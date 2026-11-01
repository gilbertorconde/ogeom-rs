//! A ball, or a sphere face, restated on a chart whose poles stand clear of the other
//! operand: of its curved faces, and of its planes other than one through
//! a pole.
//!
//! A sphere's chart collapses at its poles and folds at its seam. Where the
//! other operand's boundary passes through a pole (a drill whose wall runs
//! along the ball's axis), every section through it starts and ends on the
//! pole, hugs the seam beside it, and leaves the chart no well-placed
//! crossing to split at. The same ball charted about another axis meets the
//! same boundary with its poles well away, and the sections are ordinary
//! loops. A ball that is one face (a seam and two poles, no trim) can be
//! restated that way exactly: the surface is the same set of points, only
//! the parameterization moves. A trimmed sphere face of any solid is
//! restated too, where the new chart's seam can leave its trim at a vertex,
//! or, for a band whose two loops go round the new axis, where the seam
//! crosses each loop once: its trim edges take fitted images in the new
//! chart, and an edge the band's seam crosses is cut there, in every face
//! of the solid that holds it. A compound operand is taken member by
//! member.

use ogeom_algo::{Built, History, make_sphere};
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Curve3d as _, Surface as _};
use ogeom_geom::{PlanarCurve, SurfaceGeometry, Transformable};
use ogeom_math::{Direction, Frame, Point, Point2, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

/// How close a pole may come to the other operand's surfaces, as a
/// fraction of the ball's radius, before the chart is turned. A section a
/// thousandth and a half of the radius beside a pole is already marched
/// across the chart's collapse and snapped onto the pole.
pub const TOUCHING: f64 = 1e-2;

/// A ball operand: the solid, its one face, and the sphere it lies on in
/// world space.
struct Ball {
    solid: Shape,
    face: Shape,
    centre: Point,
    radius: f64,
    axis: Vector,
    x: Vector,
}

/// `run` on the operands, where either is or holds a whole ball or a sphere
/// face with a pole on or beside the other's surfaces, with that ball or
/// face first charted about an axis whose poles stand clear; the history
/// then runs from the operands as given. Where that refuses, and
/// otherwise, `run` on the operands as given.
pub(crate) fn with_poles_clear(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
    mut run: impl FnMut(&mut Model, &Shape, &Shape) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    #[cfg(feature = "test-hooks")]
    if crate::work::is_unturned() {
        return run(model, a, b);
    }
    if let Some((ta, tb, before)) = turned(model, a, b, tol)? {
        match run(model, &ta, &tb) {
            Ok(built) => return Ok(Built::new(built.shape, before.then(&built.history))),
            Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => return Err(e),
            Err(_) => {}
        }
    }
    run(model, a, b)
}

/// The operands with each whole ball whose pole stands on or beside the
/// other's surfaces charted about an axis whose poles stand clear, and the
/// history from the operands as given; `None` where no ball needs it.
fn turned(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape, History)>> {
    let mut before = History::new();
    let mut any = false;
    let (mut ta, mut tb) = (a.clone(), b.clone());
    for (operand, other) in [(&mut ta, b), (&mut tb, a)] {
        if let Some((made, history)) = operand_turned(model, operand, other, tol)? {
            before = before.then(&history);
            *operand = made;
            any = true;
        }
    }
    Ok(any.then_some((ta, tb, before)))
}

/// `operand` with its ball, or its sphere faces, charted clear of `other`
/// as [`turned`] states, and the history from `operand`; `None` where
/// nothing needs it. A compound is taken member by member, each placed as
/// the compound places it.
fn operand_turned(
    model: &mut Model,
    operand: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, History)>> {
    if model.kind_of(operand)? == ShapeType::Compound {
        let mut before = History::new();
        let mut any = false;
        let mut members = Vec::new();
        for member in model.children_of(operand)? {
            if let Some((made, history)) = operand_turned(model, &member, other, tol)? {
                before = before.then(&history);
                members.push(made);
                any = true;
            } else {
                members.push(member);
            }
        }
        if !any {
            return Ok(None);
        }
        let made = model.add_compound(&members)?;
        before.modify(operand, made.clone());
        return Ok(Some((made, before)));
    }
    let Some(ball) = whole_ball(model, operand, tol)? else {
        return faces_turned(model, operand, other, tol);
    };
    let walls = walls_of(model, other, tol)?;
    let clearance = |p: Point| clearance_of(&walls, p, tol);
    let poles = |axis: Vector| {
        clearance(ball.centre + axis * ball.radius).min(clearance(ball.centre - axis * ball.radius))
    };
    if poles(ball.axis) > ball.radius * TOUCHING {
        return Ok(None);
    }
    // Candidate axes in the ball's own frame: its other two axes, the
    // diagonals of each pair, and the four diagonals of the cube. The
    // one whose nearer pole stands furthest from the other's surfaces.
    let candidates = candidate_axes(ball.x, ball.axis);
    let x = ball.x;
    let y = ball.axis.cross(ball.x);
    let mut best: Option<(f64, Vector)> = None;
    for candidate in candidates {
        let Ok(axis) = candidate.normalized(tol) else {
            continue;
        };
        let score = poles(axis);
        if best.is_none_or(|(s, _)| score > s) {
            best = Some((score, axis));
        }
    }
    let Some((score, axis)) = best else {
        return Ok(None);
    };
    if score <= ball.radius * TOUCHING {
        return Ok(None);
    }
    // The seam on the side of the new axis that stands furthest from
    // the other's surfaces, read at the seam's point on the equator.
    let across = if axis.cross(x).magnitude() > 0.5 {
        x
    } else {
        y
    };
    let Ok(u) = axis.cross(across).normalized(tol) else {
        return Ok(None);
    };
    let w = axis.cross(u);
    let seam = [u, -u, w, -w]
        .into_iter()
        .map(|d| (clearance(ball.centre + d * ball.radius), d))
        .fold(None, |kept: Option<(f64, Vector)>, (c, d)| match kept {
            Some((k, _)) if k >= c => kept,
            _ => Some((c, d)),
        })
        .map_or(u, |(_, d)| d);
    let frame = Frame::new(
        ball.centre,
        Direction::new(axis, tol)?,
        Direction::new(seam, tol)?,
        tol,
    )?;
    let made = make_sphere(model, frame, ball.radius, tol)?.shape;
    let Some(face) = explore_unique(model, &made, ShapeType::Face)?
        .first()
        .cloned()
    else {
        return Ok(None);
    };
    let mut before = History::new();
    before.modify(&ball.solid, made.clone());
    for shell in explore_unique(model, &ball.solid, ShapeType::Shell)? {
        for new_shell in explore_unique(model, &made, ShapeType::Shell)? {
            before.modify(&shell, new_shell);
        }
    }
    before.modify(&ball.face, face);
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for sub in explore_unique(model, &ball.solid, kind)? {
            before.delete(&sub);
        }
    }
    Ok(Some((made, before)))
}

/// The ball `shape` is, where it is a solid of one face on a sphere with no
/// trim (a seam and its two poles), its material inside, placed as its
/// placement puts it: a ball made in its stead stands in world space, so
/// a scale or a reflection is taken whole.
fn whole_ball(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Option<Ball>> {
    if model.kind_of(shape)? != ShapeType::Solid {
        return Ok(None);
    }
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let [face] = faces.as_slice() else {
        return Ok(None);
    };
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    if !data.location.is_identity() {
        return Ok(None);
    }
    let placement = face.transform(model.datums())?;
    let Some(SurfaceGeometry::Sphere(stored)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    // The boundary: one edge with extent, used twice (the seam), and
    // degenerate edges (the poles).
    let mut seams = Vec::new();
    for wire in model.children_of(face)? {
        for edge in model.children_of(&wire)? {
            let Some(NodeData::Edge(e)) = model.node(&edge).map(|n| n.data()) else {
                return Ok(None);
            };
            if matches!(e.curve3d(), Some(EdgeRepr::Curve3d { .. })) {
                seams.push(edge.node());
            }
        }
    }
    if seams.len() != 2 || seams[0] != seams[1] {
        return Ok(None);
    }
    let sphere = stored.sphere();
    let frame = sphere.frame();
    // Material inside: the face's outward side, its natural normal or the
    // reverse, points away from the centre. Read on the sphere as stored,
    // where the face's flag states it; a placement moves the material with
    // the sphere, whichever way it turns the chart.
    let (u, v) = (0.5, 0.25);
    let at = stored.point_at(u, v, tol)?;
    let normal = stored.normal_at(u, v, tol)?.vector();
    let reversed = face.orientation() == ogeom_topo::Orientation::Reversed;
    let outward = (at - sphere.centre()).dot(normal) > 0.0;
    if outward == reversed {
        return Ok(None);
    }
    Ok(Some(Ball {
        solid: shape.clone(),
        face: face.clone(),
        centre: placement.apply(sphere.centre()),
        radius: sphere.radius() * placement.scale_factor().abs(),
        axis: placement.apply_vector(frame.z().vector()).normalized(tol)?,
        x: placement.apply_vector(frame.x().vector()).normalized(tol)?,
    }))
}

/// The surfaces of `shape`'s faces, placed in world space.
fn walls_of(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<SurfaceGeometry>> {
    let mut out = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        out.push(surface.transformed(&face.transform(model.datums())?, tol)?);
    }
    Ok(out)
}

/// How far `p`, a pole, stands from `walls`: the nearest of its distances
/// to their curved surfaces, and to their planes other than one through
/// it. A plane through the pole cuts a meridian, split there exactly; one
/// a hair beside it cuts a small circle round the pole, whose image the
/// chart cannot hold.
fn clearance_of(walls: &[SurfaceGeometry], p: Point, tol: Tolerances) -> f64 {
    walls
        .iter()
        .filter_map(|wall| match wall {
            SurfaceGeometry::Plane(plane) => {
                Some(plane.plane().distance_to(p)).filter(|d| *d > tol.confusion())
            }
            _ => ogeom_algo::project_on_surface(wall, p, 8, tol)
                .ok()
                .map(|found| found.distance),
        })
        .fold(f64::INFINITY, f64::min)
}

/// Candidate axes in a ball's own frame (`x` and its axis `z`): its three
/// axes, the diagonals of each pair, and the four diagonals of the cube.
fn candidate_axes(x: Vector, z: Vector) -> Vec<Vector> {
    let y = z.cross(x);
    let mut candidates = vec![x, y];
    for (p, q) in [(x, y), (y, z), (z, x)] {
        candidates.push(p + q);
        candidates.push(p - q);
    }
    for (s, t) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
        candidates.push(x * s + y * t + z);
    }
    candidates
}

/// How near a new pole may come to the trim of the face it is put on, as a
/// fraction of the radius, and how wide a stretch of longitude the trim
/// must leave free for a seam it does not cross. The trim's image beside a
/// pole swings fast; one well clear of it is a tame curve.
const CLEAR_OF_TRIM: f64 = 0.1;

/// Samples taken along each trim edge.
const EDGE_SAMPLES: usize = 128;

/// Samples each trim edge's image is fitted through, and the fit's target
/// at them in radians of the chart. The image is then held to a tenth of
/// its edge's tolerance in space between the samples.
const FIT_SAMPLES: usize = 512;
const FIT_IN_CHART: f64 = 1e-11;

/// How far in longitude a trim sample must stand from the seam, other than
/// at the vertex the seam leaves the face at.
const OFF_SEAM: f64 = 1e-6;

/// Meridians a band's seam is tried on, evenly spaced round the new axis,
/// and the points along each between the two loops its clearance of the
/// other operand's walls is read at.
const BAND_SEAMS: usize = 8;
const BAND_PROBES: usize = 8;

/// A shape the restating replaces and what it became, or `None` where it
/// is gone.
type Record = (Shape, Option<Shape>);

/// A trim edge cut where a new chart's seam crosses it: the edge, forward,
/// and its pieces in its own direction.
type Split = (Shape, Vec<Shape>);

/// A loop started where a band's seam crosses it: its edges, the vertex
/// it starts at, and the edge the seam cut there, if any.
type Started = (Vec<Shape>, Shape, Vec<Split>);

/// `solid` with each sphere face whose pole stands on or beside `other`'s
/// surfaces restated on a chart whose poles stand clear, and what
/// became of the shapes it replaced; `None` where no face needs it or none
/// can be restated, or the restated solid does not hold.
///
/// The face keeps its boundary edges, each given an image in the new
/// chart; its old seam and pole edges are artefacts of the old chart and
/// go. The new chart's seam crosses the trim only at a vertex: a pole
/// inside the face alone takes a seam from that vertex up to it and a pole
/// edge, two poles inside take a seam between them, and a seam the face
/// does not reach takes no edge at all. A band, no pole inside and its two
/// loops once round the axis each way, takes a seam between them on a
/// meridian crossing each loop once and reaching no other, the one of
/// [`BAND_SEAMS`] standing furthest from `other`'s walls; a loop edge it
/// crosses inside is cut there, and every face holding that edge takes the
/// pieces. A face whose trim the new chart would have to split otherwise
/// is left as it is.
///
/// A solid at a rigid placement is restated in its own frame and placed
/// as it was. One placed with a scale or a reflection, or with a pinned
/// face on a left-handed sphere, is baked into world space first, which
/// states each sphere on a right-handed frame, and the baked solid is
/// restated: baking a face already restated would re-derive images that
/// end on its seam.
fn faces_turned(
    model: &mut Model,
    solid: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, History)>> {
    if model.kind_of(solid)? != ShapeType::Solid {
        return Ok(None);
    }
    let bare = solid
        .located(ogeom_topo::Location::identity())
        .oriented(Orientation::Forward);
    let mut pinned_any = false;
    let mut left_handed = false;
    for face in explore_unique(model, &bare, ShapeType::Face)? {
        if let Some(pinned) = pinned_sphere(model, &face)? {
            pinned_any = true;
            left_handed |= pinned.sphere.frame().handedness() != ogeom_math::Handedness::Right;
        }
    }
    if !pinned_any {
        return Ok(None);
    }
    let placement = solid.transform(model.datums())?;
    let rigid =
        (placement.scale_factor().abs() - 1.0).abs() <= 1e-9 && placement.preserves_handedness();
    let (target, baked) = if rigid && !left_handed {
        (solid.clone(), None)
    } else {
        let baked = ogeom_algo::baked_shape(model, solid, tol)?;
        (baked.shape, Some(baked.history))
    };
    let Some((made, records)) = rigid_faces_turned(model, &target, other, tol)? else {
        return Ok(None);
    };
    let mut history = History::new();
    for (old, new) in records {
        match new {
            Some(new) => history.modify(&old, new),
            None => history.delete(&old),
        }
    }
    Ok(Some(match baked {
        Some(baked) => (made, baked.then(&history)),
        None => (made, history),
    }))
}

/// [`faces_turned`] on a solid at a rigid placement, with what became of
/// each shape it replaced.
fn rigid_faces_turned(
    model: &mut Model,
    solid: &Shape,
    other: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Record>)>> {
    // The solid's own node is restated in its own frame, `other`'s walls
    // brought into that frame, and the result placed as the solid is.
    let placement = solid.location().clone();
    let bare = solid
        .located(ogeom_topo::Location::identity())
        .oriented(Orientation::Forward);
    let to_local = solid.transform(model.datums())?.inverse()?;
    let mut walls: Option<Vec<SurfaceGeometry>> = None;
    let mut records = Vec::new();
    // Each trim edge a restated face cut at its seam, with the face that
    // cut it.
    let mut splits: Vec<(Split, ogeom_topo::TShapeId)> = Vec::new();
    let shells = model.children_of(&bare)?;
    let mut turned_in: Vec<Vec<Option<Shape>>> = Vec::with_capacity(shells.len());
    for shell in &shells {
        let faces = model.children_of(shell)?;
        let mut turned_here = Vec::with_capacity(faces.len());
        for face in &faces {
            let mut turned = None;
            if let Some(pinned) = pinned_sphere(model, face)?
                && pinned.sphere.frame().handedness() == ogeom_math::Handedness::Right
            {
                let walls = match &walls {
                    Some(w) => w,
                    None => walls.insert(
                        walls_of(model, other, tol)?
                            .iter()
                            .map(|w| w.transformed(&to_local, tol))
                            .collect::<OgeomResult<_>>()?,
                    ),
                };
                if let Some((new, cut)) = face_turned(model, face, &pinned, walls, tol)? {
                    records.push((face.clone(), Some(new.clone())));
                    for gone in &pinned.dropped {
                        records.push((gone.clone(), None));
                    }
                    for (old, pieces) in &cut {
                        for piece in pieces {
                            records.push((old.clone(), Some(piece.clone())));
                        }
                    }
                    splits.extend(cut.into_iter().map(|s| (s, face.node())));
                    turned = Some(new);
                }
            }
            turned_here.push(turned);
        }
        turned_in.push(turned_here);
    }
    // A restated face reads its trim edges uncut: where another restated
    // face cut one of them, the solid is left as it is.
    for (shell, turned_here) in shells.iter().zip(&turned_in) {
        for (face, turned) in model.children_of(shell)?.iter().zip(turned_here) {
            if turned.is_none() {
                continue;
            }
            for edge in explore_unique(model, face, ShapeType::Edge)? {
                if splits
                    .iter()
                    .any(|((old, _), by)| old.node() == edge.node() && *by != face.node())
                {
                    return Ok(None);
                }
            }
        }
    }
    let cut: Vec<Split> = splits.into_iter().map(|(s, _)| s).collect();
    let mut changed = false;
    let mut rebuilt_shells = Vec::with_capacity(shells.len());
    for (shell, turned_here) in shells.into_iter().zip(turned_in) {
        let faces = model.children_of(&shell)?;
        let mut kept = Vec::with_capacity(faces.len());
        let mut touched = false;
        for (face, turned) in faces.into_iter().zip(turned_here) {
            if let Some(new) = turned {
                kept.push(new);
                touched = true;
            } else if let Some(new) = with_pieces(model, &face, &cut)? {
                records.push((face, Some(new.clone())));
                kept.push(new);
                touched = true;
            } else {
                kept.push(face);
            }
        }
        if touched {
            let new_shell = model.add_shell(&kept)?;
            records.push((shell, Some(new_shell.clone())));
            rebuilt_shells.push(new_shell);
            changed = true;
        } else {
            rebuilt_shells.push(shell);
        }
    }
    if !changed {
        return Ok(None);
    }
    let shells = rebuilt_shells;
    let made = model.add_solid(&shells)?;
    if !ogeom_algo::check(model, &made, tol)?.is_valid() {
        return Ok(None);
    }
    let placed = made.moved(&placement).composed(solid.orientation());
    let mut records: Vec<Record> = records
        .into_iter()
        .map(|(old, new)| (old.moved(&placement), new.map(|n| n.moved(&placement))))
        .collect();
    records.push((solid.clone(), Some(placed.clone())));
    Ok(Some((placed, records)))
}

/// `face` on its own surface with each edge `splits` cut replaced by its
/// pieces, run the way the face runs the edge; `None` where it holds none
/// of them, or holds one at a placement.
fn with_pieces(model: &mut Model, face: &Shape, splits: &[Split]) -> OgeomResult<Option<Shape>> {
    if splits.is_empty() {
        return Ok(None);
    }
    let Some(node) = model.node(face) else {
        ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        return Ok(None);
    };
    let mut data = (**data).clone();
    data.triangulation = None;
    let wires = node.children().to_vec();
    let mut any = false;
    let mut made = Vec::with_capacity(wires.len());
    for wire in wires {
        let Some(held) = model.node(&wire) else {
            ogeom_core::ogeom_bail!(Dangling, "wire is not in this model");
        };
        let mut edges = Vec::with_capacity(held.children().len() + 1);
        for edge in held.children() {
            let Some((_, pieces)) = splits.iter().find(|(old, _)| old.node() == edge.node()) else {
                edges.push(edge.clone());
                continue;
            };
            if !edge.location().is_identity() {
                return Ok(None);
            }
            any = true;
            if edge.orientation() == Orientation::Reversed {
                edges.extend(pieces.iter().rev().map(Shape::reversed));
            } else {
                edges.extend(pieces.iter().cloned());
            }
        }
        let rebuilt = model.add_wire(&edges)?;
        made.push(
            rebuilt
                .located(wire.location().clone())
                .oriented(wire.orientation()),
        );
    }
    if !any {
        return Ok(None);
    }
    let new = model.add_face(data, &made)?;
    Ok(Some(
        new.located(face.location().clone())
            .oriented(face.orientation()),
    ))
}

/// A sphere face with a pole in its chart, its trim chained into loops.
struct Pinned {
    sphere: ogeom_math::Sphere,
    /// The face's poles in space.
    poles: Vec<Point>,
    /// The loops of its trim, seam and pole edges left out: each the edges
    /// in the order the face runs them, oriented as it uses them.
    rings: Vec<Vec<Shape>>,
    /// The seam and pole edges, and the pole vertices no loop keeps.
    dropped: Vec<Shape>,
}

/// The face as a [`Pinned`] sphere face, where it lies on a sphere at no
/// placement, has a pole, and its other edges chain head to
/// tail into closed loops one way only.
fn pinned_sphere(model: &Model, face: &Shape) -> OgeomResult<Option<Pinned>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    if !face.location().is_identity() || !data.location.is_identity() {
        return Ok(None);
    }
    let Some(SurfaceGeometry::Sphere(stored)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let sphere = stored.sphere();
    // The wires as the face stores them, whichever way the shell uses it.
    let forward = face.oriented(Orientation::Forward);
    let mut uses: ogeom_core::FastMap<ogeom_topo::TShapeId, usize> = ogeom_core::FastMap::default();
    let mut edges = Vec::new();
    for wire in model.children_of(&forward)? {
        for edge in model.children_of(&wire)? {
            *uses.entry(edge.node()).or_default() += 1;
            edges.push(edge);
        }
    }
    let mut poles = Vec::new();
    let mut dropped: Vec<Shape> = Vec::new();
    let mut kept: Vec<(Shape, Shape, Shape)> = Vec::new();
    for edge in edges {
        let Some(NodeData::Edge(e)) = model.node(&edge).map(|n| n.data()) else {
            return Ok(None);
        };
        if !edge.location().is_identity() {
            return Ok(None);
        }
        let Some((start, end)) = ogeom_algo::edge_vertices(model, &edge)? else {
            return Ok(None);
        };
        if e.degenerate {
            let Some(NodeData::Vertex(v)) = model.node(&start).map(|n| n.data()) else {
                return Ok(None);
            };
            if !start.location().is_identity() {
                return Ok(None);
            }
            poles.push(v.point);
            if !dropped.iter().any(|d| d.node() == edge.node()) {
                dropped.push(edge.clone());
                dropped.push(start);
            }
            continue;
        }
        match uses.get(&edge.node()) {
            Some(1) => {}
            // The seam, up one side of the chart and down the other.
            Some(2) => {
                if !dropped.iter().any(|d| d.node() == edge.node()) {
                    dropped.push(edge.clone());
                }
                continue;
            }
            _ => return Ok(None),
        }
        if !matches!(e.curve3d(), Some(EdgeRepr::Curve3d { location, .. }) if location.is_identity())
        {
            return Ok(None);
        }
        kept.push((edge, start, end));
    }
    if poles.is_empty() || kept.is_empty() {
        return Ok(None);
    }
    // A pole vertex a kept edge still ends at is a corner of the trim.
    dropped.retain(|d| {
        !kept
            .iter()
            .any(|(_, s, e)| s.node() == d.node() || e.node() == d.node())
    });
    let mut rings = Vec::new();
    let mut used = vec![false; kept.len()];
    while let Some(first) = used.iter().position(|u| !u) {
        used[first] = true;
        let origin = kept[first].1.node();
        let mut ring = vec![kept[first].0.clone()];
        let mut at = kept[first].2.node();
        while at != origin {
            let next: Vec<usize> = (0..kept.len())
                .filter(|&i| !used[i] && kept[i].1.node() == at)
                .collect();
            let [i] = next.as_slice() else {
                return Ok(None);
            };
            used[*i] = true;
            ring.push(kept[*i].0.clone());
            at = kept[*i].2.node();
        }
        rings.push(ring);
    }
    Ok(Some(Pinned {
        sphere,
        poles,
        rings,
        dropped,
    }))
}

/// The face restated on its sphere charted about an axis whose poles stand
/// clear of `walls` and of the face's own trim, where a pole of its chart
/// stands on or beside `walls`, and the trim edges its new seam cut;
/// `None` otherwise, or where no candidate axis lays a chart the trim
/// fits.
fn face_turned(
    model: &mut Model,
    face: &Shape,
    pinned: &Pinned,
    walls: &[SurfaceGeometry],
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Split>)>> {
    // The curved walls rank the candidate axes, as for a face no plane
    // stands beside; the planes, a pole beside which turns the chart as a
    // curved wall's does, only bar an axis whose pole they pass beside.
    let (flats, curved): (Vec<SurfaceGeometry>, Vec<SurfaceGeometry>) = walls
        .iter()
        .cloned()
        .partition(|w| matches!(w, SurfaceGeometry::Plane(_)));
    let clearance = |p: Point| clearance_of(&curved, p, tol);
    let (centre, radius) = (pinned.sphere.centre(), pinned.sphere.radius());
    let beside = |p: Point| clearance_of(walls, p, tol) <= radius * TOUCHING;
    if !pinned.poles.iter().any(|&p| beside(p)) {
        return Ok(None);
    }
    let mut trim: Vec<Point> = Vec::new();
    for ring in &pinned.rings {
        trim.extend(ring_points(model, ring, tol)?);
    }
    let near_trim = |p: Point| -> f64 {
        trim.iter()
            .map(|q| q.distance(p))
            .fold(f64::INFINITY, f64::min)
    };
    let frame = pinned.sphere.frame();
    let mut scored: Vec<(f64, Vector)> = Vec::new();
    for candidate in candidate_axes(frame.x().vector(), frame.z().vector()) {
        let Ok(axis) = candidate.normalized(tol) else {
            continue;
        };
        let ends = [centre + axis * radius, centre - axis * radius];
        if ends.iter().any(|&p| near_trim(p) <= radius * CLEAR_OF_TRIM) {
            continue;
        }
        if ends
            .iter()
            .any(|&p| clearance_of(&flats, p, tol) <= radius * TOUCHING)
        {
            continue;
        }
        let score = clearance(ends[0]).min(clearance(ends[1]));
        if score > radius * TOUCHING {
            scored.push((score, axis));
        }
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (_, axis) in scored {
        if let Some(made) = restated(model, face, pinned, axis, &clearance, tol)? {
            return Ok(Some(made));
        }
    }
    Ok(None)
}

/// `count` steps of an edge's parameters and points, in the order the edge
/// as given runs.
fn edge_samples(
    model: &Model,
    edge: &Shape,
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, Point)>> {
    let Some(NodeData::Edge(e)) = model.node(edge).map(|n| n.data()) else {
        return Ok(Vec::new());
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = e.curve3d() else {
        return Ok(Vec::new());
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(count + 1);
    for k in 0..=count {
        #[allow(clippy::cast_precision_loss)]
        let s = k as f64 / count as f64;
        let s = if edge.orientation() == Orientation::Reversed {
            1.0 - s
        } else {
            s
        };
        let t = range.0 + (range.1 - range.0) * s;
        out.push((t, curve.point_at(t, tol)?));
    }
    Ok(out)
}

/// A loop's points in the order it runs.
fn ring_points(model: &Model, ring: &[Shape], tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let mut out = Vec::new();
    for edge in ring {
        out.extend(
            edge_samples(model, edge, EDGE_SAMPLES, tol)?
                .into_iter()
                .map(|(_, p)| p),
        );
    }
    Ok(out)
}

/// Longitude and latitude of `p` in a sphere chart about `frame`.
fn chart_of(frame: &Frame, p: Point) -> (f64, f64) {
    let local = frame.to_local(p);
    (
        local.y.atan2(local.x),
        local.z.atan2(local.x.hypot(local.y)),
    )
}

/// `u` moved by whole turns to within half a turn of `near`.
fn unwrapped(u: f64, near: f64) -> f64 {
    let turn = core::f64::consts::TAU;
    u - turn * ((u - near) / turn).round()
}

/// The longitudes of a loop of points in a chart about `frame`, carried
/// along the loop from the first one.
fn longitudes(frame: &Frame, ring: &[Point]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(ring.len());
    for &p in ring {
        let u = chart_of(frame, p).0;
        out.push(out.last().map_or(u, |&before| unwrapped(u, before)));
    }
    out
}

/// Whole turns a loop of longitudes makes about the axis.
fn winding(us: &[f64]) -> i64 {
    let (Some(first), Some(last)) = (us.first(), us.last()) else {
        return 0;
    };
    #[allow(clippy::cast_possible_truncation)]
    let turns = ((last - first) / core::f64::consts::TAU).round() as i64;
    turns
}

/// The least and greatest of some longitudes.
fn span(us: &[f64]) -> (f64, f64) {
    us.iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &u| {
            (lo.min(u), hi.max(u))
        })
}

/// Whether some turn of longitude `u` lies within `(lo, hi)` widened by
/// [`OFF_SEAM`].
fn reaches(lo: f64, hi: f64, u: f64) -> bool {
    let turn = core::f64::consts::TAU;
    let first = u + turn * ((lo - OFF_SEAM - u) / turn).ceil();
    first <= hi + OFF_SEAM
}

/// The index of the sample beside the one place a loop's carried
/// longitudes pass longitude `u` (some whole turn of it), or `None` where
/// they pass it other than once.
fn crossed_once(us: &[f64], u: f64) -> Option<usize> {
    let turn = core::f64::consts::TAU;
    let mut count = 0.0;
    let mut found = None;
    for (i, w) in us.windows(2).enumerate() {
        let (lo, hi) = if w[0] <= w[1] {
            (w[0], w[1])
        } else {
            (w[1], w[0])
        };
        // The turns of `u` in `(lo, hi]`.
        let passed = (((hi - u) / turn).floor() - ((lo - u) / turn).floor()).max(0.0);
        if passed > 0.0 {
            count += passed;
            found = Some(i);
        }
    }
    if count == 1.0 { found } else { None }
}

/// The surface a face lies on.
fn face_surface(model: &Model, face: &Shape) -> Option<ogeom_topo::SurfaceId> {
    match model.node(face).map(|n| n.data()) {
        Some(NodeData::Face(data)) => Some(data.surface),
        _ => None,
    }
}

/// `ring`, a loop once round `chart`'s axis, started at the vertex where it
/// crosses the chart's seam meridian (longitude 0), the edge it crosses
/// inside cut there: the loop, that vertex, and the cut. `None` where it
/// crosses that meridian other than once, or the edge will not cut.
fn started_on_seam(
    model: &mut Model,
    ring: &[Shape],
    chart: &Frame,
    old_surface: ogeom_topo::SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Option<Started>> {
    use core::f64::consts::FRAC_PI_2;
    let mut samples: Vec<(usize, f64, f64)> = Vec::new();
    for (k, edge) in ring.iter().enumerate() {
        for (t, p) in edge_samples(model, edge, EDGE_SAMPLES, tol)? {
            samples.push((k, t, chart_of(chart, p).0));
        }
    }
    // Longitude changes sign at the seam and, a half turn round, where it
    // wraps; only the seam's is wanted.
    let crossings: Vec<(usize, f64, f64, f64)> = samples
        .windows(2)
        .filter(|w| {
            w[0].0 == w[1].0
                && (w[0].2 < 0.0) != (w[1].2 < 0.0)
                && w[0].2.abs() < FRAC_PI_2
                && w[1].2.abs() < FRAC_PI_2
        })
        .map(|w| (w[0].0, w[0].1, w[1].1, w[0].2))
        .collect();
    let [(k, lo, hi, u_lo)] = crossings[..] else {
        return Ok(None);
    };
    let (mut lo, mut hi) = (lo, hi);
    let Some(NodeData::Edge(data)) = model.node(&ring[k]).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
        return Ok(None);
    };
    let Some(curve) = model.geometry().curve(*curve).cloned() else {
        return Ok(None);
    };
    let below = u_lo < 0.0;
    for _ in 0..200 {
        let mid = f64::midpoint(lo, hi);
        if mid == lo || mid == hi {
            break;
        }
        if (chart_of(chart, curve.point_at(mid, tol)?).0 < 0.0) == below {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = f64::midpoint(lo, hi);
    let at = curve.point_at(t, tol)?;
    let Some((start, end)) = ogeom_algo::edge_vertices(model, &ring[k])? else {
        return Ok(None);
    };
    for (vertex, from) in [(&start, k), (&end, k + 1)] {
        let Some(NodeData::Vertex(v)) = model.node(vertex).map(|n| n.data()) else {
            return Ok(None);
        };
        if v.point.distance(at) <= v.tolerance.get().max(tol.confusion()) {
            let mut started = ring.to_vec();
            started.rotate_left(from % ring.len());
            return Ok(Some((started, vertex.clone(), Vec::new())));
        }
    }
    let edge = ring[k].oriented(Orientation::Forward);
    let Some((pieces, vertex)) = cut_edge(model, &edge, t, old_surface, tol)? else {
        return Ok(None);
    };
    let [head, tail] = pieces.as_slice() else {
        return Ok(None);
    };
    let run = if ring[k].orientation() == Orientation::Reversed {
        [tail.reversed(), head.reversed()]
    } else {
        [head.clone(), tail.clone()]
    };
    let mut started: Vec<Shape> = ring[..k].to_vec();
    started.extend(run);
    started.extend(ring[k + 1..].iter().cloned());
    started.rotate_left(k + 1);
    Ok(Some((started, vertex, vec![(edge, pieces)])))
}

/// `edge` (forward) cut at `t` on its curve: its two pieces in its own
/// direction, each with the edge's pcurves on surfaces other than `skip`
/// cut where they pass the cut point, and the vertex between them. `None`
/// where the edge is degenerate or placed, it is a seam of another
/// surface, or a pcurve passes nowhere within its tolerance of that point.
fn cut_edge(
    model: &mut Model,
    edge: &Shape,
    t: f64,
    skip: ogeom_topo::SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Shape>, Shape)>> {
    let Some(NodeData::Edge(data)) = model.node(edge).map(|n| n.data()) else {
        return Ok(None);
    };
    let data = data.clone();
    if data.degenerate {
        return Ok(None);
    }
    let Some(EdgeRepr::Curve3d {
        curve,
        location,
        range,
    }) = data.curve3d()
    else {
        return Ok(None);
    };
    if !location.is_identity() || !(range.0 < t && t < range.1) {
        return Ok(None);
    }
    let range = *range;
    let Some(geometry) = model.geometry().curve(*curve).cloned() else {
        return Ok(None);
    };
    let Some((start, end)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    let at = geometry.point_at(t, tol)?;
    let reach = data.tolerance.get().max(tol.confusion());
    let mut cuts = Vec::new();
    for repr in &data.representations {
        match repr {
            EdgeRepr::PCurve {
                curve,
                surface,
                location,
                range,
            } => {
                if *surface == skip {
                    continue;
                }
                if !location.is_identity() {
                    return Ok(None);
                }
                let (Some(pcurve), Some(on)) = (
                    model.geometry().pcurve(*curve),
                    model.geometry().surface(*surface),
                ) else {
                    return Ok(None);
                };
                let Some(s) = pcurve_parameter_at(on, pcurve, *range, at, reach, tol)? else {
                    return Ok(None);
                };
                if !(range.0 < s && s < range.1) {
                    return Ok(None);
                }
                cuts.push((pcurve.clone(), *surface, (range.0, s), (s, range.1)));
            }
            EdgeRepr::Seam { surface, .. } if *surface != skip => return Ok(None),
            _ => {}
        }
    }
    let vertex = model.add_vertex(ogeom_topo::VertexData::new(at));
    let head =
        ogeom_algo::make_edge_between(model, geometry.clone(), (range.0, t), &start, &vertex, tol)?
            .shape;
    let tail =
        ogeom_algo::make_edge_between(model, geometry, (t, range.1), &vertex, &end, tol)?.shape;
    for (pcurve, surface, first, second) in cuts {
        let identity = ogeom_topo::Location::identity();
        ogeom_algo::attach_pcurve(
            model,
            &head,
            pcurve.clone(),
            surface,
            identity.clone(),
            first,
        )?;
        ogeom_algo::attach_pcurve(model, &tail, pcurve, surface, identity, second)?;
    }
    model.widen(&head, data.tolerance)?;
    model.widen(&tail, data.tolerance)?;
    Ok(Some((vec![head, tail], vertex)))
}

/// The parameter in `range` at which `pcurve` on `surface` passes nearest
/// `p`, where it passes within `reach` of it: the nearest of evenly spaced
/// samples, refined by golden section between its neighbours.
fn pcurve_parameter_at(
    surface: &SurfaceGeometry,
    pcurve: &PlanarCurve,
    range: (f64, f64),
    p: Point,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let gap = |s: f64| -> OgeomResult<f64> {
        let uv = pcurve.point_at(s, tol)?;
        Ok(surface.point_at(uv.x, uv.y, tol)?.distance(p))
    };
    #[allow(clippy::cast_precision_loss)]
    let step = |k: usize| range.0 + (range.1 - range.0) * k as f64 / FIT_SAMPLES as f64;
    let mut nearest = (f64::INFINITY, 0);
    for k in 0..=FIT_SAMPLES {
        let d = gap(step(k))?;
        if d < nearest.0 {
            nearest = (d, k);
        }
    }
    let (mut lo, mut hi) = (
        step(nearest.1.saturating_sub(1)),
        step((nearest.1 + 1).min(FIT_SAMPLES)),
    );
    let golden = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..100 {
        let (a, b) = (hi - golden * (hi - lo), lo + golden * (hi - lo));
        if gap(a)? < gap(b)? {
            hi = b;
        } else {
            lo = a;
        }
    }
    let s = f64::midpoint(lo, hi);
    Ok((gap(s)? <= reach).then_some(s))
}

/// The face restated on its sphere charted about `axis`, and the trim
/// edges a band's seam cut, where the trim fits that chart crossing its
/// seam only at a vertex, or once across each of a band's two loops;
/// `None` where it does not, or an image of the trim misses its edge.
fn restated(
    model: &mut Model,
    face: &Shape,
    pinned: &Pinned,
    axis: Vector,
    clearance: &dyn Fn(Point) -> f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Split>)>> {
    use core::f64::consts::{FRAC_PI_2, TAU};
    use ogeom_algo::Containment;
    let (centre, radius) = (pinned.sphere.centre(), pinned.sphere.radius());
    let deflection = ogeom_mesh::Deflection::default();
    let mut inside = [false; 2];
    for (k, sign) in [1.0, -1.0].into_iter().enumerate() {
        let pole = centre + axis * (radius * sign);
        match ogeom_algo::classify_on_face(model, face, pole, deflection, tol)? {
            Containment::In => inside[k] = true,
            Containment::Out => {}
            Containment::On => return Ok(None),
        }
    }
    // A pole inside alone is put north.
    let axis = if inside == [false, true] { -axis } else { axis };
    let poles_in = inside.iter().filter(|&&i| i).count();
    let reference = pinned.sphere.frame();
    let Ok(across) = axis
        .cross(reference.x().vector())
        .normalized(tol)
        .or_else(|_| axis.cross(reference.y().vector()).normalized(tol))
    else {
        return Ok(None);
    };
    let provisional = Frame::new(
        centre,
        Direction::new(axis, tol)?,
        Direction::new(across, tol)?,
        tol,
    )?;
    let mut points: Vec<Vec<Point>> = Vec::with_capacity(pinned.rings.len());
    for ring in &pinned.rings {
        points.push(ring_points(model, ring, tol)?);
    }
    let turns: Vec<i64> = points
        .iter()
        .map(|r| winding(&longitudes(&provisional, r)))
        .collect();
    let around: Vec<usize> = (0..turns.len()).filter(|&i| turns[i] != 0).collect();
    // The seam's longitude in the provisional chart, and where a pole is
    // inside alone, the loop round it and the edge whose start the seam
    // leaves the face at.
    let (seam_u, through) = if poles_in == 1 {
        let [outer] = around.as_slice() else {
            return Ok(None);
        };
        let outer = *outer;
        if turns[outer].abs() != 1 {
            return Ok(None);
        }
        let mut found = None;
        for k in 0..pinned.rings[outer].len() {
            // The loop from this edge's start round to it again crosses
            // the meridian there only, and no other loop reaches it.
            let mut ring = pinned.rings[outer].clone();
            ring.rotate_left(k);
            let us = longitudes(&provisional, &ring_points(model, &ring, tol)?);
            let Some((&first, rest)) = us.split_first() else {
                return Ok(None);
            };
            let inner = &rest[..rest.len().saturating_sub(1)];
            let (lo, hi) = span(inner);
            let (below, above) = if turns[outer] > 0 {
                (first, first + TAU)
            } else {
                (first - TAU, first)
            };
            let alone = lo > below + OFF_SEAM && hi < above - OFF_SEAM;
            let clear = alone
                && points.iter().enumerate().all(|(i, r)| {
                    if i == outer {
                        return true;
                    }
                    let (lo, hi) = span(&longitudes(&provisional, r));
                    !reaches(lo, hi, first)
                });
            if clear {
                found = Some((first, Some((outer, k))));
                break;
            }
        }
        let Some(found) = found else {
            return Ok(None);
        };
        found
    } else if poles_in == 0 && !around.is_empty() {
        // A band: the two loops round the axis, the seam on the meridian
        // standing furthest from the walls of those that cross each loop
        // once and reach no other loop.
        let [first, second] = around.as_slice() else {
            return Ok(None);
        };
        if turns[*first].abs() != 1 || turns[*second] != -turns[*first] {
            return Ok(None);
        }
        let carried: Vec<Vec<f64>> = points.iter().map(|r| longitudes(&provisional, r)).collect();
        let mut best: Option<(f64, f64)> = None;
        for k in 0..BAND_SEAMS {
            #[allow(clippy::cast_precision_loss)]
            let u = TAU * k as f64 / BAND_SEAMS as f64;
            let mut latitudes = Vec::with_capacity(2);
            let mut fits = true;
            for (i, us) in carried.iter().enumerate() {
                if !around.contains(&i) {
                    let (lo, hi) = span(us);
                    fits &= !reaches(lo, hi, u);
                    continue;
                }
                match crossed_once(us, u) {
                    Some(at) => latitudes.push(chart_of(&provisional, points[i][at]).1),
                    None => fits = false,
                }
            }
            let [a, b] = latitudes.as_slice() else {
                continue;
            };
            if !fits {
                continue;
            }
            let meridian = across * u.cos() + axis.cross(across) * u.sin();
            let score = (0..=BAND_PROBES)
                .map(|j| {
                    #[allow(clippy::cast_precision_loss)]
                    let v = a + (b - a) * j as f64 / BAND_PROBES as f64;
                    clearance(centre + (meridian * v.cos() + axis * v.sin()) * radius)
                })
                .fold(f64::INFINITY, f64::min);
            if best.is_none_or(|(s, _)| score > s) {
                best = Some((score, u));
            }
        }
        let Some((_, u)) = best else {
            return Ok(None);
        };
        (u, None)
    } else {
        if !around.is_empty() {
            return Ok(None);
        }
        // The middle of the widest stretch of longitude no loop reaches.
        let mut spans: Vec<(f64, f64)> = points
            .iter()
            .map(|r| {
                let (lo, hi) = span(&longitudes(&provisional, r));
                let start = lo.rem_euclid(TAU);
                (start, start + (hi - lo))
            })
            .collect();
        spans.sort_by(|a, b| a.0.total_cmp(&b.0));
        let Some(&(first, _)) = spans.first() else {
            return Ok(None);
        };
        let mut best: Option<(f64, f64)> = None;
        let mut reach = f64::NEG_INFINITY;
        for &(lo, hi) in &spans {
            if lo > reach && reach > f64::NEG_INFINITY && best.is_none_or(|(w, _)| lo - reach > w) {
                best = Some((lo - reach, f64::midpoint(lo, reach)));
            }
            reach = reach.max(hi);
        }
        let wrap = first + TAU - reach;
        if best.is_none_or(|(w, _)| wrap > w) {
            best = Some((wrap, f64::midpoint(reach, first + TAU)));
        }
        let Some((width, middle)) = best else {
            return Ok(None);
        };
        if width <= CLEAR_OF_TRIM {
            return Ok(None);
        }
        (middle, None)
    };
    let seam_x = across * seam_u.cos() + axis.cross(across) * seam_u.sin();
    let chart = Frame::new(
        centre,
        Direction::new(axis, tol)?,
        Direction::new(seam_x, tol)?,
        tol,
    )?;
    let sphere = ogeom_math::Sphere::new(chart, radius, tol)?;
    let geometry: SurfaceGeometry = ogeom_geom::SphereSurface::new(sphere).into();
    let surface = model.geometry_mut().add_surface(geometry.clone());

    // A band's two loops cut where the seam crosses them, each then
    // starting there, with the vertex and its latitude.
    let mut rings = pinned.rings.clone();
    let mut splits: Vec<Split> = Vec::new();
    let mut band_ends: Vec<(usize, Shape, f64)> = Vec::new();
    if poles_in == 0 && !around.is_empty() {
        let Some(old_surface) = face_surface(model, face) else {
            return Ok(None);
        };
        for &i in &around {
            let Some((ring, vertex, cut)) =
                started_on_seam(model, &rings[i], &chart, old_surface, tol)?
            else {
                return Ok(None);
            };
            let Some(NodeData::Vertex(v)) = model.node(&vertex).map(|n| n.data()) else {
                return Ok(None);
            };
            band_ends.push((i, vertex.clone(), chart_of(&chart, v.point).1));
            rings[i] = ring;
            splits.extend(cut);
        }
    }

    // Each loop's edges given their image, the longitude carried along the
    // loop from where it starts: the vertex on the seam for a loop round
    // the axis, else wherever puts the whole loop inside one turn.
    let mut wires: Vec<(Vec<Shape>, f64)> = Vec::new();
    let mut seam_ends: Option<(f64, f64)> = None;
    for (i, ring) in rings.iter().enumerate() {
        let mut ring = ring.clone();
        if let Some((outer, k)) = through
            && outer == i
        {
            ring.rotate_left(k);
        }
        let mut samples: Vec<Vec<(f64, Point)>> = Vec::with_capacity(ring.len());
        for edge in &ring {
            samples.push(edge_samples(model, edge, FIT_SAMPLES, tol)?);
        }
        let flat: Vec<Point> = samples.iter().flatten().map(|(_, p)| *p).collect();
        let us = longitudes(&chart, &flat);
        let (Some(&first), Some(&last)) = (us.first(), us.last()) else {
            return Ok(None);
        };
        let turn = winding(&us);
        let shift = match turn {
            0 => -TAU * (span(&us).0 / TAU).floor(),
            1 => -first,
            -1 => TAU - first,
            _ => return Ok(None),
        };
        if turn != 0 {
            seam_ends = Some((first + shift, last + shift));
        }
        let mut area = 0.0;
        let mut at = 0;
        for (edge, along) in ring.iter().zip(&samples) {
            let mut params = Vec::with_capacity(along.len());
            let mut image_points = Vec::with_capacity(along.len());
            for (j, &(t, p)) in along.iter().enumerate() {
                params.push(t);
                image_points.push(Point2::new(us[at + j] + shift, chart_of(&chart, p).1));
            }
            at += along.len();
            for w in image_points.windows(2) {
                area += (w[0].x * w[1].y - w[1].x * w[0].y) / 2.0;
            }
            if edge.orientation() == Orientation::Reversed {
                params.reverse();
                image_points.reverse();
            }
            let (Some(&lo), Some(&hi)) = (params.first(), params.last()) else {
                return Ok(None);
            };
            let fitted =
                ogeom_geom::fit::fit_points_2d_at(&params, &image_points, 3, FIT_IN_CHART, tol)?;
            let image: PlanarCurve = fitted.curve.into();
            // The image is honest in space between its samples too.
            let Some(NodeData::Edge(e)) = model.node(edge).map(|n| n.data()) else {
                return Ok(None);
            };
            let Some(EdgeRepr::Curve3d { curve, .. }) = e.curve3d() else {
                return Ok(None);
            };
            let Some(curve) = model.geometry().curve(*curve) else {
                return Ok(None);
            };
            let budget = e.tolerance.get().max(tol.confusion()) * 0.1;
            for pair in params.windows(2) {
                let t = f64::midpoint(pair[0], pair[1]);
                let uv = image.point_at(t, tol)?;
                let Ok(on) = geometry.point_at(uv.x, uv.y, tol) else {
                    return Ok(None);
                };
                if on.distance(curve.point_at(t, tol)?) > budget {
                    return Ok(None);
                }
            }
            ogeom_algo::attach_pcurve(
                model,
                edge,
                image,
                surface,
                ogeom_topo::Location::identity(),
                (lo, hi),
            )?;
        }
        wires.push((ring, area));
    }

    // The seam and pole edges the new chart needs, and the outer loop
    // first.
    let (outer, holes): (Vec<Shape>, Vec<Vec<Shape>>) = match (poles_in, through) {
        (0, None) if !band_ends.is_empty() => {
            let [low, high] = band_ends.as_slice() else {
                return Ok(None);
            };
            let (low, high) = if low.2 <= high.2 {
                (low, high)
            } else {
                (high, low)
            };
            // Up the seam where the lower loop ends, down where it starts.
            let sides = if turns[low.0] > 0 {
                (TAU, 0.0)
            } else {
                (0.0, TAU)
            };
            let seam = seam_between(
                model,
                surface,
                &chart,
                radius,
                (&low.1, low.2),
                (&high.1, high.2),
                sides,
                tol,
            )?;
            let mut edges = wires[low.0].0.clone();
            edges.push(seam.clone());
            edges.extend(wires[high.0].0.iter().cloned());
            edges.push(seam.reversed());
            let holes = wires
                .into_iter()
                .enumerate()
                .filter(|(i, _)| *i != low.0 && *i != high.0)
                .map(|(_, (w, _))| w)
                .collect();
            (edges, holes)
        }
        (1, Some((ring_index, _))) => {
            let Some((u_start, u_end)) = seam_ends else {
                return Ok(None);
            };
            let ring = wires[ring_index].0.clone();
            let Some((start, _)) = ogeom_algo::edge_vertices(model, &ring[0])? else {
                return Ok(None);
            };
            let Some(NodeData::Vertex(v)) = model.node(&start).map(|n| n.data()) else {
                return Ok(None);
            };
            let low = chart_of(&chart, v.point).1;
            let north = model.add_vertex(ogeom_topo::VertexData::new(centre + axis * radius));
            let seam = seam_between(
                model,
                surface,
                &chart,
                radius,
                (&start, low),
                (&north, FRAC_PI_2),
                (u_end, u_start),
                tol,
            )?;
            let top = pole_at(model, surface, &north, FRAC_PI_2, tol)?;
            let mut edges = ring;
            edges.push(seam.clone());
            // Across the top from where the loop ended back to where it
            // started.
            edges.push(if u_end > u_start { top.reversed() } else { top });
            edges.push(seam.reversed());
            let holes = wires
                .into_iter()
                .enumerate()
                .filter(|(i, _)| *i != ring_index)
                .map(|(_, (w, _))| w)
                .collect();
            (edges, holes)
        }
        (2, None) => {
            let Some(hole_area) = wires.first().map(|(_, a)| *a) else {
                return Ok(None);
            };
            // Counter-clockwise round the chart where the holes run the
            // other way.
            let ccw = hole_area < 0.0;
            let south = model.add_vertex(ogeom_topo::VertexData::new(centre - axis * radius));
            let north = model.add_vertex(ogeom_topo::VertexData::new(centre + axis * radius));
            let sides = if ccw { (TAU, 0.0) } else { (0.0, TAU) };
            let seam = seam_between(
                model,
                surface,
                &chart,
                radius,
                (&south, -FRAC_PI_2),
                (&north, FRAC_PI_2),
                sides,
                tol,
            )?;
            let top = pole_at(model, surface, &north, FRAC_PI_2, tol)?;
            let bottom = pole_at(model, surface, &south, -FRAC_PI_2, tol)?;
            let edges = if ccw {
                vec![bottom, seam.clone(), top.reversed(), seam.reversed()]
            } else {
                vec![seam.clone(), top, seam.reversed(), bottom.reversed()]
            };
            (edges, wires.into_iter().map(|(w, _)| w).collect())
        }
        (0, None) => {
            let Some(widest) =
                (0..wires.len()).max_by(|&a, &b| wires[a].1.abs().total_cmp(&wires[b].1.abs()))
            else {
                return Ok(None);
            };
            let mut outer = Vec::new();
            let mut holes = Vec::new();
            for (i, (w, _)) in wires.into_iter().enumerate() {
                if i == widest {
                    outer = w;
                } else {
                    holes.push(w);
                }
            }
            (outer, holes)
        }
        _ => return Ok(None),
    };
    let mut made_wires = vec![ogeom_algo::make_wire(model, &outer, tol)?.shape];
    for hole in &holes {
        made_wires.push(ogeom_algo::make_wire(model, hole, tol)?.shape);
    }
    let Some(NodeData::Face(old)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let held = old.tolerance;
    let data = ogeom_topo::FaceData::new(surface, ogeom_topo::Location::identity());
    let made = model.add_face(data, &made_wires)?;
    model.widen(&made, held)?;
    Ok(Some((made.oriented(face.orientation()), splits)))
}

/// A pole at `at`, latitude `v` of a sphere chart: an edge with no length,
/// its image across the chart from longitude 0 to a whole turn.
fn pole_at(
    model: &mut Model,
    surface: ogeom_topo::SurfaceId,
    at: &Shape,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut data = ogeom_topo::EdgeData::new();
    data.degenerate = true;
    let edge = model.add_edge(data, &[at.clone(), at.clone()])?;
    let (from, to) = (Point2::new(0.0, v), Point2::new(core::f64::consts::TAU, v));
    ogeom_algo::attach_pcurve(
        model,
        &edge,
        ogeom_geom::Line2d::segment(from, to, tol)?.into(),
        surface,
        ogeom_topo::Location::identity(),
        (0.0, from.distance(to)),
    )?;
    Ok(edge)
}

/// The seam of a sphere chart about `chart` from `low` up to `high` (each
/// a vertex and its latitude), its image up the side of the chart at
/// longitude `sides.0` where the edge runs forward and at `sides.1` where
/// it runs back.
#[allow(clippy::too_many_arguments)]
fn seam_between(
    model: &mut Model,
    surface: ogeom_topo::SurfaceId,
    chart: &Frame,
    radius: f64,
    low: (&Shape, f64),
    high: (&Shape, f64),
    sides: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let ((from, bottom), (to, top)) = (low, high);
    // Its plane is spanned by the chart's x and axis, so the circle's
    // normal is -y and its angle is the latitude.
    let meridian = ogeom_math::Circle::new(
        Frame::new(chart.origin(), -chart.y(), chart.x(), tol)?,
        radius,
        tol,
    )?;
    let seam = ogeom_algo::make_edge_between(
        model,
        ogeom_geom::CircleCurve::new(meridian).into(),
        (bottom, top),
        from,
        to,
        tol,
    )?
    .shape;
    let side = |u: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::segment(Point2::new(u, bottom), Point2::new(u, top), tol)?.into())
    };
    ogeom_algo::attach_seam(
        model,
        &seam,
        side(sides.0)?,
        side(sides.1)?,
        surface,
        ogeom_topo::Location::identity(),
        (0.0, top - bottom),
    )?;
    Ok(seam)
}
