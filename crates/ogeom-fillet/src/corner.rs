//! Rounding a vertex: the ball-and-block tool at a planar corner.
//!
//! The ball's centre may sit anywhere a radius in from every host plane.
//! At a corner one ball touches, that region's tip is a point and the
//! rounded corner is one spherical patch; at a corner no single ball
//! touches, the tip is a few points joined by short ridges, and the
//! rounded corner is a sphere at each and a cylinder along each ridge:
//! the exact envelope of the rolling ball, where the setback family fits
//! a plate through the bands' ends instead.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// [`crate::round_vertex`], with the one-ball corner held to the labelling
/// `forced` names (by index, the ring's start then its direction) instead
/// of offered each in turn.
pub(crate) fn round_vertex_with(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    forced: Option<usize>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = crate::fillet::to_current_path("vertex:round_vertex_with");
    let mut built = match corner_at(model, solid, vertex, radius, tol)? {
        Corner::Curved(corner) => {
            crate::corner_curved::curved_corner(model, solid, vertex, corner, radius, tol)
        }
        Corner::Setback(corner, m) => setback_corner(model, solid, vertex, corner, &m, radius, tol),
        Corner::OneBall(one) => one_ball_corner(model, solid, vertex, &one, radius, forced, tol),
    }?;
    // A vertex's solid image supports the corner tool's callers, but it does
    // not identify the patch for a client naming faces. Record the new faces
    // separately, before later bands split or consume them. This uses boolean
    // ancestry, not the patch's surface kind or distance from the vertex.
    for face in corner_faces(model, solid, &built)? {
        built.history.generate(vertex, face);
    }
    Ok(built)
}

/// Faces introduced by a corner, excluding every image of an input face.
fn corner_faces(model: &Model, solid: &Shape, built: &Built) -> OgeomResult<Vec<Shape>> {
    let inherited: Vec<Shape> = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?
        .iter()
        .flat_map(|face| built.history.trace(face).to_vec())
        .collect();
    Ok(
        ogeom_topo::explore_unique(model, &built.shape, ShapeType::Face)?
            .into_iter()
            .filter(|face| !inherited.iter().any(|image| image.is_same(face)))
            .collect(),
    )
}

use crate::engine::vertex_read::{Corner, OneBall, corner_at};
pub(crate) use crate::engine::vertex_read::{CornerFrame, corner_frame};

/// [`round_vertex_with`] at a corner one ball touches every face of.
fn one_ball_corner(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    one: &OneBall,
    radius: f64,
    forced: Option<usize>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let OneBall {
        corner,
        far,
        ref m,
        ref ring,
        ref d,
        ..
    } = *one;
    let n = m.len();

    // A band left short of the block by another band's run-out is carried
    // back to it first, so the block's cut ends it as it ends a flush one.
    // The wedge carrying it starts at the vertex or on the block's plane:
    // the same material either way, which the boolean resolves in one or
    // the other. Where neither carries it the corner is cut as it stands.
    let carry = |model: &mut Model, from_vertex: bool| -> OgeomResult<Built> {
        let mut carried = Built::from_nothing(solid.clone());
        for k in 0..n {
            let planes = [m[ring[k]], m[ring[(k + 1) % n]]];
            if let Some(band) = band_gap(
                model,
                &carried.shape,
                corner,
                far,
                d,
                k,
                planes,
                radius,
                from_vertex,
                tol,
            )? {
                carried = Built {
                    shape: band.shape,
                    history: carried.history.then(&band.history),
                };
            }
        }
        Ok(carried)
    };

    // Which edge is first and which way the ring runs is the tool's
    // labelling ([`one_ball_tool`]), and the solid it builds is the same
    // for all 2N. The tool offers each labelling in turn and the first that
    // closes stands, so a corner no labelling closes is refused by name. A
    // failed attempt's nodes stay in the model unreferenced, under their
    // own operation.
    let attempt = |model: &mut Model, held: &Built, start: usize, reverse| -> OgeomResult<Built> {
        let tool = one_ball_tool(model, one, radius, start, reverse, tol)?;
        let rounded = ogeom_bool::cut(model, &held.shape, &tool.shape, tol)?;
        Ok(Built {
            shape: rounded.shape,
            history: held.history.then(&tool.history).then(&rounded.history),
        })
    };
    let mut outcome: Option<Built> = None;
    let mut last: Option<ogeom_core::OgeomError> = None;
    let mut tried: Vec<Shape> = Vec::new();
    'carries: for from_vertex in [Some(true), Some(false), None] {
        let carried = match from_vertex {
            Some(from_vertex) => match carry(model, from_vertex) {
                Ok(carried) => carried,
                Err(err @ (OgeomError::NotDone(_) | OgeomError::Construction(_))) => {
                    last = Some(err);
                    continue;
                }
                Err(err) => return Err(err),
            },
            None => Built::from_nothing(solid.clone()),
        };
        // A corner with no band to carry is cut once.
        if tried.iter().any(|held| held.is_same(&carried.shape)) {
            continue;
        }
        tried.push(carried.shape.clone());
        let labellings = (0..2 * n).map(|index| (index % n, index >= n));
        for (index, (start, reverse)) in labellings.enumerate() {
            if forced.is_some_and(|f| f != index) {
                continue;
            }
            match attempt(model, &carried, start, reverse) {
                Ok(built) => {
                    outcome = Some(built);
                    break 'carries;
                }
                Err(err) => last = Some(err),
            }
        }
    }
    let Some(rounded) = outcome else {
        ogeom_bail!(
            NotDone,
            "the corner tool's cut closed on none of the corner's {} \
             labellings; the last said: {}",
            2 * n,
            last.map_or_else(String::new, |e| e.to_string())
        );
    };
    let mut built = rounded;
    built.history.modify(vertex, built.shape.clone());
    Ok(built)
}

/// The corners at `vertices` of `solid` rounded in one cut, each as
/// [`round_vertex_with`] rounds it alone; `None` where that is not shown
/// to round them as rounding them one after another would, and the caller
/// rounds them so.
///
/// Every face of the solid must be planar, every corner one a single ball
/// touches every face of, and every corner's tool clear of every other
/// corner. With no curved face no band stands to be carried back to a
/// block. A tool lies within its block, and a corner whose block keeps a
/// margin from every other corner's leaves what that corner is read from
/// as it found it: the faces through its vertex keep their pieces round
/// it, the points its convexity was probed at stay as far inside or
/// outside, and no plane square to an edge of it passes through another
/// vertex to add a face there. A face of the solid on a plane through the
/// vertex that no face at the vertex lies on could be cut away whole; such
/// a corner is left to the caller. Each corner is then the same corner,
/// and its tool the same solid, before and after the others are cut, and
/// cutting the solid by every tool at once takes away the union of the
/// tools, as cutting by them in turn does. The margin is the chord the
/// convexity probes' classifier draws trims at, which also stands well
/// clear of the band in which a plane counts as passing through a vertex.
///
/// Cut in turn, a corner's patch meets the balls of the corners cut
/// before it, and a patch whose pole comes within the boolean's touching
/// distance of another ball is recharted, seams and all, which one cut of
/// a solid with no curved face never does. Every two balls therefore
/// stand further apart than that: no point of one comes that close to the
/// other.
///
/// The tool of each corner is built under the first labelling that builds
/// one; which labelling holds decides only the charts of the patch, as it
/// does when the corners are rounded in turn. Where the one cut does not
/// close, the tolerances it widened are put back.
pub(crate) fn round_vertices_together(
    model: &mut Model,
    solid: &Shape,
    vertices: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = crate::fillet::to_current_path("vertex:round_vertices_together");
    let mark = model.note_widened();
    match vertices_together(model, solid, vertices, radius, tol) {
        Ok(Some(built)) => {
            model.keep_widened(mark);
            Ok(Some(built))
        }
        Err(OgeomError::Cancelled) => {
            model.keep_widened(mark);
            Err(OgeomError::Cancelled)
        }
        Ok(None) | Err(_) => {
            model.undo_widened(mark);
            Ok(None)
        }
    }
}

/// [`round_vertices_together`], before what a cut that does not close
/// widened is put back.
fn vertices_together(
    model: &mut Model,
    solid: &Shape,
    vertices: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let margin = tol.confusion() * 1e4;
    let faces = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?;
    let planar = |face: &Shape| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .is_some_and(|data| {
                matches!(
                    model.geometry().surface(data.surface),
                    Some(ogeom_geom::SurfaceGeometry::Plane(_))
                )
            })
    };
    if !faces.iter().all(planar) {
        return Ok(None);
    }
    let mut corners: Vec<OneBall> = Vec::with_capacity(vertices.len());
    for vertex in vertices {
        let Corner::OneBall(one) = corner_at(model, solid, vertex, radius, tol)? else {
            return Ok(None);
        };
        // Every plane through the vertex is one a face at the vertex lies on.
        let mut at_vertex: Vec<Vector> = Vec::new();
        for face in &faces {
            let holds = ogeom_topo::explore_unique(model, face, ShapeType::Vertex)?
                .iter()
                .any(|v| v.is_same(vertex));
            if holds {
                at_vertex.push(ogeom_algo::face_normal(model, face, tol)?.1);
            }
        }
        let held = one.m.iter().all(|mk| {
            at_vertex
                .iter()
                .any(|n| n.cross(*mk).magnitude() < tol.angular() * 10.0)
        });
        if !held {
            return Ok(None);
        }
        corners.push(one);
    }
    let bounds: Vec<ogeom_math::Aabb> = corners.iter().map(|c| c.block_bounds(radius)).collect();
    for (i, one) in corners.iter().enumerate() {
        for (j, other) in corners.iter().enumerate() {
            if i == j {
                continue;
            }
            let apart = one.far.distance(other.far) - radius * 2.0;
            let clear = apart > radius * ogeom_bool::POLE_TOUCHING + margin
                && bounds[i].distance_to_box(&bounds[j]) > margin
                && other
                    .probes
                    .iter()
                    .all(|p| bounds[i].distance_to(*p) > margin)
                && one.d.iter().all(|dk| {
                    let along = dk.vector();
                    (other.corner - one.foot(along)).dot(along).abs() > margin
                });
            if !clear {
                return Ok(None);
            }
        }
    }
    let mut tools = Vec::with_capacity(corners.len());
    let mut tool_faces = Vec::with_capacity(corners.len());
    let mut history = ogeom_algo::History::new();
    for one in &corners {
        let n = one.m.len();
        let tool = (0..2 * n).find_map(|index| {
            match one_ball_tool(model, one, radius, index % n, index >= n, tol) {
                Ok(tool) => Some(Ok(tool)),
                Err(OgeomError::Cancelled) => Some(Err(OgeomError::Cancelled)),
                Err(_) => None,
            }
        });
        let Some(tool) = tool.transpose()? else {
            return Ok(None);
        };
        tool_faces.push(ogeom_topo::explore_unique(
            model,
            &tool.shape,
            ShapeType::Face,
        )?);
        history = history.then(&tool.history);
        tools.push(tool.shape);
    }
    let compound = model.add_compound(&tools)?;
    let rounded = ogeom_bool::cut(model, solid, &compound, tol)?;
    let mut history = history.then(&rounded.history);
    let patches = corner_faces(
        model,
        solid,
        &Built {
            shape: rounded.shape.clone(),
            history: history.clone(),
        },
    )?;
    for (vertex, faces) in vertices.iter().zip(tool_faces) {
        // Every tool owns its own corner, even when several are cut in one
        // boolean. A whole result solid cannot distinguish those patches.
        for patch in &patches {
            if faces
                .iter()
                .any(|face| history.trace(face).iter().any(|f| f.is_same(patch)))
            {
                history.generate(vertex, patch.clone());
            }
        }
        history.modify(vertex, rounded.shape.clone());
    }
    Ok(Some(Built {
        shape: rounded.shape,
        history,
    }))
}

/// The tool that rounds a one-ball corner: its block less its ball, under
/// the labelling whose first edge is `start` and whose ring runs backwards
/// where `reverse`.
///
/// The block: the corner bounded by its N host planes and, through the
/// ball's centre, the N planes square to its edges, where each band's
/// circle and the ball's own rim coincide, so the cut ends the band and
/// starts the patch on one curve. On a square corner it is the box of
/// side `radius`; on an oblique one a hexahedron; at a pyramid's apex a
/// polyhedron of eight faces. Its corners: the vertex, the foot of each
/// edge on its cutting plane, on each host plane the point where the
/// ball touches it, and the centre itself.
///
/// The ball on the corner's axes: a pole at one corner of the patch it
/// leaves and its seam meridian out past the block through the first
/// edge: the pole axis is the inward normal of the host plane holding
/// the first two edges, the third edge itself on a square corner. The
/// charts the block's faces and the ball wear differ between labellings
/// (where a rim is exact and where fitted, where a seam falls against a
/// patch arc); the solid does not.
fn one_ball_tool(
    model: &mut Model,
    one: &OneBall,
    radius: f64,
    start: usize,
    reverse: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let OneBall {
        corner,
        far,
        ref m,
        ref ring,
        ref d,
        ..
    } = *one;
    let n = m.len();
    let inward_of = |k: usize| m[ring[k % n]];
    // Edge t of the labelling and the host plane holding edges t and t+1.
    let edge_at = |t: usize| -> usize {
        if reverse {
            (start + n - t % n) % n
        } else {
            (start + t) % n
        }
    };
    let dl: Vec<Direction> = (0..n).map(|t| d[edge_at(t)]).collect();
    let ml: Vec<Vector> = (0..n)
        .map(|t| {
            // Forward: edges t, t+1 are E_s+t, E_s+t+1, both on plane
            // ring[s+t+1]. Reverse: E_s−t, E_s−t−1, both on ring[s−t].
            if reverse {
                inward_of((start + n - t % n) % n)
            } else {
                inward_of(edge_at(t) + 1)
            }
        })
        .collect();
    let along = |t: usize| -> Point { one.foot(dl[t].vector()) };
    let touch = |t: usize| -> Point { far - ml[t] * radius };
    if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
        eprintln!(
            "CORNER attempt start {start} reverse {reverse} far {far:?} d {:?} m {ml:?}",
            dl.iter().map(|x| x.vector()).collect::<Vec<_>>()
        );
    }
    model.begin_operation();
    let mut points: Vec<Point> = vec![corner];
    points.extend((0..n).map(along));
    points.extend((0..n).map(touch));
    points.push(far);
    let mut rings: Vec<Vec<usize>> = Vec::with_capacity(2 * n);
    for t in 0..n {
        // The host face on plane t: vertex, edge t's foot, the touch
        // point, edge t+1's foot.
        rings.push(vec![0, 1 + t, 1 + n + t, 1 + (t + 1) % n]);
        // The face square to edge t: its foot, the touch points either
        // side, the centre.
        rings.push(vec![1 + t, 1 + n + (t + n - 1) % n, 1 + 2 * n, 1 + n + t]);
    }
    let block = ogeom_algo::make_polyhedron(model, &points, &rings, tol)?.shape;
    let ball_frame = Frame::new(far, Direction::new(ml[0], tol)?, dl[0], tol)?;
    let ball = ogeom_algo::make_sphere(model, ball_frame, radius, tol)?.shape;
    let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
    if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
        let faces = ogeom_topo::explore_unique(model, &tool.shape, ShapeType::Face)?;
        let kinds: Vec<String> = faces
            .iter()
            .map(|f| {
                let kind = model
                    .node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .map_or("?", |sg| match sg {
                        ogeom_geom::SurfaceGeometry::Plane(_) => "plane",
                        ogeom_geom::SurfaceGeometry::Sphere(_) => "sphere",
                        _ => "other",
                    });
                let edges =
                    ogeom_topo::explore_unique(model, f, ShapeType::Edge).map_or(0, |e| e.len());
                format!("{kind}/{edges}")
            })
            .collect();
        eprintln!(
            "CORNER tool for start {start} reverse {reverse}: {} faces {kinds:?}",
            faces.len()
        );
    }
    Ok(tool)
}

/// The band along edge `k` of a one-ball corner carried back to the
/// corner's block, where a band another edge of the corner ran out left it
/// short.
///
/// Edge `k`'s band is the cylinder of `radius` about the line through the
/// ball's centre `far` along the edge's direction, between the two host
/// planes with inward normals `planes`; the block ends it on the plane
/// through `far` square to the edge. Rounded one at a time, another edge's
/// band may run out past the vertex and across this edge's crease, and
/// this band then stands flush against that one, further out than the
/// block reaches. In the section square to the edge, the material outside
/// the cylinder and within the two host planes is what no ball a radius in
/// from both planes reaches, so the band's wedge is cut back to the block,
/// `from_vertex` starting it at the vertex rather than on the block's
/// plane: behind that plane the block's cut takes the stretch anyway.
///
/// Only that band is carried: its near end stands against the band of
/// another edge of the corner (the crease's point there lies on that
/// band's cylinder), and no stretch of the crease between the vertex and
/// the band is still a sharp edge of the solid. A band ending against
/// anything else (a notch across the crease, a piece of the crease nobody
/// rounded) is left as it stands. `None` where no band is carried.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn band_gap(
    model: &mut Model,
    solid: &Shape,
    corner: Point,
    far: Point,
    axes: &[Direction],
    k: usize,
    planes: [Vector; 2],
    radius: f64,
    from_vertex: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let along = axes[k].vector();
    let on_line = tol.confusion() * 100.0;
    let point_of = |model: &Model, v: &Shape| -> OgeomResult<Option<Point>> {
        let Some(p) = model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
        else {
            return Ok(None);
        };
        Ok(Some(v.transform(model.datums())?.apply(p)))
    };
    // Where the band starts along the edge, from the block's plane: the
    // nearest of its vertices.
    let mut ends: Vec<Point> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        if !on_band(model, &face, far, along, radius, tol) {
            continue;
        }
        for v in ogeom_topo::explore_unique(model, &face, ShapeType::Vertex)? {
            if let Some(p) = point_of(model, &v)? {
                ends.push(p);
            }
        }
    }
    let Some(start) = ends
        .iter()
        .map(|p| (*p - far).dot(along))
        .min_by(f64::total_cmp)
    else {
        return Ok(None);
    };
    if start <= tol.confusion() {
        return Ok(None);
    }
    // The band's near end stands against another edge's band of this
    // corner: where its section there meets the crease, a point on that
    // band's cylinder.
    let from = (corner - far).dot(along);
    let gap = |p: Point, axis: Vector| {
        let w = p - far;
        ((w - axis * w.dot(axis)).magnitude() - radius).abs()
    };
    let crease_end = corner + along * (start - from);
    if !axes
        .iter()
        .enumerate()
        .any(|(j, axis)| j != k && gap(crease_end, axis.vector()) <= on_line)
    {
        return Ok(None);
    }
    // No stretch of the crease between the vertex and the band is still a
    // sharp edge of the solid. The crease leaves the other band at its rail,
    // where the band is tangent to the host plane the crease lies in, so
    // the two stand apart only quadratically there and the boolean places
    // the crossing within its vertices' tolerances. A piece each of whose
    // ends lies on one other band's cylinder within the end's own
    // tolerance lies on it within those tolerances all along (a line's
    // distance from a cylinder it stays outside is convex along it): it is
    // that band's boundary, not a stretch of sharp crease.
    for edge in ogeom_topo::explore_unique(model, solid, ShapeType::Edge)? {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &edge)? else {
            continue;
        };
        let (Some(pa), Some(pb)) = (point_of(model, &a)?, point_of(model, &b)?) else {
            continue;
        };
        let crease = |p: Point| {
            let w = p - corner;
            (w - along * w.dot(along)).magnitude() <= on_line
        };
        if !crease(pa) || !crease(pb) {
            continue;
        }
        let reach = |v: &Shape| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex())
                .map_or(0.0, |d| d.tolerance.get())
        };
        let (reach_a, reach_b) = (reach(&a), reach(&b));
        let on_band_seam = axes.iter().enumerate().any(|(j, axis)| {
            j != k && gap(pa, axis.vector()) <= reach_a && gap(pb, axis.vector()) <= reach_b
        });
        if on_band_seam {
            continue;
        }
        let (ta, tb) = ((pa - far).dot(along), (pb - far).dot(along));
        if ta.max(tb) > from + on_line && ta.min(tb) < start - on_line {
            return Ok(None);
        }
    }
    let face_on = |model: &Model, inward: Vector| -> OgeomResult<Option<Shape>> {
        for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
            let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                continue;
            };
            if !matches!(
                model.geometry().surface(data.surface),
                Some(ogeom_geom::SurfaceGeometry::Plane(_))
            ) {
                continue;
            }
            let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
            if (corner - origin).dot(outward).abs() <= on_line
                && outward.cross(inward).magnitude() < tol.angular() * 10.0
                && outward.dot(inward) < 0.0
            {
                return Ok(Some(face));
            }
        }
        Ok(None)
    };
    let (Some(first), Some(second)) = (face_on(model, planes[0])?, face_on(model, planes[1])?)
    else {
        return Ok(None);
    };
    let foot = corner + along * (far - corner).dot(along);
    let seat = crate::support::Seat {
        start: if from_vertex { corner } else { foot },
        end: foot + along * start,
        along,
        normals: [-planes[0], -planes[1]],
        faces: [first, second],
        convex: true,
    };
    let cut = crate::fillet::seated_fillet(model, solid, &seat, radius, None, &[], tol)?;
    // The carried band and the band it reaches lie on one cylinder: one
    // face, as the band the corner closes in one call is.
    let on_axis: Vec<Shape> = ogeom_topo::explore_unique(model, &cut.shape, ShapeType::Face)?
        .into_iter()
        .filter(|face| on_band(model, face, far, along, radius, tol))
        .collect();
    let (unified, _) = ogeom_heal::unify_same_domain_around(model, &cut.shape, &on_axis, tol)?;
    Ok(Some(Built {
        shape: unified.shape,
        history: cut.history.then(&unified.history),
    }))
}

/// Whether `face` lies on the cylinder of `radius` about the line through
/// `far` along `along`.
fn on_band(
    model: &Model,
    face: &Shape,
    far: Point,
    along: Vector,
    radius: f64,
    tol: Tolerances,
) -> bool {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return false;
    };
    let Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) = model.geometry().surface(data.surface)
    else {
        return false;
    };
    let frame = c.cylinder().frame();
    let axis = frame.z().vector();
    let off = far - frame.origin();
    (c.cylinder().radius() - radius).abs() <= tol.confusion()
        && axis.cross(along).magnitude() <= tol.angular()
        && (off - axis * off.dot(axis)).magnitude() <= tol.confusion()
}

/// A vertex of the region the ball's centre may occupy: a point a radius
/// in from three or more of the host planes and at least a radius from the
/// rest, with the planes it touches.
struct TipVertex {
    centre: Point,
    planes: Vec<usize>,
}

/// An edge of that region's tip between two of its vertices: the ball
/// rolling from one to the other touches two planes all the way, and
/// sweeps a cylinder.
struct Ridge {
    from: usize,
    to: usize,
    planes: [usize; 2],
}

/// Round a convex planar vertex whose faces no single ball touches.
///
/// The ball's centre may sit anywhere a radius in from every host plane:
/// a convex region whose tip, at a vertex one ball touches, is a single
/// point, and otherwise a few points joined by short edges: a rectangular
/// pyramid's apex has two, joined along the two long slopes. The rounded
/// corner is the envelope of every ball centred in that region: a sphere
/// at each tip vertex, a cylinder along each edge between them, and the
/// host planes themselves elsewhere. Exact, constant-radius, and what the
/// rolling ball leaves, where a plate would be fitted through the bands'
/// ends instead.
///
/// Each piece is cut with its own block: at a tip vertex, the corner
/// bounded by its three planes and the three planes through the centre
/// square to its edges, less the ball (the one-ball tool exactly) and
/// along an edge between two centres, the prism over the kite of the
/// virtual crease, the two touch points and the centre, between the two
/// planes square to the edge, less the cylinder. The compartments tile
/// the corner and meet on the planes square to the edges, where each
/// sphere's rim and the cylinder's end coincide, so the cuts consume one
/// another's flush faces cap to cap. An original edge leaves its tip
/// vertex along a ray, and the plane square to it there is where the
/// edge's band ends when the flush fillets follow.
fn setback_corner(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    corner: Point,
    m: &[Vector],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = m.len();
    let slack = tol.confusion() * 10.0;
    // The tip's vertices: every triple that meets at a point a radius in
    // from all the planes, merged where several triples name one point.
    let mut tips: Vec<TipVertex> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(centre) = ball_centre(corner, &[m[i], m[j], m[k]], radius) else {
                    continue;
                };
                let feasible = m.iter().all(|mk| mk.dot(centre - corner) >= radius - slack);
                if !feasible {
                    continue;
                }
                if tips.iter().any(|t| t.centre.distance(centre) <= slack) {
                    continue;
                }
                let planes: Vec<usize> = (0..n)
                    .filter(|&l| (m[l].dot(centre - corner) - radius).abs() <= slack)
                    .collect();
                tips.push(TipVertex { centre, planes });
            }
        }
    }
    if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
        for tip in &tips {
            eprintln!("TIP {:?} planes {:?}", tip.centre, tip.planes);
        }
    }
    if tips.is_empty() {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex hold no ball a radius in from \
             all of them"
        );
    }
    // Each tip vertex's edges: along every pair of its planes, the way the
    // rest of its planes allow; bounded where another plane becomes
    // tangent (a ridge to the tip vertex there) and a ray otherwise, the
    // way an original edge's band runs. A ray that lies along no original
    // edge, or a vertex with other than three edges, is a corner this tool
    // does not speak.
    let mut ridges: Vec<Ridge> = Vec::new();
    let mut rays: Vec<Vec<(Vector, [usize; 2])>> = vec![Vec::new(); tips.len()];
    for (index, tip) in tips.iter().enumerate() {
        for (a, b) in tip
            .planes
            .iter()
            .flat_map(|&a| tip.planes.iter().map(move |&b| (a, b)))
        {
            if a >= b {
                continue;
            }
            let cross = m[a].cross(m[b]);
            let length = cross.magnitude();
            if length <= tol.angular() * 10.0 {
                continue;
            }
            let Some(direction) = [1.0, -1.0]
                .into_iter()
                .map(|sign| cross * (sign / length))
                .find(|d| {
                    tip.planes
                        .iter()
                        .filter(|&&l| l != a && l != b)
                        .all(|&l| m[l].dot(*d) >= -tol.angular() * 10.0)
                })
            else {
                continue;
            };
            // The first other plane the ball meets rolling this way.
            let mut nearest: Option<(f64, usize)> = None;
            for (l, ml) in m.iter().enumerate() {
                if tip.planes.contains(&l) {
                    continue;
                }
                let rate = ml.dot(direction);
                if rate >= -tol.angular() {
                    continue;
                }
                let t = (radius - ml.dot(tip.centre - corner)) / rate;
                if t > slack && nearest.is_none_or(|(held, _)| t < held) {
                    nearest = Some((t, l));
                }
            }
            match nearest {
                Some((t, _)) => {
                    let end = tip.centre + direction * t;
                    let Some(to) = tips.iter().position(|o| o.centre.distance(end) <= slack) else {
                        ogeom_bail!(
                            Construction,
                            "the ball rolling between two of this vertex's faces meets a \
                             third where no tip vertex stands"
                        );
                    };
                    if index < to {
                        ridges.push(Ridge {
                            from: index,
                            to,
                            planes: [a, b],
                        });
                    }
                }
                None => rays[index].push((direction, [a, b])),
            }
        }
    }
    let concave = {
        let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
        let mut concave = false;
        for tip in &tips {
            if boundary.holds(model, tip.centre, tol)? != ogeom_algo::Containment::In {
                concave = true;
            }
        }
        concave
    };
    if concave {
        ogeom_bail!(
            Construction,
            "no material a radius in from every face at this vertex; a \
             concave vertex gains a ball instead of shedding one, and this \
             tool cannot round it"
        );
    }

    model.begin_operation();
    let mut rounded = solid.clone();
    let mut history = ogeom_algo::History::new();
    // The sphere at each tip vertex, with its own compartment: the wedge of
    // its planes, cut by the plane square to each of its edges through the
    // centre. A ridge's plane faces the corner, and the corner itself lies
    // in the ridge's compartment, not this one.
    for (index, tip) in tips.iter().enumerate() {
        // The vertex's edges, each on two of its planes: its rays and the
        // ridges that start or end here.
        let mut edges: Vec<(Vector, [usize; 2])> = rays[index].clone();
        for ridge in &ridges {
            let other = if ridge.from == index {
                ridge.to
            } else if ridge.to == index {
                ridge.from
            } else {
                continue;
            };
            let direction = (tips[other].centre - tip.centre).normalized(tol)?;
            edges.push((direction, ridge.planes));
        }
        let k = tip.planes.len();
        if edges.len() != k {
            ogeom_bail!(
                Construction,
                "a tip vertex of this corner touches {k} planes along {} edges, not \
                 {k}; the corner is not the simple convex one this tool speaks",
                edges.len()
            );
        }
        // In ring order, consecutive edges sharing a plane, so that host t
        // holds edges t and t+1: the labelling the ball's frame is read
        // off, its pole along a host and its seam out through an edge.
        let shared = |x: &[usize; 2], y: &[usize; 2]| -> Option<usize> {
            x.iter().copied().find(|p| y.contains(p))
        };
        let mut ring: Vec<(Vector, [usize; 2])> = vec![edges[0]];
        let mut used = vec![false; k];
        used[0] = true;
        while ring.len() < k {
            let last = ring[ring.len() - 1];
            let Some(next) = (0..k).find(|&i| {
                !used[i]
                    && shared(&edges[i].1, &last.1).is_some_and(|p| {
                        // The plane shared with the previous edge is not
                        // the one shared with the edge before that.
                        ring.len() < 2 || shared(&ring[ring.len() - 2].1, &last.1) != Some(p)
                    })
            }) else {
                ogeom_bail!(Construction, "the tip vertex's edges do not chain")
            };
            used[next] = true;
            ring.push(edges[next]);
        }
        let hosts: Vec<Vector> = (0..k)
            .map(|t| {
                shared(&ring[t].1, &ring[(t + 1) % k].1)
                    .map(|p| m[p])
                    .ok_or_else(|| {
                        ogeom_core::ogeom_err!(Construction, "the tip vertex's edges do not chain")
                    })
            })
            .collect::<OgeomResult<_>>()?;
        let directions: Vec<Vector> = ring.iter().map(|(d, _)| *d).collect();
        // Whether edge t is a ridge, whose rim plane a later cut's cap
        // stands in: a pole along a plane holding a ridge would put both
        // poles in that cap, and the cap's circle would have no chart
        // image. Labellings whose pole host holds no ridge go first.
        let is_ridge: Vec<bool> = ring
            .iter()
            .map(|(d, _)| !rays[index].iter().any(|(r, _)| r.dot(*d) > 1.0 - 1e-9))
            .collect();
        let mut walls: Vec<(Vector, Point)> = tip.planes.iter().map(|&p| (m[p], corner)).collect();
        for direction in &directions {
            walls.push((-*direction, tip.centre));
        }
        let tool = ball_block(
            model,
            &walls,
            corner,
            tip.centre,
            &hosts,
            &directions,
            &is_ridge,
            radius,
            tol,
        )?;
        let cut = ogeom_bool::cut(model, &rounded, &tool.shape, tol).map_err(|e| {
            ogeom_core::ogeom_err!(
                NotDone,
                "the cut by the ball's block at tip vertex {index} failed: {e}"
            )
        })?;
        history = history.then(&tool.history).then(&cut.history);
        rounded = cut.shape;
    }
    // The cylinder along each ridge: the flush fillet of a virtual crease (
    // the line the two planes it touches would meet along) between the
    // planes square to the ridge through its two centres, which are the
    // caps' own planes. The planar fillet builds that wedge face by face,
    // band, legs and caps, and melts it; the caps meet the spheres' rims
    // on the planes the vertex compartments already cut.
    for ridge in &ridges {
        let (v1, v2) = (tips[ridge.from].centre, tips[ridge.to].centre);
        let along = (v2 - v1).normalized(tol)?;
        let [a, c] = ridge.planes;
        let on_crease = |v: Point| corner + along * (v - corner).dot(along);
        let face_on = |model: &Model, inward: Vector| -> OgeomResult<Shape> {
            for face in ogeom_topo::explore_unique(model, &rounded, ShapeType::Face)? {
                let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                    continue;
                };
                if !matches!(
                    model.geometry().surface(data.surface),
                    Some(ogeom_geom::SurfaceGeometry::Plane(_))
                ) {
                    continue;
                }
                let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
                if (corner - origin).dot(outward).abs() <= tol.confusion() * 100.0
                    && outward.cross(inward).magnitude() < tol.angular() * 10.0
                    && outward.dot(inward) < 0.0
                {
                    return Ok(face);
                }
            }
            ogeom_bail!(
                Construction,
                "the ridge's host plane is no longer a face of the solid"
            )
        };
        let faces = [face_on(model, m[a])?, face_on(model, m[c])?];
        let seat = crate::support::Seat {
            start: on_crease(v1),
            end: on_crease(v2),
            along,
            normals: [-m[a], -m[c]],
            faces,
            convex: true,
        };
        let cut = crate::fillet::seated_fillet(model, &rounded, &seat, radius, None, &[], tol)
            .map_err(|e| ogeom_core::ogeom_err!(NotDone, "the ridge's flush fillet failed: {e}"))?;
        history = history.then(&cut.history);
        rounded = cut.shape;
    }
    history.modify(vertex, rounded.clone());
    Ok(Built {
        shape: rounded,
        history,
    })
}

use crate::engine::vertex_geometry::solve3;
pub(crate) use crate::engine::vertex_geometry::{convex_block, planes_meet};

/// The point a radius in from three planes through `corner`, or `None`
/// where they do not span.
fn ball_centre(corner: Point, m: &[Vector; 3], radius: f64) -> Option<Point> {
    solve3(m, [radius; 3]).map(|x| corner + x)
}

/// A compartment less the ball centred in it: the block from its walls,
/// cut by the sphere at `far`.
///
/// Host `t` holds edges `t` and `t + 1`. The ball's pole stands along a
/// host's normal and its seam meridian runs out through the first of that
/// host's edges (in a rim plane, so the seam doubles as a trim rather
/// than crossing a patch) and every labelling is offered in turn, the
/// first that closes standing. A cut that closes on a tool reaching past
/// its own block is a wrong tool, not a closed one, and is passed over
/// too.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn ball_block(
    model: &mut Model,
    walls: &[(Vector, Point)],
    corner: Point,
    far: Point,
    hosts: &[Vector],
    directions: &[Vector],
    is_ridge: &[bool],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = directions.len();
    let mut last: Option<ogeom_core::OgeomError> = None;
    // The ball's frames, best first. Along a ridge, the ridge's rim plane
    // is the sphere's equator, both poles stand outside the patch (the
    // patch lies on the corner's side of every rim plane, the poles on
    // the ridge's own axis either side of it) and a seam meridian turned
    // away from the corner never crosses it; the ridge fillet's cap, in
    // that same rim plane, then meets the sphere on a circle the chart
    // images exactly. Then the classic labellings: the pole along host t,
    // the seam out through edge t, hosts holding no ridge first.
    // The patch's corners are the touch points, a radius from the centre
    // against each host; a frame is judged by how far its seam meridian
    // keeps from every corner in longitude about the pole, and a pole the
    // patch itself contains (inside every rim plane's kept side) is no
    // frame at all, since a chart's degenerate point cannot stand inside
    // a face. Either sense of every ridge is offered, with the seam turned
    // from the corner or from the corners' mean, the clearest first.
    let touches: Vec<Vector> = hosts.iter().map(|m| -*m).collect();
    let inside_patch = |direction: Vector| -> bool {
        directions
            .iter()
            .all(|e| e.dot(direction) <= tol.angular() * 10.0)
    };
    let mut scored: Vec<(f64, Vector, Vector)> = Vec::new();
    for t in 0..n {
        if !is_ridge[t] {
            continue;
        }
        for pole in [directions[t], -directions[t]] {
            if inside_patch(pole) || inside_patch(-pole) {
                continue;
            }
            let flat = |v: Vector| v - pole * v.dot(pole);
            let mean = touches.iter().fold(Vector::ZERO, |acc, t| acc + *t);
            for seam in [flat(-mean), flat(far - corner)] {
                if seam.magnitude() <= tol.angular() * 10.0 {
                    continue;
                }
                let x = seam / seam.magnitude();
                let y = pole.cross(x);
                let clearance = touches
                    .iter()
                    .map(|t| t.dot(y).atan2(t.dot(x)).abs())
                    .fold(f64::INFINITY, f64::min);
                scored.push((clearance, pole, seam));
            }
        }
    }
    scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(core::cmp::Ordering::Equal));
    if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
        for (clearance, pole, seam) in &scored {
            eprintln!("FRAME clearance {clearance:.3} pole {pole:?} seam {seam:?}");
        }
    }
    let mut frames: Vec<(Vector, Vector)> = scored.into_iter().map(|(_, p, s)| (p, s)).collect();
    // Host t holds edges t and t+1; a labelling's pole host is `start`
    // forward and `start − 1` reversed.
    let pole_host = |index: usize| -> usize {
        let (start, reverse) = (index % n, index >= n);
        if reverse {
            (start + 2 * n - 1) % n
        } else {
            start
        }
    };
    let clean = |index: usize| -> bool {
        let h = pole_host(index);
        !is_ridge[h] && !is_ridge[(h + 1) % n]
    };
    let mut order: Vec<usize> = (0..2 * n).collect();
    order.sort_by_key(|&index| !clean(index));
    frames.extend(
        order
            .iter()
            .map(|&index| (hosts[pole_host(index)], directions[index % n])),
    );
    for (index, (pole, seam)) in frames.into_iter().enumerate() {
        model.begin_operation();
        let attempt = (|| -> OgeomResult<Built> {
            let block = convex_block(model, walls, tol)?;
            let block_bound = ogeom_algo::shape_bounds(model, &block, tol)?;
            let ball_frame = Frame::new(
                far,
                Direction::new(pole, tol)?,
                Direction::new(seam, tol)?,
                tol,
            )?;
            let ball = ogeom_algo::make_sphere(model, ball_frame, radius, tol)?.shape;
            let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
            let tool_bound = ogeom_algo::shape_bounds(model, &tool.shape, tol)?;
            let reach = tol.confusion() * 1e3;
            let (Some(block_lo), Some(block_hi), Some(tool_lo), Some(tool_hi)) = (
                block_bound.low(),
                block_bound.high(),
                tool_bound.low(),
                tool_bound.high(),
            ) else {
                ogeom_bail!(NotDone, "the ball's cut left no tool");
            };
            if tool_lo.x < block_lo.x - reach
                || tool_lo.y < block_lo.y - reach
                || tool_lo.z < block_lo.z - reach
                || tool_hi.x > block_hi.x + reach
                || tool_hi.y > block_hi.y + reach
                || tool_hi.z > block_hi.z + reach
            {
                ogeom_bail!(
                    NotDone,
                    "the ball's cut left a tool reaching past its own block"
                );
            }
            Ok(tool)
        })();
        match attempt {
            Ok(tool) => {
                if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
                    eprintln!("BALL frame {index} closed");
                }
                return Ok(tool);
            }
            Err(err) => {
                if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
                    eprintln!("BALL frame {index} pole {pole:?} seam {seam:?}: {err}");
                }
                last = Some(err);
            }
        }
    }
    ogeom_bail!(
        NotDone,
        "the corner tool's block closed on none of its {} frames; the last said: {}",
        2 * n,
        last.map_or_else(String::new, |e| e.to_string())
    )
}
