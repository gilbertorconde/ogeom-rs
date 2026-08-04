//! Boolean operations.
//!
//! *Elsewhere:* `BOPDS`, `BOPAlgo`, `BOPTools`, `IntTools` and `BRepAlgoAPI`.
//!
//! The architectural insight worth preserving: **fuse, common, cut, section and
//! split are one algorithm — general fuse — plus a selection predicate.** The
//! pipeline is
//!
//! 1. a data structure of indexed sub-shapes, per-type interference lists
//!    (V/V, V/E, V/F, E/E, E/F, F/F), pave blocks and common blocks;
//! 2. a pave filler running in strictly increasing dimension, ending in face/face
//!    intersection, section-edge construction and pcurve generation;
//! 3. a builder that splits faces in 2D parametric space, unifies same-domain
//!    faces, rebuilds solids from face sets, and repairs tolerances;
//! 4. filters that select from the general-fuse result.
//!
//! # The vertical slice
//!
//! What is implemented is that pipeline end to end for **planar-faced solids
//! with straight edges** — boxes, wedges, prisms, and anything else whose
//! boundary is polygons. Within that world nothing is approximated: sections
//! are exact plane/plane lines clipped to both faces' trims, faces are split
//! in their own parameter planes by an exact arrangement, each piece is
//! classified against the other solid by the exact ray classifier, and the
//! result is sewn, checked closed, and nested into solids and voids.
//!
//! A curved face, a curved edge, or a tangential same-domain contact is
//! *refused* with an error saying so, never silently mishandled — the refusal
//! is recorded in `docs/SCOPE.md`'s deferred table, and the machinery §7
//! built (marched sections with pcurves both sides) is what lifts it when the
//! slice widens.
//!
//! The pave information of stage 2 lives inside the per-face arrangements
//! here: a section's endpoint on an edge is the same world point in both
//! adjacent faces' arrangements, which is what makes the split edges sew back
//! into shared topology. The indexed, cross-referenced data structure arrives
//! when curved intersections force paves to live on curves rather than in
//! planes.

mod arrange;

use og_algo::{
    Built, Containment, History, classify_in_solid_exact, is_shell_closed, make_face, make_polygon,
    sew, shape_bounds,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d as _;
use og_geom::{Curve, PlaneSurface, SurfaceGeometry};
use og_math::{Direction, Frame, Plane, Point, Point2, Vector};
use og_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

use arrange::{Piece, inside_rings, split};

// --- gathering ---------------------------------------------------------------

/// One face of an argument, flattened into its own plane.
struct PlanarFace {
    face: Shape,
    /// The world plane, oriented so its normal points out of the solid.
    plane: Plane,
    /// The 2D frame the arrangement runs in: `z` is the outward normal, so
    /// counter-clockwise rings enclose material seen from outside.
    frame: Frame,
    /// Boundary rings as exact corner points in `frame` coordinates.
    rings: Vec<Vec<Point2>>,
}

/// An argument solid, gathered and vetted for the slice.
struct PlanarSolid {
    solid: Shape,
    faces: Vec<PlanarFace>,
}

fn gather(model: &Model, solid: &Shape, tol: Tolerances) -> OgResult<PlanarSolid> {
    if model.kind_of(solid)? != ShapeType::Solid {
        og_bail!(Construction, "boolean arguments are solids");
    }
    for shell in explore_unique(model, solid, ShapeType::Shell)? {
        if !is_shell_closed(model, &shell)? {
            og_bail!(Construction, "an open shell bounds no volume to operate on");
        }
    }

    let mut faces = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        let SurfaceGeometry::Plane(stored) = surface else {
            og_bail!(
                NotDone,
                "the boolean's vertical slice covers planar-faced solids; a \
                 curved face is refused rather than mishandled — see the \
                 deferred table in docs/SCOPE.md"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut plane = stored.plane().transformed(&placement, tol)?;
        if face.orientation() == og_topo::Orientation::Reversed {
            plane = plane.reversed();
        }
        let axis = Direction::from_cross(
            plane.normal().vector(),
            Vector::new(0.312_8, 0.549_1, 0.775_6),
            tol,
        )?;
        let frame = Frame::new(plane.origin(), plane.normal(), axis, tol)?;

        let mut rings = Vec::new();
        for wire in explore(model, &face, Filter::OfType(ShapeType::Wire))? {
            let ring = ring_corners(model, &wire, &frame, tol)?;
            if ring.len() >= 3 {
                rings.push(ring);
            }
        }
        if rings.is_empty() {
            og_bail!(Construction, "a face with no boundary ring bounds nothing");
        }
        faces.push(PlanarFace {
            face,
            plane,
            frame,
            rings,
        });
    }
    if faces.is_empty() {
        og_bail!(Construction, "a solid with no faces bounds nothing");
    }
    Ok(PlanarSolid {
        solid: solid.clone(),
        faces,
    })
}

/// A wire's corner points, chained in walking order, in `frame` coordinates.
fn ring_corners(
    model: &Model,
    wire: &Shape,
    frame: &Frame,
    tol: Tolerances,
) -> OgResult<Vec<Point2>> {
    // Each straight edge contributes its two endpoints; the chain is
    // recovered by connectivity, which for a closed wire of segments is
    // unambiguous.
    let mut spans: Vec<(Point, Point)> = Vec::new();
    for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
        let Some(node) = model.node(&edge) else {
            og_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            og_bail!(Construction, "edge node holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            og_bail!(Dangling, "curve is not in this model");
        };
        if !matches!(geometry, Curve::Line(_)) {
            og_bail!(
                NotDone,
                "the boolean's vertical slice covers straight edges; a curved \
                 edge is refused rather than mishandled — see the deferred \
                 table in docs/SCOPE.md"
            );
        }
        let placement = edge.transform(model.datums())?;
        let a = placement.apply(geometry.point_at(range.0, tol)?);
        let b = placement.apply(geometry.point_at(range.1, tol)?);
        spans.push((a, b));
    }
    if spans.is_empty() {
        return Ok(Vec::new());
    }

    let reach = tol.confusion() * 10.0;
    let mut chain: Vec<Point> = vec![spans[0].0];
    let mut tail = spans[0].1;
    chain.push(tail);
    let mut used = vec![false; spans.len()];
    used[0] = true;
    for _ in 1..spans.len() {
        let Some((index, flipped)) = spans.iter().enumerate().find_map(|(i, (a, b))| {
            if used[i] {
                return None;
            }
            if a.distance(tail) <= reach {
                return Some((i, false));
            }
            if b.distance(tail) <= reach {
                return Some((i, true));
            }
            None
        }) else {
            og_bail!(Construction, "a face boundary wire does not chain closed");
        };
        used[index] = true;
        let (a, b) = spans[index];
        tail = if flipped { a } else { b };
        chain.push(tail);
    }
    if chain
        .last()
        .is_none_or(|end| end.distance(chain[0]) > reach)
    {
        og_bail!(Construction, "a face boundary wire does not chain closed");
    }
    chain.pop();

    let mut out = Vec::with_capacity(chain.len());
    for p in chain {
        let local = frame.to_local(p);
        if local.z.abs() > reach {
            og_bail!(
                Construction,
                "a boundary corner sits {} off its face's plane",
                local.z.abs()
            );
        }
        out.push(Point2::new(local.x, local.y));
    }
    Ok(out)
}

// --- the filler: sections ----------------------------------------------------

/// One straight section segment, in world space.
#[derive(Debug, Clone, Copy)]
struct Section {
    a: Point,
    b: Point,
}

/// The line two planes share, if they are not parallel.
fn planes_line(a: &Plane, b: &Plane) -> Option<(Point, Vector)> {
    let na = a.normal().vector();
    let nb = b.normal().vector();
    let direction = na.cross(nb);
    let magnitude = direction.magnitude();
    if magnitude <= 1e-9 {
        return None;
    }
    let direction = direction / magnitude;
    // A point on both planes, found in the span of the two normals.
    let da = na.dot(a.origin().to_vector());
    let db = nb.dot(b.origin().to_vector());
    let aa = na.dot(na);
    let ab = na.dot(nb);
    let bb = nb.dot(nb);
    let determinant = ab.mul_add(-ab, aa * bb);
    let ka = db.mul_add(-ab, da * bb) / determinant;
    let kb = da.mul_add(-ab, db * aa) / determinant;
    let point = Point::ORIGIN + na * ka + nb * kb;
    Some((point, direction))
}

/// The parameter intervals of the line lying inside a face's trim.
fn line_inside_face(
    face: &PlanarFace,
    origin: Point,
    direction: Vector,
    tol: Tolerances,
) -> Vec<(f64, f64)> {
    let snap = tol.confusion();
    let l0 = face.frame.to_local(origin);
    let l1 = face.frame.to_local(origin + direction);
    let p0 = Point2::new(l0.x, l0.y);
    let e = Point2::new(l1.x, l1.y) - p0;

    let mut crossings: Vec<f64> = Vec::new();
    for ring in &face.rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let d = b - a;
            let denominator = e.cross(d);
            if denominator.abs() <= snap {
                continue;
            }
            let s = (a - p0).cross(e) / denominator;
            if !(-snap..=1.0 + snap).contains(&s) {
                continue;
            }
            crossings.push((a - p0).cross(d) / denominator);
        }
    }
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (*a - *b).abs() <= snap);

    let mut intervals = Vec::new();
    for pair in crossings.windows(2) {
        let mid = f64::midpoint(pair[0], pair[1]);
        // The parity of a sorted crossing list can lie at a grazed corner;
        // asking the rings directly cannot.
        if inside_rings(&face.rings, p0 + e * mid) {
            intervals.push((pair[0], pair[1]));
        }
    }
    intervals
}

/// Sorted-interval intersection.
fn overlap(a: &[(f64, f64)], b: &[(f64, f64)], least: f64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for &(a0, a1) in a {
        for &(b0, b1) in b {
            let lo = a0.max(b0);
            let hi = a1.min(b1);
            if hi - lo > least {
                out.push((lo, hi));
            }
        }
    }
    out
}

/// Split every section segment wherever another touches or crosses it, so
/// both faces of every pair see the same subdivision — the slice's version of
/// globally consistent paves.
fn make_consistent(sections: &[Section], tol: Tolerances) -> Vec<Section> {
    let snap = tol.confusion();
    let mut cuts: Vec<Vec<f64>> = sections.iter().map(|_| vec![0.0, 1.0]).collect();
    for i in 0..sections.len() {
        for j in 0..sections.len() {
            if i == j {
                continue;
            }
            let s = &sections[i];
            let o = &sections[j];
            let d = s.b - s.a;
            let length = d.magnitude();
            // Endpoints of the other lying on this segment.
            for p in [o.a, o.b] {
                let t = (p - s.a).dot(d) / (length * length);
                if t > snap / length && t < 1.0 - snap / length && (s.a + d * t).distance(p) <= snap
                {
                    cuts[i].push(t);
                }
            }
            // A transversal crossing in space.
            let e = o.b - o.a;
            let cross = d.cross(e);
            let denominator = cross.dot(cross);
            if denominator > snap * snap {
                let between = o.a - s.a;
                let t = between.cross(e).dot(cross) / denominator;
                let u = between.cross(d).dot(cross) / denominator;
                if (0.0..=1.0).contains(&t)
                    && (0.0..=1.0).contains(&u)
                    && (s.a + d * t).distance(o.a + e * u) <= snap
                {
                    cuts[i].push(t);
                }
            }
        }
    }
    let mut out = Vec::new();
    for (section, mut ts) in sections.iter().zip(cuts) {
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        let d = section.b - section.a;
        let length = d.magnitude();
        let mut previous = None;
        for t in ts {
            if previous.is_none_or(|p: f64| (t - p) * length > snap) {
                if let Some(p) = previous {
                    out.push(Section {
                        a: section.a + d * p,
                        b: section.a + d * t,
                    });
                }
                previous = Some(t);
            }
        }
    }
    out
}

// --- the general fuse --------------------------------------------------------

/// One piece of one argument face, classified against the other argument.
struct FacePiece {
    source: Shape,
    /// Rings in world space: outer counter-clockwise around `outward`, holes
    /// clockwise.
    rings: Vec<Vec<Point>>,
    outward: Direction,
    state: Containment,
}

/// The full result the filters select from.
struct GeneralFused {
    pieces_a: Vec<FacePiece>,
    pieces_b: Vec<FacePiece>,
    sections: Vec<Section>,
}

fn general_fuse(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<GeneralFused> {
    let ga = gather(model, a, tol)?;
    let gb = gather(model, b, tol)?;
    let snap = tol.confusion();

    // Face/face sections: the line of the two planes, clipped to both trims.
    let mut raw: Vec<Section> = Vec::new();
    for fa in &ga.faces {
        for fb in &gb.faces {
            let Some((origin, direction)) = planes_line(&fa.plane, &fb.plane) else {
                continue;
            };
            let ia = line_inside_face(fa, origin, direction, tol);
            if ia.is_empty() {
                continue;
            }
            let ib = line_inside_face(fb, origin, direction, tol);
            for (lo, hi) in overlap(&ia, &ib, snap) {
                raw.push(Section {
                    a: origin + direction * lo,
                    b: origin + direction * hi,
                });
            }
        }
    }
    let sections = make_consistent(&raw, tol);

    let pieces_a = split_side(model, &ga, &gb.solid, &sections, tol)?;
    let pieces_b = split_side(model, &gb, &ga.solid, &sections, tol)?;
    Ok(GeneralFused {
        pieces_a,
        pieces_b,
        sections,
    })
}

/// Split every face of one argument by the sections lying in it, and classify
/// each piece against the other argument.
fn split_side(
    model: &Model,
    own: &PlanarSolid,
    other: &Shape,
    sections: &[Section],
    tol: Tolerances,
) -> OgResult<Vec<FacePiece>> {
    let snap = tol.confusion();
    let mut out = Vec::new();
    for face in &own.faces {
        let mut in_plane: Vec<(Point2, Point2)> = Vec::new();
        for section in sections {
            if face.plane.distance_to(section.a) > snap || face.plane.distance_to(section.b) > snap
            {
                continue;
            }
            let la = face.frame.to_local(section.a);
            let lb = face.frame.to_local(section.b);
            in_plane.push((Point2::new(la.x, la.y), Point2::new(lb.x, lb.y)));
        }
        let pieces = split(&face.rings, &in_plane, snap)?;
        for Piece { rings, interior } in pieces {
            let world = |p: Point2| face.frame.to_world(Point::new(p.x, p.y, 0.0));
            let probe = world(interior);
            let state = classify_in_solid_exact(model, other, probe, tol)?;
            if state == Containment::On {
                og_bail!(
                    NotDone,
                    "a piece of a face lies on the other solid's boundary; \
                     same-domain contact is refused by the vertical slice \
                     rather than resolved — see the deferred table in \
                     docs/SCOPE.md"
                );
            }
            out.push(FacePiece {
                source: face.face.clone(),
                rings: rings
                    .into_iter()
                    .map(|ring| ring.into_iter().map(world).collect())
                    .collect(),
                outward: face.plane.normal(),
                state,
            });
        }
    }
    Ok(out)
}

// --- rebuilding --------------------------------------------------------------

/// Build one face from a piece, forward or with its material side flipped.
fn build_piece(
    model: &mut Model,
    piece: &FacePiece,
    reversed: bool,
    tol: Tolerances,
) -> OgResult<Shape> {
    let normal = if reversed {
        piece.outward.reversed()
    } else {
        piece.outward
    };
    let origin = piece.rings[0][0];
    let plane = Plane::through(origin, normal);

    // The surface's parameter window, from the rings themselves.
    let frame = plane.frame();
    let mut u = (f64::INFINITY, f64::NEG_INFINITY);
    let mut v = (f64::INFINITY, f64::NEG_INFINITY);
    for p in piece.rings.iter().flatten() {
        let local = frame.to_local(*p);
        u = (u.0.min(local.x), u.1.max(local.x));
        v = (v.0.min(local.y), v.1.max(local.y));
    }
    let surface = PlaneSurface::over(plane, (u.0 - 1.0, u.1 + 1.0), (v.0 - 1.0, v.1 + 1.0))?;

    let mut wires = Vec::new();
    for ring in &piece.rings {
        let points: Vec<Point> = if reversed {
            ring.iter().rev().copied().collect()
        } else {
            ring.clone()
        };
        wires.push(make_polygon(model, &points, true, tol)?.shape);
    }
    Ok(make_face(model, surface.into(), &wires, tol)?.shape)
}

/// Attach same-parameter pcurves to every edge of a sewn face.
///
/// Deliberately *after* sewing, against each surviving edge's own canonical
/// curve: sewing may keep either of two coincident edges and flip the other's
/// uses, and a pcurve computed before that decision can end up parameterized
/// backwards relative to the survivor. Computed from the survivor's endpoints
/// in the face's own plane, the pcurve agrees with the curve it annotates by
/// construction.
fn attach_pcurves(model: &mut Model, face: &Shape, tol: Tolerances) -> OgResult<()> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        og_bail!(Construction, "face node holds no face data");
    };
    let surface_id = data.surface;
    let Some(SurfaceGeometry::Plane(stored)) = model.geometry().surface(surface_id) else {
        og_bail!(Construction, "a sewn boolean face is not planar");
    };
    let frame = stored.plane().frame();
    let mut wanted: Vec<(Shape, Point2, Point2)> = Vec::new();
    for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
        let Some(node) = model.node(&edge) else {
            og_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(edge_data) = node.data() else {
            og_bail!(Construction, "edge node holds no edge data");
        };
        if edge_data
            .representations
            .iter()
            .any(|r| matches!(r, EdgeRepr::PCurve { surface, .. } if *surface == surface_id))
        {
            continue;
        }
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            og_bail!(Dangling, "curve is not in this model");
        };
        let from = geometry.point_at(range.0, tol)?;
        let to = geometry.point_at(range.1, tol)?;
        let local = |p: Point| {
            let l = frame.to_local(p);
            Point2::new(l.x, l.y)
        };
        wanted.push((edge, local(from), local(to)));
    }
    for (edge, a2, b2) in wanted {
        let pcurve: og_geom::PlanarCurve = og_geom::Line2d::segment(a2, b2, tol)?.into();
        og_algo::attach_pcurve(
            model,
            &edge,
            pcurve,
            surface_id,
            og_topo::Location::identity(),
            (0.0, a2.distance(b2)),
        )?;
    }
    Ok(())
}

/// Sew kept pieces, demand closure, and nest shells into solids and voids.
fn assemble(
    model: &mut Model,
    kept: &[(usize, bool)],
    pieces: &[&FacePiece],
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgResult<Built> {
    let mut history = History::new();

    if kept.is_empty() {
        // A legitimate answer: cutting a solid away entirely leaves nothing.
        let empty = model.add_compound(&[])?;
        for piece in pieces {
            history.delete(&piece.source);
        }
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }

    let mut faces = Vec::new();
    let mut kept_sources: Vec<Shape> = Vec::new();
    for &(index, reversed) in kept {
        let piece = pieces[index];
        let built = build_piece(model, piece, reversed, tol)?;
        history.modify(&piece.source, built.clone());
        kept_sources.push(piece.source.clone());
        faces.push(built);
    }
    for piece in pieces {
        if !kept_sources.iter().any(|s| s == &piece.source) {
            history.delete(&piece.source);
        }
    }

    let sewn = sew(model, &faces, tol)?;
    for shell in &sewn.shells {
        for face in explore(model, shell, Filter::OfType(ShapeType::Face))? {
            attach_pcurves(model, &face, tol)?;
        }
        if !is_shell_closed(model, shell)? {
            og_bail!(
                NotDone,
                "the kept pieces did not close into a shell; the configuration \
                 is beyond the vertical slice"
            );
        }
    }

    // Nest: a shell whose bound sits inside another's is that solid's void.
    let mut bounds = Vec::new();
    for shell in &sewn.shells {
        bounds.push(shape_bounds(model, shell, tol)?);
    }
    let mut solids = Vec::new();
    for (i, shell) in sewn.shells.iter().enumerate() {
        let contained = bounds
            .iter()
            .enumerate()
            .any(|(j, other)| j != i && other.contains_box(&bounds[i]));
        if contained {
            continue;
        }
        let mut group = vec![shell.clone()];
        for (j, candidate) in sewn.shells.iter().enumerate() {
            if j != i && bounds[i].contains_box(&bounds[j]) {
                group.push(candidate.clone());
            }
        }
        solids.push(model.add_solid(&group)?);
    }
    let result = if solids.len() == 1 {
        solids.remove(0)
    } else {
        model.add_compound(&solids)?
    };
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

// --- the operations ----------------------------------------------------------

/// The union of two solids.
///
/// # Errors
///
/// [`OgError::NotDone`](og_core::OgError::NotDone) for configurations beyond
/// the vertical slice — curved faces, tangential same-domain contact;
/// [`OgError::Construction`](og_core::OgError::Construction) for arguments
/// that are not closed solids.
pub fn fuse(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let pieces: Vec<&FacePiece> = fused.pieces_a.iter().chain(&fused.pieces_b).collect();
    let kept: Vec<(usize, bool)> = pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.state == Containment::Out)
        .map(|(i, _)| (i, false))
        .collect();
    assemble(model, &kept, &pieces, a, b, tol)
}

/// The intersection of two solids.
///
/// # Errors
///
/// As [`fuse`].
pub fn common(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let pieces: Vec<&FacePiece> = fused.pieces_a.iter().chain(&fused.pieces_b).collect();
    let kept: Vec<(usize, bool)> = pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.state == Containment::In)
        .map(|(i, _)| (i, false))
        .collect();
    assemble(model, &kept, &pieces, a, b, tol)
}

/// The first solid with the second removed.
///
/// The pieces of `b` that close the cut into `a`'s material bound the removed
/// volume from `b`'s side, so they join the result with their material side
/// flipped.
///
/// # Errors
///
/// As [`fuse`].
pub fn cut(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let boundary_a = fused.pieces_a.len();
    let pieces: Vec<&FacePiece> = fused.pieces_a.iter().chain(&fused.pieces_b).collect();
    let kept: Vec<(usize, bool)> = pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            let from_a = i < boundary_a;
            match (from_a, p.state) {
                (true, Containment::Out) => Some((i, false)),
                (false, Containment::In) => Some((i, true)),
                _ => None,
            }
        })
        .collect();
    assemble(model, &kept, &pieces, a, b, tol)
}

/// The edges where the two solids' boundaries cross.
///
/// # Errors
///
/// As [`fuse`].
pub fn section(model: &mut Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgResult<Built> {
    let fused = general_fuse(model, a, b, tol)?;
    let mut edges = Vec::new();
    let mut history = History::new();
    for s in &fused.sections {
        let wire = make_polygon(model, &[s.a, s.b], false, tol)?.shape;
        for edge in explore(model, &wire, Filter::OfType(ShapeType::Edge))? {
            edges.push(edge);
        }
    }
    let result = model.add_compound(&edges)?;
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_algo::{check, make_box, make_cylinder, volume_properties};
    use og_mesh::Deflection;

    const T: Tolerances = Tolerances::millimetres();

    fn boxes(model: &mut Model) -> (Shape, Shape) {
        // Overlapping in the corner cube [1,2]^3.
        let a = make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let at = Frame::new(Point::new(1.0, 1.0, 1.0), Direction::Z, Direction::X, T).unwrap();
        let b = make_box(model, at, (2.0, 2.0, 2.0), T).unwrap();
        (a.shape, b.shape)
    }

    fn volume(model: &Model, shape: &Shape) -> f64 {
        let fine = Deflection {
            chord: 1e-3,
            ..Deflection::default()
        };
        volume_properties(model, shape, fine, T).unwrap().mass
    }

    fn assert_valid(model: &Model, shape: &Shape) {
        let diagnosis = check(model, shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "the result fails validity: {:?}",
            diagnosis.problems
        );
    }

    #[test]
    fn fuse_of_overlapping_boxes_has_the_inclusion_exclusion_volume() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let fused = fuse(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &fused.shape);
        // 8 + 8 - 1, and the tessellation of planar faces is exact.
        assert!((volume(&model, &fused.shape) - 15.0).abs() < 1e-9);
        assert_eq!(
            fused.history.modified(&a),
            std::slice::from_ref(&fused.shape)
        );
    }

    #[test]
    fn common_of_overlapping_boxes_is_the_overlap_cube() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = common(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn cut_removes_the_overlap_from_the_first_argument() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = cut(&mut model, &a, &b, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 7.0).abs() < 1e-9);
    }

    #[test]
    fn section_of_overlapping_boxes_is_the_six_segment_seam() {
        // The boundaries cross along a closed loop of six unit segments
        // around the shared corner cube.
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = section(&mut model, &a, &b, T).unwrap();
        let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
        assert_eq!(edges.len(), 6);
    }

    #[test]
    fn fuse_of_disjoint_boxes_is_a_compound_of_both() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let far = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let b = make_box(&mut model, far, (2.0, 2.0, 2.0), T).unwrap();
        let fused = fuse(&mut model, &a.shape, &b.shape, T).unwrap();
        assert_eq!(
            model.kind_of(&fused.shape).unwrap(),
            ShapeType::Compound,
            "disjoint solids stay two solids, honestly"
        );
        assert!((volume(&model, &fused.shape) - 16.0).abs() < 1e-9);
    }

    #[test]
    fn cutting_a_through_post_leaves_a_slab_with_a_hole() {
        // The result's top and bottom faces have holes: the arrangement's
        // nesting, the sewing, and the checker all have to cooperate.
        let mut model = Model::new();
        let slab = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let below = Frame::new(Point::new(1.5, 1.5, -1.0), Direction::Z, Direction::X, T).unwrap();
        let post = make_box(&mut model, below, (1.0, 1.0, 3.0), T).unwrap();
        let result = cut(&mut model, &slab.shape, &post.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 15.0).abs() < 1e-9);
    }

    #[test]
    fn common_with_a_contained_box_is_that_box() {
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let at = Frame::new(Point::new(2.0, 2.0, 2.0), Direction::Z, Direction::X, T).unwrap();
        let inner = make_box(&mut model, at, (2.0, 2.0, 2.0), T).unwrap();
        let result = common(&mut model, &outer.shape, &inner.shape, T).unwrap();
        assert_valid(&model, &result.shape);
        assert!((volume(&model, &result.shape) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn cutting_everything_away_leaves_an_empty_compound() {
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let at = Frame::new(Point::new(2.0, 2.0, 2.0), Direction::Z, Direction::X, T).unwrap();
        let inner = make_box(&mut model, at, (2.0, 2.0, 2.0), T).unwrap();
        let result = cut(&mut model, &inner.shape, &outer.shape, T).unwrap();
        assert_eq!(model.kind_of(&result.shape).unwrap(), ShapeType::Compound);
        assert!(
            explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_curved_argument_is_refused_with_instructions() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let drum = make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
        assert!(fuse(&mut model, &block.shape, &drum.shape, T).is_err());
    }

    #[test]
    fn every_source_face_is_accounted_for_in_the_history() {
        let mut model = Model::new();
        let (a, b) = boxes(&mut model);
        let result = cut(&mut model, &a, &b, T).unwrap();
        for face in explore(&model, &a, Filter::OfType(ShapeType::Face)).unwrap() {
            assert!(
                result.history.is_affected(&face),
                "a face of the first argument vanished from the history"
            );
        }
    }
}
