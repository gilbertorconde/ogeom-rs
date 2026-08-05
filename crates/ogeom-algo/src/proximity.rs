//! Minimum distance between shapes.
//!
//! *Elsewhere* this is `BRepExtrema_DistShapeShape`. The geometry-level
//! extrema in `ogeom-intersect` answer where two curves or surfaces come
//! nearest; this module assembles those answers for topology, where a shape
//! is vertices, edges and faces and the nearest approach may land on any of
//! them.
//!
//! # The assembly argument
//!
//! The nearest distance between two shapes is attained either at an interior
//! stationary approach of a pair of elements, or on some element's boundary —
//! and an element's boundary is itself an element: a face's boundary is its
//! edges, an edge's boundary is its vertices. So walking every pair of
//! elements — vertex against vertex, edge and face; edge against edge and
//! face; face against face — with stationary approaches for the interiors and
//! projections for the points covers every candidate, and the geometry level
//! is allowed to answer "no interior approach" honestly because the pair that
//! owns the boundary case is in the same sweep.
//!
//! A face's interior approach is accepted only where its foot lands inside
//! the face's trimming; a foot outside or too near the boundary is dropped,
//! because the true nearest point of that configuration is on an edge and the
//! edge pairs find it exactly.
//!
//! # What the distance is between
//!
//! Boundaries. A shape strictly inside another reports the gap between their
//! boundaries, not zero — whether a point is *inside* a solid is
//! [`classify_in_solid_exact`](crate::classify_in_solid_exact)'s question,
//! and conflating the two would make this answer wrong for the shells it is
//! right for.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, SurfaceGeometry, Transformable, TrimmedCurve};
use ogeom_math::{Point, Point2, Transform};
use ogeom_mesh::{Deflection, face_boundary, inside_boundary};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

use crate::classify::{distance_to_rings, parametric_band};
use crate::measure::{project_on_curve, project_on_surface};
use ogeom_intersect::ExtremaOptions;

/// One pair of nearest points, with the elements they lie on.
#[derive(Debug, Clone)]
pub struct ClosestPair {
    /// The nearest point on the first shape.
    pub point_a: Point,
    /// The nearest point on the second.
    pub point_b: Point,
    /// The vertex, edge or face of the first shape the point lies on.
    pub support_a: Shape,
    /// The same for the second shape.
    pub support_b: Shape,
}

/// The minimum distance between two shapes, with everywhere it is attained.
#[derive(Debug, Clone)]
pub struct ShapeDistance {
    /// The distance.
    pub distance: f64,
    /// Every pair of nearest points found within tolerance of the minimum.
    /// Parallel walls meet at a representative pair, not at every point of
    /// the overlap.
    pub pairs: Vec<ClosestPair>,
}

/// One element of a shape, with its geometry carried into world space.
enum Element {
    Vertex(Shape, Point),
    Edge(Shape, Box<Curve>),
    Face(Shape, Box<Prepared>),
}

/// A face's surface twice over: in world space for the extrema, and local
/// with its placement and rings for the trim test — parameters on a
/// transformed surface need not match the rings, which live in the stored
/// surface's parameter space, so the trim question is always asked locally.
struct Prepared {
    world: SurfaceGeometry,
    local: SurfaceGeometry,
    to_local: Transform,
    rings: Vec<Vec<Point2>>,
}

/// The minimum distance between two shapes' boundaries.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either shape
/// has no vertices, edges or faces to measure to;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn distance_between_shapes(
    model: &Model,
    a: &Shape,
    b: &Shape,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<ShapeDistance> {
    let ea = elements(model, a, tol)?;
    let eb = elements(model, b, tol)?;
    if ea.is_empty() || eb.is_empty() {
        ogeom_bail!(Construction, "a shape with no elements has no distance");
    }

    let mut candidates: Vec<(f64, ClosestPair)> = Vec::new();
    for element_a in &ea {
        for element_b in &eb {
            approach(element_a, element_b, options, tol, &mut candidates)?;
        }
    }
    let Some(least) = candidates
        .iter()
        .map(|(d, _)| *d)
        .min_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal))
    else {
        ogeom_bail!(
            NotDone,
            "no candidate approach was found between these shapes"
        );
    };

    let mut pairs: Vec<ClosestPair> = Vec::new();
    for (d, pair) in candidates {
        if d - least > tol.confusion() {
            continue;
        }
        // The same nearest pair arrives from several element pairs — a corner
        // is on a vertex, three edges and three faces at once. Keep the first
        // at each location.
        if pairs.iter().any(|known| {
            known.point_a.distance(pair.point_a) <= tol.confusion() * 1e2
                && known.point_b.distance(pair.point_b) <= tol.confusion() * 1e2
        }) {
            continue;
        }
        pairs.push(pair);
    }
    Ok(ShapeDistance {
        distance: least,
        pairs,
    })
}

/// Candidate approaches between one pair of elements.
fn approach(
    a: &Element,
    b: &Element,
    options: ExtremaOptions,
    tol: Tolerances,
    out: &mut Vec<(f64, ClosestPair)>,
) -> OgeomResult<()> {
    let mut push = |distance: f64, pa: Point, pb: Point, sa: &Shape, sb: &Shape| {
        out.push((
            distance,
            ClosestPair {
                point_a: pa,
                point_b: pb,
                support_a: sa.clone(),
                support_b: sb.clone(),
            },
        ));
    };
    match (a, b) {
        (Element::Vertex(sa, pa), Element::Vertex(sb, pb)) => {
            push(pa.distance(*pb), *pa, *pb, sa, sb);
        }
        (Element::Vertex(sa, pa), Element::Edge(sb, curve)) => {
            let foot = project_on_curve(curve, *pa, 64, tol)?;
            push(foot.distance, *pa, foot.point, sa, sb);
        }
        (Element::Edge(sa, curve), Element::Vertex(sb, pb)) => {
            let foot = project_on_curve(curve, *pb, 64, tol)?;
            push(foot.distance, foot.point, *pb, sa, sb);
        }
        (Element::Vertex(sa, pa), Element::Face(sb, face)) => {
            let foot = project_on_surface(&face.world, *pa, 32, tol)?;
            if inside_trim(face, foot.point, tol)? {
                push(foot.distance, *pa, foot.point, sa, sb);
            }
        }
        (Element::Face(sa, face), Element::Vertex(sb, pb)) => {
            let foot = project_on_surface(&face.world, *pb, 32, tol)?;
            if inside_trim(face, foot.point, tol)? {
                push(foot.distance, foot.point, *pb, sa, sb);
            }
        }
        (Element::Edge(sa, ca), Element::Edge(sb, cb)) => {
            let found = ogeom_intersect::extrema_curve_curve(ca, cb, options, tol)?;
            for near in &found.approaches {
                push(near.distance, near.point_a, near.point_b, sa, sb);
            }
        }
        (Element::Edge(sa, curve), Element::Face(sb, face)) => {
            let found = ogeom_intersect::extrema_curve_surface(curve, &face.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(face, near.point_b, tol)? {
                    push(near.distance, near.point_a, near.point_b, sa, sb);
                }
            }
        }
        (Element::Face(sa, face), Element::Edge(sb, curve)) => {
            let found = ogeom_intersect::extrema_curve_surface(curve, &face.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(face, near.point_b, tol)? {
                    push(near.distance, near.point_b, near.point_a, sa, sb);
                }
            }
        }
        (Element::Face(sa, fa), Element::Face(sb, fb)) => {
            let found =
                ogeom_intersect::extrema_surface_surface(&fa.world, &fb.world, options, tol)?;
            for near in &found.approaches {
                if inside_trim(fa, near.point_a, tol)? && inside_trim(fb, near.point_b, tol)? {
                    push(near.distance, near.point_a, near.point_b, sa, sb);
                }
            }
        }
    }
    Ok(())
}

/// Whether a world-space point on a face's surface lands inside its trimming.
///
/// Asked in the stored surface's own parameter space — the point is carried
/// into the face's frame and projected there, exactly as `classify_on_face`
/// does it, because rings and world-surface parameters need not agree under a
/// placement that scales. Too near the boundary counts as outside: the edge
/// pairs own that candidate and answer it exactly.
fn inside_trim(face: &Prepared, world_point: Point, tol: Tolerances) -> OgeomResult<bool> {
    let local = face.to_local.apply(world_point);
    let projection = project_on_surface(&face.local, local, 32, tol)?;
    let (u, v) = projection.parameters;
    let at = Point2::new(u, v);
    let band = parametric_band(&face.local, (u, v), tol.confusion() + RING_CHORD, tol);
    if distance_to_rings(&face.rings, at) <= band {
        return Ok(false);
    }
    Ok(inside_boundary(&face.rings, at))
}

/// The rings' polylining error, spatially — the trim test's uncertainty band.
const RING_CHORD: f64 = 1e-3;

/// Every element of a shape, with world-space geometry.
fn elements(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Element>> {
    let mut out = Vec::new();
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let Some(node) = model.node(&vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let placed = vertex.transform(model.datums())?.apply(data.point);
        out.push(Element::Vertex(vertex, placed));
    }
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(node) = model.node(&edge) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            // A degenerate edge has no extent of its own; its vertex and its
            // face carry its geometry.
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placement = edge.transform(model.datums())?;
        let trimmed: Curve = if (range.0, range.1) == {
            use ogeom_geom::Curve3d as _;
            geometry.domain()
        } {
            geometry.clone()
        } else {
            TrimmedCurve::new(geometry.clone(), range.0, range.1, tol)?.into()
        };
        out.push(Element::Edge(
            edge,
            Box::new(trimmed.transformed(&placement, tol)?),
        ));
    }
    let ring_deflection = Deflection {
        chord: RING_CHORD,
        angular: 0.05,
        ..Deflection::default()
    };
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        let rings = face_boundary(model, &face, ring_deflection, tol)?;
        // The stored surface may declare an enormous domain — a plane spans
        // ±1e9 — and the extrema layer rightly refuses to sample that. The
        // face only uses what its rings enclose, so the surface handed over
        // is trimmed to their parameter bound, with a margin for the rings'
        // own polylining, before being carried into world space.
        let restricted = restrict_to_rings(surface, &rings, tol)?;
        out.push(Element::Face(
            face,
            Box::new(Prepared {
                world: restricted.transformed(&placement, tol)?,
                local: surface.clone(),
                to_local: placement.inverse()?,
                rings,
            }),
        ));
    }
    Ok(out)
}

/// The surface restricted to the parameter rectangle its rings enclose.
///
/// The margin is proportional to the used span: the exact boundary lies
/// within the rings' polylining of it, and `inside_trim` already treats the
/// near-boundary band as the edges' territory, so the margin only has to
/// keep the whole face inside the restriction — it does not have to be
/// tight.
fn restrict_to_rings(
    surface: &SurfaceGeometry,
    rings: &[Vec<Point2>],
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    use ogeom_geom::Surface as _;
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut u = (f64::INFINITY, f64::NEG_INFINITY);
    let mut v = (f64::INFINITY, f64::NEG_INFINITY);
    for ring in rings {
        for p in ring {
            u = (u.0.min(p.x), u.1.max(p.x));
            v = (v.0.min(p.y), v.1.max(p.y));
        }
    }
    if u.0 > u.1 || v.0 > v.1 {
        // No rings — a naturally closed face uses its whole domain.
        return Ok(surface.clone());
    }
    let margin_u = (u.1 - u.0).mul_add(0.05, tol.parametric());
    let margin_v = (v.1 - v.0).mul_add(0.05, tol.parametric());
    let lo_u = (u.0 - margin_u).max(ua);
    let hi_u = (u.1 + margin_u).min(ub);
    let lo_v = (v.0 - margin_v).max(va);
    let hi_v = (v.1 + margin_v).min(vb);
    if lo_u >= hi_u || lo_v >= hi_v {
        return Ok(surface.clone());
    }
    Ok(ogeom_geom::TrimmedSurface::new(surface.clone(), (lo_u, hi_u), (lo_v, hi_v), tol)?.into())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{make_box, make_cylinder, make_sphere};
    use ogeom_math::{Direction, Frame, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn frame_at(origin: Point) -> Frame {
        Frame::new(origin, Direction::Z, Direction::X, T).unwrap()
    }

    #[test]
    fn parallel_box_walls_meet_at_the_gap_between_them() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(5.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!((found.distance - 3.0).abs() < 1e-9, "{}", found.distance);
        assert!(!found.pairs.is_empty());
        for pair in &found.pairs {
            assert!((pair.point_a.distance(pair.point_b) - found.distance).abs() < 1e-9);
        }
    }

    #[test]
    fn diagonal_boxes_meet_corner_to_corner() {
        // Offset along all three axes: the nearest points are two vertices,
        // exactly the candidates the geometry level declines to invent.
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(3.0, 3.0, 3.0)),
            (1.0, 1.0, 1.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        let exact = (3.0_f64 * 4.0).sqrt(); // corner (1,1,1) to corner (3,3,3)
        assert!((found.distance - exact).abs() < 1e-9);
        let pair = &found.pairs[0];
        assert!(pair.point_a.is_equal(Point::new(1.0, 1.0, 1.0), T));
        assert!(pair.point_b.is_equal(Point::new(3.0, 3.0, 3.0), T));
        assert_eq!(model.kind_of(&pair.support_a).unwrap(), ShapeType::Vertex);
        assert_eq!(model.kind_of(&pair.support_b).unwrap(), ShapeType::Vertex);
    }

    #[test]
    fn a_sphere_over_a_box_measures_to_the_top_face() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let ball = make_sphere(&mut model, frame_at(Point::new(2.0, 2.0, 4.0)), 1.0, T).unwrap();
        let found = distance_between_shapes(
            &model,
            &block.shape,
            &ball.shape,
            ExtremaOptions::default(),
            T,
        )
        .unwrap();
        assert!((found.distance - 2.0).abs() < 1e-7, "{}", found.distance);
        let pair = &found.pairs[0];
        assert!(pair.point_a.is_equal(Point::new(2.0, 2.0, 1.0), T));
        assert!(pair.point_b.is_equal(Point::new(2.0, 2.0, 3.0), T));
    }

    #[test]
    fn parallel_cylinders_meet_wall_to_wall() {
        // The nearest locus is a pair of facing rulings — a family at the
        // geometry level, a representative pair here, with the distance exact.
        let mut model = Model::new();
        let a = make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T).unwrap();
        let b =
            make_cylinder(&mut model, frame_at(Point::new(5.0, 0.0, 0.0)), 1.0, 4.0, T).unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!((found.distance - 3.0).abs() < 1e-7, "{}", found.distance);
    }

    #[test]
    fn touching_boxes_report_zero() {
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let b = make_box(
            &mut model,
            frame_at(Point::new(2.0, 0.0, 0.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        assert!(found.distance < 1e-9, "{}", found.distance);
    }

    #[test]
    fn a_box_inside_a_box_measures_boundary_to_boundary() {
        // Containment is the classifier's question. Distance is between
        // boundaries, and the gap between nested walls is what comes back.
        let mut model = Model::new();
        let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
        let inner = make_box(
            &mut model,
            frame_at(Point::new(2.0, 2.0, 2.0)),
            (2.0, 2.0, 2.0),
            T,
        )
        .unwrap();
        let found = distance_between_shapes(
            &model,
            &outer.shape,
            &inner.shape,
            ExtremaOptions::default(),
            T,
        )
        .unwrap();
        assert!((found.distance - 2.0).abs() < 1e-9, "{}", found.distance);
    }

    #[test]
    fn a_rotated_box_measures_edge_to_edge() {
        // Roll one box forty-five degrees about x and lift it: what faces the
        // top of the lower box is a single edge, and the nearest pair is that
        // edge against the top face.
        let mut model = Model::new();
        let a = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
        let tilted = Frame::new(
            Point::new(2.0, 2.0, 3.0),
            Direction::new(Vector::new(0.0, 1.0, 1.0), T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let b = make_box(&mut model, tilted, (1.0, 1.0, 1.0), T).unwrap();
        let found =
            distance_between_shapes(&model, &a.shape, &b.shape, ExtremaOptions::default(), T)
                .unwrap();
        // The tilted box's lowest feature is the edge its roll brings down to
        // z = 3 - 1/sqrt(2), facing the top face at z = 1.
        let exact = 2.0 - core::f64::consts::FRAC_1_SQRT_2;
        assert!((found.distance - exact).abs() < 1e-7, "{}", found.distance);
    }
}
