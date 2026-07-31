//! Primitive solids.
//!
//! Built from numbers rather than from existing topology, so their history has
//! nothing to report — but their *provenance* very much does. Each face gets a
//! role naming which face of the primitive it is ([`roles`]), so a reference to
//! "the top of that box" survives the box being rebuilt at a different size.
//! Without it the reference would be to a handle that no longer exists.

use og_core::{OgResult, Role, Tolerances, og_bail};
use og_geom::{Curve, LineCurve, PlanarCurve, PlaneSurface, SurfaceGeometry};
use og_math::{Direction, Frame, Plane, Point, Point2};
use og_topo::{Model, Shape};

use crate::build::{make_edge_between, make_face, make_shell, make_solid, make_wire};
use crate::history::Built;

/// Roles naming which part of a primitive an entity is.
pub mod roles {
    use og_core::Role;

    /// The face at the low end of the frame's `x` axis.
    pub const FACE_MIN_X: Role = Role::op_defined(10);
    /// The face at the high end of the frame's `x` axis.
    pub const FACE_MAX_X: Role = Role::op_defined(11);
    /// The face at the low end of the frame's `y` axis.
    pub const FACE_MIN_Y: Role = Role::op_defined(12);
    /// The face at the high end of the frame's `y` axis.
    pub const FACE_MAX_Y: Role = Role::op_defined(13);
    /// The face at the low end of the frame's `z` axis.
    pub const FACE_MIN_Z: Role = Role::op_defined(14);
    /// The face at the high end of the frame's `z` axis.
    pub const FACE_MAX_Z: Role = Role::op_defined(15);
}

/// The eight corners of a box, indexed so that bit 0 is `x`, bit 1 is `y` and
/// bit 2 is `z`.
const CORNERS: [(usize, usize, usize); 8] = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
];

/// The twelve edges, as corner index pairs in a canonical direction.
const EDGES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0), // bottom ring
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4), // top ring
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7), // verticals
];

/// The six faces, each as four corner indices wound counter-clockwise *seen
/// from outside* — which is what makes every face normal point outward and the
/// shell consistently oriented.
///
/// Getting a winding backwards produces a solid that is inside out along one
/// face. Nothing about the geometry says so; the first thing to notice is a
/// volume that comes out negative, or a boolean that keeps the wrong side.
const FACES: [([usize; 4], Role); 6] = [
    ([0, 3, 2, 1], roles::FACE_MIN_Z),
    ([4, 5, 6, 7], roles::FACE_MAX_Z),
    ([0, 1, 5, 4], roles::FACE_MIN_Y),
    ([2, 3, 7, 6], roles::FACE_MAX_Y),
    ([0, 4, 7, 3], roles::FACE_MIN_X),
    ([1, 2, 6, 5], roles::FACE_MAX_X),
];

/// Build an axis-aligned box in `frame`, spanning `size` along each of its
/// axes from the frame's origin.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if any dimension
/// is not finite and positive.
pub fn make_box(
    model: &mut Model,
    frame: Frame,
    size: (f64, f64, f64),
    tol: Tolerances,
) -> OgResult<Built> {
    let (dx, dy, dz) = size;
    for (name, value) in [("x", dx), ("y", dy), ("z", dz)] {
        if !value.is_finite() || value <= tol.confusion() {
            og_bail!(
                Construction,
                "box {name} size {value} must be finite and positive"
            );
        }
    }
    model.begin_operation();

    let extent = [dx, dy, dz];
    let corner_points: Vec<Point> = CORNERS
        .iter()
        .map(|&(i, j, k)| {
            #[allow(clippy::cast_precision_loss)]
            let local = [
                i as f64 * extent[0],
                j as f64 * extent[1],
                k as f64 * extent[2],
            ];
            frame.to_world(Point::new(local[0], local[1], local[2]))
        })
        .collect();

    let vertices: Vec<Shape> = corner_points
        .iter()
        .map(|p| model.add_vertex(og_topo::VertexData::new(*p)))
        .collect();

    // One edge per pair, shared by the two faces that meet along it. Building
    // an edge per face instead would leave every edge used once, and the shell
    // would not close.
    let mut edges = Vec::with_capacity(EDGES.len());
    for &(from, to) in &EDGES {
        let curve: Curve = LineCurve::segment(corner_points[from], corner_points[to], tol)?.into();
        let length = corner_points[from].distance(corner_points[to]);
        edges.push(
            make_edge_between(
                model,
                curve,
                (0.0, length),
                &vertices[from],
                &vertices[to],
                tol,
            )?
            .shape,
        );
    }

    let mut faces = Vec::with_capacity(FACES.len());
    for (corners, role) in FACES {
        let plane = face_plane(&corner_points, corners, tol)?;
        let mut ring = Vec::with_capacity(4);
        for step in 0..4 {
            let (from, to) = (corners[step], corners[(step + 1) % 4]);
            let (index, forward) = find_edge(from, to)?;
            ring.push(if forward {
                edges[index].clone()
            } else {
                edges[index].reversed()
            });

            // The pcurve follows the edge's own parameterization, not the
            // face's traversal of it, so the two representations agree on what
            // a parameter means. Following the traversal instead would leave a
            // reversed edge's pcurve running backwards against its curve.
            let (canonical_from, canonical_to) = EDGES[index];
            attach_plane_pcurve(
                model,
                &edges[index],
                &plane,
                corner_points[canonical_from],
                corner_points[canonical_to],
                tol,
            )?;
        }

        let wire = make_wire(model, &ring, tol)?.shape;
        let surface: SurfaceGeometry = PlaneSurface::new(plane).into();
        let face = make_face(model, surface, std::slice::from_ref(&wire), tol)?.shape;
        model.set_derived(&face, &[], role)?;
        faces.push(face);
    }

    let shell = make_shell(model, &faces)?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// The plane of a box face, with its normal pointing outward.
fn face_plane(points: &[Point], corners: [usize; 4], tol: Tolerances) -> OgResult<Plane> {
    let origin = points[corners[0]];
    // The winding is counter-clockwise seen from outside, so the right-hand
    // rule over the first three corners gives the outward normal.
    let normal = Direction::from_cross(
        points[corners[1]] - origin,
        points[corners[2]] - points[corners[1]],
        tol,
    )?;
    let x = Direction::new(points[corners[1]] - origin, tol)?;
    Ok(Plane::new(Frame::new(origin, normal, x, tol)?))
}

/// Attach a line pcurve running between two points, expressed in a plane's
/// parameter space.
fn attach_plane_pcurve(
    model: &mut Model,
    edge: &Shape,
    plane: &Plane,
    from: Point,
    to: Point,
    tol: Tolerances,
) -> OgResult<()> {
    let local = |p: Point| {
        let l = plane.frame().to_local(p);
        Point2::new(l.x, l.y)
    };
    let (a, b) = (local(from), local(to));
    let pcurve: PlanarCurve = og_geom::Line2d::segment(a, b, tol)?.into();
    let surface = model
        .geometry_mut()
        .add_surface(PlaneSurface::new(*plane).into());
    crate::build::attach_pcurve(model, edge, pcurve, surface, (0.0, a.distance(b)))
}

/// Find the canonical edge joining two corners, and whether it runs that way.
fn find_edge(from: usize, to: usize) -> OgResult<(usize, bool)> {
    for (index, &(a, b)) in EDGES.iter().enumerate() {
        if a == from && b == to {
            return Ok((index, true));
        }
        if a == to && b == from {
            return Ok((index, false));
        }
    }
    og_bail!(
        Construction,
        "corners {from} and {to} are not joined by a box edge"
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::build::is_shell_closed;
    use og_geom::Surface;
    use og_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_box_has_the_topology_a_box_should_have() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let solid = &built.shape;

        assert_eq!(model.kind_of(solid).unwrap(), ShapeType::Solid);
        assert_eq!(
            explore_unique(&model, solid, ShapeType::Shell)
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            explore_unique(&model, solid, ShapeType::Face)
                .unwrap()
                .len(),
            6
        );
        assert_eq!(
            explore_unique(&model, solid, ShapeType::Wire)
                .unwrap()
                .len(),
            6
        );
        assert_eq!(
            explore_unique(&model, solid, ShapeType::Edge)
                .unwrap()
                .len(),
            12,
            "edges are shared between adjacent faces, not duplicated per face"
        );
        assert_eq!(
            explore_unique(&model, solid, ShapeType::Vertex)
                .unwrap()
                .len(),
            8
        );
    }

    #[test]
    fn a_boxs_shell_is_closed() {
        // Every edge used by exactly two faces. Building an edge per face would
        // give twenty-four edges each used once, and nothing would enclose.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap();
        assert!(is_shell_closed(&model, &shell[0]).unwrap());
    }

    #[test]
    fn every_face_normal_points_out_of_the_box() {
        // A backwards winding gives a solid that is inside out along one face.
        // Nothing about the geometry says so, and the first symptom is usually
        // a negative volume or a boolean keeping the wrong side.
        let mut model = Model::new();
        let size = (2.0, 3.0, 4.0);
        let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
        let centre = Point::new(size.0 / 2.0, size.1 / 2.0, size.2 / 2.0);

        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 6);
        for face in &faces {
            let node = model.node(face).unwrap();
            let data = node.data().as_face().unwrap();
            let surface = model.geometry().surface(data.surface).unwrap();
            let ((ua, ub), (va, vb)) = surface.domain();
            let point = surface
                .point_at((ua + ub) / 2.0, (va + vb) / 2.0, T)
                .unwrap();
            let normal = surface
                .normal_at((ua + ub) / 2.0, (va + vb) / 2.0, T)
                .unwrap();

            // A plane's parameter origin is a box corner, so sample at the
            // corner itself and check the normal leads away from the centre.
            let outward = surface.point_at(0.0, 0.0, T).unwrap() - centre;
            assert!(
                normal.dot_vector(outward) > 0.0,
                "a face normal points inward: at {point:?}, normal {normal:?}"
            );
        }
    }

    #[test]
    fn the_six_faces_carry_distinct_roles() {
        // The roles are what a rebuild matches against: "the top of that box"
        // has to survive the box being rebuilt at a different size, and a
        // handle will not.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T).unwrap();
        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();

        let mut roles: Vec<Role> = faces
            .iter()
            .map(|f| match model.provenance_of(f).unwrap() {
                og_core::Provenance::Derived { role, .. } => *role,
                other => panic!("expected a derived face, got {other:?}"),
            })
            .collect();
        roles.sort_unstable();
        roles.dedup();
        assert_eq!(roles.len(), 6, "each face is identifiable on its own");
    }

    #[test]
    fn rebuilding_at_a_different_size_gives_the_faces_the_same_roles() {
        // The point of provenance: a reference to a face survives a parameter
        // change, because the role is what identifies it and the role does not
        // depend on the size.
        let roles_of = |size| {
            let mut model = Model::new();
            let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
            let mut roles: Vec<Role> = explore_unique(&model, &built.shape, ShapeType::Face)
                .unwrap()
                .iter()
                .map(|f| match model.provenance_of(f).unwrap() {
                    og_core::Provenance::Derived { role, .. } => *role,
                    other => panic!("expected a derived face, got {other:?}"),
                })
                .collect();
            roles.sort_unstable();
            roles
        };
        assert_eq!(roles_of((1.0, 1.0, 1.0)), roles_of((10.0, 0.5, 7.0)));
    }

    #[test]
    fn every_edge_carries_a_pcurve_for_each_face_it_bounds() {
        // Face splitting during a boolean happens in parameter space, so an
        // edge without a pcurve on a face is an edge that face cannot be split
        // along.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T).unwrap();
        for edge in explore_unique(&model, &built.shape, ShapeType::Edge).unwrap() {
            let data = model.node(&edge).unwrap().data().as_edge().unwrap();
            assert_eq!(
                data.parametric_surfaces().len(),
                2,
                "a box edge borders exactly two faces"
            );
            assert!(data.curve3d().is_some());
        }
    }

    #[test]
    fn a_box_in_a_tilted_frame_is_still_a_box() {
        let frame = Frame::new(
            Point::new(5.0, -2.0, 1.0),
            Direction::from_coords(1.0, 1.0, 1.0, T).unwrap(),
            Direction::from_coords(1.0, -1.0, 0.0, T).unwrap(),
            T,
        )
        .unwrap();
        let mut model = Model::new();
        let built = make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();

        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Face)
                .unwrap()
                .len(),
            6
        );
        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap();
        assert!(is_shell_closed(&model, &shell[0]).unwrap());

        // The corners are where the frame puts them.
        let vertices = explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap();
        let origin_corner = frame.to_world(Point::ORIGIN);
        assert!(
            vertices.iter().any(|v| {
                model
                    .node(v)
                    .unwrap()
                    .data()
                    .as_vertex()
                    .unwrap()
                    .point
                    .is_equal(origin_corner, T)
            }),
            "no vertex at the frame origin"
        );
    }

    #[test]
    fn degenerate_dimensions_are_refused() {
        let mut model = Model::new();
        for size in [
            (0.0, 1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, f64::NAN),
            (f64::INFINITY, 1.0, 1.0),
        ] {
            assert!(
                make_box(&mut model, Frame::WORLD, size, T).is_err(),
                "accepted {size:?}"
            );
        }
        assert!(make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).is_ok());
    }

    #[test]
    fn a_primitive_reports_no_history_because_it_consumed_nothing() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        assert!(
            built.history.is_empty(),
            "built from numbers, so there are no inputs to report on"
        );
    }

    #[test]
    fn every_corner_pair_of_a_face_is_a_real_edge() {
        // Guards the tables themselves: a face naming two corners that share no
        // edge would silently build a wire from the wrong ones.
        for (corners, _) in FACES {
            for step in 0..4 {
                let (from, to) = (corners[step], corners[(step + 1) % 4]);
                assert!(
                    find_edge(from, to).is_ok(),
                    "face corners {from} and {to} are not joined"
                );
            }
        }
        assert!(find_edge(0, 6).is_err(), "opposite corners share no edge");
    }

    #[test]
    fn every_edge_is_used_by_exactly_two_faces_in_the_tables() {
        // Checked against the tables directly, independently of the model: if
        // this were wrong the shell would not close, and the failure would
        // point at the topology rather than at the data that caused it.
        let mut uses = [0_usize; EDGES.len()];
        for (corners, _) in FACES {
            for step in 0..4 {
                let (from, to) = (corners[step], corners[(step + 1) % 4]);
                let (index, _) = find_edge(from, to).unwrap();
                uses[index] += 1;
            }
        }
        assert!(uses.iter().all(|&n| n == 2), "edge use counts: {uses:?}");
    }
}
