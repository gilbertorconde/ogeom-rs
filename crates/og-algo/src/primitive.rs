//! Primitive solids.
//!
//! Built from numbers rather than from existing topology, so their history has
//! nothing to report — but their *provenance* very much does. Each face gets a
//! role naming which face of the primitive it is ([`roles`]), so a reference to
//! "the top of that box" survives the box being rebuilt at a different size.
//! Without it the reference would be to a handle that no longer exists.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use og_core::{OgResult, Role, Tolerances, og_bail};
use og_geom::{CircleCurve, Curve, LineCurve, PlanarCurve, PlaneSurface, SphereSurface, Surface};
use og_math::{Circle, Cylinder, Direction, Direction2, Frame, Plane, Point, Point2, Sphere};
use og_topo::{Model, Shape};

use crate::build::{make_edge_between, make_face_on, make_shell, make_solid, make_wire};
use crate::history::Built;
use og_topo::{EdgeData, VertexData};

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
    /// The face swept around the frame's `z` axis — a cylinder's side, a
    /// sphere's whole surface, a cone's flank.
    pub const FACE_LATERAL: Role = Role::op_defined(16);
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
        // One id for this face's surface, shared by the face and by every
        // pcurve on it. Registering it per pcurve would give each its own id,
        // and the face would find no pcurve on itself.
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(plane).into());
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
                surface,
                corner_points[canonical_from],
                corner_points[canonical_to],
                tol,
            )?;
        }

        let wire = make_wire(model, &ring, tol)?.shape;
        let face = make_face_on(model, surface, std::slice::from_ref(&wire), tol)?.shape;
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
    surface: og_topo::SurfaceId,
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

/// Build a cylinder in `frame`: `radius` about its `z`, `height` along it.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if either
/// dimension is not finite and positive.
pub fn make_cylinder(
    model: &mut Model,
    frame: Frame,
    radius: f64,
    height: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    check_size("cylinder radius", radius, tol)?;
    check_size("cylinder height", height, tol)?;
    model.begin_operation();

    let top_frame = raised(frame, height, tol)?;
    let bottom_circle = Circle::new(frame, radius, tol)?;
    let top_circle = Circle::new(top_frame, radius, tol)?;

    // One vertex on each rim, where the seam meets it. A full circle has to be
    // bounded somewhere or it cannot join a wire, and putting the bound on the
    // seam is what lets the lateral face close.
    let low = model.add_vertex(VertexData::new(rim_point(bottom_circle)));
    let high = model.add_vertex(VertexData::new(rim_point(top_circle)));

    let bottom_edge = full_circle_edge(model, bottom_circle, &low, tol)?;
    let top_edge = full_circle_edge(model, top_circle, &high, tol)?;
    let seam = make_edge_between(
        model,
        LineCurve::segment(rim_point(bottom_circle), rim_point(top_circle), tol)?.into(),
        (0.0, height),
        &low,
        &high,
        tol,
    )?
    .shape;

    let lateral_id = model.geometry_mut().add_surface(
        og_geom::CylinderSurface::new(Cylinder::new(frame, radius, tol)?, (0.0, height))?.into(),
    );
    let lateral = rectangle_face(
        model,
        lateral_id,
        (TAU, height),
        [&bottom_edge, &top_edge, &seam],
        tol,
    )?;
    model.set_derived(&lateral, &[], roles::FACE_LATERAL)?;

    let bottom = cap_face(model, frame, false, &bottom_edge, bottom_circle, tol)?;
    model.set_derived(&bottom, &[], roles::FACE_MIN_Z)?;
    let top = cap_face(model, top_frame, true, &top_edge, top_circle, tol)?;
    model.set_derived(&top, &[], roles::FACE_MAX_Z)?;

    let shell = make_shell(model, &[lateral, bottom, top])?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// Build a sphere of `radius` centred at `frame`'s origin.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `radius` is not
/// finite and positive.
pub fn make_sphere(
    model: &mut Model,
    frame: Frame,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    check_size("sphere radius", radius, tol)?;
    model.begin_operation();

    let centre = frame.origin();
    let south = model.add_vertex(VertexData::new(centre - frame.z().vector() * radius));
    let north = model.add_vertex(VertexData::new(centre + frame.z().vector() * radius));

    // The seam is the meridian at longitude zero, pole to pole through the
    // frame's `x`. Its plane is spanned by `x` and `z`, so the circle's normal
    // is `-y` — which makes its angle parameter the latitude exactly, and the
    // mapping onto the surface's `v` the identity rather than a rescaling.
    let meridian = Circle::new(Frame::new(centre, -frame.y(), frame.x(), tol)?, radius, tol)?;
    let seam = make_edge_between(
        model,
        CircleCurve::new(meridian).into(),
        (-FRAC_PI_2, FRAC_PI_2),
        &south,
        &north,
        tol,
    )?
    .shape;

    // The poles bound the face in parameter space and have no length in space.
    // Dropping them would leave the boundary open along the top and bottom of
    // the parameter rectangle, with nothing for the triangulator to trim to.
    let bottom_edge = degenerate_edge(model, &south, tol)?;
    let top_edge = degenerate_edge(model, &north, tol)?;

    let surface = model
        .geometry_mut()
        .add_surface(SphereSurface::new(Sphere::new(frame, radius, tol)?).into());
    let face = rectangle_face(
        model,
        surface,
        (TAU, PI),
        [&bottom_edge, &top_edge, &seam],
        tol,
    )?;
    // The sphere's `v` runs from -pi/2, not from zero, so the rectangle's
    // corner is not at the parameter origin.
    model.set_derived(&face, &[], roles::FACE_LATERAL)?;

    let shell = make_shell(model, std::slice::from_ref(&face))?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// Where a circle's own parameterization starts: its frame's `x`, one radius
/// out. The seam has to meet the rim exactly there, not merely nearby.
fn rim_point(circle: Circle) -> Point {
    circle.centre() + circle.frame().x().vector() * circle.radius()
}

/// A frame moved along its own `z`.
fn raised(frame: Frame, distance: f64, tol: Tolerances) -> OgResult<Frame> {
    Frame::new(
        frame.to_world(Point::new(0.0, 0.0, distance)),
        frame.z(),
        frame.x(),
        tol,
    )
}

/// Reject a dimension that cannot describe a solid.
fn check_size(what: &str, value: f64, tol: Tolerances) -> OgResult<()> {
    if !value.is_finite() || value <= tol.confusion() {
        og_bail!(Construction, "{what} {value} must be finite and positive");
    }
    Ok(())
}

/// A closed circular edge, bounded twice by the same vertex.
fn full_circle_edge(
    model: &mut Model,
    circle: Circle,
    at: &Shape,
    tol: Tolerances,
) -> OgResult<Shape> {
    Ok(make_edge_between(
        model,
        CircleCurve::new(circle).into(),
        (0.0, TAU),
        at,
        at,
        tol,
    )?
    .shape)
}

/// An edge with no length, bounded twice by the same vertex.
///
/// A pole or an apex: it bounds a face in parameter space and collapses to a
/// point in space. It carries no 3D curve, because there is no curve to carry —
/// its pcurve is the whole story, and the `degenerate` flag says so rather than
/// leaving a caller to notice the missing representation.
fn degenerate_edge(model: &mut Model, at: &Shape, tol: Tolerances) -> OgResult<Shape> {
    let _ = tol;
    let mut data = EdgeData::new();
    data.degenerate = true;
    model.add_edge(data, &[at.clone(), at.clone()])
}

/// A face covering a rectangle of a surface's parameter space, bounded by two
/// edges across and one seam up both sides.
///
/// The shape every surface of revolution has: `u` closes on itself, so the face
/// is a rectangle whose left and right sides are the *same* edge seen twice.
/// That edge carries two pcurves and appears in the wire twice, once each way.
///
/// `extent` is the size of the rectangle; its lower corner comes from the
/// surface's own domain, so a sphere's `v` starting at `-pi/2` needs no special
/// case here.
fn rectangle_face(
    model: &mut Model,
    surface: og_topo::SurfaceId,
    extent: (f64, f64),
    edges: [&Shape; 3],
    tol: Tolerances,
) -> OgResult<Shape> {
    let [bottom, top, seam] = edges;
    let Some(geometry) = model.geometry().surface(surface) else {
        og_bail!(Dangling, "surface is not in this model");
    };
    let ((ua, _), (va, _)) = geometry.domain();
    let (du, dv) = extent;
    let (ub, vb) = (ua + du, va + dv);

    line_pcurve(
        model,
        bottom,
        surface,
        Point2::new(ua, va),
        Point2::new(ub, va),
        tol,
    )?;
    line_pcurve(
        model,
        top,
        surface,
        Point2::new(ua, vb),
        Point2::new(ub, vb),
        tol,
    )?;
    seam_pcurves(
        model,
        seam,
        surface,
        (Point2::new(ub, va), Point2::new(ub, vb)),
        (Point2::new(ua, va), Point2::new(ua, vb)),
        tol,
    )?;

    // Counter-clockwise around the rectangle: across the bottom, up the far
    // side of the seam, back across the top, down the near side.
    let ring = [
        bottom.clone(),
        seam.clone(),
        top.reversed(),
        seam.reversed(),
    ];
    let wire = make_wire(model, &ring, tol)?.shape;
    Ok(make_face_on(model, surface, std::slice::from_ref(&wire), tol)?.shape)
}

/// A planar cap closing one end of a solid of revolution.
///
/// `outward` says whether the cap's normal follows the frame's `z` or opposes
/// it — which is the difference between a solid and one that is inside out
/// along one face, and nothing in the geometry says which was meant.
fn cap_face(
    model: &mut Model,
    frame: Frame,
    outward: bool,
    rim: &Shape,
    circle: Circle,
    tol: Tolerances,
) -> OgResult<Shape> {
    let normal = if outward { frame.z() } else { -frame.z() };
    let plane = Plane::new(Frame::new(frame.origin(), normal, frame.x(), tol)?);
    let surface = model
        .geometry_mut()
        .add_surface(PlaneSurface::new(plane).into());

    circle_pcurve_on_plane(model, rim, surface, circle, plane, tol)?;

    // The rim runs one way round; the cap that faces the other way walks it
    // backwards, so its boundary is traversed consistently with its normal.
    let edge = if outward { rim.clone() } else { rim.reversed() };
    let wire = make_wire(model, std::slice::from_ref(&edge), tol)?.shape;
    Ok(make_face_on(model, surface, std::slice::from_ref(&wire), tol)?.shape)
}

/// Attach a straight pcurve running between two parameter points.
fn line_pcurve(
    model: &mut Model,
    edge: &Shape,
    surface: og_topo::SurfaceId,
    from: Point2,
    to: Point2,
    tol: Tolerances,
) -> OgResult<()> {
    let pcurve: PlanarCurve = og_geom::Line2d::segment(from, to, tol)?.into();
    crate::build::attach_pcurve(model, edge, pcurve, surface, (0.0, from.distance(to)))
}

/// Attach a seam edge's two pcurves, one for each side of the parameter
/// rectangle it bounds.
fn seam_pcurves(
    model: &mut Model,
    edge: &Shape,
    surface: og_topo::SurfaceId,
    forward: (Point2, Point2),
    reversed: (Point2, Point2),
    tol: Tolerances,
) -> OgResult<()> {
    let length = forward.0.distance(forward.1);
    let first = model
        .geometry_mut()
        .add_pcurve(og_geom::Line2d::segment(forward.0, forward.1, tol)?.into());
    let second = model
        .geometry_mut()
        .add_pcurve(og_geom::Line2d::segment(reversed.0, reversed.1, tol)?.into());

    let Some(node) = model.node_mut(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let og_topo::NodeData::Edge(data) = node.data_mut() else {
        og_bail!(Construction, "edge node holds no edge data");
    };
    data.add(og_topo::EdgeRepr::Seam {
        forward: first,
        reversed: second,
        surface,
        location: og_topo::Location::identity(),
        range: (0.0, length),
    });
    Ok(())
}

/// Attach the pcurve of a circle lying in a plane.
///
/// Built from the circle's own axes expressed in the plane's frame, rather than
/// from the angle alone. A cap whose normal opposes the circle's sees the same
/// circle running the other way, and taking the axes through the conversion is
/// what makes that fall out instead of needing a sign to be remembered.
fn circle_pcurve_on_plane(
    model: &mut Model,
    edge: &Shape,
    surface: og_topo::SurfaceId,
    circle: Circle,
    plane: Plane,
    tol: Tolerances,
) -> OgResult<()> {
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    let flat_direction = |d: Direction| -> OgResult<Direction2> {
        let tip = flat(frame.origin() + d.vector());
        let base = flat(frame.origin());
        Direction2::new(tip - base, tol)
    };

    let frame2 = og_math::Frame2::from_axes(
        flat(circle.centre()),
        flat_direction(circle.frame().x())?,
        flat_direction(circle.frame().y())?,
        tol,
    )?;
    let pcurve: PlanarCurve =
        og_geom::Circle2d::new(og_math::Circle2::new(frame2, circle.radius(), tol)?).into();
    crate::build::attach_pcurve(model, edge, pcurve, surface, (0.0, TAU))
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod revolution_tests {
    use super::*;
    use crate::build::is_shell_closed;
    use crate::mass::{surface_properties, volume_properties};
    use approx::assert_relative_eq;
    use og_mesh::{Deflection, triangulate};
    use og_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn deflection(chord: f64) -> Deflection {
        Deflection {
            chord,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_cylinder_has_the_topology_a_cylinder_has() {
        let mut model = Model::new();
        let built = make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 3, "a side and two caps");
        assert_eq!(counts(ShapeType::Edge), 3, "two rims and one seam");
        assert_eq!(counts(ShapeType::Vertex), 2, "one on each rim");

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(
            is_shell_closed(&model, &shell).unwrap(),
            "every edge should be used an even number of times"
        );
    }

    #[test]
    fn a_cylinder_tessellates_into_a_closed_mesh_of_the_right_size() {
        // The seam is what this proves. Without both of its pcurves the lateral
        // face's boundary does not close in parameter space, and the mesh comes
        // out with a slot down one side.
        let (radius, height) = (2.0_f64, 5.0);
        let exact = PI * radius * radius * height;
        let mut model = Model::new();
        let built = make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();

        let mut previous = 0.0;
        for chord in [0.1_f64, 0.02, 0.005] {
            let mesh = triangulate(&model, &built.shape, deflection(chord), T).unwrap();
            assert!(mesh.is_closed(), "the mesh has a hole at chord {chord}");
            let props = volume_properties(&model, &built.shape, deflection(chord), T).unwrap();
            assert!(props.mass < exact, "an inscribed volume cannot exceed it");
            assert!(props.mass > previous, "refining lost volume");
            previous = props.mass;
        }
        assert!(previous > exact * 0.995, "{previous} against {exact}");
    }

    #[test]
    fn a_cylinders_caps_face_outward() {
        // A cap wound the wrong way makes the solid inside out along one face,
        // and the volume comes out short by exactly that cap's contribution
        // rather than obviously wrong.
        let mut model = Model::new();
        let built = make_cylinder(&mut model, Frame::WORLD, 1.0, 3.0, T).unwrap();
        let props = volume_properties(&model, &built.shape, deflection(0.005), T).unwrap();

        assert!(
            props.centre.distance(Point::new(0.0, 0.0, 1.5)) < 1e-3,
            "the centre of a cylinder is halfway up its axis, got {:?}",
            props.centre
        );
    }

    #[test]
    fn a_sphere_has_the_topology_a_sphere_has() {
        let mut model = Model::new();
        let built = make_sphere(&mut model, Frame::WORLD, 3.0, T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 1, "one surface covers a sphere");
        assert_eq!(counts(ShapeType::Edge), 3, "a seam and two poles");
        assert_eq!(counts(ShapeType::Vertex), 2, "the two poles");

        // The poles have no length and say so, rather than leaving a caller to
        // discover it by dividing by their length.
        let degenerate = explore_unique(&model, &built.shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
            })
            .count();
        assert_eq!(degenerate, 2);
    }

    #[test]
    fn a_sphere_converges_on_the_volume_and_area_a_sphere_has() {
        let radius = 4.0_f64;
        let volume = 4.0 / 3.0 * PI * radius.powi(3);
        let area = 4.0 * PI * radius * radius;
        let mut model = Model::new();
        let built = make_sphere(&mut model, Frame::WORLD, radius, T).unwrap();

        let props = volume_properties(&model, &built.shape, deflection(0.01), T).unwrap();
        assert!(props.mass < volume);
        assert!(
            props.mass > volume * 0.995,
            "{} against {volume}",
            props.mass
        );
        assert!(
            props.centre.distance(Point::ORIGIN) < 1e-3,
            "got {:?}",
            props.centre
        );

        let surface = surface_properties(&model, &built.shape, deflection(0.01), T).unwrap();
        assert!(surface.mass < area);
        assert!(surface.mass > area * 0.995);
    }

    #[test]
    fn a_placed_primitive_lands_where_it_was_placed() {
        let frame = Frame::new(Point::new(10.0, -5.0, 2.0), Direction::X, Direction::Y, T).unwrap();
        let mut model = Model::new();
        let built = make_cylinder(&mut model, frame, 1.0, 4.0, T).unwrap();
        let props = volume_properties(&model, &built.shape, deflection(0.005), T).unwrap();

        // Half way along the frame's own z, which here is world +x.
        assert!(
            props.centre.distance(Point::new(12.0, -5.0, 2.0)) < 1e-3,
            "got {:?}",
            props.centre
        );
        // Inscribed, so a little under the exact pi r^2 h.
        assert_relative_eq!(props.mass, PI * 4.0, max_relative = 0.01);
    }

    #[test]
    fn dimensions_that_describe_no_solid_are_refused() {
        let mut model = Model::new();
        for (r, h) in [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (f64::NAN, 1.0)] {
            assert!(make_cylinder(&mut model, Frame::WORLD, r, h, T).is_err());
        }
        for r in [0.0, -1.0, f64::INFINITY] {
            assert!(make_sphere(&mut model, Frame::WORLD, r, T).is_err());
        }
    }
}
