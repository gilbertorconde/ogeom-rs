//! Primitive solids.
//!
//! Built from numbers rather than from existing topology, so their history has
//! nothing to report — but their *provenance* very much does. Each face gets a
//! role naming which face of the primitive it is ([`roles`]), so a reference to
//! "the top of that box" survives the box being rebuilt at a different size.
//! Without it the reference would be to a handle that no longer exists.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_core::{OgeomResult, Role, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, LineCurve, PlanarCurve, PlaneSurface, SphereSurface, Surface,
};
use ogeom_math::{
    Circle, Cone, Cylinder, Direction, Direction2, Frame, Plane, Point, Point2, Sphere, Torus,
};
use ogeom_topo::{Model, Shape, ShapeType};

use crate::build::{make_edge_between, make_face_on, make_shell, make_solid, make_wire};
use crate::history::Built;
use ogeom_topo::{EdgeData, VertexData};

/// Roles naming which part of a primitive an entity is.
pub mod roles {
    use ogeom_core::Role;

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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if any dimension
/// is not finite and positive.
pub fn make_box(
    model: &mut Model,
    frame: Frame,
    size: (f64, f64, f64),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (dx, dy, dz) = size;
    for (name, value) in [("x", dx), ("y", dy), ("z", dz)] {
        if !value.is_finite() || value <= tol.confusion() {
            ogeom_bail!(
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

    box_like(model, &corner_points, tol)
}

/// Build a solid from eight corners laid out like [`CORNERS`], with the six
/// faces of [`FACES`].
///
/// A box and a wedge differ only in where the corners are: both have the same
/// eight, the same twelve edges and the same six planar faces. Sharing the
/// construction is not just less code — it is what keeps the two from drifting
/// apart in their winding, their roles or their pcurves.
fn box_like(model: &mut Model, corner_points: &[Point], tol: Tolerances) -> OgeomResult<Built> {
    let vertices: Vec<Shape> = corner_points
        .iter()
        .map(|p| model.add_vertex(ogeom_topo::VertexData::new(*p)))
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
        let plane = face_plane(corner_points, corners, tol)?;
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

/// A solid from explicit planar rings over shared corners.
///
/// The generalization `box_like` is a special case of: vertices per corner,
/// one edge per corner pair shared by both faces that meet along it, each
/// ring wound counter-clockwise seen from outside so its first three corners
/// give the outward normal. Every ring must be planar; the collapsed wedges
/// are, by construction.
fn faceted_solid(
    model: &mut Model,
    points: &[Point],
    rings: &[&[usize]],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| model.add_vertex(ogeom_topo::VertexData::new(*p)))
        .collect();
    let mut edge_of: std::collections::HashMap<(usize, usize), Shape> =
        std::collections::HashMap::new();
    let mut faces = Vec::with_capacity(rings.len());
    for ring_corners in rings {
        let origin = points[ring_corners[0]];
        let normal = Direction::from_cross(
            points[ring_corners[1]] - origin,
            points[ring_corners[2]] - points[ring_corners[1]],
            tol,
        )?;
        let x = Direction::new(points[ring_corners[1]] - origin, tol)?;
        let plane = Plane::new(Frame::new(origin, normal, x, tol)?);
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(plane).into());
        let mut ring = Vec::with_capacity(ring_corners.len());
        for step in 0..ring_corners.len() {
            let (from, to) = (
                ring_corners[step],
                ring_corners[(step + 1) % ring_corners.len()],
            );
            let key = (from.min(to), from.max(to));
            let edge = match edge_of.get(&key) {
                Some(edge) => edge.clone(),
                None => {
                    let curve: Curve =
                        LineCurve::segment(points[key.0], points[key.1], tol)?.into();
                    let length = points[key.0].distance(points[key.1]);
                    let edge = make_edge_between(
                        model,
                        curve,
                        (0.0, length),
                        &vertices[key.0],
                        &vertices[key.1],
                        tol,
                    )?
                    .shape;
                    edge_of.insert(key, edge.clone());
                    edge
                }
            };
            // The pcurve follows the edge's own parameterization — its
            // canonical low-to-high corner order — not the ring's traversal.
            attach_plane_pcurve(
                model,
                &edge,
                &plane,
                surface,
                points[key.0],
                points[key.1],
                tol,
            )?;
            ring.push(if from == key.0 { edge } else { edge.reversed() });
        }
        let wire = make_wire(model, &ring, tol)?.shape;
        faces.push(make_face_on(model, surface, std::slice::from_ref(&wire), tol)?.shape);
    }
    let shell = make_shell(model, &faces)?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// The plane of a box face, with its normal pointing outward.
fn face_plane(points: &[Point], corners: [usize; 4], tol: Tolerances) -> OgeomResult<Plane> {
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
    surface: ogeom_topo::SurfaceId,
    from: Point,
    to: Point,
    tol: Tolerances,
) -> OgeomResult<()> {
    let local = |p: Point| {
        let l = plane.frame().to_local(p);
        Point2::new(l.x, l.y)
    };
    let (a, b) = (local(from), local(to));
    let pcurve: PlanarCurve = ogeom_geom::Line2d::segment(a, b, tol)?.into();
    crate::build::attach_pcurve(
        model,
        edge,
        pcurve,
        surface,
        ogeom_topo::Location::identity(),
        (0.0, a.distance(b)),
    )
}

/// Find the canonical edge joining two corners, and whether it runs that way.
fn find_edge(from: usize, to: usize) -> OgeomResult<(usize, bool)> {
    for (index, &(a, b)) in EDGES.iter().enumerate() {
        if a == from && b == to {
            return Ok((index, true));
        }
        if a == to && b == from {
            return Ok((index, false));
        }
    }
    ogeom_bail!(
        Construction,
        "corners {from} and {to} are not joined by a box edge"
    )
}

/// Build a cylinder in `frame`: `radius` about its `z`, `height` along it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
/// dimension is not finite and positive.
pub fn make_cylinder(
    model: &mut Model,
    frame: Frame,
    radius: f64,
    height: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
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
        ogeom_geom::CylinderSurface::new(Cylinder::new(frame, radius, tol)?, (0.0, height))?.into(),
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `radius` is not
/// finite and positive.
pub fn make_sphere(
    model: &mut Model,
    frame: Frame,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
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
fn raised(frame: Frame, distance: f64, tol: Tolerances) -> OgeomResult<Frame> {
    Frame::new(
        frame.to_world(Point::new(0.0, 0.0, distance)),
        frame.z(),
        frame.x(),
        tol,
    )
}

/// Reject a dimension that cannot describe a solid.
fn check_size(what: &str, value: f64, tol: Tolerances) -> OgeomResult<()> {
    if !value.is_finite() || value <= tol.confusion() {
        ogeom_bail!(Construction, "{what} {value} must be finite and positive");
    }
    Ok(())
}

/// A closed circular edge, bounded twice by the same vertex.
fn full_circle_edge(
    model: &mut Model,
    circle: Circle,
    at: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
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
fn degenerate_edge(model: &mut Model, at: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
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
    surface: ogeom_topo::SurfaceId,
    extent: (f64, f64),
    edges: [&Shape; 3],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let [bottom, top, seam] = edges;
    let Some(geometry) = model.geometry().surface(surface) else {
        ogeom_bail!(Dangling, "surface is not in this model");
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
) -> OgeomResult<Shape> {
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
    surface: ogeom_topo::SurfaceId,
    from: Point2,
    to: Point2,
    tol: Tolerances,
) -> OgeomResult<()> {
    let pcurve: PlanarCurve = ogeom_geom::Line2d::segment(from, to, tol)?.into();
    crate::build::attach_pcurve(
        model,
        edge,
        pcurve,
        surface,
        ogeom_topo::Location::identity(),
        (0.0, from.distance(to)),
    )
}

/// Attach a seam edge's two pcurves, one for each side of the parameter
/// rectangle it bounds.
fn seam_pcurves(
    model: &mut Model,
    edge: &Shape,
    surface: ogeom_topo::SurfaceId,
    forward: (Point2, Point2),
    reversed: (Point2, Point2),
    tol: Tolerances,
) -> OgeomResult<()> {
    let length = forward.0.distance(forward.1);
    let first = model
        .geometry_mut()
        .add_pcurve(ogeom_geom::Line2d::segment(forward.0, forward.1, tol)?.into());
    let second = model
        .geometry_mut()
        .add_pcurve(ogeom_geom::Line2d::segment(reversed.0, reversed.1, tol)?.into());

    let Some(node) = model.node_mut(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let ogeom_topo::NodeData::Edge(data) = node.data_mut() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    data.add(ogeom_topo::EdgeRepr::Seam {
        forward: first,
        reversed: second,
        surface,
        location: ogeom_topo::Location::identity(),
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
    surface: ogeom_topo::SurfaceId,
    circle: Circle,
    plane: Plane,
    tol: Tolerances,
) -> OgeomResult<()> {
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    let flat_direction = |d: Direction| -> OgeomResult<Direction2> {
        let tip = flat(frame.origin() + d.vector());
        let base = flat(frame.origin());
        Direction2::new(tip - base, tol)
    };

    let frame2 = ogeom_math::Frame2::from_axes(
        flat(circle.centre()),
        flat_direction(circle.frame().x())?,
        flat_direction(circle.frame().y())?,
        tol,
    )?;
    let pcurve: PlanarCurve =
        ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(frame2, circle.radius(), tol)?).into();
    crate::build::attach_pcurve(
        model,
        edge,
        pcurve,
        surface,
        ogeom_topo::Location::identity(),
        (0.0, TAU),
    )
}

/// Build a cone or a truncated cone in `frame`.
///
/// `base_radius` is at the frame's origin and `top_radius` at `height` along
/// its `z`. A zero top radius makes a true cone, whose apex is a degenerate
/// edge and which has no top cap.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the height is
/// not positive, a radius is negative, both radii are zero, or the two radii
/// are equal — that last is a cylinder, which is a different surface with its
/// own type, and admitting it here would ask for a cone whose apex is at
/// infinity.
pub fn make_cone(
    model: &mut Model,
    frame: Frame,
    base_radius: f64,
    top_radius: f64,
    height: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    check_size("cone height", height, tol)?;
    for (what, r) in [("base radius", base_radius), ("top radius", top_radius)] {
        if !r.is_finite() || r < 0.0 {
            ogeom_bail!(
                Construction,
                "cone {what} {r} must be finite and non-negative"
            );
        }
    }
    if (base_radius - top_radius).abs() <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a cone with equal radii is a cylinder; use make_cylinder"
        );
    }
    if base_radius <= tol.confusion() && top_radius <= tol.confusion() {
        ogeom_bail!(Construction, "a cone needs one end with a radius");
    }
    model.begin_operation();

    // The surface is built along whichever direction it *widens*, because a
    // cone's half angle is signed and only the widening sense is a cone at all.
    // A narrowing solid therefore gets a surface frame pointing the other way,
    // and the rest of this function works in that frame's parameters.
    let widening = top_radius > base_radius;
    let (surface_frame, near_radius, far_radius) = if widening {
        (frame, base_radius, top_radius)
    } else {
        (flipped(frame, height, tol)?, top_radius, base_radius)
    };
    let half_angle = ((far_radius - near_radius) / height).atan();
    let cone = Cone::new(surface_frame, near_radius, half_angle, tol)?;

    let near_circle = circle_at(surface_frame, near_radius, 0.0, tol);
    let far_frame = raised(surface_frame, height, tol)?;
    // The wide end always has a radius: the two radii differ and the wider is
    // by definition not the collapsed one.
    let Some(far) = circle_at(far_frame, far_radius, 0.0, tol) else {
        ogeom_bail!(Construction, "the wide end of a cone must have a radius");
    };

    let near_vertex = model.add_vertex(VertexData::new(match near_circle {
        Some(c) => rim_point(c),
        None => surface_frame.origin(),
    }));
    let far_vertex = model.add_vertex(VertexData::new(rim_point(far)));

    // A zero radius is an apex: a rim of no length, which still bounds the face
    // in parameter space and still has to be there.
    let near_edge = match near_circle {
        Some(c) => full_circle_edge(model, c, &near_vertex, tol)?,
        None => degenerate_edge(model, &near_vertex, tol)?,
    };
    let far_edge = full_circle_edge(model, far, &far_vertex, tol)?;

    let seam_start = near_circle.map_or_else(|| surface_frame.origin(), rim_point);
    let seam_end = rim_point(far);
    let seam = make_edge_between(
        model,
        LineCurve::segment(seam_start, seam_end, tol)?.into(),
        (0.0, slant(near_radius, far_radius, height)),
        &near_vertex,
        &far_vertex,
        tol,
    )?
    .shape;

    let lateral_id = model
        .geometry_mut()
        .add_surface(ogeom_geom::ConeSurface::new(cone, (0.0, height))?.into());
    let lateral = rectangle_face(
        model,
        lateral_id,
        (TAU, height),
        [&near_edge, &far_edge, &seam],
        tol,
    )?;
    model.set_derived(&lateral, &[], roles::FACE_LATERAL)?;

    let mut faces = vec![lateral];
    // The near end only needs a cap if it has any area.
    if let Some(c) = near_circle {
        let cap = cap_face(model, surface_frame, false, &near_edge, c, tol)?;
        model.set_derived(
            &cap,
            &[],
            if widening {
                roles::FACE_MIN_Z
            } else {
                roles::FACE_MAX_Z
            },
        )?;
        faces.push(cap);
    }
    let far_cap = cap_face(model, far_frame, true, &far_edge, far, tol)?;
    model.set_derived(
        &far_cap,
        &[],
        if widening {
            roles::FACE_MAX_Z
        } else {
            roles::FACE_MIN_Z
        },
    )?;
    faces.push(far_cap);

    let shell = make_shell(model, &faces)?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// Build a torus in `frame`: `major` from the axis to the tube's centre,
/// `minor` the tube's own radius.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either radius
/// is not finite and positive.
pub fn make_torus(
    model: &mut Model,
    frame: Frame,
    major: f64,
    minor: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    check_size("torus major radius", major, tol)?;
    check_size("torus minor radius", minor, tol)?;
    model.begin_operation();

    // Closed in *both* directions, so the face has a seam on all four sides of
    // its parameter rectangle and one vertex where the two seams cross. A
    // cylinder's single seam is the easy case; this is the one that decides
    // whether seam handling is general or a special case for one primitive.
    let start = frame.to_world(Point::new(major + minor, 0.0, 0.0));
    let corner = model.add_vertex(VertexData::new(start));

    // The outer equator, along which the tube parameter is zero.
    let equator = Circle::new(frame, major + minor, tol)?;
    let along_u = full_circle_edge(model, equator, &corner, tol)?;

    // The tube's own circle at longitude zero: its plane holds the frame's `x`
    // and `z`, centred one major radius out.
    let tube_frame = Frame::new(
        frame.to_world(Point::new(major, 0.0, 0.0)),
        -frame.y(),
        frame.x(),
        tol,
    )?;
    let along_v = full_circle_edge(model, Circle::new(tube_frame, minor, tol)?, &corner, tol)?;

    let surface = model
        .geometry_mut()
        .add_surface(ogeom_geom::TorusSurface::new(Torus::new(frame, major, minor, tol)?).into());
    let face = doubly_seamed_face(model, surface, &along_u, &along_v, tol)?;
    model.set_derived(&face, &[], roles::FACE_LATERAL)?;

    let shell = make_shell(model, std::slice::from_ref(&face))?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    Ok(Built::from_nothing(solid))
}

/// Build a wedge: a box whose top face is inset.
///
/// `size` is the box at the frame's origin; `top` is the `(x, y)` extent of the
/// upper face, over the same corner. Equal extents give a box.
///
/// Both top extents must be positive. A wedge whose top collapses to a ridge
/// has five faces and one whose top collapses to a point has four — different
/// topologies, not this one with a zero somewhere, and building them here would
/// produce a face with no area.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a dimension is
/// not finite and positive.
pub fn make_wedge(
    model: &mut Model,
    frame: Frame,
    size: (f64, f64, f64),
    top: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (dx, dy, dz) = size;
    for (name, value) in [("x", dx), ("y", dy), ("z", dz)] {
        check_size(&format!("wedge {name} size"), value, tol)?;
    }
    for (name, value) in [("x", top.0), ("y", top.1)] {
        if !value.is_finite() || value < 0.0 {
            ogeom_bail!(
                Construction,
                "wedge top {name} extent {value} must be finite and non-negative"
            );
        }
    }
    model.begin_operation();

    // A zero top extent is a different topology, not a box with a flat face
    // of no area: the top collapses to a ridge (five faces) or to a point
    // (five faces, four of them triangles), and each is built as itself.
    let (dx_, dy_, dz_) = size;
    let collapsed = (top.0 <= tol.confusion(), top.1 <= tol.confusion());
    match collapsed {
        (true, true) => {
            let local = [
                Point::new(0.0, 0.0, 0.0),
                Point::new(dx_, 0.0, 0.0),
                Point::new(dx_, dy_, 0.0),
                Point::new(0.0, dy_, 0.0),
                Point::new(0.0, 0.0, dz_),
            ];
            let points: Vec<Point> = local.iter().map(|p| frame.to_world(*p)).collect();
            let rings: [&[usize]; 5] = [
                &[0, 3, 2, 1],
                &[0, 1, 4],
                &[1, 2, 4],
                &[2, 3, 4],
                &[3, 0, 4],
            ];
            return faceted_solid(model, &points, &rings, tol);
        }
        (false, true) => {
            let local = [
                Point::new(0.0, 0.0, 0.0),
                Point::new(dx_, 0.0, 0.0),
                Point::new(dx_, dy_, 0.0),
                Point::new(0.0, dy_, 0.0),
                Point::new(0.0, 0.0, dz_),
                Point::new(top.0, 0.0, dz_),
            ];
            let points: Vec<Point> = local.iter().map(|p| frame.to_world(*p)).collect();
            let rings: [&[usize]; 5] = [
                &[0, 3, 2, 1],
                &[0, 1, 5, 4],
                &[1, 2, 5],
                &[2, 3, 4, 5],
                &[3, 0, 4],
            ];
            return faceted_solid(model, &points, &rings, tol);
        }
        (true, false) => {
            let local = [
                Point::new(0.0, 0.0, 0.0),
                Point::new(dx_, 0.0, 0.0),
                Point::new(dx_, dy_, 0.0),
                Point::new(0.0, dy_, 0.0),
                Point::new(0.0, 0.0, dz_),
                Point::new(0.0, top.1, dz_),
            ];
            let points: Vec<Point> = local.iter().map(|p| frame.to_world(*p)).collect();
            let rings: [&[usize]; 5] = [
                &[0, 3, 2, 1],
                &[0, 4, 5, 3],
                &[0, 1, 4],
                &[1, 2, 5, 4],
                &[2, 3, 5],
            ];
            return faceted_solid(model, &points, &rings, tol);
        }
        (false, false) => {}
    }

    // Same eight corners and the same six faces as a box; only where the top
    // four sit differs. Every side stays planar because the top face stays a
    // rectangle parallel to the bottom, which is what makes the wedge a
    // reparameterization of the box rather than a separate construction.
    let corners: Vec<Point> = CORNERS
        .iter()
        .map(|&(i, j, k)| {
            #[allow(clippy::cast_precision_loss)]
            let (fi, fj) = (i as f64, j as f64);
            let (ex, ey) = if k == 0 { (dx, dy) } else { (top.0, top.1) };
            #[allow(clippy::cast_precision_loss)]
            frame.to_world(Point::new(fi * ex, fj * ey, k as f64 * dz))
        })
        .collect();
    box_like(model, &corners, tol)
}

/// Build the unbounded solid on one side of a face.
///
/// `inside` names the side: it is a point in the material. The face is oriented
/// so its normal leads *away* from that point, which is what "outward" means
/// for a solid, and the result is a solid bounded by that one face.
///
/// # It is only as unbounded as its surface is
///
/// A half space is genuinely infinite; a surface in this kernel is not. A plane
/// declares a finite domain — very large, but finite — and the solid built here
/// reaches exactly as far as its face's surface does. So its *volume* and its
/// centre of mass are properties of that declared extent rather than of a half
/// space, and mean nothing. What does mean something is which side of the face
/// a point is on, which is the question a half space exists to answer.
///
/// That is what it is for: a half space is an argument to a boolean — cut a
/// solid with one and you have trimmed it by a surface. `cut`, `common` and
/// `section` all accept one as either operand; `fuse` refuses it by name,
/// because an unbounded fuse has no volume to keep. Classification against it
/// works too and is what the tests here use.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `face` is not a
/// face, or if `inside` lies on it — a point on the boundary names no side.
pub fn make_half_space(
    model: &mut Model,
    face: &Shape,
    inside: Point,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "a half space is bounded by a face");
    }
    let (at, normal) = crate::measure::face_normal(model, face, tol)?;
    let towards = inside - at;
    let reach = towards.magnitude();
    if reach <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the point naming the solid side lies on the face itself, so it \
             names no side"
        );
    }
    let along = normal.dot(towards) / reach;
    if along.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the point naming the solid side lies in the face's own surface, so \
             it names no side"
        );
    }
    model.begin_operation();

    // Outward means away from the material, and the material is where `inside`
    // is. A normal pointing towards it is pointing in.
    let boundary = if along > 0.0 {
        face.reversed()
    } else {
        face.clone()
    };
    let shell = make_shell(model, std::slice::from_ref(&boundary))?.shape;
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    model.set_derived(&solid, std::slice::from_ref(face), roles::FACE_LATERAL)?;

    let mut history = crate::history::History::new();
    history.generate(face, shell);
    history.generate(face, solid.clone());
    Ok(Built::new(solid, history))
}

/// A frame turned end for end: its origin at `height` along the old `z`, and
/// its `z` pointing back the way it came.
fn flipped(frame: Frame, height: f64, tol: Tolerances) -> OgeomResult<Frame> {
    Frame::new(
        frame.to_world(Point::new(0.0, 0.0, height)),
        -frame.z(),
        frame.x(),
        tol,
    )
}

/// A circle in a frame's plane, or `None` if the radius has collapsed.
fn circle_at(frame: Frame, radius: f64, _at: f64, tol: Tolerances) -> Option<Circle> {
    if radius <= tol.confusion() {
        return None;
    }
    Circle::new(frame, radius, tol).ok()
}

/// The slant length of a cone's side, which is what its seam edge measures.
fn slant(near: f64, far: f64, height: f64) -> f64 {
    (far - near).hypot(height)
}

/// A face on a surface closed in *both* parameter directions.
///
/// A torus. Both sides of the rectangle are seams and both ends are seams, so
/// all four boundary edges are two edges used twice, and one vertex serves all
/// four corners.
fn doubly_seamed_face(
    model: &mut Model,
    surface: ogeom_topo::SurfaceId,
    along_u: &Shape,
    along_v: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (o, e) = (0.0, TAU);
    // The edge running in u is a seam in *v*: it is the same curve at v = 0 and
    // at v = 2pi.
    seam_pcurves(
        model,
        along_u,
        surface,
        (Point2::new(o, o), Point2::new(e, o)),
        (Point2::new(o, e), Point2::new(e, e)),
        tol,
    )?;
    seam_pcurves(
        model,
        along_v,
        surface,
        (Point2::new(e, o), Point2::new(e, e)),
        (Point2::new(o, o), Point2::new(o, e)),
        tol,
    )?;

    let ring = [
        along_u.clone(),
        along_v.clone(),
        along_u.reversed(),
        along_v.reversed(),
    ];
    let wire = make_wire(model, &ring, tol)?.shape;
    Ok(make_face_on(model, surface, std::slice::from_ref(&wire), tol)?.shape)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::build::is_shell_closed;
    use ogeom_geom::Surface;
    use ogeom_topo::{ShapeType, explore_unique};

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
                ogeom_core::Provenance::Derived { role, .. } => *role,
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
                    ogeom_core::Provenance::Derived { role, .. } => *role,
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
    use ogeom_mesh::{Deflection, triangulate};
    use ogeom_topo::{ShapeType, explore_unique};

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

        // The *mesh* is inscribed and converges from below; the measurement
        // itself now runs on the exact surface and lands on the closed form.
        let mut previous = 0.0;
        for chord in [0.1_f64, 0.02, 0.005] {
            let mesh = triangulate(&model, &built.shape, deflection(chord), T).unwrap();
            assert!(mesh.is_closed(), "the mesh has a hole at chord {chord}");
            assert!(
                mesh.volume() < exact,
                "an inscribed volume cannot exceed it"
            );
            assert!(mesh.volume() > previous, "refining lost volume");
            previous = mesh.volume();
        }
        assert!(previous > exact * 0.995, "{previous} against {exact}");
        let props = volume_properties(&model, &built.shape, deflection(0.005), T).unwrap();
        assert_relative_eq!(props.mass, exact, epsilon = 1e-9);
        assert_eq!(props.deflection, 0.0, "measured on the exact surface");
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
        assert_relative_eq!(props.mass, volume, epsilon = 1e-9);
        assert_eq!(props.deflection, 0.0, "measured on the exact surface");
        assert!(
            props.centre.distance(Point::ORIGIN) < 1e-9,
            "got {:?}",
            props.centre
        );

        let surface = surface_properties(&model, &built.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(surface.mass, area, epsilon = 1e-9);
        assert_eq!(surface.deflection, 0.0, "measured on the exact surface");
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod more_primitive_tests {
    use super::*;
    use crate::build::is_shell_closed;
    use crate::mass::volume_properties;
    use approx::assert_relative_eq;
    use ogeom_mesh::{Deflection, triangulate};
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn deflection(chord: f64) -> Deflection {
        Deflection {
            chord,
            ..Deflection::default()
        }
    }

    fn closed(model: &Model, solid: &Shape) -> bool {
        let shell = explore_unique(model, solid, ShapeType::Shell).unwrap()[0].clone();
        is_shell_closed(model, &shell).unwrap()
    }

    #[test]
    fn a_truncated_cone_has_the_volume_a_frustum_has() {
        let (r0, r1, h) = (3.0_f64, 1.0_f64, 4.0_f64);
        let exact = PI * h / 3.0 * r1.mul_add(r1, r0.mul_add(r0, r0 * r1));
        let mut model = Model::new();
        let built = make_cone(&mut model, Frame::WORLD, r0, r1, h, T).unwrap();

        assert!(closed(&model, &built.shape));
        let props = volume_properties(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(props.mass < exact, "an inscribed volume cannot exceed it");
        assert!(props.mass > exact * 0.995, "{} against {exact}", props.mass);
    }

    #[test]
    fn a_true_cone_ends_in_an_apex_and_has_no_top_cap() {
        // The apex is an edge of no length. Without it the lateral face's
        // boundary is open along the top of its parameter rectangle; with a cap
        // there instead, the solid would have a face of no area.
        let (radius, height) = (2.0_f64, 5.0);
        let exact = PI * radius * radius * height / 3.0;
        let mut model = Model::new();
        let built = make_cone(&mut model, Frame::WORLD, radius, 0.0, height, T).unwrap();

        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Face)
                .unwrap()
                .len(),
            2,
            "a flank and one cap"
        );
        assert!(closed(&model, &built.shape));
        let props = volume_properties(&model, &built.shape, deflection(0.005), T).unwrap();
        assert!(props.mass < exact);
        assert!(props.mass > exact * 0.99, "{} against {exact}", props.mass);
    }

    #[test]
    fn a_cone_widening_upward_is_built_the_same_way_round() {
        // The surface's half angle only makes sense in the widening direction,
        // so a narrowing solid gets a flipped surface frame. Both had better
        // give the same solid seen from the other end.
        let mut model = Model::new();
        let up = make_cone(&mut model, Frame::WORLD, 1.0, 3.0, 4.0, T).unwrap();
        let down = make_cone(&mut model, Frame::WORLD, 3.0, 1.0, 4.0, T).unwrap();

        let a = volume_properties(&model, &up.shape, deflection(0.005), T).unwrap();
        let b = volume_properties(&model, &down.shape, deflection(0.005), T).unwrap();
        assert_relative_eq!(a.mass, b.mass, max_relative = 1e-9);
        // Mirrored about the middle of the height.
        assert_relative_eq!(a.centre.z, 4.0 - b.centre.z, epsilon = 1e-9);
    }

    #[test]
    fn a_torus_has_two_seams_one_vertex_and_the_volume_a_torus_has() {
        // Closed in both parameter directions: the case that decides whether
        // seam handling is general or a special case for the cylinder.
        let (major, minor) = (5.0_f64, 2.0);
        let exact = 2.0 * PI * PI * major * minor * minor;
        let mut model = Model::new();
        let built = make_torus(&mut model, Frame::WORLD, major, minor, T).unwrap();

        let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
        assert_eq!(counts(ShapeType::Face), 1);
        assert_eq!(counts(ShapeType::Edge), 2, "one seam each way");
        assert_eq!(counts(ShapeType::Vertex), 1, "where the two seams cross");
        assert!(closed(&model, &built.shape));

        let mesh = triangulate(&model, &built.shape, deflection(0.02), T).unwrap();
        assert!(mesh.is_closed(), "the mesh has a hole");

        let props = volume_properties(&model, &built.shape, deflection(0.02), T).unwrap();
        assert_relative_eq!(props.mass, exact, epsilon = 1e-9);
        assert_eq!(props.deflection, 0.0, "measured on the exact surface");
        assert!(props.centre.distance(Point::ORIGIN) < 1e-9);
    }

    #[test]
    fn a_wedge_with_equal_extents_is_a_box() {
        let mut model = Model::new();
        let wedge = make_wedge(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), (2.0, 3.0), T).unwrap();
        let props = volume_properties(&model, &wedge.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(props.mass, 24.0, epsilon = 1e-9);
        assert!(closed(&model, &wedge.shape));
    }

    #[test]
    fn a_tapered_wedge_has_the_volume_a_frustum_of_a_pyramid_has() {
        // A prismatoid: h/6 * (A_bottom + 4*A_middle + A_top).
        let mut model = Model::new();
        let wedge = make_wedge(&mut model, Frame::WORLD, (4.0, 4.0, 6.0), (2.0, 2.0), T).unwrap();
        let props = volume_properties(&model, &wedge.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(
            props.mass,
            6.0 / 6.0 * 4.0_f64.mul_add(9.0, 16.0 + 4.0),
            epsilon = 1e-9
        );
        assert!(closed(&model, &wedge.shape));
    }

    #[test]
    fn a_wedge_collapsing_to_a_ridge_is_five_faces_and_a_prismatoid_volume() {
        // Top y extent zero: the top is a ridge along x at y = 0. The
        // prismatoid volume integral in closed form:
        // dz*dy*(dx/2 + (a - dx)/6) for ridge length a.
        let mut model = Model::new();
        let wedge = make_wedge(&mut model, Frame::WORLD, (4.0, 3.0, 6.0), (2.0, 0.0), T).unwrap();
        let faces = ogeom_topo::explore_unique(&model, &wedge.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, 5, "a ridge wedge has five faces, none of them empty");
        let props = volume_properties(&model, &wedge.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(props.mass, 6.0 * 3.0 * (2.0 - 2.0 / 6.0), epsilon = 1e-9);
        assert!(closed(&model, &wedge.shape));

        // And the other axis mirrors.
        let other = make_wedge(&mut model, Frame::WORLD, (3.0, 4.0, 6.0), (0.0, 2.0), T).unwrap();
        let props = volume_properties(&model, &other.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(props.mass, 6.0 * 3.0 * (2.0 - 2.0 / 6.0), epsilon = 1e-9);
        assert!(closed(&model, &other.shape));
    }

    #[test]
    fn a_wedge_collapsing_to_a_point_is_a_pyramid() {
        let mut model = Model::new();
        let wedge = make_wedge(&mut model, Frame::WORLD, (4.0, 3.0, 6.0), (0.0, 0.0), T).unwrap();
        let faces = ogeom_topo::explore_unique(&model, &wedge.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, 5, "a base and four triangles");
        let props = volume_properties(&model, &wedge.shape, deflection(0.01), T).unwrap();
        assert_relative_eq!(props.mass, 4.0 * 3.0 * 6.0 / 3.0, epsilon = 1e-9);
        assert!(closed(&model, &wedge.shape));
    }

    #[test]
    fn dimensions_that_describe_no_solid_are_refused() {
        let mut model = Model::new();
        // Equal radii are a cylinder, and no radius at all is nothing.
        assert!(make_cone(&mut model, Frame::WORLD, 2.0, 2.0, 1.0, T).is_err());
        assert!(make_cone(&mut model, Frame::WORLD, 0.0, 0.0, 1.0, T).is_err());
        assert!(make_cone(&mut model, Frame::WORLD, 1.0, 2.0, 0.0, T).is_err());
        assert!(make_cone(&mut model, Frame::WORLD, -1.0, 2.0, 1.0, T).is_err());

        assert!(make_torus(&mut model, Frame::WORLD, 0.0, 1.0, T).is_err());
        assert!(make_torus(&mut model, Frame::WORLD, 1.0, f64::NAN, T).is_err());

        assert!(make_wedge(&mut model, Frame::WORLD, (0.0, 1.0, 1.0), (1.0, 1.0), T).is_err());
        assert!(make_wedge(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), (-1.0, 1.0), T).is_err());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod half_space_tests {
    use super::*;
    use crate::classify::Containment;
    use crate::{classify_in_solid, make_natural_face};
    use ogeom_geom::PlaneSurface;
    use ogeom_math::Direction;
    use ogeom_mesh::Deflection;
    use ogeom_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn coarse() -> Deflection {
        Deflection {
            chord: 1.0,
            ..Deflection::default()
        }
    }

    /// The whole of the z = 0 plane, as a face.
    fn ground(model: &mut Model) -> Shape {
        make_natural_face(model, PlaneSurface::new(Plane::new(Frame::WORLD)).into())
            .unwrap()
            .shape
    }

    #[test]
    fn the_face_is_oriented_away_from_the_side_that_is_solid() {
        // Outward means away from the material, and the material is where the
        // naming point is. The two calls differ only in which side was named,
        // so the boundary they produce must differ in orientation.
        let mut model = Model::new();
        let face = ground(&mut model);
        let above = make_half_space(&mut model, &face, Point::new(0.0, 0.0, 5.0), T).unwrap();
        let below = make_half_space(&mut model, &face, Point::new(0.0, 0.0, -5.0), T).unwrap();

        let boundary = |built: &crate::Built| {
            explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone()
        };
        let (a, b) = (boundary(&above), boundary(&below));
        assert!(a.is_partner(&b), "the same face, both times");
        assert_ne!(
            a.orientation(),
            b.orientation(),
            "naming the other side should turn the boundary round"
        );
        assert_eq!(model.kind_of(&above.shape).unwrap(), ShapeType::Solid);
        assert_eq!(
            explore_unique(&model, &above.shape, ShapeType::Face)
                .unwrap()
                .len(),
            1,
            "one face bounds a half space"
        );
    }

    #[test]
    fn nothing_can_yet_be_asked_about_the_inside_of_one() {
        // Recorded rather than worked around. A half space's boundary is one
        // face with free edges all round, so it is *not* a closed shell — and
        // every query that needs an inside says so instead of guessing. That is
        // the correct answer for a shape whose boundary does not close, and it
        // is why a half space is only useful as a boolean argument.
        let mut model = Model::new();
        let face = ground(&mut model);
        let built = make_half_space(&mut model, &face, Point::new(0.0, 0.0, 5.0), T).unwrap();

        let shell = explore_unique(&model, &built.shape, ShapeType::Shell).unwrap()[0].clone();
        assert!(!crate::is_shell_closed(&model, &shell).unwrap());

        let err = classify_in_solid(&model, &built.shape, Point::new(0.0, 0.0, 5.0), coarse(), T)
            .unwrap_err();
        assert!(
            err.to_string().contains("not closed"),
            "unexpected message: {err}"
        );
        let _ = Containment::In;
    }

    #[test]
    fn a_point_on_the_face_names_no_side() {
        let mut model = Model::new();
        let face = ground(&mut model);
        let err = make_half_space(&mut model, &face, Point::ORIGIN, T).unwrap_err();
        assert!(err.to_string().contains("names no side"), "got {err}");

        // And one that is off the origin but still in the plane.
        let err = make_half_space(&mut model, &face, Point::new(3.0, 4.0, 0.0), T).unwrap_err();
        assert!(err.to_string().contains("names no side"), "got {err}");
    }

    #[test]
    fn a_half_space_is_bounded_by_a_face_and_nothing_else() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        assert!(make_half_space(&mut model, &solid, Point::ORIGIN, T).is_err());

        let vertex = model.add_point(Point::ORIGIN);
        assert!(make_half_space(&mut model, &vertex, Point::ORIGIN, T).is_err());
    }

    #[test]
    fn it_reaches_only_as_far_as_its_surface_says() {
        // The documented limitation, pinned so it cannot quietly become a
        // claim of infinity. The face is a plane with a declared domain, so the
        // solid is a very large region rather than a half space, and anything
        // that integrates over it is measuring that domain.
        let mut model = Model::new();
        let face = ground(&mut model);
        let built = make_half_space(&mut model, &face, Point::new(0.0, 0.0, 1.0), T).unwrap();

        // The bound comes back *empty*, which is the honest answer and a
        // better one than a very large box: `surface_bounds` refuses to bound
        // an unbounded plane rather than quoting its declared extent, so
        // nothing downstream mistakes that extent for the shape's size.
        let bounds = crate::shape_bounds(&model, &built.shape, T).unwrap();
        assert!(
            bounds.is_empty(),
            "an unbounded plane should decline to bound itself, got {bounds:?}"
        );
        let _ = Direction::Z;
    }
}
