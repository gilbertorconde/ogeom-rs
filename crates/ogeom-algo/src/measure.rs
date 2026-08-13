//! Bounding volumes and projection.
//!
//! # Bounds must contain what they claim to
//!
//! A bounding box in a kernel is always a *rejection* test. Too large costs
//! time; too small silently drops a real intersection, and nothing downstream
//! can tell it happened. So every bound here is derived from a property that
//! guarantees containment, never from sampling:
//!
//! - a line's bound is its endpoints, which is exact;
//! - a spline's is its control points, which is guaranteed by the convex hull
//!   property — the curve never leaves the hull of its control polygon;
//! - an analytic curve or surface's is computed from its own definition;
//! - a *trimmed* piece falls back to the bound of the whole, which is loose but
//!   never wrong.
//!
//! Sampling a curve at a few parameters and taking the extremes is the obvious
//! alternative and is not sound: the curve bulges between the samples, and the
//! amount it bulges is exactly what a bound is supposed to capture.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve2d, Curve3d, PlanarCurve, Surface, SurfaceGeometry, curve::LINE_EXTENT,
};
use ogeom_math::{Aabb, Direction, Frame, Point, Point2, Vector, solve};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

/// A guaranteed bound for a space curve.
///
/// Loose for a trimmed curve, which reports the bound of the whole rather than
/// of the piece — never wrong, and tightening it would mean solving for the
/// extremes of the trimmed range, which is the same work as an intersection.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the curve cannot be
/// evaluated at its own domain ends.
pub fn curve_bounds(curve: &Curve, tol: Tolerances) -> OgeomResult<Aabb> {
    Ok(match curve {
        // Exact: a segment is the hull of its ends.
        Curve::Line(_) => Aabb::of_corners(curve.start(tol)?, curve.end(tol)?),

        // A conic's extremes are its centre displaced by the radii along each
        // frame axis, which bounds the whole conic whatever arc is in use.
        Curve::Circle(c) => {
            let circle = c.circle();
            frame_bounds(
                circle.centre(),
                circle.frame(),
                (circle.radius(), circle.radius(), 0.0),
            )
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            frame_bounds(
                ellipse.centre(),
                ellipse.frame(),
                (ellipse.major_radius(), ellipse.minor_radius(), 0.0),
            )
        }

        // A helix never leaves its cylinder, and its rise over the trimmed
        // angle interval is linear — so the cylinder's box over that rise
        // contains it, tight along the axis and whole-circle-loose across,
        // the same convention the circle uses.
        Curve::Helix(h) => {
            let rise = h.pitch() / core::f64::consts::TAU;
            let slope = h.taper() / core::f64::consts::TAU;
            let (a, b) = h.domain();
            let reach = slope
                .mul_add(a, h.radius())
                .abs()
                .max(slope.mul_add(b, h.radius()).abs());
            let mid = h.frame().origin() + h.frame().z().vector() * (rise * f64::midpoint(a, b));
            frame_bounds(
                mid,
                *h.frame(),
                (reach, reach, (rise * (b - a) / 2.0).abs()),
            )
        }

        // An offset never strays farther than its distance from the basis:
        // the basis's bound grown by |d| on every axis contains it, exactly
        // the guarantee and no tighter.
        Curve::Offset(o) => curve_bounds(o.basis(), tol)?.expanded(o.distance().abs()),

        // A surface curve never leaves its surface, whose bound is already
        // guaranteed.
        Curve::OnSurface(c) => surface_bounds(c.surface(), tol)?,

        // A hyperbola and a parabola are unbounded, so only the trimmed extent
        // has a bound at all. Both are convex in their own frame, so the hull
        // of the two ends and the vertex contains the arc between them.
        Curve::Hyperbola(_) | Curve::Parabola(_) => {
            let (a, b) = curve.domain();
            let mid = curve.point_at(f64::midpoint(a, b), tol)?;
            let ends = Aabb::of_corners(curve.start(tol)?, curve.end(tol)?);
            // The midpoint is the extreme in the frame's x direction for both,
            // and the ends bound the rest.
            ends.with_point(mid)
        }

        // The convex hull property: a B-spline never leaves the hull of its
        // control polygon, so the polygon's box contains the curve exactly.
        Curve::BSpline(s) => Aabb::of_points(
            &s.control_points()
                .iter()
                .map(|w| w.point())
                .collect::<Vec<_>>(),
        ),

        Curve::Trimmed(t) => curve_bounds(t.basis(), tol)?,
    })
}

/// A guaranteed bound for a surface.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the surface cannot be
/// evaluated over its own domain.
pub fn surface_bounds(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Aabb> {
    Ok(match surface {
        // An offset never strays farther than its distance from the basis.
        SurfaceGeometry::Offset(o) => surface_bounds(o.basis(), tol)?.expanded(o.distance().abs()),

        // A plane is unbounded; its declared domain is what there is to bound,
        // and the four corners span it exactly.
        SurfaceGeometry::Plane(p) => {
            let ((ua, ub), (va, vb)) = p.domain();
            if ua <= -LINE_EXTENT || ub >= LINE_EXTENT {
                ogeom_bail!(
                    Domain,
                    "an unbounded plane has no finite bound; trim it before asking"
                );
            }
            let mut out = Aabb::EMPTY;
            for (u, v) in [(ua, va), (ua, vb), (ub, va), (ub, vb)] {
                out = out.with_point(p.point_at(u, v, tol)?);
            }
            out
        }

        SurfaceGeometry::Cylinder(c) => {
            let cyl = c.cylinder();
            let ((_, _), (va, vb)) = c.domain();
            let frame = cyl.frame();
            let base = frame.origin() + frame.z() * va;
            let top = frame.origin() + frame.z() * vb;
            let radial = frame_bounds(base, frame, (cyl.radius(), cyl.radius(), 0.0));
            radial.union(&frame_bounds(top, frame, (cyl.radius(), cyl.radius(), 0.0)))
        }

        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let ((_, _), (va, vb)) = c.domain();
            let frame = cone.frame();
            let mut out = Aabb::EMPTY;
            for height in [va, vb] {
                let radius = cone.radius_at(height).abs();
                let centre = frame.origin() + frame.z() * height;
                out = out.union(&frame_bounds(centre, frame, (radius, radius, 0.0)));
            }
            out
        }

        // A sphere's bound is its centre plus its radius on every axis, whatever
        // patch of it is in use.
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let r = Vector::splat(sphere.radius());
            Aabb::of_corners(sphere.centre() - r, sphere.centre() + r)
        }

        SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let reach = torus.major_radius() + torus.minor_radius();
            frame_bounds(
                torus.centre(),
                torus.frame(),
                (reach, reach, torus.minor_radius()),
            )
        }

        // Convex hull property again, in two directions.
        SurfaceGeometry::BSpline(s) => Aabb::of_points(
            &s.grid()
                .points()
                .iter()
                .map(|w| w.point())
                .collect::<Vec<_>>(),
        ),

        // A revolved curve reaches at most its own furthest distance from the
        // axis, in every direction around it.
        SurfaceGeometry::Revolution(r) => {
            let curve = curve_bounds(r.curve(), tol)?;
            let axis = r.axis();
            let mut reach: f64 = 0.0;
            let mut along = Aabb::EMPTY;
            for corner in curve.corners() {
                reach = reach.max(axis.distance_to(corner));
                along = along.with_point(axis.project(corner));
            }
            let radial = Vector::splat(reach);
            along.expanded(0.0).union(&Aabb::of_corners(
                along.low().unwrap_or(axis.location) - radial,
                along.high().unwrap_or(axis.location) + radial,
            ))
        }

        // A swept curve reaches the curve's bound at each end of the sweep.
        SurfaceGeometry::Extrusion(e) => {
            let base = curve_bounds(e.curve(), tol)?;
            let ((_, _), (va, vb)) = e.domain();
            let start = base.transformed(&ogeom_math::Transform::translation(e.direction() * va));
            let end = base.transformed(&ogeom_math::Transform::translation(e.direction() * vb));
            start.union(&end)
        }

        SurfaceGeometry::Trimmed(t) => surface_bounds(t.basis(), tol)?,
    })
}

/// The box of a point displaced by `(x, y, z)` extents along a frame's axes.
///
/// Every axis of the result gets the sum of the absolute contributions from all
/// three frame directions, which is what makes it a bound rather than an
/// estimate: a tilted frame's extent projects onto every world axis at once.
fn frame_bounds(centre: Point, frame: ogeom_math::Frame, extent: (f64, f64, f64)) -> Aabb {
    let (ex, ey, ez) = extent;
    let reach = |axis: fn(&Vector) -> f64| {
        (frame.x().vector().pipe(axis) * ex).abs()
            + (frame.y().vector().pipe(axis) * ey).abs()
            + (frame.z().vector().pipe(axis) * ez).abs()
    };
    let r = Vector::new(reach(|v| v.x), reach(|v| v.y), reach(|v| v.z));
    Aabb::of_corners(centre - r, centre + r)
}

/// A tiny helper so the reach computation above reads as one expression.
trait Pipe {
    fn pipe<R>(&self, f: impl FnOnce(&Self) -> R) -> R;
}

impl Pipe for Vector {
    fn pipe<R>(&self, f: impl FnOnce(&Self) -> R) -> R {
        f(self)
    }
}

/// A guaranteed bound for a shape, including everything below it.
///
/// Vertices contribute their point widened by their own tolerance, since a
/// vertex genuinely occupies that much space. Edges and faces contribute the
/// bound of their geometry, likewise widened.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if any handle fails to
/// resolve, and whatever the geometry's own bound reports.
pub fn shape_bounds(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Aabb> {
    let Some(node) = model.node(shape) else {
        ogeom_bail!(Dangling, "shape refers to a node not in this model");
    };
    let placement = shape.transform(model.datums())?;

    let own = match node.data() {
        NodeData::Vertex(v) => Aabb::of_point(placement.apply(v.point)).expanded(v.tolerance.get()),
        NodeData::Edge(e) => {
            let mut out = Aabb::EMPTY;
            for repr in &e.representations {
                if let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = repr
                    && let Some(geometry) = model.geometry().curve(*curve)
                {
                    out = out.union(&curve_bounds(geometry, tol)?);
                }
            }
            out.transformed(&placement).expanded(e.tolerance.get())
        }
        NodeData::Face(f) => {
            // A planar face lies inside its boundary's own hull, so its
            // wires below say everything — an imported plane's carrier
            // window can span kilometres and would drown every consumer. A
            // curved face can bulge past its boundary — a dome past its
            // equator — so those keep the whole surface's bound.
            let planar = matches!(
                model.geometry().surface(f.surface),
                Some(ogeom_geom::SurfaceGeometry::Plane(_))
            );
            if planar && !model.children_of(shape)?.is_empty() {
                Aabb::EMPTY
            } else {
                match model.geometry().surface(f.surface) {
                    Some(surface) => surface_bounds(surface, tol)
                        .unwrap_or(Aabb::EMPTY)
                        .transformed(&placement)
                        .expanded(f.tolerance.get()),
                    None => Aabb::EMPTY,
                }
            }
        }
        NodeData::Container => Aabb::EMPTY,
    };

    let mut out = own;
    for child in model.children_of(shape)? {
        out = out.union(&shape_bounds(model, &child, tol)?);
    }
    Ok(out)
}

/// A bound for a shape built only from its vertices.
///
/// Tighter than [`shape_bounds`] for a solid whose faces sit on unbounded
/// surfaces, and *not* a guarantee: a curved edge bulges past its own
/// endpoints. Use it for a quick estimate, never for a rejection test.
///
/// # Errors
///
/// As [`shape_bounds`].
pub fn vertex_bounds(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Aabb> {
    let mut out = Aabb::EMPTY;
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let Some(node) = model.node(&vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        if let Some(data) = node.data().as_vertex() {
            let placed = vertex.transform(model.datums())?.apply(data.point);
            out = out.with_point(placed);
        }
    }
    Ok(out.expanded(tol.confusion()))
}

/// A box that has been turned to fit what it bounds.
///
/// An axis-aligned box around a long thin rod lying diagonally is mostly empty;
/// this one is not. The cost is that testing a point against it is a transform
/// and then a comparison, rather than six comparisons — so [`Aabb`] stays the
/// default and this is for when the emptiness matters, which is broad-phase
/// rejection and anything that quotes a shape's real extent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Obb {
    /// The centre, and the axes the half-extents are measured along.
    pub frame: Frame,
    /// Half the size along the frame's `x`, `y` and `z`.
    pub half_extent: Vector,
}

impl Obb {
    /// The eight corners, in the same order [`Aabb::corners`] uses.
    #[must_use]
    pub fn corners(&self) -> Vec<Point> {
        let (x, y, z) = (
            self.frame.x().vector() * self.half_extent.x,
            self.frame.y().vector() * self.half_extent.y,
            self.frame.z().vector() * self.half_extent.z,
        );
        let mut out = Vec::with_capacity(8);
        for k in [-1.0_f64, 1.0] {
            for j in [-1.0_f64, 1.0] {
                for i in [-1.0_f64, 1.0] {
                    out.push(self.frame.origin() + x * i + y * j + z * k);
                }
            }
        }
        out
    }

    /// The volume it encloses.
    #[must_use]
    pub fn volume(&self) -> f64 {
        8.0 * self.half_extent.x * self.half_extent.y * self.half_extent.z
    }

    /// Whether a point is inside, measured in the box's own frame.
    #[must_use]
    pub fn contains(&self, p: Point) -> bool {
        let local = self.frame.to_local(p);
        local.x.abs() <= self.half_extent.x
            && local.y.abs() <= self.half_extent.y
            && local.z.abs() <= self.half_extent.z
    }

    /// The axis-aligned box that contains this one.
    #[must_use]
    pub fn to_aabb(&self) -> Aabb {
        Aabb::of_points(&self.corners())
    }
}

/// An oriented bound for a shape, from the spread of its geometry.
///
/// The axes come from the covariance of sampled points — the directions the
/// shape is most and least spread along — and the extents are then measured
/// along those axes, so the box is tight even though the fit is not exact.
///
/// **Not a guarantee, unlike [`shape_bounds`].** It is built from samples, so a
/// curved face can bulge a little past it. Widen it before using it to reject
/// anything.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve; [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape has no geometry to bound.
pub fn oriented_bounds(
    model: &Model,
    shape: &Shape,
    deflection: ogeom_mesh::Deflection,
    tol: Tolerances,
) -> OgeomResult<Obb> {
    let mut points = Vec::new();
    // The tessellation, where there is one to build: it follows the shape's
    // real extent, where the vertices alone would miss the bulge of a cylinder.
    if let Ok(mesh) = ogeom_mesh::triangulate(model, shape, deflection, tol) {
        points.extend(mesh.positions.iter().copied());
    }
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
            points.push(vertex.transform(model.datums())?.apply(data.point));
        }
    }
    if points.is_empty() {
        ogeom_bail!(Construction, "the shape has no geometry to bound");
    }

    let frame = spread_frame(&points, tol);
    let mut low = Vector::new(f64::MAX, f64::MAX, f64::MAX);
    let mut high = Vector::new(f64::MIN, f64::MIN, f64::MIN);
    for p in &points {
        let local = frame.to_local(*p);
        low = Vector::new(low.x.min(local.x), low.y.min(local.y), low.z.min(local.z));
        high = Vector::new(
            high.x.max(local.x),
            high.y.max(local.y),
            high.z.max(local.z),
        );
    }
    // The covariance frame is centred on the mean, which is not the middle of
    // the extent — a shape with more detail at one end pulls it. Recentring is
    // what makes the half-extents symmetric and the box actually tight.
    let middle = (low + high) * 0.5;
    let centre = frame.to_world(Point::ORIGIN + middle);
    Ok(Obb {
        frame: frame.with_origin(centre),
        half_extent: (high - low) * 0.5,
    })
}

/// The direction a face presents, in space.
///
/// Sampled at the mean of its boundary in parameter space, which for a planar
/// profile is exact everywhere and for a curved one is representative: a
/// profile whose normal turns past perpendicular to the sweep somewhere across
/// its own extent sweeps into a solid that passes through itself, and one
/// sample is enough to decide which side the material lands on in every case
/// this can build. A face with no boundary at all covers its whole surface, so
/// the middle of the domain is the point to ask about.
pub fn face_normal(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<(Point, Vector)> {
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(data) = node.data().as_face() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };

    let mut sum = (0.0, 0.0);
    let mut count = 0_u32;
    // The outer wire is the first, and it alone bounds the region; a hole would
    // only pull the sample towards a point the face does not cover.
    for edge in match model.children_of(face)?.first() {
        Some(outer) => model.children_of(outer)?,
        None => Vec::new(),
    } {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let (id, range) = match edge_data.pcurve_for(data.surface, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            Some(EdgeRepr::Seam { forward, range, .. }) => (*forward, *range),
            _ => continue,
        };
        let Some(pcurve) = model.geometry().pcurve(id) else {
            ogeom_bail!(Dangling, "pcurve is not in this model");
        };
        for at in [range.0, f64::midpoint(range.0, range.1), range.1] {
            let p = pcurve.point_at(at, tol)?;
            sum = (sum.0 + p.x, sum.1 + p.y);
            count += 1;
        }
    }

    let ((ua, ub), (va, vb)) = surface.domain();
    let (u, v) = if count == 0 {
        (f64::midpoint(ua, ub), f64::midpoint(va, vb))
    } else {
        let n = f64::from(count);
        (sum.0 / n, sum.1 / n)
    };
    let normal = surface.normal_at(u, v, tol)?;
    let point = surface.point_at(u, v, tol)?;

    let placement = face.transform(model.datums())?;
    let placed = placement.apply_vector(normal.vector());
    Ok((
        placement.apply(point),
        if face.orientation() == Orientation::Reversed {
            -placed
        } else {
            placed
        },
    ))
}

/// A deflection expressed as a fraction of a shape's own size.
///
/// "A thousandth of the part" survives the part being modelled in metres rather
/// than millimetres, and being scaled after it was drawn; an absolute chord
/// does not. [`Deflection::relative`](ogeom_mesh::Deflection::relative) does the
/// arithmetic once a size is known — this is what finds the size, which is the
/// part a caller should not have to get right.
///
/// The measure is the bounding box's *diagonal*, not its longest side: a thin
/// plate and a cube of the same longest side are not equally demanding, and the
/// diagonal is the one that notices.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `fraction` is
/// not finite and positive, or the shape has nothing a deflection would
/// describe.
pub fn relative_deflection(
    model: &Model,
    shape: &Shape,
    fraction: f64,
    tol: Tolerances,
) -> OgeomResult<ogeom_mesh::Deflection> {
    // A deflection says how closely a polyline should follow a curve, or a
    // triangle a surface. A shape with neither has nothing for it to be about —
    // and it is not enough to look at the size, because a lone vertex *does*
    // have a bound: its own tolerance. A fraction of that is a chord of about
    // 1e-10, which is not a small answer, it is a meaningless one.
    if explore_unique(model, shape, ShapeType::Edge)?.is_empty()
        && explore_unique(model, shape, ShapeType::Face)?.is_empty()
    {
        ogeom_bail!(
            Construction,
            "the shape has no edges or faces, so there is nothing a deflection \
             would describe"
        );
    }
    let diagonal = shape_bounds(model, shape, tol)?.diagonal();
    if !diagonal.is_finite() || diagonal <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the shape has no extent for a deflection to be a fraction of"
        );
    }
    ogeom_mesh::Deflection::relative(diagonal, fraction)
}

/// A frame whose axes are the directions a point set is most spread along.
///
/// The eigenvectors of the covariance, largest spread first. Falls back to the
/// world frame when the points are too few or too degenerate to say — a fit
/// that cannot decide should return something usable rather than fail, since
/// the caller then measures extents along whatever axes it gets and still
/// bounds the shape.
fn spread_frame(points: &[Point], tol: Tolerances) -> Frame {
    let Some((centroid, axes)) = covariance_axes(points) else {
        return Frame::WORLD.with_origin(points.first().copied().unwrap_or(Point::ORIGIN));
    };
    // `Frame::new` takes the primary direction as `z` and a *reference* for
    // `x`, so the most-spread axis goes in the second slot: the frame's `x` is
    // the direction the shape is longest along, which is what a caller reading
    // `half_extent.x` will expect, and its `z` is the flattest.
    let [most, _, least] = axes;
    Frame::new(centroid, least, most, tol)
        .or_else(|_| Frame::new(centroid, least, Direction::X, tol))
        .or_else(|_| Frame::new(centroid, least, Direction::Y, tol))
        .unwrap_or_else(|_| Frame::WORLD.with_origin(centroid))
}

/// The plane that best fits a point set: its centroid, and the normal to it.
///
/// The eigenvector of the covariance with the *smallest* eigenvalue — the
/// direction the points vary along least. `None` when there is nothing to fit.
///
/// Fitting says nothing about whether the points are actually planar. Every set
/// of three or more points has a best-fit plane, including a set that is
/// nowhere near one, so a caller has to measure the residual and decide.
pub(crate) fn least_squares_plane(points: &[Point], tol: Tolerances) -> Option<(Point, Direction)> {
    let (centroid, axes) = covariance_axes(points)?;
    let _ = tol;
    Some((centroid, axes[2]))
}

/// The centroid of a point set and its covariance eigenvectors, most spread
/// first.
fn covariance_axes(points: &[Point]) -> Option<(Point, [Direction; 3])> {
    if points.len() < 3 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let n = points.len() as f64;
    let mut sum = Vector::ZERO;
    for p in points {
        sum += p.to_vector();
    }
    let centroid = Point::ORIGIN + sum * (1.0 / n);

    let mut c = nalgebra::Matrix3::<f64>::zeros();
    for p in points {
        let d = *p - centroid;
        let v = nalgebra::Vector3::new(d.x, d.y, d.z);
        c += v * v.transpose();
    }
    c /= n;

    // Symmetric by construction, so the eigenvalues are real and the
    // eigenvectors orthogonal.
    let eigen = nalgebra::SymmetricEigen::new(c);
    let mut order: Vec<usize> = (0..3).collect();
    order.sort_by(|a, b| {
        eigen.eigenvalues[*b]
            .partial_cmp(&eigen.eigenvalues[*a])
            .unwrap_or(core::cmp::Ordering::Equal)
    });

    let mut axes = Vec::with_capacity(3);
    for i in order {
        let column = eigen.eigenvectors.column(i);
        axes.push(
            Direction::from_coords(column[0], column[1], column[2], Tolerances::millimetres())
                .ok()?,
        );
    }
    Some((centroid, [axes[0], axes[1], axes[2]]))
}

/// Where a point projects onto a curve, and how far away it is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Projection {
    /// The parameter of the nearest point found.
    pub parameter: f64,
    /// The nearest point.
    pub point: Point,
    /// The distance to it.
    pub distance: f64,
}

/// The nearest point on a curve to `target`.
///
/// Samples the domain to bracket the minimum, then refines. The sampling is not
/// decoration: the distance function along a curve is generally multi-modal —
/// a point inside a circle is equidistant from every part of it, and a point
/// near a spline's inflection has two competing minima — so starting a local
/// method from one guess finds whichever basin it happens to land in. The
/// sample count sets how fine a feature can be resolved, and is stated rather
/// than hidden.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the curve cannot be
/// evaluated over its own domain.
pub fn project_on_curve(
    curve: &Curve,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgeomResult<Projection> {
    let (a, b) = curve.domain();
    let steps = samples.max(8);

    let distance_at = |u: f64| -> f64 {
        curve
            .point_at(u, tol)
            .map_or(f64::INFINITY, |p| p.square_distance(target))
    };

    // Coarse scan for the best bracket.
    let mut best = (a, distance_at(a));
    let mut best_index = 0_usize;
    for i in 1..=steps {
        #[allow(clippy::cast_precision_loss)]
        let u = a + (b - a) * (i as f64 / steps as f64);
        let d = distance_at(u);
        if d < best.1 {
            best = (u, d);
            best_index = i;
        }
    }

    // Refine inside the neighbouring samples, where the minimum must lie.
    #[allow(clippy::cast_precision_loss)]
    let width = (b - a) / steps as f64;
    let lo = (best.0 - width).max(a);
    let hi = (best.0 + width).min(b);
    let _ = best_index;

    let parameter = if hi > lo {
        let refined = solve::minimize(
            distance_at,
            lo,
            hi,
            solve::Criteria {
                residual: 0.0,
                step: tol.parametric(),
                max_iterations: 100,
            },
        )?;
        // The refinement may land marginally worse than the sample if the
        // bracket was already at the boundary; keep whichever is actually
        // nearer rather than trusting the method.
        if distance_at(refined.value) <= best.1 {
            refined.value
        } else {
            best.0
        }
    } else {
        best.0
    };

    let point = curve.point_at(parameter, tol)?;
    Ok(Projection {
        parameter,
        point,
        distance: point.distance(target),
    })
}

/// Where a point projects onto a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceProjection {
    /// The parameters of the nearest point found.
    pub parameters: (f64, f64),
    /// The nearest point.
    pub point: Point,
    /// The distance to it.
    pub distance: f64,
}

/// The nearest point on a surface to `target`.
///
/// A coarse grid to bracket, then Newton on the two conditions that define a
/// foot point: the displacement from the surface to the target is perpendicular
/// to both tangents. Grid resolution is `samples` per direction, for the same
/// reason as [`project_on_curve`].
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the surface cannot be
/// evaluated over its own domain.
pub fn project_on_surface(
    surface: &SurfaceGeometry,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgeomResult<SurfaceProjection> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let steps = samples.max(4);

    let mut best = (ua, va, f64::INFINITY);
    for i in 0..=steps {
        for j in 0..=steps {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                ua + (ub - ua) * (i as f64 / steps as f64),
                va + (vb - va) * (j as f64 / steps as f64),
            );
            if let Ok(p) = surface.point_at(u, v, tol) {
                let d = p.square_distance(target);
                if d < best.2 {
                    best = (u, v, d);
                }
            }
        }
    }

    // The foot point conditions: (S - target) . Su = 0 and (S - target) . Sv = 0.
    let residual = |x: &[f64]| {
        let (u, v) = (x[0], x[1]);
        let Ok(p) = surface.point_at(u, v, tol) else {
            return (vec![0.0, 0.0], vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        };
        let Ok((du, dv)) = surface.d1_at(u, v, tol) else {
            return (vec![0.0, 0.0], vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        };
        let (d2u, duv, d2v) =
            surface
                .d2_at(u, v, tol)
                .unwrap_or((Vector::ZERO, Vector::ZERO, Vector::ZERO));
        let gap = p - target;
        (
            vec![gap.dot(du), gap.dot(dv)],
            vec![
                vec![du.dot(du) + gap.dot(d2u), du.dot(dv) + gap.dot(duv)],
                vec![du.dot(dv) + gap.dot(duv), dv.dot(dv) + gap.dot(d2v)],
            ],
        )
    };

    let refined = solve::newton_system(
        residual,
        &[best.0, best.1],
        solve::Criteria {
            residual: tol.confusion(),
            step: tol.parametric(),
            max_iterations: 60,
        },
    );

    let (u, v) = match refined {
        Ok(solution) if solution.convergence.is_converged() => {
            let (u, v) = (solution.value[0], solution.value[1]);
            // Newton is free to wander outside the domain; a foot point that
            // left it is not a foot point of this surface.
            match surface.normalize_parameters(u, v, tol) {
                Ok(inside)
                    if surface
                        .point_at(inside.0, inside.1, tol)
                        .is_ok_and(|p| p.square_distance(target) <= best.2) =>
                {
                    inside
                }
                _ => (best.0, best.1),
            }
        }
        _ => (best.0, best.1),
    };

    let point = surface.point_at(u, v, tol)?;
    Ok(SurfaceProjection {
        parameters: (u, v),
        point,
        distance: point.distance(target),
    })
}

/// The nearest point on a planar curve to a point in the same parameter space.
///
/// # Errors
///
/// As [`project_on_curve`].
pub fn project_on_planar_curve(
    curve: &PlanarCurve,
    target: Point2,
    samples: usize,
    tol: Tolerances,
) -> OgeomResult<(f64, Point2, f64)> {
    use ogeom_geom::Curve2d;

    let (a, b) = curve.domain();
    let steps = samples.max(8);
    let distance_at = |u: f64| -> f64 {
        curve
            .point_at(u, tol)
            .map_or(f64::INFINITY, |p| p.square_distance(target))
    };

    let mut best = (a, distance_at(a));
    for i in 1..=steps {
        #[allow(clippy::cast_precision_loss)]
        let u = a + (b - a) * (i as f64 / steps as f64);
        let d = distance_at(u);
        if d < best.1 {
            best = (u, d);
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let width = (b - a) / steps as f64;
    let (lo, hi) = ((best.0 - width).max(a), (best.0 + width).min(b));
    let parameter = if hi > lo {
        let refined = solve::minimize(
            distance_at,
            lo,
            hi,
            solve::Criteria {
                residual: 0.0,
                step: tol.parametric(),
                max_iterations: 100,
            },
        )?;
        if distance_at(refined.value) <= best.1 {
            refined.value
        } else {
            best.0
        }
    } else {
        best.0
    };

    let point = curve.point_at(parameter, tol)?;
    Ok((parameter, point, point.distance(target)))
}

/// A surface's parameterization window widened until it holds every one of
/// `points`.
///
/// A surface's extent is a *window*, not a trim. Anything built on the surface
/// — a face, a pcurve, a projection — has to evaluate inside it, and a window
/// clamped tight around whatever was measured last will refuse the boundary of
/// the very region it was measured from. That failure is unhelpfully quiet: it
/// arrives as a domain error from an evaluation deep inside triangulation,
/// having overshot by a part in ten million.
///
/// Widening is therefore a step to take *before* building on a surface whose
/// window came from samples, and it changes nothing about the geometry: the
/// carrier is untouched and only the window moves. The margin is proportional
/// to the span measured, plus a floor in confusion tolerances, so widening is
/// not itself a tolerance question.
///
/// Only the *bounded* directions can be widened, and only they need to be: a
/// periodic direction already covers its whole turn. A surface with no bounded
/// direction — a sphere, a torus — comes back as it went in, and so does one
/// given no points.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a point cannot be
/// projected onto the surface, or the widened window is not a valid range.
pub fn widened_to_hold(
    surface: &SurfaceGeometry,
    points: &[Point],
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    use ogeom_geom::SurfaceGeometry as S;
    use ogeom_geom::{ConeSurface, CylinderSurface, PlaneSurface};
    if points.is_empty() {
        return Ok(surface.clone());
    }
    let (mut u0, mut u1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in points {
        let projected = project_on_surface(surface, *p, 16, tol)?;
        let (u, v) = projected.parameters;
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    if !u0.is_finite() || !v0.is_finite() {
        return Ok(surface.clone());
    }
    // A margin proportional to what was measured, so widening is not itself
    // a tolerance question.
    let margin = |lo: f64, hi: f64| (hi - lo).mul_add(0.05, tol.confusion() * 1e3);
    let ((du0, du1), (dv0, dv1)) = surface.domain();
    Ok(match surface {
        S::Plane(p) => {
            let (mu, mv) = (margin(u0, u1), margin(v0, v1));
            PlaneSurface::over(
                p.plane(),
                (du0.min(u0 - mu), du1.max(u1 + mu)),
                (dv0.min(v0 - mv), dv1.max(v1 + mv)),
            )?
            .into()
        }
        S::Cylinder(c) => {
            let m = margin(v0, v1);
            CylinderSurface::new(c.cylinder(), (dv0.min(v0 - m), dv1.max(v1 + m)))?.into()
        }
        S::Cone(c) => {
            let m = margin(v0, v1);
            ConeSurface::new(c.cone(), (dv0.min(v0 - m), dv1.max(v1 + m)))?.into()
        }
        other => other.clone(),
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::make_box;
    use approx::assert_relative_eq;
    use ogeom_geom::{
        BSplineCurve, CircleCurve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface,
        TorusSurface, TrimmedCurve,
    };
    use ogeom_math::{Circle, Cylinder, Direction, Frame, KnotVector, Plane, Sphere, Torus};

    const T: Tolerances = Tolerances::millimetres();

    /// A curve sampled densely — the ground truth a bound must contain.
    fn dense_points(curve: &Curve, n: usize) -> Vec<Point> {
        let (a, b) = curve.domain();
        (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let u = a + (b - a) * (i as f64 / n as f64);
                curve.point_at(u, T).unwrap()
            })
            .collect()
    }

    #[test]
    fn a_line_bounds_exactly_to_its_endpoints() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
            .unwrap()
            .into();
        let b = curve_bounds(&curve, T).unwrap();
        assert_eq!(b.low(), Some(Point::ORIGIN));
        assert_eq!(b.high(), Some(Point::new(3.0, 4.0, 0.0)));
    }

    #[test]
    fn every_curves_bound_contains_the_curve() {
        // The one property that matters. Checked against dense sampling, which
        // is fine as a *test* oracle even though it is not sound as an
        // implementation.
        let spline_control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 5.0, 0.0),
            Point::new(3.0, -4.0, 2.0),
            Point::new(5.0, 2.0, -1.0),
            Point::new(6.0, 0.0, 0.0),
        ];
        let tilted = Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();

        let curves: Vec<Curve> = vec![
            LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
                .unwrap()
                .into(),
            CircleCurve::new(Circle::new(tilted, 2.0, T).unwrap()).into(),
            ogeom_geom::EllipseCurve::new(ogeom_math::Ellipse::new(tilted, 5.0, 3.0, T).unwrap())
                .into(),
            ogeom_geom::HyperbolaCurve::new(
                ogeom_math::Hyperbola::new(tilted, 3.0, 4.0, T).unwrap(),
                1.5,
            )
            .unwrap()
            .into(),
            ogeom_geom::ParabolaCurve::new(ogeom_math::Parabola::new(tilted, 2.0, T).unwrap(), 4.0)
                .unwrap()
                .into(),
            BSplineCurve::new(
                KnotVector::clamped_uniform(3, spline_control.len()).unwrap(),
                spline_control,
                T,
            )
            .unwrap()
            .into(),
        ];

        for curve in curves {
            let bound = curve_bounds(&curve, T).unwrap().with_tolerance(T);
            for p in dense_points(&curve, 400) {
                assert!(
                    bound.contains(p),
                    "{:?} escaped its bound at {p:?}: {bound}",
                    curve.kind()
                );
            }
        }
    }

    #[test]
    fn a_splines_bound_is_its_control_hull_and_that_is_a_guarantee() {
        // Sampling would miss the bulge between samples; the convex hull
        // property does not, because the curve provably never leaves the hull.
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 10.0, 0.0),
            Point::new(2.0, 10.0, 0.0),
            Point::new(3.0, 0.0, 0.0),
        ];
        let curve: Curve = BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control.clone(),
            T,
        )
        .unwrap()
        .into();

        let bound = curve_bounds(&curve, T).unwrap();
        assert_eq!(bound, Aabb::of_points(&control));
        for p in dense_points(&curve, 200) {
            assert!(bound.contains(p));
        }
        // And the curve really does stay well inside — the bound is loose, in
        // the safe direction.
        let peak = dense_points(&curve, 200)
            .iter()
            .fold(0.0_f64, |m, p| m.max(p.y));
        assert!(peak < 10.0, "the curve should not reach its control points");
    }

    #[test]
    fn a_trimmed_curve_reports_the_bound_of_the_whole() {
        // Loose but never wrong. Tightening it means solving for the extremes
        // of the trimmed range, which is the same work as an intersection.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let quarter: Curve = TrimmedCurve::new(circle.clone(), 0.0, 1.5, T)
            .unwrap()
            .into();

        let whole = curve_bounds(&circle, T).unwrap();
        let part = curve_bounds(&quarter, T).unwrap();
        assert_eq!(part, whole);
        for p in dense_points(&quarter, 200) {
            assert!(part.contains(p));
        }
    }

    #[test]
    fn every_surfaces_bound_contains_the_surface() {
        let tilted = Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let surfaces: Vec<SurfaceGeometry> = vec![
            PlaneSurface::over(Plane::new(tilted), (-5.0, 5.0), (-3.0, 3.0))
                .unwrap()
                .into(),
            CylinderSurface::new(Cylinder::new(tilted, 2.0, T).unwrap(), (-4.0, 4.0))
                .unwrap()
                .into(),
            ogeom_geom::ConeSurface::new(
                ogeom_math::Cone::new(tilted, 3.0, 0.6, T).unwrap(),
                (-1.0, 5.0),
            )
            .unwrap()
            .into(),
            SphereSurface::new(Sphere::new(tilted, 4.0, T).unwrap()).into(),
            TorusSurface::new(Torus::new(tilted, 5.0, 2.0, T).unwrap()).into(),
        ];

        for surface in surfaces {
            let bound = surface_bounds(&surface, T).unwrap().with_tolerance(T);
            let ((ua, ub), (va, vb)) = surface.domain();
            for i in 0..=40 {
                for j in 0..=40 {
                    let u = ua + (ub - ua) * (f64::from(i) / 40.0);
                    let v = va + (vb - va) * (f64::from(j) / 40.0);
                    let p = surface.point_at(u, v, T).unwrap();
                    assert!(
                        bound.contains(p),
                        "{:?} escaped its bound at ({u}, {v}) -> {p:?}: {bound}",
                        surface.kind()
                    );
                }
            }
        }
    }

    #[test]
    fn a_spheres_bound_is_exact() {
        let s: SurfaceGeometry =
            SphereSurface::new(Sphere::centred(Point::new(1.0, 2.0, 3.0), 4.0, T).unwrap()).into();
        let b = surface_bounds(&s, T).unwrap();
        assert_eq!(b.low(), Some(Point::new(-3.0, -2.0, -1.0)));
        assert_eq!(b.high(), Some(Point::new(5.0, 6.0, 7.0)));
    }

    #[test]
    fn an_unbounded_plane_is_refused_rather_than_bounded_wrongly() {
        // Reporting the enormous default extent as a bound would make every
        // rejection test involving it useless, and silently.
        let s: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        assert!(surface_bounds(&s, T).is_err());
    }

    #[test]
    fn a_shapes_bound_contains_every_vertex_it_holds() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let bound = shape_bounds(&model, &built.shape, T).unwrap();

        for vertex in explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap() {
            let p = model
                .node(&vertex)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            assert!(bound.contains(p), "vertex {p:?} escaped {bound}");
        }
    }

    #[test]
    fn the_vertex_bound_of_a_box_is_tight_and_the_full_bound_contains_it() {
        // A box's faces sit on planes trimmed to the box, so the two agree
        // closely here — but the vertex bound is documented as an estimate, and
        // the full bound is the one that guarantees containment.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();

        let tight = vertex_bounds(&model, &built.shape, T).unwrap();
        assert_relative_eq!(tight.size().x, 2.0, epsilon = 1e-6);
        assert_relative_eq!(tight.size().y, 3.0, epsilon = 1e-6);
        assert_relative_eq!(tight.size().z, 4.0, epsilon = 1e-6);

        let full = shape_bounds(&model, &built.shape, T).unwrap();
        assert!(full.contains_box(&tight));
    }

    #[test]
    fn a_placed_shapes_bound_moves_with_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let here = vertex_bounds(&model, &built.shape, T).unwrap();

        let moved = model.placed(
            &built.shape,
            ogeom_math::Transform::translation(Vector::new(10.0, 0.0, 0.0)),
        );
        let there = vertex_bounds(&model, &moved, T).unwrap();

        assert_relative_eq!(
            there.centre().unwrap().x - here.centre().unwrap().x,
            10.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(there.size().x, here.size().x, epsilon = 1e-9);
    }

    #[test]
    fn projecting_onto_a_line_lands_on_the_foot_of_the_perpendicular() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let p = project_on_curve(&curve, Point::new(3.0, 4.0, 0.0), 32, T).unwrap();
        assert_relative_eq!(p.parameter, 3.0, epsilon = 1e-6);
        assert!(p.point.is_equal(Point::new(3.0, 0.0, 0.0), T));
        assert_relative_eq!(p.distance, 4.0, epsilon = 1e-9);
    }

    #[test]
    fn a_window_widened_to_hold_a_point_holds_it_with_room_to_spare() {
        // The failure this exists to stop: a window clamped exactly to what was
        // measured refuses the boundary of the region it was measured from,
        // by a part in ten million, as a domain error from deep inside
        // whatever was being built on it.
        let cylinder = Cylinder::new(Frame::WORLD, 5.0, T).unwrap();
        let tight: SurfaceGeometry = CylinderSurface::new(cylinder, (0.0, 10.0)).unwrap().into();
        let just_past = Point::new(5.0, 0.0, 10.0 + 1e-6);
        assert!(
            tight.domain().1.1 < just_past.z,
            "the point is outside the tight window, which is the premise"
        );

        let wide = widened_to_hold(&tight, &[just_past], T).unwrap();
        let (_, (v0, v1)) = wide.domain();
        assert!(
            v1 > just_past.z && v0 <= 0.0,
            "the window holds the point and gives nothing back: ({v0}, {v1})"
        );
        // The window grows by the floor — a thousand confusions — and this is
        // the case that says why there is a floor at all. Projection clamps to
        // the window it is measuring, so a point that overshoots by 1e-6 comes
        // back at the old edge and the proportional term sees a span of zero.
        // Only the floor stands between the overshoot and another domain
        // error, which is why it is a thousand confusions and not one.
        assert_relative_eq!(v1 - 10.0, T.confusion() * 1e3, epsilon = 1e-12);
        assert!(
            v1 - just_past.z > 0.0,
            "and it clears the point: {}",
            v1 - just_past.z
        );

        // The carrier is untouched. Widening a window is not a change of shape.
        let SurfaceGeometry::Cylinder(c) = &wide else {
            panic!("still a cylinder");
        };
        assert_relative_eq!(c.cylinder().radius(), 5.0, epsilon = 1e-15);
    }

    #[test]
    fn a_surface_with_no_bounded_direction_comes_back_unchanged() {
        // A sphere is periodic one way and bounded by its own poles the other:
        // there is no window to widen, and inventing one would put chart space
        // past the pole where the parameterization means nothing.
        let sphere: SurfaceGeometry =
            SphereSurface::new(Sphere::new(Frame::WORLD, 4.0, T).unwrap()).into();
        let widened = widened_to_hold(&sphere, &[Point::new(4.0, 0.0, 0.0)], T).unwrap();
        assert_eq!(sphere.domain(), widened.domain());

        // And no points is no information, so nothing moves either.
        let cylinder: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.0, T).unwrap(), (1.0, 3.0))
                .unwrap()
                .into();
        assert_eq!(
            cylinder.domain(),
            widened_to_hold(&cylinder, &[], T).unwrap().domain()
        );
    }

    #[test]
    fn projecting_past_the_end_of_a_segment_clamps_to_the_end() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let p = project_on_curve(&curve, Point::new(50.0, 0.0, 0.0), 32, T).unwrap();
        assert_relative_eq!(p.parameter, 10.0, epsilon = 1e-6);
        assert_relative_eq!(p.distance, 40.0, epsilon = 1e-6);
    }

    #[test]
    fn projecting_onto_a_circle_finds_the_nearest_of_many_minima() {
        // The reason for the coarse scan: from outside the circle's plane there
        // is one minimum, but a local method started at the wrong parameter
        // converges to the far side just as happily.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 5.0, T).unwrap()).into();
        for angle in [0.1_f64, 1.0, 2.5, 4.0, 6.0] {
            let outside = Point::new(8.0 * angle.cos(), 8.0 * angle.sin(), 0.0);
            let p = project_on_curve(&circle, outside, 64, T).unwrap();
            assert_relative_eq!(p.distance, 3.0, epsilon = 1e-6);
            assert_relative_eq!(p.point.to_vector().magnitude(), 5.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn projecting_onto_a_plane_gives_the_perpendicular_foot() {
        let plane: SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (-10.0, 10.0), (-10.0, 10.0))
                .unwrap()
                .into();
        let p = project_on_surface(&plane, Point::new(2.0, 3.0, 7.0), 8, T).unwrap();
        assert!(p.point.is_equal(Point::new(2.0, 3.0, 0.0), T));
        assert_relative_eq!(p.distance, 7.0, epsilon = 1e-9);
    }

    #[test]
    fn projecting_onto_a_sphere_lands_on_the_radial_line() {
        let sphere = Sphere::centred(Point::new(1.0, 1.0, 1.0), 3.0, T).unwrap();
        let surface: SurfaceGeometry = SphereSurface::new(sphere).into();
        for target in [
            Point::new(10.0, 1.0, 1.0),
            Point::new(1.0, 1.0, 9.0),
            Point::new(-4.0, -2.0, 0.0),
        ] {
            let p = project_on_surface(&surface, target, 16, T).unwrap();
            // The foot point is on the sphere, and on the line from the centre.
            assert_relative_eq!(sphere.centre().distance(p.point), 3.0, max_relative = 1e-7);
            assert_relative_eq!(
                p.distance,
                (sphere.centre().distance(target) - 3.0).abs(),
                max_relative = 1e-6
            );
        }
    }

    #[test]
    fn projecting_onto_a_cylinder_is_radial() {
        let cylinder = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        let surface: SurfaceGeometry = CylinderSurface::new(cylinder, (-5.0, 5.0)).unwrap().into();
        let p = project_on_surface(&surface, Point::new(6.0, 0.0, 1.0), 16, T).unwrap();
        assert_relative_eq!(p.distance, 4.0, max_relative = 1e-6);
        assert_relative_eq!(p.point.z, 1.0, epsilon = 1e-6);
        assert_relative_eq!(p.point.x.hypot(p.point.y), 2.0, max_relative = 1e-7);
    }

    #[test]
    fn projection_of_a_point_already_on_the_geometry_returns_zero_distance() {
        let curve: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()).into();
        let on = curve.point_at(1.2, T).unwrap();
        let p = project_on_curve(&curve, on, 64, T).unwrap();
        assert!(p.distance < 1e-7, "distance was {}", p.distance);
    }

    #[test]
    fn projecting_onto_a_planar_curve_works_in_parameter_space() {
        let curve: PlanarCurve =
            ogeom_geom::Line2d::segment(Point2::ORIGIN, Point2::new(10.0, 0.0), T)
                .unwrap()
                .into();
        let (u, point, distance) =
            project_on_planar_curve(&curve, Point2::new(3.0, 4.0), 32, T).unwrap();
        assert_relative_eq!(u, 3.0, epsilon = 1e-6);
        assert!(point.is_equal(Point2::new(3.0, 0.0), T));
        assert_relative_eq!(distance, 4.0, epsilon = 1e-9);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod oriented_bound_tests {
    use super::*;
    use crate::{make_box, make_cylinder};
    use approx::assert_relative_eq;
    use ogeom_math::Transform;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> ogeom_mesh::Deflection {
        ogeom_mesh::Deflection {
            chord: 0.01,
            ..ogeom_mesh::Deflection::default()
        }
    }

    #[test]
    fn an_oriented_box_around_a_box_is_that_box() {
        let mut model = Model::new();
        let size = (2.0, 5.0, 1.0);
        let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
        let obb = oriented_bounds(&model, &built.shape, fine(), T).unwrap();

        assert_relative_eq!(obb.volume(), size.0 * size.1 * size.2, epsilon = 1e-9);
        assert!(
            obb.frame.origin().distance(Point::new(1.0, 2.5, 0.5)) < 1e-9,
            "got {:?}",
            obb.frame.origin()
        );
        // The half-extents are the box's, in some order: the axes come from the
        // spread, which does not know or care which one we called x.
        let mut found = [obb.half_extent.x, obb.half_extent.y, obb.half_extent.z];
        found.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut want = [size.0 / 2.0, size.1 / 2.0, size.2 / 2.0];
        want.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (a, b) in found.iter().zip(&want) {
            assert_relative_eq!(a, b, epsilon = 1e-9);
        }
    }

    #[test]
    fn turning_a_box_turns_its_oriented_bound_with_it_and_not_its_volume() {
        // The whole point. An axis-aligned box around a rotated box grows; an
        // oriented one does not, and that difference is what makes it worth the
        // transform at every containment test.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 6.0, 1.0), T).unwrap();
        let turned = crate::transformed(
            &mut model,
            &built.shape,
            Transform::rotation(
                ogeom_math::Axis::new(Point::ORIGIN, Direction::Z),
                std::f64::consts::FRAC_PI_4,
            ),
        )
        .unwrap()
        .shape;

        let obb = oriented_bounds(&model, &turned, fine(), T).unwrap();
        let aabb = shape_bounds(&model, &turned, T).unwrap();
        assert_relative_eq!(obb.volume(), 6.0, max_relative = 1e-6);
        assert!(
            aabb.volume() > obb.volume() * 1.5,
            "an axis-aligned box around a diagonal bar should be much emptier: \
             {} against {}",
            aabb.volume(),
            obb.volume()
        );
        for corner in obb.corners() {
            assert!(obb.contains(corner) || obb.to_aabb().contains(corner));
        }
    }

    #[test]
    fn a_cylinders_oriented_bound_follows_its_axis() {
        let mut model = Model::new();
        let (radius, height) = (0.5_f64, 8.0);
        let built = make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();
        let obb = oriented_bounds(&model, &built.shape, fine(), T).unwrap();

        // The long axis is the cylinder's own, and it is the first the spread
        // reports.
        assert!(
            obb.frame
                .x()
                .vector()
                .cross(Direction::Z.vector())
                .magnitude()
                < 1e-6,
            "the most-spread axis should be the cylinder's, got {:?}",
            obb.frame.x()
        );
        assert_relative_eq!(obb.half_extent.x, height / 2.0, max_relative = 1e-6);
    }

    #[test]
    fn a_shape_with_nothing_to_bound_says_so() {
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);
        // One point has a bound but no spread; it must not claim a frame it
        // cannot justify, and it must not fail either.
        assert!(oriented_bounds(&model, &vertex, fine(), T).is_ok());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod deflection_tests {
    use super::*;
    use crate::make_box;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_relative_deflection_follows_the_shape_it_is_for() {
        // The property that makes it worth having: the same fraction gives the
        // same *number of segments* whatever units the part was drawn in, and
        // an absolute chord does not.
        let mut model = Model::new();
        let small = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let large = make_box(&mut model, Frame::WORLD, (1000.0, 1000.0, 1000.0), T)
            .unwrap()
            .shape;

        let a = relative_deflection(&model, &small, 1e-3, T).unwrap();
        let b = relative_deflection(&model, &large, 1e-3, T).unwrap();
        // Not exactly a thousand: `shape_bounds` is a *guaranteed* bound, so it
        // includes each entity's tolerance, and that padding is a larger share
        // of a one-unit box than of a thousand-unit one. Which is the right
        // behaviour — the padding is really there.
        assert_relative_eq!(b.chord / a.chord, 1000.0, max_relative = 1e-3);
        assert_relative_eq!(a.chord, 3.0_f64.sqrt() * 1e-3, max_relative = 1e-3);
    }

    #[test]
    fn a_shape_with_no_extent_has_no_fraction_of_itself() {
        // A lone vertex *does* have a bound — its own tolerance — so the guard
        // cannot be about size. It is about whether there is a curve or a
        // surface for a deflection to describe.
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);
        let err = relative_deflection(&model, &vertex, 1e-3, T).unwrap_err();
        assert!(
            err.to_string().contains("no edges or faces"),
            "unexpected message: {err}"
        );

        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        assert!(relative_deflection(&model, &solid, 0.0, T).is_err());
        assert!(relative_deflection(&model, &solid, -1.0, T).is_err());
        assert!(relative_deflection(&model, &solid, f64::NAN, T).is_err());
    }

    #[test]
    fn a_planar_face_is_bounded_by_its_wires_not_its_carrier() {
        // An imported plane declares a carrier window of kilometres; the
        // face on it spans a hand's width, and its bound must say so.
        let mut model = Model::new();
        let block = crate::make_box(
            &mut model,
            ogeom_math::Frame::WORLD,
            (8.0, 6.0, 4.0),
            Tolerances::millimetres(),
        )
        .unwrap();
        let bound = shape_bounds(&model, &block.shape, Tolerances::millimetres()).unwrap();
        let (Some(lo), Some(hi)) = (bound.low(), bound.high()) else {
            panic!("the box has a bound");
        };
        assert!(
            lo.distance(ogeom_math::Point::new(0.0, 0.0, 0.0)) < 1e-3,
            "{lo:?}"
        );
        assert!(
            hi.distance(ogeom_math::Point::new(8.0, 6.0, 4.0)) < 1e-3,
            "{hi:?}"
        );
    }
}
