//! Surface/surface intersection: the analytic cases.
//!
//! Where two surfaces meet has a closed form for a specific and well-known set
//! of pairs, and a general answer that needs a marching intersector with a
//! fitting stage after it. This module is the first of those. It is deliberately
//! *not* a partial implementation of the second: a pair it cannot solve exactly
//! is reported as needing the general path, never approximated.
//!
//! # Why the exact cases come first, and separately
//!
//! Three reasons, and the third is the one that matters.
//!
//! They are common. Plane against plane, plane against cylinder, sphere against
//! sphere — a mechanical part is mostly these, and running a marching
//! intersector over a pair whose answer is a circle is slower and less accurate
//! than writing down the circle.
//!
//! They are fast. No stepping, no refinement, no approximation stage.
//!
//! And they are *ground truth*. Every result here can be checked without
//! reference to anything but the two surfaces themselves: sample the curve, ask
//! each surface how far away it is, and the answer should be zero. That check is
//! the instrument the intersection gate is measured with (`docs/SCOPE.md` §7),
//! and it only exists because these cases are exact. A benchmark whose reference
//! answers came from the thing being benchmarked would measure nothing.
//!
//! # What it reports
//!
//! Not just curves. Two surfaces can miss, touch at a point, meet along curves,
//! or be the same surface — and those are four different answers that downstream
//! code has to distinguish. A boolean that treats coincidence as "no
//! intersection" produces a solid with a face missing.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, SurfaceGeometry};
use og_math::{Circle, Direction, Ellipse, Frame, Point, Vector};

/// What two surfaces do where they meet.
#[derive(Debug, Clone, PartialEq)]
pub enum Meeting {
    /// They do not meet at all.
    Apart,
    /// They touch at isolated points, without crossing.
    ///
    /// A sphere resting on a plane. Distinguished from a curve because a
    /// tangential contact has no length to walk along, and an algorithm that
    /// treated it as a degenerate curve would divide by that length.
    Touching(Vec<Point>),
    /// They meet along these curves.
    Along(Vec<Curve>),
    /// They are the same surface wherever they overlap.
    ///
    /// A separate answer from every other, because it is the one where "the
    /// intersection curve" does not exist: the overlap is two-dimensional. A
    /// boolean has to detect this and unify the faces rather than look for a
    /// seam between them.
    Same,
}

/// Where two surfaces meet, when that has a closed form.
///
/// # Errors
///
/// [`OgError::NotDone`](og_core::OgError::NotDone) if this pair has no closed
/// form — which is a statement about the pair, not a failure to compute. The
/// general marching intersector is what answers those, and it is gated on the
/// benchmark this module makes possible.
pub fn surface_surface(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> OgResult<Meeting> {
    use SurfaceGeometry as S;
    match (a, b) {
        (S::Plane(p), S::Plane(q)) => Ok(plane_plane(p.plane(), q.plane(), tol)),
        (S::Plane(p), S::Sphere(s)) => Ok(plane_sphere(p.plane(), s.sphere(), tol)),
        (S::Sphere(s), S::Plane(p)) => Ok(plane_sphere(p.plane(), s.sphere(), tol)),
        (S::Plane(p), S::Cylinder(c)) => plane_cylinder(p.plane(), c.cylinder(), tol),
        (S::Cylinder(c), S::Plane(p)) => plane_cylinder(p.plane(), c.cylinder(), tol),
        (S::Sphere(x), S::Sphere(y)) => Ok(sphere_sphere(x.sphere(), y.sphere(), tol)),
        (S::Cylinder(x), S::Cylinder(y)) => coaxial_cylinders(x.cylinder(), y.cylinder(), tol),
        (S::Cylinder(c), S::Sphere(s)) => coaxial_cylinder_sphere(c.cylinder(), s.sphere(), tol),
        (S::Sphere(s), S::Cylinder(c)) => coaxial_cylinder_sphere(c.cylinder(), s.sphere(), tol),
        _ => og_bail!(
            NotDone,
            "this pair of surfaces has no closed-form intersection; it needs \
             the general marching intersector, which is gated on the benchmark \
             these cases provide the ground truth for"
        ),
    }
}

/// Two planes: apart, the same, or a line.
fn plane_plane(a: og_math::Plane, b: og_math::Plane, tol: Tolerances) -> Meeting {
    let along = a.normal().dot(b.normal());
    if (along.abs() - 1.0).abs() <= tol.angular() {
        // Parallel. Either the same plane or two that never meet, decided by
        // whether one contains the other's origin.
        return if a.distance_to(b.origin()) <= tol.confusion() {
            Meeting::Same
        } else {
            Meeting::Apart
        };
    }
    // The line of intersection runs along both normals' cross product, and
    // passes through the point nearest the origin that satisfies both planes.
    let Ok(direction) = Direction::from_cross(a.normal().vector(), b.normal().vector(), tol) else {
        return Meeting::Apart;
    };
    let (da, db) = (
        a.normal().dot_vector(a.origin().to_vector()),
        b.normal().dot_vector(b.origin().to_vector()),
    );
    let (na, nb) = (a.normal().vector(), b.normal().vector());
    let dot = na.dot(nb);
    let denominator = dot.mul_add(-dot, 1.0);
    if denominator.abs() <= tol.angular() {
        return Meeting::Apart;
    }
    let ca = da.mul_add(1.0, -(db * dot)) / denominator;
    let cb = db.mul_add(1.0, -(da * dot)) / denominator;
    let through = Point::from_vector(na * ca + nb * cb);
    Meeting::Along(vec![line_through(through, direction)])
}

/// A plane and a sphere: apart, a point of tangency, or a circle.
fn plane_sphere(plane: og_math::Plane, sphere: og_math::Sphere, tol: Tolerances) -> Meeting {
    let gap = plane.signed_distance_to(sphere.centre());
    let reach = gap.abs();
    if reach > sphere.radius() + tol.confusion() {
        return Meeting::Apart;
    }
    let foot = plane.project(sphere.centre());
    if (reach - sphere.radius()).abs() <= tol.confusion() {
        return Meeting::Touching(vec![foot]);
    }
    // The chord half-length: the leg of a right triangle whose hypotenuse is
    // the radius and whose other leg is the distance from the centre.
    let radius = sphere
        .radius()
        .mul_add(sphere.radius(), -(gap * gap))
        .max(0.0)
        .sqrt();
    match circle_on(foot, plane.normal(), radius, tol) {
        Some(circle) => Meeting::Along(vec![circle]),
        None => Meeting::Touching(vec![foot]),
    }
}

/// A plane and a cylinder.
///
/// Three genuinely different answers depending on the angle between them, and
/// the whole reason a closed form is worth having: a circle, an ellipse, or a
/// pair of straight lines, each exact.
fn plane_cylinder(
    plane: og_math::Plane,
    cylinder: og_math::Cylinder,
    tol: Tolerances,
) -> OgResult<Meeting> {
    let axis = cylinder.axis();
    let along = plane.normal().dot(axis.direction);

    // The plane contains the axis direction: the section is straight lines,
    // one for each side the plane cuts, or none if it misses.
    if along.abs() <= tol.angular() {
        let gap = plane.signed_distance_to(axis.location);
        let reach = gap.abs();
        if reach > cylinder.radius() + tol.confusion() {
            return Ok(Meeting::Apart);
        }
        // How far along the plane, from the foot of the axis, each line sits.
        let offset = cylinder
            .radius()
            .mul_add(cylinder.radius(), -(gap * gap))
            .max(0.0)
            .sqrt();
        let foot = plane.project(axis.location);
        let sideways =
            Direction::from_cross(plane.normal().vector(), axis.direction.vector(), tol)?;
        if offset <= tol.confusion() {
            // Tangent along one line.
            return Ok(Meeting::Along(vec![line_through(foot, axis.direction)]));
        }
        return Ok(Meeting::Along(vec![
            line_through(foot + sideways.vector() * offset, axis.direction),
            line_through(foot - sideways.vector() * offset, axis.direction),
        ]));
    }

    // Perpendicular to the axis: a circle of the cylinder's own radius.
    let centre = intersect_axis_plane(axis, plane, tol)?;
    if (along.abs() - 1.0).abs() <= tol.angular() {
        return Ok(
            match circle_on(centre, plane.normal(), cylinder.radius(), tol) {
                Some(circle) => Meeting::Along(vec![circle]),
                None => Meeting::Apart,
            },
        );
    }

    // Oblique: an ellipse. Its minor axis is the cylinder's radius, across the
    // slope; its major is that divided by the cosine of the tilt, along it.
    let minor = cylinder.radius();
    let major = minor / along.abs();
    // The minor axis runs where the plane and a plane perpendicular to the axis
    // agree — the cross of the two normals.
    let minor_direction =
        Direction::from_cross(plane.normal().vector(), axis.direction.vector(), tol)?;
    let major_direction =
        Direction::from_cross(minor_direction.vector(), plane.normal().vector(), tol)?;
    let frame = Frame::from_axes(
        centre,
        major_direction,
        minor_direction,
        plane.normal(),
        tol,
    )?;
    Ok(Meeting::Along(vec![
        og_geom::EllipseCurve::new(Ellipse::new(frame, major, minor, tol)?).into(),
    ]))
}

/// Two spheres: apart, tangent at a point, the same, or a circle.
fn sphere_sphere(a: og_math::Sphere, b: og_math::Sphere, tol: Tolerances) -> Meeting {
    let between = b.centre() - a.centre();
    let distance = between.magnitude();
    if distance <= tol.confusion() {
        return if (a.radius() - b.radius()).abs() <= tol.confusion() {
            Meeting::Same
        } else {
            // Concentric and different: one inside the other, never meeting.
            Meeting::Apart
        };
    }
    let (ra, rb) = (a.radius(), b.radius());
    if distance > ra + rb + tol.confusion() || distance < (ra - rb).abs() - tol.confusion() {
        return Meeting::Apart;
    }
    let Ok(direction) = Direction::new(between, tol) else {
        return Meeting::Apart;
    };
    // Where the plane of the intersection circle crosses the line of centres.
    let reach = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let centre = a.centre() + direction.vector() * reach;
    let squared = ra.mul_add(ra, -(reach * reach));
    if squared <= tol.confusion() * tol.confusion() {
        return Meeting::Touching(vec![centre]);
    }
    match circle_on(centre, direction, squared.max(0.0).sqrt(), tol) {
        Some(circle) => Meeting::Along(vec![circle]),
        None => Meeting::Touching(vec![centre]),
    }
}

/// Two cylinders sharing an axis.
///
/// The only cylinder pair with a closed form worth writing down. Two general
/// cylinders meet in a quartic space curve, which is what the marching
/// intersector is for.
fn coaxial_cylinders(
    a: og_math::Cylinder,
    b: og_math::Cylinder,
    tol: Tolerances,
) -> OgResult<Meeting> {
    if !a.axis().is_coaxial(b.axis(), tol) {
        og_bail!(
            NotDone,
            "two cylinders that do not share an axis meet in a quartic space \
             curve, which needs the general marching intersector"
        );
    }
    Ok(if (a.radius() - b.radius()).abs() <= tol.confusion() {
        Meeting::Same
    } else {
        // Same axis, different radii: one inside the other, touching nowhere.
        Meeting::Apart
    })
}

/// A cylinder and a sphere whose centre is on the cylinder's axis.
fn coaxial_cylinder_sphere(
    cylinder: og_math::Cylinder,
    sphere: og_math::Sphere,
    tol: Tolerances,
) -> OgResult<Meeting> {
    let axis = cylinder.axis();
    if axis.distance_to(sphere.centre()) > tol.confusion() {
        og_bail!(
            NotDone,
            "a sphere off a cylinder's axis meets it in a quartic space curve, \
             which needs the general marching intersector"
        );
    }
    let (r, radius) = (cylinder.radius(), sphere.radius());
    if r > radius + tol.confusion() {
        return Ok(Meeting::Apart);
    }
    if (r - radius).abs() <= tol.confusion() {
        // The sphere's equator lies on the cylinder, and they are tangent
        // along it rather than crossing.
        let centre = sphere.centre();
        return Ok(match circle_on(centre, axis.direction, r, tol) {
            Some(circle) => Meeting::Along(vec![circle]),
            None => Meeting::Apart,
        });
    }
    // Two circles, symmetric about the sphere's centre.
    let reach = radius.mul_add(radius, -(r * r)).max(0.0).sqrt();
    let mut out = Vec::with_capacity(2);
    for side in [reach, -reach] {
        let centre = sphere.centre() + axis.direction.vector() * side;
        if let Some(circle) = circle_on(centre, axis.direction, r, tol) {
            out.push(circle);
        }
    }
    Ok(if out.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(out)
    })
}

/// Where an axis crosses a plane.
fn intersect_axis_plane(
    axis: og_math::Axis,
    plane: og_math::Plane,
    tol: Tolerances,
) -> OgResult<Point> {
    let along = plane.normal().dot(axis.direction);
    if along.abs() <= tol.angular() {
        og_bail!(Domain, "the axis runs along the plane and never crosses it");
    }
    let t = -plane.signed_distance_to(axis.location) / along;
    Ok(axis.location + axis.direction.vector() * t)
}

/// A full circle in the plane through `centre` with the given normal.
fn circle_on(centre: Point, normal: Direction, radius: f64, tol: Tolerances) -> Option<Curve> {
    if radius <= tol.confusion() {
        return None;
    }
    // Any perpendicular will do for where the parameterization starts.
    let reference = if normal.vector().cross(Vector::X).magnitude() > 0.5 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = Direction::from_cross(normal.vector(), reference, tol).ok()?;
    let frame = Frame::new(centre, normal, x, tol).ok()?;
    Some(og_geom::CircleCurve::new(Circle::new(frame, radius, tol).ok()?).into())
}

/// An unbounded line through a point.
fn line_through(through: Point, direction: Direction) -> Curve {
    og_geom::LineCurve::new(og_math::Axis::new(through, direction)).into()
}
