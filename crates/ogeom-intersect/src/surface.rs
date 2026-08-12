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
//! the instrument the intersection gate is measured with (`docs/PLAN.md`),
//! and it only exists because these cases are exact. A benchmark whose reference
//! answers came from the thing being benchmarked would measure nothing.
//!
//! # What it reports
//!
//! Not just curves. Two surfaces can miss, touch at a point, meet along curves,
//! or be the same surface — and those are four different answers that downstream
//! code has to distinguish. A boolean that treats coincidence as "no
//! intersection" produces a solid with a face missing.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Ellipse, Frame, Point, Vector};

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
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if this pair has no closed
/// form — which is a statement about the pair, not a failure to compute. The
/// general marching intersector is what answers those, and it is gated on the
/// benchmark this module makes possible.
pub fn surface_surface(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
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
        (S::Plane(p), S::Torus(t)) => axial_plane_torus(p.plane(), t.torus(), tol),
        (S::Torus(t), S::Plane(p)) => axial_plane_torus(p.plane(), t.torus(), tol),
        (S::Cylinder(c), S::Torus(t)) => coaxial_cylinder_torus(c.cylinder(), t.torus(), tol),
        (S::Torus(t), S::Cylinder(c)) => coaxial_cylinder_torus(c.cylinder(), t.torus(), tol),
        (S::Torus(x), S::Torus(y)) => coaxial_tori(x.torus(), y.torus(), tol),
        (S::Plane(p), S::Cone(c)) => plane_cone(p.plane(), c.cone(), tol),
        (S::Cone(c), S::Plane(p)) => plane_cone(p.plane(), c.cone(), tol),
        (S::Cylinder(x), S::Cone(c)) => coaxial_cylinder_cone(x.cylinder(), c.cone(), tol),
        (S::Cone(c), S::Cylinder(x)) => coaxial_cylinder_cone(x.cylinder(), c.cone(), tol),
        (S::Cone(x), S::Cone(y)) => coaxial_cones(x.cone(), y.cone(), tol),
        _ => ogeom_bail!(
            NotDone,
            "this pair of surfaces has no closed-form intersection; it needs \
             the general marching intersector, which is gated on the benchmark \
             these cases provide the ground truth for"
        ),
    }
}

/// Two planes: apart, the same, or a line.
fn plane_plane(a: ogeom_math::Plane, b: ogeom_math::Plane, tol: Tolerances) -> Meeting {
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
fn plane_sphere(plane: ogeom_math::Plane, sphere: ogeom_math::Sphere, tol: Tolerances) -> Meeting {
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
    plane: ogeom_math::Plane,
    cylinder: ogeom_math::Cylinder,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
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
        ogeom_geom::EllipseCurve::new(Ellipse::new(frame, major, minor, tol)?).into(),
    ]))
}

/// Two spheres: apart, tangent at a point, the same, or a circle.
fn sphere_sphere(a: ogeom_math::Sphere, b: ogeom_math::Sphere, tol: Tolerances) -> Meeting {
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
    a: ogeom_math::Cylinder,
    b: ogeom_math::Cylinder,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    if !a.axis().is_coaxial(b.axis(), tol) {
        ogeom_bail!(
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
    cylinder: ogeom_math::Cylinder,
    sphere: ogeom_math::Sphere,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = cylinder.axis();
    if axis.distance_to(sphere.centre()) > tol.confusion() {
        ogeom_bail!(
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

/// A plane perpendicular to a torus's axis: apart, one tangent circle, or two
/// parallels.
///
/// The perpendicular slice is the only plane/torus configuration with a
/// closed form worth the name — an oblique plane meets a torus in a quartic
/// (with Villarceau's circles at exactly one magic tilt), and that is the
/// marching intersector's business. The blend machinery lives on this case:
/// a rolling ball's toroidal envelope is tangent to the plane it rolls on
/// along a circle, and that tangency must be *reported as the circle it is*,
/// the way a tangent plane reports its line on a cylinder — a tangential
/// answer with no curve in it would send the boolean above into a refusal.
fn axial_plane_torus(
    plane: ogeom_math::Plane,
    torus: ogeom_math::Torus,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = torus.axis();
    let along = plane.normal().dot(axis.direction);
    if (along.abs() - 1.0).abs() > tol.angular() {
        ogeom_bail!(
            NotDone,
            "a plane oblique or parallel to a torus's axis meets it in a \
             quartic, which needs the general marching intersector"
        );
    }
    // The plane's height above the tube's centre plane.
    let height = -plane.signed_distance_to(axis.location) * along.signum();
    let minor = torus.minor_radius();
    if height.abs() > minor + tol.confusion() {
        return Ok(Meeting::Apart);
    }
    let centre = axis.location + axis.direction.vector() * height;
    if (height.abs() - minor).abs() <= tol.confusion() {
        // Tangent along the parallel at the tube's top or bottom.
        return Ok(
            match circle_on(centre, axis.direction, torus.major_radius(), tol) {
                Some(circle) => Meeting::Along(vec![circle]),
                None => Meeting::Apart,
            },
        );
    }
    // Two parallels, one either side of the tube — the inner one only where
    // the tube does not swallow the axis.
    let spread = minor.mul_add(minor, -(height * height)).max(0.0).sqrt();
    let circles: Vec<Curve> = [torus.major_radius() + spread, torus.major_radius() - spread]
        .into_iter()
        .filter_map(|radius| circle_on(centre, axis.direction, radius, tol))
        .collect();
    Ok(if circles.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(circles)
    })
}

/// A cylinder sharing a torus's axis: apart, one tangent circle, or two
/// parallels at mirrored heights.
fn coaxial_cylinder_torus(
    cylinder: ogeom_math::Cylinder,
    torus: ogeom_math::Torus,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    if !cylinder.axis().is_coaxial(torus.axis(), tol) {
        ogeom_bail!(
            NotDone,
            "a cylinder off a torus's axis meets it in a quartic space curve, \
             which needs the general marching intersector"
        );
    }
    let axis = torus.axis();
    let reach = (cylinder.radius() - torus.major_radius()).abs();
    let minor = torus.minor_radius();
    if reach > minor + tol.confusion() {
        return Ok(Meeting::Apart);
    }
    if (reach - minor).abs() <= tol.confusion() {
        // Tangent along the tube's inner or outer equator.
        return Ok(
            match circle_on(axis.location, axis.direction, cylinder.radius(), tol) {
                Some(circle) => Meeting::Along(vec![circle]),
                None => Meeting::Apart,
            },
        );
    }
    let rise = minor.mul_add(minor, -(reach * reach)).max(0.0).sqrt();
    let circles: Vec<Curve> = [rise, -rise]
        .into_iter()
        .filter_map(|height| {
            circle_on(
                axis.location + axis.direction.vector() * height,
                axis.direction,
                cylinder.radius(),
                tol,
            )
        })
        .collect();
    Ok(if circles.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(circles)
    })
}

/// A plane square to a cone's axis: the parallel at that height, or the apex.
///
/// The perpendicular slice is the configuration the rebuilds lean on — a
/// drafted wall's cap, a chamfer cone against the face it melts into — and
/// the answer is a circle framed on the cone's own frame, so a caller
/// re-deriving an edge finds its parameters where the old ones were. An
/// oblique plane meets a cone in a conic, which is the marching
/// intersector's business.
fn plane_cone(
    plane: ogeom_math::Plane,
    cone: ogeom_math::Cone,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = cone.axis();
    let along = plane.normal().dot(axis.direction);
    if (along.abs() - 1.0).abs() > tol.angular() {
        ogeom_bail!(
            NotDone,
            "a plane oblique to a cone's axis meets it in a conic, which \
             needs the general marching intersector"
        );
    }
    // The plane's height along the axis, from the cone frame's origin.
    let height = -plane.signed_distance_to(axis.location) * along.signum();
    let radius = cone.radius_at(height);
    if radius.abs() <= tol.confusion() {
        // The plane passes through the apex, where the parallel has no
        // length: a touch, not a curve.
        return Ok(Meeting::Touching(vec![cone.apex()]));
    }
    // A negative radius is the far nappe, where the surface's own
    // parameterization runs half a turn out of phase; the parallel is the
    // same circle either way.
    let centre = axis.location + axis.direction.vector() * height;
    Ok(match cone_parallel(&cone, centre, radius.abs(), tol) {
        Some(circle) => Meeting::Along(vec![circle]),
        None => Meeting::Apart,
    })
}

/// A cylinder sharing a cone's axis: one parallel per nappe.
///
/// The slant crosses any coaxial cylinder's radius exactly once per nappe —
/// the radius function is linear in height and spans every value — so the
/// answer is always two circles, one of them on the far nappe.
fn coaxial_cylinder_cone(
    cylinder: ogeom_math::Cylinder,
    cone: ogeom_math::Cone,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    if !cylinder.axis().is_coaxial(cone.axis(), tol) {
        ogeom_bail!(
            NotDone,
            "a cylinder off a cone's axis meets it in a curve only the \
             general marching intersector can trace"
        );
    }
    let axis = cone.axis();
    let slope = cone.half_angle().tan();
    let circles: Vec<Curve> = [cylinder.radius(), -cylinder.radius()]
        .into_iter()
        .filter_map(|radius| {
            let height = (radius - cone.reference_radius()) / slope;
            cone_parallel(
                &cone,
                axis.location + axis.direction.vector() * height,
                cylinder.radius(),
                tol,
            )
        })
        .collect();
    Ok(if circles.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(circles)
    })
}

/// Two cones sharing an axis: the same surface, the shared apex, or the
/// parallels where the slants cross.
///
/// In height–radius coordinates along the shared axis each cone is a line,
/// and each nappe pairing is one linear equation; a solution's parallel is a
/// circle unless it lands on the apex, which is a touch.
fn coaxial_cones(
    a: ogeom_math::Cone,
    b: ogeom_math::Cone,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    if !a.axis().is_coaxial(b.axis(), tol) {
        ogeom_bail!(
            NotDone,
            "two cones that do not share an axis meet in a curve only the \
             general marching intersector can trace"
        );
    }
    let axis = a.axis();
    // Both radius functions expressed against `a`'s height origin. The axes
    // share a sense — `is_coaxial` checked — so the slopes compare directly.
    let lift = (b.axis().location - a.axis().location).dot(axis.direction.vector());
    let (slope_a, slope_b) = (a.half_angle().tan(), b.half_angle().tan());
    let (ref_a, ref_b) = (
        a.reference_radius(),
        slope_b.mul_add(-lift, b.reference_radius()),
    );
    if (slope_a - slope_b).abs() <= tol.angular() {
        // Parallel slants: the same cone, or two that never meet.
        return Ok(if (ref_a - ref_b).abs() <= tol.confusion() {
            Meeting::Same
        } else {
            Meeting::Apart
        });
    }
    if (slope_a + slope_b).abs() <= tol.angular() && (ref_a + ref_b).abs() <= tol.confusion() {
        // Mirror cones: the same double cone traversed the other way. Point
        // sets coincide but the parameterizations run opposite nappes, and
        // unifying them is not this module's call to make.
        ogeom_bail!(
            NotDone,
            "two mirror cones share their double cone; resolving that \
             coincidence needs the general machinery"
        );
    }
    let mut touches = Vec::new();
    let mut circles = Vec::new();
    // One equation per nappe pairing: like signs and opposite signs.
    for (rhs_ref, rhs_slope) in [(ref_b, slope_b), (-ref_b, -slope_b)] {
        let run = slope_a - rhs_slope;
        if run.abs() <= tol.angular() {
            continue;
        }
        let height = (rhs_ref - ref_a) / run;
        let radius = a.radius_at(height);
        if radius.abs() <= tol.confusion() {
            touches.push(a.apex());
        } else if let Some(circle) = cone_parallel(
            &a,
            axis.location + axis.direction.vector() * height,
            radius.abs(),
            tol,
        ) {
            circles.push(circle);
        }
    }
    Ok(if !circles.is_empty() {
        Meeting::Along(circles)
    } else if !touches.is_empty() {
        touches.dedup_by(|p, q| (*p - *q).magnitude() <= tol.confusion());
        Meeting::Touching(touches)
    } else {
        Meeting::Apart
    })
}

/// A parallel of a cone, framed on the cone's own frame so parameters carry.
fn cone_parallel(
    cone: &ogeom_math::Cone,
    centre: Point,
    radius: f64,
    tol: Tolerances,
) -> Option<Curve> {
    if radius <= tol.confusion() {
        return None;
    }
    let frame = cone.frame();
    let placed = Frame::new(centre, frame.z(), frame.x(), tol).ok()?;
    Some(ogeom_geom::CircleCurve::new(Circle::new(placed, radius, tol).ok()?).into())
}

/// Two tori sharing an axis: the same surface, apart, or circles where the
/// tube profiles cross.
///
/// In the shared meridian half-plane the two tubes are two circles, and
/// revolving their meetings gives the answer: radical-line algebra in the
/// `(distance-from-axis, height)` plane, each solution a parallel.
fn coaxial_tori(
    a: ogeom_math::Torus,
    b: ogeom_math::Torus,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    if !a.axis().is_coaxial(b.axis(), tol) {
        ogeom_bail!(
            NotDone,
            "two tori that do not share an axis meet in a curve only the \
             general marching intersector can trace"
        );
    }
    let axis = a.axis();
    let lift = (b.axis().location - a.axis().location).dot(axis.direction.vector());
    if (a.major_radius() - b.major_radius()).abs() <= tol.confusion()
        && lift.abs() <= tol.confusion()
        && (a.minor_radius() - b.minor_radius()).abs() <= tol.confusion()
    {
        return Ok(Meeting::Same);
    }
    // Profile circles in the meridian half-plane: centres at
    // `(major, height)`, radii the minors.
    let (ca, cb) = (
        ogeom_math::Point2::new(a.major_radius(), 0.0),
        ogeom_math::Point2::new(b.major_radius(), lift),
    );
    let between = cb - ca;
    let distance = between.magnitude();
    let (ra, rb) = (a.minor_radius(), b.minor_radius());
    if distance <= tol.confusion() {
        // Concentric profiles of different tube radii never meet; the same
        // circle was the `Same` case above.
        return Ok(Meeting::Apart);
    }
    if distance > ra + rb + tol.confusion() || distance < (ra - rb).abs() - tol.confusion() {
        return Ok(Meeting::Apart);
    }
    let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let squared = ra.mul_add(ra, -(along * along));
    let direction = between * (1.0 / distance);
    let foot = ca + direction * along;
    let mut profile_points = Vec::new();
    if squared <= tol.confusion() * tol.confusion() {
        profile_points.push(foot);
    } else {
        let offset = ogeom_math::Vector2::new(-direction.y, direction.x) * squared.max(0.0).sqrt();
        profile_points.push(foot + offset);
        profile_points.push(foot - offset);
    }
    let circles: Vec<Curve> = profile_points
        .into_iter()
        .filter_map(|p| {
            circle_on(
                axis.location + axis.direction.vector() * p.y,
                axis.direction,
                p.x,
                tol,
            )
        })
        .collect();
    Ok(if circles.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(circles)
    })
}

/// Where an axis crosses a plane.
fn intersect_axis_plane(
    axis: ogeom_math::Axis,
    plane: ogeom_math::Plane,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let along = plane.normal().dot(axis.direction);
    if along.abs() <= tol.angular() {
        ogeom_bail!(Domain, "the axis runs along the plane and never crosses it");
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
    Some(ogeom_geom::CircleCurve::new(Circle::new(frame, radius, tol).ok()?).into())
}

/// An unbounded line through a point.
fn line_through(through: Point, direction: Direction) -> Curve {
    ogeom_geom::LineCurve::new(ogeom_math::Axis::new(through, direction)).into()
}
