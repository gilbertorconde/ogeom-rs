//! Surface/surface intersection: the analytic cases.
//!
//! Where two surfaces meet has a closed form for a specific and well-known set
//! of pairs, and a general answer that needs a marching intersector with a
//! fitting stage after it. This module is the first of those. It is deliberately
//! *not* a partial implementation of the second: a pair it cannot solve exactly
//! is reported as needing the general path, never approximated.
//!
//! # The exact cases
//!
//! They are common: plane against plane, plane against cylinder, sphere
//! against sphere. A mechanical part is mostly these, and running a marching
//! intersector over a pair whose answer is a circle is slower and less accurate
//! than writing down the circle.
//!
//! They are fast. No stepping, no refinement, no approximation stage.
//!
//! And they are exact. Every result here can be checked without reference to
//! anything but the two surfaces themselves: sample the curve, ask each
//! surface how far away it is, and the answer should be zero. That check is
//! the accuracy instrument in `tests/support/`.
//!
//! # What it reports
//!
//! Not just curves. Two surfaces can miss, touch at a point, meet along curves,
//! or be the same surface, and those are four different answers that downstream
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
/// form, which is a statement about the pair, not a failure to compute. The
/// general marching intersector is what answers those.
pub fn surface_surface(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    use SurfaceGeometry as S;
    match (a, b) {
        (S::Plane(p), S::Plane(q)) => Ok(plane_plane(
            p.plane(),
            q.plane(),
            window_reach(p).min(window_reach(q)),
            tol,
        )),
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
        (S::Sphere(s), S::Torus(t)) => axial_sphere_torus(s.sphere(), t.torus(), tol),
        (S::Torus(t), S::Sphere(s)) => axial_sphere_torus(s.sphere(), t.torus(), tol),
        (S::Plane(p), S::Cone(c)) => plane_cone(p.plane(), c.cone(), heights(c), tol),
        (S::Cone(c), S::Plane(p)) => plane_cone(p.plane(), c.cone(), heights(c), tol),
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

/// How far across a plane surface's window reaches: its diagonal, or a
/// billion units for one stated unbounded, which leaves the angle alone to
/// decide.
fn window_reach(plane: &ogeom_geom::PlaneSurface) -> f64 {
    let ((u0, u1), (v0, v1)) = ogeom_geom::Surface::domain(plane);
    let reach = (u1 - u0).hypot(v1 - v0);
    if reach.is_finite() && reach > 0.0 {
        reach.min(1e9)
    } else {
        1e9
    }
}

/// Whether a plane whose normal meets an axis at cosine `along` stands
/// square to it, for a section of a surface of revolution about that axis.
///
/// Looser than an angle test on purpose. Tilted by an angle t, the section
/// taken as square is off by the radius times t squared over two: a
/// micro-radian tilt costs a millionth of a micron on a unit radius, while
/// a plane a composed placement left a rounding error off square would
/// otherwise have no closed form at all. Two planes are the opposite case,
/// where the error grows with the extent, and [`plane_plane`] tests the
/// angle itself.
fn square_to_axis(along: f64, tol: Tolerances) -> bool {
    (along.abs() - 1.0).abs() <= tol.angular()
}

/// Two planes: apart, the same, or a line.
///
/// Parallel is decided across `reach`, the size of the region the two
/// stand over: planes whose normals differ by an angle t part by t times
/// the distance, so across the region they are one plane (or two parallel
/// ones) when that stays within the confusion distance. Planes built by
/// different routes to be coplanar differ by rounding, a hundred-billionth
/// of a radian, and across a part that is nothing; a microradian across a
/// hundred millimetres is a thousand times the confusion, and a line.
fn plane_plane(a: ogeom_math::Plane, b: ogeom_math::Plane, reach: f64, tol: Tolerances) -> Meeting {
    let turn = a.normal().angle(b.normal());
    let turn = turn.min(core::f64::consts::PI - turn);
    if turn <= tol.angular().max(tol.confusion() / reach) {
        // Parallel. Either the same plane or two that never meet, decided by
        // whether one contains the other's origin.
        return if a.distance_to(b.origin()) <= tol.confusion() {
            Meeting::Same
        } else {
            Meeting::Apart
        };
    }
    // The line of intersection runs along both normals' cross product, and
    // passes through the point nearest the first plane's origin that
    // satisfies both. Measured from that origin rather than from the world's:
    // the planes' own offsets from the world origin cancel to a few digits
    // when divided by the square of a small angle between them (two faces a
    // few hundredths of a milliradian apart would meet a third of a micron
    // off their true line), where the second plane's distance from a point
    // on the first is known to rounding.
    let Ok(direction) = Direction::from_cross(a.normal().vector(), b.normal().vector(), tol) else {
        return Meeting::Apart;
    };
    let (na, nb) = (a.normal().vector(), b.normal().vector());
    let dot = na.dot(nb);
    // sin² of the angle between the normals, from the cross product: one
    // minus the dot squared loses every digit of a small angle.
    let denominator = na.cross(nb).square_magnitude();
    if denominator <= 0.0 {
        return Meeting::Apart;
    }
    let origin = a.origin();
    let off = b.signed_distance_to(origin);
    let through = origin + (na * (off * dot) - nb * off) / denominator;
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
    if square_to_axis(along, tol) {
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
    // agree: the cross of the two normals.
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
        // Equal radii with intersecting axes: the one crossing whose quartic
        // factors, into the two ellipses in the axes' bisector planes, each
        // an oblique plane section the plane machinery already speaks. The
        // ellipses cross at the two points where the cylinders are tangent;
        // that is the crossing's geometry, stated exactly rather than
        // marched through.
        if (a.radius() - b.radius()).abs() <= tol.confusion() {
            let (da, db) = (a.axis().direction.vector(), b.axis().direction.vector());
            let normal = da.cross(db);
            if normal.magnitude() > tol.angular() {
                let (pa, pb) = (a.axis().location, b.axis().location);
                // Closest points of the two axis lines; coincident when the
                // axes genuinely intersect.
                let w = pb - pa;
                let dd = da.dot(db);
                let denom = dd.mul_add(-dd, 1.0);
                let s = dd.mul_add(-db.dot(w), da.dot(w)) / denom;
                let t = dd.mul_add(da.dot(w), -db.dot(w)) / denom;
                let on_a = pa + da * s;
                let on_b = pb + db * t;
                if on_a.distance(on_b) <= tol.confusion() {
                    let centre = on_a;
                    let mut curves = Vec::new();
                    for m in [da - db, da + db] {
                        if m.magnitude() <= tol.angular() {
                            continue;
                        }
                        let plane =
                            ogeom_math::Plane::through(centre, ogeom_math::Direction::new(m, tol)?);
                        if let Meeting::Along(mut found) = plane_cylinder(plane, a, tol)? {
                            curves.append(&mut found);
                        }
                    }
                    if !curves.is_empty() {
                        return Ok(Meeting::Along(curves));
                    }
                }
            }
        }
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
/// parallels. A plane through the axis: two meridians.
///
/// Those are the plane/torus configurations with a closed form worth the
/// name: an oblique plane, or one parallel to the axis and off it, meets a
/// torus in a quartic (with Villarceau's circles at exactly one magic
/// tilt), and that is the marching intersector's business. A plane through
/// the axis often holds the torus's seam, and a fitted section there never
/// meets the seam's own vertices. The blend machinery lives on this case:
/// a rolling ball's toroidal envelope is tangent to the plane it rolls on
/// along a circle, and that tangency must be *reported as the circle it is*,
/// the way a tangent plane reports its line on a cylinder; a tangential
/// answer with no curve in it would send the boolean above into a refusal.
fn axial_plane_torus(
    plane: ogeom_math::Plane,
    torus: ogeom_math::Torus,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = torus.axis();
    let along = plane.normal().dot(axis.direction);
    if along.abs() <= tol.angular()
        && plane.signed_distance_to(axis.location).abs() <= tol.confusion()
    {
        return Ok(meridians(plane, torus, tol));
    }
    if !square_to_axis(along, tol) {
        ogeom_bail!(
            NotDone,
            "a plane oblique to a torus's axis, or parallel to it and off it, \
             meets it in a quartic, which needs the general marching \
             intersector"
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
    // Two parallels, one either side of the tube. A spindle's tube swallows
    // the axis, and its inner parallel then lies on the tube's other half,
    // past the axis: the same circle, at the radius's size.
    let spread = minor.mul_add(minor, -(height * height)).max(0.0).sqrt();
    let circles: Vec<Curve> = [
        torus.major_radius() + spread,
        (torus.major_radius() - spread).abs(),
    ]
    .into_iter()
    .filter_map(|radius| circle_on(centre, axis.direction, radius, tol))
    .collect();
    Ok(if circles.is_empty() {
        Meeting::Apart
    } else {
        Meeting::Along(circles)
    })
}

/// A plane through a torus's axis: the two tube circles either side of the
/// axis, each starting on the outer equator as the torus's own meridians
/// do, so a section lying on the torus's seam starts where the seam does.
fn meridians(plane: ogeom_math::Plane, torus: ogeom_math::Torus, tol: Tolerances) -> Meeting {
    let axis = torus.axis();
    let normal = plane.normal();
    let Ok(out) = Direction::from_cross(axis.direction.vector(), normal.vector(), tol) else {
        return Meeting::Apart;
    };
    let circles: Vec<Curve> = [out.vector(), -out.vector()]
        .into_iter()
        .filter_map(|radial| {
            let centre = axis.location + radial * torus.major_radius();
            let x = Direction::new(radial, tol).ok()?;
            let frame = Frame::new(centre, normal, x, tol).ok()?;
            let circle = Circle::new(frame, torus.minor_radius(), tol).ok()?;
            Some(ogeom_geom::CircleCurve::new(circle).into())
        })
        .collect();
    Meeting::Along(circles)
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
    let minor = torus.minor_radius();
    // In a meridian plane the tube is two circles, centred the major radius
    // either side of the axis; the far one reaches this side only on a
    // spindle, whose tube swallows the axis. The cylinder's line meets each
    // it reaches: tangent at the tube's equator, or at mirrored heights.
    let mut heights: Vec<f64> = Vec::new();
    for reach in [
        (cylinder.radius() - torus.major_radius()).abs(),
        cylinder.radius() + torus.major_radius(),
    ] {
        if reach > minor + tol.confusion() {
            continue;
        }
        if (reach - minor).abs() <= tol.confusion() {
            heights.push(0.0);
            continue;
        }
        let rise = minor.mul_add(minor, -(reach * reach)).max(0.0).sqrt();
        heights.extend([rise, -rise]);
    }
    let mut kept: Vec<f64> = Vec::new();
    for height in heights {
        if kept.iter().all(|k| (k - height).abs() > tol.confusion()) {
            kept.push(height);
        }
    }
    let circles: Vec<Curve> = kept
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
/// The perpendicular slice is the configuration the rebuilds lean on (a
/// drafted wall's cap, a chamfer cone against the face it melts into), and
/// the answer is a circle framed on the cone's own frame, so a caller
/// re-deriving an edge finds its parameters where the old ones were. An
/// oblique plane meets a cone in a conic, which is the marching
/// intersector's business.
fn plane_cone(
    plane: ogeom_math::Plane,
    cone: ogeom_math::Cone,
    heights: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = cone.axis();
    let along = plane.normal().dot(axis.direction);
    if !square_to_axis(along, tol) && plane.signed_distance_to(cone.apex()).abs() <= tol.confusion()
    {
        return Ok(rulings_in(plane, &cone, heights, tol));
    }
    if !square_to_axis(along, tol) {
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
    // Past the apex, on the far nappe, the radius runs negative: the
    // parallel is the circle of its size, which the chart places half a
    // turn round.
    let centre = axis.location + axis.direction.vector() * height;
    Ok(match cone_parallel(&cone, centre, radius.abs(), tol) {
        Some(circle) => Meeting::Along(vec![circle]),
        None => Meeting::Apart,
    })
}

/// A cone surface's window along its axis.
fn heights(cone: &ogeom_geom::ConeSurface) -> (f64, f64) {
    ogeom_geom::Surface::domain(cone).1
}

/// A plane through a cone's apex: the rulings it holds. A plane steeper
/// than the cone holds two, one tangent to it holds the one it touches
/// along, and one shallower meets the cone at the apex alone. Each ruling
/// is stated over the cone's window of `heights`, one segment for each
/// nappe the window reaches.
fn rulings_in(
    plane: ogeom_math::Plane,
    cone: &ogeom_math::Cone,
    heights: (f64, f64),
    tol: Tolerances,
) -> Meeting {
    let apex = cone.apex();
    let apex_height = cone.frame().to_local(apex).z;
    // A ruling from the apex along `direction`, cut to the window: the
    // stretch on each side of the apex the window holds.
    let stated = |direction: Direction| -> Vec<Curve> {
        let climb = direction.vector().dot(cone.axis().direction.vector());
        if climb.abs() <= tol.angular() {
            return Vec::new();
        }
        let (lo, hi) = (heights.0.min(heights.1), heights.0.max(heights.1));
        let mut out = Vec::new();
        for (from, to) in [(lo.max(apex_height), hi), (lo, hi.min(apex_height))] {
            if !(from.is_finite() && to.is_finite()) || to - from <= tol.confusion() {
                continue;
            }
            let (a, b) = ((from - apex_height) / climb, (to - apex_height) / climb);
            if let Ok(line) = ogeom_geom::LineCurve::over(
                ogeom_math::Axis::new(apex, direction),
                a.min(b),
                a.max(b),
            ) {
                out.push(line.into());
            }
        }
        out
    };
    let a = cone.axis().direction.vector();
    let n = plane.normal().vector();
    let half = cone.half_angle().abs();
    let across = n.dot(a);
    // How far the plane leans off tangency: tangent planes make the cone's
    // half angle with the axis's perpendicular.
    let lean = across.abs().clamp(0.0, 1.0).asin() - half;
    if lean.abs() <= 1e-10 {
        // The ruling the plane touches along: the plane's own direction
        // nearest the axis.
        let toward = a - n * across;
        return match Direction::new(toward, tol).map(stated) {
            Ok(lines) if !lines.is_empty() => Meeting::Along(lines),
            _ => Meeting::Touching(vec![apex]),
        };
    }
    if lean > 0.0 {
        return Meeting::Touching(vec![apex]);
    }
    // Two rulings: a direction `cos a + sin (cos t e + sin t f)` lies in the
    // plane where `cos t` takes the value below.
    let side = n - a * across;
    let Ok(e) = Direction::new(side, tol) else {
        return Meeting::Touching(vec![apex]);
    };
    let e = e.vector();
    let f = a.cross(e);
    let (sin, cos) = half.sin_cos();
    let k = (-across * cos / (side.magnitude() * sin)).clamp(-1.0, 1.0);
    let s = (1.0 - k * k).max(0.0).sqrt();
    let lines: Vec<Curve> = [s, -s]
        .into_iter()
        .filter_map(|t| {
            let d = a * cos + (e * k + f * t) * sin;
            Direction::new(d, tol).ok()
        })
        .flat_map(stated)
        .collect();
    if lines.is_empty() {
        Meeting::Touching(vec![apex])
    } else {
        Meeting::Along(lines)
    }
}

/// A cylinder sharing a cone's axis: the parallels where the slant crosses
/// the cylinder's radius, one on each nappe.
///
/// The radius function is linear in height, so it reaches the cylinder's
/// radius once on the chart's own nappe and once on the far one, past the
/// apex, where it runs negative; each crossing is a parallel, and the chart
/// places the far one half a turn round. A cone face stopping short of its
/// apex keeps only the near one: the far one's image lies outside its
/// window.
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
/// parallel where the slants cross.
///
/// In height and radius coordinates along the shared axis each cone is a line,
/// and the crossing is one linear equation; the parallel there is a circle
/// unless it lands on the apex, which is a touch.
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
    // share a sense (`is_coaxial` checked), so the slopes compare directly.
    let lift = (b.axis().location - a.axis().location).dot(axis.direction.vector());
    let (slope_a, slope_b) = (a.half_angle().tan(), b.half_angle().tan());
    let (ref_a, ref_b) = (
        a.reference_radius(),
        slope_b.mul_add(-lift, b.reference_radius()),
    );
    if (slope_a - slope_b).abs() <= tol.angular() && (ref_a - ref_b).abs() <= tol.confusion() {
        return Ok(Meeting::Same);
    }
    // Where the radius lines cross, on the same nappes or (one radius the
    // other's negative) on opposite ones, past an apex: each crossing is a
    // parallel, or a shared apex where the radius there is zero.
    let mut heights: Vec<f64> = Vec::new();
    if (slope_a - slope_b).abs() > tol.angular() {
        heights.push((ref_b - ref_a) / (slope_a - slope_b));
    }
    if (slope_a + slope_b).abs() > tol.angular() {
        heights.push(-(ref_a + ref_b) / (slope_a + slope_b));
    }
    let mut circles: Vec<Curve> = Vec::new();
    let mut touches: Vec<Point> = Vec::new();
    for height in heights {
        let radius = a.radius_at(height).abs();
        let centre = axis.location + axis.direction.vector() * height;
        if radius <= tol.confusion() {
            if touches.iter().all(|t| t.distance(centre) > tol.confusion()) {
                touches.push(centre);
            }
            continue;
        }
        let repeated = circles.iter().any(|c| {
            matches!(c, Curve::Circle(k) if k.circle().centre().distance(centre) <= tol.confusion())
        });
        if !repeated && let Some(circle) = cone_parallel(&a, centre, radius, tol) {
            circles.push(circle);
        }
    }
    Ok(if !circles.is_empty() {
        Meeting::Along(circles)
    } else if !touches.is_empty() {
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
    // Profile circles in a meridian plane: each tube is two circles, centred
    // its major radius either side of the axis at its height, radii the
    // minors. A spindle's far circle reaches past the axis onto this side,
    // so every pairing is met.
    let (ra, rb) = (a.minor_radius(), b.minor_radius());
    let profile_points = profile_meetings(
        &[
            (ogeom_math::Point2::new(a.major_radius(), 0.0), ra),
            (ogeom_math::Point2::new(-a.major_radius(), 0.0), ra),
        ],
        &[
            (ogeom_math::Point2::new(b.major_radius(), lift), rb),
            (ogeom_math::Point2::new(-b.major_radius(), lift), rb),
        ],
        tol,
    );
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

/// Where circles in a meridian plane meet, on the plane's own side of the
/// axis (a positive distance from it): each meeting revolves into a
/// parallel. Every circle of `a` is met with every circle of `b`; the
/// meetings come back once each.
fn profile_meetings(
    a: &[(ogeom_math::Point2, f64)],
    b: &[(ogeom_math::Point2, f64)],
    tol: Tolerances,
) -> Vec<ogeom_math::Point2> {
    let mut points: Vec<ogeom_math::Point2> = Vec::new();
    for &(ca, ra) in a {
        for &(cb, rb) in b {
            let between = cb - ca;
            let distance = between.magnitude();
            // Concentric circles of different radii never meet; the same
            // circle is the callers' `Same`.
            if distance <= tol.confusion()
                || distance > ra + rb + tol.confusion()
                || distance < (ra - rb).abs() - tol.confusion()
            {
                continue;
            }
            let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
            let squared = ra.mul_add(ra, -(along * along));
            let direction = between * (1.0 / distance);
            let foot = ca + direction * along;
            let found = if squared <= tol.confusion() * tol.confusion() {
                vec![foot]
            } else {
                let offset =
                    ogeom_math::Vector2::new(-direction.y, direction.x) * squared.max(0.0).sqrt();
                vec![foot + offset, foot - offset]
            };
            for p in found {
                if p.x > tol.confusion() && points.iter().all(|q| q.distance(p) > tol.confusion()) {
                    points.push(p);
                }
            }
        }
    }
    points
}

/// A sphere centred on a torus's axis: in a meridian plane the sphere is a
/// circle on the axis and the tube two circles beside it, and each place
/// they meet revolves into a parallel; a sphere seated in the tube touches
/// it along one. A sphere off the axis meets the torus in a curve only the
/// marcher traces.
fn axial_sphere_torus(
    sphere: ogeom_math::Sphere,
    torus: ogeom_math::Torus,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let axis = torus.axis();
    let lift = (sphere.centre() - axis.location).dot(axis.direction.vector());
    let foot = axis.location + axis.direction.vector() * lift;
    if foot.distance(sphere.centre()) > tol.confusion() {
        ogeom_bail!(
            NotDone,
            "a sphere off a torus's axis meets it in a curve only the general \
             marching intersector can trace"
        );
    }
    let tube = torus.minor_radius();
    let points = profile_meetings(
        &[(ogeom_math::Point2::new(0.0, lift), sphere.radius())],
        &[
            (ogeom_math::Point2::new(torus.major_radius(), 0.0), tube),
            (ogeom_math::Point2::new(-torus.major_radius(), 0.0), tube),
        ],
        tol,
    );
    let circles: Vec<Curve> = points
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
