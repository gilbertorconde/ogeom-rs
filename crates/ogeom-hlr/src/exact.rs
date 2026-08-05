//! The exact half of the drawing pipeline.
//!
//! The polygonal path draws what the *mesh* says: silhouettes are interior
//! mesh edges whose triangles disagree about facing the eye, and visibility
//! is occlusion sampling against the triangles. Both are as good as the
//! chord, and no better.
//!
//! This is the other half. A silhouette is where the surface's own normal
//! turns perpendicular to the view, which for the elementary surfaces is a
//! curve with a closed form: a great circle on a sphere, a pair of rulings
//! on a cylinder or a cone. Visibility is decided by asking the *faces*
//! whether anything stands between a point and the eye — an exact
//! curve/surface interference and a trim test, not a triangle count. What is
//! still sampled is the drawing itself, because a drawing is polylines; the
//! curves it samples and the classification it carries are exact.
//!
//! A surface whose silhouette has no closed form — a torus, a spline — is
//! refused by name rather than approximated here. The polygonal path draws
//! those, and says so.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Axis, Circle, Direction, Frame, Point, Point2, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

use crate::project::{Drawing, DrawnCurve, Source, View, Visibility};

/// One silhouette curve, and the face it belongs to.
#[derive(Debug, Clone)]
pub struct Silhouette {
    /// The face whose surface turns away here.
    pub face: Shape,
    /// The curve, in space.
    pub curve: Curve,
    /// The portion of it that lies within the face's trim.
    pub range: (f64, f64),
}

/// The exact silhouettes of a shape, seen along `direction`.
///
/// A silhouette is the locus where the surface normal is perpendicular to
/// the view: for a sphere the great circle whose plane the direction is
/// normal to, for a cylinder the two rulings furthest to either side, for a
/// cone the two rulings through its apex where the same holds. Planes have
/// none — a plane either faces the eye or does not — and a face whose
/// surface has no closed-form silhouette is refused by name.
///
/// Each curve comes back trimmed to the stretch that lies within its own
/// face, decided by sampling the face's trim, so a silhouette on a face
/// that was cut away is absent rather than drawn through thin air.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// direction has no length, or a face carries a surface whose silhouette has
/// no closed form.
pub fn silhouettes(
    model: &Model,
    shape: &Shape,
    direction: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Silhouette>> {
    let magnitude = direction.magnitude();
    if magnitude <= tol.confusion() {
        ogeom_bail!(Construction, "a silhouette needs a direction to look along");
    }
    let along = direction / magnitude;
    let deflection = Deflection::default();

    let mut out = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        let world = ogeom_geom::Transformable::transformed(&surface, &placement, tol)?;
        let candidates = match &world {
            SurfaceGeometry::Plane(_) => Vec::new(),
            SurfaceGeometry::Sphere(s) => {
                // The great circle whose plane has the view as its normal:
                // every normal on it is radial, so every one is
                // perpendicular to the view.
                let sphere = s.sphere();
                let axis = Direction::new(along, tol)?;
                let frame = Frame::new(sphere.centre(), axis, perpendicular(along, tol)?, tol)?;
                vec![Curve::Circle(CircleCurve::new(Circle::new(
                    frame,
                    sphere.radius(),
                    tol,
                )?))]
            }
            SurfaceGeometry::Cylinder(c) => {
                // The two rulings where the radial direction is
                // perpendicular to the view: the axis stepped sideways by
                // the radius, either way.
                let cylinder = c.cylinder();
                let axis = cylinder.frame().z().vector();
                let sideways = axis.cross(along);
                let m = sideways.magnitude();
                if m <= tol.angular() {
                    // Looking down the axis: the whole rim is the outline,
                    // and the face's own boundary already draws it.
                    Vec::new()
                } else {
                    let sideways = sideways / m;
                    [1.0, -1.0]
                        .iter()
                        .map(|sign| {
                            let at =
                                cylinder.frame().origin() + sideways * (cylinder.radius() * sign);
                            Curve::Line(LineCurve::new(Axis::new(at, cylinder.frame().z())))
                        })
                        .collect()
                }
            }
            SurfaceGeometry::Cone(c) => {
                // The same question on a cone: the rulings whose own normal
                // is perpendicular to the view. The normal of a ruling at
                // angle u is radial tilted by the half-angle, so the
                // condition is a linear one in (cos u, sin u) and has two
                // roots — or none, when the eye is inside the cone's own
                // angle and nothing turns away.
                let cone = c.cone();
                let frame = cone.frame();
                let (x, y, z) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
                let (sin, cos) = cone.half_angle().sin_cos();
                // n(u) = cos(half) * (x cos u + y sin u) - sin(half) * z
                let (a, b) = (cos * along.dot(x), cos * along.dot(y));
                let c0 = -sin * along.dot(z);
                let r = a.hypot(b);
                if r <= tol.angular() || c0.abs() > r {
                    Vec::new()
                } else {
                    let phase = b.atan2(a);
                    let spread = (-c0 / r).acos();
                    [phase + spread, phase - spread]
                        .iter()
                        .map(|u| {
                            let radial = x * u.cos() + y * u.sin();
                            let apex = cone.apex();
                            let direction =
                                radial * cone.half_angle().sin() + z * cone.half_angle().cos();
                            Direction::new(direction, tol)
                                .map(|d| Curve::Line(LineCurve::new(Axis::new(apex, d))))
                        })
                        .collect::<OgeomResult<Vec<Curve>>>()?
                }
            }
            _ => ogeom_bail!(
                Construction,
                "the silhouette of this surface has no closed form; the \
                 polygonal drawing path draws it, and says so"
            ),
        };

        for curve in candidates {
            for range in within_trim(model, &face, &world, &curve, deflection, tol)? {
                out.push(Silhouette {
                    face: face.clone(),
                    curve: curve.clone(),
                    range,
                });
            }
        }
    }
    Ok(out)
}

/// The reflect lines of a shape under a light: where the surface turns away
/// from the *light* rather than from the eye.
///
/// The same locus as a silhouette, asked of a different direction — which is
/// what a reflect line is, and why the two share a construction. A surface
/// inspected this way shows its own creases: the lines move a long way for a
/// small change in curvature.
///
/// # Errors
///
/// As [`silhouettes`].
pub fn reflect_lines(
    model: &Model,
    shape: &Shape,
    light: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Silhouette>> {
    silhouettes(model, shape, light, tol)
}

/// The isoparametric curves of a face: `u_count` at constant `u`, `v_count`
/// at constant `v`, each trimmed to the stretches that lie on the face.
///
/// Evenly spaced across the face's own parameter window, excluding its
/// edges, because an isoparametric at the window's edge is the face's
/// boundary and the boundary is already drawn.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// face carries no surface.
pub fn iso_curves(
    model: &Model,
    face: &Shape,
    u_count: usize,
    v_count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Point>>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;
    let world = ogeom_geom::Transformable::transformed(&surface, &placement, tol)?;
    let ((u0, u1), (v0, v1)) = world.domain();
    let rings = ogeom_mesh::face_boundary(model, face, Deflection::default(), tol)?;

    const ALONG: usize = 64;
    let mut out = Vec::new();
    for (count, constant_u) in [(u_count, true), (v_count, false)] {
        for i in 1..=count {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a curve index, far below the mantissa"
            )]
            let f = i as f64 / (count + 1) as f64;
            let mut run: Vec<Point> = Vec::new();
            for k in 0..=ALONG {
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "a station index, far below the mantissa"
                )]
                let g = k as f64 / ALONG as f64;
                let (u, v) = if constant_u {
                    ((u1 - u0).mul_add(f, u0), (v1 - v0).mul_add(g, v0))
                } else {
                    ((u1 - u0).mul_add(g, u0), (v1 - v0).mul_add(f, v0))
                };
                if inside_rings(&rings, Point2::new(u, v)) {
                    run.push(world.point_at(u, v, tol)?);
                } else if run.len() >= 2 {
                    out.push(std::mem::take(&mut run));
                } else {
                    run.clear();
                }
            }
            if run.len() >= 2 {
                out.push(run);
            }
        }
    }
    Ok(out)
}

/// Project a shape into a drawing whose silhouettes and visibility are
/// exact.
///
/// The edges and silhouettes are sampled at `deflection` to become
/// polylines, because a drawing is polylines. What is not sampled is the
/// *geometry* they sample — exact curves on the surfaces, not mesh edges —
/// or the classification, which asks the faces themselves whether anything
/// stands between a point and the eye.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// face's silhouette has no closed form, or the shape has no faces.
pub fn project_exact(
    model: &Model,
    shape: &Shape,
    view: &View,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Drawing> {
    let faces = blockers(model, shape, tol)?;
    if faces.is_empty() {
        ogeom_bail!(Construction, "a shape with no faces draws nothing");
    }
    let mut drawing = Drawing::default();

    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Ok(points) = ogeom_mesh::polyline_of_edge(model, &edge, deflection, tol) else {
            continue;
        };
        classify(
            &mut drawing,
            &points,
            Source::Edge(edge.clone()),
            view,
            &faces,
            tol,
        )?;
    }

    for silhouette in silhouettes(model, shape, view.toward_eye(), tol)? {
        let points = sampled(&silhouette.curve, silhouette.range, deflection, tol)?;
        classify(&mut drawing, &points, Source::Silhouette, view, &faces, tol)?;
    }
    Ok(drawing)
}

/// A face that can stand between a point and the eye: its surface in world
/// space and its trim as chart rings.
struct Blocker {
    surface: SurfaceGeometry,
    rings: Vec<Vec<Point2>>,
}

fn blockers(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Blocker>> {
    let mut out = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        out.push(Blocker {
            surface: ogeom_geom::Transformable::transformed(&surface, &placement, tol)?,
            rings: ogeom_mesh::face_boundary(model, &face, Deflection::default(), tol)?,
        });
    }
    Ok(out)
}

/// Split a polyline into visible and hidden runs, asking the faces.
fn classify(
    drawing: &mut Drawing,
    points: &[Point],
    source: Source,
    view: &View,
    faces: &[Blocker],
    tol: Tolerances,
) -> OgeomResult<()> {
    let mut run: Vec<Point2> = Vec::new();
    let mut held: Option<Visibility> = None;
    let mut flush = |run: &mut Vec<Point2>, visibility: Option<Visibility>| {
        if run.len() < 2 {
            run.clear();
            return;
        }
        let curve = DrawnCurve {
            points: std::mem::take(run),
            visibility: visibility.unwrap_or(Visibility::Visible),
            source: source.clone(),
        };
        if visibility == Some(Visibility::Hidden) {
            drawing.hidden.push(curve);
        } else {
            drawing.visible.push(curve);
        }
    };
    for point in points {
        let visibility = if occluded(*point, view, faces, tol)? {
            Visibility::Hidden
        } else {
            Visibility::Visible
        };
        if held.is_some_and(|was| was != visibility) {
            // The change happens somewhere between the two samples; the run
            // ends at the sample that changed, so the two runs meet there.
            let last = run.last().copied();
            flush(&mut run, held);
            if let Some(last) = last {
                run.push(last);
            }
        }
        held = Some(visibility);
        run.push(view.project(*point));
    }
    flush(&mut run, held);
    Ok(())
}

/// Whether anything stands between `at` and the eye.
fn occluded(at: Point, view: &View, faces: &[Blocker], tol: Tolerances) -> OgeomResult<bool> {
    let toward = view.toward_eye();
    let magnitude = toward.magnitude();
    if magnitude <= tol.confusion() {
        return Ok(false);
    }
    let direction = toward / magnitude;
    // Far enough to leave any shape it started inside, and started far
    // enough along not to strike the surface the point is on.
    let reach = 1e6;
    let clearance = tol.confusion() * 1e3;
    let ray = Curve::Line(LineCurve::new(Axis::new(
        at,
        Direction::new(direction, tol)?,
    )));
    let options = ogeom_intersect::CurveSurfaceOptions::default();
    for face in faces {
        let found = ogeom_intersect::intersect_curve_surface(&ray, &face.surface, options, tol)?;
        for piercing in &found.crossings {
            if piercing.on_curve <= clearance || piercing.on_curve >= reach {
                continue;
            }
            let (u, v) = piercing.on_surface;
            if inside_rings(&face.rings, Point2::new(u, v)) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// The stretches of a curve that lie within a face's trim.
fn within_trim(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    curve: &Curve,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64)>> {
    let rings = ogeom_mesh::face_boundary(model, face, deflection, tol)?;
    let (t0, t1) = curve.domain();
    // A line's domain is the whole real line as far as the type is
    // concerned; a silhouette on one is only interesting where the face is,
    // so an unbounded curve is walked over the face's own reach.
    let (t0, t1) = if t0.is_finite() && t1.is_finite() && t1 - t0 < 1e6 {
        (t0, t1)
    } else {
        let mut bound = ogeom_math::Aabb::EMPTY;
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                bound = bound.with_point(vertex.transform(model.datums())?.apply(data.point));
            }
        }
        let reach = bound.diagonal().max(1.0);
        let centre = bound.centre().unwrap_or(Point::ORIGIN);
        let at = ogeom_algo::project_on_curve(curve, centre, 64, tol)?.parameter;
        (at - reach, at + reach)
    };

    const STATIONS: usize = 96;
    let held_at = |t: f64| -> OgeomResult<bool> {
        let Ok(point) = curve.point_at(t, tol) else {
            return Ok(false);
        };
        let projection = ogeom_algo::project_on_surface(surface, point, 24, tol)?;
        let (u, v) = projection.parameters;
        // The projection clamps to the surface's own window, so a point
        // just past the end of a face comes back with a foot *at* the end
        // and the overshoot as its distance. Holding that to the confusion
        // tolerance is what stops a silhouette running off its own face.
        Ok(projection.distance <= tol.confusion() && inside_rings(&rings, Point2::new(u, v)))
    };
    // Where the answer changes between two stations, the edge of the face
    // is between them; bisecting says where to a part in a million of a
    // station, so a silhouette ends *on* its face rather than a station
    // past it.
    let edge_between = |inside: f64, outside: f64| -> OgeomResult<f64> {
        let (mut lo, mut hi) = (inside, outside);
        for _ in 0..40 {
            let mid = f64::midpoint(lo, hi);
            if held_at(mid)? {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    };

    let mut out = Vec::new();
    let mut open: Option<f64> = None;
    let mut previous: Option<(f64, bool)> = None;
    for k in 0..=STATIONS {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = (t1 - t0).mul_add(k as f64 / STATIONS as f64, t0);
        let held = held_at(t)?;
        match (held, open) {
            (true, None) => {
                open = Some(match previous {
                    Some((was, false)) => edge_between(t, was)?,
                    _ => t,
                });
            }
            (false, Some(from)) => {
                let to = match previous {
                    Some((was, true)) => edge_between(was, t)?,
                    _ => t,
                };
                if to - from > tol.parametric() {
                    out.push((from, to));
                }
                open = None;
            }
            _ => {}
        }
        previous = Some((t, held));
    }
    if let Some(from) = open
        && t1 - from > tol.parametric()
    {
        out.push((from, t1));
    }
    Ok(out)
}

/// A curve's polyline over a range, at the given deflection.
fn sampled(
    curve: &Curve,
    range: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    // How finely to walk: enough steps that the chord between two of them
    // stays inside the deflection, from the curve's own length.
    let span = curve
        .point_at(range.0, tol)?
        .distance(curve.point_at(range.1, tol)?);
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a step count, clamped into range"
    )]
    let steps =
        ((span / deflection.chord.max(tol.confusion())).sqrt().ceil() as usize).clamp(8, 512);
    let mut out = Vec::with_capacity(steps + 1);
    for k in 0..=steps {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = (range.1 - range.0).mul_add(k as f64 / steps as f64, range.0);
        out.push(curve.point_at(t, tol)?);
    }
    Ok(out)
}

/// Even-odd containment against chart rings.
fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = (b.x - a.x).mul_add((p.y - a.y) / (b.y - a.y), a.x);
                if x > p.x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// Any unit vector perpendicular to `v`.
fn perpendicular(v: Vector, tol: Tolerances) -> OgeomResult<Direction> {
    let seed = if v.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    Direction::new(v.cross(seed), tol)
}
