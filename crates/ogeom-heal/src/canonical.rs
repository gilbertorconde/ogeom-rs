//! Canonical simplification: recognizing that exact geometry is secretly
//! analytic — a B-spline surface that *is* a cylinder becomes the cylinder.
//!
//! The input is exact geometry, not samples of unknown provenance: the
//! decision is made against the surface's own equation, sampled on its own
//! chart, with its own normals. That is what separates this from reverse
//! engineering — nothing is guessed about what the data means, only checked
//! against a candidate the estimators propose. A candidate is accepted when
//! **every** sample sits within the caller's stated tolerance of it, and the
//! certificate is the worst deviation actually measured; a surface that is
//! genuinely free-form at that tolerance stays what it is.
//!
//! Why it matters: exchange formats routinely spell a cylinder out pointwise,
//! and every algorithm downstream — intersection, blending, the boolean's
//! closed forms — is faster and exacter on the analytic carrier. The
//! reference keeps this in its healing layer for the same reason.
//!
//! The estimators are classical, and shared with the sample-based recognizer
//! that lives outside the kernel: a plane is the mean normal; a sphere's
//! centre is where the normal lines meet, in least squares; a cylinder's
//! axis is the direction the normals avoid — their covariance's smallest
//! eigenvector — and a cone adds the linear taper of radius against height.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{
    CircleCurve, ConeSurface, Curve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface,
    SurfaceGeometry,
};
use ogeom_math::Circle;
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

/// What one face's simplification found.
#[derive(Debug, Clone, PartialEq)]
pub enum Simplified {
    /// The surface is a plane, to the deviation measured.
    Plane {
        /// The worst deviation any sample showed.
        worst: f64,
    },
    /// A cylinder.
    Cylinder {
        /// The recognized radius.
        radius: f64,
        /// The worst deviation any sample showed.
        worst: f64,
    },
    /// A cone.
    Cone {
        /// The recognized half angle, in radians.
        half_angle: f64,
        /// The worst deviation any sample showed.
        worst: f64,
    },
    /// A sphere.
    Sphere {
        /// The recognized radius.
        radius: f64,
        /// The worst deviation any sample showed.
        worst: f64,
    },
}

/// A report of what was simplified, face by face.
#[derive(Debug, Default)]
pub struct CanonicalReport {
    /// One entry per face whose surface became analytic.
    pub simplified: Vec<Simplified>,
    /// Faces examined and left as they were — already analytic, or
    /// genuinely free-form at the stated tolerance.
    pub untouched: usize,
}

/// Replace every free-form surface in `shape` that is secretly analytic —
/// within `tolerance` — by the plane, cylinder, cone or sphere it is.
///
/// Faces are rebuilt on the recognized carrier with exact pcurves; their
/// edges' space curves are untouched, because the curves were never wrong.
/// A surface that does not verify stays exactly as it was: the decision is
/// the product, and a wrong yes is a solid that measures nearly right with
/// the wrong surface under every later operation.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `tolerance` is not a positive distance, or a recognized face cannot be
/// rebuilt on its own boundary.
pub fn canonical_simplify(
    model: &mut Model,
    shape: &Shape,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<(Built, CanonicalReport)> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    let mut report = CanonicalReport::default();
    let mut history = ogeom_algo::History::new();
    let mut replaced: Vec<(Shape, Shape)> = Vec::new();
    // Edges whose free-form curves turned out to be lines or circles,
    // rebuilt once and shared: two faces meeting along one edge must keep
    // meeting along one edge.
    let mut edge_map: std::collections::HashMap<ogeom_topo::TShapeId, Shape> =
        std::collections::HashMap::new();

    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let Some(data) = model.node(&face).and_then(|n| match n.data() {
            NodeData::Face(d) => Some(d.clone()),
            _ => None,
        }) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            continue;
        };
        // Only free-form carriers are candidates; the analytic ones already
        // are what they are.
        if !matches!(surface, SurfaceGeometry::BSpline(_)) {
            report.untouched += 1;
            continue;
        }
        let placement = face.transform(model.datums())?;
        let world = {
            use ogeom_geom::Transformable as _;
            surface.transformed(&placement, tol)?
        };
        let Some(found) = recognize_surface(&world, tolerance, tol)? else {
            report.untouched += 1;
            continue;
        };
        let carrier: SurfaceGeometry = match &found {
            Simplified::Plane { .. } => {
                let plane = fit_plane(&world, tol)?;
                PlaneSurface::over(plane, (-1e5, 1e5), (-1e5, 1e5))?.into()
            }
            Simplified::Cylinder { .. } => {
                let cylinder = fit_cylinder(&world, tol)?;
                CylinderSurface::new(cylinder, (-1e5, 1e5))?.into()
            }
            Simplified::Cone { .. } => {
                let cone = fit_cone(&world, tol)?;
                ConeSurface::new(cone, (-1e5, 1e5))?.into()
            }
            Simplified::Sphere { .. } => {
                let sphere = fit_sphere(&world, tol)?;
                SphereSurface::new(sphere).into()
            }
        };
        // Rebuild the face on the carrier from its own boundary. The edges'
        // space curves never lied about *where* they run — but a rim spelt
        // as a B-spline that is exactly a circle becomes the circle first,
        // or the exact-pcurve machinery has nothing to project in closed
        // form. Curve recognition is held to the same certificate as the
        // surfaces': every sample within tolerance, or no change.
        let mut wires = Vec::new();
        for wire in model.ordered_children_of(&face)? {
            let mut edges = Vec::new();
            for edge in model.ordered_children_of(&wire)? {
                let mapped = simplified_edge(model, &edge, tolerance, &mut edge_map, tol)?;
                edges.push(if edge.orientation() == ogeom_topo::Orientation::Reversed {
                    mapped.reversed()
                } else {
                    mapped
                });
            }
            wires.push(edges);
        }
        let rebuilt = ogeom_algo::make_face_with_pcurves(model, carrier, &wires, tol)?.shape;
        let rebuilt = if face.orientation() == ogeom_topo::Orientation::Reversed {
            rebuilt.reversed()
        } else {
            rebuilt
        };
        history.modify(&face, rebuilt.clone());
        replaced.push((face, rebuilt));
        report.simplified.push(found);
    }

    if replaced.is_empty() {
        history.generate(shape, shape.clone());
        return Ok((Built::new(shape.clone(), history), report));
    }

    // Reassemble: untouched faces as they are, recognized ones rebuilt, the
    // shell sewn back on the shared edges.
    let mut faces = Vec::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        match replaced.iter().find(|(old, _)| old.node() == face.node()) {
            Some((_, new)) => faces.push(new.clone()),
            None => faces.push(face),
        }
    }
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    let out = match sewn.shells.as_slice() {
        [shell] if ogeom_algo::is_shell_closed(model, shell)? => {
            ogeom_algo::make_solid(model, core::slice::from_ref(shell))?.shape
        }
        [shell] => shell.clone(),
        _ => ogeom_algo::make_compound(model, &sewn.shells)?.shape,
    };
    history.modify(shape, out.clone());
    Ok((Built::new(out, history), report))
}

/// Sampled positions and unit normals over the surface's own chart.
fn samples(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<(Vec<Point>, Vec<Vector>)> {
    const N: usize = 9;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut points = Vec::with_capacity(N * N);
    let mut normals = Vec::with_capacity(N * N);
    for i in 0..N {
        for j in 0..N {
            // Interior offsets: chart corners of a swept patch can be
            // degenerate, and a degenerate sample's normal is noise.
            #[allow(clippy::cast_precision_loss, reason = "a sample index")]
            let fu = (i as f64 + 0.5) / N as f64;
            #[allow(clippy::cast_precision_loss, reason = "a sample index")]
            let fv = (j as f64 + 0.5) / N as f64;
            let u = u0 + (u1 - u0) * fu;
            let v = v0 + (v1 - v0) * fv;
            points.push(surface.point_at(u, v, tol)?);
            normals.push(surface.normal_at(u, v, tol)?.vector());
        }
    }
    Ok((points, normals))
}

/// The candidate the estimators propose, verified against every sample.
///
/// # Errors
///
/// As the samplers'.
pub fn recognize_surface(
    surface: &SurfaceGeometry,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Simplified>> {
    let (points, normals) = samples(surface, tol)?;
    // Most specific last: a plane verifies on nothing curved, a sphere and
    // cylinder disagree everywhere but a torus, so the order only decides
    // ties on genuinely ambiguous windows, where any verified answer is
    // right by the stated tolerance.
    if let Ok(plane) = fit_plane(surface, tol) {
        let worst = worst_of(&points, |p| plane.signed_distance_to(p).abs());
        if worst <= tolerance {
            return Ok(Some(Simplified::Plane { worst }));
        }
    }
    if let Ok(sphere) = fit_sphere(surface, tol) {
        let worst = worst_of(&points, |p| {
            (p.distance(sphere.centre()) - sphere.radius()).abs()
        });
        if worst <= tolerance {
            return Ok(Some(Simplified::Sphere {
                radius: sphere.radius(),
                worst,
            }));
        }
    }
    if let Ok(cylinder) = fit_cylinder(surface, tol) {
        let worst = worst_of(&points, |p| cylinder.distance_to(p).abs());
        if worst <= tolerance {
            return Ok(Some(Simplified::Cylinder {
                radius: cylinder.radius(),
                worst,
            }));
        }
    }
    if let Ok(cone) = fit_cone(surface, tol) {
        let worst = worst_of(&points, |p| cone.distance_to(p).abs());
        if worst <= tolerance {
            return Ok(Some(Simplified::Cone {
                half_angle: cone.half_angle(),
                worst,
            }));
        }
    }
    let _ = normals;
    Ok(None)
}

fn worst_of(points: &[Point], f: impl Fn(Point) -> f64) -> f64 {
    points.iter().map(|p| f(*p)).fold(0.0, f64::max)
}

/// The edge, its free-form curve replaced by the line or circle it verifiably
/// is — or the edge itself when it is analytic already or genuinely free.
fn simplified_edge(
    model: &mut Model,
    edge: &Shape,
    tolerance: f64,
    cache: &mut std::collections::HashMap<ogeom_topo::TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if let Some(found) = cache.get(&edge.node()) {
        return Ok(found.clone());
    }
    let Some((curve_id, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        return Ok(edge.clone());
    };
    let Some(curve) = model.geometry().curve(curve_id).cloned() else {
        return Ok(edge.clone());
    };
    if !matches!(curve, Curve::BSpline(_)) {
        return Ok(edge.clone());
    }
    let placement = edge.transform(model.datums())?;
    let world = {
        use ogeom_geom::Transformable as _;
        curve.transformed(&placement, tol)?
    };
    const N: usize = 17;
    let mut pts = Vec::with_capacity(N);
    for i in 0..N {
        #[allow(clippy::cast_precision_loss, reason = "a sample index")]
        let t = range.0 + (range.1 - range.0) * i as f64 / (N - 1) as f64;
        pts.push(world.point_at(t, tol)?);
    }
    let bounds = model.children_of(edge)?;
    let (Some(va), Some(vb)) = (bounds.first().cloned(), bounds.last().cloned()) else {
        return Ok(edge.clone());
    };
    let start = pts[0];
    let end = pts[N - 1];
    let closed = start.distance(end) <= tol.confusion();

    // A line: every sample on the chord.
    if !closed && start.distance(end) > tol.confusion() {
        let dir = Direction::new(end - start, tol)?;
        let worst = worst_of(&pts, |p| (p - start).cross(dir.vector()).magnitude());
        if worst <= tolerance {
            let line = LineCurve::segment(start, end, tol)?;
            let built = ogeom_algo::make_edge_between(
                model,
                Curve::from(line),
                (0.0, start.distance(end)),
                &va,
                &vb,
                tol,
            )?
            .shape;
            cache.insert(edge.node(), built.clone());
            return Ok(built);
        }
    }

    // A circle: planar, and equidistant from the fitted centre.
    if let Ok(plane) = fit_plane_points(&pts, tol) {
        let planar = worst_of(&pts, |p| plane.signed_distance_to(p).abs());
        if planar <= tolerance {
            let frame = plane.frame();
            let flat: Vec<(f64, f64)> = pts
                .iter()
                .map(|p| {
                    let l = frame.to_local(*p);
                    (l.x, l.y)
                })
                .collect();
            if let Ok((cx, cy, r)) = fit_circle_2d(&flat) {
                let centre = frame.origin() + frame.x().vector() * cx + frame.y().vector() * cy;
                let round = worst_of(&pts, |p| (p.distance(centre) - r).abs());
                if round <= tolerance {
                    // The frame's x-axis runs through the start point, so a
                    // closed rim's range is (0, τ) inside the curve's own
                    // domain rather than a window shifted past it.
                    let x = Direction::new(start - centre, tol)?;
                    let cframe = Frame::new(centre, frame.z(), x, tol)?;
                    let circle = Circle::new(cframe, r, tol)?;
                    let angle = |p: Point| -> f64 {
                        let l = cframe.to_local(p);
                        l.y.atan2(l.x).rem_euclid(core::f64::consts::TAU)
                    };
                    let (t0, t1) = if closed {
                        // The frame's x-axis runs through the start by
                        // construction, so its angle is zero — computed, it
                        // is −ε, which rem_euclid folds to τ and the range
                        // overflows the domain by a full turn.
                        (0.0, core::f64::consts::TAU)
                    } else {
                        let a = angle(start);
                        let mut b = angle(end);
                        // The arc runs the way the samples do.
                        let m = angle(pts[N / 2]);
                        let fwd = (m - a).rem_euclid(core::f64::consts::TAU)
                            <= (b - a).rem_euclid(core::f64::consts::TAU);
                        if !fwd {
                            b -= core::f64::consts::TAU;
                        }
                        if b <= a {
                            b += core::f64::consts::TAU;
                        }
                        (a, b)
                    };
                    let built = ogeom_algo::make_edge_between(
                        model,
                        Curve::from(CircleCurve::new(circle)),
                        (t0, t1),
                        &va,
                        &vb,
                        tol,
                    )?
                    .shape;
                    cache.insert(edge.node(), built.clone());
                    return Ok(built);
                }
            }
        }
    }
    Ok(edge.clone())
}

/// The best plane through raw points: centroid plus the cross-product mean.
fn fit_plane_points(points: &[Point], tol: Tolerances) -> OgeomResult<Plane> {
    let c = centroid(points);
    let mut n = Vector::ZERO;
    for w in points.windows(2) {
        n += (w[0] - c).cross(w[1] - c);
    }
    Ok(Plane::through(c, Direction::new(n, tol)?))
}

/// The plane through the samples: mean point, mean normal.
fn fit_plane(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Plane> {
    let (points, normals) = samples(surface, tol)?;
    let n = mean(&normals);
    let centroid = centroid(&points);
    Ok(Plane::through(centroid, Direction::new(n, tol)?))
}

/// The sphere whose centre is where the normal lines meet, in least squares:
/// minimise Σ‖(I − nnᵀ)(c − p)‖², a 3×3 solve.
fn fit_sphere(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Sphere> {
    let (points, normals) = samples(surface, tol)?;
    let mut a = nalgebra::Matrix3::<f64>::zeros();
    let mut b = nalgebra::Vector3::<f64>::zeros();
    for (p, n) in points.iter().zip(&normals) {
        let nv = nalgebra::Vector3::new(n.x, n.y, n.z);
        let proj = nalgebra::Matrix3::identity() - nv * nv.transpose();
        a += proj;
        b += proj * nalgebra::Vector3::new(p.x, p.y, p.z);
    }
    let c = a
        .lu()
        .solve(&b)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "normal lines meet nowhere"))?;
    let centre = Point::new(c.x, c.y, c.z);
    #[allow(clippy::cast_precision_loss, reason = "a sample count")]
    let count = points.len().max(1) as f64;
    let radius = points.iter().map(|p| p.distance(centre)).sum::<f64>() / count;
    Sphere::centred(centre, radius, tol)
}

/// The direction the normals avoid: the smallest eigenvector of Σ nnᵀ. A
/// cylinder's normals are all perpendicular to its axis; a cone's make one
/// constant angle with it, which leaves the same eigenvector.
fn normal_axis(normals: &[Vector], tol: Tolerances) -> OgeomResult<Direction> {
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    for n in normals {
        let v = nalgebra::Vector3::new(n.x, n.y, n.z);
        m += v * v.transpose();
    }
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut best = 0;
    for i in 1..3 {
        if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
            best = i;
        }
    }
    let d = eigen.eigenvectors.column(best);
    Direction::from_coords(d.x, d.y, d.z, tol)
}

fn fit_cylinder(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Cylinder> {
    let (points, normals) = samples(surface, tol)?;
    let axis = normal_axis(&normals, tol)?;
    // Project onto the plane through the centroid perpendicular to the axis
    // and fit the circle there: its centre lifts to a point on the axis.
    let origin = centroid(&points);
    let frame = frame_about(origin, axis, tol)?;
    let flat: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let l = frame.to_local(*p);
            (l.x, l.y)
        })
        .collect();
    let (cx, cy, r) = fit_circle_2d(&flat)?;
    let centre = frame.origin() + frame.x().vector() * cx + frame.y().vector() * cy;
    Cylinder::new(frame_about(centre, axis, tol)?, r, tol)
}

fn fit_cone(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Cone> {
    let (points, normals) = samples(surface, tol)?;
    let axis = normal_axis(&normals, tol)?;
    // Radius against height is affine on a cone; the slope is the taper.
    let origin = centroid(&points);
    let frame = frame_about(origin, axis, tol)?;
    let hs: Vec<f64> = points.iter().map(|p| frame.to_local(*p).z).collect();
    // The section centre may be off the centroid; fit the circle at two
    // height bands to place the axis exactly.
    let (lo, hi) = hs
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), h| {
            (a.min(*h), b.max(*h))
        });
    if hi - lo <= tol.confusion() {
        ogeom_bail!(Construction, "a flat band tapers to nothing");
    }
    let band = |keep: &dyn Fn(f64) -> bool| -> Vec<(f64, f64)> {
        points
            .iter()
            .filter(|p| keep(frame.to_local(**p).z))
            .map(|p| {
                let l = frame.to_local(*p);
                (l.x, l.y)
            })
            .collect()
    };
    let mid = f64::midpoint(lo, hi);
    let lower = band(&|h| h <= mid);
    let upper = band(&|h| h > mid);
    let (ax, ay, r_lo) = fit_circle_2d(&lower)?;
    let (bx, by, r_hi) = fit_circle_2d(&upper)?;
    let h_lo = lower_mean(&points, &frame, mid, true);
    let h_hi = lower_mean(&points, &frame, mid, false);
    if (h_hi - h_lo).abs() <= tol.confusion() {
        ogeom_bail!(Construction, "the bands coincide");
    }
    let slope = (r_hi - r_lo) / (h_hi - h_lo);
    let half_angle = slope.atan();
    // Reference circle at h = 0, centre interpolated along the fitted axis.
    let t = -h_lo / (h_hi - h_lo);
    let cx = ax + (bx - ax) * t;
    let cy = ay + (by - ay) * t;
    let r0 = r_lo + (r_hi - r_lo) * t;
    let centre = frame.origin() + frame.x().vector() * cx + frame.y().vector() * cy;
    Cone::new(frame_about(centre, axis, tol)?, r0, half_angle, tol)
}

fn lower_mean(points: &[Point], frame: &Frame, mid: f64, low: bool) -> f64 {
    let hs: Vec<f64> = points
        .iter()
        .map(|p| frame.to_local(*p).z)
        .filter(|h| (*h <= mid) == low)
        .collect();
    #[allow(clippy::cast_precision_loss, reason = "a sample count")]
    let n = hs.len().max(1) as f64;
    hs.iter().sum::<f64>() / n
}

/// Least-squares circle through 2D points, by the |c|² − r² substitution.
fn fit_circle_2d(points: &[(f64, f64)]) -> OgeomResult<(f64, f64, f64)> {
    if points.len() < 3 {
        ogeom_bail!(Construction, "a circle needs three points");
    }
    // x² + y² = 2 cx x + 2 cy y + (r² − cx² − cy²): linear in (cx, cy, k).
    let (mut sxx, mut sxy, mut sx, mut syy, mut sy, mut s1) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut bx, mut by, mut b1) = (0.0, 0.0, 0.0);
    for (x, y) in points {
        let rhs = x * x + y * y;
        sxx += 2.0 * x * 2.0 * x;
        sxy += 2.0 * x * 2.0 * y;
        sx += 2.0 * x;
        syy += 2.0 * y * 2.0 * y;
        sy += 2.0 * y;
        s1 += 1.0;
        bx += 2.0 * x * rhs;
        by += 2.0 * y * rhs;
        b1 += rhs;
    }
    let a = nalgebra::Matrix3::new(sxx, sxy, sx, sxy, syy, sy, sx, sy, s1);
    let sol = a
        .lu()
        .solve(&nalgebra::Vector3::new(bx, by, b1))
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "the points close no circle"))?;
    let (cx, cy, k) = (sol.x, sol.y, sol.z);
    let r2 = k + cx * cx + cy * cy;
    if r2 <= 0.0 {
        ogeom_bail!(Construction, "the points close no circle");
    }
    Ok((cx, cy, r2.sqrt()))
}

fn centroid(points: &[Point]) -> Point {
    let mut sum = Vector::ZERO;
    for p in points {
        sum += p.to_vector();
    }
    #[allow(clippy::cast_precision_loss, reason = "a sample count")]
    let n = points.len().max(1) as f64;
    Point::ORIGIN + sum * (1.0 / n)
}

fn mean(vs: &[Vector]) -> Vector {
    let mut sum = Vector::ZERO;
    for v in vs {
        sum += *v;
    }
    #[allow(clippy::cast_precision_loss, reason = "a sample count")]
    let n = vs.len().max(1) as f64;
    sum * (1.0 / n)
}

fn frame_about(origin: Point, axis: Direction, tol: Tolerances) -> OgeomResult<Frame> {
    let seed = if axis.vector().dot(Vector::X).abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = Direction::from_cross(axis.vector(), seed, tol)?;
    Frame::new(origin, axis, x, tol)
}
