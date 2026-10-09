//! Where a curve pierces a surface.
//!
//! *Elsewhere* this is `GeomAPI_IntCS` and the line/quadric half of `IntAna`.
//! Two consumers drive it: edge/face interference in the boolean's pave
//! filler, and the exact point-in-solid classifier, which is a ray/surface
//! query per face.
//!
//! # Well-posed
//!
//! `C(t) = S(u, v)` is three equations in three unknowns; unlike the
//! surface/surface system, nothing has to be pinned for Newton to converge to
//! a point. The analytic cases are still answered in closed form first: a line
//! against a plane or a quadric is a linear or quadratic equation, and solving
//! a quadratic by iteration would be slower and less exact than writing down
//! its roots.
//!
//! # A curve lying in the surface
//!
//! A line in a plane crosses it nowhere and everywhere. That is an overlap
//! (the parameter range of the curve that lies in the surface), and it is a
//! different answer from any list of points. Detected where the analytic
//! forms can see it, and on the general path as a run of piercings between
//! which the curve never leaves the surface.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d, Surface, SurfaceGeometry};
use ogeom_math::{Point, solve};

use crate::march::{Cell, sample_by, segment_meets_triangle};

/// One piercing of a surface by a curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Piercing {
    /// The parameter on the curve.
    pub on_curve: f64,
    /// The parameters on the surface.
    pub on_surface: (f64, f64),
    /// Where, taken from the curve.
    pub point: Point,
    /// The distance between the two evaluations there.
    pub gap: f64,
}

/// What a curve does to a surface.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSurfaceIntersection {
    /// Isolated piercings, in order along the curve.
    pub crossings: Vec<Piercing>,
    /// Parameter ranges of the curve that lie *in* the surface.
    ///
    /// Detected for the analytic cases (a line in a plane) and, on the
    /// general path, as a run of piercings longer than the curve's sampling
    /// step between which the curve never leaves the surface.
    pub lying: Vec<(f64, f64)>,
}

impl CurveSurfaceIntersection {
    /// No contact found.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.crossings.is_empty() && self.lying.is_empty()
    }

    const fn empty() -> Self {
        Self {
            crossings: Vec::new(),
            lying: Vec::new(),
        }
    }
}

/// How hard the general path looks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceOptions {
    /// How many segments the curve is sampled into for seeding.
    pub samples: usize,
    /// How finely the surface is sampled, per direction.
    pub grid: usize,
    /// The widest gap that still counts as a piercing.
    pub gap: f64,
}

impl Default for CurveSurfaceOptions {
    fn default() -> Self {
        Self {
            samples: 128,
            grid: 24,
            gap: 1e-7,
        }
    }
}

/// Where a curve pierces a surface.
///
/// Analytic line/plane, line/sphere and line/cylinder are answered in closed
/// form; everything else is seeded polyhedrally and polished by Newton on the
/// well-posed three-by-three system.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the options
/// are unusable.
pub fn intersect_curve_surface(
    curve: &Curve,
    surface: &SurfaceGeometry,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> OgeomResult<CurveSurfaceIntersection> {
    intersect_with(curve, surface, options, tol, || {
        std::borrow::Cow::Owned(sample_by(surface, seeding(surface, options.grid), tol))
    })
}

/// A surface made ready for many curve queries.
///
/// The general path seeds from a sampling of the surface into flat cells,
/// and that sampling depends only on the surface and the options. Taken
/// once, on the first query that needs it, and shared by every later one,
/// so a caller asking the same face about many curves (a ray per drawn
/// point, a ray per classified point) pays for it once. Every answer is the
/// one [`intersect_curve_surface`] gives for the same curve.
#[derive(Debug)]
pub struct PreparedSurface {
    surface: SurfaceGeometry,
    options: CurveSurfaceOptions,
    tol: Tolerances,
    cells: std::sync::OnceLock<Vec<Cell>>,
}

impl PreparedSurface {
    /// Ready `surface` for queries at `options` and `tol`.
    #[must_use]
    pub const fn new(
        surface: SurfaceGeometry,
        options: CurveSurfaceOptions,
        tol: Tolerances,
    ) -> Self {
        Self {
            surface,
            options,
            tol,
            cells: std::sync::OnceLock::new(),
        }
    }

    /// The surface queried.
    #[must_use]
    pub const fn surface(&self) -> &SurfaceGeometry {
        &self.surface
    }

    /// Where `curve` pierces the surface, as [`intersect_curve_surface`].
    ///
    /// # Errors
    ///
    /// As [`intersect_curve_surface`].
    pub fn intersect(&self, curve: &Curve) -> OgeomResult<CurveSurfaceIntersection> {
        let (surface, options, tol) = (&self.surface, self.options, self.tol);
        intersect_with(curve, surface, options, tol, || {
            std::borrow::Cow::Borrowed(
                self.cells
                    .get_or_init(|| sample_by(surface, seeding(surface, options.grid), tol))
                    .as_slice(),
            )
        })
    }
}

fn intersect_with<'c>(
    curve: &Curve,
    surface: &SurfaceGeometry,
    options: CurveSurfaceOptions,
    tol: Tolerances,
    cells: impl FnOnce() -> std::borrow::Cow<'c, [Cell]>,
) -> OgeomResult<CurveSurfaceIntersection> {
    if options.samples < 2 || options.grid < 2 {
        ogeom_bail!(Construction, "seeding needs at least two steps each way");
    }
    if !options.gap.is_finite() || options.gap <= 0.0 {
        ogeom_bail!(Construction, "a gap of {} is not a distance", options.gap);
    }

    match (curve, surface) {
        (Curve::Line(line), SurfaceGeometry::Plane(p)) => {
            Ok(line_plane(line, p.plane(), curve, surface, tol))
        }
        (Curve::Line(line), SurfaceGeometry::Sphere(s)) => Ok(line_quadric(
            line,
            curve,
            surface,
            sphere_roots(line, s.sphere()),
            options,
            tol,
        )),
        (Curve::Line(line), SurfaceGeometry::Cylinder(c)) => Ok(line_quadric(
            line,
            curve,
            surface,
            cylinder_roots(line, c.cylinder()),
            options,
            tol,
        )),
        (Curve::Line(line), SurfaceGeometry::Cone(c)) => Ok(line_quadric(
            line,
            curve,
            surface,
            cone_roots(line, c.cone(), tol),
            options,
            tol,
        )),
        (Curve::Line(line), SurfaceGeometry::Torus(t)) => Ok(line_quadric(
            line,
            curve,
            surface,
            torus_roots(line, t.torus()),
            options,
            tol,
        )),
        _ => general(curve, surface, &cells(), options, tol),
    }
}

// --- analytic ----------------------------------------------------------------

fn line_plane(
    line: &ogeom_geom::LineCurve,
    plane: ogeom_math::Plane,
    curve: &Curve,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> CurveSurfaceIntersection {
    let axis = line.axis();
    let along = plane.normal().dot(axis.direction);
    let height = plane.signed_distance_to(axis.location);

    if along.abs() <= tol.angular() {
        // Parallel: in the plane, or never touching it.
        if height.abs() <= tol.confusion() {
            return CurveSurfaceIntersection {
                crossings: Vec::new(),
                lying: vec![line.domain()],
            };
        }
        return CurveSurfaceIntersection::empty();
    }

    let t = -height / along;
    let (lo, hi) = line.domain();
    if t < lo - tol.parametric() || t > hi + tol.parametric() {
        return CurveSurfaceIntersection::empty();
    }
    let point = axis.location + axis.direction.vector() * t;
    let Some(found) = invert(surface, point, curve, t, tol) else {
        return CurveSurfaceIntersection::empty();
    };
    if found.gap > tol.confusion() {
        // The crossing is real on the unbounded plane but outside this
        // surface's stated extents; the clamped polish says so as a gap.
        return CurveSurfaceIntersection::empty();
    }
    CurveSurfaceIntersection {
        crossings: vec![found],
        lying: Vec::new(),
    }
}

/// How far rounding can move a discriminant `p - q` from its true value.
fn rounding(p: f64, q: f64) -> f64 {
    8.0 * f64::EPSILON * p.abs().max(q.abs())
}

/// The line parameters at which a line meets a sphere.
fn sphere_roots(line: &ogeom_geom::LineCurve, sphere: ogeom_math::Sphere) -> Vec<f64> {
    let axis = line.axis();
    let d = axis.direction.vector();
    let m = axis.location - sphere.centre();
    // |m + t d|^2 = r^2, with |d| = 1.
    let b = m.dot(d);
    let c = sphere.radius().mul_add(-sphere.radius(), m.dot(m));
    let discriminant = b.mul_add(b, -c);
    // A tangent line has a discriminant of exactly zero, which rounding
    // leaves a few ulps either side; below zero by no more than that is the
    // tangency, and the polish that follows rejects a false one by its gap.
    if discriminant < -rounding(b * b, c) {
        return Vec::new();
    }
    let root = discriminant.max(0.0).sqrt();
    if root == 0.0 {
        vec![-b]
    } else {
        vec![-b - root, -b + root]
    }
}

/// The line parameters at which a line meets a cylinder.
/// A line's offset from a frame's origin and its direction, in the frame's
/// own axes.
fn in_frame(
    line: &ogeom_geom::LineCurve,
    frame: ogeom_math::Frame,
) -> (ogeom_math::Vector, ogeom_math::Vector) {
    let axis = line.axis();
    let local = |v: ogeom_math::Vector| {
        ogeom_math::Vector::new(
            v.dot(frame.x().vector()),
            v.dot(frame.y().vector()),
            v.dot(frame.z().vector()),
        )
    };
    (
        local(axis.location - frame.origin()),
        local(axis.direction.vector()),
    )
}

/// The line parameters at which a line meets a cone, either nappe: the
/// polish that follows keeps those on the cone's stated extent.
fn cone_roots(
    line: &ogeom_geom::LineCurve,
    cone: ogeom_math::Cone,
    tol: ogeom_core::Tolerances,
) -> Vec<f64> {
    let (m, d) = in_frame(line, cone.frame());
    // x^2 + y^2 = (r0 + k z)^2 along m + t d.
    let (r0, k) = (cone.reference_radius(), cone.half_angle().tan());
    let rim = k.mul_add(m.z, r0);
    let a = d.x.mul_add(d.x, d.y * d.y) - k * k * d.z * d.z;
    let b = 2.0 * (m.x.mul_add(d.x, m.y * d.y) - k * d.z * rim);
    let c = m.x.mul_add(m.x, m.y * m.y) - rim * rim;
    ogeom_math::solve::roots(&[c, b, a], tol.parametric()).unwrap_or_default()
}

/// The line parameters at which a line meets a torus: the real roots of
/// `(|w|^2 + R^2 - r^2)^2 = 4 R^2 (w_x^2 + w_y^2)` along `w = m + t d` in
/// the torus's frame.
fn torus_roots(line: &ogeom_geom::LineCurve, torus: ogeom_math::Torus) -> Vec<f64> {
    let (m, d) = in_frame(line, torus.frame());
    let (big, small) = (torus.major_radius(), torus.minor_radius());
    // Solved about the line's nearest approach to the centre, in units of
    // the torus's own size, so the roots stand near one: a torus hundreds
    // of millimetres down the line puts them there otherwise, and the
    // coefficients' spread costs them digits.
    let near = -m.dot(d) / d.dot(d);
    let size = big + small;
    let m = (m + d * near) / size;
    let (big, small) = (big / size, small / size);
    let a = d.dot(d);
    let b = 2.0 * m.dot(d);
    let c = m.dot(m) + big * big - small * small;
    let p = d.x.mul_add(d.x, d.y * d.y);
    let q = 2.0 * m.x.mul_add(d.x, m.y * d.y);
    let s = m.x.mul_add(m.x, m.y * m.y);
    let four = 4.0 * big * big;
    let coefficients = [
        c.mul_add(c, -four * s),
        2.0f64.mul_add(b * c, -four * q),
        b.mul_add(b, 2.0 * a * c) - four * p,
        2.0 * a * b,
        a * a,
    ];
    let (found, count) = quartic_roots(&coefficients);
    found[..count]
        .iter()
        .map(|sigma| near + sigma * size)
        .collect()
}

/// The real roots of a quartic whose roots stand near one, tangencies
/// included: a critical point where the quartic is within 1e-9 of its coefficients'
/// size of zero is a tangency, the coefficients being differences of terms
/// that size. Two roots a few roots of epsilon apart are merged into the
/// one tangency rounding split them from.
fn quartic_roots(c: &[f64; 5]) -> ([f64; 4], usize) {
    let mut found = [0.0; 4];
    let touch = 1e-9 * c.iter().map(|x| x.abs()).sum::<f64>();
    let count = ogeom_math::solve::real_roots(c, touch, &mut found);
    let mut merged = 0;
    for i in 0..count {
        let t = found[i];
        if merged > 0 && (t - found[merged - 1]).abs() <= 1e-6 * (1.0 + t.abs()) {
            found[merged - 1] = f64::midpoint(found[merged - 1], t);
        } else {
            found[merged] = t;
            merged += 1;
        }
    }
    (found, merged)
}

fn cylinder_roots(line: &ogeom_geom::LineCurve, cylinder: ogeom_math::Cylinder) -> Vec<f64> {
    let axis = line.axis();
    let w = cylinder.axis().direction.vector();
    // Strip the components along the cylinder's axis; what is left is a 2D
    // circle problem in the perpendicular plane.
    let d = axis.direction.vector();
    let m = axis.location - cylinder.axis().location;
    let d_perp = d - w * d.dot(w);
    let m_perp = m - w * m.dot(w);
    let a = d_perp.dot(d_perp);
    if a <= f64::MIN_POSITIVE {
        // The line runs along the axis direction: on the wall it would lie,
        // not pierce, and lying is not detected here.
        return Vec::new();
    }
    let b = d_perp.dot(m_perp);
    let c = cylinder
        .radius()
        .mul_add(-cylinder.radius(), m_perp.dot(m_perp));
    let discriminant = b.mul_add(b, -(a * c));
    if discriminant < -rounding(b * b, a * c) {
        return Vec::new();
    }
    let root = discriminant.max(0.0).sqrt();
    if root == 0.0 {
        vec![-b / a]
    } else {
        vec![(-b - root) / a, (-b + root) / a]
    }
}

/// Roots dressed as piercings, filtered by the line's own range and the
/// surface's extents.
fn line_quadric(
    line: &ogeom_geom::LineCurve,
    curve: &Curve,
    surface: &SurfaceGeometry,
    roots: Vec<f64>,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> CurveSurfaceIntersection {
    let axis = line.axis();
    let (lo, hi) = line.domain();
    let mut crossings = Vec::new();
    for t in roots {
        if t < lo - tol.parametric() || t > hi + tol.parametric() {
            continue;
        }
        let point = axis.location + axis.direction.vector() * t;
        let Some(found) = invert(surface, point, curve, t, tol) else {
            continue;
        };
        // The extent check is the gap. The polish clamps the surface
        // parameters into the stated domain, so a root beyond the cylinder's
        // height converges to the rim with a gap of exactly how far past it
        // was: a piercing of the unbounded geometry, not of this surface.
        if found.gap > tol.confusion() {
            continue;
        }
        let _ = options;
        crossings.push(Piercing {
            on_curve: t,
            on_surface: found.on_surface,
            point,
            gap: found.gap,
        });
    }
    crossings.sort_by(|a, b| {
        a.on_curve
            .partial_cmp(&b.on_curve)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    CurveSurfaceIntersection {
        crossings,
        lying: Vec::new(),
    }
}

/// Surface parameters of a point known to lie on an analytic surface.
///
/// Closed-form inversion for the quadrics; refined by one Newton pass so the
/// reported parameters evaluate back onto the point to rounding.
fn invert(
    surface: &SurfaceGeometry,
    point: Point,
    curve: &Curve,
    on_curve: f64,
    tol: Tolerances,
) -> Option<Piercing> {
    let guess = match surface {
        SurfaceGeometry::Plane(p) => {
            let local = p.plane().frame().to_local(point);
            (local.x, local.y)
        }
        SurfaceGeometry::Sphere(s) => {
            let local = s.sphere().frame().to_local(point);
            let latitude = (local.z / s.sphere().radius()).clamp(-1.0, 1.0).asin();
            (
                local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU),
                latitude,
            )
        }
        SurfaceGeometry::Cylinder(c) => {
            let local = c.cylinder().frame().to_local(point);
            (
                local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU),
                local.z,
            )
        }
        SurfaceGeometry::Cone(c) => {
            ogeom_math::elementary::cone_parameters(&c.cone(), point, tol).ok()?
        }
        SurfaceGeometry::Torus(t) => {
            ogeom_math::elementary::torus_parameters(&t.torus(), point, tol).ok()?
        }
        _ => return None,
    };
    // One polish step against the curve point, so parameter rounding in the
    // inversion does not survive into the result, and the gap comes with it,
    // because the polish clamps into the surface's extents and the gap is
    // what says whether the clamped answer still touches the curve.
    polish(curve, surface, on_curve, guess, tol)
}

// --- general -----------------------------------------------------------------

fn general(
    curve: &Curve,
    surface: &SurfaceGeometry,
    cells: &[Cell],
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> OgeomResult<CurveSurfaceIntersection> {
    let (lo, hi) = curve.domain();

    let mut points = Vec::with_capacity(options.samples + 1);
    for i in 0..=options.samples {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / options.samples as f64;
        if let Ok(p) = curve.point_at(t, tol) {
            points.push((t, p));
        }
    }

    // Every segment's box lies within the samples' box, so a cell that
    // misses that box at its own margin meets no segment, and only the
    // cells near the curve are carried into the segment loop, in order.
    let first = points.first().map_or(Point::ORIGIN, |p| p.1);
    let (mut low, mut high) = (first, first);
    for &(_, p) in &points {
        low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
        high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
    }
    let near: Vec<&Cell> = cells
        .iter()
        .filter(|cell| segment_near_cell(low, high, cell, options.gap.max(cell.sag)))
        .collect();

    let mut crossings: Vec<Piercing> = Vec::new();
    for pair in points.windows(2) {
        let (t0, p0) = pair[0];
        let (t1, p1) = pair[1];
        for &cell in &near {
            // Near the cell within the surface's own bow from it: a curve
            // crossing the surface in the gap between the flat cell and
            // the curved patch it stands for (a ray starting a few microns
            // from the wall it leaves by) meets no cell, and is seeded by
            // its nearness instead.
            if !segment_near_cell(p0, p1, cell, options.gap.max(cell.sag)) {
                continue;
            }
            if segment_meets_triangle(p0, p1, cell.corners).is_none()
                && !(cell.sag > options.gap && segment_near_cell(p0, p1, cell, cell.sag))
            {
                continue;
            }
            // Newton from where the segment meets the cell, not from the
            // cell's corner: on a wall bowing a few hundredths of a
            // millimetre over a cell, the corner can stand far enough off
            // that the first step leaves the chart, and clamped at its edge
            // the solve stalls there. The corner stays a second try.
            let (near_t, near_uv) = seed_in(cell, p0, p1, t0, t1);
            let Some(found) = [(near_t, near_uv), (f64::midpoint(t0, t1), cell.at)]
                .into_iter()
                .filter_map(|(t, uv)| polish(curve, surface, t, uv, tol))
                .find(|found| found.gap <= options.gap)
            else {
                continue;
            };
            let reach = tol.confusion() * 100.0;
            if !crossings
                .iter()
                .any(|c| c.point.distance(found.point) <= reach)
            {
                crossings.push(found);
            }
        }
    }
    crossings.sort_by(|a, b| {
        a.on_curve
            .partial_cmp(&b.on_curve)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    Ok(gathered(curve, surface, crossings, options, tol))
}

/// Piercings gathered into the contacts they describe.
///
/// Newton lands a little apart each time it is started near a tangency,
/// where the gap touches zero without crossing it, and all along a stretch
/// of curve lying in the surface; a radius of rounding merges neither. So
/// neighbouring piercings between which the curve never leaves the surface
/// are one contact: a run shorter than the curve's sampling step is one
/// tangency, kept once (where the gap is least), and a longer one is a
/// stretch lying in the surface.
fn gathered(
    curve: &Curve,
    surface: &SurfaceGeometry,
    crossings: Vec<Piercing>,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> CurveSurfaceIntersection {
    let stays_on = |x: &Piercing, y: &Piercing| -> bool {
        (1..=3).all(|k| {
            let t = x.on_curve + (y.on_curve - x.on_curve) * f64::from(k) / 4.0;
            let Ok(p) = curve.point_at(t, tol) else {
                return false;
            };
            // From either end: by a pole or an apex one side's chart
            // turns too fast for the other's seed.
            [x.on_surface, y.on_surface].into_iter().any(|seed| {
                crate::march::nearest_on(surface, seed, p, tol)
                    .is_some_and(|(_, q)| q.distance(p) <= options.gap)
            })
        })
    };
    let (lo, hi) = curve.domain();
    #[allow(clippy::cast_precision_loss)]
    let step = (hi - lo).abs() / options.samples.max(1) as f64;
    let mut runs: Vec<Vec<Piercing>> = Vec::new();
    for crossing in crossings {
        match runs.last_mut() {
            Some(run) if run.last().is_some_and(|last| stays_on(last, &crossing)) => {
                run.push(crossing);
            }
            _ => runs.push(vec![crossing]),
        }
    }
    let mut out = CurveSurfaceIntersection {
        crossings: Vec::new(),
        lying: Vec::new(),
    };
    for run in runs {
        let (Some(first), Some(last)) = (run.first(), run.last()) else {
            continue;
        };
        if run.len() > 1 && (last.on_curve - first.on_curve).abs() >= step {
            out.lying.push((first.on_curve, last.on_curve));
        } else if let Some(best) = run.iter().min_by(|a, b| {
            a.gap
                .partial_cmp(&b.gap)
                .unwrap_or(core::cmp::Ordering::Equal)
        }) {
            out.crossings.push(*best);
        }
    }
    out
}

/// How many seed cells to lay along each direction of a surface: the
/// asked grid, and for a spline at least two per knot span. A patch swept
/// several turns round an axis (a thread's flank) spans dozens of knots
/// along its length, and a grid of the asked size lays flat cells a turn's
/// fraction wide whose chords stand a tenth of a millimetre off the wall;
/// a curve crossing the wall inside that gap meets no cell and is missed.
fn seeding(surface: &SurfaceGeometry, grid: usize) -> (usize, usize) {
    const CAP: usize = 1024;
    let SurfaceGeometry::BSpline(spline) = surface else {
        return (grid, grid);
    };
    let spans = |knots: &ogeom_math::KnotVector| knots.distinct().len().saturating_sub(1);
    (
        grid.max(2 * spans(spline.u_knots())).min(CAP.max(grid)),
        grid.max(2 * spans(spline.v_knots())).min(CAP.max(grid)),
    )
}

/// Where to start Newton for a segment near a cell: the point of the cell
/// the segment passes through, or failing that the one nearest its middle,
/// with its parameters on the curve and on the surface read off the cell's
/// corners.
fn seed_in(cell: &Cell, p0: Point, p1: Point, t0: f64, t1: f64) -> (f64, (f64, f64)) {
    let [a, b, c] = cell.corners;
    let (t, at) = segment_meets_triangle(p0, p1, cell.corners).map_or_else(
        || (f64::midpoint(t0, t1), p0.midpoint(p1)),
        |x| {
            let length = p0.distance(p1);
            let f = if length > 0.0 {
                p0.distance(x) / length
            } else {
                0.5
            };
            (t0 + (t1 - t0) * f, x)
        },
    );
    // Barycentric weights of the point's foot in the cell's plane, pulled
    // back inside the cell.
    let (e1, e2, d) = (b - a, c - a, at - a);
    let (d11, d12, d22) = (e1.dot(e1), e1.dot(e2), e2.dot(e2));
    let (d1, d2) = (d.dot(e1), d.dot(e2));
    let det = d11 * d22 - d12 * d12;
    let (mut wb, mut wc) = if det > 0.0 {
        ((d22 * d1 - d12 * d2) / det, (d11 * d2 - d12 * d1) / det)
    } else {
        (1.0 / 3.0, 1.0 / 3.0)
    };
    wb = wb.clamp(0.0, 1.0);
    wc = wc.clamp(0.0, 1.0);
    if wb + wc > 1.0 {
        let sum = wb + wc;
        wb /= sum;
        wc /= sum;
    }
    let wa = 1.0 - wb - wc;
    let [pa, pb, pc] = cell.params;
    (
        t,
        (
            wa * pa.0 + wb * pb.0 + wc * pc.0,
            wa * pa.1 + wb * pb.1 + wc * pc.1,
        ),
    )
}

/// Whether a segment's box comes near a cell's.
fn segment_near_cell(a: Point, b: Point, cell: &Cell, margin: f64) -> bool {
    let low = Point::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z));
    let high = Point::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z));
    low.x <= cell.high.x + margin
        && cell.low.x <= high.x + margin
        && low.y <= cell.high.y + margin
        && cell.low.y <= high.y + margin
        && low.z <= cell.high.z + margin
        && cell.low.z <= high.z + margin
}

/// Newton on the well-posed system `C(t) = S(u, v)`.
fn polish(
    curve: &Curve,
    surface: &SurfaceGeometry,
    seed_t: f64,
    seed_uv: (f64, f64),
    tol: Tolerances,
) -> Option<Piercing> {
    let clamp_t = |t: f64| {
        let (lo, hi) = curve.domain();
        if curve.is_periodic() {
            let span = hi - lo;
            if span > 0.0 {
                return lo + (t - lo).rem_euclid(span);
            }
        }
        t.clamp(lo, hi)
    };
    let clamp_uv = |u: f64, v: f64| {
        let ((ua, ub), (va, vb)) = surface.domain();
        let fold = |x: f64, lo: f64, hi: f64, periodic: bool| {
            if periodic {
                let span = hi - lo;
                if span > 0.0 {
                    return lo + (x - lo).rem_euclid(span);
                }
            }
            x.clamp(lo, hi)
        };
        (
            fold(u, ua, ub, surface.is_periodic_u()),
            fold(v, va, vb, surface.is_periodic_v()),
        )
    };

    // The gap and its Jacobian apart: the damped step measures many trial
    // points and steps onto few, and a spline's derivatives cost more than
    // its point.
    let gap_at = |x: &[f64; 3]| {
        let t = clamp_t(x[0]);
        let (u, v) = clamp_uv(x[1], x[2]);
        let (Ok(pc), Ok(ps)) = (curve.point_at(t, tol), surface.point_at(u, v, tol)) else {
            // Nowhere to measure from: infinite, so the damped step backs off.
            return [f64::INFINITY; 3];
        };
        let gap = pc - ps;
        [gap.x, gap.y, gap.z]
    };
    let jacobian_at = |x: &[f64; 3]| {
        let t = clamp_t(x[0]);
        let (u, v) = clamp_uv(x[1], x[2]);
        let (Ok(dc), Ok((du, dv))) = (curve.d1_at(t, tol), surface.d1_at(u, v, tol)) else {
            return None;
        };
        Some([
            [dc.x, -du.x, -dv.x],
            [dc.y, -du.y, -dv.y],
            [dc.z, -du.z, -dv.z],
        ])
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system_fixed_lazy(
        gap_at,
        jacobian_at,
        [seed_t, seed_uv.0, seed_uv.1],
        criteria,
    )
    .ok()?;
    let t = clamp_t(found.0[0]);
    let (u, v) = clamp_uv(found.0[1], found.0[2]);
    let pc = curve.point_at(t, tol).ok()?;
    let ps = surface.point_at(u, v, tol).ok()?;
    Some(Piercing {
        on_curve: t,
        on_surface: (u, v),
        point: pc,
        gap: pc.distance(ps),
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_geom::{
        BSplineCurve, CircleCurve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface,
    };
    use ogeom_math::{Circle, Cylinder, Direction, Frame, KnotVector, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(Point::ORIGIN, radius, T).unwrap()).into()
    }

    fn cylinder(radius: f64, height: (f64, f64)) -> SurfaceGeometry {
        CylinderSurface::new(Cylinder::new(Frame::WORLD, radius, T).unwrap(), height)
            .unwrap()
            .into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-6.0, 6.0),
            (-6.0, 6.0),
        )
        .unwrap()
        .into()
    }

    fn segment(from: Point, to: Point) -> Curve {
        LineCurve::segment(from, to, T).unwrap().into()
    }

    #[test]
    fn a_line_through_a_sphere_pierces_it_where_the_quadratic_says() {
        let ball = sphere(2.0);
        let ray = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        let found =
            intersect_curve_surface(&ray, &ball, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        assert!(
            found.crossings[0]
                .point
                .is_equal(Point::new(-2.0, 0.0, 0.0), T)
        );
        assert!(
            found.crossings[1]
                .point
                .is_equal(Point::new(2.0, 0.0, 0.0), T)
        );
        for hit in &found.crossings {
            assert!(hit.gap < 1e-12);
            // The surface parameters evaluate back onto the point.
            let lifted = ball
                .point_at(hit.on_surface.0, hit.on_surface.1, T)
                .unwrap();
            assert!(lifted.is_equal(hit.point, T));
        }

        // Tangent: one root. Missing: none.
        let grazing = segment(Point::new(-5.0, 0.0, 2.0), Point::new(5.0, 0.0, 2.0));
        assert_eq!(
            intersect_curve_surface(&grazing, &ball, CurveSurfaceOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            1
        );
        let missing = segment(Point::new(-5.0, 0.0, 3.0), Point::new(5.0, 0.0, 3.0));
        assert!(
            intersect_curve_surface(&missing, &ball, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_line_through_a_cylinder_respects_its_height() {
        let drum = cylinder(2.0, (-1.0, 1.0));
        // Crosses the infinite cylinder at z = 0: inside the height, two hits.
        let level = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        assert_eq!(
            intersect_curve_surface(&level, &drum, CurveSurfaceOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            2
        );
        // Crosses at z = 3: the unbounded geometry meets it, this surface
        // does not reach there.
        let high = segment(Point::new(-5.0, 0.0, 3.0), Point::new(5.0, 0.0, 3.0));
        assert!(
            intersect_curve_surface(&high, &drum, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_line_lying_in_a_plane_is_an_overlap_not_a_crossing_list() {
        let ground = plane(Point::ORIGIN, Vector::Z);
        let lying = segment(Point::new(-3.0, 1.0, 0.0), Point::new(3.0, 1.0, 0.0));
        let found =
            intersect_curve_surface(&lying, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert!(found.crossings.is_empty());
        assert_eq!(found.lying.len(), 1);

        let crossing = segment(Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0));
        let found =
            intersect_curve_surface(&crossing, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.is_equal(Point::ORIGIN, T));

        let parallel = segment(Point::new(-3.0, 0.0, 1.0), Point::new(3.0, 0.0, 1.0));
        assert!(
            intersect_curve_surface(&parallel, &ground, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_circle_pierces_a_plane_twice_through_the_general_path() {
        // A circle in the xz plane against the ground: no analytic case
        // handles circle/plane here, so this is the seeded Newton path, and
        // the answer is known exactly anyway.
        let ring: Curve = CircleCurve::new(
            Circle::new(
                Frame::new(Point::new(0.0, 0.0, 0.0), -Direction::Y, Direction::X, T).unwrap(),
                2.0,
                T,
            )
            .unwrap(),
        )
        .into();
        let ground = plane(Point::ORIGIN, Vector::Z);
        let found =
            intersect_curve_surface(&ring, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        for hit in &found.crossings {
            assert!(hit.gap < 1e-9);
            assert!(hit.point.z.abs() < 1e-9);
            assert!((hit.point.to_vector().magnitude() - 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn a_prepared_surface_answers_as_the_one_shot_query() {
        // A half-pipe spline strip, which only the seeded path can answer,
        // asked about many lines in turn: the stored sampling must give
        // every one the answer a fresh sampling gives, piercings and all.
        use ogeom_geom::BSplineSurface;
        use ogeom_math::ControlGrid;
        let mut points = Vec::new();
        for i in 0..7 {
            let a = core::f64::consts::PI * f64::from(i) / 6.0;
            for j in 0..2 {
                points.push(Point::new(2.0 * a.cos(), 3.0 * f64::from(j), 2.0 * a.sin()));
            }
        }
        let strip: SurfaceGeometry = BSplineSurface::new(
            KnotVector::clamped_uniform(3, 7).unwrap(),
            KnotVector::clamped_uniform(1, 2).unwrap(),
            &ControlGrid::new(points, 7, 2).unwrap(),
            T,
        )
        .unwrap()
        .into();
        let options = CurveSurfaceOptions::default();
        let prepared = PreparedSurface::new(strip.clone(), options, T);
        let mut pierced = 0;
        for k in 0..12 {
            let x = -2.5 + 5.0 * f64::from(k) / 11.0;
            let lines: [Curve; 2] = [
                LineCurve::new(ogeom_math::Axis::new(
                    Point::new(x, 1.5, 0.0),
                    Direction::new(Vector::new(0.1, 0.2, 1.0), T).unwrap(),
                ))
                .into(),
                segment(Point::new(x, -1.0, -1.0), Point::new(-x, 4.0, 3.0)),
            ];
            for line in &lines {
                let once = intersect_curve_surface(line, &strip, options, T).unwrap();
                pierced += once.crossings.len();
                assert_eq!(prepared.intersect(line).unwrap(), once, "line {k}");
            }
        }
        assert!(pierced > 0, "some lines pass through the strip");
    }

    #[test]
    fn a_spline_through_a_sphere_is_found_and_polished() {
        // A spline wandering through the ball: piercings with no closed form
        // anywhere, verified implicitly: each reported point is on the
        // sphere to the gap it claims.
        let wander: Curve = BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![
                Point::new(-4.0, -1.0, -1.0),
                Point::new(-1.0, 2.0, 1.0),
                Point::new(1.0, -2.0, -1.0),
                Point::new(4.0, 1.0, 1.0),
            ],
            T,
        )
        .unwrap()
        .into();
        let ball = sphere(2.0);
        let found =
            intersect_curve_surface(&wander, &ball, CurveSurfaceOptions::default(), T).unwrap();
        assert!(!found.crossings.is_empty(), "the spline passes through");
        for hit in &found.crossings {
            assert!(hit.gap < 1e-9);
            let SurfaceGeometry::Sphere(s) = &ball else {
                unreachable!()
            };
            assert!(s.sphere().distance_to(hit.point).abs() < 1e-9);
        }
    }

    #[test]
    fn unusable_options_are_refused() {
        let ball = sphere(1.0);
        let ray = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        for options in [
            CurveSurfaceOptions {
                samples: 1,
                ..CurveSurfaceOptions::default()
            },
            CurveSurfaceOptions {
                grid: 1,
                ..CurveSurfaceOptions::default()
            },
            CurveSurfaceOptions {
                gap: 0.0,
                ..CurveSurfaceOptions::default()
            },
        ] {
            assert!(intersect_curve_surface(&ray, &ball, options, T).is_err());
        }
    }

    /// The general path (a torus has no closed form against a line) answers
    /// each tangency once, and a curve lying in the surface as the stretch
    /// it lies along rather than as hundreds of piercings.
    /// A torus far down a line is met where it is: its quartic, solved where
    /// the line starts, has roots hundreds of units out, and the eigenvalues
    /// lost them.
    #[test]
    fn a_torus_far_down_a_line_is_still_met() {
        use ogeom_geom::TorusSurface;
        use ogeom_math::Torus;
        let options = CurveSurfaceOptions::default();
        let tilted = Frame::new(
            Point::new(0.0, 785.0, -140.0),
            Direction::new(Vector::new(0.0, -0.999_390_827, 0.034_899_497), T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let torus: SurfaceGeometry =
            TorusSurface::new(Torus::new(tilted, 120.0, 12.0, T).unwrap()).into();
        let w = Vector::new(1.0, 1.0, 1.0) / 3.0_f64.sqrt();
        let from = Point::new(-160.0, 531.0, -320.0);
        let line = segment(from, from + w * 800.0);
        let found = intersect_curve_surface(&line, &torus, options, T).unwrap();
        let cells = sample_by(&torus, seeding(&torus, options.grid), T);
        let general = general(&line, &torus, &cells, options, T).unwrap();
        assert!(!general.crossings.is_empty());
        assert_eq!(found.crossings.len(), general.crossings.len());
        for (a, b) in found.crossings.iter().zip(&general.crossings) {
            assert!(
                a.point.distance(b.point) < 1e-6,
                "{:?} against {:?}",
                a.point,
                b.point
            );
        }
    }

    /// A line through a torus's middle crosses the tube four times, and one
    /// through a cone's axis crosses its wall twice, at the closed forms'
    /// points.
    #[test]
    fn lines_cross_a_torus_and_a_cone_where_the_closed_forms_say() {
        use ogeom_geom::{ConeSurface, TorusSurface};
        use ogeom_math::{Cone, Torus};
        let options = CurveSurfaceOptions::default();
        let torus: SurfaceGeometry =
            TorusSurface::new(Torus::new(Frame::WORLD, 60.0, 20.0, T).unwrap()).into();
        let across = segment(Point::new(-100.0, 0.0, 0.0), Point::new(100.0, 0.0, 0.0));
        let found = intersect_curve_surface(&across, &torus, options, T).unwrap();
        let xs: Vec<f64> = found.crossings.iter().map(|c| c.point.x).collect();
        assert_eq!(xs.len(), 4, "{xs:?}");
        for (got, want) in xs.iter().zip([-80.0, -40.0, 40.0, 80.0]) {
            assert!((got - want).abs() < 1e-9, "{xs:?}");
        }
        let cone: SurfaceGeometry = ConeSurface::new(
            Cone::new(Frame::WORLD, 10.0, core::f64::consts::FRAC_PI_4, T).unwrap(),
            (0.0, 20.0),
        )
        .unwrap()
        .into();
        let level = segment(Point::new(-50.0, 0.0, 5.0), Point::new(50.0, 0.0, 5.0));
        let found = intersect_curve_surface(&level, &cone, options, T).unwrap();
        let xs: Vec<f64> = found.crossings.iter().map(|c| c.point.x).collect();
        assert_eq!(xs.len(), 2, "{xs:?}");
        for (got, want) in xs.iter().zip([-15.0, 15.0]) {
            assert!((got - want).abs() < 1e-9, "{xs:?}");
        }
    }

    /// Lines that miss a torus through its hole or along its axis meet
    /// nothing; one grazing the outer equator touches it once, at the
    /// tangent point; a segment ending on the tube keeps its end crossings.
    #[test]
    fn lines_through_the_hole_along_the_axis_and_grazing_a_torus() {
        use ogeom_geom::TorusSurface;
        use ogeom_math::Torus;
        let options = CurveSurfaceOptions::default();
        let tilted = Frame::new(
            Point::new(3.0, -2.0, 1.0),
            Direction::new(Vector::new(0.1, 0.2, 1.0), T).unwrap(),
            Direction::new(Vector::new(1.0, 0.0, -0.1), T).unwrap(),
            T,
        )
        .unwrap();
        let torus: SurfaceGeometry =
            TorusSurface::new(Torus::new(tilted, 60.0, 20.0, T).unwrap()).into();
        let local = |x: f64, y: f64, z: f64| tilted.to_world(Point::new(x, y, z));
        let meets = |from: Point, to: Point| {
            intersect_curve_surface(&segment(from, to), &torus, options, T)
                .unwrap()
                .crossings
        };
        assert!(meets(local(0.0, 0.0, -100.0), local(0.0, 0.0, 100.0)).is_empty());
        assert!(meets(local(5.0, 3.0, -100.0), local(-15.0, 3.0, 100.0)).is_empty());
        let grazing = meets(local(-100.0, 80.0, 0.0), local(100.0, 80.0, 0.0));
        assert_eq!(grazing.len(), 1, "{grazing:?}");
        assert!(grazing[0].point.distance(local(0.0, 80.0, 0.0)) < 1e-6);
        let ends = meets(local(-80.0, 0.0, 0.0), local(80.0, 0.0, 0.0));
        let xs: Vec<f64> = ends.iter().map(|c| tilted.to_local(c.point).x).collect();
        assert_eq!(xs.len(), 4, "{xs:?}");
        for (got, want) in xs.iter().zip([-80.0, -40.0, 40.0, 80.0]) {
            assert!((got - want).abs() < 1e-9, "{xs:?}");
        }
    }

    #[test]
    fn tangencies_come_back_once_and_lying_curves_as_stretches() {
        use ogeom_geom::TorusSurface;
        use ogeom_math::Torus;
        let torus: SurfaceGeometry =
            TorusSurface::new(Torus::new(Frame::WORLD, 60.0, 20.0, T).unwrap()).into();
        let options = CurveSurfaceOptions::default();
        // Along the tube's top: tangent at two points.
        let top = segment(Point::new(-100.0, 0.0, 20.0), Point::new(100.0, 0.0, 20.0));
        let found = intersect_curve_surface(&top, &torus, options, T).unwrap();
        assert_eq!(found.crossings.len(), 2, "{:?}", found.crossings);
        assert!(found.lying.is_empty());
        // Touching the inner equator and crossing the outer twice.
        let inner = segment(Point::new(-100.0, 40.0, 0.0), Point::new(100.0, 40.0, 0.0));
        let found = intersect_curve_surface(&inner, &torus, options, T).unwrap();
        assert_eq!(found.crossings.len(), 3, "{:?}", found.crossings);
        // A parallel of the torus lies in it the whole way round.
        let parallel: Curve = CircleCurve::new(
            Circle::new(
                Frame::new(Point::new(0.0, 0.0, 20.0), Direction::Z, Direction::X, T).unwrap(),
                60.0,
                T,
            )
            .unwrap(),
        )
        .into();
        let found = intersect_curve_surface(&parallel, &torus, options, T).unwrap();
        assert!(found.crossings.is_empty(), "{:?}", found.crossings);
        assert_eq!(found.lying.len(), 1);
        let (from, to) = found.lying[0];
        assert!((to - from).abs() > 6.0, "{from} to {to}");
    }
}
