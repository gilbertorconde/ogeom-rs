//! P3's general vertex blend for rounds: a three- or four-sided
//! transfinite patch through the strips' ends (ADR 0001, §4.4, as amended).
//!
//! Each side of the hole is a strip's end curve, and the strip's surface
//! its support. Its ribbon leaves the curve along a cross tangent in the
//! support's tangent plane and bends as the support does that way, so that
//! near its side it keeps to the support's own continuation to second
//! order (a straight ribbon off a strip tangent to a face crosses that
//! face). The cross tangent's ends are fixed by the corners: at the side's
//! start the previous side's tangent turned back, at its end the next
//! side's tangent, so each ribbon agrees to first order with its
//! neighbours' curves. Between them its components are cubic, their slopes
//! at each corner making the two ribbons' twists there agree in the face
//! both strips touch; only the strips' curvatures across that face, which
//! nothing mends, still differ. The patch over the triangle of
//! barycentrics `l` is Kato's blend of the ribbons,
//!
//! ```text
//! S(l) = sum_i mu_i R_i(s_i, d_i),  mu_i = prod_{j != i} d_j^2 / sum_k prod_{j != k} d_j^2,
//! ```
//!
//! with `d_i = l_{i+2}` (nought on side `i`) and `s_i = l_{i+1} / (l_i + l_{i+1})`
//! the side's parameter. Near side `i` every other weight falls as `d_i^2`,
//! so `S = R_i + O(d_i^2)`: the patch takes each side's curve and its cross
//! tangent, so it is tangent to the support all along the side, however
//! the supports' twists at a corner disagree. This is what a least-squares
//! fill through such corners cannot do. Four sides blend over the unit
//! square with `d = (v, 1 - u, 1 - v, u)` and `s = (u, v, 1 - u, 1 - v)`,
//! their distances cubed in the weights: at a corner the two ribbons of the
//! sides not meeting there stand off it, and squared their weights let
//! that enter the patch at second order, where the two ribbons meeting
//! there keep to their strips and so to the faces those touch (at a
//! pyramid's apex the exact patch stood past a face beside a corner);
//! cubed, it enters at third.
//!
//! The patch is not a tensor product, so the face carries it interpolated:
//! a clamped bicubic B-spline taking it at [`NET`] by [`NET`] points of the
//! unit square, and its derivatives across the square's sides at each,
//! three sides' square mapped onto the triangle with its third corner as a
//! pole, four sides' square taken as it is. What the interpolated surface
//! achieves is measured on it, never assumed from the exact patch: its gap
//! to each side's curve, its angle to each side's support, and how far it
//! stands from the exact patch inside.

use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{BSpline2d, BSplineSurface, Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{ControlGrid, KnotVector, Point, Point2, Vector};
use ogeom_topo::{EdgeData, Location, Model, Shape};

/// The interpolation points along each direction of the unit square.
pub(crate) const NET: usize = 33;

/// One side of a three- or four-sided hole: an edge's curve walked from one corner
/// to the next, and the surface the patch meets tangent along it.
#[derive(Debug, Clone)]
pub(crate) struct HoleSide {
    /// The edge, as a node (its own orientation forward).
    pub edge: Shape,
    /// The edge's curve, and the range of it the edge covers.
    pub curve: Curve,
    /// The range the edge covers.
    pub range: (f64, f64),
    /// Whether the walk runs with the curve's parameter.
    pub forward: bool,
    /// The surface across the side.
    pub support: SurfaceGeometry,
}

/// A side's ribbon: its curve over `[0, 1]` and its cross tangent.
#[derive(Debug, Clone)]
struct Ribbon {
    side: HoleSide,
    /// Into the hole, the support's tangent plane's direction square to
    /// the curve, as `+1` or `-1` times `normal x tangent`.
    sign: f64,
    /// The cross tangent's components along the curve's tangent and the
    /// direction into the hole, at the side's start and end.
    ends: [(f64, f64); 2],
    /// Those components' derivatives in `s` at the side's start and end.
    slopes: [(f64, f64); 2],
}

impl Ribbon {
    fn parameter(&self, s: f64) -> f64 {
        let (a, b) = self.side.range;
        if self.side.forward {
            a + (b - a) * s
        } else {
            b - (b - a) * s
        }
    }

    fn point(&self, s: f64, tol: Tolerances) -> OgeomResult<Point> {
        self.side.curve.point_at(self.parameter(s), tol)
    }

    /// The curve's derivative in `s`.
    fn tangent(&self, s: f64, tol: Tolerances) -> OgeomResult<Vector> {
        let (a, b) = self.side.range;
        let scale = if self.side.forward { b - a } else { a - b };
        Ok(self.side.curve.d1_at(self.parameter(s), tol)? * scale)
    }

    /// The unit tangent and the unit direction into the hole in the
    /// support's tangent plane, the latter of sign `sign`.
    fn frame(&self, s: f64, sign: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let p = self.point(s, tol)?;
        let t = self.tangent(s, tol)?.normalized(tol)?;
        let n = support_normal(&self.side.support, p, tol)?;
        let across = n.cross(t).normalized(tol)? * sign;
        Ok((t, across))
    }

    /// The cross tangent at `s`: its components carried between the
    /// corners' by cubic Hermite interpolation of their values and slopes,
    /// in the frame there.
    fn cross(&self, s: f64, tol: Tolerances) -> OgeomResult<Vector> {
        let (t, across) = self.frame(s, self.sign, tol)?;
        let [(a0, b0), (a1, b1)] = self.ends;
        let [(da0, db0), (da1, db1)] = self.slopes;
        let (s2, s3) = (s * s, s * s * s);
        let h = [
            2.0 * s3 - 3.0 * s2 + 1.0,
            s3 - 2.0 * s2 + s,
            -2.0 * s3 + 3.0 * s2,
            s3 - s2,
        ];
        let along = h[0] * a0 + h[1] * da0 + h[2] * a1 + h[3] * da1;
        let into = h[0] * b0 + h[1] * db0 + h[2] * b1 + h[3] * db1;
        Ok(t * along + across * into)
    }

    /// The ribbon `d` across from the side at `s`: off the curve along the
    /// cross tangent, bending as the support does that way, so that near
    /// its side it keeps to the support's own continuation to second
    /// order.
    fn at(&self, s: f64, d: f64, tol: Tolerances) -> OgeomResult<Point> {
        let p = self.point(s, tol)?;
        let t = self.cross(s, tol)?;
        let bend = support_bend(&self.side.support, p, t, tol)?;
        Ok(p + t * d + bend * (0.5 * d * d))
    }

    /// The derivative in `s` of the frame's two directions, by differences
    /// a hundred-thousandth apart, to second order (one-sided at the ends).
    fn frame_rate(&self, s: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let h = 1e-5;
        let at = |x: f64| self.frame(x, self.sign, tol);
        let (steps, scale) = if s <= 0.0 {
            ([0.0, h, 2.0 * h], 1.0)
        } else {
            ([1.0, 1.0 - h, 1.0 - 2.0 * h], -1.0)
        };
        let (f0, f1, f2) = (at(steps[0])?, at(steps[1])?, at(steps[2])?);
        let d = |a: Vector, b: Vector, c: Vector| (b * 4.0 - a * 3.0 - c) * (scale / (2.0 * h));
        Ok((d(f0.0, f1.0, f2.0), d(f0.1, f1.1, f2.1)))
    }
}

/// How a curve on `support` through `p` with velocity `t` (tangent to it)
/// accelerates across it: the normal curvature in `t`'s direction times
/// `|t|^2`, toward the side the support curves to. A plane does not.
fn support_bend(
    support: &SurfaceGeometry,
    p: Point,
    t: Vector,
    tol: Tolerances,
) -> OgeomResult<Vector> {
    match support {
        SurfaceGeometry::Plane(_) => Ok(Vector::ZERO),
        SurfaceGeometry::Cylinder(c) => {
            let frame = c.cylinder().frame();
            let z = frame.z().vector();
            let off = p - frame.origin();
            let radial = (off - z * off.dot(z)).normalized(tol)?;
            let round = t.dot(t) - t.dot(z).powi(2);
            Ok(radial * (-round / c.cylinder().radius()))
        }
        SurfaceGeometry::Sphere(sp) => {
            let radial = (p - sp.sphere().frame().origin()).normalized(tol)?;
            Ok(radial * (-t.dot(t) / sp.sphere().radius()))
        }
        _ => ogeom_bail!(
            Construction,
            "a vertex patch's side runs along a surface whose curvature it cannot read"
        ),
    }
}

/// The unit normal of `support` at `p`, a point on it, up to sign.
fn support_normal(support: &SurfaceGeometry, p: Point, tol: Tolerances) -> OgeomResult<Vector> {
    match support {
        SurfaceGeometry::Plane(plane) => Ok(plane.plane().normal().vector()),
        SurfaceGeometry::Cylinder(c) => {
            let frame = c.cylinder().frame();
            let off = p - frame.origin();
            let z = frame.z().vector();
            (off - z * off.dot(z)).normalized(tol)
        }
        SurfaceGeometry::Sphere(s) => (p - s.sphere().frame().origin()).normalized(tol),
        _ => ogeom_bail!(
            Construction,
            "a vertex patch's side runs along a surface whose normal it cannot read"
        ),
    }
}

/// The exact three- or four-sided patch over a hole: Kato's blend of its
/// sides' ribbons.
#[derive(Debug, Clone)]
pub(crate) struct Transfinite {
    ribbons: Vec<Ribbon>,
}

impl Transfinite {
    /// The patch over the hole `sides`, walked corner to corner, each side
    /// ending where the next starts.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where there are other than three or four sides, the sides do not
    /// close, or a support's normal cannot be read.
    pub(crate) fn new(sides: &[HoleSide], tol: Tolerances) -> OgeomResult<Self> {
        let n = sides.len();
        if !(3..=4).contains(&n) {
            ogeom_bail!(
                Construction,
                "a vertex patch takes three or four sides, not {n}"
            );
        }
        let mut ribbons: Vec<Ribbon> = sides
            .iter()
            .map(|side| Ribbon {
                side: side.clone(),
                sign: 1.0,
                ends: [(0.0, 0.0); 2],
                slopes: [(0.0, 0.0); 2],
            })
            .collect();
        let reach = tol.confusion() * 100.0;
        for i in 0..n {
            let next = (i + 1) % n;
            if ribbons[i]
                .point(1.0, tol)?
                .distance(ribbons[next].point(0.0, tol)?)
                > reach
            {
                ogeom_bail!(Construction, "a vertex patch's sides do not close");
            }
        }
        // The hole's middle, for each side's way into it: read at the
        // side's middle, not at a corner, where the hole may turn straight
        // on and its cross tangent run along the side.
        let mut middle = Vector::ZERO;
        for r in &ribbons {
            middle += r.point(0.0, tol)?.to_vector();
        }
        #[expect(clippy::cast_precision_loss, reason = "three or four sides")]
        let middle = Point::from_vector(middle / n as f64);
        for i in 0..n {
            let (prev, next) = ((i + n - 1) % n, (i + 1) % n);
            let start = -ribbons[prev].tangent(1.0, tol)?;
            let end = ribbons[next].tangent(0.0, tol)?;
            let (_, inward) = ribbons[i].frame(0.5, 1.0, tol)?;
            let towards = middle - ribbons[i].point(0.5, tol)?;
            let sign = if inward.dot(towards) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let (t0, across0) = ribbons[i].frame(0.0, sign, tol)?;
            let (t1, across1) = ribbons[i].frame(1.0, sign, tol)?;
            ribbons[i].sign = sign;
            ribbons[i].ends = [
                (start.dot(t0), start.dot(across0)),
                (end.dot(t1), end.dot(across1)),
            ];
        }
        // The twists at each corner, from the side ending there and the side
        // starting there, agree in the face both strips touch: with `a`
        // along the next side and `b` back along this one, the patch's
        // mixed derivative is `-T_i'(1)` by one ribbon and `T_{i+1}'(0)` by
        // the other. Their parts in the face are the components' slopes and
        // the frames' turning; the slopes take half each of what the
        // turning leaves. What remains is the two strips' curvatures across
        // the face, which no slope mends.
        for i in 0..n {
            let next = (i + 1) % n;
            let (t_end, across_end) = ribbons[i].frame(1.0, ribbons[i].sign, tol)?;
            let (t_start, across_start) = ribbons[next].frame(0.0, ribbons[next].sign, tol)?;
            let turning = |r: &Ribbon, k: usize, rate: (Vector, Vector)| {
                let (a, b) = r.ends[k];
                rate.0 * a + rate.1 * b
            };
            let q = turning(&ribbons[i], 1, ribbons[i].frame_rate(1.0, tol)?)
                + turning(&ribbons[next], 0, ribbons[next].frame_rate(0.0, tol)?);
            let normal = t_end.cross(across_end).normalized(tol)?;
            let in_face = q - normal * normal.dot(q);
            let half = in_face * -0.5;
            ribbons[i].slopes[1] = (half.dot(t_end), half.dot(across_end));
            ribbons[next].slopes[0] = (half.dot(t_start), half.dot(across_start));
        }
        Ok(Self { ribbons })
    }

    /// A three-sided patch at barycentrics `l` (summing to one, each in
    /// `[0, 1]`).
    ///
    /// # Errors
    ///
    /// As evaluating the sides' curves.
    pub(crate) fn point(&self, l: [f64; 3], tol: Tolerances) -> OgeomResult<Point> {
        let d = [l[2], l[0], l[1]];
        let weight = |k: usize| -> f64 {
            (0..3)
                .filter(|j| *j != k)
                .map(|j| d[j] * d[j])
                .product::<f64>()
        };
        let w = [weight(0), weight(1), weight(2)];
        let total: f64 = w.iter().sum();
        if total <= f64::MIN_POSITIVE {
            // A corner: the side starting there.
            let corner = (0..3).fold(0, |best, k| if l[k] > l[best] { k } else { best });
            return self.ribbons[corner].point(0.0, tol);
        }
        let mut sum = Vector::ZERO;
        for i in 0..3 {
            if w[i] == 0.0 {
                continue;
            }
            let along = l[i] + l[(i + 1) % 3];
            let s = if along > 0.0 {
                l[(i + 1) % 3] / along
            } else {
                0.0
            };
            let ribbon = &self.ribbons[i];
            let r = ribbon.at(s, d[i], tol)?;
            sum += r.to_vector() * (w[i] / total);
        }
        Ok(Point::from_vector(sum))
    }

    /// A four-sided patch at `(u, v)` of the unit square, its sides in turn
    /// along `v = 0`, `u = 1`, `v = 1` and `u = 0`.
    ///
    /// # Errors
    ///
    /// As evaluating the sides' curves.
    pub(crate) fn point_square(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let d = [v, 1.0 - u, 1.0 - v, u];
        let s = [u, v, 1.0 - u, 1.0 - v];
        let weight = |k: usize| -> f64 {
            (0..4)
                .filter(|j| *j != k)
                .map(|j| d[j].abs().powi(3))
                .product::<f64>()
        };
        let w = [weight(0), weight(1), weight(2), weight(3)];
        let total: f64 = w.iter().sum();
        if total <= f64::MIN_POSITIVE {
            // A corner: the side starting there.
            let corner = match (u > 0.5, v > 0.5) {
                (false, false) => 0,
                (true, false) => 1,
                (true, true) => 2,
                (false, true) => 3,
            };
            return self.ribbons[corner].point(0.0, tol);
        }
        let mut sum = Vector::ZERO;
        for i in 0..4 {
            if w[i] == 0.0 {
                continue;
            }
            let r = self.ribbons[i].at(s[i], d[i], tol)?;
            sum += r.to_vector() * (w[i] / total);
        }
        Ok(Point::from_vector(sum))
    }
}

/// The barycentrics of `(u, v)` of the unit square, its top row the third
/// corner.
fn barycentrics(u: f64, v: f64) -> [f64; 3] {
    [(1.0 - u) * (1.0 - v), u * (1.0 - v), v]
}

/// A vertex patch's face and what its interpolated surface achieved.
#[derive(Debug, Clone)]
pub(crate) struct CornerPatch {
    /// The face, bounded by the sides' edges, and a pole where there are
    /// three.
    pub face: Shape,
    /// The widest gap between the surface and a side's curve.
    pub gap: f64,
    /// The widest angle between the surface and a side's support.
    pub angle: f64,
    /// The farthest the surface stands from the exact patch at the middles
    /// of its interpolation cells.
    pub off: f64,
}

/// The vertex patch over the hole `sides` (three or four, walked corner to
/// corner, wound counter-clockwise seen from outside the solid), its face
/// built on the sides' edges, with the third corner as a pole where there
/// are three, and what it achieved.
///
/// # Errors
///
/// As [`Transfinite::new`], and as building the face.
pub(crate) fn corner_patch(
    model: &mut Model,
    sides: &[HoleSide],
    tol: Tolerances,
) -> OgeomResult<CornerPatch> {
    let exact = Transfinite::new(sides, tol)?;
    let surface = interpolated(&exact, tol)?;
    let (gap, angle) = measured(&exact, &surface, tol)?;
    let mut off = 0.0f64;
    #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
    let step = 1.0 / (NET - 1) as f64;
    for i in 0..NET - 1 {
        for j in 0..NET - 1 {
            #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
            let (u, v) = ((i as f64 + 0.5) * step, (j as f64 + 0.5) * step);
            let p = surface.point_at(u, v, tol)?;
            off = off.max(p.distance(Across { exact: &exact, tol }.point(u, v)?));
        }
    }
    let face = face_on(model, sides, surface, gap, tol)?;
    Ok(CornerPatch {
        face,
        gap,
        angle,
        off,
    })
}

/// The parameters of the interpolation points along each direction:
/// spaced by the cosine of the cosine, closest at the ends. Where a
/// corner's two supports curve differently, their twists there disagree in
/// the normal, which no choice of cross tangents mends; the exact patch
/// takes both, a polynomial only one, and it misses the other's tangent
/// plane over about its first span. The spans there fall as the square of
/// the cosine's.
fn net_parameters() -> Vec<f64> {
    #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
    let last = (NET - 1) as f64;
    (0..NET)
        .map(|k| {
            #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
            let t = k as f64 / last;
            if k == 0 {
                0.0
            } else if k == NET - 1 {
                1.0
            } else {
                {
                    let once = 0.5 * (1.0 - (core::f64::consts::PI * t).cos());
                    0.5 * (1.0 - (core::f64::consts::PI * once).cos())
                }
            }
        })
        .collect()
}

/// The exact patch over the unit square, and its derivatives across the
/// square's sides, from its ribbons: `S = R_i + O(d_i^2)` near side `i`, so
/// the first derivatives there are the ribbon's. Four sides: along `v = 0`,
/// `dS/dv = T_0(u)`; along `u = 1`, `dS/du = -T_1(v)`; along `v = 1`,
/// `dS/dv = -T_2(1 - u)`; along `u = 0`, `dS/du = T_3(1 - v)`. Three sides,
/// with `l = ((1 - u)(1 - v), u (1 - v), v)`: along
/// `v = 0`, `dS/dv = T_0(u)`; along `u = 1`,
/// `dS/du = -v (1 - v) C_1'(v) - (1 - v) T_1(v)`; along `u = 0`,
/// `dS/du = -v (1 - v) C_2'(1 - v) + (1 - v) T_2(1 - v)`; and at the pole
/// `v = 1`, where the patch is to first order the plane of the two sides
/// meeting there, `dS/dv = u C_1'(1) - (1 - u) C_2'(0)`.
struct Across<'a> {
    exact: &'a Transfinite,
    tol: Tolerances,
}

impl Across<'_> {
    fn square(&self) -> bool {
        self.exact.ribbons.len() == 4
    }

    fn point(&self, u: f64, v: f64) -> OgeomResult<Point> {
        if self.square() {
            self.exact.point_square(u, v, self.tol)
        } else {
            self.exact.point(barycentrics(u, v), self.tol)
        }
    }

    fn dv_bottom(&self, u: f64) -> OgeomResult<Vector> {
        self.exact.ribbons[0].cross(u, self.tol)
    }

    fn dv_top(&self, u: f64) -> OgeomResult<Vector> {
        let r = &self.exact.ribbons;
        if self.square() {
            return Ok(-r[2].cross(1.0 - u, self.tol)?);
        }
        Ok(r[1].tangent(1.0, self.tol)? * u - r[2].tangent(0.0, self.tol)? * (1.0 - u))
    }

    fn du_right(&self, v: f64) -> OgeomResult<Vector> {
        let r = &self.exact.ribbons[1];
        if self.square() {
            return Ok(-r.cross(v, self.tol)?);
        }
        Ok(r.tangent(v, self.tol)? * (-v * (1.0 - v)) - r.cross(v, self.tol)? * (1.0 - v))
    }

    fn du_left(&self, v: f64) -> OgeomResult<Vector> {
        if self.square() {
            return self.exact.ribbons[3].cross(1.0 - v, self.tol);
        }
        let r = &self.exact.ribbons[2];
        Ok(r.tangent(1.0 - v, self.tol)? * (-v * (1.0 - v))
            + r.cross(1.0 - v, self.tol)? * (1.0 - v))
    }

    /// The twist at a corner of the square: where the two sides' data
    /// agree (a pole's corners) their common value, and at the hole's own
    /// corners, where they need not, the mean of the two.
    fn twist(&self, u: f64, v: f64) -> OgeomResult<Vector> {
        let h = 1e-6;
        let along_u = |f: &dyn Fn(f64) -> OgeomResult<Vector>, at: f64| -> OgeomResult<Vector> {
            let (a, b) = if at <= 0.0 { (0.0, h) } else { (1.0 - h, 1.0) };
            Ok((f(b)? - f(a)?) / h)
        };
        let from_v = if v <= 0.0 {
            along_u(&|x| self.dv_bottom(x), u)?
        } else {
            along_u(&|x| self.dv_top(x), u)?
        };
        let from_u = if u <= 0.0 {
            along_u(&|x| self.du_left(x), v)?
        } else {
            along_u(&|x| self.du_right(x), v)?
        };
        Ok((from_v + from_u) * 0.5)
    }
}

/// The clamped bicubic B-spline taking the exact patch at [`NET`] by
/// [`NET`] points of the unit square and its derivatives across the
/// square's sides at each of them, its twists at the corners as
/// [`Across::twist`] gives them.
fn interpolated(exact: &Transfinite, tol: Tolerances) -> OgeomResult<BSplineSurface> {
    let params = net_parameters();
    let mut knots = vec![0.0; 4];
    knots.extend_from_slice(&params[1..NET - 1]);
    knots.extend_from_slice(&[1.0; 4]);
    let knots = KnotVector::new(knots, 3)?;
    let count = NET + 2;
    // Rows: the point at the first parameter, the derivative there, the
    // points between, the derivative at the last, the point there.
    let mut a = vec![vec![0.0; count]; count];
    let mut rows: Vec<(f64, usize)> = Vec::with_capacity(count);
    rows.push((0.0, 0));
    rows.push((0.0, 1));
    for &t in &params[1..NET - 1] {
        rows.push((t, 0));
    }
    rows.push((1.0, 1));
    rows.push((1.0, 0));
    for (row, &(t, order)) in a.iter_mut().zip(&rows) {
        let span = knots.span(t, tol)?;
        let values = knots.basis_derivatives(span, t, 1);
        for k in 0..=3 {
            row[span - 3 + k] = values[order][k];
        }
    }
    let solve = Collocation::factor(a)?;
    let across = Across { exact, tol };
    // The right-hand side along one line: its points at `params` and the
    // derivatives at its two ends, in the rows' order.
    let line = |points: Vec<Point>, start: Vector, end: Vector| -> Vec<Point> {
        let mut rhs = Vec::with_capacity(count);
        rhs.push(points[0]);
        rhs.push(Point::from_vector(start));
        rhs.extend_from_slice(&points[1..NET - 1]);
        rhs.push(Point::from_vector(end));
        rhs.push(points[NET - 1]);
        solve.apply(&rhs)
    };
    // Along `u` at each sampled `v`, the derivatives across `u = 0` and
    // `u = 1` at its ends.
    let mut along_u: Vec<Vec<Point>> = Vec::with_capacity(count);
    let mut rows_u: Vec<Vec<Point>> = Vec::with_capacity(NET);
    for &v in &params {
        let points = params
            .iter()
            .map(|&u| across.point(u, v))
            .collect::<OgeomResult<Vec<_>>>()?;
        rows_u.push(line(points, across.du_left(v)?, across.du_right(v)?));
    }
    // The derivatives across `v = 0` and `v = 1`, carried along `u` the
    // same way, the twists at their ends.
    let bottom = line(
        params
            .iter()
            .map(|&u| across.dv_bottom(u).map(Point::from_vector))
            .collect::<OgeomResult<Vec<_>>>()?,
        across.twist(0.0, 0.0)?,
        across.twist(1.0, 0.0)?,
    );
    let top = line(
        params
            .iter()
            .map(|&u| across.dv_top(u).map(Point::from_vector))
            .collect::<OgeomResult<Vec<_>>>()?,
        across.twist(0.0, 1.0)?,
        across.twist(1.0, 1.0)?,
    );
    along_u.extend(rows_u);
    // Then along `v` for each control column.
    let mut control = vec![Point::ORIGIN; count * count];
    for i in 0..count {
        let column: Vec<Point> = along_u.iter().map(|row| row[i]).collect();
        let solved = line(column, bottom[i].to_vector(), top[i].to_vector());
        for (j, p) in solved.into_iter().enumerate() {
            control[i * count + j] = p;
        }
    }
    BSplineSurface::new(
        knots.clone(),
        knots,
        &ControlGrid::new(control, count, count)?,
        tol,
    )
}

/// A square collocation matrix factored by Gaussian elimination with
/// partial pivoting.
struct Collocation {
    lu: Vec<Vec<f64>>,
    pivots: Vec<usize>,
}

impl Collocation {
    fn factor(mut a: Vec<Vec<f64>>) -> OgeomResult<Self> {
        let n = a.len();
        let mut pivots = Vec::with_capacity(n);
        for k in 0..n {
            let p = (k..n)
                .max_by(|x, y| a[*x][k].abs().total_cmp(&a[*y][k].abs()))
                .unwrap_or(k);
            if a[p][k].abs() <= f64::EPSILON {
                ogeom_bail!(Construction, "a vertex patch's interpolation is singular");
            }
            a.swap(k, p);
            pivots.push(p);
            let (upper, lower) = a.split_at_mut(k + 1);
            let pivot = &upper[k];
            for row in lower.iter_mut() {
                let factor = row[k] / pivot[k];
                row[k] = factor;
                for (x, p) in row[k + 1..].iter_mut().zip(&pivot[k + 1..]) {
                    *x -= factor * p;
                }
            }
        }
        Ok(Self { lu: a, pivots })
    }

    /// The solution for the right-hand side `b`, its points read as
    /// vectors.
    fn apply(&self, b: &[Point]) -> Vec<Point> {
        let n = self.lu.len();
        let mut x: Vec<Vector> = b.iter().map(|p| p.to_vector()).collect();
        for (k, &p) in self.pivots.iter().enumerate() {
            x.swap(k, p);
        }
        for r in 0..n {
            for c in 0..r {
                let sub = x[c] * self.lu[r][c];
                x[r] -= sub;
            }
        }
        for r in (0..n).rev() {
            for c in r + 1..n {
                let sub = x[c] * self.lu[r][c];
                x[r] -= sub;
            }
            x[r] /= self.lu[r][r];
        }
        x.into_iter().map(Point::from_vector).collect()
    }
}

/// The widest gap between `surface`'s sides and the exact sides' curves,
/// and the widest angle between its normal there and the supports', read
/// at eight points per interpolation span along each side, the angle
/// read short of the corners.
fn measured(
    exact: &Transfinite,
    surface: &BSplineSurface,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    let count = 8 * (NET - 1);
    let (mut gap, mut angle) = (0.0f64, 0.0f64);
    for (i, ribbon) in exact.ribbons.iter().enumerate() {
        for k in 0..=count {
            #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
            let s = k as f64 / count as f64;
            let (u, v) = match (exact.ribbons.len(), i) {
                (_, 0) => (s, 0.0),
                (_, 1) => (1.0, s),
                (4, 2) => (1.0 - s, 1.0),
                _ => (0.0, 1.0 - s),
            };
            if exact.ribbons.len() == 3 && v >= 1.0 {
                continue;
            }
            let p = surface.point_at(u, v, tol)?;
            let on = ribbon.point(s, tol)?;
            gap = gap.max(p.distance(on));
            // At the hole's corners the surface's two tangents are the
            // sides' own, taken exactly from the derivatives it
            // interpolates, both in the face both strips touch there: its
            // tangent plane is that face's, and where the hole runs
            // straight on the chart has no normal to read. The samples
            // beside each corner read the surface's.
            if k == 0 || k == count {
                continue;
            }
            let (du, dv) = surface.d1_at(u, v, tol)?;
            let n = du.cross(dv);

            let support = support_normal(&ribbon.side.support, on, tol)?;
            let sine = n.cross(support).magnitude() / n.magnitude();
            angle = angle.max(sine.clamp(0.0, 1.0).asin());
        }
    }
    Ok((gap, angle))
}

/// The face on `surface` bounded by the sides' edges, and a pole at the
/// third corner where there are three sides, each edge's pcurve the
/// straight side of the unit square it maps onto, the edges widened to the
/// gap measured.
fn face_on(
    model: &mut Model,
    sides: &[HoleSide],
    surface: BSplineSurface,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let id = model.geometry_mut().add_surface(surface.into());
    // Each side's ends in the unit square, in the walk's order.
    let (o, x, y, xy) = (
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(1.0, 1.0),
    );
    let corners: &[(Point2, Point2)] = if sides.len() == 4 {
        &[(o, x), (x, xy), (xy, y), (y, o)]
    } else {
        &[(o, x), (x, xy), (y, o)]
    };
    let widened = Tolerance::new(gap + tol.confusion())?;
    let mut walked = Vec::with_capacity(4);
    for (side, &(from, to)) in sides.iter().zip(corners) {
        let (a, b) = side.range;
        let (start, end) = if side.forward { (from, to) } else { (to, from) };
        let line = BSpline2d::new(KnotVector::new(vec![a, a, b, b], 1)?, vec![start, end], tol)?;
        ogeom_algo::attach_pcurve(
            model,
            &side.edge,
            line.into(),
            id,
            Location::identity(),
            (a, b),
        )?;
        model.widen(&side.edge, widened)?;
        if let Some((p, q)) = ogeom_algo::edge_vertices(model, &side.edge)? {
            model.widen(&p, widened)?;
            model.widen(&q, widened)?;
        }
        walked.push(if side.forward {
            side.edge.clone()
        } else {
            side.edge.reversed()
        });
    }
    if sides.len() == 4 {
        let wire = ogeom_algo::make_wire(model, &walked, tol)?.shape;
        return Ok(ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape);
    }
    // The pole: the walk's third corner, where the second side ends.
    let Some((_, pole)) = ogeom_algo::edge_vertices(model, &walked[1])? else {
        ogeom_bail!(Construction, "a vertex patch's side has no vertices");
    };
    let collapsed = {
        let mut data = EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[pole.clone(), pole])?
    };
    let top = BSpline2d::new(
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?,
        vec![Point2::new(1.0, 1.0), Point2::new(0.0, 1.0)],
        tol,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &collapsed,
        top.into(),
        id,
        Location::identity(),
        (0.0, 1.0),
    )?;
    walked.insert(2, collapsed);
    let wire = ogeom_algo::make_wire(model, &walked, tol)?.shape;
    Ok(ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;
    use crate::round::section::{PlanarSection, RoundSpec};
    use crate::support::Seat;
    use ogeom_algo::HalfSpaces;

    const T: Tolerances = Tolerances::millimetres();
    const FIT: f64 = 2e-4;

    /// One strip of a hole: its seat, radius, the face of its seat each of
    /// its ends lies on in the walk's order, and each face's station.
    struct Strip {
        seat: Seat,
        radius: f64,
        walk_from: usize,
        stations: [f64; 2],
    }

    fn seat(m: &mut Model, end: Point, normals: [Vector; 2], convex: bool) -> Seat {
        let face = ogeom_algo::make_vertex(m, Point::ORIGIN).shape;
        Seat {
            start: Point::ORIGIN,
            end,
            along: (end - Point::ORIGIN) / end.distance(Point::ORIGIN),
            normals,
            faces: [face.clone(), face],
            convex,
        }
    }

    /// The hole bounded by `strips`' skew ends, walked through `corners`.
    fn hole(m: &mut Model, corners: [Point; 3], strips: [Strip; 3]) -> [HoleSide; 3] {
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| ogeom_algo::make_vertex(m, *p).shape)
            .collect();
        let sides: Vec<HoleSide> = strips
            .iter()
            .enumerate()
            .map(|(i, strip)| {
                let section =
                    PlanarSection::new(&strip.seat, &RoundSpec::Circular(strip.radius), T).unwrap();
                let order = section.order();
                let at = [strip.stations[order[0]], strip.stations[order[1]]];
                let (curve, range) = section.skew_end(order, at, T).unwrap();
                let forward = order[0] == strip.walk_from;
                let (from, to) = (&vertices[i], &vertices[(i + 1) % 3]);
                let (a, b) = if forward { (from, to) } else { (to, from) };
                let edge = ogeom_algo::make_edge_between(m, curve.clone(), range, a, b, T)
                    .unwrap()
                    .shape;
                let support = section.surface((-20.0, 20.0), T).unwrap();
                HoleSide {
                    edge,
                    curve,
                    range,
                    forward,
                    support,
                }
            })
            .collect();
        [sides[0].clone(), sides[1].clone(), sides[2].clone()]
    }

    /// A convex box corner at the origin, its material the positive octant,
    /// its edges along `x`, `y` and `z` rounded at `r[0]`, `r[1]`, `r[2]`.
    /// Each face's setback corner stands each of its edges' radii in from
    /// them: `(r1, r0, 0)` on `z = 0`, `(r2, 0, r0)` on `y = 0`, `(0, r2,
    /// r1)` on `x = 0`.
    fn box_corner(m: &mut Model, r: [f64; 3]) -> [HoleSide; 3] {
        let (x, y, z) = (
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        );
        let corners = [
            Point::new(r[1], r[0], 0.0),
            Point::new(r[2], 0.0, r[0]),
            Point::new(0.0, r[2], r[1]),
        ];
        let strips = [
            // Along x, faces y = 0 and z = 0: from z = 0 to y = 0.
            Strip {
                seat: seat(m, Point::new(10.0, 0.0, 0.0), [-y, -z], true),
                radius: r[0],
                walk_from: 1,
                stations: [r[2], r[1]],
            },
            // Along z, faces y = 0 and x = 0: from y = 0 to x = 0.
            Strip {
                seat: seat(m, Point::new(0.0, 0.0, 10.0), [-y, -x], true),
                radius: r[2],
                walk_from: 0,
                stations: [r[0], r[1]],
            },
            // Along y, faces x = 0 and z = 0: from x = 0 to z = 0.
            Strip {
                seat: seat(m, Point::new(0.0, 10.0, 0.0), [-x, -z], true),
                radius: r[1],
                walk_from: 0,
                stations: [r[2], r[0]],
            },
        ];
        hole(m, corners, strips)
    }

    fn octant() -> HalfSpaces {
        let side = |n: Vector| HalfSpaces::Below {
            normal: n,
            at: Point::ORIGIN,
        };
        HalfSpaces::Every(vec![
            side(Vector::new(-1.0, 0.0, 0.0)),
            side(Vector::new(0.0, -1.0, 0.0)),
            side(Vector::new(0.0, 0.0, -1.0)),
        ])
    }

    /// The exact patch takes each side's curve and is tangent to its
    /// support all along it, at a corner of three different radii whose
    /// supports' twists disagree at every corner: on each side, to
    /// rounding; its normal there, read by central differences, to a
    /// millionth of a radian.
    #[test]
    fn the_exact_patch_takes_each_side_and_its_support_s_tangent_plane() {
        let mut m = Model::new();
        let sides = box_corner(&mut m, [1.0, 1.5, 2.0]);
        let patch = Transfinite::new(&sides, T).unwrap();
        for (i, side) in sides.iter().enumerate() {
            let ribbon = &patch.ribbons[i];
            for k in 1..20 {
                let s = f64::from(k) / 20.0;
                let mut l = [0.0; 3];
                l[i] = 1.0 - s;
                l[(i + 1) % 3] = s;
                let on = patch.point(l, T).unwrap();
                assert!(on.distance(ribbon.point(s, T).unwrap()) < 1e-12);
                // Two directions at the side, a hundred-thousandth apart:
                // along it centrally, into the patch one-sided.
                let h = 1e-5;
                let shifted = |by: [f64; 3], k: f64| {
                    let q = [l[0] + by[0] * k, l[1] + by[1] * k, l[2] + by[2] * k];
                    patch.point(q, T).unwrap()
                };
                let mut along = [0.0; 3];
                along[i] = -1.0;
                along[(i + 1) % 3] = 1.0;
                let mut inward = [0.0; 3];
                inward[i] = -0.5;
                inward[(i + 1) % 3] = -0.5;
                inward[(i + 2) % 3] = 1.0;
                let a = shifted(along, h) - shifted(along, -h);
                // Inward one-sided, to second order: the curves end at the
                // side.
                let b = (shifted(inward, h) - on) * 4.0 - (shifted(inward, 2.0 * h) - on);
                let n = a.cross(b);
                let support = support_normal(&side.support, on, T).unwrap();
                let angle = (n.cross(support).magnitude() / n.magnitude()).asin();
                assert!(angle < 1e-6, "side {i} at {s}: {angle}");
            }
        }
    }

    /// Three straight sides on the plane `z = 0`, each supported by that
    /// plane: every ribbon lies in it, so the patch does, exactly.
    #[test]
    fn a_flat_hole_is_patched_flat() {
        let mut m = Model::new();
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(1.0, 3.0, 0.0),
        ];
        let plane: SurfaceGeometry =
            ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(ogeom_math::Frame::WORLD)).into();
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| ogeom_algo::make_vertex(&mut m, *p).shape)
            .collect();
        let sides: Vec<HoleSide> = (0..3)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 3]);
                let line = ogeom_geom::LineCurve::segment(a, b, T).unwrap();
                let range = (0.0, a.distance(b));
                let edge = ogeom_algo::make_edge_between(
                    &mut m,
                    Curve::Line(line),
                    range,
                    &vertices[i],
                    &vertices[(i + 1) % 3],
                    T,
                )
                .unwrap()
                .shape;
                HoleSide {
                    edge,
                    curve: Curve::Line(line),
                    range,
                    forward: true,
                    support: plane.clone(),
                }
            })
            .collect();
        let patch =
            Transfinite::new(&[sides[0].clone(), sides[1].clone(), sides[2].clone()], T).unwrap();
        for i in 0..=10 {
            for j in 0..=(10 - i) {
                let l = [
                    f64::from(i) / 10.0,
                    f64::from(j) / 10.0,
                    f64::from(10 - i - j) / 10.0,
                ];
                let p = patch.point(l, T).unwrap();
                assert!(p.z.abs() < 1e-12, "{l:?}: {p:?}");
            }
        }
    }

    /// Four straight sides of a quadrilateral on `z = 0`, each supported by
    /// that plane: the four-sided patch lies in it, exactly, and takes each
    /// side.
    #[test]
    fn a_flat_four_sided_hole_is_patched_flat() {
        let mut m = Model::new();
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(5.0, 0.0, 0.0),
            Point::new(4.0, 3.0, 0.0),
            Point::new(1.0, 4.0, 0.0),
        ];
        let plane: SurfaceGeometry =
            ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(ogeom_math::Frame::WORLD)).into();
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| ogeom_algo::make_vertex(&mut m, *p).shape)
            .collect();
        let sides: Vec<HoleSide> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                let line = ogeom_geom::LineCurve::segment(a, b, T).unwrap();
                let range = (0.0, a.distance(b));
                let edge = ogeom_algo::make_edge_between(
                    &mut m,
                    Curve::Line(line),
                    range,
                    &vertices[i],
                    &vertices[(i + 1) % 4],
                    T,
                )
                .unwrap()
                .shape;
                HoleSide {
                    edge,
                    curve: Curve::Line(line),
                    range,
                    forward: true,
                    support: plane.clone(),
                }
            })
            .collect();
        let patch = Transfinite::new(&sides, T).unwrap();
        for i in 0..=10 {
            for j in 0..=10 {
                let (u, v) = (f64::from(i) / 10.0, f64::from(j) / 10.0);
                let p = patch.point_square(u, v, T).unwrap();
                assert!(p.z.abs() < 1e-12, "{u} {v}: {p:?}");
            }
            let s = f64::from(i) / 10.0;
            let on = patch.point_square(s, 0.0, T).unwrap();
            assert!(on.distance(corners[0].lerp(corners[1], s)) < 1e-12);
        }
        let built = corner_patch(&mut m, &sides, T).unwrap();
        assert!(built.gap <= FIT && built.angle <= FIT, "{built:?}");
    }

    /// At a convex box corner of three different radii, the corner patch's
    /// interpolated surface meets each strip's end within the bands' fit
    /// target, in position and tangent plane, and is proved inside the
    /// corner's three planes.
    #[test]
    fn a_corner_of_three_radii_patches_within_the_fit_target() {
        let mut m = Model::new();
        let sides = box_corner(&mut m, [1.0, 1.5, 2.0]);
        let patch = corner_patch(&mut m, &sides, T).unwrap();
        assert!(patch.gap <= FIT, "gap {}", patch.gap);
        assert!(patch.angle <= FIT, "angle {}", patch.angle);
        let proved = ogeom_algo::certify_face_within(&m, &patch.face, &octant(), FIT, T).unwrap();
        assert_eq!(proved, Ok(()), "off the exact patch {}", patch.off);
    }

    /// A corner rounded at 1, 3 and 2 along `x`, `y` and `z`: the strips'
    /// curvatures differ threefold, and a patch of straight ribbons stands
    /// past a face beside a corner of the hole, where a strip tangent to
    /// that face curves away from the straight line. With ribbons bending as
    /// their strips do, the patch meets its strips within the fit target and
    /// is proved inside the corner's three planes.
    #[test]
    fn a_corner_of_radii_one_three_and_two_stays_inside_its_planes() {
        let mut m = Model::new();
        let sides = box_corner(&mut m, [1.0, 3.0, 2.0]);
        let patch = corner_patch(&mut m, &sides, T).unwrap();
        assert!(patch.gap <= FIT, "gap {}", patch.gap);
        assert!(patch.angle <= FIT, "angle {}", patch.angle);
        let proved = ogeom_algo::certify_face_within(&m, &patch.face, &octant(), FIT, T).unwrap();
        assert_eq!(proved, Ok(()), "off the exact patch {}", patch.off);
    }

    /// The inner vertex of an L-shaped prism, its cap `z = 0` (material
    /// above), its two sides `x = 0` (material `x >= 0` beyond `y < 0`) and
    /// `y = 0`, the edge between the sides concave and the cap's two edges
    /// convex, all rounded at 1: the convex strips' contacts on the cap run
    /// on past the vertex to the cap's setback corner `(1, 1, 0)`, their ends
    /// skew, the concave strip's end square. The patch meets every strip
    /// within the fit target and is proved inside the corner's material,
    /// the concave round's fill added.
    #[test]
    fn a_mixed_corner_patches_within_the_fit_target() {
        let mut m = Model::new();
        let r = 1.0;
        let (x, y, z) = (
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        );
        let corners = [
            Point::new(0.0, -r, r),
            Point::new(-r, 0.0, r),
            Point::new(r, r, 0.0),
        ];
        let strips = [
            // The concave edge along z, faces x = 0 and y = 0.
            Strip {
                seat: seat(&mut m, Point::new(0.0, 0.0, 10.0), [-x, -y], false),
                radius: r,
                walk_from: 0,
                stations: [r, r],
            },
            // The cap's edge along -x, faces y = 0 and z = 0.
            Strip {
                seat: seat(&mut m, Point::new(-10.0, 0.0, 0.0), [-y, -z], true),
                radius: r,
                walk_from: 0,
                stations: [r, -r],
            },
            // The cap's edge along -y, faces x = 0 and z = 0.
            Strip {
                seat: seat(&mut m, Point::new(0.0, -10.0, 0.0), [-x, -z], true),
                radius: r,
                walk_from: 1,
                stations: [r, -r],
            },
        ];
        let sides = hole(&mut m, corners, strips);
        let patch = corner_patch(&mut m, &sides, T).unwrap();
        assert!(patch.gap <= FIT, "gap {}", patch.gap);
        assert!(patch.angle <= FIT, "angle {}", patch.angle);
        let below = |n: Vector, rise: f64| HalfSpaces::Below {
            normal: n,
            at: Point::ORIGIN + n * rise,
        };
        let region = HalfSpaces::Every(vec![
            below(-z, 0.0),
            HalfSpaces::Any(vec![
                below(-x, 0.0),
                below(-y, 0.0),
                HalfSpaces::Every(vec![below(-x, r), below(-y, r)]),
            ]),
        ]);
        let proved = ogeom_algo::certify_face_within(&m, &patch.face, &region, FIT, T).unwrap();
        assert_eq!(proved, Ok(()), "off the exact patch {}", patch.off);
    }
}
