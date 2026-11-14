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
//! cubed, it enters at third. Five or more sides blend over the regular
//! polygon with as many corners, `d_i = 1 - l_i - l_{i+1}` and
//! `s_i = l_{i+1} / (l_i + l_{i+1})` from its Wachspress coordinates `l`
//! (for a triangle and a square these are the formulas above), their
//! distances cubed too.
//!
//! The patch is not a tensor product, so the face carries it interpolated:
//! a clamped bicubic B-spline taking it at [`NET`] by [`NET`] points of the
//! unit square, and its derivatives across the square's sides at each,
//! three sides' square mapped onto the triangle with its third corner as a
//! pole, four sides' square taken as it is. Five or more sides' polygon is
//! cut from its centre into a triangle on each side, each carried on its own
//! face as three sides' are, its pole at the centre; neighbours take the
//! line between them from the same points and derivatives along it, so
//! they share it exactly, and the derivatives across it and at the centre
//! come from the one exact patch. What the interpolated surface
//! achieves is measured on it, never assumed from the exact patch: its gap
//! to each side's curve, its angle to each side's support, and how far it
//! stands from the exact patch inside.
//!
//! A side may run over several edges, each on its own support (an S side:
//! a convex round's end running on tangent into a concave one's, the two
//! bending opposite ways); its normal's sign is carried across each join,
//! and its bend eased across it. Four sides whose second and fourth sides
//! each make such a join at one `v` are fitted as two faces split along
//! that line, which one smooth surface cannot follow, on a finer net.
//! Four sides' weights read each distance in the hole's extent across it
//! (the mean length of the two sides it runs between, to the fourth), so
//! a long hole's short sides reach only near their own ends: their ribbons,
//! carried half the hole's length, are not the surface there. Each cross
//! tangent's part into the hole is lifted at its side's middle to the
//! distance across to the opposite side, so the two meet rather than
//! overshoot, and a corner's slope reaches along its side only as far as
//! lowers that part by a quarter at most, so the patch never folds over
//! its side; on a square-like hole both leave the patch as it was.

use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{BSpline2d, BSplineSurface, Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{ControlGrid, KnotVector, Point, Point2, Vector};
use ogeom_topo::{EdgeData, Location, Model, Shape};

/// The interpolation points along each direction of the unit square.
pub(crate) const NET: usize = 33;

/// One side of a hole: an edge's curve walked from one corner
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
    /// Where the side runs over more than one edge (an S side, a strip's
    /// end running on into another's through a point both are tangent at):
    /// each edge in the walk's order with its own surface, and `curve` the
    /// whole side's; empty for a side of one edge, `edge` and `support`.
    pub pieces: Vec<SidePiece>,
}

/// One edge of a side running over several.
#[derive(Debug, Clone)]
pub(crate) struct SidePiece {
    /// The edge, as a node (its own orientation forward).
    pub edge: Shape,
    /// The range of its own curve it covers.
    pub range: (f64, f64),
    /// Whether the walk runs with its own curve's parameter.
    pub forward: bool,
    /// Where along the side it starts and ends, between nought and one.
    pub along: (f64, f64),
    /// The surface across it.
    pub support: SurfaceGeometry,
}

/// An edge of a side: the edge, its own range, whether the walk runs with
/// it, and where along the side it starts and ends.
type SideEdge = (Shape, (f64, f64), bool, (f64, f64));

/// How far either side of a join of two pieces on different surfaces the
/// ribbon's bend is carried from the one's to the other's, along the side.
const JOIN_REACH: f64 = 0.05;

impl HoleSide {
    /// The surface across the side at `s` along it.
    fn support_at(&self, s: f64) -> &SurfaceGeometry {
        self.pieces
            .iter()
            .find(|p| s <= p.along.1)
            .or(self.pieces.last())
            .map_or(&self.support, |p| &p.support)
    }

    /// Whether the side runs from one kind of surface onto another.
    fn bends(&self) -> bool {
        self.pieces.windows(2).any(|w| {
            core::mem::discriminant(&w[0].support) != core::mem::discriminant(&w[1].support)
        })
    }

    /// The side's edges in the walk's order: each with its own range,
    /// whether the walk runs with it, and where along the side it stands.
    fn edges(&self) -> Vec<SideEdge> {
        if self.pieces.is_empty() {
            return vec![(self.edge.clone(), self.range, self.forward, (0.0, 1.0))];
        }
        self.pieces
            .iter()
            .map(|p| (p.edge.clone(), p.range, p.forward, p.along))
            .collect()
    }
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
    /// How far along the side from each corner its slope reaches: the whole
    /// side, or, where the slopes carried over it turn the cross tangent's
    /// part into the hole back out between the corners (a long side between
    /// sharply turning corners, the patch folding over its side), as far as
    /// lowers that part by a quarter at most.
    reach: [f64; 2],
    /// For a side of four, how far the part into the hole is lowered at the
    /// side's middle, by a bump nought with its slope at both corners, where
    /// it reaches past the opposite side there: to that distance, so the two
    /// sides' ribbons meet across a long hole rather than overshoot (a short
    /// side's speed at the corners, carried along, can be far more).
    lift: f64,
    /// Each piece's normal's sign, where the side runs over several: each
    /// turned to agree with the last's at their join (a convex round's
    /// normal faces out of the solid, a concave one's torus's into it).
    turns: Vec<f64>,
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
        let n = support_normal(self.side.support_at(s), p, tol)? * self.turn(s);
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
        // A slope's basis `x (1 - x)^2` over its reach from its corner: its
        // value nought and its slope one there, both nought where it ends.
        let slope = |x: f64, reach: f64| {
            let x = x / reach;
            if x < 1.0 {
                reach * x * (1.0 - x) * (1.0 - x)
            } else {
                0.0
            }
        };
        let h = [
            2.0 * s3 - 3.0 * s2 + 1.0,
            slope(s, self.reach[0]),
            -2.0 * s3 + 3.0 * s2,
            -slope(1.0 - s, self.reach[1]),
        ];
        let along = h[0] * a0 + h[1] * da0 + h[2] * a1 + h[3] * da1;
        let bump = 16.0 * s2 * (1.0 - s) * (1.0 - s);
        let into = h[0] * b0 + h[1] * db0 + h[2] * b1 + h[3] * db1 + bump * self.lift;
        Ok(t * along + across * into)
    }

    /// The ribbon `d` across from the side at `s`: off the curve along the
    /// cross tangent, bending as the support does that way, so that near
    /// its side it keeps to the support's own continuation to second
    /// order.
    fn at(&self, s: f64, d: f64, tol: Tolerances) -> OgeomResult<Point> {
        let p = self.point(s, tol)?;
        let t = self.cross(s, tol)?;
        let bend = self.bend(s, p, t, tol)?;
        Ok(p + t * d + bend * (0.5 * d * d))
    }

    /// The sign [`Ribbon::turns`] gives the piece at `s`.
    fn turn(&self, s: f64) -> f64 {
        self.side
            .pieces
            .iter()
            .position(|p| s <= p.along.1)
            .or_else(|| self.side.pieces.len().checked_sub(1))
            .and_then(|k| self.turns.get(k).copied())
            .unwrap_or(1.0)
    }

    /// How the ribbon bends at `s`, from `p` along `t`: as the support does,
    /// and within [`JOIN_REACH`] of a join of two pieces on different
    /// surfaces carried evenly from the one's bend to the other's, so the
    /// ribbon holds together across it.
    fn bend(&self, s: f64, p: Point, t: Vector, tol: Tolerances) -> OgeomResult<Vector> {
        for pair in self.side.pieces.windows(2) {
            let at = pair[0].along.1;
            if (s - at).abs() < JOIN_REACH {
                let x = (s - at + JOIN_REACH) / (2.0 * JOIN_REACH);
                let w = x * x * (3.0 - 2.0 * x);
                let (a, b) = (
                    support_bend(&pair[0].support, p, t, tol)?,
                    support_bend(&pair[1].support, p, t, tol)?,
                );
                return Ok(a * (1.0 - w) + b * w);
            }
        }
        support_bend(self.side.support_at(s), p, t, tol)
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
        SurfaceGeometry::Torus(to) => {
            // Off the tube's core: the meridian's curvature is the tube's,
            // the parallel's `cos phi / (R + a cos phi)`, `phi` the angle
            // round the tube from the outer equator.
            let torus = to.torus();
            let frame = torus.frame();
            let (o, z) = (frame.origin(), frame.z().vector());
            let (big, small) = (torus.major_radius(), torus.minor_radius());
            let off = p - o;
            let radial = (off - z * off.dot(z)).normalized(tol)?;
            let core = o + radial * big;
            let out = (p - core).normalized(tol)?;
            let cos_phi = out.dot(radial);
            let parallel = z.cross(radial);
            let meridian = out.cross(parallel);
            let (tm, tp) = (t.dot(meridian), t.dot(parallel));
            let bend = tm * tm / small + tp * tp * cos_phi / small.mul_add(cos_phi, big);
            Ok(out * -bend)
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
        SurfaceGeometry::Cone(c) => {
            // Off the axis, leaning back against the cone's opening.
            let cone = c.cone();
            let frame = cone.frame();
            let off = p - frame.origin();
            let z = frame.z().vector();
            let radial = (off - z * off.dot(z)).normalized(tol)?;
            (radial - z * cone.half_angle().tan()).normalized(tol)
        }
        SurfaceGeometry::Torus(to) => {
            let torus = to.torus();
            let frame = torus.frame();
            let (o, z) = (frame.origin(), frame.z().vector());
            let off = p - o;
            let radial = (off - z * off.dot(z)).normalized(tol)?;
            (p - (o + radial * torus.major_radius())).normalized(tol)
        }
        _ => ogeom_bail!(
            Construction,
            "a vertex patch's side runs along a surface whose normal it cannot read"
        ),
    }
}

/// The exact patch over a hole of three sides or more: Kato's blend of its
/// sides' ribbons, or, over four sides or more, the corners blended by the
/// polygon's Wachspress coordinates ([`corner_patch_through_corners`]).
#[derive(Debug, Clone)]
pub(crate) struct Transfinite {
    ribbons: Vec<Ribbon>,
    /// The hole's corners, where the patch is their Wachspress blend.
    through: Option<Vec<Point>>,
    /// For four sides, the hole's extent along the square's `u` and `v`:
    /// the mean length of the first and third sides, and of the second and
    /// fourth.
    extent: [f64; 2],
}

impl Transfinite {
    /// The patch over the hole `sides`, walked corner to corner, each side
    /// ending where the next starts.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where there are fewer than three sides, the sides do not close, or a
    /// support's normal cannot be read.
    pub(crate) fn new(sides: &[HoleSide], tol: Tolerances) -> OgeomResult<Self> {
        let n = sides.len();
        if n < 3 {
            ogeom_bail!(
                Construction,
                "a vertex patch takes three sides or more, not {n}"
            );
        }
        let mut ribbons: Vec<Ribbon> = sides
            .iter()
            .map(|side| Ribbon {
                side: side.clone(),
                sign: 1.0,
                ends: [(0.0, 0.0); 2],
                slopes: [(0.0, 0.0); 2],
                reach: [1.0; 2],
                lift: 0.0,
                turns: Vec::new(),
            })
            .collect();
        for r in &mut ribbons {
            let mut turns = Vec::with_capacity(r.side.pieces.len());
            for (k, piece) in r.side.pieces.iter().enumerate() {
                if k == 0 {
                    turns.push(1.0);
                    continue;
                }
                let at = r.point(piece.along.0, tol)?;
                let before = support_normal(&r.side.pieces[k - 1].support, at, tol)?;
                let after = support_normal(&piece.support, at, tol)?;
                turns.push(turns[k - 1] * before.dot(after).signum());
            }
            r.turns = turns;
        }
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
        // Where the part into the hole, carried with its slopes over the
        // whole side, turns back out of the hole between the corners (the
        // patch folding over its side), each slope reaches only as far as
        // its basis, lowering the part by at most 4/27 of its reach times
        // the slope, takes a quarter of it off.
        for r in &mut ribbons {
            let (b0, b1) = (r.ends[0].1, r.ends[1].1);
            let mut folds = false;
            for i in 1..32 {
                let s = f64::from(i) / 32.0;
                let (s2, s3) = (s * s, s * s * s);
                let into = (2.0 * s3 - 3.0 * s2 + 1.0) * b0
                    + (s3 - 2.0 * s2 + s) * r.slopes[0].1
                    + (-2.0 * s3 + 3.0 * s2) * b1
                    + (s3 - s2) * r.slopes[1].1;
                if into * b0.signum() < 0.0 && b0 * b1 > 0.0 {
                    folds = true;
                }
            }
            if !folds {
                continue;
            }
            for k in 0..2 {
                let (into, rate) = (r.ends[k].1.abs(), r.slopes[k].1.abs());
                if rate > 0.0 {
                    r.reach[k] = (27.0 / 16.0 * into / rate).min(1.0);
                }
            }
        }
        if n == 4 {
            for i in 0..4 {
                let r = &ribbons[i];
                let (_, across) = r.frame(0.5, r.sign, tol)?;
                let opposite = ribbons[(i + 2) % 4].point(0.5, tol)?;
                let wanted = (opposite - r.point(0.5, tol)?).dot(across);
                let into = r.cross(0.5, tol)?.dot(across);
                if wanted > 0.0 && into > wanted {
                    ribbons[i].lift = wanted - into;
                }
            }
        }
        let mut extent = [1.0; 2];
        if n == 4 {
            let mut lengths = [0.0; 4];
            for (length, r) in lengths.iter_mut().zip(&ribbons) {
                let mut last = r.point(0.0, tol)?;
                for k in 1..=16 {
                    let p = r.point(f64::from(k) / 16.0, tol)?;
                    *length += p.distance(last);
                    last = p;
                }
            }
            extent = [
                0.5 * (lengths[0] + lengths[2]),
                0.5 * (lengths[1] + lengths[3]),
            ];
        }
        Ok(Self {
            ribbons,
            through: None,
            extent,
        })
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
        // The weights read each distance in the hole's extent across it, so
        // a long hole's short sides reach only near their own ends.
        let [across_u, across_v] = self.extent.map(|e| e.powi(4));
        let far = [
            v * across_v,
            (1.0 - u) * across_u,
            (1.0 - v) * across_v,
            u * across_u,
        ];
        let s = [u, v, 1.0 - u, 1.0 - v];
        let weight = |k: usize| -> f64 {
            (0..4)
                .filter(|j| *j != k)
                .map(|j| far[j].abs().powi(3))
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

impl Transfinite {
    /// The corners of the regular polygon a patch of five sides or more
    /// blends over, side `i` running from corner `i` to corner `i + 1`.
    fn polygon(&self) -> Vec<Point2> {
        let n = self.ribbons.len();
        (0..n)
            .map(|k| {
                #[expect(clippy::cast_precision_loss, reason = "a handful of sides")]
                let t = core::f64::consts::TAU * k as f64 / n as f64;
                Point2::new(t.cos(), t.sin())
            })
            .collect()
    }

    /// A patch of five sides or more at `p` of its polygon.
    ///
    /// # Errors
    ///
    /// As evaluating the sides' curves.
    fn point_polygon(&self, p: Point2, tol: Tolerances) -> OgeomResult<Point> {
        let corners = self.polygon();
        let n = corners.len();
        // Twice the area of the triangle `p` makes with side `j`: nought on
        // that side, positive inside.
        let area = |j: usize| {
            let (a, b) = (corners[j], corners[(j + 1) % n]);
            (a.x - p.x) * (b.y - p.y) - (a.y - p.y) * (b.x - p.x)
        };
        let areas: Vec<f64> = (0..n).map(|j| area(j).max(0.0)).collect();
        // Wachspress's coordinates, the corners' constant factors equal on
        // a regular polygon: each the product of the areas over the sides
        // not meeting at its corner.
        let mut l: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .filter(|j| *j != i && *j != (i + n - 1) % n)
                    .map(|j| areas[j])
                    .product::<f64>()
            })
            .collect();
        let sum: f64 = l.iter().sum();
        for x in &mut l {
            *x /= sum;
        }
        if let Some(corners) = &self.through {
            // The corners' blend, and each side's departure from its chord
            // weighted by its two corners' coordinates: on side `i` the
            // blend is the chord and the departure the curve's, while the
            // sides beside it stand at an end of theirs, where a curve meets
            // its chord.
            let mut blend = corners
                .iter()
                .zip(&l)
                .fold(Vector::ZERO, |acc, (c, w)| acc + c.to_vector() * *w);
            for i in 0..n {
                let j = (i + 1) % n;
                let along = l[i] + l[j];
                if along <= 0.0 {
                    continue;
                }
                let t = l[j] / along;
                let chord = corners[i].to_vector() * (1.0 - t) + corners[j].to_vector() * t;
                let curve = self.ribbons[i].point(t, tol)?.to_vector();
                blend += (curve - chord) * along;
            }
            return Ok(Point::from_vector(blend));
        }
        let d: Vec<f64> = (0..n)
            .map(|i| (1.0 - l[i] - l[(i + 1) % n]).max(0.0))
            .collect();
        let w: Vec<f64> = (0..n)
            .map(|k| {
                (0..n)
                    .filter(|j| *j != k)
                    .map(|j| d[j].powi(3))
                    .product::<f64>()
            })
            .collect();
        let total: f64 = w.iter().sum();
        if total <= f64::MIN_POSITIVE {
            // A corner: the side starting there.
            let corner = (0..n).fold(0, |best, k| if l[k] > l[best] { k } else { best });
            return self.ribbons[corner].point(0.0, tol);
        }
        let mut sum = Vector::ZERO;
        for i in 0..n {
            if w[i] == 0.0 {
                continue;
            }
            let along = l[i] + l[(i + 1) % n];
            let s = if along > 0.0 {
                l[(i + 1) % n] / along
            } else {
                0.0
            };
            let r = self.ribbons[i].at(s, d[i], tol)?;
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
    /// The faces: one bounded by the sides' edges, with a pole where there
    /// are three; for five sides or more, one on each side, bounded by it
    /// and the lines from its ends to the centre, its pole there.
    pub faces: Vec<Shape>,
    /// The widest gap between the surface and a side's curve, or between
    /// two of its faces along the line they share.
    pub gap: f64,
    /// The widest angle between the surface and a side's support, or
    /// between two of its faces along the line they share.
    pub angle: f64,
    /// The farthest the surface stands from the exact patch at the middles
    /// of its interpolation cells.
    pub off: f64,
}

/// The vertex patch over the hole `sides` (walked corner to corner, wound
/// counter-clockwise seen from outside the solid), its faces built on the
/// sides' edges (see [`CornerPatch::faces`]), and what it achieved.
///
/// # Errors
///
/// As [`Transfinite::new`], and as building the face.
pub(crate) fn corner_patch(
    model: &mut Model,
    sides: &[HoleSide],
    tol: Tolerances,
) -> OgeomResult<CornerPatch> {
    if sides.len() != 4 && sides.iter().any(|s| !s.pieces.is_empty()) {
        ogeom_bail!(
            Construction,
            "a vertex patch takes a side of several edges only among four"
        );
    }
    // Sides running from one surface onto another second and fourth,
    // where they are.
    let turned: Vec<HoleSide>;
    let sides = if sides.len() == 4 && (sides[0].bends() || sides[2].bends()) {
        turned = sides[1..].iter().chain(&sides[..1]).cloned().collect();
        &turned[..]
    } else {
        sides
    };
    let exact = Transfinite::new(sides, tol)?;
    if let Some(at) = split_at(sides) {
        return split_patch(model, sides, &exact, at, tol);
    }
    if sides.len() >= 5 {
        return sector_patch(model, sides, &exact, tol);
    }
    let surface = interpolated(&Across { exact: &exact, tol }, tol)?;
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
        faces: vec![face],
        gap,
        angle,
        off,
    })
}

/// The patch through a hole of four sides or more, its corners off one
/// plane, wound counter-clockwise seen from outside the solid: the corners
/// blended by the regular polygon's Wachspress coordinates, with each
/// side's departure from its chord added, weighted by its two corners'
/// coordinates; a face on each side's sector as [`corner_patch`] builds
/// them. On each side of the polygon the coordinates are those of its two
/// corners alone and linear between them, and the sides beside it stand at
/// an end of theirs, so the patch takes each side exactly, straight or
/// curved; it meets its sides' supports at whatever angle they stand (a
/// chamfer's bevels share no tangent plane at a corner), and lies on one
/// plane where the corners and sides do. On a square the coordinates are
/// bilinear and the blend is the bilinearly blended Coons patch: four
/// straight sides give the bilinear patch.
///
/// # Errors
///
/// As [`corner_patch`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where there are fewer than four sides.
pub(crate) fn corner_patch_through_corners(
    model: &mut Model,
    sides: &[HoleSide],
    tol: Tolerances,
) -> OgeomResult<CornerPatch> {
    if sides.len() < 4 {
        ogeom_bail!(
            Construction,
            "a patch through its corners takes four sides or more, not {}",
            sides.len()
        );
    }
    let mut exact = Transfinite::new(sides, tol)?;
    let mut corners = Vec::with_capacity(sides.len());
    for ribbon in &exact.ribbons {
        corners.push(ribbon.point(0.0, tol)?);
    }
    exact.through = Some(corners);
    sector_patch(model, sides, &exact, tol)
}

/// The parameters of the interpolation points along each direction:
/// spaced by the cosine of the cosine, closest at the ends. Where a
/// corner's two supports curve differently, their twists there disagree in
/// the normal, which no choice of cross tangents mends; the exact patch
/// takes both, a polynomial only one, and it misses the other's tangent
/// plane over about its first span. The spans there fall as the square of
/// the cosine's.
fn net_parameters() -> Vec<f64> {
    net_parameters_of(NET)
}

/// [`net_parameters`] for a net of `net` points.
fn net_parameters_of(net: usize) -> Vec<f64> {
    #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
    let last = (net - 1) as f64;
    (0..net)
        .map(|k| {
            #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
            let t = k as f64 / last;
            if k == 0 {
                0.0
            } else if k == net - 1 {
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

/// What a patch's interpolation takes over the unit square: its points and
/// its derivatives across the square's four sides.
trait Data {
    fn point(&self, u: f64, v: f64) -> OgeomResult<Point>;
    fn dv_bottom(&self, u: f64) -> OgeomResult<Vector>;
    fn dv_top(&self, u: f64) -> OgeomResult<Vector>;
    fn du_left(&self, v: f64) -> OgeomResult<Vector>;
    fn du_right(&self, v: f64) -> OgeomResult<Vector>;

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

impl Across<'_> {
    fn square(&self) -> bool {
        self.exact.ribbons.len() == 4
    }
}

impl Data for Across<'_> {
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
}

/// The clamped bicubic B-spline taking the exact patch at [`NET`] by
/// [`NET`] points of the unit square and its derivatives across the
/// square's sides at each of them, its twists at the corners as
/// [`Data::twist`] gives them.
fn interpolated(across: &dyn Data, tol: Tolerances) -> OgeomResult<BSplineSurface> {
    interpolated_on(across, NET, tol)
}

/// [`interpolated`] on a net of `net` by `net` points.
fn interpolated_on(across: &dyn Data, net: usize, tol: Tolerances) -> OgeomResult<BSplineSurface> {
    let params = net_parameters_of(net);
    let mut knots = vec![0.0; 4];
    knots.extend_from_slice(&params[1..net - 1]);
    knots.extend_from_slice(&[1.0; 4]);
    let knots = KnotVector::new(knots, 3)?;
    let count = net + 2;
    // Rows: the point at the first parameter, the derivative there, the
    // points between, the derivative at the last, the point there.
    let mut a = vec![vec![0.0; count]; count];
    let mut rows: Vec<(f64, usize)> = Vec::with_capacity(count);
    rows.push((0.0, 0));
    rows.push((0.0, 1));
    for &t in &params[1..net - 1] {
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
    // The right-hand side along one line: its points at `params` and the
    // derivatives at its two ends, in the rows' order.
    let line = |points: Vec<Point>, start: Vector, end: Vector| -> Vec<Point> {
        let mut rhs = Vec::with_capacity(count);
        rhs.push(points[0]);
        rhs.push(Point::from_vector(start));
        rhs.extend_from_slice(&points[1..net - 1]);
        rhs.push(Point::from_vector(end));
        rhs.push(points[net - 1]);
        solve.apply(&rhs)
    };
    // Along `u` at each sampled `v`, the derivatives across `u = 0` and
    // `u = 1` at its ends.
    let mut along_u: Vec<Vec<Point>> = Vec::with_capacity(count);
    let mut rows_u: Vec<Vec<Point>> = Vec::with_capacity(net);
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

/// The interpolation points along each direction of each sector of a
/// patch of five sides or more: twice as fine as [`NET`], a sector's exact
/// patch turning faster across its narrow triangle.
const SECTOR_NET: usize = 2 * NET - 1;

/// The step of the differences a sector's derivatives are read by.
const STEP: f64 = 1e-5;

/// One triangle of a patch of five sides or more, from side `k`'s corners
/// `c_k` and `c_{k+1}` of the polygon to its centre: `(u, v)` of the unit
/// square at `(1 - v) ((1 - u) c_k + u c_{k+1})`, `v = 0` along the side and
/// `v = 1` the centre, a pole. Its derivatives are read from the exact
/// patch by differences of second order, one-sided into the square, and at
/// the centre from the exact patch's own derivatives there, which every
/// sector shares.
struct Sector<'a> {
    exact: &'a Transfinite,
    corners: [Point2; 2],
    /// The exact patch's derivatives at the centre along the polygon's `x`
    /// and `y`.
    centre: (Vector, Vector),
    tol: Tolerances,
}

impl Sector<'_> {
    /// The point of the polygon at `(u, v)`. On a side of the square, a
    /// corner's share is nought exactly, so neighbours read one line alike.
    fn at(&self, u: f64, v: f64) -> Point2 {
        let [a, b] = self.corners;
        Point2::new(
            ((1.0 - u) * a.x + u * b.x) * (1.0 - v),
            ((1.0 - u) * a.y + u * b.y) * (1.0 - v),
        )
    }

    fn exact_at(&self, u: f64, v: f64) -> OgeomResult<Vector> {
        Ok(self
            .exact
            .point_polygon(self.at(u, v), self.tol)?
            .to_vector())
    }

    /// The derivative at `from` stepping `+1` or `-1` along `u` (`along`)
    /// or `v`, one-sided, to second order.
    fn rate(&self, (u, v): (f64, f64), along: bool, sign: f64) -> OgeomResult<Vector> {
        let h = STEP * sign;
        let step = |k: f64| {
            if along {
                self.exact_at(u + h * k, v)
            } else {
                self.exact_at(u, v + h * k)
            }
        };
        Ok((step(1.0)? * 4.0 - step(0.0)? * 3.0 - step(2.0)?) / (2.0 * h))
    }
}

impl Data for Sector<'_> {
    fn point(&self, u: f64, v: f64) -> OgeomResult<Point> {
        self.exact.point_polygon(self.at(u, v), self.tol)
    }

    fn dv_bottom(&self, u: f64) -> OgeomResult<Vector> {
        self.rate((u, 0.0), false, 1.0)
    }

    fn dv_top(&self, u: f64) -> OgeomResult<Vector> {
        let [a, b] = self.corners;
        let (x, y) = ((1.0 - u) * a.x + u * b.x, (1.0 - u) * a.y + u * b.y);
        Ok(-(self.centre.0 * x + self.centre.1 * y))
    }

    fn du_left(&self, v: f64) -> OgeomResult<Vector> {
        self.rate((0.0, v), true, 1.0)
    }

    fn du_right(&self, v: f64) -> OgeomResult<Vector> {
        self.rate((1.0, v), true, -1.0)
    }
}

/// The patch of five sides or more: a face on each side's sector, and
/// what they achieved, measured on the faces: each against its side's
/// curve and support, and each against its neighbour along the line they
/// share.
fn sector_patch(
    model: &mut Model,
    sides: &[HoleSide],
    exact: &Transfinite,
    tol: Tolerances,
) -> OgeomResult<CornerPatch> {
    let n = sides.len();
    let corners = exact.polygon();
    let h = STEP;
    let at = |x: f64, y: f64| exact.point_polygon(Point2::new(x, y), tol);
    let centre = (
        (at(h, 0.0)? - at(-h, 0.0)?) / (2.0 * h),
        (at(0.0, h)? - at(0.0, -h)?) / (2.0 * h),
    );
    let count = 8 * (SECTOR_NET - 1);
    #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
    let sample = |i: usize| i as f64 / count as f64;
    let (mut gap, mut angle, mut off) = (0.0f64, 0.0f64, 0.0f64);
    let sine = |a: Vector, b: Vector| {
        (a.cross(b).magnitude() / (a.magnitude() * b.magnitude()))
            .clamp(0.0, 1.0)
            .asin()
    };
    let mut surfaces = Vec::with_capacity(n);
    for k in 0..n {
        let sector = Sector {
            exact,
            corners: [corners[k], corners[(k + 1) % n]],
            centre,
            tol,
        };
        let surface = interpolated_on(&sector, SECTOR_NET, tol)?;
        // Along its side, where the exact patch is the side's curve at the
        // same parameter; the angle read short of the corners.
        let ribbon = &exact.ribbons[k];
        for i in 0..=count {
            let t = sample(i);
            let on = ribbon.point(t, tol)?;
            gap = gap.max(surface.point_at(t, 0.0, tol)?.distance(on));
            if i == 0 || i == count {
                continue;
            }
            let (du, dv) = surface.d1_at(t, 0.0, tol)?;
            let support = support_normal(ribbon.side.support_at(t), on, tol)?;
            angle = angle.max(sine(du.cross(dv), support));
        }
        #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
        let step = 1.0 / (SECTOR_NET - 1) as f64;
        for i in 0..SECTOR_NET - 1 {
            for j in 0..SECTOR_NET - 1 {
                #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
                let (u, v) = ((i as f64 + 0.5) * step, (j as f64 + 0.5) * step);
                let p = surface.point_at(u, v, tol)?;
                off = off.max(p.distance(sector.point(u, v)?));
            }
        }
        surfaces.push(surface);
    }
    // Each sector against the next along the line from their shared corner
    // to the centre, the angle read off the corner and short of the pole.
    for k in 0..n {
        let (a, b) = (&surfaces[k], &surfaces[(k + 1) % n]);
        for i in 0..count {
            let v = sample(i);
            gap = gap.max(a.point_at(1.0, v, tol)?.distance(b.point_at(0.0, v, tol)?));
            if i == 0 {
                continue;
            }
            let (du, dv) = a.d1_at(1.0, v, tol)?;
            let (eu, ev) = b.d1_at(0.0, v, tol)?;
            angle = angle.max(sine(du.cross(dv), eu.cross(ev)));
        }
    }
    let faces = sector_faces(model, sides, &surfaces, gap, tol)?;
    Ok(CornerPatch {
        faces,
        gap,
        angle,
        off,
    })
}

/// The faces on `surfaces`, one on each side's sector: bounded by the
/// side's edge, the line from its last corner to the centre, a pole there,
/// and the line back from the centre to its first corner, each line one
/// edge shared by the two sectors beside it, on the first's `u = 0` and the
/// other's `u = 1`. The edges are widened to the gap measured.
fn sector_faces(
    model: &mut Model,
    sides: &[HoleSide],
    surfaces: &[BSplineSurface],
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let n = sides.len();
    let widened = Tolerance::new(gap + tol.confusion())?;
    let line = |from: Point2, to: Point2, (a, b): (f64, f64)| -> OgeomResult<BSpline2d> {
        BSpline2d::new(KnotVector::new(vec![a, a, b, b], 1)?, vec![from, to], tol)
    };
    let (o, x, y, xy) = (
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(1.0, 1.0),
    );
    let mut starts = Vec::with_capacity(n);
    for side in sides {
        let Some((p, q)) = ogeom_algo::edge_vertices(model, &side.edge)? else {
            ogeom_bail!(Construction, "a vertex patch's side has no vertices");
        };
        starts.push(if side.forward { p } else { q });
    }
    let ids: Vec<_> = surfaces
        .iter()
        .map(|s| model.geometry_mut().add_surface(s.clone().into()))
        .collect();
    let centre = ogeom_algo::make_vertex(model, surfaces[0].point_at(0.0, 1.0, tol)?).shape;
    model.widen(&centre, widened)?;
    let mut spokes = Vec::with_capacity(n);
    for k in 0..n {
        let curve = surfaces[k].iso_u_curve(0.0, tol)?;
        let edge = ogeom_algo::make_edge_between(
            model,
            Curve::BSpline(curve),
            (0.0, 1.0),
            &starts[k],
            &centre,
            tol,
        )?
        .shape;
        let before = (k + n - 1) % n;
        for (id, (from, to)) in [(ids[k], (o, y)), (ids[before], (x, xy))] {
            ogeom_algo::attach_pcurve(
                model,
                &edge,
                line(from, to, (0.0, 1.0))?.into(),
                id,
                Location::identity(),
                (0.0, 1.0),
            )?;
        }
        model.widen(&edge, widened)?;
        spokes.push(edge);
    }
    let mut faces = Vec::with_capacity(n);
    for (k, side) in sides.iter().enumerate() {
        let (a, b) = side.range;
        let (from, to) = if side.forward { (o, x) } else { (x, o) };
        ogeom_algo::attach_pcurve(
            model,
            &side.edge,
            line(from, to, (a, b))?.into(),
            ids[k],
            Location::identity(),
            (a, b),
        )?;
        model.widen(&side.edge, widened)?;
        let pole = {
            let mut data = EdgeData::new();
            data.degenerate = true;
            model.add_edge(data, &[centre.clone(), centre.clone()])?
        };
        ogeom_algo::attach_pcurve(
            model,
            &pole,
            line(xy, y, (0.0, 1.0))?.into(),
            ids[k],
            Location::identity(),
            (0.0, 1.0),
        )?;
        let walked = vec![
            if side.forward {
                side.edge.clone()
            } else {
                side.edge.reversed()
            },
            spokes[(k + 1) % n].clone(),
            pole,
            spokes[k].reversed(),
        ];
        let wire = ogeom_algo::make_wire(model, &walked, tol)?.shape;
        faces
            .push(ogeom_algo::make_face_on(model, ids[k], std::slice::from_ref(&wire), tol)?.shape);
    }
    for v in &starts {
        model.widen(v, widened)?;
    }
    Ok(faces)
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

            let support = support_normal(ribbon.side.support_at(s), on, tol)?;
            let sine = n.cross(support).magnitude() / n.magnitude();
            angle = angle.max(sine.clamp(0.0, 1.0).asin());
        }
    }
    Ok((gap, angle))
}

/// The interpolation points along each direction of each part of a split
/// patch: twice as fine as [`NET`], each part's exact patch turning faster.
const SPLIT_NET: usize = 4 * (NET - 1) + 1;

/// Where across the square a patch of four sides is split in two: where
/// its second and fourth sides each run from one edge into another, a
/// join of two surfaces bending opposite ways (a convex round's end running
/// on into a concave one's) that one smooth surface cannot follow, both at
/// the same `v`. `None` for any other.
fn split_at(sides: &[HoleSide]) -> Option<f64> {
    if sides.len() != 4 || sides[0].bends() || sides[2].bends() {
        return None;
    }
    let join = |s: &HoleSide| match s.pieces.as_slice() {
        [a, _] if s.bends() => Some(a.along.1),
        _ => None,
    };
    let (right, left) = (join(&sides[1])?, join(&sides[3])?);
    ((right - (1.0 - left)).abs() <= 1e-9).then_some(right)
}

/// The exact patch over `lo..hi` of the square's `v`, taken onto a square
/// of its own.
struct Part<'a> {
    across: &'a Across<'a>,
    lo: f64,
    hi: f64,
}

impl Part<'_> {
    /// The exact patch's derivative in `v` at `(u, v)`: the sides' own at
    /// the square's bottom and top, by differences a millionth apart
    /// between.
    fn rate(&self, u: f64, v: f64) -> OgeomResult<Vector> {
        if v <= 0.0 {
            return self.across.dv_bottom(u);
        }
        if v >= 1.0 {
            return self.across.dv_top(u);
        }
        let h = 1e-6;
        Ok((self.across.point(u, v + h)? - self.across.point(u, v - h)?) / (2.0 * h))
    }

    fn at(&self, w: f64) -> f64 {
        self.lo + (self.hi - self.lo) * w
    }
}

impl Data for Part<'_> {
    fn point(&self, u: f64, w: f64) -> OgeomResult<Point> {
        self.across.point(u, self.at(w))
    }

    fn dv_bottom(&self, u: f64) -> OgeomResult<Vector> {
        Ok(self.rate(u, self.lo)? * (self.hi - self.lo))
    }

    fn dv_top(&self, u: f64) -> OgeomResult<Vector> {
        Ok(self.rate(u, self.hi)? * (self.hi - self.lo))
    }

    fn du_left(&self, w: f64) -> OgeomResult<Vector> {
        self.across.du_left(self.at(w))
    }

    fn du_right(&self, w: f64) -> OgeomResult<Vector> {
        self.across.du_right(self.at(w))
    }
}

/// The patch of four sides split at `at` across the square: a face on each
/// part, below and above, each fitted on its own, so the join its second
/// and fourth sides run through falls on the line between them; and what
/// they achieved, each measured against its sides and the two against each
/// other along that line.
fn split_patch(
    model: &mut Model,
    sides: &[HoleSide],
    exact: &Transfinite,
    at: f64,
    tol: Tolerances,
) -> OgeomResult<CornerPatch> {
    let across = Across { exact, tol };
    let parts = [(0.0, at), (at, 1.0)];
    let mut surfaces = Vec::with_capacity(2);
    for (lo, hi) in parts {
        surfaces.push(interpolated_on(
            &Part {
                across: &across,
                lo,
                hi,
            },
            SPLIT_NET,
            tol,
        )?);
    }
    let sine = |a: Vector, b: Vector| {
        (a.cross(b).magnitude() / (a.magnitude() * b.magnitude()))
            .clamp(0.0, 1.0)
            .asin()
    };
    let chart = |v: f64| -> (usize, f64) {
        if v <= at {
            (0, v / at)
        } else {
            (1, (v - at) / (1.0 - at))
        }
    };
    let count = 8 * (NET - 1);
    #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
    let sample = |i: usize| i as f64 / count as f64;
    let (mut gap, mut angle, mut off) = (0.0f64, 0.0f64, 0.0f64);
    for (i, ribbon) in exact.ribbons.iter().enumerate() {
        for k in 0..=count {
            let s = sample(k);
            let (u, v) = match i {
                0 => (s, 0.0),
                1 => (1.0, s),
                2 => (1.0 - s, 1.0),
                _ => (0.0, 1.0 - s),
            };
            let (part, w) = chart(v);
            let on = ribbon.point(s, tol)?;
            gap = gap.max(surfaces[part].point_at(u, w, tol)?.distance(on));
            if k == 0 || k == count {
                continue;
            }
            let (du, dv) = surfaces[part].d1_at(u, w, tol)?;
            let support = support_normal(ribbon.side.support_at(s), on, tol)?;
            angle = angle.max(sine(du.cross(dv), support));
        }
    }
    for k in 0..=count {
        let u = sample(k);
        let (a, b) = (&surfaces[0], &surfaces[1]);
        gap = gap.max(a.point_at(u, 1.0, tol)?.distance(b.point_at(u, 0.0, tol)?));
        if k == 0 || k == count {
            continue;
        }
        let (du, dv) = a.d1_at(u, 1.0, tol)?;
        let (eu, ev) = b.d1_at(u, 0.0, tol)?;
        angle = angle.max(sine(du.cross(dv), eu.cross(ev)));
    }
    #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
    let step = 1.0 / (SPLIT_NET - 1) as f64;
    for (surface, (lo, hi)) in surfaces.iter().zip(parts) {
        for i in 0..SPLIT_NET - 1 {
            for j in 0..SPLIT_NET - 1 {
                #[expect(clippy::cast_precision_loss, reason = "a net of thirty-three")]
                let (u, w) = ((i as f64 + 0.5) * step, (j as f64 + 0.5) * step);
                let p = surface.point_at(u, w, tol)?;
                off = off.max(p.distance(across.point(u, lo + (hi - lo) * w)?));
            }
        }
    }
    let faces = split_faces(model, sides, &surfaces, at, gap, tol)?;
    Ok(CornerPatch {
        faces,
        gap,
        angle,
        off,
    })
}

/// The faces of [`split_patch`] on `surfaces`, below and above `at`: the
/// lower bounded by the first side, the second's edges up to its join, the
/// line across between the joins (one edge on the lower's top and the
/// upper's bottom) and the fourth's edges from its join; the upper by that
/// line, the rest of the second's edges, the third side and the rest of the
/// fourth's. The edges are widened to the gap measured.
fn split_faces(
    model: &mut Model,
    sides: &[HoleSide],
    surfaces: &[BSplineSurface],
    at: f64,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let widened = Tolerance::new(gap + tol.confusion())?;
    let line = |from: Point2, to: Point2, (a, b): (f64, f64)| -> OgeomResult<BSpline2d> {
        BSpline2d::new(KnotVector::new(vec![a, a, b, b], 1)?, vec![from, to], tol)
    };
    let ids: Vec<_> = surfaces
        .iter()
        .map(|s| model.geometry_mut().add_surface(s.clone().into()))
        .collect();
    // Each join's vertex: where the side's first edge ends.
    let mut joins = Vec::with_capacity(2);
    for side in [&sides[1], &sides[3]] {
        let first = &side.pieces[0];
        let Some((p, q)) = ogeom_algo::edge_vertices(model, &first.edge)? else {
            ogeom_bail!(Construction, "a vertex patch's side has no vertices");
        };
        let join = if first.forward { q } else { p };
        model.widen(&join, widened)?;
        joins.push(join);
    }
    let across = ogeom_algo::make_edge_between(
        model,
        Curve::BSpline(surfaces[0].iso_v_curve(1.0, tol)?),
        (0.0, 1.0),
        &joins[1],
        &joins[0],
        tol,
    )?
    .shape;
    for (id, v) in [(ids[0], 1.0), (ids[1], 0.0)] {
        ogeom_algo::attach_pcurve(
            model,
            &across,
            line(Point2::new(0.0, v), Point2::new(1.0, v), (0.0, 1.0))?.into(),
            id,
            Location::identity(),
            (0.0, 1.0),
        )?;
    }
    model.widen(&across, widened)?;
    // Where a side's point `s` along it stands on the square.
    let square = |i: usize, s: f64| match i {
        0 => (s, 0.0),
        1 => (1.0, s),
        2 => (1.0 - s, 1.0),
        _ => (0.0, 1.0 - s),
    };
    let reach = 1e-9;
    let mut faces = Vec::with_capacity(2);
    for (k, (lo, hi)) in [(0.0, at), (at, 1.0)].into_iter().enumerate() {
        let mut walked = Vec::new();
        if k == 1 {
            walked.push(across.clone());
        }
        for (i, side) in sides.iter().enumerate() {
            for (edge, (a, b), forward, (s0, s1)) in side.edges() {
                let (v0, v1) = (square(i, s0).1, square(i, s1).1);
                if v0.min(v1) < lo - reach || v0.max(v1) > hi + reach {
                    continue;
                }
                let chart = |s: f64| {
                    let (u, v) = square(i, s);
                    Point2::new(u, (v - lo) / (hi - lo))
                };
                let (start, end) = if forward {
                    (chart(s0), chart(s1))
                } else {
                    (chart(s1), chart(s0))
                };
                ogeom_algo::attach_pcurve(
                    model,
                    &edge,
                    line(start, end, (a, b))?.into(),
                    ids[k],
                    Location::identity(),
                    (a, b),
                )?;
                model.widen(&edge, widened)?;
                if let Some((p, q)) = ogeom_algo::edge_vertices(model, &edge)? {
                    model.widen(&p, widened)?;
                    model.widen(&q, widened)?;
                }
                walked.push(if forward { edge } else { edge.reversed() });
            }
            if k == 0 && i == 1 {
                walked.push(across.reversed());
            }
        }
        let wire = ogeom_algo::make_wire(model, &walked, tol)?.shape;
        faces
            .push(ogeom_algo::make_face_on(model, ids[k], std::slice::from_ref(&wire), tol)?.shape);
    }
    Ok(faces)
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
        let at = |x: f64| Point2::new(from.x + (to.x - from.x) * x, from.y + (to.y - from.y) * x);
        for (edge, (a, b), forward, (s0, s1)) in side.edges() {
            let (start, end) = if forward {
                (at(s0), at(s1))
            } else {
                (at(s1), at(s0))
            };
            let line =
                BSpline2d::new(KnotVector::new(vec![a, a, b, b], 1)?, vec![start, end], tol)?;
            ogeom_algo::attach_pcurve(model, &edge, line.into(), id, Location::identity(), (a, b))?;
            model.widen(&edge, widened)?;
            if let Some((p, q)) = ogeom_algo::edge_vertices(model, &edge)? {
                model.widen(&p, widened)?;
                model.widen(&q, widened)?;
            }
            walked.push(if forward { edge } else { edge.reversed() });
        }
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
                    pieces: Vec::new(),
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
                    pieces: Vec::new(),
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
                    pieces: Vec::new(),
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

    /// Five straight sides of a pentagon on `z = 0`, each supported by that
    /// plane: the five-sided patch lies in it, exactly, takes each side at
    /// its own parameter, and its five faces meet their sides and each
    /// other within the fit target.
    #[test]
    fn a_flat_five_sided_hole_is_patched_flat_on_five_faces() {
        let mut m = Model::new();
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, -0.5, 0.0),
            Point::new(5.0, 2.5, 0.0),
            Point::new(2.0, 4.5, 0.0),
            Point::new(-0.5, 3.0, 0.0),
        ];
        let plane: SurfaceGeometry =
            ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(ogeom_math::Frame::WORLD)).into();
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| ogeom_algo::make_vertex(&mut m, *p).shape)
            .collect();
        let sides: Vec<HoleSide> = (0..5)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 5]);
                let line = ogeom_geom::LineCurve::segment(a, b, T).unwrap();
                let range = (0.0, a.distance(b));
                let edge = ogeom_algo::make_edge_between(
                    &mut m,
                    Curve::Line(line),
                    range,
                    &vertices[i],
                    &vertices[(i + 1) % 5],
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
                    pieces: Vec::new(),
                }
            })
            .collect();
        let patch = Transfinite::new(&sides, T).unwrap();
        let polygon = patch.polygon();
        for i in 0..=10 {
            for j in 0..=10 {
                let p = Point2::new(-1.0 + 0.2 * f64::from(i), -1.0 + 0.2 * f64::from(j));
                // Inside the polygon: on the inner side of every side.
                let inside = (0..5).all(|k| {
                    let (a, b) = (polygon[k], polygon[(k + 1) % 5]);
                    (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) >= 0.0
                });
                if inside {
                    let q = patch.point_polygon(p, T).unwrap();
                    assert!(q.z.abs() < 1e-12, "{p:?}: {q:?}");
                }
            }
        }
        for k in 0..5 {
            let (a, b) = (polygon[k], polygon[(k + 1) % 5]);
            for i in 0..=10 {
                let t = f64::from(i) / 10.0;
                let on = Point2::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
                let q = patch.point_polygon(on, T).unwrap();
                let want = corners[k].lerp(corners[(k + 1) % 5], t);
                assert!(q.distance(want) < 1e-12, "side {k} at {t}: {q:?}");
            }
        }
        let built = corner_patch(&mut m, &sides, T).unwrap();
        assert_eq!(built.faces.len(), 5);
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
        let proved =
            ogeom_algo::certify_face_within(&m, &patch.faces[0], &octant(), FIT, T).unwrap();
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
        let proved =
            ogeom_algo::certify_face_within(&m, &patch.faces[0], &octant(), FIT, T).unwrap();
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
        let proved = ogeom_algo::certify_face_within(&m, &patch.faces[0], &region, FIT, T).unwrap();
        assert_eq!(proved, Ok(()), "off the exact patch {}", patch.off);
    }
}
