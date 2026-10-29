//! Where an open wire's walls meet inside a corner against an arc.
//!
//! Where a spine of lines and arcs turns a corner in the plane of the arcs
//! beside it, the rotation-minimizing frame keeps the section's frame in
//! that plane: a profile point at height `w` above it and at offset `u`
//! outside the turn (negative inside) runs at that height along the line
//! offset by `u` from a line leg, and round the circle of radius `R + u`
//! or `R - u` about an arc leg's axis, `R` the arc's radius. Inside the
//! turn the two legs' paths for one point cross on the legs themselves,
//! at the point `X(u)` of the corner's plane lifted by `w`, where each leg
//! runs on until it meets the other leg's wall.
//!
//! The crossings of every offset are the points at one signed distance
//! from both legs' tracks: between a line and a circle a parabola with its
//! focus at the circle's centre, between two circles an ellipse or a
//! hyperbola with the centres as foci. A straight profile edge at one
//! height meets on that conic, lifted. Any other straight edge has `u`
//! affine in `w`: a line leg's wall is a plane, an arc leg's the cone or
//! cylinder of radius `R ± u` about its axis, and the equations of two
//! such cones about parallel axes, squared, differ by an affine one, since
//! their radii change with `u` at one rate. Either way the edge's meeting
//! is a plane's section of a cone or of a cylinder: a conic, or a line
//! along the corner's axis where the edge stands square to the corner's
//! plane. A conic is held exactly as a rational quadratic through its
//! ends, the meeting of its end tangents and one point between, and
//! checked against the two tracks at points between.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSplineCurve, Curve, Curve3d as _, LineCurve};
use ogeom_math::{KnotVector, Point, Vector, Weighted};

/// A leg's paths near a corner, in the plane the corner turns in.
#[derive(Debug, Clone, Copy)]
enum Track {
    /// The lines along `heading` offset along `outward`, the unit normal
    /// pointing away from the turn.
    Line { heading: Vector, outward: Vector },
    /// The circles about `centre`, in the corner's plane, of radius
    /// `radius + sign * u` for the offset `u`.
    Circle {
        centre: Point,
        radius: f64,
        sign: f64,
    },
}

/// A corner against an arc turning in the arc's plane, read for the points
/// of the section inside it.
#[derive(Debug, Clone, Copy)]
pub(super) struct InsideTurn {
    /// The corner.
    at: Point,
    /// The unit normal of the plane the corner turns in, along its axis.
    up: Vector,
    /// The unit normal of the leg into the corner there, away from the
    /// turn.
    outward: Vector,
    /// The leg into the corner and the leg out of it.
    tracks: [Track; 2],
}

impl InsideTurn {
    /// The corner at `at` turning about `up` from the heading `d` onto
    /// `n`, each leg a line or an arc about an axis through the point given.
    pub(super) fn new(
        at: Point,
        up: Vector,
        (d, n): (Vector, Vector),
        centres: [Option<Point>; 2],
    ) -> Self {
        let unit = |v: Vector| v / v.magnitude();
        let outward = [unit(d * n.dot(d) - n), unit(d - n * d.dot(n))];
        let headings = [d, n];
        let tracks = [0, 1].map(|j| match centres[j] {
            None => Track::Line {
                heading: headings[j],
                outward: outward[j],
            },
            Some(on_axis) => {
                let centre = on_axis + up * (at - on_axis).dot(up);
                let radial = at - centre;
                Track::Circle {
                    centre,
                    radius: radial.magnitude(),
                    sign: if radial.dot(outward[j]) > 0.0 {
                        1.0
                    } else {
                        -1.0
                    },
                }
            }
        });
        Self {
            at,
            up,
            outward: outward[0],
            tracks,
        }
    }

    /// A point's offset outside the turn and its height, where the section
    /// stands at the corner.
    pub(super) fn offset_and_height(&self, p: Point) -> (f64, f64) {
        let r = p - self.at;
        (r.dot(self.outward), r.dot(self.up))
    }

    /// `x`'s offset from leg `j`'s track and the unit gradient of that
    /// offset, in the corner's plane.
    fn offset_from(&self, j: usize, x: Point) -> (f64, Vector) {
        match self.tracks[j] {
            Track::Line { outward, .. } => ((x - self.at).dot(outward), outward),
            Track::Circle {
                centre,
                radius,
                sign,
            } => {
                let r = x - centre;
                let flat = r - self.up * r.dot(self.up);
                let m = flat.magnitude();
                (sign * (m - radius), flat * (sign / m))
            }
        }
    }

    /// Where the two legs' paths for the offset `u` cross, in the corner's
    /// plane: the crossing that runs into the corner as `u` runs to zero,
    /// or `None` where the paths miss each other.
    fn crossing(&self, u: f64) -> Option<Point> {
        // The root of `t² + 2bt + c = 0` that vanishes with `c`, the one
        // the corner's own point continues.
        let line_circle = |base: Point, along: Vector, centre: Point, r: f64| -> Option<Point> {
            if r <= 0.0 {
                return None;
            }
            let q = base - centre;
            let b = q.dot(along);
            let c = q.dot(q) - r * r;
            let disc = b * b - c;
            if disc < 0.0 {
                return None;
            }
            let root = b + b.signum() * disc.sqrt();
            (root != 0.0).then(|| base + along * (-c / root))
        };
        match (self.tracks[0], self.tracks[1]) {
            (
                Track::Line { heading, outward },
                Track::Circle {
                    centre,
                    radius,
                    sign,
                },
            )
            | (
                Track::Circle {
                    centre,
                    radius,
                    sign,
                },
                Track::Line { heading, outward },
            ) => line_circle(self.at + outward * u, heading, centre, radius + sign * u),
            (
                Track::Circle {
                    centre: c1,
                    radius: r1,
                    sign: s1,
                },
                Track::Circle {
                    centre: c2,
                    radius: r2,
                    sign: s2,
                },
            ) => {
                let (r1, r2) = (r1 + s1 * u, r2 + s2 * u);
                let e = c2 - c1;
                let apart = e.magnitude();
                if r1 <= 0.0 || r2 <= 0.0 || apart == 0.0 {
                    return None;
                }
                let ex = e / apart;
                let ey = self.up.cross(ex);
                let a = (apart * apart + r1 * r1 - r2 * r2) / (2.0 * apart);
                let h2 = r1 * r1 - a * a;
                if h2 < 0.0 {
                    return None;
                }
                let side = if (self.at - c1).dot(ey) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                Some(c1 + ex * a + ey * (side * h2.sqrt()))
            }
            (Track::Line { .. }, Track::Line { .. }) => None,
        }
    }

    /// Where the two legs' paths of the point `p` of the section at the
    /// corner meet, or `None` where they miss each other.
    pub(super) fn meeting(&self, p: Point) -> Option<Point> {
        let (u, w) = self.offset_and_height(p);
        self.crossing(u).map(|x| x + self.up * w)
    }

    /// The curve the two legs' walls swept by the straight profile edge
    /// from `a` to `b` (where the section stands at the corner) meet on,
    /// from `a`'s meeting to `b`'s.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where a point's paths miss each other, or the conic cannot be held.
    pub(super) fn meeting_curve(&self, a: Point, b: Point, tol: Tolerances) -> OgeomResult<Curve> {
        let (ua, wa) = self.offset_and_height(a);
        let (ub, wb) = self.offset_and_height(b);
        let (du, dw) = (ub - ua, wb - wa);
        let at = |s: f64| -> OgeomResult<Point> {
            let Some(x) = self.crossing(ua + du * s) else {
                ogeom_bail!(
                    Construction,
                    "inside a corner against an arc, a point's paths down the \
                     two legs miss each other"
                );
            };
            Ok(x + self.up * (wa + dw * s))
        };
        // The edge's crossings lie on one line along the corner's axis
        // where the edge stands square to the corner's plane.
        if du.abs() <= tol.confusion() {
            return Ok(LineCurve::segment(at(0.0)?, at(1.0)?, tol)?.into());
        }
        // `X'(u)` meets both offsets' gradients at one: `g1 · X' = g2 · X'
        // = 1`.
        let slope = |s: f64| -> OgeomResult<Vector> {
            let x = at(s)?;
            let (g1, g2) = (self.offset_from(0, x).1, self.offset_from(1, x).1);
            let across = self.up.dot(g1.cross(g2));
            if across.abs() <= tol.angular() {
                ogeom_bail!(
                    Construction,
                    "inside a corner against an arc, the two legs' paths run \
                     tangent where they cross"
                );
            }
            Ok(self.up.cross(g1 - g2) * (du / across) + self.up * dw)
        };
        let mut spans: Vec<[Weighted<Point>; 2]> = Vec::new();
        let mut pending = vec![(0.0_f64, 1.0_f64, 0_u32)];
        let mut first = None;
        while let Some((s0, s1, depth)) = pending.pop() {
            let (p0, p2) = (at(s0)?, at(s1)?);
            first.get_or_insert(p0);
            match rational_span(
                (p0, slope(s0)?),
                (p2, slope(s1)?),
                at(f64::midpoint(s0, s1))?,
                tol,
            ) {
                Some(middle) => spans.push([
                    middle,
                    Weighted {
                        scaled: p2,
                        weight: 1.0,
                    },
                ]),
                // A span turning half a revolution or more has no single
                // quadratic: halved, each half turns less.
                None if depth < 8 => {
                    let mid = f64::midpoint(s0, s1);
                    pending.push((mid, s1, depth + 1));
                    pending.push((s0, mid, depth + 1));
                }
                None => ogeom_bail!(
                    Construction,
                    "the walls' meeting inside a corner against an arc has no \
                     rational quadratic form"
                ),
            }
        }
        let Some(p0) = first else {
            ogeom_bail!(Construction, "the walls' meeting has no span");
        };
        let mut control = vec![Weighted {
            scaled: p0,
            weight: 1.0,
        }];
        let mut knots = vec![0.0; 3];
        #[allow(clippy::cast_precision_loss, reason = "a few spans")]
        let count = spans.len() as f64;
        for (k, span) in spans.iter().enumerate() {
            control.extend(*span);
            if k + 1 < spans.len() {
                #[allow(clippy::cast_precision_loss, reason = "a few spans")]
                let knot = (k + 1) as f64 / count;
                knots.extend([knot, knot]);
            }
        }
        knots.extend([1.0; 3]);
        let curve: Curve = BSplineCurve::rational(KnotVector::new(knots, 2)?, control)?.into();
        // Held to what it stands for: each point at one offset from both
        // tracks, that offset and its height a point of the edge.
        let length = du.hypot(dw);
        for k in 1..16 {
            let y = curve.point_at(f64::from(k) / 16.0, tol)?;
            let (u1, u2) = (self.offset_from(0, y).0, self.offset_from(1, y).0);
            let w = (y - self.at).dot(self.up);
            let off_edge = ((u1 - ua) * dw - (w - wa) * du).abs() / length;
            if (u1 - u2).abs() > tol.confusion() || off_edge > tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "the walls' meeting inside a corner against an arc strays \
                     from the conic it lies on"
                );
            }
        }
        Ok(curve)
    }
}

/// The middle control point of the rational quadratic from `p0` to `p2`
/// with those end tangents through `s`: the end tangents' meeting, weighted
/// by `s`'s barycentric coordinates `(b0, b1, b2)` there as `b1 / (2 √(b0
/// b2))`. `None` where the tangents do not meet ahead of `p0` and behind
/// `p2`, or `s` lies outside their triangle.
fn rational_span(
    (p0, t0): (Point, Vector),
    (p2, t2): (Point, Vector),
    s: Point,
    tol: Tolerances,
) -> Option<Weighted<Point>> {
    let w0 = p0 - p2;
    let (a, b, c) = (t0.dot(t0), t0.dot(t2), t2.dot(t2));
    let (d, e) = (t0.dot(w0), t2.dot(w0));
    let den = a * c - b * b;
    if den <= a * c * tol.angular() * tol.angular() {
        return None;
    }
    let alpha = (b * e - c * d) / den;
    let beta = (a * e - b * d) / den;
    let (from_p0, from_p2) = (p0 + t0 * alpha, p2 + t2 * beta);
    if alpha <= 0.0 || beta >= 0.0 || from_p0.distance(from_p2) > tol.confusion() {
        return None;
    }
    let p1 = from_p0.midpoint(from_p2);
    let n = (p1 - p0).cross(p2 - p0);
    let nn = n.dot(n);
    if nn <= 0.0 {
        return None;
    }
    let bary = |x: Point, y: Point| (x - s).cross(y - s).dot(n) / nn;
    let (b0, b1, b2) = (bary(p1, p2), bary(p2, p0), bary(p0, p1));
    if b0 <= 0.0 || b1 <= 0.0 || b2 <= 0.0 {
        return None;
    }
    let weight = b1 / (2.0 * (b0 * b2).sqrt());
    Some(Weighted {
        scaled: Point::from_vector(p1.to_vector() * weight),
        weight,
    })
}
