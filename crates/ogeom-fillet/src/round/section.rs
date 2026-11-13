//! P1 for a round between two planes: the section of a blend on a straight
//! crease, solved on the original faces alone.
//!
//! A section is read off the crease's two planes and nothing else, so no
//! chain's section depends on another's (the order independence of ADR
//! 0001, §4.2). It answers what P2 and P4 ask of a strip: where it touches
//! each face along the crease (its contact curves), its section in the
//! plane square to the crease at any station (the arc a strip ends on at a
//! vertex blend), and the strip's surface over a window of the crease,
//! which may run past the edge's ends ("way long", §4.2), its true ends
//! being P2's.
//!
//! Three sections, their maths the current path's, moved over unchanged:
//!
//! - the rolling ball of radius `r`: [`crate::fillet::ball_seat`], the
//!   arithmetic `seated_fillet` builds its cylinder with, so the two paths'
//!   cylinders agree to the bit;
//! - a conic of two distances, elliptical or of shape `rho`:
//!   [`crate::asym`]'s rational arc, swept along the crease;
//! - distances following a law each: [`crate::asym_law::varying_net`].

use crate::asym::{RoundProfile, conic_arc, conic_pieces, legs};
use crate::fillet::{BallSeat, ball_seat};
use crate::law::Profile;
use crate::support::Seat;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSplineCurve, CircleCurve, Curve, CylinderSurface, ExtrusionSurface, LineCurve, SurfaceGeometry,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, KnotVector, Point, Vector};

/// What a round's section is, as a request states it.
#[derive(Debug, Clone)]
pub(crate) enum RoundSpec {
    /// The rolling ball of this radius.
    Circular(f64),
    /// A conic meeting each face of the seat at its distance, in the seat's
    /// face order.
    Conic {
        /// The distance along each face.
        d: [f64; 2],
        /// The section's shape.
        profile: RoundProfile,
    },
    /// Distances following a law each, by the fraction of the edge from
    /// the seat's start, in the seat's face order.
    Varying {
        /// The law on each face.
        laws: [Profile; 2],
        /// The section's shape.
        profile: RoundProfile,
    },
}

/// A round's section on a straight seat between two planes.
#[derive(Debug, Clone)]
pub(crate) struct PlanarSection {
    /// The seat, its crease from `start` to `end`.
    seat: Seat,
    /// On each face, the unit leg square to the crease, away from it.
    legs: [Vector; 2],
    /// The section.
    kind: Kind,
}

/// The solved section.
#[derive(Debug, Clone)]
enum Kind {
    /// The ball and its radius.
    Ball(BallSeat, f64),
    /// A conic of constant distances and weight.
    Conic([f64; 2], f64),
    /// Varying distances, and the section's weight.
    Varying(Box<[Profile; 2]>, f64),
}

impl PlanarSection {
    /// The section `spec` asks on `seat`; refused where the faces meet too
    /// sharply or too flat, or a distance is not a usable length.
    pub(crate) fn new(seat: &Seat, spec: &RoundSpec, tol: Tolerances) -> OgeomResult<Self> {
        let legs = legs(seat, tol)?;
        let cos = legs[0].dot(legs[1]).clamp(-1.0, 1.0);
        let usable = |d: f64| d.is_finite() && d > tol.confusion();
        let kind = match spec {
            RoundSpec::Circular(r) => {
                if !usable(*r) {
                    ogeom_bail!(Construction, "a round of radius {r} rounds nothing");
                }
                Kind::Ball(ball_seat(seat, *r, tol)?, *r)
            }
            RoundSpec::Conic { d, profile } => {
                if !d.iter().all(|d| usable(*d)) {
                    ogeom_bail!(
                        Construction,
                        "a round set back {d:?} from its edge rounds nothing"
                    );
                }
                Kind::Conic(*d, profile.weight(cos)?)
            }
            RoundSpec::Varying { laws, profile } => {
                Kind::Varying(Box::new(laws.clone()), profile.weight(cos)?)
            }
        };
        Ok(Self {
            seat: seat.clone(),
            legs,
            kind,
        })
    }

    /// The seat.
    pub(crate) const fn seat(&self) -> &Seat {
        &self.seat
    }

    /// The crease's length, from the seat's start to its end.
    pub(crate) fn length(&self) -> f64 {
        self.seat.start.distance(self.seat.end)
    }

    /// The point of the crease `s` along it from the seat's start.
    pub(crate) fn crease(&self, s: f64) -> Point {
        self.seat.start + self.seat.along * s
    }

    /// How far from the crease the section touches each face, in the
    /// seat's face order, `s` along the crease. A varying section holds
    /// its laws' end distances past the edge's ends.
    pub(crate) fn distances(&self, s: f64) -> [f64; 2] {
        match &self.kind {
            Kind::Ball(ball, r) => {
                let mut d = [0.0; 2];
                for (i, n) in [(ball.first, ball.n1), (ball.second, ball.n2)] {
                    d[i] = (ball.offset + n * (r * ball.sign)).dot(self.legs[i]);
                }
                d
            }
            Kind::Conic(d, _) => *d,
            Kind::Varying(laws, _) => {
                let t = (s / self.length()).clamp(0.0, 1.0);
                [laws[0].at(t).0, laws[1].at(t).0]
            }
        }
    }

    /// Where the section touches face `i` of the seat, `s` along the
    /// crease.
    pub(crate) fn contact(&self, i: usize, s: f64) -> Point {
        match &self.kind {
            Kind::Ball(ball, r) => {
                let n = if i == ball.first { ball.n1 } else { ball.n2 };
                self.crease(s) + ball.offset + n * (r * ball.sign)
            }
            _ => self.crease(s) + self.legs[i] * self.distances(s)[i],
        }
    }

    /// The section's contact curve on face `i` over the window `(s0, s1)`
    /// of the crease, and its range: a line parallel to the crease for a
    /// constant section, the law's cubics carried along the face for a
    /// varying one (within the edge, its laws held level past its ends).
    pub(crate) fn contact_curve(
        &self,
        i: usize,
        (s0, s1): (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        match &self.kind {
            Kind::Varying(laws, w) => {
                let length = self.length();
                if s0 < -tol.confusion() || s1 > length + tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a varying round's contact is asked past its edge's ends"
                    );
                }
                let net = crate::asym_law::varying_net(
                    &self.seat,
                    [&laws[0], &laws[1]],
                    self.legs,
                    *w,
                    tol,
                )?;
                let column = if i == net.first {
                    0
                } else {
                    net.rows[0].len() - 1
                };
                let control: Vec<Point> = net.rows.iter().map(|r| r[column]).collect();
                let curve = Curve::BSpline(BSplineCurve::new(
                    KnotVector::new(net.u_knots, 3)?,
                    control,
                    tol,
                )?);
                Ok((curve, (s0 / length, s1 / length)))
            }
            _ => {
                let from = self.contact(i, 0.0);
                let line = LineCurve::new(ogeom_math::Axis::new(
                    from,
                    Direction::new(self.seat.along, tol)?,
                ));
                Ok((Curve::Line(line), (s0, s1)))
            }
        }
    }

    /// The face of the seat the section's arc starts on, and the one it
    /// ends on: the order that turns positively about the crease.
    pub(crate) fn order(&self) -> [usize; 2] {
        let positive = self
            .seat
            .along
            .dot(self.seat.normals[0].cross(self.seat.normals[1]))
            .is_sign_positive();
        if positive { [0, 1] } else { [1, 0] }
    }

    /// The strip's end through its contacts on the seat's faces `faces[0]`
    /// and `faces[1]` at stations `at[0]` and `at[1]` along the crease: the
    /// strip's surface cut by the plane through both that holds square to
    /// the crease the chord of the sections there. That is the arc at the
    /// first station carried along the crease onto the plane, an affine map
    /// that keeps a rational quadratic arc's control points' weights, so the
    /// end is exact; and over `(0, 1)`. For a circular section only.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where the contacts coincide or the section is not circular.
    pub(crate) fn skew_end(
        &self,
        faces: [usize; 2],
        at: [f64; 2],
        tol: Tolerances,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        if !matches!(self.kind, Kind::Ball(..)) {
            ogeom_bail!(Construction, "a skew end of a section that is not circular");
        }
        let along = self.seat.along;
        let (p0, p1) = (self.contact(faces[0], at[0]), self.contact(faces[1], at[1]));
        // The arc at the first station: from `p0` to `q`, its middle control
        // point on the crease, weighted by the cosine of half its turn.
        let q = self.contact(faces[1], at[0]);
        let apex = self.crease(at[0]);
        let (u0, u1) = ((apex - p0).normalized(tol)?, (q - apex).normalized(tol)?);
        let weight = ((1.0 + u0.dot(u1)) / 2.0).max(0.0).sqrt();
        // The plane through `p0` and `p1` holding the chord's square in it.
        let square = along.cross(q - p0);
        let normal = (p1 - p0).cross(square);
        let onto = |x: Point| x + along * (normal.dot(p0 - x) / normal.dot(along));
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
        let control = vec![
            ogeom_math::Weighted::new(p0, 1.0, tol)?,
            ogeom_math::Weighted::new(onto(apex), weight, tol)?,
            ogeom_math::Weighted::new(p1, 1.0, tol)?,
        ];
        Ok((
            Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?),
            (0.0, 1.0),
        ))
    }

    /// The section in the plane square to the crease `s` along it, from
    /// its contact on the first face of [`PlanarSection::order`] to its
    /// contact on the second, and its range.
    pub(crate) fn arc_at(&self, s: f64, tol: Tolerances) -> OgeomResult<(Curve, (f64, f64))> {
        let [first, second] = self.order();
        match &self.kind {
            Kind::Ball(ball, r) => {
                let along = Direction::new(self.seat.along, tol)?;
                let radial = Direction::new(ball.n1 * ball.sign, tol)?;
                let frame = Frame::new(self.crease(s) + ball.offset, along, radial, tol)?;
                let circle = Circle::new(frame, *r, tol)?;
                Ok((Curve::Circle(CircleCurve::new(circle)), (0.0, ball.sweep)))
            }
            Kind::Conic(_, w) | Kind::Varying(_, w) => {
                let p = [
                    self.contact(first, s),
                    self.crease(s),
                    self.contact(second, s),
                ];
                Ok((conic_arc(p, *w, tol)?, (0.0, 1.0)))
            }
        }
    }

    /// The strip's surface over the window `(s0, s1)` of the crease: the
    /// ball's cylinder (its axis's origin `s0` along the crease, its chart
    /// across from the first face's contact), the conic's extrusion, or
    /// the varying band within the edge. Every one is exact.
    pub(crate) fn surface(
        &self,
        (s0, s1): (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<SurfaceGeometry> {
        if s1 - s0 <= tol.confusion() {
            ogeom_bail!(Construction, "a strip over no length of its crease");
        }
        let along = Direction::new(self.seat.along, tol)?;
        match &self.kind {
            Kind::Ball(ball, r) => {
                // The seam half a turn from the strip's middle, so that the
                // strip's chart lies inside the cylinder's own, clear of it.
                let radial = Direction::new((ball.n1 + ball.n2) * -ball.sign, tol)?;
                let frame = Frame::new(self.crease(s0) + ball.offset, along, radial, tol)?;
                Ok(CylinderSurface::new(Cylinder::new(frame, *r, tol)?, (0.0, s1 - s0))?.into())
            }
            Kind::Conic(_, _) => {
                let (arc, _) = self.arc_at(s0, tol)?;
                Ok(ExtrusionSurface::new(arc, along, s1 - s0)?.into())
            }
            Kind::Varying(laws, w) => {
                let length = self.length();
                if s0.abs() > tol.confusion() || (s1 - length).abs() > tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a varying round's strip is built over its edge, end to end"
                    );
                }
                let net = crate::asym_law::varying_net(
                    &self.seat,
                    [&laws[0], &laws[1]],
                    self.legs,
                    *w,
                    tol,
                )?;
                Ok(net.surface)
            }
        }
    }

    /// The control points of the section `s` along the crease, its middle
    /// weight and knots across, where it is a conic: the pieces the strip
    /// surface's rows are made of.
    pub(crate) fn conic_net_at(&self, s: f64) -> Option<(Vec<Point>, f64, Vec<f64>)> {
        let [first, second] = self.order();
        match &self.kind {
            Kind::Ball(..) => None,
            Kind::Conic(_, w) | Kind::Varying(_, w) => Some(conic_pieces(
                [
                    self.contact(first, s),
                    self.crease(s),
                    self.contact(second, s),
                ],
                *w,
            )),
        }
    }

    /// The section's area square to the crease `s` along it: the corner
    /// between the faces and the arc, which a convex strip takes off and a
    /// concave one adds.
    pub(crate) fn corner_area(&self, s: f64) -> f64 {
        let d = self.distances(s);
        let cos = self.legs[0].dot(self.legs[1]).clamp(-1.0, 1.0);
        let sin = self.legs[0].cross(self.legs[1]).magnitude();
        let triangle = 0.5 * d[0] * d[1] * sin;
        let w = match &self.kind {
            Kind::Ball(..) => ((1.0 - cos) * 0.5).sqrt(),
            Kind::Conic(_, w) | Kind::Varying(_, w) => *w,
        };
        triangle * (1.0 - segment_share(w))
    }
}

/// The share of its control triangle the segment between a rational
/// quadratic arc of middle weight `w` and its chord takes: affine-invariant,
/// read on the unit circle, the parabola or the unit hyperbola.
pub(crate) fn segment_share(w: f64) -> f64 {
    if (w - 1.0).abs() < 1e-12 {
        2.0 / 3.0
    } else if w < 1.0 {
        let t = w.acos();
        (t - t.sin() * t.cos()) * t.cos() / t.sin().powi(3)
    } else {
        let p = w.acosh();
        (p.sinh() * p.cosh() - p) * p.cosh() / p.sinh().powi(3)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::{PlanarSection, RoundSpec};
    use crate::asym::RoundProfile;
    use crate::law::{Profile, RadiusLaw};
    use crate::support::Seat;
    use ogeom_core::Tolerances;
    use ogeom_geom::Curve3d as _;
    use ogeom_math::{Point, Vector};
    use ogeom_topo::Shape;

    const T: Tolerances = Tolerances::millimetres();

    /// A box's top front edge along x from the origin, 20 long: the top
    /// (normal +z) and the front (normal -y), convex; or, flipped, the
    /// floor and wall of a re-entrant edge.
    fn seat(convex: bool) -> Seat {
        let mut m = ogeom_topo::Model::new();
        let face: Shape = ogeom_algo::make_vertex(&mut m, Point::ORIGIN).shape;
        let (top, front) = if convex {
            (Vector::new(0.0, 0.0, 1.0), Vector::new(0.0, -1.0, 0.0))
        } else {
            (Vector::new(0.0, 0.0, 1.0), Vector::new(0.0, 1.0, 0.0))
        };
        Seat {
            start: Point::new(0.0, 0.0, 0.0),
            end: Point::new(20.0, 0.0, 0.0),
            along: Vector::new(1.0, 0.0, 0.0),
            normals: [top, front],
            faces: [face.clone(), face],
            convex,
        }
    }

    /// The ball touches each face its radius from the crease on a square
    /// edge, its arc runs a quarter turn between the contacts, every point
    /// of it a radius off the axis, and the cylinder holds the arc.
    #[test]
    fn the_ball_s_section_touches_each_face_its_radius_off_a_square_crease() {
        for convex in [true, false] {
            let section = PlanarSection::new(&seat(convex), &RoundSpec::Circular(2.0), T).unwrap();
            let d = section.distances(5.0);
            assert!(
                (d[0] - 2.0).abs() < 1e-12 && (d[1] - 2.0).abs() < 1e-12,
                "{d:?}"
            );
            let (arc, (a, b)) = section.arc_at(5.0, T).unwrap();
            assert!((b - a - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
            let [first, second] = section.order();
            assert!(
                arc.point_at(a, T)
                    .unwrap()
                    .distance(section.contact(first, 5.0))
                    < 1e-12
            );
            assert!(
                arc.point_at(b, T)
                    .unwrap()
                    .distance(section.contact(second, 5.0))
                    < 1e-12
            );
            let cylinder = section.surface((-3.0, 25.0), T).unwrap();
            let p = arc.point_at((a + b) * 0.5, T).unwrap();
            let found = ogeom_algo::project_on_surface(&cylinder, p, 16, T).unwrap();
            assert!(found.distance < 1e-12);
            let area = section.corner_area(5.0);
            assert!(
                (area - (4.0 - core::f64::consts::PI)).abs() < 1e-12,
                "{area}"
            );
        }
    }

    /// A conic's contacts stand at its distances, its arc leaves each
    /// contact along the face's leg toward the crease, and its corner's
    /// area is the control triangle less the segment.
    #[test]
    fn a_conic_section_meets_each_face_tangent_at_its_distance() {
        for profile in [RoundProfile::Elliptical, RoundProfile::Conic(0.7)] {
            let spec = RoundSpec::Conic {
                d: [3.0, 1.5],
                profile,
            };
            let section = PlanarSection::new(&seat(true), &spec, T).unwrap();
            assert!(section.contact(0, 4.0).distance(Point::new(4.0, 3.0, 0.0)) < 1e-12);
            assert!(section.contact(1, 4.0).distance(Point::new(4.0, 0.0, -1.5)) < 1e-12);
            let (arc, (a, b)) = section.arc_at(4.0, T).unwrap();
            let [first, second] = section.order();
            let start = arc.d1_at(a, T).unwrap();
            let toward = section.crease(4.0) - section.contact(first, 4.0);
            assert!(
                start.cross(toward).magnitude() <= 1e-9 * start.magnitude() * toward.magnitude()
            );
            let end = arc.d1_at(b, T).unwrap();
            let from = section.contact(second, 4.0) - section.crease(4.0);
            assert!(end.cross(from).magnitude() <= 1e-9 * end.magnitude() * from.magnitude());
            let surface = section.surface((0.0, 20.0), T).unwrap();
            let p = arc.point_at(0.5 * (a + b), T).unwrap();
            assert!(
                ogeom_algo::project_on_surface(&surface, p, 16, T)
                    .unwrap()
                    .distance
                    < 1e-10
            );
        }
    }

    /// A varying section's contact stands at its law's distance at each
    /// fraction of the edge, on the contact curve and on the strip's border.
    #[test]
    fn a_varying_section_follows_its_laws_along_the_crease() {
        let laws = [
            Profile::new(
                &RadiusLaw::Points(vec![(0.0, 1.0), (0.4, 2.5), (1.0, 1.5)]),
                false,
                T,
            )
            .unwrap(),
            Profile::new(
                &RadiusLaw::Linear {
                    start: 2.0,
                    end: 1.0,
                },
                false,
                T,
            )
            .unwrap(),
        ];
        let spec = RoundSpec::Varying {
            laws,
            profile: RoundProfile::Elliptical,
        };
        let section = PlanarSection::new(&seat(true), &spec, T).unwrap();
        for (t, want) in [(0.0, 1.0), (0.4, 2.5), (1.0, 1.5)] {
            let p = section.contact(0, 20.0 * t);
            assert!((p.y - want).abs() < 1e-12, "{p:?} at {t}");
            let (curve, _) = section.contact_curve(0, (0.0, 20.0), T).unwrap();
            assert!(curve.point_at(t, T).unwrap().distance(p) < 1e-12);
        }
        let back = section.contact(1, 10.0);
        assert!((back.z + 1.5).abs() < 1e-12, "{back:?}");
        assert!(section.surface((0.0, 20.0), T).is_ok());
        assert!(section.surface((0.0, 25.0), T).is_err());
    }
}
