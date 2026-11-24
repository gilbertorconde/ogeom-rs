//! The sheet's inverse map: from a point of the surface back to the sheet
//! point that lands on it.
//!
//! The point is first taken to its nearest point on the surface. A plane
//! and a cylinder are developed affinely from the origin's chart
//! derivatives, a cone through its development about the apex, and a
//! sphere by the great circle from the origin, all in closed form. Every
//! other surface is reached by shooting: the geodesic from the origin is
//! integrated with its Jacobi field, which gives the derivative of where
//! it lands with respect to the sheet point, and Gauss-Newton steps on the
//! sheet point close the gap in space. A flattened sheet then takes the
//! geodesic polar position back through its polynomial by Newton steps.
//!
//! Where several sheet points land on the point (a cylinder or a cone
//! unrolled past a whole turn), the one nearest the origin is taken: in
//! closed form among every candidate within the reach, and by shooting the
//! one the iteration reaches from the chart's straight step.

use core::f64::consts::PI;

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_geom::Surface as _;
use ogeom_math::{Point, Point2, Vector, Vector2};

use super::{Chart, SurfaceSheet};
use crate::sheet_flat::Flattening;

/// The most Newton or Gauss-Newton steps taken.
const MOST_STEPS: usize = 60;

/// The most halvings of one damped step.
const MOST_HALVINGS: usize = 40;

impl SurfaceSheet {
    /// The sheet point that lands on `p`, a point on or near the surface:
    /// the inverse of [`SurfaceSheet::chart_point`] composed with the
    /// surface. `p` is first taken to its nearest point on the surface,
    /// and the sheet point found lands on that point to the model's
    /// confusion distance. `None` where no sheet point within the sheet's
    /// reach lands there.
    ///
    /// A plane, a cylinder and a cone are unrolled in closed form, and a
    /// sphere's geodesic polar coordinates are read in closed form. Every
    /// other surface is reached by shooting geodesics from the origin, so
    /// the sheet point is found to the accuracy the geodesics are
    /// integrated to, a thousandth of the confusion distance in space. A
    /// flattened sheet takes the geodesic polar position back through its
    /// flattening.
    ///
    /// Where several sheet points land on `p`, as on a cylinder or a cone
    /// whose sheet reaches round a whole turn, the one nearest the origin
    /// is returned. On a surface reached by shooting it is the one the
    /// iteration reaches from the chart's straight step, the nearest where
    /// the geodesics from the origin do not wrap round the surface within
    /// the reach.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where `p` is not a point, its nearest point on the surface is a pole
    /// of the surface's chart (where [`SurfaceSheet::chart_point`] refuses
    /// to land too), it is the point of a sphere opposite the origin, or
    /// where the geodesics from the origin cross before reaching it; where
    /// the shooting or the flattening's inverse does not settle on it.
    pub fn sheet_point(&self, p: Point) -> OgeomResult<Option<Point2>> {
        let tol = self.tol;
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            ogeom_bail!(Construction, "{p:?} is not a point");
        }
        let found = ogeom_geom::project_on_surface(&self.surface, p, 24, tol)?;
        let uv = Point2::new(found.parameters.0, found.parameters.1);
        if self.surface.is_degenerate_at(uv.x, uv.y, tol)? {
            ogeom_bail!(
                Construction,
                "{p:?} stands on a pole of the surface's chart, where no sheet point lands"
            );
        }
        let Some(q) = self.geodesic_back(uv, p)? else {
            return Ok(None);
        };
        let s = match &self.flat {
            None => q,
            Some(flat) => match self.unflattened(flat, q, p)? {
                Some(s) => s,
                None => return Ok(None),
            },
        };
        if s.x.hypot(s.y) > self.reach + tol.confusion() {
            return Ok(None);
        }
        Ok(Some(s))
    }

    /// The geodesic polar position that lands on the chart position `uv`,
    /// `None` past the geodesic reach.
    fn geodesic_back(&self, uv: Point2, p: Point) -> OgeomResult<Option<Point2>> {
        let q = match self.chart {
            Chart::Plane | Chart::Cylinder => self.affine_back(uv)?,
            Chart::Cone {
                sin,
                tan,
                reference,
            } => match self.cone_back(uv, sin, tan, reference)? {
                Some(q) => q,
                None => return Ok(None),
            },
            Chart::Sphere { centre, radius, .. } => self.sphere_back(uv, centre, radius, p)?,
            Chart::Integrated => match self.shot_back(uv, p)? {
                Some(q) => q,
                None => return Ok(None),
            },
        };
        if q.x.hypot(q.y) > self.geodesic_reach + self.tol.confusion() {
            return Ok(None);
        }
        Ok(Some(q))
    }

    /// The sheet point of a tangent vector at the origin.
    fn on_axes(&self, t: Vector) -> Point2 {
        Point2::new(t.dot(self.axes.0), t.dot(self.axes.1))
    }

    /// A plane's or a cylinder's development read back: the chart is
    /// affine on the sheet, so the chart step from the origin is the
    /// sheet point through the origin's chart derivatives.
    fn affine_back(&self, uv: Point2) -> OgeomResult<Point2> {
        let (su, sv) = self.surface.d1_at(self.origin.x, self.origin.y, self.tol)?;
        let step = uv - self.origin;
        let mut q = self.on_axes(su * step.x + sv * step.y);
        // A whole period of a periodic parameter moves the sheet point by
        // a fixed step `w`; the sheet point's length is least at the whole
        // number of them nearest `-(q . w) / (w . w)`.
        let ((ua, ub), (va, vb)) = self.surface.domain();
        for (periodic, along) in [
            (self.surface.is_periodic_u(), su * (ub - ua)),
            (self.surface.is_periodic_v(), sv * (vb - va)),
        ] {
            let w = self.on_axes(along);
            let length = w.x.mul_add(w.x, w.y * w.y);
            if periodic && length > 0.0 {
                let k = (-q.x.mul_add(w.x, q.y * w.y) / length).round();
                q = Point2::new(k.mul_add(w.x, q.x), k.mul_add(w.y, q.y));
            }
        }
        Ok(q)
    }

    /// A cone's development read back: the chart position's developed
    /// point about the apex, at each developed angle within half a turn of
    /// the origin's that lands on it, carried back to a tangent step at the
    /// origin through the development's derivative there. `None` on the
    /// cone's other nappe, which the development from the origin never
    /// reaches.
    fn cone_back(
        &self,
        uv: Point2,
        sin: f64,
        tan: f64,
        reference: f64,
    ) -> OgeomResult<Option<Point2>> {
        let (u0, v0) = (self.origin.x, self.origin.y);
        let (su, sv) = self.surface.d1_at(u0, v0, self.tol)?;
        let rho0 = tan.mul_add(v0, reference) / sin;
        let rho = tan.mul_add(uv.y, reference) / sin;
        if rho * rho0 <= 0.0 {
            return Ok(None);
        }
        let phi = u0 * sin;
        let (s, c) = phi.sin_cos();
        let start = Vector2::new(rho0 * c, rho0 * s);
        let du = Vector2::new(-rho0 * sin * s, rho0 * sin * c);
        let dv = Vector2::new(c, s) * (tan / sin);
        let det = du.x.mul_add(dv.y, -du.y * dv.x);
        let ((ua, ub), _) = self.surface.domain();
        let period = if self.surface.is_periodic_u() {
            ub - ua
        } else {
            0.0
        };
        // The developed angle within half a turn of the origin's, at each
        // whole number of periods of the chart's angle.
        let offset = uv.x - u0;
        let mut k = if period > 0.0 {
            ((-PI / sin - offset) / period).floor()
        } else {
            0.0
        };
        let mut best: Option<Point2> = None;
        loop {
            let turn = k.mul_add(period, offset) * sin;
            k += 1.0;
            if turn > PI {
                break;
            }
            if turn <= -PI {
                if period > 0.0 {
                    continue;
                }
                break;
            }
            let (ts, tc) = (phi + turn).sin_cos();
            let d = Vector2::new(rho * tc, rho * ts) - start;
            let a = d.x.mul_add(dv.y, -d.y * dv.x) / det;
            let b = du.x.mul_add(d.y, -du.y * d.x) / det;
            let q = self.on_axes(su * a + sv * b);
            if best.is_none_or(|w| q.x.hypot(q.y) < w.x.hypot(w.y)) {
                best = Some(q);
            }
            if period <= 0.0 {
                break;
            }
        }
        Ok(best)
    }

    /// A sphere's geodesic polar coordinates of the chart position `uv`:
    /// its angle from the origin about the centre, as a length, along the
    /// tangent at the origin of the great circle through it.
    fn sphere_back(&self, uv: Point2, centre: Point, radius: f64, p: Point) -> OgeomResult<Point2> {
        let tol = self.tol;
        let start = self.surface.point_at(self.origin.x, self.origin.y, tol)?;
        let foot = self.surface.point_at(uv.x, uv.y, tol)?;
        let out = (start - centre) * (1.0 / radius);
        let to = (foot - centre) * (1.0 / radius);
        let across = to - out * out.dot(to);
        let angle = across.magnitude().atan2(out.dot(to));
        if radius * (PI - angle) <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "{p:?} is the point of the sphere opposite the sheet's origin, where \
                 every geodesic from the origin meets"
            );
        }
        let width = across.magnitude();
        if width == 0.0 {
            return Ok(Point2::new(0.0, 0.0));
        }
        Ok(self.on_axes(across * (radius * angle / width)))
    }

    /// Where the geodesic from the origin to the sheet point `s` lands in
    /// space, and the derivatives of that point with respect to `s`'s two
    /// coordinates, from the geodesic's unit tangent and its Jacobi field
    /// there; with `J(r) / r`.
    fn shot(&self, s: Point2, accuracy: f64) -> OgeomResult<(Point, [Vector; 2], f64)> {
        let tol = self.tol;
        let r = s.x.hypot(s.y);
        if r == 0.0 {
            let start = self.surface.point_at(self.origin.x, self.origin.y, tol)?;
            return Ok((start, [self.axes.0, self.axes.1], 1.0));
        }
        let along = self.axes.0 * s.x + self.axes.1 * s.y;
        let state = self.integrated_state(along, r, accuracy, s)?;
        let (su, sv) = self.surface.d1_at(state[0], state[1], tol)?;
        let landed = self.surface.point_at(state[0], state[1], tol)?;
        let tangent = su * state[2] + sv * state[3];
        // The side of the chart's normal the sheet's axes turn about: its
        // `+y` is that normal turned against `+x`.
        let (su0, sv0) = self.surface.d1_at(self.origin.x, self.origin.y, tol)?;
        let side = self.axes.1.dot(su0.cross(sv0).cross(self.axes.0)).signum();
        let normal = su.cross(sv);
        let width = normal.magnitude();
        if width.is_nan() || width <= 0.0 {
            ogeom_bail!(
                Construction,
                "the geodesic to the sheet point {s:?} lands on a pole of the surface's chart"
            );
        }
        let across = normal.cross(tangent) * (side / width);
        let stretch = state[4] / r;
        let (c, n) = (s.x / r, s.y / r);
        Ok((
            landed,
            [
                tangent * c - across * (n * stretch),
                tangent * n + across * (c * stretch),
            ],
            stretch,
        ))
    }

    /// The geodesic polar position that lands nearest `p`, whose
    /// projection on the surface is the chart position `uv`, on a surface
    /// whose geodesics are integrated: Gauss-Newton steps on the sheet
    /// point, from the chart's straight step to `uv`, each halved until it
    /// shortens the miss along the surface. `None` where the surface's
    /// distance in space from the origin already passes the reach.
    fn shot_back(&self, uv: Point2, p: Point) -> OgeomResult<Option<Point2>> {
        let tol = self.tol;
        let accuracy = tol.confusion() * 1e-3;
        let start = self.surface.point_at(self.origin.x, self.origin.y, tol)?;
        let foot = self.surface.point_at(uv.x, uv.y, tol)?;
        // A geodesic is no shorter than the chord.
        if foot.distance(start) > self.geodesic_reach + tol.confusion() {
            return Ok(None);
        }
        let (su, sv) = self.surface.d1_at(self.origin.x, self.origin.y, tol)?;
        let periods = [
            self.surface.is_periodic_u().then(|| {
                let ((ua, ub), _) = self.surface.domain();
                ub - ua
            }),
            self.surface.is_periodic_v().then(|| {
                let (_, (va, vb)) = self.surface.domain();
                vb - va
            }),
        ];
        let near = super::nearest(uv, self.origin, periods);
        let step = near - self.origin;
        let mut s = self.on_axes(su * step.x + sv * step.y);
        // The straight step may run off a bounded surface: it is shortened
        // until a geodesic follows it.
        let mut now = None;
        for _ in 0..MOST_HALVINGS {
            if let Ok(shot) = self.shot(s, accuracy) {
                now = Some(shot);
                break;
            }
            s = Point2::new(s.x / 2.0, s.y / 2.0);
        }
        let Some(mut now) = now else {
            ogeom_bail!(
                Construction,
                "no geodesic from the sheet's origin toward {p:?} could be followed"
            );
        };
        // The miss measured along the surface: `p`'s offset from where the
        // geodesic lands, its part across the surface taken off, which
        // vanishes at the landed point nearest `p`. Aiming at `p` itself,
        // not at the foot its projection found, keeps the projection's own
        // tolerance out of the result.
        let along = |shot: &(Point, [Vector; 2], f64)| {
            let miss = p - shot.0;
            let normal = shot.1[0].cross(shot.1[1]);
            let area = normal.dot(normal);
            if area > 0.0 {
                miss - normal * (normal.dot(miss) / area)
            } else {
                miss
            }
        };
        let mut gap = along(&now).magnitude();
        for _ in 0..MOST_STEPS {
            if gap <= accuracy {
                break;
            }
            let [a, b] = now.1;
            let miss = p - now.0;
            let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
            let det = aa.mul_add(bb, -ab * ab);
            if det.is_nan() || det <= 0.0 {
                ogeom_bail!(
                    Construction,
                    "the geodesics from the sheet's origin cross before reaching {p:?}"
                );
            }
            let (ra, rb) = (a.dot(miss), b.dot(miss));
            let mut delta = Vector2::new(
                bb.mul_add(ra, -ab * rb) / det,
                aa.mul_add(rb, -ab * ra) / det,
            );
            let mut moved = false;
            for _ in 0..MOST_HALVINGS {
                let next = s + delta;
                if let Ok(shot) = self.shot(next, accuracy) {
                    let g = along(&shot).magnitude();
                    if g < gap {
                        s = next;
                        now = shot;
                        gap = g;
                        moved = true;
                        break;
                    }
                }
                delta *= 0.5;
            }
            if !moved {
                break;
            }
        }
        if gap > tol.confusion() {
            ogeom_bail!(
                Construction,
                "the geodesics from the sheet's origin do not settle on {p:?}: the nearest \
                 misses it by {gap:.3e} along the surface"
            );
        }
        if now.2.is_nan() || now.2 <= 0.0 {
            ogeom_bail!(
                Construction,
                "the geodesics from the sheet's origin cross before reaching {p:?}: the \
                 sheet folds there"
            );
        }
        Ok(Some(s))
    }

    /// The flattened sheet point the flattening carries to the geodesic
    /// polar position `q`: Newton steps on its polynomial from `q` itself,
    /// each halved until it brings the image nearer. `None` where the
    /// iteration runs past the sheet's reach, where the polynomial stands
    /// for no flattening.
    fn unflattened(&self, flat: &Flattening, q: Point2, p: Point) -> OgeomResult<Option<Point2>> {
        let target = q.x.hypot(q.y).max(1.0) * f64::EPSILON * 16.0;
        let mut s = q;
        let (mut image, mut d) = flat.at(s);
        let mut gap = (q - image).magnitude();
        for _ in 0..MOST_STEPS {
            if gap <= target {
                break;
            }
            let det = d[0][0].mul_add(d[1][1], -d[0][1] * d[1][0]);
            if det.is_nan() || det <= 0.0 {
                ogeom_bail!(
                    Construction,
                    "the flattened sheet folds over on its way back to {p:?}"
                );
            }
            let miss = q - image;
            let mut delta = Vector2::new(
                d[1][1].mul_add(miss.x, -d[0][1] * miss.y) / det,
                d[0][0].mul_add(miss.y, -d[1][0] * miss.x) / det,
            );
            let mut moved = false;
            for _ in 0..MOST_HALVINGS {
                let next = s + delta;
                let (i, dn) = flat.at(next);
                let g = (q - i).magnitude();
                if g < gap {
                    (s, image, d, gap) = (next, i, dn, g);
                    moved = true;
                    break;
                }
                delta *= 0.5;
            }
            if !moved || s.x.hypot(s.y) > 2.0 * self.reach {
                break;
            }
        }
        if s.x.hypot(s.y) > self.reach + self.tol.confusion() {
            return Ok(None);
        }
        if gap > self.tol.confusion() * 1e-3 {
            ogeom_bail!(
                Construction,
                "the flattened sheet's map back does not settle on {p:?}: the nearest \
                 lands {gap:.3e} from it"
            );
        }
        Ok(Some(s))
    }
}
