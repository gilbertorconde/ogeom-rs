//! A flat sheet laid onto a face: the chart from the plane to the face's
//! surface that carries a sketch onto it.
//!
//! The sheet is laid by geodesic polar coordinates, the exponential map at
//! its origin: the sheet point at distance `r` and angle `t` from the
//! origin lands where the geodesic leaving the origin's point at angle `t`
//! from the sheet's `+x` is after length `r`. Lengths along every straight
//! line through the origin are kept; across those lines a length is
//! stretched by `J(r) / r`, `J` the Jacobi field along the geodesic
//! (`J'' = -K J`, `J(0) = 0`, `J'(0) = 1`, `K` the Gaussian curvature), and
//! so is an area, since the radial length is kept. Where the surface is
//! flat (a plane, a cylinder, a cone) the map is the surface's development
//! and stretches nothing.
//!
//! A plane, a cylinder, a cone and a sphere are laid in closed form. Every
//! other surface has its geodesics integrated in its own chart, the
//! geodesic equation and the Jacobi equation together, by classical
//! Runge-Kutta steps whose error is estimated by step doubling and held,
//! per unit of length, to the share of the accuracy asked that the whole
//! length leaves it. A chart position is reached continuously along the
//! sheet's straight line from the origin, so it runs on past a periodic
//! seam unwrapped; a carried curve's chart curve is continued along the
//! curve itself instead, so one going round a pole of the chart winds a
//! whole period, as a parallel does.

use core::f64::consts::PI;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail, ogeom_err};
use ogeom_geom::{
    BSpline2d, Curve, Curve2d as _, CurveOnSurface, PlanarCurve, Surface as _, SurfaceGeometry,
    Transformable as _,
};
use ogeom_math::{Direction, KnotVector, Point, Point2, Vector, Vector2};
use ogeom_topo::{Location, Model, Orientation, Shape, ShapeType};

use crate::Built;
use crate::build::{attach_pcurve, make_edge};

/// How the flat sheet is carried away from its origin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SheetMethod {
    /// Geodesic polar coordinates: the sheet point at distance `r` and
    /// angle `t` from the origin lands where the geodesic leaving the
    /// origin at angle `t` from the sheet's `+x` is after length `r` (the
    /// exponential map). Lengths along lines through the origin are kept.
    Geodesic,
}

/// How a flat sheet is laid on a face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SheetPlacement {
    /// The point of the face the sheet's origin lands on.
    pub at: Point,
    /// The direction the sheet's `+x` runs along at `at`, its part in the
    /// tangent plane taken.
    pub x_direction: Direction,
    /// How the sheet is carried away from `at`.
    pub method: SheetMethod,
}

/// How far a laid sheet stretches, over the sheet points it has carried.
///
/// Each is a pair: the least and the most ratio of a length (or an area)
/// on the surface to the same on the sheet. Geodesic polar coordinates
/// keep the radial length, so the length's range always holds 1, and the
/// area's ratio is the across-ray stretch alone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SheetDistortion {
    /// The least and the most length ratio.
    pub length: (f64, f64),
    /// The least and the most area ratio.
    pub area: (f64, f64),
}

/// A flat sheet laid on a face's surface, as placed.
#[derive(Debug, Clone)]
pub struct SurfaceSheet {
    surface: SurfaceGeometry,
    /// The face's own surface, where the face stands unplaced, for an edge
    /// to carry its pcurve on.
    own: Option<ogeom_topo::SurfaceId>,
    chart: Chart,
    origin: Point2,
    /// The sheet's `+x` and `+y` at the origin, unit tangents in space, `y`
    /// the face's outward normal turned against `x`.
    axes: (Vector, Vector),
    reach: f64,
    tol: Tolerances,
    /// The least and the most across-ray stretch carried so far.
    stretch: core::cell::Cell<(f64, f64)>,
}

/// How a surface's chart is reached from the sheet.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Chart {
    /// A plane's chart, isometric to the sheet.
    Plane,
    /// A cylinder's development.
    Cylinder,
    /// A cone's development: slant `rho(v) = radius_at(v) / sin(half)`
    /// from the apex, chart angle `u` developed to `u sin(half)`.
    Cone { sin: f64, tan: f64, reference: f64 },
    /// A sphere's exponential map in closed form.
    Sphere {
        frame: ogeom_math::Frame,
        centre: Point,
        radius: f64,
    },
    /// Geodesics integrated in the surface's chart.
    Integrated,
}

impl SurfaceSheet {
    /// The sheet laid on `face` of `model`, its surface as placed.
    ///
    /// `reach` is the radius of the sheet region that will be asked for:
    /// a sheet point farther from the origin is refused.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where `face` is not a face, `at` stands off its surface by more
    /// than the face's tolerance, `x_direction` is along the normal there,
    /// `at` is a pole of the surface's chart, `reach` is not a positive
    /// length, or a sphere's sheet reaches its origin's antipode.
    pub fn new(
        model: &Model,
        face: &Shape,
        placement: &SheetPlacement,
        reach: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        match placement.method {
            SheetMethod::Geodesic => {}
        }
        if model.kind_of(face)? != ShapeType::Face {
            ogeom_bail!(Construction, "a sheet is laid on a face");
        }
        if !(reach > 0.0 && reach.is_finite()) {
            ogeom_bail!(Construction, "a sheet's reach of {reach} is not a length");
        }
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "surface is not in this model");
        };
        let unplaced = face.location().is_identity() && data.location.is_identity();
        let surface = if unplaced {
            stored.clone()
        } else {
            let placement =
                face.transform(model.datums())? * data.location.composed(model.datums())?;
            stored.transformed(&placement, tol)?
        };
        let found = ogeom_geom::project_on_surface(&surface, placement.at, 24, tol)?;
        let off = data.tolerance.get().max(tol.confusion());
        if found.distance > off {
            ogeom_bail!(
                Construction,
                "the sheet's origin {:?} stands {:.3e} off the face's surface",
                placement.at,
                found.distance
            );
        }
        let origin = Point2::new(found.parameters.0, found.parameters.1);
        let (su, sv) = surface.d1_at(origin.x, origin.y, tol)?;
        if surface.is_degenerate_at(origin.x, origin.y, tol)? {
            ogeom_bail!(
                Construction,
                "the sheet's origin {:?} is a pole of the surface's chart",
                placement.at
            );
        }
        let mut normal = Direction::new(su.cross(sv), tol)?.vector();
        if face.orientation() == Orientation::Reversed {
            normal = -normal;
        }
        let along = placement.x_direction.vector();
        let x = along - normal * along.dot(normal);
        let Ok(x) = Direction::new(x, tol) else {
            ogeom_bail!(
                Construction,
                "the sheet's x direction runs along the face's normal at its origin"
            );
        };
        let x = x.vector();
        let y = normal.cross(x);
        let chart = match &surface {
            SurfaceGeometry::Plane(_) => Chart::Plane,
            SurfaceGeometry::Cylinder(_) => Chart::Cylinder,
            SurfaceGeometry::Cone(c) => {
                let cone = c.cone();
                let (sin, cos) = cone.half_angle().sin_cos();
                Chart::Cone {
                    sin,
                    tan: sin / cos,
                    reference: cone.reference_radius(),
                }
            }
            SurfaceGeometry::Sphere(s) => {
                let sphere = s.sphere();
                if reach >= PI * sphere.radius() {
                    ogeom_bail!(
                        Construction,
                        "a sheet of reach {reach} on a sphere of radius {} reaches the \
                         point opposite its origin, where every line through the origin \
                         meets",
                        sphere.radius()
                    );
                }
                Chart::Sphere {
                    frame: sphere.frame(),
                    centre: sphere.centre(),
                    radius: sphere.radius(),
                }
            }
            _ => Chart::Integrated,
        };
        Ok(Self {
            surface,
            own: unplaced.then_some(data.surface),
            chart,
            origin,
            axes: (x, y),
            reach,
            tol,
            stretch: core::cell::Cell::new((1.0, 1.0)),
        })
    }

    /// The surface the sheet is laid on, as placed: the one the chart
    /// positions are on.
    #[must_use]
    pub const fn surface(&self) -> &SurfaceGeometry {
        &self.surface
    }

    /// The `(u, v)` of the surface a sheet point lands on, to the model's
    /// confusion distance in space.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where the point lies past the sheet's reach, lands on a pole of the
    /// chart, runs off a bounded surface's domain, or where the geodesics
    /// from the origin cross before reaching it.
    pub fn chart_point(&self, p: Point2) -> OgeomResult<Point2> {
        self.landed(p, self.tol.confusion())
    }

    /// A planar curve over `range` of its parameter carried onto the
    /// surface: its chart curve, at the same parameter, whose lift onto
    /// the surface stands within `tolerance` of the carried curve, ready
    /// for [`CurveOnSurface::new`] with [`SurfaceSheet::surface`].
    ///
    /// A line on a plane or a cylinder is carried exactly; anything else
    /// is fitted, half of `tolerance` left to the fit and half to the
    /// geodesics, and the fit is measured in space between its samples.
    /// The chart curve is continued along the curve across every period
    /// of the chart, so a curve going round a pole of the chart (a circle
    /// about a sphere's pole) is one chart curve whose ends stand a whole
    /// turn apart, as a parallel's do.
    ///
    /// # Errors
    ///
    /// As [`SurfaceSheet::chart_point`] at any point of the curve;
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where `tolerance` is not a length, where the curve runs through a
    /// pole of the chart, so its chart curve cannot be continued along
    /// it, or where the fit cannot be held to `tolerance`.
    pub fn chart_curve(
        &self,
        curve: &PlanarCurve,
        range: (f64, f64),
        tolerance: f64,
    ) -> OgeomResult<PlanarCurve> {
        let tol = self.tol;
        if !(tolerance > 0.0 && tolerance.is_finite()) {
            ogeom_bail!(Construction, "a tolerance of {tolerance} is not a length");
        }
        if range.0.is_nan() || range.1.is_nan() || range.1 <= range.0 {
            ogeom_bail!(Construction, "a curve is carried over an increasing range");
        }
        if let PlanarCurve::Line(_) = curve
            && matches!(self.chart, Chart::Plane | Chart::Cylinder)
        {
            // The chart is affine on the sheet: a line's image is the line
            // through its ends' images, at the same parameter.
            let ends = [curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?];
            let mut control = Vec::with_capacity(2);
            for end in ends {
                control.push(self.chart_point(end)?);
            }
            let knots = KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?;
            return Ok(BSpline2d::new(knots, control, tol)?.into());
        }
        let half = tolerance / 2.0;
        let surface = &self.surface;
        let lift = Unwrapped::along(self, curve, range, half)?;
        let (fitted, image) = ogeom_geom::fit::fit_trace_sampled(
            |t| {
                let uv = lift.at(self, curve, t, half)?;
                Ok((surface.point_at(uv.x, uv.y, tol)?, uv))
            },
            |uv| surface.point_at(uv.x, uv.y, tol),
            &lift.ts,
            3,
            half,
            tol,
        )?;
        if !fitted.met {
            ogeom_bail!(
                Construction,
                "the carried curve's chart curve misses it by {:.3e}, over the {half:.3e} \
                 asked",
                fitted.error
            );
        }
        Ok(image.into())
    }

    /// An edge on the surface along a planar curve carried by
    /// [`SurfaceSheet::chart_curve`]: its curve the chart curve on the
    /// surface, and, on a face standing unplaced, the chart curve its
    /// pcurve on the face's surface.
    ///
    /// # Errors
    ///
    /// As [`SurfaceSheet::chart_curve`] and
    /// [`make_edge`].
    pub fn chart_edge(
        &self,
        model: &mut Model,
        curve: &PlanarCurve,
        range: (f64, f64),
        tolerance: f64,
    ) -> OgeomResult<Built> {
        let pcurve = self.chart_curve(curve, range, tolerance)?;
        let on = Curve::OnSurface(Box::new(CurveOnSurface::new(
            pcurve.clone(),
            self.surface.clone(),
        )));
        let built = make_edge(model, on, range, self.tol)?;
        if let Some(surface) = self.own {
            attach_pcurve(
                model,
                &built.shape,
                pcurve,
                surface,
                Location::identity(),
                range,
            )?;
        }
        Ok(built)
    }

    /// How far the sheet has stretched, over every sheet point carried so
    /// far.
    #[must_use]
    pub fn distortion(&self) -> SheetDistortion {
        let (least, most) = self.stretch.get();
        SheetDistortion {
            length: (least.min(1.0), most.max(1.0)),
            area: (least, most),
        }
    }

    /// The chart position a sheet point lands on, to `accuracy` in space,
    /// its across-ray stretch recorded.
    fn landed(&self, p: Point2, accuracy: f64) -> OgeomResult<Point2> {
        let r = p.x.hypot(p.y);
        if !r.is_finite() {
            ogeom_bail!(Construction, "the sheet point {p:?} is not a point");
        }
        if r > self.reach + self.tol.confusion() {
            ogeom_bail!(
                Construction,
                "the sheet point {p:?} lies past the sheet's reach of {}",
                self.reach
            );
        }
        let along = self.axes.0 * p.x + self.axes.1 * p.y;
        let (uv, stretch) = match self.chart {
            Chart::Plane | Chart::Cylinder => (self.developed(along)?, 1.0),
            Chart::Cone {
                sin,
                tan,
                reference,
            } => (self.developed_cone(along, sin, tan, reference, p)?, 1.0),
            Chart::Sphere {
                frame,
                centre,
                radius,
            } => self.on_sphere(along, r, &frame, centre, radius, p)?,
            Chart::Integrated => self.integrated(along, r, accuracy, p)?,
        };
        if stretch.is_nan() || stretch <= 0.0 {
            ogeom_bail!(
                Construction,
                "the geodesics from the sheet's origin cross before the sheet point \
                 {p:?}: the sheet folds there"
            );
        }
        let (least, most) = self.stretch.get();
        self.stretch.set((least.min(stretch), most.max(stretch)));
        Ok(uv)
    }

    /// The chart step from the origin along the tangent vector `t` there,
    /// read in the origin's chart derivatives: where the chart is affine
    /// on the development (a plane, a cylinder), the position itself.
    fn developed(&self, t: Vector) -> OgeomResult<Point2> {
        let (su, sv) = self.surface.d1_at(self.origin.x, self.origin.y, self.tol)?;
        let step = in_chart(su, sv, t)?;
        Ok(self.origin + step)
    }

    /// A cone's development: the origin's developed point stepped along
    /// `t`, read back as slant and angle from the apex.
    fn developed_cone(
        &self,
        t: Vector,
        sin: f64,
        tan: f64,
        reference: f64,
        p: Point2,
    ) -> OgeomResult<Point2> {
        let (u0, v0) = (self.origin.x, self.origin.y);
        let (su, sv) = self.surface.d1_at(u0, v0, self.tol)?;
        let step = in_chart(su, sv, t)?;
        // Slant `rho = radius / sin`, signed, and developed angle `u sin`.
        let rho = tan.mul_add(v0, reference) / sin;
        let phi = u0 * sin;
        let (s, c) = phi.sin_cos();
        let start = Vector2::new(rho * c, rho * s);
        let du = Vector2::new(-rho * sin * s, rho * sin * c);
        let dv = Vector2::new(c, s) * (tan / sin);
        let d = start + du * step.x + dv * step.y;
        let slant = d.magnitude();
        if slant <= self.tol.confusion() {
            ogeom_bail!(
                Construction,
                "the sheet point {p:?} lands on the cone's apex, a pole of its chart"
            );
        }
        // The developed point read on the origin's own side of the apex.
        let sign = rho.signum();
        let turn = wrapped((d.y * sign).atan2(d.x * sign) - phi);
        let rho = slant * sign;
        Ok(Point2::new(u0 + turn / sin, (rho * sin - reference) / tan))
    }

    /// A sphere's exponential map: the great circle leaving the origin
    /// along `t`, after length `r`, read in the sphere's chart with its
    /// longitude continued from the origin's; the stretch is
    /// `sin(r / R) / (r / R)`.
    fn on_sphere(
        &self,
        t: Vector,
        r: f64,
        frame: &ogeom_math::Frame,
        centre: Point,
        radius: f64,
        p: Point2,
    ) -> OgeomResult<(Point2, f64)> {
        let tol = self.tol;
        let start = self.surface.point_at(self.origin.x, self.origin.y, tol)?;
        let out = (start - centre) * (1.0 / radius);
        let angle = r / radius;
        let landed = if r > 0.0 {
            let (s, c) = angle.sin_cos();
            centre + (out * c + t * (s / r)) * radius
        } else {
            start
        };
        let local = landed - centre;
        let (x, y, z) = (
            local.dot(frame.x().vector()),
            local.dot(frame.y().vector()),
            local.dot(frame.z().vector()),
        );
        if x.hypot(y) <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the sheet point {p:?} lands on a pole of the sphere's chart, at {landed:?}"
            );
        }
        let latitude = (z / radius).clamp(-1.0, 1.0).asin();
        let longitude = self.origin.x + wrapped(y.atan2(x) - self.origin.x);
        let stretch = if angle > 0.0 {
            angle.sin() / angle
        } else {
            1.0
        };
        Ok((Point2::new(longitude, latitude), stretch))
    }

    /// The geodesic leaving the origin along the unit-length direction
    /// `t / r`, after length `r`, integrated in the chart with its Jacobi
    /// field, to `accuracy` in space: its chart position and `J(r) / r`.
    fn integrated(
        &self,
        t: Vector,
        r: f64,
        accuracy: f64,
        p: Point2,
    ) -> OgeomResult<(Point2, f64)> {
        if r == 0.0 {
            return Ok((self.origin, 1.0));
        }
        let tol = self.tol;
        let (su, sv) = self.surface.d1_at(self.origin.x, self.origin.y, tol)?;
        let velocity = in_chart(su, sv, t * (1.0 / r))?;
        let mut state = [
            self.origin.x,
            self.origin.y,
            velocity.x,
            velocity.y,
            0.0,
            1.0,
        ];
        let mut s = 0.0;
        let mut h = r;
        while s < r {
            h = h.min(r - s);
            if h <= tol.confusion() && s + h < r {
                ogeom_bail!(
                    Construction,
                    "the geodesic to the sheet point {p:?} cannot be followed past \
                     chart position ({}, {}): the chart degenerates there",
                    state[0],
                    state[1]
                );
            }
            let whole = self.rk4(&state, h, p)?;
            let first = self.rk4(&state, h / 2.0, p)?;
            let both = self.rk4(&first, h / 2.0, p)?;
            let (du, dv) = self.surface.d1_at(both[0], both[1], tol)?;
            // Step doubling: the two halves' error is a fifteenth of their
            // difference from the whole step, the classical fourth order's
            // Richardson estimate.
            let miss = (du * (both[0] - whole[0]) + dv * (both[1] - whole[1]))
                .magnitude()
                .max((both[4] - whole[4]).abs())
                / 15.0;
            let allowed = accuracy * h / r;
            // A last remnant no longer than the confusion distance is taken
            // as it is: it cannot move the point by more than its length.
            if miss <= allowed || h <= tol.confusion() {
                state = both;
                // The extrapolated fifth-order value.
                for (k, value) in state.iter_mut().enumerate() {
                    *value += (both[k] - whole[k]) / 15.0;
                }
                s += h;
            }
            // The step that would just meet the allowance, its error
            // growing with its fifth power, its allowance with its first;
            // the clamp and the margin on a refused step only pace the
            // search, the allowance alone decides what is kept.
            let grow = if miss > 0.0 {
                (allowed / miss).powf(0.25)
            } else {
                2.0
            };
            h *= grow.clamp(0.25, 2.0) * if miss > allowed { 0.9 } else { 1.0 };
        }
        Ok((Point2::new(state[0], state[1]), state[4] / r))
    }

    /// One classical Runge-Kutta step of length `h`.
    fn rk4(&self, y: &[f64; 6], h: f64, p: Point2) -> OgeomResult<[f64; 6]> {
        let k1 = self.slope(y, p)?;
        let k2 = self.slope(&ahead(y, &k1, h / 2.0), p)?;
        let k3 = self.slope(&ahead(y, &k2, h / 2.0), p)?;
        let k4 = self.slope(&ahead(y, &k3, h), p)?;
        let mut out = *y;
        for k in 0..6 {
            out[k] += h / 6.0 * (k1[k] + 2.0 * (k2[k] + k3[k]) + k4[k]);
        }
        Ok(out)
    }

    /// The geodesic and Jacobi equations at the state
    /// `(u, v, u', v', J, J')`: the chart acceleration that keeps the
    /// curve's acceleration along the normal, and `J'' = -K J`.
    fn slope(&self, y: &[f64; 6], p: Point2) -> OgeomResult<[f64; 6]> {
        let tol = self.tol;
        let (u, v) = (y[0], y[1]);
        let ((ua, ub), (va, vb)) = self.surface.domain();
        let slack = tol.parametric();
        let off_u = !self.surface.is_periodic_u() && (u < ua - slack || u > ub + slack);
        let off_v = !self.surface.is_periodic_v() && (v < va - slack || v > vb + slack);
        if off_u || off_v {
            ogeom_bail!(
                Construction,
                "the sheet point {p:?} runs off the surface's domain at chart \
                 position ({u}, {v})"
            );
        }
        let jet = self.surface.jet_at(u, v, tol)?;
        let (e, f, g) = (jet.du.dot(jet.du), jet.du.dot(jet.dv), jet.dv.dot(jet.dv));
        let det = e.mul_add(g, -f * f);
        let area = jet.du.cross(jet.dv);
        if det.is_nan() || det <= 0.0 || area.magnitude() <= tol.angular() * e.max(g) {
            return Err(ogeom_err!(
                Construction,
                "the geodesic to the sheet point {p:?} reaches a pole of the surface's \
                 chart at ({u}, {v})"
            ));
        }
        let (a, b) = (y[2], y[3]);
        let bend = jet.d2u * (a * a) + jet.duv * (2.0 * a * b) + jet.d2v * (b * b);
        let (ru, rv) = (-jet.du.dot(bend), -jet.dv.dot(bend));
        let accel_u = g.mul_add(ru, -f * rv) / det;
        let accel_v = e.mul_add(rv, -f * ru) / det;
        let normal = area * (1.0 / area.magnitude());
        let (l, m, n) = (
            jet.d2u.dot(normal),
            jet.duv.dot(normal),
            jet.d2v.dot(normal),
        );
        let gaussian = l.mul_add(n, -m * m) / det;
        Ok([a, b, accel_u, accel_v, y[5], -gaussian * y[4]])
    }
}

/// A carried curve's chart positions continued along it: each periodic
/// chart coordinate taken at the period nearest its neighbours', read
/// from a table of positions dense enough that no step between two of
/// them turns by as much as a quarter of a period. A position away from
/// the table is its ray's position moved by the whole periods that put it
/// nearest the table's line between the samples either side.
struct Unwrapped {
    periods: [Option<f64>; 2],
    ts: Vec<f64>,
    uvs: Vec<Point2>,
}

impl Unwrapped {
    /// The table along `curve` over `range`, positions to `accuracy`.
    fn along(
        sheet: &SurfaceSheet,
        curve: &PlanarCurve,
        range: (f64, f64),
        accuracy: f64,
    ) -> OgeomResult<Self> {
        let surface = &sheet.surface;
        let ((ua, ub), (va, vb)) = surface.domain();
        let periods = [
            surface.is_periodic_u().then_some(ub - ua),
            surface.is_periodic_v().then_some(vb - va),
        ];
        let mut ts: Vec<f64> = (0..=64)
            .map(|k| (range.1 - range.0).mul_add(f64::from(k) / 64.0, range.0))
            .collect();
        let mut seen: ogeom_core::FastMap<u64, Point2> = ogeom_core::FastMap::default();
        loop {
            let mut uvs: Vec<Point2> = Vec::with_capacity(ts.len());
            for &t in &ts {
                let ray = match seen.get(&t.to_bits()) {
                    Some(uv) => *uv,
                    None => {
                        let uv = sheet.landed(curve.point_at(t, sheet.tol)?, accuracy)?;
                        seen.insert(t.to_bits(), uv);
                        uv
                    }
                };
                uvs.push(match uvs.last() {
                    Some(previous) => nearest(ray, *previous, periods),
                    None => ray,
                });
            }
            let mut next = Vec::with_capacity(ts.len() * 2);
            for (k, pair) in ts.windows(2).enumerate() {
                next.push(pair[0]);
                let step = uvs[k + 1] - uvs[k];
                let wide = [step.x, step.y]
                    .iter()
                    .zip(periods)
                    .any(|(d, p)| p.is_some_and(|p| d.abs() > p / 4.0));
                if wide {
                    let middle = f64::midpoint(pair[0], pair[1]);
                    if middle <= pair[0] || middle >= pair[1] {
                        ogeom_bail!(
                            Construction,
                            "the carried curve runs through a pole of the surface's chart \
                             near its parameter {middle}: its chart curve cannot be \
                             continued along it"
                        );
                    }
                    next.push(middle);
                }
            }
            next.push(ts[ts.len() - 1]);
            if next.len() == ts.len() {
                return Ok(Self { periods, ts, uvs });
            }
            if next.len() > ogeom_geom::fit::SAMPLED_SPANS + 1 {
                ogeom_bail!(
                    Construction,
                    "the carried curve runs through a pole of the surface's chart: its \
                     chart curve turns too fast to be continued along it"
                );
            }
            ts = next;
        }
    }

    /// The continued chart position at `t`.
    fn at(
        &self,
        sheet: &SurfaceSheet,
        curve: &PlanarCurve,
        t: f64,
        accuracy: f64,
    ) -> OgeomResult<Point2> {
        let ray = sheet.landed(curve.point_at(t, sheet.tol)?, accuracy)?;
        let last = self.ts.len() - 1;
        let k = self.ts.partition_point(|s| *s <= t).clamp(1, last);
        let (t0, t1) = (self.ts[k - 1], self.ts[k]);
        let f = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);
        let guide = self.uvs[k - 1] + (self.uvs[k] - self.uvs[k - 1]) * f;
        let out = nearest(ray, guide, self.periods);
        Ok(out)
    }
}

/// `uv` moved by the whole periods that bring each periodic coordinate
/// nearest `guide`'s.
fn nearest(uv: Point2, guide: Point2, periods: [Option<f64>; 2]) -> Point2 {
    let near = |x: f64, g: f64, p: Option<f64>| match p {
        Some(p) => ((g - x) / p).round().mul_add(p, x),
        None => x,
    };
    Point2::new(
        near(uv.x, guide.x, periods[0]),
        near(uv.y, guide.y, periods[1]),
    )
}

/// The state `y` stepped `h` along the slope `k`.
fn ahead(y: &[f64; 6], k: &[f64; 6], h: f64) -> [f64; 6] {
    let mut out = *y;
    for i in 0..6 {
        out[i] += h * k[i];
    }
    out
}

/// The chart step `(du, dv)` whose image under the chart derivatives
/// `su`, `sv` is the tangent vector `t`: the normal equations of the
/// first fundamental form.
fn in_chart(su: Vector, sv: Vector, t: Vector) -> OgeomResult<Vector2> {
    let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let det = e.mul_add(g, -f * f);
    if det.is_nan() || det <= 0.0 {
        ogeom_bail!(
            Construction,
            "the surface's chart degenerates at the sheet's origin"
        );
    }
    let (a, b) = (t.dot(su), t.dot(sv));
    Ok(Vector2::new(
        g.mul_add(a, -f * b) / det,
        e.mul_add(b, -f * a) / det,
    ))
}

/// An angle brought into `(-pi, pi]`.
fn wrapped(angle: f64) -> f64 {
    let turned = (angle + PI).rem_euclid(2.0 * PI) - PI;
    if turned <= -PI {
        turned + 2.0 * PI
    } else {
        turned
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

    use super::*;
    use ogeom_math::Frame;
    use ogeom_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn integrated_geodesics_on_a_sphere_meet_its_closed_form() {
        // The sphere laid in closed form and the same sphere laid by
        // integrating its geodesics in its chart: where every ray lands,
        // and how far across the rays it stretches, sin(r / R) / (r / R).
        let mut model = Model::new();
        let sphere = crate::make_sphere(&mut model, Frame::WORLD, 10.0, T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &sphere, ShapeType::Face).unwrap()[0].clone();
        // Off the equator, where the chart's metric varies along every ray.
        let at = Point::new(3.0, -9.0, 1.0);
        let placement = SheetPlacement {
            at: Point::from_vector(at.to_vector() * (10.0 / at.to_vector().magnitude())),
            x_direction: Direction::X,
            method: SheetMethod::Geodesic,
        };
        let exact = SurfaceSheet::new(&model, &face, &placement, 9.0, T).unwrap();
        let mut integrated = exact.clone();
        integrated.chart = Chart::Integrated;
        let mut worst: f64 = 0.0;
        for k in 0..24 {
            let angle = f64::from(k) * PI / 12.0;
            for r in [0.5, 2.0, 5.0, 9.0] {
                let p = Point2::new(r * angle.cos(), r * angle.sin());
                let (a, b) = (
                    exact.chart_point(p).unwrap(),
                    integrated.chart_point(p).unwrap(),
                );
                let (pa, pb) = (
                    exact.surface.point_at(a.x, a.y, T).unwrap(),
                    integrated.surface.point_at(b.x, b.y, T).unwrap(),
                );
                worst = worst.max(pa.distance(pb));
            }
        }
        assert!(worst <= 1e-6, "{worst}");
        let (e, i) = (exact.distortion(), integrated.distortion());
        assert!((e.area.0 - 0.9_f64.sin() / 0.9).abs() <= 1e-12, "{e:?}");
        assert!((i.area.0 - e.area.0).abs() <= 1e-6, "{i:?} {e:?}");
    }
}
