//! The sketch: entities over a flat parameter vector, constraints as
//! residual equations.
//!
//! A sketch is a system of equations in disguise. Every point contributes two
//! parameters, every circle one more for its radius; lines and arcs borrow
//! their points rather than owning coordinates, so coincidence never has to
//! be simulated by equations when the model can simply share. Each constraint
//! contributes one or two residual rows — zero exactly when the constraint
//! holds — and the solver's whole job is driving the stacked residual vector
//! to zero while the diagnosis reads the system's shape.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point2;

/// A point in the sketch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PointId(pub(crate) usize);

/// A line segment between two sketch points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LineId(pub(crate) usize);

/// A circle around a sketch point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CircleId(pub(crate) usize);

/// An arc: a centre and two rim points, coupled to one radius.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ArcId(pub(crate) usize);

/// A constraint in the sketch, in the order it was added.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ConstraintId(pub(crate) usize);

impl ConstraintId {
    /// The position of this constraint in the sketch's own ordering.
    #[must_use]
    pub fn index(&self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PointData {
    /// Offset of `(x, y)` in the parameter vector.
    pub at: usize,
    pub construction: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct LineData {
    pub a: PointId,
    pub b: PointId,
    pub construction: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct CircleData {
    pub centre: PointId,
    /// Offset of the radius in the parameter vector.
    pub radius_at: usize,
    pub construction: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct ArcData {
    pub centre: PointId,
    pub start: PointId,
    pub end: PointId,
    #[allow(dead_code)]
    pub construction: bool,
}

/// Which side a circle–circle tangency touches on.
///
/// Chosen when the constraint is added, from the geometry as it stands:
/// circles overlapping a shared interior are tangent internally, circles
/// apart tangent externally. Recording the side keeps the residual smooth —
/// a constraint that re-decided its own meaning every iteration would be a
/// discontinuity the solver falls into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangencySide {
    /// The circles touch from outside: centre distance equals the radii's sum.
    External,
    /// One circle touches from within: centre distance equals the radii's difference.
    Internal,
}

/// One constraint: the full 2D vocabulary.
///
/// Driving dimensions are constraints carrying their value — [`Constraint::Distance`],
/// [`Constraint::Angle`], [`Constraint::Radius`]. Driven dimensions are the
/// same quantities *read* instead of imposed: [`Sketch::measure_distance`]
/// and its siblings evaluate without entering the system.
#[derive(Debug, Clone, PartialEq)]
pub enum Constraint {
    /// Two points at the same place. Two rows.
    Coincident(PointId, PointId),
    /// A point pinned where it stands. Two rows.
    Fixed(PointId, Point2),
    /// The distance between two points.
    Distance(PointId, PointId, f64),
    /// The unsigned distance from a point to a line's carrier.
    PointLineDistance(PointId, LineId, f64),
    /// A point on a line's carrier.
    PointOnLine(PointId, LineId),
    /// A point on a circle.
    PointOnCircle(PointId, CircleId),
    /// A line parallel to the sketch x axis.
    Horizontal(LineId),
    /// A line parallel to the sketch y axis.
    Vertical(LineId),
    /// Two lines with parallel carriers.
    Parallel(LineId, LineId),
    /// Two lines with perpendicular carriers.
    Perpendicular(LineId, LineId),
    /// The counterclockwise angle from the first line to the second, radians.
    Angle(LineId, LineId, f64),
    /// Two lines of the same length.
    EqualLength(LineId, LineId),
    /// Two circles of the same radius.
    EqualRadius(CircleId, CircleId),
    /// A circle's radius.
    Radius(CircleId, f64),
    /// A line tangent to a circle.
    TangentLineCircle(LineId, CircleId),
    /// A line tangent to an arc's carrier circle.
    TangentLineArc(LineId, ArcId),
    /// Two circles tangent, on the recorded side.
    TangentCircles(CircleId, CircleId, TangencySide),
    /// Two points mirror images across a line. Two rows.
    Symmetric(PointId, PointId, LineId),
    /// An arc's two rim points at one radius from its centre.
    ///
    /// Added automatically by [`Sketch::add_arc`] — an arc with two radii is
    /// not an arc — and returned from it, so a diagnosis that names it names
    /// something the caller has seen.
    ArcRadii(ArcId),
}

/// How much of its weight a soft objective keeps: light enough that every
/// real constraint wins the fight, heavy enough to pick among the
/// configurations they leave open.
pub(crate) const SOFT_WEIGHT: f64 = 1e-3;

/// A 2D sketch: geometry, constraints, and the parameter vector under both.
#[derive(Debug, Clone, Default)]
pub struct Sketch {
    pub(crate) params: Vec<f64>,
    pub(crate) points: Vec<PointData>,
    pub(crate) lines: Vec<LineData>,
    pub(crate) circles: Vec<CircleData>,
    pub(crate) arcs: Vec<ArcData>,
    pub(crate) constraints: Vec<Constraint>,
    /// Whether the last constraint is a soft objective — a drag target —
    /// whose rows enter at a whisper of the weight.
    pub(crate) soft_last: bool,
}

impl Sketch {
    /// An empty sketch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a point at `at`.
    pub fn add_point(&mut self, at: Point2) -> PointId {
        let offset = self.params.len();
        self.params.push(at.x);
        self.params.push(at.y);
        self.points.push(PointData {
            at: offset,
            construction: false,
        });
        PointId(self.points.len() - 1)
    }

    /// Add a line between two existing points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// endpoints are the same point — a line needs two.
    pub fn add_line(&mut self, a: PointId, b: PointId) -> OgeomResult<LineId> {
        self.check_point(a)?;
        self.check_point(b)?;
        if a == b {
            ogeom_bail!(Construction, "a line needs two distinct points");
        }
        self.lines.push(LineData {
            a,
            b,
            construction: false,
        });
        Ok(LineId(self.lines.len() - 1))
    }

    /// Add a circle around an existing centre point.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// radius is not finite and positive.
    pub fn add_circle(&mut self, centre: PointId, radius: f64) -> OgeomResult<CircleId> {
        self.check_point(centre)?;
        if !radius.is_finite() || radius <= 0.0 {
            ogeom_bail!(Construction, "a circle radius of {radius} is not a size");
        }
        let radius_at = self.params.len();
        self.params.push(radius);
        self.circles.push(CircleData {
            centre,
            radius_at,
            construction: false,
        });
        Ok(CircleId(self.circles.len() - 1))
    }

    /// Add an arc from `start` to `end` around `centre`.
    ///
    /// The rim points must stand at one radius, so the arc brings its own
    /// [`Constraint::ArcRadii`] coupling; the returned constraint is that
    /// coupling, and a diagnosis may name it like any other.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// three points are not distinct.
    pub fn add_arc(
        &mut self,
        centre: PointId,
        start: PointId,
        end: PointId,
    ) -> OgeomResult<(ArcId, ConstraintId)> {
        self.check_point(centre)?;
        self.check_point(start)?;
        self.check_point(end)?;
        if centre == start || centre == end || start == end {
            ogeom_bail!(Construction, "an arc needs three distinct points");
        }
        self.arcs.push(ArcData {
            centre,
            start,
            end,
            construction: false,
        });
        let arc = ArcId(self.arcs.len() - 1);
        let coupling = self.constrain(Constraint::ArcRadii(arc))?;
        Ok((arc, coupling))
    }

    /// Mark a point as construction geometry.
    ///
    /// Construction geometry constrains and is constrained exactly like the
    /// rest — the flag says only that it is scaffolding, not profile, for
    /// whoever consumes the solved sketch.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn set_point_construction(&mut self, id: PointId, construction: bool) -> OgeomResult<()> {
        self.check_point(id)?;
        self.points[id.0].construction = construction;
        Ok(())
    }

    /// As [`Sketch::set_point_construction`], for a line.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn set_line_construction(&mut self, id: LineId, construction: bool) -> OgeomResult<()> {
        self.check_line(id)?;
        self.lines[id.0].construction = construction;
        Ok(())
    }

    /// As [`Sketch::set_point_construction`], for a circle.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn set_circle_construction(&mut self, id: CircleId, construction: bool) -> OgeomResult<()> {
        self.check_circle(id)?;
        self.circles[id.0].construction = construction;
        Ok(())
    }

    /// Whether a line is construction geometry.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn is_line_construction(&self, id: LineId) -> OgeomResult<bool> {
        self.check_line(id)?;
        Ok(self.lines[id.0].construction)
    }

    /// Add a constraint.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a referenced id
    /// is not from this sketch;
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
    /// dimension's value is not finite.
    pub fn constrain(&mut self, constraint: Constraint) -> OgeomResult<ConstraintId> {
        self.check_constraint(&constraint)?;
        self.constraints.push(constraint);
        Ok(ConstraintId(self.constraints.len() - 1))
    }

    /// The current position of a point.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn point(&self, id: PointId) -> OgeomResult<Point2> {
        self.check_point(id)?;
        let at = self.points[id.0].at;
        Ok(Point2::new(self.params[at], self.params[at + 1]))
    }

    /// The current radius of a circle.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn radius(&self, id: CircleId) -> OgeomResult<f64> {
        self.check_circle(id)?;
        Ok(self.params[self.circles[id.0].radius_at])
    }

    /// The endpoints of a line.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn line_ends(&self, id: LineId) -> OgeomResult<(Point2, Point2)> {
        self.check_line(id)?;
        let data = &self.lines[id.0];
        Ok((self.point(data.a)?, self.point(data.b)?))
    }

    /// The centre, start and end of an arc.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn arc_points(&self, id: ArcId) -> OgeomResult<(Point2, Point2, Point2)> {
        self.check_arc(id)?;
        let data = &self.arcs[id.0];
        Ok((
            self.point(data.centre)?,
            self.point(data.start)?,
            self.point(data.end)?,
        ))
    }

    /// A driven distance: the measurement, not the imposition.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if an id is not
    /// from this sketch.
    pub fn measure_distance(&self, a: PointId, b: PointId) -> OgeomResult<f64> {
        Ok(self.point(a)?.distance(self.point(b)?))
    }

    /// A driven angle between two lines' carriers, counterclockwise from the
    /// first, in `(-pi, pi]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if an id is not
    /// from this sketch.
    pub fn measure_angle(&self, a: LineId, b: LineId) -> OgeomResult<f64> {
        let (a0, a1) = self.line_ends(a)?;
        let (b0, b1) = self.line_ends(b)?;
        let u = (a1.x - a0.x, a1.y - a0.y);
        let v = (b1.x - b0.x, b1.y - b0.y);
        Ok((u.0 * v.1 - u.1 * v.0).atan2(u.0 * v.0 + u.1 * v.1))
    }

    /// A driven radius.
    ///
    /// # Errors
    ///
    /// As [`Sketch::radius`].
    pub fn measure_radius(&self, id: CircleId) -> OgeomResult<f64> {
        self.radius(id)
    }

    /// The constraints, in the order their ids name.
    #[must_use]
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// A one-line human description of a constraint, for diagnoses.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the id is not
    /// from this sketch.
    pub fn describe(&self, id: ConstraintId) -> OgeomResult<String> {
        let Some(constraint) = self.constraints.get(id.0) else {
            ogeom_bail!(Dangling, "constraint {} is not in this sketch", id.0);
        };
        Ok(match constraint {
            Constraint::Coincident(a, b) => format!("coincident(P{}, P{})", a.0, b.0),
            Constraint::Fixed(p, at) => format!("fixed(P{} at {:.3}, {:.3})", p.0, at.x, at.y),
            Constraint::Distance(a, b, d) => format!("distance(P{}, P{}) = {d}", a.0, b.0),
            Constraint::PointLineDistance(p, l, d) => {
                format!("distance(P{}, L{}) = {d}", p.0, l.0)
            }
            Constraint::PointOnLine(p, l) => format!("on(P{}, L{})", p.0, l.0),
            Constraint::PointOnCircle(p, c) => format!("on(P{}, C{})", p.0, c.0),
            Constraint::Horizontal(l) => format!("horizontal(L{})", l.0),
            Constraint::Vertical(l) => format!("vertical(L{})", l.0),
            Constraint::Parallel(a, b) => format!("parallel(L{}, L{})", a.0, b.0),
            Constraint::Perpendicular(a, b) => format!("perpendicular(L{}, L{})", a.0, b.0),
            Constraint::Angle(a, b, t) => format!("angle(L{}, L{}) = {t}", a.0, b.0),
            Constraint::EqualLength(a, b) => format!("equal-length(L{}, L{})", a.0, b.0),
            Constraint::EqualRadius(a, b) => format!("equal-radius(C{}, C{})", a.0, b.0),
            Constraint::Radius(c, r) => format!("radius(C{}) = {r}", c.0),
            Constraint::TangentLineCircle(l, c) => format!("tangent(L{}, C{})", l.0, c.0),
            Constraint::TangentLineArc(l, a) => format!("tangent(L{}, A{})", l.0, a.0),
            Constraint::TangentCircles(a, b, side) => {
                format!("tangent(C{}, C{}, {side:?})", a.0, b.0)
            }
            Constraint::Symmetric(a, b, l) => {
                format!("symmetric(P{}, P{} across L{})", a.0, b.0, l.0)
            }
            Constraint::ArcRadii(a) => format!("arc-radii(A{})", a.0),
        })
    }

    // --- residual machinery, used by the solver ------------------------------

    /// How many residual rows a constraint contributes.
    pub(crate) fn rows_of(constraint: &Constraint) -> usize {
        match constraint {
            Constraint::Coincident(..) | Constraint::Fixed(..) | Constraint::Symmetric(..) => 2,
            _ => 1,
        }
    }

    /// Evaluate every constraint's residuals at `params` into `out`.
    ///
    /// Direction residuals are sines of angles — dimensionless, bounded — so
    /// they are multiplied by `scale`, the sketch's characteristic length,
    /// to stand in the same units as the distance rows. Without that, the
    /// rank analysis would weigh a metre of error in one row against a
    /// radian's sine in another and call them comparable.
    pub(crate) fn residuals(&self, params: &[f64], scale: f64, out: &mut Vec<f64>) {
        out.clear();
        for (i, constraint) in self.constraints.iter().enumerate() {
            let start = out.len();
            self.constraint_residuals(constraint, params, scale, out);
            if self.soft_last && i + 1 == self.constraints.len() {
                for r in &mut out[start..] {
                    *r *= SOFT_WEIGHT;
                }
            }
        }
    }

    /// Append one constraint's residual rows.
    #[allow(clippy::too_many_lines)]
    pub(crate) fn constraint_residuals(
        &self,
        constraint: &Constraint,
        params: &[f64],
        scale: f64,
        out: &mut Vec<f64>,
    ) {
        let point = |id: PointId| -> (f64, f64) {
            let at = self.points[id.0].at;
            (params[at], params[at + 1])
        };
        let line = |id: LineId| -> ((f64, f64), (f64, f64)) {
            let data = &self.lines[id.0];
            (point(data.a), point(data.b))
        };
        // Direction of a line, normalized with a floor so a degenerate line
        // yields a large, smooth residual instead of a division blow-up.
        let unit = |a: (f64, f64), b: (f64, f64)| -> (f64, f64, f64) {
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let len = dx.hypot(dy).max(1e-12);
            (dx / len, dy / len, len)
        };
        let dist = |a: (f64, f64), b: (f64, f64)| (b.0 - a.0).hypot(b.1 - a.1);

        {
            match constraint {
                Constraint::Coincident(a, b) => {
                    let (pa, pb) = (point(*a), point(*b));
                    out.push(pb.0 - pa.0);
                    out.push(pb.1 - pa.1);
                }
                Constraint::Fixed(p, at) => {
                    let pp = point(*p);
                    out.push(pp.0 - at.x);
                    out.push(pp.1 - at.y);
                }
                Constraint::Distance(a, b, d) => {
                    out.push(dist(point(*a), point(*b)) - d);
                }
                Constraint::PointLineDistance(p, l, d) => {
                    let (a, b) = line(*l);
                    let (ux, uy, _) = unit(a, b);
                    let pp = point(*p);
                    let signed = ux * (pp.1 - a.1) - uy * (pp.0 - a.0);
                    // Unsigned distance, squared against the target so the
                    // residual is smooth on both sides of the carrier.
                    out.push((signed * signed - d * d) / d.abs().max(1.0));
                }
                Constraint::PointOnLine(p, l) => {
                    let (a, b) = line(*l);
                    let (ux, uy, _) = unit(a, b);
                    let pp = point(*p);
                    out.push(ux * (pp.1 - a.1) - uy * (pp.0 - a.0));
                }
                Constraint::PointOnCircle(p, c) => {
                    let data = &self.circles[c.0];
                    out.push(dist(point(*p), point(data.centre)) - params[data.radius_at]);
                }
                Constraint::Horizontal(l) => {
                    let (a, b) = line(*l);
                    out.push(b.1 - a.1);
                }
                Constraint::Vertical(l) => {
                    let (a, b) = line(*l);
                    out.push(b.0 - a.0);
                }
                Constraint::Parallel(l1, l2) => {
                    let (a1, b1) = line(*l1);
                    let (a2, b2) = line(*l2);
                    let (ux, uy, _) = unit(a1, b1);
                    let (vx, vy, _) = unit(a2, b2);
                    out.push((ux * vy - uy * vx) * scale);
                }
                Constraint::Perpendicular(l1, l2) => {
                    let (a1, b1) = line(*l1);
                    let (a2, b2) = line(*l2);
                    let (ux, uy, _) = unit(a1, b1);
                    let (vx, vy, _) = unit(a2, b2);
                    out.push((ux * vx + uy * vy) * scale);
                }
                Constraint::Angle(l1, l2, target) => {
                    let (a1, b1) = line(*l1);
                    let (a2, b2) = line(*l2);
                    let (ux, uy, _) = unit(a1, b1);
                    let (vx, vy, _) = unit(a2, b2);
                    let angle = (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
                    // Wrapped to the nearest turn, so an angle constraint
                    // pulls toward the closer of the two representations.
                    let mut delta = angle - target;
                    while delta > core::f64::consts::PI {
                        delta -= core::f64::consts::TAU;
                    }
                    while delta < -core::f64::consts::PI {
                        delta += core::f64::consts::TAU;
                    }
                    out.push(delta * scale);
                }
                Constraint::EqualLength(l1, l2) => {
                    let (a1, b1) = line(*l1);
                    let (a2, b2) = line(*l2);
                    out.push(dist(a1, b1) - dist(a2, b2));
                }
                Constraint::EqualRadius(c1, c2) => {
                    out.push(
                        params[self.circles[c1.0].radius_at] - params[self.circles[c2.0].radius_at],
                    );
                }
                Constraint::Radius(c, r) => {
                    out.push(params[self.circles[c.0].radius_at] - r);
                }
                Constraint::TangentLineCircle(l, c) => {
                    let (a, b) = line(*l);
                    let (ux, uy, _) = unit(a, b);
                    let data = &self.circles[c.0];
                    let cc = point(data.centre);
                    let signed = ux * (cc.1 - a.1) - uy * (cc.0 - a.0);
                    let r = params[data.radius_at];
                    out.push((signed * signed - r * r) / r.abs().max(1.0));
                }
                Constraint::TangentLineArc(l, arc) => {
                    let (a, b) = line(*l);
                    let (ux, uy, _) = unit(a, b);
                    let data = &self.arcs[arc.0];
                    let cc = point(data.centre);
                    let signed = ux * (cc.1 - a.1) - uy * (cc.0 - a.0);
                    let r = dist(cc, point(data.start));
                    out.push((signed * signed - r * r) / r.abs().max(1.0));
                }
                Constraint::TangentCircles(c1, c2, side) => {
                    let d1 = &self.circles[c1.0];
                    let d2 = &self.circles[c2.0];
                    let centres = dist(point(d1.centre), point(d2.centre));
                    let (r1, r2) = (params[d1.radius_at], params[d2.radius_at]);
                    let target = match side {
                        TangencySide::External => r1 + r2,
                        TangencySide::Internal => (r1 - r2).abs().max(1e-12),
                    };
                    out.push(centres - target);
                }
                Constraint::Symmetric(p, q, l) => {
                    let (a, b) = line(*l);
                    let (ux, uy, _) = unit(a, b);
                    let (pp, qq) = (point(*p), point(*q));
                    let mid = ((pp.0 + qq.0) / 2.0, (pp.1 + qq.1) / 2.0);
                    // The midpoint on the axis, and the join perpendicular
                    // to it.
                    out.push(ux * (mid.1 - a.1) - uy * (mid.0 - a.0));
                    out.push(ux * (qq.0 - pp.0) + uy * (qq.1 - pp.1));
                }
                Constraint::ArcRadii(arc) => {
                    let data = &self.arcs[arc.0];
                    let cc = point(data.centre);
                    out.push(dist(cc, point(data.start)) - dist(cc, point(data.end)));
                }
            }
        }
    }

    /// The parameter indices one constraint's residuals can touch —
    /// structural sparsity, exact by construction, because a constraint
    /// only ever reads the entities it names.
    pub(crate) fn parameters_of(&self, constraint: &Constraint) -> Vec<usize> {
        let mut out = Vec::new();
        let point = |id: PointId, out: &mut Vec<usize>| {
            let at = self.points[id.0].at;
            out.push(at);
            out.push(at + 1);
        };
        let line = |id: LineId, out: &mut Vec<usize>| {
            let data = self.lines[id.0].clone();
            point(data.a, out);
            point(data.b, out);
        };
        match constraint {
            Constraint::Coincident(a, b) | Constraint::Distance(a, b, _) => {
                point(*a, &mut out);
                point(*b, &mut out);
            }
            Constraint::Fixed(p, _) => point(*p, &mut out),
            Constraint::PointLineDistance(p, l, _) | Constraint::PointOnLine(p, l) => {
                point(*p, &mut out);
                line(*l, &mut out);
            }
            Constraint::PointOnCircle(p, c) => {
                point(*p, &mut out);
                let data = self.circles[c.0].clone();
                point(data.centre, &mut out);
                out.push(data.radius_at);
            }
            Constraint::Horizontal(l) | Constraint::Vertical(l) => line(*l, &mut out),
            Constraint::Parallel(a, b)
            | Constraint::Perpendicular(a, b)
            | Constraint::Angle(a, b, _)
            | Constraint::EqualLength(a, b) => {
                line(*a, &mut out);
                line(*b, &mut out);
            }
            Constraint::EqualRadius(a, b) => {
                out.push(self.circles[a.0].radius_at);
                out.push(self.circles[b.0].radius_at);
            }
            Constraint::Radius(c, _) => out.push(self.circles[c.0].radius_at),
            Constraint::TangentLineCircle(l, c) => {
                line(*l, &mut out);
                let data = self.circles[c.0].clone();
                point(data.centre, &mut out);
                out.push(data.radius_at);
            }
            Constraint::TangentLineArc(l, arc) => {
                line(*l, &mut out);
                let data = self.arcs[arc.0].clone();
                point(data.centre, &mut out);
                point(data.start, &mut out);
            }
            Constraint::TangentCircles(a, b, _) => {
                for c in [a, b] {
                    let data = self.circles[c.0].clone();
                    point(data.centre, &mut out);
                    out.push(data.radius_at);
                }
            }
            Constraint::Symmetric(p, q, l) => {
                point(*p, &mut out);
                point(*q, &mut out);
                line(*l, &mut out);
            }
            Constraint::ArcRadii(arc) => {
                let data = self.arcs[arc.0].clone();
                point(data.centre, &mut out);
                point(data.start, &mut out);
                point(data.end, &mut out);
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// The sketch's characteristic length: the parameter spread, floored at
    /// one, so dimensionless residual rows can stand in length units.
    pub(crate) fn characteristic_scale(&self) -> f64 {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for data in &self.points {
            for k in 0..2 {
                lo = lo.min(self.params[data.at + k]);
                hi = hi.max(self.params[data.at + k]);
            }
        }
        if hi > lo { (hi - lo).max(1.0) } else { 1.0 }
    }

    fn check_point(&self, id: PointId) -> OgeomResult<()> {
        if id.0 >= self.points.len() {
            ogeom_bail!(Dangling, "point {} is not in this sketch", id.0);
        }
        Ok(())
    }

    fn check_line(&self, id: LineId) -> OgeomResult<()> {
        if id.0 >= self.lines.len() {
            ogeom_bail!(Dangling, "line {} is not in this sketch", id.0);
        }
        Ok(())
    }

    fn check_circle(&self, id: CircleId) -> OgeomResult<()> {
        if id.0 >= self.circles.len() {
            ogeom_bail!(Dangling, "circle {} is not in this sketch", id.0);
        }
        Ok(())
    }

    fn check_arc(&self, id: ArcId) -> OgeomResult<()> {
        if id.0 >= self.arcs.len() {
            ogeom_bail!(Dangling, "arc {} is not in this sketch", id.0);
        }
        Ok(())
    }

    fn check_constraint(&self, constraint: &Constraint) -> OgeomResult<()> {
        let finite = |v: f64, what: &str| -> OgeomResult<()> {
            if !v.is_finite() {
                ogeom_bail!(Construction, "a {what} of {v} is not a dimension");
            }
            Ok(())
        };
        match constraint {
            Constraint::Coincident(a, b) => {
                self.check_point(*a)?;
                self.check_point(*b)
            }
            Constraint::Fixed(p, at) => {
                self.check_point(*p)?;
                finite(at.x, "position")?;
                finite(at.y, "position")
            }
            Constraint::Distance(a, b, d) => {
                self.check_point(*a)?;
                self.check_point(*b)?;
                finite(*d, "distance")
            }
            Constraint::PointLineDistance(p, l, d) => {
                self.check_point(*p)?;
                self.check_line(*l)?;
                finite(*d, "distance")
            }
            Constraint::PointOnLine(p, l) => {
                self.check_point(*p)?;
                self.check_line(*l)
            }
            Constraint::PointOnCircle(p, c) => {
                self.check_point(*p)?;
                self.check_circle(*c)
            }
            Constraint::Horizontal(l) | Constraint::Vertical(l) => self.check_line(*l),
            Constraint::Parallel(a, b)
            | Constraint::Perpendicular(a, b)
            | Constraint::EqualLength(a, b) => {
                self.check_line(*a)?;
                self.check_line(*b)
            }
            Constraint::Angle(a, b, t) => {
                self.check_line(*a)?;
                self.check_line(*b)?;
                finite(*t, "angle")
            }
            Constraint::EqualRadius(a, b) => {
                self.check_circle(*a)?;
                self.check_circle(*b)
            }
            Constraint::Radius(c, r) => {
                self.check_circle(*c)?;
                finite(*r, "radius")
            }
            Constraint::TangentLineCircle(l, c) => {
                self.check_line(*l)?;
                self.check_circle(*c)
            }
            Constraint::TangentLineArc(l, a) => {
                self.check_line(*l)?;
                self.check_arc(*a)
            }
            Constraint::TangentCircles(a, b, _) => {
                self.check_circle(*a)?;
                self.check_circle(*b)
            }
            Constraint::Symmetric(a, b, l) => {
                self.check_point(*a)?;
                self.check_point(*b)?;
                self.check_line(*l)
            }
            Constraint::ArcRadii(a) => self.check_arc(*a),
        }
    }
}
