//! The general surface/surface intersector: seed, then walk.
//!
//! Where two surfaces meet has no closed form in general: the curve is
//! transcendental, and `docs/DATA_MODEL.md` §9 is blunt about the consequence:
//! there is no exact answer to be exact about, which is why the topology carries
//! tolerances. What there is instead is a curve that can be *followed*, one
//! corrected step at a time, to a stated accuracy.
//!
//! # Two problems, kept apart on purpose
//!
//! **Finding a branch** and **following one** fail in completely different ways,
//! and lumping them together is how an intersector comes to look better than it
//! is. A tracer that follows one branch beautifully while never noticing the
//! second reports a smooth, accurate, *wrong* answer, and the obvious accuracy
//! measure, "is every point on both surfaces", scores it perfectly.
//!
//! So [`seeds`] and [`trace`] are separate, separately testable, and separately
//! measured. Seeding is polyhedral: both surfaces are sampled into triangles and
//! the triangle pairs that cross give starting points. A branch thinner than
//! the grid crosses no pair of triangles; where it meets a spline's border it
//! is found as a piercing, and a closed loop inside both surfaces is found
//! from the approach of the surfaces inside it. A branch running along a
//! border of either chart, where one surface ends in the other, is found by
//! sampling the border itself. Where the surfaces lie close
//! along a whole curve only part of that search is made, and there the grid
//! is the knob.
//!
//! # Following the curve
//!
//! At a point on both surfaces the intersection runs along the cross product of
//! the two normals: the one direction that stays in both tangent planes. Step
//! along it and you leave both surfaces slightly; a Newton correction brings you
//! back.
//!
//! The correction has four unknowns (two parameters on each surface) and three
//! equations, `A(u1,v1) = B(u2,v2)`. That is deliberately one short, because the
//! solution set *is* the curve and pinning it to a point needs one more
//! condition. The fourth is a plane across the direction of travel: it says how
//! far along to land, and it is what turns "somewhere on the curve" into "the
//! next point".
//!
//! # What it reports about itself
//!
//! Whether the curve closed, and whether it ran out of steps. A polyline that
//! stopped because it hit a limit is not the same answer as one that stopped
//! because the curve ended, and a caller that cannot tell them apart will treat
//! a truncated branch as a complete one.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface, SurfaceGeometry};
use ogeom_math::{Point, Vector, solve};

mod borders;
mod branching;
mod cells;
mod cut;
mod follow;
mod middles;
mod poles;
mod seeding;
mod settle;
mod split;
mod stitch;
mod touches;

use borders::{border_stretches, extrusion_borders, rims, spline_borders};
pub(crate) use branching::branches_from_touches;
pub use branching::{branches, branches_cut_at_poles};
pub(crate) use cells::{Cell, sample, sample_by, segment_meets_triangle};
use cells::{CellBins, overlap, triangles_cross};
use cut::cut_at_touches;
pub use follow::trace;
use follow::{
    Measured, SurfacePair, correct, normal_at, outside, tangent_at, tangent_of_normals,
    touch_on_edge, unit_normal, walk,
};
pub(crate) use follow::{clamp, reversed, span};
use middles::{carry_onto, fill};
pub(crate) use middles::{nearest_on, with_middles};
pub(crate) use poles::{BRANCH_POINT_SINE, distance_to_segment, joining_reach};
pub use poles::{chart_poles, cut_at_poles, passes_a_pole};
use poles::{crossing_sine, first_lap, is_fragment, passes_near, turning};
pub use seeding::seeds;
use seeding::{BOW, FAR, Seeded, point_to_triangle, seeded};
pub(crate) use seeding::{
    STATIONARY_STEPS, TANGENTIAL_SINE, crosses_measurably, leaves_at_an_angle,
};
use settle::{
    Landing, along_arm, arc_length, arm_point, arm_to_border, find_touches, settle_touches,
    settled_from, straightened,
};
use split::{
    freed_from_borders, half_arm, halves, off_the_section, on_arm, split_at_touches,
    straight_through,
};
use stitch::{interior_is_transversal, stitch_stalled};
use touches::{
    Touch, apart, beside, blur_of, crosses, same_touch, subdivided_touches, touch_arms,
    touching_point,
};

/// How hard to look, and how closely to follow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Marching {
    /// How far the polyline may sit from the true curve, in space.
    pub chord: f64,
    /// How finely each surface is sampled when looking for branches.
    ///
    /// A branch narrower than one cell crosses no pair of cells, and is found
    /// only from a border it pierces or from the approach of the surfaces
    /// near it. Raising this costs time quadratically and makes a thin branch
    /// cross cells.
    pub grid: usize,
    /// A ceiling on the points in one branch, so a curve that will not close
    /// cannot run forever.
    pub max_points: usize,
}

impl Default for Marching {
    fn default() -> Self {
        Self {
            chord: 1e-4,
            grid: 24,
            max_points: 20_000,
        }
    }
}

impl Marching {
    /// Check the settings are usable.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the chord is
    /// not positive, the grid is too coarse to hold a triangle, or no points are
    /// allowed.
    pub fn validate(&self) -> OgeomResult<()> {
        if !self.chord.is_finite() || self.chord <= 0.0 {
            ogeom_bail!(Construction, "a chord of {} is not a distance", self.chord);
        }
        if self.grid < 2 {
            ogeom_bail!(Construction, "a sampling grid needs at least two steps");
        }
        if self.max_points < 2 {
            ogeom_bail!(Construction, "a branch needs at least two points");
        }
        Ok(())
    }
}

/// A point that lies on both surfaces, with where it is on each.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// Parameters on the first surface.
    pub on_a: (f64, f64),
    /// Parameters on the second.
    pub on_b: (f64, f64),
    /// Where that is in space.
    pub point: Point,
}

/// Why a traced branch stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stopped {
    /// It came back to where it started.
    Closed,
    /// It reached the edge of one surface's domain.
    LeftTheDomain,
    /// The correction stopped converging: a tangency or a singular point.
    ///
    /// Reported rather than pushed through. Marching past a point where the two
    /// normals are parallel is how a tracer jumps onto the wrong branch, and a
    /// wrong branch is a plausible answer to a different question.
    Stalled,
    /// It hit [`Marching::max_points`].
    ///
    /// Distinct from every other reason, because this one means the answer is
    /// *incomplete* rather than finished.
    RanOut,
}

/// One traced branch.
#[derive(Debug, Clone, PartialEq)]
pub struct Traced {
    /// The points along it, in order.
    pub points: Vec<Point>,
    /// Where each point is on the first surface.
    pub on_a: Vec<(f64, f64)>,
    /// Where each point is on the second.
    pub on_b: Vec<(f64, f64)>,
    /// Why it stopped.
    pub stopped: Stopped,
}

impl Traced {
    /// Whether the branch is finished rather than truncated.
    #[must_use]
    pub const fn complete(&self) -> bool {
        !matches!(self.stopped, Stopped::RanOut)
    }

    /// Whether it closed on itself.
    #[must_use]
    pub const fn closed(&self) -> bool {
        matches!(self.stopped, Stopped::Closed)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests;
