//! Tolerances.
//!
//! See `docs/DATA_MODEL.md` §5 and §9. Two separate things live here:
//!
//! - **Global constants**: the thresholds at which the kernel decides two
//!   things are "the same". [`Tolerances`] carries them, parameterised by the
//!   model's unit scale.
//! - **Per-entity tolerances**: [`Tolerance`], the radius attached to an
//!   individual vertex, edge or face, together with the containment rule that
//!   relates them.
//!
//! The unit scale is **explicit**. A confusion tolerance hard-coded at `1e-7`
//! carries an unstated assumption that models are in millimetres, and
//! misbehaves silently on a model authored in metres or inches.

use crate::{OgeomResult, ogeom_bail};

/// The threshold constants the kernel decides identity by, for a given model
/// scale.
///
/// `linear_scale` is the length of one model unit in millimetres: `1.0` for a
/// model in millimetres, `1000.0` for metres, `25.4` for inches. Linear
/// tolerances scale with it. Angular and parametric ones do not.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerances {
    linear_scale: f64,
}

/// Two points closer than this are the same point, at unit scale.
pub const CONFUSION: f64 = 1e-7;
/// Two directions closer than this in angle are parallel. Dimensionless.
///
/// Deliberately tight: near the limit of what `f64` can distinguish. It is the
/// right threshold for comparing directions that were *stored*, such as two
/// surface axes read back from a model. It is the wrong one for comparing
/// directions *computed* through subtraction of nearby coordinates, where the
/// input's own rounding is already larger. Such a comparison needs a bound
/// derived from the magnitudes involved.
pub const ANGULAR: f64 = 1e-12;
/// Convergence target for intersection algorithms, at unit scale.
pub const INTERSECTION: f64 = CONFUSION * 1e-2;
/// Accuracy target when fitting a curve or surface to data, at unit scale.
pub const APPROXIMATION: f64 = CONFUSION * 1e1;
/// Confusion in parametric space. Dimensionless.
pub const P_CONFUSION: f64 = CONFUSION * 1e-2;
/// A length below which an entity is degenerate rather than merely small, at
/// unit scale.
pub const DEGENERATE: f64 = CONFUSION * 1e-1;

impl Default for Tolerances {
    fn default() -> Self {
        Self::millimetres()
    }
}

impl Tolerances {
    /// Tolerances for a model whose unit is one millimetre.
    #[must_use]
    pub const fn millimetres() -> Self {
        Self { linear_scale: 1.0 }
    }

    /// Tolerances for a model whose unit is one metre.
    #[must_use]
    pub const fn metres() -> Self {
        Self {
            linear_scale: 1000.0,
        }
    }

    /// Tolerances for a model whose unit is one inch.
    #[must_use]
    pub const fn inches() -> Self {
        Self { linear_scale: 25.4 }
    }

    /// Tolerances for a model whose unit is `mm_per_unit` millimetres.
    ///
    /// Every linear tolerance derived from the scale, and its square, must be
    /// a normal `f64`: comparisons are made both on distances and on squared
    /// distances, and an infinite, zero or subnormal threshold would decide
    /// them by overflow or underflow. That accepts scales from about `7.5e-161`
    /// to about `6.7e144` millimetres per unit.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](crate::OgeomError::Construction) if the scale is
    /// not finite and positive, or is outside that range.
    pub fn with_scale(mm_per_unit: f64) -> OgeomResult<Self> {
        if !mm_per_unit.is_finite() || mm_per_unit <= 0.0 {
            ogeom_bail!(
                Construction,
                "unit scale {mm_per_unit} must be finite and positive"
            );
        }
        let scaled = Self {
            linear_scale: mm_per_unit,
        };
        let linear = [
            scaled.confusion(),
            scaled.intersection(),
            scaled.approximation(),
            scaled.degenerate(),
        ];
        if !linear.iter().all(|t| t.is_normal() && (t * t).is_normal()) {
            ogeom_bail!(
                Construction,
                "unit scale {mm_per_unit} puts a linear tolerance or its square \
                 outside the normal range of f64"
            );
        }
        Ok(scaled)
    }

    /// Millimetres per model unit.
    #[must_use]
    pub const fn scale(self) -> f64 {
        self.linear_scale
    }

    /// Distance below which two points are the same point.
    #[must_use]
    pub fn confusion(self) -> f64 {
        CONFUSION / self.linear_scale
    }

    /// Angle below which two directions are parallel. Independent of scale.
    #[must_use]
    pub const fn angular(self) -> f64 {
        ANGULAR
    }

    /// Convergence target for intersection algorithms.
    #[must_use]
    pub fn intersection(self) -> f64 {
        INTERSECTION / self.linear_scale
    }

    /// Accuracy target when fitting curves and surfaces.
    #[must_use]
    pub fn approximation(self) -> f64 {
        APPROXIMATION / self.linear_scale
    }

    /// Confusion in parametric space. Independent of scale.
    #[must_use]
    pub const fn parametric(self) -> f64 {
        P_CONFUSION
    }

    /// Length below which an entity is degenerate.
    #[must_use]
    pub fn degenerate(self) -> f64 {
        DEGENERATE / self.linear_scale
    }

    /// Whether two lengths are indistinguishable.
    #[must_use]
    pub fn same_length(self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.confusion()
    }

    /// Whether two parameters are indistinguishable.
    #[must_use]
    pub fn same_parameter(self, a: f64, b: f64) -> bool {
        (a - b).abs() <= self.parametric()
    }
}

/// The tolerance carried by a single vertex, edge or face: the radius within
/// which the entity is considered to lie.
///
/// Always finite and non-negative. Operations may only widen it (see
/// [`Tolerance::widen`]) because narrowing a tolerance asserts an accuracy the
/// geometry does not have.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Tolerance(f64);

impl Tolerance {
    /// The smallest meaningful tolerance, at unit scale.
    pub const MIN: Self = Self(CONFUSION);

    /// A tolerance of `value`, clamped up to [`Tolerance::MIN`].
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](crate::OgeomError::Construction) if `value` is not
    /// finite or is negative. A NaN tolerance poisons every comparison it
    /// reaches, so it is rejected at the boundary rather than propagated.
    pub fn new(value: f64) -> OgeomResult<Self> {
        if !value.is_finite() || value < 0.0 {
            ogeom_bail!(
                Construction,
                "tolerance {value} must be finite and non-negative"
            );
        }
        Ok(Self(value.max(CONFUSION)))
    }

    /// The tolerance as a length.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }

    /// This tolerance widened to at least `other`.
    ///
    /// The only sanctioned way to change a tolerance. Boolean operations grow
    /// tolerances as they go. Nothing shrinks them.
    #[must_use]
    pub fn widen(self, other: Self) -> Self {
        Self(self.0.max(other.0))
    }

    /// This tolerance widened to at least `value`, ignoring non-finite input.
    #[must_use]
    pub fn widen_to(self, value: f64) -> Self {
        if value.is_finite() {
            Self(self.0.max(value))
        } else {
            self
        }
    }

    /// Whether a separation of `distance` is within this tolerance.
    #[must_use]
    pub fn covers(self, distance: f64) -> bool {
        distance.abs() <= self.0
    }
}

impl Default for Tolerance {
    fn default() -> Self {
        Self::MIN
    }
}

/// Check the containment rule `tol(vertex) >= tol(edge) >= tol(face)` for a
/// boundary relationship.
///
/// `docs/DATA_MODEL.md` §5. A parent entity's boundary must be at least as
/// uncertain as the entity it bounds, or the boundary does not reliably lie on
/// it.
///
/// # Errors
///
/// [`OgeomError::Invariant`](crate::OgeomError::Invariant) if `bounding` is tighter
/// than `bounded`.
pub fn check_containment(bounding: Tolerance, bounded: Tolerance) -> OgeomResult<()> {
    if bounding.get() < bounded.get() {
        ogeom_bail!(
            Invariant,
            "tolerance containment violated: bounding {} < bounded {}",
            bounding.get(),
            bounded.get()
        );
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn linear_tolerances_scale_but_angular_ones_do_not() {
        let mm = Tolerances::millimetres();
        let m = Tolerances::metres();
        // A model in metres has 1000x coarser numbers for the same physical
        // distance, so the tolerance expressed in model units is 1000x smaller.
        assert!((m.confusion() * 1000.0 - mm.confusion()).abs() < 1e-18);
        assert_eq!(m.angular(), mm.angular());
        assert_eq!(m.parametric(), mm.parametric());
    }

    #[test]
    fn rejects_nonsense_scales() {
        assert!(Tolerances::with_scale(0.0).is_err());
        assert!(Tolerances::with_scale(-1.0).is_err());
        assert!(Tolerances::with_scale(f64::NAN).is_err());
        assert!(Tolerances::with_scale(f64::INFINITY).is_err());
        assert!(Tolerances::with_scale(25.4).is_ok());
    }

    /// Micrometres to kilometres, every derived linear tolerance is the
    /// constant over the scale.
    #[test]
    fn every_working_unit_derives_its_tolerances() {
        for scale in [1e-3, 1.0, 25.4, 1e3, 1e6] {
            let t = Tolerances::with_scale(scale).unwrap();
            assert_eq!(t.confusion(), CONFUSION / scale);
            assert_eq!(t.intersection(), INTERSECTION / scale);
            assert_eq!(t.approximation(), APPROXIMATION / scale);
            assert_eq!(t.degenerate(), DEGENERATE / scale);
            assert!(t.same_length(1.0, 1.0 + 0.5 * CONFUSION / scale));
            assert!(!t.same_length(1.0, 1.0 + 2.0 * CONFUSION / scale));
        }
    }

    /// A scale whose tolerances overflow, underflow or square out of range
    /// is refused by name; the extremes it accepts still compare distances.
    #[test]
    fn a_scale_whose_tolerances_leave_the_normal_range_is_refused() {
        for scale in [f64::MIN_POSITIVE, 1e-200, 1e-161, 1e145, 1e200, f64::MAX] {
            assert!(Tolerances::with_scale(scale).is_err(), "{scale}");
        }
        for scale in [1e-159, 1e143] {
            let t = Tolerances::with_scale(scale).unwrap();
            let c = t.confusion();
            assert!(t.same_length(0.0, 0.5 * c) && !t.same_length(0.0, 2.0 * c));
            let (near, far) = (0.5 * c, 2.0 * c);
            assert!(near * near <= c * c && far * far > c * c);
        }
    }

    #[test]
    fn tolerance_never_drops_below_min() {
        assert_eq!(Tolerance::new(0.0).unwrap(), Tolerance::MIN);
        assert_eq!(Tolerance::new(1e-30).unwrap(), Tolerance::MIN);
    }

    #[test]
    fn nan_tolerance_is_rejected_not_propagated() {
        assert!(Tolerance::new(f64::NAN).is_err());
        assert!(Tolerance::new(-1.0).is_err());
        assert!(Tolerance::new(f64::INFINITY).is_err());
    }

    #[test]
    fn widening_is_monotone() {
        let a = Tolerance::new(1e-4).unwrap();
        let b = Tolerance::new(1e-2).unwrap();
        assert_eq!(a.widen(b), b);
        assert_eq!(b.widen(a), b, "widen must never shrink");
        assert_eq!(a.widen_to(f64::NAN), a, "non-finite input must not poison");
    }

    #[test]
    fn containment_rule() {
        let vertex = Tolerance::new(1e-3).unwrap();
        let edge = Tolerance::new(1e-4).unwrap();
        let face = Tolerance::new(1e-5).unwrap();
        assert!(check_containment(vertex, edge).is_ok());
        assert!(check_containment(edge, face).is_ok());
        assert!(check_containment(face, vertex).is_err());
    }
}
