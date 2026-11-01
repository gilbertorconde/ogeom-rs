//! Which points a fit keeps: the rule that drops outliers between rounds of
//! a fit, shared by every fit that iterates over a point set.

use ogeom_core::{OgeomResult, ogeom_bail};

/// Which points a fit keeps.
///
/// A fit keeps the points whose residual is no larger than a threshold, fits
/// again on those, and repeats until the kept set stops changing. The
/// threshold is measured over the points the last round kept.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Cutoff {
    /// Every point, in one round.
    #[default]
    All,
    /// The points within `k` times the RMS residual.
    Sigma(f64),
    /// The points within `k` times the median residual size.
    Median(f64),
    /// The points within a fixed distance.
    Distance(f64),
}

impl Cutoff {
    /// The rule's number checked: finite and positive.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the factor or
    /// the distance is not finite and positive.
    pub fn check(self) -> OgeomResult<()> {
        let (name, value) = match self {
            Self::All => return Ok(()),
            Self::Sigma(k) => ("sigma factor", k),
            Self::Median(k) => ("median factor", k),
            Self::Distance(d) => ("cutoff distance", d),
        };
        if !(value.is_finite() && value > 0.0) {
            ogeom_bail!(Domain, "the {name} {value} must be finite and positive");
        }
        Ok(())
    }

    /// The largest residual size kept, measured over the residuals of the
    /// points the last round kept; infinite for [`Cutoff::All`] and where
    /// there are none.
    #[must_use]
    pub fn threshold(self, kept: &[f64]) -> f64 {
        if kept.is_empty() {
            return f64::INFINITY;
        }
        #[allow(
            clippy::cast_precision_loss,
            reason = "point counts are far below 2^52"
        )]
        let count = kept.len() as f64;
        match self {
            Self::All => f64::INFINITY,
            Self::Sigma(k) => k * (kept.iter().map(|r| r * r).sum::<f64>() / count).sqrt(),
            Self::Median(k) => {
                let mut sizes: Vec<f64> = kept.iter().map(|r| r.abs()).collect();
                sizes.sort_by(f64::total_cmp);
                let mid = sizes.len() / 2;
                let median = if sizes.len().is_multiple_of(2) {
                    f64::midpoint(sizes[mid - 1], sizes[mid])
                } else {
                    sizes[mid]
                };
                k * median
            }
            Self::Distance(d) => d,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_rule_measures_its_threshold_over_the_kept_residuals() {
        let kept = [3.0, -4.0, 1.0, -2.0];
        assert!(Cutoff::All.threshold(&kept).is_infinite());
        let rms = (30.0_f64 / 4.0).sqrt();
        assert!((Cutoff::Sigma(2.0).threshold(&kept) - 2.0 * rms).abs() < 1e-15);
        assert!((Cutoff::Median(3.0).threshold(&kept) - 7.5).abs() < 1e-15);
        assert!((Cutoff::Distance(0.5).threshold(&kept) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn a_rule_without_a_positive_number_is_refused() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Cutoff::Sigma(bad).check().is_err());
            assert!(Cutoff::Median(bad).check().is_err());
            assert!(Cutoff::Distance(bad).check().is_err());
        }
        assert!(Cutoff::All.check().is_ok());
        assert!(Cutoff::Sigma(3.0).check().is_ok());
    }
}
