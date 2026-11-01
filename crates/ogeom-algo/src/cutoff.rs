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

    /// As [`Cutoff::threshold`], each residual counting by its weight: the
    /// RMS and the median of a distribution of weight rather than of
    /// count, so a fit over triangles weighs them by area. Equal weights
    /// give the unweighted threshold.
    pub(crate) fn threshold_weighted(self, kept: &[(f64, f64)]) -> f64 {
        let total: f64 = kept.iter().map(|(_, w)| w).sum();
        if kept.is_empty() || total.is_nan() || total <= 0.0 {
            return f64::INFINITY;
        }
        match self {
            Self::All => f64::INFINITY,
            Self::Sigma(k) => k * (kept.iter().map(|(r, w)| w * r * r).sum::<f64>() / total).sqrt(),
            Self::Median(k) => {
                let mut sizes: Vec<(f64, f64)> = kept.iter().map(|(r, w)| (r.abs(), *w)).collect();
                sizes.sort_by(|p, q| p.0.total_cmp(&q.0));
                let half = total * 0.5;
                let mut below = 0.0;
                for (i, (size, w)) in sizes.iter().enumerate() {
                    below += w;
                    if below > half {
                        return k * size;
                    }
                    if below == half {
                        return k * sizes
                            .get(i + 1)
                            .map_or(*size, |next| f64::midpoint(*size, next.0));
                    }
                }
                k * sizes.last().map_or(0.0, |s| s.0)
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
    fn equal_weights_give_the_unweighted_threshold_and_weights_move_it() {
        let kept = [3.0, -4.0, 1.0, -2.0];
        let even: Vec<(f64, f64)> = kept.iter().map(|r| (*r, 0.5)).collect();
        for rule in [
            Cutoff::Sigma(2.0),
            Cutoff::Median(3.0),
            Cutoff::Distance(0.5),
        ] {
            let (a, b) = (rule.threshold(&kept), rule.threshold_weighted(&even));
            assert!((a - b).abs() <= 1e-15 * a, "{rule:?}: {a} and {b}");
        }
        let odd = [(1.0, 1.0), (2.0, 1.0), (9.0, 1.0)];
        assert!((Cutoff::Median(1.0).threshold_weighted(&odd) - 2.0).abs() < 1e-15);
        // The large residual weighs as much as the others together and
        // more: the median moves to it, and the RMS follows.
        let heavy = [(1.0, 1.0), (2.0, 1.0), (9.0, 3.0)];
        assert!((Cutoff::Median(1.0).threshold_weighted(&heavy) - 9.0).abs() < 1e-15);
        let rms = ((1.0 + 4.0 + 3.0 * 81.0) / 5.0_f64).sqrt();
        assert!((Cutoff::Sigma(1.0).threshold_weighted(&heavy) - rms).abs() < 1e-14);
        assert!(Cutoff::All.threshold_weighted(&heavy).is_infinite());
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
