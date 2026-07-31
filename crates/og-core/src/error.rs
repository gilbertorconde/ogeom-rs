//! Errors as values.
//!
//! See `docs/DATA_MODEL.md` §12. The variants cover the failure vocabulary a
//! kernel needs, chosen to line up with the categories applications already
//! handle — but they are returned, not thrown, and no hardware signal is ever
//! converted into one of them.
//!
//! The rule this encodes: an algorithm that did not converge says so. It does
//! not return a null shape and set a flag for the caller to forget to check.

use core::fmt;

/// The result of any fallible kernel operation.
pub type OgResult<T> = Result<T, OgError>;

/// A kernel failure.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum OgError {
    /// Arguments cannot produce the requested entity — three collinear points
    /// for a circle, a zero-length direction, a self-intersecting wire.
    ///
    #[error("construction failed: {0}")]
    Construction(Cause),

    /// An argument lies outside the domain the operation is defined on.
    #[error("domain error: {0}")]
    Domain(Cause),

    /// A parameter or index lies outside its valid range.
    #[error("out of range: {0}")]
    Range(Cause),

    /// Collections or geometries that had to agree in size or dimension did not.
    #[error("dimension mismatch: {0}")]
    Dimension(Cause),

    /// An operation was handed a null or empty shape where one was required.
    #[error("null object: {0}")]
    NullObject(Cause),

    /// A key did not resolve — most often a stale arena key, or a shape used
    /// with an arena that does not own it.
    #[error("dangling reference: {0}")]
    Dangling(Cause),

    /// A numerical method failed: no convergence, a singular system, a step
    /// that could not be taken.
    #[error("numeric failure: {0}")]
    Numeric(Cause),

    /// The algorithm ran but could not produce a result. Distinct from
    /// [`OgError::Construction`]: the inputs were legitimate and the failure is
    /// the algorithm's.
    #[error("not done: {0}")]
    NotDone(Cause),

    /// The result exists but violates an invariant from `docs/DATA_MODEL.md` —
    /// tolerance containment, orientation consistency, edge representation
    /// agreement. Never returned silently; producing invalid topology is worse
    /// than failing.
    #[error("invariant violated: {0}")]
    Invariant(Cause),

    /// Cancelled through a progress sink.
    #[error("cancelled")]
    Cancelled,

    /// Reached a path that exists but has not been implemented yet.
    #[error("not implemented: {0}")]
    Unimplemented(Cause),
}

impl OgError {
    /// Whether retrying with a looser tolerance could plausibly succeed.
    ///
    /// Used by the fuzzy-tolerance escape hatch: numerical failures are worth
    /// retrying, malformed input is not.
    #[must_use]
    pub const fn is_tolerance_sensitive(&self) -> bool {
        matches!(
            self,
            Self::Numeric(_) | Self::NotDone(_) | Self::Invariant(_)
        )
    }
}

/// A short, cheap explanation attached to an [`OgError`].
///
/// Static in the common case so that returning an error costs no allocation on
/// paths that are hit often — failed intersections inside a boolean, for
/// instance, are routine control flow rather than exceptional.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Cause {
    /// A compile-time message.
    Static(&'static str),
    /// A message built at runtime.
    Owned(String),
}

impl fmt::Display for Cause {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Static(s) => f.write_str(s),
            Self::Owned(s) => f.write_str(s),
        }
    }
}

impl From<&'static str> for Cause {
    fn from(s: &'static str) -> Self {
        Self::Static(s)
    }
}

impl From<String> for Cause {
    fn from(s: String) -> Self {
        Self::Owned(s)
    }
}

impl From<fmt::Arguments<'_>> for Cause {
    fn from(a: fmt::Arguments<'_>) -> Self {
        a.as_str()
            .map_or_else(|| Self::Owned(a.to_string()), Self::Static)
    }
}

/// Build an [`OgError`] with a formatted cause.
///
/// ```
/// # use og_core::{og_err, OgError};
/// let e = og_err!(Construction, "radius {} is not positive", -1.0);
/// assert!(matches!(e, OgError::Construction(_)));
/// ```
#[macro_export]
macro_rules! og_err {
    ($variant:ident, $msg:literal) => {
        $crate::OgError::$variant($crate::Cause::Static($msg))
    };
    ($variant:ident, $fmt:literal, $($arg:tt)*) => {
        $crate::OgError::$variant($crate::Cause::Owned(format!($fmt, $($arg)*)))
    };
}

/// Return early with an [`OgError`] built by [`og_err!`].
#[macro_export]
macro_rules! og_bail {
    ($variant:ident, $($arg:tt)*) => {
        return Err($crate::og_err!($variant, $($arg)*))
    };
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn static_cause_does_not_allocate() {
        let e = og_err!(Domain, "parameter outside curve range");
        assert!(matches!(e, OgError::Domain(Cause::Static(_))));
    }

    #[test]
    fn formatted_cause_renders() {
        let e = og_err!(Range, "index {} of {}", 7, 3);
        assert_eq!(e.to_string(), "out of range: index 7 of 3");
    }

    #[test]
    fn tolerance_sensitivity_splits_input_errors_from_algorithm_errors() {
        assert!(og_err!(Numeric, "no convergence").is_tolerance_sensitive());
        assert!(og_err!(NotDone, "could not close shell").is_tolerance_sensitive());
        // Retrying a degenerate construction with a looser tolerance is pointless.
        assert!(!og_err!(Construction, "collinear points").is_tolerance_sensitive());
        assert!(!OgError::Cancelled.is_tolerance_sensitive());
    }

    #[test]
    fn bail_returns_early() {
        fn f(ok: bool) -> OgResult<u8> {
            if !ok {
                og_bail!(NullObject, "no shape");
            }
            Ok(1)
        }
        assert_eq!(f(true).unwrap(), 1);
        assert!(f(false).is_err());
    }
}
