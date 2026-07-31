//! Value-semantics geometric primitives and the numerical machinery under them.
//!
//! *Elsewhere:* the `gp` primitives, the B-spline basis libraries, elementary
//! curve/surface parameterization, and the solver half of a `math` package.
//!
//! Two things here are load-bearing and easy to get wrong:
//!
//! - `Dir` enforces its unit-length invariant in the type, not by convention.
//! - `Trsf` carries a **form classification** (identity / translation / rotation /
//!   scale / …), used to short-circuit transform application. Skipping it costs
//!   real performance across the whole kernel.
