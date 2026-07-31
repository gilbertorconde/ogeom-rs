//! Parametric geometry: the `Curve3d`, `Curve2d` and `Surface` traits plus every
//! concrete type behind them.
//!
//! *Elsewhere:* the `Geom` and `Geom2d` type vocabularies and the `Adaptor`
//! family.
//!
//! The trait layer is the point, and it is the best idea in the conventional
//! design: every intersection, projection and extrema algorithm is written against
//! an *adaptor*, so nothing downstream cares whether the input is an analytic
//! cylinder or a NURBS patch. Concrete types implement the traits; algorithms only
//! ever see the traits.
