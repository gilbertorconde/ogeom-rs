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

pub mod convert;
pub mod curve;
pub mod curve2d;
pub mod fit;
pub mod surface;
pub mod traits;

pub use curve::{
    BSplineCurve, CircleCurve, Curve, EllipseCurve, HelixCurve, HyperbolaCurve, LineCurve,
    ParabolaCurve, TrimmedCurve,
};
pub use curve2d::{BSpline2d, Circle2d, Ellipse2d, Line2d, PlanarCurve, Trimmed2d, tangent_angle};
pub use surface::{
    BSplineSurface, ConeSurface, CylinderSurface, ExtrusionSurface, PlaneSurface,
    RevolutionSurface, SphereSurface, SurfaceGeometry, TorusSurface, TrimmedSurface,
};
pub use traits::{
    Continuity, Curve2d, Curve3d, CurveKind, Reversible, Surface, SurfaceKind, Transformable,
};
