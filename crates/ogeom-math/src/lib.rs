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

pub mod bounds;
pub mod bspline;
pub mod conic;
pub mod direction;
pub mod elementary;
pub mod frame;
pub mod integrate;
pub mod interval;
pub mod knots;
pub mod matrix;
pub mod point;
pub mod quadric;
pub mod quaternion;
pub mod solve;
pub mod sparse;
pub mod transform;
pub mod vector;

pub use bounds::Aabb;
pub use bspline::{Blend, ControlGrid, Weighted};
pub use conic::{Circle, Circle2, Ellipse, Ellipse2, Hyperbola, Hyperbola2, Parabola, Parabola2};
pub use direction::{Direction, Direction2};
pub use elementary::{CurvePoint, SurfacePoint};
pub use frame::{Axis, Axis2, Frame, Frame2, Handedness};
pub use integrate::{gauss_legendre, integrate};
pub use interval::Interval;
pub use knots::{BasisValues, KnotVector};
pub use matrix::{Matrix2, Matrix3};
pub use point::{Point, Point2};
pub use quadric::{Cone, Cylinder, Plane, Sphere, Torus, TorusKind};
pub use quaternion::Quaternion;
pub use solve::{Convergence, Criteria, Solution, SystemSolution};
pub use sparse::{SparseMatrix, least_squares_cgnr};
pub use transform::{GeneralTransform, Transform, Transform2, TransformKind};
pub use vector::{Vector, Vector2};
