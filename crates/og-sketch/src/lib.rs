//! 2D geometric constraint solving.
//!
//! Coincidence, distance, angle, parallel, perpendicular, tangent, symmetry,
//! equality, radius, horizontal and vertical; construction geometry; driving and
//! driven dimensions; degree-of-freedom analysis.
//!
//! Diagnosis matters as much as solving. An over-constrained sketch must report
//! *which* constraints conflict, and an under-constrained one *which* degrees of
//! freedom remain — a solver that only reports success or failure pushes the real
//! work back onto the user.
//!
//! No conventional CAD kernel ships this; applications built on them each supply
//! their own. A kernel meant to be a complete foundation should not make every
//! consumer solve it again.
