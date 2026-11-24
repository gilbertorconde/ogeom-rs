//! How far a mesh's storage may have moved its vertices.

use ogeom_topo::Triangulation;

/// How far single-precision storage may have moved a mesh's vertices: half
/// the spacing of `f32` values at its largest coordinate on each axis, over
/// the three axes at once. The quantum of a mesh that went through `f32`
/// (a render mesh, a binary STL).
#[must_use]
pub fn single_precision_quantum(mesh: &Triangulation) -> f64 {
    let largest = mesh
        .positions
        .iter()
        .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
        .fold(0.0_f64, f64::max);
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the magnitude is an f32 coordinate's"
    )]
    let x = (largest as f32).max(f32::MIN_POSITIVE);
    let step = f64::from(f32::from_bits(x.to_bits() + 1) - x);
    step * 3.0_f64.sqrt() / 2.0
}
