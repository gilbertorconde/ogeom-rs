//! The generalized winding number of a point with respect to a boundary
//! (Jacobson, Kavan and Sorkine-Hornung, "Robust inside-outside
//! segmentation using generalized winding numbers", 2013).
//!
//! The solid angle a closed, outward-wound boundary subtends at a point,
//! over a full sphere's, is one inside and zero outside, and it does not
//! depend on any ray: a point every ray from which grazes something still
//! has one. Over a mesh it is a sum over triangles, exact for the mesh. A
//! mesh drawn from faces stands off them by up to its chord, and a closed
//! boundary deformed without passing through the point keeps its winding
//! number, so the mesh's answer is the solid's wherever the point stands
//! clear of the band between each face and its mesh. A plane face and its
//! mesh lie in one plane, and the band does too.
//!
//! An open or self-crossing boundary reads off the integers or past one,
//! and is not answered.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Aabb, Point, Vector};
use ogeom_mesh::{Deflection, FaceMeshCache, face_meshes_for};
use ogeom_topo::{Model, NodeData, Shape};

use crate::classify::Containment;

/// How far from zero or one a winding number may read and still answer.
///
/// A closed mesh reads an integer to rounding; what moves it is a seam
/// between two faces' meshes sampled a tolerance apart, which costs a
/// small fraction of this. A boundary missing a face, or one crossing
/// itself, reads further off, or past one.
pub(crate) const WINDING_MARGIN: f64 = 1e-2;

/// The signed solid angle the triangle `t` subtends at `p`, positive where
/// `p` stands behind it, on the side its winding faces away from (Van
/// Oosterom and Strackee 1983).
pub(crate) fn solid_angle(p: Point, t: [Point; 3]) -> f64 {
    let (a, b, c) = (t[0] - p, t[1] - p, t[2] - p);
    let (la, lb, lc) = (a.magnitude(), b.magnitude(), c.magnitude());
    let numerator = a.dot(b.cross(c));
    let denominator = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
    2.0 * numerator.atan2(denominator)
}

/// The winding number of `p` with respect to `triangles`, wound outward.
pub(crate) fn winding_number<'a>(p: Point, triangles: impl Iterator<Item = &'a [Point; 3]>) -> f64 {
    triangles.map(|t| solid_angle(p, *t)).sum::<f64>() / (4.0 * core::f64::consts::PI)
}

/// What a winding number says: In within [`WINDING_MARGIN`] of one, Out
/// within it of zero, nothing otherwise.
pub(crate) fn reading(winding: f64) -> Option<Containment> {
    if (winding - 1.0).abs() <= WINDING_MARGIN {
        Some(Containment::In)
    } else if winding.abs() <= WINDING_MARGIN {
        Some(Containment::Out)
    } else {
        None
    }
}

/// One face's mesh, as the winding number reads it.
#[derive(Debug)]
struct MeshedFace {
    triangles: Vec<[Point; 3]>,
    bound: Aabb,
    /// For a plane face, a point of its plane and the plane's unit normal.
    plane: Option<(Point, Vector)>,
    /// How far off the surface the face's own vertices and edges may lie.
    reach: f64,
}

/// A solid's faces meshed at one chord, for winding numbers.
#[derive(Debug)]
pub(crate) struct WindingMesh {
    faces: Vec<MeshedFace>,
    chord: f64,
}

/// What a point's winding number against a [`WindingMesh`] gives.
pub(crate) enum Wound {
    /// The point stands clear of every face's band and reads near an
    /// integer it can answer with.
    Reads(Containment),
    /// The point stands clear but the number is not near zero or one.
    Unread,
    /// The point stands within a curved face's band, this far from its
    /// mesh: a finer chord may read it.
    Near(f64),
}

impl WindingMesh {
    /// Mesh every face of `solid` at `chord`, edges agreed between faces.
    ///
    /// # Errors
    ///
    /// Where a face will not mesh.
    pub(crate) fn of(
        model: &Model,
        solid: &Shape,
        chord: f64,
        kept: Option<&FaceMeshCache>,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let deflection = Deflection::with_chord(chord)?;
        let (meshes, _) = face_meshes_for(model, solid, deflection, kept, tol)?;
        let mut faces = Vec::with_capacity(meshes.len());
        for (face, mesh) in meshes {
            let mesh = mesh?;
            let triangles: Vec<[Point; 3]> = mesh
                .triangles
                .iter()
                .map(|t| t.map(|i| mesh.positions[i as usize]))
                .collect();
            let bound = Aabb::of_points(&mesh.positions);
            let (reach, flat) = match model.node(&face).map(|n| n.data()) {
                Some(NodeData::Face(data)) => (
                    data.tolerance.get().max(tol.confusion()),
                    matches!(
                        model.geometry().surface(data.surface),
                        Some(ogeom_geom::SurfaceGeometry::Plane(_))
                    ),
                ),
                _ => (tol.confusion(), false),
            };
            let plane = if flat { plane_of(&triangles) } else { None };
            faces.push(MeshedFace {
                triangles,
                bound,
                plane,
                reach,
            });
        }
        Ok(Self { faces, chord })
    }

    /// The chord the faces were meshed at.
    pub(crate) const fn chord(&self) -> f64 {
        self.chord
    }

    /// `p`'s winding number against the faces, where `p` stands clear of
    /// the band between each face and its mesh.
    pub(crate) fn wind(&self, p: Point) -> Wound {
        let mut nearest = f64::INFINITY;
        let mut sum = 0.0;
        for face in &self.faces {
            // The band about a curved face's mesh, and about a plane face's
            // within its plane.
            let band = 2.0 * self.chord + face.reach;
            let clear = face
                .plane
                .is_some_and(|(o, n)| (p - o).dot(n).abs() > face.reach)
                || face.bound.distance_to(p) > band
                || {
                    let d = face
                        .triangles
                        .iter()
                        .map(|t| crate::classify::distance_to_triangle(p, *t))
                        .fold(f64::INFINITY, f64::min);
                    if d <= band {
                        nearest = nearest.min(d);
                    }
                    d > band
                };
            if clear {
                sum += face
                    .triangles
                    .iter()
                    .map(|t| solid_angle(p, *t))
                    .sum::<f64>();
            }
        }
        if nearest.is_finite() {
            return Wound::Near(nearest);
        }
        reading(sum / (4.0 * core::f64::consts::PI)).map_or(Wound::Unread, Wound::Reads)
    }
}

/// The plane of a plane face's mesh: its largest triangle's.
fn plane_of(triangles: &[[Point; 3]]) -> Option<(Point, Vector)> {
    let (t, normal) = triangles
        .iter()
        .map(|t| (t, (t[1] - t[0]).cross(t[2] - t[0])))
        .max_by(|a, b| a.1.magnitude().total_cmp(&b.1.magnitude()))?;
    let length = normal.magnitude();
    (length > 0.0).then(|| (t[0], normal / length))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tetrahedron_winds_once_about_its_inside_and_not_about_its_outside() {
        let [a, b, c, d] = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        // Wound outward.
        let tetrahedron = [[a, c, b], [a, b, d], [a, d, c], [b, c, d]];
        let inside = winding_number(Point::new(0.2, 0.2, 0.2), tetrahedron.iter());
        let outside = winding_number(Point::new(1.0, 1.0, 1.0), tetrahedron.iter());
        assert!((inside - 1.0).abs() < 1e-12, "{inside}");
        assert!(outside.abs() < 1e-12, "{outside}");
        // Without its slanted face, the shell reads neither.
        let open = winding_number(Point::new(0.2, 0.2, 0.2), tetrahedron[..3].iter());
        assert_eq!(reading(open), None, "{open}");
    }
}
