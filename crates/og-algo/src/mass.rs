//! Mass properties: how much there is, where its centre is, and how it resists
//! being spun.
//!
//! Three measures, one for each dimension a shape can have — the length of its
//! edges, the area of its faces, the volume it encloses — each with the centre
//! of that measure and the inertia tensor about that centre.
//!
//! # Computed from the tessellation, and why
//!
//! A polyhedron's volume has a closed form and a cylinder's does not — not
//! once it is trimmed by arbitrary wires. An implementation exact for planar
//! faces would be right for a box, right for a wedge, and quietly wrong for
//! anything curved, with nothing in the answer to say which case it was. So
//! every measure here comes from the same tessellation, and every result
//! carries the deflection it was computed at.
//!
//! That makes the error bounded and stated rather than hidden. Halving the
//! deflection and seeing the answer move tells a caller exactly how much to
//! trust it — [`MassProperties::deflection`] is what makes that check possible.
//!
//! # The one formula
//!
//! Length, area and volume all reduce to summing over simplices — segments,
//! triangles, tetrahedra — and the second moment of a simplex has the same
//! shape in every dimension:
//!
//! ```text
//! ∫ x_i x_j  =  m / (n(n+1)) · [ Σ_k p_k p_kᵀ + (Σ_k p_k)(Σ_k p_k)ᵀ ]
//! ```
//!
//! for `n` vertices and measure `m`. Barycentric integration gives it: the
//! integral of `λ_a λ_b` over a simplex is `m·d!·(1+δ_ab)/(d+2)!`, and
//! `n(n+1)` is what that collapses to. One function serves all three, which is
//! also why the three agree with each other rather than drifting apart.

use og_core::{OgResult, Tolerances, og_bail};
use og_math::{Direction, Matrix3, Point, Vector};
use og_mesh::{Deflection, discretize};
use og_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

/// How much of something there is, and how it is distributed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// The measure: length, area or volume, depending on what was asked for.
    ///
    /// Never negative. A volume computed from an inward-wound shell would come
    /// out negative, which says the shell is inside out rather than that the
    /// solid has negative volume, so that case is an error instead.
    pub mass: f64,
    /// The centre of the measure — the centroid, or centre of mass at uniform
    /// density.
    pub centre: Point,
    /// The inertia tensor about [`MassProperties::centre`], at unit density.
    ///
    /// About the centre, not the origin: an inertia about the origin says as
    /// much about where the part happens to sit as about the part.
    /// [`MassProperties::inertia_about`] moves it elsewhere.
    pub inertia: Matrix3,
    /// The chord deflection the tessellation was built to.
    ///
    /// The honest statement of accuracy. For a shape with only planar faces
    /// and straight edges the result is exact whatever this says, because the
    /// tessellation is exact.
    pub deflection: f64,
}

impl MassProperties {
    /// Nothing: no mass, at the origin, resisting nothing.
    #[must_use]
    pub const fn none(deflection: f64) -> Self {
        Self {
            mass: 0.0,
            centre: Point::ORIGIN,
            inertia: Matrix3::ZERO,
            deflection,
        }
    }

    /// The inertia tensor about some other point, by the parallel axis theorem.
    #[must_use]
    pub fn inertia_about(&self, point: Point) -> Matrix3 {
        let d = self.centre - point;
        // Moving *away* from the centre can only increase inertia, which is the
        // sign convention here: the centre is the minimum.
        add(self.inertia, displacement_term(self.mass, d))
    }

    /// The radius of gyration about an axis through the centre.
    ///
    /// The distance at which a point of the same mass would have the same
    /// inertia. Zero mass has no such distance, so this returns `None` rather
    /// than dividing by it.
    #[must_use]
    pub fn radius_of_gyration(&self, axis: Direction) -> Option<f64> {
        if self.mass <= 0.0 {
            return None;
        }
        let v = axis.vector();
        let i = quadratic_form(self.inertia, v);
        Some((i / self.mass).max(0.0).sqrt())
    }

    /// The principal moments, smallest first, with the axes they act about.
    ///
    /// The eigenvectors of a symmetric tensor, so the axes are orthogonal. A
    /// shape with rotational symmetry has repeated moments and the axes in that
    /// plane are arbitrary but still orthogonal — which is correct, not a
    /// failure: any pair of perpendicular axes in that plane is principal.
    ///
    /// # Errors
    ///
    /// [`OgError::NotDone`](og_core::OgError::NotDone) if the eigen-solver does
    /// not converge, which for a symmetric 3×3 means the tensor was not finite.
    pub fn principal_axes(&self, tol: Tolerances) -> OgResult<[(f64, Direction); 3]> {
        let m = nalgebra::Matrix3::from_row_slice(&[
            self.inertia.rows[0][0],
            self.inertia.rows[0][1],
            self.inertia.rows[0][2],
            self.inertia.rows[1][0],
            self.inertia.rows[1][1],
            self.inertia.rows[1][2],
            self.inertia.rows[2][0],
            self.inertia.rows[2][1],
            self.inertia.rows[2][2],
        ]);
        if !m.iter().all(|x| x.is_finite()) {
            og_bail!(NotDone, "the inertia tensor is not finite");
        }
        // Symmetric by construction, so the eigenvalues are real and this
        // always converges; the general solver would return complex pairs.
        let eigen = nalgebra::SymmetricEigen::new(m);

        let mut out: Vec<(f64, Direction)> = Vec::with_capacity(3);
        for i in 0..3 {
            let column = eigen.eigenvectors.column(i);
            let axis = Direction::new(Vector::new(column[0], column[1], column[2]), tol)?;
            out.push((eigen.eigenvalues[i], axis));
        }
        out.sort_by(|a, b| a.0.total_cmp(&b.0));
        Ok([out[0], out[1], out[2]])
    }
}

/// The length of a shape's edges, and how it is distributed.
///
/// Every distinct edge counts once, however many faces it bounds — the wire
/// frame of the shape, not a tally weighted by use.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the deflection
/// settings are unusable, or a curve is missing from the model.
pub fn linear_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<MassProperties> {
    deflection.validate()?;
    let mut acc = Accumulator::new();

    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(node) = model.node(&edge) else {
            og_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            og_bail!(Construction, "edge node holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            og_bail!(Dangling, "curve is not in this model");
        };
        let placement = edge.transform(model.datums())?;
        let line = discretize(geometry, *range, deflection, tol)?;
        for w in line.points.windows(2) {
            let (a, b) = (placement.apply(w[0]), placement.apply(w[1]));
            acc.add(&[a, b], a.distance(b));
        }
    }
    Ok(acc.finish(deflection.chord))
}

/// The area of a shape's faces, and how it is distributed.
///
/// # Errors
///
/// As [`og_mesh::triangulate_face`].
pub fn surface_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<MassProperties> {
    deflection.validate()?;
    let mut acc = Accumulator::new();

    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let mesh = og_mesh::triangulate_face(model, &face, deflection, tol)?;
        for triangle in &mesh.triangles {
            let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
            // Unsigned: a reversed face still has the same area, and summing
            // signed areas would cancel a solid's own surface to nothing.
            let area = (b - a).cross(c - a).magnitude() * 0.5;
            acc.add(&[a, b, c], area);
        }
    }
    Ok(acc.finish(deflection.chord))
}

/// The volume a shape encloses, and how it is distributed.
///
/// # Errors
///
/// As [`surface_properties`], plus
/// [`OgError::Construction`](og_core::OgError::Construction) if the boundary is
/// not closed, or is wound inward so that the volume comes out negative. Both
/// mean the answer would be meaningless rather than merely inaccurate: the
/// divergence theorem needs a closed, outward-oriented boundary, and without
/// one the sum is a number with no interpretation.
pub fn volume_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<MassProperties> {
    deflection.validate()?;
    let mesh = og_mesh::triangulate(model, shape, deflection, tol)?;
    if mesh.is_empty() {
        return Ok(MassProperties::none(deflection.chord));
    }
    if !mesh.is_closed() {
        og_bail!(
            Construction,
            "the boundary is not closed, so it encloses no volume to measure"
        );
    }

    // The apex every tetrahedron is built on. Any point serves — the signs
    // cancel outside the enclosed region wherever it sits — so it is a point on
    // the mesh, which keeps the tetrahedra the size of the shape instead of the
    // size of its distance from the world origin.
    let apex = mesh.positions[0];
    let mut acc = Accumulator::new();
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
        // The signed volume of the tetrahedron on the apex. The cancellation is
        // the divergence theorem doing the work, and why the winding has to be
        // outward.
        let volume = (a - apex).dot((b - apex).cross(c - apex)) / 6.0;
        acc.add(&[apex, a, b, c], volume);
    }

    if acc.mass < 0.0 {
        og_bail!(
            Construction,
            "the boundary is wound inward, so the volume came out negative"
        );
    }
    Ok(acc.finish(deflection.chord))
}

/// Running totals over simplices, measured from a fixed reference point.
///
/// The moments are accumulated about a *reference near the shape*, not about
/// the world origin, and that is a numerical decision rather than a stylistic
/// one. The inertia about the centre is a difference of two second moments, so
/// for a part sitting a million units from the origin the two terms agree to
/// twelve digits and their difference keeps four. Referencing the shape's own
/// bounding box keeps every intermediate the size of the shape.
///
/// The moments are about a fixed point rather than a running centre because the
/// centre is not known until the last simplex is in, and a moment about a
/// moving point is not a sum of anything.
struct Accumulator {
    /// Where the moments are measured from: the first point seen, so it is
    /// always somewhere on the shape.
    reference: Option<Point>,
    mass: f64,
    /// The first moment `∫ (x − r)`, which divided by the mass gives the centre
    /// relative to the reference.
    first: Vector,
    /// The second moment `∫ (x − r)(x − r)ᵀ`.
    second: Matrix3,
}

impl Accumulator {
    const fn new() -> Self {
        Self {
            reference: None,
            mass: 0.0,
            first: Vector::ZERO,
            second: Matrix3::ZERO,
        }
    }

    /// Add one simplex: 2 points for a segment, 3 for a triangle, 4 for a
    /// tetrahedron, with its signed or unsigned measure.
    fn add(&mut self, points: &[Point], measure: f64) {
        if points.is_empty() || measure == 0.0 || !measure.is_finite() {
            return;
        }
        let n = points.len();
        #[allow(clippy::cast_precision_loss)]
        let count = n as f64;
        let reference = *self.reference.get_or_insert(points[0]);
        let local: Vec<Vector> = points.iter().map(|p| *p - reference).collect();
        let sum: Vector = local.iter().fold(Vector::ZERO, |a, v| a + *v);

        self.mass += measure;
        self.first += sum * (measure / count);

        // ∫ x_i x_j = m/(n(n+1)) · [ Σ p p_ᵀ + (Σ p)(Σ p)ᵀ ] — see the module
        // docs. The n(n+1) is the barycentric integral collapsing.
        let scale = measure / (count * (count + 1.0));
        let mut term = outer(sum, sum);
        for v in &local {
            term = add(term, outer(*v, *v));
        }
        self.second = add(self.second, scale_matrix(term, scale));
    }

    /// Turn the running totals into the answer.
    fn finish(self, deflection: f64) -> MassProperties {
        if self.mass.abs() <= f64::MIN_POSITIVE {
            return MassProperties::none(deflection);
        }
        let offset = self.first / self.mass;
        let centre = self.reference.unwrap_or(Point::ORIGIN) + offset;

        // Inertia about the reference from the second moment: I = tr(S)·1 − S.
        let trace = self.second.rows[0][0] + self.second.rows[1][1] + self.second.rows[2][2];
        let about_reference = add(
            scale_matrix(Matrix3::IDENTITY, trace),
            scale_matrix(self.second, -1.0),
        );
        // Then shift to the centre — the reverse of `inertia_about`.
        let inertia = add(
            about_reference,
            scale_matrix(displacement_term(self.mass, offset), -1.0),
        );

        MassProperties {
            mass: self.mass.abs(),
            centre,
            inertia,
            deflection,
        }
    }
}

/// The parallel-axis contribution of a mass displaced by `d`.
fn displacement_term(mass: f64, d: Vector) -> Matrix3 {
    let squared = d.dot(d);
    add(
        scale_matrix(Matrix3::IDENTITY, mass * squared),
        scale_matrix(outer(d, d), -mass),
    )
}

/// The outer product `a bᵀ`.
fn outer(a: Vector, b: Vector) -> Matrix3 {
    Matrix3::new([
        [a.x * b.x, a.x * b.y, a.x * b.z],
        [a.y * b.x, a.y * b.y, a.y * b.z],
        [a.z * b.x, a.z * b.y, a.z * b.z],
    ])
}

/// Element-wise sum.
fn add(a: Matrix3, b: Matrix3) -> Matrix3 {
    let mut rows = a.rows;
    for (row, other) in rows.iter_mut().zip(b.rows) {
        for (value, addend) in row.iter_mut().zip(other) {
            *value += addend;
        }
    }
    Matrix3::new(rows)
}

/// Element-wise scaling.
fn scale_matrix(m: Matrix3, s: f64) -> Matrix3 {
    let mut rows = m.rows;
    for row in &mut rows {
        for value in row {
            *value *= s;
        }
    }
    Matrix3::new(rows)
}

/// `vᵀ M v`.
fn quadratic_form(m: Matrix3, v: Vector) -> f64 {
    let c = [v.x, v.y, v.z];
    let mut sum = 0.0;
    for (i, ci) in c.iter().enumerate() {
        for (j, cj) in c.iter().enumerate() {
            sum += ci * m.rows[i][j] * cj;
        }
    }
    sum
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::make_box;
    use approx::assert_relative_eq;
    use og_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_box_has_the_volume_centre_and_inertia_a_box_has() {
        // Every number here is one a textbook states, which is the point: the
        // simplex formula is general, and a general formula that gets the one
        // case everybody knows wrong is worth nothing.
        let (dx, dy, dz) = (2.0, 3.0, 4.0);
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (dx, dy, dz), T).unwrap();

        let props = volume_properties(&model, &built.shape, fine(), T).unwrap();
        assert_relative_eq!(props.mass, dx * dy * dz, epsilon = 1e-9);
        assert!(
            props
                .centre
                .is_equal(Point::new(dx / 2.0, dy / 2.0, dz / 2.0), T),
            "the centre of a box is its middle, got {:?}",
            props.centre
        );

        // I_xx = m(dy² + dz²)/12, and so round.
        let m = dx * dy * dz;
        assert_relative_eq!(
            props.inertia.rows[0][0],
            m * dz.mul_add(dz, dy * dy) / 12.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            props.inertia.rows[1][1],
            m * dz.mul_add(dz, dx * dx) / 12.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            props.inertia.rows[2][2],
            m * dy.mul_add(dy, dx * dx) / 12.0,
            epsilon = 1e-9
        );
        // A box is symmetric about its own axes, so the products vanish.
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            assert_relative_eq!(props.inertia.rows[i][j], 0.0, epsilon = 1e-9);
            assert_relative_eq!(props.inertia.rows[j][i], 0.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn a_box_has_the_area_and_edge_length_a_box_has() {
        let (dx, dy, dz) = (2.0, 3.0, 4.0);
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (dx, dy, dz), T).unwrap();

        let area = surface_properties(&model, &built.shape, fine(), T).unwrap();
        assert_relative_eq!(
            area.mass,
            2.0 * dz.mul_add(dx, dx.mul_add(dy, dy * dz)),
            epsilon = 1e-9
        );
        assert!(
            area.centre
                .is_equal(Point::new(dx / 2.0, dy / 2.0, dz / 2.0), T)
        );

        // Four edges in each direction, counted once each however many faces
        // they bound.
        let length = linear_properties(&model, &built.shape, fine(), T).unwrap();
        assert_relative_eq!(length.mass, 4.0 * (dx + dy + dz), epsilon = 1e-9);
        assert!(
            length
                .centre
                .is_equal(Point::new(dx / 2.0, dy / 2.0, dz / 2.0), T)
        );
    }

    #[test]
    fn the_answer_does_not_depend_on_where_the_shape_sits() {
        // The inertia is about the centre, so translating the box must leave it
        // alone and move only the centre. An inertia accidentally left about
        // the origin would grow with the distance.
        let mut model = Model::new();
        let here = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T).unwrap();
        let far = Frame::new(
            Point::new(100.0, -50.0, 25.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let there = make_box(&mut model, far, (1.0, 2.0, 3.0), T).unwrap();

        let a = volume_properties(&model, &here.shape, fine(), T).unwrap();
        let b = volume_properties(&model, &there.shape, fine(), T).unwrap();

        assert_relative_eq!(a.mass, b.mass, epsilon = 1e-9);
        assert!(
            b.centre
                .is_equal(a.centre + Vector::new(100.0, -50.0, 25.0), T)
        );
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(a.inertia.rows[i][j], b.inertia.rows[i][j], epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn a_part_a_long_way_from_the_origin_keeps_its_precision() {
        // The reason the moments are accumulated about a point on the shape.
        // About the world origin the two terms of the inertia agree to twelve
        // digits at this distance and their difference keeps four — the answer
        // would come back with a few percent of noise in it, or negative.
        let mut model = Model::new();
        let far = Frame::new(
            Point::new(1.0e6, -2.0e6, 5.0e5),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let built = make_box(&mut model, far, (2.0, 3.0, 4.0), T).unwrap();
        let props = volume_properties(&model, &built.shape, fine(), T).unwrap();

        assert_relative_eq!(props.mass, 24.0, epsilon = 1e-6);
        assert_relative_eq!(
            props.inertia.rows[0][0],
            24.0 * 4.0_f64.mul_add(4.0, 3.0 * 3.0) / 12.0,
            epsilon = 1e-6
        );
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            assert_relative_eq!(props.inertia.rows[i][j], 0.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn moving_the_inertia_off_the_centre_agrees_with_the_parallel_axis_theorem() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let props = volume_properties(&model, &built.shape, fine(), T).unwrap();

        // A cube of side a about a face-centre axis: I = m(a²/6 + a²/4).
        let m = 8.0;
        let corner = props.inertia_about(Point::ORIGIN);
        assert_relative_eq!(
            corner.rows[0][0],
            2.0_f64.mul_add(2.0, 2.0 * 2.0).mul_add(m / 12.0, m * 2.0),
            epsilon = 1e-9
        );
        // And about its own centre it is the smallest it can be.
        assert!(corner.rows[0][0] > props.inertia.rows[0][0]);
    }

    #[test]
    fn a_cubes_principal_moments_are_all_the_same() {
        // Full rotational symmetry: every axis is principal, so the three
        // moments must agree. Axes that came back non-orthogonal would mean the
        // solver was handed a non-symmetric tensor, which would itself be a bug.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let props = volume_properties(&model, &built.shape, fine(), T).unwrap();

        let axes = props.principal_axes(T).unwrap();
        let expected = 8.0 * 2.0_f64.mul_add(2.0, 2.0 * 2.0) / 12.0;
        for (moment, _) in &axes {
            assert_relative_eq!(*moment, expected, epsilon = 1e-6);
        }
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            assert_relative_eq!(
                axes[i].1.vector().dot(axes[j].1.vector()),
                0.0,
                epsilon = 1e-9
            );
        }
    }

    #[test]
    fn a_long_box_spins_most_easily_about_its_length() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (10.0, 1.0, 1.0), T).unwrap();
        let props = volume_properties(&model, &built.shape, fine(), T).unwrap();

        let axes = props.principal_axes(T).unwrap();
        // The smallest moment is about the long axis.
        assert!(axes[0].1.vector().x.abs() > 0.99, "got {:?}", axes[0].1);
        assert!(axes[0].0 < axes[1].0 && axes[1].0 <= axes[2].0);

        let along = props.radius_of_gyration(Direction::X).unwrap();
        let across = props.radius_of_gyration(Direction::Y).unwrap();
        assert!(along < across, "{along} should be less than {across}");
    }

    #[test]
    fn a_sphere_converges_on_the_volume_a_sphere_has() {
        // The case a planar-exact implementation would get quietly wrong. The
        // tessellation inscribes the sphere, so the volume comes in under the
        // truth and climbs as the deflection tightens — and the deflection is
        // reported, so a caller can see how far under.
        use crate::build::make_natural_face;
        use og_geom::SphereSurface;
        use og_math::Sphere;

        let radius = 5.0_f64;
        let exact = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
        let mut previous = 0.0;

        for chord in [0.5_f64, 0.1, 0.02] {
            let mut model = Model::new();
            let surface = SphereSurface::new(Sphere::new(Frame::WORLD, radius, T).unwrap());
            let face = make_natural_face(&mut model, surface.into()).unwrap().shape;
            let shell = crate::build::make_shell(&mut model, std::slice::from_ref(&face))
                .unwrap()
                .shape;

            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            let props = volume_properties(&model, &shell, deflection, T).unwrap();
            assert_relative_eq!(props.deflection, chord);
            assert!(props.mass < exact, "an inscribed volume cannot exceed it");
            assert!(
                props.mass > previous,
                "tightening the chord lost volume: {} after {previous}",
                props.mass
            );
            assert!(
                props
                    .centre
                    .is_equal(Point::ORIGIN, Tolerances::with_scale(1e4).unwrap()),
                "a sphere's centre is its centre, got {:?}",
                props.centre
            );
            previous = props.mass;
        }
        assert!(
            previous > exact * 0.99,
            "{previous} should be within a percent of {exact}"
        );
    }

    #[test]
    fn an_open_shell_is_refused_rather_than_measured() {
        // Half a boundary encloses nothing, and the divergence theorem applied
        // to it returns a number that looks like a volume and is not one.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();

        assert!(volume_properties(&model, &face, fine(), T).is_err());
        // Its area, though, is perfectly well defined.
        let area = surface_properties(&model, &face, fine(), T).unwrap();
        assert_relative_eq!(area.mass, 1.0, epsilon = 1e-9);
    }

    #[test]
    fn an_inward_shell_is_refused_rather_than_reported_as_negative() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        assert!(volume_properties(&model, &built.shape.reversed(), fine(), T).is_err());
    }

    #[test]
    fn a_shape_with_nothing_to_measure_says_so_rather_than_dividing_by_zero() {
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);

        for props in [
            volume_properties(&model, &vertex, fine(), T).unwrap(),
            surface_properties(&model, &vertex, fine(), T).unwrap(),
            linear_properties(&model, &vertex, fine(), T).unwrap(),
        ] {
            assert_relative_eq!(props.mass, 0.0);
            assert!(props.centre.is_equal(Point::ORIGIN, T));
            assert!(props.radius_of_gyration(Direction::Z).is_none());
        }
    }

    #[test]
    fn an_unusable_deflection_is_refused() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let bad = Deflection {
            chord: -1.0,
            ..Deflection::default()
        };
        assert!(volume_properties(&model, &built.shape, bad, T).is_err());
        assert!(surface_properties(&model, &built.shape, bad, T).is_err());
        assert!(linear_properties(&model, &built.shape, bad, T).is_err());
    }
}
