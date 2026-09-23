//! Mass properties: how much there is, where its centre is, and how it resists
//! being spun.
//!
//! Three measures, one for each dimension a shape can have (the length of its
//! edges, the area of its faces, the volume it encloses), each with the centre
//! of that measure and the inertia tensor about that centre.
//!
//! # Computed from the tessellation, and why
//!
//! A polyhedron's volume has a closed form and a cylinder's does not, not
//! once it is trimmed by arbitrary wires. An implementation exact for planar
//! faces would be right for a box, right for a wedge, and quietly wrong for
//! anything curved, with nothing in the answer to say which case it was. So
//! every measure here comes from the same tessellation, and every result
//! carries the deflection it was computed at.
//!
//! That makes the error bounded and stated rather than hidden. Halving the
//! deflection and seeing the answer move tells a caller exactly how much to
//! trust it; [`MassProperties::deflection`] is what makes that check possible.
//!
//! # The one formula
//!
//! Length, area and volume all reduce to summing over simplices (segments,
//! triangles, tetrahedra), and the second moment of a simplex has the same
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

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Transformable as _;
use ogeom_math::{Direction, Matrix3, Point, Vector};
use ogeom_mesh::{Deflection, discretize};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

/// How much of something there is, and how it is distributed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// The measure: length, area or volume, depending on what was asked for.
    ///
    /// Never negative. A volume computed from an inward-wound shell would come
    /// out negative, which says the shell is inside out rather than that the
    /// solid has negative volume, so that case is an error instead.
    pub mass: f64,
    /// The centre of the measure: the centroid, or centre of mass at uniform
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
    /// plane are arbitrary but still orthogonal, which is correct, not a
    /// failure: any pair of perpendicular axes in that plane is principal.
    ///
    /// # Errors
    ///
    /// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the eigen-solver does
    /// not converge, which for a symmetric 3×3 means the tensor was not finite.
    pub fn principal_axes(&self, tol: Tolerances) -> OgeomResult<[(f64, Direction); 3]> {
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
            ogeom_bail!(NotDone, "the inertia tensor is not finite");
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
/// Every distinct edge counts once, however many faces it bounds: the wire
/// frame of the shape, not a tally weighted by use.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the deflection
/// settings are unusable, or a curve is missing from the model.
pub fn linear_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<MassProperties> {
    deflection.validate()?;
    let mut acc = Accumulator::new();

    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(node) = model.node(&edge) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
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
/// As [`ogeom_mesh::triangulate_face`].
pub fn surface_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<MassProperties> {
    deflection.validate()?;
    if let Some(exact) = exact_surface_properties(model, shape, tol)? {
        return Ok(exact);
    }
    let mut acc = Accumulator::new();

    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let mesh = ogeom_mesh::triangulate_face(model, &face, deflection, tol)?;
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the boundary is
/// not closed, or is wound inward so that the volume comes out negative. Both
/// mean the answer would be meaningless rather than merely inaccurate: the
/// divergence theorem needs a closed, outward-oriented boundary, and without
/// one the sum is a number with no interpretation.
pub fn volume_properties(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<MassProperties> {
    deflection.validate()?;
    if let Some(exact) = exact_volume_properties(model, shape, tol)? {
        return Ok(exact);
    }
    let mesh = ogeom_mesh::triangulate(model, shape, deflection, tol)?;
    if mesh.is_empty() {
        return Ok(MassProperties::none(deflection.chord));
    }
    if !mesh.is_closed() {
        ogeom_bail!(
            Construction,
            "the boundary is not closed, so it encloses no volume to measure"
        );
    }

    // The apex every tetrahedron is built on. Any point serves (the signs
    // cancel outside the enclosed region wherever it sits), so it is a point on
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
        ogeom_bail!(
            Construction,
            "the boundary is wound inward, so the volume came out negative"
        );
    }
    Ok(acc.finish(deflection.chord))
}

// --- exact properties on the exact surfaces ----------------------------------

/// A face whose trim the exact integrator can walk: an analytic surface
/// trimmed to a chart rectangle, or a plane trimmed to a full disc.
enum ExactFace {
    /// `[u0, u1] x [v0, v1]` on the (placed) surface.
    ChartRectangle {
        surface: ogeom_geom::SurfaceGeometry,
        rect: (f64, f64, f64, f64),
        sign: f64,
        share: f64,
    },
    /// A full circular disc on a plane.
    Disc {
        centre: Point,
        e1: Vector,
        e2: Vector,
        normal: Vector,
        radius: f64,
        sign: f64,
        share: f64,
    },
}

impl ExactFace {
    /// Whether this region is part of its face or taken out of it: `1` for
    /// the outer boundary, `-1` for a hole. The volume integral could carry
    /// it in the normal's sign, but the area integral takes a magnitude and
    /// would hand a hole's area back as more face.
    const fn share(&self) -> f64 {
        match self {
            Self::ChartRectangle { share, .. } | Self::Disc { share, .. } => *share,
        }
    }

    /// How much of its chart the region covers, for telling a face's outer
    /// boundary from its holes. A wire's place in the face's list does not
    /// say which it is (a ring's annulus arrives inner ring first), and a
    /// hole is inside the boundary it is a hole in, so it covers less.
    fn chart_area(&self) -> f64 {
        match self {
            Self::Disc { radius, .. } => core::f64::consts::PI * radius * radius,
            Self::ChartRectangle { rect, .. } => (rect.1 - rect.0) * (rect.3 - rect.2),
        }
    }

    fn take_away(&mut self) {
        match self {
            Self::ChartRectangle { share, .. } | Self::Disc { share, .. } => *share = -1.0,
        }
    }
}

/// Mass properties integrated on the exact surfaces, when every face allows.
///
/// The integrands over an analytic surface's chart are trigonometric
/// polynomials, and panels no wider than a quarter turn under the ten-point
/// Gauss rule integrate them to rounding, exact in every sense that
/// matters, with `deflection` reported as zero. The first face that resists
/// (a non-analytic surface, a trim that is not a chart rectangle or a disc)
/// returns `None`, and the caller falls back to the tessellation with its
/// stated chord.
fn exact_volume_properties(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<MassProperties>> {
    let faces = explore(model, shape, Filter::OfType(ShapeType::Face))?;
    if faces.is_empty() {
        return Ok(None);
    }
    let mut exact = Vec::with_capacity(faces.len());
    for face in &faces {
        match exact_face(model, face, tol)? {
            Some(found) => exact.extend(found),
            None => {
                if std::env::var_os("OGEOM_DEBUG_MASS").is_some() {
                    eprintln!(
                        "MASS face {} is not exactly integrable",
                        face.node().index()
                    );
                }
                return Ok(None);
            }
        }
    }
    // And the faces must agree with each other about which way is out.
    if !flags_agree(model, shape, tol)? {
        return Ok(None);
    }
    // The divergence theorem needs a closed boundary; topology says whether
    // it has one. A shape with no shell at all (a bare face) has nothing
    // to close, and falls back to the mesh path, which refuses it properly.
    let shells = explore_unique(model, shape, ShapeType::Shell)?;
    if shells.is_empty() {
        return Ok(None);
    }
    for shell in shells {
        if !crate::build::is_shell_closed(model, &shell)? {
            ogeom_bail!(
                Construction,
                "the boundary is not closed, so it encloses no volume to measure"
            );
        }
    }

    let reference = reference_point(&exact, tol)?;
    let mut mass = 0.0;
    let mut first = Vector::ZERO;
    let mut second = Matrix3::ZERO;
    for face in &exact {
        integrate_face(face, reference, tol, &mut |p, n_da, share| {
            let n_da = n_da * share;
            let q = p - reference;
            mass += q.dot(n_da) / 3.0;
            first += Vector::new(
                q.x * q.x * n_da.x / 2.0,
                q.y * q.y * n_da.y / 2.0,
                q.z * q.z * n_da.z / 2.0,
            );
            let d = [q.x, q.y, q.z];
            let nd = [n_da.x, n_da.y, n_da.z];
            for i in 0..3 {
                // Diagonal: int q_i^2 dV = surface int q_i^3 n_i / 3.
                second.rows[i][i] += d[i] * d[i] * d[i] * nd[i] / 3.0;
                // Off-diagonal: int q_i q_j dV = surface int q_i^2 q_j n_i / 2.
                for j in 0..3 {
                    if i != j {
                        second.rows[i][j] += d[i] * d[i] * d[j] * nd[i] / 2.0;
                    }
                }
            }
        })?;
    }
    // The off-diagonal identity fills each pair twice, once from each axis;
    // average them, which also symmetrizes rounding.
    for i in 0..3 {
        for j in (i + 1)..3 {
            let mean = f64::midpoint(second.rows[i][j], second.rows[j][i]);
            second.rows[i][j] = mean;
            second.rows[j][i] = mean;
        }
    }
    if mass < 0.0 {
        ogeom_bail!(
            Construction,
            "the boundary is wound inward, so the volume came out negative"
        );
    }
    let acc = Accumulator {
        reference: Some(reference),
        mass,
        first,
        second,
    };
    Ok(Some(acc.finish(0.0)))
}

/// Surface area and its distribution, on the exact surfaces.
fn exact_surface_properties(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<MassProperties>> {
    let faces = explore(model, shape, Filter::OfType(ShapeType::Face))?;
    if faces.is_empty() {
        return Ok(None);
    }
    let mut exact = Vec::with_capacity(faces.len());
    for face in &faces {
        match exact_face(model, face, tol)? {
            Some(found) => exact.extend(found),
            None => return Ok(None),
        }
    }
    let reference = reference_point(&exact, tol)?;
    let mut mass = 0.0;
    let mut first = Vector::ZERO;
    let mut second = Matrix3::ZERO;
    for face in &exact {
        integrate_face(face, reference, tol, &mut |p, n_da, share| {
            let da = n_da.magnitude() * share;
            let q = p - reference;
            mass += da;
            first += q * da;
            for (i, qi) in [q.x, q.y, q.z].iter().enumerate() {
                for (j, qj) in [q.x, q.y, q.z].iter().enumerate() {
                    second.rows[i][j] += qi * qj * da;
                }
            }
        })?;
    }
    let acc = Accumulator {
        reference: Some(reference),
        mass,
        first,
        second,
    };
    Ok(Some(acc.finish(0.0)))
}

/// Somewhere on the shape to measure moments from.
fn reference_point(faces: &[ExactFace], tol: Tolerances) -> OgeomResult<Point> {
    use ogeom_geom::Surface as _;
    match &faces[0] {
        ExactFace::ChartRectangle { surface, rect, .. } => surface.point_at(rect.0, rect.2, tol),
        ExactFace::Disc { centre, .. } => Ok(*centre),
    }
}

/// Drive the callback over every quadrature sample of a face.
///
/// The callback receives the world point and the outward-signed `n dA`
/// already weighted; summing the callback's contributions *is* the
/// integral.
fn integrate_face(
    face: &ExactFace,
    _reference: Point,
    tol: Tolerances,
    contribute: &mut dyn FnMut(Point, Vector, f64),
) -> OgeomResult<()> {
    let share = face.share();
    use ogeom_geom::Surface as _;
    const QUARTER: f64 = core::f64::consts::FRAC_PI_2;
    match face {
        ExactFace::ChartRectangle {
            surface,
            rect,
            sign,
            ..
        } => {
            let (u0, u1, v0, v1) = *rect;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let u_panels = (((u1 - u0) / QUARTER).ceil() as usize).max(1);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let v_panels = (((v1 - v0) / QUARTER).ceil() as usize).max(1);
            let mut failure = None;
            for iu in 0..u_panels {
                #[allow(clippy::cast_precision_loss)]
                let (ua, ub) = (
                    u0 + (u1 - u0) * iu as f64 / u_panels as f64,
                    u0 + (u1 - u0) * (iu + 1) as f64 / u_panels as f64,
                );
                for iv in 0..v_panels {
                    #[allow(clippy::cast_precision_loss)]
                    let (va, vb) = (
                        v0 + (v1 - v0) * iv as f64 / v_panels as f64,
                        v0 + (v1 - v0) * (iv + 1) as f64 / v_panels as f64,
                    );
                    // Nested Gauss with the callback fed directly: the outer
                    // integrand returns 0 and the samples carry the payload,
                    // with the weights recovered from unit integrands.
                    gauss2(ua, ub, va, vb, &mut |u, v, weight| {
                        if failure.is_some() {
                            return;
                        }
                        let sample = (|| -> OgeomResult<()> {
                            let p = surface.point_at(u, v, tol)?;
                            let (du, dv) = surface.d1_at(u, v, tol)?;
                            contribute(p, du.cross(dv) * (sign * weight), share);
                            Ok(())
                        })();
                        if let Err(e) = sample {
                            failure = Some(e);
                        }
                    });
                }
            }
            match failure {
                Some(e) => Err(e),
                None => Ok(()),
            }
        }
        ExactFace::Disc {
            centre,
            e1,
            e2,
            normal,
            radius,
            sign,
            ..
        } => {
            let failure: Option<ogeom_core::OgeomError> = None;
            let turns = 4;
            for k in 0..turns {
                #[allow(clippy::cast_precision_loss)]
                let (ta, tb) = (
                    core::f64::consts::TAU * k as f64 / turns as f64,
                    core::f64::consts::TAU * (k + 1) as f64 / turns as f64,
                );
                gauss2(0.0, *radius, ta, tb, &mut |rho, theta, weight| {
                    if failure.is_some() {
                        return;
                    }
                    let p = *centre + (*e1 * theta.cos() + *e2 * theta.sin()) * rho;
                    contribute(p, *normal * (sign * rho * weight), share);
                });
            }
            match failure {
                Some(e) => Err(e),
                None => Ok(()),
            }
        }
    }
}

/// A tensor-product ten-by-ten Gauss rule over `[a,b] x [c,d]`, feeding each
/// sample and its weight to the callback.
fn gauss2(a: f64, b: f64, c: f64, d: f64, f: &mut dyn FnMut(f64, f64, f64)) {
    // The rule's nodes recovered through the public one-dimensional
    // integrator: integrating a delta-free payload is not possible, so the
    // nodes are collected by integrating an indicator that records them.
    let mut us: Vec<(f64, f64)> = Vec::with_capacity(10);
    ogeom_math::gauss_legendre(
        |u| {
            us.push((u, 0.0));
            1.0
        },
        a,
        b,
    );
    // Weight of node i: integrate a basis that is 1 at that sample order.
    // Simpler: the rule is linear, so the weight is the integral of the
    // indicator sequence, recovered by a second pass per node.
    for (i, entry) in us.iter_mut().enumerate() {
        let mut k = 0;
        let w = ogeom_math::gauss_legendre(
            |_| {
                let value = if k == i { 1.0 } else { 0.0 };
                k += 1;
                value
            },
            a,
            b,
        );
        entry.1 = w;
    }
    let mut vs: Vec<(f64, f64)> = Vec::with_capacity(10);
    ogeom_math::gauss_legendre(
        |v| {
            vs.push((v, 0.0));
            1.0
        },
        c,
        d,
    );
    for (j, entry) in vs.iter_mut().enumerate() {
        let mut k = 0;
        let w = ogeom_math::gauss_legendre(
            |_| {
                let value = if k == j { 1.0 } else { 0.0 };
                k += 1;
                value
            },
            c,
            d,
        );
        entry.1 = w;
    }
    for &(u, wu) in &us {
        for &(v, wv) in &vs {
            f(u, v, wu * wv);
        }
    }
}

/// The exact-integrable regions of one face, or `None` where there are
/// none.
///
/// One region per wire, and the integral is their sum: a face's outer
/// boundary carries its own sign and every inner one the opposite, which
/// is what a hole *is* under the divergence theorem. So a plate with a
/// bore in it is a rectangle less a disc, and a tube's end face a disc
/// less a disc, neither of which had to be meshed, and both of which
/// were.
fn exact_face(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Vec<ExactFace>>> {
    let Some(node) = model.node(face) else {
        return Ok(None);
    };
    let NodeData::Face(data) = node.data() else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let analytic = matches!(
        surface,
        ogeom_geom::SurfaceGeometry::Plane(_)
            | ogeom_geom::SurfaceGeometry::Cylinder(_)
            | ogeom_geom::SurfaceGeometry::Cone(_)
            | ogeom_geom::SurfaceGeometry::Sphere(_)
            | ogeom_geom::SurfaceGeometry::Torus(_)
    );
    if !analytic {
        return Ok(None);
    }
    let placement = face.transform(model.datums())?;
    // The chart rectangle comes from the pcurves, whose windows are the
    // *unscaled* surface's; a scaling placement changes the chart's metric
    // and the windows with it, so only rigid placements take the exact path.
    if !matches!(
        placement.kind(),
        ogeom_math::TransformKind::Identity
            | ogeom_math::TransformKind::Translation
            | ogeom_math::TransformKind::Rotation
    ) {
        return Ok(None);
    }
    let placed = surface.clone().transformed(&placement, tol)?;
    let sign = if face.orientation() == ogeom_topo::Orientation::Reversed {
        -1.0
    } else {
        1.0
    };

    let wires = model.ordered_children_of(face)?;
    // One region per wire, and the integral is their sum: a face's outer
    // boundary carries its own sign and every inner one the opposite, which
    // is what a hole *is* under the divergence theorem. So a plate with a
    // bore is a rectangle less a disc, and a tube's end face a disc less a
    // disc. Which wire is the boundary and which the holes is settled by
    // the chart each covers: a hole is inside the boundary it is a hole in,
    // so it covers less.
    let mut regions = Vec::with_capacity(wires.len());
    for wire in &wires {
        let Some(region) = exact_wire(model, data, &placed, wire, sign, 1.0, tol)? else {
            return Ok(None);
        };
        regions.push(region);
    }
    let Some(outer) = (0..regions.len()).max_by(|a, b| {
        regions[*a]
            .chart_area()
            .total_cmp(&regions[*b].chart_area())
    }) else {
        return Ok(None);
    };
    for (index, region) in regions.iter_mut().enumerate() {
        if index != outer {
            region.take_away();
        }
    }
    Ok(Some(regions))
}

/// Whether the faces agree with each other about which way is out.
///
/// The flag on a face is the only thing that says which side of its surface
/// the material is on; no winding in this kernel says it, and the wires
/// are wound however their builder wound them. But the flags can be asked
/// *about each other*: an edge between two faces is walked by one of them
/// with the material on its left and by the other with the material on its
/// left too, so the two walks run opposite ways along it. Each face's walk
/// is its outward normal crossed into the direction the material lies from
/// the edge, and both of those are had for the asking: the normal from the
/// flag, the material's direction from the chart, since a face's region
/// lies around the middle of the boundary that encloses it.
///
/// A part in the corpus has a bore wall whose flag points into the solid.
/// The tessellator repairs such a shell, flipping whichever side of the
/// disagreement is in the minority, and the closed-form integral cannot:
/// it would hand the bore back as material, a third of that part's volume.
/// So where the flags disagree this says so and the mesh is asked instead.
///
/// An instanced solid says nothing here: one edge stands in several places
/// and nothing in a name tells them apart, so its flags are taken as they
/// come, which is what they were before there was anything to ask.
fn flags_agree(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let placed_at = ogeom_topo::Location::default();
    let mut walks: std::collections::HashMap<
        ogeom_topo::TShapeId,
        Vec<(bool, ogeom_topo::TShapeId, bool)>,
    > = std::collections::HashMap::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        if face.location() != &placed_at {
            return Ok(true);
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()).cloned() else {
            return Ok(true);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(true);
        };
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let flag = if face.orientation() == ogeom_topo::Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        // Each wire's middle in the chart, and which wire is the boundary:
        // the one covering the most of it, since a hole is inside what it
        // is a hole in.
        let wires = model.ordered_children_of(&face)?;
        let mut middles: Vec<(ogeom_math::Point2, f64)> = Vec::with_capacity(wires.len());
        let mut stations: Vec<Vec<(Shape, ogeom_math::Point2, f64)>> =
            Vec::with_capacity(wires.len());
        for wire in &wires {
            let mut here = Vec::new();
            let mut sum = ogeom_math::Vector2::new(0.0, 0.0);
            let (mut lo, mut hi) = (
                ogeom_math::Point2::new(f64::INFINITY, f64::INFINITY),
                ogeom_math::Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            );
            for edge in model.ordered_children_of(wire)? {
                if edge.location() != &placed_at {
                    return Ok(true);
                }
                let Some(repr) = model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .and_then(|d| d.pcurve_for(data.surface, edge.location()))
                else {
                    return Ok(true);
                };
                // A seam bounds its face twice, once down each column.
                let sides: Vec<(ogeom_topo::PCurveId, (f64, f64))> = match repr {
                    EdgeRepr::PCurve { curve, range, .. } => vec![(*curve, *range)],
                    EdgeRepr::Seam {
                        forward,
                        reversed,
                        range,
                        ..
                    } => vec![(*forward, *range), (*reversed, *range)],
                    _ => return Ok(true),
                };
                for (id, range) in sides {
                    let Some(pcurve) = model.geometry().pcurve(id) else {
                        return Ok(true);
                    };
                    // Several stations along each edge, not one: a wire of
                    // a single closed edge has its own midpoint for a
                    // middle, and nothing lies from a point toward itself.
                    const STATIONS: usize = 4;
                    for step in 1..=STATIONS {
                        #[allow(clippy::cast_precision_loss)]
                        let t =
                            range.0 + (range.1 - range.0) * (step as f64 / (STATIONS + 1) as f64);
                        let at = pcurve.point_at(t, tol)?;
                        sum += at.to_vector();
                        lo = ogeom_math::Point2::new(lo.x.min(at.x), lo.y.min(at.y));
                        hi = ogeom_math::Point2::new(hi.x.max(at.x), hi.y.max(at.y));
                        here.push((edge.clone(), at, t));
                    }
                }
            }
            if here.is_empty() {
                return Ok(true);
            }
            #[allow(clippy::cast_precision_loss)]
            let middle = ogeom_math::Point2::ORIGIN + sum / here.len() as f64;
            middles.push((middle, (hi.x - lo.x) * (hi.y - lo.y)));
            stations.push(here);
        }
        let Some(outer) = (0..middles.len()).max_by(|a, b| middles[*a].1.total_cmp(&middles[*b].1))
        else {
            return Ok(true);
        };
        for (index, here) in stations.into_iter().enumerate() {
            let (middle, _) = middles[index];
            for (edge, at, t) in here {
                let (du, dv) = placed.d1_at(at.x, at.y, tol)?;
                let raw = du.cross(dv);
                if raw.magnitude() <= tol.angular() {
                    return Ok(true);
                }
                let out = raw / raw.magnitude() * flag;
                // Which way the material lies from this point of the
                // boundary: toward the wire's middle for the face's outer
                // wire, away from it for a hole.
                let toward = middle - at;
                let toward = if index == outer { toward } else { -toward };
                let inward = du * toward.x + dv * toward.y;
                if inward.magnitude() <= tol.angular() {
                    return Ok(true);
                }
                let walk = out.cross(inward / inward.magnitude());
                // Against the edge's own direction, so the two faces'
                // answers can be compared without comparing vectors.
                let Some(along) = edge_direction(model, &edge, t, tol)? else {
                    return Ok(true);
                };
                walks.entry(edge.node()).or_default().push((
                    walk.dot(along) > 0.0,
                    face.node(),
                    index == outer,
                ));
            }
        }
    }
    if std::env::var_os("OGEOM_DEBUG_MASS").is_some() {
        eprintln!("MASS flags_agree walked {} edges", walks.len());
    }
    for (edge, uses) in &walks {
        // An edge one face walks twice is that face's own seam, however it
        // is written down (a canonicalised drum keeps its as an ordinary
        // pcurve used twice), and one face's seam says nothing about
        // whether two faces agree.
        if uses.iter().all(|(_, owner, _)| *owner == uses[0].1) {
            continue;
        }
        let ahead = uses.iter().filter(|(ahead, ..)| *ahead).count();
        if ahead * 2 != uses.len() {
            if std::env::var_os("OGEOM_DEBUG_MASS").is_some() {
                eprintln!("MASS edge {} is walked {uses:?}", edge.index());
            }
            return Ok(false);
        }
    }
    Ok(true)
}

/// An edge's own direction in space at the parameter `t` of its curve.
fn edge_direction(
    model: &Model,
    edge: &Shape,
    t: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    use ogeom_geom::Curve3d as _;
    use ogeom_geom::Transformable as _;
    let Some(curve) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            EdgeRepr::Curve3d { curve, .. } => Some(*curve),
            _ => None,
        })
        .and_then(|id| model.geometry().curve(id).cloned())
    else {
        return Ok(None);
    };
    let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
    let along = curve.d1_at(t, tol)?;
    Ok((along.magnitude() > tol.angular()).then(|| along / along.magnitude()))
}

/// The region one of a face's wires bounds, read off its pcurves.
fn exact_wire(
    model: &Model,
    data: &ogeom_topo::FaceData,
    placed: &ogeom_geom::SurfaceGeometry,
    wire: &Shape,
    sign: f64,
    share: f64,
    tol: Tolerances,
) -> OgeomResult<Option<ExactFace>> {
    use ogeom_geom::Surface as _;
    // Gather each boundary edge's chart segments on this face.
    let mut segments: Vec<(ogeom_math::Point2, ogeom_math::Point2)> = Vec::new();
    // Where a seam bounds the chart: a column at that `u`, a row at that `v`.
    let mut columns: Vec<f64> = Vec::new();
    let mut rows: Vec<f64> = Vec::new();
    let mut circle: Option<(ogeom_geom::Circle2d, f64)> = None;
    // Where the circle's arcs start and stop, for asking whether they tile
    // its turn or merely add up to one.
    let mut arc_ends: Vec<ogeom_math::Point2> = Vec::new();
    let mut pieces = 0_usize;
    // A seam bounds the face twice; its two chart sides are gathered once.
    let mut seams_seen: Vec<ogeom_topo::TShapeId> = Vec::new();
    for edge in model.ordered_children_of(wire)? {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        let Some(repr) = edge_data.pcurve_for(data.surface, edge.location()) else {
            return Ok(None);
        };
        pieces += 1;
        match repr {
            EdgeRepr::PCurve { curve, range, .. } => {
                let Some(pcurve) = model.geometry().pcurve(*curve) else {
                    return Ok(None);
                };
                match pcurve {
                    ogeom_geom::PlanarCurve::Line(_) => {
                        use ogeom_geom::Curve2d as _;
                        let a = pcurve.point_at(range.0, tol)?;
                        let b = pcurve.point_at(range.1, tol)?;
                        segments.push((a, b));
                    }
                    ogeom_geom::PlanarCurve::Circle(arc) => {
                        use ogeom_geom::Curve2d as _;
                        // One circle's arcs, however many pieces the
                        // boundary arrives in. A boolean splits a closed rim
                        // to give the arrangement's walker somewhere to
                        // start, even a rim it never touched, and the disc
                        // those arcs bound is the same disc the whole turn
                        // bounded. The spans are summed and the total asked
                        // for a turn, so a fan of arcs that does not close
                        // is still no disc.
                        let span = (range.1 - range.0).abs();
                        arc_ends.push(pcurve.point_at(range.0, tol)?);
                        arc_ends.push(pcurve.point_at(range.1, tol)?);
                        match &mut circle {
                            None => circle = Some((*arc, span)),
                            Some((held, total)) => {
                                let (a, b) = (held.circle(), arc.circle());
                                if a.centre().distance(b.centre()) > tol.confusion()
                                    || (a.radius() - b.radius()).abs() > tol.confusion()
                                {
                                    return Ok(None);
                                }
                                *total += span;
                            }
                        }
                    }
                    _ => return Ok(None),
                }
            }
            EdgeRepr::Seam {
                forward,
                reversed,
                range,
                ..
            } => {
                use ogeom_geom::Curve2d as _;
                if seams_seen.contains(&edge.node()) {
                    pieces -= 1;
                    continue;
                }
                seams_seen.push(edge.node());
                // A seam says where the chart's edge stands, not how far
                // along it the face reaches. A reader pads its pcurve's
                // domain and a boolean shortens the face without shortening
                // the seam, so its own extent is worth nothing; the rims
                // are what say how tall the chart is, and they are ordinary
                // edges with ordinary ranges.
                for id in [forward, reversed] {
                    let Some(pcurve) = model.geometry().pcurve(*id) else {
                        return Ok(None);
                    };
                    let ogeom_geom::PlanarCurve::Line(_) = pcurve else {
                        return Ok(None);
                    };
                    let (lo, hi) = pcurve.domain();
                    let at = pcurve.point_at(range.0.clamp(lo, hi), tol)?;
                    let far = pcurve.point_at(range.1.clamp(lo, hi), tol)?;
                    if (at.x - far.x).abs() <= (at.y - far.y).abs() {
                        columns.push(at.x);
                    } else {
                        rows.push(at.y);
                    }
                }
            }
            _ => return Ok(None),
        }
    }

    if let Some((arc, span)) = circle {
        // The disc: one circle's arcs and nothing else, closing a turn, on
        // a plane.
        let _ = pieces;
        if !segments.is_empty()
            || !columns.is_empty()
            || !rows.is_empty()
            || (span - core::f64::consts::TAU).abs() > tol.parametric().max(1e-9)
        {
            return Ok(None);
        }
        // And the arcs must *tile* the turn rather than add up to one. A
        // reader that re-bases each edge's range onto its own curve can
        // leave two arcs both starting at the circle's zero, one a quarter
        // of it and one three quarters, which sums to a turn while
        // covering a quarter of the circle twice and half of it never. Each
        // arc end meets exactly one other where they genuinely chain.
        let reach = tol.confusion() * 10.0;
        for (index, at) in arc_ends.iter().enumerate() {
            let met = arc_ends
                .iter()
                .enumerate()
                .filter(|(other, q)| *other != index && q.distance(*at) <= reach)
                .count();
            if met != 1 {
                return Ok(None);
            }
        }
        let ogeom_geom::SurfaceGeometry::Plane(plane) = placed else {
            return Ok(None);
        };
        let frame = plane.plane().frame();
        let centre2 = arc.circle().centre();
        let centre = placed.point_at(centre2.x, centre2.y, tol)?;
        let radius = arc.circle().radius();
        let normal = frame.z().vector();
        return Ok(Some(ExactFace::Disc {
            centre,
            e1: frame.x().vector(),
            e2: frame.y().vector(),
            normal,
            radius,
            sign,
            share,
        }));
    }

    // A chart rectangle: every segment axis-aligned and on the hull's edge.
    // A torus's face is all seam and has no segments at all: two columns
    // and two rows, which are the rectangle.
    if segments.is_empty() && columns.is_empty() && rows.is_empty() {
        return Ok(None);
    }
    let (mut u0, mut u1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY);
    for (a, b) in &segments {
        for p in [a, b] {
            u0 = u0.min(p.x);
            u1 = u1.max(p.x);
            v0 = v0.min(p.y);
            v1 = v1.max(p.y);
        }
    }
    for u in &columns {
        u0 = u0.min(*u);
        u1 = u1.max(*u);
    }
    for v in &rows {
        v0 = v0.min(*v);
        v1 = v1.max(*v);
    }
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        return Ok(None);
    }
    if u1 - u0 <= tol.confusion() || v1 - v0 <= tol.confusion() {
        return Ok(None);
    }
    let eps = tol.confusion().max(1e-9 * (u1 - u0).max(v1 - v0));
    let on_side =
        |value: f64, lo: f64, hi: f64| (value - lo).abs() <= eps || (value - hi).abs() <= eps;
    let mut perimeter = 0.0;
    for (a, b) in &segments {
        let horizontal = (a.y - b.y).abs() <= eps;
        let vertical = (a.x - b.x).abs() <= eps;
        if !(horizontal ^ vertical) {
            return Ok(None);
        }
        if horizontal && !on_side(a.y, v0, v1) {
            return Ok(None);
        }
        if vertical && !on_side(a.x, u0, u1) {
            return Ok(None);
        }
        perimeter += a.distance(*b);
    }
    #[allow(clippy::cast_precision_loss)]
    for (values, lo, hi, span) in [(&columns, u0, u1, v1 - v0), (&rows, v0, v1, u1 - u0)] {
        for value in values {
            if !on_side(*value, lo, hi) {
                return Ok(None);
            }
            perimeter += span;
        }
    }
    let expected = 2.0 * ((u1 - u0) + (v1 - v0));
    if (perimeter - expected).abs() > 1e-6 * expected {
        return Ok(None);
    }
    Ok(Some(ExactFace::ChartRectangle {
        surface: placed.clone(),
        rect: (u0, u1, v0, v1),
        sign,
        share,
    }))
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

        // ∫ x_i x_j = m/(n(n+1)) · [ Σ p p_ᵀ + (Σ p)(Σ p)ᵀ ]; see the module
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
        // Then shift to the centre, the reverse of `inertia_about`.
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
    use ogeom_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn analytic_primitives_measure_exactly_on_their_own_surfaces() {
        // The exact path reports zero deflection and machine-precision
        // numbers: no chord band, no inscribed deficit.
        let mut model = Model::new();
        let pi = core::f64::consts::PI;

        let cylinder = crate::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        let props = volume_properties(&model, &cylinder.shape, fine(), T).unwrap();
        assert_eq!(props.deflection, 0.0, "the exact path was taken");
        assert_relative_eq!(props.mass, pi * 4.0 * 5.0, epsilon = 1e-10);
        assert!(props.centre.is_equal(Point::new(0.0, 0.0, 2.5), T));
        // I_zz of a solid cylinder: m r^2 / 2.
        let m = pi * 4.0 * 5.0;
        assert_relative_eq!(props.inertia.rows[2][2], m * 4.0 / 2.0, epsilon = 1e-8);

        let sphere = crate::make_sphere(&mut model, Frame::WORLD, 3.0, T).unwrap();
        let props = volume_properties(&model, &sphere.shape, fine(), T).unwrap();
        assert_eq!(props.deflection, 0.0);
        assert_relative_eq!(props.mass, 4.0 / 3.0 * pi * 27.0, epsilon = 1e-10);
        // I = 2/5 m r^2 about any axis through the centre.
        let m = 4.0 / 3.0 * pi * 27.0;
        assert_relative_eq!(props.inertia.rows[0][0], 0.4 * m * 9.0, epsilon = 1e-8);

        let torus = crate::make_torus(&mut model, Frame::WORLD, 5.0, 1.5, T).unwrap();
        let props = volume_properties(&model, &torus.shape, fine(), T).unwrap();
        assert_eq!(props.deflection, 0.0);
        assert_relative_eq!(props.mass, 2.0 * pi * pi * 5.0 * 1.5 * 1.5, epsilon = 1e-10);

        let cone = crate::make_cone(&mut model, Frame::WORLD, 3.0, 1.0, 4.0, T).unwrap();
        let props = volume_properties(&model, &cone.shape, fine(), T).unwrap();
        assert_eq!(props.deflection, 0.0);
        // A frustum: pi h (R^2 + R r + r^2) / 3.
        assert_relative_eq!(
            props.mass,
            pi * 4.0 * (9.0 + 3.0 + 1.0) / 3.0,
            epsilon = 1e-10
        );

        // Areas ride the same path: a sphere's is 4 pi r^2, exactly.
        let props = surface_properties(&model, &sphere.shape, fine(), T).unwrap();
        assert_eq!(props.deflection, 0.0);
        assert_relative_eq!(props.mass, 4.0 * pi * 9.0, epsilon = 1e-10);
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
        // digits at this distance and their difference keeps four, so the answer
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
        // truth and climbs as the deflection tightens, and the deflection is
        // reported, so a caller can see how far under.
        use crate::build::make_natural_face;
        use ogeom_geom::SphereSurface;
        use ogeom_math::Sphere;

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
