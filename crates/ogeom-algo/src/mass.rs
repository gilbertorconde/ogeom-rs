//! Mass properties: how much there is, where its centre is, and how it resists
//! being spun.
//!
//! Three measures, one for each dimension a shape can have (the length of its
//! edges, the area of its faces, the volume it encloses), each with the centre
//! of that measure and the inertia tensor about that centre.
//!
//! # Integrated on the surfaces where possible, meshed where not
//!
//! Area and volume are integrated on the exact surfaces first. A face on a
//! plane, cylinder, cone, sphere or torus bounded by a chart rectangle or a
//! full circle has a closed form; any other face with pcurves is integrated
//! round its chart boundary by Green's theorem. Either way the result
//! reports a deflection of zero.
//!
//! A shape with a face neither can take (no pcurves, a scaling placement,
//! a boundary that does not close in the chart) is measured on its
//! tessellation instead, and the result carries the deflection it was
//! computed at. Halving the deflection and seeing the answer move tells a
//! caller how much to trust it; [`MassProperties::deflection`] is what
//! makes that check possible. Lengths are always measured on a
//! discretization.
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
/// Each shell of a solid with several is weighed facing out of the solid's
/// material, as a probe off its faces finds: a void whose faces all point
/// into the material (turned inside out as a whole, which no edge between
/// its faces shows) is still taken away, not added.
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
    if let Some(probed) = probed_mesh_volume(model, shape, deflection, tol)? {
        return Ok(probed);
    }
    let (mut mesh, chords) = ogeom_mesh::triangulate_with_chords(model, shape, deflection, tol)?;
    if mesh.is_empty() {
        return Ok(MassProperties::none(deflection.chord));
    }
    if !mesh.is_closed() {
        // The flux through a boundary is the sum of the flux through each
        // face, and a face's share needs only that face's own mesh. The
        // welded mesh can stay open where faces meet on edges looser than
        // the weld dares reach, a mesh converted back holding threaded
        // curves a few hundredths off; closure is then asked of the
        // topology, which is what the divergence theorem is about.
        let shells = explore_unique(model, shape, ShapeType::Shell)?;
        let mut closed = !shells.is_empty();
        for shell in &shells {
            closed &= crate::build::is_shell_closed(model, shell)?;
        }
        if !closed {
            ogeom_bail!(
                Construction,
                "the boundary is not closed, so it encloses no volume to measure"
            );
        }
        mesh = ogeom_topo::Triangulation::new();
        let faces = explore(model, shape, Filter::OfType(ShapeType::Face))?;
        let meshes = ogeom_core::parallel::map_ordered(&faces, |_, face| {
            ogeom_mesh::triangulate_face_with(model, face, deflection, &chords, tol)
        });
        for face_mesh in meshes {
            mesh.append(&face_mesh?);
        }
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

/// The meshed volume of a shape holding a solid with more than one shell,
/// each face turned out by probing it against its own solid; `None` where
/// no solid of `shape` has more than one shell.
///
/// The whole-shape mesh mends a minority of turned faces within each
/// connected piece, and a void is a piece of its own: turned inside out as
/// a whole it agrees with itself and is weighed as material. So every face
/// under a solid is asked which way it faces, as `check` asks, and one
/// facing in is counted turned over. A face the probe cannot settle is
/// counted as its flag says.
fn probed_mesh_volume(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Option<MassProperties>> {
    let solids = explore(model, shape, Filter::OfType(ShapeType::Solid))?;
    let mut several = false;
    for solid in &solids {
        several |= explore(model, solid, Filter::OfType(ShapeType::Shell))?.len() > 1;
    }
    if !several {
        return Ok(None);
    }
    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        if !crate::build::is_shell_closed(model, &shell)? {
            ogeom_bail!(
                Construction,
                "the boundary is not closed, so it encloses no volume to measure"
            );
        }
    }
    // Every face with the solid it bounds, then any face under none.
    let mut faces: Vec<(Shape, Option<usize>)> = Vec::new();
    let mut seen: std::collections::HashSet<Occurrence> = std::collections::HashSet::new();
    for (at, solid) in solids.iter().enumerate() {
        for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
            seen.insert((face.node(), face.location().clone()));
            faces.push((face, Some(at)));
        }
    }
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        if !seen.contains(&(face.node(), face.location().clone())) {
            faces.push((face, None));
        }
    }
    let mut boundaries: Vec<Option<crate::SolidBoundary>> = Vec::with_capacity(solids.len());
    for solid in &solids {
        boundaries.push(or_mesh(
            crate::check::probe_boundary(model, solid, tol).map(Some),
            None,
        )?);
    }
    let chords = ogeom_mesh::edge_chords_for(model, shape, deflection, tol)?;
    let meshes = ogeom_core::parallel::map_ordered(&faces, |_, (face, at)| {
        let mesh = ogeom_mesh::triangulate_face_with(model, face, deflection, &chords, tol)?;
        let boundary = at.and_then(|at| boundaries[at].as_ref());
        let inward = match boundary {
            Some(boundary) => {
                or_mesh(crate::check::faces_inward(model, face, boundary, tol), None)?
            }
            None => None,
        };
        Ok((mesh, inward == Some(true)))
    });
    let mut acc = Accumulator::new();
    let mut apex: Option<Point> = None;
    for one in meshes {
        let (mesh, turned) = one?;
        let Some(&first) = mesh.positions.first() else {
            continue;
        };
        // Any apex serves, as for the whole-shape mesh.
        let apex = *apex.get_or_insert(first);
        let sign = if turned { -1.0 } else { 1.0 };
        for triangle in &mesh.triangles {
            let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
            let volume = (a - apex).dot((b - apex).cross(c - apex)) / 6.0;
            acc.add(&[apex, a, b, c], volume * sign);
        }
    }
    if apex.is_none() {
        return Ok(Some(MassProperties::none(deflection.chord)));
    }
    if acc.mass < 0.0 {
        ogeom_bail!(
            Construction,
            "the boundary is wound inward, so the volume came out negative"
        );
    }
    Ok(Some(acc.finish(deflection.chord)))
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
    /// Any other face, integrated round its chart loops.
    Chart(Box<crate::mass_chart::ChartFace>),
}

impl ExactFace {
    /// Whether this region is part of its face or taken out of it: `1` for
    /// the outer boundary, `-1` for a hole. The volume integral could carry
    /// it in the normal's sign, but the area integral takes a magnitude and
    /// would hand a hole's area back as more face.
    const fn share(&self) -> f64 {
        match self {
            Self::ChartRectangle { share, .. } | Self::Disc { share, .. } => *share,
            Self::Chart(_) => 1.0,
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
            Self::Chart(_) => 0.0,
        }
    }

    fn take_away(&mut self) {
        match self {
            Self::ChartRectangle { share, .. } | Self::Disc { share, .. } => *share = -1.0,
            Self::Chart(_) => {}
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
    // The face each region of `exact` belongs to.
    let mut region_of: Vec<usize> = Vec::with_capacity(faces.len());
    for (at, face) in faces.iter().enumerate() {
        // A face the closed forms cannot evaluate (a chart point a hair off
        // its surface's domain) is left to the mesh, like one they do not
        // speak at all.
        match or_mesh(integrable_face(model, face, tol), None)? {
            Some(found) => {
                region_of.extend(std::iter::repeat_n(at, found.len()));
                exact.extend(found);
            }
            None => {
                if *DEBUG_MASS {
                    eprintln!(
                        "MASS face {} is not exactly integrable",
                        face.node().index()
                    );
                }
                return Ok(None);
            }
        }
    }
    // And the faces must agree with each other about which way is out. A
    // boundary that cannot be walked to ask (a pcurve whose domain falls
    // short of its edge's range) is left to the mesh.
    let Some(flags) = or_mesh(flags_agree(model, shape, tol).map(Some), None)? else {
        return Ok(None);
    };
    if !flags.agree {
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
    // Sets of faces the walks could not tie together are asked of the
    // solid itself, a face of each probed off both its sides. A set facing
    // in that is a whole shell of a solid with several (a void turned
    // inside out) is turned over as a whole. Any other set turned in, or
    // one no face of which the probe can settle, is left to the mesh,
    // which mends a minority of turned faces.
    let mut turned: std::collections::HashSet<Occurrence> = std::collections::HashSet::new();
    if !flags.sets.is_empty() {
        let Some(inward) = sets_facing_in(model, shape, &flags.sets, tol)? else {
            return Ok(None);
        };
        let shells = shells_of_several(model, shape)?;
        for (set, inward) in flags.sets.iter().zip(inward) {
            if !inward {
                continue;
            }
            let held: std::collections::HashSet<Occurrence> = set
                .iter()
                .map(|f| (f.node(), f.location().clone()))
                .collect();
            if !shells.contains(&held) {
                return Ok(None);
            }
            turned.extend(held);
        }
    }
    let sign: Vec<f64> = region_of
        .iter()
        .map(|&at| {
            let face = &faces[at];
            if turned.contains(&(face.node(), face.location().clone())) {
                -1.0
            } else {
                1.0
            }
        })
        .collect();

    let reference = reference_point(&exact, tol)?;
    // Each face's moments are summed on their own, and added in the faces'
    // order.
    let summed = ogeom_core::parallel::map_ordered(&exact, |_, face| {
        integrate_face(
            face,
            crate::mass_chart::Measure::Volume,
            reference,
            tol,
            Moments::zero,
            |sums, p, n_da, share| {
                let n_da = n_da * share;
                let q = p - reference;
                sums.mass += q.dot(n_da) / 3.0;
                sums.first += Vector::new(
                    q.x * q.x * n_da.x / 2.0,
                    q.y * q.y * n_da.y / 2.0,
                    q.z * q.z * n_da.z / 2.0,
                );
                let d = [q.x, q.y, q.z];
                let nd = [n_da.x, n_da.y, n_da.z];
                for i in 0..3 {
                    // Diagonal: int q_i^2 dV = surface int q_i^3 n_i / 3.
                    sums.second.rows[i][i] += d[i] * d[i] * d[i] * nd[i] / 3.0;
                    // Off-diagonal: int q_i q_j dV = surface int q_i^2 q_j n_i / 2.
                    for j in 0..3 {
                        if i != j {
                            sums.second.rows[i][j] += d[i] * d[i] * d[j] * nd[i] / 2.0;
                        }
                    }
                }
            },
        )
    });
    let mut total = Moments::zero();
    for (face, sign) in summed.into_iter().zip(sign) {
        match or_mesh(face, None)? {
            Some(sums) => total.add_signed(&sums, sign),
            None => return Ok(None),
        }
    }
    let Moments {
        mass,
        first,
        mut second,
    } = total;
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
        match or_mesh(integrable_face(model, face, tol), None)? {
            Some(found) => exact.extend(found),
            None => return Ok(None),
        }
    }
    let reference = reference_point(&exact, tol)?;
    let mut total = Moments::zero();
    for face in &exact {
        let found = integrate_face(
            face,
            crate::mass_chart::Measure::Area,
            reference,
            tol,
            Moments::zero,
            |sums, p, n_da, share| {
                let da = n_da.magnitude() * share;
                let q = p - reference;
                sums.mass += da;
                sums.first += q * da;
                for (i, qi) in [q.x, q.y, q.z].iter().enumerate() {
                    for (j, qj) in [q.x, q.y, q.z].iter().enumerate() {
                        sums.second.rows[i][j] += qi * qj * da;
                    }
                }
            },
        );
        match or_mesh(found, None)? {
            Some(sums) => total.add(&sums),
            None => return Ok(None),
        }
    }
    let Moments {
        mass,
        first,
        second,
    } = total;
    let acc = Accumulator {
        reference: Some(reference),
        mass,
        first,
        second,
    };
    Ok(Some(acc.finish(0.0)))
}

/// A measure and its first and second moments about the reference, summed
/// over one face's samples or over the faces.
struct Moments {
    mass: f64,
    first: Vector,
    second: Matrix3,
}

impl Moments {
    const fn zero() -> Self {
        Self {
            mass: 0.0,
            first: Vector::ZERO,
            second: Matrix3::ZERO,
        }
    }

    fn add(&mut self, other: &Self) {
        self.add_signed(other, 1.0);
    }

    /// Add `other` times `sign`: `-1` for a face turned over.
    fn add_signed(&mut self, other: &Self, sign: f64) {
        self.mass += other.mass * sign;
        self.first += other.first * sign;
        self.second = add(self.second, scale_matrix(other.second, sign));
    }
}

/// An exact path's answer, or `fallback` where the closed forms could not
/// evaluate (a chart point a hair off its surface's domain): the mesh then
/// answers, for volumes and areas alike. A cancelled watch and a broken
/// model are errors whichever path meets them.
fn or_mesh<T>(result: OgeomResult<T>, fallback: T) -> OgeomResult<T> {
    match result {
        Ok(value) => Ok(value),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            Err(e)
        }
        Err(_) => Ok(fallback),
    }
}

/// Whether the exact path says which faces it leaves to the mesh. Read once.
static DEBUG_MASS: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var_os("OGEOM_DEBUG_MASS").is_some());

/// Somewhere on the shape to measure moments from.
fn reference_point(faces: &[ExactFace], tol: Tolerances) -> OgeomResult<Point> {
    use ogeom_geom::Surface as _;
    match &faces[0] {
        ExactFace::ChartRectangle { surface, rect, .. } => surface.point_at(rect.0, rect.2, tol),
        ExactFace::Disc { centre, .. } => Ok(*centre),
        ExactFace::Chart(chart) => chart.anchor(tol),
    }
}

/// Sum every quadrature sample of a face into an accumulator from `fresh`;
/// `None` where the face's rule did not settle on what `measure` sums.
///
/// `contribute` receives the world point, the outward-signed `n dA`
/// already weighted, and the region's share; summing those contributions
/// *is* the integral.
fn integrate_face<A>(
    face: &ExactFace,
    measure: crate::mass_chart::Measure,
    reference: Point,
    tol: Tolerances,
    fresh: impl Fn() -> A,
    contribute: impl Fn(&mut A, Point, Vector, f64),
) -> OgeomResult<Option<A>> {
    let share = face.share();
    use ogeom_geom::Surface as _;
    const QUARTER: f64 = core::f64::consts::FRAC_PI_2;
    let mut sums = fresh();
    match face {
        ExactFace::Chart(chart) => Ok(chart.integrate(measure, reference, tol, fresh, contribute)),
        ExactFace::ChartRectangle {
            surface,
            rect,
            sign,
            ..
        } => {
            let (u0, u1, v0, v1) = *rect;
            // Panels no wider than a quarter turn along a parameter that may
            // be an angle, and a spline's also cut at its knots. A length
            // parameter (a plane's, a drum's, cone's or extrusion's height)
            // carries a polynomial integrand the rule takes in one panel
            // however long it runs.
            let (angular_u, angular_v) = match surface {
                ogeom_geom::SurfaceGeometry::Plane(_) => (false, false),
                ogeom_geom::SurfaceGeometry::Cylinder(_)
                | ogeom_geom::SurfaceGeometry::Cone(_)
                | ogeom_geom::SurfaceGeometry::Extrusion(_) => (true, false),
                _ => (true, true),
            };
            let breaks =
                |lo: f64, hi: f64, angular: bool, knots: Option<&ogeom_math::KnotVector>| {
                    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                    let panels = if angular {
                        ((hi - lo) / QUARTER).ceil().clamp(1.0, 64.0) as usize
                    } else {
                        1
                    };
                    #[allow(clippy::cast_precision_loss)]
                    let mut out: Vec<f64> = (0..=panels)
                        .map(|i| lo + (hi - lo) * i as f64 / panels as f64)
                        .collect();
                    if let Some(knots) = knots {
                        out.extend(
                            knots
                                .distinct()
                                .into_iter()
                                .map(|(k, _)| k)
                                .filter(|k| *k > lo && *k < hi),
                        );
                        out.sort_by(f64::total_cmp);
                        out.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
                    }
                    out
                };
            // A swept curve's knots stand across its sweep: in `u` for an
            // extrusion, in `v` for a revolution.
            fn curve_knots(curve: &ogeom_geom::Curve) -> Option<&ogeom_math::KnotVector> {
                match curve {
                    ogeom_geom::Curve::BSpline(b) => Some(b.knots()),
                    ogeom_geom::Curve::Trimmed(t) => curve_knots(t.basis()),
                    _ => None,
                }
            }
            let (u_knots, v_knots) = match surface {
                ogeom_geom::SurfaceGeometry::BSpline(b) => (Some(b.u_knots()), Some(b.v_knots())),
                ogeom_geom::SurfaceGeometry::Extrusion(e) => (curve_knots(e.curve()), None),
                ogeom_geom::SurfaceGeometry::Revolution(r) => (None, curve_knots(r.curve())),
                _ => (None, None),
            };
            let (u_breaks, v_breaks) = (
                breaks(u0, u1, angular_u, u_knots),
                breaks(v0, v1, angular_v, v_knots),
            );
            let mut failure = None;
            for uw in u_breaks.windows(2) {
                let (ua, ub) = (uw[0], uw[1]);
                for vw in v_breaks.windows(2) {
                    let (va, vb) = (vw[0], vw[1]);
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
                            contribute(&mut sums, p, du.cross(dv) * (sign * weight), share);
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
                None => Ok(Some(sums)),
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
                    contribute(&mut sums, p, *normal * (sign * rho * weight), share);
                });
            }
            match failure {
                Some(e) => Err(e),
                None => Ok(Some(sums)),
            }
        }
    }
}

/// A tensor-product ten-by-ten Gauss rule over `[a,b] x [c,d]`, feeding each
/// sample and its weight to the callback.
fn gauss2(a: f64, b: f64, c: f64, d: f64, f: &mut dyn FnMut(f64, f64, f64)) {
    let us = ogeom_math::gauss_legendre_rule(a, b);
    let vs = ogeom_math::gauss_legendre_rule(c, d);
    for &(u, wu) in &us {
        for &(v, wv) in &vs {
            f(u, v, wu * wv);
        }
    }
}

/// A face's regions in closed form where its surface and trim allow, and
/// otherwise its chart loops for integrating round; `None` where neither
/// can be had.
fn integrable_face(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<ExactFace>>> {
    if let Some(found) = exact_face(model, face, tol)? {
        return Ok(Some(found));
    }
    Ok(crate::mass_chart::chart_face(model, face, tol)
        .map(|chart| vec![ExactFace::Chart(Box::new(chart))]))
}

/// The exact-integrable regions of one face, or `None` where there are
/// none.
///
/// One region per wire, and the integral is their sum: a face's outer
/// boundary carries its own sign and every inner one the opposite, which
/// is what a hole *is* under the divergence theorem. So a plate with a
/// bore in it is a rectangle less a disc, and a tube's end face a disc
/// less a disc, neither of which needs meshing.
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
            // A spline trimmed by its chart's own borders: a rectangle,
            // integrated knot span by knot span, where each span is one
            // polynomial piece the Gauss rule takes exactly.
            | ogeom_geom::SurfaceGeometry::BSpline(_)
            | ogeom_geom::SurfaceGeometry::Extrusion(_)
            | ogeom_geom::SurfaceGeometry::Revolution(_)
    );
    if !analytic {
        return Ok(None);
    }
    let placement = face.transform(model.datums())?;
    // The chart rectangle comes from the pcurves, whose windows are the
    // *unscaled* surface's; a scaling placement changes the chart's metric
    // and the windows with it, so only rigid placements take the exact path.
    // A reflecting one turns the placed chart's normal against the face's.
    let Some(handedness) = crate::mass_chart::rigid_handedness(&placement) else {
        return Ok(None);
    };
    let placed = surface.clone().transformed(&placement, tol)?;
    let sign = handedness
        * if face.orientation() == ogeom_topo::Orientation::Reversed {
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

/// An edge or face occurrence: its node at its placement.
type Occurrence = (ogeom_topo::TShapeId, ogeom_topo::Location);

/// Edge curves placed once for [`flags_agree`], with their ranges, by edge
/// occurrence: each is read at several stations of each face it bounds.
type PlacedCurves = std::collections::HashMap<Occurrence, (ogeom_geom::Curve, (f64, f64))>;

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
/// A bore wall whose flag points into the solid walks its edges the same
/// way as its neighbours. The tessellator repairs such a shell, flipping
/// whichever side of the disagreement is in the minority, and the
/// closed-form integral cannot: it would hand the bore back as material.
/// So where the flags disagree this says so and the mesh is asked instead.
///
/// The comparison ties faces into sets that agree among themselves, and
/// says nothing about one set against another: a face the walks cannot be
/// read on (a station at a cone's apex, where the surface has no normal, a
/// trim that is no closed loop) is a set of its own, and so is each
/// closed shell. Turning a whole set over turns its share of the volume
/// over without any edge noticing, so where there is more than one set
/// [`Flags::sets`] lists the faces of each, for the caller to ask which
/// way the set faces. A single set turned over as a whole is the whole
/// boundary wound inward, and the volume's sign says so.
///
/// An edge occurrence is its node at its placement: the top and bottom of
/// a prism are one edge moved, and each is compared only with the faces
/// that meet it where it stands.
pub(crate) fn flags_agree(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Flags> {
    // Moving the whole shape moves every face with it and changes nothing
    // about whether they agree.
    let unplaced = shape.located(ogeom_topo::Location::default());
    let faces = explore(model, &unplaced, Filter::OfType(ShapeType::Face))?;
    let mut walks: std::collections::HashMap<Occurrence, Vec<(bool, usize, bool)>> =
        std::collections::HashMap::new();
    let mut curves: PlacedCurves = std::collections::HashMap::new();
    for (index, face) in faces.iter().enumerate() {
        match face_walks(model, face, &mut curves, tol)? {
            Some(found) => {
                for (edge, ahead, outer) in found {
                    walks.entry(edge).or_default().push((ahead, index, outer));
                }
            }
            None => {
                if *DEBUG_MASS {
                    eprintln!("MASS face {:?} cannot be walked", face.node());
                }
            }
        }
    }
    if *DEBUG_MASS {
        eprintln!("MASS flags_agree walked {} edges", walks.len());
    }
    // Which set each face is in, as a forest of parents.
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    for ((edge, _), uses) in &walks {
        // An edge one face walks twice is that face's own seam, however it
        // is written down (a canonicalised drum keeps its as an ordinary
        // pcurve used twice), and one face's seam says nothing about
        // whether two faces agree.
        if uses.iter().all(|(_, owner, _)| *owner == uses[0].1) {
            continue;
        }
        let ahead = uses.iter().filter(|(ahead, ..)| *ahead).count();
        if ahead * 2 != uses.len() {
            if *DEBUG_MASS {
                eprintln!("MASS edge {} is walked {uses:?}", edge.index());
            }
            return Ok(Flags {
                agree: false,
                sets: Vec::new(),
            });
        }
        let first = set_of(&mut parent, uses[0].1);
        for (_, owner, _) in uses {
            let other = set_of(&mut parent, *owner);
            parent[other] = first;
        }
    }
    // Each set's faces in the order they were explored, named as the
    // caller's shape holds them, placement and all.
    let mut sets: Vec<Vec<Shape>> = Vec::new();
    let mut set_at: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for (index, face) in faces.iter().enumerate() {
        let root = set_of(&mut parent, index);
        let at = *set_at.entry(root).or_insert_with(|| {
            sets.push(Vec::new());
            sets.len() - 1
        });
        sets[at].push(face.moved(shape.location()));
    }
    if sets.len() < 2 {
        sets.clear();
    }
    Ok(Flags { agree: true, sets })
}

/// The face standing for the set `at` is in, shortening the way there.
fn set_of(parent: &mut [usize], mut at: usize) -> usize {
    while parent[at] != at {
        parent[at] = parent[parent[at]];
        at = parent[at];
    }
    at
}

/// What [`flags_agree`] found.
pub(crate) struct Flags {
    /// Whether every edge the walks could read is walked opposite ways by
    /// the faces either side of it.
    pub agree: bool,
    /// Where the faces fall into more than one set that agrees within
    /// itself, the faces of each set; empty where they are one set.
    pub sets: Vec<Vec<Shape>>,
}

/// One station of a face's walk round its boundary: the edge occurrence,
/// whether the walk runs with the edge's own direction there, and whether
/// the station is on the face's outer boundary.
type EdgeWalk = (Occurrence, bool, bool);

/// Each station of `face`'s walk along its edges; `None` where the walk
/// cannot be read on it.
fn face_walks(
    model: &Model,
    face: &Shape,
    curves: &mut PlacedCurves,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<EdgeWalk>>> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let placed_at = ogeom_topo::Location::default();
    let mut walks = Vec::new();
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()).cloned() else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placed = surface
        .clone()
        .transformed(&face.transform(model.datums())?, tol)?;
    let flag = if face.orientation() == ogeom_topo::Orientation::Reversed {
        -1.0
    } else {
        1.0
    };
    // Where the boundary walks into closed chart loops, their windings
    // say which side the face lies on at every point, concave or not.
    if let Some(stations) = crate::mass_chart::material_sides(model, face, tol) {
        for (edge, at, toward) in stations {
            let at = crate::mass_chart::into_domain(&placed, at);
            let (du, dv) = placed.d1_at(at.x, at.y, tol)?;
            let raw = du.cross(dv);
            let inward = du * toward.x + dv * toward.y;
            if raw.magnitude() <= tol.angular() || inward.magnitude() <= tol.angular() {
                return Ok(None);
            }
            let out = raw / raw.magnitude() * flag;
            let walk = out.cross(inward / inward.magnitude());
            let station = placed.point_at(at.x, at.y, tol)?;
            let Some(along) = edge_heading(model, &edge, station, curves, tol)? else {
                return Ok(None);
            };
            walks.push((
                (edge.node(), edge.location().clone()),
                walk.dot(along) > 0.0,
                true,
            ));
        }
        return Ok(Some(walks));
    }
    // The pcurves read below are the unplaced face's.
    if face.location() != &placed_at {
        return Ok(None);
    }
    // Otherwise each wire's middle in the chart stands in for the side
    // the face lies on, and which wire is the boundary:
    // the one covering the most of it, since a hole is inside what it
    // is a hole in.
    let wires = model.ordered_children_of(face)?;
    let mut middles: Vec<(ogeom_math::Point2, f64)> = Vec::with_capacity(wires.len());
    let mut stations: Vec<Vec<(Shape, ogeom_math::Point2)>> = Vec::with_capacity(wires.len());
    // How often each edge bounds this face: a seam the face uses once
    // (a half band, cut along its seam) bounds it down one column only.
    let mut uses: std::collections::HashMap<ogeom_topo::TShapeId, usize> =
        std::collections::HashMap::new();
    for wire in &wires {
        for edge in model.ordered_children_of(wire)? {
            *uses.entry(edge.node()).or_default() += 1;
        }
    }
    for wire in &wires {
        let mut here = Vec::new();
        let mut sum = ogeom_math::Vector2::new(0.0, 0.0);
        let (mut lo, mut hi) = (
            ogeom_math::Point2::new(f64::INFINITY, f64::INFINITY),
            ogeom_math::Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        // Seams used once, their column chosen once the rest of the
        // wire says where the face lies: by the ends of the other
        // pieces, which meet the used column and not the other.
        let mut once: Vec<(Shape, [ogeom_topo::PCurveId; 2], (f64, f64))> = Vec::new();
        let mut ends: Vec<ogeom_math::Point2> = Vec::new();
        for edge in model.ordered_children_of(wire)? {
            if edge.location() != &placed_at {
                return Ok(None);
            }
            let Some(repr) = model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(data.surface, edge.location()))
            else {
                return Ok(None);
            };
            // A seam bounds its face twice, once down each column.
            let sides: Vec<(ogeom_topo::PCurveId, (f64, f64))> = match repr {
                EdgeRepr::PCurve { curve, range, .. } => vec![(*curve, *range)],
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                } if uses.get(&edge.node()) == Some(&1) => {
                    once.push((edge.clone(), [*forward, *reversed], *range));
                    continue;
                }
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                } => vec![(*forward, *range), (*reversed, *range)],
                _ => return Ok(None),
            };
            for (id, range) in sides {
                let Some(pcurve) = model.geometry().pcurve(id) else {
                    return Ok(None);
                };
                ends.push(pcurve.point_at(range.0, tol)?);
                ends.push(pcurve.point_at(range.1, tol)?);
                // Several stations along each edge, not one: a wire of
                // a single closed edge has its own midpoint for a
                // middle, and nothing lies from a point toward itself.
                const STATIONS: usize = 4;
                for step in 1..=STATIONS {
                    #[allow(clippy::cast_precision_loss)]
                    let t = range.0 + (range.1 - range.0) * (step as f64 / (STATIONS + 1) as f64);
                    let at = pcurve.point_at(t, tol)?;
                    sum += at.to_vector();
                    lo = ogeom_math::Point2::new(lo.x.min(at.x), lo.y.min(at.y));
                    hi = ogeom_math::Point2::new(hi.x.max(at.x), hi.y.max(at.y));
                    here.push((edge.clone(), at));
                }
            }
        }
        if here.is_empty() && once.is_empty() {
            return Ok(None);
        }
        // A seam used once runs down the column its neighbours meet.
        for (edge, sides, range) in once {
            let mut best: Option<(f64, ogeom_topo::PCurveId)> = None;
            for id in sides {
                let Some(pcurve) = model.geometry().pcurve(id) else {
                    return Ok(None);
                };
                let mut d = f64::INFINITY;
                for t in [range.0, range.1] {
                    let at = pcurve.point_at(t, tol)?;
                    for end in &ends {
                        d = d.min(at.distance(*end));
                    }
                }
                if best.is_none_or(|(held, _)| d < held) {
                    best = Some((d, id));
                }
            }
            let Some((_, id)) = best else {
                return Ok(None);
            };
            let Some(pcurve) = model.geometry().pcurve(id) else {
                return Ok(None);
            };
            const STATIONS: usize = 4;
            for step in 1..=STATIONS {
                #[allow(clippy::cast_precision_loss)]
                let t = range.0 + (range.1 - range.0) * (step as f64 / (STATIONS + 1) as f64);
                let at = pcurve.point_at(t, tol)?;
                sum += at.to_vector();
                lo = ogeom_math::Point2::new(lo.x.min(at.x), lo.y.min(at.y));
                hi = ogeom_math::Point2::new(hi.x.max(at.x), hi.y.max(at.y));
                here.push((edge.clone(), at));
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let middle = ogeom_math::Point2::ORIGIN + sum / here.len() as f64;
        middles.push((middle, (hi.x - lo.x) * (hi.y - lo.y)));
        stations.push(here);
    }
    let Some(outer) = (0..middles.len()).max_by(|a, b| middles[*a].1.total_cmp(&middles[*b].1))
    else {
        return Ok(None);
    };
    for (index, here) in stations.into_iter().enumerate() {
        let (middle, _) = middles[index];
        for (edge, at) in here {
            let at = crate::mass_chart::into_domain(&placed, at);
            let (du, dv) = placed.d1_at(at.x, at.y, tol)?;
            let raw = du.cross(dv);
            if raw.magnitude() <= tol.angular() {
                return Ok(None);
            }
            let out = raw / raw.magnitude() * flag;
            // Which way the material lies from this point of the
            // boundary: toward the wire's middle for the face's outer
            // wire, away from it for a hole.
            let toward = middle - at;
            let toward = if index == outer { toward } else { -toward };
            let inward = du * toward.x + dv * toward.y;
            if inward.magnitude() <= tol.angular() {
                return Ok(None);
            }
            let walk = out.cross(inward / inward.magnitude());
            // Against the edge's own direction, so the two faces'
            // answers can be compared without comparing vectors. The
            // direction is read where the edge's curve passes the
            // station, since neither a pcurve's parameter nor its sense
            // need be its curve's.
            let station = placed.point_at(at.x, at.y, tol)?;
            let Some(along) = edge_heading(model, &edge, station, curves, tol)? else {
                return Ok(None);
            };
            walks.push((
                (edge.node(), edge.location().clone()),
                walk.dot(along) > 0.0,
                index == outer,
            ));
        }
    }
    Ok(Some(walks))
}

/// Whether each of `sets` faces into the solid it bounds; `None` where no
/// face of some set settles it.
///
/// A set agrees within itself, so any one of its faces answers for all of
/// them, and a face the probe cannot settle hands the question to the
/// next: the top of a sheet thinner than the probe's shortest step has
/// nothing but outside either side of it, while the sheet's rim answers.
/// One face of each set is asked first and more of each still unsettled
/// at every round after, in parallel. Each face is probed against its own
/// solid, not the whole shape: a compound's lumps may overlap, and a face
/// of one inside another has material on both its sides.
fn sets_facing_in(
    model: &Model,
    shape: &Shape,
    sets: &[Vec<Shape>],
    tol: Tolerances,
) -> OgeomResult<Option<Vec<bool>>> {
    let solids = explore(model, shape, Filter::OfType(ShapeType::Solid))?;
    // Which solid each face bounds, by occurrence.
    let mut owner: std::collections::HashMap<Occurrence, usize> = std::collections::HashMap::new();
    let solids = if solids.is_empty() {
        vec![shape.clone()]
    } else {
        solids
    };
    for (at, solid) in solids.iter().enumerate() {
        for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
            owner
                .entry((face.node(), face.location().clone()))
                .or_insert(at);
        }
    }
    let mut boundaries: Vec<Option<crate::SolidBoundary>> = solids.iter().map(|_| None).collect();
    // Each unsettled set with its faces still to be asked and their solids.
    let mut pending: Vec<(usize, std::collections::VecDeque<(usize, &Shape)>)> =
        Vec::with_capacity(sets.len());
    for (set, faces_of) in sets.iter().enumerate() {
        let mut faces = std::collections::VecDeque::with_capacity(faces_of.len());
        for face in faces_of {
            let Some(&at) = owner.get(&(face.node(), face.location().clone())) else {
                return Ok(None);
            };
            faces.push_back((at, face));
        }
        pending.push((set, faces));
    }
    let mut inward = vec![false; sets.len()];
    let mut take = 1;
    loop {
        let mut asked: Vec<(usize, usize, &Shape)> = Vec::new();
        for (set, faces) in &mut pending {
            for _ in 0..take {
                let Some((at, face)) = faces.pop_front() else {
                    break;
                };
                asked.push((*set, at, face));
            }
        }
        for &(_, at, _) in &asked {
            if boundaries[at].is_none() {
                let Some(boundary) = or_mesh(
                    crate::check::probe_boundary(model, &solids[at], tol).map(Some),
                    None,
                )?
                else {
                    return Ok(None);
                };
                boundaries[at] = Some(boundary);
            }
        }
        let facing = ogeom_core::parallel::map_ordered(&asked, |_, &(_, at, face)| {
            let Some(boundary) = boundaries[at].as_ref() else {
                return Ok(None);
            };
            crate::check::faces_inward(model, face, boundary, tol)
        });
        let mut settled = vec![false; sets.len()];
        for (&(set, ..), facing) in asked.iter().zip(facing) {
            let facing = or_mesh(facing, None)?;
            if *DEBUG_MASS {
                eprintln!("MASS set {set} probed {facing:?}");
            }
            if let Some(facing) = facing {
                inward[set] = facing;
                settled[set] = true;
            }
        }
        // A set out of faces with none settled is left to the mesh.
        for (set, faces) in &pending {
            if !settled[*set] && faces.is_empty() {
                if *DEBUG_MASS {
                    eprintln!("MASS set {set}: no face settles it");
                }
                return Ok(None);
            }
        }
        pending.retain(|(set, _)| !settled[*set]);
        if pending.is_empty() {
            return Ok(Some(inward));
        }
        take *= 2;
    }
}

/// The faces of every shell of a solid with more than one shell, each
/// shell's faces by occurrence.
fn shells_of_several(
    model: &Model,
    shape: &Shape,
) -> OgeomResult<Vec<std::collections::HashSet<Occurrence>>> {
    let mut out = Vec::new();
    for solid in explore(model, shape, Filter::OfType(ShapeType::Solid))? {
        let shells = explore(model, &solid, Filter::OfType(ShapeType::Shell))?;
        if shells.len() < 2 {
            continue;
        }
        for shell in shells {
            out.push(
                explore(model, &shell, Filter::OfType(ShapeType::Face))?
                    .iter()
                    .map(|f| (f.node(), f.location().clone()))
                    .collect(),
            );
        }
    }
    Ok(out)
}

/// An edge's own direction in space where its curve passes nearest `at`:
/// the best of a sampling over the edge's range, narrowed by golden
/// sections, and the curve's tangent there.
fn edge_heading(
    model: &Model,
    edge: &Shape,
    at: Point,
    curves: &mut PlacedCurves,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    use ogeom_geom::Curve3d as _;
    let key = (edge.node(), edge.location().clone());
    if !curves.contains_key(&key) {
        let Some((curve, range)) = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| match d.curve3d()? {
                EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
                _ => None,
            })
            .and_then(|(id, range)| Some((model.geometry().curve(id)?.clone(), range)))
        else {
            return Ok(None);
        };
        let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
        curves.insert(key.clone(), (placed, range));
    }
    let Some((curve, range)) = curves.get(&key) else {
        return Ok(None);
    };
    let range = *range;
    let gap = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(at)) };
    const SAMPLES: u32 = 32;
    let step = (range.1 - range.0) / f64::from(SAMPLES);
    let mut best = (range.0, gap(range.0)?);
    for k in 1..=SAMPLES {
        let t = range.0 + step * f64::from(k);
        let d = gap(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (mut a, mut b) = (
        (best.0 - step.abs()).max(range.0.min(range.1)),
        (best.0 + step.abs()).min(range.0.max(range.1)),
    );
    // Golden section, each round keeping one of its two interior points.
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    let (mut c, mut d) = (b - (b - a) * ratio, a + (b - a) * ratio);
    let (mut fc, mut fd) = (gap(c)?, gap(d)?);
    for _ in 0..60 {
        if fc < fd {
            b = d;
            (d, fd) = (c, fc);
            c = b - (b - a) * ratio;
            fc = gap(c)?;
        } else {
            a = c;
            (c, fc) = (d, fd);
            d = a + (b - a) * ratio;
            fd = gap(d)?;
        }
    }
    let along = curve.d1_at(f64::midpoint(a, b), tol)?;
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
                    // A trim fitted along a chart column or row (a rail a
                    // strip adopted from its neighbour): every control
                    // point on the one line makes it that segment.
                    ogeom_geom::PlanarCurve::BSpline(spline)
                        if along_one_chart_line(spline.control_points()) =>
                    {
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
                    // A chart column or row, stated as a line or as a
                    // spline whose control points all stand on one (a
                    // converted face's seams).
                    let straight = match pcurve {
                        ogeom_geom::PlanarCurve::Line(_) => true,
                        ogeom_geom::PlanarCurve::BSpline(spline) => {
                            along_one_chart_line(spline.control_points())
                        }
                        _ => false,
                    };
                    if !straight {
                        return Ok(None);
                    }
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
        // The chart's own normal, as every other region takes it: under a
        // reflection the frame's `z` is the image of the unplaced normal,
        // and the sign already carries the turn.
        let normal = frame.x().vector().cross(frame.y().vector());
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
    // A seam cut into pieces stands each piece on the same chart side; the
    // side is one side of the rectangle however many pieces it took.
    for values in [&mut columns, &mut rows] {
        values.sort_by(f64::total_cmp);
        values.dedup_by(|a, b| (*a - *b).abs() <= eps);
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
        // A simplex has at most four corners: no heap for them, once per
        // simplex of a whole mesh.
        let local: smallvec::SmallVec<[Vector; 4]> =
            points.iter().map(|p| *p - reference).collect();
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

/// Whether a chart curve's control points all stand on one column or one
/// row of the chart, to rounding against its extent.
fn along_one_chart_line(control: &[ogeom_math::Weighted<ogeom_math::Point2>]) -> bool {
    let points: Vec<ogeom_math::Point2> = control.iter().map(|w| w.point()).collect();
    let Some(first) = points.first() else {
        return false;
    };
    let extent = points
        .iter()
        .map(|p| p.distance(*first))
        .fold(0.0_f64, f64::max);
    let eps = 1e-9 * extent.max(1.0);
    points.iter().all(|p| (p.x - first.x).abs() <= eps)
        || points.iter().all(|p| (p.y - first.y).abs() <= eps)
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
