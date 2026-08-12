//! What hangs off a topology node: geometry, tolerances, edge representations.
//!
//! `docs/DATA_MODEL.md` §5 and §6.
//!
//! # Tolerances are per entity
//!
//! Every vertex, edge and face carries its own [`Tolerance`] — the radius
//! within which it is considered to lie. Operations may only widen one, and the
//! containment rule `tol(vertex) >= tol(edge) >= tol(face)` holds between
//! entities in a boundary relationship ([`check_containment`]).
//!
//! This is not a workaround for imprecise code. Exact arithmetic cannot
//! represent the intersection curve of two curved surfaces, so a
//! tolerance-carrying topology is the only known way to build a kernel whose
//! results close up. Every production kernel works this way.
//!
//! # An edge carries a list of representations
//!
//! Not one curve — a list. A single edge holds a 3D curve, *one pcurve per
//! adjacent face*, two pcurves where it is a seam on a closed surface, and
//! cached polylines. Face splitting during a boolean happens in a surface's
//! 2D parameter space, so without a pcurve on each face there is nothing to
//! split with; and since surfaces are parameterized differently, one 2D curve
//! cannot serve two faces.

use ogeom_core::{Arena, Key, OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry};
use ogeom_math::Point;
use smallvec::SmallVec;

use crate::location::Location;
use crate::tessellation::Triangulation;

/// A handle to a space curve.
pub type CurveId = Key<Curve>;
/// A handle to a curve in a surface's parameter space.
pub type PCurveId = Key<PlanarCurve>;
/// A handle to a surface.
pub type SurfaceId = Key<SurfaceGeometry>;

/// A handle to a cached triangulation.
pub type TriangulationId = Key<Triangulation>;

/// The geometry a model's topology refers into.
#[derive(Debug, Clone, Default)]
pub struct GeometryStore {
    curves: Arena<Curve>,
    pcurves: Arena<PlanarCurve>,
    surfaces: Arena<SurfaceGeometry>,
    triangulations: Arena<Triangulation>,
}

impl GeometryStore {
    /// An empty store.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            curves: Arena::new(),
            pcurves: Arena::new(),
            surfaces: Arena::new(),
            triangulations: Arena::new(),
        }
    }

    /// Add a space curve.
    pub fn add_curve(&mut self, curve: Curve) -> CurveId {
        self.curves.insert(curve)
    }

    /// Add a curve in parameter space.
    pub fn add_pcurve(&mut self, curve: PlanarCurve) -> PCurveId {
        self.pcurves.insert(curve)
    }

    /// Add a surface.
    pub fn add_surface(&mut self, surface: SurfaceGeometry) -> SurfaceId {
        self.surfaces.insert(surface)
    }

    /// Add a cached triangulation.
    pub fn add_triangulation(&mut self, mesh: Triangulation) -> TriangulationId {
        self.triangulations.insert(mesh)
    }

    /// The space curve behind `id`.
    #[must_use]
    pub fn curve(&self, id: CurveId) -> Option<&Curve> {
        self.curves.get(id)
    }

    /// The parameter-space curve behind `id`.
    #[must_use]
    pub fn pcurve(&self, id: PCurveId) -> Option<&PlanarCurve> {
        self.pcurves.get(id)
    }

    /// The surface behind `id`.
    #[must_use]
    pub fn surface(&self, id: SurfaceId) -> Option<&SurfaceGeometry> {
        self.surfaces.get(id)
    }

    /// The cached triangulation behind `id`.
    #[must_use]
    pub fn triangulation(&self, id: TriangulationId) -> Option<&Triangulation> {
        self.triangulations.get(id)
    }

    /// How many curves, pcurves and surfaces are held.
    #[must_use]
    pub fn counts(&self) -> (usize, usize, usize) {
        (self.curves.len(), self.pcurves.len(), self.surfaces.len())
    }

    /// The identifiers this store's arenas issue keys under.
    ///
    /// For [`Model::from_parts`](crate::Model::from_parts), which has to bind
    /// handles rebuilt by a reader to the arenas they will actually live in.
    pub(crate) const fn scopes(&self) -> GeometryScopes {
        GeometryScopes {
            curves: self.curves.scope(),
            pcurves: self.pcurves.scope(),
            surfaces: self.surfaces.scope(),
            triangulations: self.triangulations.scope(),
        }
    }

    /// Whether every arena has only ever been appended to.
    ///
    /// The precondition for extending the store by offset — see
    /// [`Arena::is_dense`](ogeom_core::Arena::is_dense).
    pub(crate) fn is_dense(&self) -> bool {
        self.curves.is_dense()
            && self.pcurves.is_dense()
            && self.surfaces.is_dense()
            && self.triangulations.is_dense()
    }

    /// Append another store's contents, returning where each kind landed.
    ///
    /// The returned offsets are this store's lengths before the append: an
    /// entry that sat at index `i` in `other` now sits at `i + offset`. The
    /// receiving arenas hand out the keys, so the values travel bare.
    pub(crate) fn append(&mut self, other: Self) -> GeometryOffsets {
        let offsets = GeometryOffsets {
            curves: arena_len(&self.curves),
            pcurves: arena_len(&self.pcurves),
            surfaces: arena_len(&self.surfaces),
            triangulations: arena_len(&self.triangulations),
        };
        for curve in other.curves.into_values() {
            self.curves.insert(curve);
        }
        for pcurve in other.pcurves.into_values() {
            self.pcurves.insert(pcurve);
        }
        for surface in other.surfaces.into_values() {
            self.surfaces.insert(surface);
        }
        for mesh in other.triangulations.into_values() {
            self.triangulations.insert(mesh);
        }
        offsets
    }

    /// Whether every piece of geometry a representation names is held here.
    ///
    /// The check a restored model needs: a handle that does not resolve is not
    /// a finding, it is a document that does not describe itself.
    #[must_use]
    pub fn holds(&self, repr: &EdgeRepr) -> bool {
        match repr {
            EdgeRepr::Curve3d { curve, .. } => self.curve(*curve).is_some(),
            EdgeRepr::PCurve { curve, surface, .. } => {
                self.pcurve(*curve).is_some() && self.surface(*surface).is_some()
            }
            EdgeRepr::Seam {
                forward,
                reversed,
                surface,
                ..
            } => {
                self.pcurve(*forward).is_some()
                    && self.pcurve(*reversed).is_some()
                    && self.surface(*surface).is_some()
            }
            // Carries its points itself, so there is nothing to resolve.
            EdgeRepr::Polyline { .. } => true,
            EdgeRepr::PolygonOnTriangulation { triangulation, .. } => {
                self.triangulation(*triangulation).is_some()
            }
        }
    }

    /// Every space curve, with its handle, in arena order.
    pub fn curves(&self) -> impl Iterator<Item = (CurveId, &Curve)> {
        self.curves.iter()
    }

    /// Every parameter-space curve, with its handle, in arena order.
    pub fn pcurves(&self) -> impl Iterator<Item = (PCurveId, &PlanarCurve)> {
        self.pcurves.iter()
    }

    /// Every surface, with its handle, in arena order.
    pub fn surfaces(&self) -> impl Iterator<Item = (SurfaceId, &SurfaceGeometry)> {
        self.surfaces.iter()
    }

    /// Every cached triangulation, with its handle, in arena order.
    pub fn triangulations(&self) -> impl Iterator<Item = (TriangulationId, &Triangulation)> {
        self.triangulations.iter()
    }

    /// How many cached triangulations are held.
    #[must_use]
    pub fn triangulation_count(&self) -> usize {
        self.triangulations.len()
    }
}

/// Which arena issues each kind of geometry handle in one store.
#[derive(Debug, Clone, Copy)]
pub(crate) struct GeometryScopes {
    pub curves: u32,
    pub pcurves: u32,
    pub surfaces: u32,
    pub triangulations: u32,
}

/// Where each kind of geometry landed in an append — the lengths of the
/// receiving arenas before it.
#[derive(Debug, Clone, Copy)]
pub(crate) struct GeometryOffsets {
    pub curves: u32,
    pub pcurves: u32,
    pub surfaces: u32,
    pub triangulations: u32,
}

/// An arena's length as the index its next append will land at.
///
/// # Panics
///
/// If the arena exceeds `u32::MAX` slots, which [`Arena::insert`] already
/// refuses to reach.
#[allow(clippy::expect_used, reason = "documented panic; see # Panics")]
pub(crate) fn arena_len<T>(arena: &ogeom_core::Arena<T>) -> u32 {
    u32::try_from(arena.len()).expect("arena exceeded u32::MAX slots")
}

/// Whether a key is unscoped and at generation zero — the state a reader
/// leaves handles in, and the only state an absorb accepts.
pub(crate) fn key_is_unbound<T>(key: ogeom_core::Key<T>) -> bool {
    key.scope() == ogeom_core::UNSCOPED && key.generation() == 0
}

/// A key shifted `offset` slots along, for landing parts in a live model.
///
/// The result stays unscoped: shifting says where an entry will sit, binding
/// says which arena it sits in, and the two are separate steps on purpose.
///
/// # Panics
///
/// If the shifted index exceeds `u32::MAX`, which
/// [`Arena::insert`](ogeom_core::Arena::insert) already refuses to reach.
#[allow(clippy::expect_used, reason = "documented panic; see # Panics")]
pub(crate) fn shifted_key<T>(key: ogeom_core::Key<T>, offset: u32) -> ogeom_core::Key<T> {
    ogeom_core::Key::from_parts(
        key.index()
            .checked_add(offset)
            .expect("arena exceeded u32::MAX slots"),
        key.generation(),
    )
}

/// One way of describing where an edge runs.
///
/// An edge holds several at once, and they must agree — see
/// [`EdgeData::same_parameter`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EdgeRepr {
    /// The edge as a curve in space.
    Curve3d {
        /// The curve.
        curve: CurveId,
        /// Where the curve sits.
        location: Location,
        /// The portion of the curve this edge covers.
        range: (f64, f64),
    },
    /// The edge as a curve in one surface's parameter space.
    PCurve {
        /// The curve in `(u, v)`.
        curve: PCurveId,
        /// The surface whose parameter space it lives in.
        surface: SurfaceId,
        /// Where that surface sits.
        location: Location,
        /// The portion of the curve this edge covers.
        range: (f64, f64),
    },
    /// The edge as a seam on a closed surface, needing two pcurves.
    ///
    /// A seam runs along a surface's closure — a cylinder's join, a sphere's
    /// date line — where the same points have two parameter values. One pcurve
    /// per side; a single one could not express both, and using only one leaves
    /// the face split open along the seam.
    Seam {
        /// The pcurve on the side the face's forward orientation sees.
        forward: PCurveId,
        /// The pcurve on the other side.
        reversed: PCurveId,
        /// The surface.
        surface: SurfaceId,
        /// Where that surface sits.
        location: Location,
        /// The portion both curves cover.
        range: (f64, f64),
    },
    /// The edge's path through one face's cached triangulation, as indices.
    ///
    /// What a renderer wants to draw a shared edge without hunting for
    /// coincident vertices: consecutive indices are consecutive mesh nodes
    /// along the edge, in the edge's own direction. Only meaningful while
    /// the named triangulation is the one stored; retessellating replaces
    /// both together.
    PolygonOnTriangulation {
        /// The face triangulation the indices point into.
        triangulation: TriangulationId,
        /// Node indices along the edge, in curve order.
        indices: Vec<u32>,
        /// Where the edge occurrence sits.
        location: Location,
    },
    /// A cached polyline approximation, with the deflection it was built to.
    ///
    /// Carried alongside the exact geometry rather than replacing it: display
    /// and coarse spatial queries want it, and rebuilding it on every frame
    /// would dominate their cost.
    Polyline {
        /// The points, in order along the edge.
        points: Vec<Point>,
        /// The curve parameter each point came from.
        ///
        /// Kept, not discarded. A face's cached triangulation has to place its
        /// boundary vertices where this polyline puts them, and it reaches them
        /// through its own pcurve — so it needs the parameters, not just the
        /// points. Without them the two caches drift and the stored mesh has
        /// gaps that the exact geometry does not.
        parameters: Vec<f64>,
        /// Where they sit.
        location: Location,
        /// The maximum distance from the exact curve.
        deflection: f64,
    },
}

impl EdgeRepr {
    /// Bind this representation's handles to the arenas that hold them.
    ///
    /// For reading a document back. Only the arena identifier changes; index
    /// and generation came from the file and are already right.
    pub(crate) fn rebind(&mut self, geometry: &GeometryScopes, datums: u32) {
        match self {
            Self::Curve3d {
                curve, location, ..
            } => {
                *curve = curve.with_scope(geometry.curves);
                *location = location.with_datum_scope(datums);
            }
            Self::PCurve {
                curve,
                surface,
                location,
                ..
            } => {
                *curve = curve.with_scope(geometry.pcurves);
                *surface = surface.with_scope(geometry.surfaces);
                *location = location.with_datum_scope(datums);
            }
            Self::Seam {
                forward,
                reversed,
                surface,
                location,
                ..
            } => {
                *forward = forward.with_scope(geometry.pcurves);
                *reversed = reversed.with_scope(geometry.pcurves);
                *surface = surface.with_scope(geometry.surfaces);
                *location = location.with_datum_scope(datums);
            }
            Self::Polyline { location, .. } => {
                *location = location.with_datum_scope(datums);
            }
            Self::PolygonOnTriangulation {
                triangulation,
                location,
                ..
            } => {
                *triangulation = triangulation.with_scope(geometry.triangulations);
                *location = location.with_datum_scope(datums);
            }
        }
    }

    /// Shift this representation's handles for landing in a live model.
    ///
    /// [`rebind`](Self::rebind)'s sibling for absorbing parts: indices move by
    /// where the source document's geometry and datums land in the target,
    /// and the handles stay unscoped for the binding pass that follows.
    pub(crate) fn shift(&mut self, geometry: &GeometryOffsets, datums: u32) {
        match self {
            Self::Curve3d {
                curve, location, ..
            } => {
                *curve = shifted_key(*curve, geometry.curves);
                *location = location.with_datum_offset(datums);
            }
            Self::PCurve {
                curve,
                surface,
                location,
                ..
            } => {
                *curve = shifted_key(*curve, geometry.pcurves);
                *surface = shifted_key(*surface, geometry.surfaces);
                *location = location.with_datum_offset(datums);
            }
            Self::Seam {
                forward,
                reversed,
                surface,
                location,
                ..
            } => {
                *forward = shifted_key(*forward, geometry.pcurves);
                *reversed = shifted_key(*reversed, geometry.pcurves);
                *surface = shifted_key(*surface, geometry.surfaces);
                *location = location.with_datum_offset(datums);
            }
            Self::Polyline { location, .. } => {
                *location = location.with_datum_offset(datums);
            }
            Self::PolygonOnTriangulation {
                triangulation,
                location,
                ..
            } => {
                *triangulation = shifted_key(*triangulation, geometry.triangulations);
                *location = location.with_datum_offset(datums);
            }
        }
    }

    /// Whether every handle here is unscoped and at generation zero — the
    /// state a reader leaves them in, and the only state an absorb accepts.
    pub(crate) fn is_unbound(&self) -> bool {
        let local_location = |location: &Location| {
            location
                .chain()
                .iter()
                .all(|&(datum, _)| key_is_unbound(datum))
        };
        match self {
            Self::Curve3d {
                curve, location, ..
            } => key_is_unbound(*curve) && local_location(location),
            Self::PCurve {
                curve,
                surface,
                location,
                ..
            } => key_is_unbound(*curve) && key_is_unbound(*surface) && local_location(location),
            Self::Seam {
                forward,
                reversed,
                surface,
                location,
                ..
            } => {
                key_is_unbound(*forward)
                    && key_is_unbound(*reversed)
                    && key_is_unbound(*surface)
                    && local_location(location)
            }
            Self::Polyline { location, .. } => local_location(location),
            Self::PolygonOnTriangulation {
                triangulation,
                location,
                ..
            } => key_is_unbound(*triangulation) && local_location(location),
        }
    }

    /// The surface this representation belongs to, if any.
    #[must_use]
    pub const fn surface(&self) -> Option<SurfaceId> {
        match self {
            Self::PCurve { surface, .. } | Self::Seam { surface, .. } => Some(*surface),
            Self::Curve3d { .. } | Self::Polyline { .. } | Self::PolygonOnTriangulation { .. } => {
                None
            }
        }
    }

    /// Where this representation sits.
    #[must_use]
    pub const fn location(&self) -> Option<&Location> {
        match self {
            Self::Curve3d { location, .. }
            | Self::PCurve { location, .. }
            | Self::Seam { location, .. }
            | Self::Polyline { location, .. }
            | Self::PolygonOnTriangulation { location, .. } => Some(location),
        }
    }

    /// The parameter range this representation covers, if it has one.
    #[must_use]
    pub const fn range(&self) -> Option<(f64, f64)> {
        match self {
            Self::Curve3d { range, .. } | Self::PCurve { range, .. } | Self::Seam { range, .. } => {
                Some(*range)
            }
            Self::Polyline { .. } | Self::PolygonOnTriangulation { .. } => None,
        }
    }

    /// Whether this is the edge's curve in space.
    #[must_use]
    pub const fn is_curve3d(&self) -> bool {
        matches!(self, Self::Curve3d { .. })
    }

    /// Whether this describes the edge in a surface's parameter space.
    #[must_use]
    pub const fn is_parametric(&self) -> bool {
        matches!(self, Self::PCurve { .. } | Self::Seam { .. })
    }
}

/// A vertex: a point and how far it may be from where it claims to be.
#[derive(Debug, Clone, PartialEq)]
pub struct VertexData {
    /// The position.
    pub point: Point,
    /// The radius within which the vertex lies.
    pub tolerance: Tolerance,
}

impl VertexData {
    /// A vertex at `point` with the minimum tolerance.
    #[must_use]
    pub fn new(point: Point) -> Self {
        Self {
            point,
            tolerance: Tolerance::MIN,
        }
    }

    /// A vertex with an explicit tolerance.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// tolerance is not finite and non-negative, or the point is non-finite.
    pub fn with_tolerance(point: Point, tolerance: f64) -> OgeomResult<Self> {
        if !point.is_finite() {
            ogeom_bail!(Construction, "vertex position is not finite");
        }
        Ok(Self {
            point,
            tolerance: Tolerance::new(tolerance)?,
        })
    }

    /// Widen this vertex's tolerance.
    pub fn widen(&mut self, to: Tolerance) {
        self.tolerance = self.tolerance.widen(to);
    }
}

/// An edge: a tolerance, a set of representations, and the flags that say
/// whether they agree.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeData {
    /// The radius within which the edge lies.
    pub tolerance: Tolerance,
    /// Every way this edge is described. See the module documentation.
    pub representations: SmallVec<[EdgeRepr; 3]>,
    /// Whether the representations agree on parameterization.
    ///
    /// A *claim*, and one that can be false — see [`EdgeData::same_parameter`].
    same_parameter: bool,
    /// Whether the edge has no length: a cone's apex, a sphere's pole.
    ///
    /// A degenerate edge still bounds a face in parameter space even though it
    /// covers no distance in space, which is why it exists rather than being
    /// dropped.
    pub degenerate: bool,
}

impl EdgeData {
    /// An edge with no representations yet.
    #[must_use]
    pub fn new() -> Self {
        Self {
            tolerance: Tolerance::MIN,
            representations: SmallVec::new(),
            same_parameter: true,
            degenerate: false,
        }
    }

    /// An edge on a curve in space.
    #[must_use]
    pub fn on_curve(curve: CurveId, location: Location, range: (f64, f64)) -> Self {
        let mut edge = Self::new();
        edge.representations.push(EdgeRepr::Curve3d {
            curve,
            location,
            range,
        });
        edge
    }

    /// Add a representation.
    ///
    /// Adding one invalidates the [`EdgeData::same_parameter`] claim: the new
    /// representation has not been shown to agree with the others. Re-establish
    /// it deliberately, with [`EdgeData::assert_same_parameter`].
    pub fn add(&mut self, repr: EdgeRepr) {
        self.representations.push(repr);
        self.same_parameter = false;
    }

    /// Whether every representation agrees on parameterization.
    ///
    /// The claim is that `curve3d(t)` and `surface(pcurve(t))` are the same
    /// point, within the edge's tolerance, for the same `t`. It matters because
    /// nearly every algorithm evaluates whichever representation is convenient
    /// and assumes the answer is interchangeable.
    ///
    /// It can be false — an imported edge whose pcurve was fitted independently
    /// of its 3D curve routinely is — which is why it is a flag to be checked
    /// rather than an invariant to be assumed.
    #[must_use]
    pub const fn same_parameter(&self) -> bool {
        self.same_parameter
    }

    /// Record that the representations have been checked and agree.
    ///
    /// Only for a caller that has actually verified it: evaluate each
    /// representation at the same parameters and confirm they land within the
    /// edge's tolerance of one another. Setting it without checking is how an
    /// edge ends up with a pcurve that does not follow its own curve, and
    /// nothing downstream will notice until a face fails to close.
    pub const fn assert_same_parameter(&mut self, agrees: bool) {
        self.same_parameter = agrees;
    }

    /// The representation on a given curve in space, if any.
    #[must_use]
    pub fn curve3d(&self) -> Option<&EdgeRepr> {
        self.representations.iter().find(|r| r.is_curve3d())
    }

    /// The representation in `surface`'s parameter space, if any.
    #[must_use]
    pub fn pcurve_on(&self, surface: SurfaceId) -> Option<&EdgeRepr> {
        self.representations
            .iter()
            .find(|r| r.surface() == Some(surface))
    }

    /// The representation in `surface`'s parameter space for an occurrence at
    /// `location`.
    ///
    /// One edge node can bound one face at more than one placement — the top
    /// and bottom of a prism are the same edge, moved — and those two
    /// occurrences run along different lines of the same parameter space. Asked
    /// by surface alone, the lookup returns whichever was attached first and
    /// both ends of the prism collapse onto one.
    ///
    /// Falls back to a representation attached without a placement, which is
    /// what every unplaced edge has and what keeps the simple case simple.
    #[must_use]
    pub fn pcurve_for(&self, surface: SurfaceId, location: &Location) -> Option<&EdgeRepr> {
        self.representations
            .iter()
            .find(|r| r.surface() == Some(surface) && r.location() == Some(location))
            .or_else(|| {
                self.representations.iter().find(|r| {
                    r.surface() == Some(surface) && r.location().is_some_and(Location::is_identity)
                })
            })
    }

    /// Every surface this edge has a pcurve on.
    #[must_use]
    pub fn parametric_surfaces(&self) -> Vec<SurfaceId> {
        self.representations
            .iter()
            .filter_map(EdgeRepr::surface)
            .collect()
    }

    /// Widen this edge's tolerance.
    pub fn widen(&mut self, to: Tolerance) {
        self.tolerance = self.tolerance.widen(to);
    }
}

impl Default for EdgeData {
    fn default() -> Self {
        Self::new()
    }
}

/// A face: a surface, where it sits, and how far the face may stray from it.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceData {
    /// The surface.
    pub surface: SurfaceId,
    /// Where the surface sits.
    pub location: Location,
    /// The radius within which the face lies.
    pub tolerance: Tolerance,
    /// Whether the face covers its surface's whole domain, with no trimming
    /// wires of its own.
    ///
    /// Worth knowing: a face with natural restriction needs no point-in-face
    /// classification at all, and that test is one of the costliest in a
    /// boolean.
    pub natural_restriction: bool,
    /// The cached triangulation of this face, if one has been built.
    ///
    /// A representation like a pcurve, not a replacement for the surface: it
    /// answers display and coarse queries, and is rebuilt when a finer
    /// deflection is asked for.
    pub triangulation: Option<TriangulationId>,
}

impl FaceData {
    /// A face on `surface`, trimmed by its own wires.
    #[must_use]
    pub fn new(surface: SurfaceId, location: Location) -> Self {
        Self {
            surface,
            location,
            tolerance: Tolerance::MIN,
            natural_restriction: false,
            triangulation: None,
        }
    }

    /// A face covering the whole of `surface`.
    #[must_use]
    pub fn natural(surface: SurfaceId, location: Location) -> Self {
        Self {
            natural_restriction: true,
            ..Self::new(surface, location)
        }
    }

    /// Widen this face's tolerance.
    pub fn widen(&mut self, to: Tolerance) {
        self.tolerance = self.tolerance.widen(to);
    }
}

/// What a topology node holds, beyond its children.
///
/// Stored inline in the node rather than in side tables keyed by handle: a
/// traversal that has the node has the data, without a second lookup, and there
/// is no way for the two to fall out of step.
#[derive(Debug, Clone, PartialEq)]
pub enum NodeData {
    /// A vertex.
    Vertex(VertexData),
    /// An edge.
    Edge(Box<EdgeData>),
    /// A face.
    Face(Box<FaceData>),
    /// A wire, shell, solid, compsolid or compound: structure, no geometry of
    /// its own.
    Container,
}

impl NodeData {
    /// The tolerance this node carries, if it carries one.
    ///
    /// Containers have none: a wire's uncertainty is that of the edges in it,
    /// and inventing a separate number for it would be a second source of truth.
    #[must_use]
    pub fn tolerance(&self) -> Option<Tolerance> {
        match self {
            Self::Vertex(v) => Some(v.tolerance),
            Self::Edge(e) => Some(e.tolerance),
            Self::Face(f) => Some(f.tolerance),
            Self::Container => None,
        }
    }

    /// Widen this node's tolerance, if it has one.
    pub fn widen(&mut self, to: Tolerance) {
        match self {
            Self::Vertex(v) => v.widen(to),
            Self::Edge(e) => e.widen(to),
            Self::Face(f) => f.widen(to),
            Self::Container => {}
        }
    }

    /// The vertex data, if this is a vertex.
    #[must_use]
    pub const fn as_vertex(&self) -> Option<&VertexData> {
        match self {
            Self::Vertex(v) => Some(v),
            _ => None,
        }
    }

    /// The edge data, if this is an edge.
    #[must_use]
    pub const fn as_edge(&self) -> Option<&EdgeData> {
        match self {
            Self::Edge(e) => Some(e),
            _ => None,
        }
    }

    /// The face data, if this is a face.
    #[must_use]
    pub const fn as_face(&self) -> Option<&FaceData> {
        match self {
            Self::Face(f) => Some(f),
            _ => None,
        }
    }
}

/// Check the containment rule between a bounding entity and what it bounds.
///
/// `docs/DATA_MODEL.md` §5. A vertex must be at least as uncertain as the edge
/// it caps, and an edge at least as uncertain as the face it borders. If it
/// were not, the boundary would not reliably lie on the thing it bounds, and
/// every containment test built on it would be answering about geometry that
/// does not quite meet.
///
/// # Errors
///
/// [`OgeomError::Invariant`](ogeom_core::OgeomError::Invariant) if `bounding` is tighter
/// than `bounded`.
pub fn check_containment(bounding: Tolerance, bounded: Tolerance) -> OgeomResult<()> {
    ogeom_core::check_containment(bounding, bounded)
}

/// Widen `bounding` just enough to satisfy the containment rule against
/// `bounded`.
///
/// The sanctioned repair: tolerances only ever grow, so restoring the rule
/// means raising the bounding entity rather than lowering what it bounds.
#[must_use]
pub fn enforce_containment(bounding: Tolerance, bounded: Tolerance) -> Tolerance {
    bounding.widen(bounded)
}

/// Whether a parameter range is usable.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the range is
/// empty or non-finite.
pub fn check_range(range: (f64, f64), tol: Tolerances) -> OgeomResult<()> {
    let (a, b) = range;
    if !a.is_finite() || !b.is_finite() || b <= a + tol.parametric() {
        ogeom_bail!(
            Construction,
            "parameter range [{a}, {b}] is empty or non-finite"
        );
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_geom::{CircleCurve, LineCurve, PlaneSurface};
    use ogeom_math::{Circle, Direction, Frame, Plane};

    const T: Tolerances = Tolerances::millimetres();

    fn store() -> (GeometryStore, CurveId, SurfaceId, PCurveId) {
        let mut s = GeometryStore::new();
        let curve = s.add_curve(
            LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
                .unwrap()
                .into(),
        );
        let surface = s.add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let pcurve = s.add_pcurve(
            ogeom_geom::Line2d::segment(
                ogeom_math::Point2::ORIGIN,
                ogeom_math::Point2::new(10.0, 0.0),
                T,
            )
            .unwrap()
            .into(),
        );
        (s, curve, surface, pcurve)
    }

    #[test]
    fn the_geometry_store_hands_back_what_it_was_given() {
        let (s, curve, surface, pcurve) = store();
        assert!(s.curve(curve).is_some());
        assert!(s.surface(surface).is_some());
        assert!(s.pcurve(pcurve).is_some());
        assert_eq!(s.counts(), (1, 1, 1));

        // Handles are typed, so a curve handle cannot address a surface. That
        // is a compile-time guarantee rather than a runtime check.
        // A plain loop, not a lazy iterator: `map(..).next_back()` would run
        // the closure once, so nothing would actually be inserted past index 0
        // and the handle would resolve after all.
        let mut other = GeometryStore::new();
        let mut beyond = curve;
        for _ in 0..5 {
            beyond = other.add_curve(
                LineCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0), T)
                    .unwrap()
                    .into(),
            );
        }
        assert!(
            s.curve(beyond).is_none(),
            "a handle past the end does not resolve"
        );
    }

    #[test]
    fn tolerances_start_at_the_minimum_and_only_widen() {
        let mut v = VertexData::new(Point::ORIGIN);
        assert_eq!(v.tolerance, Tolerance::MIN);

        let wide = Tolerance::new(1e-3).unwrap();
        v.widen(wide);
        assert_eq!(v.tolerance, wide);

        // Widening to something tighter leaves it alone: tolerances never
        // shrink, since narrowing one asserts an accuracy the geometry does not
        // have.
        v.widen(Tolerance::new(1e-9).unwrap());
        assert_eq!(v.tolerance, wide);
    }

    #[test]
    fn degenerate_vertex_data_is_refused() {
        assert!(VertexData::with_tolerance(Point::ORIGIN, -1.0).is_err());
        assert!(VertexData::with_tolerance(Point::ORIGIN, f64::NAN).is_err());
        assert!(VertexData::with_tolerance(Point::new(f64::INFINITY, 0.0, 0.0), 1e-6).is_err());
        assert!(VertexData::with_tolerance(Point::ORIGIN, 1e-3).is_ok());
    }

    #[test]
    fn the_containment_rule_holds_downward_and_is_repaired_upward() {
        let fine = Tolerance::new(1e-6).unwrap();
        let coarse = Tolerance::new(1e-3).unwrap();

        assert!(
            check_containment(coarse, fine).is_ok(),
            "vertex coarser than edge"
        );
        assert!(
            check_containment(fine, coarse).is_err(),
            "and not the other way"
        );

        // The repair raises the bounding entity, never lowers what it bounds.
        let repaired = enforce_containment(fine, coarse);
        assert_eq!(repaired, coarse);
        assert!(check_containment(repaired, coarse).is_ok());
    }

    #[test]
    fn an_edge_holds_several_representations_at_once() {
        // The whole point of §6: one edge, a curve in space and a pcurve per
        // adjacent face. A single curve could not serve two differently
        // parameterized surfaces.
        let (mut s, curve, surface, pcurve) = store();
        let other_surface =
            s.add_surface(PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Y)).into());
        let other_pcurve = s.add_pcurve(
            ogeom_geom::Line2d::segment(
                ogeom_math::Point2::ORIGIN,
                ogeom_math::Point2::new(0.0, 10.0),
                T,
            )
            .unwrap()
            .into(),
        );

        let mut edge = EdgeData::on_curve(curve, Location::identity(), (0.0, 10.0));
        edge.add(EdgeRepr::PCurve {
            curve: pcurve,
            surface,
            location: Location::identity(),
            range: (0.0, 10.0),
        });
        edge.add(EdgeRepr::PCurve {
            curve: other_pcurve,
            surface: other_surface,
            location: Location::identity(),
            range: (0.0, 10.0),
        });

        assert_eq!(edge.representations.len(), 3);
        assert!(edge.curve3d().is_some());
        assert!(edge.pcurve_on(surface).is_some());
        assert!(edge.pcurve_on(other_surface).is_some());
        assert_eq!(edge.parametric_surfaces().len(), 2);
    }

    #[test]
    fn adding_a_representation_withdraws_the_same_parameter_claim() {
        // The claim is that every representation lands on the same point for
        // the same parameter. A newly added one has not been shown to, and
        // quietly keeping the claim is how an edge ends up with a pcurve that
        // does not follow its own curve.
        let (_, curve, surface, pcurve) = store();
        let mut edge = EdgeData::on_curve(curve, Location::identity(), (0.0, 10.0));
        assert!(
            edge.same_parameter(),
            "a lone curve trivially agrees with itself"
        );

        edge.add(EdgeRepr::PCurve {
            curve: pcurve,
            surface,
            location: Location::identity(),
            range: (0.0, 10.0),
        });
        assert!(
            !edge.same_parameter(),
            "the new representation is unverified"
        );

        edge.assert_same_parameter(true);
        assert!(edge.same_parameter());
    }

    #[test]
    fn a_seam_carries_two_pcurves_because_one_cannot_express_both_sides() {
        // On a closed surface the same points have two parameter values. A
        // single pcurve names one of them, which leaves the face open along the
        // seam.
        let mut s = GeometryStore::new();
        let cylinder = s.add_surface(
            ogeom_geom::CylinderSurface::new(
                ogeom_math::Cylinder::new(Frame::WORLD, 2.0, T).unwrap(),
                (0.0, 5.0),
            )
            .unwrap()
            .into(),
        );
        let at_zero = s.add_pcurve(
            ogeom_geom::Line2d::segment(
                ogeom_math::Point2::ORIGIN,
                ogeom_math::Point2::new(0.0, 5.0),
                T,
            )
            .unwrap()
            .into(),
        );
        let at_tau = s.add_pcurve(
            ogeom_geom::Line2d::segment(
                ogeom_math::Point2::new(core::f64::consts::TAU, 0.0),
                ogeom_math::Point2::new(core::f64::consts::TAU, 5.0),
                T,
            )
            .unwrap()
            .into(),
        );

        let mut edge = EdgeData::new();
        edge.add(EdgeRepr::Seam {
            forward: at_zero,
            reversed: at_tau,
            surface: cylinder,
            location: Location::identity(),
            range: (0.0, 5.0),
        });

        let seam = &edge.representations[0];
        assert!(seam.is_parametric());
        assert_eq!(seam.surface(), Some(cylinder));
        assert_eq!(seam.range(), Some((0.0, 5.0)));
        assert!(matches!(seam, EdgeRepr::Seam { .. }));
    }

    #[test]
    fn a_polyline_representation_records_what_it_was_built_to() {
        // Without the deflection the cache cannot be judged: a polyline good
        // enough for a thumbnail is not good enough for a machining path, and
        // nothing else recorded says which this is.
        let repr = EdgeRepr::Polyline {
            points: vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)],
            parameters: vec![0.0, 1.0],
            location: Location::identity(),
            deflection: 1e-3,
        };
        assert!(!repr.is_curve3d() && !repr.is_parametric());
        assert_eq!(repr.surface(), None);
        assert_eq!(repr.range(), None, "a polyline has no parameter range");
    }

    #[test]
    fn a_degenerate_edge_is_marked_rather_than_dropped() {
        // A cone's apex has no length in space but still bounds the face in
        // parameter space. Dropping it leaves the face's boundary open.
        let mut edge = EdgeData::new();
        edge.degenerate = true;
        assert!(edge.degenerate);
        assert!(edge.representations.is_empty());
    }

    #[test]
    fn a_natural_face_covers_its_whole_surface() {
        let (_, _, surface, _) = store();
        let trimmed = FaceData::new(surface, Location::identity());
        let whole = FaceData::natural(surface, Location::identity());
        assert!(!trimmed.natural_restriction);
        assert!(whole.natural_restriction);
        assert_eq!(whole.tolerance, Tolerance::MIN);
    }

    #[test]
    fn node_data_exposes_only_what_it_holds() {
        let vertex = NodeData::Vertex(VertexData::new(Point::ORIGIN));
        let edge = NodeData::Edge(Box::default());
        let container = NodeData::Container;

        assert!(vertex.as_vertex().is_some());
        assert!(vertex.as_edge().is_none());
        assert!(edge.as_edge().is_some());
        assert!(edge.as_face().is_none());

        // A container has no tolerance of its own: a wire's uncertainty is that
        // of its edges, and a second number for it would be a second source of
        // truth to keep in step.
        assert_eq!(container.tolerance(), None);
        assert_eq!(vertex.tolerance(), Some(Tolerance::MIN));
    }

    #[test]
    fn widening_a_container_is_a_no_op_rather_than_an_error() {
        // Traversal widens whatever it walks over; a container simply has
        // nothing to widen, and making that an error would push the case
        // analysis onto every caller.
        let mut container = NodeData::Container;
        container.widen(Tolerance::new(1e-3).unwrap());
        assert_eq!(container.tolerance(), None);
    }

    #[test]
    fn empty_parameter_ranges_are_refused() {
        assert!(check_range((0.0, 1.0), T).is_ok());
        assert!(check_range((1.0, 0.0), T).is_err());
        assert!(check_range((1.0, 1.0), T).is_err());
        assert!(check_range((0.0, f64::NAN), T).is_err());
    }

    #[test]
    fn a_circular_edge_can_carry_its_own_curve() {
        let mut s = GeometryStore::new();
        let circle =
            s.add_curve(CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()).into());
        let edge = EdgeData::on_curve(circle, Location::identity(), (0.0, core::f64::consts::TAU));
        assert!(edge.curve3d().is_some());
        assert_eq!(
            edge.curve3d().unwrap().range(),
            Some((0.0, core::f64::consts::TAU))
        );
        assert!(edge.parametric_surfaces().is_empty());
    }
}
