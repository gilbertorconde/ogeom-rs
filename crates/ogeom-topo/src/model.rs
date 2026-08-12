//! The model: the arenas a shape's handles refer into, and the builder that is
//! the only way to mutate them.
//!
//! A [`Model`] owns the topology nodes, the placement datums and the geometry.
//! A [`Shape`] is meaningless without one — its handles index into these arenas
//! and nothing else (`docs/DATA_MODEL.md` §11).
//!
//! # One mutation path
//!
//! Every structural change goes through [`Model`]'s builder methods. That is
//! not ceremony: the invariants of `docs/DATA_MODEL.md` — a wire holds edges, a
//! face's tolerance does not exceed its edges', a node's kind matches its data
//! — are checkable in one place only if there is one place. Handing out
//! `&mut TShape` would scatter them across every caller, and the failures they
//! guard against are silent ones.

use std::collections::HashMap;

use ogeom_core::{
    Arena, EntityId, OgeomResult, OpId, Provenance, ProvenanceTable, Role, Tolerance, Tolerances,
    ogeom_bail,
};
use ogeom_math::{Point, Transform};

use crate::entity::{EdgeData, EdgeRepr, FaceData, NodeData, VertexData};
use crate::location::{DatumId, DatumStore, Location};
use crate::shape::{Orientation, Shape, ShapeType, TShape, TShapeId};

pub use crate::entity::GeometryStore;

/// A document: topology, placements, geometry, and where every entity came
/// from.
#[derive(Debug, Clone, Default)]
pub struct Model {
    nodes: Arena<TShape>,
    datums: DatumStore,
    geometry: GeometryStore,
    provenance: ProvenanceTable,
    identity: HashMap<TShapeId, EntityId>,
    current_op: OpId,
    tolerances: Tolerances,
}

impl Model {
    /// An empty model, in millimetres.
    #[must_use]
    pub fn new() -> Self {
        Self::with_tolerances(Tolerances::millimetres())
    }

    /// An empty model at a given unit scale.
    ///
    /// A document has a scale, and it is the document's rather than each
    /// call's: a model authored in metres does not become a model in
    /// millimetres because one caller passed the default. Algorithms still take
    /// a [`Tolerances`] argument — that is deliberate, since a caller may want
    /// to work coarser or finer than the document's own setting for one
    /// operation — but the document now says what it was built at, so a
    /// mismatch is visible rather than assumed away.
    #[must_use]
    pub fn with_tolerances(tolerances: Tolerances) -> Self {
        Self {
            nodes: Arena::new(),
            datums: DatumStore::new(),
            geometry: GeometryStore::new(),
            provenance: ProvenanceTable::new(),
            identity: HashMap::new(),
            current_op: OpId(0),
            tolerances,
        }
    }

    /// The tolerances this document was built at.
    #[must_use]
    pub const fn tolerances(&self) -> Tolerances {
        self.tolerances
    }

    /// Assemble a model from parts read back from a file.
    ///
    /// The one way into a [`Model`] that does not go through its builders, and
    /// it exists for one reason: a builder *mints* an identity for every node
    /// it makes (`docs/DATA_MODEL.md` §8). A document rebuilt through the
    /// builders is therefore a different document from the one that was
    /// written — every [`EntityId`] renumbered, every provenance record
    /// replaced by a fresh `Primitive` one — and every reference into it,
    /// which is the thing provenance exists to keep alive, is dead. Reading a
    /// file has to reproduce the document it describes, identities and all.
    ///
    /// This is not a hole in "the builder is the sole mutation path". Nothing
    /// here mutates an existing model; it assembles a new one, and it checks
    /// the structural invariants the builders check before handing it back —
    /// so a corrupt file is an error, not a model that answers wrongly.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a node names a
    /// child, a datum, a piece of geometry or an identity that is not there;
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a node's
    /// data does not match its kind, or a child is of the wrong kind for its
    /// parent.
    pub fn from_parts(parts: ModelParts) -> OgeomResult<Self> {
        let mut model = Self::with_tolerances(parts.tolerances);
        model.current_op = parts.current_op;
        // Absorbing into an empty model is exactly restoration: every offset
        // is zero, so each handle keeps the index and generation the file
        // recorded, which is what makes a document's handles survive a round
        // trip.
        model.absorb_core(parts)?;
        Ok(model)
    }

    /// Absorb another document's parts into this model.
    ///
    /// [`Model::from_parts`] lands a document in a *fresh* model; this lands
    /// one in a model that already has things in it — the operation behind
    /// bringing a serialized tool body into a live document so a boolean can
    /// use it. Everything the parts carry is appended: nodes, datums,
    /// geometry, provenance and identities all keep their relative structure,
    /// shifted past what the model already holds. The result behaves exactly
    /// as if it had been built here, because after the shift it is
    /// indistinguishable from having been.
    ///
    /// `roots` are the shapes the source document named, in the unbound state
    /// a reader leaves them; they come back bound to this model. The returned
    /// [`Absorbed::entities`] table says where every source identity landed,
    /// which is what a caller holding references against the source document
    /// resolves them through.
    ///
    /// Three deliberate refusals, each an error rather than a guess:
    ///
    /// - **Units.** A document authored at another scale is refused, not
    ///   rescaled — rescaling is a real feature with real decisions in it,
    ///   and silently absorbing metres into millimetres is a wrong model.
    /// - **Bound handles.** Parts whose handles already name an arena did not
    ///   come from a reader; absorbing them would alias whatever those
    ///   handles meant elsewhere. Serialization is the one road in.
    /// - **A model with holes.** Absorbing appends by offset, which is only
    ///   sound while the target's arenas have only ever been appended to.
    ///   Nothing in this crate removes, so this cannot trigger today; it is
    ///   checked so a future that removes gets an error, not aliasing.
    ///
    /// The current operation is left alone: absorb mints no identities, it
    /// transplants a table, and the absorbed provenance keeps its source
    /// [`OpId`]s verbatim — meaningful in the source document's rebuild, kept
    /// because renumbering them would orphan the source's own references.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for
    /// the three refusals above;
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the parts
    /// do not describe themselves — a handle that does not resolve, a
    /// derivation from an identity no entry issued, a root naming a node that
    /// is not there.
    ///
    /// On an error past the up-front gates, the model's prior contents are
    /// untouched and still fully usable: an absorbed subgraph is
    /// self-contained — its shifted handles cannot reach below the append
    /// line — and no identity is committed until every check has passed. What
    /// a failed absorb can leave behind is unreachable appended entries,
    /// which cost memory and mean nothing.
    pub fn absorb(&mut self, parts: ModelParts, roots: &[Shape]) -> OgeomResult<Absorbed> {
        #[allow(clippy::float_cmp, reason = "scales are copied, never computed")]
        if parts.tolerances.scale() != self.tolerances.scale() {
            ogeom_bail!(
                Construction,
                "these parts were authored at {} mm per unit and this model at \
                 {}; absorbing across scales needs a rescale, which is its own \
                 operation",
                parts.tolerances.scale(),
                self.tolerances.scale()
            );
        }
        if !self.nodes.is_dense() || !self.datums.is_dense() || !self.geometry.is_dense() {
            ogeom_bail!(
                Construction,
                "absorb appends by offset, and this model's arenas have holes; \
                 something removed entries, which nothing in this crate does"
            );
        }
        Self::check_parts_unbound(&parts, roots)?;

        let absorbed_entities = parts.provenance.len() as u64;
        let (node_offset, datum_offset, entity_offset) = self.absorb_core(parts)?;

        let shapes = roots
            .iter()
            .map(|root| self.bind(&root.shifted(node_offset, datum_offset)))
            .collect::<OgeomResult<Vec<Shape>>>()?;
        let entities = (1..=absorbed_entities)
            .filter_map(|raw| {
                Some((
                    EntityId::from_raw(raw)?,
                    EntityId::from_raw(raw + entity_offset)?,
                ))
            })
            .collect();
        Ok(Absorbed { shapes, entities })
    }

    /// Append parts onto this model's arenas, shifting every handle past what
    /// is already here. The shared engine of [`Model::from_parts`] (offsets
    /// all zero) and [`Model::absorb`].
    ///
    /// Returns `(node, datum, entity)` offsets — where the parts landed.
    fn absorb_core(&mut self, parts: ModelParts) -> OgeomResult<(u32, u32, u64)> {
        let ModelParts {
            mut nodes,
            datums,
            geometry,
            provenance,
            identity,
            current_op: _,
            tolerances: _,
        } = parts;

        // Every id an entry names must already have been issued, which is
        // what makes the derivation graph acyclic by construction. Checked
        // against the parts' own issue order, before any shift.
        for (issued, entry) in provenance.iter().enumerate() {
            for source in entry.inputs() {
                if source.get() > issued as u64 {
                    ogeom_bail!(
                        Dangling,
                        "an entity is derived from identity {}, which no entry \
                         before it issued",
                        source.get()
                    );
                }
            }
        }

        let node_offset = crate::entity::arena_len(&self.nodes);
        let datum_offset = u32::try_from(self.datums.len()).unwrap_or(u32::MAX);
        let entity_offset = self.provenance.len() as u64;
        let geometry_offsets = self.geometry.append(geometry);

        for node in &mut nodes {
            for child in node.children_mut() {
                *child = child.shifted(node_offset, datum_offset);
            }
            match node.data_mut() {
                NodeData::Edge(edge) => {
                    for repr in &mut edge.representations {
                        repr.shift(&geometry_offsets, datum_offset);
                    }
                }
                NodeData::Face(face) => {
                    face.surface =
                        crate::entity::shifted_key(face.surface, geometry_offsets.surfaces);
                    face.triangulation = face.triangulation.map(|mesh| {
                        crate::entity::shifted_key(mesh, geometry_offsets.triangulations)
                    });
                    face.location = face.location.with_datum_offset(datum_offset);
                }
                NodeData::Vertex(_) | NodeData::Container => {}
            }
        }

        for datum in datums {
            self.datums.insert(datum);
        }
        for mut entry in provenance {
            if let Provenance::Derived { from, .. } = &mut entry {
                for source in from.iter_mut() {
                    let Some(shifted) = EntityId::from_raw(source.get() + entity_offset) else {
                        ogeom_bail!(Construction, "an entity id overflowed in the shift");
                    };
                    *source = shifted;
                }
            }
            self.provenance.record(entry);
        }
        for node in nodes {
            self.nodes.insert(node);
        }

        // Every handle in `parts` was rebuilt by a reader that had no arenas
        // to bind them to, so they name no arena at all and resolve nowhere.
        // Bind the appended subrange now that the arenas exist; what was here
        // before is already bound.
        self.bind_handles(node_offset);
        let identity: Vec<(TShapeId, EntityId)> = identity
            .into_iter()
            .map(|(node, entity)| {
                let Some(shifted) = EntityId::from_raw(entity.get() + entity_offset) else {
                    ogeom_bail!(Construction, "an entity id overflowed in the shift");
                };
                Ok((
                    crate::entity::shifted_key(node, node_offset).with_scope(self.nodes.scope()),
                    shifted,
                ))
            })
            .collect::<OgeomResult<_>>()?;

        self.check_restored(&identity, node_offset)?;
        for (node, entity) in identity {
            self.identity.insert(node, entity);
        }
        Ok((node_offset, datum_offset, entity_offset))
    }

    /// Refuse parts whose handles are already bound to an arena or carry a
    /// non-zero generation — the state no reader produces, and the state an
    /// offset shift would silently mangle.
    fn check_parts_unbound(parts: &ModelParts, roots: &[Shape]) -> OgeomResult<()> {
        use crate::entity::key_is_unbound;

        let local_location = |location: &Location| {
            location
                .chain()
                .iter()
                .all(|&(datum, _)| key_is_unbound(datum))
        };
        let local_shape =
            |shape: &Shape| key_is_unbound(shape.node()) && local_location(shape.location());

        let mut sound = parts.identity.iter().all(|&(node, _)| key_is_unbound(node))
            && roots.iter().all(local_shape);
        for node in &parts.nodes {
            sound = sound && node.children().iter().all(local_shape);
            match node.data() {
                NodeData::Edge(edge) => {
                    sound = sound && edge.representations.iter().all(EdgeRepr::is_unbound);
                }
                NodeData::Face(face) => {
                    sound = sound
                        && key_is_unbound(face.surface)
                        && face.triangulation.is_none_or(key_is_unbound)
                        && local_location(&face.location);
                }
                NodeData::Vertex(_) | NodeData::Container => {}
            }
        }
        if !sound {
            ogeom_bail!(
                Construction,
                "these parts carry handles already bound to an arena, or at a \
                 recycled generation; absorb takes parts exactly as a reader \
                 rebuilt them, and serialization is the one road in"
            );
        }
        Ok(())
    }

    /// Bind an unscoped shape to this model.
    ///
    /// A handle rebuilt by a reader names no arena, so it resolves nowhere
    /// until it is told which document it belongs to. This is how a reader says
    /// so — and it verifies the answer, so a file naming a node that is not
    /// there is an error rather than a shape that fails mysteriously later.
    ///
    /// It will not re-home a shape that already belongs to *another* model.
    /// That is exactly the mistake scoping exists to catch, and quietly
    /// relabelling it would hand back a shape that resolves and answers about
    /// the wrong entity.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the shape
    /// already belongs to a different model;
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it does not resolve
    /// here once bound.
    pub fn bind(&self, shape: &Shape) -> OgeomResult<Shape> {
        if shape.node().scope() != ogeom_core::UNSCOPED && !self.nodes.issued(shape.node()) {
            ogeom_bail!(
                Construction,
                "this shape belongs to another model; binding it here would \
                 make it resolve and answer about a different entity"
            );
        }
        let bound = shape.rebound(self.nodes.scope(), self.datums.scope());
        if self.nodes.get(bound.node()).is_none() {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        }
        Ok(bound)
    }

    /// Bind a bare location to this model's datum store.
    ///
    /// The persistence path's sibling of [`Model::bind`]: a location read
    /// from a file names datum handles that are unscoped until the store
    /// that holds them exists. Binding checks them too — a chain naming a
    /// datum not in this model is an error here rather than wherever it is
    /// first resolved.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the chain
    /// names a datum not in this model.
    pub fn bind_location(&self, location: &Location) -> OgeomResult<Location> {
        let bound = location.with_datum_scope(self.datums.scope());
        for &(datum, _) in bound.chain() {
            if self.datums.get(datum).is_none() {
                ogeom_bail!(Dangling, "location refers to a datum not in this model");
            }
        }
        Ok(bound)
    }

    /// Bind every handle in freshly restored nodes to the arena that holds it.
    ///
    /// `from` bounds the pass to nodes at that index and past it: restoration
    /// binds everything from zero, absorption only what it appended.
    fn bind_handles(&mut self, from: u32) {
        let nodes = self.nodes.scope();
        let datums = self.datums.scope();
        let geometry = self.geometry.scopes();

        for (_, node) in self.nodes.iter_mut().filter(|(id, _)| id.index() >= from) {
            for child in node.children_mut() {
                *child = child.rebound(nodes, datums);
            }
            match node.data_mut() {
                NodeData::Edge(edge) => {
                    for repr in &mut edge.representations {
                        repr.rebind(&geometry, datums);
                    }
                }
                NodeData::Face(face) => {
                    face.surface = face.surface.with_scope(geometry.surfaces);
                    face.triangulation = face
                        .triangulation
                        .map(|mesh| mesh.with_scope(geometry.triangulations));
                    face.location = face.location.with_datum_scope(datums);
                }
                NodeData::Vertex(_) | NodeData::Container => {}
            }
        }
    }

    /// Verify that restored nodes' handles all resolve and their children are
    /// of the kinds their parents admit.
    ///
    /// `from` bounds the pass the way [`Model::bind_handles`]' does. Children
    /// still resolve through the full arena, so nothing is under-checked — an
    /// absorbed subgraph is self-contained and can only point at itself.
    fn check_restored(&self, identity: &[(TShapeId, EntityId)], from: u32) -> OgeomResult<()> {
        for (id, node) in self.nodes.iter().filter(|(id, _)| id.index() >= from) {
            let kind = node.kind();
            match (kind, node.data()) {
                (ShapeType::Vertex, NodeData::Vertex(_))
                | (ShapeType::Edge, NodeData::Edge(_))
                | (ShapeType::Face, NodeData::Face(_)) => {}
                (
                    ShapeType::Wire
                    | ShapeType::Shell
                    | ShapeType::Solid
                    | ShapeType::CompSolid
                    | ShapeType::Compound,
                    NodeData::Container,
                ) => {}
                (kind, data) => {
                    ogeom_bail!(Construction, "node {id:?} is a {kind:?} and holds {data:?}")
                }
            }
            self.check_node_geometry(id, node)?;

            // A compound may hold anything; everything else admits exactly one
            // kind of child, which is what makes traversal's assumptions safe.
            let expected = kind.child_type();
            for child in node.children() {
                let Some(below) = self.nodes.get(child.node()) else {
                    ogeom_bail!(Dangling, "node {id:?} names a child that is not there");
                };
                if let Some(expected) = expected
                    && kind != ShapeType::Compound
                    && below.kind() != expected
                {
                    ogeom_bail!(
                        Construction,
                        "a {kind:?} takes {expected:?} children; node {id:?} \
                         names a {:?}",
                        below.kind()
                    );
                }
                self.check_location(child.location())?;
            }
        }
        for (node, entity) in identity {
            if self.nodes.get(*node).is_none() {
                ogeom_bail!(Dangling, "an identity is bound to a node that is not there");
            }
            if entity.get() > self.provenance.len() as u64 {
                ogeom_bail!(
                    Dangling,
                    "node {node:?} claims identity {}, which was never issued",
                    entity.get()
                );
            }
        }
        Ok(())
    }

    /// Verify that a node's geometry handles resolve.
    fn check_node_geometry(&self, id: TShapeId, node: &TShape) -> OgeomResult<()> {
        match node.data() {
            NodeData::Edge(data) => {
                for repr in &data.representations {
                    if let Some(location) = repr.location() {
                        self.check_location(location)?;
                    }
                    if !self.geometry.holds(repr) {
                        ogeom_bail!(
                            Dangling,
                            "edge {id:?} names geometry that is not in this model"
                        );
                    }
                }
            }
            NodeData::Face(data) => {
                self.check_location(&data.location)?;
                if self.geometry.surface(data.surface).is_none() {
                    ogeom_bail!(Dangling, "face {id:?} names a surface that is not there");
                }
                if let Some(mesh) = data.triangulation
                    && self.geometry.triangulation(mesh).is_none()
                {
                    ogeom_bail!(
                        Dangling,
                        "face {id:?} names a triangulation that is not there"
                    );
                }
            }
            NodeData::Vertex(_) | NodeData::Container => {}
        }
        Ok(())
    }

    /// Verify that every datum a placement names is interned.
    fn check_location(&self, location: &Location) -> OgeomResult<()> {
        for &(datum, _) in location.chain() {
            if self.datums.get(datum).is_none() {
                ogeom_bail!(Dangling, "a placement names a datum that is not there");
            }
        }
        Ok(())
    }

    /// Begin a new operation, and return its identifier.
    ///
    /// Every node created from here on is attributed to it until the next call.
    /// The counter is deterministic — the third operation in a rebuild is
    /// `OpId(3)` every time — which is what lets provenance survive a parameter
    /// change (`docs/DATA_MODEL.md` §8).
    pub const fn begin_operation(&mut self) -> OpId {
        self.current_op = OpId(self.current_op.0 + 1);
        self.current_op
    }

    /// The operation nodes are currently attributed to.
    #[must_use]
    pub const fn current_operation(&self) -> OpId {
        self.current_op
    }

    /// The stable identity of a shape's node.
    ///
    /// Distinct from its arena handle: the handle says where the data is and
    /// dies when the shape is rebuilt, while this says what the entity *is* and
    /// survives.
    #[must_use]
    pub fn identity_of(&self, shape: &Shape) -> Option<EntityId> {
        self.identity.get(&shape.node()).copied()
    }

    /// Where a shape's node came from.
    #[must_use]
    pub fn provenance_of(&self, shape: &Shape) -> Option<&Provenance> {
        self.provenance.get(self.identity_of(shape)?)
    }

    /// The provenance table.
    #[must_use]
    pub const fn provenance(&self) -> &ProvenanceTable {
        &self.provenance
    }

    /// Trace a shape back to the entities it ultimately came from.
    ///
    /// How a reference into a rebuilt model is resolved: find what the user
    /// originally picked, then find what that became.
    #[must_use]
    pub fn roots_of(&self, shape: &Shape) -> Vec<EntityId> {
        self.identity_of(shape)
            .map(|id| self.provenance.roots(id))
            .unwrap_or_default()
    }

    /// The shape carrying a given identity, if this document has one.
    ///
    /// The inverse of [`Model::identity_of`], and the answer to "I kept a
    /// reference and the document has been saved and reloaded since". A raw
    /// [`Shape`] cannot survive that — the reloaded document is a new set of
    /// arenas and [`Model::bind`] refuses a handle from another one, on
    /// purpose. An [`EntityId`] can, because it names *what the entity is*
    /// rather than where it sits (`docs/DATA_MODEL.md` §8), and that is the
    /// whole reason it exists.
    ///
    /// Returns the shape in its default placement and orientation. A caller
    /// that wants a particular occurrence explores from here.
    #[must_use]
    pub fn shape_of(&self, id: EntityId) -> Option<Shape> {
        self.identity
            .iter()
            .find(|(_, entity)| **entity == id)
            .map(|(node, _)| Shape::of(*node))
    }

    /// Record that a node was derived from other entities.
    ///
    /// Overwrites the `Primitive` attribution a builder assigns by default.
    /// An operation that splits or reshapes existing topology calls this, and
    /// what it records is what a later rebuild will match against.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape does not
    /// resolve in this model.
    pub fn set_derived(
        &mut self,
        shape: &Shape,
        from: &[Shape],
        role: Role,
    ) -> OgeomResult<EntityId> {
        if self.node(shape).is_none() {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        }
        let sources: Vec<EntityId> = from.iter().filter_map(|s| self.identity_of(s)).collect();
        let id = self.provenance.derived(self.current_op, sources, role);
        self.identity.insert(shape.node(), id);
        Ok(id)
    }

    /// Record a node's identity as it is created.
    fn record_primitive(&mut self, node: TShapeId, role: Role) {
        let id = self.provenance.primitive(self.current_op, role);
        self.identity.insert(node, id);
    }

    /// The placement datums.
    #[must_use]
    pub const fn datums(&self) -> &DatumStore {
        &self.datums
    }

    /// The geometry.
    #[must_use]
    pub const fn geometry(&self) -> &GeometryStore {
        &self.geometry
    }

    /// Mutable access to the geometry, for adding curves and surfaces.
    #[must_use]
    pub const fn geometry_mut(&mut self) -> &mut GeometryStore {
        &mut self.geometry
    }

    /// Intern a transform for use in placements.
    pub fn add_datum(&mut self, transform: Transform) -> DatumId {
        self.datums.insert(transform)
    }

    /// The node behind a shape's handle.
    #[must_use]
    pub fn node(&self, shape: &Shape) -> Option<&TShape> {
        self.nodes.get(shape.node())
    }

    /// The node behind a handle.
    #[must_use]
    pub fn node_by_id(&self, id: TShapeId) -> Option<&TShape> {
        self.nodes.get(id)
    }

    /// Mutable access to the node behind a shape's handle.
    ///
    /// For attaching geometry to an entity that already exists — a pcurve
    /// joining an edge to a face it has just come to bound. Structural change
    /// still goes through the builders; this reaches the node's *data*, which
    /// no invariant here constrains on its own.
    #[must_use]
    pub fn node_mut(&mut self, shape: &Shape) -> Option<&mut TShape> {
        self.nodes.get_mut(shape.node())
    }

    /// What kind of shape this is.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the handle does not
    /// resolve in this model.
    pub fn kind_of(&self, shape: &Shape) -> OgeomResult<ShapeType> {
        let Some(node) = self.node(shape) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        Ok(node.kind())
    }

    /// The tolerance a shape carries, if it carries one.
    ///
    /// # Errors
    ///
    /// As [`Model::kind_of`].
    pub fn tolerance_of(&self, shape: &Shape) -> OgeomResult<Option<Tolerance>> {
        let Some(node) = self.node(shape) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        Ok(node.data().tolerance())
    }

    /// Number of topology nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the model holds no topology.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Every topology node, with its handle, in arena order.
    ///
    /// For writing a document out. Traversal from a root shape reaches only
    /// what that root bounds; a document is everything in it.
    pub fn nodes(&self) -> impl Iterator<Item = (TShapeId, &TShape)> {
        self.nodes.iter()
    }

    /// Every node that has been given an identity, with it.
    pub fn identities(&self) -> impl Iterator<Item = (TShapeId, EntityId)> {
        self.nodes
            .iter()
            .filter_map(|(id, _)| self.identity.get(&id).map(|entity| (id, *entity)))
    }

    /// Add a vertex.
    pub fn add_vertex(&mut self, data: VertexData) -> Shape {
        Shape::of(
            self.nodes
                .insert(TShape::leaf(ShapeType::Vertex, NodeData::Vertex(data))),
        )
    }

    /// Add a vertex at `point` with the minimum tolerance.
    pub fn add_point(&mut self, point: Point) -> Shape {
        self.add_vertex(VertexData::new(point))
    }

    /// Add an edge bounded by the given vertices.
    ///
    /// The vertices are the edge's ends, in order. A closed edge — a full
    /// circle — names the same vertex twice rather than once, so that walking
    /// its boundary yields a start and an end as every other edge does.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a bound is
    /// not a vertex, or if there are more than two;
    /// [`OgeomError::Invariant`](ogeom_core::OgeomError::Invariant) if a vertex's
    /// tolerance is tighter than the edge's, breaking the containment rule.
    pub fn add_edge(&mut self, data: EdgeData, bounds: &[Shape]) -> OgeomResult<Shape> {
        if bounds.len() > 2 {
            ogeom_bail!(
                Construction,
                "an edge has at most two bounding vertices, got {}",
                bounds.len()
            );
        }
        self.check_children(ShapeType::Vertex, bounds)?;
        // A vertex caps an edge, so it must be at least as uncertain as the
        // edge is — otherwise the cap does not reliably sit on what it caps.
        for bound in bounds {
            self.widen(bound, data.tolerance)?;
        }
        let node = self.nodes.insert(TShape::new(
            ShapeType::Edge,
            NodeData::Edge(Box::new(data)),
            bounds.to_vec(),
        ));
        self.record_primitive(node, Role::SOLE);
        Ok(Shape::of(node))
    }

    /// Add a wire from a sequence of edges.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is
    /// not an edge, or the wire is empty.
    pub fn add_wire(&mut self, edges: &[Shape]) -> OgeomResult<Shape> {
        if edges.is_empty() {
            ogeom_bail!(Construction, "a wire needs at least one edge");
        }
        self.check_children(ShapeType::Edge, edges)?;
        Ok(Shape::of(self.nodes.insert(TShape::container(
            ShapeType::Wire,
            edges.to_vec(),
        ))))
    }

    /// Add a face bounded by the given wires.
    ///
    /// A face with no wires covers its surface's whole domain, and is recorded
    /// as naturally restricted.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a bound is
    /// not a wire; [`OgeomError::Invariant`](ogeom_core::OgeomError::Invariant) if the
    /// containment rule is broken.
    pub fn add_face(&mut self, mut data: FaceData, wires: &[Shape]) -> OgeomResult<Shape> {
        self.check_children(ShapeType::Wire, wires)?;
        if wires.is_empty() {
            data.natural_restriction = true;
        }
        // An edge borders a face, so it must be at least as uncertain as the
        // face. Walk the wires' edges and widen them where they are not.
        let face_tolerance = data.tolerance;
        for wire in wires {
            let edges = self.children_of(wire)?;
            for edge in &edges {
                self.widen(edge, face_tolerance)?;
            }
        }
        let node = self.nodes.insert(TShape::new(
            ShapeType::Face,
            NodeData::Face(Box::new(data)),
            wires.to_vec(),
        ));
        self.record_primitive(node, Role::SOLE);
        Ok(Shape::of(node))
    }

    /// Add a shell from a set of faces.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is
    /// not a face, or the shell is empty.
    pub fn add_shell(&mut self, faces: &[Shape]) -> OgeomResult<Shape> {
        if faces.is_empty() {
            ogeom_bail!(Construction, "a shell needs at least one face");
        }
        self.check_children(ShapeType::Face, faces)?;
        Ok(Shape::of(self.nodes.insert(TShape::container(
            ShapeType::Shell,
            faces.to_vec(),
        ))))
    }

    /// Add a solid bounded by the given shells.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is
    /// not a shell, or the solid is empty.
    pub fn add_solid(&mut self, shells: &[Shape]) -> OgeomResult<Shape> {
        if shells.is_empty() {
            ogeom_bail!(Construction, "a solid needs at least one shell");
        }
        self.check_children(ShapeType::Shell, shells)?;
        Ok(Shape::of(self.nodes.insert(TShape::container(
            ShapeType::Solid,
            shells.to_vec(),
        ))))
    }

    /// Add a compsolid from solids sharing faces.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is
    /// not a solid, or it is empty.
    pub fn add_compsolid(&mut self, solids: &[Shape]) -> OgeomResult<Shape> {
        if solids.is_empty() {
            ogeom_bail!(Construction, "a compsolid needs at least one solid");
        }
        self.check_children(ShapeType::Solid, solids)?;
        Ok(Shape::of(self.nodes.insert(TShape::container(
            ShapeType::CompSolid,
            solids.to_vec(),
        ))))
    }

    /// Add a compound of arbitrary shapes.
    ///
    /// The one container with no type constraint — that is what a compound is
    /// for. It may be empty, since an empty result is a legitimate answer from
    /// a boolean and needs somewhere to live.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a child does not
    /// resolve in this model.
    pub fn add_compound(&mut self, shapes: &[Shape]) -> OgeomResult<Shape> {
        for shape in shapes {
            if self.node(shape).is_none() {
                ogeom_bail!(Dangling, "compound member is not in this model");
            }
        }
        Ok(Shape::of(self.nodes.insert(TShape::container(
            ShapeType::Compound,
            shapes.to_vec(),
        ))))
    }

    /// The direct children of a shape, with this shape's placement and
    /// orientation composed onto each.
    ///
    /// The single most important method on the model, and the reason traversal
    /// is correct by default rather than by discipline: a child's placement in
    /// the world is its parent's composed with its own, and its orientation is
    /// its parent's composed with its own. Returning raw children would leave
    /// every caller to remember both, and the failure is silent — face normals
    /// that flip inconsistently, sub-shapes drawn at the origin.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape does not
    /// resolve in this model.
    pub fn children_of(&self, shape: &Shape) -> OgeomResult<Vec<Shape>> {
        let Some(node) = self.node(shape) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        Ok(node
            .children()
            .iter()
            .map(|child| child.moved(shape.location()).composed(shape.orientation()))
            .collect())
    }

    /// A shape's children in *traversal* order.
    ///
    /// The same shapes as [`Model::children_of`], but with the list reversed
    /// when the parent is reversed. Order carries meaning for a wire — its
    /// edges run head to tail — and reversing a wire has to reverse the walk as
    /// well as each edge, or consecutive edges stop sharing a vertex and the
    /// boundary comes apart. For a shell or a solid the order means nothing and
    /// the reversal is invisible.
    ///
    /// [`Model::children_of`] stays the raw accessor: it returns what is
    /// stored, which is what a rebuild or a comparison wants.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape does not
    /// resolve in this model.
    pub fn ordered_children_of(&self, shape: &Shape) -> OgeomResult<Vec<Shape>> {
        let mut children = self.children_of(shape)?;
        if shape.orientation() == Orientation::Reversed {
            children.reverse();
        }
        Ok(children)
    }

    /// Widen a shape's tolerance, and every sub-shape's with it.
    ///
    /// The cascade is the point. The containment rule is transitive: a face's
    /// edges must be no tighter than the face, *and* those edges' vertices no
    /// tighter than the edges. Widening only one level leaves the rule broken
    /// two levels down, where nothing will notice until a containment test
    /// quietly answers about geometry that does not meet.
    ///
    /// Tolerances only ever grow, so this is the sanctioned repair: raise what
    /// bounds, never lower what is bounded.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape, or
    /// anything below it, does not resolve in this model.
    pub fn widen(&mut self, shape: &Shape, to: Tolerance) -> OgeomResult<()> {
        let mut affected = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut stack = vec![shape.node()];
        while let Some(id) = stack.pop() {
            if !seen.insert(id) {
                continue;
            }
            let Some(node) = self.nodes.get(id) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            affected.push(id);
            stack.extend(node.children().iter().map(Shape::node));
        }
        for id in affected {
            if let Some(node) = self.nodes.get_mut(id) {
                node.data_mut().widen(to);
            }
        }
        Ok(())
    }

    /// Check that every child resolves and is of the expected type.
    fn check_children(&self, expected: ShapeType, children: &[Shape]) -> OgeomResult<()> {
        for child in children {
            let Some(node) = self.node(child) else {
                ogeom_bail!(Dangling, "child refers to a node not in this model");
            };
            if node.kind() != expected {
                ogeom_bail!(
                    Construction,
                    "expected a {expected:?} child, got a {:?}",
                    node.kind()
                );
            }
        }
        Ok(())
    }

    /// Verify the containment rule across a whole shape tree.
    ///
    /// `docs/DATA_MODEL.md` §5. Walks parent to child and checks that whatever
    /// bounds is no tighter than what it bounds.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Invariant`](ogeom_core::OgeomError::Invariant) at the first
    /// violation, naming the two shape types involved.
    pub fn check_tolerances(&self, root: &Shape) -> OgeomResult<()> {
        let Some(node) = self.node(root) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        let own = node.data().tolerance();
        for child in self.children_of(root)? {
            // A boundary is *contained by* what it bounds, so the child — the
            // boundary — must be the looser of the two.
            if let (Some(parent), Some(child_tolerance)) = (own, self.tolerance_of(&child)?)
                && child_tolerance < parent
            {
                ogeom_bail!(
                    Invariant,
                    "a {:?} at tolerance {} bounds a {:?} at {}, which is tighter",
                    self.kind_of(&child)?,
                    child_tolerance.get(),
                    node.kind(),
                    parent.get()
                );
            }
            self.check_tolerances(&child)?;
        }
        Ok(())
    }

    /// Whether two shapes coincide in position, comparing composed transforms.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`](crate::Location::composed).
    pub fn same_position(&self, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<bool> {
        a.is_same_position(b, &self.datums, tol)
    }

    /// A shape placed by an additional transform.
    ///
    /// Interns the transform and composes it onto the shape's placement, so the
    /// underlying node — and all its geometry — is shared rather than copied.
    /// Placing ten thousand instances of a part costs ten thousand short chains
    /// and one copy of the geometry.
    pub fn placed(&mut self, shape: &Shape, transform: Transform) -> Shape {
        let datum = self.add_datum(transform);
        shape.moved(&Location::of(datum))
    }
}

/// What an absorb produced: the transplanted roots, bound to the model that
/// absorbed them, and where every absorbed identity ended up.
#[derive(Debug)]
pub struct Absorbed {
    /// The roots handed in, shifted and bound to the absorbing model.
    pub shapes: Vec<Shape>,
    /// Source-document identity → identity in the absorbing model, for every
    /// entity the parts carried. A caller holding references recorded against
    /// the source document resolves them through this.
    pub entities: HashMap<EntityId, EntityId>,
}

/// A model's contents, laid out the way a file holds them.
///
/// Handed to [`Model::from_parts`]. Every list is in arena order, and a node's
/// children name other nodes by their position in `nodes` — so the order is
/// load-bearing rather than incidental, and a reader has to preserve it.
#[derive(Debug, Default)]
pub struct ModelParts {
    /// The topology nodes.
    pub nodes: Vec<TShape>,
    /// The placement datums.
    pub datums: Vec<crate::location::Datum>,
    /// The geometry, already assembled.
    pub geometry: GeometryStore,
    /// Every entity's provenance, in the order identities were issued: the
    /// first entry is `EntityId(1)`.
    pub provenance: Vec<Provenance>,
    /// Which identity each node carries.
    pub identity: Vec<(TShapeId, EntityId)>,
    /// The operation the document was left in.
    pub current_op: OpId,
    /// The unit scale the document was authored at.
    pub tolerances: Tolerances,
}

/// Which sub-shapes a traversal should yield.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    /// Every shape of one type.
    OfType(ShapeType),
    /// Every shape, at every level.
    All,
}

/// Walk a shape's tree, composing placement and orientation on descent.
///
/// Yields each matching sub-shape with its *effective* placement and
/// orientation — the composition of everything from the root down. That is the
/// only form in which a sub-shape means anything outside the tree it came from.
///
/// Sub-shapes reached by more than one route — an edge shared by two faces —
/// are yielded once per route, since each occurrence has its own orientation
/// and that is usually the point. Deduplicate with
/// [`SameKey`](crate::SameKey) when it is not.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if any handle fails to
/// resolve in `model`.
pub fn explore(model: &Model, root: &Shape, filter: Filter) -> OgeomResult<Vec<Shape>> {
    let mut out = Vec::new();
    let mut stack = vec![root.clone()];
    while let Some(shape) = stack.pop() {
        let matches = match filter {
            Filter::OfType(want) => model.kind_of(&shape)? == want,
            Filter::All => true,
        };
        if matches {
            out.push(shape.clone());
        }
        // `children_of` is what does the composing; nothing here re-derives it.
        for child in model.children_of(&shape)?.into_iter().rev() {
            stack.push(child);
        }
    }
    Ok(out)
}

/// Walk a shape's tree, yielding every sub-shape of `want` exactly once.
///
/// Deduplicated by [`Shape::is_same`] — node and placement, ignoring
/// orientation — which is what "the distinct edges of this solid" means.
///
/// # Errors
///
/// As [`explore`].
pub fn explore_unique(model: &Model, root: &Shape, want: ShapeType) -> OgeomResult<Vec<Shape>> {
    use std::collections::HashSet;

    use crate::shape::SameKey;

    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for shape in explore(model, root, Filter::OfType(want))? {
        if seen.insert(SameKey(shape.clone())) {
            out.push(shape);
        }
    }
    Ok(out)
}

/// Every shape in `root` that has `target` among its sub-shapes.
///
/// The inverse of traversal: "which faces meet at this edge?" is how a boolean
/// decides what a split affects, and how a fillet finds what it is blending.
///
/// # Errors
///
/// As [`explore`].
pub fn ancestors_of(
    model: &Model,
    root: &Shape,
    target: &Shape,
    want: ShapeType,
) -> OgeomResult<Vec<Shape>> {
    let mut out = Vec::new();
    for candidate in explore(model, root, Filter::OfType(want))? {
        if explore(model, &candidate, Filter::All)?
            .iter()
            .any(|s| s.is_same(target))
        {
            out.push(candidate);
        }
    }
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_geom::PlaneSurface;
    use ogeom_math::{Direction, Frame, Plane, Vector};

    const T: Tolerances = Tolerances::millimetres();

    /// A single square face: four vertices, four edges, one wire, one face.
    fn square(model: &mut Model) -> Shape {
        let corners = [
            model.add_point(Point::new(0.0, 0.0, 0.0)),
            model.add_point(Point::new(1.0, 0.0, 0.0)),
            model.add_point(Point::new(1.0, 1.0, 0.0)),
            model.add_point(Point::new(0.0, 1.0, 0.0)),
        ];
        let mut edges = Vec::new();
        for i in 0..4 {
            let bounds = [corners[i].clone(), corners[(i + 1) % 4].clone()];
            edges.push(model.add_edge(EdgeData::new(), &bounds).unwrap());
        }
        let wire = model.add_wire(&edges).unwrap();
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        model
            .add_face(FaceData::new(surface, Location::identity()), &[wire])
            .unwrap()
    }

    /// Parts describing a two-vertex edge on a line, placed by one datum,
    /// with a primitive-and-derived provenance chain: every kind of handle a
    /// shift must move, in miniature, exactly as a reader would rebuild them.
    fn edge_parts() -> (ModelParts, Vec<Shape>) {
        use ogeom_core::Role;

        use crate::entity::CurveId;
        use crate::location::DatumId;

        let mut geometry = GeometryStore::new();
        geometry.add_curve(
            ogeom_geom::LineCurve::new(ogeom_math::Axis {
                location: Point::new(0.0, 0.0, 0.0),
                direction: Direction::X,
            })
            .into(),
        );
        let ends = [
            Shape::of(TShapeId::from_parts(0, 0)),
            Shape::of(TShapeId::from_parts(1, 0)).reversed(),
        ];
        let edge = TShape::new(
            ShapeType::Edge,
            NodeData::Edge(Box::new(EdgeData::on_curve(
                CurveId::from_parts(0, 0),
                Location::of(DatumId::from_parts(0, 0)),
                (0.0, 2.0),
            ))),
            ends.to_vec(),
        );
        let one = EntityId::from_raw(1).unwrap();
        let two = EntityId::from_raw(2).unwrap();
        let parts = ModelParts {
            nodes: vec![
                TShape::leaf(
                    ShapeType::Vertex,
                    NodeData::Vertex(VertexData::new(Point::new(0.0, 0.0, 0.0))),
                ),
                TShape::leaf(
                    ShapeType::Vertex,
                    NodeData::Vertex(VertexData::new(Point::new(2.0, 0.0, 0.0))),
                ),
                edge,
            ],
            datums: vec![Transform::translation(Vector::new(0.0, 0.0, 1.0))],
            geometry,
            provenance: vec![
                Provenance::Primitive {
                    op: OpId(1),
                    role: Role::SOLE,
                },
                Provenance::Derived {
                    op: OpId(1),
                    from: [one].into_iter().collect(),
                    role: Role::SOLE,
                },
            ],
            identity: vec![(TShapeId::from_parts(2, 0), two)],
            current_op: OpId(1),
            tolerances: T,
        };
        (parts, vec![Shape::of(TShapeId::from_parts(2, 0))])
    }

    #[test]
    fn absorbing_parts_into_a_live_model_offsets_every_handle() {
        let mut model = Model::new();
        let face = square(&mut model);

        let (parts, roots) = edge_parts();
        let absorbed = model.absorb(parts, &roots).unwrap();
        assert_eq!(absorbed.shapes.len(), 1);
        let edge = &absorbed.shapes[0];

        // The absorbed root resolves here and answers as an edge.
        assert_eq!(model.kind_of(edge).unwrap(), ShapeType::Edge);
        let vertices = explore(&model, edge, Filter::OfType(ShapeType::Vertex)).unwrap();
        assert_eq!(vertices.len(), 2);
        let points: Vec<Point> = vertices
            .iter()
            .map(|v| model.node(v).unwrap().data().as_vertex().unwrap().point)
            .collect();
        assert!(points.contains(&Point::new(2.0, 0.0, 0.0)), "{points:?}");

        // Its curve and datum handles were shifted onto this model's arenas.
        let node = model.node(edge).unwrap();
        let repr = &node.data().as_edge().unwrap().representations[0];
        assert!(
            model.geometry().holds(repr),
            "the absorbed edge's curve did not land"
        );
        let EdgeRepr::Curve3d { location, .. } = repr else {
            panic!("the representation changed kind in the shift");
        };
        assert!(
            location.composed(model.datums()).is_ok(),
            "the absorbed edge's datum did not land"
        );

        // What was here before is untouched.
        assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
    }

    #[test]
    fn absorbed_identities_keep_their_provenance_under_new_ids() {
        let mut model = Model::new();
        square(&mut model);
        let issued_before = model.provenance().len() as u64;
        assert!(
            issued_before > 0,
            "the square should have minted identities"
        );

        let (parts, roots) = edge_parts();
        let absorbed = model.absorb(parts, &roots).unwrap();

        // The remap table is exactly old-plus-offset.
        let old = EntityId::from_raw(2).unwrap();
        let new = absorbed.entities[&old];
        assert_eq!(new.get(), 2 + issued_before);
        assert_eq!(model.identity_of(&absorbed.shapes[0]), Some(new));

        // The derived entry still points at its shifted source, and walking
        // back lands on the shifted primitive.
        let entry = model.provenance().get(new).unwrap();
        let source = EntityId::from_raw(1 + issued_before).unwrap();
        assert_eq!(entry.inputs(), &[source]);
        assert_eq!(model.provenance().roots(new), vec![source]);
    }

    #[test]
    fn absorbing_into_an_empty_model_matches_from_parts() {
        let (parts, roots) = edge_parts();
        let restored = Model::from_parts(parts).unwrap();
        let bound = restored.bind(&roots[0]).unwrap();

        let (parts, roots) = edge_parts();
        let mut empty = Model::new();
        let absorbed = empty.absorb(parts, &roots).unwrap();
        let shape = &absorbed.shapes[0];

        // Zero offsets: the handles come out where the file put them, and the
        // identities are the file's own.
        assert_eq!(shape.node().index(), bound.node().index());
        assert_eq!(shape.node().generation(), bound.node().generation());
        assert_eq!(restored.identity_of(&bound), empty.identity_of(shape));
        assert_eq!(restored.provenance().len(), empty.provenance().len());
    }

    #[test]
    fn parts_with_scoped_or_generation_bearing_keys_are_refused() {
        let mut model = Model::new();

        let (mut parts, roots) = edge_parts();
        let child = parts.nodes[2].children()[0].clone();
        parts.nodes[2].children_mut()[0] = Shape::new(
            child.node().with_scope(7),
            Location::identity(),
            Orientation::Forward,
        );
        assert!(
            model.absorb(parts, &roots).is_err(),
            "a scoped child key should be refused"
        );

        let (mut parts, roots) = edge_parts();
        parts.identity[0].0 = TShapeId::from_parts(2, 1);
        assert!(
            model.absorb(parts, &roots).is_err(),
            "a recycled-generation key should be refused"
        );
    }

    #[test]
    fn parts_in_other_units_are_refused() {
        let mut model = Model::new();
        square(&mut model);
        let issued_before = model.provenance().len();

        let (mut parts, roots) = edge_parts();
        parts.tolerances = Tolerances::metres();
        assert!(model.absorb(parts, &roots).is_err());
        assert_eq!(
            model.provenance().len(),
            issued_before,
            "a refused absorb should leave the model alone"
        );
    }

    #[test]
    fn absorb_leaves_the_current_operation_alone() {
        let mut model = Model::new();
        model.begin_operation();
        model.begin_operation();
        let op = model.begin_operation();

        let (parts, roots) = edge_parts();
        model.absorb(parts, &roots).unwrap();
        assert_eq!(model.current_operation(), op);
    }

    #[test]
    fn an_absorbed_root_that_names_a_missing_node_dangles() {
        let mut model = Model::new();
        let (parts, _) = edge_parts();
        let stray = vec![Shape::of(TShapeId::from_parts(99, 0))];
        assert!(model.absorb(parts, &stray).is_err());
    }

    #[test]
    fn absorbing_empty_parts_is_a_no_op() {
        let mut model = Model::new();
        square(&mut model);
        let issued_before = model.provenance().len();

        let parts = ModelParts {
            tolerances: T,
            ..ModelParts::default()
        };
        let absorbed = model.absorb(parts, &[]).unwrap();
        assert!(absorbed.shapes.is_empty());
        assert!(absorbed.entities.is_empty());
        assert_eq!(model.provenance().len(), issued_before);
    }

    #[test]
    fn reversing_a_wire_reverses_the_walk_as_well_as_each_edge() {
        // Reversing each edge without reversing the order breaks the chain:
        // edge 1 would end where it used to start while edge 2 still starts
        // where it used to, so consecutive edges stop meeting and a face built
        // on the wire comes apart along its boundary.
        let mut model = Model::new();
        let face = square(&mut model);
        let wire = model.children_of(&face).unwrap()[0].clone();

        let forward = model.ordered_children_of(&wire).unwrap();
        let backward = model.ordered_children_of(&wire.reversed()).unwrap();

        assert_eq!(forward.len(), 4);
        assert_eq!(backward.len(), 4);
        for (i, edge) in backward.iter().enumerate() {
            let partner = &forward[3 - i];
            assert!(edge.is_same(partner), "the order did not reverse");
            assert_eq!(
                edge.orientation(),
                Orientation::Reversed.compose(partner.orientation()),
                "each edge should also flip"
            );
        }

        // The raw accessor keeps the stored order, which is what a rebuild
        // wants and what a traversal must not use.
        let raw = model.children_of(&wire.reversed()).unwrap();
        assert!(raw[0].is_same(&forward[0]));
    }

    #[test]
    fn a_built_tree_has_the_expected_shape() {
        let mut model = Model::new();
        let face = square(&mut model);

        assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
        assert_eq!(
            explore_unique(&model, &face, ShapeType::Wire)
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            explore_unique(&model, &face, ShapeType::Edge)
                .unwrap()
                .len(),
            4
        );
        assert_eq!(
            explore_unique(&model, &face, ShapeType::Vertex)
                .unwrap()
                .len(),
            4
        );
        // Four edges of two vertices each, but only four distinct vertices:
        // consecutive edges share them.
        assert_eq!(
            explore(&model, &face, Filter::OfType(ShapeType::Vertex))
                .unwrap()
                .len(),
            8
        );
    }

    #[test]
    fn children_are_returned_with_the_parents_placement_composed() {
        // The invariant traversal exists to guarantee. A vertex reported
        // without its parent's placement is a vertex at the wrong point, and
        // nothing about the value says so.
        let mut model = Model::new();
        let face = square(&mut model);
        let moved = model.placed(&face, Transform::translation(Vector::new(10.0, 0.0, 0.0)));

        let vertices = explore_unique(&model, &moved, ShapeType::Vertex).unwrap();
        assert_eq!(vertices.len(), 4);
        for v in &vertices {
            let node = model.node(v).unwrap();
            let local = node.data().as_vertex().unwrap().point;
            let world = v.transform(model.datums()).unwrap().apply(local);
            assert!(world.x >= 10.0 - 1e-12, "vertex at {world:?} was not moved");
        }
    }

    #[test]
    fn children_are_returned_with_the_parents_orientation_composed() {
        let mut model = Model::new();
        let face = square(&mut model);
        let reversed = face.reversed();

        let forward_edges = model
            .children_of(&model.children_of(&face).unwrap()[0])
            .unwrap();
        let reversed_edges = model
            .children_of(&model.children_of(&reversed).unwrap()[0])
            .unwrap();

        for (a, b) in forward_edges.iter().zip(&reversed_edges) {
            assert_eq!(
                b.orientation(),
                a.orientation().reversed(),
                "reversing a face must reverse what its edges present, \
                 without touching a single stored child"
            );
        }
    }

    #[test]
    fn reversing_a_shape_touches_no_stored_child() {
        // The point of composing on descent: the reversal lives entirely in the
        // handle, so a shared sub-tree is not disturbed for other users of it.
        let mut model = Model::new();
        let face = square(&mut model);
        let before = model.node(&face).unwrap().clone();
        let _ = face.reversed();
        assert_eq!(model.node(&face).unwrap(), &before);
    }

    #[test]
    fn placement_composes_through_nesting() {
        let mut model = Model::new();
        let vertex = model.add_point(Point::new(1.0, 0.0, 0.0));
        let edge = model
            .add_edge(EdgeData::new(), &[vertex.clone(), vertex.clone()])
            .unwrap();
        let moved_edge = model.placed(&edge, Transform::translation(Vector::new(10.0, 0.0, 0.0)));
        let compound = model.add_compound(&[moved_edge]).unwrap();
        let moved_compound = model.placed(
            &compound,
            Transform::translation(Vector::new(100.0, 0.0, 0.0)),
        );

        let found = explore_unique(&model, &moved_compound, ShapeType::Vertex).unwrap();
        assert_eq!(found.len(), 1);
        let local = model
            .node(&found[0])
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        let world = found[0].transform(model.datums()).unwrap().apply(local);
        assert!(
            world.is_equal(Point::new(111.0, 0.0, 0.0), T),
            "expected 1 + 10 + 100, got {world:?}"
        );
    }

    #[test]
    fn the_builder_refuses_children_of_the_wrong_type() {
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);
        let edge = model
            .add_edge(EdgeData::new(), std::slice::from_ref(&vertex))
            .unwrap();

        assert!(
            model.add_wire(std::slice::from_ref(&vertex)).is_err(),
            "a wire holds edges"
        );
        assert!(
            model.add_shell(std::slice::from_ref(&edge)).is_err(),
            "a shell holds faces"
        );
        assert!(
            model.add_solid(std::slice::from_ref(&edge)).is_err(),
            "a solid holds shells"
        );
        assert!(model.add_wire(&[edge]).is_ok());

        // A compound is the exception, and deliberately so.
        assert!(model.add_compound(&[vertex]).is_ok());
    }

    #[test]
    fn empty_containers_are_refused_except_a_compound() {
        let mut model = Model::new();
        assert!(model.add_wire(&[]).is_err());
        assert!(model.add_shell(&[]).is_err());
        assert!(model.add_solid(&[]).is_err());
        assert!(model.add_compsolid(&[]).is_err());
        // An empty result is a legitimate answer from a boolean and needs
        // somewhere to live.
        assert!(model.add_compound(&[]).is_ok());
    }

    #[test]
    fn an_edge_takes_at_most_two_vertices() {
        let mut model = Model::new();
        let v = model.add_point(Point::ORIGIN);
        assert!(model.add_edge(EdgeData::new(), &[]).is_ok(), "unbounded");
        assert!(
            model
                .add_edge(EdgeData::new(), std::slice::from_ref(&v))
                .is_ok()
        );
        assert!(
            model
                .add_edge(EdgeData::new(), &[v.clone(), v.clone()])
                .is_ok()
        );
        assert!(
            model
                .add_edge(EdgeData::new(), &[v.clone(), v.clone(), v])
                .is_err(),
            "three ends is not an edge"
        );
    }

    #[test]
    fn building_enforces_the_containment_rule_upward() {
        // A coarse edge must not be capped by a finer vertex. Rather than
        // refusing, the builder widens the vertex — tolerances only ever grow,
        // so the repair goes upward.
        let mut model = Model::new();
        let vertex = model.add_point(Point::ORIGIN);
        assert_eq!(model.tolerance_of(&vertex).unwrap(), Some(Tolerance::MIN));

        let mut edge_data = EdgeData::new();
        edge_data.widen(Tolerance::new(1e-3).unwrap());
        let edge = model
            .add_edge(edge_data, std::slice::from_ref(&vertex))
            .unwrap();

        assert_eq!(
            model.tolerance_of(&vertex).unwrap(),
            Some(Tolerance::new(1e-3).unwrap()),
            "the vertex was widened to contain its edge"
        );
        assert!(model.check_tolerances(&edge).is_ok());
    }

    #[test]
    fn a_face_widens_the_edges_it_borders() {
        let mut model = Model::new();
        let a = model.add_point(Point::ORIGIN);
        let b = model.add_point(Point::new(1.0, 0.0, 0.0));
        let edge = model.add_edge(EdgeData::new(), &[a.clone(), b]).unwrap();
        let wire = model.add_wire(std::slice::from_ref(&edge)).unwrap();

        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let mut face_data = FaceData::new(surface, Location::identity());
        face_data.widen(Tolerance::new(1e-2).unwrap());
        let face = model.add_face(face_data, &[wire]).unwrap();

        assert_eq!(
            model.tolerance_of(&edge).unwrap(),
            Some(Tolerance::new(1e-2).unwrap())
        );
        // And the cascade reached the vertices under those edges. Stopping one
        // level down would leave the rule broken where nothing looks.
        assert_eq!(
            model.tolerance_of(&a).unwrap(),
            Some(Tolerance::new(1e-2).unwrap()),
            "widening a face must reach its edges' vertices, not just its edges"
        );
        assert!(model.check_tolerances(&face).is_ok());
    }

    #[test]
    fn check_tolerances_catches_a_violation_the_builder_would_never_make() {
        // The builder maintains the rule, so a violation has to be assembled
        // around it — which is exactly what happens when topology arrives from
        // a file. The check has to stand on its own, or imported geometry sails
        // past it.
        let mut model = Model::new();
        let vertex = Shape::of(model.nodes.insert(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        )));

        let mut edge_data = EdgeData::new();
        edge_data.widen(Tolerance::new(1e-1).unwrap());
        let edge = Shape::of(model.nodes.insert(TShape::new(
            ShapeType::Edge,
            NodeData::Edge(Box::new(edge_data)),
            vec![vertex.clone()],
        )));

        let err = model.check_tolerances(&edge).unwrap_err();
        assert!(
            err.to_string().contains("tighter"),
            "unexpected message: {err}"
        );

        // And the sanctioned repair fixes it, cascading to the vertex.
        model.widen(&edge, Tolerance::new(1e-1).unwrap()).unwrap();
        assert!(model.check_tolerances(&edge).is_ok());
        assert_eq!(
            model.tolerance_of(&vertex).unwrap(),
            Some(Tolerance::new(1e-1).unwrap())
        );
    }

    #[test]
    fn a_face_with_no_wires_is_naturally_restricted() {
        let mut model = Model::new();
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let face = model
            .add_face(FaceData::new(surface, Location::identity()), &[])
            .unwrap();
        assert!(
            model
                .node(&face)
                .unwrap()
                .data()
                .as_face()
                .unwrap()
                .natural_restriction,
            "an untrimmed face needs no point-in-face test at all"
        );
    }

    #[test]
    fn a_shared_sub_shape_is_yielded_once_per_route_and_deduplicated_on_request() {
        // Two faces meeting at an edge. Each occurrence carries its own
        // orientation, which is usually the point; asking for distinct edges is
        // a separate question.
        let mut model = Model::new();
        let a = model.add_point(Point::ORIGIN);
        let b = model.add_point(Point::new(1.0, 0.0, 0.0));
        let shared = model.add_edge(EdgeData::new(), &[a, b]).unwrap();

        let wire_one = model.add_wire(std::slice::from_ref(&shared)).unwrap();
        let wire_two = model.add_wire(&[shared.reversed()]).unwrap();
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let face_one = model
            .add_face(FaceData::new(surface, Location::identity()), &[wire_one])
            .unwrap();
        let face_two = model
            .add_face(FaceData::new(surface, Location::identity()), &[wire_two])
            .unwrap();
        let shell = model.add_shell(&[face_one, face_two]).unwrap();

        let all = explore(&model, &shell, Filter::OfType(ShapeType::Edge)).unwrap();
        assert_eq!(all.len(), 2, "one occurrence per route");
        assert_ne!(all[0].orientation(), all[1].orientation());

        let distinct = explore_unique(&model, &shell, ShapeType::Edge).unwrap();
        assert_eq!(distinct.len(), 1, "one edge, seen from two sides");
    }

    #[test]
    fn ancestors_answers_which_faces_meet_at_an_edge() {
        let mut model = Model::new();
        let a = model.add_point(Point::ORIGIN);
        let b = model.add_point(Point::new(1.0, 0.0, 0.0));
        let shared = model.add_edge(EdgeData::new(), &[a, b]).unwrap();
        let isolated = model.add_point(Point::new(5.0, 5.0, 5.0));
        let lone = model.add_edge(EdgeData::new(), &[isolated]).unwrap();

        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let mut faces = Vec::new();
        for _ in 0..2 {
            let wire = model.add_wire(std::slice::from_ref(&shared)).unwrap();
            faces.push(
                model
                    .add_face(FaceData::new(surface, Location::identity()), &[wire])
                    .unwrap(),
            );
        }
        let third_wire = model.add_wire(std::slice::from_ref(&lone)).unwrap();
        faces.push(
            model
                .add_face(FaceData::new(surface, Location::identity()), &[third_wire])
                .unwrap(),
        );
        let shell = model.add_shell(&faces).unwrap();

        let meeting = ancestors_of(&model, &shell, &shared, ShapeType::Face).unwrap();
        assert_eq!(meeting.len(), 2, "two faces meet at the shared edge");
        let alone = ancestors_of(&model, &shell, &lone, ShapeType::Face).unwrap();
        assert_eq!(alone.len(), 1);
    }

    #[test]
    fn handles_from_another_model_are_reported_rather_than_resolved() {
        let mut model = Model::new();
        let mut other = Model::new();
        // Past the end of `model`, so it genuinely fails to resolve.
        let mut foreign = other.add_point(Point::ORIGIN);
        for _ in 0..5 {
            foreign = other.add_point(Point::ORIGIN);
        }
        assert!(model.kind_of(&foreign).is_err());
        assert!(model.children_of(&foreign).is_err());
        assert!(model.add_wire(std::slice::from_ref(&foreign)).is_err());
        assert!(model.add_compound(&[foreign]).is_err());
    }

    #[test]
    fn placing_a_shape_shares_its_geometry_rather_than_copying_it() {
        // Ten thousand fasteners cost ten thousand short chains and one copy of
        // the geometry. That is the whole reason placement is a chain.
        let mut model = Model::new();
        let face = square(&mut model);
        let before = model.node_count();

        let mut instances = Vec::new();
        for i in 0..100 {
            instances.push(model.placed(
                &face,
                Transform::translation(Vector::new(f64::from(i), 0.0, 0.0)),
            ));
        }
        assert_eq!(model.node_count(), before, "no topology was duplicated");
        assert!(instances.iter().all(|s| s.is_partner(&face)));
        assert!(instances.iter().all(|s| !s.is_same(&face)));

        // And they are all in different places.
        let a = instances[0].transform(model.datums()).unwrap();
        let b = instances[99].transform(model.datums()).unwrap();
        assert!(!a.is_equal(&b, T));
    }

    #[test]
    fn an_empty_model_reports_itself_as_empty() {
        let model = Model::new();
        assert!(model.is_empty());
        assert_eq!(model.node_count(), 0);
        assert_eq!(model.geometry().counts(), (0, 0, 0));
        assert!(model.datums().is_empty());
    }

    #[test]
    fn a_vertex_has_no_children_and_traversal_stops_there() {
        let mut model = Model::new();
        let v = model.add_point(Point::new(1.0, 2.0, 3.0));
        assert!(model.children_of(&v).unwrap().is_empty());
        assert_eq!(explore(&model, &v, Filter::All).unwrap().len(), 1);
        assert!(model.check_tolerances(&v).is_ok());
    }

    #[test]
    fn shapes_at_different_places_are_not_the_same_position() {
        let mut model = Model::new();
        let v = model.add_point(Point::ORIGIN);
        let moved = model.placed(&v, Transform::translation(Vector::X));
        assert!(!model.same_position(&v, &moved, T).unwrap());
        assert!(model.same_position(&v, &v.clone(), T).unwrap());

        // The same displacement reached twice is the same position, even
        // through two different datums.
        let again = model.placed(&v, Transform::translation(Vector::X));
        assert!(!moved.is_same(&again), "structurally different chains");
        assert!(model.same_position(&moved, &again, T).unwrap());
    }

    #[test]
    fn a_direction_is_needed_to_build_a_non_trivial_plane() {
        // Guards the test helper itself: a face built on a degenerate plane
        // would make every other assertion here meaningless.
        let mut model = Model::new();
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into());
        assert!(model.geometry().surface(surface).is_some());
    }
}
