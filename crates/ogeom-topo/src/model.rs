//! The model: the arenas a shape's handles refer into, and the builder that is
//! the only way to mutate them.
//!
//! A [`Model`] owns the topology nodes, the placement datums and the geometry.
//! A [`Shape`] is meaningless without one: its handles index into these arenas
//! and nothing else (`docs/DATA_MODEL.md` §11).
//!
//! # One mutation path
//!
//! Every structural change goes through [`Model`]'s builder methods. That is
//! not ceremony: the invariants of `docs/DATA_MODEL.md` (a wire holds edges, a
//! face's tolerance does not exceed its edges', a node's kind matches its data)
//! are checkable in one place only if there is one place. Handing out
//! `&mut TShape` would scatter them across every caller, and the failures they
//! guard against are silent ones.

use ogeom_core::{
    Arena, EntityId, Key, OgeomResult, OpId, Provenance, ProvenanceTable, Role, Tolerance,
    Tolerances, ogeom_bail,
};
use ogeom_math::{Aabb, Point, Transform, TransformKind};

use smallvec::SmallVec;

use crate::entity::{
    CurveId, EdgeData, EdgeRepr, FaceData, NodeData, PCurveId, SurfaceId, TriangulationId,
    VertexData,
};
use crate::kept::FaceBoxes;
use crate::location::{DatumId, DatumStore, Location};
use crate::shape::{Orientation, Shape, ShapeType, TShape, TShapeId};

pub use crate::entity::GeometryStore;

/// A document: topology, placements, geometry, and where every entity came
/// from.
#[derive(Debug, Clone, Default)]
pub struct Model {
    nodes: Nodes,
    datums: DatumStore,
    geometry: GeometryStore,
    provenance: ProvenanceTable,
    current_op: OpId,
    tolerances: Tolerances,
    face_boxes: FaceBoxes,
    /// Whether any node is marked held ([`Model::note_held`]).
    any_held: bool,
    /// While a journal is open, each node [`Model::widen`] grew or
    /// [`Model::node_mut`] handed out, and its tolerance before: what
    /// [`Model::undo_widened`] puts back.
    widened: Option<Vec<(TShapeId, Tolerance)>>,
}

/// A topology node and what the model keeps beside it.
#[derive(Debug, Clone)]
struct Node {
    shape: TShape,
    /// The identity the node carries.
    identity: Option<EntityId>,
    /// The node's slot among the kept face boxes, if it is a face.
    face_box: Option<u32>,
    /// The node indices of the edges, wires and faces that hold it: the way
    /// up from a vertex, an edge or a wire to the faces it bounds.
    held_by: SmallVec<[u32; 2]>,
    /// Whether an operation passed the node through into its result from
    /// the shapes it was given, which those shapes still hold. What lies
    /// below a held node is held with it.
    held: bool,
}

/// The topology nodes, looked up by the handles shapes carry.
#[derive(Debug, Clone, Default)]
struct Nodes {
    arena: Arena<Node>,
}

/// A node's key in the arena, from its handle.
#[inline]
const fn arena_key(id: TShapeId) -> Key<Node> {
    Key::from_parts(id.index(), id.generation()).with_scope(id.scope())
}

/// A node's handle, from its key in the arena.
#[inline]
const fn handle(key: Key<Node>) -> TShapeId {
    Key::from_parts(key.index(), key.generation()).with_scope(key.scope())
}

impl Nodes {
    #[inline]
    fn entry(&self, id: TShapeId) -> Option<&Node> {
        self.arena.get(arena_key(id))
    }

    #[inline]
    fn entry_mut(&mut self, id: TShapeId) -> Option<&mut Node> {
        self.arena.get_mut(arena_key(id))
    }

    #[inline]
    fn get(&self, id: TShapeId) -> Option<&TShape> {
        self.entry(id).map(|node| &node.shape)
    }

    #[inline]
    fn get_mut(&mut self, id: TShapeId) -> Option<&mut TShape> {
        self.entry_mut(id).map(|node| &mut node.shape)
    }

    fn insert(&mut self, node: Node) -> TShapeId {
        handle(self.arena.insert(node))
    }

    fn iter(&self) -> impl Iterator<Item = (TShapeId, &Node)> {
        self.arena.iter().map(|(key, node)| (handle(key), node))
    }

    fn iter_mut(&mut self) -> impl Iterator<Item = (TShapeId, &mut Node)> {
        self.arena.iter_mut().map(|(key, node)| (handle(key), node))
    }

    /// The live handle at node index `index`.
    fn at(&self, index: u32) -> Option<TShapeId> {
        self.arena.key_at(index).map(handle)
    }

    const fn scope(&self) -> u32 {
        self.arena.scope()
    }

    const fn issued(&self, id: TShapeId) -> bool {
        id.scope() == self.arena.scope()
    }

    const fn len(&self) -> usize {
        self.arena.len()
    }

    const fn is_empty(&self) -> bool {
        self.arena.is_empty()
    }
}

/// What the roots of [`Model::retain_reachable`] reach.
#[derive(Default)]
struct Reach {
    /// Node indices.
    nodes: ogeom_core::FastSet<u32>,
    curves: ogeom_core::FastSet<CurveId>,
    pcurves: ogeom_core::FastSet<PCurveId>,
    surfaces: ogeom_core::FastSet<SurfaceId>,
    meshes: ogeom_core::FastSet<TriangulationId>,
    entities: ogeom_core::FastSet<EntityId>,
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
    /// a [`Tolerances`] argument (that is deliberate, since a caller may want
    /// to work coarser or finer than the document's own setting for one
    /// operation), but the document says what it was built at, so a
    /// mismatch is visible rather than assumed away.
    #[must_use]
    pub fn with_tolerances(tolerances: Tolerances) -> Self {
        Self {
            nodes: Nodes::default(),
            datums: DatumStore::new(),
            geometry: GeometryStore::new(),
            provenance: ProvenanceTable::new(),
            current_op: OpId(0),
            tolerances,
            face_boxes: FaceBoxes::default(),
            any_held: false,
            widened: None,
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
    /// written (every [`EntityId`] renumbered, every provenance record
    /// replaced by a fresh `Primitive` one), and every reference into it,
    /// which is the thing provenance exists to keep alive, is dead. Reading a
    /// file has to reproduce the document it describes, identities and all.
    ///
    /// This is not a hole in "the builder is the sole mutation path". Nothing
    /// here mutates an existing model; it assembles a new one, and it checks
    /// the structural invariants the builders check before handing it back,
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
    /// one in a model that already has things in it: the operation behind
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
    ///   rescaled; rescaling is a real feature with real decisions in it,
    ///   and silently absorbing metres into millimetres is a wrong model.
    /// - **Bound handles.** Parts whose handles already name an arena did not
    ///   come from a reader; absorbing them would alias whatever those
    ///   handles meant elsewhere. Serialization is the one road in.
    /// - **A model that reuses slots.** Absorbing appends by offset, which
    ///   is only sound while each of the target's arenas hands out fresh
    ///   slots in order. [`Model::retain_reachable`] leaves slots empty and
    ///   never hands them out again, so a retained model absorbs; an arena
    ///   waiting to refill a freed slot is refused, and the check is what
    ///   keeps that refill from becoming aliasing.
    ///
    /// The current operation is left alone: absorb mints no identities, it
    /// transplants a table, and the absorbed provenance keeps its source
    /// [`OpId`]s verbatim, meaningful in the source document's rebuild, kept
    /// because renumbering them would orphan the source's own references.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for
    /// the three refusals above;
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the parts
    /// do not describe themselves: a handle that does not resolve, a
    /// derivation from an identity no entry issued, a root naming a node that
    /// is not there.
    ///
    /// On an error past the up-front gates, the model's prior contents are
    /// untouched and still fully usable: an absorbed subgraph is
    /// self-contained (its shifted handles cannot reach below the append
    /// line), and no identity is committed until every check has passed. What
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
        if self.nodes.arena.next_index().is_none()
            || !self.datums.is_dense()
            || !self.geometry.appends()
        {
            ogeom_bail!(
                Construction,
                "absorb appends by offset, and this model's arenas are set to \
                 refill freed slots, which nothing in this crate does"
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
    /// Returns `(node, datum, entity)` offsets: where the parts landed.
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

        let node_offset = crate::entity::arena_len(&self.nodes.arena);
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
        self.sync_face_boxes();
        let appended: Vec<TShapeId> = nodes
            .into_iter()
            .map(|node| self.insert_unlinked(node))
            .collect();

        // Every handle in `parts` was rebuilt by a reader that had no arenas
        // to bind them to, so they name no arena at all and resolve nowhere.
        // Bind the appended subrange now that the arenas exist; what was here
        // before is already bound. A child may come after its parent in the
        // parts, so the links up from children are made once all are bound.
        self.bind_handles(node_offset);
        for id in appended {
            self.link_children(id);
        }
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
            if let Some(node) = self.nodes.entry_mut(node) {
                node.identity = Some(entity);
            }
        }
        Ok((node_offset, datum_offset, entity_offset))
    }

    /// Refuse parts whose handles are already bound to an arena or carry a
    /// non-zero generation: the state no reader produces, and the state an
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
    /// so, and it verifies the answer, so a file naming a node that is not
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
    /// that holds them exists. Binding checks them too: a chain naming a
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
            let node = &mut node.shape;
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
    /// still resolve through the full arena, so nothing is under-checked; an
    /// absorbed subgraph is self-contained and can only point at itself.
    fn check_restored(&self, identity: &[(TShapeId, EntityId)], from: u32) -> OgeomResult<()> {
        for (id, node) in self.nodes.iter().filter(|(id, _)| id.index() >= from) {
            let node = &node.shape;
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
    /// The counter is deterministic (the third operation in a rebuild is
    /// `OpId(3)` every time), which is what lets provenance survive a parameter
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
        self.nodes.entry(shape.node())?.identity
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

    /// Drop the provenance entries of `ids`, keeping the ids issued.
    ///
    /// What a document reader does with the entities its file left out: no
    /// shape the file carries names them and no kept entry derives from
    /// them, so they answer `None` as an id from another document does, and
    /// no later entity can take their numbers. Ids without an entry are
    /// passed over.
    pub fn forget_provenance(&mut self, ids: impl IntoIterator<Item = EntityId>) {
        for id in ids {
            self.provenance.forget(id);
        }
    }

    /// Drop everything `roots` do not reach, in place.
    ///
    /// Kept: every node below a root, the curves, pcurves, surfaces and
    /// triangulations those nodes name, and the provenance entries of the
    /// identities they carry with every entry those derive from. Every
    /// handle into what is kept resolves as before and means the same
    /// thing, and a lineage query on a kept shape answers as before. A
    /// handle into what is dropped fails to resolve, and no later entity
    /// takes its slot or its id, so it never answers about something else.
    /// Datums are kept whole: a placement a caller holds keeps resolving.
    ///
    /// What is dropped stops costing memory, and a clone costs what is
    /// kept. For a model that runs many operations and keeps one result:
    /// call it with that result once the model has grown well past what
    /// the result reaches.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a root
    /// does not resolve in this model; nothing is dropped then.
    pub fn retain_reachable(&mut self, roots: &[Shape]) -> OgeomResult<()> {
        if roots.iter().any(|root| self.node(root).is_none()) {
            ogeom_bail!(Dangling, "a root to retain is not in this model");
        }
        if self.widened.is_some() {
            ogeom_bail!(
                Construction,
                "an operation is noting the tolerances it widens; a model is \
                 retained between operations"
            );
        }
        let reach = self.reach(roots);

        self.nodes
            .arena
            .retain(|key, _| reach.nodes.contains(&key.index()));
        let mut kept_boxes: Vec<u32> = Vec::new();
        self.any_held = false;
        for (_, node) in self.nodes.iter_mut() {
            self.any_held |= node.held;
            node.held_by.retain(|above| reach.nodes.contains(above));
            if let Some(slot) = node.face_box.as_mut() {
                kept_boxes.push(*slot);
                *slot = u32::try_from(kept_boxes.len() - 1).unwrap_or(u32::MAX);
            }
        }
        self.face_boxes.keep_slots(&kept_boxes);
        self.geometry.retain(
            |id| reach.curves.contains(&id),
            |id| reach.pcurves.contains(&id),
            |id| reach.surfaces.contains(&id),
            |id| reach.meshes.contains(&id),
        );
        self.provenance.retain(|id, _| reach.entities.contains(&id));
        Ok(())
    }

    /// Everything `roots` reach: their nodes, the geometry those name, and
    /// the identities those carry with their ancestry.
    fn reach(&self, roots: &[Shape]) -> Reach {
        let mut reach = Reach::default();
        let mut queue: Vec<TShapeId> = Vec::new();
        for root in roots {
            if reach.nodes.insert(root.node().index()) {
                queue.push(root.node());
            }
        }
        let mut ancestry: Vec<EntityId> = Vec::new();
        while let Some(id) = queue.pop() {
            let Some(node) = self.nodes.entry(id) else {
                continue;
            };
            if let Some(entity) = node.identity
                && reach.entities.insert(entity)
            {
                ancestry.push(entity);
            }
            for child in node.shape.children() {
                if self.nodes.entry(child.node()).is_some()
                    && reach.nodes.insert(child.node().index())
                {
                    queue.push(child.node());
                }
            }
            match node.shape.data() {
                NodeData::Edge(edge) => {
                    for repr in &edge.representations {
                        match repr {
                            EdgeRepr::Curve3d { curve, .. } => {
                                reach.curves.insert(*curve);
                            }
                            EdgeRepr::PCurve { curve, surface, .. } => {
                                reach.pcurves.insert(*curve);
                                reach.surfaces.insert(*surface);
                            }
                            EdgeRepr::Seam {
                                forward,
                                reversed,
                                surface,
                                ..
                            } => {
                                reach.pcurves.insert(*forward);
                                reach.pcurves.insert(*reversed);
                                reach.surfaces.insert(*surface);
                            }
                            EdgeRepr::Polyline { .. } => {}
                            EdgeRepr::PolygonOnTriangulation { triangulation, .. } => {
                                reach.meshes.insert(*triangulation);
                            }
                        }
                    }
                }
                NodeData::Face(face) => {
                    reach.surfaces.insert(face.surface);
                    if let Some(mesh) = face.triangulation {
                        reach.meshes.insert(mesh);
                    }
                }
                NodeData::Vertex(_) | NodeData::Container => {}
            }
        }
        while let Some(entity) = ancestry.pop() {
            if let Some(entry) = self.provenance.get(entity) {
                for &from in entry.inputs() {
                    if reach.entities.insert(from) {
                        ancestry.push(from);
                    }
                }
            }
        }
        reach
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
    /// [`Shape`] cannot survive that: the reloaded document is a new set of
    /// arenas and [`Model::bind`] refuses a handle from another one, on
    /// purpose. An [`EntityId`] can, because it names *what the entity is*
    /// rather than where it sits (`docs/DATA_MODEL.md` §8), and that is the
    /// whole reason it exists.
    ///
    /// Returns the shape in its default placement and orientation. A caller
    /// that wants a particular occurrence explores from here.
    #[must_use]
    pub fn shape_of(&self, id: EntityId) -> Option<Shape> {
        self.nodes
            .iter()
            .find(|(_, node)| node.identity == Some(id))
            .map(|(node, _)| Shape::of(node))
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
        if let Some(node) = self.nodes.entry_mut(shape.node()) {
            node.identity = Some(id);
        }
        Ok(id)
    }

    /// Record a node's identity as it is created.
    fn record_primitive(&mut self, node: TShapeId, role: Role) {
        let id = self.provenance.primitive(self.current_op, role);
        if let Some(node) = self.nodes.entry_mut(node) {
            node.identity = Some(id);
        }
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
    ///
    /// A surface rewritten in place through it forgets every face's kept
    /// box ([`Model::face_bounds`]).
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
    /// For attaching geometry to an entity that already exists: a pcurve
    /// joining an edge to a face it has just come to bound. Structural change
    /// still goes through the builders; this reaches the node's *data*, which
    /// no invariant here constrains on its own.
    ///
    /// The node may be changed in any way, so the box kept for every face
    /// it is part of ([`Model::face_bounds`]) is forgotten. Every shape
    /// holding the node sees the change; an editor starting from a shape it
    /// was given copies what other shapes hold first ([`Model::unshare`]).
    #[must_use]
    pub fn node_mut(&mut self, shape: &Shape) -> Option<&mut TShape> {
        self.sync_face_boxes();
        if let Some(node) = self.nodes.get(shape.node()) {
            // The tolerance as it stands, which an edit through the node
            // may grow: what an attempt that fails puts back.
            if let Some(journal) = &mut self.widened
                && let Some(was) = node.data().tolerance()
            {
                journal.push((shape.node(), was));
            }
        }
        self.forget_boxes_above(shape.node());
        self.nodes.get_mut(shape.node())
    }

    /// Insert a node, with a slot for its box if it is a face and its links
    /// up from the children it holds.
    fn insert_node(&mut self, node: TShape) -> TShapeId {
        self.sync_face_boxes();
        let id = self.insert_unlinked(node);
        self.link_children(id);
        id
    }

    /// Insert a node, with a slot for its box if it is a face.
    fn insert_unlinked(&mut self, shape: TShape) -> TShapeId {
        let face_box = (shape.kind() == ShapeType::Face).then(|| self.face_boxes.add());
        self.nodes.insert(Node {
            shape,
            identity: None,
            face_box,
            held_by: SmallVec::new(),
            held: false,
        })
    }

    /// Record node `id` with each child it holds, where it is part of a
    /// face: the link [`Model::forget_boxes_above`] climbs.
    fn link_children(&mut self, id: TShapeId) {
        let Some(node) = self.nodes.get(id) else {
            return;
        };
        if !matches!(
            node.kind(),
            ShapeType::Edge | ShapeType::Wire | ShapeType::Face
        ) {
            return;
        }
        let children: SmallVec<[TShapeId; 8]> = node.children().iter().map(Shape::node).collect();
        for child in children {
            if let Some(below) = self.nodes.entry_mut(child)
                && !below.held_by.contains(&id.index())
            {
                below.held_by.push(id.index());
            }
        }
    }

    /// Forget the box of every face that holds node `id`, the node itself
    /// included where it is a face.
    fn forget_boxes_above(&mut self, id: TShapeId) {
        if self.nodes.get(id).is_none() {
            return;
        }
        let mut stack: SmallVec<[u32; 8]> = smallvec::smallvec![id.index()];
        while let Some(at) = stack.pop() {
            let Some(node) = self.nodes.at(at).and_then(|id| self.nodes.entry(id)) else {
                continue;
            };
            if let Some(slot) = node.face_box {
                self.face_boxes.forget(slot);
            }
            stack.extend(node.held_by.iter().copied());
        }
    }

    /// Forget every kept face box once the geometry has been rewritten in
    /// place since they were found.
    fn sync_face_boxes(&mut self) {
        let revision = self.geometry.revision();
        if self.face_boxes.revision != revision {
            self.face_boxes.forget_all();
            self.face_boxes.revision = revision;
        }
    }

    /// The slot holding a face's box, while the geometry is as it was found
    /// against; `None` for anything else.
    fn face_box_slot(&self, face: &Shape) -> Option<&std::sync::OnceLock<Aabb>> {
        if self.face_boxes.revision != self.geometry.revision() {
            return None;
        }
        self.nodes
            .entry(face.node())
            .and_then(|node| node.face_box)
            .and_then(|slot| self.face_boxes.slot(slot))
    }

    /// The box kept for a face: every point of the face, its tolerance
    /// included, in the face's placement. `None` where `face` is not a face
    /// of this model or no box has been kept for it yet; a kept box is
    /// found by `ogeom_algo::face_bounds`, and by the operations that read
    /// face boxes, and is kept until the face, its wires, edges or vertices,
    /// or the geometry in place change.
    ///
    /// The box is kept in the face's own frame. Under a placement that
    /// turns the face off the axes it is the box of the turned box, which
    /// holds the face but stands clear of it.
    #[must_use]
    pub fn face_bounds(&self, face: &Shape) -> Option<Aabb> {
        let kept = *self.face_box_slot(face)?.get()?;
        let tolerance = self.node(face)?.data().tolerance()?.get();
        let placement = face.transform(&self.datums).ok()?;
        Some(placed_box(&kept, &placement).expanded(tolerance))
    }

    /// The box kept for a face in its own frame, without its placement or
    /// its tolerance; where none is kept, `find` is asked for it with the
    /// face unplaced, and its answer is kept.
    ///
    /// `find` must hold every point of the face it is given, and is what
    /// [`Model::face_bounds`] reports from then on. Where the face cannot
    /// keep a box (the geometry is being rewritten), `find`'s answer is
    /// returned and not kept.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// if `face` is not a face of this model; as `find` reports.
    pub fn kept_face_box(
        &self,
        face: &Shape,
        find: impl FnOnce(&Shape) -> OgeomResult<Aabb>,
    ) -> OgeomResult<Aabb> {
        if self.node(face).map(TShape::kind) != Some(ShapeType::Face) {
            ogeom_bail!(Construction, "only a face of this model keeps a box");
        }
        let bare = Shape::of(face.node());
        let Some(slot) = self.face_box_slot(face) else {
            return find(&bare);
        };
        if let Some(kept) = slot.get() {
            return Ok(*kept);
        }
        let found = find(&bare)?;
        Ok(self.face_boxes.keep(slot, found))
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
        self.nodes.iter().map(|(id, node)| (id, &node.shape))
    }

    /// Every node that has been given an identity, with it.
    pub fn identities(&self) -> impl Iterator<Item = (TShapeId, EntityId)> {
        self.nodes
            .iter()
            .filter_map(|(id, node)| Some((id, node.identity?)))
    }

    /// Add a vertex.
    pub fn add_vertex(&mut self, data: VertexData) -> Shape {
        Shape::of(self.insert_node(TShape::leaf(ShapeType::Vertex, NodeData::Vertex(data))))
    }

    /// Add a vertex at `point` with the minimum tolerance.
    pub fn add_point(&mut self, point: Point) -> Shape {
        self.add_vertex(VertexData::new(point))
    }

    /// Add an edge bounded by the given vertices.
    ///
    /// The vertices are the edge's ends, in order. A closed edge (a full
    /// circle) names the same vertex twice rather than once, so that walking
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
        // edge is; otherwise the cap does not reliably sit on what it caps.
        for bound in bounds {
            self.widen(bound, data.tolerance)?;
        }
        let node = self.insert_node(TShape::new(
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
        Ok(Shape::of(self.insert_node(TShape::container(
            ShapeType::Wire,
            edges.to_vec(),
        ))))
    }

    /// Add a face bounded by the given wires.
    ///
    /// The first wire is the outer boundary and any others are holes, and
    /// the face keeps them in that order whichever way it is later used; see
    /// [`Model::outer_wire`].
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
        // face. Walk the wires' edges and widen them where they are not; no
        // edge is tighter than the least tolerance.
        let face_tolerance = data.tolerance;
        let walk = if face_tolerance.get() > Tolerance::MIN.get() {
            wires
        } else {
            &[]
        };
        for wire in walk {
            let edges = self.children_of(wire)?;
            for edge in &edges {
                self.widen(edge, face_tolerance)?;
            }
        }
        let node = self.insert_node(TShape::new(
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
        Ok(Shape::of(self.insert_node(TShape::container(
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
        Ok(Shape::of(self.insert_node(TShape::container(
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
        Ok(Shape::of(self.insert_node(TShape::container(
            ShapeType::CompSolid,
            solids.to_vec(),
        ))))
    }

    /// Add a compound of arbitrary shapes.
    ///
    /// The one container with no type constraint; that is what a compound is
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
        Ok(Shape::of(self.insert_node(TShape::container(
            ShapeType::Compound,
            shapes.to_vec(),
        ))))
    }

    /// The direct children of a shape, with this shape's placement and
    /// orientation composed onto each.
    ///
    /// The reason traversal is correct by default rather than by discipline:
    /// a child's placement in the world is its parent's composed with its
    /// own, and its orientation is its parent's composed with its own.
    /// Returning raw children would leave every caller to remember both, and
    /// the failure is silent: face normals that flip inconsistently,
    /// sub-shapes drawn at the origin.
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
            .map(|child| child.beneath(shape.location(), shape.orientation()))
            .collect())
    }

    /// A shape's children in *traversal* order.
    ///
    /// The same shapes as [`Model::children_of`], but with the list reversed
    /// when the parent is reversed. Order carries meaning for a wire (its
    /// edges run head to tail), and reversing a wire has to reverse the walk as
    /// well as each edge, or consecutive edges stop sharing a vertex and the
    /// boundary comes apart. For a shell or a solid the order means nothing and
    /// the reversal is invisible.
    ///
    /// [`Model::children_of`] stays the raw accessor: it returns what is
    /// stored, which is what a rebuild or a comparison wants.
    ///
    /// A reversed face's wires come back reversed too, so its holes are
    /// listed before its outer wire. A face stores its outer wire first, and
    /// that is the only place the order means anything: code that wants the
    /// outer wire asks [`Model::outer_wire`], and code that wants every wire
    /// walked under the face's sense, outer first, asks
    /// [`Model::children_of`], whose wires carry the face's orientation in
    /// stored order.
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

    /// A face's outer wire, carrying the face's placement and orientation,
    /// or `None` for a face with no wires (one covering its whole surface).
    ///
    /// The outer wire is the first one the face stores, whichever way the
    /// face is used. [`Model::ordered_children_of`] on a reversed face lists
    /// the holes first, so its first wire is not this one.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the face does not
    /// resolve in this model; [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// if the shape is not a face.
    pub fn outer_wire(&self, face: &Shape) -> OgeomResult<Option<Shape>> {
        if self.kind_of(face)? != ShapeType::Face {
            ogeom_bail!(Construction, "only a face has an outer wire");
        }
        Ok(self.children_of(face)?.into_iter().next())
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
    /// The nodes are widened where they stand, held by other shapes as well
    /// or not ([`Model::note_held`]): a boolean widens a vertex its rebuilt
    /// edges end on though a face it set aside, and so the operand, holds
    /// it too. Only a tolerance grows: every shape holding the node keeps
    /// its geometry and its pcurves, and a looser bound still bounds what
    /// it bounded, so each stays valid. An editor starting from a shape it
    /// was given copies what is held first ([`Model::unshare`]), and widens
    /// no other shape.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape, or
    /// anything below it, does not resolve in this model.
    pub fn widen(&mut self, shape: &Shape, to: Tolerance) -> OgeomResult<()> {
        // An edge and its vertices, the usual call, are searched in a short
        // list; a larger shape's nodes are hashed.
        const SHORT: usize = 32;
        self.sync_face_boxes();
        let mut affected: smallvec::SmallVec<[TShapeId; 8]> = smallvec::SmallVec::new();
        let mut seen: Option<ogeom_core::FastSet<TShapeId>> = None;
        let mut stack: smallvec::SmallVec<[TShapeId; 8]> = smallvec::smallvec![shape.node()];
        while let Some(id) = stack.pop() {
            let fresh = match &mut seen {
                Some(seen) => seen.insert(id),
                None => !affected.contains(&id),
            };
            if !fresh {
                continue;
            }
            if seen.is_none() && affected.len() >= SHORT {
                seen = Some(affected.iter().copied().chain([id]).collect());
            }
            let Some(node) = self.nodes.get(id) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            affected.push(id);
            stack.extend(node.children().iter().map(Shape::node));
        }
        for id in affected {
            if let Some(node) = self.nodes.get_mut(id) {
                if let Some(journal) = &mut self.widened
                    && let Some(was) = node.data().tolerance()
                    && was.get() < to.get()
                {
                    journal.push((id, was));
                }
                node.data_mut().widen(to);
            }
        }
        Ok(())
    }

    /// Start noting the tolerances [`Model::widen`] grows and edits through
    /// [`Model::node_mut`] may grow, so an attempt that fails can put them
    /// back with [`Model::undo_widened`]. Where a journal is
    /// open already, it goes on and the mark is where this attempt began.
    pub fn note_widened(&mut self) -> WidenMark {
        match &self.widened {
            Some(journal) => WidenMark {
                from: journal.len(),
                opened: false,
            },
            None => {
                self.widened = Some(Vec::new());
                WidenMark {
                    from: 0,
                    opened: true,
                }
            }
        }
    }

    /// Put back every tolerance noted since `mark`, and
    /// close the journal where `mark` opened it.
    pub fn undo_widened(&mut self, mark: WidenMark) {
        let Some(journal) = &mut self.widened else {
            return;
        };
        let undone: Vec<(TShapeId, Tolerance)> =
            journal.drain(mark.from.min(journal.len())..).collect();
        // Latest first, so a node grown twice ends at what it was first.
        for (id, was) in undone.into_iter().rev() {
            if let Some(node) = self.nodes.get_mut(id) {
                node.data_mut().set_tolerance(was);
            }
        }
        if mark.opened {
            self.widened = None;
        }
    }

    /// Keep the tolerances noted since `mark` as they now stand, and close
    /// the journal where `mark` opened it.
    pub fn keep_widened(&mut self, mark: WidenMark) {
        if mark.opened {
            self.widened = None;
        }
    }

    /// Whether node `id` is marked held ([`Model::note_held`]).
    fn is_held(&self, id: TShapeId) -> bool {
        self.nodes.entry(id).is_some_and(|node| node.held)
    }

    /// Record that `result` holds nodes made before the model held `since`
    /// nodes: what an operation passed through into its result from the
    /// shapes it was given, which those shapes still hold.
    ///
    /// The walk goes down from `result` through the nodes made since, and
    /// marks each older node where it first meets one; what lies below a
    /// marked node is held with it. It costs what the operation made.
    /// [`Model::unshare`] reads the marks.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the
    /// shape, or anything below it, does not resolve in this model.
    pub fn note_held(&mut self, result: &Shape, since: usize) -> OgeomResult<()> {
        // Nodes are only appended between counting `since` and here, so the
        // ones made since are the last `made` slots handed out.
        let made = self.nodes.len().saturating_sub(since);
        let next = self.nodes.arena.next_index().unwrap_or(u32::MAX);
        let since = next.saturating_sub(u32::try_from(made).unwrap_or(u32::MAX));
        let mut seen = vec![false; made];
        let mut stack: smallvec::SmallVec<[TShapeId; 16]> = smallvec::smallvec![result.node()];
        while let Some(id) = stack.pop() {
            let Some(node) = self.nodes.get(id) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            let Some(at) = id.index().checked_sub(since) else {
                if let Some(node) = self.nodes.entry_mut(id) {
                    node.held = true;
                    self.any_held = true;
                }
                continue;
            };
            let Some(slot) = seen.get_mut(at as usize) else {
                continue;
            };
            if core::mem::replace(slot, true) {
                continue;
            }
            stack.extend(node.children().iter().map(Shape::node));
        }
        Ok(())
    }

    /// Make every node below `root` the root's alone, so that an edit in
    /// place below it changes no other shape.
    ///
    /// Each node an operation passed through into `root` from the shapes
    /// it was given ([`Model::note_held`]), and everything below it, is
    /// copied: the same kind, data and identity, and the same kept face
    /// box. The nodes `root` holds alone that held one now hold its copy,
    /// at the same placement and orientation. What is held is found from
    /// `root` down; a `root` that is itself held is left as it stands,
    /// with everything below it.
    ///
    /// Returns each node copied, with its copy, in the order copied.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the
    /// shape, or anything below it, does not resolve in this model.
    pub fn unshare(&mut self, root: &Shape) -> OgeomResult<Vec<(TShapeId, TShapeId)>> {
        if !self.any_held || self.is_held(root.node()) {
            return Ok(Vec::new());
        }
        let mut pairs = Vec::new();
        self.unshare_below(std::slice::from_ref(root), true, &mut pairs)?;
        Ok(pairs)
    }

    /// As [`Model::unshare`] for several shapes at once, for an operation
    /// that takes shapes and returns what it makes of them, and leaves the
    /// shapes it was given as they are.
    ///
    /// Nothing given is edited: each node held below a shape, and each
    /// node on the way down to one, the shape itself included, is copied
    /// once for them all, so a shape with anything held below it comes
    /// back as a copy. What a copy shares with the shape it was copied
    /// from is noted as held by both ([`Model::note_held`]).
    ///
    /// # Errors
    ///
    /// As [`Model::unshare`].
    pub fn unshare_each(&mut self, shapes: &[Shape]) -> OgeomResult<Unshared> {
        if !self.any_held {
            return Ok(Unshared {
                shapes: shapes.to_vec(),
                copies: Vec::new(),
            });
        }
        let since = self.nodes.len();
        let mut pairs = Vec::new();
        let copies = self.unshare_below(shapes, false, &mut pairs)?;
        let mut out = Vec::with_capacity(shapes.len());
        for shape in shapes {
            match copies.get(&shape.node()) {
                Some(&copy) => {
                    let copy = Shape::new(copy, shape.location().clone(), shape.orientation());
                    self.note_held(&copy, since)?;
                    out.push(copy);
                }
                None => out.push(shape.clone()),
            }
        }
        Ok(Unshared {
            shapes: out,
            copies: pairs,
        })
    }

    /// Copy every node held below `roots`; where `repoint`, point the nodes
    /// the roots hold alone at the copies, and otherwise copy those on the
    /// way down to a held node too, the roots included. The copies by
    /// original.
    fn unshare_below(
        &mut self,
        roots: &[Shape],
        repoint: bool,
        pairs: &mut Vec<(TShapeId, TShapeId)>,
    ) -> OgeomResult<ogeom_core::FastMap<TShapeId, TShapeId>> {
        // Every node reached, in the order first reached, and whether it is
        // held: marked, or below a held node on some path.
        let mut held: ogeom_core::FastMap<TShapeId, bool> = ogeom_core::FastMap::default();
        let mut order: Vec<TShapeId> = Vec::new();
        let mut stack: Vec<(TShapeId, bool)> =
            roots.iter().rev().map(|r| (r.node(), false)).collect();
        let mut any = false;
        while let Some((id, above)) = stack.pop() {
            let is_held = above || self.is_held(id);
            match held.get(&id) {
                Some(&was) if was || !is_held => continue,
                Some(_) => {}
                None => order.push(id),
            }
            held.insert(id, is_held);
            any |= is_held;
            let Some(node) = self.nodes.get(id) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            stack.extend(node.children().iter().rev().map(|c| (c.node(), is_held)));
        }
        let mut copies: ogeom_core::FastMap<TShapeId, TShapeId> = ogeom_core::FastMap::default();
        if !any {
            return Ok(copies);
        }
        // The nodes to copy: the held ones, and where nothing is pointed
        // anew, every node with one of them below it.
        let mut copied: ogeom_core::FastSet<TShapeId> = held
            .iter()
            .filter(|(_, h)| **h)
            .map(|(id, _)| *id)
            .collect();
        if !repoint {
            loop {
                let mut grew = false;
                for &id in order.iter().rev() {
                    if copied.contains(&id) {
                        continue;
                    }
                    let below = self
                        .nodes
                        .get(id)
                        .is_some_and(|n| n.children().iter().any(|c| copied.contains(&c.node())));
                    if below {
                        copied.insert(id);
                        grew = true;
                    }
                }
                if !grew {
                    break;
                }
            }
        }
        self.sync_face_boxes();
        for root in roots {
            if copied.contains(&root.node()) {
                self.copy_below(root.node(), &copied, &mut copies, pairs)?;
            }
        }
        if !repoint {
            return Ok(copies);
        }
        for parent in order {
            if copied.contains(&parent) {
                continue;
            }
            let Some(node) = self.nodes.get(parent) else {
                continue;
            };
            let children = node.children().to_vec();
            let mut changed = false;
            let mut repointed = Vec::with_capacity(children.len());
            for child in children {
                if copied.contains(&child.node()) {
                    let copy = self.copy_below(child.node(), &copied, &mut copies, pairs)?;
                    repointed.push(Shape::new(
                        copy,
                        child.location().clone(),
                        child.orientation(),
                    ));
                    changed = true;
                } else {
                    repointed.push(child);
                }
            }
            if changed && let Some(node) = self.nodes.get_mut(parent) {
                *node.children_mut() = repointed;
                self.link_children(parent);
            }
        }
        Ok(copies)
    }

    /// The copy of `id`, with every node below it that `copied` names
    /// copied too, each once.
    fn copy_below(
        &mut self,
        id: TShapeId,
        copied: &ogeom_core::FastSet<TShapeId>,
        copies: &mut ogeom_core::FastMap<TShapeId, TShapeId>,
        pairs: &mut Vec<(TShapeId, TShapeId)>,
    ) -> OgeomResult<TShapeId> {
        // Children before parents: a node is copied once every child to
        // copy has its copy.
        let mut stack: Vec<(TShapeId, bool)> = vec![(id, false)];
        while let Some((at, ready)) = stack.pop() {
            if copies.contains_key(&at) {
                continue;
            }
            let Some(node) = self.nodes.get(at) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            if !ready {
                stack.push((at, true));
                stack.extend(
                    node.children()
                        .iter()
                        .rev()
                        .filter(|c| copied.contains(&c.node()) && !copies.contains_key(&c.node()))
                        .map(|c| (c.node(), false)),
                );
                continue;
            }
            let children = node
                .children()
                .iter()
                .map(|c| {
                    let copy = copies.get(&c.node()).copied().unwrap_or(c.node());
                    Shape::new(copy, c.location().clone(), c.orientation())
                })
                .collect();
            let copy = TShape::new(node.kind(), node.data().clone(), children);
            let kept = self
                .face_box_slot(&Shape::of(at))
                .and_then(|slot| slot.get().copied());
            let made = self.insert_node(copy);
            let identity = self.nodes.entry(at).and_then(|node| node.identity);
            if let Some(node) = self.nodes.entry_mut(made) {
                node.identity = identity;
            }
            if let (Some(kept), Some(slot)) = (kept, self.face_box_slot(&Shape::of(made))) {
                self.face_boxes.keep(slot, kept);
            }
            copies.insert(at, made);
            pairs.push((at, made));
        }
        copies.get(&id).copied().ok_or_else(|| {
            ogeom_core::ogeom_err!(Dangling, "shape refers to a node not in this model")
        })
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
        // Containers carry no tolerance of their own, so each node is held
        // to the nearest tolerance above it, whatever lies between: a face's
        // edges answer to the face across the wire that holds them. Each
        // node is visited once per bound it is reached under.
        let mut stack: Vec<(Shape, Option<(ogeom_core::Tolerance, ShapeType)>)> =
            vec![(root.clone(), None)];
        let mut seen = ogeom_core::FastSet::default();
        while let Some((shape, bound)) = stack.pop() {
            let Some(node) = self.node(&shape) else {
                ogeom_bail!(Dangling, "shape refers to a node not in this model");
            };
            if !seen.insert((shape.node(), bound.map(|(t, _)| t.get().to_bits()))) {
                continue;
            }
            let own = node.data().tolerance();
            // A boundary is *contained by* what it bounds, so the child (the
            // boundary) must be the looser of the two.
            if let (Some((parent, parent_kind)), Some(child_tolerance)) = (bound, own)
                && child_tolerance < parent
            {
                ogeom_bail!(
                    Invariant,
                    "a {:?} at tolerance {} bounds a {:?} at {}, which is tighter",
                    node.kind(),
                    child_tolerance.get(),
                    parent_kind,
                    parent.get()
                );
            }
            let below = own.map(|t| (t, node.kind())).or(bound);
            for child in node.children() {
                stack.push((child.clone(), below));
            }
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
    /// underlying node (and all its geometry) is shared rather than copied.
    /// Placing ten thousand instances of a part costs ten thousand short chains
    /// and one copy of the geometry.
    pub fn placed(&mut self, shape: &Shape, transform: Transform) -> Shape {
        let datum = self.add_datum(transform);
        shape.moved(&Location::of(datum))
    }
}

/// A box carried by a placement: exactly, where the placement keeps the
/// axes; as the box of the carried box otherwise.
fn placed_box(kept: &Aabb, placement: &Transform) -> Aabb {
    match placement.kind() {
        TransformKind::Identity => *kept,
        TransformKind::Translation => match (kept.low(), kept.high()) {
            (Some(low), Some(high)) => {
                Aabb::of_corners(placement.apply(low), placement.apply(high))
            }
            _ => *kept,
        },
        _ => kept.transformed(placement),
    }
}

/// Where an attempt began in the journal [`Model::note_widened`] opens.
#[derive(Debug, Clone, Copy)]
#[must_use]
pub struct WidenMark {
    from: usize,
    opened: bool,
}

/// What [`Model::unshare_each`] made of the shapes it was given.
#[derive(Debug, Clone)]
pub struct Unshared {
    /// The shapes to work on, in the order given: each as given, or its
    /// copy at the same placement and orientation.
    pub shapes: Vec<Shape>,
    /// Each node copied, with its copy, in the order copied.
    pub copies: Vec<(TShapeId, TShapeId)>,
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
    #[allow(
        clippy::disallowed_types,
        reason = "public field since 0.1; a std map keeps the API"
    )]
    pub entities: std::collections::HashMap<EntityId, EntityId>,
}

/// A model's contents, laid out the way a file holds them.
///
/// Handed to [`Model::from_parts`]. Every list is in arena order, and a node's
/// children name other nodes by their position in `nodes`, so the order is
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
/// orientation: the composition of everything from the root down. That is the
/// only form in which a sub-shape means anything outside the tree it came from.
///
/// Sub-shapes reached by more than one route (an edge shared by two faces)
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
        // The node is fetched once and answers both questions. `kind_of` would
        // repeat the arena lookup that `children_of` is about to do anyway,
        // and this walk is the kernel's most travelled road.
        let Some(node) = model.node(&shape) else {
            ogeom_bail!(Dangling, "shape refers to a node not in this model");
        };
        let (matches, descend) = match filter {
            // Every child of anything but a compound is of a lower type (the
            // builders refuse anything else), so below a match, or below a
            // shape already lower than the one wanted, nothing can match:
            // a walk for faces never enters a wire.
            Filter::OfType(want) => (
                node.kind() == want,
                node.kind() == ShapeType::Compound || node.kind() > want,
            ),
            Filter::All => (true, true),
        };
        let children = if descend { node.children() } else { &[] };
        stack.reserve(children.len());
        // `children_of`'s composition, inline: the parent's placement and
        // sense onto each child. Done here so the walk does not allocate a
        // `Vec` per node only to drain it.
        for child in children.iter().rev() {
            stack.push(child.beneath(shape.location(), shape.orientation()));
        }
        if matches {
            // `shape` is owned and finished with; cloning it to keep it would
            // copy a location chain for nothing.
            out.push(shape);
        }
    }
    Ok(out)
}

/// Walk a shape's tree, yielding every sub-shape of `want` exactly once.
///
/// Deduplicated by [`Shape::is_same`] (node and placement, ignoring
/// orientation), which is what "the distinct edges of this solid" means.
///
/// # Errors
///
/// As [`explore`].
pub fn explore_unique(model: &Model, root: &Shape, want: ShapeType) -> OgeomResult<Vec<Shape>> {
    use crate::shape::SameKey;

    let found = explore(model, root, Filter::OfType(want))?;
    let mut seen = ogeom_core::FastSet::with_capacity_and_hasher(found.len(), Default::default());
    let mut out = Vec::with_capacity(found.len());
    for shape in found {
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
    // Each candidate is searched for the target's own type only: a walk
    // for an edge never enters a vertex list.
    let kind = model.kind_of(target)?;
    let mut out = Vec::new();
    for candidate in explore(model, root, Filter::OfType(want))? {
        if explore(model, &candidate, Filter::OfType(kind))?
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
        // edge 1 would end where it started while edge 2 still starts where
        // it did, so consecutive edges stop meeting and a face built
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
        // refusing, the builder widens the vertex; tolerances only ever grow,
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
    fn a_kept_face_box_is_forgotten_when_what_the_face_is_made_of_changes() {
        let unit = Aabb::of_corners(Point::ORIGIN, Point::new(1.0, 1.0, 0.0));
        let keep = |model: &Model, face: &Shape| {
            model.kept_face_box(face, |_| Ok(unit)).unwrap();
        };
        let mut model = Model::new();
        let face = square(&mut model);
        assert_eq!(model.face_bounds(&face), None);
        keep(&model, &face);
        let kept = model.face_bounds(&face).unwrap();
        assert!(kept.contains_box(&unit));

        // A placement carries the box; a tolerance widens it.
        let lift = model.placed(&face, Transform::translation(ogeom_math::Vector::Z));
        let lifted = model.face_bounds(&lift).unwrap();
        assert!((lifted.low().unwrap().z - 1.0).abs() < 1e-6);
        model.widen(&face, Tolerance::new(0.25).unwrap()).unwrap();
        let wide = model.face_bounds(&face).unwrap();
        assert!((wide.low().unwrap().x + 0.25).abs() < 1e-12);

        // A vertex handed out for editing forgets the box of the face above.
        let vertex = explore_unique(&model, &face, ShapeType::Vertex).unwrap()[0].clone();
        let _ = model.node_mut(&vertex);
        assert_eq!(model.face_bounds(&face), None);

        // So does a face handed out itself, and a surface rewritten in place,
        // which reaches every face.
        keep(&model, &face);
        let _ = model.node_mut(&face);
        assert_eq!(model.face_bounds(&face), None);
        keep(&model, &face);
        let Some(NodeData::Face(data)) = model.node(&face).map(TShape::data) else {
            panic!("a face holds face data");
        };
        let surface = data.surface;
        let _ = model.geometry_mut().surface_mut(surface);
        assert_eq!(model.face_bounds(&face), None);
        let _ = model.add_point(Point::ORIGIN);
        keep(&model, &face);
        assert!(model.face_bounds(&face).is_some());

        // A node no face holds forgets nothing.
        let alone = model.add_point(Point::ORIGIN);
        let _ = model.node_mut(&alone);
        assert!(model.face_bounds(&face).is_some());
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
        // around it, which is exactly what happens when topology arrives from
        // a file. The check has to stand on its own, or imported geometry sails
        // past it.
        let mut model = Model::new();
        let vertex = Shape::of(model.insert_node(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        )));

        let mut edge_data = EdgeData::new();
        edge_data.widen(Tolerance::new(1e-1).unwrap());
        let edge = Shape::of(model.insert_node(TShape::new(
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

    /// A face looser than its edges breaks the rule across the wire between
    /// them, which carries no tolerance of its own to be compared with.
    #[test]
    fn check_tolerances_holds_a_face_to_its_edges_across_the_wire() {
        let mut model = Model::new();
        let vertex = Shape::of(model.insert_node(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        )));
        let edge = Shape::of(model.insert_node(TShape::new(
            ShapeType::Edge,
            NodeData::Edge(Box::default()),
            vec![vertex],
        )));
        let wire =
            Shape::of(model.insert_node(TShape::container(ShapeType::Wire, vec![edge.clone()])));
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let mut face_data = FaceData::new(surface, Location::identity());
        face_data.widen(Tolerance::new(1e-2).unwrap());
        let face = Shape::of(model.insert_node(TShape::new(
            ShapeType::Face,
            NodeData::Face(Box::new(face_data)),
            vec![wire],
        )));
        let err = model.check_tolerances(&face).unwrap_err();
        assert!(err.to_string().contains("tighter"), "{err}");
        model.widen(&face, Tolerance::new(1e-2).unwrap()).unwrap();
        assert!(model.check_tolerances(&face).is_ok());
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
    fn what_a_failed_attempt_widened_is_put_back() {
        let mut model = Model::new();
        let tolerance =
            |model: &Model, v: &Shape| model.node(v).unwrap().data().tolerance().unwrap().get();
        let v = model.add_point(Point::new(0.0, 0.0, 0.0));
        let w = model.add_point(Point::new(1.0, 0.0, 0.0));
        let was = tolerance(&model, &v);
        let outer = model.note_widened();
        model.widen(&v, Tolerance::new(0.5).unwrap()).unwrap();
        // An attempt within another puts back only its own.
        let inner = model.note_widened();
        model.widen(&v, Tolerance::new(2.0).unwrap()).unwrap();
        model.widen(&w, Tolerance::new(3.0).unwrap()).unwrap();
        model.undo_widened(inner);
        assert_eq!(tolerance(&model, &v), 0.5);
        assert_eq!(tolerance(&model, &w), was);
        model.undo_widened(outer);
        assert_eq!(tolerance(&model, &v), was);
        // An edit through the node itself is put back as well.
        let edit = model.note_widened();
        if let Some(node) = model.node_mut(&w)
            && let NodeData::Vertex(data) = node.data_mut()
        {
            data.tolerance = Tolerance::new(4.0).unwrap();
        }
        model.undo_widened(edit);
        assert_eq!(tolerance(&model, &w), was);
        // What an attempt that holds keeps stays, and nothing is noted after.
        let kept = model.note_widened();
        model.widen(&v, Tolerance::new(0.25).unwrap()).unwrap();
        model.keep_widened(kept);
        let after = model.note_widened();
        model.undo_widened(after);
        assert_eq!(tolerance(&model, &v), 0.25);
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
