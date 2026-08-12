//! The shape triple, orientation, and the identity trichotomy.
//!
//! `docs/DATA_MODEL.md` §1, §3 and §4. A [`Shape`] is a *(topology node,
//! placement, orientation)* triple: cheap to copy, with the heavy data — the
//! children, the geometry, the tolerances — living once in an arena behind the
//! node handle. That separation is why boundary representation scales: the same
//! node appears at many placements and orientations without a byte of geometry
//! being copied.
//!
//! Two consequences run through everything downstream.
//!
//! **Traversal composes.** A sub-shape's effective placement is the product of
//! every placement from the root down, and its effective orientation is the
//! composition of every orientation on that path. An explorer that yields
//! sub-shapes without composing both is wrong in a way that produces
//! plausible-looking garbage rather than an error.
//!
//! **Identity is three questions, not one.** See [`Shape::is_partner`],
//! [`Shape::is_same`] and [`Shape::is_equal`].

use core::hash::{Hash, Hasher};

use ogeom_core::{Key, OgeomResult, Tolerances};
use ogeom_math::Transform;

use crate::entity::NodeData;
use crate::location::{DatumStore, Location};

/// What a topology node is.
///
/// Ordered by dimension, so `>=` asks a meaningful question — "is this at least
/// a face?" — and sorting a mixed collection groups it sensibly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ShapeType {
    /// A point.
    Vertex,
    /// A curve bounded by vertices.
    Edge,
    /// A connected sequence of edges.
    Wire,
    /// A surface bounded by wires.
    Face,
    /// A connected set of faces.
    Shell,
    /// A volume bounded by shells.
    Solid,
    /// A connected set of solids sharing faces.
    CompSolid,
    /// An arbitrary collection, of any dimensions.
    Compound,
}

impl ShapeType {
    /// The type of the sub-shapes this type is built from, if any.
    ///
    /// A compound has no fixed child type — it holds anything — so this reports
    /// `None` for it rather than guessing.
    #[must_use]
    pub const fn child_type(self) -> Option<Self> {
        match self {
            Self::Vertex | Self::Compound => None,
            Self::Edge => Some(Self::Vertex),
            Self::Wire => Some(Self::Edge),
            Self::Face => Some(Self::Wire),
            Self::Shell => Some(Self::Face),
            Self::Solid => Some(Self::Shell),
            Self::CompSolid => Some(Self::Solid),
        }
    }

    /// The topological dimension: 0 for a vertex, 3 for a solid.
    ///
    /// A compound reports `None`, since it may mix dimensions.
    #[must_use]
    pub const fn dimension(self) -> Option<u8> {
        match self {
            Self::Vertex => Some(0),
            Self::Edge | Self::Wire => Some(1),
            Self::Face | Self::Shell => Some(2),
            Self::Solid | Self::CompSolid => Some(3),
            Self::Compound => None,
        }
    }
}

/// Which side of a boundary the material is on.
///
/// `docs/DATA_MODEL.md` §3.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Orientation {
    /// The material is on the surface's default side.
    #[default]
    Forward,
    /// The material is on the other side.
    Reversed,
    /// The boundary lies *inside* the material — a stiffener edge embedded in a
    /// face, an edge that does not separate anything.
    Internal,
    /// The boundary lies outside the material: reference geometry, carried
    /// along but bounding nothing.
    External,
}

impl Orientation {
    /// Compose an outer orientation with an inner one.
    ///
    /// Applied at *every* level of descent: an edge's orientation within a face
    /// depends on that face's orientation within its shell, and so on to the
    /// root. Reversing a solid must therefore not touch a single child.
    ///
    /// `Internal` and `External` absorb: a boundary that lies inside the
    /// material stays inside it however the shape around it is turned.
    #[must_use]
    pub const fn compose(self, inner: Self) -> Self {
        match self {
            Self::Forward => inner,
            Self::Reversed => match inner {
                Self::Forward => Self::Reversed,
                Self::Reversed => Self::Forward,
                other => other,
            },
            Self::Internal => Self::Internal,
            Self::External => Self::External,
        }
    }

    /// This orientation reversed.
    ///
    /// `Internal` and `External` are unaffected — neither names a side, so
    /// neither has one to swap.
    #[must_use]
    pub const fn reversed(self) -> Self {
        match self {
            Self::Forward => Self::Reversed,
            Self::Reversed => Self::Forward,
            other => other,
        }
    }

    /// Whether this orientation names a side of the material at all.
    #[must_use]
    pub const fn is_boundary(self) -> bool {
        matches!(self, Self::Forward | Self::Reversed)
    }
}

/// A topology node: the shared, placeless, orientationless part of a shape.
///
/// Held in an arena; a [`Shape`] refers to one by handle. The geometry,
/// tolerances and edge representations that hang off a node live alongside it,
/// keyed by the same handle.
#[derive(Debug, Clone, PartialEq)]
pub struct TShape {
    kind: ShapeType,
    data: NodeData,
    children: Vec<Shape>,
}

/// A handle to a [`TShape`].
pub type TShapeId = Key<TShape>;

impl TShape {
    /// A node of the given type, with its data and children.
    #[must_use]
    pub const fn new(kind: ShapeType, data: NodeData, children: Vec<Shape>) -> Self {
        Self {
            kind,
            data,
            children,
        }
    }

    /// A container node — wire, shell, solid, compsolid or compound.
    #[must_use]
    pub const fn container(kind: ShapeType, children: Vec<Shape>) -> Self {
        Self {
            kind,
            data: NodeData::Container,
            children,
        }
    }

    /// A node with data and no children.
    #[must_use]
    pub const fn leaf(kind: ShapeType, data: NodeData) -> Self {
        Self {
            kind,
            data,
            children: Vec::new(),
        }
    }

    /// What kind of node this is.
    #[must_use]
    pub const fn kind(&self) -> ShapeType {
        self.kind
    }

    /// The geometry and tolerance this node carries.
    #[must_use]
    pub const fn data(&self) -> &NodeData {
        &self.data
    }

    /// Mutable access to the node's data.
    #[must_use]
    pub const fn data_mut(&mut self) -> &mut NodeData {
        &mut self.data
    }

    /// The direct children, in order.
    #[must_use]
    pub fn children(&self) -> &[Shape] {
        &self.children
    }

    /// The children, mutably, for binding restored handles to their arena.
    pub(crate) const fn children_mut(&mut self) -> &mut Vec<Shape> {
        &mut self.children
    }

    /// Number of direct children.
    #[must_use]
    pub fn child_count(&self) -> usize {
        self.children.len()
    }
}

/// A shape: a topology node, a placement, and an orientation.
///
/// Cheap to copy — a key, a small chain and an enum — so it is passed by value
/// everywhere.
#[derive(Debug, Clone)]
pub struct Shape {
    node: TShapeId,
    location: Location,
    orientation: Orientation,
}

impl Shape {
    /// A shape from its three parts.
    #[must_use]
    pub const fn new(node: TShapeId, location: Location, orientation: Orientation) -> Self {
        Self {
            node,
            location,
            orientation,
        }
    }

    /// A shape at the identity placement, oriented forward.
    #[must_use]
    pub fn of(node: TShapeId) -> Self {
        Self::new(node, Location::identity(), Orientation::Forward)
    }

    /// The topology node.
    #[must_use]
    pub const fn node(&self) -> TShapeId {
        self.node
    }

    /// The placement.
    #[must_use]
    pub const fn location(&self) -> &Location {
        &self.location
    }

    /// The orientation.
    #[must_use]
    pub const fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// This shape with a different placement.
    #[must_use]
    pub fn located(&self, location: Location) -> Self {
        Self {
            location,
            ..self.clone()
        }
    }

    /// This shape with its handles bound to the arenas that hold them.
    ///
    /// For reading a document back, where the handles are rebuilt before the
    /// arenas exist. Nothing else should need it.
    pub(crate) fn rebound(&self, nodes: u32, datums: u32) -> Self {
        Self {
            node: self.node.with_scope(nodes),
            location: self.location.with_datum_scope(datums),
            orientation: self.orientation,
        }
    }

    /// This shape with its handles shifted for landing in a live model.
    ///
    /// [`rebound`](Self::rebound)'s sibling for absorbing parts: the indices
    /// were local to the source document, and its nodes and datums are about
    /// to land `nodes` and `datums` slots into the target's arenas. The
    /// handles stay unscoped — binding is a separate, later step.
    pub(crate) fn shifted(&self, nodes: u32, datums: u32) -> Self {
        Self {
            node: crate::entity::shifted_key(self.node, nodes),
            location: self.location.with_datum_offset(datums),
            orientation: self.orientation,
        }
    }

    /// This shape moved by `outer`, applied before its own placement.
    ///
    /// The operation traversal uses on descent: a child's placement within the
    /// world is its parent's composed with its own.
    #[must_use]
    pub fn moved(&self, outer: &Location) -> Self {
        Self {
            location: outer.then(&self.location),
            ..self.clone()
        }
    }

    /// This shape with a different orientation.
    #[must_use]
    pub fn oriented(&self, orientation: Orientation) -> Self {
        Self {
            orientation,
            ..self.clone()
        }
    }

    /// This shape with its orientation reversed.
    #[must_use]
    pub fn reversed(&self) -> Self {
        self.oriented(self.orientation.reversed())
    }

    /// This shape's orientation composed under `outer`.
    ///
    /// The other half of what traversal does on descent.
    #[must_use]
    pub fn composed(&self, outer: Orientation) -> Self {
        self.oriented(outer.compose(self.orientation))
    }

    /// The composed placement as a transform.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`].
    pub fn transform(&self, store: &DatumStore) -> OgeomResult<Transform> {
        self.location.composed(store)
    }

    /// Whether two shapes share a topology node, ignoring placement and
    /// orientation.
    ///
    /// "Is this the same underlying topology, anywhere, any way round?" — the
    /// question to ask when relating a shape to another instance of itself
    /// elsewhere in an assembly.
    #[must_use]
    pub fn is_partner(&self, other: &Self) -> bool {
        self.node == other.node
    }

    /// Whether two shapes share a node *and* a placement, ignoring orientation.
    ///
    /// The common case, and the one most algorithms want: an edge and its
    /// reverse are the same edge in the same place, and a set of edges should
    /// hold one of them, not two.
    #[must_use]
    pub fn is_same(&self, other: &Self) -> bool {
        self.node == other.node && self.location == other.location
    }

    /// Whether two shapes agree in all three parts.
    ///
    /// Exact identity. An edge and its reverse are *not* equal, which is what
    /// makes a wire's direction of travel expressible.
    #[must_use]
    pub fn is_equal(&self, other: &Self) -> bool {
        self.is_same(other) && self.orientation == other.orientation
    }

    /// Whether two shapes occupy the same place, comparing composed transforms
    /// rather than chains.
    ///
    /// Costlier than [`Shape::is_same`] and answers a different question: two
    /// placements built by different routes can land in the same position.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`].
    pub fn is_same_position(
        &self,
        other: &Self,
        store: &DatumStore,
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        Ok(self.node == other.node
            && self
                .location
                .is_same_placement(&other.location, store, tol)?)
    }
}

/// Equality by [`Shape::is_equal`] — node, placement *and* orientation.
///
/// The strictest of the three, chosen as the derive-shaped default so that a
/// plain `==` never silently means something looser than the reader expects.
/// Code that wants a weaker equivalence says so, through [`SameKey`] or
/// [`PartnerKey`].
impl PartialEq for Shape {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Eq for Shape {}

impl Hash for Shape {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node.hash(state);
        self.location.hash(state);
        self.orientation.hash(state);
    }
}

/// A key that hashes and compares by [`Shape::is_same`] — node and placement,
/// ignoring orientation.
///
/// Wrapping rather than offering a custom hasher, because the danger being
/// guarded against is a map whose comparison and hash disagree: a set keyed on
/// "same" semantics but hashing the orientation as well will silently hold both
/// an edge and its reverse. Making the equivalence part of the *type* means the
/// two can never drift apart.
#[derive(Debug, Clone)]
pub struct SameKey(pub Shape);

impl PartialEq for SameKey {
    fn eq(&self, other: &Self) -> bool {
        self.0.is_same(&other.0)
    }
}

impl Eq for SameKey {}

impl Hash for SameKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.node.hash(state);
        self.0.location.hash(state);
    }
}

/// A key that hashes and compares by [`Shape::is_partner`] — the node alone.
#[derive(Debug, Clone)]
pub struct PartnerKey(pub Shape);

impl PartialEq for PartnerKey {
    fn eq(&self, other: &Self) -> bool {
        self.0.is_partner(&other.0)
    }
}

impl Eq for PartnerKey {}

impl Hash for PartnerKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.node.hash(state);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::entity::VertexData;
    use ogeom_core::Arena;
    use ogeom_math::{Point, Vector};
    use std::collections::HashSet;

    fn setup() -> (Arena<TShape>, DatumStore, TShapeId, TShapeId, Location) {
        let mut arena = Arena::new();
        let a = arena.insert(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        ));
        let b = arena.insert(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        ));
        let mut store = DatumStore::new();
        let datum = store.insert(Transform::translation(Vector::new(1.0, 0.0, 0.0)));
        (arena, store, a, b, Location::of(datum))
    }

    #[test]
    fn orientation_composition_is_a_monoid_with_forward_as_identity() {
        use Orientation::{External, Forward, Internal, Reversed};
        for o in [Forward, Reversed, Internal, External] {
            assert_eq!(Forward.compose(o), o, "forward is a left identity");
            assert_eq!(o.compose(Forward), o, "and a right identity");
        }
        assert_eq!(Reversed.compose(Reversed), Forward);
        assert_eq!(Reversed.compose(Forward), Reversed);
    }

    #[test]
    fn composition_is_associative() {
        use Orientation::{External, Forward, Internal, Reversed};
        let all = [Forward, Reversed, Internal, External];
        for a in all {
            for b in all {
                for c in all {
                    assert_eq!(
                        a.compose(b).compose(c),
                        a.compose(b.compose(c)),
                        "{a:?} {b:?} {c:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn internal_and_external_absorb_and_do_not_reverse() {
        use Orientation::{External, Forward, Internal, Reversed};
        // A boundary inside the material stays inside it however the shape
        // around it is turned.
        for o in [Forward, Reversed, Internal, External] {
            assert_eq!(Internal.compose(o), Internal);
            assert_eq!(External.compose(o), External);
        }
        assert_eq!(Internal.reversed(), Internal);
        assert_eq!(External.reversed(), External);
        assert!(!Internal.is_boundary());
        assert!(!External.is_boundary());
        assert!(Forward.is_boundary() && Reversed.is_boundary());
    }

    #[test]
    fn reversal_is_an_involution() {
        use Orientation::{External, Forward, Internal, Reversed};
        for o in [Forward, Reversed, Internal, External] {
            assert_eq!(o.reversed().reversed(), o);
        }
    }

    #[test]
    fn the_identity_trichotomy_distinguishes_three_questions() {
        let (_, _, a, b, loc) = setup();

        let base = Shape::of(a);
        let reversed = base.reversed();
        let moved = base.located(loc.clone());
        let other_node = Shape::of(b);

        // Same node, same place, opposite orientation.
        assert!(base.is_partner(&reversed));
        assert!(base.is_same(&reversed));
        assert!(
            !base.is_equal(&reversed),
            "orientation must distinguish them"
        );

        // Same node, different place.
        assert!(base.is_partner(&moved));
        assert!(!base.is_same(&moved), "placement must distinguish them");
        assert!(!base.is_equal(&moved));

        // Different node.
        assert!(!base.is_partner(&other_node));
        assert!(!base.is_same(&other_node));
        assert!(!base.is_equal(&other_node));

        // And each is reflexive.
        assert!(base.is_equal(&base) && base.is_same(&base) && base.is_partner(&base));
    }

    #[test]
    fn each_key_type_hashes_consistently_with_its_own_equality() {
        // The failure this guards: a set whose comparison and hash disagree
        // holds duplicates it believes it has excluded. Making the equivalence
        // part of the type is what keeps them together.
        let (_, _, a, _, loc) = setup();
        let base = Shape::of(a);
        let reversed = base.reversed();
        let moved = base.located(loc);

        let equal: HashSet<Shape> = [base.clone(), reversed.clone(), moved.clone()]
            .into_iter()
            .collect();
        assert_eq!(equal.len(), 3, "all three differ under strict equality");

        let same: HashSet<SameKey> = [base.clone(), reversed.clone(), moved.clone()]
            .into_iter()
            .map(SameKey)
            .collect();
        assert_eq!(same.len(), 2, "orientation is ignored, placement is not");

        let partner: HashSet<PartnerKey> = [base, reversed, moved]
            .into_iter()
            .map(PartnerKey)
            .collect();
        assert_eq!(partner.len(), 1, "only the node matters");
    }

    #[test]
    fn a_shape_placed_and_composed_reports_the_expected_position() {
        let (_, store, a, _, loc) = setup();
        let shape = Shape::of(a).located(loc.clone());
        assert!(
            shape
                .transform(&store)
                .unwrap()
                .apply(Point::ORIGIN)
                .is_equal(Point::new(1.0, 0.0, 0.0), Tolerances::millimetres())
        );

        // Moving under an outer placement composes rather than replaces.
        let nested = shape.moved(&loc);
        assert!(
            nested
                .transform(&store)
                .unwrap()
                .apply(Point::ORIGIN)
                .is_equal(Point::new(2.0, 0.0, 0.0), Tolerances::millimetres())
        );
        assert_eq!(nested.location().depth(), 1, "same datum, merged power");
    }

    #[test]
    fn positions_reached_by_different_routes_compare_equal() {
        let mut store = DatumStore::new();
        let a = store.insert(Transform::translation(Vector::X));
        let b = store.insert(Transform::translation(Vector::X));
        let mut arena = Arena::new();
        let node = arena.insert(TShape::leaf(
            ShapeType::Vertex,
            NodeData::Vertex(VertexData::new(Point::ORIGIN)),
        ));

        let via_a = Shape::of(node).located(Location::of(a));
        let via_b = Shape::of(node).located(Location::of(b));
        let tol = Tolerances::millimetres();

        assert!(!via_a.is_same(&via_b), "structurally different chains");
        assert!(
            via_a.is_same_position(&via_b, &store, tol).unwrap(),
            "but the same position"
        );
    }

    #[test]
    fn shape_types_report_their_children_and_dimensions() {
        use ShapeType::{CompSolid, Compound, Edge, Face, Shell, Solid, Vertex, Wire};
        assert_eq!(Vertex.child_type(), None);
        assert_eq!(Edge.child_type(), Some(Vertex));
        assert_eq!(Face.child_type(), Some(Wire));
        assert_eq!(Solid.child_type(), Some(Shell));
        assert_eq!(Compound.child_type(), None, "a compound holds anything");

        assert_eq!(Vertex.dimension(), Some(0));
        assert_eq!(Wire.dimension(), Some(1));
        assert_eq!(Shell.dimension(), Some(2));
        assert_eq!(CompSolid.dimension(), Some(3));
        assert_eq!(Compound.dimension(), None, "a compound may mix dimensions");

        // The ordering makes "at least a face" expressible.
        assert!(Face > Edge && Solid > Face);
        assert!(Face >= Face);
    }

    #[test]
    fn a_node_carries_its_children_in_order() {
        let (mut arena, _, a, b, _) = setup();
        let wire = arena.insert(TShape::container(
            ShapeType::Wire,
            vec![Shape::of(a), Shape::of(b).reversed()],
        ));
        let node = arena.get(wire).unwrap();
        assert_eq!(node.kind(), ShapeType::Wire);
        assert_eq!(node.child_count(), 2);
        assert_eq!(node.children()[1].orientation(), Orientation::Reversed);
    }
}
