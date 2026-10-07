//! The faces of each solid the other cannot reach, set aside before the
//! general fuse and passed through it untouched.
//!
//! A face whose box, widened by every tolerance it states and by the
//! margin the pair filters allow, misses the other solid's box lies
//! outside the other solid, and no section, contact or junction can reach
//! it. It is not gathered, split, classified or rebuilt. The faces the
//! fuse does gather keep the edges they share with it as they are, and a
//! hole of a gathered plane whose every edge is shared with such faces,
//! and whose box misses the other solid's, is left out of the arrangement
//! and put back into the piece holding it.
//!
//! A closed edge (a full circle) is split in two by the arrangement of
//! every face holding it, so a face set aside sharing one with a ring the
//! arrangement takes is gathered after all, and the edge is rebuilt in two
//! as it would be with every face taken.

use hashbrown::{HashMap, HashSet};
use smallvec::SmallVec;

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::Aabb;
use ogeom_topo::{Model, NodeData, SameKey, Shape, ShapeType, TShape, TShapeId, explore_unique};

use crate::EdgeKey;

/// What one solid sets aside.
#[derive(Default)]
pub(crate) struct Aside {
    /// The faces set aside, as the solid holds them, each with its place
    /// among the solid's faces.
    pub(crate) faces: Vec<(Shape, usize)>,
    /// The faces the fuse gathers, by node and placement, each with its
    /// place among the solid's faces.
    pub(crate) gathered: HashMap<SameKey, usize>,
    /// The edges the gathered faces share with the faces set aside.
    pub(crate) edges: HashSet<EdgeKey>,
    /// The holes of each gathered plane left out of its arrangement, by
    /// the face's node and placement.
    pub(crate) holes: HashMap<SameKey, std::sync::Arc<[Hole]>>,
    /// Whether the solid's every edge is walked an even number of times,
    /// read once for its faces: for a solid of one shell, that it is closed.
    pub(crate) closed: bool,
    /// The faces holding each edge, and which faces are set aside.
    pub(crate) beside: std::sync::Arc<Beside>,
}

/// A hole of a gathered plane left out of its arrangement.
#[derive(Clone)]
pub(crate) struct Hole {
    /// The wire, as the face's node holds it.
    pub(crate) wire: Shape,
    /// Its place among the face's wires.
    pub(crate) at: usize,
    /// The faces set aside holding its edges, by their place in the solid.
    pub(crate) beside: SmallVec<[u32; 2]>,
    /// The widest doubt its edges carry to their face's arrangement: twice
    /// an edge's tolerance, or its tolerance and its ends'.
    pub(crate) doubt: f64,
}

/// The faces holding each edge of a solid, and which of its faces are set
/// aside: what joins the faces set aside to each other and to the faces
/// rebuilt, read without walking the faces set aside again.
#[derive(Default)]
pub(crate) struct Beside {
    holders: HashMap<EdgeKey, Holders>,
    aside: Vec<bool>,
    /// The edges the gathered faces share with the faces set aside, which
    /// the faces made must walk; a hole's edges are left out, walked by the
    /// face taking it back.
    shared: Vec<EdgeKey>,
}

impl Aside {
    /// Whether anything is set aside.
    pub(crate) fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }
}

/// One face as the setting aside reads it.
struct Read {
    face: Shape,
    /// The face's box widened by its margin.
    bound: Aabb,
    /// Whether the face stands where its node puts it, unmoved.
    unmoved: bool,
    plane: bool,
    /// A plane's wires as walked, where it has holes.
    rings: Option<Rings>,
}

/// A face's wires as walked: each one's edges and the vertices they end
/// on.
#[derive(Default)]
struct Rings {
    keys: Vec<EdgeKey>,
    ends: Vec<TShapeId>,
    /// Each wire's run of `keys` and of `ends`.
    spans: Vec<((usize, usize), (usize, usize))>,
}

impl Rings {
    /// The wires of `face` walked.
    fn of(model: &Model, face: &Shape) -> OgeomResult<Self> {
        let mut rings = Self::default();
        let Some(node) = model.node(face) else {
            ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
        };
        for wire in node.children() {
            let wire = placed_wire(face, wire);
            let (k0, v0) = (rings.keys.len(), rings.ends.len());
            for_each_wire_edge(model, &wire, |edge, node| {
                rings.keys.push(EdgeKey::of(edge));
                rings.ends.extend(node.children().iter().map(Shape::node));
            })?;
            rings
                .spans
                .push(((k0, rings.keys.len()), (v0, rings.ends.len())));
        }
        Ok(rings)
    }
}

/// The faces holding an edge, and whether it is closed: one vertex at
/// both ends.
#[derive(Clone, Copy)]
struct Holders {
    faces: [u32; 2],
    count: u8,
    closed: bool,
    /// The looser of the edge's tolerance and its vertices'.
    margin: f64,
    /// Twice the edge's tolerance, or its tolerance and its ends', the
    /// wider.
    doubt: f64,
    /// How many times the faces walk it, a pole never.
    uses: u32,
    /// How many times each face holding it walks it.
    walks: [u32; 2],
}

impl Holders {
    fn held(&self) -> &[u32] {
        &self.faces[..usize::from(self.count.min(2))]
    }
}

/// One solid as the setting aside reads it: its faces, the faces holding
/// each edge, and its box.
pub(crate) struct Solid {
    faces: Vec<Read>,
    holders: HashMap<EdgeKey, Holders>,
    pub(crate) bound: Aabb,
}

impl Solid {
    /// Whether every edge but a pole is walked an even number of times.
    fn closed(&self) -> bool {
        self.holders
            .values()
            .all(|h| h.uses == u32::MAX || h.uses % 2 == 0)
    }
}

/// The faces of each operand set aside, or `None` where neither sets any
/// aside.
///
/// A solid sets none aside where a face the fuse gathers stands under a
/// placement: a gathered face is rebuilt on its own surface around the
/// edges it shares with the faces set aside, which only a face standing
/// where its node puts it can be. Nor does a solid wholly clear of the
/// other, which the general fuse settles as apart, nor one with an edge
/// more than two faces hold.
pub(crate) fn set_aside(
    model: &Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<[Aside; 2]>> {
    let read_a = read_solid(model, a, tol)?;
    let read_b = read_solid(model, b, tol)?;
    let (box_a, box_b) = (read_a.bound, read_b.bound);
    let sides = [
        side(model, read_a, &box_b, tol)?,
        side(model, read_b, &box_a, tol)?,
    ];
    if sides.iter().all(Aside::is_empty) {
        return Ok(None);
    }
    Ok(Some(sides))
}

/// One solid read: its faces, the faces holding each edge, and its box.
pub(crate) fn read_solid(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<Solid> {
    let found = explore_unique(model, solid, ShapeType::Face)?;
    let mut faces = Vec::with_capacity(found.len());
    // A closed solid has about one and a half edges per face.
    let mut holders: HashMap<EdgeKey, Holders> = HashMap::with_capacity(found.len() * 2);
    let mut whole = Aabb::EMPTY;
    for (index, face) in found.into_iter().enumerate() {
        let Some(face_node) = model.node(&face) else {
            ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = face_node.data() else {
            ogeom_core::ogeom_bail!(Construction, "face node holds no face data");
        };
        let surface = model.geometry().surface(data.surface);
        // A box found on a mesh may fall short of the surface by the
        // mesh's chord, which the margin allows for.
        let exact = matches!(
            surface,
            Some(
                SurfaceGeometry::Plane(_)
                    | SurfaceGeometry::Cylinder(_)
                    | SurfaceGeometry::Cone(_)
                    | SurfaceGeometry::Sphere(_)
                    | SurfaceGeometry::Torus(_)
            )
        );
        let plane = matches!(surface, Some(SurfaceGeometry::Plane(_)));
        let mut margin = data.tolerance.get();
        if !exact {
            margin += ogeom_mesh::Deflection::default().chord;
        }
        let index = u32::try_from(index).unwrap_or(u32::MAX);
        // A plane with holes keeps its wires as walked, for the holes it
        // may leave out.
        let mut rings = (plane && face_node.children().len() > 1).then(Rings::default);
        for wire in face_node.children() {
            let wire = placed_wire(&face, wire);
            let (k0, v0) = rings
                .as_ref()
                .map_or((0, 0), |r| (r.keys.len(), r.ends.len()));
            for_each_wire_edge(model, &wire, |edge, node| {
                let key = EdgeKey::of(edge);
                if let Some(rings) = &mut rings {
                    rings.keys.push(key);
                    rings.ends.extend(node.children().iter().map(Shape::node));
                }
                // Each edge and its ends read the first time a face walks it.
                let held = holders.entry(key).or_insert_with(|| {
                    let mut own = 0.0_f64;
                    let mut degenerate = false;
                    if let Some(e) = node.data().as_edge() {
                        own = e.tolerance.get();
                        degenerate = e.degenerate;
                    }
                    let ends = node.children();
                    let mut loosest_end = 0.0_f64;
                    for vertex in ends {
                        if let Some(v) = model.node(vertex).and_then(|n| n.data().as_vertex()) {
                            loosest_end = loosest_end.max(v.tolerance.get());
                        }
                    }
                    let reach = own.max(loosest_end);
                    let closed = match ends {
                        [_] => true,
                        [first, last] => first.node() == last.node(),
                        _ => false,
                    };
                    Holders {
                        faces: [index, u32::MAX],
                        count: 0,
                        closed,
                        margin: reach,
                        doubt: (own * 2.0).max(loosest_end + own),
                        uses: if degenerate { u32::MAX } else { 0 },
                        walks: [0, 0],
                    }
                });
                margin = margin.max(held.margin);
                // A pole stands at the most, which it keeps.
                held.uses = held.uses.saturating_add(1);
                // A seam is walked twice by its one face.
                if held.count == 0 || held.held().last() != Some(&index) {
                    if held.count < 2 {
                        held.faces[usize::from(held.count)] = index;
                    }
                    held.count = held.count.saturating_add(1);
                }
                if let Some(slot) = held.held().iter().position(|&f| f == index) {
                    held.walks[slot] += 1;
                }
            })?;
            if let Some(rings) = &mut rings {
                rings
                    .spans
                    .push(((k0, rings.keys.len()), (v0, rings.ends.len())));
            }
        }
        let bound = match model.face_bounds(&face) {
            Some(kept) => kept,
            None => ogeom_algo::face_bounds(model, &face)?,
        }
        .expanded(margin * 2.0 + tol.confusion() * 1e3);
        whole = whole.union(&bound);
        // A face under any placement is taken as moved.
        let unmoved = face.location().is_identity();
        faces.push(Read {
            face,
            bound,
            unmoved,
            plane,
            rings,
        });
    }
    Ok(Solid {
        faces,
        holders,
        bound: whole,
    })
}

/// A wire of a gathered face: its run of edges, its run of vertices, and
/// where it is a hole left out, the faces set aside beside it and its
/// edges' doubt.
type Span = (
    (usize, usize),
    (usize, usize),
    Option<(SmallVec<[u32; 2]>, f64)>,
);

/// What one solid sets aside against the other's box `other`.
pub(crate) fn side(
    model: &Model,
    solid: Solid,
    other: &Aabb,
    tol: Tolerances,
) -> OgeomResult<Aside> {
    let closed = solid.closed();
    let Solid { faces, holders, .. } = solid;
    let mut clear: Vec<bool> = faces.iter().map(|r| !r.bound.intersects(other)).collect();
    if !clear.iter().any(|c| *c) || clear.iter().all(|c| *c) {
        return Ok(Aside::default());
    }
    if faces.iter().zip(&clear).any(|(r, c)| !c && !r.unmoved)
        || holders.values().any(|h| h.count > 2)
    {
        return Ok(Aside::default());
    }
    let index = |j: &u32| *j as usize;
    let mut holes: HashMap<usize, Vec<Hole>> = HashMap::new();
    let mut shared: HashSet<EdgeKey> = HashSet::new();
    loop {
        let mut gathered_again = Vec::new();
        shared.clear();
        holes.clear();
        for (i, r) in faces.iter().enumerate() {
            if clear[i] {
                continue;
            }
            let Some(node) = model.node(&r.face) else {
                ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
            };
            let forward = ogeom_topo::Orientation::Forward;
            let wires = node.children();
            // Each wire's edges and the vertices they end on, as read.
            let walked;
            let Rings {
                keys,
                ends,
                spans: runs,
            } = match &r.rings {
                Some(rings) => rings,
                None => {
                    walked = Rings::of(model, &r.face)?;
                    &walked
                }
            };
            let mut spans: Vec<Span> = Vec::with_capacity(runs.len());
            // The ring stored first is the outer one, which a file need not
            // honour: it is where its box is the face's, and the others are
            // then holes, inside it.
            let outer_first = r.rings.is_some()
                && runs.first().is_some_and(|&((k0, k1), (v0, v1))| {
                    let face = model.face_bounds(&r.face);
                    let tolerance = model
                        .node(&r.face)
                        .and_then(|n| n.data().tolerance())
                        .map_or(0.0, |t| t.get());
                    let ring = ring_box(model, &keys[k0..k1], &ends[v0..v1], tol);
                    face.zip(ring).is_some_and(|(face, ring)| {
                        same_box(&ring, &face, tol.confusion() * 1e3 + tolerance * 2.0)
                    })
                });
            for (w, &((k0, k1), (v0, v1))) in runs.iter().enumerate() {
                // A hole of a plane whose every edge a face set aside holds
                // lies within those faces' boxes, which miss the other
                // solid's together.
                let mut beside: SmallVec<[u32; 2]> = SmallVec::new();
                let mut bound = Aabb::EMPTY;
                let mut doubt = 0.0_f64;
                let mut hole = r.plane && w > 0;
                for key in &keys[k0..k1] {
                    if !hole {
                        break;
                    }
                    let found = holders.get(key);
                    doubt = doubt.max(found.map_or(0.0, |h| h.doubt));
                    let held = found.map_or(&[][..], |h| h.held());
                    hole = false;
                    for &j in held {
                        if clear[index(&j)] {
                            hole = true;
                            bound = bound.union(&faces[index(&j)].bound);
                            if !beside.contains(&j) {
                                beside.push(j);
                            }
                        }
                    }
                }
                let hole =
                    (outer_first && hole && !bound.intersects(other)).then_some((beside, doubt));
                spans.push(((k0, k1), (v0, v1), hole));
            }
            // A hole touching a ring the arrangement takes is taken with it.
            if spans.iter().any(|(_, _, hole)| hole.is_some()) {
                let kept: HashSet<TShapeId> = spans
                    .iter()
                    .filter(|(_, _, hole)| hole.is_none())
                    .flat_map(|(_, (v0, v1), _)| ends[*v0..*v1].iter().copied())
                    .collect();
                for span in &mut spans {
                    if span.2.is_some() && ends[span.1.0..span.1.1].iter().any(|v| kept.contains(v))
                    {
                        span.2 = None;
                    }
                }
            }
            let mut arranged: Vec<&[EdgeKey]> = Vec::new();
            let mut left = Vec::new();
            for (at, (wire, ((k0, k1), _, hole))) in wires.iter().zip(spans).enumerate() {
                match hole {
                    Some((beside, doubt)) => left.push(Hole {
                        wire: if r.face.location().is_identity() {
                            wire.clone()
                        } else {
                            wire.beneath(r.face.location(), forward)
                        },
                        at,
                        beside,
                        doubt,
                    }),
                    None => arranged.push(&keys[k0..k1]),
                }
            }
            if !left.is_empty() {
                holes.insert(i, left);
            }
            // The edges the arrangement takes that a face set aside holds
            // are shared with it; a closed one gathers that face after all.
            for edges in arranged {
                for key in edges {
                    if let Some(h) = holders.get(key)
                        && h.held().iter().any(|&j| clear[index(&j)])
                    {
                        if h.closed {
                            gathered_again.extend(h.held().iter().map(index).filter(|&j| clear[j]));
                        }
                        shared.insert(*key);
                    }
                }
            }
        }
        if gathered_again.is_empty() {
            break;
        }
        for j in gathered_again {
            if !faces[j].unmoved {
                return Ok(Aside::default());
            }
            clear[j] = false;
        }
        if !clear.iter().any(|c| *c) {
            return Ok(Aside::default());
        }
    }
    let mut aside = Aside {
        closed,
        ..Aside::default()
    };
    aside.beside = std::sync::Arc::new(Beside {
        holders,
        aside: clear.clone(),
        shared: shared.iter().copied().collect(),
    });
    aside.edges = shared;
    for (i, r) in faces.into_iter().enumerate() {
        if clear[i] {
            aside.faces.push((r.face, i));
            continue;
        }
        if let Some(found) = holes.remove(&i) {
            aside.holes.insert(SameKey(r.face.clone()), found.into());
        }
        aside.gathered.insert(SameKey(r.face), i);
    }
    Ok(aside)
}

/// The box of a ring of straight and circular edges where their nodes put
/// them: its vertices, and each arc's points leading along an axis; `None`
/// for a ring with another kind of edge or one under a placement.
fn ring_box(model: &Model, keys: &[EdgeKey], ends: &[TShapeId], tol: Tolerances) -> Option<Aabb> {
    let mut bound = Aabb::EMPTY;
    for vertex in ends {
        bound = bound.with_point(model.node_by_id(*vertex)?.data().as_vertex()?.point);
    }
    for key in keys {
        if key.placement != 0 {
            return None;
        }
        let data = model.node_by_id(key.node)?.data().as_edge()?;
        let ogeom_topo::EdgeRepr::Curve3d {
            curve,
            location,
            range,
        } = data.curve3d()?
        else {
            return None;
        };
        if !location.is_identity() {
            return None;
        }
        let geometry = model.geometry().curve(*curve)?;
        match geometry {
            ogeom_geom::Curve::Line(_) => {}
            ogeom_geom::Curve::Circle(c) => {
                // An arc leads along an axis at the angle the circle's frame
                // names for it, or half a turn on, where the arc reaches it.
                let frame = c.circle().frame();
                let (x, y) = (frame.x().vector(), frame.y().vector());
                let tau = core::f64::consts::TAU;
                for axis in [
                    ogeom_math::Vector::X,
                    ogeom_math::Vector::Y,
                    ogeom_math::Vector::Z,
                ] {
                    let theta = axis.dot(y).atan2(axis.dot(x));
                    for angle in [theta, theta + core::f64::consts::PI] {
                        let t = ((range.0 - angle) / tau).ceil().mul_add(tau, angle);
                        if t <= range.1 {
                            bound = bound.with_point(geometry.point_at(t, tol).ok()?);
                        }
                    }
                }
            }
            _ => return None,
        }
    }
    (!bound.is_empty()).then_some(bound)
}

/// Whether two boxes have their sides within `room` of each other.
fn same_box(a: &Aabb, b: &Aabb, room: f64) -> bool {
    let (Some(al), Some(ah), Some(bl), Some(bh)) = (a.low(), a.high(), b.low(), b.high()) else {
        return false;
    };
    al.distance(bl) <= room * 2.0 && ah.distance(bh) <= room * 2.0
}

/// `each` of the edges of `wire`, placed as the wire places them, with its
/// node. The edge's orientation is not read here, and is left as the wire
/// stores it where the wire stands where its node puts it.
fn for_each_wire_edge(
    model: &Model,
    wire: &Shape,
    mut each: impl FnMut(&Shape, &TShape),
) -> OgeomResult<()> {
    let Some(node) = model.node(wire) else {
        ogeom_core::ogeom_bail!(Dangling, "wire is not in this model");
    };
    let unplaced = wire.location().is_identity();
    for edge in node.children() {
        let Some(edge_node) = model.node(edge) else {
            ogeom_core::ogeom_bail!(Dangling, "edge is not in this model");
        };
        if unplaced {
            each(edge, edge_node);
        } else {
            each(
                &edge.beneath(wire.location(), wire.orientation()),
                edge_node,
            );
        }
    }
    Ok(())
}

/// A wire of `face` placed as the face places it, its orientation not read:
/// the wire as the face stores it where the face stands unplaced.
fn placed_wire(face: &Shape, wire: &Shape) -> Shape {
    if face.location().is_identity() {
        wire.clone()
    } else {
        wire.beneath(face.location(), face.orientation())
    }
}

/// A face of the result, as [`shells_of`] reads it.
pub(crate) struct Member {
    pub(crate) face: Shape,
    /// For a face set aside, its solid and its place there.
    pub(crate) aside: Option<(usize, usize)>,
    /// For a face the boolean made, how many of its wires to walk: a face
    /// taking holes back holds them after its own wires.
    pub(crate) walk: usize,
    /// The faces set aside, by solid and place, holding the edges of the
    /// holes it took back, which it is joined to without walking them.
    pub(crate) joins: Vec<(usize, u32)>,
}

/// The faces given, grouped into shells by the edges they share, and
/// whether each shell is closed: every edge of it but a pole walked an even
/// number of times, as [`ogeom_algo::is_shell_closed`] counts.
///
/// A face set aside is joined to the others through what `beside` read of
/// its solid, and a hole taken back through the faces it names; only the
/// faces the boolean made are walked, their holes left out.
pub(crate) fn shells_of(
    model: &mut Model,
    faces: &[Member],
    beside: [&Beside; 2],
) -> OgeomResult<(Vec<Shape>, Vec<bool>)> {
    let mut group: Vec<usize> = (0..faces.len()).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    fn join(group: &mut [usize], a: usize, b: usize) {
        let (p, q) = (root(group, a), root(group, b));
        if p != q {
            group[p.max(q)] = p.min(q);
        }
    }
    // Each face set aside where it stands among the faces given.
    let mut at: [Vec<usize>; 2] = [
        vec![usize::MAX; beside[0].aside.len()],
        vec![usize::MAX; beside[1].aside.len()],
    ];
    for (i, member) in faces.iter().enumerate() {
        if let Some((side, place)) = member.aside
            && let Some(slot) = at[side].get_mut(place)
        {
            *slot = i;
        }
    }
    let given = |side: usize, j: u32| {
        at[side]
            .get(j as usize)
            .copied()
            .filter(|&i| i != usize::MAX)
    };
    let any_given = [0, 1].map(|side| {
        faces
            .iter()
            .any(|m| m.aside.is_some_and(|(s, _)| s == side))
    });
    for (side, read) in beside.iter().enumerate() {
        if !any_given[side] {
            continue;
        }
        for h in read.holders.values() {
            if let [p, q] = h.held()
                && let (Some(p), Some(q)) = (given(side, *p), given(side, *q))
            {
                join(&mut group, p, q);
            }
        }
    }
    // The faces made, walked: each edge's first face and its walks.
    let mut made: HashMap<EdgeKey, (usize, u32)> = HashMap::new();
    let mut edges = Vec::new();
    for (i, member) in faces.iter().enumerate() {
        if member.aside.is_some() {
            continue;
        }
        for &(side, j) in &member.joins {
            if let Some(other) = given(side, j) {
                join(&mut group, other, i);
            }
        }
        edges.clear();
        let face = &member.face;
        let Some(node) = model.node(face) else {
            ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
        };
        for wire in node.children().iter().take(member.walk) {
            let wire = placed_wire(face, wire);
            for_each_wire_edge(model, &wire, |edge, node| {
                edges.push((
                    EdgeKey::of(edge),
                    !node.data().as_edge().is_some_and(|e| e.degenerate),
                ));
            })?;
        }
        for (key, counted) in edges.drain(..) {
            let entry = made.entry(key).or_insert((i, 0));
            entry.1 += u32::from(counted);
            let first = entry.0;
            join(&mut group, first, i);
            for (side, read) in beside.iter().enumerate() {
                if let Some(h) = read.holders.get(&key) {
                    for j in h.held() {
                        if let Some(other) = given(side, *j) {
                            join(&mut group, other, i);
                        }
                    }
                }
            }
        }
    }
    // The walks of an edge by the faces set aside among those given.
    let aside_walks = |side: usize, h: &Holders| -> (u32, Option<usize>) {
        let mut walks = 0;
        let mut any = None;
        for (slot, j) in h.held().iter().enumerate() {
            if let Some(i) = given(side, *j) {
                walks += h.walks[slot];
                any = Some(i);
            }
        }
        (walks, any)
    };
    let mut open: Vec<usize> = Vec::new();
    for (key, (first, walks)) in &made {
        let mut total = *walks;
        let mut pole = false;
        for (side, read) in beside.iter().enumerate() {
            if let Some(h) = read.holders.get(key) {
                pole |= h.uses == u32::MAX;
                total += aside_walks(side, h).0;
            }
        }
        if !pole && total % 2 == 1 {
            open.push(*first);
        }
    }
    // Two faces set aside meet as they met in their solid, and the edge of
    // a hole taken back is walked once by the face taking it and once by
    // the face set aside beside it. An edge shared with a face the boolean
    // gathered is walked by a face made, or the shell is open there.
    for (side, read) in beside.iter().enumerate() {
        if !any_given[side] {
            continue;
        }
        for key in &read.shared {
            if made.contains_key(key) {
                continue;
            }
            if let Some(h) = read.holders.get(key)
                && h.uses != u32::MAX
                && let (walks, Some(i)) = aside_walks(side, h)
                && walks % 2 == 1
            {
                open.push(i);
            }
        }
    }
    let mut members: Vec<Vec<Shape>> = Vec::new();
    let mut slot_of: HashMap<usize, usize> = HashMap::new();
    for (i, Member { face, .. }) in faces.iter().enumerate() {
        let r = root(&mut group, i);
        let slot = *slot_of.entry(r).or_insert_with(|| {
            members.push(Vec::new());
            members.len() - 1
        });
        members[slot].push(face.clone());
    }
    let mut closed = vec![true; members.len()];
    for i in open {
        closed[slot_of[&root(&mut group, i)]] = false;
    }
    let shells = members
        .iter()
        .map(|faces| model.add_shell(faces))
        .collect::<OgeomResult<Vec<_>>>()?;
    Ok((shells, closed))
}
