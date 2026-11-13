//! P0: the solid read once, as given, and the request's vertices classed.
//!
//! Every element is keyed by its node in the solid, never by where the
//! request names it, and every later phase walks the elements in key
//! order: a result is the same for every order of the request.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::{Id, id};
use crate::setback::point_of;

/// One face of the solid: its outward plane, its corners in turn about its
/// outward normal, and each hole's corners in turn the other way.
pub(crate) struct Face {
    pub shape: Shape,
    pub normal: Vector,
    pub corners: Vec<Id>,
    pub holes: Vec<Vec<Id>>,
    /// The face's own plane, where the face stands unplaced on it.
    pub plane: Option<Plane>,
}

impl Face {
    /// The face's rings: the outer one, then each hole's.
    pub fn loops(&self) -> impl Iterator<Item = &[Id]> {
        std::iter::once(self.corners.as_slice()).chain(self.holes.iter().map(Vec::as_slice))
    }

    /// The corners before and after `v` on the ring holding it.
    pub fn beside(&self, v: Id) -> Option<(Id, Id)> {
        self.loops().find_map(|ring| {
            let k = ring.iter().position(|c| *c == v)?;
            let n = ring.len();
            Some((ring[(k + n - 1) % n], ring[(k + 1) % n]))
        })
    }
}

/// The solid read as planes bounded by straight edges: the domain of the
/// polyhedral classes (chamfers, M1; rounds between planes, M2).
pub(crate) struct Planar {
    pub faces: Vec<Face>,
    /// Each edge by its two vertices, the lower first.
    pub edges: FastMap<(Id, Id), Id>,
    pub points: FastMap<Id, Point>,
    pub shapes: FastMap<Id, Shape>,
}

/// One edge as a blend sees it: its two faces, each walking it its own
/// way round, its direction, and whether it is convex.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Seat {
    /// The edge.
    pub edge: Id,
    /// The face walking the edge from `ends.0` to `ends.1`, then the other.
    pub faces: [usize; 2],
    /// The edge's ends, in the first face's turn.
    pub ends: (Id, Id),
    /// Outward unit normals of the two faces.
    pub normals: [Vector; 2],
    /// Unit direction from `ends.0` to `ends.1`.
    pub along: Vector,
    /// Material inside the dihedral: a blend there takes material away.
    pub convex: bool,
    /// The dihedral's opening inside the material, in radians: a quarter
    /// turn at a box's edge, three quarters at an inside corner's.
    pub opening: f64,
}

/// A vertex's class: how many edges meet there, how many the request
/// blends, and of those how many are convex and how many concave. It picks
/// the vertex method in P3 and decides no geometry itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct VertexClass {
    pub edges: usize,
    pub blended: usize,
    pub convex: usize,
    pub concave: usize,
}

/// Each edge's uses by the faces walking it: the face and the edge's ends
/// in its turn.
type EdgeUses = FastMap<(Id, Id), Vec<(usize, Id, Id)>>;

impl Planar {
    /// Every edge's seat, keyed by its two vertices, the lower first.
    pub fn seats(&self) -> Option<FastMap<(Id, Id), Seat>> {
        let mut uses: EdgeUses = FastMap::default();
        for (fi, face) in self.faces.iter().enumerate() {
            for ring in face.loops() {
                let n = ring.len();
                for k in 0..n {
                    let (a, b) = (ring[k], ring[(k + 1) % n]);
                    uses.entry((a.min(b), a.max(b)))
                        .or_default()
                        .push((fi, a, b));
                }
            }
        }
        let mut out = FastMap::default();
        for (key, used) in uses {
            let [(f1, a, b), (f2, a2, b2)] = used.as_slice() else {
                return None;
            };
            if (*a2, *b2) != (*b, *a) {
                return None;
            }
            let along = self.points[b] - self.points[a];
            let along = along / along.magnitude();
            let (n1, n2) = (self.faces[*f1].normal, self.faces[*f2].normal);
            let convex = n1.cross(n2).dot(along) > 0.0;
            let between = n1.dot(n2).clamp(-1.0, 1.0).acos();
            let opening = if convex {
                std::f64::consts::PI - between
            } else {
                std::f64::consts::PI + between
            };
            out.insert(
                key,
                Seat {
                    edge: *self.edges.get(&key)?,
                    faces: [*f1, *f2],
                    ends: (*a, *b),
                    normals: [n1, n2],
                    along,
                    convex,
                    opening,
                },
            );
        }
        Some(out)
    }

    /// Every vertex's class under the blended edges `blended`.
    pub fn classes(
        &self,
        seats: &FastMap<(Id, Id), Seat>,
        blended: &ogeom_core::FastSet<Id>,
    ) -> FastMap<Id, VertexClass> {
        let mut out: FastMap<Id, VertexClass> = FastMap::default();
        for ((a, b), seat) in seats {
            for v in [a, b] {
                let class = out.entry(*v).or_insert(VertexClass {
                    edges: 0,
                    blended: 0,
                    convex: 0,
                    concave: 0,
                });
                class.edges += 1;
                if blended.contains(&seat.edge) {
                    class.blended += 1;
                    if seat.convex {
                        class.convex += 1;
                    } else {
                        class.concave += 1;
                    }
                }
            }
        }
        out
    }
}

/// `solid` read as planes bounded by one ring of straight edges each, or
/// `None` where it is not one.
pub(crate) fn read_planar(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Planar>> {
    if model.kind_of(solid)? != ShapeType::Solid
        || explore_unique(model, solid, ShapeType::Shell)?.len() != 1
    {
        return Ok(None);
    }
    let mut points: FastMap<Id, Point> = FastMap::default();
    let mut shapes: FastMap<Id, Shape> = FastMap::default();
    for v in explore_unique(model, solid, ShapeType::Vertex)? {
        points.insert(id(&v), point_of(model, &v)?);
        shapes.insert(id(&v), v);
    }
    let mut edges: FastMap<(Id, Id), Id> = FastMap::default();
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
            return Ok(None);
        };
        let (a, b) = (id(&a), id(&b));
        if a == b {
            return Ok(None);
        }
        let Some((curve, _)) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
            return Ok(None);
        };
        if !matches!(curve, Curve::Line(_)) || edges.insert((a.min(b), a.max(b)), id(&e)).is_some()
        {
            return Ok(None);
        }
        shapes.insert(id(&e), e);
    }
    let mut faces = Vec::new();
    for shape in explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&shape).and_then(|n| n.data().as_face()) else {
            return Ok(None);
        };
        let Some(SurfaceGeometry::Plane(surface)) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        let unplaced = shape.location().is_identity() && data.location.is_identity();
        let plane = unplaced.then(|| surface.plane());
        let (_, normal) = ogeom_algo::face_normal(model, &shape, tol)?;
        // Each wire's edges, joined into one cycle by their shared vertices;
        // the outer ring the one enclosing the most, turned about the
        // outward normal, the holes the other way.
        let mut loops: Vec<(f64, Vec<Id>)> = Vec::new();
        for wire in model.children_of(&shape)? {
            let mut pairs: Vec<(Id, Id)> = Vec::new();
            for e in explore_unique(model, &wire, ShapeType::Edge)? {
                let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
                    return Ok(None);
                };
                pairs.push((id(&a), id(&b)));
            }
            let Some(corners) = cycle(&pairs) else {
                return Ok(None);
            };
            let mut newell = Vector::ZERO;
            for k in 0..corners.len() {
                let p = points[&corners[k]].to_vector();
                let q = points[&corners[(k + 1) % corners.len()]].to_vector();
                newell += p.cross(q);
            }
            loops.push((newell.dot(normal), corners));
        }
        loops.sort_by(|a, b| b.0.abs().total_cmp(&a.0.abs()));
        let Some((turn, mut corners)) = loops.first().cloned() else {
            return Ok(None);
        };
        if turn < 0.0 {
            corners.reverse();
        }
        let mut holes = Vec::with_capacity(loops.len() - 1);
        for (turn, mut hole) in loops.into_iter().skip(1) {
            if turn > 0.0 {
                hole.reverse();
            }
            holes.push(hole);
        }
        // The class is exact planes: a corner standing off its face's plane
        // (held there only by its vertex's tolerance) is outside it.
        let at = points[&corners[0]];
        if corners
            .iter()
            .chain(holes.iter().flatten())
            .any(|c| (points[c] - at).dot(normal).abs() > tol.confusion())
        {
            return Ok(None);
        }
        faces.push(Face {
            shape,
            normal,
            corners,
            holes,
            plane,
        });
    }
    Ok(Some(Planar {
        faces,
        edges,
        points,
        shapes,
    }))
}

/// The vertices of `pairs` joined into one cycle, each pair an edge's two
/// ends in either order, or `None` where they do not make exactly one.
fn cycle(pairs: &[(Id, Id)]) -> Option<Vec<Id>> {
    let mut next: FastMap<Id, Vec<Id>> = FastMap::default();
    for (a, b) in pairs {
        next.entry(*a).or_default().push(*b);
        next.entry(*b).or_default().push(*a);
    }
    if next.values().any(|n| n.len() != 2) {
        return None;
    }
    let start = pairs.first()?.0;
    let mut out = vec![start];
    let mut prev = start;
    let mut here = next[&start][0];
    while here != start {
        out.push(here);
        let n = &next[&here];
        let step = if n[0] == prev { n[1] } else { n[0] };
        prev = here;
        here = step;
        if out.len() > pairs.len() {
            return None;
        }
    }
    (out.len() == pairs.len()).then_some(out)
}

/// The faces round vertex `v` in turn, each followed by the face across its
/// edge arriving at `v`, or `None` where they do not close round it once.
pub(crate) fn faces_round(faces: &[Face], v: Id) -> Option<Vec<usize>> {
    let holding: Vec<usize> = (0..faces.len())
        .filter(|fi| faces[*fi].beside(v).is_some())
        .collect();
    let first = *holding.first()?;
    let mut order = vec![first];
    loop {
        let (prev, _) = faces[*order.last()?].beside(v)?;
        // The face walking the edge from `v` to `prev`.
        let next = holding
            .iter()
            .copied()
            .find(|g| faces[*g].beside(v).is_some_and(|(_, after)| after == prev))?;
        if next == first {
            break;
        }
        if order.contains(&next) {
            return None;
        }
        order.push(next);
    }
    (order.len() == holding.len()).then_some(order)
}
