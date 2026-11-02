//! Booleans between closed triangle meshes, decided exactly.
//!
//! [`boolean`] takes two closed, outward-wound, 2-manifold meshes that do
//! not cross themselves and answers their union, intersection or
//! difference as one more such mesh. [`boolean_with_shape`] meshes a solid
//! first and does the same.
//!
//! # Exact decisions, rounded new vertices
//!
//! Every decision is the exact sign of a polynomial in the input
//! coordinates: which side of a plane a point lies on, whether two
//! triangles meet, the order of points along a line, whether a point lies
//! inside a triangle, the winding number of a mesh about a point. A point
//! the operands make, where an edge of one crosses a triangle of the other
//! or two coplanar edges cross, is never rounded while decisions are made
//! on it: it is kept as the input coordinates that define it, and every
//! test on it is evaluated in interval arithmetic first and exactly when
//! the interval cannot tell. So the topology of the answer is the topology
//! of the exact boolean of the two polyhedra.
//!
//! Only at the end are the new vertices rounded to `f64`, each to the
//! coordinates nearest its exact position within a few units in the last
//! place. [`MeshBoolean::rounding`] bounds how far any moved. The input
//! vertices kept are not moved at all. Rounding leaves the connectivity as
//! decided; where two operands' faces cross at a vanishing angle it can
//! tilt a sliver piece, a few units in the last place across, past its
//! neighbour.
//!
//! # How
//!
//! 1. Every pair of triangles whose boxes touch is met exactly: the points
//!    and segments they share.
//! 2. Points found more than once under different names, as where an edge
//!    passes exactly through a vertex, are merged by exact comparison.
//! 3. Each triangle that met the other operand is triangulated afresh with
//!    those points and segments as constraints, the points on its edges
//!    shared with its neighbours, so each cut operand is again closed.
//! 4. The pieces are grouped into regions that no segment separates, and
//!    each region is classified once: on the other operand, where a piece
//!    lies in a coplanar triangle of it, and otherwise inside or outside it
//!    by an exact winding number.
//! 5. The operation keeps the regions it selects, and the pieces are welded
//!    at the shared points.
//!
//! # Coplanar faces
//!
//! Where the operands share a plane, the shared area is cut out of both
//! along the other's edges, and each piece there is *on* the other
//! operand, its normal the same as the other's or opposite. One copy of a
//! shared face with the same normal is kept by a union or an intersection,
//! and none by a difference; a shared face with opposite normals, where
//! the solids touch from outside, is dropped by a union and an
//! intersection and kept, from the first operand, by a difference. These
//! are the boundaries of the regularised booleans.
//!
//! # Refused
//!
//! An operand that is open, not 2-manifold, has a degenerate triangle,
//! crosses itself, or is wound inside out, by name. Operands that do not
//! meet, neither touching nor one holding the other. A result that would
//! not be a 2-manifold, where the operands touch along an edge or at a
//! point without overlapping there.
//!
//! # References
//!
//! The implicit points and the filtered predicates over them follow
//! M. Attene, *Indirect Predicates for Geometric Constructions*, and the
//! pipeline of pairs, per-triangle cuts and region labels follows
//! G. Cherchi, M. Livesu, R. Scateni and M. Attene, *Fast and Robust Mesh
//! Arrangements using Floating-point Arithmetic*, and G. Cherchi,
//! F. Pellacini, M. Attene and M. Livesu, *Interactive and Robust Mesh
//! Booleans*. Rays through edges are resolved as in H. Edelsbrunner and
//! E. P. Mücke, *Simulation of Simplicity*; points are found by the
//! stochastic walk of O. Devillers, S. Pion and M. Teillaud, *Walking in a
//! Triangulation*.

mod classify;
mod exact;
mod pairs;
mod points;
mod split;

use ogeom_core::{
    Exact, FastMap, FastSet, OgeomResult, Predicates, Sign, Tolerances, ogeom_bail, parallel,
};
use ogeom_math::{Aabb, Point};
use ogeom_topo::{Model, Shape, Triangulation};

use self::pairs::{Key, Side, Touch, Tri};
use self::points::{Def, Mid, P3, Pt, by_weight, normal_approx, normal_sign, orient2d};
use crate::{Bvh, Deflection};

/// A boolean operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BoolKind {
    /// Everything in either operand.
    Union,
    /// Everything in both operands.
    Intersection,
    /// Everything in the first operand and not in the second.
    Difference,
}

/// The answer of a mesh boolean.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshBoolean {
    /// The result: closed, 2-manifold and wound outward, or empty where
    /// nothing is left. Its vertices are the input vertices kept, first
    /// operand first, each at its input position exactly and in input
    /// order, then the new ones from [`MeshBoolean::first_new_vertex`]
    /// on. It carries no normals and no parameters.
    pub mesh: Triangulation,
    /// For each triangle of the result, whether it is new: a piece of an
    /// input triangle that the other operand cut. A triangle not flagged is
    /// an input triangle, kept whole, possibly wound the other way.
    pub cut: Vec<bool>,
    /// The index of the first vertex made where the operands meet.
    pub first_new_vertex: usize,
    /// The farthest any new vertex lies from the exact point it rounds:
    /// a bound, of the order of the last place of its coordinates.
    pub rounding: f64,
}

/// The boolean `op` of two closed meshes.
///
/// Each operand must be closed and 2-manifold: every edge has exactly two
/// triangles running it opposite ways, and the triangles about each vertex
/// form one fan. It must have no triangle whose corners are collinear, must
/// not cross or touch itself, and must be wound so that its triangles face
/// out of the solid it bounds, which may have cavities.
///
/// Every decision is exact; only the new vertices are rounded, see the
/// [module documentation](self).
///
/// # Errors
///
/// [`OgeomError::NullObject`](ogeom_core::OgeomError::NullObject) if an
/// operand has no triangles.
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle names
/// a vertex its mesh does not have, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a corner is
/// not finite.
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if an
/// operand is open, is not 2-manifold, has a degenerate triangle, crosses
/// itself or is wound inside out, if the operands do not meet, or if the
/// result would not be 2-manifold.
pub fn boolean(a: &Triangulation, b: &Triangulation, op: BoolKind) -> OgeomResult<MeshBoolean> {
    let [first, second]: [OgeomResult<Operand<'_>>; 2] = parallel::map_ordered(
        &[(a, 0, "first"), (b, 1, "second")],
        |_, &(m, side, name)| Operand::new(m, side, name),
    )
    .try_into()
    .map_err(|_| ogeom_core::ogeom_err!(Invariant, "two operands made other than two"))?;
    let ops = [first?, second?];
    for o in &ops {
        o.check_embedded()?;
    }
    for o in &ops {
        o.check_wound()?;
    }
    let touches = meet_all(&ops[0], &ops[1]);
    let mut table = Table::new(&ops);
    let mut marks: [Vec<Mark>; 2] = [
        vec![Mark::default(); ops[0].tris.len()],
        vec![Mark::default(); ops[1].tris.len()],
    ];
    for (ta, tb, touch) in &touches {
        table.take(&mut marks, *ta, *tb, touch);
    }
    table.merge(&ops);
    let flat: [Vec<usize>; 2] = [ops[0].drops()?, ops[1].drops()?];
    for side in 0..2 {
        table.settle(&ops[side], &flat[side], &mut marks[side]);
    }
    let mut fixed: FastSet<(u32, u32)> = FastSet::default();
    let cuts = [
        table.cut(&ops[0], &flat[0], &marks[0], &mut fixed)?,
        table.cut(&ops[1], &flat[1], &marks[1], &mut fixed)?,
    ];
    let mut labels: [Vec<(Vec<usize>, Label)>; 2] = [Vec::new(), Vec::new()];
    for side in 0..2 {
        let regions = regions(&ops[side], &cuts[side], &fixed)?;
        labels[side] = table.label(
            &ops,
            side,
            &flat[side],
            &marks[side],
            &cuts[side].pieces,
            regions,
        )?;
    }
    let meets = !touches.is_empty()
        || labels
            .iter()
            .flatten()
            .any(|(_, label)| *label != Label::Out);
    if !meets {
        ogeom_bail!(
            Construction,
            "the operands do not meet: neither touches or holds the other"
        );
    }
    assemble(
        &table,
        &cuts,
        &labels,
        op,
        a.deflection_met && b.deflection_met,
    )
}

/// The boolean `op` of a closed mesh and a solid, the solid meshed first
/// to `chord`: [`triangulate`](crate::triangulate()) with
/// [`Deflection::with_chord`], then [`boolean`] with the mesh as the
/// first operand.
///
/// # Errors
///
/// As [`Deflection::with_chord`], [`triangulate`](crate::triangulate()) and
/// [`boolean`]; the solid's mesh is the second operand in [`boolean`]'s
/// messages.
pub fn boolean_with_shape(
    mesh: &Triangulation,
    model: &Model,
    shape: &Shape,
    op: BoolKind,
    chord: f64,
    tol: Tolerances,
) -> OgeomResult<MeshBoolean> {
    let solid = crate::triangulate(model, shape, Deflection::with_chord(chord)?, tol)?;
    boolean(mesh, &solid, op)
}

/// One operand, checked and indexed.
struct Operand<'a> {
    side: Side,
    name: &'static str,
    mesh: &'a Triangulation,
    tris: Vec<Tri>,
    bvh: Bvh,
    /// For each triangle, the triangle across each side, side `k` running
    /// from corner `k` to the next.
    across: Vec<[u32; 3]>,
}

fn coords(p: Point) -> P3 {
    // `+ 0.0` turns a negative zero positive, so equal points have equal
    // bits.
    [p.x + 0.0, p.y + 0.0, p.z + 0.0]
}

fn id(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(u32::MAX)
}

impl<'a> Operand<'a> {
    fn new(mesh: &'a Triangulation, side: Side, name: &'static str) -> OgeomResult<Self> {
        if mesh.triangles.is_empty() {
            ogeom_bail!(NullObject, "the {name} mesh has no triangles");
        }
        if u32::try_from(mesh.positions.len()).is_err() {
            ogeom_bail!(Range, "the {name} mesh has more than four billion vertices");
        }
        let bvh = Bvh::build(mesh)?;
        let mut tris = Vec::with_capacity(mesh.triangles.len());
        let mut edges: FastMap<(u32, u32), u32> = FastMap::default();
        edges.reserve(3 * mesh.triangles.len());
        for (t, &v) in mesh.triangles.iter().enumerate() {
            let p = v.map(|k| coords(mesh.positions[k as usize]));
            let flat = (0..3).all(|k| normal_sign(&p, k).is_zero());
            if flat {
                ogeom_bail!(
                    Construction,
                    "triangle {t} of the {name} mesh is degenerate: its corners are collinear"
                );
            }
            for k in 0..3 {
                if edges.insert((v[k], v[(k + 1) % 3]), id(t)).is_some() {
                    ogeom_bail!(
                        Construction,
                        "the {name} mesh is not 2-manifold: edge {}-{} is run the same way twice",
                        v[k],
                        v[(k + 1) % 3]
                    );
                }
            }
            tris.push(Tri {
                side,
                id: id(t),
                v,
                p,
            });
        }
        let chunks: Vec<std::ops::Range<usize>> = (0..tris.len())
            .step_by(CHUNK)
            .map(|s| s..(s + CHUNK).min(tris.len()))
            .collect();
        let across: Vec<Vec<Option<[u32; 3]>>> = parallel::map_ordered(&chunks, |_, range| {
            range
                .clone()
                .map(|t| {
                    let v: [u32; 3] = tris[t].v;
                    let mut out = [0; 3];
                    for k in 0..3 {
                        out[k] = *edges.get(&(v[(k + 1) % 3], v[k]))?;
                    }
                    Some(out)
                })
                .collect()
        });
        let across: Vec<[u32; 3]> = across
            .into_iter()
            .flatten()
            .collect::<Option<_>>()
            .map_or_else(
                || {
                    let (lo, hi) = edges
                        .keys()
                        .filter(|&&(u, w)| !edges.contains_key(&(w, u)))
                        .map(|&(u, w)| (u.min(w), u.max(w)))
                        .min()
                        .unwrap_or((0, 0));
                    Err(ogeom_core::ogeom_err!(
                        Construction,
                        "the {name} mesh is open: edge {lo}-{hi} has one triangle"
                    ))
                },
                Ok,
            )?;
        drop(edges);
        let operand = Self {
            side,
            name,
            mesh,
            tris,
            bvh,
            across,
        };
        operand.check_fans()?;
        Ok(operand)
    }

    /// Every vertex's triangles form one fan.
    fn check_fans(&self) -> OgeomResult<()> {
        let mut corners = vec![0_u32; self.mesh.positions.len()];
        let mut first = vec![u32::MAX; self.mesh.positions.len()];
        for (t, tri) in self.tris.iter().enumerate() {
            for &v in &tri.v {
                corners[v as usize] += 1;
                if first[v as usize] == u32::MAX {
                    first[v as usize] = id(t);
                }
            }
        }
        for (v, &start) in first.iter().enumerate() {
            if start == u32::MAX {
                continue;
            }
            let v = id(v);
            let mut at = start;
            let mut count = 0;
            loop {
                count += 1;
                let t = self.tris[at as usize].v;
                let k = t.iter().position(|&x| x == v).unwrap_or(0);
                // The next triangle about `v` is across the side into `v`.
                let next = self.across[at as usize][(k + 2) % 3];
                if next == start || count > corners[v as usize] {
                    break;
                }
                at = next;
            }
            if count != corners[v as usize] {
                ogeom_bail!(
                    Construction,
                    "the {} mesh is not 2-manifold at vertex {v}: its triangles form more than one fan",
                    self.name
                );
            }
        }
        Ok(())
    }

    fn bounds(t: &Tri) -> (P3, P3) {
        let mut lo = t.p[0];
        let mut hi = t.p[0];
        for q in &t.p[1..] {
            for k in 0..3 {
                lo[k] = lo[k].min(q[k]);
                hi[k] = hi[k].max(q[k]);
            }
        }
        (lo, hi)
    }

    /// The triangles of this operand whose boxes touch the box `lo..hi`.
    fn near(&self, lo: P3, hi: P3) -> Vec<u32> {
        self.bvh.triangles_where(|b: &Aabb| {
            let (Some(bl), Some(bh)) = (b.low(), b.high()) else {
                return false;
            };
            let (bl, bh) = (coords(bl), coords(bh));
            (0..3).all(|k| bl[k] <= hi[k] && lo[k] <= bh[k])
        })
    }

    /// No two triangles meet but where they share corners: along their
    /// common edge, or at their common vertex.
    fn check_embedded(&self) -> OgeomResult<()> {
        let chunks: Vec<std::ops::Range<usize>> = (0..self.tris.len())
            .step_by(CHUNK)
            .map(|s| s..(s + CHUNK).min(self.tris.len()))
            .collect();
        let found = parallel::map_ordered(&chunks, |_, range| {
            for t in range.clone() {
                let a = &self.tris[t];
                let (lo, hi) = Self::bounds(a);
                for u in self.near(lo, hi) {
                    if u as usize <= t {
                        continue;
                    }
                    let b = Tri {
                        side: 1 - self.side,
                        ..self.tris[u as usize]
                    };
                    let shared: Vec<u32> =
                        a.v.iter().copied().filter(|v| b.v.contains(v)).collect();
                    let clean = shared.len() < 3
                        && (apart_but_shared(a, &b, &shared)
                            || pairs::meet(a, &b).is_none_or(|touch| {
                                touch.points.iter().all(|(key, _)| {
                                    matches!(key, Key::Vertex(_, v) if shared.contains(v))
                                })
                            }));
                    if !clean {
                        return Some((t, u));
                    }
                }
            }
            None
        });
        if let Some((t, u)) = found.into_iter().flatten().next() {
            ogeom_bail!(
                Construction,
                "the {} mesh crosses itself: triangles {t} and {u} meet away from their shared corners",
                self.name
            );
        }
        Ok(())
    }

    /// Each connected piece faces out of the solid: from the middle of one
    /// of its triangles, a ray leaving the way it faces finds the winding
    /// number of the whole mesh 0, and behind it 1.
    fn check_wound(&self) -> OgeomResult<()> {
        let mut parent: Vec<u32> = (0..id(self.tris.len())).collect();
        for t in 0..self.tris.len() {
            for &u in &self.across[t] {
                union(&mut parent, id(t), u);
            }
        }
        let mut seen = FastSet::default();
        for t in 0..self.tris.len() {
            if !seen.insert(find(&mut parent, id(t))) {
                continue;
            }
            let tri = &self.tris[t];
            let axis = by_weight(normal_approx(&tri.p))
                .into_iter()
                .find(|&k| !normal_sign(&tri.p, k).is_zero())
                .unwrap_or(0);
            let corner = tri.p.map(|q| Pt::new(Def::Vertex(q)));
            let mid = Mid([&corner[0], &corner[1], &corner[2]]);
            let out = normal_sign(&tri.p, axis);
            let w = classify::winding(&self.tris, &self.bvh, &mid, axis, out, Some(id(t)))?;
            if w != 0 {
                ogeom_bail!(
                    Construction,
                    "the {} mesh is not wound outward: triangle {t} faces into the solid",
                    self.name
                );
            }
        }
        Ok(())
    }

    /// For each triangle, a coordinate its plane is not parallel to, the
    /// one its normal is largest along.
    fn drops(&self) -> OgeomResult<Vec<usize>> {
        self.tris
            .iter()
            .map(|t| {
                by_weight(normal_approx(&t.p))
                    .into_iter()
                    .find(|&k| !normal_sign(&t.p, k).is_zero())
                    .ok_or_else(|| ogeom_core::ogeom_err!(Invariant, "a triangle has no normal"))
            })
            .collect()
    }
}

/// Whether two triangles of one mesh sharing the corners `shared` meet
/// only there, decided from the planes alone where that is enough.
///
/// Sharing an edge and not coplanar, they meet in the line their planes
/// meet in, which each crosses only along that edge. Sharing a corner, they
/// meet only there when the other two corners of one lie strictly on one
/// side of the other's plane.
fn apart_but_shared(a: &Tri, b: &Tri, shared: &[u32]) -> bool {
    let others = |t: &Tri| -> Vec<P3> {
        (0..3)
            .filter(|&k| !shared.contains(&t.v[k]))
            .map(|k| t.p[k])
            .collect()
    };
    let against = |t: &Tri, q: P3| Exact::orient3d(t.p[0], t.p[1], t.p[2], q);
    let apart = match shared.len() {
        2 => others(b).first().is_some_and(|&w| !against(a, w).is_zero()),
        1 => [(a, b), (b, a)].iter().any(|(x, y)| {
            let o = others(y);
            let (s, t) = (against(x, o[0]), against(x, o[1]));
            !s.is_zero() && s == t
        }),
        _ => false,
    };
    apart || (others(b).iter().all(|&q| against(a, q).is_zero()) && flat_apart(a, b, shared))
}

/// As [`apart_but_shared`], for two coplanar triangles: some line through
/// a shared corner along an edge has one triangle strictly on each side,
/// but for the shared corners.
fn flat_apart(a: &Tri, b: &Tri, shared: &[u32]) -> bool {
    let Some(drop) = by_weight(normal_approx(&a.p))
        .into_iter()
        .find(|&k| !normal_sign(&a.p, k).is_zero())
    else {
        return false;
    };
    let (i, j) = points::axes(drop);
    let flat = |q: P3| [q[i], q[j]];
    let others = |t: &Tri| -> Vec<[f64; 2]> {
        (0..3)
            .filter(|&k| !shared.contains(&t.v[k]))
            .map(|k| flat(t.p[k]))
            .collect()
    };
    let at = |t: &Tri, v: u32| t.v.iter().position(|&x| x == v).map(|k| flat(t.p[k]));
    let splits = |p: [f64; 2], q: [f64; 2], mine: &[[f64; 2]], theirs: &[[f64; 2]]| {
        let side = |r: [f64; 2]| Exact::orient2d(p, q, r);
        let Some(&first) = mine.first() else {
            return false;
        };
        let s = side(first);
        !s.is_zero()
            && mine.iter().all(|&r| side(r) == s)
            && theirs.iter().all(|&r| side(r) == s.reversed())
    };
    match shared {
        [u, w] => match (at(a, *u), at(a, *w)) {
            (Some(p), Some(q)) => splits(p, q, &others(a), &others(b)),
            _ => false,
        },
        [v] => [(a, b), (b, a)].iter().any(|(x, y)| {
            let Some(p) = at(x, *v) else {
                return false;
            };
            let mine = others(x);
            (0..2).any(|k| splits(p, mine[k], &mine[1 - k..2 - k], &others(y)))
        }),
        _ => false,
    }
}

/// How many triangles a parallel task takes.
const CHUNK: usize = 256;

fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let up = parent[parent[x as usize] as usize];
        parent[x as usize] = up;
        x = up;
    }
    x
}

/// Join two sets, the lower index the root.
fn union(parent: &mut [u32], a: u32, b: u32) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra < rb {
        parent[rb as usize] = ra;
    } else if rb < ra {
        parent[ra as usize] = rb;
    }
}

/// Every pair of triangles, one from each operand, that meet, and what
/// they share, in the order of the first operand's triangles.
fn meet_all(a: &Operand<'_>, b: &Operand<'_>) -> Vec<(u32, u32, Touch)> {
    let chunks: Vec<std::ops::Range<usize>> = (0..a.tris.len())
        .step_by(CHUNK)
        .map(|s| s..(s + CHUNK).min(a.tris.len()))
        .collect();
    parallel::map_ordered(&chunks, |_, range| {
        let mut out = Vec::new();
        for t in range.clone() {
            let ta = &a.tris[t];
            let (lo, hi) = Operand::bounds(ta);
            for u in b.near(lo, hi) {
                if let Some(touch) = pairs::meet(ta, &b.tris[u as usize]) {
                    out.push((id(t), u, touch));
                }
            }
        }
        out
    })
    .into_iter()
    .flatten()
    .collect()
}

/// What the other operand left on one triangle.
#[derive(Debug, Clone, Default)]
struct Mark {
    /// Points on the triangle, as indices into the table.
    points: Vec<u32>,
    segments: Vec<[u32; 2]>,
    /// The other operand's triangles coplanar with this one that it meets.
    coplanar: Vec<u32>,
}

/// Every point, by index: the first operand's vertices, the second's, then
/// the points made where they meet.
struct Table {
    pts: Vec<Pt>,
    names: FastMap<Key, u32>,
    /// The name of each point made, in order after the input vertices.
    keys: Vec<Key>,
    /// The vertex counts of the two operands.
    counts: [usize; 2],
    /// The index each point is merged into: the lowest of the points equal
    /// to it.
    rep: Vec<u32>,
}

impl Table {
    fn new(ops: &[Operand<'_>; 2]) -> Self {
        let mut pts = Vec::with_capacity(ops[0].mesh.positions.len() + ops[1].mesh.positions.len());
        for o in ops {
            pts.extend(
                o.mesh
                    .positions
                    .iter()
                    .map(|&p| Pt::new(Def::Vertex(coords(p)))),
            );
        }
        Self {
            pts,
            names: FastMap::default(),
            keys: Vec::new(),
            counts: [ops[0].mesh.positions.len(), ops[1].mesh.positions.len()],
            rep: Vec::new(),
        }
    }

    fn index(&mut self, key: Key, pt: &Pt) -> u32 {
        match key {
            Key::Vertex(0, v) => v,
            Key::Vertex(_, v) => id(self.counts[0]) + v,
            _ => {
                if let Some(&at) = self.names.get(&key) {
                    return at;
                }
                let at = id(self.pts.len());
                self.pts.push(Pt::new(pt.def));
                self.names.insert(key, at);
                self.keys.push(key);
                at
            }
        }
    }

    fn take(&mut self, marks: &mut [Vec<Mark>; 2], ta: u32, tb: u32, touch: &Touch) {
        let ids: Vec<u32> = touch
            .points
            .iter()
            .map(|(key, pt)| self.index(*key, pt))
            .collect();
        for (side, t, other) in [(0, ta, tb), (1, tb, ta)] {
            let m = &mut marks[side][t as usize];
            m.points.extend(&ids);
            m.segments
                .extend(touch.segments.iter().map(|s| [ids[s[0]], ids[s[1]]]));
            if touch.coplanar {
                m.coplanar.push(other);
            }
        }
    }

    /// Merge every point with the lowest-indexed point equal to it.
    fn merge(&mut self, ops: &[Operand<'_>; 2]) {
        let n = self.pts.len();
        let mut parent: Vec<u32> = (0..id(n)).collect();
        let explicit = self.counts[0] + self.counts[1];
        let mut used = vec![false; explicit];
        for (side, o) in ops.iter().enumerate() {
            let base = if side == 0 { 0 } else { self.counts[0] };
            for t in &o.tris {
                for &v in &t.v {
                    used[base + v as usize] = true;
                }
            }
        }
        // Input vertices, equal when their coordinates are.
        let bits = |p: &P3| p.map(f64::to_bits);
        let mut at: FastMap<[u64; 3], u32> = FastMap::default();
        let mut sorted: Vec<(P3, u32)> = Vec::new();
        for (i, _) in used.iter().enumerate().filter(|(_, u)| **u) {
            let Def::Vertex(p) = self.pts[i].def else {
                continue;
            };
            match at.get(&bits(&p)) {
                Some(&j) => union(&mut parent, j, id(i)),
                None => {
                    at.insert(bits(&p), id(i));
                    sorted.push((p, id(i)));
                }
            }
        }
        sorted.sort_by(|x, y| {
            x.0[0]
                .total_cmp(&y.0[0])
                .then(x.0[1].total_cmp(&y.0[1]))
                .then(x.0[2].total_cmp(&y.0[2]))
                .then(x.1.cmp(&y.1))
        });
        // Made points against the input vertices inside their boxes.
        type Bounds = Option<[(f64, f64); 3]>;
        let made: Vec<(u32, Bounds)> = (explicit..n)
            .map(|i| {
                let b = self.pts[i].bounds().map(|b| b.map(|v| (v.lo, v.hi)));
                (id(i), b)
            })
            .collect();
        let pairs_found: Vec<Vec<(u32, u32)>> = parallel::map_ordered(&made, |_, &(i, b)| {
            let mut out = Vec::new();
            let pt = &self.pts[i as usize];
            match b {
                Some(b) => {
                    let start = sorted.partition_point(|e| e.0[0] < b[0].0);
                    for e in &sorted[start..] {
                        if e.0[0] > b[0].1 {
                            break;
                        }
                        let inside = (1..3).all(|k| b[k].0 <= e.0[k] && e.0[k] <= b[k].1);
                        if inside && points::same(pt, &self.pts[e.1 as usize]) {
                            out.push((e.1, i));
                        }
                    }
                }
                None => {
                    for e in &sorted {
                        if points::same(pt, &self.pts[e.1 as usize]) {
                            out.push((e.1, i));
                        }
                    }
                }
            }
            out
        });
        for (x, y) in pairs_found.into_iter().flatten() {
            union(&mut parent, x, y);
        }
        // Made points against each other, swept along x.
        let mut boxed: Vec<(u32, [(f64, f64); 3])> = Vec::new();
        let mut loose: Vec<u32> = Vec::new();
        for &(i, b) in &made {
            match b {
                Some(b) => boxed.push((i, b)),
                None => loose.push(i),
            }
        }
        boxed.sort_by(|x, y| x.1[0].0.total_cmp(&y.1[0].0).then(x.0.cmp(&y.0)));
        let mut pairs_found = Vec::new();
        for (k, (i, b)) in boxed.iter().enumerate() {
            for (j, c) in &boxed[k + 1..] {
                if c[0].0 > b[0].1 {
                    break;
                }
                let overlap = (1..3).all(|d| c[d].0 <= b[d].1 && b[d].0 <= c[d].1);
                if overlap && points::same(&self.pts[*i as usize], &self.pts[*j as usize]) {
                    pairs_found.push((*i, *j));
                }
            }
        }
        for &i in &loose {
            for &(j, _) in &made {
                if j != i && points::same(&self.pts[i as usize], &self.pts[j as usize]) {
                    pairs_found.push((i, j));
                }
            }
        }
        for (x, y) in pairs_found {
            union(&mut parent, x, y);
        }
        self.rep = (0..id(n)).map(|i| find(&mut parent, i)).collect();
    }

    /// Rename every mark's points to their merged index, drop repeats and
    /// segments merged to a point, and share the points on each edge with
    /// the triangle across it.
    fn settle(&self, o: &Operand<'_>, flat: &[usize], marks: &mut [Mark]) {
        let base = if o.side == 0 { 0 } else { self.counts[0] };
        let corners = |t: &Tri| t.v.map(|v| self.rep[base + v as usize]);
        type Settled = (Vec<u32>, Vec<[u32; 2]>, Vec<(u32, u32)>);
        let which: Vec<usize> = (0..marks.len())
            .filter(|&t| !(marks[t].points.is_empty() && marks[t].segments.is_empty()))
            .collect();
        let settled: Vec<Settled> = parallel::map_ordered(&which, |_, &t| {
            let m = &marks[t];
            let tri = &o.tris[t];
            let c = corners(tri);
            let mut raw: Vec<(u32, u32)> = m
                .points
                .iter()
                .map(|&x| (self.rep[x as usize], x))
                .filter(|(r, _)| !c.contains(r))
                .collect();
            raw.sort_unstable();
            raw.dedup();
            let mut points = Vec::new();
            let mut shares = Vec::new();
            for group in raw.chunk_by(|x, y| x.0 == y.0) {
                let r = group[0].0;
                points.push(r);
                let on = group
                    .iter()
                    .find_map(|&(_, x)| self.edge_of(o.side, tri, x))
                    .or_else(|| {
                        (0..3).find(|&k| {
                            orient2d(
                                flat[t],
                                &self.pts[c[k] as usize],
                                &self.pts[c[(k + 1) % 3] as usize],
                                &self.pts[r as usize],
                            )
                            .is_zero()
                        })
                    });
                if let Some(k) = on {
                    shares.push((o.across[t][k], r));
                }
            }
            let mut segments: Vec<[u32; 2]> = m
                .segments
                .iter()
                .map(|s| {
                    let (a, b) = (self.rep[s[0] as usize], self.rep[s[1] as usize]);
                    if a < b { [a, b] } else { [b, a] }
                })
                .filter(|s| s[0] != s[1])
                .collect();
            segments.sort_unstable();
            segments.dedup();
            (points, segments, shares)
        });
        let mut shares = Vec::new();
        for (&t, (points, segments, share)) in which.iter().zip(settled) {
            let m = &mut marks[t];
            m.points = points;
            m.segments = segments;
            m.coplanar.sort_unstable();
            m.coplanar.dedup();
            shares.extend(share);
        }
        for (u, x) in shares {
            let m = &mut marks[u as usize];
            if let Err(at) = m.points.binary_search(&x) {
                m.points.insert(at, x);
            }
        }
    }

    /// The side of `tri` that point `x` lies on by its name: where an edge
    /// of `tri` crosses the other operand.
    fn edge_of(&self, side: Side, tri: &Tri, x: u32) -> Option<usize> {
        let made = (x as usize).checked_sub(self.counts[0] + self.counts[1])?;
        let edge = match self.keys[made] {
            Key::Lpi { side: s, edge, .. } if s == side => edge,
            Key::Ssi { a, b } => {
                if side == 0 {
                    a
                } else {
                    b
                }
            }
            _ => return None,
        };
        (0..3).find(|&k| {
            let (u, w) = (tri.v[k], tri.v[(k + 1) % 3]);
            [u.min(w), u.max(w)] == edge
        })
    }

    /// The pieces of every triangle of one operand: whole where nothing
    /// marked it, cut where something did. The constraint edges of the cuts
    /// are added to `fixed`.
    fn cut(
        &self,
        o: &Operand<'_>,
        flat: &[usize],
        marks: &[Mark],
        fixed: &mut FastSet<(u32, u32)>,
    ) -> OgeomResult<Cut> {
        let marked: Vec<bool> = marks
            .iter()
            .map(|m| !(m.points.is_empty() && m.segments.is_empty()))
            .collect();
        let which: Vec<usize> = (0..marks.len()).filter(|&t| marked[t]).collect();
        let cuts = parallel::map_ordered(&which, |_, &t| {
            split::split(
                &self.pts,
                flat[t],
                self.corners(o, t),
                &marks[t].points,
                &marks[t].segments,
            )
        });
        let mut done = which.iter().zip(cuts).peekable();
        let mut pieces = Vec::with_capacity(o.tris.len());
        let mut start = Vec::with_capacity(o.tris.len() + 1);
        for t in 0..o.tris.len() {
            start.push(pieces.len());
            if done.peek().is_some_and(|(u, _)| **u == t) {
                let Some((_, c)) = done.next() else {
                    continue;
                };
                let c = c?;
                fixed.extend(c.fixed.iter().map(|&[a, b]| (a, b)));
                let cut = c.triangles.len() > 1;
                pieces.extend(c.triangles.into_iter().map(|v| Piece {
                    v,
                    parent: id(t),
                    cut,
                }));
            } else {
                pieces.push(Piece {
                    v: self.corners(o, t),
                    parent: id(t),
                    cut: false,
                });
            }
        }
        start.push(pieces.len());
        Ok(Cut {
            pieces,
            start,
            marked,
        })
    }

    /// The corners of triangle `t` of operand `o`, as merged points.
    fn corners(&self, o: &Operand<'_>, t: usize) -> [u32; 3] {
        let base = if o.side == 0 { 0 } else { self.counts[0] };
        o.tris[t].v.map(|v| self.rep[base + v as usize])
    }

    /// Classify each region of one operand's pieces against the other
    /// operand, from one piece of it.
    fn label(
        &self,
        ops: &[Operand<'_>; 2],
        side: usize,
        flat: &[usize],
        marks: &[Mark],
        pieces: &[Piece],
        regions: Vec<Vec<usize>>,
    ) -> OgeomResult<Vec<(Vec<usize>, Label)>> {
        let (own, other) = (&ops[side], &ops[1 - side]);
        let labels = parallel::map_ordered(&regions, |_, region| -> OgeomResult<Label> {
            let piece = &pieces[region[0]];
            let t = piece.parent as usize;
            let drop = flat[t];
            let corner = piece.v.map(|v| &self.pts[v as usize]);
            let mid = Mid(corner);
            let parent = &own.tris[t];
            for &u in &marks[t].coplanar {
                let q = &other.tris[u as usize];
                let qc = q.p.map(|p| Pt::new(Def::Vertex(p)));
                let facing = normal_sign(&q.p, drop);
                let o = [0, 1, 2].map(|k| orient2d(drop, &qc[k], &qc[(k + 1) % 3], &mid));
                // A piece lies wholly inside or wholly outside each coplanar
                // triangle, whose edges cut it; its centroid on an edge would
                // mean a cut was missed.
                if o.contains(&facing.reversed()) {
                    continue;
                }
                if o.contains(&Sign::Zero) {
                    ogeom_bail!(Invariant, "a piece straddles a coplanar triangle's edge");
                }
                return Ok(if facing == normal_sign(&parent.p, drop) {
                    Label::OnSame
                } else {
                    Label::OnOpposite
                });
            }
            let w = classify::winding(&other.tris, &other.bvh, &mid, drop, Sign::Positive, None)?;
            match w {
                0 => Ok(Label::Out),
                1 => Ok(Label::In),
                _ => ogeom_bail!(
                    Invariant,
                    "the {} mesh winds {w} times about a point",
                    other.name
                ),
            }
        });
        regions
            .into_iter()
            .zip(labels)
            .map(|(r, l)| Ok((r, l?)))
            .collect()
    }
}

/// A triangle of a cut operand, by point indices.
#[derive(Debug, Clone, Copy)]
struct Piece {
    v: [u32; 3],
    parent: u32,
    cut: bool,
}

/// Where a region of one operand lies against the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Label {
    In,
    Out,
    /// On a face of the other, both facing the same way.
    OnSame,
    /// On a face of the other, facing opposite ways.
    OnOpposite,
}

/// The pieces of one cut operand.
struct Cut {
    pieces: Vec<Piece>,
    /// Where each triangle's pieces start, and the end of the last.
    start: Vec<usize>,
    /// Whether each triangle was triangulated afresh.
    marked: Vec<bool>,
}

/// The pieces of one cut operand grouped into regions: joined across every
/// edge that is not a constraint, each region's pieces in order and the
/// regions in the order of their first piece.
///
/// A triangle left whole is joined to its neighbours left whole through
/// the operand's own edges: no constraint runs along an edge with a whole
/// triangle on both sides, since a constraint lies on a triangle that is
/// cut. The pieces of cut triangles, and the whole triangles beside them,
/// are joined through their edges as merged points.
fn regions(
    o: &Operand<'_>,
    cut: &Cut,
    fixed: &FastSet<(u32, u32)>,
) -> OgeomResult<Vec<Vec<usize>>> {
    let pieces = &cut.pieces;
    let near_cut = |t: usize| o.across[t].iter().any(|&u| cut.marked[u as usize]);
    let mut edges: FastMap<(u32, u32), u32> = FastMap::default();
    for t in 0..o.tris.len() {
        if !(cut.marked[t] || near_cut(t)) {
            continue;
        }
        for (i, p) in pieces
            .iter()
            .enumerate()
            .take(cut.start[t + 1])
            .skip(cut.start[t])
        {
            for k in 0..3 {
                if edges.insert((p.v[k], p.v[(k + 1) % 3]), id(i)).is_some() {
                    ogeom_bail!(Invariant, "the cut {} mesh runs an edge twice", o.name);
                }
            }
        }
    }
    let mut parent: Vec<u32> = (0..id(pieces.len())).collect();
    for t in 0..o.tris.len() {
        for (i, p) in pieces
            .iter()
            .enumerate()
            .take(cut.start[t + 1])
            .skip(cut.start[t])
        {
            for k in 0..3 {
                if !cut.marked[t] {
                    let u = o.across[t][k];
                    if !cut.marked[u as usize] {
                        union(&mut parent, id(i), id(cut.start[u as usize]));
                        continue;
                    }
                }
                let (a, b) = (p.v[k], p.v[(k + 1) % 3]);
                let Some(&j) = edges.get(&(b, a)) else {
                    ogeom_bail!(Invariant, "the cut {} mesh is open", o.name);
                };
                if !fixed.contains(&(a.min(b), a.max(b))) {
                    union(&mut parent, id(i), j);
                }
            }
        }
    }
    let mut slot = vec![usize::MAX; pieces.len()];
    let mut out: Vec<Vec<usize>> = Vec::new();
    for i in 0..pieces.len() {
        let r = find(&mut parent, id(i)) as usize;
        if slot[r] == usize::MAX {
            slot[r] = out.len();
            out.push(Vec::new());
        }
        out[slot[r]].push(i);
    }
    Ok(out)
}

/// Whether `op` keeps a region of operand `side` labelled `label`, and
/// whether it turns it over.
const fn keep(op: BoolKind, side: usize, label: Label) -> Option<bool> {
    let kept = matches!(
        (op, side, label),
        (BoolKind::Union, _, Label::Out)
            | (BoolKind::Union | BoolKind::Intersection, 0, Label::OnSame)
            | (BoolKind::Intersection, _, Label::In)
            | (BoolKind::Difference, 0, Label::Out | Label::OnOpposite)
            | (BoolKind::Difference, 1, Label::In)
    );
    if kept {
        Some(matches!(op, BoolKind::Difference) && side == 1)
    } else {
        None
    }
}

fn assemble(
    table: &Table,
    cuts: &[Cut; 2],
    labels: &[Vec<(Vec<usize>, Label)>; 2],
    op: BoolKind,
    deflection_met: bool,
) -> OgeomResult<MeshBoolean> {
    let mut kept: Vec<([u32; 3], bool, (usize, usize))> = Vec::new();
    for side in 0..2 {
        for (region, label) in &labels[side] {
            let Some(flip) = keep(op, side, *label) else {
                continue;
            };
            for &i in region {
                let p = &cuts[side].pieces[i];
                let v = if flip { [p.v[0], p.v[2], p.v[1]] } else { p.v };
                kept.push((v, p.cut, (side, i)));
            }
        }
    }
    kept.sort_unstable_by_key(|k| k.2);
    // Where the result could fail to be 2-manifold: at the points of cut
    // triangles, and where the operands share a point. Elsewhere it is an
    // operand's own surface, kept or dropped a whole region at a time.
    let mut focus = vec![false; table.pts.len()];
    for cut in cuts {
        for p in &cut.pieces {
            if cut.marked[p.parent as usize] {
                for &v in &p.v {
                    focus[v as usize] = true;
                }
            }
        }
    }
    for (i, &r) in table.rep.iter().enumerate() {
        if r as usize != i {
            focus[r as usize] = true;
        }
    }
    let mut slot = vec![u32::MAX; table.pts.len()];
    for k in &kept {
        for &v in &k.0 {
            slot[v as usize] = 0;
        }
    }
    let explicit = table.counts[0] + table.counts[1];
    let mut mesh = Triangulation::new();
    let mut rounding = 0.0_f64;
    let mut first_new_vertex = 0;
    let mut near = Vec::new();
    for (v, s) in slot.iter_mut().enumerate() {
        if *s == u32::MAX {
            continue;
        }
        *s = id(mesh.positions.len());
        if v < explicit {
            first_new_vertex += 1;
        }
        near.push(focus[v]);
        let (p, err) = table.pts[v].rounded();
        rounding = rounding.max(err);
        mesh.positions.push(Point::new(p[0], p[1], p[2]));
    }
    mesh.triangles = kept.iter().map(|k| k.0.map(|v| slot[v as usize])).collect();
    mesh.deflection_met = deflection_met;
    check_result(&mesh, &near)?;
    Ok(MeshBoolean {
        cut: kept.iter().map(|k| k.1).collect(),
        mesh,
        first_new_vertex,
        rounding,
    })
}

/// The result is closed and 2-manifold about the vertices in `focus`: each
/// edge at one of them run once each way, one fan about each.
fn check_result(mesh: &Triangulation, focus: &[bool]) -> OgeomResult<()> {
    let at = |v: u32| focus[v as usize];
    let near: Vec<usize> = (0..mesh.triangles.len())
        .filter(|&t| mesh.triangles[t].iter().any(|&v| at(v)))
        .collect();
    let mut edges: FastMap<(u32, u32), u32> = FastMap::default();
    for &t in &near {
        let v = mesh.triangles[t];
        for k in 0..3 {
            if edges.insert((v[k], v[(k + 1) % 3]), id(t)).is_some() {
                ogeom_bail!(
                    Construction,
                    "the result would not be 2-manifold: the operands touch along an edge without overlapping there"
                );
            }
        }
    }
    if edges
        .keys()
        .any(|&(a, b)| (at(a) || at(b)) && !edges.contains_key(&(b, a)))
    {
        ogeom_bail!(Invariant, "the result is open");
    }
    let mut corners = vec![0_u32; mesh.positions.len()];
    for &t in &near {
        for &v in &mesh.triangles[t] {
            corners[v as usize] += 1;
        }
    }
    let mut seen = vec![false; mesh.positions.len()];
    for &t in &near {
        for &v in &mesh.triangles[t] {
            if !at(v) || seen[v as usize] {
                continue;
            }
            seen[v as usize] = true;
            let (start, mut here, mut count) = (id(t), id(t), 0);
            loop {
                count += 1;
                let w = mesh.triangles[here as usize];
                let k = w.iter().position(|&x| x == v).unwrap_or(0);
                match edges.get(&(v, w[(k + 2) % 3])) {
                    Some(&next) if next != start && count <= corners[v as usize] => here = next,
                    _ => break,
                }
            }
            if count != corners[v as usize] {
                ogeom_bail!(
                    Construction,
                    "the result would not be 2-manifold: the operands touch at a point without overlapping there"
                );
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
