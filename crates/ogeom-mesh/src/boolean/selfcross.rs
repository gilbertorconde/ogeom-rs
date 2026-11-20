//! Where one mesh crosses itself, and the mesh made embedded.
//!
//! [`self_intersections`] lists the pairs of triangles that meet anywhere
//! but at the corners they share, decided by the predicate the boolean
//! refuses a self-crossing operand with. [`resolve_self_intersections`]
//! cuts a closed mesh along every such crossing and keeps the boundary of
//! the points it winds about at least once: the regularised union of the
//! mesh with itself.
//!
//! # How the cut differs from a boolean's
//!
//! In a boolean the segments left on a triangle come from the other
//! operand, which does not cross itself, so they never cross each other.
//! Here they come from the triangle's own mesh: where three triangles
//! cross at one point, the segments the other two leave on each cross
//! there. That point is named by the three triangles and held as the
//! meeting of their planes; where a segment runs along an edge of a
//! triangle coplanar with the one cut, the crossing is where that edge
//! meets a plane or another edge, named as the boolean names those. Each
//! triangle finds the crossings of its own segments, so a point gets the
//! same name on every triangle it lies on, and the boolean's exact merge
//! joins the names that turn out to be one point.
//!
//! # Labels
//!
//! Each region of pieces is labelled by the exact winding number of the
//! whole mesh just in front of it and just behind it. The one in front is
//! counted along a ray from the region's first piece, leaving out the
//! triangles that hold that piece; behind it is that plus one for each of
//! those triangles facing the same way as the piece and minus one for each
//! facing the other way. A region is kept where one side winds at least
//! once and the other does not, facing the side that does not. Where
//! coplanar triangles overlap, the lowest-numbered of them gives the piece.

use ogeom_core::{FastMap, FastSet, OgeomResult, Sign, ogeom_bail, parallel};
use ogeom_topo::Triangulation;

use super::pairs::{self, Key, Touch, Tri};
use super::points::{Def, Mid, P3, Pt, normal_sign, orient2d};
use super::{
    CHUNK, Cut, Mark, MeshBoolean, Operand, Piece, Table, classify, coords, find, finish, id,
    meet_at_shared_only, near, union,
};
use crate::Bvh;

/// The answer of [`resolve_self_intersections`].
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// Closed, 2-manifold, wound outward and embedded: the boundary of
    /// every point the input winds about at least once. Empty where it
    /// winds about none. Its vertices are the input vertices kept, in
    /// input order and unmoved, then the new ones; vertices with equal
    /// coordinates are one. It carries no normals and no parameters.
    pub mesh: Triangulation,
    /// The pairs of input triangles that crossed, as
    /// [`self_intersections`] counts them.
    pub crossings: usize,
    /// For each output triangle, whether it is a piece of a cut input
    /// triangle, as [`MeshBoolean::cut`](super::MeshBoolean::cut). A
    /// triangle not flagged is an input triangle kept whole, possibly
    /// wound the other way.
    pub cut: Vec<bool>,
    /// The farthest any new vertex lies from the exact point it rounds, as
    /// [`MeshBoolean::rounding`](super::MeshBoolean::rounding), summed over
    /// the passes where the result was resolved again.
    pub rounding: f64,
}

/// What the mesh's triangles are called in messages.
const NAME: &str = "input";

/// Every pair of triangles of `mesh` that meet anywhere but at the corners
/// they share (their common edge or vertex), decided exactly as the
/// boolean's predicates decide, lower number first, sorted.
///
/// Any triangulation is taken, open or closed, 2-manifold or not.
/// Triangles share a corner by vertex number: two vertices at one point
/// are different corners, so triangles meeting there are listed. A
/// triangle repeated is listed against its copy. Pairs are found through
/// a hierarchy of boxes, then decided exactly.
///
/// # Errors
///
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle
/// names a vertex the mesh does not have,
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a corner is
/// not finite, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// triangle is degenerate, its corners collinear.
pub fn self_intersections(mesh: &Triangulation) -> OgeomResult<Vec<[u32; 2]>> {
    if mesh.triangles.is_empty() {
        return Ok(Vec::new());
    }
    if u32::try_from(mesh.positions.len()).is_err() {
        ogeom_bail!(Range, "the {NAME} mesh has more than four billion vertices");
    }
    let bvh = Bvh::build(mesh)?;
    let tris = (0..mesh.triangles.len())
        .map(|t| Operand::tri(mesh, t, 0, NAME))
        .collect::<OgeomResult<Vec<Tri>>>()?;
    Ok(crossing_pairs(&tris, &bvh))
}

/// A closed mesh that crosses itself made embedded: overlapping shells
/// merged, surfaces inside dropped, as a union of the mesh with itself.
///
/// The mesh is cut along every crossing, and each piece is labelled by the
/// exact winding number of the whole mesh just in front of it and just
/// behind it. A piece is kept where one side winds at least once and the
/// other does not, facing the side that does not: a shell nested in
/// another wound the same way is dropped, one wound the other way is kept
/// as a cavity, a fold keeps what winds positively, and a mesh wound
/// inside out everywhere gives an empty mesh. Coplanar faces facing the
/// same way are kept once, and facing each other, where shells touch, not
/// at all.
///
/// A mesh that does not cross itself and faces outward comes back with
/// the same triangles in the same order and the same vertex positions,
/// less any vertex no triangle uses; normals and parameters are dropped.
///
/// Every decision is exact and only the new vertices are rounded, as in
/// [`boolean`](super::boolean()). Where rounding tilts a sliver piece past
/// its neighbour, the result is checked by [`self_intersections`] and
/// resolved again, at most three more times; a triangle rounding leaves
/// flat is first flipped into the triangle across its longest side, which
/// moves no vertex. A result without rounded vertices is embedded as
/// decided and is not checked again.
///
/// # Errors
///
/// [`OgeomError::NullObject`](ogeom_core::OgeomError::NullObject) if the
/// mesh has no triangles; as [`self_intersections`] for its vertices and
/// triangles; and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// it is open or not 2-manifold, if the result would not be 2-manifold,
/// where parts of the mesh touch along an edge or at a point without
/// overlapping there, or if after rounding the result still crosses itself
/// or keeps a flat triangle that cannot be flipped away.
pub fn resolve_self_intersections(mesh: &Triangulation) -> OgeomResult<Resolved> {
    let (first, crossings) = resolve_once(mesh)?;
    let mut out = Resolved {
        mesh: first.mesh,
        crossings,
        cut: first.cut,
        rounding: first.rounding,
    };
    let mut rounded = out.rounding > 0.0;
    for pass in 0..=PASSES {
        if !rounded {
            break;
        }
        flip_needles(&mut out.mesh, &mut out.cut)?;
        if self_intersections(&out.mesh)?.is_empty() {
            break;
        }
        if pass == PASSES {
            ogeom_bail!(
                Construction,
                "rounding the new vertices to floating point leaves the result crossing itself"
            );
        }
        let (next, _) = resolve_once(&out.mesh).map_err(|e| {
            ogeom_core::ogeom_err!(
                Construction,
                "rounding the new vertices to floating point leaves the result crossing itself, and resolving it again fails: {e}"
            )
        })?;
        out.cut = carried(&out.mesh, &out.cut, &next);
        out.rounding += next.rounding;
        rounded = next.rounding > 0.0;
        out.mesh = next.mesh;
    }
    Ok(out)
}

/// How many times a result that rounding left crossing itself is resolved
/// again.
const PASSES: usize = 3;

/// Remove every triangle of the closed, 2-manifold `mesh` whose corners
/// are collinear: with the triangle across its longest side it makes a
/// quadrilateral whose other diagonal runs from its middle corner, and
/// that diagonal is flipped in. The two new triangles cover the triangle
/// across exactly, so no vertex moves and the surface is the same point
/// set; both are flagged as pieces.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// two corners of such a triangle coincide, or the diagonal is already an
/// edge.
fn flip_needles(mesh: &mut Triangulation, cut: &mut [bool]) -> OgeomResult<()> {
    let corners =
        |m: &Triangulation, t: usize| m.triangles[t].map(|v| coords(m.positions[v as usize]));
    let needle = |m: &Triangulation, t: usize| {
        let p = corners(m, t);
        (0..3).all(|k| normal_sign(&p, k).is_zero())
    };
    let mut todo: Vec<usize> = (0..mesh.triangles.len())
        .filter(|&t| needle(mesh, t))
        .collect();
    if todo.is_empty() {
        return Ok(());
    }
    let mut edges: FastMap<(u32, u32), usize> = FastMap::default();
    for (t, v) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            edges.insert((v[k], v[(k + 1) % 3]), t);
        }
    }
    let mut flips = 0;
    while let Some(t) = todo.pop() {
        if !needle(mesh, t) {
            continue;
        }
        flips += 1;
        if flips > mesh.triangles.len() {
            ogeom_bail!(Invariant, "flipping flat triangles does not end");
        }
        let (v, p) = (mesh.triangles[t], corners(mesh, t));
        // The corner strictly between the other two, along an axis they
        // differ on.
        let middle = (0..3).find(|&k| {
            let (a, b, c) = (p[(k + 1) % 3], p[(k + 2) % 3], p[k]);
            (0..3)
                .find(|&x| a[x] != b[x])
                .is_some_and(|x| (a[x] < c[x] && c[x] < b[x]) || (b[x] < c[x] && c[x] < a[x]))
        });
        let Some(k) = middle else {
            ogeom_bail!(
                Construction,
                "rounding the new vertices to floating point makes two corners of a triangle one"
            );
        };
        let (c, a, b) = (v[k], v[(k + 1) % 3], v[(k + 2) % 3]);
        let Some(&u) = edges.get(&(b, a)) else {
            ogeom_bail!(Invariant, "the resolved mesh is open");
        };
        let w = mesh.triangles[u];
        let d = w.into_iter().find(|&x| x != a && x != b).unwrap_or(a);
        if d == c || edges.contains_key(&(c, d)) || edges.contains_key(&(d, c)) {
            ogeom_bail!(
                Construction,
                "rounding the new vertices to floating point leaves a flat triangle that cannot be removed"
            );
        }
        for (x, y) in [(c, a), (a, b), (b, c), (b, a), (a, d), (d, b)] {
            edges.remove(&(x, y));
        }
        mesh.triangles[t] = [b, c, d];
        mesh.triangles[u] = [c, a, d];
        for (x, y, at) in [
            (b, c, t),
            (c, d, t),
            (d, b, t),
            (c, a, u),
            (a, d, u),
            (d, c, u),
        ] {
            edges.insert((x, y), at);
        }
        cut[t] = true;
        cut[u] = true;
        todo.extend([t, u]);
    }
    Ok(())
}

/// The flags of `next`, a resolution of `mesh` whose triangles carry
/// `cut`: a triangle `next` cut is a piece, and one it kept whole has the
/// flag of the triangle of `mesh` at the same corners.
fn carried(mesh: &Triangulation, cut: &[bool], next: &MeshBoolean) -> Vec<bool> {
    let corners = |m: &Triangulation, t: [u32; 3]| {
        let mut c = t.map(|v| coords(m.positions[v as usize]).map(f64::to_bits));
        c.sort_unstable();
        c
    };
    let before: FastMap<[[u64; 3]; 3], bool> = mesh
        .triangles
        .iter()
        .zip(cut)
        .map(|(&t, &c)| (corners(mesh, t), c))
        .collect();
    next.mesh
        .triangles
        .iter()
        .zip(&next.cut)
        .map(|(&t, &c)| c || before.get(&corners(&next.mesh, t)).copied().unwrap_or(true))
        .collect()
}

/// One resolution of `mesh`, its new vertices rounded once, and how many
/// pairs of its triangles crossed.
fn resolve_once(mesh: &Triangulation) -> OgeomResult<(MeshBoolean, usize)> {
    let o = Operand::new(mesh, 0, NAME)?;
    let crossings = crossing_pairs(&o.tris, &o.bvh);
    let touches = parallel::map_ordered(&crossings, |_, &[t, u]| {
        pairs::meet(&o.tris[t as usize], &o.tris[u as usize]).map(|touch| (t, u, named(touch)))
    });
    let flat = o.drops()?;
    let mut table = Table::new(std::slice::from_ref(&o));
    let mut marks = vec![Mark::default(); o.tris.len()];
    let mut lines: Vec<Vec<Line>> = vec![Vec::new(); o.tris.len()];
    for (t, u, touch) in touches.iter().flatten() {
        take(&o, &flat, &mut table, &mut marks, &mut lines, *t, *u, touch)?;
    }
    let found = parallel::map_ordered(&lines, |t, own| crossings_on(&o, &table, flat[t], t, own));
    for (t, points) in found.into_iter().enumerate() {
        for (key, def) in points {
            let at = table.index(key, &Pt::new(def));
            marks[t].points.push(at);
        }
    }
    table.merge(std::slice::from_ref(&o));
    table.settle(&o, &flat, &mut marks);
    let mut fixed = FastSet::default();
    let cut = table.cut(&o, &flat, &marks, &mut fixed)?;
    let regions = regions(&o, &table, &flat, &cut, &fixed)?;
    let kept = keep(&o, &table, &flat, &marks, &cut, &regions)?;
    let out = finish(
        &table,
        std::slice::from_ref(&cut),
        &kept,
        mesh.deflection_met,
        "parts of the mesh",
    )?;
    Ok((out, crossings.len()))
}

/// Every pair of `tris`, indexed by `bvh`, that meets away from its shared
/// corners, lower number first, sorted.
fn crossing_pairs(tris: &[Tri], bvh: &Bvh) -> Vec<[u32; 2]> {
    let chunks: Vec<std::ops::Range<usize>> = (0..tris.len())
        .step_by(CHUNK)
        .map(|s| s..(s + CHUNK).min(tris.len()))
        .collect();
    parallel::map_ordered(&chunks, |_, range| {
        let mut out = Vec::new();
        for t in range.clone() {
            let a = &tris[t];
            let (lo, hi) = Operand::bounds(a);
            let mut met: Vec<u32> = near(bvh, lo, hi)
                .into_iter()
                .filter(|&u| u as usize > t)
                .filter(|&u| {
                    let b = Tri {
                        side: 1,
                        ..tris[u as usize]
                    };
                    !meet_at_shared_only(a, &b)
                })
                .collect();
            met.sort_unstable();
            out.extend(met.into_iter().map(|u| [id(t), u]));
        }
        out
    })
    .into_iter()
    .flatten()
    .collect()
}

/// `touch` with each crossing of two coplanar edges named with the lower
/// edge first, as two edges of one mesh can be found in either order.
fn named(mut touch: Touch) -> Touch {
    for (key, _) in &mut touch.points {
        if let Key::Ssi { a, b } = *key
            && b < a
        {
            *key = Key::Ssi { a: b, b: a };
        }
    }
    touch
}

/// What a segment left on a triangle runs along.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Along {
    /// The line where the triangle's plane meets triangle `u`'s.
    Plane(u32),
    /// Edge `edge`, lower vertex first, of a triangle coplanar with it.
    Edge([u32; 2], [P3; 2]),
    /// One of the triangle's own edges: a segment there crosses nothing.
    Own,
}

/// A segment on a triangle, as table points, and what it runs along.
#[derive(Debug, Clone, Copy)]
struct Line {
    ends: [u32; 2],
    along: Along,
}

/// Record what triangles `t` and `u` share on both: the points, the
/// segments and what each runs along, and each as coplanar with the other
/// where they are.
#[allow(clippy::too_many_arguments)]
fn take(
    o: &Operand<'_>,
    flat: &[usize],
    table: &mut Table,
    marks: &mut [Mark],
    lines: &mut [Vec<Line>],
    t: u32,
    u: u32,
    touch: &Touch,
) -> OgeomResult<()> {
    let ids: Vec<u32> = touch
        .points
        .iter()
        .map(|(key, pt)| table.index(*key, pt))
        .collect();
    for (host, other) in [(t, u), (u, t)] {
        let m = &mut marks[host as usize];
        m.points.extend(&ids);
        if touch.coplanar {
            m.coplanar.push(other);
        }
        for s in &touch.segments {
            let ends = [ids[s[0]], ids[s[1]]];
            m.segments.push(ends);
            let along = if touch.coplanar {
                let ends = [&touch.points[s[0]].1, &touch.points[s[1]].1];
                let drop = flat[host as usize];
                if edge_under(&o.tris[host as usize], drop, ends).is_some() {
                    Along::Own
                } else if let Some(k) = edge_under(&o.tris[other as usize], drop, ends) {
                    let (edge, line) = edge_of(&o.tris[other as usize], k);
                    Along::Edge(edge, line)
                } else {
                    ogeom_bail!(Invariant, "a coplanar segment runs along no edge");
                }
            } else {
                Along::Plane(other)
            };
            lines[host as usize].push(Line { ends, along });
        }
    }
    Ok(())
}

/// The side of `tri` that both `ends`, points in its plane, lie on the
/// line of, in the projection dropping `drop`.
fn edge_under(tri: &Tri, drop: usize, ends: [&Pt; 2]) -> Option<usize> {
    let c = tri.p.map(|q| Pt::new(Def::Vertex(q)));
    (0..3).find(|&k| {
        ends.iter()
            .all(|e| orient2d(drop, &c[k], &c[(k + 1) % 3], *e).is_zero())
    })
}

/// Side `k` of `tri` as its vertex numbers and corners, lower number first.
fn edge_of(tri: &Tri, k: usize) -> ([u32; 2], [P3; 2]) {
    let (i, j) = (k, (k + 1) % 3);
    if tri.v[i] < tri.v[j] {
        ([tri.v[i], tri.v[j]], [tri.p[i], tri.p[j]])
    } else {
        ([tri.v[j], tri.v[i]], [tri.p[j], tri.p[i]])
    }
}

/// The points where the segments on triangle `t` cross each other, each
/// strictly inside both, named and defined.
fn crossings_on(
    o: &Operand<'_>,
    table: &Table,
    drop: usize,
    t: usize,
    own: &[Line],
) -> Vec<(Key, Def)> {
    let pt = |i: u32| &table.pts[i as usize];
    let strictly = |a: u32, b: u32, c: u32, d: u32| {
        let (s, r) = (
            orient2d(drop, pt(a), pt(b), pt(c)),
            orient2d(drop, pt(a), pt(b), pt(d)),
        );
        !s.is_zero() && s == r.reversed()
    };
    let mut out = Vec::new();
    for (i, x) in own.iter().enumerate() {
        for y in &own[i + 1..] {
            if x.along == Along::Own || y.along == Along::Own || x.along == y.along {
                continue;
            }
            let ([a, b], [c, d]) = (x.ends, y.ends);
            if !(strictly(a, b, c, d) && strictly(c, d, a, b)) {
                continue;
            }
            out.push(match (x.along, y.along) {
                (Along::Plane(u), Along::Plane(w)) => {
                    let mut tris = [id(t), u, w];
                    tris.sort_unstable();
                    (
                        Key::Tpi { tris },
                        Def::Tpi {
                            planes: tris.map(|k| o.tris[k as usize].p),
                        },
                    )
                }
                (Along::Plane(u), Along::Edge(edge, line))
                | (Along::Edge(edge, line), Along::Plane(u)) => (
                    Key::Lpi {
                        side: 0,
                        edge,
                        tri: u,
                    },
                    Def::Lpi {
                        line,
                        plane: o.tris[u as usize].p,
                    },
                ),
                (Along::Edge(e, le), Along::Edge(f, lf)) => {
                    let ((a, first), (b, second)) = if e < f {
                        ((e, le), (f, lf))
                    } else {
                        ((f, lf), (e, le))
                    };
                    (
                        Key::Ssi { a, b },
                        Def::Ssi {
                            first,
                            second,
                            drop,
                        },
                    )
                }
                _ => continue,
            });
        }
    }
    out
}

/// The pieces of the cut mesh grouped into regions, each region's pieces
/// in order and the regions in the order of their first piece.
///
/// Where triangles cross, four pieces meet along one edge, so pieces are
/// joined only to the piece across each of their edges within the same
/// triangle, or, along a side of the triangle, to the piece of the
/// triangle across that side; never across a constraint.
fn regions(
    o: &Operand<'_>,
    table: &Table,
    flat: &[usize],
    cut: &Cut,
    fixed: &FastSet<(u32, u32)>,
) -> OgeomResult<Vec<Vec<usize>>> {
    let pieces = &cut.pieces;
    let near_cut = |t: usize| o.across[t].iter().any(|&u| cut.marked[u as usize]);
    let mut edges: FastMap<(u32, u32, u32), u32> = FastMap::default();
    for t in 0..o.tris.len() {
        if !(cut.marked[t] || near_cut(t)) {
            continue;
        }
        let range = cut.start[t]..cut.start[t + 1];
        for (i, piece) in pieces[range.clone()].iter().enumerate() {
            let (i, v) = (range.start + i, piece.v);
            for k in 0..3 {
                edges.insert((id(t), v[k], v[(k + 1) % 3]), id(i));
            }
        }
    }
    let mut parent: Vec<u32> = (0..id(pieces.len())).collect();
    for t in 0..o.tris.len() {
        let corner = table.corners(o, t).map(|c| &table.pts[c as usize]);
        let range = cut.start[t]..cut.start[t + 1];
        for (i, piece) in pieces[range.clone()].iter().enumerate() {
            let (i, v) = (range.start + i, piece.v);
            for k in 0..3 {
                let (a, b) = (v[k], v[(k + 1) % 3]);
                if !cut.marked[t] {
                    let u = o.across[t][k] as usize;
                    if !cut.marked[u] {
                        union(&mut parent, id(i), id(cut.start[u]));
                        continue;
                    }
                }
                let j = match edges.get(&(id(t), b, a)) {
                    Some(&j) => j,
                    None => {
                        let on = |side: usize, x: u32| {
                            orient2d(
                                flat[t],
                                corner[side],
                                corner[(side + 1) % 3],
                                &table.pts[x as usize],
                            )
                            .is_zero()
                        };
                        let side = if cut.marked[t] {
                            (0..3).find(|&side| on(side, a) && on(side, b))
                        } else {
                            Some(k)
                        };
                        let Some(&j) = side.and_then(|side| edges.get(&(o.across[t][side], b, a)))
                        else {
                            ogeom_bail!(Invariant, "the cut {NAME} mesh is open");
                        };
                        j
                    }
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

/// The triangles coplanar with piece `piece`'s triangle that hold it, each
/// with whether it faces the same way.
///
/// # Errors
///
/// [`ogeom_core::OgeomError::Invariant`] if the piece straddles an edge of
/// one: the cut along that edge was missed.
fn stack(
    o: &Operand<'_>,
    table: &Table,
    flat: &[usize],
    marks: &[Mark],
    piece: &Piece,
) -> OgeomResult<Vec<(u32, bool)>> {
    let t = piece.parent as usize;
    let drop = flat[t];
    let corner = piece.v.map(|v| &table.pts[v as usize]);
    let mid = Mid(corner);
    let own = normal_sign(&o.tris[t].p, drop);
    let mut out = Vec::new();
    for &u in &marks[t].coplanar {
        let q = &o.tris[u as usize];
        let qc = q.p.map(|p| Pt::new(Def::Vertex(p)));
        let facing = normal_sign(&q.p, drop);
        let s = [0, 1, 2].map(|k| orient2d(drop, &qc[k], &qc[(k + 1) % 3], &mid));
        if s.contains(&facing.reversed()) {
            continue;
        }
        if s.contains(&Sign::Zero) {
            ogeom_bail!(Invariant, "a piece straddles a coplanar triangle's edge");
        }
        out.push((u, facing == own));
    }
    Ok(out)
}

/// The pieces kept, each as its corners, turned to face out of the region
/// the mesh winds about, and whether it is a piece of a cut triangle, in
/// the order of the pieces.
fn keep(
    o: &Operand<'_>,
    table: &Table,
    flat: &[usize],
    marks: &[Mark],
    cut: &Cut,
    regions: &[Vec<usize>],
) -> OgeomResult<Vec<([u32; 3], bool)>> {
    let pieces = &cut.pieces;
    let turned = parallel::map_ordered(regions, |_, region| -> OgeomResult<Option<bool>> {
        let piece = &pieces[region[0]];
        let t = piece.parent as usize;
        let drop = flat[t];
        let on = stack(o, table, flat, marks, piece)?;
        let mut skip: Vec<u32> = on.iter().map(|&(u, _)| u).collect();
        skip.push(piece.parent);
        let corner = piece.v.map(|v| &table.pts[v as usize]);
        let mid = Mid(corner);
        let ahead = normal_sign(&o.tris[t].p, drop);
        let front = classify::winding(&o.tris, &o.bvh, &mid, drop, ahead, &skip)?;
        let behind = front
            + 1
            + on.iter()
                .map(|&(_, same)| if same { 1 } else { -1 })
                .sum::<i64>();
        let (inside_front, inside_behind) = (front >= 1, behind >= 1);
        Ok((inside_front != inside_behind).then_some(inside_front))
    });
    let mut chosen = vec![None; pieces.len()];
    for (region, turn) in regions.iter().zip(turned) {
        if let Some(flip) = turn? {
            for &i in region {
                chosen[i] = Some(flip);
            }
        }
    }
    let which: Vec<usize> = (0..pieces.len()).filter(|&i| chosen[i].is_some()).collect();
    let lowest = parallel::map_ordered(&which, |_, &i| -> OgeomResult<bool> {
        let piece = &pieces[i];
        if marks[piece.parent as usize].coplanar.is_empty() {
            return Ok(true);
        }
        let on = stack(o, table, flat, marks, piece)?;
        Ok(on.iter().all(|&(u, _)| piece.parent < u))
    });
    let mut kept = Vec::with_capacity(which.len());
    for (&i, lowest) in which.iter().zip(lowest) {
        if !lowest? {
            continue;
        }
        let p = &pieces[i];
        let v = if chosen[i] == Some(true) {
            [p.v[0], p.v[2], p.v[1]]
        } else {
            p.v
        };
        kept.push((v, p.cut));
    }
    Ok(kept)
}

#[cfg(test)]
mod tests;
