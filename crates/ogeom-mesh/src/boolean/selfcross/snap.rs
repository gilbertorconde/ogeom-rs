//! A resolved mesh rounded to the grid its caller stores it on.
//!
//! Every vertex is rounded to the nearest point of the grid, and the
//! result is checked on the rounded coordinates with the same exact
//! predicates that resolve it. Rounding moves a vertex by at most half a
//! grid step along each axis, which can do three things to a surface it
//! leaves embedded in `f64`, each found and mended here:
//!
//! - The two ends of an edge round to one point. The edge is collapsed,
//!   the triangles on it dropped: they had no area left.
//! - A triangle's corners round onto one line. It is flipped into the
//!   triangle across its longest side, which moves nothing.
//! - A sliver tilts through its neighbour, a shell narrower than a step
//!   turns inside out, or two vertices with no edge between them round to
//!   one point. The rounded mesh is resolved again from its rounded
//!   coordinates, and the new vertices rounded in turn.
//!
//! Where resolving the rounded mesh would leave parts touching along an
//! edge or at a point, the whole resolution is begun again from the input
//! with its vertices rounded, so the cuts fall on the grid's own points.
//!
//! This is iterated snap rounding (D. Halperin and E. Packer, *Iterated
//! snap rounding*; S. Fortune, *Vertex-rounding a three-dimensional
//! polyhedral subdivision*): rounding, then cutting again where rounding
//! made crossings, until none is left. Each round's pieces lie on the
//! triangles of the round before, so no point of the result is farther
//! from the input surface than the rounding of every round added up.

use ogeom_core::{FastMap, FastSet, OgeomError, OgeomResult, ogeom_bail, ogeom_err};
use ogeom_topo::Triangulation;

use super::super::points::normal_sign;
use super::super::{Operand, coords, find, id, union};
use super::{NAME, PASSES, Resolved, carried, crossing_pairs, resolve_once};
use ogeom_math::Point;

/// The grid the vertices of a resolved mesh are stored on, for
/// [`resolve_self_intersections_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Quantum {
    /// Double precision: only the new vertices are rounded, to the nearest
    /// `f64`, as [`resolve_self_intersections`](super::resolve_self_intersections)
    /// does.
    F64,
    /// Single precision: every coordinate is an `f32`, the nearest to the
    /// `f64` it rounds.
    F32,
    /// Every coordinate a whole multiple of the given spacing, which must
    /// be positive and finite.
    Grid(f64),
}

impl Quantum {
    /// The point of the grid nearest to `x`, as an `f64`; `None` where the
    /// grid holds no such point: outside the range of `f32`, or more than
    /// 2⁵⁰ steps of a `Grid` from zero.
    #[must_use]
    pub fn round(self, x: f64) -> Option<f64> {
        if !x.is_finite() {
            return None;
        }
        match self {
            Self::F64 => Some(x + 0.0),
            #[allow(clippy::cast_possible_truncation)]
            Self::F32 => {
                let y = x as f32;
                y.is_finite().then(|| f64::from(y) + 0.0)
            }
            Self::Grid(h) => {
                let k = (x / h).round();
                (h > 0.0 && h.is_finite() && k.abs() <= GRID_STEPS).then_some(k * h + 0.0)
            }
        }
    }

    /// The farthest [`round`](Self::round) moves a coordinate no larger in
    /// magnitude than `magnitude`: half the spacing of the grid there. A
    /// vertex moves at most `√3` times this.
    #[must_use]
    pub fn half_step(self, magnitude: f64) -> f64 {
        let m = magnitude.abs();
        match self {
            Self::F64 => {
                let up = f64::from_bits(m.to_bits() + 1);
                (up - m) / 2.0
            }
            #[allow(clippy::cast_possible_truncation)]
            Self::F32 => {
                // Rounded to the nearest, `v` is in `m`'s binade or the
                // next one up, so its spacing is never the smaller.
                let v = m as f32;
                let up = f32::from_bits(v.to_bits() + 1);
                (f64::from(up) - f64::from(v)) / 2.0
            }
            Self::Grid(h) => h / 2.0,
        }
    }
}

impl std::fmt::Display for Quantum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::F64 => f.write_str("double precision"),
            Self::F32 => f.write_str("single precision"),
            Self::Grid(h) => write!(f, "a grid of spacing {h}"),
        }
    }
}

/// How many steps of a `Grid` from zero a coordinate may lie: well inside
/// the integers an `f64` holds, so `k · h` rounds back to `k`.
const GRID_STEPS: f64 = 1_125_899_906_842_624.0;

/// [`resolve_self_intersections`](super::resolve_self_intersections), its
/// result's vertices on the grid `quantum`.
///
/// With [`Quantum::F64`] this is `resolve_self_intersections`, bit for
/// bit. Otherwise every vertex of the result, kept or new, is rounded to
/// the grid, and the result is checked on the rounded coordinates by the
/// exact predicates of [`self_intersections`](super::self_intersections):
/// no two triangles meet but at the corners they share, no triangle is
/// flat, and every shell faces out of what the mesh bounds. An edge whose
/// ends round to one point is collapsed and a triangle rounded flat is
/// flipped away; where the rounded mesh still crosses itself, or a shell
/// rounded inside out, it is resolved again from the rounded coordinates
/// and rounded again, at most three more times. Where that would leave
/// parts touching along an edge or at a point, it is all begun again once
/// from the input with its vertices rounded.
///
/// No point of the result lies farther from the input surface than
/// [`Resolved::rounding`], which is at most five times `√3 ·
/// quantum.half_step(m)`, `m` the largest coordinate of the input and the
/// result, plus the `f64` rounding of the new vertices.
///
/// # Errors
///
/// As [`resolve_self_intersections`](super::resolve_self_intersections);
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a `Grid`
/// spacing is not positive and finite, or a vertex lies outside what the
/// grid holds; and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// rounding joins vertices so that the surface is not 2-manifold, if
/// resolving the rounded result again would leave parts touching, or if
/// after the last pass it still crosses itself or faces inward; these
/// messages begin "rounding the result to" and name the grid.
pub fn resolve_self_intersections_with(
    mesh: &Triangulation,
    quantum: Quantum,
) -> OgeomResult<Resolved> {
    if quantum == Quantum::F64 {
        return super::resolve_f64(mesh);
    }
    if let Quantum::Grid(h) = quantum
        && !(h > 0.0 && h.is_finite())
    {
        ogeom_bail!(Range, "a grid spacing must be positive and finite, not {h}");
    }
    let (first, crossings) = resolve_once(mesh)?;
    let start = Resolved {
        mesh: first.mesh,
        crossings,
        cut: first.cut,
        rounding: first.rounding,
    };
    match rounds(start, quantum) {
        Err(OgeomError::Construction(cause)) => from_rounded_input(mesh, crossings, quantum)
            .map_err(|_| OgeomError::Construction(cause)),
        done => done,
    }
}

/// `out` rounded to `quantum` and resolved again from the rounded
/// coordinates until it is embedded on them, at most [`PASSES`] times.
fn rounds(mut out: Resolved, quantum: Quantum) -> OgeomResult<Resolved> {
    for pass in 0..=PASSES {
        out.rounding += snap(&mut out.mesh, &mut out.cut, quantum)?;
        if settled(&out.mesh, quantum)? {
            break;
        }
        if pass == PASSES {
            ogeom_bail!(
                Construction,
                "rounding the result to {quantum} leaves it crossing itself"
            );
        }
        let (next, _) = resolve_once(&out.mesh).map_err(|e| {
            ogeom_err!(
                Construction,
                "rounding the result to {quantum} leaves it crossing itself, and resolving it again fails: {e}"
            )
        })?;
        out.cut = carried(&out.mesh, &out.cut, &next);
        out.rounding += next.rounding;
        out.mesh = next.mesh;
    }
    Ok(out)
}

/// The resolution of `mesh` begun again from its vertices rounded to
/// `quantum`, for where rounding the resolution of the vertices as given
/// leaves parts touching: the cuts then fall on the grid's own points.
/// Triangles that are not input triangles, by their corners, are flagged
/// cut.
fn from_rounded_input(
    mesh: &Triangulation,
    crossings: usize,
    quantum: Quantum,
) -> OgeomResult<Resolved> {
    let mut start = Triangulation::new();
    start.positions.clone_from(&mesh.positions);
    start.triangles.clone_from(&mesh.triangles);
    let mut cut = vec![false; start.triangles.len()];
    let moved = snap(&mut start, &mut cut, quantum)?;
    let (first, _) = resolve_once(&start)?;
    let mut out = rounds(
        Resolved {
            mesh: first.mesh,
            crossings,
            cut: first.cut,
            rounding: moved + first.rounding,
        },
        quantum,
    )?;
    let key = |m: &Triangulation, t: [u32; 3]| {
        let mut c = t.map(|v| coords(m.positions[v as usize]).map(f64::to_bits));
        c.sort_unstable();
        c
    };
    let input: FastSet<[[u64; 3]; 3]> = mesh.triangles.iter().map(|&t| key(mesh, t)).collect();
    for (t, c) in out.mesh.triangles.iter().zip(out.cut.iter_mut()) {
        *c = !input.contains(&key(&out.mesh, *t));
    }
    Ok(out)
}

/// Whether the closed `mesh`, its vertices on the grid, is embedded and
/// faces outward, decided exactly on its coordinates.
fn settled(mesh: &Triangulation, quantum: Quantum) -> OgeomResult<bool> {
    if mesh.triangles.is_empty() {
        return Ok(true);
    }
    let o = Operand::new(mesh, 0, NAME).map_err(|e| {
        ogeom_err!(
            Construction,
            "rounding the result to {quantum} leaves it not 2-manifold: {e}"
        )
    })?;
    Ok(crossing_pairs(&o.tris, &o.bvh).is_empty() && o.check_wound(NAME, None, 0).is_ok())
}

/// Round every vertex of the closed, 2-manifold `mesh` to `quantum`,
/// collapse each edge whose ends round to one point, drop pairs of
/// triangles that collapsing left back to back, flip away triangles left
/// flat, and drop vertices no triangle uses. Answers the farthest a vertex
/// moved.
fn snap(mesh: &mut Triangulation, cut: &mut Vec<bool>, quantum: Quantum) -> OgeomResult<f64> {
    let mut moved = 0.0_f64;
    for p in &mut mesh.positions {
        let q = coords(*p).map(|x| quantum.round(x));
        let [Some(x), Some(y), Some(z)] = q else {
            ogeom_bail!(
                Range,
                "a vertex of the result, ({}, {}, {}), lies outside what {quantum} holds",
                p.x,
                p.y,
                p.z
            );
        };
        let r = Point::new(x, y, z);
        moved = moved.max(p.distance(r));
        *p = r;
    }
    let at = |v: u32| coords(mesh.positions[v as usize]).map(f64::to_bits);
    let mut parent: Vec<u32> = (0..id(mesh.positions.len())).collect();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if at(a) == at(b) {
                union(&mut parent, a, b);
            }
        }
    }
    let mut kept: Vec<([u32; 3], bool)> = Vec::with_capacity(mesh.triangles.len());
    for (t, &c) in mesh.triangles.iter().zip(cut.iter()) {
        let v = t.map(|x| find(&mut parent, x));
        if v[0] != v[1] && v[1] != v[2] && v[2] != v[0] {
            kept.push((v, c));
        }
    }
    flatten(&mesh.positions, &mut kept, quantum)?;
    let mut slot = vec![u32::MAX; mesh.positions.len()];
    for (v, _) in &kept {
        for &x in v {
            slot[x as usize] = 0;
        }
    }
    let mut positions = Vec::with_capacity(mesh.positions.len());
    for (x, s) in slot.iter_mut().enumerate() {
        if *s == 0 {
            *s = id(positions.len());
            positions.push(mesh.positions[x]);
        }
    }
    mesh.positions = positions;
    mesh.triangles = kept
        .iter()
        .map(|(v, _)| v.map(|x| slot[x as usize]))
        .collect();
    *cut = kept.iter().map(|&(_, c)| c).collect();
    closed(mesh, quantum)?;
    Ok(moved)
}

/// Remove from `kept` every triangle whose corners, at `positions`, are
/// collinear, and every pair left back to back.
///
/// A flat triangle is flipped into the triangle across its longest side,
/// as [`flip_needles`](super::flip_needles) does, which moves nothing.
/// Where the new diagonal is already an edge, the flip is made anyway:
/// what it doubles is a pair back to back, collapsed out of a sliver, and
/// is dropped with the rest.
fn flatten(
    positions: &[Point],
    kept: &mut Vec<([u32; 3], bool)>,
    quantum: Quantum,
) -> OgeomResult<()> {
    let corners = |v: [u32; 3]| v.map(|x| coords(positions[x as usize]));
    let needle = |v: [u32; 3]| {
        let p = corners(v);
        (0..3).all(|k| normal_sign(&p, k).is_zero())
    };
    let mut flips = 0;
    loop {
        drop_back_to_back(kept);
        let todo: Vec<usize> = (0..kept.len()).filter(|&t| needle(kept[t].0)).collect();
        if todo.is_empty() {
            return Ok(());
        }
        let mut edges: FastMap<(u32, u32), Vec<usize>> = FastMap::default();
        for (t, (v, _)) in kept.iter().enumerate() {
            for k in 0..3 {
                edges.entry((v[k], v[(k + 1) % 3])).or_default().push(t);
            }
        }
        for t in todo {
            let v = kept[t].0;
            if !needle(v) {
                continue;
            }
            flips += 1;
            if flips > 4 * positions.len() + kept.len() {
                ogeom_bail!(Invariant, "flipping flat triangles does not end");
            }
            let p = corners(v);
            let middle = (0..3).find(|&k| {
                let (a, b, c) = (p[(k + 1) % 3], p[(k + 2) % 3], p[k]);
                (0..3)
                    .find(|&x| a[x] != b[x])
                    .is_some_and(|x| (a[x] < c[x] && c[x] < b[x]) || (b[x] < c[x] && c[x] < a[x]))
            });
            let Some(k) = middle else {
                ogeom_bail!(
                    Invariant,
                    "a triangle of the rounded result has two corners at one point"
                );
            };
            let (c, a, b) = (v[k], v[(k + 1) % 3], v[(k + 2) % 3]);
            let Some(&u) = edges.get(&(b, a)).and_then(|across| across.first()) else {
                ogeom_bail!(
                    Construction,
                    "rounding the result to {quantum} joins its vertices so that it is not 2-manifold"
                );
            };
            let w = kept[u].0;
            let d = w.into_iter().find(|&x| x != a && x != b).unwrap_or(a);
            if d == c {
                ogeom_bail!(
                    Construction,
                    "rounding the result to {quantum} leaves a flat triangle that cannot be removed"
                );
            }
            for (x, at) in [(v, t), (w, u)] {
                for k in 0..3 {
                    if let Some(list) = edges.get_mut(&(x[k], x[(k + 1) % 3])) {
                        list.retain(|&e| e != at);
                    }
                }
            }
            kept[t] = ([b, c, d], true);
            kept[u] = ([c, a, d], true);
            for (x, at) in [([b, c, d], t), ([c, a, d], u)] {
                for k in 0..3 {
                    edges.entry((x[k], x[(k + 1) % 3])).or_default().push(at);
                }
            }
        }
    }
}

/// Remove each pair of triangles on the same three vertices wound
/// opposite ways: back to back, they bound nothing.
fn drop_back_to_back(kept: &mut Vec<([u32; 3], bool)>) {
    let key = |v: [u32; 3]| {
        let mut s = v;
        s.sort_unstable();
        s
    };
    // A triangle's winding, as whether its corners run in the order of
    // their sorted numbers up to rotation.
    let even = |v: [u32; 3]| {
        let s = key(v);
        (0..3).any(|r| [v[r], v[(r + 1) % 3], v[(r + 2) % 3]] == s)
    };
    let mut open: FastMap<[u32; 3], [Vec<usize>; 2]> = FastMap::default();
    let mut gone = vec![false; kept.len()];
    for (i, &(v, _)) in kept.iter().enumerate() {
        let side = usize::from(even(v));
        let slots = open.entry(key(v)).or_default();
        if let Some(j) = slots[1 - side].pop() {
            gone[i] = true;
            gone[j] = true;
        } else {
            slots[side].push(i);
        }
    }
    let mut i = 0;
    kept.retain(|_| {
        i += 1;
        !gone[i - 1]
    });
}

/// Every edge of `mesh` is run once each way.
fn closed(mesh: &Triangulation, quantum: Quantum) -> OgeomResult<()> {
    let mut edges: FastMap<(u32, u32), ()> = FastMap::default();
    for v in &mesh.triangles {
        for k in 0..3 {
            if edges.insert((v[k], v[(k + 1) % 3]), ()).is_some() {
                ogeom_bail!(
                    Construction,
                    "rounding the result to {quantum} joins its vertices so that it is not 2-manifold"
                );
            }
        }
    }
    if edges.keys().any(|&(a, b)| !edges.contains_key(&(b, a))) {
        ogeom_bail!(
            Construction,
            "rounding the result to {quantum} joins its vertices so that it is not 2-manifold"
        );
    }
    Ok(())
}
