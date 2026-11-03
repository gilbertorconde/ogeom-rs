//! How a network's patches, curves and corners connect: the vertices the
//! curves meet at, the strips that share knots, and the chains of tangent
//! curves that share a ratio across them.

use ogeom_core::{OgeomResult, ogeom_bail};

use super::{PatchDef, PatchSide};

/// A union-find over `0..n`, its classes named by their first member's
/// order of appearance, and an optional parity between members.
struct Classes {
    parent: Vec<usize>,
    /// The parity of each member relative to its parent.
    parity: Vec<bool>,
}

impl Classes {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            parity: vec![false; n],
        }
    }

    /// The class's root and the member's parity relative to it.
    fn find(&mut self, x: usize) -> (usize, bool) {
        let mut root = x;
        let mut parity = false;
        while self.parent[root] != root {
            parity ^= self.parity[root];
            root = self.parent[root];
        }
        // Compress, keeping each member's parity to the root.
        let mut at = x;
        let mut left = parity;
        while self.parent[at] != root && self.parent[at] != at {
            let next = self.parent[at];
            let own = self.parity[at];
            self.parent[at] = root;
            self.parity[at] = left;
            left ^= own;
            at = next;
        }
        (root, parity)
    }

    /// Join the classes of `a` and `b`, `b`'s parity relative to `a` being
    /// `odd`. False where they were already joined at the other parity.
    fn join(&mut self, a: usize, b: usize, odd: bool) -> bool {
        let (ra, pa) = self.find(a);
        let (rb, pb) = self.find(b);
        if ra == rb {
            return (pa ^ pb) == odd;
        }
        let (keep, drop) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.parent[drop] = keep;
        self.parity[drop] = pa ^ pb ^ odd;
        true
    }

    /// Each member's class, numbered from nought in order of first
    /// appearance over `members`, and the number of classes.
    fn numbered(&mut self, members: impl Iterator<Item = usize>) -> (Vec<Option<usize>>, usize) {
        let n = self.parent.len();
        let mut name = vec![None; n];
        let mut of = vec![None; n];
        let mut count = 0;
        for m in members {
            let (root, _) = self.find(m);
            let id = *name[root].get_or_insert_with(|| {
                count += 1;
                count - 1
            });
            of[m] = Some(id);
        }
        (of, count)
    }
}

/// The connections of a network.
#[derive(Debug, Clone)]
pub(super) struct Topology {
    /// The number of vertices.
    pub vertices: usize,
    /// Each patch's corners' vertices: corner `k` begins side `k`.
    pub corner_vertex: Vec<[usize; 4]>,
    /// Each curve's start and end vertex.
    pub curve_ends: Vec<(usize, usize)>,
    /// The patch sides each curve bounds, `(patch, side)`: one for a curve
    /// on the network's border, two for a curve between patches.
    pub users: Vec<Vec<(usize, usize)>>,
    /// The number of strips.
    pub strips: usize,
    /// Each patch's strips, along `u` then `v`.
    pub patch_strips: Vec<(usize, usize)>,
    /// Each curve's strip.
    pub curve_strip: Vec<usize>,
    /// Whether each vertex is the point a patch's missing side shrinks to.
    pub apex: Vec<bool>,
    /// Each vertex's patch count where every curve at it lies between two
    /// patches; `None` on the network's border.
    pub interior_valence: Vec<Option<usize>>,
    /// Whether each curve joins its patches smoothly: between two patches
    /// and not a crease.
    pub tangent: Vec<bool>,
    /// Whether each smooth curve is held tangent by the pole triples:
    /// neither of its ends is an apex or a corner of other than four
    /// patches.
    pub held: Vec<bool>,
    /// Each held curve's chain, and whether its ratio is the chain's
    /// inverse.
    pub chain: Vec<Option<(usize, bool)>>,
    /// The number of chains.
    pub chains: usize,
    /// The chains whose curves meet at the other parity somewhere, which
    /// only a ratio of one serves.
    pub conflicted: Vec<bool>,
}

/// A side's start and end, as curve ends: `2c` the start, `2c + 1` the end.
fn walk_ends(side: PatchSide) -> (usize, usize) {
    let (s, e) = (2 * side.curve, 2 * side.curve + 1);
    if side.reversed { (e, s) } else { (s, e) }
}

impl Topology {
    /// Read the connections off the patches.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// as [`fit_network`](super::fit_network) describes.
    pub fn of(patches: &[PatchDef], curves: usize, creases: &[usize]) -> OgeomResult<Self> {
        if patches.is_empty() {
            ogeom_bail!(Construction, "a network needs at least one patch");
        }
        let mut users: Vec<Vec<(usize, usize)>> = vec![Vec::new(); curves];
        for (p, patch) in patches.iter().enumerate() {
            let missing = patch.sides.iter().filter(|s| s.is_none()).count();
            if missing > 1 {
                ogeom_bail!(
                    Construction,
                    "patch {p} has {missing} missing sides; at most one may shrink to a point"
                );
            }
            for (k, side) in patch.sides.iter().enumerate() {
                let Some(side) = side else { continue };
                if side.curve >= curves {
                    ogeom_bail!(
                        Construction,
                        "patch {p} names curve {}, of {curves}",
                        side.curve
                    );
                }
                if users[side.curve].iter().any(|&(q, _)| q == p) {
                    ogeom_bail!(
                        Construction,
                        "patch {p} is bounded twice by curve {}",
                        side.curve
                    );
                }
                users[side.curve].push((p, k));
            }
        }
        for &c in creases {
            if c >= curves {
                ogeom_bail!(Construction, "crease {c} names no curve, of {curves}");
            }
        }
        let side_of = |(p, k): (usize, usize)| -> PatchSide {
            patches[p].sides[k].unwrap_or(PatchSide {
                curve: 0,
                reversed: false,
            })
        };
        for (c, used) in users.iter().enumerate() {
            match used.as_slice() {
                [] => ogeom_bail!(Construction, "curve {c} bounds no patch"),
                [_] => {}
                [a, b] => {
                    if side_of(*a).reversed == side_of(*b).reversed {
                        ogeom_bail!(
                            Construction,
                            "patches {} and {} walk curve {c} the same way; \
                             neighbours walk a shared curve in opposite directions",
                            a.0,
                            b.0
                        );
                    }
                }
                _ => ogeom_bail!(
                    Construction,
                    "curve {c} bounds {} patches; a curve bounds one or two",
                    used.len()
                ),
            }
        }

        // Vertices: the curve ends each corner joins.
        let mut ends = Classes::new(2 * curves);
        for patch in patches {
            for k in 0..4 {
                let before = patch.sides[(k + 3) % 4].map(|s| walk_ends(s).1);
                let after = patch.sides[k].map(|s| walk_ends(s).0);
                if let (Some(a), Some(b)) = (before, after) {
                    ends.join(a, b, false);
                }
            }
            if let Some(k) = patch.sides.iter().position(Option::is_none) {
                // The missing side's two corners are one point.
                let before = patch.sides[(k + 3) % 4].map(|s| walk_ends(s).1);
                let after = patch.sides[(k + 1) % 4].map(|s| walk_ends(s).0);
                if let (Some(a), Some(b)) = (before, after) {
                    ends.join(a, b, false);
                }
            }
        }
        let (vertex_of_end, vertices) = ends.numbered(0..2 * curves);
        let vertex = |end: usize| vertex_of_end[end].unwrap_or(0);
        let curve_ends: Vec<(usize, usize)> = (0..curves)
            .map(|c| (vertex(2 * c), vertex(2 * c + 1)))
            .collect();
        if let Some(c) = curve_ends.iter().position(|(a, b)| a == b) {
            ogeom_bail!(Construction, "curve {c} closes on itself");
        }
        let mut corner_vertex = Vec::with_capacity(patches.len());
        let mut apex = vec![false; vertices];
        for (p, patch) in patches.iter().enumerate() {
            let mut corners = [0; 4];
            for (k, corner) in corners.iter_mut().enumerate() {
                *corner = match (patch.sides[k], patch.sides[(k + 3) % 4]) {
                    (Some(s), _) => vertex(walk_ends(s).0),
                    (None, Some(s)) => vertex(walk_ends(s).1),
                    (None, None) => 0,
                };
            }
            if let Some(k) = patch.sides.iter().position(Option::is_none) {
                apex[corners[k]] = true;
            }
            let distinct = {
                let mut seen = corners.to_vec();
                seen.sort_unstable();
                seen.dedup();
                seen.len()
            };
            let missing = patch.sides.iter().filter(|s| s.is_none()).count();
            if distinct != 4 - missing {
                ogeom_bail!(Construction, "patch {p} meets itself at a corner");
            }
            corner_vertex.push(corners);
        }

        // Strips: the curves opposite each other in a patch share knots.
        let mut strips = Classes::new(curves + 2 * patches.len());
        for (p, patch) in patches.iter().enumerate() {
            for dir in 0..2 {
                for k in [dir, dir + 2] {
                    if let Some(s) = patch.sides[k] {
                        strips.join(curves + 2 * p + dir, s.curve, false);
                    }
                }
            }
        }
        let (strip_of, strip_count) = strips.numbered(curves..curves + 2 * patches.len());
        let patch_strips = (0..patches.len())
            .map(|p| {
                (
                    strip_of[curves + 2 * p].unwrap_or(0),
                    strip_of[curves + 2 * p + 1].unwrap_or(0),
                )
            })
            .collect();
        let curve_strip: Vec<usize> = users
            .iter()
            .map(|u| {
                let (p, k) = u[0];
                strip_of[curves + 2 * p + k % 2].unwrap_or(0)
            })
            .collect();

        let tangent: Vec<bool> = (0..curves)
            .map(|c| users[c].len() == 2 && !creases.contains(&c))
            .collect();

        // Each vertex's corners and whether it lies inside the network.
        let mut corners_at: Vec<Vec<(usize, usize)>> = vec![Vec::new(); vertices];
        for (p, corners) in corner_vertex.iter().enumerate() {
            for (k, &v) in corners.iter().enumerate() {
                // A patch counts once at the point its missing side
                // shrinks to.
                if !corners_at[v].iter().any(|&(q, _)| q == p) {
                    corners_at[v].push((p, k));
                }
            }
        }
        let mut border = vec![false; vertices];
        for (c, &(a, b)) in curve_ends.iter().enumerate() {
            if users[c].len() < 2 {
                border[a] = true;
                border[b] = true;
            }
        }
        let interior_valence: Vec<Option<usize>> = (0..vertices)
            .map(|v| (!border[v] && !apex[v]).then_some(corners_at[v].len()))
            .collect();

        // A corner of other than four patches, or an apex, turns the
        // parameter lines that cross a curve one way on one side of it and
        // the other way on the other: no pole triple of the curve can hold
        // across it.
        let free_at = |v: usize| apex[v] || interior_valence[v].is_some_and(|n| n != 4);
        let held: Vec<bool> = (0..curves)
            .map(|c| {
                let (a, b) = curve_ends[c];
                tangent[c] && !free_at(a) && !free_at(b)
            })
            .collect();

        // Chains: held curves running straight on through a corner of four
        // patches share their ratio.
        let mut chains = Classes::new(curves);
        let mut conflicted_curves = Vec::new();
        let curve_of = |p: usize, k: usize| patches[p].sides[k].map(|s| s.curve);
        for v in 0..vertices {
            if interior_valence[v] != Some(4) {
                continue;
            }
            let at: Vec<(usize, [Option<usize>; 2])> = corners_at[v]
                .iter()
                .map(|&(p, k)| (p, [curve_of(p, (k + 3) % 4), curve_of(p, k)]))
                .collect();
            let Some(ring) = ring_round(&at) else {
                continue;
            };
            for i in 0..2 {
                let (c, far) = (ring[i].0, ring[i + 2].0);
                if !(held[c] && held[far]) {
                    continue;
                }
                // The patch between `c` and the curve after it, and the one
                // between that curve and `far`: the sides of `c` and `far`
                // facing the curve between them.
                let near_side = ring[i].1;
                let far_side = ring[i + 1].1;
                let flip_c = users[c][0].0 != near_side;
                let flip_far = users[far][0].0 != far_side;
                if !chains.join(c, far, flip_c ^ flip_far) {
                    conflicted_curves.push(c);
                }
            }
        }
        let (chain_of, chain_count) = chains.numbered((0..curves).filter(|&c| held[c]));
        let mut conflicted = vec![false; chain_count];
        for c in conflicted_curves {
            if let Some(id) = chain_of[c] {
                conflicted[id] = true;
            }
        }
        let chain = (0..curves)
            .map(|c| chain_of[c].map(|id| (id, chains.find(c).1)))
            .collect();

        Ok(Self {
            vertices,
            corner_vertex,
            curve_ends,
            users,
            strips: strip_count,
            patch_strips,
            curve_strip,
            apex,
            interior_valence,
            tangent,
            held,
            chain,
            chains: chain_count,
            conflicted,
        })
    }
}

/// The curves round a corner of four patches in turn, each with the patch
/// that follows it: `(curve, patch between it and the next)`, the first
/// repeated at the end twice over so `ring[i + 2]` always reads. `None`
/// where the corners do not close into one ring.
fn ring_round(at: &[(usize, [Option<usize>; 2])]) -> Option<Vec<(usize, usize)>> {
    let (first_patch, [Some(a), Some(b)]) = at[0] else {
        return None;
    };
    let mut ring = vec![(a, first_patch)];
    let mut used = vec![false; at.len()];
    used[0] = true;
    let mut current = b;
    while ring.len() < at.len() {
        let (i, other, patch) = at.iter().enumerate().find_map(|(i, &(p, [x, y]))| {
            if used[i] {
                return None;
            }
            match (x, y) {
                (Some(x), Some(y)) if x == current => Some((i, y, p)),
                (Some(x), Some(y)) if y == current => Some((i, x, p)),
                _ => None,
            }
        })?;
        used[i] = true;
        ring.push((current, patch));
        current = other;
    }
    if current != a {
        return None;
    }
    ring.push(ring[0]);
    ring.push(ring[1]);
    Some(ring)
}
