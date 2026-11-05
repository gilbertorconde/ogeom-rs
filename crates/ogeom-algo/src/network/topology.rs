//! How a network's patches, curves and corners connect: the vertices the
//! curves meet at, the strips that share knots, and the chains of tangent
//! curves that share a ratio across them.

use ogeom_core::{OgeomResult, ogeom_bail};

use super::{ChainedPatch, PatchSide};

/// A union-find over `0..n`, its classes named by their first member's
/// order of appearance, and an optional parity between members.
pub(super) struct Classes {
    parent: Vec<usize>,
    /// The parity of each member relative to its parent.
    parity: Vec<bool>,
}

impl Classes {
    pub(super) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            parity: vec![false; n],
        }
    }

    /// The class's root and the member's parity relative to it.
    pub(super) fn find(&mut self, x: usize) -> (usize, bool) {
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
    pub(super) fn join(&mut self, a: usize, b: usize, odd: bool) -> bool {
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
    pub(super) fn numbered(
        &mut self,
        members: impl Iterator<Item = usize>,
    ) -> (Vec<Option<usize>>, usize) {
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

/// Where a patch touches a vertex: at its corner `k`, or straight through
/// it, at the junction after curve `j` of side `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Touch {
    Corner(usize),
    Through(usize, usize),
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
    /// The number of strips: runs of knot spans shared by every curve and
    /// patch that crosses them.
    pub strips: usize,
    /// Each patch's strips along `u` and along `v`, in order: one, or one
    /// per curve of a side split at junctions.
    pub patch_strips: Vec<[Vec<usize>; 2]>,
    /// Each curve's strips, in its own order.
    pub curve_strip: Vec<Vec<usize>>,
    /// Whether each vertex is the point a patch's missing side shrinks to.
    pub apex: Vec<bool>,
    /// Each vertex's patch count where every curve at it lies between two
    /// patches, a patch it lies straight on the side of counting twice;
    /// `None` on the network's border.
    pub interior_valence: Vec<Option<usize>>,
    /// The patches at each vertex and where they touch it.
    pub touching: Vec<Vec<(usize, Touch)>>,
    /// Whether each curve joins its patches smoothly: between two patches
    /// and not a crease.
    pub tangent: Vec<bool>,
    /// Whether each patch's strips along `u` and along `v` run against
    /// their strips' own direction.
    pub patch_flip: Vec<[Vec<bool>; 2]>,
    /// Whether each curve's strips run against their own direction.
    pub curve_flip: Vec<Vec<bool>>,
    /// The curves held tangent by the pole triples where no corner is held
    /// by its spans: neither of a curve's ends is an apex or a corner of
    /// other than four patches.
    pub plain: Holding,
    /// The curves held tangent by the pole triples where the
    /// [`span_corners`](Self::span_corners) are held: as
    /// [`plain`](Self::plain), and every curve into such a corner, from its
    /// span next to the corner on.
    pub spanned: Holding,
    /// The corners of an odd number of patches, other than three-sided
    /// patches' points, whose curves are all smooth and all end on the
    /// network's border at their other ends: each with its curves.
    pub odd_corners: Vec<(usize, Vec<usize>)>,
    /// The corners of other than four patches, more than four where even,
    /// other than three-sided patches' points and junctions, whose curves
    /// are all smooth, each crossing one strip that runs one way
    /// throughout, and each ending at its other end on the network's
    /// border, at a corner of four patches or a junction, or at another
    /// such corner: each with its curves.
    pub span_corners: Vec<(usize, Vec<usize>)>,
}

/// Which curves the pole triples hold tangent, and how their ratios are
/// shared.
#[derive(Debug, Clone)]
pub(super) struct Holding {
    /// Whether each smooth curve is held tangent by the pole triples.
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

/// A side's first curve's start and last curve's end, as curve ends.
fn side_ends(side: &[PatchSide]) -> Option<(usize, usize)> {
    Some((walk_ends(*side.first()?).0, walk_ends(*side.last()?).1))
}

/// How many strips each direction of patch `p` crosses, given each curve's:
/// those its side's curves cross, the most of its two sides'.
fn crossed(patch: &ChainedPatch, dir: usize, parts: &[usize]) -> usize {
    [dir, dir + 2]
        .iter()
        .map(|&k| side_parts(&patch.sides[k], parts))
        .max()
        .unwrap_or(1)
}

/// How many strips a side crosses: those of its curves, one where it is
/// missing.
fn side_parts(side: &[PatchSide], parts: &[usize]) -> usize {
    side.iter().map(|s| parts[s.curve]).sum::<usize>().max(1)
}

impl Topology {
    /// Read the connections off the patches.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// as [`fit_network`](super::fit_network) describes.
    #[allow(clippy::too_many_lines, reason = "one pass per connection")]
    pub fn of(patches: &[ChainedPatch], curves: usize, creases: &[usize]) -> OgeomResult<Self> {
        if patches.is_empty() {
            ogeom_bail!(Construction, "a network needs at least one patch");
        }
        let mut users: Vec<Vec<(usize, usize)>> = vec![Vec::new(); curves];
        for (p, patch) in patches.iter().enumerate() {
            let missing = patch.sides.iter().filter(|s| s.is_empty()).count();
            if missing > 1 {
                ogeom_bail!(
                    Construction,
                    "patch {p} has {missing} missing sides; at most one may shrink to a point"
                );
            }
            for (k, side) in patch.sides.iter().enumerate() {
                for piece in side {
                    if piece.curve >= curves {
                        ogeom_bail!(
                            Construction,
                            "patch {p} names curve {}, of {curves}",
                            piece.curve
                        );
                    }
                    if users[piece.curve].iter().any(|&(q, _)| q == p) {
                        ogeom_bail!(
                            Construction,
                            "patch {p} is bounded twice by curve {}",
                            piece.curve
                        );
                    }
                    users[piece.curve].push((p, k));
                }
            }
        }
        for &c in creases {
            if c >= curves {
                ogeom_bail!(Construction, "crease {c} names no curve, of {curves}");
            }
        }
        let piece_of = |(p, k): (usize, usize), c: usize| -> PatchSide {
            patches[p].sides[k]
                .iter()
                .copied()
                .find(|s| s.curve == c)
                .unwrap_or(PatchSide {
                    curve: c,
                    reversed: false,
                })
        };
        for (c, used) in users.iter().enumerate() {
            match used.as_slice() {
                [] => ogeom_bail!(Construction, "curve {c} bounds no patch"),
                [_] => {}
                [a, b] => {
                    if piece_of(*a, c).reversed == piece_of(*b, c).reversed {
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

        // Vertices: the curve ends each corner joins, and the junctions
        // between the curves of a split side.
        let mut ends = Classes::new(2 * curves);
        for patch in patches {
            for k in 0..4 {
                let before = side_ends(&patch.sides[(k + 3) % 4]).map(|e| e.1);
                let after = side_ends(&patch.sides[k]).map(|e| e.0);
                if let (Some(a), Some(b)) = (before, after) {
                    ends.join(a, b, false);
                }
                for w in patch.sides[k].windows(2) {
                    ends.join(walk_ends(w[0]).1, walk_ends(w[1]).0, false);
                }
            }
            if let Some(k) = patch.sides.iter().position(Vec::is_empty) {
                // The missing side's two corners are one point.
                let before = side_ends(&patch.sides[(k + 3) % 4]).map(|e| e.1);
                let after = side_ends(&patch.sides[(k + 1) % 4]).map(|e| e.0);
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
                *corner = match (
                    side_ends(&patch.sides[k]),
                    side_ends(&patch.sides[(k + 3) % 4]),
                ) {
                    (Some(s), _) => vertex(s.0),
                    (None, Some(s)) => vertex(s.1),
                    (None, None) => 0,
                };
            }
            if let Some(k) = patch.sides.iter().position(Vec::is_empty) {
                apex[corners[k]] = true;
            }
            let mut seen = corners.to_vec();
            for side in &patch.sides {
                seen.extend(side.windows(2).map(|w| vertex(walk_ends(w[0]).1)));
            }
            let junctions = seen.len() - 4;
            seen.sort_unstable();
            seen.dedup();
            let missing = patch.sides.iter().filter(|s| s.is_empty()).count();
            if seen.len() != 4 - missing + junctions {
                ogeom_bail!(Construction, "patch {p} meets itself at a corner");
            }
            corner_vertex.push(corners);
        }

        // How many strips each curve crosses: one, or as many as the side
        // opposite it crosses, carried on from patch to patch through the
        // curves they share. Each strip a curve crosses past its first
        // begins on a line from a junction between the curves of a split
        // side; a curve crossing more than there are such junctions meets
        // one line twice, which the junctions carry round a loop.
        let windows: usize = patches
            .iter()
            .flat_map(|patch| patch.sides.iter())
            .map(|side| side.len().saturating_sub(1))
            .sum();
        let mut parts = vec![1usize; curves];
        loop {
            if let Some(c) = parts.iter().position(|&n| n > windows + 1) {
                ogeom_bail!(
                    Construction,
                    "curve {c} is crossed by more junction lines than there are junctions; \
                     the split sides carry their junctions round a loop"
                );
            }
            let mut changed = false;
            for patch in patches {
                for dir in 0..2 {
                    let m = crossed(patch, dir, &parts);
                    for k in [dir, dir + 2] {
                        if let [one] = patch.sides[k].as_slice()
                            && parts[one.curve] < m
                        {
                            parts[one.curve] = m;
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
        for (p, patch) in patches.iter().enumerate() {
            for dir in 0..2 {
                let m = crossed(patch, dir, &parts);
                for k in [dir, dir + 2] {
                    let side = &patch.sides[k];
                    let own = side_parts(side, &parts);
                    if side.len() > 1 && own != m {
                        ogeom_bail!(
                            Construction,
                            "patch {p}'s side {k} is split into curves crossing {own} strips \
                             and the side opposite into {m}; opposite sides split at \
                             junctions cross as many strips"
                        );
                    }
                }
            }
        }

        // Strips: the curves opposite each other in a patch share knots, a
        // split side's curves each with the strip of the patch it crosses.
        let mut part_base = Vec::with_capacity(curves);
        let mut next = 0;
        for &n in &parts {
            part_base.push(next);
            next += n;
        }
        let mut dir_base = Vec::with_capacity(patches.len());
        for patch in patches {
            let base = next;
            next += crossed(patch, 0, &parts) + crossed(patch, 1, &parts);
            dir_base.push(base);
        }
        let mut strips = Classes::new(next);
        let mut turned = Vec::new();
        let segment = |p: usize, dir: usize, i: usize| {
            dir_base[p]
                + if dir == 0 {
                    0
                } else {
                    crossed(&patches[p], 0, &parts)
                }
                + i
        };
        for (p, patch) in patches.iter().enumerate() {
            for dir in 0..2 {
                let m = crossed(patch, dir, &parts);
                for k in [dir, dir + 2] {
                    // Sides 0 and 1 walk with their direction, 2 and 3
                    // against it.
                    let along = |w: usize| if k < 2 { w } else { m - 1 - w };
                    // Each curve crosses as many strips as it has parts,
                    // from where the curves before it along the side end; a
                    // curve runs against the patch's direction where it runs
                    // against a walk that runs against it.
                    let mut w = 0;
                    for piece in &patch.sides[k] {
                        let n = parts[piece.curve];
                        for x in 0..n {
                            let part = if piece.reversed { n - 1 - x } else { x };
                            let (a, b) =
                                (segment(p, dir, along(w + x)), part_base[piece.curve] + part);
                            if !strips.join(a, b, (k >= 2) != piece.reversed) {
                                turned.push(a);
                            }
                        }
                        w += n;
                    }
                }
            }
        }
        // Numbered over the patches' strips first, every curve's part
        // sharing one of theirs.
        let (strip_of, strip_count) = strips.numbered((dir_base[0]..next).chain(0..dir_base[0]));
        let patch_strips: Vec<[Vec<usize>; 2]> = (0..patches.len())
            .map(|p| {
                core::array::from_fn(|dir| {
                    (0..crossed(&patches[p], dir, &parts))
                        .map(|i| strip_of[segment(p, dir, i)].unwrap_or(0))
                        .collect()
                })
            })
            .collect();
        let curve_strip: Vec<Vec<usize>> = (0..curves)
            .map(|c| {
                (0..parts[c])
                    .map(|i| strip_of[part_base[c] + i].unwrap_or(0))
                    .collect()
            })
            .collect();
        let patch_flip: Vec<[Vec<bool>; 2]> = (0..patches.len())
            .map(|p| {
                core::array::from_fn(|dir| {
                    (0..crossed(&patches[p], dir, &parts))
                        .map(|i| strips.find(segment(p, dir, i)).1)
                        .collect()
                })
            })
            .collect();
        let curve_flip: Vec<Vec<bool>> = (0..curves)
            .map(|c| {
                (0..parts[c])
                    .map(|i| strips.find(part_base[c] + i).1)
                    .collect()
            })
            .collect();
        // A strip met one way by one patch and the other way by another,
        // as round a band with a half twist, has no one direction.
        let mut two_way = vec![false; strip_count];
        for a in turned {
            if let Some(s) = strip_of[a] {
                two_way[s] = true;
            }
        }

        let tangent: Vec<bool> = (0..curves)
            .map(|c| users[c].len() == 2 && !creases.contains(&c))
            .collect();

        // Each vertex's patches and whether it lies inside the network.
        let mut touching: Vec<Vec<(usize, Touch)>> = vec![Vec::new(); vertices];
        for (p, corners) in corner_vertex.iter().enumerate() {
            for (k, &v) in corners.iter().enumerate() {
                // A patch counts once at the point its missing side
                // shrinks to.
                if !touching[v].iter().any(|&(q, _)| q == p) {
                    touching[v].push((p, Touch::Corner(k)));
                }
            }
            for (k, side) in patches[p].sides.iter().enumerate() {
                for (j, w) in side.windows(2).enumerate() {
                    touching[vertex(walk_ends(w[0]).1)].push((p, Touch::Through(k, j)));
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
            .map(|v| {
                (!border[v] && !apex[v]).then(|| {
                    touching[v]
                        .iter()
                        .map(|(_, t)| match t {
                            Touch::Corner(_) => 1,
                            Touch::Through(..) => 2,
                        })
                        .sum()
                })
            })
            .collect();

        // A corner of other than four patches, or an apex, turns the
        // parameter lines that cross a curve one way on one side of it and
        // the other way on the other: no pole triple of the curve can hold
        // across it.
        let free_at = |v: usize| apex[v] || interior_valence[v].is_some_and(|n| n != 4);
        let at_vertex = |v: usize| -> Vec<usize> {
            (0..curves)
                .filter(|&c| curve_ends[c].0 == v || curve_ends[c].1 == v)
                .collect()
        };
        let through = |v: usize| {
            touching[v]
                .iter()
                .any(|(_, t)| matches!(t, Touch::Through(..)))
        };
        let far_of = |c: usize, v: usize| {
            let (a, b) = curve_ends[c];
            if a == v { b } else { a }
        };
        // A corner of an odd number of patches (or, where `even`, of an even
        // number above four), not at a junction, whose curves are all smooth
        // and pass `reaches` at their other ends; with its curves.
        let odd_with =
            |even: bool, reaches: &dyn Fn(usize, usize) -> bool| -> Vec<(usize, Vec<usize>)> {
                (0..vertices)
                    .filter_map(|v| {
                        let valence = interior_valence[v]?;
                        if (valence % 2 == 0 && !(even && valence > 4)) || apex[v] || through(v) {
                            return None;
                        }
                        let at = at_vertex(v);
                        (at.len() == valence && at.iter().all(|&c| tangent[c] && reaches(c, v)))
                            .then_some((v, at))
                    })
                    .collect()
            };
        let odd_corners = odd_with(false, &|c, v| {
            let far = far_of(c, v);
            border[far] && !apex[far]
        });
        // A curve into a corner held by its spans takes a knot of
        // multiplicity two at the end of its span next to the corner, in
        // its one strip, which must then run one way for every patch and
        // curve along it; past that span it is held by the pole triples,
        // as a curve between corners of four patches is. Its far end may be
        // on the border, a corner of four patches, a junction, or another
        // corner held by its spans: the corners are pared down until every
        // curve's far end is one of them.
        let mut span_at = vec![true; vertices];
        let span_corners = loop {
            let found = odd_with(true, &|c, v| {
                let far = far_of(c, v);
                let ends_held = (border[far] && !apex[far])
                    || interior_valence[far] == Some(4)
                    || (span_at[far] && interior_valence[far].is_some());
                ends_held && curve_strip[c].len() == 1 && !two_way[curve_strip[c][0]]
            });
            let mut kept = vec![false; vertices];
            for (v, _) in &found {
                kept[*v] = true;
            }
            if kept == span_at {
                break found;
            }
            span_at = kept;
        };

        // Chains: held curves running straight on through a corner of four
        // patches, or through a junction of a split side, share their
        // ratio.
        let curve_of = |p: usize, k: usize, last: bool| {
            let side = &patches[p].sides[k];
            if last { side.last() } else { side.first() }.map(|s| s.curve)
        };
        let holding = |held: Vec<bool>| -> Holding {
            let mut chains = Classes::new(curves);
            let mut conflicted_curves = Vec::new();
            for v in 0..vertices {
                if interior_valence[v] != Some(4) || through(v) {
                    continue;
                }
                let at: Vec<(usize, [Option<usize>; 2])> = touching[v]
                    .iter()
                    .filter_map(|&(p, t)| match t {
                        Touch::Corner(k) => {
                            Some((p, [curve_of(p, (k + 3) % 4, true), curve_of(p, k, false)]))
                        }
                        Touch::Through(..) => None,
                    })
                    .collect();
                let Some(ring) = ring_round(&at) else {
                    continue;
                };
                for i in 0..2 {
                    let (c, far) = (ring[i].0, ring[i + 2].0);
                    if !(held[c] && held[far]) {
                        continue;
                    }
                    // The patch between `c` and the curve after it, and the
                    // one between that curve and `far`: the sides of `c` and
                    // `far` facing the curve between them.
                    let near_side = ring[i].1;
                    let far_side = ring[i + 1].1;
                    let flip_c = users[c][0].0 != near_side;
                    let flip_far = users[far][0].0 != far_side;
                    if !chains.join(c, far, flip_c ^ flip_far) {
                        conflicted_curves.push(c);
                    }
                }
            }
            for (p, patch) in patches.iter().enumerate() {
                for side in &patch.sides {
                    for w in side.windows(2) {
                        let (c, next) = (w[0].curve, w[1].curve);
                        let v = vertex(walk_ends(w[0]).1);
                        if interior_valence[v] != Some(4) || !(held[c] && held[next]) {
                            continue;
                        }
                        // Patch `p` lies on the same side of both: their
                        // ratios agree where `p` is the first patch of both
                        // or of neither.
                        let flip = (users[c][0].0 == p) != (users[next][0].0 == p);
                        if !chains.join(c, next, flip) {
                            conflicted_curves.push(c);
                        }
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
            Holding {
                held,
                chain,
                chains: chain_count,
                conflicted,
            }
        };
        // A corner of other than four patches, or an apex, turns the
        // parameter lines that cross a curve one way on one side of it and
        // the other way on the other: no pole triple of the curve can hold
        // across it, but where the corner is held by its spans.
        let held_where = |spanned: &[bool]| -> Vec<bool> {
            (0..curves)
                .map(|c| {
                    let (a, b) = curve_ends[c];
                    let ok = |v: usize| !free_at(v) || spanned[v];
                    tangent[c] && ok(a) && ok(b)
                })
                .collect()
        };
        let plain = holding(held_where(&vec![false; vertices]));
        let mut spanned_at = vec![false; vertices];
        for (v, _) in &span_corners {
            spanned_at[*v] = true;
        }
        let spanned = holding(held_where(&spanned_at));

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
            touching,
            tangent,
            patch_flip,
            curve_flip,
            plain,
            spanned,
            odd_corners,
            span_corners,
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
