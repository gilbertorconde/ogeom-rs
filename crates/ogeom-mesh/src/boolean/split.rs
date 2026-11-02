//! A triangle cut along the segments the other operand left on it.
//!
//! The triangle is triangulated afresh in its own plane: its corners, then
//! every point found on it inserted one at a time, then every segment forced
//! in as a constraint. The work is done in a projection that keeps the
//! plane's area, and every orientation test is exact on the implicit points,
//! so the triangulation is valid however close the points lie. No triangle
//! quality is sought: the pieces are classified and kept or dropped, not
//! analysed.
//!
//! The kernel's planar region engine ([`crate::planar`]) is exact too, but
//! over points it first snaps to a grid. A snapped point leaves the edge
//! and the plane it was found on, so the triangle across that edge, cut on
//! its own, would split the edge somewhere else, and the operand would no
//! longer close. Here no point moves until the answer is decided.

use ogeom_core::{FastMap, FastSet, OgeomResult, Sign, ogeom_bail};

use super::points::{Pt, axes, cmp, orient2d};

/// The pieces of one cut triangle.
#[derive(Debug, Default)]
pub(super) struct Pieces {
    /// Counter-clockwise about the triangle's own normal, as global point
    /// indices.
    pub triangles: Vec<[u32; 3]>,
    /// The constraint edges, split at every point on them.
    pub fixed: Vec<[u32; 2]>,
}

/// Where a point lies in the triangulation.
enum Loc {
    Inside(u32),
    /// On side `k` of the triangle, from its corner `k` to the next.
    Edge(u32, usize),
    Vertex,
}

struct Cdt<'a> {
    pts: &'a [Pt],
    drop: usize,
    /// Whether orientations in the projection are reversed from the
    /// triangle's own.
    flip: bool,
    g: Vec<u32>,
    tris: Vec<[u32; 3]>,
    alive: Vec<bool>,
    edges: FastMap<(u32, u32), u32>,
    fixed: FastSet<(u32, u32)>,
    around: Vec<u32>,
    last: u32,
    seed: u64,
}

const fn key(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

impl Cdt<'_> {
    fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
        let s = orient2d(
            self.drop,
            &self.pts[self.g[a as usize] as usize],
            &self.pts[self.g[b as usize] as usize],
            &self.pts[self.g[c as usize] as usize],
        );
        if self.flip { s.reversed() } else { s }
    }

    fn add(&mut self, t: [u32; 3]) {
        let at = u32::try_from(self.tris.len()).unwrap_or(u32::MAX);
        self.tris.push(t);
        self.alive.push(true);
        for k in 0..3 {
            self.edges.insert((t[k], t[(k + 1) % 3]), at);
            self.around[t[k] as usize] = at;
        }
        self.last = at;
    }

    fn kill(&mut self, at: u32) {
        self.alive[at as usize] = false;
        let t = self.tris[at as usize];
        for k in 0..3 {
            let e = (t[k], t[(k + 1) % 3]);
            if self.edges.get(&e) == Some(&at) {
                self.edges.remove(&e);
            }
        }
    }

    /// The triangle across directed edge `a -> b` of a triangle: the one
    /// holding `b -> a`.
    fn across(&self, a: u32, b: u32) -> Option<u32> {
        self.edges.get(&(b, a)).copied()
    }

    fn third(&self, at: u32, a: u32, b: u32) -> u32 {
        let t = self.tris[at as usize];
        t.into_iter().find(|&v| v != a && v != b).unwrap_or(a)
    }

    fn random(&mut self) -> usize {
        self.seed ^= self.seed << 13;
        self.seed ^= self.seed >> 7;
        self.seed ^= self.seed << 17;
        usize::try_from(self.seed % 3).unwrap_or(0)
    }

    fn classify(at: u32, o: [Sign; 3]) -> Loc {
        let mut zeros = (0..3).filter(|&k| o[k].is_zero());
        match (zeros.next(), zeros.next()) {
            (None, _) => Loc::Inside(at),
            (Some(k), None) => Loc::Edge(at, k),
            _ => Loc::Vertex,
        }
    }

    /// Find point `p` by a walk that steps across an edge `p` lies beyond,
    /// chosen at random among them, so it cannot cycle; a walk that runs
    /// long falls back to a scan.
    fn locate(&mut self, p: u32) -> OgeomResult<Loc> {
        let mut at = self.last;
        let cap = 4 * self.tris.len() + 16;
        'walk: for _ in 0..cap {
            let t = self.tris[at as usize];
            let start = self.random();
            let mut o = [Sign::Zero; 3];
            for i in 0..3 {
                let k = (start + i) % 3;
                o[k] = self.orient(t[k], t[(k + 1) % 3], p);
                if o[k] == Sign::Negative {
                    let Some(next) = self.across(t[k], t[(k + 1) % 3]) else {
                        ogeom_bail!(Invariant, "a point found on a triangle lies outside it");
                    };
                    at = next;
                    continue 'walk;
                }
            }
            return Ok(Self::classify(at, o));
        }
        for at in 0..self.tris.len() {
            if !self.alive[at] {
                continue;
            }
            let t = self.tris[at];
            let o = [0, 1, 2].map(|k| self.orient(t[k], t[(k + 1) % 3], p));
            if !o.contains(&Sign::Negative) {
                return Ok(Self::classify(u32::try_from(at).unwrap_or(u32::MAX), o));
            }
        }
        ogeom_bail!(Invariant, "a point found on a triangle lies outside it")
    }

    fn insert_point(&mut self, p: u32) -> OgeomResult<()> {
        match self.locate(p)? {
            Loc::Inside(at) => {
                let [a, b, c] = self.tris[at as usize];
                self.kill(at);
                self.add([a, b, p]);
                self.add([b, c, p]);
                self.add([c, a, p]);
            }
            Loc::Edge(at, k) => {
                let t = self.tris[at as usize];
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let twin = self.across(a, b);
                self.kill(at);
                self.add([a, p, c]);
                self.add([p, b, c]);
                if let Some(tw) = twin {
                    let d = self.third(tw, a, b);
                    self.kill(tw);
                    self.add([b, p, d]);
                    self.add([p, a, d]);
                }
                if self.fixed.remove(&key(a, b)) {
                    self.fixed.insert(key(a, p));
                    self.fixed.insert(key(p, b));
                }
            }
            Loc::Vertex => ogeom_bail!(Invariant, "two points found on a triangle coincide"),
        }
        Ok(())
    }

    /// The live triangles about vertex `u`, each turned to start at `u`.
    fn fan(&self, u: u32) -> Vec<[u32; 3]> {
        let turn = |at: u32| {
            let t = self.tris[at as usize];
            let k = t.iter().position(|&v| v == u).unwrap_or(0);
            [t[k], t[(k + 1) % 3], t[(k + 2) % 3]]
        };
        let first = self.around[u as usize];
        let mut out = vec![turn(first)];
        let mut at = first;
        let mut closed = false;
        while let Some(next) = self.across(turn(at)[2], u) {
            if next == first {
                closed = true;
                break;
            }
            out.push(turn(next));
            at = next;
        }
        if !closed {
            at = first;
            while let Some(prev) = self.across(u, turn(at)[1]) {
                out.push(turn(prev));
                at = prev;
            }
        }
        out
    }

    /// Whether `x`, collinear with `u` and `v`, lies strictly between them.
    fn between(&self, u: u32, x: u32, v: u32) -> bool {
        let pt = |i: u32| &self.pts[self.g[i as usize] as usize];
        let (i, j) = axes(self.drop);
        let Some(axis) = [i, j].into_iter().find(|&k| cmp(k, pt(u), pt(v)).is_ne()) else {
            return false;
        };
        let way = cmp(axis, pt(u), pt(v));
        cmp(axis, pt(u), pt(x)) == way && cmp(axis, pt(x), pt(v)) == way
    }

    fn insert_segment(&mut self, mut u: u32, v: u32) -> OgeomResult<()> {
        'piece: while u != v {
            if self.edges.contains_key(&(u, v)) || self.edges.contains_key(&(v, u)) {
                self.fixed.insert(key(u, v));
                return Ok(());
            }
            let mut start = None;
            for [_, x, y] in self.fan(u) {
                let (ox, oy) = (self.orient(u, v, x), self.orient(u, v, y));
                for (o, w) in [(ox, x), (oy, y)] {
                    if o.is_zero() && self.between(u, w, v) {
                        self.fixed.insert(key(u, w));
                        u = w;
                        continue 'piece;
                    }
                }
                if ox == Sign::Negative && oy == Sign::Positive {
                    start = Some((x, y));
                }
            }
            let Some((x, y)) = start else {
                ogeom_bail!(Invariant, "a segment on a triangle leaves it");
            };
            let first = self.edges.get(&(x, y)).copied();
            let mut crossed: Vec<u32> = first.into_iter().collect();
            let (mut r, mut l) = (x, y);
            let mut right = vec![x];
            let mut left = vec![y];
            let end = loop {
                if self.fixed.contains(&key(r, l)) {
                    ogeom_bail!(Invariant, "two segments on a triangle cross");
                }
                let Some(next) = self.across(r, l) else {
                    ogeom_bail!(Invariant, "a segment on a triangle leaves it");
                };
                crossed.push(next);
                let z = self.third(next, r, l);
                if z == v {
                    break v;
                }
                match self.orient(u, v, z) {
                    Sign::Zero => break z,
                    Sign::Positive => {
                        left.push(z);
                        l = z;
                    }
                    Sign::Negative => {
                        right.push(z);
                        r = z;
                    }
                }
            };
            for at in crossed {
                self.kill(at);
            }
            let mut upper = vec![u, end];
            upper.extend(left.iter().rev());
            let mut lower = vec![end, u];
            lower.extend(right.iter());
            self.fill(upper)?;
            self.fill(lower)?;
            self.fixed.insert(key(u, end));
            u = end;
        }
        Ok(())
    }

    /// Triangulate a counter-clockwise simple polygon by cutting ears: a
    /// strictly convex corner whose triangle holds no other corner, not
    /// even on its sides.
    fn fill(&mut self, mut poly: Vec<u32>) -> OgeomResult<()> {
        while poly.len() > 3 {
            let n = poly.len();
            let ear = (0..n).find(|&i| {
                let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
                self.orient(a, b, c) == Sign::Positive
                    && poly.iter().all(|&w| {
                        w == a
                            || w == b
                            || w == c
                            || self.orient(a, b, w) == Sign::Negative
                            || self.orient(b, c, w) == Sign::Negative
                            || self.orient(c, a, w) == Sign::Negative
                    })
            });
            let Some(i) = ear else {
                ogeom_bail!(Invariant, "a cavity about a segment has no ear");
            };
            let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
            self.add([a, b, c]);
            poly.remove(i);
        }
        if let [a, b, c] = poly[..] {
            if self.orient(a, b, c) != Sign::Positive {
                ogeom_bail!(Invariant, "a cavity about a segment folds over");
            }
            self.add([a, b, c]);
        }
        Ok(())
    }
}

/// Cut the triangle with `corners` at the points `inner` and along
/// `segments`, all global indices into `pts`. `drop` is a coordinate the
/// triangle's plane is not parallel to.
///
/// # Errors
///
/// [`ogeom_core::OgeomError::Invariant`] if the points and segments are not
/// a planar arrangement inside the triangle: a point outside it, two
/// coincident points, or two segments crossing. The callers find them so
/// that none of these can happen.
pub(super) fn split(
    pts: &[Pt],
    drop: usize,
    corners: [u32; 3],
    inner: &[u32],
    segments: &[[u32; 2]],
) -> OgeomResult<Pieces> {
    let mut g: Vec<u32> = corners.to_vec();
    let mut local: FastMap<u32, u32> = FastMap::default();
    for (k, &c) in corners.iter().enumerate() {
        local.insert(c, u32::try_from(k).unwrap_or(0));
    }
    for &p in inner {
        if !local.contains_key(&p) {
            local.insert(p, u32::try_from(g.len()).unwrap_or(u32::MAX));
            g.push(p);
        }
    }
    let n = g.len();
    let mut cdt = Cdt {
        pts,
        drop,
        flip: false,
        g,
        tris: Vec::with_capacity(2 * n + 1),
        alive: Vec::with_capacity(2 * n + 1),
        edges: FastMap::default(),
        fixed: FastSet::default(),
        around: vec![0; n],
        last: 0,
        seed: 0x9e37_79b9_7f4a_7c15,
    };
    cdt.flip = cdt.orient(0, 1, 2) == Sign::Negative;
    if cdt.orient(0, 1, 2).is_zero() {
        ogeom_bail!(Invariant, "a triangle to cut is flat in its projection");
    }
    cdt.add([0, 1, 2]);
    for p in 3..n {
        cdt.insert_point(u32::try_from(p).unwrap_or(u32::MAX))?;
    }
    for s in segments {
        let (Some(&a), Some(&b)) = (local.get(&s[0]), local.get(&s[1])) else {
            ogeom_bail!(Invariant, "a segment's end is not a point on its triangle");
        };
        cdt.insert_segment(a, b)?;
    }
    let mut out = Pieces::default();
    for (t, &alive) in cdt.tris.iter().zip(&cdt.alive) {
        if alive {
            out.triangles.push(t.map(|v| cdt.g[v as usize]));
        }
    }
    let mut fixed: Vec<[u32; 2]> = cdt
        .fixed
        .iter()
        .map(|&(a, b)| {
            let (a, b) = (cdt.g[a as usize], cdt.g[b as usize]);
            if a < b { [a, b] } else { [b, a] }
        })
        .collect();
    fixed.sort_unstable();
    out.fixed = fixed;
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::points::Def;
    use super::*;

    fn pts(list: &[[f64; 3]]) -> Vec<Pt> {
        list.iter().map(|&p| Pt::new(Def::Vertex(p))).collect()
    }

    fn area(p: &[Pt], t: [u32; 3]) -> f64 {
        let at = |i: u32| match p[i as usize].def {
            Def::Vertex(q) => q,
            _ => [0.0; 3],
        };
        let [a, b, c] = t.map(at);
        ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.0
    }

    #[test]
    fn points_and_a_segment_split_a_triangle_exactly() {
        let p = pts(&[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [2.0, 1.0, 0.0],
        ]);
        let cut = split(&p, 2, [0, 1, 2], &[3, 4, 5, 6], &[[3, 6], [5, 4]]).unwrap();
        let total: f64 = cut.triangles.iter().map(|&t| area(&p, t)).sum();
        assert_eq!(total, 8.0);
        assert!(cut.triangles.iter().all(|&t| area(&p, t) > 0.0));
        // The segment 5 -> 4 runs through point 6, and is split there.
        assert!(cut.fixed.contains(&[5, 6]) && cut.fixed.contains(&[4, 6]));
        assert!(cut.fixed.contains(&[3, 6]));
        for [a, b] in &cut.fixed {
            let has = |t: &[u32; 3]| t.contains(a) && t.contains(b);
            assert_eq!(cut.triangles.iter().filter(|t| has(t)).count(), 2);
        }
    }

    #[test]
    fn clockwise_projection_keeps_the_triangles_winding() {
        let p = pts(&[
            [0.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [4.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ]);
        let cut = split(&p, 2, [0, 1, 2], &[3], &[]).unwrap();
        assert_eq!(cut.triangles.len(), 3);
        assert!(cut.triangles.iter().all(|&t| area(&p, t) < 0.0));
    }
}
