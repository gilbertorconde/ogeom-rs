//! Where two triangles meet, exactly: the points and segments they share.
//!
//! Two triangles not in one plane meet in a segment, a point, or nothing,
//! all on the line their planes meet in. Each triangle's own crossing of the
//! other's plane is found first, from the signs of its corners against that
//! plane, and the two crossings are then overlapped along the line.
//!
//! Two triangles in one plane meet in a convex polygon. Its corners are the
//! corners of each triangle inside the other and the points where their
//! edges cross; its sides are the parts of each triangle's edges inside the
//! other.
//!
//! Every point is named by a [`Key`] built from the input indices that
//! define it, so the same point found from two pairs of triangles has the
//! same name.

use std::cmp::Ordering;

use ogeom_core::{Exact, Predicates, Sign};

use super::points::{self, Def, P3, Pt, axes, by_weight, cmp, normal_sign, same};

/// Which operand.
pub(super) type Side = u8;

/// The name of a point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(super) enum Key {
    /// Vertex `v` of operand `side`.
    Vertex(Side, u32),
    /// Where the edge `edge` of operand `side`, lower index first, crosses
    /// the plane of triangle `tri` of the other operand.
    Lpi {
        side: Side,
        edge: [u32; 2],
        tri: u32,
    },
    /// Where edge `a` of the first operand crosses edge `b` of the second,
    /// in their common plane.
    Ssi { a: [u32; 2], b: [u32; 2] },
}

/// A triangle of one operand.
#[derive(Debug, Clone, Copy)]
pub(super) struct Tri {
    pub side: Side,
    pub id: u32,
    pub v: [u32; 3],
    pub p: [P3; 3],
}

impl Tri {
    fn corner(&self, k: usize) -> (Key, Def) {
        (Key::Vertex(self.side, self.v[k]), Def::Vertex(self.p[k]))
    }

    /// Edge `k`, from corner `k` to the next, lower vertex index first.
    fn edge(&self, k: usize) -> ([u32; 2], [P3; 2]) {
        let (i, j) = (k, (k + 1) % 3);
        if self.v[i] < self.v[j] {
            ([self.v[i], self.v[j]], [self.p[i], self.p[j]])
        } else {
            ([self.v[j], self.v[i]], [self.p[j], self.p[i]])
        }
    }
}

/// What two triangles share.
#[derive(Debug, Default)]
pub(super) struct Touch {
    pub points: Vec<(Key, Pt)>,
    /// Pairs of indices into `points`.
    pub segments: Vec<[usize; 2]>,
    /// Whether the triangles lie in one plane.
    pub coplanar: bool,
}

impl Touch {
    fn add(&mut self, key: Key, def: Def) -> usize {
        if let Some(at) = self.points.iter().position(|(k, _)| *k == key) {
            return at;
        }
        self.points.push((key, Pt::new(def)));
        self.points.len() - 1
    }

    fn segment(&mut self, a: usize, b: usize) {
        if a != b && !self.segments.contains(&[a, b]) && !self.segments.contains(&[b, a]) {
            self.segments.push([a, b]);
        }
    }
}

/// What triangles `a` and `b` share, or `None` when they do not meet.
pub(super) fn meet(a: &Tri, b: &Tri) -> Option<Touch> {
    let sb = b.p.map(|q| Exact::orient3d(a.p[0], a.p[1], a.p[2], q));
    if strictly_one_side(sb) {
        return None;
    }
    if sb.iter().all(|s| s.is_zero()) {
        return coplanar(a, b);
    }
    let sa = a.p.map(|q| Exact::orient3d(b.p[0], b.p[1], b.p[2], q));
    if strictly_one_side(sa) {
        return None;
    }
    crossing(a, sa, b, sb)
}

fn strictly_one_side(s: [Sign; 3]) -> bool {
    s.iter().all(|&x| x == Sign::Positive) || s.iter().all(|&x| x == Sign::Negative)
}

/// The points where triangle `t` meets the plane of `other`, from its
/// corners' signs against that plane: corners on it, and edges crossing it.
fn cut(t: &Tri, s: [Sign; 3], other: &Tri) -> Vec<(Key, Def)> {
    let mut out = Vec::with_capacity(2);
    for k in 0..3 {
        if s[k].is_zero() {
            out.push(t.corner(k));
        }
        let next = s[(k + 1) % 3];
        if !s[k].is_zero() && s[k] == next.reversed() {
            let (edge, line) = t.edge(k);
            out.push((
                Key::Lpi {
                    side: t.side,
                    edge,
                    tri: other.id,
                },
                Def::Lpi {
                    line,
                    plane: other.p,
                },
            ));
        }
    }
    out
}

/// Two triangles in different planes.
fn crossing(a: &Tri, sa: [Sign; 3], b: &Tri, sb: [Sign; 3]) -> Option<Touch> {
    let on_a: Vec<(Key, Pt)> = cut(a, sa, b)
        .into_iter()
        .map(|(k, d)| (k, Pt::new(d)))
        .collect();
    let on_b: Vec<(Key, Pt)> = cut(b, sb, a)
        .into_iter()
        .map(|(k, d)| (k, Pt::new(d)))
        .collect();
    if on_a.is_empty() || on_b.is_empty() {
        return None;
    }
    // Both crossings lie on the line the planes meet in; order along an
    // axis that line is not square to.
    let axis = by_weight(points::crossing_approx(&a.p, &b.p))
        .into_iter()
        .find(|&k| !points::crossing_sign(&a.p, &b.p, k).is_zero())?;
    let ends = |set: &[(Key, Pt)]| -> (usize, usize) {
        if set.len() == 2 && cmp(axis, &set[1].1, &set[0].1) == Ordering::Less {
            (1, 0)
        } else {
            (0, set.len() - 1)
        }
    };
    let (alo, ahi) = ends(&on_a);
    let (blo, bhi) = ends(&on_b);
    // The later start and the earlier end; between equals, a vertex.
    let pick = |x: &(Key, Pt), y: &(Key, Pt), later: bool| -> Option<(Key, Def)> {
        let o = cmp(axis, &x.1, &y.1);
        let first_wins = match o {
            Ordering::Equal => matches!(x.0, Key::Vertex(..)) || !matches!(y.0, Key::Vertex(..)),
            Ordering::Greater => later,
            Ordering::Less => !later,
        };
        Some(if first_wins {
            (x.0, x.1.def)
        } else {
            (y.0, y.1.def)
        })
    };
    let lo = pick(&on_a[alo], &on_b[blo], true)?;
    let hi = pick(&on_a[ahi], &on_b[bhi], false)?;
    let (lo_pt, hi_pt) = (Pt::new(lo.1), Pt::new(hi.1));
    let mut touch = Touch::default();
    match cmp(axis, &lo_pt, &hi_pt) {
        Ordering::Greater => return None,
        Ordering::Equal => {
            touch.add(lo.0, lo.1);
        }
        Ordering::Less => {
            let i = touch.add(lo.0, lo.1);
            let j = touch.add(hi.0, hi.1);
            touch.segment(i, j);
        }
    }
    Some(touch)
}

/// Two triangles in one plane.
fn coplanar(a: &Tri, b: &Tri) -> Option<Touch> {
    let drop = by_weight(points::normal_approx(&a.p))
        .into_iter()
        .find(|&k| !normal_sign(&a.p, k).is_zero())?;
    let (i, j) = axes(drop);
    let flat = |t: &Tri| t.p.map(|q| [q[i], q[j]]);
    let (fa, fb) = (flat(a), flat(b));
    let oa = Exact::orient2d(fa[0], fa[1], fa[2]);
    let ob = Exact::orient2d(fb[0], fb[1], fb[2]);
    let inside = |f: &[[f64; 2]; 3], o: Sign, q: [f64; 2]| {
        (0..3).all(|k| Exact::orient2d(f[k], f[(k + 1) % 3], q) != o.reversed())
    };
    let mut touch = Touch {
        coplanar: true,
        ..Touch::default()
    };
    // For each edge of either triangle, the points found on it.
    let mut on_edge: [[Vec<usize>; 3]; 2] = Default::default();
    for (t, ft, f, o, at) in [(a, &fa, &fb, ob, 0), (b, &fb, &fa, oa, 1)] {
        for k in 0..3 {
            if inside(f, o, ft[k]) {
                let (key, def) = t.corner(k);
                let n = touch.add(key, def);
                on_edge[at][k].push(n);
                on_edge[at][(k + 2) % 3].push(n);
                // On the other triangle's edges too, where it lies on them.
                for e in 0..3 {
                    if Exact::orient2d(f[e], f[(e + 1) % 3], ft[k]).is_zero() {
                        on_edge[1 - at][e].push(n);
                    }
                }
            }
        }
    }
    for ka in 0..3 {
        for kb in 0..3 {
            let (pa, qa) = (fa[ka], fa[(ka + 1) % 3]);
            let (pb, qb) = (fb[kb], fb[(kb + 1) % 3]);
            let (s1, s2) = (Exact::orient2d(pa, qa, pb), Exact::orient2d(pa, qa, qb));
            let (s3, s4) = (Exact::orient2d(pb, qb, pa), Exact::orient2d(pb, qb, qa));
            let proper =
                !s1.is_zero() && s1 == s2.reversed() && !s3.is_zero() && s3 == s4.reversed();
            if proper {
                let (ea, la) = a.edge(ka);
                let (eb, lb) = b.edge(kb);
                let n = touch.add(
                    Key::Ssi { a: ea, b: eb },
                    Def::Ssi {
                        first: la,
                        second: lb,
                        drop,
                    },
                );
                on_edge[0][ka].push(n);
                on_edge[1][kb].push(n);
            }
        }
    }
    if touch.points.is_empty() {
        return None;
    }
    // Each edge's part inside the other triangle runs between the extreme
    // points found on it.
    for (at, t) in [(0, a), (1, b)] {
        for (k, found) in on_edge[at].iter().enumerate() {
            let (p, q) = (t.p[k], t.p[(k + 1) % 3]);
            let Some(axis) = [i, j].into_iter().find(|&x| p[x] != q[x]) else {
                continue;
            };
            let order =
                |x: &&usize, y: &&usize| cmp(axis, &touch.points[**x].1, &touch.points[**y].1);
            let (Some(&lo), Some(&hi)) = (found.iter().min_by(order), found.iter().max_by(order))
            else {
                continue;
            };
            if !same(&touch.points[lo].1, &touch.points[hi].1) {
                touch.segment(lo, hi);
            }
        }
    }
    Some(touch)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn tri(side: Side, v: [u32; 3], p: [P3; 3]) -> Tri {
        Tri { side, id: 0, v, p }
    }

    #[test]
    fn crossing_triangles_share_a_segment_between_edge_crossings() {
        let a = tri(
            0,
            [0, 1, 2],
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
        );
        let b = tri(
            1,
            [0, 1, 2],
            [[1.0, 1.0, -1.0], [1.0, 1.0, 1.0], [3.0, -1.0, 0.0]],
        );
        let t = meet(&a, &b).unwrap();
        assert!(!t.coplanar);
        assert_eq!(t.segments.len(), 1);
        let keys: Vec<Key> = t.points.iter().map(|p| p.0).collect();
        // One end where b's edge 0-1 pierces a; the other where b leaves a
        // across a's edge 0-1.
        assert!(keys.contains(&Key::Lpi {
            side: 1,
            edge: [0, 1],
            tri: 0
        }));
        assert!(keys.contains(&Key::Lpi {
            side: 0,
            edge: [0, 1],
            tri: 0
        }));
    }

    #[test]
    fn apart_and_touching_at_a_corner() {
        let a = tri(
            0,
            [0, 1, 2],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        );
        let far = tri(
            1,
            [0, 1, 2],
            [[2.0, 2.0, -1.0], [2.0, 2.0, 1.0], [3.0, 2.0, 0.0]],
        );
        assert!(meet(&a, &far).is_none());
        let corner = tri(
            1,
            [5, 6, 7],
            [[0.0, 0.0, 0.0], [-1.0, 0.0, 1.0], [0.0, -1.0, 1.0]],
        );
        let t = meet(&a, &corner).unwrap();
        assert!(t.segments.is_empty());
        assert_eq!(t.points.len(), 1);
    }

    #[test]
    fn coplanar_triangles_share_the_sides_of_their_overlap() {
        let a = tri(
            0,
            [0, 1, 2],
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        );
        let b = tri(
            1,
            [0, 1, 2],
            [[1.0, -1.0, 0.0], [1.0, 3.0, 0.0], [-1.0, 1.0, 0.0]],
        );
        let t = meet(&a, &b).unwrap();
        assert!(t.coplanar);
        // The overlap is a quadrilateral with a corner at a's corner 0.
        assert_eq!(t.points.len(), 4);
        assert_eq!(t.segments.len(), 4);
        let shared = tri(
            0,
            [0, 2, 9],
            [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [-1.0, 1.0, 0.0]],
        );
        let t = meet(&a, &shared).unwrap();
        assert_eq!(t.points.len(), 2);
        assert_eq!(t.segments.len(), 1);
    }
}
