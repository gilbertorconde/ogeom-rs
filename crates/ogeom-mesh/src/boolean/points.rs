//! Points held implicitly, by the input coordinates that define them, and
//! the exact predicates over them.
//!
//! A point where an edge of one mesh crosses a triangle of the other has
//! rational coordinates that `f64` cannot hold. Rounding it would move it off
//! the edge or off the triangle's plane, and the next orientation test made
//! on the rounded point could contradict the one that created it. So the
//! point is kept as the inputs that define it, and every test on it is the
//! sign of a polynomial in those inputs, decided exactly.
//!
//! Each point is written in homogeneous form `(X, Y, Z, W)`, the point
//! being `(X/W, Y/W, Z/W)` with `W` not zero, so that no test divides.

use std::cmp::Ordering;
use std::sync::OnceLock;

use ogeom_core::{Exact, Predicates, Sign};

use super::exact::{Big, Iv, Ring, times};

/// A point in space, as three coordinates.
pub(super) type P3 = [f64; 3];

/// How a point is defined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Def {
    /// An input vertex.
    Vertex(P3),
    /// Where the line through `line` crosses the plane through `plane`,
    /// the line's two points lying strictly on opposite sides of it.
    Lpi { line: [P3; 2], plane: [P3; 3] },
    /// Where two coplanar lines cross, found in the projection that drops
    /// coordinate `drop`, which must keep their common plane's area.
    Ssi {
        first: [P3; 2],
        second: [P3; 2],
        drop: usize,
    },
}

fn diff<N: Ring>(a: &[N; 3], b: &[N; 3]) -> [N; 3] {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

fn cross<N: Ring>(a: &[N; 3], b: &[N; 3]) -> [N; 3] {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

fn dot<N: Ring>(a: &[N; 3], b: &[N; 3]) -> N {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn lift<N: Ring>(p: &P3) -> [N; 3] {
    [N::of(p[0]), N::of(p[1]), N::of(p[2])]
}

/// The normal `(b - a) × (c - a)` of a triangle.
fn normal<N: Ring>(t: &[P3; 3]) -> [N; 3] {
    let [a, b, c] = t.map(|p| lift::<N>(&p));
    cross(&diff(&b, &a), &diff(&c, &a))
}

impl Def {
    fn hom<N: Ring>(&self) -> [N; 4] {
        match self {
            Self::Vertex(p) => [N::of(p[0]), N::of(p[1]), N::of(p[2]), N::of(1.0)],
            Self::Lpi { line, plane } => {
                let n = normal::<N>(plane);
                let r = lift::<N>(&plane[0]);
                let (p, q) = (lift::<N>(&line[0]), lift::<N>(&line[1]));
                let d1 = dot(&n, &diff(&p, &r));
                let d2 = dot(&n, &diff(&q, &r));
                let x = |k: usize| p[k].mul(&d2).sub(&q[k].mul(&d1));
                [x(0), x(1), x(2), d2.sub(&d1)]
            }
            Self::Ssi {
                first,
                second,
                drop,
            } => {
                let (i, j) = axes(*drop);
                let (p, q) = (lift::<N>(&first[0]), lift::<N>(&first[1]));
                let (r, s) = (lift::<N>(&second[0]), lift::<N>(&second[1]));
                let a = diff(&q, &p);
                let b = diff(&s, &r);
                let rp = diff(&r, &p);
                let den = a[i].mul(&b[j]).sub(&a[j].mul(&b[i]));
                let num = rp[i].mul(&b[j]).sub(&rp[j].mul(&b[i]));
                let x = |k: usize| p[k].mul(&den).add(&a[k].mul(&num));
                [x(0), x(1), x(2), den]
            }
        }
    }
}

/// The two coordinates kept when `drop` is dropped, in the cyclic order
/// that makes a triangle's orientation there the sign of its normal's
/// `drop` component.
pub(super) const fn axes(drop: usize) -> (usize, usize) {
    ((drop + 1) % 3, (drop + 2) % 3)
}

/// A point in homogeneous form, enclosed in intervals and exactly.
pub(super) trait Hom {
    fn iv(&self) -> [Iv; 4];
    fn big(&self) -> [Big; 4];
}

/// A defined point, with its enclosure computed once and its exact form
/// computed when first asked for.
#[derive(Debug)]
pub(super) struct Pt {
    pub def: Def,
    iv: [Iv; 4],
    exact: OnceLock<[Big; 4]>,
}

impl Pt {
    pub fn new(def: Def) -> Self {
        Self {
            def,
            iv: def.hom::<Iv>(),
            exact: OnceLock::new(),
        }
    }

    /// The coordinates as `f64`, nearest the exact point, and a bound on
    /// the distance from them to it.
    pub fn rounded(&self) -> (P3, f64) {
        if let Def::Vertex(p) = self.def {
            return (p, 0.0);
        }
        let quotient = |h: &[Iv; 4]| -> Option<[Iv; 3]> {
            Some([h[0].over(&h[3])?, h[1].over(&h[3])?, h[2].over(&h[3])?])
        };
        // The enclosure carried through the predicates is tight unless the
        // point is nearly undefined, as where a line nearly parallels the
        // plane it crosses; then the exact coordinates are enclosed afresh.
        let tight = |q: &[Iv; 3]| {
            q.iter()
                .all(|v| v.hi <= v.lo.next_up().next_up().next_up().next_up())
        };
        let q = quotient(&self.iv)
            .filter(tight)
            .or_else(|| quotient(&self.big().map(|b| b.enclose())));
        let Some(q) = q else {
            return ([f64::NAN; 3], f64::INFINITY);
        };
        let mut p = [0.0; 3];
        let mut off = 0.0_f64;
        for k in 0..3 {
            let (mid, err) = q[k].middle();
            p[k] = mid;
            off = err.mul_add(err, off);
        }
        (
            p,
            if off == 0.0 {
                0.0
            } else {
                off.sqrt().next_up()
            },
        )
    }

    /// A box certainly holding the point, or `None` where the enclosure is
    /// too loose to bound it.
    pub fn bounds(&self) -> Option<[Iv; 3]> {
        let h = &self.iv;
        Some([h[0].over(&h[3])?, h[1].over(&h[3])?, h[2].over(&h[3])?])
    }
}

impl Hom for Pt {
    fn iv(&self) -> [Iv; 4] {
        self.iv
    }

    fn big(&self) -> [Big; 4] {
        if let Def::Vertex(_) = self.def {
            return self.def.hom::<Big>();
        }
        self.exact.get_or_init(|| self.def.hom::<Big>()).clone()
    }
}

/// The centroid of three points: a point certainly inside the triangle
/// they span, when they are not collinear.
#[derive(Debug, Clone, Copy)]
pub(super) struct Mid<'a>(pub [&'a Pt; 3]);

fn centroid<N: Ring>(h: [[N; 4]; 3]) -> [N; 4] {
    let w12 = h[1][3].mul(&h[2][3]);
    let w02 = h[0][3].mul(&h[2][3]);
    let w01 = h[0][3].mul(&h[1][3]);
    let x = |k: usize| {
        h[0][k]
            .mul(&w12)
            .add(&h[1][k].mul(&w02))
            .add(&h[2][k].mul(&w01))
    };
    [x(0), x(1), x(2), N::of(3.0).mul(&w01).mul(&h[2][3])]
}

impl Hom for Mid<'_> {
    fn iv(&self) -> [Iv; 4] {
        centroid(self.0.map(Hom::iv))
    }

    fn big(&self) -> [Big; 4] {
        centroid(self.0.map(Hom::big))
    }
}

/// A number type a homogeneous point can be read in.
pub(super) trait Read: Ring {
    fn read(p: &(impl Hom + ?Sized)) -> [Self; 4];
}

impl Read for Iv {
    fn read(p: &(impl Hom + ?Sized)) -> [Self; 4] {
        p.iv()
    }
}

impl Read for Big {
    fn read(p: &(impl Hom + ?Sized)) -> [Self; 4] {
        p.big()
    }
}

/// Decide a sign: in intervals where they suffice, exactly where not.
macro_rules! decide {
    ($f:ident($($arg:expr),*)) => {
        match $f::<Iv>($($arg),*) {
            Some(s) => s,
            None => $f::<Big>($($arg),*).unwrap_or(Sign::Zero),
        }
    };
}

fn orient2d_in<N: Read>(drop: usize, a: &impl Hom, b: &impl Hom, c: &impl Hom) -> Option<Sign> {
    let (i, j) = axes(drop);
    let [a, b, c] = [N::read(a), N::read(b), N::read(c)];
    let minor = |p: &[N; 4], q: &[N; 4], u: usize, v: usize| p[u].mul(&q[v]).sub(&p[v].mul(&q[u]));
    let det = a[i]
        .mul(&minor(&b, &c, j, 3))
        .sub(&a[j].mul(&minor(&b, &c, i, 3)))
        .add(&a[3].mul(&minor(&b, &c, i, j)));
    let w = times(times(a[3].sign()?, b[3].sign()?), c[3].sign()?);
    Some(times(det.sign()?, w))
}

/// The orientation of three points projected by dropping coordinate
/// `drop`: positive when counter-clockwise in the kept coordinates taken
/// in [`axes`] order.
pub(super) fn orient2d(drop: usize, a: &impl Hom, b: &impl Hom, c: &impl Hom) -> Sign {
    decide!(orient2d_in(drop, a, b, c))
}

fn side_in<N: Read>(t: &[P3; 3], p: &impl Hom) -> Option<Sign> {
    let n = normal::<N>(t);
    let h = N::read(p);
    let a = lift::<N>(&t[0]);
    let off = [
        h[0].sub(&a[0].mul(&h[3])),
        h[1].sub(&a[1].mul(&h[3])),
        h[2].sub(&a[2].mul(&h[3])),
    ];
    Some(times(dot(&n, &off).sign()?, h[3].sign()?))
}

/// Which side of a triangle's plane a point lies on: the sign of
/// `n · (p - t[0])`, `n` the triangle's counter-clockwise normal.
pub(super) fn side(t: &[P3; 3], p: &impl Hom) -> Sign {
    if let Some(s) = side_in::<Iv>(t, p) {
        return s;
    }
    side_in::<Big>(t, p).unwrap_or(Sign::Zero)
}

fn cmp_in<N: Read>(axis: usize, p: &impl Hom, q: &impl Hom) -> Option<Sign> {
    let (p, q) = (N::read(p), N::read(q));
    let d = p[axis].mul(&q[3]).sub(&q[axis].mul(&p[3]));
    Some(times(d.sign()?, times(p[3].sign()?, q[3].sign()?)))
}

/// The order of two points' coordinates along `axis`.
pub(super) fn cmp(axis: usize, p: &impl Hom, q: &impl Hom) -> Ordering {
    match decide!(cmp_in(axis, p, q)) {
        Sign::Negative => Ordering::Less,
        Sign::Zero => Ordering::Equal,
        Sign::Positive => Ordering::Greater,
    }
}

/// Whether two points coincide.
pub(super) fn same(p: &impl Hom, q: &impl Hom) -> bool {
    (0..3).all(|k| cmp(k, p, q) == Ordering::Equal)
}

/// The sign of a triangle's normal along `axis`, exactly.
pub(super) fn normal_sign(t: &[P3; 3], axis: usize) -> Sign {
    let (i, j) = axes(axis);
    let [a, b, c] = t.map(|p| [p[i], p[j]]);
    Exact::orient2d(a, b, c)
}

fn crossing_in<N: Read>(a: &[P3; 3], b: &[P3; 3], axis: usize) -> Option<Sign> {
    cross(&normal::<N>(a), &normal::<N>(b))[axis].sign()
}

/// The sign along `axis` of the direction of the line two triangles'
/// planes meet in, `n_a × n_b`.
pub(super) fn crossing_sign(a: &[P3; 3], b: &[P3; 3], axis: usize) -> Sign {
    if let Some(s) = crossing_in::<Iv>(a, b, axis) {
        return s;
    }
    crossing_in::<Big>(a, b, axis).unwrap_or(Sign::Zero)
}

/// An approximate magnitude of `n_a × n_b` along each axis, to pick the
/// axis the exact tests are made along.
pub(super) fn crossing_approx(a: &[P3; 3], b: &[P3; 3]) -> [f64; 3] {
    cross(&normal::<Iv>(a), &normal::<Iv>(b)).map(|v| v.lo.abs().max(v.hi.abs()))
}

/// An approximate magnitude of a triangle's normal along each axis.
pub(super) fn normal_approx(t: &[P3; 3]) -> [f64; 3] {
    normal::<Iv>(t).map(|v| v.lo.abs().max(v.hi.abs()))
}

/// The axes in decreasing order of `weights`, ties to the lower axis.
pub(super) fn by_weight(weights: [f64; 3]) -> [usize; 3] {
    let mut order = [0, 1, 2];
    order.sort_by(|&x, &y| weights[y].total_cmp(&weights[x]).then(x.cmp(&y)));
    order
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn v(p: P3) -> Pt {
        Pt::new(Def::Vertex(p))
    }

    #[test]
    fn homogeneous_predicates_agree_with_the_explicit_ones() {
        let t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        for d in [[0.2, 0.3, 1.0], [0.2, 0.3, -1.0], [5.0, -3.0, 0.0]] {
            let robust = Exact::orient3d(t[0], t[1], t[2], d);
            assert_eq!(side(&t, &v(d)), robust.reversed());
        }
        for (a, b, c) in [
            ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
        ] {
            for drop in 0..3 {
                let (i, j) = axes(drop);
                let want = Exact::orient2d([a[i], a[j]], [b[i], b[j]], [c[i], c[j]]);
                assert_eq!(orient2d(drop, &v(a), &v(b), &v(c)), want);
            }
        }
    }

    #[test]
    fn a_crossing_lies_on_its_line_and_its_plane() {
        // The edge from (0.1, 0.2, -1) to (0.3, 0.1, 2) crosses z = 0 at a
        // point no f64 holds; it is on the plane exactly, and on the line.
        let plane = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let line = [[0.1, 0.2, -1.0], [0.3, 0.1, 2.0]];
        let x = Pt::new(Def::Lpi { line, plane });
        assert_eq!(side(&plane, &x), Sign::Zero);
        let [p, q] = line.map(v);
        for drop in 0..3 {
            assert_eq!(orient2d(drop, &p, &q, &x), Sign::Zero);
        }
        let (r, err) = x.rounded();
        assert!(err < 1e-15 && r[2].abs() <= err);
        // Through a vertex of a tilted plane, the same point again.
        let tilted = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.3], [0.0, 1.0, 0.7]];
        let y = Pt::new(Def::Lpi {
            line: [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]],
            plane: tilted,
        });
        assert!(same(&y, &v([0.0, 0.0, 0.0])));
        assert!(!same(&x, &v([0.0, 0.0, 0.0])));
    }

    #[test]
    fn coplanar_lines_cross_where_both_run() {
        let first = [[0.0, 0.0, 1.0], [1.0, 0.3, 1.0]];
        let second = [[0.7, -0.1, 1.0], [0.1, 0.9, 1.0]];
        let x = Pt::new(Def::Ssi {
            first,
            second,
            drop: 2,
        });
        let [p, q] = first.map(v);
        let [r, s] = second.map(v);
        assert_eq!(orient2d(2, &p, &q, &x), Sign::Zero);
        assert_eq!(orient2d(2, &r, &s, &x), Sign::Zero);
        assert_eq!(cmp(2, &x, &v([0.0, 0.0, 1.0])), Ordering::Equal);
        let mid = Mid([&p, &q, &r]);
        assert_eq!(orient2d(2, &p, &q, &mid), orient2d(2, &p, &q, &r));
    }
}
