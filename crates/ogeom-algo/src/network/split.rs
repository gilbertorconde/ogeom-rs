//! A three-sided patch beside four-sided ones, split into three.

use ogeom_math::Point;

use super::topology::Topology;
use super::{ChainedPatch, CurveDef, PatchSide};

/// The network with each three-sided patch beside four-sided ones split
/// into three four-sided patches, and the creases on its curves.
pub(super) struct Split {
    pub patches: Vec<ChainedPatch>,
    pub curves: Vec<CurveDef>,
    pub creases: Vec<usize>,
}

/// The points of a drawn curve at its chord parameters over the unit
/// interval.
fn chords(drawn: &[Point]) -> Vec<f64> {
    let mut at = vec![0.0];
    for w in drawn.windows(2) {
        at.push(at[at.len() - 1] + w[0].distance(w[1]));
    }
    let total = at[at.len() - 1];
    at.iter()
        .map(|t| if total > 0.0 { t / total } else { 0.0 })
        .collect()
}

/// The point of a drawn curve at chord parameter `t`, and the index of the
/// first of its points past it.
fn along(drawn: &[Point], at: &[f64], t: f64) -> (Point, usize) {
    let k = at.partition_point(|x| *x <= t).clamp(1, drawn.len() - 1);
    let (a, b) = (at[k - 1], at[k]);
    let f = if b > a {
        ((t - a) / (b - a)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (drawn[k - 1] + (drawn[k] - drawn[k - 1]) * f, k)
}

/// The patches eligible: three-sided, each other side one curve, the
/// point inside the network, and some patch there four-sided at it.
fn eligible(topo: &Topology, patches: &[ChainedPatch]) -> Vec<(usize, usize)> {
    let missing_at = |q: usize, v: usize| {
        patches[q]
            .sides
            .iter()
            .enumerate()
            .any(|(k, side)| side.is_empty() && topo.corner_vertex[q][k] == v)
    };
    patches
        .iter()
        .enumerate()
        .filter_map(|(p, patch)| {
            let empty: Vec<usize> = (0..4).filter(|&k| patch.sides[k].is_empty()).collect();
            let [m] = empty.as_slice() else {
                return None;
            };
            if (0..4).any(|k| k != *m && patch.sides[k].len() != 1) {
                return None;
            }
            let v = topo.corner_vertex[p][*m];
            let inside = topo
                .curve_ends
                .iter()
                .enumerate()
                .filter(|&(_, &(a, b))| a == v || b == v)
                .all(|(c, _)| topo.users[c].len() == 2);
            let beside = topo.touching[v]
                .iter()
                .any(|&(q, _)| q != p && !missing_at(q, v));
            (inside && beside).then_some((p, *m))
        })
        .collect()
}

/// The points along the straight line from `a` to `b` in a patch's
/// parameters, `count` of them, mapped by `surface`.
fn drawn_line(
    surface: &dyn Fn(f64, f64) -> Point,
    a: (f64, f64),
    b: (f64, f64),
    count: usize,
) -> Vec<Point> {
    let last = super::precise(count - 1);
    (0..count)
        .map(|i| {
            let f = super::precise(i) / last;
            surface((b.0 - a.0).mul_add(f, a.0), (b.1 - a.1).mul_add(f, a.1))
        })
        .collect()
}

/// How many points each curve drawn inside a split patch has.
const INNER_POINTS: usize = 33;

/// Split each three-sided patch beside four-sided ones, `None` where none
/// is.
///
/// The patch's three sides are each split at their middles by chord, the
/// halves two curves, and a curve drawn from each middle to the patch's
/// centre: in its parameters, a third of the way from the middle of the
/// side opposite its missing one towards that side, where the triangle's
/// centroid lies when the parameters are even. The curves inside are
/// drawn on the bilinearly blended Coons patch of its sides, the surface
/// its seed is. The patch is then three: one at its point, between the
/// halves next to it, and one at each of its other two corners. Each
/// patch beside a split curve takes the two halves in its place, a chain
/// meeting at the middle.
///
/// A split patch's own place takes the patch at its point and a split
/// curve's its first half; the other two patches of each split patch, in
/// order, follow the caller's patches, and the second halves, then the
/// three curves inside each split patch, follow the caller's curves. A
/// split crease is creased in both halves.
pub(super) fn three_sided(
    topo: &Topology,
    patches: &[ChainedPatch],
    curves: &[CurveDef],
    drawn: &[Vec<Point>],
    creases: &[usize],
) -> Option<Split> {
    let split = eligible(topo, patches);
    if split.is_empty() {
        return None;
    }
    let at: Vec<Vec<f64>> = drawn.iter().map(|d| chords(d)).collect();
    let mut out_curves = curves.to_vec();
    let mut out_creases = creases.to_vec();
    // Each split curve's second half.
    let mut second = vec![None; curves.len()];
    for &(p, m) in &split {
        for k in (0..4).filter(|&k| k != m) {
            let c = patches[p].sides[k][0].curve;
            if second[c].is_some() {
                continue;
            }
            let (middle, past) = along(&drawn[c], &at[c], 0.5);
            let d = &drawn[c];
            let mut first: Vec<Point> = d[..past].to_vec();
            if first.last() != Some(&middle) {
                first.push(middle);
            }
            let mut rest = vec![middle];
            rest.extend(d[past..].iter().copied().filter(|q| *q != middle));
            out_curves[c] = CurveDef::Polyline(first);
            second[c] = Some(out_curves.len());
            if creases.contains(&c) {
                out_creases.push(out_curves.len());
            }
            out_curves.push(CurveDef::Polyline(rest));
        }
    }
    // A side's halves in its walk's order.
    let halves = |side: PatchSide| -> [PatchSide; 2] {
        let (a, b) = (side.curve, second[side.curve].unwrap_or(side.curve));
        if side.reversed {
            [
                PatchSide {
                    curve: b,
                    reversed: true,
                },
                PatchSide {
                    curve: a,
                    reversed: true,
                },
            ]
        } else {
            [
                PatchSide {
                    curve: a,
                    reversed: false,
                },
                PatchSide {
                    curve: b,
                    reversed: false,
                },
            ]
        }
    };
    let mut out_patches: Vec<ChainedPatch> = patches
        .iter()
        .map(|patch| ChainedPatch {
            sides: patch.sides.clone().map(|side| {
                side.into_iter()
                    .flat_map(|piece| {
                        if second[piece.curve].is_some() {
                            halves(piece).to_vec()
                        } else {
                            vec![piece]
                        }
                    })
                    .collect()
            }),
        })
        .collect();
    for &(p, m) in &split {
        let patch = &patches[p];
        let point = |k: usize, f: f64| -> Point {
            let side = patch.sides[k][0];
            let c = side.curve;
            along(&drawn[c], &at[c], if side.reversed { 1.0 - f } else { f }).0
        };
        let v = {
            let (k, f) = if m == 0 { (3, 1.0) } else { (m - 1, 1.0) };
            point(k, f)
        };
        // Side `k` at the share `f` of its walk; the missing side at the point.
        let walk = |k: usize, f: f64| if k == m { v } else { point(k, f) };
        let coons = |x: f64, y: f64| -> Point {
            let (s0, s1, s2, s3) = (walk(0, x), walk(1, y), walk(2, 1.0 - x), walk(3, 1.0 - y));
            let (c00, c10, c11, c01) = (walk(0, 0.0), walk(0, 1.0), walk(1, 1.0), walk(2, 1.0));
            let ruled = s0.to_vector() * (1.0 - y)
                + s2.to_vector() * y
                + s3.to_vector() * (1.0 - x)
                + s1.to_vector() * x;
            let bilinear = c00.to_vector() * ((1.0 - x) * (1.0 - y))
                + c10.to_vector() * (x * (1.0 - y))
                + c11.to_vector() * (x * y)
                + c01.to_vector() * ((1.0 - x) * y);
            Point::ORIGIN + (ruled - bilinear)
        };
        // Each side's middle in the patch's parameters.
        let middle = |k: usize| -> (f64, f64) {
            match k {
            0 => (0.5, 0.0),
            1 => (1.0, 0.5),
            2 => (0.5, 1.0),
            _ => (0.0, 0.5),
            }
        };
        let (o, d) = (middle((m + 2) % 4), middle(m));
        let centre = (
            (d.0 - o.0).mul_add(1.0 / 3.0, o.0),
            (d.1 - o.1).mul_add(1.0 / 3.0, o.1),
        );
        let [k1, k2, k3] = [(m + 1) % 4, (m + 2) % 4, (m + 3) % 4];
        let inner: Vec<usize> = [k1, k2, k3]
            .iter()
            .map(|&k| {
                let mut line = drawn_line(&coons, centre, middle(k), INNER_POINTS);
                let last = line.len() - 1;
                line[0] = coons(centre.0, centre.1);
                line[last] = point(k, 0.5);
                out_curves.push(CurveDef::Polyline(line));
                out_curves.len() - 1
            })
            .collect();
        let [ia, ib, ic] = [inner[0], inner[1], inner[2]];
        let (h1, h2, h3) = (
            halves(patch.sides[k1][0]),
            halves(patch.sides[k2][0]),
            halves(patch.sides[k3][0]),
        );
        let to = |curve: usize, reversed: bool| PatchSide { curve, reversed };
        let quad = |sides: [PatchSide; 4]| ChainedPatch {
            sides: sides.map(|s| vec![s]),
        };
        out_patches[p] = quad([h1[0], to(ia, true), to(ic, false), h3[1]]);
        out_patches.push(quad([h1[1], h2[0], to(ib, true), to(ia, false)]));
        out_patches.push(quad([h2[1], h3[0], to(ic, true), to(ib, false)]));
    }
    Some(Split {
        patches: out_patches,
        curves: out_curves,
        creases: out_creases,
    })
}
