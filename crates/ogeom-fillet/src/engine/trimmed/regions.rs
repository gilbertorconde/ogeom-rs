//! The sheet's points and rings, and the pieces of a curve between them.

use super::*;

/// The sheet's points, each held once within the confusion distance.
pub(super) struct Points {
    pub(super) points: Vec<Point>,
    /// How close two points are one: the confusion distance, or a fitted
    /// contact's width, which one corner reached by two routes (on an exact
    /// edge and at the contact's own end) may stand apart by.
    pub(super) merge: f64,
}

impl Points {
    pub(super) fn index(&mut self, p: Point) -> usize {
        if let Some(i) = self.points.iter().position(|q| q.distance(p) <= self.merge) {
            return i;
        }
        self.points.push(p);
        self.points.len() - 1
    }
}

/// The width an edge of a rolled run's fitted piece is built to: its chart's
/// curve is fitted within it.
pub(super) const ROLLED_EDGE_WIDTH: f64 = 2e-5;

/// A ring of sheet edges, each with whether it is walked forward.
type Ring = Vec<(usize, bool)>;

/// A part of a piece a rolled round cuts: from and to (each a parameter
/// and its point), and what the round does to it, if anything.
pub(super) type Part<'a> = ((f64, Point), (f64, Point), Option<&'a Fate>);

/// What a rolled round does to a piece of a contact or an edge on a face.
#[derive(Debug, Clone)]
pub(super) enum Fate {
    /// Kept, cut where the piece ends (an edge the ball rests on, on the
    /// faces it does not leave).
    Keep,
    /// Gone: under the round.
    Drop,
    /// Taken by the curve given, run from the span's low end to its high
    /// end over the curve's domain (the contact along a run).
    Replace(Curve),
}

/// A span of a curve's parameters a rolled round cuts out, `lo` below `hi`
/// (past the domain's end for a periodic curve), with the points there.
#[derive(Debug, Clone)]
pub(super) struct Span {
    pub(super) lo: (f64, Point),
    pub(super) hi: (f64, Point),
    pub(super) fate: Fate,
}

/// The material the edges a rolled run rests on bound, as a region the
/// certify reads: the face its contact left and each other face of the
/// first edge, all of them at a convex edge (where the face left turns
/// down past the other's plane), any at a concave one; `None` off planes
/// and walls of revolution.
///
/// # Errors
///
/// As evaluating the edge's curve.
pub(super) fn resting_region(
    curved: &Curved,
    r: &RolledRun,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_algo::HalfSpaces>> {
    let Some(&e) = r.run.over.first() else {
        return Ok(None);
    };
    let left = r.enter.face;
    let edge = &curved.edges[&e];
    let q = r.enter.point;
    let tangent = edge
        .curve
        .d1_at(param_on(&edge.curve, edge.range, q, tol)?, tol)?;
    let Some(n_left) = outward(curved, left, q) else {
        return Ok(None);
    };
    let Some(m_left) = material(curved, left) else {
        return Ok(None);
    };
    let mut convex = true;
    let mut pieces = vec![m_left];
    for g in faces_of(curved, e) {
        if g == left {
            continue;
        }
        let (Some(n_g), Some(m_g)) = (outward(curved, g, q), material(curved, g)) else {
            return Ok(None);
        };
        // Into the other face from the edge: on the left of the edge as
        // that face walks it, about its outward normal.
        let Some(&(_, w)) = curved.faces[g]
            .loops
            .iter()
            .flatten()
            .find(|(x, _)| *x == e)
        else {
            return Ok(None);
        };
        let walk = if w == curved.faces[g].reversed {
            -tangent
        } else {
            tangent
        };
        let into = n_g.cross(walk);
        if n_left.dot(into) > 0.0 {
            convex = false;
        }
        pieces.push(m_g);
    }
    Ok(Some(if convex {
        ogeom_algo::HalfSpaces::Every(pieces)
    } else {
        ogeom_algo::HalfSpaces::Any(pieces)
    }))
}

/// The ring from `first` (a sheet edge and whether walked forward) on
/// through each of `rest` once, each walked from where the last ended;
/// `None` where they do not close one ring.
pub(super) fn chained_ring(
    edges: &[SheetCurve],
    first: (usize, bool),
    rest: &[usize],
) -> Option<Vec<(usize, bool)>> {
    let end = |(k, f): (usize, bool)| if f { edges[k].ends.1 } else { edges[k].ends.0 };
    let start = if first.1 {
        edges[first.0].ends.0
    } else {
        edges[first.0].ends.1
    };
    let mut ring = vec![first];
    let mut left: Vec<usize> = rest.to_vec();
    let mut at = end(first);
    while !left.is_empty() {
        let n = left
            .iter()
            .position(|k| edges[*k].ends.0 == at || edges[*k].ends.1 == at)?;
        let k = left.remove(n);
        let forward = edges[k].ends.0 == at;
        ring.push((k, forward));
        at = end((k, forward));
    }
    (at == start).then_some(ring)
}

/// The pieces of `rings` (each a sheet edge and whether walked forward)
/// rejoined end to end into closed rings; `None` where a point starts two
/// pieces or none.
pub(super) fn rejoined(
    edges: &[SheetCurve],
    rings: &[Vec<(usize, bool)>],
) -> Option<Vec<Vec<(usize, bool)>>> {
    let ends = |(k, f): (usize, bool)| {
        let (a, z) = edges[k].ends;
        if f { (a, z) } else { (z, a) }
    };
    let pieces: Vec<(usize, bool)> = rings.iter().flatten().copied().collect();
    let mut from: FastMap<usize, usize> = FastMap::default();
    for (n, p) in pieces.iter().enumerate() {
        if from.insert(ends(*p).0, n).is_some() {
            return None;
        }
    }
    let mut used = vec![false; pieces.len()];
    let mut out = Vec::new();
    for n in 0..pieces.len() {
        if used[n] {
            continue;
        }
        let start = ends(pieces[n]).0;
        let mut ring = Vec::new();
        let mut at = n;
        loop {
            if used[at] {
                return None;
            }
            used[at] = true;
            ring.push(pieces[at]);
            let next = ends(pieces[at]).1;
            if next == start {
                break;
            }
            at = *from.get(&next)?;
        }
        out.push(ring);
    }
    Some(out)
}

/// `rings` on a plane parted into faces: each ring running counter-clockwise
/// about the plane's own normal with the rings inside it running the other
/// way, the largest first. On any other surface the rings as one face,
/// where none splits off (as many outer rings as one face can hold), else
/// `None`.
///
/// # Errors
///
/// As evaluating the edges' curves.
pub(super) fn grouped(
    surface: &SurfaceGeometry,
    edges: &[SheetCurve],
    rings: Vec<Ring>,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Vec<Ring>>>> {
    let SurfaceGeometry::Plane(plane) = surface else {
        return Ok(None);
    };
    let frame = plane.plane().frame();
    let (o, x, y) = (frame.origin(), frame.x().vector(), frame.y().vector());
    // Each ring's outline in the plane, eight points an edge.
    let mut outlines: Vec<Vec<(f64, f64)>> = Vec::with_capacity(rings.len());
    for ring in &rings {
        let mut outline = Vec::new();
        for (k, f) in ring {
            let e = &edges[*k];
            for j in 0..8 {
                let s = f64::from(j) / 8.0;
                let s = if *f { s } else { 1.0 - s };
                let t = e.range.0 + (e.range.1 - e.range.0) * s;
                let p = e.curve.point_at(t, tol)?;
                outline.push(((p - o).dot(x), (p - o).dot(y)));
            }
        }
        outlines.push(outline);
    }
    let area = |o: &[(f64, f64)]| {
        let mut a = 0.0;
        for k in 0..o.len() {
            let (p, q) = (o[k], o[(k + 1) % o.len()]);
            a += p.0 * q.1 - q.0 * p.1;
        }
        0.5 * a
    };
    let inside = |p: (f64, f64), o: &[(f64, f64)]| {
        let mut odd = false;
        for k in 0..o.len() {
            let (a, b) = (o[k], o[(k + 1) % o.len()]);
            if (a.1 > p.1) != (b.1 > p.1) && p.0 < a.0 + (p.1 - a.1) * (b.0 - a.0) / (b.1 - a.1) {
                odd = !odd;
            }
        }
        odd
    };
    let areas: Vec<f64> = outlines.iter().map(|o| area(o)).collect();
    let mut outer: Vec<usize> = (0..rings.len()).filter(|k| areas[*k] > 0.0).collect();
    outer.sort_by(|a, b| areas[*b].total_cmp(&areas[*a]));
    if outer.is_empty() {
        return Ok(None);
    }
    let mut groups: Vec<Vec<usize>> = outer.iter().map(|k| vec![*k]).collect();
    for k in (0..rings.len()).filter(|k| areas[*k] <= 0.0) {
        // The smallest outer ring holding it.
        let holder = outer
            .iter()
            .enumerate()
            .rev()
            .find(|(_, g)| inside(outlines[k][0], &outlines[**g]))
            .map(|(n, _)| n);
        let Some(n) = holder else {
            return Ok(None);
        };
        groups[n].push(k);
    }
    let mut rings: Vec<Option<Ring>> = rings.into_iter().map(Some).collect();
    Ok(Some(
        groups
            .into_iter()
            .map(|g| g.into_iter().filter_map(|k| rings[k].take()).collect())
            .collect(),
    ))
}

/// The period of `curve`, where it closes.
pub(super) fn period_of(curve: &Curve, tol: Tolerances) -> Option<f64> {
    (curve.is_periodic() || curve.is_closed(tol)).then(|| {
        let (lo, hi) = curve.domain();
        hi - lo
    })
}

/// `t` moved by whole periods to stand at or past `from`, below `from` plus
/// the period.
pub(super) fn wrapped(t: f64, from: f64, period: Option<f64>) -> f64 {
    let Some(p) = period else {
        return t;
    };
    let mut t = t;
    while t < from {
        t += p;
    }
    while t >= from + p {
        t -= p;
    }
    t
}

/// The span from `a` to `b` on `curve` (each a parameter and its point),
/// on a periodic curve the way round `pick` takes (given its parameters'
/// low and high ends).
pub(super) fn span_between(
    curve: &Curve,
    a: (f64, Point),
    b: (f64, Point),
    pick: impl Fn(f64, f64) -> OgeomResult<bool>,
    fate: Fate,
    tol: Tolerances,
) -> OgeomResult<Span> {
    let period = period_of(curve, tol);
    let (lo, hi) = if a.0 <= b.0 { (a, b) } else { (b, a) };
    let Some(p) = period else {
        return Ok(Span { lo, hi, fate });
    };
    let hi_up = (wrapped(hi.0, lo.0, Some(p)), hi.1);
    if pick(lo.0, hi_up.0)? {
        Ok(Span {
            lo,
            hi: hi_up,
            fate,
        })
    } else {
        let lo_up = (wrapped(lo.0, hi.0, Some(p)), lo.1);
        Ok(Span {
            lo: hi,
            hi: lo_up,
            fate,
        })
    }
}

/// A piece of `curve` from parameter `from` to `to`, running the way the
/// curve's parameter does when `rising`: its stored range (always
/// increasing) and whether the walk runs against it. A closed curve's piece
/// from a parameter to itself is its whole turn.
pub(super) fn piece(
    curve: &Curve,
    from: f64,
    to: f64,
    rising: bool,
    tol: Tolerances,
) -> ((f64, f64), bool) {
    let period = if curve.is_periodic() || curve.is_closed(tol) {
        let (lo, hi) = curve.domain();
        Some(hi - lo)
    } else {
        None
    };
    let same = (from - to).abs() <= tol.parametric();
    // A whole turn of a closed curve that does not repeat past its domain
    // (a fitted closed contact starting at its corner) is that domain.
    if same && period.is_some() && !curve.is_periodic() {
        return (curve.domain(), !rising);
    }
    let ((lo, hi), against) = match (rising, period) {
        (true, Some(p)) => {
            let mut to = to;
            while to <= from + if same { 0.5 * p } else { 0.0 } {
                to += p;
            }
            ((from, to), false)
        }
        (false, Some(p)) => {
            let mut from = from;
            while from <= to + if same { 0.5 * p } else { 0.0 } {
                from += p;
            }
            ((to, from), true)
        }
        (true, None) => ((from, to), false),
        (false, None) => ((to, from), true),
    };
    // A closed curve that does not repeat: an end within a ten-thousandth
    // of a turn past its domain (a corner at its start read a turn on) is
    // its domain's end.
    match period {
        Some(p) if !curve.is_periodic() => {
            let (a, b) = curve.domain();
            let near = |x: f64, y: f64| (x - y).abs() <= 1e-4 * p;
            let lo = if near(lo, a + p) || near(lo, a) {
                a
            } else {
                lo
            };
            let hi = if near(hi, b + p) || near(hi, b) {
                b
            } else {
                hi
            };
            ((lo, hi), against)
        }
        _ => ((lo, hi), against),
    }
}

/// Whether `curve` runs at parameter `t` the way `along` points.
pub(super) fn runs(curve: &Curve, t: f64, along: Vector, tol: Tolerances) -> OgeomResult<bool> {
    Ok(curve.d1_at(t, tol)?.dot(along) > 0.0)
}

/// The sheet's edges, each made once and keyed by what it is.
#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub(super) enum Key {
    /// A blended edge's contact on a face.
    Contact(usize, Id),
    /// An edge left sharp, between the corners at its ends.
    Sharp(Id),
    /// The piece of an edge left sharp at a vertex between its two faces'
    /// corners, where they stand apart on it: the face that runs on
    /// further shares it with the corner closing the vertex.
    Tail(Id, Id),
    /// A bevel's end, a segment between two of the sheet's points.
    Section(usize, usize),
    /// A piece of a contact or an edge a rolled round cuts: on face `.0`
    /// the contact of edge `.1` (or, `.0` [`usize::MAX`], edge `.1` itself,
    /// shared by its faces), between the sheet's points `.2` and `.3`, the
    /// lesser first.
    Piece(usize, Id, usize, usize),
    /// A contact on a consumed face, between two of the sheet's points,
    /// shared by the strips meeting along it; the last its middle's index
    /// among the middles, which tells apart the two ways round a loop.
    Shared(usize, usize, usize, usize),
}
