//! P2: where blends end, by intersection, not by case.
//!
//! On each face, the contact lines of the blends on its two edges at a
//! vertex meet at the face's setback corner there: three planes, the face's
//! and the two lines' (each square to the face). An edge not blended keeps
//! its own line, so a corner of no setback stays where it is.

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};

use super::classify::Planar;
use super::{Id, id};

/// Every face corner's setback corner (holes' corners too), keyed by face
/// index and vertex, or `None` where two edges of a face run on in one line
/// at a vertex and are set back by different distances there.
pub(crate) fn setback_corners(
    planar: &Planar,
    setback: &FastMap<(Id, Id), f64>,
    tol: Tolerances,
) -> OgeomResult<Option<FastMap<(usize, Id), Point>>> {
    let mut out = FastMap::default();
    for (fi, face) in planar.faces.iter().enumerate() {
        for ring in face.loops() {
            let n = ring.len();
            for k in 0..n {
                let (prev, v, next) = (ring[(k + n - 1) % n], ring[k], ring[(k + 1) % n]);
                // An edge's direction along the ring, the unit step into the
                // face square to it, and the face's setback from it.
                let leg = |a: Id, b: Id| -> OgeomResult<(Vector, Vector, f64)> {
                    let along = planar.points[&b] - planar.points[&a];
                    let length = along.magnitude();
                    if length <= tol.confusion() {
                        ogeom_bail!(Construction, "a face has an edge of no length");
                    }
                    let along = along / length;
                    let edge = planar.edges[&(a.min(b), a.max(b))];
                    let d = setback
                        .get(&(id(&face.shape), edge))
                        .copied()
                        .unwrap_or(0.0);
                    Ok((along, face.normal.cross(along), d))
                };
                let (u1, m1, d1) = leg(prev, v)?;
                let (u2, m2, d2) = leg(v, next)?;
                let at = planar.points[&v];
                let Some(offset) = corner_offset((u1, m1, d1), (u2, m2, d2), tol) else {
                    return Ok(None);
                };
                out.insert((fi, v), at + offset);
            }
        }
    }
    Ok(Some(out))
}

/// The step from a face's corner to its setback corner, in the face's
/// plane: `d1` in from the edge arriving along `u1` (stepping in along
/// `m1`) and `d2` in from the edge leaving along `u2` (along `m2`). An edge
/// not set back keeps the corner on its own line, so a corner of no setback
/// stays exactly where it is. `None` where the two edges run on in one line
/// with different setbacks.
fn corner_offset(
    (u1, m1, d1): (Vector, Vector, f64),
    (u2, m2, d2): (Vector, Vector, f64),
    tol: Tolerances,
) -> Option<Vector> {
    let turn = u1.cross(u2).magnitude();
    if turn <= tol.angular() {
        // Straight on: one line, one setback, or none to give.
        return ((d1 - d2).abs() <= tol.confusion()).then(|| m1 * d1);
    }
    if d1 == 0.0 && d2 == 0.0 {
        return Some(Vector::ZERO);
    }
    if d1 == 0.0 {
        // Along the arriving edge, back from the corner, to where the
        // leaving edge's setback line crosses it.
        return Some(u1 * (d2 / m2.dot(u1)));
    }
    if d2 == 0.0 {
        return Some(u2 * (d1 / m1.dot(u2)));
    }
    // In the plane: x = a u2 + b m2 with m1·x = d1 and m2·x = d2.
    let b = d2;
    let a = (d1 - b * m1.dot(m2)) / m1.dot(u2);
    Some(u2 * a + m2 * b)
}

/// A face's setback corner at a vertex, on a solid with curved faces: where
/// it stands, and where on the contact of each edge meeting there (two, or
/// more where the face's loop passes the vertex twice, as at a seam's foot).
#[derive(Debug, Clone)]
#[cfg_attr(
    not(test),
    allow(
        dead_code,
        reason = "blends of curved solids end through it from the M3 switch on"
    )
)]
pub(crate) struct Corner {
    pub point: Point,
    pub on: Vec<(Id, f64)>,
}

impl Corner {
    /// The corner's parameter on `edge`'s contact, `None` where `edge` does
    /// not meet there.
    #[cfg_attr(not(test), allow(dead_code, reason = "read by the trimmed assembly"))]
    pub fn on(&self, edge: Id) -> Option<f64> {
        self.on.iter().find(|(e, _)| *e == edge).map(|(_, t)| *t)
    }
}

/// The key a face's corner at vertex `v` is kept under where the face's
/// loop passes `v` again and puts it elsewhere (a face pinched at `v`, as
/// a top narrowing to nothing where a boss runs tangent into its side):
/// `v` with its top bit set, which no vertex's own key has.
#[cfg_attr(not(test), allow(dead_code, reason = "read by the trimmed assembly"))]
pub(crate) const fn pinched(v: Id) -> Id {
    v | 1 << 31
}

/// The corner of `face` at `v` an end of `edge` meets: the face's one
/// corner there; where the face passes `v` twice with a corner each, the
/// one on `edge`'s contact; and where both are (an edge closing on `v`),
/// the one ahead of `v` along the edge from its start where `start`, the
/// one behind it from its end otherwise. `None` where there is none.
///
/// # Errors
///
/// As evaluating the edge's curve.
#[cfg_attr(not(test), allow(dead_code, reason = "read by the trimmed assembly"))]
pub(crate) fn corner_for<'a>(
    corners: &'a FastMap<(usize, Id), Corner>,
    curved: &super::curved::Curved,
    face: usize,
    v: Id,
    edge: Id,
    start: bool,
    tol: Tolerances,
) -> OgeomResult<Option<&'a Corner>> {
    use ogeom_geom::Curve3d as _;
    let Some(second) = corners.get(&(face, pinched(v))) else {
        return Ok(corners.get(&(face, v)));
    };
    let on: Vec<&Corner> = [corners.get(&(face, v)), Some(second)]
        .into_iter()
        .flatten()
        .filter(|c| c.on(edge).is_some())
        .collect();
    match on.as_slice() {
        [one] => Ok(Some(*one)),
        [a, b] => {
            let e = &curved.edges[&edge];
            let t = if start { e.range.0 } else { e.range.1 };
            let ahead = e.curve.d1_at(t, tol)? * if start { 1.0 } else { -1.0 };
            let at = curved.points[&v];
            Ok(Some(
                if (a.point - at).dot(ahead) >= (b.point - at).dot(ahead) {
                    *a
                } else {
                    *b
                },
            ))
        }
        _ => Ok(None),
    }
}

/// Record the corner at `key`, standing at `point`, with parameters `on`; a
/// vertex the loop passes again keeps one corner with every edge's
/// parameter, or where it puts it elsewhere a second, under the vertex's
/// [`pinched`] key. `false` where a third pass puts it elsewhere again.
fn put(
    out: &mut FastMap<(usize, Id), Corner>,
    key: (usize, Id),
    point: Point,
    on: [(Id, f64); 2],
    tol: Tolerances,
) -> bool {
    match out.get_mut(&key) {
        Some(held) => {
            if held.point.distance(point) > tol.confusion() {
                return key.1 & pinched(0) == 0
                    && put(out, (key.0, pinched(key.1)), point, on, tol);
            }
            for (e, t) in on {
                if held.on(e).is_none() {
                    held.on.push((e, t));
                }
            }
            true
        }
        None => {
            let mut list = vec![on[0]];
            if on[1].0 != on[0].0 {
                list.push(on[1]);
            }
            out.insert(key, Corner { point, on: list });
            true
        }
    }
}

/// How far apart two edges' tangents may turn at a vertex and still run
/// on into each other: a millionth of a radian, looser than the angular
/// tolerance, since a vertex a boolean or a healing pass left a few
/// picometres off turns them by as much as that over its curve.
pub(crate) const TANGENT: f64 = 1e-6;

/// Every face corner's setback corner on a solid with curved faces, keyed
/// by face index and vertex: where the two contacts meeting at the vertex
/// cross on the face. A contact is the blend's, given in `contacts` by face
/// and edge as a bounded curve lying on the face; an edge not blended on
/// that face keeps its own curve, so a corner of no blend stays at its
/// vertex.
///
/// Where the two edges run on into each other at the vertex (a straight
/// round running tangent into a rim's), the contacts run on too, and the
/// corner is where each crosses the plane through the vertex square to the
/// edges; the two must agree there. A face's loop of one closed edge (a
/// rim) has its corner where its contact crosses that plane.
///
/// `None` where two contacts do not cross near their vertex (a blend
/// running off its face), or two that run on into each other do not agree
/// at the vertex.
///
/// # Errors
///
/// As evaluating the curves.
#[cfg_attr(
    not(test),
    allow(
        dead_code,
        reason = "blends of curved solids end through it from the M3 switch on"
    )
)]
pub(crate) fn curved_corners(
    curved: &super::curved::Curved,
    contacts: &FastMap<(usize, Id), ogeom_geom::Curve>,
    tol: Tolerances,
) -> OgeomResult<Option<FastMap<(usize, Id), Corner>>> {
    use ogeom_geom::Curve3d as _;
    let mut out = FastMap::default();
    for (fi, face) in curved.faces.iter().enumerate() {
        for walk in &face.loops {
            let n = walk.len();
            for k in 0..n {
                let (e_in, forward_in) = walk[k];
                let (e_out, forward_out) = walk[(k + 1) % n];
                let (arriving, leaving) = (&curved.edges[&e_in], &curved.edges[&e_out]);
                // The vertex between them: the end they share, read from
                // the vertices themselves (a face reversed in the solid
                // walks its stored loop the other way); the walk's sense
                // decides only where they share both ends.
                let ends_in = [arriving.ends.0, arriving.ends.1];
                let shared: Vec<Id> = ends_in
                    .iter()
                    .copied()
                    .filter(|v| *v == leaving.ends.0 || *v == leaving.ends.1)
                    .collect();
                let v = match shared.as_slice() {
                    [v] => *v,
                    [a, b] if a == b => *a,
                    [_, _] => {
                        if forward_in != face.reversed {
                            arriving.ends.1
                        } else {
                            arriving.ends.0
                        }
                    }
                    _ => return Ok(None),
                };
                // Which way each edge runs through `v`: the arriving one
                // ends there, the leaving one starts there.
                let closed = |e: &super::curved::CurvedEdge| e.ends.0 == e.ends.1;
                let forward_in = if closed(arriving) {
                    forward_in != face.reversed
                } else {
                    arriving.ends.1 == v
                };
                let forward_out = if closed(leaving) {
                    forward_out != face.reversed
                } else {
                    leaving.ends.0 == v
                };
                let at = curved.points[&v];
                let own = |edge: &super::curved::CurvedEdge| -> OgeomResult<ogeom_geom::Curve> {
                    let (lo, hi) = edge.curve.domain();
                    if (edge.range.0 - lo).abs() <= tol.parametric()
                        && (edge.range.1 - hi).abs() <= tol.parametric()
                    {
                        return Ok(edge.curve.clone());
                    }
                    Ok(ogeom_geom::TrimmedCurve::new(
                        edge.curve.clone(),
                        edge.range.0,
                        edge.range.1,
                        tol,
                    )?
                    .into())
                };
                let (c_in, blend_in) = match contacts.get(&(fi, e_in)) {
                    Some(c) => (c.clone(), true),
                    None => (own(arriving)?, false),
                };
                let (c_out, blend_out) = match contacts.get(&(fi, e_out)) {
                    Some(c) => (c.clone(), true),
                    None => (own(leaving)?, false),
                };
                if !blend_in && !blend_out {
                    let on = [
                        (e_in, end_param(arriving, forward_in)),
                        (e_out, start_param(leaving, forward_out)),
                    ];
                    if !put(&mut out, (fi, v), at, on, tol) {
                        return Ok(None);
                    }
                    continue;
                }
                // The edges' tangents at the vertex, in the walk's sense.
                let t_in = arriving.curve.d1_at(end_param(arriving, forward_in), tol)?
                    * if forward_in { 1.0 } else { -1.0 };
                let t_out = leaving
                    .curve
                    .d1_at(start_param(leaving, forward_out), tol)?
                    * if forward_out { 1.0 } else { -1.0 };
                let straight = n == 1
                    || t_in.cross(t_out).magnitude()
                        <= TANGENT * t_in.magnitude() * t_out.magnitude()
                        && t_in.dot(t_out) > 0.0;
                // Two fitted contacts running on into each other, each
                // through stations of its own, nearly touch at the vertex
                // (their stations there agree, or a crease's fitted seam
                // meets itself at a slight corner, where run on past their
                // ends they may part rather than cross): the corner is
                // midway between their points nearest the vertex, each
                // contact's edge widened past half the gap.
                if blend_in
                    && blend_out
                    && matches!(c_in, ogeom_geom::Curve::BSpline(_))
                    && matches!(c_out, ogeom_geom::Curve::BSpline(_))
                {
                    let nearest = |c: &ogeom_geom::Curve| -> OgeomResult<(f64, Point)> {
                        let s = super::trimmed::param_on(c, c.domain(), at, tol)?;
                        Ok((s, c.point_at(s, tol)?))
                    };
                    let ((s, p), (u, q)) = (nearest(&c_in)?, nearest(&c_out)?);
                    if p.distance(q) <= crate::marched::band_fit_target(tol) {
                        if !put(
                            &mut out,
                            (fi, v),
                            p.midpoint(q),
                            [(e_in, s), (e_out, u)],
                            tol,
                        ) {
                            return Ok(None);
                        }
                        continue;
                    }
                }
                if straight {
                    let t = t_in / t_in.magnitude();
                    // A contact running on into an edge left sharp along it:
                    // the bevel ends on the plane square to the edge there,
                    // its corner the contact's station on it, the face's
                    // own edge kept from the vertex on.
                    if blend_in != blend_out {
                        let (blend, e_blend, e_own, own_at) = if blend_in {
                            (&c_in, e_in, e_out, start_param(leaving, forward_out))
                        } else {
                            (&c_out, e_out, e_in, end_param(arriving, forward_in))
                        };
                        let Some((s, p)) = station(blend, at, t, tol)? else {
                            return Ok(None);
                        };
                        if !put(&mut out, (fi, v), p, [(e_blend, s), (e_own, own_at)], tol) {
                            return Ok(None);
                        }
                        continue;
                    }
                    let (Some((s, p)), Some((u, q))) =
                        (station(&c_in, at, t, tol)?, station(&c_out, at, t, tol)?)
                    else {
                        return Ok(None);
                    };
                    if p.distance(q) > tol.confusion() {
                        return Ok(None);
                    }
                    if !put(&mut out, (fi, v), p, [(e_in, s), (e_out, u)], tol) {
                        return Ok(None);
                    }
                    continue;
                }
                let crossings = |gap: f64| {
                    ogeom_intersect::curves::intersect_curves(
                        &c_in,
                        &c_out,
                        ogeom_intersect::curves::CurveCurveOptions {
                            gap,
                            ..ogeom_intersect::curves::CurveCurveOptions::default()
                        },
                        tol,
                    )
                };
                let mut found = crossings(tol.confusion())?;
                // A fitted contact stands within the bands' fit target of
                // its face, and may pass an edge on it that far off: its
                // edges are widened to how far it stands.
                if found.crossings.is_empty() {
                    found = crossings(crate::marched::band_fit_target(tol))?;
                }
                let Some(nearest) = found
                    .crossings
                    .iter()
                    .min_by(|a, b| a.point.distance(at).total_cmp(&b.point.distance(at)))
                else {
                    return Ok(None);
                };
                let on = [(e_in, nearest.on_a), (e_out, nearest.on_b)];
                // Where the two pass apart, the corner stands on the face's
                // own edge, which is exact: the fitted contact is what stands
                // off, within its widened edge.
                let point = if !blend_in {
                    c_in.point_at(nearest.on_a, tol)?
                } else if !blend_out {
                    c_out.point_at(nearest.on_b, tol)?
                } else {
                    nearest.point
                };
                if !put(&mut out, (fi, v), point, on, tol) {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(out))
}

/// The parameter an edge walked `forward` (or not) ends its walk at.
#[cfg_attr(not(test), allow(dead_code, reason = "used by curved_corners"))]
fn end_param(edge: &super::curved::CurvedEdge, forward: bool) -> f64 {
    if forward { edge.range.1 } else { edge.range.0 }
}

/// The parameter an edge walked `forward` (or not) starts its walk at.
#[cfg_attr(not(test), allow(dead_code, reason = "used by curved_corners"))]
fn start_param(edge: &super::curved::CurvedEdge, forward: bool) -> f64 {
    if forward { edge.range.0 } else { edge.range.1 }
}

/// Where `curve` crosses the plane through `at` square to the unit `t`: the
/// crossing nearest `at`, its parameter and point. `None` where it does not
/// cross within its domain.
#[cfg_attr(not(test), allow(dead_code, reason = "used by curved_corners"))]
fn station(
    curve: &ogeom_geom::Curve,
    at: Point,
    t: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, Point)>> {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = curve.domain();
    if !(lo.is_finite() && hi.is_finite()) {
        return Ok(None);
    }
    let side = |s: f64| -> OgeomResult<f64> { Ok((curve.point_at(s, tol)? - at).dot(t)) };
    // Every root over a fine sampling: a sample on the plane, or a sign
    // change narrowed by bisection; the one nearest `at` is kept.
    const SAMPLES: usize = 256;
    let mut roots: Vec<f64> = Vec::new();
    let mut prev = (lo, side(lo)?);
    if prev.1.abs() <= tol.confusion() {
        roots.push(lo);
    }
    for i in 1..=SAMPLES {
        #[expect(clippy::cast_precision_loss, reason = "a small sample count")]
        let s = lo + (hi - lo) * (i as f64) / (SAMPLES as f64);
        let here = (s, side(s)?);
        if here.1.abs() <= tol.confusion() {
            roots.push(s);
        } else if prev.1.abs() > tol.confusion() && prev.1.signum() != here.1.signum() {
            let (mut a, mut b) = (prev.0, here.0);
            let fa = prev.1;
            for _ in 0..200 {
                let m = f64::midpoint(a, b);
                if m <= a || m >= b {
                    break;
                }
                let fm = side(m)?;
                if fm.signum() == fa.signum() {
                    a = m;
                } else {
                    b = m;
                }
            }
            roots.push(f64::midpoint(a, b));
        }
        prev = here;
    }
    let mut best: Option<(f64, Point, f64)> = None;
    for s in roots {
        let p = curve.point_at(s, tol)?;
        let d = p.distance(at);
        if best.is_none_or(|b| d < b.2) {
            best = Some((s, p, d));
        }
    }
    Ok(best.map(|(s, p, _)| (s, p)))
}

#[cfg(test)]
mod curved_tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use crate::engine::curved::read;
    use ogeom_geom::{CircleCurve, Curve, LineCurve};
    use ogeom_math::{Circle, Direction, Frame};
    use ogeom_topo::{Model, Shape};

    const T: Tolerances = Tolerances::millimetres();

    fn up_at(x: f64, y: f64, z: f64) -> Frame {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    }

    fn segment(p: Point, q: Point) -> Curve {
        LineCurve::segment(p, q, T).unwrap().into()
    }

    /// The edge of `read` both of whose ends `keep`.
    fn edge_where(read: &super::super::curved::Curved, keep: impl Fn(Point, Point) -> bool) -> Id {
        let found: Vec<Id> = read
            .edges
            .iter()
            .filter(|(_, e)| keep(read.points[&e.ends.0], read.points[&e.ends.1]))
            .map(|(id, _)| *id)
            .collect();
        assert_eq!(found.len(), 1, "edges found: {:?}", found);
        found[0]
    }

    /// The faces of `read` walking `edge`.
    fn faces_of(read: &super::super::curved::Curved, edge: Id) -> Vec<usize> {
        (0..read.faces.len())
            .filter(|fi| {
                read.faces[*fi]
                    .loops
                    .iter()
                    .flatten()
                    .any(|(e, _)| *e == edge)
            })
            .collect()
    }

    /// On a box, two bottom edges meeting at a corner beside its sharp
    /// upright, bevelled unequally: the curved corners, from each bevel's
    /// contact lines given as segments, are the planar engine's
    /// (`corner_offset`'s 3x3 solves) at every vertex of every face.
    #[test]
    fn on_planes_the_curved_corners_are_the_planar_ones() {
        let mut m = Model::new();
        let s: Shape = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let planar = super::super::classify::read_planar(&m, &s, T)
            .unwrap()
            .unwrap();
        let curved = read(&m, &s, T).unwrap().unwrap();
        let wants = [
            (
                edge_where(&curved, |p, q| p.y.max(q.y) < 1e-9 && p.z.max(q.z) < 1e-9),
                1.5,
            ),
            (
                edge_where(&curved, |p, q| p.x.max(q.x) < 1e-9 && p.z.max(q.z) < 1e-9),
                0.7,
            ),
        ];
        let mut setback: FastMap<(Id, Id), f64> = FastMap::default();
        let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
        for (edge, d) in wants {
            let e = &curved.edges[&edge];
            let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
            let u = (b - a) / (b - a).magnitude();
            for fi in faces_of(&curved, edge) {
                let face = &curved.faces[fi];
                let shape = &face.shape;
                setback.insert((id(shape), edge), d);
                let (_, n) = ogeom_algo::face_normal(&m, shape, T).unwrap();
                let other = faces_of(&curved, edge)
                    .into_iter()
                    .find(|g| *g != fi)
                    .unwrap();
                let (_, n_other) =
                    ogeom_algo::face_normal(&m, &curved.faces[other].shape, T).unwrap();
                let mut inward = n.cross(u);
                if inward.dot(n_other) > 0.0 {
                    inward = -inward;
                }
                let shift = inward * d;
                contacts.insert(
                    (fi, edge),
                    segment(a + shift - u * 5.0, b + shift + u * 5.0),
                );
            }
        }
        let want = setback_corners(&planar, &setback, T).unwrap().unwrap();
        let got = curved_corners(&curved, &contacts, T).unwrap().unwrap();
        assert_eq!(got.len(), want.len());
        for ((pfi, v), p) in &want {
            let shape = &planar.faces[*pfi].shape;
            let cfi = curved
                .faces
                .iter()
                .position(|f| f.shape.is_same(shape))
                .unwrap();
            let q = got[&(cfi, *v)].point;
            assert!(p.distance(q) < 1e-12, "{p:?} against {q:?}");
        }
    }

    /// The corpus's D-shaped block (a drum of radius 3 cut in half along
    /// its axis, 10 long), its two straight edges chamfered 0.5: on the
    /// flat the contact is the line 0.5 in from the edge, on the curved
    /// face the ruling 0.5 along it. At each end the flat's corner is where
    /// its contact meets the cap's straight edge, the curved face's where
    /// its ruling meets the cap's arc, and the cap's the vertex itself.
    #[test]
    fn a_d_block_s_chamfer_contacts_end_on_its_caps() {
        let mut m = Model::new();
        let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
            .unwrap()
            .shape;
        let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
            .unwrap()
            .shape;
        let s = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
        let curved = read(&m, &s, T).unwrap().unwrap();
        let (r, d) = (3.0, 0.5);
        let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
        let mut expect: Vec<(usize, Id, Point, Point)> = Vec::new();
        for side in [1.0, -1.0] {
            let edge = edge_where(&curved, |p, q| {
                p.x.abs() < 1e-9
                    && q.x.abs() < 1e-9
                    && (p.y - 3.0 * side).abs() < 1e-9
                    && (p.z - q.z).abs() > 1.0
            });
            for fi in faces_of(&curved, edge) {
                let flat = matches!(
                    curved.faces[fi].surface,
                    ogeom_geom::SurfaceGeometry::Plane(_)
                );
                let (x, y) = if flat {
                    (0.0, side * (r - d))
                } else {
                    (r * (d / r).sin(), side * r * (d / r).cos())
                };
                contacts.insert(
                    (fi, edge),
                    segment(Point::new(x, y, -1.0), Point::new(x, y, 11.0)),
                );
                for z in [0.0, 10.0] {
                    expect.push((
                        fi,
                        edge,
                        Point::new(0.0, 3.0 * side, z),
                        Point::new(x, y, z),
                    ));
                }
            }
        }
        let got = curved_corners(&curved, &contacts, T).unwrap().unwrap();
        let blended: Vec<usize> = contacts.keys().map(|(fi, _)| *fi).collect();
        for (fi, edge, vertex, want) in expect {
            let (v, _) = curved
                .points
                .iter()
                .find(|(_, p)| p.distance(vertex) < 1e-9)
                .unwrap();
            let corner = &got[&(fi, *v)];
            let q = corner.point;
            assert!(q.distance(want) < 1e-9, "{q:?} against {want:?}");
            // The corner's parameters place it on the face's contact, the
            // arriving edge's or the leaving one's.
            let contact = &contacts[&(fi, edge)];
            {
                use ogeom_geom::Curve3d as _;
                let on: Vec<_> = corner
                    .on
                    .iter()
                    .map(|(_, t)| contact.point_at(*t, T))
                    .collect();
                assert!(
                    on.iter()
                        .any(|p| p.as_ref().is_ok_and(|p| p.distance(q) < 1e-9)),
                    "the corner is not where its parameters put it"
                );
            }
            // Each cap keeps its corner at the vertex: neither of its edges
            // there is blended.
            for cap in (0..curved.faces.len()).filter(|c| !blended.contains(c)) {
                if let Some(c) = got.get(&(cap, *v)) {
                    assert!(c.point.distance(vertex) < 1e-12, "a cap's corner moved");
                }
            }
        }
    }

    /// The corpus's stadium (a 10 by 4 box between two drums of radius 2,
    /// 3 tall), its top loop given a contact 0.5 in on the top face: along
    /// each straight edge the line 0.5 in, round each half rim the circle
    /// of radius 1.5. Where a straight edge runs tangent into a rim the two
    /// contacts run on into each other, and the corner is their common
    /// point on the plane through the vertex square to the loop: 0.5 in
    /// from the vertex. A rim split at its own seam is such a junction too.
    #[test]
    fn a_stadium_s_top_contacts_run_on_into_each_other() {
        let mut m = Model::new();
        let core = ogeom_algo::make_box(&mut m, up_at(-5.0, -2.0, 0.0), (10.0, 4.0, 3.0), T)
            .unwrap()
            .shape;
        let left = ogeom_algo::make_cylinder(&mut m, up_at(-5.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let right = ogeom_algo::make_cylinder(&mut m, up_at(5.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let st = ogeom_bool::fuse(&mut m, &core, &left, T).unwrap().shape;
        let s = ogeom_bool::fuse(&mut m, &st, &right, T).unwrap().shape;
        // The fuse leaves the top split along seams; read as the engine
        // reads it, with its seams dissolved.
        let s = ogeom_heal::unify_same_domain(&mut m, &s, T)
            .unwrap()
            .0
            .shape;
        let s = ogeom_heal::merge_edges(&mut m, &s, T).unwrap().0.shape;
        let curved = read(&m, &s, T).unwrap().unwrap();
        let top = curved
            .faces
            .iter()
            .position(|f| {
                matches!(&f.surface, ogeom_geom::SurfaceGeometry::Plane(_))
                    && f.loops.iter().flatten().all(|(e, _)| {
                        let e = &curved.edges[e];
                        (curved.points[&e.ends.0].z - 3.0).abs() < 1e-9
                    })
            })
            .unwrap();
        let d = 0.5;
        let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
        for (edge, _) in curved.faces[top].loops.iter().flatten() {
            let e = &curved.edges[edge];
            let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
            let contact = match &e.curve {
                Curve::Line(_) => {
                    let side = a.y.signum();
                    let shift = Vector::new(0.0, -side * d, 0.0);
                    segment(
                        a + shift + Vector::new(-6.0, 0.0, 0.0) * a.x.signum().abs(),
                        b + shift,
                    )
                }
                Curve::Circle(c) => {
                    let centre = c.circle().centre();
                    CircleCurve::new(
                        Circle::new(up_at(centre.x, centre.y, 3.0), 2.0 - d, T).unwrap(),
                    )
                    .into()
                }
                other => panic!("an edge of the top on {other:?}"),
            };
            let _ = (a, b);
            contacts.insert((top, *edge), contact);
        }
        // Lines run past their ends on both sides, so they reach the
        // stations at the junctions.
        for (edge, _) in curved.faces[top].loops.iter().flatten() {
            let e = &curved.edges[edge];
            if let Curve::Line(_) = e.curve {
                let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
                let side = a.y.signum();
                let shift = Vector::new(0.0, -side * d, 0.0);
                let u = (b - a) / (b - a).magnitude();
                contacts.insert((top, *edge), segment(a + shift - u, b + shift + u));
            }
        }
        let got = curved_corners(&curved, &contacts, T).unwrap().unwrap();
        let mut junctions = 0;
        for (v, p) in &curved.points {
            if (p.z - 3.0).abs() > 1e-9 {
                continue;
            }
            let Some(c) = got.get(&(top, *v)) else {
                continue;
            };
            // At a straight edge's end the corner is 0.5 in across the
            // loop; where a rim is split at its own seam (y = 0), 0.5 in
            // towards the rim's axis.
            let want = if p.y.abs() > 1e-9 {
                Point::new(p.x, p.y - p.y.signum() * d, 3.0)
            } else {
                Point::new(p.x - p.x.signum() * d, 0.0, 3.0)
            };
            assert!(
                c.point.distance(want) < 1e-9,
                "{:?} against {want:?}",
                c.point
            );
            junctions += 1;
        }
        assert!(junctions >= 4, "the top's junctions: {junctions}");
    }
}
