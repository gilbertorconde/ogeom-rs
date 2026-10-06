//! The middle path of a pipe-like solid: the curve its cross-sections'
//! centroids trace from one end face to the other.
//!
//! The solid is meshed once, and the march cuts that mesh with planes.
//! Each station is found by predictor and corrector: a plane square to the
//! current tangent one step ahead gives a first centroid, the chord to it
//! turns the tangent (on a circular spine the chord's direction is the mean
//! of its end tangents), and the plane square to the turned tangent through
//! the centroid gives the next, until the station stops moving. A plane
//! square to a tube's own spine cuts it in a section centred on the spine,
//! so the stations sit on the spine to within the mesh's chord.
//!
//! The spine is fitted through the stations and then measured: planes
//! square to the fitted curve between the stations are cut again, and the
//! largest distance from a cut's centroid to the curve is the deviation
//! reported. A deviation over the tolerance halves the step and marches
//! again.
//!
//! A sharp (mitred) corner stops the march: past the point where a section
//! square to the leg reaches the mitre, it takes in the next leg too and
//! its centroid jumps aside. The march then runs back from the end face as
//! well, each leg is kept up to the last station whose section cannot
//! reach the mitre, and the two legs, carried on straight along their last
//! tangents, must meet: there is the corner. The corner is measured too,
//! by the section on the plane halving it, whose centroid a mitred tube
//! has on the corner.

use ogeom_algo::{Built, History, make_edge_between, make_vertex, make_wire};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, LineCurve};
use ogeom_math::{Point, Vector};
use ogeom_mesh::{Deflection, triangulate, triangulate_face};
use ogeom_topo::{Filter, Model, Shape, ShapeType, Triangulation, explore};
use std::collections::HashMap;

/// A middle path and how closely it follows the solid's sections.
#[derive(Debug, Clone)]
pub struct MiddlePath {
    /// The wire of the path, from the start face's centroid to the end
    /// face's, and its history: both end faces generate it.
    pub built: Built,
    /// The largest distance, measured, from the centroid of a section cut
    /// square to the path to the path itself.
    pub deviation: f64,
}

/// The steps the march may halve before it gives up on the tolerance.
const MAX_REFINEMENTS: usize = 5;

/// The stations one march may place before it is taken to be lost.
const MAX_STATIONS: usize = 4000;

/// The most a station may turn the march's tangent, as the cosine of the
/// angle: a turn sharper than this is a section that has run into the
/// next leg past a sharp corner. A smooth spine turns by the step over its
/// radius of curvature, which a halved step brings under it.
const SHARPEST_TURN: f64 = core::f64::consts::FRAC_1_SQRT_2;

/// The centre line of a pipe-like `solid`, from its face `start` to its
/// face `end`.
///
/// Each point of the path is the centroid of the solid's section by a
/// plane square to the path there. The path starts at the centroid of
/// `start` and ends at the centroid of `end`, and it is a straight edge
/// where every station lies on one line and a fitted spline otherwise.
/// A tube with one sharp (mitred) corner has a path of two edges, one per
/// leg, meeting where the legs carried on straight meet; the sections
/// near the corner, which reach into the other leg, are not measured, and
/// the section on the plane halving the corner is. `tolerance` bounds the
/// reported deviation; the solid is meshed at a quarter of it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `solid` is not a solid, either face is not one of its faces, the two
/// faces are the same, the tolerance is not a positive distance, or the
/// mesh does not close. [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone)
/// if a section square to the march finds no material (the solid is not
/// pipe-like between the faces), the march never reaches `end`, it turns
/// back on itself where its legs do not meet at one corner, or the
/// deviation stays over the tolerance after every refinement.
pub fn middle_path(
    model: &mut Model,
    solid: &Shape,
    start: &Shape,
    end: &Shape,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<MiddlePath> {
    if !tolerance.is_finite() || tolerance <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a middle path to {tolerance} is not a distance"
        );
    }
    if model.kind_of(solid)? != ShapeType::Solid {
        ogeom_bail!(Construction, "a middle path runs through a solid");
    }
    let faces = explore(model, solid, Filter::OfType(ShapeType::Face))?;
    for face in [start, end] {
        if !faces.iter().any(|f| f.is_same(face)) {
            ogeom_bail!(Construction, "the end faces must be faces of the solid");
        }
    }
    if start.is_same(end) {
        ogeom_bail!(Construction, "a middle path needs two different end faces");
    }

    let deflection = Deflection {
        chord: tolerance * 0.25,
        ..Deflection::default()
    };
    let mesh = Slicer::new(triangulate(model, solid, deflection, tol)?);
    if !mesh.closed {
        ogeom_bail!(
            Construction,
            "the solid's mesh does not close, so its sections are not regions"
        );
    }
    let (c0, n0, area0) = face_centroid(model, start, deflection, tol)?;
    let (ce, ne, _) = face_centroid(model, end, deflection, tol)?;
    let size = (area0 / core::f64::consts::PI).sqrt();
    if size <= tol.confusion() {
        ogeom_bail!(Construction, "the start face has no area to start from");
    }

    // Into the solid: the side where a nearby plane cuts material round
    // an end face's centroid.
    let probe = size * 0.05;
    let inward = |c: Point, n: Vector| {
        if mesh.section(c + n * probe, n, tol).is_some() {
            Some(n)
        } else if mesh.section(c - n * probe, n, tol).is_some() {
            Some(-n)
        } else {
            None
        }
    };
    let Some(t0) = inward(c0, n0) else {
        ogeom_bail!(NotDone, "no material lies behind the start face");
    };

    let mut step = size;
    let mut last = Err("the middle path was not marched");
    for _ in 0..=MAX_REFINEMENTS {
        let forward = march(&mesh, c0, t0, ce, step, tol)?;
        let path = if forward.stopped.is_none() {
            let curve = fit_stations(&forward.stations, tolerance, tol)?;
            let deviation = measure(&mesh, &curve, forward.stations.len(), tol)?;
            Some((vec![curve], deviation))
        } else if let Some(te) = inward(ce, ne) {
            let backward = march(&mesh, ce, te, c0, step, tol)?;
            if backward.stopped.is_none() {
                let mut stations = backward.stations;
                stations.reverse();
                let curve = fit_stations(&stations, tolerance, tol)?;
                let deviation = measure(&mesh, &curve, stations.len(), tol)?;
                Some((vec![curve], deviation))
            } else {
                cornered(&mesh, &forward, &backward, tolerance, tol)?
            }
        } else {
            None
        };
        match path {
            Some((curves, deviation)) if deviation <= tolerance => {
                return build(model, curves, start, end, deviation, tol);
            }
            Some((_, deviation)) => last = Ok(deviation),
            None => {
                if let Some(stopped) = forward.stopped {
                    last = Err(stopped);
                }
            }
        }
        step *= 0.5;
    }
    match last {
        Ok(deviation) => ogeom_bail!(
            NotDone,
            "the middle path reached a deviation of {deviation} against a tolerance of {tolerance}"
        ),
        Err(stopped) => ogeom_bail!(NotDone, "{stopped}"),
    }
}

/// The stations of one march, each with the march's tangent there and the
/// reach of its section (the farthest its outline stands from its
/// centroid), and why the march stopped short of the end, if it did.
struct Marched {
    stations: Vec<Point>,
    tangents: Vec<Vector>,
    reaches: Vec<f64>,
    stopped: Option<&'static str>,
}

/// March from `c0` along `t0` until the end centroid `ce` is within a step
/// and ahead. Near is not enough: a ring bent almost shut starts beside its
/// own end face. A section that finds no material, a station that does not
/// advance, or one that turns the tangent sharper than a smooth spine
/// would, stops the march there: past a sharp corner each of these is a
/// section run into the next leg, and the stations before it are kept.
fn march(
    mesh: &Slicer,
    c0: Point,
    t0: Vector,
    ce: Point,
    step: f64,
    tol: Tolerances,
) -> OgeomResult<Marched> {
    let first = mesh
        .section(c0 + t0 * (step * 1e-3), t0, tol)
        .map_or(0.0, |s| s.reach);
    let mut marched = Marched {
        stations: vec![c0],
        tangents: vec![t0],
        reaches: vec![first],
        stopped: None,
    };
    let (mut c, mut t) = (c0, t0);
    while c.distance(ce) > step * 1.25 || (ce - c).dot(t) < 0.5 * c.distance(ce) {
        if marched.stations.len() > MAX_STATIONS {
            ogeom_bail!(NotDone, "the middle path never reached the end face");
        }
        let Some(mut q) = mesh.section(c + t * step, t, tol) else {
            marched.stopped = Some("a section square to the middle path finds no material");
            return Ok(marched);
        };
        let mut turned = t;
        for _ in 0..8 {
            let chord = (q.centre - c).normalized(tol)?;
            turned = (chord * (2.0 * t.dot(chord)) - t).normalized(tol)?;
            let Some(next) = mesh.section(q.centre, turned, tol) else {
                break;
            };
            let moved = next.centre.distance(q.centre);
            q = next;
            if moved <= tol.confusion() {
                break;
            }
        }
        if (q.centre - c).dot(t) <= 0.0 || turned.dot(t) < SHARPEST_TURN {
            marched.stopped = Some("the middle path turned back on itself");
            return Ok(marched);
        }
        marched.stations.push(q.centre);
        marched.tangents.push(turned);
        marched.reaches.push(q.reach);
        (c, t) = (q.centre, turned);
    }
    marched.stations.push(ce);
    marched.tangents.push(t);
    marched.reaches.push(0.0);
    Ok(marched)
}

/// The path of two legs meeting at a sharp corner: the stations marched
/// from the start (`forward`) and from the end (`backward`), each cut back
/// to those whose section cannot reach the mitre, and the corner where
/// the two legs carried on straight along their last tangents meet. The
/// two curves, start to corner and corner to end, and the deviation
/// measured on both and at the corner; `None` where the legs do not meet
/// within `tolerance` ahead of both.
fn cornered(
    mesh: &Slicer,
    forward: &Marched,
    backward: &Marched,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Curve>, f64)>> {
    let (mut a, mut b) = (forward.stations.len(), backward.stations.len());
    let mut corner = None;
    for _ in 0..8 {
        if a == 0 || b == 0 {
            return Ok(None);
        }
        let (pa, da) = leg_end(forward, a, tolerance, tol)?;
        let (pb, db) = leg_end(backward, b, tolerance, tol)?;
        let Some(k) = meeting(pa, da, pb, db, tolerance) else {
            return Ok(None);
        };
        // A plane square to a leg at a distance s before the corner reaches
        // the plane halving it where the section reaches past s over the
        // tangent of half the turn.
        let cos_turn = da.dot(-db).clamp(-1.0, 1.0);
        if cos_turn <= -0.99 {
            return Ok(None);
        }
        let half_turn = ((1.0 - cos_turn) / (1.0 + cos_turn)).sqrt();
        let clear = |marched: &Marched, upto: usize| {
            (0..upto)
                .take_while(|&i| {
                    marched.stations[i].distance(k)
                        > marched.reaches[i] * half_turn * 1.05 + tolerance
                })
                .count()
        };
        let (na, nb) = (clear(forward, a).max(1), clear(backward, b).max(1));
        corner = Some((k, da, db));
        if (na, nb) == (a, b) {
            break;
        }
        (a, b) = (na, nb);
    }
    let Some((k, da, db)) = corner else {
        return Ok(None);
    };

    let mut first: Vec<Point> = forward.stations[..a].to_vec();
    first.push(k);
    let mut second: Vec<Point> = backward.stations[..b].to_vec();
    second.push(k);
    second.reverse();
    let legs = [
        fit_stations(&first, tolerance, tol)?,
        fit_stations(&second, tolerance, tol)?,
    ];
    // Each leg is measured where its sections stand clear of the corner;
    // the corner by the section halving it.
    let reach = forward.reaches[..a]
        .iter()
        .chain(&backward.reaches[..b])
        .fold(0.0_f64, |m, r| m.max(*r));
    let cos_turn = da.dot(-db).clamp(-1.0, 1.0);
    let clear = reach * ((1.0 - cos_turn) / (1.0 + cos_turn)).sqrt() * 1.05 + tolerance;
    let mut worst = 0.0_f64;
    for (leg, stations) in legs.iter().zip([a, b]) {
        worst = worst.max(measure_where(
            mesh,
            leg,
            stations,
            |p| p.distance(k) > clear,
            tol,
        )?);
    }
    let Some(halving) = mesh.section(k, da - db, tol) else {
        return Ok(None);
    };
    worst = worst.max(halving.centre.distance(k));
    Ok(Some((legs.to_vec(), worst)))
}

/// Where a march's first `upto` stations leave off, and the direction they
/// leave in: the end of the curve fitted through them, which on a straight
/// leg is the line from the first station to the last.
fn leg_end(
    marched: &Marched,
    upto: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<(Point, Vector)> {
    if upto < 2 {
        return Ok((marched.stations[0], marched.tangents[0]));
    }
    let curve = fit_stations(&marched.stations[..upto], tolerance, tol)?;
    let (_, hi) = curve.domain();
    Ok((
        curve.point_at(hi, tol)?,
        curve.d1_at(hi, tol)?.normalized(tol)?,
    ))
}

/// Where the lines through `pa` along `da` and through `pb` along `db`
/// meet: the middle of their closest points, both ahead of their starts
/// and within `tolerance` of each other.
fn meeting(pa: Point, da: Vector, pb: Point, db: Vector, tolerance: f64) -> Option<Point> {
    let w = pa - pb;
    let (aa, ab, bb) = (da.dot(da), da.dot(db), db.dot(db));
    let (aw, bw) = (da.dot(w), db.dot(w));
    let det = aa * bb - ab * ab;
    if det <= 1e-12 * aa * bb {
        return None;
    }
    let s = (ab * bw - bb * aw) / det;
    let u = (aa * bw - ab * aw) / det;
    if s < 0.0 || u < 0.0 {
        return None;
    }
    let (qa, qb) = (pa + da * s, pb + db * u);
    (qa.distance(qb) <= tolerance).then(|| qa.lerp(qb, 0.5))
}

/// A straight line where the stations are collinear, else a fitted spline.
fn fit_stations(stations: &[Point], tolerance: f64, tol: Tolerances) -> OgeomResult<Curve> {
    let (first, last) = (stations[0], stations[stations.len() - 1]);
    let axis = (last - first).normalized(tol)?;
    let off_line = stations
        .iter()
        .map(|p| {
            let d = *p - first;
            (d - axis * d.dot(axis)).magnitude()
        })
        .fold(0.0_f64, f64::max);
    if off_line <= tolerance * 0.05 {
        return Ok(Curve::Line(LineCurve::segment(first, last, tol)?));
    }
    let fitted = ogeom_geom::fit::fit_points(stations, 3, tolerance * 0.1, tol)?;
    Ok(Curve::BSpline(fitted.curve))
}

/// The largest distance from a section's centroid to the curve, cut square
/// to the curve at points between the stations.
fn measure(mesh: &Slicer, curve: &Curve, stations: usize, tol: Tolerances) -> OgeomResult<f64> {
    measure_where(mesh, curve, stations, |_| true, tol)
}

/// As [`measure`], at the points of the curve `keep` holds.
fn measure_where(
    mesh: &Slicer,
    curve: &Curve,
    stations: usize,
    keep: impl Fn(Point) -> bool,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let (lo, hi) = curve.domain();
    let samples = (stations * 2).max(8);
    let mut worst = 0.0_f64;
    // The ends stand on the end faces, whose centroids are the path's ends
    // by construction; only the inside is measured.
    for i in 1..samples {
        #[allow(clippy::cast_precision_loss)]
        let u = lo + (hi - lo) * (i as f64) / (samples as f64);
        let p = curve.point_at(u, tol)?;
        if !keep(p) {
            continue;
        }
        let d = curve.d1_at(u, tol)?.normalized(tol)?;
        let Some(q) = mesh.section(p, d, tol) else {
            ogeom_bail!(
                NotDone,
                "a section square to the fitted path finds no material"
            );
        };
        worst = worst.max(q.centre.distance(p));
    }
    Ok(worst)
}

/// The path's wire: its curves in order, each an edge from the end of the
/// one before.
fn build(
    model: &mut Model,
    curves: Vec<Curve>,
    start: &Shape,
    end: &Shape,
    deviation: f64,
    tol: Tolerances,
) -> OgeomResult<MiddlePath> {
    let mut edges = Vec::with_capacity(curves.len());
    let mut from: Option<Shape> = None;
    for curve in curves {
        let domain = curve.domain();
        let from_vertex = match from.take() {
            Some(v) => v,
            None => make_vertex(model, curve.point_at(domain.0, tol)?).shape,
        };
        let to_vertex = make_vertex(model, curve.point_at(domain.1, tol)?).shape;
        let edge = make_edge_between(model, curve, domain, &from_vertex, &to_vertex, tol)?.shape;
        edges.push(edge);
        from = Some(to_vertex);
    }
    let wire = make_wire(model, &edges, tol)?.shape;
    let mut history = History::new();
    history.generate(start, wire.clone());
    history.generate(end, wire.clone());
    Ok(MiddlePath {
        built: Built::new(wire, history),
        deviation,
    })
}

/// A face's area centroid, its mean outward normal and its area, off its
/// mesh.
fn face_centroid(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<(Point, Vector, f64)> {
    let mesh = triangulate_face(model, face, deflection, tol)?;
    let (mut weighted, mut normal, mut area) = (Vector::ZERO, Vector::ZERO, 0.0);
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let n = (b - a).cross(c - a) * 0.5;
        let da = n.magnitude();
        weighted += (a.to_vector() + b.to_vector() + c.to_vector()) * (da / 3.0);
        normal += n;
        area += da;
    }
    if area <= 0.0 {
        ogeom_bail!(Construction, "an end face has no area");
    }
    Ok((
        Point::from_vector(weighted * (1.0 / area)),
        normal.normalized(tol)?,
        area,
    ))
}

/// A closed triangle mesh, cut by planes.
struct Slicer {
    mesh: Triangulation,
    closed: bool,
}

impl Slicer {
    fn new(mesh: Triangulation) -> Self {
        let closed = !mesh.triangles.is_empty() && mesh.is_closed();
        Self { mesh, closed }
    }

    /// The region the plane through `origin` square to `normal` cuts from
    /// the solid, taking the outer loop round `origin` with the loops
    /// inside that as holes: its centroid, and the farthest the outer loop
    /// stands from it. `None` where no loop holds `origin`: the point is
    /// not inside the material.
    fn section(&self, origin: Point, normal: Vector, tol: Tolerances) -> Option<Section> {
        let n = normal.normalized(tol).ok()?;
        let helper = if n.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let e1 = n.cross(helper).normalized(tol).ok()?;
        let e2 = n.cross(e1);
        let flat = |p: Point| {
            let d = p - origin;
            (d.dot(e1), d.dot(e2))
        };

        let positions = &self.mesh.positions;
        let side: Vec<f64> = positions.iter().map(|p| (*p - origin).dot(n)).collect();
        // A vertex exactly on the plane counts as above it, so every
        // crossed triangle has exactly two crossed sides.
        let above = |i: u32| side[i as usize] >= 0.0;
        let crossing = |a: u32, b: u32| {
            let (da, db) = (side[a as usize], side[b as usize]);
            let s = da / (da - db);
            positions[a as usize].lerp(positions[b as usize], s)
        };
        let key = |a: u32, b: u32| if a < b { (a, b) } else { (b, a) };

        // Each crossed triangle gives a segment from one crossed side to
        // the other, oriented so that loops run consistently.
        let mut next: HashMap<(u32, u32), (u32, u32)> = HashMap::new();
        let mut points: HashMap<(u32, u32), Point> = HashMap::new();
        for t in &self.mesh.triangles {
            let ups = t.iter().filter(|&&i| above(i)).count();
            if ups == 0 || ups == 3 {
                continue;
            }
            let mut sides = [(0u32, 0u32); 2];
            let mut found = 0;
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if above(a) != above(b) && found < 2 {
                    // The side leaving the upper half first, so the segment
                    // runs the same way round every loop.
                    sides[found] = (a, b);
                    found += 1;
                }
            }
            let (s0, s1) = if above(sides[0].0) {
                (sides[0], sides[1])
            } else {
                (sides[1], sides[0])
            };
            points.insert(key(s0.0, s0.1), crossing(s0.0, s0.1));
            points.insert(key(s1.0, s1.1), crossing(s1.0, s1.1));
            next.insert(key(s0.0, s0.1), key(s1.0, s1.1));
        }
        if next.is_empty() {
            return None;
        }

        // Chain into loops.
        let mut loops: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut seen: HashMap<(u32, u32), ()> = HashMap::new();
        let mut starts: Vec<(u32, u32)> = next.keys().copied().collect();
        starts.sort_unstable();
        for s in starts {
            if seen.contains_key(&s) {
                continue;
            }
            let mut ring = Vec::new();
            let mut at = s;
            loop {
                seen.insert(at, ());
                ring.push(flat(points[&at]));
                match next.get(&at) {
                    Some(&n) if n == s => break,
                    Some(&n) if !seen.contains_key(&n) => at = n,
                    _ => return None,
                }
            }
            if ring.len() >= 3 {
                loops.push(ring);
            }
        }

        let measured: Vec<((f64, f64), f64)> = loops.iter().map(|l| area_centroid(l)).collect();
        // The outer loop: the largest holding the origin. A plane through
        // a point outside the solid may still cut it elsewhere, and that
        // cut is not this station's section.
        let outer = (0..loops.len())
            .filter(|&i| contains(&loops[i], (0.0, 0.0)))
            .max_by(|&a, &b| measured[a].1.abs().total_cmp(&measured[b].1.abs()))?;
        let (mut sx, mut sy, mut sa) = {
            let ((x, y), a) = measured[outer];
            (x * a.abs(), y * a.abs(), a.abs())
        };
        for (i, ring) in loops.iter().enumerate() {
            if i != outer && contains(&loops[outer], ring[0]) {
                let ((x, y), a) = measured[i];
                sx -= x * a.abs();
                sy -= y * a.abs();
                sa -= a.abs();
            }
        }
        if sa <= 0.0 {
            return None;
        }
        let (cx, cy) = (sx / sa, sy / sa);
        let reach = loops[outer]
            .iter()
            .fold(0.0_f64, |m, &(x, y)| m.max((x - cx).hypot(y - cy)));
        Some(Section {
            centre: origin + e1 * cx + e2 * cy,
            reach,
        })
    }
}

/// A plane section's centroid and the farthest its outline stands from it.
struct Section {
    centre: Point,
    reach: f64,
}

/// A polygon's area centroid and signed area.
fn area_centroid(ring: &[(f64, f64)]) -> ((f64, f64), f64) {
    let (mut a, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for i in 0..ring.len() {
        let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
        let w = p.0 * q.1 - q.0 * p.1;
        a += w;
        cx += (p.0 + q.0) * w;
        cy += (p.1 + q.1) * w;
    }
    a *= 0.5;
    if a == 0.0 {
        return (ring[0], 0.0);
    }
    ((cx / (6.0 * a), cy / (6.0 * a)), a)
}

/// Whether a point lies inside a polygon, by crossing parity.
fn contains(ring: &[(f64, f64)], p: (f64, f64)) -> bool {
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a.1 > p.1) != (b.1 > p.1) {
            let x = a.0 + (p.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
            if x > p.0 {
                inside = !inside;
            }
        }
    }
    inside
}
