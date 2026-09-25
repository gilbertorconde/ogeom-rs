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

use ogeom_algo::{Built, History, make_edge, make_wire};
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

/// The centre line of a pipe-like `solid`, from its face `start` to its
/// face `end`.
///
/// Each point of the path is the centroid of the solid's section by a
/// plane square to the path there. The path starts at the centroid of
/// `start` and ends at the centroid of `end`, and it is a straight edge
/// where every station lies on one line and a fitted spline otherwise.
/// `tolerance` bounds the reported deviation; the solid is meshed at a
/// quarter of it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `solid` is not a solid, either face is not one of its faces, the two
/// faces are the same, the tolerance is not a positive distance, or the
/// mesh does not close. [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone)
/// if a section square to the march finds no material (the solid is not
/// pipe-like between the faces), the march never reaches `end`, or the
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
    let (ce, _, _) = face_centroid(model, end, deflection, tol)?;
    let size = (area0 / core::f64::consts::PI).sqrt();
    if size <= tol.confusion() {
        ogeom_bail!(Construction, "the start face has no area to start from");
    }

    // Into the solid: the side where a nearby plane cuts material round
    // the start face's centroid.
    let probe = size * 0.05;
    let t0 = if mesh.section(c0 + n0 * probe, n0, tol).is_some() {
        n0
    } else if mesh.section(c0 - n0 * probe, n0, tol).is_some() {
        -n0
    } else {
        ogeom_bail!(NotDone, "no material lies behind the start face");
    };

    let mut step = size;
    let mut last = None;
    for _ in 0..=MAX_REFINEMENTS {
        let stations = march(&mesh, c0, t0, ce, step, tol)?;
        let curve = fit_stations(&stations, tolerance, tol)?;
        let deviation = measure(&mesh, &curve, stations.len(), tol)?;
        if deviation <= tolerance {
            return build(model, curve, start, end, deviation, tol);
        }
        last = Some(deviation);
        step *= 0.5;
    }
    ogeom_bail!(
        NotDone,
        "the middle path reached a deviation of {} against a tolerance of {tolerance}",
        last.unwrap_or(f64::INFINITY)
    )
}

/// March from `c0` along `t0` until the end centroid `ce` is within a step
/// and ahead. Near is not enough: a ring bent almost shut starts beside its
/// own end face.
fn march(
    mesh: &Slicer,
    c0: Point,
    t0: Vector,
    ce: Point,
    step: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    let mut stations = vec![c0];
    let (mut c, mut t) = (c0, t0);
    while c.distance(ce) > step * 1.25 || (ce - c).dot(t) < 0.5 * c.distance(ce) {
        if stations.len() > MAX_STATIONS {
            ogeom_bail!(NotDone, "the middle path never reached the end face");
        }
        let Some(mut q) = mesh.section(c + t * step, t, tol) else {
            ogeom_bail!(
                NotDone,
                "a section square to the path finds no material {} along from {:?}",
                step,
                c
            );
        };
        let mut turned = t;
        for _ in 0..8 {
            let chord = (q - c).normalized(tol)?;
            turned = (chord * (2.0 * t.dot(chord)) - t).normalized(tol)?;
            let Some(next) = mesh.section(q, turned, tol) else {
                break;
            };
            let moved = next.distance(q);
            q = next;
            if moved <= tol.confusion() {
                break;
            }
        }
        // A station that does not advance toward the end means the march
        // turned back on itself.
        if (q - c).dot(t) <= 0.0 {
            ogeom_bail!(NotDone, "the middle path turned back on itself");
        }
        stations.push(q);
        (c, t) = (q, turned);
    }
    stations.push(ce);
    Ok(stations)
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
    let (lo, hi) = curve.domain();
    let samples = (stations * 2).max(8);
    let mut worst = 0.0_f64;
    // The ends stand on the end faces, whose centroids are the path's ends
    // by construction; only the inside is measured.
    for i in 1..samples {
        #[allow(clippy::cast_precision_loss)]
        let u = lo + (hi - lo) * (i as f64) / (samples as f64);
        let p = curve.point_at(u, tol)?;
        let d = curve.d1_at(u, tol)?.normalized(tol)?;
        let Some(q) = mesh.section(p, d, tol) else {
            ogeom_bail!(
                NotDone,
                "a section square to the fitted path finds no material"
            );
        };
        worst = worst.max(q.distance(p));
    }
    Ok(worst)
}

fn build(
    model: &mut Model,
    curve: Curve,
    start: &Shape,
    end: &Shape,
    deviation: f64,
    tol: Tolerances,
) -> OgeomResult<MiddlePath> {
    let domain = curve.domain();
    let edge = make_edge(model, curve, domain, tol)?.shape;
    let wire = make_wire(model, &[edge], tol)?.shape;
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

    /// The centroid of the region the plane through `origin` square to
    /// `normal` cuts from the solid, taking the outer loop round `origin`
    /// with the loops inside that as holes. `None` where no loop holds
    /// `origin`: the point is not inside the material.
    fn section(&self, origin: Point, normal: Vector, tol: Tolerances) -> Option<Point> {
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
        Some(origin + e1 * (sx / sa) + e2 * (sy / sa))
    }
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
