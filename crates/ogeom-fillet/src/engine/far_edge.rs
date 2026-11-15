//! A concave round whose ball rests on the far edge of the narrow face
//! beside it, all along its seam (ADR 0001, section 13, question 8, and
//! amendment "M5d, a round resting on a far edge").
//!
//! Where two solids are fused at a slant, the seam between them is a
//! concave straight edge `E` between planar faces `F1` and `F2`, and `F2`
//! is often a strip running from `E` to its own convex far edge `E'`.
//! Where the strip is narrower than the ball's setback, the ball touching
//! `F1` cannot touch `F2`: it rests on `E'`. Its centre stands a radius off
//! `F1`'s plane and a radius from `E'`'s line: in the plane square to `E`
//! at each station, on the line a radius over `F1` where it meets the
//! cylinder about `E'`, the meeting nearer the ball's side (coming in over
//! `F1` toward the seam, it meets `E'` there first). The ball touches `F1`
//! at its centre's foot on the plane and `E'` at its centre's foot on the
//! line, both on its characteristic circle, so the round is the balls'
//! arcs between the two, fitted along the stations.
//!
//! The fill is the region between the seam and the round, bounded by the
//! two faces' planes (`E'` lies on `F2`'s), closed in the arcs' planes past
//! the seam's ends and cut back to the planes of the faces the seam ends
//! at, as the engine ends every round on its end faces; it is fused into
//! the solid in one boolean (ADR 0001, amendment "M5c"). Every case this
//! does not prove (the ball touching the strip anywhere, a face within the
//! ball's reach, an end met by other than one planar face) is handed back.

use ogeom_algo::{Built, PointDistance};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, LineCurve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, read};
use super::fallback::{Sense, Tool};
use super::trimmed::{EndCap, Section, Strip, faces_of, open_tool_over, outward};
use super::{Id, id};

/// How close the fit holds the round, its contact and its spine to their
/// exact stations.
const FIT: f64 = 1e-8;

/// The stations the run is first fitted at; doubled while the fit misses.
const FIRST: usize = 33;

/// The most stations the run is fitted at.
const MOST: usize = 2049;

/// The points each ball's arc is read at.
const ACROSS: usize = 65;

/// The stations along the seam at which the ball is proved to rest on the
/// far edge alone, clear of every other face.
const PROOF: usize = 129;

/// The rounds of `asked` (edge and radius) on `solid` where each is a
/// concave straight seam whose ball rests on the far edge of the narrow
/// face beside it all along; `None` where any is not, two share a vertex,
/// or the tools or the boolean will not build.
///
/// # Errors
///
/// As reading the solid and building the tools, other than refusals,
/// which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if asked.is_empty() {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut ends: Vec<Id> = Vec::with_capacity(2 * asked.len());
    let mut runs = Vec::with_capacity(asked.len());
    for (edge, radius) in asked {
        let e = id(edge);
        let Some(c) = curved.edges.get(&e) else {
            return Ok(None);
        };
        // Two runs meeting at a vertex close a corner not built here.
        if ends.contains(&c.ends.0) || ends.contains(&c.ends.1) {
            return Ok(None);
        }
        ends.extend([c.ends.0, c.ends.1]);
        let Some(run) = run_of(model, &curved, e, *radius, tol)? else {
            return Ok(None);
        };
        runs.push(run);
    }
    let mut tools = Vec::with_capacity(runs.len());
    for run in &runs {
        let Some(tool) = tool_of(model, &curved, run, tol)? else {
            return Ok(None);
        };
        tools.push(tool);
    }
    let edges: Vec<Shape> = asked.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(model, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    match super::guards::solid_fits(model, &built.shape, tol) {
        Ok(()) => Ok(Some(built)),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// One seam's run: the seam, its two faces (`F1` the touched, `F2` the
/// narrow), the far edge's line, and what places the ball.
struct Run {
    edge: Id,
    faces: [usize; 2],
    radius: f64,
    /// The seam's start and unit direction, and its length.
    start: Point,
    along: Vector,
    length: f64,
    /// `F1`'s outward normal, and the unit direction into `F1` from the
    /// seam, square to it.
    n1: Vector,
    into1: Vector,
    /// A point of the far edge's line and its unit direction.
    far_at: Point,
    far_along: Vector,
    /// Each end face's plane, as its outward normal and offset: the fill
    /// keeps `n . p <= d`.
    caps: [(Vector, f64); 2],
    /// The stations whose arcs' planes cap the tool, past the end faces'
    /// planes.
    span: (f64, f64),
}

impl Run {
    /// Whether `p` stands within both end faces' planes, on the fill's
    /// side.
    fn inside(&self, p: Point, eps: f64) -> bool {
        self.caps
            .iter()
            .all(|(n, d)| n.dot(p - Point::ORIGIN) <= d + eps)
    }

    /// The ball's centre at station `s` along the seam: on the line a
    /// radius over `F1` in the plane square to the seam there, where it
    /// meets the cylinder about the far edge, the meeting further into
    /// `F1`. `None` where they do not meet.
    fn centre(&self, s: f64) -> Option<Point> {
        let o = self.start + self.along * s + self.n1 * self.radius;
        let d = self.far_along;
        let g = o - self.far_at;
        let (gp, hp) = (g - d * g.dot(d), self.into1 - d * self.into1.dot(d));
        let a = hp.dot(hp);
        let b = gp.dot(hp);
        let c = gp.dot(gp) - self.radius * self.radius;
        let disc = b * b - a * c;
        if a <= f64::EPSILON || disc < 0.0 {
            return None;
        }
        Some(o + self.into1 * ((-b + disc.sqrt()) / a))
    }

    /// The ball's touch on `F1` and on the far edge, for centre `c`.
    fn feet(&self, c: Point) -> (Point, Point) {
        let on_far = self.far_at + self.far_along * (c - self.far_at).dot(self.far_along);
        (c - self.n1 * self.radius, on_far)
    }

    /// The centre's path's tangent at station `s`.
    fn tangent(&self, s: f64) -> Option<Vector> {
        let h = self.length * 1e-5;
        let (p, q) = (self.centre(s - h)?, self.centre(s + h)?);
        let d = (q - p) / (2.0 * h);
        (d.magnitude() > 0.0).then_some(d)
    }

    /// Where the arc's plane at station `s` crosses the seam's line.
    fn crossing(&self, s: f64) -> Option<Point> {
        let (c, d) = (self.centre(s)?, self.tangent(s)?);
        let den = d.dot(self.along);
        (den.abs() > f64::EPSILON * d.magnitude())
            .then(|| self.start + self.along * ((c - self.start).dot(d) / den))
    }

    /// The points of the ball's arc at station `s`, from `F1` to the far
    /// edge.
    fn arc(&self, s: f64, tol: Tolerances) -> OgeomResult<Option<Vec<Point>>> {
        let Some(c) = self.centre(s) else {
            return Ok(None);
        };
        let (p, q) = self.feet(c);
        let r = self.radius;
        let (p1, q1) = ((p - c) / r, (q - c) / r);
        let axis = p1.cross(q1);
        if axis.magnitude() <= tol.angular() {
            return Ok(None);
        }
        let axis = axis / axis.magnitude();
        let sweep = p1.dot(q1).clamp(-1.0, 1.0).acos();
        Ok(Some(
            (0..ACROSS)
                .map(|k| {
                    #[expect(clippy::cast_precision_loss, reason = "a few points")]
                    let t = sweep * (k as f64) / ((ACROSS - 1) as f64);
                    c + (p1 * t.cos() + axis.cross(p1) * t.sin()) * r
                })
                .collect(),
        ))
    }
}

/// The deflection faces are classified at: their trims drawn within the
/// confusion, so a probe a thousandth of an edge off a side reads clear of
/// it.
fn fine(tol: Tolerances) -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: tol.confusion(),
        ..Default::default()
    }
}

/// The unit direction into planar face `face` from the middle `at` of its
/// straight edge of length `length` along `along`, square to it in the
/// face's plane: the side a probe a thousandth of the edge's length off it
/// lands on the face. `None` where neither side does, or both.
pub(super) fn into_face(
    model: &Model,
    curved: &Curved,
    face: usize,
    (at, along, length): (Point, Vector, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    let Some(n) = outward(curved, face, at) else {
        return Ok(None);
    };
    let w = n.cross(along);
    let w = w / w.magnitude();
    let probe = length * 1e-3;
    let on = |p: Point| -> OgeomResult<bool> {
        Ok(
            ogeom_algo::classify_on_face(model, &curved.faces[face].shape, p, fine(tol), tol)?
                == ogeom_algo::Containment::In,
        )
    };
    Ok(match (on(at + w * probe)?, on(at - w * probe)?) {
        (true, false) => Some(w),
        (false, true) => Some(-w),
        _ => None,
    })
}

/// The run of `edge` rounded at `radius` with its ball resting on the far
/// edge of the narrow face beside it all along; `None` where it is not.
fn run_of(
    model: &Model,
    curved: &Curved,
    edge: Id,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let eps = tol.confusion();
    let c = &curved.edges[&edge];
    if c.ends.0 == c.ends.1 || !matches!(c.curve, Curve::Line(_)) {
        return Ok(None);
    }
    let (a, b) = (curved.points[&c.ends.0], curved.points[&c.ends.1]);
    let length = a.distance(b);
    if length <= eps {
        return Ok(None);
    }
    let along = (b - a) / length;
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let faces = [fa, fb];
    if faces
        .iter()
        .any(|f| !matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
    {
        return Ok(None);
    }
    let mid = a.midpoint(b);
    let (Some(na), Some(nb)) = (outward(curved, fa, mid), outward(curved, fb, mid)) else {
        return Ok(None);
    };
    // Concave: the second face rises out of the first's plane.
    let Some(into_b) = into_face(model, curved, fb, (mid, along, length), tol)? else {
        return Ok(None);
    };
    if into_b.dot(na) <= tol.angular() {
        return Ok(None);
    }
    for narrow in 0..2 {
        let (f1, f2) = (faces[1 - narrow], faces[narrow]);
        let (n1, n2) = if narrow == 0 { (nb, na) } else { (na, nb) };
        let (Some(into1), Some(into2)) = (
            into_face(model, curved, f1, (mid, along, length), tol)?,
            into_face(model, curved, f2, (mid, along, length), tol)?,
        ) else {
            continue;
        };
        // The ball's setback on the narrow face, where it would touch it.
        let turn = n1.dot(n2).clamp(-1.0, 1.0).acos();
        let setback = radius * (0.5 * turn).tan();
        let Some((far, at, far_along)) = far_edge(curved, f2, edge, mid, into2, tol)? else {
            continue;
        };
        if (at - mid).magnitude() >= setback - eps {
            continue;
        }
        let far_faces = faces_of(curved, far);
        let [ga, gb] = far_faces[..] else {
            continue;
        };
        let f3 = if ga == f2 { gb } else { ga };
        // A convex far edge: the face beyond it falls behind the strip's
        // plane.
        let g = &curved.edges[&far];
        let (gp, gq) = (
            g.curve.point_at(g.range.0, tol)?,
            g.curve.point_at(g.range.1, tol)?,
        );
        let Some(into3) = into_face(
            model,
            curved,
            f3,
            (gp.midpoint(gq), far_along, gp.distance(gq)),
            tol,
        )?
        else {
            continue;
        };
        if into3.dot(n2) >= -tol.angular() {
            continue;
        }
        // Each end of the seam met by one planar face besides the two.
        let mut caps = [(Vector::ZERO, 0.0); 2];
        let mut capped = true;
        for (k, v) in [c.ends.0, c.ends.1].into_iter().enumerate() {
            let at_v: Vec<usize> = (0..curved.faces.len())
                .filter(|fi| {
                    !faces.contains(fi)
                        && curved.faces[*fi].loops.iter().flatten().any(|(x, _)| {
                            let g = &curved.edges[x];
                            g.ends.0 == v || g.ends.1 == v
                        })
                })
                .collect();
            let p = curved.points[&v];
            match at_v[..] {
                [end] if matches!(curved.faces[end].surface, SurfaceGeometry::Plane(_)) => {
                    let Some(n) = outward(curved, end, p) else {
                        capped = false;
                        break;
                    };
                    caps[k] = (n, n.dot(p - Point::ORIGIN));
                }
                _ => {
                    capped = false;
                    break;
                }
            }
        }
        if !capped {
            continue;
        }
        let run = Run {
            edge,
            faces: [f1, f2],
            radius,
            start: a,
            along,
            length,
            n1,
            into1,
            far_at: at,
            far_along,
            caps,
            span: (0.0, 0.0),
        };
        let Some(span) = capped_span(&run, tol)? else {
            continue;
        };
        let run = Run { span, ..run };
        if proved(model, curved, &run, far, f2, f3, tol)? {
            return Ok(Some(run));
        }
    }
    Ok(None)
}

/// The far edge of the planar face `face`: the first of its straight edges
/// met crossing it from `mid` on `edge`, along `into`; with where it is
/// met and its unit direction. `None` where the first edge met is not
/// straight, or none is.
pub(super) fn far_edge(
    curved: &Curved,
    face: usize,
    edge: Id,
    mid: Point,
    into: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Id, Point, Vector)>> {
    let Some(normal) = outward(curved, face, mid) else {
        return Ok(None);
    };
    let mut best: Option<(f64, Id, Point, Vector)> = None;
    for (g, _) in curved.faces[face].loops.iter().flatten() {
        if *g == edge {
            continue;
        }
        let c = &curved.edges[g];
        let (p, q) = (
            c.curve.point_at(c.range.0, tol)?,
            c.curve.point_at(c.range.1, tol)?,
        );
        let w = q - p;
        let den = into.cross(w).dot(normal);
        if den.abs() <= tol.angular() * w.magnitude() {
            continue;
        }
        let s = (p - mid).cross(w).dot(normal) / den;
        let u = (p - mid).cross(into).dot(normal) / den;
        if s > tol.confusion() && (0.0..=1.0).contains(&u) && best.is_none_or(|b| s < b.0) {
            if !matches!(c.curve, Curve::Line(_)) {
                best = Some((s, *g, Point::ORIGIN, Vector::ZERO));
                continue;
            }
            best = Some((s, *g, mid + into * s, w / w.magnitude()));
        }
    }
    Ok(best.and_then(|(_, g, at, d)| (d.magnitude() > 0.0).then_some((g, at, d))))
}

/// Whether the ball rests on the far edge alone all along the run: at each
/// of [`PROOF`] stations between its caps whose arc the end faces' planes
/// keep any of, its touch on `F1` lies on that face and its touch on the
/// far edge within the edge (where the planes keep them), its centre
/// between the narrow face's and the beyond face's normals (resting on the
/// edge, touching neither face), and no face of the solid but `F1` nearer
/// than the radius, every face but the two at the far edge and the planar
/// faces wholly behind `F1`'s plane by the distance's Lipschitz bound up to
/// the next station (ADR 0001, section 13, question 8); and the arcs'
/// planes cross the seam in order.
fn proved(
    model: &Model,
    curved: &Curved,
    run: &Run,
    far: Id,
    f2: usize,
    f3: usize,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let far_edge = &curved.edges[&far];
    let (fp, fq) = (
        far_edge.curve.point_at(far_edge.range.0, tol)?,
        far_edge.curve.point_at(far_edge.range.1, tol)?,
    );
    let (flo, fhi) = {
        let (x, y) = (
            (fp - run.far_at).dot(run.far_along),
            (fq - run.far_at).dot(run.far_along),
        );
        (x.min(y), x.max(y))
    };
    let Some(n3) = outward(curved, f3, fp) else {
        return Ok(false);
    };
    let Some(n2) = outward(curved, f2, fp) else {
        return Ok(false);
    };
    // A quick look at three stations before the faces are read.
    for share in [0.0, 0.5, 1.0] {
        let Some(c) = run.centre(run.length * share) else {
            return Ok(false);
        };
        let (_, on_far) = run.feet(c);
        let w = c - on_far;
        if w.dot(n2) < -eps || w.dot(n3) < -eps {
            return Ok(false);
        }
    }
    let others: Vec<(usize, PointDistance)> = (0..curved.faces.len())
        .filter(|fi| *fi != run.faces[0])
        .map(|fi| Ok((fi, PointDistance::of(model, &curved.faces[fi].shape, tol)?)))
        .collect::<OgeomResult<_>>()?;
    let touched = &curved.faces[run.faces[0]].shape;
    // A planar face wholly behind the touched face's plane stands a radius
    // or more from every centre, which stands a radius in front of it.
    let behind: Vec<bool> = curved
        .faces
        .iter()
        .map(|f| {
            matches!(f.surface, SurfaceGeometry::Plane(_))
                && f.loops.iter().flatten().all(|(x, _)| {
                    let g = &curved.edges[x];
                    [g.ends.0, g.ends.1]
                        .iter()
                        .all(|v| (curved.points[v] - run.start).dot(run.n1) <= eps)
                })
        })
        .collect();
    let (lo, hi) = run.span;
    #[expect(clippy::cast_precision_loss, reason = "a few stations")]
    let step = (hi - lo) / (PROOF as f64);
    let mut last: Option<f64> = None;
    for i in 0..=PROOF {
        #[expect(clippy::cast_precision_loss, reason = "a few stations")]
        let s = lo + step * (i as f64);
        let (Some(c), Some(arc), Some(x)) = (run.centre(s), run.arc(s, tol)?, run.crossing(s))
        else {
            return Ok(false);
        };
        // The arcs' planes cross the seam's line in order.
        let at = (x - run.start).dot(run.along);
        if last.is_some_and(|l| at <= l) {
            return Ok(false);
        }
        last = Some(at);
        // Only what the end faces' planes keep is built.
        if !arc.iter().any(|p| run.inside(*p, eps)) {
            continue;
        }
        let (foot, on_far) = run.feet(c);
        if run.inside(on_far, eps) {
            let along_far = (on_far - run.far_at).dot(run.far_along);
            if along_far < flo + eps || along_far > fhi - eps {
                return Ok(false);
            }
        }
        // Resting on the edge: the centre in the edge's normal cone, off
        // the strip's plane and off the beyond face's, each by a margin.
        let w = c - on_far;
        let margin = run.radius * 1e-6;
        if w.dot(n2) < margin || w.dot(n3) < margin {
            return Ok(false);
        }
        if run.inside(foot, eps)
            && ogeom_algo::classify_on_face(model, touched, foot, fine(tol), tol)?
                == ogeom_algo::Containment::Out
        {
            return Ok(false);
        }
        // Every face but the two the ball rests against at the far edge
        // clear by the radius at every centre between this station and the
        // next: the distance moves by no more than the centre does.
        let next = run.centre(s + step).unwrap_or(c);
        let reach = 0.5 * c.distance(next) * 1.01;
        for (fi, faces) in &others {
            let held = if *fi == f2 || *fi == f3 || behind[*fi] {
                0.0
            } else {
                reach
            };
            if faces.nearest(c, tol)?.0 < run.radius - eps + held {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The stations whose arcs' planes cap the tool of `run`: at each end, the
/// first, stepping out from the seam's end half a radius at a time, whose
/// arc and seam crossing stand wholly outside that end face's plane.
/// `None` where none does within twenty radii.
fn capped_span(run: &Run, tol: Tolerances) -> OgeomResult<Option<(f64, f64)>> {
    let eps = tol.confusion();
    let r = run.radius;
    let outside = |s: f64, (n, d): (Vector, f64)| -> OgeomResult<bool> {
        let (Some(arc), Some(x)) = (run.arc(s, tol)?, run.crossing(s)) else {
            return Ok(false);
        };
        Ok(arc
            .iter()
            .chain(core::iter::once(&x))
            .all(|p| n.dot(*p - Point::ORIGIN) > d + eps))
    };
    let mut stations = [0.0; 2];
    for (k, sign) in [(0, -1.0), (1, 1.0)] {
        let base = if k == 0 { 0.0 } else { run.length };
        let mut found = None;
        for step in 1..=40 {
            let s = base + sign * r * 0.5 * f64::from(step);
            if outside(s, run.caps[k])? {
                found = Some(s);
                break;
            }
        }
        let Some(s) = found else {
            return Ok(None);
        };
        stations[k] = s;
    }
    Ok(Some((stations[0], stations[1])))
}

/// The fill tool of `run`: the strip's closed tool between the arcs'
/// planes past the seam's ends, cut back to the end faces' planes.
fn tool_of(
    model: &mut Model,
    curved: &Curved,
    run: &Run,
    tol: Tolerances,
) -> OgeomResult<Option<Tool>> {
    let eps = tol.confusion();
    let r = run.radius;
    let (lo, hi) = run.span;
    let (Some(x0), Some(x1)) = (run.crossing(lo), run.crossing(hi)) else {
        return Ok(None);
    };
    let (Some(d0), Some(d1)) = (run.tangent(lo), run.tangent(hi)) else {
        return Ok(None);
    };
    // The stations run a little past the caps, for the tool to find each
    // curve's crossing there.
    let past = 0.05 * (hi - lo);
    let mut n = FIRST;
    let (surface, spine, contact) = loop {
        let mut rows = Vec::with_capacity(n);
        let mut centres = Vec::with_capacity(n);
        let mut feet = Vec::with_capacity(n);
        for i in 0..n {
            #[expect(clippy::cast_precision_loss, reason = "a few stations")]
            let s = lo - past + (hi - lo + 2.0 * past) * (i as f64) / ((n - 1) as f64);
            let (Some(c), Some(arc)) = (run.centre(s), run.arc(s, tol)?) else {
                return Ok(None);
            };
            centres.push(c);
            feet.push(run.feet(c).0);
            rows.push(arc);
        }
        let surface = ogeom_geom::fit::fit_surface_grid_chordal(&rows, 3, FIT, tol)?;
        let spine = ogeom_geom::fit::fit_points(&centres, 3, FIT, tol)?;
        let contact = ogeom_geom::fit::fit_points(&feet, 3, FIT, tol)?;
        if surface.met && spine.met && contact.met {
            break (surface.curve, spine.curve, contact.curve);
        }
        if n >= MOST {
            return Ok(None);
        }
        n = 2 * n - 1;
    };
    // The far edge's contact, its line from the first station's foot to
    // the last's.
    let (Some(c0), Some(c1)) = (run.centre(lo - past), run.centre(hi + past)) else {
        return Ok(None);
    };
    let (q0, q1) = (run.feet(c0).1, run.feet(c1).1);
    if q0.distance(q1) <= eps {
        return Ok(None);
    }
    let strip = Strip {
        edge: run.edge,
        faces: run.faces,
        contacts: [
            Curve::BSpline(contact),
            LineCurve::segment(q0, q1, tol)?.into(),
        ],
        surface: SurfaceGeometry::BSpline(surface),
        section: Section::Arc {
            spine: Curve::BSpline(spine),
        },
    };
    let line: Curve = LineCurve::segment(x0, x1, tol)?.into();
    let range = (0.0, x0.distance(x1));
    let caps = [
        EndCap::Plane { at: x0, normal: d0 },
        EndCap::Plane { at: x1, normal: d1 },
    ];
    let Some(tool) = open_tool_over(model, curved, &strip, Some((&line, range)), caps, tol)? else {
        return Ok(None);
    };
    if tool.sense != Sense::Fuse {
        return Ok(None);
    }
    // Cut back to the end faces' planes, within a box well past the tool.
    let mut planes: Vec<(Vector, f64)> = run.caps.to_vec();
    let mut box_lo = Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut box_hi = -box_lo;
    for p in [x0, x1, c0, c1, q0, q1] {
        let v = p - Point::ORIGIN;
        box_lo = Vector::new(box_lo.x.min(v.x), box_lo.y.min(v.y), box_lo.z.min(v.z));
        box_hi = Vector::new(box_hi.x.max(v.x), box_hi.y.max(v.y), box_hi.z.max(v.z));
    }
    let margin = 4.0 * r + 10.0 * eps;
    planes.extend([
        (Vector::X, box_hi.x + margin),
        (-Vector::X, -(box_lo.x - margin)),
        (Vector::Y, box_hi.y + margin),
        (-Vector::Y, -(box_lo.y - margin)),
        (Vector::Z, box_hi.z + margin),
        (-Vector::Z, -(box_lo.z - margin)),
    ]);
    let Some(slab) = super::fallback::polytope(model, &planes, tol)? else {
        return Ok(None);
    };
    let Ok(cut) = ogeom_bool::common(model, &tool.solid, &slab, tol) else {
        return Ok(None);
    };
    let mut faces = Vec::new();
    for (face, from) in &tool.faces {
        for made in cut
            .history
            .modified(face)
            .iter()
            .chain(cut.history.generated(face))
        {
            faces.push((made.clone(), from.clone()));
        }
        if !cut.history.is_deleted(face)
            && cut.history.modified(face).is_empty()
            && cut.history.generated(face).is_empty()
        {
            faces.push((face.clone(), from.clone()));
        }
    }
    Ok(Some(Tool {
        solid: cut.shape,
        sense: Sense::Fuse,
        faces,
    }))
}
