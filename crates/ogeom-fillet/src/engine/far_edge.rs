//! A concave round whose ball, touching one face of its seam, rests across
//! from it on a far edge or a face past the narrow face beside it, all
//! along the seam (ADR 0001, section 13, question 8, and amendment "M5d, a
//! round resting on a far edge").
//!
//! Where two solids are fused at a slant, the seam between them is a
//! concave straight edge `E` between planar faces `F1` and `F2`, and `F2`
//! is often a strip running from `E` to its own convex far edge `E'`.
//! Where the strip is narrower than the ball's setback, the ball touching
//! `F1` cannot touch `F2`: it rests on `E'`, or, past it, on the far edge
//! of a strip beyond (a staircase), or comes down on the face past the
//! strips (a floor). Its centre stands a radius off `F1`'s plane and a
//! radius from the support: in the plane square to the seam at each
//! station, on the line a radius over `F1` where it meets the cylinder
//! about the far edge (the meeting nearer the ball's side: coming in over
//! `F1` toward the seam, it meets the edge there first) or the plane a
//! radius over the floor. The ball touches `F1` and its support at its
//! centre's feet, both on its characteristic circle, so the round is the
//! balls' arcs between the two, fitted along the stations.
//!
//! The fill is the region between the round and the crease where `F1`'s
//! plane meets the plane of the strip whose far edge holds the ball (or of
//! the floor), bounded by those two planes: the strips before it lie in
//! the solid or under the round, and go into the fill whole. It is closed
//! in the arcs' planes past the seam's ends and cut back to the planes of
//! the faces the seam ends at, as the engine ends every round on its end
//! faces, and fused into the solid in one boolean (ADR 0001, amendment
//! "M5c"). Every case this does not prove (the ball touching the strip
//! anywhere, a face within the ball's reach, an end met by other than one
//! planar face, a strip left standing) is handed back.

use ogeom_algo::{Built, PointDistance};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, CylinderSurface, LineCurve, SurfaceGeometry};
use ogeom_math::{Cylinder, Direction, Frame, Point, Vector};
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

/// How many times an interval between two proof stations is halved, at
/// most, to show a face clear of the balls over it.
const HALVINGS: u32 = 24;

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
    // The strips each fill takes in are gone whole.
    if runs.iter().any(|run| {
        run.taken
            .iter()
            .any(|f| !built.history.trace(&curved.faces[*f].shape).is_empty())
    }) {
        return Ok(None);
    }
    match super::guards::solid_fits(model, &built.shape, tol) {
        Ok(()) => Ok(Some(built)),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// What the ball rests on across the seam from the face it touches: the
/// line of a far edge, or a face's plane.
#[derive(Debug, Clone, Copy)]
enum Support {
    /// A far edge's line: a point of it and its unit direction.
    Edge { at: Point, along: Vector },
    /// A plane: a point of it and its outward unit normal.
    Plane { at: Point, normal: Vector },
}

/// One seam's run: the seam, the touched face `F1` and the face `F2` whose
/// plane closes the fill on the other side, the support the ball rests on,
/// the strips it takes into its fill, and what places the ball.
struct Run {
    edge: Id,
    faces: [usize; 2],
    radius: f64,
    /// The crease the fill closes on, where `F1`'s plane meets `F2`'s (the
    /// seam itself where `F2` is the seam's narrow face): its start and unit
    /// direction, and its length between the seam's ends.
    start: Point,
    along: Vector,
    length: f64,
    /// `F1`'s outward normal, and the unit direction into `F1` from the
    /// seam, square to it.
    n1: Vector,
    into1: Vector,
    support: Support,
    /// The faces between the seam and `F2` (the strips of a staircase),
    /// which the fill takes in whole.
    taken: Vec<usize>,
    /// Each end face's plane, as its outward normal and offset: the fill
    /// keeps `n . p <= d`; and the end faces themselves (a plane a boolean
    /// left as several faces is all of them).
    caps: [(Vector, f64); 2],
    end_faces: Vec<usize>,
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

    /// The ball's centre at station `s` along the crease: on the line a
    /// radius over `F1` in the plane square to the crease there, where it
    /// meets the cylinder about the far edge (the meeting further into
    /// `F1`) or the plane a radius over the support's. `None` where they
    /// do not meet.
    fn centre(&self, s: f64) -> Option<Point> {
        let o = self.start + self.along * s + self.n1 * self.radius;
        match self.support {
            Support::Edge { at, along: d } => {
                let g = o - at;
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
            Support::Plane { at, normal } => {
                let den = self.into1.dot(normal);
                if den.abs() <= f64::EPSILON {
                    return None;
                }
                Some(o + self.into1 * ((self.radius - (o - at).dot(normal)) / den))
            }
        }
    }

    /// The ball's touch on `F1` and on its support, for centre `c`.
    fn feet(&self, c: Point) -> (Point, Point) {
        let other = match self.support {
            Support::Edge { at, along } => at + along * (c - at).dot(along),
            Support::Plane { normal, .. } => c - normal * self.radius,
        };
        (c - self.n1 * self.radius, other)
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
/// What a run's ball rests against across from `F1`, for its proof: the far
/// edge with the strip's and the beyond face's outward normals, or the face
/// whose plane it touches.
#[derive(Debug, Clone, Copy)]
enum Held {
    Edge {
        far: Id,
        faces: [usize; 2],
        normals: [Vector; 2],
    },
    Face {
        face: usize,
    },
}

/// The most strips walked out from the seam.
const MOST_STRIPS: usize = 6;

/// The run of `edge` rounded at `radius` where its ball, touching one of
/// the seam's faces, cannot touch the other, narrower than its setback:
/// walking the strips out from the seam across each one's far edge, the
/// first the ball is proved to rest on all along, a strip's convex far
/// edge (the strips before it taken into the fill) or the plane of the
/// face past a strip (the strips up to it taken). `None` where none is.
fn run_of(
    model: &Model,
    curved: &Curved,
    edge: Id,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Run>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
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
    let planar = |f: usize| matches!(curved.faces[f].surface, SurfaceGeometry::Plane(_));
    if !planar(fa) || !planar(fb) {
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
    if into_b.dot(na) <= ang {
        return Ok(None);
    }
    // Each end of the seam met by one planar face besides the two.
    let mut caps = [(Vector::ZERO, 0.0); 2];
    let mut end_faces: Vec<usize> = Vec::new();
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
        // One plane: a boolean may leave it as several faces.
        let p = curved.points[&v];
        let Some(&first) = at_v.first() else {
            return Ok(None);
        };
        let (true, Some(n)) = (planar(first), outward(curved, first, p)) else {
            return Ok(None);
        };
        for end in &at_v {
            if !planar(*end)
                || outward(curved, *end, p)
                    .is_none_or(|m| m.cross(n).magnitude() > ang || m.dot(n) < 0.0)
            {
                return Ok(None);
            }
        }
        caps[k] = (n, n.dot(p - Point::ORIGIN));
        end_faces.extend(at_v);
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
        // A run closing on the plane of `flap`, held by `support`, the
        // strips `taken` going into its fill.
        let attempt = |flap: usize, support: Support, held: Held, taken: &[usize]| {
            let Some(nf) = outward(curved, flap, mid) else {
                return Ok(None);
            };
            // The crease where `F1`'s plane meets the flap's, along the seam.
            let line = n1.cross(nf);
            if line.magnitude() <= ang || line.cross(along).magnitude() > ang * line.magnitude() {
                return Ok(None);
            }
            let den = into1.dot(nf);
            if den.abs() <= ang {
                return Ok(None);
            }
            let q = match support {
                Support::Edge { at, .. } | Support::Plane { at, .. } => at,
            };
            let start = a + into1 * ((q - a).dot(nf) / den);
            let run = Run {
                edge,
                faces: [f1, flap],
                radius,
                start,
                along,
                length,
                n1,
                into1,
                support,
                taken: taken.to_vec(),
                caps,
                end_faces: end_faces.clone(),
                span: (0.0, 0.0),
            };
            let Some(span) = capped_span(&run, tol)? else {
                return Ok(None);
            };
            let run = Run { span, ..run };
            Ok::<_, ogeom_core::OgeomError>(proved(model, curved, &run, held, tol)?.then_some(run))
        };
        let (mut face, mut near, mut at, mut into) = (f2, edge, mid, into2);
        let mut taken: Vec<usize> = Vec::new();
        while taken.len() <= MOST_STRIPS {
            let Some((far, hit, far_along)) = far_edge(curved, face, near, at, into, tol)? else {
                break;
            };
            if taken.is_empty() && (hit - mid).magnitude() >= setback - eps {
                break;
            }
            let [ga, gb] = faces_of(curved, far)[..] else {
                break;
            };
            let g = if ga == face { gb } else { ga };
            if g == f1 || g == face || taken.contains(&g) {
                break;
            }
            let (Some(nf), Some(ng)) = (outward(curved, face, hit), outward(curved, g, hit)) else {
                break;
            };
            let gc = &curved.edges[&far];
            let (gp, gq) = (
                gc.curve.point_at(gc.range.0, tol)?,
                gc.curve.point_at(gc.range.1, tol)?,
            );
            let into_g = into_face(
                model,
                curved,
                g,
                (gp.midpoint(gq), far_along, gp.distance(gq)),
                tol,
            )?;
            // Resting on a convex far edge: the face beyond falls behind the
            // strip's plane.
            if into_g.is_some_and(|i| i.dot(nf) < -ang) {
                let held = Held::Edge {
                    far,
                    faces: [face, g],
                    normals: [nf, ng],
                };
                let support = Support::Edge {
                    at: hit,
                    along: far_along,
                };
                if let Some(run) = attempt(face, support, held, &taken)? {
                    return Ok(Some(run));
                }
            }
            taken.push(face);
            let (true, Some(into_g)) = (planar(g), into_g) else {
                break;
            };
            // Coming down on the face past the strips.
            let support = Support::Plane {
                at: hit,
                normal: ng,
            };
            if let Some(run) = attempt(g, support, Held::Face { face: g }, &taken)? {
                return Ok(Some(run));
            }
            (face, near, at, into) = (g, far, hit, into_g);
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

/// Whether the ball rests on its support alone all along the run: at each
/// of [`PROOF`] stations between its caps whose arc the end faces' planes
/// keep any of, its touch on `F1` lies on that face and its other touch
/// within the far edge or on the face it comes down on (where the planes
/// keep them); resting on a far edge, its centre between the strip's and
/// the beyond face's normals (touching neither face); and no face of the
/// solid but `F1` and the end faces (the fill is cut to their planes)
/// nearer than the radius, every face but those it rests against and the
/// planar faces wholly behind a plane it touches by the distance's
/// Lipschitz bound up to the next station, halving where it is close
/// (ADR 0001, section 13, question 8); and the arcs' planes cross the
/// crease in order.
fn proved(
    model: &Model,
    curved: &Curved,
    run: &Run,
    held: Held,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let margin = run.radius * 1e-6;
    // Within the far edge, and resting on it alone.
    let on_edge = |c: Point, on_far: Point| -> OgeomResult<bool> {
        let Held::Edge { far, normals, .. } = held else {
            return Ok(true);
        };
        let Support::Edge { at, along } = run.support else {
            return Ok(false);
        };
        let g = &curved.edges[&far];
        let (p, q) = (
            g.curve.point_at(g.range.0, tol)?,
            g.curve.point_at(g.range.1, tol)?,
        );
        let (x, y) = ((p - at).dot(along), (q - at).dot(along));
        let t = (on_far - at).dot(along);
        let w = c - on_far;
        // Within the edge, its ends included: an end ball may touch it at
        // its end, on an end face's plane, where the fill is cut.
        Ok(
            (!run.inside(on_far, eps) || (t >= x.min(y) - eps && t <= x.max(y) + eps))
                && w.dot(normals[0]) >= margin
                && w.dot(normals[1]) >= margin,
        )
    };
    // A quick look at three stations before the faces are read.
    for share in [0.0, 0.5, 1.0] {
        let Some(c) = run.centre(run.length * share) else {
            return Ok(false);
        };
        if !on_edge(c, run.feet(c).1)? {
            return Ok(false);
        }
    }
    let resting: Vec<usize> = match held {
        Held::Edge { faces, .. } => faces.to_vec(),
        Held::Face { face } => vec![face],
    };
    let others: Vec<(usize, PointDistance)> = (0..curved.faces.len())
        .filter(|fi| *fi != run.faces[0] && !run.end_faces.contains(fi))
        .map(|fi| Ok((fi, PointDistance::of(model, &curved.faces[fi].shape, tol)?)))
        .collect::<OgeomResult<_>>()?;
    // A planar face wholly behind a plane the ball touches stands a radius
    // or more from every centre, which stands a radius in front of it.
    let mut planes = vec![(run.start, run.n1)];
    if let Support::Plane { at, normal } = run.support {
        planes.push((at, normal));
    }
    let behind: Vec<bool> = curved
        .faces
        .iter()
        .map(|f| {
            matches!(f.surface, SurfaceGeometry::Plane(_))
                && planes.iter().any(|(o, n)| {
                    f.loops.iter().flatten().all(|(x, _)| {
                        let g = &curved.edges[x];
                        [g.ends.0, g.ends.1]
                            .iter()
                            .all(|v| (curved.points[v] - *o).dot(*n) <= eps)
                    })
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
        // The arcs' planes cross the crease's line in order.
        let at = (x - run.start).dot(run.along);
        if last.is_some_and(|l| at <= l) {
            return Ok(false);
        }
        last = Some(at);
        // Only what the end faces' planes keep is built.
        if !arc.iter().any(|p| run.inside(*p, eps)) {
            continue;
        }
        let (foot, other) = run.feet(c);
        if !on_edge(c, other)? {
            return Ok(false);
        }
        let on = |face: usize, p: Point| -> OgeomResult<bool> {
            Ok(!run.inside(p, eps)
                || ogeom_algo::classify_on_face(
                    model,
                    &curved.faces[face].shape,
                    p,
                    fine(tol),
                    tol,
                )? != ogeom_algo::Containment::Out)
        };
        if !on(run.faces[0], foot)? {
            return Ok(false);
        }
        if let Held::Face { face } = held
            && !on(face, other)?
        {
            return Ok(false);
        }
        // Every face but those the ball rests against clear by the radius
        // at every centre between this station and the next: the distance
        // moves by no more than the centre does.
        for (fi, faces) in &others {
            let d = faces.nearest(c, tol)?.0;
            if d < run.radius - eps {
                return Ok(false);
            }
            if !resting.contains(fi)
                && !behind[*fi]
                && i < PROOF
                && !clear_between(run, faces, (s, s + step), d, HALVINGS, tol)?
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Whether every centre of `run` between stations `s0` and `s1` stands a
/// radius or more from `faces`, `d0` its distance at `s0` (the one at `s1`
/// is checked at that station): the distance moves by no more than the
/// centre does, so a margin of half the chord between them clears the
/// interval; where it does not, the interval is halved, at most `depth`
/// times.
fn clear_between(
    run: &Run,
    faces: &PointDistance,
    (s0, s1): (f64, f64),
    d0: f64,
    depth: u32,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let (Some(c0), Some(c1)) = (run.centre(s0), run.centre(s1)) else {
        return Ok(false);
    };
    let d1 = faces.nearest(c1, tol)?.0;
    if d0.min(d1) < run.radius - eps {
        return Ok(false);
    }
    let room = 0.5 * c0.distance(c1) * 1.01;
    if d0.min(d1) >= run.radius - eps + room {
        return Ok(true);
    }
    if depth == 0 {
        return Ok(false);
    }
    let m = 0.5 * (s0 + s1);
    let Some(cm) = run.centre(m) else {
        return Ok(false);
    };
    let dm = faces.nearest(cm, tol)?.0;
    Ok(clear_between(run, faces, (s0, m), d0, depth - 1, tol)?
        && clear_between(run, faces, (m, s1), dm, depth - 1, tol)?)
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
    let (Some(c0), Some(c1)) = (run.centre(lo - past), run.centre(hi + past)) else {
        return Ok(None);
    };
    // The support's contact, its line from the first station's foot to
    // the last's.
    let (q0, q1) = (run.feet(c0).1, run.feet(c1).1);
    if q0.distance(q1) <= eps {
        return Ok(None);
    }
    let segment =
        |p: Point, q: Point| -> OgeomResult<Curve> { Ok(LineCurve::segment(p, q, tol)?.into()) };
    // A centre running on a line (a plane for support, or a far edge along
    // the crease) sweeps a cylinder, exact; any other, the balls' arcs
    // fitted along the stations.
    let straight = match run.support {
        Support::Plane { .. } => true,
        Support::Edge { along, .. } => along.cross(run.along).magnitude() <= tol.angular(),
    };
    let (surface, spine, contact): (SurfaceGeometry, Curve, Curve) = if straight {
        let u = (c1 - c0) / c0.distance(c1);
        let x = Direction::new(run.feet(c0).0 - c0, tol)?;
        let frame = Frame::new(c0, Direction::new(u, tol)?, x, tol)?;
        (
            CylinderSurface::new(Cylinder::new(frame, r, tol)?, (0.0, c0.distance(c1)))?.into(),
            segment(c0, c1)?,
            segment(run.feet(c0).0, run.feet(c1).0)?,
        )
    } else {
        let mut n = FIRST;
        loop {
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
                break (
                    SurfaceGeometry::BSpline(surface.curve),
                    Curve::BSpline(spine.curve),
                    Curve::BSpline(contact.curve),
                );
            }
            if n >= MOST {
                return Ok(None);
            }
            n = 2 * n - 1;
        }
    };
    let strip = Strip {
        edge: run.edge,
        faces: run.faces,
        contacts: [contact, segment(q0, q1)?],
        surface,
        section: Section::Arc { spine },
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
