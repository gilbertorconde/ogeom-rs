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

/// How many distances one interval's proof may measure. A face standing a
/// radius from every centre (tangent to all the balls) leaves no margin
/// anywhere, so without a budget the halving would walk its whole tree of
/// depth [`HALVINGS`]; past the budget the face is not shown clear.
const MEASURES: u32 = 1 << 17;

/// The rounds of `asked` (edge and radius) on `solid` where each is a
/// concave straight seam whose ball rests on the far edge of the narrow
/// face beside it all along (or past it, or carried onto it); `None` where
/// any is not, two share a vertex, or the tools or the boolean will not
/// build. Where the balls reach a face standing in the fill, the faces come
/// off a copy and the rounds are built there by `round`.
///
/// # Errors
///
/// As reading the solid and building the tools, other than refusals,
/// which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    round: super::defeatured::Rounder,
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
    // The asked edges that are no such run: rounded first by the engine,
    // where they share no face and no vertex with the runs.
    let mut others: Vec<(Shape, f64)> = Vec::new();
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
        let Some(pieces) = run_of(model, &curved, e, *radius, tol)? else {
            others.push((edge.clone(), *radius));
            continue;
        };
        runs.extend(pieces);
    }
    if runs.is_empty() {
        return Ok(None);
    }
    if !others.is_empty() {
        return beside_others(model, solid, &curved, &runs, &others, round, tol);
    }
    // Faces standing in the fill: the general route, those faces off a copy
    // and the rounds built there by the engine (the round keeps its
    // surface, the faces it reaches are trimmed or covered).
    let mut far: Vec<usize> = runs.iter().flat_map(|r| r.reaches.clone()).collect();
    if !far.is_empty() {
        far.sort_unstable();
        far.dedup();
        let mut near: Vec<usize> = runs
            .iter()
            .flat_map(|r| {
                r.faces
                    .iter()
                    .chain(&r.taken)
                    .chain(&r.end_faces)
                    .copied()
                    .collect::<Vec<_>>()
            })
            .filter(|f| !far.contains(f))
            .collect();
        near.sort_unstable();
        near.dedup();
        let shape = |f: &usize| curved.faces[*f].shape.clone();
        let (far, near): (Vec<Shape>, Vec<Shape>) = (
            far.iter().map(shape).collect(),
            near.iter().map(shape).collect(),
        );
        let mut scratch = model.clone();
        let Some(tools) =
            super::defeatured::tools(&mut scratch, solid, asked, &far, &near, round, tol)?
        else {
            return Ok(None);
        };
        let edges: Vec<Shape> = asked.iter().map(|(e, _)| e.clone()).collect();
        let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)?
        else {
            return Ok(None);
        };
        return match super::guards::solid_fits(&scratch, &built.shape, tol) {
            Ok(()) => {
                *model = scratch;
                Ok(Some(built))
            }
            Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
            Err(_) => Ok(None),
        };
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

/// The far-edge `runs` with the `others` of the request: the others rounded
/// first by `round`, then each run's fill fused into what that leaves.
/// Only where no run reaches a face and the others share no face with a
/// run (the vertices were kept apart as read): the fills then stand where
/// the others' rounds change nothing. `None` where any step declines.
fn beside_others(
    model: &mut Model,
    solid: &Shape,
    curved: &Curved,
    runs: &[Run],
    others: &[(Shape, f64)],
    round: super::defeatured::Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if runs.iter().any(|r| !r.reaches.is_empty()) {
        return Ok(None);
    }
    let held: Vec<usize> = runs
        .iter()
        .flat_map(|r| r.faces.iter().chain(&r.taken).copied().collect::<Vec<_>>())
        .collect();
    if others
        .iter()
        .any(|(e, _)| faces_of(curved, id(e)).iter().any(|f| held.contains(f)))
    {
        return Ok(None);
    }
    let mut scratch = model.clone();
    let Some(first) = round(&mut scratch, solid, others, tol)? else {
        return Ok(None);
    };
    let mut tools = Vec::with_capacity(runs.len());
    for run in runs {
        let Some(tool) = tool_of(&mut scratch, curved, run, tol)? else {
            return Ok(None);
        };
        tools.push(tool);
    }
    let Some(step) = super::fallback::apply(&mut scratch, &first.shape, &tools, tol)? else {
        return Ok(None);
    };
    // The fills' faces named for their seams through both steps.
    let mut history = first.history.then(&step.history);
    for tool in &tools {
        for (_, from) in &tool.faces {
            for made in step.history.generated(from) {
                if !history.generated(from).iter().any(|m| m.is_same(made)) {
                    history.generate(from, made.clone());
                }
            }
        }
    }
    let built = Built::new(step.shape, history);
    if others
        .iter()
        .map(|(e, _)| e)
        .chain(runs.iter().filter_map(|r| curved.shapes.get(&r.edge)))
        .any(|e| built.history.generated(e).is_empty())
        || runs.iter().any(|run| {
            run.taken
                .iter()
                .any(|f| !built.history.trace(&curved.faces[*f].shape).is_empty())
        })
    {
        return Ok(None);
    }
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
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
    /// A vertex the far edge ends at: the ball rolls round it.
    Point { at: Point },
    /// A far edge's line read round it: the centre on the cylinder about
    /// it where it meets the plane a radius over `F1`, `s` a radius times
    /// the angle round the line, for an edge steep enough that the centre
    /// runs back along the seam.
    Around { at: Point, along: Vector },
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
    /// planes, or where the ball changes support.
    span: (f64, f64),
    /// Faces of the solid the balls come within their radius of, standing
    /// in the fill (a post): the round keeps its surface and goes through
    /// the general route, the faces off a copy (ADR 0001, amendment
    /// "M5c").
    reaches: Vec<usize>,
    /// At an end where the ball changes support, the cap's normal both
    /// pieces share there (the crease's direction): the piece's own
    /// tangent there is the same, but read apart it would leave two caps a
    /// rounding apart.
    joint: [Option<Vector>; 2],
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
            // A radius from the vertex: the sphere about it.
            Support::Around { at, along: d } => {
                let den = self.n1.dot(d);
                if den.abs() <= f64::EPSILON {
                    return None;
                }
                let u = Direction::new(d, Tolerances::default())
                    .ok()?
                    .any_perpendicular()
                    .vector();
                let w = d.cross(u);
                let turn = s / self.radius;
                let out = (u * turn.cos() + w * turn.sin()) * self.radius;
                let level = (self.start + self.n1 * self.radius - at).dot(self.n1);
                let l = (level - self.n1.dot(out)) / den;
                Some(at + d * l + out)
            }
            // A radius from the vertex and a radius over `F1`: the circle
            // about the vertex's foot, `s` its arc length from `along`.
            Support::Point { at } => {
                let lift = (self.start + self.n1 * self.radius - at).dot(self.n1);
                let rho2 = self.radius * self.radius - lift * lift;
                if rho2 <= 0.0 {
                    return None;
                }
                let rho = rho2.sqrt();
                let middle = at + self.n1 * lift;
                let turn = s / rho;
                Some(middle + (self.along * turn.cos() + self.into1 * turn.sin()) * rho)
            }
        }
    }

    /// The ball's touch on `F1` and on its support, for centre `c`.
    fn feet(&self, c: Point) -> (Point, Point) {
        let other = match self.support {
            Support::Edge { at, along } | Support::Around { at, along } => {
                at + along * (c - at).dot(along)
            }
            Support::Plane { normal, .. } => c - normal * self.radius,
            Support::Point { at } => at,
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
    /// The vertex the ball rolls round, touching every face there at it.
    Vertex {
        vertex: Id,
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
) -> OgeomResult<Option<Vec<Run>>> {
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
                Support::Edge { at, .. }
                | Support::Plane { at, .. }
                | Support::Point { at }
                | Support::Around { at, .. } => at,
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
                joint: [None, None],
                reaches: Vec::new(),
            };
            let Some(span) = capped_span(&run, tol)? else {
                return Ok(None);
            };
            let run = Run { span, ..run };
            Ok::<_, ogeom_core::OgeomError>(
                proved(model, curved, &run, held, tol)?.map(|reaches| Run { reaches, ..run }),
            )
        };
        let (mut face, mut near, mut at, mut into) = (f2, edge, mid, into2);
        let mut taken: Vec<usize> = Vec::new();
        while taken.len() <= MOST_STRIPS {
            let Some((far, hit, far_along)) = far_edge(curved, face, near, at, into, tol)? else {
                break;
            };
            // The narrow face widening past the setback along the seam: the
            // ball rests on its far edge up to the station where the strip
            // is as wide as the setback, and on the strip past it (ADR
            // 0001, section 13, question 8: carried across the change).
            if taken.is_empty() {
                let rate = far_along.dot(into2);
                if rate.abs() > ang {
                    let t = (setback - (hit - a).dot(into2)) / rate;
                    let switch = (hit + far_along * t - a).dot(along);
                    if switch > eps && switch < length - eps {
                        // The strip's width at the seam's start: narrower
                        // than the setback there, the ball rests on the edge
                        // first.
                        let den = far_along.dot(along);
                        if den.abs() <= ang {
                            break;
                        }
                        let t0 = (a - hit).dot(along) / den;
                        let narrow_first = (hit + far_along * t0 - a).dot(into2) < setback;
                        return carried(
                            model,
                            curved,
                            (edge, [f1, f2], radius),
                            (a, along, length, n1, into1, n2),
                            (far, hit, far_along),
                            (caps, &end_faces),
                            (switch, !narrow_first),
                            tol,
                        );
                    }
                }
            }
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
                    return Ok(Some(vec![run]));
                }
                // The far edge ending at a vertex short of the seam's end:
                // the ball rests on the edge up to where its touch reaches
                // the vertex, and rolls round the vertex past it.
                if taken.is_empty()
                    && let Some(runs) = round_the_vertex(
                        model,
                        curved,
                        (edge, [f1, face], radius),
                        (a, along, length, n1, into1),
                        (caps, &end_faces),
                        tol,
                    )?
                {
                    return Ok(Some(runs));
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
                return Ok(Some(vec![run]));
            }
            (face, near, at, into) = (g, far, hit, into_g);
        }
    }
    // The ball touching both faces: the plain round, whose fill the tool
    // ends flush on the end faces' planes where the planar engine leaves it
    // (an end standing in the open, past a corner of the face it ends on).
    let Some(into_a) = into_face(model, curved, fa, (mid, along, length), tol)? else {
        return Ok(None);
    };
    let run = Run {
        edge,
        faces: [fa, fb],
        radius,
        start: a,
        along,
        length,
        n1: na,
        into1: into_a,
        support: Support::Plane { at: a, normal: nb },
        taken: Vec::new(),
        caps,
        end_faces,
        span: (0.0, 0.0),
        joint: [None, None],
        reaches: Vec::new(),
    };
    // Along the seam itself both touches stand within the end faces' planes
    // and on their faces: an end face leaning over the seam (a corner's
    // other faces over a narrow face) would cut off a ball that touches
    // neither face there.
    for share in [0.25, 0.5, 0.75] {
        let Some(c) = run.centre(length * share) else {
            return Ok(None);
        };
        let (p, q) = run.feet(c);
        for (face, at) in [(fa, p), (fb, q)] {
            if !run.inside(at, eps)
                || ogeom_algo::classify_on_face(
                    model,
                    &curved.faces[face].shape,
                    at,
                    fine(tol),
                    tol,
                )? == ogeom_algo::Containment::Out
            {
                return Ok(None);
            }
        }
    }
    let Some(span) = capped_span(&run, tol)? else {
        return Ok(None);
    };
    let run = Run { span, ..run };
    Ok(proved(model, curved, &run, Held::Face { face: fb }, tol)?
        .map(|reaches| vec![Run { reaches, ..run }]))
}

/// The two pieces of a run carried across a change of support at station
/// `switch` along the seam: resting on the narrow face's far edge where the
/// face is narrower than the setback (before the switch where `narrowing`
/// is false, after it otherwise), touching the face past it. Each piece is
/// capped at the switch on the plane square to the seam there, which holds
/// both pieces' arcs (the ball touches the face and the edge at once there,
/// its path running along the seam). `None` where either piece is not
/// proved.
#[allow(clippy::too_many_arguments, reason = "the run's whole reading")]
fn carried(
    model: &Model,
    curved: &Curved,
    (edge, [f1, f2], radius): (Id, [usize; 2], f64),
    (a, along, length, n1, into1, n2): (Point, Vector, f64, Vector, Vector, Vector),
    (far, hit, far_along): (Id, Point, Vector),
    (caps, end_faces): ([(Vector, f64); 2], &[usize]),
    (switch, narrowing): (f64, bool),
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Run>>> {
    let ang = tol.angular();
    let [ga, gb] = faces_of(curved, far)[..] else {
        return Ok(None);
    };
    let g = if ga == f2 { gb } else { ga };
    let (Some(nf), Some(ng)) = (outward(curved, f2, hit), outward(curved, g, hit)) else {
        return Ok(None);
    };
    let gc = &curved.edges[&far];
    let (gp, gq) = (
        gc.curve.point_at(gc.range.0, tol)?,
        gc.curve.point_at(gc.range.1, tol)?,
    );
    let Some(into_g) = into_face(
        model,
        curved,
        g,
        (gp.midpoint(gq), far_along, gp.distance(gq)),
        tol,
    )?
    else {
        return Ok(None);
    };
    if into_g.dot(nf) >= -ang {
        return Ok(None);
    }
    let piece = |support: Support| Run {
        edge,
        faces: [f1, f2],
        radius,
        start: a,
        along,
        length,
        n1,
        into1,
        support,
        taken: Vec::new(),
        caps,
        end_faces: end_faces.to_vec(),
        span: (0.0, 0.0),
        joint: [None, None],
        reaches: Vec::new(),
    };
    let on_edge = piece(Support::Edge {
        at: hit,
        along: far_along,
    });
    let on_face = piece(Support::Plane { at: a, normal: n2 });
    let (Some(edge_span), Some(face_span)) =
        (capped_span(&on_edge, tol)?, capped_span(&on_face, tol)?)
    else {
        return Ok(None);
    };
    // The edge's piece where the strip is narrow, the face's past it.
    let (edge_span, face_span, edge_joint, face_joint) = if narrowing {
        (
            (switch, edge_span.1),
            (face_span.0, switch),
            [Some(along), None],
            [None, Some(along)],
        )
    } else {
        (
            (edge_span.0, switch),
            (switch, face_span.1),
            [None, Some(along)],
            [Some(along), None],
        )
    };
    let on_edge = Run {
        span: edge_span,
        joint: edge_joint,
        ..on_edge
    };
    let on_face = Run {
        span: face_span,
        joint: face_joint,
        ..on_face
    };
    let held = Held::Edge {
        far,
        faces: [f2, g],
        normals: [nf, ng],
    };
    let (Some(edge_reaches), Some(face_reaches)) = (
        proved(model, curved, &on_edge, held, tol)?,
        proved(model, curved, &on_face, Held::Face { face: f2 }, tol)?,
    ) else {
        return Ok(None);
    };
    Ok(Some(vec![
        Run {
            reaches: edge_reaches,
            ..on_edge
        },
        Run {
            reaches: face_reaches,
            ..on_face
        },
    ]))
}

/// The two pieces of a run whose far edge ends at a vertex between the
/// seam's ends: resting on the edge where the ball's touch on it stands
/// within it, rolling round the vertex past the station where the touch
/// reaches it (the ball's centre on the sphere a radius about the vertex,
/// its arc ending on the vertex, the pieces' centres' paths tangent there).
/// Each piece is capped at that station on the plane square to their
/// shared path. `None` where neither of the edge's ends is met so, or
/// either piece is not proved.
#[allow(clippy::too_many_arguments, reason = "the run's whole reading")]
fn round_the_vertex(
    model: &Model,
    curved: &Curved,
    (edge, faces, radius): (Id, [usize; 2], f64),
    (a, along, length, n1, into1): (Point, Vector, f64, Vector, Vector),
    (caps, end_faces): ([(Vector, f64); 2], &[usize]),
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Run>>> {
    let eps = tol.confusion();
    let narrow = faces[1];
    let piece = |support: Support| Run {
        edge,
        faces,
        radius,
        start: a,
        along,
        length,
        n1,
        into1,
        support,
        taken: Vec::new(),
        caps,
        end_faces: end_faces.to_vec(),
        span: (0.0, 0.0),
        joint: [None, None],
        reaches: Vec::new(),
    };
    // Each straight edge of the narrow face but the seam, convex, may hold
    // the ball up to one of its ends.
    for (far, _) in curved.faces[narrow].loops.iter().flatten() {
        let far = *far;
        if far == edge {
            continue;
        }
        let g = &curved.edges[&far];
        if !matches!(g.curve, Curve::Line(_)) || g.ends.0 == g.ends.1 {
            continue;
        }
        let (gp, gq) = (curved.points[&g.ends.0], curved.points[&g.ends.1]);
        let far_along = (gq - gp) / gp.distance(gq);
        let [ga, gb] = faces_of(curved, far)[..] else {
            continue;
        };
        let beyond = if ga == narrow { gb } else { ga };
        let (Some(nf), Some(ng)) = (outward(curved, narrow, gp), outward(curved, beyond, gp))
        else {
            continue;
        };
        let Some(into_g) = into_face(
            model,
            curved,
            beyond,
            (gp.midpoint(gq), far_along, gp.distance(gq)),
            tol,
        )?
        else {
            continue;
        };
        // Convex: the face beyond falls behind the narrow face's plane.
        if into_g.dot(nf) >= -tol.angular() {
            continue;
        }
        let held = Held::Edge {
            far,
            faces: [narrow, beyond],
            normals: [nf, ng],
        };
        let on_edge = piece(Support::Around {
            at: gp,
            along: far_along,
        });
        for v in [g.ends.0, g.ends.1] {
            let vertex = curved.points[&v];
            let other = if v == g.ends.0 { g.ends.1 } else { g.ends.0 };
            let toward_other = (curved.points[&other] - vertex).dot(far_along);
            // How far past the vertex, away from the edge, the touch stands.
            let past = |s: f64| -> Option<f64> {
                let c = on_edge.centre(s)?;
                Some(-(on_edge.feet(c).1 - vertex).dot(far_along) * toward_other.signum())
            };
            // Searched once round the edge.
            let (s_lo, s_hi) = (0.0, core::f64::consts::TAU * radius);
            let count = 400_u32;
            let mut brackets = Vec::new();
            let mut prev: Option<(f64, f64)> = None;
            for i in 0..=count {
                let s = s_lo + (s_hi - s_lo) * f64::from(i) / f64::from(count);
                let Some(p) = past(s) else {
                    prev = None;
                    continue;
                };
                if let Some((sp, pp)) = prev
                    && (pp > 0.0) != (p > 0.0)
                {
                    brackets.push((sp, s, pp));
                }
                prev = Some((s, p));
            }
            for (mut lo, mut hi, p0) in brackets {
                for _ in 0..200 {
                    let mid = f64::midpoint(lo, hi);
                    if past(mid).is_some_and(|p| (p > 0.0) == (p0 > 0.0)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let switch = f64::midpoint(lo, hi);
                // The switch's arc plane crosses the seam within it.
                let Some(x) = on_edge.crossing(switch) else {
                    continue;
                };
                let across = (x - a).dot(along);
                if across <= eps || across >= length - eps {
                    continue;
                }
                let round = piece(Support::Point { at: vertex });
                // The ball there rests on the edge's end: its centre a radius
                // from the vertex, where the vertex's piece starts.
                let Some(cv) = on_edge.centre(switch) else {
                    continue;
                };
                if (cv.distance(vertex) - radius).abs() > eps * 1e3 {
                    continue;
                }
                // The edge's piece runs from the switch the way its touch stays
                // within the edge; the vertex's piece the other way, its arc
                // length read round its circle from the switch's centre.
                let dir_e = if p0 > 0.0 { 1.0 } else { -1.0 };
                let turn_at = {
                    let lift = (a + n1 * radius - vertex).dot(n1);
                    let middle = vertex + n1 * lift;
                    let rho = (radius * radius - lift * lift).max(0.0).sqrt();
                    let d = cv - middle;
                    rho * d.dot(into1).atan2(d.dot(along))
                };
                let (Some(tv), Some(te)) = (round.tangent(turn_at), on_edge.tangent(switch)) else {
                    continue;
                };
                let joint = te / te.magnitude();
                let Some(edge_cap) = cap_toward(&on_edge, switch, dir_e, None, tol)? else {
                    continue;
                };
                // The vertex's piece runs on to the other seam end: the way its
                // arcs' planes cross the seam toward it.
                let across = |run: &Run, s: f64| -> Option<f64> {
                    Some((run.crossing(s)? - run.start).dot(run.along))
                };
                let (Some(x_switch), Some(x_cap)) =
                    (across(&on_edge, switch), across(&on_edge, edge_cap))
                else {
                    continue;
                };
                let base_v = if x_cap > x_switch { 0.0 } else { length };
                let h = radius * 1e-4;
                let (Some(xp), Some(xm)) =
                    (across(&round, turn_at + h), across(&round, turn_at - h))
                else {
                    continue;
                };
                let dir_v = if (xp - xm) * (base_v - x_switch) > 0.0 {
                    1.0
                } else {
                    -1.0
                };
                let _ = tv;
                let Some(round_cap) = cap_toward(&round, turn_at, dir_v, Some(base_v), tol)? else {
                    continue;
                };
                let span = |from: f64, to: f64| if from <= to { (from, to) } else { (to, from) };
                // The joint at whichever end of each span the switch is.
                let joint_at = |from_switch_low: bool| {
                    if from_switch_low {
                        [Some(joint), None]
                    } else {
                        [None, Some(joint)]
                    }
                };
                let (edge_span, edge_joint) = (span(switch, edge_cap), joint_at(dir_e > 0.0));
                let (round_span, round_joint) = (span(turn_at, round_cap), joint_at(dir_v > 0.0));
                let resting = Run {
                    span: edge_span,
                    joint: edge_joint,
                    ..piece(on_edge.support)
                };
                let rolling = Run {
                    span: round_span,
                    joint: round_joint,
                    ..piece(round.support)
                };
                let (Some(edge_reaches), Some(round_reaches)) = (
                    proved(model, curved, &resting, held, tol)?,
                    proved(model, curved, &rolling, Held::Vertex { vertex: v }, tol)?,
                ) else {
                    continue;
                };
                return Ok(Some(vec![
                    Run {
                        reaches: edge_reaches,
                        ..resting
                    },
                    Run {
                        reaches: round_reaches,
                        ..rolling
                    },
                ]));
            }
        }
    }
    Ok(None)
}

/// The fill tool of a run rolling round a vertex or the edge out of it:
/// a band between the caps `(x0, d0)` and `(x1, d1)` (each a point of the
/// crease and the cap's normal), the arc about the vertex turned through
/// the run or the arcs to the edge fitted; closed by its flap on the
/// touched face's plane, the flap's triangle on the far face's plane from
/// the crease to the support, and the two caps.
fn vertex_tool(
    model: &mut Model,
    curved: &Curved,
    run: &Run,
    (x0, d0): (Point, Vector),
    (x1, d1): (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<Option<Tool>> {
    use ogeom_geom::{BSpline2d, CircleCurve, PlaneSurface};
    use ogeom_math::{Circle, KnotVector, Plane, Point2};
    let round_vertex = match run.support {
        Support::Point { at } => Some(at),
        Support::Around { .. } => None,
        _ => return Ok(None),
    };
    let eps = tol.confusion();
    let r = run.radius;
    let (lo, hi) = run.span;
    // Each cap's arc lies in the cap's plane: the stations there.
    // Rolling round a vertex, the balls' arcs all lie on the spindle torus
    // their centres' circle sweeps about the vertex's axis: the band is that
    // torus, turned from the first arc to the last. Resting on an edge, it
    // is fitted to the arcs.
    let mut n = FIRST;
    let fitted = if round_vertex.is_some() {
        None
    } else {
        Some(loop {
            let mut rows = Vec::with_capacity(n);
            let mut feet = Vec::with_capacity(n);
            for i in 0..n {
                #[expect(clippy::cast_precision_loss, reason = "a few stations")]
                let s = lo + (hi - lo) * (i as f64) / ((n - 1) as f64);
                let (Some(c), Some(mut arc)) = (run.centre(s), run.arc(s, tol)?) else {
                    return Ok(None);
                };
                if let (Some(last), Some(vertex)) = (arc.last_mut(), round_vertex) {
                    *last = vertex;
                }
                feet.push(run.feet(c).0);
                rows.push(arc);
            }
            let surface = ogeom_geom::fit::fit_surface_grid_chordal(&rows, 3, FIT, tol)?;
            let contact = ogeom_geom::fit::fit_points(&feet, 3, FIT, tol)?;
            if surface.met && contact.met {
                break (surface.curve, contact.curve);
            }
            if n >= MOST {
                return Ok(None);
            }
            n = 2 * n - 1;
        })
    };
    let (Some(c0), Some(c1)) = (run.centre(lo), run.centre(hi)) else {
        return Ok(None);
    };
    let (f0, f1) = (run.feet(c0).0, run.feet(c1).0);
    // The support's side: the vertex, or the edge's touches at the caps.
    let (q0, q1) = (run.feet(c0).1, run.feet(c1).1);
    let vertex_shape = ogeom_algo::make_vertex(model, q0).shape;
    let far_shape = match round_vertex {
        Some(_) => vertex_shape.clone(),
        None => ogeom_algo::make_vertex(model, q1).shape,
    };
    let [v_f0, v_f1, v_x0, v_x1] =
        [f0, f1, x0, x1].map(|p| ogeom_algo::make_vertex(model, p).shape);
    let line = |model: &mut Model, (p, vp): (Point, &Shape), (q, vq): (Point, &Shape)| {
        let curve: Curve = LineCurve::segment(p, q, tol)?.into();
        let range = curve.domain();
        Ok::<_, ogeom_core::OgeomError>(
            ogeom_algo::make_edge_between(model, curve, range, vp, vq, tol)?.shape,
        )
    };
    // Each cap's arc, from the foot to the support's touch about the ball's
    // centre.
    let arc_edge =
        |model: &mut Model, c: Point, (f, vf): (Point, &Shape), (q, vq): (Point, &Shape)| {
            let x = Direction::new(f - c, tol)?;
            let z = Direction::new((f - c).cross(q - c), tol)?;
            let turn = (f - c).cross(q - c).magnitude().atan2((f - c).dot(q - c));
            let circle = Circle::new(Frame::new(c, z, x, tol)?, r, tol)?;
            let curve = Curve::Circle(CircleCurve::new(circle));
            Ok::<_, ogeom_core::OgeomError>(
                ogeom_algo::make_edge_between(model, curve, (0.0, turn), vf, vq, tol)?.shape,
            )
        };
    let a0 = arc_edge(model, c0, (f0, &v_f0), (q0, &vertex_shape))?;
    let (a0, e_contact, a1, revolved) = match (&fitted, round_vertex) {
        (Some((_, contact)), _) => {
            let contact_curve = Curve::BSpline(contact.clone());
            let contact_range = contact_curve.domain();
            let e_contact = ogeom_algo::make_edge_between(
                model,
                contact_curve,
                contact_range,
                &v_f0,
                &v_f1,
                tol,
            )?
            .shape;
            let a1 = arc_edge(model, c1, (f1, &v_f1), (q1, &far_shape))?;
            (a0, e_contact, a1, None)
        }
        (None, Some(vertex)) => {
            let Some(turned) = revolved_band(model, run, vertex, &a0, (f0, f1), tol)? else {
                return Ok(None);
            };
            (
                turned.a0.clone(),
                turned.contact.clone(),
                turned.a1.clone(),
                Some(turned),
            )
        }
        (None, None) => return Ok(None),
    };
    // The turned band's own corners, where it is one.
    let (v_f0, v_f1, vertex_shape, far_shape) = match &revolved {
        Some(turned) => (
            turned.f0.clone(),
            turned.f1.clone(),
            turned.vertex.clone(),
            turned.vertex.clone(),
        ),
        None => (v_f0, v_f1, vertex_shape, far_shape),
    };
    let foot0 = line(model, (f0, &v_f0), (x0, &v_x0))?;
    let foot1 = line(model, (f1, &v_f1), (x1, &v_x1))?;
    let crease = line(model, (x0, &v_x0), (x1, &v_x1))?;
    let rise0 = line(model, (x0, &v_x0), (q0, &vertex_shape))?;
    let rise1 = line(model, (x1, &v_x1), (q1, &far_shape))?;
    // The support's contact along the edge, where the ball rests on one.
    let support = match round_vertex {
        Some(_) => None,
        None => Some(line(model, (q0, &vertex_shape), (q1, &far_shape))?),
    };
    let band = if let Some(turned) = revolved {
        turned.band
    } else {
        let Some((surface, _)) = fitted else {
            return Ok(None);
        };
        // The band: `u` across each arc from the foot to the vertex, `v` along
        // the stations; its side `u = 1` the vertex.
        let id = model
            .geometry_mut()
            .add_surface(SurfaceGeometry::BSpline(surface));
        let pcurve =
            |model: &mut Model, edge: &Shape, from: Point2, to: Point2, range: (f64, f64)| {
                let line = BSpline2d::new(
                    KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?,
                    vec![from, to],
                    tol,
                )?;
                ogeom_algo::attach_pcurve(
                    model,
                    edge,
                    line.into(),
                    id,
                    ogeom_topo::Location::identity(),
                    range,
                )
            };
        let range_of = |model: &Model, e: &Shape| -> OgeomResult<(f64, f64)> {
            let Some((_, range)) = ogeom_algo::placed_edge_curve(model, e, tol)? else {
                ogeom_core::ogeom_bail!(Construction, "an edge of a vertex fill has no curve");
            };
            Ok(range)
        };
        let (o, x, y, xy) = (
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 1.0),
        );
        let rc = range_of(model, &e_contact)?;
        pcurve(model, &e_contact, o, y, rc)?;
        let r1 = range_of(model, &a1)?;
        pcurve(model, &a1, y, xy, r1)?;
        let r0 = range_of(model, &a0)?;
        pcurve(model, &a0, o, x, r0)?;
        let far_side = match &support {
            Some(edge) => {
                let range = range_of(model, edge)?;
                pcurve(model, edge, x, xy, range)?;
                edge.reversed()
            }
            None => {
                let collapsed = {
                    let mut data = ogeom_topo::EdgeData::new_in(tol);
                    data.degenerate = true;
                    model.add_edge(data, &[vertex_shape.clone(), vertex_shape.clone()])?
                };
                pcurve(model, &collapsed, xy, x, (0.0, 1.0))?;
                collapsed
            }
        };
        let fit_gap = ogeom_core::Tolerance::new_in(FIT.max(eps) * 10.0, tol)?;
        for e in [&e_contact, &a0, &a1].into_iter().chain(support.as_ref()) {
            model.widen(e, fit_gap)?;
        }
        // The ring winds about the band's own normal, as the planes' rings
        // wind about theirs: listed the other way round where it does not.
        let mut edges = vec![e_contact.clone(), a1.clone(), far_side, a0.reversed()];
        let ring = ogeom_algo::make_wire(model, &edges, tol)?.shape;
        let mut band = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&ring), tol)?.shape;
        let (_, n) = ogeom_algo::face_normal(model, &band, tol)?;
        if wound(model, &edges, n, tol)?.dot(n) < 0.0 {
            edges.reverse();
            let edges: Vec<Shape> = edges.iter().map(Shape::reversed).collect();
            let ring = ogeom_algo::make_wire(model, &edges, tol)?.shape;
            band = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&ring), tol)?.shape;
        }
        band
    };
    // The fill lies between the band and the crease, the balls outside it:
    // its band faces the balls' centres.
    let band = {
        let mid = 0.5 * (lo + hi);
        let Some(c) = run.centre(mid) else {
            return Ok(None);
        };
        let (p, n) = ogeom_algo::face_normal(model, &band, tol)?;
        if n.dot(c - p) < 0.0 {
            band.reversed()
        } else {
            band
        }
    };
    let plane_through = |p: Point, normal: Vector| -> OgeomResult<SurfaceGeometry> {
        Ok(PlaneSurface::new(Plane::through(p, Direction::new(normal, tol)?)).into())
    };
    let [f_1, f_2] = run.faces;
    let (Some(n_1), Some(n_2)) = (outward(curved, f_1, run.start), outward(curved, f_2, q0)) else {
        return Ok(None);
    };
    // Every planar face facing out of the fill: into the faces the fill
    // stands on, and away from the run past its caps.
    let along_run = c1 - c0;
    let out0 = if d0.dot(along_run) > 0.0 { -d0 } else { d0 };
    let out1 = if d1.dot(along_run) > 0.0 { d1 } else { -d1 };
    // A ring is walked as it is listed, so each plane faces the way its
    // ring winds about it; which faces turn is settled after.
    let planar = |model: &mut Model, at: Point, normal: Vector, ring: Vec<Shape>| {
        let normal = wound(model, &ring, normal, tol)?;
        crate::support::face_from_edges(model, plane_through(at, normal)?, &ring, tol)
    };
    let flap = planar(
        model,
        x0,
        -n_1,
        vec![
            e_contact.reversed(),
            foot0.clone(),
            crease.clone(),
            foot1.reversed(),
        ],
    )?;
    let mut on_far: Vec<Shape> = vec![crease.reversed(), rise0.clone()];
    if let Some(edge) = &support {
        on_far.push(edge.clone());
    }
    on_far.push(rise1.reversed());
    let triangle = planar(model, x0, -n_2, on_far)?;
    let cap0 = planar(
        model,
        x0,
        out0,
        vec![a0.clone(), rise0.reversed(), foot0.reversed()],
    )?;
    let cap1 = planar(
        model,
        x1,
        out1,
        vec![a1.clone(), rise1.reversed(), foot1.reversed()],
    )?;
    // The faces share their edges already, so the topology says which of
    // them to turn for every edge to be walked once each way; the band
    // keeps its way, and the whole shell then faces out by its material.
    let mut faces = vec![band.clone(), flap, triangle, cap0, cap1];
    let Some(agreement) = ogeom_algo::agreeing_orientation(model, &faces)? else {
        return Ok(None);
    };
    let band_turned = agreement.turn[0];
    for (face, turn) in faces.iter_mut().zip(&agreement.turn) {
        if *turn != band_turned {
            *face = face.reversed();
        }
    }
    let band = faces[0].clone();
    // Built on shared edges the faces close as they are; sewn where the
    // fitted band's edges only meet theirs within its fit.
    let made = ogeom_algo::make_shell(model, &faces)?.shape;
    let (shell, band_now) = if ogeom_algo::is_shell_closed(model, &made)? {
        (made, vec![band])
    } else {
        let sewn = ogeom_algo::sew(model, &faces, tol)?;
        let [shell] = sewn.shells.as_slice() else {
            return Ok(None);
        };
        if !ogeom_algo::is_shell_closed(model, shell)? {
            return Ok(None);
        }
        (shell.clone(), sewn.history.trace(&band).to_vec())
    };
    let solid = ogeom_algo::make_solid(model, std::slice::from_ref(&shell))?.shape;
    let tool = Tool {
        solid,
        sense: Sense::Fuse,
        faces: band_now
            .into_iter()
            .map(|f| (f, curved.shapes[&run.edge].clone()))
            .collect(),
    };
    // Its caps are the arcs' planes through the seam's ends (or the switch),
    // where the sections end: the fill is those sections whole, not cut
    // flush to the end faces, which here meet at the vertex the balls
    // roll round and pass under them.
    Ok(Some(tool))
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
/// keep them); resting on a far edge, its centre within the edge's normal
/// cone (on its side only where it comes down on a face); and no face of the
/// solid but `F1` and the end faces (the fill is cut to their planes;
/// round a vertex, where it is not, the end faces too) nearer than the
/// radius, every face but those it rests against and the
/// planar faces wholly behind a plane it touches by the distance's
/// Lipschitz bound up to the next station, halving where it is close
/// (ADR 0001, section 13, question 8); and the arcs' planes cross the
/// crease in order. `Some` with the faces the balls reach (none of those
/// the run rests on or takes in), `None` where it does not hold.
fn proved(
    model: &Model,
    curved: &Curved,
    run: &Run,
    held: Held,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<usize>>> {
    let eps = tol.confusion();
    let margin = run.radius * 1e-6;
    // A run round a vertex or the edge out of it is not cut to the end
    // faces' planes: all of it is asked, the end faces with the rest.
    let uncut = matches!(run.support, Support::Point { .. } | Support::Around { .. });
    let kept = |p: Point| uncut || run.inside(p, eps);
    // Within the far edge, and resting on it alone.
    let on_edge = |c: Point, on_far: Point| -> OgeomResult<bool> {
        let Held::Edge { far, normals, .. } = held else {
            return Ok(true);
        };
        let (Support::Edge { at, along } | Support::Around { at, along }) = run.support else {
            return Ok(false);
        };
        let g = &curved.edges[&far];
        let (p, q) = (
            g.curve.point_at(g.range.0, tol)?,
            g.curve.point_at(g.range.1, tol)?,
        );
        let (x, y) = ((p - at).dot(along), (q - at).dot(along));
        let t = (on_far - at).dot(along);
        // Resting on the edge: from the touch to the centre within the
        // edge's normal cone, `w = u n0 + v n1` with `u` and `v` not below
        // nought (on its side, where the ball comes down on a face there).
        let w = c - on_far;
        let g01 = normals[0].dot(normals[1]);
        let det = 1.0 - g01 * g01;
        if det <= f64::EPSILON {
            return Ok(false);
        }
        let (w0, w1) = (w.dot(normals[0]), w.dot(normals[1]));
        let (u, v) = ((w0 - g01 * w1) / det, (w1 - g01 * w0) / det);
        // Within the edge, its ends included: an end ball may touch it at
        // its end, on an end face's plane, where the fill is cut.
        Ok(
            (!kept(on_far) || (t >= x.min(y) - eps && t <= x.max(y) + eps))
                && u >= -margin
                && v >= -margin,
        )
    };
    // A quick look at three stations of the run along the seam before the
    // faces are read.
    // A run read round a line or a vertex has its own parameter, not the
    // seam's stations.
    let (from, to) = match run.support {
        Support::Around { .. } | Support::Point { .. } => run.span,
        _ => (run.span.0.max(0.0), run.span.1.min(run.length)),
    };
    for share in [0.0, 0.5, 1.0] {
        let Some(c) = run.centre(from + (to - from) * share) else {
            return Ok(None);
        };
        if !on_edge(c, run.feet(c).1)? {
            return Ok(None);
        }
    }
    let resting: Vec<usize> = match held {
        Held::Edge { faces, .. } => faces.to_vec(),
        Held::Face { face } => vec![face],
        Held::Vertex { vertex } => (0..curved.faces.len())
            .filter(|fi| {
                curved.faces[*fi].loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == vertex || g.ends.1 == vertex
                })
            })
            .collect(),
    };
    let others: Vec<(usize, PointDistance)> = (0..curved.faces.len())
        .filter(|fi| *fi != run.faces[0] && (uncut || !run.end_faces.contains(fi)))
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
    let mut reaches: Vec<usize> = Vec::new();
    // The faces the run holds, and the edges of their outlines (and of the
    // touched face's).
    let holds: Vec<usize> = resting
        .iter()
        .chain(&run.taken)
        .chain(core::iter::once(&run.faces[1]))
        .copied()
        .collect();
    let mut outlines: Vec<Id> = holds
        .iter()
        .flat_map(|f| curved.faces[*f].loops.iter().flatten().map(|(e, _)| *e))
        .collect();
    outlines.extend(
        curved.faces[run.faces[0]]
            .loops
            .first()
            .into_iter()
            .flatten()
            .map(|(e, _)| *e),
    );
    let apart = |f: usize| {
        !holds.contains(&f)
            && !curved.faces[f]
                .loops
                .iter()
                .flatten()
                .any(|(e, _)| outlines.contains(e))
    };
    let mut last: Option<f64> = None;
    let mut order: Option<f64> = None;
    for i in 0..=PROOF {
        #[expect(clippy::cast_precision_loss, reason = "a few stations")]
        let s = lo + step * (i as f64);
        let (Some(c), Some(arc), Some(x)) = (run.centre(s), run.arc(s, tol)?, run.crossing(s))
        else {
            return Ok(None);
        };
        // The arcs' planes cross the crease's line in order, one way or
        // the other along the run.
        let at = (x - run.start).dot(run.along);
        if let Some(l) = last {
            let step = at - l;
            if step == 0.0 || order.is_some_and(|o: f64| o * step <= 0.0) {
                return Ok(None);
            }
            order = Some(step.signum());
        }
        last = Some(at);
        // Only what the end faces' planes keep is built.
        if !arc.iter().any(|p| kept(*p)) {
            continue;
        }
        let (foot, other) = run.feet(c);
        if !on_edge(c, other)? {
            return Ok(None);
        }
        let on = |face: usize, p: Point| -> OgeomResult<bool> {
            Ok(!kept(p)
                || ogeom_algo::classify_on_face(
                    model,
                    &curved.faces[face].shape,
                    p,
                    fine(tol),
                    tol,
                )? != ogeom_algo::Containment::Out)
        };
        if !on(run.faces[0], foot)? {
            return Ok(None);
        }
        if let Held::Face { face } = held
            && !on(face, other)?
        {
            return Ok(None);
        }
        // Every face but those the ball rests against clear by the radius
        // at every centre between this station and the next: the distance
        // moves by no more than the centre does.
        for (fi, faces) in &others {
            if reaches.contains(fi) {
                continue;
            }
            let d = faces.nearest(c, tol)?.0;
            let clear = d >= run.radius - eps
                && (resting.contains(fi) || behind[*fi] || i == PROOF || {
                    let mut budget = MEASURES;
                    clear_between(run, faces, (s, s + step), d, (HALVINGS, &mut budget), tol)?
                });
            if !clear {
                // A face the run rests on, takes in or closes on must stand
                // clear, and so must one joined to those or to the touched
                // face's outline (a strip further on, a wall): the run is not
                // this one. Any other (a post on the touched face, joined
                // to it only through a hole) is one the round reaches.
                if !apart(*fi) {
                    return Ok(None);
                }
                reaches.push(*fi);
            }
        }
    }
    Ok(Some(reaches))
}

/// Whether every centre of `run` between stations `s0` and `s1` stands a
/// radius or more from `faces`, `d0` its distance at `s0` (the one at `s1`
/// is checked at that station): the distance moves by no more than the
/// centre does, so a margin of half the chord between them clears the
/// interval; where it does not, the interval is halved, at most `depth`
/// times and while `budget` distances are left to measure.
fn clear_between(
    run: &Run,
    faces: &PointDistance,
    (s0, s1): (f64, f64),
    d0: f64,
    (depth, budget): (u32, &mut u32),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let eps = tol.confusion();
    if *budget == 0 {
        return Ok(false);
    }
    *budget -= 1;
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
    Ok(
        clear_between(run, faces, (s0, m), d0, (depth - 1, budget), tol)?
            && clear_between(run, faces, (m, s1), dm, (depth - 1, budget), tol)?,
    )
}

/// The station that caps `run` stepping from `from` the way `dir` says:
/// the first, half a radius at a time, whose arc and seam crossing stand
/// wholly outside either end face's plane; else, and for a run read round
/// a line or a vertex, where its arc's plane crosses a seam end (the balls
/// run out round a vertex, or a short seam's arcs never clear its end
/// face).
/// Where `only` names a seam end, only that one caps. `None` where neither
/// is met within twenty radii.
fn cap_toward(
    run: &Run,
    from: f64,
    dir: f64,
    only: Option<f64>,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let eps = tol.confusion();
    let r = run.radius;
    let outside = |s: f64| -> OgeomResult<bool> {
        let (Some(arc), Some(x)) = (run.arc(s, tol)?, run.crossing(s)) else {
            return Ok(false);
        };
        Ok(run.caps.iter().any(|(n, d)| {
            arc.iter()
                .chain(core::iter::once(&x))
                .all(|p| n.dot(*p - Point::ORIGIN) > d + eps)
        }))
    };
    if !matches!(run.support, Support::Point { .. } | Support::Around { .. }) {
        for step in 1..=40 {
            let s = (r * 0.5 * f64::from(step)).mul_add(dir, from);
            if outside(s)? {
                return Ok(Some(s));
            }
        }
    }
    let across = |s: f64| -> Option<f64> { Some((run.crossing(s)? - run.start).dot(run.along)) };
    let Some(x0) = across(from) else {
        return Ok(None);
    };
    let count = 800_u32;
    let mut prev = (from, x0);
    for i in 1..=count {
        let s = (20.0 * r * f64::from(i) / f64::from(count)).mul_add(dir, from);
        let Some(x) = across(s) else {
            return Ok(None);
        };
        for base in [0.0, run.length] {
            if only.is_some_and(|o| o != base) {
                continue;
            }
            if (prev.1 - base > 0.0) != (x - base > 0.0) {
                let (mut lo, mut hi) = (prev.0, s);
                let low = prev.1 - base > 0.0;
                for _ in 0..200 {
                    let mid = f64::midpoint(lo, hi);
                    if across(mid).is_some_and(|m| (m - base > 0.0) == low) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return Ok(Some(f64::midpoint(lo, hi)));
            }
        }
        prev = (s, x);
    }
    Ok(None)
}

/// A band turned about a vertex, with its edges and corners.
struct Turned {
    band: Shape,
    /// The first arc, the last, and the contact circle on the touched face.
    a0: Shape,
    a1: Shape,
    contact: Shape,
    /// The arcs' feet and the vertex they all end on.
    f0: Shape,
    f1: Shape,
    vertex: Shape,
}

/// The band of `run` rolling round `vertex`: the arc `arc` from the
/// touched face to the vertex, turned about the line through the vertex
/// square to that face through the run's span of its centres' circle, a
/// patch of the spindle torus the balls sweep. `None` where it does not
/// come round from the foot `feet.0` to `feet.1`, or its edges are not the
/// two arcs and the contact circle.
fn revolved_band(
    model: &mut Model,
    run: &Run,
    vertex: Point,
    arc: &Shape,
    feet: (Point, Point),
    tol: Tolerances,
) -> OgeomResult<Option<Turned>> {
    use ogeom_topo::ShapeType;
    let near = tol.confusion();
    let (lo, hi) = run.span;
    let lift = (run.start + run.n1 * run.radius - vertex).dot(run.n1);
    let rho2 = run.radius * run.radius - lift * lift;
    if rho2 <= 0.0 {
        return Ok(None);
    }
    let turn = (hi - lo) / rho2.sqrt();
    if turn <= tol.angular() || turn >= std::f64::consts::TAU {
        return Ok(None);
    }
    // The centres turn from `along` toward `into1`.
    let axis = ogeom_math::Axis {
        location: vertex,
        direction: Direction::new(run.along.cross(run.into1), tol)?,
    };
    let swept = ogeom_algo::make_revolution(model, arc, axis, turn, tol)?.shape;
    let [band] = &ogeom_topo::explore_unique(model, &swept, ShapeType::Face)?[..] else {
        return Ok(None);
    };
    let band = band.clone();
    let (mut a0, mut a1, mut contact) = (None, None, None);
    let (mut f0, mut f1, mut at_vertex) = (None, None, None);
    for edge in ogeom_topo::explore_unique(model, &band, ShapeType::Edge)? {
        if model
            .node(&edge)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.degenerate)
        {
            continue;
        }
        let Some((s, t)) = ogeom_algo::edge_vertices(model, &edge)? else {
            return Ok(None);
        };
        let (ps, pt) = (
            crate::setback::point_of(model, &s)?,
            crate::setback::point_of(model, &t)?,
        );
        // Each arc walked from its foot to the vertex, the contact from
        // the first foot to the last.
        let (tip, foot, p, edge) = if ps.distance(vertex) <= near {
            (Some(s), t, pt, edge.reversed())
        } else if pt.distance(vertex) <= near {
            (Some(t), s, ps, edge)
        } else if ps.distance(feet.0) <= near {
            (None, s, ps, edge)
        } else {
            (None, t, pt, edge.reversed())
        };
        match tip {
            Some(tip) if p.distance(feet.0) <= near => {
                (a0, f0, at_vertex) = (Some(edge), Some(foot), Some(tip));
            }
            Some(_) if p.distance(feet.1) <= near => (a1, f1) = (Some(edge), Some(foot)),
            None if [ps, pt]
                .iter()
                .all(|q| q.distance(feet.0) <= near || q.distance(feet.1) <= near) =>
            {
                contact = Some(edge);
            }
            _ => return Ok(None),
        }
    }
    let (Some(a0), Some(a1), Some(contact), Some(f0), Some(f1), Some(vertex)) =
        (a0, a1, contact, f0, f1, at_vertex)
    else {
        return Ok(None);
    };
    Ok(Some(Turned {
        band,
        a0,
        a1,
        contact,
        f0,
        f1,
        vertex,
    }))
}

/// `normal` or its opposite, whichever `ring` winds about by the right
/// hand: Newell's sum over points along its edges as they are walked.
fn wound(model: &Model, ring: &[Shape], normal: Vector, tol: Tolerances) -> OgeomResult<Vector> {
    const STEPS: usize = 8;
    let mut points = Vec::with_capacity(ring.len() * STEPS);
    for edge in ring {
        let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            continue;
        };
        let forward = edge.orientation() == ogeom_topo::Orientation::Forward;
        for k in 0..STEPS {
            #[expect(clippy::cast_precision_loss, reason = "a few points an edge")]
            let share = k as f64 / STEPS as f64;
            let t = if forward {
                a + (b - a) * share
            } else {
                b - (b - a) * share
            };
            points.push(curve.point_at(t, tol)?);
        }
    }
    let mut sum = Vector::ZERO;
    for (k, &p) in points.iter().enumerate() {
        let q = points[(k + 1) % points.len()];
        sum += (p - Point::ORIGIN).cross(q - Point::ORIGIN);
    }
    Ok(if sum.dot(normal) < 0.0 {
        -normal
    } else {
        normal
    })
}

/// The stations whose arcs' planes cap the tool of `run`: at each end, the
/// first, stepping out from the seam's end half a radius at a time, whose
/// arc and seam crossing stand wholly outside that end face's plane.
/// `None` where none does within twenty radii.
fn capped_span(run: &Run, tol: Tolerances) -> OgeomResult<Option<(f64, f64)>> {
    let (Some(s0), Some(s1)) = (capped_end(run, 0, tol)?, capped_end(run, 1, tol)?) else {
        return Ok(None);
    };
    Ok(Some((s0, s1)))
}

/// The station whose arc's plane caps the tool of `run` at its end `k`
/// (0 at the seam's start, 1 at its end): the first, stepping out from the
/// seam's end half a radius at a time, whose arc and seam crossing stand
/// wholly outside that end face's plane. `None` where none does within
/// twenty radii.
fn capped_end(run: &Run, k: usize, tol: Tolerances) -> OgeomResult<Option<f64>> {
    let eps = tol.confusion();
    let r = run.radius;
    let (n, d) = run.caps[k];
    let outside = |s: f64| -> OgeomResult<bool> {
        let (Some(arc), Some(x)) = (run.arc(s, tol)?, run.crossing(s)) else {
            return Ok(false);
        };
        Ok(arc
            .iter()
            .chain(core::iter::once(&x))
            .all(|p| n.dot(*p - Point::ORIGIN) > d + eps))
    };
    let (base, sign) = if k == 0 {
        (0.0, -1.0)
    } else {
        (run.length, 1.0)
    };
    let mut ran_out = false;
    for step in 1..=40 {
        let s = base + sign * r * 0.5 * f64::from(step);
        if outside(s)? {
            return Ok(Some(s));
        }
        // The line of centres at a station draws away from the edge's line
        // linearly: once the ball cannot reach the edge, no later one can.
        if matches!(run.support, Support::Edge { .. }) && run.centre(s).is_none() {
            ran_out = true;
            break;
        }
    }
    // Rolling round a vertex, or resting on a far edge that leaves `F1`
    // past the seam's end, the balls run out before their arcs clear the
    // end face's plane: the run ends where its arc's plane crosses the
    // seam's end (ADR 0001, amendment "a far edge leaving the seam").
    if matches!(run.support, Support::Point { .. })
        || (ran_out && matches!(run.support, Support::Edge { .. }))
    {
        let across =
            |s: f64| -> Option<f64> { Some((run.crossing(s)? - run.start).dot(run.along) - base) };
        let count = 400_u32;
        let reach = 20.0 * r;
        let mut prev: Option<(f64, f64)> = None;
        for i in 0..=count {
            let s = (2.0 * reach * f64::from(i) / f64::from(count) - reach).mul_add(sign, base);
            let Some(x) = across(s) else {
                prev = None;
                continue;
            };
            if let Some((sp, xp)) = prev
                && (xp > 0.0) != (x > 0.0)
            {
                let (mut lo, mut hi) = (sp, s);
                for _ in 0..200 {
                    let mid = f64::midpoint(lo, hi);
                    if across(mid).is_some_and(|m| (m > 0.0) == (xp > 0.0)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return Ok(Some(f64::midpoint(lo, hi)));
            }
            prev = Some((s, x));
        }
    }
    Ok(None)
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
    // Each cap on the arc's plane there, or on the plane both pieces share
    // at a change of support.
    let cap = |s: f64, joint: Option<Vector>| -> Option<(Point, Vector)> {
        let normal = match joint {
            Some(n) => n,
            None => run.tangent(s)?,
        };
        let c = run.centre(s)?;
        let den = normal.dot(run.along);
        (den.abs() > f64::EPSILON * normal.magnitude()).then(|| {
            (
                run.start + run.along * ((c - run.start).dot(normal) / den),
                normal,
            )
        })
    };
    let (Some((x0, d0)), Some((x1, d1))) = (cap(lo, run.joint[0]), cap(hi, run.joint[1])) else {
        return Ok(None);
    };
    // Read round a line or a vertex: built face by face from its stations.
    if matches!(run.support, Support::Point { .. } | Support::Around { .. }) {
        return vertex_tool(model, curved, run, (x0, d0), (x1, d1), tol);
    }
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
        Support::Point { .. } | Support::Around { .. } => false,
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
