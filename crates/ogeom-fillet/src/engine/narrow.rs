//! Narrow passages (ADR 0001, amendment "overflow handling"): concave
//! rounds of one floor's rim where the zone between them is tighter than
//! the radius, so the ball cannot roll through touching its seats on both
//! sides: two walls no farther apart than twice the radius, or a corner
//! face whose curvature is tighter than the radius.
//!
//! Every round keeps its own constant section. Each wall's band runs on
//! past its edge's ends through the tight zone, to the plane of the next
//! wall or to the far side of the face closing the zone, and the rounds
//! meet each other where they cross. The fill is the union of the bands'
//! fills, cut to the solid's air: one spandrel prism per wall, each fused
//! into the solid. The faces of the tight zone (a corner face, or a floor
//! narrower than twice the radius) are buried where the fill covers them.
//!
//! Built for a planar floor with planar walls square to it. Neighbouring
//! walls meet at a sharp corner of any angle (a concave one: the band runs
//! on to the next wall's plane and is cut there; a convex one: the ball
//! pivots on the walls' edge and the spandrel is turned about it), at a
//! square corner through a corner face tighter than the radius all along
//! (a circular or elliptic cylinder square to the floor), or across such a
//! face closing two facing walls, tangent to both; open ends are
//! flush on end faces square to their edge. Facing walls may stand closer
//! than the radius, each band then crossing the wall it faces. Pieces are
//! kept only where nothing else of the solid stands in their room, a
//! section entering the solid through a wall never comes out into air
//! again, the fuses add no more than the pieces hold, and the result is
//! valid.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Aabb, Axis, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{Curved, read};
use super::overflow::{arc, segment};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// Which narrow-passage strategy a passage needs at the least. Each
/// strategy takes what the ones before it take.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Needs {
    /// The tight zone is made of constant-radius rounds only: the
    /// request's own rounds crossing over a narrow floor, and corner faces
    /// that are cylinders tangent to both their walls.
    BlendMade,
    /// The tight zone spans one whole corner face that no constant-radius
    /// round made: a cylinder meeting its walls at kinks, or an elliptic
    /// cylinder.
    WholeFace,
    /// The tight zone covers part of a wall's face: walls facing each
    /// other closer than twice the radius over part of one's length only.
    PartialFace,
}

/// One wall's band, run on through the tight zone at either end: its
/// prism starts at `start` and runs `len` along `along`, `inward` the
/// wall's outward normal (into the passage). `lead` is how far it runs on
/// past its edge at the start and at the end.
#[derive(Debug, Clone)]
struct Band {
    edge: Id,
    start: Point,
    along: Vector,
    len: f64,
    inward: Vector,
    lead: [f64; 2],
    /// The faces of the passage its prism may reach: its own wall, the
    /// faces at its ends, those its run on reaches, and walls facing it.
    reaches: Vec<usize>,
    /// The planes its prism is cut by, at the corners where it runs on to
    /// a wall not square to it: a point and the wall's normal into the air.
    trims: Vec<(Point, Vector)>,
}

/// A convex corner of the walls (the floor's rim turning outward): the
/// rolling ball pivots on the vertical edge at `at`, so the fill there is
/// the spandrel turned about that edge from the direction `from` (the
/// first wall's normal into the air) to `to` (the second's).
#[derive(Debug, Clone)]
struct Pivot {
    vertex: Shape,
    at: Point,
    from: Vector,
    to: Vector,
    reaches: Vec<usize>,
}

/// A narrow passage the union of its bands' fills rounds.
#[derive(Debug, Clone)]
pub(crate) struct Narrow {
    pub(crate) needs: Needs,
    up: Vector,
    radius: f64,
    bands: Vec<Band>,
    pivots: Vec<Pivot>,
    /// Each corner face's rim edge, with the bands either side of it.
    corners: Vec<(Id, [usize; 2])>,
    shapes: Vec<(Id, Shape)>,
}

/// A rim edge of the floor: along a planar wall, or round a corner face
/// square to the floor whose air lies on its concave side.
#[derive(Debug, Clone)]
enum Rim {
    Wall { ends: (Id, Id), inward: Vector },
    Corner { ends: (Id, Id), bend: Bend },
}

impl Rim {
    const fn ends(&self) -> (Id, Id) {
        match self {
            Self::Wall { ends, .. } | Self::Corner { ends, .. } => *ends,
        }
    }
}

/// A corner face's rim edge: the arc of `o + a cos(t) x + b sin(t) y`
/// over `range`, a circle where `a == b` (a corner a constant-radius round
/// can make) and an ellipse otherwise.
#[derive(Debug, Clone, Copy)]
struct Bend {
    o: Point,
    x: Vector,
    y: Vector,
    a: f64,
    b: f64,
    range: (f64, f64),
    circle: bool,
}

impl Bend {
    fn at(&self, t: f64) -> Point {
        self.o + self.x * (self.a * t.cos()) + self.y * (self.b * t.sin())
    }

    /// The unit tangent at `t`.
    fn tangent(&self, t: f64) -> Vector {
        let d = self.y * (self.b * t.cos()) - self.x * (self.a * t.sin());
        d / d.magnitude()
    }

    /// The radius of curvature at `t`.
    fn bend_radius(&self, t: f64) -> f64 {
        let (s, c) = t.sin_cos();
        let q = self.a * self.a * s * s + self.b * self.b * c * c;
        q * q.sqrt() / (self.a * self.b)
    }

    /// The greatest radius of curvature over the arc: `a^2 sin^2 t +
    /// b^2 cos^2 t` is greatest at an end or at a quarter turn between.
    fn loosest(&self) -> f64 {
        let mut worst = self
            .bend_radius(self.range.0)
            .max(self.bend_radius(self.range.1));
        let mut t = (self.range.0 / FRAC_PI_2).ceil() * FRAC_PI_2;
        while t < self.range.1 {
            worst = worst.max(self.bend_radius(t));
            t += FRAC_PI_2;
        }
        worst
    }

    /// The least and greatest of `(p - from) . d` over the arc.
    fn extent(&self, from: Point, d: Vector) -> (f64, f64) {
        let (cx, cy) = (self.a * self.x.dot(d), self.b * self.y.dot(d));
        let at = |t: f64| (self.o - from).dot(d) + cx * t.cos() + cy * t.sin();
        let (lo, hi) = self.range;
        let mut values = vec![at(lo), at(hi)];
        let peak = cy.atan2(cx);
        for base in [peak, peak + PI] {
            let mut t = base;
            while t < lo {
                t += TAU;
            }
            while t - TAU >= lo {
                t -= TAU;
            }
            if t <= hi {
                values.push(at(t));
            }
        }
        let least = values.iter().copied().fold(f64::INFINITY, f64::min);
        let most = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (least, most)
    }

    /// Whether the arc leaves its end at `p` square to `inward`: its
    /// tangent's lean toward `inward`, scaled by its radius of curvature
    /// there, within the length `slack`.
    fn tangent_at(&self, p: Point, inward: Vector, slack: f64) -> bool {
        let (lo, hi) = self.range;
        let t = if self.at(lo).distance(p) <= self.at(hi).distance(p) {
            lo
        } else {
            hi
        };
        self.tangent(t).dot(inward).abs() * self.bend_radius(t) <= slack
    }
}

/// The faces at vertex `v` other than those listed.
fn faces_at(curved: &Curved, v: Id, but: &[usize]) -> Vec<usize> {
    let mut found = Vec::new();
    for (fi, face) in curved.faces.iter().enumerate() {
        if but.contains(&fi) || found.contains(&fi) {
            continue;
        }
        let touches = face.loops.iter().flatten().any(|(e, _)| {
            curved
                .edges
                .get(e)
                .is_some_and(|edge| edge.ends.0 == v || edge.ends.1 == v)
        });
        if touches {
            found.push(fi);
        }
    }
    found
}

/// Reads a rim edge on `wall`, or `None` where it is neither kind.
fn rim_of(
    curved: &Curved,
    eid: Id,
    wall: usize,
    up: Vector,
    radius: f64,
    tol: Tolerances,
) -> Option<Rim> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let e = curved.edges.get(&eid)?;
    let (pa, pb) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    // A corner face tighter than the radius all along, its air on the
    // arc's concave side, where the arc's centre lies.
    let corner = |ends: (Id, Id), range: (f64, f64), f: Frame, (a, b): (f64, f64), circle| {
        let range = if range.0 <= range.1 {
            range
        } else {
            (range.1, range.0)
        };
        let bend = Bend {
            o: f.origin(),
            x: f.x().vector(),
            y: f.y().vector(),
            a,
            b,
            range,
            circle,
        };
        if bend.loosest() >= radius - eps {
            return None;
        }
        let p = bend.at(0.5 * (range.0 + range.1));
        let n = outward(curved, wall, p)?;
        (n.dot(bend.o - p) > 0.0).then_some(Rim::Corner { ends, bend })
    };
    match (&curved.faces[wall].surface, &e.curve) {
        (SurfaceGeometry::Plane(_), Curve::Line(_)) => {
            let n = outward(curved, wall, pa)?;
            let len = pa.distance(pb);
            if len <= eps {
                return None;
            }
            let along = (pb - pa) / len;
            (n.dot(up).abs() <= ang && along.dot(up).abs() <= ang && along.dot(n).abs() <= ang)
                .then_some(Rim::Wall {
                    ends: e.ends,
                    inward: n,
                })
        }
        (SurfaceGeometry::Cylinder(c), Curve::Circle(circle)) => {
            let circle = circle.forward(tol).ok()?.circle();
            let rc = circle.radius();
            if c.cylinder().axis().direction.vector().cross(up).magnitude() > ang
                || (rc - c.cylinder().radius()).abs() > eps
            {
                return None;
            }
            corner(e.ends, e.range, circle.frame(), (rc, rc), true)
        }
        (SurfaceGeometry::Extrusion(x), Curve::Ellipse(ellipse)) => {
            let ellipse = ellipse.forward(tol).ok()?.ellipse();
            let f = ellipse.frame();
            if x.direction().vector().cross(up).magnitude() > ang
                || f.z().vector().cross(up).magnitude() > ang
            {
                return None;
            }
            let axes = (ellipse.major_radius(), ellipse.minor_radius());
            corner(e.ends, e.range, f, axes, false)
        }
        _ => None,
    }
}

/// The narrow passage `edges` round at `radius`, or `None` where they are
/// not one this module builds.
///
/// # Errors
///
/// As reading the solid.
#[allow(
    clippy::too_many_lines,
    reason = "one walk round the rim, each end of each band read in turn"
)]
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Narrow>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    if edges.len() < 2 || !radius.is_finite() || radius <= eps {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let ids: Vec<Id> = edges.iter().map(id).collect();
    if ids.iter().any(|e| !curved.edges.contains_key(e)) {
        return Ok(None);
    }
    // The floor: the one face every edge of the request bounds.
    let mut common = faces_of(&curved, ids[0]);
    for e in &ids[1..] {
        let here = faces_of(&curved, *e);
        common.retain(|f| here.contains(f));
    }
    let [floor] = common[..] else {
        return Ok(None);
    };
    let SurfaceGeometry::Plane(_) = &curved.faces[floor].surface else {
        return Ok(None);
    };
    let Some(up) = outward(&curved, floor, curved.points[&curved.edges[&ids[0]].ends.0]) else {
        return Ok(None);
    };
    let mut rims = Vec::with_capacity(ids.len());
    let mut walls = Vec::with_capacity(ids.len());
    for (k, e) in ids.iter().enumerate() {
        let [f0, f1] = faces_of(&curved, *e)[..] else {
            return Ok(None);
        };
        let wall = if f0 == floor { f1 } else { f0 };
        if crate::setback::edge_is_convex(model, solid, &edges[k], radius, tol)? != Some(false) {
            return Ok(None);
        }
        let Some(rim) = rim_of(&curved, *e, wall, up, radius, tol) else {
            return Ok(None);
        };
        rims.push(rim);
        walls.push(wall);
    }
    // Which edges of the request meet at each vertex.
    let mut at: ogeom_core::FastMap<Id, Vec<usize>> = ogeom_core::FastMap::default();
    for (k, rim) in rims.iter().enumerate() {
        let (a, b) = rim.ends();
        at.entry(a).or_default().push(k);
        at.entry(b).or_default().push(k);
    }
    if at.values().any(|v| v.len() > 2) {
        return Ok(None);
    }
    let other = |v: Id, k: usize| -> Option<usize> {
        at.get(&v)
            .and_then(|list| list.iter().copied().find(|j| *j != k))
    };
    let mut needs = Needs::BlendMade;
    let mut tight = false;
    let mut bands: Vec<Band> = Vec::new();
    let mut walls_of: Vec<usize> = Vec::new();
    let mut pivots: Vec<Pivot> = Vec::new();
    let mut band_of: Vec<Option<usize>> = vec![None; rims.len()];
    let near = |v: Id| faces_at(&curved, v, &[floor]);
    for (k, rim) in rims.iter().enumerate() {
        let Rim::Wall { ends, inward } = rim else {
            continue;
        };
        let (pa, pb) = (curved.points[&ends.0], curved.points[&ends.1]);
        let along = (pb - pa) / pa.distance(pb);
        let mut ext = [0.0; 2];
        let mut reaches = vec![floor, walls[k]];
        // The planes the band runs on to at concave corners that are not
        // square: a point and the next wall's normal into the air.
        let mut trims: Vec<(Point, Vector)> = Vec::new();
        for (side, (v, out)) in [(ends.0, -along), (ends.1, along)].into_iter().enumerate() {
            let pv = curved.points[&v];
            reaches.extend(near(v));
            match other(v, k) {
                None => {
                    // An open end, flush on an end face square to the edge.
                    let [end] = faces_at(&curved, v, &[floor, walls[k]])[..] else {
                        return Ok(None);
                    };
                    let SurfaceGeometry::Plane(_) = &curved.faces[end].surface else {
                        return Ok(None);
                    };
                    if !outward(&curved, end, pv).is_some_and(|n| n.dot(out) >= 1.0 - ang) {
                        return Ok(None);
                    }
                }
                Some(j) => match &rims[j] {
                    Rim::Wall {
                        ends: next_ends,
                        inward: next,
                    } => {
                        let turn = next.dot(out);
                        if turn > ang {
                            // A convex corner: the ball pivots on the
                            // walls' edge there, and the band stops flush.
                            if k < j {
                                let Some(vertex) = curved.shapes.get(&v) else {
                                    return Ok(None);
                                };
                                pivots.push(Pivot {
                                    vertex: vertex.clone(),
                                    at: pv,
                                    from: *inward,
                                    to: *next,
                                    reaches: [vec![floor], near(v)].concat(),
                                });
                            }
                        } else if turn < -ang {
                            // A concave corner at any angle: the band runs
                            // on to the next wall's plane.
                            ext[side] = (-radius * next.dot(*inward) / turn).max(0.0);
                            if next.dot(*inward).abs() > ang {
                                trims.push((pv, *next));
                            }
                            // Run on past the corner, the band reaches past
                            // the next wall's far end only where the corner
                            // there is concave too: the run on is then behind
                            // the wall beyond, in the solid. Otherwise it must
                            // stand against the next wall's own face.
                            let w = if next_ends.0 == v {
                                next_ends.1
                            } else {
                                next_ends.0
                            };
                            let pw = curved.points[&w];
                            let past = other(w, j).map(|m| &rims[m]);
                            if ext[side] > eps {
                                if let Some(Rim::Wall { inward: beyond, .. }) = past
                                    && beyond.dot(pw - pv) < 0.0
                                {
                                    reaches.extend(near(w));
                                } else {
                                    let reach = pv + out * ext[side] + *inward * radius;
                                    if (reach - pv).dot(pw - pv) > pv.distance(pw).powi(2) + eps {
                                        return Ok(None);
                                    }
                                }
                            }
                        } else if next.dot(*inward) < 1.0 - ang {
                            return Ok(None);
                        }
                    }
                    Rim::Corner {
                        ends: arc_ends,
                        bend,
                    } => {
                        let w = if arc_ends.0 == v {
                            arc_ends.1
                        } else {
                            arc_ends.0
                        };
                        let Some(m) = other(w, j) else {
                            return Ok(None);
                        };
                        let Rim::Wall { inward: beyond, .. } = &rims[m] else {
                            return Ok(None);
                        };
                        reaches.extend(near(w));
                        let pw = curved.points[&w];
                        let tangent = bend.tangent_at(pv, *inward, eps * 10.0)
                            && bend.tangent_at(pw, *beyond, eps * 10.0);
                        // No constant-radius round makes an elliptic corner.
                        if !bend.circle {
                            needs = needs.max(Needs::WholeFace);
                        }
                        let square = beyond.dot(out) <= -1.0 + ang;
                        let facing = beyond.dot(*inward) <= -1.0 + ang;
                        if square {
                            // The corner face stands within the square
                            // the two bands cover, on their air side.
                            let (lo_i, hi_i) = bend.extent(pv, *inward);
                            let (lo_m, hi_m) = bend.extent(pw, *beyond);
                            if lo_i < -eps || hi_i > radius || lo_m < -eps || hi_m > radius {
                                return Ok(None);
                            }
                            ext[side] = (pw - pv).dot(out);
                            if !tangent {
                                needs = needs.max(Needs::WholeFace);
                            }
                        } else if facing && tangent {
                            // The band runs on to the far side of the
                            // face closing the zone.
                            ext[side] = bend.extent(pv, out).1;
                        } else {
                            return Ok(None);
                        }
                        tight = true;
                    }
                },
            }
        }
        let band = Band {
            edge: ids[k],
            start: pa - along * ext[0],
            along,
            len: pa.distance(pb) + ext[0] + ext[1],
            inward: *inward,
            lead: ext,
            reaches,
            trims,
        };
        band_of[k] = Some(bands.len());
        bands.push(band);
        walls_of.push(walls[k]);
    }
    // Facing walls: twice the radius apart or closer is a narrow passage.
    // Where they face each other over part of a wall's edge only, the zone
    // covers part of that wall's face.
    let mut facing = Vec::new();
    for (i, a) in bands.iter().enumerate() {
        for (j, b) in bands.iter().enumerate().skip(i + 1) {
            if a.inward.dot(b.inward) > -1.0 + ang {
                continue;
            }
            let span = |band: &Band| {
                let (s0, s1) = (
                    (band.start + band.along * band.lead[0] - a.start).dot(a.along),
                    (band.start + band.along * (band.len - band.lead[1]) - a.start).dot(a.along),
                );
                (s0.min(s1), s0.max(s1))
            };
            let (u0, u1) = span(a);
            let (w0, w1) = span(b);
            let overlap = u1.min(w1) - u0.max(w0);
            if overlap <= eps {
                continue;
            }
            let d = (b.start - a.start).dot(a.inward);
            if d <= 0.0 {
                continue;
            }
            if d <= eps {
                return Ok(None);
            }
            facing.push((i, j));
            // Exactly twice the radius apart, the contacts meet on one
            // line of the floor: the trough the two rounds make is one half
            // circle, and the floor between is taken whole. Closer than the
            // radius, each band crosses the wall it faces and is cut there.
            if d < 2.0 * radius + eps {
                tight = true;
                if overlap < u1 - u0 - eps || overlap < w1 - w0 - eps {
                    needs = needs.max(Needs::PartialFace);
                }
            }
        }
    }
    for (i, j) in facing {
        bands[i].reaches.push(walls_of[j]);
        bands[j].reaches.push(walls_of[i]);
    }
    if !tight || bands.is_empty() {
        return Ok(None);
    }
    let mut corners = Vec::new();
    for (k, rim) in rims.iter().enumerate() {
        if let Rim::Corner { ends, .. } = rim {
            let beside = [other(ends.0, k), other(ends.1, k)];
            let [Some(i), Some(j)] = beside.map(|b| b.and_then(|x| band_of[x])) else {
                return Ok(None);
            };
            corners.push((ids[k], [i, j]));
        }
    }
    let narrow = Narrow {
        needs,
        up,
        radius,
        bands,
        pivots,
        corners,
        shapes: ids.iter().copied().zip(edges.iter().cloned()).collect(),
    };
    let ours: Vec<usize> = {
        let mut ours = vec![floor];
        ours.extend(walls.iter().copied());
        for v in at.keys() {
            if at[v].len() == 1 {
                ours.extend(faces_at(&curved, *v, &[floor]));
            }
        }
        ours.sort_unstable();
        ours.dedup();
        ours
    };
    if !clear(model, &curved, &narrow, &ours, &ids, tol)? {
        return Ok(None);
    }
    Ok(Some(narrow))
}

impl Band {
    /// The eight corners of the box the band's prism stands in, the bit
    /// `4` set at its end, `2` across the floor and `1` up the wall.
    fn corners(&self, up: Vector, radius: f64) -> [Point; 8] {
        core::array::from_fn(|i| {
            let t = if i & 4 == 0 { 0.0 } else { self.len };
            let s = if i & 2 == 0 { 0.0 } else { radius };
            let z = if i & 1 == 0 { 0.0 } else { radius };
            self.start + self.along * t + self.inward * s + up * z
        })
    }

    /// The box the band's prism stands in.
    fn room(&self, up: Vector, radius: f64, margin: f64) -> Option<Aabb> {
        padded(&self.corners(up, radius), margin)
    }
}

impl Pivot {
    /// The cube about the corner the turned spandrel stands in.
    fn room(&self, radius: f64, margin: f64) -> Option<Aabb> {
        let h = Vector::new(radius, radius, radius);
        padded(&[self.at - h, self.at + h], margin)
    }
}

/// The box of `points` padded by `margin` all round.
fn padded(points: &[Point], margin: f64) -> Option<Aabb> {
    let room = Aabb::of_points(points);
    let (lo, hi) = (room.low()?, room.high()?);
    let pad = Vector::new(margin, margin, margin);
    Some(Aabb::of_points(&[lo - pad, hi + pad]))
}

/// Whether face `fi` stands in `room`: its box meets it, and where it is
/// planar, its plane passes through it.
fn stands_in(
    model: &Model,
    curved: &Curved,
    fi: usize,
    room: &Aabb,
    tol: Tolerances,
) -> OgeomResult<bool> {
    if !room.intersects(&ogeom_algo::shape_bounds(
        model,
        &curved.faces[fi].shape,
        tol,
    )?) {
        return Ok(false);
    }
    let SurfaceGeometry::Plane(plane) = &curved.faces[fi].surface else {
        return Ok(true);
    };
    let (Some(lo), Some(hi)) = (room.low(), room.high()) else {
        return Ok(true);
    };
    let (o, n) = (
        plane.plane().frame().origin(),
        plane.plane().normal().vector(),
    );
    let sides: Vec<f64> = (0..8)
        .map(|i| {
            let p = Point::new(
                if i & 1 == 0 { lo.x } else { hi.x },
                if i & 2 == 0 { lo.y } else { hi.y },
                if i & 4 == 0 { lo.z } else { hi.z },
            );
            n.dot(p - o)
        })
        .collect();
    Ok(sides.iter().any(|v| *v > 0.0) && sides.iter().any(|v| *v < 0.0))
}

/// Whether nothing of the solid but the passage's own faces (`ours`)
/// stands in any band's or pivot's room, none of their edges but the
/// request's and those between two of them, and of the passage's own faces
/// only those each piece may reach: a prism passes into the solid only
/// through the walls it reaches, never on into air beyond.
fn clear(
    model: &Model,
    curved: &Curved,
    narrow: &Narrow,
    ours: &[usize],
    asked: &[Id],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let margin = narrow.radius * 1e-3 + tol.confusion() * 10.0;
    // Each piece's room, the faces it may reach, and the directions its
    // sections run from the wall into the air.
    let mut rooms: Vec<(Aabb, &[usize], Vec<Vector>)> = Vec::new();
    for band in &narrow.bands {
        let Some(room) = band.room(narrow.up, narrow.radius, margin) else {
            return Ok(false);
        };
        rooms.push((room, &band.reaches, vec![band.inward]));
    }
    for pivot in &narrow.pivots {
        let Some(room) = pivot.room(narrow.radius, margin) else {
            return Ok(false);
        };
        rooms.push((room, &pivot.reaches, vec![pivot.from, pivot.to]));
    }
    let meets = |b: &Aabb| rooms.iter().any(|(r, _, _)| r.intersects(b));
    for (fi, face) in curved.faces.iter().enumerate() {
        if ours.contains(&fi) {
            for (room, reaches, rays) in &rooms {
                if reaches.contains(&fi) || !stands_in(model, curved, fi, room, tol)? {
                    continue;
                }
                // Another wall of the passage in the room: a section that
                // passes into the solid through a wall facing it must not
                // come out into air again, so no such wall may face the
                // way the section runs.
                let SurfaceGeometry::Plane(plane) = &face.surface else {
                    return Ok(false);
                };
                let o = plane.plane().frame().origin();
                let Some(n) = outward(curved, fi, o) else {
                    return Ok(false);
                };
                if rays.iter().any(|d| n.dot(*d) > tol.angular()) {
                    return Ok(false);
                }
            }
            continue;
        }
        if meets(&ogeom_algo::shape_bounds(model, &face.shape, tol)?) {
            return Ok(false);
        }
    }
    for &fi in ours {
        for (e, _) in curved.faces[fi].loops.iter().flatten() {
            if asked.contains(e) {
                continue;
            }
            let between = faces_of(curved, *e);
            if between.len() == 2 && between.iter().all(|f| ours.contains(f)) {
                continue;
            }
            if meets(&ogeom_algo::shape_bounds(
                model,
                &curved.edges[e].shape,
                tol,
            )?) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The band's fill: its spandrel (the region between the wall, the floor
/// and the round's quarter circle) swept along the band; the prism and
/// its round face.
fn prism(
    model: &mut Model,
    band: &Band,
    up: Vector,
    r: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape)>> {
    let e = band.start;
    let on_floor = e + band.inward * r;
    let on_wall = e + up * r;
    let centre = e + (up + band.inward) * r;
    let corner = ogeom_algo::make_vertex(model, e).shape;
    let foot = ogeom_algo::make_vertex(model, on_floor).shape;
    let top = ogeom_algo::make_vertex(model, on_wall).shape;
    let (c0, r0) = segment(e, on_floor, tol)?;
    let first = ogeom_algo::make_edge_between(model, c0, r0, &corner, &foot, tol)?.shape;
    let (c1, r1) = arc(centre, up.cross(band.inward), -up, r, FRAC_PI_2, tol)?;
    let round = ogeom_algo::make_edge_between(model, c1, r1, &foot, &top, tol)?.shape;
    let (c2, r2) = segment(on_wall, e, tol)?;
    let last = ogeom_algo::make_edge_between(model, c2, r2, &top, &corner, tol)?.shape;
    let wire = ogeom_algo::make_wire(model, &[first, round, last], tol)?.shape;
    let plane = PlaneSurface::new(Plane::new(Frame::new(
        e,
        Direction::new(band.along, tol)?,
        Direction::new(up, tol)?,
        tol,
    )?));
    let face = ogeom_algo::make_face(model, SurfaceGeometry::Plane(plane), &[wire], tol)?.shape;
    let mut solid = ogeom_algo::make_prism(model, &face, band.along * band.len, tol)?.shape;
    // Each trim takes off what lies behind the next wall's plane: a box
    // on that side, wider and deeper than the prism.
    let span = band.len + 4.0 * r;
    for (at, n) in &band.trims {
        let into = Direction::new(-*n, tol)?;
        // The walls stand square to the floor: `up` lies in each plane.
        let x = Direction::new(up, tol)?;
        let frame = Frame::new(*at, into, x, tol)?;
        let (fx, fy) = (frame.x().vector(), frame.y().vector());
        let corner = *at - (fx + fy) * span;
        let frame = Frame::new(corner, into, x, tol)?;
        let block = ogeom_algo::make_box(model, frame, (2.0 * span, 2.0 * span, span), tol)?.shape;
        match ogeom_bool::cut(model, &solid, &block, tol) {
            Ok(cut) => solid = cut.shape,
            Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
            Err(_) => return Ok(None),
        }
    }
    let rounds = round_faces(model, &solid, |g| matches!(g, SurfaceGeometry::Cylinder(_)))?;
    let [round_face] = &rounds[..] else {
        return Ok(None);
    };
    Ok(Some((solid, round_face.clone())))
}

/// The faces of `solid` whose surface `round` picks.
fn round_faces(
    model: &Model,
    solid: &Shape,
    round: impl Fn(&SurfaceGeometry) -> bool,
) -> OgeomResult<Vec<Shape>> {
    let mut found = Vec::new();
    for f in explore_unique(model, solid, ShapeType::Face)? {
        let data = model.node(&f).and_then(|n| n.data().as_face());
        if data
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(&round)
        {
            found.push(f);
        }
    }
    Ok(found)
}

/// The fill at a convex corner: the band's spandrel turned about the
/// walls' vertical edge from `from` to `to`; the solid and its round face.
fn turned(
    model: &mut Model,
    pivot: &Pivot,
    up: Vector,
    r: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape)>> {
    let e = pivot.at;
    let (on_floor, on_axis) = (e + pivot.from * r, e + up * r);
    let centre = e + (up + pivot.from) * r;
    let corner = ogeom_algo::make_vertex(model, e).shape;
    let foot = ogeom_algo::make_vertex(model, on_floor).shape;
    let top = ogeom_algo::make_vertex(model, on_axis).shape;
    let (c0, r0) = segment(e, on_floor, tol)?;
    let first = ogeom_algo::make_edge_between(model, c0, r0, &corner, &foot, tol)?.shape;
    let (c1, r1) = arc(centre, up.cross(pivot.from), -up, r, FRAC_PI_2, tol)?;
    let round = ogeom_algo::make_edge_between(model, c1, r1, &foot, &top, tol)?.shape;
    let (c2, r2) = segment(on_axis, e, tol)?;
    let last = ogeom_algo::make_edge_between(model, c2, r2, &top, &corner, tol)?.shape;
    let wire = ogeom_algo::make_wire(model, &[first, round, last], tol)?.shape;
    let plane = PlaneSurface::new(Plane::new(Frame::new(
        e,
        Direction::new(up.cross(pivot.from), tol)?,
        Direction::new(pivot.from, tol)?,
        tol,
    )?));
    let face = ogeom_algo::make_face(model, SurfaceGeometry::Plane(plane), &[wire], tol)?.shape;
    let turn = pivot.from.angle(pivot.to, tol)?;
    let axis = if pivot.from.cross(pivot.to).dot(up) >= 0.0 {
        up
    } else {
        -up
    };
    let solid = ogeom_algo::make_revolution(
        model,
        &face,
        Axis {
            location: e,
            direction: Direction::new(axis, tol)?,
        },
        turn,
        tol,
    )?;
    let rounds = round_faces(model, &solid.shape, |g| {
        !matches!(g, SurfaceGeometry::Plane(_))
    })?;
    let [round_face] = &rounds[..] else {
        return Ok(None);
    };
    let round_face = round_face.clone();
    Ok(Some((solid.shape, round_face)))
}

/// The passage rounded: each band's prism fused into `solid`; `None`
/// where a prism will not build or the result will not prove.
///
/// # Errors
///
/// As the booleans, other than their refusals.
pub(crate) fn built(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    narrow: &Narrow,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let deflection = ogeom_mesh::Deflection::default();
    let shape_of = |e: Id| {
        narrow
            .shapes
            .iter()
            .find(|(i, _)| *i == e)
            .map(|(_, s)| s.clone())
    };
    let mut held = 0.0;
    let mut tools = Vec::with_capacity(narrow.bands.len());
    for (k, band) in narrow.bands.iter().enumerate() {
        let Some((fill, round)) = prism(model, band, narrow.up, narrow.radius, tol)? else {
            return Ok(None);
        };
        held += ogeom_algo::volume_properties(model, &fill, deflection, tol)?.mass;
        let mut faces = Vec::new();
        if let Some(edge) = shape_of(band.edge) {
            faces.push((round.clone(), edge));
        }
        for (corner, beside) in &narrow.corners {
            if beside.contains(&k)
                && let Some(edge) = shape_of(*corner)
            {
                faces.push((round.clone(), edge));
            }
        }
        tools.push(super::fallback::Tool {
            solid: fill,
            sense: super::fallback::Sense::Fuse,
            faces,
        });
    }
    for pivot in &narrow.pivots {
        let Some((fill, round)) = turned(model, pivot, narrow.up, narrow.radius, tol)? else {
            return Ok(None);
        };
        held += ogeom_algo::volume_properties(model, &fill, deflection, tol)?.mass;
        tools.push(super::fallback::Tool {
            solid: fill,
            sense: super::fallback::Sense::Fuse,
            faces: vec![(round, pivot.vertex.clone())],
        });
    }
    let before = ogeom_algo::volume_properties(model, solid, deflection, tol)?.mass;
    let Some(built) = super::fallback::applied(model, solid, edges, &tools, tol)? else {
        return Ok(None);
    };
    let built = super::merged(model, solid, built, tol)?;
    // The fill adds material and no more than its prisms hold.
    let after = ogeom_algo::volume_properties(model, &built.shape, deflection, tol)?.mass;
    let slack = 1e-6 * (1.0 + before.abs());
    if after - before <= slack || after - before > held + slack {
        return Ok(None);
    }
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    Ok(Some(built))
}
