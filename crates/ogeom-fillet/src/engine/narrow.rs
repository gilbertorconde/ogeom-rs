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
//! fills: one spandrel prism per wall, each fused into the solid. The
//! faces of the tight zone (a corner face, or a floor narrower than twice
//! the radius) are buried where the fill covers them.
//!
//! Built for a planar floor with planar walls square to it, joined at
//! square corners (sharp, or through a cylindrical corner face square to
//! the floor tighter than the radius) or across a half cylinder tangent to
//! two facing walls, and for open ends flush on end faces square to their
//! edge. Prisms are kept only where nothing else of the solid stands in
//! their room, the fuses add no more than the prisms hold, and the result
//! is valid.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Aabb, Circle, Direction, Frame, Plane, Point, Vector};
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
    /// round made: a cylinder meeting its walls at kinks.
    WholeFace,
}

/// One wall's band, run on through the tight zone at either end: its
/// prism starts at `start` and runs `len` along `along`, `inward` the
/// wall's outward normal (into the passage).
#[derive(Debug, Clone)]
struct Band {
    edge: Id,
    start: Point,
    along: Vector,
    len: f64,
    inward: Vector,
}

/// A narrow passage the union of its bands' fills rounds.
#[derive(Debug, Clone)]
pub(crate) struct Narrow {
    pub(crate) needs: Needs,
    up: Vector,
    radius: f64,
    bands: Vec<Band>,
    /// Each corner face's rim edge, with the bands either side of it.
    corners: Vec<(Id, [usize; 2])>,
    shapes: Vec<(Id, Shape)>,
}

/// A rim edge of the floor: along a planar wall, or round a cylindrical
/// corner face whose air lies inside it.
#[derive(Debug, Clone)]
enum Rim {
    Wall {
        ends: (Id, Id),
        inward: Vector,
    },
    Corner {
        ends: (Id, Id),
        circle: Circle,
        range: (f64, f64),
    },
}

impl Rim {
    const fn ends(&self) -> (Id, Id) {
        match self {
            Self::Wall { ends, .. } | Self::Corner { ends, .. } => *ends,
        }
    }
}

/// The least and greatest of `(p - from) . d` over the arc of `circle`
/// across `range`.
fn arc_extent(circle: &Circle, range: (f64, f64), from: Point, d: Vector) -> (f64, f64) {
    let f = circle.frame();
    let (o, x, y) = (f.origin(), f.x().vector(), f.y().vector());
    let rc = circle.radius();
    let at = |t: f64| (o - from).dot(d) + rc * (t.cos() * x.dot(d) + t.sin() * y.dot(d));
    let mut values = vec![at(range.0), at(range.1)];
    let peak = y.dot(d).atan2(x.dot(d));
    for base in [peak, peak + PI] {
        let mut t = base;
        while t < range.0 {
            t += TAU;
        }
        while t - TAU >= range.0 {
            t -= TAU;
        }
        if t <= range.1 {
            values.push(at(t));
        }
    }
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (lo, hi)
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
            let circle = circle.circle();
            let rc = circle.radius();
            if c.cylinder().axis().direction.vector().cross(up).magnitude() > ang
                || (rc - c.cylinder().radius()).abs() > eps
                || rc >= radius - eps
            {
                return None;
            }
            // The air lies inside the cylinder: its outward normal there
            // points to the axis.
            let mid = 0.5 * (e.range.0 + e.range.1);
            let f = circle.frame();
            let p = f.origin() + (f.x().vector() * mid.cos() + f.y().vector() * mid.sin()) * rc;
            let n = outward(curved, wall, p)?;
            if n.dot(f.origin() - p) <= 0.0 {
                return None;
            }
            let (lo, hi) = if e.range.0 <= e.range.1 {
                e.range
            } else {
                (e.range.1, e.range.0)
            };
            Some(Rim::Corner {
                ends: e.ends,
                circle,
                range: (lo, hi),
            })
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
    let mut band_of: Vec<Option<usize>> = vec![None; rims.len()];
    for (k, rim) in rims.iter().enumerate() {
        let Rim::Wall { ends, inward } = rim else {
            continue;
        };
        let (pa, pb) = (curved.points[&ends.0], curved.points[&ends.1]);
        let along = (pb - pa) / pa.distance(pb);
        let mut ext = [0.0; 2];
        for (side, (v, out)) in [(ends.0, -along), (ends.1, along)].into_iter().enumerate() {
            let pv = curved.points[&v];
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
                    Rim::Wall { inward: next, .. } => {
                        if next.dot(out) > -1.0 + ang {
                            return Ok(None);
                        }
                    }
                    Rim::Corner {
                        ends: arc_ends,
                        circle,
                        range,
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
                        let pw = curved.points[&w];
                        let rc = circle.radius();
                        let o = circle.frame().origin();
                        let tangent = o.distance(pv + *inward * rc) <= eps * 10.0
                            && o.distance(pw + *beyond * rc) <= eps * 10.0;
                        let square = beyond.dot(out) <= -1.0 + ang;
                        let facing = beyond.dot(*inward) <= -1.0 + ang;
                        if square {
                            // The corner face stands within the square
                            // the two bands cover, on their air side.
                            let (lo_i, hi_i) = arc_extent(circle, *range, pv, *inward);
                            let (lo_m, hi_m) = arc_extent(circle, *range, pw, *beyond);
                            if lo_i < -eps || hi_i > radius || lo_m < -eps || hi_m > radius {
                                return Ok(None);
                            }
                            ext[side] = (pw - pv).dot(out);
                            if !tangent {
                                needs = needs.max(Needs::WholeFace);
                            }
                        } else if facing && tangent {
                            ext[side] = rc;
                        } else {
                            return Ok(None);
                        }
                        tight = true;
                    }
                },
            }
        }
        band_of[k] = Some(bands.len());
        bands.push(Band {
            edge: ids[k],
            start: pa - along * ext[0],
            along,
            len: pa.distance(pb) + ext[0] + ext[1],
            inward: *inward,
        });
    }
    // Facing walls: twice the radius apart or closer is a narrow passage, and
    // a band must stop short of the wall it faces.
    for (i, a) in bands.iter().enumerate() {
        for b in &bands[i + 1..] {
            if a.inward.dot(b.inward) > -1.0 + ang {
                continue;
            }
            let span = |band: &Band| {
                let (s0, s1) = (
                    (band.start - a.start).dot(a.along),
                    (band.start + band.along * band.len - a.start).dot(a.along),
                );
                (s0.min(s1), s0.max(s1))
            };
            let (u0, u1) = span(a);
            let (w0, w1) = span(b);
            if u1.min(w1) - u0.max(w0) <= eps {
                continue;
            }
            let d = (b.start - a.start).dot(a.inward);
            if d <= 0.0 {
                continue;
            }
            if d <= radius + eps {
                return Ok(None);
            }
            // Exactly twice the radius apart, the contacts meet on one
            // line of the floor: the trough the two rounds make is one half
            // circle, and the floor between is taken whole.
            if d < 2.0 * radius + eps {
                tight = true;
            }
        }
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
    /// The box the band's prism stands in.
    fn room(&self, up: Vector, radius: f64, margin: f64) -> Option<Aabb> {
        let end = self.start + self.along * self.len;
        let mut points = Vec::with_capacity(8);
        for p in [self.start, end] {
            points.extend([
                p,
                p + self.inward * radius,
                p + up * radius,
                p + (self.inward + up) * radius,
            ]);
        }
        let room = Aabb::of_points(&points);
        let (lo, hi) = (room.low()?, room.high()?);
        let pad = Vector::new(margin, margin, margin);
        Some(Aabb::of_points(&[lo - pad, hi + pad]))
    }
}

/// Whether nothing of the solid but the passage's own faces (`ours`)
/// stands in any band's room, and none of their edges but the request's
/// and those between two of them.
fn clear(
    model: &Model,
    curved: &Curved,
    narrow: &Narrow,
    ours: &[usize],
    asked: &[Id],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let margin = narrow.radius * 1e-3 + tol.confusion() * 10.0;
    let mut rooms = Vec::with_capacity(narrow.bands.len());
    for band in &narrow.bands {
        let Some(room) = band.room(narrow.up, narrow.radius, margin) else {
            return Ok(false);
        };
        rooms.push(room);
    }
    let meets = |b: &Aabb| rooms.iter().any(|r| r.intersects(b));
    for (fi, face) in curved.faces.iter().enumerate() {
        if ours.contains(&fi) {
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
    let solid = ogeom_algo::make_prism(model, &face, band.along * band.len, tol)?.shape;
    let mut rounds = Vec::new();
    for f in explore_unique(model, &solid, ShapeType::Face)? {
        let data = model.node(&f).and_then(|n| n.data().as_face());
        let cylinder = data
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Cylinder(_)));
        if cylinder {
            rounds.push(f);
        }
    }
    let [round_face] = &rounds[..] else {
        return Ok(None);
    };
    Ok(Some((solid, round_face.clone())))
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
