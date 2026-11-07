//! The largest radius a round, or distance a chamfer, seats at: read from
//! where the ball's contacts walk across the faces as it grows, without
//! building a band.
//!
//! A ball of radius `r` on an edge whose faces turn by `φ` touches each
//! face `r · tan(φ/2)` back from the edge, square to it, at the stations
//! the request check holds it to. As `r` grows each contact walks across
//! its face, and the largest radius whose contacts all stay on their faces
//! is where the first walk leaves its face. On a plane, and along a ruling
//! of a cylinder or a cone, the walk is straight in the face's chart, and
//! where it leaves the region the check reads the face by is exact; on
//! other surfaces the walk curves, and is followed by stations along it.
//! Two bands overlap where a contact of one comes within the other's
//! setback of its edge: the contact's distance to that edge changes no
//! faster than the contact walks, so the radius where it first does is
//! approached from below in steps that cannot pass it.

use crate::support::{chart_holds, distance_to_edge, edge_curve, face_mesh, same_occurrence};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// Stations along the walk of a contact on a surface where it curves.
const WALK_STATIONS: u32 = 256;

/// Where one contact of an edge's ball walks as the radius grows: from the
/// edge's point `from`, along the unit `step` square to the edge into
/// `face`, `factor` (the tangent of half the faces' turn) per unit of
/// radius, the point found on `surface` (the face's, placed).
struct Walk {
    face: Shape,
    surface: SurfaceGeometry,
    from: Point,
    step: Vector,
    factor: f64,
}

impl Walk {
    /// The contact a setback `sigma` back from the edge.
    fn contact(&self, sigma: f64, tol: Tolerances) -> OgeomResult<(Point, Point2)> {
        let found =
            ogeom_algo::project_on_surface(&self.surface, self.from + self.step * sigma, 16, tol)?;
        let (u, v) = found.parameters;
        Ok((found.point, Point2::new(u, v)))
    }
}

/// How far a walk stays on its face: the setback where it first leaves
/// the face's region, at most the `up` it was followed to, and how fast
/// its contact moves per unit of setback.
struct Reach {
    setback: f64,
    speed: f64,
}

/// The walks of `edge`'s ball, two at each station [`ball_fits`] reads
/// (stations a quarter, half and three quarters along it), on each of the
/// edge's two faces of `solid`. None for an edge not between two faces;
/// stations where the faces meet tangent are left out, as the check leaves
/// them.
///
/// [`ball_fits`]: crate::support::ball_fits
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where the side
/// of the edge a face lies on reads differently at the hair's steps the
/// check reads it at: that walk's direction depends on the radius.
fn walks(model: &Model, solid: &Shape, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<Walk>> {
    let mut found = Vec::new();
    let Ok((curve, range)) = edge_curve(model, edge, tol) else {
        return Ok(found);
    };
    let mut sides: Vec<(Shape, bool)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        for e in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if same_occurrence(model, &e, edge, tol) {
                sides.push((face.clone(), e.orientation() == Orientation::Reversed));
                break;
            }
        }
    }
    let [(face0, rev0), (face1, rev1)] = sides.as_slice() else {
        return Ok(found);
    };
    let normal_at = |face: &Shape, p: Point| -> OgeomResult<Option<(SurfaceGeometry, Vector)>> {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return Ok(None);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let (u, v) = ogeom_algo::project_on_surface(&placed, p, 16, tol)?.parameters;
        let Ok(n) = placed.normal_at(u, v, tol) else {
            return Ok(None);
        };
        let n = if face.orientation() == Orientation::Reversed {
            -n.vector()
        } else {
            n.vector()
        };
        Ok(Some((placed, n)))
    };
    for fraction in [0.25, 0.5, 0.75] {
        let t = range.0 + (range.1 - range.0) * fraction;
        let p = curve.point_at(t, tol)?;
        let d = curve.d1_at(t, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            continue;
        }
        let tangent = d / m;
        let (Some((s0, n0)), Some((s1, n1))) = (normal_at(face0, p)?, normal_at(face1, p)?) else {
            continue;
        };
        let turn = n0.dot(n1).clamp(-1.0, 1.0).acos();
        if turn <= tol.angular() * 1e3 {
            continue;
        }
        let factor = (turn * 0.5).tan();
        // Into each face: square to the edge and the normal, on the side
        // a hair's step lands in the face, for every hair the check steps.
        let into = |face: &Shape, surface: &SurfaceGeometry, n: Vector, reversed: bool| {
            let walk = if reversed { -tangent } else { tangent };
            let d = n.cross(walk);
            let l = d.magnitude();
            if l <= 0.0 {
                return Ok(None);
            }
            let d = d / l;
            let lands = |step: Vector| -> OgeomResult<bool> {
                let (u, v) = ogeom_algo::project_on_surface(surface, p + step, 16, tol)?.parameters;
                chart_holds(model, face, Point2::new(u, v), tol)
            };
            // The check's hair grows with the setback, from ten confusion
            // distances to 1e5. Landing on both sides, it takes the first
            // guess, and holds to it while the contact that way is on the
            // face, as it is wherever the walk below holds.
            let mut side: Option<Vector> = None;
            for decade in 1..=5 {
                let hair = tol.confusion() * 10_f64.powi(decade);
                let here = match [lands(d * hair)?, lands(-d * hair)?] {
                    [true, _] => d,
                    [false, true] => -d,
                    [false, false] => continue,
                };
                if side.is_some_and(|s| s.dot(here) < 0.0) {
                    ogeom_bail!(
                        NotDone,
                        "the face beside the edge through {p:?} reads on different sides of \
                         it at the steps the check takes, which grow with the radius"
                    );
                }
                side = Some(here);
            }
            Ok(side)
        };
        let (Some(d0), Some(d1)) = (into(face0, &s0, n0, *rev0)?, into(face1, &s1, n1, *rev1)?)
        else {
            continue;
        };
        for (face, surface, step) in [(face0, s0, d0), (face1, s1, d1)] {
            found.push(Walk {
                face: face.clone(),
                surface,
                from: p,
                step,
                factor,
            });
        }
    }
    Ok(found)
}

/// How far `walk` stays on its face, followed up to the setback `up`.
///
/// On a plane, cylinder or cone, a walk that lies on the surface at its
/// start, middle and end is a line on it (a line meets a quadric in two
/// points unless it lies on it) and runs straight in the chart, and where
/// it leaves the face's region is found exactly, triangle by triangle of
/// the mesh the check reads the face by. Otherwise it is followed by
/// [`WALK_STATIONS`] stations and the first it leaves at is bisected.
fn reach(model: &Model, walk: &Walk, up: f64, tol: Tolerances) -> OgeomResult<Reach> {
    let conf = tol.confusion();
    let (p0, a) = walk.contact(0.0, tol)?;
    let (pm, m) = walk.contact(up * 0.5, tol)?;
    let (p1, b) = walk.contact(up, tol)?;
    let quadric = matches!(
        walk.surface,
        SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)
    );
    let on = |p: Point, sigma: f64| p.distance(walk.from + walk.step * sigma) <= conf;
    let mid = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
    let chart_span =
        a.x.abs()
            .max(a.y.abs())
            .max(b.x.abs())
            .max(b.y.abs())
            .max(1.0);
    if quadric
        && on(p0, 0.0)
        && on(pm, up * 0.5)
        && on(p1, up)
        && (m.x - mid.x).abs().max((m.y - mid.y).abs()) <= chart_span * 1e-12
    {
        let mesh = face_mesh(model, &walk.face, tol)?;
        return Ok(Reach {
            setback: up * straight_exit(&mesh, a, b),
            speed: 1.0,
        });
    }
    // A contact found where the surface's window clamps the search, off
    // square from the point walked to, stands past the surface's edge.
    let holds = |sigma: f64| -> OgeomResult<bool> {
        let target = walk.from + walk.step * sigma;
        let (at, uv) = walk.contact(sigma, tol)?;
        let Ok(n) = walk.surface.normal_at(uv.x, uv.y, tol) else {
            return Ok(false);
        };
        let off = target - at;
        let n = n.vector();
        if (off - n * n.dot(off)).magnitude() > conf {
            return Ok(false);
        }
        chart_holds(model, &walk.face, uv, tol)
    };
    let mut speed = 1.0_f64;
    let mut last = (0.0, p0);
    for k in 1..=WALK_STATIONS {
        let sigma = up * f64::from(k) / f64::from(WALK_STATIONS);
        if !holds(sigma)? {
            let (mut lo, mut hi) = (last.0, sigma);
            for _ in 0..48 {
                let mid = f64::midpoint(lo, hi);
                if holds(mid)? {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return Ok(Reach { setback: lo, speed });
        }
        let (at, _) = walk.contact(sigma, tol)?;
        speed = speed.max(at.distance(last.1) / (sigma - last.0) * 2.0);
        last = (sigma, at);
    }
    Ok(Reach { setback: up, speed })
}

/// Where the chart segment from `a` to `b` first leaves the region a
/// face's mesh lays over its chart, as a fraction of the way along it: 0
/// where it starts outside, 1 where it never leaves. A point is in the
/// region where it lies in a triangle or on one's side, as
/// [`chart_holds`] reads it.
fn straight_exit(mesh: &ogeom_topo::Triangulation, a: Point2, b: Point2) -> f64 {
    let chart = |i: u32| {
        let (u, v) = mesh.parameters[i as usize];
        Point2::new(u, v)
    };
    let span = mesh
        .parameters
        .iter()
        .fold(0.0_f64, |m, (u, v)| m.max(u.abs()).max(v.abs()))
        .max(1.0);
    let eps = span * 1e-9;
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    // Each side's cross product along the segment, alpha + beta s.
    let side = |p: Point2, q: Point2| {
        let (ex, ey) = (q.x - p.x, q.y - p.y);
        (ex * (a.y - p.y) - ey * (a.x - p.x), ex * dy - ey * dx)
    };
    // Where alpha + beta s >= bound over [0, 1].
    let at_least = |(alpha, beta): (f64, f64), bound: f64, (lo, hi): (f64, f64)| {
        let excess = alpha - bound;
        if beta == 0.0 {
            return if excess >= 0.0 { (lo, hi) } else { (1.0, 0.0) };
        }
        let root = -excess / beta;
        if beta > 0.0 {
            (lo.max(root), hi)
        } else {
            (lo, hi.min(root))
        }
    };
    let mut spans: Vec<(f64, f64)> = Vec::new();
    for t in &mesh.triangles {
        let [p, q, r] = t.map(chart);
        let sides = [side(p, q), side(q, r), side(r, p)];
        // Every side's cross product at least -eps, or every one at most
        // eps: in the triangle either way it winds.
        let up = sides
            .iter()
            .fold((0.0, 1.0), |span, s| at_least(*s, -eps, span));
        let down = sides
            .iter()
            .fold((0.0, 1.0), |span, s| at_least((-s.0, -s.1), -eps, span));
        for span in [up, down] {
            if span.0 <= span.1 {
                spans.push(span);
            }
        }
    }
    spans.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut reached: Option<f64> = None;
    for (lo, hi) in spans {
        match reached {
            None if lo <= 0.0 => reached = Some(hi),
            Some(r) if lo <= r => reached = Some(r.max(hi)),
            _ => {}
        }
    }
    reached.unwrap_or(0.0)
}

/// The largest measure up to `up_to` at which every contact of `edge`'s
/// ball is on its face at every station, where a contact stands `per` of
/// the measure back from the edge, with each walk and its reach.
fn edge_largest(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    up_to: f64,
    per: impl Fn(&Walk) -> f64,
    tol: Tolerances,
) -> OgeomResult<(f64, Vec<(Walk, Reach)>)> {
    let conf = tol.confusion();
    let mut largest = up_to;
    let mut found = Vec::new();
    for walk in walks(model, solid, edge, tol)? {
        let rate = per(&walk);
        if rate <= 0.0 {
            continue;
        }
        let bounds = ogeom_algo::shape_bounds(model, &walk.face, tol)?;
        // Past its face's diagonal from the edge, a contact is off it.
        let up = (up_to * rate).min(bounds.diagonal() + conf);
        let reached = reach(model, &walk, up, tol)?;
        if reached.setback < up || up < up_to * rate {
            largest = largest.min((reached.setback - conf) / rate);
        }
        found.push((walk, reached));
    }
    Ok((largest, found))
}

/// The largest measure up to `up_to` at which every contact of a bevel on
/// `edge` is on its face at every station, where the bevel's leg on a
/// face runs `per` of the measure along it.
pub(crate) fn legs_largest(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    up_to: f64,
    per: impl Fn(ogeom_topo::TShapeId) -> f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    crate::support::keeping_face_meshes(|| {
        Ok(edge_largest(model, solid, edge, up_to, |w| per(w.face.node()), tol)?.0)
    })
}

/// The largest radius up to `up_to` at which a ball rolling along `edge`
/// of `solid` sits with both contacts on their faces at every station the
/// request check reads, so that [`fillet_edges`](crate::fillet_edges) on
/// that edge alone is not refused for setback at it or at any radius
/// below it. None where no radius above the confusion distance fits.
///
/// Read from the seat, without building a band: each contact walks square
/// to the edge, `r · tan(φ/2)` into its face, and the answer is where the
/// first walk leaves its face, a confusion distance short. On a plane,
/// and along a ruling of a cylinder or a cone, that is exact against the
/// region the check reads the face by; on other surfaces the walk is
/// followed by stations, and a sliver of the face's outside narrower than
/// their spacing can pass between two. A contact the check lets run off a
/// face out of the solid, where another station holds the ball, bounds the
/// answer here all the same, so the answer can lie under the largest
/// radius the check lets through. Refusals other than for setback (a seat
/// the round does not speak, a corner) are outside its promise.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to` is not a usable length;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where which
/// way a contact walks into its face depends on the radius, or the radius
/// found does not pass the check (a station whose answer the walk does not
/// predict).
pub fn largest_fillet_radius(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    largest_fillet_radius_all(model, solid, std::slice::from_ref(edge), up_to, tol)
}

/// The largest radius up to `up_to` at which the bands of `edges` all sit
/// with their contacts on their faces, as [`largest_fillet_radius`] reads
/// each, and no band's contact comes within another's setback of its edge
/// across a face they share, the overlap
/// [`fillet_edges`](crate::fillet_edges) refuses. Edges meeting at a
/// vertex are corner mates, whose bands meet there by design, and are not
/// held apart. None where no radius above the confusion distance fits.
///
/// A contact's distance to another edge changes no faster than the
/// contact walks, and the edge's setback grows linearly, so the radius
/// where the two first meet is approached from below in steps that cannot
/// pass it.
///
/// # Errors
///
/// As [`largest_fillet_radius`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `edges` is empty.
pub fn largest_fillet_radius_all(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let conf = tol.confusion();
    if !up_to.is_finite() || up_to <= conf {
        ogeom_bail!(
            Construction,
            "a largest radius up to {up_to} asks for nothing"
        );
    }
    if edges.is_empty() {
        ogeom_bail!(Construction, "a request of no edges has no largest radius");
    }
    crate::support::keeping_face_meshes(|| -> OgeomResult<Option<f64>> {
        let mut largest = up_to;
        let mut seats = Vec::with_capacity(edges.len());
        for edge in edges {
            let (radius, walks) = edge_largest(model, solid, edge, up_to, |w| w.factor, tol)?;
            largest = largest.min(radius);
            seats.push(walks);
        }
        if largest <= conf {
            return Ok(None);
        }
        let ends = |edge: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
            Ok(ogeom_algo::edge_vertices(model, edge)?
                .map(|(a, b)| vec![a.node(), b.node()])
                .unwrap_or_default())
        };
        for (i, first) in edges.iter().enumerate() {
            for (j, second) in edges.iter().enumerate() {
                if i == j {
                    continue;
                }
                let (own, other) = (ends(first)?, ends(second)?);
                if own.iter().any(|v| other.contains(v)) {
                    continue;
                }
                let Ok((curve, range)) = edge_curve(model, second, tol) else {
                    continue;
                };
                for (walk, reached) in &seats[i] {
                    // The other edge's widest setback per unit radius on
                    // this face, if it has a contact there.
                    let Some(strip) = seats[j]
                        .iter()
                        .filter(|(w, _)| w.face.node() == walk.face.node())
                        .map(|(w, _)| w.factor)
                        .reduce(f64::max)
                    else {
                        continue;
                    };
                    // The gap between the contact and the other edge's
                    // strip, shrinking at most this fast per unit radius.
                    let rate = walk.factor * reached.speed + strip;
                    let gap = |r: f64| -> OgeomResult<f64> {
                        let (at, _) = walk.contact(r * walk.factor, tol)?;
                        Ok(distance_to_edge(&curve, range, at, tol)? - strip * r)
                    };
                    let mut r = 0.0_f64;
                    for _ in 0..100_000 {
                        if r >= largest {
                            break;
                        }
                        let g = gap(r)?;
                        // Where the gap closes to a confusion distance, a
                        // confusion distance short.
                        if g <= conf {
                            largest = (r - conf / rate).max(0.0);
                            break;
                        }
                        r += g / rate;
                    }
                    largest = largest.min(r);
                }
            }
        }
        if largest <= conf {
            return Ok(None);
        }
        // The answer passes the request's own check, or it is not one.
        if let Err(refusal) = crate::fillet::request_fits(model, solid, edges, largest, tol) {
            ogeom_bail!(
                NotDone,
                "the largest radius read from the seats, {largest}, does not pass the \
                 request's check: {refusal}"
            );
        }
        Ok(Some(largest))
    })
}
