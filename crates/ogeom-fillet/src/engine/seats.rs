//! A round's ball re-seated where its contact leaves its face (ADR 0001,
//! section 13, question 8).
//!
//! The ball touches the solid's boundary on each side at a seat: a face, an
//! edge of that face, or a vertex. Its centre runs where the two seats'
//! offsets meet, and the round is the envelope of those balls. Where a
//! contact crosses its face's edge `E` (a [`RollOver`]), the ball's seat on
//! that side moves from the face to `E`: its centre stays a radius off the
//! other face and comes a radius from `E`, on the pipe about it. The face's
//! offset and `E`'s pipe are tangent along the curve a radius in from `E`,
//! so the centre's path runs on tangent across the change, and the two
//! pieces of the round share the section there.
//!
//! [`edge_run`] reads one such run, from the crossing where the contact
//! leaves its face to the one where it comes back, station by station on
//! the strip's own spine: each station's plane (square to the spine there)
//! holds the centre, found where both distances are the radius. The run's
//! surface is the balls' arcs fitted along the stations; the contact on the
//! other face and the centre's path are fitted through the stations too; the
//! seat on `E` is `E` itself between the crossings, exact.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};

use super::Id;
use super::curved::Curved;
use super::setbacks::RollOver;

/// One run of a round's ball resting on an edge `over` of the face its
/// contact left, between the two crossings: what [`edge_run`] reads.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRun {
    /// The edges the ball rests on, in the run's order: one, or several
    /// running on tangent into each other (a circle a seam splits).
    pub over: Vec<Id>,
    /// Where on them each ball touches: the edge's place in `over` and its
    /// parameter on that edge's curve.
    pub on_over: Vec<(usize, f64)>,
    /// The ball's centre at each station, from the first crossing's to the
    /// second's.
    pub centres: Vec<Point>,
    /// Where each ball touches the other face.
    pub feet: Vec<Point>,
    /// Where each ball touches `over`.
    pub on_edge: Vec<Point>,
    /// The balls' arcs fitted along the stations, from the other face's
    /// contact across to `over`.
    pub surface: SurfaceGeometry,
    /// The contact on the other face, fitted through the feet.
    pub contact: Curve,
    /// How far the fitted surface and contact stand off the exact
    /// stations at most.
    pub off: f64,
}

/// Why a run is not built: a seat change this milestone refuses by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NotSeated {
    /// The two crossings are over different edges: the ball would pass a
    /// vertex, a seat of its own.
    PastAVertex,
    /// The ball's touch on the edge reaches the edge's end.
    OffTheEdge,
    /// The ball comes down on the face beyond the edge.
    OntoTheFaceBeyond,
    /// No centre a radius from both seats at a station.
    NoCentre,
}

impl NotSeated {
    /// The refusal's name.
    pub(crate) const fn name(self) -> &'static str {
        match self {
            Self::PastAVertex => {
                "a round's ball would roll past a vertex of the edge it rests on, a seat not built"
            }
            Self::OffTheEdge => {
                "a round's ball would run off the end of the edge it rests on, a seat not built"
            }
            Self::OntoTheFaceBeyond => {
                "a round's ball would come down on the face beyond the edge it rests on"
            }
            Self::NoCentre => "a round's ball cannot rest on the edge its contact leaves over",
        }
    }
}

/// How close the fit holds the run to its exact stations.
const RUN_FIT: f64 = 2e-7;

/// The stations a run is first read at; doubled while the fit misses.
const FIRST_STATIONS: usize = 33;

/// The most stations a run is read at.
const MOST_STATIONS: usize = 1025;

/// The points each ball's arc is read at across the run.
const ACROSS: usize = 17;

/// The run of the round whose face-to-face centre path is `spine` (its
/// station planes square to it), radius `radius`, its contact on face
/// `left` leaving over `enter.over` at `enter` and coming back at `leave`,
/// its other face `other`: the ball resting on the edge between.
///
/// `Ok(Err(..))` names a seat change not built here.
///
/// # Errors
///
/// As evaluating and fitting the curves and surfaces.
#[allow(clippy::too_many_arguments, reason = "the run's whole reading")]
pub(crate) fn edge_run(
    model: &ogeom_topo::Model,
    curved: &Curved,
    spine: &Curve,
    radius: f64,
    left: usize,
    other: usize,
    enter: &RollOver,
    leave: &RollOver,
    tol: Tolerances,
) -> OgeomResult<Result<EdgeRun, NotSeated>> {
    let Some(chain) = seat_chain(curved, left, enter.over, leave.over, tol)? else {
        return Ok(Err(NotSeated::PastAVertex));
    };
    let curves: Vec<Curve> = chain
        .iter()
        .map(|e| {
            let c = &curved.edges[e];
            trimmed(&c.curve, c.range, tol)
        })
        .collect::<OgeomResult<_>>()?;
    let nearest = |p: Point| -> OgeomResult<(usize, ogeom_algo::Projection)> {
        let mut best: Option<(usize, ogeom_algo::Projection)> = None;
        for (k, c) in curves.iter().enumerate() {
            let near = ogeom_algo::project_on_curve(c, p, 32, tol)?;
            if best
                .as_ref()
                .is_none_or(|(_, b)| near.distance < b.distance)
            {
                best = Some((k, near));
            }
        }
        best.ok_or_else(|| ogeom_core::OgeomError::Construction("a seat of no edge".into()))
    };
    let face = &curved.faces[other].surface;
    let beyond: Vec<(&SurfaceGeometry, ogeom_algo::PointDistance)> = (0..curved.faces.len())
        .filter(|fi| *fi != left)
        .filter(|fi| {
            curved.faces[*fi]
                .loops
                .iter()
                .flatten()
                .any(|(e, _)| chain.contains(e))
        })
        .map(|fi| {
            let f = &curved.faces[fi];
            Ok((
                &f.surface,
                ogeom_algo::PointDistance::of(model, &f.shape, tol)?,
            ))
        })
        .collect::<OgeomResult<_>>()?;
    let (s0, s1) = (
        super::trimmed::param_on(spine, spine.domain(), enter.point, tol)?,
        super::trimmed::param_on(spine, spine.domain(), leave.point, tol)?,
    );
    let mut n = FIRST_STATIONS;
    loop {
        let mut centres = Vec::with_capacity(n);
        let mut feet = Vec::with_capacity(n);
        let mut on_edge = Vec::with_capacity(n);
        let mut on_params = Vec::with_capacity(n);
        let mut guess: Option<Point> = None;
        for i in 0..n {
            #[expect(clippy::cast_precision_loss, reason = "a few stations")]
            let s = s0 + (s1 - s0) * (i as f64) / ((n - 1) as f64);
            let o = spine.point_at(s, tol)?;
            let t = spine.d1_at(s, tol)?.normalized(tol)?;
            let start = guess.unwrap_or(o);
            let Some(c) = centre_in(face, &nearest, o, t, start, radius, tol)? else {
                return Ok(Err(NotSeated::NoCentre));
            };
            let foot = ogeom_geom::project::project_on_surface(face, c, 16, tol)?;
            let (k, near) = nearest(c)?;
            // The touch on the chain inside it, short of its outer ends.
            let (lo, hi) = curves[k].domain();
            let first_end = if forward_of(curved, &chain, 0) {
                lo
            } else {
                hi
            };
            let last_end = if forward_of(curved, &chain, chain.len() - 1) {
                hi
            } else {
                lo
            };
            if i > 0
                && i < n - 1
                && ((k == 0 && (near.parameter - first_end).abs() <= tol.parametric())
                    || (k == chain.len() - 1
                        && (near.parameter - last_end).abs() <= tol.parametric()))
            {
                return Ok(Err(NotSeated::OffTheEdge));
            }
            // The ball comes down on a face beyond the edge where that
            // face's own foot of the centre stands on the face, nearer than
            // the radius.
            for (surface, face) in &beyond {
                let at = ogeom_geom::project::project_on_surface(surface, c, 16, tol)?;
                if at.distance < radius - tol.confusion()
                    && face.nearest(at.point, tol)?.0 <= tol.confusion()
                {
                    return Ok(Err(NotSeated::OntoTheFaceBeyond));
                }
            }
            centres.push(c);
            feet.push(foot.point);
            on_edge.push(near.point);
            on_params.push((k, near.parameter));
            guess = Some(c);
        }
        let target = RUN_FIT;
        let rows: Vec<Vec<Point>> = (0..n)
            .map(|i| arc_points(centres[i], feet[i], on_edge[i], radius, tol))
            .collect::<OgeomResult<_>>()?;
        let surface = ogeom_geom::fit::fit_surface_grid_chordal(&rows, 3, target, tol)?;
        let contact_fit = ogeom_geom::fit::fit_points(&feet, 3, target, tol)?;
        let met = surface.met && contact_fit.met;
        if met || n >= MOST_STATIONS {
            if !met {
                ogeom_core::ogeom_bail!(
                    NotDone,
                    "a round's ball resting on an edge fits within {} of its stations, against {target}",
                    surface.error.max(contact_fit.error)
                );
            }
            let off = surface.error.max(contact_fit.error);
            return Ok(Ok(EdgeRun {
                over: chain,
                on_over: on_params,
                centres,
                feet,
                on_edge,
                surface: SurfaceGeometry::BSpline(surface.curve),
                contact: Curve::BSpline(contact_fit.curve),
                off,
            }));
        }
        n = 2 * n - 1;
    }
}

/// The edges of face `face`'s loop from `first` to `last`, each running on
/// tangent into the next at the vertex they share, the shorter way round
/// the loop that does: the seat a ball rests on across a seam. `None`
/// where neither way runs on tangent (the ball would pass a vertex).
///
/// # Errors
///
/// As evaluating the edges' curves.
fn seat_chain(
    curved: &Curved,
    face: usize,
    first: Id,
    last: Id,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Id>>> {
    if first == last {
        return Ok(Some(vec![first]));
    }
    let Some(walk) = curved.faces[face]
        .loops
        .iter()
        .find(|w| w.iter().any(|(e, _)| *e == first))
    else {
        return Ok(None);
    };
    let n = walk.len();
    let start = walk.iter().position(|(e, _)| *e == first).unwrap_or(0);
    let mut best: Option<Vec<Id>> = None;
    for step in [1, n - 1] {
        let mut chain = vec![first];
        let mut k = start;
        let mut ok = true;
        while chain.last() != Some(&last) {
            let next = (k + step) % n;
            if next == start {
                ok = false;
                break;
            }
            let (a, b) = (walk[k].0, walk[next].0);
            if !runs_on(curved, a, b, tol)? {
                ok = false;
                break;
            }
            chain.push(b);
            k = next;
        }
        if ok && best.as_ref().is_none_or(|c| chain.len() < c.len()) {
            best = Some(chain);
        }
    }
    Ok(best)
}

/// Whether edges `a` and `b` share a vertex where they run on tangent into
/// each other.
fn runs_on(curved: &Curved, a: Id, b: Id, tol: Tolerances) -> OgeomResult<bool> {
    let (ea, eb) = (&curved.edges[&a], &curved.edges[&b]);
    for (va, ta) in [(ea.ends.0, ea.range.0), (ea.ends.1, ea.range.1)] {
        for (vb, tb) in [(eb.ends.0, eb.range.0), (eb.ends.1, eb.range.1)] {
            if va != vb {
                continue;
            }
            let (da, db) = (ea.curve.d1_at(ta, tol)?, eb.curve.d1_at(tb, tol)?);
            if da.cross(db).magnitude()
                <= super::setbacks::TANGENT * da.magnitude() * db.magnitude()
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Whether the chain's `k`th edge runs, along its curve, away from the
/// edge before it (toward the one after): read at the vertex it shares
/// with its neighbour.
fn forward_of(curved: &Curved, chain: &[Id], k: usize) -> bool {
    let e = &curved.edges[&chain[k]];
    if chain.len() == 1 {
        return true;
    }
    if k + 1 < chain.len() {
        let next = &curved.edges[&chain[k + 1]];
        e.ends.1 == next.ends.0 || e.ends.1 == next.ends.1
    } else {
        let prev = &curved.edges[&chain[k - 1]];
        e.ends.0 == prev.ends.0 || e.ends.0 == prev.ends.1
    }
}

/// `curve` over `range`, as a bounded curve: itself where the range is its
/// whole domain.
fn trimmed(curve: &Curve, range: (f64, f64), tol: Tolerances) -> OgeomResult<Curve> {
    let (lo, hi) = curve.domain();
    if (range.0 - lo).abs() <= tol.parametric() && (range.1 - hi).abs() <= tol.parametric() {
        return Ok(curve.clone());
    }
    Ok(ogeom_geom::TrimmedCurve::new(curve.clone(), range.0, range.1, tol)?.into())
}

/// The point of the plane through `o` square to the unit `t` a `radius` from
/// both `face` and the seat `edge` measures to, by Newton from `start`;
/// `None` where it does not converge.
fn centre_in(
    face: &SurfaceGeometry,
    edge: &impl Fn(Point) -> OgeomResult<(usize, ogeom_algo::Projection)>,
    o: Point,
    t: Vector,
    start: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    let u = ogeom_math::Direction::new(t, tol)?
        .any_perpendicular()
        .vector();
    let w = t.cross(u);
    let at = |x: f64, y: f64| o + u * x + w * y;
    let residual = |p: Point| -> OgeomResult<[f64; 2]> {
        let f = ogeom_geom::project::project_on_surface(face, p, 16, tol)?;
        let (_, e) = edge(p)?;
        Ok([f.distance - radius, e.distance - radius])
    };
    let (mut x, mut y) = ((start - o).dot(u), (start - o).dot(w));
    let h = radius * 1e-7;
    for _ in 0..50 {
        let r = residual(at(x, y))?;
        if r[0].abs().max(r[1].abs()) <= tol.confusion() * 1e-3 {
            return Ok(Some(at(x, y)));
        }
        let (rx, ry) = (residual(at(x + h, y))?, residual(at(x, y + h))?);
        let j = [
            [(rx[0] - r[0]) / h, (ry[0] - r[0]) / h],
            [(rx[1] - r[1]) / h, (ry[1] - r[1]) / h],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() <= f64::EPSILON {
            return Ok(None);
        }
        let dx = (r[0] * j[1][1] - r[1] * j[0][1]) / det;
        let dy = (j[0][0] * r[1] - j[1][0] * r[0]) / det;
        x -= dx;
        y -= dy;
        if dx.abs().max(dy.abs()) <= radius * 1e-14 {
            return Ok(Some(at(x, y)));
        }
    }
    Ok(None)
}

/// The ball's arc from `a` to `b` about `centre`, at [`ACROSS`] points.
fn arc_points(
    centre: Point,
    a: Point,
    b: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    let (p, q) = ((a - centre) / radius, (b - centre) / radius);
    let axis = p.cross(q).normalized(tol)?;
    let sweep = p.dot(q).clamp(-1.0, 1.0).acos();
    Ok((0..ACROSS)
        .map(|k| {
            #[expect(clippy::cast_precision_loss, reason = "a few points")]
            let theta = sweep * (k as f64) / ((ACROSS - 1) as f64);
            centre + (p * theta.cos() + axis.cross(p) * theta.sin()) * radius
        })
        .collect())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use crate::engine::curved::read;
    use crate::engine::setbacks::{curved_corners, roll_overs};
    use ogeom_core::FastMap;
    use ogeom_geom::{LineCurve, Surface as _};
    use ogeom_math::{Direction, Frame};
    use ogeom_topo::Model;

    const T: Tolerances = Tolerances::millimetres();

    /// #974's shelf rounded at 3 on the box's top front edge: between the
    /// crossings at y = ±2.8 the ball rests on the cylinder's rim, its
    /// centre (9.6, y, 21 - sqrt(9 - (10 - rho)²)), rho = |(9.6, y)|, its
    /// arc from (12.6, y, z) to the rim's point (96 / rho, 10 y / rho, 21).
    /// The run's surface holds every such arc within 1e-6 and the contact
    /// its feet.
    #[test]
    fn a_shelf_s_ball_rests_on_the_cylinder_s_rim_between_the_crossings() {
        let mut model = Model::new();
        let up = |x: f64, y: f64, z: f64| {
            Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
        };
        let block = ogeom_algo::make_box(&mut model, up(0.0, -5.0, 0.0), (12.6, 10.0, 21.0), T)
            .unwrap()
            .shape;
        let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 26.0, T)
            .unwrap()
            .shape;
        let solid = ogeom_bool::fuse(&mut model, &cylinder, &block, T)
            .unwrap()
            .shape;
        let curved = read(&model, &solid, T).unwrap().unwrap();
        let edge = *curved
            .edges
            .iter()
            .find(|(_, e)| {
                [e.ends.0, e.ends.1].iter().all(|v| {
                    let p = curved.points[v];
                    (p.x - 12.6).abs() < 1e-7 && (p.z - 21.0).abs() < 1e-7
                })
            })
            .unwrap()
            .0;
        let faces = crate::engine::trimmed::faces_of(&curved, edge);
        let line = |a: Point, b: Point| -> Curve { LineCurve::segment(a, b, T).unwrap().into() };
        let top_contact = line(Point::new(9.6, -10.0, 21.0), Point::new(9.6, 10.0, 21.0));
        let front_contact = line(Point::new(12.6, -10.0, 18.0), Point::new(12.6, 10.0, 18.0));
        let is_top = |fi: usize| {
            matches!(&curved.faces[fi].surface, SurfaceGeometry::Plane(p)
                if p.plane().normal().vector().z.abs() > 0.5)
        };
        let (top, front) = if is_top(faces[0]) {
            (faces[0], faces[1])
        } else {
            (faces[1], faces[0])
        };
        let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
        contacts.insert((top, edge), top_contact);
        contacts.insert((front, edge), front_contact);
        let corners = curved_corners(&curved, &contacts, T).unwrap().unwrap();
        let rolls = roll_overs(&curved, &contacts, &corners, T).unwrap();
        assert_eq!(rolls.len(), 2);
        let spine = line(Point::new(9.6, -10.0, 18.0), Point::new(9.6, 10.0, 18.0));
        let run = edge_run(
            &model, &curved, &spine, 3.0, top, front, &rolls[0], &rolls[1], T,
        )
        .unwrap()
        .unwrap();
        assert!(run.off <= 1e-6, "{}", run.off);
        let SurfaceGeometry::BSpline(surface) = &run.surface else {
            panic!("fitted");
        };
        let _ = surface.domain();
        for k in 0..=40 {
            let y = -2.8 + 5.6 * f64::from(k) / 40.0;
            let rho = 9.6_f64.hypot(y);
            let over = 10.0 - rho;
            let centre = Point::new(9.6, y, 21.0 - (9.0 - over * over).sqrt());
            let rim = Point::new(96.0 / rho, 10.0 * y / rho, 21.0);
            let foot = Point::new(12.6, y, centre.z);
            let arc = arc_points(centre, foot, rim, 3.0, T).unwrap();
            for p in arc {
                let d = ogeom_geom::project::project_on_surface(&run.surface, p, 16, T)
                    .unwrap()
                    .distance;
                assert!(d < 1e-6, "y {y}: {p:?} is {d} off the run");
            }
            let on_contact = ogeom_algo::project_on_curve(&run.contact, foot, 64, T).unwrap();
            assert!(
                on_contact.distance < 1e-6,
                "y {y}: foot {}",
                on_contact.distance
            );
        }
    }
}
