//! Rounds of rims between a planar cap and a coaxial wall (a drum's edge, a
//! bore's, a boss's foot), built as one sheet over the solid read by
//! [`curved`](super::curved).
//!
//! A rim is a circle, in one closed edge or in arcs, running the whole turn
//! between a planar cap and a coaxial cylinder or cone, every arc asked at
//! one radius. Its round has no vertex blend. The cap is set back to the
//! ring the ball touches it along and the wall to the ring it touches the
//! wall along, each arc of the rim becoming an arc of each ring over the
//! same angles; between them the strip is the torus the ball sweeps, a face
//! per arc. At each of the rim's vertices the meridian arc from the one
//! ring to the other bounds the torus's faces, and the wall's straight edge
//! through the vertex is cut back to the wall's ring. Every other face,
//! edge and vertex is built again as it was.

use ogeom_algo::{Built, HalfSpaces, History};
use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, LineCurve, SurfaceGeometry, TorusSurface};
use ogeom_math::{Circle, Direction, Frame, Point, Torus};
use ogeom_topo::{Model, Shape};

use super::Id;
use super::curved::{Curved, SheetCurve, SheetRings, build_with_widths, read};

/// One arc of a rim, solved: its cap and wall faces, its arcs on the cap's
/// and the wall's rings, and its vertices from first to last.
struct RimArc {
    edge: Id,
    cap: usize,
    wall: usize,
    on_cap: (Curve, (f64, f64)),
    on_wall: (Curve, (f64, f64)),
    ends: (Id, Id),
}

/// One vertex of a rim, solved: where the ball touches the cap and the wall
/// on its meridian, and the meridian arc from the one to the other.
struct RimVertex {
    vertex: Id,
    on_cap: Point,
    on_wall: Point,
    meridian: (Curve, (f64, f64)),
}

/// One rim, solved: its arcs, its vertices, its torus, and whether the
/// torus's own normal faces into the solid (a round that fills, a boss's
/// foot, rather than cuts).
struct Rim {
    arcs: Vec<RimArc>,
    vertices: Vec<RimVertex>,
    torus: SurfaceGeometry,
    fills: bool,
    material: HalfSpaces,
    /// Outside the box the round can take up: past the cap's plane or the
    /// wall's ring along the axis, or nearer or farther off it than the
    /// rim and both rings.
    clear: HalfSpaces,
}

/// The points each other edge of a rim's cap or wall is read at, to keep
/// the ring clear of it.
const CAP_SAMPLES: u32 = 32;

/// A rim's section: the rolling ball's, or a conic meeting the cap and the
/// wall each its own distance from the rim (an asymmetric round, ADR 0001,
/// amendment "round section laws, conic and elliptical").
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RimSection {
    Ball(f64),
    Conic {
        cap: f64,
        wall: f64,
        profile: crate::asym::RoundProfile,
    },
}

/// Where a rim's round meets its cap and its wall, in the meridian: the cap
/// ring's radius, the wall ring's height along `up` and radius.
#[derive(Debug, Clone, Copy)]
struct Placed {
    rho: f64,
    wall_height: f64,
    wall_rho: f64,
}

/// The rounded solid of `rounds` (each a rim's arc and its radius) on
/// `solid`, or `None` where nothing is asked, `solid` is not one [`read`]
/// takes or the edges asked are not whole rims, each between a planar cap and a coaxial wall,
/// asked at one radius, each of its vertices meeting nothing but its two
/// arcs and one straight edge of the wall, with a ball that fits between
/// its cap and wall.
///
/// # Errors
///
/// As building the sheet and checking the solid.
pub(crate) fn assemble(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let sections: Vec<(Shape, RimSection)> = rounds
        .iter()
        .map(|(e, r)| (e.clone(), RimSection::Ball(*r)))
        .collect();
    assemble_with(model, solid, &sections, tol)
}

/// [`assemble`] with each rim's own [`RimSection`]: a conic round's strip
/// is its meridian conic turned about the axis.
///
/// # Errors
///
/// As [`assemble`].
pub(crate) fn assemble_with(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, RimSection)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(read) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let Some(rims) = rims(model, solid, &read, rounds, tol)? else {
        return Ok(None);
    };
    // No rim, nothing rounded: the sheet would only rebuild the solid as
    // given, which is not a round of anything asked.
    if rims.is_empty() {
        return Ok(None);
    }
    // A round reaching a face it does not blend (a cavity's wall under a
    // thin cap) is not built here: every face but the rim's caps and walls
    // proved outside the box the round can take up, to the bands' fit
    // target. A face whose guaranteed bounding box the region's bound
    // already holds is proved without drawing it: the box's corners span a
    // hull holding the face, and the bound is at least the height over it.
    let boxes: Vec<Vec<Point>> = read
        .faces
        .iter()
        .map(|f| Ok(ogeom_algo::shape_bounds(model, &f.shape, tol)?.corners()))
        .collect::<OgeomResult<_>>()?;
    for rim in &rims {
        let own: Vec<usize> = rim.arcs.iter().flat_map(|a| [a.cap, a.wall]).collect();
        for (fi, face) in read.faces.iter().enumerate() {
            if own.contains(&fi) {
                continue;
            }
            if !boxes[fi].is_empty() && rim.clear.bound_over(&boxes[fi]) <= 0.0 {
                continue;
            }
            if ogeom_algo::certify_face_within(
                model,
                &face.shape,
                &rim.clear,
                crate::marched::band_fit_target(tol),
                tol,
            )?
            .is_err()
            {
                return Ok(None);
            }
        }
    }

    // The sheet's points: every vertex but the rims', and each rim
    // vertex's two contact points.
    let gone: Vec<Id> = rims
        .iter()
        .flat_map(|r| r.vertices.iter().map(|v| v.vertex))
        .collect();
    let mut ids: Vec<Id> = read
        .points
        .keys()
        .copied()
        .filter(|v| !gone.contains(v))
        .collect();
    ids.sort_unstable();
    let index: FastMap<Id, usize> = ids.iter().enumerate().map(|(i, v)| (*v, i)).collect();
    let mut points: Vec<Point> = ids.iter().map(|v| read.points[v]).collect();
    // Each rim vertex's points on the cap and on the wall.
    let mut contact: FastMap<Id, (usize, usize)> = FastMap::default();
    for v in rims.iter().flat_map(|r| &r.vertices) {
        points.push(v.on_cap);
        points.push(v.on_wall);
        contact.insert(v.vertex, (points.len() - 2, points.len() - 1));
    }

    // The sheet's edges: every edge but the rims', a wall's straight edge
    // through a rim vertex cut back to the wall's ring there; then each
    // rim arc's two ring arcs and each rim vertex's meridian.
    let rim_edges: Vec<Id> = rims
        .iter()
        .flat_map(|r| r.arcs.iter().map(|a| a.edge))
        .collect();
    let mut edge_ids: Vec<Id> = read
        .edges
        .keys()
        .copied()
        .filter(|e| !rim_edges.contains(e))
        .collect();
    edge_ids.sort_unstable();
    let mut edges: Vec<SheetCurve> = Vec::new();
    // Each edge carried over as it was keeps the tolerance it stood within.
    let mut widened: Vec<Option<ogeom_core::Tolerance>> = Vec::new();
    let mut edge_index: FastMap<Id, usize> = FastMap::default();
    for e in &edge_ids {
        let edge = &read.edges[e];
        let end = |v: Id| {
            contact
                .get(&v)
                .map(|c| c.1)
                .or_else(|| index.get(&v).copied())
        };
        let (Some(a), Some(b)) = (end(edge.ends.0), end(edge.ends.1)) else {
            return Ok(None);
        };
        let moved = contact.contains_key(&edge.ends.0) || contact.contains_key(&edge.ends.1);
        let (curve, range) = if moved {
            let line = LineCurve::segment(points[a], points[b], tol)?;
            widened.push(None);
            (Curve::Line(line), (0.0, points[a].distance(points[b])))
        } else {
            widened.push(model.tolerance_of(&edge.shape)?);
            (edge.curve.clone(), edge.range)
        };
        edge_index.insert(*e, edges.len());
        edges.push(SheetCurve {
            ends: (a, b),
            curve,
            range,
        });
    }
    // Each rim arc's edges on the cap's ring and on the wall's, by its edge.
    let mut ring_arcs: FastMap<Id, (usize, usize)> = FastMap::default();
    let mut meridians: FastMap<Id, usize> = FastMap::default();
    for rim in &rims {
        for a in &rim.arcs {
            let (c0, w0) = contact[&a.ends.0];
            let (c1, w1) = contact[&a.ends.1];
            edges.push(SheetCurve {
                ends: (c0, c1),
                curve: a.on_cap.0.clone(),
                range: a.on_cap.1,
            });
            edges.push(SheetCurve {
                ends: (w0, w1),
                curve: a.on_wall.0.clone(),
                range: a.on_wall.1,
            });
            ring_arcs.insert(a.edge, (edges.len() - 2, edges.len() - 1));
        }
        // The torus's seam: the meridian at the rim's first vertex.
        for v in rim.vertices.iter().take(1) {
            let (c, w) = contact[&v.vertex];
            edges.push(SheetCurve {
                ends: (c, w),
                curve: v.meridian.0.clone(),
                range: v.meridian.1,
            });
            meridians.insert(v.vertex, edges.len() - 1);
        }
    }

    // The faces: each as it was, a rim arc on its cap swapped for its arc
    // of the cap's ring and on its wall for its arc of the wall's, each
    // walked as the rim arc was.
    let arc_of = |fi: usize, e: Id| -> Option<usize> {
        let a = rims.iter().flat_map(|r| &r.arcs).find(|a| a.edge == e)?;
        let (cap, wall) = ring_arcs[&e];
        if a.cap == fi {
            Some(cap)
        } else if a.wall == fi {
            Some(wall)
        } else {
            None
        }
    };
    let mut faces: Vec<SheetRings> = Vec::with_capacity(read.faces.len());
    for (fi, f) in read.faces.iter().enumerate() {
        let mut rings = Vec::with_capacity(f.loops.len());
        for l in &f.loops {
            let mut ring = Vec::with_capacity(l.len());
            for (e, walk) in l {
                let Some(at) = arc_of(fi, *e).or_else(|| edge_index.get(e).copied()) else {
                    return Ok(None);
                };
                ring.push((at, *walk != f.reversed));
            }
            rings.push(ring);
        }
        faces.push(SheetRings {
            surface: f.surface.clone(),
            rings,
            reversed: f.reversed,
        });
    }
    // Each rim arc's piece of torus: its cap ring arc walked against the
    // cap's walk of the rim arc, on along the meridian at that arc's end to
    // the wall's ring, its wall ring arc against the wall's walk, and back
    // along the meridian at its start.
    // Each rim's torus, one face: its cap ring's arcs walked against the
    // cap's walk of the rim all the way round from the seam, on along the
    // seam to the wall's ring, its wall ring's arcs back round against the
    // wall's walk, and back along the seam.
    let mut torus_of: Vec<Vec<Id>> = Vec::new();
    for rim in &rims {
        let mut on_cap_ring = Vec::with_capacity(rim.arcs.len());
        let mut on_wall_ring = Vec::with_capacity(rim.arcs.len());
        for a in &rim.arcs {
            let walked_on = |fi: usize| -> Option<bool> {
                read.faces[fi]
                    .loops
                    .iter()
                    .flatten()
                    .find(|(e, _)| *e == a.edge)
                    .map(|(_, w)| *w)
            };
            let (Some(on_cap), Some(on_wall)) = (walked_on(a.cap), walked_on(a.wall)) else {
                return Ok(None);
            };
            if on_cap == on_wall {
                return Ok(None);
            }
            let (cap, wall) = ring_arcs[&a.edge];
            on_cap_ring.push((cap, !on_cap));
            on_wall_ring.push((wall, !on_wall));
        }
        let Some(first) = rim.vertices.first() else {
            return Ok(None);
        };
        let seam = meridians[&first.vertex];
        let (c0, w0) = edges[seam].ends;
        let (Some(cap_round), Some(wall_round)) = (
            chained(&edges, &on_cap_ring, c0),
            chained(&edges, &on_wall_ring, w0),
        ) else {
            return Ok(None);
        };
        let mut standing = cap_round;
        standing.push((seam, true));
        standing.extend(wall_round);
        standing.push((seam, false));
        let ring: Vec<(usize, bool)> = if rim.fills {
            standing.iter().rev().map(|(e, w)| (*e, !w)).collect()
        } else {
            standing
        };
        faces.push(SheetRings {
            surface: rim.torus.clone(),
            rings: vec![ring],
            reversed: rim.fills,
        });
        torus_of.push(rim.arcs.iter().map(|a| a.edge).collect());
    }

    let widths: Vec<f64> = widened
        .iter()
        .map(|w| w.map_or(0.0, ogeom_core::Tolerance::get))
        .collect();
    let sheet = build_with_widths(model, &points, &edges, &faces, &widths, tol)?;
    let shell = ogeom_algo::make_shell(model, &sheet.faces)?.shape;
    let rounded = ogeom_algo::make_solid(model, &[shell])?.shape;
    // A wall's or a torus's seam, walked twice, as a seam's pair of pcurves.
    ogeom_algo::join_seam_columns(model, &rounded, tol)?;
    super::guards::solid_fits(model, &rounded, tol)?;
    // P5: each piece of torus proved within its rim's material, to the
    // bands' fit target.
    let fit = crate::marched::band_fit_target(tol);
    let materials = rims.iter().map(|r| &r.material);
    for (face, region) in sheet.faces[read.faces.len()..].iter().zip(materials) {
        if let Err(out) = ogeom_algo::certify_face_within(model, face, region, fit, tol)? {
            ogeom_core::ogeom_bail!(
                Construction,
                "a rim's torus stands {:.1e} past the material it blends, at {:?}{}",
                out.past,
                out.at,
                if out.found { "" } else { " (not cleared)" }
            );
        }
    }

    let mut history = History::new();
    for (f, b) in read.faces.iter().zip(&sheet.faces) {
        history.modify(&f.shape, b.clone());
    }
    for (arcs, torus) in torus_of.iter().zip(&sheet.faces[read.faces.len()..]) {
        for e in arcs {
            history.delete(&read.shapes[e]);
            history.generate(&read.shapes[e], torus.clone());
        }
    }
    for v in &gone {
        history.delete(&read.shapes[v]);
    }
    for e in &edge_ids {
        history.modify(&read.edges[e].shape, sheet.edges[edge_index[e]].clone());
    }
    for (i, v) in ids.iter().enumerate() {
        history.modify(&read.shapes[v], sheet.vertices[i].clone());
    }
    Ok(Some(Built::new(rounded, history)))
}

/// A rim as asked: its seat, its section and its arcs with their circles.
type RimAsked = (crate::support::RevolvedSeat, RimSection, Vec<(Id, Circle)>);

/// The rims `rounds` asks, each solved, or `None` where they are not whole
/// rims of the class [`assemble`] builds.
fn rims(
    model: &Model,
    solid: &Shape,
    read: &Curved,
    rounds: &[(Shape, RimSection)],
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Rim>>> {
    // Each arc asked with its seat, grouped by its circle.
    let mut groups: Vec<RimAsked> = Vec::new();
    for (edge, radius) in rounds {
        let e = super::id(edge);
        let Some(arc) = read.edges.get(&e) else {
            return Ok(None);
        };
        let Curve::Circle(circle) = &arc.curve else {
            return Ok(None);
        };
        let Ok(seat) = crate::support::revolved_seat(model, solid, &arc.shape, circle, tol) else {
            return Ok(None);
        };
        let c = circle.circle();
        let same = |g: &Circle| {
            g.centre().distance(c.centre()) <= tol.confusion()
                && (g.radius() - c.radius()).abs() <= tol.confusion()
        };
        match groups.iter_mut().find(|(_, _, arcs)| same(&arcs[0].1)) {
            Some((_, r, arcs)) => {
                if *r != *radius || arcs.iter().any(|(x, _)| *x == e) {
                    return Ok(None);
                }
                arcs.push((e, c));
            }
            None => groups.push((seat, radius.clone(), vec![(e, c)])),
        }
    }
    // Built in the solid's order, not the request's: each rim's arcs by
    // their edges, the rims by their first.
    for (_, _, arcs) in &mut groups {
        arcs.sort_by_key(|(e, _)| *e);
    }
    groups.sort_by_key(|(_, _, arcs)| arcs[0].0);
    let face_of = |s: &Shape| read.faces.iter().position(|f| f.shape.is_same(s));
    let tau = core::f64::consts::TAU;
    let mut out = Vec::with_capacity(groups.len());
    for (seat, radius, arcs) in groups {
        // The arcs run the whole turn.
        let turn: f64 = arcs
            .iter()
            .map(|(e, _)| read.edges[e].range.1 - read.edges[e].range.0)
            .sum();
        if (turn - tau).abs() > tol.angular() {
            return Ok(None);
        }
        // Where the round meets the cap and the wall: the ball's rings, or
        // a conic's at its two distances (the cap's along it, the wall's
        // along its generator, both away from the rim into the material).
        let (ball, weight) = match &radius {
            RimSection::Ball(r) => {
                let b = seat.ball(*r);
                (
                    Placed {
                        rho: b.rho,
                        wall_height: b.wall_height,
                        wall_rho: b.wall_rho,
                    },
                    None,
                )
            }
            RimSection::Conic {
                cap: d_cap,
                wall: d_wall,
                profile,
            } => {
                let (d_cap, d_wall) = (*d_cap, *d_wall);
                let w = seat.slope.hypot(1.0);
                let wall_height = -seat.tau * d_wall / w;
                (
                    Placed {
                        rho: seat.radius - seat.sigma * seat.tau * d_cap,
                        wall_height,
                        wall_rho: seat.radius_at(wall_height),
                    },
                    Some(*profile),
                )
            }
        };
        // The round would touch its cap at or past the axis: no ring is
        // left on the cap for it to stand on, the ball as wide as the cap
        // or wider. Refused by name, as mainstream refuses a round that
        // consumes the face it rounds onto.
        if ball.rho <= tol.confusion() || ball.wall_rho <= tol.confusion() {
            ogeom_core::ogeom_bail!(
                Construction,
                "a rim's round meets its cap at or past the axis: the round is as \
                 wide as the cap it rounds onto, or wider; see docs/CAPABILITIES.md, \
                 fillet.edge-blends"
            );
        }
        let mut solved_arcs = Vec::with_capacity(arcs.len());
        let mut vertex_ids: Vec<Id> = Vec::new();
        for (e, circle) in &arcs {
            let edge = &read.edges[e];
            let Ok(arc_seat) = crate::support::revolved_seat(
                model,
                solid,
                &edge.shape,
                &CircleCurve::new(*circle),
                tol,
            ) else {
                return Ok(None);
            };
            let (Some(cap), Some(wall)) =
                (face_of(&arc_seat.cap_face), face_of(&arc_seat.wall_face))
            else {
                return Ok(None);
            };
            // Each ring's arc over the rim arc's own angles, its frame the
            // rim's moved along the axis.
            let frame = circle.frame();
            let ring = |level: f64, r: f64| -> OgeomResult<Curve> {
                let origin = seat.centre + seat.up * level;
                Ok(Curve::Circle(CircleCurve::new(Circle::new(
                    Frame::new(origin, frame.z(), frame.x(), tol)?,
                    r,
                    tol,
                )?)))
            };
            solved_arcs.push(RimArc {
                edge: *e,
                cap,
                wall,
                on_cap: (ring(0.0, ball.rho)?, edge.range),
                on_wall: (ring(ball.wall_height, ball.wall_rho)?, edge.range),
                ends: edge.ends,
            });
            for v in [edge.ends.0, edge.ends.1] {
                if !vertex_ids.contains(&v) {
                    vertex_ids.push(v);
                }
            }
        }
        let rim_edges: Vec<Id> = arcs.iter().map(|(e, _)| *e).collect();
        // The wall's ring stands on every wall face, short of its far side:
        // every other edge of the wall farther along the axis from the rim
        // than the ring, read along its length, or at its far end for a
        // straight edge from a rim vertex, which is cut back to the ring.
        for a in &solved_arcs {
            for (x, _) in read.faces[a.wall].loops.iter().flatten() {
                if rim_edges.contains(x) {
                    continue;
                }
                let edge = &read.edges[x];
                let from_rim = [edge.ends.0, edge.ends.1]
                    .iter()
                    .any(|w| vertex_ids.contains(w));
                let (t0, t1) = edge.range;
                let samples: Vec<Point> = if from_rim {
                    [edge.ends.0, edge.ends.1]
                        .iter()
                        .filter(|w| !vertex_ids.contains(w))
                        .map(|w| read.points[w])
                        .collect()
                } else {
                    (0..=CAP_SAMPLES)
                        .map(|k| {
                            let t = t0 + (t1 - t0) * f64::from(k) / f64::from(CAP_SAMPLES);
                            ogeom_geom::Curve3d::point_at(&edge.curve, t, tol)
                        })
                        .collect::<OgeomResult<_>>()?
                };
                for p in samples {
                    let along = (p - seat.centre).dot(seat.up).abs();
                    if ball.wall_height.abs() >= along - tol.confusion() {
                        return Ok(None);
                    }
                }
            }
        }
        // The cap's ring stands on every cap face, nothing else of the cap
        // between it and the rim: every other edge of the cap, read along
        // its length, stays farther from the rim than the ring. A cap
        // narrower than the setback is left to the current path.
        let inward = ball.rho < seat.radius;
        for a in &solved_arcs {
            for (x, _) in read.faces[a.cap].loops.iter().flatten() {
                if rim_edges.contains(x) {
                    continue;
                }
                let edge = &read.edges[x];
                let (t0, t1) = edge.range;
                for k in 0..=CAP_SAMPLES {
                    let t = t0 + (t1 - t0) * f64::from(k) / f64::from(CAP_SAMPLES);
                    let off = ogeom_geom::Curve3d::point_at(&edge.curve, t, tol)? - seat.centre;
                    let rho = (off - seat.up * off.dot(seat.up)).magnitude();
                    let clear = if inward {
                        rho < ball.rho - tol.confusion()
                    } else {
                        rho > ball.rho + tol.confusion()
                    };
                    if !clear {
                        return Ok(None);
                    }
                }
            }
        }
        let mut vertices = Vec::with_capacity(vertex_ids.len());
        for v in vertex_ids {
            // The vertex meets its rim's arcs and at most one more edge, a
            // straight one of the wall.
            let others: Vec<Id> = read
                .edges
                .iter()
                .filter(|(id, x)| !rim_edges.contains(id) && (x.ends.0 == v || x.ends.1 == v))
                .map(|(id, _)| *id)
                .collect();
            // Beside a straight wall edge, or at a mere split of the rim.
            match others.as_slice() {
                [] => {}
                [wall_edge] if matches!(read.edges[wall_edge].curve, Curve::Line(_)) => {}
                _ => return Ok(None),
            }
            let off = read.points[&v] - seat.centre;
            let radial = (off - seat.up * off.dot(seat.up)).normalized(tol)?;
            let meridian = meridian_at(&seat, &radius, ball, weight, radial, tol)?;
            vertices.push(RimVertex {
                vertex: v,
                on_cap: meridian.0,
                on_wall: meridian.1,
                meridian: meridian.2,
            });
        }
        let torus: SurfaceGeometry = match &radius {
            RimSection::Ball(r) => TorusSurface::new(Torus::new(
                seat.frame_at(seat.centre + seat.up * seat.ball(*r).height, tol)?,
                ball.rho,
                *r,
                tol,
            )?)
            .into(),
            RimSection::Conic { .. } => {
                // The meridian conic at the rim frame's own `x`, turned the
                // whole way about the axis.
                let frame = seat.frame_at(seat.centre, tol)?;
                let (_, _, (curve, _)) =
                    meridian_at(&seat, &radius, ball, weight, frame.x().vector(), tol)?;
                // Run from the wall to the cap, as the torus's own tube
                // turns, so the surface faces the way the torus would.
                let curve = ogeom_geom::Reversible::reversed(&curve);
                ogeom_geom::RevolutionSurface::new(
                    curve,
                    ogeom_math::Axis {
                        location: seat.centre,
                        direction: frame.z(),
                    },
                    core::f64::consts::TAU,
                )?
                .into()
            }
        };
        // The material round the rim once rounded: below the cap and on the
        // wall's side of it at a convex rim; at a concave one below the cap,
        // or on the wall's side, or in the fill, the wedge between the cap
        // raised to the wall's contact ring and the wall moved out to the
        // cap's.
        let below = |rise: f64| HalfSpaces::Below {
            normal: seat.up,
            at: seat.centre + seat.up * rise,
        };
        let wall = |radius: f64| HalfSpaces::Revolved {
            at: seat.centre,
            axis: seat.up,
            radius,
            slope: seat.slope,
            inside: seat.sigma > 0.0,
        };
        let material = if seat.additive() {
            HalfSpaces::Any(vec![
                below(0.0),
                wall(seat.radius),
                HalfSpaces::Every(vec![below(ball.wall_height), wall(ball.rho)]),
            ])
        } else {
            HalfSpaces::Every(vec![below(0.0), wall(seat.radius)])
        };
        let (h_lo, h_hi) = (ball.wall_height.min(0.0), ball.wall_height.max(0.0));
        let radii = [seat.radius, ball.rho, ball.wall_rho];
        let (r_lo, r_hi) = (
            radii.iter().copied().fold(f64::INFINITY, f64::min),
            radii.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        );
        let off_axis = |radius: f64, inside: bool| HalfSpaces::Revolved {
            at: seat.centre,
            axis: seat.up,
            radius,
            slope: 0.0,
            inside,
        };
        let clear = HalfSpaces::Any(vec![
            HalfSpaces::Below {
                normal: -seat.up,
                at: seat.centre + seat.up * h_hi,
            },
            below(h_lo),
            off_axis(r_lo, true),
            off_axis(r_hi, false),
        ]);
        out.push(Rim {
            arcs: solved_arcs,
            vertices,
            torus,
            fills: seat.additive(),
            material,
            clear,
        });
    }
    Ok(Some(out))
}

/// The oriented `pieces` (each an edge and whether it is walked from its
/// first end) chained into one loop from `start` back to it, or `None`
/// where they do not close.
fn chained(
    edges: &[SheetCurve],
    pieces: &[(usize, bool)],
    start: usize,
) -> Option<Vec<(usize, bool)>> {
    let ends = |(e, forward): (usize, bool)| {
        let (a, b) = edges[e].ends;
        if forward { (a, b) } else { (b, a) }
    };
    let mut left: Vec<(usize, bool)> = pieces.to_vec();
    let mut at = start;
    let mut out = Vec::with_capacity(pieces.len());
    while !left.is_empty() {
        let k = left.iter().position(|p| ends(*p).0 == at)?;
        let p = left.remove(k);
        at = ends(p).1;
        out.push(p);
    }
    (at == start).then_some(out)
}

/// A meridian as [`meridian_at`] reads it: the touches on the cap and the
/// wall, and the section's curve between them with its range.
type Meridian = (Point, Point, (Curve, (f64, f64)));

/// A rim's meridian at the radial direction `radial`: where its round
/// meets the cap and the wall there, and the section between them (the
/// ball's arc, or the conic whose control point is the rim's own point,
/// tangent to the cap and the wall at its contacts, of `weight`'s
/// profile).
fn meridian_at(
    seat: &crate::support::RevolvedSeat,
    section: &RimSection,
    at: Placed,
    weight: Option<crate::asym::RoundProfile>,
    radial: ogeom_math::Vector,
    tol: Tolerances,
) -> OgeomResult<Meridian> {
    let on_cap = seat.centre + radial * at.rho;
    let on_wall = seat.centre + seat.up * at.wall_height + radial * at.wall_rho;
    match (section, weight) {
        (RimSection::Ball(r), _) => {
            let ball = seat.ball(*r);
            let centre = seat.centre + seat.up * ball.height + radial * ball.rho;
            let from = (on_cap - centre).normalized(tol)?;
            let to = on_wall - centre;
            let normal = Direction::new(from.cross(to), tol)?;
            let circle = Circle::new(
                Frame::new(centre, normal, Direction::new(from, tol)?, tol)?,
                *r,
                tol,
            )?;
            let y = normal.vector().cross(from);
            Ok((
                on_cap,
                on_wall,
                (
                    Curve::Circle(CircleCurve::new(circle)),
                    (0.0, to.dot(y).atan2(to.dot(from))),
                ),
            ))
        }
        (RimSection::Conic { .. }, profile) => {
            let rim = seat.centre + radial * seat.radius;
            let (a, b) = (
                (on_cap - rim).normalized(tol)?,
                (on_wall - rim).normalized(tol)?,
            );
            let w = profile
                .unwrap_or(crate::asym::RoundProfile::Elliptical)
                .weight(a.dot(b))?;
            // The conic in quarters, each exact (its homogeneous control
            // polygon halved twice): a strongly weighted conic is one span
            // too few for the measures that integrate along it.
            let mut pieces: Vec<[(ogeom_math::Vector, f64); 3]> = vec![[
                (on_cap - Point::ORIGIN, 1.0),
                ((rim - Point::ORIGIN) * w, w),
                (on_wall - Point::ORIGIN, 1.0),
            ]];
            for _ in 0..2 {
                pieces = pieces
                    .iter()
                    .flat_map(|[a, b, c]| {
                        let mid = |p: &(ogeom_math::Vector, f64), q: &(ogeom_math::Vector, f64)| {
                            ((p.0 + q.0) * 0.5, 0.5 * (p.1 + q.1))
                        };
                        let (ab, bc) = (mid(a, b), mid(b, c));
                        let abc = mid(&ab, &bc);
                        [[*a, ab, abc], [abc, bc, *c]]
                    })
                    .collect();
            }
            let n = pieces.len();
            let mut knots = vec![0.0; 3];
            #[expect(clippy::cast_precision_loss, reason = "four pieces")]
            for k in 1..n {
                let at = k as f64 / n as f64;
                knots.extend_from_slice(&[at, at]);
            }
            knots.extend_from_slice(&[1.0; 3]);
            let knots = ogeom_math::KnotVector::new(knots, 2)?;
            let mut control = Vec::with_capacity(2 * n + 1);
            for (k, piece) in pieces.iter().enumerate() {
                for (h, wt) in &piece[usize::from(k > 0)..] {
                    control.push(ogeom_math::Weighted::new(
                        Point::ORIGIN + *h / *wt,
                        *wt,
                        tol,
                    )?);
                }
            }
            Ok((
                on_cap,
                on_wall,
                (
                    Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?),
                    (0.0, 1.0),
                ),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::assemble;
    use core::f64::consts::PI;
    use ogeom_core::Tolerances;
    use ogeom_math::{Direction, Frame, Point, Vector};
    use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn volume(m: &Model, s: &Shape) -> f64 {
        ogeom_algo::volume_properties(
            m,
            s,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass
    }

    /// The area a round of radius `r` takes off a square corner or adds to
    /// one, and how far its centroid stands from the corner along each leg.
    fn corner(r: f64) -> (f64, f64) {
        (
            (1.0 - PI / 4.0) * r * r,
            r * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI),
        )
    }

    fn up_at(x: f64, y: f64, z: f64) -> Frame {
        Frame::new(
            Point::new(x, y, z),
            Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
            Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
            T,
        )
        .unwrap()
    }

    /// The circle edges of `s` at height `z`.
    fn rims_at(m: &Model, s: &Shape, z: f64) -> Vec<Shape> {
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (curve, _) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                matches!(curve, ogeom_geom::Curve::Circle(c) if (c.circle().centre().z - z).abs() < 1e-9)
            })
            .collect()
    }

    fn valid(m: &Model, s: &Shape) {
        let found = ogeom_algo::check(m, s, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
    }

    /// A drum of radius 2, 3 tall, its top rim rounded at 0.5: the drum
    /// less the corner's area swept round its centroid, Pappus' ring at
    /// `2 - x`.
    #[test]
    fn a_drum_s_rim_is_its_torus() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let rim = rims_at(&m, &drum, 3.0);
        assert_eq!(rim.len(), 1);
        let built = assemble(&mut m, &drum, &[(rim[0].clone(), 0.5)], T)
            .unwrap()
            .unwrap();
        valid(&m, &built.shape);
        let (area, x) = corner(0.5);
        let want = PI * 4.0 * 3.0 - area * 2.0 * PI * (2.0 - x);
        let v = volume(&m, &built.shape);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A round reported built changes the solid: the drum's rim takes off
    /// its torus's spandrel, and a request of nothing rounds nothing and is
    /// not built (the sheet rebuilt the drum as given, reported as the
    /// round). Each request through the rims alone and through the
    /// engine's whole dispatch.
    #[test]
    fn a_round_reported_built_changes_the_solid() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let rim = rims_at(&m, &drum, 3.0);
        let before = volume(&m, &drum);
        let requests: [Vec<(Shape, f64)>; 2] = [Vec::new(), vec![(rim[0].clone(), 0.5)]];
        for asked in &requests {
            let by_rims = assemble(&mut m, &drum, asked, T).unwrap();
            let by_engine = crate::fillet::rounds_by_engine(&mut m, &drum, asked, T).unwrap();
            for built in [by_rims, by_engine].into_iter().flatten() {
                let after = volume(&m, &built.shape);
                assert!(
                    (before - after).abs() > 1e-6 * before,
                    "{} edges asked: reported built, yet {after} against {before}",
                    asked.len()
                );
            }
        }
        assert!(assemble(&mut m, &drum, &requests[1], T).unwrap().is_some());
    }

    /// The drum's wall and the rim's torus each walk their seam twice: each
    /// seam carries its pair of pcurves, one per walk, so a later walk of
    /// the face in its chart closes; nothing is left for
    /// `join_seam_columns` to join.
    #[test]
    fn a_rounded_drum_s_seams_carry_their_pairs_of_pcurves() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let top = rims_at(&m, &drum, 3.0);
        let built = assemble(&mut m, &drum, &[(top[0].clone(), 0.5)], T)
            .unwrap()
            .unwrap()
            .shape;
        assert_eq!(ogeom_algo::join_seam_columns(&mut m, &built, T).unwrap(), 0);
    }

    /// Both of the drum's rims at once: each takes its own ring.
    #[test]
    fn both_rims_of_a_drum_round_at_once() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let mut asked: Vec<(Shape, f64)> = rims_at(&m, &drum, 3.0)
            .into_iter()
            .map(|e| (e, 0.5))
            .collect();
        asked.extend(rims_at(&m, &drum, 0.0).into_iter().map(|e| (e, 0.75)));
        let built = assemble(&mut m, &drum, &asked, T).unwrap().unwrap();
        valid(&m, &built.shape);
        let (a1, x1) = corner(0.5);
        let (a2, x2) = corner(0.75);
        let want = PI * 4.0 * 3.0 - a1 * 2.0 * PI * (2.0 - x1) - a2 * 2.0 * PI * (2.0 - x2);
        let v = volume(&m, &built.shape);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A block with a bore of radius 1.5 through it, the bore's top rim
    /// rounded at 0.5: the corner's area swept round outside the bore, at
    /// `1.5 + x`.
    #[test]
    fn a_bore_s_rim_rounds_outward() {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 4.0), T)
            .unwrap()
            .shape;
        let bore = ogeom_algo::make_cylinder(&mut m, up_at(5.0, 4.0, -1.0), 1.5, 6.0, T)
            .unwrap()
            .shape;
        let holed = ogeom_bool::cut(&mut m, &block, &bore, T).unwrap().shape;
        let before = volume(&m, &holed);
        let rim: Vec<(Shape, f64)> = rims_at(&m, &holed, 4.0)
            .into_iter()
            .map(|e| (e, 0.5))
            .collect();
        assert!(!rim.is_empty());
        let built = assemble(&mut m, &holed, &rim, T).unwrap().unwrap();
        valid(&m, &built.shape);
        let (area, x) = corner(0.5);
        let want = before - area * 2.0 * PI * (1.5 + x);
        let v = volume(&m, &built.shape);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A boss of radius 1.5 on a block, its foot rounded at 0.5: the round
    /// fills the concave corner, the corner's area swept round outside the
    /// boss, at `1.5 + x`.
    #[test]
    fn a_boss_s_foot_rounds_by_filling() {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 4.0), T)
            .unwrap()
            .shape;
        let boss = ogeom_algo::make_cylinder(&mut m, up_at(5.0, 4.0, 4.0), 1.5, 2.0, T)
            .unwrap()
            .shape;
        let bossed = ogeom_bool::fuse(&mut m, &block, &boss, T).unwrap().shape;
        let before = volume(&m, &bossed);
        let foot: Vec<(Shape, f64)> = rims_at(&m, &bossed, 4.0)
            .into_iter()
            .map(|e| (e, 0.5))
            .collect();
        assert!(!foot.is_empty());
        let built = assemble(&mut m, &bossed, &foot, T).unwrap().unwrap();
        valid(&m, &built.shape);
        let (area, x) = corner(0.5);
        let want = before + area * 2.0 * PI * (1.5 + x);
        let v = volume(&m, &built.shape);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A cone 4 tall narrowing from radius 3 to 1, its top rim rounded at
    /// 0.4: the cap and wall meet at the cone's angle, not square. The
    /// round takes the corner's section in the meridian (the kite from the
    /// corner to the ball's two contacts and its centre, less the ball's
    /// sector between the contacts) turned about the axis: Pappus' volume
    /// from the section's area and its centroid's radius, each in closed
    /// form.
    #[test]
    fn a_cone_s_rim_is_its_torus_by_pappus() {
        let mut m = Model::new();
        let cone = ogeom_algo::make_cone(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 1.0, 4.0, T)
            .unwrap()
            .shape;
        let rim = rims_at(&m, &cone, 4.0);
        let asked: Vec<(Shape, f64)> = rim.iter().map(|e| (e.clone(), 0.4)).collect();
        let built = assemble(&mut m, &cone, &asked, T).unwrap().unwrap();
        valid(&m, &built.shape);
        // In the meridian (radius, height): the corner, and the unit ways
        // into the material along the cap and down the wall.
        let r = 0.4;
        let corner = (1.0_f64, 4.0_f64);
        let along_cap = (-1.0_f64, 0.0_f64);
        let down_wall = (2.0 / 20.0_f64.sqrt(), -4.0 / 20.0_f64.sqrt());
        let angle = (along_cap.0 * down_wall.0 + along_cap.1 * down_wall.1).acos();
        let half = 0.5 * angle;
        let reach = r / half.tan();
        let bisector = {
            let (x, y) = (along_cap.0 + down_wall.0, along_cap.1 + down_wall.1);
            let n = x.hypot(y);
            (x / n, y / n)
        };
        let at = |(x, y): (f64, f64), k: f64| (corner.0 + x * k, corner.1 + y * k);
        let (cap_touch, wall_touch) = (at(along_cap, reach), at(down_wall, reach));
        let centre = at(bisector, r / half.sin());
        // Each right triangle corner-touch-centre, and the sector the ball
        // spans between its touches, its centroid toward the corner.
        let triangle = 0.5 * reach * r;
        let rho = |p: (f64, f64), q: (f64, f64)| (corner.0 + p.0 + q.0) / 3.0;
        let sweep = PI - angle;
        let sector = 0.5 * r * r * sweep;
        let sector_rho = centre.0 - bisector.0 * (4.0 * r * (0.5 * sweep).sin() / (3.0 * sweep));
        let area = 2.0 * triangle - sector;
        let centroid = (triangle * rho(cap_touch, centre) + triangle * rho(wall_touch, centre)
            - sector * sector_rho)
            / area;
        let want = volume(&m, &cone) - 2.0 * PI * centroid * area;
        let v = volume(&m, &built.shape);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// Both rims of a drum asked in either order: the same solid to the bit.
    #[test]
    fn the_rims_order_does_not_change_the_drum() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let mut asked: Vec<(Shape, f64)> = rims_at(&m, &drum, 3.0)
            .into_iter()
            .map(|e| (e, 0.5))
            .collect();
        asked.extend(rims_at(&m, &drum, 0.0).into_iter().map(|e| (e, 0.75)));
        let bits = |m: &Model, s: &Shape| {
            let mut points: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
                .unwrap()
                .iter()
                .map(|v| {
                    let p = m
                        .node(v)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .unwrap();
                    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
                })
                .collect();
            points.sort_unstable();
            (volume(m, s).to_bits(), points)
        };
        let forward = assemble(&mut m, &drum, &asked, T).unwrap().unwrap().shape;
        let reversed: Vec<(Shape, f64)> = asked.iter().rev().cloned().collect();
        let backward = assemble(&mut m, &drum, &reversed, T)
            .unwrap()
            .unwrap()
            .shape;
        assert_eq!(bits(&m, &forward), bits(&m, &backward));
    }

    /// A rim's round comes through the history as the torus it generates,
    /// the rim deleted; every other face and edge as one of the result.
    #[test]
    fn a_rim_s_round_is_its_torus_in_the_history() {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 3.0, T)
            .unwrap()
            .shape;
        let rim = rims_at(&m, &drum, 3.0);
        let built = assemble(&mut m, &drum, &[(rim[0].clone(), 0.5)], T)
            .unwrap()
            .unwrap();
        assert!(built.history.is_deleted(&rim[0]));
        let made = built.history.generated(&rim[0]);
        assert_eq!(made.len(), 1);
        let data = m.node(&made[0]).unwrap().data().as_face().unwrap();
        assert!(matches!(
            m.geometry().surface(data.surface).unwrap(),
            ogeom_geom::SurfaceGeometry::Torus(_)
        ));
        let faces = explore_unique(&m, &built.shape, ShapeType::Face).unwrap();
        for f in explore_unique(&m, &drum, ShapeType::Face).unwrap() {
            let images = built.history.trace(&f);
            assert_eq!(images.len(), 1);
            assert!(faces.iter().any(|g| g.is_same(&images[0])));
        }
    }

    /// A conical tube's four rims (a frustum bored through: two caps each
    /// walking two rims, the bore's cut in arcs by the boolean) asked in
    /// either order: the same solid to the bit, its rims built in the
    /// solid's order rather than the request's.
    #[test]
    fn a_conical_tube_s_rims_build_alike_whatever_the_order() {
        let mut m = Model::new();
        let outer = ogeom_algo::make_cone(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 2.0, 4.0, T)
            .unwrap()
            .shape;
        let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 1.0, 6.0, T)
            .unwrap()
            .shape;
        let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
        let mut asked: Vec<(Shape, f64)> = rims_at(&m, &tube, 4.0)
            .into_iter()
            .map(|e| (e, 0.3))
            .collect();
        asked.extend(rims_at(&m, &tube, 0.0).into_iter().map(|e| (e, 0.3)));
        let bits = |m: &Model, s: &Shape| {
            let mut points: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
                .unwrap()
                .iter()
                .map(|v| {
                    let p = m
                        .node(v)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .unwrap();
                    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
                })
                .collect();
            points.sort_unstable();
            (volume(m, s).to_bits(), points)
        };
        let forward = assemble(&mut m, &tube, &asked, T).unwrap().unwrap().shape;
        valid(&m, &forward);
        let reversed: Vec<(Shape, f64)> = asked.iter().rev().cloned().collect();
        let backward = assemble(&mut m, &tube, &reversed, T)
            .unwrap()
            .unwrap()
            .shape;
        assert_eq!(bits(&m, &forward), bits(&m, &backward));
    }
}
