//! Constant-radius fillets, where the rolling ball's envelope is closed-form.
//!
//! A ball of constant radius rolling along a straight edge between two planar
//! faces stays tangent to both, and its envelope is exactly a cylinder; the
//! same ball rolling around the rim where a cylindrical wall meets its
//! perpendicular cap traces a torus. These are the blends whose *new* surface
//! costs nothing, and both are built the same way: the chamfer's wedge with
//! the bevel exchanged for the envelope (legs running along the faces to the
//! tangency lines, and the envelope between them, every pcurve exact). The
//! boolean then does what it does for the chamfer: the legs melt into the
//! solid's own faces by same-domain resolution, and the envelope stays as
//! the blend.

use crate::support::{
    Seat, apply_wedge, apply_wedge_solid, edge_curve, face_from_edges, planar_face, planar_seat,
    revolved_flanks, revolved_seat, segment_between, trimmed_wedge, wedge_solid,
};
use ogeom_algo::{Built, make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Torus, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// Round a convex edge of a solid with a constant-radius blend tangent to
/// both of its faces.
///
/// Two seats are recognised, and each has a closed-form envelope: a straight
/// edge between two planar faces, where the rolling ball traces a cylinder,
/// and a circular rim where a cylindrical wall meets a perpendicular planar
/// cap, where it traces a torus. In both, the result is the boolean
/// difference with a wedge whose curved face is that envelope, so the history
/// reads as a cut: the faces are modified into their trimmed pieces, and the
/// edge's neighbourhood gains the blend face.
///
/// The blend is the one [`fillet_edges`] builds on this one edge. Where a
/// face cuts the corner away past an end (a floor giving way to a ramp
/// under a wall that carries on, or a boss's cylinder past the end of a
/// straight edge beside it), both run it out onto that face until its whole
/// section has left the material. Where a
/// host face narrows along the edge and the ball's contact runs off it
/// out of the solid, both cut the band off on the face that follows. One
/// thing differs: where the edge ends against a neighbouring blend, this
/// blend stops flush there instead of running on through it. Three such
/// flush blends and their vertex are what
/// [`round_vertex`](crate::round_vertex) closes.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`] if the edge is
/// neither of the seats above, is concave, or `radius` is not a usable
/// length, or the ball sets back past the far side of a face into the
/// material, or past it all along the edge.
pub fn fillet_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| fillet_edge_own(model, solid, &parts[0], radius, tol),
    )
}

/// [`fillet_edge`] on a shape every node below which is its own.
fn fillet_edge_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    // The ball touches each face a setback back from the edge; set back
    // past a face's far side into the material, or all along the edge,
    // the band would cut through the face.
    crate::support::ball_fits(model, solid, edge, radius, true, tol)?;
    fillet_edge_meeting(model, solid, edge, radius, None, tol)
}

/// An edge of a multi-edge request, by its ends: where each end stands and
/// which way the edge leaves it, so a later edge can tell a chain mate from
/// a corner mate at a shared vertex.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Mate {
    /// `(point, unit tangent leaving the point along the edge)` per end.
    pub ends: [(Point, Vector); 2],
    /// Per end, whether the corner tool has already closed the vertex there:
    /// the band stops flush against the corner's patch rather than running
    /// out through it.
    pub settled: [bool; 2],
    /// Per end, whether a junction of the request meets the band there in
    /// a plane: the band stops at the vertex and the junction closes the
    /// corner.
    pub meets: [bool; 2],
    /// Per end, whether a crossed junction of the request across a convex
    /// turn stands there: the band, built from the solid as given, runs out
    /// there as it does alone.
    pub alone: [bool; 2],
}

impl Mate {
    /// Whether any mate's end at `at` is a vertex the corner tool closed,
    /// or `at` is an end of the request's own edge `index` whose vertex the
    /// corner tool closed.
    ///
    /// The second reading is the one a band asks of itself: a corner tool
    /// trims the edge back to the plane square to it through the ball's
    /// centre, so the edge's end now stands a setback from the vertex it
    /// was named at, and the nearer of its original ends says which one
    /// it is.
    pub(crate) fn settled_at(mates: &[Self], index: usize, at: Point, tol: Tolerances) -> bool {
        let other = mates.iter().any(|mate| {
            mate.ends
                .iter()
                .zip(mate.settled)
                .any(|((p, _), settled)| settled && p.distance(at) <= tol.confusion() * 1e3)
        });
        let own = mates.get(index).is_some_and(|mate| {
            let nearer = usize::from(mate.ends[1].0.distance(at) < mate.ends[0].0.distance(at));
            mate.settled[nearer]
        });
        other || own
    }

    /// Whether the band of edge `index` runs out at its end nearest `at`
    /// as it does alone: no other edge of the request ends there, or a
    /// crossed junction across a convex turn stands there.
    pub(crate) fn alone_at(mates: &[Self], index: usize, at: Point, tol: Tolerances) -> bool {
        let near = |p: &Point| p.distance(at) <= tol.confusion() * 1e3;
        let shared = mates
            .iter()
            .enumerate()
            .any(|(i, mate)| i != index && mate.ends.iter().any(|(p, _)| near(p)));
        let crossed = mates.get(index).is_some_and(|mate| {
            mate.ends
                .iter()
                .zip(mate.alone)
                .any(|((p, _), alone)| alone && near(p))
        });
        !shared || crossed
    }

    /// Whether a junction of the request meets edge `index` at its end
    /// nearest `at`.
    pub(crate) fn meets_at(mates: &[Self], index: usize, at: Point, tol: Tolerances) -> bool {
        mates.get(index).is_some_and(|mate| {
            mate.ends
                .iter()
                .zip(mate.meets)
                .any(|((p, _), meets)| meets && p.distance(at) <= tol.confusion() * 1e3)
        })
    }
}

/// [`fillet_edge`], with the straight seat told what else was asked for.
///
/// One edge at a time stops flush against a neighbouring blend it ends on:
/// that is the state the corner tool (`round_vertex`) is built for, three
/// flush blends and a vertex. Where no blend follows, the seat is the same
/// alone and in a request, and runs out the same way onto a face cutting
/// its corner away. Edges asked for together meet: where one of them ends at
/// a vertex another leaves *tangentially*, the two are a chain and their
/// caps stand flush in one plane; where it ends against the band of one
/// that left the vertex some other way, the later seat runs on through
/// that band and the cut trims the two against each other: two fillets
/// meeting at a corner.
fn fillet_edge_meeting(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut built = fillet_edge_seat(model, solid, edge, radius, mates, tol)?;
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(built)
}

/// [`fillet_edge_meeting`] by seat, before the blend is credited to the
/// edge.
fn fillet_edge_seat(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    match route(model, solid, edge, tol)? {
        Route::Planar => planar_fillet(model, solid, edge, radius, mates, tol),
        Route::Ruled(seat) => {
            crate::ruled::ruled_fillet(model, solid, edge, radius, mates, seat, tol)
        }
        Route::Revolved(c, crange, seat) => {
            let closed =
                ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
            if closed {
                return revolved_fillet(model, solid, edge, &c, radius, tol);
            }
            let whole_turn = crate::support::runs_whole_turn(model, &seat, &c, crange, tol)?;
            if whole_turn {
                revolved_fillet(model, solid, edge, &c, radius, tol)
            } else {
                revolved_arc_fillet(model, solid, edge, &c, radius, tol)
            }
        }
        Route::Marched => crate::marched::marched_fillet(model, solid, edge, radius, mates, tol),
    }
}

/// The construction that rounds an edge, chosen by its curve and the two
/// faces it lies between.
enum Route {
    /// A straight edge between two planes.
    Planar,
    /// A straight edge between planes or drums parallel to it.
    Ruled(crate::ruled::RuledSeat),
    /// A circle between a planar cap and its coaxial wall, with the
    /// edge's range on it.
    Revolved(CircleCurve, (f64, f64), crate::support::RevolvedSeat),
    /// Any other seat: the march.
    Marched,
}

/// How `edge` of `solid` is rounded. A rim on a planar cap and its
/// coaxial wall is the revolved blend's, exact; a circle on any other pair
/// of hosts (a bore down a ball's axis, a ring's own rim) is a seat like
/// any other for the march.
fn route(model: &Model, solid: &Shape, edge: &Shape, tol: Tolerances) -> OgeomResult<Route> {
    let (curve, crange) = edge_curve(model, edge, tol)?;
    Ok(match curve {
        Curve::Line(_) if crate::support::hosts_planar(model, solid, edge, tol)? => Route::Planar,
        Curve::Line(_) => match crate::ruled::ruled_seat(model, solid, edge, tol)? {
            Some(seat) => Route::Ruled(seat),
            None => Route::Marched,
        },
        Curve::Circle(c) => match revolved_seat(model, solid, edge, &c, tol) {
            Ok(seat) => Route::Revolved(c, crange, seat),
            Err(OgeomError::Construction(_)) => Route::Marched,
            Err(e) => return Err(e),
        },
        _ => Route::Marched,
    })
}

/// Round a chain of edges in one call, each blend applied in turn and every
/// later edge resolved through the history the earlier blends wrote.
///
/// At a tangent junction the neighbouring wedges' end caps stand in one
/// plane with one cross-section, and the boolean's same-domain resolution
/// melts them: the blends join without a seam face between them. Two
/// convex edges of the chain that share a face and turn at a vertex meet
/// there, the way one edge at a time ([`fillet_edge`]) deliberately does
/// not, since the flush-ended state is what the corner tool
/// ([`round_vertex`](crate::round_vertex)) is built for. Carried on past
/// the vertex, each band is a straight band of its end section, and two of
/// one radius cross in the plane through the meeting of their centre lines
/// that bisects them, where the reflection in that plane carries the one
/// band's end section onto the other's (the shared host flat, the two
/// walls' tangent planes at the vertex reflections of each other, to the
/// fit a marched band is held to). Across a re-entrant turn (an L-shaped
/// top's inside corner) both then run on past their caps to that plane,
/// taking the corner cell they both hold, where the material holds that
/// cell; across a convex turn both stop in it where the bands' own
/// curvature keeps their sections there within that fit. Either way the
/// two are built from the solid as it stands and joined across the plane
/// into one wedge, each band's faces credited to its own edge. Where the
/// end sections differ (a drafted wall beside a square one) and both
/// creases are straight, across a re-entrant turn the corner takes what
/// lies above both bands carried on, bounded by the boolean's own
/// intersection of the two, where the material holds it all; across a
/// convex turn each band runs out past the vertex as it does alone and the
/// cuts take their union, as they do where one of the two is a drum's rim
/// (a D-shaped block's arc beside its chord). Across a re-entrant turn one
/// of the two may be a drum's rim, the material within the drum, where the
/// straight band carried on stands nowhere off the drum before the rim's
/// cap. Such
/// junctions stand beside flush ones in a loop or chain, a band crossed at
/// both ends cut in two halfway along, and two bands meeting at both their
/// ends meet at each.
/// Elsewhere (a wider convex turn of two curved bands, or of a curved band
/// other than a drum's rim beside a straight one)
/// and wherever the meeting does not close, the later seat runs on
/// through the earlier blend's band and the cut trims the two bands
/// against each other along their own intersection, or the caps stand, as
/// without a junction. Where three or
/// more edges of the chain meet at one vertex, the corner tool closes it
/// with the rolling ball's own patch (the octant of a sphere at a box
/// corner, the envelope of spheres and cylinders at a vertex no single
/// ball touches) rather than leaving the bands' caps standing. The corner
/// goes first and the bands stop flush against its patch: bands built
/// first crash into each other at an apex. A corner the tool does not
/// speak (a concave vertex, more than three surfaces with no shared
/// ball when one of them is curved) keeps its caps,
/// which is the honest picture of a corner no ball rolls around. Two
/// creases meeting on no shared face leave the wedges' caps standing
/// likewise.
///
/// Only two blends that round the same way trim each other. A wedge's cut
/// would eat a fill, and a fill cannot run on through a wedge's band, so
/// where a convex edge meets a concave one (an L-bracket's front edge at
/// its re-entrant edge) whichever is asked first takes the corner and the
/// other stops flush against its rail. The two orders then land on
/// different solids, each exact, because the corner is genuinely one
/// rounding or the other and no ball rolls round both.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`fillet_edge`] per edge, and additionally if an earlier blend
/// consumed or split a later edge: a chain whose members interfere is
/// refused rather than guessed at. A request of two or more edges also
/// refuses a ball that sets back past the far side of a face anywhere
/// along an edge, since each band's contacts are checked against the
/// others' and a contact off its face holds none to check.
pub fn fillet_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(model, solid, edges, |model, solid, parts| {
        fillet_edges_meeting(model, solid, parts, radius, tol)
    })
}

/// [`fillet_edges`] on a shape every node below which is its own, its
/// bands meeting at the junctions where two of its creases turn. Where
/// that meeting refuses, the request is rounded again from the solid as it
/// was given (its tolerances and surfaces put back), every band capped or
/// run on at its ends as it is without junctions: a junction never turns a
/// rounding into a refusal.
fn fillet_edges_meeting(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut ends: Vec<Shape> = Vec::new();
    let mut shared = false;
    for edge in edges {
        if let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? {
            for v in [a, b] {
                if ends.iter().any(|w| w.is_same(&v)) {
                    shared = true;
                }
                ends.push(v);
            }
        }
    }
    if !shared {
        return fillet_edges_own(model, solid, edges, radius, None, tol);
    }
    // What an attempt changes of the given solid where it stands: the
    // tolerances it widens, which the model's journal notes, and, once a
    // junction is taken, the surfaces a marched seat continues past its
    // faces and writes back, which the attempt notes before it builds.
    let mark = model.note_widened();
    // Crossed junctions first; where one does not close, the request is
    // rounded again with the flush junctions alone, and where that does
    // not close either, with none.
    for crossed in [true, false] {
        let mut taken = Taken {
            crossed,
            ..Taken::default()
        };
        match fillet_edges_own(model, solid, edges, radius, Some(&mut taken), tol) {
            Ok(built) => {
                model.keep_widened(mark);
                return Ok(built);
            }
            // A cancelled request stops, and a key that does not resolve
            // says the caller's model is not the one the shapes are in:
            // neither is the junctions' to answer for.
            Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => {
                model.keep_widened(mark);
                return Err(e);
            }
            // With no junction taken the attempt was the rounding without
            // one.
            Err(e) if !taken.any => {
                model.keep_widened(mark);
                return Err(e);
            }
            Err(_) => {
                model.undo_widened(mark);
                for (id, surface) in taken.surfaces {
                    if let Some(held) = model.geometry_mut().surface_mut(id) {
                        *held = surface;
                    }
                }
                // Without a crossed junction taken, the next attempt would
                // take the same flush ones again.
                if !taken.any_crossed {
                    break;
                }
            }
        }
    }
    fillet_edges_own(model, solid, edges, radius, None, tol)
}

/// What a rounding with junctions may take and took: whether crossed
/// junctions may be taken, whether any junction and any crossed one was,
/// and the solid's face surfaces as they stood before it built a band.
#[derive(Default)]
struct Taken {
    crossed: bool,
    any: bool,
    any_crossed: bool,
    surfaces: Vec<(ogeom_topo::SurfaceId, SurfaceGeometry)>,
}

/// [`fillet_edges`] on a shape every node below which is its own, with
/// junctions where `meeting` asks, which notes there what it took.
fn fillet_edges_own(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    meeting: Option<&mut Taken>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
    // The ball touches each face a setback back from the edge; set back
    // past a face's far side, the band would cut through the face, and set
    // back into another edge's band, the two balls overlap.
    let mut contacts = Vec::with_capacity(edges.len());
    crate::support::keeping_face_meshes(|| -> OgeomResult<()> {
        // A lone edge's band is checked against no other, and is held to
        // what one edge alone is.
        let alone = edges.len() == 1;
        for edge in edges {
            contacts.push(crate::support::ball_fits(
                model, solid, edge, radius, alone, tol,
            )?);
        }
        Ok(())
    })?;
    crate::support::bands_clear(model, edges, &contacts, tol)?;
    use ogeom_geom::Curve3d as _;
    // The vertices three or more of the chain's edges meet at: the corners
    // the ball rolls round, closed by the corner tool.
    let mut corners: Vec<(Shape, usize)> = Vec::new();
    for edge in edges {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        for v in [a, b] {
            if let Some(slot) = corners.iter_mut().find(|(held, _)| held.is_same(&v)) {
                slot.1 += 1;
            } else {
                corners.push((v, 1));
            }
        }
    }
    corners.retain(|(_, count)| *count >= 3);
    // A corner the tool rounds, or the solid unchanged where it does not
    // speak the corner, which keeps the bands' caps.
    let round = |model: &mut Model, built: Built, vertex: &Shape| -> OgeomResult<Built> {
        match crate::corner::round_vertex_with(model, &built.shape, vertex, radius, None, tol) {
            Ok(rounded) => Ok(Built {
                shape: rounded.shape,
                history: built.history.then(&rounded.history),
            }),
            Err(OgeomError::Construction(_)) => Ok(built),
            Err(other) => Err(other),
        }
    };
    // The corner first, the bands after it: bands built first crash into
    // each other at an apex, and the corner tool wants the sharp vertex's
    // planes, not the bands' mitre. Corners that stand clear of each other
    // are cut in one boolean, the rest one after another.
    let mut built: Option<Built> = None;
    let mut rounded: Vec<Point> = Vec::new();
    let together = if corners.len() > 1 {
        let vertices: Vec<Shape> = corners.iter().map(|(v, _)| v.clone()).collect();
        crate::corner::round_vertices_together(model, solid, &vertices, radius, tol)?
    } else {
        None
    };
    let one_by_one: &[(Shape, usize)] = if together.is_some() { &[] } else { &corners };
    if let Some(all) = together {
        for (vertex, _) in &corners {
            if let Some(point) = model
                .node(vertex)
                .and_then(|n| n.data().as_vertex())
                .map(|d| d.point)
            {
                rounded.push(vertex.transform(model.datums())?.apply(point));
            }
        }
        built = Some(all);
    }
    for (vertex, _) in one_by_one {
        let start = built
            .take()
            .unwrap_or_else(|| Built::from_nothing(solid.clone()));
        let before = start.shape.clone();
        let after = round(model, start, vertex)?;
        if !after.shape.is_same(&before)
            && let Some(point) = model
                .node(vertex)
                .and_then(|n| n.data().as_vertex())
                .map(|d| d.point)
        {
            rounded.push(vertex.transform(model.datums())?.apply(point));
        }
        built = Some(after);
    }
    let mut mates: Vec<Mate> = edges
        .iter()
        .map(|e| -> OgeomResult<Mate> {
            let (curve, range) = edge_curve(model, e, tol)?;
            let end_of = |t: f64, leaving: f64| -> OgeomResult<(Point, Vector)> {
                let d = curve.d1_at(t, tol)?;
                let m = d.magnitude();
                let unit = if m > tol.confusion() {
                    d / m
                } else {
                    Vector::ZERO
                };
                Ok((curve.point_at(t, tol)?, unit * leaving))
            };
            let ends = [end_of(range.0, 1.0)?, end_of(range.1, -1.0)?];
            let settled = ends.map(|(p, _)| {
                rounded
                    .iter()
                    .any(|q| q.distance(p) <= tol.confusion() * 1e3)
            });
            Ok(Mate {
                ends,
                settled,
                meets: [false, false],
                alone: [false, false],
            })
        })
        .collect::<OgeomResult<_>>()?;
    // The circles this chain has rounded the whole turn round. A piece of
    // a rim that closes on the solid in parts is blended the whole turn
    // round, so the rim's other pieces, asked for later or split off by the
    // boolean, are rounded already: done, not consumed. A piece of a rim
    // that stops is rounded over its own arc only, and the other pieces
    // still want theirs.
    let mut round_rims: Vec<ogeom_math::Circle> = Vec::new();
    let on_round_rim = |model: &Model, edge: &Shape, rims: &[ogeom_math::Circle]| {
        rim_circle(model, edge, tol).is_some_and(|c| rims.iter().any(|r| same_circle(r, &c, tol)))
    };
    // The two faces each edge of the chain lies between on the solid as
    // given. A blend along a crease runs the whole of it where the crease
    // closes on itself (a junction loop the march follows round), so a
    // later edge of the same crease that the blend consumed is rounded
    // already, not interfered with.
    let mut host_faces: ogeom_core::FastMap<ogeom_topo::TShapeId, Vec<Shape>> =
        ogeom_core::FastMap::default();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        for e in explore_unique(model, &face, ShapeType::Edge)? {
            let held = host_faces.entry(e.node()).or_default();
            if !held.iter().any(|f| f.is_same(&face)) {
                held.push(face.clone());
            }
        }
    }
    // Where two of the request's creases turn at a shared vertex, their
    // bands meet there: built from the solid as it stands, each stops in
    // the plane the two meet in (carried on to it past its cap across a
    // re-entrant turn), and the two are joined across it into one wedge.
    let junctions = if meeting.is_some() {
        let given: Vec<Vec<Shape>> = edges
            .iter()
            .map(|e| host_faces.get(&e.node()).cloned().unwrap_or_default())
            .collect();
        let mut found =
            crate::junction::junctions(model, solid, edges, &given, &mates, radius, tol)
                .unwrap_or_default();
        if meeting.as_ref().is_some_and(|taken| !taken.crossed) {
            found.retain(|junction| !junction.crossed);
        }
        joinable(model, edges, found)?
    } else {
        Vec::new()
    };
    if let Some(taken) = meeting
        && !junctions.is_empty()
    {
        taken.any = true;
        taken.any_crossed = junctions.iter().any(|junction| junction.crossed);
        for face in explore_unique(model, solid, ShapeType::Face)? {
            let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                continue;
            };
            if taken.surfaces.iter().any(|(id, _)| *id == data.surface) {
                continue;
            }
            if let Some(surface) = model.geometry().surface(data.surface) {
                taken.surfaces.push((data.surface, surface.clone()));
            }
        }
    }
    // A band stops at the vertex where a junction closes the corner; across
    // a crossed convex turn each runs out as it does alone.
    for junction in &junctions {
        for &index in &junction.pair {
            let mate = &mut mates[index];
            let end = usize::from(
                mate.ends[1].0.distance(junction.at) < mate.ends[0].0.distance(junction.at),
            );
            if junction.reentrant || !junction.crossed {
                mate.meets[end] = true;
            } else {
                mate.alone[end] = true;
            }
        }
    }
    let hosts_of = |edge: &Shape| -> Vec<ogeom_topo::TShapeId> {
        let mut out: Vec<_> = host_faces
            .get(&edge.node())
            .map(|faces| faces.iter().map(Shape::node).collect())
            .unwrap_or_default();
        out.sort_unstable();
        out
    };
    // The edges of the solid a round's blends are built on, by occurrence:
    // until a round is applied that solid is unchanged, and an edge of it
    // is found there by name rather than re-found by its geometry.
    let mut occurrences: Option<(Shape, ogeom_core::FastMap<ogeom_topo::SameKey, Shape>)> = None;
    let mut refind =
        |model: &Model, on: &Shape, edge: &Shape, hosts: &[Shape]| -> OgeomResult<Vec<Shape>> {
            if !occurrences
                .as_ref()
                .is_some_and(|(held, _)| held.is_same(on))
            {
                let mut by_name = ogeom_core::FastMap::default();
                for e in explore_unique(model, on, ShapeType::Edge)? {
                    by_name.entry(ogeom_topo::SameKey(e.clone())).or_insert(e);
                }
                occurrences = Some((on.clone(), by_name));
            }
            let named = occurrences
                .as_ref()
                .and_then(|(_, by_name)| by_name.get(&ogeom_topo::SameKey(edge.clone())));
            match named {
                Some(own) => Ok(vec![own.clone()]),
                None => refind_edges(model, on, edge, hosts, tol),
            }
        };
    let mut blended_creases: Vec<Vec<ogeom_topo::TShapeId>> = Vec::new();
    // Edges that share no vertex round without meeting: each blend is built
    // from the faces as they stand, and theirs can be applied together.
    // Edges that do share one are rounded in turn, the later built against
    // what the earlier left. So the chain is taken in rounds, no two edges
    // of a round sharing a vertex (pieces of one circle are built from the
    // same faces: the first rounds the whole turn where the rim closes, and
    // where it stops each piece is capped at its ends, flush with its
    // neighbour's cap), each round's wedges set
    // aside and applied together, one boolean each way.
    let rims: Vec<Option<Circle>> = edges.iter().map(|e| rim_circle(model, e, tol)).collect();
    let mut ends_of: Vec<Vec<Shape>> = Vec::with_capacity(edges.len());
    for edge in edges {
        let mut ends: Vec<Shape> = Vec::new();
        if let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? {
            for v in [a, b] {
                if !ends.iter().any(|held| held.is_same(&v)) {
                    ends.push(v);
                }
            }
        }
        ends_of.push(ends);
    }
    let meet = |i: usize, j: usize| {
        ends_of[i]
            .iter()
            .any(|v| ends_of[j].iter().any(|w| w.is_same(v)))
            && !matches!((&rims[i], &rims[j]), (Some(a), Some(b)) if same_circle(a, b, tol))
            // Two bands meeting in a junction's plane are built from the
            // solid as it stands and joined into one wedge.
            && !junctions
                .iter()
                .any(|junction| junction.slot(i).is_some() && junction.slot(j).is_some())
    };
    // Where a marched band and an analytic one meet at a vertex, the
    // marched band is built first, against the faces as they stand: its
    // run-out then melts with untouched hosts, and the analytic band built
    // after it runs on through it, which the melt resolves. Built the other
    // way round, the marched band's run-out lands on the analytic band's
    // rail along a stretch, a contact the melt does not resolve. Each edge
    // is routed on the solid as given: earlier rounds trim an edge's hosts
    // but leave their surfaces, which are what the route reads. The order
    // the edges were asked in decides everything else.
    let mut marched: Vec<Option<bool>> = vec![None; edges.len()];
    for i in 0..edges.len() {
        if (0..edges.len()).any(|j| j != i && meet(i, j)) {
            marched[i] = Some(matches!(
                route(model, solid, &edges[i], tol)?,
                Route::Marched
            ));
        }
    }
    let first: Vec<bool> = (0..edges.len())
        .map(|i| {
            marched[i] == Some(true)
                && (0..edges.len()).any(|j| meet(i, j) && marched[j] == Some(false))
        })
        .collect();
    let mut taken: Vec<usize> = (0..edges.len()).collect();
    taken.sort_by_key(|&i| (!first[i], i));
    // Bands joined across a junction's plane are built in one round: each
    // joined group takes the first round in which none of its members
    // meets an edge already placed there.
    let mut group: Vec<usize> = (0..edges.len()).collect();
    for junction in &junctions {
        let [a, b] = junction.pair;
        let (from, to) = (group[a], group[b]);
        for g in &mut group {
            if *g == from {
                *g = to;
            }
        }
    }
    let mut round_of: Vec<Option<usize>> = vec![None; edges.len()];
    for &i in &taken {
        if round_of[i].is_some() {
            continue;
        }
        let members: Vec<usize> = (0..edges.len()).filter(|&m| group[m] == group[i]).collect();
        let mut round = 0;
        while members
            .iter()
            .any(|&m| (0..edges.len()).any(|j| round_of[j] == Some(round) && meet(m, j)))
        {
            round += 1;
        }
        for m in members {
            round_of[m] = Some(round);
        }
    }
    let round_of: Vec<usize> = round_of.into_iter().map(|r| r.unwrap_or(0)).collect();
    let rank: Vec<usize> = {
        let mut rank = vec![0; edges.len()];
        for (k, &i) in taken.iter().enumerate() {
            rank[i] = k;
        }
        rank
    };
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_by_key(|&i| (round_of[i], rank[i]));
    let mut set_aside: Vec<(Vec<usize>, Shape, crate::support::Wedge)> = Vec::new();
    let mut this_round = 0;
    for index in order {
        let edge = &edges[index];
        // A round's wedges are applied before the next round is built.
        if round_of[index] != this_round {
            this_round = round_of[index];
            if let Some(b) = built.take() {
                let halves = with_corners(model, &mut set_aside, &junctions, &b.shape, tol)?;
                let joined = join_flush(
                    model,
                    core::mem::take(&mut set_aside),
                    &junctions,
                    &halves,
                    tol,
                )?;
                built = Some(apply_set_aside(model, b, joined, tol)?);
            }
        }
        // A piece of a circle whose whole turn another piece's blend rounds.
        if on_round_rim(model, edge, &round_rims) {
            if let Some(b) = built.as_mut() {
                b.history.delete(edge);
            }
            continue;
        }
        // The faces the edge lies between on the solid as given: what a
        // piece of it re-found later still lies between.
        let hosts: Vec<Shape> = host_faces.get(&edge.node()).cloned().unwrap_or_default();
        // The edge as it stands on the current solid: itself on the first
        // step, and afterwards whatever the earlier blends left of it: one
        // re-found stand-in, or the pieces a blend running out across it
        // split it into, each of which is a seat of its own ending against
        // that blend's band, which is exactly the corner it then meets.
        // The pieces are re-found from the edge itself, between its hosts:
        // the history also traces it to stubs along its line where an
        // earlier blend's cap meets a host tangentially, no corner to round.
        let (current, targets): (Shape, Vec<Shape>) = match &built {
            None => (solid.clone(), vec![edge.clone()]),
            Some(b) => {
                let gone = b.history.trace(edge).is_empty();
                // An edge both of whose ends this chain's corners rounded,
                // and which they left nothing of, ran its whole length
                // inside their balls: the corners are its blend.
                let on_blended_crease = blended_creases.contains(&hosts_of(edge))
                    || mates[index].settled == [true, true];
                let rounded = on_blended_crease || on_round_rim(model, edge, &round_rims);
                if gone && rounded {
                    continue;
                }
                if gone {
                    ogeom_bail!(
                        Construction,
                        "an earlier blend in the chain consumed this edge; \
                         the chain's members interfere"
                    );
                }
                let found = match refind(model, &b.shape, edge, &hosts) {
                    Ok(live) => live,
                    Err(_) if rounded => Vec::new(),
                    Err(e) => return Err(e),
                };
                (b.shape.clone(), found)
            }
        };
        // Nothing left of it: the crease it lies on is rounded, and the
        // edge with it.
        if targets.is_empty() {
            if let Some(b) = built.as_mut() {
                b.history.delete(edge);
            }
            continue;
        }
        let mut current = current;
        for target in &targets {
            // Each piece's blend replaces the solid; the next piece is
            // re-found on what that blend left.
            let live = match refind(model, &current, target, &hosts) {
                Ok(live) => live,
                Err(_) if on_round_rim(model, target, &round_rims) => continue,
                Err(e) => return Err(e),
            };
            let [target] = live.as_slice() else {
                ogeom_bail!(
                    Construction,
                    "a blend in the chain split a piece of a later edge again; \
                     the chain's members interfere"
                );
            };
            let whole_rim = match rim_circle(model, target, tol) {
                Some(rim) if rounds_whole_turn(model, &current, target, tol)? => Some(rim),
                _ => None,
            };
            let mut step = {
                let (step, wedges) = crate::support::collecting_wedges(|| {
                    fillet_edge_meeting(model, &current, target, radius, Some((index, &mates)), tol)
                });
                let step = step?;
                for mut wedge in wedges {
                    if !wedge.additive {
                        for junction in &junctions {
                            let Some(slot) = junction.slot(index) else {
                                continue;
                            };
                            // A crossed junction's bands stand as built:
                            // a re-entrant one's corner is built from the
                            // two once the round has set both aside.
                            if junction.crossed {
                                continue;
                            }
                            let plane = junction.keeps(slot);
                            // Across a re-entrant turn the band is carried
                            // on past its cap to the plane it meets the
                            // other in; across a convex one it is trimmed
                            // back to it. Nothing of the band but its end
                            // stands within two radii of the vertex across
                            // the plane: a band of its own edge coming back
                            // that close would leave no room for its ball.
                            if junction.reentrant {
                                let Some(cap) = crate::junction::cap_at(
                                    model,
                                    &wedge.solid,
                                    junction.at,
                                    junction.outward[slot],
                                    tol,
                                )?
                                else {
                                    ogeom_bail!(
                                        Construction,
                                        "a band meeting another across a re-entrant turn has \
                                         no cap at the vertex to carry on from"
                                    );
                                };
                                wedge.solid = crate::junction::carried_into(
                                    model,
                                    &wedge.solid,
                                    &cap,
                                    (edge, &hosts),
                                    junction.outward[slot],
                                    plane,
                                    tol,
                                )?;
                            } else {
                                wedge.solid = crate::support::trimmed_wedge(
                                    model,
                                    &wedge.solid,
                                    plane,
                                    radius * 2.0,
                                    tol,
                                )?;
                            }
                        }
                    }
                    set_aside.push((vec![index], edge.clone(), wedge));
                }
                step
            };
            if let Some(rim) = whole_rim {
                round_rims.push(rim);
            }
            let hosts = hosts_of(edge);
            if hosts.len() == 2 && !blended_creases.contains(&hosts) {
                blended_creases.push(hosts);
            }
            // The blend consumed the re-found stand-in; the caller's edge is
            // the same fact under its original name.
            if !target.is_same(edge) {
                step.history.delete(edge);
            }
            current = step.shape.clone();
            built = Some(match built {
                None => step,
                Some(prev) => Built {
                    shape: step.shape.clone(),
                    history: prev.history.then(&step.history),
                },
            });
        }
    }
    let built = built.unwrap_or_else(|| unreachable!());
    let halves = with_corners(model, &mut set_aside, &junctions, &built.shape, tol)?;
    let joined = join_flush(model, set_aside, &junctions, &halves, tol)?;
    apply_set_aside(model, built, joined, tol)
}

/// The junctions of `found` whose bands can be joined: taken in turn, each
/// joins the groups of bands already joined at its two edges, and is
/// dropped where two edges of those groups also meet at another vertex,
/// which joined into one round would lose the later band's running on
/// through the earlier there.
fn joinable(
    model: &Model,
    edges: &[Shape],
    found: Vec<crate::junction::Junction>,
) -> OgeomResult<Vec<crate::junction::Junction>> {
    let mut ends: Vec<Vec<Shape>> = Vec::with_capacity(edges.len());
    for edge in edges {
        ends.push(
            ogeom_algo::edge_vertices(model, edge)?
                .map(|(a, b)| vec![a, b])
                .unwrap_or_default(),
        );
    }
    // Flush junctions are taken first: a crossed junction joins its two
    // bands only where that leaves every flush one standing.
    let (flush, crossed): (Vec<_>, Vec<_>) = found.into_iter().partition(|j| !j.crossed);
    let mut found = flush.clone();
    found.extend(crossed);
    let pending: Vec<[usize; 2]> = found.iter().map(|j| j.pair).collect();
    let mut group: Vec<usize> = (0..edges.len()).collect();
    let mut kept: Vec<crate::junction::Junction> = Vec::new();
    for junction in found {
        let [a, b] = junction.pair;
        if group[a] == group[b] {
            // A junction closing a loop of joined bands joins no groups,
            // and is taken for its corner. Each junction stands at a vertex
            // of its own, so two bands meeting at both their ends (a pocket
            // of a chord and its arc) close a loop of two there.
            kept.push(junction);
            continue;
        }
        // Two edges of the groups this joins that meet elsewhere must meet
        // there in a junction of their own: one taken, this one, or a
        // crossed one still to come.
        let joined = |m: usize, n: usize| {
            kept.iter()
                .chain(std::iter::once(&junction))
                .any(|j| j.slot(m).is_some() && j.slot(n).is_some())
                || pending
                    .iter()
                    .any(|pair| pair.contains(&m) && pair.contains(&n))
        };
        let of_a: Vec<usize> = (0..edges.len()).filter(|&m| group[m] == group[a]).collect();
        let of_b: Vec<usize> = (0..edges.len()).filter(|&m| group[m] == group[b]).collect();
        let clash = of_a.iter().any(|&m| {
            of_b.iter().any(|&n| {
                !joined(m, n) && ends[m].iter().any(|v| ends[n].iter().any(|w| w.is_same(v)))
            })
        });
        if clash {
            continue;
        }
        let (from, to) = (group[a], group[b]);
        for g in &mut group {
            if *g == from {
                *g = to;
            }
        }
        kept.push(junction);
    }
    // Every two edges of one joined group that share a vertex meet there
    // in a kept junction; where one was dropped, the request is taken with
    // its flush junctions alone.
    let meets_unjoined = (0..edges.len()).any(|m| {
        (m + 1..edges.len()).any(|n| {
            group[m] == group[n]
                && ends[m].iter().any(|v| ends[n].iter().any(|w| w.is_same(v)))
                && !kept
                    .iter()
                    .any(|j| j.slot(m).is_some() && j.slot(n).is_some())
        })
    });
    if meets_unjoined && kept.iter().any(|j| j.crossed) {
        return joinable(model, edges, flush);
    }
    Ok(kept)
}

/// Join, among the wedges a round set aside, the two bands of each crossed
/// re-entrant junction the round built into the one wedge that also takes
/// their corner, in the material of `solid`, the solid the round's wedges
/// apply to. The joined wedge is set aside with the first band's edge, its
/// faces from the second band credited to the second's.
///
/// A band with a crossed junction at each end is first cut in two in the
/// plane square to its crease halfway along: each half is the straight
/// band of its own end section, joined into its own corner, and the two
/// corners' wedges are joined again across that plane, where both hold the
/// band's one section. Returns those planes, each with the band's index.
fn with_corners(
    model: &mut Model,
    set_aside: &mut Vec<(Vec<usize>, Shape, crate::support::Wedge)>,
    junctions: &[crate::junction::Junction],
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<(usize, (Point, Vector))>> {
    let crossed: Vec<&crate::junction::Junction> = junctions
        .iter()
        .filter(|j| j.crossed && j.reentrant)
        .collect();
    let holding = |set_aside: &[(Vec<usize>, Shape, crate::support::Wedge)], index: usize| {
        (0..set_aside.len())
            .filter(|&k| set_aside[k].0.contains(&index))
            .collect::<Vec<usize>>()
    };
    let mut halves = Vec::new();
    let indices: Vec<usize> = crossed.iter().flat_map(|j| j.pair).collect();
    for &index in &indices {
        let ends: Vec<Point> = crossed
            .iter()
            .filter(|j| j.slot(index).is_some())
            .map(|j| j.at)
            .collect();
        if ends.len() != 2 || halves.iter().any(|(i, _)| *i == index) {
            continue;
        }
        let held = holding(set_aside, index);
        let [k] = held.as_slice() else {
            continue;
        };
        let k = *k;
        let along = ends[1] - ends[0];
        let length = along.magnitude();
        if length <= tol.confusion() {
            ogeom_bail!(Construction, "a band crossed at both ends has no length");
        }
        let toward = along / length;
        let middle = ends[0].midpoint(ends[1]);
        let span = length + crate::marched::band_fit_target(tol).max(length);
        let (indices_k, edge, wedge) = set_aside.remove(k);
        for keep in [-toward, toward] {
            let half =
                crate::support::trimmed_wedge(model, &wedge.solid, (middle, keep), span, tol)?;
            set_aside.push((
                indices_k.clone(),
                edge.clone(),
                crate::support::Wedge {
                    solid: half,
                    additive: wedge.additive,
                    edges: wedge.edges.clone(),
                    melt_refusal: wedge.melt_refusal,
                    credit: Vec::new(),
                },
            ));
        }
        halves.push((index, (middle, toward)));
    }
    for junction in crossed {
        // Each side's band: the one wedge holding it, or of a band cut in
        // two, the half ending at this junction's vertex.
        let mut sides = [None, None];
        for (slot, side) in sides.iter_mut().enumerate() {
            // A band not yet joined into another corner holds no other.
            let held = holding(set_aside, junction.pair[slot])
                .into_iter()
                .filter(|&k| set_aside[k].0.len() == 1);
            let mut here = Vec::new();
            for k in held {
                if crate::junction::cap_at(
                    model,
                    &set_aside[k].2.solid,
                    junction.at,
                    junction.outward[slot],
                    tol,
                )?
                .is_some()
                {
                    here.push(k);
                }
            }
            *side = match here.as_slice() {
                [] => None,
                [k] => Some(*k),
                _ => ogeom_bail!(
                    Construction,
                    "a crossed junction's band is more than one cut wedge in its round"
                ),
            };
        }
        // Built in the same round, or not yet.
        let (a, b) = match sides {
            [None, None] => continue,
            [Some(a), Some(b)] if a != b => (a, b),
            _ => ogeom_bail!(
                Construction,
                "a crossed junction's bands are not one cut wedge each in their round"
            ),
        };
        if set_aside[a].2.additive || set_aside[b].2.additive {
            ogeom_bail!(Construction, "a crossed junction's band is a fill");
        }
        if !set_aside[a].2.credit.is_empty() || !set_aside[b].2.credit.is_empty() {
            ogeom_bail!(Construction, "a crossed junction's band is already joined");
        }
        let (taken, of_second) = crate::junction::corner(
            model,
            [&set_aside[a].2.solid, &set_aside[b].2.solid],
            junction,
            solid,
            tol,
        )?;
        let second = set_aside[b].1.clone();
        let later = set_aside.remove(a.max(b));
        let earlier = set_aside.remove(a.min(b));
        let ((mut indices, edge, first), (others, other, last)) = if a < b {
            (earlier, later)
        } else {
            (later, earlier)
        };
        indices.extend(others);
        let mut rounded = first.edges;
        rounded.extend(last.edges);
        rounded.push(other);
        set_aside.push((
            indices,
            edge,
            crate::support::Wedge {
                solid: taken,
                additive: false,
                edges: rounded,
                melt_refusal: first.melt_refusal.or(last.melt_refusal),
                credit: of_second
                    .into_iter()
                    .map(|face| (face, second.clone()))
                    .collect(),
            },
        ));
    }
    Ok(halves)
}

/// The wedges a round set aside, those of two bands meeting in a
/// junction's plane joined into one across it, each wedge with the edge it
/// rounds.
fn join_flush(
    model: &mut Model,
    set_aside: Vec<(Vec<usize>, Shape, crate::support::Wedge)>,
    junctions: &[crate::junction::Junction],
    halves: &[(usize, (Point, Vector))],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, crate::support::Wedge)>> {
    // Each entry with every edge of the request whose band it holds.
    let mut held: Vec<Option<(Vec<usize>, Shape, crate::support::Wedge)>> =
        set_aside.into_iter().map(Some).collect();
    // The flush junctions, and the planes a band crossed at both ends was
    // cut in, each as the plane, a point of it the meeting stands at, and
    // the edges on either side.
    let meetings: Vec<((Point, Vector), Point, [usize; 2])> = junctions
        .iter()
        .filter(|j| !j.crossed)
        .map(|j| (j.flush, j.at, j.pair))
        .chain(
            halves
                .iter()
                .map(|&(index, (on, normal))| ((on, normal), on, [index, index])),
        )
        .collect();
    for ((on, normal), at, pair) in meetings {
        // The wedge holding each side's piece that ends in the plane, the
        // two sides two wedges.
        let mut sides = [None, None];
        for slot in 0..2 {
            for (k, entry) in held.iter().enumerate() {
                let Some((indices, _, wedge)) = entry else {
                    continue;
                };
                if !indices.contains(&pair[slot]) || (slot == 1 && sides[0] == Some(k)) {
                    continue;
                }
                let gap = crate::marched::band_fit_target(tol) * 2.0;
                if crate::junction::ends_in(model, &wedge.solid, (on, normal), at, gap, tol)? {
                    sides[slot] = Some(k);
                    if slot == 0 {
                        break;
                    }
                }
            }
        }
        let [Some(a), Some(b)] = sides else {
            continue;
        };
        if a == b {
            // A loop closing on itself: its wedge's two ends in the plane
            // are sewn to each other.
            let Some((indices, edge, ring)) = held[a].take() else {
                continue;
            };
            let (solid, traced) = crate::junction::joined(
                model,
                [&ring.solid, &ring.solid],
                (on, normal),
                at,
                crate::marched::band_fit_target(tol) * 2.0,
                tol,
            )?;
            let mut credit = Vec::new();
            for (given, image) in &traced {
                let owner = ring
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(given))
                    .map_or_else(|| edge.clone(), |(_, e)| e.clone());
                credit.push((image.clone(), owner));
            }
            held[a] = Some((
                indices,
                edge,
                crate::support::Wedge {
                    solid,
                    credit,
                    ..ring
                },
            ));
            continue;
        }
        let (Some((mut indices, edge, first)), Some((others, other, second))) =
            (held[a].take(), held[b].take())
        else {
            continue;
        };
        let (solid, traced) = crate::junction::joined(
            model,
            [&first.solid, &second.solid],
            (on, normal),
            at,
            crate::marched::band_fit_target(tol) * 2.0,
            tol,
        )?;
        // Each face of the joined wedge stays credited to the edge whose
        // band it came from.
        let mut credit: Vec<(Shape, Shape)> = Vec::new();
        for (from, rounded) in [(&first, &edge), (&second, &other)] {
            for face in ogeom_topo::explore_unique(model, &from.solid, ShapeType::Face)? {
                let owner = from
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or_else(|| rounded.clone(), |(_, e)| e.clone());
                for (given, image) in &traced {
                    if given.is_same(&face) {
                        credit.push((image.clone(), owner.clone()));
                    }
                }
            }
        }
        indices.extend(others);
        indices.dedup();
        let mut rounded: Vec<Shape> = first.edges;
        rounded.extend(second.edges);
        if !other.is_same(&edge) {
            rounded.push(other);
        }
        held[a] = Some((
            indices,
            edge,
            crate::support::Wedge {
                solid,
                additive: false,
                edges: rounded,
                melt_refusal: first.melt_refusal.or(second.melt_refusal),
                credit,
            },
        ));
    }
    Ok(held
        .into_iter()
        .flatten()
        .map(|(_, edge, wedge)| (edge, wedge))
        .collect())
}

/// Apply the wedges a chain set aside: the cut ones together and the
/// fused ones together, each blend's faces credited to its edge through
/// the boolean's history. Wedges whose boxes meet go in different
/// booleans, as a compound of lumps that overlap would count their overlap
/// twice. A cut is [`ogeom_bool::cut_seated`], which trims the faces
/// straight bands and rims sit on where nothing else reaches into them, and
/// a fill [`ogeom_bool::fuse_seated`], which does so for rims.
fn apply_set_aside(
    model: &mut Model,
    mut built: Built,
    set_aside: Vec<(Shape, crate::support::Wedge)>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if set_aside.is_empty() {
        return Ok(built);
    }
    let mut boxes = Vec::with_capacity(set_aside.len());
    for (_, wedge) in &set_aside {
        boxes.push(ogeom_algo::shape_bounds(model, &wedge.solid, tol)?);
    }
    let (fused, cut): (Vec<usize>, Vec<usize>) =
        (0..set_aside.len()).partition(|&i| set_aside[i].1.additive);
    // A cut and a fill whose boxes meet are applied in the order they were
    // built, one by one: which goes first decides what the corner keeps.
    let crossed = cut
        .iter()
        .any(|&i| fused.iter().any(|&j| boxes[i].intersects(&boxes[j])));
    let groups: Vec<Vec<usize>> = if crossed {
        (0..set_aside.len()).map(|i| vec![i]).collect()
    } else {
        // Each side's wedges in as few groups as keep every group's boxes
        // apart, each group one boolean: blends close side by side (rims on
        // a tight grid) take a handful of booleans, not one each.
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for side in [cut, fused] {
            let first = groups.len();
            for i in side {
                let slot = groups[first..]
                    .iter()
                    .position(|g| g.iter().all(|&j| !boxes[i].intersects(&boxes[j])));
                match slot {
                    Some(k) => groups[first + k].push(i),
                    None => groups.push(vec![i]),
                }
            }
        }
        groups
    };
    let given = built.clone();
    let mut groups: std::collections::VecDeque<Vec<usize>> = groups.into();
    while let Some(group) = groups.pop_front() {
        let additive = set_aside[group[0]].1.additive;
        let lumps: Vec<Shape> = group
            .iter()
            .map(|&i| set_aside[i].1.solid.clone())
            .collect();
        let tool = if let [one] = lumps.as_slice() {
            one.clone()
        } else {
            model.add_compound(&lumps)?
        };
        let applied = if additive {
            ogeom_bool::fuse_seated(model, &built.shape, &tool, tol)
        } else {
            ogeom_bool::cut_seated(model, &built.shape, &tool, tol)
        };
        // A group the boolean does not resolve at once is no refusal: the
        // wedges are then applied one at a time, in the order they were
        // built, to the solid as it was given.
        if matches!(applied, Err(OgeomError::NotDone(_))) && group.len() > 1 {
            built = given.clone();
            groups = (0..set_aside.len()).map(|i| vec![i]).collect();
            continue;
        }
        let refusal = group.iter().find_map(|&i| set_aside[i].1.melt_refusal);
        let step = match (applied, refusal) {
            (Err(OgeomError::NotDone(cause)), Some(refusal)) => {
                ogeom_bail!(NotDone, "{refusal} ({cause})")
            }
            (applied, _) => applied?,
        };
        // The chain's history already says each edge is gone; the blends
        // it generated are recorded on the composed history, where an edge
        // consumed before this step still answers for them.
        let mut history = built.history.then(&step.history);
        let result = ogeom_topo::explore_unique(model, &step.shape, ShapeType::Face)?;
        for &i in &group {
            let (edge, wedge) = &set_aside[i];
            for face in ogeom_topo::explore_unique(model, &wedge.solid, ShapeType::Face)? {
                let owner = wedge
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or(edge, |(_, e)| e);
                for image in step.history.trace(&face) {
                    if result.iter().any(|f| f.is_same(image)) {
                        history.generate(owner, image.clone());
                    }
                }
            }
            for rounded in &wedge.edges {
                history.delete(rounded);
            }
            history.delete(edge);
        }
        built = Built {
            shape: step.shape.clone(),
            history,
        };
    }
    Ok(built)
}

/// Round a straight convex edge with a blend whose radius runs linearly from
/// `start_radius` at the edge's start to `end_radius` at its end.
///
/// For a linear law on a straight edge between planes the rolling ball's
/// envelope is *exactly* a rational B-spline surface (degree one along the
/// edge, a rational quadratic arc across it, the control net affine in the
/// radius), so nothing here is fitted. The tangency lines are straight, the
/// legs stay planar, and the wedge subtracts through the boolean like its
/// constant-radius siblings.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`fillet_edge`], for the straight planar seat only.
pub fn fillet_edge_variable(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    start_radius: f64,
    end_radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| {
            fillet_edge_variable_own(model, solid, &parts[0], start_radius, end_radius, tol)
        },
    )
}

/// [`fillet_edge_variable`] on a shape every node below which is its own.
fn fillet_edge_variable_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    start_radius: f64,
    end_radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for radius in [start_radius, end_radius] {
        if !radius.is_finite() || radius <= tol.confusion() {
            ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
        }
    }
    if (start_radius - end_radius).abs() <= tol.confusion() {
        return fillet_edge_own(model, solid, edge, start_radius, tol);
    }
    let seat = planar_seat(model, solid, edge, tol)?;
    let (first, second) = if seat
        .along
        .dot(seat.normals[0].cross(seat.normals[1]))
        .is_sign_positive()
    {
        (0, 1)
    } else {
        (1, 0)
    };
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let n1 = seat.normals[first] * sign;
    let n2 = seat.normals[second] * sign;
    let a = seat.leg(first, tol)? * sign;
    let b = seat.leg(second, tol)? * sign;
    let bisector = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
        }
        u / m
    };
    let depth = -n1.dot(bisector);
    if depth <= tol.angular() {
        ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
    }
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2)) * sign;
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }

    // The sections at the two ends: everything else is affine between them.
    // On a concave edge every direction above is already mirrored: the
    // arithmetic below cannot tell which case it serves.
    let radii = [start_radius, end_radius];
    let apex = [seat.start, seat.end];
    let centre = [
        apex[0] + bisector * (radii[0] / depth),
        apex[1] + bisector * (radii[1] / depth),
    ];
    let contact_a = [centre[0] + n1 * radii[0], centre[1] + n1 * radii[1]];
    let contact_b = [centre[0] + n2 * radii[0], centre[1] + n2 * radii[1]];

    // The unit arc from n1 to n2 as a rational quadratic: the middle control
    // point is where the end tangents meet, its weight the half-angle cosine.
    let half_cos = (sweep / 2.0).cos();
    let arc_mid = (n1 + n2) / (1.0 + n1.dot(n2));
    let weights = [1.0, half_cos, 1.0];

    // The blend: degree one by rational quadratic, the control net
    // c_i + r_i * V_j, weights carried across the rows unchanged. Every
    // section of the surface is that section's exact blend arc.
    let u_knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
    let v_knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
    let mut net: Vec<ogeom_math::Weighted<Point>> = Vec::with_capacity(6);
    for i in 0..2 {
        let arms = [n1, arc_mid, n2];
        for (j, arm) in arms.iter().enumerate() {
            net.push(ogeom_math::Weighted::new(
                centre[i] + *arm * radii[i],
                weights[j],
                tol,
            )?);
        }
    }
    let grid = ogeom_math::ControlGrid::new(net, 2, 3)?;
    let surface: SurfaceGeometry =
        ogeom_geom::BSplineSurface::rational(u_knots, v_knots, grid)?.into();
    let surface_id = model.geometry_mut().add_surface(surface.clone());

    // Shared vertices at the four tangency corners.
    let va = [
        make_vertex(model, contact_a[0]).shape,
        make_vertex(model, contact_a[1]).shape,
    ];
    let vb = [
        make_vertex(model, contact_b[0]).shape,
        make_vertex(model, contact_b[1]).shape,
    ];

    // The tangency lines, straight because the law is linear, with
    // degree-one pcurves mapping their arc length onto the chart's rows.
    let rail = |model: &mut Model,
                from: (&Shape, Point),
                to: (&Shape, Point),
                row: f64|
     -> OgeomResult<Shape> {
        let line = ogeom_geom::LineCurve::segment(from.1, to.1, tol)?;
        let curve = Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        let built = make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape;
        let knots = ogeom_math::KnotVector::new(vec![domain.0, domain.0, domain.1, domain.1], 1)?;
        let pcurve = ogeom_geom::BSpline2d::new(
            knots,
            vec![
                ogeom_math::Point2::new(0.0, row),
                ogeom_math::Point2::new(1.0, row),
            ],
            tol,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            pcurve.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            domain,
        )?;
        Ok(built)
    };
    let rail_a = rail(model, (&va[0], contact_a[0]), (&va[1], contact_a[1]), 0.0)?;
    let rail_b = rail(model, (&vb[0], contact_b[0]), (&vb[1], contact_b[1]), 1.0)?;

    // The end arcs as the surface's own rational sections, so the cap edges
    // and the blend's chart speak the same parameter.
    let arc_edge = |model: &mut Model, i: usize| -> OgeomResult<Shape> {
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
        let control = vec![
            ogeom_math::Weighted::new(contact_a[i], 1.0, tol)?,
            ogeom_math::Weighted::new(centre[i] + arc_mid * radii[i], half_cos, tol)?,
            ogeom_math::Weighted::new(contact_b[i], 1.0, tol)?,
        ];
        let curve = Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?);
        let built = make_edge_between(model, curve, (0.0, 1.0), &va[i], &vb[i], tol)?.shape;
        #[allow(clippy::cast_precision_loss)]
        let column = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                ogeom_math::Point2::new(i as f64, 0.0),
                ogeom_math::Direction2::Y,
            ),
            0.0,
            1.0,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            column.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, 1.0),
        )?;
        Ok(built)
    };
    let arc0 = arc_edge(model, 0)?;
    let arc1 = arc_edge(model, 1)?;

    // The blend face on the registered surface, oriented so its outward side
    // leaves the wedge, decided by measurement at the middle rather than by
    // convention.
    let blend = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                rail_a.clone(),
                arc1.clone(),
                rail_b.reversed(),
                arc0.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        use ogeom_geom::Surface as _;
        let s_mid = surface.point_at(0.5, 0.5, tol)?;
        let (du, dv) = surface.d1_at(0.5, 0.5, tol)?;
        let apex_mid = apex[0].midpoint(apex[1]);
        let outward = s_mid - apex_mid.midpoint(s_mid);
        if du.cross(dv).dot(outward) >= 0.0 {
            face
        } else {
            face.reversed()
        }
    };

    // Legs along the two faces, and caps closed by the end arcs. The caps go
    // through the generic builder: the rational arc projects into a plane
    // exactly, control point by control point.
    let leg_a = planar_face(
        model,
        &[apex[0], contact_a[0], contact_a[1], apex[1]],
        n1 * sign,
        tol,
    )?;
    let leg_b = planar_face(
        model,
        &[apex[0], contact_b[0], contact_b[1], apex[1]],
        n2 * sign,
        tol,
    )?;
    let cap = |model: &mut Model, i: usize, outward: Vector| -> OgeomResult<Built> {
        let plane = Plane::through(apex[i], Direction::new(outward, tol)?);
        let reach = (apex[i].distance(contact_a[i]) + apex[i].distance(contact_b[i]) + 1.0) * 2.0;
        let cap_surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let apex_v = make_vertex(model, apex[i]).shape;
        let edges = vec![
            crate::support::segment_between(
                model,
                (&apex_v, apex[i]),
                (&va[i], contact_a[i]),
                tol,
            )?,
            if i == 0 { arc0.clone() } else { arc1.clone() },
            crate::support::segment_between(
                model,
                (&vb[i], contact_b[i]),
                (&apex_v, apex[i]),
                tol,
            )?,
        ];
        // Apex, contact on `a`, contact on `b`: walked back where that
        // turns against `outward`, so the cap keeps itself on the left.
        let turn = (contact_a[i] - apex[i]).cross(contact_b[i] - apex[i]);
        let edges: Vec<Shape> = if turn.dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        ogeom_algo::make_face_with_pcurves(model, cap_surface, &[edges], tol)
    };
    let cap0 = cap(model, 0, -seat.along)?.shape;
    let cap1 = cap(model, 1, seat.along)?.shape;

    let faces = [leg_a, leg_b, cap0, cap1, blend];
    apply_wedge(model, solid, Some(edge), &faces, !seat.convex, tol)
}

/// The straight-edge fillet: the rolling ball's envelope is a cylinder.
fn planar_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut seat = planar_seat(model, solid, edge, tol)?;
    // Where a convex crease runs into a face overhanging its corner (a rib
    // running up under a slanted roof), the band is trimmed by that face's
    // plane: carried on until the plane has passed its whole section, and
    // what lies behind the plane cut off it, so the band stops on the face
    // rather than cutting into what it bounds. A face square to the crease
    // meets the band's cap in the cap's own plane, and the cap stands.
    let mut trims: Vec<(Point, Vector)> = Vec::new();
    let mut trimmed = [false; 2];
    for (slot, end) in [true, false].into_iter().enumerate() {
        if !seat.convex {
            break;
        }
        let at = if end { seat.start } else { seat.end };
        // A junction of the request closes the corner at this end.
        if mates.is_some_and(|(index, mates)| Mate::meets_at(mates, index, at, tol)) {
            continue;
        }
        let outward = if end { -seat.along } else { seat.along };
        let hosts = [&seat.faces[0], &seat.faces[1]];
        let Some((origin, normal)) =
            crate::support::overhang_at(model, solid, edge, hosts, at, outward, tol)?
        else {
            continue;
        };
        let Some(inward) = seat_centre(&seat, radius, tol)? else {
            continue;
        };
        // Where each point of the section crosses the plane, carried along
        // `outward`: a plane passes a section last at one of its corners.
        let rate = -outward.dot(normal);
        let crossings =
            section_at(&seat, at, inward, radius).map(|p| (origin - p).dot(normal) / rate);
        if crossings.iter().all(|c| c.abs() <= tol.confusion() * 1e3) {
            continue;
        }
        let reach = crossings.iter().fold(0.0_f64, |m, c| m.max(-*c)) + radius * 0.25;
        if end {
            seat.start -= seat.along * reach;
        } else {
            seat.end += seat.along * reach;
        }
        trims.push((at - normal * (at - origin).dot(normal), normal));
        trimmed[slot] = true;
    }
    // A straight seat runs out. Where a convex edge ends at a face that
    // cuts the corner away obliquely (a floor giving way to a ramp under a
    // wall that carries on), the seat runs on along the host it stays
    // tangent to until that face has cut its whole section away, and the
    // cut trims it there. At a plain wall the cap already stands in the
    // wall's own plane and stays; at a step the seat stops where the edge
    // does. Where a seat asked to meet ends against a neighbouring blend
    // (a band tangent to one of the hosts at the end vertex), the material
    // past the end is that blend's own rounding, and the seat carries on
    // through it until the ball has left the solid; the cut then trims the
    // two bands against each other along their own intersection, which is
    // what two fillets meeting at a corner are.
    let mut probes = crate::support::LazyProbe::default();
    for (slot, end) in [true, false].into_iter().enumerate() {
        // A band trimmed by a face it runs into stops on that face.
        if trimmed[slot] {
            continue;
        }
        let at = if end { seat.start } else { seat.end };
        let outward = if end { -seat.along } else { seat.along };
        if let Some((index, mates)) = mates {
            // A mate leaving this vertex the way this edge arrives is a
            // chain mate: the caps stand flush and the melt joins the
            // blends.
            let chain = mates.iter().enumerate().any(|(i, mate)| {
                i != index
                    && mate.ends.iter().any(|(p, leaving)| {
                        p.distance(at) <= tol.confusion() * 1e3
                            && leaving.cross(outward).magnitude() <= 1e-2
                            && leaving.dot(outward) > 0.0
                    })
            });
            if chain {
                continue;
            }
            // A corner the corner tool closed, or a junction of the request
            // that meets the band here.
            if Mate::settled_at(mates, index, at, tol) || Mate::meets_at(mates, index, at, tol) {
                continue;
            }
        }
        let hosts = [&seat.faces[0], &seat.faces[1]];
        if !crate::marched::crease_terminates_at(model, solid, edge, hosts, at, tol)? {
            continue;
        }
        let through_blend =
            crate::marched::neighbour_blend_at(model, solid, hosts, at, seat.convex, radius, tol)?;
        // One edge at a time stops flush against a neighbouring blend: the
        // state the corner tool closes. A fill runs on only through one.
        if through_blend && mates.is_none() || !through_blend && !seat.convex {
            continue;
        }
        // The ball's centre line, a radius in from both hosts.
        let Some(inward) = seat_centre(&seat, radius, tol)? else {
            continue;
        };
        let centre_at = |s: f64| at + outward * s + inward;
        let step = radius / 8.0;
        let mut reach: Option<f64> = None;
        if through_blend {
            let Some(probe) = probes.get(model, solid, tol)? else {
                continue;
            };
            for k in 1..=32 {
                let s = step * f64::from(k);
                let inside = probe.holds(model, centre_at(s), tol)? == ogeom_algo::Containment::In;
                if inside != seat.convex {
                    reach = Some(s + radius * 0.25);
                    break;
                }
            }
        } else {
            // No blend follows: the band runs out across whatever face cuts
            // the corner away past the end, probed over its section there,
            // the crease, the two contacts, each leg's middle and the arc's.
            // A boss's cylinder is run out across only where the band runs
            // out as it does alone: where another band of the request ends
            // at the vertex with no junction there, the two meet as before,
            // the later running on through the earlier.
            let section = section_at(&seat, at, inward, radius);
            let boss = mates.is_none_or(|(index, mates)| Mate::alone_at(mates, index, at, tol));
            reach = crate::support::run_out_reach(
                model,
                solid,
                edge,
                hosts,
                (at, outward),
                &section,
                (radius, boss),
                &mut probes,
                tol,
            )?;
        }
        if std::env::var_os("OGEOM_DEBUG_RUNOUT").is_some() {
            eprintln!("RUNOUT straight seat end {end}: reach {reach:?}");
        }
        let Some(mut reach) = reach else {
            continue;
        };
        // The edge as asked for ends at its own vertex, where its blend
        // alone caps flush. Where the ball's centre is still in the material
        // a step past that vertex, the solid runs on there (a curved face
        // leaving the vertex obliquely), and a seat carried further would
        // round material no ball on this edge touches: it stops at the
        // vertex, as rounding this edge before the neighbour does.
        if let Some(mate) = mates.and_then(|(index, mates)| mates.get(index)) {
            let own = mate
                .ends
                .iter()
                .map(|(p, _)| *p)
                .min_by(|p, q| p.distance(at).total_cmp(&q.distance(at)));
            let to_own = own.map(|own| (own - at).dot(outward));
            if let Some(to_own) = to_own.filter(|&t| t > tol.confusion() && t < reach)
                && let Some(probe) = probes.get(model, solid, tol)?
                && probe.holds(model, centre_at(to_own + step), tol)? == ogeom_algo::Containment::In
            {
                reach = to_own;
            }
        }
        if end {
            seat.start -= seat.along * reach;
        } else {
            seat.end += seat.along * reach;
        }
    }
    seated_fillet(model, solid, &seat, radius, Some(edge), &trims, tol)
}

/// From a point of a straight seat's crease to the ball's centre in the
/// section square to the crease there; `None` where the faces meet too
/// sharply or too flat to seat a ball.
fn seat_centre(seat: &Seat, radius: f64, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let bisector = seat.leg(0, tol)? * sign + seat.leg(1, tol)? * sign;
    let m = bisector.magnitude();
    if m <= tol.angular() {
        return Ok(None);
    }
    let depth = seat.normals[0].dot(bisector / m).abs();
    if depth <= tol.angular() {
        return Ok(None);
    }
    Ok(Some(bisector / m * (radius / depth)))
}

/// Points of a straight seat's band section square to the crease at `at`,
/// the ball's centre `inward` of it: the section's corners (the crease and
/// the two contacts) first, then each leg's middle and the arc's.
fn section_at(seat: &Seat, at: Point, inward: Vector, radius: f64) -> [Point; 6] {
    let centre = at + inward;
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let touch = [0, 1].map(|i| centre + seat.normals[i] * (radius * sign));
    let arc = centre - inward * (radius / inward.magnitude());
    [
        at,
        touch[0],
        touch[1],
        at.midpoint(touch[0]),
        at.midpoint(touch[1]),
        arc,
    ]
}

/// The cylindrical blend of a seat, whatever found the seat.
///
/// A blend's construction cares about the seat (where the ball rolls and
/// which way the two planes face) and not at all about whether an edge of
/// the solid runs along it. An edge gives one; two faces that share nothing
/// give the same one through their planes' own intersection, and everything
/// from here down is common to both.
pub(crate) fn seated_fillet(
    model: &mut Model,
    solid: &Shape,
    seat: &Seat,
    radius: f64,
    edge: Option<&Shape>,
    trims: &[(Point, Vector)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    // Order the two faces so the blend arc sweeps positively about the edge
    // direction: the cylinder's parameterization and every arc below then
    // run from the first face's tangency line to the second's.
    let (first, second) = if seat
        .along
        .dot(seat.normals[0].cross(seat.normals[1]))
        .is_sign_positive()
    {
        (0, 1)
    } else {
        (1, 0)
    };
    let n1 = seat.normals[first];
    let n2 = seat.normals[second];
    // On a concave edge the legs mirror into the open dihedral, the ball
    // rolls there instead of in the material, and the contacts sit on the
    // faces' other side of the centre.
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(first, tol)? * sign;
    let b = seat.leg(second, tol)? * sign;

    // The rolling ball's centre line: along the bisector, at the distance
    // where the ball touches both planes.
    let bisector = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
        }
        u / m
    };
    let depth = (n1.dot(bisector)).abs();
    if depth <= tol.angular() {
        ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
    }
    let centre = seat.start + bisector * (radius / depth);
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2));
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }

    let travel = seat.end - seat.start;
    let length = travel.magnitude();
    let apex0 = seat.start;
    let apex1 = seat.end;
    let contact_a0 = centre + n1 * (radius * sign);
    let contact_b0 = centre + n2 * (radius * sign);
    let contact_a1 = contact_a0 + travel;
    let contact_b1 = contact_b0 + travel;

    let along_dir = Direction::new(seat.along, tol)?;
    let radial_dir = Direction::new(n1 * sign, tol)?;

    // An arc of the blend circle at height `h` along the edge, running from
    // the first face's contact vertex to the second's. Its frame is the
    // cylinder's frame translated up the axis, which is what makes its pcurve
    // on the cylinder a straight line in `(u, v)` with the arc's own
    // parameter.
    let arc = |model: &mut Model, h: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let frame = Frame::new(centre + seat.along * h, along_dir, radial_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };

    // The legs, coplanar with the solid's own faces: aligned with them on a
    // convex edge, where the melt keeps one copy, opposed on a concave one,
    // where the fuse cancels both.
    let leg_a = planar_face(
        model,
        &[apex0, contact_a0, contact_a1, apex1],
        n1 * sign,
        tol,
    )?;
    let leg_b = planar_face(
        model,
        &[apex0, contact_b0, contact_b1, apex1],
        n2 * sign,
        tol,
    )?;

    // The caps: two segments to the tangency points, closed by the arc.
    let cap = |model: &mut Model,
               apex: Point,
               ca: Point,
               cb: Point,
               h: f64,
               outward: Vector|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex, Direction::new(outward, tol)?);
        let reach = (apex.distance(ca).max(apex.distance(cb)) + radius).max(1.0) * 2.0;
        let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
        let apex_v = make_vertex(model, apex).shape;
        let ca_v = make_vertex(model, ca).shape;
        let cb_v = make_vertex(model, cb).shape;
        let edges = vec![
            segment_between(model, (&apex_v, apex), (&ca_v, ca), tol)?,
            arc(model, h, &ca_v, &cb_v)?,
            segment_between(model, (&cb_v, cb), (&apex_v, apex), tol)?,
        ];
        // Walked back where apex, `ca`, `cb` turn against `outward`, so the
        // cap keeps itself on the left.
        let edges: Vec<Shape> = if (ca - apex).cross(cb - apex).dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        face_from_edges(model, surface.into(), &edges, tol)
    };
    let cap0 = cap(model, apex0, contact_a0, contact_b0, 0.0, -seat.along)?;
    let cap1 = cap(model, apex1, contact_a1, contact_b1, length, seat.along)?;

    // The blend face itself. A cylinder's natural normal points away from its
    // axis (into this wedge, whose material lies between the cylinder and
    // the apex), so the face enters the shell reversed.
    let blend = {
        let frame = Frame::new(centre, along_dir, radial_dir, tol)?;
        let surface = CylinderSurface::new(Cylinder::new(frame, radius, tol)?, (0.0, length))?;
        let a0_v = make_vertex(model, contact_a0).shape;
        let b0_v = make_vertex(model, contact_b0).shape;
        let a1_v = make_vertex(model, contact_a1).shape;
        let b1_v = make_vertex(model, contact_b1).shape;
        let edges = vec![
            arc(model, 0.0, &a0_v, &b0_v)?,
            segment_between(model, (&b0_v, contact_b0), (&b1_v, contact_b1), tol)?,
            arc(model, length, &a1_v, &b1_v)?.reversed(),
            segment_between(model, (&a0_v, contact_a0), (&a1_v, contact_a1), tol)?.reversed(),
        ];
        face_from_edges(model, surface.into(), &edges, tol)?
    };

    let faces = [leg_a, leg_b, cap0, cap1, blend.reversed()];
    let mut wedge = wedge_solid(model, &faces, false, tol)?;
    // Wide and deep enough to hold all of the wedge behind each plane.
    let span = length + (centre.distance(apex0) + radius) * 4.0;
    for &plane in trims {
        wedge = trimmed_wedge(model, &wedge, plane, span, tol)?;
    }
    apply_wedge_solid(model, solid, edge, &wedge, !seat.convex, tol)
}

/// The circular-rim fillet: a ball rolling around a rim where a cylindrical
/// wall meets a perpendicular planar cap traces a torus.
///
/// Four seats, one parameterization. With `sigma` the wall's outward radial
/// sign and `tau` telling whether the wall extends away from the cap's
/// outward side, the tube's centre circle sits at radius `R − sigma tau r`,
/// lifted `tau r` against the cap's normal, and `tau` alone decides whether
/// the wedge subtracts (the external rim and the hole's rim, both convex) or
/// fuses (the boss base and the blind hole's floor, both concave). The wedge
/// is always the same three revolved faces: a band of the wall, an annulus
/// of the cap, and the quarter-tube between them.
fn revolved_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    let tube_rho = seat.sigma.mul_add(-(seat.tau * radius), seat.radius);
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre - seat.up * (seat.tau * radius);
    let flanks = revolved_flanks(model, &seat, radius, tube_rho, tol)?;

    // The blend: the quarter-tube between the tangencies, its natural normal
    // away from the tube's centre (into the wedge), so always reversed.
    let blend_band = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            seat.frame_at(tube_level, tol)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?.reversed()
    };

    let faces = [flanks.wall_band, flanks.annulus, blend_band];
    apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol)
}

/// Whether the blend along `edge` of `solid`, a circle, runs the circle's
/// whole turn: the edge is closed, or it is a piece of a rim the solid
/// holds whole on its revolved seat, or of a crease that closes on itself
/// on any other pair of hosts.
fn rounds_whole_turn(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (Curve::Circle(c), range) = edge_curve(model, edge, tol)? else {
        return Ok(false);
    };
    if ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b)) {
        return Ok(true);
    }
    match revolved_seat(model, solid, edge, &c, tol) {
        Ok(seat) => crate::support::runs_whole_turn(model, &seat, &c, range, tol),
        Err(OgeomError::Construction(_)) => crate::marched::crease_closes(model, solid, edge, tol),
        Err(e) => Err(e),
    }
}

/// The circle an edge runs on, where it runs on one.
fn rim_circle(model: &Model, edge: &Shape, tol: Tolerances) -> Option<ogeom_math::Circle> {
    match edge_curve(model, edge, tol).ok()?.0 {
        Curve::Circle(c) => Some(c.circle()),
        _ => None,
    }
}

/// Whether two circles are one: the same centre, radius and plane.
fn same_circle(a: &ogeom_math::Circle, b: &ogeom_math::Circle, tol: Tolerances) -> bool {
    let reach = tol.confusion() * 100.0;
    a.centre().distance(b.centre()) <= reach
        && (a.radius() - b.radius()).abs() <= reach
        && a.frame()
            .z()
            .vector()
            .cross(b.frame().z().vector())
            .magnitude()
            <= tol.angular()
}

/// The edge of `solid` standing where `edge` stood.
///
/// A boolean rebuilds every face it splits with fresh edges and its history
/// speaks of faces, not of them, so a chain's later edge is re-found by
/// geometry: the solid's edges whose own samples all lie on the sought
/// edge's curve. A trimmed survivor qualifies, and so does each piece an
/// earlier blend running out across the edge left of it; every one is a
/// part of the edge still there to blend.
///
/// A corner, though, and not merely a line. An earlier blend's end cap can
/// stand in one of the hosts' own planes and leave a fresh edge along the
/// very line the sought one runs on (an L-bracket's front blend caps in
/// the wall's plane, along the wall's own top line) and there the two
/// faces are tangent: a seam across one flat, no corner to round, and no
/// part of what was asked for. Nor is a corner between other surfaces: a
/// cap standing square to the line where a curved host leaves it meets
/// that host along the very line, a corner of the cap's and none of the
/// edge's. A piece lies between the surfaces of `hosts`, the faces the
/// edge lay between where it was asked for.
fn refind_edges(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    hosts: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    use ogeom_geom::Curve3d as _;
    use ogeom_topo::ShapeType;
    let (curve, range) = edge_curve(model, edge, tol)?;
    // The edge's own stretch of its curve, for a candidate to lie within:
    // the wall's bottom edge is two pieces either side of a scoop, on one
    // line, and asking for the left one must not blend the right.
    let within = |t: f64| -> bool {
        let slack = (range.1 - range.0).abs() * 1e-3;
        let t = if curve.is_periodic() {
            let (lo, hi) = curve.domain();
            let period = hi - lo;
            if period > 0.0 {
                range.0 + (t - range.0).rem_euclid(period)
            } else {
                t
            }
        } else {
            t
        };
        t >= range.0.min(range.1) - slack && t <= range.0.max(range.1) + slack
    };
    // A candidate lying on the edge's curve lies in its box: one that does
    // not reach the box is not projected at all.
    let reach = ogeom_algo::shape_bounds(model, edge, tol)?.expanded(tol.confusion() * 100.0);
    let mut matches: Vec<Shape> = Vec::new();
    for candidate in ogeom_topo::explore_unique(model, solid, ShapeType::Edge)? {
        if candidate.is_same(edge) {
            return Ok(vec![candidate]);
        }
        if !ogeom_algo::shape_bounds(model, &candidate, tol)?.intersects(&reach) {
            continue;
        }
        let Ok((c_curve, c_range)) = edge_curve(model, &candidate, tol) else {
            continue;
        };
        let mut all_on = true;
        for i in 0..3 {
            let t = c_range.0 + (c_range.1 - c_range.0) * f64::from(i) / 2.0;
            let p = c_curve.point_at(t, tol)?;
            let projected = ogeom_algo::project_on_curve(&curve, p, 64, tol)?;
            if projected.distance > tol.confusion() * 100.0
                || (i == 1 && !within(projected.parameter))
            {
                all_on = false;
                break;
            }
        }
        if !all_on || matches.iter().any(|m| m.is_same(&candidate)) {
            continue;
        }
        let mid = c_curve.point_at(f64::midpoint(c_range.0, c_range.1), tol)?;
        let mut normals: Vec<Vector> = Vec::new();
        for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
            let touches = ogeom_topo::explore_unique(model, &face, ShapeType::Edge)?
                .iter()
                .any(|e| crate::support::same_occurrence(model, e, &candidate, tol));
            if let (true, Some(normal)) = (
                touches,
                crate::marched::face_normal_near(model, &face, mid, tol)?,
            ) {
                normals.push(normal);
            }
        }
        // Two faces that agree on their normal there are one smooth
        // surface; anything else the seat machinery speaks to in its own
        // words.
        if matches!(normals.as_slice(), [first, second]
            if first.cross(*second).magnitude() <= 1e-2)
        {
            continue;
        }
        let mut host_normals: Vec<Vector> = Vec::with_capacity(hosts.len());
        for host in hosts {
            if let Some(normal) = crate::marched::face_normal_near(model, host, mid, tol)? {
                host_normals.push(normal);
            }
        }
        let along = |a: &[Vector], b: &[Vector]| {
            a.iter()
                .all(|n| b.iter().any(|m| n.cross(*m).magnitude() <= 1e-4))
        };
        if let ([_, _], [_, _]) = (normals.as_slice(), host_normals.as_slice())
            && !(along(&normals, &host_normals) && along(&host_normals, &normals))
        {
            continue;
        }
        matches.push(candidate);
    }
    if matches.is_empty() {
        ogeom_bail!(
            Construction,
            "an earlier blend in the chain consumed this edge; the chain's \
             members interfere"
        );
    }
    Ok(matches)
}

/// The rim fillet over an arc of the rim rather than its whole turn.
///
/// The same four-seat torus geometry as [`revolved_fillet`], restricted to
/// the arc's angular window and closed by two planar end caps in the
/// meridian half-planes at its ends. Every face is a chart rectangle or a
/// planar triangle, every edge an arc or a segment, every pcurve exact,
/// which is what lets a chain's neighbouring wedges melt at their shared
/// caps.
fn revolved_arc_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    let tube_rho = seat.sigma.mul_add(-(seat.tau * radius), seat.radius);
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre - seat.up * (seat.tau * radius);

    // The arc's window in the seat's own parameterization: the edge's stored
    // circle may be framed anyhow, so the ends are read back through the
    // seat frame and the window is the one that holds the arc's midpoint.
    let seat_frame = seat.frame_at(seat.centre, tol)?;
    let (curve, crange) = edge_curve(model, edge, tol)?;
    let angle_of = |t: f64| -> OgeomResult<f64> {
        let local = seat_frame.to_local(curve.point_at(t, tol)?);
        Ok(local.y.atan2(local.x))
    };
    let tau_turn = core::f64::consts::TAU;
    let raw0 = angle_of(crange.0)?;
    let raw1 = angle_of(crange.1)?;
    let mid = angle_of(f64::midpoint(crange.0, crange.1))?;
    let (theta0, theta1) = {
        let span = (raw1 - raw0).rem_euclid(tau_turn);
        let inside = (mid - raw0).rem_euclid(tau_turn) <= span;
        if inside {
            (raw0, raw0 + span)
        } else {
            (raw1, raw1 + (tau_turn - span))
        }
    };
    if theta1 - theta0 <= tol.angular() {
        ogeom_bail!(Construction, "the rim arc has no length to blend along");
    }
    let sweep = theta1 - theta0;

    // Corner points at each end: A on the rim, W down the wall, C in on the
    // cap; and the frames the arcs run on.
    let dir_at = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    // Every wedge chart is anchored at the arc's start: the window then
    // lives in (0, sweep), clear of the chart's own seam whatever part of
    // the rim the arc covers.
    let wedge_frame_at = |origin: Point| -> OgeomResult<Frame> {
        Frame::new(
            origin,
            Direction::new(seat.up, tol)?,
            Direction::new(dir_at(theta0), tol)?,
            tol,
        )
    };
    let apex_at = |theta: f64| seat.centre + dir_at(theta) * seat.radius;
    let wall_at = |theta: f64| tube_level + dir_at(theta) * seat.radius;
    let cap_at = |theta: f64| seat.centre + dir_at(theta) * tube_rho;

    // Six shared corner vertices, nine shared edges: a wire chains through
    // vertex objects, and neighbouring faces sew through shared edges.
    let va0 = make_vertex(model, apex_at(theta0)).shape;
    let va1 = make_vertex(model, apex_at(theta1)).shape;
    let vw0 = make_vertex(model, wall_at(theta0)).shape;
    let vw1 = make_vertex(model, wall_at(theta1)).shape;
    let vc0 = make_vertex(model, cap_at(theta0)).shape;
    let vc1 = make_vertex(model, cap_at(theta1)).shape;

    let ring_arc = |model: &mut Model,
                    origin: Point,
                    r: f64,
                    from: &Shape,
                    to: &Shape|
     -> OgeomResult<Shape> {
        let circle = Circle::new(wedge_frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };
    // The quarter meridian at `theta`, from the wall contact to the cap
    // contact: the tube's own circle, in the meridian plane.
    let quarter = |model: &mut Model, theta: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let dir = dir_at(theta);
        let centre_m = tube_level + dir * tube_rho;
        let x_dir = Direction::new(dir * (seat.sigma * seat.tau), tol)?;
        let z_dir = Direction::new(dir.cross(seat.up) * seat.sigma, tol)?;
        let frame = Frame::new(centre_m, z_dir, x_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(
            model,
            curve,
            (0.0, core::f64::consts::FRAC_PI_2),
            from,
            to,
            tol,
        )?
        .shape)
    };

    let arc_apex = ring_arc(model, seat.centre, seat.radius, &va0, &va1)?;
    let arc_wall = ring_arc(model, tube_level, seat.radius, &vw0, &vw1)?;
    let arc_cap = ring_arc(model, seat.centre, tube_rho, &vc0, &vc1)?;
    let ruling0 = segment_between(model, (&va0, apex_at(theta0)), (&vw0, wall_at(theta0)), tol)?;
    let ruling1 = segment_between(model, (&vw1, wall_at(theta1)), (&va1, apex_at(theta1)), tol)?;
    let radial0 = segment_between(model, (&vc0, cap_at(theta0)), (&va0, apex_at(theta0)), tol)?;
    let radial1 = segment_between(model, (&va1, apex_at(theta1)), (&vc1, cap_at(theta1)), tol)?;
    let quarter0 = quarter(model, theta0, &vw0, &vc0)?;
    let quarter1 = quarter(model, theta1, &vw1, &vc1)?;
    // Each face keeps itself on the left of its ring in its own chart: a
    // ring given turning the other way is walked back.
    let walked = |edges: Vec<Shape>, back: bool| -> Vec<Shape> {
        if back {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        }
    };

    // The wall patch: the solid wall's own chart between the contact ring
    // and the rim, over the window, its chart reaching past the rim on
    // both sides whichever side of the cap the wall stands on. The ring
    // runs the contact ring forward and the rim back, anticlockwise in the
    // chart where the contact ring stands below the rim.
    let wall_patch = {
        let below = radius + 1.0;
        let origin = seat.centre - seat.up * below;
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(wedge_frame_at(origin)?, seat.radius, tol)?,
            (0.0, 2.0 * below),
        )?
        .into();
        let edges = walked(
            vec![
                arc_wall.clone(),
                ruling1.clone(),
                arc_apex.reversed(),
                ruling0.clone(),
            ],
            seat.tau < 0.0,
        );
        let face = face_from_edges(model, surface, &edges, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            face.reversed()
        } else {
            face
        }
    };

    // The cap sector: the cap plane between the rim and the contact ring.
    // The ring runs the rim forward about `up`, which turns it about the
    // plane's normal `up * tau` where the rim is the wider ring and the cap
    // faces up, or the narrower and it faces down.
    let cap_sector = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let reach = (seat.radius + radius + 1.0) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let rim_wider = seat.radius > tube_rho;
        let edges = walked(
            vec![
                arc_apex.clone(),
                radial1.clone(),
                arc_cap.reversed(),
                radial0.clone(),
            ],
            rim_wider != (seat.tau > 0.0),
        );
        face_from_edges(model, surface, &edges, tol)?
    };

    // The blend: the quarter-tube patch over the window, reversed as the
    // full rim's band is: its natural normal points into the wedge. The
    // ring runs the wall's contact forward and the cap's back, anticlockwise
    // in the torus chart where the tube's angle grows from the wall's contact
    // to the cap's, which is where `sigma` is positive.
    let blend_patch = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            wedge_frame_at(tube_level)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        let edges = walked(
            vec![
                arc_wall.clone(),
                quarter1.clone(),
                arc_cap.reversed(),
                quarter0.reversed(),
            ],
            seat.sigma < 0.0,
        );
        face_from_edges(model, surface, &edges, tol)?.reversed()
    };

    // The end caps: planar triangles in the meridian half-planes, outward
    // along the rim's travel: behind it at the start, ahead at the end.
    // Each ring is walked back where the corners it passes in turn
    // against `outward`, so the cap keeps itself on the left: the quarter
    // arc bulges toward the rim corner and never past it, so the corners'
    // turn is the ring's.
    let tangent_at = |theta: f64| seat.up.cross(dir_at(theta));
    let end_cap = |model: &mut Model,
                   theta: f64,
                   outward: Vector,
                   edges: Vec<Shape>,
                   [p, q, r]: [Point; 3]|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex_at(theta), Direction::new(outward, tol)?);
        let reach = (seat.radius + radius + 1.0) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let edges = walked(edges, (q - p).cross(r - p).dot(outward) < 0.0);
        face_from_edges(model, surface, &edges, tol)
    };
    let cap_start = end_cap(
        model,
        theta0,
        -tangent_at(theta0),
        vec![ruling0.clone(), quarter0.clone(), radial0.clone()],
        [apex_at(theta0), wall_at(theta0), cap_at(theta0)],
    )?;
    let cap_end = end_cap(
        model,
        theta1,
        tangent_at(theta1),
        vec![ruling1.clone(), radial1.clone(), quarter1.reversed()],
        [wall_at(theta1), apex_at(theta1), cap_at(theta1)],
    )?;

    let faces = [wall_patch, cap_sector, blend_patch, cap_start, cap_end];
    apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol)
}
