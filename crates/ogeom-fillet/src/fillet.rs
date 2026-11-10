//! Constant-radius fillets, where the rolling ball's envelope is closed-form.
//!
//! A ball of constant radius rolling along a straight edge between two planar
//! faces stays tangent to both, and its envelope is exactly a cylinder; the
//! same ball rolling around the rim where a cylindrical or conical wall meets its
//! perpendicular cap traces a torus. These are the blends whose *new* surface
//! costs nothing, and both are built the same way: the chamfer's wedge with
//! the bevel exchanged for the envelope (legs running along the faces to the
//! tangency lines, and the envelope between them, every pcurve exact). The
//! boolean then does what it does for the chamfer: the legs melt into the
//! solid's own faces by same-domain resolution, and the envelope stays as
//! the blend.

use crate::law::{Profile, RadiusLaw};
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
/// and a circular rim where a cylindrical or conical wall meets a
/// perpendicular planar cap, where it traces a torus. In both, the result is the boolean
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
/// out of the solid, both cut the band off on the face that follows. Where
/// a convex edge's face is a planar strip narrower than the setback all
/// along, its far edge convex too and its ends on faces square to the edge
/// (a chamfer, a narrow top), the ball rolls between the planes either
/// side of the strip, and where the strip lies wholly in the corner that
/// ball takes off, the round is that ball's cylinder and the strip goes
/// with the corner. Where a concave edge's face is such a strip and its
/// far edge is straight and convex (the seam where two solids are fused
/// at a slant), the ball touches the other face and rests on that far
/// edge, and where it does so alone over the whole crease, clear of every
/// other face, the round is the band fitted through those balls' arcs.
/// Where it rests there over part of the crease, coming down on a face
/// either side or rolling round a corner of the strip onto its next edge,
/// the round is a chain of such bands meeting tangent where the ball
/// changes what it holds. One thing differs: where the edge ends against a neighbouring blend, this
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
/// material, or past it all along the edge where no strip is taken off and
/// no far edge holds the ball as above; naming why where the ball rests
/// on such a far edge but comes onto the face beyond it, rests on far
/// edges either side at once, or leaves a face or comes within its radius
/// of another along the crease.
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
    // the band would cut through the face. Past a strip's far side all
    // along, the ball may roll between the faces either side of the strip
    // and take it off whole.
    if let Err(refusal) = crate::support::ball_fits(model, solid, edge, radius, false, tol) {
        return consuming(model, solid, edge, radius, refusal, tol);
    }
    let edges = std::slice::from_ref(edge);
    or_rim_consumed(
        model,
        solid,
        edge,
        radius,
        |model| {
            with_fitted_rims(
                model,
                solid,
                edges,
                |model, rims| fillet_edge_meeting(model, solid, edge, radius, None, rims, tol),
                tol,
            )
        },
        tol,
    )
}

/// The round of `edge` that takes a narrow strip beside it off whole, where
/// [`crate::consumed::consumed_seat`] finds one, and `refusal` otherwise.
fn consuming(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    refusal: OgeomError,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !matches!(refusal, OgeomError::Construction(_)) {
        return Err(refusal);
    }
    let Some(seat) = crate::consumed::consumed_seat(model, solid, edge, radius, tol)? else {
        let rolled = match crate::rolled::rolled_fillet(model, solid, edge, radius, tol)? {
            Some(built) => Some(built),
            None => crate::consumed::consumed_slanted(model, solid, edge, radius, tol)?,
        };
        let Some(mut built) = rolled else {
            return match crate::consumed_rim::consumed_rim(model, solid, edge, radius, tol)? {
                Some(rim) => rim_consuming(model, solid, edge, rim, radius, tol),
                None => crate::consumed_marched::consumed_marched(model, solid, edge, radius, tol)?
                    .ok_or(refusal),
            };
        };
        crate::support::credit_new_faces(model, solid, edge, &mut built)?;
        return Ok(built);
    };
    let mut built = seated_fillet(model, solid, &seat, radius, Some(edge), &[], tol)?;
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(built)
}

/// The round of a rim on `rim`, the seat of the cap and wall either side of
/// a ring beside it that the round takes off whole.
fn rim_consuming(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: crate::support::RevolvedSeat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut built = revolved_fillet_on(model, solid, edge, rim, radius, tol)?;
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(built)
}

/// `round`, the rounding of the one `edge`, and where it refuses, the round
/// that takes a ring beside a rim off whole, where
/// [`crate::consumed_rim::consumed_rim`] finds one: a rim beside a small
/// round meets it tangentially, so no seat of the rim's own builds. What
/// the refused attempt widened is put back first; elsewhere the attempt
/// stands as it is.
fn or_rim_consumed(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    round: impl FnOnce(&mut Model) -> OgeomResult<Built>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mark = model.note_widened();
    match round(model) {
        Err(refusal @ OgeomError::Construction(_)) => {
            match crate::consumed_rim::consumed_rim(model, solid, edge, radius, tol)? {
                Some(rim) => {
                    model.undo_widened(mark);
                    rim_consuming(model, solid, edge, rim, radius, tol)
                }
                None => {
                    // The tolerances the refused attempt widened stand, as
                    // they do where the refusal stands.
                    model.keep_widened(mark);
                    crate::consumed_marched::consumed_marched(model, solid, edge, radius, tol)?
                        .ok_or(refusal)
                }
            }
        }
        // A strip beside a marched seat can leave the edge's own round
        // nothing to stand on short of a construction refusal.
        Err(refusal @ OgeomError::NotDone(_)) => {
            model.keep_widened(mark);
            crate::consumed_marched::consumed_marched(model, solid, edge, radius, tol)?
                .ok_or(refusal)
        }
        other => {
            model.keep_widened(mark);
            other
        }
    }
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
    rims: Rims,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut built = fillet_edge_seat(model, solid, edge, radius, mates, rims, tol)?;
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
    rims: Rims,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    match route(model, solid, edge, rims, tol)? {
        Route::Planar => planar_fillet(model, solid, edge, radius, mates, tol),
        Route::Ruled(seat) => {
            crate::ruled::ruled_fillet(model, solid, edge, radius, mates, seat, tol)
        }
        // A cone's ball whose tube would reach its own axis takes the march:
        // the torus would cross the axis.
        Route::Revolved(_, _, seat)
            if seat.slope != 0.0 && seat.ball(radius).rho - radius <= tol.confusion() =>
        {
            crate::marched::marched_fillet(model, solid, edge, radius, mates, tol)
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
                revolved_arc_fillet(model, solid, edge, &c, radius, mates, tol)
            }
        }
        Route::Marched => crate::marched::marched_fillet(model, solid, edge, radius, mates, tol),
    }
}

/// Whether an open rim arc's torus band, capped in the meridian half-plane
/// through its end `at` (the crease leaving along `ahead` past it), would
/// stop short of a face the arc ends on: the crease terminates there, no
/// mate of the request ends there, and a planar face there stands oblique
/// to the crease, its normal neither square to the crease's direction (a
/// wall continuing the crease tangentially, which meets the cap square) nor
/// along it (a wall square to the crease, which is the cap). Past such an
/// end the ball still rolls through material between the meridian cap and
/// the face while its contact on the cap stands on it. The face's plane,
/// facing back over the arc, where it does.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn rim_end_stops_short(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &crate::support::RevolvedSeat,
    mates: Option<(usize, &[Mate])>,
    at: Point,
    ahead: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Vector)>> {
    let hosts = [&seat.cap_face, &seat.wall_face];
    let met = mates.is_some_and(|(index, mates)| {
        mates.iter().enumerate().any(|(i, mate)| {
            i != index
                && mate
                    .ends
                    .iter()
                    .any(|(p, _)| p.distance(at) <= tol.confusion() * 1e3)
        })
    });
    if met || !crate::marched::crease_terminates_at(model, solid, edge, hosts, at, tol)? {
        return Ok(None);
    }
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if hosts.iter().any(|h| h.is_same(&face)) {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            continue;
        };
        use ogeom_geom::Transformable as _;
        let SurfaceGeometry::Plane(placed) =
            SurfaceGeometry::Plane(*plane).transformed(&face.transform(model.datums())?, tol)?
        else {
            continue;
        };
        let plane = placed.plane();
        if plane.signed_distance_to(at).abs() > tol.confusion() * 10.0 {
            continue;
        }
        let meets = explore_unique(model, &face, ShapeType::Vertex)?
            .iter()
            .any(|v| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .zip(v.transform(model.datums()).ok())
                    .is_some_and(|(d, place)| {
                        place.apply(d.point).distance(at) <= tol.confusion() * 1e3
                    })
            });
        let normal = plane.normal().vector();
        let lean = normal.dot(ahead).abs();
        if meets && lean > 1e-6 && lean < 1.0 - 1e-6 {
            // Facing back over the arc: the side of the face the arc lies on.
            let back = if normal.dot(ahead) > 0.0 {
                -normal
            } else {
                normal
            };
            return Ok(Some((plane.origin(), back)));
        }
    }
    Ok(None)
}

/// Which construction rounds a rim on a planar cap and a coaxial cone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Rims {
    /// The exact torus band.
    Exact,
    /// The marched band, fitted through the ball's arcs.
    Fitted,
}

/// The construction that rounds an edge, chosen by its curve and the two
/// faces it lies between.
enum Route {
    /// A straight edge between two planes.
    Planar,
    /// A straight edge between planes or drums parallel to it.
    Ruled(crate::ruled::RuledSeat),
    /// A circle between a planar cap and its coaxial wall, a cylinder or
    /// a cone, with the edge's range on it.
    Revolved(CircleCurve, (f64, f64), crate::support::RevolvedSeat),
    /// Any other seat: the march.
    Marched,
}

/// How `edge` of `solid` is rounded. A rim on a planar cap and its
/// coaxial wall is the revolved blend's, exact, but where `rims` asks a
/// cone's for the march; a circle on any other pair of hosts (a bore down
/// a ball's axis, a ring's own rim) is a seat like any other for the march.
fn route(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rims: Rims,
    tol: Tolerances,
) -> OgeomResult<Route> {
    let (curve, crange) = edge_curve(model, edge, tol)?;
    Ok(match curve {
        Curve::Line(_) if crate::support::hosts_planar(model, solid, edge, tol)? => Route::Planar,
        Curve::Line(_) => match crate::ruled::ruled_seat(model, solid, edge, tol)? {
            Some(seat) => Route::Ruled(seat),
            None => Route::Marched,
        },
        Curve::Circle(c) => match revolved_seat(model, solid, edge, &c, tol) {
            Ok(seat) if seat.slope != 0.0 && rims == Rims::Fitted => Route::Marched,
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
/// cell; across a convex turn of a curved crease both stop in it where the
/// bands' own curvature keeps their sections there within that fit. Either
/// way the two are built from the solid as it stands and joined across the
/// plane into one wedge, each band's faces credited to its own edge. Where
/// the end sections differ (a drafted wall beside a square one) and both
/// creases are straight, across a re-entrant turn the corner takes what
/// lies above both bands carried on, bounded by the boolean's own
/// intersection of the two, where the material holds it all. Across a
/// convex turn of two straight creases, whatever their sections, each band
/// runs out past the vertex as it does alone and the cuts take their union
/// (walls leaning alike with their feet out each reach past the bisecting
/// plane under the other's band, which a meeting there would leave
/// standing), as they do where one of the two is a drum's rim (a D-shaped
/// block's arc beside its chord, and past where they meet flush, a rib's
/// spline wall running into a boss's rim). Across a re-entrant turn one
/// of the two may be a drum's rim, upright or drafted, the material within
/// the drum, where the
/// straight band carried on stands nowhere off the drum before the rim's
/// cap. Such
/// junctions stand beside flush ones in a loop or chain, a band crossed at
/// both ends cut in two halfway along, and two bands meeting at both their
/// ends meet at each.
/// Elsewhere (a wider convex turn of two curved bands, or of a curved band
/// other than a drum's rim beside one of another kind)
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
/// its re-entrant edge) no ball rolls round both. There the request is
/// rounded in one order whatever its own, and a vertex of three planar
/// faces and three straight edges, all of them asked, is closed as
/// [`round_vertex_setback`](crate::round_vertex_setback) closes a corner:
/// each band stopped one of its own widths past the least setback at which
/// its section clears the others, and one patch G1 to every band end and
/// face, the same solid in any order.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`fillet_edge`] per edge, and additionally if an earlier blend
/// consumed or split a later edge: a chain whose members interfere is
/// refused rather than guessed at. Where a convex and a concave edge of the
/// request meet at a vertex that is not three planar faces and three
/// straight edges all asked, the request is refused by name in any order:
/// with an edge of it left sharp, the opening between the band ends would
/// cross that edge, and no patch is tangent to both its faces at one
/// point. A request of two or more edges also
/// refuses a ball that sets back past the far side of a face anywhere
/// along an edge that shares a face or a vertex with another of the
/// request, since each band's contacts are checked against the others'
/// and a contact off its face holds none to check. An edge sharing
/// neither is held to what [`fillet_edge`] holds it to.
pub fn fillet_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let radii = vec![radius; edges.len()];
    ogeom_algo::on_own_nodes(model, solid, edges, |model, solid, parts| {
        fillet_edges_meeting(model, solid, parts, &radii, tol)
    })
}

/// [`fillet_edges`] with a radius per edge.
///
/// With every radius the same, exactly [`fillet_edges`] at that radius.
/// Each band is built at its own radius and checked against the others as
/// [`fillet_edges`] checks them. Two edges of one radius meet as there.
/// Two edges of different radii turning at a vertex are rounded in turn,
/// the later band running on through the earlier and the cut trimming
/// the two against each other along their own intersection. A tangent
/// chain whose radius steps at a junction is rounded likewise, the later
/// band's end standing against the earlier band.
///
/// # Errors
///
/// As [`fillet_edges`], and as [`round_vertex_setback`](crate::round_vertex_setback)
/// where three or more of the edges meeting at a vertex do not share one
/// radius. No ball rolls round such a corner: each band is stopped one of
/// its own widths past the least setback at which its section clears every
/// other band there, and one patch closes the opening G1 to every band end
/// and face. [`OgeomError::NotDone`] where that patch misses the bands' fit
/// target, and [`OgeomError::Construction`] where the corner is not a
/// convex planar one every edge of which the request rounds. An edge asked
/// twice at two radii is refused likewise.
pub fn fillet_edges_with(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let shapes: Vec<Shape> = edges.iter().map(|(edge, _)| edge.clone()).collect();
    let radii: Vec<f64> = edges.iter().map(|&(_, radius)| radius).collect();
    for (i, (edge, radius)) in edges.iter().enumerate() {
        if !radius.is_finite() || *radius <= tol.confusion() {
            ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
        }
        if edges[..i]
            .iter()
            .any(|(other, r)| other.is_same(edge) && r != radius)
        {
            ogeom_bail!(
                Construction,
                "an edge asked at two radii, {radius} and another, rounds neither"
            );
        }
    }
    ogeom_algo::on_own_nodes(model, solid, &shapes, |model, solid, parts| {
        fillet_edges_meeting(model, solid, parts, &radii, tol)
    })
}

/// [`fillet_edges_with`] on a shape every node below which is its own, its
/// cone rims rounded exactly where that builds and fitted where it does
/// not, as [`with_fitted_rims`] says.
fn fillet_edges_meeting(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    tol: Tolerances,
) -> OgeomResult<Built> {
    // Where a convex and a concave band meet at a vertex, the bands are
    // built in one order whatever the request's, and the corner closed
    // after.
    let request: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii.iter().copied()).collect();
    let signed = crate::setback::signed_corners(model, solid, &request, tol)?;
    let (edges, radii): (Vec<Shape>, Vec<f64>) = if signed.is_empty() {
        (edges.to_vec(), radii.to_vec())
    } else {
        let mut sorted = request.clone();
        sorted.sort_by_key(|(edge, _)| edge.node());
        sorted.into_iter().unzip()
    };
    let (edges, radii) = (edges.as_slice(), radii.as_slice());
    let round = |model: &mut Model| {
        with_fitted_rims(
            model,
            solid,
            edges,
            |model, rims| fillet_edges_meeting_with(model, solid, edges, radii, rims, tol),
            tol,
        )
    };
    let mut mixed = mixed_corners(model, edges, radii)?;
    mixed.retain(|(vertex, _)| !signed.iter().any(|c| c.vertex.is_same(vertex)));
    // A request of one edge is that edge alone, a ring beside a rim taken
    // off whole as there.
    let mut built = match (edges, radii) {
        ([edge], [radius]) => or_rim_consumed(model, solid, edge, *radius, round, tol)?,
        _ => round(model)?,
    };
    for (vertex, request) in &mixed {
        built = crate::setback::close_mixed_corner(model, built, vertex, request, tol)?;
    }
    for corner in &signed {
        built = crate::setback::close_signed_corner(model, built, corner, tol)?;
    }
    Ok(built)
}

/// A vertex of the request and the request's edges there with their radii.
type Corner = (Shape, Vec<(Shape, f64)>);

/// The vertices where three or more edges of the request meet at
/// different radii, each with the request's edges there and their radii.
fn mixed_corners(model: &Model, edges: &[Shape], radii: &[f64]) -> OgeomResult<Vec<Corner>> {
    let mut corners: Vec<Corner> = Vec::new();
    for (edge, &radius) in edges.iter().zip(radii) {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        for v in [a, b] {
            if let Some(slot) = corners.iter_mut().find(|(held, _)| held.is_same(&v)) {
                slot.1.push((edge.clone(), radius));
            } else {
                corners.push((v, vec![(edge.clone(), radius)]));
            }
        }
    }
    corners.retain(|(_, at)| at.len() >= 3 && at.iter().any(|&(_, r)| r != at[0].1));
    Ok(corners)
}

/// Round with every cone rim's band an exact torus band, and where that
/// refuses, again from the solid as it was given (its tolerances put back)
/// with those bands fitted through the ball's arcs, as the march builds
/// them. The fitted band stands within its fit of the torus, and the
/// boolean resolves some meetings of it that it does not of the exact one,
/// so an exact rim never turns a rounding into a refusal.
fn with_fitted_rims(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    mut round: impl FnMut(&mut Model, Rims) -> OgeomResult<Built>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mark = model.note_widened();
    let exact = round(model, Rims::Exact);
    if matches!(
        exact,
        Ok(_) | Err(OgeomError::Cancelled | OgeomError::Dangling(_))
    ) || !has_cone_rim(model, solid, edges, tol)?
    {
        model.keep_widened(mark);
        return exact;
    }
    model.undo_widened(mark);
    round(model, Rims::Fitted)
}

/// Whether any of `edges` is a rim of `solid` on a planar cap and a
/// coaxial cone.
fn has_cone_rim(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<bool> {
    for edge in edges {
        if let Route::Revolved(_, _, seat) = route(model, solid, edge, Rims::Exact, tol)?
            && seat.slope != 0.0
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// [`fillet_edges_meeting`] with cone rims rounded as `rims` says, its
/// bands meeting at the junctions where two of its creases turn. Where
/// that meeting refuses, the request is rounded again from the solid as it
/// was given (its tolerances and surfaces put back), every band capped or
/// run on at its ends as it is without junctions: a junction never turns a
/// rounding into a refusal.
fn fillet_edges_meeting_with(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    rims: Rims,
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
        return fillet_edges_own(model, solid, edges, radii, None, rims, None, tol);
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
        match fillet_edges_own(
            model,
            solid,
            edges,
            radii,
            Some(&mut taken),
            rims,
            None,
            tol,
        ) {
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
    fillet_edges_own(model, solid, edges, radii, None, rims, None, tol)
}

/// The checks a request's bands pass before any is built: each ball sits
/// on both faces of its edge, and no band overlaps another across a face
/// they share. Which edges, sharing nothing with the rest of the request,
/// had their ball rest on the far edge of a narrow face beside them, to
/// be rounded as alone, comes back.
///
/// # Errors
///
/// [`OgeomError::Construction`] where a ball does not sit on its faces or
/// two bands overlap, naming the edges refused.
pub(crate) fn request_fits(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    tol: Tolerances,
) -> Result<Vec<bool>, crate::support::Blamed> {
    use crate::support::Blamed;
    let blame = |i: usize| {
        move |error: OgeomError| Blamed {
            edges: vec![i],
            error,
        }
    };
    // The ball touches each face a setback back from the edge; set back
    // past a face's far side, the band would cut through the face, and set
    // back into another edge's band, the two balls overlap.
    let mut contacts = Vec::with_capacity(edges.len());
    // An edge sharing no face and no vertex with another of the request
    // stands alone: no corner or junction takes its ends, and its band
    // is built as it is alone.
    let mut touched: Vec<(Vec<ogeom_topo::TShapeId>, Vec<ogeom_topo::TShapeId>)> =
        Vec::with_capacity(edges.len());
    for edge in edges {
        let faces = ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?
            .iter()
            .map(Shape::node)
            .collect();
        let ends = ogeom_algo::edge_vertices(model, edge)?
            .map(|(a, b)| vec![a.node(), b.node()])
            .unwrap_or_default();
        touched.push((faces, ends));
    }
    let alone = |i: usize| {
        let (faces, ends) = &touched[i];
        touched.iter().enumerate().all(|(j, (f, e))| {
            j == i || (!f.iter().any(|x| faces.contains(x)) && !e.iter().any(|x| ends.contains(x)))
        })
    };
    // Edges sharing nothing with the rest of the request whose ball rests
    // on the far edge of a narrow face beside them, rounded as alone.
    let mut rolled = vec![false; edges.len()];
    crate::support::keeping_face_meshes(|| -> Result<(), Blamed> {
        // A band checked against no other is held to what one edge alone
        // is.
        for (i, edge) in edges.iter().enumerate() {
            match crate::support::ball_fits(model, solid, edge, radii[i], !alone(i), tol) {
                Ok(found) => contacts.push(found),
                Err(refusal @ OgeomError::Construction(_)) if edges.len() > 1 && alone(i) => {
                    if crate::rolled::rolled_run(model, solid, edge, radii[i], tol)
                        .map_err(blame(i))?
                        .is_none()
                    {
                        return Err(blame(i)(refusal));
                    }
                    rolled[i] = true;
                    contacts.push(Vec::new());
                }
                Err(other) => return Err(blame(i)(other)),
            }
        }
        // A band lies within its setback of its edge, its contacts its
        // farthest points. A band standing alone whose ball runs off a face
        // is cut off on the face across that face's side, which may be one
        // another band stands on; where the two edges stand within their
        // setbacks' sum, it is checked against the others as a neighbour's
        // is.
        let reach =
            |found: &[crate::support::Contact]| found.iter().map(|c| c.setback).reduce(f64::max);
        let mut boxes = Vec::with_capacity(edges.len());
        for edge in edges {
            boxes.push(ogeom_algo::shape_bounds(model, edge, tol)?);
        }
        let reaches: Vec<Option<f64>> = contacts.iter().map(|found| reach(found)).collect();
        for (i, edge) in edges.iter().enumerate() {
            let Some(own) = reaches[i].filter(|_| alone(i)) else {
                continue;
            };
            let near = reaches.iter().enumerate().any(|(j, other)| {
                j != i
                    && other.is_some_and(|other| boxes[i].distance_to_box(&boxes[j]) < own + other)
            });
            if near {
                contacts[i] = crate::support::ball_fits(model, solid, edge, radii[i], true, tol)
                    .map_err(blame(i))?;
            }
        }
        Ok(())
    })?;
    crate::support::bands_clear(model, edges, &contacts, tol)?;
    Ok(rolled)
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
#[allow(clippy::too_many_arguments, reason = "one request's settings")]
fn fillet_edges_own(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    meeting: Option<&mut Taken>,
    cones: Rims,
    laws: Option<&[Profile]>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
    let rolled = match request_fits(model, solid, edges, radii, tol) {
        Ok(rolled) => rolled,
        Err(crate::support::Blamed { error: refusal, .. }) if laws.is_some() => {
            return Err(refusal);
        }
        // A request of one edge is that edge alone, a strip beside it taken
        // off whole as there.
        Err(crate::support::Blamed { error: refusal, .. }) => {
            return match edges {
                [edge] if meeting.is_none() => {
                    consuming(model, solid, edge, radii[0], refusal, tol)
                }
                _ => Err(refusal),
            };
        }
    };
    use ogeom_geom::Curve3d as _;
    // The vertices three or more of the chain's edges meet at: the corners
    // the ball rolls round, closed by the corner tool.
    // Each holds the radii of the edges ending there.
    let mut corners: Vec<(Shape, Vec<f64>)> = Vec::new();
    for (edge, &radius) in edges.iter().zip(radii) {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        for v in [a, b] {
            if let Some(slot) = corners.iter_mut().find(|(held, _)| held.is_same(&v)) {
                slot.1.push(radius);
            } else {
                corners.push((v, vec![radius]));
            }
        }
    }
    corners.retain(|(_, at)| at.len() >= 3);
    if laws.is_some() && (!corners.is_empty() || rolled.iter().any(|&r| r)) {
        ogeom_bail!(
            Construction,
            "a variable radius over a chain three of whose edges meet at a \
             vertex, or one of whose balls rests on a narrow face's far edge, \
             is still owed; see docs/CAPABILITIES.md, fillet.variable"
        );
    }
    // No ball rolls round a corner whose bands differ in radius: the
    // corner tool's ball would touch a band of another radius nowhere
    // along its rim. Such a corner is left to the bands here and closed by
    // the setback patch once they stand ([`mixed_corners`]).
    corners.retain(|(_, at)| at.iter().all(|&r| r == at[0]));
    // A corner the tool rounds, or the solid unchanged where it does not
    // speak the corner, which keeps the bands' caps.
    let round =
        |model: &mut Model, built: Built, vertex: &Shape, radius: f64| -> OgeomResult<Built> {
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
    // Corners of different radii are cut one after another.
    let together = if corners.len() > 1 && corners.iter().all(|(_, at)| at[0] == corners[0].1[0]) {
        let vertices: Vec<Shape> = corners.iter().map(|(v, _)| v.clone()).collect();
        crate::corner::round_vertices_together(model, solid, &vertices, corners[0].1[0], tol)?
    } else {
        None
    };
    let one_by_one: &[(Shape, Vec<f64>)] = if together.is_some() { &[] } else { &corners };
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
    for (vertex, at) in one_by_one {
        let start = built
            .take()
            .unwrap_or_else(|| Built::from_nothing(solid.clone()));
        let before = start.shape.clone();
        let after = round(model, start, vertex, at[0])?;
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
        let mut found = crate::junction::junctions(model, solid, edges, &given, &mates, radii, tol)
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
                route(model, solid, &edges[i], cones, tol)?,
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
                let halves = with_corners(model, &mut set_aside, &junctions, &b.shape, false, tol)?;
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
                    if rolled[index] {
                        let Some(mut built) = crate::rolled::rolled_fillet(
                            model,
                            &current,
                            target,
                            radii[index],
                            tol,
                        )?
                        else {
                            ogeom_bail!(
                                Construction,
                                "the ball no longer rests on the far edge beside an edge \
                                 of the chain"
                            );
                        };
                        crate::support::credit_new_faces(model, &current, target, &mut built)?;
                        return Ok(built);
                    }
                    if let Some(laws) = laws {
                        return law_piece(
                            model,
                            &current,
                            &edges[index],
                            target,
                            &laws[index],
                            (index, &mates),
                            tol,
                        );
                    }
                    fillet_edge_meeting(
                        model,
                        &current,
                        target,
                        radii[index],
                        Some((index, &mates)),
                        cones,
                        tol,
                    )
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
                                    radii[index] * 2.0,
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
    let halves = with_corners(model, &mut set_aside, &junctions, &built.shape, false, tol)?;
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
/// `law` says the bands are under a radius law, their sections changing up
/// to the vertex, as [`crate::junction::corner`] takes them.
fn with_corners(
    model: &mut Model,
    set_aside: &mut Vec<(Vec<usize>, Shape, crate::support::Wedge)>,
    junctions: &[crate::junction::Junction],
    solid: &Shape,
    law: bool,
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
            law,
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
/// twice. A cut is [`ogeom_bool::cut_seated`] and a fill
/// [`ogeom_bool::fuse_seated`], which trim the faces straight bands, rims
/// and the pieces of rims a tangent chain's arcs round sit on, where nothing
/// else reaches into them, and straight bands running at a block's corner
/// into a band an earlier group rounded.
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
    let law = Profile::new(
        &RadiusLaw::Linear {
            start: start_radius,
            end: end_radius,
        },
        false,
        tol,
    )?;
    if let Some(built) = consumed_law(model, solid, edge, &law, false, tol)? {
        return Ok(built);
    }
    law_ball_on_faces(model, solid, edge, &law, Some(false), tol)?;
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

/// Round a convex or concave edge, or a closed chain of curved edges, with
/// a blend whose radius follows `law` along `chain`, taken as one
/// parameter in arc length.
///
/// At arc length `s` along the chain the blend's section, in the plane
/// square to the edge there, is the arc of radius `law(s)` tangent to both
/// faces: the rolling ball of that radius seated in the section. The
/// chain starts where its first edge does as that edge is oriented (a
/// reversed edge starts at its stored curve's far end) and runs the way
/// that edge does.
///
/// Rounded:
///
/// - a straight edge between two planes, for any law. The linear law is
///   [`fillet_edge_variable`]'s blend; under [`RadiusLaw::Points`] the
///   envelope is again exact: a cubic B-spline along the edge, the
///   rational quadratic arc across it, the control net the edge's
///   Greville points pushed out by the law's own control radii.
/// - a closed loop of curved edges between the same two faces, given as
///   one closed edge or as the arcs it is split into, end to end round
///   the loop: a drum's rim, a cut drum's oblique rim, a branch's seam. The
///   march carries the law's radius station by station round the whole
///   loop, and the band is fitted through each station's own arc.
/// - an open curved edge whose ball runs out past its ends (a half drum's
///   rim): the law is held at its end radius past each end, so the ball
///   runs out there at that one radius, as a constant round's does.
/// - an open chain of arcs of one curve between the same two faces, end
///   to end (a crease split at a vertex): the one seat over the window of
///   that curve its arcs cover.
/// - a straight edge along a drum's ruling (a flat or a wall meeting a
///   drum, or two drums): marched with the law's radius station by station.
/// - a seat pinched where its hosts turn tangent at an end (two equal
///   drums crossing): each station solved at the law's radius for its
///   place, the band pinching to the pole.
/// - a chain of edges of different seats continuing each other end to
///   end, open or closed (the top loop of a rounded-rectangle pad: its
///   sides between two planes, its corners' arcs between the top and a
///   drum). The run is walked from the first edge the way it runs; each
///   edge takes the law over its own stretch of the run, its band stopped
///   in its end sections, and the bands' wedges are joined across each
///   junction, where one ball touches the same two points from either
///   side, and applied in one boolean.
/// - a chain of straight edges between planes turning where they meet on
///   one shared planar host (a block's top loop, an L-block's). Past a
///   turning vertex each band goes on at the law's radius there, as the
///   law holds level past an edge's end: across a convex turn it runs out
///   as the constant band of that radius alone does and the cuts take the
///   union; across a re-entrant one it is carried to the plane bisecting
///   the turn, where the two are joined. A loop the tangent junctions join
///   round to one convex turn (a pad with one vertical left sharp) is
///   opened at a straight band between two of them, built in two halves.
/// - a chain turning where marched bands meet (a drum's ruling turning
///   into its rim), met in rounds as [`fillet_edges`] meets it: each band
///   a later round builds stands over what the earlier left of its edge,
///   under its own window of the law, and runs on through the earlier
///   band at the law's radius there.
/// - a straight edge alone beside a narrow strip the ball takes off whole
///   at every radius of the law: the exact envelope over the seat between
///   the faces either side of the strip.
///
/// A law holding one radius all along is [`fillet_edges`]' round.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`] where the law states no usable radius
/// (no points, a fraction off the chain or out of order, a radius not a
/// positive length), where a linear law changes radius round a closed
/// chain, where `chain` is empty, a chain of more than one edge whose
/// edges do not run end to end, a turn where a marched band meets a
/// straight one between planes, a loop joined
/// round a convex turn with no straight band between two of its junctions
/// to open it at, a ball resting
/// on a narrow face's far edge or two bands overlapping across a face
/// they share at an edge's largest radius, a
/// straight edge between faces that are neither two planes nor a drum and
/// a face along its ruling (only two planes inside a chain of different
/// seats), where two bands' end sections at a junction part by more than
/// a marched band's fit, where a law is given over one arc of a
/// crease that closes through other edges, or where the ball at the law's
/// largest radius does not fit the edge's faces; and as [`fillet_edge`]
/// for the seat itself.
pub fn fillet_edges_variable(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if chain.is_empty() {
        ogeom_bail!(
            Construction,
            "a variable radius over no edges rounds nothing"
        );
    }
    ogeom_algo::on_own_nodes(model, solid, chain, |model, solid, parts| {
        fillet_chain_law_own(model, solid, parts, law, tol)
    })
}

/// Whether `chain` runs end to end and back to its start, each edge
/// between the same two faces of `solid` as the first, none twice.
fn closed_chain(model: &Model, solid: &Shape, chain: &[Shape]) -> OgeomResult<bool> {
    if !runs_between_one_pair(model, solid, chain)? {
        return Ok(false);
    }
    let (Some(first), Some(last)) = (
        ogeom_algo::edge_vertices(model, &chain[0])?,
        ogeom_algo::edge_vertices(model, &chain[chain.len() - 1])?,
    ) else {
        return Ok(false);
    };
    Ok(last.1.is_same(&first.0))
}

/// The window of one curve an open `chain` covers, where its edges run
/// end to end between the same two faces of `solid` and are each a piece
/// of that curve (the same geometry), one after the other along it: in the curve's own
/// parameters, a piece past a periodic curve's seam taken a period on.
/// `None` for any other chain.
fn one_curve_span(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    if !runs_between_one_pair(model, solid, chain)? {
        return Ok(None);
    }
    let mut pieces = Vec::with_capacity(chain.len());
    for edge in chain {
        let Some(ogeom_topo::EdgeRepr::Curve3d {
            curve,
            location,
            range,
        }) = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.curve3d())
            .cloned()
        else {
            return Ok(None);
        };
        if edge.location() != &ogeom_topo::Location::default() {
            return Ok(None);
        }
        pieces.push((curve, location, range, edge.orientation()));
    }
    // One curve: the same stored geometry, whether or not under one handle
    // (a boolean that splits a crease gives each piece its own copy).
    let (curve, location, _, orientation) = pieces[0].clone();
    let Some(geometry) = model.geometry().curve(curve) else {
        return Ok(None);
    };
    if pieces.iter().any(|(c, l, _, o)| {
        *l != location || *o != orientation || model.geometry().curve(*c) != Some(geometry)
    }) {
        return Ok(None);
    }
    let period = {
        use ogeom_geom::Curve3d as _;
        geometry.is_periodic().then(|| {
            let (lo, hi) = geometry.domain();
            hi - lo
        })
    };
    // Along the curve the way the chain runs: forward pieces end where the
    // next begins, reversed ones begin where the next ends.
    let forward = orientation != ogeom_topo::Orientation::Reversed;
    let ordered: Vec<(f64, f64)> = if forward {
        pieces.iter().map(|p| p.2).collect()
    } else {
        pieces.iter().rev().map(|p| p.2).collect()
    };
    let reach = tol.parametric();
    let mut span = ordered[0];
    for &(lo, hi) in &ordered[1..] {
        let shift = match period {
            Some(period) if (lo + period - span.1).abs() <= reach => period,
            _ => 0.0,
        };
        if (lo + shift - span.1).abs() > reach {
            return Ok(None);
        }
        span.1 = hi + shift;
    }
    Ok(Some(span))
}

/// Whether `chain` runs end to end, each edge between the same two faces of
/// `solid` as the first, none twice.
fn runs_between_one_pair(model: &Model, solid: &Shape, chain: &[Shape]) -> OgeomResult<bool> {
    let faces_of = |edge: &Shape| -> OgeomResult<Vec<Shape>> {
        let mut faces: Vec<Shape> = Vec::new();
        for f in ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)? {
            if !faces.iter().any(|g| g.is_same(&f)) {
                faces.push(f);
            }
        }
        Ok(faces)
    };
    let hosts = faces_of(&chain[0])?;
    let mut ends = Vec::with_capacity(chain.len());
    for (i, edge) in chain.iter().enumerate() {
        if chain[..i].iter().any(|e| e.is_same(edge)) {
            return Ok(false);
        }
        let faces = faces_of(edge)?;
        let same =
            faces.len() == hosts.len() && faces.iter().all(|f| hosts.iter().any(|h| h.is_same(f)));
        if !same {
            return Ok(false);
        }
        let Some(pair) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(false);
        };
        ends.push(pair);
    }
    Ok((1..ends.len()).all(|i| ends[i - 1].1.is_same(&ends[i].0)))
}

/// [`fillet_edges_variable`] on a shape every node below which is its own.
fn fillet_chain_law_own(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let edge = &chain[0];
    let closed = closed_chain(model, solid, chain)?;
    // An open chain of arcs of one curve (a crease the boolean split at a
    // vertex) is rounded as the one seat it is, over the window of that
    // curve its arcs cover end to end.
    let span = if chain.len() > 1 && !closed {
        let Some(span) = one_curve_span(model, solid, chain, tol)? else {
            return chain_law(model, solid, chain, law, tol);
        };
        Some(span)
    } else {
        None
    };
    let profile = Profile::new(law, closed, tol)?;
    if profile.level(tol) {
        let radii = vec![profile.largest(); chain.len()];
        return fillet_edges_meeting(model, solid, chain, &radii, tol);
    }
    let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
    if chain.len() == 1
        && let RadiusLaw::Linear { start, end } = law
    {
        let (from, to) = if reversed {
            (*end, *start)
        } else {
            (*start, *end)
        };
        return fillet_edge_variable_own(model, solid, edge, from, to, tol);
    }
    // A straight edge alone: between two planes the envelope is exact;
    // along a drum's ruling (a wall meeting a drum, or two drums) every
    // section changes with the radius, and the march carries it.
    let mut planar = false;
    for e in chain {
        let (curve, _) = edge_curve(model, e, tol)?;
        if matches!(curve, Curve::Line(_)) {
            let seat = if chain.len() == 1 {
                Some(route(model, solid, e, Rims::Exact, tol)?)
            } else {
                None
            };
            planar = matches!(seat, Some(Route::Planar));
            if !planar && !matches!(seat, Some(Route::Ruled(_))) {
                ogeom_bail!(
                    Construction,
                    "a variable radius on a straight edge other than one \
                     between two planes or along a drum's ruling, alone, is \
                     still owed; see docs/CAPABILITIES.md, fillet.variable"
                );
            }
        }
    }
    // A narrow strip beside a straight edge alone, taken off whole at every
    // radius of the law, as the constant round takes it off at one.
    if planar
        && chain.len() == 1
        && let Some(built) = consumed_law(model, solid, edge, &profile, reversed, tol)?
    {
        return Ok(built);
    }
    // Each ball sits on both faces of its edge all along the law.
    let stations = (chain.len() == 1 && span.is_none()).then_some(reversed);
    if closed && chain.len() > 1 && !law_loop_fits(model, solid, chain, &profile, tol)? {
        for e in chain {
            law_ball_on_faces(model, solid, e, &profile, None, tol)?;
        }
    } else if !(closed && chain.len() > 1) {
        for e in chain {
            law_ball_on_faces(model, solid, e, &profile, stations, tol)?;
        }
    }
    let mut built = if planar {
        straight_law_fillet(model, solid, edge, &profile, reversed, (0.0, 1.0), tol)?
    } else {
        crate::marched::marched_fillet_law(model, solid, edge, &profile, span, None, tol)?
    };
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    // The chain's other arcs went with the first: the one seat over them.
    if span.is_some() {
        for other in &chain[1..] {
            built.history.delete(other);
        }
    }
    Ok(built)
}

/// The band over `piece`, what an earlier round of a chain left of the
/// chain's edge `whole`, under `law`: the law over `whole` from its own
/// start the way it runs, read over the stretch of it the piece covers.
///
/// A straight piece between planes is the exact envelope of that stretch;
/// where the constant band of the radius at an end runs on past it
/// (through a blend an earlier round built, as `fillet_edges` runs a
/// later band through an earlier one), what of that band stands past the
/// end is joined on, the law held level there at the end's radius. Any
/// other piece is marched under the law with the request's mates, which
/// run it out past an end as the constant band would.
fn law_piece(
    model: &mut Model,
    solid: &Shape,
    whole: &Shape,
    piece: &Shape,
    law: &Profile,
    (index, mates): (usize, &[Mate]),
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let (curve, range) = edge_curve(model, whole, tol)?;
    let ruler = crate::law::Ruler::new(curve.clone(), range, false, tol)?;
    let whole_reversed = whole.orientation() == ogeom_topo::Orientation::Reversed;
    // The fraction of the whole edge, from its own start the way it runs,
    // at a point of it.
    let fraction = |p: Point| -> OgeomResult<f64> {
        let t = ogeom_algo::project_on_curve(&curve, p, 64, tol)?.parameter;
        let s = ruler.fraction(t, tol)?.0;
        Ok(if whole_reversed { 1.0 - s } else { s })
    };
    let (piece_curve, (lo, hi)) = edge_curve(model, piece, tol)?;
    let reversed = piece.orientation() == ogeom_topo::Orientation::Reversed;
    let (start, end) = if reversed {
        (
            piece_curve.point_at(hi, tol)?,
            piece_curve.point_at(lo, tol)?,
        )
    } else {
        (
            piece_curve.point_at(lo, tol)?,
            piece_curve.point_at(hi, tol)?,
        )
    };
    let (from, to) = (fraction(start)?, fraction(end)?);
    if (to - from).abs() <= tol.parametric() {
        ogeom_bail!(
            Construction,
            "a piece of a chain's edge under a law has no length"
        );
    }
    let own = law.window(from, to);
    let planar = matches!(piece_curve, Curve::Line(_))
        && matches!(route(model, solid, piece, Rims::Exact, tol)?, Route::Planar);
    if !planar {
        return crate::marched::marched_fillet_law(
            model,
            solid,
            piece,
            &own,
            None,
            Some((index, mates)),
            tol,
        );
    }
    let (built, wedges) = crate::support::collecting_wedges(|| {
        straight_law_fillet(model, solid, piece, &own, reversed, (0.0, 1.0), tol)
    });
    built?;
    let [band] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|w| {
        ogeom_core::ogeom_err!(
            Construction,
            "a piece of a chain's edge under a law built {} wedges, not one",
            w.len()
        )
    })?;
    let mut wedge = band.solid;
    let gap = crate::marched::band_fit_target(tol) * 2.0;
    let seat = planar_seat(model, solid, piece, tol)?;
    for (at, outward, radius) in [
        (
            seat.start,
            -seat.along,
            if reversed {
                own.at(1.0).0
            } else {
                own.at(0.0).0
            },
        ),
        (
            seat.end,
            seat.along,
            if reversed {
                own.at(0.0).0
            } else {
                own.at(1.0).0
            },
        ),
    ] {
        let (alone, wedges) = crate::support::collecting_wedges(|| {
            planar_fillet(model, solid, piece, radius, Some((index, mates)), tol)
        });
        alone?;
        let Some(alone) = wedges.into_iter().next() else {
            continue;
        };
        // A straight band's sections stand square to its crease, so its
        // vertices reach as far past the end as any point of it.
        let mut past = f64::NEG_INFINITY;
        for vertex in explore_unique(model, &alone.solid, ShapeType::Vertex)? {
            if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                let p = vertex.transform(model.datums())?.apply(data.point);
                past = past.max((p - at).dot(outward));
            }
        }
        if past <= gap {
            continue;
        }
        let span = ogeom_algo::shape_bounds(model, &alone.solid, tol)?.diagonal() * 2.0;
        let run_out = crate::support::trimmed_wedge(model, &alone.solid, (at, outward), span, tol)?;
        wedge = crate::junction::joined(model, [&wedge, &run_out], (at, outward), at, gap, tol)?.0;
    }
    crate::support::apply_wedge_solid(model, solid, Some(piece), &wedge, band.additive, tol)
}

/// The straight edge's band under `law` where a narrow strip beside it is
/// taken off whole at every radius of the law, as the constant round
/// takes it off at one ([`crate::consumed::consumed_seat`]): the exact
/// envelope over the seat between the faces either side of the strip.
/// The constant round takes it off where its ball does not sit on the
/// edge's faces and the seat beyond the strip holds; each of that seat's
/// conditions grows with the radius or is linear in it, and a ball that
/// does not sit on the strip at one radius sits on it at no larger one, so
/// where the law's least and largest radius both take the strip off,
/// every radius between does. `None` where neither does; a law taking it
/// off at one and not the other is refused by name.
fn consumed_law(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    law: &Profile,
    reversed: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    // The constant round takes the strip off only where its ball does not
    // sit on the edge's faces: at a radius the strip holds, it rounds the
    // edge as any other.
    let consumes = |radius: f64| -> OgeomResult<Option<crate::support::Seat>> {
        match crate::support::ball_fits(model, solid, edge, radius, false, tol) {
            Ok(_) => Ok(None),
            Err(OgeomError::Construction(_)) => {
                crate::consumed::consumed_seat(model, solid, edge, radius, tol)
            }
            Err(e) => Err(e),
        }
    };
    let (least, largest) = (consumes(law.least())?, consumes(law.largest())?);
    let seat = match (least, largest) {
        (Some(_), Some(seat)) => seat,
        (None, None) => return Ok(None),
        _ => ogeom_bail!(
            Construction,
            "a variable radius whose ball sits on a narrow face beside the \
             edge at some radii of the law and takes it off whole at others \
             is still owed: the constant round changes construction between \
             them; see docs/CAPABILITIES.md, fillet.variable"
        ),
    };
    let mut built = straight_law_on(model, solid, edge, seat, law, reversed, (0.0, 1.0), tol)?;
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(Some(built))
}

/// The law's radius at a fraction of `edge`'s curve range, start to end:
/// the law laid by arc length from where the edge starts, the way it runs
/// (from the curve's end where `reversed`), as the band is built under it.
fn law_on_edge<'a>(
    model: &Model,
    edge: &Shape,
    law: &'a Profile,
    reversed: bool,
    tol: Tolerances,
) -> OgeomResult<Box<dyn Fn(f64) -> f64 + 'a>> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let flip = move |s: f64| if reversed { 1.0 - s } else { s };
    if matches!(curve, Curve::Line(_)) {
        return Ok(Box::new(move |f: f64| law.at(flip(f)).0));
    }
    let ruler = crate::law::Ruler::new(curve, range, false, tol)?;
    Ok(Box::new(move |f: f64| {
        let t = (range.1 - range.0).mul_add(f, range.0);
        let s = ruler.fraction(t, tol).map_or(f, |(s, _)| s);
        law.at(flip(s)).0
    }))
}

/// Whether the ball under `law` sits on its faces round a closed loop of
/// arcs of one curve, the one seat the march walks: each arc read station
/// by station at the law's radius there, as [`law_ball_on_faces`] reads one
/// edge alone, the law laid by arc length round the first arc's curve from
/// where that arc starts, the way it runs, as the march lays it.
/// `Ok(false)` where the first arc's curve does not close on itself, for
/// the loop to be read as [`law_ball_on_faces`] reads any other.
fn law_loop_fits(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    law: &Profile,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Curve3d as _;
    let (stored, (lo, hi)) = edge_curve(model, &chain[0], tol)?;
    if !stored.is_periodic() {
        return Ok(false);
    }
    let ruler = crate::law::Ruler::new(stored.clone(), stored.domain(), true, tol)?;
    let (start, heading) = if chain[0].orientation() == ogeom_topo::Orientation::Reversed {
        (stored.point_at(hi, tol)?, -stored.d1_at(hi, tol)?)
    } else {
        (stored.point_at(lo, tol)?, stored.d1_at(lo, tol)?)
    };
    let at = ogeom_algo::project_on_curve(&stored, start, 256, tol)?.parameter;
    let against = stored.d1_at(at, tol)?.dot(heading) < 0.0;
    let origin = ruler.fraction(at, tol)?.0;
    let fractions: Vec<f64> = (1..LAW_STATIONS)
        .map(|k| f64::from(k) / f64::from(LAW_STATIONS))
        .collect();
    for edge in chain {
        let (curve, (a, b)) = edge_curve(model, edge, tol)?;
        let radius_at = |f: f64| -> f64 {
            let read = || -> OgeomResult<f64> {
                let p = curve.point_at((b - a).mul_add(f, a), tol)?;
                let t = ogeom_algo::project_on_curve(&stored, p, 64, tol)?.parameter;
                let s = ruler.fraction(t, tol)?.0 - origin;
                Ok(law.at((if against { -s } else { s }).rem_euclid(1.0)).0)
            };
            read().unwrap_or_else(|_| law.largest())
        };
        crate::support::ball_fits_by(model, solid, edge, &radius_at, &fractions, false, tol)?;
    }
    Ok(true)
}

/// Whether the balls of a chain under a law sit on their faces, and no
/// band overlaps another across a face they share.
///
/// A straight edge between two planes is read station by station at its
/// own law's radius as [`law_ball_on_faces`] reads one alone, its contacts
/// carried past a face's side where its ball runs off a narrow face
/// ([`crate::support::ball_fits_by`] with `carry`), as a request's band
/// sharing a vertex with another is read; any other edge is read as
/// [`law_ball_on_faces`] reads it, every contact on its own face at the
/// law's least and largest radius. The contacts of every station are then
/// checked against the other edges' bands
/// ([`crate::support::bands_clear`]).
fn law_chain_fits(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    laws: &[&Profile],
    tol: Tolerances,
) -> OgeomResult<()> {
    let fractions: Vec<f64> = (1..LAW_STATIONS)
        .map(|k| f64::from(k) / f64::from(LAW_STATIONS))
        .collect();
    let mut contacts = Vec::with_capacity(chain.len());
    for (edge, law) in chain.iter().zip(laws) {
        let (curve, _) = edge_curve(model, edge, tol)?;
        let planar = matches!(curve, Curve::Line(_))
            && matches!(route(model, solid, edge, Rims::Exact, tol)?, Route::Planar);
        if planar {
            let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let radius_at = law_on_edge(model, edge, law, reversed, tol)?;
            contacts.push(crate::support::ball_fits_by(
                model,
                solid,
                edge,
                &*radius_at,
                &fractions,
                true,
                tol,
            )?);
        } else {
            law_ball_on_faces(model, solid, edge, law, None, tol)?;
            contacts.push(crate::support::ball_fits(
                model,
                solid,
                edge,
                law.largest(),
                true,
                tol,
            )?);
        }
    }
    crate::support::bands_clear(model, chain, &contacts, tol).map_err(|b| b.error)
}

/// How many pieces a straight edge under a law is read in for whether
/// its ball sits on both faces: a station at each inner end.
const LAW_STATIONS: u32 = 32;

/// Whether the ball under `law` sits on both faces of `edge` all along it.
///
/// Where `stations` is given (a straight edge between two planes rounded
/// alone, the law laid from the curve's start, or from its end where
/// `reversed`), the ball at each station has the law's radius there and is
/// read as the constant round reads its own
/// ([`crate::support::ball_fits_by`]): a contact on its face or on any face
/// of the solid holds; one past its face outside the solid is cut off
/// there, the band running on over the open as the constant round's does;
/// one standing in the material refuses, and so does a ball held at no
/// station. Where the law is level that is the constant round's own
/// reading at that radius.
///
/// Elsewhere [`crate::support::ball_fits`] is read at the law's least and
/// largest radius and every contact must stand on its own face: a chain's
/// bands and a marched band are not cut off where a contact runs off.
fn law_ball_on_faces(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    law: &Profile,
    stations: Option<bool>,
    tol: Tolerances,
) -> OgeomResult<()> {
    if let Some(reversed) = stations {
        let fractions: Vec<f64> = (1..LAW_STATIONS)
            .map(|k| f64::from(k) / f64::from(LAW_STATIONS))
            .collect();
        let radius_at = law_on_edge(model, edge, law, reversed, tol)?;
        crate::support::ball_fits_by(model, solid, edge, &*radius_at, &fractions, false, tol)?;
        return Ok(());
    }
    for radius in [law.least(), law.largest()] {
        let deflection = ogeom_mesh::Deflection {
            chord: (radius * 1e-3).max(tol.confusion() * 1e2),
            ..ogeom_mesh::Deflection::default()
        };
        for contact in crate::support::ball_fits(model, solid, edge, radius, false, tol)? {
            let on =
                ogeom_algo::classify_on_face(model, &contact.face, contact.at, deflection, tol)?;
            if on == ogeom_algo::Containment::Out {
                ogeom_bail!(
                    Construction,
                    "a variable radius whose ball at radius {radius} runs off a \
                     face beside the edge (a narrow face's far edge) is still \
                     owed; see docs/CAPABILITIES.md, fillet.variable"
                );
            }
        }
    }
    Ok(())
}

/// The request's edges by their ends, as a band asks of its mates: each
/// end's point and the unit tangent leaving it along the edge.
fn chain_mates(model: &Model, chain: &[Shape], tol: Tolerances) -> OgeomResult<Vec<Mate>> {
    use ogeom_geom::Curve3d as _;
    chain
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
            Ok(Mate {
                ends: [end_of(range.0, 1.0)?, end_of(range.1, -1.0)?],
                settled: [false, false],
                meets: [false, false],
                alone: [false, false],
            })
        })
        .collect()
}

/// A chain of edges of different seats running end to end, rounded under
/// `law` as one run: edges continuing each other (the top loop of a
/// rounded-rectangle pad: its sides between two planes, its corners' arcs
/// between the top and a drum), or straight edges between planes turning
/// where they meet on one shared planar host (a block's top loop).
///
/// The run is walked from the first edge the way it runs, on through the
/// edges sharing its ends, and back from its start where it is open; the
/// law's fraction is the run's arc length from where the walk starts, and
/// the run closes where it comes back to its own start. Each edge takes
/// the law over its own stretch of the run ([`Profile::window`], the law
/// itself there), so the ball at a junction is one radius from either
/// side. Each band is built against the solid as given and stopped in its
/// end sections: a straight edge's are square to it by construction, and a
/// marched band stops flush where a mate continues its crease. Where the
/// creases continue each other and the hosts either side are tangent, one
/// ball touches the same two points from both sides, so the two end
/// sections are one; the bands' wedges are joined across each junction's
/// plane into one wedge (sewn within the fit a marched band is held to,
/// which refuses ends that part by more) and applied in one boolean.
/// Where they turn, each band goes on past the vertex at the law's radius
/// there, which the law holds level past an edge's end: the constant band
/// of that radius, read as `fillet_edges` reads it, run out across a
/// convex turn and carried to the plane bisecting a re-entrant one.
fn chain_law(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mates = chain_mates(model, chain, tol)?;
    // Each edge's ends the way it runs: its oriented start, then its end.
    let runs: Vec<[Point; 2]> = chain
        .iter()
        .zip(&mates)
        .map(|(e, m)| {
            let [(a, _), (b, _)] = m.ends;
            if e.orientation() == ogeom_topo::Orientation::Reversed {
                [b, a]
            } else {
                [a, b]
            }
        })
        .collect();
    let near = |p: Point, q: Point| p.distance(q) <= tol.confusion() * 1e3;
    // The walk: each step an edge and whether it runs with the walk.
    let mut walk: Vec<(usize, bool)> = vec![(0, true)];
    let mut used = vec![false; chain.len()];
    used[0] = true;
    let mut head = runs[0][1];
    while let Some(k) =
        (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, head)))
    {
        let with = near(runs[k][0], head);
        used[k] = true;
        head = if with { runs[k][1] } else { runs[k][0] };
        walk.push((k, with));
    }
    let closed = near(head, runs[0][0]);
    if !closed {
        let mut tail = runs[0][0];
        while let Some(k) =
            (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, tail)))
        {
            let with = near(runs[k][1], tail);
            used[k] = true;
            tail = if with { runs[k][0] } else { runs[k][1] };
            walk.insert(0, (k, with));
        }
    }
    if used.iter().any(|u| !u) {
        ogeom_bail!(
            Construction,
            "a variable radius over a chain whose edges do not run end to end \
             rounds nothing; see docs/CAPABILITIES.md, fillet.variable"
        );
    }
    let tangent = crate::junction::tangent_junctions(chain, &mates, tol);
    let profile = Profile::new(law, closed, tol)?;
    if profile.level(tol) {
        let radii = vec![profile.largest(); chain.len()];
        return fillet_edges_meeting(model, solid, chain, &radii, tol);
    }
    // The run's arc length, edge by edge along the walk.
    let mut lengths = Vec::with_capacity(walk.len());
    for &(k, _) in &walk {
        let (curve, range) = edge_curve(model, &chain[k], tol)?;
        lengths.push(crate::law::Ruler::new(curve, range, false, tol)?.length());
    }
    let total: f64 = lengths.iter().sum();
    // Each meeting of the run: the two edges, where they meet, and the
    // law's radius there, one radius from either side.
    let mut meetings: Vec<([usize; 2], Point, f64)> = Vec::new();
    let mut before = 0.0;
    for (step, &(k, with)) in walk.iter().enumerate() {
        before += lengths[step];
        let Some(&(next, _)) = walk
            .get(step + 1)
            .or(if closed { walk.first() } else { None })
        else {
            continue;
        };
        let at = if with { runs[k][1] } else { runs[k][0] };
        meetings.push(([k, next], at, profile.at(before / total).0));
    }
    // A meeting where the creases continue each other is a tangent
    // junction; any other is the junction two bands of the meeting's
    // radius make there, read as `fillet_edges` reads it.
    let mut junctions: Vec<crate::junction::Junction> = Vec::new();
    // Whether a turn is met in rounds, as `fillet_edges` meets a curved
    // crease's band and the band beside it: the one built first, the other
    // against what it left, running on through it.
    let mut in_rounds = false;
    for &([i, j], at, radius) in &meetings {
        if let Some(found) = tangent
            .iter()
            .find(|t| t.slot(i).is_some() && t.slot(j).is_some() && near(t.at, at))
        {
            junctions.push(found.clone());
            continue;
        }
        let hosts: Vec<Vec<Shape>> = [i, j]
            .iter()
            .map(|&k| -> OgeomResult<Vec<Shape>> {
                let mut faces: Vec<Shape> = Vec::new();
                for f in ogeom_topo::ancestors_of(model, solid, &chain[k], ShapeType::Face)? {
                    if !faces.iter().any(|g| g.is_same(&f)) {
                        faces.push(f);
                    }
                }
                Ok(faces)
            })
            .collect::<OgeomResult<_>>()?;
        let turning = crate::junction::junctions(
            model,
            solid,
            &[chain[i].clone(), chain[j].clone()],
            &hosts,
            &[mates[i], mates[j]],
            &[radius, radius],
            tol,
        )?
        .into_iter()
        .find(|t| near(t.at, at));
        let planar = |model: &Model, k: usize| -> OgeomResult<bool> {
            Ok(
                edge_curve(model, &chain[k], tol).is_ok_and(|(c, _)| matches!(c, Curve::Line(_)))
                    && matches!(
                        route(model, solid, &chain[k], Rims::Exact, tol)?,
                        Route::Planar
                    ),
            )
        };
        let straight = planar(model, i)? && planar(model, j)?;
        // Met in rounds, a marched band and a straight one between planes
        // are refused: under a law the marched band has no exact form (its
        // law runs by arc length, which no rational parameter of a circle
        // does), and its fitted end meets the straight band beside it in
        // what the boolean does not melt, whichever is built first and
        // whether it stops in the wall it ends on or runs out past it. The
        // constant round's band there is an exact torus band.
        let mixed = planar(model, i)? != planar(model, j)?;
        let in_turn = |mixed: bool| -> OgeomResult<()> {
            if mixed {
                ogeom_bail!(
                    Construction,
                    "a variable radius over a chain turning where a marched \
                     band meets a straight one between planes is still owed: \
                     under a law the marched band is fitted, and its end meets \
                     the straight band in what the boolean does not melt; see \
                     docs/CAPABILITIES.md, fillet.variable"
                );
            }
            Ok(())
        };
        let Some(mut turning) = turning else {
            // No junction on one shared planar host: the two bands meet as
            // `fillet_edges`' do, in rounds.
            in_turn(mixed)?;
            in_rounds = true;
            continue;
        };
        turning.pair = [i, j];
        if !straight && turning.crossed && !turning.reentrant {
            in_turn(mixed)?;
            in_rounds = true;
            continue;
        }
        let convex = turning.crossed && !turning.reentrant;
        let reentrant = turning.reentrant;
        if !straight || !(convex || reentrant) {
            ogeom_bail!(
                Construction,
                "a variable radius over a chain turning where it meets is \
                 rounded between straight creases across a convex or a \
                 re-entrant turn, and a curved crease across a convex turn; \
                 any other turn is still owed; see docs/CAPABILITIES.md, \
                 fillet.variable"
            );
        }
        junctions.push(turning);
    }
    // Each edge's law over its own stretch of the run, laid from its own
    // start the way it runs.
    let mut owns: Vec<Profile> = Vec::with_capacity(walk.len());
    let mut before = 0.0;
    for (&(_, with), length) in walk.iter().zip(&lengths) {
        let (from, to) = (before / total, (before + length) / total);
        before += length;
        owns.push(if with {
            profile.window(from, to)
        } else {
            profile.window(to, from)
        });
    }
    // A chain met in rounds goes to `fillet_edges`' own rounds, each edge
    // under its own window of the law; its other meetings must then be
    // turns too, met the same way.
    if in_rounds {
        if !junctions.is_empty() {
            ogeom_bail!(
                Construction,
                "a variable radius over a chain that turns where a curved \
                 crease meets and also joins bands across another junction is \
                 still owed; see docs/CAPABILITIES.md, fillet.variable"
            );
        }
        let mut laws: Vec<Profile> = owns.clone();
        let mut radii: Vec<f64> = vec![0.0; chain.len()];
        for (&(k, _), own) in walk.iter().zip(&owns) {
            laws[k] = own.clone();
            radii[k] = own.largest();
        }
        for (edge, own) in chain.iter().zip(&laws) {
            law_ball_on_faces(model, solid, edge, own, None, tol)?;
        }
        return fillet_edges_own(
            model,
            solid,
            chain,
            &radii,
            None,
            Rims::Exact,
            Some(&laws),
            tol,
        );
    }
    // The request's checks at each edge's largest radius, where its ball
    // sets back furthest: each ball on its faces, and no band over another
    // across a face they share. A ball resting on a narrow face's far
    // edge, which `fillet_edges` takes off whole, is refused here.
    {
        let mut edges: Vec<Shape> = Vec::with_capacity(walk.len());
        let mut radii: Vec<f64> = Vec::with_capacity(walk.len());
        for (&(k, _), own) in walk.iter().zip(&owns) {
            edges.push(chain[k].clone());
            radii.push(own.largest());
        }
        let rests = request_fits(model, solid, &edges, &radii, tol).map_err(|b| b.error)?;
        if rests.iter().any(|&r| r) {
            ogeom_bail!(
                Construction,
                "a variable radius over a chain one of whose balls rests on \
                 the far edge of a narrow face beside it is still owed; see \
                 docs/CAPABILITIES.md, fillet.variable"
            );
        }
    }
    // Bands the junctions join into one wedge, by the group each edge's
    // band falls in. A group holding both bands of a convex turn (a loop of
    // tangent junctions closing on one sharp corner) would be one wedge
    // overlapping itself there; one straight band of it between two
    // junctions it joins is built in two halves instead, the halves'
    // shared section the same arc, so the group opens into two wedges that
    // meet flush there and overlap only at the corner, where the cuts take
    // their union.
    let mut group: Vec<usize> = (0..chain.len()).collect();
    for junction in junctions.iter().filter(|j| !j.crossed) {
        let [a, b] = junction.pair;
        let (from, to) = (group[a], group[b]);
        for g in &mut group {
            if *g == from {
                *g = to;
            }
        }
    }
    let mut halved: Vec<usize> = Vec::new();
    for junction in junctions.iter().filter(|j| j.crossed && !j.reentrant) {
        let [a, b] = junction.pair;
        if group[a] != group[b] || halved.iter().any(|&h| group[h] == group[a]) {
            continue;
        }
        let between = (0..chain.len()).find(|&m| {
            group[m] == group[a]
                && junctions
                    .iter()
                    .filter(|j| !j.crossed && j.slot(m).is_some())
                    .count()
                    == 2
                && edge_curve(model, &chain[m], tol).is_ok_and(|(c, _)| matches!(c, Curve::Line(_)))
        });
        let Some(m) = between else {
            ogeom_bail!(
                Construction,
                "a variable radius over a loop joined end to end round a convex \
                 turn with no straight edge between two of its junctions to \
                 open it at is still owed; see docs/CAPABILITIES.md, \
                 fillet.variable"
            );
        };
        halved.push(m);
    }
    {
        let mut laws: Vec<Option<&Profile>> = vec![None; chain.len()];
        for (&(k, _), own) in walk.iter().zip(&owns) {
            laws[k] = Some(own);
        }
        let laws: Vec<&Profile> = laws.into_iter().flatten().collect();
        if laws.len() != chain.len() {
            ogeom_bail!(Invariant, "a chain's walk takes every edge once");
        }
        law_chain_fits(model, solid, chain, &laws, tol)?;
    }
    let mut set_aside: Vec<(Vec<usize>, Shape, crate::support::Wedge)> = Vec::new();
    for ((&(k, _), own), _) in walk.iter().zip(&owns).zip(&lengths) {
        let edge = &chain[k];
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let (curve, _) = edge_curve(model, edge, tol)?;
        if halved.contains(&k) {
            // Each half of the edge by its stored parameter, its law laid
            // from its own start the way the edge runs.
            for (part, (from, to)) in [((0.0, 0.5), (0.0, 0.5)), ((0.5, 1.0), (0.5, 1.0))] {
                let law = if reversed {
                    own.window(1.0 - to, 1.0 - from)
                } else {
                    own.window(from, to)
                };
                let (built, wedges) = crate::support::collecting_wedges(|| {
                    straight_law_fillet(model, solid, edge, &law, reversed, part, tol)
                });
                built?;
                let [wedge] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|w| {
                    ogeom_core::ogeom_err!(
                        Construction,
                        "a half band of a chain under a law built {} wedges, not one",
                        w.len()
                    )
                })?;
                set_aside.push((vec![k], edge.clone(), wedge));
            }
            continue;
        }
        let (built, wedges) = crate::support::collecting_wedges(|| -> OgeomResult<Built> {
            if matches!(curve, Curve::Line(_)) {
                if !matches!(route(model, solid, edge, Rims::Exact, tol)?, Route::Planar) {
                    ogeom_bail!(
                        Construction,
                        "a variable radius on a straight edge of a chain other \
                         than one between two planes is still owed; see \
                         docs/CAPABILITIES.md, fillet.variable"
                    );
                }
                straight_law_fillet(model, solid, edge, own, reversed, (0.0, 1.0), tol)
            } else {
                crate::marched::marched_fillet_law(
                    model,
                    solid,
                    edge,
                    own,
                    None,
                    Some((k, &mates)),
                    tol,
                )
            }
        });
        built?;
        let [wedge] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|w| {
            ogeom_core::ogeom_err!(
                Construction,
                "a band of a chain under a law built {} wedges, not one",
                w.len()
            )
        })?;
        set_aside.push((vec![k], edge.clone(), wedge));
    }
    // At a turning junction each band goes on past the vertex at the
    // meeting's radius, where the law holds level past an edge's end: the
    // constant band of that radius, which is what the junction is read for.
    // Across a convex turn each band runs out as the constant band alone
    // does and the cuts take their union; across a re-entrant one each is
    // carried on to the plane the two meet in, where they are joined.
    let gap = crate::marched::band_fit_target(tol) * 2.0;
    let run_on = run_on_past_turns(
        model,
        solid,
        chain,
        (&walk, &owns, &halved),
        (&junctions, &meetings),
        &mut set_aside,
        tol,
    )?;
    for (junction, &(_, _, radius)) in junctions.iter().zip(&meetings) {
        if junction.tangent {
            continue;
        }
        for slot in 0..2 {
            let index = junction.pair[slot];
            let edge = chain[index].clone();
            let Some(entry) = set_aside
                .iter_mut()
                .find(|(held, _, _)| held.contains(&index))
            else {
                continue;
            };
            let outward = junction.outward[slot];
            // Across a re-entrant turn whose sections differ, the corner
            // joins the two bands, each followed under the law to the
            // vertex and carried on level past it.
            if junction.crossed && junction.reentrant {
                continue;
            }
            if junction.reentrant {
                let Some(cap) =
                    crate::junction::cap_at(model, &entry.2.solid, junction.at, outward, tol)?
                else {
                    ogeom_bail!(
                        Construction,
                        "a band meeting another across a re-entrant turn has no \
                         cap at the vertex to carry on from"
                    );
                };
                let hosts: Vec<Shape> =
                    ogeom_topo::ancestors_of(model, solid, &edge, ShapeType::Face)?;
                entry.2.solid = crate::junction::carried_into(
                    model,
                    &entry.2.solid,
                    &cap,
                    (&edge, &hosts),
                    outward,
                    junction.keeps(slot),
                    tol,
                )?;
                continue;
            }
            // A marched band ran out past the vertex as it was built, and a
            // band run on past it is built so.
            if !edge_curve(model, &edge, tol).is_ok_and(|(c, _)| matches!(c, Curve::Line(_)))
                || run_on
                    .iter()
                    .any(|&(k, at)| k == index && at.distance(junction.at) <= tol.confusion())
            {
                continue;
            }
            // The constant band alone, and what of it stands past the
            // vertex: its run-out, where it runs out there.
            let (built, wedges) = crate::support::collecting_wedges(|| {
                planar_fillet(model, solid, &edge, radius, None, tol)
            });
            built?;
            let Some(alone) = wedges.into_iter().next() else {
                continue;
            };
            // A straight band's sections stand square to its crease, so
            // its vertices reach as far past the vertex as any point of it.
            let mut past = f64::NEG_INFINITY;
            for vertex in explore_unique(model, &alone.solid, ShapeType::Vertex)? {
                if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                    let p = vertex.transform(model.datums())?.apply(data.point);
                    past = past.max((p - junction.at).dot(outward));
                }
            }
            if past <= gap {
                continue;
            }
            let run_out = crate::support::trimmed_wedge(
                model,
                &alone.solid,
                (junction.at, outward),
                ogeom_algo::shape_bounds(model, &alone.solid, tol)?.diagonal() * 2.0,
                tol,
            )?;
            let (solid_joined, _) = crate::junction::joined(
                model,
                [&entry.2.solid, &run_out],
                (junction.at, outward),
                junction.at,
                gap,
                tol,
            )?;
            entry.2.solid = solid_joined;
        }
    }
    let halves = with_corners(model, &mut set_aside, &junctions, solid, true, tol)?;
    let joined = join_flush(model, set_aside, &junctions, &halves, tol)?;
    let mut built = apply_set_aside(model, Built::from_nothing(solid.clone()), joined, tol)?;
    for edge in chain {
        built.history.delete(edge);
    }
    Ok(built)
}

/// A meeting of a chain under a law: its two edges, where they meet, and
/// the law's radius there.
type Meeting = ([usize; 2], Point, f64);

/// The straight bands of a chain under a law that run on past a convex
/// turn rather than meet a constant run-out there: where the two end
/// sections are no reflection of each other in the plane bisecting the turn
/// (a drafted wall beside a square one), the other band crosses this one's
/// run-out near the vertex, and a seam there between the law's band and the
/// constant run-out is what the cut leaves open. Each such band is built
/// again over its edge and past it as far as the constant band of the
/// meeting's radius runs out alone, the law held level past the edge's
/// end: one surface through the vertex, which past the vertex is that
/// constant band itself where it ends square to its crease. A constant band
/// trimmed by a face it runs into is left to the join. Returns each band's
/// index with the vertex it runs on past.
fn run_on_past_turns(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    (walk, owns, halved): (&[(usize, bool)], &[Profile], &[usize]),
    (junctions, meetings): (&[crate::junction::Junction], &[Meeting]),
    set_aside: &mut [(Vec<usize>, Shape, crate::support::Wedge)],
    tol: Tolerances,
) -> OgeomResult<Vec<(usize, Point)>> {
    let gap = crate::marched::band_fit_target(tol) * 2.0;
    // How far past each end of its oriented run, as fractions of its
    // length, and each turn it runs on past.
    let mut ends: Vec<[f64; 2]> = vec![[0.0; 2]; chain.len()];
    let mut alone_at: Vec<(usize, Point)> = Vec::new();
    for (junction, &(_, _, radius)) in junctions.iter().zip(meetings) {
        if junction.tangent || junction.reentrant || !junction.crossed {
            continue;
        }
        let (meet, bisector) = junction.flush;
        let reflect = |v: Vector| v - bisector * (2.0 * v.dot(bisector));
        let miss = radius.mul_add(
            (reflect(junction.sides[0]) - junction.sides[1]).magnitude(),
            2.0 * (junction.at - meet).dot(bisector).abs(),
        );
        if miss <= crate::marched::band_fit_target(tol) {
            continue;
        }
        for slot in 0..2 {
            let index = junction.pair[slot];
            let edge = &chain[index];
            let (curve, range) = edge_curve(model, edge, tol)?;
            if !matches!(curve, Curve::Line(_)) || halved.contains(&index) {
                continue;
            }
            let outward = junction.outward[slot];
            let (built, wedges) = crate::support::collecting_wedges(|| {
                planar_fillet(model, solid, edge, radius, None, tol)
            });
            built?;
            let Some(alone) = wedges.into_iter().next() else {
                continue;
            };
            let mut past = f64::NEG_INFINITY;
            for vertex in explore_unique(model, &alone.solid, ShapeType::Vertex)? {
                if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                    let p = vertex.transform(model.datums())?.apply(data.point);
                    past = past.max((p - junction.at).dot(outward));
                }
            }
            // The constant band runs out as the straight band of its
            // section, ended square to its crease where it has left the
            // solid; one trimmed by a face's plane it runs into is met as
            // before.
            let far = junction.at + outward * past;
            if past <= gap
                || !crate::junction::ends_in(model, &alone.solid, (far, outward), far, gap, tol)?
            {
                continue;
            }
            let first = ogeom_geom::Curve3d::point_at(&curve, range.0, tol)?;
            let last = ogeom_geom::Curve3d::point_at(&curve, range.1, tol)?;
            let start = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                last
            } else {
                first
            };
            let end = usize::from(start.distance(junction.at) > first.distance(last) / 2.0);
            ends[index][end] = ends[index][end].max(past / first.distance(last));
            alone_at.push((index, junction.at));
        }
    }
    let mut run_on = Vec::new();
    for (index, &[before, after]) in ends.iter().enumerate() {
        if before == 0.0 && after == 0.0 {
            continue;
        }
        let Some(slot) = walk.iter().position(|&(k, _)| k == index) else {
            continue;
        };
        let edge = &chain[index];
        let law = owns[slot].stretched(before, after);
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let part = if reversed {
            (-after, 1.0 + before)
        } else {
            (-before, 1.0 + after)
        };
        let (built, wedges) = crate::support::collecting_wedges(|| {
            straight_law_fillet(model, solid, edge, &law, reversed, part, tol)
        });
        built?;
        let [wedge] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|w| {
            ogeom_core::ogeom_err!(
                Construction,
                "a band run on past its turns built {} wedges, not one",
                w.len()
            )
        })?;
        let Some(entry) = set_aside
            .iter_mut()
            .find(|(held, _, _)| held.contains(&index))
        else {
            continue;
        };
        entry.2 = wedge;
        for (_, at) in alone_at.iter().filter(|(k, _)| *k == index) {
            run_on.push((index, *at));
        }
    }
    Ok(run_on)
}

/// The straight planar seat under a law of points: the envelope exact, a
/// cubic B-spline along the edge by the rational quadratic arc across it.
///
/// Along a straight edge the arc length is affine in the parameter, so the
/// law's cubic B-spline is the radius at every section; the apex line in
/// the same basis has the Greville points for control points. Each row of
/// the net is that apex control point moved out along the section's arm
/// by the law's control radius, which is the section's arc at every
/// parameter, since both are affine in the net. `part` is the stretch of
/// the edge built, by fractions of its stored parameter, `(0, 1)` for the
/// whole: the law is laid over that stretch, and the band is capped
/// square to the edge at its ends.
fn straight_law_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    law: &Profile,
    reversed: bool,
    part: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let seat = planar_seat(model, solid, edge, tol)?;
    straight_law_on(model, solid, edge, seat, law, reversed, part, tol)
}

/// [`straight_law_fillet`] on `seat`, the edge's own or one its band takes
/// a narrow strip beside it off whole on, running the way the edge does.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
#[allow(clippy::too_many_lines, reason = "one wedge, assembled end to end")]
fn straight_law_on(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    mut seat: crate::support::Seat,
    law: &Profile,
    reversed: bool,
    part: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (whole_start, whole_end) = (seat.start, seat.end);
    seat.start = whole_start + (whole_end - whole_start) * part.0;
    seat.end = whole_start + (whole_end - whole_start) * part.1;
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

    // The law over the edge's own parameter, start to end: the chain's
    // fraction where the edge runs with it, and its complement where
    // against.
    let (knots, radii) = {
        let (knots, values) = law.spline(tol);
        if reversed {
            (
                knots.iter().rev().map(|k| 1.0 - k).collect::<Vec<f64>>(),
                values.into_iter().rev().collect::<Vec<f64>>(),
            )
        } else {
            (knots, values)
        }
    };
    let greville: Vec<f64> = (0..radii.len())
        .map(|i| (knots[i + 1] + knots[i + 2] + knots[i + 3]) / 3.0)
        .collect();
    let apex_at = |u: f64| seat.start + (seat.end - seat.start) * u;
    let apex = [seat.start, seat.end];
    let end_radii = [radii[0], radii[radii.len() - 1]];
    let centre = [
        apex[0] + bisector * (end_radii[0] / depth),
        apex[1] + bisector * (end_radii[1] / depth),
    ];
    let contact_a = [centre[0] + n1 * end_radii[0], centre[1] + n1 * end_radii[1]];
    let contact_b = [centre[0] + n2 * end_radii[0], centre[1] + n2 * end_radii[1]];

    let half_cos = (sweep / 2.0).cos();
    let arc_mid = (n1 + n2) / (1.0 + n1.dot(n2));
    let weights = [1.0, half_cos, 1.0];
    let arms = [n1, arc_mid, n2];
    let lift = bisector / depth;

    let u_knots = ogeom_math::KnotVector::new(knots.clone(), 3)?;
    let v_knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
    let mut net: Vec<ogeom_math::Weighted<Point>> = Vec::with_capacity(radii.len() * 3);
    for (g, r) in greville.iter().zip(&radii) {
        for (arm, w) in arms.iter().zip(weights) {
            net.push(ogeom_math::Weighted::new(
                apex_at(*g) + (lift + *arm) * *r,
                w,
                tol,
            )?);
        }
    }
    let grid = ogeom_math::ControlGrid::new(net, radii.len(), 3)?;
    let surface: SurfaceGeometry =
        ogeom_geom::BSplineSurface::rational(u_knots, v_knots, grid)?.into();
    let surface_id = model.geometry_mut().add_surface(surface.clone());

    let va = [
        make_vertex(model, contact_a[0]).shape,
        make_vertex(model, contact_a[1]).shape,
    ];
    let vb = [
        make_vertex(model, contact_b[0]).shape,
        make_vertex(model, contact_b[1]).shape,
    ];

    // The tangency lines: the net's border rows, each in its own face's
    // plane, since every arm's end lies in its face.
    let rail = |model: &mut Model,
                arm: Vector,
                ends: [&Shape; 2],
                row: f64|
     -> OgeomResult<Shape> {
        let control: Vec<Point> = greville
            .iter()
            .zip(&radii)
            .map(|(g, r)| apex_at(*g) + (lift + arm) * *r)
            .collect();
        let knot_vector = ogeom_math::KnotVector::new(knots.clone(), 3)?;
        let curve = Curve::BSpline(ogeom_geom::BSplineCurve::new(knot_vector, control, tol)?);
        let built = make_edge_between(model, curve, (0.0, 1.0), ends[0], ends[1], tol)?.shape;
        let image = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(ogeom_math::Point2::new(0.0, row), ogeom_math::Direction2::X),
            0.0,
            1.0,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            image.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, 1.0),
        )?;
        Ok(built)
    };
    let rail_a = rail(model, n1, [&va[0], &va[1]], 0.0)?;
    let rail_b = rail(model, n2, [&vb[0], &vb[1]], 1.0)?;

    let arc_edge = |model: &mut Model, i: usize| -> OgeomResult<Shape> {
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
        let control = vec![
            ogeom_math::Weighted::new(contact_a[i], 1.0, tol)?,
            ogeom_math::Weighted::new(centre[i] + arc_mid * end_radii[i], half_cos, tol)?,
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

    // A planar face on `corners`' plane bounded by `edges` in ring order,
    // walked back where the corners turn against `outward`.
    let plane_face = |model: &mut Model,
                      corners: [Point; 3],
                      edges: Vec<Shape>,
                      outward: Vector|
     -> OgeomResult<Shape> {
        let plane = Plane::through(corners[0], Direction::new(outward, tol)?);
        let reach = corners
            .iter()
            .fold(1.0_f64, |m, p| m.max(p.distance(corners[0])))
            * 2.0
            + seat.start.distance(seat.end) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let turn = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
        let edges: Vec<Shape> = if turn.dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[edges], tol)?.shape)
    };
    let apex_v = [
        make_vertex(model, apex[0]).shape,
        make_vertex(model, apex[1]).shape,
    ];
    let leg = |model: &mut Model,
               rail: &Shape,
               ends: &[Shape; 2],
               contact: [Point; 2],
               outward: Vector|
     -> OgeomResult<Shape> {
        let edges = vec![
            segment_between(model, (&apex_v[0], apex[0]), (&ends[0], contact[0]), tol)?,
            rail.clone(),
            segment_between(model, (&ends[1], contact[1]), (&apex_v[1], apex[1]), tol)?,
            segment_between(model, (&apex_v[1], apex[1]), (&apex_v[0], apex[0]), tol)?,
        ];
        plane_face(model, [apex[0], contact[0], contact[1]], edges, outward)
    };
    let leg_a = leg(model, &rail_a, &va, contact_a, n1 * sign)?;
    let leg_b = leg(model, &rail_b, &vb, contact_b, n2 * sign)?;
    let cap = |model: &mut Model, i: usize, outward: Vector| -> OgeomResult<Shape> {
        let edges = vec![
            segment_between(model, (&apex_v[i], apex[i]), (&va[i], contact_a[i]), tol)?,
            if i == 0 { arc0.clone() } else { arc1.clone() },
            segment_between(model, (&vb[i], contact_b[i]), (&apex_v[i], apex[i]), tol)?,
        ];
        plane_face(model, [apex[i], contact_a[i], contact_b[i]], edges, outward)
    };
    let cap0 = cap(model, 0, -seat.along)?;
    let cap1 = cap(model, 1, seat.along)?;

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
/// or conical wall meets a perpendicular planar cap traces a torus.
///
/// Four seats, one parameterization. With `sigma` the wall's outward radial
/// sign and `tau` telling whether the wall extends away from the cap's
/// outward side, the tube's centre circle sits at radius `R − sigma tau r`
/// beside a cylinder (where a cone's wall leans, as
/// [`RevolvedSeat::ball`](crate::support::RevolvedSeat::ball) derives),
/// lifted `tau r` against the cap's normal, and `tau` alone decides whether
/// the wedge subtracts (the external rim and the hole's rim, both convex) or
/// fuses (the boss base and the blind hole's floor, both concave). The wedge
/// is always the same three revolved faces: a band of the wall, an annulus
/// of the cap, and the stretch of tube between them, a quarter beside a
/// cylinder.
fn revolved_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let seat = revolved_seat(model, solid, edge, rim, tol)?;
    revolved_fillet_on(model, solid, edge, seat, radius, tol)
}

/// [`revolved_fillet`] on `seat`, read from `edge` or from the faces either
/// side of a ring beside it.
fn revolved_fillet_on(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    seat: crate::support::RevolvedSeat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let ball = seat.ball(radius);
    let tube_rho = ball.rho;
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre + seat.up * ball.height;
    let flanks = revolved_flanks(model, &seat, -seat.tau * ball.wall_height, tube_rho, tol)?;

    // The blend: the tube between the tangencies, its natural normal away
    // from the tube's centre (into the wedge), so always reversed.
    let blend_band = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            seat.frame_at(tube_level, tol)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        if seat.slope == 0.0 {
            make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?
        } else {
            tube_band(model, &surface, &seat, &flanks, tol)?
        }
        .reversed()
    };

    let faces = [flanks.wall_band, flanks.annulus, blend_band];
    apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol)
}

/// The tube between a cone seat's tangency rings, over the whole turn.
///
/// In the tube's chart the angle about the tube is counted from its outer
/// equator, and beside a cone the ball's arc can cross it: where the wall
/// leans out over the ball, it touches below the equator. So the rings'
/// rows are given, not read back: the cap's straight along `tau up` from
/// the tube's centre and the wall's within a half turn of it, the arc
/// between them less than one.
fn tube_band(
    model: &mut Model,
    surface: &SurfaceGeometry,
    seat: &crate::support::RevolvedSeat,
    flanks: &crate::support::RevolvedFlanks,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use core::f64::consts::{PI, TAU};
    let cap_row = seat.tau.atan2(0.0);
    let lean = seat.sigma * seat.tau;
    let wall_row = {
        let row = (-lean * seat.slope).atan2(lean);
        row + ((cap_row - row) / TAU).round() * TAU
    };
    debug_assert!((cap_row - wall_row).abs() < PI);
    let row = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                ogeom_math::Point2::new(0.0, v),
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(1.0, 0.0), tol)?,
            ),
            0.0,
            TAU,
        )?
        .into())
    };
    ogeom_algo::make_band_between(
        model,
        surface,
        [
            (&flanks.wall_ring, row(wall_row)?),
            (&flanks.cap_ring, row(cap_row)?),
        ],
        tol,
    )
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
pub(crate) fn refind_edges(
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
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    let ball = seat.ball(radius);
    let tube_rho = ball.rho;
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre + seat.up * ball.height;
    let wall_level = seat.centre + seat.up * ball.wall_height;

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
    let dir_of = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    // At an end where the meridian cap would stop short of an oblique face
    // the arc ends on, the band runs on past the end, its window widened
    // until the ball's section stands clear of the material (a plane leaves
    // the section's triangle last at a corner: the rim, the wall's contact
    // and the cap's), and a quarter radius further; the cut trims it
    // against the face, as a straight band's run-out is.
    let mut trims: Vec<(Point, Vector)> = Vec::new();
    let (theta0, theta1) = {
        let mut window = [theta0, theta1];
        let step = (radius / seat.radius) / 8.0;
        for (slot, sign) in [(0_usize, -1.0), (1, 1.0)] {
            let end = window[slot];
            let at = seat.centre + dir_of(end) * seat.radius;
            let ahead = seat.up.cross(dir_of(end)) * sign;
            let Some(face) = rim_end_stops_short(model, solid, edge, &seat, mates, at, ahead, tol)?
            else {
                continue;
            };
            let probe = crate::support::material_probe(model, solid, tol)?;
            let (on, back) = face;
            let mut past = None;
            let mut k = 1.0_f64;
            // Over four radii along the rim, as a straight band's run-out
            // is read: a section that has not left the material by then, or
            // whose corner leaves it other than across the face the arc ends
            // on, meets a step or a ledge past the end rather than a ramp
            // across it, and the band stops at the end.
            while k * step * seat.radius <= radius * 4.0 {
                let theta = end + sign * k * step;
                // The corners stand on the hosts; each is read a hundredth
                // of the way in towards the others, inside the corner.
                let (apex, wall, cap) = (
                    seat.centre + dir_of(theta) * seat.radius,
                    wall_level + dir_of(theta) * ball.wall_rho,
                    seat.centre + dir_of(theta) * tube_rho,
                );
                let section = [
                    apex + ((wall - apex) + (cap - apex)) * 0.01,
                    wall + (apex - wall) * 0.01 + (cap - wall) * 0.01,
                    cap + (apex - cap) * 0.01 + (wall - cap) * 0.01,
                ];
                let mut clear = true;
                let mut ledge = false;
                for p in section {
                    if probe.holds(model, p, tol)? == ogeom_algo::Containment::In {
                        clear = false;
                    } else if (p - on).dot(back) > tol.confusion() {
                        ledge = true;
                    }
                }
                if ledge {
                    break;
                }
                if clear {
                    past = Some(k * step + 2.0 * step);
                    break;
                }
                k += 1.0;
            }
            if let Some(past) = past {
                window[slot] = end + sign * past;
                trims.push(face);
            }
        }
        if window[1] - window[0] >= tau_turn - tol.angular() {
            ogeom_bail!(
                Construction,
                "a rim arc's band run on past both its ends would close the turn"
            );
        }
        (window[0], window[1])
    };
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
    let wall_at = |theta: f64| wall_level + dir_at(theta) * ball.wall_rho;
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
    // The ball's arc in the meridian at `theta`, from the wall contact to
    // the cap contact: the tube's own circle, in the meridian plane, a
    // quarter of it beside a cylinder.
    let quarter = |model: &mut Model, theta: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let dir = dir_at(theta);
        let centre_m = tube_level + dir * tube_rho;
        let x_dir = Direction::new((dir - seat.up * seat.slope) * (seat.sigma * seat.tau), tol)?;
        let z_dir = Direction::new(dir.cross(seat.up) * seat.sigma, tol)?;
        let frame = Frame::new(centre_m, z_dir, x_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, ball.sweep), from, to, tol)?.shape)
    };

    let arc_apex = ring_arc(model, seat.centre, seat.radius, &va0, &va1)?;
    let arc_wall = ring_arc(model, wall_level, ball.wall_rho, &vw0, &vw1)?;
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
        let surface = seat.wall_surface(
            (-below, below),
            (ball.wall_height.min(0.0), ball.wall_height.max(0.0)),
            Direction::new(dir_at(theta0), tol)?,
            tol,
        )?;
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
    if !seat.additive() || trims.is_empty() {
        return apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol);
    }
    let mut wedge = crate::support::wedge_solid(model, &faces, false, tol)?;
    let span = (seat.radius + radius) * 4.0;
    for plane in trims {
        wedge = crate::support::trimmed_wedge(model, &wedge, plane, span, tol)?;
    }
    crate::support::apply_wedge_solid(model, solid, Some(edge), &wedge, true, tol)
}
