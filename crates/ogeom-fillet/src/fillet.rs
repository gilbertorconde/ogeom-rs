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

use crate::asym::Section;
use crate::law::{Profile, RadiusLaw};
use crate::support::{
    Seat, apply_wedge_solid, edge_curve, face_from_edges, planar_face, planar_seat, revolved_seat,
    segment_between, trimmed_wedge, wedge_solid,
};
use ogeom_algo::{Built, make_edge_between, make_vertex};
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
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
/// The blend engine takes the edge first, as a request of one, as it takes
/// many in [`fillet_edges`]; what it hands back is built as above.
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edge");
    // One edge is a request of one: the engine first, as for many.
    if let Some(built) = rounds_by_engine(model, solid, &[(edge.clone(), radius)], tol)? {
        return Ok(built);
    }
    refused_alone(model, solid, edge, radius, Check::Ball, tol)?;
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("fillet_edge");
    fillet_edge_wedge(model, solid, edge, radius, tol)
}

/// A request of two edges or more that no engine route builds, checked as
/// it is checked before any band of it is built: its corners of both
/// kinds and of mixed radii read, then each ball seated on both faces of
/// its edge and no two bands overlapping across a face they share
/// ([`crate::engine::request_fits`]), a concave band past its face's side kept for the
/// cut-off as there. One failing is refused here, by name, in the words
/// the current path refuses it in, for it refuses there before it builds
/// anything (ADR 0001, amendment "what M6 refuses by name"); one passing
/// goes on. A request of one edge is [`refused_alone`]'s: a strip beside
/// it may be taken off whole.
///
/// # Errors
///
/// The refusal, as the current path would raise it.
fn refused_as_checked(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    tol: Tolerances,
) -> OgeomResult<()> {
    if let [edge] = edges {
        return refused_alone(model, solid, edge, radii[0], Check::Request, tol);
    }
    // In the current path's order: its corners first, then the bands, in
    // node order where a corner of both kinds or a cusp is met.
    let request: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii.iter().copied()).collect();
    let signed = crate::setback::signed_corners(model, solid, &request, tol)?;
    let cusps = crate::setback::cusp_corners(model, solid, &request, tol)?;
    let (edges, radii): (Vec<Shape>, Vec<f64>) = if signed.is_empty() && cusps.is_empty() {
        (edges.to_vec(), radii.to_vec())
    } else {
        let mut sorted = request;
        sorted.sort_by_key(|(edge, _)| edge.node());
        sorted.into_iter().unzip()
    };
    mixed_corners(model, &edges, &radii)?;
    // As the current path checks them, within its cut-offs: a concave band
    // running past its face's side there is cut off at that side after it
    // is built ([`cut_off_bands`]), not refused.
    let (fits, _) = crate::support::with_cut_offs(|| {
        crate::engine::request_fits(model, solid, &edges, &radii, tol)
    });
    fits.map_err(|b| b.error)?;
    Ok(())
}

/// Which of the current path's first checks a request of one edge meets:
/// [`crate::engine::request_fits`] (a request through [`fillet_edges`]), or the ball on
/// its two faces alone ([`fillet_edge`]).
#[derive(Clone, Copy)]
enum Check {
    Request,
    Ball,
}

/// [`refused_as_checked`] for one edge: checked as the current path checks
/// it, and where that refuses, the strip or ring beside it read as taken
/// off whole as there ([`consuming`], on a copy of the model). Refused, in
/// the current path's words, only where that reads nothing either: the
/// current path then refuses before it builds anything.
///
/// # Errors
///
/// The refusal, as the current path would raise it.
fn refused_alone(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    check: Check,
    tol: Tolerances,
) -> OgeomResult<()> {
    if !radius.is_finite() || radius <= tol.confusion() {
        return Ok(());
    }
    let (checked, _) = crate::support::with_cut_offs(|| match check {
        Check::Request => {
            crate::engine::request_fits(model, solid, std::slice::from_ref(edge), &[radius], tol)
                .map(|_| ())
                .map_err(|b| b.error)
        }
        Check::Ball => {
            crate::support::ball_fits(model, solid, edge, radius, false, tol).map(|_| ())
        }
    });
    let refusal = match checked {
        Err(refusal @ OgeomError::Construction(_)) => refusal,
        _ => return Ok(()),
    };
    let mut scratch = model.clone();
    let (taken, _) = crate::support::with_cut_offs(|| {
        consuming(&mut scratch, solid, edge, radius, refusal, tol)
    });
    taken.map(|_| ())
}

#[cfg(feature = "old-path-log")]
pub(crate) use crate::old_path::entry as to_current_path;

/// [`fillet_edge`] on the current path alone: the band and its wedge, as
/// built before the blend engine took the edge.
///
/// # Errors
///
/// As [`fillet_edge`].
pub(crate) fn fillet_edge_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("round:fillet_edge_wedge");
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| {
            cut_off_bands(
                model,
                solid,
                |model| fillet_edge_own(model, solid, &parts[0], radius, tol),
                tol,
            )
        },
    )
}

/// Build with `build`, and cut each concave band on a curved edge that ran
/// past the side of the face it rolls on off at that side's plane: what
/// stands in front of the plane, within twice the band's radius of where
/// it ran past, goes, so the band's fill ends flush with the face's side.
/// Refused by name where `solid` itself holds anything there.
fn cut_off_bands(
    model: &mut Model,
    solid: &Shape,
    build: impl FnOnce(&mut Model) -> OgeomResult<Built>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (built, cuts) = crate::support::with_cut_offs(|| build(model));
    let mut built = built?;
    for cut in cuts {
        let reach = 2.0 * cut.radius;
        let (mut lo, mut hi) = (cut.points[0], cut.points[0]);
        for p in &cut.points {
            lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let (lo, hi) = (
            lo - Vector::new(reach, reach, reach),
            hi + Vector::new(reach, reach, reach),
        );
        let walls = [
            (cut.normal, cut.origin),
            (Vector::new(1.0, 0.0, 0.0), lo),
            (Vector::new(0.0, 1.0, 0.0), lo),
            (Vector::new(0.0, 0.0, 1.0), lo),
            (Vector::new(-1.0, 0.0, 0.0), hi),
            (Vector::new(0.0, -1.0, 0.0), hi),
            (Vector::new(0.0, 0.0, -1.0), hi),
        ];
        let tool = crate::engine::vertex_geometry::convex_block(model, &walls, tol)?;
        // The solid itself stands behind the plane there: the cut takes the
        // band's overhang and nothing it was given.
        let shared = ogeom_bool::common(model, solid, &tool, tol)?.shape;
        let held = if explore_unique(model, &shared, ShapeType::Face)?.is_empty() {
            0.0
        } else {
            ogeom_algo::volume_properties(model, &shared, ogeom_mesh::Deflection::default(), tol)?
                .mass
                .abs()
        };
        if held > tol.confusion() {
            ogeom_bail!(
                Construction,
                "a concave blend of radius {} runs past the side of the face it rolls on \
                 where the solid stands beyond that side, so its fill cannot be cut off \
                 there; see docs/CAPABILITIES.md, fillet.edge-blends",
                cut.radius
            );
        }
        let trimmed = ogeom_bool::cut(model, &built.shape, &tool, tol)?;
        built = Built::new(trimmed.shape, built.history.then(&trimmed.history));
    }
    Ok(built)
}

/// [`fillet_edge`] on a shape every node below which is its own.
fn fillet_edge_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("fillet_edge");
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
/// [`crate::engine::consumed_seat::consumed_seat`] finds one, and `refusal` otherwise.
fn consuming(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    refusal: OgeomError,
    tol: Tolerances,
) -> OgeomResult<Built> {
    crate::engine::old_marched_constructor("consuming");
    if !matches!(refusal, OgeomError::Construction(_)) {
        return Err(refusal);
    }
    let Some(seat) = crate::engine::consumed_seat::consumed_seat(model, solid, edge, radius, tol)?
    else {
        let rolled = match crate::rolled::rolled_fillet(model, solid, edge, radius, tol)? {
            Some(built) => Some(built),
            None => crate::consumed::consumed_slanted(model, solid, edge, radius, tol)?,
        };
        let Some(mut built) = rolled else {
            return match crate::consumed_rim::consumed_rim(model, solid, edge, radius, tol)? {
                Some(_) => Err(declined()),
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
    crate::engine::old_marched_constructor("or_rim_consumed");
    let mark = model.note_widened();
    match round(model) {
        Err(refusal @ OgeomError::Construction(_)) => {
            match crate::consumed_rim::consumed_rim(model, solid, edge, radius, tol)? {
                Some(_) => {
                    model.undo_widened(mark);
                    Err(declined())
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
    crate::engine::old_marched_constructor("fillet_edge_meeting");
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
        Route::Planar => planar_fillet(model, solid, edge, Section::Ball(radius), mates, tol),
        Route::Ruled(seat) => {
            crate::ruled::ruled_fillet(model, solid, edge, radius, mates, seat, tol)
        }
        // A cone's ball whose tube would reach its own axis is marched by
        // the engine.
        Route::Revolved(_, _, seat)
            if seat.slope != 0.0 && seat.ball(radius).rho - radius <= tol.confusion() =>
        {
            Err(declined())
        }
        Route::Revolved(c, crange, seat) => {
            let closed =
                ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
            // A whole rim's torus band is the engine's.
            if closed || crate::support::runs_whole_turn(model, &seat, &c, crange, tol)? {
                Err(declined())
            } else {
                crate::rim_arc::revolved_arc_fillet(
                    model, solid, edge, &c, radius, None, mates, None, tol,
                )
            }
        }
        // A seat with no closed form is marched by the engine.
        Route::Marched => Err(declined()),
    }
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
    Ruled(crate::engine::ruled_seat::RuledSeat),
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
        Curve::Line(_) => match crate::engine::ruled_seat::ruled_seat(model, solid, edge, tol)? {
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
/// A solid bounded by planes, each face one ring of straight edges, goes
/// first through the blend engine (docs/adr/0001-blend-architecture.md),
/// which builds the whole request in one pass, the same solid in any order.
/// Every face is set back to where its rounds' contact lines meet. Every
/// round is its exact cylinder between them, ended through its faces'
/// setback corners (skew on the cylinder where they stand at different
/// stations). A vertex closes as a sphere where one ball touches its three
/// faces there, flush in its third face where one round ends square to it,
/// and otherwise with one transfinite patch tangent to the three rounds'
/// ends within the bands' fit target. Every face is proved within the
/// material it rounds. A request the engine does not take (curved faces,
/// a face set back to a line, a vertex of four edges, an edge's two ends
/// set back past each other) is rounded as follows.
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
/// straight edges, all asked or one left sharp and convex, the request is
/// refused by name in any order, and so where a band's setback reaches its
/// edge's far end or the corner's cell reaches through the solid behind
/// the vertex. A request of two or more edges also
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edges");
    let asked: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), radius)).collect();
    if let Some(built) = rounds_by_engine(model, solid, &asked, tol)? {
        return Ok(built);
    }
    let radii = vec![radius; edges.len()];
    refused_as_checked(model, solid, edges, &radii, tol)?;
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("fillet_edges");
    ogeom_algo::on_own_nodes(model, solid, edges, |model, solid, parts| {
        cut_off_bands(
            model,
            solid,
            |model| fillet_edges_meeting(model, solid, parts, &radii, tol),
            tol,
        )
    })
}

/// How [`fillet_edges_corner`] closes a vertex where rounds meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletCorner {
    /// Today's behaviour: bands set back, one G1 patch between them.
    #[default]
    Setback,
    /// The ball rolls through the corner: no setback, the bands meet along
    /// the curves where they intersect, and any opening left is closed by
    /// the piece of the ball's sweep (sphere or torus) that rolls round it.
    Rolling,
}

/// [`fillet_edges`], each vertex where rounds meet closed as `corner` says.
///
/// [`FilletCorner::Setback`] is [`fillet_edges`] itself.
/// [`FilletCorner::Rolling`] closes every vertex where two or more of the
/// rounds meet by the ball rolling through it, every new face a plane, a
/// cylinder, a torus or a sphere, every contact G1:
///
/// - three convex rounds, or three concave: the ball's sphere where it
///   touches all three faces;
/// - two convex and one concave (a step's inside corner), or two concave
///   and one convex (a boss's foot at its upright edge): the lone round
///   runs to the vertex unbroken, to where its ball touches the face across
///   from it, and the other two run on round its end as the torus their
///   ball sweeps rolling round it, meeting it, that face and each other's
///   straight bands along circles;
/// - two rounds of one sense beside an edge left sharp: the bands meet
///   along the curve where they cross, as [`fillet_edges`] builds them.
///
/// The result does not depend on the order of `edges`.
///
/// # Errors
///
/// As [`fillet_edges`] for `Setback`. For `Rolling`,
/// [`OgeomError::Construction`], naming the vertex, where the ball's
/// sweep does not close it: a vertex
/// where rounds meet touching a curved face or joining other than three
/// edges, a convex and a concave round meeting beside an edge left sharp,
/// or the lone round's edge not square to the face across it (the sweep
/// round its end then no torus). A request the engine's rounds between
/// planes do not build is refused, never closed by setbacks instead.
pub fn fillet_edges_corner(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: FilletCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edges_corner");
    match corner {
        FilletCorner::Setback => fillet_edges(model, solid, edges, radius, tol),
        FilletCorner::Rolling => crate::engine::rolling_rounds(model, solid, edges, radius, tol),
    }
}

/// [`fillet_edges`] with a radius per edge.
///
/// With every radius the same, exactly [`fillet_edges`] at that radius. A
/// solid bounded by planes goes through the blend engine as there, a corner
/// of three radii closed with its transfinite patch; the rest is rounded as
/// follows.
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edges_with");
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
    if let Some(built) = rounds_by_engine(model, solid, edges, tol)? {
        return Ok(built);
    }
    refused_as_checked(model, solid, &shapes, &radii, tol)?;
    Err(declined())
}

/// The refusal of a round the blend engine does not build: by name, not
/// handed to another construction.
pub(crate) fn declined() -> OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "the blend engine builds no round of this kind here: the edges' \
         faces, their seats or the request lie outside what it rounds; see \
         docs/CAPABILITIES.md"
    )
}

/// The rounds of `edges` built by the blend engine in one pass, where the
/// solid is bounded by planes and every vertex is one the engine closes, or
/// where every edge asked is a rim between a planar cap and a coaxial wall,
/// or a straight edge along planes and cylinders or such a rim beside them;
/// `None` leaves the request to the current path. An edge asked twice at
/// one radius is asked once. Faces a route leaves in pieces on one carrier
/// are merged after ([`crate::engine::merged`]).
pub(crate) fn rounds_by_engine(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    match rounds_by_engine_pieces(model, solid, edges, tol)? {
        Some(built) => crate::engine::merged(model, solid, built, tol).map(Some),
        None => Ok(None),
    }
}

/// [`rounds_by_engine`] before the faces a route left in pieces on one
/// carrier are merged.
fn rounds_by_engine_pieces(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_topo::ShapeType;
    match model.kind_of(solid)? {
        ShapeType::Compound => return rounds_by_engine_compound(model, solid, edges, tol),
        ShapeType::Solid if model.children_of(solid)?.len() > 1 => {
            return rounds_by_engine_voided(model, solid, edges, tol);
        }
        _ => {}
    }
    // A solid with placed shapes in it (a prism's far end is its profile
    // moved, the same nodes standing twice) is read with its placements
    // baked in, the request traced onto the baked shapes and the history
    // carried back.
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    if baked.shape.is_same(solid) {
        return rounds_by_engine_unplaced(model, solid, edges, tol);
    }
    let mut traced = Vec::with_capacity(edges.len());
    for (edge, radius) in edges {
        let [only] = baked.history.trace(edge) else {
            return Ok(None);
        };
        traced.push((only.clone(), *radius));
    }
    let Some(built) = rounds_by_engine_unplaced(model, &baked.shape, &traced, tol)? else {
        return Ok(None);
    };
    Ok(Some(Built::new(
        built.shape,
        baked.history.then(&built.history),
    )))
}

/// [`rounds_by_engine`] on a compound of solids: each solid with the edges
/// it holds, the rest kept as they are, the results gathered into a
/// compound again. `None` where the compound holds anything but solids (a
/// sheet's rounds are the current path's), an edge lies in no solid or in
/// two, or any solid's rounds are handed back.
fn rounds_by_engine_compound(
    model: &mut Model,
    compound: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_topo::ShapeType;
    let parts = model.children_of(compound)?;
    if parts.is_empty()
        || parts
            .iter()
            .any(|p| !matches!(model.kind_of(p), Ok(ShapeType::Solid)))
    {
        return Ok(None);
    }
    let held: Vec<Vec<Shape>> = parts
        .iter()
        .map(|p| ogeom_topo::explore_unique(model, p, ShapeType::Edge))
        .collect::<OgeomResult<_>>()?;
    let mut asked: Vec<Vec<(Shape, f64)>> = vec![Vec::new(); parts.len()];
    for (edge, radius) in edges {
        let holders: Vec<usize> = (0..parts.len())
            .filter(|k| held[*k].iter().any(|e| e.is_same(edge)))
            .collect();
        let [k] = holders[..] else {
            return Ok(None);
        };
        asked[k].push((edge.clone(), *radius));
    }
    let mut out = Vec::with_capacity(parts.len());
    let mut history = ogeom_algo::History::new();
    for (part, own) in parts.iter().zip(&asked) {
        if own.is_empty() {
            out.push(part.clone());
            continue;
        }
        let Some(built) = rounds_by_engine(model, part, own, tol)? else {
            return Ok(None);
        };
        history = history.then(&built.history);
        out.push(built.shape);
    }
    let gathered = ogeom_algo::make_compound(model, &out)?.shape;
    history.modify(compound, gathered.clone());
    Ok(Some(Built::new(gathered, history)))
}

/// [`rounds_by_engine`] on a solid with voids: the outer shell rounded as a
/// solid alone where it holds edges of the request, and each void holding
/// some rounded as the solid it bounds turned inside out (a void's edges
/// are concave in the solid, and filling its corner is taking it off that
/// solid), then the shells put together again. `None` where an edge lies on
/// no shell or two, or a shell's rounds are handed back.
///
/// # Errors
///
/// [`OgeomError::Construction`],
/// refused by name, where a round reaches into a void: what a shell's
/// rounds leave of it touches another shell, or a void then stands outside
/// the outer shell (in what a round took off).
fn rounds_by_engine_voided(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_topo::ShapeType;
    let shells = model.children_of(solid)?;
    // The outer shell: the one whose bounds hold every other's.
    let mut sized = Vec::with_capacity(shells.len());
    for shell in &shells {
        sized.push(ogeom_algo::shape_bounds(model, shell, tol)?.diagonal());
    }
    let outer = (0..shells.len())
        .max_by(|a, b| sized[*a].total_cmp(&sized[*b]))
        .unwrap_or(0);
    let held: Vec<Vec<Shape>> = shells
        .iter()
        .map(|shell| ogeom_topo::explore_unique(model, shell, ShapeType::Edge))
        .collect::<OgeomResult<_>>()?;
    let mut asked: Vec<Vec<(Shape, f64)>> = vec![Vec::new(); shells.len()];
    for (edge, radius) in edges {
        let holders: Vec<usize> = (0..shells.len())
            .filter(|k| held[*k].iter().any(|e| e.is_same(edge)))
            .collect();
        let [k] = holders[..] else {
            return Ok(None);
        };
        asked[k].push((edge.clone(), *radius));
    }
    let mut history = ogeom_algo::History::new();
    let mut order: Vec<usize> = vec![outer];
    order.extend((0..shells.len()).filter(|k| *k != outer));
    let mut kept: Vec<Shape> = Vec::with_capacity(shells.len());
    for k in order {
        if asked[k].is_empty() {
            kept.push(shells[k].clone());
            continue;
        }
        // The shell as a solid of its own, the outer one as it stands and
        // a void turned inside out: whichever way bounds a valid solid.
        let mut alone = None;
        for shell in [shells[k].clone(), shells[k].reversed()] {
            let candidate = ogeom_algo::make_solid(model, &[shell])?.shape;
            if ogeom_algo::check(model, &candidate, tol)?.is_valid()
                && ogeom_algo::volume_properties(
                    model,
                    &candidate,
                    ogeom_mesh::Deflection::default(),
                    tol,
                )?
                .mass
                    > 0.0
            {
                alone = Some(candidate);
                break;
            }
        }
        let Some(alone) = alone else {
            return Ok(None);
        };
        let Some(built) = rounds_by_engine(model, &alone, &asked[k], tol)? else {
            return Ok(None);
        };
        let [shell] = &model.children_of(&built.shape)?[..] else {
            return Ok(None);
        };
        let shell = shell.clone();
        let turned = model
            .children_of(&alone)?
            .first()
            .is_some_and(|s| s.orientation() != shells[k].orientation());
        kept.push(if turned { shell.reversed() } else { shell });
        history = history.then(&built.history);
    }
    // No round reaches into a void: every void apart from every other shell
    // and wholly inside the outer one.
    let outside = ogeom_algo::make_solid(model, &[kept[0].clone()])?.shape;
    for i in 1..kept.len() {
        for j in 0..i {
            let apart = ogeom_algo::distance_between_shapes(
                model,
                &kept[i],
                &kept[j],
                ogeom_intersect::ExtremaOptions::default(),
                tol,
            )?;
            if apart.distance <= tol.confusion() {
                ogeom_core::ogeom_bail!(
                    Construction,
                    "a round reaches into a void: what it leaves of the solid meets the void"
                );
            }
        }
        let faces = ogeom_topo::explore_unique(model, &kept[i], ShapeType::Face)?;
        let Some(face) = faces.first() else {
            return Ok(None);
        };
        // A point of the void's: a triangle of its first face, at its
        // centroid (within the mesh's chord of the face).
        let mesh =
            ogeom_mesh::triangulate_face(model, face, ogeom_mesh::Deflection::default(), tol)?;
        let Some(&[i, j, k]) = mesh.triangles.first() else {
            return Ok(None);
        };
        let at = |n: u32| mesh.positions[n as usize].to_vector();
        let probe = ogeom_math::Point::from_vector((at(i) + at(j) + at(k)) / 3.0);
        if ogeom_algo::classify_in_solid_exact(model, &outside, probe, tol)?
            != ogeom_algo::Containment::In
        {
            ogeom_core::ogeom_bail!(
                Construction,
                "a round reaches into a void: the void stands in what the round takes off"
            );
        }
    }
    let whole = ogeom_algo::make_solid(model, &kept)?.shape;
    history.modify(solid, whole.clone());
    Ok(Some(Built::new(whole, history)))
}

/// [`rounds_by_engine`] on one solid of one shell with nothing placed in
/// it: the engine's routes, and where none takes the request, the creases
/// carried out of the solid, then the marched seats' strips (ADR 0001,
/// amendments "M5e" and "marched seats").
fn rounds_by_engine_unplaced(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if let Some(built) = rounds_by_engine_routes(model, solid, edges, tol)? {
        return Ok(Some(built));
    }
    let mut once: Vec<(Shape, f64)> = Vec::with_capacity(edges.len());
    for (edge, radius) in edges {
        if !once.iter().any(|(e, _)| e.is_same(edge)) {
            once.push((edge.clone(), *radius));
        }
    }
    if once.is_empty() {
        return Ok(None);
    }
    // Rounds either side of a cusp, each its own tool (amendment "rounds
    // either side of a cusp").
    if let Some(built) = crate::engine::cusp_side_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // A rim arc ending obliquely on flats, rolled on past its ends (amendment
    // "a rim arc ending obliquely on a flat").
    if let Some(built) = crate::engine::oblique_rim_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // Then straight creases no route closes, each carried on past both ends
    // out of the solid (amendment "M5e"'s free ends).
    if once
        .iter()
        .all(|(_, r)| r.is_finite() && *r > tol.confusion())
        && let Some(built) = crate::engine::carried_out_rounds(model, solid, &once, false, tol)?
    {
        return Ok(Some(built));
    }
    // A rim whose band runs on past an end of its arc, round its circle,
    // until the end face trims it, or else ends in the plane through the
    // axis at its vertex; others asked with it rounded alone (amendment
    // "corners where a round's rail misses").
    if let Some(built) =
        crate::engine::rim_runs_on_rounds(model, solid, &once, rounds_by_engine, tol)?
    {
        return Ok(Some(built));
    }
    // Creases meeting only across convex turns, each built as it is alone,
    // the cuts taking their union (amendment "M5e").
    if let Some(built) = crate::engine::alone_rounds(model, solid, &once, rounds_by_engine, tol)? {
        return Ok(Some(built));
    }
    // Three convex creases at a vertex one ball rounds, each built as it is
    // alone, the corner's cell beyond the ball taken off too (amendment "a
    // corner one ball rounds, from its bands alone").
    if let Some(built) =
        crate::engine::ball_corner_rounds(model, solid, &once, rounds_by_engine, tol)?
    {
        return Ok(Some(built));
    }
    if let Some(built) = crate::engine::marched_rounds(model, solid, &once, false, tol)? {
        return Ok(Some(built));
    }
    // A cone rim asked alone whose exact torus band nothing above closes,
    // its band fitted through the ball's arcs as the march builds it, as
    // the current path retries it: the boolean resolves some meetings of
    // the fitted band it does not of the exact one. With other edges the
    // rim's band meets theirs at their corners, which the march's
    // junctions do not read as the exact rounds do.
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    if edges.len() == 1
        && has_cone_rim(model, solid, &edges, tol)?
        && let Some(built) = crate::engine::marched_rounds(model, solid, &once, true, tol)?
    {
        return Ok(Some(built));
    }
    // Arcs of one closed crease asked together: one round over the loop,
    // the march's band from the first arc.
    if let Some(built) = crate::engine::crease_loop_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // A contour on a planar host beside swept or circular walls that no
    // route above closes, its bands joined at their junctions (amendment
    // "contour chains on swept walls").
    if let Some(built) = crate::engine::contour_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // Last, groups of edges whose rounds stand apart, each built as it is
    // alone (amendment "rounds standing apart").
    crate::engine::apart_rounds(model, solid, &once, rounds_by_engine, tol)
}

/// The engine's routes for [`rounds_by_engine_unplaced`], in turn, but the
/// carried-out creases and the marched seats, tried last.
fn rounds_by_engine_routes(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut once: Vec<(Shape, f64)> = Vec::with_capacity(edges.len());
    for (edge, radius) in edges {
        if !once.iter().any(|(e, _)| e.is_same(edge)) {
            once.push((edge.clone(), *radius));
        }
    }
    // Nothing asked, or a radius that rounds nothing, is refused by the
    // current path, in its words.
    if once.is_empty()
        || once
            .iter()
            .any(|(_, r)| !(r.is_finite() && *r > tol.confusion()))
    {
        return Ok(None);
    }
    // A ball resting on a far edge changes its seat. Read that before a
    // cylindrical strip can extend across the narrow face (ADR 0001, M5d).
    if let Some(built) = crate::engine::far_edge_rounds(model, solid, &once, rounds_by_engine, tol)?
    {
        return Ok(Some(built));
    }
    let request = crate::engine::Request::from_rounds(&once, tol)?;
    if let Some(built) = crate::engine::run(model, solid, &request, tol)? {
        return Ok(Some(built));
    }
    // A solid bounded by planes holding seams (two flats a boolean joined
    // into one): read with its seams dissolved, as the chamfers read it, on
    // a copy of the model, the history told of the solid as given.
    {
        let mut scratch = model.clone();
        if let crate::engine::Dissolved::Into(seamless) =
            crate::engine::dissolve(&mut scratch, solid, tol)?
            && let Some(rounds) = seamless.rounds(&once)
            && !rounds.is_empty()
            && let Some(built) = crate::engine::run(
                &mut scratch,
                &seamless.solid,
                &crate::engine::Request::from_rounds(&rounds, tol)?,
                tol,
            )?
        {
            *model = scratch;
            return Ok(Some(Built::new(
                built.shape,
                seamless.carry(&built.history),
            )));
        }
    }
    // A round wider than a narrow face beside its edge: the face taken off
    // whole, the crease its neighbours close in rounded (amendment "M5d").
    if let Some(built) =
        crate::engine::taken_whole_rounds(model, solid, &once, rounds_by_engine, tol)?
    {
        return Ok(Some(built));
    }
    // Rounds meeting at inside corners, each carried on past the other
    // (amendment "M5e").
    if let Some(built) = crate::engine::meeting_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // The same where a seam splits the face the creases share at their
    // corner: read with the seams dissolved, on a copy of the model.
    {
        let mut scratch = model.clone();
        if let Some(seamless) = crate::engine::dissolve_curved(&mut scratch, solid, tol)?
            && let Some(rounds) = seamless.rounds(&once)
            && !rounds.is_empty()
            && let Some(built) =
                crate::engine::meeting_rounds(&mut scratch, &seamless.solid, &rounds, tol)?
        {
            *model = scratch;
            return Ok(Some(Built::new(
                built.shape,
                seamless.carry(&built.history),
            )));
        }
    }
    // A straight crease ending at a band an earlier round left, seated on
    // planes: its band runs on into that band, which trims it (amendment
    // "a later round reaching an earlier band at its end"), before the
    // curved routes would stop it there.
    if crate::engine::reaches_earlier_band(model, solid, &once, false, tol)?
        && let Some(built) = crate::engine::carried_out_rounds(model, solid, &once, false, tol)?
    {
        return Ok(Some(built));
    }
    // The same where the earlier band is a rim's torus, read as a band
    // only here: where the crease's ball does not seat along its faces up
    // to the torus, the torus stays the face it ends flush on, as a free
    // end (the straight creases carried out, below).
    if crate::engine::reaches_earlier_band(model, solid, &once, true, tol)?
        && let Some(built) = crate::engine::carried_out_rounds(model, solid, &once, true, tol)?
    {
        return Ok(Some(built));
    }
    // A drum's ruling ending at an earlier band runs on into it, trimmed by
    // it (amendment "a later round reaching a band seated on a curved
    // face"), before the curved routes would stop it at that band.
    if let Some(built) = crate::engine::ruling_out_rounds(model, solid, &once, tol)? {
        return Ok(Some(built));
    }
    // A solid with curved faces, read with its seams dissolved where it
    // holds any: rims between a planar cap and a coaxial cylinder or cone,
    // each the torus its ball sweeps; else straight edges along planes and
    // cylinders, and rims beside them, through the trimmed assembly.
    if crate::engine::bounded_by_planes(model, solid, tol)? {
        return Ok(None);
    }
    let on_curved = |model: &mut Model, solid: &Shape, once: &[(Shape, f64)]| {
        // Every edge asked a seam of the dissolved solid: nothing to round.
        if once.is_empty() {
            return Ok(None);
        }
        let once = whole_rims(model, solid, once, tol)?;
        if let Some(built) = crate::engine::round_rims(model, solid, &once, tol)? {
            return Ok(Some(built));
        }
        crate::engine::round_curved(model, solid, &once, tol)
    };
    // A seam (an edge between faces of one carrier a boolean split) names
    // nothing, as for chamfers.
    let mut kept = Vec::with_capacity(once.len());
    for (edge, radius) in &once {
        // A face's own seam (an edge one face walks both ways) too.
        let own =
            ogeom_topo::ancestors_of(model, solid, edge, ogeom_topo::ShapeType::Face)?.len() == 1;
        if !own && !crate::engine::is_seam(model, solid, edge, tol)? {
            kept.push((edge.clone(), *radius));
        }
    }
    if kept.is_empty() {
        return Ok(None);
    }
    let once = kept;
    // The seams dissolved on a copy of the model: dissolving builds over
    // the solid's own nodes, and a request the engine hands back must
    // reach the current path as it came.
    let mut scratch = model.clone();
    match crate::engine::dissolve_curved(&mut scratch, solid, tol)? {
        Some(seamless) => {
            let Some(rounds) = seamless.rounds(&once) else {
                return Ok(None);
            };
            if let Some(built) = on_curved(&mut scratch, &seamless.solid, &rounds)? {
                *model = scratch;
                return Ok(Some(Built::new(
                    built.shape,
                    seamless.carry(&built.history),
                )));
            }
            // A round reaching a face it does not end at: its tools, read
            // on the dissolved solid and named for the edges asked, applied
            // by booleans on the solid as given.
            let Some(tools) =
                crate::engine::round_non_local_tools(&mut scratch, &seamless.solid, &rounds, tol)?
            else {
                // The general route: the far feature taken off a copy of
                // the solid as given, the rounds built there.
                let Some((far, near)) =
                    crate::engine::round_far_faces(&mut scratch, &seamless.solid, &rounds, tol)?
                else {
                    return Ok(None);
                };
                let inputs = |faces: &[Shape]| -> Vec<Shape> {
                    faces
                        .iter()
                        .flat_map(|f| seamless.face_inputs_of(f))
                        .collect()
                };
                let Some(built) = rounds_defeatured(
                    &mut scratch,
                    solid,
                    &once,
                    &inputs(&far),
                    &inputs(&near),
                    tol,
                )?
                else {
                    return Ok(None);
                };
                *model = scratch;
                return Ok(Some(built));
            };
            let tools: Vec<crate::engine::Tool> = tools
                .into_iter()
                .map(|t| crate::engine::Tool {
                    faces: t
                        .faces
                        .iter()
                        .flat_map(|(face, from)| {
                            seamless
                                .inputs_of(from)
                                .into_iter()
                                .map(|e| (face.clone(), e))
                        })
                        .collect(),
                    ..t
                })
                .collect();
            let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
            let Some(built) = crate::engine::apply_tools(&mut scratch, solid, &edges, &tools, tol)?
            else {
                return Ok(None);
            };
            *model = scratch;
            Ok(Some(built))
        }
        None => {
            if let Some(built) = on_curved(model, solid, &once)? {
                return Ok(Some(built));
            }
            let mut scratch = model.clone();
            let Some(tools) =
                crate::engine::round_non_local_tools(&mut scratch, solid, &once, tol)?
            else {
                let Some((far, near)) =
                    crate::engine::round_far_faces(&mut scratch, solid, &once, tol)?
                else {
                    return Ok(None);
                };
                let Some(built) = rounds_defeatured(&mut scratch, solid, &once, &far, &near, tol)?
                else {
                    return Ok(None);
                };
                *model = scratch;
                return Ok(Some(built));
            };
            let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
            let Some(built) = crate::engine::apply_tools(&mut scratch, solid, &edges, &tools, tol)?
            else {
                return Ok(None);
            };
            *model = scratch;
            Ok(Some(built))
        }
    }
}

/// `once` on `solid` where the rounds reach the faces `far` and end at the
/// faces `near`, by the general non-local route (ADR 0001, M5c): the far
/// feature taken off a copy, the rounds built there by the engine, the
/// differences applied by booleans on `solid`. `None` where any step
/// declines.
fn rounds_defeatured(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    far: &[Shape],
    near: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(tools) =
        crate::engine::defeatured_tools(model, solid, once, far, near, rounds_by_engine, tol)?
    else {
        return Ok(None);
    };
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    crate::engine::apply_tools(model, solid, &edges, &tools, tol)
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
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("round:fillet_edges_meeting");
    // A corner where a convex and a concave band meet, a cusp, or three
    // bands of different radii is the engine's or refused by name.
    let request: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii.iter().copied()).collect();
    if !crate::setback::signed_corners(model, solid, &request, tol)?.is_empty()
        || !crate::setback::cusp_corners(model, solid, &request, tol)?.is_empty()
        || !mixed_corners(model, edges, radii)?.is_empty()
    {
        return Err(declined());
    }
    let round = |model: &mut Model| {
        with_fitted_rims(
            model,
            solid,
            edges,
            |model, rims| fillet_edges_meeting_with(model, solid, edges, radii, rims, tol),
            tol,
        )
    };
    // A request of one edge is that edge alone, a ring beside a rim taken
    // off whole as there.
    match (edges, radii) {
        ([edge], [radius]) => or_rim_consumed(model, solid, edge, *radius, round, tol),
        _ => round(model),
    }
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
        return fillet_edges_own(model, solid, edges, radii, None, rims, tol);
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
        match fillet_edges_own(model, solid, edges, radii, Some(&mut taken), rims, tol) {
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
    fillet_edges_own(model, solid, edges, radii, None, rims, tol)
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
    tol: Tolerances,
) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
    let rolled = match crate::engine::request_fits(model, solid, edges, radii, tol) {
        Ok(rolled) => rolled,
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
    // No ball rolls round a corner whose bands differ in radius: the
    // corner tool's ball would touch a band of another radius nowhere
    // along its rim. Such a corner is left to the bands here and closed by
    // the setback patch once they stand ([`mixed_corners`]).
    corners.retain(|(_, at)| at.iter().all(|&r| r == at[0]));
    // A corner the tool rounds, or the solid unchanged where it does not
    // speak the corner, which keeps the bands' caps; a corner of planes is
    // the engine's, refused here by name.
    let round =
        |model: &mut Model, built: Built, vertex: &Shape, radius: f64| -> OgeomResult<Built> {
            match crate::corner::round_vertex_with(model, &built.shape, vertex, radius, tol) {
                Ok(Some(rounded)) => Ok(Built {
                    shape: rounded.shape,
                    history: built.history.then(&rounded.history),
                }),
                Ok(None) => Err(declined()),
                Err(OgeomError::Construction(_)) => Ok(built),
                Err(other) => Err(other),
            }
        };
    // The corner first, the bands after it: bands built first crash into
    // each other at an apex, and the corner tool wants the sharp vertex's
    // planes, not the bands' mitre. Corners are cut one after another.
    let mut built: Option<Built> = None;
    let mut rounded: Vec<Point> = Vec::new();
    for (vertex, at) in &corners {
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
    let mut mates = crate::engine::request_mates(model, edges, &rounded, tol)?;
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
            crate::engine::junction::junctions(model, solid, edges, &given, &mates, radii, tol)
                .unwrap_or_default();
        if meeting.as_ref().is_some_and(|taken| !taken.crossed) {
            found.retain(|junction| !junction.crossed);
        }
        crate::engine::joinable(model, edges, found)?
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
                let halves = crate::engine::with_corners(
                    model,
                    &mut set_aside,
                    &junctions,
                    &b.shape,
                    false,
                    tol,
                )?;
                let joined = crate::engine::join_flush(
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
                    crate::engine::fit_to_junctions(
                        model,
                        &mut wedge,
                        (index, edge, &hosts),
                        &junctions,
                        radii[index],
                        tol,
                    )?;
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
    let halves =
        crate::engine::with_corners(model, &mut set_aside, &junctions, &built.shape, false, tol)?;
    let joined = crate::engine::join_flush(model, set_aside, &junctions, &halves, tol)?;
    apply_set_aside(model, built, joined, tol)
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
pub(crate) fn apply_set_aside(
    model: &mut Model,
    mut built: Built,
    set_aside: Vec<(Shape, crate::support::Wedge)>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = to_current_path("marched-backedge:apply_set_aside");
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
        // A concave ruled wedge adds planar end patches beside the cap it
        // extends. Dissolve only these same-domain seams before assigning
        // provenance: the cap remains the cap, and its outer boundary is
        // the new circular rim, not several ambiguous cap/patch edges.
        let mut planar_images = Vec::new();
        for &i in &group {
            for face in explore_unique(model, &set_aside[i].1.solid, ShapeType::Face)? {
                for image in step.history.trace(&face) {
                    let planar = model
                        .node(image)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| model.geometry().surface(d.surface))
                        .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
                    if planar {
                        planar_images.push(image.clone());
                    }
                }
            }
        }
        let (unified, _) =
            ogeom_heal::unify_same_domain_around(model, &step.shape, &planar_images, tol)?;
        let step = Built {
            shape: unified.shape,
            history: step.history.then(&unified.history),
        };
        // The chain's history already says each edge is gone; the blends
        // it generated are recorded on the composed history, where an edge
        // consumed before this step still answers for them.
        let mut history = built.history.then(&step.history);
        let result = ogeom_topo::explore_unique(model, &step.shape, ShapeType::Face)?;
        // A wedge's legs and caps can melt into the solid's existing cap.
        // That cap is an image of its original face, not a new blend: giving
        // it the blend's name too makes their common circular rim disappear
        // from a consumer's disjoint-face edge names (a following Face fillet).
        let inherited: Vec<_> = explore_unique(model, &built.shape, ShapeType::Face)?
            .iter()
            .flat_map(|face| step.history.trace(face).to_vec())
            .collect();
        for &i in &group {
            let (edge, wedge) = &set_aside[i];
            for face in ogeom_topo::explore_unique(model, &wedge.solid, ShapeType::Face)? {
                let owner = wedge
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or(edge, |(_, e)| e);
                for image in step.history.trace(&face) {
                    if result.iter().any(|f| f.is_same(image))
                        && !inherited.iter().any(|f| f.is_same(image))
                    {
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edge_variable");
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| {
            fillet_edge_variable_own(model, solid, &parts[0], start_radius, end_radius, tol)
        },
    )
}

/// A linear law on a straight planar seat, including its level limit.
fn fillet_edge_variable_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    start_radius: f64,
    end_radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let law = RadiusLaw::Linear {
        start: start_radius,
        end: end_radius,
    };
    let forward = edge.oriented(ogeom_topo::Orientation::Forward);
    if let Some(request) = crate::engine::Request::from_varying(model, &forward, &law)?
        && let Some(built) = crate::engine::run(model, solid, &request, tol)?
    {
        return Ok(built);
    }
    let profile = Profile::new(&law, false, tol)?;
    crate::engine::linear_round(model, solid, edge, &profile, tol)
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edges_variable");
    if chain.is_empty() {
        ogeom_bail!(
            Construction,
            "a variable radius over no edges rounds nothing"
        );
    }
    // One straight edge between planes: the blend engine first.
    if let [edge] = chain
        && let Some(request) = crate::engine::Request::from_varying(model, edge, law)?
        && let Some(built) = crate::engine::run(model, solid, &request, tol)?
    {
        return Ok(built);
    }
    // A chain of straight edges between planes: each edge rounded under its
    // own stretch of the law, the bands meeting at its turns as the engine
    // meets rounds (amendment "M5e").
    if chain.len() > 1
        && let Some(specs) = straight_chain_laws(model, chain, law, tol)?
        && let Some(built) = crate::engine::meeting_rounds_with(model, solid, &specs, tol)?
    {
        return Ok(built);
    }
    // Any other chain: the bands the march builds under the law, applied
    // by the engine as one-boolean tools.
    if let Some(built) = crate::engine::varying_rounds(model, solid, chain, law, tol)? {
        return Ok(built);
    }
    ogeom_bail!(Construction, "the variable-law request leaves no blend")
}

/// Each edge of `chain` (straight, end to end, open or closed) with the
/// stretch of `law` it covers by arc length, and the point the stretch
/// starts at; `None` where an edge is not straight or the edges do not run
/// end to end.
///
/// # Errors
///
/// As reading the edges and laying the law.
fn straight_chain_laws(
    model: &Model,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<(Shape, crate::engine::MeetingSpec)>>> {
    let point = |v: &Shape| {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
    };
    // Each edge's ends, the way it runs (its orientation's way).
    let mut ends: Vec<[Point; 2]> = Vec::with_capacity(chain.len());
    for e in chain {
        let (curve, _) = edge_curve(model, e, tol)?;
        if !matches!(curve, Curve::Line(_)) {
            return Ok(None);
        }
        let Some((a, b)) = ogeom_algo::edge_vertices(model, e)? else {
            return Ok(None);
        };
        let (Some(a), Some(b)) = (point(&a), point(&b)) else {
            return Ok(None);
        };
        ends.push([a, b]);
    }
    // `edge_vertices` reads them the way each edge runs already.
    let runs = ends;
    let near = |p: Point, q: Point| p.distance(q) <= tol.confusion() * 1e3;
    // The walk from the first edge's start, each step an edge and the
    // point it starts at.
    let mut walk: Vec<(usize, Point, Point)> = vec![(0, runs[0][0], runs[0][1])];
    let mut used = vec![false; chain.len()];
    used[0] = true;
    let mut head = runs[0][1];
    while let Some(k) =
        (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, head)))
    {
        let with = near(runs[k][0], head);
        used[k] = true;
        let (from, to) = if with {
            (runs[k][0], runs[k][1])
        } else {
            (runs[k][1], runs[k][0])
        };
        walk.push((k, from, to));
        head = to;
    }
    let closed = near(head, runs[0][0]);
    // An open chain asked from its middle: on back from the first edge's
    // start, the chain starting where that walk ends.
    if !closed {
        let mut tail = runs[0][0];
        while let Some(k) =
            (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, tail)))
        {
            let with = near(runs[k][1], tail);
            used[k] = true;
            let (from, to) = if with {
                (runs[k][0], runs[k][1])
            } else {
                (runs[k][1], runs[k][0])
            };
            walk.insert(0, (k, from, to));
            tail = from;
        }
    }
    if used.iter().any(|u| !u) {
        return Ok(None);
    }
    let profile = crate::law::Profile::new(law, closed, tol)?;
    let total: f64 = walk.iter().map(|(_, a, b)| a.distance(*b)).sum();
    if total <= tol.confusion() {
        return Ok(None);
    }
    let mut out: Vec<(Shape, crate::engine::MeetingSpec)> = Vec::with_capacity(chain.len());
    let mut before = 0.0;
    for (k, from, to) in walk {
        let length = from.distance(to);
        let (f0, f1) = (before / total, (before + length) / total);
        before += length;
        out.push((
            chain[k].clone(),
            crate::engine::MeetingSpec::Law(profile.window(f0, f1), from),
        ));
    }
    Ok(Some(out))
}

/// Whether `chain` runs end to end and back to its start, each edge
/// between the same two faces of `solid` as the first, none twice.
pub(crate) fn closed_chain(model: &Model, solid: &Shape, chain: &[Shape]) -> OgeomResult<bool> {
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
pub(crate) fn one_curve_span(
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
pub(crate) fn law_loop_fits(
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
        // As one edge alone is read: the constant round at the law's
        // largest radius holds the ball at some station.
        crate::support::ball_fits(model, solid, edge, law.largest(), false, tol)?;
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
pub(crate) fn law_chain_fits(
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
pub(crate) fn law_ball_on_faces(
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
        // The constant round at the law's largest radius, as a chain's
        // request is read before its stations: a contact that runs off a
        // strip past the solid at one radius runs off at every larger one,
        // and a ball held at no station there is one the constant round
        // takes the strip off for, or refuses, not a band cut off.
        crate::support::ball_fits(model, solid, edge, law.largest(), false, tol)?;
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
pub(crate) fn chain_mates(
    model: &Model,
    chain: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<Mate>> {
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

/// The straight-edge fillet: the rolling ball's envelope is a cylinder.
pub(crate) fn planar_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    section: Section,
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let radius = section.size();
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
        let Some(inward) = section.centre(&seat, tol)? else {
            continue;
        };
        // Where each point of the section crosses the plane, carried along
        // `outward`: a plane passes a section last at one of its corners.
        let rate = -outward.dot(normal);
        let crossings = section
            .points(&seat, at, inward, tol)?
            .map(|p| (origin - p).dot(normal) / rate);
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
        let Some(inward) = section.centre(&seat, tol)? else {
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
            let points = section.points(&seat, at, inward, tol)?;
            let boss = mates.is_none_or(|(index, mates)| Mate::alone_at(mates, index, at, tol));
            reach = crate::support::run_out_reach(
                model,
                solid,
                edge,
                hosts,
                (at, outward),
                &points,
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
    match section {
        Section::Ball(radius) => {
            seated_fillet(model, solid, &seat, radius, Some(edge), &trims, tol)
        }
        Section::Conic { .. } => Err(declined()),
    }
}

/// From a point of a straight seat's crease to the ball's centre in the
/// section square to the crease there; `None` where the faces meet too
/// sharply or too flat to seat a ball.
pub(crate) fn seat_centre(
    seat: &Seat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
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
pub(crate) fn section_at(seat: &Seat, at: Point, inward: Vector, radius: f64) -> [Point; 6] {
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

/// Where the rolling ball of `radius` seats between a straight seat's two
/// planes: the faces in the order that turns positively about the crease,
/// so the ball's arc and the cylinder's chart run from the first face's
/// contact to the second's; their normals; the sign that mirrors every
/// direction onto a concave seat's open side; the ball's centre from the
/// crease, square to it; and the arc's sweep.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BallSeat {
    /// The seat's face the arc starts on.
    pub first: usize,
    /// The face it ends on.
    pub second: usize,
    /// The first face's outward normal.
    pub n1: Vector,
    /// The second's.
    pub n2: Vector,
    /// `1` on a convex seat, `-1` on a concave one.
    pub sign: f64,
    /// From a point of the crease to the ball's centre.
    pub offset: Vector,
    /// The angle the arc turns through, about the crease.
    pub sweep: f64,
}

/// [`BallSeat`] of `seat` at `radius`; refused where the faces meet too
/// sharply or are parallel.
pub(crate) fn ball_seat(seat: &Seat, radius: f64, tol: Tolerances) -> OgeomResult<BallSeat> {
    // Order the two faces so the blend arc sweeps positively about the edge
    // direction: the cylinder's parameterization and every arc then run
    // from the first face's tangency line to the second's.
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
    let offset = bisector * (radius / depth);
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2));
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }
    Ok(BallSeat {
        first,
        second,
        n1,
        n2,
        sign,
        offset,
        sweep,
    })
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
    let BallSeat {
        n1,
        n2,
        sign,
        offset,
        sweep,
        ..
    } = ball_seat(seat, radius, tol)?;
    let centre = seat.start + offset;

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

/// `asked` with each rim closed: where the pieces of a circle asked all
/// take one radius and the solid's arcs of that circle run the whole turn,
/// the arcs not asked join at that radius. One value rounds a tangent chain
/// (ADR 0001, "what M6 refuses by name"), and the current path rounds the
/// whole turn where a rim closes; a rim that stops, or pieces asked at
/// different radii, are as asked.
///
/// # Errors
///
/// As reading the solid's edges.
fn whole_rims(
    model: &Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, f64)>> {
    let mut out = asked.to_vec();
    let circles: Vec<Option<ogeom_math::Circle>> = asked
        .iter()
        .map(|(e, _)| rim_circle(model, e, tol))
        .collect();
    if circles.iter().all(Option::is_none) {
        return Ok(out);
    }
    let all = ogeom_topo::explore_unique(model, solid, ShapeType::Edge)?;
    let on: Vec<(Shape, Option<ogeom_math::Circle>)> = all
        .into_iter()
        .map(|e| {
            let c = rim_circle(model, &e, tol);
            (e, c)
        })
        .collect();
    let mut done: Vec<ogeom_math::Circle> = Vec::new();
    for (k, circle) in circles.iter().enumerate() {
        let Some(circle) = circle else {
            continue;
        };
        if done.iter().any(|d| same_circle(d, circle, tol)) {
            continue;
        }
        done.push(*circle);
        let radius = asked[k].1;
        let alike = asked.iter().zip(&circles).all(|((_, r), c)| {
            c.as_ref().is_none_or(|c| !same_circle(c, circle, tol)) || *r == radius
        });
        if !alike {
            continue;
        }
        let pieces: Vec<&Shape> = on
            .iter()
            .filter(|(_, c)| c.as_ref().is_some_and(|c| same_circle(c, circle, tol)))
            .map(|(e, _)| e)
            .collect();
        let mut turn = 0.0;
        for piece in &pieces {
            let (_, (a, b)) = edge_curve(model, piece, tol)?;
            turn += (b - a).abs();
        }
        if (turn - core::f64::consts::TAU).abs() > tol.angular() {
            continue;
        }
        for piece in pieces {
            if !out.iter().any(|(e, _)| e.is_same(piece)) {
                out.push((piece.clone(), radius));
            }
        }
    }
    Ok(out)
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
