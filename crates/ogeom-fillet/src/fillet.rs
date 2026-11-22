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
    Seat, apply_wedge_solid, edge_curve, face_from_edges, planar_face, revolved_seat,
    segment_between, trimmed_wedge, wedge_solid,
};
use ogeom_algo::{Built, make_edge_between, make_vertex};
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, CylinderSurface, PlaneSurface};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

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
/// The blend engine builds the edge, as a request of one, as it takes many
/// in [`fillet_edges`]; a request it does not build is refused by name.
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
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    // One edge is a request of one: the engine first, as for many.
    if let Some(built) = rounds_by_engine(model, solid, &[(edge.clone(), radius)], tol)? {
        return Ok(built);
    }
    refused_alone(model, solid, edge, radius, Check::Ball, tol)?;
    Err(declined())
}

/// A request of two edges or more that no engine route builds, checked as
/// it is checked before any band of it is built: its corners of both
/// kinds and of mixed radii read, then each ball seated on both faces of
/// its edge and no two bands overlapping across a face they share
/// ([`crate::engine::request_fits`]), a concave band past its face's side kept for the
/// cut-off as there. One failing is refused here, by name, in that
/// check's words (ADR 0001, amendment "what M6 refuses by name"). A request
/// of one edge is [`refused_alone`]'s.
///
/// # Errors
///
/// The refusal.
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
    // The corners first, then the bands, in node order where a corner of
    // both kinds or a cusp is met.
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
    // Within the cut-offs: a concave band running past its face's side is
    // read as cut off at that side, not refused for it.
    let (fits, _) = crate::support::with_cut_offs(|| {
        crate::engine::request_fits(model, solid, &edges, &radii, tol)
    });
    fits.map_err(|b| b.error)?;
    Ok(())
}

/// Which check a request of one edge the engine does not build meets:
/// [`crate::engine::request_fits`] (a request through [`fillet_edges`]), or the ball on
/// its two faces alone ([`fillet_edge`]).
#[derive(Clone, Copy)]
enum Check {
    Request,
    Ball,
}

/// [`refused_as_checked`] for one edge the engine does not build: its ball
/// checked on its two faces, or its request as [`crate::engine::request_fits`]
/// reads it, and where that fails the refusal raised in those words, or in
/// the words of the engine's reading of a ball on a narrow face's far edge
/// where that names why its run does not hold.
///
/// # Errors
///
/// The refusal.
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
    match checked {
        // A ball resting on a narrow face's far edge names why its run does
        // not hold, where the engine's reading of it says.
        Err(refusal @ OgeomError::Construction(_)) => {
            match crate::engine::rolled::rolled_run(model, solid, edge, radius, tol) {
                Err(named @ OgeomError::Construction(_)) => Err(named),
                Err(other) => Err(other),
                Ok(_) => Err(refusal),
            }
        }
        _ => Ok(()),
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
    Ruled,
    /// A circle between a planar cap and its coaxial wall, a cylinder or
    /// a cone.
    Revolved(Box<crate::support::RevolvedSeat>),
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
    let (curve, _) = edge_curve(model, edge, tol)?;
    Ok(match curve {
        Curve::Line(_) if crate::support::hosts_planar(model, solid, edge, tol)? => Route::Planar,
        Curve::Line(_) => match crate::engine::ruled_seat::ruled_seat(model, solid, edge, tol)? {
            Some(_) => Route::Ruled,
            None => Route::Marched,
        },
        Curve::Circle(c) => match revolved_seat(model, solid, edge, &c, tol) {
            Ok(seat) if seat.slope != 0.0 && rims == Rims::Fitted => Route::Marched,
            Ok(seat) => Route::Revolved(Box::new(seat)),
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
/// A compound of solids is blended solid by solid, as mainstream CAD
/// blends each body on its own: each solid with the edges it holds, a
/// solid with none handed back as it was, the results gathered into a
/// compound again. An edge two of its solids share is refused by name.
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
///
/// Where two rounds of the request meet tangent at a vertex on one shared
/// face, their other faces meeting G1 along a seam there, the ball changes
/// walls where its contact crosses the seam at its own height. Where that
/// stands off the section square to the chain at the vertex (a band
/// leaning to the face its foot is rounded on, an arc's torus band meeting
/// a plane along a circle), the rounds joined in that section would leave
/// a gap or a step. Where a link's crossing at one end stands so over its
/// neighbour, and at its other end so or in that vertex's section (a half
/// frustum's foot or top round its leaning uprights' bands, a half
/// cylinder's foot round its end arcs' bands, a rib's root on a drum round
/// its rounded top), the link is rounded as the pipe about the ball's
/// centre curve between the two crossings, fitted through the ball's exact
/// sections, and its neighbours end on those sections; elsewhere the
/// request is refused by name, naming the vertex.
///
/// Where a round reaches geometry other than its own two faces, the
/// behaviours of [`Overflow`] are tried in its order, every one on:
/// [`fillet_edges_overflow`] chooses among them.
pub fn fillet_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    fillet_edges_overflow(model, solid, edges, radius, Overflow::default(), tol)
}

/// What a round does where it reaches geometry other than its own two
/// faces (ADR 0001, amendment "overflow handling"). The behaviours are
/// tried in the order of the fields, and the first that builds wins; all
/// are on by default.
///
/// - `blend_into_round`: where the round reaches an existing round of the
///   body, its surface bends to meet that round tangentially, with no step
///   and no crease. Built for a boss's foot: the closed rim between a
///   planar floor and a wall about an axis square to it (a cylinder, or a
///   cone narrowing away from the floor), whose round reaches straight
///   rounds tangent to the floor, convex or concave, over arcs apart. The
///   rim may be an open arc whose ends stand on planes through the axis,
///   the rounds it reaches inside it: each end's plane caps the round.
/// - `hold_edge`: where the ball would roll over a sharp edge of the body,
///   the edge is kept and the ball rests on it: it touches the other face
///   and the round there is the canal surface it sweeps, meeting the band
///   G1. Where one held edge stops and the next goes on, the ball rests on
///   the vertex between them. `hold_edge_at_end` applies the same rule
///   where a stretch held reaches an end of the band. Built for a round of
///   a straight edge, or a chain of edges along one line, between two
///   planes, ending on end faces within 60 degrees of square to it, whose
///   contact runs past lines or circles bounding one plane's faces; any
///   other request goes on to `trim_through`.
/// - `trim_through`: the round keeps its own surface, and the faces it
///   reaches are trimmed or extended to meet it.
///
/// A behaviour switched off passes the request to the next. Where none
/// that is on builds a round reaching other geometry, the request is
/// refused by name.
///
/// A **narrow passage** is a concave zone tighter than the radius, which
/// the ball cannot roll through touching its seats on both sides: two
/// walls of the floor no farther apart than twice the radius, or a corner
/// face whose curvature is tighter than the radius. Each round there keeps
/// its constant section and runs on through the zone, and the rounds meet
/// each other where they cross. Three strategies, cheapest first, say what
/// the tight zone may be made of; each takes what the ones before it take:
///
/// - `narrow_blend_made`: constant-radius rounds only (the request's own
///   rounds over a narrow floor, corner faces tangent to both walls);
/// - `narrow_whole_face`: also one whole corner face no round made (a
///   cylinder meeting its walls at kinks);
/// - `narrow_partial_face`: also a zone covering part of a face. Nothing
///   needs it alone yet: such a zone goes on to the engine, which refuses
///   it by name.
///
/// Built for a planar floor's rim with planar walls square to it, joined
/// at square corners (sharp or through a corner face) or by a half
/// cylinder between facing walls, and open ends flush on end faces. A
/// passage whose strategies are all off is refused by name.
///
/// A chamfer takes the same options
/// ([`chamfer_edges_overflow`](crate::chamfer_edges_overflow)): blend into
/// round and trim through apply to a flat bevel, and `hold_edge` and
/// `hold_edge_at_end` are ignored, a bevel having no ball to rest on an
/// edge. The narrow-passage switches are ignored too: a bevel in a narrow
/// passage goes to the engine as before.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct Overflow {
    /// Bend the round into an existing round it reaches.
    pub blend_into_round: bool,
    /// Keep a sharp edge the ball would roll over.
    pub hold_edge: bool,
    /// Keep a sharp edge where only a band's end runs over it.
    pub hold_edge_at_end: bool,
    /// Keep the round's surface and trim or extend what it reaches.
    pub trim_through: bool,
    /// Round a narrow passage made of constant-radius rounds only.
    pub narrow_blend_made: bool,
    /// Round a narrow passage spanning one whole face.
    pub narrow_whole_face: bool,
    /// Round a narrow passage covering part of a face.
    pub narrow_partial_face: bool,
}

impl Default for Overflow {
    fn default() -> Self {
        Self {
            blend_into_round: true,
            hold_edge: true,
            hold_edge_at_end: true,
            trim_through: true,
            narrow_blend_made: true,
            narrow_whole_face: true,
            narrow_partial_face: true,
        }
    }
}

impl Overflow {
    /// Every behaviour off: a round reaching other geometry is refused.
    #[must_use]
    pub const fn none() -> Self {
        Self {
            blend_into_round: false,
            hold_edge: false,
            hold_edge_at_end: false,
            trim_through: false,
            narrow_blend_made: false,
            narrow_whole_face: false,
            narrow_partial_face: false,
        }
    }
}

/// [`fillet_edges`], a round reaching geometry other than its own two
/// faces handled as `overflow` says.
///
/// With [`Overflow::default`] this is [`fillet_edges`]. Blend into round
/// is tried first where it is on, then hold edge. A request neither builds
/// goes on to the engine's own construction, which trims through where a
/// round reaches another face. With `trim_through` off, a result that
/// changes any face other than the edges' own and those at their ends is
/// refused by name, as is a request that blend into round or hold edge
/// would take while it is off.
///
/// # Errors
///
/// As [`fillet_edges`], and [`OgeomError::Construction`], by name, where a
/// round reaches other geometry and no behaviour that is on builds it.
pub fn fillet_edges_overflow(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let asked: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), radius)).collect();
    if let Some(built) = piped_or_refused(model, solid, &asked, tol)? {
        return Ok(built);
    }
    if let Some(built) = before_trim_through(model, solid, edges, radius, overflow, tol)? {
        return Ok(built);
    }
    let built = rounded_as_today(model, solid, edges, radius, tol);
    trimmed_through(model, solid, edges, built, overflow, tol)
}

/// The behaviours of `overflow` tried before trim through on a round of
/// `edges` at one `radius`: a narrow passage, blend into round, hold edge;
/// the first that builds, or `None` where the request goes on to trim
/// through.
///
/// # Errors
///
/// [`OgeomError::Construction`], by name, where one of them would take the
/// request, it is off, and so is every behaviour after it that could.
fn before_trim_through(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if model.kind_of(solid)? == ShapeType::Solid
        && let Some(narrow) = crate::engine::narrow_reach(model, solid, edges, radius, tol)?
    {
        let blend_made = narrow.needs == crate::engine::NarrowNeeds::BlendMade;
        let general = overflow.narrow_whole_face || overflow.narrow_partial_face;
        if (blend_made && overflow.narrow_blend_made) || general {
            if let Some(built) = crate::engine::narrow_built(model, solid, edges, &narrow, tol)? {
                return Ok(Some(built));
            }
        } else if blend_made {
            ogeom_bail!(
                Construction,
                "the round passes a narrow passage, a concave zone tighter than \
                 the radius, and every narrow-passage strategy is off; see \
                 docs/CAPABILITIES.md, fillet.overflow"
            );
        } else {
            ogeom_bail!(
                Construction,
                "the round passes a narrow passage spanning a face no \
                 constant-radius round made, and whole face and partial face \
                 are off; see docs/CAPABILITIES.md, fillet.overflow"
            );
        }
    }
    if let Some(built) = blend_into_round_step(model, solid, edges, radius, overflow, tol)? {
        return Ok(Some(built));
    }
    if model.kind_of(solid)? == ShapeType::Solid
        && let Some(hold) = crate::engine::hold_edge_reach(model, solid, edges, radius, tol)?
    {
        let wanted = overflow.hold_edge && (overflow.hold_edge_at_end || !hold.at_end);
        if wanted {
            if let Some(built) = crate::engine::hold_edge_held(model, solid, edges, &hold, tol)? {
                return Ok(Some(built));
            }
        } else if !overflow.trim_through {
            let which = if overflow.hold_edge {
                "hold edge at the end is off"
            } else {
                "hold edge is off"
            };
            ogeom_bail!(
                Construction,
                "the round's ball would roll over a sharp edge of the body, \
                 {which} and trim through is off; see docs/CAPABILITIES.md, \
                 fillet.overflow"
            );
        }
    }
    Ok(None)
}

/// Blend into round alone of `overflow` on a round of `edges` at one
/// `radius`: the foot blended where it is on and builds; `None` where the
/// request goes on.
///
/// # Errors
///
/// [`OgeomError::Construction`], by name, where blend into round would
/// take the request and both it and trim through are off.
pub(crate) fn blend_into_round_step(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if overflow.blend_into_round {
        if let Some(built) = blended_into_round(model, solid, edges, radius, tol)? {
            return Ok(Some(built));
        }
    } else if !overflow.trim_through
        && crate::engine::overflow_reach(
            model,
            solid,
            edges,
            crate::engine::OverflowSize::Round(radius),
            tol,
        )?
    {
        ogeom_bail!(
            Construction,
            "the round reaches an existing round of the body, and both blend into \
             round and trim through are off; see docs/CAPABILITIES.md, \
             fillet.overflow"
        );
    }
    Ok(None)
}

/// Trim through's result `built` for the request of `edges`: as built
/// where trim through is on; with it off, refused by name where it changed
/// a face other than the edges' own and those at their ends.
///
/// # Errors
///
/// The error `built` holds, and that refusal.
pub(crate) fn trimmed_through(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    built: OgeomResult<Built>,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let built = built?;
    if !overflow.trim_through
        && crate::engine::overflow_changes_others(model, solid, edges, &built, tol)?
    {
        ogeom_bail!(
            Construction,
            "the round reaches a face it is not seated on, which only trim \
             through builds here, and trim through is off; see \
             docs/CAPABILITIES.md, fillet.overflow"
        );
    }
    Ok(built)
}

/// The foot round of [`Overflow::blend_into_round`], on a solid or on
/// each solid of a compound; `None` where it does not take the request.
fn blended_into_round(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if let Some(built) = crate::compound::per_solid(
        model,
        solid,
        edges,
        |edge| edge,
        |model, part, own| {
            crate::engine::overflow_blended(
                model,
                part,
                own,
                crate::engine::OverflowSize::Round(radius),
                tol,
            )
        },
    )? {
        return Ok(Some(built));
    }
    if model.kind_of(solid)? != ShapeType::Solid {
        return Ok(None);
    }
    crate::engine::overflow_blended(
        model,
        solid,
        edges,
        crate::engine::OverflowSize::Round(radius),
        tol,
    )
}

/// The engine's rounds of `edges`, trimming through where they reach
/// another face, or the refusal by name of what it does not build.
fn rounded_as_today(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let asked: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), radius)).collect();
    if let Some(built) = rounds_by_engine(model, solid, &asked, tol)? {
        return Ok(built);
    }
    let radii = vec![radius; edges.len()];
    refused_as_checked(model, solid, edges, &radii, tol)?;
    // Concave rounds the engine declines on the solid's side are the
    // void's convex rounds, which it may build.
    if let Some(built) = crate::engine::void_side_rounds(model, solid, &asked, tol)? {
        return crate::engine::merged(model, solid, built, tol);
    }
    Err(declined())
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
/// Where the rounds meeting at a vertex are all of one sense, on planes or
/// curved faces, the corner [`fillet_edges`] builds is already the ball's
/// and is taken: bands running on through a tangent junction with no
/// corner, the sphere where one ball touches every face (a pyramid's apex
/// of four rounds too), bands crossing beside an edge left sharp. The
/// torus round a lone band's end is built between planes, and on a curved
/// face where the lone edge is a line and the faces beside it planes,
/// cylinders along it or cones it is a ruling of: there the lone bands are
/// rounded first and the rest roll round their feet. Square to a plane
/// across (a D-shaped boss's foot with its uprights) that roll is the
/// torus about the band's axis; leaning to it (a half frustum's uprights
/// on a plate), or with a curved face across (a rib's top under a drum's),
/// the ball changes walls where its contact crosses the band's rails, and
/// the round between is the pipe [`fillet_edges`] builds there.
///
/// The result does not depend on the order of `edges`.
///
/// # Errors
///
/// As [`fillet_edges`] for `Setback`, and for `Rolling` where every vertex
/// is of one sense. For `Rolling`, [`OgeomError::Construction`], naming
/// the vertex, where the ball's sweep does not close it: a vertex mixing
/// convex and concave rounds that joins other than three edges, or that
/// touches a curved face where the lone round's edge is curved or its
/// faces are not planes, cylinders along it or cones it is a ruling of;
/// a convex and a concave round
/// meeting beside an edge left sharp; and as [`fillet_edges`] where the
/// rounds rolling round a leaning band's foot are refused. A request
/// mixing senses that the engine's rounds between planes do not build is
/// refused, never closed by setbacks instead.
///
/// Where a round reaches geometry other than its own two faces, the
/// behaviours of [`Overflow`] are tried in its order, every one on:
/// [`fillet_edges_corner_overflow`] chooses among them.
pub fn fillet_edges_corner(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: FilletCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    fillet_edges_corner_overflow(
        model,
        solid,
        edges,
        radius,
        corner,
        Overflow::default(),
        tol,
    )
}

/// [`fillet_edges_corner`], a round reaching geometry other than its own
/// two faces handled as `overflow` says.
///
/// With [`Overflow::default`] this is [`fillet_edges_corner`], and
/// [`FilletCorner::Setback`] is [`fillet_edges_overflow`]. The corner and
/// the overflow are two settings: the corner says how rounds meeting at a
/// vertex close it, the overflow what a round does where it reaches other
/// geometry. For [`FilletCorner::Rolling`], every request whose corners
/// [`fillet_edges`] builds already (a vertex of one sense, or none) is
/// [`fillet_edges_overflow`] with `overflow`, each behaviour off passing
/// to the next as there; at a vertex mixing convex and concave rounds,
/// each pass of rounds is that call too, and the corner closed between
/// planes by the ball's sphere or torus is trim through's construction:
/// blend into round, hold edge and the narrow passages are not built
/// there, and with trim through off, a result changing any face other than
/// the edges' own and those at their ends is refused by name.
///
/// # Errors
///
/// As [`fillet_edges_corner`], and as [`fillet_edges_overflow`] where a
/// round reaches other geometry and no behaviour that is on builds it.
pub fn fillet_edges_corner_overflow(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: FilletCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    match corner {
        FilletCorner::Setback => fillet_edges_overflow(model, solid, edges, radius, overflow, tol),
        FilletCorner::Rolling => {
            if let Some(built) = crate::compound::each_solid(
                model,
                solid,
                edges,
                |edge| edge,
                |model, part, own| {
                    crate::engine::rolling_rounds(model, part, own, radius, overflow, tol)
                },
            )? {
                return Ok(built);
            }
            crate::engine::rolling_rounds(model, solid, edges, radius, overflow, tol)
        }
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
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
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
    if let Some(built) = piped_or_refused(model, solid, edges, tol)? {
        return Ok(built);
    }
    if let Some(built) = rounds_by_engine(model, solid, edges, tol)? {
        return Ok(built);
    }
    refused_as_checked(model, solid, &shapes, &radii, tol)?;
    Err(declined())
}

/// A chain whose ball changes walls off a vertex where its rounds meet
/// (ADR 0001, amendment "a ball changing walls off the vertex"), its
/// links between two such crossings built as pipes; `None` where no
/// crossing stands off its vertex.
///
/// # Errors
///
/// [`OgeomError::Construction`], by name, where one does and no pipe
/// closes it: the rounds, joined in the section square to the chain at
/// the vertex, would leave a gap or an overlap beside the seam.
fn piped_or_refused(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if let Some(off) = crate::engine::switch_over::off_vertex(model, solid, edges, tol)? {
        if let Some(built) = crate::engine::piped_rounds(model, solid, edges, true, tol)? {
            return crate::engine::merged(model, solid, built, tol).map(Some);
        }
        let p = off.at;
        ogeom_bail!(
            Construction,
            "the ball changes walls off the vertex at ({:.6}, {:.6}, {:.6}): its contact \
             crosses the seam between the walls {:.3} along the chain from the vertex's \
             section, and the rounds joined there do not close the corner",
            p.x,
            p.y,
            p.z,
            off.by
        );
    }
    Ok(None)
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
/// `None` where no route builds it. An edge asked twice at
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
/// sheet's rounds are not the engine's), an edge lies in no solid, or any
/// solid's rounds are handed back.
///
/// # Errors
///
/// [`OgeomError::Construction`] where an edge lies in two of its solids
/// ([`crate::compound::per_solid`]).
fn rounds_by_engine_compound(
    model: &mut Model,
    compound: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    crate::compound::per_solid(
        model,
        compound,
        edges,
        |(edge, _)| edge,
        |model, part, own| rounds_by_engine(model, part, own, tol),
    )
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
    // Straight concave creases, each fill run on past both ends until the
    // end faces' surfaces trim it (amendment "a concave band ends against
    // its end face").
    if let Some(built) = crate::engine::concave_end_rounds(model, solid, &once, tol)? {
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
    // Straight creases whose marched bands stop short of an end face they
    // lean into: carried past their ends, each wedge clipped by that face
    // (amendment "a marched band leaning into its end face").
    if lines_only(model, &once, tol)?
        && let Some(built) = crate::marched::carried_past_ends(|| {
            crate::engine::marched_rounds(model, solid, &once, false, tol)
        })?
    {
        return Ok(Some(built));
    }
    // A cone rim asked alone whose exact torus band nothing above closes,
    // its band fitted through the ball's arcs as the march builds it: the
    // boolean resolves some meetings of
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

/// Whether every edge of `once` lies on a line.
fn lines_only(model: &Model, once: &[(Shape, f64)], tol: Tolerances) -> OgeomResult<bool> {
    for (edge, _) in once {
        let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
        if !matches!(curve, ogeom_geom::Curve::Line(_)) {
            return Ok(false);
        }
    }
    Ok(true)
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
    // Nothing asked, or a radius that rounds nothing, is the entry's to
    // refuse by name.
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
        let whole = whole_rims(model, solid, once, tol)?;
        if let Some(built) = crate::engine::round_rims(model, solid, &whole, tol)? {
            return Ok(Some(built));
        }
        if let Some(built) = crate::engine::round_curved(model, solid, &whole, tol)? {
            return Ok(Some(built));
        }
        // A cusp's S asked on one side of the tangent point, the rim's arc
        // past it left sharp: the arc as asked (amendment "a cusp's S, one
        // side or two").
        if whole.len() > once.len()
            && !crate::setback::cusp_corners(model, solid, once, tol)?.is_empty()
        {
            return crate::engine::round_curved(model, solid, once, tol);
        }
        Ok(None)
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
    // the solid's own nodes, and a request the engine hands back leaves
    // the model as it came.
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

/// Whether any of `edges` is a rim of `solid` on a planar cap and a
/// coaxial cone.
fn has_cone_rim(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<bool> {
    for edge in edges {
        if let Route::Revolved(seat) = route(model, solid, edge, Rims::Exact, tol)?
            && seat.slope != 0.0
        {
            return Ok(true);
        }
    }
    Ok(false)
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
    if let Some(built) = crate::compound::each_solid(
        model,
        solid,
        std::slice::from_ref(edge),
        |edge| edge,
        |model, part, _| fillet_edge_variable(model, part, edge, start_radius, end_radius, tol),
    )? {
        return Ok(built);
    }
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
///
/// Where a round reaches geometry other than its own two faces, the
/// behaviours of [`Overflow`] are tried in its order, every one on:
/// [`fillet_edges_variable_overflow`] chooses among them.
pub fn fillet_edges_variable(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Built> {
    fillet_edges_variable_overflow(model, solid, chain, law, Overflow::default(), tol)
}

/// [`fillet_edges_variable`], a round reaching geometry other than its
/// own two faces handled as `overflow` says.
///
/// With [`Overflow::default`] this is [`fillet_edges_variable`].
///
/// A law holding one radius all along is the constant round of that
/// radius, and takes every behaviour of [`fillet_edges_overflow`] in its
/// order: a narrow passage, blend into round (a boss's foot bending into
/// the edge round it reaches), hold edge, then trim through, which is the
/// variable round's own construction here, or, where that refuses, the
/// constant round's.
///
/// A law whose radius changes along the chain takes trim through alone:
/// its round keeps its own surface where it reaches other geometry, as it
/// did before this option set. Blend into round, hold edge and the narrow
/// passages are not built for it, so switching them on changes nothing and
/// switching them off refuses nothing; with trim through off, a result
/// changing any face other than the chain's own and those at its ends is
/// refused by name.
///
/// # Errors
///
/// As [`fillet_edges_variable`], and [`OgeomError::Construction`], by
/// name, where a round reaches other geometry and no behaviour that is on
/// builds it.
pub fn fillet_edges_variable_overflow(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if chain.is_empty() {
        ogeom_bail!(
            Construction,
            "a variable radius over no edges rounds nothing"
        );
    }
    let level = Profile::new(law, false, tol)
        .ok()
        .filter(|p| p.level(tol))
        .map(|p| p.largest());
    let Some(radius) = level else {
        let built = varying_as_today(model, solid, chain, law, tol);
        return trimmed_through(model, solid, chain, built, overflow, tol);
    };
    if let Some(built) = before_trim_through(model, solid, chain, radius, overflow, tol)? {
        return Ok(built);
    }
    let built = varying_as_today(model, solid, chain, law, tol);
    trimmed_through(model, solid, chain, built, overflow, tol)
}

/// [`fillet_edges_variable`]'s construction, the round keeping its own
/// surface where it reaches other geometry.
fn varying_as_today(
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
    // A chain runs on one solid: on a compound, the solid holding it.
    let mut solids = 0;
    if let Some(built) = crate::compound::each_solid(
        model,
        solid,
        chain,
        |edge| edge,
        |model, part, own| {
            solids += 1;
            if solids > 1 {
                ogeom_bail!(
                    Construction,
                    "a variable radius runs along one chain on one solid; these edges \
                     lie on two solids of the compound"
                );
            }
            varying_as_today(model, part, own, law, tol)
        },
    )? {
        return Ok(built);
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
/// (ADR 0001, "what M6 refuses by name"), and a rim that closes is rounded
/// the whole turn; a rim that stops, or pieces asked at
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
pub(crate) fn same_circle(a: &ogeom_math::Circle, b: &ogeom_math::Circle, tol: Tolerances) -> bool {
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
