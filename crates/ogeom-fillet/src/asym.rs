//! Asymmetric rounds: a blend that meets one face a distance `d1` from the
//! edge and the other `d2`, its section a conic tangent to both faces there.
//!
//! In the plane square to a straight crease, with `A` the crease and `P1`,
//! `P2` the two contacts, the section is the rational quadratic Bézier arc
//! with control points `P1`, `A`, `P2` and weights `1`, `w`, `1`. It is
//! tangent to `P1A` at `P1` and to `AP2` at `P2`, so the band meets both
//! faces with their tangent planes along its contact lines; and the arc's
//! middle stands on the line from the chord's middle `M` to `A`, a share
//! `rho = w / (1 + w)` of the way from `M`. A conic of shape `rho` takes it
//! as asked; the elliptical profile takes the weight of the circular round
//! of the same dihedral, `w = sin(alpha / 2)` with `alpha` the angle at `A`
//! between the faces, which makes the arc the image of that round's arc
//! under the affine map fixing `A` and stretching each face's leg to its
//! own distance: an ellipse, and with equal distances the circle itself.
//!
//! The band along a straight crease is that arc swept along the crease: an
//! extrusion of a rational B-spline, exact, which a plane crosses in closed
//! form. Everything else about the seat (where its ends run out, a face
//! overhanging its corner, the wedge cut or fused) is the circular round's
//! own, with the section's corners and middle standing for the ball's.

use crate::support::Seat;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSplineCurve, Curve};
use ogeom_math::bspline::Weighted;
use ogeom_math::{Direction, KnotVector, Point, Vector};
use ogeom_topo::{Model, Shape};

/// The shape of an asymmetric round's section, between its two contacts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RoundProfile {
    /// The ellipse tangent to each face at its distance: the circular
    /// round's section stretched along each face to that face's distance.
    /// With equal distances it is the circular round of
    /// [`fillet_edge`](crate::fillet_edge).
    Elliptical,
    /// The conic tangent to each face at its distance whose middle stands
    /// the share `rho` of the way from the middle of the chord between the
    /// contacts to the crease: near `0` a bevel, `0.5` a parabola, near `1`
    /// the sharp edge. `rho` lies strictly between `0` and `1`.
    Conic(f64),
}

impl RoundProfile {
    /// The middle control point's weight for a section whose faces' legs
    /// from the crease make the cosine `cos_alpha` with each other.
    pub(crate) fn weight(self, cos_alpha: f64) -> OgeomResult<f64> {
        match self {
            Self::Elliptical => Ok(((1.0 - cos_alpha) * 0.5).max(0.0).sqrt()),
            Self::Conic(rho) => {
                if !(rho > 0.0 && rho < 1.0) {
                    ogeom_bail!(
                        Construction,
                        "a conic round's rho of {rho} is not strictly between 0 and 1"
                    );
                }
                Ok(rho / (1.0 - rho))
            }
        }
    }
}

/// The conic arc a straight seat's band sweeps along its crease, tangent
/// to each face of the seat at its distance from the crease.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Section {
    /// The distance along each face of the seat, in the seat's face order.
    pub(crate) d: [f64; 2],
}

/// On each face of the seat, the unit direction square to the crease that
/// runs along the face away from it.
pub(crate) fn legs(seat: &Seat, tol: Tolerances) -> OgeomResult<[Vector; 2]> {
    let sign = if seat.convex { 1.0 } else { -1.0 };
    Ok([seat.leg(0, tol)? * sign, seat.leg(1, tol)? * sign])
}

/// How far from 1 a piece of a conic arc's weight may stand before the arc
/// is split there: a piece nearer the polynomial arc of weight 1 keeps the
/// rules that integrate and intersect it near their polynomial accuracy.
/// The circle's quarter, of weight `1/sqrt(2)`, stays whole.
const WEIGHT_RANGE: (f64, f64) = (0.7, 1.25);

/// The conic arc with control points `p`, weight `w` on the middle one, as
/// a rational B-spline on `[0, 1]`: the arc split at its middle, each half
/// again, until every piece's weight in standard form lies within
/// [`WEIGHT_RANGE`]. Splitting the standard arc at its middle gives two
/// standard arcs of weight `sqrt((1 + w) / 2)` each, the same curve; the
/// pieces meet at double knots.
pub(crate) fn conic_arc(p: [Point; 3], w: f64, tol: Tolerances) -> OgeomResult<Curve> {
    let (control, weight, knots) = conic_pieces(p, w);
    let control = control
        .into_iter()
        .enumerate()
        .map(|(i, q)| Weighted::new(q, if i % 2 == 1 { weight } else { 1.0 }, tol))
        .collect::<OgeomResult<Vec<_>>>()?;
    Ok(Curve::BSpline(BSplineCurve::rational(
        KnotVector::new(knots, 2)?,
        control,
    )?))
}

/// The pieces [`conic_arc`] splits its arc into: their control points end
/// to end (each piece's first, middle and, for the last, its end), the
/// middle points' common weight, the end points' being 1, and the knots
/// on `[0, 1]`. Each control point is an affine combination of `p` whose
/// coefficients depend on `w` alone, so a family of arcs whose control
/// points move affinely splits into control points that move alike.
pub(crate) fn conic_pieces(p: [Point; 3], w: f64) -> (Vec<Point>, f64, Vec<f64>) {
    let mut pieces = vec![p];
    let mut weight = w;
    for _ in 0..8 {
        if weight >= WEIGHT_RANGE.0 && weight <= WEIGHT_RANGE.1 {
            break;
        }
        pieces = pieces
            .into_iter()
            .flat_map(|[p0, a, p2]| {
                let near =
                    Point::from_vector((p0.to_vector() + a.to_vector() * weight) / (1.0 + weight));
                let far =
                    Point::from_vector((a.to_vector() * weight + p2.to_vector()) / (1.0 + weight));
                let middle = near.midpoint(far);
                [[p0, near, middle], [middle, far, p2]]
            })
            .collect();
        weight = ((1.0 + weight) * 0.5).sqrt();
    }
    let n = pieces.len();
    let mut knots = vec![0.0; 3];
    let mut control = Vec::with_capacity(2 * n + 1);
    for (i, [p0, a, _]) in pieces.iter().enumerate() {
        control.push(*p0);
        control.push(*a);
        if i + 1 < n {
            #[allow(clippy::cast_precision_loss, reason = "at most 256 pieces")]
            let k = (i + 1) as f64 / n as f64;
            knots.extend([k, k]);
        }
    }
    control.push(pieces[n - 1][2]);
    knots.extend([1.0; 3]);
    (control, weight, knots)
}

/// Round a straight edge between two planar faces with a band that meets
/// `face` `on_face` from the edge and the edge's other face `on_other`
/// from it, its section `profile` between the two contacts.
///
/// The band is tangent to each face along its contact line. Naming the
/// other face with the distances swapped builds the same solid, to the
/// bit. With equal distances the elliptical profile is the circular round
/// [`fillet_edge`](crate::fillet_edge) builds on the edge, by the radius
/// whose ball touches each face that far from the edge: the same call.
/// Otherwise the band is the conic arc swept along the edge, exact, and
/// its ends are the circular round's: square caps where the edge ends at
/// a wall, run out across a face that cuts the corner away past an end,
/// trimmed by a face overhanging the corner.
///
/// What `solid` shares with other shapes is copied before it is worked on
/// ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where `face` is not one of the edge's two faces, a distance is not a
/// usable length, a conic's `rho` is not strictly between 0 and 1, the
/// edge is not straight between two planes, a contact sets back past
/// its face's side, or on a drum's ruling the section's conic crosses the
/// drum between its contacts; with equal distances and the elliptical
/// profile, as
/// [`fillet_edge`](crate::fillet_edge).
#[allow(
    clippy::too_many_arguments,
    reason = "the chamfer's distances' own, and the profile"
)]
pub fn fillet_edge_distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    fillet_edge_distances_overflow(
        model,
        solid,
        edge,
        face,
        on_face,
        on_other,
        profile,
        crate::Overflow::default(),
        tol,
    )
}

/// [`fillet_edge_distances`], a round reaching geometry other than its own
/// two faces handled as `overflow` says.
///
/// With [`Overflow::default`](crate::Overflow::default) this is
/// [`fillet_edge_distances`]. Two equal distances of the elliptical
/// profile are the circular round, which takes blend into round as
/// [`fillet_edges_overflow`](crate::fillet_edges_overflow) does where the
/// edge is a boss's whole foot rim. Every other request takes trim through
/// alone, the construction [`fillet_edge_distances`] built before this
/// option set: blend into round is not built for two different distances
/// or a conic profile, nor hold edge and the narrow passages for this
/// call, so switching them on changes nothing and switching them off
/// refuses nothing. With trim through off, a result where the round
/// changes a face other than the edge's own and those at its ends, or
/// carries a face on past its side over nothing, is refused by name, as is
/// a foot blend into round would take while both are off.
///
/// # Errors
///
/// As [`fillet_edge_distances`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where the round reaches other geometry and no behaviour that is
/// on builds it.
#[allow(
    clippy::too_many_arguments,
    reason = "the chamfer's distances' own, the profile and the overflow"
)]
pub fn fillet_edge_distances_overflow(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let usable = |d: f64| d.is_finite() && d > tol.confusion();
    if profile == RoundProfile::Elliptical
        && on_face == on_other
        && usable(on_face)
        && model.kind_of(solid)? == ogeom_topo::ShapeType::Solid
        && constant_dihedral(model, solid, edge, tol)?
    {
        let radius = equal_round(model, solid, edge, face, on_face, tol)?;
        let one = std::slice::from_ref(edge);
        if let Some(built) =
            crate::fillet::blend_into_round_step(model, solid, one, radius, overflow, tol)?
        {
            return Ok(built);
        }
    }
    let built = distances_as_today(model, solid, edge, face, on_face, on_other, profile, tol);
    let sized = [(edge.clone(), on_face.max(on_other))];
    crate::fillet::trimmed_through(model, solid, &sized, built, overflow, tol)
}

/// [`fillet_edge_distances`]' construction, the round keeping its own
/// surface where it reaches other geometry.
#[allow(
    clippy::too_many_arguments,
    reason = "the chamfer's distances' own, and the profile"
)]
fn distances_as_today(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for d in [on_face, on_other] {
        if !d.is_finite() || d <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "a round set back {d} from its edge rounds nothing"
            );
        }
    }
    if let RoundProfile::Conic(_) = profile {
        profile.weight(0.0)?;
    }
    if let Some(built) = crate::compound::each_solid(
        model,
        solid,
        std::slice::from_ref(edge),
        |edge| edge,
        |model, part, _| {
            distances_as_today(model, part, edge, face, on_face, on_other, profile, tol)
        },
    )? {
        return Ok(built);
    }
    if profile == RoundProfile::Elliptical
        && on_face == on_other
        && constant_dihedral(model, solid, edge, tol)?
    {
        // The circular round, by the same call.
        let radius = equal_round(model, solid, edge, face, on_face, tol)?;
        return crate::fillet_edge(model, solid, edge, radius, tol);
    }
    // The blend engine first, as for a round of one radius: a straight edge
    // between two planes on a solid bounded by planes.
    let asked = Asymmetric {
        edge: edge.clone(),
        face: face.clone(),
        on_face,
        on_other,
        profile,
    };
    if let Some(request) = crate::engine::Request::from_round_distances(model, solid, &asked, tol)?
        && let Some(built) = crate::engine::run(model, solid, &request, tol)?
    {
        return Ok(built);
    }
    // A rim between a planar cap and a coaxial wall: its conic turned about
    // the axis, through the engine's rims, over every arc of its circle; a
    // rim it does not build is handed on.
    let along = ogeom_topo::ancestors_of(model, solid, edge, ogeom_topo::ShapeType::Face)?;
    if along.iter().any(|f| f.is_same(face))
        && let Some(arcs) = rim_arcs(model, solid, edge, tol)?
    {
        let cap_named = ogeom_topo::Model::node(model, face)
            .and_then(|n| n.data().as_face().map(|d| d.surface))
            .and_then(|s| model.geometry().surface(s))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        let (cap, wall) = if cap_named {
            (on_face, on_other)
        } else {
            (on_other, on_face)
        };
        let section = crate::engine::RimSection::Conic { cap, wall, profile };
        let asked: Vec<(Shape, crate::engine::RimSection)> =
            arcs.iter().map(|e| (e.clone(), section.clone())).collect();
        match crate::engine::round_rims_with(model, solid, &asked, tol) {
            Ok(Some(built)) => return Ok(built),
            Ok(None) | Err(ogeom_core::OgeomError::Construction(_)) => {}
            Err(e) => return Err(e),
        }
        // A rim the seams split (a wall a boolean left in pieces): its pieces
        // and faces joined, rounded there, the history carried back.
        if let Some(built) = conic_rim_dissolved(model, solid, &arcs, &section, tol)? {
            return Ok(built);
        }
    }
    // Any other seat: the band applied by the engine as a one-boolean tool.
    if let Some(built) =
        crate::engine::distances_rounds(model, solid, edge, face, on_face, on_other, profile, tol)?
    {
        return Ok(built);
    }
    Err(crate::fillet::declined())
}

/// The conic section of a straight edge between planes asked as `face`,
/// `on_face`, `on_other` and `profile`, its distances in the seat's face
/// order; `None` where it is the circular round's (equal distances and the
/// elliptical profile).
#[allow(clippy::too_many_arguments, reason = "the request's own")]
pub(crate) fn straight_section(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Option<Section>> {
    let (seat, named) = straight_seat(model, solid, edge, face, tol)?;
    if profile == RoundProfile::Elliptical && on_face == on_other {
        return Ok(None);
    }
    let mut d = [0.0; 2];
    d[named] = on_face;
    d[1 - named] = on_other;
    let f = legs(&seat, tol)?;
    // The conic's weight reads at the seat's opening, or the round is refused.
    profile.weight(f[0].dot(f[1]).clamp(-1.0, 1.0))?;
    Ok(Some(Section { d }))
}

/// The straight seat of an asymmetric round, and which of its faces is
/// `face`; refused by name where the edge is not straight between two
/// planes or `face` is not beside it.
pub(crate) fn straight_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Seat, usize)> {
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    if !matches!(curve, Curve::Line(_)) || !crate::support::hosts_planar(model, solid, edge, tol)? {
        ogeom_bail!(
            Construction,
            "an asymmetric round is built on a straight edge between two planes or a \
             rim where a planar cap meets a coaxial cylinder or cone; this edge is \
             neither; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
        );
    }
    let seat = crate::support::planar_seat(model, solid, edge, tol)?;
    let named = if seat.faces[0].is_same(face) {
        0
    } else if seat.faces[1].is_same(face) {
        1
    } else {
        ogeom_bail!(
            Construction,
            "the face an asymmetric round's first distance runs along is not one of the \
             edge's two faces"
        );
    };
    Ok((seat, named))
}

/// The radius of the circular round an elliptical round of `distance` on
/// each face is, refused as [`fillet_edge_distances`] refuses the request.
pub(crate) fn equal_round(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    if !distance.is_finite() || distance <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a round set back {distance} from its edge rounds nothing"
        );
    }
    let profile = RoundProfile::Elliptical;
    if let Some(rim) = rim_conic(model, solid, edge, face, distance, distance, profile, tol)? {
        return Ok(rim.equal_radius(distance));
    }
    straight_section(model, solid, edge, face, distance, distance, profile, tol)?;
    equal_radius(model, solid, edge, distance, tol)
}

/// The radius of the circular round whose ball touches each face of the
/// straight edge `distance` from it.
fn equal_radius(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let seat = crate::support::planar_seat(model, solid, edge, tol)?;
    let f = legs(&seat, tol)?;
    let cos = f[0].dot(f[1]).clamp(-1.0, 1.0);
    let sin = f[0].cross(f[1]).magnitude();
    if sin <= tol.angular() || 1.0 + cos <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the faces meet too sharply or too flat to round"
        );
    }
    Ok(distance * sin / (1.0 + cos))
}

/// Refuses by name where a contact of the conic band, read at a quarter,
/// the middle and three quarters of the edge, stands off its face and off
/// every other face of the solid: the band would set back past the face's
/// side.
pub(crate) fn contacts_fit(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    section: Section,
    tol: Tolerances,
) -> OgeomResult<()> {
    let Section { d } = section;
    let seat = crate::support::planar_seat(model, solid, edge, tol)?;
    let f = legs(&seat, tol)?;
    for fraction in [0.25, 0.5, 0.75] {
        let p = seat.start + (seat.end - seat.start) * fraction;
        for i in 0..2 {
            let contact = p + f[i] * d[i];
            if !crate::support::on_face_or_boundary(model, solid, &seat.faces[i], contact, tol)? {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round sets back {} along a face of the edge through \
                     {p:?}, past that face's side; see docs/CAPABILITIES.md, \
                     fillet.asymmetric-rounds",
                    d[i]
                );
            }
        }
    }
    Ok(())
}

/// One asymmetric round of a [`blend_edges`](crate::blend_edges) request.
#[derive(Clone)]
pub(crate) struct Asymmetric {
    /// The edge.
    pub edge: Shape,
    /// The face the first distance runs along.
    pub face: Shape,
    /// The distance along `face`.
    pub on_face: f64,
    /// The distance along the other face.
    pub on_other: f64,
    /// The section.
    pub profile: RoundProfile,
}

/// Whether two asymmetric rounds meeting at `vertex` continue each other:
/// the same face named, the same distances and section, and their edges
/// tangent there.
pub(crate) fn continues(
    model: &Model,
    a: &Asymmetric,
    b: &Asymmetric,
    vertex: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    if !a.face.is_same(&b.face)
        || a.on_face != b.on_face
        || a.on_other != b.on_other
        || a.profile != b.profile
    {
        return Ok(false);
    }
    let at = {
        let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
            return Ok(false);
        };
        vertex.transform(model.datums())?.apply(data.point)
    };
    let tangent = |edge: &Shape| -> OgeomResult<Vector> {
        use ogeom_geom::Curve3d as _;
        let (curve, (t0, t1)) = crate::support::edge_curve(model, edge, tol)?;
        let t = if curve.point_at(t0, tol)?.distance(at) <= curve.point_at(t1, tol)?.distance(at) {
            t0
        } else {
            t1
        };
        let d = curve.d1_at(t, tol)?;
        Ok(d / d.magnitude())
    };
    Ok(tangent(&a.edge)?.cross(tangent(&b.edge)?).magnitude() <= tol.angular() * 10.0)
}

/// An asymmetric round's section on a rim where a planar cap meets a
/// coaxial cylinder or cone, in the meridian half-plane through the seat's
/// reference direction.
pub(crate) struct RimConic {
    /// The rim.
    pub(crate) seat: crate::support::RevolvedSeat,
    /// The edge's circle.
    pub(crate) circle: ogeom_geom::CircleCurve,
    /// Whether the rims on the circle run its whole turn.
    pub(crate) whole: bool,
    /// The wall's contact's depth along the axis from the cap's plane.
    pub(crate) wall_depth: f64,
    /// The cap's contact ring's radius.
    pub(crate) cap_rho: f64,
    /// The section's control points: the wall's contact, the rim, the
    /// cap's contact.
    pub(crate) points: [Point; 3],
    /// The weight of the rim's control point.
    pub(crate) w: f64,
}

impl RimConic {
    /// How far the section reaches from the rim: the farther distance.
    pub(crate) fn size(&self) -> f64 {
        let [wall, rim, cap] = self.points;
        wall.distance(rim).max(cap.distance(rim))
    }

    /// The contacts in the meridian, as a ball's would stand: the cap's
    /// ring's radius and the wall's ring's radius and height. The ball's
    /// own centre height and arc sweep stand for nothing here.
    pub(crate) fn contacts(
        &self,
        seat: &crate::support::RevolvedSeat,
    ) -> crate::support::RevolvedBall {
        let wall_height = -seat.tau * self.wall_depth;
        crate::support::RevolvedBall {
            rho: self.cap_rho,
            height: 0.0,
            wall_rho: seat.radius_at(wall_height),
            wall_height,
            sweep: 0.0,
        }
    }

    /// The section in the meridian half-plane at `theta` from the seat's
    /// reference direction, from the wall's contact to the cap's.
    pub(crate) fn meridian(
        &self,
        seat: &crate::support::RevolvedSeat,
        theta: f64,
        tol: Tolerances,
    ) -> OgeomResult<Curve> {
        let turn = ogeom_math::Transform::rotation(
            ogeom_math::Axis {
                location: seat.centre,
                direction: Direction::new(seat.up, tol)?,
            },
            theta,
        );
        conic_arc(self.points.map(|p| turn.apply(p)), self.w, tol)
    }

    /// The band over the window of `sweep` from `theta0`: the section at
    /// `theta0` revolved through `sweep`, bounded by the wall's and the
    /// cap's contact arcs (framed from `theta0`, on `(0, sweep)`) and the
    /// meridian sections at both ends, each given its line in the band's
    /// chart, facing away from the rim.
    pub(crate) fn band_sector(
        &self,
        model: &mut Model,
        seat: &crate::support::RevolvedSeat,
        (theta0, sweep): (f64, f64),
        [wall, cap, start, end]: [&Shape; 4],
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Surface as _;
        let axis = ogeom_math::Axis {
            location: seat.centre,
            direction: Direction::new(seat.up, tol)?,
        };
        let profile = self.meridian(seat, theta0, tol)?;
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::RevolutionSurface::over(profile, axis, (0.0, sweep))?.into();
        let surface_id = model.geometry_mut().add_surface(surface.clone());
        let line = |from: (f64, f64),
                    along: (f64, f64),
                    to: f64|
         -> OgeomResult<ogeom_geom::PlanarCurve> {
            Ok(ogeom_geom::Line2d::over(
                ogeom_math::Axis2::new(
                    ogeom_math::Point2::new(from.0, from.1),
                    ogeom_math::Direction2::new(ogeom_math::Vector2::new(along.0, along.1), tol)?,
                ),
                0.0,
                to,
            )?
            .into())
        };
        for (edge, pcurve, range) in [
            (wall, line((0.0, 0.0), (1.0, 0.0), sweep)?, (0.0, sweep)),
            (cap, line((0.0, 1.0), (1.0, 0.0), sweep)?, (0.0, sweep)),
            (start, line((0.0, 0.0), (0.0, 1.0), 1.0)?, (0.0, 1.0)),
            (end, line((sweep, 0.0), (0.0, 1.0), 1.0)?, (0.0, 1.0)),
        ] {
            ogeom_algo::attach_pcurve(
                model,
                edge,
                pcurve,
                surface_id,
                ogeom_topo::Location::identity(),
                range,
            )?;
        }
        let wire = ogeom_algo::make_wire(
            model,
            &[wall.clone(), end.clone(), cap.reversed(), start.reversed()],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        let middle = surface.point_at(sweep * 0.5, 0.5, tol)?;
        let (du, dv) = surface.d1_at(sweep * 0.5, 0.5, tol)?;
        let rim = ogeom_math::Transform::rotation(axis, theta0 + sweep * 0.5).apply(self.points[1]);
        Ok(if du.cross(dv).dot(middle - rim) >= 0.0 {
            face
        } else {
            face.reversed()
        })
    }
    /// The radius of the circular round whose ball touches cap and wall
    /// `distance` from the rim: the section's legs meet at the angle
    /// `alpha` there, and the ball's tangents from the rim are
    /// `r cot(alpha / 2)` long.
    fn equal_radius(&self, distance: f64) -> f64 {
        let [wall, rim, cap] = self.points;
        let (a, b) = (wall - rim, cap - rim);
        let (a, b) = (a / a.magnitude(), b / b.magnitude());
        let cos = a.dot(b).clamp(-1.0, 1.0);
        let sin = a.cross(b).magnitude();
        distance * sin / (1.0 + cos)
    }
}

/// The asymmetric round's section on a rim, where `edge` is one: `None`
/// where it is not a circle; refused by name where it is a circle but not
/// a whole turn of a rim on a cap and its coaxial wall, `face` is neither,
/// or a contact would reach the axis.
#[allow(clippy::too_many_arguments, reason = "the request's own")]
pub(crate) fn rim_conic(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Option<RimConic>> {
    let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
    let Curve::Circle(circle) = curve else {
        return Ok(None);
    };
    let seat = match crate::support::revolved_seat(model, solid, edge, &circle, tol) {
        Ok(seat) => seat,
        // A circle that is no rim of a cap and its coaxial wall is a
        // marched crease.
        Err(ogeom_core::OgeomError::Construction(_)) => return Ok(None),
        Err(e) => return Err(e),
    };
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    let whole = closed || crate::support::runs_whole_turn(model, &seat, &circle, range, tol)?;
    let (on_cap, on_wall) = if seat.cap_face.is_same(face) {
        (on_face, on_other)
    } else if seat.wall_face.is_same(face) {
        (on_other, on_face)
    } else {
        ogeom_bail!(
            Construction,
            "the face an asymmetric round's first distance runs along is not one of the \
             edge's two faces"
        );
    };
    // Down the wall `on_wall` along its generator, `wall_depth` along the
    // axis; in along the cap `on_cap`, toward the material.
    let wall_depth = on_wall / seat.slope.hypot(1.0);
    let wall_height = -seat.tau * wall_depth;
    let wall_rho = seat.radius_at(wall_height);
    let cap_rho = seat.radius - seat.sigma * seat.tau * on_cap;
    if wall_rho <= tol.confusion() || cap_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "an asymmetric round of {on_cap} on the cap and {on_wall} on the wall reaches \
             the axis of a rim of radius {}",
            seat.radius
        );
    }
    let x = seat.x_ref.vector();
    let points = [
        seat.centre + seat.up * wall_height + x * wall_rho,
        seat.centre + x * seat.radius,
        seat.centre + x * cap_rho,
    ];
    let (a, b) = (points[0] - points[1], points[2] - points[1]);
    let (a, b) = (a / a.magnitude(), b / b.magnitude());
    let w = profile.weight(a.dot(b).clamp(-1.0, 1.0))?;
    Ok(Some(RimConic {
        seat,
        circle,
        whole,
        wall_depth,
        cap_rho,
        points,
        w,
    }))
}

/// Whether `edge` is straight between two planes.
fn straight_between_planes(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    Ok(matches!(curve, Curve::Line(_)) && crate::support::hosts_planar(model, solid, edge, tol)?)
}

/// Whether the faces beside `edge` meet at one angle all along it, so an
/// elliptical round of equal distances is the circular round: a straight
/// edge between planes, or a rim on a cap and its coaxial wall.
fn constant_dihedral(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    if straight_between_planes(model, solid, edge, tol)? {
        return Ok(true);
    }
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    let Curve::Circle(circle) = curve else {
        return Ok(false);
    };
    Ok(crate::support::revolved_seat(model, solid, edge, &circle, tol).is_ok())
}

/// Refuses by name where the rim's contacts stand off their faces, read
/// across from a quarter, the middle and three quarters of the edge.
pub(crate) fn rim_fits(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &RimConic,
    tol: Tolerances,
) -> OgeomResult<()> {
    use ogeom_geom::Curve3d as _;
    let seat = &rim.seat;
    let (curve, (t0, t1)) = crate::support::edge_curve(model, edge, tol)?;
    for fraction in [0.25, 0.5, 0.75] {
        let on_rim = curve.point_at(t0 + (t1 - t0) * fraction, tol)?;
        let radial = {
            let off = on_rim - seat.centre;
            let r = off - seat.up * off.dot(seat.up);
            r / r.magnitude()
        };
        let across = |p: Point| {
            let off = p - seat.centre;
            let height = off.dot(seat.up);
            let r = (off - seat.up * height).magnitude();
            seat.centre + seat.up * height + radial * r
        };
        for (face, p) in [
            (&seat.wall_face, rim.points[0]),
            (&seat.cap_face, rim.points[2]),
        ] {
            if !crate::support::on_face_or_boundary(model, solid, face, across(p), tol)? {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round on a rim sets back past the side of a face beside \
                     it; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                );
            }
        }
    }
    Ok(())
}

/// The largest common scale up to `up_to` of an asymmetric round's two
/// distances at which [`fillet_edge_distances`] on `edge` builds, the
/// round meeting `face` `s * on_face` from the edge and the other face
/// `s * on_other` from it, its section elliptical: every scale probed
/// below it builds and a scale [`LIMIT_TOLERANCE`](crate::LIMIT_TOLERANCE)
/// above it is refused, whatever refuses it. None where no scale leaving
/// both distances longer than the confusion distance builds.
///
/// The seat is read first, without building a band, as
/// [`largest_chamfer_distance`](crate::largest_chamfer_distance) reads a
/// bevel's legs: each contact walks square to the edge, its own distance
/// per unit of scale into its face, and the seat is where the first walk
/// leaves its face, a confusion distance short; on a rim the cap's contact
/// stops short of the axis. On any other seat (a curved crease, an edge
/// along a drum's ruling) each leg walks its face's section in the plane
/// square to the crease, as the band's own legs walk it, at stations over
/// the edge, and the reading also stops where the faces' tangents at a
/// station's contacts no longer meet in front of the crease: the shoulder
/// the conic stands on, gone to infinity where a leg turns a quarter round
/// a drum, and on a concave crease past the other face's contact once the
/// drum's tangent crosses that face's line beyond it. That reading only
/// starts the search: the answer is [`blend_limit`](crate::blend_limit)'s
/// for [`BlendSizing::RoundDistances`](crate::BlendSizing::RoundDistances)
/// with the elliptical section, confirmed by building the round on a copy
/// of the model at scales bracketing it; another section is read there.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to` or a distance is not a usable length, and as
/// [`fillet_edge_distances`] where the seat is refused and no scale builds.
#[allow(
    clippy::too_many_arguments,
    reason = "fillet_edge_distances' own, and the bound"
)]
pub fn largest_fillet_distances(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let sizing = crate::BlendSizing::RoundDistances {
        edge,
        face,
        on_face,
        on_other,
        profile: RoundProfile::Elliptical,
    };
    Ok(crate::blend_limit(model, solid, &sizing, up_to, tol)?.size)
}

/// The seat reading [`largest_fillet_distances`] starts from: the largest
/// common scale up to `up_to` at which each contact of the round stays on
/// its face, read without building a band.
///
/// # Errors
///
/// As [`fillet_edge_distances`] where the seat is refused at any scale.
#[allow(
    clippy::too_many_arguments,
    reason = "fillet_edge_distances' own, and the bound"
)]
pub(crate) fn seat_fillet_distances(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let conf = tol.confusion();
    let profile = RoundProfile::Elliptical;
    let rim = rim_conic(model, solid, edge, face, on_face, on_other, profile, tol)?;
    if rim.is_none() && straight_seat(model, solid, edge, face, tol).is_err() {
        return largest_marched(model, solid, edge, face, [on_face, on_other], up_to, tol);
    }
    let mut largest = crate::fits::legs_largest(
        model,
        solid,
        edge,
        up_to,
        |f| if f.is_same(face) { on_face } else { on_other },
        tol,
    )?;
    if let Some(rim) = &rim {
        // The cap's contact walking in toward the axis reaches it at the
        // rim's radius.
        let seat = &rim.seat;
        if seat.sigma * seat.tau > 0.0 {
            let on_cap = if seat.cap_face.is_same(face) {
                on_face
            } else {
                on_other
            };
            largest = largest.min((seat.radius - conf) / on_cap);
        }
    }
    Ok((largest * on_face.min(on_other) > conf).then_some(largest))
}

/// How many steps [`largest_marched`] walks the scale up in before it
/// bisects the first that does not fit.
const SCALE_STEPS: u32 = 64;

/// How many pieces [`largest_marched`] splits the edge into for its
/// stations.
const STATIONS: u32 = 32;

/// [`seat_fillet_distances`] on a seat neither a straight edge between
/// planes nor a rim: each leg walked along its face's section square to
/// the crease, as the marched band walks it, at stations over the edge,
/// the scale stepped up until a contact first leaves its face or a
/// station's shoulder leaves the front of its crease, and that step
/// bisected.
fn largest_marched(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    [on_face, on_other]: [f64; 2],
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let conf = tol.confusion();
    let mut scratch = model.clone();
    let mut fits = |scale: f64| -> OgeomResult<bool> {
        let spec = crate::Chamfer::Distances {
            face: face.clone(),
            on_face: scale * on_face,
            on_other: scale * on_other,
        };
        let Some(found) = crate::chamfer_marched::crease_contacts(
            &mut scratch,
            solid,
            edge,
            &spec,
            STATIONS,
            None,
            tol,
        )?
        else {
            return Ok(false);
        };
        // Each section's conic stands on the shoulder its faces' tangents
        // meet on in front of the crease: a leg turned a quarter round a
        // drum leaves it no corner.
        if found.shouldered.iter().any(|s| !s) {
            return Ok(false);
        }
        // On a ruled seat the section is read whole, and its conic must not
        // cross a drum it is tangent to, which it does well short of the
        // quarter turn where the drum's leg turns far past the flat's.
        let asked = Asymmetric {
            edge: edge.clone(),
            face: face.clone(),
            on_face: scale * on_face,
            on_other: scale * on_other,
            profile: RoundProfile::Elliptical,
        };
        if crate::engine::ruled_fits(&scratch, solid, &asked, tol)? == Some(false) {
            return Ok(false);
        }
        for (points, f) in found.points.iter().zip(&found.faces) {
            for p in points {
                if !crate::support::on_face_or_boundary(&scratch, solid, f, *p, tol)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    };
    let least = conf / on_face.min(on_other);
    let lo = crate::support::keeping_face_meshes(|| -> OgeomResult<f64> {
        let mut lo = 0.0;
        let mut hi = None;
        for k in 1..=SCALE_STEPS {
            let scale = up_to * f64::from(k) / f64::from(SCALE_STEPS);
            if scale > least && !fits(scale)? {
                hi = Some(scale);
                break;
            }
            lo = scale;
        }
        if let Some(mut hi) = hi {
            for _ in 0..48 {
                let mid = f64::midpoint(lo, hi);
                if fits(mid)? {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            // The contacts are read on their faces within the reach the
            // check reads them by, a hundred confusion distances: back off
            // twice that, so the limiting contact stands on its face.
            lo -= 200.0 * conf / on_face.min(on_other);
        }
        Ok(lo)
    })?;
    Ok((lo * on_face.min(on_other) > conf).then_some(lo))
}

/// The edges of `solid` on `edge`'s circle (the arcs of one rim), `edge`
/// first; `None` where `edge` is not an arc of a circle.
///
/// # Errors
///
/// As reading the solid.
fn rim_arcs(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Shape>>> {
    let circle_of = |e: &Shape| -> OgeomResult<Option<ogeom_math::Circle>> {
        Ok(match ogeom_algo::placed_edge_curve(model, e, tol)? {
            Some((Curve::Circle(c), _)) => Some(c.circle()),
            _ => None,
        })
    };
    let Some(own) = circle_of(edge)? else {
        return Ok(None);
    };
    let mut arcs = vec![edge.clone()];
    for e in ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Edge)? {
        if e.is_same(edge) {
            continue;
        }
        if let Some(c) = circle_of(&e)?
            && c.centre().distance(own.centre()) <= tol.confusion()
            && (c.radius() - own.radius()).abs() <= tol.confusion()
            && c.frame()
                .z()
                .vector()
                .cross(own.frame().z().vector())
                .magnitude()
                <= tol.angular()
        {
            arcs.push(e);
        }
    }
    Ok(Some(arcs))
}

/// The conic round of the rim of `arcs` on `solid` with its seams
/// dissolved: the arcs named on the dissolved solid, rounded by the
/// engine's rims there, the history told of the solid as given; `None`
/// where nothing dissolves or the rims decline.
///
/// # Errors
///
/// As dissolving; a construction failure of the rims declines.
fn conic_rim_dissolved(
    model: &mut Model,
    solid: &Shape,
    arcs: &[Shape],
    section: &crate::engine::RimSection,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut scratch = model.clone();
    let Some(seamless) = crate::engine::dissolve_curved(&mut scratch, solid, tol)? else {
        return Ok(None);
    };
    let asked: Vec<(Shape, f64)> = arcs.iter().map(|e| (e.clone(), 1.0)).collect();
    let Some(rounds) = seamless.rounds(&asked) else {
        return Ok(None);
    };
    let asked: Vec<(Shape, crate::engine::RimSection)> = rounds
        .into_iter()
        .map(|(e, _)| (e, section.clone()))
        .collect();
    match crate::engine::round_rims_with(&mut scratch, &seamless.solid, &asked, tol) {
        Ok(Some(built)) => {
            *model = scratch;
            Ok(Some(Built::new(
                built.shape,
                seamless.carry(&built.history),
            )))
        }
        Ok(None) | Err(ogeom_core::OgeomError::Construction(_)) => Ok(None),
        Err(e) => Err(e),
    }
}
