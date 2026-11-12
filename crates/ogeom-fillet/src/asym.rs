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

use crate::support::{Seat, face_from_edges, planar_face, segment_between};
use ogeom_algo::{Built, make_edge_between, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSplineCurve, Curve, ExtrusionSurface, PlaneSurface};
use ogeom_math::bspline::Weighted;
use ogeom_math::{Direction, KnotVector, Plane, Point, Vector};
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
    fn weight(self, cos_alpha: f64) -> OgeomResult<f64> {
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

/// What a straight seat's band sweeps along its crease.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Section {
    /// The rolling ball's arc, of this radius.
    Ball(f64),
    /// A conic arc tangent to each face of the seat at its distance from
    /// the crease, in the seat's face order, of weight `w` on the crease.
    Conic {
        /// The distance along each face of the seat.
        d: [f64; 2],
        /// The weight of the crease's control point.
        w: f64,
    },
}

impl Section {
    /// The section's reach from the crease: the ball's radius, or the
    /// farther distance, which scales the steps an end's run-out takes.
    pub(crate) fn size(self) -> f64 {
        match self {
            Self::Ball(radius) => radius,
            Self::Conic { d, .. } => d[0].max(d[1]),
        }
    }

    /// From a point of the seat's crease to a point well inside the band's
    /// section: the ball's centre, or the conic arc's middle. `None` where
    /// the faces meet too sharply or too flat to seat a ball.
    pub(crate) fn centre(self, seat: &Seat, tol: Tolerances) -> OgeomResult<Option<Vector>> {
        match self {
            Self::Ball(radius) => crate::fillet::seat_centre(seat, radius, tol),
            Self::Conic { d, w } => {
                let f = legs(seat, tol)?;
                Ok(Some((f[0] * d[0] + f[1] * d[1]) / (2.0 * (1.0 + w))))
            }
        }
    }

    /// Points of the band's section square to the crease at `at`, with
    /// `inward` from [`Section::centre`]: the crease and the two contacts
    /// first, then each leg's middle and the arc's.
    pub(crate) fn points(
        self,
        seat: &Seat,
        at: Point,
        inward: Vector,
        tol: Tolerances,
    ) -> OgeomResult<[Point; 6]> {
        match self {
            Self::Ball(radius) => Ok(crate::fillet::section_at(seat, at, inward, radius)),
            Self::Conic { d, .. } => {
                let f = legs(seat, tol)?;
                let touch = [at + f[0] * d[0], at + f[1] * d[1]];
                Ok([
                    at,
                    touch[0],
                    touch[1],
                    at.midpoint(touch[0]),
                    at.midpoint(touch[1]),
                    at + inward,
                ])
            }
        }
    }
}

/// On each face of the seat, the unit direction square to the crease that
/// runs along the face away from it.
fn legs(seat: &Seat, tol: Tolerances) -> OgeomResult<[Vector; 2]> {
    let sign = if seat.convex { 1.0 } else { -1.0 };
    Ok([seat.leg(0, tol)? * sign, seat.leg(1, tol)? * sign])
}

/// The conic band of a straight seat: the wedge between the faces' legs
/// and the swept arc, cut from the solid on a convex crease and fused on a
/// concave one, trimmed first by each of `trims` as the circular round's.
pub(crate) fn seated_conic(
    model: &mut Model,
    solid: &Shape,
    seat: &Seat,
    section: Section,
    edge: Option<&Shape>,
    trims: &[(Point, Vector)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let Section::Conic { d, w } = section else {
        ogeom_bail!(Construction, "a conic band asked of a ball's section");
    };
    // The faces in the order that turns positively about the crease, as
    // the circular round orders them, so the arc and the extrusion's chart
    // run from the first face's contact to the second's.
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
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2));
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let f = legs(seat, tol)?;
    let (fa, fb) = (f[first], f[second]);
    let (da, db) = (d[first], d[second]);

    let travel = seat.end - seat.start;
    let length = travel.magnitude();
    let apex0 = seat.start;
    let apex1 = seat.end;
    let contact_a0 = apex0 + fa * da;
    let contact_b0 = apex0 + fb * db;
    let contact_a1 = contact_a0 + travel;
    let contact_b1 = contact_b0 + travel;
    let along_dir = Direction::new(seat.along, tol)?;

    // The section's arc at height `h` along the crease, from the first
    // face's contact to the second's.
    let arc_curve = |h: f64| -> OgeomResult<Curve> {
        let lift = seat.along * h;
        conic_arc([contact_a0 + lift, apex0 + lift, contact_b0 + lift], w, tol)
    };
    let arc = |model: &mut Model, h: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        Ok(make_edge_between(model, arc_curve(h)?, (0.0, 1.0), from, to, tol)?.shape)
    };

    // The legs, in the solid's own faces' planes.
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

    let reach = da.max(db);
    let cap = |model: &mut Model,
               apex: Point,
               ca: Point,
               cb: Point,
               h: f64,
               outward: Vector|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex, Direction::new(outward, tol)?);
        let span = (apex.distance(ca).max(apex.distance(cb)) + reach).max(1.0) * 2.0;
        let surface = PlaneSurface::over(plane, (-span, span), (-span, span))?;
        let apex_v = make_vertex(model, apex).shape;
        let ca_v = make_vertex(model, ca).shape;
        let cb_v = make_vertex(model, cb).shape;
        let edges = vec![
            segment_between(model, (&apex_v, apex), (&ca_v, ca), tol)?,
            arc(model, h, &ca_v, &cb_v)?,
            segment_between(model, (&cb_v, cb), (&apex_v, apex), tol)?,
        ];
        let edges: Vec<Shape> = if (ca - apex).cross(cb - apex).dot(outward) < 0.0 {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        };
        face_from_edges(model, surface.into(), &edges, tol)
    };
    let cap0 = cap(model, apex0, contact_a0, contact_b0, 0.0, -seat.along)?;
    let cap1 = cap(model, apex1, contact_a1, contact_b1, length, seat.along)?;

    // The band: the arc swept along the crease. Its chart's normal, the
    // arc's tangent across the sweep, faces the crease as a cylinder's
    // faces away from its axis, so the face enters the wedge reversed as
    // the circular round's does.
    let band = {
        let surface = ExtrusionSurface::new(arc_curve(0.0)?, along_dir, length)?;
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

    let faces = [leg_a, leg_b, cap0, cap1, band.reversed()];
    let mut wedge = crate::support::wedge_solid(model, &faces, false, tol)?;
    let span = length + (reach * 2.0) * 4.0;
    for &plane in trims {
        wedge = crate::support::trimmed_wedge(model, &wedge, plane, span, tol)?;
    }
    crate::support::apply_wedge_solid(model, solid, edge, &wedge, !seat.convex, tol)
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
fn conic_arc(p: [Point; 3], w: f64, tol: Tolerances) -> OgeomResult<Curve> {
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
        control.push(Weighted::new(*p0, 1.0, tol)?);
        control.push(Weighted::new(*a, weight, tol)?);
        if i + 1 < n {
            #[allow(clippy::cast_precision_loss, reason = "at most 256 pieces")]
            let k = (i + 1) as f64 / n as f64;
            knots.extend([k, k]);
        }
    }
    control.push(Weighted::new(pieces[n - 1][2], 1.0, tol)?);
    knots.extend([1.0; 3]);
    Ok(Curve::BSpline(BSplineCurve::rational(
        KnotVector::new(knots, 2)?,
        control,
    )?))
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
/// edge is not straight between two planes, or a contact sets back past
/// its face's side; with equal distances and the elliptical profile, as
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
    if profile == RoundProfile::Elliptical && on_face == on_other {
        // The circular round, by the same call.
        let radius = equal_round(model, solid, edge, face, on_face, tol)?;
        return crate::fillet_edge(model, solid, edge, radius, tol);
    }
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            fillet_edge_distances_own(
                model, solid, &parts[0], &parts[1], on_face, on_other, profile, tol,
            )
        },
    )
}

/// [`fillet_edge_distances`] on a shape every node below which is its own.
#[allow(clippy::too_many_arguments, reason = "fillet_edge_distances' own")]
fn fillet_edge_distances_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let Some(section) =
        straight_section(model, solid, edge, face, on_face, on_other, profile, tol)?
    else {
        ogeom_bail!(
            Construction,
            "the circular round asked as an asymmetric one"
        );
    };
    contacts_fit(model, solid, edge, section, tol)?;
    let mut built = crate::fillet::planar_fillet(model, solid, edge, section, None, tol)?;
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(built)
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
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    if !matches!(curve, Curve::Line(_)) || !crate::support::hosts_planar(model, solid, edge, tol)? {
        ogeom_bail!(
            Construction,
            "an asymmetric round is built on a straight edge between two planes; this \
             edge is curved or a face beside it is not planar; see docs/CAPABILITIES.md, \
             fillet.asymmetric-rounds"
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
    if profile == RoundProfile::Elliptical && on_face == on_other {
        return Ok(None);
    }
    let mut d = [0.0; 2];
    d[named] = on_face;
    d[1 - named] = on_other;
    let f = legs(&seat, tol)?;
    let w = profile.weight(f[0].dot(f[1]).clamp(-1.0, 1.0))?;
    Ok(Some(Section::Conic { d, w }))
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
fn contacts_fit(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    section: Section,
    tol: Tolerances,
) -> OgeomResult<()> {
    let Section::Conic { d, .. } = section else {
        return Ok(());
    };
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

/// The request's rounds and bevels of one radius each, built by
/// `symmetric`, then each asymmetric round in the order of its edge's node
/// on what they leave, its edge and face traced through them.
///
/// An asymmetric round meeting another edge of the request at a vertex is
/// refused by name: the corner where it meets another band is not closed.
pub(crate) fn after_symmetric(
    model: &mut Model,
    others: &[Shape],
    mut asymmetric: Vec<Asymmetric>,
    symmetric: impl FnOnce(&mut Model) -> OgeomResult<Built>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    asymmetric.sort_by_key(|a| a.edge.node());
    let ends = |model: &Model, edge: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
        Ok(match ogeom_algo::edge_vertices(model, edge)? {
            Some((a, b)) => vec![a.node(), b.node()],
            None => Vec::new(),
        })
    };
    for (i, a) in asymmetric.iter().enumerate() {
        let own = ends(model, &a.edge)?;
        let other_edges = others.iter().chain(
            asymmetric
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, b)| &b.edge),
        );
        for other in other_edges {
            if ends(model, other)?.iter().any(|v| own.contains(v)) {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round meets another blended edge at a vertex, and the \
                     corner where they meet is not closed yet; see docs/CAPABILITIES.md, \
                     fillet.asymmetric-rounds"
                );
            }
        }
    }
    let mut built = symmetric(model)?;
    for a in &asymmetric {
        let [edge] = built.history.trace(&a.edge) else {
            ogeom_bail!(
                Construction,
                "an asymmetric round's edge does not come through the request's other \
                 blends as one edge"
            );
        };
        let edge = edge.clone();
        // The named face's image beside the edge, on what the request's
        // other blends left.
        let beside = crate::support::planar_seat(model, &built.shape, &edge, tol)?.faces;
        let images = built.history.trace(&a.face);
        let Some(face) = beside
            .iter()
            .find(|f| images.iter().any(|g| g.is_same(f)))
            .cloned()
        else {
            ogeom_bail!(
                Construction,
                "the face an asymmetric round's first distance runs along is not one of the \
                 edge's two faces"
            );
        };
        let next = fillet_edge_distances(
            model,
            &built.shape,
            &edge,
            &face,
            a.on_face,
            a.on_other,
            a.profile,
            tol,
        )?;
        built = Built::new(next.shape.clone(), built.history.then(&next.history));
    }
    Ok(built)
}
