//! Chamfers: the bevel that replaces an edge.
//!
//! A chamfer along a straight edge between two planar faces is a wedge
//! subtracted, and the wedge's own faces lie *exactly* on the solid's
//! (coplanar, materials aligned), which is the boolean's same-domain case.
//!
//! Three spellings, one construction. The symmetric chamfer cuts the same
//! distance along both faces; the distance-distance form cuts a named
//! distance along a named face and another along its neighbour; the
//! distance-angle form cuts a distance along the named face and leaves it at
//! an angle, with the second distance derived where that bevel meets the
//! other face. All three end in the same wedge subtraction.
//!
//! Two seats per spelling. A straight edge between planes takes the
//! triangular prism above; the circular rim where a coaxial cylinder or cone
//! meets a perpendicular planar cap takes a revolved wedge whose bevel is a
//! *cone*, its leg on a cone's wall running along the slant:
//! the same flanks the rim fillet builds, with the quarter-tube exchanged
//! for the slant, and the same melt taking the legs away.

use crate::support::{
    Blamed, RevolvedSeat, Seat, apply_wedge, edge_curve, planar_face, planar_seat, revolved_flanks,
    revolved_seat,
};
use ogeom_algo::{Built, make_revolution_band};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{ConeSurface, Curve, SurfaceGeometry};
use ogeom_math::Cone;
use ogeom_topo::{Model, Shape};

/// Bevel a straight edge of a solid, cutting `distance` back along each of
/// its two faces.
///
/// The edge must be straight, convex, and shared by exactly two planar faces;
/// the distances are equal (the symmetric chamfer). The result is the boolean
/// difference with a wedge whose legs run along the two faces, so the
/// history reads as a cut: the two faces are modified into their trimmed
/// pieces, the edge's neighbourhood gains the bevel face.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, not convex, not shared by exactly two planar faces of
/// `solid`, or `distance` is not a usable length.
pub fn chamfer_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| chamfer_edge_own(model, solid, &parts[0], distance, tol),
    )
}

/// [`chamfer_edge`] on a shape every node below which is its own.
fn chamfer_edge_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    wedge_for(model, solid, edge, &Chamfer::Symmetric(distance), tol)?
        .apply(model, solid, edge, tol)
}

/// Bevel a straight edge, cutting `on_face` back along `face` and `on_other`
/// along the edge's other face.
///
/// The asymmetric chamfer: `face` names which side the first distance applies
/// to, and must be one of the two faces meeting at the edge.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge`], and additionally if `face` is not one of the edge's
/// two faces.
pub fn chamfer_edge_distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            chamfer_edge_distances_own(model, solid, &parts[0], &parts[1], on_face, on_other, tol)
        },
    )
}

/// [`chamfer_edge_distances`] on a shape every node below which is its own.
fn chamfer_edge_distances_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spec = Chamfer::Distances {
        face: face.clone(),
        on_face,
        on_other,
    };
    wedge_for(model, solid, edge, &spec, tol)?.apply(model, solid, edge, tol)
}

/// Bevel a straight edge, cutting `distance` back along `face` and leaving it
/// at `angle` radians from that face.
///
/// The distance-angle chamfer: the second distance is where the bevel,
/// departing the named face at the given angle, meets the other face. An
/// angle of `π/4` on a square edge reproduces the symmetric chamfer.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge_distances`], and additionally if the bevel at that angle
/// never reaches the other face.
pub fn chamfer_edge_angle(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            chamfer_edge_angle_own(model, solid, &parts[0], &parts[1], distance, angle, tol)
        },
    )
}

/// [`chamfer_edge_angle`] on a shape every node below which is its own.
fn chamfer_edge_angle_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spec = Chamfer::Angle {
        face: face.clone(),
        distance,
        angle,
    };
    wedge_for(model, solid, edge, &spec, tol)?.apply(model, solid, edge, tol)
}

/// The shapes chamfer specs name: each edge, then the face its chamfer
/// runs along where it names one.
fn chamfer_parts(specs: &[(Shape, Chamfer)]) -> Vec<Shape> {
    let mut parts = Vec::with_capacity(specs.len() * 2);
    for (edge, spec) in specs {
        parts.push(edge.clone());
        match spec {
            Chamfer::Symmetric(_) => {}
            Chamfer::Distances { face, .. } | Chamfer::Angle { face, .. } => {
                parts.push(face.clone());
            }
        }
    }
    parts
}

/// The specs again, naming `parts` in the order [`chamfer_parts`] lists
/// them.
fn chamfer_specs(specs: &[(Shape, Chamfer)], parts: &[Shape]) -> Vec<(Shape, Chamfer)> {
    let mut parts = parts.iter().cloned();
    let mut next = |was: &Shape| parts.next().unwrap_or_else(|| was.clone());
    specs
        .iter()
        .map(|(edge, spec)| {
            let edge = next(edge);
            let spec = match spec {
                Chamfer::Symmetric(d) => Chamfer::Symmetric(*d),
                Chamfer::Distances {
                    face,
                    on_face,
                    on_other,
                } => Chamfer::Distances {
                    face: next(face),
                    on_face: *on_face,
                    on_other: *on_other,
                },
                Chamfer::Angle {
                    face,
                    distance,
                    angle,
                } => Chamfer::Angle {
                    face: next(face),
                    distance: *distance,
                    angle: *angle,
                },
            };
            (edge, spec)
        })
        .collect()
}

/// One edge's chamfer, in any of the three spellings.
#[derive(Debug, Clone)]
pub enum Chamfer {
    /// The same distance back along both faces, as [`chamfer_edge`].
    Symmetric(f64),
    /// `on_face` along `face` and `on_other` along the other face, as
    /// [`chamfer_edge_distances`].
    Distances {
        /// The face the first distance runs along.
        face: Shape,
        /// The distance along `face`.
        on_face: f64,
        /// The distance along the edge's other face.
        on_other: f64,
    },
    /// `distance` along `face`, leaving it at `angle` radians, as
    /// [`chamfer_edge_angle`].
    Angle {
        /// The face the distance runs along.
        face: Shape,
        /// The distance along `face`.
        distance: f64,
        /// The bevel's angle from `face`, in radians.
        angle: f64,
    },
}

/// Bevel several edges of a solid as one operation, each by `distance`.
///
/// As [`chamfer_edges_with`] with the symmetric chamfer on every edge.
///
/// # Errors
///
/// As [`chamfer_edges_with`].
pub fn chamfer_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let specs: Vec<(Shape, Chamfer)> = edges
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(distance)))
        .collect();
    chamfer_edges_with(model, solid, &specs, tol)
}

/// Bevel several edges of a solid as one operation, each with its own
/// chamfer.
///
/// Every wedge is built on the solid as it stands before the call, and then
/// every wedge is applied. That is what makes the bevels mitre: where two
/// meet at a vertex, each wedge still reaches the corner, the two cut the
/// region above either bevel plane, and the planes meet along their own
/// line with neither a step nor a cap. One edge at a time
/// ([`chamfer_edge`] in a loop) builds the second wedge on an edge the
/// first bevel has already shortened, and each corner keeps a small
/// tetrahedron and two extra faces. Where three bevels meet at a convex
/// vertex the three planes meet at one point, the mitre of three planes.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge`], [`chamfer_edge_distances`] and
/// [`chamfer_edge_angle`] per edge, judged on the solid as it stands before
/// the call; and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if no edges are given.
pub fn chamfer_edges_with(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    chamfer_edges_blamed(model, solid, specs, tol).map_err(|blamed| blamed.error)
}

/// [`chamfer_edges_with`], its refusal naming the specs it concerns.
pub(crate) fn chamfer_edges_blamed(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, Blamed> {
    let mut blame = Vec::new();
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &chamfer_parts(specs),
        |model, solid, parts| {
            chamfer_edges_with_own(model, solid, &chamfer_specs(specs, parts), tol).map_err(
                |blamed| {
                    blame = blamed.edges;
                    blamed.error
                },
            )
        },
    )
    .map_err(|error| Blamed {
        edges: blame,
        error,
    })
}

/// [`chamfer_edges_with`] on a shape every node below which is its own.
///
/// A wedge refused names its spec; one that does not apply names its
/// spec, then the specs before it whose edges share a vertex with its
/// edge, whose bevels it was applied after.
fn chamfer_edges_with_own(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, Blamed> {
    if specs.is_empty() {
        return Err(
            ogeom_core::ogeom_err!(Construction, "a chain of no edges bevels nothing").into(),
        );
    }
    let mut wedges = Vec::with_capacity(specs.len());
    for (i, (edge, spec)) in specs.iter().enumerate() {
        wedges.push(
            wedge_for(model, solid, edge, spec, tol).map_err(|error| Blamed {
                edges: vec![i],
                error,
            })?,
        );
    }
    let mut built = Built::from_nothing(solid.clone());
    for (i, ((edge, _), wedge)) in specs.iter().zip(wedges).enumerate() {
        let step = match wedge.apply(model, &built.shape, edge, tol) {
            Ok(step) => step,
            Err(error) => {
                let mut edges = vec![i];
                let ends = |e: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
                    Ok(ogeom_algo::edge_vertices(model, e)?
                        .map(|(a, b)| vec![a.node(), b.node()])
                        .unwrap_or_default())
                };
                let own = ends(edge)?;
                for (j, (other, _)) in specs[..i].iter().enumerate() {
                    if ends(other)?.iter().any(|v| own.contains(v)) {
                        edges.push(j);
                    }
                }
                return Err(Blamed { edges, error });
            }
        };
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// The seat of a rim's bevel: a revolved seat, its wall a cylinder or a
/// cone.
fn rim_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &ogeom_geom::CircleCurve,
    tol: Tolerances,
) -> OgeomResult<RevolvedSeat> {
    revolved_seat(model, solid, edge, rim, tol)
}

/// How far along the axis a leg of `on_wall` down the wall's ruling
/// reaches: the leg itself on a cylinder, its run along the slant on a
/// cone.
fn wall_depth(seat: &RevolvedSeat, on_wall: f64) -> f64 {
    if seat.slope == 0.0 {
        on_wall
    } else {
        on_wall / seat.slope.hypot(1.0)
    }
}

/// The angle at the rim between the cap's leg and the wall's, in the
/// rim's section: square on a cylinder, and on a cone the cone's lean
/// from it, obtuse where the wall leans out under the cap.
fn rim_corner(seat: &RevolvedSeat) -> f64 {
    // In the section, the cap's leg runs `-sigma tau` radially and the
    // wall's runs `(-slope tau, -tau)`, both from the rim.
    (seat.sigma * seat.slope / seat.slope.hypot(1.0)).acos()
}

/// The wedge one chamfer cuts, built on `solid` as it stands.
fn wedge_for(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    match legs(model, solid, edge, spec, tol)? {
        Legs::Straight(seat, distances) => bevel(model, &seat, distances, tol),
        Legs::Rim(seat, on_wall, on_cap) => {
            revolved_bevel(model, edge, &seat, on_wall, on_cap, tol)
        }
    }
}

/// The largest measure up to `up_to` at which a chamfer of `spec`'s
/// spelling on `edge` of `solid` is not refused for running past its
/// faces: for [`Chamfer::Symmetric`] the distance, for
/// [`Chamfer::Distances`] the largest common scale of its two distances,
/// and for [`Chamfer::Angle`] the distance along its face at its angle;
/// the measure `spec` itself names is not read. None where no measure
/// leaving both legs longer than the confusion distance fits.
///
/// Read from the seat, without building a wedge. On a straight edge
/// between planes, a convex bevel's leg runs at most as far as its face
/// reaches back from the edge, the bound [`chamfer_edges_with`] refuses
/// past, and the answer is exact against it; a concave one fuses and runs
/// past no face. On a rim, the cap's leg stops short of the axis, and
/// each leg's contact stays on its face at the stations
/// [`largest_fillet_radius`](crate::largest_fillet_radius) reads, a
/// confusion distance short.
///
/// # Errors
///
/// As [`chamfer_edges_with`] where the seat or spelling is refused at any
/// measure, and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if `up_to` is not a usable length; as
/// [`largest_fillet_radius`](crate::largest_fillet_radius) where a
/// contact's walk on a rim depends on the measure.
pub fn largest_chamfer_distance(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let conf = tol.confusion();
    if !up_to.is_finite() || up_to <= conf {
        ogeom_bail!(
            Construction,
            "a largest distance up to {up_to} asks for nothing"
        );
    }
    // The legs at a measure of one, which scale with it.
    let unit = match spec {
        Chamfer::Symmetric(_) => Chamfer::Symmetric(1.0),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => Chamfer::Distances {
            face: face.clone(),
            on_face: *on_face,
            on_other: *on_other,
        },
        Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
            face: face.clone(),
            distance: 1.0,
            angle: *angle,
        },
    };
    let mut largest = up_to;
    let shortest = match legs(model, solid, edge, &unit, tol)? {
        Legs::Straight(seat, per) => {
            if seat.convex {
                for (i, leg) in per.into_iter().enumerate() {
                    if leg > 0.0 {
                        let along = seat.leg(i, tol)?;
                        let reach = crate::support::face_reach(
                            model,
                            &seat.faces[i],
                            seat.start,
                            along,
                            tol,
                        )?;
                        largest = largest.min(reach / leg);
                    }
                }
            }
            per[0].min(per[1])
        }
        Legs::Rim(seat, on_wall, on_cap) => {
            // The cap's leg ends at the radius `radius - sigma tau on_cap`,
            // which stays above the confusion distance.
            let inward = seat.sigma * seat.tau * on_cap;
            if inward > 0.0 {
                largest = largest.min((seat.radius - 2.0 * conf) / inward);
            }
            // A cone's wall leg ends at the radius the wall has a depth
            // down, which narrows toward the apex with the measure.
            let narrowing = seat.slope * seat.tau * wall_depth(&seat, on_wall);
            if narrowing > 0.0 {
                largest = largest.min((seat.radius - 2.0 * conf) / narrowing);
            }
            let per = |face: &Shape| {
                if face.is_same(&seat.wall_face) {
                    on_wall
                } else if face.is_same(&seat.cap_face) {
                    on_cap
                } else {
                    0.0
                }
            };
            largest = largest.min(crate::fits::legs_largest(
                model, solid, edge, largest, per, tol,
            )?);
            on_wall.min(on_cap)
        }
    };
    Ok((largest * shortest > conf).then_some(largest))
}

/// The seat a chamfer stands on and how far its legs run along each face.
pub(crate) enum Legs {
    /// A straight edge between planes, its legs along the seat's faces in
    /// order.
    Straight(Seat, [f64; 2]),
    /// A rim on a cylinder or a cone: its legs down the wall and in along the
    /// cap.
    Rim(RevolvedSeat, f64, f64),
}

/// Read one chamfer's seat and legs from `solid` as it stands.
pub(crate) fn legs(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Legs> {
    match spec {
        Chamfer::Symmetric(distance) => match seat_kind(model, edge, tol)? {
            SeatKind::Straight => {
                let seat = planar_seat(model, solid, edge, tol)?;
                Ok(Legs::Straight(seat, [*distance, *distance]))
            }
            SeatKind::Rim(rim) => {
                let seat = rim_seat(model, solid, edge, &rim, tol)?;
                Ok(Legs::Rim(seat, *distance, *distance))
            }
        },
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match seat_kind(model, edge, tol)? {
            SeatKind::Straight => {
                let seat = planar_seat(model, solid, edge, tol)?;
                let i = seat_side(&seat, face)?;
                let mut distances = [0.0; 2];
                distances[i] = *on_face;
                distances[1 - i] = *on_other;
                Ok(Legs::Straight(seat, distances))
            }
            SeatKind::Rim(rim) => {
                let seat = rim_seat(model, solid, edge, &rim, tol)?;
                let (on_wall, on_cap) = revolved_side(&seat, face, *on_face, *on_other)?;
                Ok(Legs::Rim(seat, on_wall, on_cap))
            }
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            let (distance, angle) = (*distance, *angle);
            if !angle.is_finite() || angle <= tol.angular() {
                ogeom_bail!(
                    Construction,
                    "a chamfer at an angle of {angle} cuts nothing"
                );
            }
            match seat_kind(model, edge, tol)? {
                SeatKind::Straight => {
                    let seat = planar_seat(model, solid, edge, tol)?;
                    let i = seat_side(&seat, face)?;
                    // In the cross-section: from the contact on the named
                    // face, the bevel leaves at `angle` into the wedge's own
                    // side: the material on a convex edge, the open dihedral
                    // on a concave one. Where it crosses the other leg's ray
                    // is the derived distance; no crossing, no chamfer.
                    let sign = if seat.convex { 1.0 } else { -1.0 };
                    let a = seat.leg(i, tol)? * sign;
                    let b = seat.leg(1 - i, tol)? * sign;
                    let inward = -seat.normals[i] * sign;
                    let denominator = angle.sin().mul_add(b.dot(a), angle.cos() * b.dot(inward));
                    if denominator <= tol.angular() {
                        ogeom_bail!(
                            Construction,
                            "the bevel at that angle never meets the edge's other face"
                        );
                    }
                    let derived = distance * angle.sin() / denominator;
                    let mut distances = [0.0; 2];
                    distances[i] = distance;
                    distances[1 - i] = derived;
                    Ok(Legs::Straight(seat, distances))
                }
                SeatKind::Rim(rim) => {
                    // In the rim's section the legs make the triangle with
                    // the bevel: `angle` at the named face's contact and the
                    // corner's own angle at the rim. A cylinder's seat is
                    // square, the derived distance the plain tangent; past
                    // the two angles' sum of a half turn the bevel walks
                    // away from the other face instead of toward it.
                    let seat = rim_seat(model, solid, edge, &rim, tol)?;
                    let corner = rim_corner(&seat);
                    let past = if seat.slope == 0.0 {
                        angle >= core::f64::consts::FRAC_PI_2 - tol.angular()
                    } else {
                        angle + corner >= core::f64::consts::PI - tol.angular()
                    };
                    if past {
                        ogeom_bail!(
                            Construction,
                            "the bevel at that angle never meets the edge's other face"
                        );
                    }
                    let derived = if seat.slope == 0.0 {
                        distance * angle.tan()
                    } else {
                        distance * angle.sin() / (angle + corner).sin()
                    };
                    let (on_wall, on_cap) = revolved_side(&seat, face, distance, derived)?;
                    Ok(Legs::Rim(seat, on_wall, on_cap))
                }
            }
        }
    }
}

/// A chamfer's wedge, built and not yet applied: its faces, and whether it
/// fuses (a concave edge) or cuts (a convex one).
struct Wedge {
    faces: Vec<Shape>,
    /// The wedge already closed into a solid, where it was built as one.
    solid: Option<Shape>,
    additive: bool,
}

impl Wedge {
    fn apply(
        self,
        model: &mut Model,
        solid: &Shape,
        edge: &Shape,
        tol: Tolerances,
    ) -> OgeomResult<Built> {
        let mut built = match self.solid {
            Some(wedge) => {
                let mut result = if self.additive {
                    ogeom_bool::fuse(model, solid, &wedge, tol)?
                } else {
                    ogeom_bool::cut(model, solid, &wedge, tol)?
                };
                result.history.delete(edge);
                result
            }
            None => apply_wedge(model, solid, Some(edge), &self.faces, self.additive, tol)?,
        };
        crate::support::credit_new_faces(model, solid, edge, &mut built)?;
        Ok(built)
    }
}

/// Which seat a chamfer is standing on, read from the edge's curve.
enum SeatKind {
    /// A straight edge between planes: the triangular-prism wedge.
    Straight,
    /// A circular rim: the revolved wedge with a conical bevel.
    Rim(ogeom_geom::CircleCurve),
}

fn seat_kind(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<SeatKind> {
    let (curve, _) = edge_curve(model, edge, tol)?;
    match curve {
        Curve::Line(_) => Ok(SeatKind::Straight),
        Curve::Circle(c) => Ok(SeatKind::Rim(c)),
        _ => ogeom_bail!(
            Construction,
            "chamfering an edge that is neither straight nor circular needs \
             the marching blend machinery"
        ),
    }
}

/// Assign a named face's distance to the wall or the cap.
///
/// On the wall the distance runs axially down from the rim; on the cap it
/// runs radially in from it.
fn revolved_side(
    seat: &RevolvedSeat,
    face: &Shape,
    on_face: f64,
    on_other: f64,
) -> OgeomResult<(f64, f64)> {
    if seat.wall_face.is_same(face) {
        Ok((on_face, on_other))
    } else if seat.cap_face.is_same(face) {
        Ok((on_other, on_face))
    } else {
        ogeom_bail!(
            Construction,
            "the named face does not meet the edge being chamfered"
        )
    }
}

/// The revolved wedge with the quarter-tube exchanged for a slant: legs
/// `on_wall` down the wall along its ruling and `on_cap` radially along the
/// cap, and the cone between the two tangency rings as the bevel.
fn revolved_bevel(
    model: &mut Model,
    edge: &Shape,
    seat: &RevolvedSeat,
    on_wall: f64,
    on_cap: f64,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    for distance in [on_wall, on_cap] {
        if !distance.is_finite() || distance <= tol.confusion() {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    let cap_rho = seat.sigma.mul_add(-(seat.tau * on_cap), seat.radius);
    if cap_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_cap} along the cap swallows the axis of a rim \
             of radius {}",
            seat.radius
        );
    }
    // A cone's wall leg runs its distance along the slant, a depth down
    // the axis, to the wall's radius there.
    let depth = wall_depth(seat, on_wall);
    let wall_rho = seat.radius_at(-seat.tau * depth);
    if wall_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_wall} down the wall runs past the apex of its cone"
        );
    }
    if let Some((start, sweep)) = arc_span(model, edge, seat, tol)? {
        return sector_bevel(model, seat, depth, cap_rho, start, sweep, tol);
    }
    let flanks = revolved_flanks(model, seat, depth, cap_rho, tol)?;

    // The bevel: the cone through both tangency rings (reference radius
    // `cap_rho` at the cap's level, the wall's radius a depth below), a
    // cylinder where the two radii agree. Unlike the fillet's quarter-tube,
    // whose away-from-the-tube normal tracks the wedge seat by seat, the
    // cone's natural normal always points away from the axis, and the
    // wedge sits on the axis side of the slant exactly when `sigma` and
    // `tau` agree.
    let bevel_band = {
        let slope = (cap_rho - wall_rho) / (seat.tau * depth);
        // The domain covers the band's two rows (the cap ring at zero and
        // the wall ring a depth away) with a margin that stays clear of the
        // apex, where the surface degenerates.
        let rows = (
            0.0_f64.min(-seat.tau * depth),
            0.0_f64.max(-seat.tau * depth),
        );
        let frame = seat.frame_at(seat.centre, tol)?;
        let surface: SurfaceGeometry = if slope.atan().abs() <= tol.angular() {
            let pad = 0.1 * depth;
            ogeom_geom::CylinderSurface::new(
                ogeom_math::Cylinder::new(frame, cap_rho, tol)?,
                (rows.0 - pad, rows.1 + pad),
            )?
            .into()
        } else {
            // Beside a cone's wall the bevel's apex can come near its
            // rows; the margin then stops halfway to it.
            let pad = if seat.slope == 0.0 {
                0.1 * depth
            } else {
                let apex = -cap_rho / slope;
                let clear = (apex - rows.0).abs().min((apex - rows.1).abs());
                (0.1 * depth).min(0.5 * clear)
            };
            let cone = Cone::new(frame, cap_rho, slope.atan(), tol)?;
            ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into()
        };
        let band = make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?;
        if seat.sigma * seat.tau > 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    Ok(Wedge {
        faces: vec![flanks.wall_band, flanks.annulus, bevel_band],
        solid: None,
        additive: seat.additive(),
    })
}

/// Where a rim that is only an arc of its circle starts, and how far it
/// turns, as angles about the seat's axis from its `x_ref`, turning
/// positively about `up`. `None` for a whole circle, or an arc of one whose
/// wall runs on round the rest of the turn.
fn arc_span(
    model: &Model,
    edge: &Shape,
    seat: &RevolvedSeat,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    use ogeom_geom::Curve3d as _;
    let (curve, range) = edge_curve(model, edge, tol)?;
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    let Curve::Circle(arc) = &curve else {
        return Ok(None);
    };
    if closed || crate::support::runs_whole_turn(model, seat, arc, range, tol)? {
        return Ok(None);
    }
    let frame = seat.frame_at(seat.centre, tol)?;
    let angle = |t: f64| -> OgeomResult<f64> {
        let d = curve.point_at(t, tol)? - seat.centre;
        Ok(d.dot(frame.y().vector()).atan2(d.dot(frame.x().vector())))
    };
    let turn = core::f64::consts::TAU;
    let (a, b, m) = (
        angle(range.0)?,
        angle(range.1)?,
        angle(0.5 * (range.0 + range.1))?,
    );
    // The arc runs from `a` to `b` one way round or the other: the way
    // that passes through its middle.
    let forward = (b - a).rem_euclid(turn);
    let start = if (m - a).rem_euclid(turn) <= forward {
        a
    } else {
        b
    };
    let sweep = if start == a { forward } else { turn - forward };
    Ok(Some((start, sweep)))
}

/// The wedge of a chamfer along an arc of a rim: the bevel's cross-section,
/// the triangle between the rim and the two tangency rings, turned about the
/// rim's axis through the arc alone. Its ends are the triangles a straight
/// bevel meeting the arc at a tangent ends in, so the bevels of a chain of
/// lines and arcs abut exactly; a whole ring would cut the material on the
/// far side of the axis as well.
fn sector_bevel(
    model: &mut Model,
    seat: &RevolvedSeat,
    depth: f64,
    cap_rho: f64,
    start: f64,
    sweep: f64,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    let frame = seat.frame_at(seat.centre, tol)?;
    let radial = frame.x().vector() * start.cos() + frame.y().vector() * start.sin();
    let rim = seat.centre + radial * seat.radius;
    let wall = if seat.slope == 0.0 {
        rim - seat.up * (seat.tau * depth)
    } else {
        let height = -seat.tau * depth;
        seat.centre + seat.up * height + radial * seat.radius_at(height)
    };
    let cap = seat.centre + radial * cap_rho;
    let across = seat.up.cross(radial);
    let profile = crate::support::planar_face(model, &[rim, wall, cap], across, tol)?;
    let axis = ogeom_math::Axis::new(seat.centre, frame.z());
    let solid = ogeom_algo::make_revolution(model, &profile, axis, sweep, tol)?.shape;
    Ok(Wedge {
        faces: Vec::new(),
        solid: Some(solid),
        additive: seat.additive(),
    })
}

/// Which side of the seat a named face is, by identity.
fn seat_side(seat: &Seat, face: &Shape) -> OgeomResult<usize> {
    if seat.faces[0].is_same(face) {
        Ok(0)
    } else if seat.faces[1].is_same(face) {
        Ok(1)
    } else {
        ogeom_bail!(
            Construction,
            "the named face does not meet the edge being chamfered"
        )
    }
}

/// The one construction under all three spellings: the wedge with legs
/// `distances[i]` along face `i`, subtracted.
fn bevel(
    model: &mut Model,
    seat: &Seat,
    distances: [f64; 2],
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    for distance in distances {
        if !distance.is_finite() || distance <= tol.confusion() {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    // On a concave edge every leg mirrors: the wedge sits in the open
    // dihedral, its legs walk the faces' planes into it, and its strips face
    // the material they will melt against with *opposed* orientation, which
    // is exactly what a fuse cancels.
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(0, tol)? * sign;
    let b = seat.leg(1, tol)? * sign;
    // A setback runs back across each face; past the face's far side it
    // would cut through the face and on into whatever lies beyond it.
    if seat.convex {
        for (i, leg) in [a, b].into_iter().enumerate() {
            let reach = crate::support::face_reach(model, &seat.faces[i], seat.start, leg, tol)?;
            if distances[i] > reach + tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "a chamfer of {} on the edge from {:?} to {:?} runs past its face, \
                     which reaches {reach} back from the edge",
                    distances[i],
                    seat.start,
                    seat.end
                );
            }
        }
    }

    let travel = seat.end - seat.start;
    let apex0 = seat.start;
    let apex1 = seat.end;
    let a0 = apex0 + a * distances[0];
    let b0 = apex0 + b * distances[1];
    let a1 = a0 + travel;
    let b1 = b0 + travel;

    // The bevel's outward normal: perpendicular to the cut line and the edge,
    // pointing from the apex toward the cut. For equal distances this is the
    // leg bisector exactly.
    let bevel_out = {
        let across = b0 - a0;
        let mut n = seat.along.cross(across);
        let m = n.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "the chamfer's cut line has no direction");
        }
        n /= m;
        if n.dot(a0 - apex0) < 0.0 {
            n = -n;
        }
        n
    };

    // The wedge: a triangular prism whose apex line is the edge and whose
    // legs run the distances along each face. Built from five explicit planar
    // faces rather than swept, because a sweep's walls are extrusion
    // surfaces even when they are geometrically planes, and the boolean's
    // same-domain resolution (which is what makes the coplanar legs melt
    // into the solid's own faces) recognises coincidence between *planes*.
    let faces = [
        planar_face(model, &[apex0, a0, b0], -seat.along, tol)?,
        planar_face(model, &[apex1, a1, b1], seat.along, tol)?,
        planar_face(model, &[apex0, a0, a1, apex1], seat.normals[0] * sign, tol)?,
        planar_face(model, &[apex0, b0, b1, apex1], seat.normals[1] * sign, tol)?,
        planar_face(model, &[a0, b0, b1, a1], bevel_out, tol)?,
    ];
    Ok(Wedge {
        faces: faces.to_vec(),
        solid: None,
        additive: !seat.convex,
    })
}
