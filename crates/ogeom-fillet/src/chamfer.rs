//! Chamfers: the bevel that replaces an edge.
//!
//! P4's opening stone, built deliberately on M3's shoulders: a chamfer along
//! a straight edge between two planar faces is a wedge subtracted, and the
//! wedge's own faces lie *exactly* on the solid's — coplanar, materials
//! aligned — which is the same-domain case the boolean learned to resolve.
//!
//! Three spellings, one construction. The symmetric chamfer cuts the same
//! distance along both faces; the distance-distance form cuts a named
//! distance along a named face and another along its neighbour; the
//! distance-angle form cuts a distance along the named face and leaves it at
//! an angle, with the second distance derived where that bevel meets the
//! other face. All three end in the same wedge subtraction.
//!
//! Two seats per spelling. A straight edge between planes takes the
//! triangular prism above; the circular rim where a cylindrical wall meets a
//! perpendicular planar cap takes a revolved wedge whose bevel is a *cone* —
//! the same flanks the rim fillet builds, with the quarter-tube exchanged
//! for the slant, and the same melt taking the legs away.

use crate::support::{
    RevolvedSeat, Seat, apply_wedge, edge_curve, planar_face, planar_seat, revolved_flanks,
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
/// difference with a wedge whose legs run along the two faces — so the
/// history reads as a cut: the two faces are modified into their trimmed
/// pieces, the edge's neighbourhood gains the bevel face.
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
    wedge_for(model, solid, edge, &Chamfer::Symmetric(distance), tol)?
        .apply(model, solid, edge, tol)
}

/// Bevel a straight edge, cutting `on_face` back along `face` and `on_other`
/// along the edge's other face.
///
/// The asymmetric chamfer: `face` names which side the first distance applies
/// to, and must be one of the two faces meeting at the edge.
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
    let spec = Chamfer::Angle {
        face: face.clone(),
        distance,
        angle,
    };
    wedge_for(model, solid, edge, &spec, tol)?.apply(model, solid, edge, tol)
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
    if specs.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges bevels nothing");
    }
    let wedges: Vec<Wedge> = specs
        .iter()
        .map(|(edge, spec)| wedge_for(model, solid, edge, spec, tol))
        .collect::<OgeomResult<_>>()?;
    let mut built = Built::from_nothing(solid.clone());
    for ((edge, _), wedge) in specs.iter().zip(wedges) {
        let step = wedge.apply(model, &built.shape, edge, tol)?;
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// The wedge one chamfer cuts, built on `solid` as it stands.
fn wedge_for(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    match spec {
        Chamfer::Symmetric(distance) => match seat_kind(model, edge, tol)? {
            SeatKind::Straight => {
                let seat = planar_seat(model, solid, edge, tol)?;
                bevel(model, &seat, [*distance, *distance], tol)
            }
            SeatKind::Rim(rim) => {
                let seat = revolved_seat(model, solid, edge, &rim, tol)?;
                revolved_bevel(model, &seat, *distance, *distance, tol)
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
                bevel(model, &seat, distances, tol)
            }
            SeatKind::Rim(rim) => {
                let seat = revolved_seat(model, solid, edge, &rim, tol)?;
                let (on_wall, on_cap) = revolved_side(&seat, face, *on_face, *on_other)?;
                revolved_bevel(model, &seat, on_wall, on_cap, tol)
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
                    // side — the material on a convex edge, the open dihedral
                    // on a concave one. Where it crosses the other leg's ray
                    // is the derived distance — no crossing, no chamfer.
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
                    bevel(model, &seat, distances, tol)
                }
                SeatKind::Rim(rim) => {
                    // The rim's seat is square by construction — the cap is
                    // perpendicular to the wall — so the derived distance is
                    // the plain tangent, and past a right angle the bevel
                    // walks away from the other face instead of toward it.
                    if angle >= core::f64::consts::FRAC_PI_2 - tol.angular() {
                        ogeom_bail!(
                            Construction,
                            "the bevel at that angle never meets the edge's other face"
                        );
                    }
                    let seat = revolved_seat(model, solid, edge, &rim, tol)?;
                    let derived = distance * angle.tan();
                    let (on_wall, on_cap) = revolved_side(&seat, face, distance, derived)?;
                    revolved_bevel(model, &seat, on_wall, on_cap, tol)
                }
            }
        }
    }
}

/// A chamfer's wedge, built and not yet applied: its faces, and whether it
/// fuses (a concave edge) or cuts (a convex one).
struct Wedge {
    faces: Vec<Shape>,
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
        apply_wedge(model, solid, Some(edge), &self.faces, self.additive, tol)
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
    if seat.wall_face.node() == face.node() {
        Ok((on_face, on_other))
    } else if seat.cap_face.node() == face.node() {
        Ok((on_other, on_face))
    } else {
        ogeom_bail!(
            Construction,
            "the named face does not meet the edge being chamfered"
        )
    }
}

/// The revolved wedge with the quarter-tube exchanged for a slant: legs
/// `on_wall` axially down the wall and `on_cap` radially along the cap, and
/// the cone between the two tangency rings as the bevel.
fn revolved_bevel(
    model: &mut Model,
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
    let flanks = revolved_flanks(model, seat, on_wall, cap_rho, tol)?;

    // The bevel: the cone through both tangency rings — reference radius
    // `cap_rho` at the cap's level, the rim's radius a wall-depth below.
    // Unlike the fillet's quarter-tube, whose away-from-the-tube normal
    // tracks the wedge seat by seat, the cone's natural normal always points
    // away from the axis — and the wedge sits on the axis side of the slant
    // exactly when `sigma` and `tau` agree.
    let bevel_band = {
        let slope = (cap_rho - seat.radius) / (seat.tau * on_wall);
        let cone = Cone::new(seat.frame_at(seat.centre, tol)?, cap_rho, slope.atan(), tol)?;
        // The domain covers the band's two rows — the cap ring at zero and
        // the wall ring a depth away — with a margin that stays clear of the
        // apex, where the surface degenerates.
        let rows = (
            0.0_f64.min(-seat.tau * on_wall),
            0.0_f64.max(-seat.tau * on_wall),
        );
        let pad = 0.1 * on_wall;
        let surface: SurfaceGeometry = ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into();
        let band = make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?;
        if seat.sigma * seat.tau > 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    Ok(Wedge {
        faces: vec![flanks.wall_band, flanks.annulus, bevel_band],
        additive: seat.additive(),
    })
}

/// Which side of the seat a named face is, by identity.
fn seat_side(seat: &Seat, face: &Shape) -> OgeomResult<usize> {
    if seat.faces[0].node() == face.node() {
        Ok(0)
    } else if seat.faces[1].node() == face.node() {
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
    // the material they will melt against with *opposed* orientation — which
    // is exactly what a fuse cancels.
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(0, tol)? * sign;
    let b = seat.leg(1, tol)? * sign;

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
    // same-domain resolution — which is what makes the coplanar legs melt
    // into the solid's own faces — recognises coincidence between *planes*.
    let faces = [
        planar_face(model, &[apex0, a0, b0], -seat.along, tol)?,
        planar_face(model, &[apex1, a1, b1], seat.along, tol)?,
        planar_face(model, &[apex0, a0, a1, apex1], seat.normals[0] * sign, tol)?,
        planar_face(model, &[apex0, b0, b1, apex1], seat.normals[1] * sign, tol)?,
        planar_face(model, &[a0, b0, b1, a1], bevel_out, tol)?,
    ];
    Ok(Wedge {
        faces: faces.to_vec(),
        additive: !seat.convex,
    })
}
