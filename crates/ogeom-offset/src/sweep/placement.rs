//! Bringing a profile onto a spine before it is swept.
//!
//! A pipe shell sweeps a profile standing square to the spine at its
//! start, through the start point. A profile drawn elsewhere, or leaning,
//! is brought there by a rigid motion about its centroid (the centre of
//! the area it bounds): turned by the least rotation that stands its plane
//! square to the spine's start tangent, and moved so its centroid lands on
//! the spine's start point.
//!
//! A pipe shell may instead bring its profile into contact itself, as
//! [`PipeContact`] asks: moved along its plane's normal onto the spine's
//! start, and turned about that point square to the start tangent. That
//! needs only the profile's plane, so it serves open wires as well as
//! faces and closed wires, and it keeps the profile's offset from the
//! spine within its plane.

use ogeom_algo::{Built, find_plane, surface_properties};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, PlaneSurface};
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType};

use super::{moved_ring, realized_profile};

/// How a profile is brought onto a spine by [`place_on_spine`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PipePlacement {
    /// Move the profile, unturned, so its centroid stands on the spine's
    /// start point.
    pub contact: bool,
    /// Turn the profile about its centroid, by the least rotation, until
    /// its plane stands square to the spine's tangent at its start.
    pub correction: bool,
}

/// How a pipe shell brings its profile onto the spine before sweeping it,
/// a rigid motion of the profile defined by the profile's plane and the
/// spine's start alone.
///
/// The contact point is the foot of the spine's start point on the
/// profile's plane. It is one point for any planar profile, an open wire
/// as much as a face, and a face and its outer wire share it; the
/// profile's own point nearest the start would not be (every point of a
/// ring about the start is nearest at once, and a face's nearest point
/// can lie inside it where its boundary's cannot).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PipeContact {
    /// The profile is swept where it stands: square to the spine at its
    /// start and through the start point, or refused.
    #[default]
    AsPlaced,
    /// The profile is moved, unturned, along its plane's normal by the
    /// distance from its plane to the spine's start point, so that the
    /// contact point lands on the start point. A profile already within
    /// the confusion distance of the start is swept as it stands.
    Contact,
    /// As `Contact`, then the profile is turned about the spine's start
    /// point by the least rotation standing its plane square to the
    /// spine's start tangent: about the axis square to both, through the
    /// angle between the tangent and the plane's normal taken whichever
    /// way round the normal is nearer the tangent. A plane has no side, so
    /// a normal pointing back along the tangent is already square and the
    /// profile is not turned. A plane holding the tangent (to an angle's
    /// tolerance) is a quarter turn from square either way round, and is
    /// refused.
    Correction,
}

/// A planar profile (a face, or a closed wire) moved rigidly onto the start
/// of `spine` as `placement` asks: turned about its centroid until it
/// stands square to the spine's start tangent where `correction`, and
/// moved so its centroid is the spine's start point where `contact`. The
/// result is a fresh profile, ready for [`make_pipe_shell`](super::make_pipe_shell)
/// and the sweeps like it; the history records it as the profile
/// modified.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the profile is not a planar face or closed wire, or the spine is not an
/// edge or a wire with a tangent at its start.
pub fn place_on_spine(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    placement: PipePlacement,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let kind = model.kind_of(profile)?;
    if !matches!(kind, ShapeType::Face | ShapeType::Wire) {
        ogeom_bail!(
            Construction,
            "a profile is a planar face or closed wire, not a {kind:?}"
        );
    }
    let Some(plane) = find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let (start, tangent) = spine_start(model, spine, tol)?;
    let centre = centroid(model, profile, kind, &plane, tol)?;
    let mut motion = Transform::IDENTITY;
    if placement.correction
        && let Some(turn) = least_turn(plane.normal().vector(), tangent, centre, tol)?
    {
        motion = turn;
    }
    if placement.contact {
        motion = Transform::translation(start - centre) * motion;
    }
    let placed = match kind {
        ShapeType::Face => realized_profile(model, profile, &motion, tol)?,
        _ => moved_ring(model, profile, &motion, tol)?,
    };
    let mut built = Built::new(placed.clone(), ogeom_algo::History::new());
    built.history.modify(profile, placed);
    Ok(built)
}

/// The least rotation about `about` standing a plane of normal `normal`
/// square to the unit `tangent`, or `None` where it stands square already
/// to an angle's tolerance. The plane faces whichever way round is nearer
/// the tangent.
fn least_turn(
    normal: Vector,
    tangent: Vector,
    about: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Transform>> {
    let n = if normal.dot(tangent) < 0.0 {
        -normal
    } else {
        normal
    };
    let axis = n.cross(tangent);
    if axis.magnitude() <= tol.angular() {
        return Ok(None);
    }
    let about = Axis {
        location: about,
        direction: Direction::new(axis, tol)?,
    };
    Ok(Some(Transform::rotation(
        about,
        axis.magnitude().atan2(n.dot(tangent)),
    )))
}

/// The profile a pipe shell sweeps under `contact`: `profile` itself where
/// it is swept as placed or the motion moves nothing, else a fresh profile
/// moved rigidly as [`PipeContact`] defines.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the profile is not a planar wire or face, if the spine has no tangent
/// at its start, or if a correction is asked of a profile whose plane
/// holds the start tangent.
pub(crate) fn profile_in_contact(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    contact: PipeContact,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if contact == PipeContact::AsPlaced {
        return Ok(profile.clone());
    }
    let kind = model.kind_of(profile)?;
    if !matches!(kind, ShapeType::Face | ShapeType::Wire) {
        ogeom_bail!(
            Construction,
            "a pipe shell sweeps a planar wire or face, not a {kind:?}"
        );
    }
    let Some(plane) = find_plane(model, profile, tol)? else {
        ogeom_bail!(
            Construction,
            "a profile is brought into contact along its plane's normal; this \
             one is not planar"
        );
    };
    let (start, tangent) = spine_start(model, spine, tol)?;
    let normal = plane.normal().vector();
    let shift = normal * (start - plane.origin()).dot(normal);
    let mut motion = Transform::translation(shift);
    let mut moves = shift.magnitude() > tol.confusion();
    if contact == PipeContact::Correction {
        if normal.dot(tangent).abs() <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the profile's plane holds the spine's start tangent: a quarter \
                 turn either way round stands it square, so its correction is \
                 not one turn"
            );
        }
        if let Some(turn) = least_turn(normal, tangent, start, tol)? {
            motion = turn * motion;
            moves = true;
        }
    }
    if !moves {
        return Ok(profile.clone());
    }
    match kind {
        ShapeType::Face => realized_profile(model, profile, &motion, tol),
        _ => moved_ring(model, profile, &motion, tol),
    }
}

/// The spine's start point and its unit tangent there, in the direction
/// of travel, read exactly off its first edge.
fn spine_start(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<(Point, Vector)> {
    let first = match model.kind_of(spine)? {
        ShapeType::Edge => spine.clone(),
        ShapeType::Wire => match model.ordered_children_of(spine)?.first() {
            Some(edge) => edge.clone(),
            None => ogeom_bail!(Construction, "the spine has no edge"),
        },
        other => ogeom_bail!(
            Construction,
            "a spine is an edge or a wire, not a {other:?}"
        ),
    };
    let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, &first, tol)? else {
        ogeom_bail!(Construction, "an edge has no curve");
    };
    let reversed = first.orientation() == ogeom_topo::Orientation::Reversed;
    let t = if reversed { range.1 } else { range.0 };
    let d = curve.d1_at(t, tol)?;
    let d = if reversed { -d } else { d };
    if d.magnitude() <= tol.confusion() {
        ogeom_bail!(Construction, "the spine has no tangent at its start");
    }
    Ok((curve.point_at(t, tol)?, d / d.magnitude()))
}

/// The centroid of the area a planar profile bounds: a face's own, or a
/// closed wire's as the face it bounds on `plane`.
fn centroid(
    model: &mut Model,
    profile: &Shape,
    kind: ShapeType,
    plane: &ogeom_math::Plane,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let face = if kind == ShapeType::Face {
        profile.clone()
    } else {
        ogeom_algo::make_face(
            model,
            PlaneSurface::new(*plane).into(),
            std::slice::from_ref(profile),
            tol,
        )?
        .shape
    };
    // Rebuilt with its pcurves on the plane, the face is integrated exactly
    // where its edges allow; else over a mesh as fine as the model's
    // confusion.
    let charted = realized_profile(model, &face, &Transform::IDENTITY, tol)?;
    let chord = Deflection::with_chord(tol.confusion())?;
    Ok(surface_properties(model, &charted, chord, tol)?.centre)
}
