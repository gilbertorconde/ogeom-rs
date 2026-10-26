//! Bringing a profile onto a spine before it is swept.
//!
//! A pipe shell sweeps a profile standing square to the spine at its
//! start, through the start point. A profile drawn elsewhere, or leaning,
//! is brought there by a rigid motion about its centroid (the centre of
//! the area it bounds): turned by the least rotation that stands its plane
//! square to the spine's start tangent, and moved so its centroid lands on
//! the spine's start point.

use ogeom_algo::{Built, find_plane, surface_properties};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, PlaneSurface};
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType};

use super::{moved_ring, realized_profile, spine_curve_of};

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
    if placement.correction {
        // The least rotation standing the plane square to the tangent turns
        // whichever way round the plane faces nearer it.
        let n = plane.normal().vector();
        let n = if n.dot(tangent) < 0.0 { -n } else { n };
        let axis = n.cross(tangent);
        if axis.magnitude() > tol.angular() {
            let about = Axis {
                location: centre,
                direction: Direction::new(axis, tol)?,
            };
            motion = Transform::rotation(about, axis.magnitude().atan2(n.dot(tangent)));
        }
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
    let (curve, range) = spine_curve_of(model, &first)?;
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
