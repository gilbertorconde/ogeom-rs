//! Draft: turning faces about a neutral plane so a part can leave its mould.
//!
//! A drafted face is the same face on a *tilted* support. It keeps the line
//! where it crosses the neutral plane — that line does not move, which is
//! what makes the draft measurable from a datum — and turns about it by the
//! draft angle. Everything else follows: the neighbouring faces re-meet the
//! tilted plane, the vertices re-solve, and the solid comes back with the
//! same topology on new geometry.
//!
//! That last part is not this module's work. It is the offset's rebuild,
//! which already puts a solid back together on moved supports; a draft
//! hands it turned surfaces instead of translated ones.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{PlaneSurface, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Plane, Point, Transform, Vector};
use ogeom_topo::{Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore_unique};

use crate::shape::rebuilt;

/// Draft the named faces of a solid about a neutral plane.
///
/// Each face turns about its own intersection with `neutral` by `angle`,
/// in the sense that leans the face inwards as it goes: a positive angle
/// narrows the solid in the `pull` direction — the way the part leaves its
/// mould — and a negative one widens it. Leaning inwards tilts the face's
/// outward normal *towards* the pull, which is how the sense is picked,
/// measured rather than assumed from a convention nobody can check. A face parallel
/// to the neutral plane has no line to turn about and is refused by name,
/// as is a face the rebuild cannot re-meet.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// named face is not a planar face of `solid`, is parallel to the neutral
/// plane, or the angle is not a usable one; plus whatever the rebuild
/// refuses.
pub fn apply_draft(
    model: &mut Model,
    solid: &Shape,
    faces: &[Shape],
    neutral: Plane,
    pull: Direction,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !angle.is_finite() || angle.abs() >= core::f64::consts::FRAC_PI_2 {
        ogeom_bail!(
            Construction,
            "a draft of {angle} radians turns the face past its own plane"
        );
    }
    if faces.is_empty() {
        ogeom_bail!(Construction, "a draft of no faces drafts nothing");
    }
    let own: Vec<TShapeId> = explore_unique(model, solid, ShapeType::Face)?
        .iter()
        .map(Shape::node)
        .collect();

    // The turned surface for each named face, worked out before the
    // rebuild, so a face that cannot be drafted says so here rather than
    // half-way through a solid.
    let mut turned: Vec<(TShapeId, SurfaceGeometry)> = Vec::with_capacity(faces.len());
    for face in faces {
        if !own.contains(&face.node()) {
            ogeom_bail!(Construction, "a drafted face is not a face of the solid");
        }
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
            ogeom_bail!(Construction, "expected a face");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let SurfaceGeometry::Plane(p) = surface else {
            ogeom_bail!(
                Construction,
                "drafting a curved face needs a rebuild that can turn its \
                 support; this is the planar form"
            );
        };
        let plane = p.plane();
        let ((u0, u1), (v0, v1)) = surface.domain();
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        let outward = plane.normal().vector() * sign;

        // The hinge: the line where this face crosses the neutral plane.
        let along = plane.normal().vector().cross(neutral.normal().vector());
        let magnitude = along.magnitude();
        if magnitude <= tol.angular() {
            ogeom_bail!(
                Construction,
                "a face parallel to the neutral plane has no line to turn \
                 about"
            );
        }
        let along = along / magnitude;
        let hinge = meet(plane, neutral, along, tol)?;

        // Which way to turn: the sense whose outward normal ends up leaning
        // furthest towards the pull, which is the face leaning inwards and
        // the solid narrowing as it leaves.
        let axis = ogeom_math::Axis::new(hinge, Direction::new(along, tol)?);
        let mut candidates = Vec::with_capacity(2);
        for sense in [1.0, -1.0] {
            let turn = Transform::rotation(axis, angle * sense);
            candidates.push((sense, turn.apply_vector(outward).dot(pull.vector())));
        }
        let leaning = candidates
            .iter()
            .copied()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(core::cmp::Ordering::Equal))
            .map_or(1.0, |(sense, _)| sense);
        let turn = Transform::rotation(axis, angle * leaning);
        let moved_normal = Direction::new(turn.apply_vector(plane.normal().vector()), tol)?;
        let tilted = Plane::new(Frame::new(
            hinge,
            moved_normal,
            Direction::new(along, tol)?,
            tol,
        )?);
        // The window grows with the turn: a tilted plane reaches further
        // across the same solid than the one it replaces.
        let grow = (u1 - u0).abs().max((v1 - v0).abs()).mul_add(0.5, 1.0) * angle.abs().tan()
            + tol.confusion();
        turned.push((
            face.node(),
            PlaneSurface::over(tilted, (u0 - grow, u1 + grow), (v0 - grow, v1 + grow))?.into(),
        ));
    }

    rebuilt(
        model,
        solid,
        &|_| 0.0,
        &|face| {
            turned
                .iter()
                .find(|(node, _)| *node == face.node())
                .map(|(_, surface)| surface.clone())
        },
        tol,
    )
}

/// A point on the line where two planes meet, nearest their origins.
fn meet(a: Plane, b: Plane, along: Vector, tol: Tolerances) -> OgeomResult<Point> {
    let rows = [a.normal().vector(), b.normal().vector(), along];
    let rhs = [
        rows[0].dot(a.origin().to_vector()),
        rows[1].dot(b.origin().to_vector()),
        along.dot(Point::midpoint(a.origin(), b.origin()).to_vector()),
    ];
    let det = rows[0].dot(rows[1].cross(rows[2]));
    if det.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "the two planes do not meet in a line");
    }
    Ok(Point::ORIGIN
        + (rows[1].cross(rows[2]) * rhs[0]
            + rows[2].cross(rows[0]) * rhs[1]
            + rows[0].cross(rows[1]) * rhs[2])
            / det)
}
