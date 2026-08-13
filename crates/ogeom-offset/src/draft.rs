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
    let (canonical, mapped, prefix) = crate::shape::canonical_input(model, solid, faces, tol)?;
    if let Some(prefix) = prefix {
        let mut out = apply_draft(model, &canonical, &mapped, neutral, pull, angle, tol)?;
        out.history = prefix.then(&out.history);
        return Ok(out);
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
        let sign_of = |face: &Shape| {
            if face.orientation() == Orientation::Reversed {
                -1.0
            } else {
                1.0
            }
        };
        // A wall of revolution drafts about its neutral *circle*: the same
        // axis, the radius at the neutral plane held, the slant turned.
        match surface {
            SurfaceGeometry::Cylinder(c) => {
                let cylinder = c.cylinder();
                let (_, (v0, v1)) = surface.domain();
                turned.push((
                    face.node(),
                    revolved_draft(
                        cylinder.frame(),
                        cylinder.radius(),
                        0.0,
                        (v0, v1),
                        sign_of(face),
                        neutral,
                        pull,
                        angle,
                        tol,
                    )?,
                ));
                continue;
            }
            SurfaceGeometry::Cone(co) => {
                let cone = co.cone();
                let (_, (v0, v1)) = surface.domain();
                turned.push((
                    face.node(),
                    revolved_draft(
                        cone.frame(),
                        cone.reference_radius(),
                        cone.half_angle(),
                        (v0, v1),
                        sign_of(face),
                        neutral,
                        pull,
                        angle,
                        tol,
                    )?,
                ));
                continue;
            }
            SurfaceGeometry::Plane(_) => {}
            _ => ogeom_bail!(
                Construction,
                "drafting a face that is neither planar nor a wall of \
                 revolution needs a support the rebuild cannot turn — \
                 docs/PARITY.md, offset.draft"
            ),
        }
        let SurfaceGeometry::Plane(p) = surface else {
            unreachable!("the match above let only planes through");
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

        // Which way to turn: probed at the angle's *magnitude*, so the
        // sense names the inward lean — outward normal furthest towards the
        // pull, the solid narrowing as it leaves — and the angle's sign
        // stays the caller's: positive drafts inward, negative outward.
        let axis = ogeom_math::Axis::new(hinge, Direction::new(along, tol)?);
        let mut candidates = Vec::with_capacity(2);
        for sense in [1.0, -1.0] {
            let turn = Transform::rotation(axis, angle.abs() * sense);
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

/// The turned support for a drafted wall of revolution: a cone about the
/// same axis, holding the radius at the neutral plane and leaning the slant
/// by the draft, in the sense that tips the outward normal towards the pull.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn revolved_draft(
    frame: Frame,
    reference_radius: f64,
    half_angle: f64,
    window: (f64, f64),
    sign: f64,
    neutral: Plane,
    pull: Direction,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    use ogeom_geom::ConeSurface;

    let axis_dir = frame.z().vector();
    let along = axis_dir.dot(neutral.normal().vector());
    if (along.abs() - 1.0).abs() > tol.angular().max(1e-9) {
        ogeom_bail!(
            Construction,
            "a wall of revolution drafts about a neutral plane square to \
             its axis; the oblique neutral needs the general machinery — \
             docs/PARITY.md, offset.draft"
        );
    }
    // The neutral circle: where the axis meets the plane, and the radius
    // the wall holds there.
    let height = -neutral.signed_distance_to(frame.origin()) * along.signum();
    let neutral_point = frame.origin() + axis_dir * height;
    let neutral_radius = half_angle.tan().mul_add(height, reference_radius);
    if neutral_radius <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the wall has no radius left at the neutral plane to hold"
        );
    }
    let hinge_frame = Frame::new(neutral_point, frame.z(), frame.x(), tol)?;

    // Which way to lean, by measurement: of the two candidate slants, keep
    // the one whose outward normal — probed a little above the neutral
    // circle — ends up leaning furthest towards the pull.
    let mut best: Option<(f64, f64)> = None;
    for sense in [1.0_f64, -1.0] {
        // Probed at the magnitude: the sense names the inward lean, and the
        // caller's sign then picks inward or outward through it.
        let probe = half_angle + angle.abs() * sense;
        let candidate = half_angle + angle * sense;
        if probe.abs() <= tol.angular()
            || probe.abs() >= core::f64::consts::FRAC_PI_2 - tol.angular()
            || candidate.abs() <= tol.angular()
            || candidate.abs() >= core::f64::consts::FRAC_PI_2 - tol.angular()
        {
            continue;
        }
        let cone = ogeom_math::Cone::new(hinge_frame, neutral_radius, probe, tol)?;
        let surface: SurfaceGeometry = ConeSurface::new(cone, (-1.0, 1.0))?.into();
        let (du, dv) = surface.d1_at(0.0, 1.0, tol)?;
        let n = du.cross(dv);
        let outward = n / n.magnitude() * sign;
        let lean = outward.dot(pull.vector());
        if best.as_ref().is_none_or(|(_, held)| lean > *held) {
            best = Some((candidate, lean));
        }
    }
    let Some((leaned, _)) = best else {
        ogeom_bail!(
            Construction,
            "a draft of {angle} radians flattens the wall or swallows it"
        );
    };
    let cone = ogeom_math::Cone::new(hinge_frame, neutral_radius, leaned, tol)?;

    // The old window, re-expressed against the neutral origin and grown a
    // little; refused when the slant runs out of radius inside it.
    let shift = height;
    let grow = (window.1 - window.0).abs().mul_add(0.1, 1.0);
    let (w0, w1) = (window.0 - shift - grow, window.1 - shift + grow);
    let apex_height = -neutral_radius / leaned.tan();
    if apex_height > w0 && apex_height < w1 {
        ogeom_bail!(
            Construction,
            "the draft swallows the drafted face's own apex"
        );
    }
    Ok(ConeSurface::new(cone, (w0, w1))?.into())
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
