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
use ogeom_topo::{Model, NodeData, Shape, ShapeType, TShapeId};

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
    // The solid's own face occurrences, orientation and all: a handle a
    // caller got from a canonical exploration carries no use-orientation,
    // and the draft's sense probe needs the true outward.
    let own: Vec<Shape> = {
        let mut seen: Vec<Shape> = Vec::new();
        for f in ogeom_topo::explore(model, solid, ogeom_topo::Filter::OfType(ShapeType::Face))? {
            if !seen.iter().any(|s| s.node() == f.node()) {
                seen.push(f);
            }
        }
        seen
    };

    // The turned surface for each named face, worked out before the
    // rebuild, so a face that cannot be drafted says so here rather than
    // half-way through a solid.
    let mut turned: Vec<(TShapeId, SurfaceGeometry)> = Vec::with_capacity(faces.len());
    for face in faces {
        let Some(used) = own.iter().find(|f| f.node() == face.node()).cloned() else {
            ogeom_bail!(Construction, "a drafted face is not a face of the solid");
        };
        let face = &used;
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
            ogeom_bail!(Construction, "expected a face");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        // Which sign turns the raw surface normal outward, read from the
        // material itself rather than from an orientation flag a handle may
        // or may not carry: a step along the raw normal that lands inside
        // the solid means the raw normal points inward.
        let sign = outward_sign(model, solid, face, surface, tol)?;
        let sign_of = |_: &Shape| sign;
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
            SurfaceGeometry::Extrusion(e) => {
                turned.push((
                    face.node(),
                    extruded_draft(
                        e,
                        surface.domain(),
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
                "drafting a face that is neither planar, a wall of \
                 revolution nor an extruded wall needs a support the \
                 rebuild cannot turn — docs/PARITY.md, offset.draft"
            ),
        }
        let SurfaceGeometry::Plane(p) = surface else {
            unreachable!("the match above let only planes through");
        };
        let plane = p.plane();
        let ((u0, u1), (v0, v1)) = surface.domain();
        // The plane's raw normal is its du x dv; the measured sign turns it
        // outward.
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

/// The turned support for a drafted extruded wall: every ruling rotated
/// about the hinge curve's own tangent by the draft, the result re-fitted.
///
/// The hinge is where the wall crosses the neutral plane — one closed-form
/// height per profile parameter — and it does not move, exactly as a planar
/// draft's hinge line does not. Each ruling turns about the hinge's local
/// tangent in the sense that leans the outward normal towards the pull,
/// probed at the profile's midpoint the way the planar draft probes its
/// candidates. The turned rulings are sampled on a grid and fitted; a draft
/// whose rulings cross inside the drafted window — a concave profile turned
/// far enough to fold — is refused by name before anything is fitted.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn extruded_draft(
    extrusion: &ogeom_geom::ExtrusionSurface,
    window: ((f64, f64), (f64, f64)),
    sign: f64,
    neutral: Plane,
    pull: Direction,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    use ogeom_geom::Curve3d as _;
    let ((u0, u1), (v0, v1)) = window;
    let d = extrusion.direction().vector();
    let n = neutral.normal().vector();
    let den = n.dot(d);
    if den.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the neutral plane runs along the wall's rulings; there is no \
             hinge to turn about"
        );
    }
    let curve = extrusion.curve();
    let o = neutral.origin().to_vector();
    // Where the ruling through C(u) crosses the neutral plane, and which way
    // the hinge runs there.
    let height_at = |c: Point| n.dot(o - c.to_vector()) / den;
    let hinge_tangent = |cd: Vector| cd - d * (n.dot(cd) / den);

    // The sense, probed at the profile's midpoint exactly as the planar
    // draft probes its two candidates: the turn whose outward normal leans
    // furthest towards the pull is the inward one, and the caller's sign
    // picks inward or outward through it.
    let um = f64::midpoint(u0, u1);
    let cm = curve.point_at(um, tol)?;
    let cdm = curve.d1_at(um, tol)?;
    let hinge_m = cm + d * height_at(cm);
    let tangent_m = Direction::new(hinge_tangent(cdm), tol)?;
    let outward = {
        let nw = cdm.cross(d);
        nw / nw.magnitude() * sign
    };
    let axis_m = ogeom_math::Axis::new(hinge_m, tangent_m);
    let mut leaning = 1.0;
    let mut best = f64::NEG_INFINITY;
    for sense in [1.0_f64, -1.0] {
        let turn = Transform::rotation(axis_m, angle.abs() * sense);
        let lean = turn.apply_vector(outward).dot(pull.vector());
        if lean > best {
            best = lean;
            leaning = sense;
        }
    }
    let theta = angle * leaning;

    // One ruling per sample: the hinge point, the hinge tangent, and the
    // extrusion direction turned about it.
    const ALONG: usize = 65;
    let mut hinges: Vec<Point> = Vec::with_capacity(ALONG);
    let mut rulings: Vec<Vector> = Vec::with_capacity(ALONG);
    let (mut s_lo, mut s_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..ALONG {
        #[allow(clippy::cast_precision_loss)]
        let u = u0 + (u1 - u0) * (i as f64) / ((ALONG - 1) as f64);
        let c = curve.point_at(u, tol)?;
        let cd = curve.d1_at(u, tol)?;
        let h = height_at(c);
        let hinge = c + d * h;
        let tangent = Direction::new(hinge_tangent(cd), tol)?;
        let turn = Transform::rotation(ogeom_math::Axis::new(hinge, tangent), theta);
        hinges.push(hinge);
        rulings.push(turn.apply_vector(d));
        s_lo = s_lo.min(v0 - h);
        s_hi = s_hi.max(v1 - h);
    }
    // The window grows with the turn, as the planar draft's does. Along the
    // profile the curve itself ends, so the growth is a tangent-line
    // continuation at each end: the wall must still reach the neighbours it
    // re-meets, and the continuation exists only to be trimmed away there.
    let grow = (u1 - u0).abs().max((v1 - v0).abs()).mul_add(0.5, 1.0) * angle.abs().tan()
        + tol.confusion();
    let (s_lo, s_hi) = (s_lo - grow, s_hi + grow);
    {
        // Quadratic continuation, so the fitted wall keeps its end
        // curvature across the join instead of kinking straight.
        let extend = |hinges: &mut Vec<Point>, rulings: &mut Vec<Vector>, front: bool| {
            let (i0, i1, i2) = if front {
                (0, 1, 2)
            } else {
                let n = hinges.len();
                (n - 1, n - 2, n - 3)
            };
            let d1 = hinges[i0] - hinges[i1];
            let d2 = (hinges[i0] - hinges[i1]) - (hinges[i1] - hinges[i2]);
            let r1 = rulings[i0] - rulings[i1];
            let steps = (grow / d1.magnitude().max(tol.confusion())).ceil().max(2.0);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let steps = (steps as usize).min(16);
            for k in 1..=steps {
                #[allow(clippy::cast_precision_loss)]
                let k = k as f64;
                let station = (
                    hinges[i0] + d1 * k + d2 * (k * (k + 1.0) / 2.0),
                    rulings[i0] + r1 * k,
                );
                if front {
                    hinges.insert(0, station.0);
                    rulings.insert(0, station.1);
                } else {
                    hinges.push(station.0);
                    rulings.push(station.1);
                }
            }
        };
        extend(&mut hinges, &mut rulings, true);
        extend(&mut hinges, &mut rulings, false);
    }
    let along_total = hinges.len();

    // A fold is two rulings crossing inside the window: walking the wall at
    // either extreme height must still advance the way the hinge advances.
    for edge in [s_lo, s_hi] {
        for i in 0..along_total - 1 {
            let step = (hinges[i + 1] + rulings[i + 1] * edge) - (hinges[i] + rulings[i] * edge);
            if step.dot(hinges[i + 1] - hinges[i]) <= 0.0 {
                ogeom_bail!(
                    Construction,
                    "the draft folds the wall onto itself inside the drafted \
                     window; refused — docs/PARITY.md, offset.draft"
                );
            }
        }
    }

    // Rulings are straight, so a handful of rows fits them exactly; the
    // profile direction carries the shape.
    const ACROSS: usize = 9;
    let rows: Vec<Vec<Point>> = (0..ACROSS)
        .map(|j| {
            #[allow(clippy::cast_precision_loss)]
            let s = s_lo + (s_hi - s_lo) * (j as f64) / ((ACROSS - 1) as f64);
            (0..along_total)
                .map(|i| hinges[i] + rulings[i] * s)
                .collect()
        })
        .collect();
    let fit_target = (tol.confusion() * 1e3).max(1e-4);
    let fitted = ogeom_geom::fit::fit_surface_grid(&rows, 3, fit_target, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the drafted wall's fit reached {} against a target of {fit_target}",
            fitted.error
        );
    }
    Ok(fitted.curve.into())
}

/// Which sign turns a surface's raw normal (du x dv) outward, read from
/// the solid itself: probed a step off the face midpoint on both sides, at
/// growing steps until one side is material and the other is not.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if no
/// probe separates the sides — a wall thinner than the probe can resolve.
fn outward_sign(
    model: &Model,
    solid: &Shape,
    face: &Shape,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_algo::Containment;
    // A point genuinely on the face — the surface's domain midpoint may lie
    // outside the trim — from the face's own triangulation, at its largest
    // triangle's centre.
    let mesh = ogeom_mesh::triangulate_face(model, face, ogeom_mesh::Deflection::default(), tol)?;
    let mut at = None;
    let mut largest = 0.0_f64;
    for t in &mesh.triangles {
        let [a, b, c] = [
            mesh.positions[t[0] as usize],
            mesh.positions[t[1] as usize],
            mesh.positions[t[2] as usize],
        ];
        let area = (b - a).cross(c - a).magnitude();
        if area > largest {
            largest = area;
            let params = [
                mesh.parameters[t[0] as usize],
                mesh.parameters[t[1] as usize],
                mesh.parameters[t[2] as usize],
            ];
            at = Some((
                (params[0].0 + params[1].0 + params[2].0) / 3.0,
                (params[0].1 + params[1].1 + params[2].1) / 3.0,
            ));
        }
    }
    let Some((um, vm)) = at else {
        ogeom_bail!(Construction, "the drafted face has no interior to probe");
    };
    let p = surface.point_at(um, vm, tol)?;
    let (du, dv) = surface.d1_at(um, vm, tol)?;
    let n = du.cross(dv);
    let m = n.magnitude();
    if m <= tol.confusion() {
        ogeom_bail!(Construction, "the face has no normal at its midpoint");
    }
    let n = n / m;
    let scale = largest.sqrt().max(tol.confusion() * 1e3);
    for eps_scale in [1e-3, 1e-2, 5e-2] {
        let eps = scale * eps_scale;
        let deflection = ogeom_mesh::Deflection {
            chord: (eps * 0.1).max(1e-4),
            ..ogeom_mesh::Deflection::default()
        };
        let ahead = ogeom_algo::classify_in_solid(model, solid, p + n * eps, deflection, tol)?;
        let behind = ogeom_algo::classify_in_solid(model, solid, p - n * eps, deflection, tol)?;
        match (ahead, behind) {
            (Containment::Out, Containment::In) => return Ok(1.0),
            (Containment::In, Containment::Out) => return Ok(-1.0),
            _ => {}
        }
    }
    ogeom_bail!(
        Construction,
        "cannot read which side of the drafted face holds material; the          wall is thinner than the probe can resolve"
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
