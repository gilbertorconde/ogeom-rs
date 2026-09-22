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
        let axial = |frame: Frame| -> bool {
            (frame.z().vector().dot(neutral.normal().vector()).abs() - 1.0).abs()
                <= tol.angular().max(1e-9)
        };
        match surface {
            SurfaceGeometry::Cylinder(c) if axial(c.cylinder().frame()) => {
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
            SurfaceGeometry::Cone(co) if axial(co.cone().frame()) => {
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
            // Everything else — a raw fitted patch, a wall of revolution
            // about an oblique neutral — is drafted the way a mould-maker
            // drafts: along the pull, tilted by the angle, through the
            // line where the face crosses the neutral plane.
            _ => {
                turned.push((
                    face.node(),
                    general_draft(model, face, surface, sign, neutral, pull, angle, tol)?,
                ));
                continue;
            }
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

/// How many stations a general draft's hinge is sampled at, all round.
const HINGE_STATIONS: usize = 256;

/// The turned support for any face: the ruled surface through the face's
/// crossing with the neutral plane, its rulings the pull direction turned
/// about the crossing's tangent by the draft — what a mould-maker means by
/// a draft, and what the planar and revolved paths are the closed forms
/// of. The crossing is read off the face's own mesh and corrected onto the
/// surface, so it lies inside the face whatever the surface's chart does.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn general_draft(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    sign: f64,
    neutral: Plane,
    pull: Direction,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let n = neutral.normal().vector();
    let mesh = ogeom_mesh::triangulate_face(
        model,
        face,
        ogeom_mesh::Deflection {
            chord: 0.05,
            angular: 0.2,
            ..ogeom_mesh::Deflection::default()
        },
        tol,
    )?;
    let extent = {
        let b = mesh
            .positions
            .iter()
            .fold(ogeom_math::Aabb::EMPTY, |acc, p| acc.with_point(*p));
        match (b.low(), b.high()) {
            (Some(lo), Some(hi)) => (hi - lo).magnitude(),
            _ => ogeom_bail!(Construction, "the drafted face has no extent"),
        }
    };
    if extent <= tol.confusion() {
        ogeom_bail!(Construction, "the drafted face has no extent");
    }

    // The crossing, one segment per triangle the plane cuts, in the chart
    // with its ends in space. A vertex the plane passes through — the hinge
    // running along a rim, as a draft about a base does — is a crossing in
    // itself, not a sign to read; read as one, rounding gives every rim
    // triangle a hair of a segment pointing anywhere.
    let on = tol.confusion() * 10.0;
    let side: Vec<f64> = mesh
        .positions
        .iter()
        .map(|p| neutral.signed_distance_to(*p))
        .collect();
    let mut segments: Vec<[((f64, f64), Point); 2]> = Vec::new();
    for t in &mesh.triangles {
        let mut ends: Vec<((f64, f64), Point)> = Vec::with_capacity(2);
        for &corner in t {
            let i = corner as usize;
            if side[i].abs() <= on {
                ends.push((mesh.parameters[i], mesh.positions[i]));
            }
        }
        for k in 0..3 {
            let (i, j) = (t[k] as usize, t[(k + 1) % 3] as usize);
            let (a, b) = (side[i], side[j]);
            if a.abs() <= on || b.abs() <= on || (a < 0.0) == (b < 0.0) {
                continue;
            }
            let f = a / (a - b);
            let (pa, pb) = (mesh.parameters[i], mesh.parameters[j]);
            let (qa, qb) = (mesh.positions[i], mesh.positions[j]);
            ends.push((
                (pa.0 + (pb.0 - pa.0) * f, pa.1 + (pb.1 - pa.1) * f),
                qa + (qb - qa) * f,
            ));
        }
        ends.dedup_by(|a, b| a.1.distance(b.1) <= on);
        if ends.len() == 2 && ends[0].1.distance(ends[1].1) > on {
            segments.push([ends[0], ends[1]]);
        }
    }
    if segments.is_empty() {
        ogeom_bail!(
            Construction,
            "the neutral plane does not cross the drafted face; there is no \
             hinge to turn about"
        );
    }
    // Chained end to end into one run, closed or open, by where the ends
    // are in space — across a seam the chart says two things and space one.
    // Two runs is a face the plane crosses twice, which has no one hinge.
    // A rim's segments come once per triangle on either side of it, so a
    // segment already covered by the chain is dropped rather than chained.
    let same = |a: Point, b: Point| a.distance(b) <= on;
    let mut chain: Vec<((f64, f64), Point)> = vec![segments[0][0], segments[0][1]];
    let mut used = vec![false; segments.len()];
    used[0] = true;
    loop {
        let tail = chain[chain.len() - 1].1;
        let head = chain[0].1;
        let mut grew = false;
        for (k, seg) in segments.iter().enumerate() {
            if used[k] {
                continue;
            }
            // The segment that joins the tail back to the head closes the
            // run; it is kept, and the run is closed by it.
            if (same(seg[0].1, tail) && same(seg[1].1, head))
                || (same(seg[1].1, tail) && same(seg[0].1, head))
            {
                chain.push(chain[0]);
                used[k] = true;
                grew = true;
                break;
            }
            let covered = |p: Point| chain.iter().any(|c| same(c.1, p));
            if covered(seg[0].1) && covered(seg[1].1) {
                used[k] = true;
                continue;
            }
            if same(seg[0].1, tail) {
                chain.push(seg[1]);
            } else if same(seg[1].1, tail) {
                chain.push(seg[0]);
            } else if same(seg[0].1, head) {
                chain.insert(0, seg[1]);
            } else if same(seg[1].1, head) {
                chain.insert(0, seg[0]);
            } else {
                continue;
            }
            used[k] = true;
            grew = true;
        }
        if !grew {
            break;
        }
    }
    if used.iter().any(|u| !u) {
        ogeom_bail!(
            Construction,
            "the neutral plane crosses the drafted face more than once; \
             there is no one hinge to turn about"
        );
    }
    let closed = chain.len() > 3 && same(chain[0].1, chain[chain.len() - 1].1);
    if closed {
        chain.pop();
    }
    if chain.len() < 2 {
        ogeom_bail!(
            Construction,
            "the neutral plane touches the drafted face at a point; there is \
             no hinge to turn about"
        );
    }
    let chain: Vec<(f64, f64)> = chain.into_iter().map(|c| c.0).collect();

    // The mesh's stations are a chord apart; a cubic fitted through them
    // sits a fraction of that chord off the true hinge. Resampled between
    // them in the chart — across a periodic seam by the short way — and
    // corrected onto the surface below, the stations are as dense as the
    // fit's target wants.
    let chain: Vec<(f64, f64)> = {
        let ((ua, ub), (va, vb)) = surface.domain();
        // Closed counts as periodic here: a fitted tube closes on itself
        // without repeating, and a chain crossing its join must still take
        // the short way round.
        let period = (
            (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(ub - ua),
            (surface.is_periodic_v() || surface.is_closed_v(tol)).then_some(vb - va),
        );
        let short = |a: f64, b: f64, period: Option<f64>| -> f64 {
            let d = b - a;
            match period {
                Some(p) if d.abs() > p * 0.5 => d - p * d.signum(),
                _ => d,
            }
        };
        // A closed hinge starts where it crosses the chart's own seam
        // column, so the drafted support's seam stands where the old one
        // did and the rebuild finds it there.
        let chain: Vec<(f64, f64)> = if closed {
            // The station exactly on the column: where the chain's segments
            // cross `u = ua`, interpolated there; the nearest chain point
            // otherwise.
            let n = chain.len();
            let mut exact: Option<(usize, (f64, f64))> = None;
            for i in 0..n {
                let (a, b) = (chain[i], chain[(i + 1) % n]);
                let (da, db) = (short(ua, a.0, period.0), short(ua, b.0, period.0));
                if da == 0.0 {
                    exact = Some((i, a));
                    break;
                }
                if (da < 0.0) != (db < 0.0) && (da - db).abs() > 0.0 {
                    let f = da / (da - db);
                    let dv = short(a.1, b.1, period.1);
                    exact = Some((i + 1, (ua, a.1 + dv * f)));
                    break;
                }
            }
            let (start, inserted) = exact.unwrap_or_else(|| {
                let mut best = (0usize, f64::INFINITY);
                for (i, c) in chain.iter().enumerate() {
                    let d = short(ua, c.0, period.0).abs();
                    if d < best.1 {
                        best = (i, d);
                    }
                }
                (best.0, chain[best.0])
            });
            let mut rotated: Vec<(f64, f64)> = Vec::with_capacity(n + 1);
            rotated.push(inserted);
            for k in 0..n {
                let c = chain[(start + k) % n];
                if rotated.len() == 1 && c == inserted {
                    continue;
                }
                rotated.push(c);
            }
            rotated
        } else {
            chain
        };
        // Evenly by arc length, so the fit's parameter is the hinge's own
        // length: the mesh's segments run from a hair to a chord, and a
        // cubic through stations that uneven wanders between them.
        let pairs = if closed { chain.len() } else { chain.len() - 1 };
        let mut lengths = Vec::with_capacity(pairs);
        let mut total = 0.0;
        for i in 0..pairs {
            let (a, b) = (chain[i], chain[(i + 1) % chain.len()]);
            let step = surface
                .point_at(a.0, a.1, tol)?
                .distance(surface.point_at(b.0, b.1, tol)?);
            lengths.push(step);
            total += step;
        }
        let count = if closed {
            HINGE_STATIONS
        } else {
            HINGE_STATIONS + 1
        };
        let mut dense = Vec::with_capacity(count);
        let (mut pair, mut walked) = (0usize, 0.0_f64);
        for k in 0..count {
            #[allow(clippy::cast_precision_loss)]
            let target = total * k as f64 / HINGE_STATIONS as f64;
            while pair + 1 < pairs && walked + lengths[pair] < target {
                walked += lengths[pair];
                pair += 1;
            }
            let (a, b) = (chain[pair], chain[(pair + 1) % chain.len()]);
            let (du, dv) = (short(a.0, b.0, period.0), short(a.1, b.1, period.1));
            let f = if lengths[pair] > 0.0 {
                ((target - walked) / lengths[pair]).clamp(0.0, 1.0)
            } else {
                0.0
            };
            dense.push((a.0 + du * f, a.1 + dv * f));
        }
        dense
    };

    // Each station corrected onto the surface's crossing with the plane —
    // the mesh's chord is not the surface — and read for its tangent and
    // outward normal.
    let mut hinges: Vec<Point> = Vec::with_capacity(chain.len());
    let mut tangents: Vec<Vector> = Vec::with_capacity(chain.len());
    let mut outwards: Vec<Vector> = Vec::with_capacity(chain.len());
    for (k, &(mut u, mut v)) in chain.iter().enumerate() {
        // The first station of a closed hinge is the seam column's, and
        // stays on it: corrected along `v` alone.
        let pinned = closed && k == 0;
        for _ in 0..8 {
            let p = surface.point_at(u, v, tol)?;
            let f = neutral.signed_distance_to(p);
            if f.abs() <= tol.confusion() * 1e-2 {
                break;
            }
            let (du, dv) = surface.d1_at(u, v, tol)?;
            let g = (if pinned { 0.0 } else { n.dot(du) }, n.dot(dv));
            let g2 = g.0 * g.0 + g.1 * g.1;
            if g2 <= 0.0 {
                break;
            }
            u -= f * g.0 / g2;
            v -= f * g.1 / g2;
        }
        let p = surface.point_at(u, v, tol)?;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let raw = du.cross(dv);
        if raw.magnitude() <= tol.confusion() {
            ogeom_bail!(Construction, "the drafted face has no normal on its hinge");
        }
        let outward = raw / raw.magnitude() * sign;
        // Along the crossing: square to both normals, run the way the chain
        // runs.
        let mut t = outward.cross(n);
        let next = chain[(k + 1) % chain.len()];
        let prev = chain[(k + chain.len() - 1) % chain.len()];
        let ahead =
            surface.point_at(next.0, next.1, tol)? - surface.point_at(prev.0, prev.1, tol)?;
        if t.dot(ahead) < 0.0 {
            t = -t;
        }
        if t.magnitude() <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the neutral plane is tangent to the drafted face; there is no \
                 hinge to turn about"
            );
        }
        hinges.push(p);
        tangents.push(t / t.magnitude());
        outwards.push(outward);
    }

    // The sense, probed at the middle station as the other paths probe:
    // the turn whose outward normal leans furthest towards the pull is the
    // inward one, and the caller's sign picks inward or outward through it.
    let m = hinges.len() / 2;
    let axis_m = ogeom_math::Axis::new(hinges[m], Direction::new(tangents[m], tol)?);
    let mut leaning = 1.0;
    let mut best = f64::NEG_INFINITY;
    for sense in [1.0_f64, -1.0] {
        let turn = Transform::rotation(axis_m, angle.abs() * sense);
        let lean = turn.apply_vector(outwards[m]).dot(pull.vector());
        if lean > best {
            best = lean;
            leaning = sense;
        }
    }
    let theta = angle * leaning;

    // One ruling per station: the pull turned about the hinge's tangent.
    let mut rulings: Vec<Vector> = Vec::with_capacity(hinges.len());
    for (hinge, tangent) in hinges.iter().zip(&tangents) {
        let turn = Transform::rotation(
            ogeom_math::Axis::new(*hinge, Direction::new(*tangent, tol)?),
            theta,
        );
        rulings.push(turn.apply_vector(pull.vector()));
    }
    // How far the face reaches along the pull either side of its hinge,
    // grown with the turn so the wall still reaches the neighbours it
    // re-meets.
    // Measured from every hinge station, not the middle one: an oblique
    // hinge rises and falls along the pull, and a ruling from its low end
    // must still reach the face's top.
    let (mut s_lo, mut s_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut h_lo, mut h_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for h in &hinges {
        let s = h.to_vector().dot(pull.vector());
        h_lo = h_lo.min(s);
        h_hi = h_hi.max(s);
    }
    for p in &mesh.positions {
        let s = p.to_vector().dot(pull.vector());
        s_lo = s_lo.min(s - h_hi);
        s_hi = s_hi.max(s - h_lo);
    }
    let grow = extent.mul_add(0.5, 1.0) * angle.abs().tan() + tol.confusion();
    let (s_lo, s_hi) = (s_lo - grow, s_hi + grow);
    if !closed {
        // An open hinge is continued straight past both ends, for the
        // same reason.
        let extend = |hinges: &mut Vec<Point>, rulings: &mut Vec<Vector>, front: bool| {
            let (i0, i1) = if front {
                (0, 1)
            } else {
                (hinges.len() - 1, hinges.len() - 2)
            };
            let d1 = hinges[i0] - hinges[i1];
            let steps = (grow / d1.magnitude().max(tol.confusion())).ceil().max(2.0);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let steps = (steps as usize).min(16);
            for k in 1..=steps {
                #[allow(clippy::cast_precision_loss)]
                let station = (hinges[i0] + d1 * k as f64, rulings[i0]);
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
    } else {
        hinges.push(hinges[0]);
        rulings.push(rulings[0]);
    }
    for edge in [s_lo, s_hi] {
        for i in 0..hinges.len() - 1 {
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
    // The ruled surface itself, exactly: the hinge fitted as a cubic, the
    // rulings' tips fitted as another at the *same* parameters, the two
    // brought onto one knot vector, and the surface linear between them —
    // degree one along the ruling, so a straight line is a straight line
    // and the fit's only error is the two curves' own. A grid fit through
    // rows at several heights parameterizes each row by its own chord and
    // averages, and rows that converge along their rulings disagree by
    // enough for a cubic across them to wander.
    // Chord-length parameters, not centripetal: the stations are as evenly
    // spaced as the mesh's segments let them be, and centripetal
    // parameters kink wherever the spacing changes, which a cubic then
    // cannot follow. By chord the parameter is the arc length whatever
    // the spacing.
    let params: Vec<f64> = {
        let mut out = Vec::with_capacity(hinges.len());
        let mut total = 0.0;
        out.push(0.0);
        for pair in hinges.windows(2) {
            total += pair[0].distance(pair[1]);
            out.push(total);
        }
        if total > 0.0 {
            for t in &mut out {
                *t /= total;
            }
        }
        if let Some(last) = out.last_mut() {
            *last = 1.0;
        }
        out
    };
    let reach = s_lo.abs().max(s_hi.abs()).max(1.0);
    let fit_target = (tol.confusion() * 1e3).max(1e-4);
    let hinge_fit = ogeom_geom::fit::fit_points_at(&params, &hinges, 3, fit_target, tol)?;
    let tips: Vec<Point> = hinges.iter().zip(&rulings).map(|(h, r)| *h + *r).collect();
    let tip_fit = ogeom_geom::fit::fit_points_at(&params, &tips, 3, fit_target / reach, tol)?;
    if !hinge_fit.met || !tip_fit.met {
        ogeom_bail!(
            NotDone,
            "the drafted wall's hinge fit reached {} and its rulings' {} against \
             a target of {fit_target}",
            hinge_fit.error,
            tip_fit.error * reach
        );
    }
    let (mut hinge_curve, mut tip_curve) = (hinge_fit.curve, tip_fit.curve);
    for (value, count) in tip_curve.knots().distinct() {
        let have = hinge_curve.knots().multiplicity_of(value);
        if count > have {
            hinge_curve = hinge_curve.with_knot_inserted(value, count - have, tol)?;
        }
    }
    for (value, count) in hinge_curve.knots().distinct() {
        let have = tip_curve.knots().multiplicity_of(value);
        if count > have {
            tip_curve = tip_curve.with_knot_inserted(value, count - have, tol)?;
        }
    }
    let (hc, tc) = (hinge_curve.control_points(), tip_curve.control_points());
    if hc.len() != tc.len() {
        ogeom_bail!(
            Construction,
            "the hinge and its rulings did not share a knot vector"
        );
    }
    let mut net: Vec<Point> = Vec::with_capacity(hc.len() * 2);
    for (h, t) in hc.iter().zip(tc) {
        let (h, d) = (h.point(), t.point() - h.point());
        net.push(h + d * s_lo);
        net.push(h + d * s_hi);
    }
    let grid = ogeom_math::ControlGrid::new(net, hc.len(), 2)?;
    // The chart: `u` over the old surface's own `u` domain for a closed
    // hinge, so the seam column carries over; `v` the height along the
    // ruling in the model's own units.
    let u_knots = if closed {
        let ((ua, ub), _) = surface.domain();
        hinge_curve.knots().reparameterized(ua, ub)?
    } else {
        hinge_curve.knots().clone()
    };
    let v_knots = ogeom_math::KnotVector::clamped_uniform(1, 2)?.reparameterized(s_lo, s_hi)?;
    Ok(ogeom_geom::BSplineSurface::new(u_knots, v_knots, &grid, tol)?.into())
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
