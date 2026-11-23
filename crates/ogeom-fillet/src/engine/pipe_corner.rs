//! The pipe where a chain's ball changes walls off its vertices (ADR 0001,
//! amendments "a ball changing walls off the vertex" and "the pipe's other
//! corners"): between the two sections where its contact crosses the seams
//! at a link's ends, the ball touches the cap and that link's wall alone,
//! and its centre runs on the curve at the radius from both, outside a
//! concave chain's faces and within a convex one's material. The round
//! there is the pipe of the radius about that curve, fitted through the
//! ball's exact sections; the neighbouring rounds end on those sections.
//! A crossing may stand in its vertex's own section, and the pipe then
//! ends there.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::SheetCurve;
use super::overflow::{ToolFace, boundary_from, closed_within, segment};
use super::switch_over::{Crossing, Faces, from_end, outward_at};
use crate::support::Wedge;

/// Sections fitted along the pipe.
const STATIONS: usize = 64;
/// Points across each section.
const ACROSS: usize = 33;
/// Steps of a two-hundredth of the radius the centre's march may take:
/// a hundred radii.
const MARCH: usize = 20_000;

/// A curve and the parameters of a piece of it.
type Crease = (Curve, (f64, f64));

/// One end of a pipe: where its link's neighbour hands the ball over.
pub(crate) struct End {
    /// The neighbour, by its index in the request.
    pub(crate) neighbour: usize,
    /// The ball's centre in the section the two meet in.
    pub(crate) centre: Point,
    /// The section's unit normal, from the neighbour into the pipe.
    pub(crate) into: Vector,
    /// Where the section crosses the neighbour's crease.
    pub(crate) foot: Point,
    /// The ball's contacts there: on the cap, and on the seam.
    pub(crate) on_cap: Point,
    pub(crate) touch: Point,
    /// The vertex.
    pub(crate) at: Point,
    /// The crossing stands in the vertex's own section: `foot` is the
    /// vertex, and the leg is the seam.
    in_section: bool,
    /// The section's leg on the neighbour's wall, from `foot` to `touch`.
    leg: Crease,
    /// The section's leg on the cap, from `on_cap` to `foot`.
    cap_leg: Crease,
    /// The fractions of the neighbour's curve range, start to end, its
    /// band keeps: from `foot` to its far end.
    pub(crate) kept: (f64, f64),
    /// The neighbour's crease, its parameters at the vertex and at `foot`.
    crease: Crease,
    /// The seam, its parameters at the vertex and at `touch`.
    seam: Crease,
    /// The neighbour's wall, by its index among the solid's faces.
    wall: usize,
}

/// The pipe of the request's `link`-th edge between the two crossings at
/// its ends, the neighbours there `ends[k].neighbour`; `None` where a
/// crossing off its vertex's section stands over the link itself, the two
/// differ in sense, a section's legs do not fit on their faces, the ball's
/// centre does not run between the two sections, or the fit or the tool
/// will not close.
///
/// # Errors
///
/// As reading the faces and edges and building the tool.
#[allow(clippy::too_many_lines, reason = "one tool, its pieces in turn")]
pub(crate) fn pipe(
    model: &mut Model,
    faces: &Faces,
    edges: &[(Shape, f64)],
    link: usize,
    crossings: [&Crossing; 2],
    tol: Tolerances,
) -> OgeomResult<Option<(Wedge, [End; 2])>> {
    let radius = edges[link].1;
    let cap = crossings[0].cap;
    let convex = crossings[0].convex;
    if crossings[1].cap != cap || crossings[1].convex != convex {
        return Ok(None);
    }
    // The ball stands off the faces outward of a concave chain, within the
    // material of a convex one.
    let side = if convex { -1.0 } else { 1.0 };
    let slot_of = |c: &Crossing| usize::from(c.pair[1] == link);
    let wall = crossings[0].walls[slot_of(crossings[0])];
    if crossings[1].walls[slot_of(crossings[1])] != wall {
        return Ok(None);
    }
    let host = |k: usize| crate::engine::flat_faces::host(model, &faces.faces[k], tol);
    let (Some(cap_host), Some(wall_host)) = (host(cap)?, host(wall)?) else {
        return Ok(None);
    };
    let mut ends = Vec::with_capacity(2);
    for crossing in crossings {
        let slot = slot_of(crossing);
        // The crossing stands past the vertex over the neighbour, the link's
        // ball running on there and the neighbour's ending at the crossing,
        // or in the vertex's own section, both ending there.
        let in_section = crossing.by <= super::switch_over::off_by(tol);
        if !in_section && crossing.over == slot {
            return Ok(None);
        }
        match end(
            model,
            faces,
            edges,
            crossing,
            slot,
            &cap_host,
            side * radius,
            tol,
        )? {
            Some(e) => ends.push(e),
            None => return Ok(None),
        }
    }
    let [e0, e1]: [End; 2] = match ends.try_into() {
        Ok(two) => two,
        Err(_) => return Ok(None),
    };

    // The centre's run from the one section to the other, on the cap's and
    // the wall's offsets at the radius.
    let Some(run) = centres(&cap_host, &wall_host, &e0, &e1, side * radius, tol)? else {
        return Ok(None);
    };
    let length: f64 = run.windows(2).map(|w| w[0].distance(w[1])).sum();
    let at_length = |s: f64| -> Point {
        let mut walked = 0.0;
        for w in run.windows(2) {
            let step = w[0].distance(w[1]);
            if walked + step >= s {
                let f = if step > 0.0 { (s - walked) / step } else { 0.0 };
                return w[0] + (w[1] - w[0]) * f;
            }
            walked += step;
        }
        run[run.len() - 1]
    };
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(STATIONS + 1);
    let mut mid = None;
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a few dozen stations")]
        let s = i as f64 / STATIONS as f64;
        let (c, b, w) = if i == 0 {
            (e0.centre, e0.on_cap, e0.touch)
        } else if i == STATIONS {
            (e1.centre, e1.on_cap, e1.touch)
        } else {
            let guess = at_length(length * 0.5 * (1.0 - (core::f64::consts::PI * s).cos()));
            let Some(found) = settled(&cap_host, &wall_host, guess, side * radius, tol)? else {
                return Ok(None);
            };
            found
        };
        let row = section(c, b, w, radius, tol)?;
        if i == STATIONS / 2 {
            mid = Some((c, b, w, row[ACROSS / 2]));
        }
        rows.push(row);
    }
    let Some((c_mid, b_mid, w_mid, p_mid)) = mid else {
        return Ok(None);
    };
    let target = tol.confusion() * 0.5;
    let Ok(fitted) = ogeom_geom::fit::fit_surface_grid_chordal(&rows, 3, target, tol) else {
        return Ok(None);
    };
    if fitted.error > tol.confusion() {
        return Ok(None);
    }
    let surface = fitted.curve;
    // Between its sections the fit holds the ball's radius about the run.
    if !holds(
        &surface,
        &run,
        radius,
        crate::marched::band_fit_target(tol),
        tol,
    )? {
        return Ok(None);
    }

    // Points: each end's contact on the cap and on the seam, its vertex,
    // and its foot on the neighbour's crease where that is not the vertex.
    let mut points = vec![e0.on_cap, e0.touch, e1.on_cap, e1.touch, e0.at, e1.at];
    let (b0, q0, b1, q1, v0, v1) = (0, 1, 2, 3, 4, 5);
    let mut foot_of = |e: &End, v: usize| {
        if e.in_section {
            v
        } else {
            points.push(e.foot);
            points.len() - 1
        }
    };
    let (f0, f1) = (foot_of(&e0, v0), foot_of(&e1, v1));
    let mut curves: Vec<SheetCurve> = Vec::new();
    let from_surface =
        |curves: &mut Vec<SheetCurve>, ends: (usize, usize)| -> OgeomResult<Option<usize>> {
            Ok(boundary_from(&surface, &points, ends, tol)?.map(|c| {
                curves.push(c);
                curves.len() - 1
            }))
        };
    let (Some(bottom), Some(top), Some(arc0), Some(arc1)) = (
        from_surface(&mut curves, (b0, b1))?,
        from_surface(&mut curves, (q0, q1))?,
        from_surface(&mut curves, (b0, q0))?,
        from_surface(&mut curves, (b1, q1))?,
    ) else {
        return Ok(None);
    };
    let push = |curves: &mut Vec<SheetCurve>,
                ends: (usize, usize),
                (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let bf0 = push(&mut curves, (b0, f0), e0.cap_leg.clone());
    let bf1 = push(&mut curves, (b1, f1), e1.cap_leg.clone());
    // A curve between two points by its parameters there.
    let piece = |(curve, (ta, tb)): &(Curve, (f64, f64)),
                 a: usize,
                 b: usize|
     -> (usize, usize, Curve, (f64, f64)) {
        if ta <= tb {
            (a, b, curve.clone(), (*ta, *tb))
        } else {
            (b, a, curve.clone(), (*tb, *ta))
        }
    };
    let push_piece =
        |curves: &mut Vec<SheetCurve>, (a, b, curve, range): (usize, usize, Curve, (f64, f64))| {
            curves.push(SheetCurve {
                ends: (a, b),
                curve,
                range,
            });
            curves.len() - 1
        };
    let seam0 = push_piece(&mut curves, piece(&e0.seam, v0, q0));
    let seam1 = push_piece(&mut curves, piece(&e1.seam, v1, q1));
    // An end in its vertex's section has no crease piece, and its leg on
    // the neighbour's wall is the seam.
    let leg_and_crease = |curves: &mut Vec<SheetCurve>, e: &End, f, q, v, seam| {
        if e.in_section {
            (seam, None)
        } else {
            let fq = push(curves, (f, q), e.leg.clone());
            (fq, Some(push_piece(curves, piece(&e.crease, v, f))))
        }
    };
    let (fq0, crease0) = leg_and_crease(&mut curves, &e0, f0, q0, v0, seam0);
    let (fq1, crease1) = leg_and_crease(&mut curves, &e1, f1, q1, v1, seam1);
    let Some((link_curve, link_range)) = ogeom_algo::placed_edge_curve(model, &edges[link].0, tol)?
    else {
        return Ok(None);
    };
    let (t_v0, t_v1) = from_end(&link_curve, link_range, e0.at, tol)?;
    let foot_link = push_piece(&mut curves, piece(&(link_curve, (t_v0, t_v1)), v0, v1));

    // The legs lie on the cap and the walls, the sections on their planes:
    // the tool is the fill outside a concave chain's faces, and the cut
    // within a convex one's.
    let section_plane = |e: &End| -> OgeomResult<SurfaceGeometry> {
        let reach = (radius * 4.0).max(1.0);
        let plane = Plane::through(e.centre, Direction::new(e.into, tol)?);
        Ok(ogeom_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into())
    };
    let (_, n_cap) = outward_at(&cap_host, b_mid, tol)?;
    let (_, n_wall) = outward_at(&wall_host, w_mid, tol)?;
    let mut tool_faces = vec![
        ToolFace {
            surface: surface.clone().into(),
            curves: vec![bottom, arc1, top, arc0],
            at: p_mid,
            out: (c_mid - p_mid) / radius,
        },
        ToolFace {
            surface: cap_host.0.clone(),
            curves: [
                Some(bottom),
                Some(bf1),
                crease1,
                Some(foot_link),
                crease0,
                Some(bf0),
            ]
            .into_iter()
            .flatten()
            .collect(),
            at: b_mid,
            out: n_cap * -side,
        },
        ToolFace {
            surface: wall_host.0.clone(),
            curves: vec![foot_link, seam1, top, seam0],
            at: w_mid,
            out: n_wall * -side,
        },
    ];
    for (e, crease, seam, fq, bf, arc) in [
        (&e0, crease0, seam0, fq0, bf0, arc0),
        (&e1, crease1, seam1, fq1, bf1, arc1),
    ] {
        if let Some(crease) = crease {
            let Some(neighbour) = host(e.wall)? else {
                return Ok(None);
            };
            let between = e.foot.midpoint(e.touch);
            let (_, n) = outward_at(&neighbour, between, tol)?;
            tool_faces.push(ToolFace {
                // Run on past its height range: the vertex stands within its
                // tolerance of the wall's end.
                surface: crate::marched::ruled_on(&neighbour.0, radius)?,
                curves: vec![crease, seam, fq],
                at: between,
                out: n * -side,
            });
        }
        tool_faces.push(ToolFace {
            surface: section_plane(e)?,
            curves: vec![bf, fq, arc],
            at: e.on_cap.midpoint(e.touch),
            out: -e.into,
        });
    }
    // The solid's own edges, and the fitted pipe's, meet within the bands'
    // fit.
    let widths = vec![crate::marched::band_fit_target(tol); curves.len()];
    let Some((solid, _)) = closed_within(model, &points, &curves, &tool_faces, &widths, tol)?
    else {
        return Ok(None);
    };
    Ok(Some((
        Wedge {
            solid,
            additive: !convex,
            edges: vec![edges[link].0.clone()],
            melt_refusal: None,
            credit: Vec::new(),
        },
        [e0, e1],
    )))
}

/// The end of a pipe at `crossing`, the link in `slot`: the ball's centre
/// off the neighbour's wall at the seam's crossing (`offset` outward, the
/// radius signed), the section through it square to the neighbour's run,
/// and where it cuts the neighbour's crease (the vertex, where the
/// crossing stands in its section). `None` where the section's legs on the
/// cap and the neighbour's wall will not fit, or an in-section crossing's
/// seam leaves the section.
#[allow(clippy::too_many_arguments, reason = "one end, read whole")]
#[allow(clippy::too_many_lines, reason = "one end, read whole")]
fn end(
    model: &Model,
    faces: &Faces,
    edges: &[(Shape, f64)],
    crossing: &Crossing,
    slot: usize,
    cap_host: &(SurfaceGeometry, f64),
    offset: f64,
    tol: Tolerances,
) -> OgeomResult<Option<End>> {
    let radius = offset.abs();
    let side = offset.signum();
    let other = 1 - slot;
    let neighbour = crossing.pair[other];
    let wall = crossing.walls[other];
    let Some(wall_host) = crate::engine::flat_faces::host(model, &faces.faces[wall], tol)? else {
        return Ok(None);
    };
    let Some(link_host) =
        crate::engine::flat_faces::host(model, &faces.faces[crossing.walls[slot]], tol)?
    else {
        return Ok(None);
    };
    // The crossing again, off one wall alone: its ball stands the radius
    // from the cap and from the seam's point. The two walls meet G1 at the
    // seam only within a fitted band's fit, so the wall read is the one
    // that is no fitted band: the neighbour's, unless only the link's is.
    let fitted = |h: &(SurfaceGeometry, f64)| matches!(h.0, SurfaceGeometry::BSpline(_));
    let solve_host = if fitted(&wall_host) && !fitted(&link_host) {
        &link_host
    } else {
        &wall_host
    };
    let (seam_curve, (t_v, t_q)) = &crossing.seam;
    let off = |t: f64| -> OgeomResult<(Point, Point, Point, Vector, f64)> {
        let touch = seam_curve.point_at(t, tol)?;
        let (_, n) = outward_at(solve_host, touch, tol)?;
        let centre = touch + n * offset;
        let (on_cap, up) = outward_at(cap_host, centre, tol)?;
        Ok((
            touch,
            centre,
            on_cap,
            up,
            side * (centre - on_cap).dot(up) - radius,
        ))
    };
    let mut t = *t_q;
    let mut h = off(t)?.4;
    let mut dt = (t_q - t_v) * 1e-3;
    for _ in 0..40 {
        if h.abs() <= tol.confusion() * 1e-3 || dt == 0.0 {
            break;
        }
        let h2 = off(t + dt)?.4;
        if h2 == h {
            break;
        }
        dt = -h * dt / (h2 - h);
        t += dt;
        h = off(t)?.4;
    }
    let (touch, centre, on_cap, up, _) = off(t)?;
    let seam = (seam_curve.clone(), (*t_v, t));
    if (side * (centre - on_cap).dot(up) - radius).abs() > crate::marched::band_fit_target(tol) {
        return Ok(None);
    }
    // The neighbour's run there is square to both contacts.
    let run = (centre - on_cap).cross(centre - touch);
    let m = run.magnitude();
    if m <= tol.confusion() {
        return Ok(None);
    }
    let near = crate::marched::band_fit_target(tol);
    let mut into = run / m;
    let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, &edges[neighbour].0, tol)?
    else {
        return Ok(None);
    };
    let (t_at, t_far) = from_end(&curve, range, crossing.at, tol)?;
    // From the neighbour into the pipe: from the crease's far end towards
    // the vertex.
    if (crossing.at - curve.point_at(t_far, tol)?).dot(into) < 0.0 {
        into = -into;
    }
    // A crossing in the vertex's own section ends the pipe there, its foot
    // the vertex and its leg on the neighbour's wall the seam, which lies in
    // that section.
    if crossing.by <= super::switch_over::off_by(tol) {
        let (_, (t_v, t_q)) = &seam;
        for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let q = seam.0.point_at(t_v + (t_q - t_v) * f, tol)?;
            if (q - centre).dot(into).abs() > near {
                return Ok(None);
            }
        }
        let Some(cap_leg) = leg(&cap_host.0, (centre, into), on_cap, crossing.at, tol)? else {
            return Ok(None);
        };
        return Ok(Some(End {
            neighbour,
            centre,
            into,
            foot: crossing.at,
            on_cap,
            touch,
            at: crossing.at,
            in_section: true,
            leg: seam.clone(),
            cap_leg,
            kept: (0.0, 1.0),
            crease: (curve, (t_at, t_at)),
            seam,
            wall,
        }));
    }
    // Where the section cuts the neighbour's crease, between the vertex
    // and its far end.
    let g = |t: f64| -> OgeomResult<f64> { Ok((curve.point_at(t, tol)? - centre).dot(into)) };
    if g(t_at)? <= 0.0 {
        return Ok(None);
    }
    let steps = 64;
    let mut last = t_at;
    let mut foot_t = None;
    for k in 1..=steps {
        let t = t_at + (t_far - t_at) * f64::from(k) / f64::from(steps);
        if g(t)? <= 0.0 {
            let (mut lo, mut hi) = (last, t);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if g(mid)? > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            foot_t = Some(0.5 * (lo + hi));
            break;
        }
        last = t;
    }
    let Some(foot_t) = foot_t else {
        return Ok(None);
    };
    let foot = curve.point_at(foot_t, tol)?;
    let Some(cap_leg) = leg(&cap_host.0, (centre, into), on_cap, foot, tol)? else {
        return Ok(None);
    };
    let Some(leg) = leg(&wall_host.0, (centre, into), foot, touch, tol)? else {
        return Ok(None);
    };
    Ok(Some(End {
        neighbour,
        centre,
        into,
        foot,
        on_cap,
        touch,
        at: crossing.at,
        in_section: false,
        leg,
        cap_leg,
        kept: {
            let f = (foot_t - range.0) / (range.1 - range.0);
            if t_at <= t_far { (f, 1.0) } else { (0.0, f) }
        },
        crease: (curve, (t_at, foot_t)),
        seam,
        wall,
    }))
}

/// The ball's section at centre `c` touching the cap at `b` and the wall
/// at `w`: `ACROSS` points of its arc from `b` to `w`.
fn section(c: Point, b: Point, w: Point, radius: f64, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let u = b - c;
    let v = w - c;
    let z = u.cross(v);
    let m = z.magnitude();
    if m <= tol.confusion() {
        ogeom_core::ogeom_bail!(Construction, "a pipe's section has no arc");
    }
    let z = z / m;
    let x = u / u.magnitude();
    let y = z.cross(x);
    let sweep = (x.dot(v) / v.magnitude()).clamp(-1.0, 1.0).acos();
    Ok((0..ACROSS)
        .map(|j| {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen points")]
            let phi = sweep * j as f64 / (ACROSS - 1) as f64;
            c + (x * phi.cos() + y * phi.sin()) * radius
        })
        .collect())
}

/// The ball's centre nearest `guess` at `offset` outward off the cap (the
/// radius signed) and the radius off the wall, with its contacts on each;
/// `None` where it does not settle.
fn settled(
    cap: &(SurfaceGeometry, f64),
    wall: &(SurfaceGeometry, f64),
    guess: Point,
    offset: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point, Point)>> {
    let radius = offset.abs();
    let mut c = guess;
    for _ in 0..40 {
        let (on_cap, up) = outward_at(cap, c, tol)?;
        c = on_cap + up * offset;
        let found = ogeom_geom::project_on_surface(&wall.0, c, 16, tol)?;
        let d = found.distance;
        if d <= tol.confusion() {
            return Ok(None);
        }
        let g = (c - found.point) / d;
        let along = g - up * g.dot(up);
        let m2 = along.dot(along);
        if m2 <= 1e-6 {
            return Ok(None);
        }
        if (d - radius).abs() <= tol.confusion() * 1e-3 {
            let (on_cap, _) = outward_at(cap, c, tol)?;
            return Ok(Some((c, on_cap, found.point)));
        }
        c -= along * ((d - radius) / m2);
    }
    Ok(None)
}

/// The centres from `e0`'s section to `e1`'s, marched in short steps on
/// the cap's and the wall's offsets (`offset` outward, the radius
/// signed); `None` where the march loses them.
fn centres(
    cap: &(SurfaceGeometry, f64),
    wall: &(SurfaceGeometry, f64),
    e0: &End,
    e1: &End,
    offset: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Point>>> {
    let step = offset.abs() / 200.0;
    let mut run = vec![e0.centre];
    let mut c = e0.centre;
    let mut heading = e0.into;
    for _ in 0..MARCH {
        let (_, up) = outward_at(cap, c, tol)?;
        let found = ogeom_geom::project_on_surface(&wall.0, c, 16, tol)?;
        let g = c - found.point;
        let mut t = up.cross(g);
        let m = t.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        t /= m;
        if t.dot(heading) < 0.0 {
            t = -t;
        }
        heading = t;
        let Some((next, _, _)) = settled(cap, wall, c + t * step, offset, tol)? else {
            return Ok(None);
        };
        // Past the far section: the run ends on its centre.
        if (next - e1.centre).dot(-e1.into) >= 0.0 {
            if next.distance(e1.centre) > step * 4.0 {
                return Ok(None);
            }
            run.push(e1.centre);
            return Ok(Some(run));
        }
        run.push(next);
        c = next;
    }
    Ok(None)
}

/// Whether `surface` stands within `fit` of the ball of `radius` about
/// `run` between its sections: read at points halfway between the grid's.
fn holds(
    surface: &ogeom_geom::BSplineSurface,
    run: &[Point],
    radius: f64,
    fit: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let to_run = |p: Point| -> f64 {
        run.windows(2)
            .map(|w| {
                let d = w[1] - w[0];
                let l2 = d.dot(d);
                let f = if l2 > 0.0 {
                    ((p - w[0]).dot(d) / l2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                p.distance(w[0] + d * f)
            })
            .fold(f64::INFINITY, f64::min)
    };
    for i in 0..2 * STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred samples")]
        let v = v0 + (v1 - v0) * (i as f64 + 0.5) / (2 * STATIONS) as f64;
        for j in 0..=8 {
            let u = u0 + (u1 - u0) * f64::from(j) / 8.0;
            let p = surface.point_at(u, v, tol)?;
            if (to_run(p) - radius).abs() > fit {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The leg from `foot` to `touch` where the section plane through `on`
/// square to `normal` cuts `wall`: the segment where it lies on the wall
/// within the bands' fit, and otherwise the cut itself, fitted; `None`
/// where the cut does not settle or will not fit.
pub(super) fn leg(
    wall: &SurfaceGeometry,
    (on, normal): (Point, Vector),
    foot: Point,
    touch: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Crease>> {
    let near = crate::marched::band_fit_target(tol);
    let mut straight = true;
    for f in [0.25, 0.5, 0.75] {
        let p = foot + (touch - foot) * f;
        straight &= ogeom_geom::project_on_surface(wall, p, 16, tol)?.distance <= near;
    }
    if straight {
        return segment(foot, touch, tol).map(Some);
    }
    // Each point of the chord taken to the cut: on the wall, then stepped
    // across it square to the cut onto the plane, and again.
    let cut = |f: f64| -> OgeomResult<Point> {
        if f <= 0.0 {
            return Ok(foot);
        }
        if f >= 1.0 {
            return Ok(touch);
        }
        let mut p = foot + (touch - foot) * f;
        for _ in 0..32 {
            let found = ogeom_geom::project_on_surface(wall, p, 16, tol)?;
            let (u, v) = found.parameters;
            let n = wall.normal_at(u, v, tol)?.vector();
            let across = n.cross(normal).cross(n);
            let off = (found.point - on).dot(normal);
            let lean = across.dot(normal);
            if lean.abs() <= 1e-9 {
                ogeom_core::ogeom_bail!(Construction, "a pipe's leg runs along its section");
            }
            p = found.point - across * (off / lean);
            if off.abs() <= tol.confusion() * 1e-3 && found.distance <= tol.confusion() * 1e-3 {
                break;
            }
        }
        let found = ogeom_geom::project_on_surface(wall, p, 16, tol)?;
        if found.distance > tol.confusion()
            || (found.point - on).dot(normal).abs() > tol.confusion()
        {
            ogeom_core::ogeom_bail!(Construction, "a pipe's leg does not settle on its wall");
        }
        Ok(found.point)
    };
    // Fitted through the cut's points, and held to it halfway between.
    let n = 32;
    let points = (0..=n)
        .map(|i| cut(f64::from(i) / f64::from(n)))
        .collect::<OgeomResult<Vec<_>>>();
    let Ok(points) = points else {
        return Ok(None);
    };
    let Ok(fitted) = ogeom_geom::fit::fit_points(&points, 3, tol.confusion() * 0.5, tol) else {
        return Ok(None);
    };
    if !fitted.met {
        return Ok(None);
    }
    for i in 0..n {
        let Ok(p) = cut((f64::from(i) + 0.5) / f64::from(n)) else {
            return Ok(None);
        };
        if ogeom_algo::project_on_curve(&Curve::BSpline(fitted.curve.clone()), p, 16, tol)?.distance
            > tol.confusion()
        {
            return Ok(None);
        }
    }
    let range = fitted.curve.domain();
    Ok(Some((Curve::BSpline(fitted.curve), range)))
}
