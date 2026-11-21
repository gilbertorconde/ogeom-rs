//! The pipe where a concave chain's ball changes walls off its vertices
//! (ADR 0001, amendment "a ball changing walls off the vertex"): between
//! the two sections where its contact crosses the seams at a link's ends,
//! the ball touches the cap and that link's wall alone, and its centre
//! runs on the curve at the radius from both. The round there is the pipe
//! of the radius about that curve, fitted through the ball's exact
//! sections; the neighbouring rounds end on those sections.

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
    /// The neighbour's crease, its parameters at the vertex and at `foot`.
    crease: Crease,
    /// The seam, its parameters at the vertex and at `touch`.
    seam: Crease,
    /// The neighbour's wall, by its index among the solid's faces.
    wall: usize,
}

/// The pipe of the request's `link`-th edge between the two crossings at
/// its ends, the neighbours there `ends[k].neighbour`; `None` where a
/// crossing stands over the link itself, a section's legs do not lie on
/// their faces, the ball's centre does not run between the two sections,
/// or the fit or the tool will not close.
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
    if crossings[1].cap != cap || crossings.iter().any(|c| c.convex) {
        return Ok(None);
    }
    let slot_of = |c: &Crossing| usize::from(c.pair[1] == link);
    let wall = crossings[0].walls[slot_of(crossings[0])];
    if crossings[1].walls[slot_of(crossings[1])] != wall {
        return Ok(None);
    }
    let host = |k: usize| crate::engine::rolled::host(model, &faces.faces[k], tol);
    let (Some(cap_host), Some(wall_host)) = (host(cap)?, host(wall)?) else {
        return Ok(None);
    };
    let mut ends = Vec::with_capacity(2);
    for crossing in crossings {
        let slot = slot_of(crossing);
        // The crossing stands past the vertex over the neighbour: the link's
        // ball runs on there, and the neighbour's ends at the crossing.
        if crossing.by <= super::switch_over::off_by(tol) || crossing.over == slot {
            return Ok(None);
        }
        match end(model, faces, edges, crossing, slot, &cap_host, radius, tol)? {
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
    let Some(run) = centres(&cap_host, &wall_host, &e0, &e1, radius, tol)? else {
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
            let Some(found) = settled(&cap_host, &wall_host, guess, radius, tol)? else {
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

    // Points: each end's contact on the cap and on the seam, its foot on
    // the neighbour's crease, and its vertex.
    let points = vec![
        e0.on_cap, e0.touch, e1.on_cap, e1.touch, e0.foot, e1.foot, e0.at, e1.at,
    ];
    let (b0, q0, b1, q1, f0, f1, v0, v1) = (0, 1, 2, 3, 4, 5, 6, 7);
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
    let bf0 = push(&mut curves, (b0, f0), segment(points[b0], points[f0], tol)?);
    let bf1 = push(&mut curves, (b1, f1), segment(points[b1], points[f1], tol)?);
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
    let fq0 = push(&mut curves, (f0, q0), segment(points[f0], points[q0], tol)?);
    let fq1 = push(&mut curves, (f1, q1), segment(points[f1], points[q1], tol)?);
    let crease0 = push_piece(&mut curves, piece(&e0.crease, v0, f0));
    let crease1 = push_piece(&mut curves, piece(&e1.crease, v1, f1));
    let seam0 = push_piece(&mut curves, piece(&e0.seam, v0, q0));
    let seam1 = push_piece(&mut curves, piece(&e1.seam, v1, q1));
    let Some((link_curve, link_range)) = ogeom_algo::placed_edge_curve(model, &edges[link].0, tol)?
    else {
        return Ok(None);
    };
    let (t_v0, t_v1) = from_end(&link_curve, link_range, e0.at, tol)?;
    let foot_link = push_piece(&mut curves, piece(&(link_curve, (t_v0, t_v1)), v0, v1));

    // The legs lie on the cap and the walls, the sections on their planes.
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
            curves: vec![bottom, bf1, crease1, foot_link, crease0, bf0],
            at: b_mid,
            out: -n_cap,
        },
        ToolFace {
            surface: wall_host.0.clone(),
            curves: vec![foot_link, seam1, top, seam0],
            at: w_mid,
            out: -n_wall,
        },
    ];
    for (e, crease, seam, fq, bf, arc) in [
        (&e0, crease0, seam0, fq0, bf0, arc0),
        (&e1, crease1, seam1, fq1, bf1, arc1),
    ] {
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
            out: -n,
        });
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
            additive: true,
            edges: vec![edges[link].0.clone()],
            melt_refusal: None,
            credit: Vec::new(),
        },
        [e0, e1],
    )))
}

/// The end of a pipe at `crossing`, the link in `slot`: the ball's centre
/// off the neighbour's wall at the seam's crossing, the section through
/// it square to the neighbour's run, and where it cuts the neighbour's
/// crease. `None` where a leg of that section does not lie on its face.
#[allow(clippy::too_many_arguments, reason = "one end, read whole")]
#[allow(clippy::too_many_lines, reason = "one end, read whole")]
fn end(
    model: &Model,
    faces: &Faces,
    edges: &[(Shape, f64)],
    crossing: &Crossing,
    slot: usize,
    cap_host: &(SurfaceGeometry, f64),
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<End>> {
    let other = 1 - slot;
    let neighbour = crossing.pair[other];
    let wall = crossing.walls[other];
    let Some(wall_host) = crate::engine::rolled::host(model, &faces.faces[wall], tol)? else {
        return Ok(None);
    };
    // The crossing again, off the neighbour's wall alone: its ball stands
    // the radius from the cap and from the seam's point.
    let (seam_curve, (t_v, t_q)) = &crossing.seam;
    let off = |t: f64| -> OgeomResult<(Point, Point, Point, Vector, f64)> {
        let touch = seam_curve.point_at(t, tol)?;
        let (_, n) = outward_at(&wall_host, touch, tol)?;
        let centre = touch + n * radius;
        let (on_cap, up) = outward_at(cap_host, centre, tol)?;
        Ok((
            touch,
            centre,
            on_cap,
            up,
            (centre - on_cap).dot(up) - radius,
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
    if ((centre - on_cap).dot(up) - radius).abs() > crate::marched::band_fit_target(tol) {
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
    // From the neighbour into the pipe: back towards the vertex.
    if (crossing.at - centre).dot(into) < 0.0 {
        into = -into;
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
    // The section's legs: straight from the crease to each contact, on the
    // cap and on the neighbour's wall.
    for (a, b, on) in [(foot, on_cap, cap_host), (foot, touch, &wall_host)] {
        for f in [0.25, 0.5, 0.75] {
            let p = a + (b - a) * f;
            let dd = ogeom_geom::project_on_surface(&on.0, p, 16, tol)?.distance;
            if dd > near {
                return Ok(None);
            }
        }
    }
    Ok(Some(End {
        neighbour,
        centre,
        into,
        foot,
        on_cap,
        touch,
        at: crossing.at,
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

/// The ball's centre nearest `guess` at `radius` off the cap and the
/// wall, with its contacts on each; `None` where it does not settle.
fn settled(
    cap: &(SurfaceGeometry, f64),
    wall: &(SurfaceGeometry, f64),
    guess: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point, Point)>> {
    let mut c = guess;
    for _ in 0..40 {
        let (on_cap, up) = outward_at(cap, c, tol)?;
        c = on_cap + up * radius;
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
/// the cap's and the wall's offsets; `None` where the march loses them.
fn centres(
    cap: &(SurfaceGeometry, f64),
    wall: &(SurfaceGeometry, f64),
    e0: &End,
    e1: &End,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Point>>> {
    let step = radius / 200.0;
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
        let Some((next, _, _)) = settled(cap, wall, c + t * step, radius, tol)? else {
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
