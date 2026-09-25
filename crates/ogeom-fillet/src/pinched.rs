//! The pinched fillet: a seat whose two hosts turn tangent at an end of
//! the crease, so the ball's section shrinks to nothing there.
//!
//! Two equal drums crossing square are the model case. Their seam is two
//! ellipses meeting at the two points where the drums touch. Along the
//! seam the drums' normals open an angle that falls to zero at each such
//! pole, and the ball's arc, whose sweep is that angle, falls with it: the
//! band pinches to a point there.
//!
//! The walker cannot march into a pole (the seat's equations turn singular
//! as the two touch points merge), so the stations are solved one by one at
//! guide parameters crowded toward the pole, each warm-started from its
//! neighbour, until the shrinking section stops converging. The pole's own
//! section closes the run exactly: both touch points on the pole, the
//! ball's centre a radius off it along the normal the hosts share. The band
//! and wedge are then the open band's, with no cap at a pinched end.

use crate::march::{BlendStop, MarchedBlend, Sides, seat_section};
use crate::support::edge_curve;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

/// Stations across the edge, crowded toward its ends by a cosine spacing,
/// where a pinched section changes fastest relative to its size.
const STATIONS: usize = 96;

/// Whether the two hosts are tangent at `p`: their unit normals parallel to
/// within what a seat solve can separate.
pub(crate) fn hosts_tangent_at(
    first: &SurfaceGeometry,
    second: &SurfaceGeometry,
    p: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (n1, _) = unit_normal_near(first, p, tol)?;
    let (n2, _) = unit_normal_near(second, p, tol)?;
    Ok(n1.cross(n2).magnitude() <= 1e-6)
}

fn unit_normal_near(
    surface: &SurfaceGeometry,
    p: Point,
    tol: Tolerances,
) -> OgeomResult<(Vector, (f64, f64))> {
    let (u, v) = ogeom_algo::project_on_surface(surface, p, 32, tol)?.parameters;
    Ok((unit_normal(surface, u, v, tol)?, (u, v)))
}

fn unit_normal(surface: &SurfaceGeometry, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Vector> {
    let (du, dv) = surface.d1_at(u, v, tol)?;
    let n = du.cross(dv);
    let m = n.magnitude();
    if m <= tol.angular() {
        ogeom_bail!(Construction, "a host has no normal at the seat");
    }
    Ok(n / m)
}

/// One station: guide parameter, contact parameters, centre, touches.
type Station = (f64, [f64; 4], Point, Point, Point);

/// The fillet of a crease pinched at the ends marked in `pinched`, built as
/// the open band over the edge's own window.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
pub(crate) fn pinched_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    hosts: [(&SurfaceGeometry, f64); 2],
    sides: Sides,
    convex: bool,
    radius: f64,
    pinched: [bool; 2],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (guide, (lo, hi)) = edge_curve(model, edge, tol)?;
    let [(first, _), (second, _)] = hosts;

    // One solved section, or `None` where the seat does not converge to
    // two distinct touch points.
    let solve_at = |t: f64, near: [f64; 4]| -> Option<Station> {
        let x = seat_section(first, second, radius, &guide, sides, t, near, tol).ok()?;
        let p1 = first.point_at(x[0], x[1], tol).ok()?;
        let p2 = second.point_at(x[2], x[3], tol).ok()?;
        let n1 = unit_normal(first, x[0], x[1], tol).ok()?;
        let n2 = unit_normal(second, x[2], x[3], tol).ok()?;
        let c1 = p1 + n1 * (f64::from(sides.first) * radius);
        let c2 = p2 + n2 * (f64::from(sides.second) * radius);
        if c1.distance(c2) > tol.confusion() * 10.0 || p1.distance(p2) <= tol.confusion() * 10.0 {
            return None;
        }
        Some((t, [x[0], x[1], x[2], x[3]], c1.midpoint(c2), p1, p2))
    };
    #[allow(clippy::cast_precision_loss)]
    let parameter = |i: usize| -> f64 {
        let s = (1.0 - (core::f64::consts::PI * (i as f64) / (STATIONS as f64)).cos()) * 0.5;
        lo + (hi - lo) * s
    };

    // Solved from the middle out, so each Newton starts in the basin of
    // this seating as the section shrinks toward a pole.
    let middle = STATIONS / 2;
    let seed = {
        let p = guide.point_at(parameter(middle), tol)?;
        let a = ogeom_algo::project_on_surface(first, p, 32, tol)?.parameters;
        let b = ogeom_algo::project_on_surface(second, p, 32, tol)?.parameters;
        [a.0, a.1, b.0, b.1]
    };
    let Some(centre) = solve_at(parameter(middle), seed) else {
        ogeom_bail!(
            NotDone,
            "no ball of radius {radius} seats at the middle of the pinched crease"
        );
    };
    let mut run: [Vec<Station>; 2] = [Vec::new(), Vec::new()];
    for (side, indices) in [
        (0, (0..middle).rev().collect::<Vec<_>>()),
        (1, (middle + 1..=STATIONS).collect::<Vec<_>>()),
    ] {
        let mut near = centre.1;
        for i in indices {
            // A pinched end's own station is the pole's, placed below.
            if (i == 0 || i == STATIONS) && pinched[side] {
                break;
            }
            match solve_at(parameter(i), near) {
                Some(found) => {
                    near = found.1;
                    run[side].push(found);
                }
                None if pinched[side] => break,
                None => ogeom_bail!(
                    NotDone,
                    "the ball does not seat at the unpinched end of the crease"
                ),
            }
        }
    }
    let [mut stations, toward_hi] = run;
    stations.reverse();
    stations.push(centre);
    stations.extend(toward_hi);
    // Between a pinched run's last station and its pole the band is only
    // fitted; that gap is held to a small share of the crease.
    let gap = (stations[0].0 - lo).max(hi - stations[stations.len() - 1].0);
    if gap > (hi - lo) * 0.02 {
        ogeom_bail!(
            NotDone,
            "the ball stops seating {gap} short of a pole of the pinched crease"
        );
    }

    let pole = |t: f64| -> OgeomResult<Station> {
        let at = guide.point_at(t, tol)?;
        let (n, a) = unit_normal_near(first, at, tol)?;
        let (_, b) = unit_normal_near(second, at, tol)?;
        let centre = at + n * (f64::from(sides.first) * radius);
        Ok((t, [a.0, a.1, b.0, b.1], centre, at, at))
    };
    if pinched[0] {
        stations.insert(0, pole(lo)?);
    }
    if pinched[1] {
        stations.push(pole(hi)?);
    }

    let mut blend = MarchedBlend {
        spine: Vec::with_capacity(stations.len()),
        on_first: Vec::with_capacity(stations.len()),
        on_second: Vec::with_capacity(stations.len()),
        touch_first: Vec::with_capacity(stations.len()),
        touch_second: Vec::with_capacity(stations.len()),
        along: Vec::with_capacity(stations.len()),
        sides,
        stopped: BlendStop::RanPastTheGuide,
    };
    for (t, x, centre, p1, p2) in stations {
        blend.along.push(t);
        blend.on_first.push((x[0], x[1]));
        blend.on_second.push((x[2], x[3]));
        blend.spine.push(centre);
        blend.touch_first.push(p1);
        blend.touch_second.push(p2);
    }
    crate::marched::build_open_band(
        model, solid, edge, &blend, &guide, hosts, radius, convex, pinched, tol,
    )
}
