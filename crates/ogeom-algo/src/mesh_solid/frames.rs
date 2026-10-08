//! Frames for curved regions, chosen so every ring of a region's boundary
//! reads in one piece in its chart: border loops walked from a region's
//! half-edges, frames for closed surfaces the boundary only makes holes in
//! and for sphere caps, bands cut open along a slit, and coaxial surfaces
//! put on one axis and angular origin.

use ogeom_core::{FastMap, Tolerances};
use ogeom_math::{Direction, Frame, Point, Sphere, Vector};

use super::seams::{free_angle, loop_turns, widest_gap_holding};
use super::segment::{angular_spread, axis_frame, chart, on_frame, periodic};
use super::snap::plane_through;
use super::weld::Adjacency;
use super::{Carrier, Curved, Groups};
use crate::recognize::{Canonical, worst_deviation};

/// A region's boundary, walked into loops of mesh vertices: each border
/// half-edge leads to the one leaving its end. `None` where a vertex has
/// more than one way on.
pub(super) fn border_loops(
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    of: &[usize],
    g: usize,
) -> Option<Vec<Vec<u32>>> {
    let mut leaving: FastMap<u32, Vec<u32>> = FastMap::default();
    for (t, tri) in triangles.iter().enumerate() {
        if of[t] != g {
            continue;
        }
        for k in 0..3 {
            let inside = adjacency.twin[3 * t + k].is_some_and(|o| of[o / 3] == g);
            if !inside {
                leaving.entry(tri[k]).or_default().push(tri[(k + 1) % 3]);
            }
        }
    }
    if leaving.is_empty() || leaving.values().any(|to| to.len() != 1) {
        return None;
    }
    let mut loops: Vec<Vec<u32>> = Vec::new();
    let mut done: ogeom_core::FastSet<u32> = ogeom_core::FastSet::default();
    let mut starts: Vec<u32> = leaving.keys().copied().collect();
    starts.sort_unstable();
    for start in starts {
        if !done.insert(start) {
            continue;
        }
        let mut ring = vec![start];
        let mut at = leaving[&start][0];
        while at != start && ring.len() <= leaving.len() {
            done.insert(at);
            ring.push(at);
            at = leaving.get(&at).map_or(start, |to| to[0]);
        }
        loops.push(ring);
    }
    Some(loops)
}

/// A frame for each closed surface its boundary only makes holes in: a
/// sphere round its axis whose rings are not the parallels of one axis, or
/// a torus round both ways. Its seams are placed clear of every ring, the
/// sphere's poles as far from them as any axis puts them, so the whole
/// surface's face can carry the rings as holes. A torus whose rings go
/// round it one way is a band instead, and goes round the other way no
/// more. A sphere bounded by a rim or two as well as holes gets its frame
/// from [`cap_frame`].
pub(super) fn hole_frames(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    tol: Tolerances,
) {
    for g in 0..groups.carriers.len() {
        hole_frame(points, triangles, adjacency, groups, g, tol);
        cap_frame(points, triangles, adjacency, groups, g, tol);
    }
}

/// Whether the triangle whose centroid stands nearest `p` is in the region
/// `g`.
fn on_region(points: &[Point], triangles: &[[u32; 3]], of: &[usize], g: usize, p: Point) -> bool {
    triangles
        .iter()
        .enumerate()
        .map(|(t, tri)| {
            let c = Point::from_vector(
                (points[tri[0] as usize].to_vector()
                    + points[tri[1] as usize].to_vector()
                    + points[tri[2] as usize].to_vector())
                    / 3.0,
            );
            (c.distance(p), t)
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .is_some_and(|(_, t)| of[t] == g)
}

/// A frame for a sphere bounded by a rim and holes, or by two rims and
/// holes, that [`hole_frame`] left without one: an axis one ring goes
/// round (or two do) and every other ring does not. With one rim the face
/// is a cap about the pole on it, clear of every ring; with two, a zone
/// whose poles stand inside the rims, clear of them. The axes tried first
/// are the region's own and each ring's plane normal and mean direction,
/// so a rim in a plane is a latitude where that leaves the pole clear;
/// then the spread axes, the pole farthest from every ring.
fn cap_frame(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    g: usize,
    tol: Tolerances,
) {
    let Carrier::Curved(curved) = &groups.carriers[g] else {
        return;
    };
    let Canonical::Sphere(sphere) = curved.shape else {
        return;
    };
    if curved.fixed {
        return;
    }
    let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
        return;
    };
    if loops.len() < 2 {
        return;
    }
    let unit = |p: Point| {
        let d = p - sphere.centre();
        let m = d.magnitude();
        (m > 0.0).then(|| d / m)
    };
    let rings: Vec<Vec<Vector>> = loops
        .iter()
        .map(|ring| {
            ring.iter()
                .filter_map(|&v| unit(points[v as usize]))
                .collect()
        })
        .collect();
    let own = sphere.frame().z().vector();
    let mut preferred: Vec<Vector> = vec![own, -own];
    for (ring, directions) in loops.iter().zip(&rings) {
        let pts: Vec<Point> = ring.iter().map(|&v| points[v as usize]).collect();
        if let Some((_, n)) = (pts.len() >= 3).then(|| plane_through(&pts, tol)).flatten() {
            // A normal all but the region's own axis is that axis: the
            // axis may have been shared with the surfaces round it.
            let n = if n.vector().cross(own).magnitude() <= 1e-3 {
                own * n.vector().dot(own).signum()
            } else {
                n.vector()
            };
            preferred.extend([n, -n]);
        }
        let mean = directions.iter().fold(Vector::ZERO, |a, d| a + *d);
        let m = mean.magnitude();
        #[allow(clippy::cast_precision_loss, reason = "ring lengths are small")]
        if m > directions.len() as f64 * 0.1 {
            preferred.extend([mean / m, -mean / m]);
        }
    }
    // How near a ring the pole the face keeps comes (between two rims,
    // either pole), as a cosine; `None` for an axis that makes neither a
    // cap nor a zone, or a cap whose pole is off the region.
    let near = |z: Vector| -> Option<f64> {
        let turns: Vec<i32> = rings.iter().map(|ring| turns_about(ring, z)).collect();
        if turns.iter().any(|t| t.abs() > 1) {
            return None;
        }
        match turns.iter().filter(|t| **t != 0).count() {
            1 => Some(
                rings
                    .iter()
                    .flatten()
                    .map(|d| d.dot(z))
                    .fold(-1.0, f64::max),
            ),
            2 => Some(
                rings
                    .iter()
                    .flatten()
                    .map(|d| d.dot(z).abs())
                    .fold(0.0, f64::max),
            ),
            _ => None,
        }
    };
    let rims = |z: Vector| rings.iter().filter(|r| turns_about(r, z) != 0).count();
    let pole_on = |z: Vector| {
        rims(z) == 2
            || on_region(
                points,
                triangles,
                &groups.of,
                g,
                sphere.centre() + z * sphere.radius(),
            )
    };
    let clear = POLE_CLEARANCE.cos();
    let fits = |z: &Vector| near(*z).is_some_and(|n| n <= clear) && pole_on(*z);
    let chosen = preferred.iter().copied().find(fits).or_else(|| {
        let mut ranked: Vec<(f64, Vector)> = spread_directions(POLE_CANDIDATES)
            .into_iter()
            .filter_map(|z| near(z).filter(|n| *n <= clear).map(|n| (n, z)))
            .collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
        ranked.into_iter().map(|(_, z)| z).find(|z| pole_on(*z))
    });
    let Some(z) = chosen else {
        return;
    };
    let Ok(z) = Direction::new(z, tol) else {
        return;
    };
    let x = if z.vector().cross(own).magnitude() <= 1e-12 {
        sphere.frame().x()
    } else {
        z.any_perpendicular()
    };
    let Ok(frame) = Frame::new(sphere.centre(), z, x, tol) else {
        return;
    };
    let Ok(turned) = Sphere::new(frame, sphere.radius(), tol) else {
        return;
    };
    if let Carrier::Curved(curved) = &mut groups.carriers[g] {
        curved.shape = Canonical::Sphere(turned);
        curved.fixed = true;
        curved.wraps = true;
        curved.wraps_v = false;
        curved.centre = (core::f64::consts::PI, curved.centre.1);
    }
}

/// [`hole_frames`] for the region `g`.
pub(super) fn hole_frame(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    g: usize,
    tol: Tolerances,
) {
    let Carrier::Curved(curved) = &groups.carriers[g] else {
        return;
    };
    let wanted = match curved.shape {
        Canonical::Sphere(_) => curved.wraps && !curved.fixed,
        Canonical::Torus(_) => curved.wraps && curved.wraps_v,
        _ => false,
    };
    if !wanted {
        return;
    }
    let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
        return;
    };
    // A torus whose boundary goes round one way only is a band, however
    // narrow the strip its rims leave out: it goes round the other way
    // no more, and its chart is centred on its own middle across the rims.
    let turns: Option<Vec<(i32, i32)>> = loops
        .iter()
        .map(|ring| loop_turns(&curved.shape, ring, points, tol))
        .collect();
    if matches!(curved.shape, Canonical::Torus(_))
        && let Some(turns) = turns
    {
        let round_axis = turns.iter().any(|t| t.0 != 0);
        let round_tube = turns.iter().any(|t| t.1 != 0);
        if round_axis != round_tube {
            let across =
                |p: Point| chart(&curved.shape, p, tol).map(|c| if round_axis { c.1 } else { c.0 });
            let on_rims: ogeom_core::FastSet<u32> = loops.iter().flatten().copied().collect();
            let mut rims: Vec<f64> = on_rims
                .iter()
                .filter_map(|&v| across(points[v as usize]))
                .collect();
            let inside: Vec<f64> = curved
                .vertices
                .iter()
                .filter(|v| !on_rims.contains(*v))
                .filter_map(|&v| across(points[v as usize]))
                .collect();
            let Some(middle) = widest_gap_holding(&mut rims, &inside) else {
                return;
            };
            if let Carrier::Curved(curved) = &mut groups.carriers[g] {
                if round_axis {
                    curved.wraps_v = false;
                    curved.centre.1 = middle;
                } else {
                    curved.wraps = false;
                    curved.centre.0 = middle;
                }
            }
            return;
        }
    }
    let ring_points: Vec<Point> = loops
        .iter()
        .flatten()
        .map(|&v| points[v as usize])
        .collect();
    let shape = match curved.shape.clone() {
        Canonical::Sphere(sphere) => {
            // The axis whose poles stand farthest from every ring point,
            // of those no ring goes round: a pole inside a hole is off
            // the face, however far it stands from the hole's edge.
            let unit = |p: Point| {
                let d = p - sphere.centre();
                let m = d.magnitude();
                (m > 0.0).then(|| d / m)
            };
            let directions: Vec<Vector> = ring_points.iter().filter_map(|p| unit(*p)).collect();
            let rings: Vec<Vec<Vector>> = loops
                .iter()
                .map(|ring| {
                    ring.iter()
                        .filter_map(|&v| unit(points[v as usize]))
                        .collect()
                })
                .collect();
            let mut ranked: Vec<(f64, Vector)> = spread_directions(POLE_CANDIDATES)
                .into_iter()
                .map(|z| {
                    let nearest = directions
                        .iter()
                        .map(|d| d.dot(z).abs())
                        .fold(0.0_f64, f64::max);
                    (nearest, z)
                })
                .collect();
            ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
            // And both poles on the region itself: a loop no axis goes
            // round has both poles to one side of it, which may be the
            // hole's.
            let on_region = |p: Point| on_region(points, triangles, &groups.of, g, p);
            let best = ranked.into_iter().find(|(_, z)| {
                rings.iter().all(|ring| turns_about(ring, *z) == 0)
                    && on_region(sphere.centre() + *z * sphere.radius())
                    && on_region(sphere.centre() - *z * sphere.radius())
            });
            let Some((nearest, z)) = best else {
                return;
            };
            // A pole within a few degrees of a ring has no room round it.
            if nearest > POLE_CLEARANCE.cos() {
                return;
            }
            let Ok(z) = Direction::new(z, tol) else {
                return;
            };
            let Ok(frame) = Frame::new(sphere.centre(), z, z.any_perpendicular(), tol) else {
                return;
            };
            let Ok(turned) = Sphere::new(frame, sphere.radius(), tol) else {
                return;
            };
            Canonical::Sphere(turned)
        }
        other => other,
    };
    // Then the seam, turned about the axis into the widest angle the
    // rings leave free.
    let angles: Vec<Vec<(f64, f64)>> = vec![
        ring_points
            .iter()
            .filter_map(|p| chart(&shape, *p, tol))
            .collect(),
    ];
    let Some(free) = free_angle(&angles) else {
        return;
    };
    let Some(frame) = axis_frame(&shape) else {
        return;
    };
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let Ok(x) = Direction::new(x * free.cos() + y * free.sin(), tol) else {
        return;
    };
    let Ok(turned) = Frame::new(frame.origin(), frame.z(), x, tol) else {
        return;
    };
    let Some(shape) = on_frame(&shape, turned, tol) else {
        return;
    };
    // A torus's seam round its axis is a parallel, placed at the tube
    // angle the rings leave widest free; its chart is centred half a
    // turn on from it.
    let centre_v = match shape {
        Canonical::Torus(_) => {
            let across: Vec<Vec<(f64, f64)>> =
                vec![angles[0].iter().map(|&(u, v)| (v, u)).collect()];
            let Some(free_v) = free_angle(&across) else {
                return;
            };
            free_v + core::f64::consts::PI
        }
        _ => curved.centre.1,
    };
    if let Carrier::Curved(curved) = &mut groups.carriers[g] {
        curved.shape = shape;
        curved.fixed = true;
        curved.centre = (core::f64::consts::PI, centre_v);
    }
}

/// How many times a loop of directions goes round an axis.
fn turns_about(ring: &[Vector], z: Vector) -> i32 {
    let x = if z.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = x - z * x.dot(z);
    let y = z.cross(x);
    let angle = |d: &Vector| d.dot(y).atan2(d.dot(x));
    let mut turned = 0.0;
    for (i, d) in ring.iter().enumerate() {
        let next = &ring[(i + 1) % ring.len()];
        turned += ogeom_math::elementary::wrap_signed_angle(angle(next) - angle(d));
    }
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    let turns = (turned / core::f64::consts::TAU).round() as i32;
    turns
}

/// A region round its axis whose boundary is one loop going round it no
/// times is a band cut open along a slit (a few triangles its growth left
/// across it): no seam can join rims it does not have. It is built as a
/// patch instead, its chart cut placed in the slit, the widest gap in the
/// angles its vertices stand at.
pub(super) fn slit_bands(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    tol: Tolerances,
) {
    for g in 0..groups.carriers.len() {
        slit_band(points, triangles, adjacency, groups, g, tol);
    }
}

/// [`slit_bands`] for the region `g`.
pub(super) fn slit_band(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    g: usize,
    tol: Tolerances,
) {
    let Carrier::Curved(curved) = &groups.carriers[g] else {
        return;
    };
    if !curved.wraps
        || curved.wraps_v
        || !matches!(curved.shape, Canonical::Cylinder(_) | Canonical::Cone(_))
    {
        return;
    }
    let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
        return;
    };
    let [ring] = &loops[..] else {
        return;
    };
    let Some(frame) = axis_frame(&curved.shape) else {
        return;
    };
    let directions: Vec<Vector> = ring
        .iter()
        .map(|&v| points[v as usize] - frame.origin())
        .collect();
    if turns_about(&directions, frame.z().vector()) != 0 {
        return;
    }
    let mut angles: Vec<f64> = curved
        .vertices
        .iter()
        .filter_map(|&v| chart(&curved.shape, points[v as usize], tol).map(|c| c.0))
        .collect();
    let Some(gap) = widest_gap(&mut angles) else {
        return;
    };
    if let Carrier::Curved(curved) = &mut groups.carriers[g] {
        curved.wraps = false;
        curved.centre = (
            ogeom_math::elementary::wrap_angle(gap + core::f64::consts::PI),
            curved.centre.1,
        );
    }
}

/// The middle of the widest gap between angles.
fn widest_gap(angles: &mut [f64]) -> Option<f64> {
    let tau = core::f64::consts::TAU;
    if angles.is_empty() {
        return None;
    }
    for a in angles.iter_mut() {
        *a = a.rem_euclid(tau);
    }
    angles.sort_by(f64::total_cmp);
    let mut best = (
        angles[0] + tau - angles[angles.len() - 1],
        angles[angles.len() - 1],
    );
    for pair in angles.windows(2) {
        if pair[1] - pair[0] > best.0 {
            best = (pair[1] - pair[0], pair[0]);
        }
    }
    Some(best.1 + best.0 / 2.0)
}

/// How many axes a sphere's poles are tried along.
const POLE_CANDIDATES: usize = 400;

/// How near a ring a sphere's pole may stand: five degrees.
const POLE_CLEARANCE: f64 = 0.087;

/// Directions spread evenly over the sphere: a Fibonacci lattice.
fn spread_directions(count: usize) -> Vec<Vector> {
    let golden = core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            #[allow(clippy::cast_precision_loss, reason = "a few hundred directions")]
            let (i, n) = (i as f64, count as f64);
            let z = 1.0 - 2.0 * (i + 0.5) / n;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let a = golden * i;
            Vector::new(r * a.cos(), r * a.sin(), z)
        })
        .collect()
}

/// A curved region put on `frame`, whose axis is nearly its own: its
/// surface carried there if it still holds the region's vertices within
/// `flat`, or a cylinder or cone fitted again about the new axis if that
/// does. Otherwise the region keeps its own fit.
fn onto_axis(curved: &mut Curved, points: &[Point], frame: Frame, flat: f64, tol: Tolerances) {
    let Some(shape) = on_frame(&curved.shape, frame, tol) else {
        return;
    };
    let pts: Vec<Point> = curved
        .vertices
        .iter()
        .map(|&v| points[v as usize])
        .collect();
    let deviation = worst_deviation(&shape, &pts);
    if deviation <= flat {
        curved.shape = shape;
        curved.deviation = deviation;
    } else if matches!(curved.shape, Canonical::Cylinder(_) | Canonical::Cone(_))
        && let Some(refitted) =
            crate::recognize::ruled_about(&pts, frame.origin(), frame.z(), flat, tol)
        && let Some(refitted) = on_frame(&refitted, frame, tol)
    {
        // The fit's own axis stood a few slops off this one, and its
        // radius and lean with it: fitted again about this axis, it meets
        // its neighbours on their circles.
        let deviation = worst_deviation(&refitted, &pts);
        if deviation <= flat {
            curved.shape = refitted;
            curved.deviation = deviation;
        }
    }
}

/// Put coaxial surfaces on one axis and one angular origin, so bands that
/// meet along a circle meet at their seams too, and an axis all but square
/// to a plane square to it; then fix each curved region's chart branch and
/// whether it wraps.
pub(super) fn align_axes(
    points: &[Point],
    groups: &mut Groups,
    normals: &[Direction],
    flat: f64,
    tol: Tolerances,
) {
    let mut leaders: Vec<Frame> = Vec::new();
    // The best determined axis leads: a cylinder's before a cone's, whose
    // lean trades against its axis on a short band, and a larger region
    // before a smaller.
    let rank = |c: &Carrier| match c {
        Carrier::Curved(curved) => (
            match curved.shape {
                Canonical::Cylinder(_) => 0,
                Canonical::Cone(_) => 1,
                _ => 2,
            },
            usize::MAX - curved.vertices.len(),
        ),
        _ => (3, 0),
    };
    let mut order: Vec<usize> = (0..groups.carriers.len()).collect();
    order.sort_by_key(|&i| rank(&groups.carriers[i]));
    for i in order {
        let Carrier::Curved(curved) = &mut groups.carriers[i] else {
            continue;
        };
        align_one(points, curved, &mut leaders, normals, flat, tol);
    }
}

/// [`align_axes`] for one curved region, after the regions `leaders`
/// were taken from.
pub(super) fn align_one(
    points: &[Point],
    curved: &mut Curved,
    leaders: &mut Vec<Frame>,
    normals: &[Direction],
    flat: f64,
    tol: Tolerances,
) {
    let Some(frame) = axis_frame(&curved.shape) else {
        return;
    };
    let sphere = matches!(curved.shape, Canonical::Sphere(_));
    // An axis all but square to a plane of the solid is square to it:
    // the mesh's slop leans the fit by a few millionths, and a leaning
    // axis meets the plane in an ellipse where the part has a circle.
    let frame = match normals.iter().find(|n| {
        let lean = n.vector().cross(frame.z().vector()).magnitude();
        lean > 0.0 && lean <= 1e-3
    }) {
        Some(&normal) if !sphere => {
            let axis = if normal.vector().dot(frame.z().vector()) >= 0.0 {
                normal
            } else {
                -normal
            };
            if let Ok(square) = Frame::new(frame.origin(), axis, frame.x(), tol) {
                onto_axis(curved, points, square, flat, tol);
            }
            axis_frame(&curved.shape).unwrap_or(frame)
        }
        _ => frame,
    };
    let lead = leaders.iter().find(|l| {
        let parallel = l.z().vector().cross(frame.z().vector()).magnitude() <= 1e-3;
        let w = frame.origin() - l.origin();
        let off = (w - l.z().vector() * w.dot(l.z().vector())).magnitude();
        // A sphere centred on the axis takes its frame too: any frame
        // through its centre is exact, and seams meeting on the circle
        // it shares with the axis's other surfaces must start from one
        // angle.
        parallel && off <= flat * 10.0
    });
    if let Some(lead) = lead {
        let z = lead.z().vector();
        let w = frame.origin() - lead.origin();
        let origin = lead.origin() + z * w.dot(z);
        let axis = if frame.z().vector().dot(z) >= 0.0 {
            lead.z()
        } else {
            -lead.z()
        };
        if let Ok(snapped) = Frame::new(origin, axis, lead.x(), tol) {
            onto_axis(curved, points, snapped, flat, tol);
        }
    } else if !sphere {
        leaders.push(frame);
    }
    // The branch: the region's mean angle, and whether it wraps.
    let charts: Vec<(f64, f64)> = curved
        .vertices
        .iter()
        .filter_map(|&v| chart(&curved.shape, points[v as usize], tol))
        .collect();
    let mut us: Vec<f64> = charts.iter().map(|c| c.0).collect();
    let (_, gap_u) = angular_spread(&mut us);
    let (_, pv) = periodic(&curved.shape);
    if pv {
        let mut vs: Vec<f64> = charts.iter().map(|c| c.1).collect();
        curved.wraps_v = angular_spread(&mut vs).1 < core::f64::consts::FRAC_PI_2;
    }
    let wraps_u = gap_u < core::f64::consts::FRAC_PI_2;
    curved.wraps = wraps_u;
    if !curved.wraps
        && !curved.wraps_v
        && !curved.fixed
        && let Some(reframed) = away_from(curved, points, tol)
    {
        curved.shape = reframed;
    }
    // The branch every pcurve is read on: the region's own mean chart
    // point, which after the re-framing sits half a turn from the cut.
    let charts: Vec<(f64, f64)> = curved
        .vertices
        .iter()
        .filter_map(|&v| chart(&curved.shape, points[v as usize], tol))
        .collect();
    let mut us: Vec<f64> = charts.iter().map(|c| c.0).collect();
    let (mean_u, _) = angular_spread(&mut us);
    let mean_v = if curved.wraps_v {
        core::f64::consts::PI
    } else if pv {
        let mut vs: Vec<f64> = charts.iter().map(|c| c.1).collect();
        angular_spread(&mut vs).0
    } else {
        #[allow(
            clippy::cast_precision_loss,
            reason = "vertex counts are far below 2^52"
        )]
        let count = charts.len().max(1) as f64;
        charts.iter().map(|c| c.1).sum::<f64>() / count
    };
    let centre_u = if wraps_u {
        core::f64::consts::PI
    } else {
        ogeom_math::elementary::wrap_angle(mean_u)
    };
    curved.centre = (centre_u, mean_v);
}

/// The surface on a frame that puts the region half a turn from its
/// chart's cut (and a sphere's poles a quarter turn to either side of it),
/// so every image of its boundary reads in one piece.
fn away_from(curved: &Curved, points: &[Point], tol: Tolerances) -> Option<Canonical> {
    let mut mean = Vector::ZERO;
    match curved.shape {
        Canonical::Sphere(s) => {
            for &v in &curved.vertices {
                let d = points[v as usize] - s.centre();
                let m = d.magnitude();
                if m > 0.0 {
                    mean += d / m;
                }
            }
            let facing = Direction::new(mean, tol).ok()?;
            let frame = Frame::new(s.centre(), facing.any_perpendicular(), -facing, tol).ok()?;
            Some(Canonical::Sphere(Sphere::new(frame, s.radius(), tol).ok()?))
        }
        _ => {
            let frame = axis_frame(&curved.shape)?;
            let z = frame.z().vector();
            for &v in &curved.vertices {
                let w = points[v as usize] - frame.origin();
                let r = w - z * w.dot(z);
                let m = r.magnitude();
                if m > 0.0 {
                    mean += r / m;
                }
            }
            let facing = Direction::new(mean, tol).ok()?;
            on_frame(
                &curved.shape,
                Frame::new(frame.origin(), frame.z(), -facing, tol).ok()?,
                tol,
            )
        }
    }
}
