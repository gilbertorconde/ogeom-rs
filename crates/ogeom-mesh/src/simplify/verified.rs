//! Decimation checked by measuring, both ways, and repeated until it passes.
//!
//! The quadric cost that orders collapses is a mean over planes, and a mean
//! can be small where one point has moved far. So the result is measured
//! against the source rather than trusted: each source vertex to the result,
//! and each result vertex and triangle centre to the source. Where a point is
//! too far, the source vertices responsible are pinned and the decimation runs
//! again from the source.
//!
//! # Why it stops
//!
//! A pinned vertex never moves and nothing merges into it, so a result vertex
//! made only of pinned source vertices is one of them, where it was, and a
//! triangle whose corners are all pinned is a source triangle: neither is
//! farther from the source than rounding. A source vertex that is pinned and
//! still far from the result has lost every triangle around it, which needs an
//! edge between two unpinned neighbours, so its neighbours are pinned in turn. Every point
//! that fails by more than rounding therefore names a source vertex not yet
//! pinned, each round pins at least one more, and with every vertex pinned the
//! result is the source. A round that pins nothing new ends the loop, the
//! shortfall reported.

use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::Triangulation;

use super::{Deviation, Simplified, State};
use crate::Bvh;

/// [`super::simplify`] for [`super::Target::Verified`], on a mesh whose
/// triangles name only vertices it has.
pub(super) fn simplify(
    mesh: &Triangulation,
    tolerance: f64,
    share: f64,
    max_triangles: Option<usize>,
    tol: Tolerances,
) -> OgeomResult<Simplified> {
    let source = Bvh::build(mesh)?;
    let mut used: Vec<u32> = mesh.triangles.iter().flatten().copied().collect();
    used.sort_unstable();
    used.dedup();

    let budget = tolerance * tolerance;
    let mut held: FastSet<u32> = FastSet::default();
    let mut rounds = 0;
    let (mut state, mut collapsed, mut measured) = loop {
        rounds += 1;
        let mut state = State::new(mesh, tol);
        state.held.clone_from(&held);
        let collapsed = state.decimate_queued(budget, 1);
        let measured = Measured::of(&source, mesh, &used, &state, share)?;
        if measured.distance_at_share <= tolerance {
            break (state, collapsed, measured);
        }
        let before = held.len();
        pin(mesh, &used, &state, &measured, tolerance, &mut held);
        if held.len() == before {
            // Every far point is made of pinned vertices only, so rounding,
            // not a collapse, put it past a tolerance below it.
            break (state, collapsed, measured);
        }
    };

    if let Some(most) = max_triangles
        && state.live_triangles() > most
    {
        state.held.clear();
        collapsed += state.decimate_queued(f64::INFINITY, most);
        measured = Measured::of(&source, mesh, &used, &state, share)?;
    }

    let target_met = measured.distance_at_share <= tolerance
        && max_triangles.is_none_or(|most| state.live_triangles() <= most);
    let within_share = measured.within(tolerance);
    Ok(Simplified {
        mesh: measured.mesh,
        error: measured.furthest.1,
        collapsed,
        target_met,
        deviation: Some(Deviation {
            within_share,
            distance_at_share: measured.distance_at_share,
            furthest: measured.furthest,
            rms: measured.rms,
            pinned: held.len(),
            rounds,
        }),
    })
}

/// A decimation's result and its distances from the source, both ways.
struct Measured {
    mesh: Triangulation,
    /// Each surviving vertex's number in [`Measured::mesh`].
    renumbered: FastMap<u32, u32>,
    /// Each used source vertex's distance to the result, in `used`'s order.
    source: Vec<f64>,
    /// Each result vertex's distance to the source, then each result
    /// triangle centre's.
    result: Vec<f64>,
    distance_at_share: f64,
    furthest: (Point, f64),
    rms: f64,
}

impl Measured {
    fn of(
        index: &Bvh,
        mesh: &Triangulation,
        used: &[u32],
        state: &State,
        share: f64,
    ) -> OgeomResult<Self> {
        let (result, renumbered) = state.harvest();
        let mut out = result.positions.clone();
        out.extend(result.triangles.iter().map(|t| centre(&result, *t)));
        let back = distances(index, &out);
        let sources: Vec<Point> = used.iter().map(|&s| mesh.positions[s as usize]).collect();
        let there = distances(&Bvh::build(&result)?, &sources);

        let mut furthest = (Point::ORIGIN, 0.0);
        let mut squares = 0.0;
        for (p, d) in sources.iter().zip(&there).chain(out.iter().zip(&back)) {
            if *d > furthest.1 {
                furthest = (*p, *d);
            }
            squares += d * d;
        }
        let count = there.len() + back.len();
        #[allow(clippy::cast_precision_loss)]
        let rms = if count == 0 {
            0.0
        } else {
            (squares / count as f64).sqrt()
        };
        Ok(Self {
            distance_at_share: at_share(&there, share).max(at_share(&back, share)),
            mesh: result,
            renumbered,
            source: there,
            result: back,
            furthest,
            rms,
        })
    }

    /// The smaller of the two sides' fractions within `tolerance`.
    fn within(&self, tolerance: f64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let fraction = |side: &[f64]| {
            if side.is_empty() {
                1.0
            } else {
                side.iter().filter(|d| **d <= tolerance).count() as f64 / side.len() as f64
            }
        };
        fraction(&self.source).min(fraction(&self.result))
    }
}

/// The distance from each point to the nearest point on the mesh.
fn distances(index: &Bvh, points: &[Point]) -> Vec<f64> {
    index
        .nearest_many(points, f64::INFINITY)
        .into_iter()
        .map(|hit| hit.map_or(f64::INFINITY, |hit| hit.distance))
        .collect()
}

/// The smallest distance within which at least `share` of the distances lie.
fn at_share(distances: &[f64], share: f64) -> f64 {
    if distances.is_empty() {
        return 0.0;
    }
    let mut sorted = distances.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let needed = ((share * sorted.len() as f64).ceil() as usize).clamp(1, sorted.len());
    sorted[needed - 1]
}

fn centre(mesh: &Triangulation, triangle: [u32; 3]) -> Point {
    let [a, b, c] = triangle.map(|v| mesh.positions[v as usize].to_vector());
    Point::from_vector((a + b + c) * (1.0 / 3.0))
}

/// Pin the source vertices behind every point farther than `tolerance`.
///
/// A far result vertex or triangle centre holds the source vertex nearest it
/// among those merged into it or into the triangle's corners; none of them is
/// pinned, or the point would not have moved. A far source vertex holds
/// itself, or, already held, its neighbours: it is far only because every
/// triangle around it was collapsed away.
fn pin(
    mesh: &Triangulation,
    used: &[u32],
    state: &State,
    measured: &Measured,
    tolerance: f64,
    held: &mut FastSet<u32>,
) {
    let vertices = measured.mesh.positions.len();
    let mut members: Vec<Vec<u32>> = vec![Vec::new(); vertices];
    for &s in used {
        if let Some(&r) = measured.renumbered.get(&state.survivor(s)) {
            members[r as usize].push(s);
        }
    }

    for (k, d) in measured.result.iter().enumerate() {
        if *d <= tolerance {
            continue;
        }
        let (point, corners) = if k < vertices {
            #[allow(clippy::cast_possible_truncation)]
            (measured.mesh.positions[k], vec![k as u32])
        } else {
            let triangle = measured.mesh.triangles[k - vertices];
            (centre(&measured.mesh, triangle), triangle.to_vec())
        };
        let nearest = corners
            .iter()
            .flat_map(|r| &members[*r as usize])
            .filter(|&&s| !state.pinned.contains(&s) && !state.held.contains(&s))
            .map(|&s| ((mesh.positions[s as usize] - point).magnitude(), s))
            .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        if let Some((_, s)) = nearest {
            held.insert(s);
        }
    }

    let mut stuck: FastSet<u32> = FastSet::default();
    for (&s, d) in used.iter().zip(&measured.source) {
        if *d <= tolerance {
            continue;
        }
        if state.pinned.contains(&s) || state.held.contains(&s) {
            stuck.insert(s);
        } else {
            held.insert(s);
        }
    }
    if !stuck.is_empty() {
        for triangle in &mesh.triangles {
            if triangle.iter().any(|v| stuck.contains(v)) {
                held.extend(triangle.iter().copied());
            }
        }
    }
    held.retain(|v| !state.pinned.contains(v));
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::{Target, simplify};
    use crate::{Deflection, triangulate};
    use ogeom_algo::make_sphere;
    use ogeom_core::Tolerances;
    use ogeom_math::{Frame, Point};
    use ogeom_topo::{Model, Triangulation};

    const T: Tolerances = Tolerances::millimetres();

    fn mesh_of(shape: impl FnOnce(&mut Model) -> ogeom_topo::Shape, chord: f64) -> Triangulation {
        let mut model = Model::new();
        let shape = shape(&mut model);
        triangulate(
            &model,
            &shape,
            Deflection {
                chord,
                ..Deflection::default()
            },
            T,
        )
        .unwrap()
    }

    /// A sphere of radius 10 with every vertex moved along its radius by a
    /// reproducible noise of about `sigma`: the sum of three uniform draws,
    /// near enough a normal distribution for a scan.
    fn noisy_sphere(chord: f64, sigma: f64) -> Triangulation {
        let mut mesh = mesh_of(
            |m| make_sphere(m, Frame::WORLD, 10.0, T).unwrap().shape,
            chord,
        );
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut uniform = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            #[allow(clippy::cast_precision_loss)]
            let unit = (state >> 11) as f64 / (1_u64 << 53) as f64;
            unit - 0.5
        };
        for p in &mut mesh.positions {
            // Three draws on [-1/2, 1/2] have variance 1/4, so twice their
            // sum has a standard deviation of one.
            let noise = 2.0 * (uniform() + uniform() + uniform()) * sigma;
            let radial = p.to_vector() * (1.0 / p.to_vector().magnitude());
            *p += radial * noise;
        }
        mesh
    }

    /// The distance from `p` to triangle `abc`, from the closest point found
    /// by the triangle's Voronoi regions, independently of the hierarchy.
    fn to_triangle(p: Point, [a, b, c]: [Point; 3]) -> f64 {
        let (ab, ac, ap) = (b - a, c - a, p - a);
        let (d1, d2) = (ab.dot(ap), ac.dot(ap));
        if d1 <= 0.0 && d2 <= 0.0 {
            return (p - a).magnitude();
        }
        let bp = p - b;
        let (d3, d4) = (ab.dot(bp), ac.dot(bp));
        if d3 >= 0.0 && d4 <= d3 {
            return (p - b).magnitude();
        }
        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return (p - (a + ab * v)).magnitude();
        }
        let cp = p - c;
        let (d5, d6) = (ab.dot(cp), ac.dot(cp));
        if d6 >= 0.0 && d5 <= d6 {
            return (p - c).magnitude();
        }
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return (p - (a + ac * w)).magnitude();
        }
        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return (p - (b + (c - b) * w)).magnitude();
        }
        let denom = 1.0 / (va + vb + vc);
        let (v, w) = (vb * denom, vc * denom);
        (p - (a + ab * v + ac * w)).magnitude()
    }

    fn brute(p: Point, mesh: &Triangulation) -> f64 {
        mesh.triangles
            .iter()
            .map(|t| to_triangle(p, t.map(|v| mesh.positions[v as usize])))
            .fold(f64::INFINITY, f64::min)
    }

    /// Both ways by brute force: the source vertices' distances to the
    /// result, and the result vertices' and centres' to the source.
    fn both_ways(source: &Triangulation, result: &Triangulation) -> (Vec<f64>, Vec<f64>) {
        let there = source.positions.iter().map(|p| brute(*p, result)).collect();
        let mut points = result.positions.clone();
        for t in &result.triangles {
            let [a, b, c] = t.map(|v| result.positions[v as usize].to_vector());
            points.push(Point::from_vector((a + b + c) * (1.0 / 3.0)));
        }
        let back = points.iter().map(|p| brute(*p, source)).collect();
        (there, back)
    }

    #[allow(clippy::cast_precision_loss)]
    fn share_within(side: &[f64], tolerance: f64) -> f64 {
        side.iter().filter(|d| **d <= tolerance).count() as f64 / side.len() as f64
    }

    /// The reported deviation against an independent measurement.
    fn check_readout(source: &Triangulation, done: &super::Simplified, tolerance: f64) {
        let deviation = done.deviation.unwrap();
        let (there, back) = both_ways(source, &done.mesh);
        let worst = there.iter().chain(&back).fold(0.0_f64, |m, d| m.max(*d));
        #[allow(clippy::cast_precision_loss)]
        let rms = (there.iter().chain(&back).map(|d| d * d).sum::<f64>()
            / (there.len() + back.len()) as f64)
            .sqrt();
        let within = share_within(&there, tolerance).min(share_within(&back, tolerance));
        assert!(
            (deviation.furthest.1 - worst).abs() < 1e-9,
            "furthest reported {} measured {worst}",
            deviation.furthest.1
        );
        assert!((done.error - worst).abs() < 1e-9);
        assert!(
            (deviation.rms - rms).abs() < 1e-9,
            "rms reported {} measured {rms}",
            deviation.rms
        );
        assert!(
            (deviation.within_share - within).abs() < 1e-12,
            "within reported {} measured {within}",
            deviation.within_share
        );
    }

    #[test]
    fn a_noisy_scan_holds_the_share_both_ways() {
        let source = noisy_sphere(0.02, 0.01);
        let target = Target::Verified {
            tolerance: 0.05,
            share: 0.999,
            max_triangles: None,
        };
        let done = simplify(&source, target, T).unwrap();
        let deviation = done.deviation.unwrap();
        eprintln!(
            "{} -> {} triangles, {deviation:?}",
            source.triangle_count(),
            done.mesh.triangle_count()
        );
        assert!(done.target_met);
        assert!(done.mesh.is_closed());

        let (there, back) = both_ways(&source, &done.mesh);
        assert!(share_within(&there, 0.05) >= 0.999);
        assert!(share_within(&back, 0.05) >= 0.999);
        check_readout(&source, &done, 0.05);
    }

    #[test]
    fn a_triangle_cap_forcing_more_reduction_reports_what_it_cost() {
        let source = noisy_sphere(0.02, 0.01);
        let target = Target::Verified {
            tolerance: 0.05,
            share: 0.999,
            max_triangles: Some(500),
        };
        let done = simplify(&source, target, T).unwrap();
        let deviation = done.deviation.unwrap();
        eprintln!(
            "capped: {} triangles, {deviation:?}",
            done.mesh.triangle_count()
        );
        assert!(done.mesh.triangle_count() <= 500);
        // The cap is far below what the tolerance allows, so the deviation
        // is past it and the target is not claimed.
        assert!(deviation.distance_at_share > 0.05);
        assert!(!done.target_met);
        check_readout(&source, &done, 0.05);

        let (there, back) = both_ways(&source, &done.mesh);
        let at_share = |side: &[f64]| {
            let mut sorted = side.to_vec();
            sorted.sort_by(f64::total_cmp);
            #[allow(
                clippy::cast_precision_loss,
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss
            )]
            let needed = (0.999 * sorted.len() as f64).ceil() as usize;
            sorted[needed - 1]
        };
        let measured = at_share(&there).max(at_share(&back));
        assert!(
            (deviation.distance_at_share - measured).abs() < 1e-9,
            "reported {} measured {measured}",
            deviation.distance_at_share
        );
    }

    /// A closed box of `cells` squares of side `size` along each axis, every
    /// face a grid of them, each square two triangles wound outward.
    fn plate(cells: [u32; 3], size: f64) -> Triangulation {
        let mut mesh = Triangulation::new();
        let mut numbered: ogeom_core::FastMap<[u32; 3], u32> = ogeom_core::FastMap::default();
        let mut vertex = |at: [u32; 3], mesh: &mut Triangulation| {
            *numbered.entry(at).or_insert_with(|| {
                mesh.positions.push(Point::new(
                    f64::from(at[0]) * size,
                    f64::from(at[1]) * size,
                    f64::from(at[2]) * size,
                ));
                u32::try_from(mesh.positions.len() - 1).unwrap()
            })
        };
        for axis in 0..3 {
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            for side in [0, cells[axis]] {
                for i in 0..cells[u] {
                    for j in 0..cells[v] {
                        let corner = |di: u32, dj: u32| {
                            let mut at = [0; 3];
                            at[axis] = side;
                            at[u] = i + di;
                            at[v] = j + dj;
                            at
                        };
                        let quad = [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)]
                            .map(|at| vertex(at, &mut mesh));
                        // (u, v, axis) is right-handed, so counter-clockwise in
                        // (u, v) faces along +axis: outward on the far side.
                        let [a, b, c, d] = if side == 0 {
                            [quad[0], quad[3], quad[2], quad[1]]
                        } else {
                            quad
                        };
                        mesh.triangles.push([a, b, c]);
                        mesh.triangles.push([a, c, d]);
                    }
                }
            }
        }
        mesh
    }

    #[test]
    fn a_flat_plate_comes_down_to_a_few_triangles_exactly() {
        let source = plate([40, 30, 2], 0.5);
        assert_eq!(source.triangle_count(), 4 * (40 * 30 + 30 * 2 + 2 * 40));
        let target = Target::Verified {
            tolerance: 1e-9,
            share: 1.0,
            max_triangles: None,
        };
        let done = simplify(&source, target, T).unwrap();
        let deviation = done.deviation.unwrap();
        eprintln!(
            "plate: {} -> {} triangles, {deviation:?}",
            source.triangle_count(),
            done.mesh.triangle_count()
        );
        assert!(done.target_met);
        assert!(done.mesh.is_closed());
        assert!(deviation.furthest.1 <= 1e-9);
        assert!((done.mesh.volume() - 20.0 * 15.0).abs() < 1e-9);
        // Midpoint collapses cannot always pair off an odd run of cells, so
        // not the box's twelve, but a hundredth of the source at most.
        assert!(
            done.mesh.triangle_count() * 100 <= source.triangle_count(),
            "{} triangles left",
            done.mesh.triangle_count()
        );
        check_readout(&source, &done, 1e-9);
    }

    #[test]
    fn the_same_mesh_decimates_the_same_way() {
        let source = noisy_sphere(0.05, 0.01);
        let target = Target::Verified {
            tolerance: 0.05,
            share: 0.999,
            max_triangles: None,
        };
        let one = simplify(&source, target, T).unwrap();
        let two = simplify(&source, target, T).unwrap();
        assert_eq!(one.mesh.triangles, two.mesh.triangles);
        assert_eq!(one.mesh.positions, two.mesh.positions);
        assert_eq!(one.deviation, two.deviation);
    }

    #[test]
    fn a_verified_target_that_means_nothing_is_refused() {
        let source = noisy_sphere(0.2, 0.0);
        let verified = |tolerance, share, max_triangles| Target::Verified {
            tolerance,
            share,
            max_triangles,
        };
        for target in [
            verified(0.0, 0.999, None),
            verified(-1.0, 0.999, None),
            verified(f64::NAN, 0.999, None),
            verified(0.05, 0.0, None),
            verified(0.05, 1.5, None),
            verified(0.05, f64::NAN, None),
            verified(0.05, 0.999, Some(0)),
        ] {
            assert!(simplify(&source, target, T).is_err(), "{target:?}");
        }
        let empty = simplify(&Triangulation::new(), verified(0.05, 0.999, None), T).unwrap();
        assert_eq!(empty.mesh.triangle_count(), 0);
        assert!(empty.target_met);
    }

    #[test]
    #[ignore = "heavy"]
    fn a_scan_of_over_two_hundred_thousand_triangles_holds_the_share() {
        let source = noisy_sphere(0.0007, 0.01);
        let target = Target::Verified {
            tolerance: 0.05,
            share: 0.999,
            max_triangles: None,
        };
        let start = ogeom_core::clock::Instant::now();
        let done = simplify(&source, target, T).unwrap();
        let took = start.elapsed();
        let deviation = done.deviation.unwrap();
        eprintln!(
            "scan: {} -> {} triangles in {took:?}, {deviation:?}",
            source.triangle_count(),
            done.mesh.triangle_count()
        );
        assert!(done.target_met);
        assert!(source.triangle_count() > 200_000);
        assert!(took < std::time::Duration::from_secs(80), "took {took:?}");

        // Brute force on a sample, both ways: nothing farther than the
        // reported furthest point, and that point where it says.
        let reach = deviation.furthest.1 + 1e-12;
        for p in source.positions.iter().step_by(211) {
            assert!(brute(*p, &done.mesh) <= reach);
        }
        for p in done.mesh.positions.iter().step_by(13) {
            assert!(brute(*p, &source) <= reach);
        }
        let at = deviation.furthest.0;
        // On its own side it is at no distance; on the other, at the reported one.
        let there = brute(at, &done.mesh).max(brute(at, &source));
        assert!(
            (there - deviation.furthest.1).abs() < 1e-9,
            "the furthest point is {there} from the other mesh"
        );
    }
}
