//! Recovering topology from a triangle soup.
//!
//! The inverse of tessellation, and not its mirror image: tessellation throws
//! information away, and no amount of care getting it back invents what was
//! lost. A mesh does not say which triangles were one face, what surface they
//! were cut from, or which chains of edges were one curve. Those have to be
//! *decided*, and the decisions are what this module is.
//!
//! # It recognises before it rebuilds
//!
//! Two stages, deliberately separable. [`planar_regions`] groups triangles that
//! lie in a common plane and reports what it found — which is useful on its own,
//! and is the honest answer for a mesh whose regions are not planar.
//! [`to_brep`] turns those regions into faces, and refuses a mesh it cannot
//! account for completely.
//!
//! # What it refuses, and why refusing is the point
//!
//! A curved region is refused rather than approximated. Fitting a cylinder to a
//! band of triangles is not hard; deciding that a band of triangles *is* a
//! cylinder rather than a smooth patch that resembles one is the whole problem,
//! and getting it wrong produces a solid that looks right, measures nearly
//! right, and has the wrong surface underneath every operation that follows.
//! [`crate::canonical`] is where that decision is made, when the caller
//! states a tolerance to make it at.
//!
//! So: planes are recovered exactly, because a plane through a set of coplanar
//! triangles is the plane they are in and there is nothing to decide. Anything
//! else is reported as what it is — a region this cannot name — and the caller
//! is told rather than handed a guess.

use std::collections::{HashMap, HashSet};

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Plane, Point, Point2, Vector};
use ogeom_topo::{Location, Model, Shape, Triangulation};

use ogeom_algo::{Built, History, make_face_on, make_shell, make_solid, make_wire};

/// Roles reconstruction assigns.
pub mod roles {
    use ogeom_core::Role;

    /// A face recovered from a group of coplanar triangles.
    pub const RECOVERED_FACE: Role = Role::op_defined(50);
}

/// A group of triangles found to lie in one plane.
#[derive(Debug, Clone)]
pub struct Region {
    /// The plane they lie in.
    pub plane: Plane,
    /// Which triangles of the mesh belong to it.
    pub triangles: Vec<usize>,
    /// Its boundary, as loops of mesh vertex indices.
    ///
    /// The first loop is the outer one; any others are holes. Each names its
    /// vertices in order and does not repeat the first at the end.
    pub loops: Vec<Vec<u32>>,
}

/// The smallest dihedral angle that counts as a real edge rather than a
/// tessellation seam.
///
/// Thirty degrees. See [`planar_regions`] for why this number is the whole
/// difficulty and not an implementation detail.
pub const CREASE: f64 = core::f64::consts::FRAC_PI_6;

/// Group a mesh's triangles into planar regions.
///
/// A triangle joins a region when its plane agrees with the region's — same
/// normal to within `tol.angular()`, and every one of its vertices within
/// `tol.confusion()` of the plane. Both tests are needed: normals alone would
/// merge two parallel faces on opposite sides of a slab into one region, and
/// distance alone would merge a face with a coplanar face pointing the other
/// way.
///
/// # Coplanar is not the same as flat
///
/// Every triangle is planar, and that is the trap. A cylinder tessellates into
/// quads between rulings, and a ruled quad *is* planar — so its two triangles
/// are genuinely coplanar and grow into a perfectly good two-triangle region.
/// Accepting those would rebuild a cylinder as a faceted drum, silently, and
/// every later operation would work on the facets.
///
/// So a region is only evidence of a plane if its *boundary is a boundary*:
/// every edge of it either has no triangle on the other side, or meets one
/// across a dihedral angle of at least `crease`. Below that the surface simply
/// continues, and a seam in the tessellation is not an edge of the model.
///
/// **That threshold is the difficulty, not a detail.** No angle is right for
/// every mesh: a coarse tessellation of a large cylinder turns by more than a
/// fine one of a small one, so a threshold that rejects the first accepts the
/// second. Deciding it from the mesh's own deflection is what canonical
/// recognition does, and `docs/PLAN.md` carries it. Until then the number is
/// a parameter with a stated default, so a caller can see it and choose.
///
/// Returns the regions found and the triangles that joined none. A caller that
/// needs the whole mesh accounted for checks that the leftovers are empty; one
/// that only wants the flat parts does not.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the mesh names
/// a vertex it does not have, or `crease` is not a positive angle.
pub fn planar_regions(
    mesh: &Triangulation,
    crease: f64,
    tol: Tolerances,
) -> OgeomResult<(Vec<Region>, Vec<usize>)> {
    if !crease.is_finite() || crease <= 0.0 {
        ogeom_bail!(Construction, "a crease angle of {crease} names no edge");
    }
    for triangle in &mesh.triangles {
        for index in triangle {
            if *index as usize >= mesh.positions.len() {
                ogeom_bail!(
                    Construction,
                    "a triangle names vertex {index}, and the mesh has {}",
                    mesh.positions.len()
                );
            }
        }
    }

    let planes: Vec<Option<Plane>> = mesh
        .triangles
        .iter()
        .map(|t| plane_of(mesh, *t, tol))
        .collect();
    let neighbours = adjacency(mesh);

    let mut taken = vec![false; mesh.triangles.len()];
    let mut regions = Vec::new();
    let mut leftover = Vec::new();

    for seed in 0..mesh.triangles.len() {
        if taken[seed] {
            continue;
        }
        let Some(plane) = planes[seed] else {
            // A sliver with no normal joins nothing; it is not evidence of a
            // surface, only of a tessellation that produced a degenerate.
            taken[seed] = true;
            leftover.push(seed);
            continue;
        };

        // Grow across shared edges only. Two coplanar patches that do not touch
        // are two faces, not one with a disconnected boundary.
        let mut members = vec![seed];
        taken[seed] = true;
        let mut frontier = vec![seed];
        while let Some(current) = frontier.pop() {
            for next in neighbours.get(&current).into_iter().flatten() {
                if taken[*next] {
                    continue;
                }
                if !belongs(mesh, mesh.triangles[*next], &plane, tol) {
                    continue;
                }
                taken[*next] = true;
                members.push(*next);
                frontier.push(*next);
            }
        }

        // Every edge out of the region has to be a real edge. One that is not
        // means the surface continues, and this region is a piece of something
        // curved rather than a face.
        if !bounded_by_creases(&members, &planes, &neighbours, crease) {
            leftover.extend(members.iter().copied());
            continue;
        }
        let loops = boundary_loops(mesh, &members);
        regions.push(Region {
            plane,
            triangles: members,
            loops,
        });
    }

    // A region of one degenerate triangle is not a region.
    regions.retain(|r| {
        if r.loops.is_empty() {
            leftover.extend(r.triangles.iter().copied());
            false
        } else {
            true
        }
    });
    leftover.sort_unstable();
    Ok((regions, leftover))
}

/// A curved region: triangles grown across smooth (sub-crease) edges, with
/// its boundary loops — the raw material canonical recognition decides on.
#[derive(Debug, Clone)]
pub struct CurvedRegion {
    /// Which triangles of the mesh belong to it.
    pub triangles: Vec<usize>,
    /// Its boundary, as loops of mesh vertex indices, outer first.
    pub loops: Vec<Vec<u32>>,
    /// One sample per triangle: centroid and unit normal.
    pub samples: Vec<(Point, Vector)>,
}

/// One chord of a recovered boundary: the edge, its chart line, and the sag
/// its own tolerance records.
///
/// Shared by index pair, so two regions meeting along a chord meet along one
/// edge — which is what makes the recovered shell close.
#[allow(clippy::too_many_arguments)]
fn chart_edge(
    model: &mut Model,
    mesh: &ogeom_topo::Triangulation,
    (from, to): (u32, u32),
    (pa2, pb2): (Point2, Point2),
    surface: &ogeom_geom::SurfaceGeometry,
    surface_id: ogeom_topo::SurfaceId,
    vertices: &mut HashMap<u32, Shape>,
    shared: &mut HashMap<(u32, u32), Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::Surface as _;
    let (a, b) = (mesh.positions[from as usize], mesh.positions[to as usize]);
    for (index, at) in [(from, a), (to, b)] {
        vertices
            .entry(index)
            .or_insert_with(|| model.add_vertex(ogeom_topo::VertexData::new(at)));
    }
    let key = (from.min(to), from.max(to));
    let edge = match shared.get(&key) {
        Some(existing) => existing.reversed(),
        None => {
            let fresh = ogeom_algo::make_edge_between(
                model,
                ogeom_geom::LineCurve::segment(a, b, tol)?.into(),
                (0.0, a.distance(b)),
                &vertices[&from].clone(),
                &vertices[&to].clone(),
                tol,
            )?
            .shape;
            shared.insert(key, fresh.clone());
            fresh
        }
    };
    // The pcurve is stated in the *edge's* own direction, which a reversed
    // occurrence traverses backwards.
    let (pa2, pb2) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
        (pb2, pa2)
    } else {
        (pa2, pb2)
    };
    let length = a.distance(b);
    let span = pb2 - pa2;
    if span.magnitude() <= f64::MIN_POSITIVE || length <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a curved boundary chord collapsed in the chart"
        );
    }
    let towards = ogeom_math::Direction2::new(span / length, tol).map_err(|_| {
        ogeom_core::ogeom_err!(
            Construction,
            "a curved boundary chord has no chart direction"
        )
    })?;
    ogeom_algo::attach_pcurve(
        model,
        &edge,
        ogeom_geom::Line2d::over(ogeom_math::Axis2::new(pa2, towards), 0.0, length)?.into(),
        surface_id,
        Location::identity(),
        (0.0, length),
    )?;
    // The chord stands off the recognized surface by its sag; the edge's own
    // tolerance records it.
    let mid2 = pa2 + span * 0.5;
    if let Ok(on_surface) = surface.point_at(mid2.x, mid2.y, tol) {
        let chord_mid = a + (b - a) * 0.5;
        let sag = on_surface.distance(chord_mid);
        if let Some(node) = model.node_mut(&edge)
            && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
        {
            data.tolerance = data.tolerance.widen_to(sag + tol.confusion());
        }
    }
    Ok(edge)
}

/// The one wire two full-period rings bound, closed by a seam.
///
/// Two rings that each wrap the whole chart bound a band, and in the chart a
/// band is not closed by its rims alone: the walk goes along the lower rim
/// for a turn, *up a seam*, back along the upper rim for a turn, and down the
/// seam again. The seam is one edge occurring twice, carrying the two chart
/// columns a turn apart, which is exactly what a cylinder built from first
/// principles has and what a face recovered from a mesh was missing.
#[allow(clippy::too_many_arguments)]
fn seamed_band(
    model: &mut Model,
    mesh: &ogeom_topo::Triangulation,
    first: (&[u32], &[Point2]),
    second: (&[u32], &[Point2]),
    advances: (f64, f64),
    period: f64,
    surface: &ogeom_geom::SurfaceGeometry,
    surface_id: ogeom_topo::SurfaceId,
    vertices: &mut HashMap<u32, Shape>,
    shared: &mut HashMap<(u32, u32), Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    // The lower rim is the one at the smaller chart height; the walk runs it
    // with `u` increasing and the upper one with `u` decreasing, which is the
    // sense that encloses the band rather than its complement.
    let height = |params: &[Point2]| -> f64 {
        #[expect(clippy::cast_precision_loss, reason = "a vertex count")]
        let count = params.len().max(1) as f64;
        params.iter().map(|p| p.y).sum::<f64>() / count
    };
    let (low, low_advance, high, high_advance) = if height(first.1) <= height(second.1) {
        (first, advances.0, second, advances.1)
    } else {
        (second, advances.1, first, advances.0)
    };
    let ordered = |ring: &[u32], params: &[Point2], forward: bool| -> (Vec<u32>, Vec<Point2>) {
        if forward {
            (ring.to_vec(), params.to_vec())
        } else {
            let mut r: Vec<u32> = ring.to_vec();
            let mut p: Vec<Point2> = params.to_vec();
            r.reverse();
            p.reverse();
            (r, p)
        }
    };
    let (low_ring, low_params) = ordered(low.0, low.1, low_advance > 0.0);
    let (high_ring, high_params) = ordered(high.0, high.1, high_advance < 0.0);
    if low_ring.len() < 3 || high_ring.len() < 3 {
        ogeom_bail!(
            Construction,
            "a rim of fewer than three chords bounds nothing"
        );
    }

    // Where the seam stands: at the lower rim's own start, and at whichever
    // vertex of the upper rim is nearest it round the chart — the most nearly
    // vertical seam there is, which is the one least likely to cross a rim.
    let start = low_params[0].x;
    let mut at = 0;
    let mut best = f64::INFINITY;
    for (k, p) in high_params.iter().enumerate() {
        let mut gap = (p.x - start).rem_euclid(period);
        if gap > period / 2.0 {
            gap -= period;
        }
        if gap.abs() < best {
            best = gap.abs();
            at = k;
        }
    }
    // The upper rim, rotated to start at the seam and shifted so its walk
    // begins a whole turn along and ends where the seam comes back down.
    let mut high_walk: Vec<u32> = high_ring[at..].to_vec();
    high_walk.extend_from_slice(&high_ring[..at]);
    let mut high_chart: Vec<Point2> = Vec::with_capacity(high_params.len() + 1);
    {
        let base = high_params[at].x;
        let mut previous = base;
        for k in 0..high_params.len() {
            let mut u = high_params[(at + k) % high_params.len()].x;
            while u - previous > period / 2.0 {
                u -= period;
            }
            while previous - u > period / 2.0 {
                u += period;
            }
            previous = u;
            high_chart.push(Point2::new(
                u + period,
                high_params[(at + k) % high_params.len()].y,
            ));
        }
        // Closing back onto the seam, a turn lower than it started.
        high_chart.push(Point2::new(base, high_params[at].y));
    }

    // The lower rim's walk: a whole turn, closing a period along.
    let mut low_chart: Vec<Point2> = low_params.clone();
    low_chart.push(Point2::new(low_params[0].x + period, low_params[0].y));

    let mut edges = Vec::with_capacity(low_ring.len() + high_ring.len() + 2);
    for i in 0..low_ring.len() {
        edges.push(chart_edge(
            model,
            mesh,
            (low_ring[i], low_ring[(i + 1) % low_ring.len()]),
            (low_chart[i], low_chart[i + 1]),
            surface,
            surface_id,
            vertices,
            shared,
            tol,
        )?);
    }

    // The seam itself: one edge from the lower rim's start to the upper rim's,
    // carrying both chart columns.
    let (from, to) = (low_ring[0], high_walk[0]);
    let (a, b) = (mesh.positions[from as usize], mesh.positions[to as usize]);
    for (index, point) in [(from, a), (to, b)] {
        vertices
            .entry(index)
            .or_insert_with(|| model.add_vertex(ogeom_topo::VertexData::new(point)));
    }
    let length = a.distance(b);
    if length <= tol.confusion() {
        ogeom_bail!(Construction, "the band's two rims meet, so it has no seam");
    }
    let seam = ogeom_algo::make_edge_between(
        model,
        ogeom_geom::LineCurve::segment(a, b, tol)?.into(),
        (0.0, length),
        &vertices[&from].clone(),
        &vertices[&to].clone(),
        tol,
    )?
    .shape;
    let column = |from2: Point2, to2: Point2| -> OgeomResult<ogeom_geom::PlanarCurve> {
        let span = to2 - from2;
        let towards = ogeom_math::Direction2::new(span / length, tol)
            .map_err(|_| ogeom_core::ogeom_err!(Construction, "the seam has no chart direction"))?;
        Ok(ogeom_geom::Line2d::over(ogeom_math::Axis2::new(from2, towards), 0.0, length)?.into())
    };
    // Forward is the walk *up* the far column; reversed is the near one, and
    // both are stated in the edge's own direction.
    let far = column(low_chart[low_chart.len() - 1], high_chart[0])?;
    let near = column(low_chart[0], high_chart[high_chart.len() - 1])?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        far,
        near,
        surface_id,
        Location::identity(),
        (0.0, length),
    )?;
    edges.push(seam.clone());

    for i in 0..high_walk.len() {
        edges.push(chart_edge(
            model,
            mesh,
            (high_walk[i], high_walk[(i + 1) % high_walk.len()]),
            (high_chart[i], high_chart[i + 1]),
            surface,
            surface_id,
            vertices,
            shared,
            tol,
        )?);
    }
    edges.push(seam.reversed());
    Ok(make_wire(model, &edges, tol)?.shape)
}

/// Distance to a surface's *unbounded* carrier, where it has one — the
/// adoption test must not be defeated by a window clamped tight around a
/// half-grown region.
fn carrier_distance(surface: &ogeom_geom::SurfaceGeometry, p: Point) -> Option<f64> {
    use ogeom_geom::SurfaceGeometry as S;
    Some(match surface {
        S::Plane(pl) => pl.plane().signed_distance_to(p).abs(),
        S::Cylinder(c) => c.cylinder().distance_to(p),
        S::Cone(c) => c.cone().distance_to(p),
        S::Sphere(s) => (p.distance(s.sphere().centre()) - s.sphere().radius()).abs(),
        S::Torus(t) => t.torus().distance_to(p),
        _ => return None,
    })
}

/// Segment the triangles no planar region claimed, *recognition-driven*:
/// a region grows only while the recognizer still accounts for it.
///
/// Tangent-smooth junctions — a fillet meeting its wall — are invisible to
/// the crease test, so a smooth blob can span several true surfaces. The
/// segmentation therefore grows greedily but validates: below the sample
/// floor it takes the smoothest continuation, and from there on a triangle
/// joins only if the grown set still recognizes as one surface. A blob no
/// segmentation of which recognizes refuses the rebuild, by name.
fn segmented_regions(
    mesh: &Triangulation,
    leftover: &[usize],
    crease: f64,
    recognition_tolerance: f64,
    tol: Tolerances,
    planar: &mut [Region],
    recognize: &SurfaceRecognizer<'_>,
) -> OgeomResult<Vec<(CurvedRegion, ogeom_geom::SurfaceGeometry)>> {
    let in_leftover: HashSet<usize> = leftover.iter().copied().collect();
    let neighbours = adjacency(mesh);
    let sample_of = |t: usize| -> Option<(Point, Vector)> {
        let [a, b, c] = mesh.triangles[t];
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        let n = (pb - pa).cross(pc - pa);
        let m = n.magnitude();
        if m <= tol.confusion() {
            return None;
        }
        let centroid = Point::new(
            (pa.x + pb.x + pc.x) / 3.0,
            (pa.y + pb.y + pc.y) / 3.0,
            (pa.z + pb.z + pc.z) / 3.0,
        );
        Some((centroid, n / m))
    };
    let try_recognize = |members: &[usize]| -> OgeomResult<Option<ogeom_geom::SurfaceGeometry>> {
        let mut points = Vec::with_capacity(members.len());
        let mut normals = Vec::with_capacity(members.len());
        for t in members {
            let Some((p, n)) = sample_of(*t) else {
                continue;
            };
            points.push(p);
            normals.push(n);
        }
        if points.len() < 3 {
            return Ok(None);
        }
        recognize(&points, &normals)
    };

    let mut taken: HashSet<usize> = HashSet::new();
    let mut grown: Vec<(Vec<usize>, ogeom_geom::SurfaceGeometry)> = Vec::new();
    // Seeds re-run until nothing changes: a triangle peeled back by one
    // region's validation gets its own chance to seed the next.
    let mut pool: Vec<usize> = leftover.to_vec();
    pool.sort_unstable();
    let mut passes = 0;
    loop {
        passes += 1;
        let claimed_before = taken.len();
        let mut seeds: Vec<usize> = pool
            .iter()
            .copied()
            .filter(|t| !taken.contains(t))
            .collect();
        // Interiors first: a seed whose neighbourhood is flattest sits deep
        // inside one surface, so the blind opening of its growth cannot
        // straddle a junction. Junction triangles join late, when validation
        // is already deciding.
        let interiorness = |t: &usize| -> u64 {
            let Some((_, n)) = sample_of(*t) else {
                return u64::MAX;
            };
            let worst = neighbours
                .get(t)
                .into_iter()
                .flatten()
                .filter_map(|other| {
                    sample_of(*other).map(|(_, m)| n.dot(m).clamp(-1.0, 1.0).acos())
                })
                .fold(0.0f64, f64::max);
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "a bounded angle key"
            )]
            let key = (worst * 1e9) as u64;
            key
        };
        seeds.sort_by_key(|t| (interiorness(t), *t));
        for seed in seeds {
            if taken.contains(&seed) || sample_of(seed).is_none() {
                taken.insert(seed);
                continue;
            }
            let mut members = vec![seed];
            taken.insert(seed);
            let mut surface: Option<ogeom_geom::SurfaceGeometry> = None;
            // A flat seed bootstraps smoothest-first and stays on its plane; a
            // curved seed bootstraps roughest-first so its opening samples span
            // the curvature instead of running along a coplanar strip.
            let seed_is_flat = {
                let local = sample_of(seed).map(|(_, n)| {
                    neighbours
                        .get(&seed)
                        .into_iter()
                        .flatten()
                        .filter_map(|other| {
                            sample_of(*other).map(|(_, m)| n.dot(m).clamp(-1.0, 1.0).acos())
                        })
                        .fold(0.0f64, f64::max)
                });
                local.is_some_and(|worst| worst < 1e-6)
            };
            loop {
                ogeom_core::progress::checkpoint()?;
                // Frontier: unclaimed leftover triangles smooth against some
                // member, ordered smoothest-first so the bootstrap stays on one
                // surface.
                let mut frontier: Vec<(f64, usize)> = Vec::new();
                for m in &members {
                    let Some((_, nm)) = sample_of(*m) else {
                        continue;
                    };
                    for next in neighbours.get(m).into_iter().flatten() {
                        if taken.contains(next) || !in_leftover.contains(next) {
                            continue;
                        }
                        let Some((_, nn)) = sample_of(*next) else {
                            continue;
                        };
                        let turn = nm.dot(nn).clamp(-1.0, 1.0).acos();
                        if turn < crease && !frontier.iter().any(|(_, t)| t == next) {
                            frontier.push((turn, *next));
                        }
                    }
                }
                // Bootstrap roughest-first: the opening samples must *span* the
                // region's curvature, or a curved band bootstraps along its own
                // coplanar strips and locks as the plane those strips are.
                // Once the kind is locked, validation decides and the order
                // stops mattering.
                if surface.is_none() && !seed_is_flat {
                    frontier.sort_by(|x, y| {
                        y.0.partial_cmp(&x.0).unwrap_or(core::cmp::Ordering::Equal)
                    });
                } else {
                    frontier.sort_by(|x, y| {
                        x.0.partial_cmp(&y.0).unwrap_or(core::cmp::Ordering::Equal)
                    });
                }
                let mut grew = false;
                // The first validation happens on the members alone, before any
                // ninth candidate — so the lock is set by what the seed's own
                // neighbourhood is, not by what a mixture happens to fit.
                if members.len() >= 8 && surface.is_none() {
                    surface = try_recognize(&members)?;
                    if surface.is_none() {
                        break;
                    }
                }
                for (_, candidate) in frontier {
                    if members.len() < 8 {
                        members.push(candidate);
                        taken.insert(candidate);
                        grew = true;
                        break;
                    }
                    let mut trial = members.clone();
                    trial.push(candidate);
                    if let Some(fit) = try_recognize(&trial)? {
                        // Kind lock: a region that has recognized as one kind
                        // may not morph into another as it grows — a nearly
                        // planar sample set also fits a giant cylinder, and
                        // without the lock a plane eats its neighbouring blend
                        // through exactly that loophole.
                        if let Some(held) = &surface
                            && core::mem::discriminant(held) != core::mem::discriminant(&fit)
                        {
                            continue;
                        }
                        members = trial;
                        taken.insert(candidate);
                        surface = Some(fit);
                        grew = true;
                        break;
                    }
                }
                if !grew {
                    break;
                }
            }
            let mut surface = match surface {
                Some(s) => Some(s),
                None => try_recognize(&members)?,
            };
            // The blind bootstrap can straddle a tangent junction: peel the
            // most recent additions back until what remains recognizes, and
            // return the peeled triangles to the pool for their own region.
            while surface.is_none() && members.len() > 1 {
                let popped = members.pop().unwrap_or(seed);
                taken.remove(&popped);
                if members.len() >= 3 {
                    surface = try_recognize(&members)?;
                }
            }
            if let Some(surface) = surface {
                grown.push((members, surface));
            } else if members.len() >= 3 {
                ogeom_bail!(
                    Construction,
                    "{} of {} triangles form a curved region canonical recognition \
                 does not account for. Fitting a surface anyway would give a \
                 solid that looks right and has the wrong surface underneath \
                 every later operation",
                    members.len(),
                    mesh.triangles.len()
                );
            } else {
                // Too small to recognize alone; the adoption pass below places
                // it with a neighbour whose surface accounts for it.
                for t in members {
                    taken.remove(&t);
                }
            }
        }

        if taken.len() == claimed_before || passes >= 8 {
            break;
        }
    }

    // Adoption: a leftover triangle joins an adjacent recognized region if
    // it lies on that region's own surface.
    let mut orphans: Vec<usize> = leftover
        .iter()
        .copied()
        .filter(|t| !taken.contains(t))
        .collect();
    orphans.sort_unstable();
    let mut settled = true;
    let mut dirty_planar = false;
    while settled {
        settled = false;
        orphans.retain(|orphan| {
            let Some((centroid, _)) = sample_of(*orphan) else {
                return false;
            };
            // Every vertex-sharing region is a candidate; the orphan goes
            // to the one whose carrier it actually lies on — best fit, not
            // first fit, because a corner triangle can graze several
            // planes within any reasonable gate.
            let orphan_vertices: Vec<u32> = mesh.triangles[*orphan].to_vec();
            let touches = |t: &usize| {
                mesh.triangles[*t]
                    .iter()
                    .any(|v| orphan_vertices.contains(v))
            };
            let mut candidates: Vec<(f64, usize, bool)> = Vec::new();
            for (index, region) in planar.iter().enumerate() {
                if region.triangles.iter().any(touches) {
                    candidates.push((region.plane.signed_distance_to(centroid).abs(), index, true));
                }
            }
            for (index, (members, surface)) in grown.iter().enumerate() {
                if members.iter().any(touches) {
                    let d = carrier_distance(surface, centroid)
                        .or_else(|| {
                            ogeom_algo::project_on_surface(surface, centroid, 12, tol)
                                .ok()
                                .map(|projection| projection.point.distance(centroid))
                        })
                        .unwrap_or(f64::INFINITY);
                    candidates.push((d, index, false));
                }
            }
            candidates.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(core::cmp::Ordering::Equal));
            eprintln!("DBGX orphan {orphan} at {centroid:?}");
            for (d, index, is_planar) in candidates {
                if d > recognition_tolerance.max(tol.confusion() * 10.0) {
                    continue;
                }
                if is_planar {
                    planar[index].triangles.push(*orphan);
                    dirty_planar = true;
                } else {
                    grown[index].0.push(*orphan);
                }
                settled = true;
                return false;
            }
            true
        });
    }
    if !orphans.is_empty() {
        ogeom_bail!(
            Construction,
            "{} of {} triangles form a curved region canonical recognition \
             does not account for. Fitting a surface anyway would give a \
             solid that looks right and has the wrong surface underneath \
             every later operation",
            orphans.len(),
            mesh.triangles.len()
        );
    }
    if dirty_planar {
        for region in planar.iter_mut() {
            region.loops = boundary_loops(mesh, &region.triangles);
        }
    }

    let mut out = Vec::new();
    for (members, surface) in grown {
        let loops = boundary_loops(mesh, &members);
        if loops.is_empty() {
            continue;
        }
        let samples = members.iter().filter_map(|t| sample_of(*t)).collect();
        out.push((
            CurvedRegion {
                triangles: members,
                loops,
                samples,
            },
            surface,
        ));
    }
    Ok(out)
}

/// Rebuild a mesh as topology.
///
/// Every triangle must belong to a planar region: a mesh with a curved part is
/// refused, naming how many triangles it could not account for. See the module
/// documentation for why that is a refusal rather than a fit.
///
/// The result is a solid when the recovered shell closes and a shell when it
/// does not, which is the same distinction the mesh itself makes.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the mesh is
/// empty, has triangles outside any planar region, or a region's boundary
/// cannot be built into a wire.
pub fn to_brep(
    model: &mut Model,
    mesh: &Triangulation,
    crease: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    to_brep_with(model, mesh, crease, tol.confusion() * 1e3, tol, &|_, _| {
        Ok(None)
    })
}

/// The kind of recognizer [`to_brep_with`] consumes: samples with normals
/// in, a recognized surface out — or `None`, which is a refusal to guess.
pub type SurfaceRecognizer<'a> =
    dyn Fn(&[Point], &[Vector]) -> OgeomResult<Option<ogeom_geom::SurfaceGeometry>> + 'a;

/// As [`to_brep`], with a canonical recognizer deciding the curved regions.
///
/// Triangles no planar region claims are grown into smooth curved regions
/// and handed to `recognize`; a region the recognizer declines still
/// refuses the whole rebuild, because fitting a wrong surface underneath is
/// worse than failing. Curved boundary edges keep the mesh's own chords,
/// shared with their planar neighbours — which is what closes the shell —
/// and each records the measured sag between chord and surface on its own
/// tolerance.
///
/// # Errors
///
/// As [`to_brep`].
pub fn to_brep_with(
    model: &mut Model,
    mesh: &Triangulation,
    crease: f64,
    recognition_tolerance: f64,
    tol: Tolerances,
    recognize: &SurfaceRecognizer<'_>,
) -> OgeomResult<Built> {
    if mesh.triangles.is_empty() {
        ogeom_bail!(Construction, "there are no triangles to rebuild from");
    }
    let (mut regions, leftover) = planar_regions(mesh, crease, tol)?;
    let recognized = segmented_regions(
        mesh,
        &leftover,
        crease,
        recognition_tolerance,
        tol,
        &mut regions,
        recognize,
    )?;
    model.begin_operation();

    // One vertex per mesh vertex that any boundary uses, shared between every
    // face that reaches it — which is what makes the shell close rather than
    // merely abut.
    let mut vertices: HashMap<u32, Shape> = HashMap::new();
    // And one *edge* per mesh edge, shared between the two faces that meet
    // along it. Building an edge per face instead leaves every one used once,
    // and the shell does not close — it abuts.
    let mut shared: HashMap<(u32, u32), Shape> = HashMap::new();
    let mut faces = Vec::with_capacity(regions.len());
    let mut history = History::new();

    for region in &regions {
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::PlaneSurface::new(region.plane).into());
        let mut wires = Vec::with_capacity(region.loops.len());
        for ring in &region.loops {
            let mut edges = Vec::with_capacity(ring.len());
            for i in 0..ring.len() {
                let (from, to) = (ring[i], ring[(i + 1) % ring.len()]);
                let (a, b) = (mesh.positions[from as usize], mesh.positions[to as usize]);
                for (index, at) in [(from, a), (to, b)] {
                    vertices
                        .entry(index)
                        .or_insert_with(|| model.add_vertex(ogeom_topo::VertexData::new(at)));
                }
                let key = (from.min(to), from.max(to));
                let edge = match shared.get(&key) {
                    // The neighbouring face built it, walking it the other way.
                    // The occurrence turns round; the edge does not.
                    Some(existing) => existing.reversed(),
                    None => {
                        let fresh = ogeom_algo::make_edge_between(
                            model,
                            ogeom_geom::LineCurve::segment(a, b, tol)?.into(),
                            (0.0, a.distance(b)),
                            &vertices[&from].clone(),
                            &vertices[&to].clone(),
                            tol,
                        )?
                        .shape;
                        shared.insert(key, fresh.clone());
                        fresh
                    }
                };
                // A pcurve per face, whichever way this face walks the edge:
                // the pcurve describes the edge in *this* surface's parameter
                // space, and without one the face cannot be triangulated.
                let (pa, pb) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                    (b, a)
                } else {
                    (a, b)
                };
                attach(model, &edge, region.plane, surface, pa, pb, tol)?;
                edges.push(edge);
            }
            wires.push(make_wire(model, &edges, tol)?.shape);
        }
        let face = make_face_on(model, surface, &wires, tol)?.shape;
        model.set_derived(&face, &[], roles::RECOVERED_FACE)?;
        faces.push(face);
    }

    for (region, curved_surface) in &recognized {
        use ogeom_geom::Surface as _;
        // The recognizer sizes its window from the samples it recognized on,
        // and a sample is a triangle's *centroid* — half a chord inside the
        // region's own rim. The face is trimmed to that rim, so its chart
        // reaches past the window the surface was given, and evaluating it
        // there is a domain error rather than a wrong answer. Widen the
        // window to hold the boundary before anything is built on it: the
        // extent of a surface is a parameterization window and not a trim,
        // so making it larger changes nothing about the face.
        let boundary: Vec<Point> = region
            .loops
            .iter()
            .flatten()
            .map(|index| mesh.positions[*index as usize])
            .collect();
        let curved_surface = &ogeom_algo::widened_to_hold(curved_surface, &boundary, tol)?;
        let surface_id = model.geometry_mut().add_surface(curved_surface.clone());
        let periodic_u = curved_surface.is_periodic_u();
        let u_period = {
            let ((ua, ub), _) = curved_surface.domain();
            ub - ua
        };
        // Chart parameters per ring vertex, the angle unwrapped so a ring
        // that crosses the chart's period stays continuous.
        let mut charts: Vec<Vec<Point2>> = Vec::with_capacity(region.loops.len());
        for ring in &region.loops {
            let mut params: Vec<Point2> = Vec::with_capacity(ring.len());
            for (k, index) in ring.iter().enumerate() {
                let p = mesh.positions[*index as usize];
                let projected = ogeom_algo::project_on_surface(curved_surface, p, 16, tol)?;
                let (mut u, v) = projected.parameters;
                if periodic_u && k > 0 {
                    let prev = params[k - 1].x;
                    while u - prev > u_period / 2.0 {
                        u -= u_period;
                    }
                    while prev - u > u_period / 2.0 {
                        u += u_period;
                    }
                }
                params.push(Point2::new(u, v));
            }
            charts.push(params);
        }

        // How far round the chart each ring travels. A ring that comes back
        // where it started travels nothing; one that wraps the whole period
        // travels a turn — and that is a face with no boundary in `u` at
        // all, whose chart region two such rings bound only once a *seam*
        // closes them. A drum recovered from its own mesh is exactly this.
        let wrap = |params: &[Point2]| -> f64 {
            if !periodic_u || params.len() < 2 {
                return 0.0;
            }
            let mut advance = params[params.len() - 1].x - params[0].x;
            // The closing step, unwrapped like the rest.
            let mut last = params[0].x - params[params.len() - 1].x;
            while last > u_period / 2.0 {
                last -= u_period;
            }
            while last < -u_period / 2.0 {
                last += u_period;
            }
            advance += last;
            advance
        };
        let wrapping: Vec<usize> = (0..charts.len())
            .filter(|&i| (wrap(&charts[i]).abs() - u_period).abs() < u_period * 1e-3)
            .collect();

        let mut wires = Vec::with_capacity(region.loops.len());
        if wrapping.len() == 2 {
            wires.push(seamed_band(
                model,
                mesh,
                (&region.loops[wrapping[0]], &charts[wrapping[0]]),
                (&region.loops[wrapping[1]], &charts[wrapping[1]]),
                (wrap(&charts[wrapping[0]]), wrap(&charts[wrapping[1]])),
                u_period,
                curved_surface,
                surface_id,
                &mut vertices,
                &mut shared,
                tol,
            )?);
        }
        for (index, ring) in region.loops.iter().enumerate() {
            if wrapping.len() == 2 && wrapping.contains(&index) {
                continue;
            }
            let params = &charts[index];
            let mut edges = Vec::with_capacity(ring.len());
            for i in 0..ring.len() {
                let (from, to) = (ring[i], ring[(i + 1) % ring.len()]);
                edges.push(chart_edge(
                    model,
                    mesh,
                    (from, to),
                    (params[i], params[(i + 1) % ring.len()]),
                    curved_surface,
                    surface_id,
                    &mut vertices,
                    &mut shared,
                    tol,
                )?);
            }
            wires.push(make_wire(model, &edges, tol)?.shape);
        }
        let face = make_face_on(model, surface_id, &wires, tol)?.shape;
        model.set_derived(&face, &[], roles::RECOVERED_FACE)?;
        faces.push(face);
    }

    let shell = make_shell(model, &faces)?.shape;
    for face in &faces {
        history.generate(face, shell.clone());
    }
    // Closed shells bound solids; open ones are sheets and stay sheets.
    if ogeom_algo::is_shell_closed(model, &shell)? {
        let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
        history.generate(&shell, solid.clone());
        return Ok(Built::new(solid, history));
    }
    Ok(Built::new(shell, history))
}

/// Attach a straight pcurve for an edge lying in a plane.
fn attach(
    model: &mut Model,
    edge: &Shape,
    plane: Plane,
    surface: ogeom_topo::SurfaceId,
    from: Point,
    to: Point,
    tol: Tolerances,
) -> OgeomResult<()> {
    let flat = |p: Point| {
        let local = plane.frame().to_local(p);
        Point2::new(local.x, local.y)
    };
    let (a, b) = (flat(from), flat(to));
    ogeom_algo::attach_pcurve(
        model,
        edge,
        ogeom_geom::Line2d::segment(a, b, tol)?.into(),
        surface,
        Location::identity(),
        (0.0, a.distance(b)),
    )
}

/// The plane a triangle lies in, or `None` if it has no area.
fn plane_of(mesh: &Triangulation, triangle: [u32; 3], tol: Tolerances) -> Option<Plane> {
    let [a, b, c] = triangle.map(|v| mesh.positions[v as usize]);
    let normal = Direction::from_cross(b - a, c - a, tol).ok()?;
    Some(Plane::new(Frame::about(a, normal)))
}

/// Whether a triangle belongs in a region's plane.
fn belongs(mesh: &Triangulation, triangle: [u32; 3], plane: &Plane, tol: Tolerances) -> bool {
    let Some(own) = plane_of(mesh, triangle, tol) else {
        return false;
    };
    // Same *direction*, not merely the same line: a coplanar face pointing the
    // other way is the far side of a zero-thickness sheet, not this face.
    if own.normal().dot(plane.normal()) < 1.0 - tol.angular() {
        return false;
    }
    triangle
        .iter()
        .all(|v| plane.distance_to(mesh.positions[*v as usize]) <= tol.confusion())
}

/// Whether every edge leaving a region is a real edge of the model.
fn bounded_by_creases(
    members: &[usize],
    planes: &[Option<Plane>],
    neighbours: &HashMap<usize, Vec<usize>>,
    crease: f64,
) -> bool {
    let inside: HashSet<usize> = members.iter().copied().collect();
    for index in members {
        let Some(own) = planes[*index] else {
            return false;
        };
        for next in neighbours.get(index).into_iter().flatten() {
            if inside.contains(next) {
                continue;
            }
            let Some(other) = planes[*next] else {
                continue;
            };
            // The angle between the two facets. A tessellation seam turns by a
            // little; an edge of the model turns by a lot.
            let turn = own.normal().dot(other.normal()).clamp(-1.0, 1.0).acos();
            if turn < crease {
                return false;
            }
        }
    }
    true
}

/// Which triangles share an edge with which.
fn adjacency(mesh: &Triangulation) -> HashMap<usize, Vec<usize>> {
    let mut across: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (i, triangle) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
            across.entry((a.min(b), a.max(b))).or_default().push(i);
        }
    }
    let mut out: HashMap<usize, Vec<usize>> = HashMap::new();
    for sharing in across.values() {
        for a in sharing {
            for b in sharing {
                if a != b {
                    out.entry(*a).or_default().push(*b);
                }
            }
        }
    }
    out
}

/// The boundary of a set of triangles, as ordered loops of vertices.
///
/// An edge is on the boundary when exactly one triangle *of the region* uses
/// it. Chaining those gives the region's outline, and a region with a hole
/// gives more than one chain — the outer one first, decided by which encloses
/// the most area in the region's own plane.
fn boundary_loops(mesh: &Triangulation, members: &[usize]) -> Vec<Vec<u32>> {
    let mut uses: HashMap<(u32, u32), (usize, (u32, u32))> = HashMap::new();
    for index in members {
        let triangle = mesh.triangles[*index];
        for k in 0..3 {
            let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
            let entry = uses.entry((a.min(b), a.max(b))).or_insert((0, (a, b)));
            entry.0 += 1;
            // Keep the direction the *first* user walked it, which is the
            // region's own winding.
            if entry.0 == 1 {
                entry.1 = (a, b);
            }
        }
    }

    // Directed boundary edges, keyed by where each starts.
    let mut onward: HashMap<u32, Vec<u32>> = HashMap::new();
    for (count, (a, b)) in uses.values() {
        if *count == 1 {
            onward.entry(*a).or_default().push(*b);
        }
    }

    let mut loops = Vec::new();
    let mut walked: HashSet<(u32, u32)> = HashSet::new();
    let starts: Vec<u32> = {
        let mut keys: Vec<u32> = onward.keys().copied().collect();
        // Deterministic: a boundary that comes back in a different order each
        // run gives a different model each run.
        keys.sort_unstable();
        keys
    };
    for start in starts {
        while let Some(first) = onward
            .get(&start)
            .and_then(|next| next.iter().find(|to| !walked.contains(&(start, **to))))
            .copied()
        {
            let mut ring = vec![start];
            walked.insert((start, first));
            let mut at = first;
            while at != start {
                ring.push(at);
                let Some(next) = onward
                    .get(&at)
                    .and_then(|next| next.iter().find(|to| !walked.contains(&(at, **to))))
                    .copied()
                else {
                    break;
                };
                walked.insert((at, next));
                at = next;
            }
            if ring.len() >= 3 {
                loops.push(ring);
            }
        }
    }

    // The outer loop first: it is the one enclosing the most area.
    loops.sort_by(|a, b| {
        enclosed(mesh, b)
            .partial_cmp(&enclosed(mesh, a))
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    loops
}

/// Twice the area a loop encloses, as a magnitude.
fn enclosed(mesh: &Triangulation, ring: &[u32]) -> f64 {
    let mut total = Vector::ZERO;
    for i in 0..ring.len() {
        let a = mesh.positions[ring[i] as usize];
        let b = mesh.positions[ring[(i + 1) % ring.len()] as usize];
        total += a.to_vector().cross(b.to_vector());
    }
    total.magnitude()
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_algo::{make_box, make_cylinder, make_wedge, shape_bounds};
    use ogeom_mesh::{Deflection, triangulate};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 0.01,
            ..Deflection::default()
        }
    }

    fn mesh_of(model: &Model, shape: &Shape) -> Triangulation {
        triangulate(model, shape, fine(), T).unwrap()
    }

    #[test]
    fn a_boxs_mesh_is_recognised_as_six_planar_regions() {
        // Twelve triangles, six faces: the pairs that share a diagonal have to
        // merge, and the pairs that merely touch at an edge must not.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let mesh = mesh_of(&model, &built.shape);

        let (regions, leftover) = planar_regions(&mesh, CREASE, T).unwrap();
        assert!(
            leftover.is_empty(),
            "{} triangles unaccounted for",
            leftover.len()
        );
        assert_eq!(regions.len(), 6);
        for region in &regions {
            assert_eq!(region.triangles.len(), 2);
            assert_eq!(region.loops.len(), 1, "a box face has no holes");
            assert_eq!(region.loops[0].len(), 4, "and four corners");
        }
    }

    #[test]
    fn a_box_survives_the_round_trip_through_a_mesh() {
        // The end-to-end claim: tessellate, throw the topology away, and get a
        // solid back that measures the same.
        let mut model = Model::new();
        let size = (2.0, 3.0, 4.0);
        let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
        let mesh = mesh_of(&model, &built.shape);

        let recovered = to_brep(&mut model, &mesh, CREASE, T).unwrap();
        assert_eq!(
            model.kind_of(&recovered.shape).unwrap(),
            ogeom_topo::ShapeType::Solid,
            "a closed mesh should come back a solid"
        );
        let counts = |kind| {
            ogeom_topo::explore_unique(&model, &recovered.shape, kind)
                .unwrap()
                .len()
        };
        assert_eq!(counts(ogeom_topo::ShapeType::Face), 6);
        assert_eq!(counts(ogeom_topo::ShapeType::Edge), 12);
        assert_eq!(counts(ogeom_topo::ShapeType::Vertex), 8);

        assert!(
            ogeom_algo::check(&model, &recovered.shape, T)
                .unwrap()
                .is_valid(),
            "{}",
            ogeom_algo::check(&model, &recovered.shape, T).unwrap()
        );
        let again = triangulate(&model, &recovered.shape, fine(), T).unwrap();
        assert!(again.is_closed());
        assert_relative_eq!(again.volume(), size.0 * size.1 * size.2, epsilon = 1e-9);
    }

    #[test]
    fn a_wedges_slanted_faces_come_back_slanted() {
        // Not axis-aligned, so a region grower that leaned on the axes would
        // merge the slanted sides with something.
        let mut model = Model::new();
        let built = make_wedge(&mut model, Frame::WORLD, (4.0, 4.0, 6.0), (2.0, 2.0), T).unwrap();
        let mesh = mesh_of(&model, &built.shape);
        let before = mesh.volume();

        let (regions, leftover) = planar_regions(&mesh, CREASE, T).unwrap();
        assert!(leftover.is_empty());
        assert_eq!(regions.len(), 6);

        let recovered = to_brep(&mut model, &mesh, CREASE, T).unwrap();
        let again = triangulate(&model, &recovered.shape, fine(), T).unwrap();
        assert_relative_eq!(again.volume(), before, epsilon = 1e-9);
    }

    #[test]
    fn a_curved_mesh_is_refused_rather_than_fitted() {
        // The decision this module exists to avoid making badly. A band of
        // triangles that resembles a cylinder is not evidence that it *is*
        // one, and a solid built on a guessed surface looks right and measures
        // nearly right while being wrong underneath every later operation.
        let mut model = Model::new();
        let built = make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        let mesh = mesh_of(&model, &built.shape);

        let (regions, leftover) = planar_regions(&mesh, CREASE, T).unwrap();
        assert_eq!(regions.len(), 2, "the two flat caps are recognised");
        assert!(!leftover.is_empty(), "the curved side is not");

        let err = to_brep(&mut model, &mesh, CREASE, T).unwrap_err();
        assert!(
            err.to_string().contains("canonical recognition"),
            "unexpected message: {err}"
        );
    }

    #[test]
    fn two_coplanar_patches_that_do_not_touch_are_two_regions() {
        // Coplanarity is not adhesion. Merging them would give one face with a
        // boundary in two pieces, which is a face nothing can triangulate.
        let mut model = Model::new();
        let here = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let there = ogeom_algo::transformed(
            &mut model,
            &here.shape,
            ogeom_math::Transform::translation(Vector::new(5.0, 0.0, 0.0)),
        )
        .unwrap();

        let mut both = mesh_of(&model, &here.shape);
        both.append(&mesh_of(&model, &there.shape));
        let (regions, _) = planar_regions(&both, CREASE, T).unwrap();
        assert_eq!(regions.len(), 12, "six faces each, none merged");
    }

    #[test]
    fn an_empty_or_broken_mesh_is_refused() {
        let mut model = Model::new();
        assert!(to_brep(&mut model, &Triangulation::new(), CREASE, T).is_err());

        let mut broken = Triangulation::new();
        broken.positions = vec![Point::ORIGIN];
        broken.triangles = vec![[0, 1, 2]];
        assert!(planar_regions(&broken, CREASE, T).is_err());
        assert!(to_brep(&mut model, &broken, CREASE, T).is_err());
        let _ = shape_bounds(&model, &Shape::of(ogeom_core::Key::from_parts(0, 0)), T);
    }
}
