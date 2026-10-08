//! Built faces against the triangles they replace: each face's own mesh,
//! the free edges of the built faces, and the regions whose face strays
//! past their triangles.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::Curve;
use ogeom_math::{Point, Vector};
use ogeom_topo::{EdgeRepr, Model, Shape, Triangulation};

use super::{Carrier, Groups};

/// Each built face's mesh drawn to the shape's agreed `chords`, taken from
/// the shape-wide pass's meshes (`drawn`) where it holds the face, and drawn
/// on its own where it does not; `None` for a face that does not draw.
pub(super) fn own_meshes(
    model: &Model,
    faces: &[(usize, &Shape)],
    drawn: Vec<ogeom_mesh::FaceMesh>,
    chords: &ogeom_mesh::EdgeChords,
    deflection: ogeom_mesh::Deflection,
    tol: Tolerances,
) -> Vec<Option<Triangulation>> {
    let key = |face: &Shape| (face.node(), face.orientation(), face.location().clone());
    let mut by_face: FastMap<_, Option<Triangulation>> =
        FastMap::with_capacity_and_hasher(drawn.len(), Default::default());
    for (face, mesh) in drawn {
        by_face.entry(key(&face)).or_insert_with(|| mesh.ok());
    }
    let missing: Vec<usize> = (0..faces.len())
        .filter(|&i| !by_face.contains_key(&key(faces[i].1)))
        .collect();
    let mut own = ogeom_core::parallel::map_ordered(&missing, |_, &i| {
        ogeom_mesh::triangulate_face_with(model, faces[i].1, deflection, chords, tol).ok()
    })
    .into_iter();
    faces
        .iter()
        .map(|&(_, face)| match by_face.get_mut(&key(face)) {
            Some(mesh) => mesh.take(),
            None => own.next().flatten(),
        })
        .collect()
}

/// A free edge as [`unmatched_faces`] reads it: its trimmed curve, how
/// far from it a point counts as on it, and the box holding every such
/// point.
///
/// [`unmatched_faces`]: super::checks::unmatched_faces
type FreeEdge = (Curve, f64, (Point, Point));

/// The edges of the built faces that bound one face only and are no seam
/// of it: each one's curve, trimmed to its range, its tolerance with the
/// confusion distance's margin, and a box holding every point within that
/// of it (its samples' box, widened by the tolerance and the longest step
/// between samples).
pub(super) fn free_edges(
    model: &Model,
    faces: &[(usize, &Shape)],
    tol: Tolerances,
) -> OgeomResult<Vec<FreeEdge>> {
    use ogeom_geom::Curve3d as _;
    let mut count: FastMap<ogeom_topo::SameKey, usize> = FastMap::default();
    for &(_, face) in faces {
        for edge in ogeom_topo::explore(
            model,
            face,
            ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Edge),
        )? {
            *count.entry(ogeom_topo::SameKey(edge)).or_default() += 1;
        }
    }
    let mut out = Vec::new();
    for (key, n) in count {
        if n != 1 {
            continue;
        }
        let Some(data) = model.node(&key.0).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        if data.degenerate {
            continue;
        }
        let Some(curve) = model.geometry().curve(*curve) else {
            continue;
        };
        let trimmed = ogeom_geom::TrimmedCurve::new(curve.clone(), range.0, range.1, tol)
            .map_or_else(|_| curve.clone(), |t| Curve::Trimmed(Box::new(t)));
        let reach = data.tolerance.get() + tol.confusion();
        let mut samples = Vec::with_capacity(33);
        for k in 0..=32 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 32.0;
            samples.push(curve.point_at(t, tol)?);
        }
        let step = samples
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0_f64, f64::max);
        let (mut lo, mut hi) = (samples[0], samples[0]);
        for p in &samples {
            lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let margin = Vector::new(1.0, 1.0, 1.0) * (reach + step);
        out.push((trimmed, reach, (lo - margin, hi + margin)));
    }
    Ok(out)
}

/// The most mesh triangles a planar face between curved ones may hold for
/// its seams with them to be threaded straight.
pub(super) const SLIVER_FACETS: usize = 2;

/// How much more a planar face's drawn triangles may cover than its exact
/// area before its trim is taken to cross itself: a planar face's mesh
/// covers its area exactly, but for the rounding of its boundary's chords.
pub(super) const OVERRUN: f64 = 1.2;

/// The recognized regions whose built face strays from the triangles it
/// replaces: its bounds stand beyond theirs by more than the surface bulges
/// past its facets and its boundary can run on to meet its neighbours, or
/// some point of it lies away from all of them.
#[allow(clippy::too_many_arguments, reason = "the build's state, read")]
pub(super) fn astray_faces(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    flat: f64,
    kept: &ogeom_mesh::FaceMeshCache,
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    let mut reach: Vec<Option<(Point, Point)>> = vec![None; groups.carriers.len()];
    let mut bulge = vec![0.0_f64; groups.carriers.len()];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, tri) in triangles.iter().enumerate() {
        let g = groups.of[t];
        let Some(Carrier::Curved(curved)) = groups.carriers.get(g) else {
            continue;
        };
        members[g].push(t);
        // The surface stands off a facet's edges by the sagitta, and the
        // face bulges past the facets' corners by as much.
        let [a, b, c] = tri.map(|v| points[v as usize]);
        for (p, q) in [(a, b), (b, c), (c, a)] {
            bulge[g] = bulge[g].max(curved.shape.distance_to(p + (q - p) * 0.5));
        }
        for &v in tri {
            let p = points[v as usize];
            let (lo, hi) = reach[g].get_or_insert((p, p));
            *lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            *hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); points.len()];
    for (t, tri) in triangles.iter().enumerate() {
        for &v in tri {
            around[v as usize].push(t);
        }
    }
    // Each planar region's triangles, and its face's drawn triangles'
    // middles, drawn once for every curved face beside it to ask.
    let mut plane_members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Plane(_))) {
            plane_members[g].push(t);
        }
    }
    let plane_middles: Vec<std::sync::OnceLock<Option<Vec<Point>>>> = (0..groups.carriers.len())
        .map(|_| std::sync::OnceLock::new())
        .collect();
    // Each face is judged on its own, and the answers taken in face order.
    let judged: Vec<(usize, &Shape, (Point, Point))> = built
        .iter()
        .enumerate()
        .filter_map(|(g, face)| Some((g, face.as_ref()?, reach[g]?)))
        .collect();
    let strays = ogeom_core::parallel::map_ordered(
        &judged,
        |_, &(g, face, (lo, hi))| -> OgeomResult<bool> {
            // Its boundary follows its neighbours' surfaces, which may meet it a
            // little past the last row of the mesh; a face closed the wrong way
            // round reaches a fair fraction of the region's size past it.
            let margin = bulge[g] * 2.0 + flat * 20.0 + lo.distance(hi) * 0.1;
            let bounds = crate::tight_bounds(model, face, tol)?;
            let (Some(flo), Some(fhi)) = (bounds.low(), bounds.high()) else {
                return Ok(false);
            };
            if flo.x < lo.x - margin
                || flo.y < lo.y - margin
                || flo.z < lo.z - margin
                || fhi.x > hi.x + margin
                || fhi.y > hi.y + margin
                || fhi.z > hi.z + margin
            {
                return Ok(true);
            }
            // Nor may it cover another part of its surface than they do: every
            // point of the face lies near one of them, by the surface's bulge
            // over them and the chord it is sampled at.
            let chord = (bulge[g] * 2.0).max(flat * 10.0);
            // A face that cannot be drawn cannot be vouched for either.
            let drawn = ogeom_mesh::triangulate_face_kept(
                model,
                face,
                ogeom_mesh::Deflection::with_chord(chord)?,
                kept,
                tol,
            )
            .or_else(|_| {
                ogeom_mesh::triangulate_face_kept(
                    model,
                    face,
                    ogeom_mesh::Deflection::default(),
                    kept,
                    tol,
                )
            });
            let Ok(mesh) = drawn else {
                return Ok(true);
            };
            if mesh.triangles.is_empty() {
                return Ok(true);
            }
            // Each triangle allows what the surface rises over it: a large facet
            // on a gentle curve lets the face stand well off its middle, and
            // lends nothing to the face anywhere else.
            let Some(Carrier::Curved(curved)) = groups.carriers.get(g) else {
                return Ok(false);
            };
            // The region's triangles and those touching them: the face's
            // boundary runs on to meet its neighbours' surfaces, a little past
            // its own last row.
            let near: Vec<usize> = {
                let mut near: Vec<usize> = members[g]
                    .iter()
                    .flat_map(|&t| triangles[t])
                    .flat_map(|v| around[v as usize].iter().copied())
                    .collect();
                near.sort_unstable();
                near.dedup();
                near
            };
            let rise: Vec<f64> = near
                .iter()
                .map(|&t| {
                    let [a, b, c] = triangles[t].map(|v| points[v as usize]);
                    [(a, b), (b, c), (c, a)]
                        .into_iter()
                        .map(|(p, q)| curved.shape.distance_to(p + (q - p) * 0.5))
                        .fold(0.0_f64, f64::max)
                })
                .collect();
            // A triangle standing far off the surface (a flat neighbour's facet
            // the region took in along a tangent) lends no more than a typical
            // one of the region does.
            let typical = {
                let mut own: Vec<f64> = near
                    .iter()
                    .zip(&rise)
                    .filter(|(t, _)| groups.of[**t] == g)
                    .map(|(_, &r)| r)
                    .collect();
                own.sort_by(f64::total_cmp);
                own.get(own.len() / 2).copied().unwrap_or(0.0)
            };
            let rise: Vec<f64> = rise.iter().map(|&r| r.min(typical * 4.0)).collect();
            let count = mesh.triangles.len();
            let samples = count.min(32);
            // A face running on past its last row where it meets a neighbour
            // strays at a few points by its end; one covering another part of
            // its surface strays over much of it.
            let wandering = (0..samples)
                .filter(|&k| {
                    let [a, b, c] =
                        mesh.triangles[k * count / samples].map(|i| mesh.positions[i as usize]);
                    let middle =
                        Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0);
                    near.iter().zip(&rise).all(|(&t, &rise)| {
                        let [p, q, r] = triangles[t].map(|v| points[v as usize]);
                        distance_to_triangle(middle, p, q, r)
                            > rise * 2.0 + chord * 2.0 + flat * 20.0
                    })
                })
                .count();
            if wandering * 4 > samples {
                return Ok(true);
            }
            // A plane beside the face takes its boundary from it, and a face a
            // little off its triangles carries that plane off its own: sampled
            // where it meets this region, each plane must lie on the triangles
            // too.
            let mut planes: Vec<usize> = near
                .iter()
                .map(|&t| groups.of[t])
                .filter(|&o| matches!(groups.carriers.get(o), Some(Carrier::Plane(_))))
                .collect();
            planes.sort_unstable();
            planes.dedup();
            // A gross check: a plane carried off by a face placed a little
            // wrong stands off by a fair part of the region's size, while one
            // meeting a rough neighbour a little past the last row stands off
            // by the roughness.
            let allowance = typical * 8.0 + chord * 2.0 + flat * 20.0 + lo.distance(hi) * 0.01;
            let reaches = |p: Point| {
                p.x >= lo.x - margin
                    && p.y >= lo.y - margin
                    && p.z >= lo.z - margin
                    && p.x <= hi.x + margin
                    && p.y <= hi.y + margin
                    && p.z <= hi.z + margin
            };
            for plane in planes {
                let Some(beside) = built.get(plane).and_then(Option::as_ref) else {
                    continue;
                };
                let Some(drawn) = plane_middles[plane].get_or_init(|| {
                    let drawn = ogeom_mesh::triangulate_face_kept(
                        model,
                        beside,
                        ogeom_mesh::Deflection::default(),
                        kept,
                        tol,
                    )
                    .ok()?;
                    Some(
                        drawn
                            .triangles
                            .iter()
                            .map(|t| {
                                let [a, b, c] = t.map(|i| drawn.positions[i as usize]);
                                Point::from_vector(
                                    (a.to_vector() + b.to_vector() + c.to_vector()) / 3.0,
                                )
                            })
                            .collect(),
                    )
                }) else {
                    continue;
                };
                let middles: Vec<Point> = drawn.iter().copied().filter(|&m| reaches(m)).collect();
                if middles.is_empty() {
                    continue;
                }
                let own = &plane_members[plane];
                let samples = middles.len().min(32);
                let off = (0..samples)
                    .filter(|&k| {
                        let m = middles[k * middles.len() / samples];
                        near.iter().chain(own).all(|&t| {
                            let [p, q, r] = triangles[t].map(|v| points[v as usize]);
                            distance_to_triangle(m, p, q, r) > allowance
                        })
                    })
                    .count();
                if off * 4 > samples {
                    return Ok(true);
                }
            }
            Ok(false)
        },
    );
    let mut astray = Vec::new();
    for (&(g, _, _), stray) in judged.iter().zip(strays) {
        if stray? {
            astray.push(g);
        }
    }
    Ok(astray)
}

/// The distance from `x` to the triangle `a b c`.
pub(crate) fn distance_to_triangle(x: Point, a: Point, b: Point, c: Point) -> f64 {
    let (ab, ac, ax) = (b - a, c - a, x - a);
    let n = ab.cross(ac);
    let area = n.magnitude();
    if area <= f64::MIN_POSITIVE {
        return x.distance(a).min(x.distance(b)).min(x.distance(c));
    }
    // Inside the prism over the triangle, the distance is the height above
    // its plane; outside, it is to the nearest side.
    let inside = [(a, b), (b, c), (c, a)]
        .iter()
        .all(|&(p, q)| (q - p).cross(x - p).dot(n) >= 0.0);
    if inside {
        return (ax.dot(n) / area).abs();
    }
    [(a, b), (b, c), (c, a)]
        .into_iter()
        .map(|(p, q)| {
            let d = q - p;
            let t = ((x - p).dot(d) / d.dot(d).max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
            x.distance(p + d * t)
        })
        .fold(f64::INFINITY, f64::min)
}
