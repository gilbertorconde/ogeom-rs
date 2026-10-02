//! A B-rep solid from a triangle mesh.
//!
//! A mesh already says which triangles meet along which edges: once its
//! repeated vertices are welded, two triangles that share two vertices
//! share an edge. So the topology is built from that connectivity
//! directly (one vertex per mesh vertex, one edge per mesh edge), and no
//! geometric sewing search is needed, which is what keeps a mesh of
//! hundreds of thousands of triangles cheap to convert.
//!
//! Adjacent triangles that lie in one plane merge into one planar face,
//! bounded by the region's outer loop and its holes, and a run of boundary
//! segments along one straight line between the same two faces becomes
//! one edge: a cube's twelve triangles become six faces and twelve edges,
//! and the faces take fillets and chamfers as a modelled box's do.
//!
//! Curved regions are recognized: triangles across which the surface
//! turns smoothly are grown into regions for as long as their vertices lie
//! on one cylinder, cone, sphere or torus, verified at the stated
//! tolerance, and the region is rebuilt on that surface. Its boundary with
//! each neighbour is placed on the surface exactly (a parallel circle or
//! a ruling of it), so a meshed bore comes back as a cylinder between two
//! circles, and a band all the way round its axis gets a seam as a
//! modelled one has. A sphere or torus with no boundary is one face, a
//! sphere bounded by one circle a cap, and a band round a torus's tube
//! runs between two of its meridians. A region whose boundary is no such
//! curve, or that nothing canonical fits, stays faceted.
//!
//! Windings are made consistent across each connected piece and turned
//! outward. A mesh that does not close (a hole, an edge three triangles
//! share, a piece that cannot be oriented) does not fail: it comes back
//! as open shells, with the report saying why.

use std::collections::HashMap;

use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, LineCurve, PlanarCurve, PlaneSurface};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Point2, Sphere, Torus, Vector};

use crate::recognize::{Canonical, recognize_curved, worst_deviation};
use ogeom_topo::{EdgeData, EdgeRepr, FaceData, Location, Model, Shape, Triangulation, VertexData};

/// How [`solid_from_mesh`] builds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshSolidOptions {
    /// Merge adjacent coplanar triangles into one planar face, and collinear
    /// boundary segments into one edge. Without it every triangle is a face
    /// and every mesh edge an edge.
    pub merge_coplanar: bool,
    /// The angle, in radians, below which two triangles' normals count as
    /// parallel for merging.
    pub coplanar_angle: f64,
    /// How far a triangle's corner may sit off a face's plane, or a
    /// merged edge's intermediate vertex off its line, and still be merged.
    /// `None` takes the largest of twice the mesh's `quantum`, a millionth
    /// of its bounding box's diagonal, one and a half times how far its
    /// vertices stand off the flat faces they lie on (measured across
    /// nearly coplanar edges), and the weld distance. With `recognize`, a
    /// distance given below two and a half times that scatter (or twice
    /// the `quantum`) is raised to it, and the report says so: no surface
    /// fitted to the vertices holds them closer. The
    /// vertices and edges of a merged face widen their tolerances to cover
    /// what they stand off it.
    pub coplanar_distance: Option<f64>,
    /// How far the mesh's own encoding may have moved a vertex: the
    /// rounding of its coordinates (see
    /// [`single_precision_quantum`] for a mesh held in `f32`, and the STL
    /// reader's for a file). A vertex that far off its surface is on it;
    /// `None` knows nothing of the encoding.
    pub quantum: Option<f64>,
    /// Vertices closer than this are welded into one before building; an
    /// STL repeats every vertex once per triangle that uses it. `None` is
    /// the confusion tolerance.
    pub weld: Option<f64>,
    /// Recognize curved regions as cylinders, cones, spheres and tori, at
    /// the coplanar distance. Needs `merge_coplanar`.
    pub recognize: bool,
    /// The angle, in radians, from which a turn between two triangles is a
    /// crease (an edge of the model) rather than the surface curving on:
    /// thirty degrees by default. A mesh drawn coarser than this round a
    /// curve reads as facets.
    pub crease: f64,
    /// Keep every mesh vertex a vertex of the solid: collinear boundary
    /// segments stay edges of their own. A solid built without recognition
    /// to be refined later ([`refine_solid`]) keeps the samples that
    /// recognition reads, where merged edges along a curved stretch's
    /// facets would have dropped all but their ends.
    pub keep_vertices: bool,
    /// After the canonical surfaces, rebuild a smooth region none of them
    /// fits as an extrusion of a fitted profile (its normals all square to
    /// one direction) or a surface of revolution (its normals all meeting
    /// one axis), verified at every vertex. Needs `recognize`.
    pub sweeps: bool,
    /// After the sweeps, rebuild a smooth region nothing else fits as one
    /// fitted B-spline patch, where the region is one disk bounded by one
    /// loop and the patch verifies: every vertex within the coplanar
    /// distance of it, every triangle's interior within that and its own
    /// sag, both ways, and the patch's normal regular and agreeing with
    /// every triangle's. Its chart is a nearly fitting canonical surface's
    /// where one fits within ten times the distance, and otherwise the
    /// mean-value map of the region onto a square. A region that is not a
    /// disk, is too narrow to hold a patch across it, or whose patch does
    /// not verify stays faceted and is counted
    /// ([`MeshSolidReport::patches_not_disk`],
    /// [`MeshSolidReport::patches_narrow`],
    /// [`MeshSolidReport::patches_unverified`]). A patch whose seams cannot
    /// be built falls back to facets as any recognized face does. Needs
    /// `recognize`.
    pub patches: bool,
}

impl Default for MeshSolidOptions {
    fn default() -> Self {
        Self {
            merge_coplanar: true,
            coplanar_angle: 1e-3,
            coplanar_distance: None,
            quantum: None,
            weld: None,
            recognize: true,
            crease: core::f64::consts::FRAC_PI_6,
            keep_vertices: false,
            sweeps: true,
            patches: true,
        }
    }
}

/// How far single-precision storage may have moved a mesh's vertices: half
/// the spacing of `f32` values at its largest coordinate on each axis, over
/// the three axes at once. The quantum of a mesh that went through `f32`
/// (a render mesh, a binary STL), for [`MeshSolidOptions::quantum`].
#[must_use]
pub fn single_precision_quantum(mesh: &Triangulation) -> f64 {
    let largest = mesh
        .positions
        .iter()
        .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
        .fold(0.0_f64, f64::max);
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the magnitude is an f32 coordinate's"
    )]
    let x = (largest as f32).max(f32::MIN_POSITIVE);
    let step = f64::from(f32::from_bits(x.to_bits() + 1) - x);
    step * 3.0_f64.sqrt() / 2.0
}

/// What [`solid_from_mesh`] did, and where the mesh does not close.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MeshSolidReport {
    /// Triangles built from, after the dropped ones.
    pub triangles: usize,
    /// Faces built: the triangles, or their coplanar groups and recognized
    /// regions.
    pub faces: usize,
    /// Faces built on a recognized curved surface.
    pub curved_faces: usize,
    /// Regions recognized as curved whose boundary could not be placed on
    /// the surface exactly, and which were faceted instead.
    pub curved_faceted: usize,
    /// Mesh vertices that welded onto another.
    pub vertices_welded: usize,
    /// Triangles dropped for having no area.
    pub degenerate_dropped: usize,
    /// Triangles dropped for repeating another's three vertices.
    pub duplicates_dropped: usize,
    /// Triangles whose winding was reversed to agree with their neighbours
    /// and face outward.
    pub windings_flipped: usize,
    /// Mesh edges only one triangle uses: the rims of holes.
    pub edges_used_once: usize,
    /// Mesh edges three or more triangles use: non-manifold.
    pub edges_used_more: usize,
    /// Shared edges whose windings cannot be made to agree, in a piece that
    /// is not orientable.
    pub orientation_conflicts: usize,
    /// Connected pieces, each one shell.
    pub shells: usize,
    /// Whether recognition was withdrawn from a body because the solid it
    /// built did not hold that body's mesh volume, and the body was built
    /// faceted instead.
    pub recognition_withdrawn: bool,
    /// Whether the coplanar distance asked for was below what the mesh's
    /// own scatter allows recognition to hold to, and was raised to that
    /// ([`MeshSolid::coplanar_distance`] is the one used).
    pub coplanar_distance_raised: bool,
    /// Faces built on a fitted B-spline patch, among the curved faces.
    pub patch_faces: usize,
    /// Of those, the ones whose chart is the mean-value map onto a square
    /// rather than a nearly fitting canonical surface's.
    pub patch_charts_mapped: usize,
    /// Smooth regions nothing else fitted that were left faceted because
    /// they are not one disk bounded by one loop (a region with a hole, or
    /// one that touches itself).
    pub patches_not_disk: usize,
    /// Smooth regions nothing else fitted that were left faceted because
    /// too few of their vertices lie inside their boundary for a patch to be
    /// held across them (a row or two of facets along a seam).
    pub patches_narrow: usize,
    /// Smooth regions nothing else fitted whose fitted patch failed its
    /// verification, left faceted.
    pub patches_unverified: usize,
    /// Free edges of curved faces (a run of the mesh's boundary with one
    /// face and nothing across it) that no parallel or ruling of the face
    /// holds, built as a curve fitted onto the face through the feet of
    /// the run's points, its tolerance measured against every vertex of the
    /// run. The others are parallels and rulings, placed exactly.
    pub free_edges_fitted: usize,
}

/// The result of [`solid_from_mesh`].
#[derive(Debug, Clone)]
pub struct MeshSolid {
    /// A solid when the mesh closes (a compound of solids when it is
    /// several disjoint closed pieces), and otherwise a shell, or a
    /// compound of shells.
    pub shape: Shape,
    /// Whether every piece closed and became a solid.
    pub closed: bool,
    /// The coplanar distance the faces were built to: the one asked for,
    /// or the default, as widened where the recognized surfaces pressed
    /// against it or raised to the mesh's scatter (see
    /// [`MeshSolidReport::coplanar_distance_raised`]).
    pub coplanar_distance: f64,
    /// What was built, and where the mesh does not close.
    pub report: MeshSolidReport,
}

/// Seal the cracks flat slivers leave when they are dropped.
///
/// A triangle with three distinct corners and no area has one corner on the
/// segment between the other two. In a closed mesh the triangle across that
/// long side fills the other side of it, and the two short sides are shared
/// with the neighbours beyond: dropped alone, the sliver leaves all three
/// used once. The triangle across the long side is split at the middle
/// corner instead, so each short side is shared again. Slivers against
/// slivers resolve as their neighbours are split, round by round.
fn split_across_slivers(points: &[Point], triangles: &mut Vec<[u32; 3]>, slivers: &[[u32; 3]]) {
    let mut pending: Vec<[u32; 3]> = slivers
        .iter()
        .map(|&[a, b, c]| {
            // The middle corner is the one opposite the longest side.
            let long = |x: u32, y: u32| points[x as usize].distance(points[y as usize]);
            let sides = [
                (long(b, c), a, [b, c]),
                (long(c, a), b, [c, a]),
                (long(a, b), c, [a, b]),
            ];
            let (_, middle, [p, q]) =
                sides.into_iter().fold(
                    sides[0],
                    |best, side| if side.0 > best.0 { side } else { best },
                );
            [p, middle, q]
        })
        .collect();
    // The triangles on each undirected edge, by index, lowest first: the
    // triangle across a sliver's long side is found from its edge, not by
    // scanning the mesh, and the lowest index is the one a scan would meet
    // first.
    let key = |x: u32, y: u32| (x.min(y), x.max(y));
    let mut on_edge: std::collections::HashMap<(u32, u32), Vec<usize>> =
        std::collections::HashMap::new();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            on_edge
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(i);
        }
    }
    let unlink =
        |on_edge: &mut std::collections::HashMap<(u32, u32), Vec<usize>>, t: [u32; 3], i: usize| {
            for k in 0..3 {
                if let Some(list) = on_edge.get_mut(&key(t[k], t[(k + 1) % 3]))
                    && let Ok(at) = list.binary_search(&i)
                {
                    list.remove(at);
                }
            }
        };
    let link =
        |on_edge: &mut std::collections::HashMap<(u32, u32), Vec<usize>>, t: [u32; 3], i: usize| {
            for k in 0..3 {
                let list = on_edge.entry(key(t[k], t[(k + 1) % 3])).or_default();
                if let Err(at) = list.binary_search(&i) {
                    list.insert(at, i);
                }
            }
        };
    loop {
        let mut progress = false;
        let mut left = Vec::new();
        for [p, middle, q] in pending {
            let across = on_edge.get(&key(p, q)).and_then(|list| {
                list.iter()
                    .copied()
                    .find(|&i| !triangles[i].contains(&middle))
            });
            let Some(index) = across else {
                left.push([p, middle, q]);
                continue;
            };
            let t = triangles[index];
            let Some(k) = (0..3).find(|&k| {
                let (x, y) = (t[k], t[(k + 1) % 3]);
                (x == p && y == q) || (x == q && y == p)
            }) else {
                continue;
            };
            let (x, y, z) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            unlink(&mut on_edge, t, index);
            triangles[index] = [x, middle, z];
            link(&mut on_edge, triangles[index], index);
            triangles.push([middle, y, z]);
            link(&mut on_edge, [middle, y, z], triangles.len() - 1);
            progress = true;
        }
        pending = left;
        if pending.is_empty() || !progress {
            break;
        }
    }
}

/// Rebuild a solid's faceted stretches on the surfaces they approximate.
///
/// The second of two steps: [`solid_from_mesh`] without recognition keeps
/// a mesh as it is, every flat stretch one face and every curved one its
/// facets; this finds the cylinders, cones, spheres and tori among the
/// facets of such a solid (or of any solid) and builds them anew as
/// `solid_from_mesh` with recognition would. The solid is triangulated,
/// its flat faces exactly and any curved ones to a chord of a
/// ten-thousandth of its size, and the triangles read as a mesh. The first
/// step reads best with [`MeshSolidOptions::keep_vertices`], so the
/// solid's edges keep the mesh's samples along curved stretches; pass its
/// `quantum` in `options` where the mesh came from a file.
///
/// # Errors
///
/// As [`solid_from_mesh`], and where the solid cannot be triangulated.
pub fn refine_solid(
    model: &mut Model,
    shape: &Shape,
    options: &MeshSolidOptions,
    tol: Tolerances,
) -> OgeomResult<MeshSolid> {
    let bounds = crate::tight_bounds(model, shape, tol)?;
    let size = match (bounds.low(), bounds.high()) {
        (Some(lo), Some(hi)) => lo.distance(hi),
        _ => ogeom_bail!(Construction, "the shape has no extent to refine"),
    };
    let chord = (size * 1e-4).max(tol.confusion() * 10.0);
    let mesh = ogeom_mesh::triangulate(
        model,
        shape,
        ogeom_mesh::Deflection::with_chord(chord)?,
        tol,
    )?;
    let options = MeshSolidOptions {
        recognize: true,
        keep_vertices: false,
        ..*options
    };
    solid_from_mesh(model, &mesh, &options, tol)
}

/// How near its distance a recognized fit must come for the distance, not
/// the surface, to be what bounds it.
const PRESSED: f64 = 0.7;

/// Build a B-rep from a triangle mesh.
///
/// See the module documentation for the construction. A closed piece
/// inside another becomes a void of the solid around it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// mesh has no triangle with area, an index names a vertex that is not
/// there, or an option is not finite and positive.
pub fn solid_from_mesh(
    model: &mut Model,
    mesh: &Triangulation,
    options: &MeshSolidOptions,
    tol: Tolerances,
) -> OgeomResult<MeshSolid> {
    let weld = options.weld.unwrap_or_else(|| tol.confusion());
    if !(weld.is_finite() && weld > 0.0) {
        ogeom_bail!(
            Construction,
            "the weld distance {weld} is not finite and positive"
        );
    }
    if !(options.coplanar_angle.is_finite() && options.coplanar_angle >= 0.0) {
        ogeom_bail!(Construction, "the coplanar angle is not finite");
    }
    if let Some(distance) = options.coplanar_distance
        && !(distance.is_finite() && distance > 0.0)
    {
        ogeom_bail!(
            Construction,
            "the coplanar distance {distance} is not finite and positive"
        );
    }
    if let Some(quantum) = options.quantum
        && !(quantum.is_finite() && quantum >= 0.0)
    {
        ogeom_bail!(
            Construction,
            "the quantum {quantum} is not finite and non-negative"
        );
    }
    // A crease turns by more than nothing and no more than a half turn: at
    // zero every turn is a crease, and nothing turns further than back on
    // itself.
    if !(options.crease.is_finite()
        && options.crease > 0.0
        && options.crease <= core::f64::consts::PI)
    {
        ogeom_bail!(
            Construction,
            "the crease angle {} is not a turn between nothing and a half turn",
            options.crease
        );
    }
    let count = mesh.positions.len();
    if mesh
        .triangles
        .iter()
        .flatten()
        .any(|&v| v as usize >= count)
    {
        ogeom_bail!(
            Construction,
            "a triangle names a vertex past the {count} the mesh has"
        );
    }
    let mut report = MeshSolidReport::default();

    let (points, remap) = weld_points(&mesh.positions, weld);
    report.vertices_welded = count - points.len();
    let mut triangles = Vec::with_capacity(mesh.triangles.len());
    let mut seen: HashMap<[u32; 3], ()> = HashMap::with_capacity(mesh.triangles.len());
    let mut flat_slivers: Vec<[u32; 3]> = Vec::new();
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|v| remap[v as usize]);
        if a == b || b == c || c == a || !has_area(&points, [a, b, c], weld) {
            report.degenerate_dropped += 1;
            if a != b && b != c && c != a {
                flat_slivers.push([a, b, c]);
            }
            continue;
        }
        let mut key = [a, b, c];
        key.sort_unstable();
        if seen.insert(key, ()).is_some() {
            report.duplicates_dropped += 1;
            continue;
        }
        triangles.push([a, b, c]);
    }
    if triangles.is_empty() {
        ogeom_bail!(Construction, "the mesh has no triangle with area");
    }
    split_across_slivers(&points, &mut triangles, &flat_slivers);
    report.triangles = triangles.len();

    // Orient each piece consistently, then outward where it closes.
    let adjacency = Adjacency::new(&triangles);
    report.edges_used_once = adjacency.used_once;
    report.edges_used_more = adjacency.used_more;
    let pieces = orient(&points, &mut triangles, &adjacency, &mut report);
    // A closed piece nested an odd number of pieces deep bounds a void of
    // the one around it, and faces into the void, out of the material.
    let all_closed = pieces.iter().all(|p| p.closed);
    let depth: Vec<usize> = if all_closed {
        (0..pieces.len())
            .map(|i| {
                (0..pieces.len())
                    .filter(|&j| j != i && inside(&points, &triangles, &pieces[j], &pieces[i]))
                    .count()
            })
            .collect()
    } else {
        vec![0; pieces.len()]
    };
    for (piece, d) in pieces.iter().zip(&depth) {
        if d % 2 == 1 {
            for &t in &piece.triangles {
                triangles[t as usize].swap(1, 2);
                report.windings_flipped += 1;
            }
        }
    }
    // Flipping renumbers every triangle's edges; the twins are found again.
    let mut adjacency = Adjacency::new(&triangles);
    report.shells = pieces.len();
    for _ in 0..3 {
        if unfold(&points, &mut triangles, &adjacency) == 0 {
            break;
        }
        adjacency = Adjacency::new(&triangles);
    }

    let diagonal = diagonal(&points);
    let noise = flat_noise(&points, &triangles, &adjacency);
    let quantum = options.quantum.unwrap_or(0.0);
    let mut flat = options
        .coplanar_distance
        .unwrap_or_else(|| (1e-6 * diagonal).max(2.0 * quantum).max(1.5 * noise))
        .max(weld);
    // A distance asked below what the mesh's own scatter allows cannot be
    // met: no surface fitted to the vertices holds them closer, regions
    // break into fragments, and the seams between the fragments open. A
    // seam's points lie within the distance of two surfaces each fitted
    // to vertices that scatter by the noise, and the noise is measured at
    // its ninetieth percentile, so recognition holds to at least
    // [`NOISE_FLOOR`] times it, and says so.
    if options.recognize && options.coplanar_distance.is_some() {
        let floor = (2.0 * quantum).max(NOISE_FLOOR * noise);
        if flat < floor {
            flat = floor;
            report.coplanar_distance_raised = true;
        }
    }
    let mut groups = segment(&points, &triangles, &adjacency, options, flat, tol)?;
    // The default distance is what single precision resolves, and some
    // exporters place their vertices a few times farther off their own
    // surfaces than that. The recognized surfaces say so: where their fits
    // press against the distance, it is the distance that stops them, and
    // their rims' last triangles stay facets. Unless the caller chose the
    // distance, it widens while that holds, twice at most.
    if options.coplanar_distance.is_none() && options.recognize {
        for _ in 0..2 {
            let pressed = groups
                .carriers
                .iter()
                .filter_map(|c| match c {
                    Carrier::Curved(curved) => Some(curved.fitted),
                    _ => None,
                })
                .fold(0.0_f64, f64::max);
            if pressed <= flat * PRESSED {
                break;
            }
            flat *= 2.0;
            groups = segment(&points, &triangles, &adjacency, options, flat, tol)?;
        }
    }
    // Facets left between curved faces go to them first. Whether that
    // leaves every face facing out is known only once the solid is built;
    // where one faces into the material, the conversion is made again
    // without them.
    let unabsorbed = (groups.clone(), report.clone());
    let mut absorbed = if options.recognize {
        absorb_facets(&points, &triangles, &adjacency, &mut groups, flat)
    } else {
        HashMap::new()
    };
    let mut absorbing = !absorbed.is_empty();
    let snaps = std::sync::Mutex::new(SnapCache::new());
    let shape = 'attempt: loop {
        // Plan until every curved face's boundary is exact, faceting the ones
        // whose boundary is not; then build, and facet any recognized face that
        // reaches past the triangles it replaces (a boundary placed on the wrong
        // turn of its surface closes a face of the wrong extent) and build again.
        let mut pinned: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut straight: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
        // What stood before any seam was threaded straight: a straightened
        // build with a face turned into its material is set aside for it.
        let mut unthreaded: Option<(Groups, MeshSolidReport)> = None;
        let mut threading_refused = false;
        // Whether faces were faceted after the first seam was threaded: the
        // build then differs from the unthreaded one beyond those seams.
        let mut refaceted = false;
        let shape = loop {
            let planner = Planner {
                points: &points,
                triangles: &triangles,
                adjacency: &adjacency,
                groups: &groups,
                merge: options.merge_coplanar && !options.keep_vertices,
                pinned: &pinned,
                straight: &straight,
                crease: options.crease,
                flat,
                tol,
                snaps: &snaps,
            };
            let failed = match planner.plan()? {
                Err(Replan::Pin(vertices)) => {
                    pinned.extend(vertices);
                    continue;
                }
                Err(Replan::Facet(failed)) => failed,
                Ok(plan) => {
                    model.begin_operation();
                    let built = Builder {
                        model,
                        points: &points,
                        triangles: &triangles,
                        groups: &groups,
                        plan: &plan,
                        tol,
                    }
                    .build()?;
                    let astray =
                        astray_faces(model, &points, &triangles, &groups, &built, flat, tol)?;
                    if astray.is_empty() {
                        let (shape, bodies) = assemble(
                            model, &points, &triangles, &pieces, &depth, all_closed, &groups,
                            &built,
                        )?;
                        let mut culprits = if options.recognize && all_closed {
                            body_culprits(
                                model,
                                &points,
                                &triangles,
                                &groups,
                                &built,
                                &bodies,
                                flat,
                                tol,
                                &mut report,
                            )?
                        } else {
                            Vec::new()
                        };
                        if culprits.is_empty() && options.recognize && !threading_refused {
                            let crossed = crossed_seams(
                                model, &shape, &triangles, &adjacency, &groups, &built, tol,
                            )?;
                            let before = straight.len();
                            if before == 0 && !crossed.is_empty() {
                                unthreaded = Some((groups.clone(), report.clone()));
                            }
                            straight.extend(crossed);
                            if straight.len() > before {
                                continue;
                            }
                            if let Some((was, then)) = unthreaded.take()
                                && (turned_over(
                                    model, &shape, &triangles, &groups, &built, &straight, tol,
                                )? || (refaceted && any_turned_in(model, &shape, tol)?))
                            {
                                groups = was;
                                report = then;
                                straight.clear();
                                threading_refused = true;
                                continue;
                            }
                        }
                        if culprits.is_empty() && options.recognize {
                            culprits = folded_seams(model, &groups, &built, tol)?;
                        }
                        if culprits.is_empty() && options.recognize {
                            culprits =
                                overlapping_faces(model, &shape, &adjacency, &groups, &built, tol)?;
                        }
                        if culprits.is_empty() && options.recognize {
                            culprits = unmatched_faces(model, &shape, &groups, &built, flat, tol)?;
                        }
                        if culprits.is_empty() {
                            report.faces = built.iter().flatten().count();
                            report.curved_faces = groups
                                .carriers
                                .iter()
                                .zip(&built)
                                .filter(|(c, b)| matches!(c, Carrier::Curved(_)) && b.is_some())
                                .count();
                            let patches: Vec<&Curved> = groups
                                .carriers
                                .iter()
                                .zip(&built)
                                .filter_map(|(c, b)| match (c, b) {
                                    (Carrier::Curved(curved), Some(_))
                                        if curved.patch.is_some() =>
                                    {
                                        Some(curved)
                                    }
                                    _ => None,
                                })
                                .collect();
                            report.patch_faces = patches.len();
                            report.patch_charts_mapped =
                                patches.iter().filter(|c| c.patch == Some(true)).count();
                            report.patches_not_disk = groups.refused.not_disk;
                            report.patches_narrow = groups.refused.narrow;
                            report.patches_unverified = groups.refused.unverified;
                            report.free_edges_fitted = plan.free_fitted;
                            break shape;
                        }
                        culprits
                    } else {
                        astray
                    }
                }
            };
            refaceted |= unthreaded.is_some() && !failed.is_empty();
            // A region that took facets gives them back and is tried again
            // without them before it is faceted.
            let returned: Vec<usize> = failed
                .iter()
                .copied()
                .filter(|g| absorbed.contains_key(g))
                .collect();
            if !returned.is_empty() {
                for g in returned {
                    if let Some(record) = absorbed.remove(&g) {
                        give_back(&mut groups, g, record);
                    }
                }
                continue;
            }
            for g in failed {
                groups.carriers[g] = Carrier::Gone;
                report.curved_faceted += 1;
                for of in &mut groups.of {
                    if *of == g {
                        *of = usize::MAX;
                    }
                }
            }
            coplanar_groups(
                &points,
                &triangles,
                &adjacency,
                options,
                flat,
                &mut groups,
                tol,
            )?;
        };
        if absorbing && any_turned_in(model, &shape, tol)? {
            absorbing = false;
            absorbed.clear();
            (groups, report) = unabsorbed.clone();
            continue 'attempt;
        }
        break shape;
    };
    Ok(MeshSolid {
        shape,
        closed: all_closed,
        coplanar_distance: flat,
        report,
    })
}

/// The recognized faces that overlap a face beside them: a fitted face
/// running past the facet it should end on encloses a sliver outside that
/// facet, and the facet then faces into material on its outer side while
/// every check on either face alone passes. Every recognized face and every
/// face beside one is probed as the validity check probes orientation; a
/// recognized face found turned in is a culprit, and for any other face
/// so found the recognized faces beside it are.
fn overlapping_faces(
    model: &Model,
    shape: &Shape,
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    let curved = |g: usize| {
        matches!(groups.carriers.get(g), Some(Carrier::Curved(_)))
            && built.get(g).is_some_and(Option::is_some)
    };
    let mut beside: HashMap<usize, std::collections::BTreeSet<usize>> = HashMap::new();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(t) = *twin else {
            continue;
        };
        let (mine, theirs) = (groups.of[h / 3], groups.of[t / 3]);
        if mine != theirs && curved(theirs) {
            beside.entry(mine).or_default().insert(theirs);
        }
    }
    // Each built face's box; a face a fitted face runs through need not
    // share an edge with it.
    let mut boxes: HashMap<usize, ogeom_math::Aabb> = HashMap::new();
    for (g, face) in built.iter().enumerate() {
        if let Some(face) = face {
            boxes.insert(
                g,
                crate::shape_bounds(model, face, tol)?.expanded(tol.confusion()),
            );
        }
    }
    let curved_boxes: Vec<(usize, ogeom_math::Aabb)> = boxes
        .iter()
        .filter(|(g, _)| curved(**g))
        .map(|(g, b)| (*g, *b))
        .collect();
    let crossing = |g: usize| -> Vec<usize> {
        boxes.get(&g).map_or_else(Vec::new, |own| {
            curved_boxes
                .iter()
                .filter(|(c, b)| *c != g && b.intersects(own))
                .map(|(c, _)| *c)
                .collect()
        })
    };
    let mut near: std::collections::BTreeSet<usize> = beside.keys().copied().collect();
    near.extend((0..groups.carriers.len()).filter(|&g| curved(g) || !crossing(g).is_empty()));
    if near.is_empty() {
        return Ok(Vec::new());
    }
    let mut of_face: HashMap<ogeom_topo::TShapeId, usize> = HashMap::new();
    let mut faces = Vec::new();
    for &g in &near {
        if let Some(face) = built.get(g).and_then(Option::as_ref) {
            of_face.insert(face.node(), g);
            faces.push(face.clone());
        }
    }
    let mut culprits = std::collections::BTreeSet::new();
    for solid in ogeom_topo::explore_unique(model, shape, ogeom_topo::ShapeType::Solid)? {
        for face in crate::check::faces_turned_in(model, &solid, &faces, tol)? {
            let Some(&g) = of_face.get(&face.node()) else {
                continue;
            };
            if curved(g) {
                culprits.insert(g);
            } else if let Some(next) = beside.get(&g) {
                culprits.extend(next.iter().copied());
            } else {
                culprits.extend(crossing(g));
            }
        }
    }
    Ok(culprits.into_iter().collect())
}

/// Whether any face of the shape's solids faces into their material.
fn any_turned_in(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    for solid in ogeom_topo::explore_unique(model, shape, ogeom_topo::ShapeType::Solid)? {
        if !crate::check::inside_out_faces(model, &solid, tol)?.is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether a face beside a seam threaded straight faces into the material
/// of the shape's solids, or has collapsed: thinner than the confusion
/// distance, its seams threaded onto one line.
fn turned_over(
    model: &Model,
    shape: &Shape,
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    straight: &std::collections::HashSet<(u32, u32)>,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mut beside: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for h in 0..triangles.len() * 3 {
        let (a, b) = from_to(triangles, h);
        if straight.contains(&(a.min(b), a.max(b))) {
            beside.insert(groups.of[h / 3]);
        }
    }
    let faces: Vec<Shape> = beside
        .iter()
        .filter_map(|&g| built.get(g).and_then(Option::as_ref).cloned())
        .collect();
    let fine = ogeom_mesh::Deflection::with_chord(ogeom_mesh::Deflection::default().chord * 1e-2)?;
    for face in &faces {
        let area = crate::surface_properties(model, face, fine, tol)?.mass;
        let reach = crate::shape_bounds(model, face, tol)?.diagonal();
        if area <= reach * tol.confusion() {
            return Ok(true);
        }
    }
    for solid in ogeom_topo::explore_unique(model, shape, ogeom_topo::ShapeType::Solid)? {
        if !crate::check::faces_turned_in(model, &solid, &faces, tol)?.is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The seams between a planar face and a curved one that bend the planar
/// face's trim back across itself: a facet thinner than the curve the two
/// surfaces meet along bulges, or a sliver between two curved faces whose
/// seams cross in its plane. Every face is meshed at the default deflection
/// on the shape's agreed edge chords, as a caller meshing the shape would.
/// A facet or two left between curved faces whose triangles cover a fifth
/// more than its exact area, or which draws a mesh edge no other face
/// draws (or that two others draw too), has its seams with curved faces
/// returned, to be threaded straight.
fn crossed_seams(
    model: &Model,
    shape: &Shape,
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<Vec<(u32, u32)>> {
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    // Only a facet or two left between curved faces: a larger planar face
    // holds its own shape, and straightening its seams can turn it over.
    let mut size: HashMap<usize, usize> = HashMap::new();
    for &g in &groups.of {
        *size.entry(g).or_insert(0) += 1;
    }
    // Each such face's mesh edges against a curved face.
    let mut seams: HashMap<usize, Vec<(u32, u32)>> = HashMap::new();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(g) = *twin else {
            continue;
        };
        let (mine, theirs) = (groups.of[h / 3], groups.of[g / 3]);
        if mine == theirs
            || !curved(theirs)
            || !matches!(groups.carriers.get(mine), Some(Carrier::Plane(_)))
            || size.get(&mine).copied().unwrap_or(0) > SLIVER_FACETS
        {
            continue;
        }
        let (a, b) = from_to(triangles, h);
        seams.entry(mine).or_default().push((a.min(b), a.max(b)));
    }
    if seams.is_empty() {
        return Ok(Vec::new());
    }
    let deflection = ogeom_mesh::Deflection::default();
    let fine = ogeom_mesh::Deflection::with_chord(deflection.chord * 1e-2)?;
    let chords = ogeom_mesh::edge_chords_for(model, shape, deflection, tol)?;
    let faces: Vec<(usize, &Shape)> = built
        .iter()
        .enumerate()
        .filter_map(|(g, b)| b.as_ref().map(|f| (g, f)))
        .collect();
    let meshes = ogeom_core::parallel::map_ordered(&faces, |_, &(_, face)| {
        ogeom_mesh::triangulate_face_with(model, face, deflection, &chords, tol).ok()
    });
    // The faces' points welded within the confusion distance: a shared
    // edge's points come from the same chords, and a seam's two columns
    // land a rounding apart.
    let cell = tol.confusion();
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    let cell_of = |p: Point| {
        (
            (p.x / cell).round() as i64,
            (p.y / cell).round() as i64,
            (p.z / cell).round() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    let mut welded: Vec<Point> = Vec::new();
    let mut weld = |p: Point| -> usize {
        let (x, y, z) = cell_of(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(found) = grid.get(&(x + dx, y + dy, z + dz)).and_then(|list| {
                        list.iter()
                            .copied()
                            .find(|&i| welded[i].distance(p) <= cell)
                    }) {
                        return found;
                    }
                }
            }
        }
        welded.push(p);
        grid.entry((x, y, z)).or_default().push(welded.len() - 1);
        welded.len() - 1
    };
    // The exact area of each meshed face with seams, which its drawn mesh
    // may overrun.
    let areas = ogeom_core::parallel::map_ordered(&faces, |i, &(g, face)| {
        (meshes[i].is_some() && seams.contains_key(&g))
            .then(|| crate::surface_properties(model, face, fine, tol))
    });
    let mut flagged: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut uses: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for ((&(g, _), mesh), area) in faces.iter().zip(&meshes).zip(areas) {
        let Some(mesh) = mesh else {
            if seams.contains_key(&g) {
                flagged.insert(g);
            }
            continue;
        };
        let at: Vec<usize> = mesh.positions.iter().map(|&p| weld(p)).collect();
        let mut drawn = 0.0;
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            drawn += (b - a).cross(c - a).magnitude() * 0.5;
            let corners = t.map(|i| at[i as usize]);
            if corners[0] == corners[1] || corners[1] == corners[2] || corners[2] == corners[0] {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (corners[k], corners[(k + 1) % 3]);
                uses.entry((a.min(b), a.max(b))).or_default().push(g);
            }
        }
        if let Some(area) = area
            && drawn > area?.mass * OVERRUN + tol.confusion()
        {
            flagged.insert(g);
        }
    }
    for users in uses.values() {
        if users.len() != 2 {
            flagged.extend(users.iter().copied().filter(|g| seams.contains_key(g)));
        }
    }
    let mut out: Vec<(u32, u32)> = flagged
        .iter()
        .flat_map(|g| seams[g].iter().copied())
        .collect();
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// The curved faces whose seams fold another face's boundary over itself
/// by more than those seams' tolerance.
///
/// Each built face's boundary is drawn in its surface's chart as the
/// tessellator trims it, and two of its chords that cross without sharing
/// an end are a fold: the face encloses some of its area twice, or none.
/// A fold no deeper than the tolerance of the edges that cross, or inside
/// the tolerance of a corner of the face, is within what those claim, and
/// stands. A deeper one names its two edges by
/// the points they were drawn through, and the curved faces across them
/// are the culprits; where neither is curved and the folded face is, it
/// is.
fn folded_seams(
    model: &Model,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    use ogeom_geom::Surface as _;
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    let deflection = ogeom_mesh::Deflection::default();
    let fine = ogeom_mesh::Deflection::with_chord(deflection.chord * 1e-2)?;
    // Each edge's faces and drawn points, by its node.
    let mut owners: HashMap<ogeom_topo::TShapeId, Vec<usize>> = HashMap::new();
    let mut drawn: HashMap<ogeom_topo::TShapeId, (Vec<Point>, f64)> = HashMap::new();
    for (g, face) in built.iter().enumerate() {
        let Some(face) = face else {
            continue;
        };
        for edge in ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)? {
            owners.entry(edge.node()).or_default().push(g);
            if let std::collections::hash_map::Entry::Vacant(slot) = drawn.entry(edge.node()) {
                let tolerance = model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .map_or(0.0, |d| d.tolerance.get());
                let polyline = ogeom_mesh::polyline_of_edge(model, &edge, fine, tol)?;
                slot.insert((polyline, tolerance));
            }
        }
    }
    let nearest = |at: Point, face_edges: &[ogeom_topo::TShapeId]| {
        face_edges
            .iter()
            .map(|id| {
                let d = drawn[id]
                    .0
                    .windows(2)
                    .map(|w| distance_to_segment(at, w[0], w[1]))
                    .fold(f64::INFINITY, f64::min);
                (d, *id)
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
    };
    let mut out: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for (g, face) in built.iter().enumerate() {
        let Some(face) = face else {
            continue;
        };
        let Some(surface) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        else {
            continue;
        };
        // The face's own edges; a seam bounds the face twice, up one side
        // of the chart and down the other, and is no fold.
        let Some(surface_id) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            continue;
        };
        let face_edges: Vec<ogeom_topo::TShapeId> =
            ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)?
                .iter()
                .filter(|e| {
                    !model
                        .node(e)
                        .and_then(|n| n.data().as_edge())
                        .and_then(|d| d.pcurve_for(surface_id, e.location()))
                        .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
                })
                .map(Shape::node)
                .collect();
        // The face's corners and how far each claims.
        let balls: Vec<(Point, f64)> =
            ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Vertex)?
                .iter()
                .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
                .map(|v| (v.point, v.tolerance.get()))
                .collect();
        // Drawn as the tessellator draws it first; a face whose chords
        // cross there is drawn a hundred times finer, where chords that only
        // crossed for their sag (two long arcs meeting at a sharp corner)
        // part, and what still crosses is the boundary's own fold.
        if chart_crossings(&ogeom_mesh::face_boundary(model, face, deflection, tol)?).is_empty() {
            continue;
        }
        let rings = ogeom_mesh::face_boundary(model, face, fine, tol)?;
        for [a, b, p, q] in chart_crossings(&rings) {
            let lift = |x: Point2| surface.point_at(x.x, x.y, tol);
            // Where the chords cross, along the first.
            let side = |o: Point2, x: Point2, y: Point2| {
                (x.x - o.x).mul_add(y.y - o.y, -((x.y - o.y) * (y.x - o.x)))
            };
            let (da, db) = (side(p, q, a), side(p, q, b));
            let t = da / (da - db);
            let Ok(at) = lift(Point2::new(
                (b.x - a.x).mul_add(t, a.x),
                (b.y - a.y).mul_add(t, a.y),
            )) else {
                continue;
            };
            // A crossing inside a corner's tolerance is the corner's slack:
            // its edges' ends stand off it by up to that much.
            if balls.iter().any(|&(v, r)| at.distance(v) <= r) {
                continue;
            }
            let (Ok(a), Ok(b), Ok(p), Ok(q)) = (lift(a), lift(b), lift(p), lift(q)) else {
                continue;
            };
            let depth = distance_to_line(a, p, q)
                .min(distance_to_line(b, p, q))
                .min(distance_to_line(p, a, b))
                .min(distance_to_line(q, a, b));
            // Each chord named by the edge it was drawn from; a chord no
            // edge was drawn through (a side of the chart's window, or a
            // seam) folds nothing.
            let (Some((d1, first)), Some((d2, second))) = (
                nearest(a.lerp(b, 0.5), &face_edges),
                nearest(p.lerp(q, 0.5), &face_edges),
            ) else {
                continue;
            };
            let drawn_within = fine.chord * 2.0;
            if d1 > drawn[&first].1 + drawn_within || d2 > drawn[&second].1 + drawn_within {
                continue;
            }
            if depth <= drawn[&first].1.max(drawn[&second].1) {
                continue;
            }
            let across: Vec<usize> = [first, second]
                .iter()
                .flat_map(|id| owners[id].iter().copied())
                .filter(|&o| o != g && curved(o))
                .collect();
            if across.is_empty() {
                if curved(g) {
                    out.insert(g);
                }
            } else {
                out.extend(across);
            }
        }
    }
    Ok(out.into_iter().collect())
}

/// Each pair of a boundary's chords that cross without sharing an end, as
/// the two chords' ends.
fn chart_crossings(rings: &[Vec<Point2>]) -> Vec<[Point2; 4]> {
    // Each chord, with its ring and the index of its first point, swept in
    // order of its least `u`.
    let mut chords: Vec<(usize, usize, Point2, Point2)> = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        for i in 0..ring.len() {
            chords.push((r, i, ring[i], ring[(i + 1) % ring.len()]));
        }
    }
    chords.sort_by(|a, b| a.2.x.min(a.3.x).total_cmp(&b.2.x.min(b.3.x)));
    let mut out = Vec::new();
    let mut active: Vec<usize> = Vec::new();
    for k in 0..chords.len() {
        let (r, i, a, b) = chords[k];
        let low = a.x.min(b.x);
        active.retain(|&j| chords[j].2.x.max(chords[j].3.x) >= low);
        for &j in &active {
            let (rj, ij, p, q) = chords[j];
            if r == rj {
                let n = rings[r].len();
                if (i + 1) % n == ij || (ij + 1) % n == i {
                    continue;
                }
            }
            if a.y.max(b.y) < p.y.min(q.y) || p.y.max(q.y) < a.y.min(b.y) {
                continue;
            }
            if segments_cross((a.x, a.y), (b.x, b.y), (p.x, p.y), (q.x, q.y)) {
                out.push([a, b, p, q]);
            }
        }
        active.push(k);
    }
    out
}

/// Facets given to a curved region, and what it was before, so they can
/// be handed back.
struct Absorbed {
    /// Each facet: its own group, its plane, its triangles.
    facets: Vec<(usize, Carrier, Vec<usize>)>,
    /// The region's vertices before it took any.
    vertices: Vec<u32>,
}

/// Give each facet beside a curved region to that region.
///
/// Where curved regions meet at an angle, or where a fit stops a row short,
/// the mesh leaves planar facets of a triangle or two that no surface
/// claimed, between curved faces. Built as faces of their own they are
/// bounded by seams whose tolerance (a chord threaded through the mesh's
/// vertices strays up to a twentieth of its span) is as wide as they are,
/// so their trims fold and they mesh over themselves. A facet of at most
/// [`SLIVER_FACETS`] triangles whose corners all lie within the reach of a
/// curved neighbour's surface goes to the neighbour they lie nearest. Such
/// a facet has no vertex inside it: each of its corners ends on a seam,
/// held to the reach as any seam's points are, so the face it joins is
/// held to what it was. What each region took is returned, so a region
/// that cannot be built with its facets gives them back before it is
/// faceted itself.
fn absorb_facets(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    flat: f64,
) -> HashMap<usize, Absorbed> {
    let reach = flat * REACH;
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Plane(_))) {
            members.entry(g).or_default().push(t);
        }
    }
    let mut facets: Vec<(usize, Vec<usize>)> = members
        .into_iter()
        .filter(|(_, ts)| ts.len() <= SLIVER_FACETS)
        .collect();
    facets.sort_unstable();
    let mut absorbed: HashMap<usize, Absorbed> = HashMap::new();
    for (g, ts) in facets {
        let mut beside: Vec<usize> = ts
            .iter()
            .flat_map(|&t| {
                (3 * t..3 * t + 3).filter_map(|h| adjacency.twin[h].map(|o| groups.of[o / 3]))
            })
            .filter(|&o| o != g && matches!(groups.carriers.get(o), Some(Carrier::Curved(_))))
            .collect();
        beside.sort_unstable();
        beside.dedup();
        let mut corners: Vec<u32> = ts.iter().flat_map(|&t| triangles[t]).collect();
        corners.sort_unstable();
        corners.dedup();
        let at: Vec<Point> = corners.iter().map(|&v| points[v as usize]).collect();
        let nearest = beside
            .iter()
            .filter_map(|&o| match &groups.carriers[o] {
                Carrier::Curved(c) => Some((worst_deviation(&c.shape, &at), o)),
                _ => None,
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        let Some((deviation, to)) = nearest else {
            continue;
        };
        if deviation > reach {
            continue;
        }
        let Carrier::Curved(region) = &mut groups.carriers[to] else {
            continue;
        };
        let record = absorbed.entry(to).or_insert_with(|| Absorbed {
            facets: Vec::new(),
            vertices: region.vertices.clone(),
        });
        region.vertices.extend(corners);
        region.vertices.sort_unstable();
        region.vertices.dedup();
        for &t in &ts {
            groups.of[t] = to;
        }
        let carrier = std::mem::replace(&mut groups.carriers[g], Carrier::Gone);
        record.facets.push((g, carrier, ts));
    }
    absorbed
}

/// Hand a region's facets back to themselves, as they were before it took
/// them.
fn give_back(groups: &mut Groups, to: usize, record: Absorbed) {
    for (g, carrier, ts) in record.facets {
        for t in ts {
            groups.of[t] = g;
        }
        groups.carriers[g] = carrier;
    }
    if let Carrier::Curved(region) = &mut groups.carriers[to] {
        region.vertices = record.vertices;
    }
}

/// The curved faces whose meshes do not meet their neighbours', where the
/// solid does not mesh closed.
///
/// Each built face is meshed on its own, its edges drawn alike for every
/// face, and the meshes are joined where their points coincide. A mesh
/// edge used other than twice is a gap or an overlap: a face whose trim
/// folds within its seams' tolerance, which every check on the solid
/// passes, still meshes over itself. The curved faces among those using
/// such an edge are the culprits; where only planes use it, the curved
/// faces with an edge drawn within the reach of it are.
fn unmatched_faces(
    model: &Model,
    shape: &Shape,
    groups: &Groups,
    built: &[Option<Shape>],
    flat: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    type Key = (u64, u64, u64);
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    let deflection = ogeom_mesh::Deflection::default();
    // What the solid is drawn as decides: a crack between two faces' own
    // meshes that the drawing welds shut costs nothing.
    let (drawn, chords) = ogeom_mesh::triangulate_with_chords(model, shape, deflection, tol)?;
    if drawn.is_closed() {
        return Ok(Vec::new());
    }
    let faces: Vec<(usize, &Shape)> = built
        .iter()
        .enumerate()
        .filter_map(|(g, b)| b.as_ref().map(|f| (g, f)))
        .collect();
    let meshes = ogeom_core::parallel::map_ordered(&faces, |_, &(_, face)| {
        ogeom_mesh::triangulate_face_with(model, face, deflection, &chords, tol).ok()
    });
    // The faces' points welded within the confusion distance: a shared
    // edge's points come from the same chords, but a seam drawn from either
    // side of a closed chart lands a rounding apart.
    let cell = tol.confusion();
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    let cell_of = |p: Point| -> Key {
        (
            (p.x / cell).round() as i64 as u64,
            (p.y / cell).round() as i64 as u64,
            (p.z / cell).round() as i64 as u64,
        )
    };
    let mut grid: HashMap<Key, Vec<usize>> = HashMap::new();
    let mut welded: Vec<Point> = Vec::new();
    let mut weld = |p: Point| -> Key {
        let c = cell_of(p);
        for dx in [0u64, 1, u64::MAX] {
            for dy in [0u64, 1, u64::MAX] {
                for dz in [0u64, 1, u64::MAX] {
                    let near = (
                        c.0.wrapping_add(dx),
                        c.1.wrapping_add(dy),
                        c.2.wrapping_add(dz),
                    );
                    if let Some(&i) = grid
                        .get(&near)
                        .and_then(|list| list.iter().find(|&&i| welded[i].distance(p) <= cell))
                    {
                        return cell_of(welded[i]);
                    }
                }
            }
        }
        welded.push(p);
        grid.entry(c).or_default().push(welded.len() - 1);
        c
    };
    let mut uses: HashMap<(Key, Key), (Vec<usize>, [Point; 2])> = HashMap::new();
    let mut out: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for (&(g, _), mesh) in faces.iter().zip(&meshes) {
        let Some(mesh) = mesh else {
            if curved(g) {
                out.insert(g);
            }
            continue;
        };
        for t in &mesh.triangles {
            for k in 0..3 {
                let (p, q) = (
                    mesh.positions[t[k] as usize],
                    mesh.positions[t[(k + 1) % 3] as usize],
                );
                let (a, b) = (weld(p), weld(q));
                if a == b {
                    continue;
                }
                let entry = uses
                    .entry((a.min(b), a.max(b)))
                    .or_insert_with(|| (Vec::new(), [p, q]));
                entry.0.push(g);
            }
        }
    }
    let free = free_edges(model, &faces, tol)?;
    let mut planar: Vec<Point> = Vec::new();
    for (users, ends) in uses.values() {
        // One face's own seam, drawn from both sides of its closed chart,
        // meets itself there, and a pole's fans with it: an even count
        // from that face alone is its own closure, not a gap.
        if users.len() == 2 || (users.len() % 2 == 0 && users.iter().all(|&g| g == users[0])) {
            continue;
        }
        // A free edge (one face, nothing across) is drawn by that face
        // alone.
        if users.len() == 1
            && free.iter().any(|(curve, reach, (lo, hi))| {
                let inside = |p: Point| {
                    (lo.x..=hi.x).contains(&p.x)
                        && (lo.y..=hi.y).contains(&p.y)
                        && (lo.z..=hi.z).contains(&p.z)
                };
                ends.iter().all(|&p| {
                    inside(p)
                        && crate::measure::project_on_curve(curve, p, 64, tol)
                            .is_ok_and(|near| near.distance <= *reach)
                })
            })
        {
            continue;
        }
        let at = ends[0].lerp(ends[1], 0.5);
        let mine: Vec<usize> = users.iter().copied().filter(|&g| curved(g)).collect();
        if mine.is_empty() {
            planar.push(at);
        } else {
            out.extend(mine);
        }
    }
    if !planar.is_empty() {
        let reach = flat * REACH;
        for &(g, face) in &faces {
            if !curved(g) || out.contains(&g) {
                continue;
            }
            'edges: for edge in
                ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)?
            {
                let drawn = ogeom_mesh::polyline_of_edge(model, &edge, deflection, tol)?;
                for w in drawn.windows(2) {
                    if planar
                        .iter()
                        .any(|&p| distance_to_segment(p, w[0], w[1]) <= reach)
                    {
                        out.insert(g);
                        break 'edges;
                    }
                }
            }
        }
    }
    Ok(out.into_iter().collect())
}

/// A free edge as [`unmatched_faces`] reads it: its trimmed curve, how
/// far from it a point counts as on it, and the box holding every such
/// point.
type FreeEdge = (Curve, f64, (Point, Point));

/// The edges of the built faces that bound one face only and are no seam
/// of it: each one's curve, trimmed to its range, its tolerance with the
/// confusion distance's margin, and a box holding every point within that
/// of it (its samples' box, widened by the tolerance and the longest step
/// between samples).
fn free_edges(
    model: &Model,
    faces: &[(usize, &Shape)],
    tol: Tolerances,
) -> OgeomResult<Vec<FreeEdge>> {
    use ogeom_geom::Curve3d as _;
    let mut count: HashMap<ogeom_topo::SameKey, usize> = HashMap::new();
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
const SLIVER_FACETS: usize = 2;

/// How much more a planar face's drawn triangles may cover than its exact
/// area before its trim is taken to cross itself: a planar face's mesh
/// covers its area exactly, but for the rounding of its boundary's chords.
const OVERRUN: f64 = 1.2;

/// The recognized regions whose built face strays from the triangles it
/// replaces: its bounds stand beyond theirs by more than the surface bulges
/// past its facets and its boundary can run on to meet its neighbours, or
/// some point of it lies away from all of them.
fn astray_faces(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    flat: f64,
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
            let drawn = ogeom_mesh::triangulate_face(
                model,
                face,
                ogeom_mesh::Deflection::with_chord(chord)?,
                tol,
            )
            .or_else(|_| {
                ogeom_mesh::triangulate_face(model, face, ogeom_mesh::Deflection::default(), tol)
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
                let Ok(drawn) = ogeom_mesh::triangulate_face(
                    model,
                    beside,
                    ogeom_mesh::Deflection::default(),
                    tol,
                ) else {
                    continue;
                };
                let middles: Vec<Point> = drawn
                    .triangles
                    .iter()
                    .map(|t| {
                        let [a, b, c] = t.map(|i| drawn.positions[i as usize]);
                        Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0)
                    })
                    .filter(|&m| reaches(m))
                    .collect();
                if middles.is_empty() {
                    continue;
                }
                let own: Vec<usize> = (0..triangles.len())
                    .filter(|&t| groups.of[t] == plane)
                    .collect();
                let samples = middles.len().min(32);
                let off = (0..samples)
                    .filter(|&k| {
                        let m = middles[k * middles.len() / samples];
                        near.iter().chain(&own).all(|&t| {
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

/// One closed body of the result: the solid, and the mesh triangles of its
/// shell and of the voids inside it.
struct Body {
    solid: Shape,
    triangles: Vec<u32>,
}

/// Faces into one shell per piece; closed pieces into solids, a piece
/// nested an odd number of times deep being a void of the one around it.
/// The result, and each solid with its triangles.
#[allow(
    clippy::too_many_arguments,
    reason = "the build's own state, passed through"
)]
fn assemble(
    model: &mut Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    pieces: &[Piece],
    depth: &[usize],
    all_closed: bool,
    groups: &Groups,
    built: &[Option<Shape>],
) -> OgeomResult<(Shape, Vec<Body>)> {
    let mut shells = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let mut faces: Vec<Shape> = Vec::new();
        let mut taken = vec![false; groups.carriers.len()];
        for &t in &piece.triangles {
            let g = groups.of[t as usize];
            if !taken[g] {
                taken[g] = true;
                if let Some(face) = &built[g] {
                    faces.push(face.clone());
                }
            }
        }
        shells.push(model.add_shell(&faces)?);
    }
    if !all_closed {
        let shape = if shells.len() == 1 {
            shells.swap_remove(0)
        } else {
            model.add_compound(&shells)?
        };
        return Ok((shape, Vec::new()));
    }
    let mut bodies = Vec::new();
    for (i, shell) in shells.iter().enumerate() {
        if depth[i] % 2 == 1 {
            continue;
        }
        let mut members = vec![shell.clone()];
        let mut mine = pieces[i].triangles.clone();
        for (j, void) in shells.iter().enumerate() {
            if depth[j] == depth[i] + 1 && inside(points, triangles, &pieces[i], &pieces[j]) {
                members.push(void.clone());
                mine.extend(&pieces[j].triangles);
            }
        }
        bodies.push(Body {
            solid: model.add_solid(&members)?,
            triangles: mine,
        });
    }
    let shape = if bodies.len() == 1 {
        bodies[0].solid.clone()
    } else {
        let solids: Vec<Shape> = bodies.iter().map(|b| b.solid.clone()).collect();
        model.add_compound(&solids)?
    };
    Ok((shape, bodies))
}

/// The recognized regions to facet so every body of the result is sound: a
/// recognized face turned into the material gives up its region, and a
/// body whose volume is not
/// the mesh's, within what its recognized surfaces may add over their
/// facets, gives up all of its recognized regions (a face closed over the
/// wrong part of its surface passes every local test and is caught only
/// there). Other bodies keep theirs.
#[allow(
    clippy::too_many_arguments,
    reason = "the build's own state, passed through"
)]
fn body_culprits(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    bodies: &[Body],
    flat: f64,
    tol: Tolerances,
    report: &mut MeshSolidReport,
) -> OgeomResult<Vec<usize>> {
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    let own = |body: &Body| -> Vec<usize> {
        let mut own: Vec<usize> = body
            .triangles
            .iter()
            .map(|&t| groups.of[t as usize])
            .filter(|&g| curved(g))
            .collect();
        own.sort_unstable();
        own.dedup();
        own
    };
    // Orientation first: a recognized face turned against the triangles it
    // replaces faces into the material, and is facetted whatever the
    // volumes say.
    let culprits = inverted_faces(model, points, triangles, groups, built, tol)?;
    if !culprits.is_empty() {
        return Ok(culprits);
    }
    let mut culprits: Vec<usize> = Vec::new();
    for body in bodies {
        let own = own(body);
        if own.is_empty() {
            continue;
        }
        let (mesh_volume, allowance, area) =
            volume_allowance(points, triangles, &body.triangles, groups, flat);
        let diagonal = body_diagonal(points, triangles, &body.triangles);
        // Measured no finer than the allowance needs: at the facets' mean
        // offset from their surfaces, whose error over the area is the
        // allowance again, counted in the slack.
        let chord = (allowance / area.max(f64::MIN_POSITIVE)).max(flat);
        let deflection = ogeom_mesh::Deflection::with_chord(chord)?;
        // A body with a swept face is measured on its tessellation at that
        // chord: the exact integral over a fitted profile is held to far
        // finer than the slack asks, at far greater cost.
        let swept = own
            .iter()
            .any(|&g| matches!(&groups.carriers[g], Carrier::Curved(c) if matches!(c.shape, Canonical::Swept(_))));
        let measured = if swept {
            let mesh = ogeom_mesh::triangulate(model, &body.solid, deflection, tol)?;
            mesh.triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| mesh.positions[i as usize].to_vector());
                    a.dot(b.cross(c)) / 6.0
                })
                .sum::<f64>()
        } else {
            crate::volume_properties(model, &body.solid, deflection, tol)?.mass
        };
        // Measured at the facets' corners and edge middles, the allowance
        // misses the surface's rise inside a facet and a fitted boundary's
        // wander between the rows; twice over covers both.
        let slack = allowance * 2.0 + area * chord * 2.0 + diagonal.powi(3) * 1e-12;
        if (measured - mesh_volume).abs() > slack {
            // The sweeps answer first: the least proven of the fits. Only
            // where none is left does the body give up its recognition.
            let sweeps: Vec<usize> = own
                .iter()
                .copied()
                .filter(|&g| matches!(&groups.carriers[g], Carrier::Curved(c) if matches!(c.shape, Canonical::Swept(_))))
                .collect();
            if sweeps.is_empty() {
                report.recognition_withdrawn = true;
                culprits.extend(own);
            } else {
                culprits.extend(sweeps);
            }
        }
    }
    culprits.sort_unstable();
    culprits.dedup();
    Ok(culprits)
}

/// The recognized regions whose built face points against the triangles it
/// replaces. The mesh is oriented outward before anything is built, so each
/// triangle's normal is the side the material is not on; the face's own
/// outward normal, at the foot of the triangle's middle on its surface,
/// must agree with it. A few triangles vote, and the face is turned only
/// where most of them say so.
fn inverted_faces(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    use ogeom_geom::Surface as _;
    const VOTES: usize = 9;
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Curved(_))) {
            members[g].push(t);
        }
    }
    let mut inverted = Vec::new();
    for (g, tris) in members.iter().enumerate() {
        let Some(face) = built.get(g).and_then(Option::as_ref) else {
            continue;
        };
        if tris.is_empty() {
            continue;
        }
        let Some(surface) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|data| model.geometry().surface(data.surface))
        else {
            continue;
        };
        let turned = face.orientation() == ogeom_topo::Orientation::Reversed;
        let votes = VOTES.min(tris.len());
        let (mut against, mut asked) = (0, 0);
        for k in 0..votes {
            let [a, b, c] = triangles[tris[k * tris.len() / votes]].map(|v| points[v as usize]);
            let normal = (b - a).cross(c - a);
            let middle = Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0);
            let Ok(foot) = crate::measure::project_on_surface(surface, middle, 8, tol) else {
                continue;
            };
            let Ok(outward) = surface.normal_at(foot.parameters.0, foot.parameters.1, tol) else {
                continue;
            };
            let outward = if turned {
                -outward.vector()
            } else {
                outward.vector()
            };
            asked += 1;
            if outward.dot(normal) < 0.0 {
                against += 1;
            }
        }
        if asked > 0 && against * 2 > asked {
            inverted.push(g);
        }
    }
    Ok(inverted)
}

/// The diagonal of the box around some of the mesh's triangles.
fn body_diagonal(points: &[Point], triangles: &[[u32; 3]], mine: &[u32]) -> f64 {
    let mut lo = Point::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Point::new(f64::MIN, f64::MIN, f64::MIN);
    for &t in mine {
        for &v in &triangles[t as usize] {
            let p = points[v as usize];
            lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    lo.distance(hi)
}

/// Some triangles' own volume, how far the solid built on their groups may
/// differ from it, and their area. Each triangle on a recognized surface may
/// stand off it by its worst corner or edge middle, over its whole area; a
/// triangle on a plane by the coplanar distance.
fn volume_allowance(
    points: &[Point],
    triangles: &[[u32; 3]],
    mine: &[u32],
    groups: &Groups,
    flat: f64,
) -> (f64, f64, f64) {
    let (mut volume, mut allowance, mut area) = (0.0, 0.0, 0.0);
    for &t in mine {
        let t = t as usize;
        let tri = &triangles[t];
        let [a, b, c] = tri.map(|v| points[v as usize]);
        let (va, vb, vc) = (a - Point::ORIGIN, b - Point::ORIGIN, c - Point::ORIGIN);
        volume += va.dot(vb.cross(vc)) / 6.0;
        let size = (b - a).cross(c - a).magnitude() / 2.0;
        area += size;
        let off = match groups.carriers.get(groups.of[t]) {
            Some(Carrier::Curved(curved)) => {
                let middle = |p: Point, q: Point| p + (q - p) * 0.5;
                [
                    a,
                    b,
                    c,
                    middle(a, b),
                    middle(b, c),
                    middle(c, a),
                    middle(a, middle(b, c)),
                ]
                .into_iter()
                .map(|p| curved.shape.distance_to(p))
                .fold(0.0_f64, f64::max)
            }
            _ => flat,
        };
        allowance += size * off;
    }
    (volume, allowance, area)
}

/// Weld points on a grid of the weld distance, looking in the neighbouring
/// cells too, so two points either side of a cell wall still meet.
fn weld_points(positions: &[Point], weld: f64) -> (Vec<Point>, Vec<u32>) {
    // A cast saturates, so a coordinate past `i64`'s cells shares the last
    // one and is still compared by distance.
    #[allow(clippy::cast_possible_truncation, reason = "saturating")]
    let cell = |p: Point| {
        (
            (p.x / weld).floor() as i64,
            (p.y / weld).floor() as i64,
            (p.z / weld).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::with_capacity(positions.len());
    let mut kept: Vec<Point> = Vec::with_capacity(positions.len());
    let mut remap = Vec::with_capacity(positions.len());
    for p in positions {
        let (x, y, z) = cell(*p);
        let mut found = None;
        'search: for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let key = (
                        x.saturating_add(dx),
                        y.saturating_add(dy),
                        z.saturating_add(dz),
                    );
                    if let Some(list) = grid.get(&key)
                        && let Some(&k) = list
                            .iter()
                            .find(|&&k| kept[k as usize].distance(*p) <= weld)
                    {
                        found = Some(k);
                        break 'search;
                    }
                }
            }
        }
        let index = found.unwrap_or_else(|| {
            let k = u32::try_from(kept.len()).unwrap_or(u32::MAX);
            kept.push(*p);
            grid.entry((x, y, z)).or_default().push(k);
            k
        });
        remap.push(index);
    }
    (kept, remap)
}

/// Whether a triangle stands more than the weld distance tall over its
/// longest side.
fn has_area(points: &[Point], [a, b, c]: [u32; 3], weld: f64) -> bool {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    let twice_area = (b - a).cross(c - a).magnitude();
    let longest = a.distance(b).max(b.distance(c)).max(c.distance(a));
    twice_area > weld * longest
}

fn diagonal(points: &[Point]) -> f64 {
    let (lo, hi) = points
        .iter()
        .fold(([f64::MAX; 3], [f64::MIN; 3]), |(lo, hi), p| {
            (
                [lo[0].min(p.x), lo[1].min(p.y), lo[2].min(p.z)],
                [hi[0].max(p.x), hi[1].max(p.y), hi[2].max(p.z)],
            )
        });
    Point::new(lo[0], lo[1], lo[2]).distance(Point::new(hi[0], hi[1], hi[2]))
}

/// A half-edge: side `k` of triangle `t`, as `3 t + k`, running from the
/// triangle's corner `k` to corner `k + 1`.
type Half = usize;

/// One way to unfold a sliver: the neighbour across the swapped edge, the
/// two triangles that replace the pair, the edge given up and the one taken,
/// and how well shaped the smaller of the two is.
struct Unfolding {
    shape: f64,
    neighbour: usize,
    pair: [[u32; 3]; 2],
    old: (u32, u32),
    new: (u32, u32),
}

/// Swap the diagonal under each fold of the mesh: a sliver facing against
/// all three of its neighbours, which agree among themselves, is the
/// surface folded back over itself, as an exporter leaves where it moved a
/// vertex across a thin triangle. Consistently wound, the fold survives
/// orientation and becomes a face pointing into the material. Across one
/// of its edges the sliver and its neighbour make a quadrilateral whose
/// other diagonal gives two triangles facing with the neighbours; of the
/// diagonals that do and are not already edges of the mesh, the one whose
/// smaller triangle is largest is taken. Returns how many folds were
/// swapped.
fn unfold(points: &[Point], triangles: &mut [[u32; 3]], adjacency: &Adjacency) -> usize {
    let normal = |t: [u32; 3]| -> Vector {
        let [a, b, c] = t.map(|v| points[v as usize]);
        (b - a).cross(c - a)
    };
    let unit = |v: Vector| {
        let m = v.magnitude();
        if m > 0.0 { v / m } else { v }
    };
    let mut edges: std::collections::HashSet<(u32, u32)> = triangles
        .iter()
        .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
        .map(|(a, b)| (a.min(b), a.max(b)))
        .collect();
    let mut touched = vec![false; triangles.len()];
    let mut swapped = 0;
    for t in 0..triangles.len() {
        if touched[t] {
            continue;
        }
        let Some(twins) = (0..3)
            .map(|k| adjacency.twin[3 * t + k])
            .collect::<Option<Vec<Half>>>()
        else {
            continue;
        };
        if twins.iter().any(|g| touched[g / 3]) {
            continue;
        }
        let own = unit(normal(triangles[t]));
        let around: Vec<Vector> = twins
            .iter()
            .map(|g| unit(normal(triangles[g / 3])))
            .collect();
        let folded = around.iter().all(|n| n.dot(own) < -0.5)
            && around
                .iter()
                .enumerate()
                .all(|(i, n)| around[i + 1..].iter().all(|m| m.dot(*n) > 0.5));
        if !folded {
            continue;
        }
        let facing = around
            .iter()
            .fold(Vector::new(0.0, 0.0, 0.0), |acc, n| acc + *n);
        // Of the diagonals that unfold it, the one whose smaller triangle
        // is the larger: the best-shaped pair.
        let mut best: Option<Unfolding> = None;
        for (k, &g) in twins.iter().enumerate() {
            let u = g / 3;
            let (a, b) = from_to(triangles, 3 * t + k);
            let c = triangles[t][(k + 2) % 3];
            let d = triangles[u][(g % 3 + 2) % 3];
            if c == d || edges.contains(&(c.min(d), c.max(d))) {
                continue;
            }
            let (first, second) = ([c, a, d], [d, b, c]);
            let (n1, n2) = (normal(first), normal(second));
            if n1.dot(facing) <= 0.0 || n2.dot(facing) <= 0.0 {
                continue;
            }
            let shape = n1.magnitude().min(n2.magnitude());
            if best.as_ref().is_none_or(|b| shape > b.shape) {
                best = Some(Unfolding {
                    shape,
                    neighbour: u,
                    pair: [first, second],
                    old: (a.min(b), a.max(b)),
                    new: (c.min(d), c.max(d)),
                });
            }
        }
        if let Some(swap) = best {
            edges.remove(&swap.old);
            edges.insert(swap.new);
            triangles[t] = swap.pair[0];
            triangles[swap.neighbour] = swap.pair[1];
            touched[t] = true;
            touched[swap.neighbour] = true;
            swapped += 1;
        }
    }
    swapped
}

fn from_to(triangles: &[[u32; 3]], h: Half) -> (u32, u32) {
    let t = triangles[h / 3];
    (t[h % 3], t[(h % 3 + 1) % 3])
}

fn next(h: Half) -> Half {
    h - h % 3 + (h % 3 + 1) % 3
}

/// Which half-edges pair across a mesh edge exactly two triangles share.
struct Adjacency {
    /// The other triangle's half-edge on the same mesh edge, where exactly
    /// one other triangle uses it.
    twin: Vec<Option<Half>>,
    used_once: usize,
    used_more: usize,
}

impl Adjacency {
    fn new(triangles: &[[u32; 3]]) -> Self {
        let mut keyed: Vec<(u64, Half)> = (0..triangles.len() * 3)
            .map(|h| {
                let (a, b) = from_to(triangles, h);
                ((u64::from(a.min(b)) << 32) | u64::from(a.max(b)), h)
            })
            .collect();
        keyed.sort_unstable();
        let mut twin = vec![None; keyed.len()];
        let (mut used_once, mut used_more) = (0, 0);
        let mut crowded: Vec<std::ops::Range<usize>> = Vec::new();
        let mut i = 0;
        while i < keyed.len() {
            let mut j = i;
            while j < keyed.len() && keyed[j].0 == keyed[i].0 {
                j += 1;
            }
            let n = j - i;
            match n {
                1 => used_once += 1,
                2 => {
                    twin[keyed[i].1] = Some(keyed[i + 1].1);
                    twin[keyed[i + 1].1] = Some(keyed[i].1);
                }
                _ => {
                    used_more += 1;
                    crowded.push(i..j);
                }
            }
            i = j;
        }
        // An edge more than two triangles use is where bodies meet (two
        // blocks sharing an edge, as an exporter writes glued parts). The
        // bodies are what the edges used exactly twice join; at a crowded
        // edge each body's own two triangles are each other's twins, and each
        // body closes on its own.
        if !crowded.is_empty() {
            let mut body: Vec<usize> = (0..triangles.len()).collect();
            fn root(body: &mut [usize], mut t: usize) -> usize {
                while body[t] != t {
                    body[t] = body[body[t]];
                    t = body[t];
                }
                t
            }
            for (h, g) in twin.iter().enumerate() {
                if let Some(g) = g {
                    let (a, b) = (root(&mut body, h / 3), root(&mut body, g / 3));
                    body[a.max(b)] = a.min(b);
                }
            }
            for range in crowded {
                let halves: Vec<Half> = keyed[range].iter().map(|&(_, h)| h).collect();
                let bodies: Vec<usize> = halves.iter().map(|&h| root(&mut body, h / 3)).collect();
                for (k, &h) in halves.iter().enumerate() {
                    let mine: Vec<usize> = (0..halves.len())
                        .filter(|&m| bodies[m] == bodies[k])
                        .collect();
                    if let [x, y] = mine[..]
                        && x == k
                        && halves[x] / 3 != halves[y] / 3
                    {
                        twin[h] = Some(halves[y]);
                        twin[halves[y]] = Some(h);
                    }
                }
            }
        }
        Self {
            twin,
            used_once,
            used_more,
        }
    }
}

/// One connected piece: its triangles, and whether it closes.
struct Piece {
    triangles: Vec<u32>,
    closed: bool,
}

/// Make windings agree across every shared edge, piece by piece; turn a
/// closed piece outward and leave an open one the way most of its
/// triangles already were.
fn orient(
    points: &[Point],
    triangles: &mut [[u32; 3]],
    adjacency: &Adjacency,
    report: &mut MeshSolidReport,
) -> Vec<Piece> {
    let n = triangles.len();
    let mut flip: Vec<Option<bool>> = vec![None; n];
    let mut pieces = Vec::new();
    for seed in 0..n {
        if flip[seed].is_some() {
            continue;
        }
        flip[seed] = Some(false);
        let mut members = vec![seed];
        let mut stack = vec![seed];
        let mut conflicts = 0;
        let mut open = false;
        while let Some(t) = stack.pop() {
            let mine = flip[t].unwrap_or(false);
            for h in 3 * t..3 * t + 3 {
                let Some(g) = adjacency.twin[h] else {
                    open = true;
                    continue;
                };
                let other = g / 3;
                // Agreeing windings run a shared edge opposite ways.
                let same_way = from_to(triangles, h) == from_to(triangles, g);
                let wanted = mine ^ same_way;
                match flip[other] {
                    None => {
                        flip[other] = Some(wanted);
                        members.push(other);
                        stack.push(other);
                    }
                    Some(have) if have != wanted => conflicts += 1,
                    Some(_) => {}
                }
            }
        }
        // Each conflicting edge was met from both sides.
        conflicts /= 2;
        report.orientation_conflicts += conflicts;
        let flipped = members.iter().filter(|&&t| flip[t] == Some(true)).count();
        let closed = !open && conflicts == 0;
        let turn_all = if closed {
            let volume: f64 = members
                .iter()
                .map(|&t| {
                    let mut tri = triangles[t];
                    if flip[t] == Some(true) {
                        tri.swap(1, 2);
                    }
                    signed_volume(points, tri)
                })
                .sum();
            volume < 0.0
        } else {
            flipped * 2 > members.len()
        };
        for &t in &members {
            if flip[t].unwrap_or(false) ^ turn_all {
                triangles[t].swap(1, 2);
                report.windings_flipped += 1;
            }
        }
        let mut members: Vec<u32> = members
            .into_iter()
            .map(|t| u32::try_from(t).unwrap_or(u32::MAX))
            .collect();
        members.sort_unstable();
        pieces.push(Piece {
            triangles: members,
            closed,
        });
    }
    pieces
}

fn signed_volume(points: &[Point], [a, b, c]: [u32; 3]) -> f64 {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize] - Point::ORIGIN);
    a.dot(b.cross(c)) / 6.0
}

/// Whether piece `inner` lies inside closed piece `outer`: a ray from one
/// of its vertices crosses `outer` an odd number of times.
///
/// A ray through a shared edge or vertex of `outer` meets two triangles
/// there, or none, and its count says nothing; so a ray passing that close
/// to any triangle's boundary is set aside and the next direction tried.
/// Where every direction grazes something, the parities they read vote.
fn inside(points: &[Point], triangles: &[[u32; 3]], outer: &Piece, inner: &Piece) -> bool {
    let Some(&first) = inner.triangles.first() else {
        return false;
    };
    let origin = points[triangles[first as usize][0] as usize];
    // Off-axis directions, so a ray along a mesh's grid lines is unlikely;
    // the crossing test needs no unit length.
    const DIRECTIONS: [[f64; 3]; 6] = [
        [0.577_215_664_9, 0.618_033_988_7, 0.533_751_168_7],
        [-0.412_310_562_6, 0.723_606_797_7, 0.553_574_358_9],
        [0.682_384_715_9, -0.291_637_412_8, 0.669_740_133_4],
        [0.267_949_192_4, 0.414_213_562_4, -0.869_565_217_4],
        [-0.713_825_491_7, -0.327_419_853_1, 0.618_574_239_6],
        [0.229_416_525_6, -0.881_784_197_0, -0.412_310_562_6],
    ];
    // How near a triangle's boundary, in its own barycentric terms, a hit
    // stands too close to say which triangle it belongs to.
    const GRAZE: f64 = 1e-9;
    let mut odd = 0_usize;
    let mut read = 0_usize;
    for direction in DIRECTIONS {
        let direction = Vector::new(direction[0], direction[1], direction[2]);
        let mut crossings = 0_usize;
        let mut grazed = false;
        for &t in &outer.triangles {
            let [a, b, c] = triangles[t as usize].map(|i| points[i as usize]);
            let (e1, e2) = (b - a, c - a);
            let p = direction.cross(e2);
            let det = e1.dot(p);
            if det.abs() < 1e-300 {
                continue;
            }
            let s = origin - a;
            let u = s.dot(p) / det;
            if !(-GRAZE..=1.0 + GRAZE).contains(&u) {
                continue;
            }
            let q = s.cross(e1);
            let v = direction.dot(q) / det;
            if v < -GRAZE || u + v > 1.0 + GRAZE {
                continue;
            }
            if e2.dot(q) / det <= 0.0 {
                continue;
            }
            if u < GRAZE || v < GRAZE || u + v > 1.0 - GRAZE {
                grazed = true;
                break;
            }
            crossings += 1;
        }
        read += 1;
        if crossings % 2 == 1 {
            odd += 1;
        }
        if !grazed {
            return crossings % 2 == 1;
        }
    }
    odd * 2 > read
}

/// What a face is built on.
#[derive(Debug, Clone)]
enum Carrier {
    /// A plane, its normal outward.
    Plane(Plane),
    /// A surface recognition decided the region is.
    Curved(Curved),
    /// A region demoted to faceted, whose triangles went to other faces.
    Gone,
}

#[derive(Debug, Clone)]
struct Curved {
    shape: Canonical,
    /// How far the region's vertices stand off it.
    deviation: f64,
    /// How far they stood off the region's own fit, before a shared axis
    /// or a sphere's frame moved it: how hard the fit pressed against the
    /// distance it was grown within.
    fitted: f64,
    /// The chart point every pcurve on it is unwrapped around, so they all
    /// read on one branch.
    centre: (f64, f64),
    /// Whether the region runs all the way round its axis, and so needs a
    /// seam across it.
    wraps: bool,
    /// Whether it runs all the way round a torus's tube, and so needs a
    /// seam along it.
    wraps_v: bool,
    /// Whether the frame was set from the region's own boundary: a sphere
    /// bounded by circles in parallel planes turns about their normal, so
    /// they are its latitudes.
    fixed: bool,
    /// The region's mesh vertices.
    vertices: Vec<u32>,
    /// For a fitted patch, whether its chart is the mean-value map.
    patch: Option<bool>,
}

/// How a curved face's boundary lies in its chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// A patch inside one branch of the chart.
    Open,
    /// A band between two full circles, joined by a seam: round the axis,
    /// or round a torus's tube.
    Band { round_tube: bool },
    /// A sphere's cap: one latitude circle, a seam to the pole, and the pole.
    Cap,
    /// The whole surface, with no boundary of its own.
    Whole,
    /// Round the axis between two rims of any shape, with holes: a seam
    /// joins a vertex of each rim, straight in the chart and clear of the
    /// holes.
    Wrapped,
    /// A sphere or a torus whole but for holes, its seams and poles clear
    /// of them: the whole surface's face, the holes its inner wires.
    Holed,
}

/// Triangles gathered into faces: which face each triangle is in, and what
/// each face is built on.
#[derive(Clone)]
struct Groups {
    of: Vec<usize>,
    carriers: Vec<Carrier>,
    /// The smooth regions refused a patch, by why.
    refused: PatchRefusals,
}

/// How many smooth regions were refused a patch, by why.
#[derive(Debug, Clone, Copy, Default)]
struct PatchRefusals {
    not_disk: usize,
    narrow: usize,
    unverified: usize,
}

/// The plane of a triangle, its normal by the winding, its `x` axis along
/// the first side. The triangle has area, so both are defined.
fn plane_of(points: &[Point], [a, b, c]: [u32; 3], tol: Tolerances) -> OgeomResult<Plane> {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    // Scaled to unit length before the kernel's own normalization, which
    // refuses a vector shorter than the confusion distance, as the cross
    // product of a small triangle's sides may be.
    let normal = (b - a).cross(c - a);
    let z = Direction::new(normal / normal.magnitude(), tol)?;
    let x = Direction::new((b - a) / (b - a).magnitude(), tol)?;
    Ok(Plane::new(Frame::new(a, z, x, tol)?))
}

fn unit_normal(points: &[Point], [a, b, c]: [u32; 3]) -> Vector {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    let n = (b - a).cross(c - a);
    n / n.magnitude()
}

/// Grow planar faces across shared edges from the largest triangles down,
/// among the triangles no face holds yet, taking a neighbour whose normal
/// is within the angle of the face's and whose corners all lie within the
/// distance of its plane. Measured against the face's own plane, not the
/// last triangle's, so a gently curved surface does not drift into one
/// face: the seed's plane at first, then the plane fitted to what the face
/// has gathered.
fn coplanar_groups(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    let area = |t: usize| {
        let [a, b, c] = triangles[t].map(|i| points[i as usize]);
        (b - a).cross(c - a).magnitude()
    };
    let mut order: Vec<usize> = (0..triangles.len())
        .filter(|&t| groups.of[t] == usize::MAX)
        .collect();
    order.sort_by(|&x, &y| area(y).total_cmp(&area(x)));
    let cos = options.coplanar_angle.cos();
    let quantum = options.quantum.unwrap_or(0.0);
    let mut stack = Vec::new();
    for seed in order {
        if groups.of[seed] != usize::MAX {
            continue;
        }
        let g = groups.carriers.len();
        let mut plane = plane_of(points, triangles[seed], tol)?;
        groups.of[seed] = g;
        let mut members = vec![seed];
        let mut vertices: Vec<u32> = triangles[seed].to_vec();
        let mut rim = vec![seed];
        loop {
            let (origin, normal) = (plane.frame().origin(), plane.frame().z().vector());
            stack.clone_from(&rim);
            rim.clear();
            let grown = members.len();
            while let Some(t) = stack.pop() {
                for h in 3 * t..3 * t + 3 {
                    let Some(twin) = adjacency.twin[h] else {
                        continue;
                    };
                    let other = twin / 3;
                    if groups.of[other] != usize::MAX {
                        continue;
                    }
                    let [a, b, c] = triangles[other].map(|i| points[i as usize]);
                    let n = (b - a).cross(c - a);
                    // A small triangle's corners, each rounded by the
                    // mesh's encoding, tilt its normal by up to twice the
                    // rounding over its smallest altitude: on a mesh stored
                    // far from its origin, a tilt that is no turn of the
                    // surface. Its corners' distance decides it then.
                    let longest = (b - a)
                        .magnitude()
                        .max((c - b).magnitude())
                        .max((a - c).magnitude());
                    let altitude = n.magnitude() / longest.max(f64::MIN_POSITIVE);
                    let slop = (2.0 * quantum / altitude.max(f64::MIN_POSITIVE))
                        .min(1.0)
                        .asin();
                    // A triangle no higher than the distance (three corners
                    // all but on a line) has no plane of its own at that
                    // distance: built alone it is a face of almost no area,
                    // turned whichever way its rounding points.
                    let cos = if altitude <= flat {
                        -1.0
                    } else {
                        cos.min(slop.cos())
                    };
                    if n.dot(normal) < cos * n.magnitude()
                        || [a, b, c]
                            .iter()
                            .any(|p| (*p - origin).dot(normal).abs() > flat)
                    {
                        rim.push(t);
                        continue;
                    }
                    groups.of[other] = g;
                    members.push(other);
                    vertices.extend(triangles[other]);
                    stack.push(other);
                }
            }
            // A plane through one triangle leans with its corners' slop, and
            // across a large face that lean carries the far side past the
            // distance. Refitted to everything gathered, where the fit still
            // holds it all, the face grows again from its rim.
            if members.len() == grown || rim.is_empty() {
                break;
            }
            vertices.sort_unstable();
            vertices.dedup();
            let at: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
            let Some((centre, fitted)) = plane_through(&at, tol) else {
                break;
            };
            let fitted = if fitted.vector().dot(normal) < 0.0 {
                fitted.reversed()
            } else {
                fitted
            };
            if at
                .iter()
                .any(|p| (*p - centre).dot(fitted.vector()).abs() > flat)
            {
                break;
            }
            let x = plane.frame().x().vector();
            let x = x - fitted.vector() * x.dot(fitted.vector());
            plane = Plane::new(Frame::new(centre, fitted, Direction::new(x, tol)?, tol)?);
            rim.sort_unstable();
            rim.dedup();
        }
        rim.clear();
        groups.carriers.push(Carrier::Plane(plane));
    }
    Ok(())
}

/// One face per triangle, for the triangles no face holds yet.
fn one_each(
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    for (t, triangle) in triangles.iter().enumerate() {
        if groups.of[t] == usize::MAX {
            groups.of[t] = groups.carriers.len();
            groups
                .carriers
                .push(Carrier::Plane(plane_of(points, *triangle, tol)?));
        }
    }
    Ok(())
}

/// The direction of a canonical surface's own normal at a point near it,
/// not normalized: the gradient of its distance.
/// How far a sample's normals must turn before the smallest sample stage
/// is fitted: twenty degrees.
const SMALL_STAGE_TURN: f64 = 0.35;

/// The widest angle between any two of a sample's normals.
fn turn_of(normals: &[Vector]) -> f64 {
    let mut widest: f64 = 0.0;
    for (i, a) in normals.iter().enumerate() {
        for b in &normals[i + 1..] {
            widest = widest.max(a.dot(*b).clamp(-1.0, 1.0).acos());
        }
    }
    widest
}

/// Whether a facet whose corners lie on `shape` leans only as the surface
/// turns under it: its normal within the spread of the surface's normals
/// at its corners (and [`FACET_LEAN`]) of their mean, the surface turning
/// no more than sixty degrees under it. A facet of a coarse mesh over a
/// tight bend passes however far it turns from its neighbours; a flat
/// cap's triangle with its corners on a cylinder's rim, square to every
/// one of those normals, does not.
fn leans_as_the_surface(shape: &Canonical, corners: [Point; 3], normal: Vector) -> bool {
    let mut at = [Vector::ZERO; 3];
    for (n, p) in at.iter_mut().zip(corners) {
        let g = gradient(shape, p);
        let m = g.magnitude();
        if m == 0.0 {
            return false;
        }
        *n = if g.dot(normal) < 0.0 { -g / m } else { g / m };
    }
    let angle = |a: Vector, b: Vector| a.dot(b).clamp(-1.0, 1.0).acos();
    let spread = angle(at[0], at[1])
        .max(angle(at[1], at[2]))
        .max(angle(at[0], at[2]));
    // A facet of any mesh turns far less than this across itself; one that
    // spans more is a chord across the surface, not a piece of it.
    if spread > FACET_TURN {
        return false;
    }
    let mean = at[0] + at[1] + at[2];
    let m = mean.magnitude();
    m > 0.0 && angle(mean / m, normal) <= spread + FACET_LEAN
}

/// How far a facet's normal may lean past the spread of the surface's
/// normals at its corners: a skewed triangle across a cylinder, joining
/// points at different heights and angles, tilts along the axis out of the
/// plane its corners' normals span. About eleven degrees.
const FACET_LEAN: f64 = 0.2;

/// The most a surface may turn under one facet: sixty degrees.
const FACET_TURN: f64 = core::f64::consts::FRAC_PI_3;

fn gradient(shape: &Canonical, p: Point) -> Vector {
    let radial = |o: Point, z: Vector| {
        let w = p - o;
        let r = w - z * w.dot(z);
        let m = r.magnitude();
        (if m > 0.0 { r / m } else { Vector::ZERO }, w.dot(z))
    };
    match shape {
        Canonical::Plane(plane) => plane.frame().z().vector(),
        Canonical::Cylinder(c) => radial(c.frame().origin(), c.frame().z().vector()).0,
        Canonical::Cone(c) => {
            let z = c.frame().z().vector();
            let (out, _) = radial(c.frame().origin(), z);
            out - z * c.half_angle().tan()
        }
        Canonical::Sphere(s) => p - s.centre(),
        Canonical::Torus(t) => {
            let z = t.frame().z().vector();
            let (out, _) = radial(t.frame().origin(), z);
            p - (t.frame().origin() + out * t.major_radius())
        }
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            s.foot(p)
                .and_then(|f| {
                    s.surface
                        .normal_at(f.parameters.0, f.parameters.1, s.tol)
                        .ok()
                })
                .map_or(Vector::ZERO, |n| n.vector())
        }
    }
}

/// The axis of a surface of revolution: its frame.
fn axis_frame(shape: &Canonical) -> Option<Frame> {
    match shape {
        Canonical::Cylinder(c) => Some(c.frame()),
        Canonical::Cone(c) => Some(c.frame()),
        Canonical::Torus(t) => Some(t.frame()),
        Canonical::Sphere(s) => Some(s.frame()),
        Canonical::Plane(_) | Canonical::Swept(_) => None,
    }
}

/// The same surface on another frame whose axis is the same line: its
/// radii kept, a cone's reference radius carried to the new origin.
fn on_frame(shape: &Canonical, frame: Frame, tol: Tolerances) -> Option<Canonical> {
    Some(match shape {
        Canonical::Cylinder(c) => Canonical::Cylinder(Cylinder::new(frame, c.radius(), tol).ok()?),
        Canonical::Cone(c) => {
            // The new origin's height on the old axis, and which way the
            // new axis runs against the old.
            let old = c.frame();
            let shift = (frame.origin() - old.origin()).dot(old.z().vector());
            let same = frame.z().vector().dot(old.z().vector()) > 0.0;
            if !same {
                return None;
            }
            let r0 = c.radius_at(shift);
            Canonical::Cone(Cone::new(frame, r0.max(tol.confusion()), c.half_angle(), tol).ok()?)
        }
        Canonical::Torus(t) => {
            Canonical::Torus(Torus::new(frame, t.major_radius(), t.minor_radius(), tol).ok()?)
        }
        Canonical::Sphere(s) => Canonical::Sphere(Sphere::new(frame, s.radius(), tol).ok()?),
        Canonical::Plane(_) | Canonical::Swept(_) => return None,
    })
}

/// A point's raw chart coordinates on a canonical surface.
fn chart(shape: &Canonical, p: Point, tol: Tolerances) -> Option<(f64, f64)> {
    use ogeom_math::elementary as e;
    match shape {
        Canonical::Plane(plane) => Some(e::plane_parameters(plane, p)),
        Canonical::Cylinder(c) => e::cylinder_parameters(c, p, tol).ok(),
        Canonical::Cone(c) => e::cone_parameters(c, p, tol).ok(),
        Canonical::Sphere(s) => e::sphere_parameters(s, p, tol).ok(),
        Canonical::Torus(t) => e::torus_parameters(t, p, tol).ok(),
        Canonical::Swept(s) => s.foot(p).map(|f| f.parameters),
    }
}

fn evaluate(shape: &Canonical, (u, v): (f64, f64)) -> Point {
    use ogeom_math::elementary as e;
    match shape {
        Canonical::Plane(plane) => e::plane_at(plane, u, v).point,
        Canonical::Cylinder(c) => e::cylinder_at(c, u, v).point,
        Canonical::Cone(c) => e::cone_at(c, u, v).point,
        Canonical::Sphere(s) => e::sphere_at(s, u, v).point,
        Canonical::Torus(t) => e::torus_at(t, u, v).point,
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            s.surface.point_at(u, v, s.tol).unwrap_or(Point::ORIGIN)
        }
    }
}

/// Which chart directions are angles, and so wrap.
fn periodic(shape: &Canonical) -> (bool, bool) {
    match shape {
        Canonical::Plane(_) => (false, false),
        Canonical::Cylinder(_) | Canonical::Cone(_) | Canonical::Sphere(_) => (true, false),
        Canonical::Torus(_) => (true, true),
        Canonical::Swept(s) => {
            use ogeom_geom::Surface as _;
            (s.surface.is_periodic_u(), s.surface.is_periodic_v())
        }
    }
}

/// Chart coordinates on the branch around `centre`.
fn unwrapped(curved: &Curved, p: Point, tol: Tolerances) -> Option<(f64, f64)> {
    let (u, v) = chart(&curved.shape, p, tol)?;
    let (pu, pv) = periodic(&curved.shape);
    let near = |x: f64, c: f64, wraps: bool| {
        if wraps {
            c + ogeom_math::elementary::wrap_signed_angle(x - c)
        } else {
            x
        }
    };
    Some((near(u, curved.centre.0, pu), near(v, curved.centre.1, pv)))
}

/// The circular mean of angles, and the widest gap between them.
fn angular_spread(angles: &mut [f64]) -> (f64, f64) {
    let (s, c) = angles
        .iter()
        .fold((0.0, 0.0), |(s, c), a| (s + a.sin(), c + a.cos()));
    angles.sort_by(f64::total_cmp);
    let mut gap: f64 = 0.0;
    for w in angles.windows(2) {
        gap = gap.max(w[1] - w[0]);
    }
    if let (Some(first), Some(last)) = (angles.first(), angles.last()) {
        gap = gap.max(first + core::f64::consts::TAU - last);
    }
    (s.atan2(c), gap)
}

/// Grow regions of triangles that recognition says lie on one curved
/// canonical surface, then planar faces over the rest.
///
/// A region starts at a triangle with a curved edge (one across which the
/// surface turns by less than the crease angle but more than the coplanar
/// angle) and first grows across such edges only, which keeps a flat face
/// tangent to a fillet out of the fillet's first samples. Once it holds
/// enough vertices it is recognized, and from then on grows across any
/// smooth edge to a triangle whose corners lie on the surface and whose
/// normal agrees with it, the surface refitted to everything held as the
/// region doubles. A region whose samples are free-form, or also flat,
/// returns its triangles to the planar pass.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
fn segment(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    tol: Tolerances,
) -> OgeomResult<Groups> {
    let n = triangles.len();
    let mut groups = Groups {
        of: vec![usize::MAX; n],
        carriers: Vec::new(),
        refused: PatchRefusals::default(),
    };
    if !options.merge_coplanar {
        one_each(points, triangles, &mut groups, tol)?;
        return Ok(groups);
    }
    if options.recognize {
        recognized_regions(
            points,
            triangles,
            adjacency,
            options,
            flat,
            &mut groups,
            tol,
        );
        if options.sweeps {
            swept_regions(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut groups,
                tol,
            );
        }
        if options.patches {
            patch_regions(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut groups,
                tol,
            )?;
        }
        split_disconnected(triangles, adjacency, &mut groups);
        merge_same_surface(points, triangles, adjacency, &mut groups, flat);
        sphere_axes(points, triangles, adjacency, &mut groups, flat, tol);
        let mut planes = groups.clone();
        coplanar_groups(
            points,
            triangles,
            adjacency,
            options,
            flat,
            &mut planes,
            tol,
        )?;
        tangent_rounds(
            points,
            triangles,
            adjacency,
            &mut groups,
            &planes,
            options.crease.cos(),
            flat,
            tol,
        );
        let normals = plane_normals(&planes);
        align_axes(points, &mut groups, &normals, flat, tol);
        tangent_blends(
            points,
            triangles,
            adjacency,
            &mut groups,
            &planes,
            options.crease.cos(),
            flat,
            tol,
        );
        hole_frames(points, triangles, adjacency, &mut groups, tol);
        slit_bands(points, triangles, adjacency, &mut groups, tol);
    }
    coplanar_groups(
        points,
        triangles,
        adjacency,
        options,
        flat,
        &mut groups,
        tol,
    )?;
    Ok(groups)
}

/// The normals of the planes the triangles left over from recognition
/// gather into (`planes`, the groups with those planes grown), the largest
/// plane first: a small plane's fit leans with its few vertices' slop.
fn plane_normals(planes: &Groups) -> Vec<Direction> {
    let mut size = vec![0_usize; planes.carriers.len()];
    for &g in &planes.of {
        if let Some(n) = size.get_mut(g) {
            *n += 1;
        }
    }
    let mut normals: Vec<(usize, Direction)> = planes
        .carriers
        .iter()
        .zip(&size)
        .filter_map(|(c, &n)| match c {
            Carrier::Plane(plane) => Some((n, plane.frame().z())),
            _ => None,
        })
        .collect();
    normals.sort_by_key(|&(n, _)| core::cmp::Reverse(n));
    normals.into_iter().map(|(_, d)| d).collect()
}

/// Put a round between two flat faces on the cylinder tangent to both.
///
/// A round only a row or two of facets across has its vertices on a few
/// lines, which a cone, or a cylinder leaning from the faces, fits as well
/// as the round's own cylinder; the one fitted then meets the faces along
/// lines that are not where the facets end. The two faces fix the round's
/// axis (along their line of meeting) and leave only the radius, which each
/// vertex gives: the circle through it tangent to both faces. A curved
/// region meeting two non-parallel planes (of `planes`, the groups with the
/// planes grown) across smooth edges, and no other plane so, takes the
/// cylinder at the vertices' median radius where it holds every vertex
/// within the distance.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
fn tangent_rounds(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    planes: &Groups,
    cos_crease: f64,
    flat: f64,
    tol: Tolerances,
) {
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, &g) in groups.of.iter().enumerate() {
        if let Some(list) = members.get_mut(g) {
            list.push(t);
        }
    }
    for (i, region) in members.iter().enumerate() {
        let Carrier::Curved(curved) = &groups.carriers[i] else {
            continue;
        };
        if matches!(curved.shape, Canonical::Sphere(_) | Canonical::Torus(_)) {
            continue;
        }
        let mut beside: Vec<usize> = Vec::new();
        for &t in region {
            for h in 3 * t..3 * t + 3 {
                let Some(g) = adjacency.twin[h] else {
                    continue;
                };
                let other = g / 3;
                if groups.of[other] == i
                    || unit_normal(points, triangles[t]).dot(unit_normal(points, triangles[other]))
                        < cos_crease
                {
                    continue;
                }
                beside.push(planes.of[other]);
            }
        }
        beside.sort_unstable();
        beside.dedup();
        let flanks: Vec<&Plane> = beside
            .iter()
            .filter_map(|&j| match planes.carriers.get(j) {
                Some(Carrier::Plane(plane)) => Some(plane),
                _ => None,
            })
            .collect();
        if flanks.len() != 2 {
            continue;
        }
        let pts: Vec<Point> = curved
            .vertices
            .iter()
            .map(|&v| points[v as usize])
            .collect();
        let Some(shape) = tangent_cylinder(flanks[0], flanks[1], &pts, None, tol) else {
            continue;
        };
        let deviation = worst_deviation(&shape, &pts);
        if deviation > flat {
            continue;
        }
        if let Carrier::Curved(curved) = &mut groups.carriers[i] {
            curved.shape = shape;
            curved.deviation = deviation;
        }
    }
}

/// The cylinder tangent to two planes, on the side of them the points are
/// on, at `radius` or else through the points at their median radius.
///
/// With the planes' unit normals `a` and `b` and `w = (a + b) / (1 + a·b)`,
/// the axis of a circle of radius `r` tangent to both runs through
/// `c0 + s·r·w`, `c0` on both planes and `s` the side (-1 within both, +1
/// beyond both). A point `q` from `c0` (square to the axis) is on that
/// circle where `r²(|w|² - 1) - 2s(q·w)r + |q|² = 0`, the larger root.
fn tangent_cylinder(
    a: &Plane,
    b: &Plane,
    pts: &[Point],
    radius: Option<f64>,
    tol: Tolerances,
) -> Option<Canonical> {
    let (na, nb) = (a.frame().z().vector(), b.frame().z().vector());
    let g = na.dot(nb);
    // Nearly parallel planes leave the axis to their slop; nearly opposite
    // ones a radius the points barely fix.
    if g.abs() > 5.0_f64.to_radians().cos() {
        return None;
    }
    let axis = Direction::new(na.cross(nb), tol).ok()?;
    let (pa, pb) = (
        a.frame().origin().to_vector(),
        b.frame().origin().to_vector(),
    );
    // The point on both planes nearest the origin of the axis's normal
    // plane: solve c·na = pa·na, c·nb = pb·nb, c·axis = 0.
    let (ha, hb) = (pa.dot(na), pb.dot(nb));
    let c0 = (nb.cross(axis.vector()) * ha + axis.vector().cross(na) * hb)
        / na.cross(nb).dot(axis.vector());
    let w = (na + nb) / (1.0 + g);
    let side = pts
        .iter()
        .map(|p| (p.to_vector() - c0).dot(na) + (p.to_vector() - c0).dot(nb))
        .sum::<f64>()
        .signum();
    let k = w.dot(w) - 1.0;
    let mut radii: Vec<f64> = pts
        .iter()
        .filter_map(|p| {
            let d = p.to_vector() - c0;
            let q = d - axis.vector() * d.dot(axis.vector());
            let qw = side * q.dot(w);
            let disc = qw * qw - k * q.dot(q);
            (qw > 0.0).then(|| (qw + disc.max(0.0).sqrt()) / k)
        })
        .collect();
    if radii.len() < pts.len() / 2 + 1 {
        return None;
    }
    let radius = match radius {
        Some(radius) => radius,
        None => {
            let at = radii.len() / 2;
            *radii.select_nth_unstable_by(at, f64::total_cmp).1
        }
    };
    let through = Point::from_vector(c0 + w * (side * radius));
    Some(Canonical::Cylinder(
        Cylinder::new(Frame::about(through, axis), radius, tol).ok()?,
    ))
}

/// Put fillets and corner balls on the surfaces their neighbours fix.
///
/// A fillet is fitted as freely as any region, so it meets the faces it
/// blends into at a slight angle or a slight gap, and the seam solved
/// between two nearly tangent surfaces wanders along them. A fillet's
/// supports fix it but for its radius:
///
/// - rounds between two planes that meet at corner balls are one rolling
///   ball's: they take their median radius together where their own
///   radii agree to within `flat`, and each ball is centred a radius off
///   the planes its rounds run between, the point all their axes pass
///   through;
/// - a torus between a plane and a cylinder or a cone whose axis is square
///   to the plane sits on that axis, its tube's centre a radius off both;
/// - a sphere where cylinders of its own radius meet otherwise is centred
///   nearest their axes.
///
/// The surface derived so is tangent to its supports by construction, and
/// replaces the fitted one where it holds every vertex of the region within
/// `flat`, on the fitted one's frame so the chart branch fixed for it still
/// holds. Other regions keep their fits. `planes` are the groups with the
/// planes grown, as for [`tangent_rounds`].
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
fn tangent_blends(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    planes: &Groups,
    cos_crease: f64,
    flat: f64,
    tol: Tolerances,
) {
    let count = groups.carriers.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (t, &g) in groups.of.iter().enumerate() {
        if let Some(list) = members.get_mut(g) {
            list.push(t);
        }
    }
    // Each curved region's planes and curved regions across its smooth
    // edges, and its vertices.
    let mut flanks: Vec<Vec<Plane>> = vec![Vec::new(); count];
    let mut supports: Vec<Vec<usize>> = vec![Vec::new(); count];
    let mut samples: Vec<Vec<Point>> = vec![Vec::new(); count];
    for (i, region) in members.iter().enumerate() {
        let Carrier::Curved(curved) = &groups.carriers[i] else {
            continue;
        };
        let mut beside: Vec<usize> = Vec::new();
        for &t in region {
            for h in 3 * t..3 * t + 3 {
                let Some(g) = adjacency.twin[h] else {
                    continue;
                };
                let other = g / 3;
                if groups.of[other] == i
                    || unit_normal(points, triangles[t]).dot(unit_normal(points, triangles[other]))
                        < cos_crease
                {
                    continue;
                }
                beside.push(planes.of[other]);
            }
        }
        beside.sort_unstable();
        beside.dedup();
        for &j in &beside {
            match (planes.carriers.get(j), groups.carriers.get(j)) {
                (Some(Carrier::Plane(plane)), _) => flanks[i].push(*plane),
                (_, Some(Carrier::Curved(_))) if j != i => supports[i].push(j),
                _ => {}
            }
        }
        samples[i] = curved
            .vertices
            .iter()
            .map(|&v| points[v as usize])
            .collect();
    }
    let shape_of = |groups: &Groups, i: usize| match &groups.carriers[i] {
        Carrier::Curved(c) => Some(c.shape.clone()),
        _ => None,
    };
    let is_round = |groups: &Groups, i: usize| {
        flanks[i].len() == 2 && matches!(shape_of(groups, i), Some(Canonical::Cylinder(_)))
    };
    let is_ball =
        |groups: &Groups, i: usize| matches!(shape_of(groups, i), Some(Canonical::Sphere(_)));
    // Rounds joined through the balls they meet at.
    let mut chain: Vec<usize> = (0..count).collect();
    fn root(chain: &mut [usize], mut i: usize) -> usize {
        while chain[i] != i {
            chain[i] = chain[chain[i]];
            i = chain[i];
        }
        i
    }
    for (i, beside) in supports.iter().enumerate() {
        if !is_ball(groups, i) {
            continue;
        }
        for &j in beside {
            if is_round(groups, j) {
                let (a, b) = (root(&mut chain, i), root(&mut chain, j));
                chain[a] = b;
            }
        }
    }
    let mut parts: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..count {
        if is_round(groups, i) || is_ball(groups, i) {
            parts.entry(root(&mut chain, i)).or_default().push(i);
        }
    }
    let mut settled = vec![false; count];
    for part in parts.values() {
        let rounds: Vec<usize> = part
            .iter()
            .copied()
            .filter(|&j| is_round(groups, j))
            .collect();
        let balls: Vec<usize> = part
            .iter()
            .copied()
            .filter(|&j| is_ball(groups, j))
            .collect();
        if balls.is_empty() || rounds.is_empty() {
            continue;
        }
        let mut radii: Vec<f64> = rounds
            .iter()
            .filter_map(|&j| match shape_of(groups, j) {
                Some(Canonical::Cylinder(c)) => Some(c.radius()),
                _ => None,
            })
            .collect();
        let (lo, hi) = radii
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &r| {
                (lo.min(r), hi.max(r))
            });
        if hi - lo > flat {
            continue;
        }
        let at = radii.len() / 2;
        let radius = *radii.select_nth_unstable_by(at, f64::total_cmp).1;
        for &j in &rounds {
            let Some(old) = shape_of(groups, j) else {
                continue;
            };
            let derived =
                tangent_cylinder(&flanks[j][0], &flanks[j][1], &samples[j], Some(radius), tol)
                    .and_then(|d| on_frame_of(&d, &old, tol));
            if let Some(shape) = derived {
                settled[j] = put(groups, j, shape, &samples[j], flat);
            }
        }
        for &j in &balls {
            let Some(Canonical::Sphere(old)) = shape_of(groups, j) else {
                continue;
            };
            let mut around: Vec<Plane> = Vec::new();
            for &k in &supports[j] {
                if rounds.contains(&k) {
                    around.extend(flanks[k].iter().copied());
                }
            }
            if let Some(shape) = ball_off_planes(&old, &around, &samples[j], radius, tol) {
                settled[j] = put(groups, j, shape, &samples[j], flat);
            }
        }
    }
    for i in 0..count {
        if settled[i] {
            continue;
        }
        let derived = match shape_of(groups, i) {
            Some(Canonical::Torus(torus)) => {
                // The support is the band on the torus's own axis; other
                // fillets running into it smoothly are its neighbours, not
                // its supports.
                let (o, z) = (torus.frame().origin(), torus.frame().z().vector());
                let coaxial: Vec<Canonical> = supports[i]
                    .iter()
                    .filter_map(|&j| shape_of(groups, j))
                    .filter(|s| {
                        axis_frame(s).is_some_and(|f| {
                            let w = f.origin() - o;
                            f.z().vector().cross(z).magnitude() <= 1e-3
                                && (w - z * w.dot(z)).magnitude() <= flat * 10.0
                        })
                    })
                    .collect();
                match (&flanks[i][..], &coaxial[..]) {
                    ([plane], [support]) => tangent_torus(plane, support, &samples[i], tol),
                    _ => None,
                }
            }
            Some(Canonical::Sphere(sphere)) => {
                let around: Vec<Canonical> = supports[i]
                    .iter()
                    .filter_map(|&j| shape_of(groups, j))
                    .collect();
                corner_ball(&sphere, &around, flat, tol)
            }
            _ => None,
        };
        if let Some(shape) = derived {
            put(groups, i, shape, &samples[i], flat);
        }
    }
}

/// `shape` as group `i`'s surface where it holds every sample within
/// `flat`; whether it does.
fn put(groups: &mut Groups, i: usize, shape: Canonical, samples: &[Point], flat: f64) -> bool {
    let deviation = worst_deviation(&shape, samples);
    if deviation > flat {
        return false;
    }
    if let Carrier::Curved(curved) = &mut groups.carriers[i] {
        curved.shape = shape;
        curved.deviation = deviation;
    }
    true
}

/// A cylinder put on the frame of the one it replaces: its axis turned to
/// run the same way, its angle measured from the same direction and its
/// height from the same level, so the chart branch fixed for the old one
/// reads the same on the new.
fn on_frame_of(shape: &Canonical, old: &Canonical, tol: Tolerances) -> Option<Canonical> {
    let (Canonical::Cylinder(new), Some(was)) = (shape, axis_frame(old)) else {
        return None;
    };
    let z = new.frame().z();
    let z = if z.vector().dot(was.z().vector()) >= 0.0 {
        z
    } else {
        -z
    };
    let x = was.x().vector() - z.vector() * was.x().vector().dot(z.vector());
    let origin =
        new.frame().origin() + z.vector() * (was.origin() - new.frame().origin()).dot(z.vector());
    let frame = Frame::new(origin, z, Direction::new(x, tol).ok()?, tol).ok()?;
    Some(Canonical::Cylinder(
        Cylinder::new(frame, new.radius(), tol).ok()?,
    ))
}

/// The ball a radius off the planes its rounds run between, on the side of
/// each the samples are on, on the sphere's own frame. `None` unless the
/// planes are three that meet in a point.
fn ball_off_planes(
    sphere: &Sphere,
    around: &[Plane],
    samples: &[Point],
    radius: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let mut distinct: Vec<Plane> = Vec::new();
    for plane in around {
        let n = plane.frame().z().vector();
        let same = distinct.iter().any(|q| {
            let m = q.frame().z().vector();
            n.cross(m).magnitude() <= 1e-9
                && (plane.frame().origin() - q.frame().origin()).dot(m).abs() <= tol.confusion()
        });
        if !same {
            distinct.push(*plane);
        }
    }
    let [a, b, c] = distinct[..] else {
        return None;
    };
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let count = samples.len().max(1) as f64;
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut rhs = nalgebra::Vector3::<f64>::zeros();
    for (row, plane) in [a, b, c].iter().enumerate() {
        let (o, n) = (plane.frame().origin(), plane.frame().z().vector());
        let side = (samples.iter().map(|p| (*p - o).dot(n)).sum::<f64>() / count).signum();
        m.set_row(row, &nalgebra::RowVector3::new(n.x, n.y, n.z));
        rhs[row] = o.to_vector().dot(n) + side * radius;
    }
    let solved = m.lu().solve(&rhs)?;
    let centre = Point::new(solved[0], solved[1], solved[2]);
    let frame = sphere.frame();
    let on = Frame::new(centre, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Sphere(Sphere::new(on, radius, tol).ok()?))
}

/// The torus tangent to a plane and to a cylinder or cone whose axis is
/// square to it, on the side of each the points are on, through the
/// points at their median tube radius.
///
/// In the support's axial half-plane, with `rho` the distance from the axis
/// and `z` the height along it, the support's line is `rho = a + b z` and
/// the plane is `z = h`. A tube circle of radius `r` tangent to both has its
/// centre at `A + B r`, with `A = (a + b h, h)` and
/// `B = (s2 √(1 + b²) + s1 b, s1)`, `s1` and `s2` the sides of the plane and
/// the support the points are on. A point `q` lies on that circle where
/// `(|B|² - 1) r² - 2 (q - A)·B r + |q - A|² = 0`, the larger root.
fn tangent_torus(
    plane: &Plane,
    support: &Canonical,
    pts: &[Point],
    tol: Tolerances,
) -> Option<Canonical> {
    let (frame, a, b) = match support {
        Canonical::Cylinder(c) => (c.frame(), c.radius(), 0.0),
        Canonical::Cone(c) => (
            c.frame(),
            c.radius_at(0.0),
            c.radius_at(1.0) - c.radius_at(0.0),
        ),
        _ => return None,
    };
    let (o, z) = (frame.origin(), frame.z().vector());
    let n = plane.frame().z().vector();
    if n.cross(z).magnitude() > 1e-9 {
        return None;
    }
    let h = (plane.frame().origin() - o).dot(z);
    let local: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let w = *p - o;
            let along = w.dot(z);
            ((w - z * along).magnitude(), along)
        })
        .collect();
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let count = local.len().max(1) as f64;
    let s1 = (local.iter().map(|q| q.1 - h).sum::<f64>() / count).signum();
    let s2 = (local.iter().map(|q| q.0 - a - b * q.1).sum::<f64>() / count).signum();
    let base = (a + b * h, h);
    let step = (s2 * b.hypot(1.0) + s1 * b, s1);
    let k = step.0.mul_add(step.0, step.1 * step.1) - 1.0;
    if k <= 0.0 {
        return None;
    }
    let mut radii: Vec<f64> = local
        .iter()
        .filter_map(|q| {
            let d = (q.0 - base.0, q.1 - base.1);
            let along = d.0.mul_add(step.0, d.1 * step.1);
            let disc = along * along - k * d.0.mul_add(d.0, d.1 * d.1);
            (along > 0.0).then(|| (along + disc.max(0.0).sqrt()) / k)
        })
        .collect();
    if radii.len() < local.len() / 2 + 1 {
        return None;
    }
    let at = radii.len() / 2;
    let (_, radius, _) = radii.select_nth_unstable_by(at, f64::total_cmp);
    let radius = *radius;
    let centre = (
        step.0.mul_add(radius, base.0),
        step.1.mul_add(radius, base.1),
    );
    if centre.0 <= radius {
        return None;
    }
    let on = Frame::new(o + z * centre.1, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Torus(
        Torus::new(on, centre.0, radius, tol).ok()?,
    ))
}

/// The ball where cylinders of its own radius meet: centred at the point
/// nearest their axes, its radius theirs, on the sphere's own frame. The
/// cylinders taken are those within a hundredth of the fitted sphere's
/// radius; `None` unless there are at least two, their radii agree to
/// within `flat`, their axes are not parallel, and every axis passes within
/// `flat` of the centre.
fn corner_ball(
    sphere: &Sphere,
    supports: &[Canonical],
    flat: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let rounds: Vec<&Cylinder> = supports
        .iter()
        .filter_map(|s| match s {
            Canonical::Cylinder(c)
                if (c.radius() - sphere.radius()).abs() <= sphere.radius() * 0.01 =>
            {
                Some(c)
            }
            _ => None,
        })
        .collect();
    if rounds.len() < 2 {
        return None;
    }
    let (lo, hi) = rounds
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), c| {
            (lo.min(c.radius()), hi.max(c.radius()))
        });
    if hi - lo > flat {
        return None;
    }
    // Least squares: the sum over the axes of the projections square to
    // each, applied to the centre, equals the same applied to their origins.
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut rhs = nalgebra::Vector3::<f64>::zeros();
    for c in &rounds {
        let z = c.frame().z().vector();
        let zv = nalgebra::Vector3::new(z.x, z.y, z.z);
        let across = nalgebra::Matrix3::identity() - zv * zv.transpose();
        let o = c.frame().origin();
        m += across;
        rhs += across * nalgebra::Vector3::new(o.x, o.y, o.z);
    }
    let eigen = m.symmetric_eigen();
    if eigen.eigenvalues.min() < 1e-3 {
        return None;
    }
    let solved = m.lu().solve(&rhs)?;
    let centre = Point::new(solved[0], solved[1], solved[2]);
    let met = rounds.iter().all(|c| {
        let z = c.frame().z().vector();
        let w = centre - c.frame().origin();
        (w - z * w.dot(z)).magnitude() <= flat
    });
    if !met {
        return None;
    }
    #[allow(clippy::cast_precision_loss, reason = "a few cylinders")]
    let radius = rounds.iter().map(|c| c.radius()).sum::<f64>() / rounds.len() as f64;
    let frame = sphere.frame();
    let on = Frame::new(centre, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Sphere(Sphere::new(on, radius, tol).ok()?))
}

/// How much larger than a curved region's own facets a flat patch must be
/// to be a face of its own rather than facets of the curve.
const PATCH_SCALE: f64 = 20.0;

/// How much larger than a seed's flat patch another must be for a first
/// sample to take its facets without growing on from them.
const LEAF_SCALE: f64 = 16.0;

/// The mesh cut into flat patches: each grown from its largest facet across
/// edges to neighbours whose normals stay within a couple of degrees of
/// that facet's and whose corners stay on its plane. A curved region takes
/// a large patch whole or not at all.
struct FlatPatches {
    /// The patch each triangle is in.
    of: Vec<usize>,
    /// Each patch's triangles.
    members: Vec<Vec<usize>>,
    /// Each patch's area.
    area: Vec<f64>,
    /// Each triangle's area.
    facet: Vec<f64>,
}

impl FlatPatches {
    fn of(
        points: &[Point],
        triangles: &[[u32; 3]],
        adjacency: &Adjacency,
        normals: &[Vector],
        flat: f64,
    ) -> Self {
        // Two degrees: the folds a curve's facets make are larger except on
        // a very finely drawn one, whose facets' corners then leave the seed's
        // plane within a few rows.
        let cos_fold = 2.0_f64.to_radians().cos();
        let n = triangles.len();
        let facet: Vec<f64> = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| points[v as usize]);
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&x, &y| facet[y].total_cmp(&facet[x]).then(x.cmp(&y)));
        let mut of = vec![usize::MAX; n];
        let (mut members, mut area) = (Vec::new(), Vec::new());
        for seed in order {
            if of[seed] != usize::MAX {
                continue;
            }
            let id = members.len();
            let normal = normals[seed];
            let origin = points[triangles[seed][0] as usize];
            let mut patch = vec![seed];
            of[seed] = id;
            let mut i = 0;
            while i < patch.len() {
                let t = patch[i];
                i += 1;
                for h in 3 * t..3 * t + 3 {
                    let Some(g) = adjacency.twin[h] else {
                        continue;
                    };
                    let u = g / 3;
                    if of[u] != usize::MAX || normals[u].dot(normal) < cos_fold {
                        continue;
                    }
                    if triangles[u]
                        .iter()
                        .all(|&v| (points[v as usize] - origin).dot(normal).abs() <= flat)
                    {
                        of[u] = id;
                        patch.push(u);
                    }
                }
            }
            area.push(patch.iter().map(|&t| facet[t]).sum());
            members.push(patch);
        }
        Self {
            of,
            members,
            area,
            facet,
        }
    }
}

/// How far the mesh's vertices stand off the flat faces they lie on: across
/// every edge between two triangles as good as coplanar (a turn of under a
/// twentieth of a degree, finer than any curve is drawn), the height of the
/// one's far corner over the other's plane, at the ninetieth percentile.
/// Zero where the mesh has no such edges.
fn flat_noise(points: &[Point], triangles: &[[u32; 3]], adjacency: &Adjacency) -> f64 {
    let cos_level = 0.05_f64.to_radians().cos();
    let mut heights: Vec<f64> = Vec::new();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(g) = *twin else {
            continue;
        };
        if g < h {
            continue;
        }
        let (t, u) = (h / 3, g / 3);
        let a = unit_normal(points, triangles[t]);
        let b = unit_normal(points, triangles[u]);
        if a.dot(b) < cos_level {
            continue;
        }
        let far = triangles[u][(g % 3 + 2) % 3];
        let base = points[triangles[t][0] as usize];
        heights.push((points[far as usize] - base).dot(a).abs());
    }
    if heights.is_empty() {
        return 0.0;
    }
    let at = (heights.len() * 9 / 10).min(heights.len() - 1);
    let (_, value, _) = heights.select_nth_unstable_by(at, f64::total_cmp);
    *value
}

/// Keep only the largest edge-connected piece of a set of triangles.
fn keep_largest_piece(region: &mut Vec<usize>, adjacency: &Adjacency) {
    let members: std::collections::HashSet<usize> = region.iter().copied().collect();
    let mut piece: HashMap<usize, usize> = HashMap::with_capacity(region.len());
    let mut sizes: Vec<usize> = Vec::new();
    for &start in region.iter() {
        if piece.contains_key(&start) {
            continue;
        }
        let id = sizes.len();
        piece.insert(start, id);
        let mut stack = vec![start];
        let mut size = 0;
        while let Some(t) = stack.pop() {
            size += 1;
            for h in 3 * t..3 * t + 3 {
                if let Some(twin) = adjacency.twin[h] {
                    let other = twin / 3;
                    if members.contains(&other) && !piece.contains_key(&other) {
                        piece.insert(other, id);
                        stack.push(other);
                    }
                }
            }
        }
        sizes.push(size);
    }
    if sizes.len() < 2 {
        return;
    }
    let largest = (0..sizes.len()).max_by_key(|&i| sizes[i]).unwrap_or(0);
    region.retain(|t| piece.get(t) == Some(&largest));
}

/// Recognized regions of one kind sharing a mesh edge and lying on one
/// surface are one region: each within twice the distance of the other's
/// (a band peeled of a flat face's facets can leave its two ends to be
/// fitted apart), or the smaller within the distance of the larger's, whose
/// surface the merged region keeps.
fn merge_same_surface(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    flat: f64,
) {
    loop {
        let mut pair: Option<(usize, usize)> = None;
        'find: for (t, tri) in triangles.iter().enumerate() {
            let g = groups.of[t];
            let Some(Carrier::Curved(a)) = groups.carriers.get(g) else {
                continue;
            };
            for h in 3 * t..3 * t + 3 {
                let Some(twin) = adjacency.twin[h] else {
                    continue;
                };
                let o = groups.of[twin / 3];
                if o == g {
                    continue;
                }
                let Some(Carrier::Curved(b)) = groups.carriers.get(o) else {
                    continue;
                };
                if core::mem::discriminant(&a.shape) != core::mem::discriminant(&b.shape) {
                    continue;
                }
                // Each was fitted to its own vertices within the distance,
                // so the other's lie within twice it where the two are one.
                let on = |shape: &Canonical, vertices: &[u32], within: f64| {
                    vertices
                        .iter()
                        .all(|&v| shape.distance_to(points[v as usize]) <= within)
                };
                if on(&a.shape, &b.vertices, flat * 2.0) && on(&b.shape, &a.vertices, flat * 2.0) {
                    pair = Some((g.min(o), g.max(o)));
                    break 'find;
                }
                // A small piece fitted on its own leans with its few
                // vertices' slop, and the larger's fit may miss it by more;
                // where its vertices lie on the larger's surface, it is part
                // of the larger.
                let (large, small) = if a.vertices.len() >= b.vertices.len() {
                    ((g, a), (o, b))
                } else {
                    ((o, b), (g, a))
                };
                if on(&large.1.shape, &small.1.vertices, flat) {
                    pair = Some((large.0, small.0));
                    break 'find;
                }
            }
            let _ = tri;
        }
        let Some((keep, gone)) = pair else {
            return;
        };
        let Carrier::Curved(absorbed) =
            core::mem::replace(&mut groups.carriers[gone], Carrier::Gone)
        else {
            return;
        };
        for of in &mut groups.of {
            if *of == gone {
                *of = keep;
            }
        }
        if let Carrier::Curved(kept) = &mut groups.carriers[keep] {
            let pts: Vec<Point> = absorbed
                .vertices
                .iter()
                .map(|&v| points[v as usize])
                .collect();
            kept.vertices.extend(absorbed.vertices);
            kept.vertices.sort_unstable();
            kept.vertices.dedup();
            kept.deviation = kept
                .deviation
                .max(absorbed.deviation)
                .max(worst_deviation(&kept.shape, &pts));
            kept.fitted = kept.fitted.max(absorbed.fitted);
        }
    }
}

/// One recognized region per connected patch. Triangles dropped from a
/// region for touching a vertex its fit refused can take with them the only
/// triangles joining the rest, and a face is built from one region's
/// boundary: the patches past the first would be lost from the shell. Each
/// patch past the first becomes a region of its own on the same surface.
fn split_disconnected(triangles: &[[u32; 3]], adjacency: &Adjacency, groups: &mut Groups) {
    let count = groups.carriers.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (t, &g) in groups.of.iter().enumerate() {
        if g < count && matches!(groups.carriers[g], Carrier::Curved(_)) {
            members[g].push(t);
        }
    }
    for (g, tris) in members.into_iter().enumerate() {
        if tris.len() < 2 {
            continue;
        }
        let mut patch = vec![usize::MAX; triangles.len()];
        let mut patches = 0_usize;
        for &start in &tris {
            if patch[start] != usize::MAX {
                continue;
            }
            patch[start] = patches;
            let mut stack = vec![start];
            while let Some(t) = stack.pop() {
                for h in 3 * t..3 * t + 3 {
                    if let Some(twin) = adjacency.twin[h] {
                        let other = twin / 3;
                        if groups.of[other] == g && patch[other] == usize::MAX {
                            patch[other] = patches;
                            stack.push(other);
                        }
                    }
                }
            }
            patches += 1;
        }
        if patches < 2 {
            continue;
        }
        let Carrier::Curved(template) = &groups.carriers[g] else {
            continue;
        };
        let template = template.clone();
        let base = groups.carriers.len();
        for k in 1..patches {
            let mut vertices: Vec<u32> = tris
                .iter()
                .filter(|&&t| patch[t] == k)
                .flat_map(|&t| triangles[t])
                .collect();
            vertices.sort_unstable();
            vertices.dedup();
            groups.carriers.push(Carrier::Curved(Curved {
                vertices,
                ..template.clone()
            }));
        }
        let mut first: Vec<u32> = tris
            .iter()
            .filter(|&&t| patch[t] == 0)
            .flat_map(|&t| triangles[t])
            .collect();
        first.sort_unstable();
        first.dedup();
        if let Carrier::Curved(curved) = &mut groups.carriers[g] {
            curved.vertices = first;
        }
        for &t in &tris {
            if patch[t] > 0 {
                groups.of[t] = base + patch[t] - 1;
            }
        }
    }
}

/// What a seed's first samples came to: the triangles and vertices
/// gathered, those fitted, how many triangles the smallest sample held, the
/// fit if one held, with which of the fitted vertices it kept.
struct FirstFit {
    region: Vec<usize>,
    vertices: Vec<u32>,
    shared: Vec<u32>,
    first_sample: usize,
    /// Whether the first sample already turned through a tight bend: on a
    /// coarse mesh it can span several surfaces, and a failed fit says
    /// nothing about the triangles in it but the seed.
    wide: bool,
    found: Option<(crate::recognize::Recognized, Vec<bool>)>,
}

/// The mesh as recognition reads it.
struct Surfaces<'a> {
    points: &'a [Point],
    triangles: &'a [[u32; 3]],
    adjacency: &'a Adjacency,
    normals: Vec<Vector>,
    patches: &'a FlatPatches,
    cos_crease: f64,
    cos_flat: f64,
    flat: f64,
    tol: Tolerances,
}

impl Surfaces<'_> {
    fn turn(&self, h: Half) -> Option<f64> {
        self.adjacency.twin[h].map(|g| self.normals[h / 3].dot(self.normals[g / 3]))
    }

    fn curved(&self, h: Half) -> bool {
        self.turn(h)
            .is_some_and(|c| c >= self.cos_crease && c < self.cos_flat)
    }

    fn smooth(&self, h: Half) -> bool {
        self.turn(h).is_some_and(|c| c >= self.cos_crease)
    }

    fn bends(&self, t: usize) -> bool {
        (3 * t..3 * t + 3).any(|h| self.curved(h))
    }

    /// Sample points with normals averaged over the region's triangles.
    fn samples(&self, vertices: &[u32], region: &[usize]) -> (Vec<Point>, Vec<Vector>) {
        let mut sum: HashMap<u32, Vector> = HashMap::with_capacity(vertices.len());
        for &t in region {
            for &v in &self.triangles[t] {
                *sum.entry(v).or_insert(Vector::ZERO) += self.normals[t];
            }
        }
        let pts = vertices.iter().map(|&v| self.points[v as usize]).collect();
        let nrm = vertices
            .iter()
            .map(|v| {
                let s = sum.get(v).copied().unwrap_or(Vector::Z);
                let m = s.magnitude();
                if m > 0.0 { s / m } else { Vector::Z }
            })
            .collect();
        (pts, nrm)
    }

    /// The region's edges as segments, a few hundred at most, evenly
    /// through the region.
    fn chords(&self, region: &[usize]) -> Vec<(Point, Point)> {
        let stride = region.len().div_ceil(100).max(1);
        region
            .iter()
            .step_by(stride)
            .flat_map(|&t| {
                let [a, b, c] = self.triangles[t].map(|v| self.points[v as usize]);
                [(a, b), (b, c), (c, a)]
            })
            .collect()
    }

    /// A seed's first samples and their fit, read from where the faces and
    /// the retired seeds stand; it changes nothing.
    ///
    /// Samples grow across smooth edges into triangles where the surface
    /// is seen to bend, each with an edge across which it turns. Within a
    /// strip of a cylinder or a torus the two triangles of a cell meet
    /// flat, and each still bends across its other side; a flat face met
    /// tangentially (a fillet's run-out) lends only its triangles along
    /// the tangent line, whose far corners the trimmed fit drops. The fit
    /// is tried as the sample grows: a narrow fillet fits from a few dozen
    /// vertices and would take in its neighbours' by a hundred; a patch of a
    /// thick torus needs the hundred to show its tube.
    fn first_fit(&self, seed: usize, of: &[usize], tried: &[bool]) -> FirstFit {
        const STAGES: [usize; 4] = [12, 24, 60, 150];
        let mut held: std::collections::HashSet<usize> = std::collections::HashSet::from([seed]);
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut region = vec![seed];
        let mut vertices: Vec<u32> = Vec::new();
        let take =
            |t: usize, vertices: &mut Vec<u32>, seen: &mut std::collections::HashSet<u32>| {
                for &v in &self.triangles[t] {
                    if seen.insert(v) {
                        vertices.push(v);
                    }
                }
            };
        take(seed, &mut vertices, &mut seen);
        let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::from([seed]);
        let mut found = None;
        let mut shared = Vec::new();
        let mut first_sample = usize::MAX;
        let mut wide = false;
        for target in STAGES {
            while vertices.len() < target {
                let Some(next) = queue.pop_front() else {
                    break;
                };
                for h in 3 * next..3 * next + 3 {
                    let Some(g) = self.adjacency.twin[h] else {
                        continue;
                    };
                    let other = g / 3;
                    if !self.smooth(h) || held.contains(&other) {
                        continue;
                    }
                    if of[other] != usize::MAX
                        || tried[other]
                        || !(self.curved(h) || self.bends(other))
                    {
                        continue;
                    }
                    held.insert(other);
                    region.push(other);
                    take(other, &mut vertices, &mut seen);
                    // A facet of a flat patch much larger than the seed's
                    // own (a flat face beside a fillet, where the seed's
                    // patch is a row of it) lends its corners but leads
                    // nowhere: on a coarse mesh it borders other fillets
                    // than the seed's.
                    let (patch, own) = (self.patches.of[other], self.patches.of[seed]);
                    if self.patches.area[patch] <= self.patches.area[own] * LEAF_SCALE {
                        queue.push_back(other);
                    }
                }
            }
            // Fitted on the vertices two or more of the sampled triangles
            // share: a corner of a flat face across a tangent line is
            // touched by the one triangle that reached it, and a single
            // such corner far off the surface decides a least-squares axis.
            shared = {
                let mut count: HashMap<u32, u32> = HashMap::with_capacity(vertices.len());
                for &t in &region {
                    for &v in &self.triangles[t] {
                        *count.entry(v).or_insert(0) += 1;
                    }
                }
                vertices
                    .iter()
                    .copied()
                    .filter(|v| count[v] >= 2)
                    .collect::<Vec<u32>>()
            };
            first_sample = first_sample.min(region.len());
            if target == STAGES[0] {
                let normals: Vec<Vector> = region.iter().map(|&t| self.normals[t]).collect();
                wide = turn_of(&normals) >= SMALL_STAGE_TURN;
            }
            if shared.len() < 8 {
                if queue.is_empty() {
                    break;
                }
                continue;
            }
            let (pts, nrm) = self.samples(&shared, &region);
            // The first, smallest stage is for a coarse mesh, where a dozen
            // vertices already span a tight bend; on a fine one they lie
            // nearly flat and say little about the surface they are on.
            if target == STAGES[0] && !wide {
                continue;
            }
            let chords = self.chords(&region);
            match crate::recognize::recognize_trimmed(&pts, &nrm, &chords, self.flat, self.tol) {
                Ok(fit) => {
                    found = Some(fit);
                    break;
                }
                // A larger sample is worth fitting only where this one came
                // near for its size (within a hundredth of its own span),
                // as a small patch of a thick torus does, its tube not yet
                // seen; a free-form patch misses by more, and is left.
                Err(closest) if closest > span(&pts) * 1e-2 || queue.is_empty() => break,
                Err(_) => {}
            }
        }
        FirstFit {
            region,
            vertices,
            shared,
            first_sample,
            wide,
            found,
        }
    }
}

/// Grow the recognized regions, seed by seed in triangle order.
///
/// A seed's first fit is the costly part, and it only reads: so seeds are
/// fitted a batch at a time in parallel against where things stand, then
/// taken in order, and a seed whose gathering took a triangle an earlier
/// seed of the same batch has since claimed or retired is fitted again.
/// The answer is the one seed-by-seed order gives, at any thread count.
#[allow(clippy::too_many_lines, reason = "one growth, read in one place")]
fn recognized_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) {
    let n = triangles.len();
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    let patches = FlatPatches::of(points, triangles, adjacency, &normals, flat);
    let mesh = Surfaces {
        points,
        triangles,
        adjacency,
        normals,
        patches: &patches,
        cos_crease: options.crease.cos(),
        cos_flat: options.coplanar_angle.cos(),
        flat,
        tol,
    };
    // A facet of a coarse mesh leans from the surface at its centre by up
    // to half the turn between facets, which the crease bounds.
    let agree = options.crease.cos();
    let mut tried = vec![false; n];
    // The batch in which each triangle's standing last changed.
    let mut changed = vec![0_u32; n];
    let batch_size = ogeom_core::parallel::threads().max(1) * 4;
    let mut batch = 0_u32;
    let eligible =
        |t: usize, of: &[usize], tried: &[bool]| of[t] == usize::MAX && !tried[t] && mesh.bends(t);
    // Seeds are taken in a fixed stride through the triangles rather than
    // one after another: a mesh lists neighbours together, and a batch of
    // neighbouring seeds would mostly gather what the first of them took.
    // The stride is a prime not dividing the count, so every triangle comes
    // up once.
    let stride = [7919_usize, 7907, 7901]
        .into_iter()
        .find(|p| !n.is_multiple_of(*p))
        .unwrap_or(1);
    let mut order: Vec<usize> = (0..n).map(|i| (i * stride) % n).collect();
    // Seeds whose wide first sample fitted nothing, with that sample. Once
    // every seed has been taken, those whose sample has since lost
    // triangles to other regions are taken once more, the rest of their
    // samples free again: a fillet sampled together with its unclaimed
    // neighbours fits nothing, and its triangles, retired one by one, leave
    // ever smaller samples; on its own, its neighbours claimed, it fits.
    let mut straddled: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut again = true;
    let mut next = 0;
    loop {
        if next == order.len() {
            let samples: Vec<(usize, Vec<usize>)> = std::mem::take(&mut straddled)
                .into_iter()
                .filter(|(seed, sample)| {
                    groups.of[*seed] == usize::MAX
                        && sample.iter().any(|&t| groups.of[t] != usize::MAX)
                })
                .collect();
            if !std::mem::take(&mut again) || samples.is_empty() {
                break;
            }
            for (seed, sample) in samples {
                for &t in &sample {
                    if groups.of[t] == usize::MAX {
                        tried[t] = false;
                    }
                }
                order.push(seed);
            }
        }
        batch += 1;
        let mut seeds = Vec::with_capacity(batch_size);
        while next < order.len() && seeds.len() < batch_size {
            if eligible(order[next], &groups.of, &tried) {
                seeds.push(order[next]);
            }
            next += 1;
        }
        let fits = {
            let (of, tried) = (&groups.of, &tried);
            ogeom_core::parallel::map_ordered(&seeds, |_, &seed| mesh.first_fit(seed, of, tried))
        };
        for (seed, fit) in seeds.into_iter().zip(fits) {
            if !eligible(seed, &groups.of, &tried) {
                continue;
            }
            // A triangle's standing only moves one way (free to taken), so
            // a triangle the gathering passed over stays passed over, and
            // only one it took can have changed what it gathers.
            let fit = if fit.region.iter().any(|&t| changed[t] == batch) {
                mesh.first_fit(seed, &groups.of, &tried)
            } else {
                fit
            };
            let FirstFit {
                mut region,
                mut vertices,
                shared,
                first_sample,
                wide,
                found,
                ..
            } = fit;
            let Some((found, keep)) = found else {
                // A fine sample that fits nothing lies on nothing canonical,
                // and its triangles are not seeded again; a wide one may
                // have straddled a fillet and its neighbours, and only the
                // seed is retired.
                let retired = if wide {
                    straddled.push((seed, region.clone()));
                    1
                } else {
                    first_sample
                };
                for &t in &region[..retired.min(region.len())] {
                    tried[t] = true;
                    changed[t] = batch;
                }
                continue;
            };
            // What the fit dropped, and what it never saw but misses the
            // fit: those vertices go, and the triangles that brought them.
            let kept: HashMap<u32, bool> = shared.iter().copied().zip(keep).collect();
            let dropped: std::collections::HashSet<u32> = vertices
                .iter()
                .copied()
                .filter(|v| {
                    !kept
                        .get(v)
                        .copied()
                        .unwrap_or_else(|| found.surface.distance_to(points[*v as usize]) <= flat)
                })
                .collect();
            if !dropped.is_empty() {
                let gathered = region.clone();
                region.retain(|&t| triangles[t].iter().all(|v| !dropped.contains(v)));
                if region.is_empty() {
                    for &t in &gathered[..first_sample.min(gathered.len())] {
                        tried[t] = true;
                        changed[t] = batch;
                    }
                    continue;
                }
                // What is left may have come apart where the dropped
                // triangles joined it. The largest piece is the claim; the
                // others stay free, to seed regions of their own or fall to
                // the planes.
                keep_largest_piece(&mut region, adjacency);
                let mut seen = std::collections::HashSet::new();
                vertices.clear();
                for &t in &region {
                    for &v in &triangles[t] {
                        if seen.insert(v) {
                            vertices.push(v);
                        }
                    }
                }
            }
            let mut shape = found.surface;
            // A facet of the region as it was first fitted, for telling a
            // flat face's facets from its own.
            let typical = {
                let mut areas: Vec<f64> = region.iter().map(|&t| patches.facet[t]).collect();
                areas.sort_by(f64::total_cmp);
                areas.get(areas.len() / 2).copied().unwrap_or(0.0)
            };
            let mut mine: std::collections::HashSet<usize> = region.iter().copied().collect();
            let mut seen: std::collections::HashSet<u32> = vertices.iter().copied().collect();

            // Then across any smooth edge, while the surface holds. A fit
            // from a small patch extrapolates only so far; when the growth
            // stalls with vertices gained since the last fit, the surface is
            // refitted to everything held and the rim tried again.
            let mut fitted_at = vertices.len();
            loop {
                let mut i = 0;
                while i < region.len() {
                    let t = region[i];
                    i += 1;
                    for h in 3 * t..3 * t + 3 {
                        let Some(g) = adjacency.twin[h] else {
                            continue;
                        };
                        let other = g / 3;
                        if mine.contains(&other) || groups.of[other] != usize::MAX {
                            continue;
                        }
                        let corners = triangles[other].map(|v| points[v as usize]);
                        if corners.iter().any(|p| shape.distance_to(*p) > flat) {
                            continue;
                        }
                        // Across a smooth edge, the facet's normal agrees
                        // with the surface's. Across a sharper one (a coarse
                        // mesh spanning two rows of a small fillet in one
                        // triangle), the surface must account for the whole
                        // lean.
                        if mesh.smooth(h) {
                            let centroid = Point::from_vector(
                                (corners[0].to_vector()
                                    + corners[1].to_vector()
                                    + corners[2].to_vector())
                                    / 3.0,
                            );
                            let direction = gradient(&shape, centroid);
                            let m = direction.magnitude();
                            if m == 0.0 || (direction.dot(mesh.normals[other]) / m).abs() < agree {
                                continue;
                            }
                        } else if !leans_as_the_surface(&shape, corners, mesh.normals[other])
                            || !sags_as_the_surface(&shape, corners, flat)
                        {
                            continue;
                        }
                        // A facet of a flat face much larger than the
                        // region's own (a plane tangent to it) comes with its
                        // whole face or not at all: its corners by the
                        // tangent line lie on the surface, its far ones not.
                        let patch = patches.of[other];
                        let taken: Vec<usize> = if patches.area[patch] > typical * PATCH_SCALE {
                            let whole = &patches.members[patch];
                            let on = whole.iter().all(|&t| {
                                mine.contains(&t)
                                    || (groups.of[t] == usize::MAX
                                        && triangles[t].iter().all(|&v| {
                                            shape.distance_to(points[v as usize]) <= flat
                                        }))
                            });
                            if !on {
                                continue;
                            }
                            whole
                                .iter()
                                .copied()
                                .filter(|t| !mine.contains(t))
                                .collect()
                        } else {
                            vec![other]
                        };
                        for other in taken {
                            mine.insert(other);
                            region.push(other);
                            for &v in &triangles[other] {
                                if seen.insert(v) {
                                    vertices.push(v);
                                }
                            }
                        }
                    }
                }
                if vertices.len() <= fitted_at {
                    break;
                }
                fitted_at = vertices.len();
                let (pts, nrm) = mesh.samples(&vertices, &region);
                match recognize_curved(&pts, &nrm, &mesh.chords(&region), flat, tol) {
                    Some(better) => shape = better.surface,
                    None => break,
                }
            }
            // A few triangles the region surrounds on every side, with their
            // corners on the surface, belong to it whatever their normals:
            // a sliver's plane through three nearly collinear points on the
            // surface tilts far from the surface's normal, as in the fans
            // round a sphere's pole, where several such slivers touch.
            let mut claimed: std::collections::HashSet<usize> = std::collections::HashSet::new();
            let rim: Vec<usize> = region
                .iter()
                .flat_map(|&t| (3 * t..3 * t + 3).filter_map(|h| adjacency.twin[h]))
                .map(|g| g / 3)
                .filter(|&other| !mine.contains(&other) && groups.of[other] == usize::MAX)
                .collect();
            for start in rim {
                if claimed.contains(&start) {
                    continue;
                }
                let mut cluster = vec![start];
                let mut inside: std::collections::HashSet<usize> =
                    std::collections::HashSet::from([start]);
                let mut surrounded = true;
                let mut i = 0;
                while i < cluster.len() && surrounded {
                    let t = cluster[i];
                    i += 1;
                    for h in 3 * t..3 * t + 3 {
                        let Some(g) = adjacency.twin[h] else {
                            surrounded = false;
                            break;
                        };
                        let next = g / 3;
                        if mine.contains(&next) || inside.contains(&next) {
                            continue;
                        }
                        if groups.of[next] != usize::MAX || cluster.len() >= ENCLOSED_CLUSTER {
                            surrounded = false;
                            break;
                        }
                        inside.insert(next);
                        cluster.push(next);
                    }
                }
                let on_surface = cluster.iter().all(|&t| {
                    let corners = triangles[t].map(|v| points[v as usize]);
                    corners.iter().all(|p| shape.distance_to(*p) <= flat)
                        && sags_as_the_surface(&shape, corners, flat)
                        && (is_sliver(corners)
                            || leans_as_the_surface(&shape, corners, mesh.normals[t]))
                });
                if surrounded && on_surface {
                    for t in cluster {
                        claimed.insert(t);
                        if mine.insert(t) {
                            region.push(t);
                        }
                    }
                }
            }
            // The vertices on the surface do not put the triangles on it: a
            // long facet across a flat stretch has its corners where a
            // surface through both ends of the stretch passes, and its
            // middle far from it. Where that facet is a large flat face's (a
            // plane meeting the surface along a tangent circle, triangulated
            // across it), it is the face's and is peeled off, the largest
            // piece left being the region; any other such facet says the
            // surface is wrong.
            let peel: Vec<usize> = region
                .iter()
                .copied()
                .filter(|&t| {
                    patches.area[patches.of[t]] > typical * PATCH_SCALE
                        && !sags_as_the_surface(
                            &shape,
                            triangles[t].map(|v| points[v as usize]),
                            flat,
                        )
                })
                .collect();
            if !peel.is_empty() {
                region.retain(|t| !peel.contains(t));
                keep_largest_piece(&mut region, adjacency);
                let mut seen = std::collections::HashSet::new();
                vertices.clear();
                for &t in &region {
                    for &v in &triangles[t] {
                        if seen.insert(v) {
                            vertices.push(v);
                        }
                    }
                }
                // What is left is asked again what it is: a band a row or
                // two high lies on a whole family of surfaces, and the one
                // chosen with the flat facets in is not the one without.
                // Where nothing is found, the surface held so far stands if
                // it still holds what is left: the flat facets gone, the
                // normals at the rim they shared lean to one side, and can
                // mislead the fit that the surface already answers.
                let (pts, nrm) = mesh.samples(&vertices, &region);
                match recognize_curved(&pts, &nrm, &mesh.chords(&region), flat, tol) {
                    Some(better) => shape = better.surface,
                    None if worst_deviation(&shape, &pts) <= flat => {}
                    None => {
                        for &t in &region {
                            tried[t] = true;
                            changed[t] = batch;
                        }
                        continue;
                    }
                }
            }
            let spans = |t: usize| {
                !sags_as_the_surface(&shape, triangles[t].map(|v| points[v as usize]), flat)
            };
            let spans_off = region.iter().any(|&t| spans(t));
            let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
            let deviation = worst_deviation(&shape, &pts);
            let flat_too = crate::recognize::is_flat(&pts, flat, tol);
            // A sphere or a torus through a handful of noisy vertices is
            // as likely the noise's as the part's: twice its unknowns are
            // asked for.
            let few = match shape {
                Canonical::Sphere(_) => vertices.len() < 8,
                Canonical::Torus(_) => vertices.len() < 14,
                _ => false,
            };
            if deviation > flat || flat_too || spans_off || region.len() < 2 || few {
                for &t in &region {
                    tried[t] = true;
                    changed[t] = batch;
                }
                continue;
            }
            let g = groups.carriers.len();
            for &t in &region {
                groups.of[t] = g;
                changed[t] = batch;
            }
            groups.carriers.push(Carrier::Curved(Curved {
                shape,
                deviation,
                fitted: deviation,
                centre: (0.0, 0.0),
                wraps: false,
                wraps_v: false,
                fixed: false,
                vertices,
                patch: None,
            }));
        }
    }
}

/// The smooth regions recognition left free, rebuilt as extrusions or
/// surfaces of revolution where one fits.
///
/// A region is the free triangles joined across edges that turn but do not
/// crease. Its vertices, with the normals averaged over its triangles, are
/// fitted as an extrusion and then as a surface of revolution, each held
/// to the coplanar distance at every vertex; a region that goes all the way
/// round (a closed profile, a whole turn) is left to its facets, as is one
/// whose triangles do not sag as the surface does.
fn swept_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) {
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    let cos_crease = options.crease.cos();
    let turn = |h: Half| adjacency.twin[h].map(|g| normals[h / 3].dot(normals[g / 3]));
    for region in smooth_regions(triangles, adjacency, &normals, options, groups) {
        // The recognized bands the region runs into smoothly, offered with
        // it: a stretch of a wavy profile straight enough to pass for a
        // cone is part of the one sweep, and the sweep takes it where the
        // whole fits.
        let mut bands: Vec<usize> = Vec::new();
        for &t in &region {
            for h in 3 * t..3 * t + 3 {
                let (Some(c), Some(g)) = (turn(h), adjacency.twin[h]) else {
                    continue;
                };
                let other = groups.of[g / 3];
                if c >= cos_crease
                    && other != usize::MAX
                    && matches!(groups.carriers.get(other), Some(Carrier::Curved(_)))
                    && !bands.contains(&other)
                {
                    bands.push(other);
                }
            }
        }
        if !bands.is_empty() {
            let mut widened = region.clone();
            widened.extend((0..triangles.len()).filter(|&t| bands.contains(&groups.of[t])));
            if let Some(claim) = swept_claim(points, triangles, &normals, &widened, flat, tol) {
                for &b in &bands {
                    groups.carriers[b] = Carrier::Gone;
                }
                claim_swept(groups, &widened, claim);
                continue;
            }
        }
        if let Some(claim) = swept_claim(points, triangles, &normals, &region, flat, tol) {
            claim_swept(groups, &region, claim);
        }
    }
}

/// The free triangles joined across edges that turn but do not crease,
/// each such region holding at least [`SWEPT_TRIANGLES`] and an edge that
/// turns: the smooth regions no canonical surface took.
fn smooth_regions(
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    normals: &[Vector],
    options: &MeshSolidOptions,
    groups: &Groups,
) -> Vec<Vec<usize>> {
    let (cos_crease, cos_flat) = (options.crease.cos(), options.coplanar_angle.cos());
    let turn = |h: Half| adjacency.twin[h].map(|g| normals[h / 3].dot(normals[g / 3]));
    let mut seen = vec![false; triangles.len()];
    let mut regions = Vec::new();
    for seed in 0..triangles.len() {
        if seen[seed] || groups.of[seed] != usize::MAX {
            continue;
        }
        let mut region = vec![seed];
        seen[seed] = true;
        let mut bends = false;
        let mut i = 0;
        while i < region.len() {
            let t = region[i];
            i += 1;
            for h in 3 * t..3 * t + 3 {
                let (Some(c), Some(g)) = (turn(h), adjacency.twin[h]) else {
                    continue;
                };
                if c < cos_crease {
                    continue;
                }
                let other = g / 3;
                if c < cos_flat {
                    bends = true;
                }
                if !seen[other] && groups.of[other] == usize::MAX {
                    seen[other] = true;
                    region.push(other);
                }
            }
        }
        if bends && region.len() >= SWEPT_TRIANGLES {
            regions.push(region);
        }
    }
    regions
}

/// A sweep fitted to a region's vertices, with the normals averaged over
/// its triangles: held at every vertex, its triangles sagging as it does,
/// and its chart centre the vertices' mean. `None` where none fits or the
/// region goes all the way round.
fn swept_claim(
    points: &[Point],
    triangles: &[[u32; 3]],
    normals: &[Vector],
    region: &[usize],
    flat: f64,
    tol: Tolerances,
) -> Option<Curved> {
    let mut sums: HashMap<u32, Vector> = HashMap::new();
    for &t in region {
        for &v in &triangles[t] {
            *sums.entry(v).or_insert(Vector::ZERO) += normals[t];
        }
    }
    let mut vertices: Vec<u32> = sums.keys().copied().collect();
    vertices.sort_unstable();
    let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
    let nrm: Vec<Vector> = vertices
        .iter()
        .map(|v| {
            let s = sums[v];
            let m = s.magnitude();
            if m > 0.0 { s / m } else { Vector::Z }
        })
        .collect();
    let found = crate::recognize_swept::fit_extrusion(&pts, &nrm, flat, tol)
        .or_else(|| crate::recognize_swept::fit_revolution(&pts, &nrm, flat, tol))?;
    let shape = Canonical::Swept(Box::new(crate::recognize::SweptShape::new(
        found.surface,
        tol,
    )));
    if !region
        .iter()
        .all(|&t| sags_as_the_surface(&shape, triangles[t].map(|v| points[v as usize]), flat))
    {
        return None;
    }
    let charts: Vec<(f64, f64)> = pts.iter().filter_map(|p| chart(&shape, *p, tol)).collect();
    if charts.len() != pts.len() {
        return None;
    }
    let (pu, pv) = periodic(&shape);
    // Each chart direction's centre, and whether the region goes all the
    // way round it: then its centre is half a turn from the seam, as a
    // canonical band's is.
    let centre_of = |mut values: Vec<f64>, periodic: bool| -> (f64, bool) {
        if periodic {
            let (mean, gap) = angular_spread(&mut values);
            if gap < core::f64::consts::FRAC_PI_2 {
                (core::f64::consts::PI, true)
            } else {
                (mean, false)
            }
        } else {
            #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
            (
                values.iter().sum::<f64>() / values.len().max(1) as f64,
                false,
            )
        }
    };
    let (cu, wraps) = centre_of(charts.iter().map(|c| c.0).collect(), pu);
    let (cv, wraps_v) = centre_of(charts.iter().map(|c| c.1).collect(), pv);
    Some(Curved {
        shape,
        deviation: found.deviation,
        fitted: found.deviation,
        centre: (cu, cv),
        wraps,
        wraps_v,
        fixed: true,
        vertices,
        patch: None,
    })
}

/// The smooth regions recognition and the sweeps left free, each rebuilt
/// on one fitted B-spline patch where the region is one disk and the patch
/// verifies.
///
/// A region is the free triangles joined across edges that turn but do not
/// crease, as for the sweeps. The square its chart maps onto takes its
/// corners first where the face across the region's boundary changes; the
/// faces there are read as the planar pass would make them, every region a
/// face of its own and the other free triangles grouped into planes. A
/// region refused a patch keeps its triangles for the planes, and is
/// counted by why.
fn patch_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    let regions = smooth_regions(triangles, adjacency, &normals, options, groups);
    if regions.is_empty() {
        return Ok(());
    }
    let mut across = groups.clone();
    for region in &regions {
        let g = across.carriers.len();
        across.carriers.push(Carrier::Gone);
        for &t in region {
            across.of[t] = g;
        }
    }
    coplanar_groups(
        points,
        triangles,
        adjacency,
        options,
        flat,
        &mut across,
        tol,
    )?;
    for region in regions {
        // The faces across the boundary on either side of each boundary
        // vertex: a corner where they differ.
        let own = across.of[region[0]];
        let mut sides: HashMap<u32, [Vec<usize>; 2]> = HashMap::new();
        for &t in &region {
            for h in 3 * t..3 * t + 3 {
                let other = adjacency.twin[h].map_or(usize::MAX, |g| across.of[g / 3]);
                if other == own {
                    continue;
                }
                let (a, b) = from_to(triangles, h);
                sides.entry(a).or_default()[0].push(other);
                sides.entry(b).or_default()[1].push(other);
            }
        }
        let corners: std::collections::HashSet<u32> = sides
            .into_iter()
            .filter(|(_, [leaving, arriving])| {
                let mut faces = leaving.iter().chain(arriving);
                faces.next().is_some_and(|first| faces.any(|f| f != first))
            })
            .map(|(v, _)| v)
            .collect();
        let found = crate::recognize_patch::fit_patch(
            &crate::recognize_patch::Region {
                points,
                triangles,
                members: &region,
                corners: &corners,
            },
            flat,
            tol,
        );
        match found {
            Ok(patch) => {
                let mut vertices: Vec<u32> = region.iter().flat_map(|&t| triangles[t]).collect();
                vertices.sort_unstable();
                vertices.dedup();
                let claim = Curved {
                    shape: Canonical::Swept(Box::new(crate::recognize::SweptShape::new(
                        patch.surface,
                        tol,
                    ))),
                    deviation: patch.deviation,
                    fitted: patch.deviation,
                    centre: (0.5, 0.5),
                    wraps: false,
                    wraps_v: false,
                    fixed: true,
                    vertices,
                    patch: Some(patch.mapped),
                };
                claim_swept(groups, &region, claim);
            }
            Err(crate::recognize_patch::Refused::NotDisk) => groups.refused.not_disk += 1,
            Err(crate::recognize_patch::Refused::Narrow) => groups.refused.narrow += 1,
            Err(crate::recognize_patch::Refused::Unverified) => groups.refused.unverified += 1,
        }
    }
    Ok(())
}

/// A region's triangles given to a new swept face.
fn claim_swept(groups: &mut Groups, region: &[usize], claim: Curved) {
    let g = groups.carriers.len();
    for &t in region {
        groups.of[t] = g;
    }
    groups.carriers.push(Carrier::Curved(claim));
}

/// The fewest triangles a smooth region needs to be tried as a sweep.
const SWEPT_TRIANGLES: usize = 32;

/// The largest distance between two of the points, from the first.
fn span(points: &[Point]) -> f64 {
    points.first().map_or(0.0, |a| {
        points.iter().map(|p| p.distance(*a)).fold(0.0, f64::max)
    })
}

/// Turn each recognized sphere whose boundary is circles in parallel planes
/// about their common normal, so the circles are latitudes and the face a
/// cap or a zone about a pole; the axis points into a cap, so its pole is
/// the frame's north one. A sphere with no boundary is whole and keeps its
/// frame, and one bounded otherwise keeps it for the re-framing that moves
/// its poles away.
fn sphere_axes(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    flat: f64,
    tol: Tolerances,
) {
    let reach = flat * REACH;
    for g in 0..groups.carriers.len() {
        let Carrier::Curved(curved) = &groups.carriers[g] else {
            continue;
        };
        let Canonical::Sphere(sphere) = curved.shape else {
            continue;
        };
        let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
            continue;
        };
        let mut axis: Option<Vector> = None;
        let mut planar = true;
        for ring in &loops {
            let pts: Vec<Point> = ring.iter().map(|&v| points[v as usize]).collect();
            let Some((through, normal)) =
                (pts.len() >= 3).then(|| plane_through(&pts, tol)).flatten()
            else {
                planar = false;
                break;
            };
            let n = normal.vector();
            if pts.iter().any(|p| (*p - through).dot(n).abs() > reach) {
                planar = false;
                break;
            }
            match axis {
                None => axis = Some(n),
                Some(a) if a.cross(n).magnitude() <= 1e-3 => {}
                Some(_) => {
                    planar = false;
                    break;
                }
            }
        }
        let (true, Some(mut z)) = (planar, axis) else {
            continue;
        };
        // Into the region: the pole a cap covers is the north one.
        let side: f64 = curved
            .vertices
            .iter()
            .map(|&v| (points[v as usize] - sphere.centre()).dot(z))
            .sum();
        if side < 0.0 {
            z = -z;
        }
        let Ok(z) = Direction::new(z, tol) else {
            continue;
        };
        let Ok(frame) = Frame::new(sphere.centre(), z, z.any_perpendicular(), tol) else {
            continue;
        };
        let Ok(turned) = Sphere::new(frame, sphere.radius(), tol) else {
            continue;
        };
        if let Carrier::Curved(curved) = &mut groups.carriers[g] {
            curved.shape = Canonical::Sphere(turned);
            curved.fixed = true;
        }
    }
}

/// A region's boundary, walked into loops of mesh vertices: each border
/// half-edge leads to the one leaving its end. `None` where a vertex has
/// more than one way on.
fn border_loops(
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    of: &[usize],
    g: usize,
) -> Option<Vec<Vec<u32>>> {
    let mut leaving: HashMap<u32, Vec<u32>> = HashMap::new();
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
    let mut done: std::collections::HashSet<u32> = std::collections::HashSet::new();
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
/// surface's face can carry the rings as holes.
fn hole_frames(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    tol: Tolerances,
) {
    for g in 0..groups.carriers.len() {
        let Carrier::Curved(curved) = &groups.carriers[g] else {
            continue;
        };
        let wanted = match curved.shape {
            Canonical::Sphere(_) => curved.wraps && !curved.fixed,
            Canonical::Torus(_) => curved.wraps && curved.wraps_v,
            _ => false,
        };
        if !wanted {
            continue;
        }
        let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
            continue;
        };
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
                let on_region = |p: Point| {
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
                        .is_some_and(|(_, t)| groups.of[t] == g)
                };
                let best = ranked.into_iter().find(|(_, z)| {
                    rings.iter().all(|ring| turns_about(ring, *z) == 0)
                        && on_region(sphere.centre() + *z * sphere.radius())
                        && on_region(sphere.centre() - *z * sphere.radius())
                });
                let Some((nearest, z)) = best else {
                    continue;
                };
                // A pole within a few degrees of a ring has no room round it.
                if nearest > POLE_CLEARANCE.cos() {
                    continue;
                }
                let Ok(z) = Direction::new(z, tol) else {
                    continue;
                };
                let Ok(frame) = Frame::new(sphere.centre(), z, z.any_perpendicular(), tol) else {
                    continue;
                };
                let Ok(turned) = Sphere::new(frame, sphere.radius(), tol) else {
                    continue;
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
            continue;
        };
        let Some(frame) = axis_frame(&shape) else {
            continue;
        };
        let (x, y) = (frame.x().vector(), frame.y().vector());
        let Ok(x) = Direction::new(x * free.cos() + y * free.sin(), tol) else {
            continue;
        };
        let Ok(turned) = Frame::new(frame.origin(), frame.z(), x, tol) else {
            continue;
        };
        let Some(shape) = on_frame(&shape, turned, tol) else {
            continue;
        };
        // A torus's seam round its axis is a parallel, placed at the tube
        // angle the rings leave widest free; its chart is centred half a
        // turn on from it.
        let centre_v = match shape {
            Canonical::Torus(_) => {
                let across: Vec<Vec<(f64, f64)>> =
                    vec![angles[0].iter().map(|&(u, v)| (v, u)).collect()];
                let Some(free_v) = free_angle(&across) else {
                    continue;
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
fn slit_bands(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    tol: Tolerances,
) {
    for g in 0..groups.carriers.len() {
        let Carrier::Curved(curved) = &groups.carriers[g] else {
            continue;
        };
        if !curved.wraps
            || curved.wraps_v
            || !matches!(curved.shape, Canonical::Cylinder(_) | Canonical::Cone(_))
        {
            continue;
        }
        let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
            continue;
        };
        let [ring] = &loops[..] else {
            continue;
        };
        let Some(frame) = axis_frame(&curved.shape) else {
            continue;
        };
        let directions: Vec<Vector> = ring
            .iter()
            .map(|&v| points[v as usize] - frame.origin())
            .collect();
        if turns_about(&directions, frame.z().vector()) != 0 {
            continue;
        }
        let mut angles: Vec<f64> = curved
            .vertices
            .iter()
            .filter_map(|&v| chart(&curved.shape, points[v as usize], tol).map(|c| c.0))
            .collect();
        let Some(gap) = widest_gap(&mut angles) else {
            continue;
        };
        if let Carrier::Curved(curved) = &mut groups.carriers[g] {
            curved.wraps = false;
            curved.centre = (
                ogeom_math::elementary::wrap_angle(gap + core::f64::consts::PI),
                curved.centre.1,
            );
        }
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
fn align_axes(
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
        let Some(frame) = axis_frame(&curved.shape) else {
            continue;
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

/// A vertex as built: a mesh vertex, or a point the construction placed,
/// where a closed circle is bounded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Corner {
    Mesh(u32),
    Placed(usize),
}

/// An edge as planned.
struct EdgeSpec {
    curve: Curve,
    range: (f64, f64),
    ends: [Corner; 2],
    tolerance: f64,
    closed_circle: bool,
}

/// Everything the build needs, decided before anything is built.
struct Plan {
    edges: Vec<EdgeSpec>,
    /// Each boundary mesh edge, by its lower vertex first: its planned
    /// edge, and whether walking it from the lower vertex runs the edge
    /// forward.
    edge_of: HashMap<(u32, u32), (usize, bool)>,
    /// Each face's loops as half-edges, walked with the face on the left.
    loops: Vec<Vec<Vec<Half>>>,
    placed: Vec<Point>,
    /// Each curved face's surface, windowed to its boundary.
    surfaces: Vec<Option<ogeom_geom::SurfaceGeometry>>,
    /// Each edge's image on each curved face it bounds, by (edge, face),
    /// with how far the image strays from the edge.
    pcurves: HashMap<(usize, usize), (PlanarCurve, f64)>,
    /// How each curved face's boundary lies in its chart.
    layouts: Vec<Layout>,
    /// How many free edges of curved faces are fitted curves.
    free_fitted: usize,
}

/// Plans the edges and loops.
struct Planner<'a> {
    points: &'a [Point],
    triangles: &'a [[u32; 3]],
    adjacency: &'a Adjacency,
    groups: &'a Groups,
    merge: bool,
    /// Vertices kept as edge ends whatever lies either side of them.
    pinned: &'a std::collections::HashSet<u32>,
    /// Mesh edges whose seam is threaded straight through the chain's own
    /// vertices: the solved curve bent the trim of a face beside it back
    /// across itself.
    straight: &'a std::collections::HashSet<(u32, u32)>,
    /// The turn, in radians, from which a free boundary of a curved face
    /// is cut into separate edges at a vertex.
    crease: f64,
    flat: f64,
    tol: Tolerances,
    /// Curves snapped by earlier plans with the same points and distance.
    snaps: &'a std::sync::Mutex<SnapCache>,
}

/// A chain of boundary vertices between two kept ones, and its mesh edges.
type Chain = (Vec<u32>, Vec<(u32, u32)>);

/// Why a plan was refused: curved faces to facet, or vertices to keep.
enum Replan {
    Facet(Vec<usize>),
    Pin(Vec<u32>),
}

/// Whether a triangle with its corners on `shape` stands off it at its
/// centroid no farther than a facet of its size over that much turn of the
/// surface would: its longest side times the turn of the surface's normals
/// across its corners, over six. A triangle spanning a recess, its
/// corners on the rim and its middle over the hollow, stands off more.
fn sags_as_the_surface(shape: &Canonical, corners: [Point; 3], flat: f64) -> bool {
    let unit = |p: Point| {
        let g = gradient(shape, p);
        let m = g.magnitude();
        (m > 0.0).then(|| g / m)
    };
    let (Some(a), Some(b), Some(c)) = (unit(corners[0]), unit(corners[1]), unit(corners[2])) else {
        return false;
    };
    let angle = |x: Vector, y: Vector| x.dot(y).clamp(-1.0, 1.0).acos();
    let turn = angle(a, b).max(angle(b, c)).max(angle(a, c));
    let longest = corners[0]
        .distance(corners[1])
        .max(corners[1].distance(corners[2]))
        .max(corners[0].distance(corners[2]));
    let centroid = Point::from_vector(
        (corners[0].to_vector() + corners[1].to_vector() + corners[2].to_vector()) / 3.0,
    );
    shape.distance_to(centroid) <= longest * turn / 6.0 + flat
}

/// Whether a triangle is a sliver: no taller across its longest side than
/// a tenth of that side, so its normal says little.
fn is_sliver(corners: [Point; 3]) -> bool {
    let [a, b, c] = corners;
    let sides = [(a, b, c), (b, c, a), (c, a, b)];
    let (p, q, r) = sides
        .into_iter()
        .max_by(|x, y| x.0.distance(x.1).total_cmp(&y.0.distance(y.1)))
        .unwrap_or((a, b, c));
    let base = p.distance(q);
    base > 0.0 && distance_to_line(r, p, q) <= base * 0.1
}

/// How far a chord taken for an edge may stand off the curved face it
/// bounds, against its length: a twentieth.
const CHORD_SAG: f64 = 0.05;

/// How much wider than the gap it was measured from an edge's tolerance is
/// recorded: a millionth.
const TOLERANCE_MARGIN: f64 = 1e-6;

/// The most triangles a cluster the region surrounds may hold and still
/// be taken into it whatever its normals.
const ENCLOSED_CLUSTER: usize = 16;

/// The least coplanar distance recognition holds to, against the mesh's
/// measured scatter: twice it, for a seam between two surfaces each fitted
/// to scattered vertices, and a quarter more for the tail its ninetieth
/// percentile leaves out.
const NOISE_FLOOR: f64 = 2.5;

/// How far a snapped curve may stand off its chain or its faces.
const REACH: f64 = 20.0;

impl Planner<'_> {
    fn border(&self, h: Half) -> bool {
        match self.adjacency.twin[h] {
            None => true,
            Some(g) => self.groups.of[g / 3] != self.groups.of[h / 3],
        }
    }

    fn curved(&self, g: usize) -> Option<&Curved> {
        match &self.groups.carriers[g] {
            Carrier::Curved(c) => Some(c),
            _ => None,
        }
    }

    /// The plan, or the curved faces whose boundary could not be built
    /// exactly and are to be faceted instead.
    #[allow(clippy::too_many_lines, reason = "one pass over the boundary")]
    fn plan(&self) -> OgeomResult<Result<Plan, Replan>> {
        let halves = self.triangles.len() * 3;
        let mut edge_faces: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for h in 0..halves {
            if self.border(h) {
                let (a, b) = from_to(self.triangles, h);
                edge_faces
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push(self.groups.of[h / 3]);
            }
        }
        for faces in edge_faces.values_mut() {
            faces.sort_unstable();
        }
        let mut incident: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
        for &(a, b) in edge_faces.keys() {
            incident.entry(a).or_default().push((a, b));
            incident.entry(b).or_default().push((a, b));
        }
        let any_curved = |faces: &[usize]| faces.iter().any(|&g| self.curved(g).is_some());
        // A vertex between two boundary edges with the same faces either
        // side is inside one edge: between planes only where the two are
        // on one line; along a curved face always, the curve deciding,
        // but for a free boundary (one face, nothing across), which ends
        // an edge where it turns by the crease angle or more.
        let removable = |v: u32| -> bool {
            if !self.merge || self.pinned.contains(&v) {
                return false;
            }
            let Some(list) = incident.get(&v) else {
                return false;
            };
            let [e1, e2] = list[..] else {
                return false;
            };
            let faces = &edge_faces[&e1];
            if *faces != edge_faces[&e2] {
                return false;
            }
            let far = |(a, b): (u32, u32)| if a == v { b } else { a };
            if any_curved(faces) {
                if faces.len() > 1 {
                    return true;
                }
                let at = self.points[v as usize];
                let (inward, onward) = (
                    at - self.points[far(e1) as usize],
                    self.points[far(e2) as usize] - at,
                );
                let lengths = inward.magnitude() * onward.magnitude();
                return lengths > 0.0 && inward.dot(onward) / lengths > self.crease.cos();
            }
            let (p, q) = (self.points[far(e1) as usize], self.points[far(e2) as usize]);
            let at = self.points[v as usize];
            (at - p).dot(q - at) > 0.0 && distance_to_line(at, p, q) <= self.flat
        };
        let is_kept: HashMap<u32, bool> = incident.keys().map(|&v| (v, !removable(v))).collect();

        let mut plan = Plan {
            edges: Vec::new(),
            edge_of: HashMap::new(),
            loops: vec![Vec::new(); self.groups.carriers.len()],
            placed: Vec::new(),
            surfaces: vec![None; self.groups.carriers.len()],
            pcurves: HashMap::new(),
            layouts: vec![Layout::Open; self.groups.carriers.len()],
            free_fitted: 0,
        };
        let mut failed: Vec<usize> = Vec::new();
        let mut keys: Vec<(u32, u32)> = edge_faces.keys().copied().collect();
        keys.sort_unstable();
        // The chain of boundary edges from `key` to the next kept vertex,
        // with its edges; the first pass starts only from a kept vertex.
        let walk = |key: (u32, u32), pass: usize| -> Option<Chain> {
            let (a, b) = key;
            let start = if is_kept[&a] {
                a
            } else if is_kept[&b] {
                b
            } else if pass == 1 {
                a
            } else {
                return None;
            };
            let mut chain = vec![start];
            let mut edges = vec![key];
            let mut at = if start == a { b } else { a };
            chain.push(at);
            while !is_kept[&at] && at != start {
                let Some(&following) = incident[&at].iter().find(|e| !edges.contains(e)) else {
                    break;
                };
                edges.push(following);
                at = if following.0 == at {
                    following.1
                } else {
                    following.0
                };
                chain.push(at);
            }
            Some((chain, edges))
        };
        // Every chain is walked as below, and the curves no earlier plan
        // snapped are solved side by side before the edges are planned in
        // order.
        {
            let mut taken: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
            let mut wanted: Vec<(Vec<u32>, Vec<usize>)> = Vec::new();
            for pass in 0..2 {
                for &key in &keys {
                    if taken.contains(&key) {
                        continue;
                    }
                    let Some((chain, edges)) = walk(key, pass) else {
                        continue;
                    };
                    let faces = &edge_faces[&key];
                    if any_curved(faces)
                        && !edges.iter().any(|e| self.straight.contains(e))
                        && !self.snap_known(&chain, faces)
                    {
                        wanted.push((chain, faces.clone()));
                    }
                    taken.extend(edges);
                }
            }
            self.snap_all(wanted);
        }
        for pass in 0..2 {
            for &key in &keys {
                if plan.edge_of.contains_key(&key) {
                    continue;
                }
                let Some((chain, edges)) = walk(key, pass) else {
                    continue;
                };
                let (start, at) = (chain[0], chain[chain.len() - 1]);
                let faces = edge_faces[&key].clone();
                if any_curved(&faces) {
                    let threaded =
                        faces.len() > 1 && edges.iter().any(|e| self.straight.contains(e));
                    let snapped = if threaded {
                        let closed = chain.len() > 2 && chain[0] == chain[chain.len() - 1];
                        let pts: Vec<Point> = chain[..chain.len() - usize::from(closed)]
                            .iter()
                            .map(|&v| self.points[v as usize])
                            .collect();
                        self.chord(&pts, closed, &faces, self.flat * REACH)
                    } else {
                        self.snapped(&chain, &faces)
                    };
                    match snapped {
                        Some(spec) => {
                            let index = plan.edges.len();
                            let (spec, forward, images) = spec;
                            if faces.len() == 1 && !images.is_empty() {
                                plan.free_fitted += 1;
                            }
                            for (g, pcurve, deviation) in images {
                                plan.pcurves.insert((index, g), (pcurve, deviation));
                            }
                            let spec = match spec {
                                Snapped::Open(curve, range, tolerance) => EdgeSpec {
                                    curve,
                                    range,
                                    ends: if forward {
                                        [
                                            Corner::Mesh(chain[0]),
                                            Corner::Mesh(*chain.last().unwrap_or(&chain[0])),
                                        ]
                                    } else {
                                        [
                                            Corner::Mesh(*chain.last().unwrap_or(&chain[0])),
                                            Corner::Mesh(chain[0]),
                                        ]
                                    },
                                    tolerance,
                                    closed_circle: false,
                                },
                                Snapped::Loop(curve, range, tolerance) => EdgeSpec {
                                    curve,
                                    range,
                                    ends: [Corner::Mesh(chain[0]), Corner::Mesh(chain[0])],
                                    tolerance,
                                    closed_circle: false,
                                },
                                Snapped::Closed(curve, tolerance) => {
                                    use ogeom_geom::Curve3d as _;
                                    let at = curve.point_at(0.0, self.tol)?;
                                    plan.placed.push(at);
                                    let corner = Corner::Placed(plan.placed.len() - 1);
                                    EdgeSpec {
                                        curve,
                                        range: (0.0, core::f64::consts::TAU),
                                        ends: [corner, corner],
                                        tolerance,
                                        closed_circle: true,
                                    }
                                }
                            };
                            plan.edges.push(spec);
                            for (i, e) in edges.iter().enumerate() {
                                plan.edge_of
                                    .insert(*e, (index, (chain[i] == e.0) == forward));
                            }
                        }
                        None => {
                            for g in faces {
                                if self.curved(g).is_some() && !failed.contains(&g) {
                                    failed.push(g);
                                }
                            }
                            for e in edges {
                                plan.edge_of.insert(e, (usize::MAX, true));
                            }
                        }
                    }
                    continue;
                }
                let end = at;
                let (p, q) = (self.points[start as usize], self.points[end as usize]);
                let straight = start != end
                    && chain[1..chain.len() - 1]
                        .iter()
                        .all(|&v| distance_to_line(self.points[v as usize], p, q) <= self.flat);
                // Each piece: its vertices in order, and its mesh edges.
                type Piece = (Vec<u32>, Vec<(u32, u32)>);
                let pieces: Vec<Piece> = if straight {
                    vec![(chain.clone(), edges.clone())]
                } else {
                    edges.iter().map(|&e| (vec![e.0, e.1], vec![e])).collect()
                };
                for (piece, piece_edges) in pieces {
                    let (from, to) = (piece[0], piece[piece.len() - 1]);
                    let (p, q) = (self.points[from as usize], self.points[to as usize]);
                    let mut reach = self.tol.confusion();
                    for &v in &piece {
                        let at = self.points[v as usize];
                        reach = reach.max(distance_to_line(at, p, q));
                        for &g in &faces {
                            if let Carrier::Plane(plane) = &self.groups.carriers[g] {
                                reach = reach.max(plane.signed_distance_to(at).abs());
                            }
                        }
                    }
                    let index = plan.edges.len();
                    plan.edges.push(EdgeSpec {
                        curve: LineCurve::segment(p, q, self.tol)?.into(),
                        range: (0.0, p.distance(q)),
                        ends: [Corner::Mesh(from), Corner::Mesh(to)],
                        tolerance: reach,
                        closed_circle: false,
                    });
                    for (i, e) in piece_edges.iter().enumerate() {
                        plan.edge_of.insert(*e, (index, piece[i] == e.0));
                    }
                }
            }
        }

        // Each face's loops, walked with the face on the left: from a
        // boundary half-edge to the next one at its end, turning through the
        // face's own triangles around the vertex, which keeps a loop that
        // touches itself at a vertex on its own side.
        let mut walked = vec![false; halves];
        for h in 0..halves {
            if walked[h] || !self.border(h) {
                continue;
            }
            let mut ring = Vec::new();
            let mut at = h;
            loop {
                walked[at] = true;
                ring.push(at);
                let mut step = next(at);
                let mut guard = 0;
                while !self.border(step) {
                    let Some(twin) = self.adjacency.twin[step] else {
                        break;
                    };
                    step = next(twin);
                    guard += 1;
                    if guard > halves {
                        ogeom_bail!(Construction, "a face's boundary does not close");
                    }
                }
                if step == h {
                    break;
                }
                if walked[step] {
                    ogeom_bail!(Construction, "a face's boundary runs into itself");
                }
                at = step;
            }
            plan.loops[self.groups.of[h / 3]].push(ring);
        }

        // A face whose boundary would run out and back along one line (a
        // strip of slivers merged into two straight edges between the same
        // two vertices) encloses nothing; its runs keep their vertices. Two
        // edges of which one is curved (a flat face cut from a ball by a
        // second plane: an arc and a line) enclose a face like any other.
        let mut pin: Vec<u32> = Vec::new();
        for (g, rings) in plan.loops.iter().enumerate() {
            if !matches!(self.groups.carriers[g], Carrier::Plane(_)) {
                continue;
            }
            for ring in rings {
                let mut entries: Vec<usize> = Vec::new();
                for &h in ring {
                    let (edge, _) = self.entry(&plan, h);
                    if entries.last() != Some(&edge) {
                        entries.push(edge);
                    }
                }
                if entries.len() > 1 && entries.first() == entries.last() {
                    entries.pop();
                }
                if entries.len() < 3
                    && entries.iter().all(|&e| {
                        e != usize::MAX
                            && !plan.edges[e].closed_circle
                            && matches!(plan.edges[e].curve, Curve::Line(_))
                    })
                {
                    for &h in ring {
                        let (a, b) = from_to(self.triangles, h);
                        for v in [a, b] {
                            if !pin.contains(&v) && !self.pinned.contains(&v) {
                                pin.push(v);
                            }
                        }
                    }
                }
            }
        }
        if !pin.is_empty() && failed.is_empty() {
            return Ok(Err(Replan::Pin(pin)));
        }

        // A face round its axis, or round a torus's tube, is built as a band
        // between two full circles joined by a seam; a sphere's cap as one
        // circle, a seam and the pole; a sphere or torus with no boundary
        // whole. A face round its axis between two rims of any other shape,
        // holed or not, gets a seam of its own between them.
        for (g, carrier) in self.groups.carriers.iter().enumerate() {
            let Carrier::Curved(curved) = carrier else {
                continue;
            };
            if failed.contains(&g) || !(curved.wraps || curved.wraps_v) {
                continue;
            }
            let sphere = matches!(curved.shape, Canonical::Sphere(_));
            let torus = matches!(curved.shape, Canonical::Torus(_));
            let rings = &plan.loops[g];
            let circles = rings.iter().all(|ring| {
                let first = self.entry(&plan, ring[0]);
                first.0 != usize::MAX
                    && plan.edges[first.0].closed_circle
                    && ring.iter().all(|&h| self.entry(&plan, h).0 == first.0)
            });
            let wrapped = || {
                let resolved = rings
                    .iter()
                    .all(|ring| ring.iter().all(|&h| self.entry(&plan, h).0 != usize::MAX));
                if sphere || !curved.wraps || curved.wraps_v || !resolved {
                    return None;
                }
                let windings: Option<Vec<i32>> = rings
                    .iter()
                    .map(|ring| winding(&curved.shape, ring, self.triangles, self.points, self.tol))
                    .collect();
                let windings = windings?;
                let rims = windings.iter().filter(|w| w.abs() == 1).count();
                let holes = windings.iter().filter(|w| **w == 0).count();
                (rims == 2 && rims + holes == windings.len()).then_some(Layout::Wrapped)
            };
            let holed = || {
                let closed_round = sphere || (torus && curved.wraps && curved.wraps_v);
                let resolved = rings
                    .iter()
                    .all(|ring| ring.iter().all(|&h| self.entry(&plan, h).0 != usize::MAX));
                if !closed_round || !resolved {
                    return false;
                }
                let tau = core::f64::consts::TAU;
                let (_, wraps_v) = periodic(&curved.shape);
                rings.iter().all(|ring| {
                    let turns =
                        windings(&curved.shape, ring, self.triangles, self.points, self.tol);
                    // Round neither way, and clear of the seams.
                    let Some((0, 0)) = turns else {
                        return false;
                    };
                    let Some(polygon) = hole_polygons(
                        &curved.shape,
                        &[ring.as_slice()],
                        self.triangles,
                        self.points,
                        self.tol,
                    )
                    .pop() else {
                        return false;
                    };
                    let clear = |values: &mut dyn Iterator<Item = f64>, seam: f64| {
                        let (lo, hi) = values
                            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| {
                                (a.min(x), b.max(x))
                            });
                        ((lo - seam) / tau).floor() == ((hi - seam) / tau).floor()
                    };
                    clear(&mut polygon.iter().map(|p| p.0), 0.0)
                        && (!wraps_v
                            || clear(&mut polygon.iter().map(|p| p.1), torus_seam_v(curved)))
                })
            };
            let layout = if rings.is_empty() {
                (sphere || torus).then_some(Layout::Whole)
            } else if holed() {
                Some(Layout::Holed)
            } else if (curved.wraps && curved.wraps_v) || !circles {
                wrapped()
            } else if sphere && rings.len() == 1 && curved.fixed {
                Some(Layout::Cap)
            } else if rings.len() == 2
                && (!sphere || curved.fixed)
                && !matches!(curved.shape, Canonical::Swept(_))
            {
                // A band's seam is a ruling, a meridian or a tube's circle;
                // a sweep's would be its profile, and it is laid out as
                // wrapped instead, its seam seated between its rims.
                Some(Layout::Band {
                    round_tube: curved.wraps_v,
                })
            } else {
                wrapped()
            };
            let layout = match layout {
                Some(Layout::Wrapped) => {
                    let rings = plan.loops[g].clone();
                    self.seat_seam(&mut plan, curved, &rings)
                        .then_some(Layout::Wrapped)
                }
                other => other,
            };
            match layout {
                Some(layout) => plan.layouts[g] = layout,
                None => failed.push(g),
            }
        }

        // Every curved face's images of its edges, held to the reach: the
        // straight image in the chart where a parallel or a ruling has one,
        // a fit by projection where it has not.
        let reach = self.flat * REACH;
        for (g, carrier) in self.groups.carriers.iter().enumerate() {
            let Carrier::Curved(curved) = carrier else {
                continue;
            };
            if failed.contains(&g) || plan.loops[g].is_empty() {
                continue;
            }
            let surface = surface_of(curved, self.points, self.tol)?;
            if matches!(
                plan.layouts[g],
                Layout::Open | Layout::Wrapped | Layout::Holed
            ) {
                let mut held = true;
                'rings: for ring in &plan.loops[g] {
                    for &h in ring {
                        let (edge, _) = self.entry(&plan, h);
                        if edge == usize::MAX {
                            held = false;
                            break 'rings;
                        }
                        if plan.pcurves.contains_key(&(edge, g)) {
                            continue;
                        }
                        let spec = &plan.edges[edge];
                        // A chord taken for an edge is imaged as loosely as
                        // it stands off the face.
                        let reach = reach.max(spec.tolerance * 2.0);
                        let image =
                            image_on(curved, &surface, &spec.curve, spec.range, reach, self.tol);
                        match image {
                            Some(found) => {
                                plan.pcurves.insert((edge, g), found);
                            }
                            None => {
                                held = false;
                                break 'rings;
                            }
                        }
                    }
                }
                if !held {
                    failed.push(g);
                    continue;
                }
            }
            plan.surfaces[g] = Some(surface);
        }
        if failed.is_empty() {
            Ok(Ok(plan))
        } else {
            Ok(Err(Replan::Facet(failed)))
        }
    }

    /// Whether a face round its axis has a seam clear of its holes, turning
    /// a rim that is one full circle so its vertex stands in the widest
    /// stretch the holes leave free where the rims' own vertices offer none.
    fn seat_seam(&self, plan: &mut Plan, curved: &Curved, rings: &[Vec<Half>]) -> bool {
        use ogeom_geom::Curve3d as _;
        let shape = &curved.shape;
        let windings: Vec<i32> = rings
            .iter()
            .map(|ring| winding(shape, ring, self.triangles, self.points, self.tol).unwrap_or(0))
            .collect();
        let rims: Vec<usize> = (0..rings.len())
            .filter(|&k| windings[k].abs() == 1)
            .collect();
        let [low, high] = rims[..] else {
            return false;
        };
        let holes: Vec<&[Half]> = (0..rings.len())
            .filter(|&k| windings[k] == 0)
            .map(|k| rings[k].as_slice())
            .collect();
        let hole_rings = hole_polygons(shape, &holes, self.triangles, self.points, self.tol);
        let entries = |plan: &Plan, ring: &[Half]| -> Vec<(usize, bool)> {
            let mut out: Vec<(usize, bool)> = Vec::new();
            for &h in ring {
                let entry = self.entry(plan, h);
                if out.last() != Some(&entry) {
                    out.push(entry);
                }
            }
            if out.len() > 1 && out.first() == out.last() {
                out.pop();
            }
            out
        };
        let starts = |plan: &Plan, ring: &[Half]| -> Vec<Point> {
            entries(plan, ring)
                .into_iter()
                .map(
                    |(edge, forward)| match plan.edges[edge].ends[usize::from(!forward)] {
                        Corner::Mesh(v) => self.points[v as usize],
                        Corner::Placed(k) => plan.placed[k],
                    },
                )
                .collect()
        };
        let clear = |plan: &Plan| {
            choose_seam(
                shape,
                &starts(plan, &rings[low]),
                &starts(plan, &rings[high]),
                &hole_rings,
                self.tol,
            )
            .is_some()
        };
        if clear(plan) {
            return true;
        }
        let Some(angle) = free_angle(&hole_rings) else {
            return false;
        };
        for rim in [low, high] {
            let list = entries(plan, &rings[rim]);
            let [(edge, _)] = list[..] else {
                continue;
            };
            let spec = &plan.edges[edge];
            let (Curve::Circle(c), Corner::Placed(k)) = (&spec.curve, spec.ends[0]) else {
                continue;
            };
            let circle = c.circle();
            let Some((_, v)) = chart(shape, plan.placed[k], self.tol) else {
                continue;
            };
            let target = evaluate(shape, (angle, v));
            let Ok(x) = Direction::new(target - circle.centre(), self.tol) else {
                continue;
            };
            let Ok(frame) = Frame::new(circle.centre(), circle.frame().z(), x, self.tol) else {
                continue;
            };
            let Ok(turned) = ogeom_math::Circle::new(frame, circle.radius(), self.tol) else {
                continue;
            };
            let curve: Curve = ogeom_geom::CircleCurve::new(turned).into();
            let Ok(at) = curve.point_at(0.0, self.tol) else {
                continue;
            };
            plan.edges[edge].curve = curve;
            plan.placed[k] = at;
        }
        clear(plan)
    }

    /// A half-edge's planned edge, and whether the half-edge runs it
    /// forward.
    fn entry(&self, plan: &Plan, h: Half) -> (usize, bool) {
        let (a, b) = from_to(self.triangles, h);
        let (edge, along) = plan.edge_of[&(a.min(b), a.max(b))];
        (edge, along == (a < b))
    }

    /// The exact curve a chain between two faces, one of them curved,
    /// lies on: a parallel circle or a ruling of a curved face, placed on
    /// that face itself so its pcurve there is exact. Failing that, a
    /// section solved onto both, or a chord. A free chain (one curved
    /// face, nothing across) is held to its own face's parallels and
    /// rulings, and failing them is fitted onto that face as a section is.
    /// `None` when no curve holds the chain and the faces.
    fn snapped(&self, chain: &[u32], faces: &[usize]) -> Option<(Snapped, bool, Images)> {
        let read = self.snap_read(faces);
        let key = (chain.to_vec(), faces.to_vec());
        if let Some(known) = self.snaps_held().get(&key)
            && let Some((_, found)) = known.iter().find(|(was, _)| *was == read)
        {
            return found.clone();
        }
        let found = self.snap(chain, faces);
        self.snaps_held()
            .entry(key)
            .or_default()
            .push((read, found.clone()));
        found
    }

    /// What a snapped curve between `faces` reads of them.
    fn snap_read(&self, faces: &[usize]) -> Vec<SnapFace> {
        faces
            .iter()
            .map(|&g| SnapFace::of(&self.groups.carriers[g]))
            .collect()
    }

    /// The curves snapped so far.
    fn snaps_held(&self) -> std::sync::MutexGuard<'_, SnapCache> {
        self.snaps
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Whether the curve for `chain` between `faces` as they stand is known.
    fn snap_known(&self, chain: &[u32], faces: &[usize]) -> bool {
        let read = self.snap_read(faces);
        self.snaps_held()
            .get(&(chain.to_vec(), faces.to_vec()))
            .is_some_and(|known| known.iter().any(|(was, _)| *was == read))
    }

    /// The curves for chains between faces, each solved on its own, kept.
    fn snap_all(&self, wanted: Vec<(Vec<u32>, Vec<usize>)>) {
        let found =
            ogeom_core::parallel::map_ordered(&wanted, |_, (chain, faces)| self.snap(chain, faces));
        for ((chain, faces), found) in wanted.into_iter().zip(found) {
            let read = self.snap_read(&faces);
            self.snaps_held()
                .entry((chain, faces))
                .or_default()
                .push((read, found));
        }
    }

    /// As [`Self::snapped`], solved afresh.
    fn snap(&self, chain: &[u32], faces: &[usize]) -> Option<(Snapped, bool, Images)> {
        if faces.len() > 2 {
            return None;
        }
        let closed = chain.len() > 2 && chain[0] == chain[chain.len() - 1];
        let pts: Vec<Point> = chain[..chain.len() - usize::from(closed)]
            .iter()
            .map(|&v| self.points[v as usize])
            .collect();
        let reach = self.flat * REACH;
        // A band round its axis bounds itself with its own rims, so its
        // seam meets their vertices; its candidates go first.
        let mut order: Vec<usize> = faces.to_vec();
        order.sort_by_key(|&g| !self.curved(g).is_some_and(|c| c.wraps || c.wraps_v));
        // A recognized plane across the chain is where the chain lies,
        // exactly: better than a plane fitted to the chain's own points.
        let across = faces.iter().find_map(|&g| match &self.groups.carriers[g] {
            Carrier::Plane(plane) => Some(*plane),
            _ => None,
        });
        for &g in &order {
            let Some(curved) = self.curved(g) else {
                continue;
            };
            for candidate in candidates(&curved.shape, &pts, across, reach, self.tol) {
                if let Some((snapped, forward)) = self.fitted(candidate, &pts, closed, faces, reach)
                {
                    return Some((snapped, forward, Vec::new()));
                }
            }
        }
        self.section(&pts, closed, faces, reach)
            .or_else(|| self.chord(&pts, closed, faces, reach))
    }

    /// The last resort between two faces that meet all but tangentially (a
    /// fillet running on into a corner ball, or into a patch the mesh
    /// leaves faceted): no curve the two surfaces share can be solved for
    /// along the chain. The chain's own vertices lie on both, and a curve
    /// is threaded through them, a line for a single span; its tolerance is
    /// how far it strays from either surface between them, up to a
    /// twentieth of its longest span. A face bounded by such a curve is
    /// good to its tolerance, where it would otherwise fall to facets
    /// whole.
    fn chord(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool, Images)> {
        use ogeom_geom::Curve3d as _;
        let [a, b] = faces[..] else {
            return None;
        };
        let (fa, fb) = (self.signed(a)?, self.signed(b)?);
        let mut on: Vec<Point> = if closed {
            pts.to_vec()
        } else {
            // Along a tangency the chain wanders, and can step past its end
            // and back; a curve through such a step hooks back on itself.
            // Only the points that close in on the end are threaded.
            let end = pts[pts.len() - 1];
            let mut kept = vec![pts[0]];
            for &p in &pts[1..pts.len() - 1] {
                let last = kept[kept.len() - 1];
                if p.distance(end) < last.distance(end) && (p - last).dot(end - last) > 0.0 {
                    kept.push(p);
                }
            }
            kept.push(end);
            kept
        };
        if closed {
            on.push(pts[0]);
        }
        let longest = on
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0_f64, f64::max);
        if longest <= self.tol.confusion() {
            return None;
        }
        // A cubic through the points first; where it overshoots between
        // them (a short span beside long ones), the polyline through them.
        for degree in [3, 1] {
            let (curve, samples): (Curve, Vec<f64>) = if on.len() == 2 {
                let length = on[0].distance(on[1]);
                let line: Curve = LineCurve::segment(on[0], on[1], self.tol).ok()?.into();
                (
                    line,
                    (0..=16).map(|k| length * f64::from(k) / 16.0).collect(),
                )
            } else {
                // A loop is carried on past its ends and cut back, as a
                // fitted section is.
                let n = on.len();
                let pad = if closed { 3.min(n / 3) } else { 0 };
                let padded: Vec<Point> = if pad > 0 {
                    on[n - 1 - pad..n - 1]
                        .iter()
                        .chain(&on)
                        .chain(&on[1..=pad])
                        .copied()
                        .collect()
                } else {
                    on.clone()
                };
                let parameters =
                    crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
                let mut spline =
                    crate::fit::interpolate_at(&padded, &parameters, degree, self.tol).ok()?;
                if pad > 0 {
                    spline = spline.split_at(parameters[pad], self.tol).ok()?.1;
                    spline = spline.split_at(parameters[pad + n - 1], self.tol).ok()?.0;
                }
                let own = &parameters[pad..pad + n];
                let samples = own
                    .windows(2)
                    .flat_map(|w| {
                        [
                            w[0],
                            w[0] + (w[1] - w[0]) * 0.25,
                            w[0] + (w[1] - w[0]) * 0.5,
                            w[0] + (w[1] - w[0]) * 0.75,
                        ]
                    })
                    .chain(std::iter::once(own[n - 1]))
                    .collect();
                (spline.into(), samples)
            };
            let range = curve.domain();
            let mut tolerance = self.tol.confusion();
            // The curve runs the way the points it was threaded through do,
            // or it hooks back on itself between them.
            let mut runs = true;
            for (i, &t) in samples.iter().enumerate() {
                let t = t.clamp(range.0, range.1);
                let p = curve.point_at(t, self.tol).ok()?;
                tolerance = tolerance.max(fa(p).abs()).max(fb(p).abs());
                let span = (i / 4).min(on.len() - 2);
                runs &= curve.d1_at(t, self.tol).ok()?.dot(on[span + 1] - on[span]) > 0.0;
            }
            if tolerance > reach.max(longest * CHORD_SAG) {
                return None;
            }
            if !runs {
                continue;
            }
            return Some((
                if closed {
                    Snapped::Loop(curve, range, tolerance)
                } else {
                    Snapped::Open(curve, range, tolerance)
                },
                true,
                Vec::new(),
            ));
        }
        None
    }

    /// The curve two faces meet along, where it is no parallel or ruling of
    /// either: vertices of the chain, and points between them, solved onto
    /// both surfaces, and the curve interpolated through them. Its image on
    /// each curved face is interpolated through the same points' chart
    /// positions at the same parameters, so the two run together. `None`
    /// where the surfaces meet tangentially, the solve does not settle, or
    /// the curve strays from either face past the reach.
    ///
    /// A free chain, with one curved face and nothing across, is fitted
    /// the same way through the feet of its vertices on that face's
    /// surface; its tolerance also holds every vertex of the chain, so the
    /// edge keeps to the mesh's boundary within it.
    fn section(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool, Images)> {
        // Made through more points until it keeps to the surfaces to a
        // fiftieth of a micron's worth of confusion distances, as close as
        // an exact edge's image would, or the points run out.
        let close = self.tol.confusion() * 50.0;
        // A free chain's points are its own vertices, however many are
        // asked for.
        if faces.len() == 1 {
            return self
                .section_through(pts, closed, faces, reach, SECTION_POINTS)
                .map(|(_, found)| found);
        }
        let mut best: Option<(f64, (Snapped, bool, Images))> = None;
        let mut count = SECTION_POINTS;
        while count <= SECTION_POINTS * 8 {
            let Some((worst, found)) = self.section_through(pts, closed, faces, reach, count)
            else {
                break;
            };
            let done = worst <= close;
            if best.as_ref().is_none_or(|(held, _)| worst < *held) {
                best = Some((worst, found));
            }
            if done {
                break;
            }
            count *= 2;
        }
        best.map(|(_, found)| found)
    }

    /// One fitted section through at least `count` points, with the worst
    /// of its own stray and its images'.
    fn section_through(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
        count: usize,
    ) -> Option<(f64, (Snapped, bool, Images))> {
        use ogeom_geom::Curve3d as _;
        let (fa, fb) = match faces[..] {
            [a, b] => (self.signed(a)?, Some(self.signed(b)?)),
            [a] => (self.signed(a)?, None),
            _ => return None,
        };
        // A point on the faces near a start: solved onto both, or for a
        // free chain the foot on its face.
        let onto = |start: Point, limit: f64| -> Option<Point> {
            match &fb {
                Some(fb) => onto_both(&fa, fb, start, reach, limit),
                None => {
                    let shape = &self.curved(faces[0])?.shape;
                    let foot = evaluate(shape, chart(shape, start, self.tol)?);
                    (fa(foot).abs() <= reach * 1e-3 && foot.distance(start) <= limit)
                        .then_some(foot)
                }
            }
        };
        let off = |p: Point| fa(p).abs().max(fb.as_ref().map_or(0.0, |fb| fb(p).abs()));
        let steps = if closed { pts.len() } else { pts.len() - 1 };
        // A long chain is taken a few vertices at a time: the curve needs
        // its shape, not every vertex. A free chain is taken through every
        // vertex and nothing between: the foot of a point on a chord of a
        // curved boundary lies off the boundary, across the surface, by as
        // much as the chord cuts off it. A free chain of one or two spans
        // has too few vertices for a cubic, and its spans are split.
        let free = fb.is_none();
        let stride = if free {
            1
        } else {
            steps.div_ceil(SECTION_SPANS).max(1)
        };
        let split = if free {
            if steps >= 3 { 1 } else { SECTION_SPLIT }
        } else {
            count.div_ceil(steps.div_ceil(stride)).max(SECTION_SPLIT)
        };
        // A chord of the mesh lies across a face recognized from it, and a
        // face recognized from a mesh turns through well under a right
        // angle from one end of one of its chords to the other. Where the
        // fitted surface turns farther, it is not the surface the chord
        // lies on, and the curve solved on it runs where the mesh does not.
        let turns_away = |p: Point, q: Point| {
            faces.iter().filter_map(|&g| self.curved(g)).any(|curved| {
                let (gp, gq) = (gradient(&curved.shape, p), gradient(&curved.shape, q));
                gp.dot(gq) <= 0.0
            })
        };
        let mut on = Vec::new();
        let mut i = 0;
        while i < steps {
            let next = (i + stride).min(steps);
            let (p, q) = (pts[i], pts[next % pts.len()]);
            if turns_away(p, q) {
                return None;
            }
            for k in 0..split {
                #[allow(clippy::cast_precision_loss, reason = "a handful of splits")]
                let f = k as f64 / split as f64;
                let limit = if k == 0 { reach } else { p.distance(q) + reach };
                on.push(onto(p + (q - p) * f, limit)?);
            }
            i = next;
        }
        on.push(if closed {
            on[0]
        } else {
            onto(pts[pts.len() - 1], reach)?
        });
        // A loop is interpolated with a few of its points carried on past
        // each end, then cut back to its own: an open interpolation left
        // free at its ends wanders where the loop meets itself.
        let pad = if closed {
            SECTION_PAD.min(on.len() / 4)
        } else {
            0
        };
        let n = on.len();
        let padded: Vec<Point> = if pad > 0 {
            on[n - 1 - pad..n - 1]
                .iter()
                .chain(&on)
                .chain(&on[1..=pad])
                .copied()
                .collect()
        } else {
            on.clone()
        };
        let parameters =
            crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
        let cut = (parameters[pad], parameters[pad + n - 1]);
        let trimmed = |spline: ogeom_geom::BSplineCurve| -> Option<ogeom_geom::BSplineCurve> {
            if pad == 0 {
                return Some(spline);
            }
            let (_, after) = spline.split_at(cut.0, self.tol).ok()?;
            Some(after.split_at(cut.1, self.tol).ok()?.0)
        };
        let curve: Curve =
            trimmed(crate::fit::interpolate_at(&padded, &parameters, 3, self.tol).ok()?)?.into();
        let range = curve.domain();
        let parameters = parameters[pad..pad + n].to_vec();
        // Between the points it was made through, how far it strays from
        // either surface.
        let between = |k: usize, f: f64| parameters[k] + (parameters[k + 1] - parameters[k]) * f;
        // Its ends are the chain's ends solved onto both surfaces, a hair
        // from the vertices they meet.
        let last = if closed { pts[0] } else { pts[pts.len() - 1] };
        let mut tolerance = self
            .tol
            .confusion()
            .max(on[0].distance(pts[0]))
            .max(on[on.len() - 1].distance(last));
        for k in 0..on.len() - 1 {
            for f in [0.25, 0.5, 0.75] {
                let p = curve.point_at(between(k, f), self.tol).ok()?;
                tolerance = tolerance.max(off(p));
            }
        }
        if free {
            // Every vertex of a free chain to the foot the curve runs
            // through.
            for (i, p) in pts.iter().enumerate() {
                tolerance = tolerance.max(p.distance(on[i * split]));
            }
        }
        if tolerance > reach {
            return None;
        }
        let mut images = Vec::new();
        for &g in faces {
            let Some(curved) = self.curved(g) else {
                continue;
            };
            let (pu, pv) = periodic(&curved.shape);
            let mut uv: Vec<Point> = Vec::with_capacity(padded.len());
            for p in &padded {
                let (u, v) = match uv.last() {
                    None => unwrapped(curved, *p, self.tol)?,
                    Some(last) => {
                        let (u, v) = chart(&curved.shape, *p, self.tol)?;
                        let near = |x: f64, c: f64, wraps: bool| {
                            if wraps {
                                c + ogeom_math::elementary::wrap_signed_angle(x - c)
                            } else {
                                x
                            }
                        };
                        (near(u, last.x, pu), near(v, last.y, pv))
                    }
                };
                uv.push(Point::new(u, v, 0.0));
            }
            let all =
                crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
            let flat = trimmed(crate::fit::interpolate_at(&uv, &all, 3, self.tol).ok()?)?;
            let control: Vec<Point2> = flat
                .control_points()
                .iter()
                .map(|c| Point2::new(c.scaled.x, c.scaled.y))
                .collect();
            let pcurve: PlanarCurve =
                ogeom_geom::BSpline2d::new(flat.knots().clone(), control, self.tol)
                    .ok()?
                    .into();
            let mut deviation = self.tol.confusion();
            for k in 0..on.len() - 1 {
                for f in [0.0, 0.25, 0.5, 0.75] {
                    use ogeom_geom::Curve2d as _;
                    let t = between(k, f);
                    let at = pcurve.point_at(t, self.tol).ok()?;
                    let lifted = evaluate(&curved.shape, (at.x, at.y));
                    deviation = deviation.max(lifted.distance(curve.point_at(t, self.tol).ok()?));
                }
            }
            if deviation > reach {
                return None;
            }
            images.push((g, pcurve, deviation));
        }
        let worst = images.iter().map(|(_, _, d)| *d).fold(tolerance, f64::max);
        Some((
            worst,
            (
                if closed {
                    Snapped::Loop(curve, range, tolerance)
                } else {
                    Snapped::Open(curve, range, tolerance)
                },
                true,
                images,
            ),
        ))
    }

    /// A face's surface as a signed distance, where it has one.
    fn signed(&self, g: usize) -> Option<Box<dyn Fn(Point) -> f64>> {
        match &self.groups.carriers[g] {
            Carrier::Plane(plane) => {
                let plane = *plane;
                Some(Box::new(move |p: Point| plane.signed_distance_to(p)))
            }
            Carrier::Curved(c) => {
                let shape = c.shape.clone();
                Some(Box::new(move |p: Point| shape.signed_distance_to(p)))
            }
            Carrier::Gone => None,
        }
    }

    /// A candidate curve held against the chain and the faces: its range,
    /// which way it runs along the chain, and its tolerance.
    fn fitted(
        &self,
        curve: Curve,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool)> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let parameter = |p: Point| -> Option<f64> {
            match &curve {
                Curve::Line(l) => {
                    let axis = l.axis();
                    Some((p - axis.location).dot(axis.direction.vector()))
                }
                Curve::Circle(c) => {
                    ogeom_math::elementary::circle_parameter(&c.circle(), p, self.tol).ok()
                }
                _ => None,
            }
        };
        let mut tolerance = self.tol.confusion();
        let mut ts = Vec::with_capacity(pts.len());
        for p in pts {
            let t = parameter(*p)?;
            let on = curve.point_at(t, self.tol).ok()?;
            tolerance = tolerance.max(on.distance(*p));
            ts.push(t);
        }
        // The sweep along the chain, unwrapped for a circle.
        let mut sweep = 0.0;
        let steps = if closed { ts.len() } else { ts.len() - 1 };
        for i in 0..steps {
            let (a, b) = (ts[i], ts[(i + 1) % ts.len()]);
            let d = b - a;
            sweep += if matches!(curve, Curve::Circle(_)) {
                ogeom_math::elementary::wrap_signed_angle(d)
            } else {
                d
            };
        }
        let (snapped, forward) = if closed {
            if !matches!(curve, Curve::Circle(_)) || (sweep.abs() - tau).abs() > 1e-3 {
                return None;
            }
            (None, sweep > 0.0)
        } else if sweep > 0.0 {
            (Some((ts[0], ts[0] + sweep)), true)
        } else {
            let last = ts[ts.len() - 1];
            (Some((last, last - sweep)), false)
        };
        let range = snapped.unwrap_or((0.0, tau));
        if range.1 - range.0 <= self.tol.parametric() {
            return None;
        }
        // Drawn in every plane it bounds: a circle is imaged in a plane by
        // projection, which a circle standing across the plane has not.
        for &g in faces {
            if let Carrier::Plane(plane) = &self.groups.carriers[g] {
                let surface: ogeom_geom::SurfaceGeometry = PlaneSurface::new(*plane).into();
                ogeom_intersect::exact_pcurve_of(&curve, &surface, self.tol)?;
            }
        }
        // Standing on every face it bounds.
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            let p = curve.point_at(t, self.tol).ok()?;
            for &g in faces {
                let off = match &self.groups.carriers[g] {
                    Carrier::Plane(plane) => plane.signed_distance_to(p).abs(),
                    Carrier::Curved(c) => c.shape.distance_to(p),
                    Carrier::Gone => return None,
                };
                tolerance = tolerance.max(off);
            }
        }
        if tolerance > reach {
            return None;
        }
        Some((
            match snapped {
                Some(range) => Snapped::Open(curve, range, tolerance),
                None => Snapped::Closed(curve, tolerance),
            },
            forward,
        ))
    }
}

#[derive(Clone)]
enum Snapped {
    Open(Curve, (f64, f64), f64),
    Closed(Curve, f64),
    /// A closed curve that is no circle, starting and ending at the
    /// chain's first vertex.
    Loop(Curve, (f64, f64), f64),
}

/// How many pieces each span of a chain is cut into for a fitted section:
/// its vertices alone leave a coarse chain's curve unconstrained between
/// them.
const SECTION_SPLIT: usize = 4;

/// The most spans a chain is taken in for a fitted section.
const SECTION_SPANS: usize = 40;

/// How many points a closed section is carried on past each end while it
/// is interpolated.
const SECTION_PAD: usize = 8;

/// How many points a fitted section is interpolated through, at least:
/// a short chain's spans are cut finer to reach it.
const SECTION_POINTS: usize = 160;

/// A fitted section's images on the curved faces it bounds: the face, the
/// image, and how far the image strays from the curve.
type Images = Vec<(usize, PlanarCurve, f64)>;

/// What a snapped curve reads of one face beside its chain: the plane, or
/// the curved surface with the chart point its images unwrap around and
/// whether it wraps.
#[derive(PartialEq)]
enum SnapFace {
    Plane(Plane),
    Curved(Canonical, (f64, f64), bool, bool),
    Gone,
}

impl SnapFace {
    fn of(carrier: &Carrier) -> Self {
        match carrier {
            Carrier::Plane(plane) => Self::Plane(*plane),
            Carrier::Curved(c) => Self::Curved(c.shape.clone(), c.centre, c.wraps, c.wraps_v),
            Carrier::Gone => Self::Gone,
        }
    }
}

/// The curves snapped to chains, kept across the plans of one conversion:
/// by chain and faces, each with what it read of the faces. The coplanar
/// distance and the points are those of the whole conversion.
type SnapCache =
    HashMap<(Vec<u32>, Vec<usize>), Vec<(Vec<SnapFace>, Option<(Snapped, bool, Images)>)>>;

/// A point solved onto where two surfaces meet, from a start near both:
/// Newton's step, the least one that zeroes both signed distances to first
/// order. `None` where the surfaces meet tangentially there, where the
/// solve does not settle within a thousandth of `reach` of both, or where
/// it lands farther than `limit` from its start.
fn onto_both(
    fa: &dyn Fn(Point) -> f64,
    fb: &dyn Fn(Point) -> f64,
    start: Point,
    reach: f64,
    limit: f64,
) -> Option<Point> {
    let mut p = start;
    let h = 1e-7 * (1.0 + p.to_vector().magnitude());
    let gradient = |f: &dyn Fn(Point) -> f64, p: Point| {
        let d = |v: Vector| (f(p + v * h) - f(p - v * h)) / (2.0 * h);
        Vector::new(d(Vector::X), d(Vector::Y), d(Vector::Z))
    };
    let (mut last, mut stalled) = (f64::INFINITY, 0);
    for _ in 0..40 {
        let (va, vb) = (fa(p), fb(p));
        let residual = va.abs().max(vb.abs());
        if residual <= 1e-13 * (1.0 + p.to_vector().magnitude()) {
            break;
        }
        // A surface whose distance is itself a search answers to its
        // search's precision, short of the residual asked: two steps that
        // do not halve it have reached that floor.
        stalled = if residual > last * 0.5 {
            stalled + 1
        } else {
            0
        };
        if stalled >= 2 {
            break;
        }
        last = residual;
        let (ga, gb) = (gradient(fa, p), gradient(fb, p));
        let (aa, ab, bb) = (ga.dot(ga), ga.dot(gb), gb.dot(gb));
        let det = aa.mul_add(bb, -(ab * ab));
        if det <= 1e-12 * aa * bb {
            return None;
        }
        let la = (va * bb - vb * ab) / det;
        let lb = (vb * aa - va * ab) / det;
        p = p - ga * la - gb * lb;
    }
    (fa(p).abs().max(fb(p).abs()) <= reach * 1e-3 && p.distance(start) <= limit).then_some(p)
}

/// The curves a chain on a curved surface may be: the parallel circle
/// through its mean height, when every point sits at one height and one
/// distance from the axis; and, on a cylinder or a cone, the ruling
/// through its mean angle, when the chain is straight along it. Each is
/// placed on the surface exactly, its angle measured from the surface's
/// own origin.
fn candidates(
    shape: &Canonical,
    pts: &[Point],
    across: Option<Plane>,
    reach: f64,
    tol: Tolerances,
) -> Vec<Curve> {
    let mut out = Vec::new();
    let plane_of = |pts: &[Point]| match across {
        Some(plane) => Some((plane.project(pts[0]), plane.normal())),
        None => plane_through(pts, tol),
    };
    if let Canonical::Sphere(sphere) = shape {
        // The section of the sphere by the chain's plane. Where that plane
        // is square to the sphere's axis the section is a latitude, and is
        // built on the sphere's own frame, so its parameter is the sphere's
        // angle and it starts where the sphere's seam does.
        let axis = sphere.frame().z();
        if pts.len() >= 3
            && let Some((centre, normal)) = plane_of(pts)
            && normal.vector().cross(axis.vector()).magnitude()
                <= if across.is_some() { 1e-12 } else { 1e-3 }
        {
            let h = (centre - sphere.centre()).dot(axis.vector());
            let r2 = sphere.radius().powi(2) - h * h;
            if r2 > 0.0
                && let Ok(frame) = Frame::new(
                    sphere.centre() + axis.vector() * h,
                    axis,
                    sphere.frame().x(),
                    tol,
                )
                && let Ok(circle) = ogeom_math::Circle::new(frame, r2.sqrt(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
                return out;
            }
        }
        if pts.len() >= 3
            && let Some((centre, normal)) = plane_of(pts)
        {
            let n = normal.vector();
            let d = (centre - sphere.centre()).dot(n);
            let r2 = sphere.radius().powi(2) - d * d;
            if r2 > 0.0
                && let Ok(frame) = Frame::new(
                    sphere.centre() + n * d,
                    normal,
                    normal.any_perpendicular(),
                    tol,
                )
                && let Ok(circle) = ogeom_math::Circle::new(frame, r2.sqrt(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
            }
        }
        return out;
    }
    let Some(frame) = axis_frame(shape) else {
        return out;
    };
    let (o, z) = (frame.origin(), frame.z().vector());
    let heights: Vec<f64> = pts.iter().map(|p| (*p - o).dot(z)).collect();
    let radii: Vec<f64> = pts
        .iter()
        .zip(&heights)
        .map(|(p, h)| ((*p - o) - z * *h).magnitude())
        .collect();
    #[allow(
        clippy::cast_precision_loss,
        reason = "chain lengths are far below 2^52"
    )]
    let count = pts.len() as f64;
    let mean_h = heights.iter().sum::<f64>() / count;
    let mean_r = radii.iter().sum::<f64>() / count;
    let level = heights.iter().all(|h| (h - mean_h).abs() <= reach)
        && radii.iter().all(|r| (r - mean_r).abs() <= reach);
    if level {
        let radius = match shape {
            Canonical::Cylinder(c) => c.radius(),
            Canonical::Cone(c) => c.radius_at(mean_h),
            Canonical::Torus(t) => {
                let (big, small) = (t.major_radius(), t.minor_radius());
                let off = (small * small - mean_h * mean_h).max(0.0).sqrt();
                // At the tube's crest the circle's radius turns infinitely
                // fast with its height, and a few slops of height move it by
                // many: a ring within the reach of the crest is the crest.
                if (small - mean_h.abs()).abs() <= reach {
                    big
                } else if (big + off - mean_r).abs() <= (big - off - mean_r).abs() {
                    big + off
                } else {
                    big - off
                }
            }
            _ => mean_r,
        };
        if let Ok(at) = Frame::new(o + z * mean_h, frame.z(), frame.x(), tol)
            && let Ok(circle) = ogeom_math::Circle::new(at, radius, tol)
        {
            out.push(ogeom_geom::CircleCurve::new(circle).into());
        }
    }
    // A circle of a torus's tube: every point in one plane through the
    // axis, at the angle the chain stands at. Built with the tube's own
    // angle as its parameter, starting at the outer equator.
    if let Canonical::Torus(torus) = shape
        && pts.len() >= 2
    {
        let mut angles: Vec<f64> = pts
            .iter()
            .filter_map(|x| chart(shape, *x, tol).map(|c| c.0))
            .collect();
        if !angles.is_empty() {
            let (u, _) = angular_spread(&mut angles);
            let (x, y) = (frame.x().vector(), frame.y().vector());
            let out_u = x * u.cos() + y * u.sin();
            let across = y * u.cos() - x * u.sin();
            let in_plane = pts.iter().all(|p| (*p - o).dot(across).abs() <= reach);
            if in_plane
                && let Ok(radial) = Direction::new(out_u, tol)
                && let Ok(normal) = Direction::new(out_u.cross(z), tol)
                && let Ok(at) = Frame::new(o + out_u * torus.major_radius(), normal, radial, tol)
                && let Ok(circle) = ogeom_math::Circle::new(at, torus.minor_radius(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
            }
        }
    }
    if matches!(shape, Canonical::Cylinder(_) | Canonical::Cone(_)) && pts.len() >= 2 {
        let (p, q) = (pts[0], pts[pts.len() - 1]);
        let straight = pts.iter().all(|x| distance_to_line(*x, p, q) <= reach);
        let mut angles: Vec<f64> = pts
            .iter()
            .filter_map(|x| chart(shape, *x, tol).map(|c| c.0))
            .collect();
        if straight && !angles.is_empty() {
            let (u, _) = angular_spread(&mut angles);
            let (lo, hi) = (
                heights.iter().copied().fold(f64::INFINITY, f64::min),
                heights.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            );
            let (a, b) = (evaluate(shape, (u, lo)), evaluate(shape, (u, hi)));
            if let Ok(line) = LineCurve::segment(a, b, tol) {
                // Measured from the chain's first end, so its parameter is
                // the chain's distance along it.
                let _ = line;
                let (a, b) = if (pts[0] - a).magnitude() <= (pts[0] - b).magnitude() {
                    (a, b)
                } else {
                    (b, a)
                };
                if let Ok(line) = LineCurve::segment(a, b, tol) {
                    out.push(line.into());
                }
            }
        }
    }
    out
}

/// The plane nearest a set of points: their centroid and the covariance's
/// smallest direction.
fn plane_through(points: &[Point], tol: Tolerances) -> Option<(Point, Direction)> {
    #[allow(
        clippy::cast_precision_loss,
        reason = "chain lengths are far below 2^52"
    )]
    let count = points.len() as f64;
    let c = points.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / count;
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    for p in points {
        let d = p.to_vector() - c;
        let v = nalgebra::Vector3::new(d.x, d.y, d.z);
        m += v * v.transpose();
    }
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut best = 0;
    for i in 1..3 {
        if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
            best = i;
        }
    }
    let v = eigen.eigenvectors.column(best);
    Some((
        Point::from_vector(c),
        Direction::new(Vector::new(v[0], v[1], v[2]), tol).ok()?,
    ))
}

/// The surface a curved face is built on, windowed along its axis to hold
/// every point its boundary reaches.
fn surface_of(
    curved: &Curved,
    points: &[Point],
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::SurfaceGeometry> {
    use ogeom_geom::{ConeSurface, CylinderSurface, SphereSurface, TorusSurface};
    let heights = || {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &curved.vertices {
            if let Some((_, h)) = chart(&curved.shape, points[v as usize], tol) {
                lo = lo.min(h);
                hi = hi.max(h);
            }
        }
        let margin = (hi - lo).mul_add(0.25, tol.confusion() * 10.0);
        (lo - margin, hi + margin)
    };
    Ok(match &curved.shape {
        Canonical::Cylinder(c) => CylinderSurface::new(*c, heights())?.into(),
        Canonical::Cone(c) => {
            // Short of the apex, where the cone's radius runs out.
            let (lo, hi) = heights();
            let apex = -c.reference_radius() / c.half_angle().tan();
            let lo = lo.max(apex + (hi - apex) * 1e-6);
            ConeSurface::new(*c, (lo, hi))?.into()
        }
        Canonical::Sphere(s) => SphereSurface::new(*s).into(),
        Canonical::Torus(t) => TorusSurface::new(*t).into(),
        Canonical::Plane(p) => PlaneSurface::new(*p).into(),
        Canonical::Swept(s) => s.surface.clone(),
    })
}

/// An edge's image on a curved face, and how far it strays from the edge:
/// the straight chart segment where the edge is a parallel or a ruling of
/// the face, otherwise a fit by projection; `None` past the reach.
fn image_on(
    curved: &Curved,
    surface: &ogeom_geom::SurfaceGeometry,
    curve: &Curve,
    range: (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    if let Some((pcurve, deviation)) = straight_image(curved, curve, range, tol)
        && deviation <= reach
    {
        return Some((pcurve, deviation));
    }
    if let Some((pcurve, deviation)) = interpolated_image(curved, curve, range, tol)
        && deviation <= reach
    {
        return Some((pcurve, deviation));
    }
    let (pcurve, error, _, off, _) =
        crate::pcurve_fit::fit_projected_pcurve_capped(curve, range, surface, reach, tol).ok()?;
    let deviation = error.max(off);
    (deviation <= reach).then_some((pcurve, deviation))
}

/// An edge's image interpolated through the chart positions of points
/// along it, at the edge's own parameters: the curve and its image agree
/// at every sample by construction, and between them to the fourth power
/// of the spacing, which is doubled until the image keeps to the curve
/// within a hundredth of the confusion distance. `None` where a sample
/// has no chart position (a pole).
fn interpolated_image(
    curved: &Curved,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    use ogeom_geom::{Curve2d as _, Curve3d as _};
    let (pu, pv) = periodic(&curved.shape);
    let near = |x: f64, c: f64, wraps: bool| {
        if wraps {
            c + ogeom_math::elementary::wrap_signed_angle(x - c)
        } else {
            x
        }
    };
    let mut best: Option<(PlanarCurve, f64)> = None;
    let mut count: u32 = 32;
    while count <= 1024 {
        let parameters: Vec<f64> = (0..=count)
            .map(|k| range.0 + (range.1 - range.0) * f64::from(k) / f64::from(count))
            .collect();
        let mut uv: Vec<Point> = Vec::with_capacity(parameters.len());
        for &t in &parameters {
            let p = curve.point_at(t, tol).ok()?;
            let (u, v) = match uv.last() {
                None => unwrapped(curved, p, tol)?,
                Some(last) => {
                    let (u, v) = chart(&curved.shape, p, tol)?;
                    (near(u, last.x, pu), near(v, last.y, pv))
                }
            };
            uv.push(Point::new(u, v, 0.0));
        }
        let flat = crate::fit::interpolate_at(&uv, &parameters, 3, tol).ok()?;
        let control: Vec<Point2> = flat
            .control_points()
            .iter()
            .map(|c| Point2::new(c.scaled.x, c.scaled.y))
            .collect();
        let pcurve: PlanarCurve = ogeom_geom::BSpline2d::new(flat.knots().clone(), control, tol)
            .ok()?
            .into();
        let mut deviation = tol.confusion() * 1e-2;
        for pair in parameters.windows(2) {
            for f in [0.25, 0.5, 0.75] {
                let t = pair[0] + (pair[1] - pair[0]) * f;
                let at = pcurve.point_at(t, tol).ok()?;
                let lifted = evaluate(&curved.shape, (at.x, at.y));
                deviation = deviation.max(lifted.distance(curve.point_at(t, tol).ok()?));
            }
        }
        let done = deviation <= tol.confusion() * 1e-2;
        if best.as_ref().is_none_or(|(_, held)| deviation < *held) {
            best = Some((pcurve, deviation));
        }
        if done {
            break;
        }
        count *= 2;
    }
    best
}

/// A degree-one pcurve over the edge's range, from the chart points of its
/// ends on the face's branch, and how far it strays along its length.
fn straight_image(
    curved: &Curved,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    use ogeom_geom::Curve3d as _;
    let at =
        |t: f64| -> Option<(f64, f64)> { unwrapped(curved, curve.point_at(t, tol).ok()?, tol) };
    let start = at(range.0)?;
    // Carried along the edge a quarter at a time, so a full turn ends a
    // whole period from where it began.
    let (pu, pv) = periodic(&curved.shape);
    let step = |a: f64, b: f64, wraps: bool| {
        if wraps {
            a + ogeom_math::elementary::wrap_signed_angle(b - a)
        } else {
            b
        }
    };
    let mut end = start;
    for k in 1..=4 {
        let next = at(range.0 + (range.1 - range.0) * f64::from(k) / 4.0)?;
        end = (step(end.0, next.0, pu), step(end.1, next.1, pv));
    }
    let pcurve = linear(start, end, range, tol).ok()?;
    let mut deviation: f64 = 0.0;
    for k in 0..=16 {
        let f = f64::from(k) / 16.0;
        let t = range.0 + (range.1 - range.0) * f;
        let uv = (
            (end.0 - start.0).mul_add(f, start.0),
            (end.1 - start.1).mul_add(f, start.1),
        );
        deviation =
            deviation.max(evaluate(&curved.shape, uv).distance(curve.point_at(t, tol).ok()?));
    }
    Some((pcurve, deviation))
}

/// Builds the planned vertices, edges and faces.
struct Builder<'a> {
    model: &'a mut Model,
    points: &'a [Point],
    triangles: &'a [[u32; 3]],
    groups: &'a Groups,
    plan: &'a Plan,
    tol: Tolerances,
}

impl Builder<'_> {
    fn entry(&self, h: Half) -> (usize, bool) {
        let (a, b) = from_to(self.triangles, h);
        let (edge, along) = self.plan.edge_of[&(a.min(b), a.max(b))];
        (edge, along == (a < b))
    }

    /// The faces, one per live group, by group index.
    fn build(mut self) -> OgeomResult<Vec<Option<Shape>>> {
        let mut corners: HashMap<Corner, Shape> = HashMap::new();
        let mut edges: Vec<Shape> = Vec::with_capacity(self.plan.edges.len());
        for spec in &self.plan.edges {
            for corner in spec.ends {
                corners.entry(corner).or_insert_with(|| {
                    let at = match corner {
                        Corner::Mesh(v) => self.points[v as usize],
                        Corner::Placed(i) => self.plan.placed[i],
                    };
                    self.model.add_vertex(VertexData::new(at))
                });
            }
            let id = self.model.geometry_mut().add_curve(spec.curve.clone());
            let mut data = EdgeData::on_curve(id, Location::identity(), spec.range);
            // A tolerance measured as a gap is held a millionth wider: the
            // checker measures the same gap again, a rounding apart.
            data.tolerance = Tolerance::new(spec.tolerance * (1.0 + TOLERANCE_MARGIN))?;
            let bounds = if spec.ends[0] == spec.ends[1] {
                vec![
                    corners[&spec.ends[0]].clone(),
                    corners[&spec.ends[0]].clone(),
                ]
            } else {
                vec![
                    corners[&spec.ends[0]].clone(),
                    corners[&spec.ends[1]].clone(),
                ]
            };
            edges.push(self.model.add_edge(data, &bounds)?);
        }

        let mut faces = Vec::with_capacity(self.groups.carriers.len());
        for (g, carrier) in self.groups.carriers.iter().enumerate() {
            let rings = &self.plan.loops[g];
            let face = match carrier {
                Carrier::Gone => None,
                Carrier::Curved(curved) if self.plan.layouts[g] == Layout::Whole => {
                    Some(self.whole_face(curved, g, 0.0)?)
                }
                _ if rings.is_empty() => None,
                Carrier::Plane(plane) => Some(self.plane_face(*plane, rings, &edges)?),
                Carrier::Curved(curved) => match self.plan.layouts[g] {
                    Layout::Band { round_tube } => {
                        Some(self.band_face(curved, rings, &edges, round_tube)?)
                    }
                    Layout::Cap => Some(self.cap_face(curved, rings, &edges)?),
                    Layout::Wrapped => Some(self.wrapped_face(curved, g, rings, &edges)?),
                    Layout::Holed => Some(self.holed_face(curved, g, rings, &edges)?),
                    Layout::Open | Layout::Whole => {
                        Some(self.curved_face(curved, g, rings, &edges)?)
                    }
                },
            };
            faces.push(face);
        }
        Ok(faces)
    }

    /// A ring's planned edges in walking order, repeats run together.
    fn entries(&self, ring: &[Half]) -> Vec<(usize, bool)> {
        let mut entries: Vec<(usize, bool)> = Vec::new();
        for &h in ring {
            let entry = self.entry(h);
            if entries.last() != Some(&entry) {
                entries.push(entry);
            }
        }
        if entries.len() > 1 && entries.first() == entries.last() {
            entries.pop();
        }
        entries
    }

    fn has_pcurve(&self, edge: &Shape, surface: ogeom_topo::SurfaceId) -> bool {
        self.model.node(edge).and_then(|n| n.data().as_edge()).is_some_and(|d| {
            d.representations.iter().any(|rep| {
                matches!(rep, EdgeRepr::PCurve { surface: s, .. } | EdgeRepr::Seam { surface: s, .. } if *s == surface)
            })
        })
    }

    fn plane_face(
        &mut self,
        plane: Plane,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let geometry: ogeom_geom::SurfaceGeometry = PlaneSurface::new(plane).into();
        let surface = self.model.geometry_mut().add_surface(geometry.clone());
        let local = |p: Point| {
            let l = plane.frame().to_local(p);
            Point2::new(l.x, l.y)
        };
        let mut wires: Vec<(f64, Shape)> = Vec::with_capacity(rings.len());
        for ring in rings {
            let area = ring_area(ring, self.triangles, |p| Some(local(p)), self.points);
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
                let spec = &self.plan.edges[edge];
                if !self.has_pcurve(&edges[edge], surface) {
                    let Some(pcurve) =
                        ogeom_intersect::exact_pcurve_of(&spec.curve, &geometry, self.tol)
                    else {
                        ogeom_bail!(Construction, "an edge has no image in its face's plane");
                    };
                    crate::build::attach_pcurve(
                        self.model,
                        &edges[edge],
                        pcurve,
                        surface,
                        Location::identity(),
                        spec.range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            wires.push((area, self.model.add_wire(&ring_edges)?));
        }
        wires.sort_by(|a, b| b.0.total_cmp(&a.0));
        let wires: Vec<Shape> = wires.into_iter().map(|(_, w)| w).collect();
        self.model
            .add_face(FaceData::new(surface, Location::identity()), &wires)
    }

    /// Whether the region's outward side is the surface's own normal side.
    fn outward(&self, curved: &Curved, g: usize) -> bool {
        let mut vote = 0.0;
        for (t, tri) in self.triangles.iter().enumerate() {
            if self.groups.of[t] != g {
                continue;
            }
            let corners = tri.map(|v| self.points[v as usize]);
            let centroid = Point::from_vector(
                (corners[0].to_vector() + corners[1].to_vector() + corners[2].to_vector()) / 3.0,
            );
            let n = unit_normal(self.points, *tri);
            vote += n.dot(self.surface_normal(curved, centroid));
        }
        vote >= 0.0
    }

    /// The surface's own normal (the chart's `du × dv`) near a point.
    fn surface_normal(&self, curved: &Curved, p: Point) -> Vector {
        let Some(at) = chart(&curved.shape, p, self.tol) else {
            return Vector::ZERO;
        };
        let e = 1e-6;
        let o = evaluate(&curved.shape, at);
        let du = evaluate(&curved.shape, (at.0 + e, at.1)) - o;
        let dv = evaluate(&curved.shape, (at.0, at.1 + e)) - o;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m > 0.0 { n / m } else { Vector::ZERO }
    }

    fn curved_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(
                Construction,
                "a curved face was planned without its surface"
            );
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let mut wires: Vec<(f64, Shape)> = Vec::with_capacity(rings.len());
        for ring in rings {
            let area = ring_area(
                ring,
                self.triangles,
                |p| unwrapped(curved, p, self.tol).map(|(u, v)| Point2::new(u, v)),
                self.points,
            );
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
                let spec = &self.plan.edges[edge];
                if !self.has_pcurve(&edges[edge], surface) {
                    let Some((pcurve, deviation)) = self.plan.pcurves.get(&(edge, g)).cloned()
                    else {
                        ogeom_bail!(
                            Construction,
                            "an edge was planned without its image on a face"
                        );
                    };
                    self.model.widen(
                        &edges[edge],
                        Tolerance::new(deviation.max(self.tol.confusion()))?,
                    )?;
                    crate::build::attach_pcurve(
                        self.model,
                        &edges[edge],
                        pcurve,
                        surface,
                        Location::identity(),
                        spec.range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            let sign = if outward { 1.0 } else { -1.0 };
            wires.push((area * sign, self.model.add_wire(&ring_edges)?));
        }
        wires.sort_by(|a, b| b.0.total_cmp(&a.0));
        let wires: Vec<Shape> = wires.into_iter().map(|(_, w)| w).collect();
        // Fitted images answer on whichever branch their projection chose.
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A face round its axis between two rims of any shape, with holes.
    ///
    /// The seam joins a vertex of one rim to a vertex of the other along a
    /// straight line in the chart: the pair turning least between them
    /// whose line crosses no hole. It is a ruling where the pair stand at
    /// one angle on a cylinder or a cone, and otherwise the curve that line
    /// traces on the surface, interpolated. The outer wire walks the seam
    /// down, one rim round, the seam up a whole turn over, and the other rim
    /// back; each hole is its own wire.
    fn wrapped_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(
                Construction,
                "a face round its axis was planned without its surface"
            );
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        // Every edge's image, as the plan made it.
        for ring in rings {
            for (edge, _) in self.entries(ring) {
                if self.has_pcurve(&edges[edge], surface) {
                    continue;
                }
                let Some((pcurve, deviation)) = self.plan.pcurves.get(&(edge, g)).cloned() else {
                    ogeom_bail!(
                        Construction,
                        "an edge was planned without its image on a face"
                    );
                };
                self.model.widen(
                    &edges[edge],
                    Tolerance::new(deviation.max(self.tol.confusion()))?,
                )?;
                crate::build::attach_pcurve(
                    self.model,
                    &edges[edge],
                    pcurve,
                    surface,
                    Location::identity(),
                    self.plan.edges[edge].range,
                )?;
            }
        }
        let windings: Vec<i32> = rings
            .iter()
            .map(|ring| {
                winding(&curved.shape, ring, self.triangles, self.points, self.tol).unwrap_or(0)
            })
            .collect();
        let rims: Vec<usize> = (0..rings.len())
            .filter(|&k| windings[k].abs() == 1)
            .collect();
        let [low, high] = rims[..] else {
            ogeom_bail!(Construction, "a face round its axis has two rims");
        };
        if windings[low] != -windings[high] {
            ogeom_bail!(
                Construction,
                "a face round its axis has its rims turning one way"
            );
        }
        let holes: Vec<usize> = (0..rings.len()).filter(|&k| windings[k] == 0).collect();

        // Each rim's entries, and the vertex each one starts from.
        let starts = |this: &Self, ring: &[Half]| -> OgeomResult<Vec<RimStart>> {
            let mut out = Vec::new();
            for (edge, forward) in this.entries(ring) {
                let ends = this.model.children_of(&edges[edge])?;
                let vertex = if forward { ends.first() } else { ends.last() };
                let Some(vertex) = vertex.cloned() else {
                    ogeom_bail!(Construction, "a rim edge has no vertex");
                };
                let Some(ogeom_topo::NodeData::Vertex(data)) =
                    this.model.node(&vertex).map(|n| n.data())
                else {
                    ogeom_bail!(Construction, "a rim vertex has no position");
                };
                out.push(((edge, forward), vertex, data.point));
            }
            Ok(out)
        };
        let from = starts(self, &rings[low])?;
        let to = starts(self, &rings[high])?;
        let hole_rings = hole_polygons(
            &curved.shape,
            &holes
                .iter()
                .map(|&k| rings[k].as_slice())
                .collect::<Vec<_>>(),
            self.triangles,
            self.points,
            self.tol,
        );
        let from_at: Vec<Point> = from.iter().map(|x| x.2).collect();
        let to_at: Vec<Point> = to.iter().map(|x| x.2).collect();
        let Some((i, j, a, b)) =
            choose_seam(&curved.shape, &from_at, &to_at, &hole_rings, self.tol)
        else {
            ogeom_bail!(Construction, "no seam joins the rims clear of the holes");
        };
        let (pa, pb) = (from[i].2, to[j].2);
        let straight = (b.0 - a.0).abs() <= 1e-12
            && matches!(curved.shape, Canonical::Cylinder(_) | Canonical::Cone(_));
        let (seam_curve, range, deviation): (Curve, (f64, f64), f64) = if straight {
            let line = LineCurve::segment(pa, pb, self.tol)?;
            (line.into(), (0.0, pa.distance(pb)), self.tol.confusion())
        } else {
            const SAMPLES: u32 = 96;
            let along = |f: f64| (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
            let mut pts: Vec<Point> = (0..=SAMPLES)
                .map(|k| evaluate(&curved.shape, along(f64::from(k) / f64::from(SAMPLES))))
                .collect();
            pts[0] = pa;
            pts[SAMPLES as usize] = pb;
            let curve: Curve =
                crate::fit::interpolate(&pts, 3, crate::fit::Spacing::Uniform, self.tol)?.into();
            let range = curve.domain();
            let mut deviation = self.tol.confusion();
            for k in 0..=(SAMPLES * 4) {
                let f = f64::from(k) / f64::from(SAMPLES * 4);
                let t = range.0 + (range.1 - range.0) * f;
                deviation = deviation.max(
                    curve
                        .point_at(t, self.tol)?
                        .distance(evaluate(&curved.shape, along(f))),
                );
            }
            (curve, range, deviation)
        };
        let id = self.model.geometry_mut().add_curve(seam_curve);
        let mut data = EdgeData::on_curve(id, Location::identity(), range);
        data.tolerance = Tolerance::new(deviation)?;
        let seam = self
            .model
            .add_edge(data, &[from[i].1.clone(), to[j].1.clone()])?;
        // Down its near side where the wire leaves the high rim, up its far
        // side a whole turn over, where the low rim's walk comes round to.
        let over = tau * f64::from(windings[low]);
        let back = linear(a, b, range, self.tol)?;
        let forward = linear((a.0 + over, a.1), (b.0 + over, b.1), range, self.tol)?;
        crate::build::attach_seam(
            self.model,
            &seam,
            forward,
            back,
            surface,
            Location::identity(),
            range,
        )?;
        let rotated = |list: &[RimStart], at: usize| -> Vec<Shape> {
            (0..list.len())
                .map(|k| {
                    let ((edge, forward), _, _) = list[(at + k) % list.len()];
                    oriented(&edges[edge], forward)
                })
                .collect()
        };
        let mut outer = vec![seam.reversed()];
        outer.extend(rotated(&from, i));
        outer.push(seam.clone());
        outer.extend(rotated(&to, j));
        let mut wires = vec![self.model.add_wire(&outer)?];
        for &k in &holes {
            let ring_edges: Vec<Shape> = self
                .entries(&rings[k])
                .into_iter()
                .map(|(edge, forward)| oriented(&edges[edge], forward))
                .collect();
            wires.push(self.model.add_wire(&ring_edges)?);
        }
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A band round the axis: the two rim circles, and a seam at the
    /// surface's angle zero joining their vertices.
    fn band_face(
        &mut self,
        curved: &Curved,
        rings: &[Vec<Half>],
        edges: &[Shape],
        round_tube: bool,
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let g = self.groups.of[rings[0][0] / 3];
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a band was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let Some(frame) = axis_frame(&curved.shape) else {
            ogeom_bail!(Construction, "a band has no axis");
        };
        // Each rim: its edge, its chart height, and whether its parameter
        // runs with the surface's angle.
        let mut rims = Vec::with_capacity(2);
        for ring in rings {
            let (edge, _) = self.entry(ring[0]);
            let spec = &self.plan.edges[edge];
            let Curve::Circle(c) = &spec.curve else {
                ogeom_bail!(Construction, "a band's rim is not a circle");
            };
            let circle = c.circle();
            let start = spec.curve.point_at(0.0, self.tol)?;
            let Some((u, v)) = unwrapped(curved, start, self.tol) else {
                ogeom_bail!(Construction, "a band's rim has no chart position");
            };
            // Whether the rim's own parameter runs with the chart's angle
            // it goes round: the axis's for a band round it, the tube's for
            // a band round the tube.
            let turning = if round_tube {
                let radial = circle.centre() - frame.origin();
                radial.cross(frame.z().vector())
            } else {
                frame.z().vector()
            };
            let with = circle.frame().z().vector().dot(turning) > 0.0;
            rims.push((edge, if round_tube { u } else { v }, with, start));
        }
        rims.sort_by(|a, b| a.1.total_cmp(&b.1));
        let [
            (low, v_low, low_with, low_at),
            (high, v_high, high_with, high_at),
        ] = rims[..]
        else {
            ogeom_bail!(Construction, "a band has two rims");
        };
        for (edge, at, with) in [(low, v_low, low_with), (high, v_high, high_with)] {
            let (a, b) = if with { (0.0, tau) } else { (tau, 0.0) };
            let (from, to) = if round_tube {
                ((at, a), (at, b))
            } else {
                ((a, at), (b, at))
            };
            let pcurve = linear(from, to, (0.0, tau), self.tol)?;
            crate::build::attach_pcurve(
                self.model,
                &edges[edge],
                pcurve,
                surface,
                Location::identity(),
                (0.0, tau),
            )?;
        }
        // The seam, low rim to high, on the surface along angle zero: a
        // ruling, a meridian of a sphere, a circle of a torus's tube, or
        // for a band round the tube, an arc of its outer equator.
        let seam_curve: Curve = match curved.shape {
            Canonical::Torus(t) if round_tube => ogeom_geom::CircleCurve::new(
                ogeom_math::Circle::new(frame, t.major_radius() + t.minor_radius(), self.tol)?,
            )
            .into(),
            Canonical::Sphere(s) => {
                let normal =
                    Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
                ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
                    Frame::new(s.centre(), normal, frame.x(), self.tol)?,
                    s.radius(),
                    self.tol,
                )?)
                .into()
            }
            Canonical::Torus(t) => {
                let spine = frame.origin() + frame.x().vector() * t.major_radius();
                let normal =
                    Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
                ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
                    Frame::new(spine, normal, frame.x(), self.tol)?,
                    t.minor_radius(),
                    self.tol,
                )?)
                .into()
            }
            _ => LineCurve::segment(low_at, high_at, self.tol)?.into(),
        };
        let seam_range = match curved.shape {
            Canonical::Torus(_) | Canonical::Sphere(_) => (v_low, v_high),
            _ => (0.0, low_at.distance(high_at)),
        };
        // The seam runs between the rims' own vertices, placed where each
        // circle's parameter starts: the surface's angle zero.
        let placed = |at: Point| {
            self.plan
                .placed
                .iter()
                .any(|p| p.distance(at) <= self.tol.confusion())
        };
        if !placed(low_at) || !placed(high_at) {
            ogeom_bail!(
                Construction,
                "a band's rim vertex is not where its seam starts"
            );
        }
        let vertex = |edge: usize| -> OgeomResult<Shape> {
            match self.model.children_of(&edges[edge])?.first() {
                Some(v) => Ok(v.clone()),
                None => ogeom_bail!(Construction, "a rim has no vertex"),
            }
        };
        let (from_vertex, to_vertex) = (vertex(low)?, vertex(high)?);
        let id = self.model.geometry_mut().add_curve(seam_curve);
        let data = EdgeData::on_curve(id, Location::identity(), seam_range);
        let seam = self.model.add_edge(data, &[from_vertex, to_vertex])?;
        // The seam's two images: where the wire walks it forward, then back.
        let (forward, back) = if round_tube {
            (
                linear((v_low, 0.0), (v_high, 0.0), seam_range, self.tol)?,
                linear((v_low, tau), (v_high, tau), seam_range, self.tol)?,
            )
        } else {
            (
                linear((tau, v_low), (tau, v_high), seam_range, self.tol)?,
                linear((0.0, v_low), (0.0, v_high), seam_range, self.tol)?,
            )
        };
        crate::build::attach_seam(
            self.model,
            &seam,
            forward,
            back,
            surface,
            Location::identity(),
            seam_range,
        )?;
        // Counter-clockwise in the chart, then the whole ring the other way
        // for a face that faces against the surface. Round the axis: along
        // the low rim, up the seam's far side, back along the high rim, down
        // its near side. Round the tube: along the seam's near side, up the
        // high rim, back along the seam's far side, down the low rim.
        let mut ring = if round_tube {
            vec![
                seam.clone(),
                oriented(&edges[high], high_with),
                seam.reversed(),
                oriented(&edges[low], !low_with),
            ]
        } else {
            vec![
                oriented(&edges[low], low_with),
                seam.clone(),
                oriented(&edges[high], !high_with),
                seam.reversed(),
            ]
        };
        if !outward {
            ring.reverse();
            ring = ring.iter().map(Shape::reversed).collect();
        }
        let wire = self.model.add_wire(&ring)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, std::slice::from_ref(&wire))?;
        Ok(if outward { face } else { face.reversed() })
    }
}

impl Builder<'_> {
    /// A sphere's cap: its rim, a meridian seam from the rim to the pole,
    /// and the pole as an edge of no length, the frame's axis pointing into
    /// the cap so the pole is the north one.
    fn cap_face(
        &mut self,
        curved: &Curved,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let north = core::f64::consts::FRAC_PI_2;
        let g = self.groups.of[rings[0][0] / 3];
        let Canonical::Sphere(sphere) = curved.shape else {
            ogeom_bail!(Construction, "a cap is a sphere's");
        };
        let Some(geometry) = self.plan.surfaces[g].clone() else {
            ogeom_bail!(Construction, "a cap was planned without its surface");
        };
        let surface = self.model.geometry_mut().add_surface(geometry);
        let outward = self.outward(curved, g);
        let frame = sphere.frame();
        let (rim, _) = self.entry(rings[0][0]);
        let spec = &self.plan.edges[rim];
        let Curve::Circle(c) = &spec.curve else {
            ogeom_bail!(Construction, "a cap's rim is not a circle");
        };
        let with = c.circle().frame().z().vector().dot(frame.z().vector()) > 0.0;
        let start = spec.curve.point_at(0.0, self.tol)?;
        let Some((_, v_rim)) = unwrapped(curved, start, self.tol) else {
            ogeom_bail!(Construction, "a cap's rim has no chart position");
        };
        let (a, b) = if with { (0.0, tau) } else { (tau, 0.0) };
        crate::build::attach_pcurve(
            self.model,
            &edges[rim],
            linear((a, v_rim), (b, v_rim), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        let Some(rim_vertex) = self.model.children_of(&edges[rim])?.first().cloned() else {
            ogeom_bail!(Construction, "a rim has no vertex");
        };
        let pole = self.model.add_vertex(VertexData::new(
            sphere.centre() + frame.z().vector() * sphere.radius(),
        ));
        let normal = Direction::new(frame.x().vector().cross(frame.z().vector()), self.tol)?;
        let meridian: Curve = ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(
            Frame::new(sphere.centre(), normal, frame.x(), self.tol)?,
            sphere.radius(),
            self.tol,
        )?)
        .into();
        let id = self.model.geometry_mut().add_curve(meridian);
        let seam_range = (v_rim, north);
        let seam = self.model.add_edge(
            EdgeData::on_curve(id, Location::identity(), seam_range),
            &[rim_vertex, pole.clone()],
        )?;
        crate::build::attach_seam(
            self.model,
            &seam,
            linear((tau, v_rim), (tau, north), seam_range, self.tol)?,
            linear((0.0, v_rim), (0.0, north), seam_range, self.tol)?,
            surface,
            Location::identity(),
            seam_range,
        )?;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let tip = self.model.add_edge(data, &[pole.clone(), pole])?;
        crate::build::attach_pcurve(
            self.model,
            &tip,
            linear((0.0, north), (tau, north), (0.0, tau), self.tol)?,
            surface,
            Location::identity(),
            (0.0, tau),
        )?;
        // Counter-clockwise in the chart: along the rim, up the seam's far
        // side, back along the pole, down the seam's near side.
        let mut ring = vec![
            oriented(&edges[rim], with),
            seam.clone(),
            tip.reversed(),
            seam.reversed(),
        ];
        if !outward {
            ring.reverse();
            ring = ring.iter().map(Shape::reversed).collect();
        }
        let wire = self.model.add_wire(&ring)?;
        let mut data = FaceData::new(surface, Location::identity());
        data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        let face = self.model.add_face(data, std::slice::from_ref(&wire))?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A sphere or torus whole but for holes: the whole surface's face, a
    /// torus's seam round its axis on the parallel the plan placed clear of
    /// the holes, with each ring an inner wire. The rings
    /// are turned to run as the primitive's own wires do, so the face
    /// flips as one where the region faces against its surface.
    fn holed_face(
        &mut self,
        curved: &Curved,
        g: usize,
        rings: &[Vec<Half>],
        edges: &[Shape],
    ) -> OgeomResult<Shape> {
        let outward = self.outward(curved, g);
        let whole = self.whole_face(curved, g, torus_seam_v(curved))?;
        let whole = if outward { whole } else { whole.reversed() };
        let Some(ogeom_topo::NodeData::Face(data)) =
            self.model.node(&whole).map(|n| n.data().clone())
        else {
            ogeom_bail!(Construction, "a whole surface's face has no data");
        };
        let surface = data.surface;
        let mut wires = self.model.ordered_children_of(&whole)?;
        for ring in rings {
            let mut ring_edges = Vec::new();
            for (edge, forward) in self.entries(ring) {
                if !self.has_pcurve(&edges[edge], surface) {
                    let Some((pcurve, deviation)) = self.plan.pcurves.get(&(edge, g)).cloned()
                    else {
                        ogeom_bail!(
                            Construction,
                            "an edge was planned without its image on a face"
                        );
                    };
                    self.model.widen(
                        &edges[edge],
                        Tolerance::new(deviation.max(self.tol.confusion()))?,
                    )?;
                    crate::build::attach_pcurve(
                        self.model,
                        &edges[edge],
                        pcurve,
                        surface,
                        Location::identity(),
                        self.plan.edges[edge].range,
                    )?;
                }
                ring_edges.push(oriented(&edges[edge], forward));
            }
            if !outward {
                ring_edges.reverse();
                ring_edges = ring_edges.iter().map(Shape::reversed).collect();
            }
            wires.push(self.model.add_wire(&ring_edges)?);
        }
        crate::build::chain_wire_branches(self.model, surface, &wires, self.tol)?;
        let mut face_data = FaceData::new(surface, Location::identity());
        face_data.tolerance = data.tolerance;
        let face = self.model.add_face(face_data, &wires)?;
        Ok(if outward { face } else { face.reversed() })
    }

    /// A whole sphere or torus: the face the primitive builds on the
    /// recognized surface, a torus's seam round its axis on the parallel at
    /// tube angle `seam_v`.
    fn whole_face(&mut self, curved: &Curved, g: usize, seam_v: f64) -> OgeomResult<Shape> {
        let outward = self.outward(curved, g);
        let face = match curved.shape {
            Canonical::Sphere(s) => {
                let built =
                    crate::primitive::make_sphere(self.model, s.frame(), s.radius(), self.tol)?;
                let Some(face) = ogeom_topo::explore_unique(
                    self.model,
                    &built.shape,
                    ogeom_topo::ShapeType::Face,
                )?
                .into_iter()
                .next() else {
                    ogeom_bail!(Construction, "a primitive came back with no face");
                };
                face
            }
            Canonical::Torus(t) => crate::primitive::torus_face_seamed_at(
                self.model,
                t.frame(),
                t.major_radius(),
                t.minor_radius(),
                seam_v,
                self.tol,
            )?,
            _ => ogeom_bail!(Construction, "only a sphere or a torus is whole"),
        };
        // The facets stand off the surface by up to the fit's deviation, and
        // so do the seam and the poles built on it: the face's tolerance
        // reaches down to every edge and vertex it holds, as the checker
        // requires of a face looser than its boundary.
        self.model.widen(
            &face,
            Tolerance::new(curved.deviation.max(self.tol.confusion()))?,
        )?;
        Ok(if outward { face } else { face.reversed() })
    }
}

fn oriented(edge: &Shape, forward: bool) -> Shape {
    if forward {
        edge.clone()
    } else {
        edge.reversed()
    }
}

/// A straight chart segment from `a` to `b`, linear in the edge's own
/// parameter over `range`.
fn linear(
    a: (f64, f64),
    b: (f64, f64),
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let knots = ogeom_math::KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?;
    Ok(ogeom_geom::BSpline2d::new(
        knots,
        vec![Point2::new(a.0, a.1), Point2::new(b.0, b.1)],
        tol,
    )?
    .into())
}

/// Where a rim's entry starts: the entry, its vertex, and where that
/// stands.
type RimStart = ((usize, bool), Shape, Point);

/// A seam's two rim vertices, by their places in each rim, and the ends of
/// its line in the chart.
type SeamChoice = (usize, usize, (f64, f64), (f64, f64));

/// The tube angle of a torus's seam round its axis: half a turn from its
/// chart's centre.
fn torus_seam_v(curved: &Curved) -> f64 {
    curved.centre.1 - core::f64::consts::PI
}

/// A face's holes in its chart: each a polygon of its mesh vertices,
/// carried round continuously.
fn hole_polygons(
    shape: &Canonical,
    holes: &[&[Half]],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Vec<Vec<(f64, f64)>> {
    holes
        .iter()
        .filter_map(|ring| {
            let mut out: Vec<(f64, f64)> = Vec::new();
            for &h in *ring {
                let (a, _) = from_to(triangles, h);
                let (u, v) = chart(shape, points[a as usize], tol)?;
                // Unwrapped along the ring both ways the chart closes on
                // itself: a hole across a torus's outer equator, where the
                // tube's angle starts, runs on past a whole turn rather
                // than jumping back to nought, so a seam there is seen to
                // cross it.
                let (u, v) = match out.last() {
                    Some(&(last_u, last_v)) => (
                        last_u + ogeom_math::elementary::wrap_signed_angle(u - last_u),
                        if matches!(shape, Canonical::Torus(_)) {
                            last_v + ogeom_math::elementary::wrap_signed_angle(v - last_v)
                        } else {
                            v
                        },
                    ),
                    None => (u, v),
                };
                out.push((u, v));
            }
            Some(out)
        })
        .collect()
}

/// The seam for a face round its axis: of the pairs of a vertex on one
/// rim and a vertex on the other, the one turning least between them whose
/// straight chart line crosses no hole, a whole turn either way included.
/// Its indices and the line's ends in the chart.
fn choose_seam(
    shape: &Canonical,
    from: &[Point],
    to: &[Point],
    holes: &[Vec<(f64, f64)>],
    tol: Tolerances,
) -> Option<SeamChoice> {
    let tau = core::f64::consts::TAU;
    let crosses = |a: (f64, f64), b: (f64, f64)| {
        holes.iter().any(|ring| {
            [-tau, 0.0, tau].iter().any(|shift| {
                (0..ring.len()).any(|i| {
                    let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                    segments_cross(a, b, (p.0 + shift, p.1), (q.0 + shift, q.1))
                })
            })
        })
    };
    let mut best: Option<(f64, SeamChoice)> = None;
    for (i, pa) in from.iter().enumerate() {
        let Some((ua, va)) = chart(shape, *pa, tol) else {
            continue;
        };
        for (j, pb) in to.iter().enumerate() {
            let Some((ub, vb)) = chart(shape, *pb, tol) else {
                continue;
            };
            let turn = ogeom_math::elementary::wrap_signed_angle(ub - ua);
            let (a, b) = ((ua, va), (ua + turn, vb));
            if crosses(a, b) {
                continue;
            }
            let score = turn.abs() * 1e3 + pa.distance(*pb);
            if best.is_none_or(|held| score < held.0) {
                best = Some((score, (i, j, a, b)));
            }
        }
    }
    best.map(|(_, choice)| choice)
}

/// The middle of the widest stretch of angle no hole covers.
fn free_angle(holes: &[Vec<(f64, f64)>]) -> Option<f64> {
    let tau = core::f64::consts::TAU;
    let mut angles: Vec<f64> = holes
        .iter()
        .flatten()
        .map(|(u, _)| u.rem_euclid(tau))
        .collect();
    if angles.is_empty() {
        return None;
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

/// Whether two chart segments cross, each at a point strictly inside both.
fn segments_cross(a: (f64, f64), b: (f64, f64), p: (f64, f64), q: (f64, f64)) -> bool {
    let side = |o: (f64, f64), x: (f64, f64), y: (f64, f64)| {
        (x.0 - o.0).mul_add(y.1 - o.1, -((x.1 - o.1) * (y.0 - o.0)))
    };
    let (d1, d2) = (side(p, q, a), side(p, q, b));
    let (d3, d4) = (side(a, b, p), side(a, b, q));
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

/// How many times a ring of half-edges goes round a surface's angle, with
/// the sign of its sense; `None` where a vertex has no chart position.
fn winding(
    shape: &Canonical,
    ring: &[Half],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<i32> {
    let mut turned = 0.0;
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        let (ua, _) = chart(shape, points[a as usize], tol)?;
        let (ub, _) = chart(shape, points[b as usize], tol)?;
        turned += ogeom_math::elementary::wrap_signed_angle(ub - ua);
    }
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    Some((turned / core::f64::consts::TAU).round() as i32)
}

/// How many times a ring goes round each of a surface's chart directions,
/// as [`winding`] counts the first.
fn windings(
    shape: &Canonical,
    ring: &[Half],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<(i32, i32)> {
    let (_, wraps_v) = periodic(shape);
    let mut turned = (0.0, 0.0);
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        let (ua, va) = chart(shape, points[a as usize], tol)?;
        let (ub, vb) = chart(shape, points[b as usize], tol)?;
        turned.0 += ogeom_math::elementary::wrap_signed_angle(ub - ua);
        if wraps_v {
            turned.1 += ogeom_math::elementary::wrap_signed_angle(vb - va);
        }
    }
    let tau = core::f64::consts::TAU;
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    Some((
        (turned.0 / tau).round() as i32,
        (turned.1 / tau).round() as i32,
    ))
}

/// Twice the area a ring of half-edges encloses in a chart, signed.
fn ring_area(
    ring: &[Half],
    triangles: &[[u32; 3]],
    chart: impl Fn(Point) -> Option<Point2>,
    points: &[Point],
) -> f64 {
    let mut area = 0.0;
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        if let (Some(a), Some(b)) = (chart(points[a as usize]), chart(points[b as usize])) {
            area += a.x * b.y - b.x * a.y;
        }
    }
    area
}

/// The distance from `p` to the segment from `a` to `b`.
fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.dot(along);
    if length <= 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

fn distance_to_line(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let m = d.magnitude();
    if m == 0.0 {
        return p.distance(a);
    }
    (p - a).cross(d).magnitude() / m
}
