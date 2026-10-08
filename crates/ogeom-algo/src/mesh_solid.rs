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

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Plane, Point};
use ogeom_topo::{Model, Shape, Triangulation};

use crate::recognize::Canonical;

mod astray;
mod bodies;
mod builder;
mod caches;
mod checks;
mod frames;
mod holed;
mod planner;
mod regions;
mod rounds;
mod seams;
mod segment;
mod snap;
mod steps;
mod weld;

pub use steps::{
    FallbackReason, FitConstraints, MeshRegion, MeshRegions, RegionFallback, RegionId,
    RegionRefusal, SurfaceKind,
};

use astray::astray_faces;
pub(crate) use astray::distance_to_triangle;
use bodies::{assemble, body_culprits};
use builder::Builder;
use caches::{AreaCache, ImageCache, SnapCache};
use checks::{
    Regrouped, absorb_facets, any_turned_in, backwards_facets, collapsed_beside, crossed_seams,
    failing_faces, folded_seams, give_back, overlapping_faces, planes_off_their_vertices,
    unbuilt_culprits, unmatched_faces,
};
use planner::{NOISE_FLOOR, Planner, Replan};
use regions::{flat_noise, split_disconnected};
use segment::{coplanar_groups, found_within, one_each, plane_of, segment};
use weld::{Adjacency, Piece, diagonal, has_area, inside, orient, unfold, weld_points};

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
    /// mean-value map of the region onto a square. Curved regions the
    /// smooth area encloses (meeting nothing outside it except across
    /// creases), as a fine mesh of a free-form surface is cut into, are
    /// tried with the smooth regions they join as one region first, and
    /// kept apart where that patch does not verify. A region that is not a
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
    /// Regions whose face could not be built on their surface, and which
    /// were faceted instead: curved regions, and the few planar ones whose
    /// face failed or whose plane does not hold their vertices, their
    /// triangles gathered again.
    pub curved_faceted: usize,
    /// Those regions by name, each with why it was faceted, in the order
    /// they fell back: as many as [`MeshSolidReport::curved_faceted`]
    /// counts. The names are the regions' in
    /// [`MeshRegions`], whose [`MeshRegions::find`] names them alike for
    /// the same mesh and options.
    pub fallbacks: Vec<RegionFallback>,
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
    /// [`MeshSolidReport::coplanar_distance_raised`]). Where the build put
    /// back every region the default found, it is the default.
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
    let mut on_edge: ogeom_core::FastMap<(u32, u32), Vec<usize>> = ogeom_core::FastMap::default();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            on_edge
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(i);
        }
    }
    let unlink =
        |on_edge: &mut ogeom_core::FastMap<(u32, u32), Vec<usize>>, t: [u32; 3], i: usize| {
            for k in 0..3 {
                if let Some(list) = on_edge.get_mut(&key(t[k], t[(k + 1) % 3]))
                    && let Ok(at) = list.binary_search(&i)
                {
                    list.remove(at);
                }
            }
        };
    let link =
        |on_edge: &mut ogeom_core::FastMap<(u32, u32), Vec<usize>>, t: [u32; 3], i: usize| {
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
/// inside another becomes a void of the solid around it. It is
/// [`MeshRegions::find`] followed by [`MeshRegions::build`], with no region
/// changed between them.
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
    MeshRegions::find(mesh, options, tol)?.build(model)
}

/// A mesh welded, cleaned and oriented, its pieces found, and its triangles
/// gathered into regions: what a solid is built from.
#[derive(Clone)]
struct Found {
    points: Vec<Point>,
    /// The welded vertex each of the input mesh's vertices became.
    remap: Vec<u32>,
    triangles: Vec<[u32; 3]>,
    adjacency: Adjacency,
    pieces: Vec<Piece>,
    /// How many closed pieces each piece lies inside.
    depth: Vec<usize>,
    all_closed: bool,
    /// The coplanar distance the regions were found to.
    flat: f64,
    /// The default distance, where the regions were found at a wider one.
    tight: f64,
    /// The regions found at the default, where they were found again at a
    /// wider distance.
    first: Option<Box<Groups>>,
    groups: Groups,
    report: MeshSolidReport,
}

/// Weld, clean and orient the mesh, and gather its triangles into regions.
fn find(mesh: &Triangulation, options: &MeshSolidOptions, tol: Tolerances) -> OgeomResult<Found> {
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
    let mut seen: FastMap<[u32; 3], ()> =
        FastMap::with_capacity_and_hasher(mesh.triangles.len(), Default::default());
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
    // distance, it widens while that holds, twice at most. A wider distance
    // also lets a fit spread over triangles of other surfaces, so each
    // curved region found at it keeps the regions the default found among
    // its triangles, which the build puts back where it cannot place the
    // wider one.
    let found_at = flat;
    let mut first = None;
    if options.coplanar_distance.is_none() && options.recognize {
        first = Some(Box::new(groups.clone()));
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
        match &first {
            Some(narrow) if flat > found_at => found_within(&mut groups, narrow, flat),
            _ => first = None,
        }
    }
    Ok(Found {
        points,
        remap,
        triangles,
        adjacency,
        pieces,
        depth,
        all_closed,
        flat,
        tight: found_at,
        first,
        groups,
        report,
    })
}

/// The solid built from found regions.
///
/// A region whose carrier is gone and still holds triangles gives them to
/// planar faces first, as the planar pass gathers them. A facet holding a
/// triangle marked in `protected` (one per triangle, or empty for none) is
/// never given to a curved neighbour. A region found at a wider distance
/// than the default is tried again as the default found it before it is
/// faceted: the build starts over from the regions found, with it put
/// back.
fn build(
    model: &mut Model,
    found: &Found,
    options: &MeshSolidOptions,
    protected: &[bool],
    tol: Tolerances,
) -> OgeomResult<MeshSolid> {
    let (triangles, adjacency) = (&found.triangles, &found.adjacency);
    let mut start = found.groups.clone();
    let mut flat = found.flat;
    let mut first = found.first.clone();
    loop {
        let narrowing = match build_from(model, found, &start, flat, options, protected, tol)? {
            Ok(solid) => return Ok(solid),
            Err(narrowing) => narrowing,
        };
        for g in narrowing {
            narrow(&mut start, triangles, g);
        }
        // With every region put back, the regions are the ones the default
        // found.
        if start.narrower.is_empty()
            && let Some(regions) = first.take()
        {
            start = *regions;
        } else {
            split_disconnected(triangles, adjacency, &mut start);
            coplanar_groups(
                &found.points,
                triangles,
                adjacency,
                options,
                found.tight,
                &mut start,
                tol,
            )?;
        }
        flat = start
            .narrower
            .values()
            .map(|n| n.distance)
            .fold(found.tight, f64::max);
    }
}

/// [`build`] from regions `start`, to the coplanar distance `flat`: the
/// solid, or the regions found at a wider distance than the default that
/// it would facet.
#[allow(clippy::too_many_arguments, reason = "the build's inputs")]
fn build_from(
    model: &mut Model,
    found: &Found,
    start: &Groups,
    flat: f64,
    options: &MeshSolidOptions,
    protected: &[bool],
    tol: Tolerances,
) -> OgeomResult<Result<MeshSolid, Vec<usize>>> {
    let (points, triangles, adjacency) = (&found.points, &found.triangles, &found.adjacency);
    let (pieces, depth, all_closed) = (&found.pieces, &found.depth, found.all_closed);
    let mut groups = start.clone();
    let mut report = found.report.clone();
    // A plane that does not hold its region's vertices (one a caller put
    // there) can place none of its boundary; its triangles are gathered
    // again.
    let off = planes_off_their_vertices(points, triangles, &groups, flat);
    let mut released = !off.is_empty();
    for g in off {
        withdraw(
            &mut groups,
            &mut report,
            g,
            FallbackReason::BoundaryNotPlaced,
        );
    }
    for t in 0..triangles.len() {
        if matches!(groups.carriers.get(groups.of[t]), Some(Carrier::Gone)) {
            groups.of[t] = usize::MAX;
            released = true;
        }
    }
    if released {
        if options.merge_coplanar {
            coplanar_groups(
                points,
                triangles,
                adjacency,
                options,
                flat,
                &mut groups,
                tol,
            )?;
        } else {
            one_each(points, triangles, &mut groups, tol)?;
        }
    }
    // Facets left between curved faces go to them first. Whether that
    // leaves every face facing out is known only once the solid is built;
    // where one faces into the material, the conversion is made again
    // without them.
    let unabsorbed = (groups.clone(), report.clone());
    let mut absorbed = if options.recognize {
        absorb_facets(points, triangles, adjacency, &mut groups, protected, flat)
    } else {
        FastMap::default()
    };
    let mut absorbing = !absorbed.is_empty();
    let mut regrouped = Regrouped::default();
    let snaps = std::sync::Mutex::new(SnapCache::default());
    let images = std::sync::Mutex::new(ImageCache::default());
    // Face meshes and areas kept across the builds: most faces recur
    // unchanged.
    let kept = crate::mass::VolumeKept::default();
    let areas = AreaCache::default();
    let shape = 'attempt: loop {
        // Plan until every curved face's boundary is exact, faceting the ones
        // whose boundary is not; then build, and facet any recognized face that
        // reaches past the triangles it replaces (a boundary placed on the wrong
        // turn of its surface closes a face of the wrong extent) and build again.
        let mut pinned: ogeom_core::FastSet<u32> = ogeom_core::FastSet::default();
        let mut straight: ogeom_core::FastSet<(u32, u32)> = ogeom_core::FastSet::default();
        // What stood before any seam was threaded straight: a straightened
        // build with a face collapsed beside a threaded seam is set aside for
        // it. A face turned into the material beside one is an overlap like
        // any other, and the curved faces beside it are faceted.
        let mut unthreaded: Option<(Groups, MeshSolidReport)> = None;
        let mut threading_refused = false;
        let shape = 'build: loop {
            let fans = groups.fans_by_group(triangles, adjacency);
            let planner = Planner {
                points,
                triangles,
                adjacency,
                groups: &groups,
                merge: options.merge_coplanar && !options.keep_vertices,
                pinned: &pinned,
                straight: &straight,
                crease: options.crease,
                flat,
                tol,
                snaps: &snaps,
                images: &images,
                fans: &fans,
            };
            // A step that fails on the built faces withdraws the faces that
            // cannot be meshed or measured on their own (or, where none is
            // found, every recognized one), and the error stands only where
            // nothing is left to withdraw.
            let mut step_error = None;
            macro_rules! or_withdraw {
                ($label:lifetime, $built:expr, $step:expr) => {
                    match $step {
                        Ok(value) => value,
                        Err(error) => {
                            step_error = Some(error);
                            break $label(
                                failing_faces(model, &groups, $built, tol),
                                FallbackReason::BuildFailed,
                            );
                        }
                    }
                };
            }
            let failed = match planner.plan() {
                Err(error) => {
                    step_error = Some(error);
                    (
                        failing_faces(model, &groups, &[], tol),
                        FallbackReason::BuildFailed,
                    )
                }
                Ok(Err(Replan::Pin(vertices))) => {
                    pinned.extend(vertices);
                    continue;
                }
                Ok(Err(Replan::Fan(facets))) => {
                    groups.fans.extend(facets);
                    continue;
                }
                Ok(Err(Replan::Facet(failed))) => (failed, FallbackReason::BoundaryNotPlaced),
                Ok(Ok(plan)) => 'checks: {
                    model.begin_operation();
                    let (built, unbuilt) = Builder {
                        model,
                        points,
                        triangles,
                        groups: &groups,
                        plan: &plan,
                        fans: &fans,
                        tol,
                    }
                    .build();
                    if !unbuilt.is_empty() {
                        break 'checks (
                            unbuilt_culprits(&unbuilt, &groups, adjacency, &fans),
                            FallbackReason::BuildFailed,
                        );
                    }
                    let astray = or_withdraw!(
                        'checks,
                        &built,
                        astray_faces(model, points, triangles, &groups, &built, flat, kept.meshes(), tol)
                    );
                    if astray.is_empty() {
                        let (shape, bodies) = or_withdraw!(
                            'checks,
                            &built,
                            assemble(
                                model, points, triangles, pieces, depth, all_closed, &groups,
                                &built,
                            )
                        );
                        let (mut culprits, mut reason) = if options.recognize && all_closed {
                            or_withdraw!(
                                'checks,
                                &built,
                                body_culprits(
                                    model,
                                    points,
                                    triangles,
                                    &groups,
                                    &built,
                                    &bodies,
                                    flat,
                                    &kept,
                                    tol,
                                    &mut report,
                                )
                            )
                        } else {
                            (Vec::new(), FallbackReason::TurnedIn)
                        };
                        if culprits.is_empty() && options.recognize && !threading_refused {
                            let crossed = or_withdraw!(
                                'checks,
                                &built,
                                crossed_seams(
                                    model, &shape, triangles, adjacency, &groups, &built,
                                    (kept.meshes(), &areas), tol,
                                )
                            );
                            let before = straight.len();
                            if before == 0 && !crossed.is_empty() {
                                unthreaded = Some((groups.clone(), report.clone()));
                            }
                            straight.extend(crossed);
                            if straight.len() > before {
                                continue 'build;
                            }
                            if let Some((was, then)) = unthreaded.take()
                                && or_withdraw!(
                                    'checks,
                                    &built,
                                    collapsed_beside(
                                        model, triangles, &groups, &built, &straight, &areas,
                                        tol,
                                    )
                                )
                            {
                                groups = was;
                                report = then;
                                straight.clear();
                                threading_refused = true;
                                continue 'build;
                            }
                        }
                        if culprits.is_empty() && options.recognize {
                            let backwards = or_withdraw!(
                                'checks,
                                &built,
                                backwards_facets(
                                    model, points, triangles, adjacency, &groups, &built, &fans,
                                    tol,
                                )
                            );
                            if !backwards.is_empty() {
                                groups.fans.extend(backwards);
                                continue 'build;
                            }
                        }
                        if culprits.is_empty() && options.recognize {
                            let fanned;
                            (culprits, fanned) = or_withdraw!(
                                'checks,
                                &built,
                                folded_seams(model, triangles, adjacency, &groups, &built, tol)
                            );
                            // A facet folding a curved face beside it is
                            // built as a fan before that face is faceted.
                            if !fanned.is_empty() {
                                groups.fans.extend(fanned);
                                continue 'build;
                            }
                            reason = FallbackReason::FoldedSeam;
                        }
                        if culprits.is_empty() && options.recognize {
                            let (turned, alone, covering);
                            (culprits, turned, alone, covering) = or_withdraw!(
                                'checks,
                                &built,
                                overlapping_faces(model, &shape, adjacency, &groups, &built, tol)
                            );
                            reason = FallbackReason::Overlaps;
                            // A facet turned in beside one curved face is
                            // built as a fan before that face is faceted.
                            let fanned: Vec<usize> = turned
                                .iter()
                                .filter(|&&g| !fans.contains_key(&g))
                                .filter_map(|&g| {
                                    let t = groups.of.iter().position(|&o| o == g)?;
                                    groups.fan_at(t, triangles, adjacency).map(|_| t)
                                })
                                .collect();
                            if !fanned.is_empty() {
                                groups.fans.extend(fanned);
                                continue 'build;
                            }
                            // A curved face that took facets in can run
                            // over the planar faces near it, which then face
                            // into it: every such face across a turned-in
                            // one gives its facets back before the curved
                            // face beside it is blamed.
                            let covered: Vec<usize> = covering
                                .into_iter()
                                .filter(|g| absorbed.contains_key(g))
                                .collect();
                            if !covered.is_empty() {
                                for g in covered {
                                    if let Some(record) = absorbed.remove(&g) {
                                        give_back(&mut groups, g, record);
                                    }
                                }
                                continue 'build;
                            }
                            // A planar face turned in with no recognized
                            // face near it is gathered again, then built a
                            // face per triangle, while that is left to try.
                            if culprits.is_empty() {
                                culprits = alone
                                    .into_iter()
                                    .filter(|&g| regrouped.can_withdraw(&groups, g))
                                    .collect();
                                reason = FallbackReason::TurnedIn;
                            }
                        }
                        if culprits.is_empty() && options.recognize {
                            culprits = or_withdraw!(
                                'checks,
                                &built,
                                unmatched_faces(model, &shape, &groups, &built, flat, kept.meshes(), tol)
                            );
                            reason = FallbackReason::MeshesOpen;
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
                            break 'build shape;
                        }
                        (culprits, reason)
                    } else {
                        // A facet all but tangent to a curved face beside
                        // it meets the face's surface along a seam that can
                        // lie well past the facet, and the face then reaches
                        // past its triangles: each facet beside such a face
                        // that can be a fan is built as one onto it before
                        // the face is faceted.
                        let fanned: Vec<usize> = (0..triangles.len())
                            .filter(|&t| !fans.contains_key(&groups.of[t]))
                            .filter(|&t| {
                                groups.fan_at(t, triangles, adjacency).is_some_and(|fan| {
                                    fan.seams().any(|(_, c)| astray.contains(&c))
                                })
                            })
                            .collect();
                        if !fanned.is_empty() {
                            groups.fans.extend(fanned);
                            continue 'build;
                        }
                        (astray, FallbackReason::ReachesPast)
                    }
                }
            };
            let (failed, reason) = failed;
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
            let narrowing: Vec<usize> = failed
                .iter()
                .copied()
                .filter(|g| groups.narrower.contains_key(g))
                .collect();
            if !narrowing.is_empty() {
                return Ok(Err(narrowing));
            }
            let mut alone = Vec::new();
            let mut withdrawn = false;
            for g in failed {
                if !regrouped.can_withdraw(&groups, g) {
                    continue;
                }
                if matches!(groups.carriers[g], Carrier::Plane(_)) {
                    alone.extend(regrouped.withdraw_plane(&groups, g));
                }
                withdraw(&mut groups, &mut report, g, reason);
                withdrawn = true;
            }
            if !withdrawn {
                if let Some(error) = step_error {
                    return Err(error);
                }
                ogeom_bail!(Construction, "a facet of the mesh could not be built");
            }
            for t in alone {
                groups.of[t] = groups.carriers.len();
                groups
                    .carriers
                    .push(Carrier::Plane(plane_of(points, triangles[t], tol)?));
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
        };
        if absorbing && any_turned_in(model, &shape, tol)? {
            absorbing = false;
            absorbed.clear();
            (groups, report) = unabsorbed.clone();
            continue 'attempt;
        }
        break shape;
    };
    crate::pcurve_gap::state_pcurve_gaps(model, &shape, tol)?;
    Ok(Ok(MeshSolid {
        shape,
        closed: all_closed,
        coplanar_distance: flat,
        report,
    }))
}

/// Facet region `g`: its face is not built on its surface, its triangles
/// are left for the planar faces to gather, and the report names it.
fn withdraw(groups: &mut Groups, report: &mut MeshSolidReport, g: usize, reason: FallbackReason) {
    groups.carriers[g] = Carrier::Gone;
    report.curved_faceted += 1;
    report.fallbacks.push(RegionFallback {
        region: RegionId(g),
        reason,
    });
    for of in &mut groups.of {
        if *of == g {
            *of = usize::MAX;
        }
    }
}

/// Put back in place of region `g`, found at a wider distance, the curved
/// regions the default found within it: each takes back its triangles,
/// from `g`, from a planar face (gathered again round it) or from another
/// curved region found wider, and what `g` held beyond them is left for
/// the planar faces. Nothing changes where `g` was not found wider.
fn narrow(groups: &mut Groups, triangles: &[[u32; 3]], g: usize) {
    let Some(record) = groups.narrower.remove(&g) else {
        return;
    };
    groups.carriers[g] = Carrier::Gone;
    let mut touched = vec![g];
    for (_, held) in &record.regions {
        for &t in held {
            let p = groups.of[t];
            if matches!(groups.carriers.get(p), Some(Carrier::Plane(_))) {
                groups.carriers[p] = Carrier::Gone;
            }
        }
    }
    for of in &mut groups.of {
        if *of == g || matches!(groups.carriers.get(*of), Some(Carrier::Gone)) {
            *of = usize::MAX;
        }
    }
    for (k, (carrier, held)) in record.regions.into_iter().enumerate() {
        let at = if k == 0 {
            groups.carriers[g] = carrier;
            g
        } else {
            groups.carriers.push(carrier);
            groups.carriers.len() - 1
        };
        touched.push(at);
        for t in held {
            let p = groups.of[t];
            if p == usize::MAX || groups.narrower.contains_key(&p) {
                touched.push(p);
                groups.of[t] = at;
            }
        }
    }
    touched.sort_unstable();
    touched.dedup();
    for c in touched {
        if c == usize::MAX {
            continue;
        }
        let mut vertices: Vec<u32> = groups
            .of
            .iter()
            .enumerate()
            .filter(|&(_, &o)| o == c)
            .flat_map(|(t, _)| triangles[t])
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        if let Carrier::Curved(curved) = &mut groups.carriers[c] {
            curved.vertices = vertices;
        }
    }
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
    /// A sphere's cap: one latitude circle, a seam to the pole, and the
    /// pole; or a cone closing at its apex: one rim circle, a ruling to the
    /// apex, and the apex.
    Cap,
    /// A sphere's cap with holes: one rim of any shape going round the
    /// axis, a seam from a vertex of it to the pole, straight in the chart
    /// and clear of the holes, the pole, and the holes as inner wires.
    HoledCap,
    /// The whole surface, with no boundary of its own.
    Whole,
    /// Round the axis (or round a torus's tube) between two rims of any
    /// shape, with holes: a seam joins a vertex of each rim, straight in
    /// the chart and clear of the holes.
    Wrapped,
    /// A sphere or a torus whole but for holes, its seams and poles clear
    /// of them: the whole surface's face, the holes its inner wires.
    Holed,
    /// A torus whole but for holes that between them leave no parallel
    /// (or no meridian) free: its seam that way runs through some of
    /// them (see [`Thread`]).
    Threaded,
}

/// The seam of a torus whose holes leave no parallel (or no meridian)
/// free, read in a working chart whose first angle is the one the seam
/// goes round (the axis's, or the tube's where `round_tube`). The seam is a
/// chain: from a point on the free circle the other way, straight in the
/// chart to a hole's vertex furthest back along it, round the hole to its
/// vertex furthest on, straight to the next hole, and so on back to the
/// start. Each straight piece is a seam edge, met once from either side.
/// Each hole it passes through is split at those two vertices: one arc
/// bounds the face along the chain, the other a whole turn across.
#[derive(Debug, Clone)]
struct Thread {
    round_tube: bool,
    /// The free circle's first angle, where the chain starts and ends.
    start: (f64, f64),
    /// The holes the chain passes through, in order: the ring, the mesh
    /// vertices it enters and leaves at and where those stand in the
    /// working chart.
    stops: Vec<Stop>,
}

/// One hole a [`Thread`] passes through.
#[derive(Debug, Clone, Copy)]
struct Stop {
    ring: usize,
    enter: u32,
    leave: u32,
    enter_at: (f64, f64),
    leave_at: (f64, f64),
}

/// Triangles gathered into faces: which face each triangle is in, and what
/// each face is built on.
#[derive(Clone)]
struct Groups {
    of: Vec<usize>,
    carriers: Vec<Carrier>,
    /// The smooth regions refused a patch, by why.
    refused: PatchRefusals,
    /// Triangles each built as a fan (see [`Fan`]) rather than flat, while
    /// they still stand alone in their group beside one curved face.
    fans: std::collections::BTreeSet<usize>,
    /// The curved regions found at a wider distance than the default, each
    /// with what the default found among its triangles.
    narrower: FastMap<usize, Narrower>,
}

/// What the default distance found among the triangles of a region found
/// at a wider one.
#[derive(Clone)]
struct Narrower {
    /// The distance the region was found at.
    distance: f64,
    /// The curved regions found at the default holding any of its
    /// triangles, each with all of its own.
    regions: Vec<(Carrier, Vec<usize>)>,
}

/// A planar facet of one triangle beside one curved face, built as the
/// ruled surface from its far corner to its seam lifted onto the curved
/// surface. Where a facet meets a curved face all but tangentially (a chord
/// of a neighbour left faceted), the plane through it meets the curved
/// surface along a curve that bulges past the facet's third corner, or not
/// along the seam at all; a seam threaded straight instead trims the curved
/// face short of it by the seam's sag. The fan shares the lifted seam with
/// the curved face and its two straight sides with the faces beside it,
/// exactly, and closes at its apex on an edge with no length.
///
/// A facet with two sides on curved faces and the third on a planar one is
/// a wedge: built as the ruled surface between its two seams, each lifted
/// onto its curved surface, from the corner they share (where it closes on
/// an edge with no length) to its third side, the ruling at their far ends.
#[derive(Debug, Clone, Copy)]
struct Fan {
    /// The triangle's corner off the seam.
    apex: u32,
    /// The seam's ends, in the triangle's winding.
    seam: (u32, u32),
    /// The curved face across the seam.
    curved: usize,
    /// For a wedge, the curved face across its second seam, from the
    /// seam's second end to the apex.
    across: Option<usize>,
}

impl Fan {
    /// Its seams, each by its ends and the curved face across it.
    fn seams(&self) -> impl Iterator<Item = ((u32, u32), usize)> + use<> {
        let second = self.across.map(|c| ((self.seam.1, self.apex), c));
        core::iter::once((self.seam, self.curved)).chain(second)
    }

    /// The curved face across the seam between `p` and `q`, either way.
    fn seam_between(&self, p: u32, q: u32) -> Option<usize> {
        self.seams()
            .find(|&((a, b), _)| (a, b) == (p, q) || (b, a) == (p, q))
            .map(|(_, c)| c)
    }
}

impl Groups {
    /// The fan a triangle would be built as: alone in its planar group, one
    /// side against a curved face and the other two against planar ones, or
    /// a wedge, two sides against curved faces and the third against a
    /// planar one.
    fn fan_at(&self, t: usize, triangles: &[[u32; 3]], adjacency: &Adjacency) -> Option<Fan> {
        let g = *self.of.get(t)?;
        if !matches!(self.carriers.get(g), Some(Carrier::Plane(_))) {
            return None;
        }
        let mut sides = [None; 3];
        for (k, side) in sides.iter_mut().enumerate() {
            let o = self.of[adjacency.twin[3 * t + k]? / 3];
            match self.carriers.get(o) {
                _ if o == g => return None,
                Some(Carrier::Curved(_)) => *side = Some(o),
                Some(Carrier::Plane(_)) => {}
                _ => return None,
            }
        }
        // The seam is the curved side; a wedge's, the curved side whose
        // successor in the winding is curved too.
        let k = match sides.iter().flatten().count() {
            1 => sides.iter().position(Option::is_some)?,
            2 => (0..3).find(|&k| sides[k].is_some() && sides[(k + 1) % 3].is_some())?,
            _ => return None,
        };
        let tri = triangles[t];
        Some(Fan {
            apex: tri[(k + 2) % 3],
            seam: (tri[k], tri[(k + 1) % 3]),
            curved: sides[k]?,
            across: sides[(k + 1) % 3],
        })
    }

    /// The fans standing, by group.
    fn fans_by_group(&self, triangles: &[[u32; 3]], adjacency: &Adjacency) -> FastMap<usize, Fan> {
        self.fans
            .iter()
            .filter_map(|&t| Some((self.of[t], self.fan_at(t, triangles, adjacency)?)))
            .collect()
    }
}

/// How many smooth regions were refused a patch, by why.
#[derive(Debug, Clone, Copy, Default)]
struct PatchRefusals {
    not_disk: usize,
    narrow: usize,
    unverified: usize,
}
