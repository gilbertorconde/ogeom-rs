//! A B-rep solid from a triangle mesh.
//!
//! A mesh already says which triangles meet along which edges: once its
//! repeated vertices are welded, two triangles that share two vertices
//! share an edge. So the topology is built from that connectivity
//! directly (one vertex per mesh vertex, one edge per mesh edge), and no
//! geometric sewing search is needed, which is what lets a printable file
//! of hundreds of thousands of triangles convert in seconds.
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
    /// `None` takes the resolution a mesh stored in single precision has:
    /// a millionth of its bounding box's diagonal, never below the weld
    /// distance. The vertices and edges of a merged face widen their
    /// tolerances to cover what they stand off it.
    pub coplanar_distance: Option<f64>,
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
}

impl Default for MeshSolidOptions {
    fn default() -> Self {
        Self {
            merge_coplanar: true,
            coplanar_angle: 1e-3,
            coplanar_distance: None,
            weld: None,
            recognize: true,
            crease: core::f64::consts::FRAC_PI_6,
        }
    }
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
    /// What was built, and where the mesh does not close.
    pub report: MeshSolidReport,
}

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
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|v| remap[v as usize]);
        if a == b || b == c || c == a || !has_area(&points, [a, b, c], weld) {
            report.degenerate_dropped += 1;
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
    let adjacency = Adjacency::new(&triangles);
    report.shells = pieces.len();

    let diagonal = diagonal(&points);
    let flat = options
        .coplanar_distance
        .unwrap_or(1e-6 * diagonal)
        .max(weld);
    let mut groups = segment(&points, &triangles, &adjacency, options, flat, tol)?;
    // Plan until every curved face's boundary is exact, faceting the ones
    // whose boundary is not.
    let mut pinned: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let plan = loop {
        let planner = Planner {
            points: &points,
            triangles: &triangles,
            adjacency: &adjacency,
            groups: &groups,
            merge: options.merge_coplanar,
            pinned: &pinned,
            flat,
            tol,
        };
        let planned = planner.plan()?;

        match planned {
            Ok(plan) => break plan,
            Err(Replan::Pin(vertices)) => pinned.extend(vertices),
            Err(Replan::Facet(failed)) => {
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
                    options.coplanar_angle,
                    flat,
                    &mut groups,
                    tol,
                )?;
            }
        }
    };

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
    report.faces = built.iter().flatten().count();
    report.curved_faces = groups
        .carriers
        .iter()
        .zip(&built)
        .filter(|(c, b)| matches!(c, Carrier::Curved(_)) && b.is_some())
        .count();

    // Faces into one shell per piece; closed pieces into solids, a piece
    // nested an odd number of times deep being a void of the one around it.
    let mut shells = Vec::with_capacity(pieces.len());
    for piece in &pieces {
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
    let shape = if all_closed {
        let mut solids = Vec::new();
        for (i, shell) in shells.iter().enumerate() {
            if depth[i] % 2 == 1 {
                continue;
            }
            let mut members = vec![shell.clone()];
            for (j, void) in shells.iter().enumerate() {
                if depth[j] == depth[i] + 1 && inside(&points, &triangles, &pieces[i], &pieces[j]) {
                    members.push(void.clone());
                }
            }
            solids.push(model.add_solid(&members)?);
        }
        if solids.len() == 1 {
            solids.swap_remove(0)
        } else {
            model.add_compound(&solids)?
        }
    } else if shells.len() == 1 {
        shells.swap_remove(0)
    } else {
        model.add_compound(&shells)?
    };
    Ok(MeshSolid {
        shape,
        closed: all_closed,
        report,
    })
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
                _ => used_more += 1,
            }
            i = j;
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
fn inside(points: &[Point], triangles: &[[u32; 3]], outer: &Piece, inner: &Piece) -> bool {
    let Some(&first) = inner.triangles.first() else {
        return false;
    };
    let origin = points[triangles[first as usize][0] as usize];
    // An off-axis direction, so a ray along a mesh's grid lines is unlikely.
    let direction = Vector::new(0.577_215_664_9, 0.618_033_988_7, 0.533_751_168_7);
    let mut crossings = 0;
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
        if !(0.0..=1.0).contains(&u) {
            continue;
        }
        let q = s.cross(e1);
        let v = direction.dot(q) / det;
        if v < 0.0 || u + v > 1.0 {
            continue;
        }
        if e2.dot(q) / det > 0.0 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
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
struct Groups {
    of: Vec<usize>,
    carriers: Vec<Carrier>,
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
/// face.
fn coplanar_groups(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    angle: f64,
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
    let cos = angle.cos();
    let mut stack = Vec::new();
    for seed in order {
        if groups.of[seed] != usize::MAX {
            continue;
        }
        let g = groups.carriers.len();
        let plane = plane_of(points, triangles[seed], tol)?;
        let (origin, normal) = (plane.frame().origin(), plane.frame().z().vector());
        groups.carriers.push(Carrier::Plane(plane));
        groups.of[seed] = g;
        stack.push(seed);
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
                if n.dot(normal) < cos * n.magnitude() {
                    continue;
                }
                if [a, b, c]
                    .iter()
                    .all(|p| (*p - origin).dot(normal).abs() <= flat)
                {
                    groups.of[other] = g;
                    stack.push(other);
                }
            }
        }
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
    }
}

/// The axis of a surface of revolution: its frame.
fn axis_frame(shape: &Canonical) -> Option<Frame> {
    match shape {
        Canonical::Cylinder(c) => Some(c.frame()),
        Canonical::Cone(c) => Some(c.frame()),
        Canonical::Torus(t) => Some(t.frame()),
        Canonical::Sphere(s) => Some(s.frame()),
        Canonical::Plane(_) => None,
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
        Canonical::Plane(_) => return None,
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
    }
}

/// Which chart directions are angles, and so wrap.
fn periodic(shape: &Canonical) -> (bool, bool) {
    match shape {
        Canonical::Plane(_) => (false, false),
        Canonical::Cylinder(_) | Canonical::Cone(_) | Canonical::Sphere(_) => (true, false),
        Canonical::Torus(_) => (true, true),
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
        sphere_axes(points, triangles, adjacency, &mut groups, flat, tol);
        align_axes(points, &mut groups, flat, tol);
        hole_frames(points, triangles, adjacency, &mut groups, tol);
        slit_bands(points, triangles, adjacency, &mut groups, tol);
    }
    coplanar_groups(
        points,
        triangles,
        adjacency,
        options.coplanar_angle,
        flat,
        &mut groups,
        tol,
    )?;
    Ok(groups)
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
                    queue.push_back(other);
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
    let mesh = Surfaces {
        points,
        triangles,
        adjacency,
        normals: triangles.iter().map(|t| unit_normal(points, *t)).collect(),
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
    let order: Vec<usize> = (0..n).map(|i| (i * stride) % n).collect();
    let mut next = 0;
    while next < n {
        batch += 1;
        let mut seeds = Vec::with_capacity(batch_size);
        while next < n && seeds.len() < batch_size {
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
                let retired = if wide { 1 } else { first_sample };
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
                        mine.insert(other);
                        region.push(other);
                        for &v in &triangles[other] {
                            if seen.insert(v) {
                                vertices.push(v);
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
            let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
            let deviation = worst_deviation(&shape, &pts);
            let flat_too = crate::recognize::is_flat(&pts, flat, tol);
            if deviation > flat || flat_too || region.len() < 2 {
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
                centre: (0.0, 0.0),
                wraps: false,
                wraps_v: false,
                fixed: false,
                vertices,
            }));
        }
    }
}

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
        let shape = match curved.shape {
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
        if let Carrier::Curved(curved) = &mut groups.carriers[g] {
            curved.shape = shape;
            curved.fixed = true;
            curved.centre = (core::f64::consts::PI, curved.centre.1);
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

/// Put coaxial surfaces on one axis and one angular origin, so bands that
/// meet along a circle meet at their seams too; then fix each curved
/// region's chart branch and whether it wraps.
fn align_axes(points: &[Point], groups: &mut Groups, flat: f64, tol: Tolerances) {
    let mut leaders: Vec<Frame> = Vec::new();
    for carrier in &mut groups.carriers {
        let Carrier::Curved(curved) = carrier else {
            continue;
        };
        let Some(frame) = axis_frame(&curved.shape) else {
            continue;
        };
        let sphere = matches!(curved.shape, Canonical::Sphere(_));
        let lead = leaders.iter().find(|l| {
            let parallel = l.z().vector().cross(frame.z().vector()).magnitude() <= 1e-3;
            let w = frame.origin() - l.origin();
            let off = (w - l.z().vector() * w.dot(l.z().vector())).magnitude();
            !sphere && parallel && off <= flat * 10.0
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
            if let Ok(snapped) = Frame::new(origin, axis, lead.x(), tol)
                && let Some(shape) = on_frame(&curved.shape, snapped, tol)
            {
                let pts: Vec<Point> = curved
                    .vertices
                    .iter()
                    .map(|&v| points[v as usize])
                    .collect();
                let deviation = worst_deviation(&shape, &pts);
                if deviation <= flat {
                    curved.shape = shape;
                    curved.deviation = deviation;
                }
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
    flat: f64,
    tol: Tolerances,
}

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
        // on one line; along a curved face always, the curve deciding.
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
            if any_curved(faces) {
                return true;
            }
            let far = |(a, b): (u32, u32)| if a == v { b } else { a };
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
        };
        let mut failed: Vec<usize> = Vec::new();
        let mut keys: Vec<(u32, u32)> = edge_faces.keys().copied().collect();
        keys.sort_unstable();
        for pass in 0..2 {
            for &key in &keys {
                if plan.edge_of.contains_key(&key) {
                    continue;
                }
                let (a, b) = key;
                let start = if is_kept[&a] {
                    a
                } else if is_kept[&b] {
                    b
                } else if pass == 1 {
                    a
                } else {
                    continue;
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
                let faces = edge_faces[&key].clone();
                if any_curved(&faces) {
                    match self.snapped(&chain, &faces) {
                        Some(spec) => {
                            let index = plan.edges.len();
                            let (spec, forward, images) = spec;
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
                    let clear = |values: &mut dyn Iterator<Item = f64>| {
                        let (lo, hi) = values
                            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| {
                                (a.min(x), b.max(x))
                            });
                        (lo / tau).floor() == (hi / tau).floor()
                    };
                    clear(&mut polygon.iter().map(|p| p.0))
                        && (!wraps_v || clear(&mut polygon.iter().map(|p| p.1)))
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
            } else if rings.len() == 2 && (!sphere || curved.fixed) {
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
    /// that face itself so its pcurve there is exact. `None` when no such
    /// curve holds the chain and the faces both.
    fn snapped(&self, chain: &[u32], faces: &[usize]) -> Option<(Snapped, bool, Images)> {
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
        let mut on: Vec<Point> = pts.to_vec();
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
        let (curve, samples): (Curve, Vec<f64>) = if on.len() == 2 {
            let length = on[0].distance(on[1]);
            let line: Curve = LineCurve::segment(on[0], on[1], self.tol).ok()?.into();
            (
                line,
                (0..=16).map(|k| length * f64::from(k) / 16.0).collect(),
            )
        } else {
            // A loop is carried on past its ends and cut back, as a fitted
            // section is.
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
            let mut spline = crate::fit::interpolate_at(&padded, &parameters, 3, self.tol).ok()?;
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
        for t in samples {
            let p = curve.point_at(t.clamp(range.0, range.1), self.tol).ok()?;
            tolerance = tolerance.max(fa(p).abs()).max(fb(p).abs());
        }
        if tolerance > reach.max(longest * CHORD_SAG) {
            return None;
        }
        Some((
            if closed {
                Snapped::Loop(curve, range, tolerance)
            } else {
                Snapped::Open(curve, range, tolerance)
            },
            true,
            Vec::new(),
        ))
    }

    /// The curve two faces meet along, where it is no parallel or ruling of
    /// either: vertices of the chain, and points between them, solved onto
    /// both surfaces, and the curve interpolated through them. Its image on
    /// each curved face is interpolated through the same points' chart
    /// positions at the same parameters, so the two run together. `None`
    /// where the surfaces meet tangentially, the solve does not settle, or
    /// the curve strays from either face past the reach.
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
        let [a, b] = faces[..] else {
            return None;
        };
        let (fa, fb) = (self.signed(a)?, self.signed(b)?);
        let steps = if closed { pts.len() } else { pts.len() - 1 };
        // A long chain is taken a few vertices at a time: the curve needs
        // its shape, not every vertex.
        let stride = steps.div_ceil(SECTION_SPANS).max(1);
        let split = count.div_ceil(steps.div_ceil(stride)).max(SECTION_SPLIT);
        let mut on = Vec::new();
        let mut i = 0;
        while i < steps {
            let next = (i + stride).min(steps);
            let (p, q) = (pts[i], pts[next % pts.len()]);
            for k in 0..split {
                #[allow(clippy::cast_precision_loss, reason = "a handful of splits")]
                let f = k as f64 / split as f64;
                let limit = if k == 0 { reach } else { p.distance(q) + reach };
                on.push(onto_both(&fa, &fb, p + (q - p) * f, reach, limit)?);
            }
            i = next;
        }
        on.push(if closed {
            on[0]
        } else {
            onto_both(&fa, &fb, pts[pts.len() - 1], reach, reach)?
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
                tolerance = tolerance.max(fa(p).abs()).max(fb(p).abs());
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
                let shape = c.shape;
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
    for _ in 0..40 {
        let (va, vb) = (fa(p), fb(p));
        if va.abs().max(vb.abs()) <= 1e-13 * (1.0 + p.to_vector().magnitude()) {
            break;
        }
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
                if (big + off - mean_r).abs() <= (big - off - mean_r).abs() {
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
    Ok(match curved.shape {
        Canonical::Cylinder(c) => CylinderSurface::new(c, heights())?.into(),
        Canonical::Cone(c) => {
            // Short of the apex, where the cone's radius runs out.
            let (lo, hi) = heights();
            let apex = -c.reference_radius() / c.half_angle().tan();
            let lo = lo.max(apex + (hi - apex) * 1e-6);
            ConeSurface::new(c, (lo, hi))?.into()
        }
        Canonical::Sphere(s) => SphereSurface::new(s).into(),
        Canonical::Torus(t) => TorusSurface::new(t).into(),
        Canonical::Plane(p) => PlaneSurface::new(p).into(),
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
                    Some(self.whole_face(curved, g)?)
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

    /// A whole sphere or torus, a piece with no boundary: the face the
    /// primitive builds on the recognized surface.
    /// A sphere or torus whole but for holes: the whole surface's face, as
    /// the primitive builds it, with each ring an inner wire. The rings
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
        let whole = self.whole_face(curved, g)?;
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

    fn whole_face(&mut self, curved: &Curved, g: usize) -> OgeomResult<Shape> {
        let outward = self.outward(curved, g);
        let built = match curved.shape {
            Canonical::Sphere(s) => {
                crate::primitive::make_sphere(self.model, s.frame(), s.radius(), self.tol)?
            }
            Canonical::Torus(t) => crate::primitive::make_torus(
                self.model,
                t.frame(),
                t.major_radius(),
                t.minor_radius(),
                self.tol,
            )?,
            _ => ogeom_bail!(Construction, "only a sphere or a torus is whole"),
        };
        let Some(face) =
            ogeom_topo::explore_unique(self.model, &built.shape, ogeom_topo::ShapeType::Face)?
                .into_iter()
                .next()
        else {
            ogeom_bail!(Construction, "a primitive came back with no face");
        };
        if let Some(node) = self.model.node_mut(&face)
            && let ogeom_topo::NodeData::Face(data) = node.data_mut()
        {
            data.tolerance = Tolerance::new(curved.deviation.max(self.tol.confusion()))?;
        }
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
                let u = match out.last() {
                    Some(&(last, _)) => last + ogeom_math::elementary::wrap_signed_angle(u - last),
                    None => u,
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

fn distance_to_line(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let m = d.magnitude();
    if m == 0.0 {
        return p.distance(a);
    }
    (p - a).cross(d).magnitude() / m
}
