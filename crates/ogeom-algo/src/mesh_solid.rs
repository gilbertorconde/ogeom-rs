//! A B-rep solid from a triangle mesh.
//!
//! A mesh already says which triangles meet along which edges: once its
//! repeated vertices are welded, two triangles that share two vertices
//! share an edge. So the topology is built from that connectivity
//! directly — one vertex per mesh vertex, one edge per mesh edge — and no
//! geometric sewing search is needed, which is what lets a printable file
//! of hundreds of thousands of triangles convert in seconds.
//!
//! Adjacent triangles that lie in one plane merge into one planar face,
//! bounded by the region's outer loop and its holes, and a run of boundary
//! segments along one straight line between the same two faces becomes
//! one edge: a cube's twelve triangles become six faces and twelve edges,
//! and the faces take fillets and chamfers as a modelled box's do. A
//! curved region stays faceted; recognising cylinders and the like is a
//! separate step.
//!
//! Windings are made consistent across each connected piece and turned
//! outward. A mesh that does not close — a hole, an edge three triangles
//! share, a piece that cannot be oriented — does not fail: it comes back
//! as open shells, with the report saying why.

use std::collections::HashMap;

use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, LineCurve, PlanarCurve, PlaneSurface};
use ogeom_math::{Direction, Frame, Plane, Point, Point2, Vector};
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
}

impl Default for MeshSolidOptions {
    fn default() -> Self {
        Self {
            merge_coplanar: true,
            coplanar_angle: 1e-3,
            coplanar_distance: None,
            weld: None,
        }
    }
}

/// What [`solid_from_mesh`] did, and where the mesh does not close.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MeshSolidReport {
    /// Triangles built from, after the dropped ones.
    pub triangles: usize,
    /// Faces built: the triangles, or their coplanar groups.
    pub faces: usize,
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
    /// A solid when the mesh closes — a compound of solids when it is
    /// several disjoint closed pieces — and otherwise a shell, or a
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
    let groups = if options.merge_coplanar {
        coplanar_groups(
            &points,
            &triangles,
            &adjacency,
            options.coplanar_angle,
            flat,
            tol,
        )?
    } else {
        Groups::one_each(&points, &triangles, tol)?
    };
    report.faces = groups.planes.len();

    model.begin_operation();
    let built = Builder {
        model,
        points: &points,
        triangles: &triangles,
        adjacency: &adjacency,
        groups: &groups,
        merge: options.merge_coplanar,
        flat,
        tol,
    }
    .build()?;

    // Faces into one shell per piece; closed pieces into solids, a piece
    // nested an odd number of times deep being a void of the one around it.
    let mut shells = Vec::with_capacity(pieces.len());
    for piece in &pieces {
        let mut faces: Vec<Shape> = Vec::new();
        let mut taken = vec![false; groups.planes.len()];
        for &t in &piece.triangles {
            let g = groups.of[t as usize];
            if !taken[g] {
                taken[g] = true;
                faces.push(built[g].clone());
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

/// Triangles gathered into faces: which face each triangle is in, and each
/// face's plane.
struct Groups {
    of: Vec<usize>,
    planes: Vec<Plane>,
}

impl Groups {
    fn one_each(points: &[Point], triangles: &[[u32; 3]], tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self {
            of: (0..triangles.len()).collect(),
            planes: triangles
                .iter()
                .map(|t| plane_of(points, *t, tol))
                .collect::<OgeomResult<_>>()?,
        })
    }
}

/// The plane of a triangle, its normal by the winding, its `x` axis along
/// the first side. The triangle has area, so both are defined.
fn plane_of(points: &[Point], [a, b, c]: [u32; 3], tol: Tolerances) -> OgeomResult<Plane> {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    // Scaled to unit length before the kernel's own normalization, which
    // refuses a vector shorter than the confusion distance — as the cross
    // product of a small triangle's sides may be.
    let normal = (b - a).cross(c - a);
    let z = Direction::new(normal / normal.magnitude(), tol)?;
    let x = Direction::new((b - a) / (b - a).magnitude(), tol)?;
    Ok(Plane::new(Frame::new(a, z, x, tol)?))
}

/// Grow faces across shared edges from the largest triangles down, taking
/// a neighbour whose normal is within the angle of the face's and whose
/// corners all lie within the distance of its plane. Measured against the
/// face's own plane, not the last triangle's, so a gently curved surface
/// does not drift into one face.
fn coplanar_groups(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    angle: f64,
    flat: f64,
    tol: Tolerances,
) -> OgeomResult<Groups> {
    let area = |t: usize| {
        let [a, b, c] = triangles[t].map(|i| points[i as usize]);
        (b - a).cross(c - a).magnitude()
    };
    let mut order: Vec<usize> = (0..triangles.len()).collect();
    order.sort_by(|&x, &y| area(y).total_cmp(&area(x)));
    let cos = angle.cos();
    let mut of = vec![usize::MAX; triangles.len()];
    let mut planes = Vec::new();
    let mut stack = Vec::new();
    for seed in order {
        if of[seed] != usize::MAX {
            continue;
        }
        let g = planes.len();
        let plane = plane_of(points, triangles[seed], tol)?;
        let (origin, normal) = (plane.frame().origin(), plane.frame().z().vector());
        planes.push(plane);
        of[seed] = g;
        stack.push(seed);
        while let Some(t) = stack.pop() {
            for h in 3 * t..3 * t + 3 {
                let Some(twin) = adjacency.twin[h] else {
                    continue;
                };
                let other = twin / 3;
                if of[other] != usize::MAX {
                    continue;
                }
                let tri = triangles[other];
                let [a, b, c] = tri.map(|i| points[i as usize]);
                let n = (b - a).cross(c - a);
                if n.dot(normal) < cos * n.magnitude() {
                    continue;
                }
                if [a, b, c]
                    .iter()
                    .all(|p| (*p - origin).dot(normal).abs() <= flat)
                {
                    of[other] = g;
                    stack.push(other);
                }
            }
        }
    }
    Ok(Groups { of, planes })
}

/// Builds the vertices, edges and faces.
struct Builder<'a> {
    model: &'a mut Model,
    points: &'a [Point],
    triangles: &'a [[u32; 3]],
    adjacency: &'a Adjacency,
    groups: &'a Groups,
    merge: bool,
    flat: f64,
    tol: Tolerances,
}

/// A boundary edge as built: from one kept vertex to another along a run
/// of mesh edges.
struct Run {
    from: u32,
    to: u32,
    /// The points in between, for the tolerance.
    through: Vec<u32>,
    /// The faces either side.
    faces: Vec<usize>,
}

impl Builder<'_> {
    /// Whether a half-edge bounds its face: nothing across it, or another
    /// face.
    fn border(&self, h: Half) -> bool {
        match self.adjacency.twin[h] {
            None => true,
            Some(g) => self.groups.of[g / 3] != self.groups.of[h / 3],
        }
    }

    /// The faces, one per group, in group order.
    fn build(mut self) -> OgeomResult<Vec<Shape>> {
        let halves = self.triangles.len() * 3;
        // The boundary's mesh edges, each once, by its lower vertex first,
        // with the faces on either side.
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
        // A vertex between two boundary edges on one line, with the same
        // faces either side, is dropped: the two become one edge.
        let removable = |v: u32| -> bool {
            if !self.merge {
                return false;
            }
            let Some(list) = incident.get(&v) else {
                return false;
            };
            let [e1, e2] = list[..] else {
                return false;
            };
            if edge_faces[&e1] != edge_faces[&e2] {
                return false;
            }
            let far = |(a, b): (u32, u32)| if a == v { b } else { a };
            let (p, q) = (self.points[far(e1) as usize], self.points[far(e2) as usize]);
            let at = self.points[v as usize];
            (at - p).dot(q - at) > 0.0 && distance_to_line(at, p, q) <= self.flat
        };

        // Runs between kept vertices.
        let mut run_of: HashMap<(u32, u32), (usize, bool)> = HashMap::new();
        let mut runs: Vec<Run> = Vec::new();
        let mut keys: Vec<(u32, u32)> = edge_faces.keys().copied().collect();
        keys.sort_unstable();
        let is_kept: HashMap<u32, bool> = incident.keys().map(|&v| (v, !removable(v))).collect();
        for pass in 0..2 {
            for &key in &keys {
                if run_of.contains_key(&key) {
                    continue;
                }
                let (a, b) = key;
                // Runs start at a kept vertex; a loop of removable ones —
                // which a straight line cannot close — is cut at its first.
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
                let end = at;
                // Held to the straight line from end to end as a whole,
                // not only vertex by vertex, or a slow curve would pass.
                let (p, q) = (self.points[start as usize], self.points[end as usize]);
                let straight = start != end
                    && chain[1..chain.len() - 1]
                        .iter()
                        .all(|&v| distance_to_line(self.points[v as usize], p, q) <= self.flat);
                if straight {
                    let index = runs.len();
                    for (i, e) in edges.iter().enumerate() {
                        run_of.insert(*e, (index, chain[i] == e.0));
                    }
                    runs.push(Run {
                        from: start,
                        to: end,
                        through: chain[1..chain.len() - 1].to_vec(),
                        faces: edge_faces[&key].clone(),
                    });
                } else {
                    for e in edges {
                        run_of.insert(e, (runs.len(), true));
                        runs.push(Run {
                            from: e.0,
                            to: e.1,
                            through: Vec::new(),
                            faces: edge_faces[&e].clone(),
                        });
                    }
                }
            }
        }

        // Each vertex that ends a run stands off the planes of the faces
        // around it by as much as it does; each edge takes the most any of
        // its points stands off its faces' planes or its line.
        let off_plane = |p: Point, g: usize| {
            let frame = self.groups.planes[g].frame();
            (p - frame.origin()).dot(frame.z().vector()).abs()
        };
        let confusion = self.tol.confusion();
        let mut vertices: HashMap<u32, Shape> = HashMap::new();
        let mut edges: Vec<Shape> = Vec::with_capacity(runs.len());
        for run in &runs {
            let (p, q) = (self.points[run.from as usize], self.points[run.to as usize]);
            let mut reach = confusion;
            for &v in [run.from, run.to].iter().chain(&run.through) {
                let at = self.points[v as usize];
                reach = reach.max(distance_to_line(at, p, q));
                for &g in &run.faces {
                    reach = reach.max(off_plane(at, g));
                }
            }
            for v in [run.from, run.to] {
                vertices.entry(v).or_insert_with(|| {
                    self.model
                        .add_vertex(VertexData::new(self.points[v as usize]))
                });
            }
            let curve: Curve = LineCurve::segment(p, q, self.tol)?.into();
            let id = self.model.geometry_mut().add_curve(curve);
            let mut data = EdgeData::on_curve(id, Location::identity(), (0.0, p.distance(q)));
            data.tolerance = Tolerance::new(reach)?;
            let edge = self.model.add_edge(
                data,
                &[vertices[&run.from].clone(), vertices[&run.to].clone()],
            )?;
            edges.push(edge);
        }

        // Each face's loops, walked with the face on the left: from a
        // boundary half-edge to the next one at its end, turning through the
        // face's own triangles around the vertex, which keeps a loop that
        // touches itself at a vertex on its own side.
        let group_count = self.groups.planes.len();
        let mut loops: Vec<Vec<Vec<Half>>> = vec![Vec::new(); group_count];
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
            loops[self.groups.of[h / 3]].push(ring);
        }

        let mut faces = Vec::with_capacity(group_count);
        for (g, rings) in loops.iter().enumerate() {
            let plane = self.groups.planes[g];
            let surface = self
                .model
                .geometry_mut()
                .add_surface(PlaneSurface::new(plane).into());
            let local = |p: Point| {
                let l = plane.frame().to_local(p);
                Point2::new(l.x, l.y)
            };
            // Outer loop first: the one enclosing positive area about the
            // face's normal.
            let mut wires: Vec<(f64, Shape)> = Vec::with_capacity(rings.len());
            for ring in rings {
                let mut entries: Vec<(usize, bool)> = Vec::new();
                for &h in ring {
                    let (a, b) = from_to(self.triangles, h);
                    let (run, along) = run_of[&(a.min(b), a.max(b))];
                    // Along the run if the mesh edge's own direction
                    // agrees with the run's and this half-edge runs that way.
                    let forward = along == (a < b);
                    if entries.last() != Some(&(run, forward)) {
                        entries.push((run, forward));
                    }
                }
                if entries.len() > 1 && entries.first() == entries.last() {
                    entries.pop();
                }
                let mut area = 0.0;
                let mut ring_edges = Vec::with_capacity(entries.len());
                for &(run, forward) in &entries {
                    let r = &runs[run];
                    let (from, to) = if forward {
                        (r.from, r.to)
                    } else {
                        (r.to, r.from)
                    };
                    let (a, b) = (
                        local(self.points[from as usize]),
                        local(self.points[to as usize]),
                    );
                    area += a.x * b.y - b.x * a.y;
                    self.attach(&edges[run], r, surface, local)?;
                    ring_edges.push(if forward {
                        edges[run].clone()
                    } else {
                        edges[run].reversed()
                    });
                }
                wires.push((area, self.model.add_wire(&ring_edges)?));
            }
            wires.sort_by(|a, b| b.0.total_cmp(&a.0));
            let wires: Vec<Shape> = wires.into_iter().map(|(_, w)| w).collect();
            faces.push(
                self.model
                    .add_face(FaceData::new(surface, Location::identity()), &wires)?,
            );
        }
        Ok(faces)
    }

    /// The edge's line in the face's plane, once per face.
    fn attach(
        &mut self,
        edge: &Shape,
        run: &Run,
        surface: ogeom_topo::SurfaceId,
        local: impl Fn(Point) -> Point2,
    ) -> OgeomResult<()> {
        let Some(node) = self.model.node(edge) else {
            ogeom_bail!(Dangling, "an edge just built is not in this model");
        };
        if node.data().as_edge().is_some_and(|d| {
            d.representations
                .iter()
                .any(|rep| matches!(rep, EdgeRepr::PCurve { surface: s, .. } if *s == surface))
        }) {
            return Ok(());
        }
        let (a, b) = (
            local(self.points[run.from as usize]),
            local(self.points[run.to as usize]),
        );
        let pcurve: PlanarCurve = ogeom_geom::Line2d::segment(a, b, self.tol)?.into();
        crate::build::attach_pcurve(
            self.model,
            edge,
            pcurve,
            surface,
            Location::identity(),
            (0.0, a.distance(b)),
        )
    }
}

fn distance_to_line(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let m = d.magnitude();
    if m == 0.0 {
        return p.distance(a);
    }
    (p - a).cross(d).magnitude() / m
}
