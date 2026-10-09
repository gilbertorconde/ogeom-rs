//! The mesh made ready to segment: its vertices welded, the half-edges
//! across each shared edge paired, folds unfolded, and each connected piece
//! wound consistently, a closed one outward, with which pieces lie inside
//! which.

use ogeom_core::FastMap;
use ogeom_math::{Point, Vector};

use super::MeshSolidReport;

/// Weld points on a grid of the weld distance, looking in the neighbouring
/// cells too, so two points either side of a cell wall still meet.
pub(super) fn weld_points(positions: &[Point], weld: f64) -> (Vec<Point>, Vec<u32>) {
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
    let mut grid: FastMap<(i64, i64, i64), Vec<u32>> =
        FastMap::with_capacity_and_hasher(positions.len(), Default::default());
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
pub(super) fn has_area(points: &[Point], [a, b, c]: [u32; 3], weld: f64) -> bool {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
    let twice_area = (b - a).cross(c - a).magnitude();
    let longest = a.distance(b).max(b.distance(c)).max(c.distance(a));
    twice_area > weld * longest
}

pub(super) fn diagonal(points: &[Point]) -> f64 {
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
pub(super) type Half = usize;

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
pub(super) fn unfold(points: &[Point], triangles: &mut [[u32; 3]], adjacency: &Adjacency) -> usize {
    let normal = |t: [u32; 3]| -> Vector {
        let [a, b, c] = t.map(|v| points[v as usize]);
        (b - a).cross(c - a)
    };
    let unit = |v: Vector| {
        let m = v.magnitude();
        if m > 0.0 { v / m } else { v }
    };
    let mut edges: ogeom_core::FastSet<(u32, u32)> = triangles
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

pub(super) fn from_to(triangles: &[[u32; 3]], h: Half) -> (u32, u32) {
    let t = triangles[h / 3];
    (t[h % 3], t[(h % 3 + 1) % 3])
}

pub(super) fn next(h: Half) -> Half {
    h - h % 3 + (h % 3 + 1) % 3
}

/// Which half-edges pair across a mesh edge exactly two triangles share.
#[derive(Clone)]
pub(super) struct Adjacency {
    /// The other triangle's half-edge on the same mesh edge, where exactly
    /// one other triangle uses it.
    pub(super) twin: Vec<Option<Half>>,
    pub(super) used_once: usize,
    pub(super) used_more: usize,
}

impl Adjacency {
    pub(super) fn new(triangles: &[[u32; 3]]) -> Self {
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
#[derive(Clone)]
pub(super) struct Piece {
    pub(super) triangles: Vec<u32>,
    pub(super) closed: bool,
}

/// Make windings agree across every shared edge, piece by piece; turn a
/// closed piece outward and leave an open one the way most of its
/// triangles already were.
pub(super) fn orient(
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
            // About a vertex of the piece: the sum is the same about any
            // point, and about a far one it is a difference of huge
            // products whose sign is rounding.
            let apex = points[triangles[seed][0] as usize];
            let volume: f64 = members
                .iter()
                .map(|&t| {
                    let mut tri = triangles[t];
                    if flip[t] == Some(true) {
                        tri.swap(1, 2);
                    }
                    signed_volume(points, tri, apex)
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

/// The signed volume of the tetrahedron a triangle forms with `apex`.
fn signed_volume(points: &[Point], [a, b, c]: [u32; 3], apex: Point) -> f64 {
    let [a, b, c] = [a, b, c].map(|i| points[i as usize] - apex);
    a.dot(b.cross(c)) / 6.0
}

/// Whether piece `inner` lies inside closed piece `outer`: a ray from one
/// of its vertices crosses `outer` an odd number of times.
///
/// A ray through a shared edge or vertex of `outer` meets two triangles
/// there, or none, and its count says nothing; so a ray passing that close
/// to any triangle's boundary is set aside and the next direction tried.
/// Where every direction grazes something, the parities they read vote.
pub(super) fn inside(
    points: &[Point],
    triangles: &[[u32; 3]],
    outer: &Piece,
    inner: &Piece,
) -> bool {
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
