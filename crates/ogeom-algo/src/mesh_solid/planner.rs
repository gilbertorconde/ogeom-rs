//! The plan of the build, decided before anything is built: the edges with
//! their ends and tolerances, each face's loops, each curved face's surface
//! and boundary layout in its chart, the seams of faces round an axis or
//! threaded through their holes, and the fans and wedges between facets
//! and curved faces. A plan that cannot be made asks for vertices to keep,
//! facets to build as fans, or curved faces to facet.

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, LineCurve, PlanarCurve};
use ogeom_math::{Direction, Frame, Point, Point2, Vector};

use super::caches::{ImageCache, SnapCache};
use super::seams::{
    cap_rings, cap_seam, centred_rings, choose_seam, distance_to_line, free_angle, hole_polygons,
    round_the_tube, segments_cross, swapped, swapped_rings, thread_pieces, torus_seam_v,
    turns_along, windings,
};
use super::segment::{
    chart, cone_rim_start, cone_seamed_through, evaluate, gradient, periodic, unwrapped,
    walk_entries,
};
use super::snap::{Images, Snapped, surface_of};
use super::weld::{Adjacency, Half, from_to, next};
use super::{Carrier, Curved, Fan, Groups, Layout, Stop, Thread};
use crate::recognize::Canonical;

/// A vertex as built: a mesh vertex, or a point the construction placed,
/// where a closed circle is bounded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Corner {
    Mesh(u32),
    Placed(usize),
}

/// An edge as planned.
pub(super) struct EdgeSpec {
    pub(super) curve: Curve,
    pub(super) range: (f64, f64),
    pub(super) ends: [Corner; 2],
    pub(super) tolerance: f64,
    pub(super) closed_circle: bool,
}

/// Everything the build needs, decided before anything is built.
pub(super) struct Plan {
    pub(super) edges: Vec<EdgeSpec>,
    /// Each boundary mesh edge, by its lower vertex first: its planned
    /// edge, and whether walking it from the lower vertex runs the edge
    /// forward.
    pub(super) edge_of: FastMap<(u32, u32), (usize, bool)>,
    /// Each face's loops as half-edges, walked with the face on the left.
    pub(super) loops: Vec<Vec<Vec<Half>>>,
    pub(super) placed: Vec<Point>,
    /// Each curved face's surface, windowed to its boundary.
    pub(super) surfaces: Vec<Option<ogeom_geom::SurfaceGeometry>>,
    /// Each edge's image on each curved face it bounds, by (edge, face),
    /// with how far the image strays from the edge.
    pub(super) pcurves: FastMap<(usize, usize), (PlanarCurve, f64)>,
    /// How each curved face's boundary lies in its chart.
    pub(super) layouts: Vec<Layout>,
    /// The seam chain of each face laid out [`Layout::Threaded`].
    pub(super) threads: FastMap<usize, Thread>,
    /// How many free edges of curved faces are fitted curves.
    pub(super) free_fitted: usize,
}

/// Plans the edges and loops.
pub(super) struct Planner<'a> {
    pub(super) points: &'a [Point],
    pub(super) triangles: &'a [[u32; 3]],
    pub(super) adjacency: &'a Adjacency,
    pub(super) groups: &'a Groups,
    pub(super) merge: bool,
    /// Vertices kept as edge ends whatever lies either side of them.
    pub(super) pinned: &'a ogeom_core::FastSet<u32>,
    /// Mesh edges whose seam is threaded straight through the chain's own
    /// vertices: the solved curve bent the trim of a face beside it back
    /// across itself.
    pub(super) straight: &'a ogeom_core::FastSet<(u32, u32)>,
    /// The turn, in radians, from which a free boundary of a curved face
    /// is cut into separate edges at a vertex.
    pub(super) crease: f64,
    pub(super) flat: f64,
    pub(super) tol: Tolerances,
    /// Curves snapped by earlier plans with the same points and distance.
    pub(super) snaps: &'a std::sync::Mutex<SnapCache>,
    /// Edge images found by earlier plans.
    pub(super) images: &'a std::sync::Mutex<ImageCache>,
    /// The facets built as fans, by group.
    pub(super) fans: &'a FastMap<usize, Fan>,
}

/// A chain of boundary vertices between two kept ones, and its mesh edges.
type Chain = (Vec<u32>, Vec<(u32, u32)>);

/// Why a plan was refused: curved faces to facet, vertices to keep, or
/// facets to build as fans.
pub(super) enum Replan {
    Facet(Vec<usize>),
    Pin(Vec<u32>),
    Fan(Vec<usize>),
}

/// Whether a triangle with its corners on `shape` stands off it at its
/// centroid no farther than a facet of its size over that much turn of the
/// surface would: its longest side times the turn of the surface's normals
/// across its corners, over six. A triangle spanning a recess, its
/// corners on the rim and its middle over the hollow, stands off more.
pub(super) fn sags_as_the_surface(shape: &Canonical, corners: [Point; 3], flat: f64) -> bool {
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
pub(super) fn is_sliver(corners: [Point; 3]) -> bool {
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
pub(super) const CHORD_SAG: f64 = 0.05;

/// How much wider than the gap it was measured from an edge's tolerance is
/// recorded: a millionth.
pub(super) const TOLERANCE_MARGIN: f64 = 1e-6;

/// The most triangles a cluster the region surrounds may hold and still
/// be taken into it whatever its normals.
pub(super) const ENCLOSED_CLUSTER: usize = 16;

/// The least coplanar distance recognition holds to, against the mesh's
/// measured scatter: twice it, for a seam between two surfaces each fitted
/// to scattered vertices, and a quarter more for the tail its ninetieth
/// percentile leaves out.
pub(super) const NOISE_FLOOR: f64 = 2.5;

/// How far a snapped curve may stand off its chain or its faces.
pub(super) const REACH: f64 = 20.0;

impl Planner<'_> {
    fn border(&self, h: Half) -> bool {
        match self.adjacency.twin[h] {
            None => true,
            Some(g) => self.groups.of[g / 3] != self.groups.of[h / 3],
        }
    }

    pub(super) fn curved(&self, g: usize) -> Option<&Curved> {
        match &self.groups.carriers[g] {
            Carrier::Curved(c) => Some(c),
            _ => None,
        }
    }

    /// The plan, or the curved faces whose boundary could not be built
    /// exactly and are to be faceted instead.
    #[allow(clippy::too_many_lines, reason = "one pass over the boundary")]
    pub(super) fn plan(&self) -> OgeomResult<Result<Plan, Replan>> {
        let halves = self.triangles.len() * 3;
        let mut edge_faces: FastMap<(u32, u32), Vec<usize>> = FastMap::default();
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
        let mut incident: FastMap<u32, Vec<(u32, u32)>> = FastMap::default();
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
        let is_kept: FastMap<u32, bool> = incident.keys().map(|&v| (v, !removable(v))).collect();

        let mut plan = Plan {
            edges: Vec::new(),
            edge_of: FastMap::default(),
            loops: vec![Vec::new(); self.groups.carriers.len()],
            placed: Vec::new(),
            surfaces: vec![None; self.groups.carriers.len()],
            pcurves: FastMap::default(),
            layouts: vec![Layout::Open; self.groups.carriers.len()],
            threads: FastMap::default(),
            free_fitted: 0,
        };
        let mut failed: Vec<usize> = Vec::new();
        let mut fan_wanted: Vec<usize> = Vec::new();
        let mut wedge_corners: Vec<u32> = Vec::new();
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
            let mut taken: ogeom_core::FastSet<(u32, u32)> = ogeom_core::FastSet::default();
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
                        && self.fan_seam_of(&chain, faces).is_none()
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
                    let fan = self
                        .fan_seam_of(&chain, &faces)
                        .and_then(|(g, wedge)| self.lifted_seam(&chain, g, wedge));
                    let snapped = if fan.is_some() {
                        fan
                    } else if threaded {
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
                            if let Some(t) = self.fan_for(&chain, &faces) {
                                fan_wanted.push(t);
                            } else if let Some(v) = self.wedge_corner(&chain, &faces) {
                                wedge_corners.push(v);
                            }
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
        // holed or not, gets a seam of its own between them. A torus whose
        // holes leave no seam position free has its seam run through them.
        let mut thread_pins: Vec<u32> = Vec::new();
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
                if (sphere && !curved.fixed) || curved.wraps == curved.wraps_v || !resolved {
                    return None;
                }
                let windings: Option<Vec<i32>> = rings
                    .iter()
                    .map(|ring| turns_along(curved, ring, self.triangles, self.points, self.tol))
                    .collect();
                let windings = windings?;
                let rims = windings.iter().filter(|w| w.abs() == 1).count();
                let holes = windings.iter().filter(|w| **w == 0).count();
                if rims + holes != windings.len() {
                    return None;
                }
                match rims {
                    2 => Some(Layout::Wrapped),
                    1 if sphere && holes > 0 => Some(Layout::HoledCap),
                    _ => None,
                }
            };
            // A sphere's circles all latitudes of its frame.
            let latitudes = || {
                rings.iter().all(|ring| {
                    let (edge, _) = self.entry(&plan, ring[0]);
                    match (&plan.edges[edge].curve, &curved.shape) {
                        (Curve::Circle(c), Canonical::Sphere(s)) => {
                            c.circle()
                                .frame()
                                .z()
                                .vector()
                                .cross(s.frame().z().vector())
                                .magnitude()
                                <= 1e-2
                        }
                        _ => false,
                    }
                })
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
            let holed = !rings.is_empty() && holed();
            let thread = if !holed && !rings.is_empty() && torus && curved.wraps && curved.wraps_v {
                self.seam_thread(&plan, curved, rings)
            } else {
                None
            };
            let layout = if rings.is_empty() {
                (sphere || torus).then_some(Layout::Whole)
            } else if holed {
                Some(Layout::Holed)
            } else if thread.is_some() {
                Some(Layout::Threaded)
            } else if (curved.wraps && curved.wraps_v) || !circles {
                // A drill's point whose rim is arcs of one circle, met by
                // several faces, closes at its apex as one met by one.
                wrapped().or_else(|| {
                    (rings.len() == 1 && self.cone_tip(&plan, curved, &rings[0]))
                        .then_some(Layout::Cap)
                })
            } else if rings.len() == 1
                && ((sphere && curved.fixed) || self.cone_tip(&plan, curved, &rings[0]))
            {
                Some(Layout::Cap)
            } else if rings.len() == 2
                && (!sphere || (curved.fixed && latitudes()))
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
                Some(Layout::HoledCap) => self
                    .cap_seam_clear(&plan, curved, &plan.loops[g])
                    .then_some(Layout::HoledCap),
                other => other,
            };
            if let Some(thread) = thread {
                for stop in &thread.stops {
                    for v in [stop.enter, stop.leave] {
                        if !self.pinned.contains(&v) && !thread_pins.contains(&v) {
                            thread_pins.push(v);
                        }
                    }
                }
                plan.threads.insert(g, thread);
            }
            match layout {
                Some(layout) => plan.layouts[g] = layout,
                None => failed.push(g),
            }
        }
        if !thread_pins.is_empty() && failed.is_empty() {
            return Ok(Err(Replan::Pin(thread_pins)));
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
            let cap = plan.layouts[g] == Layout::Cap;
            let turned;
            let curved = match &curved.shape {
                Canonical::Cone(cone) if cap => {
                    let entries = walk_entries(&plan.loops[g][0], |h| self.entry(&plan, h));
                    let Some(cone) = cone_rim_start(cone, &plan.edges, &entries, reach, self.tol)
                        .and_then(|start| cone_seamed_through(cone, start, self.tol))
                    else {
                        failed.push(g);
                        continue;
                    };
                    turned = Curved {
                        shape: Canonical::Cone(cone),
                        ..curved.clone()
                    };
                    &turned
                }
                _ => curved,
            };
            let Ok(surface) = surface_of(curved, self.points, cap, self.tol) else {
                failed.push(g);
                continue;
            };
            if matches!(
                plan.layouts[g],
                Layout::Open
                    | Layout::Wrapped
                    | Layout::Holed
                    | Layout::HoledCap
                    | Layout::Threaded
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
                        let image = self.image(curved, &surface, spec, reach);
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
        if !fan_wanted.is_empty() {
            Ok(Err(Replan::Fan(fan_wanted)))
        } else if !wedge_corners.is_empty() {
            Ok(Err(Replan::Pin(wedge_corners)))
        } else if failed.is_empty() {
            Ok(Ok(plan))
        } else {
            Ok(Err(Replan::Facet(failed)))
        }
    }

    /// Whether a cone's face closes at its apex inside one rim: the rim one
    /// full circle starting off the axis, or arcs of one parallel of the
    /// cone going once round (see [`cone_rim_start`]), a vertex of the
    /// region at the apex within the reach, and every vertex on the apex's
    /// own nappe. A drill's point is such a face; its surface is the cone
    /// turned about its axis so its seam runs through the rim's start.
    fn cone_tip(&self, plan: &Plan, curved: &Curved, ring: &[Half]) -> bool {
        let Canonical::Cone(cone) = &curved.shape else {
            return false;
        };
        let reach = self.flat * REACH;
        let entries = walk_entries(ring, |h| self.entry(plan, h));
        let Some(start) = cone_rim_start(cone, &plan.edges, &entries, reach, self.tol) else {
            return false;
        };
        if cone_seamed_through(cone, start, self.tol).is_none() {
            return false;
        }
        let frame = cone.frame();
        let axis = frame.z().vector();
        let apex = frame.origin() + axis * (-cone.reference_radius() / cone.half_angle().tan());
        let mut nearest = f64::INFINITY;
        for &v in &curved.vertices {
            let p = self.points[v as usize];
            if (p - apex).dot(axis) < -reach {
                return false;
            }
            nearest = nearest.min(p.distance(apex));
        }
        nearest <= reach
    }

    /// Whether a face round its axis has a seam clear of its holes, turning
    /// a rim that is one full circle so its vertex stands in the widest
    /// stretch the holes leave free where the rims' own vertices offer none.
    fn seat_seam(&self, plan: &mut Plan, curved: &Curved, rings: &[Vec<Half>]) -> bool {
        use ogeom_geom::Curve3d as _;
        let shape = &curved.shape;
        let round_tube = round_the_tube(curved);
        let windings: Vec<i32> = rings
            .iter()
            .map(|ring| {
                turns_along(curved, ring, self.triangles, self.points, self.tol).unwrap_or(0)
            })
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
        let hole_rings = swapped_rings(
            centred_rings(
                curved,
                hole_polygons(shape, &holes, self.triangles, self.points, self.tol),
            ),
            round_tube,
        );
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
                curved,
                round_tube,
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
            let Some(at) = chart(shape, plan.placed[k], self.tol) else {
                continue;
            };
            let (_, across) = swapped(at, round_tube);
            let target = evaluate(shape, swapped((angle, across), round_tube));
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

    /// Whether a sphere's cap with holes has a seam from a vertex of its
    /// rim to the pole clear of the holes (see [`cap_seam`]).
    fn cap_seam_clear(&self, plan: &Plan, curved: &Curved, rings: &[Vec<Half>]) -> bool {
        let Some((rim, _, holes)) = cap_rings(curved, rings, self.triangles, self.points, self.tol)
        else {
            return false;
        };
        let mut entries: Vec<(usize, bool)> = Vec::new();
        for &h in &rings[rim] {
            let entry = self.entry(plan, h);
            if entries.last() != Some(&entry) {
                entries.push(entry);
            }
        }
        if entries.len() > 1 && entries.first() == entries.last() {
            entries.pop();
        }
        let starts: Vec<Point> = entries
            .into_iter()
            .map(
                |(edge, forward)| match plan.edges[edge].ends[usize::from(!forward)] {
                    Corner::Mesh(v) => self.points[v as usize],
                    Corner::Placed(k) => plan.placed[k],
                },
            )
            .collect();
        let holes: Vec<&[Half]> = holes.iter().map(|&k| rings[k].as_slice()).collect();
        cap_seam(
            curved,
            &rings[rim],
            &holes,
            &starts,
            self.triangles,
            self.points,
            self.tol,
        )
        .is_some()
    }

    /// The seam chain of a torus whole but for holes that between them
    /// cross every parallel, or every meridian, while leaving a circle the
    /// other way free (see [`Thread`]). Of the chains at a few levels
    /// across, the one through fewest holes, then turning least across;
    /// `None` where every chain crosses a hole.
    fn seam_thread(&self, plan: &Plan, curved: &Curved, rings: &[Vec<Half>]) -> Option<Thread> {
        use core::f64::consts::{PI, TAU};
        const LEVELS: u32 = 32;
        let resolved = rings
            .iter()
            .all(|ring| ring.iter().all(|&h| self.entry(plan, h).0 != usize::MAX));
        if !resolved {
            return None;
        }
        for ring in rings {
            let turns = windings(&curved.shape, ring, self.triangles, self.points, self.tol);
            if turns != Some((0, 0)) {
                return None;
            }
        }
        let slices: Vec<&[Half]> = rings.iter().map(Vec::as_slice).collect();
        let polygons = centred_rings(
            curved,
            hole_polygons(
                &curved.shape,
                &slices,
                self.triangles,
                self.points,
                self.tol,
            ),
        );
        if polygons.len() != rings.len() || polygons.iter().any(Vec::is_empty) {
            return None;
        }
        let span = |ring: &[(f64, f64)], k: usize| {
            ring.iter()
                .map(|p| if k == 0 { p.0 } else { p.1 })
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| {
                    (a.min(x), b.max(x))
                })
        };
        let mut best: Option<((usize, f64), Thread)> = None;
        for round_tube in [false, true] {
            let (centre_s, centre_t) = swapped(curved.centre, round_tube);
            let side = centre_s - PI;
            // Every hole clear of the free circle, moved to lie on past it.
            let mut work: Vec<Vec<(f64, f64)>> = Vec::with_capacity(polygons.len());
            for ring in swapped_rings(polygons.clone(), round_tube) {
                let (lo, hi) = span(&ring, 0);
                let turn = ((lo - side) / TAU).floor();
                if ((hi - side) / TAU).floor() != turn {
                    work.clear();
                    break;
                }
                work.push(ring.into_iter().map(|(a, b)| (a - turn * TAU, b)).collect());
            }
            if work.len() != polygons.len() {
                continue;
            }
            for level in 0..LEVELS {
                let at = centre_t - PI + (f64::from(level) + 0.5) * TAU / f64::from(LEVELS);
                let mut stops: Vec<Stop> = Vec::new();
                for (k, ring) in work.iter().enumerate() {
                    let (lo, hi) = span(ring, 1);
                    let shift = ((at - lo) / TAU).floor() * TAU;
                    if lo + shift >= at || hi + shift <= at {
                        continue;
                    }
                    let extreme = |further: fn(f64, f64) -> bool| {
                        (0..ring.len())
                            .fold(0, |m, i| if further(ring[i].0, ring[m].0) { i } else { m })
                    };
                    let (first, last) = (extreme(|a, b| a < b), extreme(|a, b| a > b));
                    let vertex = |i: usize| from_to(self.triangles, rings[k][i]).0;
                    stops.push(Stop {
                        ring: k,
                        enter: vertex(first),
                        leave: vertex(last),
                        enter_at: (ring[first].0, ring[first].1 + shift),
                        leave_at: (ring[last].0, ring[last].1 + shift),
                    });
                }
                if stops.is_empty() {
                    continue;
                }
                stops.sort_by(|a, b| a.enter_at.0.total_cmp(&b.enter_at.0));
                if stops.windows(2).any(|w| w[0].leave_at.0 >= w[1].enter_at.0) {
                    continue;
                }
                let (first, last) = (stops[0], stops[stops.len() - 1]);
                let from = (last.leave_at.0 - TAU, last.leave_at.1);
                let f = (side - from.0) / (first.enter_at.0 - from.0);
                let start = (side, from.1 + (first.enter_at.1 - from.1) * f);
                let pieces = thread_pieces(start, &stops);
                let crosses = pieces.iter().any(|&(a, b)| {
                    work.iter().any(|ring| {
                        [-TAU, 0.0, TAU].iter().any(|ds| {
                            [-TAU, 0.0, TAU].iter().any(|dt| {
                                (0..ring.len()).any(|i| {
                                    let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                                    segments_cross(a, b, (p.0 + ds, p.1 + dt), (q.0 + ds, q.1 + dt))
                                })
                            })
                        })
                    })
                });
                if crosses {
                    continue;
                }
                let rise: f64 = pieces.iter().map(|(a, b)| (b.1 - a.1).abs()).sum();
                let score = (stops.len(), rise);
                if best.as_ref().is_none_or(|(held, _)| {
                    score.0 < held.0 || (score.0 == held.0 && score.1 < held.1)
                }) {
                    best = Some((
                        score,
                        Thread {
                            round_tube,
                            start,
                            stops,
                        },
                    ));
                }
            }
            if best.is_some() {
                break;
            }
        }
        best.map(|(_, thread)| thread)
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
    /// The facet a seam that no curve places could be built as a fan
    /// from: one triangle alone in its planar group, its side along the
    /// chain's single span against the curved face across it, its other two
    /// against planar faces, and not a fan already. Its seam is then the
    /// span drawn straight in the curved face's chart, where no curve
    /// through the span's ends keeps to both the facet's plane and the
    /// curved surface (a chord of a neighbour left faceted, cutting across
    /// the curved surface's curvature).
    fn fan_for(&self, chain: &[u32], faces: &[usize]) -> Option<usize> {
        let [a, b] = faces[..] else {
            return None;
        };
        let [p, q] = chain[..] else {
            return None;
        };
        let (facet, curved) = match (self.curved(a), self.curved(b)) {
            (Some(_), None) => (b, a),
            (None, Some(_)) => (a, b),
            _ => return None,
        };
        let mut members = self
            .groups
            .of
            .iter()
            .enumerate()
            .filter(|&(_, &g)| g == facet)
            .map(|(t, _)| t);
        let t = members.next()?;
        if members.next().is_some() || self.groups.fans.contains(&t) {
            return None;
        }
        let fan = self.groups.fan_at(t, self.triangles, self.adjacency)?;
        (fan.seam_between(p, q) == Some(curved)).then_some(t)
    }

    /// The corner to pin where a chain no curve places runs along both
    /// curved sides of a wedge onto one curved face: the two sides meet at
    /// the chain's middle vertex, which split there leaves each side a span
    /// of its own, a seam the wedge can be built from.
    fn wedge_corner(&self, chain: &[u32], faces: &[usize]) -> Option<u32> {
        let [a, b] = faces[..] else {
            return None;
        };
        let [p, middle, q] = chain[..] else {
            return None;
        };
        let (facet, curved) = match (self.curved(a), self.curved(b)) {
            (Some(_), None) => (b, a),
            (None, Some(_)) => (a, b),
            _ => return None,
        };
        let mut members = self
            .groups
            .of
            .iter()
            .enumerate()
            .filter(|&(_, &g)| g == facet)
            .map(|(t, _)| t);
        let t = members.next()?;
        if members.next().is_some() || self.pinned.contains(&middle) {
            return None;
        }
        let fan = self.groups.fan_at(t, self.triangles, self.adjacency)?;
        (fan.seam.1 == middle
            && fan.seam_between(p, middle) == Some(curved)
            && fan.seam_between(middle, q) == Some(curved))
        .then_some(middle)
    }

    /// The curved face a chain is a fan's seam with, if it is one, and
    /// whether the fan is a wedge.
    fn fan_seam_of(&self, chain: &[u32], faces: &[usize]) -> Option<(usize, bool)> {
        let [a, b] = faces[..] else {
            return None;
        };
        let [p, q] = chain[..] else {
            return None;
        };
        [(a, b), (b, a)].into_iter().find_map(|(fan, curved)| {
            let f = self.fans.get(&fan)?;
            (f.seam_between(p, q) == Some(curved)).then_some((curved, f.across.is_some()))
        })
    }

    /// A fan's seam: the straight line between its two ends in the chart
    /// of the curved face `g`, lifted onto it and interpolated, with its
    /// image there. Its tolerance is how far the image's lift strays from
    /// the curve. `None` where `g` has no chart to draw the line in. A
    /// wedge's seams are interpolated at evenly spaced parameters, so the
    /// two share their knots either way along.
    fn lifted_seam(&self, chain: &[u32], g: usize, even: bool) -> Option<(Snapped, bool, Images)> {
        use ogeom_geom::{Curve2d as _, Curve3d as _};
        /// Points the lifted line is interpolated through.
        const SAMPLES: u32 = 16;
        let curved = self.curved(g)?;
        if curved.patch.is_some()
            || matches!(curved.shape, Canonical::Swept(_) | Canonical::Plane(_))
        {
            return None;
        }
        let [p, q] = chain[..] else {
            return None;
        };
        let (p, q) = (self.points[p as usize], self.points[q as usize]);
        let (pu, pv) = periodic(&curved.shape);
        let from = unwrapped(curved, p, self.tol)?;
        let (u, v) = chart(&curved.shape, q, self.tol)?;
        let near = |x: f64, c: f64, wraps: bool| {
            if wraps {
                c + ogeom_math::elementary::wrap_signed_angle(x - c)
            } else {
                x
            }
        };
        let to = (near(u, from.0, pu), near(v, from.1, pv));
        let at: Vec<Point> = (0..=SAMPLES)
            .map(|k| {
                let f = f64::from(k) / f64::from(SAMPLES);
                Point::new(
                    (to.0 - from.0).mul_add(f, from.0),
                    (to.1 - from.1).mul_add(f, from.1),
                    0.0,
                )
            })
            .collect();
        let mut lifted: Vec<Point> = at
            .iter()
            .map(|c| evaluate(&curved.shape, (c.x, c.y)))
            .collect();
        // The ends are the mesh's own vertices, which the edge's vertices
        // stand at.
        lifted[0] = p;
        lifted[SAMPLES as usize] = q;
        let parameters = if even {
            (0..=SAMPLES)
                .map(|k| f64::from(k) / f64::from(SAMPLES))
                .collect()
        } else {
            crate::fit::spaced(&lifted, crate::fit::Spacing::Centripetal, self.tol).ok()?
        };
        let curve = crate::fit::interpolate_at(&lifted, &parameters, 3, self.tol).ok()?;
        let image = crate::fit::interpolate_at(&at, &parameters, 3, self.tol).ok()?;
        let image: PlanarCurve = ogeom_geom::BSpline2d::new(
            image.knots().clone(),
            image
                .control_points()
                .iter()
                .map(|c| {
                    let c = c.point();
                    Point2::new(c.x, c.y)
                })
                .collect(),
            self.tol,
        )
        .ok()?
        .into();
        let curve: Curve = curve.into();
        let range = curve.domain();
        let mut deviation = self.tol.confusion() * 1e-2;
        for w in parameters.windows(2) {
            for f in [0.25, 0.5, 0.75] {
                let t = (w[1] - w[0]).mul_add(f, w[0]);
                let c = image.point_at(t, self.tol).ok()?;
                let on = evaluate(&curved.shape, (c.x, c.y));
                deviation = deviation.max(on.distance(curve.point_at(t, self.tol).ok()?));
            }
        }
        Some((
            Snapped::Open(curve, range, deviation.max(self.tol.confusion())),
            true,
            vec![(g, image, deviation)],
        ))
    }
}
