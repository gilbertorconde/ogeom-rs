//! Sections of a triangle mesh by a plane: polylines that know which side
//! the material is on.
//!
//! Which side of the plane each vertex lies on is decided once per vertex,
//! by [`Exact::orient3d`] against three points spanning the plane, and every
//! triangle reads its corners' sides from that one table. A triangle is cut
//! exactly when its corners disagree, the cut runs between the two edges
//! whose ends disagree, and two triangles sharing an edge agree on whether
//! it is cut because they share both ends' sides. The chains therefore
//! neither break nor double, whatever the rounding of the points on them.
//!
//! A vertex lying exactly in the plane counts as above it, on the side the
//! normal points to: the section is the one the plane gives moved an
//! infinitesimal distance against its normal. A cut edge whose upper end
//! lies in the plane is cut at that vertex itself, and the cuts through one
//! such vertex become one point of the chain. A triangle lying in the plane
//! is not cut, and is listed in [`Section::coplanar`]: a face in the plane
//! with the material below gives its outline, one with the material above
//! gives nothing, as the moved plane would.
//!
//! Chains are welded by vertex index, through the edges two triangles share,
//! never by position. Each runs with the material on its left seen from the
//! side the normal points to, so an outer boundary runs counter-clockwise
//! and a hole clockwise, and `outward` is the direction of travel crossed
//! with the plane's normal.

use ogeom_core::collections::FastMap;
use ogeom_core::parallel::map_ordered;
use ogeom_core::{Exact, OgeomResult, Predicates, Sign, ogeom_bail};
use ogeom_math::{Plane, Point, Vector};
use ogeom_topo::Triangulation;

/// One connected polyline of a section.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Chain {
    /// The points, in order along the chain. A closed chain does not repeat
    /// its first point at its end.
    pub points: Vec<Point>,
    /// At each point, the unit direction in the plane pointing out of the
    /// material, from the normals of the triangles cut there projected into
    /// the plane.
    pub outward: Vec<Vector>,
    /// Whether the chain closes on itself: every edge it crosses is shared
    /// by two triangles. An open chain ends on the mesh's boundary or on an
    /// edge listed in [`Section::unlinked`].
    pub closed: bool,
}

/// A mesh's section by one plane.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Section {
    /// The chains, in the order of the lowest-numbered triangle each cuts.
    pub chains: Vec<Chain>,
    /// The triangles lying in the plane, every corner exactly on it, in
    /// index order. None of them is cut.
    pub coplanar: Vec<u32>,
    /// The cut edges, as their two vertex indices smaller first, that are
    /// not crossed once each way by the triangles using them: an edge of
    /// more than two triangles, or of two wound the same way. A chain ends
    /// on each rather than joining across it, and the edges are sorted.
    pub unlinked: Vec<[u32; 2]>,
}

/// Section a mesh by a plane.
///
/// `outward` and the chains' direction come from the triangles' winding,
/// counter-clockwise about their outward normals as
/// [`Triangulation::triangles`] are. A triangle naming one vertex twice has
/// no area and is passed over.
///
/// # Errors
///
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle names
/// a vertex the mesh does not have, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a position or
/// the plane's origin is not finite.
pub fn section(mesh: &Triangulation, plane: &Plane) -> OgeomResult<Section> {
    checked(mesh, plane)?;
    Ok(cut(mesh, &Cutter::new(plane, 0.0)))
}

/// Section a mesh by a stack of planes, `plane` moved along its normal by
/// each of `offsets` in turn.
///
/// Each section is [`section`]'s by the moved plane, computed in parallel
/// and returned in the order of `offsets`, the same at any thread count.
///
/// # Errors
///
/// As [`section`], and [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric)
/// if an offset is not finite.
pub fn sections(mesh: &Triangulation, plane: &Plane, offsets: &[f64]) -> OgeomResult<Vec<Section>> {
    checked(mesh, plane)?;
    if let Some(d) = offsets.iter().find(|d| !d.is_finite()) {
        ogeom_bail!(Numeric, "a section offset of {d} is not a distance");
    }
    Ok(map_ordered(offsets, |_, &d| {
        cut(mesh, &Cutter::new(plane, d))
    }))
}

fn checked(mesh: &Triangulation, plane: &Plane) -> OgeomResult<()> {
    if !plane.origin().is_finite() {
        ogeom_bail!(Numeric, "the section plane's origin is not finite");
    }
    if let Some(k) = mesh.positions.iter().position(|p| !p.is_finite()) {
        ogeom_bail!(Numeric, "mesh vertex {k} is not finite");
    }
    let count = mesh.positions.len();
    if let Some(&k) = mesh
        .triangles
        .iter()
        .flatten()
        .find(|&&k| k as usize >= count)
    {
        ogeom_bail!(Range, "a triangle names vertex {k} of {count}");
    }
    Ok(())
}

/// Below, in, or above the plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Below,
    On,
    Above,
}

impl Side {
    const fn up(self) -> bool {
        !matches!(self, Self::Below)
    }
}

/// The plane as the side test sees it.
#[derive(Debug, Clone, Copy)]
struct Cutter {
    origin: Point,
    normal: Vector,
    /// Three points spanning the plane, counter-clockwise about `normal`.
    span: [[f64; 3]; 3],
}

impl Cutter {
    fn new(plane: &Plane, offset: f64) -> Self {
        let z = plane.normal();
        let normal = z.vector();
        let origin = plane.origin() + normal * offset;
        let u = z.any_perpendicular().vector();
        let v = normal.cross(u);
        Self {
            origin,
            normal,
            span: [
                origin.to_array(),
                (origin + u).to_array(),
                (origin + v).to_array(),
            ],
        }
    }

    fn side(&self, p: Point) -> Side {
        let [a, b, c] = self.span;
        // Positive below: seen from below, the counter-clockwise span turns
        // clockwise.
        match Exact::orient3d(a, b, c, p.to_array()) {
            Sign::Positive => Side::Below,
            Sign::Zero => Side::On,
            Sign::Negative => Side::Above,
        }
    }

    fn height(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.normal)
    }
}

/// A cut through one triangle, from the edge it leaves the lower side by
/// to the edge it returns by, in the triangle's winding.
#[derive(Debug, Clone, Copy)]
struct Segment {
    from: Edge,
    to: Edge,
    triangle: u32,
}

/// The cuts that start and end on one edge.
#[derive(Debug, Clone, Copy)]
struct Slot {
    starts: u32,
    ends: u32,
    start: u32,
    end: u32,
}

const NONE: u32 = u32::MAX;

/// An edge, as its two vertex indices smaller first.
type Edge = [u32; 2];

fn key(a: u32, b: u32) -> Edge {
    [a.min(b), a.max(b)]
}

/// One point of a chain: a vertex in the plane, or the crossing on an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Node {
    Vertex(u32),
    Edge(Edge),
}

fn cut(mesh: &Triangulation, cutter: &Cutter) -> Section {
    let sides: Vec<Side> = mesh.positions.iter().map(|&p| cutter.side(p)).collect();

    let mut segments = Vec::new();
    let mut coplanar = Vec::new();
    for (index, t) in mesh.triangles.iter().enumerate() {
        let triangle = u32::try_from(index).unwrap_or(NONE);
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            continue;
        }
        let s = t.map(|k| sides[k as usize]);
        let ups = s.iter().filter(|s| s.up()).count();
        if ups == 0 || ups == 3 {
            if s.iter().all(|&s| s == Side::On) {
                coplanar.push(triangle);
            }
            continue;
        }
        let (mut from, mut to) = ([0, 0], [0, 0]);
        for k in 0..3 {
            let (p, q) = (k, (k + 1) % 3);
            match (s[p].up(), s[q].up()) {
                (true, false) => from = key(t[p], t[q]),
                (false, true) => to = key(t[p], t[q]),
                _ => {}
            }
        }
        segments.push(Segment { from, to, triangle });
    }

    let mut slots: FastMap<Edge, Slot> = FastMap::default();
    let empty = Slot {
        starts: 0,
        ends: 0,
        start: NONE,
        end: NONE,
    };
    for (i, s) in segments.iter().enumerate() {
        let i = u32::try_from(i).unwrap_or(NONE);
        let slot = slots.entry(s.from).or_insert(empty);
        slot.starts += 1;
        slot.start = i;
        let slot = slots.entry(s.to).or_insert(empty);
        slot.ends += 1;
        slot.end = i;
    }
    let linked = |slot: &Slot| slot.starts <= 1 && slot.ends <= 1;
    let mut unlinked: Vec<[u32; 2]> = slots
        .iter()
        .filter(|(_, slot)| !linked(slot))
        .map(|(&k, _)| k)
        .collect();
    unlinked.sort_unstable();

    let next: Vec<u32> = segments
        .iter()
        .map(|s| {
            slots
                .get(&s.to)
                .filter(|x| linked(x))
                .map_or(NONE, |x| x.start)
        })
        .collect();
    let prev: Vec<u32> = segments
        .iter()
        .map(|s| {
            slots
                .get(&s.from)
                .filter(|x| linked(x))
                .map_or(NONE, |x| x.end)
        })
        .collect();

    let mut seen = vec![false; segments.len()];
    let mut chains = Vec::new();
    for first in 0..segments.len() {
        if seen[first] {
            continue;
        }
        let mut head = first;
        let mut closed = false;
        loop {
            let p = prev[head];
            if p == NONE {
                break;
            }
            if p as usize == first {
                closed = true;
                head = first;
                break;
            }
            head = p as usize;
        }
        let mut run = vec![head];
        seen[head] = true;
        let mut at = head;
        loop {
            let n = next[at];
            if n == NONE || n as usize == head {
                break;
            }
            at = n as usize;
            seen[at] = true;
            run.push(at);
        }
        if let Some(chain) = chain(mesh, cutter, &sides, &segments, &run, closed) {
            chains.push(chain);
        }
    }

    Section {
        chains,
        coplanar,
        unlinked,
    }
}

/// The chain through `run`, consecutive segments in order, or `None` where
/// it has fewer than two distinct points.
fn chain(
    mesh: &Triangulation,
    cutter: &Cutter,
    sides: &[Side],
    segments: &[Segment],
    run: &[usize],
    closed: bool,
) -> Option<Chain> {
    let node = |edge: Edge| {
        let [a, b] = edge;
        let up = if sides[a as usize].up() { a } else { b };
        if sides[up as usize] == Side::On {
            Node::Vertex(up)
        } else {
            Node::Edge(edge)
        }
    };
    // Each segment's lean out of the material, none for one that a vertex
    // in the plane shrinks to that vertex. A node takes the leans of the
    // segments either side of it: in a closed chain the node a segment
    // starts on joins it to the one before, in an open one the last
    // segment's end is a node of its own.
    let lean = |s: &Segment| {
        if node(s.from) == node(s.to) {
            Vector::ZERO
        } else {
            in_plane(mesh, cutter, s.triangle)
        }
    };
    let leans: Vec<Vector> = run.iter().map(|&s| lean(&segments[s])).collect();
    let mut nodes: Vec<(Node, Vector)> = Vec::with_capacity(run.len() + 1);
    for (i, &s) in run.iter().enumerate() {
        let before = if i > 0 {
            leans[i - 1]
        } else if closed && let Some(&last) = leans.last() {
            last
        } else {
            Vector::ZERO
        };
        push(&mut nodes, node(segments[s].from), leans[i] + before);
    }
    if !closed && let (Some(&last), Some(&w)) = (run.last(), leans.last()) {
        push(&mut nodes, node(segments[last].to), w);
    }
    if closed
        && nodes.len() > 1
        && let (Some(&(first, _)), Some(&(last, w))) = (nodes.first(), nodes.last())
        && first == last
    {
        nodes.pop();
        nodes[0].1 += w;
    }
    if nodes.len() < 2 {
        return None;
    }

    let points: Vec<Point> = nodes
        .iter()
        .map(|&(n, _)| point(mesh, cutter, sides, n))
        .collect();
    let count = points.len();
    let outward = (0..count)
        .map(|i| {
            let w = nodes[i].1;
            let length = w.magnitude();
            if length > 0.0 && length.is_finite() {
                return w / length;
            }
            // No cut triangle here leans out of the plane by a representable
            // amount: the chord through the neighbours, crossed with the
            // normal, is the direction the material's edge runs.
            let before = if i > 0 || closed {
                (i + count - 1) % count
            } else {
                i
            };
            let after = if i + 1 < count || closed {
                (i + 1) % count
            } else {
                i
            };
            let c = (points[after] - points[before]).cross(cutter.normal);
            let length = c.magnitude();
            if length > 0.0 && length.is_finite() {
                c / length
            } else {
                Vector::ZERO
            }
        })
        .collect();
    Some(Chain {
        points,
        outward,
        closed,
    })
}

/// Append a node, merged into the last one where both are the same vertex.
fn push(nodes: &mut Vec<(Node, Vector)>, node: Node, w: Vector) {
    if let Some(last) = nodes.last_mut()
        && last.0 == node
    {
        last.1 += w;
        return;
    }
    nodes.push((node, w));
}

/// A triangle's unit normal with its component along the plane's normal
/// removed: the sine of its tilt from the plane long, so the triangles
/// nearest lying in the plane, whose projected normals are least well
/// determined, weigh least.
fn in_plane(mesh: &Triangulation, cutter: &Cutter, triangle: u32) -> Vector {
    let [a, b, c] = mesh.triangles[triangle as usize].map(|k| mesh.positions[k as usize]);
    let n = (b - a).cross(c - a);
    let length = n.magnitude();
    if !(length > 0.0 && length.is_finite()) {
        return Vector::ZERO;
    }
    let n = n / length;
    n - cutter.normal * n.dot(cutter.normal)
}

/// Where a node lies: the vertex itself, or the crossing on the edge,
/// interpolated from the lower end so that every triangle sharing the edge
/// computes the same bits.
fn point(mesh: &Triangulation, cutter: &Cutter, sides: &[Side], node: Node) -> Point {
    match node {
        Node::Vertex(v) => mesh.positions[v as usize],
        Node::Edge(edge) => {
            let [a, b] = edge;
            let (low, high) = if sides[a as usize].up() {
                (b, a)
            } else {
                (a, b)
            };
            let (p, q) = (mesh.positions[low as usize], mesh.positions[high as usize]);
            let (hp, hq) = (cutter.height(p), cutter.height(q));
            // The exact sides put `p` strictly below and `q` strictly above;
            // rounded heights a hair the wrong way are held to the edge, and
            // two that round equal put the crossing anywhere on it, its
            // middle as well as any.
            let t = if hp == hq {
                0.5
            } else {
                (hp / (hp - hq)).clamp(0.0, 1.0)
            };
            p.lerp(q, t)
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]
mod tests {
    use super::*;
    use core::f64::consts::TAU;
    use ogeom_core::Tolerances;
    use ogeom_math::Direction;
    use proptest::prelude::*;

    const T: Tolerances = Tolerances::millimetres();

    fn mesh(positions: Vec<Point>, triangles: Vec<[u32; 3]>) -> Triangulation {
        Triangulation {
            positions,
            triangles,
            ..Triangulation::new()
        }
    }

    fn index(i: usize) -> u32 {
        u32::try_from(i).unwrap()
    }

    /// A wall of `sides` facets about the `z` axis at `radius`, `rings`
    /// rings of vertices from `z = 0` to `z = height`, wound outward when
    /// `outward` and inward otherwise, its vertices appended to `points`.
    fn wall(
        points: &mut Vec<Point>,
        triangles: &mut Vec<[u32; 3]>,
        radius: f64,
        height: f64,
        sides: usize,
        rings: usize,
        outward: bool,
    ) -> usize {
        let base = points.len();
        for j in 0..rings {
            let z = height * j as f64 / (rings - 1) as f64;
            for i in 0..sides {
                let a = TAU * i as f64 / sides as f64;
                points.push(Point::new(radius * a.cos(), radius * a.sin(), z));
            }
        }
        let at = |i: usize, j: usize| index(base + j * sides + i % sides);
        for j in 0..rings - 1 {
            for i in 0..sides {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                if outward {
                    triangles.extend([[a, b, c], [a, c, d]]);
                } else {
                    triangles.extend([[a, c, b], [a, d, c]]);
                }
            }
        }
        base
    }

    /// A closed cylinder about the `z` axis, capped by fans about the axis.
    fn cylinder(radius: f64, height: f64, sides: usize, rings: usize) -> Triangulation {
        let (mut points, mut triangles) = (Vec::new(), Vec::new());
        let base = wall(
            &mut points,
            &mut triangles,
            radius,
            height,
            sides,
            rings,
            true,
        );
        let top = base + (rings - 1) * sides;
        let (low, high) = (index(points.len()), index(points.len() + 1));
        points.extend([Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, height)]);
        for i in 0..sides {
            let j = (i + 1) % sides;
            triangles.push([low, index(base + j), index(base + i)]);
            triangles.push([high, index(top + i), index(top + j)]);
        }
        mesh(points, triangles)
    }

    /// A closed tube about the `z` axis between two radii, with annular caps.
    fn tube(inner: f64, outer: f64, height: f64, sides: usize) -> Triangulation {
        let (mut points, mut triangles) = (Vec::new(), Vec::new());
        let o = wall(&mut points, &mut triangles, outer, height, sides, 2, true);
        let i = wall(&mut points, &mut triangles, inner, height, sides, 2, false);
        for k in 0..sides {
            let m = (k + 1) % sides;
            let (ok, om, ik, im) = (index(o + k), index(o + m), index(i + k), index(i + m));
            triangles.extend([[ok, ik, im], [ok, im, om]]);
            let s = index(sides);
            triangles.extend([[ok + s, om + s, im + s], [ok + s, im + s, ik + s]]);
        }
        mesh(points, triangles)
    }

    /// The unit cube, each face a grid of `n` by `n` squares, welded and
    /// wound outward.
    fn cube(n: u32) -> Triangulation {
        let mut ids: FastMap<[u32; 3], u32> = FastMap::default();
        let (mut points, mut triangles) = (Vec::new(), Vec::new());
        let scale = f64::from(n);
        let mut vertex = |c: [u32; 3], points: &mut Vec<Point>| {
            *ids.entry(c).or_insert_with(|| {
                points.push(Point::new(
                    f64::from(c[0]) / scale,
                    f64::from(c[1]) / scale,
                    f64::from(c[2]) / scale,
                ));
                index(points.len() - 1)
            })
        };
        for axis in 0..3 {
            for side in [0, n] {
                let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                let at = |a: u32, b: u32| {
                    let mut c = [0; 3];
                    c[axis] = side;
                    c[u] = a;
                    c[v] = b;
                    c
                };
                for a in 0..n {
                    for b in 0..n {
                        let q = [at(a, b), at(a + 1, b), at(a + 1, b + 1), at(a, b + 1)]
                            .map(|c| vertex(c, &mut points));
                        // `u` then `v` turns counter-clockwise about `+axis`.
                        if side == n {
                            triangles.extend([[q[0], q[1], q[2]], [q[0], q[2], q[3]]]);
                        } else {
                            triangles.extend([[q[0], q[2], q[1]], [q[0], q[3], q[2]]]);
                        }
                    }
                }
            }
        }
        mesh(points, triangles)
    }

    fn plane(origin: Point, normal: Vector) -> Plane {
        Plane::through(origin, Direction::new(normal, T).unwrap())
    }

    fn horizontal(z: f64) -> Plane {
        plane(Point::new(0.0, 0.0, z), Vector::Z)
    }

    /// The area a chain encloses, positive counter-clockwise about `normal`.
    fn area(chain: &Chain, normal: Vector) -> f64 {
        let p0 = chain.points[0];
        let mut sum = 0.0;
        for w in chain.points.windows(2).skip(1) {
            sum += (w[0] - p0).cross(w[1] - p0).dot(normal);
        }
        sum * 0.5
    }

    fn radial(p: Point) -> Vector {
        Vector::new(p.x, p.y, 0.0) / p.x.hypot(p.y)
    }

    #[test]
    fn a_closed_cylinder_cut_across_is_one_circle_facing_away_from_its_axis() {
        let sides = 48;
        let m = cylinder(2.0, 3.0, sides, 4);
        // Between rings, where the facets' diagonals are cut as well as
        // their vertical edges, then through the second ring's vertices.
        let polygon = 0.5 * sides as f64 * 4.0 * (TAU / sides as f64).sin();
        let half = (TAU / sides as f64 / 2.0).cos();
        for (z, count) in [(0.7, 2 * sides), (1.0, sides)] {
            let s = section(&m, &horizontal(z)).unwrap();
            assert!(s.unlinked.is_empty() && s.coplanar.is_empty());
            assert_eq!(s.chains.len(), 1);
            let c = &s.chains[0];
            assert!(c.closed);
            assert_eq!(c.points.len(), count);
            assert_eq!(c.outward.len(), count);
            let mut on_edges = 0;
            for (p, w) in c.points.iter().zip(&c.outward) {
                assert!((p.z - z).abs() < 1e-15, "{p:?}");
                assert!(w.z.abs() < 1e-15);
                if (p.x.hypot(p.y) - 2.0).abs() < 1e-14 {
                    // On a vertical edge the two facets lean out equally.
                    assert!(w.dot(radial(*p)) > 1.0 - 1e-12, "{w:?} at {p:?}");
                    on_edges += 1;
                } else {
                    // Inside one facet, its own normal.
                    assert!(w.dot(radial(*p)) > half - 1e-12, "{w:?} at {p:?}");
                }
            }
            assert_eq!(on_edges, sides);
            // Counter-clockwise about the normal: the regular polygon's area.
            assert!((area(c, Vector::Z) - polygon).abs() < 1e-12);
        }
    }

    #[test]
    fn a_thin_wall_gives_two_chains_with_the_material_between_them() {
        let sides = 64;
        let (inner, outer) = (0.95, 1.0);
        let s = section(&tube(inner, outer, 1.0, sides), &horizontal(0.5)).unwrap();
        assert!(s.unlinked.is_empty());
        assert_eq!(s.chains.len(), 2);
        let polygon = |r: f64| 0.5 * sides as f64 * r * r * (TAU / sides as f64).sin();
        // No direction out of a facet is further than half a facet's turn
        // from the radius.
        let half = (TAU / sides as f64 / 2.0).cos();
        let mut faces = Vec::new();
        for c in &s.chains {
            assert!(c.closed);
            assert_eq!(c.points.len(), 2 * sides);
            let r = c.points.iter().map(|p| p.x.hypot(p.y)).fold(0.0, f64::max);
            let facing: Vec<f64> = c
                .points
                .iter()
                .zip(&c.outward)
                .map(|(p, w)| w.dot(radial(*p)))
                .collect();
            if (r - outer).abs() < 1e-12 {
                assert!(facing.iter().all(|&f| f > half - 1e-12));
                assert!((area(c, Vector::Z) - polygon(outer)).abs() < 1e-12);
                faces.push(1);
            } else {
                assert!((r - inner).abs() < 1e-12);
                // The bore's side faces the axis, and runs clockwise.
                assert!(facing.iter().all(|&f| f < -half + 1e-12));
                assert!((area(c, Vector::Z) + polygon(inner)).abs() < 1e-12);
                faces.push(-1);
            }
        }
        faces.sort_unstable();
        assert_eq!(faces, [-1, 1]);
    }

    #[test]
    fn an_open_sheet_gives_an_open_chain_from_edge_to_edge() {
        // A height field over [0, 4]², its normals up out of the material.
        let n = 40;
        let mut points = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                let (x, y) = (4.0 * i as f64 / n as f64, 4.0 * j as f64 / n as f64);
                points.push(Point::new(x, y, 0.1 * x.sin() * y.cos()));
            }
        }
        let at = |i: usize, j: usize| index(j * (n + 1) + i);
        let mut triangles = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                triangles.extend([[a, b, c], [a, c, d]]);
            }
        }
        let m = mesh(points, triangles);
        let s = section(&m, &plane(Point::new(2.05, 0.0, 0.0), Vector::X)).unwrap();
        assert!(s.unlinked.is_empty());
        assert_eq!(s.chains.len(), 1);
        let c = &s.chains[0];
        assert!(!c.closed);
        // One point on each row of vertical edges and one on each diagonal.
        assert_eq!(c.points.len(), 2 * n + 1);
        // Material on the left seen from `+x`, below the sheet: it runs to `-y`.
        assert_eq!(c.points[0].y, 4.0);
        assert_eq!(c.points[2 * n].y, 0.0);
        for (p, w) in c.points.iter().zip(&c.outward) {
            assert!((p.x - 2.05).abs() < 1e-15);
            assert!((p.z - 0.1 * p.x.sin() * p.y.cos()).abs() < 1e-3);
            assert!(w.x.abs() < 1e-15 && w.z > 0.99, "{w:?}");
        }
    }

    #[test]
    fn a_plane_through_vertices_and_edges_passes_them_once() {
        let m = cube(2);
        // Through the ring of mid vertices: every edge of the ring lies in
        // the plane, and each vertex is one point of the chain.
        let s = section(&m, &horizontal(0.5)).unwrap();
        assert!(s.unlinked.is_empty() && s.coplanar.is_empty());
        assert_eq!(s.chains.len(), 1);
        let c = &s.chains[0];
        assert!(c.closed);
        assert_eq!(c.points.len(), 8);
        assert!(c.points.iter().all(|p| m.positions.contains(p)));
        assert!((area(c, Vector::Z) - 1.0).abs() < 1e-15);
        for (p, w) in c.points.iter().zip(&c.outward) {
            let away = Vector::new(p.x - 0.5, p.y - 0.5, 0.0);
            assert!(w.dot(away) > 0.0, "{w:?} at {p:?}");
            assert!(w.z == 0.0);
        }

        // Along the top face: the plane moved down a hair cuts the sides
        // just under it, so the outline of the face, which lies in it.
        let s = section(&m, &horizontal(1.0)).unwrap();
        assert_eq!(s.coplanar.len(), 8);
        for &t in &s.coplanar {
            assert!(
                m.triangles[t as usize]
                    .iter()
                    .all(|&k| m.positions[k as usize].z == 1.0)
            );
        }
        assert_eq!(s.chains.len(), 1);
        assert_eq!(s.chains[0].points.len(), 8);
        assert!((area(&s.chains[0], Vector::Z) - 1.0).abs() < 1e-15);

        // Along the bottom face the moved plane misses the cube; facing the
        // other way it moves up into it and finds the outline again.
        let s = section(&m, &horizontal(0.0)).unwrap();
        assert_eq!(s.coplanar.len(), 8);
        assert!(s.chains.is_empty());
        let s = section(&m, &plane(Point::ORIGIN, -Vector::Z)).unwrap();
        assert_eq!(s.chains.len(), 1);
        let c = &s.chains[0];
        assert!((area(c, -Vector::Z) - 1.0).abs() < 1e-15);
        for (p, w) in c.points.iter().zip(&c.outward) {
            assert!(w.dot(Vector::new(p.x - 0.5, p.y - 0.5, 0.0)) > 0.0);
        }

        // Diagonally through two opposite vertical edges and the face
        // diagonals between them.
        let s = section(&m, &plane(Point::ORIGIN, Vector::new(1.0, -1.0, 0.0))).unwrap();
        assert!(s.unlinked.is_empty());
        assert_eq!(s.chains.len(), 1);
        assert!(s.chains[0].closed);
        assert!(
            (area(&s.chains[0], Vector::new(1.0, -1.0, 0.0) / 2f64.sqrt()) - 2f64.sqrt()).abs()
                < 1e-14
        );
    }

    #[test]
    fn every_lattice_plane_of_a_cube_cuts_it_once_or_not_at_all() {
        // Each axis both ways, through every layer of the 4-by-4 grids:
        // the faces the plane lies in are coplanar, and the plane moved a
        // hair against its normal is inside the cube or not.
        let m = cube(4);
        for axis in 0..3 {
            for sign in [1.0, -1.0] {
                for layer in 0..=4 {
                    let mut normal = [0.0; 3];
                    normal[axis] = sign;
                    let normal = Vector::from_array(normal);
                    let mut origin = [0.5; 3];
                    origin[axis] = f64::from(layer) / 4.0;
                    let s = section(&m, &plane(Point::from_array(origin), normal)).unwrap();
                    let inside = if sign > 0.0 { layer > 0 } else { layer < 4 };
                    let face = layer == 0 || layer == 4;
                    assert_eq!(s.coplanar.len(), if face { 32 } else { 0 });
                    assert!(s.unlinked.is_empty());
                    assert_eq!(s.chains.len(), usize::from(inside), "{axis} {sign} {layer}");
                    if let Some(c) = s.chains.first() {
                        assert!(c.closed);
                        assert_eq!(c.points.len(), 16);
                        assert!((area(c, normal) - 1.0).abs() < 1e-15);
                    }
                }
            }
        }
    }

    proptest! {
        /// Any plane through the middle of a cube cuts it in one convex
        /// loop, counter-clockwise about the normal, whose outward
        /// directions all point away from its inside.
        #[test]
        fn a_plane_through_a_cube_cuts_one_convex_loop_facing_out(
            (a, b) in (-3.2..3.2f64, 0.0..3.2f64),
            at in (0.3..0.7f64, 0.3..0.7f64, 0.3..0.7f64),
        ) {
            let m = cube(3);
            let normal = Vector::new(a.cos() * b.sin(), a.sin() * b.sin(), b.cos());
            let s = section(&m, &plane(Point::new(at.0, at.1, at.2), normal)).unwrap();
            prop_assert!(s.unlinked.is_empty() && s.coplanar.is_empty());
            prop_assert_eq!(s.chains.len(), 1);
            let c = &s.chains[0];
            prop_assert!(c.closed);
            prop_assert!(area(c, normal) > 0.0);
            let centre = Point::centroid(&c.points).unwrap();
            for (p, w) in c.points.iter().zip(&c.outward) {
                prop_assert!(w.dot(*p - centre) > 0.0);
                prop_assert!((w.magnitude() - 1.0).abs() < 1e-15);
                prop_assert!(w.dot(normal).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn a_vertex_touching_the_plane_from_below_adds_no_chain() {
        // A square pyramid with its apex in the plane.
        let points = vec![
            Point::new(-1.0, -1.0, 0.0),
            Point::new(1.0, -1.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(-1.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        let triangles = vec![
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
            [0, 2, 1],
            [0, 3, 2],
        ];
        let s = section(&mesh(points, triangles), &horizontal(1.0)).unwrap();
        assert_eq!(s, Section::default());
    }

    #[test]
    fn an_edge_of_three_sheets_is_reported_and_every_chain_ends_on_it() {
        // Three fins about the `z` axis, sharing the edge from vertex 0 to 1.
        let mut points = vec![Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0)];
        let mut triangles = Vec::new();
        for k in 0..3 {
            let a = TAU * f64::from(k) / 3.0;
            let (b, t) = (index(points.len()), index(points.len() + 1));
            points.extend([
                Point::new(a.cos(), a.sin(), -1.0),
                Point::new(a.cos(), a.sin(), 1.0),
            ]);
            triangles.extend([[0, b, t], [0, t, 1]]);
        }
        let s = section(&mesh(points, triangles), &horizontal(0.0)).unwrap();
        assert_eq!(s.unlinked, [[0, 1]]);
        assert_eq!(s.chains.len(), 3);
        for c in &s.chains {
            assert!(!c.closed);
            assert_eq!(c.points[0], Point::ORIGIN);
            assert!((c.points.last().unwrap().to_vector().magnitude() - 1.0).abs() < 1e-15);
        }

        // Two triangles wound against each other across their shared edge.
        let points = vec![
            Point::new(0.0, 0.0, -1.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
        ];
        let s = section(&mesh(points, vec![[0, 2, 1], [0, 3, 1]]), &horizontal(0.0)).unwrap();
        assert_eq!(s.unlinked, [[0, 1]]);
        assert_eq!(s.chains.len(), 2);
    }

    #[test]
    fn a_single_precision_mesh_is_cut_through_its_rings_exactly() {
        // Rounding to single precision is the point: the mesh a binary STL holds.
        #[allow(clippy::cast_possible_truncation)]
        let single = |v: f64| f64::from(v as f32);
        let mut m = cylinder(1.3, 1.0, 40, 11);
        for p in &mut m.positions {
            *p = Point::new(single(p.x), single(p.y), single(p.z));
        }
        let ring = single(0.3);
        let s = section(&m, &horizontal(ring)).unwrap();
        assert_eq!(s.chains.len(), 1);
        let c = &s.chains[0];
        assert!(c.closed && c.points.len() == 40);
        assert!(
            c.points
                .iter()
                .all(|p| p.z == ring && m.positions.contains(p))
        );
        // A hair off the ring, where the double-precision height falls.
        let s = section(&m, &horizontal(0.3)).unwrap();
        assert_eq!(s.chains.len(), 1);
        assert!(s.chains[0].closed && s.chains[0].points.len() == 80);
    }

    #[test]
    fn a_stack_is_the_sections_of_the_moved_planes_in_order() {
        let m = cylinder(1.0, 2.0, 32, 5);
        let base = plane(Point::new(0.1, 0.2, 0.0), Vector::new(0.1, 0.2, 1.0));
        let offsets = [1.5, 0.25, 0.5, 1.0, -3.0];
        let stack = sections(&m, &base, &offsets).unwrap();
        assert_eq!(stack.len(), offsets.len());
        for (s, &d) in stack.iter().zip(&offsets) {
            let moved = Plane::through(base.origin() + base.normal().vector() * d, base.normal());
            assert_eq!(*s, section(&m, &moved).unwrap());
        }
        assert_eq!(stack[1].chains.len(), 1);
        assert!(stack[4].chains.is_empty());
        assert!(sections(&m, &base, &[f64::NAN]).is_err());
    }

    #[test]
    fn a_bad_index_or_position_is_refused() {
        let mut m = cube(1);
        m.triangles.push([0, 1, 99]);
        assert!(section(&m, &horizontal(0.5)).is_err());
        let mut m = cube(1);
        m.positions[0].x = f64::NAN;
        assert!(section(&m, &horizontal(0.5)).is_err());
    }
}
