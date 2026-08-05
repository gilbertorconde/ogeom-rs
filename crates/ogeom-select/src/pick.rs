//! The pickable scene: a shape flattened for questions from the screen.
//!
//! Selection answers geometry's questions backwards. Modelling asks "where
//! is this face"; a pointer asks "which face is *here*" — and wants the
//! answer at sub-shape granularity, nearest first, thousands of times a
//! second. So the shape is flattened once: every face triangulated with a
//! record of which face owns which triangles, every edge discretized, every
//! vertex kept, and a bounding-volume hierarchy built over the triangles so
//! a ray touches `log` of the scene instead of all of it.
//!
//! The triangle-to-topology mapping is *stable*: triangle indices are
//! assigned in face traversal order at build time and never reshuffled, so
//! an application can hold a triangle index across frames and still name
//! the face it came from.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Aabb, Frame, Point, Point2, Vector};
use ogeom_mesh::{Deflection, polyline_of_edge, triangulate_face};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore, explore_unique};

/// What kind of sub-shape a hit resolved to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PickKind {
    /// A vertex within the aperture of the ray.
    Vertex,
    /// An edge within the aperture of the ray.
    Edge,
    /// The face the ray struck.
    Face,
}

/// One hit along a pick ray.
#[derive(Debug, Clone)]
pub struct Hit {
    /// The sub-shape the hit resolved to: the struck face, or an edge or
    /// vertex of it that passes within the aperture.
    pub shape: Shape,
    /// What kind of sub-shape [`Hit::shape`] is.
    pub kind: PickKind,
    /// Where the ray met the surface, in world coordinates.
    pub position: Point,
    /// The distance along the ray, for depth ordering.
    pub distance: f64,
    /// The index of the struck triangle — stable across the scene's life,
    /// resolvable back to its face through [`Pickable::triangle_face`].
    pub triangle: usize,
}

/// Whether a marquee keeps only what it swallows or everything it touches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Marquee {
    /// Select only shapes entirely inside the region.
    Inside,
    /// Select shapes the region touches at all.
    Crossing,
}

/// A ray into the scene.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Where the ray starts.
    pub origin: Point,
    /// Where it points; need not be unit length.
    pub direction: Vector,
}

#[derive(Debug, Clone)]
struct FaceRecord {
    shape: Shape,
    /// Triangle range `[start, end)` in the flat arrays.
    triangles: (usize, usize),
}

#[derive(Debug, Clone)]
struct EdgeRecord {
    shape: Shape,
    polyline: Vec<Point>,
}

#[derive(Debug, Clone)]
struct VertexRecord {
    shape: Shape,
    at: Point,
}

#[derive(Debug, Clone)]
struct BvhNode {
    bounds: Aabb,
    /// Leaf: range into `order`. Interior: children indices.
    content: BvhContent,
}

#[derive(Debug, Clone)]
enum BvhContent {
    Leaf(usize, usize),
    Split(usize, usize),
}

/// A shape flattened for picking, with its acceleration structure.
#[derive(Debug, Clone)]
pub struct Pickable {
    positions: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    faces: Vec<FaceRecord>,
    /// For each triangle, the index into `faces` — the stable mapping.
    owner: Vec<u32>,
    edges: Vec<EdgeRecord>,
    vertices: Vec<VertexRecord>,
    nodes: Vec<BvhNode>,
    /// Triangle indices as the BVH leaves order them.
    order: Vec<u32>,
}

/// How many triangles a BVH leaf holds before it splits.
const LEAF_SIZE: usize = 8;

impl Pickable {
    /// Flatten a shape for picking at the given deflection.
    ///
    /// The deflection is the level of detail: a scene for a distant view
    /// can be built coarse and cheap, one for close work fine — several
    /// scenes over the same shape are independent and each keeps its own
    /// stable triangle indices.
    ///
    /// # Errors
    ///
    /// As [`triangulate_face`] — a face that cannot be triangulated cannot
    /// be picked, and saying so beats silently ignoring it.
    pub fn build(
        model: &Model,
        shape: &Shape,
        deflection: Deflection,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let mut positions = Vec::new();
        let mut triangles: Vec<[u32; 3]> = Vec::new();
        let mut faces = Vec::new();
        let mut owner = Vec::new();

        for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
            let mesh = triangulate_face(model, &face, deflection, tol)?;
            let base = u32::try_from(positions.len()).map_err(|_| {
                ogeom_core::ogeom_err!(Construction, "a scene of four billion vertices")
            })?;
            let start = triangles.len();
            positions.extend_from_slice(&mesh.positions);
            for t in &mesh.triangles {
                triangles.push([t[0] + base, t[1] + base, t[2] + base]);
                owner.push(u32::try_from(faces.len()).unwrap_or(u32::MAX));
            }
            faces.push(FaceRecord {
                shape: face,
                triangles: (start, triangles.len()),
            });
        }

        let mut edges = Vec::new();
        for edge in explore_unique(model, shape, ShapeType::Edge)? {
            // A degenerate edge has no curve and no length; there is
            // nothing to point at.
            if let Ok(polyline) = polyline_of_edge(model, &edge, deflection, tol)
                && polyline.len() >= 2
            {
                edges.push(EdgeRecord {
                    shape: edge,
                    polyline,
                });
            }
        }

        let mut vertices = Vec::new();
        for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
            let Some(node) = model.node(&vertex) else {
                ogeom_bail!(Dangling, "vertex is not in this model");
            };
            if let Some(data) = node.data().as_vertex() {
                let at = vertex.transform(model.datums())?.apply(data.point);
                vertices.push(VertexRecord { shape: vertex, at });
            }
        }

        let (nodes, order) = build_bvh(&positions, &triangles);
        Ok(Self {
            positions,
            triangles,
            faces,
            owner,
            edges,
            vertices,
            nodes,
            order,
        })
    }

    /// The face that produced a triangle — the stable mapping back from
    /// tessellation to topology.
    #[must_use]
    pub fn triangle_face(&self, triangle: usize) -> Option<&Shape> {
        let face = *self.owner.get(triangle)?;
        self.faces.get(face as usize).map(|record| &record.shape)
    }

    /// How many triangles the scene holds.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Every hit along a ray, nearest first.
    ///
    /// `aperture` is the pick radius in world units — how close the ray
    /// must pass to an edge or vertex of the struck face to resolve to it
    /// instead of the face. Zero apertures pick faces only.
    #[must_use]
    pub fn pick(&self, ray: Ray, aperture: f64) -> Vec<Hit> {
        let length = ray.direction.magnitude();
        if !(length.is_finite()) || length <= 0.0 {
            return Vec::new();
        }
        let direction = ray.direction / length;

        let mut hits = Vec::new();
        if !self.nodes.is_empty() {
            let mut stack = vec![0_usize];
            while let Some(at) = stack.pop() {
                let node = &self.nodes[at];
                if !ray_touches_box(&ray.origin, &direction, &node.bounds) {
                    continue;
                }
                match node.content {
                    BvhContent::Split(left, right) => {
                        stack.push(left);
                        stack.push(right);
                    }
                    BvhContent::Leaf(start, end) => {
                        for &index in &self.order[start..end] {
                            let triangle = index as usize;
                            let [a, b, c] =
                                self.triangles[triangle].map(|k| self.positions[k as usize]);
                            if let Some(t) = ray_triangle(&ray.origin, &direction, a, b, c) {
                                let position = ray.origin + direction * t;
                                hits.push((t, triangle, position));
                            }
                        }
                    }
                }
            }
        }
        hits.sort_by(|x, y| x.0.total_cmp(&y.0));
        // A ray along a triangulation diagonal strikes both triangles of
        // one face at one depth; the pick reports the face once.
        hits.dedup_by(|next, held| {
            self.owner[next.1] == self.owner[held.1]
                && (next.0 - held.0).abs() <= 1e-9 * (1.0 + held.0.abs())
        });

        hits.into_iter()
            .map(|(t, triangle, position)| {
                let (shape, kind) = self.resolve(triangle, &ray.origin, &direction, t, aperture);
                Hit {
                    shape,
                    kind,
                    position,
                    distance: t,
                    triangle,
                }
            })
            .collect()
    }

    /// The nearest hit along a ray, if any.
    #[must_use]
    pub fn pick_first(&self, ray: Ray, aperture: f64) -> Option<Hit> {
        self.pick(ray, aperture).into_iter().next()
    }

    /// Sub-shape resolution: the vertex, else the edge, else the face.
    fn resolve(
        &self,
        triangle: usize,
        origin: &Point,
        direction: &Vector,
        along: f64,
        aperture: f64,
    ) -> (Shape, PickKind) {
        let face = &self.faces[self.owner[triangle] as usize];
        if aperture > 0.0 {
            let struck = *origin + *direction * along;
            let mut best_vertex: Option<(f64, usize)> = None;
            for (i, vertex) in self.vertices.iter().enumerate() {
                // Near the strike point, not merely near the infinite ray.
                if vertex.at.distance(struck) <= aperture {
                    let off = ray_point_distance(origin, direction, vertex.at);
                    if off <= aperture && best_vertex.is_none_or(|(held, _)| off < held) {
                        best_vertex = Some((off, i));
                    }
                }
            }
            if let Some((_, i)) = best_vertex {
                return (self.vertices[i].shape.clone(), PickKind::Vertex);
            }
            let mut best_edge: Option<(f64, usize)> = None;
            for (i, edge) in self.edges.iter().enumerate() {
                for pair in edge.polyline.windows(2) {
                    let (off, near) = ray_segment_distance(origin, direction, pair[0], pair[1]);
                    if off <= aperture
                        && near.distance(struck) <= aperture * 4.0
                        && best_edge.is_none_or(|(held, _)| off < held)
                    {
                        best_edge = Some((off, i));
                    }
                }
            }
            if let Some((_, i)) = best_edge {
                return (self.edges[i].shape.clone(), PickKind::Edge);
            }
        }
        (face.shape.clone(), PickKind::Face)
    }

    /// Marquee selection: the faces, edges and vertices a screen-space
    /// rectangle selects under an orthographic view.
    ///
    /// The view is a frame whose `z` looks along the view direction and
    /// whose `x`/`y` are the screen axes; everything is projected onto its
    /// `xy` plane. [`Marquee::Inside`] keeps sub-shapes wholly within the
    /// rectangle; [`Marquee::Crossing`] keeps anything touching it.
    #[must_use]
    pub fn select_rectangle(
        &self,
        view: &Frame,
        low: Point2,
        high: Point2,
        mode: Marquee,
    ) -> Vec<Shape> {
        let inside =
            |p: Point2| -> bool { p.x >= low.x && p.x <= high.x && p.y >= low.y && p.y <= high.y };
        self.select_where(view, mode, &inside)
    }

    /// Marquee selection by an arbitrary closed polygon in screen space.
    ///
    /// Even–odd containment, so the polygon may be concave; the last point
    /// closes to the first.
    #[must_use]
    pub fn select_polygon(&self, view: &Frame, polygon: &[Point2], mode: Marquee) -> Vec<Shape> {
        if polygon.len() < 3 {
            return Vec::new();
        }
        let inside = |p: Point2| -> bool {
            let mut odd = false;
            let n = polygon.len();
            for i in 0..n {
                let (a, b) = (polygon[i], polygon[(i + 1) % n]);
                if (a.y > p.y) != (b.y > p.y) {
                    let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                    if p.x < x {
                        odd = !odd;
                    }
                }
            }
            odd
        };
        self.select_where(view, mode, &inside)
    }

    fn select_where(
        &self,
        view: &Frame,
        mode: Marquee,
        inside: &dyn Fn(Point2) -> bool,
    ) -> Vec<Shape> {
        let project = |p: Point| -> Point2 {
            let local = view.to_local(p);
            Point2::new(local.x, local.y)
        };
        let mut selected: Vec<Shape> = Vec::new();
        let mut keep = |shape: &Shape| {
            if !selected.iter().any(|held| held.is_same(shape)) {
                selected.push(shape.clone());
            }
        };

        for face in &self.faces {
            let (start, end) = face.triangles;
            let mut all = true;
            let mut any = false;
            for triangle in &self.triangles[start..end] {
                for &k in triangle {
                    let hit = inside(project(self.positions[k as usize]));
                    all &= hit;
                    any |= hit;
                }
            }
            if start == end {
                continue;
            }
            let wanted = match mode {
                Marquee::Inside => all,
                Marquee::Crossing => any,
            };
            if wanted {
                keep(&face.shape);
            }
        }
        for edge in &self.edges {
            let mut all = true;
            let mut any = false;
            for p in &edge.polyline {
                let hit = inside(project(*p));
                all &= hit;
                any |= hit;
            }
            let wanted = match mode {
                Marquee::Inside => all,
                Marquee::Crossing => any,
            };
            if wanted {
                keep(&edge.shape);
            }
        }
        for vertex in &self.vertices {
            if inside(project(vertex.at)) {
                keep(&vertex.shape);
            }
        }
        selected
    }
}

/// Build the BVH: median split on the longest axis of the centroid bounds.
fn build_bvh(positions: &[Point], triangles: &[[u32; 3]]) -> (Vec<BvhNode>, Vec<u32>) {
    let mut order: Vec<u32> = (0..triangles.len())
        .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
        .collect();
    let mut nodes = Vec::new();
    if triangles.is_empty() {
        return (nodes, order);
    }
    let centroid = |i: u32| -> Point {
        let [a, b, c] = triangles[i as usize].map(|k| positions[k as usize]);
        Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        )
    };
    let bounds_of = |slice: &[u32]| -> Aabb {
        slice.iter().fold(Aabb::EMPTY, |acc, &i| {
            triangles[i as usize]
                .iter()
                .fold(acc, |acc, &k| acc.with_point(positions[k as usize]))
        })
    };

    // An explicit stack of (range, parent slot) builds the tree without
    // recursion; each pending range claims its node index up front so the
    // parent can point at it.
    let mut pending = vec![(0_usize, order.len())];
    let mut slots: Vec<Option<(usize, bool)>> = vec![None];
    while let Some(((start, end), slot)) = pending.pop().zip(slots.pop()) {
        let index = nodes.len();
        if let Some((parent, right)) = slot
            && let BvhContent::Split(l, r) = &mut nodes[parent].content
        {
            if right {
                *r = index;
            } else {
                *l = index;
            }
        }
        let bounds = bounds_of(&order[start..end]);
        if end - start <= LEAF_SIZE {
            nodes.push(BvhNode {
                bounds,
                content: BvhContent::Leaf(start, end),
            });
            continue;
        }
        // Median split along the widest centroid spread.
        let size = bounds.size();
        let axis = if size.x >= size.y && size.x >= size.z {
            0
        } else if size.y >= size.z {
            1
        } else {
            2
        };
        let key = |i: u32| -> f64 {
            let c = centroid(i);
            match axis {
                0 => c.x,
                1 => c.y,
                _ => c.z,
            }
        };
        order[start..end].sort_by(|&a, &b| key(a).total_cmp(&key(b)));
        let mid = start + (end - start) / 2;
        nodes.push(BvhNode {
            bounds,
            content: BvhContent::Split(0, 0),
        });
        pending.push((mid, end));
        slots.push(Some((index, true)));
        pending.push((start, mid));
        slots.push(Some((index, false)));
    }
    (nodes, order)
}

/// Whether a ray touches an axis-aligned box: the slab test.
fn ray_touches_box(origin: &Point, direction: &Vector, bounds: &Aabb) -> bool {
    let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
        return false;
    };
    let mut t_near = 0.0_f64;
    let mut t_far = f64::INFINITY;
    for axis in 0..3 {
        let (o, d, lo, hi) = match axis {
            0 => (origin.x, direction.x, low.x, high.x),
            1 => (origin.y, direction.y, low.y, high.y),
            _ => (origin.z, direction.z, low.z, high.z),
        };
        if d.abs() < f64::MIN_POSITIVE {
            if o < lo || o > hi {
                return false;
            }
            continue;
        }
        let (mut t0, mut t1) = ((lo - o) / d, (hi - o) / d);
        if t0 > t1 {
            core::mem::swap(&mut t0, &mut t1);
        }
        t_near = t_near.max(t0);
        t_far = t_far.min(t1);
        if t_near > t_far {
            return false;
        }
    }
    true
}

/// Möller–Trumbore, front and back faces alike; `None` off the triangle or
/// behind the origin.
fn ray_triangle(origin: &Point, direction: &Vector, a: Point, b: Point, c: Point) -> Option<f64> {
    let e1 = b - a;
    let e2 = c - a;
    let h = direction.cross(e2);
    let det = e1.dot(h);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv = 1.0 / det;
    let s = *origin - a;
    let u = s.dot(h) * inv;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = direction.dot(q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(q) * inv;
    (t > 0.0).then_some(t)
}

/// Distance from a ray's carrier line to a point.
fn ray_point_distance(origin: &Point, direction: &Vector, p: Point) -> f64 {
    let to = p - *origin;
    to.cross(*direction).magnitude()
}

/// Distance from a ray's carrier line to a segment, and the nearest point
/// on the segment.
fn ray_segment_distance(origin: &Point, direction: &Vector, a: Point, b: Point) -> (f64, Point) {
    // Closest points between two lines, the segment's parameter clamped.
    let u = *direction;
    let v = b - a;
    let w = *origin - a;
    let uu = u.dot(u);
    let uv = u.dot(v);
    let vv = v.dot(v);
    let uw = u.dot(w);
    let vw = v.dot(w);
    let denominator = uu * vv - uv * uv;
    let s = if denominator.abs() < f64::MIN_POSITIVE {
        0.0
    } else {
        ((uu * vw - uv * uw) / denominator).clamp(0.0, 1.0)
    };
    let on_segment = a + v * s;
    (
        ray_point_distance(origin, direction, on_segment),
        on_segment,
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_math::Direction;

    const T: Tolerances = Tolerances::millimetres();

    fn box_scene() -> (Model, Pickable) {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let scene = Pickable::build(&model, &solid, Deflection::default(), T).unwrap();
        (model, scene)
    }

    fn down_at(x: f64, y: f64) -> Ray {
        Ray {
            origin: Point::new(x, y, 10.0),
            direction: Vector::new(0.0, 0.0, -1.0),
        }
    }

    #[test]
    fn a_ray_into_a_face_hits_it_nearest_first_and_names_the_face() {
        let (model, scene) = box_scene();
        let hits = scene.pick(down_at(1.0, 1.0), 0.0);
        // In through the top, out through the bottom: two faces, in depth
        // order.
        assert_eq!(hits.len(), 2);
        assert!(hits[0].distance < hits[1].distance);
        assert_eq!(hits[0].kind, PickKind::Face);
        assert_relative_eq!(hits[0].position.z, 2.0, epsilon = 1e-12);
        assert_relative_eq!(hits[1].position.z, 0.0, epsilon = 1e-12);
        // The stable mapping agrees with the resolved hit.
        let owner = scene.triangle_face(hits[0].triangle).unwrap();
        assert!(owner.is_same(&hits[0].shape));
        assert_eq!(model.kind_of(owner).unwrap(), ShapeType::Face);
    }

    #[test]
    fn an_aperture_resolves_edges_and_vertices_before_faces() {
        let (model, scene) = box_scene();
        // Straight down the box's corner: with an aperture, the corner
        // vertex wins.
        let corner = scene.pick_first(down_at(0.001, 0.001), 0.05).unwrap();
        assert_eq!(corner.kind, PickKind::Vertex);
        assert_eq!(model.kind_of(&corner.shape).unwrap(), ShapeType::Vertex);
        // Down the middle of the x = 0..2, y = 0 top edge: the edge wins.
        let edge = scene.pick_first(down_at(1.0, 0.001), 0.05).unwrap();
        assert_eq!(edge.kind, PickKind::Edge);
        assert_eq!(model.kind_of(&edge.shape).unwrap(), ShapeType::Edge);
        // Mid-face, the same aperture is nowhere near anything lower.
        let face = scene.pick_first(down_at(1.0, 1.0), 0.05).unwrap();
        assert_eq!(face.kind, PickKind::Face);
    }

    #[test]
    fn a_missing_ray_hits_nothing() {
        let (_, scene) = box_scene();
        assert!(scene.pick(down_at(5.0, 5.0), 0.1).is_empty());
        assert!(scene.pick_first(down_at(-1.0, -1.0), 0.0).is_none());
    }

    #[test]
    fn the_triangle_mapping_is_stable_and_total() {
        let (_, scene) = box_scene();
        assert_eq!(scene.triangle_count(), 12, "a box is twelve triangles");
        for t in 0..scene.triangle_count() {
            assert!(scene.triangle_face(t).is_some());
        }
        // Rebuilding the same shape gives the same mapping: stability is
        // determinism plus traversal order.
        let (_, again) = box_scene();
        assert_eq!(again.triangle_count(), scene.triangle_count());
    }

    #[test]
    fn rectangle_selection_distinguishes_inside_from_crossing() {
        let (_, scene) = box_scene();
        let view = Frame::new(Point::new(0.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        // The whole footprint: everything is inside.
        let all = scene.select_rectangle(
            &view,
            Point2::new(-0.1, -0.1),
            Point2::new(2.1, 2.1),
            Marquee::Inside,
        );
        // 6 faces + 12 edges + 8 vertices.
        assert_eq!(all.len(), 26);
        // Half the footprint: nothing whole fits inside, plenty crosses.
        let inside = scene.select_rectangle(
            &view,
            Point2::new(-0.1, -0.1),
            Point2::new(1.0, 2.1),
            Marquee::Inside,
        );
        let crossing = scene.select_rectangle(
            &view,
            Point2::new(-0.1, -0.1),
            Point2::new(1.0, 2.1),
            Marquee::Crossing,
        );
        assert!(inside.len() < crossing.len());
        // The x = 0 side's faces and edges cross; the x = 2 side's do not.
        assert!(!crossing.is_empty());
    }

    #[test]
    fn polygon_selection_honours_concavity() {
        let (_, scene) = box_scene();
        let view = Frame::new(Point::new(0.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        // An L-shaped marquee that covers the (0, 0) corner but leaves the
        // (2, 2) corner out.
        let marquee = [
            Point2::new(-0.5, -0.5),
            Point2::new(1.0, -0.5),
            Point2::new(1.0, 1.0),
            Point2::new(-0.5, 1.0),
        ];
        let picked = scene.select_polygon(&view, &marquee, Marquee::Crossing);
        assert!(!picked.is_empty());
        let none = scene.select_polygon(&view, &marquee[..2], Marquee::Crossing);
        assert!(none.is_empty(), "two points are not a polygon");
    }
}
