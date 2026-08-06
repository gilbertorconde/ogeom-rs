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
    /// Whether [`Hit::position`] and [`Hit::distance`] were refined onto
    /// the face's exact surface, or stand on the tessellation alone.
    pub refined: bool,
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
    /// The face's surface, for exact refinement.
    surface: ogeom_geom::SurfaceGeometry,
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

/// One face's draft reading: the signed angle range against the pull.
#[derive(Debug, Clone)]
pub struct FaceDraft {
    /// The face.
    pub face: Shape,
    /// The smallest signed draft angle sampled, radians.
    pub min: f64,
    /// The largest.
    pub max: f64,
}

/// One face's least material depth.
#[derive(Debug, Clone)]
pub struct FaceThickness {
    /// The face.
    pub face: Shape,
    /// The least inward distance to the opposite wall; infinity where the
    /// rays found none.
    pub least: f64,
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
            let surface = {
                let Some(node) = model.node(&face) else {
                    ogeom_bail!(Dangling, "face is not in this model");
                };
                let ogeom_topo::NodeData::Face(data) = node.data() else {
                    ogeom_bail!(Construction, "face node holds no face data");
                };
                let Some(surface) = model.geometry().surface(data.surface) else {
                    ogeom_bail!(Dangling, "face refers to a surface not in this model");
                };
                surface.clone()
            };
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
                surface,
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

    /// Draft analysis: each face's signed angle range against a pull
    /// direction, sampled at this scene's triangles — the deflection the
    /// scene was built at is the sampling, stated by construction.
    ///
    /// The angle at a sample is `asin(n · pull)`: positive where the face
    /// tilts its outward normal along the pull — drafted — negative where
    /// it undercuts, zero on a straight wall.
    #[must_use]
    pub fn draft_analysis(&self, pull: ogeom_math::Direction) -> Vec<FaceDraft> {
        let pull = pull.vector();
        self.faces
            .iter()
            .map(|record| {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for t in record.triangles.0..record.triangles.1 {
                    let [a, b, c] = self.triangles[t];
                    let (pa, pb, pc) = (
                        self.positions[a as usize],
                        self.positions[b as usize],
                        self.positions[c as usize],
                    );
                    let n = (pb - pa).cross(pc - pa);
                    let m = n.magnitude();
                    if m <= f64::MIN_POSITIVE {
                        continue;
                    }
                    let angle = (n.dot(pull) / m).clamp(-1.0, 1.0).asin();
                    lo = lo.min(angle);
                    hi = hi.max(angle);
                }
                FaceDraft {
                    face: record.shape.clone(),
                    min: lo,
                    max: hi,
                }
            })
            .filter(|d| d.min.is_finite())
            .collect()
    }

    /// Thickness analysis: for each face, the least material depth found by
    /// casting inward from its triangle centroids to the first opposite
    /// wall — sampled at this scene's own triangles, as fine as its
    /// deflection.
    ///
    /// Faces whose inward rays strike nothing — an open sheet — report
    /// infinity, which is the honest answer for material with no far side.
    #[must_use]
    pub fn thickness_analysis(&self) -> Vec<FaceThickness> {
        self.faces
            .iter()
            .map(|record| {
                let mut least = f64::INFINITY;
                for t in record.triangles.0..record.triangles.1 {
                    let [a, b, c] = self.triangles[t];
                    let (pa, pb, pc) = (
                        self.positions[a as usize],
                        self.positions[b as usize],
                        self.positions[c as usize],
                    );
                    let n = (pb - pa).cross(pc - pa);
                    let m = n.magnitude();
                    if m <= f64::MIN_POSITIVE {
                        continue;
                    }
                    let inward = n / -m;
                    let centroid = Point::new(
                        (pa.x + pb.x + pc.x) / 3.0,
                        (pa.y + pb.y + pc.y) / 3.0,
                        (pa.z + pb.z + pc.z) / 3.0,
                    );
                    // Step off the surface so the ray does not strike home.
                    let skin = 1e-9 * (1.0 + centroid.to_vector().magnitude());
                    let hits = self.pick(
                        Ray {
                            origin: centroid + inward * skin,
                            direction: inward,
                        },
                        0.0,
                    );
                    for hit in hits {
                        if hit.triangle != t && hit.distance > skin {
                            least = least.min(hit.distance + skin);
                            break;
                        }
                    }
                }
                FaceThickness {
                    face: record.shape.clone(),
                    least,
                }
            })
            .collect()
    }

    /// How many faces the scene holds.
    ///
    /// Faces are numbered in traversal order at build time, so the same index
    /// names the same face in every scene built over one shape — which is what
    /// lets a [`PickHierarchy`] carry an answer from one level to another.
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// The index of the face a triangle came from.
    #[must_use]
    pub fn face_index(&self, triangle: usize) -> Option<usize> {
        self.owner.get(triangle).map(|&f| f as usize)
    }

    /// Which faces this scene's tessellation brings within `slack` of the
    /// ray — a *conservative* answer, by which nothing that could be hit is
    /// left out.
    ///
    /// Conservative is the whole point. A coarse tessellation stands within
    /// its own chord of the true surface, so a ray that strikes the surface
    /// passes within that chord of the coarse mesh; testing the mesh's boxes
    /// widened by the chord therefore cannot miss a face a finer scene would
    /// find. The reverse is not true, and does not need to be: a face this
    /// returns and a finer level rejects simply costs one more test.
    #[must_use]
    pub fn faces_near(&self, ray: Ray, slack: f64) -> Vec<bool> {
        let mut near = vec![false; self.faces.len()];
        let length = ray.direction.magnitude();
        if !length.is_finite() || length <= 0.0 || self.nodes.is_empty() {
            return near;
        }
        let direction = ray.direction / length;
        let mut stack = vec![0_usize];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at];
            if !ray_touches_box(&ray.origin, &direction, &node.bounds.expanded(slack)) {
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
                        let owner = self.owner[triangle] as usize;
                        if near[owner] {
                            continue;
                        }
                        let bounds = self.triangles[triangle]
                            .iter()
                            .fold(Aabb::EMPTY, |acc, &k| {
                                acc.with_point(self.positions[k as usize])
                            })
                            .expanded(slack);
                        if ray_touches_box(&ray.origin, &direction, &bounds) {
                            near[owner] = true;
                        }
                    }
                }
            }
        }
        near
    }

    /// Every hit along a ray, nearest first.
    ///
    /// `aperture` is the pick radius in world units — how close the ray
    /// must pass to an edge or vertex of the struck face to resolve to it
    /// instead of the face. Zero apertures pick faces only.
    #[must_use]
    pub fn pick(&self, ray: Ray, aperture: f64) -> Vec<Hit> {
        self.pick_within(ray, aperture, None)
    }

    /// As [`Pickable::pick`], considering only the faces `allowed` marks.
    ///
    /// The mask is indexed by face, in the order [`Pickable::face_count`]
    /// counts them.
    #[must_use]
    pub fn pick_within(&self, ray: Ray, aperture: f64, allowed: Option<&[bool]>) -> Vec<Hit> {
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
                            if let Some(mask) = allowed
                                && !mask
                                    .get(self.owner[triangle] as usize)
                                    .copied()
                                    .unwrap_or(false)
                            {
                                continue;
                            }
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
                    refined: false,
                }
            })
            .collect()
    }

    /// The nearest hit along a ray, if any.
    #[must_use]
    pub fn pick_first(&self, ray: Ray, aperture: f64) -> Option<Hit> {
        self.pick(ray, aperture).into_iter().next()
    }

    /// As [`Pickable::pick`], with each hit refined onto its face's exact
    /// surface.
    ///
    /// The tessellation finds the hits and orders them; the analytic
    /// surface then answers *where*, exactly: the ray is intersected with
    /// the struck face's own geometry and the crossing nearest the mesh
    /// answer replaces it. A hit whose exact refinement resolves nothing —
    /// a grazing ray, a surface kind the intersector seeds poorly — keeps
    /// the tessellated answer and says so through [`Hit::refined`].
    #[must_use]
    pub fn pick_refined(&self, ray: Ray, aperture: f64, tol: Tolerances) -> Vec<Hit> {
        let hits = self.pick(ray, aperture);
        self.refine(hits, ray, tol)
    }

    /// Put already-found hits onto their faces' exact surfaces.
    ///
    /// Separated from [`Pickable::pick_refined`] because refining is per hit
    /// and knows nothing about how the hits were found — which is what lets a
    /// [`PickHierarchy`] narrow the search and refine the same way.
    #[must_use]
    pub fn refine(&self, hits: Vec<Hit>, ray: Ray, tol: Tolerances) -> Vec<Hit> {
        let mut hits = hits;
        let length = ray.direction.magnitude();
        if length <= 0.0 || !length.is_finite() {
            return hits;
        }
        let direction = ray.direction / length;
        for hit in &mut hits {
            let face = &self.faces[self.owner[hit.triangle] as usize];
            // A segment comfortably bracketing the mesh answer.
            let margin = (hit.distance * 0.5).max(1.0);
            let Ok(segment) = ogeom_geom::LineCurve::segment(
                ray.origin + direction * (hit.distance - margin).max(0.0),
                ray.origin + direction * (hit.distance + margin),
                tol,
            ) else {
                continue;
            };
            let curve: ogeom_geom::Curve = segment.into();
            let Ok(pierced) = ogeom_intersect::intersect_curve_surface(
                &curve,
                &face.surface,
                ogeom_intersect::CurveSurfaceOptions::default(),
                tol,
            ) else {
                continue;
            };
            let Some(best) = pierced
                .crossings
                .iter()
                .min_by(|x, y| {
                    x.point
                        .distance(hit.position)
                        .total_cmp(&y.point.distance(hit.position))
                })
                .filter(|c| c.point.distance(hit.position) < margin)
            else {
                continue;
            };
            hit.position = best.point;
            hit.distance = (best.point - ray.origin).dot(direction);
            hit.refined = true;
        }
        hits
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

/// One structure over several deflections, and a pick that descends it.
///
/// A scene built fine answers precisely and costs what fine costs; a scene
/// built coarse costs little and answers coarsely. Held separately, an
/// application picking against both pays for both and has to reconcile two
/// numberings. Held here they are one thing: built in one call over one
/// shape, so face index `i` names the same face at every level, and queried
/// coarse to fine, so the fine level is asked only about the faces the coarse
/// one could not rule out.
///
/// **The answer does not change.** Descending is a way of *not asking* the
/// fine level about most of the scene; what it is asked, it answers exactly as
/// it would alone. That holds because each level's mesh stands within its own
/// stated chord of the true surface, so widening the coarse level's boxes by
/// the two chords cannot rule out a face the fine level would hit. The
/// equality is pinned by test, not asserted here.
#[derive(Debug, Clone)]
pub struct PickHierarchy {
    /// Coarsest first, finest last.
    levels: Vec<Pickable>,
    /// Each level's chord, in the same order.
    chords: Vec<f64>,
}

impl PickHierarchy {
    /// Build one scene per deflection, over one shape.
    ///
    /// The deflections are sorted coarsest first; duplicates are kept, since
    /// a caller that asked for two identical levels asked for two.
    ///
    /// # Errors
    ///
    /// As [`Pickable::build`], plus
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if no
    /// deflection is given — a hierarchy of no levels answers nothing.
    pub fn build(
        model: &Model,
        shape: &Shape,
        deflections: &[Deflection],
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        if deflections.is_empty() {
            ogeom_bail!(Construction, "a pick hierarchy needs at least one level");
        }
        let mut ordered: Vec<Deflection> = deflections.to_vec();
        ordered.sort_by(|a, b| b.chord.total_cmp(&a.chord));
        let mut levels = Vec::with_capacity(ordered.len());
        let mut chords = Vec::with_capacity(ordered.len());
        for deflection in ordered {
            levels.push(Pickable::build(model, shape, deflection, tol)?);
            chords.push(deflection.chord);
        }
        Ok(Self { levels, chords })
    }

    /// How many levels the hierarchy holds.
    #[must_use]
    pub fn level_count(&self) -> usize {
        self.levels.len()
    }

    /// One level, coarsest first.
    #[must_use]
    pub fn level(&self, index: usize) -> Option<&Pickable> {
        self.levels.get(index)
    }

    /// The chord a level was built at.
    #[must_use]
    pub fn chord(&self, index: usize) -> Option<f64> {
        self.chords.get(index).copied()
    }

    /// The finest level: the one every descent ends at.
    #[must_use]
    pub fn finest(&self) -> &Pickable {
        &self.levels[self.levels.len() - 1]
    }

    /// The coarsest level.
    #[must_use]
    pub fn coarsest(&self) -> &Pickable {
        &self.levels[0]
    }

    /// The coarsest level no coarser than `chord` — what a view at that
    /// detail should be asked, without building anything new.
    #[must_use]
    pub fn for_chord(&self, chord: f64) -> &Pickable {
        for (index, held) in self.chords.iter().enumerate() {
            if *held <= chord {
                return &self.levels[index];
            }
        }
        self.finest()
    }

    /// Every hit along a ray, nearest first — the finest level's answer,
    /// reached by ruling faces out on the coarser ones first.
    #[must_use]
    pub fn pick(&self, ray: Ray, aperture: f64) -> Vec<Hit> {
        match self.narrowed(ray) {
            None => self.finest().pick(ray, aperture),
            Some(mask) => self.finest().pick_within(ray, aperture, Some(&mask)),
        }
    }

    /// The nearest hit along a ray, if any.
    #[must_use]
    pub fn pick_first(&self, ray: Ray, aperture: f64) -> Option<Hit> {
        self.pick(ray, aperture).into_iter().next()
    }

    /// As [`PickHierarchy::pick`], refined onto the exact surfaces.
    #[must_use]
    pub fn pick_refined(&self, ray: Ray, aperture: f64, tol: Tolerances) -> Vec<Hit> {
        // The refinement is the finest level's own, over the hits the descent
        // left: refining is per hit, so narrowing changes what is refined and
        // never how.
        let hits = self.pick(ray, aperture);
        self.finest().refine(hits, ray, tol)
    }

    /// The faces the coarse levels could not rule out, or `None` where there
    /// is nothing to rule out with.
    fn narrowed(&self, ray: Ray) -> Option<Vec<bool>> {
        if self.levels.len() < 2 {
            return None;
        }
        let fine = self.chords[self.chords.len() - 1];
        let mut mask: Option<Vec<bool>> = None;
        for (index, level) in self.levels[..self.levels.len() - 1].iter().enumerate() {
            let near = level.faces_near(ray, self.chords[index] + fine);
            mask = Some(match mask {
                None => near,
                Some(held) => held.into_iter().zip(near).map(|(a, b)| a && b).collect(),
            });
        }
        mask
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

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod refine_tests {
    use super::*;
    use ogeom_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    /// On a coarsely meshed cylinder the tessellated hit sits a visible
    /// sagitta off the true wall; the refined hit stands on the exact
    /// radius.
    #[test]
    fn refinement_lands_on_the_exact_surface() {
        let mut model = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 8.0, 10.0, T)
            .unwrap()
            .shape;
        let coarse = Deflection {
            chord: 0.5,
            ..Deflection::default()
        };
        let scene = Pickable::build(&model, &drum, coarse, T).unwrap();
        let ray = Ray {
            origin: Point::new(20.0, 1.7, 5.0),
            direction: Vector::new(-1.0, 0.0, 0.0),
        };
        let mesh_hit = scene.pick_first(ray, 0.0).unwrap();
        let refined = scene.pick_refined(ray, 0.0, T).into_iter().next().unwrap();
        assert!(refined.refined);
        let radius_of = |p: Point| p.x.hypot(p.y);
        assert!(
            (radius_of(refined.position) - 8.0).abs() < 1e-9,
            "refined radius {}",
            radius_of(refined.position)
        );
        assert!(
            (radius_of(mesh_hit.position) - 8.0).abs() > 1e-4,
            "the coarse mesh visibly sags: {}",
            radius_of(mesh_hit.position)
        );
        // Depth order and ownership survive refinement.
        assert_eq!(refined.triangle, mesh_hit.triangle);
        assert!(refined.distance < mesh_hit.distance + 1.0);
    }
}
