//! Polygonal hidden line removal: project, classify, draw.
//!
//! The mesh does the occlusion work. The drawing's curves come from two
//! places — the model's own edges, discretized by the same machinery every
//! face boundary uses, and the tessellation's silhouettes, the mesh edges
//! where the surface turns away from the eye. Every sampled segment is
//! classified by casting its midpoint toward the eye against the whole
//! mesh: a triangle strictly in front hides it. Runs of same-classified
//! segments merge back into polylines, so a curve that dips behind a boss
//! comes out as visible, hidden, visible — three curves, which is what a
//! drawing shows.
//!
//! Polygonal, not exact: the classification is as fine as the tessellation
//! and the sampling. That is the honest half of `HLRBRep`; the exact half —
//! curve/surface interference resolved analytically — is deferred and the
//! deferred table says so.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Point2, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, Shape, ShapeType, Triangulation, explore};

/// Which side of the pencil a curve lands on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Visibility {
    /// Nothing stands between the curve and the eye.
    Visible,
    /// Something does: drawn dashed, or not at all.
    Hidden,
}

/// Where a drawn curve came from.
#[derive(Debug, Clone)]
pub enum Source {
    /// A model edge, with the occurrence that produced it.
    Edge(Shape),
    /// A silhouette: the tessellation turning away from the eye.
    Silhouette,
}

/// One polyline of the drawing, in view-plane coordinates.
#[derive(Debug, Clone)]
pub struct DrawnCurve {
    /// The projected points, in order.
    pub points: Vec<Point2>,
    /// Visible or hidden.
    pub visibility: Visibility,
    /// What it is a picture of.
    pub source: Source,
}

/// A 2D drawing: the classified projection of a shape.
#[derive(Debug, Clone, Default)]
pub struct Drawing {
    /// Curves nothing occludes.
    pub visible: Vec<DrawnCurve>,
    /// Curves something does.
    pub hidden: Vec<DrawnCurve>,
}

impl Drawing {
    /// Every curve, visible first.
    pub fn curves(&self) -> impl Iterator<Item = &DrawnCurve> {
        self.visible.iter().chain(self.hidden.iter())
    }
}

/// The view for a drawing: an orthographic camera looking along `-z` of the
/// frame it carries, with `x` right and `y` up on the sheet.
#[derive(Debug, Clone, Copy)]
pub struct View {
    frame: Frame,
}

impl View {
    /// A view looking along `direction`, with `up` steadying the sheet.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `up` is parallel to the view direction.
    pub fn looking(direction: Vector, up: Vector, tol: Tolerances) -> OgeomResult<Self> {
        let toward_eye = Direction::new(-direction, tol)?;
        let right = Direction::new(up.cross(toward_eye.vector()), tol)?;
        Ok(Self {
            frame: Frame::new(Point::ORIGIN, toward_eye, right, tol)?,
        })
    }

    /// Sheet coordinates of a world point: `x` right, `y` up.
    #[must_use]
    pub fn project(&self, p: Point) -> Point2 {
        let local = self.frame.to_local(p);
        Point2::new(local.x, local.y)
    }

    /// Depth of a world point: greater is nearer the eye.
    #[must_use]
    pub fn depth(&self, p: Point) -> f64 {
        self.frame.to_local(p).z
    }

    /// The world direction toward the eye.
    #[must_use]
    pub fn toward_eye(&self) -> Vector {
        self.frame.z().vector()
    }
}

/// Project a shape into a classified 2D drawing.
///
/// The model's edges and the tessellation's silhouettes, each split into
/// visible and hidden runs by occlusion against the shape's own mesh.
/// Segments that project to nothing — an edge running straight along the
/// view direction — are dropped: a point is not a line in a drawing.
///
/// # Errors
///
/// As [`ogeom_mesh::triangulate()`]; and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape has no faces to draw.
pub fn project(
    model: &Model,
    shape: &Shape,
    view: &View,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Drawing> {
    let mesh = ogeom_mesh::triangulate(model, shape, deflection, tol)?;
    if mesh.is_empty() {
        ogeom_bail!(Construction, "the shape tessellates to nothing to draw");
    }
    // A segment on the exact surface sits up to a sagitta outside the
    // inscribed mesh, and a sample on a front face must not be occluded by
    // that face's own triangles: the clearance a hit must beat.
    let clearance = deflection.chord.max(tol.confusion() * 1e3) * 4.0;

    let mut drawing = Drawing::default();

    // The model's own edges.
    let mut seen = std::collections::HashSet::new();
    for edge in explore(model, shape, Filter::OfType(ShapeType::Edge))? {
        let key = (edge.node(), edge.location().clone());
        if !seen.insert(key) {
            continue;
        }
        let Ok(points) = ogeom_mesh::polyline_of_edge(model, &edge, deflection, tol) else {
            continue;
        };
        classify_into(
            &mut drawing,
            &points,
            Source::Edge(edge.clone()),
            view,
            &mesh,
            clearance,
            tol,
        );
    }

    // Silhouettes: interior mesh edges whose triangles disagree about facing
    // the eye, and border edges, which are their own outline.
    let toward_eye = view.toward_eye();
    let mut uses: std::collections::HashMap<(u32, u32), Vec<usize>> =
        std::collections::HashMap::new();
    for (t, triangle) in mesh.triangles.iter().enumerate() {
        for i in 0..3 {
            let (a, b) = (triangle[i], triangle[(i + 1) % 3]);
            uses.entry((a.min(b), a.max(b))).or_default().push(t);
        }
    }
    let facing = |t: usize| -> f64 {
        let [a, b, c] = mesh.triangles[t];
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        (pb - pa).cross(pc - pa).dot(toward_eye)
    };
    let mut edges: Vec<(&(u32, u32), &Vec<usize>)> = uses.iter().collect();
    edges.sort_by_key(|&(&(a, b), _)| (a, b));
    for (&(a, b), triangles) in edges {
        let silhouette = match triangles.as_slice() {
            [t] => facing(*t) > 0.0,
            [s, t] => (facing(*s) > 0.0) != (facing(*t) > 0.0),
            _ => false,
        };
        if !silhouette {
            continue;
        }
        let points = [mesh.positions[a as usize], mesh.positions[b as usize]];
        classify_into(
            &mut drawing,
            &points,
            Source::Silhouette,
            view,
            &mesh,
            clearance,
            tol,
        );
    }
    Ok(drawing)
}

/// Split a polyline into visible and hidden runs against the mesh.
fn classify_into(
    drawing: &mut Drawing,
    points: &[Point],
    source: Source,
    view: &View,
    mesh: &Triangulation,
    clearance: f64,
    tol: Tolerances,
) {
    let mut run: Vec<Point2> = Vec::new();
    let mut run_visibility: Option<Visibility> = None;
    let mut flush = |run: &mut Vec<Point2>, visibility: Option<Visibility>| {
        if run.len() < 2 {
            run.clear();
            return;
        }
        let curve = DrawnCurve {
            points: std::mem::take(run),
            visibility: visibility.unwrap_or(Visibility::Visible),
            source: source.clone(),
        };
        match visibility {
            Some(Visibility::Hidden) => drawing.hidden.push(curve),
            _ => drawing.visible.push(curve),
        }
    };
    for pair in points.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let (pa, pb) = (view.project(a), view.project(b));
        if pa.distance(pb) <= tol.confusion() {
            // Projects to a point: not a line in a drawing.
            flush(&mut run, run_visibility);
            run_visibility = None;
            continue;
        }
        let mid = Point::new(
            f64::midpoint(a.x, b.x),
            f64::midpoint(a.y, b.y),
            f64::midpoint(a.z, b.z),
        );
        let visibility = if occluded(mesh, mid, view, clearance) {
            Visibility::Hidden
        } else {
            Visibility::Visible
        };
        if run_visibility != Some(visibility) {
            flush(&mut run, run_visibility);
            run_visibility = Some(visibility);
        }
        if run.is_empty() {
            run.push(pa);
        }
        run.push(pb);
    }
    flush(&mut run, run_visibility);
}

/// Whether anything in the mesh stands between the point and the eye.
fn occluded(mesh: &Triangulation, p: Point, view: &View, clearance: f64) -> bool {
    let toward_eye = view.toward_eye();
    let depth = view.depth(p);
    for triangle in &mesh.triangles {
        let [a, b, c] = *triangle;
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        // Möller–Trumbore, orthographic: the ray from p toward the eye.
        let (e1, e2) = (pb - pa, pc - pa);
        let h = toward_eye.cross(e2);
        let det = e1.dot(h);
        if det.abs() < 1e-14 {
            continue;
        }
        let inv = 1.0 / det;
        let s = p - pa;
        let u = s.dot(h) * inv;
        if !(0.0..=1.0).contains(&u) {
            continue;
        }
        let q = s.cross(e1);
        let v = toward_eye.dot(q) * inv;
        if v < 0.0 || u + v > 1.0 {
            continue;
        }
        let t = e2.dot(q) * inv;
        if t <= clearance {
            continue;
        }
        let hit_depth = depth + t;
        if hit_depth > depth + clearance {
            return true;
        }
    }
    false
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_math::Frame as MFrame;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-2,
            ..Deflection::default()
        }
    }

    fn edge_curves(drawing: &Drawing, visibility: Visibility) -> usize {
        drawing
            .curves()
            .filter(|c| c.visibility == visibility && matches!(c.source, Source::Edge(_)))
            .count()
    }

    #[test]
    fn a_box_face_on_shows_its_front_and_hides_its_back() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, MFrame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        // Looking down -z: the eye is above, the top face is the front.
        let view =
            View::looking(Vector::new(0.0, 0.0, -1.0), Vector::new(0.0, 1.0, 0.0), T).unwrap();
        let drawing = super::project(&model, &solid.shape, &view, fine(), T).unwrap();

        // Four top edges visible, four bottom edges hidden behind the top
        // face; the four vertical edges project to points and are dropped.
        assert_eq!(edge_curves(&drawing, Visibility::Visible), 4);
        assert_eq!(edge_curves(&drawing, Visibility::Hidden), 4);
    }

    #[test]
    fn a_box_in_three_quarter_view_shows_nine_and_hides_three() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, MFrame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        // The classic drawing-class view: three faces show, nine edges
        // visible, the three edges meeting at the far corner hidden.
        let view =
            View::looking(Vector::new(-1.0, -1.2, -0.9), Vector::new(0.0, 0.0, 1.0), T).unwrap();
        let drawing = super::project(&model, &solid.shape, &view, fine(), T).unwrap();
        assert_eq!(edge_curves(&drawing, Visibility::Visible), 9);
        assert_eq!(edge_curves(&drawing, Visibility::Hidden), 3);
    }

    #[test]
    fn a_cylinder_from_the_side_draws_its_silhouette() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_cylinder(&mut model, MFrame::WORLD, 3.0, 8.0, T).unwrap();
        let view =
            View::looking(Vector::new(-1.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0), T).unwrap();
        let drawing = super::project(&model, &solid.shape, &view, fine(), T).unwrap();

        // Silhouette runs exist, and the visible drawing spans the
        // cylinder's height and diameter.
        let silhouettes = drawing
            .visible
            .iter()
            .filter(|c| matches!(c.source, Source::Silhouette))
            .count();
        assert!(silhouettes > 0, "a curved side draws by its silhouette");
        let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
        for curve in &drawing.visible {
            for p in &curve.points {
                min_x = min_x.min(p.x);
                max_x = max_x.max(p.x);
                min_y = min_y.min(p.y);
                max_y = max_y.max(p.y);
            }
        }
        assert!(
            (max_x - min_x - 6.0).abs() < 0.1,
            "diameter across the sheet"
        );
        assert!((max_y - min_y - 8.0).abs() < 0.1, "height up the sheet");
    }

    #[test]
    fn a_small_box_behind_a_large_one_is_entirely_hidden() {
        let mut model = Model::new();
        let front = ogeom_algo::make_box(&mut model, MFrame::WORLD, (20.0, 20.0, 2.0), T).unwrap();
        let behind_frame =
            MFrame::new(Point::new(8.0, 8.0, -10.0), Direction::Z, Direction::X, T).unwrap();
        let back = ogeom_algo::make_box(&mut model, behind_frame, (4.0, 4.0, 2.0), T).unwrap();
        let both = ogeom_algo::build::make_compound(
            &mut model,
            &[front.shape.clone(), back.shape.clone()],
        )
        .unwrap();
        let view =
            View::looking(Vector::new(0.0, 0.0, -1.0), Vector::new(0.0, 1.0, 0.0), T).unwrap();
        let drawing = super::project(&model, &both.shape, &view, fine(), T).unwrap();

        // Every drawable edge of the back box is hidden by the front plate.
        let back_edges_visible = drawing
            .visible
            .iter()
            .filter_map(|c| match &c.source {
                Source::Edge(e) => Some(e),
                Source::Silhouette => None,
            })
            .filter(|e| {
                ogeom_topo::explore(&model, &back.shape, Filter::OfType(ShapeType::Edge))
                    .unwrap()
                    .iter()
                    .any(|be| be.node() == e.node())
            })
            .count();
        assert_eq!(back_edges_visible, 0, "the plate hides the block");
    }
}
