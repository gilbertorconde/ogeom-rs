//! Polygonal hidden line removal: project, classify, draw.
//!
//! The mesh does the occlusion work. The drawing's curves come from two
//! places: the model's own edges, discretized by the same machinery every
//! face boundary uses, and the tessellation's silhouettes, the mesh edges
//! where the surface turns away from the eye. Every sampled segment is
//! split where its projection crosses a contour of the mesh (a silhouette
//! or a free border), the only places its visibility can change, and each
//! piece is classified by casting its midpoint toward the eye against the
//! whole mesh: a triangle strictly in front hides it. Runs of
//! same-classified pieces merge back into polylines, so a curve that dips
//! behind a boss comes out as visible, hidden, visible: three curves, which
//! is what a drawing shows.
//!
//! Polygonal, not exact: the classification is as fine as the tessellation
//! and the sampling. The exact half (silhouettes in closed form, visibility
//! asked of the faces themselves) is the `exact` module.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Point2, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, Shape, ShapeType, Triangulation, explore};

use crate::crossings::{Contours, locate};

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
/// Segments that project to nothing (an edge running straight along the
/// view direction) are dropped: a point is not a line in a drawing.
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
    let occluders = Occluders::over(&mesh, view);

    // The contours: interior mesh edges whose triangles disagree about
    // facing the eye, border edges, and edges more than two triangles
    // share. Faces are meshed apart, so the mesh is welded first: otherwise
    // every face's border would read as one. The silhouettes drawn are the
    // contours that turn away from the eye and the borders of front faces.
    let welded = mesh.welded(tol);
    let toward_eye = view.toward_eye();
    let mut uses: ogeom_core::FastMap<(u32, u32), Vec<usize>> = ogeom_core::FastMap::default();
    for (t, triangle) in welded.triangles.iter().enumerate() {
        for i in 0..3 {
            let (a, b) = (triangle[i], triangle[(i + 1) % 3]);
            uses.entry((a.min(b), a.max(b))).or_default().push(t);
        }
    }
    let facing = |t: usize| -> f64 {
        let [a, b, c] = welded.triangles[t];
        let (pa, pb, pc) = (
            welded.positions[a as usize],
            welded.positions[b as usize],
            welded.positions[c as usize],
        );
        (pb - pa).cross(pc - pa).dot(toward_eye)
    };
    let mut edges: Vec<(&(u32, u32), &Vec<usize>)> = uses.iter().collect();
    edges.sort_by_key(|&(&(a, b), _)| (a, b));
    let mut contours: Vec<(Point2, Point2)> = Vec::new();
    let mut silhouettes: Vec<[Point; 2]> = Vec::new();
    for (&(a, b), triangles) in edges {
        let (contour, silhouette) = match triangles.as_slice() {
            [t] => (true, facing(*t) > 0.0),
            [s, t] => {
                let turns = (facing(*s) > 0.0) != (facing(*t) > 0.0);
                (turns, turns)
            }
            _ => (true, false),
        };
        let points = [welded.positions[a as usize], welded.positions[b as usize]];
        if contour {
            contours.push((view.project(points[0]), view.project(points[1])));
        }
        if silhouette {
            silhouettes.push(points);
        }
    }
    let contours = Contours::new(contours);
    let classifier = Classifier {
        view,
        occluders: &occluders,
        contours: &contours,
        clearance,
        tol,
    };

    let mut drawing = Drawing::default();

    // The model's own edges, their segments kept so a silhouette along one
    // is not drawn over it.
    let mut seen = ogeom_core::FastSet::default();
    let mut drawn = DrawnSegments::new(tol);
    for edge in explore(model, shape, Filter::OfType(ShapeType::Edge))? {
        let key = (edge.node(), edge.location().clone());
        if !seen.insert(key) {
            continue;
        }
        let Ok(points) = ogeom_mesh::polyline_of_edge(model, &edge, deflection, tol) else {
            continue;
        };
        drawn.add(&points);
        classifier.classify_into(&mut drawing, &points, Source::Edge(edge.clone()));
    }

    for points in silhouettes {
        if drawn.holds(points[0], points[1]) {
            continue;
        }
        classifier.classify_into(&mut drawing, &points, Source::Silhouette);
    }
    Ok(drawing)
}

/// A grid cell of [`DrawnSegments`].
type Cell = (i64, i64, i64);

/// The segments of the model edges already drawn, found by their ends.
struct DrawnSegments {
    cell: f64,
    tol: Tolerances,
    /// Each polyline point by its grid cell: which polyline, and where on
    /// it.
    points: ogeom_core::FastMap<Cell, Vec<(usize, usize, Point)>>,
    lines: usize,
}

impl DrawnSegments {
    fn new(tol: Tolerances) -> Self {
        Self {
            cell: tol.confusion().max(f64::MIN_POSITIVE) * 10.0,
            tol,
            points: ogeom_core::FastMap::default(),
            lines: 0,
        }
    }

    fn key(&self, p: Point) -> Cell {
        #[allow(clippy::cast_possible_truncation)]
        (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
            (p.z / self.cell).floor() as i64,
        )
    }

    fn add(&mut self, polyline: &[Point]) {
        for (index, p) in polyline.iter().enumerate() {
            let key = self.key(*p);
            self.points
                .entry(key)
                .or_default()
                .push((self.lines, index, *p));
        }
        self.lines += 1;
    }

    /// Where `p` stands on the drawn polylines.
    fn at(&self, p: Point) -> Vec<(usize, usize)> {
        let (x, y, z) = self.key(p);
        let mut out = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(found) = self.points.get(&(x + dx, y + dy, z + dz)) else {
                        continue;
                    };
                    out.extend(
                        found
                            .iter()
                            .filter(|(_, _, q)| q.is_equal(p, self.tol))
                            .map(|&(line, index, _)| (line, index)),
                    );
                }
            }
        }
        out
    }

    /// Whether `a` to `b` is a segment of a drawn polyline.
    fn holds(&self, a: Point, b: Point) -> bool {
        let at_b = self.at(b);
        self.at(a).into_iter().any(|(line, index)| {
            at_b.iter()
                .any(|&(other, j)| other == line && index.abs_diff(j) == 1)
        })
    }
}

/// What a polyline is classified against: the mesh that hides, and the
/// contours where hiding can begin or end.
struct Classifier<'a> {
    view: &'a View,
    occluders: &'a Occluders<'a>,
    contours: &'a Contours,
    clearance: f64,
    tol: Tolerances,
}

impl Classifier<'_> {
    /// Split a polyline into visible and hidden runs against the mesh.
    ///
    /// Each segment is cut where its projection crosses a contour, and each
    /// piece is classified by one ray from its middle: between two
    /// crossings nothing can pass in front of it.
    fn classify_into(&self, drawing: &mut Drawing, points: &[Point], source: Source) {
        let (view, tol) = (self.view, self.tol);
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
        let projected: Vec<Point2> = points.iter().map(|p| view.project(*p)).collect();
        let at = |position: f64| -> Point {
            let (k, f) = locate(position, points.len());
            points[k] + (points[k + 1] - points[k]) * f
        };
        for pair in self.contours.split(&projected, tol).windows(2) {
            let (from, to) = (pair[0], pair[1]);
            let (k, _) = locate(from, points.len());
            if projected[k].distance(projected[k + 1]) <= tol.confusion() {
                // Projects to a point: not a line in a drawing.
                flush(&mut run, run_visibility);
                run_visibility = None;
                continue;
            }
            let middle = at(f64::midpoint(from, to));
            let visibility = if self.occluders.occlude(middle, view, self.clearance) {
                Visibility::Hidden
            } else {
                Visibility::Visible
            };
            if run_visibility != Some(visibility) {
                flush(&mut run, run_visibility);
                run_visibility = Some(visibility);
            }
            if run.is_empty() {
                run.push(view.project(at(from)));
            }
            run.push(view.project(at(to)));
        }
        flush(&mut run, run_visibility);
    }
}

/// A mesh's triangles binned by where they project in the view.
///
/// The view is orthographic, so the ray from a point toward the eye can
/// only hit a triangle whose projection covers the point's: a sample is
/// tested against the triangles of its own cell of a grid over the drawing,
/// not the whole mesh.
struct Occluders<'m> {
    mesh: &'m Triangulation,
    cells: Vec<Vec<u32>>,
    low: Point2,
    size: f64,
    columns: usize,
    rows: usize,
}

impl<'m> Occluders<'m> {
    fn over(mesh: &'m Triangulation, view: &View) -> Self {
        let projected: Vec<Point2> = mesh.positions.iter().map(|p| view.project(*p)).collect();
        let (mut low, mut high) = (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for p in &projected {
            low = Point2::new(low.x.min(p.x), low.y.min(p.y));
            high = Point2::new(high.x.max(p.x), high.y.max(p.y));
        }
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let side = ((mesh.triangles.len() as f64).sqrt().ceil() as usize).clamp(1, 512);
        let span = (high.x - low.x).max(high.y - low.y);
        #[allow(clippy::cast_precision_loss)]
        let size = if span.is_finite() && span > 0.0 {
            span / side as f64
        } else {
            1.0
        };
        let (columns, rows) = (side, side);
        let mut cells: Vec<Vec<u32>> = vec![Vec::new(); columns * rows];
        let cell = |x: f64, lo: f64, n: usize| -> usize {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let k = ((x - lo) / size).floor().max(0.0) as usize;
            k.min(n - 1)
        };
        for (t, triangle) in mesh.triangles.iter().enumerate() {
            let corners = triangle.map(|i| projected[i as usize]);
            let (x0, x1) = (
                corners.iter().map(|p| p.x).fold(f64::INFINITY, f64::min),
                corners
                    .iter()
                    .map(|p| p.x)
                    .fold(f64::NEG_INFINITY, f64::max),
            );
            let (y0, y1) = (
                corners.iter().map(|p| p.y).fold(f64::INFINITY, f64::min),
                corners
                    .iter()
                    .map(|p| p.y)
                    .fold(f64::NEG_INFINITY, f64::max),
            );
            if !(x0.is_finite() && x1.is_finite() && y0.is_finite() && y1.is_finite()) {
                continue;
            }
            // A cell's worth of margin each way: the test's own tolerance
            // lets a hit land a hair outside the exact projection.
            let (c0, c1) = (
                cell(x0 - size, low.x, columns),
                cell(x1 + size, low.x, columns),
            );
            let (r0, r1) = (cell(y0 - size, low.y, rows), cell(y1 + size, low.y, rows));
            #[allow(clippy::cast_possible_truncation)]
            for r in r0..=r1 {
                for c in c0..=c1 {
                    cells[r * columns + c].push(t as u32);
                }
            }
        }
        Self {
            mesh,
            cells,
            low,
            size,
            columns,
            rows,
        }
    }

    /// Whether anything in the mesh stands between `p` and the eye.
    fn occlude(&self, p: Point, view: &View, clearance: f64) -> bool {
        let q = view.project(p);
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let at = |x: f64, lo: f64, n: usize| -> Option<usize> {
            // A point on the drawing's outer edge falls exactly at the
            // grid's end: within a cell of it, it belongs to the last cell.
            let k = ((x - lo) / self.size).floor();
            (k >= -1.0 && k <= n as f64).then(|| (k.max(0.0) as usize).min(n - 1))
        };
        let (Some(c), Some(r)) = (
            at(q.x, self.low.x, self.columns),
            at(q.y, self.low.y, self.rows),
        ) else {
            return false;
        };
        occluded(
            self.mesh,
            self.cells[r * self.columns + c].iter().map(|t| *t as usize),
            p,
            view,
            clearance,
        )
    }
}

/// Whether any of the given triangles stands between the point and the eye.
fn occluded(
    mesh: &Triangulation,
    candidates: impl Iterator<Item = usize>,
    p: Point,
    view: &View,
    clearance: f64,
) -> bool {
    let toward_eye = view.toward_eye();
    let depth = view.depth(p);
    for t in candidates {
        let [a, b, c] = mesh.triangles[t];
        let (pa, pb, pc) = (
            mesh.positions[a as usize],
            mesh.positions[b as usize],
            mesh.positions[c as usize],
        );
        // Möller-Trumbore, orthographic: the ray from p toward the eye.
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
        // face. The four vertical edges project to points and are dropped.
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

    fn silhouette_length(drawing: &Drawing, visibility: Visibility) -> f64 {
        drawing
            .curves()
            .filter(|c| c.visibility == visibility && matches!(c.source, Source::Silhouette))
            .map(|c| {
                c.points
                    .windows(2)
                    .map(|w| w[0].distance(w[1]))
                    .sum::<f64>()
            })
            .sum()
    }

    /// Each outline is drawn once: a box's outline is its own edges, with
    /// no silhouette over them, and a drum's side silhouette is its two
    /// generators, not its face borders as well.
    #[test]
    fn every_outline_is_drawn_once() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, MFrame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        let view =
            View::looking(Vector::new(-1.0, -1.2, -0.9), Vector::new(0.0, 0.0, 1.0), T).unwrap();
        let drawing = super::project(&model, &solid.shape, &view, fine(), T).unwrap();
        assert!(silhouette_length(&drawing, Visibility::Visible) < 1e-9);
        assert!(silhouette_length(&drawing, Visibility::Hidden) < 1e-9);

        let drum = ogeom_algo::make_cylinder(&mut model, MFrame::WORLD, 3.0, 8.0, T).unwrap();
        let side =
            View::looking(Vector::new(-1.0, 0.3, 0.0), Vector::new(0.0, 0.0, 1.0), T).unwrap();
        let drawing = super::project(&model, &drum.shape, &side, fine(), T).unwrap();
        let visible = silhouette_length(&drawing, Visibility::Visible);
        assert!((visible - 16.0).abs() < 1e-6, "two generators: {visible}");
    }

    /// A bar under a block that covers its far end, seen from above: the
    /// bar's top edge is one straight segment, visible up to the block's
    /// side and hidden past it.
    #[test]
    fn a_straight_edge_half_under_a_block_is_split_at_the_blocks_outline() {
        let mut model = Model::new();
        let bar = ogeom_algo::make_box(&mut model, MFrame::WORLD, (20.0, 1.0, 1.0), T).unwrap();
        let at = MFrame::new(Point::new(12.0, -5.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let block = ogeom_algo::make_box(&mut model, at, (20.0, 10.0, 2.0), T).unwrap();
        let both = ogeom_algo::build::make_compound(&mut model, &[bar.shape.clone(), block.shape])
            .unwrap();
        let view =
            View::looking(Vector::new(0.0, 0.0, -1.0), Vector::new(0.0, 1.0, 0.0), T).unwrap();
        let drawing = super::project(&model, &both.shape, &view, fine(), T).unwrap();

        let top = explore(&model, &bar.shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .into_iter()
            .find(|e| {
                let points = ogeom_mesh::polyline_of_edge(&model, e, fine(), T).unwrap();
                points.len() == 2
                    && points
                        .iter()
                        .all(|p| p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9)
            })
            .unwrap();
        let spans = |curves: &[DrawnCurve]| -> Vec<(f64, f64)> {
            curves
                .iter()
                .filter(|c| matches!(&c.source, Source::Edge(e) if e.node() == top.node()))
                .map(|c| {
                    let xs = c.points.iter().map(|p| p.x);
                    (
                        xs.clone().fold(f64::INFINITY, f64::min),
                        xs.fold(f64::NEG_INFINITY, f64::max),
                    )
                })
                .collect()
        };
        let visible = spans(&drawing.visible);
        let hidden = spans(&drawing.hidden);
        assert_eq!(visible.len(), 1, "{visible:?} {hidden:?}");
        assert_eq!(hidden.len(), 1, "{hidden:?}");
        assert!(visible[0].0.abs() < 1e-9 && (visible[0].1 - 12.0).abs() < 1e-9);
        assert!((hidden[0].0 - 12.0).abs() < 1e-9 && (hidden[0].1 - 20.0).abs() < 1e-9);
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
