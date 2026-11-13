//! Proving that a face lies within a region bounded by planes.
//!
//! A [`HalfSpaces`] is a tree of half-spaces under maxima and minima: a convex
//! corner is the maximum of its faces' heights, a concave edge the minimum
//! of its two, and so on. [`certify_face_within`] proves that every point
//! of a face stands no farther past the region than a stated fit, or finds
//! a point that does.
//!
//! The proof rests on the convex hull. Every surface kind converts exactly
//! to a B-spline over any rectangle of its own parameters
//! ([`SurfaceGeometry::to_bspline`] of a [`TrimmedSurface`]), and a B-spline
//! whose weights are all positive lies in the convex hull of its control
//! points. A plane's height is affine, so its greatest over the hull is its
//! greatest at a control point. The maximum of several pieces is bounded by
//! the maximum of their bounds, and the minimum by the minimum of theirs. A
//! rectangle whose bound is within the fit is cleared; one wholly outside
//! the face's trim is passed over; the rest is split in four and tried
//! again, down to [`DEPTH`] halvings each way and [`BUDGET`] rectangles in
//! all.
//!
//! The trim is read from the face's rings in its chart, each drawn within a
//! chord of its edges in space through the surface's tangent scales
//! ([`ogeom_mesh::face_boundary_in_chart_met`]). A rectangle is passed over
//! only where its centre lies outside the rings and every ring stands
//! farther from that centre than the rectangle's half diagonal and a slack
//! of [`SLACK`] chords read in the chart at that centre. That reading of
//! the trim is the one input the proof takes from the mesher.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{BSplineSurface, Surface as _, SurfaceGeometry, TrimmedSurface};
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{Model, Shape};

/// A region of space bounded by planes, as a tree of half-spaces.
#[derive(Debug, Clone, PartialEq)]
pub enum HalfSpaces {
    /// The points on the side of the plane through `at` that `normal`
    /// points away from: a point's height is `normal · (p - at)`, with
    /// `normal` of unit length.
    Below {
        /// The plane's outward unit normal.
        normal: Vector,
        /// A point on the plane.
        at: Point,
    },
    /// The points inside every piece: a height the greatest of the
    /// pieces'.
    Every(Vec<HalfSpaces>),
    /// The points inside any piece: a height the least of the pieces'.
    Any(Vec<HalfSpaces>),
}

impl HalfSpaces {
    /// How far `p` stands past the region (negative inside), and the
    /// outward normal of the plane that decides it.
    #[must_use]
    pub fn height(&self, p: Point) -> (f64, Vector) {
        match self {
            Self::Below { normal, at } => (normal.dot(p - *at), *normal),
            Self::Every(pieces) => pieces.iter().map(|r| r.height(p)).fold(
                (f64::NEG_INFINITY, Vector::ZERO),
                |a, b| {
                    if b.0 > a.0 { b } else { a }
                },
            ),
            Self::Any(pieces) => {
                pieces
                    .iter()
                    .map(|r| r.height(p))
                    .fold(
                        (f64::INFINITY, Vector::ZERO),
                        |a, b| {
                            if b.0 < a.0 { b } else { a }
                        },
                    )
            }
        }
    }

    /// A bound on the height over the convex hull of `hull`: exact for a
    /// plane and for the greatest of pieces, at least the least of pieces'
    /// heights for the least.
    #[must_use]
    pub fn bound_over(&self, hull: &[Point]) -> f64 {
        match self {
            Self::Below { normal, at } => hull
                .iter()
                .map(|p| normal.dot(*p - *at))
                .fold(f64::NEG_INFINITY, f64::max),
            Self::Every(pieces) => pieces
                .iter()
                .map(|r| r.bound_over(hull))
                .fold(f64::NEG_INFINITY, f64::max),
            Self::Any(pieces) => pieces
                .iter()
                .map(|r| r.bound_over(hull))
                .fold(f64::INFINITY, f64::min),
        }
    }
}

/// What kept a face from being proved within its region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Uncertified {
    /// How far past the region: the height at `at` where `found`, else the
    /// bound over the rectangle the proof could not clear.
    pub past: f64,
    /// A point of the face standing `past` outside, or the centre of the
    /// rectangle the proof could not clear.
    pub at: Point,
    /// Whether `at` is a point of the face found standing past the fit, as
    /// against one the proof ran out of halvings over.
    pub found: bool,
}

/// How many times the face's parameter rectangle is halved each way before
/// a part of it no bound clears is taken as standing outside.
pub const DEPTH: u32 = 20;

/// How many rectangles a proof may try before what it has not cleared is
/// taken as standing outside.
pub const BUDGET: usize = 1 << 20;

/// The chords of the trim's drawing allowed between a rectangle and the
/// rings before it is passed over as outside the trim.
pub const SLACK: f64 = 2.0;

/// Whether every point of `face` stands within `fit` of `region`: `Ok(())`
/// when proved, else what kept it from being proved (see the module
/// documentation for the argument and its one input from the mesher).
///
/// # Errors
///
/// As reading the face's surface and its rings, or converting a rectangle
/// of the surface to a B-spline.
pub fn certify_face_within(
    model: &Model,
    face: &Shape,
    region: &HalfSpaces,
    fit: f64,
    tol: Tolerances,
) -> OgeomResult<Result<(), Uncertified>> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_core::ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let place = face.transform(model.datums())?;
    let chord = fit / 4.0;
    let deflection = ogeom_mesh::Deflection {
        chord,
        ..ogeom_mesh::Deflection::default()
    };
    let Some(rings) = ogeom_mesh::face_boundary_in_chart_met(model, face, deflection, tol)? else {
        // A trim drawn short of its chord says nothing the proof can use.
        let ((u, _), (v, _)) = surface.domain();
        return Ok(Err(Uncertified {
            past: f64::INFINITY,
            at: place.apply(surface.point_at(u, v, tol)?),
            found: false,
        }));
    };
    let Some(trim) = Trim::of(&rings) else {
        return Ok(Ok(()));
    };
    // A B-spline's pieces are cut from their parent's, keeping its
    // parameters; any other surface's converted afresh over each rectangle.
    let spline = match &surface {
        SurfaceGeometry::BSpline(s) => {
            let ((u0, u1), (v0, v1)) = s.domain();
            let (lo, hi) = (trim.lo, trim.hi);
            Some(s.segment(
                (lo.x.max(u0), hi.x.min(u1)),
                (lo.y.max(v0), hi.y.min(v1)),
                tol,
            )?)
        }
        _ => None,
    };
    let start = match &spline {
        Some(s) => {
            let ((u0, u1), (v0, v1)) = s.domain();
            (Point2::new(u0, v0), Point2::new(u1, v1))
        }
        None => (trim.lo, trim.hi),
    };
    let mut pending: Vec<(Point2, Point2, u32, Option<BSplineSurface>)> =
        vec![(start.0, start.1, 0, spline)];
    let mut tried = 0usize;
    while let Some((a, b, depth, parent)) = pending.pop() {
        tried += 1;
        let centre = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
        let inside = trim.inside(centre);
        if !inside {
            // The chart's slack at the centre: chords of space read back
            // through the tangent scales there.
            let (du, dv) = surface.d1_at(centre.x, centre.y, tol)?;
            let scale = du.magnitude().min(dv.magnitude());
            let slack = if scale > 0.0 {
                SLACK * chord / scale
            } else {
                f64::INFINITY
            };
            let half = (b.x - a.x).hypot(b.y - a.y) * 0.5;
            if !trim.near(centre, half + slack) {
                continue;
            }
        }
        let piece = match &parent {
            Some(p) => p.segment((a.x, b.x), (a.y, b.y), tol)?,
            None => SurfaceGeometry::Trimmed(Box::new(TrimmedSurface::new(
                surface.clone(),
                (a.x, b.x),
                (a.y, b.y),
                tol,
            )?))
            .to_bspline(tol)?,
        };
        let net = piece.grid().points();
        let bound = if net.iter().all(|w| w.weight > 0.0) {
            let hull: Vec<Point> = net.iter().map(|w| place.apply(w.point())).collect();
            region.bound_over(&hull)
        } else {
            f64::INFINITY
        };
        if bound <= fit {
            continue;
        }
        let at = place.apply(surface.point_at(centre.x, centre.y, tol)?);
        if inside {
            let (past, _) = region.height(at);
            if past > fit {
                return Ok(Err(Uncertified {
                    past,
                    at,
                    found: true,
                }));
            }
        }
        if depth >= DEPTH || tried >= BUDGET {
            return Ok(Err(Uncertified {
                past: bound,
                at,
                found: false,
            }));
        }
        let carried = parent.is_some().then_some(piece);
        for (x0, x1) in [(a.x, centre.x), (centre.x, b.x)] {
            for (y0, y1) in [(a.y, centre.y), (centre.y, b.y)] {
                pending.push((
                    Point2::new(x0, y0),
                    Point2::new(x1, y1),
                    depth + 1,
                    carried.clone(),
                ));
            }
        }
    }
    Ok(Ok(()))
}

/// A face's trim rings in its chart, bucketed: their segments on a grid of
/// cells for the nearest within a reach, and by bands of height for the
/// crossings a ray to `+x` makes.
struct Trim {
    segments: Vec<(Point2, Point2)>,
    lo: Point2,
    hi: Point2,
    step: (f64, f64),
    cells: Vec<Vec<usize>>,
    bands: Vec<Vec<usize>>,
}

/// The cells of a [`Trim`] each way.
const CELLS: usize = 64;

impl Trim {
    /// The rings bucketed; `None` where they hold no segment.
    fn of(rings: &[Vec<Point2>]) -> Option<Self> {
        let segments: Vec<(Point2, Point2)> = rings
            .iter()
            .flat_map(|ring| (0..ring.len()).map(move |k| (ring[k], ring[(k + 1) % ring.len()])))
            .collect();
        if segments.is_empty() {
            return None;
        }
        let (mut lo, mut hi) = (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for (a, _) in &segments {
            lo = Point2::new(lo.x.min(a.x), lo.y.min(a.y));
            hi = Point2::new(hi.x.max(a.x), hi.y.max(a.y));
        }
        #[expect(clippy::cast_precision_loss, reason = "a cell count of sixty-four")]
        let step = (
            ((hi.x - lo.x) / CELLS as f64).max(f64::MIN_POSITIVE),
            ((hi.y - lo.y) / CELLS as f64).max(f64::MIN_POSITIVE),
        );
        let mut trim = Self {
            segments,
            lo,
            hi,
            step,
            cells: vec![Vec::new(); CELLS * CELLS],
            bands: vec![Vec::new(); CELLS],
        };
        for k in 0..trim.segments.len() {
            let (a, b) = trim.segments[k];
            let (i0, j0) = trim.cell(Point2::new(a.x.min(b.x), a.y.min(b.y)));
            let (i1, j1) = trim.cell(Point2::new(a.x.max(b.x), a.y.max(b.y)));
            for i in i0..=i1 {
                for j in j0..=j1 {
                    trim.cells[i * CELLS + j].push(k);
                }
            }
            for band in &mut trim.bands[j0..=j1] {
                band.push(k);
            }
        }
        Some(trim)
    }

    /// The cell holding `p`, clamped to the grid.
    fn cell(&self, p: Point2) -> (usize, usize) {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "clamped to the grid's sixty-four cells"
        )]
        let at = |x: f64, lo: f64, step: f64| {
            (((x - lo) / step).floor().max(0.0) as usize).min(CELLS - 1)
        };
        (
            at(p.x, self.lo.x, self.step.0),
            at(p.y, self.lo.y, self.step.1),
        )
    }

    /// Whether a segment of the rings comes within `reach` of `p`.
    fn near(&self, p: Point2, reach: f64) -> bool {
        if !reach.is_finite() {
            return true;
        }
        if p.x + reach < self.lo.x
            || p.x - reach > self.hi.x
            || p.y + reach < self.lo.y
            || p.y - reach > self.hi.y
        {
            return false;
        }
        let (i0, j0) = self.cell(Point2::new(p.x - reach, p.y - reach));
        let (i1, j1) = self.cell(Point2::new(p.x + reach, p.y + reach));
        (i0..=i1).any(|i| {
            (j0..=j1).any(|j| {
                self.cells[i * CELLS + j]
                    .iter()
                    .any(|&k| distance_to(p, self.segments[k]) <= reach)
            })
        })
    }

    /// Whether `p` lies inside the rings, by the parity of the crossings
    /// its ray to `+x` makes: only segments of its band can cross it.
    fn inside(&self, p: Point2) -> bool {
        if p.y < self.lo.y || p.y > self.hi.y {
            return false;
        }
        let (_, j) = self.cell(p);
        let mut crossings = 0usize;
        for &k in &self.bands[j] {
            let (a, b) = self.segments[k];
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > p.x {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    }
}

/// The distance from `p` to the segment `(a, b)`.
fn distance_to(p: Point2, (a, b): (Point2, Point2)) -> f64 {
    let (ux, uy) = (b.x - a.x, b.y - a.y);
    let length = ux.mul_add(ux, uy * uy);
    let t = if length > 0.0 {
        (((p.x - a.x) * ux + (p.y - a.y) * uy) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p.x - (a.x + ux * t)).hypot(p.y - (a.y + uy * t))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_geom::BSplineSurface;
    use ogeom_math::{ControlGrid, Frame, KnotVector};
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();
    const FIT: f64 = 2e-4;

    fn below(normal: Vector, at: Point) -> HalfSpaces {
        HalfSpaces::Below { normal, at }
    }

    /// The square prism `|x|, |y| <= half`, as the maximum of its four sides.
    fn square(half: f64) -> HalfSpaces {
        let side =
            |x: f64, y: f64| below(Vector::new(x, y, 0.0), Point::new(x * half, y * half, 0.0));
        HalfSpaces::Every(vec![
            side(1.0, 0.0),
            side(-1.0, 0.0),
            side(0.0, 1.0),
            side(0.0, -1.0),
        ])
    }

    /// The faces of `solid` on a surface of `kind`.
    fn faces_of(model: &Model, solid: &Shape, kind: ogeom_geom::SurfaceKind) -> Vec<Shape> {
        explore_unique(model, solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                let data = model.node(f).unwrap().data().as_face().unwrap();
                model.geometry().surface(data.surface).unwrap().kind() == kind
            })
            .collect()
    }

    /// A unit cylinder's side touches the square prism round it along four
    /// lines and lies inside it: proved, the tangency cleared by the arcs'
    /// middle control points standing on the sides. Drawn in by a tenth,
    /// the sides cut it, and a point standing out is found.
    #[test]
    fn a_cylinder_inside_its_square_is_proved_and_one_poking_out_is_found() {
        let mut m = Model::new();
        let solid = crate::make_cylinder(&mut m, Frame::WORLD, 1.0, 3.0, T)
            .unwrap()
            .shape;
        let side = faces_of(&m, &solid, ogeom_geom::SurfaceKind::Cylinder);
        assert_eq!(side.len(), 1);
        assert_eq!(
            certify_face_within(&m, &side[0], &square(1.0), FIT, T).unwrap(),
            Ok(())
        );
        let out = certify_face_within(&m, &side[0], &square(0.9), FIT, T)
            .unwrap()
            .unwrap_err();
        assert!(out.found, "{out:?}");
        assert!(out.past > FIT && out.past <= 0.1 + 1e-12, "{out:?}");
    }

    /// A sphere of radius 2 inside the cube it touches, the cube as the
    /// maximum of six planes: proved. With one side drawn in by 0.05, the
    /// point found stands past it, by no more than 0.05.
    #[test]
    fn a_sphere_in_its_cube_is_proved_and_past_a_drawn_in_side_is_found() {
        let mut m = Model::new();
        let solid = crate::make_sphere(&mut m, Frame::WORLD, 2.0, T)
            .unwrap()
            .shape;
        let faces = faces_of(&m, &solid, ogeom_geom::SurfaceKind::Sphere);
        let side = |n: Vector, d: f64| below(n, Point::ORIGIN + n * d);
        let axes = [
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(-1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, -1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
            Vector::new(0.0, 0.0, -1.0),
        ];
        let cube = HalfSpaces::Every(axes.iter().map(|n| side(*n, 2.0)).collect());
        for f in &faces {
            assert_eq!(certify_face_within(&m, f, &cube, FIT, T).unwrap(), Ok(()));
        }
        let mut pieces: Vec<HalfSpaces> = axes.iter().map(|n| side(*n, 2.0)).collect();
        pieces[4] = side(axes[4], 1.95);
        let drawn = HalfSpaces::Every(pieces);
        let worst = faces
            .iter()
            .filter_map(|f| certify_face_within(&m, f, &drawn, FIT, T).unwrap().err())
            .collect::<Vec<_>>();
        assert!(!worst.is_empty());
        assert!(
            worst
                .iter()
                .all(|u| u.found && u.past > FIT && u.past <= 0.05 + 1e-12)
        );
    }

    /// The soundness hole a guard reading only vertices leaves: a
    /// biquadratic patch over the unit square whose corners and edges lie
    /// on the plane `z = 0` but whose middle control point stands 1 above
    /// it. Every vertex is on the plane; the patch stands a quarter past it
    /// at its middle (the middle basis function is a half squared there),
    /// and that is what is found, to the bit.
    #[test]
    fn a_patch_bulging_past_its_plane_between_its_vertices_is_found() {
        let mut m = Model::new();
        let mut points = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                let z = if i == 1 && j == 1 { 1.0 } else { 0.0 };
                points.push(Point::new(f64::from(i) * 0.5, f64::from(j) * 0.5, z));
            }
        }
        let knots = || KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap();
        let grid = ControlGrid::new(points, 3, 3).unwrap();
        let surface = BSplineSurface::new(knots(), knots(), &grid, T).unwrap();
        let face = crate::make_natural_face(&mut m, surface.into())
            .unwrap()
            .shape;
        let plane = below(Vector::new(0.0, 0.0, 1.0), Point::ORIGIN);
        let out = certify_face_within(&m, &face, &plane, FIT, T)
            .unwrap()
            .unwrap_err();
        assert!(out.found, "{out:?}");
        assert_eq!(out.past, 0.25, "{out:?}");
        // Raised by the quarter, the plane holds it.
        let raised = below(Vector::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, 0.25));
        assert_eq!(
            certify_face_within(&m, &face, &raised, FIT, T).unwrap(),
            Ok(())
        );
    }

    /// The least of pieces: a point is inside where it is inside any. The
    /// unit cylinder's side lies inside the union of the half-spaces
    /// `x <= 1` and `x >= -1` everywhere, so proved; and inside the union of
    /// `x <= -0.5` and `x >= 0.5` nowhere near `x = 0`, so a point is found.
    #[test]
    fn the_least_of_pieces_proves_a_union_and_finds_what_falls_between() {
        let mut m = Model::new();
        let solid = crate::make_cylinder(&mut m, Frame::WORLD, 1.0, 3.0, T)
            .unwrap()
            .shape;
        let side = &faces_of(&m, &solid, ogeom_geom::SurfaceKind::Cylinder)[0];
        let x = Vector::new(1.0, 0.0, 0.0);
        let both = HalfSpaces::Any(vec![
            below(x, Point::new(1.0, 0.0, 0.0)),
            below(-x, Point::new(-1.0, 0.0, 0.0)),
        ]);
        assert_eq!(
            certify_face_within(&m, side, &both, FIT, T).unwrap(),
            Ok(())
        );
        let apart = HalfSpaces::Any(vec![
            below(x, Point::new(-0.5, 0.0, 0.0)),
            below(-x, Point::new(0.5, 0.0, 0.0)),
        ]);
        let out = certify_face_within(&m, side, &apart, FIT, T)
            .unwrap()
            .unwrap_err();
        assert!(
            out.found && out.past > FIT && out.past <= 0.5 + 1e-12,
            "{out:?}"
        );
    }
}
