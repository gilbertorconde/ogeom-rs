//! The exact half of the drawing pipeline.
//!
//! The polygonal path draws what the *mesh* says: silhouettes are interior
//! mesh edges whose triangles disagree about facing the eye, and visibility
//! is occlusion sampling against the triangles. Both are as good as the
//! chord, and no better.
//!
//! This is the other half. A silhouette is where the surface's own normal
//! turns perpendicular to the view, which for the elementary surfaces is a
//! curve with a closed form: a great circle on a sphere, a pair of rulings
//! on a cylinder or a cone. Visibility is decided by asking the *faces*
//! whether anything stands between a point and the eye: an exact
//! curve/surface interference and a trim test, not a triangle count. What is
//! still sampled is the drawing itself, because a drawing is polylines. The
//! curves it samples and the classification it carries are exact.
//!
//! A surface whose silhouette has no closed form (a torus, a spline) has it
//! marched: the silhouette is one equation on the surface's own chart, and
//! the shared walker follows it.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Axis, Circle, Direction, Frame, Point, Point2, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

use crate::crossings::{Contours, locate};
use crate::project::{Drawing, DrawnCurve, Source, View, Visibility};

/// One silhouette curve, and the face it belongs to.
#[derive(Debug, Clone)]
pub struct Silhouette {
    /// The face whose surface turns away here.
    pub face: Shape,
    /// The curve, in space.
    pub curve: Curve,
    /// The portion of it that lies within the face's trim.
    pub range: (f64, f64),
}

/// The exact silhouettes of a shape, seen along `direction`.
///
/// A silhouette is the locus where the surface normal is perpendicular to
/// the view: for a sphere the great circle whose plane the direction is
/// normal to, for a cylinder the two rulings furthest to either side, for a
/// cone the two rulings through its apex where the same holds. Planes have
/// none (a plane either faces the eye or does not), and a face whose
/// surface has no closed-form silhouette has it marched and fitted.
///
/// Each curve comes back trimmed to the stretch that lies within its own
/// face, decided by sampling the face's trim, so a silhouette on a face
/// that was cut away is absent rather than drawn through thin air.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// direction has no length. A marched silhouette whose walk or fit fails
/// contributes no curve.
pub fn silhouettes(
    model: &Model,
    shape: &Shape,
    direction: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Silhouette>> {
    let magnitude = direction.magnitude();
    if magnitude <= tol.confusion() {
        ogeom_bail!(Construction, "a silhouette needs a direction to look along");
    }
    let along = direction / magnitude;
    let deflection = Deflection::default();

    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let per_face = ogeom_core::parallel::map_ordered(&faces, |_, face| {
        face_silhouettes(model, face, along, deflection, tol)
    });
    let mut out = Vec::new();
    for found in per_face {
        out.extend(found?);
    }
    Ok(out)
}

/// The silhouettes of one face, trimmed to it.
fn face_silhouettes(
    model: &Model,
    face: &Shape,
    along: Vector,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Silhouette>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(Vec::new());
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(Vec::new());
    };
    let placement = face.transform(model.datums())?;
    let world = ogeom_geom::Transformable::transformed(surface, &placement, tol)?;
    let candidates = match &world {
        SurfaceGeometry::Plane(_) => Vec::new(),
        SurfaceGeometry::Sphere(s) => {
            // The great circle whose plane has the view as its normal:
            // every normal on it is radial, so every one is
            // perpendicular to the view.
            let sphere = s.sphere();
            let axis = Direction::new(along, tol)?;
            let frame = Frame::new(sphere.centre(), axis, perpendicular(along, tol)?, tol)?;
            vec![Curve::Circle(CircleCurve::new(Circle::new(
                frame,
                sphere.radius(),
                tol,
            )?))]
        }
        SurfaceGeometry::Cylinder(c) => {
            // The two rulings where the radial direction is
            // perpendicular to the view: the axis stepped sideways by
            // the radius, either way.
            let cylinder = c.cylinder();
            let axis = cylinder.frame().z().vector();
            let sideways = axis.cross(along);
            let m = sideways.magnitude();
            if m <= tol.angular() {
                // Looking down the axis: the whole rim is the outline,
                // and the face's own boundary already draws it.
                Vec::new()
            } else {
                let sideways = sideways / m;
                [1.0, -1.0]
                    .iter()
                    .map(|sign| {
                        let at = cylinder.frame().origin() + sideways * (cylinder.radius() * sign);
                        Curve::Line(LineCurve::new(Axis::new(at, cylinder.frame().z())))
                    })
                    .collect()
            }
        }
        SurfaceGeometry::Cone(c) => {
            // The same question on a cone: the rulings whose own normal
            // is perpendicular to the view. The normal of a ruling at
            // angle u is radial tilted by the half-angle, so the
            // condition is a linear one in (cos u, sin u) and has two
            // roots, or none, when the eye is inside the cone's own
            // angle and nothing turns away.
            let cone = c.cone();
            let frame = cone.frame();
            let (x, y, z) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
            let (sin, cos) = cone.half_angle().sin_cos();
            // n(u) = cos(half) * (x cos u + y sin u) - sin(half) * z
            let (a, b) = (cos * along.dot(x), cos * along.dot(y));
            let c0 = -sin * along.dot(z);
            let r = a.hypot(b);
            if r <= tol.angular() || c0.abs() > r {
                Vec::new()
            } else {
                let phase = b.atan2(a);
                let spread = (-c0 / r).acos();
                [phase + spread, phase - spread]
                    .iter()
                    .map(|u| {
                        let radial = x * u.cos() + y * u.sin();
                        let apex = cone.apex();
                        let direction =
                            radial * cone.half_angle().sin() + z * cone.half_angle().cos();
                        Direction::new(direction, tol)
                            .map(|d| Curve::Line(LineCurve::new(Axis::new(apex, d))))
                    })
                    .collect::<OgeomResult<Vec<Curve>>>()?
            }
        }
        // No closed form: a torus, a spline. The silhouette is still
        // one equation on the surface's own chart, and one equation in
        // two unknowns is a curve, so it is *walked* rather than refused.
        other => marched_silhouettes(other, along, tol)?,
    };

    if candidates.is_empty() {
        return Ok(Vec::new());
    }
    let rings = ogeom_mesh::face_boundary(model, face, deflection, tol)?;
    // Every trim station projects onto the same surface, so its seeding
    // grid is evaluated once for the face.
    let seeds = ogeom_algo::SurfaceSeeds::over(&world, 24, tol)?;
    let mut out = Vec::new();
    for curve in candidates {
        for range in within_trim(model, face, &world, &seeds, &rings, &curve, tol)? {
            out.push(Silhouette {
                face: face.clone(),
                curve: curve.clone(),
                range,
            });
        }
    }
    Ok(out)
}

/// The reflect lines of a shape under a light: where the surface turns away
/// from the *light* rather than from the eye.
///
/// The same locus as a silhouette, asked of a different direction, which is
/// what a reflect line is, and why the two share a construction. A surface
/// inspected this way shows its own creases: the lines move a long way for a
/// small change in curvature.
///
/// # Errors
///
/// As [`silhouettes`].
pub fn reflect_lines(
    model: &Model,
    shape: &Shape,
    light: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Silhouette>> {
    silhouettes(model, shape, light, tol)
}

/// The isoparametric curves of a face: `u_count` at constant `u`, `v_count`
/// at constant `v`, each trimmed to the stretches that lie on the face.
///
/// Evenly spaced across the face's own parameter window, excluding its
/// edges, because an isoparametric at the window's edge is the face's
/// boundary and the boundary is already drawn.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// face carries no surface.
pub fn iso_curves(
    model: &Model,
    face: &Shape,
    u_count: usize,
    v_count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Point>>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;
    let world = ogeom_geom::Transformable::transformed(&surface, &placement, tol)?;
    let ((u0, u1), (v0, v1)) = world.domain();
    let rings = ogeom_mesh::face_boundary(model, face, Deflection::default(), tol)?;

    const ALONG: usize = 64;
    let mut out = Vec::new();
    for (count, constant_u) in [(u_count, true), (v_count, false)] {
        for i in 1..=count {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a curve index, far below the mantissa"
            )]
            let f = i as f64 / (count + 1) as f64;
            let mut run: Vec<Point> = Vec::new();
            for k in 0..=ALONG {
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "a station index, far below the mantissa"
                )]
                let g = k as f64 / ALONG as f64;
                let (u, v) = if constant_u {
                    ((u1 - u0).mul_add(f, u0), (v1 - v0).mul_add(g, v0))
                } else {
                    ((u1 - u0).mul_add(g, u0), (v1 - v0).mul_add(f, v0))
                };
                if inside_rings(&rings, Point2::new(u, v)) {
                    run.push(world.point_at(u, v, tol)?);
                } else if run.len() >= 2 {
                    out.push(std::mem::take(&mut run));
                } else {
                    run.clear();
                }
            }
            if run.len() >= 2 {
                out.push(run);
            }
        }
    }
    Ok(out)
}

/// Project a shape into a drawing whose silhouettes and visibility are
/// exact.
///
/// The edges and silhouettes are sampled at `deflection` to become
/// polylines, because a drawing is polylines. What is not sampled is the
/// *geometry* they sample (exact curves on the surfaces, not mesh edges),
/// or the classification, which asks the faces themselves whether anything
/// stands between a point and the eye.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape has no faces, or the direction has no length. A marched silhouette
/// whose walk or fit fails contributes no curve.
pub fn project_exact(
    model: &Model,
    shape: &Shape,
    view: &View,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Drawing> {
    let faces = Blockers::new(model, shape, view, tol)?;
    if faces.list.is_empty() {
        ogeom_bail!(Construction, "a shape with no faces draws nothing");
    }

    // Every curve is sampled before any is classified: the projection of
    // each is a contour the others may pass behind.
    let edges = explore_unique(model, shape, ShapeType::Edge)?;
    let mut traced: Vec<Traced> = ogeom_core::parallel::map_ordered(&edges, |_, edge| {
        traced_edge(model, edge, deflection, tol)
    })
    .into_iter()
    .flatten()
    .collect();
    for silhouette in silhouettes(model, shape, view.toward_eye(), tol)? {
        let line = ogeom_mesh::discretize(&silhouette.curve, silhouette.range, deflection, tol)?;
        traced.push(Traced {
            curve: silhouette.curve,
            placement: Transform::IDENTITY,
            parameters: line.parameters,
            points: line.points,
            source: Source::Silhouette,
        });
    }
    let projected: Vec<Vec<Point2>> = traced
        .iter()
        .map(|t| t.points.iter().map(|p| view.project(*p)).collect())
        .collect();
    let contours = Contours::new(
        projected
            .iter()
            .flat_map(|line| line.windows(2).map(|w| (w[0], w[1])))
            .collect(),
    );

    // Each curve is classified on its own, and the runs are gathered in
    // curve order.
    let classified = ogeom_core::parallel::map_ordered(&traced, |i, curve| {
        let mut drawing = Drawing::default();
        classify(
            &mut drawing,
            curve,
            &projected[i],
            view,
            &faces,
            &contours,
            tol,
        )
        .map(|()| drawing)
    });
    let mut drawing = Drawing::default();
    for one in classified {
        let one = one?;
        drawing.visible.extend(one.visible);
        drawing.hidden.extend(one.hidden);
    }
    Ok(drawing)
}

/// A curve to draw: its exact geometry, and the samples its polyline is
/// drawn through.
struct Traced {
    curve: Curve,
    /// Where the curve's own coordinates stand in the world.
    placement: Transform,
    /// The curve parameter at each sample, in drawing order.
    parameters: Vec<f64>,
    /// The samples, in the world.
    points: Vec<Point>,
    source: Source,
}

impl Traced {
    /// The exact point at a position along the samples (`k + f`: the
    /// fraction `f` of the way from sample `k` to the next in parameter).
    fn at(&self, position: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (k, f) = locate(position, self.parameters.len());
        let (t0, t1) = (self.parameters[k], self.parameters[k + 1]);
        Ok(self
            .placement
            .apply(self.curve.point_at((t1 - t0).mul_add(f, t0), tol)?))
    }
}

/// An edge's curve and samples, in the edge's own direction; `None` for an
/// edge with no curve to draw.
fn traced_edge(
    model: &Model,
    edge: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> Option<Traced> {
    let data = model.node(edge)?.data().as_edge()?;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return None;
    };
    let curve = model.geometry().curve(*curve)?.clone();
    let placement = edge.transform(model.datums()).ok()?;
    let line = ogeom_mesh::discretize(&curve, *range, deflection, tol).ok()?;
    let mut points: Vec<Point> = line.points.iter().map(|p| placement.apply(*p)).collect();
    let mut parameters = line.parameters;
    if edge.orientation() == Orientation::Reversed {
        points.reverse();
        parameters.reverse();
    }
    Some(Traced {
        curve,
        placement,
        parameters,
        points,
        source: Source::Edge(edge.clone()),
    })
}

/// A face that can stand between a point and the eye: its surface in world
/// space and its trim as chart rings.
struct Blocker {
    surface: ogeom_intersect::PreparedSurface,
    rings: Vec<Vec<Point2>>,
    /// The projection of the corners of the face's box: a face whose
    /// projected box misses a point's projection cannot hide it, and is not
    /// intersected.
    low: Point2,
    high: Point2,
    /// The depth of the box corner nearest the eye: a face wholly behind a
    /// point cannot hide it either.
    front: f64,
}

/// The faces of a shape seen from one view, binned by their projected
/// boxes so a point is tested only against the faces whose boxes cover it.
struct Blockers {
    list: Vec<Blocker>,
    /// For each cell of a square grid over the drawing, the faces whose
    /// projected box meets it, ascending.
    cells: Vec<Vec<u32>>,
    low: Point2,
    size: f64,
    side: usize,
}

impl Blockers {
    fn new(model: &Model, shape: &Shape, view: &View, tol: Tolerances) -> OgeomResult<Self> {
        let faces = explore_unique(model, shape, ShapeType::Face)?;
        let built =
            ogeom_core::parallel::map_ordered(&faces, |_, face| blocker(model, face, view, tol));
        let mut list = Vec::with_capacity(built.len());
        for one in built {
            list.extend(one?);
        }

        let (mut low, mut high) = (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for b in &list {
            low = Point2::new(low.x.min(b.low.x), low.y.min(b.low.y));
            high = Point2::new(high.x.max(b.high.x), high.y.max(b.high.y));
        }
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let side = ((list.len() as f64).sqrt().ceil() as usize * 2).clamp(1, 256);
        let span = (high.x - low.x).max(high.y - low.y);
        #[allow(clippy::cast_precision_loss)]
        let size = if span.is_finite() && span > 0.0 {
            span / side as f64
        } else {
            1.0
        };
        if !low.x.is_finite() || !low.y.is_finite() {
            low = Point2::new(0.0, 0.0);
        }
        let mut out = Self {
            list: Vec::new(),
            cells: vec![Vec::new(); side * side],
            low,
            size,
            side,
        };
        for (index, b) in list.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation)]
            let index = index as u32;
            let (c0, c1) = (out.cell(b.low.x, low.x), out.cell(b.high.x, low.x));
            let (r0, r1) = (out.cell(b.low.y, low.y), out.cell(b.high.y, low.y));
            for r in r0..=r1 {
                for c in c0..=c1 {
                    out.cells[r * side + c].push(index);
                }
            }
        }
        out.list = list;
        Ok(out)
    }

    /// The grid column (or row) of a coordinate, clamped to the grid; a
    /// box's ends land in the cells bracketing every point it covers.
    fn cell(&self, x: f64, lo: f64) -> usize {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = ((x - lo) / self.size).floor().max(0.0) as usize;
        k.min(self.side - 1)
    }

    /// The faces whose projected box covers `seen` and that stand at least
    /// partly in front of `depth`, by index, in face order.
    fn covering(&self, seen: Point2, depth: f64) -> impl Iterator<Item = usize> + '_ {
        let cell = self.cell(seen.y, self.low.y) * self.side + self.cell(seen.x, self.low.x);
        self.cells[cell]
            .iter()
            .map(|&i| i as usize)
            .filter(move |&i| self.list[i].covers(seen, depth))
    }
}

impl Blocker {
    /// Whether the face's projected box covers `seen` with some of the box
    /// in front of `depth`.
    fn covers(&self, seen: Point2, depth: f64) -> bool {
        self.low.x <= seen.x
            && seen.x <= self.high.x
            && self.low.y <= seen.y
            && seen.y <= self.high.y
            && self.front > depth
    }
}

/// One face as a blocker; `None` for a face with no surface or an empty box.
fn blocker(
    model: &Model,
    face: &Shape,
    view: &View,
    tol: Tolerances,
) -> OgeomResult<Option<Blocker>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placement = face.transform(model.datums())?;
    let surface = ogeom_geom::Transformable::transformed(surface, &placement, tol)?;
    let rings = ogeom_mesh::face_boundary(model, face, Deflection::default(), tol)?;
    let bound = ogeom_algo::shape_bounds(model, face, tol)?.expanded(tol.confusion() * 1e3);
    let corners = bound.corners();
    if corners.is_empty() {
        return Ok(None);
    }
    let (mut low, mut high) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    let mut front = f64::NEG_INFINITY;
    for c in &corners {
        let q = view.project(*c);
        low = Point2::new(low.x.min(q.x), low.y.min(q.y));
        high = Point2::new(high.x.max(q.x), high.y.max(q.y));
        front = front.max(view.depth(*c));
    }
    Ok(Some(Blocker {
        surface: ogeom_intersect::PreparedSurface::new(
            surface,
            ogeom_intersect::CurveSurfaceOptions::default(),
            tol,
        ),
        rings,
        low,
        high,
        front,
    }))
}

/// Split a curve into visible and hidden runs, asking the faces.
///
/// The samples are cut where the projection crosses another drawn curve,
/// and each piece is classified at its middle on the exact curve. Where two
/// neighbouring pieces disagree, the change is bisected on the curve until
/// its two sides are the confusion tolerance apart, so the runs meet where
/// the curve passes behind an outline rather than at a sample.
fn classify(
    drawing: &mut Drawing,
    traced: &Traced,
    projected: &[Point2],
    view: &View,
    faces: &Blockers,
    contours: &Contours,
    tol: Tolerances,
) -> OgeomResult<()> {
    if traced.points.len() < 2 {
        return Ok(());
    }
    // Neighbouring points along a curve are mostly hidden by the same face,
    // so the face that hid the last one is asked first.
    let last = std::cell::Cell::new(None);
    let seen = |position: f64| -> OgeomResult<Visibility> {
        Ok(
            if occluded(traced.at(position, tol)?, view, faces, &last, tol)? {
                Visibility::Hidden
            } else {
                Visibility::Visible
            },
        )
    };
    let change_between = |mut lo: f64, mut hi: f64, was: Visibility| -> OgeomResult<f64> {
        for _ in 0..64 {
            if traced.at(lo, tol)?.distance(traced.at(hi, tol)?) <= tol.confusion() {
                break;
            }
            let mid = f64::midpoint(lo, hi);
            if seen(mid)? == was {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(f64::midpoint(lo, hi))
    };
    let drawn_at =
        |position: f64| -> OgeomResult<Point2> { Ok(view.project(traced.at(position, tol)?)) };
    let mut flush = |run: &mut Vec<Point2>, visibility: Visibility| {
        if run.len() < 2 {
            run.clear();
            return;
        }
        let curve = DrawnCurve {
            points: std::mem::take(run),
            visibility,
            source: traced.source.clone(),
        };
        if visibility == Visibility::Hidden {
            drawing.hidden.push(curve);
        } else {
            drawing.visible.push(curve);
        }
    };

    let positions = contours.split(projected, tol);
    let middles: Vec<f64> = positions
        .windows(2)
        .map(|w| f64::midpoint(w[0], w[1]))
        .collect();
    let mut verdicts = Vec::with_capacity(middles.len());
    for middle in &middles {
        verdicts.push(seen(*middle)?);
    }
    let Some(&first) = verdicts.first() else {
        return Ok(());
    };
    let mut run = vec![drawn_at(positions[0])?];
    let mut held = first;
    for (i, &verdict) in verdicts.iter().enumerate() {
        if verdict != held {
            // The run so far ends at the change and the next starts there;
            // the piece boundary goes to whichever side it is on.
            let boundary = positions[i];
            let change = change_between(middles[i - 1], middles[i], held)?;
            let at_change = drawn_at(change)?;
            if change < boundary {
                run.pop();
                run.push(at_change);
                flush(&mut run, held);
                run = vec![at_change, drawn_at(boundary)?];
            } else {
                run.push(at_change);
                flush(&mut run, held);
                run = vec![at_change];
            }
            held = verdict;
        }
        run.push(drawn_at(positions[i + 1])?);
    }
    flush(&mut run, held);
    Ok(())
}

/// Whether anything stands between `at` and the eye.
///
/// The answer is whether any face hides the point, so the order the faces
/// are asked in changes only the time: the face in `last` is asked first,
/// and the face that hides the point is left there.
fn occluded(
    at: Point,
    view: &View,
    faces: &Blockers,
    last: &std::cell::Cell<Option<usize>>,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let toward = view.toward_eye();
    let magnitude = toward.magnitude();
    if magnitude <= tol.confusion() {
        return Ok(false);
    }
    let direction = toward / magnitude;
    // Far enough to leave any shape it started inside, and started far
    // enough along not to strike the surface the point is on.
    let reach = 1e6;
    let clearance = tol.confusion() * 1e3;
    let ray = Curve::Line(LineCurve::new(Axis::new(
        at,
        Direction::new(direction, tol)?,
    )));
    let hides = |face: &Blocker| -> OgeomResult<bool> {
        let found = face.surface.intersect(&ray)?;
        Ok(found.crossings.iter().any(|piercing| {
            piercing.on_curve > clearance
                && piercing.on_curve < reach
                && inside_rings(
                    &face.rings,
                    Point2::new(piercing.on_surface.0, piercing.on_surface.1),
                )
        }))
    };
    let (seen, depth) = (view.project(at), view.depth(at));
    let first = last.get().filter(|&i| faces.list[i].covers(seen, depth));
    if let Some(i) = first
        && hides(&faces.list[i])?
    {
        return Ok(true);
    }
    for i in faces.covering(seen, depth) {
        if Some(i) != first && hides(&faces.list[i])? {
            last.set(Some(i));
            return Ok(true);
        }
    }
    Ok(false)
}

/// The stretches of a curve that lie within a face's trim.
fn within_trim(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    seeds: &ogeom_algo::SurfaceSeeds,
    rings: &[Vec<Point2>],
    curve: &Curve,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64)>> {
    let (t0, t1) = curve.domain();
    // A line's domain is the whole real line as far as the type is
    // concerned. A silhouette on one is only interesting where the face is,
    // so an unbounded curve is walked over the face's own reach.
    let (t0, t1) = if t0.is_finite() && t1.is_finite() && t1 - t0 < 1e6 {
        (t0, t1)
    } else {
        let mut bound = ogeom_math::Aabb::EMPTY;
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                bound = bound.with_point(vertex.transform(model.datums())?.apply(data.point));
            }
        }
        let reach = bound.diagonal().max(1.0);
        let centre = bound.centre().unwrap_or(Point::ORIGIN);
        let at = ogeom_algo::project_on_curve(curve, centre, 64, tol)?.parameter;
        (at - reach, at + reach)
    };

    const STATIONS: usize = 96;
    let held_at = |t: f64| -> OgeomResult<bool> {
        let Ok(point) = curve.point_at(t, tol) else {
            return Ok(false);
        };
        let projection = seeds.project(surface, point, tol)?;
        let (u, v) = projection.parameters;
        // The projection clamps to the surface's own window, so a point
        // just past the end of a face comes back with a foot *at* the end
        // and the overshoot as its distance. Holding that to the confusion
        // tolerance is what stops a silhouette running off its own face.
        Ok(projection.distance <= tol.confusion() && inside_rings(rings, Point2::new(u, v)))
    };
    // Where the answer changes between two stations, the edge of the face
    // is between them. Bisecting says where to a part in a million of a
    // station, so a silhouette ends *on* its face rather than a station
    // past it.
    let edge_between = |inside: f64, outside: f64| -> OgeomResult<f64> {
        let (mut lo, mut hi) = (inside, outside);
        for _ in 0..40 {
            let mid = f64::midpoint(lo, hi);
            if held_at(mid)? {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    };

    let mut out = Vec::new();
    let mut open: Option<f64> = None;
    let mut previous: Option<(f64, bool)> = None;
    for k in 0..=STATIONS {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = (t1 - t0).mul_add(k as f64 / STATIONS as f64, t0);
        let held = held_at(t)?;
        match (held, open) {
            (true, None) => {
                open = Some(match previous {
                    Some((was, false)) => edge_between(t, was)?,
                    _ => t,
                });
            }
            (false, Some(from)) => {
                let to = match previous {
                    Some((was, true)) => edge_between(was, t)?,
                    _ => t,
                };
                if to - from > tol.parametric() {
                    out.push((from, to));
                }
                open = None;
            }
            _ => {}
        }
        previous = Some((t, held));
    }
    if let Some(from) = open
        && t1 - from > tol.parametric()
    {
        out.push((from, t1));
    }
    Ok(out)
}

/// Even-odd containment against chart rings.
fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = (b.x - a.x).mul_add((p.y - a.y) / (b.y - a.y), a.x);
                if x > p.x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// Any unit vector perpendicular to `v`.
fn perpendicular(v: Vector, tol: Tolerances) -> OgeomResult<Direction> {
    let seed = if v.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    Direction::new(v.cross(seed), tol)
}

// --- the marched silhouette --------------------------------------------------

/// A surface's silhouette, as a condition the shared walker can follow.
///
/// The whole content of a silhouette is one equation on the surface's own
/// chart: the normal is square to the view,
///
/// > `n(u, v) · d = 0`
///
/// which is one equation in two unknowns, and one equation in two unknowns
/// is a curve. So a torus's silhouette needs no machinery a surface
/// intersection did not already need: it is the same walk, following a
/// different condition.
///
/// Stated with the **unit** normal. The unnormalized `Sᵤ × Sᵥ` has the same
/// zero set and a simpler derivative, but its residual carries the surface's
/// own scale: a correction that drives `|Sᵤ × Sᵥ| · d` below a *length*
/// tolerance demands the angle tighter by the surface's own size, and the
/// walk halves its step until it crawls. A dimensionless residual is a
/// dimensionless tolerance.
struct SilhouetteOn<'s> {
    surface: &'s SurfaceGeometry,
    along: Vector,
    /// A length scale, for the walker's step control.
    reach: f64,
}

impl ogeom_intersect::walk::Condition for SilhouetteOn<'_> {
    fn unknowns(&self) -> usize {
        2
    }

    fn position(&self, x: &[f64], tol: Tolerances) -> Option<Point> {
        self.surface.point_at(x[0], x[1], tol).ok()
    }

    fn position_gradient(&self, x: &[f64], tol: Tolerances) -> Option<Vec<Vector>> {
        let (du, dv) = self.surface.d1_at(x[0], x[1], tol).ok()?;
        Some(vec![du, dv])
    }

    fn system(&self, x: &[f64], tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
        let (su, sv) = self.surface.d1_at(x[0], x[1], tol).ok()?;
        let (suu, suv, svv) = self.surface.d2_at(x[0], x[1], tol).ok()?;
        let cross = su.cross(sv);
        let length = cross.magnitude();
        if length <= tol.confusion() {
            return None;
        }
        let normal = cross / length;
        // The unit normal's own derivative: the part of the unnormalized
        // one's across the normal, over the length. The projection is what
        // keeps a unit vector unit.
        let across = |d: Vector| (d - normal * d.dot(normal)) / length;
        let du = across(suu.cross(sv) + su.cross(suv));
        let dv = across(suv.cross(sv) + su.cross(svv));
        Some((
            vec![normal.dot(self.along)],
            vec![vec![du.dot(self.along), dv.dot(self.along)]],
        ))
    }

    fn clamp(&self, x: &mut [f64]) {
        let ((ua, ub), (va, vb)) = self.surface.domain();
        let hold = |value: f64, lo: f64, hi: f64, periodic: bool| {
            if periodic && hi > lo {
                lo + (value - lo).rem_euclid(hi - lo)
            } else {
                value.clamp(lo, hi)
            }
        };
        x[0] = hold(x[0], ua, ub, self.surface.is_periodic_u());
        x[1] = hold(x[1], va, vb, self.surface.is_periodic_v());
    }

    fn outside(&self, x: &[f64], tol: Tolerances) -> bool {
        let ((ua, ub), (va, vb)) = self.surface.domain();
        let band = tol.parametric();
        (!self.surface.is_periodic_u() && (x[0] < ua - band || x[0] > ub + band))
            || (!self.surface.is_periodic_v() && (x[1] < va - band || x[1] > vb + band))
    }

    fn near_edge(&self, x: &[f64]) -> bool {
        let ((ua, ub), (va, vb)) = self.surface.domain();
        let near = |value: f64, lo: f64, hi: f64| {
            let reach = (hi - lo) * 1e-6;
            value <= lo + reach || value >= hi - reach
        };
        (!self.surface.is_periodic_u() && near(x[0], ua, ub))
            || (!self.surface.is_periodic_v() && near(x[1], va, vb))
    }

    fn extent(&self) -> f64 {
        self.reach
    }
}

/// How far a point stands from a polyline, segment by segment.
fn on_polyline(line: &[Point], p: Point) -> f64 {
    let mut best = f64::INFINITY;
    for pair in line.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let d = b - a;
        let len2 = d.dot(d);
        let t = if len2 > 0.0 {
            ((p - a).dot(d) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        best = best.min(p.distance(a + d * t));
    }
    best
}

/// Silhouettes of a surface with no closed form, marched.
///
/// Seeded the way the intersector seeds: the chart is sampled on a grid and
/// every sign change of `n · d` between neighbours is a starting point,
/// refined onto the condition before the walk begins. A silhouette loop
/// smaller than one grid cell is stepped over.
///
/// The walked polylines are fitted to curves at `chord` through the walked
/// points, which lie on the silhouette, so what comes back is worth the
/// chord rather than exactness; nothing is stated with it.
fn marched_silhouettes(
    surface: &SurfaceGeometry,
    along: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Curve>> {
    use ogeom_intersect::walk::Condition as _;
    let options = ogeom_intersect::Marching {
        chord: tol.confusion() * 1e2,
        ..ogeom_intersect::Marching::default()
    };
    let ((ua, ub), (va, vb)) = surface.domain();
    // The step control wants a *length*, and it must be the surface's own:
    // a torus's face is bounded by a seam and one vertex, so its vertices'
    // bounding box is a point, and a step control fed that walks the whole
    // ring in steps of a ten-thousandth.
    let reach = {
        let mut bound = ogeom_math::Aabb::EMPTY;
        for i in 0..=8 {
            for j in 0..=8 {
                let u = (ub - ua).mul_add(f64::from(i) / 8.0, ua);
                let v = (vb - va).mul_add(f64::from(j) / 8.0, va);
                if let Ok(p) = surface.point_at(u, v, tol) {
                    bound = bound.with_point(p);
                }
            }
        }
        bound.diagonal().max(tol.confusion() * 1e3)
    };
    let condition = SilhouetteOn {
        surface,
        along,
        reach,
    };
    let value = |u: f64, v: f64| -> Option<f64> {
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        Some(su.cross(sv).dot(along))
    };

    // Seeds: a sign change between grid neighbours, bisected to the crossing
    // and then corrected onto the condition by the walker's own solve.
    let mut seeds: Vec<[f64; 2]> = Vec::new();
    let steps = options.grid;
    #[expect(clippy::cast_precision_loss, reason = "a grid index")]
    let at = |i: usize, n: usize, lo: f64, hi: f64| lo + (hi - lo) * (i as f64) / (n as f64);
    // Each grid point is evaluated once and read by up to three cell edges.
    let grid: Vec<Option<f64>> = (0..=steps)
        .flat_map(|i| (0..=steps).map(move |j| (i, j)))
        .map(|(i, j)| value(at(i, steps, ua, ub), at(j, steps, va, vb)))
        .collect();
    for i in 0..=steps {
        for j in 0..=steps {
            let (u, v) = (at(i, steps, ua, ub), at(j, steps, va, vb));
            let Some(here) = grid[i * (steps + 1) + j] else {
                continue;
            };
            for (du, dv) in [(1_usize, 0_usize), (0, 1)] {
                if i + du > steps || j + dv > steps {
                    continue;
                }
                let (u2, v2) = (at(i + du, steps, ua, ub), at(j + dv, steps, va, vb));
                let Some(there) = grid[(i + du) * (steps + 1) + j + dv] else {
                    continue;
                };
                if here.signum() == there.signum() || here == 0.0 {
                    continue;
                }
                // Bisect to the crossing along the cell edge.
                let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
                for _ in 0..40 {
                    let mid = f64::midpoint(lo, hi);
                    let (um, vm) = (u + (u2 - u) * mid, v + (v2 - v) * mid);
                    let Some(m) = value(um, vm) else { break };
                    if m.signum() == here.signum() {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let mid = f64::midpoint(lo, hi);
                seeds.push([u + (u2 - u) * mid, v + (v2 - v) * mid]);
            }
        }
    }

    let mut out = Vec::new();
    let mut walked_points: Vec<Vec<Point>> = Vec::new();
    for seed in seeds {
        let mut start = seed;
        condition.clamp(&mut start);
        // A seed already covered by a curve found earlier is the same branch
        // met in another cell, not a new one.
        let Some(here) = condition.position(&start, tol) else {
            continue;
        };
        // Against the polyline, not its vertices: the walk's own step is
        // hundreds of times the chord, so a seed landing between two points
        // of a curve already found would otherwise be walked all over again.
        // A torus seen down its axis has a seed in every grid column of both
        // equators.
        if walked_points
            .iter()
            .any(|line| on_polyline(line, here) <= options.chord * 32.0)
        {
            continue;
        }
        let Ok(walked) = ogeom_intersect::walk::follow(&condition, &start, options, tol) else {
            continue;
        };
        if walked.points.len() < 4 {
            continue;
        }
        // Fitted through the walked points at the walk's chord; the fit is
        // measured at those points only.
        let Ok(fitted) = ogeom_geom::fit::fit_points(&walked.points, 3, options.chord, tol) else {
            continue;
        };
        walked_points.push(walked.points);
        out.push(Curve::BSpline(fitted.curve));
    }
    Ok(out)
}
