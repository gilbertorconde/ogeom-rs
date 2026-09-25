//! Dividing a shape into smaller pieces: faces cut along lines of their own
//! parameters, edges cut at parameters of their own.
//!
//! Every piece keeps the curve or surface it was cut from, and so its
//! parameters: nothing is refitted, and a pcurve that held on the whole
//! holds on each piece. A cut face's new boundary is the surface's own
//! iso-curve, exact for every surface whose iso-curves have a closed form
//! (planes, drums, cones, balls, tori, extrusions, revolutions and
//! splines).
//!
//! The divisions on top pick where to cut: at the knots where a spline is
//! less smooth than asked ([`divide_by_continuity`]), wherever a face turns
//! further than an angle ([`divide_by_angle`]), until no face is larger than
//! an area ([`divide_by_area`]), and at every knot, each piece's geometry
//! then restated as the single Bézier span it covers ([`to_bezier`]).

use std::collections::{HashMap, HashSet};

use ogeom_algo::{Built, History, edge_vertices, surface_properties};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSpline2d, CircleCurve, Continuity, Curve, Curve2d as _, Curve3d as _, LineCurve, PlanarCurve,
    Surface as _, SurfaceGeometry, Transformable as _,
};
use ogeom_math::{Axis, Circle, Direction, Frame, KnotVector, Point, Point2, Transform};
use ogeom_mesh::Deflection;
use ogeom_topo::{
    EdgeData, EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId,
    TShapeId, VertexData, explore_unique,
};

use crate::Reshape;

/// A line of a face's parameters: `u = at` or `v = at`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsoLine {
    /// The line `u = at`, running along `v`.
    U(f64),
    /// The line `v = at`, running along `u`.
    V(f64),
}

impl IsoLine {
    const fn at(self) -> f64 {
        match self {
            Self::U(c) | Self::V(c) => c,
        }
    }

    /// The parameter the line holds fixed, read off a chart point.
    const fn fixed(self, p: Point2) -> f64 {
        match self {
            Self::U(_) => p.x,
            Self::V(_) => p.y,
        }
    }

    /// The parameter the line runs along.
    const fn free(self, p: Point2) -> f64 {
        match self {
            Self::U(_) => p.y,
            Self::V(_) => p.x,
        }
    }

    const fn point(self, free: f64) -> Point2 {
        match self {
            Self::U(c) => Point2::new(c, free),
            Self::V(c) => Point2::new(free, c),
        }
    }
}

/// Cut `face` of `shape` in two or more along `line`, and rebuild `shape`
/// around the pieces. The edges the line crosses are cut where it crosses
/// them, in every face that holds them.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the line does not cross the face, runs along one of its edges, or the
/// surface has no closed-form iso-curve.
pub fn divide_face(
    model: &mut Model,
    shape: &Shape,
    face: &Shape,
    line: IsoLine,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let Some(reshape) = cut_face(model, face, line, tol)? else {
        ogeom_bail!(Construction, "the line {line:?} does not cross the face");
    };
    reshape.apply(model, shape)
}

/// Cut every edge and face of `shape` where its spline is less smooth than
/// `at_least`, so each piece is at least that smooth throughout. A knot is
/// judged by its multiplicity, so G1 and G2 ask what C1 and C2 ask.
///
/// Edges are cut at their curves' knots; faces along their surfaces' knot
/// lines, for a spline, an extrusion or revolution of a spline, and a
/// trimmed or offset surface over one.
///
/// # Errors
///
/// As [`divide_face`], for a cut that cannot be made.
pub fn divide_by_continuity(
    model: &mut Model,
    shape: &Shape,
    at_least: Continuity,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let order = match at_least {
        Continuity::C0 => return Ok(Built::new(shape.clone(), History::identity())),
        Continuity::G1 | Continuity::C1 => 1,
        Continuity::G2 | Continuity::C2 => 2,
        Continuity::CInfinity => usize::MAX,
    };
    let start = placed_baked(model, shape, tol)?;
    let edges = divide_edges(model, &start.shape, order, tol)?;
    let faces = divide_faces(model, &edges.shape, tol, |_, _, surface, _| {
        Ok(knot_lines(surface, order))
    })?;
    Ok(Built::new(
        faces.shape,
        start.history.then(&edges.history).then(&faces.history),
    ))
}

/// Cut every face of `shape` whose angular parameters (a drum's, cone's,
/// ball's, torus's or revolution's turn) sweep more than `max_angle`, into
/// equal pieces sweeping no more than it. A closed face comes apart at its
/// seam into open ones.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `max_angle` is not a positive angle, and as [`divide_face`].
pub fn divide_by_angle(
    model: &mut Model,
    shape: &Shape,
    max_angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !max_angle.is_finite() || max_angle <= tol.angular() {
        ogeom_bail!(Construction, "{max_angle} is not an angle to divide by");
    }
    let start = placed_baked(model, shape, tol)?;
    let faces = divide_faces(model, &start.shape, tol, |_, _, surface, bounds| {
        let (u_turns, v_turns) = angular(surface);
        let mut lines = Vec::new();
        for (turns, (lo, hi), iso) in [
            (u_turns, bounds.0, IsoLine::U as fn(f64) -> IsoLine),
            (v_turns, bounds.1, IsoLine::V),
        ] {
            let span = hi - lo;
            if turns && span > max_angle * (1.0 + 1e-9) {
                let pieces = (span / max_angle - 1e-6).ceil();
                lines.push(iso(lo + span / pieces));
            }
        }
        Ok(lines)
    })?;
    Ok(Built::new(faces.shape, start.history.then(&faces.history)))
}

/// Cut every face of `shape` larger than `max_area` in half, across its
/// longer extent, until none is.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `max_area` is not a positive area, and as [`divide_face`].
pub fn divide_by_area(
    model: &mut Model,
    shape: &Shape,
    max_area: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !max_area.is_finite() || max_area <= tol.confusion() * tol.confusion() {
        ogeom_bail!(Construction, "{max_area} is not an area to divide by");
    }
    let start = placed_baked(model, shape, tol)?;
    let faces = divide_faces(model, &start.shape, tol, |model, face, surface, bounds| {
        let area = surface_properties(model, face, Deflection::default(), tol)?.mass;
        if area <= max_area {
            return Ok(Vec::new());
        }
        let ((u0, u1), (v0, v1)) = bounds;
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        // How long the face runs each way, along the lines through its
        // middle.
        let length = |from: (f64, f64), to: (f64, f64)| -> OgeomResult<f64> {
            let mut total = 0.0;
            let mut last = surface.point_at(from.0, from.1, tol)?;
            for i in 1..=16 {
                let s = f64::from(i) / 16.0;
                let p = surface.point_at(
                    from.0 + (to.0 - from.0) * s,
                    from.1 + (to.1 - from.1) * s,
                    tol,
                )?;
                total += p.distance(last);
                last = p;
            }
            Ok(total)
        };
        let across_u = length((u0, vm), (u1, vm))?;
        let across_v = length((um, v0), (um, v1))?;
        Ok(if across_u >= across_v {
            vec![IsoLine::U(um), IsoLine::V(vm)]
        } else {
            vec![IsoLine::V(vm), IsoLine::U(um)]
        })
    })?;
    Ok(Built::new(faces.shape, start.history.then(&faces.history)))
}

/// `shape` with every curve and surface a single Bézier span: converted to
/// splines, cut at every knot, and each piece's geometry restated as the
/// span it covers, parameters kept.
///
/// # Errors
///
/// As [`ogeom_algo::to_nurbs`] and [`divide_by_continuity`].
pub fn to_bezier(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let nurbs = ogeom_algo::to_nurbs(model, shape, tol)?;
    let divided = divide_by_continuity(model, &nurbs.shape, Continuity::CInfinity, tol)?;

    // Each edge's curve restated as its span.
    let mut reshape = Reshape::new();
    for edge in explore_unique(model, &divided.shape, ShapeType::Edge)? {
        let data = edge_data(model, &edge)?;
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
            continue;
        };
        let Some(Curve::BSpline(spline)) = model.geometry().curve(curve).cloned() else {
            continue;
        };
        if spline.knots().distinct().len() <= 2 {
            continue;
        }
        let span = model
            .geometry_mut()
            .add_curve(Curve::BSpline(spline.segment(range, tol)?));
        let mut fresh = data.clone();
        for repr in &mut fresh.representations {
            if let EdgeRepr::Curve3d { curve, .. } = repr {
                *curve = span;
            }
        }
        let bounds = model.children_of(&forward(&edge))?;
        let rebuilt = model.add_edge(fresh, &bounds)?;
        reshape.replace(&forward(&edge), rebuilt);
    }
    let spans = if reshape.is_empty() {
        Built::new(divided.shape.clone(), History::identity())
    } else {
        reshape.apply(model, &divided.shape)?
    };

    // Each face's surface restated as its patch, its edges' trims kept
    // under the patch's id: the patch keeps the parameters they are in.
    let mut reshape = Reshape::new();
    for face in explore_unique(model, &spans.shape, ShapeType::Face)? {
        let data = face_data(model, &face)?;
        let Some(SurfaceGeometry::BSpline(spline)) =
            model.geometry().surface(data.surface).cloned()
        else {
            continue;
        };
        if spline.u_knots().distinct().len() <= 2 && spline.v_knots().distinct().len() <= 2 {
            continue;
        }
        let rings = rings(model, &forward(&face), data.surface, tol)?;
        let ((u0, u1), (v0, v1)) = chart_bounds(&rings, tol);
        let ((du0, du1), (dv0, dv1)) = spline.domain();
        let patch = spline.segment((u0.max(du0), u1.min(du1)), (v0.max(dv0), v1.min(dv1)), tol)?;
        let id = model
            .geometry_mut()
            .add_surface(SurfaceGeometry::BSpline(patch));
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let repr = edge_data(model, &edge)?
                .pcurve_for(data.surface, edge.location())
                .cloned();
            let Some(mut repr) = repr else { continue };
            match &mut repr {
                EdgeRepr::PCurve { surface, .. } | EdgeRepr::Seam { surface, .. } => *surface = id,
                _ => continue,
            }
            let Some(node) = model.node_mut(&edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            if let NodeData::Edge(e) = node.data_mut() {
                let agreed = e.same_parameter();
                e.add(repr);
                e.assert_same_parameter(agreed);
            }
        }
        let mut fresh = data.clone();
        fresh.surface = id;
        fresh.triangulation = None;
        let wires = model.children_of(&forward(&face))?;
        let rebuilt = model.add_face(fresh, &wires)?;
        reshape.replace(&forward(&face), rebuilt);
    }
    let patches = if reshape.is_empty() {
        Built::new(spans.shape.clone(), History::identity())
    } else {
        reshape.apply(model, &spans.shape)?
    };
    Ok(Built::new(
        patches.shape,
        nurbs
            .history
            .then(&divided.history)
            .then(&spans.history)
            .then(&patches.history),
    ))
}

// --- drivers ---------------------------------------------------------------

type Bounds = ((f64, f64), (f64, f64));

/// Cut faces until `lines` asks for no cut any face takes: each round asks
/// every face not yet settled for the lines it wants, in order, and makes
/// the first that crosses it.
fn divide_faces(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
    mut lines: impl FnMut(&Model, &Shape, &SurfaceGeometry, Bounds) -> OgeomResult<Vec<IsoLine>>,
) -> OgeomResult<Built> {
    let mut current = Built::new(shape.clone(), History::identity());
    let mut settled: HashSet<TShapeId> = HashSet::new();
    'rounds: for _ in 0..100_000 {
        for face in explore_unique(model, &current.shape, ShapeType::Face)? {
            if settled.contains(&face.node()) {
                continue;
            }
            let data = face_data(model, &face)?;
            if data.natural_restriction {
                settled.insert(face.node());
                continue;
            }
            let Some(surface) = model.geometry().surface(data.surface).cloned() else {
                ogeom_bail!(Dangling, "surface is not in this model");
            };
            let rings = rings(model, &forward(&face), data.surface, tol)?;
            let bounds = chart_bounds(&rings, tol);
            for line in lines(&*model, &face, &surface, bounds)? {
                let (lo, hi) = match line {
                    IsoLine::U(_) => bounds.0,
                    IsoLine::V(_) => bounds.1,
                };
                let margin = (hi - lo).abs() * 1e-7 + tol.parametric();
                if line.at() <= lo + margin || line.at() >= hi - margin {
                    continue;
                }
                if let Some(reshape) = cut_face(model, &face, line, tol)? {
                    let next = reshape.apply(model, &current.shape)?;
                    current = Built::new(next.shape, current.history.then(&next.history));
                    continue 'rounds;
                }
            }
            settled.insert(face.node());
        }
        return Ok(current);
    }
    ogeom_bail!(Construction, "the division did not settle");
}

/// Cut every edge whose spline is less than `order` times differentiable at
/// an interior knot, there.
fn divide_edges(
    model: &mut Model,
    shape: &Shape,
    order: usize,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut reshape = Reshape::new();
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let edge = forward(&edge);
        let data = edge_data(model, &edge)?;
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let margin = (range.1 - range.0) * 1e-9 + tol.parametric();
        let at: Vec<f64> = curve_knots(geometry)
            .map(|knots| weak_knots(&knots, order))
            .unwrap_or_default()
            .into_iter()
            .filter(|t| *t > range.0 + margin && *t < range.1 - margin)
            .collect();
        if at.is_empty() {
            continue;
        }
        let fractions: Vec<f64> = at
            .iter()
            .map(|t| (t - range.0) / (range.1 - range.0))
            .collect();
        let (pieces, _) = split_edge(model, &edge, &fractions, tol)?;
        reshape.split(&edge, pieces);
    }
    if reshape.is_empty() {
        return Ok(Built::new(shape.clone(), History::identity()));
    }
    reshape.apply(model, shape)
}

/// The shape, its placements baked in where any edge or face is placed:
/// cuts make new nodes in the chart of the node they cut, and a node placed
/// twice has two.
fn placed_baked(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let mut placed = false;
    for kind in [ShapeType::Edge, ShapeType::Face] {
        for s in explore_unique(model, shape, kind)? {
            placed |= !s.location().is_identity();
        }
    }
    for face in explore_unique(model, shape, ShapeType::Face)? {
        placed |= !face_data(model, &face)?.location.is_identity();
    }
    if placed {
        ogeom_algo::baked_shape(model, shape, tol)
    } else {
        Ok(Built::new(shape.clone(), History::identity()))
    }
}

// --- where to cut -----------------------------------------------------------

fn curve_knots(curve: &Curve) -> Option<KnotVector> {
    match curve {
        Curve::BSpline(b) => Some(b.knots().clone()),
        Curve::Trimmed(t) => curve_knots(t.basis()),
        _ => None,
    }
}

/// Interior knots where the spline is less than `order` times
/// differentiable.
fn weak_knots(knots: &KnotVector, order: usize) -> Vec<f64> {
    let (a, b) = knots.domain();
    let p = knots.degree();
    knots
        .distinct()
        .into_iter()
        .filter(|(value, m)| *value > a && *value < b && p.saturating_sub(*m) < order)
        .map(|(value, _)| value)
        .collect()
}

/// The knot lines of a surface where it is less than `order` times
/// differentiable.
fn knot_lines(surface: &SurfaceGeometry, order: usize) -> Vec<IsoLine> {
    match surface {
        SurfaceGeometry::BSpline(b) => weak_knots(b.u_knots(), order)
            .into_iter()
            .map(IsoLine::U)
            .chain(weak_knots(b.v_knots(), order).into_iter().map(IsoLine::V))
            .collect(),
        SurfaceGeometry::Extrusion(e) => curve_knots(e.curve())
            .map(|k| weak_knots(&k, order).into_iter().map(IsoLine::U).collect())
            .unwrap_or_default(),
        SurfaceGeometry::Revolution(r) => curve_knots(r.curve())
            .map(|k| weak_knots(&k, order).into_iter().map(IsoLine::V).collect())
            .unwrap_or_default(),
        SurfaceGeometry::Trimmed(t) => knot_lines(t.basis(), order),
        // An offset is one order rougher than its basis at a knot.
        SurfaceGeometry::Offset(o) => knot_lines(o.basis(), order.saturating_add(1)),
        _ => Vec::new(),
    }
}

/// Whether each parameter of a surface is an angle.
const fn angular(surface: &SurfaceGeometry) -> (bool, bool) {
    match surface {
        SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Revolution(_) => (true, false),
        SurfaceGeometry::Sphere(_) | SurfaceGeometry::Torus(_) => (true, true),
        _ => (false, false),
    }
}

// --- one cut ------------------------------------------------------------------

/// One traversal of an edge around a face, in the chart.
#[derive(Debug, Clone)]
struct Occurrence {
    edge: Shape,
    pcurve: PlanarCurve,
    range: (f64, f64),
}

impl Occurrence {
    fn reversed(&self) -> bool {
        self.edge.orientation() == Orientation::Reversed
    }

    /// The chart point a fraction `s` of the way along the traversal.
    fn at(&self, s: f64, tol: Tolerances) -> OgeomResult<Point2> {
        let s = if self.reversed() { 1.0 - s } else { s };
        self.pcurve
            .point_at(self.range.0 + (self.range.1 - self.range.0) * s, tol)
    }

    fn polyline(&self, samples: usize, tol: Tolerances) -> OgeomResult<Vec<Point2>> {
        (0..=samples)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let s = i as f64 / samples as f64;
                self.at(s, tol)
            })
            .collect()
    }
}

/// An edge's trim on a face: one pcurve, or a seam's two, over a range.
type Trims = (PlanarCurve, Option<PlanarCurve>, (f64, f64));

/// A ring's edges, its chart outline, and the holes it holds.
type Outer = (Vec<Shape>, Vec<Point2>, Vec<Vec<Shape>>);

/// A face's rings as traversals, each seam traversal on the side of the
/// chart its ring continues on.
fn rings(
    model: &Model,
    face: &Shape,
    surface: SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Occurrence>>> {
    let mut out = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let edges = model.ordered_children_of(&wire)?;
        let mut choices: Vec<Option<Trims>> = Vec::new();
        for edge in &edges {
            let data = edge_data(model, edge)?;
            let pcurve =
                |id| -> OgeomResult<PlanarCurve> {
                    model.geometry().pcurve(id).cloned().ok_or_else(|| {
                        ogeom_core::ogeom_err!(Dangling, "pcurve is not in this model")
                    })
                };
            choices.push(match data.pcurve_for(surface, edge.location()) {
                Some(EdgeRepr::PCurve { curve, range, .. }) => {
                    Some((pcurve(*curve)?, None, *range))
                }
                Some(EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                }) => Some((pcurve(*forward)?, Some(pcurve(*reversed)?), *range)),
                _ => None,
            });
        }
        // Walked from an edge off the seam where there is one; a ring of
        // seams alone (a torus's) starts on the side its first occurrence's
        // orientation names.
        let first_plain = choices
            .iter()
            .position(|c| c.as_ref().is_some_and(|c| c.1.is_none()))
            .unwrap_or(0);
        let n = edges.len();
        let mut ring: Vec<Option<Occurrence>> = vec![None; n];
        let mut last: Option<Point2> = None;
        for k in 0..n {
            let i = (first_plain + k) % n;
            let Some((a, b, range)) = choices[i].clone() else {
                ogeom_bail!(Construction, "an edge has no trim on this face");
            };
            let mut occurrence = Occurrence {
                edge: edges[i].clone(),
                pcurve: a,
                range,
            };
            if let Some(b) = b {
                let other = Occurrence {
                    pcurve: b,
                    ..occurrence.clone()
                };
                let take_other = match last {
                    Some(last) => {
                        other.at(0.0, tol)?.distance(last) < occurrence.at(0.0, tol)?.distance(last)
                    }
                    None => occurrence.reversed(),
                };
                if take_other {
                    occurrence = other;
                }
            }
            last = Some(occurrence.at(1.0, tol)?);
            ring[i] = Some(occurrence);
        }
        out.push(ring.into_iter().flatten().collect());
    }
    Ok(out)
}

fn chart_bounds(rings: &[Vec<Occurrence>], tol: Tolerances) -> Bounds {
    let mut b = (
        (f64::INFINITY, f64::NEG_INFINITY),
        (f64::INFINITY, f64::NEG_INFINITY),
    );
    for o in rings.iter().flatten() {
        for i in 0..=32 {
            if let Ok(p) = o.at(f64::from(i) / 32.0, tol) {
                b.0 = (b.0.0.min(p.x), b.0.1.max(p.x));
                b.1 = (b.1.0.min(p.y), b.1.1.max(p.y));
            }
        }
    }
    b
}

/// Where a line meets a face's boundary.
#[derive(Debug, Clone, Copy)]
enum Meeting {
    /// At an existing vertex.
    Vertex(TShapeId),
    /// Inside an edge, a fraction of the way along it (forward).
    Inside(TShapeId, u64),
}

/// The substitutions that cut `face` along `line`, or `None` where the
/// line does not cross it.
#[allow(clippy::too_many_lines)]
fn cut_face(
    model: &mut Model,
    face: &Shape,
    line: IsoLine,
    tol: Tolerances,
) -> OgeomResult<Option<Reshape>> {
    let face_fwd = forward(face);
    let data = face_data(model, &face_fwd)?;
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(Dangling, "surface is not in this model");
    };
    let rings = rings(model, &face_fwd, data.surface, tol)?;
    let bounds = chart_bounds(&rings, tol);
    let span = (bounds.0.1 - bounds.0.0).max(bounds.1.1 - bounds.1.0);
    let eps = span * 1e-9 + tol.parametric();
    let c = line.at();

    // Where the line crosses each edge, as fractions of the edge's forward
    // range; and every point where it meets the boundary, along the line.
    let mut cuts: HashMap<TShapeId, Vec<f64>> = HashMap::new();
    let mut meetings: Vec<(f64, Meeting)> = Vec::new();
    for occurrence in rings.iter().flatten() {
        let f = |s: f64| -> OgeomResult<f64> { Ok(line.fixed(occurrence.at(s, tol)?) - c) };
        const N: usize = 64;
        #[allow(clippy::cast_precision_loss)]
        let values: Vec<f64> = (0..=N)
            .map(|i| f(i as f64 / N as f64))
            .collect::<OgeomResult<_>>()?;
        if values.iter().all(|v| v.abs() <= eps) {
            ogeom_bail!(Construction, "an edge of the face runs along {line:?}");
        }
        let end = occurrence.at(1.0, tol)?;
        if (line.fixed(end) - c).abs() <= eps {
            let Some((_, v)) = edge_vertices(model, &occurrence.edge)? else {
                continue;
            };
            meetings.push((line.free(end), Meeting::Vertex(v.node())));
        }
        for i in 0..N {
            let (a, b) = (values[i], values[i + 1]);
            if a.abs() <= eps || b.abs() <= eps || a.signum() == b.signum() {
                // Through a sample exactly: a crossing there only where the
                // sign changes across it.
                if i + 1 < N
                    && b.abs() <= eps
                    && a.abs() > eps
                    && values[i + 2].abs() > eps
                    && a.signum() != values[i + 2].signum()
                {
                    #[allow(clippy::cast_precision_loss)]
                    let s = (i + 1) as f64 / N as f64;
                    record(occurrence, s, line, &mut cuts, &mut meetings, tol)?;
                }
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let (mut lo, mut hi) = (i as f64 / N as f64, (i + 1) as f64 / N as f64);
            let mut f_lo = a;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let f_mid = f(mid)?;
                if f_mid.signum() == f_lo.signum() {
                    (lo, f_lo) = (mid, f_mid);
                } else {
                    hi = mid;
                }
            }
            record(
                occurrence,
                0.5 * (lo + hi),
                line,
                &mut cuts,
                &mut meetings,
                tol,
            )?;
        }
    }
    meetings.sort_by(|a, b| a.0.total_cmp(&b.0));
    meetings.dedup_by(|b, a| {
        (a.0 - b.0).abs() <= eps
            && match (a.1, b.1) {
                (Meeting::Vertex(x), Meeting::Vertex(y)) => x == y,
                (Meeting::Inside(x, s), Meeting::Inside(y, t)) => x == y && s == t,
                _ => false,
            }
    });

    // The stretches of the line inside the face.
    let outlines: Vec<Vec<Point2>> = rings
        .iter()
        .map(|ring| -> OgeomResult<Vec<Point2>> {
            let mut out = Vec::new();
            for o in ring {
                out.extend(o.polyline(32, tol)?);
            }
            Ok(out)
        })
        .collect::<OgeomResult<_>>()?;
    let mut stretches: Vec<(usize, usize)> = Vec::new();
    for i in 0..meetings.len().saturating_sub(1) {
        let (w0, w1) = (meetings[i].0, meetings[i + 1].0);
        if w1 - w0 <= eps * 10.0 {
            continue;
        }
        if inside(&outlines, line.point(0.5 * (w0 + w1))) {
            stretches.push((i, i + 1));
        }
    }
    if stretches.is_empty() {
        return Ok(None);
    }
    let outer_sign = outlines
        .iter()
        .map(|o| signed_area(o))
        .max_by(|a, b| a.abs().total_cmp(&b.abs()))
        .unwrap_or(1.0)
        .signum();

    // Cut the edges, and every meeting inside one becomes its new vertex.
    let mut reshape = Reshape::new();
    let mut pieces_of: HashMap<TShapeId, Vec<(Shape, f64, f64)>> = HashMap::new();
    let mut vertex_at: HashMap<(TShapeId, u64), Shape> = HashMap::new();
    let mut vertex_of_node: HashMap<TShapeId, Shape> = HashMap::new();
    for v in explore_unique(model, &face_fwd, ShapeType::Vertex)? {
        vertex_of_node.insert(v.node(), v);
    }
    for occurrence in rings.iter().flatten() {
        let edge = forward(&occurrence.edge);
        if pieces_of.contains_key(&edge.node()) {
            continue;
        }
        let Some(fractions) = cuts.get(&edge.node()) else {
            continue;
        };
        let mut fractions = fractions.clone();
        fractions.sort_by(f64::total_cmp);
        fractions.dedup_by(|b, a| (*a - *b).abs() <= 1e-12);
        let (pieces, vertices) = split_edge(model, &edge, &fractions, tol)?;
        for (s, v) in fractions.iter().zip(&vertices) {
            vertex_at.insert((edge.node(), s.to_bits()), v.clone());
        }
        let mut bounds = vec![0.0];
        bounds.extend(&fractions);
        bounds.push(1.0);
        pieces_of.insert(
            edge.node(),
            pieces
                .iter()
                .zip(bounds.windows(2))
                .map(|(p, w)| (p.clone(), w[0], w[1]))
                .collect(),
        );
        reshape.split(&edge, pieces);
    }
    let meeting_vertex = |m: Meeting| -> OgeomResult<Shape> {
        match m {
            Meeting::Vertex(node) => vertex_of_node
                .get(&node)
                .cloned()
                .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a meeting vertex is lost")),
            Meeting::Inside(edge, s) => vertex_at
                .get(&(edge, s))
                .cloned()
                .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a cut vertex is lost")),
        }
    };

    // The boundary's pieces, each on its side of the line.
    let mut sides: [Vec<Item>; 2] = [Vec::new(), Vec::new()];
    for occurrence in rings.iter().flatten() {
        let edge = forward(&occurrence.edge);
        let parts: Vec<(Shape, f64, f64)> = match pieces_of.get(&edge.node()) {
            Some(parts) => parts.clone(),
            None => vec![(edge.clone(), 0.0, 1.0)],
        };
        let mut ordered: Vec<(Shape, f64, f64)> = if occurrence.reversed() {
            parts
                .into_iter()
                .rev()
                .map(|(p, a, b)| (p.reversed(), 1.0 - b, 1.0 - a))
                .collect()
        } else {
            parts
        };
        for (piece, a, b) in ordered.drain(..) {
            let sub = |s: f64| occurrence.at(a + (b - a) * s, tol);
            let polyline: Vec<Point2> = (0..=16)
                .map(|i| sub(f64::from(i) / 16.0))
                .collect::<OgeomResult<_>>()?;
            let side = usize::from(line.fixed(polyline[8]) > c);
            let Some((from, to)) = edge_vertices(model, &piece)? else {
                ogeom_bail!(Construction, "a boundary piece has no vertices");
            };
            sides[side].push(Item {
                edge: piece,
                from: from.node(),
                to: to.node(),
                polyline,
            });
        }
    }

    // The new edges along the line, each walked one way by each side.
    let probes: Vec<f64> = stretches
        .iter()
        .map(|&(i, j)| 0.5 * (meetings[i].0 + meetings[j].0))
        .collect();
    let Some((curve, scale, offset)) = iso_curve(&surface, line, &probes, tol)? else {
        ogeom_bail!(
            Construction,
            "the surface has no closed-form iso-curve to cut along"
        );
    };
    for (i, j) in stretches {
        let (w0, w1) = (meetings[i].0, meetings[j].0);
        let (a, b) = (
            meeting_vertex(meetings[i].1)?,
            meeting_vertex(meetings[j].1)?,
        );
        let (t0, t1) = (scale * w0 + offset, scale * w1 + offset);
        let rising = t1 > t0;
        let (range, (start, end), (p0, p1)) = if rising {
            ((t0, t1), (&a, &b), (line.point(w0), line.point(w1)))
        } else {
            ((t1, t0), (&b, &a), (line.point(w1), line.point(w0)))
        };
        let curve_id = model.geometry_mut().add_curve(curve.clone());
        let mut edge_data = EdgeData::on_curve(curve_id, Location::identity(), range);
        let trim = BSpline2d::new(
            KnotVector::new(vec![range.0, range.0, range.1, range.1], 1)?,
            vec![p0, p1],
            tol,
        )?;
        let trim_id = model.geometry_mut().add_pcurve(PlanarCurve::BSpline(trim));
        edge_data.add(EdgeRepr::PCurve {
            curve: trim_id,
            surface: data.surface,
            location: Location::identity(),
            range,
        });
        edge_data.assert_same_parameter(true);
        // The vertices were placed by the edges they cut; the new edge's
        // ends reach them within the gap.
        let reach = [(start, range.0), (end, range.1)]
            .iter()
            .map(|(v, t)| -> OgeomResult<f64> {
                Ok(vertex_point(model, v)?.distance(curve.point_at(*t, tol)?))
            })
            .collect::<OgeomResult<Vec<f64>>>()?
            .into_iter()
            .fold(0.0_f64, f64::max);
        edge_data.widen(ogeom_core::Tolerance::new(reach + tol.confusion())?);
        let edge = model.add_edge(edge_data, &[start.clone(), end.clone()])?;
        let polyline: Vec<Point2> = (0..=16)
            .map(|k| {
                let w = w0 + (w1 - w0) * f64::from(k) / 16.0;
                line.point(w)
            })
            .collect();
        // Walked up the line (increasing free parameter), the material on
        // the left is the low side of a `u` line and the high side of a `v`
        // line, for a boundary wound counter-clockwise.
        let left = match line {
            IsoLine::U(_) => 0,
            IsoLine::V(_) => 1,
        };
        let up_side = if outer_sign >= 0.0 { left } else { 1 - left };
        let up = if rising {
            edge.clone()
        } else {
            edge.reversed()
        };
        let down = up.reversed();
        sides[up_side].push(Item {
            edge: up,
            from: a.node(),
            to: b.node(),
            polyline: polyline.clone(),
        });
        sides[1 - up_side].push(Item {
            edge: down,
            from: b.node(),
            to: a.node(),
            polyline: polyline.into_iter().rev().collect(),
        });
    }

    // Each side's pieces chained into rings, and the rings into faces.
    let mut faces: Vec<Shape> = Vec::new();
    for items in sides {
        let loops = chain(items, eps.max(tol.parametric() * 1e3))?;
        let mut outers: Vec<Outer> = Vec::new();
        let mut holes: Vec<(Vec<Shape>, Vec<Point2>)> = Vec::new();
        for (edges, outline) in loops {
            if signed_area(&outline) * outer_sign > 0.0 {
                outers.push((edges, outline, Vec::new()));
            } else {
                holes.push((edges, outline));
            }
        }
        for (edges, outline) in holes {
            let probe = outline[outline.len() / 2];
            let Some(host) = outers
                .iter_mut()
                .filter(|o| inside(std::slice::from_ref(&o.1), probe))
                .min_by(|a, b| signed_area(&a.1).abs().total_cmp(&signed_area(&b.1).abs()))
            else {
                ogeom_bail!(Construction, "a hole of the cut face lies in no piece");
            };
            host.2.push(edges);
        }
        for (edges, _, hole_edges) in outers {
            let mut wires = vec![model.add_wire(&edges)?];
            for h in hole_edges {
                wires.push(model.add_wire(&h)?);
            }
            let mut fresh = data.clone();
            fresh.triangulation = None;
            fresh.natural_restriction = false;
            faces.push(model.add_face(fresh, &wires)?);
        }
    }
    if faces.len() < 2 {
        return Ok(None);
    }
    reshape.split(&face_fwd, faces);
    Ok(Some(reshape))
}

/// Record a crossing a fraction `s` along `occurrence`.
fn record(
    occurrence: &Occurrence,
    s: f64,
    line: IsoLine,
    cuts: &mut HashMap<TShapeId, Vec<f64>>,
    meetings: &mut Vec<(f64, Meeting)>,
    tol: Tolerances,
) -> OgeomResult<()> {
    let at = occurrence.at(s, tol)?;
    // As a fraction of the forward edge, rounded so two traversals of one
    // seam name the same point.
    let forward_s = if occurrence.reversed() { 1.0 - s } else { s };
    let forward_s = (forward_s * 1e12).round() / 1e12;
    let node = occurrence.edge.node();
    cuts.entry(node).or_default().push(forward_s);
    meetings.push((line.free(at), Meeting::Inside(node, forward_s.to_bits())));
    Ok(())
}

/// One directed piece of a ring being assembled.
#[derive(Debug, Clone)]
struct Item {
    edge: Shape,
    from: TShapeId,
    to: TShapeId,
    polyline: Vec<Point2>,
}

/// Chain directed pieces end to start into closed rings: at a vertex, the
/// piece whose chart start is nearest where the ring stands.
fn chain(mut items: Vec<Item>, snap: f64) -> OgeomResult<Vec<(Vec<Shape>, Vec<Point2>)>> {
    let mut loops = Vec::new();
    while let Some(first) = items.pop() {
        let start_node = first.from;
        let start_at = first.polyline[0];
        let mut here = *first.polyline.last().unwrap_or(&start_at);
        let mut at_node = first.to;
        let mut edges = vec![first.edge];
        let mut outline = first.polyline;
        loop {
            if at_node == start_node && here.distance(start_at) <= snap {
                break;
            }
            let next = items
                .iter()
                .enumerate()
                .filter(|(_, it)| it.from == at_node)
                .min_by(|a, b| {
                    a.1.polyline[0]
                        .distance(here)
                        .total_cmp(&b.1.polyline[0].distance(here))
                })
                .map(|(i, _)| i);
            let Some(i) = next else {
                ogeom_bail!(Construction, "a piece of the cut face does not close");
            };
            let it = items.swap_remove(i);
            here = *it.polyline.last().unwrap_or(&here);
            at_node = it.to;
            edges.push(it.edge);
            outline.extend(it.polyline);
        }
        loops.push((edges, outline));
    }
    Ok(loops)
}

/// Split `edge` (forward) at fractions of its range, in order: the pieces
/// in order, and the new vertex at each fraction. Every description the
/// edge has is cut at the same fraction of its own range.
fn split_edge(
    model: &mut Model,
    edge: &Shape,
    fractions: &[f64],
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, Vec<Shape>)> {
    let data = edge_data(model, edge)?;
    let Some((start, end)) = edge_vertices(model, edge)? else {
        ogeom_bail!(Construction, "an edge with no vertices cannot be cut");
    };
    let mut vertices = Vec::with_capacity(fractions.len());
    for &s in fractions {
        let vertex = match data.curve3d() {
            Some(EdgeRepr::Curve3d { curve, range, .. }) if !data.degenerate => {
                let Some(curve) = model.geometry().curve(*curve) else {
                    ogeom_bail!(Dangling, "curve is not in this model");
                };
                let p = curve.point_at(range.0 + (range.1 - range.0) * s, tol)?;
                model.add_vertex(VertexData::with_tolerance(p, data.tolerance.get())?)
            }
            // A degenerate edge is one point, however it is cut.
            _ => start.clone(),
        };
        vertices.push(vertex);
    }
    let mut bounds = vec![0.0];
    bounds.extend(fractions);
    bounds.push(1.0);
    let mut ends = vec![start];
    ends.extend(vertices.iter().cloned());
    ends.push(end);
    let mut pieces = Vec::with_capacity(bounds.len() - 1);
    for k in 0..bounds.len() - 1 {
        let (a, b) = (bounds[k], bounds[k + 1]);
        let sub = |range: (f64, f64)| {
            (
                range.0 + (range.1 - range.0) * a,
                range.0 + (range.1 - range.0) * b,
            )
        };
        let mut piece = data.clone();
        piece.representations.retain(|r| {
            matches!(
                r,
                EdgeRepr::Curve3d { .. } | EdgeRepr::PCurve { .. } | EdgeRepr::Seam { .. }
            )
        });
        for repr in &mut piece.representations {
            match repr {
                EdgeRepr::Curve3d { range, .. }
                | EdgeRepr::PCurve { range, .. }
                | EdgeRepr::Seam { range, .. } => *range = sub(*range),
                _ => {}
            }
        }
        pieces.push(model.add_edge(piece, &[ends[k].clone(), ends[k + 1].clone()])?);
    }
    Ok((pieces, vertices))
}

/// The surface's iso-curve along `line`, and the map from the line's free
/// parameter `w` to the curve's own: `t = scale * w + offset`.
fn iso_curve(
    surface: &SurfaceGeometry,
    line: IsoLine,
    probes: &[f64],
    tol: Tolerances,
) -> OgeomResult<Option<(Curve, f64, f64)>> {
    let c = line.at();
    let circle = |centre: Point, z: Direction, x: Point, radius: f64| -> OgeomResult<Curve> {
        let x = Direction::new(x - centre, tol)?;
        let frame = Frame::new(centre, z, x, tol)?;
        Ok(CircleCurve::new(Circle::new(frame, radius, tol)?).into())
    };
    let straight = |at: Point, along: ogeom_math::Vector| -> OgeomResult<(Curve, f64)> {
        let length = along.magnitude();
        let direction = Direction::new(along, tol)?;
        Ok((
            LineCurve::new(Axis {
                location: at,
                direction,
            })
            .into(),
            length,
        ))
    };
    let found: (Curve, f64, f64) = match (surface, line) {
        (SurfaceGeometry::BSpline(b), IsoLine::U(_)) => {
            (Curve::BSpline(b.iso_u_curve(c, tol)?), 1.0, 0.0)
        }
        (SurfaceGeometry::BSpline(b), IsoLine::V(_)) => {
            (Curve::BSpline(b.iso_v_curve(c, tol)?), 1.0, 0.0)
        }
        (
            SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_),
            IsoLine::U(_),
        )
        | (SurfaceGeometry::Plane(_) | SurfaceGeometry::Extrusion(_), _) => {
            // Straight in the free parameter: through the point at w = 0,
            // along the rate the surface moves with w.
            let at = surface.point_at(line.point(0.0).x, line.point(0.0).y, tol)?;
            let next = surface.point_at(line.point(1.0).x, line.point(1.0).y, tol)?;
            match (surface, line) {
                (SurfaceGeometry::Extrusion(e), IsoLine::V(_)) => {
                    let shift = Transform::translation(e.direction().vector() * c);
                    (e.curve().transformed(&shift, tol)?, 1.0, 0.0)
                }
                _ => {
                    let (curve, speed) = straight(at, next - at)?;
                    (curve, speed, 0.0)
                }
            }
        }
        (
            SurfaceGeometry::Cylinder(_)
            | SurfaceGeometry::Cone(_)
            | SurfaceGeometry::Sphere(_)
            | SurfaceGeometry::Torus(_)
            | SurfaceGeometry::Revolution(_),
            IsoLine::V(_),
        ) => {
            // A parallel: the circle the surface's `u = 0` point turns on.
            let axis = match surface {
                SurfaceGeometry::Cylinder(s) => frame_axis(s.cylinder().frame()),
                SurfaceGeometry::Cone(s) => frame_axis(s.cone().frame()),
                SurfaceGeometry::Sphere(s) => frame_axis(s.sphere().frame()),
                SurfaceGeometry::Torus(s) => frame_axis(s.torus().frame()),
                SurfaceGeometry::Revolution(s) => s.axis(),
                _ => return Ok(None),
            };
            let start = surface.point_at(0.0, c, tol)?;
            let centre = axis.project(start);
            let radius = start.distance(centre);
            if radius <= tol.confusion() {
                return Ok(None);
            }
            (circle(centre, axis.direction, start, radius)?, 1.0, 0.0)
        }
        (SurfaceGeometry::Sphere(_) | SurfaceGeometry::Torus(_), IsoLine::U(_)) => {
            // A meridian or tube circle, its angle the surface's `v`.
            let at0 = surface.point_at(c, 0.0, tol)?;
            let quarter = surface.point_at(c, core::f64::consts::FRAC_PI_2, tol)?;
            let opposite = surface.point_at(c, -core::f64::consts::FRAC_PI_2, tol)?;
            let centre = Point::from_vector((quarter.to_vector() + opposite.to_vector()) * 0.5);
            let radius = at0.distance(centre);
            let z = Direction::new((at0 - centre).cross(quarter - centre), tol)?;
            (circle(centre, z, at0, radius)?, 1.0, 0.0)
        }
        (SurfaceGeometry::Revolution(r), IsoLine::U(_)) => {
            let turn = Transform::rotation(r.axis(), c);
            (r.curve().transformed(&turn, tol)?, 1.0, 0.0)
        }
        (SurfaceGeometry::Trimmed(t), _) => return iso_curve(t.basis(), line, probes, tol),
        _ => return Ok(None),
    };
    // The closed forms above assume the surfaces' own conventions; hold
    // them to it.
    let (curve, scale, offset) = &found;
    for &w in probes {
        let on = line.point(w);
        let expected = surface.point_at(on.x, on.y, tol)?;
        let got = curve.point_at(scale * w + offset, tol)?;
        if expected.distance(got) > tol.confusion() {
            return Ok(None);
        }
    }
    Ok(Some(found))
}

fn frame_axis(frame: Frame) -> Axis {
    Axis {
        location: frame.origin(),
        direction: frame.z(),
    }
}

// --- helpers -------------------------------------------------------------------

fn forward(shape: &Shape) -> Shape {
    if shape.orientation() == Orientation::Reversed {
        shape.reversed()
    } else {
        shape.clone()
    }
}

fn edge_data(model: &Model, edge: &Shape) -> OgeomResult<EdgeData> {
    match model.node(edge).map(ogeom_topo::TShape::data) {
        Some(NodeData::Edge(data)) => Ok((**data).clone()),
        Some(_) => ogeom_bail!(Construction, "expected an edge"),
        None => ogeom_bail!(Dangling, "edge is not in this model"),
    }
}

fn face_data(model: &Model, face: &Shape) -> OgeomResult<ogeom_topo::FaceData> {
    match model.node(face).map(ogeom_topo::TShape::data) {
        Some(NodeData::Face(data)) => Ok((**data).clone()),
        Some(_) => ogeom_bail!(Construction, "expected a face"),
        None => ogeom_bail!(Dangling, "face is not in this model"),
    }
}

fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Construction, "a vertex holds no point");
    };
    Ok(data.point)
}

fn signed_area(ring: &[Point2]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            a.x * b.y - b.x * a.y
        })
        .sum::<f64>()
        * 0.5
}

/// Even-odd containment against every ring.
fn inside(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut odd = false;
    for ring in rings {
        let n = ring.len();
        for i in 0..n {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > p.x {
                    odd = !odd;
                }
            }
        }
    }
    odd
}
