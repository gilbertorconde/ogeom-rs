//! Splitting a face along curves: edges lying on it, or dropped onto it
//! along its normals or along a direction.
//!
//! Each curve is first given its trace in the face's chart: its own pcurve
//! where it carries one on the surface, the exact or fitted one where it
//! lies on the surface, or, for a projection, the curve through its landed
//! points fitted together with their chart positions, so curve and trace
//! share a parameter. A trace leaving the face's window of a periodic chart
//! is cut there and brought back by a period, and the pieces of all the
//! curves are chained end to end into paths.
//!
//! Each path is then cut where it meets the face's boundary, and each
//! stretch of it inside the face divides the face it runs through: a
//! stretch from one point of a boundary loop to another point of the same
//! loop cuts that loop in two, and a path closing on itself inside the face
//! cuts out the region it encloses, leaving a hole in the rest. A closed
//! path that meets the boundary is walked from a meeting, so it is a
//! stretch from the boundary to the boundary. Stretches are cut one at a
//! time, so a later path meets the pieces of the earlier ones as boundary
//! and may end on them.
//!
//! The boundary edges a stretch ends on are cut there in every face that
//! holds them, so the rest of the shape stays joined to the pieces.

use ogeom_core::FastSet;

use ogeom_algo::{Built, edge_vertices, project_on_surface, project_on_surface_from};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail, ogeom_err};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry,
    Transformable as _,
};
use ogeom_math::{Direction, Point, Point2, Transform2, Vector2};
use ogeom_topo::{
    EdgeData, EdgeRepr, Location, Model, Shape, ShapeType, SurfaceId, TShapeId, VertexData,
    explore_unique,
};

use crate::Reshape;
use crate::divide::{
    Bounds, Occurrence, bounded_natural, chart_bounds, edge_data, face_data, forward, inside,
    placed_baked, rings, signed_area, split_edge, vertex_point,
};

/// How the curves to split along reach the face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Projection {
    /// The curves lie on the face's surface already.
    OnFace,
    /// Each point of a curve is taken to its nearest point on the surface.
    AlongNormals,
    /// Each point of a curve is carried along the direction, either way,
    /// to where it meets the surface.
    Along(Direction),
}

/// Cut `face` of `shape` along `edges` and rebuild `shape` around the
/// pieces.
///
/// `edges` are edges, or wires and compounds of them, taken onto the face
/// as `projection` says; curves meeting end to end are followed as one.
/// Every stretch of a curve inside the face must run from the face's
/// boundary to its boundary, the boundary including the stretches cut
/// before it, or close on itself inside the face; a curve that runs past
/// the face is cut where it leaves it. The pieces share the new edges, and
/// the boundary edges a cut ends on are cut there in every face that holds
/// them. A face given on its own comes back as a shell of its pieces.
///
/// The history records `face` modified into its pieces and each edge in
/// `edges` generating the edges cut along it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// reason, where `face` is not a face of `shape`, no edge is given, an edge
/// asked to lie on the face does not, a curve or its projection misses the
/// face, a curve ends inside the face, runs along its boundary, or runs from
/// one boundary loop of it to another.
pub fn split_face(
    model: &mut Model,
    shape: &Shape,
    face: &Shape,
    edges: &[Shape],
    projection: Projection,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "only a face can be split along curves");
    }
    if !explore_unique(model, shape, ShapeType::Face)?
        .iter()
        .any(|f| f.node() == face.node())
    {
        ogeom_bail!(Construction, "the face to split is not a face of the shape");
    }
    let mut sources: Vec<Shape> = Vec::new();
    for given in edges {
        for edge in explore_unique(model, given, ShapeType::Edge)? {
            if !edge_data(model, &edge)?.degenerate {
                sources.push(edge);
            }
        }
    }
    if sources.is_empty() {
        ogeom_bail!(Construction, "a face is split along at least one edge");
    }

    // A face given alone is held in a shell, so its pieces have a home.
    let root = if model.kind_of(shape)? == ShapeType::Face {
        model.add_shell(std::slice::from_ref(shape))?
    } else {
        shape.clone()
    };
    let start = placed_baked(model, &root, tol)?;
    let mut current = Built::new(start.shape.clone(), start.history);
    let Some(mut target) = current.history.trace(face).first().cloned() else {
        ogeom_bail!(Construction, "the face to split is lost in its placement");
    };
    if face_data(model, &target)?.natural_restriction {
        let Some(bounded) = bounded_natural(model, &target, tol)? else {
            ogeom_bail!(
                Construction,
                "a face covering an unbounded surface has no boundary to cut to"
            );
        };
        let mut reshape = Reshape::new();
        reshape.replace(&forward(&target), bounded.clone());
        let next = reshape.apply(model, &current.shape)?;
        current = Built::new(next.shape, current.history.then(&next.history));
        target = bounded;
    }
    let data = face_data(model, &target)?;
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(Dangling, "surface is not in this model");
    };
    ogeom_algo::build::trimmed_where_bare(model, &target, tol)?;
    let window = chart_bounds(&rings(model, &forward(&target), data.surface, tol)?, tol);

    let mut segments: Vec<Segment> = Vec::new();
    for source in &sources {
        let found = segments_of(model, source, &surface, data.surface, projection, tol)?;
        if found.is_empty() {
            ogeom_bail!(Construction, "the projection of an edge misses the face");
        }
        for segment in found {
            segments.extend(unfold(segment, &surface, window, tol)?);
        }
    }
    let mut cutters = chain(segments, &surface, tol)?;

    let mut targets: Vec<Shape> = vec![target];
    let mut pending: Vec<usize> = (0..cutters.len()).collect();
    loop {
        let mut progressed = false;
        let mut held: Vec<usize> = Vec::new();
        for &k in &pending {
            let dangling = loop {
                let mut found = None;
                let mut dangling = false;
                for face in &targets {
                    match next_stretch(model, face, &cutters[k], &surface, tol)? {
                        Found::Nothing => {}
                        Found::Dangling => dangling = true,
                        Found::Rebase(t) => {
                            found = Some(Next::Rebase(t));
                            break;
                        }
                        Found::Cut(stretch) => {
                            found = Some(Next::Cut(face.clone(), Box::new(stretch)));
                            break;
                        }
                    }
                }
                let (face, stretch) = match found {
                    None => break dangling,
                    Some(Next::Rebase(t)) => {
                        cutters[k] = cutters[k].rebased(t, tol)?;
                        continue;
                    }
                    Some(Next::Cut(face, stretch)) => (face, *stretch),
                };
                let made = cut(model, &face, &cutters[k], &stretch, tol)?;
                let mut step = made.reshape.apply(model, &current.shape)?;
                for (source, edge) in &made.edges {
                    step.history.generate(source, edge.clone());
                }
                let mut next_targets = Vec::new();
                for t in &targets {
                    next_targets.extend(step.history.trace(t).iter().cloned());
                }
                targets = next_targets;
                for cutter in &mut cutters {
                    for (old, pieces) in &made.splits {
                        if cutter.own.remove(old) {
                            cutter.own.extend(pieces.iter().map(Shape::node));
                        }
                    }
                }
                cutters[k]
                    .own
                    .extend(made.edges.iter().map(|(_, e)| e.node()));
                cutters[k].consumed.push(stretch.span);
                current = Built::new(step.shape, current.history.then(&step.history));
                progressed = true;
            };
            if dangling {
                held.push(k);
            } else if cutters[k].consumed.is_empty() {
                ogeom_bail!(
                    Construction,
                    "a curve to split along does not cross the face"
                );
            }
        }
        if held.is_empty() {
            break;
        }
        if !progressed {
            ogeom_bail!(
                Construction,
                "a curve to split along ends inside the face; it must run from \
                 boundary to boundary"
            );
        }
        pending = held;
    }
    Ok(current)
}

/// What the search for a path's next stretch turned up.
enum Next {
    Rebase(f64),
    Cut(Shape, Box<Planned>),
}

// --- paths ---------------------------------------------------------------------

/// A piece of a path: a curve over a range of its parameter, with its trace
/// on the surface at the same parameter.
#[derive(Debug, Clone)]
struct Segment {
    /// The edge it came from.
    source: Shape,
    curve: Curve,
    pcurve: PlanarCurve,
    /// The range it spans, increasing.
    range: (f64, f64),
    /// Whether the path walks it from `range.1` to `range.0`.
    backwards: bool,
    /// How far the curve and the surface along its trace part.
    tolerance: f64,
}

impl Segment {
    /// The curve's parameter a fraction `f` along the walk.
    fn local(&self, f: f64) -> f64 {
        let (a, b) = if self.backwards {
            (self.range.1, self.range.0)
        } else {
            self.range
        };
        (b - a).mul_add(f, a)
    }

    fn chart(&self, f: f64, tol: Tolerances) -> OgeomResult<Point2> {
        self.pcurve.point_at(self.local(f), tol)
    }

    fn point(&self, f: f64, tol: Tolerances) -> OgeomResult<Point> {
        self.curve.point_at(self.local(f), tol)
    }

    /// The part between fractions `f0 < f1` of the walk.
    fn part(&self, f0: f64, f1: f64) -> Self {
        let (a, b) = (self.local(f0), self.local(f1));
        Self {
            range: (a.min(b), a.max(b)),
            ..self.clone()
        }
    }

    fn reversed(mut self) -> Self {
        self.backwards = !self.backwards;
        self
    }

    /// The segment with its trace moved across the chart by `by`.
    fn shifted(mut self, by: Vector2, tol: Tolerances) -> OgeomResult<Self> {
        if by.x != 0.0 || by.y != 0.0 {
            self.pcurve = self.pcurve.transformed(&Transform2::translation(by), tol)?;
        }
        Ok(self)
    }
}

/// A path to cut along: segments end to end, walked by a parameter that
/// runs from `k` to `k + 1` along segment `k`.
struct Cutter {
    path: Vec<Segment>,
    closed: bool,
    /// The edges already cut along it: its own stretches, never a boundary
    /// it crosses.
    own: FastSet<TShapeId>,
    /// The parameter spans already cut.
    consumed: Vec<(f64, f64)>,
}

impl Cutter {
    #[expect(clippy::cast_precision_loss, reason = "a segment count")]
    fn range(&self) -> (f64, f64) {
        (0.0, self.path.len() as f64)
    }

    /// The segment the parameter `t` falls in, and how far along it.
    fn locate(&self, t: f64) -> (&Segment, f64) {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a parameter clamped to the segment count"
        )]
        let k = (t.floor().max(0.0) as usize).min(self.path.len() - 1);
        #[expect(clippy::cast_precision_loss, reason = "a segment index")]
        let f = t - k as f64;
        (&self.path[k], f)
    }

    fn chart(&self, t: f64, tol: Tolerances) -> OgeomResult<Point2> {
        let (s, f) = self.locate(t);
        s.chart(f, tol)
    }

    fn point(&self, t: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (s, f) = self.locate(t);
        s.point(f, tol)
    }

    /// The trace's derivative by the path's parameter.
    fn chart_d1(&self, t: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        let (s, f) = self.locate(t);
        let rate = if s.backwards {
            s.range.0 - s.range.1
        } else {
            s.range.1 - s.range.0
        };
        Ok(s.pcurve.d1_at(s.local(f), tol)? * rate)
    }

    fn tolerance(&self) -> f64 {
        self.path.iter().map(|s| s.tolerance).fold(0.0, f64::max)
    }

    fn resolution(&self, tol: Tolerances) -> f64 {
        (self.range().1 - self.range().0) * 1e-9 + tol.parametric()
    }

    fn is_consumed(&self, a: f64, b: f64) -> bool {
        let mid = 0.5 * (a + b);
        self.consumed.iter().any(|(lo, hi)| mid > *lo && mid < *hi)
    }

    /// The closed path walked from `t` round to `t`, open.
    fn rebased(&self, t: f64, tol: Tolerances) -> OgeomResult<Self> {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a parameter clamped to the segment count"
        )]
        let k = (t.floor().max(0.0) as usize).min(self.path.len() - 1);
        #[expect(clippy::cast_precision_loss, reason = "a segment index")]
        let f = t - k as f64;
        let eps = self.resolution(tol);
        let mut path = Vec::with_capacity(self.path.len() + 1);
        if f < 1.0 - eps {
            path.push(self.path[k].part(f, 1.0));
        }
        path.extend(self.path[k + 1..].iter().cloned());
        path.extend(self.path[..k].iter().cloned());
        if f > eps {
            path.push(self.path[k].part(0.0, f));
        }
        Ok(Self {
            path,
            closed: false,
            own: FastSet::default(),
            consumed: Vec::new(),
        })
    }
}

/// The segments `edge` gives on the surface.
fn segments_of(
    model: &Model,
    edge: &Shape,
    surface: &SurfaceGeometry,
    surface_id: SurfaceId,
    projection: Projection,
    tol: Tolerances,
) -> OgeomResult<Vec<Segment>> {
    let data = edge_data(model, edge)?;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
        ogeom_bail!(Construction, "an edge to split along has no curve in space");
    };
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let unplaced = edge.location().is_identity();
    let geometry = if unplaced {
        geometry
    } else {
        geometry.transformed(&edge.transform(model.datums())?, tol)?
    };
    let mut out = Vec::new();
    match projection {
        Projection::OnFace => {
            let cap = data.tolerance.get().max(tol.confusion());
            let mut worst: f64 = 0.0;
            for k in 0..=32 {
                let t = (range.1 - range.0).mul_add(f64::from(k) / 32.0, range.0);
                let p = geometry.point_at(t, tol)?;
                worst = worst.max(project_on_surface(surface, p, 24, tol)?.distance);
            }
            if worst > cap {
                ogeom_bail!(
                    Construction,
                    "an edge to split along does not lie on the face: it stands \
                     {worst:.3e} off the surface; project it instead"
                );
            }
            let own = if unplaced && data.same_parameter() {
                match data.pcurve_for(surface_id, edge.location()) {
                    Some(EdgeRepr::PCurve {
                        curve, range: r, ..
                    }) if *r == range => model.geometry().pcurve(*curve).cloned(),
                    _ => None,
                }
            } else {
                None
            };
            let (pcurve, fit) = match own
                .or_else(|| ogeom_intersect::exact_pcurve_over(&geometry, range, surface, tol))
            {
                Some(p) => (p, 0.0),
                None => {
                    let (p, error, _, off, _) =
                        ogeom_algo::pcurve_fit::fit_projected_pcurve_capped(
                            &geometry, range, surface, cap, tol,
                        )?;
                    (p, error.max(off))
                }
            };
            out.push(segment(edge, geometry, pcurve, range, fit, surface, tol)?);
        }
        Projection::AlongNormals | Projection::Along(_) => {
            for (span, seed) in landed_runs(&geometry, range, surface, projection, tol)? {
                out.push(projected(
                    edge, &geometry, span, seed, surface, projection, tol,
                )?);
            }
        }
    }
    Ok(out)
}

/// A segment, its tolerance measured: how far the curve and the surface
/// along its trace part, and no less than the fit's own error.
fn segment(
    edge: &Shape,
    curve: Curve,
    pcurve: PlanarCurve,
    range: (f64, f64),
    fit: f64,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Segment> {
    let mut gap = fit;
    for k in 0..=64 {
        let t = (range.1 - range.0).mul_add(f64::from(k) / 64.0, range.0);
        let uv = pcurve.point_at(t, tol)?;
        let on = surface.point_at(uv.x, uv.y, tol)?;
        gap = gap.max(on.distance(curve.point_at(t, tol)?));
    }
    Ok(Segment {
        source: edge.clone(),
        curve,
        pcurve,
        range,
        backwards: false,
        tolerance: gap,
    })
}

/// The periods of the surface's chart, zero where a direction is not
/// periodic.
fn periods(surface: &SurfaceGeometry) -> (f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.domain();
    (
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    )
}

/// The segment cut where its trace crosses a line a whole number of periods
/// from the face's window, each part moved by the periods that bring it
/// into the window.
fn unfold(
    segment: Segment,
    surface: &SurfaceGeometry,
    window: Bounds,
    tol: Tolerances,
) -> OgeomResult<Vec<Segment>> {
    let (pu, pv) = periods(surface);
    let mut parts = vec![segment];
    for (axis, period, low) in [(0, pu, window.0.0), (1, pv, window.1.0)] {
        if period <= 0.0 || !low.is_finite() {
            continue;
        }
        // Which copy of the window a chart point is in.
        let copy = |p: Point2| -> f64 {
            let x = if axis == 0 { p.x } else { p.y };
            ((x - low) / period + 1e-9).floor()
        };
        let mut next = Vec::new();
        for part in parts {
            const N: usize = 128;
            let mut cuts = vec![0.0];
            let mut last = copy(part.chart(0.0, tol)?);
            for i in 1..=N {
                #[expect(clippy::cast_precision_loss, reason = "a sample index")]
                let f = i as f64 / N as f64;
                let here = copy(part.chart(f, tol)?);
                if here != last {
                    #[expect(clippy::cast_precision_loss, reason = "a sample index")]
                    let (mut a, mut b) = ((i - 1) as f64 / N as f64, f);
                    for _ in 0..60 {
                        let mid = 0.5 * (a + b);
                        if copy(part.chart(mid, tol)?) == last {
                            a = mid;
                        } else {
                            b = mid;
                        }
                    }
                    cuts.push(0.5 * (a + b));
                    last = here;
                }
            }
            cuts.push(1.0);
            // A crossing within the confusion distance of a neighbouring
            // one, or of an end (a closed curve starting on a seam the
            // window's edge stands a hair off), cuts no piece: the part
            // beside it takes the sliver.
            let mut kept = vec![0.0];
            for &c in &cuts[1..] {
                let from = kept[kept.len() - 1];
                let sliver = c - from <= 1e-6
                    && part.point(from, tol)?.distance(part.point(c, tol)?) <= tol.confusion();
                if !sliver {
                    kept.push(c);
                } else if c >= 1.0 {
                    if kept.len() > 1 {
                        let last = kept.len() - 1;
                        kept[last] = 1.0;
                    } else {
                        kept.push(1.0);
                    }
                }
            }
            for w in kept.windows(2) {
                if w[1] - w[0] <= 1e-9 {
                    continue;
                }
                let piece = part.part(w[0], w[1]);
                let turns = copy(piece.chart(0.5, tol)?);
                let by = if axis == 0 {
                    Vector2::new(-turns * period, 0.0)
                } else {
                    Vector2::new(0.0, -turns * period)
                };
                next.push(piece.shifted(by, tol)?);
            }
        }
        parts = next;
    }
    Ok(parts)
}

/// The segments chained end to end into paths: two ends join where they
/// meet in space and their traces meet in the chart, not a period apart.
fn chain(
    mut segments: Vec<Segment>,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Vec<Cutter>> {
    let (pu, pv) = periods(surface);
    let meet = |a: &Segment, fa: f64, b: &Segment, fb: f64| -> OgeomResult<bool> {
        let reach = a.tolerance + b.tolerance + tol.confusion();
        if a.point(fa, tol)?.distance(b.point(fb, tol)?) > reach {
            return Ok(false);
        }
        let (p, q) = (a.chart(fa, tol)?, b.chart(fb, tol)?);
        Ok((pu <= 0.0 || (p.x - q.x).abs() < 0.25 * pu)
            && (pv <= 0.0 || (p.y - q.y).abs() < 0.25 * pv))
    };
    let mut out = Vec::new();
    while let Some(first) = segments.pop() {
        let mut path = std::collections::VecDeque::from([first]);
        loop {
            let mut grown = false;
            let mut k = 0;
            while k < segments.len() {
                let (head, tail) = (&path[0], &path[path.len() - 1]);
                let s = &segments[k];
                if meet(tail, 1.0, s, 0.0)? {
                    path.push_back(segments.swap_remove(k));
                } else if meet(tail, 1.0, s, 1.0)? {
                    path.push_back(segments.swap_remove(k).reversed());
                } else if meet(head, 0.0, s, 1.0)? {
                    path.push_front(segments.swap_remove(k));
                } else if meet(head, 0.0, s, 0.0)? {
                    path.push_front(segments.swap_remove(k).reversed());
                } else {
                    k += 1;
                    continue;
                }
                grown = true;
            }
            if !grown {
                break;
            }
        }
        let path: Vec<Segment> = path.into_iter().collect();
        let closed = meet(&path[path.len() - 1], 1.0, &path[0], 0.0)?;
        out.push(Cutter {
            path,
            closed,
            own: FastSet::default(),
            consumed: Vec::new(),
        });
    }
    Ok(out)
}

// --- projection --------------------------------------------------------------

/// Where `p` lands on the surface, by `projection`, from a nearby landing
/// where one is known: the point and its chart position. `None` where it
/// does not land: no hit along the direction, or a nearest point held at the
/// edge of the surface's window, which is not where the normal through `p`
/// meets the surface.
fn land(
    surface: &SurfaceGeometry,
    p: Point,
    seed: Option<Point2>,
    projection: Projection,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point2)>> {
    match projection {
        Projection::OnFace | Projection::AlongNormals => {
            let found = match seed
                .and_then(|s| project_on_surface_from(surface, p, (s.x, s.y), tol).ok())
            {
                Some(found) => found,
                None => project_on_surface(surface, p, 24, tol)?,
            };
            let (u, v) = found.parameters;
            let gap = p - found.point;
            let d = gap.magnitude();
            if d > tol.confusion() {
                let (su, sv) = surface.d1_at(u, v, tol)?;
                for t in [su, sv] {
                    let m = t.magnitude();
                    if m > 0.0 && gap.dot(t).abs() > 1e-7 * m * d {
                        return Ok(None);
                    }
                }
            }
            Ok(Some((found.point, Point2::new(u, v))))
        }
        Projection::Along(direction) => {
            let start = match seed {
                Some(s) => (s.x, s.y),
                None => project_on_surface(surface, p, 24, tol)?.parameters,
            };
            Ok(ray_hit(surface, p, direction.vector(), start, tol)?
                .map(|(at, u, v)| (at, Point2::new(u, v))))
        }
    }
}

/// Newton on `S(u, v) = p + s d` from `start`, the parameters held to the
/// surface's window. `None` where it does not converge to a hit inside it.
fn ray_hit(
    surface: &SurfaceGeometry,
    p: Point,
    d: ogeom_math::Vector,
    start: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<(Point, f64, f64)>> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let hold = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic && hi > lo {
            lo + (x - lo).rem_euclid(hi - lo)
        } else {
            x.clamp(lo, hi)
        }
    };
    let (mut u, mut v) = start;
    let mut s = (surface.point_at(u, v, tol)? - p).dot(d);
    for _ in 0..64 {
        let (at, su, sv) = surface.point_d1_at(u, v, tol)?;
        let f = at - (p + d * s);
        if f.magnitude() <= tol.confusion() * 1e-3 {
            return Ok(Some((at, u, v)));
        }
        let j = nalgebra::Matrix3::new(su.x, sv.x, -d.x, su.y, sv.y, -d.y, su.z, sv.z, -d.z);
        let Some(step) = j.lu().solve(&nalgebra::Vector3::new(-f.x, -f.y, -f.z)) else {
            return Ok(None);
        };
        u = hold(u + step[0], u0, u1, surface.is_periodic_u());
        v = hold(v + step[1], v0, v1, surface.is_periodic_v());
        s += step[2];
    }
    let at = surface.point_at(u, v, tol)?;
    if at.distance(p + d * s) <= tol.confusion() {
        return Ok(Some((at, u, v)));
    }
    Ok(None)
}

/// The parameter spans of the curve that land on the surface, each run's
/// ends pressed to where the landing stops, with the chart point its start
/// landed on.
fn landed_runs(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    projection: Projection,
    tol: Tolerances,
) -> OgeomResult<Vec<((f64, f64), Point2)>> {
    const STATIONS: usize = 128;
    let mut stations: Vec<(f64, Option<Point2>)> = Vec::with_capacity(STATIONS + 1);
    let mut seed = None;
    for k in 0..=STATIONS {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = (range.1 - range.0).mul_add(k as f64 / STATIONS as f64, range.0);
        let landed = land(surface, curve.point_at(t, tol)?, seed, projection, tol)?;
        seed = landed.map(|(_, uv)| uv);
        stations.push((t, seed));
    }
    // Between a station that lands and one that does not, bisected to
    // where the landing stops.
    let edge_of = |inside: (f64, Point2), outside: f64| -> OgeomResult<(f64, Point2)> {
        let (mut good, mut seed) = inside;
        let mut bad = outside;
        for _ in 0..60 {
            let mid = 0.5 * (good + bad);
            match land(
                surface,
                curve.point_at(mid, tol)?,
                Some(seed),
                projection,
                tol,
            )? {
                Some((_, uv)) => (good, seed) = (mid, uv),
                None => bad = mid,
            }
        }
        Ok((good, seed))
    };
    let mut runs = Vec::new();
    let mut k = 0;
    while k < stations.len() {
        let Some(uv) = stations[k].1 else {
            k += 1;
            continue;
        };
        let first = k;
        while k + 1 < stations.len() && stations[k + 1].1.is_some() {
            k += 1;
        }
        let (lo, seed) = if first == 0 {
            (range.0, uv)
        } else {
            edge_of((stations[first].0, uv), stations[first - 1].0)?
        };
        let hi = match (stations[k].1, stations.get(k + 1)) {
            (Some(last), Some(next)) => edge_of((stations[k].0, last), next.0)?.0,
            _ => range.1,
        };
        if hi - lo > (range.1 - range.0) * 1e-9 {
            runs.push(((lo, hi), seed));
        }
        k += 1;
    }
    Ok(runs)
}

/// The curve through the landings of `span` of the edge, fitted together
/// with its chart trace, the first landing sought from `start`.
fn projected(
    edge: &Shape,
    curve: &Curve,
    span: (f64, f64),
    start: Point2,
    surface: &SurfaceGeometry,
    projection: Projection,
    tol: Tolerances,
) -> OgeomResult<Segment> {
    const SAMPLES: usize = 96;
    let mut points = Vec::with_capacity(SAMPLES + 1);
    let mut chart: Vec<Point2> = Vec::with_capacity(SAMPLES + 1);
    let mut seed = Some(start);
    for k in 0..=SAMPLES {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a sample index, far below the mantissa"
        )]
        let t = (span.1 - span.0).mul_add(k as f64 / SAMPLES as f64, span.0);
        // A run's end is where the landing stops, found to rounding; it is
        // read a hair inside where it falls on the wrong side of that line.
        let inward = (span.1 - span.0) * if k == 0 { 1e-12 } else { -1e-12 };
        let tries = if k == 0 || k == SAMPLES { 16 } else { 0 };
        let mut landed = None;
        for step in 0..=tries {
            let at = curve.point_at(f64::from(step).mul_add(inward, t), tol)?;
            landed = land(surface, at, seed, projection, tol)?;
            if landed.is_some() {
                break;
            }
        }
        let Some((at, uv)) = landed else {
            ogeom_bail!(
                Construction,
                "the projection of an edge breaks off between two landings"
            );
        };
        points.push(at);
        chart.push(uv);
        seed = Some(uv);
    }
    unwrap(&mut chart, surface);
    let (fitted, trace, _) =
        ogeom_geom::fit::fit_points_joint(&points, &chart, &chart, 3, tol.confusion(), tol)?;
    let error = fitted.error;
    let fitted: Curve = Curve::BSpline(fitted.curve);
    let range = fitted.domain();
    segment(
        edge,
        fitted,
        PlanarCurve::BSpline(trace),
        range,
        error,
        surface,
        tol,
    )
}

/// Undo the folding a periodic chart applies: every step longer than half a
/// period is that period the other way.
fn unwrap(chart: &mut [Point2], surface: &SurfaceGeometry) {
    let (pu, pv) = periods(surface);
    for k in 1..chart.len() {
        let previous = chart[k - 1];
        let here = &mut chart[k];
        if pu > 0.0 {
            here.x -= ((here.x - previous.x) / pu).round() * pu;
        }
        if pv > 0.0 {
            here.y -= ((here.y - previous.y) / pv).round() * pv;
        }
    }
}

// --- where a path meets a face -----------------------------------------------

/// Where a path meets a face's boundary.
#[derive(Debug, Clone)]
struct Meeting {
    /// The path's parameter.
    t: f64,
    /// The path's chart point there.
    at: Point2,
    on: On,
}

#[derive(Debug, Clone)]
enum On {
    /// At a vertex of the boundary.
    Vertex(Shape),
    /// Inside a boundary edge (forward), a fraction of the way along it.
    Inside(Shape, f64),
}

/// One cut to make.
#[derive(Debug, Clone)]
enum Stretch {
    /// From one boundary meeting to the next along the path.
    Open(Meeting, Meeting),
    /// The whole closed path, inside the face.
    Loop,
}

/// A stretch with its path's parameter span.
struct Planned {
    stretch: Stretch,
    span: (f64, f64),
}

enum Found {
    Nothing,
    /// The path ends inside the face.
    Dangling,
    /// The closed path meets the boundary at this parameter: walk it from
    /// there.
    Rebase(f64),
    Cut(Planned),
}

/// The face's boundary in its chart.
struct Chart {
    rings: Vec<Vec<Occurrence>>,
    outlines: Vec<Vec<Point2>>,
    /// Chart resolution.
    eps: f64,
}

fn chart_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Chart> {
    let face_fwd = forward(face);
    let data = face_data(model, &face_fwd)?;
    let rings = rings(model, &face_fwd, data.surface, tol)?;
    let bounds = chart_bounds(&rings, tol);
    let span = (bounds.0.1 - bounds.0.0).max(bounds.1.1 - bounds.1.0);
    let outlines = rings
        .iter()
        .map(|ring| -> OgeomResult<Vec<Point2>> {
            let mut out = Vec::new();
            for o in ring {
                let mut line = o.polyline(64, tol)?;
                line.pop();
                out.extend(line);
            }
            Ok(out)
        })
        .collect::<OgeomResult<_>>()?;
    Ok(Chart {
        rings,
        outlines,
        eps: span * 1e-9 + tol.parametric(),
    })
}

/// The next stretch of `cutter` to cut in `face`, if any.
fn next_stretch(
    model: &Model,
    face: &Shape,
    cutter: &Cutter,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Found> {
    let chart = chart_of(model, face, tol)?;
    let meetings = meetings(model, &chart, cutter, surface, tol)?;
    let t_eps = cutter.resolution(tol);
    let (lo, hi) = cutter.range();

    // Whether the stretch `(a, b)` lies inside the face, refusing one that
    // runs along its boundary.
    let within = |a: f64, b: f64| -> OgeomResult<bool> {
        let mid = 0.5 * (a + b);
        let at = cutter.chart(mid, tol)?;
        let p = cutter.point(mid, tol)?;
        for o in chart.rings.iter().flatten() {
            if cutter.own.contains(&o.edge.node()) {
                continue;
            }
            let (_, near) = nearest_on(o, at, tol)?;
            let reach = boundary_reach(model, o, cutter, tol)?;
            if surface.point_at(near.x, near.y, tol)?.distance(p) <= reach {
                ogeom_bail!(
                    Construction,
                    "a curve to split along runs along the face's boundary"
                );
            }
        }
        Ok(inside(&chart.outlines, at))
    };

    if meetings.is_empty() {
        if cutter.is_consumed(lo, hi) || !within(lo, hi)? {
            return Ok(Found::Nothing);
        }
        if cutter.closed {
            return Ok(Found::Cut(Planned {
                stretch: Stretch::Loop,
                span: (lo, hi),
            }));
        }
        return Ok(Found::Dangling);
    }
    if cutter.closed {
        return Ok(Found::Rebase(meetings[0].t));
    }
    let (first, last) = (&meetings[0], &meetings[meetings.len() - 1]);
    for (a, b) in [(lo, first.t), (last.t, hi)] {
        if b - a > t_eps && !cutter.is_consumed(a, b) && within(a, b)? {
            return Ok(Found::Dangling);
        }
    }
    for pair in meetings.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        if b.t - a.t <= t_eps || cutter.is_consumed(a.t, b.t) {
            continue;
        }
        if within(a.t, b.t)? {
            return Ok(Found::Cut(Planned {
                stretch: Stretch::Open(a.clone(), b.clone()),
                span: (a.t, b.t),
            }));
        }
    }
    Ok(Found::Nothing)
}

/// How near a path must come to a boundary edge to be on it.
fn boundary_reach(
    model: &Model,
    o: &Occurrence,
    cutter: &Cutter,
    tol: Tolerances,
) -> OgeomResult<f64> {
    Ok(edge_data(model, &o.edge)?.tolerance.get() + cutter.tolerance() + tol.confusion())
}

/// The fraction along `o` nearest the chart point `p`, and its point.
fn nearest_on(o: &Occurrence, p: Point2, tol: Tolerances) -> OgeomResult<(f64, Point2)> {
    const N: usize = 64;
    let line = o.polyline(N, tol)?;
    let mut best = (0, f64::INFINITY);
    for (i, q) in line.iter().enumerate() {
        let d = q.distance(p);
        if d < best.1 {
            best = (i, d);
        }
    }
    #[expect(clippy::cast_precision_loss, reason = "a sample index")]
    let step = 1.0 / N as f64;
    #[expect(clippy::cast_precision_loss, reason = "a sample index")]
    let centre = best.0 as f64 * step;
    let (mut a, mut b) = ((centre - step).max(0.0), (centre + step).min(1.0));
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    let far = |s: f64| -> OgeomResult<f64> { Ok(o.at(s, tol)?.distance(p)) };
    for _ in 0..80 {
        let (x, y) = (b - g * (b - a), a + g * (b - a));
        if far(x)? <= far(y)? {
            b = y;
        } else {
            a = x;
        }
    }
    let s = 0.5 * (a + b);
    Ok((s, o.at(s, tol)?))
}

/// Every meeting of the path with the face's boundary, in order along the
/// path, one per point.
fn meetings(
    model: &Model,
    chart: &Chart,
    cutter: &Cutter,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Vec<Meeting>> {
    const ALONG: usize = 128;
    const ACROSS: usize = 64;
    let (lo, hi) = cutter.range();
    let samples = ALONG * cutter.path.len();
    let track: Vec<(f64, Point2)> = (0..=samples)
        .map(|i| {
            #[expect(clippy::cast_precision_loss, reason = "a sample index")]
            let t = (hi - lo).mul_add(i as f64 / samples as f64, lo);
            Ok((t, cutter.chart(t, tol)?))
        })
        .collect::<OgeomResult<_>>()?;
    let mut found: Vec<Meeting> = Vec::new();
    for o in chart.rings.iter().flatten() {
        if cutter.own.contains(&o.edge.node()) {
            continue;
        }
        let line = o.polyline(ACROSS, tol)?;
        for j in 1..line.len() {
            for i in 1..track.len() {
                let Some((x, y)) =
                    segments_cross((track[i - 1].1, track[i].1), (line[j - 1], line[j]))
                else {
                    continue;
                };
                let t = (track[i].0 - track[i - 1].0).mul_add(x, track[i - 1].0);
                #[expect(clippy::cast_precision_loss, reason = "a sample index")]
                let s = ((j - 1) as f64 + y) / ACROSS as f64;
                if let Some((t, s)) = polish(cutter, o, t, s, chart.eps, tol)? {
                    found.push(meeting(model, o, cutter, t, s, tol)?);
                }
            }
        }
        // An open path's ends, where they stand on this edge.
        if !cutter.closed {
            let reach = boundary_reach(model, o, cutter, tol)?;
            for t in [lo, hi] {
                let at = cutter.chart(t, tol)?;
                let (s, near) = nearest_on(o, at, tol)?;
                let p = cutter.point(t, tol)?;
                if surface.point_at(near.x, near.y, tol)?.distance(p) <= reach {
                    found.push(meeting(model, o, cutter, t, s, tol)?);
                }
            }
        }
    }
    found.sort_by(|a, b| a.t.total_cmp(&b.t));
    // One meeting per point: a crossing found on two edges at their shared
    // vertex is that vertex.
    let t_eps = cutter.resolution(tol);
    let mut out: Vec<Meeting> = Vec::new();
    for m in found {
        match out.last_mut() {
            Some(held) if (m.t - held.t).abs() <= t_eps * 1e3 => {
                if matches!(m.on, On::Vertex(_)) {
                    *held = m;
                }
            }
            _ => out.push(m),
        }
    }
    Ok(out)
}

/// Newton on `path(t) = o(s)` from a seed, in the chart.
fn polish(
    cutter: &Cutter,
    o: &Occurrence,
    mut t: f64,
    mut s: f64,
    eps: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    let (lo, hi) = cutter.range();
    let h = 1e-7;
    for _ in 0..40 {
        let f = cutter.chart(t, tol)? - o.at(s, tol)?;
        if f.magnitude() <= eps * 1e-3 {
            break;
        }
        let dc = cutter.chart_d1(t, tol)?;
        let (s0, s1) = ((s - h).max(0.0), (s + h).min(1.0));
        let dq = (o.at(s1, tol)? - o.at(s0, tol)?) / (s1 - s0);
        // Solve [dc, -dq] (dt, ds) = -f.
        let det = dc.x.mul_add(-dq.y, dq.x * dc.y);
        if det.abs() <= f64::MIN_POSITIVE {
            return Ok(None);
        }
        let dt = (-f.x).mul_add(-dq.y, dq.x * (-f.y)) / det;
        let ds = dc.x.mul_add(-f.y, f.x * dc.y) / det;
        t = (t + dt).clamp(lo, hi);
        s = (s + ds).clamp(0.0, 1.0);
    }
    let f = cutter.chart(t, tol)? - o.at(s, tol)?;
    Ok((f.magnitude() <= eps).then_some((t, s)))
}

/// The meeting a fraction `s` along `o`, at the path's `t`: a vertex where
/// it stands on one.
fn meeting(
    model: &Model,
    o: &Occurrence,
    cutter: &Cutter,
    t: f64,
    s: f64,
    tol: Tolerances,
) -> OgeomResult<Meeting> {
    let at = cutter.chart(t, tol)?;
    let edge = forward(&o.edge);
    let forward_s = if o.reversed() { 1.0 - s } else { s };
    // Rounded so two traversals of one seam name the same point.
    let forward_s = (forward_s * 1e12).round() / 1e12;
    let p = cutter.point(t, tol)?;
    if let Some((first, last)) = edge_vertices(model, &edge)? {
        for v in [first, last] {
            let held = model
                .node(&v)
                .and_then(|n| n.data().as_vertex())
                .map_or(0.0, |d| d.tolerance.get());
            if vertex_point(model, &v)?.distance(p) <= held + cutter.tolerance() + tol.confusion() {
                return Ok(Meeting {
                    t,
                    at,
                    on: On::Vertex(v),
                });
            }
        }
    }
    Ok(Meeting {
        t,
        at,
        on: On::Inside(edge, forward_s),
    })
}

/// Where two chart segments cross, as fractions along each, their ends
/// included.
fn segments_cross(a: (Point2, Point2), b: (Point2, Point2)) -> Option<(f64, f64)> {
    let da = a.1 - a.0;
    let db = b.1 - b.0;
    let cross = da.x.mul_add(db.y, -(da.y * db.x));
    if cross.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let w = b.0 - a.0;
    let x = w.x.mul_add(db.y, -(w.y * db.x)) / cross;
    let y = w.x.mul_add(da.y, -(w.y * da.x)) / cross;
    let slack = 1e-9;
    ((-slack..=1.0 + slack).contains(&x) && (-slack..=1.0 + slack).contains(&y))
        .then_some((x.clamp(0.0, 1.0), y.clamp(0.0, 1.0)))
}

// --- one cut -----------------------------------------------------------------

/// What one cut made.
struct Made {
    reshape: Reshape,
    /// The edges cut along, each with the edge it came from.
    edges: Vec<(Shape, Shape)>,
    /// The boundary edges cut, and their pieces.
    splits: Vec<(TShapeId, Vec<Shape>)>,
}

/// A boundary edge's pieces, each with the fractions of the edge it spans.
type Parts = Vec<(Shape, f64, f64)>;

/// One directed piece of a ring.
#[derive(Debug, Clone)]
struct Item {
    edge: Shape,
    from: TShapeId,
    polyline: Vec<Point2>,
}

/// Cut `face` along one stretch of `cutter`.
#[allow(clippy::too_many_lines)]
fn cut(
    model: &mut Model,
    face: &Shape,
    cutter: &Cutter,
    planned: &Planned,
    tol: Tolerances,
) -> OgeomResult<Made> {
    let face_fwd = forward(face);
    let data = face_data(model, &face_fwd)?;
    let chart = chart_of(model, face, tol)?;
    let mut reshape = Reshape::new();
    let mut splits: Vec<(TShapeId, Vec<Shape>)> = Vec::new();

    // The boundary edges the stretch ends inside, cut there.
    let ends: Vec<Meeting> = match &planned.stretch {
        Stretch::Open(a, b) => vec![a.clone(), b.clone()],
        Stretch::Loop => Vec::new(),
    };
    let mut fractions: Vec<(Shape, Vec<f64>)> = Vec::new();
    for m in &ends {
        if let On::Inside(edge, s) = &m.on {
            match fractions.iter_mut().find(|(e, _)| e.node() == edge.node()) {
                Some((_, list)) => list.push(*s),
                None => fractions.push((edge.clone(), vec![*s])),
            }
        }
    }
    let mut pieces_of: Vec<(TShapeId, Parts)> = Vec::new();
    let mut vertex_at: Vec<((TShapeId, u64), Shape)> = Vec::new();
    for (edge, mut list) in fractions {
        list.sort_by(f64::total_cmp);
        list.dedup_by(|b, a| (*a - *b).abs() <= 1e-12);
        let (pieces, vertices) = split_edge(model, &edge, &list, tol)?;
        for (s, v) in list.iter().zip(&vertices) {
            vertex_at.push(((edge.node(), s.to_bits()), v.clone()));
        }
        let mut bounds = vec![0.0];
        bounds.extend(&list);
        bounds.push(1.0);
        pieces_of.push((
            edge.node(),
            pieces
                .iter()
                .zip(bounds.windows(2))
                .map(|(p, w)| (p.clone(), w[0], w[1]))
                .collect(),
        ));
        reshape.split(&edge, pieces.clone());
        splits.push((edge.node(), pieces));
    }
    let vertex_of = |m: &Meeting| -> OgeomResult<Shape> {
        match &m.on {
            On::Vertex(v) => Ok(v.clone()),
            On::Inside(edge, s) => vertex_at
                .iter()
                .find(|(key, _)| *key == (edge.node(), s.to_bits()))
                .map(|(_, v)| v.clone())
                .ok_or_else(|| ogeom_err!(Construction, "a cut vertex is lost")),
        }
    };

    // The rings as directed pieces.
    let mut items: Vec<Vec<Item>> = Vec::new();
    for ring in &chart.rings {
        let mut out = Vec::new();
        for o in ring {
            let edge = forward(&o.edge);
            let parts = pieces_of
                .iter()
                .find(|(node, _)| *node == edge.node())
                .map_or_else(|| vec![(edge.clone(), 0.0, 1.0)], |(_, p)| p.clone());
            let ordered: Parts = if o.reversed() {
                parts
                    .into_iter()
                    .rev()
                    .map(|(p, a, b)| (p.reversed(), 1.0 - b, 1.0 - a))
                    .collect()
            } else {
                parts
            };
            for (piece, a, b) in ordered {
                let polyline: Vec<Point2> = (0..=16)
                    .map(|i| o.at((b - a).mul_add(f64::from(i) / 16.0, a), tol))
                    .collect::<OgeomResult<_>>()?;
                let Some((from, _)) = edge_vertices(model, &piece)? else {
                    ogeom_bail!(Construction, "a boundary piece has no vertices");
                };
                out.push(Item {
                    edge: piece,
                    from: from.node(),
                    polyline,
                });
            }
        }
        items.push(out);
    }
    let areas: Vec<f64> = items
        .iter()
        .map(|ring| signed_area(&outline(ring)))
        .collect();
    let outer = (0..areas.len())
        .max_by(|a, b| areas[*a].abs().total_cmp(&areas[*b].abs()))
        .unwrap_or(0);
    let outer_sign = areas.get(outer).copied().unwrap_or(1.0).signum();

    // The new edges, one per segment the stretch runs through.
    let (va, vb) = match &planned.stretch {
        Stretch::Open(a, b) => (vertex_of(a)?, vertex_of(b)?),
        Stretch::Loop => {
            let p = cutter.point(planned.span.0, tol)?;
            let v = model.add_vertex(VertexData::new(p));
            (v.clone(), v)
        }
    };
    let (along, edges) = cut_edges(model, cutter, planned.span, &va, &vb, data.surface, tol)?;
    let back: Vec<Item> = along
        .iter()
        .rev()
        .map(|item| -> OgeomResult<Item> {
            let Some((from, _)) = edge_vertices(model, &item.edge.reversed())? else {
                ogeom_bail!(Construction, "a cut edge has no vertices");
            };
            Ok(Item {
                edge: item.edge.reversed(),
                from: from.node(),
                polyline: item.polyline.iter().rev().copied().collect(),
            })
        })
        .collect::<OgeomResult<_>>()?;

    // The pieces' rings: (outer, holes).
    let mut pieces: Vec<(Vec<Item>, Vec<Vec<Item>>)> = Vec::new();
    let probe = |ring: &[Item]| ring[0].polyline[ring[0].polyline.len() / 2];
    match &planned.stretch {
        Stretch::Open(a, b) => {
            let (ra, ia) = junction(&items, va.node(), a.at)?;
            let (rb, ib) = junction(&items, vb.node(), b.at)?;
            if ra != rb {
                ogeom_bail!(
                    Construction,
                    "a curve running from one boundary loop of the face to \
                     another does not divide it"
                );
            }
            let ring = &items[ra];
            let n = ring.len();
            // The ring from one junction up to the other; where both ends
            // are one junction, nothing one way and the whole ring the other.
            let arc = |from: usize, to: usize, whole: bool| -> Vec<Item> {
                let mut out = Vec::new();
                let mut k = from;
                while out.len() < n && (k != to || (whole && out.is_empty())) {
                    out.push(ring[k].clone());
                    k = (k + 1) % n;
                }
                out
            };
            let mut first = arc(ia, ib, false);
            first.extend(back.iter().cloned());
            let mut second = arc(ib, ia, true);
            second.extend(along.iter().cloned());
            let (s1, s2) = (
                signed_area(&outline(&first)) * outer_sign,
                signed_area(&outline(&second)) * outer_sign,
            );
            if ra == outer {
                if s1 <= 0.0 || s2 <= 0.0 {
                    ogeom_bail!(Construction, "the cut does not divide the face in two");
                }
                let (mut h1, mut h2) = (Vec::new(), Vec::new());
                let line1 = outline(&first);
                for (k, hole) in items.iter().enumerate() {
                    if k == ra {
                        continue;
                    }
                    if inside(std::slice::from_ref(&line1), probe(hole)) {
                        h1.push(hole.clone());
                    } else {
                        h2.push(hole.clone());
                    }
                }
                pieces.push((first, h1));
                pieces.push((second, h2));
            } else {
                // A cut from a hole back to it: one loop is a new piece,
                // the other the hole the rest keeps.
                let (fresh, hole) = match (s1 > 0.0, s2 > 0.0) {
                    (true, false) => (first, second),
                    (false, true) => (second, first),
                    _ => ogeom_bail!(Construction, "the cut does not divide the face in two"),
                };
                let line = outline(&fresh);
                let (mut held, mut rest) = (Vec::new(), vec![hole]);
                for (k, ring) in items.iter().enumerate() {
                    if k == ra || k == outer {
                        continue;
                    }
                    if inside(std::slice::from_ref(&line), probe(ring)) {
                        held.push(ring.clone());
                    } else {
                        rest.push(ring.clone());
                    }
                }
                pieces.push((items[outer].clone(), rest));
                pieces.push((fresh, held));
            }
        }
        Stretch::Loop => {
            let (disk, hole) = if signed_area(&outline(&along)) * outer_sign > 0.0 {
                (along, back)
            } else {
                (back, along)
            };
            let line = outline(&disk);
            let (mut held, mut rest) = (Vec::new(), vec![hole]);
            for (k, ring) in items.iter().enumerate() {
                if k == outer {
                    continue;
                }
                if inside(std::slice::from_ref(&line), probe(ring)) {
                    held.push(ring.clone());
                } else {
                    rest.push(ring.clone());
                }
            }
            pieces.push((items[outer].clone(), rest));
            pieces.push((disk, held));
        }
    }

    let mut faces = Vec::with_capacity(pieces.len());
    for (outer_ring, holes) in pieces {
        let mut wires = vec![model.add_wire(&edges_of(&outer_ring))?];
        for hole in holes {
            wires.push(model.add_wire(&edges_of(&hole))?);
        }
        let mut fresh = data.clone();
        fresh.triangulation = None;
        fresh.natural_restriction = false;
        faces.push(model.add_face(fresh, &wires)?);
    }
    reshape.split(&face_fwd, faces);
    Ok(Made {
        reshape,
        edges,
        splits,
    })
}

fn edges_of(ring: &[Item]) -> Vec<Shape> {
    ring.iter().map(|i| i.edge.clone()).collect()
}

fn outline(ring: &[Item]) -> Vec<Point2> {
    ring.iter()
        .flat_map(|i| i.polyline.iter().copied())
        .collect()
}

/// The ring and the index of the piece leaving the vertex `v` nearest the
/// chart point `at`.
fn junction(items: &[Vec<Item>], v: TShapeId, at: Point2) -> OgeomResult<(usize, usize)> {
    let mut best: Option<(usize, usize, f64)> = None;
    for (r, ring) in items.iter().enumerate() {
        for (k, item) in ring.iter().enumerate() {
            if item.from != v {
                continue;
            }
            let d = item.polyline[0].distance(at);
            if best.is_none_or(|(.., held)| d < held) {
                best = Some((r, k, d));
            }
        }
    }
    best.map(|(r, k, _)| (r, k)).ok_or_else(|| {
        ogeom_err!(
            Construction,
            "a cut ends at a vertex the face does not have"
        )
    })
}

/// A stretch's new edges as directed pieces, and each edge with the edge
/// it came from.
type Walked = (Vec<Item>, Vec<(Shape, Shape)>);

/// The edges along `cutter` over `span`, from `from` to `to`, one per
/// segment the span runs through, each with its trace on the surface.
fn cut_edges(
    model: &mut Model,
    cutter: &Cutter,
    span: (f64, f64),
    from: &Shape,
    to: &Shape,
    surface: SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Walked> {
    let eps = cutter.resolution(tol);
    let mut stops = vec![span.0];
    let mut joint = span.0.floor() + 1.0;
    while joint < span.1 - eps {
        if joint > span.0 + eps {
            stops.push(joint);
        }
        joint += 1.0;
    }
    stops.push(span.1);
    // A vertex at each joint, where the two segments meeting there end.
    let mut vertices = vec![from.clone()];
    for &t in &stops[1..stops.len() - 1] {
        let before = cutter.locate(t - 0.5).0.point(1.0, tol)?;
        let after = cutter.point(t, tol)?;
        let at = Point::from_vector((before.to_vector() + after.to_vector()) * 0.5);
        vertices.push(model.add_vertex(VertexData::new(at)));
    }
    vertices.push(to.clone());

    let mut items = Vec::with_capacity(stops.len() - 1);
    let mut edges = Vec::with_capacity(stops.len() - 1);
    for (k, w) in stops.windows(2).enumerate() {
        let (segment, f0) = cutter.locate(0.5 * (w[0] + w[1]));
        let base = 0.5f64.mul_add(-(w[1] - w[0]), f0);
        let (f0, f1) = (base, base + (w[1] - w[0]));
        let (t0, t1) = (segment.local(f0), segment.local(f1));
        let range = (t0.min(t1), t0.max(t1));
        let (start, end) = if t0 <= t1 {
            (&vertices[k], &vertices[k + 1])
        } else {
            (&vertices[k + 1], &vertices[k])
        };
        let curve_id = model.geometry_mut().add_curve(segment.curve.clone());
        let mut data = EdgeData::on_curve(curve_id, Location::identity(), range);
        let trace = model.geometry_mut().add_pcurve(segment.pcurve.clone());
        data.add(EdgeRepr::PCurve {
            curve: trace,
            surface,
            location: Location::identity(),
            range,
        });
        data.assert_same_parameter(true);
        // The vertices were placed by the edges they cut or as the joints'
        // middles; the new edge's ends reach them within the gap.
        let reach = vertex_point(model, start)?
            .distance(segment.curve.point_at(range.0, tol)?)
            .max(vertex_point(model, end)?.distance(segment.curve.point_at(range.1, tol)?));
        let edge = model.add_edge(data, &[start.clone(), end.clone()])?;
        model.widen(
            &edge,
            Tolerance::new(reach.max(segment.tolerance) + tol.confusion())?,
        )?;
        let walked = if t0 <= t1 {
            edge.clone()
        } else {
            edge.reversed()
        };
        let polyline: Vec<Point2> = (0..=32)
            .map(|i| segment.chart((f1 - f0).mul_add(f64::from(i) / 32.0, f0), tol))
            .collect::<OgeomResult<_>>()?;
        items.push(Item {
            edge: walked,
            from: vertices[k].node(),
            polyline,
        });
        edges.push((segment.source.clone(), edge));
    }
    Ok((items, edges))
}
