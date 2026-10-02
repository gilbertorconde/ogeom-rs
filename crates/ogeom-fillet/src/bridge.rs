//! Bridging a gap: a blend surface between edges of two faces, and a blend
//! curve between the ends of two edges.
//!
//! Both are Hermite constructions across the gap. Each side contributes its
//! position, a crossing tangent pointing off its face (or along its curve)
//! for G1, and the second derivative that matches its curvature for G2; the
//! polynomial across is the Bézier of the least degree holding all of them.
//!
//! The surface runs along the two edges in `u` and across the gap in `v`,
//! with `v = 0` on the first edge and `v = 1` on the second. Its rows next
//! to each edge are that edge's own B-spline control polygon, made
//! compatible with the other's (the same degree and knots), so the edges
//! lie on the surface exactly and are shared as they are, gaining a trim on
//! the blend. The crossing derivatives along each edge are B-splines on the
//! same knots, interpolated at the Greville abscissae. On a plane every
//! interpolated tangent lies in the plane, so the join is tangent (and for
//! G2 curvature continuous) exactly; on a curved face the join is measured
//! along the edge and the knots refined until it is within budget.

use ogeom_algo::{
    Built, History, attach_pcurve, edge_vertices, make_edge_between, make_face_on, make_wire,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSpline2d, BSplineCurve, BSplineSurface, Continuity, Curve, Curve2d as _, Curve3d as _, Line2d,
    PlanarCurve, Surface as _, SurfaceGeometry,
};
use ogeom_math::{
    Axis2, ControlGrid, Direction, Direction2, KnotVector, Point, Point2, Vector, Vector2, Weighted,
};
use ogeom_topo::{
    EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, explore_unique,
};

/// One end of an edge, in the direction the edge is traversed: the start of
/// a reversed edge is the end of its curve's range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
    /// Where the edge begins.
    Start,
    /// Where the edge finishes.
    End,
}

/// The largest angle, in radians, between the blend's normal and a face's
/// at a sampled station on their shared edge, for a G1 or G2 join.
const TANGENCY_BUDGET: f64 = 1e-5;

/// The largest difference between the blend's and a face's normal
/// curvatures square to the edge at a station, for a G2 join, as a fraction
/// of the inverse of the gap's width there.
const CURVATURE_BUDGET: f64 = 1e-4;

/// How many times the knots along the edges are halved before a join that
/// stays outside its budget is refused.
const REFINEMENTS: usize = 8;

/// Stations measured in each knot span along the edges.
const STATIONS_PER_SPAN: usize = 4;

/// Bridge the gap between two edges with a blend face.
///
/// `first` and `second` are each an edge and the face it bounds. The blend
/// is a B-spline face whose boundary is the two edges, shared as they are,
/// and two new side edges joining their ends, paired so the sides do not
/// cross. It meets each face with the continuity asked of that side:
///
/// - [`Continuity::C0`]: along the edge only;
/// - [`Continuity::G1`]: tangent to the face, its normal within `1e-5`
///   radians of the face's at every sampled station;
/// - [`Continuity::G2`]: tangent as above, and bending alike square to the
///   edge, the two normal curvatures there apart by at most `1e-4` divided
///   by the gap's width.
///
/// The crossing tangent's length on each side is the gap's width, the
/// distance between the two edges at the same `u`, so a blend across a
/// wider stretch of gap swings out further. Joins to planar faces are exact; on a
/// curved face the knots along the edges are refined until the measured join
/// is within budget.
///
/// Each edge gains a trim on the blend's surface. History: each edge
/// generates the blend face and the two side edges.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by name,
/// where:
///
/// - a continuity is [`Continuity::C1`], [`Continuity::C2`] or
///   [`Continuity::CInfinity`]: a blend's joins are geometric, to G2 at most;
/// - a shape is not an edge or a face, an edge is not on the face given with
///   it, or both are the same edge;
/// - an edge or a face is placed;
/// - an edge is degenerate, closed, a seam of its face, or has no trim on its
///   face's surface;
/// - an edge's curve has no exact B-spline form on its own parameter (a
///   circle or an ellipse, whose B-spline form runs on another), so it
///   cannot be shared with the blend exactly;
/// - the face's side of an edge cannot be told from its boundary, or the face
///   has no normal at a point of the edge;
/// - the edges meet at an end, where the blend would have no side;
/// - on a curved face, the join stays outside its budget after refinement.
pub fn make_blend_surface(
    model: &mut Model,
    first: (&Shape, &Shape),
    second: (&Shape, &Shape),
    continuity: (Continuity, Continuity),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let orders = (order_of(continuity.0)?, order_of(continuity.1)?);
    if first.0.is_same(second.0) {
        ogeom_bail!(
            Construction,
            "a blend bridges two different edges; both are the same one"
        );
    }
    let mut near = read_side(model, first.0, first.1, orders.0, tol)?;
    let mut far = read_side(model, second.0, second.1, orders.1, tol)?;

    // The second edge runs whichever way pairs its ends with the first's
    // without the sides crossing.
    let (a0, a1) = (near.spline.start(tol)?, near.spline.end(tol)?);
    let (b0, b1) = (far.spline.start(tol)?, far.spline.end(tol)?);
    far.reversed = a0.distance(b1) + a1.distance(b0) < a0.distance(b0) + a1.distance(b1);
    if far.reversed {
        let (knots, control) =
            ogeom_math::bspline::reverse(far.spline.knots(), far.spline.control_points());
        far.spline = BSplineCurve::rational(knots, control)?;
    }
    let (b0, b1) = (far.spline.start(tol)?, far.spline.end(tol)?);
    if a0.distance(b0).min(a1.distance(b1)) <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the edges meet at an end; a blend between them would have no side there"
        );
    }

    compatible(&mut near.spline, &mut far.spline, tol)?;
    let mut rounds = 0;
    let patch = loop {
        let patch = hermite(&near, &far, tol)?;
        let near_worst = measure(&patch, &near, &far, 0.0, tol)?;
        let far_worst = measure(&patch, &near, &far, 1.0, tol)?;
        let angle = near_worst.0.max(far_worst.0);
        let bend = near_worst.1.max(far_worst.1);
        if angle <= TANGENCY_BUDGET && bend <= CURVATURE_BUDGET {
            break patch;
        }
        if rounds == REFINEMENTS {
            ogeom_bail!(
                Construction,
                "the blend meets its faces {angle} radians off tangent and {bend} off their \
                 curvature (as a fraction of the gap's inverse width) after refining its knots \
                 {REFINEMENTS} times, outside the budget of {TANGENCY_BUDGET} and \
                 {CURVATURE_BUDGET}"
            );
        }
        refine(&mut near.spline, &mut far.spline, tol)?;
        rounds += 1;
    };

    let sides_of = [
        Curve::BSpline(patch.iso_u_curve(0.0, tol)?),
        Curve::BSpline(patch.iso_u_curve(1.0, tol)?),
    ];
    let surface_id = model
        .geometry_mut()
        .add_surface(SurfaceGeometry::BSpline(patch));
    let ends_of = |side: &Side| {
        if side.reversed {
            (side.vertices.1.clone(), side.vertices.0.clone())
        } else {
            side.vertices.clone()
        }
    };
    let (near_ends, far_ends) = (ends_of(&near), ends_of(&far));
    let [at_start, at_end] = sides_of;
    let side_start = side_edge(
        model,
        at_start,
        0.0,
        (&near_ends.0, &far_ends.0),
        surface_id,
        tol,
    )?;
    let side_end = side_edge(
        model,
        at_end,
        1.0,
        (&near_ends.1, &far_ends.1),
        surface_id,
        tol,
    )?;
    for (side, v) in [(&near, 0.0), (&far, 1.0)] {
        let (s0, s1) = if side.reversed {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };
        let knots = KnotVector::new(
            vec![side.range.0, side.range.0, side.range.1, side.range.1],
            1,
        )?;
        let trim = BSpline2d::new(knots, vec![Point2::new(s0, v), Point2::new(s1, v)], tol)?;
        attach_pcurve(
            model,
            &side.edge,
            trim.into(),
            surface_id,
            Location::identity(),
            side.range,
        )?;
    }

    let ring = [
        near.edge.clone(),
        side_end.clone(),
        if far.reversed {
            far.edge.clone()
        } else {
            far.edge.reversed()
        },
        side_start.reversed(),
    ];
    let wire = make_wire(model, &ring, tol)?.shape;
    let face = make_face_on(model, surface_id, &[wire], tol)?.shape;

    let mut history = History::new();
    for edge in [first.0, second.0] {
        history.generate(edge, face.clone());
        history.generate(edge, side_start.clone());
        history.generate(edge, side_end.clone());
    }
    Ok(Built::new(face, history))
}

/// Bridge the gap between the ends of two edges with a blend edge.
///
/// The blend is a Bézier curve from the end of `a` to the end of `b`,
/// joining their vertices. On each side it continues the edge with the
/// continuity asked: [`Continuity::C0`] meets it, [`Continuity::G1`] also
/// leaves along its tangent, [`Continuity::G2`] also with its curvature
/// vector. The tangent's length on each side is the distance between the
/// two ends.
///
/// History: each edge generates the blend edge.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by name,
/// where a continuity is [`Continuity::C1`], [`Continuity::C2`] or
/// [`Continuity::CInfinity`]; a shape is not an edge, is placed or
/// degenerate, or has no curve; an edge stands still at the end asked, so
/// it has no tangent there; or the two ends are one point.
pub fn make_blend_curve(
    model: &mut Model,
    a: (&Shape, End),
    b: (&Shape, End),
    continuity: (Continuity, Continuity),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let orders = (order_of(continuity.0)?, order_of(continuity.1)?);
    let (from, near) = joint(model, a.0, a.1, orders.0, tol)?;
    let (to, far) = joint(model, b.0, b.1, orders.1, tol)?;
    let width = near.point.distance(far.point);
    if width <= tol.confusion() || from.is_same(&to) {
        ogeom_bail!(
            Construction,
            "the two ends are one point; there is no gap to bridge"
        );
    }
    let control = near.column(&far, width);
    let degree = control.len() - 1;
    let mut knots = vec![0.0; degree + 1];
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    let curve = BSplineCurve::new(KnotVector::new(knots, degree)?, control, tol)?;
    let edge = make_edge_between(model, Curve::BSpline(curve), (0.0, 1.0), &from, &to, tol)?.shape;
    let mut history = History::new();
    history.generate(a.0, edge.clone());
    history.generate(b.0, edge.clone());
    Ok(Built::new(edge, history))
}

/// How many derivatives across a join `continuity` asks to match.
fn order_of(continuity: Continuity) -> OgeomResult<usize> {
    match continuity {
        Continuity::C0 => Ok(0),
        Continuity::G1 => Ok(1),
        Continuity::G2 => Ok(2),
        Continuity::C1 | Continuity::C2 => ogeom_bail!(
            Construction,
            "a blend's join is geometric: parametric continuity would tie it to the two \
             sides' parameterizations; ask G1 or G2"
        ),
        Continuity::CInfinity => ogeom_bail!(
            Construction,
            "a blend is built to G2 at most; a join smooth to every order is not offered"
        ),
    }
}

/// One side of the blend: an edge, its face, and the edge's curve as a
/// B-spline on `[0, 1]` that runs on the edge's own parameter, scaled.
struct Side {
    /// The edge, forward.
    edge: Shape,
    /// Its vertices at the start and the end of its curve's range.
    vertices: (Shape, Shape),
    /// The face's surface.
    surface: SurfaceGeometry,
    /// The edge's trim on it, with the trim's range.
    trim: (PlanarCurve, (f64, f64)),
    /// The edge's curve and range.
    curve: Curve,
    range: (f64, f64),
    /// The curve as a B-spline on `[0, 1]`, reversed where `reversed`.
    spline: BSplineCurve,
    /// Whether the blend's `u` runs against the edge.
    reversed: bool,
    /// The sign that turns the face's normal crossed with the edge's
    /// tangent into the direction off the face.
    off: f64,
    /// How many derivatives across the join are matched.
    order: usize,
}

/// Where the blend's crossing tangent starts on one side, at one station.
struct Station {
    /// The face's unit normal.
    normal: Vector,
    /// The unit direction square to the edge, in the face's tangent plane,
    /// pointing off the face.
    off: Vector,
    /// The face's normal curvature along `off` times its normal: the second
    /// derivative a unit-speed crossing curve has across the edge.
    bend: Vector,
    /// The edge's tangent.
    along: Vector,
}

/// The side `edge` of `face` gives the blend, checked.
fn read_side(
    model: &Model,
    edge: &Shape,
    face: &Shape,
    order: usize,
    tol: Tolerances,
) -> OgeomResult<Side> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(
            Construction,
            "a blend surface bridges edges; got a {:?}",
            model.kind_of(edge)?
        );
    }
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "each edge is given with the face it bounds; got a {:?}",
            model.kind_of(face)?
        );
    }
    let Some(NodeData::Face(face_data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let surface_id = face_data.surface;
    if !face.location().is_identity()
        || !face_data.location.is_identity()
        || !edge.location().is_identity()
    {
        ogeom_bail!(
            Construction,
            "a placed edge or face is not blended; bake its placement into its geometry first"
        );
    }
    if !explore_unique(model, face, ShapeType::Edge)?
        .iter()
        .any(|e| e.is_same(edge))
    {
        ogeom_bail!(Construction, "the edge is not on the face given with it");
    }
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let Some(NodeData::Edge(data)) = model.node(edge).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    if data.degenerate {
        ogeom_bail!(
            Construction,
            "a degenerate edge has no extent to blend from"
        );
    }
    let Some(EdgeRepr::Curve3d {
        curve,
        location,
        range,
    }) = data.curve3d()
    else {
        ogeom_bail!(Construction, "the edge has no curve in space");
    };
    if !location.is_identity() {
        ogeom_bail!(
            Construction,
            "a placed edge is not blended; bake its placement into its geometry first"
        );
    }
    let range = *range;
    let Some(curve) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let trim = match data.pcurve_for(surface_id, edge.location()) {
        Some(EdgeRepr::PCurve { curve, range, .. }) => {
            let Some(pcurve) = model.geometry().pcurve(*curve).cloned() else {
                ogeom_bail!(Dangling, "pcurve is not in this model");
            };
            (pcurve, *range)
        }
        Some(EdgeRepr::Seam { .. }) => ogeom_bail!(
            Construction,
            "the edge is a seam of its face; there is no one side to blend from"
        ),
        _ => ogeom_bail!(
            Construction,
            "the edge has no trim on its face's surface, so the face's side of it is unknown"
        ),
    };
    let edge = edge.oriented(Orientation::Forward);
    let Some((start, end)) = edge_vertices(model, &edge)? else {
        ogeom_bail!(Construction, "the edge has no vertices");
    };
    if start.is_same(&end)
        || curve
            .point_at(range.0, tol)?
            .distance(curve.point_at(range.1, tol)?)
            <= tol.confusion()
    {
        ogeom_bail!(
            Construction,
            "the edge is closed; a blend's sides run from the ends of its edges"
        );
    }

    let spline = match curve.to_bspline_over(range, tol) {
        Ok(spline) if !spline.is_periodic() => spline,
        _ => ogeom_bail!(
            Construction,
            "the edge's curve has no exact B-spline form to share with the blend"
        ),
    };
    for k in 0..=8 {
        let f = f64::from(k) / 8.0;
        let on_edge = curve.point_at((range.1 - range.0).mul_add(f, range.0), tol)?;
        if spline.point_at(f, tol)?.distance(on_edge) > tol.confusion() {
            ogeom_bail!(
                Construction,
                "the edge's curve runs on a parameter its B-spline form does not share (as a \
                 circle's or an ellipse's does), so the blend cannot share the edge exactly"
            );
        }
    }

    let mut side = Side {
        edge,
        vertices: (start, end),
        surface,
        trim,
        curve,
        range,
        spline,
        reversed: false,
        off: 1.0,
        order,
    };
    side.off = off_sign(model, face, &side, tol)?;
    Ok(side)
}

/// Which way the face's normal crossed with the edge's tangent points
/// relative to the face: `1` where it points off it, `-1` where onto it.
/// Read at the edge's middle by stepping either way across its trim.
fn off_sign(model: &Model, face: &Shape, side: &Side, tol: Tolerances) -> OgeomResult<f64> {
    let rings = ogeom_mesh::face_boundary(model, face, ogeom_mesh::Deflection::default(), tol)?;
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in rings.iter().flatten() {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let size = lo.distance(hi);
    let (pcurve, prange) = &side.trim;
    let p = f64::midpoint(prange.0, prange.1);
    let uv = pcurve.point_at(p, tol)?;
    let heading = pcurve.d1_at(p, tol)?;
    let across = Vector2::new(-heading.y, heading.x);
    let Ok(across) = across.normalized(tol) else {
        ogeom_bail!(
            Construction,
            "the edge's trim stands still at its middle, so the face's side of it is unknown"
        );
    };
    let mut outward = None;
    for step in [1e-2, 1e-3, 1e-4] {
        let h = across * (step * size);
        let (left, right) = (
            ogeom_mesh::inside_boundary(&rings, uv + h),
            ogeom_mesh::inside_boundary(&rings, uv - h),
        );
        if left != right {
            outward = Some(if left { -1.0 } else { 1.0 });
            break;
        }
    }
    let Some(outward) = outward else {
        ogeom_bail!(
            Construction,
            "the face's side of the edge cannot be told from its boundary"
        );
    };
    let (du, dv) = side.surface.d1_at(uv.x, uv.y, tol)?;
    let off_chart = (du * across.x + dv * across.y) * outward;
    let t = (side.range.0 + side.range.1) * 0.5;
    let normal = face_normal(&side.surface, uv, tol)?;
    let crossed = normal.cross(side.curve.d1_at(t, tol)?);
    let agree = crossed.dot(off_chart);
    if agree.abs() <= f64::EPSILON * crossed.magnitude() * off_chart.magnitude() {
        ogeom_bail!(
            Construction,
            "the face's side of the edge cannot be told from its boundary"
        );
    }
    Ok(agree.signum())
}

/// The unit normal of a face's surface at a point of the edge.
fn face_normal(surface: &SurfaceGeometry, uv: Point2, tol: Tolerances) -> OgeomResult<Vector> {
    match surface.normal_at(uv.x, uv.y, tol) {
        Ok(n) => Ok(n.vector()),
        Err(_) => ogeom_bail!(
            Construction,
            "the face has no normal at a point of the edge, so there is no tangent plane to \
             meet"
        ),
    }
}

impl Side {
    /// The edge's parameter and its trim's at the blend's `u = s`.
    fn parameters(&self, s: f64) -> (f64, f64) {
        let f = if self.reversed { 1.0 - s } else { s };
        let (t0, t1) = self.range;
        let (p0, p1) = self.trim.1;
        ((t1 - t0).mul_add(f, t0), (p1 - p0).mul_add(f, p0))
    }

    /// The face's frame at the blend's `u = s`.
    fn station(&self, s: f64, tol: Tolerances) -> OgeomResult<Station> {
        let (t, p) = self.parameters(s);
        let uv = self.trim.0.point_at(p, tol)?;
        let normal = face_normal(&self.surface, uv, tol)?;
        let along = self.curve.d1_at(t, tol)?;
        let Ok(off) = (normal.cross(along) * self.off).normalized(tol) else {
            ogeom_bail!(
                Construction,
                "the edge stands still at a point, so it has no direction to blend square to"
            );
        };
        let bend = if self.order >= 2 {
            let Ok(curvature) = self.surface.curvature_at(uv.x, uv.y, tol) else {
                ogeom_bail!(
                    Construction,
                    "the face has no curvature at a point of the edge to match"
                );
            };
            let Some(k) = curvature.normal_curvature(off) else {
                ogeom_bail!(
                    Construction,
                    "the face has no curvature at a point of the edge to match"
                );
            };
            curvature.normal.vector() * k
        } else {
            Vector::ZERO
        };
        Ok(Station {
            normal,
            off,
            bend,
            along,
        })
    }
}

/// The derivatives one side of a Hermite column holds.
struct Hermite {
    point: Point,
    /// The unit direction the column leaves `point` in.
    out: Vector,
    /// The second derivative of a unit-speed crossing curve at `point`.
    bend: Vector,
    order: usize,
}

impl Hermite {
    /// The Bézier control points from `self` to `far`, of the least degree
    /// matching both sides' derivatives, the crossing tangent `width` long.
    fn column(&self, far: &Self, width: f64) -> Vec<Point> {
        let n = self.order + far.order + 1;
        #[expect(clippy::cast_precision_loss, reason = "a degree of at most five")]
        let (n1, n2) = (n as f64, (n * n.saturating_sub(1)) as f64);
        let mut rows = vec![Point::ORIGIN; n + 1];
        rows[0] = self.point;
        rows[n] = far.point;
        for (side, at, step) in [(self, 0usize, 1isize), (far, n, -1isize)] {
            let index = |k: isize| at.wrapping_add_signed(step * k);
            if side.order >= 1 {
                rows[index(1)] = side.point + side.out * (width / n1);
            }
            if side.order >= 2 {
                rows[index(2)] = rows[index(1)]
                    + (rows[index(1)] - side.point)
                    + side.bend * (width * width / n2);
            }
        }
        rows
    }
}

/// Give two B-splines on `[0, 1]` the same degree and knots without moving
/// either. Knots closer than the parametric tolerance are one knot.
fn compatible(a: &mut BSplineCurve, b: &mut BSplineCurve, tol: Tolerances) -> OgeomResult<()> {
    while a.degree() < b.degree() {
        *a = a.elevated(tol)?;
    }
    while b.degree() < a.degree() {
        *b = b.elevated(tol)?;
    }
    let merge = |into: &BSplineCurve, from: &BSplineCurve| -> OgeomResult<BSplineCurve> {
        let mut out = into.clone();
        for (k, m) in from.knots().distinct() {
            let have = out.knots().distinct();
            let (value, held) = have
                .iter()
                .find(|(h, _)| (h - k).abs() <= tol.parametric())
                .map_or((k, 0), |(h, c)| (*h, *c));
            if held < m {
                out = out.with_knot_inserted(value, m - held, tol)?;
            }
        }
        Ok(out)
    };
    *b = merge(b, a)?;
    *a = merge(a, b)?;
    let (ka, kb) = (a.knots().knots(), b.knots().knots());
    if ka.len() != kb.len()
        || ka
            .iter()
            .zip(kb)
            .any(|(x, y)| (x - y).abs() > tol.parametric())
    {
        ogeom_bail!(
            Invariant,
            "the two edges' B-spline forms did not come to share their knots"
        );
    }
    Ok(())
}

/// Halve every knot span of both B-splines.
fn refine(a: &mut BSplineCurve, b: &mut BSplineCurve, tol: Tolerances) -> OgeomResult<()> {
    let distinct = a.knots().distinct();
    for pair in distinct.windows(2) {
        let middle = f64::midpoint(pair[0].0, pair[1].0);
        *a = a.with_knot_inserted(middle, 1, tol)?;
        *b = b.with_knot_inserted(middle, 1, tol)?;
    }
    Ok(())
}

/// The Greville abscissae of a knot vector: where each control point
/// stands along the curve.
fn greville(knots: &KnotVector) -> Vec<f64> {
    let p = knots.degree();
    let k = knots.knots();
    #[expect(
        clippy::cast_precision_loss,
        reason = "a degree, far below the mantissa"
    )]
    let scale = 1.0 / p as f64;
    (0..knots.control_point_count())
        .map(|i| k[i + 1..=i + p].iter().sum::<f64>() * scale)
        .collect()
}

/// The blend patch: the two splines' control polygons as its outer rows,
/// each column a Bézier across the gap.
///
/// Next to each edge, the differences between the rows are the control
/// vectors of the crossing derivatives, which share the edge's weights, so
/// each derivative along the edge is a B-spline in its own right. They are
/// interpolated at the Greville abscissae: the first to the direction off
/// the face, the gap's width long, and the second to the face's normal
/// curvature along it.
fn hermite(near: &Side, far: &Side, tol: Tolerances) -> OgeomResult<BSplineSurface> {
    let n = near.order + far.order + 1;
    let knots = near.spline.knots().clone();
    let abscissae = greville(&knots);
    let (ca, cb) = (near.spline.control_points(), far.spline.control_points());
    let count = ca.len();
    let mut widths = Vec::with_capacity(count);
    for s in &abscissae {
        widths.push(
            near.spline
                .point_at(*s, tol)?
                .distance(far.spline.point_at(*s, tol)?),
        );
    }
    #[expect(clippy::cast_precision_loss, reason = "a degree of at most five")]
    let (n1, n2) = (n as f64, (n * (n - 1)) as f64);
    // Per side, the control vectors of the first and second derivatives
    // across, each divided by the factor the Bézier form puts on it.
    let crossing = |side: &Side, control: &[Weighted<Point>]| -> OgeomResult<[Vec<Vector>; 2]> {
        if side.order == 0 {
            return Ok([vec![Vector::ZERO; count], vec![Vector::ZERO; count]]);
        }
        let mut first = Vec::with_capacity(count);
        let mut second = Vec::with_capacity(count);
        for (s, width) in abscissae.iter().zip(&widths) {
            let station = side.station(*s, tol)?;
            first.push(station.off * (width / n1));
            second.push(if side.order >= 2 {
                station.bend * (width * width / n2)
            } else {
                Vector::ZERO
            });
        }
        let weights: Vec<f64> = control.iter().map(|c| c.weight).collect();
        let matrix = collocation(&knots, &weights, &abscissae, tol)?;
        Ok([solve(&matrix, &first)?, solve(&matrix, &second)?])
    };
    let [near_first, near_second] = crossing(near, ca)?;
    let [far_first, far_second] = crossing(far, cb)?;
    let mut grid = Vec::with_capacity(count * (n + 1));
    for i in 0..count {
        let (pa, pb) = (ca[i].point(), cb[i].point());
        let mut column = vec![Point::ORIGIN; n + 1];
        column[0] = pa;
        column[n] = pb;
        if near.order >= 1 {
            column[1] = pa + near_first[i];
        }
        if near.order >= 2 {
            column[2] = pa + near_first[i] * 2.0 + near_second[i];
        }
        if far.order >= 1 {
            column[n - 1] = pb + far_first[i];
        }
        if far.order >= 2 {
            column[n - 2] = pb + far_first[i] * 2.0 + far_second[i];
        }
        for (j, point) in column.into_iter().enumerate() {
            // The rows a side's derivatives are read from share its
            // weights, so the weight function is flat across the edge.
            let weight = if j <= near.order {
                ca[i].weight
            } else {
                cb[i].weight
            };
            grid.push(Weighted::new(point, weight, tol)?);
        }
    }
    let mut v_knots = vec![0.0; n + 1];
    v_knots.extend(std::iter::repeat_n(1.0, n + 1));
    BSplineSurface::rational(
        knots,
        KnotVector::new(v_knots, n)?,
        ControlGrid::new(grid, count, n + 1)?,
    )
}

/// The rational basis functions with `weights` over `knots`, row `k` read
/// at `at[k]`.
fn collocation(
    knots: &KnotVector,
    weights: &[f64],
    at: &[f64],
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<f64>>> {
    let p = knots.degree();
    let mut rows = Vec::with_capacity(at.len());
    for u in at {
        let span = knots.span(*u, tol)?;
        let basis = knots.basis(span, *u);
        let mut row = vec![0.0; weights.len()];
        let mut total = 0.0;
        for (r, value) in basis.iter().enumerate() {
            let i = span - p + r;
            row[i] = value * weights[i];
            total += row[i];
        }
        for value in &mut row {
            *value /= total;
        }
        rows.push(row);
    }
    Ok(rows)
}

/// `matrix · x = rhs` for vectors `x`, by elimination with partial
/// pivoting.
fn solve(matrix: &[Vec<f64>], rhs: &[Vector]) -> OgeomResult<Vec<Vector>> {
    let size = rhs.len();
    let mut a: Vec<Vec<f64>> = matrix.to_vec();
    let mut b: Vec<Vector> = rhs.to_vec();
    for col in 0..size {
        let pivot = (col..size)
            .max_by(|x, y| a[*x][col].abs().total_cmp(&a[*y][col].abs()))
            .unwrap_or(col);
        if a[pivot][col].abs() <= f64::EPSILON {
            ogeom_bail!(
                Numeric,
                "the crossing derivatives cannot be interpolated along the edge"
            );
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (above, below) = a.split_at_mut(col + 1);
        let leading = &above[col];
        for (offset, target) in below.iter_mut().enumerate() {
            let factor = target[col] / leading[col];
            if factor == 0.0 {
                continue;
            }
            for (t, l) in target[col..].iter_mut().zip(&leading[col..]) {
                *t -= factor * l;
            }
            let pivot_value = b[col];
            b[col + 1 + offset] -= pivot_value * factor;
        }
    }
    let mut x = vec![Vector::ZERO; size];
    for row in (0..size).rev() {
        let mut sum = b[row];
        for k in row + 1..size {
            sum -= x[k] * a[row][k];
        }
        x[row] = sum * (1.0 / a[row][row]);
    }
    Ok(x)
}

/// The worst join along the edge at `v`: the largest angle between the two
/// normals, and the largest difference of normal curvatures square to the
/// edge times the gap's width there. Each is zero where its side does not
/// ask for it.
fn measure(
    patch: &BSplineSurface,
    near: &Side,
    far: &Side,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    let side = if v == 0.0 { near } else { far };
    if side.order == 0 {
        return Ok((0.0, 0.0));
    }
    let (mut angle, mut bend) = (0.0f64, 0.0f64);
    let knots = side.spline.knots().distinct();
    for pair in knots.windows(2) {
        for k in 0..=STATIONS_PER_SPAN {
            #[expect(clippy::cast_precision_loss, reason = "a station index")]
            let f = k as f64 / STATIONS_PER_SPAN as f64;
            let s = (pair[1].0 - pair[0].0).mul_add(f, pair[0].0);
            let station = side.station(s, tol)?;
            let Ok(blend_normal) = patch.normal_at(s, v, tol) else {
                ogeom_bail!(
                    Construction,
                    "the blend folds onto the edge, so it has no tangent plane there"
                );
            };
            let blend_normal = blend_normal.vector();
            let along = blend_normal.dot(station.normal);
            angle = angle.max(along.abs().min(1.0).acos());
            if side.order >= 2 {
                let width = near
                    .spline
                    .point_at(s, tol)?
                    .distance(far.spline.point_at(s, tol)?);
                let Ok(curvature) = patch.curvature_at(s, v, tol) else {
                    ogeom_bail!(
                        Construction,
                        "the blend has no curvature at a point of the edge to match"
                    );
                };
                let kf = station.bend.dot(station.normal);
                let Some(kb) =
                    curvature.normal_curvature(curvature.normal.vector().cross(station.along))
                else {
                    ogeom_bail!(
                        Construction,
                        "the blend has no curvature at a point of the edge to match"
                    );
                };
                let kb = if curvature.normal.vector().dot(station.normal) < 0.0 {
                    -kb
                } else {
                    kb
                };
                bend = bend.max((kb - kf).abs() * width);
            }
        }
    }
    Ok((angle, bend))
}

/// A new side edge along the patch's parameter line `u = at`, from the
/// first edge's end to the second's, trimmed on the patch.
fn side_edge(
    model: &mut Model,
    curve: Curve,
    at: f64,
    (from, to): (&Shape, &Shape),
    surface: ogeom_topo::SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let edge = make_edge_between(model, curve, (0.0, 1.0), from, to, tol)?.shape;
    let trim = Line2d::over(
        Axis2 {
            location: Point2::new(at, 0.0),
            direction: Direction2::new(Vector2::new(0.0, 1.0), tol)?,
        },
        0.0,
        1.0,
    )?;
    attach_pcurve(
        model,
        &edge,
        trim.into(),
        surface,
        Location::identity(),
        (0.0, 1.0),
    )?;
    if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
        data.assert_same_parameter(true);
    }
    Ok(edge)
}

/// The end `end` of `edge`, read for a blend curve: its vertex, and the
/// derivatives the curve matches there to `order`, leaving along the
/// direction that continues the edge past its end.
fn joint(
    model: &Model,
    edge: &Shape,
    end: End,
    order: usize,
    tol: Tolerances,
) -> OgeomResult<(Shape, Hermite)> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(
            Construction,
            "a blend curve bridges the ends of edges; got a {:?}",
            model.kind_of(edge)?
        );
    }
    if !edge.location().is_identity() {
        ogeom_bail!(
            Construction,
            "a placed edge is not blended; bake its placement into its geometry first"
        );
    }
    let Some(NodeData::Edge(data)) = model.node(edge).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    if data.degenerate {
        ogeom_bail!(Construction, "a degenerate edge has no tangent to continue");
    }
    let Some(EdgeRepr::Curve3d {
        curve,
        location,
        range,
    }) = data.curve3d()
    else {
        ogeom_bail!(Construction, "the edge has no curve in space");
    };
    if !location.is_identity() {
        ogeom_bail!(
            Construction,
            "a placed edge is not blended; bake its placement into its geometry first"
        );
    }
    let Some(curve) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let Some((start, finish)) = edge_vertices(model, edge)? else {
        ogeom_bail!(Construction, "the edge has no vertices");
    };
    let backwards = edge.orientation() == Orientation::Reversed;
    // At the curve's range end the edge continues forward along its
    // parameter; at the range start it continues backward.
    let (vertex, at_range_end) = match end {
        End::Start => (start, backwards),
        End::End => (finish, !backwards),
    };
    let t = if at_range_end { range.1 } else { range.0 };
    let d = curve.derivatives_at(t, 2, tol)?;
    let speed = d[1].magnitude();
    let Ok(tangent) = Direction::new(d[1], tol) else {
        ogeom_bail!(
            Construction,
            "the edge stands still at that end, so it has no tangent to continue"
        );
    };
    let tangent = tangent.vector();
    let bend = (d[2] - tangent * d[2].dot(tangent)) * (1.0 / (speed * speed));
    Ok((
        vertex,
        Hermite {
            point: curve.point_at(t, tol)?,
            out: if at_range_end { tangent } else { -tangent },
            bend,
            order,
        },
    ))
}
