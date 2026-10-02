//! The N-sided filling: one face over a hole bounded by any number of
//! edges, meeting each side's face at G0, G1 or G2.

use ogeom_algo::{Built, History, attach_pcurve, edge_vertices, make_face, make_wire};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSpline2d, BSplineSurface, Continuity, Curve, Curve2d as _, Curve3d as _, CurveKind,
    PlanarCurve, Surface as _, SurfaceCurvature, SurfaceGeometry, Trig2d,
};
use ogeom_math::{Direction, Point, Point2, Vector, Vector2, Weighted};
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Orientation, Shape, ShapeType, explore,
};

use crate::fill_patch::{Condition, DEGREE, PlaneFrame, fit_height};

/// One side of an N-sided filling.
#[derive(Debug, Clone)]
pub struct FillBoundary {
    /// The boundary edge. The filling's face is bounded by this edge node
    /// itself, not by a copy of it.
    pub edge: Shape,
    /// The face the edge belongs to, which the filling meets across it.
    /// Required for [`Continuity::G1`] and [`Continuity::G2`]; on a
    /// [`Continuity::C0`] side it is measured against and helps decide which
    /// way the filling faces.
    pub support: Option<Shape>,
    /// How the filling meets `support`: [`Continuity::C0`] (position),
    /// [`Continuity::G1`] (tangent plane) or [`Continuity::G2`] (tangent
    /// plane and normal curvature).
    pub continuity: Continuity,
}

/// What a filling achieved along one side, measured at stations spread
/// over the edge.
#[derive(Debug, Clone, PartialEq)]
pub struct FillSide {
    /// The side's edge.
    pub edge: Shape,
    /// The largest distance, in model units, between the edge's curve and
    /// the filling's surface read through the edge's pcurve on it.
    pub gap: f64,
    /// The largest angle, in radians, between the filling's normal and the
    /// support's; `None` on a side with no support.
    pub angle: Option<f64>,
    /// The largest difference, in inverse model units, between the
    /// filling's and the support's normal curvatures square to the edge,
    /// signed against one shared normal; `None` on a side with no support
    /// or where no station gave a curvature on both surfaces.
    pub curvature: Option<f64>,
    /// How many stations were measured.
    pub stations: usize,
}

/// A filling and what it achieved.
#[derive(Debug, Clone)]
pub struct Filled {
    /// The face; every boundary edge and constraint generates it.
    pub built: Built,
    /// One report per side, in the order the sides were given.
    pub sides: Vec<FillSide>,
    /// The largest distance, along the filling plane's normal, from a
    /// constraint point or a sampled point of a constraint curve to the
    /// surface; zero with no constraints.
    pub constraint_gap: f64,
}

/// Samples per side for the loop's outline.
const OUTLINE_SAMPLES: usize = 64;
/// The margin round the hole, as a share of its larger extent.
const MARGIN: f64 = 0.05;
/// The control counts across the larger extent, round by round.
const NETS: [usize; 6] = [8, 12, 18, 27, 40, 60];
/// The bending energy's weight against the conditions: small enough that
/// the conditions win wherever they reach, leaving the energy to settle
/// the controls they do not (the hole's interior, the margin round it).
const SMOOTHING: f64 = 1e-12;
/// The least cosine between a support's normal and the plane's normal:
/// about 84 degrees.
const MIN_LIFT: f64 = 0.1;
/// A constraint point's weight, as the share of the hole's size a boundary
/// sample of that length would carry.
const POINT_SHARE: f64 = 0.05;
/// Samples per constraint curve.
const PER_CURVE: usize = 32;

/// Fill the hole a loop of edges bounds with one face meeting each side's
/// support at the continuity asked.
///
/// The sides may come in any order and either direction; they must chain
/// into one simple closed loop. `constraints` are vertices and edges inside
/// the hole that the surface passes through. The surface is a cubic
/// B-spline height patch over the plane the loop spans, fitted by least
/// squares to the sides' positions, to the tangent planes of G1 and G2
/// sides' supports and to the normal curvatures of G2 sides' supports, with
/// a thin-plate bending energy settling the rest; the control net is
/// refined until every side meets `tolerance` or the refinement runs out.
/// The face is trimmed by the given edges themselves, each given a pcurve
/// on the patch, and faces the way the supports say: across each edge it
/// runs opposite to its support's use of the edge, so sewing it to the
/// supports gives a consistently oriented shell. With no supports it faces
/// the side the loop turns counter-clockwise about, walked from the first
/// side in that edge's own direction.
///
/// `tolerance` is the target for every measured deviation: a distance in
/// model units for each side's gap and each constraint, an angle in radians
/// for a G1 or G2 side's tangency, and a curvature difference in inverse
/// model units for a G2 side. An edge the surface stands further from than
/// the edge's own tolerance has that tolerance, and its vertices', widened
/// to the measured gap.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where:
///
/// - `tolerance` is not positive and finite, or there are no sides;
/// - a side is not an edge, has no 3D curve or no vertices, or is placed;
/// - the same edge is given twice;
/// - a side asks [`Continuity::C1`], [`Continuity::C2`] or
///   [`Continuity::CInfinity`] (parametric continuity between two
///   surfaces' charts), or G1 or G2 with no support;
/// - a support is not a face, is placed, does not hold its edge, or the
///   edge does not lie on its surface within `tolerance`;
/// - the sides do not chain into one closed loop;
/// - the loop encloses no area, or crosses itself seen along the normal of
///   the plane it spans, so the hole is not a height field over that plane;
/// - a G1 or G2 side's support stands within about 6 degrees of square to
///   that plane;
/// - a constraint is neither a vertex nor an edge, or lies outside the hole
///   seen along the plane's normal.
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the finest
/// control net still misses `tolerance`, naming the side and the deviation.
pub fn make_filling_n(
    model: &mut Model,
    boundary: &[FillBoundary],
    constraints: &[Shape],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Filled> {
    if !(tolerance.is_finite() && tolerance > 0.0) {
        ogeom_bail!(
            Construction,
            "a filling's tolerance is positive and finite; got {tolerance}"
        );
    }
    if boundary.is_empty() {
        ogeom_bail!(Construction, "a filling needs at least one boundary edge");
    }
    let mut sides = Vec::with_capacity(boundary.len());
    for (i, entry) in boundary.iter().enumerate() {
        sides.push(read_side(model, i, entry, tolerance, tol)?);
    }
    for i in 0..sides.len() {
        for j in (i + 1)..sides.len() {
            if sides[i].edge.node() == sides[j].edge.node() {
                ogeom_bail!(Construction, "side {j} is side {i}'s edge again");
            }
        }
    }
    let mut order = chain(model, &mut sides, tol)?;
    face_the_supports(model, &mut sides, &mut order)?;

    // The plane the loop spans, and the loop seen along its normal.
    let outline = loop_points(&sides, &order, tol)?;
    let frame = frame_of(&outline, tol)?;
    let chart_outline: Vec<Point2> = outline.iter().map(|p| frame.chart(*p)).collect();
    if crosses_itself(&chart_outline) {
        ogeom_bail!(
            Construction,
            "the boundary loop crosses itself seen along the normal of the \
             plane it spans; the hole is not a height field over that plane"
        );
    }
    let interior = constraint_points(model, constraints, tol)?;
    for (k, (p, _)) in interior.iter().enumerate() {
        if !inside(&chart_outline, frame.chart(*p)) {
            ogeom_bail!(
                Construction,
                "a point of constraint {k} at {p:?} lies outside the hole seen \
                 along the normal of the plane the boundary spans"
            );
        }
    }

    // The rectangle the patch covers: the hole and a margin round it.
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for q in &chart_outline {
        lo = Point2::new(lo.x.min(q.x), lo.y.min(q.y));
        hi = Point2::new(hi.x.max(q.x), hi.y.max(q.y));
    }
    let margin = MARGIN * (hi.x - lo.x).max(hi.y - lo.y);
    let domain = (
        (lo.x - margin, hi.x + margin),
        (lo.y - margin, hi.y + margin),
    );
    let (du, dv) = (domain.0.1 - domain.0.0, domain.1.1 - domain.1.0);
    let size = du.max(dv);

    let mut traces = Vec::with_capacity(sides.len());
    let mut lengths = Vec::with_capacity(sides.len());
    for side in &sides {
        traces.push(trace(side, &frame, tolerance, tol)?);
        lengths.push(side_length(side, tol)?);
    }

    let mut last_miss = String::new();
    for base in NETS {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            reason = "a control count of a few dozen, from a positive ratio"
        )]
        let count = |extent: f64| -> usize {
            ((base as f64 * extent / size).round() as usize).max(DEGREE + 2)
        };
        let controls = (count(du), count(dv));
        let samples = 4 * controls.0.max(controls.1) + 8;

        let mut conditions = Vec::new();
        for (side, length) in sides.iter().zip(&lengths) {
            side_conditions(
                side,
                *length / size,
                size,
                samples,
                &frame,
                &mut conditions,
                tol,
            )?;
        }
        for (p, share) in &interior {
            conditions.push(Condition {
                at: frame.chart(*p),
                order: (0, 0),
                target: frame.height(*p),
                weight: share.sqrt() / size,
            });
        }
        let surface = fit_height(&frame, domain, controls, &conditions, SMOOTHING, tol)?;

        let mut reports = Vec::with_capacity(sides.len());
        for (side, pcurve) in sides.iter().zip(&traces) {
            reports.push(measure_side(side, pcurve, &surface, 3 * samples + 1, tol)?);
        }
        let mut constraint_gap = 0.0f64;
        for (p, _) in &interior {
            let q = frame.chart(*p);
            constraint_gap = constraint_gap.max(surface.point_at(q.x, q.y, tol)?.distance(*p));
        }

        match first_miss(&sides, &reports, constraint_gap, tolerance) {
            Some(miss) => {
                last_miss = format!("at {}x{} controls, {miss}", controls.0, controls.1);
            }
            None => {
                let face = build(model, &sides, &order, &traces, &reports, surface, tol)?;
                let mut history = History::new();
                for side in &sides {
                    history.generate(&side.edge, face.clone());
                }
                for constraint in constraints {
                    history.generate(constraint, face.clone());
                }
                return Ok(Filled {
                    built: Built::new(face, history),
                    sides: reports,
                    constraint_gap,
                });
            }
        }
    }
    ogeom_bail!(
        NotDone,
        "the filling misses its tolerance of {tolerance} on the finest net: {last_miss}"
    )
}

/// The surface a side's support lies on and the side's pcurve there.
struct Support {
    face: Shape,
    surface: SurfaceGeometry,
    pcurve: PlanarCurve,
    prange: (f64, f64),
}

/// One side, read and checked.
struct Side {
    entry: usize,
    /// The edge node, forward.
    edge: Shape,
    curve: Curve,
    range: (f64, f64),
    edge_tolerance: f64,
    /// 0, 1 or 2: positions, tangent planes, normal curvatures.
    order: usize,
    support: Option<Support>,
    /// Whether the loop runs against the edge.
    reversed: bool,
}

impl Side {
    fn parameter(&self, f: f64) -> f64 {
        (self.range.1 - self.range.0).mul_add(f, self.range.0)
    }

    /// The support's chart position and unit normal at the edge's
    /// parameter `t`.
    fn support_at(&self, t: f64, tol: Tolerances) -> OgeomResult<Option<(Point2, Vector)>> {
        let Some(support) = &self.support else {
            return Ok(None);
        };
        let f = (t - self.range.0) / (self.range.1 - self.range.0);
        let pt = (support.prange.1 - support.prange.0).mul_add(f, support.prange.0);
        let uv = support.pcurve.point_at(pt, tol)?;
        let normal = support.surface.normal_at(uv.x, uv.y, tol)?.vector();
        Ok(Some((uv, normal)))
    }
}

fn read_side(
    model: &Model,
    i: usize,
    entry: &FillBoundary,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Side> {
    let edge = &entry.edge;
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "side {i} is not an edge");
    }
    if !edge.location().is_identity() {
        ogeom_bail!(
            Construction,
            "side {i}'s edge is placed; bake its placement into its geometry first"
        );
    }
    let order = match entry.continuity {
        Continuity::C0 => 0,
        Continuity::G1 => 1,
        Continuity::G2 => 2,
        other => ogeom_bail!(
            Construction,
            "side {i} asks {other:?}: parametric continuity between two \
             surfaces' charts is not what a filling meets; ask for G1 or G2"
        ),
    };
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "side {i}'s edge holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "side {i}'s edge has no 3D curve");
    };
    let Some(curve) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "side {i}'s curve is not in this model");
    };
    let range = *range;
    let edge_tolerance = data.tolerance.get();
    let edge = edge.oriented(Orientation::Forward);

    let support = match &entry.support {
        None if order > 0 => ogeom_bail!(
            Construction,
            "side {i} asks {:?} but names no support face to meet",
            entry.continuity
        ),
        None => None,
        Some(face) => Some(read_support(
            model,
            i,
            face,
            &edge,
            &curve,
            range,
            tolerance.max(edge_tolerance),
            tol,
        )?),
    };
    Ok(Side {
        entry: i,
        edge,
        curve,
        range,
        edge_tolerance,
        order,
        support,
        reversed: false,
    })
}

#[expect(
    clippy::too_many_arguments,
    reason = "the side's edge, curve and range are read once by the caller"
)]
fn read_support(
    model: &Model,
    i: usize,
    face: &Shape,
    edge: &Shape,
    curve: &Curve,
    range: (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Support> {
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "side {i}'s support is not a face");
    }
    let Some(NodeData::Face(face_data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "side {i}'s support is not in this model");
    };
    if !face.location().is_identity() || !face_data.location.is_identity() {
        ogeom_bail!(
            Construction,
            "side {i}'s support is placed; bake its placement into its geometry first"
        );
    }
    let holds = explore(model, face, Filter::OfType(ShapeType::Edge))?
        .iter()
        .any(|e| e.node() == edge.node());
    if !holds {
        ogeom_bail!(
            Construction,
            "side {i}'s support face does not hold the side's edge"
        );
    }
    let surface_id = face_data.surface;
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "side {i}'s support surface is not in this model");
    };
    let Some(edge_data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "side {i}'s edge holds no edge data");
    };
    let stored = match edge_data.pcurve_for(surface_id, edge.location()) {
        Some(
            EdgeRepr::PCurve { curve, range, .. }
            | EdgeRepr::Seam {
                forward: curve,
                range,
                ..
            },
        ) => model
            .geometry()
            .pcurve(*curve)
            .cloned()
            .map(|c| (c, *range)),
        _ => None,
    };
    let (pcurve, prange) = if let Some(found) = stored {
        found
    } else {
        let (fitted, _, _, _, _) =
            ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, range, &surface, tol)?;
        (fitted, range)
    };
    // The edge must lie on the surface it is said to bound.
    let mut worst = 0.0f64;
    for k in 0..=16 {
        let f = f64::from(k) / 16.0;
        let t = (range.1 - range.0).mul_add(f, range.0);
        let pt = (prange.1 - prange.0).mul_add(f, prange.0);
        let uv = pcurve.point_at(pt, tol)?;
        let on = surface.point_at(uv.x, uv.y, tol)?;
        worst = worst.max(on.distance(curve.point_at(t, tol)?));
    }
    if worst > reach {
        ogeom_bail!(
            Construction,
            "side {i}'s edge stands {worst} off its support face, past {reach}"
        );
    }
    Ok(Support {
        face: face.clone(),
        surface,
        pcurve,
        prange,
    })
}

/// Chain the sides into one loop: the order they are walked in, with each
/// side's `reversed` set to the direction it is walked.
fn chain(model: &Model, sides: &mut [Side], tol: Tolerances) -> OgeomResult<Vec<usize>> {
    let mut ends = Vec::with_capacity(sides.len());
    for side in sides.iter() {
        let Some(pair) = edge_vertices(model, &side.edge)? else {
            ogeom_bail!(
                Construction,
                "side {} has no vertices, so it cannot be shown to join the loop",
                side.entry
            );
        };
        ends.push(pair);
    }
    let meets = |a: &Shape, b: &Shape| -> OgeomResult<bool> {
        Ok(a.is_same(b) || model.same_position(a, b, tol)?)
    };
    let n = sides.len();
    let mut used = vec![false; n];
    used[0] = true;
    let mut order = vec![0];
    let first = ends[0].0.clone();
    let mut cursor = ends[0].1.clone();
    while order.len() < n {
        let last = sides[order[order.len() - 1]].entry;
        let mut found: Vec<(usize, bool)> = Vec::new();
        for j in 0..n {
            if used[j] {
                continue;
            }
            if meets(&ends[j].0, &cursor)? {
                found.push((j, false));
            } else if meets(&ends[j].1, &cursor)? {
                found.push((j, true));
            }
        }
        let (j, reversed) = match found.as_slice() {
            [one] => *one,
            [] => ogeom_bail!(
                Construction,
                "the boundary does not close: no side continues where side {last} ends"
            ),
            _ => ogeom_bail!(
                Construction,
                "{} sides continue where side {last} ends; a filling's boundary \
                 is one simple loop",
                found.len()
            ),
        };
        used[j] = true;
        sides[j].reversed = reversed;
        cursor = if reversed {
            ends[j].0.clone()
        } else {
            ends[j].1.clone()
        };
        order.push(j);
    }
    if !meets(&cursor, &first)? {
        ogeom_bail!(
            Construction,
            "the boundary does not close: the last side ends away from where the first begins"
        );
    }
    Ok(order)
}

/// Turn the loop to run against its supports' uses of the shared edges,
/// the majority deciding where they disagree.
fn face_the_supports(model: &Model, sides: &mut [Side], order: &mut [usize]) -> OgeomResult<()> {
    let (mut agree, mut disagree) = (0usize, 0usize);
    for side in sides.iter() {
        let Some(support) = &side.support else {
            continue;
        };
        let uses: Vec<Orientation> =
            explore(model, &support.face, Filter::OfType(ShapeType::Edge))?
                .iter()
                .filter(|e| e.node() == side.edge.node())
                .map(Shape::orientation)
                .collect();
        let Some(&first) = uses.first() else {
            continue;
        };
        if uses.iter().any(|o| *o != first)
            || !matches!(first, Orientation::Forward | Orientation::Reversed)
        {
            continue;
        }
        if (first == Orientation::Forward) == !side.reversed {
            disagree += 1;
        } else {
            agree += 1;
        }
    }
    if disagree > agree {
        order.reverse();
        for side in sides.iter_mut() {
            side.reversed = !side.reversed;
        }
    }
    Ok(())
}

/// Points round the loop in its walking order.
fn loop_points(sides: &[Side], order: &[usize], tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let mut out = Vec::with_capacity(order.len() * OUTLINE_SAMPLES);
    for &i in order {
        let side = &sides[i];
        for k in 0..OUTLINE_SAMPLES {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a sample index, far below the mantissa"
            )]
            let mut f = k as f64 / OUTLINE_SAMPLES as f64;
            if side.reversed {
                f = 1.0 - f;
            }
            out.push(side.curve.point_at(side.parameter(f), tol)?);
        }
    }
    Ok(out)
}

/// The plane a closed loop spans: its vector area's normal, which the loop
/// turns counter-clockwise about, through its centroid, with `e1` along the
/// loop's widest spread in the plane.
fn frame_of(outline: &[Point], tol: Tolerances) -> OgeomResult<PlaneFrame> {
    let Ok(origin) = Point::centroid(outline) else {
        ogeom_bail!(Construction, "the boundary loop has no points");
    };
    let mut area = Vector::ZERO;
    let mut reach = 0.0f64;
    for (k, p) in outline.iter().enumerate() {
        let q = outline[(k + 1) % outline.len()];
        area += (*p - origin).cross(q - origin) * 0.5;
        reach = reach.max(p.distance(origin));
    }
    let magnitude = area.magnitude();
    if magnitude <= 1e-9 * reach * reach || magnitude <= tol.confusion() * tol.confusion() {
        ogeom_bail!(
            Construction,
            "the boundary loop encloses no area seen from any plane"
        );
    }
    let n = area * (1.0 / magnitude);
    let b1 = Direction::new(n, tol)?.any_perpendicular().vector();
    let b2 = n.cross(b1);
    let (mut sxx, mut syy, mut sxy) = (0.0f64, 0.0f64, 0.0f64);
    for p in outline {
        let d = *p - origin;
        let (x, y) = (d.dot(b1), d.dot(b2));
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
    }
    let theta = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let e1 = b1 * theta.cos() + b2 * theta.sin();
    let e2 = n.cross(e1);
    Ok(PlaneFrame { origin, e1, e2, n })
}

/// Whether a closed polygon crosses itself: any two segments that are not
/// neighbours crossing.
fn crosses_itself(polygon: &[Point2]) -> bool {
    let n = polygon.len();
    let orient = |a: Point2, b: Point2, c: Point2| (b - a).cross(c - a);
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (c, d) = (polygon[j], polygon[(j + 1) % n]);
            if orient(a, b, c) * orient(a, b, d) < 0.0 && orient(c, d, a) * orient(c, d, b) < 0.0 {
                return true;
            }
        }
    }
    false
}

/// Whether a point lies inside a closed polygon, by its winding number.
fn inside(polygon: &[Point2], q: Point2) -> bool {
    let n = polygon.len();
    let mut winding = 0i32;
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        let side = (b - a).cross(q - a);
        if a.y <= q.y {
            if b.y > q.y && side > 0.0 {
                winding += 1;
            }
        } else if b.y <= q.y && side < 0.0 {
            winding -= 1;
        }
    }
    winding != 0
}

/// The points the constraints put inside the hole, each with its share of
/// the fit's weight.
fn constraint_points(
    model: &Model,
    constraints: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<(Point, f64)>> {
    let mut out = Vec::new();
    for (k, shape) in constraints.iter().enumerate() {
        let placement = shape.transform(model.datums())?;
        match model.kind_of(shape)? {
            ShapeType::Vertex => {
                let Some(data) = model.node(shape).and_then(|n| n.data().as_vertex()) else {
                    ogeom_bail!(Construction, "constraint {k} holds no vertex data");
                };
                out.push((placement.apply(data.point), POINT_SHARE));
            }
            ShapeType::Edge => {
                let Some(data) = model.node(shape).and_then(|n| n.data().as_edge()) else {
                    ogeom_bail!(Construction, "constraint {k} holds no edge data");
                };
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    ogeom_bail!(Construction, "constraint {k} has no 3D curve");
                };
                let Some(curve) = model.geometry().curve(*curve) else {
                    ogeom_bail!(Dangling, "constraint {k}'s curve is not in this model");
                };
                for i in 0..=PER_CURVE {
                    #[expect(
                        clippy::cast_precision_loss,
                        reason = "a sample index, far below the mantissa"
                    )]
                    let f = i as f64 / PER_CURVE as f64;
                    let t = (range.1 - range.0).mul_add(f, range.0);
                    out.push((placement.apply(curve.point_at(t, tol)?), POINT_SHARE / 4.0));
                }
            }
            other => ogeom_bail!(
                Construction,
                "constraint {k} is a {other:?}; a filling passes through vertices and edges"
            ),
        }
    }
    Ok(out)
}

/// The basis curve of a forward trim, whose parameter is the trim's own.
fn untrimmed(curve: &Curve) -> &Curve {
    match curve {
        Curve::Trimmed(t) if !t.is_reversed() => untrimmed(t.basis()),
        other => other,
    }
}

/// The side's curve seen in the plane, at the curve's own parameter: exact
/// where the curve is a line, a circle or ellipse, or a B-spline (an affine
/// image of each is a curve of the same form), checked against the curve
/// before it is trusted; fitted at the curve's parameters otherwise.
fn trace(
    side: &Side,
    frame: &PlaneFrame,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let truth = |t: f64| -> OgeomResult<Point2> { Ok(frame.chart(side.curve.point_at(t, tol)?)) };
    let (t0, t1) = side.range;
    let exact: Option<PlanarCurve> = match untrimmed(&side.curve) {
        Curve::BSpline(spline) => {
            let mut control = Vec::with_capacity(spline.control_points().len());
            for w in spline.control_points() {
                control.push(Weighted::new(frame.chart(w.point()), w.weight, tol)?);
            }
            BSpline2d::rational(spline.knots().clone(), control)
                .ok()
                .map(PlanarCurve::BSpline)
        }
        _ => match side.curve.kind() {
            CurveKind::Line => {
                let (q0, q1) = (truth(t0)?, truth(t1)?);
                let d = (q1 - q0) * (1.0 / (t1 - t0));
                let c = q0 - d * t0;
                Trig2d::new(c, d, Vector2::ZERO, Vector2::ZERO, side.range)
                    .ok()
                    .map(PlanarCurve::Trig)
            }
            CurveKind::Circle | CurveKind::Ellipse => {
                // Thirds rather than ends and middle: a full circle's ends
                // are one point.
                let ts = [
                    t0,
                    (t1 - t0).mul_add(1.0 / 3.0, t0),
                    (t1 - t0).mul_add(2.0 / 3.0, t0),
                ];
                let rows = ts.map(|t| [1.0, t.cos(), t.sin()]);
                let values = [truth(ts[0])?, truth(ts[1])?, truth(ts[2])?];
                match (
                    solve3(rows, values.map(|q| q.x)),
                    solve3(rows, values.map(|q| q.y)),
                ) {
                    (Some(x), Some(y)) => Trig2d::new(
                        Point2::new(x[0], y[0]),
                        Vector2::ZERO,
                        Vector2::new(x[1], y[1]),
                        Vector2::new(x[2], y[2]),
                        side.range,
                    )
                    .ok()
                    .map(PlanarCurve::Trig),
                    _ => None,
                }
            }
            _ => None,
        },
    };
    if let Some(curve) = exact {
        let mut worst = 0.0f64;
        let mut scale = 1.0f64;
        for k in 0..=32 {
            let t = (t1 - t0).mul_add(f64::from(k) / 32.0, t0);
            let q = truth(t)?;
            scale = scale.max(q.to_vector().magnitude());
            worst = worst.max(curve.point_at(t, tol)?.distance(q));
        }
        if worst <= 1e-12 * scale + 1e-3 * tol.confusion() {
            return Ok(curve);
        }
    }
    const SAMPLES: usize = 96;
    let mut parameters = Vec::with_capacity(SAMPLES + 1);
    let mut points = Vec::with_capacity(SAMPLES + 1);
    for k in 0..=SAMPLES {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a sample index, far below the mantissa"
        )]
        let t = (t1 - t0).mul_add(k as f64 / SAMPLES as f64, t0);
        parameters.push(t);
        points.push(truth(t)?);
    }
    let fitted = ogeom_geom::fit::fit_points_2d_at(&parameters, &points, 3, tolerance * 1e-2, tol)?;
    Ok(PlanarCurve::BSpline(fitted.curve))
}

/// Solve a 3x3 system by Cramer's rule; `None` where it is singular.
fn solve3(rows: [[f64; 3]; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * m[1][1].mul_add(m[2][2], -m[1][2] * m[2][1])
            - m[0][1] * m[1][0].mul_add(m[2][2], -m[1][2] * m[2][0])
            + m[0][2] * m[1][0].mul_add(m[2][1], -m[1][1] * m[2][0])
    };
    let d = det(rows);
    if d.abs() < 1e-12 {
        return None;
    }
    let mut out = [0.0; 3];
    for (c, slot) in out.iter_mut().enumerate() {
        let mut m = rows;
        for r in 0..3 {
            m[r][c] = rhs[r];
        }
        *slot = det(m) / d;
    }
    Some(out)
}

/// A side's length, by chords.
fn side_length(side: &Side, tol: Tolerances) -> OgeomResult<f64> {
    let mut length = 0.0;
    let mut previous = side.curve.point_at(side.range.0, tol)?;
    for k in 1..=OUTLINE_SAMPLES {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a sample index, far below the mantissa"
        )]
        let f = k as f64 / OUTLINE_SAMPLES as f64;
        let p = side.curve.point_at(side.parameter(f), tol)?;
        length += p.distance(previous);
        previous = p;
    }
    Ok(length)
}

/// The fit's conditions along one side: positions always, tangent planes
/// from G1 on, normal curvatures at G2. Every row carries the side's share
/// of the boundary per sample (`share`, its length against the hole's
/// `size`) and is scaled to a residual free of units.
fn side_conditions(
    side: &Side,
    share: f64,
    size: f64,
    samples: usize,
    frame: &PlaneFrame,
    out: &mut Vec<Condition>,
    tol: Tolerances,
) -> OgeomResult<()> {
    #[expect(
        clippy::cast_precision_loss,
        reason = "a sample count, far below the mantissa"
    )]
    let root = (share / samples as f64).max(f64::MIN_POSITIVE).sqrt();
    for k in 0..=samples {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a sample index, far below the mantissa"
        )]
        let t = side.parameter(k as f64 / samples as f64);
        let p = side.curve.point_at(t, tol)?;
        let at = frame.chart(p);
        out.push(Condition {
            at,
            order: (0, 0),
            target: frame.height(p),
            weight: root / size,
        });
        if side.order == 0 {
            continue;
        }
        let (Some((uv, normal)), Some(support)) = (side.support_at(t, tol)?, &side.support) else {
            continue;
        };
        let lift = normal.dot(frame.n);
        if lift.abs() < MIN_LIFT {
            ogeom_bail!(
                Construction,
                "side {}'s support stands within {:.1} degrees of square to the \
                 plane the boundary spans; the hole is not a height field over it",
                side.entry,
                lift.abs().asin().to_degrees()
            );
        }
        // The support's normal on the plane's side, and the slopes that
        // put the patch's tangent plane on it.
        let normal = if lift < 0.0 { normal * -1.0 } else { normal };
        let lift = lift.abs();
        let (hu, hv) = (-normal.dot(frame.e1) / lift, -normal.dot(frame.e2) / lift);
        out.push(Condition {
            at,
            order: (1, 0),
            target: hu,
            weight: root,
        });
        out.push(Condition {
            at,
            order: (0, 1),
            target: hv,
            weight: root,
        });
        if side.order < 2 {
            continue;
        }
        // With the tangent plane fixed, the patch's second fundamental form
        // is `n · S_ab = h_ab (normal · n)`: the support's, read off its
        // principal curvatures, fixes the three second derivatives.
        let curvature = support.surface.curvature_at(uv.x, uv.y, tol)?;
        let sense = if curvature.normal.dot_vector(normal) < 0.0 {
            -1.0
        } else {
            1.0
        };
        let (dmax, dmin) = (
            curvature.max_direction.vector(),
            curvature.min_direction.vector(),
        );
        let second = |a: Vector, b: Vector| {
            sense
                * curvature.max.mul_add(
                    a.dot(dmax) * b.dot(dmax),
                    curvature.min * a.dot(dmin) * b.dot(dmin),
                )
                / lift
        };
        let su = frame.e1 + frame.n * hu;
        let sv = frame.e2 + frame.n * hv;
        for (order, target) in [
            ((2, 0), second(su, su)),
            ((1, 1), second(su, sv)),
            ((0, 2), second(sv, sv)),
        ] {
            out.push(Condition {
                at,
                order,
                target,
                weight: root * size,
            });
        }
    }
    Ok(())
}

/// Measure one side at `stations` spread over its edge: the gap between
/// the edge's curve and the surface through the pcurve, and against the
/// support the angle between normals and the difference of normal
/// curvatures square to the edge.
fn measure_side(
    side: &Side,
    pcurve: &PlanarCurve,
    surface: &BSplineSurface,
    stations: usize,
    tol: Tolerances,
) -> OgeomResult<FillSide> {
    let (mut gap, mut angle, mut curvature) = (0.0f64, 0.0f64, None::<f64>);
    for k in 0..stations {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = side.parameter(k as f64 / (stations - 1) as f64);
        let p = side.curve.point_at(t, tol)?;
        let uv = pcurve.point_at(t, tol)?;
        gap = gap.max(surface.point_at(uv.x, uv.y, tol)?.distance(p));
        let (Some((suv, theirs)), Some(support)) = (side.support_at(t, tol)?, &side.support) else {
            continue;
        };
        let ours = surface.normal_at(uv.x, uv.y, tol)?.vector();
        angle = angle.max(ours.cross(theirs).magnitude().atan2(ours.dot(theirs).abs()));
        let across = ours.cross(side.curve.d1_at(t, tol)?);
        let signed = |c: SurfaceCurvature| -> Option<f64> {
            let k = c.normal_curvature(across)?;
            Some(if c.normal.dot_vector(ours) < 0.0 {
                -k
            } else {
                k
            })
        };
        let a = surface.curvature_at(uv.x, uv.y, tol).ok().and_then(signed);
        let b = support
            .surface
            .curvature_at(suv.x, suv.y, tol)
            .ok()
            .and_then(signed);
        if let (Some(a), Some(b)) = (a, b) {
            curvature = Some(curvature.unwrap_or(0.0).max((a - b).abs()));
        }
    }
    let supported = side.support.is_some();
    Ok(FillSide {
        edge: side.edge.clone(),
        gap,
        angle: supported.then_some(angle),
        curvature: curvature.filter(|_| supported),
        stations,
    })
}

/// The first deviation past `tolerance`, described, or `None` when every
/// side and constraint meets it.
fn first_miss(
    sides: &[Side],
    reports: &[FillSide],
    constraint_gap: f64,
    tolerance: f64,
) -> Option<String> {
    for (side, report) in sides.iter().zip(reports) {
        let i = side.entry;
        if report.gap > tolerance {
            return Some(format!("side {i} stands {} off its edge", report.gap));
        }
        if side.order >= 1 {
            let angle = report.angle.unwrap_or(f64::INFINITY);
            if angle > tolerance {
                return Some(format!(
                    "side {i} meets its support {angle} radians from tangent"
                ));
            }
        }
        if side.order >= 2 {
            match report.curvature {
                Some(c) if c <= tolerance => {}
                Some(c) => {
                    return Some(format!(
                        "side {i}'s normal curvature differs from its support's by {c}"
                    ));
                }
                None => {
                    return Some(format!(
                        "side {i}'s curvature could not be read on both surfaces"
                    ));
                }
            }
        }
    }
    (constraint_gap > tolerance)
        .then(|| format!("a constraint stands {constraint_gap} off the surface"))
}

/// Build the face on the fitted patch, bounded by the sides' own edges,
/// each given its pcurve on it and a tolerance that holds its gap.
fn build(
    model: &mut Model,
    sides: &[Side],
    order: &[usize],
    traces: &[PlanarCurve],
    reports: &[FillSide],
    surface: BSplineSurface,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let walk: Vec<Shape> = order
        .iter()
        .map(|&i| {
            let side = &sides[i];
            if side.reversed {
                side.edge.reversed()
            } else {
                side.edge.clone()
            }
        })
        .collect();
    let wire = make_wire(model, &walk, tol)?.shape;
    let face = make_face(model, SurfaceGeometry::BSpline(surface), &[wire], tol)?.shape;
    let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "the face just built is not in this model");
    };
    let surface_id = data.surface;
    for ((side, pcurve), report) in sides.iter().zip(traces).zip(reports) {
        attach_pcurve(
            model,
            &side.edge,
            pcurve.clone(),
            surface_id,
            Location::identity(),
            side.range,
        )?;
        if report.gap > side.edge_tolerance {
            // The edge owns the gap, and so do the vertices bounding it.
            let widened = Tolerance::new(report.gap + tol.confusion())?;
            model.widen(&side.edge, widened)?;
            if let Some((a, b)) = edge_vertices(model, &side.edge)? {
                model.widen(&a, widened)?;
                model.widen(&b, widened)?;
            }
        }
    }
    Ok(face)
}
