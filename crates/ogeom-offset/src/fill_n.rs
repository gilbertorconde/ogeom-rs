//! The N-sided filling: one face over a hole bounded by any number of
//! edges, meeting each side's face at G0, G1 or G2.

use ogeom_algo::{
    Built, History, attach_pcurve, edge_vertices, make_edge_between, make_face, make_vertex,
    make_wire,
};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSpline2d, BSplineSurface, Continuity, Curve, Curve2d as _, Curve3d as _, CurveKind,
    PlanarCurve, Surface as _, SurfaceCurvature, SurfaceGeometry, Transformable as _, Trig2d,
};
use ogeom_math::{Direction, Point, Point2, Vector, Vector2, Weighted};
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Orientation, Shape, ShapeType, explore,
};

use crate::fill_patch::{Condition, DEGREE, PlaneFrame, fit_free, fit_height};

/// One side of an N-sided filling.
#[derive(Debug, Clone)]
pub struct FillBoundary {
    /// The boundary edge, placed or not. The filling's face is bounded by
    /// this edge node itself where the edge is not placed and its ends are
    /// the vertex nodes its neighbours' ends are; otherwise by a new edge
    /// on the edge's curve, where it stands, between vertices shared with
    /// the neighbouring sides, which the history records as generated from
    /// this edge. Either way the face sews to the edge's own faces.
    pub edge: Shape,
    /// The face the edge belongs to, placed or not, which the filling
    /// meets across it; the edge as an edge of this face, at the same
    /// placement.
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
    /// The edge bounding the filling along this side: the side's own edge,
    /// or the new edge standing in for it.
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
/// The share of the tolerance the bending energy may cost the conditions;
/// see [`smoothing_for`].
const BUDGET: f64 = 0.1;
/// The least weight of the bending energy: below it rounding rather than
/// the energy would settle the controls the conditions leave free.
const LEAST_SMOOTHING: f64 = 1e-12;
/// The least cosine between a support's normal and the plane's normal:
/// about 84 degrees.
const MIN_LIFT: f64 = 0.1;
/// A constraint point's weight, as the share of the hole's size a boundary
/// sample of that length would carry.
const POINT_SHARE: f64 = 0.05;
/// Samples per constraint curve.
const PER_CURVE: usize = 32;
/// The least share of each corner's turn the loop keeps seen along its
/// plane's normal for a height over that plane to fill it; below it the
/// free patch does.
const KEEP: f64 = 0.25;
/// Samples per side for the chart the free patch is drawn over.
const CHART_SAMPLES: usize = 128;
/// The largest share of a half turn a corner of that chart turns.
const MOST_TURN: f64 = 0.9;
/// The least angle, in radians, between the tangents meeting at a corner
/// of the loop for it to count as one: about a degree.
const CORNER: f64 = 0.0175;

/// Fill the hole a loop of edges bounds with one face meeting each side's
/// support at the continuity asked.
///
/// The sides may come in any order and either direction; they must chain
/// into one simple closed loop, each side's end meeting the next one's at
/// a shared vertex node or where the two vertices lie within their own and
/// their edges' tolerances of each other. `constraints` are vertices and
/// edges inside the hole that the surface passes through. The surface is a
/// cubic B-spline height patch over the plane the loop spans, fitted by least
/// squares to the sides' positions, to the tangent planes of G1 and G2
/// sides' supports and to the normal curvatures of G2 sides' supports, with
/// a thin-plate bending energy settling the rest, weighted from
/// `tolerance` so it costs the conditions a small share of it; the control
/// net is refined until every side meets `tolerance` and the conditions'
/// residual no longer outweighs the energy, or the refinement runs out.
/// Where a corner of the loop is seen smooth along the plane's normal (two
/// sides meeting at an angle, each in a plane through that normal), no
/// height over the plane meets both sides there: the patch is then fitted
/// in all three coordinates over a chart drawn from the loop itself, whose
/// corners are the loop's, and refined until it also does not fold over
/// inside the hole.
/// The face is trimmed by the given edges themselves where they are not
/// placed and share vertex nodes end to end, each given a pcurve on the
/// patch; any other side is stood in for by a new edge on its curve
/// between vertices shared round the loop, so the face is one closed wire
/// either way. A placed edge or support is read where it stands. The face
/// faces the way the supports say: across each edge it
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
/// - a side is not an edge, or has no 3D curve or no vertices;
/// - the same edge is given twice;
/// - a side asks [`Continuity::C1`], [`Continuity::C2`] or
///   [`Continuity::CInfinity`] (parametric continuity between two
///   surfaces' charts), or G1 or G2 with no support;
/// - a support is not a face, does not hold its edge, or the
///   edge does not lie on its surface within `tolerance`;
/// - the sides do not chain into one closed loop;
/// - the loop encloses no area, or crosses itself seen along the normal of
///   the plane it spans (the plane it encloses the most area seen square
///   to), so the hole is not a height field over that plane;
/// - a G1 or G2 side's support stands within about 6 degrees of square to
///   that plane;
/// - a constraint is neither a vertex nor an edge, or lies outside the hole
///   seen along the plane's normal, or is given where a corner of the loop
///   is seen smooth along it;
/// - the loop's corners turn so far that no chart drawn from it closes.
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
            if sides[i].given.is_same(&sides[j].given) {
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

    // A corner seen smooth along the plane's normal asks a height over the
    // plane for two slopes at one point, one from each side; the free patch
    // over a chart drawn from the loop, whose corners are the loop's, takes
    // it.
    let corners = corners_of(&sides, &order, tol)?;
    let free = kept_turn(frame.n, &corners) < KEEP;
    if free && !interior.is_empty() {
        ogeom_bail!(
            Construction,
            "a corner of the boundary loop is seen smooth along the normal of \
             the plane it spans, so the filling is drawn over a chart of the \
             loop's own, which places no interior constraints"
        );
    }
    let mut traces = Vec::with_capacity(sides.len());
    let mut lengths = Vec::with_capacity(sides.len());
    for side in &sides {
        if !free {
            traces.push(trace(side, &frame, tolerance, tol)?);
        }
        lengths.push(side_length(side, tol)?);
    }
    let chart_outline = if free {
        traces = boundary_chart(&sides, &order, &corners, tolerance, tol)?;
        traced_outline(&sides, &order, &traces, tol)?
    } else {
        chart_outline
    };

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

    let smoothing = smoothing_for(&sides, tolerance, size);
    let mut last_miss = String::new();
    // The latest fit that met the tolerance.
    let mut fallback = None;
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

        let fit = if free {
            let mut conditions = Vec::new();
            for ((side, length), pcurve) in sides.iter().zip(&lengths).zip(&traces) {
                free_conditions(
                    side,
                    pcurve,
                    *length / size,
                    size,
                    samples,
                    &mut conditions,
                    tol,
                )?;
            }
            fit_free(domain, controls, &conditions, smoothing, tol)?
        } else {
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
                conditions.push(Condition::partial(
                    frame.chart(*p),
                    (0, 0),
                    [frame.height(*p)],
                    share.sqrt() / size,
                ));
            }
            fit_height(&frame, domain, controls, &conditions, smoothing, tol)?
        };
        let surface = fit.surface;

        let mut reports = Vec::with_capacity(sides.len());
        for (side, pcurve) in sides.iter().zip(&traces) {
            reports.push(measure_side(side, pcurve, &surface, 3 * samples + 1, tol)?);
        }
        let mut constraint_gap = 0.0f64;
        for (p, _) in &interior {
            let q = frame.chart(*p);
            constraint_gap = constraint_gap.max(surface.point_at(q.x, q.y, tol)?.distance(*p));
        }

        if let Some(miss) = first_miss(&sides, &reports, constraint_gap, tolerance) {
            last_miss = format!("at {}x{} controls, {miss}", controls.0, controls.1);
            continue;
        }
        if free && folds(&surface, &chart_outline, domain, tol)? {
            last_miss = format!(
                "at {}x{} controls, the patch folds over inside the hole",
                controls.0, controls.1
            );
            continue;
        }
        fallback = Some((surface, reports, constraint_gap));
        // A residual-led net is refined; see [`smoothing_for`].
        if fit.residual <= smoothing {
            break;
        }
    }
    let Some((surface, reports, constraint_gap)) = fallback else {
        ogeom_bail!(
            NotDone,
            "the filling misses its tolerance of {tolerance} on the finest net: {last_miss}"
        )
    };
    let face = build(model, &mut sides, &order, &traces, &reports, surface, tol)?;
    let mut history = History::new();
    let mut reports = reports;
    for (side, report) in sides.iter().zip(&mut reports) {
        history.generate(&side.given, face.clone());
        if !side.edge.is_same(&side.given) {
            history.generate(&side.given, side.edge.clone());
        }
        report.edge = side.edge.clone();
    }
    for constraint in constraints {
        history.generate(constraint, face.clone());
    }
    Ok(Filled {
        built: Built::new(face, history),
        sides: reports,
        constraint_gap,
    })
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
    /// The edge as given, forward.
    given: Shape,
    /// The edge bounding the face: the given one, or its stand-in once the
    /// face is built.
    edge: Shape,
    /// Whether the given edge or its curve is placed.
    placed: bool,
    /// The curve where the edge stands.
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
    let Some(EdgeRepr::Curve3d {
        curve,
        range,
        location: own,
    }) = data.curve3d()
    else {
        ogeom_bail!(Construction, "side {i}'s edge has no 3D curve");
    };
    let Some(curve) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "side {i}'s curve is not in this model");
    };
    let range = *range;
    let edge_tolerance = data.tolerance.get();
    let placed = !(edge.location().is_identity() && own.is_identity());
    let curve = if placed {
        curve
            .transformed(&own.composed(model.datums())?, tol)?
            .transformed(&edge.transform(model.datums())?, tol)?
    } else {
        curve
    };
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
        given: edge.clone(),
        edge,
        placed,
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
    let face_placed = !(face.location().is_identity() && face_data.location.is_identity());
    let holds = explore(model, face, Filter::OfType(ShapeType::Edge))?
        .iter()
        .any(|e| e.is_same(edge));
    if !holds {
        ogeom_bail!(
            Construction,
            "side {i}'s support face does not hold the side's edge at its placement"
        );
    }
    let surface_id = face_data.surface;
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "side {i}'s support surface is not in this model");
    };
    let surface = if face_placed {
        surface
            .transformed(&face_data.location.composed(model.datums())?, tol)?
            .transformed(&face.transform(model.datums())?, tol)?
    } else {
        surface
    };
    let Some(edge_data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "side {i}'s edge holds no edge data");
    };
    // A stored pcurve describes the surface's own chart, which a placed
    // face's restated surface need not share.
    let stored = match edge_data
        .pcurve_for(surface_id, edge.location())
        .filter(|_| !face_placed)
    {
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
    // The edge must lie on the surface it is said to bound. A stored image
    // that misses it (one kept for another placement of the same edge
    // node) gives way to one fitted where the edge stands.
    let off = |pcurve: &PlanarCurve, prange: (f64, f64)| -> OgeomResult<f64> {
        let mut worst = 0.0f64;
        for k in 0..=16 {
            let f = f64::from(k) / 16.0;
            let t = (range.1 - range.0).mul_add(f, range.0);
            let pt = (prange.1 - prange.0).mul_add(f, prange.0);
            let uv = pcurve.point_at(pt, tol)?;
            let on = surface.point_at(uv.x, uv.y, tol)?;
            worst = worst.max(on.distance(curve.point_at(t, tol)?));
        }
        Ok(worst)
    };
    let stored = match stored {
        Some((pcurve, prange)) => {
            let worst = off(&pcurve, prange)?;
            (worst <= reach).then_some((pcurve, prange, worst))
        }
        None => None,
    };
    let (pcurve, prange, worst) = if let Some(found) = stored {
        found
    } else {
        let (fitted, _, _, _, _) =
            ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, range, &surface, tol)?;
        let worst = off(&fitted, range)?;
        (fitted, range, worst)
    };
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

/// One end of a side: its vertex, where the vertex stands, and how far
/// it reaches (its own tolerance or its edge's, the wider).
#[derive(Clone)]
struct End {
    vertex: Shape,
    at: Point,
    reach: f64,
}

/// A side's two ends, in the edge's own direction.
fn ends_of(model: &Model, side: &Side) -> OgeomResult<(End, End)> {
    let Some((a, b)) = edge_vertices(model, &side.given)? else {
        ogeom_bail!(
            Construction,
            "side {} has no vertices, so it cannot be shown to join the loop",
            side.entry
        );
    };
    let end = |vertex: Shape| -> OgeomResult<End> {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Construction, "side {}'s vertex holds no point", side.entry);
        };
        Ok(End {
            at: vertex.transform(model.datums())?.apply(data.point),
            reach: data.tolerance.get().max(side.edge_tolerance),
            vertex,
        })
    };
    Ok((end(a)?, end(b)?))
}

/// Chain the sides into one loop: the order they are walked in, with each
/// side's `reversed` set to the direction it is walked. Ends join where
/// they are one vertex, or failing that where they lie within reach of
/// each other.
fn chain(model: &Model, sides: &mut [Side], tol: Tolerances) -> OgeomResult<Vec<usize>> {
    let mut ends = Vec::with_capacity(sides.len());
    for side in sides.iter() {
        let (a, b) = ends_of(model, side)?;
        ends.push((a, b));
    }
    let same = |a: &End, b: &End| -> OgeomResult<bool> {
        Ok(a.vertex.is_same(&b.vertex) || model.same_position(&a.vertex, &b.vertex, tol)?)
    };
    let near = |a: &End, b: &End| a.at.distance(b.at) <= a.reach + b.reach + tol.confusion();
    let meets = |a: &End, b: &End| -> OgeomResult<bool> { Ok(same(a, b)? || near(a, b)) };
    let n = sides.len();
    let mut used = vec![false; n];
    used[0] = true;
    let mut order = vec![0];
    let first = ends[0].0.clone();
    let mut cursor = ends[0].1.clone();
    while order.len() < n {
        let last = sides[order[order.len() - 1]].entry;
        // A shared vertex is the stronger word: ends merely near the
        // cursor count only where none is the cursor's own vertex.
        let mut found: Vec<(usize, bool)> = Vec::new();
        for strict in [true, false] {
            for j in 0..n {
                if used[j] {
                    continue;
                }
                let joins = |e: &End| -> OgeomResult<bool> {
                    if strict {
                        same(e, &cursor)
                    } else {
                        meets(e, &cursor)
                    }
                };
                if joins(&ends[j].0)? {
                    found.push((j, false));
                } else if joins(&ends[j].1)? {
                    found.push((j, true));
                }
            }
            if !found.is_empty() {
                break;
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
                .filter(|e| e.is_same(&side.given))
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

/// The loop's corners in its walking order, corner `k` where side
/// `order[k]` ends and side `order[k + 1]` begins: the unit tangents
/// arriving and leaving, and the angle between them, zero where the loop
/// runs on smoothly (within [`CORNER`]) or a tangent vanishes.
fn corners_of(
    sides: &[Side],
    order: &[usize],
    tol: Tolerances,
) -> OgeomResult<Vec<(Vector, Vector, f64)>> {
    let walked = |side: &Side, end: bool| -> OgeomResult<Option<Vector>> {
        let t = if end == side.reversed {
            side.range.0
        } else {
            side.range.1
        };
        let d = side.curve.d1_at(t, tol)?;
        let d = if side.reversed { d * -1.0 } else { d };
        Ok(Direction::new(d, tol).ok().map(Direction::vector))
    };
    let mut out = Vec::with_capacity(order.len());
    for (k, &i) in order.iter().enumerate() {
        let next = &sides[order[(k + 1) % order.len()]];
        let corner = match (walked(&sides[i], true)?, walked(next, false)?) {
            (Some(a), Some(b)) => {
                let turn = a.cross(b).magnitude().atan2(a.dot(b));
                (a, b, if turn > CORNER { turn } else { 0.0 })
            }
            _ => (Vector::ZERO, Vector::ZERO, 0.0),
        };
        out.push(corner);
    }
    Ok(out)
}

/// The least share of its turn any corner of the loop still turns seen
/// along the unit `n`; one where the loop has no corners.
fn kept_turn(n: Vector, corners: &[(Vector, Vector, f64)]) -> f64 {
    let flat = |v: Vector| v - n * v.dot(n);
    let mut worst = 1.0f64;
    for (a, b, turn) in corners {
        if *turn > 0.0 {
            let (a, b) = (flat(*a), flat(*b));
            let seen = a.cross(b).magnitude().atan2(a.dot(b));
            worst = worst.min((seen / turn).min(1.0));
        }
    }
    worst
}

/// A chart drawn from the loop itself, for the free patch: each side's
/// pcurve, by side index, at the side's own parameter.
///
/// The chart walks the loop at the loop's own speed, turning at each corner
/// by the corner's turn (held under [`MOST_TURN`] of a half turn, and
/// scaled down together where they would leave the sides no turn of their
/// own) and spreading the rest of a full turn evenly along its length, so
/// its corners are the loop's corners and it is smooth wherever the loop
/// is. Such a walk ends short of where it began by a little; that drift is
/// taken out evenly along it, one vector off every step's velocity, which
/// leaves the walk smooth where it was smooth and its corners corners.
fn boundary_chart(
    sides: &[Side],
    order: &[usize],
    corners: &[(Vector, Vector, f64)],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<PlanarCurve>> {
    use std::f64::consts::{PI, TAU};
    // Each side at even parameters, with the length walked from where the
    // loop enters the side.
    let mut walks: Vec<(Vec<f64>, Vec<f64>, f64)> =
        vec![(Vec::new(), Vec::new(), 0.0); sides.len()];
    let mut total = 0.0;
    for &i in order {
        let side = &sides[i];
        let mut params = Vec::with_capacity(CHART_SAMPLES + 1);
        let mut along = Vec::with_capacity(CHART_SAMPLES + 1);
        let mut previous = side.curve.point_at(side.range.0, tol)?;
        let mut length = 0.0;
        for k in 0..=CHART_SAMPLES {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a sample index, far below the mantissa"
            )]
            let t = side.parameter(k as f64 / CHART_SAMPLES as f64);
            let p = side.curve.point_at(t, tol)?;
            length += p.distance(previous);
            previous = p;
            params.push(t);
            along.push(length);
        }
        if side.reversed {
            for a in &mut along {
                *a = length - *a;
            }
        }
        total += length;
        walks[i] = (params, along, length);
    }
    if total <= tol.confusion() {
        ogeom_bail!(Construction, "the boundary loop has no length");
    }
    let turns: Vec<f64> = corners.iter().map(|c| c.2.min(MOST_TURN * PI)).collect();
    let sum: f64 = turns.iter().sum();
    let scale = if sum > MOST_TURN * TAU {
        MOST_TURN * TAU / sum
    } else {
        1.0
    };
    let bend = sum.mul_add(-scale, TAU) / total;
    // The point a walk of length `s` from `start`, heading `heading`,
    // reaches turning at `bend`.
    let arc = |start: Point2, heading: f64, s: f64| -> Point2 {
        let delta = bend * s;
        if delta.abs() < 1e-9 {
            start + Vector2::new(heading.cos(), heading.sin()) * s
        } else {
            start
                + Vector2::new(
                    (heading + delta).sin() - heading.sin(),
                    heading.cos() - (heading + delta).cos(),
                ) * (1.0 / bend)
        }
    };
    let mut starts = vec![(Point2::new(0.0, 0.0), 0.0, 0.0); sides.len()];
    let (mut at, mut heading, mut walked) = (Point2::new(0.0, 0.0), 0.0f64, 0.0f64);
    for (k, &i) in order.iter().enumerate() {
        starts[i] = (at, heading, walked);
        let length = walks[i].2;
        at = arc(at, heading, length);
        heading += bend.mul_add(length, turns[k] * scale);
        walked += length;
    }
    let drift = at - Point2::new(0.0, 0.0);
    if drift.magnitude() > 0.5 * total / TAU {
        ogeom_bail!(
            Construction,
            "the boundary loop's corners leave no chart drawn from it closing"
        );
    }
    let mut out = Vec::with_capacity(sides.len());
    for (i, (params, along, _)) in walks.iter().enumerate() {
        let (start, heading, before) = starts[i];
        let points: Vec<Point2> = along
            .iter()
            .map(|&s| arc(start, heading, s) - drift * ((before + s) / total))
            .collect();
        let fitted = ogeom_geom::fit::fit_points_2d_at(params, &points, 3, tolerance * 1e-2, tol)?;
        out.push(PlanarCurve::BSpline(fitted.curve));
    }
    Ok(out)
}
/// Whether a closed polygon crosses itself: any two segments that are not
/// neighbours crossing, each segment's ends standing clear of the other's
/// line on opposite sides.
///
/// A point within rounding of a line is on it, not on a side: samples of a
/// straight side are collinear, and two segments of one straight side
/// would otherwise cross wherever rounding scatters their ends' sides.
fn crosses_itself(polygon: &[Point2]) -> bool {
    let n = polygon.len();
    let reach = polygon
        .iter()
        .map(|q| q.x.abs().max(q.y.abs()))
        .fold(0.0f64, f64::max);
    let floor = 1e-9 * reach;
    // The side of line `ab` that `c` stands on: 1, -1, or 0 within `floor`.
    let side = |a: Point2, b: Point2, c: Point2| -> i8 {
        let ab = b - a;
        let o = ab.cross(c - a);
        if o.abs() <= floor * ab.magnitude() {
            0
        } else if o > 0.0 {
            1
        } else {
            -1
        }
    };
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (c, d) = (polygon[j], polygon[(j + 1) % n]);
            if side(a, b, c) * side(a, b, d) < 0 && side(c, d, a) * side(c, d, b) < 0 {
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

/// The bending energy's weight against the conditions.
///
/// The fit minimises `Q + λ·E`. `Q` sums the conditions' squared
/// residuals, each free of units (a position over the patch size `size`, a
/// slope as it is, a curvature times `size`) and weighted by its side's
/// share of the boundary. `E` is the thin-plate energy, also free of units:
/// about the square of the angle, in radians, the patch turns through.
/// Scaling the whole problem changes neither, so `λ` means the same at any
/// size and on any net; only the tolerance sets it.
///
/// For any surface `s` the net can carry, `Q(fit) + λ·E(fit) <= Q(s) +
/// λ·E(s)`: the energy costs the conditions at most `λ·E(s)`. With `ε` the
/// tolerance in the residuals' own units (the strictest of `tolerance /
/// size` for positions, `tolerance` for G1 slopes and `tolerance · size`
/// for G2 curvatures), `λ = (BUDGET·ε)²` holds that cost to a tenth of the
/// tolerance, root mean square, per radian of turning, and gives the energy
/// the largest weight the bound allows, so the energy rather than the
/// conditions' approximation error settles every control the conditions
/// reach only weakly.
///
/// The same bound says when a net is too coarse. Where the conditions'
/// residual `Q` exceeds `λ`, the fit would bend the interior by up to a
/// unit of energy to save that residual, so the hole's middle follows the
/// boundary's approximation error (on a coarse net every control reaches
/// the boundary, and the middle dips or bulges). Such a net is refined
/// even when its sides meet the tolerance; where no net gets clear of it,
/// the finest one that met the tolerance is taken.
fn smoothing_for(sides: &[Side], tolerance: f64, size: f64) -> f64 {
    let order = sides.iter().map(|s| s.order).max().unwrap_or(0);
    let mut eps = tolerance / size;
    if order >= 1 {
        eps = eps.min(tolerance);
    }
    if order >= 2 {
        eps = eps.min(tolerance * size);
    }
    (BUDGET * eps).powi(2).max(LEAST_SMOOTHING)
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
    out: &mut Vec<Condition<1>>,
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
        out.push(Condition::partial(
            at,
            (0, 0),
            [frame.height(p)],
            root / size,
        ));
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
        out.push(Condition::partial(at, (1, 0), [hu], root));
        out.push(Condition::partial(at, (0, 1), [hv], root));
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
            out.push(Condition::partial(at, order, [target], root * size));
        }
    }
    Ok(())
}

/// The free patch's conditions along one side: the side's points, at its
/// pcurve on the chart drawn from the loop, weighted as in
/// [`side_conditions`].
fn free_conditions(
    side: &Side,
    pcurve: &PlanarCurve,
    share: f64,
    size: f64,
    samples: usize,
    out: &mut Vec<Condition<3>>,
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
        let at = pcurve.point_at(t, tol)?;
        out.push(Condition::partial(at, (0, 0), [p.x, p.y, p.z], root / size));
    }
    Ok(())
}

/// Points round the loop on the chart, in its walking order.
fn traced_outline(
    sides: &[Side],
    order: &[usize],
    traces: &[PlanarCurve],
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
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
            out.push(traces[i].point_at(side.parameter(f), tol)?);
        }
    }
    Ok(out)
}

/// Whether a free patch folds over inside the hole: on a grid over the
/// chart, inside its outline, the normal vanishes against the largest, or
/// turns over between neighbouring samples.
fn folds(
    surface: &BSplineSurface,
    outline: &[Point2],
    domain: ((f64, f64), (f64, f64)),
    tol: Tolerances,
) -> OgeomResult<bool> {
    const GRID: usize = 40;
    let ((ua, ub), (va, vb)) = domain;
    let mut normals = vec![None; (GRID + 1) * (GRID + 1)];
    let mut largest = 0.0f64;
    for i in 0..=GRID {
        for j in 0..=GRID {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a grid index, far below the mantissa"
            )]
            let q = Point2::new(
                (ub - ua).mul_add(i as f64 / GRID as f64, ua),
                (vb - va).mul_add(j as f64 / GRID as f64, va),
            );
            if !inside(outline, q) {
                continue;
            }
            let (du, dv) = surface.d1_at(q.x, q.y, tol)?;
            let normal = du.cross(dv);
            largest = largest.max(normal.magnitude());
            normals[i * (GRID + 1) + j] = Some(normal);
        }
    }
    for i in 0..=GRID {
        for j in 0..=GRID {
            let Some(here) = normals[i * (GRID + 1) + j] else {
                continue;
            };
            if here.magnitude() <= 1e-6 * largest {
                return Ok(true);
            }
            let right = if i < GRID {
                normals[(i + 1) * (GRID + 1) + j]
            } else {
                None
            };
            let up = if j < GRID {
                normals[i * (GRID + 1) + j + 1]
            } else {
                None
            };
            if [right, up]
                .into_iter()
                .flatten()
                .any(|n| n.dot(here) <= 0.0)
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
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

/// Build the face on the fitted patch, bounded by the sides' own edges or
/// their stand-ins, each given its pcurve on it and a tolerance that holds
/// its gap.
fn build(
    model: &mut Model,
    sides: &mut [Side],
    order: &[usize],
    traces: &[PlanarCurve],
    reports: &[FillSide],
    surface: BSplineSurface,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    stand_in(model, sides, order, tol)?;
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

/// Give each side the edge the face is bounded by: its own where it is not
/// placed and both its corners are vertex nodes it shares with its
/// neighbours, otherwise a new edge on its curve where it stands, between
/// the loop's corner vertices. A corner keeps the shared vertex where there
/// is one, unplaced; elsewhere it is a new vertex midway between the two
/// ends, reaching both.
fn stand_in(
    model: &mut Model,
    sides: &mut [Side],
    order: &[usize],
    tol: Tolerances,
) -> OgeomResult<()> {
    let n = order.len();
    // Corner k is where side `order[k]` ends its walk and `order[k + 1]`
    // starts its.
    let mut corners: Vec<(Shape, bool)> = Vec::with_capacity(n);
    for k in 0..n {
        let (a, b) = (&sides[order[k]], &sides[order[(k + 1) % n]]);
        let (a0, a1) = ends_of(model, a)?;
        let (b0, b1) = ends_of(model, b)?;
        let (arrive, leave) = (
            if a.reversed { a0 } else { a1 },
            if b.reversed { b1 } else { b0 },
        );
        let (from, to) = (
            a.curve
                .point_at(if a.reversed { a.range.0 } else { a.range.1 }, tol)?,
            b.curve
                .point_at(if b.reversed { b.range.1 } else { b.range.0 }, tol)?,
        );
        if arrive.vertex.is_same(&leave.vertex) && arrive.vertex.location().is_identity() {
            corners.push((arrive.vertex, true));
            continue;
        }
        let vertex = make_vertex(model, from.midpoint(to)).shape;
        let reach = (0.5 * from.distance(to) + tol.confusion())
            .max(a.edge_tolerance)
            .max(b.edge_tolerance);
        model.widen(&vertex, Tolerance::new(reach)?)?;
        corners.push((vertex, false));
    }
    for k in 0..n {
        let side = &sides[order[k]];
        let (start, end) = (&corners[(k + n - 1) % n], &corners[k]);
        if !side.placed && start.1 && end.1 {
            continue;
        }
        let (from, to) = if side.reversed {
            (&end.0, &start.0)
        } else {
            (&start.0, &end.0)
        };
        let edge = make_edge_between(model, side.curve.clone(), side.range, from, to, tol)?.shape;
        model.widen(&edge, Tolerance::new(side.edge_tolerance)?)?;
        sides[order[k]].edge = edge;
    }
    Ok(())
}
