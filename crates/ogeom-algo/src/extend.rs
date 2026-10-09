//! Growing a face past one of its boundary edges.
//!
//! The face is rebuilt with the chosen edge replaced by three new ones: a
//! side up from each of its ends and a far edge across, `length` beyond it.
//! Every other edge of the face is the same edge node as before, so the
//! neighbours that share them are untouched; the chosen edge itself stays
//! in the model for the neighbour across it.
//!
//! A plane grows square to a straight edge in its own metric. On the other
//! surfaces the edge must run along a parameter line, and the face grows
//! along the crossing parameter: an analytic surface on its own equation, a
//! B-spline patch continued past the side of its domain the edge lies on.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSplineSurface, CircleCurve, ConeSurface, Continuity, Curve, Curve2d as _, Curve3d as _,
    CylinderSurface, Line2d, LineCurve, PlanarCurve, PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_math::{Axis, Axis2, Circle, Direction, Direction2, Frame, Point, Point2, Vector2};
use ogeom_topo::{
    EdgeRepr, FaceData, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId,
};

use crate::build::{attach_pcurve, edge_vertices, make_edge_between, make_wire};
use crate::history::{Built, History};

/// How a face continues past the edge it is extended across.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Extension {
    /// The surface continues as itself: an analytic surface along its own
    /// equation, a B-spline patch as the polynomial of the span at the edge,
    /// which is as smooth across the old edge as the patch is anywhere.
    Natural,
    /// The surface continues as the polynomial of its derivatives across the
    /// edge, carried to the least order that gives the join `continuity`:
    /// the first order (straight out along the crossing tangent) for
    /// [`Continuity::C0`], [`Continuity::G1`] and [`Continuity::C1`], the
    /// second for [`Continuity::G2`] and [`Continuity::C2`]. An analytic
    /// surface takes it only where its crossing parameter line is straight
    /// (a plane, a cylinder or cone grown along its rulings), where it is
    /// the surface itself.
    Linear {
        /// The smoothness the join across the old edge keeps.
        continuity: Continuity,
    },
}

/// Grow the face `shape` past its boundary edge `edge` by `length`.
///
/// The length is measured along the surface, square to the edge: on a
/// plane every point of the edge moves that far, square to it; elsewhere it
/// is the arc length of the crossing parameter line through the edge's
/// middle, which is square to the edge on every analytic surface. The
/// result is a new face on the same surface (an analytic one with its
/// domain widened where the extension needs it) or on the patch continued
/// as `mode` says. The edges the face kept are the same nodes, carrying a
/// trim on the new surface where it is a new one.
///
/// A placed face, or a face with placed edges (a prism's far end edges are
/// its profile's edges under the prism's translation), is first rebuilt
/// with its placements baked in, and that rebuilt face is extended: its
/// edges are new nodes on the same geometry, which sew to the old ones.
///
/// History: the face is modified into the new face; `edge` is modified into
/// the far edge and generates the two sides; where the face was baked
/// first, every other edge is modified into its baked twin.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by name,
/// where:
///
/// - `length` is not positive and finite;
/// - `shape` is not a face, or `edge` is not an edge of its outer boundary
///   (an edge on no boundary of it, an edge bounding a hole, a seam the
///   boundary runs along twice, a degenerate edge);
/// - the surface is neither a plane, a cylinder, a cone, a sphere, a torus
///   nor a B-spline patch;
/// - on a plane, the edge is not straight;
/// - elsewhere, the edge does not run along a parameter line, has no trim on
///   the surface, or the face's boundary cannot be walked in its chart;
/// - the extension would wrap the face round a periodic direction onto
///   itself, or reach a cone's apex or a sphere's pole;
/// - a B-spline edge is not on a side of the patch's domain, or the patch
///   closes on itself across it;
/// - a linear extension asks [`Continuity::CInfinity`], or crosses a curved
///   parameter line of an analytic surface.
pub fn extend_face(
    model: &mut Model,
    shape: &Shape,
    edge: &Shape,
    length: f64,
    mode: Extension,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !(length > 0.0 && length.is_finite()) {
        ogeom_bail!(
            Construction,
            "a face is extended by a positive length; got {length}"
        );
    }
    let order = order_of(mode)?;
    let kind = model.kind_of(shape)?;
    if kind != ShapeType::Face {
        ogeom_bail!(Construction, "extend_face extends a face; got a {kind:?}");
    }
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "a face is extended across one of its edges");
    }
    let baked = crate::convert::baked_where_placed(model, shape, tol)?;
    if !baked.history.is_empty() {
        let Some(twin) = baked.history.modified(edge).first().cloned() else {
            ogeom_bail!(Construction, "the edge is not on the face");
        };
        let grown = extend_face(model, &baked.shape, &twin, length, mode, tol)?;
        return Ok(Built::new(grown.shape, baked.history.then(&grown.history)));
    }
    let (surface_id, face_tolerance) = {
        let Some(NodeData::Face(data)) = model.node(shape).map(|n| n.data()) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        (data.surface, data.tolerance)
    };

    // The face's wires and edges as stored, so the rebuilt face walks its
    // boundary exactly as the old one did.
    let wires = model.children_of(&shape.oriented(Orientation::Forward))?;
    let Some(outer) = wires.first().cloned() else {
        ogeom_bail!(
            Construction,
            "the face has no boundary, so no edge to extend it across"
        );
    };
    let mut loops: Vec<Vec<Shape>> = Vec::with_capacity(wires.len());
    for wire in &wires {
        loops.push(model.children_of(&wire.oriented(Orientation::Forward))?);
    }
    let hits: Vec<usize> = loops[0]
        .iter()
        .enumerate()
        .filter(|(_, e)| e.is_same(edge))
        .map(|(i, _)| i)
        .collect();
    let at = match hits.as_slice() {
        [i] => *i,
        [] if loops[1..].iter().flatten().any(|e| e.is_same(edge)) => ogeom_bail!(
            Construction,
            "the edge bounds a hole in the face; a face is extended across its outer boundary"
        ),
        [] => ogeom_bail!(Construction, "the edge is not on the face"),
        _ => ogeom_bail!(
            Construction,
            "the face's boundary runs along the edge twice (a seam); there is no one side \
             to extend it across"
        ),
    };
    let crossed = loops[0][at].clone();
    let Some(NodeData::Edge(edge_data)) = model.node(&crossed).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    if edge_data.degenerate {
        ogeom_bail!(
            Construction,
            "a degenerate edge has no extent to extend the face across"
        );
    }
    let Some((from_vertex, to_vertex)) = edge_vertices(model, &crossed)? else {
        ogeom_bail!(Construction, "the edge has no vertices");
    };
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };

    let plan = match &surface {
        SurfaceGeometry::Plane(plane) => plan_on_plane(
            model,
            *plane,
            &loops[0],
            &crossed,
            &from_vertex,
            &to_vertex,
            length,
            tol,
        )?,
        SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Sphere(_)
        | SurfaceGeometry::Torus(_)
        | SurfaceGeometry::BSpline(_) => plan_on_chart(
            model, shape, &surface, surface_id, &crossed, length, order, tol,
        )?,
        other => ogeom_bail!(
            Construction,
            "a face on a {:?} surface is not extended; only planes, cylinders, cones, \
             spheres, tori and B-spline patches are",
            other.kind()
        ),
    };

    let new_surface_id = match plan.surface {
        Some(widened) => model.geometry_mut().add_surface(widened),
        None => surface_id,
    };
    let Some(target) = model.geometry().surface(new_surface_id).cloned() else {
        ogeom_bail!(Dangling, "the extended surface is not in this model");
    };
    if new_surface_id != surface_id {
        let kept: Vec<Shape> = loops
            .iter()
            .flatten()
            .filter(|e| !e.is_same(&crossed))
            .cloned()
            .collect();
        retrim(model, &kept, surface_id, new_surface_id)?;
    }

    let far_from = model.add_point(target.point_at(plan.far.0.x, plan.far.0.y, tol)?);
    let far_to = model.add_point(target.point_at(plan.far.1.x, plan.far.1.y, tol)?);
    let side_a = chart_edge(
        model,
        &target,
        new_surface_id,
        plan.flat,
        (plan.near.0, &from_vertex),
        (plan.far.0, &far_from),
        tol,
    )?;
    let far = chart_edge(
        model,
        &target,
        new_surface_id,
        plan.flat,
        (plan.far.0, &far_from),
        (plan.far.1, &far_to),
        tol,
    )?;
    let side_b = chart_edge(
        model,
        &target,
        new_surface_id,
        plan.flat,
        (plan.far.1, &far_to),
        (plan.near.1, &to_vertex),
        tol,
    )?;

    let mut ring = loops[0].clone();
    ring.splice(at..=at, [side_a.clone(), far.clone(), side_b.clone()]);
    let new_outer = make_wire(model, &ring, tol)?
        .shape
        .oriented(outer.orientation());
    let mut new_wires = vec![new_outer];
    new_wires.extend(wires[1..].iter().cloned());
    let mut data = FaceData::new(new_surface_id, Location::identity());
    data.tolerance = face_tolerance;
    let face = model
        .add_face(data, &new_wires)?
        .oriented(shape.orientation());

    let mut history = History::new();
    history.modify(shape, face.clone());
    history.modify(edge, far.oriented(Orientation::Forward));
    history.generate(edge, side_a.oriented(Orientation::Forward));
    history.generate(edge, side_b.oriented(Orientation::Forward));
    Ok(Built::new(face, history))
}

/// The order of derivatives a mode carries across the edge; `None` for all
/// of them.
fn order_of(mode: Extension) -> OgeomResult<Option<usize>> {
    match mode {
        Extension::Natural => Ok(None),
        Extension::Linear { continuity } => match continuity {
            Continuity::C0 | Continuity::G1 | Continuity::C1 => Ok(Some(1)),
            Continuity::G2 | Continuity::C2 => Ok(Some(2)),
            Continuity::CInfinity => ogeom_bail!(
                Construction,
                "a linear extension carries derivatives to the second order at most; \
                 a join smooth to every order is the natural extension"
            ),
        },
    }
}

/// Where the new edges go, in the chart of the surface the face ends up on.
struct Plan {
    /// The surface to build on, where it is not the face's own.
    surface: Option<SurfaceGeometry>,
    /// Whether the chart is a plane's, where the new edges are straight
    /// between any two chart points.
    flat: bool,
    /// The crossed edge's ends, in the order the face's boundary walks it.
    near: (Point2, Point2),
    /// Where each end goes.
    far: (Point2, Point2),
}

/// A plane's extension: the straight edge moved square to itself, away
/// from the face, in the plane's own isometric chart.
#[allow(clippy::too_many_arguments)]
fn plan_on_plane(
    model: &Model,
    plane: PlaneSurface,
    outer: &[Shape],
    crossed: &Shape,
    from_vertex: &Shape,
    to_vertex: &Shape,
    length: f64,
    tol: Tolerances,
) -> OgeomResult<Plan> {
    let frame = plane.plane().frame();
    let chart = |p: Point| frame.to_local(p).xy();
    let a = vertex_point(model, from_vertex)?;
    let b = vertex_point(model, to_vertex)?;
    let chord = b - a;
    let span = chord.magnitude();
    if span <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the edge closes on itself; there is no straight run to extend across"
        );
    }
    for p in edge_samples(model, crossed, 8, tol)? {
        let off = (p - a).cross(chord).magnitude() / span;
        if off > tol.confusion() {
            ogeom_bail!(
                Construction,
                "a planar face is extended across a straight edge; this one strays {off} \
                 from its chord"
            );
        }
    }
    // The outer boundary as a polygon in the chart: its winding says which
    // side of the edge the face is on.
    let mut polygon: Vec<Point2> = Vec::new();
    for e in outer {
        let mut samples = edge_samples(model, e, 8, tol)?;
        if e.orientation() == Orientation::Reversed {
            samples.reverse();
        }
        samples.pop();
        polygon.extend(samples.into_iter().map(chart));
    }
    // Summed about the polygon's first point: the same sum as about the
    // chart's origin, which may stand far off, where it is a difference of
    // huge products.
    let mut area = 0.0;
    if let Some(&anchor) = polygon.first() {
        for (i, p) in polygon.iter().enumerate() {
            let q = polygon[(i + 1) % polygon.len()];
            area += (*p - anchor).cross(q - anchor);
        }
    }
    let (pa, pb) = (chart(a), chart(b));
    let along = (pb - pa).normalized(tol)?;
    // The face lies to the left of a counter-clockwise boundary.
    let outward = if area > 0.0 {
        -along.perpendicular()
    } else {
        along.perpendicular()
    };
    let (qa, qb) = (pa + outward * length, pb + outward * length);
    let ((ua, ub), (va, vb)) = plane.domain();
    let inside = |p: Point2| p.x >= ua && p.x <= ub && p.y >= va && p.y <= vb;
    let surface = if inside(qa) && inside(qb) {
        None
    } else {
        let lo = (ua.min(qa.x).min(qb.x), va.min(qa.y).min(qb.y));
        let hi = (ub.max(qa.x).max(qb.x), vb.max(qa.y).max(qb.y));
        Some(SurfaceGeometry::Plane(PlaneSurface::over(
            plane.plane(),
            (lo.0, hi.0),
            (lo.1, hi.1),
        )?))
    };
    Ok(Plan {
        surface,
        flat: true,
        near: (pa, pb),
        far: (qa, qb),
    })
}

/// An extension along a crossing parameter line: the edge must run along
/// the other parameter, and the face grows across it.
#[allow(clippy::too_many_arguments)]
fn plan_on_chart(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    surface_id: SurfaceId,
    crossed: &Shape,
    length: f64,
    order: Option<usize>,
    tol: Tolerances,
) -> OgeomResult<Plan> {
    let Some(NodeData::Edge(data)) = model.node(crossed).map(|n| n.data()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let (pcurve, range) = match data.pcurve_for(surface_id, crossed.location()) {
        Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
        Some(EdgeRepr::Seam { .. }) => ogeom_bail!(
            Construction,
            "the edge is a seam of the surface; there is no one side to extend the face across"
        ),
        _ => ogeom_bail!(
            Construction,
            "the edge has no trim on the face's surface, so where it runs in the chart is unknown"
        ),
    };
    let Some(pcurve) = model.geometry().pcurve(pcurve) else {
        ogeom_bail!(Dangling, "pcurve is not in this model");
    };
    let samples: Vec<Point2> = (0..=8)
        .map(|k| pcurve.point_at(range.0 + (range.1 - range.0) * f64::from(k) / 8.0, tol))
        .collect::<OgeomResult<_>>()?;
    let (first, last) = (samples[0], samples[8]);
    let (near_from, near_to) = if crossed.orientation() == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    };
    let slack = tol.parametric();
    let still = |pick: fn(Point2) -> f64| {
        samples
            .iter()
            .all(|p| (pick(*p) - pick(first)).abs() <= slack)
    };
    // `across_u`: the edge holds `u` and the face grows in `u`.
    let across_u = if still(|p| p.x) {
        true
    } else if still(|p| p.y) {
        false
    } else {
        ogeom_bail!(
            Construction,
            "the edge does not run along a parameter line of the surface, so there is no \
             crossing parameter to extend the face along"
        );
    };
    let cross = |p: Point2| if across_u { p.x } else { p.y };
    let along = |p: Point2| if across_u { p.y } else { p.x };
    let place = |c: f64, a: f64| {
        if across_u {
            Point2::new(c, a)
        } else {
            Point2::new(a, c)
        }
    };
    let c = cross(first);

    // Which way is out: away from the side the face lies on.
    let Some(sides) = crate::mass_chart::material_sides(model, face, tol) else {
        ogeom_bail!(
            Construction,
            "the face's boundary cannot be walked in its chart, so its side of the edge is unknown"
        );
    };
    let mut inward = 0.0;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (e, p, d) in &sides {
        lo = lo.min(cross(*p));
        hi = hi.max(cross(*p));
        if e.is_same(crossed) {
            inward += cross(Point2::from_vector(*d));
        }
    }
    if inward == 0.0 {
        ogeom_bail!(
            Construction,
            "the face's side of the edge cannot be told from its boundary walk"
        );
    }
    let outward = -inward.signum();
    let middle = 0.5 * (along(near_from) + along(near_to));

    let periodic = if across_u {
        surface.is_periodic_u()
    } else {
        surface.is_periodic_v()
    };
    let ((ua, ub), (va, vb)) = surface.domain();
    let (dlo, dhi) = if across_u { (ua, ub) } else { (va, vb) };

    let (target, widened) = if let SurfaceGeometry::BSpline(patch) = surface {
        let (to, longer) = continue_patch(
            patch,
            across_u,
            outward > 0.0,
            c,
            middle,
            (dlo, dhi),
            length,
            order,
            tol,
        )?;
        (to, Some(SurfaceGeometry::BSpline(longer)))
    } else {
        let straight = matches!(
            surface,
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)
        ) && !across_u;
        if order.is_some() && !straight {
            ogeom_bail!(
                Construction,
                "a linear extension of an analytic surface across a curved parameter line \
                 would leave the surface; extend it naturally"
            );
        }
        let (du, dv) = surface.d1_at(place(c, middle).x, place(c, middle).y, tol)?;
        let speed = if across_u { du } else { dv }.magnitude();
        if speed <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the surface stands still across the edge; there is no direction to extend in"
            );
        }
        let step = length / speed;
        let to = c + outward * step;
        if periodic {
            let period = dhi - dlo;
            if hi - lo + step >= period - slack {
                ogeom_bail!(
                    Construction,
                    "the extension would wrap the face round the periodic surface onto itself"
                );
            }
        }
        (to, widened_analytic(surface, across_u, c, to, slack)?)
    };
    Ok(Plan {
        surface: widened,
        flat: false,
        near: (near_from, near_to),
        far: (
            place(target, along(near_from)),
            place(target, along(near_to)),
        ),
    })
}

/// The analytic surface with its domain widened to hold `to` along the
/// crossing parameter, or `None` where it already does; refused where the
/// surface ends before it.
fn widened_analytic(
    surface: &SurfaceGeometry,
    across_u: bool,
    from: f64,
    to: f64,
    slack: f64,
) -> OgeomResult<Option<SurfaceGeometry>> {
    if across_u {
        // Every analytic chart here repeats in `u`.
        return Ok(None);
    }
    let ((_, _), (va, vb)) = surface.domain();
    let wider = (va.min(to), vb.max(to));
    let holds = to >= va && to <= vb;
    match surface {
        SurfaceGeometry::Cylinder(c) => Ok((!holds).then_some(SurfaceGeometry::Cylinder(
            CylinderSurface::new(c.cylinder(), wider)?,
        ))),
        SurfaceGeometry::Cone(c) => {
            let apex = c.apex_height();
            if (to - apex) * (from - apex) <= 0.0 || (to - apex).abs() <= slack {
                ogeom_bail!(
                    Construction,
                    "the extension would reach the cone's apex, where the far edge collapses"
                );
            }
            Ok((!holds).then_some(SurfaceGeometry::Cone(ConeSurface::new(c.cone(), wider)?)))
        }
        SurfaceGeometry::Sphere(_) => {
            if to.abs() >= std::f64::consts::FRAC_PI_2 - slack {
                ogeom_bail!(
                    Construction,
                    "the extension would reach the sphere's pole, where the far edge collapses"
                );
            }
            Ok(None)
        }
        SurfaceGeometry::Torus(_) => Ok(None),
        _ => ogeom_bail!(Construction, "not an analytic surface"),
    }
}

/// A B-spline patch continued past the side of its domain at `c`, cut back
/// so its new side lies `length` beyond along the crossing parameter line
/// at `middle`; with that new side's parameter.
#[allow(clippy::too_many_arguments)]
fn continue_patch(
    patch: &BSplineSurface,
    across_u: bool,
    at_end: bool,
    c: f64,
    middle: f64,
    (dlo, dhi): (f64, f64),
    length: f64,
    order: Option<usize>,
    tol: Tolerances,
) -> OgeomResult<(f64, BSplineSurface)> {
    let slack = tol.parametric();
    let side = if at_end { dhi } else { dlo };
    if (c - side).abs() > slack {
        ogeom_bail!(
            Construction,
            "the edge is not on a side of the patch's domain; a B-spline face is extended \
             across a side of its patch"
        );
    }
    let closed = if across_u {
        patch.is_closed_u(tol)
    } else {
        patch.is_closed_v(tol)
    };
    if closed {
        ogeom_bail!(
            Construction,
            "the patch closes on itself across the edge; continuing it would overlap the face"
        );
    }
    let degree = if across_u {
        patch.u_knots().degree()
    } else {
        patch.v_knots().degree()
    };
    let order = order.unwrap_or(degree);
    // The continuation is the same polynomial whatever its span: grow it
    // until the crossing line reaches the length, then cut it back there.
    let mut reach = length;
    for _ in 0..12 {
        let longer = patch.extended(across_u, at_end, reach, order, tol)?;
        let ((ua, ub), (va, vb)) = longer.domain();
        let (lo, hi) = if across_u { (ua, ub) } else { (va, vb) };
        let line = Curve::BSpline(if across_u {
            longer.iso_v_curve(middle, tol)?
        } else {
            longer.iso_u_curve(middle, tol)?
        });
        let run = if at_end { (c, hi) } else { (lo, c) };
        let available = crate::length::curve_length(&line, run, tol)?;
        if available >= length {
            let to = if at_end {
                crate::length::parameter_at_length(&line, run, length, tol)?
            } else {
                crate::length::parameter_at_length(&line, run, available - length, tol)?
            };
            let (u, v) = match (across_u, at_end) {
                (true, true) => ((ua, to), (va, vb)),
                (true, false) => ((to, ub), (va, vb)),
                (false, true) => ((ua, ub), (va, to)),
                (false, false) => ((ua, ub), (to, vb)),
            };
            return Ok((to, longer.segment(u, v, tol)?));
        }
        reach *= 2.0;
    }
    ogeom_bail!(
        Construction,
        "the patch's continuation does not reach the length along the crossing line"
    )
}

/// A new edge between two chart points, on the surface: a straight line on
/// a plane, otherwise the parameter line through both. Its trim on
/// `surface_id` is attached and checked against its curve, and it comes
/// back oriented to run `from` to `to`.
fn chart_edge(
    model: &mut Model,
    surface: &SurfaceGeometry,
    surface_id: SurfaceId,
    flat: bool,
    (from, from_vertex): (Point2, &Shape),
    (to, to_vertex): (Point2, &Shape),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (curve, range, pcurve, reversed) = if flat {
        let a = surface.point_at(from.x, from.y, tol)?;
        let b = surface.point_at(to.x, to.y, tol)?;
        let run = from.distance(to);
        let curve: Curve = LineCurve::new(Axis {
            location: a,
            direction: Direction::new(b - a, tol)?,
        })
        .into();
        let pcurve: PlanarCurve = Line2d::over(
            Axis2 {
                location: from,
                direction: Direction2::new(to - from, tol)?,
            },
            0.0,
            run,
        )?
        .into();
        (curve, (0.0, run), pcurve, false)
    } else {
        let fixed_u = (from.x - to.x).abs() <= tol.parametric();
        let (fixed, a, b) = if fixed_u {
            (from.x, from.y, to.y)
        } else {
            (from.y, from.x, to.x)
        };
        let (lo, hi) = (a.min(b), a.max(b));
        let (curve, range, pcurve) = parameter_line(surface, fixed_u, fixed, lo, hi, tol)?;
        (curve, range, pcurve, a > b)
    };
    // The trim follows the curve: checked here, so a wrong chart is an
    // error and never a face whose edges miss its surface.
    for k in 0..=4 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 4.0;
        let on_curve = curve.point_at(t, tol)?;
        let uv = pcurve.point_at(t, tol)?;
        let on_surface = surface.point_at(uv.x, uv.y, tol)?;
        let off = on_curve.distance(on_surface);
        if off > tol.confusion() {
            ogeom_bail!(
                Invariant,
                "a new edge strays {off} from the surface it was built on"
            );
        }
    }
    let (start, end) = if reversed {
        (to_vertex, from_vertex)
    } else {
        (from_vertex, to_vertex)
    };
    let edge = make_edge_between(model, curve, range, start, end, tol)?.shape;
    attach_pcurve(
        model,
        &edge,
        pcurve,
        surface_id,
        Location::identity(),
        range,
    )?;
    if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
        data.assert_same_parameter(true);
    }
    Ok(if reversed { edge.reversed() } else { edge })
}

/// The parameter line of `surface` at `fixed` (of `u` where `fixed_u`, else
/// of `v`) from `lo` to `hi` of the other parameter: its curve, the curve's
/// range, and its trim.
fn parameter_line(
    surface: &SurfaceGeometry,
    fixed_u: bool,
    fixed: f64,
    lo: f64,
    hi: f64,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64), PlanarCurve)> {
    let place = |p: f64| {
        if fixed_u {
            Point2::new(fixed, p)
        } else {
            Point2::new(p, fixed)
        }
    };
    let heading = if fixed_u {
        Direction2::new(Vector2::new(0.0, 1.0), tol)?
    } else {
        Direction2::new(Vector2::new(1.0, 0.0), tol)?
    };
    if let SurfaceGeometry::BSpline(patch) = surface {
        // The patch's own iso-curve, parameterized as the chart is.
        let curve = Curve::BSpline(if fixed_u {
            patch.iso_u_curve(fixed, tol)?
        } else {
            patch.iso_v_curve(fixed, tol)?
        });
        let pcurve: PlanarCurve = Line2d::over(
            Axis2 {
                location: place(0.0),
                direction: heading,
            },
            lo,
            hi,
        )?
        .into();
        return Ok((curve, (lo, hi), pcurve));
    }
    let start_uv = place(lo);
    let start = surface.point_at(start_uv.x, start_uv.y, tol)?;
    let (du, dv) = surface.d1_at(start_uv.x, start_uv.y, tol)?;
    let tangent = if fixed_u { dv } else { du };
    let speed = tangent.magnitude();
    let ruling = matches!(
        surface,
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)
    ) && fixed_u;
    // `rate`: the curve's parameter per unit of the chart's.
    let (curve, rate): (Curve, f64) = if ruling {
        let line = LineCurve::new(Axis {
            location: start,
            direction: Direction::new(tangent, tol)?,
        });
        (line.into(), speed)
    } else {
        let centre = match surface {
            SurfaceGeometry::Sphere(s) if fixed_u => s.sphere().centre(),
            SurfaceGeometry::Torus(t) if fixed_u => {
                let torus = t.torus();
                let f = torus.frame();
                let out = f.x().vector() * fixed.cos() + f.y().vector() * fixed.sin();
                f.origin() + out * torus.major_radius()
            }
            _ => {
                // A parallel: centred on the axis at the start's height.
                let f = axis_frame(surface)?;
                let z = f.z().vector();
                f.origin() + z * (start - f.origin()).dot(z)
            }
        };
        let radial = start - centre;
        let radius = radial.magnitude();
        if radius <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the parameter line collapses to a point; no edge runs along it"
            );
        }
        let frame = Frame::new(
            centre,
            Direction::from_cross(radial, tangent, tol)?,
            Direction::new(radial, tol)?,
            tol,
        )?;
        let circle = CircleCurve::new(Circle::new(frame, radius, tol)?);
        (circle.into(), speed / radius)
    };
    let end = rate * (hi - lo);
    let pcurve: PlanarCurve = if (rate - 1.0).abs() <= f64::EPSILON * 16.0 {
        Line2d::over(
            Axis2 {
                location: start_uv,
                direction: heading,
            },
            0.0,
            end,
        )?
        .into()
    } else {
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, end, end], 1)?;
        ogeom_geom::BSpline2d::new(knots, vec![start_uv, place(hi)], tol)?.into()
    };
    Ok((curve, (0.0, end), pcurve))
}

/// The frame whose `z` is a surface of revolution's axis.
fn axis_frame(surface: &SurfaceGeometry) -> OgeomResult<Frame> {
    Ok(match surface {
        SurfaceGeometry::Cylinder(c) => c.cylinder().frame(),
        SurfaceGeometry::Cone(c) => c.cone().frame(),
        SurfaceGeometry::Sphere(s) => s.sphere().frame(),
        SurfaceGeometry::Torus(t) => t.torus().frame(),
        _ => ogeom_bail!(Construction, "the surface has no axis"),
    })
}

/// Give each edge the trims it has on `from` on `to` as well: the surface
/// grown keeps its parameters, so a trim on the old surface is one on the
/// new.
fn retrim(model: &mut Model, edges: &[Shape], from: SurfaceId, to: SurfaceId) -> OgeomResult<()> {
    let mut done: Vec<Shape> = Vec::new();
    for edge in edges {
        if done.iter().any(|d| d.node() == edge.node()) {
            continue;
        }
        done.push(edge.clone());
        let Some(NodeData::Edge(data)) = model.node_mut(edge).map(|n| n.data_mut()) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let agreed = data.same_parameter();
        let copies: Vec<EdgeRepr> = data
            .representations
            .iter()
            .filter_map(|r| match r {
                EdgeRepr::PCurve {
                    curve,
                    surface,
                    location,
                    range,
                } if *surface == from => Some(EdgeRepr::PCurve {
                    curve: *curve,
                    surface: to,
                    location: location.clone(),
                    range: *range,
                }),
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    surface,
                    location,
                    range,
                } if *surface == from => Some(EdgeRepr::Seam {
                    forward: *forward,
                    reversed: *reversed,
                    surface: to,
                    location: location.clone(),
                    range: *range,
                }),
                _ => None,
            })
            .collect();
        for repr in copies {
            data.add(repr);
        }
        data.assert_same_parameter(agreed);
    }
    Ok(())
}

/// A vertex's point, placed.
fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Dangling, "vertex is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// `count + 1` points along an edge's curve, placed, from its range's start
/// to its end.
fn edge_samples(
    model: &Model,
    edge: &Shape,
    count: u32,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "an edge of the face has no curve in space");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placement = edge.transform(model.datums())?;
    (0..=count)
        .map(|k| {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(count);
            Ok(placement.apply(geometry.point_at(t, tol)?))
        })
        .collect()
}
