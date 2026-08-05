//! The upgrade family: same-domain face unification, collinear edge
//! merging, and tolerance reduction — undoing the splits an operation left
//! behind without changing the shape they describe.

use std::collections::{HashMap, HashSet};

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve2d as _, Curve3d as _, PlanarCurve, SurfaceGeometry};
use ogeom_math::{Direction2, Plane, Point2, Transform};
use ogeom_topo::{
    EdgeRepr, Location, Model, NodeData, PCurveId, Shape, ShapeType, SurfaceId, TShapeId,
    explore_unique,
};

use crate::reshape::Reshape;

/// Merge adjacent faces lying on one carrier into single faces.
///
/// Two faces qualify when they share an edge and their surfaces are the
/// same plane, in the same parameterization — the split a boolean or an
/// exchange leaves behind. The merged face keeps the first face's surface;
/// the shared edges dissolve; the remaining boundary re-chains into wires.
/// Faces on other carriers pass through untouched — a curved unification
/// wants parameterization transport this deliberately does not guess at.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// merged boundary fails to chain into closed wires.
pub fn unify_same_domain(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    let faces = explore_unique(model, shape, ShapeType::Face)?;

    // The carrier of each planar face: where it lies in the world, which
    // decides the grouping, and the chart it was stored in, which is where
    // any synthesized pcurve has to land.
    let mut carriers: Vec<Option<Carrier>> = Vec::with_capacity(faces.len());
    for face in &faces {
        carriers.push(carrier_of(model, face, tol)?);
    }

    // Union-find over faces sharing an edge on one carrier.
    let mut group: Vec<usize> = (0..faces.len()).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    let mut edge_users: HashMap<TShapeId, Vec<usize>> = HashMap::new();
    for (i, face) in faces.iter().enumerate() {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            edge_users.entry(edge.node()).or_default().push(i);
        }
    }
    for users in edge_users.values() {
        for pair in users.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            // One carrier means one *plane*, not one frame: two halves of a
            // wall carry their own origins on the same flat.
            if let (Some(ca), Some(cb)) = (carriers[a].as_ref(), carriers[b].as_ref())
                && ca.world.normal().dot(cb.world.normal()) > 1.0 - tol.angular()
                && ca.world.signed_distance_to(cb.world.frame().origin()).abs() <= tol.confusion()
            {
                let (ra, rb) = (root(&mut group, a), root(&mut group, b));
                group[ra] = rb;
            }
        }
    }

    let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..faces.len() {
        let r = root(&mut group, i);
        clusters.entry(r).or_default().push(i);
    }

    let mut reshape = Reshape::new();
    for members in clusters.values() {
        if members.len() < 2 {
            continue;
        }
        // Interior edges — used by two members — dissolve; the rest chain.
        let mut counts: HashMap<TShapeId, (usize, Shape)> = HashMap::new();
        for &i in members {
            for edge in explore_unique(model, &faces[i], ShapeType::Edge)? {
                counts
                    .entry(edge.node())
                    .and_modify(|(n, _)| *n += 1)
                    .or_insert((1, edge));
            }
        }
        let boundary: Vec<Shape> = {
            let mut edges: Vec<(TShapeId, Shape)> = counts
                .into_iter()
                .filter(|(_, (n, _))| *n == 1)
                .map(|(id, (_, e))| (id, e))
                .collect();
            edges.sort_by_key(|(id, _)| *id);
            edges.into_iter().map(|(_, e)| e).collect()
        };
        if boundary.len() < 3 {
            continue;
        }
        let keeper = members[0];
        let surface_id = {
            let Some(NodeData::Face(data)) = model.node(&faces[keeper]).map(|n| n.data().clone())
            else {
                continue;
            };
            data.surface
        };
        // Every boundary edge needs a pcurve for the kept surface; a plane's
        // is direct projection.
        let Some(carrier) = carriers[keeper].clone() else {
            continue;
        };
        for edge in &boundary {
            ensure_planar_pcurve(model, edge, surface_id, &carrier, tol)?;
        }
        let wire = ogeom_algo::make_wire_unordered(model, &boundary, tol)?.shape;
        let merged = {
            let Some(NodeData::Face(data)) = model.node(&faces[keeper]).map(|n| n.data().clone())
            else {
                continue;
            };
            model.add_face(*data, std::slice::from_ref(&wire))?
        };
        reshape.replace(&faces[keeper], merged);
        for &other in &members[1..] {
            reshape.remove(&faces[other]);
        }
    }
    if reshape.is_empty() {
        let mut built = Built::from_nothing(shape.clone());
        built.history.modify(shape, shape.clone());
        return Ok(built);
    }
    reshape.apply(model, shape)
}

/// Where a planar face lies and which chart it was stored in.
#[derive(Debug, Clone)]
struct Carrier {
    /// The plane in world space — what "same carrier" is decided on.
    world: Plane,
    /// The plane as the surface stores it: the chart pcurves speak in.
    stored: Plane,
    /// What takes the stored chart to the world.
    placement: Transform,
}

/// The carrier of a planar face, if it is one.
fn carrier_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Carrier>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        return Ok(None);
    };
    let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let stored = p.plane();
    let placement = face.transform(model.datums())?;
    if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
        return Ok(None);
    }
    Ok(Some(Carrier {
        world: stored.transformed(&placement, tol)?,
        stored,
        placement,
    }))
}

/// Attach a pcurve for `surface_id` to `edge` if it does not carry one —
/// direct projection into the stored plane's own frame, exact.
///
/// The chart is the *stored* plane's, so the world point is carried back
/// through the face's placement before it is flattened; a pcurve read under
/// that placement then lands where the edge is.
fn ensure_planar_pcurve(
    model: &mut Model,
    edge: &Shape,
    surface_id: SurfaceId,
    carrier: &Carrier,
    tol: Tolerances,
) -> OgeomResult<()> {
    let (curve, range, has) = {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge holds no edge data");
        };
        let has = data.pcurve_for(surface_id, edge.location()).is_some();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "an edge with no curve cannot be unified over");
        };
        (*curve, *range, has)
    };
    if has {
        return Ok(());
    }
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        ogeom_bail!(Construction, "edge refers to a curve not in this model");
    };
    let back = carrier.placement.inverse()?;
    let flat = |p: ogeom_math::Point| {
        let local = carrier.stored.frame().to_local(back.apply(p));
        Point2::new(local.x, local.y)
    };
    let a = flat(geometry.point_at(range.0, tol)?);
    let b = flat(geometry.point_at(range.1, tol)?);
    ogeom_algo::attach_pcurve(
        model,
        edge,
        chart_line(a, b, range, tol)?,
        surface_id,
        edge.location().clone(),
        range,
    )
}

/// The chart line running from `a` at `range.0` to `b` at `range.1`.
///
/// A [`Line2d`](ogeom_geom::Line2d) reads its parameter from the axis
/// origin outward, so the origin is where the parameter would have been
/// zero — not where the range starts.
fn chart_line(
    a: Point2,
    b: Point2,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let span = range.1 - range.0;
    let direction = Direction2::new((b - a) / span, tol)?;
    let origin = a - direction.vector() * range.0;
    Ok(
        ogeom_geom::Line2d::over(ogeom_math::Axis2::new(origin, direction), range.0, range.1)?
            .into(),
    )
}

/// Merge chains of edges lying on one curve into single edges.
///
/// Within each wire, consecutive edges continuing one curve — the same
/// stored curve over contiguous ranges, or two collinear lines meeting end
/// to end — become one edge over the joined range, and the vertex between
/// them dissolves. Every pcurve the pair carried is joined with it and then
/// *measured* against the joined curve: a pair whose parameterizations do
/// not join cleanly is left split rather than merged into a face that could
/// no longer be triangulated.
///
/// # Errors
///
/// As the model's own builders.
pub fn merge_edges(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    // One pass merges disjoint pairs; a chain of three needs the next pass
    // to see the pair the first one made. Eight is well past any real chain.
    let mut current = shape.clone();
    let mut steps: Vec<History> = Vec::new();
    for _ in 0..8 {
        let Some(step) = merge_pass(model, &current, tol)? else {
            break;
        };
        current = step.shape;
        steps.push(step.history);
    }
    if steps.is_empty() {
        let mut built = Built::from_nothing(shape.clone());
        built.history.modify(shape, shape.clone());
        return Ok(built);
    }
    Ok(Built::new(current, History::chain(&steps)))
}

/// One merge pass: every disjoint joinable pair, in one rebuild. `None`
/// when nothing joins.
fn merge_pass(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Option<Built>> {
    // Which parameterizations each edge is actually *read* on: the faces
    // bounding it, not whatever reprs it accumulated along the way. An
    // earlier unification leaves pcurves on surfaces nothing carries any
    // more, and a join that tried to honour those would refuse work it can
    // do.
    let mut needed: HashMap<TShapeId, Vec<(SurfaceId, Location)>> = HashMap::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
            continue;
        };
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let slot = (data.surface, edge.location().clone());
            let slots = needed.entry(edge.node()).or_default();
            if !slots.contains(&slot) {
                slots.push(slot);
            }
        }
    }

    let wires = explore_unique(model, shape, ShapeType::Wire)?;
    let mut reshape = Reshape::new();
    let mut merged_nodes: HashSet<TShapeId> = HashSet::new();
    for wire in &wires {
        let edges = model.ordered_children_of(wire)?;
        if edges.len() < 2 {
            continue;
        }
        // A closed wire has no last edge: its end is the start again, and
        // the split sitting across that join is a split like any other.
        let closed = is_closed_wire(model, &edges, tol);
        let pairs = if closed { edges.len() } else { edges.len() - 1 };
        let mut i = 0;
        while i < pairs {
            let (a, b) = (edges[i].clone(), edges[(i + 1) % edges.len()].clone());
            if a.node() == b.node()
                || merged_nodes.contains(&a.node())
                || merged_nodes.contains(&b.node())
            {
                i += 1;
                continue;
            }
            let mut slots: Vec<(SurfaceId, Location)> =
                needed.get(&a.node()).cloned().unwrap_or_default();
            for slot in needed.get(&b.node()).into_iter().flatten() {
                if !slots.contains(slot) {
                    slots.push(slot.clone());
                }
            }
            let Some(join) = joinable(model, &a, &b, &slots, tol)? else {
                i += 1;
                continue;
            };
            let joined = ogeom_algo::make_edge(model, join.curve, join.range, tol)?.shape;
            for pcurve in join.pcurves {
                ogeom_algo::attach_pcurve(
                    model,
                    &joined,
                    pcurve.curve,
                    pcurve.surface,
                    pcurve.location,
                    pcurve.range,
                )?;
            }
            merged_nodes.insert(a.node());
            merged_nodes.insert(b.node());
            reshape.replace(&a, joined);
            reshape.remove(&b);
            i += 2;
        }
    }
    if reshape.is_empty() {
        return Ok(None);
    }
    reshape.apply(model, shape).map(Some)
}

/// Whether an ordered edge list closes back on its own start.
fn is_closed_wire(model: &Model, edges: &[Shape], tol: Tolerances) -> bool {
    let ends = |e: &Shape| -> Option<(ogeom_math::Point, ogeom_math::Point)> {
        let vertices = explore_unique(model, e, ShapeType::Vertex).ok()?;
        let first = model.node(vertices.first()?)?.data().as_vertex()?.point;
        let last = model.node(vertices.last()?)?.data().as_vertex()?.point;
        Some((first, last))
    };
    let (Some(first), Some(last)) = (ends(&edges[0]), ends(&edges[edges.len() - 1])) else {
        return false;
    };
    [first.0, first.1]
        .iter()
        .any(|p| p.distance(last.0) <= tol.confusion() || p.distance(last.1) <= tol.confusion())
}

/// One joined edge, ready to build.
struct Join {
    /// The curve the joined edge runs on.
    curve: ogeom_geom::Curve,
    /// The range it runs over.
    range: (f64, f64),
    /// The joined pcurve for every surface the pair was parameterized on.
    pcurves: Vec<PcurveJoin>,
}

/// A joined pcurve, in the slot it goes back into.
struct PcurveJoin {
    surface: SurfaceId,
    location: Location,
    curve: PlanarCurve,
    range: (f64, f64),
}

/// The three-dimensional half of a join: the curve, the range it runs
/// over, and whether `a` runs first along it.
type CurveJoin = (ogeom_geom::Curve, (f64, f64), bool);

/// The three-dimensional half of a join: the curve, its range, and whether
/// `a` runs first along it.
fn joinable_curve(
    model: &Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<CurveJoin>> {
    let read = |e: &Shape| -> Option<(ogeom_topo::CurveId, (f64, f64))> {
        let data = model.node(e)?.data().as_edge()?;
        let EdgeRepr::Curve3d { curve, range, .. } = data.curve3d()? else {
            return None;
        };
        Some((*curve, *range))
    };
    let (Some((ca, ra)), Some((cb, rb))) = (read(a), read(b)) else {
        return Ok(None);
    };
    if ca == cb {
        let Some(geometry) = model.geometry().curve(ca).cloned() else {
            return Ok(None);
        };
        if (ra.1 - rb.0).abs() <= tol.parametric() {
            return Ok(Some((geometry, (ra.0, rb.1), true)));
        }
        if (rb.1 - ra.0).abs() <= tol.parametric() {
            return Ok(Some((geometry, (rb.0, ra.1), false)));
        }
        return Ok(None);
    }
    // Distinct curves: collinear lines meeting end to end still join.
    let (Some(ga), Some(gb)) = (
        model.geometry().curve(ca).cloned(),
        model.geometry().curve(cb).cloned(),
    ) else {
        return Ok(None);
    };
    let (ogeom_geom::Curve::Line(la), ogeom_geom::Curve::Line(lb)) = (&ga, &gb) else {
        return Ok(None);
    };
    if !la.axis().is_collinear(lb.axis(), tol) {
        return Ok(None);
    }
    let (a0, a1) = (ga.point_at(ra.0, tol)?, ga.point_at(ra.1, tol)?);
    let (b0, b1) = (gb.point_at(rb.0, tol)?, gb.point_at(rb.1, tol)?);
    let fresh = |from: ogeom_math::Point,
                 to: ogeom_math::Point,
                 a_first: bool|
     -> OgeomResult<Option<CurveJoin>> {
        let segment = ogeom_geom::LineCurve::segment(from, to, tol)?;
        let range = segment.domain();
        Ok(Some((segment.into(), range, a_first)))
    };
    if a1.distance(b0) <= tol.confusion() {
        return fresh(a0, b1, true);
    }
    if b1.distance(a0) <= tol.confusion() {
        return fresh(b0, a1, false);
    }
    Ok(None)
}

/// Whether two consecutive edges continue one curve, pcurves and all.
///
/// `slots` is what the joined edge will be read on: every one of them has
/// to come out of the join, or the merge would leave a face that cannot be
/// triangulated — worse than the split it came to fix.
fn joinable(
    model: &Model,
    a: &Shape,
    b: &Shape,
    slots: &[(SurfaceId, Location)],
    tol: Tolerances,
) -> OgeomResult<Option<Join>> {
    let Some((curve, range, a_first)) = joinable_curve(model, a, b, tol)? else {
        return Ok(None);
    };
    let (first, second) = if a_first { (a, b) } else { (b, a) };
    let reprs = |e: &Shape| -> Vec<EdgeRepr> {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .map(|d| d.representations.to_vec())
            .unwrap_or_default()
    };
    let (first_reprs, second_reprs) = (reprs(first), reprs(second));
    // A seam is a face's own doubling, not a split; leave those pairs be.
    if first_reprs
        .iter()
        .chain(&second_reprs)
        .any(|r| matches!(r, EdgeRepr::Seam { .. }))
    {
        return Ok(None);
    }
    let parametric = |reprs: &[EdgeRepr]| -> Vec<(SurfaceId, Location, PCurveId, (f64, f64))> {
        reprs
            .iter()
            .filter_map(|r| match r {
                EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    location,
                } => Some((*surface, location.clone(), *curve, *range)),
                _ => None,
            })
            .collect()
    };
    let (ones, twos) = (parametric(&first_reprs), parametric(&second_reprs));
    let mut pcurves = Vec::with_capacity(slots.len());
    for (surface, location) in slots {
        let (surface, location) = (*surface, location.clone());
        let find = |from: &[(SurfaceId, Location, PCurveId, (f64, f64))]| {
            from.iter()
                .find(|(s, l, _, _)| *s == surface && *l == location)
                .cloned()
        };
        let (Some((_, _, id_one, range_one)), Some((_, _, id_two, range_two))) =
            (find(&ones), find(&twos))
        else {
            return Ok(None);
        };
        let (Some(one), Some(two)) = (
            model.geometry().pcurve(id_one).cloned(),
            model.geometry().pcurve(id_two).cloned(),
        ) else {
            return Ok(None);
        };
        let joined = if id_one == id_two && (range_one.1 - range_two.0).abs() <= tol.parametric() {
            (one, (range_one.0, range_two.1))
        } else {
            // Two chart lines meeting end to end join over the 3D range the
            // joined edge will be read at.
            let (PlanarCurve::Line(_), PlanarCurve::Line(_)) = (&one, &two) else {
                return Ok(None);
            };
            let start = one.point_at(range_one.0, tol)?;
            let end = two.point_at(range_two.1, tol)?;
            if one
                .point_at(range_one.1, tol)?
                .distance(two.point_at(range_two.0, tol)?)
                > tol.confusion()
            {
                return Ok(None);
            }
            (chart_line(start, end, range, tol)?, range)
        };
        pcurves.push(PcurveJoin {
            surface,
            location,
            curve: joined.0,
            range: joined.1,
        });
    }
    // The plan is only a plan until it measures up against the joined curve.
    for pcurve in &pcurves {
        if !agrees(model, &curve, range, pcurve, tol)? {
            return Ok(None);
        }
    }
    Ok(Some(Join {
        curve,
        range,
        pcurves,
    }))
}

/// Whether a joined pcurve, read the way a face reads it — the same
/// fraction along both ranges — lands on the joined curve.
fn agrees(
    model: &Model,
    curve: &ogeom_geom::Curve,
    range: (f64, f64),
    pcurve: &PcurveJoin,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    let Some(surface) = model.geometry().surface(pcurve.surface) else {
        return Ok(false);
    };
    let placement = pcurve.location.composed(model.datums())?;
    for k in 0..=6 {
        let f = f64::from(k) / 6.0;
        let t = (range.1 - range.0).mul_add(f, range.0);
        let pt = (pcurve.range.1 - pcurve.range.0).mul_add(f, pcurve.range.0);
        let (Ok(on_curve), Ok(chart)) = (curve.point_at(t, tol), pcurve.curve.point_at(pt, tol))
        else {
            return Ok(false);
        };
        let Ok(lifted) = surface.point_at(chart.x, chart.y, tol) else {
            return Ok(false);
        };
        if placement.apply(lifted).distance(on_curve) > tol.confusion() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Shrink every edge and vertex tolerance to what the geometry measures,
/// floored at the confusion tolerance — the inverse of the widenings repair
/// operations apply, safe because it is measured the same way.
///
/// Returns how many claims shrank.
///
/// # Errors
///
/// As evaluation.
pub fn reduce_tolerances(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<usize> {
    let mut shrunk = 0;
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let measured = {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                continue;
            };
            // The claim an edge tolerance covers: its pcurves against its
            // curve at matched parameters.
            let mut worst = 0.0f64;
            let mut measurable = false;
            for representation in &data.representations {
                let EdgeRepr::PCurve {
                    curve: pc,
                    range: prange,
                    surface,
                    location,
                } = representation
                else {
                    continue;
                };
                let Some(pcurve) = model.geometry().pcurve(*pc) else {
                    continue;
                };
                let Some(surface_geometry) = model.geometry().surface(*surface) else {
                    continue;
                };
                use ogeom_geom::Surface as _;
                for k in 0..=8 {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / 8.0;
                    let pt = prange.0 + (prange.1 - prange.0) * f64::from(k) / 8.0;
                    let (Ok(on_curve), Ok(chart)) =
                        (geometry.point_at(t, tol), pcurve.point_at(pt, tol))
                    else {
                        continue;
                    };
                    let Ok(lifted) = surface_geometry.point_at(chart.x, chart.y, tol) else {
                        continue;
                    };
                    let Ok(placement) = location.composed(model.datums()) else {
                        continue;
                    };
                    worst = worst.max(placement.apply(lifted).distance(on_curve));
                    measurable = true;
                }
            }
            measurable.then_some(worst)
        };
        let Some(worst) = measured else { continue };
        let target = ogeom_core::Tolerance::new(worst + tol.confusion())?;
        if let Some(node) = model.node_mut(&edge)
            && let NodeData::Edge(data) = node.data_mut()
            && target.get() < data.tolerance.get()
        {
            data.tolerance = target;
            shrunk += 1;
        }
    }
    Ok(shrunk)
}
