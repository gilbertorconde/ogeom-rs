//! Unification on curved carriers: cylinders, cones, spheres and tori on
//! one axis with one size, and extrusions of one basis curve along one
//! line.
//!
//! Two charts of one such carrier differ by a map that is affine in each
//! parameter alone: the angle about the axis turned and possibly run the
//! other way, the height (or latitude, or tube angle, or sweep) shifted and
//! possibly run the other way. The map is read off the two frames exactly
//! and then measured over the face it carries; a pair the frames do not
//! relate so is not one carrier here, however close the surfaces pass.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction2, Frame, Point, Point2, Transform, Transform2, Vector, Vector2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, SurfaceId, TShapeId};

/// How many points each pcurve is read at, ends included.
const SAMPLES: u32 = 17;

/// Where a curved face lies: its surface as stored, what places that in
/// the world, the side the face takes, and the box its boundary spans in
/// the stored chart.
#[derive(Debug, Clone)]
pub(super) struct Curved {
    pub(super) id: SurfaceId,
    surface: SurfaceGeometry,
    placement: Transform,
    reversed: bool,
    window: (Point2, Point2),
}

/// The carrier of a face on a cylinder, cone, sphere, torus or extrusion,
/// or `None` for any other face, a face placed by a scaling, or one whose
/// boundary is not all plain pcurves on its surface (a seam closes the face
/// round its carrier, and is never merged across).
pub(super) fn carrier_of(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Curved>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        return Ok(None);
    };
    if !matches!(
        surface,
        SurfaceGeometry::Cylinder(_)
            | SurfaceGeometry::Cone(_)
            | SurfaceGeometry::Sphere(_)
            | SurfaceGeometry::Torus(_)
            | SurfaceGeometry::Extrusion(_)
    ) {
        return Ok(None);
    }
    let placement = face.transform(model.datums())?;
    if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
        return Ok(None);
    }
    let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for wire in model.children_of(face)? {
        for edge in model.ordered_children_of(&wire)? {
            let Some((pcurve, range)) = pcurve_on(model, &edge, data.surface) else {
                return Ok(None);
            };
            for p in samples(&pcurve, range, false, tol)? {
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
    }
    if !(lo.x.is_finite() && lo.y.is_finite() && hi.x.is_finite() && hi.y.is_finite()) {
        return Ok(None);
    }
    Ok(Some(Curved {
        id: data.surface,
        surface,
        placement,
        reversed: face.orientation() == Orientation::Reversed,
        window: (lo, hi),
    }))
}

/// Whether two curved carriers are one surface facing one way.
pub(super) fn same_carrier(a: &Curved, b: &Curved, tol: Tolerances) -> bool {
    chart_map(a, b, tol).is_some_and(|map| faces_one_way(a, b, &map))
}

/// Whether two faces related by `map` take the same side of their
/// carrier. Each face's outward normal is its placement applied to its
/// stored chart's `Su × Sv`, turned back where the face is used reversed;
/// in placed tangents that is `Su × Sv` turned again by a mirroring
/// placement, and the placed tangents of one chart are those of the other
/// times the map's Jacobian, whose determinant is the sign of its flips.
fn faces_one_way(a: &Curved, b: &Curved, map: &ChartMap) -> bool {
    let turned = |c: &Curved| c.reversed != !c.placement.preserves_handedness();
    (turned(a) != (map.flip.0 != map.flip.1)) == turned(b)
}

/// A map between two charts of one carrier: each parameter run one way or
/// the other, then shifted.
#[derive(Debug, Clone, Copy)]
struct ChartMap {
    flip: (bool, bool),
    shift: Vector2,
}

impl ChartMap {
    const IDENTITY: Self = Self {
        flip: (false, false),
        shift: Vector2::ZERO,
    };

    fn apply(&self, p: Point2) -> Point2 {
        Point2::new(
            if self.flip.0 { -p.x } else { p.x } + self.shift.x,
            if self.flip.1 { -p.y } else { p.y } + self.shift.y,
        )
    }

    /// The same map as a planar transform, built from exact reflections in
    /// the chart's axes so a pcurve carried by it lands where `apply` puts
    /// its points.
    fn transform(&self) -> Transform2 {
        let mut map = Transform2::translation(self.shift);
        if self.flip.0 {
            map = map * Transform2::line_mirror(Point2::ORIGIN, Direction2::X);
        }
        if self.flip.1 {
            map = map * Transform2::line_mirror(Point2::ORIGIN, Direction2::Y);
        }
        map
    }
}

/// A stored frame's origin and axes carried into the world, as vectors: a
/// mirroring placement leaves them left-handed, which the chart formulas
/// read the same way.
struct Axes {
    origin: Point,
    x: Vector,
    y: Vector,
    z: Vector,
}

fn axes(frame: &Frame, placement: &Transform) -> Axes {
    Axes {
        origin: placement.apply(frame.origin()),
        x: placement.apply_vector(frame.x().vector()),
        y: placement.apply_vector(frame.y().vector()),
        z: placement.apply_vector(frame.z().vector()),
    }
}

/// The angle map between two frames sharing an axis line, for charts
/// whose angle `u` names the direction `cos u·x + sin u·y`: `from`'s angle
/// `u` is `to`'s `±u + turn`, run backwards where the two frames see the
/// turn about the axis opposite ways. Returns that flip, the turn, and
/// whether the axes run opposite ways. `None` where the axes are not
/// parallel.
fn about_one_axis(from: &Axes, to: &Axes, tol: Tolerances) -> Option<(bool, f64, bool)> {
    let along = from.z.dot(to.z);
    if (along.abs() - 1.0).abs() > tol.angular() {
        return None;
    }
    let turn = from.x.dot(to.y).atan2(from.x.dot(to.x));
    let backwards = from.x.cross(from.y).dot(to.x.cross(to.y)) < 0.0;
    Some((backwards, turn, along < 0.0))
}

/// How far `offset` runs along `axis`, where it runs along it alone.
fn along_axis(offset: Vector, axis: Vector, tol: Tolerances) -> Option<f64> {
    let run = offset.dot(axis);
    ((offset - axis * run).magnitude() <= tol.confusion()).then_some(run)
}

/// The map taking `from`'s chart onto `to`'s, where the two are one
/// carrier, and `None` otherwise.
fn chart_map(from: &Curved, to: &Curved, tol: Tolerances) -> Option<ChartMap> {
    use SurfaceGeometry as S;
    let close = |a: f64, b: f64| (a - b).abs() <= tol.confusion();
    let map = match (&from.surface, &to.surface) {
        (S::Cylinder(a), S::Cylinder(b)) => {
            let (a, b) = (a.cylinder(), b.cylinder());
            if !close(a.radius(), b.radius()) {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            let rise = along_axis(fa.origin - fb.origin, fb.z, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, rise),
            }
        }
        (S::Cone(a), S::Cone(b)) => {
            let (a, b) = (a.cone(), b.cone());
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            let rise = along_axis(fa.origin - fb.origin, fb.z, tol)?;
            // The radius `r + h·tan(a)` at `from`'s height `h` is the
            // radius at `to`'s height `±h + rise`.
            let (slope_a, slope_b) = (a.half_angle().tan(), b.half_angle().tan());
            let slope_b_seen = if opposed { -slope_b } else { slope_b };
            if (slope_a - slope_b_seen).abs() > tol.angular()
                || !close(
                    a.reference_radius(),
                    rise.mul_add(slope_b, b.reference_radius()),
                )
            {
                return None;
            }
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, rise),
            }
        }
        (S::Sphere(a), S::Sphere(b)) => {
            let (a, b) = (a.sphere(), b.sphere());
            if !close(a.radius(), b.radius()) {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            if fa.origin.distance(fb.origin) > tol.confusion() {
                return None;
            }
            // Latitude is the angle off the equator, so it runs the other
            // way where the axes do.
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, 0.0),
            }
        }
        (S::Torus(a), S::Torus(b)) => {
            let (a, b) = (a.torus(), b.torus());
            if !close(a.major_radius(), b.major_radius())
                || !close(a.minor_radius(), b.minor_radius())
            {
                return None;
            }
            let (fa, fb) = (
                axes(&a.frame(), &from.placement),
                axes(&b.frame(), &to.placement),
            );
            if fa.origin.distance(fb.origin) > tol.confusion() {
                return None;
            }
            // The tube angle rises along the axis, so it too runs the other
            // way where the axes do.
            let (backwards, turn, opposed) = about_one_axis(&fa, &fb, tol)?;
            ChartMap {
                flip: (backwards, opposed),
                shift: Vector2::new(turn, 0.0),
            }
        }
        (S::Extrusion(a), S::Extrusion(b)) => {
            use ogeom_geom::Transformable as _;
            let basis_a = a.curve().transformed(&from.placement, tol).ok()?;
            let basis_b = b.curve().transformed(&to.placement, tol).ok()?;
            let da = from.placement.apply_vector(a.direction().vector());
            let db = to.placement.apply_vector(b.direction().vector());
            let along = da.dot(db);
            if (along.abs() - 1.0).abs() > tol.angular() {
                return None;
            }
            // `b`'s basis is `a`'s moved `lift` along `da` at the same
            // parameter, so `a`'s `(u, v)` is `b`'s `(u, ±(v - lift))`.
            let lift = translate_of(&basis_a, &basis_b, da, tol)?;
            ChartMap {
                flip: (false, along < 0.0),
                shift: Vector2::new(0.0, if along < 0.0 { lift } else { -lift }),
            }
        }
        _ => return None,
    };
    // The map is read off the frames; it must also carry every point of
    // the face where the other chart has it.
    let (lo, hi) = from.window;
    const GRID: u32 = 4;
    for i in 0..=GRID {
        for j in 0..=GRID {
            let (fu, fv) = (
                f64::from(i) / f64::from(GRID),
                f64::from(j) / f64::from(GRID),
            );
            let at = Point2::new(
                (hi.x - lo.x).mul_add(fu, lo.x),
                (hi.y - lo.y).mul_add(fv, lo.y),
            );
            let (Some(p), Some(q)) = (
                world_point(from, at, tol),
                world_point(to, map.apply(at), tol),
            ) else {
                return None;
            };
            if p.distance(q) > tol.confusion() {
                return None;
            }
        }
    }
    Some(map)
}

/// The point at `(u, v)` of a carrier's chart, placed in the world, read
/// off the unbounded surface so a chart point past the stored window still
/// names its point.
fn world_point(c: &Curved, at: Point2, tol: Tolerances) -> Option<Point> {
    use ogeom_math::elementary as e;
    let local = match &c.surface {
        SurfaceGeometry::Cylinder(s) => e::cylinder_at(&s.cylinder(), at.x, at.y).point,
        SurfaceGeometry::Cone(s) => e::cone_at(&s.cone(), at.x, at.y).point,
        SurfaceGeometry::Sphere(s) => e::sphere_at(&s.sphere(), at.x, at.y).point,
        SurfaceGeometry::Torus(s) => e::torus_at(&s.torus(), at.x, at.y).point,
        SurfaceGeometry::Extrusion(s) => {
            s.curve().point_at(at.x, tol).ok()? + s.direction().vector() * at.y
        }
        _ => return None,
    };
    Some(c.placement.apply(local))
}

/// How far along `d` the curve `b` is `a` moved, where it is `a` moved
/// along `d` with its parameter kept: a line, circle or B-spline whose
/// defining points all moved by one multiple of `d` and whose shape and
/// parameter range are the same. Moving every defining point by at most
/// the confusion distance moves every point of the curve by no more.
fn translate_of(a: &Curve, b: &Curve, d: Vector, tol: Tolerances) -> Option<f64> {
    let close = |x: f64, y: f64| (x - y).abs() <= tol.confusion();
    let same_range = {
        let ((a0, a1), (b0, b1)) = (a.domain(), b.domain());
        (a0 - b0).abs() <= tol.parametric() && (a1 - b1).abs() <= tol.parametric()
    };
    if !same_range {
        return None;
    }
    match (a, b) {
        (Curve::Line(la), Curve::Line(lb)) => {
            let (aa, ab) = (la.axis(), lb.axis());
            let (a0, a1) = a.domain();
            let reach = a0.abs().max(a1.abs());
            let turned = (aa.direction.vector() - ab.direction.vector()).magnitude() * reach;
            if turned > tol.confusion() {
                return None;
            }
            along_axis(ab.location - aa.location, d, tol)
        }
        (Curve::Circle(ca), Curve::Circle(cb)) => {
            let (ra, rb) = (ca.circle(), cb.circle());
            if ca.is_reversed() != cb.is_reversed() || !close(ra.radius(), rb.radius()) {
                return None;
            }
            let (fa, fb) = (ra.frame(), rb.frame());
            let r = ra.radius();
            if (fa.x().vector() - fb.x().vector()).magnitude() * r > tol.confusion()
                || (fa.y().vector() - fb.y().vector()).magnitude() * r > tol.confusion()
            {
                return None;
            }
            along_axis(fb.origin() - fa.origin(), d, tol)
        }
        (Curve::BSpline(sa), Curve::BSpline(sb)) => {
            let (ka, kb) = (sa.knots(), sb.knots());
            let (pa, pb) = (sa.control_points(), sb.control_points());
            if ka.degree() != kb.degree()
                || ka.knots().len() != kb.knots().len()
                || ka
                    .knots()
                    .iter()
                    .zip(kb.knots())
                    .any(|(x, y)| (x - y).abs() > tol.parametric())
                || a.is_periodic() != b.is_periodic()
                || pa.len() != pb.len()
                || pa.is_empty()
            {
                return None;
            }
            let mut lift: Option<f64> = None;
            for (p, q) in pa.iter().zip(pb) {
                if (p.weight - q.weight).abs() > 1e-12 * p.weight.abs() || p.weight <= 0.0 {
                    return None;
                }
                let moved = along_axis(q.point() - p.point(), d, tol)?;
                match lift {
                    None => lift = Some(moved),
                    Some(first) if close(first, moved) => {}
                    Some(_) => return None,
                }
            }
            lift
        }
        _ => None,
    }
}

/// The plain pcurve `edge` carries on `surface`, where it carries one.
fn pcurve_on(model: &Model, edge: &Shape, surface: SurfaceId) -> Option<(PlanarCurve, (f64, f64))> {
    let data = model.node(edge)?.data().as_edge()?;
    let EdgeRepr::PCurve { curve, range, .. } = data.pcurve_for(surface, edge.location())? else {
        return None;
    };
    Some((model.geometry().pcurve(*curve)?.clone(), *range))
}

/// Whether `edge` carries any parameterization on `surface` already.
fn has_slot(model: &Model, edge: &Shape, surface: SurfaceId) -> bool {
    model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .is_some_and(|d| d.pcurve_for(surface, edge.location()).is_some())
}

/// `pcurve` read along `range` at evenly spaced parameters, backwards where
/// `backwards`.
fn samples(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    backwards: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    let mut out = Vec::with_capacity(SAMPLES as usize);
    for k in 0..SAMPLES {
        let f = f64::from(k) / f64::from(SAMPLES - 1);
        let f = if backwards { 1.0 - f } else { f };
        out.push(pcurve.point_at((range.1 - range.0).mul_add(f, range.0), tol)?);
    }
    Ok(out)
}

/// The chart's period in each direction, where it repeats or closes.
fn periods(surface: &SurfaceGeometry, tol: Tolerances) -> (Option<f64>, Option<f64>) {
    let ((ua, ub), (va, vb)) = surface.domain();
    (
        (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(ub - ua),
        (surface.is_periodic_v() || surface.is_closed_v(tol)).then_some(vb - va),
    )
}

/// Whether two chart points that name one point of the carrier are one
/// chart point too, and not a period apart.
fn same_sheet(a: Point2, b: Point2, periods: (Option<f64>, Option<f64>)) -> bool {
    periods.0.is_none_or(|p| (a.x - b.x).abs() < p / 4.0)
        && periods.1.is_none_or(|p| (a.y - b.y).abs() < p / 4.0)
}

/// Twice the area a chain of chart points sweeps about `anchor`: summed
/// over closed rings, what they enclose, positive counter-clockwise.
fn swept(points: &[Point2], anchor: Point2) -> f64 {
    points
        .windows(2)
        .map(|w| (w[0] - anchor).cross(w[1] - anchor))
        .sum()
}

/// The cluster's members, as indices into the faces, each with its
/// carrier; the first is the keeper, whose chart the merged face takes.
pub(super) struct Members<'a> {
    pub(super) members: &'a [usize],
    pub(super) carriers: &'a [&'a Curved],
}

/// Carry a cluster of faces on one curved carrier into its first face's
/// chart, and chain the merged boundary into rings, the outer one first,
/// each a list of indices into `walks`. `walks` holds each member's edges
/// as it walks them, seen from the side they all face, tagged with the
/// member's index into the faces; `kept` marks those no two members share.
///
/// Each member's chart is mapped onto the keeper's, the period it lands on
/// settled by an edge it shares with a member already placed; the members'
/// pcurves are carried across by that map, which is exact, and measured
/// against the members' own lifts. An edge two members share dissolves
/// where its two sides land on one sheet of the chart, and stays as a seam
/// where they land a period apart (two halves of a bore make a whole one).
/// The rings are chained through the chart, so a seam's two sides are
/// told apart. The cluster is left as it is (`None`) where a member walks
/// an edge twice, a shared edge's sides land neither together nor a period
/// apart, a ring does not close in the chart (a band round the axis with
/// no seam), the rings do not wind as one outer ring round its holes, or
/// the keeper's window cannot be grown to hold the others'. Nothing is
/// changed in the model until every check has passed.
pub(super) fn merge(
    model: &mut Model,
    cluster: &Members<'_>,
    walks: &[(usize, Shape)],
    kept: &[bool],
    reversed: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Vec<usize>>>> {
    let Members { members, carriers } = *cluster;
    let keeper = carriers[0];
    let periods = periods(&keeper.surface, tol);
    let slot_of: FastMap<usize, usize> = members.iter().enumerate().map(|(k, &i)| (i, k)).collect();

    // Each walk's member, as a slot in `members`, and its pcurve there.
    let mut by_node: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    let mut own: Vec<(usize, PlanarCurve, (f64, f64))> = Vec::with_capacity(walks.len());
    for (w, (i, edge)) in walks.iter().enumerate() {
        let slot = slot_of[i];
        let Some((pcurve, range)) = pcurve_on(model, edge, carriers[slot].id) else {
            return Ok(None);
        };
        let uses = by_node.entry(edge.node()).or_default();
        if uses.iter().any(|&o| walks[o].0 == *i) {
            return Ok(None);
        }
        uses.push(w);
        own.push((slot, pcurve, range));
    }
    let shared: Vec<(TShapeId, usize, usize)> = {
        let mut shared: Vec<(TShapeId, usize, usize)> = by_node
            .iter()
            .filter_map(|(node, uses)| match uses.as_slice() {
                &[a, b] => Some((*node, a, b)),
                _ => None,
            })
            .collect();
        shared.sort_by_key(|&(_, a, _)| a);
        shared
    };

    // Each member's chart placed on the keeper's, reached across the edges
    // the members share.
    let mut maps: Vec<Option<ChartMap>> = vec![None; members.len()];
    maps[0] = Some(ChartMap::IDENTITY);
    let mut placed = 1;
    while placed < members.len() {
        let before = placed;
        for &(_, a, b) in &shared {
            let (sa, sb) = (own[a].0, own[b].0);
            let (from, onto, unplaced) = match (maps[sa], maps[sb]) {
                (Some(m), None) => (b, (a, m), sb),
                (None, Some(m)) => (a, (b, m), sa),
                _ => continue,
            };
            let Some(base) = chart_map(carriers[unplaced], keeper, tol) else {
                return Ok(None);
            };
            // Both sides read at the middle of the edge.
            let mid = |w: usize| -> OgeomResult<Point2> {
                let (_, pcurve, range) = &own[w];
                pcurve.point_at(0.5 * (range.0 + range.1), tol)
            };
            let there = onto.1.apply(mid(onto.0)?);
            let here = base.apply(mid(from)?);
            let shift = whole_periods(there - here, periods);
            maps[unplaced] = Some(ChartMap {
                flip: base.flip,
                shift: base.shift + shift,
            });
            placed += 1;
        }
        if placed == before {
            return Ok(None);
        }
    }
    let maps: Vec<ChartMap> = maps.into_iter().flatten().collect();

    // Every walk read in the keeper's chart, in the direction it is
    // walked, each point measured against the member's own lift.
    let mut lines: Vec<Vec<Point2>> = Vec::with_capacity(walks.len());
    for (w, (_, edge)) in walks.iter().enumerate() {
        let (slot, pcurve, range) = &own[w];
        let backwards = edge.orientation() == Orientation::Reversed;
        let mut line = Vec::new();
        for p in samples(pcurve, *range, backwards, tol)? {
            let q = maps[*slot].apply(p);
            if *slot != 0 {
                let (Some(theirs), Some(ours)) = (
                    world_point(carriers[*slot], p, tol),
                    world_point(keeper, q, tol),
                ) else {
                    return Ok(None);
                };
                if theirs.distance(ours) > tol.confusion() {
                    return Ok(None);
                }
            }
            line.push(q);
        }
        lines.push(line);
    }

    // Each shared edge dissolves, or is a seam a whole period across.
    let mut kept = kept.to_vec();
    let mut seams: Vec<(TShapeId, usize, usize)> = Vec::new();
    for &(node, a, b) in &shared {
        let pairs = || lines[a].iter().zip(lines[b].iter().rev());
        if pairs().all(|(p, q)| same_sheet(*p, *q, periods)) {
            continue;
        }
        let (Some(&p), Some(&q)) = (lines[a].first(), lines[b].last()) else {
            return Ok(None);
        };
        let across = whole_periods(p - q, periods);
        if across == Vector2::ZERO || !pairs().all(|(p, q)| same_sheet(*p, *q + across, periods)) {
            return Ok(None);
        }
        kept[a] = true;
        kept[b] = true;
        seams.push((node, a, b));
    }

    // Chain the kept walks into rings, head to tail through the chart: of
    // the walks leaving a vertex, the one starting where the last ended.
    let boundary: Vec<usize> = (0..walks.len()).filter(|&w| kept[w]).collect();
    if boundary.is_empty() {
        return Ok(None);
    }
    let mut ends: Vec<(TShapeId, TShapeId)> = Vec::with_capacity(boundary.len());
    let mut leaving: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for (k, &w) in boundary.iter().enumerate() {
        let Some((start, end)) = ogeom_algo::edge_vertices(model, &walks[w].1)? else {
            return Ok(None);
        };
        leaving.entry(start.node()).or_default().push(k);
        ends.push((start.node(), end.node()));
    }
    let first_of = |k: usize| lines[boundary[k]].first().copied();
    let last_of = |k: usize| lines[boundary[k]].last().copied();
    let mut chained = vec![false; boundary.len()];
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for first in 0..boundary.len() {
        if chained[first] {
            continue;
        }
        let mut ring = Vec::new();
        let mut k = first;
        loop {
            chained[k] = true;
            ring.push(boundary[k]);
            let Some(end) = last_of(k) else {
                return Ok(None);
            };
            let mut next = leaving
                .get(&ends[k].1)
                .map_or(&[][..], Vec::as_slice)
                .iter()
                .copied()
                .filter(|&n| first_of(n).is_some_and(|start| same_sheet(start, end, periods)));
            let (Some(next), None) = (next.next(), next.next()) else {
                return Ok(None);
            };
            if next == first {
                break;
            }
            if chained[next] {
                return Ok(None);
            }
            k = next;
        }
        rings.push(ring);
    }

    // The side the members' boundaries wind about, all alike.
    let anchor = lines[0][0];
    let mut winding = 0.0f64;
    for slot in 0..members.len() {
        let area: f64 = lines
            .iter()
            .zip(&own)
            .filter(|(_, o)| o.0 == slot)
            .map(|(line, _)| swept(line, anchor))
            .sum();
        if area == 0.0 || (winding != 0.0 && area.signum() != winding) {
            return Ok(None);
        }
        winding = area.signum();
    }

    // One merged ring winds as the members do, the rest the other way.
    let mut outer = None;
    for (r, ring) in rings.iter().enumerate() {
        let mut points: Vec<Point2> = Vec::new();
        for &w in ring {
            points.extend_from_slice(&lines[w]);
        }
        points.push(points[0]);
        let area = swept(&points, anchor);
        if area == 0.0 {
            return Ok(None);
        }
        if area.signum() == winding && outer.replace(r).is_some() {
            return Ok(None);
        }
    }
    let Some(outer) = outer else {
        return Ok(None);
    };

    // The keeper's window grown to hold every member's, carried across.
    let Some(grown) = grown_window(keeper, carriers, &maps, tol) else {
        return Ok(None);
    };

    // The boundary's pcurves on the keeper's surface: its own where it has
    // them, the members' carried across where the edge has none there yet,
    // and where it has one already, that one only if it lies where the
    // carried one does.
    let carried = |w: usize| -> OgeomResult<(PlanarCurve, (f64, f64))> {
        let (slot, pcurve, range) = &own[w];
        Ok((pcurve.transformed(&maps[*slot].transform(), tol)?, *range))
    };
    let seam_nodes: Vec<TShapeId> = seams.iter().map(|s| s.0).collect();
    let mut attach: Vec<(Shape, PlanarCurve, (f64, f64))> = Vec::new();
    for (w, (_, edge)) in walks.iter().enumerate() {
        if !kept[w] || own[w].0 == 0 || seam_nodes.contains(&edge.node()) {
            continue;
        }
        let (pcurve, range) = carried(w)?;
        if has_slot(model, edge, keeper.id) {
            let Some((there, there_range)) = pcurve_on(model, edge, keeper.id) else {
                return Ok(None);
            };
            let theirs = samples(&there, there_range, false, tol)?;
            let ours = samples(&pcurve, range, false, tol)?;
            if !theirs
                .iter()
                .zip(&ours)
                .all(|(p, q)| same_sheet(*p, *q, periods))
            {
                return Ok(None);
            }
            continue;
        }
        attach.push((edge.clone(), pcurve, range));
    }
    // A seam's Forward use, as the merged face's own node walks it, reads
    // the side that use walks.
    let mut seamed: Vec<(Shape, PlanarCurve, PlanarCurve, (f64, f64))> = Vec::new();
    for &(_, a, b) in &seams {
        let ((pa, ra), (pb, rb)) = (carried(a)?, carried(b)?);
        let Some(pb) = over_range(pb, rb, ra, tol) else {
            return Ok(None);
        };
        let forward_is_a = (walks[a].1.orientation() == Orientation::Forward) != reversed;
        let (forward, backward) = if forward_is_a { (pa, pb) } else { (pb, pa) };
        seamed.push((walks[a].1.clone(), forward, backward, ra));
    }

    if let Some(grown) = grown
        && let Some(held) = model.geometry_mut().surface_mut(keeper.id)
    {
        *held = grown;
    }
    for (edge, pcurve, range) in attach {
        if !has_slot(model, &edge, keeper.id) {
            ogeom_algo::attach_pcurve(
                model,
                &edge,
                pcurve,
                keeper.id,
                edge.location().clone(),
                range,
            )?;
        }
    }
    for (edge, forward, backward, range) in seamed {
        if let Some(node) = model.node_mut(&edge)
            && let NodeData::Edge(data) = node.data_mut()
        {
            data.representations.retain(|r| {
                r.surface() != Some(keeper.id) || r.location() != Some(edge.location())
            });
        }
        ogeom_algo::attach_seam(
            model,
            &edge,
            forward,
            backward,
            keeper.id,
            edge.location().clone(),
            range,
        )?;
    }
    rings.swap(0, outer);
    Ok(Some(rings))
}

/// `pcurve`, read over `from`, read over `to` instead at the same fraction
/// along: itself where the ranges agree, a chart line slid along itself
/// where they span one length, and `None` otherwise.
fn over_range(
    pcurve: PlanarCurve,
    from: (f64, f64),
    to: (f64, f64),
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let eps = tol.parametric();
    if (from.0 - to.0).abs() <= eps && (from.1 - to.1).abs() <= eps {
        return Some(pcurve);
    }
    let PlanarCurve::Line(line) = pcurve else {
        return None;
    };
    if ((from.1 - from.0) - (to.1 - to.0)).abs() > eps {
        return None;
    }
    let axis = line.axis();
    let origin = axis.location + axis.direction.vector() * (from.0 - to.0);
    ogeom_geom::Line2d::over(ogeom_math::Axis2::new(origin, axis.direction), to.0, to.1)
        .ok()
        .map(Into::into)
}

/// The whole periods nearest `gap` in each direction the chart repeats.
fn whole_periods(gap: Vector2, periods: (Option<f64>, Option<f64>)) -> Vector2 {
    let wrap = |g: f64, period: Option<f64>| period.map_or(0.0, |p| (g / p).round() * p);
    Vector2::new(wrap(gap.x, periods.0), wrap(gap.y, periods.1))
}

/// The keeper's surface with its bounded window grown to hold every
/// member's carried into its chart: `Some(None)` where it holds them
/// already, `None` where a member reaches past what the keeper's chart
/// has (an extrusion's basis curve ends there).
fn grown_window(
    keeper: &Curved,
    carriers: &[&Curved],
    maps: &[ChartMap],
    tol: Tolerances,
) -> Option<Option<SurfaceGeometry>> {
    let ((ua, ub), (va, vb)) = keeper.surface.domain();
    let (pu, _) = periods(&keeper.surface, tol);
    let (mut lo, mut hi) = (va, vb);
    for (c, map) in carriers.iter().zip(maps).skip(1) {
        let ((a, b), (c0, c1)) = c.surface.domain();
        let corners = [map.apply(Point2::new(a, c0)), map.apply(Point2::new(b, c1))];
        if pu.is_none()
            && corners
                .iter()
                .any(|p| p.x < ua - tol.parametric() || p.x > ub + tol.parametric())
        {
            return None;
        }
        for p in corners {
            (lo, hi) = (lo.min(p.y), hi.max(p.y));
        }
    }
    if (lo, hi) == (va, vb) {
        return Some(None);
    }
    match &keeper.surface {
        SurfaceGeometry::Cylinder(s) => ogeom_geom::CylinderSurface::new(s.cylinder(), (lo, hi))
            .ok()
            .map(|s| Some(s.into())),
        SurfaceGeometry::Cone(s) => ogeom_geom::ConeSurface::new(s.cone(), (lo, hi))
            .ok()
            .map(|s| Some(s.into())),
        SurfaceGeometry::Extrusion(s) => {
            ogeom_geom::ExtrusionSurface::over(s.curve().clone(), s.direction(), (lo, hi))
                .ok()
                .map(|s| Some(s.into()))
        }
        // A sphere's latitude and a torus's tube angle span their whole
        // range already.
        _ => Some(None),
    }
}
