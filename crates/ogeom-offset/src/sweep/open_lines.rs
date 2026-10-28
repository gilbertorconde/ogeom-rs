//! An open wire swept down an open spine of straight legs, as one sheet.
//!
//! An open wire bounds no face, so there is no solid to fuse round a
//! corner and take the walls of. The sheet is built from its walls. The
//! rotation-minimizing frame carries the section to each leg: along the
//! legs before it, turned at each corner by the least rotation between the
//! headings. Each leg is the prism of its section, trimmed at a corner by
//! the plane halving it:
//!
//! - mitred, the leg runs on straight past the corner and is cut there. The
//!   turn is the reflection in that plane followed by the reflection in the
//!   plane square to the next leg, so a profile point's path down either
//!   leg meets the plane at one point: the two legs' cuts are one curve;
//! - round, the leg ends square at the corner, and the plane cuts only the
//!   part of the section inside the turn. The part outside is turned about
//!   the corner's axis through the corner's angle, onto the next leg's
//!   section there; a profile edge along the axis turns on itself and
//!   sweeps nothing.
//!
//! Each leg lies between its corners' planes and each round between the
//! planes square to its legs, so neighbours meet only on the edges they
//! share. The pieces are sewn and the sheet must come out a disc: one
//! shell facing one way, every edge of one face or two, one boundary loop
//! running only along the two end sections and the paths of the profile's
//! two ends, and no two faces meeting other than along what they share.

use ogeom_algo::{
    Built, History, check, check_self_intersection, edge_vertices, make_half_space,
    make_natural_face, make_prism, make_revolution, make_wire, sew_within, shape_bounds,
};
use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, PlaneSurface, SurfaceGeometry, Transformable as _};
use ogeom_math::{Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_topo::{Filter, Location, Model, Shape, ShapeType, TShapeId, explore, explore_unique};

use super::{PipeCorners, carried_down_lines, realized_rings, spine_curve_of, wire_ends};

/// Sweep an open wire down an open spine of straight legs under the
/// rotation-minimizing frame, each corner mitred or rounded; see the module
/// documentation. The sheet faces the way the prism of the profile along
/// the first leg does.
pub(super) fn open_wire_lines(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    corners: PipeCorners,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let legs = placed_legs(model, spine, tol)?;
    let unit = |v: Vector| v / v.magnitude();
    let sharp = |a: Vector, b: Vector| {
        let (a, b) = (unit(a), unit(b));
        a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0
    };
    if corners == PipeCorners::Extended && legs.windows(2).any(|pair| sharp(pair[0].1, pair[1].1)) {
        ogeom_bail!(
            Construction,
            "an extended corner squares the outside of a solid's section; an open \
             wire sweeps a sheet, whose corners are mitred or round"
        );
    }
    let Some(&(start, first)) = legs.first() else {
        ogeom_bail!(Construction, "the spine has no leg");
    };
    match ogeom_algo::find_plane(model, profile, tol)? {
        Some(plane) => plane_square_at_start(&plane, start, unit(first), tol)?,
        None => lines_square_at_start(model, profile, start, unit(first), tol)?,
    }
    let (carried, turns) = carried_down_lines(&legs, tol)?;
    let hull = hull_points(model, profile, tol)?;
    let width = {
        let bounds = shape_bounds(model, profile, tol)?;
        match (bounds.low(), bounds.high()) {
            (Some(a), Some(b)) => a.distance(b),
            _ => ogeom_bail!(Construction, "the profile has no extent"),
        }
    };
    // Whether the corner after leg `k` is rounded: a corner that turns
    // nothing has no round, and is cut on the plane square to both legs.
    let rounded = |k: usize| corners == PipeCorners::Round && turns[k].is_some();
    let last = legs.len() - 1;
    // A round corner is told from its square end by the arc its section's
    // farthest point outside the turn sweeps, and by the depth the plane
    // halving the corner cuts its farthest point inside the turn short. A
    // section within the confusion distance of the axis on a side has
    // nothing on that side; past it, an arc or a depth within that distance
    // is no curve or cut the kernel holds.
    for k in 0..last {
        let Some((axis, angle)) = turns[k].filter(|_| rounded(k)) else {
            continue;
        };
        let (d, next) = (unit(legs[k].1), unit(legs[k + 1].1));
        let outward = unit(d * next.dot(d) - next);
        let at_corner = Transform::translation(legs[k].1) * carried[k];
        let (mut outside, mut inside) = (0.0_f64, 0.0_f64);
        for p in &hull {
            let alpha = (at_corner.apply(*p) - axis.location).dot(outward);
            outside = outside.max(alpha);
            inside = inside.max(-alpha);
        }
        let held = |reach: f64, run: f64| reach <= tol.confusion() || run > tol.confusion();
        if !held(outside, outside * angle) || !held(inside, inside * (angle / 2.0).tan()) {
            ogeom_bail!(
                Construction,
                "a round corner turns so little that its round or its cut lies \
                 within the confusion distance of its square end; a mitred \
                 corner sweeps it"
            );
        }
    }

    // Each corner's plane halving it, its normal heading on down the spine:
    // where each leg starts and where it ends.
    let cuts: Vec<(Cut, Cut)> = (0..legs.len())
        .map(|i| {
            let (from, v) = legs[i];
            let d = unit(v);
            (
                (i > 0).then(|| (from, unit(unit(legs[i - 1].1) + d))),
                (i < last).then(|| (from + v, unit(d + unit(legs[i + 1].1)))),
            )
        })
        .collect();
    // Where a profile point's path down leg `i`, from the leg's start,
    // leaves the start's cut and meets the end's, both linear in the point.
    // A rounded corner stops the path of a point outside the turn square at
    // the corner.
    let path = |i: usize, p: Point| -> (f64, f64) {
        let (from, v) = legs[i];
        let (d, length) = (unit(v), v.magnitude());
        let r = carried[i].apply(p) - from;
        let (start_cut, end_cut) = cuts[i];
        let enters = start_cut.map_or(0.0, |(_, n)| -r.dot(n) / d.dot(n));
        let leaves = end_cut.map_or(length, |(_, n)| length - r.dot(n) / d.dot(n));
        let enters = if i > 0 && rounded(i - 1) {
            enters.max(0.0)
        } else {
            enters
        };
        let leaves = if i < last && rounded(i) {
            leaves.min(length)
        } else {
            leaves
        };
        (enters, leaves)
    };

    let mut faces: Vec<Shape> = Vec::new();
    for (i, &(_, v)) in legs.iter().enumerate() {
        let (d, length) = (unit(v), v.magnitude());
        let (start_cut, end_cut) = cuts[i];
        // The extremes of the crossings over the profile are taken at the
        // hull points.
        let mut shortest = f64::INFINITY;
        let (mut before, mut after) = (0.0_f64, 0.0_f64);
        for p in &hull {
            let (enters, leaves) = path(i, *p);
            shortest = shortest.min(leaves - enters);
            before = before.max(-enters);
            after = after.max(leaves - length);
        }
        if shortest <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "a leg is shorter than its corners' reach: the cuts at its ends \
                 cross within the section"
            );
        }
        // A mitred end runs on a section's width past its farthest
        // crossing, so the cut falls inside the prism. A rounded end stops
        // square at the corner.
        let before = if start_cut.is_some() && !rounded(i - 1) {
            before + width
        } else {
            0.0
        };
        let after = if end_cut.is_some() && !rounded(i) {
            after + width
        } else {
            0.0
        };
        let motion = Transform::translation(-d * before) * carried[i];
        let Some(ring) = realized_rings(model, profile, &motion, tol)?
            .into_iter()
            .next()
        else {
            ogeom_bail!(Construction, "the profile has no wire");
        };
        let section = make_wire(model, &ring, tol)?.shape;
        let mut piece = make_prism(model, &section, d * (length + before + after), tol)?.shape;
        if let Some((at, n)) = start_cut {
            piece = kept_side(model, &piece, at, n, tol)?;
        }
        if let Some((at, n)) = end_cut {
            piece = kept_side(model, &piece, at, -n, tol)?;
        }
        faces.extend(explore(model, &piece, Filter::OfType(ShapeType::Face))?);
        if i < last
            && rounded(i)
            && let Some((axis, angle)) = turns[i]
        {
            faces.extend(round_corner(model, &piece, axis, angle, d, tol)?);
        }
    }

    // The sheet's boundary: the profile where the spine starts and its image
    // where it ends, and the paths of the profile's two ends.
    let ends = wire_ends(model, profile, tol)?.map(|(at, _)| at);
    let end_along = legs[last].1;
    let mut boundary = Boundary {
        sections: [
            (ends[0], unit(first)),
            (carried[last].apply(ends[0]) + end_along, unit(end_along)),
        ],
        rails: Vec::new(),
    };
    for end in ends {
        for (i, &(_, v)) in legs.iter().enumerate() {
            let (d, length) = (unit(v), v.magnitude());
            let at = carried[i].apply(end);
            let (enters, leaves) = path(i, end);
            boundary
                .rails
                .push(Rail::Line(at + d * enters, at + d * leaves));
            if i < last
                && rounded(i)
                && leaves >= length - tol.confusion()
                && let Some((axis, angle)) = turns[i]
            {
                boundary.rails.push(Rail::Arc(axis, at + v, angle));
            }
        }
    }
    let shape = one_disc(model, &faces, &boundary, tol)?;
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// The legs of a spine of straight edges where they stand, each edge's
/// placement applied, in the order the spine runs them: each as its start
/// and its displacement.
fn placed_legs(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<(Point, Vector)>> {
    let edges = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        _ => model.ordered_children_of(spine)?,
    };
    let mut legs = Vec::with_capacity(edges.len());
    for edge in &edges {
        let Some((a, b)) = placed_line(model, edge, tol)? else {
            ogeom_bail!(Construction, "a leg of the spine is not straight");
        };
        let (a, b) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
            (b, a)
        } else {
            (a, b)
        };
        legs.push((a, b - a));
    }
    Ok(legs)
}

/// Refuse a profile on `plane` that does not stand square to the spine at
/// its start `at`, heading `heading`: square to the heading to an angle's
/// tolerance, and through the start, as [`profile_square_at_start`] asks.
///
/// [`profile_square_at_start`]: super::profile_square_at_start
fn plane_square_at_start(
    plane: &Plane,
    at: Point,
    heading: Vector,
    tol: Tolerances,
) -> OgeomResult<()> {
    if plane.normal().vector().cross(heading).magnitude() > tol.angular() {
        ogeom_bail!(
            Construction,
            "the profile leans along its spine; a pipe shell runs square to \
             the start"
        );
    }
    if plane.distance_to(at) > tol.confusion() * 100.0 {
        ogeom_bail!(
            Construction,
            "the profile does not sit at the spine's start"
        );
    }
    Ok(())
}

/// Refuse a profile of lines spanning no plane that does not lie in the
/// plane square to the spine at its start `at`, heading `heading`: each
/// line square to the heading to an angle's tolerance, and through that
/// plane, as [`plane_square_at_start`] asks of a planar profile.
fn lines_square_at_start(
    model: &Model,
    profile: &Shape,
    at: Point,
    heading: Vector,
    tol: Tolerances,
) -> OgeomResult<()> {
    for edge in model.ordered_children_of(profile)? {
        let (a, b) = match placed_line(model, &edge, tol)? {
            Some(ends) => ends,
            None => ogeom_bail!(Construction, "a pipe shell sweeps a planar profile"),
        };
        let along = b - a;
        if along.dot(heading).abs() > tol.angular() * along.magnitude() {
            ogeom_bail!(
                Construction,
                "the profile leans along its spine; a pipe shell runs square to \
                 the start"
            );
        }
        if (a - at).dot(heading).abs() > tol.confusion() * 100.0 {
            ogeom_bail!(
                Construction,
                "the profile does not sit at the spine's start"
            );
        }
    }
    Ok(())
}

/// A line edge's ends where it stands; `None` for any other curve.
fn placed_line(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point)>> {
    let (curve, range) = spine_curve_of(model, edge)?;
    let basis = match &curve {
        Curve::Trimmed(t) => t.basis().clone(),
        other => other.clone(),
    };
    if !matches!(basis, Curve::Line(_)) {
        return Ok(None);
    }
    let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
    Ok(Some((
        placed.point_at(range.0, tol)?,
        placed.point_at(range.1, tol)?,
    )))
}

/// Points whose convex hull holds the profile: a line's ends, and the
/// corners of any other edge's bound.
fn hull_points(model: &Model, profile: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let mut out = Vec::new();
    for edge in model.ordered_children_of(profile)? {
        match placed_line(model, &edge, tol)? {
            Some((a, b)) => out.extend([a, b]),
            None => out.extend(shape_bounds(model, &edge, tol)?.corners()),
        }
    }
    Ok(out)
}

/// The part of `piece` on the side of the plane through `at` square to
/// `toward` that `toward` points into.
fn kept_side(
    model: &mut Model,
    piece: &Shape,
    at: Point,
    toward: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    // The plane's patch reaches twice as far as anything of the piece, so
    // its rim stays clear of the piece.
    let reach = shape_bounds(model, piece, tol)?
        .corners()
        .into_iter()
        .map(|p| p.distance(at))
        .fold(0.0_f64, f64::max);
    let plane = Plane::through(at, Direction::new(toward, tol)?);
    let face = make_natural_face(
        model,
        SurfaceGeometry::Plane(PlaneSurface::over(
            plane,
            (-reach * 2.0, reach * 2.0),
            (-reach * 2.0, reach * 2.0),
        )?),
    )?
    .shape;
    let side = make_half_space(model, &face, at + toward / toward.magnitude() * reach, tol)?.shape;
    Ok(ogeom_bool::common(model, piece, &side, tol)?.shape)
}

/// The round at a corner: each edge the leg `piece` ends on square to its
/// heading `d` at the corner, turned about `axis` through `angle`. An edge
/// along the axis turns on itself and is left out.
fn round_corner(
    model: &mut Model,
    piece: &Shape,
    axis: Axis,
    angle: f64,
    d: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let c = axis.location;
    let a = axis.direction.vector();
    let off_axis = |p: Point| {
        let r = p - c;
        (r - a * r.dot(a)).magnitude()
    };
    let mut out = Vec::new();
    for edge in explore_unique(model, piece, ShapeType::Edge)? {
        let points = edge_points(model, &edge, tol)?;
        if points
            .iter()
            .any(|p| (*p - c).dot(d).abs() > tol.confusion())
            || points.iter().all(|p| off_axis(*p) <= tol.confusion())
        {
            continue;
        }
        let swept = make_revolution(model, &edge, axis, angle, tol)?.shape;
        out.extend(explore(model, &swept, Filter::OfType(ShapeType::Face))?);
    }
    Ok(out)
}

/// An edge's ends and the middle of its range, where it stands.
fn edge_points(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<[Point; 3]> {
    let (curve, (t0, t1)) = spine_curve_of(model, edge)?;
    let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
    Ok([
        placed.point_at(t0, tol)?,
        placed.point_at(0.5 * (t0 + t1), tol)?,
        placed.point_at(t1, tol)?,
    ])
}

/// A corner's plane halving it where a leg ends: a point on it and its
/// normal, heading on down the spine.
type Cut = Option<(Point, Vector)>;

/// Where a swept sheet may be open: its two end sections, each on the
/// plane through a point square to a heading, and the paths of the
/// profile's two ends.
struct Boundary {
    sections: [(Point, Vector); 2],
    rails: Vec<Rail>,
}

/// A piece of the path of one of the profile's ends: a segment down a leg,
/// or an arc turned about a round corner's axis from a point through an
/// angle.
enum Rail {
    Line(Point, Point),
    Arc(Axis, Point, f64),
}

impl Rail {
    /// Whether `x` lies on the rail, to the confusion distance.
    fn holds(&self, x: Point, tol: Tolerances) -> bool {
        match *self {
            Rail::Line(a, b) => {
                let ab = b - a;
                let t = ((x - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
                x.distance(a + ab * t) <= tol.confusion()
            }
            Rail::Arc(axis, from, angle) => {
                let (c, a) = (axis.location, axis.direction.vector());
                let (f, r) = (from - c, x - c);
                let (pf, pr) = (f - a * f.dot(a), r - a * r.dot(a));
                let radius = pf.magnitude();
                if (f.dot(a) - r.dot(a)).abs() > tol.confusion()
                    || (radius - pr.magnitude()).abs() > tol.confusion()
                {
                    return false;
                }
                if radius <= tol.confusion() {
                    return true;
                }
                let turned = pf.cross(pr).dot(a).atan2(pf.dot(pr));
                let slack = tol.confusion() / radius;
                turned >= -slack && turned <= angle + slack
            }
        }
    }
}

/// The faces sewn into one sheet that is a disc: one shell facing one way,
/// every edge of one face or two, the edges of one face forming one loop
/// that runs only along `boundary`, Euler characteristic one, no finding
/// but its being open, and no two faces meeting other than along what
/// they share.
fn one_disc(
    model: &mut Model,
    faces: &[Shape],
    boundary: &Boundary,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    // Sewn within the confusion distance, as any sewing joins, and with an
    // edge split where a neighbour's vertices meet it: the trim may split
    // one leg's cut where the next leg's runs whole.
    let sewn = sew_within(model, faces, tol.confusion(), tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_bail!(
            Construction,
            "the swept sheet's pieces do not sew into one sheet"
        );
    };
    if !sewn.edges_walked_one_way(model)?.is_empty() {
        ogeom_bail!(Construction, "the swept sheet's pieces do not face one way");
    }
    type Key = (TShapeId, Location);
    let key = |s: &Shape| (s.node(), s.location().clone());
    let mut uses: FastMap<Key, (Shape, usize)> = FastMap::default();
    let mut face_count = 0_i64;
    for face in explore(model, shell, Filter::OfType(ShapeType::Face))? {
        face_count += 1;
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|data| data.degenerate)
            {
                continue;
            }
            uses.entry(key(&edge))
                .or_insert_with(|| (edge.clone(), 0))
                .1 += 1;
        }
    }
    if uses.values().any(|(_, n)| *n > 2) {
        ogeom_bail!(
            Construction,
            "the swept sheet has an edge of more than two faces"
        );
    }
    // The boundary: each of its vertices on two of its edges, and all of
    // them one loop.
    let mut around: FastMap<Key, Vec<Key>> = FastMap::default();
    for (edge, _) in uses.values().filter(|(_, n)| *n == 1) {
        let points = edge_points(model, edge, tol)?;
        let on_section = boundary.sections.iter().any(|&(at, n)| {
            points
                .iter()
                .all(|p| (*p - at).dot(n).abs() <= tol.confusion())
        });
        let on_rails = points
            .iter()
            .all(|p| boundary.rails.iter().any(|rail| rail.holds(*p, tol)));
        if !on_section && !on_rails {
            ogeom_bail!(
                Construction,
                "the swept sheet is open where its pieces should join"
            );
        }
        let Some((a, b)) = edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "a swept sheet's edge has no vertices");
        };
        around.entry(key(&a)).or_default().push(key(&b));
        around.entry(key(&b)).or_default().push(key(&a));
    }
    let one_loop = !around.is_empty() && around.values().all(|next| next.len() == 2) && {
        let mut seen = ogeom_core::FastSet::default();
        let mut pending: Vec<&Key> = around.keys().take(1).collect();
        while let Some(at) = pending.pop() {
            if seen.insert(at.clone()) {
                pending.extend(around[at].iter());
            }
        }
        seen.len() == around.len()
    };
    let vertex_count = explore_unique(model, shell, ShapeType::Vertex)?.len();
    #[allow(clippy::cast_possible_wrap, reason = "counts of a model's shapes")]
    let euler = vertex_count as i64 - uses.len() as i64 + face_count;
    if !one_loop || euler != 1 {
        ogeom_bail!(
            Construction,
            "the swept sheet's pieces do not join into one sheet with one boundary"
        );
    }
    if check(model, shell, tol)?
        .problems
        .iter()
        .any(|problem| problem.kind != ShapeType::Shell)
    {
        ogeom_bail!(Construction, "the swept sheet does not check clean");
    }
    if !check_self_intersection(model, shell, tol)?.is_empty() {
        ogeom_bail!(Construction, "the swept sheet runs into itself");
    }
    Ok(shell.clone())
}
