//! Mixed asymmetric requests built from original-seat tools and setback openings.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::fallback::{Sense, Tool};
use crate::asym::Asymmetric;
use crate::round::section::RoundSpec;

pub(crate) fn blend(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    bevels: &[(Shape, crate::Chamfer)],
    asymmetric: &[Asymmetric],
    corner: crate::ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for a in asymmetric {
        for d in [a.on_face, a.on_other] {
            if !d.is_finite() || d <= tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "a round set back {d} from its edge rounds nothing"
                );
            }
        }
    }
    let mut scratch = model.clone();
    let corners = read_corners(&scratch, solid, rounds, bevels, asymmetric, tol)?;
    let at_corner = |e: &Shape| corners.edges.iter().any(|c| c.is_same(e));
    let free_rounds: Vec<_> = rounds
        .iter()
        .filter(|(e, _)| !at_corner(e))
        .cloned()
        .collect();
    let free_bevels: Vec<_> = bevels
        .iter()
        .filter(|(e, _)| !at_corner(e))
        .cloned()
        .collect();
    let mut built = Built::from_nothing(solid.clone());
    if !free_bevels.is_empty() {
        let step = match corner {
            crate::ChamferCorner::Point => {
                crate::chamfer::engine_chamfer_mitred(&mut scratch, solid, &free_bevels, tol)?
            }
            crate::ChamferCorner::Facet => {
                crate::chamfer::engine_chamfer(&mut scratch, solid, &free_bevels, tol)?
            }
        }
        .ok_or_else(|| {
            ogeom_core::ogeom_err!(
                Construction,
                "the blend engine could not close the mixed request's bevels"
            )
        })?;
        built = step;
    }
    if !free_rounds.is_empty() {
        let traced = free_rounds
            .iter()
            .map(|(e, r)| Ok((one_edge(&built, e)?, *r)))
            .collect::<OgeomResult<Vec<_>>>()?;
        let step = crate::fillet::rounds_by_engine(&mut scratch, &built.shape, &traced, tol)?
            .ok_or_else(|| {
                ogeom_core::ogeom_err!(
                    Construction,
                    "the blend engine could not close the mixed request's rounds"
                )
            })?;
        built = Built::new(step.shape, built.history.then(&step.history));
    }
    let mut free_asym = Vec::new();
    for a in asymmetric.iter().filter(|a| !at_corner(&a.edge)) {
        let (edge, face) = trace_request(&scratch, &built, a, tol)?;
        free_asym.push(Asymmetric {
            edge,
            face,
            ..a.clone()
        });
    }
    if !free_asym.is_empty() {
        let step = super::conic_tools::asymmetric(&mut scratch, &built.shape, &free_asym, tol)?
            .ok_or_else(|| {
                ogeom_core::ogeom_err!(
                    Construction,
                    "the blend engine could not close the mixed request's asymmetric rounds"
                )
            })?;
        built = Built::new(step.shape, built.history.then(&step.history));
    }
    // An edge with an end at a patched corner stops at its section there;
    // one only at vertices left sharp runs to its end faces.
    let patched = |model: &Model, e: &Shape| -> OgeomResult<bool> {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, e)? else {
            return Ok(false);
        };
        Ok(corners
            .corners
            .iter()
            .any(|(v, ..)| v.is_same(&a) || v.is_same(&b)))
    };
    let mut tools = Vec::new();
    for (e, radius) in rounds.iter().filter(|(e, _)| at_corner(e)) {
        if crate::support::planar_seat(&scratch, solid, e, tol).is_err()
            && let Some(a) = round_on_rim(&scratch, solid, e, *radius, tol)?
        {
            let grouped = patched(&scratch, e)?;
            tools.push((
                e.clone(),
                super::conic_tools::tool(&mut scratch, solid, &a, grouped, tol)?,
            ));
            continue;
        }
        crate::support::ball_fits(&scratch, solid, e, *radius, false, tol)?;
        let seat = crate::support::planar_seat(&scratch, solid, e, tol)?;
        let extend = !patched(&scratch, e)?;
        tools.push((
            e.clone(),
            super::conic_tools::extended_sweep(
                &mut scratch,
                solid,
                e,
                &seat,
                &RoundSpec::Circular(*radius),
                [extend; 2],
                tol,
            )?,
        ));
    }
    for (e, spec) in bevels.iter().filter(|(e, _)| at_corner(e)) {
        let tool = bevel_tool(&mut scratch, solid, e, spec, tol)?;
        tools.push((e.clone(), tool));
    }
    for a in asymmetric.iter().filter(|a| at_corner(&a.edge)) {
        let grouped = patched(&scratch, &a.edge)?;
        tools.push((
            a.edge.clone(),
            super::conic_tools::tool(&mut scratch, solid, a, grouped, tol)?,
        ));
    }
    tools.sort_by_key(|(e, _)| e.node());
    if !tools.is_empty() {
        let tools: Vec<_> = tools.into_iter().map(|(_, t)| t).collect();
        let step =
            super::fallback::applied(&mut scratch, &built.shape, &corners.edges, &tools, tol)?
                .ok_or_else(|| {
                    ogeom_core::ogeom_err!(
                        Construction,
                        "the mixed request's corner tools did not close"
                    )
                })?;
        // The corner tools' faces name the original edges, which the earlier
        // steps may have rebuilt: their generation is kept on the originals.
        let mut history = built.history.then(&step.history);
        for e in &corners.edges {
            for made in step.history.generated(e) {
                if !history.generated(e).iter().any(|f| f.is_same(made)) {
                    history.generate(e, made.clone());
                }
            }
            history.delete(e);
        }
        built = Built::new(step.shape, history);
        for (vertex, frame, setback, kinds) in &corners.corners {
            built = super::mixed_corner::close(
                &mut scratch,
                solid,
                built,
                vertex,
                frame,
                setback,
                kinds,
                ogeom_geom::Continuity::G1,
                tol,
            )?;
        }
    }
    super::guards::solid_fits(&scratch, &built.shape, tol)?;
    *model = scratch;
    Ok(built)
}

fn one_edge(built: &Built, edge: &Shape) -> OgeomResult<Shape> {
    let [one] = built.history.trace(edge) else {
        ogeom_bail!(
            Construction,
            "a round's edge does not come through the bevels as one edge"
        );
    };
    Ok(one.clone())
}

fn bevel_tool(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &crate::Chamfer,
    tol: Tolerances,
) -> OgeomResult<Tool> {
    let crate::chamfer::Legs::Straight(seat, d) =
        crate::chamfer::legs(model, solid, edge, spec, tol)?
    else {
        ogeom_bail!(
            Construction,
            "a corner with an asymmetric round takes bevels on straight edges between planes"
        );
    };
    // `leg` runs off a convex crease; a concave one's faces run on the
    // other way from it.
    let sense = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(0, tol)? * (d[0] * sense);
    let b = seat.leg(1, tol)? * (d[1] * sense);
    let points = [
        seat.start,
        seat.end,
        seat.start + a,
        seat.end + a,
        seat.start + b,
        seat.end + b,
    ];
    let mut rings = vec![
        vec![0, 1, 3, 2],
        vec![0, 4, 5, 1],
        vec![0, 2, 4],
        vec![1, 5, 3],
        vec![2, 3, 5, 4],
    ];
    let inside = seat.start + (seat.end - seat.start) * 0.5 + (a + b) / 3.0;
    for ring in &mut rings {
        let p = points[ring[0]];
        if (points[ring[1]] - p)
            .cross(points[ring[2]] - p)
            .dot(inside - p)
            > 0.0
        {
            ring.reverse();
        }
    }
    let (made, _) = ogeom_algo::make_faceted_solid(model, &points, &rings, &[None; 5], tol)?;
    let faces = explore_unique(model, &made.shape, ShapeType::Face)?
        .into_iter()
        .map(|f| (f, edge.clone()))
        .collect();
    Ok(Tool {
        solid: made.shape,
        faces,
        sense: if seat.convex { Sense::Cut } else { Sense::Fuse },
    })
}

/// An asymmetric round's edge and named face as they stand on what the
/// request's rounds and bevels of one radius built.
fn trace_request(
    model: &Model,
    built: &Built,
    a: &Asymmetric,
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let images = built.history.trace(&a.edge);
    // A rim's band takes its whole turn from any arc of it, so a rim split
    // into arcs by the others stands for itself.
    let one_rim = !images.is_empty() && {
        let circle = |e: &Shape| match crate::support::edge_curve(model, e, tol) {
            Ok((Curve::Circle(c), _)) => Some(c.circle()),
            _ => None,
        };
        match circle(&a.edge) {
            Some(rim) => images.iter().all(|e| {
                circle(e).is_some_and(|c| {
                    c.centre().distance(rim.centre()) <= tol.confusion()
                        && (c.radius() - rim.radius()).abs() <= tol.confusion()
                })
            }),
            None => false,
        }
    };
    let edge = match images {
        [edge] => edge.clone(),
        [edge, ..] if one_rim => edge.clone(),
        _ => ogeom_bail!(
            Construction,
            "an asymmetric round's edge does not come through the request's other blends \
             as one edge"
        ),
    };
    // The named face's image beside the edge, on what the request's other
    // blends left.
    let beside = faces_beside(model, &built.shape, &edge, tol)?;
    let images = built.history.trace(&a.face);
    let Some(face) = beside
        .iter()
        .find(|f| images.iter().any(|g| g.is_same(f)))
        .cloned()
    else {
        ogeom_bail!(
            Construction,
            "the face an asymmetric round's first distance runs along is not one of the \
             edge's two faces"
        );
    };
    Ok((edge, face))
}

/// The two faces beside an edge an asymmetric round takes: a rim's cap
/// and wall, or a straight edge's two planes.
fn faces_beside(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<[Shape; 2]> {
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    if let Curve::Circle(circle) = curve {
        let seat = crate::support::revolved_seat(model, solid, edge, &circle, tol)?;
        return Ok([seat.cap_face, seat.wall_face]);
    }
    if matches!(curve, Curve::Line(_)) && crate::support::hosts_planar(model, solid, edge, tol)? {
        return Ok(crate::support::planar_seat(model, solid, edge, tol)?.faces);
    }
    let mut found: Vec<Shape> = Vec::new();
    for f in ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Face)? {
        let holds = ogeom_topo::explore(
            model,
            &f,
            ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Edge),
        )?
        .iter()
        .any(|e| crate::support::same_occurrence(model, e, edge, tol));
        if holds {
            found.push(f);
        }
    }
    let [a, b] = <[Shape; 2]>::try_from(found).map_err(|_| {
        ogeom_core::ogeom_err!(Construction, "an asymmetric round's edge has not two faces")
    })?;
    Ok([a, b])
}

/// The corners where an asymmetric round meets the request's other edges:
/// each vertex, its frame, its bands' setbacks and how each band ends, and
/// every edge of the request at any of them.
struct Corners {
    /// Each corner: its vertex, frame, setbacks and band ends.
    corners: Vec<(
        Shape,
        crate::engine::vertex_read::CornerFrame,
        Vec<f64>,
        Vec<crate::setback::CornerEnd>,
    )>,
    /// The request's edges at the corners, and at vertices where two of
    /// three edges are blended and their bands meet by intersection.
    edges: Vec<Shape>,
}

/// The corners of `asymmetric`: each vertex where an asymmetric round meets
/// another edge of the request other than as one tangent chain. Refused by
/// name unless it is a convex corner of three planes whose three edges are
/// all in the request, each straight, or a vertex of three edges two of
/// which are blended and of one convexity, each a round or bevel on a
/// straight edge between planes or a round on a rim's arc ending there in
/// a plane through the axis, the third left sharp; and no
/// edge at either meets an edge of the request elsewhere other than at
/// another such vertex.
fn read_corners(
    model: &Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    bevels: &[(Shape, crate::Chamfer)],
    asymmetric: &[Asymmetric],
    tol: Tolerances,
) -> OgeomResult<Corners> {
    use crate::setback::CornerEnd;
    let ends = |edge: &Shape| -> OgeomResult<Vec<Shape>> {
        Ok(match ogeom_algo::edge_vertices(model, edge)? {
            Some((a, b)) if a.is_same(&b) => vec![a],
            Some((a, b)) => vec![a, b],
            None => Vec::new(),
        })
    };
    // Every request edge, by kind: 0 a round, 1 a bevel, 2 an asymmetric
    // round, with its index.
    let mut all: Vec<(Shape, usize, usize)> = Vec::new();
    all.extend(
        rounds
            .iter()
            .enumerate()
            .map(|(i, (e, _))| (e.clone(), 0, i)),
    );
    all.extend(
        bevels
            .iter()
            .enumerate()
            .map(|(i, (e, _))| (e.clone(), 1, i)),
    );
    all.extend(
        asymmetric
            .iter()
            .enumerate()
            .map(|(i, a)| (a.edge.clone(), 2, i)),
    );
    let mut at: std::collections::BTreeMap<ogeom_topo::TShapeId, (Shape, Vec<usize>)> =
        std::collections::BTreeMap::new();
    for (k, (edge, _, _)) in all.iter().enumerate() {
        for v in ends(edge)? {
            at.entry(v.node())
                .or_insert_with(|| (v.clone(), Vec::new()))
                .1
                .push(k);
        }
    }
    let solid_edges = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Edge)?;
    let mut corners = Vec::new();
    let mut sharp: Vec<Shape> = Vec::new();
    let mut edges: Vec<Shape> = Vec::new();
    for (vertex, here) in at.into_values() {
        let asym_here: Vec<usize> = here.iter().copied().filter(|&k| all[k].1 == 2).collect();
        if asym_here.is_empty() || here.len() < 2 {
            continue;
        }
        // Two asymmetric rounds continuing each other: a chain's junction.
        if here.len() == 2 && asym_here.len() == 2 {
            let (a, b) = (&asymmetric[all[here[0]].2], &asymmetric[all[here[1]].2]);
            if crate::asym::continues(model, a, b, &vertex, tol)? {
                continue;
            }
        }
        let mut meeting = 0;
        for e in &solid_edges {
            if ends(e)?.iter().any(|v| v.is_same(&vertex)) {
                meeting += 1;
            }
        }
        // Two of three edges blended, the third left sharp: the two bands
        // meet by intersection, each its own tool run to its end faces,
        // where both are rounds or bevels on straight edges between planes
        // and both cut or both fill. A bevel's tool is its prism between
        // its ends, so each end meets only faces square to it or along it.
        if meeting == 3 && here.len() == 2 {
            let mut convex = Vec::with_capacity(2);
            for &k in &here {
                let seat = if all[k].1 == 1 {
                    match crate::chamfer::legs(model, solid, &all[k].0, &bevels[all[k].2].1, tol) {
                        Ok(crate::chamfer::Legs::Straight(seat, _)) => {
                            square_ended(model, solid, &all[k].0, &seat, tol)?;
                            Some(seat)
                        }
                        _ => None,
                    }
                } else {
                    crate::support::planar_seat(model, solid, &all[k].0, tol).ok()
                };
                let sense = match seat {
                    Some(seat) => Some(seat.convex),
                    None if all[k].1 == 2 => {
                        rim_arc_sense(model, solid, &asymmetric[all[k].2], &vertex, tol)?
                    }
                    None if all[k].1 == 0 => {
                        match round_on_rim(model, solid, &all[k].0, rounds[all[k].2].1, tol)? {
                            Some(a) => rim_arc_sense(model, solid, &a, &vertex, tol)?,
                            None => None,
                        }
                    }
                    None => None,
                };
                convex.push(sense);
            }
            // A bevel and an asymmetric round on a straight edge between
            // planes, one cutting and the other filling: bevels first
            // (ADR 0001, section 13, question 3), the round then taken on
            // what the bevel left, its end meeting the bevel's plane.
            let kinds = [all[here[0]].1, all[here[1]].1];
            let straight_asym = here.iter().all(|&k| {
                all[k].1 != 2 || crate::support::planar_seat(model, solid, &all[k].0, tol).is_ok()
            });
            if let [Some(a), Some(b)] = convex[..]
                && a != b
                && kinds.contains(&1)
                && kinds.contains(&2)
                && straight_asym
            {
                continue;
            }
            if let [Some(a), Some(b)] = convex[..]
                && a == b
            {
                for &k in &here {
                    if !edges.iter().any(|e| e.is_same(&all[k].0)) {
                        edges.push(all[k].0.clone());
                    }
                }
                sharp.push(vertex);
                continue;
            }
            ogeom_bail!(
                Construction,
                "an asymmetric round meets another blended edge at a vertex whose third \
                 edge is left sharp, and the two bands meet there only where both cut or \
                 both fill and each is a round or a bevel on a straight edge between \
                 planes, or a round on a rim's arc ending there in a plane through the \
                 rim's axis, or where a convex and a concave round of a \
                 request of rounds alone close on one patch; see docs/CAPABILITIES.md, \
                 fillet.asymmetric-rounds"
            );
        }
        if meeting != 3 || here.len() != 3 {
            ogeom_bail!(
                Construction,
                "an asymmetric round meets another blended edge at a vertex, and such a \
                 corner is closed only where three edges meet and all three are blended; \
                 see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
            );
        }
        let largest = here
            .iter()
            .map(|&k| match all[k].1 {
                0 => rounds[all[k].2].1,
                1 => 0.0,
                _ => asymmetric[all[k].2]
                    .on_face
                    .max(asymmetric[all[k].2].on_other),
            })
            .fold(0.0, f64::max)
            .max(tol.confusion() * 10.0);
        let Ok(frame) =
            crate::engine::vertex_read::corner_frame(model, solid, &vertex, largest, tol)
        else {
            ogeom_bail!(
                Construction,
                "an asymmetric round meets other blended edges at a vertex that is not a \
                 convex corner of planes, and the corner is not closed there; see \
                 docs/CAPABILITIES.md, fillet.asymmetric-rounds"
            );
        };
        let n = frame.d.len();
        let point = |v: &Shape| -> OgeomResult<Point> {
            let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Dangling, "a vertex is not in this model");
            };
            Ok(v.transform(model.datums())?.apply(data.point))
        };
        let mut widths = Vec::with_capacity(n);
        let mut kinds = Vec::with_capacity(n);
        for k in 0..n {
            let (ma, mb) = (frame.m[frame.ring[k]], frame.m[frame.ring[(k + 1) % n]]);
            let dk = frame.d[k].vector();
            let mut found = None;
            for &j in &here {
                let edge = &all[j].0;
                let Some((a, z)) = ogeom_algo::edge_vertices(model, edge)? else {
                    continue;
                };
                let far = if a.is_same(&vertex) { z } else { a };
                let p = point(&far)? - frame.corner;
                let length = p.magnitude();
                if length <= tol.confusion() || p.dot(dk) / length < 1.0 - tol.angular().sqrt() {
                    continue;
                }
                // The width on each of the frame's two faces.
                let on = |normals: [Vector; 2], d: [f64; 2]| -> OgeomResult<(f64, f64)> {
                    let pick = |m: Vector| {
                        (0..2).find(|&i| normals[i].cross(m).magnitude() <= tol.angular() * 10.0)
                    };
                    let (Some(ia), Some(ib)) = (pick(ma), pick(mb)) else {
                        ogeom_bail!(Construction, "a band at a corner does not lie on its faces");
                    };
                    Ok((d[ia], d[ib]))
                };
                found = Some(match all[j].1 {
                    0 => {
                        let cos = -ma.dot(mb);
                        let sin = (1.0 - cos * cos).max(0.0).sqrt();
                        let w = rounds[all[j].2].1 * (1.0 + cos) / sin;
                        ((w, w), CornerEnd::Round)
                    }
                    1 => {
                        let crate::chamfer::Legs::Straight(seat, d) =
                            crate::chamfer::legs(model, solid, edge, &bevels[all[j].2].1, tol)?
                        else {
                            ogeom_bail!(
                                Construction,
                                "a corner with an asymmetric round takes bevels on straight \
                                 edges between planes"
                            );
                        };
                        (on(seat.normals, d)?, CornerEnd::Bevel)
                    }
                    _ => {
                        let a = &asymmetric[all[j].2];
                        let Some(crate::asym::Section { d }) = crate::asym::straight_section(
                            model, solid, edge, &a.face, a.on_face, a.on_other, a.profile, tol,
                        )?
                        else {
                            ogeom_bail!(Construction, "an asymmetric round asked as a round");
                        };
                        let seat = crate::support::planar_seat(model, solid, edge, tol)?;
                        (on(seat.normals, d)?, CornerEnd::Round)
                    }
                });
            }
            let Some((w, kind)) = found else {
                ogeom_bail!(
                    Construction,
                    "an edge at a corner of asymmetric rounds is not one of the request's, or \
                     is curved"
                );
            };
            widths.push(w);
            kinds.push(kind);
        }
        let setback = crate::setback::least_setbacks_widths(&frame, &widths, tol)?;
        for &k in &here {
            if !edges.iter().any(|e| e.is_same(&all[k].0)) {
                edges.push(all[k].0.clone());
            }
        }
        corners.push((vertex, frame, setback, kinds));
    }
    // A corner's edge meets the rest of the request only at corners.
    for e in &edges {
        for v in ends(e)? {
            let corner_here =
                corners.iter().any(|(c, ..)| c.is_same(&v)) || sharp.iter().any(|c| c.is_same(&v));
            let others_here = all.iter().any(|(o, _, _)| {
                !o.is_same(e) && ends(o).is_ok_and(|w| w.iter().any(|x| x.is_same(&v)))
            });
            if others_here && !corner_here {
                ogeom_bail!(
                    Construction,
                    "an edge at a corner of asymmetric rounds meets another blended edge at \
                     its far end, and that meeting is not closed yet; see \
                     docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                );
            }
        }
    }
    Ok(Corners { corners, edges })
}

/// Whether the asymmetric round `a` on an arc of a rim cuts (`Some(true)`)
/// or fills (`Some(false)`) where its band meets another at `vertex` by
/// intersection: its tool is the conic revolved over the arc's own turn,
/// capped in the planes through the axis at the arc's ends, so every face
/// following the rim at `vertex` must be such a plane. `None` for any other
/// edge or end.
fn rim_arc_sense(
    model: &Model,
    solid: &Shape,
    a: &Asymmetric,
    vertex: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    let rim = match crate::asym::rim_conic(
        model, solid, &a.edge, &a.face, a.on_face, a.on_other, a.profile, tol,
    ) {
        Ok(Some(rim)) if !rim.whole => rim,
        Ok(_) | Err(ogeom_core::OgeomError::Construction(_)) => return Ok(None),
        Err(e) => return Err(e),
    };
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Dangling, "a vertex is not in this model");
    };
    let at = vertex.transform(model.datums())?.apply(data.point);
    let hosts = [&rim.seat.cap_face, &rim.seat.wall_face];
    let Some(planes) = crate::support::following_planes(model, solid, &a.edge, hosts, at)? else {
        return Ok(None);
    };
    let radial = planes.iter().all(|(_, origin, normal)| {
        normal.dot(rim.seat.up).abs() <= tol.angular()
            && (rim.seat.centre - *origin).dot(*normal).abs() <= tol.confusion()
    });
    Ok(radial.then_some(!rim.seat.additive()))
}

/// The circular round of `radius` on an open arc of a rim as the
/// asymmetric round of equal distances with the elliptical section, which
/// is that round: each leg the ball's tangent from the rim, `r cot(alpha /
/// 2)` where cap and wall meet at `alpha`. `None` for any other edge.
fn round_on_rim(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Asymmetric>> {
    let Ok((Curve::Circle(circle), _)) = crate::support::edge_curve(model, edge, tol) else {
        return Ok(None);
    };
    let Ok(seat) = crate::support::revolved_seat(model, solid, edge, &circle, tol) else {
        return Ok(None);
    };
    let profile = crate::RoundProfile::Elliptical;
    let face = seat.cap_face.clone();
    let rim = match crate::asym::rim_conic(model, solid, edge, &face, radius, radius, profile, tol)
    {
        Ok(Some(rim)) if !rim.whole => rim,
        Ok(_) | Err(ogeom_core::OgeomError::Construction(_)) => return Ok(None),
        Err(e) => return Err(e),
    };
    let [wall, at, cap] = rim.points;
    let (a, b) = (wall - at, cap - at);
    let cos = (a / a.magnitude()).dot(b / b.magnitude()).clamp(-1.0, 1.0);
    let sin = (1.0 - cos * cos).max(0.0).sqrt();
    if sin <= tol.angular() {
        return Ok(None);
    }
    let leg = radius * (1.0 + cos) / sin;
    Ok(Some(Asymmetric {
        edge: edge.clone(),
        face,
        on_face: leg,
        on_other: leg,
        profile,
    }))
}

/// Refuses by name where an end of the straight bevel `edge` meets a face
/// of `solid`, other than its seat's two, that is not a plane square to the
/// edge or running along it: the bevel's tool stops square at its ends.
fn square_ended(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &crate::support::Seat,
    tol: Tolerances,
) -> OgeomResult<()> {
    for at in [seat.start, seat.end] {
        let hosts = [&seat.faces[0], &seat.faces[1]];
        let planes = crate::support::following_planes(model, solid, edge, hosts, at)?;
        let square = planes.is_some_and(|planes| {
            planes.iter().all(|(_, _, normal)| {
                let rate = normal.dot(seat.along).abs();
                rate <= tol.angular() || rate >= 1.0 - tol.angular()
            })
        });
        if !square {
            ogeom_bail!(
                Construction,
                "a bevel meeting an asymmetric round at a vertex whose third edge is left \
                 sharp ends at a face not square to it, and its prism stops square there; \
                 see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
            );
        }
    }
    Ok(())
}
