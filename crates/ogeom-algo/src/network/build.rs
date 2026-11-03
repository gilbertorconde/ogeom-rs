//! A fitted network as faces, and how well its patches join.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{BSpline2d, BSplineCurve, Curve, Curve3d as _, Line2d, PlanarCurve, Surface as _};
use ogeom_math::{KnotVector, Point2};
use ogeom_topo::{EdgeData, Location, Model, NodeData, Shape};

use super::topology::{Topology, Touch};
use super::{ChainedPatch, IrregularCorner, Round, corner_chart, walk_share};
use crate::build::{
    attach_pcurve, make_edge_between, make_face_on, make_shell, make_vertex, make_wire,
};

/// The faces, their shell and the curves' edges.
pub(super) struct Built {
    pub faces: Vec<Shape>,
    pub shell: Shape,
    pub edges: Vec<Shape>,
}

/// How well the patches join.
pub(super) struct Measured {
    pub worst_tangent: f64,
    pub worst_gap: f64,
    pub corners: Vec<IrregularCorner>,
}

/// Curve `c` of the fitted round, on its strips' knots.
fn curve_of(
    topo: &Topology,
    round: &Round,
    c: usize,
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let n = round.layout.along(c);
    let poles = (0..n)
        .map(|q| round.poles[round.layout.curve_slot(topo, c, q)])
        .collect();
    BSplineCurve::new(round.curve_knots[c].clone(), poles, tol)
}

/// Where the walk along side `k` stands on the chart at share `f` of it.
fn on_side(k: usize, f: f64) -> Point2 {
    let (start, end) = (corner_chart(k), corner_chart(k + 1));
    Point2::new(
        (end.0 - start.0).mul_add(f, start.0),
        (end.1 - start.1).mul_add(f, start.1),
    )
}

/// The share of side `k`'s walk of patch `p` that its `i`-th curve
/// covers in `round`.
fn share(patches: &[ChainedPatch], round: &Round, (p, k): (usize, usize), i: usize) -> (f64, f64) {
    walk_share(
        k,
        i,
        patches[p].sides[k].len(),
        &round.patch_bounds[p][k % 2],
    )
}

/// Where curve `c`'s parameter `t` lies on the chart of the patch whose
/// side `k` it is on.
fn chart_of(
    patches: &[ChainedPatch],
    round: &Round,
    (p, k): (usize, usize),
    c: usize,
    t: f64,
) -> Point2 {
    let side = &patches[p].sides[k];
    let i = side.iter().position(|s| s.curve == c).unwrap_or(0);
    let (f0, f1) = share(patches, round, (p, k), i);
    let (from, to) = (on_side(k, f0), on_side(k, f1));
    let (start, end) = if side[i].reversed {
        (to, from)
    } else {
        (from, to)
    };
    Point2::new(
        (end.x - start.x).mul_add(t, start.x),
        (end.y - start.y).mul_add(t, start.y),
    )
}

/// The faces of the round, sharing an edge per curve, in one shell.
pub(super) fn faces(
    model: &mut Model,
    topo: &Topology,
    patches: &[ChainedPatch],
    round: &Round,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let vertices: Vec<Shape> = (0..topo.vertices)
        .map(|v| make_vertex(model, round.poles[v]).shape)
        .collect();
    let mut edges = Vec::with_capacity(topo.users.len());
    for (c, &(a, b)) in topo.curve_ends.iter().enumerate() {
        let curve = Curve::BSpline(curve_of(topo, round, c, tol)?);
        edges.push(
            make_edge_between(model, curve, (0.0, 1.0), &vertices[a], &vertices[b], tol)?.shape,
        );
    }
    let mut faces = Vec::with_capacity(patches.len());
    for (p, patch) in patches.iter().enumerate() {
        let surface = model
            .geometry_mut()
            .add_surface(ogeom_geom::SurfaceGeometry::BSpline(
                round.surfaces[p].clone(),
            ));
        let mut walk = Vec::with_capacity(4);
        let mut traced: Vec<(Shape, PlanarCurve)> = Vec::with_capacity(4);
        for (k, side) in patch.sides.iter().enumerate() {
            if side.is_empty() {
                let apex = &vertices[topo.corner_vertex[p][k]];
                let mut data = EdgeData::new();
                data.degenerate = true;
                let edge = model.add_edge(data, &[apex.clone(), apex.clone()])?;
                traced.push((
                    edge.clone(),
                    Line2d::segment(on_side(k, 0.0), on_side(k, 1.0), tol)?.into(),
                ));
                walk.push(edge);
                continue;
            }
            for (i, s) in side.iter().enumerate() {
                let (f0, f1) = share(patches, round, (p, k), i);
                let (from, to) = (on_side(k, f0), on_side(k, f1));
                let edge = &edges[s.curve];
                let (start, end) = if s.reversed { (to, from) } else { (from, to) };
                // The pcurve runs over the edge's unit interval: a whole
                // side is a unit segment, by length; a part of one, a
                // segment of degree one over the unit interval.
                let pcurve: PlanarCurve = if side.len() == 1 {
                    Line2d::segment(start, end, tol)?.into()
                } else {
                    BSpline2d::new(KnotVector::clamped_uniform(1, 2)?, vec![start, end], tol)?
                        .into()
                };
                traced.push((edge.clone(), pcurve));
                walk.push(if s.reversed {
                    edge.reversed()
                } else {
                    edge.clone()
                });
            }
        }
        let wire = make_wire(model, &walk, tol)?.shape;
        let face = make_face_on(model, surface, &[wire], tol)?.shape;
        for (edge, pcurve) in traced {
            attach_pcurve(
                model,
                &edge,
                pcurve,
                surface,
                Location::identity(),
                (0.0, 1.0),
            )?;
            // The pcurve is the patch's parameter line along a boundary row
            // of poles the edge's curve is drawn on, over the same knots:
            // the same point at every parameter.
            if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
                data.assert_same_parameter(true);
            }
        }
        faces.push(face);
    }
    let shell = make_shell(model, &faces)?.shape;
    Ok(Built {
        faces,
        shell,
        edges,
    })
}

/// The parameters each curve is read at: the ends and the Gauss nodes of
/// every span.
fn nodes(curve: &BSplineCurve) -> Vec<f64> {
    let mut at = Vec::new();
    let distinct = curve.knots().distinct();
    for w in distinct.windows(2) {
        at.push(w[0].0);
        at.extend(
            ogeom_math::gauss_legendre_rule(w[0].0, w[1].0)
                .iter()
                .map(|n| n.0),
        );
    }
    if let Some(last) = distinct.last() {
        at.push(last.0);
    }
    at
}

/// The angle between two unit normals.
fn angle(a: ogeom_math::Vector, b: ogeom_math::Vector) -> f64 {
    a.cross(b).magnitude().atan2(a.dot(b))
}

/// How well the round's patches join: across every curve, the distance
/// from its edge to each patch's boundary, and across every curve held
/// tangent and at every irregular corner, the angle between the normals.
pub(super) fn measure(
    topo: &Topology,
    patches: &[ChainedPatch],
    round: &Round,
    tol: Tolerances,
) -> OgeomResult<Measured> {
    let mut worst_gap = 0.0_f64;
    let mut worst_tangent = 0.0_f64;
    for c in 0..topo.users.len() {
        let curve = curve_of(topo, round, c, tol)?;
        let (a, b) = topo.curve_ends[c];
        for t in nodes(&curve) {
            let on = curve.point_at(t, tol)?;
            let mut normals = Vec::with_capacity(2);
            for &user in &topo.users[c] {
                let chart = chart_of(patches, round, user, c, t);
                let surface = &round.surfaces[user.0];
                worst_gap = worst_gap.max(surface.point_at(chart.x, chart.y, tol)?.distance(on));
                let at_apex = (t <= 0.0 && topo.apex[a]) || (t >= 1.0 && topo.apex[b]);
                if topo.tangent[c] && !at_apex {
                    normals.push(surface.normal_at(chart.x, chart.y, tol)?.vector());
                }
            }
            if let [n1, n2] = normals.as_slice() {
                worst_tangent = worst_tangent.max(angle(*n1, *n2));
            }
        }
    }
    let mut corners = Vec::new();
    for v in 0..topo.vertices {
        let Some(valence) = topo.interior_valence[v] else {
            continue;
        };
        let all_tangent = topo
            .curve_ends
            .iter()
            .enumerate()
            .filter(|&(_, &(a, b))| a == v || b == v)
            .all(|(c, _)| topo.tangent[c]);
        if valence == 4 || !all_tangent {
            continue;
        }
        let mut normals = Vec::with_capacity(valence);
        for &(p, touch) in &topo.touching[v] {
            let chart = match touch {
                Touch::Corner(k) => {
                    let (u, w) = corner_chart(k);
                    Point2::new(u, w)
                }
                Touch::Through(k, j) => on_side(k, share(patches, round, (p, k), j).1),
            };
            normals.push(round.surfaces[p].normal_at(chart.x, chart.y, tol)?.vector());
        }
        let mut widest = 0.0_f64;
        for (i, a) in normals.iter().enumerate() {
            for b in &normals[i + 1..] {
                widest = widest.max(angle(*a, *b));
            }
        }
        worst_tangent = worst_tangent.max(widest);
        corners.push(IrregularCorner {
            point: round.poles[v],
            valence,
            angle: widest,
        });
    }
    Ok(Measured {
        worst_tangent,
        worst_gap,
        corners,
    })
}
