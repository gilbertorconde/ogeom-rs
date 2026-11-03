//! A fitted network as faces, and how well its patches join.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{BSplineCurve, Curve, Curve3d as _, Line2d, Surface as _};
use ogeom_math::Point2;
use ogeom_topo::{EdgeData, Location, Model, NodeData, Shape};

use super::topology::Topology;
use super::{IrregularCorner, PatchDef, Round, corner_chart};
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

/// Curve `c` of the fitted round, on its strip's knots.
fn curve_of(
    topo: &Topology,
    round: &Round,
    c: usize,
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let n = round.layout.along(topo, c);
    let poles = (0..n)
        .map(|q| round.poles[round.layout.curve_slot(topo, c, q)])
        .collect();
    BSplineCurve::new(round.knots[topo.curve_strip[c]].clone(), poles, tol)
}

/// Where curve `c`'s parameter `t` lies on the chart of the patch whose
/// side `k` it is.
fn chart_of(patches: &[PatchDef], (p, k): (usize, usize), c: usize, t: f64) -> Point2 {
    let reversed = patches[p].sides[k].is_some_and(|s| s.curve == c && s.reversed);
    let (start, end) = if reversed {
        (corner_chart(k + 1), corner_chart(k))
    } else {
        (corner_chart(k), corner_chart(k + 1))
    };
    Point2::new(
        (end.0 - start.0).mul_add(t, start.0),
        (end.1 - start.1).mul_add(t, start.1),
    )
}

/// The faces of the round, sharing an edge per curve, in one shell.
pub(super) fn faces(
    model: &mut Model,
    topo: &Topology,
    patches: &[PatchDef],
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
        let mut traced: Vec<(Shape, ogeom_geom::PlanarCurve)> = Vec::with_capacity(4);
        for (k, side) in patch.sides.iter().enumerate() {
            let (from, to) = (corner_chart(k), corner_chart(k + 1));
            let (from, to) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
            match side {
                Some(s) => {
                    let edge = &edges[s.curve];
                    let (start, end) = if s.reversed { (to, from) } else { (from, to) };
                    traced.push((edge.clone(), Line2d::segment(start, end, tol)?.into()));
                    walk.push(if s.reversed {
                        edge.reversed()
                    } else {
                        edge.clone()
                    });
                }
                None => {
                    let apex = &vertices[topo.corner_vertex[p][k]];
                    let mut data = EdgeData::new();
                    data.degenerate = true;
                    let edge = model.add_edge(data, &[apex.clone(), apex.clone()])?;
                    traced.push((edge.clone(), Line2d::segment(from, to, tol)?.into()));
                    walk.push(edge);
                }
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
    patches: &[PatchDef],
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
                let chart = chart_of(patches, user, c, t);
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
        for (p, corners) in topo.corner_vertex.iter().enumerate() {
            if let Some(k) = corners.iter().position(|&w| w == v) {
                let (u, w) = corner_chart(k);
                normals.push(round.surfaces[p].normal_at(u, w, tol)?.vector());
            }
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
