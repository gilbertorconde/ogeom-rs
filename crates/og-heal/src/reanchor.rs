//! Re-anchoring the rings of a seamless periodic face.
//!
//! The first healing operation an imported file demanded by name. A STEP
//! face on a periodic surface may arrive bounded by two closed rings and no
//! seam; where the rings' vertices share an angle, the reader synthesises
//! the seam itself, but where they sit at *different* angles no seam can
//! join them — the vertices are load-bearing, and one of them is in the
//! wrong place. Moving it is surgery: the ring's closed edge is rebuilt
//! anchored at the aligned angle, and every face that shares the edge — the
//! neighbour on the other side of the ring — is rebuilt to use the new one,
//! or the shell tears exactly where the file was trying to close it.

use og_algo::{
    Built, History, attach_seam, is_shell_closed, make_edge_between, make_face_on, make_shell,
    make_solid, make_vertex, make_wire, surface_iso_u_curve,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve2d as _;
use og_geom::Curve3d as _;
use og_geom::Surface as _;
use og_geom::{Curve, PlanarCurve, SurfaceGeometry};
use og_math::Point;
use og_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};
use std::collections::HashMap;

/// Re-anchor misaligned ring vertices so every periodic face can carry a seam.
///
/// Faces already whole pass through untouched; a shape with nothing to heal
/// comes back as itself with an empty history.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the shape is
/// not a solid, or a broken face's structure resists the repair — in which
/// case nothing is modified.
pub fn reanchor_periodic_rings(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgResult<Built> {
    if model.kind_of(shape)? != ShapeType::Solid {
        og_bail!(Construction, "ring re-anchoring heals solids");
    }

    // --- find the broken faces ---------------------------------------------
    struct Broken {
        face: Shape,
        surface_id: og_topo::SurfaceId,
        surface: SurfaceGeometry,
        /// Each ring: its single closed edge, that edge's vertex, curve and
        /// range.
        rings: Vec<(Shape, Shape, Curve, (f64, f64))>,
    }
    let mut revolution: Vec<Broken> = Vec::new();
    let mut any_broken = false;
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        if !surface.is_periodic_u() {
            continue;
        }
        let surface = surface.clone();
        let surface_id = data.surface;
        let wires = explore(model, &face, Filter::OfType(ShapeType::Wire))?;

        // The face's ring edges: its closed non-seam boundaries, whether it
        // arrived broken (two single-ring wires, no seam) or already seamed
        // (one wire in which the seam edge appears twice).
        let mut ring_edges: Vec<Shape> = Vec::new();
        let mut seen_twice: Vec<TShapeId> = Vec::new();
        for wire in &wires {
            let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
            for edge in &edges {
                if edges.iter().filter(|e| e.node() == edge.node()).count() == 2 {
                    if !seen_twice.contains(&edge.node()) {
                        seen_twice.push(edge.node());
                    }
                    continue;
                }
                let Some((a, b)) = og_algo::edge_vertices(model, edge)? else {
                    continue;
                };
                if a.is_same(&b) && !ring_edges.iter().any(|r| r.node() == edge.node()) {
                    ring_edges.push(edge.clone());
                }
            }
        }
        if ring_edges.len() != 2 {
            continue;
        }
        let broken_here = wires.len() == 2 && seen_twice.is_empty();
        any_broken |= broken_here;
        let mut rings = Vec::new();
        for edge in ring_edges {
            let Some((a, _)) = og_algo::edge_vertices(model, &edge)? else {
                continue;
            };
            let Some(edge_node) = model.node(&edge) else {
                continue;
            };
            let Some(edge_data) = edge_node.data().as_edge() else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                continue;
            };
            rings.push((edge, a, geometry.clone(), *range));
        }
        if rings.len() != 2 {
            continue;
        }
        revolution.push(Broken {
            face,
            surface_id,
            surface,
            rings,
        });
    }
    if !any_broken {
        return Ok(Built::new(shape.clone(), History::new()));
    }
    let broken = revolution;

    // --- one half-plane for the whole chain ----------------------------------
    //
    // Re-anchoring one ring per face is not enough: rings are shared, and a
    // neighbouring face's seam ties to the vertex being moved. For a coaxial
    // chain — which a fillet stack is — a single half-plane through the axis
    // meets every ring exactly once, and anchoring them all there satisfies
    // every face's equal-angle constraint at one stroke.
    let mut history = History::new();
    let (axis_point, axis_dir, radial) = {
        let (_, vertex, curve, _) = &broken[0].rings[0];
        let Curve::Circle(c) = curve else {
            og_bail!(
                Construction,
                "a ring to re-anchor is not a circle; the repair does not \
                 know its parameterization"
            );
        };
        let circle = c.circle();
        let Some(node) = model.node(vertex) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            og_bail!(Construction, "vertex node holds no vertex data");
        };
        let axis_dir = circle.frame().z().vector();
        let r = data.point - circle.centre();
        let radial = r - axis_dir * r.dot(axis_dir);
        (circle.centre(), axis_dir, radial / radial.magnitude())
    };

    let mut substitution: HashMap<TShapeId, Shape> = HashMap::new();
    for b in &broken {
        for (old_edge, old_vertex, curve, _) in &b.rings {
            if substitution.contains_key(&old_edge.node()) {
                continue;
            }
            let Curve::Circle(c) = curve else {
                og_bail!(
                    Construction,
                    "a ring to re-anchor is not a circle; the repair does \
                     not know its parameterization"
                );
            };
            let circle = c.circle();
            // Coaxial or nothing: the half-plane trick needs one axis.
            if circle.frame().z().vector().cross(axis_dir).magnitude() > 1e-6
                || (circle.centre() - axis_point).cross(axis_dir).magnitude()
                    > tol.confusion() * 1e3
            {
                og_bail!(
                    Construction,
                    "the rings are not coaxial; the half-plane repair does \
                     not apply"
                );
            }
            let target = circle.centre() + radial * circle.radius();
            let current = {
                let Some(node) = model.node(old_vertex) else {
                    og_bail!(Dangling, "vertex is not in this model");
                };
                let Some(data) = node.data().as_vertex() else {
                    og_bail!(Construction, "vertex node holds no vertex data");
                };
                data.point
            };
            if current.distance(target) <= tol.confusion() * 1e2 {
                continue;
            }
            let Some(t_star) = circle_parameter(curve, target) else {
                og_bail!(Construction, "an anchor point fell off its circle");
            };
            let period = {
                let (lo, hi) = curve.domain();
                hi - lo
            };
            let vertex = make_vertex(model, target).shape;
            let rebuilt = make_edge_between(
                model,
                curve.clone(),
                (t_star, t_star + period),
                &vertex,
                &vertex,
                tol,
            )?
            .shape;
            history.modify(old_edge, rebuilt.clone());
            history.delete(old_vertex);
            substitution.insert(old_edge.node(), rebuilt);
        }
    }

    // --- rebuild every face that touches a substituted edge ------------------
    let mut face_map: HashMap<TShapeId, Shape> = HashMap::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let uses_any = explore_unique(model, &face, ShapeType::Edge)?
            .iter()
            .any(|e| substitution.contains_key(&e.node()));
        if !uses_any {
            continue;
        }
        let rebuilt = if let Some(b) = broken.iter().find(|b| b.face.is_same(&face)) {
            rebuild_broken_face(
                model,
                b.surface_id,
                &b.surface,
                &b.rings,
                &substitution,
                tol,
            )?
        } else {
            rebuild_plain_face(model, &face, &substitution, tol)?
        };
        let rebuilt = if face.orientation() == og_topo::Orientation::Reversed {
            rebuilt.reversed()
        } else {
            rebuilt
        };
        history.modify(&face, rebuilt.clone());
        face_map.insert(face.node(), rebuilt);
    }

    // --- rebuild shells and the solid ----------------------------------------
    let mut shells = Vec::new();
    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        let faces: Vec<Shape> = explore(model, &shell, Filter::OfType(ShapeType::Face))?
            .into_iter()
            .map(|f| {
                face_map.get(&f.node()).map_or(f.clone(), |n| {
                    if f.orientation() == og_topo::Orientation::Reversed {
                        n.reversed()
                    } else {
                        n.clone()
                    }
                })
            })
            .collect();
        let rebuilt = make_shell(model, &faces)?.shape;
        if !is_shell_closed(model, &rebuilt)? {
            og_bail!(
                Construction,
                "re-anchoring left the shell open; the shape resists this \
                 repair"
            );
        }
        history.modify(&shell, rebuilt.clone());
        shells.push(rebuilt);
    }
    let solid = make_solid(model, &shells)?.shape;
    history.modify(shape, solid.clone());
    Ok(Built::new(solid, history))
}

/// The circle parameter of a point on a circle curve.
fn circle_parameter(curve: &Curve, p: Point) -> Option<f64> {
    let Curve::Circle(c) = curve else {
        return None;
    };
    let local = c.circle().frame().to_local(p);
    Some(local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU))
}

/// Rebuild a healed face: both rings now share an angle, so the seam the
/// reader could not synthesise can exist.
#[allow(clippy::type_complexity)]
fn rebuild_broken_face(
    model: &mut Model,
    surface_id: og_topo::SurfaceId,
    surface: &SurfaceGeometry,
    rings: &[(Shape, Shape, Curve, (f64, f64))],
    substitution: &HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgResult<Shape> {
    let resolved: Vec<Shape> = rings
        .iter()
        .map(|(e, ..)| substitution.get(&e.node()).unwrap_or(e).clone())
        .collect();
    let point_of = |model: &Model, edge: &Shape| -> OgResult<Point> {
        let Some((a, _)) = og_algo::edge_vertices(model, edge)? else {
            og_bail!(Construction, "a ring edge lost its vertex");
        };
        let Some(node) = model.node(&a) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            og_bail!(Construction, "vertex node holds no vertex data");
        };
        Ok(data.point)
    };
    let pl = point_of(model, &resolved[0])?;
    let ph = point_of(model, &resolved[1])?;
    let (ua, va) = og_algo::project_on_surface(surface, pl, 32, tol)?.parameters;
    let (_, vb) = og_algo::project_on_surface(surface, ph, 32, tol)?.parameters;
    let ((dlo, dhi), _) = surface.domain();
    let span = dhi - dlo;

    let Some(seam_curve) = surface_iso_u_curve(surface, ua, tol) else {
        og_bail!(
            Construction,
            "the healed surface's iso-curve has no closed form; no seam can \
             be built on it"
        );
    };
    let v_lo = og_algo::edge_vertices(model, &resolved[0])?
        .map(|(a, _)| a)
        .ok_or_else(|| og_core::og_err!(Construction, "a ring edge lost its vertex"))?;
    let v_hi = og_algo::edge_vertices(model, &resolved[1])?
        .map(|(a, _)| a)
        .ok_or_else(|| og_core::og_err!(Construction, "a ring edge lost its vertex"))?;
    let (range, from, to, downward) = if va <= vb {
        ((va, vb), v_lo, v_hi, false)
    } else {
        ((vb, va), v_hi, v_lo, true)
    };
    let seam = make_edge_between(model, seam_curve, range, &from, &to, tol)?.shape;
    let side = |u: f64| -> OgResult<PlanarCurve> {
        Ok(og_geom::Line2d::over(
            og_math::Axis2::new(
                og_math::Point2::new(u, 0.0),
                og_math::Direction2::new(og_math::Vector2::new(0.0, 1.0), tol)?,
            ),
            range.0 - 1.0,
            range.1 + 1.0,
        )?
        .into())
    };
    attach_seam(
        model,
        &seam,
        side(ua)?,
        side(ua + span)?,
        surface_id,
        Location::identity(),
        range,
    )?;
    attach_face_pcurves(model, surface, surface_id, &resolved, tol)?;

    let up = if downward {
        seam.reversed()
    } else {
        seam.clone()
    };
    let ring = vec![
        resolved[0].clone(),
        up.clone(),
        resolved[1].reversed(),
        up.reversed(),
    ];
    let wire = make_wire(model, &ring, tol)?.shape;
    Ok(make_face_on(model, surface_id, &[wire], tol)?.shape)
}

/// Rebuild an ordinary face around a substituted edge, pcurves recomputed.
fn rebuild_plain_face(
    model: &mut Model,
    face: &Shape,
    substitution: &HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgResult<Shape> {
    let (surface_id, surface) = {
        let Some(node) = model.node(face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            og_bail!(Dangling, "face refers to a surface not in this model");
        };
        (data.surface, surface.clone())
    };
    let mut wires = Vec::new();
    let mut replaced = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let mut edges = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            match substitution.get(&edge.node()) {
                Some(new_edge) => {
                    let placed = if edge.orientation() == og_topo::Orientation::Reversed {
                        new_edge.reversed()
                    } else {
                        new_edge.clone()
                    };
                    replaced.push(new_edge.clone());
                    edges.push(placed);
                }
                None => edges.push(edge),
            }
        }
        wires.push(make_wire(model, &edges, tol)?.shape);
    }
    attach_face_pcurves(model, &surface, surface_id, &replaced, tol)?;
    Ok(make_face_on(model, surface_id, &wires, tol)?.shape)
}

/// Attach this face's pcurves to freshly rebuilt edges, where the projection
/// has a closed form.
fn attach_face_pcurves(
    model: &mut Model,
    surface: &SurfaceGeometry,
    surface_id: og_topo::SurfaceId,
    edges: &[Shape],
    tol: Tolerances,
) -> OgResult<()> {
    for edge in edges {
        let (curve, range, already) = {
            let Some(node) = model.node(edge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let already = data
                .representations
                .iter()
                .any(|r| matches!(r, EdgeRepr::PCurve { surface: s, .. } if *s == surface_id));
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range, already)
        };
        if already {
            continue;
        }
        let Some(pcurve) = og_intersect::exact_pcurve_of(&curve, surface, tol) else {
            continue;
        };
        // A line pcurve's stated domain must cover the edge's range, which
        // for a wrapped circle runs past one period.
        let pcurve = if let PlanarCurve::Line(l) = &pcurve {
            let (lo, hi) = (l.domain().0.min(range.0), l.domain().1.max(range.1));
            og_geom::Line2d::over(l.axis(), lo, hi)
                .map(Into::into)
                .unwrap_or(pcurve)
        } else {
            pcurve
        };
        og_algo::attach_pcurve(model, edge, pcurve, surface_id, Location::identity(), range)?;
    }
    Ok(())
}
