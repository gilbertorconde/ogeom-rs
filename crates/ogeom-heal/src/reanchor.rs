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

use ogeom_algo::{
    Built, History, is_shell_closed, make_edge_between, make_face_on, make_shell, make_solid,
    make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry};
use ogeom_math::Point;
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};
use std::collections::HashMap;

/// A face bounded by two closed rings on a periodic surface, and what the
/// repair needs to know about it.
struct Broken {
    face: Shape,
    surface_id: ogeom_topo::SurfaceId,
    surface: SurfaceGeometry,
    /// Each ring: its single closed edge, that edge's vertex, curve and
    /// range.
    rings: Vec<(Shape, Shape, Curve, (f64, f64))>,
}

/// A coaxial chain of broken faces, healed with one shared half-plane.
struct Chain {
    point: Point,
    dir: ogeom_math::Vector,
    members: Vec<usize>,
}

/// Re-anchor misaligned ring vertices so every periodic face can carry a seam.
///
/// Faces already whole pass through untouched; a shape with nothing to heal
/// comes back as itself with an empty history.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the shape is
/// not a solid, or a broken face's structure resists the repair — in which
/// case nothing is modified.
pub fn reanchor_periodic_rings(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    if model.kind_of(shape)? != ShapeType::Solid {
        ogeom_bail!(Construction, "ring re-anchoring heals solids");
    }

    // --- find the broken faces ---------------------------------------------
    let mut revolution: Vec<Broken> = Vec::new();
    let mut any_broken = false;
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
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
                let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
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
            let Some((a, _)) = ogeom_algo::edge_vertices(model, &edge)? else {
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
        return Ok((Built::new(shape.clone(), History::new()), 0));
    }
    let broken = revolution;

    // --- one half-plane per coaxial chain ------------------------------------
    //
    // Re-anchoring one ring per face is not enough: rings are shared, and a
    // neighbouring face's seam ties to the vertex being moved. For a coaxial
    // chain — which a fillet stack is — a single half-plane through the axis
    // meets every ring exactly once, and anchoring them all there satisfies
    // every face's equal-angle constraint at one stroke. A real part carries
    // several such chains on several axes — every bored hole is its own —
    // so the broken faces are grouped by axis first and each group healed
    // with its own half-plane; a group whose rings resist (a non-circle
    // ring, a stray non-coaxial circle) is left alone rather than damning
    // the rest of the part.
    let mut history = History::new();
    let axis_of = |s: &SurfaceGeometry| -> Option<(Point, ogeom_math::Vector)> {
        match s {
            SurfaceGeometry::Cylinder(c) => {
                let f = c.cylinder().frame();
                Some((f.origin(), f.z().vector()))
            }
            SurfaceGeometry::Cone(c) => {
                let f = c.cone().frame();
                Some((f.origin(), f.z().vector()))
            }
            SurfaceGeometry::Sphere(c) => {
                let f = c.sphere().frame();
                Some((f.origin(), f.z().vector()))
            }
            SurfaceGeometry::Torus(c) => {
                let f = c.torus().frame();
                Some((f.origin(), f.z().vector()))
            }
            _ => None,
        }
    };
    let mut chains: Vec<Chain> = Vec::new();
    for (i, b) in broken.iter().enumerate() {
        let Some((point, dir)) = axis_of(&b.surface) else {
            continue;
        };
        let joined = chains.iter_mut().find(|c| {
            c.dir.cross(dir).magnitude() < 1e-6
                && (point - c.point).cross(c.dir).magnitude() < tol.confusion() * 1e3
        });
        match joined {
            Some(c) => c.members.push(i),
            None => chains.push(Chain {
                point,
                dir,
                members: vec![i],
            }),
        }
    }

    let mut substitution: HashMap<TShapeId, Shape> = HashMap::new();
    let mut healable: Vec<usize> = Vec::new();
    for chain in &chains {
        let anchored = anchor_chain(model, &broken, chain, tol);
        match anchored {
            Ok(subs) => {
                for (node, edge, old_edge, old_vertex) in subs {
                    history.modify(&old_edge, edge.clone());
                    history.delete(&old_vertex);
                    substitution.insert(node, edge);
                }
                healable.extend_from_slice(&chain.members);
            }
            Err(_) => {
                // This chain resists; its faces stay as imported.
            }
        }
    }
    if healable.is_empty() {
        return Ok((Built::new(shape.clone(), History::new()), 0));
    }

    // --- rebuild every face that touches a substituted edge ------------------
    let mut face_map: HashMap<TShapeId, Shape> = HashMap::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        // Every revolution face is rebuilt, moved rings or not: an
        // already-aligned band may still carry import-time pcurves whose
        // chart windows disagree, and the rebuild is what makes its
        // annotation coherent. Plain faces rebuild only when an edge they
        // use actually changed.
        let uses_any = explore_unique(model, &face, ShapeType::Edge)?
            .iter()
            .any(|e| substitution.contains_key(&e.node()));
        let is_revolution = healable.iter().any(|&i| broken[i].face.is_same(&face));
        if !uses_any && !is_revolution {
            continue;
        }
        let rebuilt = if let Some(b) = healable
            .iter()
            .map(|&i| &broken[i])
            .find(|b| b.face.is_same(&face))
        {
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
        let rebuilt = if face.orientation() == ogeom_topo::Orientation::Reversed {
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
                    if f.orientation() == ogeom_topo::Orientation::Reversed {
                        n.reversed()
                    } else {
                        n.clone()
                    }
                })
            })
            .collect();
        let rebuilt = make_shell(model, &faces)?.shape;
        if !is_shell_closed(model, &rebuilt)? {
            ogeom_bail!(
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
    let moved = substitution.len();
    Ok((Built::new(solid, history), moved))
}

/// Anchor every ring of one coaxial chain at the chain's own half-plane.
///
/// Returns the substitutions to apply — `(old node, new edge, old edge, old
/// vertex)` — or an error if any ring in the chain resists, in which case
/// nothing has been decided and the chain is left as imported. New vertices
/// and edges may have been added to the model, but nothing references them.
#[allow(clippy::type_complexity)]
fn anchor_chain(
    model: &mut Model,
    broken: &[Broken],
    chain: &Chain,
    tol: Tolerances,
) -> OgeomResult<Vec<(TShapeId, Shape, Shape, Shape)>> {
    let (axis_point, axis_dir) = (chain.point, chain.dir);
    let radial = {
        let (_, vertex, curve, _) = &broken[chain.members[0]].rings[0];
        let Curve::Circle(c) = curve else {
            ogeom_bail!(
                Construction,
                "a ring to re-anchor is not a circle; the repair does not \
                 know its parameterization"
            );
        };
        let circle = c.circle();
        let Some(node) = model.node(vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let r = data.point - circle.centre();
        let radial = r - axis_dir * r.dot(axis_dir);
        if radial.magnitude() <= tol.confusion() {
            ogeom_bail!(Construction, "an anchor vertex sits on the axis");
        }
        radial / radial.magnitude()
    };

    let mut out = Vec::new();
    let mut done: Vec<TShapeId> = Vec::new();
    for &i in &chain.members {
        for (old_edge, old_vertex, curve, _) in &broken[i].rings {
            if done.contains(&old_edge.node()) {
                continue;
            }
            done.push(old_edge.node());
            let Curve::Circle(c) = curve else {
                ogeom_bail!(
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
                ogeom_bail!(
                    Construction,
                    "the rings are not coaxial; the half-plane repair does \
                     not apply"
                );
            }
            let target = circle.centre() + radial * circle.radius();
            let current = {
                let Some(node) = model.node(old_vertex) else {
                    ogeom_bail!(Dangling, "vertex is not in this model");
                };
                let Some(data) = node.data().as_vertex() else {
                    ogeom_bail!(Construction, "vertex node holds no vertex data");
                };
                data.point
            };
            if current.distance(target) <= tol.confusion() * 1e2 {
                continue;
            }
            let Some(t_star) = circle_parameter(curve, target) else {
                ogeom_bail!(Construction, "an anchor point fell off its circle");
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
            out.push((
                old_edge.node(),
                rebuilt,
                old_edge.clone(),
                old_vertex.clone(),
            ));
        }
    }
    Ok(out)
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
/// reader could not synthesise can exist. The band construction itself is
/// ogeom-algo's `make_revolution_band` — one authority for reader and healer.
#[allow(clippy::type_complexity)]
fn rebuild_broken_face(
    model: &mut Model,
    _old_surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    rings: &[(Shape, Shape, Curve, (f64, f64))],
    substitution: &HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let resolved: Vec<Shape> = rings
        .iter()
        .map(|(e, ..)| substitution.get(&e.node()).unwrap_or(e).clone())
        .collect();
    ogeom_algo::make_revolution_band(model, surface, &resolved[0], &resolved[1], tol)
}

/// Rebuild an ordinary face around a substituted edge, pcurves recomputed.
fn rebuild_plain_face(
    model: &mut Model,
    face: &Shape,
    substitution: &HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (surface_id, surface) = {
        let Some(node) = model.node(face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
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
                    let placed = if edge.orientation() == ogeom_topo::Orientation::Reversed {
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
    surface_id: ogeom_topo::SurfaceId,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<()> {
    for edge in edges {
        let (curve, range, already) = {
            let Some(node) = model.node(edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let already = data
                .representations
                .iter()
                .any(|r| matches!(r, EdgeRepr::PCurve { surface: s, .. } if *s == surface_id));
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range, already)
        };
        if already {
            continue;
        }
        let Some(pcurve) = ogeom_intersect::exact_pcurve_over(&curve, range, surface, tol) else {
            continue;
        };
        // A line pcurve's stated domain must cover the edge's range, which
        // for a wrapped circle runs past one period.
        let pcurve = if let PlanarCurve::Line(l) = &pcurve {
            let (lo, hi) = (l.domain().0.min(range.0), l.domain().1.max(range.1));
            ogeom_geom::Line2d::over(l.axis(), lo, hi)
                .map(Into::into)
                .unwrap_or(pcurve)
        } else {
            pcurve
        };
        ogeom_algo::attach_pcurve(model, edge, pcurve, surface_id, Location::identity(), range)?;
    }
    Ok(())
}
