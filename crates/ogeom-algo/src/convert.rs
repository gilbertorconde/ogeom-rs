//! Whole-shape conversion: everything as NURBS, and general affine
//! transforms over a shape.
//!
//! The geometry-level conversions are exact and were done long ago; what a
//! *shape* needs is the operator that walks it and restates every dependent
//! description, because conversion necessarily reparameterizes — a circle's
//! parameter is its angle and a rational quadratic's is not. So every edge's
//! range moves to its converted curve's domain, and every pcurve is
//! re-derived against a surface whose parameterization has also moved.
//! Re-deriving a pcurve is a *fit* — at the edge's own parameters, so
//! same-parameter holds by construction — which is why this waited for the
//! adaptive fitting core.
//!
//! The affine operator rides on top: an affine map moves a B-spline's
//! control points and nothing else, so the parameterizations of the
//! converted shape survive the transform untouched — including every fitted
//! pcurve, which is the point of converting first.

use crate::build::{attach_pcurve, attach_seam, make_edge_between, make_face_on};
use crate::{Built, History, make_shell, make_solid, make_vertex, make_wire};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{GeneralTransform, Point, Point2, Weighted};
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};
use std::collections::HashMap;

/// One converted edge: the shape, its curve and the curve's own range.
type ConvertedEdge = (Shape, Curve, (f64, f64));

/// A surface's chart window: the `u` then `v` interval.
type ChartWindow = ((f64, f64), (f64, f64));

/// Rebuild a solid with every surface and curve in B-spline form.
///
/// The result is a new solid in world coordinates — every occurrence
/// placement baked in — whose geometry is exactly the original's wherever
/// the conversions are exact (everywhere but the fitted pcurves, whose error
/// is bounded by the fit target derived from the tolerance). History records
/// each original face modified into its converted twin.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape is not a solid or a conversion has no exact form (a trimmed
/// surface's basis is converted; nothing else refuses today);
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if a pcurve
/// refit cannot reach its target.
pub fn to_nurbs(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    rebuild(model, shape, None, true, tol)
}

/// Rebuild a solid with its placements baked into the geometry, keeping
/// every surface and curve in its own analytic vocabulary.
///
/// A uniform-scale placement carries a cylinder to a cylinder and a line to
/// a line — `Transformable` states each exactly — but the stored pcurves
/// describe the *old* parameterizations. This rebuild restates the geometry
/// in world space and re-derives every pcurve against it, exact where the
/// chart alignment allows and projection-fitted with recorded slop where
/// not, which is what lets a scaled shape enter the boolean's analytic
/// pipeline instead of refusing.
///
/// # Errors
///
/// As [`to_nurbs`].
pub fn baked_shape(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Built> {
    rebuild(model, shape, None, false, tol)
}

/// Rebuild a solid under a general affine transform.
///
/// A shear or an uneven scale is not a placement: it carries a circle to an
/// ellipse and a sphere to something with no analytic name here, so the
/// shape is converted to its exact B-spline form and the *control points*
/// are moved — which an affine map does exactly. The pcurves fitted during
/// conversion survive untouched, because an affine map does not
/// reparameterize.
///
/// # Errors
///
/// As [`to_nurbs`].
pub fn general_transformed_shape(
    model: &mut Model,
    shape: &Shape,
    transform: &GeneralTransform,
    tol: Tolerances,
) -> OgeomResult<Built> {
    rebuild(model, shape, Some(transform), true, tol)
}

/// The shared engine: convert, refit, optionally move control points.
fn rebuild(
    model: &mut Model,
    shape: &Shape,
    affine: Option<&GeneralTransform>,
    convert: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(shape)? != ShapeType::Solid {
        ogeom_bail!(Construction, "whole-shape conversion rebuilds solids");
    }
    let map = |p: Point| affine.map_or(p, |t| t.apply(p));
    // The fit target for re-derived pcurves: comfortably inside the model's
    // own working band.
    let target = tol.confusion() * 1e2;

    let mut history = History::new();
    let mut new_vertices: HashMap<(TShapeId, [u64; 3]), Shape> = HashMap::new();
    let mut new_edges: HashMap<(TShapeId, [u64; 3]), ConvertedEdge> = HashMap::new();
    let mut shells = Vec::new();

    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        let mut faces = Vec::new();
        for face in model.ordered_children_of(&shell)? {
            let placement = face.transform(model.datums())?;
            let (surface_id_old, old_surface) = {
                let Some(node) = model.node(&face) else {
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
            use ogeom_geom::Transformable as _;
            let placed = old_surface.transformed(&placement, tol)?;
            // Bounded to the face's own chart region first: a plane declares
            // an enormous domain, and a patch over all of it would map the
            // face to a dot in its chart.
            let (placed, old_window) = bounded_to_face(
                model,
                &face,
                surface_id_old,
                &placed,
                placement.scale_factor().abs(),
                tol,
            )?;
            let patch_surface: SurfaceGeometry = if convert {
                let mut patch = placed.to_bspline(tol)?;
                if let Some(t) = affine {
                    patch = transformed_patch(&patch, t)?;
                }
                patch.into()
            } else {
                placed.clone()
            };
            let surface_id = model.geometry_mut().add_surface(patch_surface.clone());

            let mut wires = Vec::new();
            let mut corner_uv: HashMap<TShapeId, Point2> = HashMap::new();
            for wire in model.ordered_children_of(&face)? {
                let mut ring = Vec::new();
                let mut seams_done: Vec<TShapeId> = Vec::new();
                for edge in model.ordered_children_of(&wire)? {
                    let edge_placement = edge.transform(model.datums())?;
                    let key = (edge.node(), placement_bits(&edge_placement));
                    let data = {
                        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                            ogeom_bail!(Construction, "edge node holds no edge data");
                        };
                        data.clone()
                    };
                    if data.degenerate {
                        // A pole or an apex: still an edge in parameter
                        // space; its pcurve is rebuilt below from the old
                        // chart row's place in the new chart.
                        let Some(vertex) = model.children_of(&edge)?.first().cloned() else {
                            ogeom_bail!(Construction, "a degenerate edge has no vertex");
                        };
                        let at = {
                            let Some(v) = model.node(&vertex).and_then(|n| n.data().as_vertex())
                            else {
                                ogeom_bail!(Construction, "vertex node holds no vertex data");
                            };
                            map(edge_placement.apply(v.point))
                        };
                        let new_vertex = new_vertices
                            .entry((vertex.node(), point_bits(at)))
                            .or_insert_with(|| make_vertex(model, at).shape)
                            .clone();
                        let mut degenerate = ogeom_topo::EdgeData::new();
                        degenerate.degenerate = true;
                        let new_edge =
                            model.add_edge(degenerate, &[new_vertex.clone(), new_vertex])?;
                        let row = degenerate_row(
                            model,
                            &data,
                            surface_id_old,
                            &edge,
                            &placed,
                            &patch_surface,
                            tol,
                        )?;
                        attach_pcurve(
                            model,
                            &new_edge,
                            row.into(),
                            surface_id,
                            Location::identity(),
                            row.domain(),
                        )?;
                        ring.push(if edge.orientation() == Orientation::Reversed {
                            new_edge.reversed()
                        } else {
                            new_edge
                        });
                        continue;
                    }

                    let (new_edge, new_curve, new_range) = match new_edges.get(&key) {
                        Some(found) => found.clone(),
                        None => {
                            let built = convert_edge(
                                model,
                                &edge,
                                &data,
                                &map,
                                convert,
                                &mut new_vertices,
                                tol,
                            )?;
                            history.modify(&edge, built.0.clone());
                            new_edges.insert(key, built.clone());
                            built
                        }
                    };

                    // The pcurve on this face: fitted at the new edge's own
                    // parameters, seam sides each fitted against their own
                    // half of the chart.
                    let old_repr = data.pcurve_for(surface_id_old, edge.location()).cloned();
                    let was_seam = matches!(old_repr, Some(EdgeRepr::Seam { .. }));
                    if was_seam {
                        if !seams_done.contains(&new_edge.node()) {
                            seams_done.push(new_edge.node());
                            let (forward, reversed) =
                                exact_seam_columns(&new_curve, new_range, &patch_surface, tol)
                                    .map_or_else(
                                        || {
                                            seam_pcurves(
                                                &new_curve,
                                                new_range,
                                                &patch_surface,
                                                target,
                                                tol,
                                            )
                                        },
                                        Ok,
                                    )?;
                            let seam_range = forward.domain();
                            attach_seam(
                                model,
                                &new_edge,
                                forward,
                                reversed,
                                surface_id,
                                Location::identity(),
                                seam_range,
                            )?;
                        }
                    } else {
                        // The exact path: an edge that ran along the old
                        // chart's iso direction runs along the new chart's,
                        // and the boundary conversion shares the patch
                        // direction's parameterization by construction.
                        if let Some(iso) = exact_iso_pcurve(
                            model,
                            old_repr.as_ref(),
                            old_window,
                            &new_curve,
                            new_range,
                            &patch_surface,
                            tol,
                        )? {
                            let iso_range = iso.domain();
                            attach_pcurve(
                                model,
                                &new_edge,
                                iso,
                                surface_id,
                                Location::identity(),
                                iso_range,
                            )?;
                            ring.push(if edge.orientation() == Orientation::Reversed {
                                new_edge.reversed()
                            } else {
                                new_edge
                            });
                            continue;
                        }
                        let uv_of = |model: &Model,
                                     cache: &mut HashMap<TShapeId, Point2>,
                                     vertex: &Shape|
                         -> OgeomResult<Point2> {
                            if let Some(&uv) = cache.get(&vertex.node()) {
                                return Ok(uv);
                            }
                            let Some(v) = model.node(vertex).and_then(|n| n.data().as_vertex())
                            else {
                                ogeom_bail!(Construction, "vertex node holds no vertex data");
                            };
                            let projection = crate::measure::project_on_surface(
                                &patch_surface,
                                v.point,
                                24,
                                tol,
                            )?;
                            let uv = Point2::new(projection.parameters.0, projection.parameters.1);
                            cache.insert(vertex.node(), uv);
                            Ok(uv)
                        };
                        let bounds = model.children_of(&new_edge)?;
                        // A closed edge — a rim — has one vertex at both
                        // ends, and near a seam its single chart image is
                        // one side's; pinning both ends there would fold the
                        // ring. Its trace closes on its own.
                        let closed_edge =
                            bounds.len() < 2 || bounds[0].node() == bounds[bounds.len() - 1].node();
                        let ends = if closed_edge {
                            (None, None)
                        } else {
                            (
                                Some(uv_of(model, &mut corner_uv, &bounds[0])?),
                                Some(uv_of(model, &mut corner_uv, &bounds[bounds.len() - 1])?),
                            )
                        };
                        let fitted = fit_pcurve(
                            &new_curve,
                            new_range,
                            &patch_surface,
                            None,
                            ends,
                            target,
                            tol,
                        )?;
                        attach_pcurve(
                            model,
                            &new_edge,
                            fitted,
                            surface_id,
                            Location::identity(),
                            new_range,
                        )?;
                    }
                    ring.push(if edge.orientation() == Orientation::Reversed {
                        new_edge.reversed()
                    } else {
                        new_edge
                    });
                }
                wires.push(make_wire(model, &ring, tol)?.shape);
            }
            let built = make_face_on(model, surface_id, &wires, tol)?.shape;
            let built = if face.orientation() == Orientation::Reversed {
                built.reversed()
            } else {
                built
            };
            history.modify(&face, built.clone());
            faces.push(built);
        }
        shells.push(make_shell(model, &faces)?.shape);
    }
    let solid = make_solid(model, &shells)?.shape;
    history.modify(shape, solid.clone());
    Ok(Built::new(solid, history))
}

/// One edge converted: exact B-spline over its range, vertices carried.
fn convert_edge(
    model: &mut Model,
    edge: &Shape,
    data: &ogeom_topo::EdgeData,
    map: &dyn Fn(Point) -> Point,
    convert: bool,
    vertices: &mut HashMap<(TShapeId, [u64; 3]), Shape>,
    tol: Tolerances,
) -> OgeomResult<(Shape, Curve, (f64, f64))> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "an edge with no curve cannot be converted");
    };
    let Some(geometry) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placement = edge.transform(model.datums())?;
    use ogeom_geom::Transformable as _;
    let placed = geometry.transformed(&placement, tol)?;
    // A placement with scale rescales a length-parameterized curve's domain
    // — the same rule `transformed` itself applies — so the edge's range
    // must move with it or the conversion covers only part of the edge.
    let stretch = placement.scale_factor().abs();
    let range_on_placed = match &geometry {
        Curve::Line(_) => (range.0 * stretch, range.1 * stretch),
        Curve::Trimmed(t) if matches!(t.basis(), Curve::Line(_)) => {
            (range.0 * stretch, range.1 * stretch)
        }
        _ => *range,
    };
    let (curve, new_range): (Curve, (f64, f64)) = if convert {
        let spline = placed.to_bspline_over(range_on_placed, tol)?;
        // The affine map moves control points; the parameterization and the
        // weights stay, which is the whole point of converting first.
        let weighted: Vec<Weighted<Point>> = spline
            .control_points()
            .iter()
            .map(|c| {
                let p = Point::from_vector(c.scaled.to_vector() / c.weight);
                Weighted::new(map(p), c.weight, tol)
            })
            .collect::<OgeomResult<_>>()?;
        let spline = ogeom_geom::BSplineCurve::rational(spline.knots().clone(), weighted)?;
        let curve: Curve = spline.into();
        let range = curve.domain();
        (curve, range)
    } else {
        (placed, range_on_placed)
    };
    // Vertices are shared across every edge that meets them: cached by the
    // old node and the old vertex's own mapped point — the same bits every
    // neighbouring edge computes, unlike each spline's evaluated end.
    let old = model.children_of(edge)?;
    if old.is_empty() {
        ogeom_bail!(Construction, "an edge with no vertices cannot be converted");
    }
    let mut resolve = |model: &mut Model, occurrence: &Shape| -> OgeomResult<Shape> {
        let Some(data) = model.node(occurrence).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let at = map(placement.apply(data.point));
        Ok(vertices
            .entry((occurrence.node(), point_bits(at)))
            .or_insert_with(|| make_vertex(model, at).shape)
            .clone())
    };
    let va = resolve(model, &old[0])?;
    let vb = if old.len() == 1 {
        va.clone()
    } else {
        resolve(model, &old[old.len() - 1])?
    };
    let built = make_edge_between(model, curve.clone(), new_range, &va, &vb, tol)?.shape;
    Ok((built, curve, new_range))
}

/// A pcurve refit at the edge's own parameters, optionally biased toward one
/// side of the chart's `u` seam.
#[allow(clippy::too_many_arguments)]
fn fit_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    side: Option<f64>,
    ends: (Option<Point2>, Option<Point2>),
    target: f64,
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::PlanarCurve> {
    const SAMPLES: usize = 48;
    let ((u0, u1), _) = surface.domain();
    let mut parameters = Vec::with_capacity(SAMPLES + 1);
    let mut trace = Vec::with_capacity(SAMPLES + 1);
    for i in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = range.0 + (range.1 - range.0) * i as f64 / SAMPLES as f64;
        let p = curve.point_at(t, tol)?;
        let projection = crate::measure::project_on_surface(surface, p, 24, tol)?;
        let (mut u, v) = projection.parameters;
        if let Some(bias) = side {
            // A point exactly on the seam projects to either side; the bias
            // names which column this pcurve is.
            if (u - u0).abs() < (u1 - u0) * 0.25 || (u - u1).abs() < (u1 - u0) * 0.25 {
                u = bias;
            }
        }
        parameters.push(t);
        trace.push(Point2::new(u, v));
    }
    // A closed patch is clamped, not periodic, but its geometry still closes:
    // a trace crossing the closure jumps by the whole u span. Unwrapped, the
    // fit sees the continuous curve the edge is; the pcurve may leave [u0,u1]
    // by a little at the join, which evaluation tolerates.
    let span = u1 - u0;
    for i in 1..trace.len() {
        let mut du = trace[i].x - trace[i - 1].x;
        while du > span * 0.5 {
            trace[i].x -= span;
            du = trace[i].x - trace[i - 1].x;
        }
        while du < -span * 0.5 {
            trace[i].x += span;
            du = trace[i].x - trace[i - 1].x;
        }
    }
    // Ring corners are shared property: every edge meeting a vertex uses the
    // vertex's own chart image, projected once per face, so the ring closes
    // exactly instead of to each fit's own rounding.
    if let Some(start) = ends.0
        && let Some(first) = trace.first_mut()
    {
        *first = start;
    }
    if let Some(end) = ends.1
        && let Some(last) = trace.last_mut()
    {
        *last = end;
    }
    let fitted = ogeom_geom::fit::fit_points_2d_at(&parameters, &trace, 3, target, tol)?;
    Ok(fitted.curve.into())
}

/// An exact iso pcurve in the new chart, when the old one was an iso line.
///
/// Both endpoints project onto the patch; the pcurve is the straight chart
/// segment between them, exact because the boundary conversion and the patch
/// direction share one parameterization. `None` for anything that was not an
/// axis-aligned chart line, which takes the fitted path.
fn exact_iso_pcurve(
    model: &Model,
    old_repr: Option<&EdgeRepr>,
    _old_window: ChartWindow,
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::PlanarCurve>> {
    let Some(EdgeRepr::PCurve { curve: old_pc, .. }) = old_repr else {
        return Ok(None);
    };
    let Some(ogeom_geom::PlanarCurve::Line(old_line)) = model.geometry().pcurve(*old_pc) else {
        return Ok(None);
    };
    let d = old_line.axis().direction;
    let axis_aligned = d.vector().x.abs() <= 1e-12 || d.vector().y.abs() <= 1e-12;
    if !axis_aligned {
        return Ok(None);
    }
    let ((nu0, nu1), (nv0, nv1)) = surface.domain();
    let start = curve.point_at(range.0, tol)?;
    let end = curve.point_at(range.1, tol)?;
    let a = crate::measure::project_on_surface(surface, start, 24, tol)?;
    let b = crate::measure::project_on_surface(surface, end, 24, tol)?;
    let mut a = Point2::new(a.parameters.0, a.parameters.1);
    let mut b = Point2::new(b.parameters.0, b.parameters.1);
    let row = d.vector().y.abs() <= 1e-12;
    if row {
        // A row's v is one value; the endpoints agree up to projection noise
        // — and a full-width row on a closed chart runs edge to edge, its
        // seam-side endpoints snapped to the chart's own edges.
        let v = f64::midpoint(a.y, b.y);
        a.y = v;
        b.y = v;
        let span = nu1 - nu0;
        for u in [&mut a.x, &mut b.x] {
            if (*u - nu0).abs() < span * 1e-6 || (*u - nu1).abs() < span * 1e-6 {
                // On the closure both images are the same point; the edge's
                // own direction decides which end this is.
            }
        }
        if (a.x - b.x).abs() < span * 1e-9 {
            // A closed row: the full width, oriented by the curve's start.
            a.x = nu0;
            b.x = nu1;
        }
    } else {
        let u = f64::midpoint(a.x, b.x);
        a.x = u;
        b.x = u;
        if (a.y - b.y).abs() < (nv1 - nv0) * 1e-9 {
            a.y = nv0;
            b.y = nv1;
        }
    }
    Ok(Some(ogeom_geom::Line2d::segment(a, b, tol)?.into()))
}

/// Exact seam columns for a closed patch: the chart's two `u` edges over the
/// seam's projected `v` span. `None` when the projections cannot say.
fn exact_seam_columns(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<(ogeom_geom::PlanarCurve, ogeom_geom::PlanarCurve)> {
    let ((u0, u1), _) = {
        use ogeom_geom::Surface as _;
        surface.domain()
    };
    let start = curve.point_at(range.0, tol).ok()?;
    let end = curve.point_at(range.1, tol).ok()?;
    let a = crate::measure::project_on_surface(surface, start, 24, tol).ok()?;
    let b = crate::measure::project_on_surface(surface, end, 24, tol).ok()?;
    let (va, vb) = (a.parameters.1, b.parameters.1);
    let forward: ogeom_geom::PlanarCurve =
        ogeom_geom::Line2d::segment(Point2::new(u1, va), Point2::new(u1, vb), tol)
            .ok()?
            .into();
    let reversed: ogeom_geom::PlanarCurve =
        ogeom_geom::Line2d::segment(Point2::new(u0, va), Point2::new(u0, vb), tol)
            .ok()?
            .into();
    Some((forward, reversed))
}

/// The two seam-side pcurves of a closed patch: the chart's first and last
/// `u` columns, each fitted at the edge's own parameters.
fn seam_pcurves(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    target: f64,
    tol: Tolerances,
) -> OgeomResult<(ogeom_geom::PlanarCurve, ogeom_geom::PlanarCurve)> {
    let ((u0, u1), _) = surface.domain();
    let forward = fit_pcurve(curve, range, surface, Some(u1), (None, None), target, tol)?;
    let reversed = fit_pcurve(curve, range, surface, Some(u0), (None, None), target, tol)?;
    Ok((forward, reversed))
}

/// The pole row of a degenerate edge, restated in the new chart.
fn degenerate_row(
    model: &Model,
    data: &ogeom_topo::EdgeData,
    old_surface: ogeom_topo::SurfaceId,
    edge: &Shape,
    old: &SurfaceGeometry,
    new: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::Line2d> {
    // Which end of the old chart the pole row sat at decides which end of
    // the new chart it sits at: conversion preserves ends.
    let Some(EdgeRepr::PCurve { curve, .. }) = data.pcurve_for(old_surface, edge.location()) else {
        ogeom_bail!(Construction, "a degenerate edge has no pcurve to restate");
    };
    let Some(row) = model.geometry().pcurve(*curve) else {
        ogeom_bail!(Dangling, "pcurve is not in this model");
    };
    let (lo, hi) = row.domain();
    let v_old = f64::midpoint(row.point_at(lo, tol)?.y, row.point_at(hi, tol)?.y);
    let (_, (va_old, vb_old)) = old.domain();
    let ((nu0, nu1), (nv0, nv1)) = new.domain();
    let v_new = if (v_old - va_old).abs() <= (v_old - vb_old).abs() {
        nv0
    } else {
        nv1
    };
    ogeom_geom::Line2d::segment(Point2::new(nu0, v_new), Point2::new(nu1, v_new), tol)
}

/// The surface shrunk to the face's own chart region, with a margin.
///
/// The face's pcurves say which part of the surface the face actually uses;
/// converting the whole declared domain would spend the patch's parameter
/// range on empty plane. Kinds whose windows are structural — a sphere's, a
/// torus's, the closed direction of a cylinder — keep them.
fn bounded_to_face(
    model: &Model,
    face: &Shape,
    surface_id: ogeom_topo::SurfaceId,
    placed: &SurfaceGeometry,
    stretch: f64,
    tol: Tolerances,
) -> OgeomResult<(SurfaceGeometry, ChartWindow)> {
    let (mut u0, mut u1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY);
    for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let mut widen = |pc: ogeom_topo::PCurveId, range: (f64, f64)| -> OgeomResult<()> {
            let Some(pcurve) = model.geometry().pcurve(pc) else {
                return Ok(());
            };
            for i in 0..=8 {
                #[allow(clippy::cast_precision_loss)]
                let t = range.0 + (range.1 - range.0) * f64::from(i) / 8.0;
                let p = pcurve.point_at(t, tol)?;
                u0 = u0.min(p.x);
                u1 = u1.max(p.x);
                v0 = v0.min(p.y);
                v1 = v1.max(p.y);
            }
            Ok(())
        };
        match data.pcurve_for(surface_id, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => widen(*curve, *range)?,
            Some(EdgeRepr::Seam {
                forward,
                reversed,
                range,
                ..
            }) => {
                widen(*forward, *range)?;
                widen(*reversed, *range)?;
            }
            _ => {}
        }
    }
    if !u0.is_finite() || u1 - u0 <= tol.confusion() || v1 - v0 <= tol.confusion() {
        return Ok((placed.clone(), placed.domain()));
    }
    // The window was read off the *old* pcurves, in the old chart's units;
    // a placement with scale stretches the placed chart's length-like
    // directions, and the trim must stretch with them. Angle directions
    // never stretch.
    let (su, sv) = match placed {
        SurfaceGeometry::Plane(_) => (stretch, stretch),
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => (1.0, stretch),
        _ => (1.0, 1.0),
    };
    let (u0, u1, v0, v1) = (u0 * su, u1 * su, v0 * sv, v1 * sv);
    let margin = ((u1 - u0) + (v1 - v0)) * 0.05 + tol.confusion();
    let bounded: SurfaceGeometry = match placed {
        SurfaceGeometry::Plane(p) => ogeom_geom::PlaneSurface::over(
            p.plane(),
            (u0 - margin, u1 + margin),
            (v0 - margin, v1 + margin),
        )?
        .into(),
        SurfaceGeometry::Cylinder(c) => {
            ogeom_geom::CylinderSurface::new(c.cylinder(), (v0 - margin, v1 + margin))?.into()
        }
        SurfaceGeometry::Cone(c) => {
            ogeom_geom::ConeSurface::new(c.cone(), (v0 - margin, v1 + margin))?.into()
        }
        other => other.clone(),
    };
    let window = bounded.domain();
    Ok((bounded, window))
}

/// A B-spline patch with its control points carried through an affine map.
fn transformed_patch(
    patch: &ogeom_geom::BSplineSurface,
    transform: &GeneralTransform,
) -> OgeomResult<ogeom_geom::BSplineSurface> {
    let grid = patch.grid();
    let mapped = grid.map(|c| {
        let p = Point::from_vector(c.scaled.to_vector() / c.weight);
        let moved = transform.apply(p);
        Weighted {
            scaled: Point::from_vector(moved.to_vector() * c.weight),
            weight: c.weight,
        }
    });
    ogeom_geom::BSplineSurface::rational(patch.u_knots().clone(), patch.v_knots().clone(), mapped)
}

/// A rigid placement quantized for deduplication keys.
fn placement_bits(t: &ogeom_math::Transform) -> [u64; 3] {
    let p = t.apply(Point::new(0.123_456_789, 9.87, -3.21));
    point_bits(p)
}

fn point_bits(p: Point) -> [u64; 3] {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::{make_box, make_cylinder, volume_properties};
    use ogeom_math::Frame;
    use ogeom_mesh::Deflection;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_box_converts_to_nurbs_and_measures_the_same() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let converted = to_nurbs(&mut model, &solid.shape, T).unwrap();

        // Every surface is now a spline.
        for face in explore(&model, &converted.shape, Filter::OfType(ShapeType::Face)).unwrap() {
            let NodeData::Face(data) = model.node(&face).unwrap().data() else {
                panic!("face data");
            };
            assert!(matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::BSpline(_))
            ));
        }
        let props = volume_properties(&model, &converted.shape, fine(), T).unwrap();
        assert!((props.mass - 24.0).abs() < 1e-6, "volume {}", props.mass);
        assert!(
            crate::check(&model, &converted.shape, T)
                .unwrap()
                .is_valid(),
            "{}",
            crate::check(&model, &converted.shape, T).unwrap()
        );
        assert_eq!(
            converted.history.modified(&solid.shape),
            std::slice::from_ref(&converted.shape)
        );
    }

    #[test]
    fn a_cylinder_converts_seams_poles_and_all() {
        let mut model = Model::new();
        let solid = make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        let converted = to_nurbs(&mut model, &solid.shape, T).unwrap();
        let exact = core::f64::consts::PI * 4.0 * 5.0;
        let props = volume_properties(&model, &converted.shape, fine(), T).unwrap();
        assert!(
            (props.mass - exact).abs() < exact * 5e-3,
            "volume {} against {exact}",
            props.mass
        );
    }

    #[test]
    fn a_shear_preserves_volume_and_a_stretch_scales_it() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();

        // A unit-determinant shear: x' = x + 0.5 y.
        let shear = GeneralTransform::new(
            ogeom_math::Matrix3 {
                rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            ogeom_math::Vector::ZERO,
        );
        let sheared = general_transformed_shape(&mut model, &solid.shape, &shear, T).unwrap();
        let props = volume_properties(&model, &sheared.shape, fine(), T).unwrap();
        assert!((props.mass - 24.0).abs() < 1e-6, "sheared {}", props.mass);

        // An uneven stretch doubles x: volume doubles.
        let stretch = GeneralTransform::new(
            ogeom_math::Matrix3 {
                rows: [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            ogeom_math::Vector::ZERO,
        );
        let stretched = general_transformed_shape(&mut model, &solid.shape, &stretch, T).unwrap();
        let props = volume_properties(&model, &stretched.shape, fine(), T).unwrap();
        assert!((props.mass - 48.0).abs() < 1e-6, "stretched {}", props.mass);
    }

    #[test]
    fn a_sheared_cylinder_still_encloses_its_volume() {
        // The case placements cannot express: the circular rims become
        // ellipses, and the volume is invariant under the unit-det shear.
        let mut model = Model::new();
        let solid = make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        let shear = GeneralTransform::new(
            ogeom_math::Matrix3 {
                rows: [[1.0, 0.0, 0.4], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            ogeom_math::Vector::ZERO,
        );
        let sheared = general_transformed_shape(&mut model, &solid.shape, &shear, T).unwrap();
        let exact = core::f64::consts::PI * 4.0 * 5.0;
        let props = volume_properties(&model, &sheared.shape, fine(), T).unwrap();
        assert!(
            (props.mass - exact).abs() < exact * 5e-3,
            "sheared cylinder {} against {exact}",
            props.mass
        );
    }
}
