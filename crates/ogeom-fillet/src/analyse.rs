//! What a blend actually achieved, measured rather than asserted.
//!
//! A blend claims two things: that it meets each of its supports along a
//! curve, and that it meets them *smoothly*: the two surfaces sharing a
//! normal there. Both are claims about geometry that construction can get
//! subtly wrong: a fitted section drifts, a marched spine carries its
//! chord budget, a rebuilt face lands on a neighbour a hair off. This
//! module reports the two numbers instead of trusting them, with a third
//! for surface work that asks more: whether the two faces also bend alike
//! across the edge. Per-point curvature over a face, for curvature maps and
//! zebra stripes, is [`face_curvature_samples`].
//!
//! The report is per shared edge, because that is where the claim lives.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceCurvature};
use ogeom_math::{Direction, Point, Point2, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

/// How well a blend face meets one neighbour along one edge.
#[derive(Debug, Clone, PartialEq)]
pub struct BlendContact {
    /// The shared edge.
    pub edge: Shape,
    /// The face on the other side of it.
    pub neighbour: Shape,
    /// The largest angle, in radians, between the two surfaces' normals at
    /// the sampled stations: zero for a tangent join.
    pub tangency_error: f64,
    /// The largest difference, in inverse model units, between the two
    /// surfaces' normal curvatures in the direction square to the edge at
    /// the sampled stations, signed against one shared sense of the
    /// normal: zero for a curvature-continuous join. It speaks of a join
    /// whose `tangency_error` is small; across a crease the two directions
    /// it compares are not one direction.
    pub curvature_error: f64,
    /// The largest distance, in model units, between the edge's curve and
    /// the two surfaces it is supposed to lie on.
    pub gap: f64,
    /// How many stations were sampled along the edge.
    pub stations: usize,
}

/// Measure a blend face against every face it shares an edge with.
///
/// Every number is sampled: `stations` points spread over each shared
/// edge's range, the surfaces read at the parameters the edge's own pcurves
/// give, so nothing is inverted and nothing is guessed. A blend that is
/// tangent everywhere except between two stations reports zero, which is
/// the honest limit of sampling and the reason the count is in the report.
///
/// An edge whose pcurve on either face is missing cannot be measured this
/// way at all: it is reported with its `gap` and `tangency_error` set to
/// infinity rather than quietly skipped, because a blend nobody can measure
/// is not a blend anybody should trust. Its `curvature_error` is infinite
/// on the same grounds, and so is that of an edge where no station gave a
/// curvature on both sides.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `blend` is not a face of `shape`, or `stations` is less than two.
pub fn analyse_blend(
    model: &Model,
    shape: &Shape,
    blend: &Shape,
    stations: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<BlendContact>> {
    if stations < 2 {
        ogeom_bail!(
            Construction,
            "a blend measured at {stations} stations is not measured"
        );
    }
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    if !faces.iter().any(|f| f.node() == blend.node()) {
        ogeom_bail!(Construction, "that face is not part of this shape");
    }
    let surface_of =
        |face: &Shape| -> Option<(ogeom_topo::SurfaceId, ogeom_geom::SurfaceGeometry)> {
            let NodeData::Face(data) = model.node(face)?.data() else {
                return None;
            };
            let geometry = model.geometry().surface(data.surface)?.clone();
            Some((data.surface, geometry))
        };
    let Some((blend_id, blend_surface)) = surface_of(blend) else {
        ogeom_bail!(Construction, "the blend face carries no surface");
    };
    let blend_edges = explore_unique(model, blend, ShapeType::Edge)?;

    let mut out = Vec::new();
    for neighbour in &faces {
        if neighbour.node() == blend.node() {
            continue;
        }
        let Some((other_id, other_surface)) = surface_of(neighbour) else {
            continue;
        };
        for edge in explore_unique(model, neighbour, ShapeType::Edge)? {
            if !blend_edges.iter().any(|e| e.node() == edge.node()) {
                continue;
            }
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve).cloned() else {
                continue;
            };
            let chart =
                |id: ogeom_topo::SurfaceId| -> Option<(ogeom_geom::PlanarCurve, (f64, f64))> {
                    match data.pcurve_for(id, edge.location())? {
                        EdgeRepr::PCurve { curve, range, .. } => {
                            Some((model.geometry().pcurve(*curve)?.clone(), *range))
                        }
                        EdgeRepr::Seam { forward, range, .. } => {
                            Some((model.geometry().pcurve(*forward)?.clone(), *range))
                        }
                        _ => None,
                    }
                };
            let (Some((pc_blend, pr_blend)), Some((pc_other, pr_other))) =
                (chart(blend_id), chart(other_id))
            else {
                out.push(BlendContact {
                    edge: edge.clone(),
                    neighbour: neighbour.clone(),
                    tangency_error: f64::INFINITY,
                    curvature_error: f64::INFINITY,
                    gap: f64::INFINITY,
                    stations: 0,
                });
                continue;
            };
            let (mut worst_angle, mut worst_gap, mut worst_curvature) = (0.0f64, 0.0f64, 0.0f64);
            let mut curvature_read = false;
            for k in 0..stations {
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "a station index, far below the mantissa"
                )]
                let f = k as f64 / (stations - 1) as f64;
                let t = (range.1 - range.0).mul_add(f, range.0);
                let on_curve = geometry.point_at(t, tol)?;
                let tangent = geometry.d1_at(t, tol)?;
                let mut normals = Vec::with_capacity(2);
                let mut across = Vec::with_capacity(2);
                for ((pcurve, prange), surface) in [
                    ((&pc_blend, pr_blend), &blend_surface),
                    ((&pc_other, pr_other), &other_surface),
                ] {
                    let pt = (prange.1 - prange.0).mul_add(f, prange.0);
                    let uv = ogeom_geom::Curve2d::point_at(pcurve, pt, tol)?;
                    // A pcurve ending on a chart's pole may stand a rounding
                    // past the chart's window there; the station is read at
                    // the window's edge, and the gap says what that costs.
                    let uv = {
                        use ogeom_geom::Surface as _;
                        let ((u0, u1), (v0, v1)) = surface.domain();
                        ogeom_math::Point2::new(
                            if surface.is_periodic_u() {
                                uv.x
                            } else {
                                uv.x.clamp(u0, u1)
                            },
                            if surface.is_periodic_v() {
                                uv.y
                            } else {
                                uv.y.clamp(v0, v1)
                            },
                        )
                    };
                    worst_gap =
                        worst_gap.max(surface.point_at(uv.x, uv.y, tol)?.distance(on_curve));
                    // A station on a chart's pole (a corner patch's own
                    // corner sits on the ball's pole by construction) has
                    // no normal from the chart; the stations beside it say
                    // what the join does there.
                    if let Ok(normal) = surface.normal_at(uv.x, uv.y, tol) {
                        normals.push(normal.vector());
                        // The direction square to the edge in this
                        // surface's own tangent plane, and how the surface
                        // bends along it.
                        across.push(
                            surface
                                .curvature_at(uv.x, uv.y, tol)
                                .ok()
                                .and_then(|c| c.normal_curvature(normal.vector().cross(tangent))),
                        );
                    }
                }
                if normals.len() < 2 {
                    continue;
                }
                // Orientation is the topology's business, not the join's:
                // two faces meeting smoothly may still be wound opposite
                // ways, so the angle is taken to the nearer sense, and the
                // curvatures are compared against one normal.
                let along = normals[0].dot(normals[1]);
                worst_angle = worst_angle.max(along.abs().min(1.0).acos());
                if let (Some(a), Some(b)) = (across[0], across[1]) {
                    let b = if along < 0.0 { -b } else { b };
                    worst_curvature = worst_curvature.max((a - b).abs());
                    curvature_read = true;
                }
            }
            out.push(BlendContact {
                edge: edge.clone(),
                neighbour: neighbour.clone(),
                tangency_error: worst_angle,
                curvature_error: if curvature_read {
                    worst_curvature
                } else {
                    f64::INFINITY
                },
                gap: worst_gap,
                stations,
            });
        }
    }
    Ok(out)
}

/// The surface's curvature over a face, on a grid: what a curvature map or
/// a zebra display draws from.
///
/// The box the face's trim spans in its surface's parameters is cut into
/// `density` by `density` cells, and each cell's centre that lies inside
/// the trim gives one sample: the point in space and the principal
/// curvatures there. The trim is read as the mesher reads it, at the
/// default deflection. A centre where the surface has no normal (a pole,
/// an apex) has no curvature and gives no sample.
///
/// Each sample is placed where the face stands and signed against the
/// normal the face's mesh carries: the placement moves the point and the
/// directions, its scale divides the curvatures, and a reversed face turns
/// the normal round. A convex face reads negative, a cavity positive.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `face` is not a face or `density` is zero;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it is not
/// in this model.
pub fn face_curvature_samples(
    model: &Model,
    face: &Shape,
    density: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<(Point, SurfaceCurvature)>> {
    if density == 0 {
        ogeom_bail!(Construction, "a face sampled on no grid is not sampled");
    }
    let rings = ogeom_mesh::face_boundary(model, face, ogeom_mesh::Deflection::default(), tol)?;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;
    let flip = face.orientation() == Orientation::Reversed;

    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in rings.iter().flatten() {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }

    let mut out = Vec::new();
    for i in 0..density {
        for j in 0..density {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a grid index, far below the mantissa"
            )]
            let (fu, fv) = (
                (i as f64 + 0.5) / density as f64,
                (j as f64 + 0.5) / density as f64,
            );
            let uv = Point2::new(
                (hi.x - lo.x).mul_add(fu, lo.x),
                (hi.y - lo.y).mul_add(fv, lo.y),
            );
            if !uv.x.is_finite() || !uv.y.is_finite() || !ogeom_mesh::inside_boundary(&rings, uv) {
                continue;
            }
            let Ok(curvature) = surface.curvature_at(uv.x, uv.y, tol) else {
                continue;
            };
            let point = placement.apply(surface.point_at(uv.x, uv.y, tol)?);
            out.push((point, placed(&curvature, &placement, flip, tol)?));
        }
    }
    Ok(out)
}

/// The curvature of a face at one point of it, as [`face_curvature_at`]
/// reads it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceCurvature {
    /// The principal curvatures and directions, placed where the face
    /// stands and signed against its outward normal there.
    pub curvature: SurfaceCurvature,
    /// The point of the face they are read at, where the face stands.
    pub point: Point,
    /// That point's parameters on the face's surface: a hint for the next
    /// query nearby.
    pub uv: Point2,
    /// How far the queried point is from it.
    pub distance: f64,
}

/// The principal curvatures of `face`, as placed and as oriented, at the
/// point of its surface nearest `p`.
///
/// `p` is taken where the face stands. It is carried into the surface's own
/// frame and projected onto the surface: in closed form where the surface
/// has one, over its spans for a B-spline patch, from a seed grid for
/// anything else. A `hint` (surface parameters, such as the
/// [`FaceCurvature::uv`] of a query nearby) is refined from as well, and
/// its foot is kept wherever it is no farther than the global one by more
/// than the confusion distance, so a probe dragged across a surface with
/// several equally near feet (a point on a cylinder's axis) stays on the
/// one it had.
///
/// The foot must lie on the face: inside its trim, or on its boundary
/// within the face's tolerance and the trim's default deflection. One
/// outside it is refused rather than moved to the nearest point of the
/// boundary: what is asked about is the surface's nearest point.
///
/// The curvatures are signed as [`face_curvature_samples`] signs them:
/// against the face's outward normal, so a convex face reads negative and
/// a cavity positive. The placement's scale divides them, a reversed face
/// turns their sign round (the largest and smallest trading places), and a
/// mirror changes nothing: it carries the outward normal to the image's.
/// A plane reads `0` twice; a cylinder of radius `r` reads `0` along its
/// axis ([`SurfaceCurvature::max_direction`]) and `-1/r` round it.
///
/// Where every direction bends alike (an umbilic: every point of a sphere
/// or a plane, and the pole of a surface of revolution) the two principal
/// curvatures are equal and the directions are any perpendicular pair of
/// tangents; [`SurfaceCurvature::is_umbilic`] says so.
///
/// Where the surface's parameterization degenerates (the pole of a sphere,
/// a B-spline patch with a collapsed row of control points) its own
/// derivatives give no normal. The curvature there is the limit of the
/// curvatures read approaching it along the parameter that does not
/// collapse, which is the curvature of the geometry whenever the geometry
/// is smooth there. Where it has no limit (a cone's apex, where one
/// principal curvature grows without bound) the point is refused.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `face` is not a face;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it is not
/// in this model;
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `p` is not
/// finite, if the surface's nearest point lies outside the face, or if the
/// surface has no curvature there (an apex).
pub fn face_curvature_at(
    model: &Model,
    face: &Shape,
    p: Point,
    hint: Option<Point2>,
    tol: Tolerances,
) -> OgeomResult<FaceCurvature> {
    if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
        ogeom_bail!(Domain, "curvature asked at a point that is not finite");
    }
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;
    let local = placement.inverse()?.apply(p);

    // The grid only seeds a surface with neither a closed form nor spans;
    // 32 a side tells the near and far sides of a folded surface apart, as
    // the face classifier finds.
    let mut foot = ogeom_algo::project_on_surface(surface, local, 32, tol)?;
    if let Some(hint) = hint
        && let Ok(near) = ogeom_algo::project_on_surface_from(surface, local, (hint.x, hint.y), tol)
        && near.distance <= foot.distance + tol.confusion()
    {
        foot = near;
    }

    let trim =
        ogeom_algo::classify::FaceTrim::of(model, face, ogeom_mesh::Deflection::default(), tol)?;
    let reach = tol.confusion().max(data.tolerance.get());
    if trim.holds(model, foot.parameters, reach, tol)? == ogeom_algo::Containment::Out {
        ogeom_bail!(
            Domain,
            "the nearest point of the surface, at ({}, {}), lies outside the face",
            foot.parameters.0,
            foot.parameters.1
        );
    }

    let (u, v) = foot.parameters;
    let curvature = surface_curvature_or_limit(surface, u, v, tol)?;
    let flip = face.orientation() == Orientation::Reversed;
    Ok(FaceCurvature {
        curvature: placed(&curvature, &placement, flip, tol)?,
        point: placement.apply(foot.point),
        uv: Point2::new(u, v),
        distance: foot.distance * placement.scale_factor(),
    })
}

/// The surface's curvature at `(u, v)`, or where its parameterization
/// degenerates there, the limit of it approaching along the parameter
/// that does not collapse.
///
/// The surface is read twice, a ten-thousandth and a hundred-thousandth of
/// that parameter's range in from the point. Where the two readings agree
/// to a hundredth (a smooth surface's curvature changes by the step's
/// order between them, an apex's by the steps' ratio) the curvatures and
/// the normal are carried on to the point along the line through them,
/// and the directions are laid in the tangent plane of that normal.
fn surface_curvature_or_limit(
    surface: &ogeom_geom::SurfaceGeometry,
    u: f64,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<SurfaceCurvature> {
    if !surface.is_degenerate_at(u, v, tol)?
        && let Ok(c) = surface.curvature_at(u, v, tol)
    {
        return Ok(c);
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let (du, dv) = surface.d1_at(u, v, tol)?;
    // The parameter whose derivative survives is the one to approach along,
    // inward from whichever end of its range is nearer.
    let along_v = du.magnitude() <= dv.magnitude();
    let (at, (lo, hi)) = if along_v {
        (v, (va, vb))
    } else {
        (u, (ua, ub))
    };
    let span = if (hi - lo).is_finite() { hi - lo } else { 1.0 };
    let inward = if at - lo <= hi - at { 1.0 } else { -1.0 };
    let read = |fraction: f64| -> Option<SurfaceCurvature> {
        let t = (inward * fraction).mul_add(span, at).clamp(lo, hi);
        let (s, w) = if along_v { (u, t) } else { (t, v) };
        surface.curvature_at(s, w, tol).ok()
    };
    let (Some(far), Some(near)) = (read(1e-4), read(1e-5)) else {
        ogeom_bail!(
            Domain,
            "the surface has no normal near ({u}, {v}), so no curvature there"
        );
    };
    let agree = |a: f64, b: f64| (a - b).abs() <= 1e-2 * a.abs().max(b.abs()) + tol.angular();
    if !(agree(far.max, near.max) && agree(far.min, near.min)) {
        ogeom_bail!(
            Domain,
            "the surface's curvature grows without bound towards ({u}, {v}), an apex"
        );
    }
    // Each reading is the limit plus a term of the step's order: on the
    // line through them, the step's ninth past the nearer.
    let on = |far: f64, near: f64| (near - far).mul_add(1.0 / 9.0, near);
    let (max, min) = (on(far.max, near.max), on(far.min, near.min));
    let towards = |far: Direction, near: Direction| -> Vector {
        (near.vector() - far.vector()) * (1.0 / 9.0) + near.vector()
    };
    let normal = Direction::new(towards(far.normal, near.normal), tol)?;
    let n = normal.vector();
    let d = near.max_direction.vector();
    let max_direction = Direction::new(d - n * d.dot(n), tol)?;
    let min_direction = Direction::new(n.cross(max_direction.vector()), tol)?;
    Ok(if max >= min {
        SurfaceCurvature {
            max,
            min,
            max_direction,
            min_direction,
            normal,
        }
    } else {
        SurfaceCurvature {
            max: min,
            min: max,
            max_direction: min_direction,
            min_direction: max_direction,
            normal,
        }
    })
}

/// A curvature read in a surface's own frame, carried to where its face
/// stands and signed against the face's normal there.
///
/// A similarity `s R` maps the normal `n` to the direction of `s R n`, and
/// every normal curvature `k` to `k / s` against it. A mirror among them
/// changes neither: it carries the outward normal to the image's outward
/// normal, as the mesher carries it. Turning the normal round for a
/// reversed face negates the curvatures, and the largest and smallest
/// trade places.
fn placed(
    c: &SurfaceCurvature,
    placement: &ogeom_math::Transform,
    flip: bool,
    tol: Tolerances,
) -> OgeomResult<SurfaceCurvature> {
    let carry = |d: Direction| Direction::new(placement.apply_vector(d.vector()), tol);
    let sign = if flip { -1.0 } else { 1.0 };
    let scale = placement.scale_factor();
    let (a, b) = (sign * c.max / scale, sign * c.min / scale);
    let normal = carry(c.normal)?;
    let normal = if flip { -normal } else { normal };
    let (max_direction, min_direction) = (carry(c.max_direction)?, carry(c.min_direction)?);
    Ok(if a >= b {
        SurfaceCurvature {
            max: a,
            min: b,
            max_direction,
            min_direction,
            normal,
        }
    } else {
        SurfaceCurvature {
            max: b,
            min: a,
            max_direction: min_direction,
            min_direction: max_direction,
            normal,
        }
    })
}
