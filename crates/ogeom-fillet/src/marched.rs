//! The marched fillet: a rolling-ball blend whose seat has no closed form,
//! carried all the way to topology.
//!
//! [`march_blend`](crate::march_blend) solves the ball's two contact points
//! station by station; this module turns those stations into the same wedge
//! the closed-form blends build. The blend face is a surface fitted through
//! the ball's own arcs, its rails are the surface's own border iso-curves,
//! and the legs are patches of the hosts' *own* surfaces — exact geometry
//! bounded by fitted rails — so the boolean's same-domain resolution melts
//! them exactly as it melts a closed-form wedge's.
//!
//! The chart discipline that makes it sound: the marcher *solves* the
//! contact parameters rather than projecting, so a pcurve fitted through
//! them at the grid's own parameters is same-parameter by construction; the
//! blend chart's rail images are straight iso rows; and the apex ring's
//! images come from closed-form inversion on the host charts, unwrapped for
//! continuity.

use crate::march::{BlendStop, Sides, march_blend_seeded};
use crate::support::{apply_wedge, edge_curve};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry};
use ogeom_intersect::Marching;
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// How many samples cross each blend arc.
const ACROSS: usize = 9;

/// Round an edge whose seat only the marching machinery can speak: the
/// intersection curve of two analytic faces, closed on itself.
pub(crate) fn marched_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (stored_guide, edge_range) = edge_curve(model, edge, tol)?;
    // The seat is the whole loop even when the boolean split it into arcs:
    // the apex ring runs the curve's full turn, and every arc of the old
    // seat melts with the legs. A conic arc is re-opened to its full
    // period; a fitted seam already spans its loop.
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    let (guide, guide_range) = if closed {
        (stored_guide, edge_range)
    } else {
        match &stored_guide {
            Curve::Ellipse(e) => {
                let full: Curve = ogeom_geom::EllipseCurve::new(e.ellipse()).into();
                let domain = full.domain();
                (full, domain)
            }
            _ => {
                let domain = stored_guide.domain();
                (stored_guide, domain)
            }
        }
    };

    // The two host faces at the edge, with their surfaces and outward signs.
    let mut hosts: Vec<(Shape, SurfaceGeometry, f64)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| crate::support::same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        // Baked into the world: the face's surface lives wherever its
        // placement puts it, and everything below — the march, the legs,
        // the melt — speaks world coordinates.
        let surface = {
            use ogeom_geom::Transformable as _;
            let placement = face.transform(model.datums())?;
            stored.transformed(&placement, tol)?
        };
        match surface {
            SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) => {}
            SurfaceGeometry::Cone(_) | SurfaceGeometry::Sphere(_) | SurfaceGeometry::Torus(_) => {
                ogeom_bail!(
                    Construction,
                    "a marched fillet on a cone, sphere or torus host needs \
                     that chart's inversion carried; planes and cylinders \
                     are what this speaks — docs/PARITY.md, fillet.edge-blends"
                )
            }
            _ => ogeom_bail!(
                Construction,
                "a marched fillet's hosts must be analytic; a fitted host \
                 has no chart the legs can melt against — docs/PARITY.md, \
                 fillet.edge-blends"
            ),
        }
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        hosts.push((face, surface, sign));
    }
    let [(face_first, first, sign_first), (_, second, sign_second)] = hosts.as_slice() else {
        ogeom_bail!(
            Construction,
            "a marched fillet needs an edge shared by exactly two faces, \
             found {}",
            hosts.len()
        );
    };
    let face_first = face_first.clone();
    let (first, second) = (first.clone(), second.clone());
    let (sign_first, sign_second) = (*sign_first, *sign_second);

    // A seat running through a point where its hosts are tangent — the
    // crossing of two equal drums — has no section there: the ball's arc
    // collapses at the pole, and the march can only stall on it. Refused by
    // name up front, sampled along the whole reconstructed loop.
    {
        use ogeom_geom::Surface as _;
        for i in 0..64 {
            #[allow(clippy::cast_precision_loss)]
            let t = guide_range.0 + (guide_range.1 - guide_range.0) * (i as f64) / 64.0;
            let p = guide.point_at(t, tol)?;
            let normal_of = |surface: &SurfaceGeometry| -> OgeomResult<Vector> {
                let projection = ogeom_algo::project_on_surface(surface, p, 32, tol)?;
                let (u, v) = projection.parameters;
                let (du, dv) = surface.d1_at(u, v, tol)?;
                let n = du.cross(dv);
                Ok(n / n.magnitude())
            };
            let (n1, n2) = (normal_of(&first)?, normal_of(&second)?);
            if n1.cross(n2).magnitude() <= 1e-2 {
                ogeom_bail!(
                    Construction,
                    "the seat passes through a point where its two hosts are \
                     tangent; the ball's section collapses at that pole and \
                     the pinched seam is refused — docs/PARITY.md, \
                     fillet.edge-blends"
                );
            }
        }
    }

    // Convexity, read from the solid itself the way the planar seat reads
    // it: which way the first face extends from the edge, leaned against
    // the second's outward normal. It decides everything downstream — which
    // of the four ball seatings is the fillet's, and whether the wedge adds
    // or removes.
    //
    // Probed at the *edge's* own midpoint, not the reconstructed loop's: a
    // conic arc re-opened to its full period runs through territory the
    // boolean cut away, and a probe standing off the solid reads nothing in
    // either direction. The edge's midpoint is on the crease by definition.
    use ogeom_geom::Surface as _;
    let mid_t = f64::midpoint(edge_range.0, edge_range.1);
    let mid = guide.point_at(mid_t, tol)?;
    let outward_at = |surface: &SurfaceGeometry, sign: f64| -> OgeomResult<Vector> {
        let projection = ogeom_algo::project_on_surface(surface, mid, 32, tol)?;
        let (u, v) = projection.parameters;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        Ok(n / n.magnitude() * sign)
    };
    let n1 = outward_at(&first, sign_first)?;
    let n2 = outward_at(&second, sign_second)?;
    let convex = {
        let tangent = {
            let d = guide.d1_at(mid_t, tol)?;
            d / d.magnitude()
        };
        let raw = {
            let t = n1.cross(tangent);
            let m = t.magnitude();
            if m <= tol.angular() {
                ogeom_bail!(Construction, "a face is tangent to its own edge");
            }
            t / m
        };
        let span = guide.point_at(edge_range.0, tol)?.distance(mid).max(radius);
        let mut extends: Option<Vector> = None;
        'scales: for scale in [1e-3, 1e-2, 5e-2] {
            let eps = span * scale;
            let deflection = ogeom_mesh::Deflection {
                chord: eps * 0.1,
                ..ogeom_mesh::Deflection::default()
            };
            for dir in [raw, -raw] {
                // The step is chordal; on a curved host it leaves the
                // surface quadratically, and an off-surface probe classifies
                // as nothing. Project it home first.
                let probe = ogeom_algo::project_on_surface(&first, mid + dir * eps, 32, tol)?.point;
                if ogeom_algo::classify_on_face(model, &face_first, probe, deflection, tol)?
                    == ogeom_algo::Containment::In
                {
                    extends = Some(dir);
                    break 'scales;
                }
            }
        }
        let Some(extends) = extends else {
            ogeom_bail!(
                Construction,
                "cannot read which way the edge's face extends; the face is \
                 thinner than the probe can resolve"
            );
        };
        let lean = extends.dot(n2);
        if lean.abs() <= tol.angular() {
            ogeom_bail!(
                Construction,
                "the edge's faces are tangent; there is no corner"
            );
        }
        lean < 0.0
    };
    // The fillet's ball rides the material's own side of each support: its
    // centre sits inside the material at a convex corner and out in the
    // notch at a concave one.
    let seat_sign = if convex { -1.0 } else { 1.0 };
    #[allow(clippy::cast_possible_truncation)]
    let sides = Sides {
        first: (seat_sign * sign_first) as i8,
        second: (seat_sign * sign_second) as i8,
    };

    // Seeded at the edge's own midpoint: on a reconstructed loop the domain
    // midpoint may stand in cut-away territory where no ball seats, but the
    // crease's own midpoint is seat by definition, and a closed loop closes
    // from wherever the walker starts.
    let blend = march_blend_seeded(
        &first,
        &second,
        radius,
        &guide,
        sides,
        mid_t,
        Marching {
            chord: 3e-6,
            ..Marching::default()
        },
        tol,
    )?;
    if blend.stopped != BlendStop::Closed {
        ogeom_bail!(
            Construction,
            "the blend did not close on itself ({:?}); the open seat's \
             run-out is still owed — docs/PARITY.md, fillet.edge-blends",
            blend.stopped
        );
    }
    if blend.len() < 8 {
        ogeom_bail!(Construction, "the march produced too few stations to fit");
    }
    // The band wants every winding rail to run its period forward; when the
    // march went the other way round, the whole loop reverses.
    let mut blend = blend;
    // A closed march may hand its first station back as its last; the loop
    // owns it once, and the grid closes itself.
    while blend.len() > 8
        && blend.spine[blend.len() - 1].distance(blend.spine[0]) <= tol.confusion() * 100.0
    {
        blend.spine.pop();
        blend.touch_first.pop();
        blend.touch_second.pop();
        blend.on_first.pop();
        blend.on_second.pop();
        blend.along.pop();
    }
    // A march may also overshoot its start before noticing closure, leaving
    // trailing stations that re-trace the loop's opening arc: the sequence
    // folds back on itself and no smooth fit can follow it. The overshoot
    // reads off the closing step — while it points *against* the march, the
    // last station is past the start and goes.
    while blend.len() > 8 {
        let last = blend.spine[blend.len() - 1];
        let prev = blend.spine[blend.len() - 2];
        if (blend.spine[0] - last).dot(last - prev) >= 0.0 {
            break;
        }
        blend.spine.pop();
        blend.touch_first.pop();
        blend.touch_second.pop();
        blend.on_first.pop();
        blend.on_second.pop();
        blend.along.pop();
    }
    {
        let winding = |on: &[(f64, f64)], surface: &SurfaceGeometry| -> f64 {
            let Some(period) = period_of(surface) else {
                return 0.0;
            };
            // Unwrapped: the walker clamps parameters into the window, and a
            // wrapped sequence reads as no winding at all.
            let mut last = on[0].0;
            let mut total = 0.0;
            for &(u, _) in &on[1..] {
                let mut step = u - last;
                while step > period / 2.0 {
                    step -= period;
                }
                while step < -period / 2.0 {
                    step += period;
                }
                total += step;
                last = u;
            }
            // The loop's closing step back to the start.
            let mut close = on[0].0 - last;
            while close > period / 2.0 {
                close -= period;
            }
            while close < -period / 2.0 {
                close += period;
            }
            ((total + close) / period).round()
        };
        let (w1, w2) = (
            winding(&blend.on_first, &first),
            winding(&blend.on_second, &second),
        );
        if w1 * w2 < 0.0 {
            ogeom_bail!(
                Construction,
                "the seat winds its two hosts in opposite senses; that \
                 configuration is still owed — docs/PARITY.md, \
                 fillet.edge-blends"
            );
        }
        if w1 < 0.0 || w2 < 0.0 {
            blend.spine.reverse();
            blend.touch_first.reverse();
            blend.touch_second.reverse();
            blend.on_first.reverse();
            blend.on_second.reverse();
            blend.along.reverse();
        }

        // A winding leg closes as a band, and the band's connector joins the
        // ring *starts*. Anchor the march so station zero stands on the
        // guide's own start column: the connector then runs (nearly) up the
        // iso between the rings, where it cannot cross either of them — a
        // connector thrown diagonally across the band does, and the
        // arrangement downstream cannot hold strands that cross mid-span.
        let anchor_host = if w1.abs() > 0.5 {
            Some((&first, &blend.on_first))
        } else if w2.abs() > 0.5 {
            Some((&second, &blend.on_second))
        } else {
            None
        };
        if let Some((host, on)) = anchor_host
            && let Some(period) = period_of(host)
        {
            let apex_u = chart_of(host, guide.point_at(guide_range.0, tol)?, None, tol)?.x;
            let circular = |u: f64| -> f64 {
                let d = (u - apex_u).rem_euclid(period);
                d.min(period - d)
            };
            let mut k = 0;
            for (i, &(u, _)) in on.iter().enumerate() {
                if circular(u) < circular(on[k].0) {
                    k = i;
                }
            }
            blend.spine.rotate_left(k);
            blend.touch_first.rotate_left(k);
            blend.touch_second.rotate_left(k);
            blend.on_first.rotate_left(k);
            blend.on_second.rotate_left(k);
            blend.along.rotate_left(k);
            // The march's closing step may be far shorter than its stride;
            // rotated into the loop's interior, that cramped pair would put
            // two grid columns nearly on top of each other and poison the
            // fits' parameterization. One sweep drops any station standing
            // within a fraction of the loop's median stride of its
            // predecessor.
            let mut steps: Vec<f64> = blend
                .spine
                .windows(2)
                .map(|w| w[0].distance(w[1]))
                .collect();
            steps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            let cramped = steps.get(steps.len() / 2).copied().unwrap_or(0.0) * 0.25;
            let mut i = 1;
            while i < blend.len() {
                if blend.spine[i].distance(blend.spine[i - 1]) <= cramped && blend.len() > 8 {
                    blend.spine.remove(i);
                    blend.touch_first.remove(i);
                    blend.touch_second.remove(i);
                    blend.on_first.remove(i);
                    blend.on_second.remove(i);
                    blend.along.remove(i);
                } else {
                    i += 1;
                }
            }
        }
    }

    // A convex corner's wedge is material removed; a concave notch's is
    // material added.
    let additive = !convex;

    let n = blend.len();
    let fit_target = (tol.confusion() * 1e3).max(1e-4);

    // The blend surface: each station's exact ball arc, the loop of stations
    // fitted *closed* — the join is C1 wherever the seam lands, so anchoring
    // the seam column costs nothing.
    let mut rows: Vec<Vec<Point>> = (0..ACROSS).map(|_| Vec::with_capacity(n + 1)).collect();
    for i in 0..=n {
        let at = i % n;
        let centre = blend.spine[at];
        let a = (blend.touch_first[at] - centre) / radius;
        let b = (blend.touch_second[at] - centre) / radius;
        let cross = a.cross(b);
        let m = cross.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(
                Construction,
                "a blend section collapsed; the radius wedges rather than \
                 seats at station {at}"
            );
        }
        let axis = cross / m;
        let sweep = a.dot(b).clamp(-1.0, 1.0).acos();
        for (k, row) in rows.iter_mut().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let theta = sweep * (k as f64) / ((ACROSS - 1) as f64);
            let dir = a * theta.cos() + axis.cross(a) * theta.sin();
            row.push(centre + dir * radius);
        }
    }
    // The closed direction is `v` by the fit's convention, so the grid goes
    // in station-major: one row per station arc, first repeated at the end.
    let loops: Vec<Vec<Point>> = (0..=n)
        .map(|i| (0..ACROSS).map(|k| rows[k][i]).collect())
        .collect();
    let fitted = ogeom_geom::fit::fit_surface_grid_closed_v_chordal(&loops, 3, fit_target, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the blend surface reached {} against a target of {fit_target}",
            fitted.error
        );
    }

    // The loop's own parameters, recomputed the way the fit computes them
    // (averaged centripetal), so the rail pcurves fitted at them are
    // same-parameter with the surface's borders.
    let u_params = averaged_chordal(&rows);

    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k_count, l_count, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l_count + j] };
    let blend_geo: SurfaceGeometry = surface.into();
    let (u_dom, v_dom) = blend_geo.domain();
    let blend_id = model.geometry_mut().add_surface(blend_geo.clone());

    // Rails and seam straight off the control net, exactly as a skin's:
    // the rails are the station loops at the arc's two ends, the seam the
    // station-zero column across the arc.
    let border = |i: usize| -> OgeomResult<Curve> {
        let control: Vec<Point> = (0..l_count).map(|j| point_at(i, j)).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?))
    };
    let seam_curve = {
        let control: Vec<Point> = (0..k_count).map(|i| point_at(i, 0)).collect();
        Curve::BSpline(ogeom_geom::BSplineCurve::new(u_knots, control, tol)?)
    };
    let rail_first = ogeom_algo::make_edge(model, border(0)?, v_dom, tol)?.shape;
    let rail_second = ogeom_algo::make_edge(model, border(k_count - 1)?, v_dom, tol)?.shape;
    // The rails carry the fit's honest slop: every downstream filter that
    // compares them against exact geometry — the melt's crossing paver above
    // all — widens by an edge's recorded tolerance, not by wishful thinking.
    for rail in [&rail_first, &rail_second] {
        if let Some(node) = model.node_mut(rail)
            && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
        {
            data.tolerance = data.tolerance.widen_to(fit_target);
        }
    }
    let anchor0 = ogeom_algo::edge_vertices(model, &rail_first)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a rail has no vertex"))?;
    let anchor1 = ogeom_algo::edge_vertices(model, &rail_second)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a rail has no vertex"))?;
    let seam =
        ogeom_algo::make_edge_between(model, seam_curve, u_dom, &anchor0, &anchor1, tol)?.shape;

    let row_line = |v: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_pcurve(
        model,
        &rail_first,
        column_line(u_dom.0)?,
        blend_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &rail_second,
        column_line(u_dom.1)?,
        blend_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        row_line(v_dom.0)?,
        row_line(v_dom.1)?,
        blend_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    let blend_face = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                seam.clone(),
                rail_second.clone(),
                seam.reversed(),
                rail_first.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, blend_id, std::slice::from_ref(&wire), tol)?.shape;
        // The wedge's outward at the blend: towards the ball's centre when
        // cutting — the wedge is the corner material the ball displaced —
        // and away from it when fusing.
        let mid_u = f64::midpoint(u_dom.0, u_dom.1);
        let mid_v = f64::midpoint(v_dom.0, v_dom.1);
        let p = blend_geo.point_at(mid_u, mid_v, tol)?;
        let (du, dv) = blend_geo.d1_at(mid_u, mid_v, tol)?;
        let towards_centre = (blend.spine[n / 2] - p).dot(du.cross(dv)) > 0.0;
        if towards_centre == !additive {
            face
        } else {
            face.reversed()
        }
    };

    // The legs: one per host, the exact host surface bounded by the fitted
    // rail and a fresh ring on the edge's own curve.
    let leg_first = host_leg(
        model,
        &first,
        &blend.on_first,
        &u_params,
        &rail_first,
        &guide,
        guide_range,
        fit_target,
        tol,
    )?;
    let leg_second = host_leg(
        model,
        &second,
        &blend.on_second,
        &u_params,
        &rail_second,
        &guide,
        guide_range,
        fit_target,
        tol,
    )?;
    // Legs coincide with the solid's own faces: aligned when subtracting,
    // opposed when fusing, which is what the melt needs either way.
    let orient = |face: Shape, host_sign: f64| -> Shape {
        let aligned = if additive { -host_sign } else { host_sign };
        if aligned > 0.0 { face } else { face.reversed() }
    };
    let faces = [
        orient(leg_first, sign_first),
        orient(leg_second, sign_second),
        blend_face,
    ];
    apply_wedge(model, solid, Some(edge), &faces, additive, tol)
}

/// One leg: the host's own surface bounded by the marched rail and the
/// edge's ring. A rail that winds the chart's period takes the band with its
/// connector seam; a contractible one takes the annular two-wire face.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn host_leg(
    model: &mut Model,
    host: &SurfaceGeometry,
    on_host: &[(f64, f64)],
    u_params: &[f64],
    rail: &Shape,
    guide: &Curve,
    guide_range: (f64, f64),
    fit_target: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let n = on_host.len();
    // The rail's chart image, closed: the marcher solved these parameters,
    // so the fit through them at the grid's own u-parameters is
    // same-parameter with the rail's curve.
    let mut rail_chart: Vec<Point2> = on_host.iter().map(|(u, v)| Point2::new(*u, *v)).collect();
    // The walker clamps parameters into the chart window, so a loop that
    // crosses the seam comes back wrapped; unwrap it before anything is
    // fitted through it.
    if let Some(period) = period_of(host) {
        for i in 1..rail_chart.len() {
            let mut u = rail_chart[i].x;
            while u - rail_chart[i - 1].x > period / 2.0 {
                u -= period;
            }
            while rail_chart[i - 1].x - u > period / 2.0 {
                u += period;
            }
            rail_chart[i].x = u;
        }
    }
    let winding = period_of(host).map_or(0.0, |period| {
        // The unwrapped loop's chart displacement over one closing step,
        // rounded to whole periods.
        let delta = rail_chart[n - 1].x - rail_chart[0].x + (rail_chart[1].x - rail_chart[0].x);
        (delta / period).round() * period
    });
    rail_chart.push(Point2::new(rail_chart[0].x + winding, rail_chart[0].y));
    let rail_pcurve = {
        let fitted =
            ogeom_geom::fit::fit_points_2d_at_closed(u_params, &rail_chart, 3, fit_target, tol)?;
        if !fitted.met {
            ogeom_bail!(
                NotDone,
                "a rail's chart image reached {} against a target of {fit_target}",
                fitted.error
            );
        }
        PlanarCurve::from(fitted.curve)
    };

    // The apex ring: a fresh closed edge on the guide's own curve — turned
    // round when its own parameterization winds the chart backward, since
    // the band runs every period forward — its chart image inverted in
    // closed form and unwrapped for continuity.
    let chart_run = |guide: &Curve| -> OgeomResult<(Vec<f64>, Vec<Point2>)> {
        let samples = 96;
        let mut params: Vec<f64> = Vec::with_capacity(samples + 1);
        let mut chart: Vec<Point2> = Vec::with_capacity(samples + 1);
        let mut prev: Option<Point2> = None;
        for i in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let t = guide_range.0 + (guide_range.1 - guide_range.0) * (i as f64) / (samples as f64);
            let p = guide.point_at(t, tol)?;
            let uv = chart_of(host, p, prev, tol)?;
            params.push(t);
            chart.push(uv);
            prev = Some(uv);
        }
        Ok((params, chart))
    };
    let (mut apex_params, mut apex_chart) = chart_run(guide)?;
    let mut apex_guide = guide.clone();
    if apex_chart[apex_chart.len() - 1].x < apex_chart[0].x - 1e-6 {
        use ogeom_geom::Reversible as _;
        apex_guide = guide.clone().reversed();
        (apex_params, apex_chart) = chart_run(&apex_guide)?;
    }
    let apex = ogeom_algo::make_edge(model, apex_guide.clone(), guide_range, tol)?.shape;
    // Exact wherever the chart has a closed form — coinciding to machine
    // precision with the strands the boolean draws from the same curve is
    // what keeps slivers out of the arrangement — fitted only past that.
    // The leg's surface, windowed to its own neighbourhood. A leg carrying
    // the host's whole domain pairs with faces far from the seat, and the
    // duplicate strands it draws there interleave with the faces' own
    // boundaries as slivers no classifier can hold.
    let host = &windowed(host, &apex_chart, &rail_chart, tol)?;
    let apex_pcurve = match ogeom_intersect::exact_pcurve_of(&apex_guide, host, tol) {
        Some(exact) => exact,
        None => {
            let fitted = ogeom_geom::fit::fit_points_2d_at_closed(
                &apex_params,
                &apex_chart,
                3,
                fit_target,
                tol,
            )?;
            if !fitted.met {
                ogeom_bail!(
                    NotDone,
                    "the edge's chart image reached {} against a target of {fit_target}",
                    fitted.error
                );
            }
            PlanarCurve::from(fitted.curve)
        }
    };

    if winding.abs() > 1e-6 {
        // Both loops wind the period; the band with its connector closes the
        // strip between them. The band wants the period run *forward*.
        if winding < 0.0 {
            ogeom_bail!(
                Construction,
                "the seat winds against its host's chart; reversing the \
                 guide is still owed — docs/PARITY.md, fillet.edge-blends"
            );
        }
        ogeom_algo::make_band_between(
            model,
            host,
            [(&apex, apex_pcurve), (rail, rail_pcurve.clone())],
            tol,
        )
    } else {
        // Contractible loops: an annular patch, outer wire first.
        let apex_area = chart_area(&apex_chart);
        let rail_area = chart_area(&rail_chart);
        let surface_id = model.geometry_mut().add_surface(host.clone());
        ogeom_algo::attach_pcurve(
            model,
            &apex,
            apex_pcurve,
            surface_id,
            ogeom_topo::Location::identity(),
            guide_range,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            rail,
            rail_pcurve,
            surface_id,
            ogeom_topo::Location::identity(),
            (u_params[0], u_params[u_params.len() - 1]),
        )?;
        let apex_wire = ogeom_algo::make_wire(model, std::slice::from_ref(&apex), tol)?.shape;
        let rail_wire = ogeom_algo::make_wire(model, std::slice::from_ref(rail), tol)?.shape;
        let wires = if apex_area.abs() >= rail_area.abs() {
            [apex_wire, rail_wire]
        } else {
            [rail_wire, apex_wire]
        };
        Ok(ogeom_algo::make_face_on(model, surface_id, &wires, tol)?.shape)
    }
}

/// The host surface cut down to the window the leg actually spans, padded a
/// little, so the leg's bound stays near the seat.
fn windowed(
    host: &SurfaceGeometry,
    apex_chart: &[Point2],
    rail_chart: &[Point2],
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let mut u = (f64::INFINITY, f64::NEG_INFINITY);
    let mut v = (f64::INFINITY, f64::NEG_INFINITY);
    for p in apex_chart.iter().chain(rail_chart) {
        u = (u.0.min(p.x), u.1.max(p.x));
        v = (v.0.min(p.y), v.1.max(p.y));
    }
    let pad = |(lo, hi): (f64, f64)| {
        let span = (hi - lo).max(tol.confusion() * 1e3);
        (span.mul_add(-0.25, lo), span.mul_add(0.25, hi))
    };
    Ok(match host {
        SurfaceGeometry::Cylinder(c) => {
            let (lo, hi) = pad(v);
            ogeom_geom::CylinderSurface::new(c.cylinder(), (lo, hi))?.into()
        }
        SurfaceGeometry::Plane(p) => {
            let (ulo, uhi) = pad(u);
            let (vlo, vhi) = pad(v);
            ogeom_geom::PlaneSurface::over(p.plane(), (ulo, uhi), (vlo, vhi))?.into()
        }
        other => other.clone(),
    })
}

/// The chart period in `u`, for surfaces that have one.
fn period_of(surface: &SurfaceGeometry) -> Option<f64> {
    match surface {
        SurfaceGeometry::Cylinder(_) => Some(core::f64::consts::TAU),
        _ => None,
    }
}

/// Closed-form chart inversion, unwrapped against the previous sample.
fn chart_of(
    surface: &SurfaceGeometry,
    p: Point,
    prev: Option<Point2>,
    tol: Tolerances,
) -> OgeomResult<Point2> {
    let _ = tol;
    let raw = match surface {
        SurfaceGeometry::Plane(pl) => {
            let local = pl.plane().frame().to_local(p);
            Point2::new(local.x, local.y)
        }
        SurfaceGeometry::Cylinder(c) => {
            let local = c.cylinder().frame().to_local(p);
            Point2::new(local.y.atan2(local.x), local.z)
        }
        _ => ogeom_bail!(
            Construction,
            "no closed-form chart inversion for this surface"
        ),
    };
    let Some(prev) = prev else {
        return Ok(raw);
    };
    let Some(period) = period_of(surface) else {
        return Ok(raw);
    };
    let mut u = raw.x;
    while u - prev.x > period / 2.0 {
        u -= period;
    }
    while prev.x - u > period / 2.0 {
        u += period;
    }
    Ok(Point2::new(u, raw.y))
}

/// The signed area a closed chart image encloses, by the shoelace.
fn chart_area(points: &[Point2]) -> f64 {
    let mut sum = 0.0;
    for pair in points.windows(2) {
        sum += pair[0].x.mul_add(pair[1].y, -(pair[1].x * pair[0].y));
    }
    sum / 2.0
}

/// The fit's own u-parameters: averaged centripetal across the rows, the
/// last pinned to one, exactly as the surface fit assigns them.
fn averaged_chordal(rows: &[Vec<Point>]) -> Vec<f64> {
    let len = rows[0].len();
    let mut sums = vec![0.0_f64; len];
    for row in rows {
        let mut partial = Vec::with_capacity(len);
        partial.push(0.0);
        let mut total = 0.0;
        for pair in row.windows(2) {
            total += pair[0].distance(pair[1]);
            partial.push(total);
        }
        if total > 0.0 {
            for p in &mut partial {
                *p /= total;
            }
        }
        if let Some(last) = partial.last_mut() {
            *last = 1.0;
        }
        for (s, p) in sums.iter_mut().zip(partial) {
            *s += p;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let count = rows.len() as f64;
    sums.iter().map(|s| s / count).collect()
}
