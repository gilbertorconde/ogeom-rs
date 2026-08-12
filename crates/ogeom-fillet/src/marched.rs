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

use crate::march::{BlendStop, march_blend};
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
    let (guide, edge_range) = edge_curve(model, edge)?;
    // The seat is the whole loop even when the boolean split it into arcs:
    // the apex ring runs the curve's full turn, and every arc of the old
    // seat melts with the legs.
    let guide_range = {
        let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
        if closed { edge_range } else { guide.domain() }
    };

    // The two host faces at the edge, with their surfaces and outward signs.
    let mut hosts: Vec<(Shape, SurfaceGeometry, f64)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| e.node() == edge.node());
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
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
    let [(_, first, sign_first), (_, second, sign_second)] = hosts.as_slice() else {
        ogeom_bail!(
            Construction,
            "a marched fillet needs an edge shared by exactly two faces, \
             found {}",
            hosts.len()
        );
    };
    let (first, second) = (first.clone(), second.clone());
    let (sign_first, sign_second) = (*sign_first, *sign_second);

    let blend = march_blend(
        &first,
        &second,
        radius,
        &guide,
        Marching {
            chord: 1e-5,
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
    }
    {
        let winding = |on: &[(f64, f64)], surface: &SurfaceGeometry| -> f64 {
            let Some(period) = period_of(surface) else {
                return 0.0;
            };
            let n = on.len();
            let delta = on[n - 1].0 - on[0].0 + (on[1].0 - on[0].0);
            (delta / period).round()
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
        }
    }

    // Whether the wedge adds or removes: where the ball's centre sits. A
    // centre inside the material is a ball rolling a concave seat, and the
    // wedge fuses.
    let additive =
        ogeom_algo::classify_in_solid_exact(model, solid, blend.spine[blend.len() / 2], tol)?
            == ogeom_algo::Containment::In;

    let n = blend.len();
    let fit_target = (tol.confusion() * 1e3).max(1e-4);

    // The blend surface: each station's exact ball arc, the grid closed
    // along the seam by the pinned row ends.
    let mut rows: Vec<Vec<Point>> = vec![Vec::with_capacity(n + 1); ACROSS];
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
    let fitted = ogeom_geom::fit::fit_surface_grid(&rows, 3, fit_target, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the blend surface reached {} against a target of {fit_target}",
            fitted.error
        );
    }

    // The grid's own u-parameters, recomputed the way the fit computes them
    // (averaged centripetal), so the rail pcurves fitted at them are
    // same-parameter with the surface's borders.
    let u_params = averaged_centripetal(&rows);

    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k_count, l_count, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l_count + j] };
    let blend_geo: SurfaceGeometry = surface.into();
    use ogeom_geom::Surface as _;
    let (u_dom, v_dom) = blend_geo.domain();
    let blend_id = model.geometry_mut().add_surface(blend_geo.clone());

    // Rails and seam straight off the control net, exactly as a skin's.
    let border = |j: usize| -> OgeomResult<Curve> {
        let control: Vec<Point> = (0..k_count).map(|i| point_at(i, j)).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?))
    };
    let seam_curve = {
        let control: Vec<Point> = (0..l_count).map(|j| point_at(0, j)).collect();
        Curve::BSpline(ogeom_geom::BSplineCurve::new(v_knots, control, tol)?)
    };
    let rail_first = ogeom_algo::make_edge(model, border(0)?, u_dom, tol)?.shape;
    let rail_second = ogeom_algo::make_edge(model, border(l_count - 1)?, u_dom, tol)?.shape;
    let anchor0 = ogeom_algo::edge_vertices(model, &rail_first)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a rail has no vertex"))?;
    let anchor1 = ogeom_algo::edge_vertices(model, &rail_second)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a rail has no vertex"))?;
    let seam =
        ogeom_algo::make_edge_between(model, seam_curve, v_dom, &anchor0, &anchor1, tol)?.shape;

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
        row_line(v_dom.0)?,
        blend_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &rail_second,
        row_line(v_dom.1)?,
        blend_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        column_line(u_dom.1)?,
        column_line(u_dom.0)?,
        blend_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    let blend_face = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                rail_first.clone(),
                seam.clone(),
                rail_second.reversed(),
                seam.reversed(),
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
    let winding = period_of(host).map_or(0.0, |period| {
        let delta = rail_chart[n - 1].x - rail_chart[0].x + (rail_chart[1].x - rail_chart[0].x);
        // The loop's chart displacement, rounded to whole periods.
        (delta / period).round() * period
    });
    rail_chart.push(Point2::new(rail_chart[0].x + winding, rail_chart[0].y));
    let rail_pcurve = {
        let fitted = ogeom_geom::fit::fit_points_2d_at(u_params, &rail_chart, 3, fit_target, tol)?;
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
    let apex = ogeom_algo::make_edge(model, apex_guide, guide_range, tol)?.shape;
    let apex_pcurve = {
        let fitted =
            ogeom_geom::fit::fit_points_2d_at(&apex_params, &apex_chart, 3, fit_target, tol)?;
        if !fitted.met {
            ogeom_bail!(
                NotDone,
                "the edge's chart image reached {} against a target of {fit_target}",
                fitted.error
            );
        }
        PlanarCurve::from(fitted.curve)
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
fn averaged_centripetal(rows: &[Vec<Point>]) -> Vec<f64> {
    let len = rows[0].len();
    let mut sums = vec![0.0_f64; len];
    for row in rows {
        let mut partial = Vec::with_capacity(len);
        partial.push(0.0);
        let mut total = 0.0;
        for pair in row.windows(2) {
            total += pair[0].distance(pair[1]).sqrt();
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
