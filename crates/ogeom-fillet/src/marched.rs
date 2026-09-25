//! The marched fillet: a rolling-ball blend whose seat has no closed form,
//! carried all the way to topology.
//!
//! [`march_blend`](crate::march_blend) solves the ball's two contact points
//! station by station; this module turns those stations into the same wedge
//! the closed-form blends build. The blend face is a surface fitted through
//! the ball's own arcs, its rails are the surface's own border iso-curves,
//! and the legs are patches of the hosts' *own* surfaces (exact geometry
//! bounded by fitted rails), so the boolean's same-domain resolution melts
//! them exactly as it melts a closed-form wedge's.
//!
//! The chart discipline that makes it sound: the marcher *solves* the
//! contact parameters rather than projecting, so a pcurve fitted through
//! them at the grid's own parameters is same-parameter by construction; the
//! blend chart's rail images are straight iso rows; and the apex ring's
//! images come from closed-form inversion on the host charts, unwrapped for
//! continuity.

use crate::march::{BlendStop, Sides, march_blend_seeded};
use crate::support::{apply_wedge, edge_curve, face_from_edges};
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
    mates: Option<(usize, &[crate::fillet::Mate])>,
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
    // A closed seat on a fitted seam (a bore's rim where it leaves a
    // sphere, the two arcs of a boolean's seam joined end to end) is a
    // loop whose join is a corner: the ends meet, the tangents do not. The
    // march evaluates the guide's derivatives at every step and cannot
    // cross a corner in them; it crawls onto the join and stalls. So the
    // march steers by a *second* guide, the loop re-fitted smooth through
    // its join. It guides the section planes only: the ball seats on the
    // exact hosts, and every edge the wedge builds rides the seat's own
    // curve, so its fit error costs the blend nothing it can measure.
    let closure_of =
        |guide: &Curve, guide_range: (f64, f64)| -> OgeomResult<(bool, Option<Curve>)> {
            let loops = closed || {
                let (lo, hi) = guide_range;
                guide
                    .point_at(lo, tol)
                    .and_then(|p| guide.point_at(hi, tol).map(|q| p.distance(q)))
                    .is_ok_and(|d| d <= tol.confusion() * 10.0)
            };
            let march_guide: Option<Curve> =
                if loops && matches!(guide, Curve::BSpline(b) if !b.is_periodic()) {
                    const SAMPLES: usize = 256;
                    let (lo, hi) = guide_range;
                    let mut points: Vec<Point> = Vec::with_capacity(SAMPLES + 1);
                    for i in 0..=SAMPLES {
                        #[allow(clippy::cast_precision_loss)]
                        let t = lo + (hi - lo) * ((i % SAMPLES) as f64) / (SAMPLES as f64);
                        points.push(guide.point_at(t, tol)?);
                    }
                    let fitted =
                        ogeom_geom::fit::fit_points_closed(&points, 3, tol.confusion() * 1e2, tol)?;
                    fitted.met.then_some(Curve::BSpline(fitted.curve))
                } else {
                    None
                };
            Ok((loops, march_guide))
        };
    let (loops, march_guide) = closure_of(&guide, guide_range)?;

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
        let surface_id = data.surface;
        let Some(stored) = model.geometry().surface(surface_id).cloned() else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        // A fitted patch ends where its face ends, and a ball rolling out
        // through a wall needs the host to go on past it: the patch is
        // continued on every open side by a few radii, as itself, and
        // *written back* as the face's own surface (the same parameters,
        // a wider window, exactly as a reader widens a window to hold a
        // face), so the legs built on it and the face they melt with
        // stand on one chart, as an analytic host's windowed copy does.
        let stored = match stored {
            SurfaceGeometry::BSpline(patch) => {
                use ogeom_geom::Surface as _;
                let reach = radius * 6.0;
                let mut longer = patch;
                for (along_u, closed) in [
                    (true, longer.is_closed_u(tol) || longer.is_periodic_u()),
                    (false, longer.is_closed_v(tol) || longer.is_periodic_v()),
                ] {
                    if closed {
                        continue;
                    }
                    for at_end in [false, true] {
                        if let Ok(grown) = longer.extended(along_u, at_end, reach, 2, tol) {
                            longer = grown;
                        }
                    }
                }
                let wider = SurfaceGeometry::BSpline(longer);
                if let Some(held) = model.geometry_mut().surface_mut(surface_id) {
                    *held = wider.clone();
                }
                wider
            }
            other => other,
        };
        // Baked into the world: the face's surface lives wherever its
        // placement puts it, and everything below (the march, the legs,
        // the melt) speaks world coordinates.
        let surface = {
            use ogeom_geom::Transformable as _;
            let placement = face.transform(model.datums())?;
            stored.transformed(&placement, tol)?
        };
        // Any host the ball can seat on: the analytics invert their charts
        // in closed form, and a fitted patch (or a swept or revolved
        // surface) inverts by projection, warm-started from the last
        // station. A trimmed or offset host is its basis with a story the
        // march does not read.
        if matches!(
            surface,
            SurfaceGeometry::Trimmed(_) | SurfaceGeometry::Offset(_)
        ) {
            ogeom_bail!(
                Construction,
                "a marched fillet's hosts must carry their own chart; a trimmed \
                 or offset host is refused; see docs/PARITY.md, fillet.edge-blends"
            );
        }
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        hosts.push((face, surface, sign));
    }
    let [
        (face_first, first, sign_first),
        (face_second, second, sign_second),
    ] = hosts.as_slice()
    else {
        ogeom_bail!(
            Construction,
            "a marched fillet needs an edge shared by exactly two faces, \
             found {}",
            hosts.len()
        );
    };
    let (face_first, face_second) = (face_first.clone(), face_second.clone());
    let (first, second) = (first.clone(), second.clone());
    let (sign_first, sign_second) = (*sign_first, *sign_second);

    // A seat the boolean split into arcs at its hosts' seams, on a solid
    // whose curves are the arcs themselves (a converted solid, an imported
    // one) has no stored loop to run the whole turn on. The loop is put
    // back through the neighbours the two hosts share, each continuing the
    // last tangentially, and the seat becomes the whole turn as it is
    // where the stored curve carries it.
    let (guide, guide_range, edge_range, loops, march_guide) =
        if !closed && ends_apart(&guide, guide_range, tol) {
            match loop_through_neighbours(
                model,
                edge,
                &guide,
                edge_range,
                [&face_first, &face_second],
                tol,
            )? {
                Some((whole, seat)) => {
                    // The whole turn is a loop where the arc was not: decided
                    // again, or the station on the seam's column is bracketed
                    // the long way round the turn and re-solved on its far side.
                    let domain = whole.domain();
                    let (loops, march_guide) = closure_of(&whole, domain)?;
                    (whole, domain, seat, loops, march_guide)
                }
                None => (guide, guide_range, edge_range, loops, march_guide),
            }
        } else {
            (guide, guide_range, edge_range, loops, march_guide)
        };
    // An open seat on a spline that ends with the edge (a converted
    // solid's every edge) leaves the ball nowhere to run out: the guide
    // is continued past both ends by a few radii, as itself. It steers
    // the section planes only, so its continuation need only be smooth.
    let (guide, guide_range) = match (&guide, closed) {
        (Curve::BSpline(spline), false)
            if !spline.is_periodic() && ends_apart(&guide, guide_range, tol) =>
        {
            let mut longer = spline.clone();
            for at_end in [false, true] {
                if let Ok(grown) = longer.extended(at_end, radius * 8.0, 2, tol) {
                    longer = grown;
                }
            }
            let domain = longer.domain();
            (Curve::BSpline(longer), domain)
        }
        _ => (guide, guide_range),
    };

    // A seat running through a point where its hosts are tangent (the
    // crossing of two equal drums) has no section there: the ball's arc
    // collapses at the pole, and the march can only stall on it. Where a
    // pole is one of the crease's own ends, the band pinches there and the
    // pinched construction builds it; a pole inside a crease with none at
    // its ends is refused by name, sampled along the whole reconstructed
    // loop.
    let pinched = [
        crate::pinched::hosts_tangent_at(
            &first,
            &second,
            stored_end(model, edge, false, tol)?,
            tol,
        )?,
        crate::pinched::hosts_tangent_at(
            &first,
            &second,
            stored_end(model, edge, true, tol)?,
            tol,
        )?,
    ];
    if pinched == [false, false] {
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
                     tangent inside the crease; the ball's section collapses \
                     at that pole and the pinched seam is refused; see \
                     docs/PARITY.md, fillet.edge-blends"
                );
            }
        }
    }

    // Convexity, read from the solid itself the way the planar seat reads
    // it: which way the first face extends from the edge, leaned against
    // the second's outward normal. It decides everything downstream: which
    // of the four ball seatings is the fillet's, and whether the wedge adds
    // or removes.
    //
    // Probed at the *edge's* own midpoint, not the reconstructed loop's: a
    // conic arc re-opened to its full period runs through territory the
    // boolean cut away, and a probe standing off the solid reads nothing in
    // either direction. The edge's midpoint is on the crease by definition.
    use ogeom_geom::Surface as _;
    let mid_t = f64::midpoint(edge_range.0, edge_range.1);
    let convex = crease_convexity(
        model,
        &face_first,
        [(&first, sign_first), (&second, sign_second)],
        &guide,
        edge_range,
        radius,
        tol,
    )?;
    // The fillet's ball rides the material's own side of each support: its
    // centre sits inside the material at a convex corner and out in the
    // notch at a concave one.
    let seat_sign = if convex { -1.0 } else { 1.0 };
    #[allow(clippy::cast_possible_truncation)]
    let sides = Sides {
        first: (seat_sign * sign_first) as i8,
        second: (seat_sign * sign_second) as i8,
    };

    if pinched != [false, false] {
        return crate::pinched::pinched_fillet(
            model,
            solid,
            edge,
            [(&first, sign_first), (&second, sign_second)],
            sides,
            convex,
            radius,
            pinched,
            tol,
        );
    }

    // Seeded at the edge's own midpoint: on a reconstructed loop the domain
    // midpoint may stand in cut-away territory where no ball seats, but the
    // crease's own midpoint is seat by definition, and a closed loop closes
    // from wherever the walker starts.
    let (steering, seed_t) = match &march_guide {
        Some(smooth) => {
            let seed = guide.point_at(mid_t, tol)?;
            let on_smooth = ogeom_algo::project_on_curve(smooth, seed, 256, tol)?;
            (smooth, on_smooth.parameter)
        }
        None => (&guide, mid_t),
    };
    let mut blend = march_blend_seeded(
        &first,
        &second,
        radius,
        steering,
        sides,
        seed_t,
        Marching {
            chord: 3e-6,
            ..Marching::default()
        },
        tol,
    )?;
    if std::env::var_os("OGEOM_DEBUG_RUNOUT").is_some() {
        eprintln!(
            "MARCH stopped {:?} with {} stations; guide domain {:?} closed {closed}; first {:?}/{:?} at {:?}; last {:?}/{:?} at {:?}; hosts {:?} {:?}",
            blend.stopped,
            blend.len(),
            guide.domain(),
            blend.on_first.first(),
            blend.on_second.first(),
            blend.spine.first(),
            blend.on_first.last(),
            blend.on_second.last(),
            blend.spine.last(),
            first.domain(),
            second.domain()
        );
    }
    // Stations steered by the smooth loop carry its parameters; the seat's
    // own curve is what the wedge measures windows on, so each station is
    // re-placed on it by projection.
    if let Some(smooth) = &march_guide {
        for t in &mut blend.along {
            let p = smooth.point_at(*t, tol)?;
            *t = ogeom_algo::project_on_curve(&guide, p, 256, tol)?.parameter;
        }
    }
    // An open seat (the ball ran off the end of a support in each
    // direction) ends in run-out caps instead of closing: the arc-restricted
    // wedge of the revolved fillets, generalised to the fitted band.
    let open_stop = matches!(
        blend.stopped,
        BlendStop::LeftTheFirstSupport
            | BlendStop::LeftTheSecondSupport
            | BlendStop::LeftBothSupports
            | BlendStop::RanPastTheGuide
    );
    if blend.stopped != BlendStop::Closed && !open_stop {
        ogeom_bail!(
            Construction,
            "the blend neither closed on itself nor ran off its supports \
             ({:?}); this stop has no construction yet; see docs/PARITY.md, \
             fillet.edge-blends",
            blend.stopped
        );
    }
    if blend.len() < 8 {
        ogeom_bail!(Construction, "the march produced too few stations to fit");
    }
    if open_stop {
        return open_runout_wedge(
            model,
            solid,
            edge,
            blend,
            &guide,
            edge_range,
            [(&first, sign_first), (&second, sign_second)],
            [&face_first, &face_second],
            mates,
            radius,
            convex,
            tol,
        );
    }
    // The band wants every winding rail to run its period forward; when the
    // march went the other way round, the whole loop reverses.
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
    // reads off the closing step: while it points *against* the march, the
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
                 configuration is still owed; see docs/PARITY.md, \
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
        // iso between the rings, where it cannot cross either of them; a
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
            // The nearest station stands a fraction of a stride off the
            // column; where the column is the host's own seam (a rim the
            // boolean opened at the sphere's meridian) that fraction is a
            // sliver between the rail's start and its seam crossing, which
            // no arrangement holds. Station zero is re-solved exactly on
            // the column, bracketed by its two neighbours.
            let host_is_first = core::ptr::eq(host, &first);
            let u_of = |x: &[f64; 5]| if host_is_first { x[0] } else { x[2] };
            let signed = |u: f64| -> f64 {
                let d = (u - apex_u).rem_euclid(period);
                if d <= period / 2.0 { d } else { d - period }
            };
            let n = blend.len();
            let (lo, hi) = guide_range;
            let span = hi - lo;
            let w0 = blend.along[0];
            let unwrap = |w: f64| -> f64 {
                if !loops {
                    return w;
                }
                let mut w = w;
                while w - w0 > span / 2.0 {
                    w -= span;
                }
                while w0 - w > span / 2.0 {
                    w += span;
                }
                w
            };
            let fold = |w: f64| -> f64 {
                if loops {
                    lo + (w - lo).rem_euclid(span)
                } else {
                    w
                }
            };
            let near = [
                blend.on_first[0].0,
                blend.on_first[0].1,
                blend.on_second[0].0,
                blend.on_second[0].1,
            ];
            let station_at = |w: f64| -> Option<([f64; 5], f64)> {
                let x = crate::march::seat_section(
                    &first,
                    &second,
                    radius,
                    &guide,
                    blend.sides,
                    fold(w),
                    near,
                    tol,
                )
                .ok()?;
                Some((x, signed(u_of(&x))))
            };
            let mut bracket = (unwrap(blend.along[n - 1]), unwrap(blend.along[1]));
            if let (Some((_, fa)), Some((_, fb))) = (station_at(bracket.0), station_at(bracket.1))
                && fa * fb < 0.0
            {
                let (mut fa, mut solved) = (fa, None);
                for _ in 0..60 {
                    let mid = f64::midpoint(bracket.0, bracket.1);
                    let Some((x, fm)) = station_at(mid) else {
                        break;
                    };
                    if fm.abs() <= tol.parametric() || bracket.1 - bracket.0 <= tol.parametric() {
                        solved = Some((x, mid));
                        break;
                    }
                    if fa * fm < 0.0 {
                        bracket.1 = mid;
                    } else {
                        bracket.0 = mid;
                        fa = fm;
                    }
                }
                if let Some((x, w)) = solved
                    && let (Ok(p1), Ok(p2), Ok((du, dv))) = (
                        first.point_at(x[0], x[1], tol),
                        second.point_at(x[2], x[3], tol),
                        first.d1_at(x[0], x[1], tol),
                    )
                {
                    let n1 = du.cross(dv);
                    let centre = p1 + n1 / n1.magnitude() * (f64::from(blend.sides.first) * radius);
                    blend.spine[0] = centre;
                    blend.touch_first[0] = p1;
                    blend.touch_second[0] = p2;
                    blend.on_first[0] = (x[0], x[1]);
                    blend.on_second[0] = (x[2], x[3]);
                    blend.along[0] = fold(w);
                }
            }
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
    // Two tenths of a micron at unit scale: what the march's own stations
    // hold to, and what the band's edges are widened to say.
    let fit_target = (tol.confusion() * 2e3).max(2e-4);

    // The blend surface: each station's exact ball arc, the loop of stations
    // fitted *closed*: the join is C1 wherever the seam lands, so anchoring
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
    // compares them against exact geometry (the melt's crossing paver above
    // all) widens by an edge's recorded tolerance, not by wishful thinking.
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
        // cutting (the wedge is the corner material the ball displaced)
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
    let fitted_host = matches!(first, SurfaceGeometry::BSpline(_))
        || matches!(second, SurfaceGeometry::BSpline(_));
    match apply_wedge(model, solid, Some(edge), &faces, additive, tol) {
        Err(e) if fitted_host => ogeom_bail!(
            NotDone,
            "the blend marched and its wedge was built, but the melt against a \
             fitted host is beyond what the boolean resolves ({e}); see \
             docs/PARITY.md, fillet.edge-blends"
        ),
        other => other,
    }
}

/// The open seat's wedge: a marched band that ran off its supports, capped
/// at both ends: the revolved fillets' arc-restricted wedge generalised to
/// the fitted band.
///
/// Five faces close it: the band fitted *open* through the ball's arcs, one
/// leg on each host between the crease and the touch rail, and one planar
/// cap in the section plane of each end station: the plane the march's own
/// fourth equation held every section to, so the end arc, both touch points
/// and the crease point all stand in it by construction.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
#[allow(clippy::too_many_lines, reason = "one wedge, assembled end to end")]
fn open_runout_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    mut blend: crate::march::MarchedBlend,
    guide: &Curve,
    edge_range: (f64, f64),
    hosts: [(&SurfaceGeometry, f64); 2],
    host_faces: [&Shape; 2],
    mates: Option<(usize, &[crate::fillet::Mate])>,
    radius: f64,
    convex: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Surface as _;
    let [(first, _), (second, _)] = hosts;

    // The walker clamps the guide parameter into its window; a run that
    // crossed the period comes back wrapped. Unwrap it into one monotonic
    // sweep, turn the whole band forward, and drop any station that fails
    // to advance: the seed join can hand back a duplicate.
    if guide.is_periodic() {
        let (lo, hi) = guide.domain();
        let period = hi - lo;
        for i in 1..blend.along.len() {
            let mut t = blend.along[i];
            while t - blend.along[i - 1] > period / 2.0 {
                t -= period;
            }
            while blend.along[i - 1] - t > period / 2.0 {
                t += period;
            }
            blend.along[i] = t;
        }
    }
    if blend.along.last() < blend.along.first() {
        blend.spine.reverse();
        blend.touch_first.reverse();
        blend.touch_second.reverse();
        blend.on_first.reverse();
        blend.on_second.reverse();
        blend.along.reverse();
    }
    let mut i = 1;
    while i < blend.len() {
        if blend.along[i] <= blend.along[i - 1] && blend.len() > 8 {
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
    // The seat the fillet owns is the *edge's* window, not everywhere the
    // supports happen to extend: the surfaces run on past the solid (a box
    // face's plane does not end at the box) and the march runs with them
    // into territory the boolean cut away. Trim the band to the edge's own
    // window, and solve the exact section at each end: the caps stand on
    // those, not on wherever the walker's last step landed.
    {
        let span = edge_range.1 - edge_range.0;
        if span <= 0.0 {
            ogeom_bail!(Construction, "the edge's window has no length");
        }
        // The unwrapped run lives on its own branch of a periodic guide;
        // shift the edge window onto it before comparing parameters.
        let mut w0 = edge_range.0;
        if guide.is_periodic() {
            let (lo, hi) = guide.domain();
            let period = hi - lo;
            let mid_run = f64::midpoint(blend.along[0], blend.along[blend.len() - 1]);
            let k = ((mid_run - (w0 + span / 2.0)) / period).round();
            w0 += k * period;
        }
        let mut w1 = w0 + span;
        // Where the crease *terminates* at the solid's own boundary (its end
        // vertex belongs to a third face, not to a continuation of the seat
        // past a seam split), the blend runs out through the wall: the band
        // carries on past the window until the ball's contacts have left
        // both host faces, and the cut trims the wedge against whatever the
        // crease ended on. Capped in its own arc plane at the crease's end
        // instead, the wedge stops short of the wall by the plane's slant
        // and leaves a sliver of sharp crease between cap and wall, the
        // remnant a second fillet then has to meet.
        if std::env::var_os("OGEOM_DEBUG_RUNOUT").is_some() {
            eprintln!(
                "RUNOUT start: edge_range {edge_range:?} w0 {w0:.5} w1 {w1:.5} run {:.5}..{:.5} guide domain {:?}",
                blend.along[0],
                blend.along[blend.len() - 1],
                guide.domain()
            );
        }
        for end in [true, false] {
            // The window may stand a period past a periodic guide's stored
            // domain; the point is the same turn either way.
            let at = {
                let mut w = if end { w0 } else { w1 };
                if guide.is_periodic() {
                    let (lo, hi) = guide.domain();
                    w = lo + (w - lo).rem_euclid(hi - lo);
                }
                guide.point_at(w, tol)?
            };
            let terminates = crease_terminates_at(model, solid, edge, host_faces, at, tol)?;
            if std::env::var_os("OGEOM_DEBUG_RUNOUT").is_some() {
                eprintln!("RUNOUT end {end} at {at:?} terminates {terminates}");
            }
            if !terminates {
                continue;
            }
            // A mate of the same request ending at this vertex some way
            // other than tangentially is a corner mate: the band runs on
            // through that blend until the ball has left the material, and
            // the cut trims the two against each other. A chain mate,
            // leaving the vertex the way this crease arrives, is a junction
            // the caps close flush.
            let outward_here = {
                let d = guide.d1_at(
                    {
                        let mut w = if end { w0 } else { w1 };
                        if guide.is_periodic() {
                            let (lo, hi) = guide.domain();
                            w = lo + (w - lo).rem_euclid(hi - lo);
                        }
                        w
                    },
                    tol,
                )?;
                let m = d.magnitude();
                let unit = if m > tol.confusion() {
                    d / m
                } else {
                    Vector::ZERO
                };
                if end { -unit } else { unit }
            };
            let (chain_mate, corner_mate) = match mates {
                None => (false, false),
                Some((index, mates)) => {
                    let mut chain = false;
                    let mut corner = false;
                    for (i, mate) in mates.iter().enumerate() {
                        if i == index {
                            continue;
                        }
                        for (p, leaving) in &mate.ends {
                            if p.distance(at) > tol.confusion() * 1e3 {
                                continue;
                            }
                            if leaving.cross(outward_here).magnitude() <= 1e-2
                                && leaving.dot(outward_here) > 0.0
                            {
                                chain = true;
                            } else {
                                corner = true;
                            }
                        }
                    }
                    (chain, corner && !chain)
                }
            };
            let settled = mates.is_some_and(|(index, mates)| {
                crate::fillet::Mate::settled_at(mates, index, at, tol)
            });
            if chain_mate || settled {
                continue;
            }
            // Where the ball's contacts stand against the host faces past
            // the window: `Out` of both is clear of the solid, `On` either
            // is a neighbouring blend's own rail; the seat goes on under
            // that blend, which a cap in the section's own plane meets.
            let standing = |i: usize| -> OgeomResult<(bool, bool)> {
                let deflection = ogeom_mesh::Deflection {
                    chord: (radius * 1e-3).max(tol.confusion() * 1e2),
                    ..ogeom_mesh::Deflection::default()
                };
                let first = ogeom_algo::classify_on_face(
                    model,
                    host_faces[0],
                    blend.touch_first[i],
                    deflection,
                    tol,
                )?;
                let second = ogeom_algo::classify_on_face(
                    model,
                    host_faces[1],
                    blend.touch_second[i],
                    deflection,
                    tol,
                )?;
                use ogeom_algo::Containment as C;
                Ok((
                    first == C::Out && second == C::Out,
                    first == C::On || second == C::On,
                ))
            };
            let n = blend.len();
            // Stations outside the window, nearest the window first. The
            // first clear one is where the wedge has left the solid; the
            // run-out carries on a radius further so the wall crosses the
            // band well inside it (a crossing a hair from the band's end
            // is one the intersector's seeding can miss) and the run's own
            // end is the honest stop where the walker quit first.
            let outside: Vec<usize> = if end {
                (0..n).rev().filter(|&i| blend.along[i] < w0).collect()
            } else {
                (0..n).filter(|&i| blend.along[i] > w1).collect()
            };
            let mut cleared: Option<Point> = None;
            // A contact standing on a host's boundary for one station is
            // the contact crossing a wall's edge; standing there station
            // after station, it is riding a neighbouring blend's rail.
            let mut on_since: Option<Point> = None;
            for i in outside {
                let Some(from) = cleared else {
                    let (is_clear, on_edge) = standing(i)?;
                    if on_edge {
                        let since = *on_since.get_or_insert(blend.spine[i]);
                        if blend.spine[i].distance(since) > radius * 0.05 {
                            break;
                        }
                    } else {
                        on_since = None;
                    }
                    if is_clear {
                        // Clear of the host faces is not clear of the solid.
                        // Past a seam vertex whose other half is blended
                        // already, the hosts are cut away too, yet the ball
                        // still sits in the material there: the seat goes
                        // on under the neighbouring blend, which is exactly
                        // where a cap in the section's own plane meets it.
                        // Only a ball that has left the material altogether
                        // (a convex seat's, past a wall) runs out.
                        let deflection = ogeom_mesh::Deflection {
                            chord: (radius * 1e-2).max(tol.confusion() * 1e3),
                            ..ogeom_mesh::Deflection::default()
                        };
                        let inside = ogeom_algo::classify_in_solid(
                            model,
                            solid,
                            blend.spine[i],
                            deflection,
                            tol,
                        )? == ogeom_algo::Containment::In;
                        if inside == convex && !corner_mate {
                            break;
                        }
                        if inside == convex {
                            // Under the corner mate's blend still: walk on.
                            continue;
                        }
                        cleared = Some(blend.spine[i]);
                    }
                    continue;
                };
                if end {
                    w0 = blend.along[i];
                } else {
                    w1 = blend.along[i];
                }
                if blend.spine[i].distance(from) >= radius {
                    break;
                }
            }
        }
        // Only cap at an end the march actually reached past; where it
        // stopped short (the true seat ended first), the walker's own last
        // station is the honest end.
        let cap0 = blend.along[0] < w0;
        let cap1 = blend.along[blend.len() - 1] > w1;
        if std::env::var_os("OGEOM_DEBUG_RUNOUT").is_some() {
            eprintln!(
                "RUNOUT edge_range {edge_range:?} window ({w0:.5}, {w1:.5}) run {:.5}..{:.5} cap0 {cap0} cap1 {cap1}",
                blend.along[0],
                blend.along[blend.len() - 1]
            );
        }
        let mut keep_from = 0;
        let mut keep_to = blend.len();
        for (i, t) in blend.along.iter().enumerate() {
            if *t <= w0 {
                keep_from = i + 1;
            }
            if *t >= w1 && keep_to == blend.len() {
                keep_to = i;
            }
        }
        if keep_from >= keep_to {
            ogeom_bail!(
                Construction,
                "the marched band and the edge's window do not overlap; the \
                 guide does not run along this seat"
            );
        }
        let cut = |v: &mut Vec<Point>, from: usize, to: usize| {
            v.truncate(to);
            v.drain(..from);
        };
        cut(&mut blend.spine, keep_from, keep_to);
        cut(&mut blend.touch_first, keep_from, keep_to);
        cut(&mut blend.touch_second, keep_from, keep_to);
        blend.on_first.truncate(keep_to);
        blend.on_first.drain(..keep_from);
        blend.on_second.truncate(keep_to);
        blend.on_second.drain(..keep_from);
        blend.along.truncate(keep_to);
        blend.along.drain(..keep_from);
        let mut end_station = |w: f64, front: bool| -> OgeomResult<()> {
            // Seeded from the adjacent kept station: the Newton must settle
            // in *this* seat's basin: a drum's far side holds a ball too.
            let i = if front { 0 } else { blend.len() - 1 };
            let near = [
                blend.on_first[i].0,
                blend.on_first[i].1,
                blend.on_second[i].0,
                blend.on_second[i].1,
            ];
            let x = crate::march::seat_section(
                first,
                second,
                radius,
                guide,
                blend.sides,
                w,
                near,
                tol,
            )?;
            let p1 = first.point_at(x[0], x[1], tol)?;
            let p2 = second.point_at(x[2], x[3], tol)?;
            let n1 = {
                let (du, dv) = first.d1_at(x[0], x[1], tol)?;
                let n = du.cross(dv);
                n / n.magnitude()
            };
            let centre = p1 + n1 * (f64::from(blend.sides.first) * radius);
            if front {
                blend.spine.insert(0, centre);
                blend.touch_first.insert(0, p1);
                blend.touch_second.insert(0, p2);
                blend.on_first.insert(0, (x[0], x[1]));
                blend.on_second.insert(0, (x[2], x[3]));
                blend.along.insert(0, w);
            } else {
                blend.spine.push(centre);
                blend.touch_first.push(p1);
                blend.touch_second.push(p2);
                blend.on_first.push((x[0], x[1]));
                blend.on_second.push((x[2], x[3]));
                blend.along.push(w);
            }
            Ok(())
        };
        if cap0 {
            end_station(w0, true)?;
        }
        if cap1 {
            end_station(w1, false)?;
        }
    }
    if blend.len() < 8 {
        ogeom_bail!(
            Construction,
            "the edge's window holds too few marched stations to fit"
        );
    }

    // The forward and backward halves join at the seed, and the walker's
    // landing on each boundary comes in short refining steps: both leave
    // stations standing nearly on top of a neighbour, and a cramped pair
    // puts two grid columns in the same place and poisons the fit's
    // parameterization. One sweep, the closed loop's own remedy.
    {
        let mut steps: Vec<f64> = blend
            .spine
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .collect();
        steps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        let cramped = steps.get(steps.len() / 2).copied().unwrap_or(0.0) * 0.25;
        let mut i = 1;
        while i < blend.len() {
            // The two boundary stations are the run-out itself; the sweep
            // drops their crowding neighbours, never the ends.
            let at_end = i == blend.len() - 1;
            let victim = if at_end { i - 1 } else { i };
            if victim > 0
                && blend.spine[i].distance(blend.spine[i - 1]) <= cramped
                && blend.len() > 8
            {
                blend.spine.remove(victim);
                blend.touch_first.remove(victim);
                blend.touch_second.remove(victim);
                blend.on_first.remove(victim);
                blend.on_second.remove(victim);
                blend.along.remove(victim);
            } else {
                i += 1;
            }
        }
    }

    build_open_band(
        model,
        solid,
        edge,
        &blend,
        guide,
        hosts,
        radius,
        convex,
        [false, false],
        tol,
    )
}

/// The band and wedge over an open run of stations that spans exactly the
/// window the wedge owns, applied to the solid.
///
/// An end that is `pinched` stands on a pole where the hosts turn tangent:
/// its station's touch points and crease point are the one pole, its row
/// collapses to it, and the wedge closes there without a cap, the band's
/// end arc a degenerate edge and both legs meeting the rails at the pole.
/// Every other end is capped in its section's own plane.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
#[allow(clippy::too_many_lines, reason = "one wedge, assembled end to end")]
pub(crate) fn build_open_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    blend: &crate::march::MarchedBlend,
    guide: &Curve,
    hosts: [(&SurfaceGeometry, f64); 2],
    radius: f64,
    convex: bool,
    pinched: [bool; 2],
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Surface as _;
    let [(first, sign_first), (second, sign_second)] = hosts;
    let additive = !convex;
    let n = blend.len();
    // Two tenths of a micron at unit scale: what the march's own stations
    // hold to, and what the band's edges are widened to say.
    let fit_target = (tol.confusion() * 2e3).max(2e-4);

    // The band: each station's exact ball arc, fitted open along the
    // stations: same arcs as the closed case, no wrap. Sampled twice as
    // finely across: the scoop's widest sections sweep well past a right
    // angle, and the knot refinement can only split spans that still hold
    // data.
    const ACROSS_OPEN: usize = 17;
    let mut rows: Vec<Vec<Point>> = (0..ACROSS_OPEN).map(|_| Vec::with_capacity(n)).collect();
    for at in 0..n {
        if (at == 0 && pinched[0]) || (at == n - 1 && pinched[1]) {
            for row in &mut rows {
                row.push(blend.touch_first[at]);
            }
            continue;
        }
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
            let theta = sweep * (k as f64) / ((ACROSS_OPEN - 1) as f64);
            let dir = a * theta.cos() + axis.cross(a) * theta.sin();
            row.push(centre + dir * radius);
        }
    }
    let stations: Vec<Vec<Point>> = (0..n)
        .map(|i| (0..ACROSS_OPEN).map(|k| rows[k][i]).collect())
        .collect();
    // Fitted at half the target: the two passes each hold their half, but
    // the assembled surface's honest error is measured across both, and an
    // open band lands near their sum where the closed band's wrap absorbs
    // it. The acceptance stays the caller's target.
    let fitted = ogeom_geom::fit::fit_surface_grid_chordal(&stations, 3, fit_target * 0.5, tol)?;
    if fitted.error > fit_target {
        ogeom_bail!(
            NotDone,
            "the blend surface reached {} against a target of {fit_target}",
            fitted.error
        );
    }
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

    // Six shared vertices: the four band corners off the control net (
    // which the clamped borders interpolate exactly) and the crease's two
    // ends off the guide itself.
    // The cap at each end stands in the end section's *own* plane: the
    // plane of the ball's arc, which holds both touch points exactly. The
    // march's guide condition holds only one point of the section to the
    // guide-normal plane, so that plane cannot close a cap. The cap's apex
    // is where the crease crosses the arc plane: a whisker off the end
    // parameter, found by one-dimensional Newton along the guide.
    let section_plane = |end: usize| -> OgeomResult<ogeom_math::Plane> {
        let centre = blend.spine[end];
        let a = (blend.touch_first[end] - centre) / radius;
        let b = (blend.touch_second[end] - centre) / radius;
        let n = a.cross(b);
        let m = n.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "a run-out section has no plane to cap in");
        }
        Ok(ogeom_math::Plane::through(
            centre,
            ogeom_math::Direction::new(n / m, tol)?,
        ))
    };
    let apex_on = |plane: &ogeom_math::Plane, near: f64| -> OgeomResult<f64> {
        let mut t = near;
        for _ in 0..40 {
            let f = plane.signed_distance_to(guide.point_at(t, tol)?);
            if f.abs() <= tol.confusion() * 0.01 {
                return Ok(t);
            }
            let df = plane.normal().vector().dot(guide.d1_at(t, tol)?);
            if df.abs() <= 1e-12 {
                break;
            }
            t -= f / df;
        }
        ogeom_bail!(
            Construction,
            "the crease does not cross a run-out cap's section plane"
        );
    };
    // A pinched end has no section plane: its apex is the pole itself,
    // at the station's own parameter.
    let plane0 = if pinched[0] {
        None
    } else {
        Some(section_plane(0)?)
    };
    let plane1 = if pinched[1] {
        None
    } else {
        Some(section_plane(n - 1)?)
    };
    let (t0, t1) = (
        match &plane0 {
            Some(plane) => apex_on(plane, blend.along[0])?,
            None => blend.along[0],
        },
        match &plane1 {
            Some(plane) => apex_on(plane, blend.along[n - 1])?,
            None => blend.along[n - 1],
        },
    );
    if t1 <= t0 {
        ogeom_bail!(Construction, "the run-out caps cross; nothing to blend");
    }
    let apex0 = guide.point_at(t0, tol)?;
    let apex1 = guide.point_at(t1, tol)?;
    let corner = |i: usize, j: usize| point_at(i, j);
    let va0 = ogeom_algo::make_vertex(model, apex0).shape;
    let va1 = ogeom_algo::make_vertex(model, apex1).shape;
    // At a pinched end the crease, both rails and the collapsed row all
    // stand on the pole: one vertex.
    let vc00 = if pinched[0] {
        va0.clone()
    } else {
        ogeom_algo::make_vertex(model, corner(0, 0)).shape
    };
    let vc01 = if pinched[1] {
        va1.clone()
    } else {
        ogeom_algo::make_vertex(model, corner(0, l_count - 1)).shape
    };
    let vc10 = if pinched[0] {
        va0.clone()
    } else {
        ogeom_algo::make_vertex(model, corner(k_count - 1, 0)).shape
    };
    let vc11 = if pinched[1] {
        va1.clone()
    } else {
        ogeom_algo::make_vertex(model, corner(k_count - 1, l_count - 1)).shape
    };
    // The corners stand a fit error from the exact touch points the
    // connectors end at; the vertices own that slop.
    for v in [&vc00, &vc01, &vc10, &vc11] {
        model.widen(v, ogeom_core::Tolerance::new(fit_target)?)?;
    }

    // The band's borders, straight off the control net.
    let border = |i: usize| -> OgeomResult<Curve> {
        let control: Vec<Point> = (0..l_count).map(|j| point_at(i, j)).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?))
    };
    let end_arc = |j: usize| -> OgeomResult<Curve> {
        let control: Vec<Point> = (0..k_count).map(|i| point_at(i, j)).collect();
        Ok(Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?))
    };
    let rail_first =
        ogeom_algo::make_edge_between(model, border(0)?, v_dom, &vc00, &vc01, tol)?.shape;
    let rail_second =
        ogeom_algo::make_edge_between(model, border(k_count - 1)?, v_dom, &vc10, &vc11, tol)?.shape;
    // A pinched end's arc is the collapsed row: a degenerate edge, which
    // carries its chart image and no curve.
    let collapsed = |model: &mut Model, at: &Shape| -> OgeomResult<Shape> {
        let mut data = ogeom_topo::EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[at.clone(), at.clone()])
    };
    let arc_start = if pinched[0] {
        collapsed(model, &vc00)?
    } else {
        ogeom_algo::make_edge_between(model, end_arc(0)?, u_dom, &vc00, &vc10, tol)?.shape
    };
    let arc_end = if pinched[1] {
        collapsed(model, &vc01)?
    } else {
        ogeom_algo::make_edge_between(model, end_arc(l_count - 1)?, u_dom, &vc01, &vc11, tol)?.shape
    };
    for rail in [&rail_first, &rail_second, &arc_start, &arc_end] {
        if let Some(node) = model.node_mut(rail)
            && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
        {
            data.tolerance = data.tolerance.widen_to(fit_target);
        }
    }

    // The crease over the marched window. On a periodic guide the unwrapped
    // window may stand past the stored domain; the curve is the same turn
    // either way, so shift it back in.
    let apex_edge = {
        let (lo, hi) = guide.domain();
        let mut window = (t0, t1);
        if guide.is_periodic() {
            let period = hi - lo;
            while window.1 > hi {
                window = (window.0 - period, window.1 - period);
            }
            while window.0 < lo {
                window = (window.0 + period, window.1 + period);
            }
        }
        ogeom_algo::make_edge_between(model, guide.clone(), window, &va0, &va1, tol)?.shape
    };

    // One connector per host per end: the host's own section in the cap's
    // plane, from the crease to the touch rail. A pinched end has none.
    let connector = |model: &mut Model,
                     host: &SurfaceGeometry,
                     plane: Option<ogeom_math::Plane>,
                     apex: (&Shape, Point),
                     touch: (&Shape, Point)|
     -> OgeomResult<Option<Shape>> {
        match plane {
            Some(plane) => Ok(Some(section_connector(
                model, host, plane, apex, touch, radius, fit_target, tol,
            )?)),
            None => Ok(None),
        }
    };
    let conn_first_0 = connector(model, first, plane0, (&va0, apex0), (&vc00, corner(0, 0)))?;
    let conn_first_1 = connector(
        model,
        first,
        plane1,
        (&va1, apex1),
        (&vc01, corner(0, l_count - 1)),
    )?;
    let conn_second_0 = connector(
        model,
        second,
        plane0,
        (&va0, apex0),
        (&vc10, corner(k_count - 1, 0)),
    )?;
    let conn_second_1 = connector(
        model,
        second,
        plane1,
        (&va1, apex1),
        (&vc11, corner(k_count - 1, l_count - 1)),
    )?;

    // The band face: same-parameter iso pcurves on its own chart, no seam.
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
    ogeom_algo::attach_pcurve(
        model,
        &arc_start,
        row_line(v_dom.0)?,
        blend_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &arc_end,
        row_line(v_dom.1)?,
        blend_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    let blend_face = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                arc_start.clone(),
                rail_second.clone(),
                arc_end.reversed(),
                rail_first.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, blend_id, std::slice::from_ref(&wire), tol)?.shape;
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

    // The legs: each host's own surface between the crease and its rail,
    // closed at the ends by the connectors.
    let leg = |model: &mut Model,
               host: &SurfaceGeometry,
               rail: &Shape,
               ends: [&Option<Shape>; 2]|
     -> OgeomResult<Shape> {
        let mut edges = vec![apex_edge.clone()];
        edges.extend(ends[1].clone());
        edges.push(rail.reversed());
        edges.extend(ends[0].as_ref().map(Shape::reversed));
        face_from_edges(model, host.clone(), &edges, tol)
    };
    let leg_first = leg(model, first, &rail_first, [&conn_first_0, &conn_first_1])?;
    let leg_second = leg(
        model,
        second,
        &rail_second,
        [&conn_second_0, &conn_second_1],
    )?;

    // The caps: the section planes, bounded by connector–arc–connector.
    // Outward is out of the marched window: against the guide at the
    // start, along it at the end.
    let cap = |model: &mut Model,
               plane: ogeom_math::Plane,
               edges: &[Shape],
               outward: Vector|
     -> OgeomResult<Shape> {
        let reach = (radius * 4.0).max(1.0);
        let surface: SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let face = face_from_edges(model, surface, edges, tol)?;
        if plane.normal().vector().dot(outward) > 0.0 {
            Ok(face)
        } else {
            Ok(face.reversed())
        }
    };
    let mut caps = Vec::new();
    if let (Some(plane), Some(c1), Some(c2)) = (plane0, &conn_first_0, &conn_second_0) {
        caps.push(cap(
            model,
            plane,
            &[c1.clone(), arc_start.clone(), c2.reversed()],
            -guide.d1_at(t0, tol)?,
        )?);
    }
    if let (Some(plane), Some(c1), Some(c2)) = (plane1, &conn_first_1, &conn_second_1) {
        caps.push(cap(
            model,
            plane,
            &[c1.clone(), arc_end.clone(), c2.reversed()],
            guide.d1_at(t1, tol)?,
        )?);
    }

    let orient = |face: Shape, host_sign: f64| -> Shape {
        let aligned = if additive { -host_sign } else { host_sign };
        if aligned > 0.0 { face } else { face.reversed() }
    };
    let mut faces = vec![
        orient(leg_first, sign_first),
        orient(leg_second, sign_second),
        blend_face,
    ];
    faces.extend(caps);
    apply_wedge(model, solid, Some(edge), &faces, additive, tol)
}

/// Whether the crease between the two hosts is convex, read from the solid
/// itself the way the planar seat reads it: which way the first face
/// extends from the edge, leaned against the second's outward normal. It
/// decides which of the four ball seatings is the fillet's, and whether
/// the wedge adds or removes.
///
/// Probed at the *edge's* own midpoint, not a reconstructed loop's: a
/// conic arc re-opened to its full period runs through territory the
/// boolean cut away, and a probe standing off the solid reads nothing in
/// either direction. The edge's midpoint is on the crease by definition.
pub(crate) fn crease_convexity(
    model: &Model,
    face_first: &Shape,
    hosts: [(&SurfaceGeometry, f64); 2],
    guide: &Curve,
    edge_range: (f64, f64),
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    let [(first, sign_first), (second, sign_second)] = hosts;
    let mid_t = f64::midpoint(edge_range.0, edge_range.1);
    let mid = guide.point_at(mid_t, tol)?;
    let outward_at = |surface: &SurfaceGeometry, sign: f64| -> OgeomResult<Vector> {
        let projection = ogeom_algo::project_on_surface(surface, mid, 32, tol)?;
        let (u, v) = projection.parameters;
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let n = du.cross(dv);
        Ok(n / n.magnitude() * sign)
    };
    let n1 = outward_at(first, sign_first)?;
    let n2 = outward_at(second, sign_second)?;
    {
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
                let probe = ogeom_algo::project_on_surface(first, mid + dir * eps, 32, tol)?.point;
                if ogeom_algo::classify_on_face(model, face_first, probe, deflection, tol)?
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
        Ok(lean < 0.0)
    }
}

/// The start or end of the edge's own stored curve over its own range,
/// placed: the crease's ends, whatever loop the seat was rebuilt on.
fn stored_end(model: &Model, edge: &Shape, end: bool, tol: Tolerances) -> OgeomResult<Point> {
    let (curve, (lo, hi)) = edge_curve(model, edge, tol)?;
    curve.point_at(if end { hi } else { lo }, tol)
}

/// Whether the crease ends at `at` because the solid does (its end vertex
/// there belongs to a third face) rather than continuing as another edge
/// on the same two hosts past a split.
pub(crate) fn crease_terminates_at(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    host_faces: [&Shape; 2],
    at: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(true);
    };
    let point_of = |v: &Shape| -> OgeomResult<Point> {
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex().cloned()) else {
            ogeom_bail!(Construction, "vertex node holds no point");
        };
        Ok(v.transform(model.datums())?.apply(data.point))
    };
    let vertex = if point_of(&a)?.distance(at) <= point_of(&b)?.distance(at) {
        a
    } else {
        b
    };
    let _ = tol;
    for other in explore(model, solid, Filter::OfType(ShapeType::Edge))? {
        if other.is_same(edge) {
            continue;
        }
        let shares = ogeom_algo::edge_vertices(model, &other)?
            .is_some_and(|(x, y)| x.is_same(&vertex) || y.is_same(&vertex));
        if !shares {
            continue;
        }
        // The crease continues where another edge on both hosts leaves the
        // vertex.
        let on_both = host_faces.iter().all(|face| {
            explore(model, face, Filter::OfType(ShapeType::Edge))
                .is_ok_and(|edges| edges.iter().any(|e| e.is_same(&other)))
        });
        if on_both {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether a face other than the hosts meets the crease's end at `at`
/// *tangentially* to a host, and rounds the way a seat of this `convex`
/// does: a neighbouring blend's band, which is tangent to the host it
/// rides. Material under such a face is the neighbour's rounding, not a
/// step: a blend meeting it runs on through it and the cut trims the two
/// bands against each other.
///
/// Only where the two round the same way. A convex blend's cut trims a
/// convex band against its own, but a fill in a re-entrant corner is
/// material that cut would eat (an L-bracket's front edge filleted after
/// its re-entrant one ran the whole length of the leg and out the far side
/// of the wall) and a fill cannot run on through a band either. `radius`
/// sets the chord the band's own side is read over.
pub(crate) fn neighbour_blend_at(
    model: &Model,
    solid: &Shape,
    host_faces: [&Shape; 2],
    at: Point,
    convex: bool,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let host_normals: Vec<Vector> = host_faces
        .iter()
        .filter_map(|f| face_normal_near(model, f, at, tol).ok().flatten())
        .collect();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        if host_faces.iter().any(|h| h.is_same(&face)) {
            continue;
        }
        let touches = explore(model, &face, Filter::OfType(ShapeType::Vertex))?
            .iter()
            .any(|v| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex().cloned())
                    .and_then(|d| v.transform(model.datums()).ok().map(|t| t.apply(d.point)))
                    .is_some_and(|p| p.distance(at) <= tol.confusion() * 1e3)
            });
        if !touches {
            continue;
        }
        let Some(normal) = face_normal_near(model, &face, at, tol)? else {
            continue;
        };
        if host_normals
            .iter()
            .any(|h| h.cross(normal).magnitude() <= 1e-2)
            && band_rounds_out_near(model, &face, at, radius * 0.5, tol)? == Some(convex)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Which way a face rounds where it comes nearest `at`.
///
/// `Some(true)` for a band whose material bulges out (a chord between two
/// of its points sinks below the surface, into the solid) and
/// `Some(false)` for one that fills a re-entrant corner, where the chord
/// stands proud of it. `None` for a face flat over `step` at `at`, which
/// has no side to be on, or one whose normal cannot be read there.
///
/// Which side is out is the face's own orientation over its surface's
/// normal, the same reading [`ogeom_algo::face_normal`] takes, so nothing
/// is classified and no ray is cast.
fn band_rounds_out_near(
    model: &Model,
    face: &Shape,
    at: Point,
    step: f64,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    use ogeom_geom::Surface as _;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(stored) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let surface = {
        use ogeom_geom::Transformable as _;
        stored.transformed(&face.transform(model.datums())?, tol)?
    };
    let projection = ogeom_algo::project_on_surface(&surface, at, 24, tol)?;
    if projection.distance > tol.confusion() * 1e3 {
        return Ok(None);
    }
    let (u, v) = projection.parameters;
    let outward_sign = if face.orientation() == ogeom_topo::Orientation::Reversed {
        -1.0
    } else {
        1.0
    };
    let ((ulo, uhi), (vlo, vhi)) = surface.domain();
    // The sagitta of a short chord each way across the surface, the larger
    // of the two answering: a plane has neither, and on a band the way
    // across it dwarfs the way along.
    let mut sagitta = 0.0_f64;
    for across in [true, false] {
        let (du, dv) = surface.d1_at(u, v, tol)?;
        let length = if across { du } else { dv }.magnitude();
        if length <= tol.angular() {
            continue;
        }
        let (lo, hi, here) = if across { (ulo, uhi, u) } else { (vlo, vhi, v) };
        let span = hi - lo;
        let mut half = step / length;
        if span.is_finite() && half * 2.0 > span {
            half = span * 0.5;
        }
        if half <= tol.parametric() {
            continue;
        }
        // The chord's own middle, slid inside the domain: at a patch's
        // corner (which is where a band's rail ends, and so where this is
        // asked), a chord centred on the projection falls half off the
        // surface, and a band read as flat is a band run on through.
        let mut middle = here;
        if lo.is_finite() {
            middle = middle.max(lo + half);
        }
        if hi.is_finite() {
            middle = middle.min(hi - half);
        }
        let at_parameter = |t: f64| -> (f64, f64) { if across { (t, v) } else { (u, t) } };
        let ends = [middle - half, middle, middle + half]
            .map(at_parameter)
            .map(|(a, b)| surface.point_at(a, b, tol));
        let [Ok(low), Ok(mid), Ok(high)] = ends else {
            continue;
        };
        let (mdu, mdv) = {
            let (a, b) = at_parameter(middle);
            surface.d1_at(a, b, tol)?
        };
        let raw = mdu.cross(mdv);
        let magnitude = raw.magnitude();
        if magnitude <= tol.angular() {
            continue;
        }
        let found = (low.midpoint(high) - mid).dot(raw / magnitude * outward_sign);
        if found.abs() > sagitta.abs() {
            sagitta = found;
        }
    }
    if sagitta.abs() <= tol.confusion() * 10.0 {
        return Ok(None);
    }
    Ok(Some(sagitta < 0.0))
}

/// A face's surface normal (either sign) at the point of it nearest `at`.
pub(crate) fn face_normal_near(
    model: &Model,
    face: &Shape,
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    use ogeom_geom::Surface as _;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(stored) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let surface = {
        use ogeom_geom::Transformable as _;
        stored.transformed(&face.transform(model.datums())?, tol)?
    };
    let projection = ogeom_algo::project_on_surface(&surface, at, 24, tol)?;
    if projection.distance > tol.confusion() * 1e3 {
        return Ok(None);
    }
    let (u, v) = projection.parameters;
    let (du, dv) = surface.d1_at(u, v, tol)?;
    let n = du.cross(dv);
    let m = n.magnitude();
    if m <= tol.angular() {
        return Ok(None);
    }
    Ok(Some(n / m))
}

/// The host's own curve in a cap's section plane, from the crease vertex to
/// the touch-rail corner: a segment on a planar host, the exact conic the
/// plane cuts from a drum, always the short way round.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn section_connector(
    model: &mut Model,
    host: &SurfaceGeometry,
    plane: ogeom_math::Plane,
    from: (&Shape, Point),
    to: (&Shape, Point),
    radius: f64,
    slack: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if matches!(host, SurfaceGeometry::Plane(_)) {
        return crate::support::segment_between(model, from, to, tol);
    }
    let reach = (radius * 8.0).max(from.1.distance(to.1) * 4.0);
    let section: SurfaceGeometry =
        ogeom_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
    // The exact conic where the plane cuts an analytic host; a marched
    // section where the host is a patch and no closed form exists.
    let curves: Vec<Curve> = match ogeom_intersect::surface::surface_surface(&section, host, tol) {
        Ok(ogeom_intersect::surface::Meeting::Along(curves)) => curves,
        Ok(_) => ogeom_bail!(
            Construction,
            "a run-out cap's plane does not cut its host in a curve"
        ),
        Err(ogeom_core::OgeomError::NotDone(_)) => {
            match ogeom_intersect::intersect_surfaces(
                &section,
                host,
                ogeom_intersect::IntersectOptions::default(),
                tol,
            )? {
                ogeom_intersect::SurfaceIntersection::Along(sections) => {
                    sections.into_iter().map(|s| s.curve).collect()
                }
                _ => ogeom_bail!(
                    Construction,
                    "a run-out cap's plane does not cut its host in a curve"
                ),
            }
        }
        Err(e) => return Err(e),
    };
    for curve in curves {
        let pa = ogeom_algo::project_on_curve(&curve, from.1, 64, tol)?;
        let pb = ogeom_algo::project_on_curve(&curve, to.1, 64, tol)?;
        // The rail corner stands a fit error off the exact touch point; the
        // section still owns it, at the band's own slack.
        let slack = slack.max(tol.confusion() * 100.0);
        if pa.distance > slack || pb.distance > slack {
            continue;
        }
        let (pa, pb) = (pa.parameter, pb.parameter);
        if curve.is_periodic() {
            let (lo, hi) = curve.domain();
            let period = hi - lo;
            let d = (pb - pa).rem_euclid(period);
            if d <= period - d {
                return Ok(ogeom_algo::make_edge_between(
                    model,
                    curve,
                    (pa, pa + d),
                    from.0,
                    to.0,
                    tol,
                )?
                .shape);
            }
            return Ok(ogeom_algo::make_edge_between(
                model,
                curve,
                (pb, pb + (period - d)),
                to.0,
                from.0,
                tol,
            )?
            .shape
            .reversed());
        }
        if pa <= pb {
            return Ok(
                ogeom_algo::make_edge_between(model, curve, (pa, pb), from.0, to.0, tol)?.shape,
            );
        }
        return Ok(
            ogeom_algo::make_edge_between(model, curve, (pb, pa), to.0, from.0, tol)?
                .shape
                .reversed(),
        );
    }
    ogeom_bail!(
        Construction,
        "a run-out cap's section does not pass through its own corner"
    )
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
    if let Some(period) = period_v_of(host) {
        for i in 1..rail_chart.len() {
            let mut v = rail_chart[i].y;
            while v - rail_chart[i - 1].y > period / 2.0 {
                v -= period;
            }
            while rail_chart[i - 1].y - v > period / 2.0 {
                v += period;
            }
            rail_chart[i].y = v;
        }
    }
    let winding = period_of(host).map_or(0.0, |period| {
        // The unwrapped loop's chart displacement over one closing step,
        // rounded to whole periods.
        let delta = rail_chart[n - 1].x - rail_chart[0].x + (rail_chart[1].x - rail_chart[0].x);
        (delta / period).round() * period
    });
    // The loop's first station stands where the march put it, which on a
    // seam's column is either side of the seam by rounding: a station
    // re-solved onto the column came back a hair under the chart's far
    // edge, the loop unwrapped forward from there ran a whole period past
    // the window, and the host face never met the rail: its arrangement
    // drops what lies outside its chart. The whole loop is slid by whole
    // periods so its first station starts inside the window, the hair
    // under the far edge read as the near one.
    if let Some(period) = period_of(host) {
        use ogeom_geom::Surface as _;
        let ((u0, _), _) = host.domain();
        let turns = ((rail_chart[0].x - u0 + period * 1e-6) / period).floor();
        if turns != 0.0 {
            for uv in &mut rail_chart {
                uv.x -= turns * period;
            }
        }
    }
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

    // The apex ring: a fresh closed edge on the guide's own curve (turned
    // round when its own parameterization winds the chart backward, since
    // the band runs every period forward), its chart image inverted in
    // closed form and unwrapped for continuity.
    let chart_run = |guide: &Curve| -> OgeomResult<(Vec<f64>, Vec<Point2>)> {
        // Dense enough for a loop round a bore on a cone, whose image bends
        // hard where the bore leaves the wall; a fit through fewer misses.
        let samples = 256;
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
    // Exact wherever the chart has a closed form (coinciding to machine
    // precision with the strands the boolean draws from the same curve is
    // what keeps slivers out of the arrangement), fitted only past that.
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

    let apex_winds = period_of(host).is_some_and(|period| {
        let run = apex_chart[apex_chart.len() - 1].x - apex_chart[0].x;
        (run / period).round().abs() >= 1.0
    });
    if (winding.abs() > 1e-6) != apex_winds {
        // One loop winds the period and the other does not: the leg holds
        // the pole between them (a rim beside a sphere's pole, its rail
        // passing over it).
        let (wound, contractible) = if apex_winds {
            (
                (&apex, apex_pcurve, guide_range),
                (
                    rail,
                    rail_pcurve,
                    (u_params[0], u_params[u_params.len() - 1]),
                    rail_chart.as_slice(),
                ),
            )
        } else {
            (
                (
                    rail,
                    rail_pcurve,
                    (u_params[0], u_params[u_params.len() - 1]),
                ),
                (&apex, apex_pcurve, guide_range, apex_chart.as_slice()),
            )
        };
        let wound_chart = if apex_winds { &apex_chart } else { &rail_chart };
        return pole_leg(model, host, wound, wound_chart, contractible, tol);
    }
    if winding.abs() > 1e-6 {
        // Both loops wind the period; the band with its connector closes the
        // strip between them. The band wants the period run *forward*.
        if winding < 0.0 {
            ogeom_bail!(
                Construction,
                "the seat winds against its host's chart; reversing the \
                 guide is still owed; see docs/PARITY.md, fillet.edge-blends"
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

/// A leg holding a pole: the band from the loop that winds the period to
/// the pole on the other loop's side, with the other loop cut from it as a
/// hole.
fn pole_leg(
    model: &mut Model,
    host: &SurfaceGeometry,
    wound: (&Shape, PlanarCurve, (f64, f64)),
    wound_chart: &[Point2],
    contractible: (&Shape, PlanarCurve, (f64, f64), &[Point2]),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::{Curve2d as _, Surface as _};
    if !matches!(host, SurfaceGeometry::Sphere(_)) {
        ogeom_bail!(
            Construction,
            "the blend's rail winds round its host on one side and its edge \
             does not, and this host has no pole to close the leg on; a ball \
             too big for the wall beside the edge walks such a rail"
        );
    }
    let (wound_edge, wound_pcurve, wound_range) = wound;
    let (hole_edge, hole_pcurve, hole_range, hole_chart) = contractible;
    let start = wound_pcurve.point_at(wound_range.0, tol)?;
    let end = wound_pcurve.point_at(wound_range.1, tol)?;
    if end.x < start.x {
        ogeom_bail!(
            Construction,
            "the seat winds against its host's chart; reversing the \
             guide is still owed; see docs/PARITY.md, fillet.edge-blends"
        );
    }
    let mean = |chart: &[Point2]| {
        #[allow(clippy::cast_precision_loss)]
        let n = chart.len() as f64;
        chart.iter().map(|p| p.y).sum::<f64>() / n
    };
    // The pole on the hole's side of the winding loop.
    let ((u0, u1), (v0, v1)) = host.domain();
    let above = mean(hole_chart) > mean(wound_chart);
    let row = if above { v1 } else { v0 };
    let pole = ogeom_algo::make_vertex(model, host.point_at(start.x, row, tol)?).shape;
    let mut data = ogeom_topo::EdgeData::new();
    data.degenerate = true;
    let pole_edge = model.add_edge(data, &[pole.clone(), pole])?;
    let pole_pcurve: PlanarCurve = ogeom_geom::Line2d::over(
        ogeom_math::Axis2::new(
            Point2::new(start.x, row),
            ogeom_math::Direction2::new(ogeom_math::Vector2::new(1.0, 0.0), tol)?,
        ),
        0.0,
        u1 - u0,
    )?
    .into();
    let band = ogeom_algo::make_band_between(
        model,
        host,
        [(wound_edge, wound_pcurve), (&pole_edge, pole_pcurve)],
        tol,
    )?;
    let Some(surface_id) = model
        .node(&band)
        .and_then(|n| n.data().as_face())
        .map(|d| d.surface)
    else {
        ogeom_bail!(Construction, "the pole band holds no face data");
    };
    let outer = model.children_of(&band)?[0].clone();
    ogeom_algo::attach_pcurve(
        model,
        hole_edge,
        hole_pcurve,
        surface_id,
        ogeom_topo::Location::identity(),
        hole_range,
    )?;
    // The band's walk runs the winding loop forward, so it turns
    // anticlockwise in the chart when the pole is above it; the hole turns
    // the other way.
    let outer_turn = if above { 1.0 } else { -1.0 };
    let hole = if chart_area(hole_chart) * outer_turn > 0.0 {
        hole_edge.reversed()
    } else {
        hole_edge.clone()
    };
    let hole_wire = ogeom_algo::make_wire(model, &[hole], tol)?.shape;
    Ok(ogeom_algo::make_face_on(model, surface_id, &[outer, hole_wire], tol)?.shape)
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
        SurfaceGeometry::Cone(c) => {
            let (lo, hi) = pad(v);
            ogeom_geom::ConeSurface::new(c.cone(), (lo, hi))?.into()
        }
        SurfaceGeometry::Plane(p) => {
            let (ulo, uhi) = pad(u);
            let (vlo, vhi) = pad(v);
            ogeom_geom::PlaneSurface::over(p.plane(), (ulo, uhi), (vlo, vhi))?.into()
        }
        // A sphere or a torus is its whole self: neither has a window to
        // cut, and a leg on one pairs with what its rings bound.
        other => other.clone(),
    })
}

/// The chart period in `u`, for surfaces that have one.
fn period_of(surface: &SurfaceGeometry) -> Option<f64> {
    use ogeom_geom::Surface as _;
    match surface {
        SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Sphere(_)
        | SurfaceGeometry::Torus(_) => Some(core::f64::consts::TAU),
        // A patch that meets itself at its seam without being periodic (a
        // converted cylinder's) wraps like one: the march wraps its
        // parameter there, and the rail's image must unwrap it back.
        other if other.is_periodic_u() || other.is_closed_u(Tolerances::millimetres()) => {
            let ((ua, ub), _) = other.domain();
            Some(ub - ua)
        }
        _ => None,
    }
}

/// The chart period in `v`, for the one surface that has one.
fn period_v_of(surface: &SurfaceGeometry) -> Option<f64> {
    use ogeom_geom::Surface as _;
    match surface {
        SurfaceGeometry::Torus(_) => Some(core::f64::consts::TAU),
        other if other.is_periodic_v() || other.is_closed_v(Tolerances::millimetres()) => {
            let (_, (va, vb)) = other.domain();
            Some(vb - va)
        }
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
    let raw = match surface {
        SurfaceGeometry::Plane(pl) => {
            let local = pl.plane().frame().to_local(p);
            Point2::new(local.x, local.y)
        }
        SurfaceGeometry::Cylinder(c) => {
            let local = c.cylinder().frame().to_local(p);
            Point2::new(local.y.atan2(local.x), local.z)
        }
        SurfaceGeometry::Cone(c) => {
            let (u, v) = ogeom_math::elementary::cone_parameters(&c.cone(), p, tol)?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Sphere(s) => {
            let (u, v) = ogeom_math::elementary::sphere_parameters(&s.sphere(), p, tol)?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Torus(t) => {
            let (u, v) = ogeom_math::elementary::torus_parameters(&t.torus(), p, tol)?;
            Point2::new(u, v)
        }
        // No closed form: the foot by projection, warm-started from the
        // last station where there is one (consecutive stations are
        // neighbours on the surface) and seeded from a grid otherwise, or
        // where the warm start wandered off.
        _ => {
            let near = prev.and_then(|q| {
                ogeom_algo::project_on_surface_from(surface, p, (q.x, q.y), tol).ok()
            });
            let found = match near {
                Some(close) if close.distance <= tol.confusion() * 1e3 => close,
                _ => ogeom_algo::project_on_surface(surface, p, 32, tol)?,
            };
            Point2::new(found.parameters.0, found.parameters.1)
        }
    };
    let Some(prev) = prev else {
        return Ok(raw);
    };
    let unwrap = |x: f64, from: f64, period: Option<f64>| -> f64 {
        let Some(period) = period else {
            return x;
        };
        let mut x = x;
        while x - from > period / 2.0 {
            x -= period;
        }
        while from - x > period / 2.0 {
            x += period;
        }
        x
    };
    Ok(Point2::new(
        unwrap(raw.x, prev.x, period_of(surface)),
        unwrap(raw.y, prev.y, period_v_of(surface)),
    ))
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

/// The next piece of a loop walk: its curve and range, whether it runs the
/// walk's way, the vertex it arrives at, and the direction it arrives in.
type NextPiece = (Curve, (f64, f64), bool, ogeom_topo::TShapeId, Vector);

/// Whether a guide's two ends stand apart: an arc, not a loop.
fn ends_apart(guide: &Curve, range: (f64, f64), tol: Tolerances) -> bool {
    let (lo, hi) = range;
    guide
        .point_at(lo, tol)
        .and_then(|p| guide.point_at(hi, tol).map(|q| p.distance(q)))
        .is_ok_and(|d| d > tol.confusion() * 10.0)
}

/// A seat's loop closed back through the edges its two host faces share,
/// each taken where it continues the last tangentially, as one spline:
/// the pieces in their exact spline forms over unit spans, turned to run
/// the walk's way, raised to one degree and joined end to end, the seat
/// itself the first span. `None` where the walk does not come back to the
/// seat's start.
fn loop_through_neighbours(
    model: &Model,
    edge: &Shape,
    guide: &Curve,
    edge_range: (f64, f64),
    hosts: [&Shape; 2],
    tol: Tolerances,
) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    let Some((start, end)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    // The edges both hosts share, other than the seat: the rest of the
    // crease, arc by arc.
    let on_second: Vec<ogeom_topo::TShapeId> =
        explore(model, hosts[1], Filter::OfType(ShapeType::Edge))?
            .iter()
            .map(Shape::node)
            .collect();
    let mut shared: Vec<Shape> = Vec::new();
    for candidate in explore(model, hosts[0], Filter::OfType(ShapeType::Edge))? {
        if candidate.node() != edge.node()
            && on_second.contains(&candidate.node())
            && !shared.iter().any(|s| s.node() == candidate.node())
        {
            shared.push(candidate);
        }
    }
    if shared.is_empty() {
        return Ok(None);
    }
    let unit = |v: Vector| -> Option<Vector> {
        let m = v.magnitude();
        (m > tol.confusion()).then(|| v / m)
    };
    // The seat first, run forward; then each neighbour that leaves the
    // current vertex the way the last piece arrived.
    let mut pieces: Vec<(Curve, (f64, f64), bool)> = vec![(guide.clone(), edge_range, true)];
    let mut at = end.node();
    let Some(mut heading) = unit(guide.d1_at(edge_range.1, tol)?) else {
        return Ok(None);
    };
    let mut used: Vec<ogeom_topo::TShapeId> = vec![edge.node()];
    let mut closed = false;
    for _ in 0..shared.len() {
        let mut next: Option<NextPiece> = None;
        for candidate in &shared {
            if used.contains(&candidate.node()) {
                continue;
            }
            let Some((a, b)) = ogeom_algo::edge_vertices(model, candidate)? else {
                continue;
            };
            let (forward, far) = if a.node() == at {
                (true, b.node())
            } else if b.node() == at {
                (false, a.node())
            } else {
                continue;
            };
            let (curve, range) = edge_curve(model, candidate, tol)?;
            let (leaving, arriving) = if forward {
                (curve.d1_at(range.0, tol)?, curve.d1_at(range.1, tol)?)
            } else {
                (-curve.d1_at(range.1, tol)?, -curve.d1_at(range.0, tol)?)
            };
            let (Some(leaving), Some(arriving)) = (unit(leaving), unit(arriving)) else {
                continue;
            };
            if leaving.dot(heading) < 0.99 {
                continue;
            }
            used.push(candidate.node());
            next = Some((curve, range, forward, far, arriving));
            break;
        }
        let Some((curve, range, forward, far, arriving)) = next else {
            break;
        };
        pieces.push((curve, range, forward));
        heading = arriving;
        at = far;
        if at == start.node() {
            closed = true;
            break;
        }
    }
    if !closed {
        return Ok(None);
    }
    // Each piece over its own arc length rather than a unit span, so the
    // joined guide's speed is continuous across the joins: a chart image
    // fitted at the guide's parameters would otherwise carry a kink at
    // every join that the fit spends its budget on.
    let mut splines: Vec<ogeom_geom::BSplineCurve> = Vec::with_capacity(pieces.len());
    for (curve, range, forward) in &pieces {
        let mut spline = curve.to_bspline_over(*range, tol)?;
        if !forward {
            let (knots, control) =
                ogeom_math::bspline::reverse(spline.knots(), spline.control_points());
            spline = ogeom_geom::BSplineCurve::rational(knots, control)?;
        }
        let (lo, hi) = spline.domain();
        let mut length = 0.0;
        let mut last = spline.point_at(lo, tol)?;
        for k in 1..=64 {
            let p = spline.point_at(lo + (hi - lo) * f64::from(k) / 64.0, tol)?;
            length += last.distance(p);
            last = p;
        }
        if length > tol.confusion() {
            spline = ogeom_geom::BSplineCurve::rational(
                spline.knots().reparameterized(0.0, length)?,
                spline.control_points().to_vec(),
            )?;
        }
        splines.push(spline);
    }
    let degree = splines
        .iter()
        .map(ogeom_geom::BSplineCurve::degree)
        .max()
        .unwrap_or(1);
    for spline in &mut splines {
        while spline.degree() < degree {
            *spline = spline.elevated(tol)?;
        }
    }
    let seat = splines[0].domain();
    let mut whole = (
        splines[0].knots().clone(),
        splines[0].control_points().to_vec(),
    );
    for spline in &splines[1..] {
        whole = ogeom_math::bspline::join(
            &whole,
            &(spline.knots().clone(), spline.control_points().to_vec()),
        )?;
    }
    Ok(Some((
        Curve::BSpline(ogeom_geom::BSplineCurve::rational(whole.0, whole.1)?),
        seat,
    )))
}
