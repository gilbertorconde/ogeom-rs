//! Remodelling: geometry restated wholesale, the solid unchanged. Swept
//! forms are made by restating primitives as revolutions, then named back
//! as what they are. High-degree splines are made by restating a box in
//! elevated form, then restricted back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::{OgeomResult, Tolerances};
use ogeom::geom::{
    BSplineCurve, BSplineSurface, CircleCurve, Curve, Curve2d as _, Curve3d as _, LineCurve,
    RevolutionSurface, Surface as _, SurfaceGeometry, TrimmedCurve,
};
use ogeom::math::bspline::ControlGrid;
use ogeom::math::{Axis, Circle, Frame, KnotVector, Point, Weighted};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn kinds(model: &Model, shape: &Shape) -> Vec<&'static str> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                panic!("not a face");
            };
            match model.geometry().surface(d.surface).unwrap() {
                SurfaceGeometry::Plane(_) => "plane",
                SurfaceGeometry::Cylinder(_) => "cylinder",
                SurfaceGeometry::Cone(_) => "cone",
                SurfaceGeometry::Sphere(_) => "sphere",
                SurfaceGeometry::Torus(_) => "torus",
                SurfaceGeometry::Revolution(_) => "revolution",
                SurfaceGeometry::BSpline(_) => "spline",
                _ => "other",
            }
        })
        .collect()
}

fn holds(model: &Model, shape: &Shape, volume: f64) {
    let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let fine = Deflection {
        chord: 1e-3,
        angular: 0.02,
        ..Deflection::default()
    };
    let measured = ogeom::algo::volume_properties(model, shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (measured - volume).abs() <= volume * 1e-3,
        "{measured} against {volume}"
    );
}

/// The elementary surfaces restated as the revolutions they are.
fn as_revolution(s: &SurfaceGeometry) -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
    let turn = core::f64::consts::TAU;
    let axis_of = |f: Frame| Axis {
        location: f.origin(),
        direction: f.z(),
    };
    let revolved = |profile: Curve, frame: Frame| -> OgeomResult<SurfaceGeometry> {
        Ok(RevolutionSurface::new(profile, axis_of(frame), turn)?.into())
    };
    Ok(Some(match s {
        SurfaceGeometry::Cylinder(c) => {
            let (lo, hi) = c.domain().1;
            let f = c.cylinder().frame();
            let ruling = Axis {
                location: f.origin() + f.x().vector() * c.cylinder().radius(),
                direction: f.z(),
            };
            (revolved(LineCurve::over(ruling, lo, hi)?.into(), f)?, false)
        }
        SurfaceGeometry::Sphere(b) => {
            let f = b.sphere().frame();
            let meridian = Frame::new(f.origin(), -f.y(), f.x(), T)?;
            let circle: Curve =
                CircleCurve::new(Circle::new(meridian, b.sphere().radius(), T)?).into();
            let half = core::f64::consts::FRAC_PI_2;
            (
                revolved(TrimmedCurve::new(circle, -half, half, T)?.into(), f)?,
                false,
            )
        }
        SurfaceGeometry::Torus(t) => {
            let f = t.torus().frame();
            let centre = f.origin() + f.x().vector() * t.torus().major_radius();
            let tube = Frame::new(centre, -f.y(), f.x(), T)?;
            let circle: Curve =
                CircleCurve::new(Circle::new(tube, t.torus().minor_radius(), T)?).into();
            (revolved(circle, f)?, false)
        }
        _ => return Ok(None),
    }))
}

fn keep_curves(_: &Curve, _: (f64, f64)) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    Ok(None)
}

fn round_trip(solid: &Shape, model: &mut Model, volume: f64, expect: &str) {
    let swept = ogeom::algo::restate_geometry(model, solid, &as_revolution, &keep_curves, T)
        .unwrap()
        .shape;
    assert!(kinds(model, &swept).contains(&"revolution"));
    holds(model, &swept, volume);
    let named = ogeom::heal::swept_to_elementary(model, &swept, T)
        .unwrap()
        .shape;
    let found = kinds(model, &named);
    assert!(!found.contains(&"revolution"), "{found:?}");
    assert!(found.contains(&expect), "{found:?}");
    holds(model, &named, volume);
}

#[test]
fn a_revolved_ruling_is_a_drum() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    round_trip(
        &drum,
        &mut model,
        core::f64::consts::PI * 4.0 * 5.0,
        "cylinder",
    );
}

#[test]
fn a_revolved_meridian_is_a_ball() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
        .unwrap()
        .shape;
    round_trip(
        &ball,
        &mut model,
        4.0 / 3.0 * core::f64::consts::PI * 27.0,
        "sphere",
    );
}

#[test]
fn a_revolved_tube_circle_is_a_torus() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 5.0, 1.0, T)
        .unwrap()
        .shape;
    round_trip(
        &ring,
        &mut model,
        2.0 * core::f64::consts::PI.powi(2) * 5.0,
        "torus",
    );
}

/// A plane as a cubic patch over the same chart window, and a line as a
/// cubic over the same range: the same points at the same parameters.
fn elevated_plane(s: &SurfaceGeometry) -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
    let SurfaceGeometry::Plane(_) = s else {
        return Ok(None);
    };
    let ((u0, u1), (v0, v1)) = s.domain();
    let mut points = Vec::with_capacity(16);
    for i in 0..4 {
        for j in 0..4 {
            let (fu, fv) = (f64::from(i) / 3.0, f64::from(j) / 3.0);
            points.push(s.point_at(u0 + (u1 - u0) * fu, v0 + (v1 - v0) * fv, T)?);
        }
    }
    let knots = |a: f64, b: f64| KnotVector::new(vec![a, a, a, a, b, b, b, b], 3);
    let patch = BSplineSurface::new(
        knots(u0, u1)?,
        knots(v0, v1)?,
        &ControlGrid::new(points, 4, 4)?,
        T,
    )?;
    Ok(Some((patch.into(), false)))
}

fn elevated_line(c: &Curve, range: (f64, f64)) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    let Curve::Line(_) = c else {
        return Ok(None);
    };
    let control = (0..4)
        .map(|k| {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 3.0;
            Weighted::new(c.point_at(t, T)?, 1.0, T)
        })
        .collect::<OgeomResult<Vec<Weighted<Point>>>>()?;
    let (a, b) = range;
    let knots = KnotVector::new(vec![a, a, a, a, b, b, b, b], 3)?;
    Ok(Some((
        BSplineCurve::rational(knots, control)?.into(),
        range,
    )))
}

#[test]
fn a_cubic_box_restricts_to_degree_one() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let cubic =
        ogeom::algo::restate_geometry(&mut model, &block, &elevated_plane, &elevated_line, T)
            .unwrap()
            .shape;
    assert!(kinds(&model, &cubic).iter().all(|k| *k == "spline"));
    holds(&model, &cubic, 1000.0);

    let linear = ogeom::heal::restrict_degree(&mut model, &cubic, 1, 1e-6, T)
        .unwrap()
        .shape;
    for face in explore_unique(&model, &linear, ShapeType::Face).unwrap() {
        let Some(NodeData::Face(d)) = model.node(&face).map(|n| n.data()) else {
            panic!("not a face");
        };
        let SurfaceGeometry::BSpline(s) = model.geometry().surface(d.surface).unwrap() else {
            panic!("still a spline");
        };
        assert_eq!((s.u_knots().degree(), s.v_knots().degree()), (1, 1));
    }
    for edge in explore_unique(&model, &linear, ShapeType::Edge).unwrap() {
        let Some(NodeData::Edge(d)) = model.node(&edge).map(|n| n.data()) else {
            panic!("not an edge");
        };
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, .. }) = d.curve3d() else {
            continue;
        };
        if let Curve::BSpline(b) = model.geometry().curve(*curve).unwrap() {
            assert_eq!(b.degree(), 1);
        }
    }
    holds(&model, &linear, 1000.0);
}

#[test]
fn a_torus_rebuilds_with_both_seams() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 5.0, 1.0, T)
        .unwrap()
        .shape;
    let keep = |_: &SurfaceGeometry| -> OgeomResult<Option<(SurfaceGeometry, bool)>> { Ok(None) };
    let same = ogeom::algo::restate_geometry(&mut model, &ring, &keep, &keep_curves, T)
        .unwrap()
        .shape;
    holds(&model, &same, 2.0 * core::f64::consts::PI.powi(2) * 5.0);
    let nurbs = ogeom::algo::to_nurbs(&mut model, &ring, T).unwrap().shape;
    holds(&model, &nurbs, 2.0 * core::f64::consts::PI.powi(2) * 5.0);
}

/// Planes and drums restated as offsets of a plane or drum set back by
/// half a unit: the same points, spelt with no exact spline form.
fn as_offset(s: &SurfaceGeometry) -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
    let back = 0.5;
    let basis: SurfaceGeometry = match s {
        SurfaceGeometry::Plane(p) => {
            let f = p.plane().frame();
            let moved = Frame::new(f.origin() - f.z().vector() * back, f.z(), f.x(), T)?;
            let (u, v) = s.domain();
            ogeom::geom::PlaneSurface::over(ogeom::math::Plane::new(moved), u, v)?.into()
        }
        SurfaceGeometry::Cylinder(c) => {
            let cylinder = c.cylinder();
            ogeom::geom::CylinderSurface::new(
                ogeom::math::Cylinder::new(cylinder.frame(), cylinder.radius() - back, T)?,
                c.domain().1,
            )?
            .into()
        }
        _ => return Ok(None),
    };
    Ok(Some((
        SurfaceGeometry::Offset(Box::new(ogeom::geom::OffsetSurface::new(basis, back)?)),
        false,
    )))
}

#[test]
fn offset_faces_convert_to_nurbs_within_a_tolerance() {
    for drum in [false, true] {
        let mut model = Model::new();
        let (solid, volume) = if drum {
            (
                ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
                    .unwrap()
                    .shape,
                core::f64::consts::PI * 4.0 * 5.0,
            )
        } else {
            (
                ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
                    .unwrap()
                    .shape,
                1000.0,
            )
        };
        let offset = ogeom::algo::restate_geometry(&mut model, &solid, &as_offset, &keep_curves, T)
            .unwrap()
            .shape;
        assert!(kinds(&model, &offset).contains(&"other"));
        assert!(ogeom::algo::to_nurbs(&mut model, &offset, T).is_err());
        let nurbs = ogeom::algo::to_nurbs_within(&mut model, &offset, 1e-4, T)
            .unwrap()
            .shape;
        assert!(kinds(&model, &nurbs).iter().all(|k| *k == "spline"));
        holds(&model, &nurbs, volume);
    }
}

/// The distance from `p` to the nearest point of `curve` over `range`: a
/// scan, then a golden-section search about the nearest sample.
fn nearest_on(curve: &Curve, range: (f64, f64), p: ogeom::math::Point) -> f64 {
    let at = |t: f64| curve.point_at(t, T).unwrap().distance(p);
    let step = (range.1 - range.0) / 400.0;
    let mut best = (range.0, at(range.0));
    for k in 1..=400 {
        let t = range.0 + step * f64::from(k);
        if at(t) < best.1 {
            best = (t, at(t));
        }
    }
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let (mut a, mut b) = ((best.0 - step.abs()).max(lo), (best.0 + step.abs()).min(hi));
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..100 {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if at(c) < at(d) {
            b = d;
        } else {
            a = c;
        }
    }
    best.1.min(at(f64::midpoint(a, b)))
}

/// The edges whose pcurves, lifted through their surfaces, leave their
/// curve by more than the edge's tolerance anywhere along it: how far, and
/// the tolerance. Each pcurve is sampled along its range, and each lifted
/// point measured against the curve's point at the same fraction of its
/// range or, where that is farther than the tolerance, the nearest point
/// of the curve. A pcurve at another placement than the curve describes
/// another occurrence and is not compared.
fn pcurves_off_their_curves(model: &Model, shape: &Shape) -> Vec<(f64, f64)> {
    let mut off = Vec::new();
    for edge in explore_unique(model, shape, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d {
            curve,
            range,
            location,
        }) = data.curve3d()
        else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let tolerance = data.tolerance.get();
        for repr in &data.representations {
            let (pcurves, pcurve_range, surface) = match repr {
                EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    location: at,
                } if at == location => (vec![*curve], *range, *surface),
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    surface,
                    location: at,
                } if at == location => (vec![*forward, *reversed], *range, *surface),
                _ => continue,
            };
            let surface = model.geometry().surface(surface).unwrap();
            for id in pcurves {
                let pcurve = model.geometry().pcurve(id).unwrap();
                let mut widest: f64 = 0.0;
                for k in 0..=64 {
                    let f = f64::from(k) / 64.0;
                    let uv = pcurve
                        .point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * f, T)
                        .unwrap();
                    let Ok(lifted) = surface.point_at(uv.x, uv.y, T) else {
                        continue;
                    };
                    let mut gap = curve
                        .point_at(range.0 + (range.1 - range.0) * f, T)
                        .unwrap()
                        .distance(lifted);
                    if gap > tolerance {
                        gap = gap.min(nearest_on(curve, *range, lifted));
                    }
                    widest = widest.max(gap);
                }
                if widest > tolerance {
                    off.push((widest, tolerance));
                }
            }
        }
    }
    off
}

/// Faces restated on offsets of a drum and of planes keep every pcurve on
/// its edge to the edge's own tolerance, and so does the same solid fitted
/// as splines, sampled along each pcurve independently of the checker.
/// The offsets are exact, so their edges keep what they stated; the fitted
/// drum misses the true one between its samples, and the rims say so.
#[test]
fn offset_faces_keep_their_pcurves_on_their_edges() {
    for drum in [false, true] {
        let mut model = Model::new();
        let solid = if drum {
            ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
                .unwrap()
                .shape
        } else {
            ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
                .unwrap()
                .shape
        };
        let offset = ogeom::algo::restate_geometry(&mut model, &solid, &as_offset, &keep_curves, T)
            .unwrap()
            .shape;
        assert_eq!(pcurves_off_their_curves(&model, &offset), []);
        for edge in explore_unique(&model, &offset, ShapeType::Edge).unwrap() {
            let Some(NodeData::Edge(d)) = model.node(&edge).map(|n| n.data()) else {
                panic!("not an edge");
            };
            assert!(d.tolerance.get() <= 1e-5, "{}", d.tolerance.get());
        }
        let nurbs = ogeom::algo::to_nurbs_within(&mut model, &offset, 1e-4, T)
            .unwrap()
            .shape;
        assert_eq!(pcurves_off_their_curves(&model, &nurbs), []);
    }
}
