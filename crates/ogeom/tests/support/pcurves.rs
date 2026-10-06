//! Pcurves compared across an exchange round trip: each pcurve a read
//! shape holds against the source's pcurve for the same edge on the same
//! surface.

use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom::math::{Point, Point2};
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

/// Parameters agree to this: the file states every coordinate in its
/// shortest round-trip form, so a pcurve read as written evaluates to the
/// same bits, and a period subtracted leaves rounding far below this. A
/// pcurve derived on reading instead stands a fit's error away.
const SAME: f64 = 1e-9;

/// How many of the read shape's pcurves came back as the source's, and
/// how many it holds.
pub fn kept(source: (&Model, &Shape), read: (&Model, &Shape), tol: Tolerances) -> (usize, usize) {
    let (mut same, mut total) = (0, 0);
    let sources: Vec<(Shape, Point)> = explore_unique(source.0, source.1, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|e| Some((e.clone(), midpoint(source.0, &e, tol)?.0)))
        .collect();
    for face in explore_unique(read.0, read.1, ShapeType::Face).unwrap() {
        let held = read
            .0
            .node(&face)
            .unwrap()
            .data()
            .as_face()
            .unwrap()
            .surface;
        let surface = read.0.geometry().surface(held).unwrap();
        for edge in explore_unique(read.0, &face, ShapeType::Edge).unwrap() {
            let data = read.0.node(&edge).unwrap().data().as_edge().unwrap();
            if data.degenerate {
                continue;
            }
            let Some((pcurves, range)) = pcurves_of(read.0, data.pcurve_for(held, edge.location()))
            else {
                continue;
            };
            total += 1;
            let Some((at, read_curve)) = midpoint(read.0, &edge, tol) else {
                continue;
            };
            let Some((twin, _)) = sources
                .iter()
                .find(|(_, p)| p.distance(at) <= tol.confusion())
            else {
                continue;
            };
            let twin_data = source.0.node(twin).unwrap().data().as_edge().unwrap();
            let (_, twin_curve) = midpoint(source.0, twin, tol).unwrap();
            let found = twin_data.parametric_surfaces().into_iter().any(|sid| {
                let theirs = source.0.geometry().surface(sid).unwrap();
                let Some((originals, held_over)) = pcurves_of(source.0, twin_data.pcurve_on(sid))
                else {
                    return false;
                };
                same_chart(surface, theirs, &pcurves[0], range, tol)
                    && pcurves.iter().all(|p| {
                        originals.iter().any(|o| {
                            apart(
                                surface,
                                (p, range, read_curve),
                                (o, held_over, twin_curve),
                                tol,
                            ) <= SAME
                        })
                    })
            });
            if found {
                same += 1;
            }
        }
    }
    (same, total)
}

/// The pcurves a representation names, one or a seam's two, and its range.
fn pcurves_of(model: &Model, repr: Option<&EdgeRepr>) -> Option<(Vec<PlanarCurve>, (f64, f64))> {
    let pcurve = |id| model.geometry().pcurve(id).cloned();
    match repr? {
        EdgeRepr::PCurve { curve, range, .. } => Some((vec![pcurve(*curve)?], *range)),
        EdgeRepr::Seam {
            forward,
            reversed,
            range,
            ..
        } => Some((vec![pcurve(*forward)?, pcurve(*reversed)?], *range)),
        _ => None,
    }
}

/// The point halfway along an edge's curve, and the edge's range on its
/// curve with the curve's period, zero where it has none.
fn midpoint(model: &Model, edge: &Shape, tol: Tolerances) -> Option<(Point, Span)> {
    let data = model.node(edge)?.data().as_edge()?;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return None;
    };
    let curve = model.geometry().curve(*curve)?;
    let period = if curve.is_periodic() {
        let (lo, hi) = curve.domain();
        hi - lo
    } else {
        0.0
    };
    let at = curve.point_at(0.5 * (range.0 + range.1), tol).ok()?;
    Some((at, (*range, period)))
}

/// An edge's range on its curve, and the curve's period or zero.
type Span = ((f64, f64), f64);

/// Whether two surfaces put the same point at a pcurve's parameters: the
/// same surface, the same chart.
fn same_chart(
    ours: &SurfaceGeometry,
    theirs: &SurfaceGeometry,
    pcurve: &PlanarCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> bool {
    (0..=4).all(|i| {
        let t = range.0 + (range.1 - range.0) * f64::from(i) / 4.0;
        let Ok(uv) = pcurve.point_at(t, tol) else {
            return false;
        };
        match (
            ours.point_at(uv.x, uv.y, tol),
            theirs.point_at(uv.x, uv.y, tol),
        ) {
            (Ok(a), Ok(b)) => a.distance(b) <= tol.confusion(),
            _ => false,
        }
    })
}

/// How far two pcurves of one edge stand apart in the chart, each over
/// its own range and compared at the same parameter of the edge's curve,
/// which each range follows in step with the curve's: whole periods of a
/// closed curve, and of a closed direction of the chart, aside.
fn apart(
    surface: &SurfaceGeometry,
    (a, over_a, (curve_a, period)): (&PlanarCurve, (f64, f64), Span),
    (b, over_b, (curve_b, _)): (&PlanarCurve, (f64, f64), Span),
    tol: Tolerances,
) -> f64 {
    let turns = if period > 0.0 {
        period * ((curve_a.0 - curve_b.0) / period).round()
    } else {
        0.0
    };
    let paced = |t: f64, (lo, hi): (f64, f64), curve: (f64, f64)| {
        lo + (t - curve.0) * (hi - lo) / (curve.1 - curve.0)
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let reduce = |d: f64, periodic: bool, width: f64| {
        if periodic {
            d - width * (d / width).round()
        } else {
            d
        }
    };
    (0..=16)
        .map(|i| {
            let t = (curve_a.1 - curve_a.0).mul_add(f64::from(i) / 16.0, curve_a.0);
            let (s, r) = (paced(t, over_a, curve_a), paced(t - turns, over_b, curve_b));
            match (a.point_at(s, tol), b.point_at(r, tol)) {
                (Ok(p), Ok(q)) => {
                    let du = reduce(p.x - q.x, surface.is_periodic_u(), u1 - u0);
                    let dv = reduce(p.y - q.y, surface.is_periodic_v(), v1 - v0);
                    Point2::new(du, dv).distance(Point2::new(0.0, 0.0))
                }
                _ => f64::INFINITY,
            }
        })
        .fold(0.0, f64::max)
}
