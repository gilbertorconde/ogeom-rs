//! Normal projection: a wire dropped onto a shape along its faces' normals.
//!
//! The projection of a point onto a surface is its foot — the point where
//! the displacement is perpendicular to both tangents — and the projection
//! of a curve is the curve through the feet. That curve almost never has a
//! closed form, so it is sampled, fitted to a stated tolerance, and carries
//! the pcurve fitted *with* it: same parameter in space and in the chart,
//! which is what makes the result an edge a face can be split along.
//!
//! Which face a point lands on is decided by measurement, not by order:
//! every face is asked, the nearest foot inside a face's own trim wins, and
//! a sample no face claims ends the run it was in. So a wire projected onto
//! a solid comes back as one edge per stretch that actually landed, and the
//! stretches that fell off the shape are simply absent.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry, Transformable as _,
};
use ogeom_math::{Point, Point2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, SurfaceId, explore_unique};

/// One projected stretch: the edge that was built, and the face it lies on.
#[derive(Debug, Clone)]
pub struct Projected {
    /// The edge, with its pcurve attached on the face's surface.
    pub edge: Shape,
    /// The face it landed on.
    pub face: Shape,
    /// How far the fitted curve may sit from the sampled feet.
    pub tolerance: f64,
}

/// Project every edge of `wire` onto the faces of `target`.
///
/// `stations` is how finely each edge is sampled; the fit is held to
/// `tolerance` against those samples. Both are the caller's, because both
/// are the answer's accuracy and this cannot guess what it is for.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `stations` is under four or `tolerance` is not a usable distance, and
/// whatever the fit refuses.
pub fn normal_projection(
    model: &mut Model,
    target: &Shape,
    wire: &Shape,
    stations: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<(Vec<Projected>, Built)> {
    if stations < 4 {
        ogeom_bail!(
            Construction,
            "a projection sampled at {stations} stations is a guess"
        );
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }

    // The faces to land on, with their surfaces in world space and their
    // trims as chart rings — a foot outside the trim is on the surface but
    // not on the face, and landing there would be a projection onto
    // geometry the shape does not have.
    let mut seats: Vec<Seat> = Vec::new();
    let deflection = ogeom_mesh::Deflection::default();
    for face in explore_unique(model, target, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
            ogeom_bail!(
                Construction,
                "a scaled placement changes a surface's parameterization out \
                 from under its pcurves; bake the scale before projecting"
            );
        }
        let rings = ogeom_mesh::face_boundary(model, &face, deflection, tol)?;
        seats.push(Seat {
            face,
            surface_id: data.surface,
            surface: surface.transformed(&placement, tol)?,
            rings,
        });
    }
    if seats.is_empty() {
        ogeom_bail!(Construction, "a shape with no faces catches nothing");
    }

    let mut out = Vec::new();
    let mut history = History::new();
    for edge in explore_unique(model, wire, ShapeType::Edge)? {
        let (curve, range) = {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve).cloned() else {
                continue;
            };
            (geometry, *range)
        };
        // Walk the edge, landing each station on the nearest face that will
        // have it, and break the run wherever the face changes or nothing
        // catches — each run becomes one projected edge.
        let mut run: Vec<(usize, Point, Point2)> = Vec::new();
        let mut runs: Vec<(usize, Vec<(Point, Point2)>)> = Vec::new();
        let mut flush = |run: &mut Vec<(usize, Point, Point2)>| {
            if run.len() >= 4 {
                let seat = run[0].0;
                runs.push((seat, run.iter().map(|(_, p, uv)| (*p, *uv)).collect()));
            }
            run.clear();
        };
        for k in 0..=stations {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a station index, far below the mantissa"
            )]
            let t = (range.1 - range.0).mul_add(k as f64 / stations as f64, range.0);
            let at = curve.point_at(t, tol)?;
            let landed = nearest_seat(&seats, at, tol)?;
            match landed {
                Some((seat, point, uv)) => {
                    if run.first().is_some_and(|(held, _, _)| *held != seat) {
                        flush(&mut run);
                    }
                    run.push((seat, point, uv));
                }
                None => flush(&mut run),
            }
        }
        flush(&mut run);

        for (seat, samples) in runs {
            let points: Vec<Point> = samples.iter().map(|(p, _)| *p).collect();
            let mut chart: Vec<Point2> = samples.iter().map(|(_, uv)| *uv).collect();
            // A periodic chart's parameters come back folded into the
            // surface's own window, so a run crossing the seam arrives torn
            // — and a fit through a tear is a fit through a jump it cannot
            // make. Unwrapped, the run is continuous again.
            unwrap(&mut chart, &seats[seat].surface);
            // Fitted together, so the two descriptions share a parameter:
            // the pcurve rides the same knots as the curve.
            let (fitted, on_face, _) =
                ogeom_geom::fit::fit_points_joint(&points, &chart, &chart, 3, tolerance, tol)?;
            let curve: Curve = fitted.curve.into();
            let built = ogeom_algo::make_edge(model, curve, (0.0, 1.0), tol)?.shape;
            let pcurve: PlanarCurve = on_face.into();
            ogeom_algo::attach_pcurve(
                model,
                &built,
                pcurve,
                seats[seat].surface_id,
                ogeom_topo::Location::identity(),
                (0.0, 1.0),
            )?;
            history.generate(&edge, built.clone());
            out.push(Projected {
                edge: built,
                face: seats[seat].face.clone(),
                tolerance: fitted.error,
            });
        }
    }

    let edges: Vec<Shape> = out.iter().map(|p| p.edge.clone()).collect();
    let result = model.add_compound(&edges)?;
    history.modify(wire, result.clone());
    Ok((out, Built::new(result, history)))
}

/// A face ready to catch a projection.
struct Seat {
    face: Shape,
    surface_id: SurfaceId,
    /// The surface in world space.
    surface: SurfaceGeometry,
    /// The face's trim, as chart rings.
    rings: Vec<Vec<Point2>>,
}

/// The nearest face whose *trim* holds the foot, with the foot and its
/// chart position.
fn nearest_seat(
    seats: &[Seat],
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Option<(usize, Point, Point2)>> {
    let mut best: Option<(usize, Point, Point2, f64)> = None;
    for (i, seat) in seats.iter().enumerate() {
        let projection = ogeom_algo::project_on_surface(&seat.surface, at, 24, tol)?;
        let (u, v) = projection.parameters;
        let uv = Point2::new(u, v);
        if !inside_rings(&seat.rings, uv) {
            continue;
        }
        let foot = seat.surface.point_at(u, v, tol)?;
        let distance = foot.distance(at);
        if best.as_ref().is_none_or(|(.., held)| distance < *held) {
            best = Some((i, foot, uv, distance));
        }
    }
    Ok(best.map(|(i, foot, uv, _)| (i, foot, uv)))
}

/// Undo the folding a periodic chart applies: every step longer than half a
/// period is that period the other way.
fn unwrap(chart: &mut [Point2], surface: &SurfaceGeometry) {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let periods = [
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    ];
    for k in 1..chart.len() {
        let previous = chart[k - 1];
        let mut here = chart[k];
        for (axis, period) in periods.iter().enumerate() {
            if *period <= 0.0 {
                continue;
            }
            let (was, is) = if axis == 0 {
                (previous.x, here.x)
            } else {
                (previous.y, here.y)
            };
            let shifted = (is - was) / period;
            let turns = shifted.round();
            if turns.abs() >= 1.0 {
                if axis == 0 {
                    here.x = turns.mul_add(-period, is);
                } else {
                    here.y = turns.mul_add(-period, is);
                }
            }
        }
        chart[k] = here;
    }
}

/// Even-odd containment against a face's chart rings, holes included by the
/// same counting.
fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = (b.x - a.x).mul_add((p.y - a.y) / (b.y - a.y), a.x);
                if x > p.x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}
