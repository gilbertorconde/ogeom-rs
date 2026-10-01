//! The smallest axis-aligned box that holds a shape: how big a body is.
//!
//! [`shape_bounds`](crate::shape_bounds) promises to hold everything and
//! so keeps each surface's carrier bound, which for a surface of
//! revolution or a trimmed cone can stand well clear of the solid. Here
//! each extreme is found where it is: on an edge, by a one-dimensional
//! search along its curve, or inside a face, at a point of the exact
//! surface found from the face's own mesh and held to the face's chart.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _};
use ogeom_math::{Aabb, Point, Point2, Vector, solve};
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const AXES: [Vector; 3] = [Vector::X, Vector::Y, Vector::Z];

/// The smallest axis-aligned box holding `shape`, to within `tol`: each of
/// its six sides where the shape actually reaches, on an edge, at a vertex
/// or inside a face.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// does not resolve; [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if the shape holds nothing to bound, or as the face meshes report.
pub fn tight_bounds(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Aabb> {
    let mut points: Vec<Point> = Vec::new();

    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        points.push(vertex.transform(model.datums())?.apply(data.point));
    }

    // Along every edge, each coordinate's least and greatest, the search
    // bracketed by the samples that led.
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placement = edge.transform(model.datums())?;
        let at = |t: f64| -> OgeomResult<Point> { Ok(placement.apply(geometry.point_at(t, tol)?)) };
        const N: usize = 64;
        #[allow(clippy::cast_precision_loss)]
        let step = (range.1 - range.0) / N as f64;
        let samples: Vec<(f64, Point)> = (0..=N)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = range.0 + step * i as f64;
                at(t).map(|p| (t, p))
            })
            .collect::<OgeomResult<_>>()?;
        for axis in AXES {
            for sense in [1.0, -1.0] {
                let score = |p: Point| p.to_vector().dot(axis) * sense;
                let Some((best, _)) = samples
                    .iter()
                    .enumerate()
                    .max_by(|a, b| score(a.1.1).total_cmp(&score(b.1.1)))
                else {
                    continue;
                };
                let lo = samples[best.saturating_sub(1)].0;
                let hi = samples[(best + 1).min(N)].0;
                points.push(samples[best].1);
                if hi > lo {
                    let found = solve::minimize(
                        |t| at(t).map_or(f64::INFINITY, |p| -score(p)),
                        lo,
                        hi,
                        solve::Criteria::default(),
                    )?;
                    points.push(at(found.value)?);
                }
            }
        }
    }

    // Inside every face: the mesh's leading vertex for each side, moved on
    // the exact surface as far as the side leads while it stays inside the
    // face's own chart triangles.
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "surface is not in this model");
        };
        let placement = face.transform(model.datums())?;
        let mesh = ogeom_mesh::triangulate_face(model, &face, Deflection::default(), tol)?;
        if mesh.positions.is_empty() {
            continue;
        }
        let chart: Vec<[Point2; 3]> = mesh
            .triangles
            .iter()
            .map(|t| {
                t.map(|i| {
                    let (u, v) = mesh.parameters[i as usize];
                    Point2::new(u, v)
                })
            })
            .collect();
        let inside = |p: Point2| chart.iter().any(|t| in_triangle(*t, p));
        let span = chart.iter().flatten().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |b, p| (b.0.min(p.x), b.1.max(p.x), b.2.min(p.y), b.3.max(p.y)),
        );
        let lower = [span.0, span.2];
        let upper = [span.1, span.3];
        if !(upper[0] > lower[0] && upper[1] > lower[1]) {
            continue;
        }
        // The chart triangles' own border: the sides only one triangle has.
        // It runs on chords of the face's boundary, and where that boundary
        // curves in the chart the chords stand outside it. A descent that
        // ends against the border was stopped by the chords, not by the
        // surface, and an extreme on the boundary is the edges' to find.
        let mut sides: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        for t in &mesh.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *sides.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let border: Vec<(Point2, Point2)> = sides
            .into_iter()
            .filter(|(_, n)| *n == 1)
            .map(|((a, b), _)| {
                let at = |i: u32| {
                    let (u, v) = mesh.parameters[i as usize];
                    Point2::new(u, v)
                };
                (at(a), at(b))
            })
            .collect();
        let against = ((upper[0] - lower[0]).min(upper[1] - lower[1])) * 1e-6;
        let on_border = |p: Point2| {
            border
                .iter()
                .any(|(a, b)| distance_to_segment(p, *a, *b) <= against)
        };
        // A sphere's and a torus's extreme along each side is closed form,
        // and where it falls inside the face it is the face's extreme; the
        // descent below can stall short of it in a general frame. Angles
        // are tried a whole turn either way, as the chart may run them.
        use ogeom_geom::Transformable as _;
        let placed = surface.transformed(&placement, tol)?;
        let tau = core::f64::consts::TAU;
        let inside_turned = |u: f64, v: f64, periodic_v: bool| {
            let turns: &[f64] = &[0.0, 1.0, -1.0];
            turns.iter().any(|&a| {
                let shifts: &[f64] = if periodic_v { turns } else { &[0.0] };
                shifts
                    .iter()
                    .any(|&b| inside(Point2::new(u + a * tau, v + b * tau)))
            })
        };
        for axis in AXES {
            for sense in [1.0, -1.0] {
                let toward = axis * sense;
                let extreme = match &placed {
                    ogeom_geom::SurfaceGeometry::Sphere(s) => {
                        let sphere = s.sphere();
                        let p = sphere.centre() + toward * sphere.radius();
                        ogeom_math::elementary::sphere_parameters(&sphere, p, tol)
                            .ok()
                            .map(|(u, v)| (p, u, v, false))
                    }
                    ogeom_geom::SurfaceGeometry::Torus(t) => {
                        let torus = t.torus();
                        let frame = torus.frame();
                        let z = frame.z().vector();
                        let along = toward.dot(z);
                        let across = toward - z * along;
                        let reach = across.magnitude();
                        // Facing straight along the axis every parallel
                        // leads equally; any one of them will do.
                        let out = if reach > 1e-12 {
                            across / reach
                        } else {
                            frame.x().vector()
                        };
                        let lean = along.hypot(reach);
                        let tilt = (out * reach + z * along) / lean;
                        let p = frame.origin()
                            + out * torus.major_radius()
                            + tilt * torus.minor_radius();
                        ogeom_math::elementary::torus_parameters(&torus, p, tol)
                            .ok()
                            .map(|(u, v)| (p, u, v, true))
                    }
                    _ => None,
                };
                if let Some((p, u, v, periodic_v)) = extreme
                    && inside_turned(u, v, periodic_v)
                {
                    points.push(p);
                }
            }
        }
        for axis in AXES {
            for sense in [1.0, -1.0] {
                let score = |p: Point| p.to_vector().dot(axis) * sense;
                let Some(best) = (0..mesh.positions.len())
                    .max_by(|a, b| score(mesh.positions[*a]).total_cmp(&score(mesh.positions[*b])))
                else {
                    continue;
                };
                points.push(mesh.positions[best]);
                let (u, v) = mesh.parameters[best];
                let refined = solve_on_face(
                    |uv: &[f64]| {
                        let p = Point2::new(uv[0], uv[1]);
                        if !inside(p) {
                            return f64::INFINITY;
                        }
                        surface
                            .point_at(uv[0], uv[1], tol)
                            .map_or(f64::INFINITY, |q| -score(placement.apply(q)))
                    },
                    [u, v],
                    lower,
                    upper,
                )?;
                if let Some([u, v]) = refined
                    && !on_border(Point2::new(u, v))
                {
                    points.push(placement.apply(surface.point_at(u, v, tol)?));
                }
            }
        }
    }

    if points.is_empty() {
        ogeom_bail!(Construction, "the shape holds nothing to bound");
    }
    Ok(Aabb::of_points(&points))
}

/// A local descent from `start` inside the chart window; `None` where it
/// found nothing better than the start.
fn solve_on_face(
    f: impl FnMut(&[f64]) -> f64,
    start: [f64; 2],
    lower: [f64; 2],
    upper: [f64; 2],
) -> OgeomResult<Option<[f64; 2]>> {
    let step = ((upper[0] - lower[0]).min(upper[1] - lower[1]) * 1e-2).max(1e-9);
    let found = ogeom_math::minimize_local(f, &start, &lower, &upper, step, 1e-15, 2000)?;
    Ok(found
        .value
        .is_finite()
        .then(|| [found.point[0], found.point[1]]))
}

fn distance_to_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let squared = dx * dx + dy * dy;
    let t = if squared > 0.0 {
        (((p.x - a.x) * dx + (p.y - a.y) * dy) / squared).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p.x - a.x - dx * t).hypot(p.y - a.y - dy * t)
}

fn in_triangle(t: [Point2; 3], p: Point2) -> bool {
    let cross =
        |a: Point2, b: Point2, c: Point2| (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    let (d1, d2, d3) = (
        cross(t[0], t[1], p),
        cross(t[1], t[2], p),
        cross(t[2], t[0], p),
    );
    let eps = 1e-12;
    let neg = d1 < -eps || d2 < -eps || d3 < -eps;
    let pos = d1 > eps || d2 > eps || d3 > eps;
    !(neg && pos)
}
