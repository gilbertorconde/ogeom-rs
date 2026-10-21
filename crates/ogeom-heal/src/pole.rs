//! Degenerate edges restored where a face's boundary crosses a pole.
//!
//! Where a surface collapses a whole row of its chart to one point (a
//! sphere's pole, a cone's apex), a boundary through that point arrives at
//! the row at one parameter and leaves it at another. In space the two are
//! the same place; in the chart they are a stretch of the row apart, and
//! the edge that spans it is a degenerate one: no curve in space, a pcurve
//! along the row. Without it the stretch is left to whoever reads the
//! chart, and on a row that repeats it is ambiguous: the boundary of a
//! sphere walking its seam up and down crosses each pole by a whole turn,
//! which a reader taking the shortest way across reads as no turn at all,
//! and the face meshes to nothing. `check` passes such a boundary, since
//! its two ends meet on the surface.

use ogeom_algo::edge_vertices;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Line2d, Surface as _, SurfaceGeometry};
use ogeom_math::{Point2, Vector2};
use ogeom_topo::{EdgeData, EdgeRepr, Location, Model, Shape, ShapeType, explore_unique};

use crate::Reshape;

/// The most joints at poles one wire is solved for; each has three
/// candidate turns to try.
const MOST_POLES: usize = 4;

/// Stage, for every wire whose crossings of poles are read wrongly the
/// shortest way across, a rebuilt wire with a degenerate edge spanning
/// each; returns how many edges were made.
///
/// A joint is at a pole where the end of one edge's pcurve and the start
/// of the next's stand apart in the chart along one parameter direction,
/// and the whole row of the chart through them, sampled across the
/// surface's domain, lifts onto the joint vertex within its tolerance: the
/// surface collapses that row. Every other joint must meet, or stand a
/// whole number of periods apart where the chart repeats (a seam's two
/// columns); a wire with any other gap is left as it is.
///
/// Across a repeating row each crossing may run the shortest way or a
/// whole turn either way more. The wire is solved for the one choice
/// under which its chart, unrolled, closes on itself and encloses area
/// the right way round: positive for a face's first wire, negative for
/// its holes. Where exactly one choice does and it is not the shortest
/// way at every pole, each pole gets a degenerate edge spanning its
/// crossing; where the shortest way is the answer, a reader of the chart
/// already reads it, and nothing is made. Only unplaced pcurves and
/// vertices are read; a placed boundary is left as it is.
pub(crate) fn stage_pole_edges(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let mut made = 0;
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let surface_id = data.surface;
        let Some(surface) = model.geometry().surface(surface_id).cloned() else {
            continue;
        };
        let stored = face.oriented(ogeom_topo::Orientation::Forward);
        for (w, wire) in model.ordered_children_of(&stored)?.into_iter().enumerate() {
            let edges = model.ordered_children_of(&wire)?;
            let Some(spans) = pole_spans(model, &surface, surface_id, &edges, w == 0, tol)? else {
                continue;
            };
            let mut ring = Vec::with_capacity(edges.len() + spans.len());
            let mut spans = spans.into_iter().peekable();
            for (i, edge) in edges.iter().enumerate() {
                ring.push(edge.clone());
                if let Some(pole) = spans.next_if(|pole| pole.after == i) {
                    ring.push(degenerate_edge(model, &pole, surface_id, tol)?);
                    made += 1;
                }
            }
            let rebuilt = model.add_wire(&ring)?;
            reshape.replace(&wire, rebuilt);
        }
    }
    Ok(made)
}

/// A degenerate edge to make: after which edge of the wire, at which
/// vertex, and its pcurve's ends.
struct Pole {
    after: usize,
    vertex: Shape,
    from: Point2,
    to: Point2,
}

/// One edge as the wire walks it in the chart: its points in order, first
/// and last its ends.
type Walked = Vec<Point2>;

/// The degenerate edges `edges` (a wire, `outer` where it is a face's
/// first) needs at its poles, or `None` where it needs none or cannot be
/// solved.
fn pole_spans(
    model: &Model,
    surface: &SurfaceGeometry,
    surface_id: ogeom_topo::SurfaceId,
    edges: &[Shape],
    outer: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Pole>>> {
    let n = edges.len();
    let mut walked: Vec<Walked> = Vec::with_capacity(n);
    for edge in edges {
        let Some(points) = walk(model, surface_id, edge, tol)? else {
            return Ok(None);
        };
        walked.push(points);
    }
    let ((u0, u1), (v0, v1)) = surface.domain();
    let periods = [
        (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(u1 - u0),
        (surface.is_periodic_v() || surface.is_closed_v(tol)).then_some(v1 - v0),
    ];
    // Each joint: met, at a pole, or a whole number of periods apart.
    let mut shifts = vec![Vector2::new(0.0, 0.0); n];
    let mut poles: Vec<(usize, Shape, usize, f64)> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let (end, start) = (walked[i][walked[i].len() - 1], walked[j][0]);
        let gap = start - end;
        if gap.x.abs() <= tol.parametric() && gap.y.abs() <= tol.parametric() {
            continue;
        }
        let axis = if gap.y.abs() <= tol.parametric() {
            Some(0)
        } else if gap.x.abs() <= tol.parametric() {
            Some(1)
        } else {
            None
        };
        if let Some(axis) = axis
            && let Some(vertex) =
                collapsed_row(model, surface, &edges[i], &edges[j], end, axis, tol)?
        {
            poles.push((i, vertex, axis, if axis == 0 { gap.x } else { gap.y }));
            continue;
        }
        let whole = |value: f64, period: Option<f64>| {
            value.abs() <= tol.parametric()
                || period
                    .is_some_and(|p| (value - (value / p).round() * p).abs() <= tol.parametric())
        };
        if whole(gap.x, periods[0]) && whole(gap.y, periods[1]) {
            shifts[i] = gap;
            continue;
        }
        return Ok(None);
    }
    if poles.is_empty() || poles.len() > MOST_POLES {
        return Ok(None);
    }
    // Each pole's crossing: the shortest way, or a turn either way more
    // where its row repeats.
    let shortest: Vec<f64> = poles
        .iter()
        .map(|(_, _, axis, gap)| periods[*axis].map_or(*gap, |p| gap - (gap / p).round() * p))
        .collect();
    let turns: Vec<&[i32]> = poles
        .iter()
        .map(|(_, _, axis, _)| -> &[i32] {
            if periods[*axis].is_some() {
                &[-1, 0, 1]
            } else {
                &[0]
            }
        })
        .collect();
    let mut solutions: Vec<Vec<f64>> = Vec::new();
    let mut choice = vec![0usize; poles.len()];
    loop {
        let crossing: Vec<f64> = (0..poles.len())
            .map(|k| {
                let p = periods[poles[k].2].unwrap_or(0.0);
                shortest[k] + f64::from(turns[k][choice[k]]) * p
            })
            .collect();
        if let Some(area) = unrolled_area(&walked, &shifts, &poles, &crossing, tol)
            && (if outer { area } else { -area }) > tol.parametric() * tol.parametric()
        {
            solutions.push(crossing);
        }
        // The next choice, odometer fashion.
        let mut k = 0;
        loop {
            if k == poles.len() {
                break;
            }
            choice[k] += 1;
            if choice[k] < turns[k].len() {
                break;
            }
            choice[k] = 0;
            k += 1;
        }
        if k == poles.len() {
            break;
        }
    }
    let [crossing] = &solutions[..] else {
        return Ok(None);
    };
    if crossing
        .iter()
        .zip(&shortest)
        .all(|(c, s)| (c - s).abs() <= tol.parametric())
    {
        return Ok(None);
    }
    Ok(Some(
        poles
            .into_iter()
            .zip(crossing)
            .map(|((i, vertex, axis, _), span)| {
                let from = walked[i][walked[i].len() - 1];
                let step = if axis == 0 {
                    Vector2::new(*span, 0.0)
                } else {
                    Vector2::new(0.0, *span)
                };
                Pole {
                    after: i,
                    vertex,
                    from,
                    to: from + step,
                }
            })
            .collect(),
    ))
}

/// The signed area the wire encloses in its chart unrolled: each edge
/// placed where the one before it, and the joint between, leave off,
/// whole-period joints closed and each pole crossed by its `crossing`.
/// `None` where the unrolled wire does not come back to where it began.
fn unrolled_area(
    walked: &[Walked],
    shifts: &[Vector2],
    poles: &[(usize, Shape, usize, f64)],
    crossing: &[f64],
    tol: Tolerances,
) -> Option<f64> {
    let mut offset = Vector2::new(0.0, 0.0);
    let mut path: Vec<Point2> = Vec::new();
    for (i, points) in walked.iter().enumerate() {
        path.extend(points.iter().map(|p| *p + offset));
        // Where the next edge starts, unrolled, against where it is stored.
        let end = points[points.len() - 1] + offset;
        let stored_next = walked[(i + 1) % walked.len()][0];
        let next = match poles.iter().position(|(at, ..)| *at == i) {
            Some(k) => {
                let step = if poles[k].2 == 0 {
                    Vector2::new(crossing[k], 0.0)
                } else {
                    Vector2::new(0.0, crossing[k])
                };
                end + step
            }
            None => end + (stored_next - points[points.len() - 1] - shifts[i]),
        };
        offset = next - stored_next;
    }
    if offset.x.abs() > tol.parametric() || offset.y.abs() > tol.parametric() {
        return None;
    }
    let mut twice = 0.0;
    for k in 0..path.len() {
        let (a, b) = (path[k], path[(k + 1) % path.len()]);
        twice += a.x.mul_add(b.y, -(b.x * a.y));
    }
    Some(twice / 2.0)
}

/// The vertex `edge` and `next` meet at, where the whole row of the chart
/// through `at` along `axis` lifts onto it within its tolerance.
fn collapsed_row(
    model: &Model,
    surface: &SurfaceGeometry,
    edge: &Shape,
    next: &Shape,
    at: Point2,
    axis: usize,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    /// Samples across the row.
    const PIECES: u32 = 16;
    let (Some((_, end)), Some((start, _))) =
        (edge_vertices(model, edge)?, edge_vertices(model, next)?)
    else {
        return Ok(None);
    };
    if !end.is_same(&start) || !end.location().is_identity() {
        return Ok(None);
    }
    let Some(vertex) = model.node(&end).and_then(|n| n.data().as_vertex()) else {
        return Ok(None);
    };
    let reach = vertex.tolerance.get().max(tol.confusion());
    let ((u0, u1), (v0, v1)) = surface.domain();
    let (lo, hi) = if axis == 0 { (u0, u1) } else { (v0, v1) };
    if !(lo.is_finite() && hi.is_finite()) {
        return Ok(None);
    }
    for k in 0..=PIECES {
        let t = (hi - lo).mul_add(f64::from(k) / f64::from(PIECES), lo);
        let uv = if axis == 0 {
            Point2::new(t, at.y)
        } else {
            Point2::new(at.x, t)
        };
        let Ok(lifted) = surface.point_at(uv.x, uv.y, tol) else {
            return Ok(None);
        };
        if lifted.distance(vertex.point) > reach {
            return Ok(None);
        }
    }
    Ok(Some(end))
}

/// `edge`'s pcurve on the surface as the wire walks it, sampled: a seam's
/// side by the occurrence's orientation. `None` for an edge with no
/// unplaced pcurve there.
fn walk(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Walked>> {
    /// Samples along each edge, for the area its wire encloses.
    const PIECES: u32 = 16;
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
    let (id, range) = match data.pcurve_for(surface_id, edge.location()) {
        Some(EdgeRepr::PCurve {
            curve,
            range,
            location,
            ..
        }) if location.is_identity() => (*curve, *range),
        Some(EdgeRepr::Seam {
            forward,
            reversed: back,
            range,
            location,
            ..
        }) if location.is_identity() => (if reversed { *back } else { *forward }, *range),
        _ => return Ok(None),
    };
    let Some(pcurve) = model.geometry().pcurve(id) else {
        return Ok(None);
    };
    let mut points = Vec::with_capacity(PIECES as usize + 1);
    for k in 0..=PIECES {
        let f = f64::from(k) / f64::from(PIECES);
        let f = if reversed { 1.0 - f } else { f };
        points.push(pcurve.point_at((range.1 - range.0).mul_add(f, range.0), tol)?);
    }
    Ok(Some(points))
}

/// A degenerate edge at `pole`'s vertex, its pcurve the chart segment it
/// spans.
fn degenerate_edge(
    model: &mut Model,
    pole: &Pole,
    surface: ogeom_topo::SurfaceId,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut data = EdgeData::new();
    data.degenerate = true;
    let edge = model.add_edge(data, &[pole.vertex.clone(), pole.vertex.clone()])?;
    ogeom_algo::attach_pcurve(
        model,
        &edge,
        Line2d::segment(pole.from, pole.to, tol)?.into(),
        surface,
        Location::identity(),
        (0.0, pole.from.distance(pole.to)),
    )?;
    Ok(edge)
}
