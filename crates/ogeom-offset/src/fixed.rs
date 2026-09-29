//! A section carried along a spine by translation alone.
//!
//! Every point of the section runs along its own copy of the spine, so each
//! edge of the section sweeps the surface `e(u) + c(v) - c(start)`. For
//! rational B-splines that sum is itself a tensor-product rational
//! B-spline, exactly: control points `E_i + C_j` less the start, weights
//! `w_i v_j`, since the product of the two denominators is the tensor
//! product of their weights. Nothing is fitted. The section's vertices run
//! along translated copies of the spine, which the side faces share, and
//! the section and its copy at the spine's end cap the ends.

use std::collections::HashMap;

use ogeom_algo::{Built, make_edge_between, make_face_on, make_shell, make_solid, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Transformable as _;
use ogeom_geom::{BSplineCurve, BSplineSurface, Curve, Curve3d as _, Line2d, PlaneSurface};
use ogeom_math::{
    Axis2, ControlGrid, Direction, Direction2, Frame, Plane, Point, Point2, Vector, Vector2,
    Weighted,
};
use ogeom_topo::{EdgeRepr, Location, Model, Orientation, Shape, ShapeType, TShapeId};

/// A section edge, placed and spelt as a B-spline over its own parameter.
struct SectionEdge {
    /// The occurrence in its ring, for the direction it is walked.
    occurrence: Shape,
    spline: BSplineCurve,
    /// The stored edge's start and end vertices.
    ends: (Shape, Shape),
}

/// Sweep `profile` along `spine`, carried by translation: the section keeps
/// the orientation it has at the spine's start all the way.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// profile is not a planar face or closed wire, the spine is not an edge or
/// wire, or the spine crosses the section's plane one way and back the
/// other, where the section would sweep through itself.
pub(crate) fn fixed_pipe(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let rings: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => model.ordered_children_of(profile)?,
        ShapeType::Wire => vec![profile.clone()],
        _ => ogeom_bail!(Construction, "a pipe sweeps a planar face or wire"),
    };
    let mut section: Vec<Vec<SectionEdge>> = Vec::with_capacity(rings.len());
    for ring in &rings {
        let mut edges = Vec::new();
        for occurrence in model.ordered_children_of(ring)? {
            edges.push(section_edge(model, &occurrence, tol)?);
        }
        section.push(edges);
    }
    let normal = section_normal(&section, tol)?;

    let legs: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        ShapeType::Wire => model.ordered_children_of(spine)?,
        _ => ogeom_bail!(Construction, "a spine is an edge or a wire"),
    };
    let mut start: Option<Point> = None;
    let mut pieces = Vec::with_capacity(legs.len());
    for leg in &legs {
        let (spline, from_start) = leg_spline(model, leg, tol)?;
        let (a, b) = spline.domain();
        let (v_start, v_end) = if from_start { (a, b) } else { (b, a) };
        let origin = *start.get_or_insert(spline.point_at(v_start, tol)?);
        pieces.push(leg_solid(
            model,
            &section,
            normal,
            &spline,
            (v_start, v_end),
            origin,
            tol,
        )?);
    }
    let mut runs = pieces;
    while runs.len() > 1 {
        let tail = runs.split_off(1);
        let mut joined = runs.remove(0);
        for next in tail {
            joined = ogeom_bool::fuse(model, &joined, &next, tol)?.shape;
        }
        runs = vec![joined];
    }
    let Some(solid) = runs.pop() else {
        ogeom_bail!(Construction, "a spine of no edges carries nothing");
    };
    let mut history = ogeom_algo::History::new();
    history.generate(spine, solid.clone());
    history.generate(profile, solid.clone());
    Ok(Built::new(solid, history))
}

/// A section edge's curve, placed where the occurrence stands, as a spline.
fn section_edge(model: &Model, occurrence: &Shape, tol: Tolerances) -> OgeomResult<SectionEdge> {
    let Some(data) = model.node(occurrence).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "a section edge has no curve");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placed = geometry.transformed(&occurrence.transform(model.datums())?, tol)?;
    let spline = placed.to_bspline_over(*range, tol)?;
    let bounds = model.children_of(&stored(occurrence))?;
    let (Some(first), Some(last)) = (bounds.first(), bounds.last()) else {
        ogeom_bail!(Construction, "a section edge has no vertices");
    };
    Ok(SectionEdge {
        occurrence: occurrence.clone(),
        spline,
        ends: (first.clone(), last.clone()),
    })
}

/// A spine edge as a spline, and whether its occurrence walks it from its
/// spline's start.
fn leg_spline(model: &Model, leg: &Shape, tol: Tolerances) -> OgeomResult<(BSplineCurve, bool)> {
    let Some(data) = model.node(leg).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "a spine edge has no curve");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placed = geometry.transformed(&leg.transform(model.datums())?, tol)?;
    Ok((
        placed.to_bspline_over(*range, tol)?,
        leg.orientation() != Orientation::Reversed,
    ))
}

/// The occurrence in its stored direction.
fn stored(shape: &Shape) -> Shape {
    if shape.orientation() == Orientation::Reversed {
        shape.reversed()
    } else {
        shape.clone()
    }
}

/// The normal about which the section's outer ring runs anticlockwise:
/// the material of the section lies to the left of every ring's walk
/// about it.
fn section_normal(section: &[Vec<SectionEdge>], tol: Tolerances) -> OgeomResult<Vector> {
    let Some(outer) = section.first() else {
        ogeom_bail!(Construction, "a section with no ring bounds nothing");
    };
    let mut points = Vec::new();
    for edge in outer {
        let (a, b) = edge.spline.domain();
        let (t0, t1) = if edge.occurrence.orientation() == Orientation::Reversed {
            (b, a)
        } else {
            (a, b)
        };
        for k in 0..16 {
            let t = t0 + (t1 - t0) * f64::from(k) / 16.0;
            points.push(edge.spline.point_at(t, tol)?);
        }
    }
    // Newell's normal: twice the ring's vector area.
    let mut n = Vector::ZERO;
    for (k, p) in points.iter().enumerate() {
        let q = points[(k + 1) % points.len()];
        n += (*p - Point::ORIGIN).cross(q - Point::ORIGIN);
    }
    if n.magnitude() <= tol.confusion() {
        ogeom_bail!(Construction, "the section encloses no area");
    }
    Ok(n / n.magnitude())
}

/// `spline` moved by `by`: the same parameter and weights.
fn moved(spline: &BSplineCurve, by: Vector, tol: Tolerances) -> OgeomResult<BSplineCurve> {
    let control = spline
        .control_points()
        .iter()
        .map(|c| Weighted::new(c.point() + by, c.weight, tol))
        .collect::<OgeomResult<Vec<_>>>()?;
    BSplineCurve::rational(spline.knots().clone(), control)
}

/// The surface `edge(u) + spine(v) - origin`, exactly.
fn swept(
    edge: &BSplineCurve,
    spine: &BSplineCurve,
    origin: Point,
    tol: Tolerances,
) -> OgeomResult<BSplineSurface> {
    let mut grid = Vec::with_capacity(edge.control_points().len() * spine.control_points().len());
    for e in edge.control_points() {
        for c in spine.control_points() {
            let at = e.point() + (c.point() - origin);
            grid.push(Weighted::new(at, e.weight * c.weight, tol)?);
        }
    }
    let (nu, nv) = (edge.control_points().len(), spine.control_points().len());
    BSplineSurface::rational(
        edge.knots().clone(),
        spine.knots().clone(),
        ControlGrid::new(grid, nu, nv)?,
    )
}

/// A straight pcurve through `from` along `along`, parameterized as the
/// edge it describes over `range`.
fn chart_line(
    from: Point2,
    along: Vector2,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::PlanarCurve> {
    Ok(Line2d::over(
        Axis2::new(from, Direction2::new(along, tol)?),
        range.0,
        range.1,
    )?
    .into())
}

/// One spine edge's solid: the section's sides swept along it, capped by
/// the section where the edge starts and where it ends.
fn leg_solid(
    model: &mut Model,
    section: &[Vec<SectionEdge>],
    normal: Vector,
    spine: &BSplineCurve,
    (v_start, v_end): (f64, f64),
    origin: Point,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let spine_range = spine.domain();
    let (to_start, to_end) = (
        spine.point_at(v_start, tol)? - origin,
        spine.point_at(v_end, tol)? - origin,
    );
    // Which way the leg carries the section through its plane: one way
    // throughout, or the section sweeps through itself. Where the spine
    // runs in the plane (its tangent square to the normal) it carries the
    // section along itself, which bounds nothing either way.
    let heading = if v_end >= v_start { 1.0 } else { -1.0 };
    let (mut ahead, mut behind) = (0.0_f64, 0.0_f64);
    for k in 0..=32 {
        let v = spine_range.0 + (spine_range.1 - spine_range.0) * f64::from(k) / 32.0;
        let d = spine.derivatives_at(v, 1, tol)?[1] * heading;
        let across = d.dot(normal) / d.magnitude().max(f64::MIN_POSITIVE);
        ahead = ahead.max(across);
        behind = behind.max(-across);
    }
    if ahead > tol.angular() && behind > tol.angular() {
        ogeom_bail!(
            Construction,
            "the spine crosses the section's plane one way and back the \
             other; the fixed section would sweep through itself"
        );
    }
    if ahead <= tol.angular() && behind <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the spine runs in the section's plane; a fixed section carried \
             along it sweeps no volume"
        );
    }
    let forward = ahead > behind;

    // The section's vertices at both ends, and the rail each runs along.
    struct Rail {
        start: Shape,
        end: Shape,
        edge: Shape,
    }
    let mut rails: HashMap<TShapeId, Rail> = HashMap::new();
    for edge in section.iter().flatten() {
        for vertex in [&edge.ends.0, &edge.ends.1] {
            if rails.contains_key(&vertex.node()) {
                continue;
            }
            let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Dangling, "vertex is not in this model");
            };
            let at = vertex.transform(model.datums())?.apply(data.point);
            let start = make_vertex(model, at + to_start).shape;
            let end = make_vertex(model, at + to_end).shape;
            let path = moved(spine, at - origin, tol)?;
            let (low, high) = if v_start <= v_end {
                (&start, &end)
            } else {
                (&end, &start)
            };
            let edge =
                make_edge_between(model, Curve::from(path), spine_range, low, high, tol)?.shape;
            rails.insert(vertex.node(), Rail { start, end, edge });
        }
    }

    // The section's edges moved to both ends, shared by caps and sides.
    let mut start_rings: Vec<Vec<Shape>> = Vec::with_capacity(section.len());
    let mut end_rings: Vec<Vec<Shape>> = Vec::with_capacity(section.len());
    let mut faces: Vec<Shape> = Vec::new();
    for ring in section {
        let (mut starts, mut ends) = (Vec::new(), Vec::new());
        for edge in ring {
            let range = edge.spline.domain();
            let rail_a = &rails[&edge.ends.0.node()];
            let rail_b = &rails[&edge.ends.1.node()];
            let at_start = make_edge_between(
                model,
                Curve::from(moved(&edge.spline, to_start, tol)?),
                range,
                &rail_a.start,
                &rail_b.start,
                tol,
            )?
            .shape;
            let at_end = make_edge_between(
                model,
                Curve::from(moved(&edge.spline, to_end, tol)?),
                range,
                &rail_a.end,
                &rail_b.end,
                tol,
            )?
            .shape;

            // The side: bounded by the edge at both ends and the two rails,
            // every pcurve a line of the surface's own chart.
            let surface = swept(&edge.spline, spine, origin, tol)?;
            let surface_id = model.geometry_mut().add_surface(surface.clone().into());
            let identity = Location::identity();
            ogeom_algo::attach_pcurve(
                model,
                &at_start,
                chart_line(
                    Point2::new(0.0, v_start),
                    Vector2::new(1.0, 0.0),
                    range,
                    tol,
                )?,
                surface_id,
                identity.clone(),
                range,
            )?;
            ogeom_algo::attach_pcurve(
                model,
                &at_end,
                chart_line(Point2::new(0.0, v_end), Vector2::new(1.0, 0.0), range, tol)?,
                surface_id,
                identity.clone(),
                range,
            )?;
            for (rail, u) in [(rail_a, range.0), (rail_b, range.1)] {
                ogeom_algo::attach_pcurve(
                    model,
                    &rail.edge,
                    chart_line(
                        Point2::new(u, 0.0),
                        Vector2::new(0.0, 1.0),
                        spine_range,
                        tol,
                    )?,
                    surface_id,
                    identity.clone(),
                    spine_range,
                )?;
            }
            // Round the side: along the edge at the start, up the far rail,
            // back along the edge at the end, down the near rail; run so it
            // goes anticlockwise in the chart, about the chart's own normal.
            let mut loop_edges = vec![
                at_start.clone(),
                rail_b.edge.clone(),
                at_end.reversed(),
                rail_a.edge.reversed(),
            ];
            if v_start > v_end {
                // The rails climb in the chart from the start row to the
                // end row; walked the other way the loop turns clockwise.
                loop_edges = vec![
                    at_start.clone(),
                    rail_b.edge.reversed(),
                    at_end.reversed(),
                    rail_a.edge.clone(),
                ];
                loop_edges = loop_edges.iter().rev().map(Shape::reversed).collect();
            }
            let wire = ogeom_algo::make_wire(model, &loop_edges, tol)?.shape;
            let side = make_face_on(model, surface_id, &[wire], tol)?.shape;
            // The chart's normal is the edge's parameter direction across the
            // spine's; it points out of the solid where the section's walk
            // and the spine's crossing of the plane agree.
            let walks_along = edge.occurrence.orientation() != Orientation::Reversed;
            let crosses_along = (heading > 0.0) == forward;
            faces.push(if walks_along == crosses_along {
                side
            } else {
                side.reversed()
            });
            let place = |shape: Shape| {
                if edge.occurrence.orientation() == Orientation::Reversed {
                    shape.reversed()
                } else {
                    shape
                }
            };
            starts.push(place(at_start));
            ends.push(place(at_end));
        }
        start_rings.push(starts);
        end_rings.push(ends);
    }

    // The caps: the section where the leg starts, facing back along it, and
    // where it ends, facing on.
    for (rings, at, out) in [
        (&start_rings, to_start, !forward),
        (&end_rings, to_end, forward),
    ] {
        let frame = Frame::new(
            origin + at,
            Direction::new(normal, tol)?,
            Direction::new(any_square_to(normal), tol)?,
            tol,
        )?;
        let plane: ogeom_geom::SurfaceGeometry = PlaneSurface::new(Plane::new(frame)).into();
        let cap = ogeom_algo::make_face_with_pcurves(model, plane, rings, tol)?.shape;
        faces.push(if out { cap } else { cap.reversed() });
    }
    let shell = make_shell(model, &faces)?.shape;
    Ok(make_solid(model, &[shell])?.shape)
}

/// A unit vector square to `n`.
fn any_square_to(n: Vector) -> Vector {
    let helper = if n.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let s = n.cross(helper);
    s / s.magnitude()
}
