//! A solid's sphere faces charted again with the sphere's pole turned away
//! from a blend's seat.
//!
//! A rail passing a sphere's pole closer than a stride of the march turns
//! half the chart round within that stride, and no image fitted through
//! the stations follows it, on the leg or in the boolean's shared chart.
//! The same sphere on a frame turned so its pole stands clear of the seat
//! is the same surface: the faces on it are rebuilt there, each bounded by
//! its own edges with the seam and pole edges of the old chart dropped and
//! the new chart's own put in, and the leg is built on that chart too.

use ogeom_algo::History;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, Curve2d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Point2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType};

/// The solid with its faces on a sphere charted afresh, and the chart.
pub(crate) struct TurnedPole {
    /// The sphere on its turned frame.
    pub(crate) host: SurfaceGeometry,
    /// The solid, its faces on the sphere rebuilt on `host`.
    pub(crate) solid: Shape,
    /// What became of the solid's faces.
    pub(crate) history: History,
}

/// Samples per edge where the edges' clearance and columns are read.
const PER_EDGE: u32 = 32;

/// The solid's faces on the sphere `host` rebuilt on the same sphere with
/// its pole on whichever axis of its own frame stands farthest from the
/// seat (`seat`, and the faces' own edges), turned about that axis so its
/// seam runs through the widest gap between them.
///
/// `None` where no other axis stands farther than the pole's own, the host
/// is no sphere or a mirrored one, no face of the solid lies on it, or a
/// face's own edges do not bound it there as a sphere less some disks, or
/// one disk less some others: an edge wound round the new pole, crossing
/// its seam, or with no image the chart holds.
pub(crate) fn turned_from_pole(
    model: &mut Model,
    solid: &Shape,
    host: &SurfaceGeometry,
    seat: &[Point],
    tol: Tolerances,
) -> OgeomResult<Option<TurnedPole>> {
    let SurfaceGeometry::Sphere(sphere) = host else {
        return Ok(None);
    };
    let sphere = sphere.sphere();
    let frame = sphere.frame();
    if frame.handedness() != ogeom_math::Handedness::Right {
        return Ok(None);
    }
    let (centre, radius) = (frame.origin(), sphere.radius());

    // The faces on the host, each with its own edges walked as loops.
    let mut faces: Vec<(Shape, Vec<Vec<Shape>>)> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let on_host = matches!(
            model.geometry().surface(data.surface),
            Some(SurfaceGeometry::Sphere(s))
                if s.sphere().frame() == frame && s.sphere().radius() == radius
        );
        if !on_host {
            continue;
        }
        if !data.location.is_identity() || !face.location().is_identity() {
            return Ok(None);
        }
        let face = face.oriented(ogeom_topo::Orientation::Forward);
        let Some(loops) = own_loops(model, &face, data.surface)? else {
            return Ok(None);
        };
        faces.push((face, loops));
    }
    if faces.is_empty() {
        return Ok(None);
    }

    let mut points: Vec<Point> = seat.to_vec();
    for (_, loops) in &faces {
        for edge in loops.iter().flatten() {
            let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
            for k in 0..=PER_EDGE {
                let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(PER_EDGE);
                points.push(curve.point_at(t, tol)?);
            }
        }
    }
    let clearance = |axis: Direction| -> f64 {
        let reach = axis.vector() * radius;
        points.iter().fold(f64::INFINITY, |near, p| {
            near.min(p.distance(centre + reach).min(p.distance(centre - reach)))
        })
    };
    let (x, y, z) = (frame.x(), frame.y(), frame.z());
    let own = clearance(z);
    let mut best: Option<(f64, [Direction; 3])> = None;
    // Turned cyclically, so the frame stays right-handed and the sphere's
    // outward sense is kept.
    for axes in [[y, z, x], [z, x, y]] {
        let far = clearance(axes[2]);
        if far > own && best.is_none_or(|(seen, _)| far > seen) {
            best = Some((far, axes));
        }
    }
    let Some((_, [x, y, z])) = best else {
        return Ok(None);
    };

    // The seam through the middle of the widest gap between the columns the
    // points stand on.
    let mut columns: Vec<f64> = points
        .iter()
        .map(|p| {
            let off = *p - centre;
            off.dot(y.vector()).atan2(off.dot(x.vector()))
        })
        .collect();
    columns.sort_by(f64::total_cmp);
    let mut gap = (
        columns[0] + core::f64::consts::TAU - columns[columns.len() - 1],
        columns[columns.len() - 1],
    );
    for w in columns.windows(2) {
        if w[1] - w[0] > gap.0 {
            gap = (w[1] - w[0], w[0]);
        }
    }
    let middle = gap.1 + gap.0 / 2.0;
    let seam = x.vector() * middle.cos() + y.vector() * middle.sin();
    let turned = Frame::new(centre, z, Direction::new(seam, tol)?, tol)?;

    let turned_host: SurfaceGeometry =
        ogeom_geom::SphereSurface::new(ogeom_math::Sphere::new(turned, radius, tol)?).into();
    let surface_id = model.geometry_mut().add_surface(turned_host.clone());
    let ((u0, u1), _) = turned_host.domain();
    let period = u1 - u0;

    let mut reshape = ogeom_heal::Reshape::new();
    for (face, loops) in faces {
        let mut outer: Option<Shape> = None;
        let mut holes = Vec::with_capacity(loops.len());
        for edges in loops {
            let mut walk: Vec<Point2> = Vec::new();
            for edge in &edges {
                let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
                let Some(image) = image_on(&curve, range, &turned_host, tol) else {
                    return Ok(None);
                };
                // Moved by whole periods into the window: the seam stands
                // clear of every edge, so each image lies in one copy.
                let middle = image.point_at(f64::midpoint(range.0, range.1), tol)?.x;
                let turns = ((middle - u0) / period).floor();
                let image = if turns == 0.0 {
                    image
                } else {
                    image.transformed(
                        &ogeom_math::Transform2::translation(ogeom_math::Vector2::new(
                            -turns * period,
                            0.0,
                        )),
                        tol,
                    )?
                };
                let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
                for k in 0..PER_EDGE {
                    let s = f64::from(k) / f64::from(PER_EDGE);
                    let t = if reversed {
                        range.1 - (range.1 - range.0) * s
                    } else {
                        range.0 + (range.1 - range.0) * s
                    };
                    let uv = image.point_at(t, tol)?;
                    if uv.x <= u0 || uv.x >= u1 {
                        return Ok(None);
                    }
                    walk.push(uv);
                }
                ogeom_algo::attach_pcurve(
                    model,
                    edge,
                    image,
                    surface_id,
                    ogeom_topo::Location::identity(),
                    range,
                )?;
            }
            let wire = ogeom_algo::make_wire(model, &edges, tol)?.shape;
            // Anticlockwise in the chart bounds the face from outside, and
            // clockwise is a hole; a face has one outside at most.
            if signed_area(&walk) > 0.0 {
                if outer.replace(wire).is_some() {
                    return Ok(None);
                }
            } else {
                holes.push(wire);
            }
        }
        let outer = match outer {
            Some(outer) => outer,
            None => rectangle(model, surface_id, turned, radius, tol)?,
        };
        let wires: Vec<Shape> = core::iter::once(outer).chain(holes).collect();
        let rebuilt = ogeom_algo::make_face_on(model, surface_id, &wires, tol)?.shape;
        if let Some(stated) = model.tolerance_of(&face)? {
            model.widen(&rebuilt, stated)?;
        }
        reshape.replace(&face, rebuilt);
    }
    let built = reshape.apply(model, solid)?;
    Ok(Some(TurnedPole {
        host: turned_host,
        solid: built.shape,
        history: built.history,
    }))
}

/// The whole sphere's boundary in its chart, for a face its own edges leave
/// open round the poles: the meridian through the frame's `x` from pole to
/// pole as the seam, and each pole as an edge with no length, walked
/// anticlockwise round the chart's rectangle.
fn rectangle(
    model: &mut Model,
    surface_id: ogeom_topo::SurfaceId,
    frame: Frame,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use core::f64::consts::{FRAC_PI_2, TAU};
    let centre = frame.origin();
    let south = ogeom_algo::make_vertex(model, centre - frame.z().vector() * radius).shape;
    let north = ogeom_algo::make_vertex(model, centre + frame.z().vector() * radius).shape;
    // In the plane of `x` and `z`, its normal `-y`: the circle's angle is
    // the latitude.
    let meridian =
        ogeom_math::Circle::new(Frame::new(centre, -frame.y(), frame.x(), tol)?, radius, tol)?;
    let rows = (-FRAC_PI_2, FRAC_PI_2);
    let seam = ogeom_algo::make_edge_between(
        model,
        ogeom_geom::CircleCurve::new(meridian).into(),
        rows,
        &south,
        &north,
        tol,
    )?
    .shape;
    let upright = |u: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            rows.0,
            rows.1,
        )?
        .into())
    };
    ogeom_algo::attach_seam(
        model,
        &seam,
        upright(TAU)?,
        upright(0.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        rows,
    )?;
    let mut poles = Vec::with_capacity(2);
    for (pole, row) in [(&south, rows.0), (&north, rows.1)] {
        let mut data = ogeom_topo::EdgeData::new();
        data.degenerate = true;
        let edge = model.add_edge(data, &[pole.clone(), pole.clone()])?;
        let image: PlanarCurve = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, row), ogeom_math::Direction2::X),
            0.0,
            TAU,
        )?
        .into();
        ogeom_algo::attach_pcurve(
            model,
            &edge,
            image,
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, TAU),
        )?;
        poles.push(edge);
    }
    let ring = [
        poles[0].clone(),
        seam.clone(),
        poles[1].reversed(),
        seam.reversed(),
    ];
    Ok(ogeom_algo::make_wire(model, &ring, tol)?.shape)
}

/// A face's own edges, its seam and pole edges dropped, walked head to tail
/// as closed loops in the face's sense; `None` where an edge bounds the
/// face twice without being its chart's seam, or the rest do not close.
fn own_loops(
    model: &Model,
    face: &Shape,
    surface: ogeom_topo::SurfaceId,
) -> OgeomResult<Option<Vec<Vec<Shape>>>> {
    let mut walked: Vec<Shape> = Vec::new();
    for wire in model.children_of(face)? {
        walked.extend(model.children_of(&wire)?);
    }
    let mut own: Vec<Shape> = Vec::with_capacity(walked.len());
    for edge in &walked {
        let Some(NodeData::Edge(data)) = model.node(edge).map(ogeom_topo::TShape::data) else {
            return Ok(None);
        };
        if data.degenerate {
            continue;
        }
        if walked.iter().filter(|e| e.node() == edge.node()).count() > 1 {
            if matches!(data.pcurve_on(surface), Some(EdgeRepr::Seam { .. })) {
                continue;
            }
            return Ok(None);
        }
        own.push(edge.clone());
    }
    let mut ends = Vec::with_capacity(own.len());
    for edge in &own {
        let Some((from, to)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(None);
        };
        ends.push((from.node(), to.node()));
    }
    let mut used = vec![false; own.len()];
    let mut loops = Vec::new();
    while let Some(first) = used.iter().position(|u| !u) {
        used[first] = true;
        let mut chain = vec![own[first].clone()];
        let start = ends[first].0;
        let mut at = ends[first].1;
        while at != start {
            let mut next = (0..own.len()).filter(|&i| !used[i] && ends[i].0 == at);
            let (Some(i), None) = (next.next(), next.next()) else {
                return Ok(None);
            };
            used[i] = true;
            chain.push(own[i].clone());
            at = ends[i].1;
        }
        loops.push(chain);
    }
    Ok(Some(loops))
}

/// An edge's image on the turned sphere: in closed form where there is
/// one, fitted to its target otherwise.
fn image_on(
    curve: &Curve,
    range: (f64, f64),
    host: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    if let Some(exact) = ogeom_intersect::exact_pcurve_of(curve, host, tol) {
        return Some(exact);
    }
    let (image, _, met, _, _) =
        ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, range, host, tol).ok()?;
    met.then_some(image)
}

/// The signed area a closed chart polyline encloses, by the shoelace about
/// its first point.
fn signed_area(points: &[Point2]) -> f64 {
    let Some(&anchor) = points.first() else {
        return 0.0;
    };
    let mut twice = 0.0;
    for (a, b) in points.iter().zip(points.iter().cycle().skip(1)) {
        let (a, b) = (*a - anchor, *b - anchor);
        twice += a.x * b.y - a.y * b.x;
    }
    twice / 2.0
}
