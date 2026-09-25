//! The ruled fillet: a straight crease whose faces are planes or drums
//! running parallel to it, blended exactly.
//!
//! A wall meeting a drum along one of its rulings (the vertical edge where
//! an extruded profile's line meets its arc at an angle) or two drums
//! meeting along a shared ruling: every section square to the edge is the
//! same, so the blend is a section swept straight. The section is a plane
//! problem: each host cuts the section plane in a line or a circle, the
//! ball in a circle tangent to both, and the wedge is the region between
//! the crease point, the two touch points and the ball's arc. Extruded
//! along the edge, its faces are planes and drums: the legs lie exactly on
//! the hosts and the band is a drum of the fillet's radius, so the melt is
//! exact and a band's cap meets a corner ball's rim on the same circle.

use crate::fillet::Mate;
use crate::march::{Sides, seat_section};
use crate::support::{edge_curve, same_occurrence};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore};

/// A host's trace in the section plane: a line or a circle about a point.
enum Trace {
    Line,
    Circle { centre_on_axis: Point, radius: f64 },
}

/// Blend a straight `edge` whose two faces are planes or drums parallel to
/// it, at least one a drum; `None` where the seat is not of that kind.
pub(crate) fn ruled_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    mates: Option<(usize, &[Mate])>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let (guide, range) = edge_curve(model, edge, tol)?;
    let Curve::Line(_) = &guide else {
        return Ok(None);
    };
    let start = guide.point_at(range.0, tol)?;
    let end = guide.point_at(range.1, tol)?;
    let length = start.distance(end);
    if length <= tol.confusion() {
        ogeom_bail!(Construction, "the edge has no length to blend along");
    }
    let along = (end - start) / length;

    let mut faces: Vec<Shape> = Vec::new();
    let mut hosts: Vec<(SurfaceGeometry, f64, Trace)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face().cloned()) else {
            ogeom_bail!(Construction, "a face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface).cloned() else {
            ogeom_bail!(Dangling, "a face's surface is not in this model");
        };
        use ogeom_geom::Transformable as _;
        let surface = stored.transformed(&face.transform(model.datums())?, tol)?;
        let trace = match &surface {
            SurfaceGeometry::Plane(_) => Trace::Line,
            SurfaceGeometry::Cylinder(c) => {
                let cylinder = c.cylinder();
                let axis = cylinder.frame().z().vector();
                if axis.cross(along).magnitude() > tol.angular() {
                    return Ok(None);
                }
                // The axis's crossing with the section plane at the start.
                let o = cylinder.frame().origin();
                Trace::Circle {
                    centre_on_axis: o + axis * (start - o).dot(axis),
                    radius: cylinder.radius(),
                }
            }
            _ => return Ok(None),
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        faces.push(face);
        hosts.push((surface, sign, trace));
    }
    if hosts.len() != 2 {
        return Ok(None);
    }
    if hosts.iter().all(|(_, _, t)| matches!(t, Trace::Line)) {
        return Ok(None);
    }

    let convex = crate::marched::crease_convexity(
        model,
        &faces[0],
        [(&hosts[0].0, hosts[0].1), (&hosts[1].0, hosts[1].1)],
        &guide,
        range,
        radius,
        tol,
    )?;
    let seat_sign = if convex { -1.0 } else { 1.0 };
    #[allow(clippy::cast_possible_truncation)]
    let sides = Sides {
        first: (seat_sign * hosts[0].1) as i8,
        second: (seat_sign * hosts[1].1) as i8,
    };
    let near = {
        let a = ogeom_algo::project_on_surface(&hosts[0].0, start, 32, tol)?.parameters;
        let b = ogeom_algo::project_on_surface(&hosts[1].0, start, 32, tol)?.parameters;
        [a.0, a.1, b.0, b.1]
    };
    let x = seat_section(
        &hosts[0].0,
        &hosts[1].0,
        radius,
        &guide,
        sides,
        range.0,
        near,
        tol,
    )?;
    // Every section is the same, so the start's is carried into the plane
    // square to the edge through the start.
    let square = |p: Point| p - along * (p - start).dot(along);
    use ogeom_geom::Surface as _;
    let touch = [
        square(hosts[0].0.point_at(x[0], x[1], tol)?),
        square(hosts[1].0.point_at(x[2], x[3], tol)?),
    ];
    let centre = {
        let (du, dv) = hosts[0].0.d1_at(x[0], x[1], tol)?;
        let n = du.cross(dv);
        square(touch[0] + n / n.magnitude() * (f64::from(sides.first) * radius))
    };

    // How far past each end the band runs: nowhere past a flush end, and
    // out through a neighbouring blend until the ball has left the solid,
    // the way the planar seat runs out.
    let mut reach = [0.0_f64; 2];
    if let Some((index, mates)) = mates {
        for (slot, at_end) in [(0, false), (1, true)] {
            let at = if at_end { end } else { start };
            let outward = if at_end { along } else { -along };
            let chain = mates.iter().enumerate().any(|(i, mate)| {
                i != index
                    && mate.ends.iter().any(|(p, leaving)| {
                        p.distance(at) <= tol.confusion() * 1e3
                            && leaving.cross(outward).magnitude() <= 1e-2
                            && leaving.dot(outward) > 0.0
                    })
            });
            if chain || Mate::settled_at(mates, index, at, tol) {
                continue;
            }
            let pair = [&faces[0], &faces[1]];
            if !crate::marched::crease_terminates_at(model, solid, edge, pair, at, tol)?
                || !crate::marched::neighbour_blend_at(model, solid, pair, at, convex, radius, tol)?
            {
                continue;
            }
            let centre_at_end = centre + along * (at - start).dot(along);
            let deflection = ogeom_mesh::Deflection {
                chord: (radius * 1e-2).max(tol.confusion() * 1e3),
                ..ogeom_mesh::Deflection::default()
            };
            let step = radius / 8.0;
            for k in 1..=32 {
                let s = step * f64::from(k);
                let inside = ogeom_algo::classify_in_solid(
                    model,
                    solid,
                    centre_at_end + outward * s,
                    deflection,
                    tol,
                )? == ogeom_algo::Containment::In;
                if inside != convex {
                    reach[slot] = s + radius * 0.25;
                    break;
                }
            }
        }
    }

    // The section at the band's own start, and the wedge swept from it.
    let shift = -along * reach[0];
    let crease = start + shift;
    let touch = [touch[0] + shift, touch[1] + shift];
    let centre = centre + shift;
    let v_crease = ogeom_algo::make_vertex(model, crease).shape;
    let v_touch = [
        ogeom_algo::make_vertex(model, touch[0]).shape,
        ogeom_algo::make_vertex(model, touch[1]).shape,
    ];
    let arc = |model: &mut Model,
               about: Point,
               r: f64,
               from: (&Shape, Point),
               to: (&Shape, Point)|
     -> OgeomResult<Shape> {
        // Swept positively about whichever of ±along turns the short way.
        let x_dir = (from.1 - about) / r;
        let y_dir = along.cross(x_dir);
        let angle = (to.1 - about).dot(y_dir).atan2((to.1 - about).dot(x_dir));
        let normal = if angle >= 0.0 { along } else { -along };
        let frame = Frame::new(
            about,
            Direction::new(normal, tol)?,
            Direction::new(x_dir, tol)?,
            tol,
        )?;
        let circle: Curve = CircleCurve::new(Circle::new(frame, r, tol)?).into();
        Ok(
            ogeom_algo::make_edge_between(model, circle, (0.0, angle.abs()), from.0, to.0, tol)?
                .shape,
        )
    };
    let leg =
        |model: &mut Model, trace: &Trace, i: usize| -> OgeomResult<Shape> {
            match trace {
                Trace::Line => {
                    let line: Curve = LineCurve::segment(crease, touch[i], tol)?.into();
                    let domain = line.domain();
                    Ok(ogeom_algo::make_edge_between(
                        model,
                        line,
                        domain,
                        &v_crease,
                        &v_touch[i],
                        tol,
                    )?
                    .shape)
                }
                Trace::Circle {
                    centre_on_axis,
                    radius: r,
                } => arc(
                    model,
                    *centre_on_axis + shift,
                    *r,
                    (&v_crease, crease),
                    (&v_touch[i], touch[i]),
                ),
            }
        };
    let leg_first = leg(model, &hosts[0].2, 0)?;
    let leg_second = leg(model, &hosts[1].2, 1)?;
    let ball_arc = arc(
        model,
        centre,
        radius,
        (&v_touch[0], touch[0]),
        (&v_touch[1], touch[1]),
    )?;
    // Wound counter-clockwise about the section plane's normal.
    let turn = (touch[0] - crease).cross(touch[1] - crease).dot(along);
    let normal = if turn >= 0.0 { along } else { -along };
    let reach_plane = (radius * 4.0).max(crease.distance(centre) * 4.0).max(1.0);
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(crease, Direction::new(normal, tol)?),
        (-reach_plane, reach_plane),
        (-reach_plane, reach_plane),
    )?
    .into();
    let loop_edges = if turn >= 0.0 {
        vec![leg_first, ball_arc, leg_second.reversed()]
    } else {
        vec![leg_second, ball_arc.reversed(), leg_first.reversed()]
    };
    let section = ogeom_algo::make_face_with_pcurves(model, plane, &[loop_edges], tol)?.shape;
    let travel: Vector = along * (length + reach[0] + reach[1]);
    let wedge = ogeom_algo::make_prism(model, &section, travel, tol)?.shape;
    let wedge = explore(model, &wedge, Filter::OfType(ShapeType::Solid))?
        .into_iter()
        .next()
        .unwrap_or(wedge);
    let mut result = if convex {
        ogeom_bool::cut(model, solid, &wedge, tol)?
    } else {
        ogeom_bool::fuse(model, solid, &wedge, tol)?
    };
    result.history.delete(edge);
    Ok(Some(result))
}
