//! Tools swept from conic sections on straight, revolved and marched seats.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, Curve3d as _, PlaneSurface, RevolutionSurface, Surface as _,
    SurfaceGeometry,
};
use ogeom_math::{Axis, Circle, Direction, Plane, Transform};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{SheetEdge, SheetFace, build_faces};
use super::fallback::{Sense, Tool};
use crate::asym::{Asymmetric, RimConic, RoundProfile};
use crate::round::section::{PlanarSection, RoundSpec};

#[allow(
    clippy::too_many_arguments,
    reason = "the edge's two distances and section"
)]
pub(crate) fn distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut scratch = model.clone();
    let tool = tool(
        &mut scratch,
        solid,
        &Asymmetric {
            edge: edge.clone(),
            face: face.clone(),
            on_face,
            on_other,
            profile,
        },
        false,
        tol,
    )?;
    let Some(built) = super::planar_tools::applied(
        &mut scratch,
        solid,
        std::slice::from_ref(edge),
        &[tool],
        tol,
    )?
    else {
        return Ok(None);
    };
    *model = scratch;
    Ok(Some(built))
}

pub(crate) fn asymmetric(
    model: &mut Model,
    solid: &Shape,
    asked: &[Asymmetric],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut asked = asked.to_vec();
    asked.sort_by_key(|a| a.edge.node());
    let mut count = std::collections::BTreeMap::new();
    for (i, a) in asked.iter().enumerate() {
        let Some((a0, a1)) = ogeom_algo::edge_vertices(model, &a.edge)? else {
            continue;
        };
        *count.entry(a0.node()).or_insert(0) += 1;
        if !a0.is_same(&a1) {
            *count.entry(a1.node()).or_insert(0) += 1;
        }
        for b in &asked[i + 1..] {
            let Some((b0, b1)) = ogeom_algo::edge_vertices(model, &b.edge)? else {
                continue;
            };
            for v in [&a0, &a1] {
                if (v.is_same(&b0) || v.is_same(&b1))
                    && !crate::asym::continues(model, a, b, v, tol)?
                {
                    ogeom_bail!(
                        Construction,
                        "two asymmetric rounds meet at a vertex other than as one tangent chain naming one face with the same distances and section, and the corner where they meet is not closed yet; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                    );
                }
            }
        }
    }
    if count.values().any(|n| *n > 2) {
        ogeom_bail!(
            Construction,
            "three asymmetric rounds meet at a vertex, and the corner where they meet is not closed yet; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
        );
    }
    let mut scratch = model.clone();
    let tools: Vec<_> = asked
        .iter()
        .map(|a| tool(&mut scratch, solid, a, true, tol))
        .collect::<OgeomResult<_>>()?;
    let edges: Vec<_> = asked.iter().map(|a| a.edge.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    *model = scratch;
    Ok(Some(built))
}

fn tool(
    model: &mut Model,
    solid: &Shape,
    a: &Asymmetric,
    grouped: bool,
    tol: Tolerances,
) -> OgeomResult<Tool> {
    if let Some(rim) = crate::asym::rim_conic(
        model, solid, &a.edge, &a.face, a.on_face, a.on_other, a.profile, tol,
    )? {
        crate::asym::rim_fits(model, solid, &a.edge, &rim, tol)?;
        return revolved(model, &a.edge, &rim, tol);
    }
    if crate::support::planar_seat(model, solid, &a.edge, tol).is_ok() {
        let (mut seat, first) = crate::asym::straight_seat(model, solid, &a.edge, &a.face, tol)?;
        let d = if first == 0 {
            [a.on_face, a.on_other]
        } else {
            [a.on_other, a.on_face]
        };
        let Some(spec) = crate::asym::straight_section(
            model, solid, &a.edge, &a.face, a.on_face, a.on_other, a.profile, tol,
        )?
        else {
            ogeom_bail!(
                Construction,
                "the circular round asked as an asymmetric one"
            );
        };
        crate::asym::contacts_fit(model, solid, &a.edge, spec, tol)?;
        let mut trims = Vec::new();
        if seat.convex && !grouped {
            let original = seat.clone();
            let section = PlanarSection::new(
                &original,
                &RoundSpec::Conic {
                    d,
                    profile: a.profile,
                },
                tol,
            )?;
            for (end, outward, at) in [(0, -seat.along, seat.start), (1, seat.along, seat.end)] {
                let hosts = [&original.faces[0], &original.faces[1]];
                let Some(planes) =
                    crate::support::following_planes(model, solid, &a.edge, hosts, at)?
                else {
                    continue;
                };
                let mut reach = 0.0_f64;
                let mut hull = section
                    .conic_net_at(if end == 0 { 0.0 } else { section.length() })
                    .map_or_else(Vec::new, |(points, _, _)| points);
                hull.push(at);
                for (_, origin, normal) in planes {
                    let rate = normal.dot(outward);
                    if rate.abs() <= tol.angular() {
                        continue;
                    }
                    for p in &hull {
                        reach = reach.max((origin - *p).dot(normal) / rate);
                    }
                    if rate < -tol.angular() {
                        trims.push((origin, normal));
                    }
                }
                if reach > tol.confusion() {
                    reach += a.on_face.max(a.on_other) * 0.25;
                }
                if end == 0 {
                    seat.start += outward * reach;
                } else {
                    seat.end += outward * reach;
                }
            }
        }
        let section = PlanarSection::new(
            &seat,
            &RoundSpec::Conic {
                d,
                profile: a.profile,
            },
            tol,
        )?;
        let mut tool = super::planar_tools::sweep(model, &a.edge, &section, tol)?;
        for plane in trims {
            let span = ogeom_algo::shape_bounds(model, &tool.solid, tol)?.diagonal() * 2.0;
            let solid = crate::support::trimmed_wedge(model, &tool.solid, plane, span, tol)?;
            tool = named(model, solid, &a.edge, tool.sense)?;
        }
        return Ok(tool);
    }
    let weight = match a.profile {
        RoundProfile::Elliptical => crate::marched::ConicWeight::Elliptical,
        RoundProfile::Conic(_) => crate::marched::ConicWeight::Fixed(a.profile.weight(0.0)?),
    };
    let spec = crate::Chamfer::Distances {
        face: a.face.clone(),
        on_face: a.on_face,
        on_other: a.on_other,
    };
    let (built, wedges) = crate::support::collecting_wedges(|| {
        crate::chamfer_marched::marched_bevel_with(
            model,
            solid,
            &a.edge,
            &spec,
            &[],
            crate::marched::Chord::Conic(weight),
            tol,
        )
    });
    built?;
    let [wedge] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|_| {
        ogeom_core::ogeom_err!(
            Construction,
            "a conic's marched strip does not make one tool"
        )
    })?;
    named(
        model,
        wedge.solid,
        &a.edge,
        if wedge.additive {
            Sense::Fuse
        } else {
            Sense::Cut
        },
    )
}

fn named(model: &Model, solid: Shape, edge: &Shape, sense: Sense) -> OgeomResult<Tool> {
    let faces = explore_unique(model, &solid, ShapeType::Face)?
        .into_iter()
        .map(|f| (f, edge.clone()))
        .collect();
    Ok(Tool {
        solid,
        faces,
        sense,
    })
}

fn revolved(model: &mut Model, edge: &Shape, rim: &RimConic, tol: Tolerances) -> OgeomResult<Tool> {
    let seat = &rim.seat;
    let axis = Axis {
        location: seat.centre,
        direction: Direction::new(seat.up, tol)?,
    };
    let sense = if seat.additive() {
        Sense::Fuse
    } else {
        Sense::Cut
    };
    if rim.whole {
        let flanks =
            crate::support::revolved_flanks(model, seat, rim.wall_depth, rim.cap_rho, tol)?;
        let arc = crate::asym::conic_arc(rim.points, rim.w, tol)?;
        let middle = arc.point_at(0.5, tol)?;
        let surface: SurfaceGeometry =
            RevolutionSurface::new(arc, axis, core::f64::consts::TAU)?.into();
        let band = ogeom_algo::make_profile_revolution_band(
            model,
            &surface,
            &flanks.wall_ring,
            &flanks.cap_ring,
            tol,
        )?;
        let (du, dv) = surface.d1_at(0.0, 0.5, tol)?;
        let band = if du.cross(dv).dot(middle - rim.points[1]) < 0.0 {
            band.reversed()
        } else {
            band
        };
        let solid = crate::support::wedge_solid(
            model,
            &[flanks.wall_band, flanks.annulus, band],
            false,
            tol,
        )?;
        return named(model, solid, edge, sense);
    }
    let (_, range) = crate::support::edge_curve(model, edge, tol)?;
    let (a, b) = if rim.circle.circle().frame().z().vector().dot(seat.up) > 0.0 {
        range
    } else {
        (-range.1, -range.0)
    };
    let rotate = |t| Transform::rotation(axis, t);
    let points = [
        rotate(a).apply(rim.points[1]),
        rotate(b).apply(rim.points[1]),
        rotate(a).apply(rim.points[0]),
        rotate(b).apply(rim.points[0]),
        rotate(a).apply(rim.points[2]),
        rotate(b).apply(rim.points[2]),
    ];
    let contact = rim.contacts(seat);
    let circle = |radius, height| -> OgeomResult<(Curve, (f64, f64))> {
        Ok((
            Curve::Circle(CircleCurve::new(Circle::new(
                seat.frame_at(seat.centre + seat.up * height, tol)?,
                radius,
                tol,
            )?)),
            (a, b),
        ))
    };
    let mut edges = vec![SheetEdge {
        ends: (0, 1),
        curve: Some(circle(seat.radius, 0.0)?),
    }];
    edges.extend(
        [(0, 2), (1, 3), (0, 4), (1, 5)]
            .into_iter()
            .map(|ends| SheetEdge { ends, curve: None }),
    );
    edges.extend([
        SheetEdge {
            ends: (2, 3),
            curve: Some(circle(contact.wall_rho, contact.wall_height)?),
        },
        SheetEdge {
            ends: (4, 5),
            curve: Some(circle(contact.rho, 0.0)?),
        },
        SheetEdge {
            ends: (2, 4),
            curve: Some((rim.meridian(seat, a, tol)?, (0.0, 1.0))),
        },
        SheetEdge {
            ends: (3, 5),
            curve: Some((rim.meridian(seat, b, tol)?, (0.0, 1.0))),
        },
    ]);
    let band: SurfaceGeometry =
        RevolutionSurface::over(rim.meridian(seat, a, tol)?, axis, (0.0, b - a))?.into();
    let wall = seat.wall_surface(
        (contact.wall_height.min(0.0), contact.wall_height.max(0.0)),
        (contact.wall_height.min(0.0), contact.wall_height.max(0.0)),
        seat.x_ref,
        tol,
    )?;
    let cap: SurfaceGeometry =
        PlaneSurface::new(Plane::new(seat.frame_at(seat.centre, tol)?)).into();
    let rings = [
        vec![(0, true), (2, true), (5, false), (1, false)],
        vec![(0, false), (3, true), (6, true), (4, false)],
        vec![(1, true), (7, true), (3, false)],
        vec![(2, false), (4, true), (8, false)],
        vec![(5, true), (8, true), (6, false), (7, false)],
    ];
    let surfaces = [Some(wall), Some(cap), None, None, Some(band)];
    let specs: Vec<_> = rings
        .into_iter()
        .zip(surfaces)
        .map(|(ring, surface)| SheetFace {
            ring,
            surface,
            holes: Vec::new(),
        })
        .collect();
    let mut faces = build_faces(model, &points, &edges, &specs, tol)?;
    let middle = rim
        .meridian(seat, f64::midpoint(a, b), tol)?
        .point_at(0.5, tol)?;
    let radial = rotate(f64::midpoint(a, b)).apply_vector(seat.x_ref.vector());
    let sign = if seat.additive() { -1.0 } else { 1.0 };
    let outward = [
        (radial - seat.up * seat.slope) * (seat.sigma * sign),
        seat.up * sign,
        -seat.up.cross(rotate(a).apply_vector(seat.x_ref.vector())),
        seat.up.cross(rotate(b).apply_vector(seat.x_ref.vector())),
        middle - rotate(f64::midpoint(a, b)).apply(rim.points[1]),
    ];
    for (face, wanted) in faces.iter_mut().zip(outward) {
        let data = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .ok_or_else(|| ogeom_core::ogeom_err!(Dangling, "a conic tool's face"))?;
        let surface = model
            .geometry()
            .surface(data.surface)
            .ok_or_else(|| ogeom_core::ogeom_err!(Dangling, "a conic tool's surface"))?;
        let p = ogeom_algo::shape_bounds(model, face, tol)?
            .centre()
            .unwrap_or(middle);
        let foot = ogeom_algo::project_on_surface(surface, p, 16, tol)?;
        let (du, dv) = surface.d1_at(foot.parameters.0, foot.parameters.1, tol)?;
        if du.cross(dv).dot(wanted) < 0.0 {
            *face = face.reversed();
        }
    }
    let shell = model.add_shell(&faces)?;
    let solid = model.add_solid(&[shell])?;
    named(model, solid, edge, sense)
}
