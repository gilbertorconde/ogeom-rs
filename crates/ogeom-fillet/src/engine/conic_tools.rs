//! Tools swept from conic sections on straight, revolved and marched seats.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, Curve3d as _, PlaneSurface, RevolutionSurface, Surface as _,
    SurfaceGeometry,
};
use ogeom_math::{Axis, Circle, Direction, Frame, Plane, Transform};
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
    let mut joints = Vec::new();
    for (i, a) in asked.iter().enumerate() {
        let Some((a0, a1)) = ogeom_algo::edge_vertices(model, &a.edge)? else {
            continue;
        };
        *count.entry(a0.node()).or_insert(0) += 1;
        if !a0.is_same(&a1) {
            *count.entry(a1.node()).or_insert(0) += 1;
        }
        for (j, b) in asked.iter().enumerate().skip(i + 1) {
            let Some((b0, b1)) = ogeom_algo::edge_vertices(model, &b.edge)? else {
                continue;
            };
            for v in [&a0, &a1] {
                if !v.is_same(&b0) && !v.is_same(&b1) {
                    continue;
                }
                if !crate::asym::continues(model, a, b, v, tol)? {
                    ogeom_bail!(
                        Construction,
                        "two asymmetric rounds meet at a vertex other than as one tangent chain naming one face with the same distances and section, and the corner where they meet is not closed yet; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                    );
                }
                joints.push((i, j));
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
        .enumerate()
        .map(|(i, a)| {
            tool(
                &mut scratch,
                solid,
                a,
                joints.iter().any(|(x, y, ..)| *x == i || *y == i),
                tol,
            )
        })
        .collect::<OgeomResult<_>>()?;
    let edges: Vec<_> = asked.iter().map(|a| a.edge.clone()).collect();
    let Some(built) = super::planar_tools::applied(&mut scratch, solid, &edges, &tools, tol)?
    else {
        return Ok(None);
    };
    *model = scratch;
    Ok(Some(built))
}

pub(super) fn tool(
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
        let (seat, first) = crate::asym::straight_seat(model, solid, &a.edge, &a.face, tol)?;
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
        return extended_sweep(
            model,
            solid,
            &a.edge,
            &seat,
            &RoundSpec::Conic {
                d,
                profile: a.profile,
            },
            [!grouped; 2],
            tol,
        );
    }
    let weight = match a.profile {
        RoundProfile::Elliptical => crate::marched::ConicWeight::Elliptical,
        RoundProfile::Conic(_) => crate::marched::ConicWeight::Fixed(a.profile.weight(0.0)?),
    };
    if let Some(tool) = ruled(model, solid, a, weight, tol)? {
        return Ok(tool);
    }
    marched(model, solid, a, None, None, weight, tol)
}

/// The tool of `a`'s band on a seat only the march speaks, its legs `a`'s
/// distances or, with `laws`, a law for the named face and for the other;
/// stopped flush where `mates` continue its crease.
pub(super) fn marched(
    model: &mut Model,
    solid: &Shape,
    a: &Asymmetric,
    laws: Option<&[crate::law::Profile; 2]>,
    mates: Option<(usize, &[crate::fillet::Mate])>,
    weight: crate::marched::ConicWeight,
    tol: Tolerances,
) -> OgeomResult<Tool> {
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
            laws,
            &[],
            mates,
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

/// An asymmetric round whose distances follow `laws` (the named face's
/// and the other's) on a seat only the march speaks: a rim, a curved
/// crease, open or closed. The band is marched under the laws and applied
/// as a one-boolean tool; `None` where the tool does not close on the
/// solid.
///
/// # Errors
///
/// As the marched band, and where a contact runs off its face.
pub(crate) fn varying_marched(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    laws: &[crate::law::Profile; 2],
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let weight = match profile {
        RoundProfile::Elliptical => crate::marched::ConicWeight::Elliptical,
        RoundProfile::Conic(_) => crate::marched::ConicWeight::Fixed(profile.weight(0.0)?),
    };
    let asked = Asymmetric {
        edge: edge.clone(),
        face: face.clone(),
        on_face: laws[0].largest(),
        on_other: laws[1].largest(),
        profile,
    };
    let mut scratch = model.clone();
    let tool = marched(&mut scratch, solid, &asked, Some(laws), None, weight, tol)?;
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

/// The tool of `spec` swept over the straight `seat` of `edge`. With
/// `extend` at an end (the seat's start, then its end), a convex seat's
/// band runs past it where planes following the crease there lean back
/// from it, far enough to clear them, and is trimmed by those overhanging
/// it: the circular round's own ends. A concave seat's fill ends in the
/// faces following its crease there, planes, cylinders or spheres: where
/// one is curved or a plane not square to the crease, the fill runs on
/// until it clears them all and is trimmed by each, kept on the side of
/// the material its face bounds. A varying law is held level over the run
/// past an end.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a fill's end meets a face that is not a plane, a cylinder or a
/// sphere, a curved face the crease runs into tangent to it, or a sphere
/// that does not span the fill's end.
pub(super) fn extended_sweep(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    seat: &crate::support::Seat,
    spec: &RoundSpec,
    extend: [bool; 2],
    tol: Tolerances,
) -> OgeomResult<Tool> {
    let mut seat = seat.clone();
    let mut trims = Vec::new();
    let mut past = [0.0_f64; 2];
    if extend.contains(&true) {
        let original = seat.clone();
        let section = PlanarSection::new(&original, spec, tol)?;
        for (end, outward, at) in [(0, -seat.along, seat.start), (1, seat.along, seat.end)] {
            if !extend[end] {
                continue;
            }
            let hosts = [&original.faces[0], &original.faces[1]];
            // A fill also ends in a cylinder following its crease; a band
            // cut away runs on only to planes.
            let ends = if original.convex {
                crate::support::following_planes(model, solid, edge, hosts, at)?.map(|planes| {
                    planes
                        .into_iter()
                        .map(|(face, origin, normal)| {
                            (face, crate::support::EndFace::Plane(origin, normal))
                        })
                        .collect::<Vec<_>>()
                })
            } else {
                crate::support::following_ends(model, solid, edge, hosts, at, tol)?
            };
            let Some(ends) = ends else {
                if !original.convex
                    && crate::support::following_faces(model, solid, edge, hosts, at)?
                        .is_some_and(|faces| !faces.is_empty())
                {
                    ogeom_bail!(
                        Construction,
                        "a round's fill on a concave edge ends in the faces at its edge's \
                         ends, and an end of this edge meets a face that is not a plane, a \
                         cylinder or a sphere; see docs/CAPABILITIES.md, \
                         fillet.asymmetric-rounds"
                    );
                }
                continue;
            };
            // A cylinder ends the fill only where it spans the fill's end:
            // each point of the end section meets it, run along the crease,
            // on the face itself, read at the vertex's height along the
            // axis. Elsewhere a constant fill stops square and a varying one
            // is refused.
            let hull_at = |s: f64| {
                let mut hull = section.conic_net_at(s).map_or_else(
                    || vec![section.contact(0, s), section.contact(1, s)],
                    |(points, _, _)| points,
                );
                hull.push(at);
                hull
            };
            let mut spans = true;
            for (face, end_face) in &ends {
                let crate::support::EndFace::Drum(axis, ..) = *end_face else {
                    continue;
                };
                let up = axis.direction.vector();
                for p in hull_at(if end == 0 { 0.0 } else { section.length() }) {
                    let Some(t) = end_face.exit(p, outward) else {
                        spans = false;
                        break;
                    };
                    let c = p + outward * t;
                    let c = c + up * (at - c).dot(up);
                    let nudge = up * (tol.confusion() * 10.0);
                    let mut on = false;
                    for q in [c, c + nudge, c - nudge] {
                        on = on || crate::support::on_face_itself(model, face, q, tol)?;
                    }
                    spans = spans && on;
                }
            }
            // A sphere ends the fill only where its continuation past the
            // face is read from the hosts the face meets: each point of the
            // end section, run along the crease to the sphere, and pulled to
            // the sphere's nearest point in the plane of each host sharing
            // an edge with the face, lands on the face itself. Elsewhere the
            // fill is refused, constant or varying.
            for (face, end_face) in &ends {
                let crate::support::EndFace::Ball(centre, radius, _) = *end_face else {
                    continue;
                };
                let mut planes = Vec::new();
                for (host, normal) in hosts.iter().zip(original.normals) {
                    if shares_edge(model, face, host)? {
                        planes.push(normal);
                    }
                }
                let mut held = !planes.is_empty();
                for p in hull_at(if end == 0 { 0.0 } else { section.length() }) {
                    let Some(t) = end_face.exit(p, outward) else {
                        held = false;
                        break;
                    };
                    let c = p + outward * t;
                    for normal in &planes {
                        let flat = c - *normal * (c - at).dot(*normal);
                        let hub = centre - *normal * (centre - at).dot(*normal);
                        let off = flat - hub;
                        let ring = radius.mul_add(radius, -(centre - hub).dot(centre - hub));
                        if ring <= 0.0 || off.magnitude() <= tol.confusion() {
                            held = false;
                            continue;
                        }
                        let q = hub + off * (ring.sqrt() / off.magnitude());
                        let nudge = *normal * (tol.confusion() * 10.0);
                        let mut on = false;
                        for q in [q, q + nudge, q - nudge] {
                            on = on || crate::support::on_face_itself(model, face, q, tol)?;
                        }
                        held = held && on;
                    }
                }
                if !held {
                    ogeom_bail!(
                        Construction,
                        "a round's fill on a concave edge ends in the faces at its edge's \
                         ends, and a sphere at an end of this edge does not span the fill's \
                         end; see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                    );
                }
            }
            if !spans {
                if matches!(spec, RoundSpec::Varying { .. }) {
                    ogeom_bail!(
                        Construction,
                        "an asymmetric round whose distances vary on a concave edge ends in \
                         the faces at its edge's ends, and a cylinder at an end of this edge \
                         does not span the fill's end; see docs/CAPABILITIES.md, \
                         fillet.asymmetric-rounds"
                    );
                }
                continue;
            }
            // A fill capped square already ends in planes square to it or
            // running along it.
            let square = |end: &crate::support::EndFace| match end {
                crate::support::EndFace::Plane(_, normal) => {
                    let rate = normal.dot(outward).abs();
                    rate <= tol.angular() || rate >= 1.0 - tol.angular()
                }
                crate::support::EndFace::Drum(..) | crate::support::EndFace::Ball(..) => false,
            };
            if !original.convex && ends.iter().all(|(_, end)| square(end)) {
                continue;
            }
            let mut reach = 0.0_f64;
            let s = if end == 0 { 0.0 } else { section.length() };
            let hull = hull_at(s);
            for (_, face) in ends {
                match face {
                    crate::support::EndFace::Plane(origin, normal) => {
                        let rate = normal.dot(outward);
                        if rate.abs() <= tol.angular() {
                            continue;
                        }
                        for p in &hull {
                            reach = reach.max((origin - *p).dot(normal) / rate);
                        }
                        if !original.convex {
                            // The fill stands only behind the face, in front
                            // of which there is no material for it to fill
                            // against.
                            trims.push(Trim::Plane(origin, -normal));
                        } else if rate < -tol.angular() {
                            trims.push(Trim::Plane(origin, normal));
                        }
                    }
                    crate::support::EndFace::Drum(..) | crate::support::EndFace::Ball(..) => {
                        // A curved face the crease runs into along its
                        // surface gives no end the fill could stop in.
                        if face.normal_at(at).dot(outward).abs() <= tol.angular() {
                            ogeom_bail!(
                                Construction,
                                "a round's fill on a concave edge ends in the face at its \
                                 end, and the edge runs into that curved face tangent to it; \
                                 see docs/CAPABILITIES.md, fillet.asymmetric-rounds"
                            );
                        }
                        for p in &hull {
                            if let Some(t) = face.exit(*p, outward) {
                                reach = reach.max(t);
                            }
                        }
                        trims.push(Trim::Curved(face));
                    }
                }
            }
            if reach > tol.confusion() {
                let [d0, d1] = section.distances(s);
                reach += d0.max(d1) * 0.25;
            }
            past[end] = reach;
            if end == 0 {
                seat.start += outward * reach;
            } else {
                seat.end += outward * reach;
            }
        }
    }
    let spec = match spec {
        RoundSpec::Varying { laws, profile } if past.iter().any(|r| *r > 0.0) => {
            let length = (seat.end - seat.start).magnitude() - past[0] - past[1];
            RoundSpec::Varying {
                laws: laws
                    .clone()
                    .map(|law| law.stretched(past[0] / length, past[1] / length)),
                profile: *profile,
            }
        }
        other => other.clone(),
    };
    let section = PlanarSection::new(&seat, &spec, tol)?;
    let mut tool = super::planar_tools::sweep(model, edge, &section, tol)?;
    for trim in trims {
        let span = ogeom_algo::shape_bounds(model, &tool.solid, tol)?.diagonal() * 2.0;
        let solid = match trim {
            Trim::Plane(origin, normal) => {
                crate::support::trimmed_wedge(model, &tool.solid, (origin, normal), span, tol)?
            }
            Trim::Curved(face @ crate::support::EndFace::Ball(..)) => {
                crate::support::trimmed_by_ball(model, &tool.solid, face, tol)?
            }
            Trim::Curved(face) => {
                crate::support::trimmed_by_drum(model, &tool.solid, face, span, tol)?
            }
        };
        tool = named(model, solid, edge, tool.sense)?;
    }
    Ok(tool)
}

/// What cuts a tool run on past an end of its seat.
enum Trim {
    /// The half-space behind this plane (a point on it and the normal
    /// pointing away from what is kept) is cut away.
    Plane(ogeom_math::Point, ogeom_math::Vector),
    /// Only what lies on the material side of this cylinder or sphere is
    /// kept.
    Curved(crate::support::EndFace),
}

/// Whether faces `a` and `b` share an edge.
fn shares_edge(model: &Model, a: &Shape, b: &Shape) -> OgeomResult<bool> {
    use ogeom_topo::{Filter, explore};
    let theirs = explore(model, b, Filter::OfType(ogeom_topo::ShapeType::Edge))?;
    Ok(
        explore(model, a, Filter::OfType(ogeom_topo::ShapeType::Edge))?
            .iter()
            .any(|e| theirs.iter().any(|f| f.is_same(e))),
    )
}

/// How a host of a ruled seat meets the plane square to its edge: a line,
/// or a circle about the foot of its cylinder's axis.
#[derive(Clone, Copy)]
enum Ruled {
    /// A plane, of this normal.
    Line(ogeom_math::Vector),
    /// A cylinder of this axis and radius.
    Circle(Axis, f64),
}

impl Ruled {
    /// The point `distance` along the host's section from the crease
    /// point `at`, leaving it along `into`, and the host's normal there.
    fn walk(
        self,
        at: ogeom_math::Point,
        into: ogeom_math::Vector,
        distance: f64,
    ) -> (ogeom_math::Point, ogeom_math::Vector) {
        match self {
            Self::Line(normal) => (at + into * distance, normal),
            Self::Circle(axis, radius) => {
                let up = axis.direction.vector();
                let foot = axis.location + up * (at - axis.location).dot(up);
                let turn = distance / radius;
                let p = foot + (at - foot) * turn.cos() + into * (radius * turn.sin());
                (p, (p - foot) / radius)
            }
        }
    }
}

/// The section of an asymmetric round on a ruled seat, as
/// [`ruled_section`] reads it.
struct RuledSection {
    start: ogeom_math::Point,
    along: ogeom_math::Vector,
    length: f64,
    sections: [Ruled; 2],
    into: [ogeom_math::Vector; 2],
    d: [f64; 2],
    control: [ogeom_math::Point; 3],
    w: f64,
    convex: bool,
}

/// The section of an asymmetric round on a straight edge where two faces,
/// a plane or a cylinder each and at least one a cylinder, meet along the
/// cylinders' ruling: the same at every station. Each leg runs its
/// distance along its face's section, the arc on a cylinder, as the
/// marched band walks it, and the conic stands on the shoulder the faces'
/// tangents meet on at the contacts. `None` for any other seat, or where
/// an end of the edge meets a face not square to it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a contact sets back past its face's side, the faces' sections
/// leave the contacts along tangents that do not meet in front of the
/// crease, or the conic crosses a cylinder it is tangent to between the
/// contacts.
fn ruled_section(
    model: &Model,
    solid: &Shape,
    a: &Asymmetric,
    weight: crate::marched::ConicWeight,
    tol: Tolerances,
) -> OgeomResult<Option<RuledSection>> {
    use ogeom_geom::Transformable as _;
    use ogeom_topo::{Filter, explore};
    let conf = tol.confusion();
    let (Curve::Line(_), _) = crate::support::edge_curve(model, &a.edge, tol)? else {
        return Ok(None);
    };
    let Some((v0, v1)) = ogeom_algo::edge_vertices(model, &a.edge)? else {
        return Ok(None);
    };
    let point = |v: &Shape| crate::setback::point_of(model, v);
    let (start, end) = (point(&v0)?, point(&v1)?);
    let length = start.distance(end);
    if length <= conf {
        return Ok(None);
    }
    let along = (end - start) / length;
    let mut hosts: Vec<(Shape, Ruled)> = Vec::new();
    let mut outward: Vec<ogeom_math::Vector> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let holds = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| crate::support::same_occurrence(model, e, &a.edge, tol));
        if !holds {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            return Ok(None);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        let placed = surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?;
        let section = match placed {
            SurfaceGeometry::Plane(p) => Ruled::Line(p.plane().normal().vector()),
            SurfaceGeometry::Cylinder(c) => {
                let axis = c.cylinder().axis();
                let up = axis.direction.vector();
                let off = start - axis.location;
                let radial = off - up * off.dot(up);
                if up.cross(along).magnitude() > tol.angular()
                    || (radial.magnitude() - c.cylinder().radius()).abs() > conf
                {
                    return Ok(None);
                }
                Ruled::Circle(axis, c.cylinder().radius())
            }
            _ => return Ok(None),
        };
        // The face's outward normal at the edge's start.
        let sign = if face.orientation() == ogeom_topo::Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        outward.push(section.walk(start, along, 0.0).1 * sign);
        hosts.push((face, section));
    }
    let Ok([(f0, s0), (f1, s1)]) = <[(Shape, Ruled); 2]>::try_from(hosts) else {
        return Ok(None);
    };
    if matches!((s0, s1), (Ruled::Line(_), Ruled::Line(_))) {
        return Ok(None);
    }
    // Each end square: the band is capped in the end faces' planes.
    for v in [&v0, &v1] {
        for face in explore_unique(model, solid, ShapeType::Face)? {
            if face.is_same(&f0) || face.is_same(&f1) {
                continue;
            }
            let meets = explore_unique(model, &face, ShapeType::Vertex)?
                .iter()
                .any(|w| crate::support::same_occurrence(model, w, v, tol));
            if !meets {
                continue;
            }
            let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                return Ok(None);
            };
            let square = match model
                .geometry()
                .surface(data.surface)
                .map(|s| s.clone().transformed(&face.transform(model.datums())?, tol))
            {
                Some(Ok(SurfaceGeometry::Plane(p))) => {
                    p.plane().normal().vector().cross(along).magnitude() <= tol.angular()
                }
                _ => false,
            };
            if !square {
                return Ok(None);
            }
        }
    }
    let faces = [f0, f1];
    let sections = [s0, s1];
    let d = if faces[0].is_same(&a.face) {
        [a.on_face, a.on_other]
    } else if faces[1].is_same(&a.face) {
        [a.on_other, a.on_face]
    } else {
        ogeom_bail!(
            Construction,
            "the face an asymmetric round's first distance runs along is not one of the \
             edge's two faces"
        );
    };
    // Into each face from the middle of the edge: square to the edge on
    // the face, the side a short step lands on it.
    let middle = start + along * (length * 0.5);
    let mut into = [ogeom_math::Vector::ZERO; 2];
    for i in 0..2 {
        let normal = sections[i].walk(middle, along, 0.0).1;
        let heading = along.cross(normal);
        let heading = heading / heading.magnitude();
        let step = (d[i] * 1e-2).max(conf * 100.0);
        let mut found = None;
        for sign in [1.0, -1.0] {
            let (p, _) = sections[i].walk(middle, heading * sign, step);
            if crate::support::on_face_or_boundary(model, solid, &faces[i], p, tol)? {
                if found.is_some() {
                    return Ok(None);
                }
                found = Some(heading * sign);
            }
        }
        let Some(found) = found else {
            return Ok(None);
        };
        into[i] = found;
    }
    // The legs' turn on a cylinder stays short of half way round.
    for i in 0..2 {
        if let Ruled::Circle(_, radius) = sections[i]
            && d[i] >= radius * core::f64::consts::PI * 0.5
        {
            return Ok(None);
        }
    }
    let walked = [
        sections[0].walk(start, into[0], d[0]),
        sections[1].walk(start, into[1], d[1]),
    ];
    for (i, (p, _)) in walked.iter().enumerate() {
        for fraction in [0.25, 0.5, 0.75] {
            let q = *p + along * (length * fraction);
            if !crate::support::on_face_or_boundary(model, solid, &faces[i], q, tol)? {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round sets back {} along a face of the edge through \
                     {q:?}, past that face's side; see docs/CAPABILITIES.md, \
                     fillet.asymmetric-rounds",
                    d[i]
                );
            }
        }
    }
    let Some((control, w)) = crate::marched::conic_section(
        [walked[0].0, walked[1].0],
        [walked[0].1, walked[1].1],
        start,
        weight,
        tol,
    ) else {
        ogeom_bail!(
            Construction,
            "an asymmetric round's faces leave their contacts at {start:?} along tangents \
             that do not meet in front of the crease; see docs/CAPABILITIES.md, \
             fillet.asymmetric-rounds"
        );
    };
    // Each cylinder's circle is touched at its contact; the conic must not
    // cross it between the contacts, or the band leaves the face it is
    // tangent to and no round of that section fits the seat.
    for (i, section) in sections.iter().enumerate() {
        if let Ruled::Circle(axis, radius) = *section
            && !keeps_to_circle(control, w, axis, radius, i)
        {
            ogeom_bail!(
                Construction,
                "an asymmetric round's section crosses the drum it is tangent to between its \
                 contacts at {start:?}: the conic standing on its faces' shoulder bulges past \
                 the drum, so no band of that section fits there; see docs/CAPABILITIES.md, \
                 fillet.asymmetric-rounds"
            );
        }
    }
    // The material side: at a convex crease a leg leaves its face's
    // crease under the other face's tangent plane.
    let convex = into[0].dot(outward[1]) < 0.0;
    Ok(Some(RuledSection {
        start,
        along,
        length,
        sections,
        into,
        d,
        control,
        w,
        convex,
    }))
}

/// Whether the rational quadratic arc about `p`, of middle weight `w`,
/// keeps to one side of the circle a cylinder of `axis` and `radius` cuts
/// in the arc's plane, square to the axis, touching it at its end `at`
/// (0 for `p[0]`, 1 for `p[2]`): the side its other end stands on.
///
/// With the arc's homogeneous point `N(t)` over its weight `S(t)`, each a
/// quadratic in Bernstein form, `G = |N − c S|² − r² S²` is a quartic of
/// the condition's sign, `c` the axis's foot. The touch is a double root
/// at `at`, divided out in closed form, which leaves a quadratic whose
/// least value over `[0, 1]` is read exactly: the arc crosses the circle
/// between its ends where that changes sign, however shallowly.
fn keeps_to_circle(p: [ogeom_math::Point; 3], w: f64, axis: Axis, radius: f64, at: usize) -> bool {
    let up = axis.direction.vector();
    let weights = [1.0, w, 1.0];
    let q: [ogeom_math::Vector; 3] = core::array::from_fn(|k| {
        let off = p[k] - axis.location;
        (off - up * off.dot(up)) * weights[k]
    });
    // The Bernstein coefficients of the product of two quadratics.
    let product = |a: [f64; 3], b: [f64; 3]| -> [f64; 5] {
        [
            a[0] * b[0],
            0.5 * a[0].mul_add(b[1], a[1] * b[0]),
            (4.0 * a[1]).mul_add(b[1], a[0].mul_add(b[2], a[2] * b[0])) / 6.0,
            0.5 * a[1].mul_add(b[2], a[2] * b[1]),
            a[2] * b[2],
        ]
    };
    let (x, y, z) = (q.map(|v| v.x), q.map(|v| v.y), q.map(|v| v.z));
    let (xx, yy, zz, ww) = (
        product(x, x),
        product(y, y),
        product(z, z),
        product(weights, weights),
    );
    let g: [f64; 5] =
        core::array::from_fn(|k| (-radius * radius).mul_add(ww[k], xx[k] + yy[k] + zz[k]));
    // `G = t² H` touching at the start, `G = (1 − t)² H` at the end; the
    // other end's value is `H`'s there.
    let (h, far) = if at == 0 {
        ([6.0 * g[2], 2.0 * g[3], g[4]], g[4])
    } else {
        ([g[0], 2.0 * g[1], 6.0 * g[2]], g[0])
    };
    // Each coefficient is a blend of products of coordinates rounded to
    // their size and the radius's; the division multiplies by six.
    let size = q
        .iter()
        .zip(weights)
        .map(|(v, wk)| (v.magnitude() + radius * wk).powi(2))
        .fold(0.0, f64::max);
    let noise = size * f64::EPSILON * 64.0 * 6.0;
    if far.abs() <= noise {
        return false;
    }
    let h = h.map(|c| c * far.signum());
    if h[0] < -noise || h[2] < -noise {
        return false;
    }
    if h[1] >= 0.0 {
        return true;
    }
    // Both ends above, the middle coefficient below: the least value is
    // inside, at `(h0 − h1) / (h0 − 2 h1 + h2)`.
    let bend = 2.0f64.mul_add(-h[1], h[0] + h[2]);
    h[0].mul_add(h[2], -h[1] * h[1]) / bend >= -noise
}

/// Whether the asymmetric round `a`, with the elliptical section, fits its
/// ruled seat by [`ruled_section`]: `None` where the seat is not ruled.
///
/// # Errors
///
/// As reading the solid; a refusal of the section is `Some(false)`.
pub(crate) fn ruled_fits(
    model: &Model,
    solid: &Shape,
    a: &Asymmetric,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    match ruled_section(
        model,
        solid,
        a,
        crate::marched::ConicWeight::Elliptical,
        tol,
    ) {
        Ok(found) => Ok(found.map(|_| true)),
        Err(ogeom_core::OgeomError::Construction(_)) => Ok(Some(false)),
        Err(e) => Err(e),
    }
}

/// The tool of an asymmetric round on a ruled seat ([`ruled_section`]):
/// the band is the section extruded along the edge, exact. At an end whose
/// plane bounds the whole solid the prism runs on past it, so the tool
/// shares no face with the solid's end there: the end cuts the band in its
/// own exact section, rather than the band's end edge being matched in
/// one plane against the solid's edge it runs beside. `None` where the
/// seat is not ruled.
///
/// # Errors
///
/// As [`ruled_section`].
fn ruled(
    model: &mut Model,
    solid: &Shape,
    a: &Asymmetric,
    weight: crate::marched::ConicWeight,
    tol: Tolerances,
) -> OgeomResult<Option<Tool>> {
    let Some(found) = ruled_section(model, solid, a, weight, tol)? else {
        return Ok(None);
    };
    let RuledSection {
        start,
        along,
        length,
        sections,
        into,
        d,
        control,
        w,
        convex,
    } = found;
    // How far past each end a cut's prism runs: a leg's length where
    // nothing of the solid stands beyond that end's plane, read on its
    // guaranteed bound (widened by its vertices' and edges' tolerances, a
    // hundred confusion distances at most), else none. A fill runs on
    // nowhere: past the solid it would add material.
    let bound = ogeom_algo::shape_bounds(model, solid, tol)?;
    let reach = 100.0 * tol.confusion();
    let past = |at: ogeom_math::Point, outward: ogeom_math::Vector| -> f64 {
        let beyond = bound
            .corners()
            .iter()
            .map(|c| (*c - at).dot(outward))
            .fold(f64::NEG_INFINITY, f64::max);
        if convex && beyond <= reach {
            d[0].max(d[1])
        } else {
            0.0
        }
    };
    let end = start + along * length;
    let (before, after) = (past(start, -along), past(end, along));
    let shift = -along * before;
    let start = start + shift;
    let length = length + before + after;
    let control = control.map(|p| p + shift);
    let vertex = |model: &mut Model, p| ogeom_algo::make_vertex(model, p).shape;
    let (vc, va, vb) = (
        vertex(model, start),
        vertex(model, control[0]),
        vertex(model, control[2]),
    );
    let leg = |model: &mut Model, i: usize, to: &Shape| -> OgeomResult<Shape> {
        let (curve, range): (Curve, (f64, f64)) = match sections[i] {
            Ruled::Line(_) => (
                ogeom_geom::LineCurve::new(Axis::new(start, Direction::new(into[i], tol)?)).into(),
                (0.0, d[i]),
            ),
            Ruled::Circle(axis, radius) => {
                let up = axis.direction.vector();
                let foot = axis.location + up * (start - axis.location).dot(up);
                let x = (start - foot) / radius;
                (
                    Curve::Circle(CircleCurve::new(Circle::new(
                        Frame::new(
                            foot,
                            Direction::new(x.cross(into[i]), tol)?,
                            Direction::new(x, tol)?,
                            tol,
                        )?,
                        radius,
                        tol,
                    )?)),
                    (0.0, d[i] / radius),
                )
            }
        };
        Ok(ogeom_algo::make_edge_between(model, curve, range, &vc, to, tol)?.shape)
    };
    let first = leg(model, 0, &va)?;
    let second = leg(model, 1, &vb)?;
    let arc = crate::asym::conic_arc(control, w, tol)?;
    let arc = ogeom_algo::make_edge_between(model, arc, (0.0, 1.0), &va, &vb, tol)?.shape;
    let wire = ogeom_algo::make_wire(model, &[first, arc, second.reversed()], tol)?.shape;
    let plane: SurfaceGeometry =
        PlaneSurface::new(Plane::through(start, Direction::new(along, tol)?)).into();
    let profile = ogeom_algo::make_face(model, plane, &[wire], tol)?.shape;
    let prism = ogeom_algo::make_prism(model, &profile, along * length, tol)?.shape;
    Ok(Some(named(
        model,
        prism,
        &a.edge,
        if convex { Sense::Cut } else { Sense::Fuse },
    )?))
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
    // Anchor periodic charts at the arc's start, clear of the circle's seam.
    let circle = |radius, height| -> OgeomResult<(Curve, (f64, f64))> {
        Ok((
            Curve::Circle(CircleCurve::new(Circle::new(
                Frame::new(
                    seat.centre + seat.up * height,
                    Direction::new(seat.up, tol)?,
                    Direction::new(rotate(a).apply_vector(seat.x_ref.vector()), tol)?,
                    tol,
                )?,
                radius,
                tol,
            )?)),
            (0.0, b - a),
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
    // Contact levels trim the wall, not its artificial chart boundary.
    let reach = rim.size() + 1.0;
    let wall = seat.wall_surface(
        (-reach, reach),
        (contact.wall_height.min(0.0), contact.wall_height.max(0.0)),
        Direction::new(rotate(a).apply_vector(seat.x_ref.vector()), tol)?,
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
