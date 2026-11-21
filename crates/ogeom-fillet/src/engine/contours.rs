//! Constant-radius contour bands on planar and swept hosts (ADR 0001,
//! amendment "contour chains on swept walls use engine-owned tools"):
//! every edge of the request rides a planar face, and its other face is a
//! wall swept from a fitted profile, or every edge is an arc. Each band is
//! built on the solid as given with the request's other edges as its
//! mates, the bands are joined at their junctions, and the joined tools
//! are applied together.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fillet::Mate;
use crate::round::section::{PlanarSection, RoundSpec};
use crate::support::Wedge;

/// `asked` rounded as one contour. `None` where the radii differ, an edge
/// is not between a planar face and one other face, no edge has a swept
/// host and not every edge is an arc, a ball rolls onto a far edge, three
/// edges meet at a vertex, or a band, junction or boolean declines.
///
/// # Errors
///
/// As the section builders and the booleans, other than their refusals,
/// which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    rounds_piped(model, solid, asked, false, tol)
}

/// [`rounds`], and where `piped`, each link whose ball changes walls off
/// both its vertices built as the pipe between the two crossings
/// ([`super::pipe_corner`]), its neighbours ending on those sections.
/// Without `piped` such a request is not this route's to build.
///
/// # Errors
///
/// As [`rounds`].
pub(crate) fn rounds_piped(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    piped: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some((_, radius)) = asked.first() else {
        return Ok(None);
    };
    if !radius.is_finite() || *radius <= tol.confusion() || asked.iter().any(|(_, r)| r != radius) {
        return Ok(None);
    }
    let mut swept = false;
    let mut circles = true;
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let face_edges = faces
        .iter()
        .map(|f| explore_unique(model, f, ShapeType::Edge))
        .collect::<OgeomResult<Vec<_>>>()?;
    for (edge, _) in asked {
        let faces: Vec<_> = faces
            .iter()
            .zip(&face_edges)
            .filter(|(_, edges)| {
                edges
                    .iter()
                    .any(|e| crate::support::same_occurrence(model, e, edge, tol))
            })
            .map(|(f, _)| f.clone())
            .collect();
        if faces.len() != 2 {
            return Ok(None);
        }
        let mut planar = false;
        for face in faces {
            let data = model
                .node(&face)
                .and_then(|n| n.data().as_face())
                .ok_or_else(|| ogeom_core::ogeom_err!(Dangling, "a contour host is missing"))?;
            match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Plane(_)) => planar = true,
                Some(SurfaceGeometry::Extrusion(_)) | Some(SurfaceGeometry::BSpline(_)) => {
                    swept = true
                }
                _ => {}
            }
        }
        if !planar {
            return Ok(None);
        }
        circles &= matches!(
            crate::support::edge_curve(model, edge, tol)?.0,
            Curve::Circle(_)
        );
    }
    if !(swept || circles || piped) {
        return Ok(None);
    }
    let mut once = asked.to_vec();
    once.sort_unstable_by_key(|(e, _)| super::id(e));
    let mut scratch = model.clone();
    let built =
        super::marched_seats::building_tools(|| build(&mut scratch, solid, &once, piped, tol));
    match built {
        Ok(Some(built)) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ OgeomError::Cancelled) => Err(e),
        Ok(None) | Err(OgeomError::Construction(_) | OgeomError::NotDone(_)) => Ok(None),
        Err(e) => Err(e),
    }
}

/// [`rounds`] on a scratch model: each band set aside unapplied, joined at
/// its junctions, then applied, the result checked against the solid.
fn build(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    piped: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let edges: Vec<_> = once.iter().map(|(e, _)| e.clone()).collect();
    let radii: Vec<_> = once.iter().map(|(_, r)| *r).collect();
    if super::request_fits(model, solid, &edges, &radii, tol)
        .map_err(|b| b.error)?
        .iter()
        .any(|r| *r)
    {
        return Ok(None);
    }
    let ends = edges
        .iter()
        .map(|edge| {
            ogeom_algo::edge_vertices(model, edge)
                .map(|p| p.map(|(a, b)| vec![a, b]).unwrap_or_default())
        })
        .collect::<OgeomResult<Vec<_>>>()?;
    if ends
        .iter()
        .flatten()
        .any(|v| ends.iter().flatten().filter(|w| w.is_same(v)).count() > 2)
    {
        return Ok(None);
    }
    let Some(mut met) =
        super::marched_junctions::met(model, solid, &edges, &radii, &ends, true, tol)?
    else {
        return Ok(None);
    };
    let pipes = if piped {
        let Some(pipes) = pipes(model, solid, once, tol)? else {
            return Ok(None);
        };
        pipes
    } else {
        Vec::new()
    };
    let mut bands = Vec::new();
    for (index, (edge, radius)) in once.iter().enumerate() {
        if pipes.iter().any(|(link, ..)| *link == index) {
            continue;
        }
        let aside = core::cell::RefCell::new(Vec::new());
        let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
        let rim = if let Curve::Circle(c) = &curve {
            match crate::support::revolved_seat(model, solid, edge, c, tol) {
                Ok(_) => true,
                Err(OgeomError::Construction(_)) => false,
                Err(e) => return Err(e),
            }
        } else {
            false
        };
        match curve {
            Curve::Line(_) if crate::support::hosts_planar(model, solid, edge, tol)? => {
                bands.push((
                    index,
                    straight(model, solid, edge, *radius, (index, &met.mates), tol)?,
                ));
                continue;
            }
            Curve::Circle(ref c) if rim => {
                crate::rim_arc::revolved_arc_fillet(
                    model,
                    solid,
                    edge,
                    c,
                    *radius,
                    None,
                    Some((index, &met.mates)),
                    Some(&aside),
                    tol,
                )?;
            }
            _ => {
                crate::marched::marched_wedge(
                    model,
                    solid,
                    edge,
                    *radius,
                    Some((index, &met.mates)),
                    &aside,
                    tol,
                )?;
            }
        }
        bands.extend(aside.into_inner().into_iter().map(|w| (index, w)));
    }
    // Each pipe's neighbours end on its sections, and meet it there.
    for (link, pipe, ends) in pipes {
        for end in &ends {
            let reach = once[end.neighbour].1 * 3.0;
            let mut trimmed = false;
            for (index, band) in &mut bands {
                if *index == end.neighbour && band.additive {
                    band.solid = crate::support::trimmed_wedge(
                        model,
                        &band.solid,
                        (end.centre, -end.into),
                        reach,
                        tol,
                    )?;
                    trimmed = true;
                }
            }
            if !trimmed || !met.flush_at([link, end.neighbour], (end.centre, end.into), end.foot) {
                return Ok(None);
            }
        }
        bands.push((link, pipe));
    }
    let wedges = super::marched_junctions::joined(model, solid, once, &met, bands, tol)?;
    let Some(built) = super::marched_seats::applied(model, solid, once, &wedges, tol)? else {
        return Ok(None);
    };
    super::guards::solid_fits(model, &built.shape, tol)?;
    Ok(Some(built))
}

/// The exact planar section swept between its end planes. An overhanging
/// face clips the tool before application, so no material behind it is cut.
fn straight(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    (index, mates): (usize, &[Mate]),
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    let mut seat = crate::support::planar_seat(model, solid, edge, tol)?;
    let mut trims = Vec::<(Point, Vector)>::new();
    for start in [true, false] {
        let at = if start { seat.start } else { seat.end };
        if !seat.convex || Mate::meets_at(mates, index, at, tol) {
            continue;
        }
        let outward = seat.along * if start { -1.0 } else { 1.0 };
        let Some((on, n)) = crate::support::overhang_at(
            model,
            solid,
            edge,
            [&seat.faces[0], &seat.faces[1]],
            at,
            outward,
            tol,
        )?
        else {
            continue;
        };
        let section = PlanarSection::new(&seat, &RoundSpec::Circular(radius), tol)?;
        let s = if start { 0.0 } else { section.length() };
        let reach = [at, section.contact(0, s), section.contact(1, s)]
            .iter()
            .map(|p| (on - *p).dot(n) / outward.dot(n))
            .fold(0.0, f64::max)
            + radius * 0.25;
        if start {
            seat.start -= seat.along * reach;
        } else {
            seat.end += seat.along * reach;
        }
        trims.push((at - n * (at - on).dot(n), n));
    }
    let section = PlanarSection::new(&seat, &RoundSpec::Circular(radius), tol)?;
    let length = section.length();
    let mut tool = super::planar_tools::sweep(model, edge, &section, tol)?.solid;
    for plane in trims {
        tool = crate::support::trimmed_wedge(model, &tool, plane, length + radius * 8.0, tol)?;
    }
    let checked = ogeom_algo::check(model, &tool, tol)?;
    if !checked.is_valid() {
        ogeom_bail!(
            Construction,
            "a contour tool fails its check: {:?}",
            checked.problems
        );
    }
    Ok(Wedge {
        solid: tool,
        additive: !seat.convex,
        edges: vec![edge.clone()],
        melt_refusal: None,
        credit: Vec::new(),
    })
}

/// A pipe: its link by its index in the request, its tool and its ends.
type Piped = (usize, Wedge, [super::pipe_corner::End; 2]);

/// The pipes of `once` on `solid`: each link whose ball changes walls off
/// both its vertices, with its tool and its two ends. `None` where an
/// off-vertex crossing is not one of a pipe's two ends, or a pipe will not
/// build.
fn pipes(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Piped>>> {
    let faces = super::switch_over::Faces::of(model, solid)?;
    let off: Vec<super::switch_over::Crossing> =
        super::switch_over::crossings(model, solid, &faces, once, tol)?
            .into_iter()
            .filter(|c| c.by > super::switch_over::off_by(tol))
            .collect();
    let mut built = Vec::new();
    for link in 0..once.len() {
        // The crossings at the link's ends that stand over its neighbours.
        let here: Vec<&super::switch_over::Crossing> = off
            .iter()
            .filter(|c| c.pair.contains(&link) && c.pair[c.over] != link)
            .collect();
        let pair: [&super::switch_over::Crossing; 2] = match here.as_slice() {
            [] => continue,
            [a, b] => [a, b],
            _ => return Ok(None),
        };
        let Some((wedge, ends)) = super::pipe_corner::pipe(model, &faces, once, link, pair, tol)?
        else {
            return Ok(None);
        };
        built.push((link, wedge, ends));
    }
    let covered = |c: &super::switch_over::Crossing| {
        built
            .iter()
            .any(|(link, ..)| c.pair.contains(link) && c.pair[c.over] != *link)
    };
    if !off.iter().all(covered) {
        return Ok(None);
    }
    Ok(Some(built))
}
