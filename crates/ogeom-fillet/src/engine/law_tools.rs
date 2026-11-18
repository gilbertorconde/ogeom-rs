//! Radius laws laid over an edge walk, with one tool per joined strip.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

use super::fallback::{Sense, Tool};
use crate::law::{Profile, RadiusLaw, Ruler};
use crate::round::section::{PlanarSection, RoundSpec};

/// Lay the law on the chain's arc-length walk, build sections on the
/// original seats, join their shared end sections, and apply their tools.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut scratch = model.clone();
    match build(&mut scratch, solid, chain, law, tol) {
        Ok(Some(built)) => {
            *model = scratch;
            Ok(Some(built))
        }
        other => other,
    }
}

fn build(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let closed = crate::fillet::closed_chain(model, solid, chain)?;
    let profile = Profile::new(law, closed, tol)?;
    if profile.level(tol) {
        let asked: Vec<_> = chain
            .iter()
            .map(|e| (e.clone(), profile.largest()))
            .collect();
        if let Some(built) = crate::fillet::rounds_by_engine(model, solid, &asked, tol)? {
            return Ok(Some(built));
        }
    }
    let span = if chain.len() > 1 && !closed {
        crate::fillet::one_curve_span(model, solid, chain, tol)?
    } else {
        None
    };
    if chain.len() == 1 || closed || span.is_some() {
        if crate::support::planar_seat(model, solid, &chain[0], tol).is_ok() {
            let reversed = chain[0].orientation() == Orientation::Reversed;
            let law = profile.scaled_by(1.0, reversed);
            return super::planar_tools::linear_round(model, solid, &chain[0], &law, tol).map(Some);
        }
        if closed && chain.len() > 1 {
            if !crate::fillet::law_loop_fits(model, solid, chain, &profile, tol)? {
                for edge in chain {
                    crate::fillet::law_ball_on_faces(model, solid, edge, &profile, None, tol)?;
                }
            }
        } else {
            crate::fillet::law_ball_on_faces(
                model,
                solid,
                &chain[0],
                &profile,
                (chain.len() == 1).then_some(chain[0].orientation() == Orientation::Reversed),
                tol,
            )?;
        }
        let tool = marched(model, solid, &chain[0], &profile, span, None, tol)?;
        let Some(mut built) =
            super::planar_tools::applied(model, solid, &chain[..1], &[tool], tol)?
        else {
            return Ok(None);
        };
        for edge in &chain[1..] {
            built.history.delete(edge);
        }
        return Ok(Some(built));
    }
    let mates = crate::fillet::chain_mates(model, chain, tol)?;
    let runs: Vec<_> = chain
        .iter()
        .zip(&mates)
        .map(|(e, m)| {
            let [a, b] = m.ends.map(|(p, _)| p);
            if e.orientation() == Orientation::Reversed {
                [b, a]
            } else {
                [a, b]
            }
        })
        .collect();
    let near = |p: Point, q: Point| p.distance(q) <= tol.confusion() * 1e3;
    let mut walk = vec![(0, true)];
    let mut used = vec![false; chain.len()];
    used[0] = true;
    let mut head = runs[0][1];
    while let Some(k) =
        (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, head)))
    {
        let forward = near(runs[k][0], head);
        used[k] = true;
        head = runs[k][usize::from(forward)];
        walk.push((k, forward));
    }
    let closed = near(head, runs[0][0]);
    if !closed {
        let mut tail = runs[0][0];
        while let Some(k) =
            (0..chain.len()).find(|&k| !used[k] && runs[k].iter().any(|p| near(*p, tail)))
        {
            let forward = near(runs[k][1], tail);
            used[k] = true;
            tail = runs[k][usize::from(!forward)];
            walk.insert(0, (k, forward));
        }
    }
    if used.iter().any(|u| !u) {
        return Ok(None);
    }
    let profile = Profile::new(law, closed, tol)?;
    let lengths: Vec<f64> = chain
        .iter()
        .map(|e| {
            let (curve, range) = crate::support::edge_curve(model, e, tol)?;
            Ok(Ruler::new(curve, range, false, tol)?.length())
        })
        .collect::<OgeomResult<_>>()?;
    let total: f64 = lengths.iter().sum();
    let mut laws = vec![profile.clone(); chain.len()];
    let mut before = 0.0;
    for &(k, forward) in &walk {
        let (a, b) = (before / total, (before + lengths[k]) / total);
        laws[k] = if forward {
            profile.window(a, b)
        } else {
            profile.window(b, a)
        };
        before += lengths[k];
    }
    let tangent = crate::engine::junction::tangent_junctions(chain, &mates, tol);
    let planar: Vec<_> = chain
        .iter()
        .map(|e| crate::support::planar_seat(model, solid, e, tol).is_ok())
        .collect();
    for i in 0..chain.len() {
        for j in i + 1..chain.len() {
            if planar[i] == planar[j] {
                continue;
            }
            for (p, _) in mates[i].ends {
                if mates[j].ends.iter().any(|(q, _)| near(p, *q))
                    && !tangent
                        .iter()
                        .any(|t| t.pair.contains(&i) && t.pair.contains(&j) && near(t.at, p))
                {
                    ogeom_core::ogeom_bail!(
                        Construction,
                        "a variable radius over a chain turning where a marched band meets a straight one between planes is still owed: under a law the marched band is fitted, and its end meets the straight band in what the boolean does not melt; see docs/CAPABILITIES.md, fillet.variable"
                    );
                }
            }
        }
    }
    crate::fillet::law_chain_fits(model, solid, chain, &laws.iter().collect::<Vec<_>>(), tol)?;
    if planar.iter().all(|p| !p) && tangent.is_empty() {
        return resolved(model, solid, chain, &laws, &mates, tol).map(Some);
    }
    let mut tools = Vec::with_capacity(chain.len());
    for (k, edge) in chain.iter().enumerate() {
        if let Ok(mut seat) = crate::support::planar_seat(model, solid, edge, tol) {
            let mut law = laws[k].scaled_by(1.0, edge.orientation() == Orientation::Reversed);
            let tangent = crate::engine::junction::tangent_junctions(chain, &mates, tol);
            let sharp = |p: Point| {
                (0..chain.len()).any(|j| j != k && mates[j].ends.iter().any(|(q, _)| near(p, *q)))
                    && !tangent.iter().any(|t| t.pair.contains(&k) && near(t.at, p))
            };
            let reach = 4.0 * law.largest();
            let before = if sharp(seat.start) { reach } else { 0.0 };
            let after = if sharp(seat.end) { reach } else { 0.0 };
            let length = seat.start.distance(seat.end);
            seat.start -= seat.along * before;
            seat.end += seat.along * after;
            law = law.stretched(before / length, after / length);
            let unit = PlanarSection::new(&seat, &RoundSpec::Circular(1.0), tol)?;
            let distance = law.scaled_by(unit.distances(0.0)[0], false);
            let section = PlanarSection::new(
                &seat,
                &RoundSpec::Varying {
                    laws: [distance.clone(), distance],
                    profile: crate::RoundProfile::Elliptical,
                },
                tol,
            )?;
            tools.push(super::planar_tools::sweep(model, edge, &section, tol)?);
        } else {
            tools.push(marched(
                model,
                solid,
                edge,
                &laws[k],
                None,
                Some((k, &mates)),
                tol,
            )?);
        }
    }
    // Join only shared tangent sections. Sharp meetings remain separate
    // tools so their overlapping cuts take the union of the removals.
    let mut groups: Vec<usize> = (0..chain.len()).collect();
    let gap = crate::marched::band_fit_target(tol) * 2.0;
    let tangent = crate::engine::junction::tangent_junctions(chain, &mates, tol);
    let mut sharp = Vec::new();
    for i in 0..chain.len() {
        for j in i + 1..chain.len() {
            for (p, _) in mates[i].ends {
                if mates[j].ends.iter().any(|(q, _)| near(p, *q))
                    && !tangent
                        .iter()
                        .any(|t| t.pair.contains(&i) && t.pair.contains(&j) && near(t.at, p))
                {
                    sharp.push([i, j]);
                }
            }
        }
    }
    for joint in tangent {
        let [i, j] = joint.pair.map(|k| groups[k]);
        if sharp.iter().any(|[a, b]| {
            (groups[*a] == i && groups[*b] == j) || (groups[*a] == j && groups[*b] == i)
        }) {
            continue;
        }
        if i == j {
            let joined = crate::engine::junction::joined(
                model,
                [&tools[i].solid, &tools[i].solid],
                joint.flush,
                joint.at,
                gap,
                tol,
            )?;
            tools[i].solid = joined.0;
        } else {
            let joined = crate::engine::junction::joined(
                model,
                [&tools[i].solid, &tools[j].solid],
                joint.flush,
                joint.at,
                gap,
                tol,
            )?;
            let source: Vec<_> = tools[i]
                .faces
                .iter()
                .chain(&tools[j].faces)
                .cloned()
                .collect();
            tools[i].faces = joined
                .1
                .iter()
                .filter_map(|(from, to)| {
                    source
                        .iter()
                        .find(|(f, _)| f.is_same(from))
                        .map(|(_, owner)| (to.clone(), owner.clone()))
                })
                .collect();
            tools[i].solid = joined.0;
            for group in &mut groups {
                if *group == j {
                    *group = i;
                }
            }
        }
    }
    let mut roots = groups;
    roots.sort_unstable();
    roots.dedup();
    let tools: Vec<_> = tools
        .into_iter()
        .enumerate()
        .filter_map(|(i, tool)| roots.contains(&i).then_some(tool))
        .collect();
    super::planar_tools::applied(model, solid, chain, &tools, tol)
}

/// Resolve intersecting curved strips in a fixed guide order. A shortened
/// crease still has its original support surfaces; only its interval is
/// trimmed by the preceding strip, and its law is restricted by arc length.
fn resolved(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    laws: &[Profile],
    mates: &[crate::fillet::Mate],
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::{Curve, Curve3d as _};
    let mut order = Vec::with_capacity(chain.len());
    for (i, edge) in chain.iter().enumerate() {
        let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
        order.push((!matches!(curve, Curve::Line(_)), super::id(edge), i));
    }
    order.sort_unstable();
    let mut built = Built::from_nothing(solid.clone());
    for (_, _, i) in order {
        let whole = &chain[i];
        let hosts = ogeom_topo::ancestors_of(model, solid, whole, ShapeType::Face)?;
        let pieces = crate::fillet::refind_edges(model, &built.shape, whole, &hosts, tol)?;
        let (curve, range) = crate::support::edge_curve(model, whole, tol)?;
        let ruler = Ruler::new(curve.clone(), range, false, tol)?;
        let reversed = whole.orientation() == Orientation::Reversed;
        for piece in pieces {
            let (own, range) = crate::support::edge_curve(model, &piece, tol)?;
            let (a, b) = if piece.orientation() == Orientation::Reversed {
                (range.1, range.0)
            } else {
                range
            };
            let fraction = |p| -> OgeomResult<f64> {
                let t = ogeom_algo::project_on_curve(&curve, p, 64, tol)?.parameter;
                let s = ruler.fraction(t, tol)?.0;
                Ok(if reversed { 1.0 - s } else { s })
            };
            let law = laws[i].window(
                fraction(own.point_at(a, tol)?)?,
                fraction(own.point_at(b, tol)?)?,
            );
            let mut tool = marched(
                model,
                &built.shape,
                &piece,
                &law,
                None,
                Some((i, mates)),
                tol,
            )?;
            for (_, owner) in &mut tool.faces {
                *owner = whole.clone();
            }
            let next = super::planar_tools::applied(
                model,
                &built.shape,
                std::slice::from_ref(whole),
                &[tool],
                tol,
            )?
            .ok_or_else(|| {
                ogeom_core::ogeom_err!(Construction, "a trimmed law strip leaves no blend")
            })?;
            built = Built::new(next.shape, built.history.then(&next.history));
        }
        built.history.delete(whole);
    }
    Ok(built)
}

fn marched(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    law: &Profile,
    span: Option<(f64, f64)>,
    mates: Option<(usize, &[crate::fillet::Mate])>,
    tol: Tolerances,
) -> OgeomResult<Tool> {
    let (built, wedges) = crate::support::collecting_wedges(|| {
        crate::marched::marched_fillet_law(model, solid, edge, law, span, mates, tol)
    });
    built?;
    let [wedge] = <[crate::support::Wedge; 1]>::try_from(wedges).map_err(|_| {
        ogeom_core::ogeom_err!(Construction, "a law's marched strip does not make one tool")
    })?;
    let faces = explore_unique(model, &wedge.solid, ShapeType::Face)?
        .into_iter()
        .map(|f| (f, edge.clone()))
        .collect();
    Ok(Tool {
        solid: wedge.solid,
        sense: if wedge.additive {
            Sense::Fuse
        } else {
            Sense::Cut
        },
        faces,
    })
}
