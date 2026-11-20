//! A tangent chain of asymmetric rounds whose distances follow laws laid
//! by arc length along the whole chain.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::Point;
use ogeom_topo::{Model, Orientation, Shape, ShapeType};

use crate::asym::{Asymmetric, RoundProfile};
use crate::law::{Profile, RadiusLaw, Ruler};
use crate::round::section::{PlanarSection, RoundSpec};

/// The chain `chain` of `solid` rounded with a band meeting `face` at the
/// distance `laws[0]` gives and the chain's other faces at `laws[1]`'s,
/// each law laid by arc length from the start of the walk along the chain
/// (the first edge's start, the way it runs, carried on back past it where
/// the chain is asked from its middle), closed round a closed chain. Each
/// straight edge between planes takes its exact sweep under its stretch of
/// the laws, any other edge the marched band, each stopped in its end
/// sections at a junction, and the bands are joined across every junction
/// into one tool per run.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where `face` is not beside every edge, the edges do not run end to end,
/// two meet other than tangent, a law states no usable distance or an
/// open law is asked round a closed chain, a contact sets back past its
/// face's side, or the bands do not join; and as the bands themselves.
#[allow(clippy::too_many_lines, reason = "one chain, walked end to end")]
pub(crate) fn varying(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    face: &Shape,
    laws: [&RadiusLaw; 2],
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for edge in chain {
        let beside = ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?;
        if !beside.iter().any(|f| f.is_same(face)) {
            ogeom_bail!(
                Construction,
                "a chain of asymmetric rounds names one face beside every edge, and an \
                 edge of this chain does not lie on the face named"
            );
        }
    }
    let mates = crate::fillet::chain_mates(model, chain, tol)?;
    // Each edge's ends the way it runs, and the walk from the first edge.
    let runs: Vec<[Point; 2]> = chain
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
        ogeom_bail!(
            Construction,
            "the edges of a chain of asymmetric rounds do not run end to end"
        );
    }
    // Every end two edges share is a tangent junction.
    let tangent = super::junction::tangent_junctions(chain, &mates, tol);
    for i in 0..chain.len() {
        for j in i + 1..chain.len() {
            for (p, _) in mates[i].ends {
                if mates[j].ends.iter().any(|(q, _)| near(p, *q))
                    && !tangent
                        .iter()
                        .any(|t| t.pair.contains(&i) && t.pair.contains(&j) && near(t.at, p))
                {
                    ogeom_bail!(
                        Construction,
                        "two edges of a chain of asymmetric rounds whose distances vary meet \
                         other than tangent, and such a chain runs only where each edge \
                         continues the last; see docs/CAPABILITIES.md, \
                         fillet.asymmetric-rounds"
                    );
                }
            }
        }
    }
    let profiles = [
        Profile::new(laws[0], closed, tol)?,
        Profile::new(laws[1], closed, tol)?,
    ];
    if profiles.iter().all(|p| p.level(tol)) {
        let request: Vec<(Shape, crate::Blend)> = chain
            .iter()
            .map(|e| {
                (
                    e.clone(),
                    crate::Blend::RoundDistances {
                        face: face.clone(),
                        on_face: profiles[0].largest(),
                        on_other: profiles[1].largest(),
                        profile,
                    },
                )
            })
            .collect();
        return crate::blend_edges(model, solid, &request, crate::ChamferCorner::Point, tol);
    }
    let weight = match profile {
        RoundProfile::Elliptical => crate::marched::ConicWeight::Elliptical,
        RoundProfile::Conic(_) => crate::marched::ConicWeight::Fixed(profile.weight(0.0)?),
    };
    let lengths: Vec<f64> = chain
        .iter()
        .map(|e| {
            let (curve, range) = crate::support::edge_curve(model, e, tol)?;
            Ok(Ruler::new(curve, range, false, tol)?.length())
        })
        .collect::<OgeomResult<_>>()?;
    let total: f64 = lengths.iter().sum();
    // Each edge's stretch of the two laws, the way the edge runs.
    let mut stretches: Vec<[Profile; 2]> = vec![profiles.clone(); chain.len()];
    let mut before = 0.0;
    for &(k, forward) in &walk {
        let (a, b) = (before / total, (before + lengths[k]) / total);
        stretches[k] = profiles.clone().map(|p| {
            if forward {
                p.window(a, b)
            } else {
                p.window(b, a)
            }
        });
        before += lengths[k];
    }

    let mut scratch = model.clone();
    let mut tools = Vec::with_capacity(chain.len());
    for (k, edge) in chain.iter().enumerate() {
        let at_junction = |p: Point| tangent.iter().any(|t| t.pair.contains(&k) && near(t.at, p));
        if let Ok((seat, named)) = crate::asym::straight_seat(&scratch, solid, edge, face, tol) {
            // The seat reads its laws from its stored start.
            let reversed = edge.orientation() == Orientation::Reversed;
            let [on_face, on_other] = stretches[k].clone().map(|p| p.scaled_by(1.0, reversed));
            let spec = RoundSpec::Varying {
                laws: if named == 0 {
                    [on_face, on_other]
                } else {
                    [on_other, on_face]
                },
                profile,
            };
            let section = PlanarSection::new(&seat, &spec, tol)?;
            for j in 1..16 {
                let s = section.length() * f64::from(j) / 16.0;
                for (i, f) in seat.faces.iter().enumerate() {
                    if !crate::support::on_face_or_boundary(
                        &scratch,
                        solid,
                        f,
                        section.contact(i, s),
                        tol,
                    )? {
                        ogeom_bail!(
                            Construction,
                            "an asymmetric round sets back {} along a face of the chain's edge, \
                             at {:.3} of it, past that face's side; see docs/CAPABILITIES.md, \
                             fillet.asymmetric-rounds",
                            section.distances(s)[i],
                            s / section.length()
                        );
                    }
                }
            }
            let extend = [!at_junction(seat.start), !at_junction(seat.end)];
            tools.push(super::conic_tools::extended_sweep(
                &mut scratch,
                solid,
                edge,
                &seat,
                &spec,
                extend,
                tol,
            )?);
        } else {
            let asked = Asymmetric {
                edge: edge.clone(),
                face: face.clone(),
                on_face: stretches[k][0].largest(),
                on_other: stretches[k][1].largest(),
                profile,
            };
            tools.push(super::conic_tools::marched(
                &mut scratch,
                solid,
                &asked,
                Some(&stretches[k]),
                Some((k, &mates)),
                weight,
                tol,
            )?);
        }
    }
    // One tool per run: the bands joined across each junction, where their
    // end sections stand in one plane square to both creases.
    let gap = crate::marched::band_fit_target(tol) * 2.0;
    let mut groups: Vec<usize> = (0..chain.len()).collect();
    for joint in &tangent {
        let [i, j] = joint.pair.map(|k| groups[k]);
        let wedges = [&tools[i].solid, &tools[j].solid];
        let (solid_joined, traced) =
            super::junction::joined(&mut scratch, wedges, joint.flush, joint.at, gap, tol)?;
        if i != j {
            let source: Vec<_> = tools[i]
                .faces
                .iter()
                .chain(&tools[j].faces)
                .cloned()
                .collect();
            tools[i].faces = traced
                .iter()
                .filter_map(|(from, to)| {
                    source
                        .iter()
                        .find(|(f, _)| f.is_same(from))
                        .map(|(_, owner)| (to.clone(), owner.clone()))
                })
                .collect();
            for group in &mut groups {
                if *group == j {
                    *group = i;
                }
            }
        } else {
            let source = tools[i].faces.clone();
            tools[i].faces = traced
                .iter()
                .filter_map(|(from, to)| {
                    source
                        .iter()
                        .find(|(f, _)| f.is_same(from))
                        .map(|(_, owner)| (to.clone(), owner.clone()))
                })
                .collect();
        }
        tools[i].solid = solid_joined;
    }
    let mut roots = groups;
    roots.sort_unstable();
    roots.dedup();
    let tools: Vec<_> = tools
        .into_iter()
        .enumerate()
        .filter_map(|(i, tool)| roots.contains(&i).then_some(tool))
        .collect();
    let built = super::planar_tools::applied(&mut scratch, solid, chain, &tools, tol)?.ok_or_else(
        || {
            ogeom_core::ogeom_err!(
                Construction,
                "a chain of asymmetric rounds whose distances vary does not close on its solid"
            )
        },
    )?;
    *model = scratch;
    Ok(built)
}
