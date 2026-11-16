//! Where marched bands meet at a junction (ADR 0001, amendment "marched
//! seats", junctions): two of a request's creases turning at a shared
//! vertex (a rib's spline crease into a boss's rim) meet as the bands of
//! one round do. Each band is fitted to the plane the two meet in (carried
//! on past its cap across a re-entrant turn, trimmed back across a convex
//! one), a flush pair joined across it into one wedge, a crossed re-entrant
//! pair's corner built from the two. The wedges stay the engine's tools;
//! the one-edge entries' bands meet through the same steps.

use core::cell::RefCell;

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fillet::Mate;
use crate::junction::Junction;
use crate::support::{Wedge, edge_curve};

/// A request's edges as their bands meet: each edge's mates by its ends,
/// the junctions the bands are joined at, and each edge's faces.
pub(crate) struct Met {
    pub(crate) mates: Vec<Mate>,
    junctions: Vec<Junction>,
    hosts: Vec<Vec<Shape>>,
}

/// How `edges` (rounded at `radii`, with `ends` their vertices) meet on
/// `solid`. `None` where two of them meet at a vertex unjoined and not
/// both are marched (a marched band meeting one of another kind there is
/// built after it, one at a time), or the junctions are not read.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled, and the reading
/// of an edge's curve.
pub(crate) fn met(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    ends: &[Vec<Shape>],
    all_marched: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Met>> {
    let mut mates = request_mates(model, edges, &[], tol)?;
    let mut host_faces: ogeom_core::FastMap<ogeom_topo::TShapeId, Vec<Shape>> =
        ogeom_core::FastMap::default();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        for e in explore_unique(model, &face, ShapeType::Edge)? {
            let held = host_faces.entry(e.node()).or_default();
            if !held.iter().any(|f| f.is_same(&face)) {
                held.push(face.clone());
            }
        }
    }
    let hosts: Vec<Vec<Shape>> = edges
        .iter()
        .map(|e| host_faces.get(&e.node()).cloned().unwrap_or_default())
        .collect();
    let found = match crate::junction::junctions(model, solid, edges, &hosts, &mates, radii, tol) {
        Ok(found) => found,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Err(_) => Vec::new(),
    };
    let junctions = joinable(model, edges, found)?;
    for m in 0..edges.len() {
        for n in m + 1..edges.len() {
            let meet = ends[m].iter().any(|v| ends[n].iter().any(|w| w.is_same(v)));
            let joined = junctions
                .iter()
                .any(|j| j.slot(m).is_some() && j.slot(n).is_some());
            // Marched bands meeting unjoined run out through each other,
            // each as it is alone.
            if meet && !joined && !all_marched {
                return Ok(None);
            }
        }
    }
    for junction in &junctions {
        for &index in &junction.pair {
            let mate = &mut mates[index];
            let end = usize::from(
                mate.ends[1].0.distance(junction.at) < mate.ends[0].0.distance(junction.at),
            );
            if junction.reentrant || !junction.crossed {
                mate.meets[end] = true;
            } else {
                mate.alone[end] = true;
            }
        }
    }
    Ok(Some(Met {
        mates,
        junctions,
        hosts,
    }))
}

/// The band of an edge of another kind meeting a marched one at a junction,
/// set `aside` on `solid` as given: a rim's arc its exact torus band (or
/// the march's where a cone's ball would reach its axis).
///
/// # Errors
///
/// As the band's builder; [`OgeomError::NotDone`] for a straight crease
/// or a whole rim, whose bands do not meet a marched one here.
pub(crate) fn other_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    mates: (usize, &[Mate]),
    aside: &RefCell<Vec<Wedge>>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (curve, crange) = crate::support::edge_curve(model, edge, tol)?;
    let Curve::Circle(c) = curve else {
        ogeom_bail!(
            NotDone,
            "a straight crease's band does not meet a marched one at a junction"
        );
    };
    let seat = crate::support::revolved_seat(model, solid, edge, &c, tol)?;
    if seat.slope != 0.0 && seat.ball(radius).rho - radius <= tol.confusion() {
        return crate::marched::marched_wedge(model, solid, edge, radius, Some(mates), aside, tol);
    }
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    if closed || crate::support::runs_whole_turn(model, &seat, &c, crange, tol)? {
        ogeom_bail!(
            NotDone,
            "a whole rim's band does not meet a marched one at a junction"
        );
    }
    crate::rim_arc::revolved_arc_fillet(
        model,
        solid,
        edge,
        &c,
        radius,
        None,
        Some(mates),
        Some(aside),
        tol,
    )
}

/// The `set_aside` wedges (each with the index of the edge of `once` it
/// rounds) fitted to `met`'s junctions and joined across them, built on
/// `solid`.
///
/// # Errors
///
/// As fitting and joining the wedges.
pub(crate) fn joined(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    met: &Met,
    set_aside: Vec<(usize, Wedge)>,
    tol: Tolerances,
) -> OgeomResult<Vec<Wedge>> {
    let mut fitted: Vec<(Vec<usize>, Shape, Wedge)> = Vec::with_capacity(set_aside.len());
    for (index, mut wedge) in set_aside {
        let (edge, radius) = &once[index];
        fit_to_junctions(
            model,
            &mut wedge,
            (index, edge, &met.hosts[index]),
            &met.junctions,
            *radius,
            tol,
        )?;
        fitted.push((vec![index], edge.clone(), wedge));
    }
    let halves = with_corners(model, &mut fitted, &met.junctions, solid, false, tol)?;
    let joined = join_flush(model, fitted, &met.junctions, &halves, tol)?;
    Ok(joined.into_iter().map(|(_, wedge)| wedge).collect())
}

/// The mates of a request's `edges` by their ends, each end settled where
/// it stands at one of the `rounded` corners.
pub(crate) fn request_mates(
    model: &Model,
    edges: &[Shape],
    rounded: &[Point],
    tol: Tolerances,
) -> OgeomResult<Vec<Mate>> {
    use ogeom_geom::Curve3d as _;
    edges
        .iter()
        .map(|e| -> OgeomResult<Mate> {
            let (curve, range) = edge_curve(model, e, tol)?;
            let end_of = |t: f64, leaving: f64| -> OgeomResult<(Point, Vector)> {
                let d = curve.d1_at(t, tol)?;
                let m = d.magnitude();
                let unit = if m > tol.confusion() {
                    d / m
                } else {
                    Vector::ZERO
                };
                Ok((curve.point_at(t, tol)?, unit * leaving))
            };
            let ends = [end_of(range.0, 1.0)?, end_of(range.1, -1.0)?];
            let settled = ends.map(|(p, _)| {
                rounded
                    .iter()
                    .any(|q| q.distance(p) <= tol.confusion() * 1e3)
            });
            Ok(Mate {
                ends,
                settled,
                meets: [false, false],
                alone: [false, false],
            })
        })
        .collect()
}

/// A convex band's wedge fitted to the junctions at its edge's ends (the
/// request's `index`-th edge `edge`, between `hosts` on the solid as
/// given): across a re-entrant flush turn carried on past its cap to the
/// plane it meets the other band in, across a convex one trimmed back to
/// it. A crossed junction's bands stand as built: a re-entrant one's corner
/// is built from the two once both are set aside.
pub(crate) fn fit_to_junctions(
    model: &mut Model,
    wedge: &mut Wedge,
    (index, edge, hosts): (usize, &Shape, &[Shape]),
    junctions: &[Junction],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<()> {
    if wedge.additive {
        return Ok(());
    }
    for junction in junctions {
        let Some(slot) = junction.slot(index) else {
            continue;
        };
        if junction.crossed {
            continue;
        }
        let plane = junction.keeps(slot);
        // Nothing of the band but its end stands within two radii of the
        // vertex across the plane: a band of its own edge coming back that
        // close would leave no room for its ball.
        if junction.reentrant {
            let Some(cap) = crate::junction::cap_at(
                model,
                &wedge.solid,
                junction.at,
                junction.outward[slot],
                tol,
            )?
            else {
                ogeom_bail!(
                    Construction,
                    "a band meeting another across a re-entrant turn has \
                     no cap at the vertex to carry on from"
                );
            };
            wedge.solid = crate::junction::carried_into(
                model,
                &wedge.solid,
                &cap,
                (edge, hosts),
                junction.outward[slot],
                plane,
                tol,
            )?;
        } else {
            wedge.solid =
                crate::support::trimmed_wedge(model, &wedge.solid, plane, radius * 2.0, tol)?;
        }
    }
    Ok(())
}

/// The junctions of `found` whose bands can be joined: taken in turn, each
/// joins the groups of bands already joined at its two edges, and is
/// dropped where two edges of those groups also meet at another vertex,
/// which joined into one round would lose the later band's running on
/// through the earlier there.
pub(crate) fn joinable(
    model: &Model,
    edges: &[Shape],
    found: Vec<Junction>,
) -> OgeomResult<Vec<Junction>> {
    let mut ends: Vec<Vec<Shape>> = Vec::with_capacity(edges.len());
    for edge in edges {
        ends.push(
            ogeom_algo::edge_vertices(model, edge)?
                .map(|(a, b)| vec![a, b])
                .unwrap_or_default(),
        );
    }
    // Flush junctions are taken first: a crossed junction joins its two
    // bands only where that leaves every flush one standing.
    let (flush, crossed): (Vec<_>, Vec<_>) = found.into_iter().partition(|j| !j.crossed);
    let mut found = flush.clone();
    found.extend(crossed);
    let pending: Vec<[usize; 2]> = found.iter().map(|j| j.pair).collect();
    let mut group: Vec<usize> = (0..edges.len()).collect();
    let mut kept: Vec<Junction> = Vec::new();
    for junction in found {
        let [a, b] = junction.pair;
        if group[a] == group[b] {
            // A junction closing a loop of joined bands joins no groups,
            // and is taken for its corner. Each junction stands at a vertex
            // of its own, so two bands meeting at both their ends (a pocket
            // of a chord and its arc) close a loop of two there.
            kept.push(junction);
            continue;
        }
        // Two edges of the groups this joins that meet elsewhere must meet
        // there in a junction of their own: one taken, this one, or a
        // crossed one still to come.
        let joined = |m: usize, n: usize| {
            kept.iter()
                .chain(std::iter::once(&junction))
                .any(|j| j.slot(m).is_some() && j.slot(n).is_some())
                || pending
                    .iter()
                    .any(|pair| pair.contains(&m) && pair.contains(&n))
        };
        let of_a: Vec<usize> = (0..edges.len()).filter(|&m| group[m] == group[a]).collect();
        let of_b: Vec<usize> = (0..edges.len()).filter(|&m| group[m] == group[b]).collect();
        let clash = of_a.iter().any(|&m| {
            of_b.iter().any(|&n| {
                !joined(m, n) && ends[m].iter().any(|v| ends[n].iter().any(|w| w.is_same(v)))
            })
        });
        if clash {
            continue;
        }
        let (from, to) = (group[a], group[b]);
        for g in &mut group {
            if *g == from {
                *g = to;
            }
        }
        kept.push(junction);
    }
    // Every two edges of one joined group that share a vertex meet there
    // in a kept junction; where one was dropped, the request is taken with
    // its flush junctions alone.
    let meets_unjoined = (0..edges.len()).any(|m| {
        (m + 1..edges.len()).any(|n| {
            group[m] == group[n]
                && ends[m].iter().any(|v| ends[n].iter().any(|w| w.is_same(v)))
                && !kept
                    .iter()
                    .any(|j| j.slot(m).is_some() && j.slot(n).is_some())
        })
    });
    if meets_unjoined && kept.iter().any(|j| j.crossed) {
        return joinable(model, edges, flush);
    }
    Ok(kept)
}

/// Join, among the wedges a round set aside, the two bands of each crossed
/// re-entrant junction the round built into the one wedge that also takes
/// their corner, in the material of `solid`, the solid the round's wedges
/// apply to. The joined wedge is set aside with the first band's edge, its
/// faces from the second band credited to the second's.
///
/// A band with a crossed junction at each end is first cut in two in the
/// plane square to its crease halfway along: each half is the straight
/// band of its own end section, joined into its own corner, and the two
/// corners' wedges are joined again across that plane, where both hold the
/// band's one section. Returns those planes, each with the band's index.
/// `law` says the bands are under a radius law, their sections changing up
/// to the vertex, as [`crate::junction::corner`] takes them.
pub(crate) fn with_corners(
    model: &mut Model,
    set_aside: &mut Vec<(Vec<usize>, Shape, Wedge)>,
    junctions: &[Junction],
    solid: &Shape,
    law: bool,
    tol: Tolerances,
) -> OgeomResult<Vec<(usize, (Point, Vector))>> {
    let crossed: Vec<&Junction> = junctions
        .iter()
        .filter(|j| j.crossed && j.reentrant)
        .collect();
    let holding = |set_aside: &[(Vec<usize>, Shape, Wedge)], index: usize| {
        (0..set_aside.len())
            .filter(|&k| set_aside[k].0.contains(&index))
            .collect::<Vec<usize>>()
    };
    let mut halves = Vec::new();
    let indices: Vec<usize> = crossed.iter().flat_map(|j| j.pair).collect();
    for &index in &indices {
        let ends: Vec<Point> = crossed
            .iter()
            .filter(|j| j.slot(index).is_some())
            .map(|j| j.at)
            .collect();
        if ends.len() != 2 || halves.iter().any(|(i, _)| *i == index) {
            continue;
        }
        let held = holding(set_aside, index);
        let [k] = held.as_slice() else {
            continue;
        };
        let k = *k;
        let along = ends[1] - ends[0];
        let length = along.magnitude();
        if length <= tol.confusion() {
            ogeom_bail!(Construction, "a band crossed at both ends has no length");
        }
        let toward = along / length;
        let middle = ends[0].midpoint(ends[1]);
        let span = length + crate::marched::band_fit_target(tol).max(length);
        let (indices_k, edge, wedge) = set_aside.remove(k);
        for keep in [-toward, toward] {
            let half =
                crate::support::trimmed_wedge(model, &wedge.solid, (middle, keep), span, tol)?;
            set_aside.push((
                indices_k.clone(),
                edge.clone(),
                Wedge {
                    solid: half,
                    additive: wedge.additive,
                    edges: wedge.edges.clone(),
                    melt_refusal: wedge.melt_refusal,
                    credit: Vec::new(),
                },
            ));
        }
        halves.push((index, (middle, toward)));
    }
    for junction in crossed {
        // Each side's band: the one wedge holding it, or of a band cut in
        // two, the half ending at this junction's vertex.
        let mut sides = [None, None];
        for (slot, side) in sides.iter_mut().enumerate() {
            // A band not yet joined into another corner holds no other.
            let held = holding(set_aside, junction.pair[slot])
                .into_iter()
                .filter(|&k| set_aside[k].0.len() == 1);
            let mut here = Vec::new();
            for k in held {
                if crate::junction::cap_at(
                    model,
                    &set_aside[k].2.solid,
                    junction.at,
                    junction.outward[slot],
                    tol,
                )?
                .is_some()
                {
                    here.push(k);
                }
            }
            *side = match here.as_slice() {
                [] => None,
                [k] => Some(*k),
                _ => ogeom_bail!(
                    Construction,
                    "a crossed junction's band is more than one cut wedge in its round"
                ),
            };
        }
        // Built in the same round, or not yet.
        let (a, b) = match sides {
            [None, None] => continue,
            [Some(a), Some(b)] if a != b => (a, b),
            _ => ogeom_bail!(
                Construction,
                "a crossed junction's bands are not one cut wedge each in their round"
            ),
        };
        if set_aside[a].2.additive || set_aside[b].2.additive {
            ogeom_bail!(Construction, "a crossed junction's band is a fill");
        }
        if !set_aside[a].2.credit.is_empty() || !set_aside[b].2.credit.is_empty() {
            ogeom_bail!(Construction, "a crossed junction's band is already joined");
        }
        let (taken, of_second) = crate::junction::corner(
            model,
            [&set_aside[a].2.solid, &set_aside[b].2.solid],
            junction,
            solid,
            law,
            tol,
        )?;
        let second = set_aside[b].1.clone();
        let later = set_aside.remove(a.max(b));
        let earlier = set_aside.remove(a.min(b));
        let ((mut indices, edge, first), (others, other, last)) = if a < b {
            (earlier, later)
        } else {
            (later, earlier)
        };
        indices.extend(others);
        let mut rounded = first.edges;
        rounded.extend(last.edges);
        rounded.push(other);
        set_aside.push((
            indices,
            edge,
            Wedge {
                solid: taken,
                additive: false,
                edges: rounded,
                melt_refusal: first.melt_refusal.or(last.melt_refusal),
                credit: of_second
                    .into_iter()
                    .map(|face| (face, second.clone()))
                    .collect(),
            },
        ));
    }
    Ok(halves)
}

/// The wedges a round set aside, those of two bands meeting in a
/// junction's plane joined into one across it, each wedge with the edge it
/// rounds.
pub(crate) fn join_flush(
    model: &mut Model,
    set_aside: Vec<(Vec<usize>, Shape, Wedge)>,
    junctions: &[Junction],
    halves: &[(usize, (Point, Vector))],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Wedge)>> {
    // Each entry with every edge of the request whose band it holds.
    let mut held: Vec<Option<(Vec<usize>, Shape, Wedge)>> =
        set_aside.into_iter().map(Some).collect();
    // The flush junctions, and the planes a band crossed at both ends was
    // cut in, each as the plane, a point of it the meeting stands at, and
    // the edges on either side.
    let meetings: Vec<((Point, Vector), Point, [usize; 2])> = junctions
        .iter()
        .filter(|j| !j.crossed)
        .map(|j| (j.flush, j.at, j.pair))
        .chain(
            halves
                .iter()
                .map(|&(index, (on, normal))| ((on, normal), on, [index, index])),
        )
        .collect();
    for ((on, normal), at, pair) in meetings {
        // The wedge holding each side's piece that ends in the plane, the
        // two sides two wedges.
        let mut sides = [None, None];
        for slot in 0..2 {
            for (k, entry) in held.iter().enumerate() {
                let Some((indices, _, wedge)) = entry else {
                    continue;
                };
                if !indices.contains(&pair[slot]) || (slot == 1 && sides[0] == Some(k)) {
                    continue;
                }
                let gap = crate::marched::band_fit_target(tol) * 2.0;
                if crate::junction::ends_in(model, &wedge.solid, (on, normal), at, gap, tol)? {
                    sides[slot] = Some(k);
                    if slot == 0 {
                        break;
                    }
                }
            }
        }
        let [Some(a), Some(b)] = sides else {
            continue;
        };
        if a == b {
            // A loop closing on itself: its wedge's two ends in the plane
            // are sewn to each other.
            let Some((indices, edge, ring)) = held[a].take() else {
                continue;
            };
            let (solid, traced) = crate::junction::joined(
                model,
                [&ring.solid, &ring.solid],
                (on, normal),
                at,
                crate::marched::band_fit_target(tol) * 2.0,
                tol,
            )?;
            let mut credit = Vec::new();
            for (given, image) in &traced {
                let owner = ring
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(given))
                    .map_or_else(|| edge.clone(), |(_, e)| e.clone());
                credit.push((image.clone(), owner));
            }
            held[a] = Some((
                indices,
                edge,
                Wedge {
                    solid,
                    credit,
                    ..ring
                },
            ));
            continue;
        }
        let (Some((mut indices, edge, first)), Some((others, other, second))) =
            (held[a].take(), held[b].take())
        else {
            continue;
        };
        let (solid, traced) = crate::junction::joined(
            model,
            [&first.solid, &second.solid],
            (on, normal),
            at,
            crate::marched::band_fit_target(tol) * 2.0,
            tol,
        )?;
        // Each face of the joined wedge stays credited to the edge whose
        // band it came from.
        let mut credit: Vec<(Shape, Shape)> = Vec::new();
        for (from, rounded) in [(&first, &edge), (&second, &other)] {
            for face in ogeom_topo::explore_unique(model, &from.solid, ShapeType::Face)? {
                let owner = from
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or_else(|| rounded.clone(), |(_, e)| e.clone());
                for (given, image) in &traced {
                    if given.is_same(&face) {
                        credit.push((image.clone(), owner.clone()));
                    }
                }
            }
        }
        indices.extend(others);
        indices.dedup();
        let mut rounded: Vec<Shape> = first.edges;
        rounded.extend(second.edges);
        if !other.is_same(&edge) {
            rounded.push(other);
        }
        held[a] = Some((
            indices,
            edge,
            Wedge {
                solid,
                additive: false,
                edges: rounded,
                melt_refusal: first.melt_refusal.or(second.melt_refusal),
                credit,
            },
        ));
    }
    Ok(held
        .into_iter()
        .flatten()
        .map(|(_, edge, wedge)| (edge, wedge))
        .collect())
}
