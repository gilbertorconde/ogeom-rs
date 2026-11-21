//! The non-local rounds' general tool (ADR 0001, §4.5, M5c): the far
//! feature the rounds reach taken off a copy of the solid, the rounds built
//! on that copy as the local request they are there, and the tools the
//! copy less its rounding (what the rounds take) and the rounding less the
//! copy (what they fill), applied by booleans on the solid as given. The
//! rounds keep their own surfaces and the feature they reach is changed,
//! as where a round's strip or corner the other tools close is not one
//! they close: a corner of mixed convexity, strips meeting at an angle on
//! a curved solid, a strip along one closed edge, one along a cone.

use ogeom_algo::Built;
use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

use super::fallback::{Sense, Tool};
use super::{Id, id};

/// How a solid is rounded where nothing else reaches the rounds: the
/// engine's own entry.
pub(crate) type Rounder =
    fn(&mut Model, &Shape, &[(Shape, f64)], Tolerances) -> OgeomResult<Option<Built>>;

/// The feature of `solid` holding `far`: those faces, and every face that
/// meets one along its outer ring, but those in `near` (the faces the
/// rounds end at) and those it stands in through a hole of theirs (the
/// faces it is cut into or stands on). `None` where it holds every face.
fn feature(
    model: &Model,
    solid: &Shape,
    far: &[Shape],
    near: &[Shape],
) -> OgeomResult<Option<Vec<Shape>>> {
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    // Each edge's faces, with whether it stands in a hole of that face.
    let mut uses: FastMap<Id, Vec<(usize, bool)>> = FastMap::default();
    for (fi, f) in faces.iter().enumerate() {
        for (k, wire) in model.children_of(f)?.iter().enumerate() {
            for e in model.children_of(wire)? {
                uses.entry(id(&e)).or_default().push((fi, k > 0));
            }
        }
    }
    let index = |s: &Shape| faces.iter().position(|f| f.is_same(s));
    let near: Vec<usize> = near.iter().filter_map(index).collect();
    let mut taken: Vec<usize> = far.iter().filter_map(index).collect();
    let mut k = 0;
    while k < taken.len() {
        let f = taken[k];
        for wire in model.children_of(&faces[f])? {
            for e in model.children_of(&wire)? {
                for (g, hole) in uses.get(&id(&e)).into_iter().flatten() {
                    if *g != f && !hole && !near.contains(g) && !taken.contains(g) {
                        taken.push(*g);
                    }
                }
            }
        }
        k += 1;
    }
    if taken.is_empty() || taken.len() >= faces.len() {
        return Ok(None);
    }
    taken.sort_unstable();
    Ok(Some(taken.into_iter().map(|f| faces[f].clone()).collect()))
}

/// The tools of `rounds` on `solid` where they reach the faces `far` and
/// end at the faces `near`: the feature holding `far` removed from a copy
/// ([`ogeom_bool::remove_faces`]), the rounds built on the copy by `round`,
/// and the copy less its rounding cut, the rounding less the copy fused,
/// each strip's faces named for its edge. `None` where the feature will not
/// come off, the copy does not round, or a boolean refuses.
///
/// # Errors
///
/// As rounding the copy.
pub(crate) fn tools(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    far: &[Shape],
    near: &[Shape],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Tool>>> {
    tools_across(model, solid, rounds, far, near, &[], round, tol)
}

/// [`tools`], each edge of `rounds` named in `across` with two faces
/// rounded on the copy as the edge where those two faces meet there (an
/// edge of a face the feature takes, which the copy no longer has, standing
/// for the crease its neighbours close in). `None` also where those faces
/// do not meet on the copy along one edge.
///
/// # Errors
///
/// As [`tools`].
#[allow(clippy::too_many_arguments, reason = "the tools' whole reading")]
pub(crate) fn tools_across(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    far: &[Shape],
    near: &[Shape],
    across: &[(Shape, [Shape; 2])],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Tool>>> {
    let Some(taken) = feature(model, solid, far, near)? else {
        return Ok(None);
    };
    let Ok(copy) = ogeom_bool::remove_faces(model, solid, &taken, tol) else {
        return Ok(None);
    };
    let mut traced: Vec<(Shape, f64)> = Vec::with_capacity(rounds.len());
    for (edge, radius) in rounds {
        let image = match across.iter().find(|(e, _)| e.is_same(edge)) {
            Some((_, [a, b])) => {
                let ([a], [b]) = (copy.history.trace(a), copy.history.trace(b)) else {
                    return Ok(None);
                };
                let edges_of = |f: &Shape| explore_unique(model, f, ShapeType::Edge);
                let of_b = edges_of(b)?;
                let shared: Vec<Shape> = edges_of(a)?
                    .into_iter()
                    .filter(|e| of_b.iter().any(|x| x.is_same(e)))
                    .collect();
                let [only] = &shared[..] else {
                    return Ok(None);
                };
                only.clone()
            }
            None => {
                let Some(image) = image_of(model, solid, &copy, edge, tol)? else {
                    return Ok(None);
                };
                image
            }
        };
        // Two edges standing for one crease are rounded once, alike.
        if let Some((_, r)) = traced.iter().find(|(e, _)| e.is_same(&image))
            && (r - radius).abs() > tol.confusion()
        {
            return Ok(None);
        }
        traced.push((image, *radius));
    }
    let mut once: Vec<(Shape, f64)> = Vec::with_capacity(traced.len());
    for (e, r) in &traced {
        if !once.iter().any(|(x, _)| x.is_same(e)) {
            once.push((e.clone(), *r));
        }
    }
    let Some(rounded) = round(model, &copy.shape, &once, tol)? else {
        return Ok(None);
    };
    let mut tools = Vec::with_capacity(2);
    for (sense, (from, less)) in [
        (Sense::Fuse, (&rounded.shape, &copy.shape)),
        (Sense::Cut, (&copy.shape, &rounded.shape)),
    ] {
        let Ok(part) = ogeom_bool::cut(model, from, less, tol) else {
            return Ok(None);
        };
        if explore_unique(model, &part.shape, ShapeType::Solid)?.is_empty() {
            continue;
        }
        let kept = explore_unique(model, &part.shape, ShapeType::Face)?;
        let mut named = Vec::new();
        for ((edge, _), (image, _)) in rounds.iter().zip(&traced) {
            for blend in rounded.history.generated(image) {
                let images = part
                    .history
                    .modified(blend)
                    .iter()
                    .chain(part.history.generated(blend))
                    .cloned()
                    .chain(std::iter::once(blend.clone()));
                for face in images {
                    if let Some(k) = kept.iter().find(|k| k.is_same(&face)) {
                        named.push((k.clone(), edge.clone()));
                    }
                }
            }
        }
        tools.push(Tool {
            solid: part.shape,
            sense,
            faces: named,
        });
    }
    Ok((!tools.is_empty()).then_some(tools))
}

/// The edge of `copy` standing for `edge` of `solid`. Taking faces off
/// records what became of faces only: an edge between two rebuilt faces
/// traces to itself though the copy holds a new one, which is found as
/// the edge the two faces' images share that runs through `edge`'s
/// middle. `None` where no one edge does.
fn image_of(
    model: &Model,
    solid: &Shape,
    copy: &Built,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let edges = explore_unique(model, &copy.shape, ShapeType::Edge)?;
    if let [only] = copy.history.trace(edge)
        && edges.iter().any(|e| e.is_same(only))
    {
        return Ok(Some(only.clone()));
    }
    let faces = ancestors_of(model, solid, edge, ShapeType::Face)?;
    let [a, b] = faces.as_slice() else {
        return Ok(None);
    };
    let ([a], [b]) = (copy.history.trace(a), copy.history.trace(b)) else {
        return Ok(None);
    };
    let Some((curve, (s, t))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    let middle = curve.point_at(0.5 * (s + t), tol)?;
    let of_b = explore_unique(model, b, ShapeType::Edge)?;
    let mut through = Vec::new();
    for e in explore_unique(model, a, ShapeType::Edge)? {
        if of_b.iter().any(|x| x.is_same(&e)) && runs_through(model, &e, middle, tol)? {
            through.push(e);
        }
    }
    let [only] = &through[..] else {
        return Ok(None);
    };
    Ok(Some(only.clone()))
}

/// Whether `edge` passes within the confusion distance of `p` inside its
/// own range.
fn runs_through(model: &Model, edge: &Shape, p: Point, tol: Tolerances) -> OgeomResult<bool> {
    let Some((curve, (s, t))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(false);
    };
    let near = ogeom_algo::project_on_curve(&curve, p, 64, tol)?;
    if near.distance > tol.confusion() {
        return Ok(false);
    }
    let (lo, hi) = (s.min(t), s.max(t));
    let within = |u: f64| u >= lo - tol.parametric() && u <= hi + tol.parametric();
    let (a, b) = curve.domain();
    let period = b - a;
    Ok(within(near.parameter)
        || (curve.is_periodic()
            && (within(near.parameter + period) || within(near.parameter - period))))
}
