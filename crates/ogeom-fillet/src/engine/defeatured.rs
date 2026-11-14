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
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

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
    let Some(taken) = feature(model, solid, far, near)? else {
        return Ok(None);
    };
    let Ok(copy) = ogeom_bool::remove_faces(model, solid, &taken, tol) else {
        return Ok(None);
    };
    let mut traced = Vec::with_capacity(rounds.len());
    for (edge, radius) in rounds {
        let [only] = copy.history.trace(edge) else {
            return Ok(None);
        };
        traced.push((only.clone(), *radius));
    }
    let Some(rounded) = round(model, &copy.shape, &traced, tol)? else {
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
