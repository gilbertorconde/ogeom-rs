//! A blend on a compound of solids, solid by solid.
//!
//! Mainstream CAD blends each body on its own: a solid of the compound is
//! blended with the requests whose edges (or vertices) it holds, a solid
//! with none comes back as it was, and the results are gathered into a
//! compound again. Nothing is handed to a boolean with the compound as one
//! operand, which would fuse solids that touch.

use ogeom_algo::{Built, History};
use ogeom_core::OgeomResult;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// `build` run on each solid of `shape`, a compound of solids, with the
/// requests of `asked` whose shape (`held`) that solid holds, the other
/// solids kept as they are, the results gathered into a compound again
/// whose history tells of each solid's and of the compound modified into
/// the new one.
///
/// `None` where `shape` is no compound, holds anything but solids or
/// nothing, an asked shape lies in none of its solids, or `build` hands a
/// solid back: the caller's own path then reads the request.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where an asked shape lies in two solids of the compound (each solid is
/// blended on its own, so a shape two share belongs to neither); and as
/// `build`.
pub(crate) fn per_solid<S: Clone>(
    model: &mut Model,
    shape: &Shape,
    asked: &[S],
    held: impl Fn(&S) -> &Shape,
    mut build: impl FnMut(&mut Model, &Shape, &[S]) -> OgeomResult<Option<Built>>,
) -> OgeomResult<Option<Built>> {
    if model.kind_of(shape)? != ShapeType::Compound {
        return Ok(None);
    }
    let parts = model.children_of(shape)?;
    if parts.is_empty()
        || parts
            .iter()
            .any(|p| !matches!(model.kind_of(p), Ok(ShapeType::Solid)))
    {
        return Ok(None);
    }
    let mut own: Vec<Vec<S>> = vec![Vec::new(); parts.len()];
    let mut seen: Vec<(ShapeType, Vec<Vec<Shape>>)> = Vec::new();
    for request in asked {
        let target = held(request);
        let kind = model.kind_of(target)?;
        let k = match seen.iter().position(|(t, _)| *t == kind) {
            Some(k) => k,
            None => {
                let of_kind = parts
                    .iter()
                    .map(|p| explore_unique(model, p, kind))
                    .collect::<OgeomResult<_>>()?;
                seen.push((kind, of_kind));
                seen.len() - 1
            }
        };
        let holders: Vec<usize> = (0..parts.len())
            .filter(|i| seen[k].1[*i].iter().any(|s| s.is_same(target)))
            .collect();
        match holders[..] {
            [] => return Ok(None),
            [i] => own[i].push(request.clone()),
            _ => return Err(shared_refused()),
        }
    }
    let mut out = Vec::with_capacity(parts.len());
    let mut history = History::new();
    for (part, mine) in parts.iter().zip(&own) {
        if mine.is_empty() {
            out.push(part.clone());
            continue;
        }
        let Some(built) = build(model, part, mine)? else {
            return Ok(None);
        };
        history = history.then(&built.history);
        out.push(built.shape);
    }
    let gathered = ogeom_algo::make_compound(model, &out)?.shape;
    history.modify(shape, gathered.clone());
    Ok(Some(Built::new(gathered, history)))
}

/// [`per_solid`] for a construction that builds or refuses, never hands
/// back: `None` only where `shape` is not a compound of solids holding
/// every asked shape.
///
/// # Errors
///
/// As [`per_solid`].
pub(crate) fn each_solid<S: Clone>(
    model: &mut Model,
    shape: &Shape,
    asked: &[S],
    held: impl Fn(&S) -> &Shape,
    mut build: impl FnMut(&mut Model, &Shape, &[S]) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    per_solid(model, shape, asked, held, |model, part, mine| {
        build(model, part, mine).map(Some)
    })
}

/// The refusal of a shape two solids of a compound share.
fn shared_refused() -> ogeom_core::OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "an edge or vertex asked lies on two solids of the compound; each solid is \
         blended on its own, so a shape two of them share belongs to neither"
    )
}
