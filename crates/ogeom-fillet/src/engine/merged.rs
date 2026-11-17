//! Faces a blend leaves in pieces on one carrier, merged after it.
//!
//! A route that applies its strips by a boolean, or closes a corner over
//! a cell, can leave two faces sharing an edge and one carrier facing one
//! way: one surface standing as two faces with a line between them that no
//! feature of the part draws. Mainstream CAD shows no such line after a
//! blend. The pieces the blend made are merged with their neighbours on
//! one carrier ([`ogeom_heal::unify_same_domain_around`]), and the merge's
//! history follows the blend's, so every input element is told of the face
//! that now holds it. Splits the input already held between faces the
//! blend left alone stay as they were.

use ogeom_algo::Built;
use ogeom_core::{FastSet, OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

/// `built`, a blend of `input`, with the faces it made merged with their
/// neighbours on one carrier where any two such share an edge. `built` as
/// it is where none do, or where the merge does not apply (a cluster whose
/// boundary does not re-chain is left split by the merge itself).
///
/// # Errors
///
/// [`OgeomError::Cancelled`] when cancelled, and as reading the solids.
pub(crate) fn merged(
    model: &mut Model,
    input: &Shape,
    built: Built,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let before: FastSet<_> = explore_unique(model, input, ShapeType::Face)?
        .iter()
        .map(Shape::node)
        .collect();
    let made: Vec<Shape> = explore_unique(model, &built.shape, ShapeType::Face)?
        .into_iter()
        .filter(|f| !before.contains(&f.node()))
        .collect();
    if made.is_empty() || !any_split(model, &built.shape, &made, tol)? {
        return Ok(built);
    }
    let unified = match ogeom_heal::unify_same_domain_around(model, &built.shape, &made, tol) {
        Ok((_, 0)) => return Ok(built),
        Ok((unified, _)) => unified,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(built),
    };
    if !ogeom_algo::check(model, &unified.shape, tol)?.is_valid() {
        return Ok(built);
    }
    Ok(Built::new(
        unified.shape,
        built.history.then(&unified.history),
    ))
}

/// Whether a face of `made` shares an edge of `solid` with another face on
/// its carrier facing its way.
fn any_split(model: &Model, solid: &Shape, made: &[Shape], tol: Tolerances) -> OgeomResult<bool> {
    let made_nodes: FastSet<_> = made.iter().map(Shape::node).collect();
    for face in made {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            let faces = ancestors_of(model, solid, &edge, ShapeType::Face)?;
            let [a, b] = faces.as_slice() else {
                continue;
            };
            if a.node() == b.node() {
                continue;
            }
            if !(made_nodes.contains(&a.node()) || made_nodes.contains(&b.node())) {
                continue;
            }
            if ogeom_heal::same_domain(model, a, b, tol)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
