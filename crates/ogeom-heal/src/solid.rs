//! Solids from closed shells, outer boundaries told from voids by which
//! shell lies inside which.
//!
//! *Elsewhere:* `ShapeFix_Solid::SolidFromShell`.
//!
//! A file or a sewing hands over closed shells with nothing saying which
//! bounds a part and which a cavity in one. Nesting says it: a shell inside
//! no other is a part's outer boundary, a shell inside one is a void of
//! it, a shell inside two is the outer boundary of a part inside that
//! cavity, and so on down.

use ogeom_algo::{Built, Containment, History, classify_in_solid_exact, is_shell_closed};
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// The solids `shells` bound: one where a single shell is outermost, else
/// a compound of them, one per outermost shell and per shell lying in a
/// void.
///
/// Each shell is placed against each other shell by its vertices that lie
/// off it, classified by ray parity, which needs no shell to face either
/// way; they must all read inside or all outside. Shells whose faces cross
/// with no vertex of either inside the other are not told apart from
/// shells side by side. A shell's parent is the deepest shell that
/// holds it. A shell held by an even number of shells is an outer
/// boundary, and the shells it is the parent of are its voids. Each solid
/// is put together by [`make_solid`](ogeom_algo::make_solid), which turns
/// an outer shell to face out of the volume it bounds and a void's to face
/// into its cavity.
///
/// # Errors
///
/// [`OgeomError::Construction`] if
/// `shells` is empty, holds something that is not a closed shell, or two
/// shells cross (vertices of one on both sides of the other, or each
/// holding the other) or one has every vertex on the other; as
/// classification otherwise.
pub fn solids_from_shells(
    model: &mut Model,
    shells: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
    if shells.is_empty() {
        ogeom_bail!(Construction, "there are no shells to make solids of");
    }
    for shell in shells {
        if model.kind_of(shell)? != ShapeType::Shell || !is_shell_closed(model, shell)? {
            ogeom_bail!(
                Construction,
                "a solid is bounded by closed shells, and this is not one"
            );
        }
    }
    // holds[i][j]: shell i holds shell j.
    let n = shells.len();
    let mut holds = vec![vec![false; n]; n];
    for j in 0..n {
        for i in (0..n).filter(|&i| i != j) {
            holds[i][j] = holds_shell(model, &shells[i], &shells[j], tol)?;
        }
    }
    if (0..n).any(|i| (0..i).any(|j| holds[i][j] && holds[j][i])) {
        ogeom_bail!(Construction, "two shells cross: each holds the other");
    }
    let depth: Vec<usize> = (0..n)
        .map(|j| (0..n).filter(|&i| holds[i][j]).count())
        .collect();
    let parent: Vec<Option<usize>> = (0..n)
        .map(|j| (0..n).filter(|&i| holds[i][j]).max_by_key(|&i| depth[i]))
        .collect();
    let mut solids = Vec::new();
    let mut history = History::new();
    for outer in (0..n).filter(|&i| depth[i].is_multiple_of(2)) {
        let mut bounds = vec![shells[outer].clone()];
        bounds.extend(
            (0..n)
                .filter(|&j| !depth[j].is_multiple_of(2) && parent[j] == Some(outer))
                .map(|j| shells[j].clone()),
        );
        let solid = ogeom_algo::make_solid(model, &bounds)?.shape;
        for shell in &bounds {
            history.generate(shell, solid.clone());
        }
        solids.push(solid);
    }
    if let [solid] = &solids[..] {
        return Ok(Built::new(solid.clone(), history));
    }
    let compound = ogeom_algo::make_compound(model, &solids)?.shape;
    Ok(Built::new(compound, history))
}

/// Whether `inner` lies inside the volume `outer` encloses, read at every
/// vertex of `inner` that is not on `outer`: all must agree.
fn holds_shell(model: &Model, outer: &Shape, inner: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    let mut read = None;
    for vertex in explore_unique(model, inner, ShapeType::Vertex)? {
        let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        let point = vertex.transform(model.datums())?.apply(data.point);
        let inside = match classify_in_solid_exact(model, outer, point, tol) {
            Ok(Containment::In) => true,
            Ok(Containment::Out) => false,
            Ok(Containment::On) | Err(OgeomError::NotDone(_)) => continue,
            Err(e) => return Err(e),
        };
        if read.is_some_and(|was| was != inside) {
            ogeom_bail!(
                Construction,
                "two shells cross: one has vertices inside the other and outside it"
            );
        }
        read = Some(inside);
    }
    let Some(inside) = read else {
        ogeom_bail!(
            Construction,
            "every vertex of one shell lies on another: neither holds the other"
        );
    };
    Ok(inside)
}
