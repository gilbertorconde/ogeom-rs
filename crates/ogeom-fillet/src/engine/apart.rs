//! Rounds whose changes stand apart, each group built as it is alone (ADR
//! 0001, amendment "rounds standing apart"). The edges asked fall into
//! groups joined at shared vertices. The engine rounds each group alone on
//! the solid as given. What a group's round takes off, and what it adds,
//! are its tools, the faces its bands leave named for their edges. Where
//! every tool stands clear of every other group's, no band meets another
//! group's band or the faces it changes, so the rounds asked together are
//! the rounds alone. Where two groups' cuts meet (convex bands crossing
//! though their edges share no vertex), the cuts take their union, as
//! convex creases meeting across a convex turn do, provided the request
//! passes the overlap guards: no two bands overlap across a face they
//! share (amendment "convex bands crossing apart"). The tools are applied
//! as the one-boolean fallback's.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};

/// `once` rounded on `solid` where its edges fall into at least two groups
/// joined at shared vertices, the engine rounds each group alone, and what
/// each group's round takes off or adds stands clear of every other
/// group's, or meets it only where both take material off and no two
/// bands overlap across a face they share; `None` otherwise, or where a
/// boolean refuses or the result is not valid.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled, and reading the
/// edges.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let groups = joined(model, once)?;
    if groups.len() < 2 {
        return Ok(None);
    }
    let mut scratch = model.clone();
    // Each group's tools, as (group, tool).
    let mut tools: Vec<(usize, Tool)> = Vec::new();
    for (g, group) in groups.iter().enumerate() {
        let asked: Vec<(Shape, f64)> = group.iter().map(|&i| once[i].clone()).collect();
        let alone = match round(&mut scratch, solid, &asked, tol) {
            Ok(Some(alone)) => alone,
            Err(e @ OgeomError::Cancelled) => return Err(e),
            _ => return Ok(None),
        };
        let (Ok(taken), Ok(added)) = (
            ogeom_bool::cut(&mut scratch, solid, &alone.shape, tol),
            ogeom_bool::cut(&mut scratch, &alone.shape, solid, tol),
        ) else {
            return Ok(None);
        };
        let mut any = false;
        for (sense, change) in [(Sense::Cut, &taken), (Sense::Fuse, &added)] {
            if explore_unique(&scratch, &change.shape, ShapeType::Solid)?.is_empty() {
                continue;
            }
            let mut faces = Vec::new();
            for (edge, _) in &asked {
                for band in alone.history.generated(edge) {
                    if scratch.kind_of(band)? != ShapeType::Face {
                        continue;
                    }
                    for face in change.history.trace(band) {
                        faces.push((face.clone(), edge.clone()));
                    }
                }
            }
            // A change no band of the group's bounds is not a round's.
            if faces.is_empty() {
                return Ok(None);
            }
            any = true;
            tools.push((
                g,
                Tool {
                    solid: change.shape.clone(),
                    sense,
                    faces,
                },
            ));
        }
        if !any {
            return Ok(None);
        }
    }
    // Every tool clear of every other group's: no band meets another
    // group's, nor a face another group's round changes. Two cuts may
    // meet: convex bands crossing take their cuts' union, where no two
    // bands overlap across a face they share.
    let mut cuts_meet = false;
    for (i, (gi, ti)) in tools.iter().enumerate() {
        for (gj, tj) in &tools[..i] {
            if gi == gj {
                continue;
            }
            let gap = ogeom_algo::distance_between_shapes(
                &scratch,
                &ti.solid,
                &tj.solid,
                ogeom_intersect::ExtremaOptions::default(),
                tol,
            )?;
            if gap.distance <= tol.confusion() {
                if ti.sense != Sense::Cut || tj.sense != Sense::Cut {
                    return Ok(None);
                }
                cuts_meet = true;
            }
        }
    }
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    if cuts_meet {
        let radii: Vec<f64> = once.iter().map(|(_, r)| *r).collect();
        match super::fits::request_fits(model, solid, &edges, &radii, tol) {
            Ok(_) => {}
            Err(blamed) if matches!(blamed.error, OgeomError::Cancelled) => {
                return Err(blamed.error);
            }
            Err(_) => return Ok(None),
        }
    }
    let tools: Vec<Tool> = tools.into_iter().map(|(_, t)| t).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(&scratch, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    *model = scratch;
    Ok(Some(built))
}

/// The edges of `once` gathered into groups, two edges in one group where a
/// chain of edges sharing vertices joins them, each group's indices in
/// request order. Empty where an edge has no vertices.
fn joined(model: &Model, once: &[(Shape, f64)]) -> OgeomResult<Vec<Vec<usize>>> {
    let mut ends = Vec::with_capacity(once.len());
    for (edge, _) in once {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(Vec::new());
        };
        ends.push([a, b]);
    }
    let meet = |i: usize, j: usize| ends[i].iter().any(|v| ends[j].iter().any(|w| w.is_same(v)));
    let mut group: Vec<usize> = (0..once.len()).collect();
    let root = |group: &[usize], mut i: usize| {
        while group[i] != i {
            i = group[i];
        }
        i
    };
    for i in 0..once.len() {
        for j in 0..i {
            if meet(i, j) {
                let (ri, rj) = (root(&group, i), root(&group, j));
                group[ri.max(rj)] = ri.min(rj);
            }
        }
    }
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut slot: Vec<Option<usize>> = vec![None; once.len()];
    for i in 0..once.len() {
        let r = root(&group, i);
        match slot[r] {
            Some(k) => out[k].push(i),
            None => {
                slot[r] = Some(out.len());
                out.push(vec![i]);
            }
        }
    }
    Ok(out)
}
