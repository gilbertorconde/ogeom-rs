//! Thickening an open mesh into a closed one: two offset layers of it
//! joined by a wall along each rim.
//!
//! The sheet's faces face the way their corners wind counter-clockwise.
//! Each layer is the sheet offset by [`Editable::offset`], folds repaired:
//! the front layer, which keeps the sheet's winding, by the larger
//! distance, and the back layer, wound the other way, by the smaller. The
//! layers' rims are the sheet's, less any rim vertex a fold repair merged
//! along the rim, so each rim of one layer is walked beside the same rim of
//! the other, by where its vertices came from on the sheet's rim, and the
//! strip between them is closed with triangles, each joining an edge of one
//! layer to a vertex of the other.
//!
//! The result is checked exactly as [`boolean`](crate::boolean()) checks
//! an operand, so layers or walls that run into each other, where the sheet
//! bends tighter than the distance or its rim nearly touches itself, are
//! refused rather than returned.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_topo::Triangulation;

use crate::boolean::check_solid;
use crate::editable::{Editable, OffsetReport};

/// Which side of a sheet a thickening grows to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Thicken {
    /// The side the sheet's faces face away from: the sheet is the front.
    Inward,
    /// The side the sheet's faces face: the sheet is the back.
    Outward,
    /// Half to each side.
    Both,
}

/// A thickened sheet.
#[derive(Debug, Clone, PartialEq)]
pub struct Thickened {
    /// The solid's boundary: the front layer, then the back layer, then the
    /// walls; closed, 2-manifold and wound outward. It carries no normals
    /// and no parameters.
    pub mesh: Triangulation,
    /// What the offset did to the front layer, the one the sheet's faces
    /// face out of the solid on.
    pub front: OffsetReport,
    /// What it did to the back layer.
    pub back: OffsetReport,
}

const NOTHING: OffsetReport = OffsetReport {
    merged: 0,
    eased: 0,
    capped: 0,
    folded: 0,
};

/// Thicken an open mesh by `distance` to `side`.
///
/// The sheet must be 2-manifold and have a rim; its vertices are taken as
/// numbered, not welded.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `distance` is not a positive distance, the sheet is not 2-manifold or is
/// closed, or the result does not bound a solid as
/// [`boolean`](crate::boolean()) asks of an operand, named "thickened";
/// [`OgeomError::NullObject`](ogeom_core::OgeomError::NullObject) if it has
/// no triangles; and as [`Editable::from`] and [`Editable::offset`].
pub fn thicken(mesh: &Triangulation, distance: f64, side: Thicken) -> OgeomResult<Thickened> {
    if !(distance.is_finite() && distance > 0.0) {
        ogeom_bail!(Construction, "a thickness of {distance} is not a distance");
    }
    if mesh.triangles.is_empty() {
        ogeom_bail!(NullObject, "the sheet has no triangles");
    }
    let sheet = Editable::from(mesh, 0.0)?;
    if let Some(v) = sheet.non_manifold().first() {
        ogeom_bail!(
            Construction,
            "the sheet is not 2-manifold at vertex {}",
            sheet.source_vertex(*v).unwrap_or(*v)
        );
    }
    let rims = sheet.rims();
    if rims.is_empty() {
        ogeom_bail!(
            Construction,
            "the sheet is closed: it bounds a solid already, to hollow rather than thicken"
        );
    }
    let (back, front) = match side {
        Thicken::Inward => (-distance, 0.0),
        Thicken::Outward => (0.0, distance),
        Thicken::Both => (-0.5 * distance, 0.5 * distance),
    };
    let layer = |d: f64| -> OgeomResult<(Editable, OffsetReport)> {
        let mut e = sheet.clone();
        let report = if d == 0.0 {
            NOTHING
        } else {
            e.offset(d, None)?
        };
        Ok((e, report))
    };
    let (front_layer, front_report) = layer(front)?;
    let (back_layer, back_report) = layer(back)?;

    // Each sheet rim vertex's rim and place along it.
    let mut place = vec![None; sheet.vertex_count()];
    let mut lengths = Vec::with_capacity(rims.len());
    for (r, rim) in rims.iter().enumerate() {
        for (at, &h) in rim.edges.iter().enumerate() {
            if let Some([v, _]) = sheet.edge(h) {
                place[v as usize] = Some((r, at));
            }
        }
        lengths.push(rim.edges.len());
    }

    let mut out = Triangulation::new();
    let front_numbers = take(&front_layer, &mut out, false);
    let back_numbers = take(&back_layer, &mut out, true);
    let front_rims = along(&front_layer, &front_numbers, &place, rims.len())?;
    let back_rims = along(&back_layer, &back_numbers, &place, rims.len())?;
    for r in 0..rims.len() {
        zip(&mut out, &back_rims[r], &front_rims[r], lengths[r])?;
    }
    out.deflection_met = false;
    check_solid(&out, "thickened", 0)?;
    Ok(Thickened {
        mesh: out,
        front: front_report,
        back: back_report,
    })
}

/// Add a layer's faces to `out`, wound the other way when `flip`; returns
/// each vertex's number in `out`.
fn take(e: &Editable, out: &mut Triangulation, flip: bool) -> Vec<Option<u32>> {
    let mut number = vec![None; e.vertex_count()];
    for f in 0..u32::try_from(e.face_count()).unwrap_or(u32::MAX) {
        let Some(corners) = e.face(f) else {
            continue;
        };
        let mapped = corners.map(|v| {
            *number[v as usize].get_or_insert_with(|| {
                let n = u32::try_from(out.positions.len()).unwrap_or(u32::MAX);
                out.positions.push(e.position(v).unwrap_or_default());
                n
            })
        });
        out.triangles.push(if flip {
            [mapped[0], mapped[2], mapped[1]]
        } else {
            mapped
        });
    }
    number
}

/// Each sheet rim's vertices on a layer, as their places along the sheet's
/// rim and their numbers in the output, in the rim's order from the lowest
/// place.
fn along(
    e: &Editable,
    number: &[Option<u32>],
    place: &[Option<(usize, usize)>],
    count: usize,
) -> OgeomResult<Vec<Vec<(usize, u32)>>> {
    let mut found: Vec<Vec<(usize, u32)>> = vec![Vec::new(); count];
    for rim in e.rims() {
        let mut walk = Vec::with_capacity(rim.edges.len());
        let mut which = None;
        for &h in &rim.edges {
            let Some([v, _]) = e.edge(h) else {
                continue;
            };
            let at = e
                .source_vertex(v)
                .and_then(|s| place.get(s as usize).copied().flatten());
            let (Some((r, at)), Some(n)) = (at, number[v as usize]) else {
                ogeom_bail!(
                    Invariant,
                    "an offset layer's rim has a vertex off the sheet's rim"
                );
            };
            if which.is_some_and(|w| w != r) {
                ogeom_bail!(
                    Invariant,
                    "an offset layer's rim joins two of the sheet's rims"
                );
            }
            which = Some(r);
            walk.push((at, n));
        }
        let Some(r) = which else { continue };
        if !found[r].is_empty() {
            ogeom_bail!(
                Invariant,
                "an offset layer parted one of the sheet's rims in two"
            );
        }
        let start = walk
            .iter()
            .enumerate()
            .min_by_key(|(_, (at, _))| *at)
            .map_or(0, |(k, _)| k);
        walk.rotate_left(start);
        if walk.windows(2).any(|w| w[0].0 >= w[1].0) {
            ogeom_bail!(
                Invariant,
                "an offset layer's rim does not keep the order of the sheet's"
            );
        }
        found[r] = walk;
    }
    if found.iter().any(Vec::is_empty) {
        ogeom_bail!(Invariant, "an offset layer lost one of the sheet's rims");
    }
    Ok(found)
}

/// Close the strip between a rim of the back layer and the same rim of the
/// front, each given by its vertices' places along the sheet's rim of
/// `length` vertices. Each triangle takes the next edge of whichever layer
/// reaches the lesser place next.
fn zip(
    out: &mut Triangulation,
    back: &[(usize, u32)],
    front: &[(usize, u32)],
    length: usize,
) -> OgeomResult<()> {
    let (nb, nf) = (back.len(), front.len());
    if nb < 2 || nf < 2 {
        ogeom_bail!(
            Invariant,
            "an offset layer's rim has fewer than two vertices"
        );
    }
    // The place of the `k`th vertex along, counted past the start.
    let at = |side: &[(usize, u32)], k: usize| {
        let n = side.len();
        side[k % n].0 + length * (k / n)
    };
    let vertex = |side: &[(usize, u32)], k: usize| side[k % side.len()].1;
    let (mut i, mut j) = (0, 0);
    while i < nb || j < nf {
        let take_back = j == nf || (i < nb && at(back, i + 1) <= at(front, j + 1));
        if take_back {
            out.triangles
                .push([vertex(back, i), vertex(back, i + 1), vertex(front, j)]);
            i += 1;
        } else {
            out.triangles
                .push([vertex(back, i), vertex(front, j + 1), vertex(front, j)]);
            j += 1;
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "thicken/tests.rs"]
mod tests;
