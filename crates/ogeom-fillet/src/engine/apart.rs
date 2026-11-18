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
//!
//! What a group's round changes is read by booleans between the solid and
//! its rounded copy, which share every face but the round's. Where such a
//! boolean refuses, the groups are rounded in turn instead, each on the
//! result of those before. Each later group's edges must come through the
//! earlier rounds unchanged, and the faces the earlier rounds made must
//! come through the later ones untouched; then no group's band met another
//! group's, and the rounds in turn are the rounds alone.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};

/// `once` rounded on `solid` where its edges fall into at least two groups
/// joined at shared vertices, the engine rounds each group alone, and what
/// each group's round takes off or adds stands clear of every other
/// group's, or meets it only where both take material off and no two
/// bands overlap across a face they share; `None` otherwise, or where a
/// boolean refuses or the result is not valid. Where the booleans reading
/// a group's changes refuse, the groups are rounded in turn ([`in_turn`]).
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
            return in_turn(model, solid, once, &groups, round, tol);
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

/// `once` rounded group by group, each group as the engine rounds it alone
/// on the result of the groups before; `None` where a group is not rounded,
/// a later group's edge does not come through the earlier rounds as one
/// edge on the same ends, a face an earlier round made does not come
/// through a later one untouched (the bands met), or the result is not
/// valid. Groups holding an edge on a marched seat go first, and a group
/// the engine does not round on the earlier rounds' result is rounded
/// first instead, once: the rounds standing apart, their order does not
/// change what each takes off or adds.
fn in_turn(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    groups: &[Vec<usize>],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    // Groups on marched seats first: the march reads the solid as given
    // more surely than one an earlier round rebuilt.
    let mut marched = Vec::with_capacity(groups.len());
    for group in groups {
        let mut on = false;
        for &i in group {
            on |= super::marched_seats::on_marched_seat(model, solid, &once[i].0, tol)?;
        }
        marched.push(on);
    }
    let mut order: Vec<&[usize]> = Vec::with_capacity(groups.len());
    for want in [true, false] {
        order.extend(
            groups
                .iter()
                .zip(&marched)
                .filter(|(_, m)| **m == want)
                .map(|(g, _)| g.as_slice()),
        );
    }
    let mut first: Vec<&[usize]> = Vec::new();
    loop {
        let mut scratch = model.clone();
        match turned(&mut scratch, solid, once, &order, round, tol)? {
            Turn::Built(built) => {
                *model = scratch;
                return Ok(Some(*built));
            }
            Turn::Declined(k) if k > 0 && !first.contains(&order[k]) => {
                let group = order.remove(k);
                first.push(group);
                order.insert(0, group);
            }
            Turn::Declined(_) | Turn::Met => return Ok(None),
        }
    }
}

/// How rounding groups in turn went.
enum Turn {
    Built(Box<Built>),
    /// The engine did not round the group at this place in the order.
    Declined(usize),
    /// Two groups' rounds met, or the result is not valid.
    Met,
}

/// The groups rounded in `order` on `model`, for [`in_turn`].
fn turned(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    order: &[&[usize]],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Turn> {
    let mut shape = solid.clone();
    let mut history = History::identity();
    // The faces the rounds so far made, as they stand now.
    let mut made: Vec<Shape> = Vec::new();
    for (k, group) in order.iter().enumerate() {
        let mut asked: Vec<(Shape, f64)> = Vec::with_capacity(group.len());
        for &i in *group {
            let (edge, radius) = &once[i];
            let [image] = history.trace(edge) else {
                return Ok(Turn::Met);
            };
            if !same_ends(model, edge, image, tol)? {
                return Ok(Turn::Met);
            }
            asked.push((image.clone(), *radius));
        }
        let step = match round(model, &shape, &asked, tol) {
            Ok(Some(step)) => step,
            Err(e @ OgeomError::Cancelled) => return Err(e),
            _ => return Ok(Turn::Declined(k)),
        };
        let mut now = Vec::with_capacity(made.len());
        for face in &made {
            match step.history.trace(face) {
                [image] if untouched(model, face, image, tol)? => now.push(image.clone()),
                _ => return Ok(Turn::Met),
            }
        }
        let mut next = history.then(&step.history);
        for (&i, (image, _)) in group.iter().zip(&asked) {
            let edge = &once[i].0;
            for band in step.history.generated(image) {
                if !next.generated(edge).iter().any(|f| f.is_same(band)) {
                    next.generate(edge, band.clone());
                }
                if model.kind_of(band)? == ShapeType::Face {
                    now.push(band.clone());
                }
            }
            next.delete(edge);
        }
        made = now;
        history = next;
        shape = step.shape;
    }
    if !ogeom_algo::check(model, &shape, tol)?.is_valid() {
        return Ok(Turn::Met);
    }
    Ok(Turn::Built(Box::new(Built::new(shape, history))))
}

/// Whether `image` is `face` as it stood: its edges as many, each running
/// between the points one of the face's did. A band that met the face
/// would have cut a new edge into it.
pub(super) fn untouched(
    model: &Model,
    face: &Shape,
    image: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (was, is) = (
        explore_unique(model, face, ShapeType::Edge)?,
        explore_unique(model, image, ShapeType::Edge)?,
    );
    if was.len() != is.len() {
        return Ok(false);
    }
    for edge in &is {
        let mut found = false;
        for old in &was {
            if same_ends(model, old, edge, tol)? {
                found = true;
                break;
            }
        }
        if !found {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether `image` runs between the points `edge` does, either way round.
pub(super) fn same_ends(
    model: &Model,
    edge: &Shape,
    image: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let ends = |e: &Shape| -> OgeomResult<Option<[ogeom_math::Point; 2]>> {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, e)? else {
            return Ok(None);
        };
        Ok(Some([
            crate::setback::point_of(model, &a)?,
            crate::setback::point_of(model, &b)?,
        ]))
    };
    let (Some([a, b]), Some([c, d])) = (ends(edge)?, ends(image)?) else {
        return Ok(false);
    };
    let near = |p: ogeom_math::Point, q: ogeom_math::Point| p.distance(q) <= tol.confusion();
    Ok((near(a, c) && near(b, d)) || (near(a, d) && near(b, c)))
}

/// The edges of `once` gathered into groups, two edges in one group where a
/// chain of edges sharing vertices joins them, each group's indices in
/// request order. Empty where an edge has no vertices.
pub(super) fn joined(model: &Model, once: &[(Shape, f64)]) -> OgeomResult<Vec<Vec<usize>>> {
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
