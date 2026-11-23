//! The union of many solids in one call.
//!
//! Each solid lies within its guaranteed bound, tolerances included, so two
//! solids whose bounds widened by the confusion distance do not overlap are
//! more than that distance apart everywhere: no contact the boolean would
//! see. A set of solids whose widened bounds are pairwise apart is then a
//! compound of disjoint lumps, which the boolean takes as one argument. The
//! solids that meet anything are coloured so that no two of one colour
//! have overlapping widened bounds, each colour becomes one such compound,
//! and the colours are fused in a balanced tree. Union is associative and
//! commutative, so the tree bounds the same material as fusing one at a
//! time. A row of solids each touching the next is two colours and one
//! boolean.

use ogeom_algo::{Built, History};
use ogeom_core::{Cause, FastMap, FastSet, OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::Aabb;
use ogeom_topo::{Model, SameKey, Shape, ShapeType, TShapeId, explore_unique};

/// One solid of the inputs: an input solid, or a solid of an input
/// compound.
struct Lump {
    /// The input it belongs to, by its place in the slice given.
    input: usize,
    /// Its place among that input's solids, where the input is a compound.
    part: Option<usize>,
    shape: Shape,
    /// Its guaranteed bound, widened by the confusion distance.
    bounds: Aabb,
}

impl Lump {
    fn name(&self) -> String {
        match self.part {
            None => format!("input {}", self.input),
            Some(part) => format!("solid {part} of input {}", self.input),
        }
    }
}

/// A shape on the way to the union, with what was done to the inputs to
/// make it.
struct Step {
    shape: Shape,
    history: History,
}

/// The union of `solids`, each a solid or a compound of solids: one solid
/// per group of solids that meet, and a compound of those where groups
/// stay apart.
///
/// - **Disjoint inputs** come back as a compound of separate solids, never
///   as one solid with several outer shells. A solid lying apart from
///   every other is handed back as it is, its nodes shared with the input.
///   A solid inside another's void lies apart from it.
/// - **Solids of one compound** need not stand apart: those that overlap
///   are fused like any other two.
/// - **Contact.** Solids that share volume or a patch of boundary join into
///   one solid. Solids meeting only along an edge or at a point are joined
///   as [`fuse`](crate::fuse) joins the two.
/// - **Faces.** Where faces of two solids that joined meet across an edge
///   and lie on one carrier facing one way, they are merged into one face,
///   and the edges continuing one curve across where they met into one
///   edge: fifty cubes in a row come back as one box of six faces. A
///   solid's own splits, between faces of that solid alone, are kept. Where
///   the merge does not leave a valid solid, the union is given unmerged.
/// - **Order.** The solids are taken in an order set by their bounds, not
///   by their place in `solids`, so a reordering gives the same result,
///   but for solids whose bounds coincide.
/// - **History.** Every face, edge and vertex of the inputs traces to its
///   images in the result, or is deleted where nothing of it is left (a
///   wall two solids shared). Each input solid and each input compound is
///   modified into the result solids holding what became of it; a solid
///   handed back as it is has no record and traces to itself.
///
/// One input that is a single solid is given back as it is.
///
/// # Errors
///
/// [`OgeomError::Construction`] for no input, or an input that is not a
/// solid or a compound of solids (a sheet, a compsolid, an empty compound),
/// naming it by its place in `solids`.
///
/// Where the boolean refuses the solids taken together, they are fused
/// again one at a time, each meeting what was fused before; where that too
/// refuses, the error is the boolean's, its message naming the two inputs
/// that do not fuse, or the input that does not fuse with the union of the
/// others it meets. [`OgeomError::Cancelled`] passes through as it is.
pub fn fuse_all(model: &mut Model, solids: &[Shape], tol: Tolerances) -> OgeomResult<Built> {
    if solids.is_empty() {
        ogeom_bail!(Construction, "a union of many takes at least one solid");
    }
    if let [only] = solids
        && model.kind_of(only)? == ShapeType::Solid
    {
        return Ok(Built::from_nothing(only.clone()));
    }
    let since = model.node_count();
    let lumps = lumps_of(model, solids, tol)?;
    let neighbours = meeting(&lumps);
    let fused = fused_tree(model, &lumps, &neighbours, tol);
    let fused = match fused {
        Err(OgeomError::Cancelled) => return Err(OgeomError::Cancelled),
        Err(_) => fused_one_by_one(model, &lumps, &neighbours, tol)?,
        Ok(fused) => fused,
    };
    let fused = match fused {
        Some(step) => Some(merged(model, &lumps, &neighbours, step, tol)?),
        None => None,
    };

    // The result: the fused solids, and each solid that meets nothing as
    // it came.
    let mut out: Vec<Shape> = Vec::new();
    let mut history = History::new();
    if let Some(step) = fused {
        out.extend(solids_of(model, &step.shape)?);
        history = step.history;
    }
    let passed: Vec<usize> = (0..lumps.len())
        .filter(|&i| neighbours[i].is_empty())
        .collect();
    out.extend(passed.iter().map(|&i| lumps[i].shape.clone()));
    if passed.len() == lumps.len() && solids.len() == 1 {
        return Ok(Built::from_nothing(solids[0].clone()));
    }
    let shape = match out.as_slice() {
        [one] => one.clone(),
        _ => model.add_compound(&out)?,
    };
    model.note_held(&shape, since)?;
    let history = with_solids(model, solids, &lumps, &neighbours, &out, history, tol)?;
    Ok(Built::new(shape, history))
}

/// Every solid of the inputs, with its widened bound, in the order their
/// bounds set.
fn lumps_of(model: &Model, solids: &[Shape], tol: Tolerances) -> OgeomResult<Vec<Lump>> {
    let mut lumps = Vec::new();
    for (input, shape) in solids.iter().enumerate() {
        let parts = match model.kind_of(shape)? {
            ShapeType::Solid => vec![(None, shape.clone())],
            ShapeType::Compound => {
                let mut parts = Vec::new();
                let mut stack = vec![shape.clone()];
                while let Some(at) = stack.pop() {
                    for child in model.children_of(&at)?.into_iter().rev() {
                        match model.kind_of(&child)? {
                            ShapeType::Solid => parts.push(child),
                            ShapeType::Compound => stack.push(child),
                            kind => ogeom_bail!(
                                Construction,
                                "input {input} holds a {kind:?}; a union of many takes \
                                 solids and compounds of solids"
                            ),
                        }
                    }
                }
                if parts.is_empty() {
                    ogeom_bail!(Construction, "input {input} holds no solid");
                }
                parts
                    .into_iter()
                    .enumerate()
                    .map(|(part, solid)| (Some(part), solid))
                    .collect()
            }
            kind => ogeom_bail!(
                Construction,
                "input {input} is a {kind:?}; a union of many takes solids and \
                 compounds of solids"
            ),
        };
        for (part, shape) in parts {
            let bounds = ogeom_algo::shape_bounds(model, &shape, tol)?.expanded(tol.confusion());
            if bounds.is_empty() {
                ogeom_bail!(Construction, "input {input} holds a solid with no boundary");
            }
            lumps.push(Lump {
                input,
                part,
                shape,
                bounds,
            });
        }
    }
    let key = |lump: &Lump| {
        let corner = |p: Option<ogeom_math::Point>| p.map_or([0.0; 3], |p| [p.x, p.y, p.z]);
        let (low, high) = (corner(lump.bounds.low()), corner(lump.bounds.high()));
        [low[0], low[1], low[2], high[0], high[1], high[2]]
    };
    lumps.sort_by(|a, b| {
        key(a)
            .iter()
            .zip(key(b).iter())
            .map(|(x, y)| x.total_cmp(y))
            .find(|o| o.is_ne())
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(a.input.cmp(&b.input))
            .then(a.part.cmp(&b.part))
    });
    Ok(lumps)
}

/// For each solid, the others whose widened bounds overlap its own: the
/// only ones it can meet, since each solid lies within its bound. Sorted.
fn meeting(lumps: &[Lump]) -> Vec<Vec<usize>> {
    let low_x = |i: usize| lumps[i].bounds.low().map_or(0.0, |p| p.x);
    let high_x = |i: usize| lumps[i].bounds.high().map_or(0.0, |p| p.x);
    let mut by_low: Vec<usize> = (0..lumps.len()).collect();
    by_low.sort_by(|&a, &b| low_x(a).total_cmp(&low_x(b)));
    let mut neighbours = vec![Vec::new(); lumps.len()];
    let mut open: Vec<usize> = Vec::new();
    for &i in &by_low {
        open.retain(|&j| high_x(j) >= low_x(i));
        for &j in &open {
            if lumps[i].bounds.intersects(&lumps[j].bounds) {
                neighbours[i].push(j);
                neighbours[j].push(i);
            }
        }
        open.push(i);
    }
    for list in &mut neighbours {
        list.sort_unstable();
    }
    neighbours
}

/// The solids that meet anything, coloured so that no two of one colour
/// meet, each colour a compound, and the colours fused in a balanced tree.
/// `None` where no solid meets another.
fn fused_tree(
    model: &mut Model,
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    tol: Tolerances,
) -> OgeomResult<Option<Step>> {
    let mut colour: Vec<Option<usize>> = vec![None; lumps.len()];
    let mut colours: Vec<Vec<usize>> = Vec::new();
    for i in 0..lumps.len() {
        if neighbours[i].is_empty() {
            continue;
        }
        let taken: FastSet<usize> = neighbours[i].iter().filter_map(|&j| colour[j]).collect();
        let c = (0..).find(|c| !taken.contains(c)).unwrap_or(0);
        colour[i] = Some(c);
        if c == colours.len() {
            colours.push(Vec::new());
        }
        colours[c].push(i);
    }
    let mut level: Vec<Step> = Vec::with_capacity(colours.len());
    for members in &colours {
        let shapes: Vec<Shape> = members.iter().map(|&i| lumps[i].shape.clone()).collect();
        let shape = match shapes.as_slice() {
            [one] => one.clone(),
            _ => model.add_compound(&shapes)?,
        };
        level.push(Step {
            shape,
            history: History::new(),
        });
    }
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        let mut steps = level.into_iter();
        while let Some(a) = steps.next() {
            let Some(b) = steps.next() else {
                next.push(a);
                break;
            };
            let fused = crate::fuse(model, &a.shape, &b.shape, tol)?;
            next.push(Step {
                shape: fused.shape,
                history: a.history.joined(b.history).followed_by(&fused.history),
            });
        }
        level = next;
    }
    Ok(level.pop())
}

/// The solids that meet anything fused one at a time, each group grown
/// from its first solid by the solids meeting what it holds; the refusal
/// named by the solids it came from.
fn fused_one_by_one(
    model: &mut Model,
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    tol: Tolerances,
) -> OgeomResult<Option<Step>> {
    let mut placed = vec![false; lumps.len()];
    let mut groups: Vec<Step> = Vec::new();
    for start in 0..lumps.len() {
        if placed[start] || neighbours[start].is_empty() {
            continue;
        }
        placed[start] = true;
        let mut held = vec![start];
        let mut queue = std::collections::VecDeque::from([start]);
        let mut step = Step {
            shape: lumps[start].shape.clone(),
            history: History::new(),
        };
        while let Some(at) = queue.pop_front() {
            for &next in &neighbours[at] {
                if placed[next] {
                    continue;
                }
                placed[next] = true;
                queue.push_back(next);
                match crate::fuse(model, &step.shape, &lumps[next].shape, tol) {
                    Ok(fused) => {
                        step = Step {
                            shape: fused.shape,
                            history: step.history.followed_by(&fused.history),
                        };
                    }
                    Err(e) => return Err(refusal(model, lumps, neighbours, &held, next, e, tol)),
                }
                held.push(next);
            }
        }
        groups.push(step);
    }
    let Some(mut step) = groups.pop() else {
        return Ok(None);
    };
    if !groups.is_empty() {
        let mut shapes = Vec::new();
        for group in groups.iter().chain(core::iter::once(&step)) {
            shapes.extend(solids_of(model, &group.shape)?);
        }
        let history = groups
            .into_iter()
            .fold(History::new(), |all, group| all.joined(group.history))
            .joined(step.history);
        step = Step {
            shape: model.add_compound(&shapes)?,
            history,
        };
    }
    Ok(Some(step))
}

/// The boolean's refusal of `lump` against the union of `held`, named: by
/// the one solid of `held` it does not fuse with, where fusing it with each
/// it meets alone finds one, and by the union otherwise.
fn refusal(
    model: &mut Model,
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    held: &[usize],
    lump: usize,
    e: OgeomError,
    tol: Tolerances,
) -> OgeomError {
    if matches!(e, OgeomError::Cancelled) {
        return e;
    }
    let name = lumps[lump].name();
    for &other in neighbours[lump].iter().filter(|j| held.contains(j)) {
        match crate::fuse(model, &lumps[other].shape, &lumps[lump].shape, tol) {
            Ok(_) => {}
            Err(OgeomError::Cancelled) => return OgeomError::Cancelled,
            Err(pair) => {
                let other = lumps[other].name();
                return renamed(pair, format!("{other} and {name} do not fuse"));
            }
        }
    }
    let met: Vec<String> = neighbours[lump]
        .iter()
        .filter(|j| held.contains(j))
        .map(|&j| lumps[j].name())
        .collect();
    renamed(
        e,
        format!(
            "{name} does not fuse with the union of {}, though it fuses with each alone",
            met.join(", ")
        ),
    )
}

/// `e` with `what` before its message.
fn renamed(e: OgeomError, what: String) -> OgeomError {
    let told = |cause: Cause| Cause::Owned(format!("{what}: {cause}"));
    match e {
        OgeomError::Construction(c) => OgeomError::Construction(told(c)),
        OgeomError::Domain(c) => OgeomError::Domain(told(c)),
        OgeomError::Range(c) => OgeomError::Range(told(c)),
        OgeomError::Dimension(c) => OgeomError::Dimension(told(c)),
        OgeomError::NullObject(c) => OgeomError::NullObject(told(c)),
        OgeomError::Dangling(c) => OgeomError::Dangling(told(c)),
        OgeomError::Numeric(c) => OgeomError::Numeric(told(c)),
        OgeomError::NotDone(c) => OgeomError::NotDone(told(c)),
        OgeomError::Invariant(c) => OgeomError::Invariant(told(c)),
        OgeomError::Unimplemented(c) => OgeomError::Unimplemented(told(c)),
        other => other,
    }
}

/// The solids of a boolean's result: itself, or the solids of its
/// compound.
fn solids_of(model: &Model, shape: &Shape) -> OgeomResult<Vec<Shape>> {
    if model.kind_of(shape)? == ShapeType::Solid {
        return Ok(vec![shape.clone()]);
    }
    explore_unique(model, shape, ShapeType::Solid)
}

/// The faces of the solids `step` fused, with which solid each came from.
fn owners(
    model: &Model,
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    history: &History,
) -> OgeomResult<FastMap<SameKey, usize>> {
    let mut owner: FastMap<SameKey, usize> = FastMap::default();
    for (i, lump) in lumps.iter().enumerate() {
        if neighbours[i].is_empty() {
            continue;
        }
        for face in explore_unique(model, &lump.shape, ShapeType::Face)? {
            for image in history.trace(&face) {
                owner.entry(SameKey(image.clone())).or_insert(i);
            }
        }
    }
    Ok(owner)
}

/// `step` with the faces of different solids lying on one carrier across
/// an edge merged, and the edges continuing one curve across where they
/// met joined; `step` as it is where nothing merges or the merge leaves no
/// valid solid.
fn merged(
    model: &mut Model,
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    step: Step,
    tol: Tolerances,
) -> OgeomResult<Step> {
    let owner = owners(model, lumps, neighbours, &step.history)?;
    let faces = explore_unique(model, &step.shape, ShapeType::Face)?;
    let mut users: FastMap<SameKey, Vec<usize>> = FastMap::default();
    for (index, face) in faces.iter().enumerate() {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            users.entry(SameKey(edge)).or_default().push(index);
        }
    }
    let mut seam = vec![false; faces.len()];
    for list in users.values() {
        let mut from = list
            .iter()
            .filter_map(|&f| owner.get(&SameKey(faces[f].clone())));
        if let Some(first) = from.next()
            && from.any(|o| o != first)
        {
            for &f in list {
                seam[f] = true;
            }
        }
    }
    let around: Vec<Shape> = faces
        .iter()
        .zip(&seam)
        .filter(|(_, s)| **s)
        .map(|(f, _)| f.clone())
        .collect();
    if around.is_empty() {
        return Ok(step);
    }
    let unified = match ogeom_algo::unify_same_domain_around(model, &step.shape, &around, tol) {
        Ok((_, 0)) => return Ok(step),
        Ok((unified, _)) => unified,
        Err(OgeomError::Cancelled) => return Err(OgeomError::Cancelled),
        Err(_) => return Ok(step),
    };
    // The vertices the merge may dissolve: the ends of the edges it took
    // out between the faces it merged.
    let kept: FastSet<SameKey> = explore_unique(model, &unified.shape, ShapeType::Edge)?
        .into_iter()
        .map(SameKey)
        .collect();
    let mut at: FastSet<TShapeId> = FastSet::default();
    for face in &around {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            let images = unified.history.trace(&edge);
            if images.iter().any(|e| kept.contains(&SameKey(e.clone()))) {
                continue;
            }
            for vertex in explore_unique(model, &edge, ShapeType::Vertex)? {
                at.extend(unified.history.trace(&vertex).iter().map(Shape::node));
            }
        }
    }
    let mut built = unified;
    if !at.is_empty() {
        match ogeom_algo::merge_edges_at(model, &built.shape, Some(&at), tol) {
            Ok((_, 0)) => {}
            Ok((joined, _)) => {
                built = Built::new(joined.shape, built.history.followed_by(&joined.history));
            }
            Err(OgeomError::Cancelled) => return Err(OgeomError::Cancelled),
            Err(_) => {}
        }
    }
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(step);
    }
    Ok(Step {
        shape: built.shape,
        history: step.history.followed_by(&built.history),
    })
}

/// `history` told only of the inputs' faces, edges and vertices, with each
/// input solid that met another, and each input compound, modified into
/// the solids of `out` holding what became of it.
fn with_solids(
    model: &Model,
    solids: &[Shape],
    lumps: &[Lump],
    neighbours: &[Vec<usize>],
    out: &[Shape],
    history: History,
    tol: Tolerances,
) -> OgeomResult<History> {
    let mut originals: Vec<Shape> = Vec::new();
    let mut held_parts: FastSet<SameKey> = FastSet::default();
    for kind in [ShapeType::Face, ShapeType::Edge, ShapeType::Vertex] {
        for (i, lump) in lumps.iter().enumerate() {
            if !neighbours[i].is_empty() {
                originals.extend(explore_unique(model, &lump.shape, kind)?);
            }
        }
        for solid in out {
            held_parts.extend(explore_unique(model, solid, kind)?.into_iter().map(SameKey));
        }
    }
    let keys: FastSet<SameKey> = originals.iter().cloned().map(SameKey).collect();
    let mut history = history.only_for(|shape| keys.contains(&SameKey(shape.clone())));
    // A step that merged or joined may leave a shape it took out unrecorded:
    // what the result does not hold is not an image.
    for part in &originals {
        let images = history.trace(part);
        let held: Vec<Shape> = images
            .iter()
            .filter(|image| held_parts.contains(&SameKey((*image).clone())))
            .cloned()
            .collect();
        if held.len() == images.len() {
            continue;
        }
        history.delete(part);
        for image in held {
            history.modify(part, image);
        }
    }

    let mut holder: FastMap<SameKey, usize> = FastMap::default();
    for (index, solid) in out.iter().enumerate() {
        for face in explore_unique(model, solid, ShapeType::Face)? {
            holder.insert(SameKey(face), index);
        }
    }
    let place: FastMap<SameKey, usize> = out
        .iter()
        .enumerate()
        .map(|(index, solid)| (SameKey(solid.clone()), index))
        .collect();
    let mut bounds: Vec<Option<Aabb>> = vec![None; out.len()];
    let mut images: Vec<Vec<usize>> = vec![Vec::new(); lumps.len()];
    for (i, lump) in lumps.iter().enumerate() {
        if neighbours[i].is_empty() {
            images[i] = place
                .get(&SameKey(lump.shape.clone()))
                .copied()
                .into_iter()
                .collect();
            continue;
        }
        let mut held: Vec<usize> = Vec::new();
        for face in explore_unique(model, &lump.shape, ShapeType::Face)? {
            for image in history.trace(&face) {
                if let Some(&at) = holder.get(&SameKey(image.clone()))
                    && !held.contains(&at)
                {
                    held.push(at);
                }
            }
        }
        // A solid wholly inside another left no face: its material is in
        // the solid holding its bound.
        if held.is_empty() {
            for (index, solid) in out.iter().enumerate() {
                if bounds[index].is_none() {
                    bounds[index] = Some(
                        ogeom_algo::shape_bounds(model, solid, tol)?.expanded(tol.confusion()),
                    );
                }
                if bounds[index]
                    .as_ref()
                    .is_some_and(|b| b.contains_box(&lump.bounds))
                {
                    held.push(index);
                }
            }
        }
        held.sort_unstable();
        for &at in &held {
            history.modify(&lump.shape, out[at].clone());
        }
        images[i] = held;
    }
    for (input, shape) in solids.iter().enumerate() {
        if model.kind_of(shape)? == ShapeType::Solid {
            continue;
        }
        let mut held: Vec<usize> = lumps
            .iter()
            .enumerate()
            .filter(|(_, lump)| lump.input == input)
            .flat_map(|(i, _)| images[i].iter().copied())
            .collect();
        held.sort_unstable();
        held.dedup();
        for at in held {
            history.modify(shape, out[at].clone());
        }
    }
    Ok(history)
}
