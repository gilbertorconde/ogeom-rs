//! Free bounds: the edges of a shape only one face walks, chained into the
//! contours they form.
//!
//! A shell open along some of its edges is open along contours: the rim of
//! a sheet, the edge of a hole where a face is missing. Each free edge is
//! walked by its one face in the direction that face's material keeps to
//! the left of, so where the faces agree in orientation, the free edges
//! chain head to tail round each opening. A contour that comes back to the
//! vertex it left from is closed: an opening a face could fill. One that
//! stops where no free edge carries on (two openings meeting at vertices
//! the model holds apart, or faces that disagree in orientation) is open.
//!
//! Nothing is joined here. Free edges within the tolerances they state are
//! already one edge to sewing, and [`fix_shape`](crate::fix_shape()) says
//! which pairs farther apart are each other's twins.

use ogeom_algo::edge_vertices;
use ogeom_core::{FastMap, OgeomResult};
use ogeom_topo::{Model, Shape, ShapeType, TShapeId, explore_unique};

/// The free edges of a shape, chained.
#[derive(Debug, Clone, Default)]
pub struct FreeBounds {
    /// Contours that come back to where they began: each the edges in the
    /// order they walk, as their faces walk them.
    pub closed: Vec<Vec<Shape>>,
    /// Chains that stop at a vertex where no free edge carries on.
    pub open: Vec<Vec<Shape>>,
}

impl FreeBounds {
    /// Whether the shape has no free edge at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.closed.is_empty() && self.open.is_empty()
    }

    /// Each contour as a wire of its edges in the order they walk, the
    /// closed contours first, then the open chains.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if an
    /// edge is not in `model`.
    pub fn wires(&self, model: &mut Model) -> OgeomResult<Vec<Shape>> {
        self.closed
            .iter()
            .chain(&self.open)
            .map(|edges| model.add_wire(edges))
            .collect()
    }
}

/// One free edge as its face walks it.
struct Free {
    edge: Shape,
    from: TShapeId,
    to: TShapeId,
    /// The free edge its wire walks next, where that is one.
    next_in_wire: Option<usize>,
}

/// The free bounds of `shape`: every edge only one of its faces walks, a
/// degenerate edge aside, chained head to tail into contours. See the
/// module documentation.
///
/// At a vertex where more than one free edge carries on, the one the same
/// wire walks next is taken first, so a face's own rim comes out whole.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// fails to resolve.
pub fn free_bounds(model: &Model, shape: &Shape) -> OgeomResult<FreeBounds> {
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let mut walks: Vec<Vec<Shape>> = Vec::new();
    let mut uses: FastMap<TShapeId, usize> = FastMap::default();
    for face in &faces {
        for wire in model.ordered_children_of(face)? {
            let edges = model.ordered_children_of(&wire)?;
            for edge in &edges {
                *uses.entry(edge.node()).or_default() += 1;
            }
            walks.push(edges);
        }
    }
    let is_free = |edge: &Shape| {
        uses.get(&edge.node()) == Some(&1)
            && !model
                .node(edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
    };
    let mut free: Vec<Free> = Vec::new();
    for edges in &walks {
        let first = free.len();
        let mut at: Vec<Option<usize>> = Vec::with_capacity(edges.len());
        for edge in edges {
            if !is_free(edge) {
                at.push(None);
                continue;
            }
            let Some((start, end)) = edge_vertices(model, edge)? else {
                at.push(None);
                continue;
            };
            let (from, to) = (start.node(), end.node());
            at.push(Some(free.len()));
            free.push(Free {
                edge: edge.clone(),
                from,
                to,
                next_in_wire: None,
            });
        }
        if free.len() == first {
            continue;
        }
        for i in 0..at.len() {
            if let (Some(k), Some(next)) = (at[i], at[(i + 1) % at.len()]) {
                free[k].next_in_wire = Some(next);
            }
        }
    }

    let mut leaving: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    let mut arriving: FastMap<TShapeId, usize> = FastMap::default();
    for (k, f) in free.iter().enumerate() {
        leaving.entry(f.from).or_default().push(k);
        *arriving.entry(f.to).or_default() += 1;
    }
    let mut used = vec![false; free.len()];
    let mut out = FreeBounds::default();
    // Chains from vertices no free edge arrives at first, so an open chain
    // is found whole from its head; what is left goes round.
    let heads: Vec<usize> = (0..free.len())
        .filter(|&k| !arriving.contains_key(&free[k].from))
        .collect();
    let rest = 0..free.len();
    for start in heads.into_iter().chain(rest) {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = vec![start];
        let mut last = start;
        while free[last].to != free[start].from {
            let next = free[last]
                .next_in_wire
                .filter(|&n| !used[n] && free[n].from == free[last].to)
                .or_else(|| {
                    leaving
                        .get(&free[last].to)
                        .and_then(|ks| ks.iter().copied().find(|&n| !used[n]))
                });
            let Some(next) = next else { break };
            used[next] = true;
            chain.push(next);
            last = next;
        }
        let edges: Vec<Shape> = chain.iter().map(|&k| free[k].edge.clone()).collect();
        if free[last].to == free[start].from {
            out.closed.push(edges);
        } else {
            out.open.push(edges);
        }
    }
    Ok(out)
}
