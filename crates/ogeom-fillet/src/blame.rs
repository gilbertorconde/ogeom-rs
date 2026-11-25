//! Refusals of a blend request that name the request's edges they
//! concern, and builds of what of a request can be built.

use crate::support::Blamed;
use crate::{Chamfer, Overflow};
use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, TShapeId};

/// A blend request's refusal, naming the edges of the request it concerns.
#[derive(Debug, Clone)]
pub struct BlendRefusal {
    /// Indices into the request's edge list: the edge refused, then any
    /// edge it collides with. Empty where the refusal is no edge's (a
    /// request of no edges, a cancelled one).
    pub edges: Vec<usize>,
    /// Why.
    pub error: OgeomError,
}

impl core::fmt::Display for BlendRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "edges {:?} of the request: {}", self.edges, self.error)
    }
}

impl std::error::Error for BlendRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

/// [`fillet_edges`](crate::fillet_edges), its refusal naming the edges of
/// `edges` it concerns, by their indices there.
///
/// On success the same as [`fillet_edges`](crate::fillet_edges). Refused,
/// the indices are read from the request's check where it refused there:
/// the edge whose ball does not sit on its faces, or the edge whose
/// contact stands in another's band, then that other. A refusal of the
/// build after the check names the first edge of the request whose
/// prefix of the request does not build (each prefix built on a copy of
/// the model), then, where that edge builds alone, the edges before it
/// sharing a vertex or a face with it.
///
/// # Errors
///
/// As [`fillet_edges`](crate::fillet_edges), in a [`BlendRefusal`].
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: `BlendRequest::new().round(edges, r).build_checked(..)`, which names every refused edge"
)]
pub fn fillet_edges_checked(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    // What a failed request leaves in the model is beside the given
    // solid, whose nodes the request copies before it works on them.
    let error = match crate::fillet::rounds(model, solid, edges, radius, Overflow::default(), tol) {
        Ok(built) => return Ok(built),
        Err(error) => error,
    };
    if matches!(error, OgeomError::Cancelled) || edges.is_empty() {
        return Err(BlendRefusal {
            edges: Vec::new(),
            error,
        });
    }
    if let Err(Blamed { edges: named, .. }) =
        crate::engine::request_fits(model, solid, edges, &vec![radius; edges.len()], tol)
        && !named.is_empty()
    {
        return Err(BlendRefusal {
            edges: named,
            error,
        });
    }
    let build = |model: &mut Model, set: &[usize]| {
        let picked: Vec<Shape> = set.iter().map(|&i| edges[i].clone()).collect();
        crate::fillet::rounds(model, solid, &picked, radius, Overflow::default(), tol)
    };
    Err(first_failing_prefix(model, solid, edges, error, build))
}

/// [`fillet_edges_corner`](crate::fillet_edges_corner), its refusal naming
/// the edges of `edges` it concerns, by their indices there, read as
/// [`fillet_edges_checked`] reads them, each prefix of the request built
/// with the same corners.
///
/// # Errors
///
/// As [`fillet_edges_corner`](crate::fillet_edges_corner), in a
/// [`BlendRefusal`].
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: `BlendRequest::new().round(edges, r).round_corner(c).build_checked(..)`, which names every refused edge"
)]
#[allow(
    deprecated,
    reason = "the setback reading is fillet_edges_checked's own"
)]
pub fn fillet_edges_corner_checked(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: crate::FilletCorner,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    if corner == crate::FilletCorner::Setback {
        return fillet_edges_checked(model, solid, edges, radius, tol);
    }
    let error = match crate::fillet::rounds_corner(
        model,
        solid,
        edges,
        radius,
        corner,
        Overflow::default(),
        tol,
    ) {
        Ok(built) => return Ok(built),
        Err(error) => error,
    };
    if matches!(error, OgeomError::Cancelled) || edges.is_empty() {
        return Err(BlendRefusal {
            edges: Vec::new(),
            error,
        });
    }
    if let Err(Blamed { edges: named, .. }) =
        crate::engine::request_fits(model, solid, edges, &vec![radius; edges.len()], tol)
        && !named.is_empty()
    {
        return Err(BlendRefusal {
            edges: named,
            error,
        });
    }
    let build = |model: &mut Model, set: &[usize]| {
        let picked: Vec<Shape> = set.iter().map(|&i| edges[i].clone()).collect();
        crate::fillet::rounds_corner(
            model,
            solid,
            &picked,
            radius,
            corner,
            Overflow::default(),
            tol,
        )
    };
    Err(first_failing_prefix(model, solid, edges, error, build))
}

/// [`chamfer_edges_with`](crate::chamfer_edges_with), its refusal naming
/// the specs it concerns, by their indices in `specs`: a wedge refused
/// names its spec, and a wedge that does not apply names its spec, then
/// the specs before it whose edges share a vertex with its edge.
///
/// A vertex of three bevels or more takes the facet
/// ([`ChamferCorner::Facet`](crate::ChamferCorner::Facet)), as
/// [`chamfer_edges_with`](crate::chamfer_edges_with) does; the mitre is
/// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered) with
/// [`ChamferCorner::Point`](crate::ChamferCorner::Point).
///
/// # Errors
///
/// As [`chamfer_edges_with`](crate::chamfer_edges_with), in a
/// [`BlendRefusal`].
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: one `.edge(e, spec)` per spec, then `.build_checked(..)`, which names every refused edge"
)]
pub fn chamfer_edges_checked(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    crate::chamfer::chamfer_edges_blamed(model, solid, specs, tol).map_err(|blamed| BlendRefusal {
        edges: blamed.edges,
        error: blamed.error,
    })
}

/// [`fillet_edges_corner_overflow`](crate::fillet_edges_corner_overflow),
/// its refusal naming every edge of `edges` it refuses, by their indices
/// there.
///
/// On success exactly the same [`Built`] as
/// [`fillet_edges_corner_overflow`](crate::fillet_edges_corner_overflow).
/// Refused, the request is read in order, each set built on a copy of the
/// model with the same `corner` and `overflow`: an edge joins the
/// accepted edges where they build with it, and is refused where they do
/// not. [`BlendRefusal::edges`] lists the refused edges in request order,
/// so that the request without them builds with the same options, and
/// [`BlendRefusal::error`] is the request's own refusal. Under overflow
/// options no check of the default reading refuses an edge unbuilt: a
/// request it refuses may build with another behaviour, so every reading
/// is a build (ADR 0001, amendment "checked blends under overflow").
///
/// Where two edges are refused only together, the later in the request is
/// named: reorder the request to keep another.
///
/// # Errors
///
/// As [`fillet_edges_corner_overflow`](crate::fillet_edges_corner_overflow),
/// in a [`BlendRefusal`]; its `edges` empty where cancelled or where
/// `edges` is empty.
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: `BlendRequest::new().round(edges, r).round_corner(c).overflow(o).build_checked(..)`"
)]
pub fn fillet_edges_corner_overflow_checked(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: crate::FilletCorner,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    crate::BlendRequest::new()
        .round(edges.iter().cloned(), radius)
        .round_corner(corner)
        .overflow(overflow)
        .build_checked(model, solid, tol)
}

/// [`fillet_edges_corner_overflow_checked`]'s reading, which a
/// [`BlendRequest`](crate::BlendRequest) of rounds of one radius reads
/// checked.
pub(crate) fn rounds_checked(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corner: crate::FilletCorner,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    let error =
        match crate::fillet::rounds_corner(model, solid, edges, radius, corner, overflow, tol) {
            Ok(built) => return Ok(built),
            Err(error) => error,
        };
    let build = |model: &mut Model, set: &[usize]| {
        let picked: Vec<Shape> = set.iter().map(|&i| edges[i].clone()).collect();
        crate::fillet::rounds_corner(model, solid, &picked, radius, corner, overflow, tol)
    };
    Err(refused_in_order(model, solid, edges, error, build))
}

/// [`chamfer_edges_overflow`](crate::chamfer_edges_overflow), its refusal
/// naming every spec it refuses, by their indices in `specs`, read as
/// [`fillet_edges_corner_overflow_checked`] reads a round's: each set built
/// with the same `corner` and `overflow`, the refused specs listed in
/// request order so that the request without them builds. Covers both
/// corners, [`ChamferCorner::Point`](crate::ChamferCorner::Point) as well
/// as [`ChamferCorner::Facet`](crate::ChamferCorner::Facet).
///
/// On success exactly the same [`Built`] as
/// [`chamfer_edges_overflow`](crate::chamfer_edges_overflow).
///
/// # Errors
///
/// As [`chamfer_edges_overflow`](crate::chamfer_edges_overflow), in a
/// [`BlendRefusal`]; its `edges` empty where cancelled or where `specs` is
/// empty.
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: one `.edge(e, spec)` per spec, then `.chamfer_corner(c).overflow(o).build_checked(..)`"
)]
pub fn chamfer_edges_overflow_checked(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: crate::ChamferCorner,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    crate::BlendRequest::bevels(specs)
        .chamfer_corner(corner)
        .overflow(overflow)
        .build_checked(model, solid, tol)
}

/// [`chamfer_edges_overflow_checked`]'s reading, which a
/// [`BlendRequest`](crate::BlendRequest) of bevels alone reads checked.
pub(crate) fn bevels_checked(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: crate::ChamferCorner,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> Result<Built, BlendRefusal> {
    let bevelled = |model: &mut Model, specs: &[(Shape, Chamfer)]| {
        crate::chamfer_overflow::bevels_overflow(model, solid, specs, corner, overflow, tol)
    };
    let error = match bevelled(model, specs) {
        Ok(built) => return Ok(built),
        Err(error) => error,
    };
    let edges: Vec<Shape> = specs.iter().map(|(edge, _)| edge.clone()).collect();
    let build = |model: &mut Model, set: &[usize]| {
        let picked: Vec<(Shape, Chamfer)> = set.iter().map(|&i| specs[i].clone()).collect();
        bevelled(model, &picked)
    };
    Err(refused_in_order(model, solid, &edges, error, build))
}

/// Round what of a request can be rounded, leaving the rest sharp.
///
/// The edges are taken in request order; an edge joins the built set when
/// [`fillet_edges`](crate::fillet_edges) builds it together with the edges
/// accepted so far, corners, junctions and meetings between them as a
/// full call builds them; an edge that does not is left sharp and
/// reported, naming it, then the accepted edges the request's check finds
/// it collides with, or, refused by the build, the accepted edges sharing
/// a vertex or a face with it. Each attempt is built on a copy of the
/// model, and the model then becomes the copy holding the last attempt
/// that built: the [`Built`] is [`fillet_edges`](crate::fillet_edges) on
/// exactly the accepted edges, and its history says nothing of the
/// refused ones. An edge the request's check refuses together with the
/// accepted ones is refused without building.
///
/// # Errors
///
/// As [`fillet_edges`](crate::fillet_edges) where no edge is accepted (the
/// first edge's refusal), and [`OgeomError::Cancelled`] when cancelled.
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: `BlendRequest::new().round(edges, r).build_partial(..)`"
)]
pub fn fillet_edges_partial(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
    crate::BlendRequest::new()
        .round(edges.iter().cloned(), radius)
        .build_partial(model, solid, tol)
}

/// [`fillet_edges_partial`]'s reading, which a
/// [`BlendRequest`](crate::BlendRequest) of rounds of one radius reads in
/// part.
pub(crate) fn rounds_partial(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
    let pick = |set: &[usize]| -> Vec<Shape> { set.iter().map(|&i| edges[i].clone()).collect() };
    partial(
        model,
        solid,
        edges,
        |model, set| {
            crate::fillet::rounds(model, solid, &pick(set), radius, Overflow::default(), tol)
        },
        |model, set| {
            if set.len() < 2 {
                return None;
            }
            match crate::engine::request_fits(
                model,
                solid,
                &pick(set),
                &vec![radius; set.len()],
                tol,
            ) {
                Err(Blamed { edges, error }) if !edges.is_empty() => {
                    Some((edges.iter().map(|&i| set[i]).collect(), error))
                }
                _ => None,
            }
        },
    )
}

/// Bevel what of a request can be bevelled, leaving the rest sharp, as
/// [`fillet_edges_partial`] rounds: the specs taken in request order, each
/// joining the built set when
/// [`chamfer_edges_with`](crate::chamfer_edges_with) builds it with the
/// specs accepted so far, a refused spec reported with the accepted specs
/// whose edges share a vertex or a face with its edge.
///
/// # Errors
///
/// As [`chamfer_edges_with`](crate::chamfer_edges_with) where no spec is
/// accepted, and [`OgeomError::Cancelled`] when cancelled.
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: one `.edge(e, spec)` per spec, then `.build_partial(..)`"
)]
pub fn chamfer_edges_partial(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
    crate::BlendRequest::bevels(specs).build_partial(model, solid, tol)
}

/// [`chamfer_edges_partial`]'s reading, which a
/// [`BlendRequest`](crate::BlendRequest) of bevels alone reads in part.
pub(crate) fn bevels_partial(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
    let edges: Vec<Shape> = specs.iter().map(|(edge, _)| edge.clone()).collect();
    partial(
        model,
        solid,
        &edges,
        |model, set| {
            let picked: Vec<(Shape, Chamfer)> = set.iter().map(|&i| specs[i].clone()).collect();
            crate::chamfer::bevels(model, solid, &picked, tol)
        },
        |_, _| None,
    )
}

/// The vertices and faces of `solid` each of `edges` touches.
fn touches(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
) -> OgeomResult<Vec<(Vec<TShapeId>, Vec<TShapeId>)>> {
    edges
        .iter()
        .map(|edge| {
            let ends = ogeom_algo::edge_vertices(model, edge)?
                .map(|(a, b)| vec![a.node(), b.node()])
                .unwrap_or_default();
            let faces = ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?
                .iter()
                .map(Shape::node)
                .collect();
            Ok((ends, faces))
        })
        .collect()
}

/// Whether two edges, as [`touches`] reads them, share a vertex or a face.
fn meet(a: &(Vec<TShapeId>, Vec<TShapeId>), b: &(Vec<TShapeId>, Vec<TShapeId>)) -> bool {
    a.0.iter().any(|v| b.0.contains(v)) || a.1.iter().any(|f| b.1.contains(f))
}

/// The refusal of a request whose check passed and whose build did not:
/// the first edge whose prefix of the request does not build, each on a
/// copy of the model, then, where that edge builds alone, the edges
/// before it it meets. The whole request, which did not build, is the last
/// prefix.
fn first_failing_prefix(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    error: OgeomError,
    build: impl Fn(&mut Model, &[usize]) -> OgeomResult<Built>,
) -> BlendRefusal {
    let refusal = |edges: Vec<usize>, error: OgeomError| BlendRefusal { edges, error };
    let n = edges.len();
    let mut culprit = (n - 1, error);
    for k in 0..n - 1 {
        let prefix: Vec<usize> = (0..=k).collect();
        match build(&mut model.clone(), &prefix) {
            Ok(_) => {}
            Err(e @ OgeomError::Cancelled) => return refusal(Vec::new(), e),
            Err(e) => {
                culprit = (k, e);
                break;
            }
        }
    }
    let (k, error) = culprit;
    let mut named = vec![k];
    match build(&mut model.clone(), &[k]) {
        Err(e @ OgeomError::Cancelled) => return refusal(Vec::new(), e),
        Err(_) => {}
        Ok(_) => {
            if let Ok(touched) = touches(model, solid, edges) {
                named.extend((0..k).filter(|&j| meet(&touched[k], &touched[j])));
            }
        }
    }
    refusal(named, error)
}

/// What [`fillet_edges_partial`] and [`chamfer_edges_partial`] share:
/// `build` on a set of the request's indices, and `check`, the request's
/// check on a set, refusing it before any build where it names edges. The
/// whole request is tried first, and where it builds nothing is refused.
fn partial(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    build: impl Fn(&mut Model, &[usize]) -> OgeomResult<Built>,
    check: impl Fn(&Model, &[usize]) -> Option<(Vec<usize>, OgeomError)>,
) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
    let all: Vec<usize> = (0..edges.len()).collect();
    if check(model, &all).is_none() {
        let mut trial = model.clone();
        match build(&mut trial, &all) {
            Ok(built) => {
                *model = trial;
                return Ok((built, Vec::new()));
            }
            Err(e @ OgeomError::Cancelled) => return Err(e),
            Err(_) => {}
        }
    }
    let (last, refused) = in_order(model, solid, edges, build, check)?;
    let Some((built_model, built)) = last else {
        return Err(match refused.into_iter().next() {
            Some(first) => first.error,
            None => ogeom_core::ogeom_err!(Construction, "a request of no edges blends nothing"),
        });
    };
    *model = built_model;
    Ok((built, refused))
}

/// The request's edges taken in order, each joining the accepted set
/// where `check` passes the set with it and `build` builds that set on a
/// copy of `model`; each refused edge reported, naming it, then the edges
/// `check` names with it or, refused by the build, the accepted edges it
/// meets. With the model holding the last set that built and its result,
/// or `None` where no edge was accepted.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] when cancelled, and as reading the edges'
/// vertices and faces.
#[allow(
    clippy::type_complexity,
    reason = "the accepted build, then the refusals; used by two callers only"
)]
fn in_order(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    build: impl Fn(&mut Model, &[usize]) -> OgeomResult<Built>,
    check: impl Fn(&Model, &[usize]) -> Option<(Vec<usize>, OgeomError)>,
) -> OgeomResult<(Option<(Model, Built)>, Vec<BlendRefusal>)> {
    let touched = touches(model, solid, edges)?;
    let mut accepted: Vec<usize> = Vec::new();
    let mut refused: Vec<BlendRefusal> = Vec::new();
    let mut last: Option<(Model, Built)> = None;
    for k in 0..edges.len() {
        let mut set = accepted.clone();
        set.push(k);
        if let Some((named, error)) = check(model, &set) {
            // The edge refused first, then those it collides with.
            let mut named_k = vec![k];
            named_k.extend(named.into_iter().filter(|&i| i != k));
            refused.push(BlendRefusal {
                edges: named_k,
                error,
            });
            continue;
        }
        let mut trial = model.clone();
        match build(&mut trial, &set) {
            Ok(built) => {
                accepted.push(k);
                last = Some((trial, built));
            }
            Err(e @ OgeomError::Cancelled) => return Err(e),
            Err(error) => {
                let mut named = vec![k];
                named.extend(
                    accepted
                        .iter()
                        .copied()
                        .filter(|&j| meet(&touched[k], &touched[j])),
                );
                refused.push(BlendRefusal {
                    edges: named,
                    error,
                });
            }
        }
    }
    Ok((last, refused))
}

/// The refusal of a request under its options that did not build with
/// `error`: every edge of it the request read in order leaves out, each
/// set built on a copy of `model` with the same options, so that the
/// request without them builds (ADR 0001, amendment "checked blends under
/// overflow"). Empty where cancelled or where the request holds no edge.
fn refused_in_order(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    error: OgeomError,
    build: impl Fn(&mut Model, &[usize]) -> OgeomResult<Built>,
) -> BlendRefusal {
    if matches!(error, OgeomError::Cancelled) || edges.is_empty() {
        return BlendRefusal {
            edges: Vec::new(),
            error,
        };
    }
    match in_order(model, solid, edges, build, |_, _| None) {
        Ok((_, refused)) if !refused.is_empty() => BlendRefusal {
            edges: refused.iter().map(|r| r.edges[0]).collect(),
            error,
        },
        Err(e @ OgeomError::Cancelled) => BlendRefusal {
            edges: Vec::new(),
            error: e,
        },
        // Every edge accepted is the whole request, which did not build;
        // or its edges did not read: the whole request is named.
        _ => BlendRefusal {
            edges: (0..edges.len()).collect(),
            error,
        },
    }
}
