//! Chamfers: the bevel that replaces an edge.
//!
//! A chamfer along a straight edge between two planar faces is a wedge
//! subtracted, and the wedge's own faces lie *exactly* on the solid's
//! (coplanar, materials aligned), which is the boolean's same-domain case.
//!
//! Three spellings, one construction. The symmetric chamfer cuts the same
//! distance along both faces; the distance-distance form cuts a named
//! distance along a named face and another along its neighbour; the
//! distance-angle form cuts a distance along the named face and leaves it at
//! an angle, with the second distance derived where that bevel meets the
//! other face. All three end in the same wedge subtraction.
//!
//! Two seats per spelling. A straight edge between planes takes the
//! triangular prism above; the circular rim where a coaxial cylinder or cone
//! meets a perpendicular planar cap takes a revolved wedge whose bevel is a
//! *cone*, its leg on a cone's wall running along the slant:
//! the same flanks the rim fillet builds, with the quarter-tube exchanged
//! for the slant, and the same melt taking the legs away.

use crate::chamfer_sphere::{SphereSeat, sphere_bevel, sphere_seat};
use crate::support::{
    Blamed, RevolvedSeat, Seat, apply_wedge, edge_curve, planar_face, planar_seat, revolved_flanks,
    revolved_seat,
};
use ogeom_algo::{Built, make_revolution_band};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{ConeSurface, Curve, SurfaceGeometry};
use ogeom_math::{Cone, Direction, Frame};
use ogeom_topo::{Model, Shape};

/// Bevel a straight edge of a solid, cutting `distance` back along each of
/// its two faces.
///
/// The edge must be straight, convex, and shared by exactly two planar faces;
/// the distances are equal (the symmetric chamfer). The result is the boolean
/// difference with a wedge whose legs run along the two faces, so the
/// history reads as a cut: the two faces are modified into their trimmed
/// pieces, the edge's neighbourhood gains the bevel face.
///
/// A seam (an edge between two faces of one shell on one surface facing
/// one way, as a boolean leaves where it joined two flats) is no edge of
/// the part's form, and naming one names nothing: the solid comes back as
/// it is. See [`chamfer_edges_with`].
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, not convex, not shared by exactly two planar faces of
/// `solid`, or `distance` is not a usable length.
pub fn chamfer_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        std::slice::from_ref(edge),
        |model, solid, parts| chamfer_edge_own(model, solid, &parts[0], distance, tol),
    )
}

/// [`chamfer_edge`] on a shape every node below which is its own.
fn chamfer_edge_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if crate::engine::is_seam(model, solid, edge, tol)? {
        return Ok(unchanged(solid));
    }
    let specs = [(edge.clone(), Chamfer::Symmetric(distance))];
    if let Some(built) = engine_chamfer(model, solid, &specs, tol)? {
        return Ok(built);
    }
    wedge_for(model, solid, edge, &specs[0].1, tol)?.apply(model, solid, edge, tol)
}

/// Bevel a straight edge, cutting `on_face` back along `face` and `on_other`
/// along the edge's other face.
///
/// The asymmetric chamfer: `face` names which side the first distance applies
/// to, and must be one of the two faces meeting at the edge.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge`], and additionally if `face` is not one of the edge's
/// two faces.
pub fn chamfer_edge_distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            chamfer_edge_distances_own(model, solid, &parts[0], &parts[1], on_face, on_other, tol)
        },
    )
}

/// [`chamfer_edge_distances`] on a shape every node below which is its own.
fn chamfer_edge_distances_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spec = Chamfer::Distances {
        face: face.clone(),
        on_face,
        on_other,
    };
    if crate::engine::is_seam(model, solid, edge, tol)? {
        return Ok(unchanged(solid));
    }
    let specs = [(edge.clone(), spec)];
    if let Some(built) = engine_chamfer(model, solid, &specs, tol)? {
        return Ok(built);
    }
    wedge_for(model, solid, edge, &specs[0].1, tol)?.apply(model, solid, edge, tol)
}

/// Bevel a straight edge, cutting `distance` back along `face` and leaving it
/// at `angle` radians from that face.
///
/// The distance-angle chamfer: the second distance is where the bevel,
/// departing the named face at the given angle, meets the other face. An
/// angle of `π/4` on a square edge reproduces the symmetric chamfer.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge_distances`], and additionally if the bevel at that angle
/// never reaches the other face.
pub fn chamfer_edge_angle(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &[edge.clone(), face.clone()],
        |model, solid, parts| {
            chamfer_edge_angle_own(model, solid, &parts[0], &parts[1], distance, angle, tol)
        },
    )
}

/// [`chamfer_edge_angle`] on a shape every node below which is its own.
fn chamfer_edge_angle_own(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let spec = Chamfer::Angle {
        face: face.clone(),
        distance,
        angle,
    };
    if crate::engine::is_seam(model, solid, edge, tol)? {
        return Ok(unchanged(solid));
    }
    let specs = [(edge.clone(), spec)];
    if let Some(built) = engine_chamfer(model, solid, &specs, tol)? {
        return Ok(built);
    }
    wedge_for(model, solid, edge, &specs[0].1, tol)?.apply(model, solid, edge, tol)
}

/// The shapes chamfer specs name: each edge, then the face its chamfer
/// runs along where it names one.
fn chamfer_parts(specs: &[(Shape, Chamfer)]) -> Vec<Shape> {
    let mut parts = Vec::with_capacity(specs.len() * 2);
    for (edge, spec) in specs {
        parts.push(edge.clone());
        match spec {
            Chamfer::Symmetric(_) => {}
            Chamfer::Distances { face, .. } | Chamfer::Angle { face, .. } => {
                parts.push(face.clone());
            }
        }
    }
    parts
}

/// The specs again, naming `parts` in the order [`chamfer_parts`] lists
/// them.
fn chamfer_specs(specs: &[(Shape, Chamfer)], parts: &[Shape]) -> Vec<(Shape, Chamfer)> {
    let mut parts = parts.iter().cloned();
    let mut next = |was: &Shape| parts.next().unwrap_or_else(|| was.clone());
    specs
        .iter()
        .map(|(edge, spec)| {
            let edge = next(edge);
            let spec = match spec {
                Chamfer::Symmetric(d) => Chamfer::Symmetric(*d),
                Chamfer::Distances {
                    face,
                    on_face,
                    on_other,
                } => Chamfer::Distances {
                    face: next(face),
                    on_face: *on_face,
                    on_other: *on_other,
                },
                Chamfer::Angle {
                    face,
                    distance,
                    angle,
                } => Chamfer::Angle {
                    face: next(face),
                    distance: *distance,
                    angle: *angle,
                },
            };
            (edge, spec)
        })
        .collect()
}

/// One edge's chamfer, in any of the three spellings.
#[derive(Debug, Clone)]
pub enum Chamfer {
    /// The same distance back along both faces, as [`chamfer_edge`].
    Symmetric(f64),
    /// `on_face` along `face` and `on_other` along the other face, as
    /// [`chamfer_edge_distances`].
    Distances {
        /// The face the first distance runs along.
        face: Shape,
        /// The distance along `face`.
        on_face: f64,
        /// The distance along the edge's other face.
        on_other: f64,
    },
    /// `distance` along `face`, leaving it at `angle` radians, as
    /// [`chamfer_edge_angle`].
    Angle {
        /// The face the distance runs along.
        face: Shape,
        /// The distance along `face`.
        distance: f64,
        /// The bevel's angle from `face`, in radians.
        angle: f64,
    },
}

/// Bevel several edges of a solid as one operation, each by `distance`.
///
/// As [`chamfer_edges_with`] with the symmetric chamfer on every edge.
///
/// # Errors
///
/// As [`chamfer_edges_with`].
pub fn chamfer_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let specs: Vec<(Shape, Chamfer)> = edges
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(distance)))
        .collect();
    chamfer_edges_with(model, solid, &specs, tol)
}

/// Bevel several edges of a solid as one operation, each with its own
/// chamfer.
///
/// The whole request is read on the solid as it stands before the call, so
/// the result is the same for every order of `specs`. Where two bevels meet
/// at a vertex they mitre along the line their planes share, with neither a
/// step nor a cap, at a face's reflex corner as at a convex one. Where three
/// or more meet, the vertex is closed by the facet through the points where
/// their contact lines cross on the faces round it
/// ([`ChamferCorner::Facet`](crate::ChamferCorner::Facet)); the mitre of
/// their planes is [`chamfer_edges_cornered`](crate::chamfer_edges_cornered)
/// with [`ChamferCorner::Point`](crate::ChamferCorner::Point).
///
/// A solid bounded by planes is chamfered by the blend engine in one pass:
/// every face set back from its bevelled edges, every bevel the plane
/// between its two setback lines, every vertex of three bevels its facet.
/// Where a face's setback lines cross (a thin wall bevelled on both edges,
/// a plate bevelled all round past half its thickness), the face keeps what
/// lies inside every bevel, or goes, and its bevels meet along the ridges
/// their planes share: round such a face, and at a convex corner of two
/// bevels, what is left is the solid cut by the planes there. A request
/// past that (a consumed face with a reflex corner or a concave edge beside
/// it), and any other solid, are chamfered by wedges built on the solid as
/// given and applied together, mitred where bevels meet.
///
/// A seam is no edge of the part's form: an edge between two faces of one
/// shell lying on one surface facing one way
/// ([`ogeom_heal::same_domain`]), as a boolean leaves where it joined two
/// flats, the surface running on across it and the material on one side.
/// Mainstream CAD merges such faces after the boolean, so the edge cannot
/// be picked; here a spec naming one names nothing and is passed over, and
/// a request naming nothing but seams returns the solid as it is. A solid
/// holding seams is read by the blend engine with them dissolved, as
/// [`ogeom_heal::unify_same_domain`] and [`ogeom_heal::merge_edges`] merge
/// them: the faces on one plane come out one face, the pieces of an edge
/// split along its length one bevel, and the history tells of the solid
/// as given (each face modified into the face holding it, each seam and
/// each vertex the seams left inside an edge deleted). Coplanar faces
/// facing opposite ways, the material changing side across the edge
/// between them, meet at a real edge.
///
/// What `solid` shares with other shapes is copied before it is
/// worked on ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`chamfer_edge`], [`chamfer_edge_distances`] and
/// [`chamfer_edge_angle`] per edge, judged on the solid as it stands before
/// the call; and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if no edges are given.
pub fn chamfer_edges_with(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    chamfer_edges_blamed(model, solid, specs, tol).map_err(|blamed| blamed.error)
}

/// [`chamfer_edges_with`], its refusal naming the specs it concerns.
///
/// A solid bounded by planes goes through the blend engine
/// ([`crate::engine`]): one pass by setbacks, each vertex of three bevels
/// closed by a facet. Any other solid is chamfered by wedges, mitred.
pub(crate) fn chamfer_edges_blamed(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, Blamed> {
    if specs.is_empty() {
        return Err(
            ogeom_core::ogeom_err!(Construction, "a chain of no edges bevels nothing").into(),
        );
    }
    let (specs, index) = without_seams(model, solid, specs, tol)?;
    if specs.is_empty() {
        return Ok(unchanged(solid));
    }
    if let Some(built) = engine_chamfer(model, solid, &specs, tol)? {
        return Ok(built);
    }
    chamfer_edges_mitred_blamed(model, solid, &specs, tol).map_err(|mut blamed| {
        for i in &mut blamed.edges {
            *i = index[*i];
        }
        blamed
    })
}

/// Chamfer specs, each edge with its chamfer.
pub(crate) type Specs = Vec<(Shape, Chamfer)>;

/// `specs` without those naming a seam of `solid` (an edge between two
/// faces on one surface facing one way, no edge of the part's form), with
/// the index in `specs` of each kept.
pub(crate) fn without_seams(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<(Specs, Vec<usize>)> {
    let mut kept = Vec::with_capacity(specs.len());
    let mut index = Vec::with_capacity(specs.len());
    for (i, spec) in specs.iter().enumerate() {
        if !crate::engine::is_seam(model, solid, &spec.0, tol)? {
            kept.push(spec.clone());
            index.push(i);
        }
    }
    Ok((kept, index))
}

/// `solid` as given: what a request naming nothing but seams makes.
pub(crate) fn unchanged(solid: &Shape) -> Built {
    let mut built = Built::from_nothing(solid.clone());
    built.history.modify(solid, solid.clone());
    built
}

/// The chamfer of `specs` by the blend engine, or `None` where the solid
/// or a spec lies outside its classes. A solid holding seams is chamfered
/// with them dissolved, as [`ogeom_heal::unify_same_domain`] merges them;
/// no spec may name a seam.
pub(crate) fn engine_chamfer(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if model.kind_of(solid)? == ogeom_topo::ShapeType::Compound {
        return engine_chamfer_compound(model, solid, specs, tol);
    }
    // A solid with placed shapes in it (a prism's far end is its profile
    // moved, the same nodes standing twice) is read with its placements
    // baked in, the request traced onto the baked shapes and the history
    // carried back.
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    if baked.shape.is_same(solid) {
        return engine_chamfer_unplaced(model, solid, specs, tol);
    }
    let Some(traced) = traced_specs(&baked.history, specs) else {
        return Ok(None);
    };
    let Some(built) = engine_chamfer_unplaced(model, &baked.shape, &traced, tol)? else {
        return Ok(None);
    };
    Ok(Some(Built::new(
        built.shape,
        baked.history.then(&built.history),
    )))
}

/// [`engine_chamfer`] closing each vertex of three chamfers by their mitre:
/// on a solid bounded by planes the bevels run on to the point their
/// planes share, or, at a convex vertex of four or more sharing none, until
/// the others cut them; on curved faces to the point their surfaces share,
/// meeting each other along their sections from there. `None` elsewhere,
/// the model then left as it was found.
///
/// # Errors
///
/// As [`engine_chamfer`].
pub(crate) fn engine_chamfer_mitred(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let before = model.clone();
    let built = engine_chamfer_mitred_attempt(model, solid, specs, tol);
    if !matches!(built, Ok(Some(_))) {
        *model = before;
    }
    built
}

/// [`engine_chamfer_mitred`]'s attempt, on the live model.
fn engine_chamfer_mitred_attempt(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    let (solid, specs) = if baked.shape.is_same(solid) {
        (solid.clone(), specs.to_vec())
    } else {
        let Some(traced) = traced_specs(&baked.history, specs) else {
            return Ok(None);
        };
        (baked.shape.clone(), traced)
    };
    if let Some(what) = crate::engine::degenerate_planar(model, &solid)? {
        return Err(degenerate_refused(&what));
    }
    if !matches!(
        crate::engine::dissolve(model, &solid, tol)?,
        crate::engine::Dissolved::Plain
    ) {
        return Ok(None);
    }
    if let Some(mut request) = crate::engine::Request::from_chamfers(model, &solid, &specs, tol)? {
        request.mitre = true;
        if let Some(built) = crate::engine::run(model, &solid, &request, tol)? {
            return Ok(Some(Built::new(
                built.shape,
                baked.history.then(&built.history),
            )));
        }
    }
    // A solid with curved faces: mitred on the curved path, its seams
    // dissolved where it holds any.
    if crate::engine::bounded_by_planes(model, &solid, tol)? {
        return Ok(None);
    }
    let built = if let Some(seamless) = crate::engine::dissolve_curved(model, &solid, tol)? {
        if seamless.asked_apart(&specs) {
            return Err(apart_refused());
        }
        let Some(specs) = seamless.specs(&specs) else {
            return Ok(None);
        };
        let Some(built) =
            crate::engine::chamfer_curved_mitred(model, &seamless.solid, &specs, tol)?
        else {
            return Ok(None);
        };
        Built::new(built.shape, seamless.carry(&built.history))
    } else {
        let Some(built) = crate::engine::chamfer_curved_mitred(model, &solid, &specs, tol)? else {
            return Ok(None);
        };
        built
    };
    Ok(Some(Built::new(
        built.shape,
        baked.history.then(&built.history),
    )))
}

/// [`engine_chamfer`] on a compound of solids: each solid with the specs
/// whose edges it holds, the rest kept as they are, the results gathered
/// into a compound again. `None` where an edge lies in no solid or in two,
/// or any solid with specs is handed back.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the compound holds anything but solids, or nothing.
fn engine_chamfer_compound(
    model: &mut Model,
    compound: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_topo::ShapeType;
    let parts = model.children_of(compound)?;
    // A chamfer cuts solids: a compound holding anything else (a shell, a
    // face, a nested compound) or nothing is refused by name (ADR 0001,
    // amendment "what M6 refuses by name").
    if parts.is_empty()
        || parts
            .iter()
            .any(|p| !matches!(model.kind_of(p), Ok(ShapeType::Solid)))
    {
        ogeom_core::ogeom_bail!(
            Construction,
            "a chamfer cuts solids; this compound holds something other than solids"
        );
    }
    let mut asked: Vec<Vec<(Shape, Chamfer)>> = vec![Vec::new(); parts.len()];
    let edges: Vec<Vec<Shape>> = parts
        .iter()
        .map(|p| ogeom_topo::explore_unique(model, p, ShapeType::Edge))
        .collect::<OgeomResult<_>>()?;
    for (edge, spec) in specs {
        let holders: Vec<usize> = (0..parts.len())
            .filter(|k| edges[*k].iter().any(|e| e.is_same(edge)))
            .collect();
        let [k] = holders[..] else {
            return Ok(None);
        };
        asked[k].push((edge.clone(), spec.clone()));
    }
    let mut out = Vec::with_capacity(parts.len());
    let mut history = ogeom_algo::History::new();
    for (part, own) in parts.iter().zip(&asked) {
        if own.is_empty() {
            out.push(part.clone());
            continue;
        }
        let Some(built) = engine_chamfer(model, part, own, tol)? else {
            return Ok(None);
        };
        history = history.then(&built.history);
        out.push(built.shape);
    }
    let gathered = ogeom_algo::make_compound(model, &out)?.shape;
    history.modify(compound, gathered.clone());
    Ok(Some(Built::new(gathered, history)))
}

/// [`engine_chamfer`] on a solid with voids: the outer shell chamfered as a
/// solid alone where it holds edges of the request, and each void holding
/// some chamfered as the solid it bounds turned inside out (a void's
/// edges are concave in the solid, and filling a wedge of the cavity is
/// taking it off that solid), then the shells put together again. `None`
/// where an edge lies on no shell or two, a shell's chamfer is handed back,
/// or a void touches the outer shell or another void or lies outside the
/// outer shell (in what a bevel takes off).
fn engine_chamfer_voided(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_topo::ShapeType;
    let shells = model.children_of(solid)?;
    // The outer shell: the one whose bounds hold every other's.
    let mut sized = Vec::with_capacity(shells.len());
    for shell in &shells {
        sized.push(ogeom_algo::shape_bounds(model, shell, tol)?.diagonal());
    }
    let outer = (0..shells.len())
        .max_by(|a, b| sized[*a].total_cmp(&sized[*b]))
        .unwrap_or(0);
    let edges: Vec<Vec<Shape>> = shells
        .iter()
        .map(|shell| ogeom_topo::explore_unique(model, shell, ShapeType::Edge))
        .collect::<OgeomResult<_>>()?;
    let mut asked: Vec<Vec<(Shape, Chamfer)>> = vec![Vec::new(); shells.len()];
    for (edge, spec) in specs {
        let holders: Vec<usize> = (0..shells.len())
            .filter(|k| edges[*k].iter().any(|e| e.is_same(edge)))
            .collect();
        let [k] = holders[..] else {
            return Ok(None);
        };
        asked[k].push((edge.clone(), spec.clone()));
    }
    let mut history = ogeom_algo::History::new();
    // Each shell as it comes out: the outer one first.
    let mut order: Vec<usize> = vec![outer];
    order.extend((0..shells.len()).filter(|k| *k != outer));
    let mut kept: Vec<Shape> = Vec::with_capacity(shells.len());
    for k in order {
        if asked[k].is_empty() {
            kept.push(shells[k].clone());
            continue;
        }
        // The shell as a solid of its own, the outer one as it stands and
        // a void turned inside out: whichever way bounds a valid solid.
        let mut alone = None;
        for shell in [shells[k].clone(), shells[k].reversed()] {
            let candidate = ogeom_algo::make_solid(model, &[shell])?.shape;
            if ogeom_algo::check(model, &candidate, tol)?.is_valid()
                && ogeom_algo::volume_properties(
                    model,
                    &candidate,
                    ogeom_mesh::Deflection::default(),
                    tol,
                )?
                .mass
                    > 0.0
            {
                alone = Some(candidate);
                break;
            }
        }
        let Some(alone) = alone else {
            return Ok(None);
        };
        let Some(built) = engine_chamfer(model, &alone, &asked[k], tol)? else {
            return Ok(None);
        };
        let [shell] = &model.children_of(&built.shape)?[..] else {
            return Ok(None);
        };
        let shell = shell.clone();
        // Back the way the solid held it.
        let turned = model
            .children_of(&alone)?
            .first()
            .is_some_and(|s| s.orientation() != shells[k].orientation());
        kept.push(if turned { shell.reversed() } else { shell });
        history = history.then(&built.history);
    }
    // Every void apart from the outer shell and from every other, and
    // inside the outer shell. A bevel reaching a void is non-local (ADR
    // 0001, amendment "M5a"): what it takes there is already empty, so the
    // solid is the outer shell's chamfered body less each void's, by one
    // boolean a void.
    let outside = ogeom_algo::make_solid(model, &[kept[0].clone()])?.shape;
    for i in 1..kept.len() {
        for j in 0..i {
            let apart = ogeom_algo::distance_between_shapes(
                model,
                &kept[i],
                &kept[j],
                ogeom_intersect::ExtremaOptions::default(),
                tol,
            )?;
            if apart.distance <= tol.confusion() {
                return voids_cut(model, solid, &kept, history, tol);
            }
        }
        let Some(corner) = ogeom_topo::explore_unique(model, &kept[i], ShapeType::Vertex)?
            .first()
            .cloned()
        else {
            return Ok(None);
        };
        let at = crate::setback::point_of(model, &corner)?;
        if ogeom_algo::classify_in_solid_exact(model, &outside, at, tol)?
            != ogeom_algo::Containment::In
        {
            return Ok(None);
        }
    }
    let whole = ogeom_algo::make_solid(model, &kept)?.shape;
    history.modify(solid, whole.clone());
    Ok(Some(Built::new(whole, history)))
}

/// The outer shell's body less each void's, one boolean a void, where a
/// bevel reaches a void (amendment "M5a"): each void's shell turned to face
/// out bounds the void's body. `history` is the shells' chamfers', carried
/// on through the booleans. `None` where a boolean refuses.
///
/// # Errors
///
/// As the solid builder and the booleans, other than their refusals.
fn voids_cut(
    model: &mut Model,
    solid: &Shape,
    kept: &[Shape],
    history: ogeom_algo::History,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut body = ogeom_algo::make_solid(model, &[kept[0].clone()])?.shape;
    let mut history = history;
    for void in &kept[1..] {
        let mut hollow = None;
        for shell in [void.reversed(), void.clone()] {
            let candidate = ogeom_algo::make_solid(model, &[shell])?.shape;
            if ogeom_algo::volume_properties(
                model,
                &candidate,
                ogeom_mesh::Deflection::default(),
                tol,
            )?
            .mass
                > 0.0
            {
                hollow = Some(candidate);
                break;
            }
        }
        let Some(hollow) = hollow else {
            return Ok(None);
        };
        let Ok(step) = ogeom_bool::cut(model, &body, &hollow, tol) else {
            return Ok(None);
        };
        body = step.shape.clone();
        history = history.then(&step.history);
    }
    history.modify(solid, body.clone());
    Ok(Some(Built::new(body, history)))
}

/// `specs` traced through `history`: each edge, and each face a spec names,
/// to the one shape it went to; `None` where any went to none or several.
fn traced_specs(
    history: &ogeom_algo::History,
    specs: &[(Shape, Chamfer)],
) -> Option<Vec<(Shape, Chamfer)>> {
    let one = |s: &Shape| match history.trace(s) {
        [only] => Some(only.clone()),
        _ => None,
    };
    specs
        .iter()
        .map(|(edge, spec)| {
            let spec = match spec {
                Chamfer::Symmetric(d) => Chamfer::Symmetric(*d),
                Chamfer::Distances {
                    face,
                    on_face,
                    on_other,
                } => Chamfer::Distances {
                    face: one(face)?,
                    on_face: *on_face,
                    on_other: *on_other,
                },
                Chamfer::Angle {
                    face,
                    distance,
                    angle,
                } => Chamfer::Angle {
                    face: one(face)?,
                    distance: *distance,
                    angle: *angle,
                },
            };
            Some((one(edge)?, spec))
        })
        .collect()
}

/// [`engine_chamfer`] on a solid with nothing placed in it. Dissolving
/// seams rewrites nodes the solid shares, so where the engine hands the
/// request back the model is left as it was found, for the caller's own
/// path to read.
fn engine_chamfer_unplaced(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let before = model.clone();
    let built = engine_chamfer_attempt(model, solid, specs, tol);
    if !matches!(built, Ok(Some(_))) {
        *model = before;
    }
    built
}

/// [`engine_chamfer_unplaced`]'s attempt, on the live model.
fn engine_chamfer_attempt(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use crate::engine::Dissolved;
    if model.kind_of(solid)? == ogeom_topo::ShapeType::Solid && model.children_of(solid)?.len() > 1
    {
        return engine_chamfer_voided(model, solid, specs, tol);
    }
    // Degenerate planar topology is refused by name, not chamfered past.
    if let Some(what) = crate::engine::degenerate_planar(model, solid)? {
        return Err(degenerate_refused(&what));
    }
    match crate::engine::dissolve(model, solid, tol)? {
        Dissolved::Plain => {
            if let Some(request) = crate::engine::Request::from_chamfers(model, solid, specs, tol)?
                && let Some(built) = crate::engine::run(model, solid, &request, tol)?
            {
                return Ok(Some(built));
            }
            // A solid with curved faces: its seats on cylinders and cones,
            // read with its seams dissolved where it holds any.
            // On planes, only where a corner's facet is off one plane: the
            // curved path closes it with a patch.
            if crate::engine::bounded_by_planes(model, solid, tol)? {
                return crate::engine::chamfer_patched(model, solid, specs, tol);
            }
            // Straight bevels between planes reaching a face they do not end
            // at: the one boolean (ADR 0001, amendment "M5a").
            if let Some(built) = crate::engine::straight_non_local(model, solid, specs, tol)? {
                return Ok(Some(built));
            }
            if let Some(seamless) = crate::engine::dissolve_curved(model, solid, tol)? {
                // Pieces of a curved edge asked apart step between them,
                // each piece's rim tool ending on the plane through the
                // joint (amendment "what M6 refuses by name").
                if seamless.asked_apart(specs) {
                    return match rims_one_boolean(model, solid, specs, tol)? {
                        Some(built) => Ok(Some(built)),
                        None => Err(apart_refused()),
                    };
                }
                let Some(asked) = seamless.specs(specs) else {
                    return Ok(None);
                };
                if let Some(built) =
                    crate::engine::chamfer_curved(model, &seamless.solid, &asked, tol)?
                {
                    return Ok(Some(Built::new(
                        built.shape,
                        seamless.carry(&built.history),
                    )));
                }
                return rims_one_boolean(model, solid, specs, tol);
            }
            if let Some(built) = crate::engine::chamfer_curved(model, solid, specs, tol)? {
                return Ok(Some(built));
            }
            rims_one_boolean(model, solid, specs, tol)
        }
        Dissolved::Outside => Ok(None),
        Dissolved::Into(seamless) => {
            if seamless.asked_apart(specs) {
                return engine_chamfer_apart(model, &seamless, specs, tol);
            }
            let Some(specs) = seamless.specs(specs) else {
                return Ok(None);
            };
            let Some(request) =
                crate::engine::Request::from_chamfers(model, &seamless.solid, &specs, tol)?
            else {
                return Ok(None);
            };
            let Some(built) = crate::engine::run(model, &seamless.solid, &request, tol)? else {
                return Ok(None);
            };
            Ok(Some(Built::new(
                built.shape,
                seamless.carry(&built.history),
            )))
        }
    }
}

/// The chamfer of `specs` where pieces of one edge the seams split are
/// asked with different chamfers: each such edge is read as its pieces,
/// each piece chamfered as asked, and the bevels step at the vertex
/// between two pieces, each ending square to the edge there (as two
/// chamfer features on the pieces would build; ADR 0001, amendment "what
/// M6 refuses by name").
fn engine_chamfer_apart(
    model: &mut Model,
    seamless: &crate::engine::Seamless,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some((whole, pieces, splits)) = seamless.apart(model, specs)? else {
        return Ok(None);
    };
    let Some(mut request) =
        crate::engine::Request::from_chamfers(model, &seamless.solid, &whole, tol)?
    else {
        return Ok(None);
    };
    for (piece, edge, spec) in &pieces {
        let Some(mut one) = crate::engine::Request::from_chamfers(
            model,
            &seamless.solid,
            &[(edge.clone(), spec.clone())],
            tol,
        )?
        else {
            return Ok(None);
        };
        let Some(mut blend) = one.blends.pop() else {
            return Ok(None);
        };
        blend.edge = piece.clone();
        request.blends.push(blend);
    }
    let named: Vec<Shape> = splits
        .iter()
        .flat_map(|s| s.pieces.iter().chain(&s.joints).cloned())
        .collect();
    request.splits = splits;
    let Some(built) = crate::engine::run(model, &seamless.solid, &request, tol)? else {
        return Ok(None);
    };
    // The pieces and joints are the input's own: their history passes
    // through as the engine made it.
    let mut history = seamless.carry(&built.history);
    for shape in &named {
        if built.history.is_deleted(shape) {
            history.delete(shape);
        }
        for made in built.history.generated(shape) {
            history.generate(shape, made.clone());
        }
        for made in built.history.modified(shape) {
            history.modify(shape, made.clone());
        }
    }
    Ok(Some(Built::new(built.shape, history)))
}

/// Rim chamfers the curved path's rings do not place (a rim rolling over
/// its wall's far side as a notch, a rim reaching a face it does not end
/// at, pieces of one rim asked apart), built by the engine's one boolean
/// (ADR 0001, amendments "M5a" and "M5b"): each rim's tool is its bevel's
/// meridian triangle turned about the axis over the rim's own span, its
/// ends the planes through the axis at an arc's ends (the step at a joint).
/// The tool's face on the bevel's cone is named for its edge. Only where
/// every asked edge is a rim between a planar cap and a coaxial cylinder
/// or cone, and no two meet but pieces of one circle; `None` otherwise or
/// where a boolean refuses or a bevel is left nothing of.
///
/// # Errors
///
/// As reading the seats (a leg past the axis or a cone's apex is refused
/// by name) and the booleans, other than their refusals.
fn rims_one_boolean(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if specs.is_empty() {
        return Ok(None);
    }
    let mut seats: Vec<RevolvedSeat> = Vec::with_capacity(specs.len());
    for (edge, spec) in specs {
        let Ok(Legs::Rim(seat, _, _)) = legs(model, solid, edge, spec, tol) else {
            return Ok(None);
        };
        seats.push(seat);
    }
    // No two asked edges meet, but pieces of one circle.
    for (i, (a, _)) in specs.iter().enumerate() {
        for (j, (b, _)) in specs.iter().enumerate().skip(i + 1) {
            let ends = |e: &Shape| -> OgeomResult<Vec<Shape>> {
                Ok(ogeom_algo::edge_vertices(model, e)?
                    .map(|(u, v)| vec![u, v])
                    .unwrap_or_default())
            };
            let (ea, eb) = (ends(a)?, ends(b)?);
            let meet = ea.iter().any(|u| eb.iter().any(|v| u.is_same(v)));
            let one_circle = seats[i].centre.distance(seats[j].centre) <= tol.confusion()
                && (seats[i].radius - seats[j].radius).abs() <= tol.confusion()
                && seats[i].up.dot(seats[j].up) >= 1.0 - tol.angular();
            if meet && !one_circle {
                return Ok(None);
            }
        }
    }
    let mut tools: Vec<crate::engine::Tool> = Vec::with_capacity(specs.len());
    for ((edge, spec), seat) in specs.iter().zip(&seats) {
        let wedge = wedge_for(model, solid, edge, spec, tol)?;
        let shape = match &wedge.solid {
            Some(s) => s.clone(),
            None => match crate::support::wedge_solid(model, &wedge.faces, true, tol) {
                Ok(s) => s,
                Err(_) => return Ok(None),
            },
        };
        // The tool's face on the bevel: neither on the cap's plane, nor on
        // the wall, nor a plane through the axis.
        let mut named = Vec::new();
        for face in ogeom_topo::explore_unique(model, &shape, ogeom_topo::ShapeType::Face)? {
            let (at, n) = ogeom_algo::face_normal(model, &face, tol)?;
            let off_axis = at - seat.centre;
            let h = off_axis.dot(seat.up);
            let rho = (off_axis - seat.up * h).magnitude();
            let on_cap = h.abs() <= tol.confusion() && n.dot(seat.up).abs() >= 1.0 - tol.angular();
            let on_wall = (rho - seat.radius_at(h)).abs() <= tol.confusion();
            let meridian = n.dot(seat.up).abs() <= tol.angular()
                && n.dot(off_axis - seat.up * h).abs() <= tol.confusion() * rho.max(1.0);
            if !(on_cap || on_wall || meridian) {
                named.push((face, edge.clone()));
            }
        }
        tools.push(crate::engine::Tool {
            solid: shape,
            sense: if wedge.additive {
                crate::engine::Sense::Fuse
            } else {
                crate::engine::Sense::Cut
            },
            faces: named,
        });
    }
    let edges: Vec<Shape> = specs.iter().map(|(e, _)| e.clone()).collect();
    crate::engine::apply_tools(model, solid, &edges, &tools, tol)
}

/// The chamfer of `specs` by wedges built on the solid as given and
/// applied together, mitred where bevels meet: three bevels at a convex
/// vertex meet at the point their planes share.
pub(crate) fn chamfer_edges_mitred(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    chamfer_edges_mitred_blamed(model, solid, specs, tol).map_err(|blamed| blamed.error)
}

/// [`chamfer_edges_mitred`], its refusal naming the specs it concerns.
fn chamfer_edges_mitred_blamed(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, Blamed> {
    let mut blame = Vec::new();
    ogeom_algo::on_own_nodes(
        model,
        solid,
        &chamfer_parts(specs),
        |model, solid, parts| {
            chamfer_edges_with_own(model, solid, &chamfer_specs(specs, parts), tol).map_err(
                |blamed| {
                    blame = blamed.edges;
                    blamed.error
                },
            )
        },
    )
    .map_err(|error| Blamed {
        edges: blame,
        error,
    })
}

/// [`chamfer_edges_with`] on a shape every node below which is its own.
///
/// A wedge refused names its spec; one that does not apply names its
/// spec, then the specs before it whose edges share a vertex with its
/// edge, whose bevels it was applied after.
fn chamfer_edges_with_own(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> Result<Built, Blamed> {
    if specs.is_empty() {
        return Err(
            ogeom_core::ogeom_err!(Construction, "a chain of no edges bevels nothing").into(),
        );
    }
    let mut wedges = Vec::with_capacity(specs.len());
    for (i, (edge, spec)) in specs.iter().enumerate() {
        wedges.push(
            wedge_for(model, solid, edge, spec, tol).map_err(|error| Blamed {
                edges: vec![i],
                error,
            })?,
        );
    }
    // A bevel on a marched seat is built against the solid as it stands
    // when its turn comes: its edge is found again there, between the faces
    // it lay between as asked.
    let marched = specs
        .iter()
        .any(|(edge, _)| matches!(seat_kind(model, edge, tol), Ok(SeatKind::Marched)));
    let given = if marched {
        let mut given = Vec::with_capacity(specs.len());
        for (edge, _) in specs {
            given.push(faces_holding(model, solid, edge)?);
        }
        given
    } else {
        Vec::new()
    };
    let mut built = Built::from_nothing(solid.clone());
    let mut rings: Vec<(usize, WholeRing)> = Vec::new();
    for (i, ((edge, _), wedge)) in specs.iter().zip(wedges).enumerate() {
        // The arcs a boolean split one sphere rim into each chamfer the
        // whole ring: the first cuts it, the others asking the same legs
        // are taken with it, and others asking other legs contradict it.
        if let Some(ring) = wedge.whole_ring {
            if let Some((j, first)) = rings.iter().find(|(_, r)| r.same_circle(&ring, tol)) {
                if first.same_legs(&ring, tol) {
                    continue;
                }
                return Err(Blamed {
                    edges: vec![i, *j],
                    error: ogeom_core::ogeom_err!(
                        Construction,
                        "two arcs of one sphere rim ask for different chamfers \
                         of its whole ring"
                    ),
                });
            }
            rings.push((i, ring));
        }
        let applied = match &wedge.marched {
            Some(spec) => marched_step(model, &built, edge, spec, &given[i], &specs[..i], tol),
            None => wedge.apply(model, &built.shape, edge, tol),
        };
        let step = match applied {
            Ok(step) => step,
            Err(error) => {
                let mut edges = vec![i];
                let ends = |e: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
                    Ok(ogeom_algo::edge_vertices(model, e)?
                        .map(|(a, b)| vec![a.node(), b.node()])
                        .unwrap_or_default())
                };
                let own = ends(edge)?;
                for (j, (other, _)) in specs[..i].iter().enumerate() {
                    if ends(other)?.iter().any(|v| own.contains(v)) {
                        edges.push(j);
                    }
                }
                return Err(Blamed { edges, error });
            }
        };
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// The seat of a rim's bevel: a revolved seat, its wall a cylinder or a
/// cone, or a sphere centred on the rim's axis.
enum RimSeat {
    Revolved(RevolvedSeat),
    Sphere(SphereSeat),
}

fn rim_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    rim: &ogeom_geom::CircleCurve,
    tol: Tolerances,
) -> OgeomResult<RimSeat> {
    match revolved_seat(model, solid, edge, rim, tol) {
        Ok(seat) => Ok(RimSeat::Revolved(seat)),
        Err(refused) => match sphere_seat(model, solid, edge, rim, tol)? {
            Some(seat) => Ok(RimSeat::Sphere(seat)),
            None => Err(refused),
        },
    }
}

/// How far along the axis a leg of `on_wall` down the wall's ruling
/// reaches: the leg itself on a cylinder, its run along the slant on a
/// cone.
pub(crate) fn wall_depth(seat: &RevolvedSeat, on_wall: f64) -> f64 {
    if seat.slope == 0.0 {
        on_wall
    } else {
        on_wall / seat.slope.hypot(1.0)
    }
}

/// The angle at the rim between the cap's leg and the wall's, in the
/// rim's section: square on a cylinder, and on a cone the cone's lean
/// from it, obtuse where the wall leans out under the cap.
fn rim_corner(seat: &RevolvedSeat) -> f64 {
    // In the section, the cap's leg runs `-sigma tau` radially and the
    // wall's runs `(-slope tau, -tau)`, both from the rim.
    (seat.sigma * seat.slope / seat.slope.hypot(1.0)).acos()
}

/// The wedge one chamfer cuts, built on `solid` as it stands.
fn wedge_for(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    match legs(model, solid, edge, spec, tol)? {
        Legs::Straight(seat, distances) => bevel(model, &seat, distances, tol),
        Legs::Rim(seat, on_wall, on_cap) => {
            revolved_bevel(model, solid, edge, &seat, on_wall, on_cap, tol)
        }
        Legs::Marched => Ok(Wedge {
            faces: Vec::new(),
            solid: None,
            additive: false,
            whole_ring: None,
            marched: Some(spec.clone()),
        }),
        Legs::Sphere(seat, on_wall, on_cap) => {
            // Each contact stays on its face at every station the walk
            // reads, or the wedge would cut past the face.
            let per = |face: &Shape| sphere_leg(&seat, face, on_wall, on_cap);
            if crate::fits::legs_largest(model, solid, edge, 1.0, per, tol)? < 1.0 {
                ogeom_bail!(
                    Construction,
                    "a chamfer on the rim between a plane and a sphere runs past \
                     one of its faces"
                );
            }
            let closed =
                ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
            let Curve::Circle(rim) = edge_curve(model, edge, tol)?.0 else {
                ogeom_bail!(Construction, "a sphere's rim is not a circle");
            };
            let whole = closed
                || crate::chamfer_sphere::rim_closes(model, solid, &seat, &rim.circle(), tol)?;
            // An arc that stops short of its circle is the section turned
            // through the arc's own sweep.
            let span = if whole {
                None
            } else {
                let frame =
                    Frame::new(seat.centre, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
                Some(arc_window(model, edge, &frame, tol)?)
            };
            Ok(Wedge {
                faces: sphere_bevel(model, &seat, on_wall, on_cap, span, tol)?,
                solid: None,
                additive: seat.tau < 0.0,
                whole_ring: whole.then_some(WholeRing {
                    centre: seat.centre,
                    radius: seat.radius,
                    up: seat.up,
                    on_wall,
                    on_cap,
                }),
                marched: None,
            })
        }
    }
}

/// How far a sphere seat's bevel runs along `face`: its leg on the wall or
/// the cap, nothing on any other face.
fn sphere_leg(seat: &SphereSeat, face: &Shape, on_wall: f64, on_cap: f64) -> f64 {
    if face.is_same(&seat.wall_face) {
        on_wall
    } else if face.is_same(&seat.cap_face) {
        on_cap
    } else {
        0.0
    }
}

/// The largest measure up to `up_to` at which a chamfer of `spec`'s
/// spelling on `edge` of `solid` stays seated on its faces, its contacts
/// short of their faces' far sides: for [`Chamfer::Symmetric`] the distance, for
/// [`Chamfer::Distances`] the largest common scale of its two distances,
/// and for [`Chamfer::Angle`] the distance along its face at its angle;
/// the measure `spec` itself names is not read. None where no measure
/// leaving both legs longer than the confusion distance fits.
///
/// Read from the seat, without building a wedge. On a straight edge
/// between planes, a convex bevel's leg runs at most as far as its face
/// reaches back from the edge, and the answer is exact against it; past
/// it [`chamfer_edges_with`] rolls the bevel over its face's far side as a
/// notch (ADR 0001, amendment "M5b"); a concave one fuses and runs past no
/// face. On a rim, the cap's leg stops short of the axis, and
/// each leg's contact stays on its face at the stations
/// [`largest_fillet_radius`](crate::largest_fillet_radius) reads, a
/// confusion distance short.
///
/// # Errors
///
/// As [`chamfer_edges_with`] where the seat or spelling is refused at any
/// measure, and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if `up_to` is not a usable length; as
/// [`largest_fillet_radius`](crate::largest_fillet_radius) where a
/// contact's walk on a rim depends on the measure.
pub fn largest_chamfer_distance(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let conf = tol.confusion();
    if !up_to.is_finite() || up_to <= conf {
        ogeom_bail!(
            Construction,
            "a largest distance up to {up_to} asks for nothing"
        );
    }
    // The legs at a measure of one, which scale with it.
    let unit = match spec {
        Chamfer::Symmetric(_) => Chamfer::Symmetric(1.0),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => Chamfer::Distances {
            face: face.clone(),
            on_face: *on_face,
            on_other: *on_other,
        },
        Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
            face: face.clone(),
            distance: 1.0,
            angle: *angle,
        },
    };
    let mut largest = up_to;
    let shortest = match legs(model, solid, edge, &unit, tol)? {
        Legs::Straight(seat, per) => {
            if seat.convex {
                for (i, leg) in per.into_iter().enumerate() {
                    if leg > 0.0 {
                        let along = seat.leg(i, tol)?;
                        let reach = crate::support::face_reach(
                            model,
                            &seat.faces[i],
                            seat.start,
                            along,
                            tol,
                        )?;
                        largest = largest.min(reach / leg);
                    }
                }
            }
            per[0].min(per[1])
        }
        Legs::Rim(seat, on_wall, on_cap) => {
            // The cap's leg ends at the radius `radius - sigma tau on_cap`,
            // which stays above the confusion distance.
            let inward = seat.sigma * seat.tau * on_cap;
            if inward > 0.0 {
                largest = largest.min((seat.radius - 2.0 * conf) / inward);
            }
            // A cone's wall leg ends at the radius the wall has a depth
            // down, which narrows toward the apex with the measure.
            let narrowing = seat.slope * seat.tau * wall_depth(&seat, on_wall);
            if narrowing > 0.0 {
                largest = largest.min((seat.radius - 2.0 * conf) / narrowing);
            }
            let per = |face: &Shape| {
                if face.is_same(&seat.wall_face) {
                    on_wall
                } else if face.is_same(&seat.cap_face) {
                    on_cap
                } else {
                    0.0
                }
            };
            largest = largest.min(crate::fits::legs_largest(
                model, solid, edge, largest, per, tol,
            )?);
            on_wall.min(on_cap)
        }
        Legs::Marched => ogeom_bail!(
            NotDone,
            "the largest distance of a chamfer on a marched seat is not read \
             from the seat"
        ),
        Legs::Sphere(seat, on_wall, on_cap) => {
            // The derived leg of an angle on a sphere's rim does not scale
            // with the measure: the wall's section is an arc.
            if matches!(spec, Chamfer::Angle { .. }) {
                ogeom_bail!(
                    NotDone,
                    "the largest distance of a chamfer at an angle on a sphere's \
                     rim is not read from the seat"
                );
            }
            let inward = seat.sigma * seat.tau * on_cap;
            if inward > 0.0 {
                largest = largest.min((seat.radius - 2.0 * conf) / inward);
            }
            largest = largest.min((seat.to_pole() - 2.0 * conf) / on_wall);
            let per = |face: &Shape| sphere_leg(&seat, face, on_wall, on_cap);
            largest = largest.min(crate::fits::legs_largest(
                model, solid, edge, largest, per, tol,
            )?);
            on_wall.min(on_cap)
        }
    };
    Ok((largest * shortest > conf).then_some(largest))
}

/// The seat a chamfer stands on and how far its legs run along each face.
pub(crate) enum Legs {
    /// A straight edge between planes, its legs along the seat's faces in
    /// order.
    Straight(Seat, [f64; 2]),
    /// A rim on a cylinder or a cone: its legs down the wall and in along the
    /// cap.
    Rim(RevolvedSeat, f64, f64),
    /// A rim on a sphere: its legs down the sphere's great circle and along
    /// the cap.
    Sphere(SphereSeat, f64, f64),
    /// A crease on a marched seat, its legs read station by station when
    /// the bevel is built.
    Marched,
}

/// Read one chamfer's seat and legs from `solid` as it stands.
pub(crate) fn legs(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Legs> {
    match spec {
        Chamfer::Symmetric(distance) => match seat_kind(model, edge, tol)? {
            SeatKind::Straight => {
                let seat = planar_seat(model, solid, edge, tol)?;
                Ok(Legs::Straight(seat, [*distance, *distance]))
            }
            SeatKind::Marched => Ok(Legs::Marched),
            SeatKind::Rim(rim) => match rim_seat(model, solid, edge, &rim, tol)? {
                RimSeat::Revolved(seat) => Ok(Legs::Rim(seat, *distance, *distance)),
                RimSeat::Sphere(seat) => Ok(Legs::Sphere(seat, *distance, *distance)),
            },
        },
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match seat_kind(model, edge, tol)? {
            SeatKind::Straight => {
                let seat = planar_seat(model, solid, edge, tol)?;
                let i = seat_side(&seat, face)?;
                let mut distances = [0.0; 2];
                distances[i] = *on_face;
                distances[1 - i] = *on_other;
                Ok(Legs::Straight(seat, distances))
            }
            SeatKind::Marched => Ok(Legs::Marched),
            SeatKind::Rim(rim) => match rim_seat(model, solid, edge, &rim, tol)? {
                RimSeat::Revolved(seat) => {
                    let (on_wall, on_cap) =
                        revolved_side(&seat.wall_face, &seat.cap_face, face, *on_face, *on_other)?;
                    Ok(Legs::Rim(seat, on_wall, on_cap))
                }
                RimSeat::Sphere(seat) => {
                    let (on_wall, on_cap) =
                        revolved_side(&seat.wall_face, &seat.cap_face, face, *on_face, *on_other)?;
                    Ok(Legs::Sphere(seat, on_wall, on_cap))
                }
            },
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            let (distance, angle) = (*distance, *angle);
            if !angle.is_finite() || angle <= tol.angular() {
                ogeom_bail!(
                    Construction,
                    "a chamfer at an angle of {angle} cuts nothing"
                );
            }
            match seat_kind(model, edge, tol)? {
                SeatKind::Straight => {
                    let seat = planar_seat(model, solid, edge, tol)?;
                    let i = seat_side(&seat, face)?;
                    // In the cross-section: from the contact on the named
                    // face, the bevel leaves at `angle` into the wedge's own
                    // side: the material on a convex edge, the open dihedral
                    // on a concave one. Where it crosses the other leg's ray
                    // is the derived distance; no crossing, no chamfer.
                    let sign = if seat.convex { 1.0 } else { -1.0 };
                    let a = seat.leg(i, tol)? * sign;
                    let b = seat.leg(1 - i, tol)? * sign;
                    let inward = -seat.normals[i] * sign;
                    let denominator = angle.sin().mul_add(b.dot(a), angle.cos() * b.dot(inward));
                    if denominator <= tol.angular() {
                        ogeom_bail!(
                            Construction,
                            "the bevel at that angle never meets the edge's other face"
                        );
                    }
                    let derived = distance * angle.sin() / denominator;
                    let mut distances = [0.0; 2];
                    distances[i] = distance;
                    distances[1 - i] = derived;
                    Ok(Legs::Straight(seat, distances))
                }
                SeatKind::Marched => Ok(Legs::Marched),
                SeatKind::Rim(rim) => {
                    // In the rim's section the legs make the triangle with
                    // the bevel: `angle` at the named face's contact and the
                    // corner's own angle at the rim. A cylinder's seat is
                    // square, the derived distance the plain tangent; past
                    // the two angles' sum of a half turn the bevel walks
                    // away from the other face instead of toward it.
                    let seat = match rim_seat(model, solid, edge, &rim, tol)? {
                        RimSeat::Revolved(seat) => seat,
                        RimSeat::Sphere(seat) => {
                            // On a sphere the section's wall is an arc: the
                            // other contact is where the line meets it.
                            let on_cap_face = if seat.cap_face.is_same(face) {
                                true
                            } else if seat.wall_face.is_same(face) {
                                false
                            } else {
                                ogeom_bail!(
                                    Construction,
                                    "the named face does not meet the edge being chamfered"
                                )
                            };
                            let (on_wall, on_cap) =
                                seat.angle_legs(on_cap_face, distance, angle, tol)?;
                            return Ok(Legs::Sphere(seat, on_wall, on_cap));
                        }
                    };
                    let corner = rim_corner(&seat);
                    let past = if seat.slope == 0.0 {
                        angle >= core::f64::consts::FRAC_PI_2 - tol.angular()
                    } else {
                        angle + corner >= core::f64::consts::PI - tol.angular()
                    };
                    if past {
                        ogeom_bail!(
                            Construction,
                            "the bevel at that angle never meets the edge's other face"
                        );
                    }
                    let derived = if seat.slope == 0.0 {
                        distance * angle.tan()
                    } else {
                        distance * angle.sin() / (angle + corner).sin()
                    };
                    let (on_wall, on_cap) =
                        revolved_side(&seat.wall_face, &seat.cap_face, face, distance, derived)?;
                    Ok(Legs::Rim(seat, on_wall, on_cap))
                }
            }
        }
    }
}

/// A chamfer's wedge, built and not yet applied: its faces, and whether it
/// fuses (a concave edge) or cuts (a convex one).
struct Wedge {
    faces: Vec<Shape>,
    /// The wedge already closed into a solid, where it was built as one.
    solid: Option<Shape>,
    additive: bool,
    /// The whole ring a sphere's rim chamfers from any of the arcs a
    /// boolean split it into.
    whole_ring: Option<WholeRing>,
    /// The spelling of a bevel on a marched seat, built against the solid
    /// as it stands when its turn in the request comes.
    marched: Option<Chamfer>,
}

/// A sphere rim's whole ring and the legs its bevel runs.
#[derive(Debug, Clone, Copy)]
struct WholeRing {
    centre: ogeom_math::Point,
    radius: f64,
    up: ogeom_math::Vector,
    on_wall: f64,
    on_cap: f64,
}

impl WholeRing {
    /// Whether `other` is a ring on the same circle and the same cap.
    fn same_circle(&self, other: &Self, tol: Tolerances) -> bool {
        self.centre.distance(other.centre) <= tol.confusion()
            && (self.radius - other.radius).abs() <= tol.confusion()
            && self.up.dot(other.up) >= 1.0 - tol.angular()
    }

    /// Whether `other`'s legs are this ring's.
    fn same_legs(&self, other: &Self, tol: Tolerances) -> bool {
        (self.on_wall - other.on_wall).abs() <= tol.confusion()
            && (self.on_cap - other.on_cap).abs() <= tol.confusion()
    }
}

impl Wedge {
    fn apply(
        self,
        model: &mut Model,
        solid: &Shape,
        edge: &Shape,
        tol: Tolerances,
    ) -> OgeomResult<Built> {
        let mut built = match self.solid {
            Some(wedge) => {
                let mut result = if self.additive {
                    ogeom_bool::fuse(model, solid, &wedge, tol)?
                } else {
                    ogeom_bool::cut(model, solid, &wedge, tol)?
                };
                result.history.delete(edge);
                result
            }
            None => apply_wedge(model, solid, Some(edge), &self.faces, self.additive, tol)?,
        };
        crate::support::credit_new_faces(model, solid, edge, &mut built)?;
        Ok(built)
    }
}

/// Which seat a chamfer is standing on, read from the edge's curve.
enum SeatKind {
    /// A straight edge between planes: the triangular-prism wedge.
    Straight,
    /// A circular rim: the revolved wedge with a conical bevel.
    Rim(ogeom_geom::CircleCurve),
    /// Any other crease: the marched seat, its bevel ruled.
    Marched,
}

fn seat_kind(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<SeatKind> {
    let (curve, _) = edge_curve(model, edge, tol)?;
    match curve {
        Curve::Line(_) => Ok(SeatKind::Straight),
        Curve::Circle(c) => Ok(SeatKind::Rim(c)),
        _ => Ok(SeatKind::Marched),
    }
}

/// Assign a named face's distance to the wall or the cap.
///
/// On the wall the distance runs axially down from the rim; on the cap it
/// runs radially in from it.
fn revolved_side(
    wall_face: &Shape,
    cap_face: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
) -> OgeomResult<(f64, f64)> {
    if wall_face.is_same(face) {
        Ok((on_face, on_other))
    } else if cap_face.is_same(face) {
        Ok((on_other, on_face))
    } else {
        ogeom_bail!(
            Construction,
            "the named face does not meet the edge being chamfered"
        )
    }
}

/// The bevel's cone (a cylinder where the two radii agree) through a rim
/// chamfer's two tangency rings: radius `cap_rho` at the cap's level and
/// the wall's `wall_rho` a `depth` down the axis, charted over both rows
/// with a margin clear of the apex.
pub(crate) fn rim_bevel_surface(
    seat: &RevolvedSeat,
    depth: f64,
    cap_rho: f64,
    wall_rho: f64,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let slope = (cap_rho - wall_rho) / (seat.tau * depth);
    // The domain covers the band's two rows (the cap ring at zero and
    // the wall ring a depth away) with a margin that stays clear of the
    // apex, where the surface degenerates.
    let rows = (
        0.0_f64.min(-seat.tau * depth),
        0.0_f64.max(-seat.tau * depth),
    );
    let frame = seat.frame_at(seat.centre, tol)?;
    let surface: SurfaceGeometry = if slope.atan().abs() <= tol.angular() {
        let pad = 0.1 * depth;
        ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(frame, cap_rho, tol)?,
            (rows.0 - pad, rows.1 + pad),
        )?
        .into()
    } else {
        // Beside a cone's wall the bevel's apex can come near its
        // rows; the margin then stops halfway to it.
        let pad = if seat.slope == 0.0 {
            0.1 * depth
        } else {
            let apex = -cap_rho / slope;
            let clear = (apex - rows.0).abs().min((apex - rows.1).abs());
            (0.1 * depth).min(0.5 * clear)
        };
        let cone = Cone::new(frame, cap_rho, slope.atan(), tol)?;
        ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into()
    };
    Ok(surface)
}

/// The revolved wedge with the quarter-tube exchanged for a slant: legs
/// `on_wall` down the wall along its ruling and `on_cap` radially along the
/// cap, and the cone between the two tangency rings as the bevel.
fn revolved_bevel(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    seat: &RevolvedSeat,
    on_wall: f64,
    on_cap: f64,
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    for distance in [on_wall, on_cap] {
        if !distance.is_finite() || distance <= tol.confusion() {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    let cap_rho = seat.sigma.mul_add(-(seat.tau * on_cap), seat.radius);
    if cap_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_cap} along the cap swallows the axis of a rim \
             of radius {}",
            seat.radius
        );
    }
    // A cone's wall leg runs its distance along the slant, a depth down
    // the axis, to the wall's radius there.
    let depth = wall_depth(seat, on_wall);
    let wall_rho = seat.radius_at(-seat.tau * depth);
    if wall_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a chamfer of {on_wall} down the wall runs past the apex of its cone"
        );
    }
    if let Some((start, sweep)) = arc_span(model, edge, seat, tol)? {
        let legs = (depth, cap_rho, on_wall.max(on_cap));
        return sector_bevel(model, solid, edge, seat, legs, (start, sweep), tol);
    }
    let flanks = revolved_flanks(model, seat, depth, cap_rho, tol)?;

    // The bevel: the cone through both tangency rings (reference radius
    // `cap_rho` at the cap's level, the wall's radius a depth below), a
    // cylinder where the two radii agree. Unlike the fillet's quarter-tube,
    // whose away-from-the-tube normal tracks the wedge seat by seat, the
    // cone's natural normal always points away from the axis, and the
    // wedge sits on the axis side of the slant exactly when `sigma` and
    // `tau` agree.
    let bevel_band = {
        let surface = rim_bevel_surface(seat, depth, cap_rho, wall_rho, tol)?;
        let band = make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?;
        if seat.sigma * seat.tau > 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    Ok(Wedge {
        faces: vec![flanks.wall_band, flanks.annulus, bevel_band],
        solid: None,
        additive: seat.additive(),
        whole_ring: None,
        marched: None,
    })
}

/// Where a rim that is only an arc of its circle starts, and how far it
/// turns, as angles about the seat's axis from its `x_ref`, turning
/// positively about `up`. `None` for a whole circle, or an arc of one whose
/// wall runs on round the rest of the turn.
fn arc_span(
    model: &Model,
    edge: &Shape,
    seat: &RevolvedSeat,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    let (curve, range) = edge_curve(model, edge, tol)?;
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    let Curve::Circle(arc) = &curve else {
        return Ok(None);
    };
    if closed || crate::support::runs_whole_turn(model, seat, arc, range, tol)? {
        return Ok(None);
    }
    let frame = seat.frame_at(seat.centre, tol)?;
    Ok(Some(arc_window(model, edge, &frame, tol)?))
}

/// Where an arc of a rim starts and how far it turns, as angles about
/// `frame`'s `z` from its `x`, whose origin is the rim's centre.
fn arc_window(
    model: &Model,
    edge: &Shape,
    frame: &Frame,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    use ogeom_geom::Curve3d as _;
    let (curve, range) = edge_curve(model, edge, tol)?;
    let centre = frame.origin();
    let angle = |t: f64| -> OgeomResult<f64> {
        let d = curve.point_at(t, tol)? - centre;
        Ok(d.dot(frame.y().vector()).atan2(d.dot(frame.x().vector())))
    };
    let turn = core::f64::consts::TAU;
    let (a, b, m) = (
        angle(range.0)?,
        angle(range.1)?,
        angle(0.5 * (range.0 + range.1))?,
    );
    // The arc runs from `a` to `b` one way round or the other: the way
    // that passes through its middle.
    let forward = (b - a).rem_euclid(turn);
    let start = if (m - a).rem_euclid(turn) <= forward {
        a
    } else {
        b
    };
    let sweep = if start == a { forward } else { turn - forward };
    Ok((start, sweep))
}

/// The wedge of a chamfer along an arc of a rim: the bevel's cross-section,
/// the triangle between the rim and the two tangency rings, turned about the
/// rim's axis through the arc alone. Its ends are the triangles a straight
/// bevel meeting the arc at a tangent ends in, so the bevels of a chain of
/// lines and arcs abut exactly; a whole ring would cut the material on the
/// far side of the axis as well.
///
/// Where the arc ends on a planar face oblique to the crease there, the
/// section in the meridian through that end would stop short of the face
/// and leave the corner between them: the bevel runs on past the end until
/// its section has left the solid across that face, within four of its
/// longer leg, as a straight band's run-out is read, and is cut back to the
/// face's plane.
fn sector_bevel(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    seat: &RevolvedSeat,
    (depth, cap_rho, leg): (f64, f64, f64),
    (start, sweep): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    let frame = seat.frame_at(seat.centre, tol)?;
    let radial_at =
        |theta: f64| frame.x().vector() * theta.cos() + frame.y().vector() * theta.sin();
    let section_at = |theta: f64| -> [ogeom_math::Point; 3] {
        let radial = radial_at(theta);
        let rim = seat.centre + radial * seat.radius;
        let wall = if seat.slope == 0.0 {
            rim - seat.up * (seat.tau * depth)
        } else {
            let height = -seat.tau * depth;
            seat.centre + seat.up * height + radial * seat.radius_at(height)
        };
        [rim, wall, seat.centre + radial * cap_rho]
    };
    let mut window = [start, start + sweep];
    let mut trims = Vec::new();
    let hosts = [&seat.cap_face, &seat.wall_face];
    for (slot, sign) in [(0_usize, -1.0), (1, 1.0)] {
        let end = window[slot];
        let at = seat.centre + radial_at(end) * seat.radius;
        let ahead = seat.up.cross(radial_at(end)) * sign;
        if !crate::marched::crease_terminates_at(model, solid, edge, hosts, at, tol)? {
            continue;
        }
        let Some(face) = crate::support::oblique_end_face(model, solid, at, ahead, tol)? else {
            continue;
        };
        let step = (leg / seat.radius) / 8.0;
        let reach = 4.0 * leg / seat.radius;
        let section = |t: f64| Ok(section_at(end + sign * t));
        if let Some(past) =
            crate::support::run_out_past(model, solid, face, step, reach, &section, tol)?
        {
            window[slot] = end + sign * past;
            trims.push(face);
        }
    }
    // A bevel that would close the turn running on keeps its ends.
    if window[1] - window[0] >= core::f64::consts::TAU - tol.angular() {
        window = [start, start + sweep];
        trims.clear();
    }
    let [rim, wall, cap] = section_at(window[0]);
    let across = seat.up.cross(radial_at(window[0]));
    let profile = crate::support::planar_face(model, &[rim, wall, cap], across, tol)?;
    let axis = ogeom_math::Axis::new(seat.centre, frame.z());
    let turn = window[1] - window[0];
    let mut wedge = ogeom_algo::make_revolution(model, &profile, axis, turn, tol)?.shape;
    let span = (seat.radius + leg) * 4.0;
    for plane in trims {
        wedge = crate::support::trimmed_wedge(model, &wedge, plane, span, tol)?;
    }
    Ok(Wedge {
        faces: Vec::new(),
        solid: Some(wedge),
        additive: seat.additive(),
        whole_ring: None,
        marched: None,
    })
}

/// Which side of the seat a named face is, by identity.
fn seat_side(seat: &Seat, face: &Shape) -> OgeomResult<usize> {
    if seat.faces[0].is_same(face) {
        Ok(0)
    } else if seat.faces[1].is_same(face) {
        Ok(1)
    } else {
        ogeom_bail!(
            Construction,
            "the named face does not meet the edge being chamfered"
        )
    }
}

/// The one construction under all three spellings: the wedge with legs
/// `distances[i]` along face `i`, subtracted.
fn bevel(
    model: &mut Model,
    seat: &Seat,
    distances: [f64; 2],
    tol: Tolerances,
) -> OgeomResult<Wedge> {
    for distance in distances {
        if !distance.is_finite() || distance <= tol.confusion() {
            ogeom_bail!(Construction, "a chamfer of {distance} cuts nothing");
        }
    }
    // On a concave edge every leg mirrors: the wedge sits in the open
    // dihedral, its legs walk the faces' planes into it, and its strips face
    // the material they will melt against with *opposed* orientation, which
    // is exactly what a fuse cancels.
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(0, tol)? * sign;
    let b = seat.leg(1, tol)? * sign;
    // A setback runs back across each face; past the face's far side it
    // would cut through the face and on into whatever lies beyond it.
    if seat.convex {
        for (i, leg) in [a, b].into_iter().enumerate() {
            let reach = crate::support::face_reach(model, &seat.faces[i], seat.start, leg, tol)?;
            if distances[i] > reach + tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "a chamfer of {} on the edge from {:?} to {:?} runs past its face, \
                     which reaches {reach} back from the edge",
                    distances[i],
                    seat.start,
                    seat.end
                );
            }
        }
    }

    let travel = seat.end - seat.start;
    let apex0 = seat.start;
    let apex1 = seat.end;
    let a0 = apex0 + a * distances[0];
    let b0 = apex0 + b * distances[1];
    let a1 = a0 + travel;
    let b1 = b0 + travel;

    // The bevel's outward normal: perpendicular to the cut line and the edge,
    // pointing from the apex toward the cut. For equal distances this is the
    // leg bisector exactly.
    let bevel_out = {
        let across = b0 - a0;
        let mut n = seat.along.cross(across);
        let m = n.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "the chamfer's cut line has no direction");
        }
        n /= m;
        if n.dot(a0 - apex0) < 0.0 {
            n = -n;
        }
        n
    };

    // The wedge: a triangular prism whose apex line is the edge and whose
    // legs run the distances along each face. Built from five explicit planar
    // faces rather than swept, because a sweep's walls are extrusion
    // surfaces even when they are geometrically planes, and the boolean's
    // same-domain resolution (which is what makes the coplanar legs melt
    // into the solid's own faces) recognises coincidence between *planes*.
    let faces = [
        planar_face(model, &[apex0, a0, b0], -seat.along, tol)?,
        planar_face(model, &[apex1, a1, b1], seat.along, tol)?,
        planar_face(model, &[apex0, a0, a1, apex1], seat.normals[0] * sign, tol)?,
        planar_face(model, &[apex0, b0, b1, apex1], seat.normals[1] * sign, tol)?,
        planar_face(model, &[a0, b0, b1, a1], bevel_out, tol)?,
    ];
    Ok(Wedge {
        faces: faces.to_vec(),
        solid: None,
        additive: !seat.convex,
        whole_ring: None,
        marched: None,
    })
}

/// The faces of `solid` holding `edge`.
fn faces_holding(model: &Model, solid: &Shape, edge: &Shape) -> Result<Vec<Shape>, Blamed> {
    let mut found = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Face)? {
        if ogeom_topo::explore_unique(model, &face, ogeom_topo::ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(edge))
        {
            found.push(face);
        }
    }
    Ok(found)
}

/// A bevel on a marched seat, applied to the solid the request's earlier
/// bevels (`earlier`) left (`so_far`, with their history): the edge found
/// again there between the faces `given` it lay between as asked, a face
/// the spelling names followed to what those bevels made of it, each piece
/// of the edge bevelled in turn.
///
/// A bevel runs on through a vertex where the crease carries on (an open
/// marched bevel runs out as the marched fillet does), so an edge an
/// earlier bevel of the same spelling meeting it at a vertex has taken
/// away is bevelled already; one an earlier bevel of another spelling has
/// taken is refused.
fn marched_step(
    model: &mut Model,
    so_far: &Built,
    edge: &Shape,
    spec: &Chamfer,
    given: &[Shape],
    earlier: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let solid = &so_far.shape;
    let faces_now = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Face)?;
    let now = |named: &Shape| -> Shape {
        if faces_now.iter().any(|f| f.is_same(named)) {
            return named.clone();
        }
        so_far
            .history
            .modified(named)
            .iter()
            .find(|f| faces_now.iter().any(|g| g.is_same(f)))
            .cloned()
            .unwrap_or_else(|| named.clone())
    };
    let asked = spec;
    let spec = &match spec {
        Chamfer::Symmetric(d) => Chamfer::Symmetric(*d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => Chamfer::Distances {
            face: now(face),
            on_face: *on_face,
            on_other: *on_other,
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => Chamfer::Angle {
            face: now(face),
            distance: *distance,
            angle: *angle,
        },
    };
    // An edge nothing of which is left is reported consumed.
    let pieces = match crate::fillet::refind_edges(model, solid, edge, given, tol) {
        Ok(pieces) => pieces,
        Err(ogeom_core::OgeomError::Construction(_)) => Vec::new(),
        Err(e) => return Err(e),
    };
    if pieces.is_empty() {
        let ends = |e: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
            Ok(ogeom_algo::edge_vertices(model, e)?
                .map(|(a, b)| vec![a.node(), b.node()])
                .unwrap_or_default())
        };
        let own = ends(edge)?;
        for (other, theirs) in earlier {
            if same_spelling(theirs, asked, tol) && ends(other)?.iter().any(|v| own.contains(v)) {
                let mut done = Built::from_nothing(solid.clone());
                done.history.delete(edge);
                return Ok(done);
            }
        }
        ogeom_bail!(
            Construction,
            "a bevel on a marched seat whose edge an earlier bevel of the request, \
             of another spelling, took away is not built"
        );
    }
    let mut built = Built::from_nothing(solid.clone());
    for piece in pieces {
        let mut step =
            crate::chamfer_marched::marched_bevel(model, &built.shape, &piece, spec, given, tol)?;
        crate::support::credit_new_faces(model, &built.shape, &piece, &mut step)?;
        step.history.delete(edge);
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// Whether two spellings bevel alike: the same kind, the same measures, the
/// same face named.
fn same_spelling(a: &Chamfer, b: &Chamfer, tol: Tolerances) -> bool {
    let near = |x: f64, y: f64| (x - y).abs() <= tol.confusion();
    match (a, b) {
        (Chamfer::Symmetric(x), Chamfer::Symmetric(y)) => near(*x, *y),
        (
            Chamfer::Distances {
                face: fa,
                on_face: oa,
                on_other: ra,
            },
            Chamfer::Distances {
                face: fb,
                on_face: ob,
                on_other: rb,
            },
        ) => fa.is_same(fb) && near(*oa, *ob) && near(*ra, *rb),
        (
            Chamfer::Angle {
                face: fa,
                distance: da,
                angle: aa,
            },
            Chamfer::Angle {
                face: fb,
                distance: db,
                angle: ab,
            },
        ) => fa.is_same(fb) && near(*da, *db) && (aa - ab).abs() <= tol.angular(),
        _ => false,
    }
}

/// The refusal of two pieces of one curved edge (split by seams) asked
/// with different chamfers: on planes the bevels step between them, and
/// the curved engine does not build that step yet (ADR 0001, amendment
/// "what M6 refuses by name").
fn apart_refused() -> ogeom_core::OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "two pieces of one curved edge are asked with different chamfers; the step \
         between them is built on planes only, so ask the edge once or give its \
         pieces one chamfer"
    )
}

/// The refusal of a solid whose planar topology is degenerate (`what`):
/// a chamfer needs every edge between two faces (ADR 0001, the M6
/// refusals).
fn degenerate_refused(what: &str) -> ogeom_core::OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "the solid's topology is degenerate ({what}); a chamfer needs every edge between \
         exactly two faces"
    )
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_math::{Direction, Frame, Point};
    use ogeom_topo::ShapeType;

    const T: Tolerances = Tolerances::millimetres();

    fn block(m: &mut Model, (x, y, z): (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
        let at = Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_box(m, at, size, T).unwrap().shape
    }

    /// A plate, a step and a block fused: faces split by seams, which the
    /// engine dissolves before it reads the solid. A chamfer wider than the
    /// part, which the engine does not build, leaves the model as it was:
    /// the dissolved solid's nodes are gone again, for the wedge path to
    /// read the solid as given.
    #[test]
    fn a_request_the_engine_does_not_build_leaves_the_model_as_it_was() {
        let mut m = Model::new();
        let plate = block(&mut m, (0.0, 0.0, 0.0), (80.0, 60.0, 6.0));
        let step = block(&mut m, (0.0, 30.0, 6.0), (80.0, 30.0, 10.0));
        let s = ogeom_bool::fuse(&mut m, &plate, &step, T).unwrap().shape;
        let edge = ogeom_topo::explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let nodes = m.node_count();
        let built = engine_chamfer(&mut m, &s, &[(edge, Chamfer::Symmetric(100.0))], T);
        assert!(!matches!(built, Ok(Some(_))), "a chamfer of 100 built");
        assert_eq!(m.node_count(), nodes, "nodes left behind");
    }
}
