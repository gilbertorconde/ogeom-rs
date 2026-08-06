//! Defeaturing: putting back what a feature took away, or taking away what
//! it added.
//!
//! A recognized feature knows its own geometry, and for the features that
//! are a *volume* — a hole, a pocket, a boss — that geometry is enough to
//! rebuild the volume and undo the operation with the boolean. A hole is
//! filled by the bore it is; a pocket by the prism its floor sweeps up to
//! the material it was cut from; a boss is shaved by the same prism run the
//! other way.
//!
//! The features that are a *shape* rather than a volume — a fillet, a
//! round, a chamfer — are refused here by name. Undoing one means restoring
//! the corner it eased, which is the blend's own wedge construction run
//! backwards, and that belongs with the blends.

use ogeom_algo::{Built, make_cylinder, make_prism};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Aabb, Direction, Frame};
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

use crate::recognize::Feature;

/// How far a filling tool runs past the opening it fills, in confusion
/// tolerances.
///
/// Not zero: a tool flush with the faces it meets is a coincidence at every
/// opening at once, and the boolean does not assemble it — recorded in
/// SCOPE's deferred table with this as its reproduction.
///
/// Not small, either, and this is the part that is not obvious. A margin
/// leaves a sliver band standing past the opening, and that band's own
/// interior probes have to be *decisively* outside the part, or the exact
/// classifier finds every ray from them grazing the face they sit against,
/// exhausts its whole fan of directions, and answers On the slow way. A
/// micron of overshoot is inside the band the classifier reads as "on the
/// boundary" and costs fifty seconds on a part that takes a fifth of one
/// otherwise; ten microns is outside it and costs nothing.
///
/// So: a hundred thousand confusions, ten microns at millimetre tolerances.
/// The restored solid is larger than the original by that times the
/// openings' area — a cubic millimetre for a ten-millimetre bore — and this
/// is the number to look at if a caller's tolerance is tighter than that.
const OVERSHOOT: f64 = 1e5;

/// Remove a recognized feature from the solid it was recognized on.
///
/// A hole is filled, a pocket is filled, a boss is shaved. The result is the
/// solid as it would have been without that operation — measured, not
/// approximated: each tool is built from the feature's own surfaces.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// feature is one whose undoing is a blend construction rather than a
/// volume, or its faces do not carry the geometry the tool is built from.
pub fn remove_feature(
    model: &mut Model,
    solid: &Shape,
    feature: &Feature,
    tol: Tolerances,
) -> OgeomResult<Built> {
    match feature {
        Feature::Hole(hole) => {
            // One tool per bore face, each the cylinder that face lies on
            // over its own reach along the axis, so a counterbore fills
            // step by step and a plain bore in one.
            if hole.faces.is_empty() {
                ogeom_bail!(Construction, "a hole with no faces fills nothing");
            }
            let mut result: Option<Built> = None;
            for face in &hole.faces {
                let Some(cylinder) = bore_of(model, face, tol)? else {
                    continue;
                };
                let filled = ogeom_bool::fuse(
                    model,
                    result.as_ref().map_or(solid, |built| &built.shape),
                    &cylinder,
                    tol,
                )?;
                result = Some(filled);
            }
            match result {
                Some(built) => Ok(built),
                None => ogeom_bail!(
                    Construction,
                    "the hole carries no cylindrical bore to fill it with"
                ),
            }
        }
        Feature::Pocket(pocket) => {
            let tool = prism_over(model, &pocket.floor, &pocket.walls, false, tol)?;
            ogeom_bool::fuse(model, solid, &tool, tol)
        }
        Feature::Boss(boss) => {
            let tool = prism_over(model, &boss.top, &boss.walls, true, tol)?;
            ogeom_bool::cut(model, solid, &tool, tol)
        }
        Feature::Fillet(_) | Feature::PartialRound(_) | Feature::Chamfer(_) => ogeom_bail!(
            Construction,
            "undoing a blend or a bevel means restoring the corner it eased, \
             which is the wedge construction run backwards and belongs with \
             the blends — see docs/SCOPE.md"
        ),
    }
}

/// The solid cylinder a bore face lies on, over that face's own reach.
fn bore_of(model: &mut Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Shape>> {
    let (cylinder, placement) = {
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
            return Ok(None);
        };
        let Some(SurfaceGeometry::Cylinder(c)) = model.geometry().surface(data.surface) else {
            return Ok(None);
        };
        (c.cylinder(), face.transform(model.datums())?)
    };
    let axis = placement.apply_vector(cylinder.frame().z().vector());
    let origin = placement.apply(cylinder.frame().origin());
    // The bore's own reference direction, carried through the placement:
    // the fill has to be the *same surface* as the bore it fills, not a
    // congruent one. A fresh perpendicular would make a cylinder equal in
    // every measurable way and different in its frame, which the
    // intersector then has to discover by marching two coaxial equal-radius
    // cylinders — minutes of work to conclude what an identical frame says
    // for free.
    let reference = Direction::new(placement.apply_vector(cylinder.frame().x().vector()), tol)?;
    let bound = extent(model, std::slice::from_ref(face))?;
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        return Ok(None);
    };
    let direction = Direction::new(axis, tol)?;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for corner in Aabb::of_corners(low, high).corners() {
        let t = (corner - origin).dot(direction.vector());
        lo = lo.min(t);
        hi = hi.max(t);
    }
    if hi - lo <= tol.confusion() {
        return Ok(None);
    }
    // Past each end by the overshoot the boolean needs to see a crossing
    // rather than a coincidence — see `OVERSHOOT`. Exactly flush is the
    // coincidence, and the arrangement will not resolve it at both ends of
    // a bore at once.
    let margin = OVERSHOOT * tol.confusion();
    let base = Frame::new(
        origin + direction.vector() * (lo - margin),
        direction,
        reference,
        tol,
    )?;
    let built = make_cylinder(
        model,
        base,
        cylinder.radius() * placement.scale_factor().abs(),
        (hi - lo) + margin * 2.0,
        tol,
    )?;
    Ok(Some(built.shape))
}

/// The prism a planar face sweeps to cover the walls standing on it.
fn prism_over(
    model: &mut Model,
    planar: &Shape,
    walls: &[Shape],
    outward: bool,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some(NodeData::Face(data)) = model.node(planar).map(|n| n.data().clone()) else {
        ogeom_bail!(Construction, "expected a face");
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        ogeom_bail!(
            Construction,
            "a pocket's floor or a boss's top is planar; this one is not"
        );
    };
    let placement = planar.transform(model.datums())?;
    let mut normal = placement.apply_vector(plane.plane().normal().vector());
    if planar.orientation() == ogeom_topo::Orientation::Reversed {
        normal = -normal;
    }
    let normal = Direction::new(normal, tol)?;
    if walls.is_empty() {
        ogeom_bail!(
            Construction,
            "a floor with no walls says nothing about how deep it sits"
        );
    }
    // How far the walls reach off the face, along its own normal.
    let origin = placement.apply(plane.plane().frame().origin());
    let bound = extent(model, walls)?;
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        ogeom_bail!(Construction, "the walls bound nothing");
    };
    let mut reach = 0.0f64;
    for corner in Aabb::of_corners(low, high).corners() {
        reach = reach.max((corner - origin).dot(normal.vector()));
    }
    if reach <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the walls do not stand off the face; there is no prism to build"
        );
    }
    // As with a bore, and for the same reason.
    let margin = OVERSHOOT * tol.confusion();
    let sense = if outward { -1.0 } else { 1.0 };
    let step = normal.vector() * sense;
    let started = ogeom_algo::transformed(
        model,
        planar,
        ogeom_math::Transform::translation(-step * margin),
    )?;
    Ok(make_prism(model, &started.shape, step * (reach + margin * 2.0), tol)?.shape)
}

/// The world bound of a set of faces.
fn extent(model: &Model, faces: &[Shape]) -> OgeomResult<Aabb> {
    let mut bound = Aabb::EMPTY;
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                bound = bound.with_point(vertex.transform(model.datums())?.apply(data.point));
            }
        }
    }
    Ok(bound)
}
