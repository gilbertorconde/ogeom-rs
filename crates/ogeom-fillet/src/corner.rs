//! Rounding a vertex: the ball-and-block tool at a planar corner.
//!
//! The ball's centre may sit anywhere a radius in from every host plane.
//! At a corner one ball touches, that region's tip is a point and the
//! rounded corner is one spherical patch; at a corner no single ball
//! touches, the tip is a few points joined by short ridges, and the
//! rounded corner is a sphere at each and a cylinder along each ridge:
//! the exact envelope of the rolling ball, where the setback family fits
//! a plate through the bands' ends instead.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType};

/// The old corner tool at `vertex` where a curved face meets it; `None`
/// at a corner of planes, which one ball or the ridge balls round and the
/// engine builds.
pub(crate) fn round_vertex_with(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    #[cfg(feature = "old-path-log")]
    let _old_path_entry = crate::fillet::to_current_path("vertex:round_vertex_with");
    let Corner::Curved(corner) = corner_at(model, solid, vertex, radius, tol)? else {
        return Ok(None);
    };
    let mut built = crate::corner_curved::curved_corner(model, solid, vertex, corner, radius, tol)?;
    // A vertex's solid image supports the corner tool's callers, but it does
    // not identify the patch for a client naming faces. Record the new faces
    // separately, before later bands split or consume them. This uses boolean
    // ancestry, not the patch's surface kind or distance from the vertex.
    for face in corner_faces(model, solid, &built)? {
        built.history.generate(vertex, face);
    }
    Ok(Some(built))
}

/// Faces introduced by a corner, excluding every image of an input face.
fn corner_faces(model: &Model, solid: &Shape, built: &Built) -> OgeomResult<Vec<Shape>> {
    let inherited: Vec<Shape> = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?
        .iter()
        .flat_map(|face| built.history.trace(face).to_vec())
        .collect();
    Ok(
        ogeom_topo::explore_unique(model, &built.shape, ShapeType::Face)?
            .into_iter()
            .filter(|face| !inherited.iter().any(|image| image.is_same(face)))
            .collect(),
    )
}

use crate::engine::vertex_read::{Corner, corner_at};
pub(crate) use crate::engine::vertex_read::{CornerFrame, corner_frame};

pub(crate) use crate::engine::vertex_geometry::{convex_block, planes_meet};
