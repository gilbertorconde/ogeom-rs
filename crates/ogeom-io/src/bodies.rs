//! A part's shape as the exchange writers see it: solids, and the sheet
//! bodies beside them.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_topo::{Model, Shape, ShapeType};

/// The bodies under a part's shape, each with its composed placement and
/// sense, in the shape's own order.
pub(crate) struct Bodies {
    /// Every solid.
    pub solids: Vec<Shape>,
    /// Every shell no solid owns.
    pub shells: Vec<Shape>,
    /// Every face no shell owns.
    pub faces: Vec<Shape>,
    /// Every wire and edge no face owns, where the caller writes wireframe.
    pub wires: Vec<Shape>,
}

/// Split a part's shape into its solids, free shells and free faces,
/// looking through compounds and compound solids.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape holds a vertex outside every face, or a wire or edge outside every
/// face where `wireframe` is false (wireframe, which `format` is not written
/// with here), or holds nothing to write at all.
pub(crate) fn bodies_of(
    model: &Model,
    shape: &Shape,
    format: &str,
    wireframe: bool,
) -> OgeomResult<Bodies> {
    let mut bodies = Bodies {
        solids: Vec::new(),
        shells: Vec::new(),
        faces: Vec::new(),
        wires: Vec::new(),
    };
    let mut stack = vec![shape.clone()];
    while let Some(shape) = stack.pop() {
        match model.kind_of(&shape)? {
            ShapeType::Compound | ShapeType::CompSolid => {
                let children = model.children_of(&shape)?;
                stack.extend(children.into_iter().rev());
            }
            ShapeType::Solid => bodies.solids.push(shape),
            ShapeType::Shell => bodies.shells.push(shape),
            ShapeType::Face => bodies.faces.push(shape),
            ShapeType::Wire | ShapeType::Edge if wireframe => bodies.wires.push(shape),
            kind @ (ShapeType::Wire | ShapeType::Edge | ShapeType::Vertex) => ogeom_bail!(
                Construction,
                "a part's shape holds a free {}, which {format} export does not write: \
                 it writes solids, shells and faces{}",
                format!("{kind:?}").to_lowercase(),
                if wireframe {
                    " and wires, not points"
                } else {
                    ", not wireframe"
                }
            ),
        }
    }
    if bodies.solids.is_empty()
        && bodies.shells.is_empty()
        && bodies.faces.is_empty()
        && bodies.wires.is_empty()
    {
        ogeom_bail!(
            Construction,
            "a part's shape holds no solid, shell, face{} to write as {format}",
            if wireframe { " or wire" } else { "" }
        );
    }
    Ok(bodies)
}

/// Whether a writer restates `surface`'s normal: a plane, drum, cone, ball
/// or torus goes out about its axis, and a reader takes the right-handed
/// frame on that axis, so where the frame is left-handed the file's normal
/// is the surface's chart normal (`Su × Sv`) turned back. With `through`,
/// a trimmed or offset surface goes out as one over its basis, which says
/// for it (STEP's spelling); without, as a spline patch, chart and all.
/// Anything else keeps its chart.
pub(crate) fn restates(surface: &ogeom_geom::SurfaceGeometry, through: bool) -> bool {
    use ogeom_geom::SurfaceGeometry as S;
    use ogeom_math::Handedness;
    let frame = match surface {
        S::Plane(s) => s.plane().frame(),
        S::Cylinder(s) => s.cylinder().frame(),
        S::Cone(s) => s.cone().frame(),
        S::Sphere(s) => s.sphere().frame(),
        S::Torus(s) => s.torus().frame(),
        S::Trimmed(t) if through => return restates(t.basis(), through),
        S::Offset(o) if through => return restates(o.basis(), through),
        _ => return false,
    };
    frame.handedness() == Handedness::Left
}

/// Whether a face faces against its surface as the file states it. The
/// side a face faces is its placement applied to its own chart's normal,
/// which a mirroring placement turns against the placed chart's; a
/// `restated` surface (see [`restates`]) turns the file's normal against
/// the placed chart's. Where exactly one of the two holds, the face's
/// sense in the file is its orientation turned.
pub(crate) fn turned_in_file(placement: &ogeom_math::Transform, restated: bool) -> bool {
    !placement.preserves_handedness() != restated
}
