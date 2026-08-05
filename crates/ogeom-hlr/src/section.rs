//! Section views: cut the part, draw what the plane reveals.
//!
//! A section view is a boolean wearing drawing clothes. The material on the
//! plane's positive side is removed — through the general cut, against a
//! proxy box sized off the shape's own bounds, so every downstream guarantee
//! the boolean makes holds here too — and the faces the cut created *on* the
//! plane become the section outline: the closed loops a draughtsman hatches.
//! The rest of the drawing is the cut solid through the projection machinery,
//! viewed straight down the plane's normal.
//!
//! A broken-out section is the same construction with the proxy box shrunk
//! to a window: the cut reveals the interior only where the break is.

use crate::project::{Drawing, View, project};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Frame, Plane, Point, Point2};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

/// A section view: the outline on the cutting plane, and the drawing of what
/// remains behind it.
#[derive(Debug, Clone)]
pub struct SectionView {
    /// The closed loops where the plane cut material, in the plane's own
    /// `(x, y)` coordinates: one outer loop per cut face, holes after it, in
    /// the face's own wire order.
    pub outline: Vec<Vec<Point2>>,
    /// The remaining solid, projected along the plane's normal.
    pub drawing: Drawing,
    /// The solid the cut produced, for measuring or further sectioning.
    pub remainder: Shape,
}

/// Cut `solid` at `plane` and draw the section.
///
/// Material on the plane's `+z` side is removed. The outline loops are the
/// cut faces' boundaries in the plane's `(x, y)`; the drawing views the
/// remainder along the plane's normal.
///
/// # Errors
///
/// As the boolean cut and [`project()`]; and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// plane misses the solid entirely.
pub fn section(
    model: &mut Model,
    solid: &Shape,
    plane: &Plane,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<SectionView> {
    let reach = reach_of(model, solid, tol)?;
    section_with_window(
        model,
        solid,
        plane,
        (-reach, -reach),
        (reach, reach),
        deflection,
        tol,
    )
}

/// Cut only within a window of the plane: the broken-out section.
///
/// The window is a rectangle in the plane's `(x, y)`, and only material on
/// the `+z` side behind that rectangle is removed.
///
/// # Errors
///
/// As [`section()`].
pub fn broken_section(
    model: &mut Model,
    solid: &Shape,
    plane: &Plane,
    window_min: (f64, f64),
    window_max: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<SectionView> {
    section_with_window(model, solid, plane, window_min, window_max, deflection, tol)
}

/// The shared construction: a proxy box over the window, cut, outline, draw.
#[allow(clippy::too_many_arguments)]
fn section_with_window(
    model: &mut Model,
    solid: &Shape,
    plane: &Plane,
    window_min: (f64, f64),
    window_max: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<SectionView> {
    let reach = reach_of(model, solid, tol)?;
    let frame = plane.frame();
    // The proxy box: one face exactly on the plane, extending along +z; its
    // footprint is the window. Slight overshoot along z keeps the far face
    // clear of the solid.
    let corner = frame.to_world(Point::new(window_min.0, window_min.1, 0.0));
    let box_frame = Frame::new(corner, frame.z(), frame.x(), tol)?;
    let sizes = (
        window_max.0 - window_min.0,
        window_max.1 - window_min.1,
        reach,
    );
    if sizes.0 <= tol.confusion() || sizes.1 <= tol.confusion() {
        ogeom_bail!(Construction, "a section window must have area");
    }
    let proxy = ogeom_algo::make_box(model, box_frame, sizes, tol)?;
    let cut = ogeom_bool::cut(model, solid, &proxy.shape, tol)?;

    // The cut faces on the plane are the section outline. Their rings come
    // from the same boundary machinery the triangulator trusts — walked in
    // order, seams and orientations resolved — and lift from the face's own
    // chart into the section plane's coordinates through the surface.
    let mut outline = Vec::new();
    for face in explore(model, &cut.shape, Filter::OfType(ShapeType::Face))? {
        let Some(surface) = face_on_plane(model, &face, plane, tol)? else {
            continue;
        };
        let placement = face.transform(model.datums())?;
        for ring in ogeom_mesh::face_boundary(model, &face, deflection, tol)? {
            let mut loop_points = Vec::with_capacity(ring.len());
            for uv in ring {
                use ogeom_geom::Surface as _;
                let world = placement.apply(surface.point_at(uv.x, uv.y, tol)?);
                let local = frame.to_local(world);
                loop_points.push(Point2::new(local.x, local.y));
            }
            if loop_points.len() >= 3 {
                outline.push(loop_points);
            }
        }
    }
    if outline.is_empty() {
        ogeom_bail!(Construction, "the section plane misses the solid");
    }

    let view = View::looking(-frame.z().vector(), frame.y().vector(), tol)?;
    let drawing = project(model, &cut.shape, &view, deflection, tol)?;
    Ok(SectionView {
        outline,
        drawing,
        remainder: cut.shape,
    })
}

/// A margin that certainly covers the solid from any plane through it.
fn reach_of(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let bounds = ogeom_algo::shape_bounds(model, solid, tol)?;
    Ok(bounds.diagonal().max(1.0) * 2.0)
}

/// The face's surface, when the face lies on the section plane.
fn face_on_plane(
    model: &Model,
    face: &Shape,
    plane: &Plane,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::SurfaceGeometry>> {
    let Some(node) = model.node(face) else {
        return Ok(None);
    };
    let NodeData::Face(data) = node.data() else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let ogeom_geom::SurfaceGeometry::Plane(planar) = surface else {
        return Ok(None);
    };
    let placement = face.transform(model.datums())?;
    let own = planar.plane().frame();
    let origin = placement.apply(own.origin());
    let normal = placement.apply_vector(own.z().vector());
    let aligned = normal.cross(plane.frame().z().vector()).magnitude() <= 1e-9;
    let on = plane.distance_to(origin).abs() <= tol.confusion() * 1e3;
    Ok((aligned && on).then(|| surface.clone()))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_math::{Direction, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            ..Deflection::default()
        }
    }

    /// Signed shoelace area of a loop.
    fn area(points: &[Point2]) -> f64 {
        let mut sum = 0.0;
        for pair in points.windows(2) {
            sum += pair[0].x.mul_add(pair[1].y, -(pair[1].x * pair[0].y));
        }
        if let (Some(first), Some(last)) = (points.first(), points.last()) {
            sum += last.x.mul_add(first.y, -(first.x * last.y));
        }
        sum / 2.0
    }

    /// Total material area: outer loops positive, holes negative, by winding.
    fn material(outline: &[Vec<Point2>]) -> f64 {
        outline.iter().map(|l| area(l)).sum::<f64>().abs()
    }

    #[test]
    fn a_box_sections_into_its_cross_section() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        // The plane x = 5, normal +x: the half x > 5 is removed.
        let plane = Plane::new(
            Frame::new(Point::new(5.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
        );
        let view = section(&mut model, &solid.shape, &plane, fine(), T).unwrap();
        assert!(
            (material(&view.outline) - 24.0).abs() < 1e-6,
            "6 x 4 revealed, got {} from {} loops {:?}",
            material(&view.outline),
            view.outline.len(),
            view.outline.iter().map(|l| area(l)).collect::<Vec<_>>()
        );

        // Half the box remains.
        let volume = ogeom_algo::volume_properties(&model, &view.remainder, fine(), T)
            .unwrap()
            .mass;
        assert!((volume - 120.0).abs() < 0.1);
        assert!(!view.drawing.visible.is_empty());
    }

    #[test]
    fn a_bored_block_sections_through_its_hole() {
        let mut model = Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        let bore_frame =
            Frame::new(Point::new(5.0, 3.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let bore = ogeom_algo::make_cylinder(&mut model, bore_frame, 1.0, 6.0, T).unwrap();
        let part = ogeom_bool::cut(&mut model, &block.shape, &bore.shape, T).unwrap();

        // Section through the bore, half a radius off its axis: the plane
        // exactly through the axis meets the wall along its rulings, a
        // configuration the boolean still refuses (honestly), so the test
        // sections where the answer is just as exact: the slot's width is
        // the chord at that offset.
        let plane = Plane::new(
            Frame::new(Point::new(0.0, 2.5, 0.0), Direction::Y, Direction::Z, T).unwrap(),
        );
        let view = section(&mut model, &part.shape, &plane, fine(), T).unwrap();
        let chord = 2.0 * (1.0_f64 - 0.25).sqrt();
        let expected = 10.0f64.mul_add(4.0, -(chord * 4.0));
        assert!(
            (material(&view.outline) - expected).abs() < 1e-3,
            "revealed {} against {expected}",
            material(&view.outline)
        );
    }

    #[test]
    fn a_broken_section_reveals_only_its_window() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        let plane = Plane::new(
            Frame::new(Point::new(5.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
        );
        // A window over y in [1, 3] and z in [1, 3]: four square units.
        let view = broken_section(
            &mut model,
            &solid.shape,
            &plane,
            (1.0, 1.0),
            (3.0, 3.0),
            fine(),
            T,
        )
        .unwrap();
        assert!((material(&view.outline) - 4.0).abs() < 1e-6);
        // Only the window's pocket is missing.
        let volume = ogeom_algo::volume_properties(&model, &view.remainder, fine(), T)
            .unwrap()
            .mass;
        let expected = 240.0 - 4.0 * 5.0;
        assert!(
            (volume - expected).abs() < 0.1,
            "{volume} against {expected}"
        );
    }

    #[test]
    fn a_plane_that_misses_the_solid_is_refused() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
        let plane = Plane::new(
            Frame::new(Point::new(50.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
        );
        assert!(section(&mut model, &solid.shape, &plane, fine(), T).is_err());
        let _ = Vector::ZERO;
    }
}
