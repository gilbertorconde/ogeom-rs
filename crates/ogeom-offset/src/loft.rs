//! Solid lofts closed from loft sheets.
//!
//! The walls are the loft sheet through the sections ([`make_loft_surface`]):
//! ruled or smooth, matched by arc length across sections of different edge
//! counts, and following guide curves where they are given. The sheet
//! passes through every section exactly, so its two open ends are the
//! first and last sections; each is closed by a face on its own plane,
//! bounded by the sheet's own edges there, and walls and caps are sewn
//! into one closed shell.

use ogeom_algo::{
    Built, attach_pcurve, edge_vertices, find_plane, is_wire_closed, make_face, make_solid,
    make_wire, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{PlaneSurface, SurfaceGeometry};
use ogeom_math::Plane;
use ogeom_topo::{Filter, Location, Model, NodeData, Shape, ShapeType, TShapeId, explore};

use crate::sheet::make_loft_surface;

/// Loft a solid through closed sections, in order, its ends capped by the
/// first and last sections' planes.
///
/// The walls are [`make_loft_surface`] through the sections, open in the
/// loft direction: with `ruled` straight between each pair of neighbouring
/// sections, otherwise the smooth skin through them all, and with `guides`
/// the Gordon surface that also follows each guide. Sections of different
/// edge counts, lines and arcs and splines mixed, are first split to match
/// by arc length from each section's own start, every piece keeping its
/// exact form, and the pieces pair point to point at one fraction of their
/// exact parameters. A ruled solid between a square and a circle is so
/// exact: its walls are the degree-one surfaces between the square's sides
/// and the quarter arcs they pair with, and its caps the two sections'
/// planes. Aligning the sections' starts is the caller's authorship.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// section is not a closed edge or wire, if the first or last section is
/// not planar, if the walls and caps do not close into one shell, and as
/// [`make_loft_surface`] refuses the sections or the guides (fewer than two
/// sections; with guides, a ruled loft or sections of more than one edge).
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) as
/// [`make_loft_surface`].
pub fn make_loft_solid(
    model: &mut Model,
    sections: &[Shape],
    guides: &[Shape],
    ruled: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for (k, section) in sections.iter().enumerate() {
        if !is_closed_section(model, section, tol)? {
            ogeom_bail!(
                Construction,
                "a solid loft runs through closed sections; section {k} is not a closed edge or \
                 wire"
            );
        }
    }
    let sheet = make_loft_surface(model, sections, false, guides, ruled, tol)?.shape;
    let mut faces = explore(model, &sheet, Filter::OfType(ShapeType::Face))?;
    let ends = free_loops(model, &faces)?;
    if ends.len() != 2 {
        ogeom_bail!(
            Construction,
            "the loft's walls are open along {} loops, not at its two ends",
            ends.len()
        );
    }
    for ring in &ends {
        let wire = make_wire(model, ring, tol)?.shape;
        let Some(plane) = find_plane(model, &wire, tol)? else {
            ogeom_bail!(
                Construction,
                "a solid loft is capped by its end sections' planes; an end section is not planar"
            );
        };
        faces.push(cap(model, plane, ring, wire, tol)?);
    }
    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !sewn.free_edges.is_empty() {
        ogeom_bail!(
            Construction,
            "the loft's walls and caps did not close into one shell"
        );
    }
    let mut built = make_solid(model, &sewn.shells)?;
    for section in sections.iter().chain(guides) {
        built.history.generate(section, built.shape.clone());
    }
    Ok(built)
}

/// The face on `plane` bounded by `wire`, the ring of `edges`, each edge
/// given its exact trim in the plane's chart.
fn cap(
    model: &mut Model,
    plane: Plane,
    edges: &[Shape],
    wire: Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let surface: SurfaceGeometry = PlaneSurface::new(plane).into();
    let face = make_face(model, surface.clone(), &[wire], tol)?.shape;
    let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "the cap holds no face data");
    };
    let surface_id = data.surface;
    for edge in edges {
        let (curve, range) = crate::sweep::spine_curve_of(model, edge)?;
        let Some(pcurve) = ogeom_intersect::exact_pcurve_over(&curve, range, &surface, tol) else {
            ogeom_bail!(
                Construction,
                "an end section's edge has no exact trim on the section's plane"
            );
        };
        attach_pcurve(model, edge, pcurve, surface_id, Location::identity(), range)?;
    }
    Ok(face)
}

/// Whether `section` is an edge or a wire that closes on itself.
fn is_closed_section(model: &Model, section: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    match model.kind_of(section)? {
        ShapeType::Wire => is_wire_closed(model, section, tol),
        ShapeType::Edge => match edge_vertices(model, section)? {
            Some((start, end)) => {
                Ok(start.is_same(&end) || model.same_position(&start, &end, tol)?)
            }
            None => Ok(false),
        },
        _ => Ok(false),
    }
}

/// The edges walked by an odd number of the faces' rings, chained end to
/// end into loops. An edge walked twice by one face, its seam, is no
/// border.
fn free_loops(model: &Model, faces: &[Shape]) -> OgeomResult<Vec<Vec<Shape>>> {
    let mut uses: Vec<(TShapeId, usize, Shape)> = Vec::new();
    for face in faces {
        for wire in model.children_of(face)? {
            for edge in model.children_of(&wire)? {
                if model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
                {
                    continue;
                }
                match uses.iter_mut().find(|(id, _, _)| *id == edge.node()) {
                    Some(entry) => entry.1 += 1,
                    None => uses.push((edge.node(), 1, edge)),
                }
            }
        }
    }
    let mut free: Vec<(Shape, TShapeId, TShapeId)> = Vec::new();
    for (_, count, edge) in uses {
        if count % 2 == 0 {
            continue;
        }
        let Some((start, end)) = edge_vertices(model, &edge)? else {
            ogeom_bail!(Construction, "a loft wall's border edge has no vertices");
        };
        free.push((edge, start.node(), end.node()));
    }
    let mut loops = Vec::new();
    while let Some((edge, first, mut at)) = free.pop() {
        let mut ring = vec![edge];
        while at != first {
            let Some(next) = free.iter().position(|(_, a, b)| *a == at || *b == at) else {
                ogeom_bail!(Construction, "a loft wall's border does not close");
            };
            let (edge, a, b) = free.swap_remove(next);
            at = if a == at { b } else { a };
            ring.push(edge);
        }
        loops.push(ring);
    }
    Ok(loops)
}
