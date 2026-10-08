//! The checks the build loop runs on what it built, and the bookkeeping they
//! drive: planes off their vertices, faces that do not build or cannot be
//! meshed and measured, faces overlapping or turned in, seams crossed or
//! folded, facets wound backwards, faces whose meshes leave edges open, and
//! facets lent to curved faces and given back.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::Curve;
use ogeom_math::{Point, Point2};
use ogeom_topo::{EdgeRepr, Model, Shape};

use super::astray::{OVERRUN, SLIVER_FACETS, free_edges, own_meshes};
use super::caches::AreaCache;
use super::planner::REACH;
use super::seams::{distance_to_line, distance_to_segment, segments_cross};
use super::weld::{Adjacency, from_to};
use super::{Carrier, Fan, Groups};
use crate::recognize::worst_deviation;

/// The planar regions some of whose vertices stand off their plane by more
/// than twice the coplanar distance. Every plane the regions are found on
/// holds its own vertices within the distance.
pub(super) fn planes_off_their_vertices(
    points: &[Point],
    triangles: &[[u32; 3]],
    groups: &Groups,
    flat: f64,
) -> Vec<usize> {
    let mut off: Vec<usize> = triangles
        .iter()
        .zip(&groups.of)
        .filter_map(|(tri, &g)| match groups.carriers.get(g) {
            Some(Carrier::Plane(plane))
                if tri
                    .iter()
                    .any(|&v| plane.distance_to(points[v as usize]) > flat * 2.0) =>
            {
                Some(g)
            }
            _ => None,
        })
        .collect();
    off.sort_unstable();
    off.dedup();
    off
}

/// The planar regions the build has given up on, by their triangles. A
/// planar face that fails, or faces into the material, is gathered again
/// from its triangles the first time; once its triangles have been
/// gathered again, it is built a face to each triangle; a face of one
/// triangle so built is left as it is. A curved region can always be
/// faceted.
#[derive(Default)]
pub(super) struct Regrouped {
    once: ogeom_core::FastSet<usize>,
}

impl Regrouped {
    pub(super) fn can_withdraw(&self, groups: &Groups, g: usize) -> bool {
        match groups.carriers.get(g) {
            Some(Carrier::Curved(_)) => true,
            Some(Carrier::Plane(_)) => {
                let mut members = groups.of.iter().enumerate().filter(|&(_, &o)| o == g);
                let Some((first, _)) = members.next() else {
                    return false;
                };
                members.next().is_some() || !self.once.contains(&first)
            }
            _ => false,
        }
    }

    /// Record planar region `g` as withdrawn; the triangles to build a
    /// face each, where it has been gathered again already.
    pub(super) fn withdraw_plane(&mut self, groups: &Groups, g: usize) -> Vec<usize> {
        let members: Vec<usize> = (0..groups.of.len())
            .filter(|&t| groups.of[t] == g)
            .collect();
        if members.iter().any(|t| self.once.contains(t)) {
            members
        } else {
            self.once.extend(members.iter().copied());
            Vec::new()
        }
    }
}

/// The groups to withdraw for faces that could not be built: a curved one
/// itself; for a fan, the curved faces its seams lie on, which leaves it a
/// flat facet; for another planar face, the curved faces beside it, whose
/// seams it could not take, or where there are none, the face itself.
pub(super) fn unbuilt_culprits(
    unbuilt: &[usize],
    groups: &Groups,
    adjacency: &Adjacency,
    fans: &FastMap<usize, Fan>,
) -> Vec<usize> {
    let mut out = std::collections::BTreeSet::new();
    for &g in unbuilt {
        if !matches!(groups.carriers.get(g), Some(Carrier::Plane(_))) {
            out.insert(g);
        } else if let Some(fan) = fans.get(&g) {
            out.extend(fan.seams().map(|(_, c)| c));
        } else {
            let beside: Vec<usize> = adjacency
                .twin
                .iter()
                .enumerate()
                .filter_map(|(h, twin)| {
                    let o = groups.of[(*twin)? / 3];
                    (groups.of[h / 3] == g
                        && matches!(groups.carriers.get(o), Some(Carrier::Curved(_))))
                    .then_some(o)
                })
                .collect();
            if beside.is_empty() {
                out.insert(g);
            } else {
                out.extend(beside);
            }
        }
    }
    out.into_iter().collect()
}

/// The groups to withdraw when a step after the build fails on the built
/// faces: those whose face cannot be meshed on its own, or whose area
/// cannot be measured, and where none is found, every recognized region
/// (with nothing built, every one).
pub(super) fn failing_faces(
    model: &Model,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> Vec<usize> {
    let faces: Vec<(usize, &Shape)> = built
        .iter()
        .enumerate()
        .filter_map(|(g, b)| b.as_ref().map(|f| (g, f)))
        .collect();
    let fine = ogeom_mesh::Deflection::with_chord(ogeom_mesh::Deflection::default().chord * 1e-2);
    let failing = ogeom_core::parallel::map_ordered(&faces, |_, &(_, face)| {
        ogeom_mesh::triangulate_face(model, face, ogeom_mesh::Deflection::default(), tol).is_err()
            || fine
                .as_ref()
                .is_ok_and(|&fine| crate::surface_properties(model, face, fine, tol).is_err())
    });
    let found: Vec<usize> = faces
        .iter()
        .zip(failing)
        .filter(|(_, failing)| *failing)
        .map(|(&(g, _), _)| g)
        .collect();
    if !found.is_empty() {
        return found;
    }
    (0..groups.carriers.len())
        .filter(|&g| {
            matches!(groups.carriers[g], Carrier::Curved(_))
                && (built.is_empty() || built.get(g).is_some_and(Option::is_some))
        })
        .collect()
}

/// The recognized faces that overlap a face beside them: a fitted face
/// running past the facet it should end on encloses a sliver outside that
/// facet, and the facet then faces into material on its outer side while
/// every check on either face alone passes. Every face is probed as the
/// validity check probes orientation; a recognized face found turned in is
/// a culprit, and for any other face so found the recognized faces beside
/// it, or whose box it meets, are. Returned with the planar faces so found
/// that share an edge with a recognized one, those with no recognized face
/// beside them or across their box, and every recognized face whose box
/// meets a planar face so found.
#[allow(clippy::type_complexity, reason = "four lists of groups")]
pub(super) fn overlapping_faces(
    model: &Model,
    shape: &Shape,
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<(Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>)> {
    let curved = |g: usize| {
        matches!(groups.carriers.get(g), Some(Carrier::Curved(_)))
            && built.get(g).is_some_and(Option::is_some)
    };
    let mut beside: FastMap<usize, std::collections::BTreeSet<usize>> = FastMap::default();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(t) = *twin else {
            continue;
        };
        let (mine, theirs) = (groups.of[h / 3], groups.of[t / 3]);
        if mine != theirs && curved(theirs) {
            beside.entry(mine).or_default().insert(theirs);
        }
    }
    // Each built face's box; a face a fitted face runs through need not
    // share an edge with it.
    let mut boxes: FastMap<usize, ogeom_math::Aabb> = FastMap::default();
    for (g, face) in built.iter().enumerate() {
        if let Some(face) = face {
            boxes.insert(
                g,
                crate::shape_bounds(model, face, tol)?.expanded(tol.confusion()),
            );
        }
    }
    let curved_boxes: Vec<(usize, ogeom_math::Aabb)> = boxes
        .iter()
        .filter(|(g, _)| curved(**g))
        .map(|(g, b)| (*g, *b))
        .collect();
    let crossing = |g: usize| -> Vec<usize> {
        boxes.get(&g).map_or_else(Vec::new, |own| {
            curved_boxes
                .iter()
                .filter(|(c, b)| *c != g && b.intersects(own))
                .map(|(c, _)| *c)
                .collect()
        })
    };
    let of_face: FastMap<ogeom_topo::TShapeId, usize> = built
        .iter()
        .enumerate()
        .filter_map(|(g, face)| Some((face.as_ref()?.node(), g)))
        .collect();
    let mut culprits = std::collections::BTreeSet::new();
    let mut turned = Vec::new();
    let mut alone = Vec::new();
    let mut covering = std::collections::BTreeSet::new();
    for solid in ogeom_topo::explore_unique(model, shape, ogeom_topo::ShapeType::Solid)? {
        for face in crate::check::inside_out_faces(model, &solid, tol)? {
            let Some(&g) = of_face.get(&face.node()) else {
                continue;
            };
            if curved(g) {
                culprits.insert(g);
                continue;
            }
            covering.extend(crossing(g));
            if let Some(next) = beside.get(&g) {
                culprits.extend(next.iter().copied());
                turned.push(g);
            } else {
                let across = crossing(g);
                if across.is_empty() {
                    alone.push(g);
                }
                culprits.extend(across);
            }
        }
    }
    Ok((
        culprits.into_iter().collect(),
        turned,
        alone,
        covering.into_iter().collect(),
    ))
}

/// The facets that could be fans (see [`Fan`]) whose boundary as built
/// winds the wrong way round their plane: the seam with the curved face
/// bulges past the facet's far corner, and the loop encloses the crescent
/// between the seam and the facet's other two sides, turned over. Each
/// boundary is drawn through its edges' curves, a few points a curved one,
/// and its area in the plane compared with the facet's own. By triangle.
#[allow(clippy::too_many_arguments, reason = "the build's state, read")]
pub(super) fn backwards_facets(
    model: &Model,
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    fans: &FastMap<usize, Fan>,
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    use ogeom_geom::Curve3d as _;
    /// Points drawn along each curved edge.
    const DRAWN: u32 = 32;
    let mut out = Vec::new();
    for t in 0..triangles.len() {
        let g = groups.of[t];
        if fans.contains_key(&g) || groups.fan_at(t, triangles, adjacency).is_none() {
            continue;
        }
        let (Some(Carrier::Plane(plane)), Some(Some(face))) =
            (groups.carriers.get(g), built.get(g))
        else {
            continue;
        };
        let local = |p: Point| {
            let l = plane.frame().to_local(p);
            Point2::new(l.x, l.y)
        };
        let shoelace = |ring: &[Point2]| {
            ring.iter()
                .zip(ring.iter().cycle().skip(1))
                .map(|(a, b)| a.x * b.y - b.x * a.y)
                .sum::<f64>()
                * 0.5
        };
        let own = shoelace(&triangles[t].map(|v| local(points[v as usize])));
        let wires = model.ordered_children_of(face)?;
        let [wire] = &wires[..] else {
            continue;
        };
        let mut ring: Vec<Point2> = Vec::new();
        for edge in model.ordered_children_of(wire)? {
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
            else {
                continue;
            };
            let Some(curve) = model.geometry().curve(*curve) else {
                continue;
            };
            let count = if matches!(curve, Curve::Line(_)) {
                1
            } else {
                DRAWN
            };
            let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
            for k in 0..count {
                let f = f64::from(k) / f64::from(count);
                let f = if reversed { 1.0 - f } else { f };
                ring.push(local(
                    curve.point_at((range.1 - range.0).mul_add(f, range.0), tol)?,
                ));
            }
        }
        if ring.len() >= 3 && shoelace(&ring) * own <= 0.0 {
            out.push(t);
        }
    }
    Ok(out)
}

/// Whether any face of the shape's solids faces into their material.
pub(super) fn any_turned_in(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    for solid in ogeom_topo::explore_unique(model, shape, ogeom_topo::ShapeType::Solid)? {
        if !crate::check::inside_out_faces(model, &solid, tol)?.is_empty() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether a face beside a seam threaded straight has collapsed: thinner
/// than the confusion distance, its seams threaded onto one line.
pub(super) fn collapsed_beside(
    model: &Model,
    triangles: &[[u32; 3]],
    groups: &Groups,
    built: &[Option<Shape>],
    straight: &ogeom_core::FastSet<(u32, u32)>,
    areas: &AreaCache,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mut beside: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for h in 0..triangles.len() * 3 {
        let (a, b) = from_to(triangles, h);
        if straight.contains(&(a.min(b), a.max(b))) {
            beside.insert(groups.of[h / 3]);
        }
    }
    let fine = ogeom_mesh::Deflection::with_chord(ogeom_mesh::Deflection::default().chord * 1e-2)?;
    for face in beside
        .iter()
        .filter_map(|&g| built.get(g).and_then(Option::as_ref))
    {
        let area = areas.area(model, face, fine, tol)?;
        let reach = crate::shape_bounds(model, face, tol)?.diagonal();
        if area <= reach * tol.confusion() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The seams between a planar face and a curved one that bend the planar
/// face's trim back across itself: a facet thinner than the curve the two
/// surfaces meet along bulges, or a sliver between two curved faces whose
/// seams cross in its plane. Every face is meshed at the default deflection
/// on the shape's agreed edge chords, as a caller meshing the shape would.
/// A facet or two left between curved faces whose triangles cover a fifth
/// more than its exact area, or which draws a mesh edge no other face
/// draws (or that two others draw too), has its seams with curved faces
/// returned, to be threaded straight.
#[allow(clippy::too_many_arguments, reason = "the build's state, read")]
pub(super) fn crossed_seams(
    model: &Model,
    shape: &Shape,
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    (kept, areas): (&ogeom_mesh::FaceMeshCache, &AreaCache),
    tol: Tolerances,
) -> OgeomResult<Vec<(u32, u32)>> {
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    // Only a facet or two left between curved faces: a larger planar face
    // holds its own shape, and straightening its seams can turn it over.
    let mut size: FastMap<usize, usize> = FastMap::default();
    for &g in &groups.of {
        *size.entry(g).or_insert(0) += 1;
    }
    // Each such face's mesh edges against a curved face.
    let mut seams: FastMap<usize, Vec<(u32, u32)>> = FastMap::default();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(g) = *twin else {
            continue;
        };
        let (mine, theirs) = (groups.of[h / 3], groups.of[g / 3]);
        if mine == theirs
            || !curved(theirs)
            || !matches!(groups.carriers.get(mine), Some(Carrier::Plane(_)))
            || size.get(&mine).copied().unwrap_or(0) > SLIVER_FACETS
        {
            continue;
        }
        let (a, b) = from_to(triangles, h);
        seams.entry(mine).or_default().push((a.min(b), a.max(b)));
    }
    if seams.is_empty() {
        return Ok(Vec::new());
    }
    let deflection = ogeom_mesh::Deflection::default();
    let fine = ogeom_mesh::Deflection::with_chord(deflection.chord * 1e-2)?;
    let (drawn, chords) = ogeom_mesh::face_meshes_for(model, shape, deflection, Some(kept), tol)?;
    let faces: Vec<(usize, &Shape)> = built
        .iter()
        .enumerate()
        .filter_map(|(g, b)| b.as_ref().map(|f| (g, f)))
        .collect();
    let meshes = own_meshes(model, &faces, drawn, &chords, deflection, tol);
    // The faces' points welded within the confusion distance: a shared
    // edge's points come from the same chords, and a seam's two columns
    // land a rounding apart.
    let cell = tol.confusion();
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    let cell_of = |p: Point| {
        (
            (p.x / cell).round() as i64,
            (p.y / cell).round() as i64,
            (p.z / cell).round() as i64,
        )
    };
    let mut grid: FastMap<(i64, i64, i64), Vec<usize>> = FastMap::default();
    let mut welded: Vec<Point> = Vec::new();
    let mut weld = |p: Point| -> usize {
        let (x, y, z) = cell_of(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(found) = grid.get(&(x + dx, y + dy, z + dz)).and_then(|list| {
                        list.iter()
                            .copied()
                            .find(|&i| welded[i].distance(p) <= cell)
                    }) {
                        return found;
                    }
                }
            }
        }
        welded.push(p);
        grid.entry((x, y, z)).or_default().push(welded.len() - 1);
        welded.len() - 1
    };
    // The exact area of each meshed face with seams, which its drawn mesh
    // may overrun.
    let areas = ogeom_core::parallel::map_ordered(&faces, |i, &(g, face)| {
        (meshes[i].is_some() && seams.contains_key(&g)).then(|| areas.area(model, face, fine, tol))
    });
    let mut flagged: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut uses: FastMap<(usize, usize), Vec<usize>> = FastMap::default();
    for ((&(g, _), mesh), area) in faces.iter().zip(&meshes).zip(areas) {
        let Some(mesh) = mesh else {
            if seams.contains_key(&g) {
                flagged.insert(g);
            }
            continue;
        };
        let at: Vec<usize> = mesh.positions.iter().map(|&p| weld(p)).collect();
        let mut drawn = 0.0;
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            drawn += (b - a).cross(c - a).magnitude() * 0.5;
            let corners = t.map(|i| at[i as usize]);
            if corners[0] == corners[1] || corners[1] == corners[2] || corners[2] == corners[0] {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (corners[k], corners[(k + 1) % 3]);
                uses.entry((a.min(b), a.max(b))).or_default().push(g);
            }
        }
        if let Some(area) = area
            && drawn > area? * OVERRUN + tol.confusion()
        {
            flagged.insert(g);
        }
    }
    for users in uses.values() {
        if users.len() != 2 {
            flagged.extend(users.iter().copied().filter(|g| seams.contains_key(g)));
        }
    }
    let mut out: Vec<(u32, u32)> = flagged
        .iter()
        .flat_map(|g| seams[g].iter().copied())
        .collect();
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// The curved faces whose seams fold another face's boundary over itself
/// by more than those seams' tolerance.
///
/// Each built face's boundary is drawn in its surface's chart as the
/// tessellator trims it, and two of its chords that cross without sharing
/// an end are a fold: the face encloses some of its area twice, or none.
/// A fold inside the tolerance of a corner of the face is within what the
/// corner claims, and stands. One no deeper than the tolerance of the
/// edges that cross stands where one edge crosses itself; where two edges
/// of a curved face cross, the face is a sliver narrower than their
/// tolerance whose sides swap along it, and is a culprit itself. A deeper
/// one names its two edges by
/// the points they were drawn through, and the curved faces across them
/// are the culprits; where neither is curved and the folded face is, it
/// is.
///
/// Returned with the facets that could be fans (see [`Fan`]) along a fold
/// that names a curved face beside them: a facet's plane meets the curved
/// surface along a curve that can bulge past the curved face's far side
/// where that face is narrow, and the fan's seam, drawn straight in the
/// curved face's chart, does not.
pub(super) fn folded_seams(
    model: &Model,
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &Groups,
    built: &[Option<Shape>],
    tol: Tolerances,
) -> OgeomResult<(Vec<usize>, Vec<usize>)> {
    use ogeom_geom::Surface as _;
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    // The triangle of facet `f` that could be built as a fan onto curved
    // face `c`, and is not one yet.
    let fan_of = |f: usize, c: usize| -> Option<usize> {
        let mut members = groups.of.iter().enumerate().filter(|&(_, &g)| g == f);
        let (t, _) = members.next()?;
        if members.next().is_some() || groups.fans.contains(&t) {
            return None;
        }
        let fan = groups.fan_at(t, triangles, adjacency)?;
        fan.seams().any(|(_, g)| g == c).then_some(t)
    };
    let mut fans: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let deflection = ogeom_mesh::Deflection::default();
    let fine = ogeom_mesh::Deflection::with_chord(deflection.chord * 1e-2)?;
    // Each edge's faces and drawn points, by its node.
    let mut owners: FastMap<ogeom_topo::TShapeId, Vec<usize>> = FastMap::default();
    let mut drawn: FastMap<ogeom_topo::TShapeId, (Vec<Point>, f64)> = FastMap::default();
    for (g, face) in built.iter().enumerate() {
        let Some(face) = face else {
            continue;
        };
        for edge in ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)? {
            owners.entry(edge.node()).or_default().push(g);
            if let ogeom_core::collections::hash_map::Entry::Vacant(slot) = drawn.entry(edge.node())
            {
                let tolerance = model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .map_or(0.0, |d| d.tolerance.get());
                let polyline = ogeom_mesh::polyline_of_edge(model, &edge, fine, tol)?;
                slot.insert((polyline, tolerance));
            }
        }
    }
    let nearest = |at: Point, face_edges: &[ogeom_topo::TShapeId]| {
        face_edges
            .iter()
            .map(|id| {
                let d = drawn[id]
                    .0
                    .windows(2)
                    .map(|w| distance_to_segment(at, w[0], w[1]))
                    .fold(f64::INFINITY, f64::min);
                (d, *id)
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
    };
    let mut out: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for (g, face) in built.iter().enumerate() {
        let Some(face) = face else {
            continue;
        };
        let Some(surface) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        else {
            continue;
        };
        // The face's own edges; a seam bounds the face twice, up one side
        // of the chart and down the other, and is no fold.
        let Some(surface_id) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            continue;
        };
        let face_edges: Vec<ogeom_topo::TShapeId> =
            ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)?
                .iter()
                .filter(|e| {
                    !model
                        .node(e)
                        .and_then(|n| n.data().as_edge())
                        .and_then(|d| d.pcurve_for(surface_id, e.location()))
                        .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
                })
                .map(Shape::node)
                .collect();
        // The face's corners and how far each claims.
        let balls: Vec<(Point, f64)> =
            ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Vertex)?
                .iter()
                .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
                .map(|v| (v.point, v.tolerance.get()))
                .collect();
        // Drawn as the tessellator draws it first; a face whose chords
        // cross there is drawn a hundred times finer, where chords that only
        // crossed for their sag (two long arcs meeting at a sharp corner)
        // part, and what still crosses is the boundary's own fold.
        if chart_crossings(&ogeom_mesh::face_boundary(model, face, deflection, tol)?).is_empty() {
            continue;
        }
        let rings = ogeom_mesh::face_boundary(model, face, fine, tol)?;
        for [a, b, p, q] in chart_crossings(&rings) {
            let lift = |x: Point2| surface.point_at(x.x, x.y, tol);
            // Where the chords cross, along the first.
            let side = |o: Point2, x: Point2, y: Point2| {
                (x.x - o.x).mul_add(y.y - o.y, -((x.y - o.y) * (y.x - o.x)))
            };
            let (da, db) = (side(p, q, a), side(p, q, b));
            let t = da / (da - db);
            let Ok(at) = lift(Point2::new(
                (b.x - a.x).mul_add(t, a.x),
                (b.y - a.y).mul_add(t, a.y),
            )) else {
                continue;
            };
            // A crossing inside a corner's tolerance is the corner's slack:
            // its edges' ends stand off it by up to that much.
            if balls.iter().any(|&(v, r)| at.distance(v) <= r) {
                continue;
            }
            let (Ok(a), Ok(b), Ok(p), Ok(q)) = (lift(a), lift(b), lift(p), lift(q)) else {
                continue;
            };
            let depth = distance_to_line(a, p, q)
                .min(distance_to_line(b, p, q))
                .min(distance_to_line(p, a, b))
                .min(distance_to_line(q, a, b));
            // Each chord named by the edge it was drawn from; a chord no
            // edge was drawn through (a side of the chart's window, or a
            // seam) folds nothing.
            let (Some((d1, first)), Some((d2, second))) = (
                nearest(a.lerp(b, 0.5), &face_edges),
                nearest(p.lerp(q, 0.5), &face_edges),
            ) else {
                continue;
            };
            let drawn_within = fine.chord * 2.0;
            if d1 > drawn[&first].1 + drawn_within || d2 > drawn[&second].1 + drawn_within {
                continue;
            }
            let fannable: Vec<usize> = [first, second]
                .iter()
                .flat_map(|id| owners[id].iter().copied())
                .filter_map(|o| {
                    if curved(g) {
                        fan_of(o, g)
                    } else if curved(o) {
                        fan_of(g, o)
                    } else {
                        None
                    }
                })
                .collect();
            if depth <= drawn[&first].1.max(drawn[&second].1) {
                // Two of a curved face's own edges crossing within their
                // tolerance make it a sliver no wider than that tolerance,
                // its sides swapping along it: the face is drawn over
                // itself, and its neighbours fold with it.
                if first != second && curved(g) {
                    out.insert(g);
                    fans.extend(fannable);
                }
                continue;
            }
            let across: Vec<usize> = [first, second]
                .iter()
                .flat_map(|id| owners[id].iter().copied())
                .filter(|&o| o != g && curved(o))
                .collect();
            if across.is_empty() {
                if curved(g) {
                    out.insert(g);
                    fans.extend(fannable);
                }
            } else {
                out.extend(across);
                fans.extend(fannable);
            }
        }
    }
    Ok((out.into_iter().collect(), fans.into_iter().collect()))
}

/// Each pair of a boundary's chords that cross without sharing an end, as
/// the two chords' ends.
fn chart_crossings(rings: &[Vec<Point2>]) -> Vec<[Point2; 4]> {
    // Each chord, with its ring and the index of its first point, swept in
    // order of its least `u`.
    let mut chords: Vec<(usize, usize, Point2, Point2)> = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        for i in 0..ring.len() {
            chords.push((r, i, ring[i], ring[(i + 1) % ring.len()]));
        }
    }
    chords.sort_by(|a, b| a.2.x.min(a.3.x).total_cmp(&b.2.x.min(b.3.x)));
    let mut out = Vec::new();
    let mut active: Vec<usize> = Vec::new();
    for k in 0..chords.len() {
        let (r, i, a, b) = chords[k];
        let low = a.x.min(b.x);
        active.retain(|&j| chords[j].2.x.max(chords[j].3.x) >= low);
        for &j in &active {
            let (rj, ij, p, q) = chords[j];
            if r == rj {
                let n = rings[r].len();
                if (i + 1) % n == ij || (ij + 1) % n == i {
                    continue;
                }
            }
            if a.y.max(b.y) < p.y.min(q.y) || p.y.max(q.y) < a.y.min(b.y) {
                continue;
            }
            if segments_cross((a.x, a.y), (b.x, b.y), (p.x, p.y), (q.x, q.y)) {
                out.push([a, b, p, q]);
            }
        }
        active.push(k);
    }
    out
}

/// Facets given to a curved region, and what it was before, so they can
/// be handed back.
pub(super) struct Absorbed {
    /// Each facet: its own group, its plane, its triangles.
    facets: Vec<(usize, Carrier, Vec<usize>)>,
    /// The region's vertices before it took any.
    vertices: Vec<u32>,
}

/// Give each facet beside a curved region to that region.
///
/// Where curved regions meet at an angle, or where a fit stops a row short,
/// the mesh leaves planar facets of a triangle or two that no surface
/// claimed, between curved faces. Built as faces of their own they are
/// bounded by seams whose tolerance (a chord threaded through the mesh's
/// vertices strays up to a twentieth of its span) is as wide as they are,
/// so their trims fold and they mesh over themselves. A facet of at most
/// [`SLIVER_FACETS`] triangles whose corners all lie within the reach of a
/// curved neighbour's surface goes to the neighbour they lie nearest. Such
/// a facet has no vertex inside it: each of its corners ends on a seam,
/// held to the reach as any seam's points are, so the face it joins is
/// held to what it was. What each region took is returned, so a region
/// that cannot be built with its facets gives them back before it is
/// faceted itself. A facet holding a `protected` triangle stays a face.
pub(super) fn absorb_facets(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    protected: &[bool],
    flat: f64,
) -> FastMap<usize, Absorbed> {
    let reach = flat * REACH;
    let mut members: FastMap<usize, Vec<usize>> = FastMap::default();
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Plane(_))) {
            members.entry(g).or_default().push(t);
        }
    }
    let mut facets: Vec<(usize, Vec<usize>)> = members
        .into_iter()
        .filter(|(_, ts)| ts.len() <= SLIVER_FACETS)
        .filter(|(_, ts)| {
            !ts.iter()
                .any(|&t| protected.get(t).copied().unwrap_or(false))
        })
        .collect();
    facets.sort_unstable();
    let mut absorbed: FastMap<usize, Absorbed> = FastMap::default();
    for (g, ts) in facets {
        let mut beside: Vec<usize> = ts
            .iter()
            .flat_map(|&t| {
                (3 * t..3 * t + 3).filter_map(|h| adjacency.twin[h].map(|o| groups.of[o / 3]))
            })
            .filter(|&o| o != g && matches!(groups.carriers.get(o), Some(Carrier::Curved(_))))
            .collect();
        beside.sort_unstable();
        beside.dedup();
        let mut corners: Vec<u32> = ts.iter().flat_map(|&t| triangles[t]).collect();
        corners.sort_unstable();
        corners.dedup();
        let at: Vec<Point> = corners.iter().map(|&v| points[v as usize]).collect();
        let nearest = beside
            .iter()
            .filter_map(|&o| match &groups.carriers[o] {
                Carrier::Curved(c) => Some((worst_deviation(&c.shape, &at), o)),
                _ => None,
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        let Some((deviation, to)) = nearest else {
            continue;
        };
        if deviation > reach {
            continue;
        }
        let Carrier::Curved(region) = &mut groups.carriers[to] else {
            continue;
        };
        let record = absorbed.entry(to).or_insert_with(|| Absorbed {
            facets: Vec::new(),
            vertices: region.vertices.clone(),
        });
        region.vertices.extend(corners);
        region.vertices.sort_unstable();
        region.vertices.dedup();
        for &t in &ts {
            groups.of[t] = to;
        }
        let carrier = std::mem::replace(&mut groups.carriers[g], Carrier::Gone);
        record.facets.push((g, carrier, ts));
    }
    absorbed
}

/// Hand a region's facets back to themselves, as they were before it took
/// them.
pub(super) fn give_back(groups: &mut Groups, to: usize, record: Absorbed) {
    for (g, carrier, ts) in record.facets {
        for t in ts {
            groups.of[t] = g;
        }
        groups.carriers[g] = carrier;
    }
    if let Carrier::Curved(region) = &mut groups.carriers[to] {
        region.vertices = record.vertices;
    }
}

/// The curved faces whose meshes do not meet their neighbours', where the
/// solid does not mesh closed.
///
/// Each built face is meshed on its own, its edges drawn alike for every
/// face, and the meshes are joined where their points coincide. A mesh
/// edge used other than twice is a gap or an overlap: a face whose trim
/// folds within its seams' tolerance, which every check on the solid
/// passes, still meshes over itself. The curved faces among those using
/// such an edge are the culprits; where only planes use it, the curved
/// faces with an edge drawn within the reach of it are.
pub(super) fn unmatched_faces(
    model: &Model,
    shape: &Shape,
    groups: &Groups,
    built: &[Option<Shape>],
    flat: f64,
    kept: &ogeom_mesh::FaceMeshCache,
    tol: Tolerances,
) -> OgeomResult<Vec<usize>> {
    type Key = (u64, u64, u64);
    let curved = |g: usize| matches!(groups.carriers.get(g), Some(Carrier::Curved(_)));
    let deflection = ogeom_mesh::Deflection::default();
    // What the solid is drawn as decides: a crack between two faces' own
    // meshes that the drawing welds shut costs nothing.
    let (drawn, chords, each) =
        ogeom_mesh::triangulate_with_face_meshes(model, shape, deflection, Some(kept), tol)?;
    if drawn.is_closed() {
        return Ok(Vec::new());
    }
    let faces: Vec<(usize, &Shape)> = built
        .iter()
        .enumerate()
        .filter_map(|(g, b)| b.as_ref().map(|f| (g, f)))
        .collect();
    let meshes = own_meshes(model, &faces, each, &chords, deflection, tol);
    // The faces' points welded within the confusion distance: a shared
    // edge's points come from the same chords, but a seam drawn from either
    // side of a closed chart lands a rounding apart.
    let cell = tol.confusion();
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    let cell_of = |p: Point| -> Key {
        (
            (p.x / cell).round() as i64 as u64,
            (p.y / cell).round() as i64 as u64,
            (p.z / cell).round() as i64 as u64,
        )
    };
    let mut grid: FastMap<Key, Vec<usize>> = FastMap::default();
    let mut welded: Vec<Point> = Vec::new();
    let mut weld = |p: Point| -> Key {
        let c = cell_of(p);
        for dx in [0u64, 1, u64::MAX] {
            for dy in [0u64, 1, u64::MAX] {
                for dz in [0u64, 1, u64::MAX] {
                    let near = (
                        c.0.wrapping_add(dx),
                        c.1.wrapping_add(dy),
                        c.2.wrapping_add(dz),
                    );
                    if let Some(&i) = grid
                        .get(&near)
                        .and_then(|list| list.iter().find(|&&i| welded[i].distance(p) <= cell))
                    {
                        return cell_of(welded[i]);
                    }
                }
            }
        }
        welded.push(p);
        grid.entry(c).or_default().push(welded.len() - 1);
        c
    };
    let mut uses: FastMap<(Key, Key), (Vec<usize>, [Point; 2])> = FastMap::default();
    let mut out: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for (&(g, _), mesh) in faces.iter().zip(&meshes) {
        let Some(mesh) = mesh else {
            if curved(g) {
                out.insert(g);
            }
            continue;
        };
        for t in &mesh.triangles {
            for k in 0..3 {
                let (p, q) = (
                    mesh.positions[t[k] as usize],
                    mesh.positions[t[(k + 1) % 3] as usize],
                );
                let (a, b) = (weld(p), weld(q));
                if a == b {
                    continue;
                }
                let entry = uses
                    .entry((a.min(b), a.max(b)))
                    .or_insert_with(|| (Vec::new(), [p, q]));
                entry.0.push(g);
            }
        }
    }
    let free = free_edges(model, &faces, tol)?;
    let mut planar: Vec<Point> = Vec::new();
    for (users, ends) in uses.values() {
        // One face's own seam, drawn from both sides of its closed chart,
        // meets itself there, and a pole's fans with it: an even count
        // from that face alone is its own closure, not a gap.
        if users.len() == 2 || (users.len() % 2 == 0 && users.iter().all(|&g| g == users[0])) {
            continue;
        }
        // A free edge (one face, nothing across) is drawn by that face
        // alone.
        if users.len() == 1
            && free.iter().any(|(curve, reach, (lo, hi))| {
                let inside = |p: Point| {
                    (lo.x..=hi.x).contains(&p.x)
                        && (lo.y..=hi.y).contains(&p.y)
                        && (lo.z..=hi.z).contains(&p.z)
                };
                ends.iter().all(|&p| {
                    inside(p)
                        && crate::measure::project_on_curve(curve, p, 64, tol)
                            .is_ok_and(|near| near.distance <= *reach)
                })
            })
        {
            continue;
        }
        let at = ends[0].lerp(ends[1], 0.5);
        let mine: Vec<usize> = users.iter().copied().filter(|&g| curved(g)).collect();
        if mine.is_empty() {
            planar.push(at);
        } else {
            out.extend(mine);
        }
    }
    if !planar.is_empty() {
        let reach = flat * REACH;
        for &(g, face) in &faces {
            if !curved(g) || out.contains(&g) {
                continue;
            }
            'edges: for edge in
                ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Edge)?
            {
                let drawn = ogeom_mesh::polyline_of_edge(model, &edge, deflection, tol)?;
                for w in drawn.windows(2) {
                    if planar
                        .iter()
                        .any(|&p| distance_to_segment(p, w[0], w[1]) <= reach)
                    {
                        out.insert(g);
                        break 'edges;
                    }
                }
            }
        }
    }
    Ok(out.into_iter().collect())
}
