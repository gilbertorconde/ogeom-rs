//! Storing a tessellation back onto the model.
//!
//! [`triangulate`](crate::triangulate::triangulate) computes a mesh and hands it to the
//! caller. This puts one on the shape itself, as a representation alongside the
//! exact geometry (`docs/DATA_MODEL.md` §6): a polyline on each edge, a
//! triangulation on each face.
//!
//! # Why cache at all
//!
//! A viewer redraws at sixty frames a second and cannot re-solve a NURBS patch
//! each time. Data exchange writes the mesh, not the surface. Both want the
//! same answer every time they ask, which a cache guarantees and recomputation
//! does not — two calls with the same deflection can differ in their last bits,
//! and a display that flickers along a shared edge is the visible result.
//!
//! # Why the polyline keeps its parameters
//!
//! An edge's cached polyline and the boundary of a face's cached triangulation
//! have to be the same points, or the stored form has gaps where the exact
//! geometry has none. They agree because both come from the edge's 3D curve at
//! one set of parameters — so the parameters are stored, not just the points.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{
    EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, Triangulation, explore_unique,
};

use crate::discretize::{Deflection, discretize};
use crate::triangulate::triangulate_face;

/// What a tessellation pass produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tessellated {
    /// How many faces received a triangulation.
    pub faces: usize,
    /// How many edges received a polyline.
    pub edges: usize,
    /// How many triangles were produced in total.
    pub triangles: usize,
    /// Whether every face and edge met the requested deflection.
    ///
    /// `false` says the stored mesh is coarser than asked for, which a caller
    /// about to quote a tolerance needs to know.
    pub deflection_met: bool,
}

/// Tessellate every face and edge below `shape`, storing the result on the
/// model.
///
/// Replaces any tessellation already stored: a cache built to a different
/// deflection is not the one that was asked for.
///
/// # Errors
///
/// As [`triangulate_face`], plus
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn tessellate(
    model: &mut Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Tessellated> {
    deflection.validate()?;
    let mut done = Tessellated {
        faces: 0,
        edges: 0,
        triangles: 0,
        deflection_met: true,
    };

    // Edges first. A face's triangulation is built from its boundary edges, so
    // doing them in the other order would store a face mesh whose boundary the
    // edge polylines then contradict.
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        if attach_polyline(model, &edge, deflection, tol)? {
            done.edges += 1;
        }
    }

    for face in ogeom_topo::explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let mesh = triangulate_face(model, &face, deflection, tol)?;
        done.triangles += mesh.triangle_count();
        done.deflection_met &= mesh.deflection_met;

        // Each boundary edge's path through this mesh, as node indices — the
        // PolygonOnTriangulation representation. Matched while the mesh is
        // still owned, attached after it is stored.
        let mut paths: Vec<(Shape, Vec<u32>)> = Vec::new();
        let mut seen: Vec<(ogeom_topo::TShapeId, ogeom_topo::Location)> = Vec::new();
        for edge in ogeom_topo::explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let key = (edge.node(), edge.location().clone());
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);
            let points = crate::triangulate::polyline_of_edge(model, &edge, deflection, tol)?;
            if points.len() < 2 {
                continue;
            }
            if let Some(indices) = index_path(&mesh, &points, edge_reach(model, &edge, tol)) {
                paths.push((edge, indices));
            }
        }

        let id = model.geometry_mut().add_triangulation(mesh);
        for (edge, indices) in paths {
            let Some(node) = model.node_mut(&edge) else {
                continue;
            };
            let NodeData::Edge(data) = node.data_mut() else {
                continue;
            };
            data.representations.push(EdgeRepr::PolygonOnTriangulation {
                triangulation: id,
                indices,
                location: edge.location().clone(),
            });
        }

        let Some(node) = model.node_mut(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data_mut() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        data.triangulation = Some(id);
        done.faces += 1;
    }
    Ok(done)
}

/// How far a polyline point may sit from its mesh node and still be it:
/// the edge's own recorded tolerance, floored at a resolution the weld uses.
fn edge_reach(model: &Model, edge: &Shape, tol: Tolerances) -> f64 {
    let recorded = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .map_or(0.0, |d| d.tolerance.get());
    recorded.max(tol.confusion() * 1e3)
}

/// The polyline's node indices in the mesh, chosen so consecutive indices
/// are triangle edges.
///
/// A position may name several nodes — a seam's two chart columns lift to
/// the same points — so matching by position alone can jump between the
/// copies. Candidates come from position (exact bits, else within `reach`),
/// and the walk picks, at each step, a candidate adjacent in the mesh to the
/// one before it; the first point tries each of its candidates as a start.
/// `None` if no adjacency-respecting path exists.
fn index_path(mesh: &Triangulation, points: &[ogeom_math::Point], reach: f64) -> Option<Vec<u32>> {
    use std::collections::{HashMap, HashSet};
    let mut by_bits: HashMap<[u64; 3], Vec<u32>> = HashMap::new();
    for (i, p) in mesh.positions.iter().enumerate() {
        #[allow(clippy::cast_possible_truncation)]
        by_bits
            .entry([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
            .or_default()
            .push(i as u32);
    }
    let mut adjacent: HashSet<(u32, u32)> = HashSet::new();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            adjacent.insert((a.min(b), a.max(b)));
        }
    }
    let candidates = |p: &ogeom_math::Point| -> Vec<u32> {
        if let Some(exact) = by_bits.get(&[p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]) {
            return exact.clone();
        }
        let mut near: Vec<(f64, u32)> = Vec::new();
        for (i, q) in mesh.positions.iter().enumerate() {
            let d = q.distance(*p);
            if d <= reach {
                #[allow(clippy::cast_possible_truncation)]
                near.push((d, i as u32));
            }
        }
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        near.into_iter().map(|(_, i)| i).collect()
    };

    let walk = |start: u32| -> Option<Vec<u32>> {
        let mut out = vec![start];
        for p in &points[1..] {
            let previous = *out.last()?;
            let next = candidates(p)
                .into_iter()
                .find(|&c| adjacent.contains(&(previous.min(c), previous.max(c))))?;
            out.push(next);
        }
        Some(out)
    };
    candidates(points.first()?).into_iter().find_map(walk)
}

/// The triangulation stored on a face, if one has been built.
#[must_use]
pub fn triangulation_of<'a>(model: &'a Model, face: &Shape) -> Option<&'a Triangulation> {
    let NodeData::Face(data) = model.node(face)?.data() else {
        return None;
    };
    model.geometry().triangulation(data.triangulation?)
}

/// The polyline stored on an edge, if one has been built.
#[must_use]
pub fn polyline_of(model: &Model, edge: &Shape) -> Option<(Vec<ogeom_math::Point>, Vec<f64>)> {
    let NodeData::Edge(data) = model.node(edge)?.data() else {
        return None;
    };
    data.representations.iter().find_map(|repr| match repr {
        EdgeRepr::Polyline {
            points, parameters, ..
        } => Some((points.clone(), parameters.clone())),
        _ => None,
    })
}

/// Discretize an edge and store the polyline on it, replacing any earlier one.
///
/// Returns whether a polyline was stored; an edge with no 3D curve — a
/// degenerate edge at a cone's apex — has nothing to discretize.
fn attach_polyline(
    model: &mut Model,
    edge: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(false);
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let line = discretize(geometry, *range, deflection, tol)?;

    let Some(node) = model.node_mut(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data_mut() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    data.representations.retain(|repr| {
        !matches!(
            repr,
            EdgeRepr::Polyline { .. } | EdgeRepr::PolygonOnTriangulation { .. }
        )
    });
    data.add(EdgeRepr::Polyline {
        points: line.points,
        parameters: line.parameters,
        location: ogeom_topo::Location::identity(),
        deflection: deflection.chord,
    });
    Ok(true)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_algo::make_box;
    use ogeom_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn tessellating_a_box_stores_a_mesh_on_every_face_and_edge() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();

        let done = tessellate(&mut model, &built.shape, fine(), T).unwrap();
        assert_eq!(done.faces, 6);
        assert_eq!(done.edges, 12);
        assert_eq!(done.triangles, 12);
        assert!(done.deflection_met);

        for face in explore_unique(&model, &built.shape, ShapeType::Face).unwrap() {
            let mesh = triangulation_of(&model, &face).expect("face has no triangulation");
            assert_eq!(mesh.triangle_count(), 2);
        }
        for edge in explore_unique(&model, &built.shape, ShapeType::Edge).unwrap() {
            let (points, parameters) = polyline_of(&model, &edge).expect("edge has no polyline");
            assert_eq!(points.len(), parameters.len());
            assert_eq!(points.len(), 2, "a straight edge is its own polyline");
        }
    }

    #[test]
    fn the_cached_boundary_agrees_with_the_cached_faces() {
        // The point of storing the polyline's parameters. If the two caches
        // were built independently they would disagree by a hair along every
        // shared edge, and the stored form would have gaps the exact geometry
        // does not.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T).unwrap();
        tessellate(&mut model, &built.shape, fine(), T).unwrap();

        for edge in explore_unique(&model, &built.shape, ShapeType::Edge).unwrap() {
            let (points, _) = polyline_of(&model, &edge).unwrap();
            for face in
                ogeom_topo::ancestors_of(&model, &built.shape, &edge, ShapeType::Face).unwrap()
            {
                let mesh = triangulation_of(&model, &face).unwrap();
                for p in &points {
                    assert!(
                        mesh.positions.iter().any(|q| q.is_equal(*p, T)),
                        "the face's mesh has no vertex at {p:?}, which its edge's \
                         polyline passes through"
                    );
                }
            }
        }
    }

    #[test]
    fn tessellating_again_replaces_rather_than_accumulates() {
        // A cache built to a different deflection is not the one that was
        // asked for, and keeping both would leave the reader picking.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();

        tessellate(&mut model, &built.shape, Deflection::default(), T).unwrap();
        tessellate(&mut model, &built.shape, fine(), T).unwrap();

        for edge in explore_unique(&model, &built.shape, ShapeType::Edge).unwrap() {
            let NodeData::Edge(data) = model.node(&edge).unwrap().data() else {
                unreachable!()
            };
            let polylines = data
                .representations
                .iter()
                .filter(|r| matches!(r, EdgeRepr::Polyline { .. }))
                .count();
            assert_eq!(polylines, 1, "the earlier polyline was left behind");
        }
    }

    #[test]
    fn a_shape_with_no_faces_tessellates_to_nothing_rather_than_failing() {
        let mut model = Model::new();
        let vertex = model.add_point(ogeom_math::Point::ORIGIN);
        let done = tessellate(&mut model, &vertex, fine(), T).unwrap();
        assert_eq!(done.faces, 0);
        assert_eq!(done.edges, 0);
        assert_eq!(done.triangles, 0);
        assert!(done.deflection_met);
        assert!(triangulation_of(&model, &vertex).is_none());
        assert!(polyline_of(&model, &vertex).is_none());
    }

    #[test]
    fn an_unusable_deflection_is_refused_before_anything_is_stored() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let bad = Deflection {
            chord: 0.0,
            ..Deflection::default()
        };
        assert!(tessellate(&mut model, &built.shape, bad, T).is_err());

        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();
        assert!(triangulation_of(&model, &face).is_none());
    }
}
#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod polygon_on_tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_math::Frame;
    use ogeom_topo::{EdgeRepr, Filter, ShapeType, explore};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-2,
            ..Deflection::default()
        }
    }

    #[test]
    fn every_edge_walks_its_faces_triangulations_by_index() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        tessellate(&mut model, &solid.shape, fine(), T).unwrap();

        let mut checked = 0;
        for face in explore(&model, &solid.shape, Filter::OfType(ShapeType::Face)).unwrap() {
            let mesh_id = {
                let ogeom_topo::NodeData::Face(data) = model.node(&face).unwrap().data() else {
                    panic!("face data");
                };
                data.triangulation.unwrap()
            };
            let mesh = model.geometry().triangulation(mesh_id).unwrap();
            // Triangle edge set for the adjacency check.
            let mut edges_of = std::collections::HashSet::new();
            for t in &mesh.triangles {
                for i in 0..3 {
                    let (a, b) = (t[i], t[(i + 1) % 3]);
                    edges_of.insert((a.min(b), a.max(b)));
                }
            }
            for edge in explore(&model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
                let data = model.node(&edge).unwrap().data().as_edge().unwrap();
                if data.degenerate {
                    continue;
                }
                let paths: Vec<&Vec<u32>> = data
                    .representations
                    .iter()
                    .filter_map(|r| match r {
                        EdgeRepr::PolygonOnTriangulation {
                            triangulation,
                            indices,
                            ..
                        } if *triangulation == mesh_id => Some(indices),
                        _ => None,
                    })
                    .collect();
                assert!(
                    !paths.is_empty(),
                    "an edge of a tessellated face walks its triangulation"
                );
                for indices in paths {
                    assert!(indices.len() >= 2);
                    for pair in indices.windows(2) {
                        let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
                        assert!(
                            edges_of.contains(&key),
                            "consecutive indices are a triangle edge of the mesh: \
                             {:?} at {:?} and {:?}",
                            pair,
                            mesh.positions[pair[0] as usize],
                            mesh.positions[pair[1] as usize]
                        );
                    }
                    checked += 1;
                }
            }
        }
        assert!(checked >= 6, "rings, seam sides and rims all walked");
    }
}
