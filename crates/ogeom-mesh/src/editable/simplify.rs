//! Decimating an editable mesh in place, its source vertices carried.
//!
//! [`Editable::simplify`] writes the mesh out as
//! [`Editable::to_triangulation`] does, decimates that with
//! [`simplify`](crate::simplify()), and builds the mesh again from the
//! result. The result is the free function's on the written mesh, triangle
//! for triangle, and so is what it reports. Each vertex of the result
//! takes the source vertices of every vertex merged into it, through
//! [`Simplified::result_vertex`](crate::Simplified::result_vertex), so a
//! marking of the first source still selects what it did, as the
//! [module documentation](super) defines selection.
//!
//! # What stays
//!
//! Nothing of the numbering: decimation is a compaction, as
//! [`Editable::to_triangulation`] is. Vertex `v` afterwards is the result's
//! vertex `v`, and face `f` its triangle `f`; a point the result's faces
//! meet in several fans is split as [`Editable::from`] splits it, each copy
//! with the point's source vertices. A vertex no face used is dropped with
//! its source vertices.

use ogeom_core::{OgeomResult, Tolerances};

use super::{Editable, index};
use crate::simplify::{Deviation, Target, simplify};

/// What decimating an editable mesh did: [`Simplified`](crate::Simplified)
/// without the mesh, which is now the editable one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimplifyReport {
    /// The worst error accepted, as a distance; as
    /// [`Simplified::error`](crate::Simplified::error).
    pub error: f64,
    /// How many edge collapses were made.
    pub collapsed: usize,
    /// Whether the target was reached; as
    /// [`Simplified::target_met`](crate::Simplified::target_met).
    pub target_met: bool,
    /// With [`Target::Verified`], what was measured against the mesh as it
    /// stood before; `None` otherwise.
    pub deviation: Option<Deviation>,
}

impl Editable {
    /// Decimate the mesh to `target`, as [`simplify`](crate::simplify())
    /// decimates [`Editable::to_triangulation`]'s mesh, keeping every
    /// vertex's source vertices: a vertex of the result has those of every
    /// vertex merged into it. See the [module
    /// documentation](crate::editable::simplify) for the numbering after.
    ///
    /// # Errors
    ///
    /// As [`simplify`](crate::simplify()); nothing changes when it
    /// refuses.
    pub fn simplify(&mut self, target: Target, tol: Tolerances) -> OgeomResult<SimplifyReport> {
        let written = self.to_triangulation();
        let done = simplify(&written, target, tol)?;
        let mut sources: Vec<Vec<u32>> = vec![Vec::new(); done.mesh.positions.len()];
        let kept = (0..self.positions.len()).filter(|&v| !self.gone[v]);
        for (v, r) in kept.zip(&done.result_vertex) {
            if let Some(r) = r {
                sources[*r as usize].extend_from_slice(self.source_vertices(index(v)));
            }
        }
        let mut next = Self::from(&done.mesh, 0.0)?;
        for set in &mut sources {
            set.sort_unstable();
            set.dedup();
        }
        // Built without a weld, each vertex's source is the result vertex
        // it is or was split from.
        for v in 0..next.positions.len() {
            let r = next.source[v] as usize;
            next.set_sources(index(v), sources[r].clone());
        }
        next.source_count = self.source_count;
        *self = next;
        Ok(SimplifyReport {
            error: done.error,
            collapsed: done.collapsed,
            target_met: done.target_met,
            deviation: done.deviation,
        })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::tests::{check, torus};
    use super::*;
    use crate::editable::{Cap, Side, Smooth};
    use crate::{Deflection, triangulate};
    use ogeom_algo::make_sphere;
    use ogeom_math::{Direction, Frame, Plane, Point, Vector};
    use ogeom_topo::{Model, Triangulation};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(chord: f64) -> Triangulation {
        let mut model = Model::new();
        let built = make_sphere(&mut model, Frame::WORLD, 10.0, T).unwrap();
        let deflection = Deflection {
            chord,
            ..Deflection::default()
        };
        triangulate(&model, &built.shape, deflection, T).unwrap()
    }

    const TO_200: Target = Target::Verified {
        tolerance: 0.2,
        share: 0.999,
        max_triangles: Some(200),
    };

    /// The vertex numbers there, with a face or not.
    fn present(e: &Editable) -> Vec<u32> {
        (0..e.vertex_count())
            .map(index)
            .filter(|&v| e.position(v).is_some())
            .collect()
    }

    #[test]
    fn marking_every_source_vertex_offsets_the_reduced_mesh_whole() {
        let source = sphere(0.005);
        assert!(source.is_closed() && source.triangle_count() > 10_000);
        let mut e = Editable::from(&source, 0.0).unwrap();
        let report = e.simplify(TO_200, T).unwrap();
        check(&e);

        // The free function's mesh and report, unchanged.
        let free = simplify(&source, TO_200, T).unwrap();
        let written = e.to_triangulation();
        assert_eq!(written.positions, free.mesh.positions);
        assert_eq!(written.triangles, free.mesh.triangles);
        assert!(free.mesh.triangle_count() <= 200);
        assert_eq!(report.deviation, free.deviation);
        assert_eq!(report.error.to_bits(), free.error.to_bits());
        assert_eq!(report.collapsed, free.collapsed);
        assert_eq!(report.target_met, free.target_met);
        assert_eq!(e.source_count(), source.positions.len());

        let region = e
            .region_from_sources(&vec![true; e.source_count()])
            .unwrap();
        assert!(present(&e).iter().all(|&v| region[v as usize]));
        let (mut marked, mut whole) = (e.clone(), e.clone());
        let by_marking = marked.offset(0.1, Some(&region)).unwrap();
        let by_whole = whole.offset(0.1, None).unwrap();
        assert_eq!(by_marking, by_whole);
        assert_eq!(marked.to_triangulation(), whole.to_triangulation());
        // Every vertex moved outward, the collapsed ones as well.
        let moved = marked.to_triangulation();
        for (p, q) in moved.positions.iter().zip(&written.positions) {
            assert!(p.to_vector().magnitude() > q.to_vector().magnitude());
        }
    }

    #[test]
    fn one_collapse_merges_the_two_source_sets() {
        // One collapse on a closed mesh removes two triangles, and puts
        // the vertex it keeps at the midpoint of the two it merged.
        let source = torus(12, 8);
        let mut e = Editable::from(&source, 0.0).unwrap();
        let report = e
            .simplify(Target::Triangles(source.triangle_count() - 2), T)
            .unwrap();
        check(&e);
        assert_eq!(report.collapsed, 1);
        assert_eq!(e.vertex_count(), source.positions.len() - 1);

        let mut pair = None;
        for v in present(&e) {
            let at = e.position(v).unwrap();
            match *e.source_vertices(v) {
                [s] => {
                    assert_eq!(e.source_vertex(v), Some(s));
                    assert_eq!(at, source.positions[s as usize]);
                }
                [a, b] => {
                    assert!(pair.is_none(), "a second merged vertex {v}");
                    assert_eq!(e.source_vertex(v), None);
                    let (pa, pb) = (source.positions[a as usize], source.positions[b as usize]);
                    assert_eq!(
                        at,
                        Point::from_vector((pa.to_vector() + pb.to_vector()) * 0.5)
                    );
                    pair = Some((v, a, b));
                }
                ref other => panic!("vertex {v} from {other:?}"),
            }
        }
        let (v, a, b) = pair.unwrap();
        // Every source vertex is behind exactly one vertex.
        let mut all: Vec<u32> = present(&e)
            .iter()
            .flat_map(|&w| e.source_vertices(w).to_vec())
            .collect();
        all.sort_unstable();
        assert_eq!(all, (0..index(source.positions.len())).collect::<Vec<_>>());

        // The merged vertex is selected only with both its sources marked;
        // the others with theirs, whatever their numbers now.
        for marks in [vec![a], vec![b], vec![a, b]] {
            let mut marked = vec![false; source.positions.len()];
            for &s in &marks {
                marked[s as usize] = true;
            }
            marked[0] = true;
            let region = e.region_from_sources(&marked).unwrap();
            assert_eq!(region[v as usize], marks.len() == 2);
            for w in present(&e).into_iter().filter(|&w| w != v) {
                let s = e.source_vertex(w).unwrap();
                assert_eq!(region[w as usize], marked[s as usize], "vertex {w}");
            }
        }
        assert!(e.region_from_sources(&[true; 3]).is_err());
    }

    #[test]
    fn a_partial_marking_follows_the_recorded_sets() {
        let source = sphere(0.005);
        let mut e = Editable::from(&source, 0.0).unwrap();
        e.simplify(TO_200, T).unwrap();
        let marked: Vec<bool> = source.positions.iter().map(|p| p.z > 0.0).collect();
        let region = e.region_from_sources(&marked).unwrap();

        let mut seen = vec![false; source.positions.len()];
        let mut mixed = 0;
        for v in present(&e) {
            let sources = e.source_vertices(v);
            assert!(!sources.is_empty() && sources.windows(2).all(|w| w[0] < w[1]));
            for &s in sources {
                assert!(!seen[s as usize], "source vertex {s} behind two vertices");
                seen[s as usize] = true;
            }
            if let [s] = *sources {
                // A vertex nothing merged into has not moved.
                assert_eq!(e.position(v), Some(source.positions[s as usize]));
            }
            let all = sources.iter().all(|&s| marked[s as usize]);
            let any = sources.iter().any(|&s| marked[s as usize]);
            mixed += usize::from(any && !all);
            assert_eq!(region[v as usize], all, "vertex {v} from {sources:?}");
        }
        assert!(seen.iter().all(|&s| s), "a source vertex behind none");
        // Not a marking every set falls on one side of.
        assert!(mixed > 0);
        assert!(region.iter().any(|&r| r) && present(&e).iter().any(|&v| !region[v as usize]));
    }

    /// Each source vertex, and how many vertices name it.
    fn uses(e: &Editable) -> Vec<usize> {
        let mut count = vec![0; e.source_count()];
        for v in present(e) {
            for &s in e.source_vertices(v) {
                count[s as usize] += 1;
            }
        }
        count
    }

    #[test]
    fn edits_then_two_reductions_carry_sources_and_make_none_up() {
        let source = sphere(0.1);
        let mut e = Editable::from(&source, 0.0).unwrap();
        // A cut makes vertices with no source, an offset of a region
        // copies the vertices on its edge with their sources, and a smooth
        // moves them all.
        let plane = Plane::through(
            Point::new(0.0, 0.0, 6.0),
            Direction::new(Vector::Z, T).unwrap(),
        );
        e.cut(&plane, Side::Below, Some(Cap::default())).unwrap();
        let marked: Vec<bool> = source.positions.iter().map(|p| p.x > 0.0).collect();
        let region = e.region_from_sources(&marked).unwrap();
        e.offset(0.5, Some(&region)).unwrap();
        e.smooth(&Smooth {
            size: 1.0,
            max_move: 0.05,
            keep_edges: Some(0.5),
            region: None,
        })
        .unwrap();
        check(&e);
        let unsourced = |e: &Editable| {
            present(e)
                .into_iter()
                .filter(|&v| e.source_vertices(v).is_empty())
                .count()
        };
        assert!(unsourced(&e) > 0);
        assert!(uses(&e).iter().any(|&n| n > 1), "the wall copied nothing");

        for target in [
            Target::Verified {
                tolerance: 0.05,
                share: 0.999,
                max_triangles: None,
            },
            Target::Triangles(400),
        ] {
            let before = uses(&e);
            let none_before = unsourced(&e);
            let report = e.simplify(target, T).unwrap();
            check(&e);
            assert!(report.collapsed > 0);
            // A collapse only merges: no source vertex behind more vertices
            // than before, none appearing, and no vertex without a source
            // that did not come only of such vertices.
            for (s, (&after, &was)) in uses(&e).iter().zip(&before).enumerate() {
                assert!(after <= was, "source vertex {s} behind {after}, was {was}");
            }
            assert!(unsourced(&e) <= none_before);
            for v in present(&e) {
                assert!(
                    e.source_vertices(v)
                        .iter()
                        .all(|&s| (s as usize) < source.positions.len())
                );
            }
        }
        assert!(e.to_triangulation().triangle_count() <= 400);
    }
}
