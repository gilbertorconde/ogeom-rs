//! Small faces and small solids, removed where a caller names the size
//! below which they are noise.
//!
//! A face the size of a tolerance is no feature: an exchange round trip or
//! a boolean at the edge of its resolution leaves spots (a face that fits
//! in a ball of that size) and strips (a face whose two long sides stand
//! that close along their length). Neither is removed by [`fix_shape`],
//! which drops nothing the model says is there; this is the step that does,
//! at a size the caller states.
//!
//! A spot collapses to a point: its edges go and its vertices become one.
//! A strip collapses to an edge: its short ends go, one long side stands
//! for both, and its neighbours meet on it. Either way the neighbours'
//! geometry does not move; the vertices and edges that close the gap widen
//! to own it, as the containment rule asks, and trims are fitted where an
//! edge now bounds a face it did not.
//!
//! [`fix_shape`]: crate::fix_shape()

use std::collections::HashMap;

use ogeom_algo::{Built, History, edge_vertices, volume_properties};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_math::Point;
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Model, Shape, ShapeType, TShapeId, explore_unique};

use crate::{Reshape, fix_face_pcurves};

/// What [`fix_small_faces`] removed.
#[derive(Debug, Clone)]
pub struct SmallFaces {
    /// The shape without them, and what became of every input.
    pub built: Built,
    /// Faces collapsed to a point.
    pub spots: usize,
    /// Faces collapsed to an edge.
    pub strips: usize,
}

/// Remove every face of `shape` smaller than `size`: a spot (it fits in a
/// ball of that diameter) collapses to a point, and a strip (two long
/// sides within `size` of each other along their length, every other side
/// shorter than `size`) collapses to one of its long sides.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `size` is not a distance or the rebuilt shape comes out empty.
pub fn fix_small_faces(
    model: &mut Model,
    shape: &Shape,
    size: f64,
    tol: Tolerances,
) -> OgeomResult<SmallFaces> {
    if !size.is_finite() || size <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a small-face size of {size} is not a distance"
        );
    }
    // Collapses merge vertices and substitute edges node by node, and a
    // node placed twice (a prism's far cap is its near cap moved) would be
    // merged at both places at once. Such a shape is baked first, every
    // occurrence its own node in world space.
    let placed_twice = explore_unique(model, shape, ShapeType::Edge)?
        .iter()
        .any(|e| !e.location().is_identity());
    let start = if placed_twice && model.kind_of(shape)? == ShapeType::Solid {
        ogeom_algo::baked_shape(model, shape, tol)?
    } else {
        Built::new(shape.clone(), History::identity())
    };
    // Passes until none removes anything: a face beside one just removed
    // waits for the next pass, where it stands as the last one left it.
    let mut total = SmallFaces {
        built: start,
        spots: 0,
        strips: 0,
    };
    for _ in 0..16 {
        let pass = one_pass(model, &total.built.shape, size, tol)?;
        if pass.spots + pass.strips == 0 {
            break;
        }
        total = SmallFaces {
            built: Built::new(
                pass.built.shape,
                total.built.history.then(&pass.built.history),
            ),
            spots: total.spots + pass.spots,
            strips: total.strips + pass.strips,
        };
    }
    Ok(total)
}

fn one_pass(
    model: &mut Model,
    shape: &Shape,
    size: f64,
    tol: Tolerances,
) -> OgeomResult<SmallFaces> {
    let mut reshape = Reshape::new();
    // Vertex merges, resolved through earlier ones: a survivor's survivor.
    let mut survivor: HashMap<TShapeId, Shape> = HashMap::new();
    let root = |survivor: &HashMap<TShapeId, Shape>, v: &Shape| -> Shape {
        let mut current = v.clone();
        while let Some(next) = survivor.get(&current.node()) {
            if next.node() == current.node() {
                break;
            }
            current = next.clone();
        }
        current
    };
    let mut gone: Vec<TShapeId> = Vec::new();
    let (mut spots, mut strips) = (0, 0);

    for face in explore_unique(model, shape, ShapeType::Face)? {
        let edges = explore_unique(model, &face, ShapeType::Edge)?;
        if edges.iter().any(|e| gone.contains(&e.node())) {
            // A neighbour already collapsed onto this face's boundary; the
            // next pass sees the rebuilt face.
            continue;
        }
        let mut samples: Vec<(Shape, Vec<Point>)> = Vec::with_capacity(edges.len());
        for edge in &edges {
            samples.push((edge.clone(), edge_points(model, edge, tol)?));
        }
        let all: Vec<Point> = samples
            .iter()
            .flat_map(|(_, p)| p.iter().copied())
            .collect();
        if all.is_empty() {
            continue;
        }
        let spread = all
            .iter()
            .flat_map(|p| all.iter().map(move |q| p.distance(*q)))
            .fold(0.0_f64, f64::max);

        if spread < size {
            // A spot: every edge goes, every vertex becomes the first.
            let mut vertices = explore_unique(model, &face, ShapeType::Vertex)?.into_iter();
            let Some(first) = vertices.next() else {
                continue;
            };
            let keep = root(&survivor, &first);
            for v in vertices {
                let drop = root(&survivor, &v);
                if !drop.is_same(&keep) {
                    survivor.insert(drop.node(), keep.clone());
                }
            }
            for edge in &edges {
                reshape.remove(edge);
                gone.push(edge.node());
            }
            reshape.remove(&face);
            spots += 1;
            continue;
        }

        // A strip: two long sides, every other side short, the long ones
        // within `size` of each other all along.
        let long: Vec<usize> = (0..samples.len())
            .filter(|&i| polyline_length(&samples[i].1) >= size)
            .collect();
        let [a, b] = long.as_slice() else {
            continue;
        };
        if explore_unique(model, &face, ShapeType::Wire)?.len() != 1 {
            continue;
        }
        let (side_a, side_b) = (&samples[*a], &samples[*b]);
        let apart = hausdorff(&side_a.1, &side_b.1);
        if apart >= size {
            continue;
        }
        let (Some((a0, a1)), Some((b0, b1))) = (
            edge_vertices(model, &side_a.0)?,
            edge_vertices(model, &side_b.0)?,
        ) else {
            continue;
        };
        let at = |v: &Shape| -> OgeomResult<Point> {
            let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Construction, "a vertex holds no point");
            };
            Ok(v.transform(model.datums())?.apply(data.point))
        };
        // Which of b's ends faces which of a's: the same way round or the
        // other.
        let same_way = at(&a0)?.distance(at(&b0)?) + at(&a1)?.distance(at(&b1)?)
            <= at(&a0)?.distance(at(&b1)?) + at(&a1)?.distance(at(&b0)?);
        let (to0, to1) = if same_way { (&b0, &b1) } else { (&b1, &b0) };
        for (from, to) in [(&a0, to0), (&a1, to1)] {
            let (keep, drop) = (root(&survivor, to), root(&survivor, from));
            if !drop.is_same(&keep) {
                survivor.insert(drop.node(), keep.clone());
            }
        }
        for (i, (edge, _)) in samples.iter().enumerate() {
            if i != *b {
                gone.push(edge.node());
            }
            if i != *a && i != *b {
                reshape.remove(edge);
            }
        }
        // Side a stands for side b wherever a neighbour held it, and owns
        // the width it spans.
        let stand_in = if same_way {
            side_b.0.clone()
        } else {
            side_b.0.reversed()
        };
        model.widen(&side_b.0, Tolerance::new(apart + tol.confusion())?)?;
        reshape.replace(&side_a.0, stand_in);
        reshape.remove(&face);
        strips += 1;
    }

    // Each absorbed vertex's survivor widens to reach where it stood.
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let to = root(&survivor, &vertex);
        if to.is_same(&vertex) {
            continue;
        }
        let from = {
            let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            (
                vertex.transform(model.datums())?.apply(data.point),
                data.tolerance.get(),
            )
        };
        let Some(data) = model.node(&to).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        let here = to.transform(model.datums())?.apply(data.point);
        model.widen(
            &to,
            Tolerance::new((from.0.distance(here) + from.1).max(tol.confusion()))?,
        )?;
        reshape.replace(&vertex, to);
    }

    if reshape.is_empty() {
        return Ok(SmallFaces {
            built: Built::new(shape.clone(), History::identity()),
            spots,
            strips,
        });
    }
    let built = reshape.apply(model, shape)?;
    // Edges now bounding faces they did not are given their trims there.
    for face in explore_unique(model, &built.shape, ShapeType::Face)? {
        fix_face_pcurves(model, &face, tol.confusion() * 1e7, tol)?;
    }
    ogeom_algo::restore_containment(model, &built.shape)?;
    Ok(SmallFaces {
        built,
        spots,
        strips,
    })
}

/// Remove every solid of `shape` enclosing less than `volume`: the debris
/// an exchange or a boolean leaves beside the part.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `volume` is not positive, a solid's volume cannot be measured, or every
/// solid would go.
pub fn remove_small_solids(
    model: &mut Model,
    shape: &Shape,
    volume: f64,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    if !volume.is_finite() || volume <= 0.0 {
        ogeom_bail!(
            Construction,
            "a small-solid volume of {volume} is not positive"
        );
    }
    let solids = explore_unique(model, shape, ShapeType::Solid)?;
    let mut reshape = Reshape::new();
    let mut removed = 0;
    for solid in &solids {
        let measured = volume_properties(model, solid, Deflection::default(), tol)?.mass;
        if measured < volume {
            reshape.remove(solid);
            removed += 1;
        }
    }
    if removed == 0 {
        return Ok((Built::new(shape.clone(), History::identity()), 0));
    }
    if removed == solids.len() {
        ogeom_bail!(
            Construction,
            "every solid is smaller than {volume}; removing them leaves nothing"
        );
    }
    Ok((reshape.apply(model, shape)?, removed))
}

/// Points along an edge, placed.
fn edge_points(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(Vec::new());
    };
    if data.degenerate {
        return Ok(Vec::new());
    }
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(Vec::new());
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        return Ok(Vec::new());
    };
    let placement = edge.transform(model.datums())?;
    let mut out = Vec::with_capacity(17);
    for i in 0..=16 {
        let t = range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
        out.push(placement.apply(geometry.point_at(t, tol)?));
    }
    Ok(out)
}

fn polyline_length(points: &[Point]) -> f64 {
    points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The largest distance from either polyline to the other.
fn hausdorff(a: &[Point], b: &[Point]) -> f64 {
    let nearest = |p: Point, line: &[Point]| -> f64 {
        line.windows(2)
            .map(|w| {
                let d = w[1] - w[0];
                let dd = d.dot(d);
                let s = if dd > 0.0 {
                    ((p - w[0]).dot(d) / dd).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                p.distance(w[0] + d * s)
            })
            .fold(f64::INFINITY, f64::min)
    };
    let one = a.iter().map(|p| nearest(*p, b)).fold(0.0_f64, f64::max);
    let other = b.iter().map(|p| nearest(*p, a)).fold(0.0_f64, f64::max);
    one.max(other)
}
