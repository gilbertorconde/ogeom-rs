//! A whole ball restated on a chart whose poles stand clear of the other
//! operand.
//!
//! A sphere's chart collapses at its poles and folds at its seam. Where the
//! other operand's boundary passes through a pole (a drill whose wall runs
//! along the ball's axis), every section through it starts and ends on the
//! pole, hugs the seam beside it, and leaves the chart no well-placed
//! crossing to split at. The same ball charted about another axis meets the
//! same boundary with its poles well away, and the sections are ordinary
//! loops. A ball that is one face (a seam and two poles, no trim) can be
//! restated that way exactly: the surface is the same set of points, only
//! the parameterization moves.

use ogeom_algo::{Built, History, make_sphere};
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Surface as _;
use ogeom_geom::{SurfaceGeometry, Transformable};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

/// How close a pole may come to the other operand's curved surfaces, as a
/// fraction of the ball's radius, before the chart is turned.
const TOUCHING: f64 = 1e-3;

/// A ball operand: the solid, its one face, and the sphere it lies on in
/// world space.
struct Ball {
    solid: Shape,
    face: Shape,
    centre: Point,
    radius: f64,
    axis: Vector,
    x: Vector,
}

/// `run` on the operands, where either is a whole ball with a pole on or
/// beside the other's curved surfaces, with that ball first charted about
/// an axis whose poles stand clear; the history then runs from the
/// operands as given. Where that refuses, and otherwise, `run` on the
/// operands as given.
pub(crate) fn with_poles_clear(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
    mut run: impl FnMut(&mut Model, &Shape, &Shape) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    if let Some((ta, tb, before)) = turned(model, a, b, tol)? {
        match run(model, &ta, &tb) {
            Ok(built) => return Ok(Built::new(built.shape, before.then(&built.history))),
            Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => return Err(e),
            Err(_) => {}
        }
    }
    run(model, a, b)
}

/// The operands with each whole ball whose pole stands on or beside the
/// other's surfaces charted about an axis whose poles stand clear, and the
/// history from the operands as given; `None` where no ball needs it.
fn turned(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape, History)>> {
    let mut before = History::new();
    let mut any = false;
    let (mut ta, mut tb) = (a.clone(), b.clone());
    for (operand, other) in [(&mut ta, b), (&mut tb, a)] {
        let Some(ball) = whole_ball(model, operand, tol)? else {
            continue;
        };
        let walls = curved_surfaces_of(model, other, tol)?;
        let clearance = |p: Point| -> f64 {
            walls
                .iter()
                .filter_map(|s| ogeom_algo::project_on_surface(s, p, 8, tol).ok())
                .map(|found| found.distance)
                .fold(f64::INFINITY, f64::min)
        };
        let poles = |axis: Vector| {
            clearance(ball.centre + axis * ball.radius)
                .min(clearance(ball.centre - axis * ball.radius))
        };
        if poles(ball.axis) > ball.radius * TOUCHING {
            continue;
        }
        // Candidate axes in the ball's own frame: its other two axes, the
        // diagonals of each pair, and the four diagonals of the cube. The
        // one whose nearer pole stands furthest from the other's surfaces.
        let y = ball.axis.cross(ball.x);
        let (x, z) = (ball.x, ball.axis);
        let mut candidates = vec![x, y];
        for (p, q) in [(x, y), (y, z), (z, x)] {
            candidates.push(p + q);
            candidates.push(p - q);
        }
        for (s, t) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
            candidates.push(x * s + y * t + z);
        }
        let mut best: Option<(f64, Vector)> = None;
        for candidate in candidates {
            let Ok(axis) = candidate.normalized(tol) else {
                continue;
            };
            let score = poles(axis);
            if best.is_none_or(|(s, _)| score > s) {
                best = Some((score, axis));
            }
        }
        let Some((score, axis)) = best else {
            continue;
        };
        if score <= ball.radius * TOUCHING {
            continue;
        }
        // The seam on the side of the new axis that stands furthest from
        // the other's surfaces, read at the seam's point on the equator.
        let across = if axis.cross(x).magnitude() > 0.5 {
            x
        } else {
            y
        };
        let Ok(u) = axis.cross(across).normalized(tol) else {
            continue;
        };
        let w = axis.cross(u);
        let seam = [u, -u, w, -w]
            .into_iter()
            .map(|d| (clearance(ball.centre + d * ball.radius), d))
            .fold(None, |kept: Option<(f64, Vector)>, (c, d)| match kept {
                Some((k, _)) if k >= c => kept,
                _ => Some((c, d)),
            })
            .map_or(u, |(_, d)| d);
        let frame = Frame::new(
            ball.centre,
            Direction::new(axis, tol)?,
            Direction::new(seam, tol)?,
            tol,
        )?;
        let made = make_sphere(model, frame, ball.radius, tol)?.shape;
        let Some(face) = explore_unique(model, &made, ShapeType::Face)?
            .first()
            .cloned()
        else {
            continue;
        };
        before.modify(&ball.solid, made.clone());
        for shell in explore_unique(model, &ball.solid, ShapeType::Shell)? {
            for new_shell in explore_unique(model, &made, ShapeType::Shell)? {
                before.modify(&shell, new_shell);
            }
        }
        before.modify(&ball.face, face);
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in explore_unique(model, &ball.solid, kind)? {
                before.delete(&sub);
            }
        }
        *operand = made;
        any = true;
    }
    Ok(any.then_some((ta, tb, before)))
}

/// The ball `shape` is, where it is a solid of one face on a sphere with no
/// trim (a seam and its two poles), its material inside, its placement
/// neither scaled nor mirrored.
fn whole_ball(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Option<Ball>> {
    if model.kind_of(shape)? != ShapeType::Solid {
        return Ok(None);
    }
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let [face] = faces.as_slice() else {
        return Ok(None);
    };
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let placement = face.transform(model.datums())?;
    if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 || !placement.preserves_handedness() {
        return Ok(None);
    }
    let Some(stored) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let SurfaceGeometry::Sphere(_) = stored else {
        return Ok(None);
    };
    let SurfaceGeometry::Sphere(placed) = stored.transformed(&placement, tol)? else {
        return Ok(None);
    };
    // The boundary: one edge with extent, used twice (the seam), and
    // degenerate edges (the poles).
    let mut seams = Vec::new();
    for wire in model.children_of(face)? {
        for edge in model.children_of(&wire)? {
            let Some(NodeData::Edge(e)) = model.node(&edge).map(|n| n.data()) else {
                return Ok(None);
            };
            if matches!(e.curve3d(), Some(EdgeRepr::Curve3d { .. })) {
                seams.push(edge.node());
            }
        }
    }
    if seams.len() != 2 || seams[0] != seams[1] {
        return Ok(None);
    }
    let sphere = placed.sphere();
    let frame = sphere.frame();
    // Material inside: the face's outward side, its natural normal or the
    // reverse, points away from the centre.
    let (u, v) = (0.5, 0.25);
    let at = placed.point_at(u, v, tol)?;
    let normal = placed.normal_at(u, v, tol)?.vector();
    let reversed = face.orientation() == ogeom_topo::Orientation::Reversed;
    let outward = (at - sphere.centre()).dot(normal) > 0.0;
    if outward == reversed {
        return Ok(None);
    }
    Ok(Some(Ball {
        solid: shape.clone(),
        face: face.clone(),
        centre: sphere.centre(),
        radius: sphere.radius(),
        axis: frame.z().vector(),
        x: frame.x().vector(),
    }))
}

/// The curved surfaces of `shape`'s faces, placed in world space. A plane
/// through a pole is left out: its section with the ball is a circle split
/// at the pole with exact pcurves on both sides, and the ball's own chart
/// arranges it as it stands.
fn curved_surfaces_of(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<SurfaceGeometry>> {
    let mut out = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        if matches!(surface, SurfaceGeometry::Plane(_)) {
            continue;
        }
        out.push(surface.transformed(&face.transform(model.datums())?, tol)?);
    }
    Ok(out)
}
