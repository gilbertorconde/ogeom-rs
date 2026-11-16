//! Rounds of creases that meet only across convex turns, each built as it
//! is alone (ADR 0001, amendment "M5e": across a convex turn each band runs
//! out past the vertex as it does alone, and the cuts take their union).
//! The engine rounds each edge alone on the solid as given; what each takes
//! off is its tool, the faces its band leaves named for its edge, and the
//! tools are applied as the one-boolean fallback's. The solid that comes
//! out is what every round alone keeps.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};

/// `once` rounded on `solid` where every edge is a convex crease the
/// engine rounds alone and the edges meet one another only at vertices of
/// three edges, across a convex turn on the face they share; `None`
/// otherwise, or where a boolean refuses or the result is not valid.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled, and reading the
/// edges.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if once.len() < 2 || !meet_across_convex_turns(model, solid, once, tol)? {
        return Ok(None);
    }
    let mut scratch = model.clone();
    let mut tools = Vec::with_capacity(once.len());
    for (edge, radius) in once {
        let alone = match round(&mut scratch, solid, &[(edge.clone(), *radius)], tol) {
            Ok(Some(alone)) => alone,
            Err(e @ OgeomError::Cancelled) => return Err(e),
            _ => return Ok(None),
        };
        // A convex crease's round takes material and adds none.
        let (Ok(taken), Ok(added)) = (
            ogeom_bool::cut(&mut scratch, solid, &alone.shape, tol),
            ogeom_bool::cut(&mut scratch, &alone.shape, solid, tol),
        ) else {
            return Ok(None);
        };
        if !explore_unique(&scratch, &added.shape, ShapeType::Solid)?.is_empty()
            || explore_unique(&scratch, &taken.shape, ShapeType::Solid)?.is_empty()
        {
            return Ok(None);
        }
        let mut faces = Vec::new();
        for band in alone.history.generated(edge) {
            if scratch.kind_of(band)? != ShapeType::Face {
                continue;
            }
            for face in taken.history.trace(band) {
                faces.push((face.clone(), edge.clone()));
            }
        }
        if faces.is_empty() {
            return Ok(None);
        }
        tools.push(Tool {
            solid: taken.shape,
            sense: Sense::Cut,
            faces,
        });
    }
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(&scratch, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    *model = scratch;
    Ok(Some(built))
}

/// Whether the edges of `once` that share a vertex meet only there, at a
/// vertex of three edges whose third is not asked, sharing one face, each
/// leaving the vertex into the material behind the other's other face (a
/// convex turn on the shared face), and at least one pair meets so.
fn meet_across_convex_turns(
    model: &Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let asked = |e: &Shape| once.iter().any(|(o, _)| o.is_same(e));
    let mut met = false;
    for (i, (a, _)) in once.iter().enumerate() {
        for (b, _) in &once[i + 1..] {
            let (Some(ea), Some(eb)) = (
                ogeom_algo::edge_vertices(model, a)?,
                ogeom_algo::edge_vertices(model, b)?,
            ) else {
                return Ok(false);
            };
            let shared: Vec<Shape> = [ea.0.clone(), ea.1.clone()]
                .into_iter()
                .filter(|v| v.is_same(&eb.0) || v.is_same(&eb.1))
                .collect();
            for v in &shared {
                let around = ogeom_topo::ancestors_of(model, solid, v, ShapeType::Edge)?;
                let mut edges: Vec<Shape> = Vec::new();
                for e in around {
                    if !edges.iter().any(|f| f.is_same(&e)) {
                        edges.push(e);
                    }
                }
                if edges.len() != 3 || edges.iter().filter(|e| asked(e)).count() != 2 {
                    return Ok(false);
                }
                if !convex_turn(model, solid, [a, b], v, tol)? {
                    return Ok(false);
                }
                met = true;
            }
        }
    }
    Ok(met)
}

/// Whether `edges` meet at `v` across a convex turn: one face shared, and
/// each leaving `v` into the material behind the other's other face.
fn convex_turn(
    model: &Model,
    solid: &Shape,
    edges: [&Shape; 2],
    v: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let faces = |e: &Shape| -> OgeomResult<Vec<Shape>> {
        let mut out: Vec<Shape> = Vec::new();
        for f in ogeom_topo::ancestors_of(model, solid, e, ShapeType::Face)? {
            if !out.iter().any(|g| g.is_same(&f)) {
                out.push(f);
            }
        }
        Ok(out)
    };
    let (fa, fb) = (faces(edges[0])?, faces(edges[1])?);
    let shared: Vec<&Shape> = fa
        .iter()
        .filter(|f| fb.iter().any(|g| g.is_same(f)))
        .collect();
    if shared.len() != 1 || fa.len() != 2 || fb.len() != 2 {
        return Ok(false);
    }
    let other = |fs: &[Shape]| fs.iter().find(|f| !f.is_same(shared[0])).cloned();
    let (Some(wall_a), Some(wall_b)) = (other(&fa), other(&fb)) else {
        return Ok(false);
    };
    let Some(at) = model
        .node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .map(|p| v.transform(model.datums()).map(|t| t.apply(p)))
        .transpose()?
    else {
        return Ok(false);
    };
    let ang = tol.angular() * 100.0;
    let (Some(la), Some(lb)) = (
        leaving(model, edges[0], at, tol)?,
        leaving(model, edges[1], at, tol)?,
    ) else {
        return Ok(false);
    };
    let (Some(na), Some(nb)) = (
        outward(model, &wall_a, at, tol)?,
        outward(model, &wall_b, at, tol)?,
    ) else {
        return Ok(false);
    };
    Ok(la.dot(nb) < -ang && lb.dot(na) < -ang)
}

/// The unit tangent of `edge` leaving `at`, one of its ends.
fn leaving(model: &Model, edge: &Shape, at: Point, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let Some((curve, (t0, t1))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    let (t, sense) =
        if curve.point_at(t0, tol)?.distance(at) <= curve.point_at(t1, tol)?.distance(at) {
            (t0, 1.0)
        } else {
            (t1, -1.0)
        };
    let d = curve.d1_at(t, tol)? * sense;
    let m = d.magnitude();
    Ok((m > tol.confusion()).then(|| d / m))
}

/// The outward unit normal of `face` at `at`, a point of it.
fn outward(model: &Model, face: &Shape, at: Point, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let Some(n) = crate::marched::face_normal_near(model, face, at, tol)? else {
        return Ok(None);
    };
    Ok(Some(
        if face.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        },
    ))
}
