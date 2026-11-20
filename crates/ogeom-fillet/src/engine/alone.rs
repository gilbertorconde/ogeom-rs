//! Rounds of creases that meet only across convex turns, each built as it
//! is alone (ADR 0001, amendment "M5e": across a convex turn each band runs
//! out past the vertex as it does alone, and the cuts take their union).
//! The engine rounds each edge alone on the solid as given; what each takes
//! off is its tool, the faces its band leaves named for its edge, and the
//! tools are applied as the one-boolean fallback's. The solid that comes
//! out is what every round alone keeps.
//!
//! Concave creases meet the same way across a turn convex in the void (a
//! pocket floor's corner, its upright left sharp): each fill runs on past
//! the vertex as it does alone, and the fills take their union (amendment
//! "a concave band ends against its end face").

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};

/// `once` rounded on `solid` where every edge is a crease the engine
/// rounds alone, all convex or all concave, and the edges meet one another
/// only at vertices of three edges, across a turn on the face they share
/// convex in the material (concave creases: in the void); `None`
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
    if once.len() < 2 {
        return Ok(None);
    }
    let mut convex = None;
    for (edge, radius) in once {
        let this = crate::setback::edge_is_convex(model, solid, edge, *radius, tol)?;
        if this.is_none() || convex.is_some_and(|c| Some(c) != this) {
            return Ok(None);
        }
        convex = this;
    }
    let Some(convex) = convex else {
        return Ok(None);
    };
    if !meet_across_convex_turns(model, solid, once, convex, tol)? {
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
        // A convex crease's round takes material and adds none; a concave
        // one's adds and takes none.
        let (Ok(taken), Ok(added)) = (
            ogeom_bool::cut(&mut scratch, solid, &alone.shape, tol),
            ogeom_bool::cut(&mut scratch, &alone.shape, solid, tol),
        ) else {
            return Ok(None);
        };
        let (taken, added) = if convex {
            (taken, added)
        } else {
            (added, taken)
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
            sense: if convex { Sense::Cut } else { Sense::Fuse },
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
/// convex turn on the shared face), or for concave creases (not `convex`)
/// into the void before it, and at least one pair meets so.
fn meet_across_convex_turns(
    model: &Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    convex: bool,
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
                if !convex_turn(model, solid, [a, b], v, convex, tol)? {
                    return Ok(false);
                }
                met = true;
            }
        }
    }
    Ok(met)
}

/// Whether `edges` meet at `v` across a convex turn: one face shared, and
/// each leaving `v` into the material behind the other's other face, or
/// where not `convex` into the void before it.
fn convex_turn(
    model: &Model,
    solid: &Shape,
    edges: [&Shape; 2],
    v: &Shape,
    convex: bool,
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
    let side = if convex { -1.0 } else { 1.0 };
    Ok(side * la.dot(nb) > ang && side * lb.dot(na) > ang)
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
