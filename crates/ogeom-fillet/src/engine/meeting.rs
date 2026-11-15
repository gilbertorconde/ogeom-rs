//! Rounds of straight convex creases meeting at inside corners, carried on
//! past each other (ADR 0001, amendment "M5e, rounds carried on through an
//! inside corner").
//!
//! Where two rounded top creases meet at a corner that is reflex in plan
//! (an inside corner of a notch, a pocket's rim), each band runs on past
//! the corner as its own crease's round, and the solid keeps whatever lies
//! under either band: a point goes only where it lies above both bands
//! carried on, or above one and past the other's wall, where the corner is
//! open. Off square in plan, or with walls leaning, the bands then meet
//! along their own intersection and the corner's cell past both walls,
//! which no flush cap or patch at the vertex holds.
//!
//! Each crease's band carried on is its spandrel prism: its two faces'
//! wedge, far longer than the solid, less that wedge rounded by the engine
//! (an exact cylinder). Each crease's tool is its spandrel where, at each
//! end shared with another asked crease, it also lies past that crease's
//! wall or within its spandrel; the solid is cut by the tools (the one
//! boolean of amendment "M5a"). A crease end with no asked neighbour must
//! carry its band out of the solid there, so the band ends flush on the
//! faces it meets.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{Curved, read};
use super::fallback::{Sense, Tool, polytope};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// One asked crease: its edge, its ends, its two faces with their outward
/// normals and offsets, and its radius.
struct Crease {
    shape: Shape,
    ends: (Id, Id),
    a: Point,
    b: Point,
    faces: [usize; 2],
    planes: [(Vector, f64); 2],
    radius: f64,
}

/// The rounds of `asked` on `solid` where every asked edge is a straight
/// convex crease between planes, two or more of them meet at inside
/// corners of a face they share, and every crease end met by no other
/// asked crease runs out of the solid; `None` otherwise, or where a tool
/// or the boolean will not build.
///
/// # Errors
///
/// As reading the solid and building the tools, other than refusals,
/// which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let ang = tol.angular();
    if asked.len() < 2 {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut creases = Vec::with_capacity(asked.len());
    for (edge, radius) in asked {
        let Some(crease) = crease_of(&curved, id(edge), edge, *radius, tol)? else {
            return Ok(None);
        };
        creases.push(crease);
    }
    // The extent of the solid, for the wedges and the walls' half-spaces.
    let (mut lo, mut hi) = (
        Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in curved.points.values() {
        let v = *p - Point::ORIGIN;
        lo = Vector::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vector::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let reach = (hi - lo).magnitude() + 4.0 * creases.iter().map(|c| c.radius).fold(0.0, f64::max);
    let mid = (lo + hi) * 0.5;
    let boxed = |planes: &mut Vec<(Vector, f64)>| {
        for axis in [Vector::X, Vector::Y, Vector::Z] {
            let m = axis.dot(mid);
            planes.push((axis, m + reach));
            planes.push((-axis, -(m - reach)));
        }
    };
    // Each crease's neighbour at each end: the asked crease sharing that
    // vertex, with the face both walk and the neighbour's other face (its
    // wall); the corner reflex in plan.
    let mut neighbours: Vec<[Option<(usize, usize)>; 2]> = vec![[None, None]; creases.len()];
    let mut met = false;
    for (i, c) in creases.iter().enumerate() {
        for (k, v) in [c.ends.0, c.ends.1].into_iter().enumerate() {
            let others: Vec<usize> = (0..creases.len())
                .filter(|j| *j != i && (creases[*j].ends.0 == v || creases[*j].ends.1 == v))
                .collect();
            match others[..] {
                [] => {
                    // Carried on past a free end, the band leaves the
                    // solid: a point just past the vertex, a little inside
                    // the wedge, stands outside it.
                    let (at, away) = if k == 0 {
                        (c.a, c.a - c.b)
                    } else {
                        (c.b, c.b - c.a)
                    };
                    let away = away / away.magnitude();
                    let inward = -(c.planes[0].0 + c.planes[1].0);
                    let inward = inward / inward.magnitude();
                    let probe = at + away * (c.radius * 0.5) + inward * (c.radius * 0.1);
                    if ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?
                        != ogeom_algo::Containment::Out
                    {
                        return Ok(None);
                    }
                }
                [j] => {
                    let n = &creases[j];
                    let shared: Vec<usize> = c
                        .faces
                        .iter()
                        .copied()
                        .filter(|f| n.faces.contains(f))
                        .collect();
                    let [top] = shared[..] else {
                        return Ok(None);
                    };
                    let Some(wall) = n.faces.iter().position(|f| *f != top) else {
                        return Ok(None);
                    };
                    // An inside corner: this crease, leaving the vertex,
                    // runs on the open side of the neighbour's wall, and the
                    // neighbour on the open side of this one's.
                    let along = |x: &Crease, at: Id| {
                        let d = if x.ends.0 == at { x.b - x.a } else { x.a - x.b };
                        d / d.magnitude()
                    };
                    let own_wall = c.planes[usize::from(c.faces[0] == top)].0;
                    if along(c, v).dot(n.planes[wall].0) <= ang || along(n, v).dot(own_wall) <= ang
                    {
                        return Ok(None);
                    }
                    neighbours[i][k] = Some((j, wall));
                    met = true;
                }
                _ => return Ok(None),
            }
        }
    }
    if !met {
        return Ok(None);
    }
    let mut scratch = model.clone();
    // Each crease's spandrel carried on, its band's face named.
    let mut spandrels: Vec<(Shape, Vec<Shape>)> = Vec::with_capacity(creases.len());
    for c in &creases {
        let Some(spandrel) = spandrel(&mut scratch, c, reach, tol)? else {
            return Ok(None);
        };
        spandrels.push(spandrel);
    }
    // Each crease's wall's far side, with its spandrel.
    let mut open: Vec<Option<Shape>> = vec![None; creases.len()];
    let mut tools = Vec::with_capacity(creases.len());
    for (i, c) in creases.iter().enumerate() {
        let (mut tool, mut named) = spandrels[i].clone();
        for (j, wall) in neighbours[i].iter().flatten() {
            let either = match &open[*j] {
                Some(s) => s.clone(),
                None => {
                    let (n, d) = creases[*j].planes[*wall];
                    let mut planes = vec![(-n, -d)];
                    boxed(&mut planes);
                    let Some(past) = polytope(&mut scratch, &planes, tol)? else {
                        return Ok(None);
                    };
                    let Ok(both) = ogeom_bool::fuse(&mut scratch, &past, &spandrels[*j].0, tol)
                    else {
                        return Ok(None);
                    };
                    open[*j] = Some(both.shape.clone());
                    both.shape
                }
            };
            let Ok(kept) = ogeom_bool::common(&mut scratch, &tool, &either, tol) else {
                return Ok(None);
            };
            named = named
                .iter()
                .flat_map(|f| kept.history.trace(f).to_vec())
                .collect();
            tool = kept.shape;
        }
        if explore_unique(&scratch, &tool, ShapeType::Solid)?.is_empty() || named.is_empty() {
            return Ok(None);
        }
        tools.push(Tool {
            solid: tool,
            sense: Sense::Cut,
            faces: named.into_iter().map(|f| (f, c.shape.clone())).collect(),
        });
    }
    let edges: Vec<Shape> = creases.iter().map(|c| c.shape.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// `edge` as a crease: straight, open, between two planar faces, convex;
/// `None` otherwise.
fn crease_of(
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Crease>> {
    let Some(c) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    if c.ends.0 == c.ends.1 || !matches!(c.curve, Curve::Line(_)) {
        return Ok(None);
    }
    let (a, b) = (curved.points[&c.ends.0], curved.points[&c.ends.1]);
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    if [fa, fb]
        .iter()
        .any(|f| !matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
    {
        return Ok(None);
    }
    let mid = a.midpoint(b);
    let (Some(na), Some(nb)) = (outward(curved, fa, mid), outward(curved, fb, mid)) else {
        return Ok(None);
    };
    // Convex: each face's points away from the edge stand behind the
    // other's plane (the material between them).
    let behind = |f: usize, n: Vector| {
        curved.faces[f].loops.iter().flatten().any(|(x, _)| {
            let g = &curved.edges[x];
            [g.ends.0, g.ends.1]
                .iter()
                .any(|v| (curved.points[v] - mid).dot(n) < -tol.confusion())
        })
    };
    if !behind(fb, na) || !behind(fa, nb) || na.cross(nb).magnitude() <= tol.angular() {
        return Ok(None);
    }
    Ok(Some(Crease {
        shape: shape.clone(),
        ends: c.ends,
        a,
        b,
        faces: [fa, fb],
        planes: [
            (na, na.dot(a - Point::ORIGIN)),
            (nb, nb.dot(a - Point::ORIGIN)),
        ],
        radius,
    }))
}

/// The spandrel of `crease` carried on along its line `reach` past both
/// ends: the wedge of its two faces' planes, capped square to the crease
/// and behind it, less that wedge with its edge rounded by the engine (an
/// exact cylinder); with the band's faces in it.
fn spandrel(
    model: &mut Model,
    crease: &Crease,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let eps = tol.confusion();
    let along = (crease.b - crease.a) / crease.a.distance(crease.b);
    let o = crease.a;
    let o_along = along.dot(o - Point::ORIGIN);
    let bisect = crease.planes[0].0 + crease.planes[1].0;
    let bisect = bisect / bisect.magnitude();
    let mut planes = crease.planes.to_vec();
    planes.push((along, o_along + crease.a.distance(crease.b) + reach));
    planes.push((-along, -(o_along - reach)));
    // Behind the crease, a few radii in: the wedge's back.
    planes.push((
        -bisect,
        -(bisect.dot(o - Point::ORIGIN) - 4.0 * crease.radius),
    ));
    let Some(wedge) = polytope(model, &planes, tol)? else {
        return Ok(None);
    };
    // The wedge's edge where its two faces' planes meet (an edge's own
    // vertices may stand off that line by a little).
    let line_point = |p: Point| {
        crease
            .planes
            .iter()
            .all(|(n, d)| (n.dot(p - Point::ORIGIN) - d).abs() <= 10.0 * eps)
    };
    let point = |model: &Model, v: &Shape| {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
    };
    let mut edge = None;
    for e in explore_unique(model, &wedge, ShapeType::Edge)? {
        if let Some((p, q)) = ogeom_algo::edge_vertices(model, &e)?
            && let (Some(p), Some(q)) = (point(model, &p), point(model, &q))
            && line_point(p)
            && line_point(q)
        {
            edge = Some(e);
        }
    }
    let Some(edge) = edge else {
        return Ok(None);
    };
    let request = super::Request {
        blends: vec![super::request::Blend {
            edge: edge.clone(),
            section: super::Section::Round {
                radius: crease.radius,
            },
        }],
        mitre: false,
        splits: Vec::new(),
        corners: Vec::new(),
    };
    let Ok(Some(rounded)) = super::run(model, &wedge, &request, tol) else {
        return Ok(None);
    };
    let band: Vec<Shape> = rounded.history.generated(&edge).to_vec();
    let Ok(cut) = ogeom_bool::cut(model, &wedge, &rounded.shape, tol) else {
        return Ok(None);
    };
    // The band's faces as the spandrel holds them: traced through the cut,
    // or where it recorded none, the faces on the band's own surface.
    let kept = explore_unique(model, &cut.shape, ShapeType::Face)?;
    let mut named: Vec<Shape> = band
        .iter()
        .flat_map(|b| cut.history.trace(b).to_vec())
        .filter(|f| kept.iter().any(|k| k.is_same(f)))
        .collect();
    if named.is_empty() {
        let surface = |s: &Shape| {
            model
                .node(s)
                .and_then(|n| n.data().as_face().map(|d| d.surface))
        };
        named = kept
            .into_iter()
            .filter(|f| {
                band.iter()
                    .any(|b| surface(b).is_some() && surface(b) == surface(f))
            })
            .collect();
    }
    Ok(Some((cut.shape, named)))
}
