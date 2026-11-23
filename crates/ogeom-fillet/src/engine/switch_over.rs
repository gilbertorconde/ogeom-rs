//! Where a chain's ball changes walls off its vertex (ADR 0001, amendment
//! "a ball changing walls off the vertex"): two rounds of a request meet
//! tangent at a vertex on one shared face, their other faces meeting G1
//! along a seam. The ball rolling along the chain changes walls where its
//! contact crosses the seam, and the chain's rounds, joined in the plane
//! square to the chain at the vertex, close the corner only where that
//! crossing stands in that plane.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// Where at a vertex the ball rolling along two rounds of a request
/// crosses the seam between their walls.
#[derive(Clone)]
pub(crate) struct Crossing {
    /// The vertex's point.
    pub(crate) at: Point,
    /// The two edges, by their index in the request.
    pub(crate) pair: [usize; 2],
    /// The face the two share, by its index among the solid's faces.
    pub(crate) cap: usize,
    /// Each edge's other face, by the same index.
    pub(crate) walls: [usize; 2],
    /// The seam's curve, and its parameters at the vertex and where the
    /// ball's contact crosses it.
    pub(crate) seam: (Curve, (f64, f64)),
    /// How far along the chain that centre stands from the vertex's
    /// section.
    pub(crate) by: f64,
    /// The slot of the edge the centre stands over, past the vertex.
    pub(crate) over: usize,
    /// The rounds are convex.
    pub(crate) convex: bool,
}

/// The solid's faces, each with its edges.
pub(crate) struct Faces {
    pub(crate) faces: Vec<Shape>,
    pub(crate) edges: Vec<Vec<Shape>>,
}

impl Faces {
    pub(crate) fn of(model: &Model, solid: &Shape) -> OgeomResult<Self> {
        let faces = explore_unique(model, solid, ShapeType::Face)?;
        let edges = faces
            .iter()
            .map(|f| explore_unique(model, f, ShapeType::Edge))
            .collect::<OgeomResult<Vec<_>>>()?;
        Ok(Self { faces, edges })
    }

    fn holding(&self, edge: &Shape) -> Vec<usize> {
        (0..self.faces.len())
            .filter(|&k| self.edges[k].iter().any(|e| e.is_same(edge)))
            .collect()
    }
}

/// How far off a vertex's section a crossing may stand and still be the
/// section's: the bands' own fit.
pub(crate) fn off_by(tol: Tolerances) -> f64 {
    crate::marched::band_fit_target(tol) * 10.0
}

/// The crossing of `edges` (each with its radius) on `solid` that stands
/// farthest off its vertex's section, where one stands off it by more than
/// [`off_by`] (of equals, the first by its vertex's point, so the same in
/// any order); `None` where every one stands in it.
///
/// # Errors
///
/// As [`crossings`].
pub(crate) fn off_vertex(
    model: &Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Crossing>> {
    if edges.len() < 2 || model.kind_of(solid)? != ShapeType::Solid {
        return Ok(None);
    }
    let faces = Faces::of(model, solid)?;
    let off_by = off_by(tol);
    let mut found: Vec<Crossing> = crossings(model, solid, &faces, edges, tol)?
        .into_iter()
        .filter(|c| c.by > off_by)
        .collect();
    let key = |c: &Crossing| (c.by / off_by).round();
    found.sort_by(|a, b| {
        key(b).total_cmp(&key(a)).then(
            a.at.x
                .total_cmp(&b.at.x)
                .then(a.at.y.total_cmp(&b.at.y))
                .then(a.at.z.total_cmp(&b.at.z)),
        )
    });
    Ok(found.into_iter().next())
}

/// Every crossing of `edges` on `solid`: each vertex where two of them of
/// one sense and radius meet tangent on one shared face, their other faces
/// meeting along a seam there that the ball's contact reaches.
///
/// # Errors
///
/// As reading the edges' curves and the faces' surfaces.
pub(crate) fn crossings(
    model: &Model,
    solid: &Shape,
    faces: &Faces,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<Crossing>> {
    let mut read = Vec::with_capacity(edges.len());
    for (edge, radius) in edges {
        let held = faces.holding(edge);
        let ends = ogeom_algo::edge_vertices(model, edge)?;
        // An edge between tangent faces has no sense, and meets no seam.
        let sense = match crate::setback::edge_is_convex(model, solid, edge, *radius, tol) {
            Ok(sense) => sense,
            Err(ogeom_core::OgeomError::Construction(_)) => None,
            Err(e) => return Err(e),
        };
        read.push((held, ends, sense));
    }
    let mut found = Vec::new();
    for i in 0..edges.len() {
        for j in i + 1..edges.len() {
            let ((fi, Some((a, b)), Some(si)), (fj, Some((c, d)), Some(sj))) = (&read[i], &read[j])
            else {
                continue;
            };
            let radius = edges[i].1;
            if si != sj || (edges[j].1 - radius).abs() > tol.confusion() {
                continue;
            }
            let ([f0, f1], [g0, g1]) = (fi.as_slice(), fj.as_slice()) else {
                continue;
            };
            let shared: Vec<usize> = [*f0, *f1]
                .into_iter()
                .filter(|f| [*g0, *g1].contains(f))
                .collect();
            let [cap] = shared.as_slice() else {
                continue;
            };
            let wall_i = if f0 == cap { *f1 } else { *f0 };
            let wall_j = if g0 == cap { *g1 } else { *g0 };
            let ends = if a.is_same(b) { vec![a] } else { vec![a, b] };
            for v in ends {
                if !(v.is_same(c) || v.is_same(d)) {
                    continue;
                }
                if let Some(x) = crossing(
                    model,
                    faces,
                    [(i, &edges[i].0), (j, &edges[j].0)],
                    v,
                    [*cap, wall_i, wall_j],
                    radius,
                    *si,
                    tol,
                )? {
                    found.push(x);
                }
            }
        }
    }
    Ok(found)
}

/// The point of `curve` at `t`, and its unit tangent.
fn point_and_tangent(
    curve: &Curve,
    t: f64,
    tol: Tolerances,
) -> OgeomResult<(Point, Option<Vector>)> {
    let p = curve.point_at(t, tol)?;
    let d = curve.d1_at(t, tol)?;
    let m = d.magnitude();
    Ok((p, (m > tol.confusion()).then(|| d / m)))
}

/// The outward unit normal of a face (its placed surface and sign) at the
/// point of it nearest `p`, and that point.
pub(crate) fn outward_at(
    (surface, sign): &(SurfaceGeometry, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<(Point, Vector)> {
    let found = ogeom_geom::project_on_surface(surface, p, 16, tol)?;
    let (u, v) = found.parameters;
    Ok((found.point, surface.normal_at(u, v, tol)?.vector() * *sign))
}

/// Of `curve` over `range`, the parameter at whichever end lies nearer
/// `at`, and the other.
pub(crate) fn from_end(
    curve: &Curve,
    (t0, t1): (f64, f64),
    at: Point,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    Ok(
        if curve.point_at(t0, tol)?.distance(at) <= curve.point_at(t1, tol)?.distance(at) {
            (t0, t1)
        } else {
            (t1, t0)
        },
    )
}

/// The crossing at vertex `v` of the ball of `radius` (convex rounds'
/// within the material) rolling along `pair`. `None` where the two do not
/// meet tangent, the walls meet along no seam at `v`, or the seam ends
/// before the ball's contact reaches it.
#[allow(clippy::too_many_arguments, reason = "one crossing, read whole")]
fn crossing(
    model: &Model,
    faces: &Faces,
    pair: [(usize, &Shape); 2],
    v: &Shape,
    [cap, wall_i, wall_j]: [usize; 3],
    radius: f64,
    convex: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Crossing>> {
    let Some(at) = model
        .node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
    else {
        return Ok(None);
    };
    // Each edge's tangent at `v`, away from it along its own edge: the two
    // must run on tangent.
    let mut into = Vec::with_capacity(2);
    for (_, edge) in pair {
        let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            return Ok(None);
        };
        let (t, other) = from_end(&curve, range, at, tol)?;
        let (_, Some(d)) = point_and_tangent(&curve, t, tol)? else {
            return Ok(None);
        };
        into.push(if other > t { d } else { -d });
    }
    // A fitted band's foot runs on within its fit, a few tenths of a
    // milliradian.
    if into[0].cross(into[1]).magnitude() > 1e-2 || into[0].dot(into[1]) > 0.0 {
        return Ok(None);
    }
    // The seam: an edge at `v` the two walls share.
    let Some(seam) = faces.edges[wall_i]
        .iter()
        .find(|e| {
            faces.edges[wall_j].iter().any(|f| f.is_same(e))
                && ogeom_algo::edge_vertices(model, e)
                    .ok()
                    .flatten()
                    .is_some_and(|(a, b)| a.is_same(v) || b.is_same(v))
        })
        .cloned()
    else {
        return Ok(None);
    };
    let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, &seam, tol)? else {
        return Ok(None);
    };
    let (from, to) = from_end(&curve, range, at, tol)?;
    let hosts =
        [cap, wall_i, wall_j].map(|k| crate::engine::flat_faces::host(model, &faces.faces[k], tol));
    let [Ok(Some(cap_host)), Ok(Some(host_i)), Ok(Some(host_j))] = hosts else {
        return Ok(None);
    };
    // The ball's centre where it touches the walls at the seam's point at
    // `t`, and how far it stands from the cap less the radius: the ball
    // crosses the seam where that is nought.
    let side = if convex { -1.0 } else { 1.0 };
    let centre = |t: f64| -> OgeomResult<(Point, Point, f64)> {
        let q = curve.point_at(t, tol)?;
        let (_, ni) = outward_at(&host_i, q, tol)?;
        let (_, nj) = outward_at(&host_j, q, tol)?;
        let n = ni + nj;
        let n = n / n.magnitude().max(tol.confusion());
        let c = q + n * (side * radius);
        let (foot, up) = outward_at(&cap_host, c, tol)?;
        Ok((q, c, side * (c - foot).dot(up) - radius))
    };
    let steps = 64;
    let mut last = (from, centre(from)?.2);
    if last.1 >= 0.0 {
        return Ok(None);
    }
    for k in 1..=steps {
        let t = from + (to - from) * f64::from(k) / f64::from(steps);
        let h = centre(t)?.2;
        if h >= 0.0 {
            let (mut lo, mut hi) = (last.0, t);
            for _ in 0..60 {
                let m = 0.5 * (lo + hi);
                if centre(m)?.2 < 0.0 {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            let t = 0.5 * (lo + hi);
            let (_, c, _) = centre(t)?;
            let over = usize::from((c - at).dot(into[1]) > (c - at).dot(into[0]));
            return Ok(Some(Crossing {
                at,
                pair: [pair[0].0, pair[1].0],
                cap,
                walls: [wall_i, wall_j],
                seam: (curve, (from, t)),
                by: (c - at).dot(into[0]).abs(),
                over,
                convex,
            }));
        }
        last = (t, h);
    }
    Ok(None)
}
