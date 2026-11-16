//! The exact envelope of the feasible ball centres at a planar vertex.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// A vertex of the region the ball's centre may occupy: a point a radius
/// in from three or more of the host planes and at least a radius from the
/// rest, with the planes it touches.
struct TipVertex {
    centre: Point,
    planes: Vec<usize>,
}

/// An edge of that region's tip between two of its vertices: the ball
/// rolling from one to the other touches two planes all the way, and
/// sweeps a cylinder.
struct Ridge {
    from: usize,
    to: usize,
    planes: [usize; 2],
}

/// Round a convex planar vertex by the offset region's spherical and cylindrical envelope.
///
/// The ball's centre may sit anywhere a radius in from every host plane:
/// a convex region whose tip, at a vertex one ball touches, is a single
/// point, and otherwise a few points joined by short edges: a rectangular
/// pyramid's apex has two, joined along the two long slopes. The rounded
/// corner is the envelope of every ball centred in that region: a sphere
/// at each tip vertex, a cylinder along each edge between them, and the
/// host planes themselves elsewhere. Exact, constant-radius, and what the
/// rolling ball leaves, where a plate would be fitted through the bands'
/// ends instead.
///
/// Each piece is cut with its own block: at a tip vertex, the corner
/// bounded by its three planes and the three planes through the centre
/// square to its edges, less the ball (the one-ball tool exactly) and
/// along an edge between two centres, the prism over the kite of the
/// virtual crease, the two touch points and the centre, between the two
/// planes square to the edge, less the cylinder. The compartments tile
/// the corner and meet on the planes square to the edges, where each
/// sphere's rim and the cylinder's end coincide, so the cuts consume one
/// another's flush faces cap to cap. An original edge leaves its tip
/// vertex along a ray, and the plane square to it there is where the
/// edge's band ends when the flush fillets follow.
pub(super) fn round(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    corner: Point,
    m: &[Vector],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = m.len();
    let slack = tol.confusion() * 10.0;
    // The tip's vertices: every triple that meets at a point a radius in
    // from all the planes, merged where several triples name one point.
    let mut tips: Vec<TipVertex> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(centre) = super::vertex_geometry::solve3(&[m[i], m[j], m[k]], [radius; 3])
                    .map(|v| corner + v)
                else {
                    continue;
                };
                let feasible = m.iter().all(|mk| mk.dot(centre - corner) >= radius - slack);
                if !feasible {
                    continue;
                }
                if tips.iter().any(|t| t.centre.distance(centre) <= slack) {
                    continue;
                }
                let planes: Vec<usize> = (0..n)
                    .filter(|&l| (m[l].dot(centre - corner) - radius).abs() <= slack)
                    .collect();
                tips.push(TipVertex { centre, planes });
            }
        }
    }
    if tips.is_empty() {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex hold no ball a radius in from \
             all of them"
        );
    }
    // Each tip vertex's edges: along every pair of its planes, the way the
    // rest of its planes allow; bounded where another plane becomes
    // tangent (a ridge to the tip vertex there) and a ray otherwise, the
    // way an original edge's band runs. A ray that lies along no original
    // edge, or a vertex with other than three edges, is a corner this tool
    // does not speak.
    let mut ridges: Vec<Ridge> = Vec::new();
    let mut rays: Vec<Vec<(Vector, [usize; 2])>> = vec![Vec::new(); tips.len()];
    for (index, tip) in tips.iter().enumerate() {
        for (a, b) in tip
            .planes
            .iter()
            .flat_map(|&a| tip.planes.iter().map(move |&b| (a, b)))
        {
            if a >= b {
                continue;
            }
            let cross = m[a].cross(m[b]);
            let length = cross.magnitude();
            if length <= tol.angular() * 10.0 {
                continue;
            }
            let Some(direction) = [1.0, -1.0]
                .into_iter()
                .map(|sign| cross * (sign / length))
                .find(|d| {
                    tip.planes
                        .iter()
                        .filter(|&&l| l != a && l != b)
                        .all(|&l| m[l].dot(*d) >= -tol.angular() * 10.0)
                })
            else {
                continue;
            };
            // The first other plane the ball meets rolling this way.
            let mut nearest: Option<(f64, usize)> = None;
            for (l, ml) in m.iter().enumerate() {
                if tip.planes.contains(&l) {
                    continue;
                }
                let rate = ml.dot(direction);
                if rate >= -tol.angular() {
                    continue;
                }
                let t = (radius - ml.dot(tip.centre - corner)) / rate;
                if t > slack && nearest.is_none_or(|(held, _)| t < held) {
                    nearest = Some((t, l));
                }
            }
            match nearest {
                Some((t, _)) => {
                    let end = tip.centre + direction * t;
                    let Some(to) = tips.iter().position(|o| o.centre.distance(end) <= slack) else {
                        ogeom_bail!(
                            Construction,
                            "the ball rolling between two of this vertex's faces meets a \
                             third where no tip vertex stands"
                        );
                    };
                    if index < to {
                        ridges.push(Ridge {
                            from: index,
                            to,
                            planes: [a, b],
                        });
                    }
                }
                None => rays[index].push((direction, [a, b])),
            }
        }
    }
    let concave = {
        let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
        let mut concave = false;
        for tip in &tips {
            if boundary.holds(model, tip.centre, tol)? != ogeom_algo::Containment::In {
                concave = true;
            }
        }
        concave
    };
    if concave {
        ogeom_bail!(
            Construction,
            "no material a radius in from every face at this vertex; a \
             concave vertex gains a ball instead of shedding one, and this \
             tool cannot round it"
        );
    }

    model.begin_operation();
    let mut rounded = solid.clone();
    let mut history = ogeom_algo::History::new();
    // The sphere at each tip vertex, with its own compartment: the wedge of
    // its planes, cut by the plane square to each of its edges through the
    // centre. A ridge's plane faces the corner, and the corner itself lies
    // in the ridge's compartment, not this one.
    for (index, tip) in tips.iter().enumerate() {
        // The vertex's edges, each on two of its planes: its rays and the
        // ridges that start or end here.
        let mut edges: Vec<(Vector, [usize; 2])> = rays[index].clone();
        for ridge in &ridges {
            let other = if ridge.from == index {
                ridge.to
            } else if ridge.to == index {
                ridge.from
            } else {
                continue;
            };
            let direction = (tips[other].centre - tip.centre).normalized(tol)?;
            edges.push((direction, ridge.planes));
        }
        let k = tip.planes.len();
        if edges.len() != k {
            ogeom_bail!(
                Construction,
                "a tip vertex of this corner touches {k} planes along {} edges, not \
                 {k}; the corner is not the simple convex one this tool speaks",
                edges.len()
            );
        }
        let mut walls: Vec<(Vector, Point)> = tip.planes.iter().map(|&p| (m[p], corner)).collect();
        walls.extend(edges.iter().map(|(direction, _)| (-*direction, tip.centre)));
        let tool = super::vertex_tools::ball_polytope(model, &walls, tip.centre, radius, tol)?;
        let cut = ogeom_bool::cut(model, &rounded, &tool, tol).map_err(|e| {
            ogeom_core::ogeom_err!(
                NotDone,
                "the cut by the ball's block at tip vertex {index} failed: {e}"
            )
        })?;
        history = history.then(&cut.history);
        rounded = cut.shape;
    }
    // Each ridge uses one section over its virtual crease between the
    // two centre planes. Its end arcs coincide with the spherical pieces'
    // rims, so their cuts consume the adjoining caps.
    for ridge in &ridges {
        let (v1, v2) = (tips[ridge.from].centre, tips[ridge.to].centre);
        let along = (v2 - v1).normalized(tol)?;
        let [a, c] = ridge.planes;
        let on_crease = |v: Point| corner + along * (v - corner).dot(along);
        let face_on = |model: &Model, inward: Vector| -> OgeomResult<Shape> {
            for face in ogeom_topo::explore_unique(model, &rounded, ShapeType::Face)? {
                let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                    continue;
                };
                if !matches!(
                    model.geometry().surface(data.surface),
                    Some(ogeom_geom::SurfaceGeometry::Plane(_))
                ) {
                    continue;
                }
                let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
                if (corner - origin).dot(outward).abs() <= tol.confusion() * 100.0
                    && outward.cross(inward).magnitude() < tol.angular() * 10.0
                    && outward.dot(inward) < 0.0
                {
                    return Ok(face);
                }
            }
            ogeom_bail!(
                Construction,
                "the ridge's host plane is no longer a face of the solid"
            )
        };
        let faces = [face_on(model, m[a])?, face_on(model, m[c])?];
        let seat = crate::support::Seat {
            start: on_crease(v1),
            end: on_crease(v2),
            along,
            normals: [-m[a], -m[c]],
            faces,
            convex: true,
        };
        let cut = super::vertex_tools::cut_strip(model, &rounded, &seat, radius, tol)
            .map_err(|e| ogeom_core::ogeom_err!(NotDone, "the ridge's flush fillet failed: {e}"))?;
        history = history.then(&cut.history);
        rounded = cut.shape;
    }
    history.modify(vertex, rounded.clone());
    Ok(Built {
        shape: rounded,
        history,
    })
}
