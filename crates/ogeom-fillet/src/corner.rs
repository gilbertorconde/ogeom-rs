//! Rounding a vertex: the ball-and-block tool at a planar corner.
//!
//! The ball's centre may sit anywhere a radius in from every host plane.
//! At a corner one ball touches, that region's tip is a point and the
//! rounded corner is one spherical patch; at a corner no single ball
//! touches, the tip is a few points joined by short ridges, and the
//! rounded corner is a sphere at each and a cylinder along each ridge —
//! the exact envelope of the rolling ball, where the setback family fits
//! a plate through the bands' ends instead.
//!
//! *Elsewhere:* the vertex blend of `ChFi3d`'s setback family.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// Round a solid's vertex with a ball of `radius`.
///
/// The construction is the corner family's centre of gravity, promoted from
/// the B2 proof: the corner block spanned by the edges less the ball seated
/// a radius in from every host plane is exactly the spike a rounded corner
/// sheds, and the general boolean does the shedding. The block's faces
/// through the ball's centre are square to the edges, so each edge's flush
/// band ends on the ball's rim: three sequential fillets at a box corner
/// followed by this call round the vertex the setback way, and the
/// `b2_three_fillets_and_the_corner_tool_round_the_vertex` pin measures
/// the result against a closed form. At a vertex of more edges the corner
/// goes first and the fillets follow — four bands built before the corner
/// crash into each other at a pyramid's apex — and the block is the
/// polyhedron of the N host planes and the N planes square to the edges.
///
/// A vertex whose planes share no tangent ball — a rectangular pyramid's
/// apex, any general N-edged vertex — is rounded by the envelope of every
/// ball a radius in from all of them: a sphere at each vertex of the
/// region the ball's centre may occupy and a cylinder along each ridge
/// between two of them, each cut with its own compartment, the spheres
/// by the one-ball tool on their three planes and the ridges by the flush
/// fillet of a virtual crease. The compartments meet on the planes
/// square to the ridges, cap to cap.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// vertex is not a vertex of the solid; if fewer than three planes pass
/// through it (a curved-edged corner is the setback family's, still
/// owed — docs/PARITY.md, fillet.edge-blends); if the region the ball's
/// centre may occupy has a tip vertex touching more than three planes
/// without one ball touching all of the corner's, or with other than three
/// edges leaving it; or if the corner turns out concave, where a ball adds
/// material instead of shedding it and a tool built from a cut cannot say
/// so. The corner may be oblique: the block is then the hexahedron bounded
/// by the host planes and the three planes through the ball's centre
/// square to the edges.
pub fn round_vertex(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(vertex)? != ShapeType::Vertex {
        ogeom_bail!(Construction, "round_vertex rounds a vertex");
    }
    if radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend radius must be a positive distance");
    }
    let Some(corner) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
    else {
        ogeom_bail!(Construction, "the vertex holds no point");
    };

    // The corner's frame comes from the planes that pass through the
    // vertex's point — not from the vertex's own adjacency, which the very
    // sequence this tool serves destroys: after three fillets the tip is
    // consumed, but the three shrunk planes still contain the corner, and
    // still say exactly which corner it was. The vertex argument may
    // therefore come from an earlier state of the solid — the sharp box's
    // corner captured before the fillets — and anchors the history. Each
    // plane's material side is read off its face, which the fillets shrink
    // but never turn over.
    let mut m: Vec<Vector> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        if !matches!(surface, ogeom_geom::SurfaceGeometry::Plane(_)) {
            continue;
        }
        let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
        if (corner - origin).dot(outward).abs() > tol.confusion() * 100.0 {
            continue;
        }
        // One vote per plane: coplanar trims share it.
        if m.iter()
            .any(|n| n.cross(outward).magnitude() < tol.angular() * 10.0)
        {
            continue;
        }
        m.push(-outward);
    }
    let n = m.len();
    if n < 3 {
        ogeom_bail!(
            Construction,
            "round_vertex speaks the planar corner: at least three planes \
             must pass through the vertex, found {n}; the curved-edged corner \
             is the setback family's, still owed — docs/PARITY.md, \
             fillet.edge-blends"
        );
    }
    // The ball's centre: the point a radius in from every plane. Three
    // planes that span always hold one; more only when they share a
    // tangent ball, which the least-squares fit's residual tells — a
    // square pyramid's apex does, a general N-edged vertex does not, and
    // that vertex is owed the general setback patch instead.
    let centre_for = |m: &[Vector]| -> Option<(Point, f64)> {
        // m_k · (c − corner) = radius over all k: solved directly for
        // three planes, through the normal equations for more. The direct
        // solve is kept for three not for speed but for its last bit — the
        // boolean's paving at the touch points is still sensitive to an
        // ulp of the centre, and the oblique corner closes on the direct
        // solve's value.
        let mut a = [[0.0_f64; 3]; 3];
        let mut b = [0.0_f64; 3];
        if m.len() == 3 {
            for (i, mk) in m.iter().enumerate() {
                a[i] = [mk.x, mk.y, mk.z];
                b[i] = radius;
            }
        } else {
            for mk in m {
                let v = [mk.x, mk.y, mk.z];
                for i in 0..3 {
                    b[i] += radius * v[i];
                    for j in 0..3 {
                        a[i][j] += v[i] * v[j];
                    }
                }
            }
        }
        let det3 = |a: &[[f64; 3]; 3]| -> f64 {
            a[0][0].mul_add(
                a[1][1].mul_add(a[2][2], -(a[1][2] * a[2][1])),
                -a[0][1].mul_add(
                    a[1][0].mul_add(a[2][2], -(a[1][2] * a[2][0])),
                    -(a[0][2] * a[1][0].mul_add(a[2][1], -(a[1][1] * a[2][0]))),
                ),
            )
        };
        let det = det3(&a);
        if !det.is_finite() || det.abs() <= 1e-12 {
            return None;
        }
        let mut x = [0.0_f64; 3];
        for (k, xk) in x.iter_mut().enumerate() {
            let mut ak = a;
            for (row, bi) in ak.iter_mut().zip(b) {
                row[k] = bi;
            }
            *xk = det3(&ak) / det;
        }
        let centre = corner + Vector::new(x[0], x[1], x[2]);
        let residual = m
            .iter()
            .map(|mk| (mk.dot(centre - corner) - radius).abs())
            .fold(0.0, f64::max);
        Some((centre, residual))
    };
    let Some((far, residual)) = centre_for(&m) else {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex do not span a corner"
        );
    };
    if residual > tol.confusion() * 10.0 {
        // No one ball touches every face: the rounded corner is the
        // envelope of every ball a radius in from all of them, which is
        // more than one sphere. Its pieces are read off the region the
        // ball's centre may occupy.
        return setback_corner(model, solid, vertex, corner, &m, radius, tol);
    }
    // A concave vertex puts the ball's centre outside the material: a ball
    // there adds material instead of shedding it, and a tool built from a
    // cut cannot say so.
    let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    if boundary.holds(model, far, tol)? != ogeom_algo::Containment::In {
        ogeom_bail!(
            Construction,
            "no material a radius in from every face at this vertex; a \
             concave vertex gains a ball instead of shedding one, and this \
             tool cannot round it"
        );
    }
    // The corner's edges: where two of the planes meet along a ray that
    // lies inside all the others, directed from the vertex into the solid.
    // A simple convex corner has one edge per plane, and the edges chain the
    // planes into a ring round the vertex.
    let mut edges: Vec<(usize, usize, Direction)> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let cross = m[i].cross(m[j]);
            let length = cross.magnitude();
            if length <= tol.angular() * 10.0 {
                continue;
            }
            for sign in [1.0, -1.0] {
                let d = cross * (sign / length);
                let inside = (0..n)
                    .filter(|&k| k != i && k != j)
                    .all(|k| m[k].dot(d) >= -tol.angular() * 10.0);
                if inside {
                    edges.push((i, j, Direction::new(d, tol)?));
                    break;
                }
            }
        }
    }
    if edges.len() != n {
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex meet along {} edges, not {n}; \
             the corner is not the simple convex one this tool speaks",
            edges.len()
        );
    }
    let mut ring: Vec<usize> = vec![0];
    let mut last_edge: Option<usize> = None;
    loop {
        let here = *ring.last().unwrap_or(&0);
        let Some((e, other)) = edges.iter().enumerate().find_map(|(e, &(i, j, _))| {
            if Some(e) == last_edge {
                None
            } else if i == here {
                Some((e, j))
            } else if j == here {
                Some((e, i))
            } else {
                None
            }
        }) else {
            ogeom_bail!(
                Construction,
                "the planes through this vertex do not chain into a ring round it"
            );
        };
        last_edge = Some(e);
        if other == 0 {
            break;
        }
        if ring.contains(&other) || ring.len() >= n {
            ogeom_bail!(
                Construction,
                "the planes through this vertex do not chain into a ring round it"
            );
        }
        ring.push(other);
    }
    if ring.len() != n {
        ogeom_bail!(
            Construction,
            "the planes through this vertex do not chain into a ring round it"
        );
    }
    // Plane `ring[k]` and plane `ring[k+1]` meet along edge k.
    let d: Vec<Direction> = (0..n)
        .map(|k| {
            let (p, q) = (ring[k], ring[(k + 1) % n]);
            edges
                .iter()
                .find(|&&(i, j, _)| (i == p && j == q) || (i == q && j == p))
                .map(|&(_, _, dir)| dir)
        })
        .collect::<Option<_>>()
        .ok_or_else(|| {
            ogeom_core::ogeom_err!(Construction, "the corner's edges lost their ring")
        })?;
    let inward_of = |k: usize| m[ring[k % n]];

    // The block: the corner bounded by its N host planes and, through the
    // ball's centre, the N planes square to its edges — where each band's
    // circle and the ball's own rim coincide, so the cut ends the band and
    // starts the patch on one curve. On a square corner it is the box of
    // side `radius`; on an oblique one a hexahedron; at a pyramid's apex a
    // polyhedron of eight faces. Its corners: the vertex, the foot of each
    // edge on its cutting plane, on each host plane the point where the
    // ball touches it, and the centre itself.
    //
    // The ball on the corner's axes: a pole at one corner of the patch it
    // leaves and its seam meridian out past the block through the first
    // edge — the pole axis is the inward normal of the host plane holding
    // the first two edges, the third edge itself on a square corner.
    //
    // Which edge is first and which way the ring runs is the tool's
    // labelling, and the solid it builds is the same for all 2N. The
    // boolean is not yet indifferent to it: the charts the block's faces
    // and the ball wear decide where a rim is exact and where fitted, where
    // a seam falls against a patch arc, and at an oblique corner two of the
    // six labellings still die in the cut. So the tool is offered on each
    // labelling in turn and the first that closes stands — every one of
    // them is the same exact construction — and the corner is refused by
    // name only when none does. A failed attempt's nodes stay in the model
    // unreferenced, under their own operation. The boolean closing all of
    // them is owed (docs/PARITY.md, fillet.edge-blends).
    let attempt = |model: &mut Model, start: usize, reverse: bool| -> OgeomResult<Built> {
        // Edge t of the labelling and the host plane holding edges t and t+1.
        let edge_at = |t: usize| -> usize {
            if reverse {
                (start + n - t % n) % n
            } else {
                (start + t) % n
            }
        };
        let dl: Vec<Direction> = (0..n).map(|t| d[edge_at(t)]).collect();
        let ml: Vec<Vector> = (0..n)
            .map(|t| {
                // Forward: edges t, t+1 are E_s+t, E_s+t+1, both on plane
                // ring[s+t+1]. Reverse: E_s−t, E_s−t−1, both on ring[s−t].
                if reverse {
                    inward_of((start + n - t % n) % n)
                } else {
                    inward_of(edge_at(t) + 1)
                }
            })
            .collect();
        let along =
            |t: usize| -> Point { corner + dl[t].vector() * (far - corner).dot(dl[t].vector()) };
        let touch = |t: usize| -> Point { far - ml[t] * radius };
        if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
            eprintln!(
                "CORNER attempt start {start} reverse {reverse} far {far:?} d {:?} m {ml:?}",
                dl.iter().map(|x| x.vector()).collect::<Vec<_>>()
            );
        }
        model.begin_operation();
        let mut points: Vec<Point> = vec![corner];
        points.extend((0..n).map(along));
        points.extend((0..n).map(touch));
        points.push(far);
        let mut rings: Vec<Vec<usize>> = Vec::with_capacity(2 * n);
        for t in 0..n {
            // The host face on plane t: vertex, edge t's foot, the touch
            // point, edge t+1's foot.
            rings.push(vec![0, 1 + t, 1 + n + t, 1 + (t + 1) % n]);
            // The face square to edge t: its foot, the touch points either
            // side, the centre.
            rings.push(vec![1 + t, 1 + n + (t + n - 1) % n, 1 + 2 * n, 1 + n + t]);
        }
        let block = ogeom_algo::make_polyhedron(model, &points, &rings, tol)?.shape;
        let ball_frame = Frame::new(far, Direction::new(ml[0], tol)?, dl[0], tol)?;
        let ball = ogeom_algo::make_sphere(model, ball_frame, radius, tol)?.shape;
        let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
        if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
            let faces = ogeom_topo::explore_unique(model, &tool.shape, ShapeType::Face)?;
            let kinds: Vec<String> = faces
                .iter()
                .map(|f| {
                    let kind = model
                        .node(f)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| model.geometry().surface(d.surface))
                        .map_or("?", |sg| match sg {
                            ogeom_geom::SurfaceGeometry::Plane(_) => "plane",
                            ogeom_geom::SurfaceGeometry::Sphere(_) => "sphere",
                            _ => "other",
                        });
                    let edges = ogeom_topo::explore_unique(model, f, ShapeType::Edge)
                        .map_or(0, |e| e.len());
                    format!("{kind}/{edges}")
                })
                .collect();
            eprintln!(
                "CORNER tool for start {start} reverse {reverse}: {} faces {kinds:?}",
                faces.len()
            );
        }
        let rounded = ogeom_bool::cut(model, solid, &tool.shape, tol)?;
        Ok(Built {
            shape: rounded.shape,
            history: tool.history.then(&rounded.history),
        })
    };
    let mut outcome: Option<Built> = None;
    let mut last: Option<ogeom_core::OgeomError> = None;
    // Forensics: one labelling only, by index, for the boolean's benefit.
    let forced: Option<usize> = std::env::var("OGEOM_CORNER_LABELLING")
        .ok()
        .and_then(|v| v.parse().ok());
    let labellings = (0..2 * n).map(|index| (index % n, index >= n));
    for (index, (start, reverse)) in labellings.enumerate() {
        if forced.is_some_and(|f| f != index) {
            continue;
        }
        match attempt(model, start, reverse) {
            Ok(built) => {
                outcome = Some(built);
                break;
            }
            Err(err) => last = Some(err),
        }
    }
    let Some(rounded) = outcome else {
        ogeom_bail!(
            NotDone,
            "the corner tool's cut closed on none of the corner's {} \
             labellings; the last said: {}",
            2 * n,
            last.map_or_else(String::new, |e| e.to_string())
        );
    };
    let mut built = rounded;
    built.history.modify(vertex, built.shape.clone());
    Ok(built)
}

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

/// Round a convex planar vertex whose faces no single ball touches.
///
/// The ball's centre may sit anywhere a radius in from every host plane:
/// a convex region whose tip, at a vertex one ball touches, is a single
/// point, and otherwise a few points joined by short edges — a rectangular
/// pyramid's apex has two, joined along the two long slopes. The rounded
/// corner is the envelope of every ball centred in that region: a sphere
/// at each tip vertex, a cylinder along each edge between them, and the
/// host planes themselves elsewhere. Exact, constant-radius, and what the
/// rolling ball leaves — where a plate would be fitted through the bands'
/// ends instead.
///
/// Each piece is cut with its own block: at a tip vertex, the corner
/// bounded by its three planes and the three planes through the centre
/// square to its edges, less the ball — the one-ball tool exactly — and
/// along an edge between two centres, the prism over the kite of the
/// virtual crease, the two touch points and the centre, between the two
/// planes square to the edge, less the cylinder. The compartments tile
/// the corner and meet on the planes square to the edges, where each
/// sphere's rim and the cylinder's end coincide, so the cuts consume one
/// another's flush faces cap to cap. An original edge leaves its tip
/// vertex along a ray, and the plane square to it there is where the
/// edge's band ends when the flush fillets follow.
fn setback_corner(
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
                let Some(centre) = ball_centre(corner, &[m[i], m[j], m[k]], radius) else {
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
    if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
        for tip in &tips {
            eprintln!("TIP {:?} planes {:?}", tip.centre, tip.planes);
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
    // tangent — a ridge to the tip vertex there — and a ray otherwise, the
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
        // In ring order, consecutive edges sharing a plane, so that host t
        // holds edges t and t+1 — the labelling the ball's frame is read
        // off, its pole along a host and its seam out through an edge.
        let shared = |x: &[usize; 2], y: &[usize; 2]| -> Option<usize> {
            x.iter().copied().find(|p| y.contains(p))
        };
        let mut ring: Vec<(Vector, [usize; 2])> = vec![edges[0]];
        let mut used = vec![false; k];
        used[0] = true;
        while ring.len() < k {
            let last = ring[ring.len() - 1];
            let Some(next) = (0..k).find(|&i| {
                !used[i]
                    && shared(&edges[i].1, &last.1).is_some_and(|p| {
                        // The plane shared with the previous edge is not
                        // the one shared with the edge before that.
                        ring.len() < 2 || shared(&ring[ring.len() - 2].1, &last.1) != Some(p)
                    })
            }) else {
                ogeom_bail!(Construction, "the tip vertex's edges do not chain")
            };
            used[next] = true;
            ring.push(edges[next]);
        }
        let hosts: Vec<Vector> = (0..k)
            .map(|t| {
                shared(&ring[t].1, &ring[(t + 1) % k].1)
                    .map(|p| m[p])
                    .ok_or_else(|| {
                        ogeom_core::ogeom_err!(Construction, "the tip vertex's edges do not chain")
                    })
            })
            .collect::<OgeomResult<_>>()?;
        let directions: Vec<Vector> = ring.iter().map(|(d, _)| *d).collect();
        // Whether edge t is a ridge, whose rim plane a later cut's cap
        // stands in: a pole along a plane holding a ridge would put both
        // poles in that cap, and the cap's circle would have no chart
        // image. Labellings whose pole host holds no ridge go first.
        let is_ridge: Vec<bool> = ring
            .iter()
            .map(|(d, _)| !rays[index].iter().any(|(r, _)| r.dot(*d) > 1.0 - 1e-9))
            .collect();
        let mut walls: Vec<(Vector, Point)> = tip.planes.iter().map(|&p| (m[p], corner)).collect();
        for direction in &directions {
            walls.push((-*direction, tip.centre));
        }
        let tool = ball_block(
            model,
            &walls,
            corner,
            tip.centre,
            &hosts,
            &directions,
            &is_ridge,
            radius,
            tol,
        )?;
        let cut = ogeom_bool::cut(model, &rounded, &tool.shape, tol).map_err(|e| {
            ogeom_core::ogeom_err!(
                NotDone,
                "the cut by the ball's block at tip vertex {index} failed: {e}"
            )
        })?;
        history = history.then(&tool.history).then(&cut.history);
        rounded = cut.shape;
    }
    // The cylinder along each ridge: the flush fillet of a virtual crease —
    // the line the two planes it touches would meet along — between the
    // planes square to the ridge through its two centres, which are the
    // caps' own planes. The planar fillet builds that wedge face by face,
    // band, legs and caps, and melts it; the caps meet the spheres' rims
    // on the planes the vertex compartments already cut.
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
        let cut = crate::fillet::seated_fillet(model, &rounded, &seat, radius, None, tol)
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

/// A convex polytope from the half-spaces that bound it, each given by
/// its inward normal and a point on its plane.
///
/// Its corners are the feasible meetings of three planes; each plane's
/// face is those of its corners that lie on it, walked round the face's
/// centroid. The polyhedron builder checks what this hands it — planarity,
/// every edge shared by two faces — so a set of half-spaces that bounds
/// nothing, or bounds a sliver, is refused rather than built.
fn convex_block(
    model: &mut Model,
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let slack = tol.confusion() * 100.0;
    let n = walls.len();
    let mut points: Vec<Point> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(p) = planes_meet(&walls[i], &walls[j], &walls[k]) else {
                    continue;
                };
                if walls
                    .iter()
                    .any(|(normal, on)| normal.dot(p - *on) < -slack)
                {
                    continue;
                }
                if points.iter().any(|q| q.distance(p) <= slack) {
                    continue;
                }
                points.push(p);
            }
        }
    }
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for (normal, on) in walls {
        let mine: Vec<usize> = (0..points.len())
            .filter(|&i| normal.dot(points[i] - *on).abs() <= slack)
            .collect();
        if mine.len() < 3 {
            continue;
        }
        let centroid = mine
            .iter()
            .fold(Vector::ZERO, |acc, &i| acc + (points[i] - Point::ORIGIN))
            * (1.0 / f64::from(u32::try_from(mine.len()).unwrap_or(u32::MAX)));
        let centroid = Point::ORIGIN + centroid;
        let axis = normal.normalized(tol)?;
        let first = (points[mine[0]] - centroid).normalized(tol)?;
        let second = axis.cross(first);
        let mut ordered: Vec<(f64, usize)> = mine
            .iter()
            .map(|&i| {
                let v = points[i] - centroid;
                (v.dot(second).atan2(v.dot(first)), i)
            })
            .collect();
        ordered.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
        rings.push(ordered.into_iter().map(|(_, i)| i).collect());
    }
    Ok(ogeom_algo::make_polyhedron(model, &points, &rings, tol)?.shape)
}

/// Where three planes meet, or `None` where they do not span.
fn planes_meet(a: &(Vector, Point), b: &(Vector, Point), c: &(Vector, Point)) -> Option<Point> {
    let m = [a.0, b.0, c.0];
    let rhs = [
        a.0.dot(a.1 - Point::ORIGIN),
        b.0.dot(b.1 - Point::ORIGIN),
        c.0.dot(c.1 - Point::ORIGIN),
    ];
    let x = solve3(&m, rhs)?;
    Some(Point::ORIGIN + x)
}

/// `m_k · x = rhs_k` for three rows, or `None` where they do not span.
fn solve3(m: &[Vector; 3], rhs: [f64; 3]) -> Option<Vector> {
    let a: [[f64; 3]; 3] = std::array::from_fn(|i| [m[i].x, m[i].y, m[i].z]);
    let det3 = |a: &[[f64; 3]; 3]| -> f64 {
        a[0][0].mul_add(
            a[1][1].mul_add(a[2][2], -(a[1][2] * a[2][1])),
            -a[0][1].mul_add(
                a[1][0].mul_add(a[2][2], -(a[1][2] * a[2][0])),
                -(a[0][2] * a[1][0].mul_add(a[2][1], -(a[1][1] * a[2][0]))),
            ),
        )
    };
    let det = det3(&a);
    if !det.is_finite() || det.abs() <= 1e-12 {
        return None;
    }
    let mut x = [0.0_f64; 3];
    for (k, xk) in x.iter_mut().enumerate() {
        let mut ak = a;
        for (row, r) in ak.iter_mut().zip(rhs) {
            row[k] = r;
        }
        *xk = det3(&ak) / det;
    }
    Some(Vector::new(x[0], x[1], x[2]))
}

/// The point a radius in from three planes through `corner`, or `None`
/// where they do not span.
fn ball_centre(corner: Point, m: &[Vector; 3], radius: f64) -> Option<Point> {
    solve3(m, [radius; 3]).map(|x| corner + x)
}

/// A compartment less the ball centred in it: the block from its walls,
/// cut by the sphere at `far`.
///
/// Host `t` holds edges `t` and `t + 1`. The ball's pole stands along a
/// host's normal and its seam meridian runs out through the first of that
/// host's edges — in a rim plane, so the seam doubles as a trim rather
/// than crossing a patch — and every labelling is offered in turn, the
/// first that closes standing. A cut that closes on a tool reaching past
/// its own block is a wrong tool, not a closed one, and is passed over
/// too.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn ball_block(
    model: &mut Model,
    walls: &[(Vector, Point)],
    corner: Point,
    far: Point,
    hosts: &[Vector],
    directions: &[Vector],
    is_ridge: &[bool],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = directions.len();
    let mut last: Option<ogeom_core::OgeomError> = None;
    // The ball's frames, best first. Along a ridge, the ridge's rim plane
    // is the sphere's equator, both poles stand outside the patch — the
    // patch lies on the corner's side of every rim plane, the poles on
    // the ridge's own axis either side of it — and a seam meridian turned
    // away from the corner never crosses it; the ridge fillet's cap, in
    // that same rim plane, then meets the sphere on a circle the chart
    // images exactly. Then the classic labellings: the pole along host t,
    // the seam out through edge t, hosts holding no ridge first.
    let mut frames: Vec<(Vector, Vector)> = Vec::new();
    for t in 0..n {
        if is_ridge[t] {
            let axis = directions[t];
            let away = far - corner;
            let seam = away - axis * away.dot(axis);
            if seam.magnitude() > tol.angular() * 10.0 {
                frames.push((axis, seam));
            }
        }
    }
    // Host t holds edges t and t+1; a labelling's pole host is `start`
    // forward and `start − 1` reversed.
    let pole_host = |index: usize| -> usize {
        let (start, reverse) = (index % n, index >= n);
        if reverse {
            (start + 2 * n - 1) % n
        } else {
            start
        }
    };
    let clean = |index: usize| -> bool {
        let h = pole_host(index);
        !is_ridge[h] && !is_ridge[(h + 1) % n]
    };
    let mut order: Vec<usize> = (0..2 * n).collect();
    order.sort_by_key(|&index| !clean(index));
    frames.extend(
        order
            .iter()
            .map(|&index| (hosts[pole_host(index)], directions[index % n])),
    );
    for (index, (pole, seam)) in frames.into_iter().enumerate() {
        model.begin_operation();
        let attempt = (|| -> OgeomResult<Built> {
            let block = convex_block(model, walls, tol)?;
            let block_bound = ogeom_algo::shape_bounds(model, &block, tol)?;
            let ball_frame = Frame::new(
                far,
                Direction::new(pole, tol)?,
                Direction::new(seam, tol)?,
                tol,
            )?;
            let ball = ogeom_algo::make_sphere(model, ball_frame, radius, tol)?.shape;
            let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
            let tool_bound = ogeom_algo::shape_bounds(model, &tool.shape, tol)?;
            let reach = tol.confusion() * 1e3;
            let (Some(block_lo), Some(block_hi), Some(tool_lo), Some(tool_hi)) = (
                block_bound.low(),
                block_bound.high(),
                tool_bound.low(),
                tool_bound.high(),
            ) else {
                ogeom_bail!(NotDone, "the ball's cut left no tool");
            };
            if tool_lo.x < block_lo.x - reach
                || tool_lo.y < block_lo.y - reach
                || tool_lo.z < block_lo.z - reach
                || tool_hi.x > block_hi.x + reach
                || tool_hi.y > block_hi.y + reach
                || tool_hi.z > block_hi.z + reach
            {
                ogeom_bail!(
                    NotDone,
                    "the ball's cut left a tool reaching past its own block"
                );
            }
            Ok(tool)
        })();
        match attempt {
            Ok(tool) => {
                if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
                    eprintln!("BALL frame {index} closed");
                }
                return Ok(tool);
            }
            Err(err) => {
                if std::env::var_os("OGEOM_DEBUG_CORNER").is_some() {
                    eprintln!("BALL frame {index} pole {pole:?} seam {seam:?}: {err}");
                }
                last = Some(err);
            }
        }
    }
    ogeom_bail!(
        NotDone,
        "the corner tool's block closed on none of its {} frames; the last said: {}",
        2 * n,
        last.map_or_else(String::new, |e| e.to_string())
    )
}
