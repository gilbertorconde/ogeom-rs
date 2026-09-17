//! Rounding a vertex: the ball-and-block tool at a corner whose faces one
//! ball touches.
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
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// vertex is not a vertex of the solid; if fewer than three planes pass
/// through it (a curved-edged corner is the N-support setback's, still
/// owed — docs/PARITY.md, fillet.edge-blends); if the planes share no
/// tangent ball — three that span always do, a square pyramid's four do,
/// a general N-edged vertex does not, and that vertex is owed the general
/// setback patch; or if the corner turns out concave, where a ball adds
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
        ogeom_bail!(
            Construction,
            "the {n} planes through this vertex share no tangent ball — the \
             nearest fit misses one by {residual}; the corner tool rounds \
             the vertex whose faces one ball touches, and the general \
             N-support setback is still owed — docs/PARITY.md, \
             fillet.edge-blends"
        );
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
