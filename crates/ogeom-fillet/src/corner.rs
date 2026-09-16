//! Rounding a vertex: the ball-and-block tool at a trihedral corner.
//!
//! *Elsewhere:* the vertex blend of `ChFi3d`'s setback family.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// Round a solid's vertex with a ball of `radius`.
///
/// The construction is the corner family's centre of gravity, promoted from
/// the B2 proof: the corner block spanned by the three edges less the ball
/// seated at the same origin is exactly the spike a rounded corner sheds,
/// and the general boolean does the shedding. Three sequential fillets at a
/// box corner followed by this call round the vertex the setback way — the
/// `b2_three_fillets_and_the_corner_tool_round_the_vertex` pin measures the
/// result against a closed form.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// vertex is not a vertex of the solid; if it is not trihedral — exactly
/// three planes must pass through it (a curved-edged corner is the
/// N-support setback's, still owed — docs/PARITY.md, fillet.edge-blends);
/// or if the corner turns out concave, where a ball adds material instead
/// of shedding it and a tool built from a cut cannot say so. The corner
/// may be oblique: the block is then the hexahedron bounded by the host
/// planes and the three planes through the ball's centre square to the
/// edges.
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

    // The corner's frame comes from the three planes that pass through the
    // vertex's point — not from the vertex's own adjacency, which the very
    // sequence this tool serves destroys: after three fillets the tip is
    // consumed, but the three shrunk planes still contain the corner, and
    // still say exactly which trihedral corner it was. The vertex argument
    // may therefore come from an earlier state of the solid — the sharp
    // box's corner captured before the fillets — and anchors the history.
    let mut normals: Vec<Vector> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        if let ogeom_geom::SurfaceGeometry::Plane(p) = surface {
            let placed = face.transform(model.datums())?;
            let origin = placed.apply(p.plane().frame().origin());
            let normal = placed.apply_vector(p.plane().frame().z().vector());
            if (corner - origin).dot(normal).abs() > tol.confusion() * 100.0 {
                continue;
            }
            // One vote per plane: coplanar trims share it.
            if normals
                .iter()
                .any(|n| n.cross(normal).magnitude() < tol.angular() * 10.0)
            {
                continue;
            }
            normals.push(normal);
        }
    }
    if normals.len() != 3 {
        ogeom_bail!(
            Construction,
            "round_vertex speaks the trihedral corner: exactly three planes \
             must pass through the vertex, found {}; the curved and N-face \
             corners are the setback family's, still owed — docs/PARITY.md, \
             fillet.edge-blends",
            normals.len()
        );
    }
    // The corner's own axes: each edge is the meet of two of the planes,
    // so its direction is their normals' cross product, signed to point
    // from the vertex into the solid. Three planes that do not span are no
    // corner.
    let span = normals[0].cross(normals[1]).dot(normals[2]).abs();
    if span <= tol.angular() * 10.0 {
        ogeom_bail!(
            Construction,
            "the three planes through this vertex do not span a corner"
        );
    }
    let edge_of = |i: usize, j: usize| -> OgeomResult<Direction> {
        Direction::new(normals[i].cross(normals[j]), tol)
    };
    let axes = [edge_of(1, 2)?, edge_of(2, 0)?, edge_of(0, 1)?];
    // Which way along each edge the material lies: the ball's centre — the
    // point a radius in from all three planes — sits inside the solid for
    // exactly one of the eight sign choices on a convex corner. It is
    // probed rather than derived from the faces' orientations because the
    // vertex may already be consumed: after three fillets the tip is gone,
    // but the shrunk planes still say where the material was.
    let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    let centre_for = |dirs: [Vector; 3]| -> Option<Point> {
        // Inward plane normals: each signed toward the diagonal.
        let diagonal = dirs[0] + dirs[1] + dirs[2];
        let m: Vec<Vector> = normals
            .iter()
            .map(|n| if n.dot(diagonal) >= 0.0 { *n } else { -*n })
            .collect();
        // Solve m_k · (c − corner) = radius for c.
        let det = m[0].cross(m[1]).dot(m[2]);
        if det.abs() <= 1e-12 {
            return None;
        }
        let rhs = Vector::new(radius, radius, radius);
        let col = |k: usize| -> f64 {
            let mut a = [m[0], m[1], m[2]];
            for (r, row) in a.iter_mut().enumerate() {
                let v = [row.x, row.y, row.z];
                let mut w = v;
                w[k] = [rhs.x, rhs.y, rhs.z][r];
                *row = Vector::new(w[0], w[1], w[2]);
            }
            a[0].cross(a[1]).dot(a[2]) / det
        };
        Some(corner + Vector::new(col(0), col(1), col(2)))
    };
    let mut inward: Option<([Vector; 3], Point)> = None;
    for signs in 0..8_u8 {
        let cand = [
            axes[0].vector() * if signs & 1 == 0 { 1.0 } else { -1.0 },
            axes[1].vector() * if signs & 2 == 0 { 1.0 } else { -1.0 },
            axes[2].vector() * if signs & 4 == 0 { 1.0 } else { -1.0 },
        ];
        let Some(centre) = centre_for(cand) else {
            continue;
        };
        if boundary.holds(model, centre, tol)? == ogeom_algo::Containment::In {
            if inward.is_some() {
                ogeom_bail!(
                    Construction,
                    "two sign choices probe inside; the corner is not the \
                     simple trihedral this tool speaks"
                );
            }
            inward = Some((cand, centre));
        }
    }
    let Some((inward, far)) = inward else {
        ogeom_bail!(
            Construction,
            "no side of this corner holds material; a concave vertex gains a \
             ball instead of shedding one, and this tool cannot round it"
        );
    };
    let d: Vec<Direction> = inward
        .iter()
        .map(|v| Direction::new(*v, tol))
        .collect::<OgeomResult<_>>()?;

    // The block: the corner bounded by its three host planes and, through
    // the ball's centre, the three planes square to its edges — where each
    // band's circle and the ball's own rim coincide, so the cut ends the
    // band and starts the patch on one curve. On a square corner it is the
    // box of side `radius`; on an oblique one a hexahedron whose planar
    // faces `make_hexahedron` checks. Its corners: the vertex, the foot of
    // each edge on its cutting plane, on each host plane the point square
    // to both of that plane's edges, and the centre itself.
    //
    // The ball on the corner's axes: a pole at one corner of the patch it
    // leaves and its seam meridian out past the block through the first
    // edge — the pole axis is the inward normal of the host plane holding
    // the first two edges, the third edge itself on a square corner.
    //
    // Which edge is first, second and third is the tool's labelling, and
    // the solid it builds is the same for all six. The boolean is not yet
    // indifferent to it: the charts the block's faces and the ball wear
    // decide where a rim is exact and where fitted, where a seam falls
    // against a patch arc, and at an oblique corner five of the six
    // labellings still die in the cut in five different ways. So the tool
    // is offered on each labelling in turn and the first that closes
    // stands — every one of them is the same exact construction — and the
    // corner is refused by name only when none does. A failed attempt's
    // nodes stay in the model unreferenced, under their own operation.
    // The boolean closing all six is owed (docs/PARITY.md,
    // fillet.edge-blends).
    let attempt = |model: &mut Model, order: [usize; 3]| -> OgeomResult<Built> {
        let d = [d[order[0]], d[order[1]], d[order[2]]];
        let along =
            |k: usize| -> Point { corner + d[k].vector() * (far - corner).dot(d[k].vector()) };
        let on_plane = |i: usize, j: usize| -> OgeomResult<Point> {
            let (ei, ej) = (d[i].vector(), d[j].vector());
            let rhs = [(far - corner).dot(ei), (far - corner).dot(ej)];
            let m = [[ei.dot(ei), ej.dot(ei)], [ei.dot(ej), ej.dot(ej)]];
            let det = m[0][0].mul_add(m[1][1], -(m[0][1] * m[1][0]));
            if det.abs() <= 1e-12 {
                ogeom_bail!(Construction, "two edges of this corner are parallel");
            }
            let alpha = rhs[0].mul_add(m[1][1], -(m[0][1] * rhs[1])) / det;
            let beta = m[0][0].mul_add(rhs[1], -(rhs[0] * m[1][0])) / det;
            Ok(corner + ei * alpha + ej * beta)
        };
        let corners = [
            corner,
            along(0),
            on_plane(0, 1)?,
            along(1),
            along(2),
            on_plane(0, 2)?,
            far,
            on_plane(1, 2)?,
        ];
        model.begin_operation();
        let block = ogeom_algo::make_hexahedron(model, corners, tol)?.shape;
        let ball_frame = {
            let plane_normal = d[0].vector().cross(d[1].vector());
            let z = if plane_normal.dot(d[2].vector()) >= 0.0 {
                plane_normal
            } else {
                -plane_normal
            };
            Frame::new(far, Direction::new(z, tol)?, d[0], tol)?
        };
        let ball = ogeom_algo::make_sphere(model, ball_frame, radius, tol)?.shape;
        let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
        let rounded = ogeom_bool::cut(model, solid, &tool.shape, tol)?;
        Ok(Built {
            shape: rounded.shape,
            history: tool.history.then(&rounded.history),
        })
    };
    const LABELLINGS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [1, 0, 2],
        [0, 2, 1],
        [2, 1, 0],
    ];
    let mut outcome: Option<Built> = None;
    let mut last: Option<ogeom_core::OgeomError> = None;
    for order in LABELLINGS {
        match attempt(model, order) {
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
            "the corner tool's cut closed on none of the corner's six \
             labellings; the last said: {}",
            last.map_or_else(String::new, |e| e.to_string())
        );
    };
    let mut built = rounded;
    built.history.modify(vertex, built.shape.clone());
    Ok(built)
}
