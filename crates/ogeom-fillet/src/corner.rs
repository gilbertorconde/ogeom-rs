//! Rounding a vertex: the ball-and-block tool at a trihedral corner.
//!
//! *Elsewhere:* the vertex blend of `ChFi3d`'s setback family.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Direction, Frame, Vector};
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
/// three edges must meet it; if the three edges do not leave it mutually
/// orthogonal and straight, which is the corner this tool speaks (the
/// oblique and curved-edged corners are the N-support setback's, still
/// owed — docs/PARITY.md, fillet.edge-blends); or if the corner turns out
/// concave, where a ball adds material instead of shedding it and a tool
/// built from a cut cannot say so.
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
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        if normals[i].dot(normals[j]).abs()
            > normals[i].magnitude() * normals[j].magnitude() * tol.angular() * 10.0
        {
            ogeom_bail!(
                Construction,
                "the faces meet this vertex obliquely; the oblique corner is \
                 the N-support setback's, still owed — docs/PARITY.md, \
                 fillet.edge-blends"
            );
        }
    }
    let n: Vec<Direction> = normals
        .iter()
        .map(|v| Direction::new(*v, tol))
        .collect::<OgeomResult<_>>()?;

    // The block spans *into* the solid: each plane normal, signed so a probe
    // just inside the would-be block lands inside the solid. The boundary is
    // asked once per sign choice; the corner that answers to none of them is
    // concave or stranger, and this construction cannot round it.
    let boundary = ogeom_algo::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)?;
    // The probe's stand-off is chosen for the corner this tool follows: on a
    // solid whose edges are already filleted at `radius`, material near the
    // corner survives only past the fillet prisms (per-axis offset over
    // r/√2) and within the coming ball's reach (under r). 0.85r sits in
    // that band at every radius, and trivially inside a sharp corner.
    let probe = radius * 0.85;
    let mut inward: Option<[Vector; 3]> = None;
    for signs in 0..8_u8 {
        let cand = [
            n[0].vector() * if signs & 1 == 0 { 1.0 } else { -1.0 },
            n[1].vector() * if signs & 2 == 0 { 1.0 } else { -1.0 },
            n[2].vector() * if signs & 4 == 0 { 1.0 } else { -1.0 },
        ];
        let at = corner + (cand[0] + cand[1] + cand[2]) * probe;
        if boundary.holds(model, at, tol)? == ogeom_algo::Containment::In {
            if inward.is_some() {
                ogeom_bail!(
                    Construction,
                    "two sign choices probe inside; the corner is not the \
                     simple trihedral this tool speaks"
                );
            }
            inward = Some(cand);
        }
    }
    let Some(inward) = inward else {
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

    // The corner block stands at the far corner — the vertex walked in by
    // `radius` along all three edges — and spans back toward the vertex, so
    // its frame's axes are the *negated* inward directions. Negating a
    // right-handed triple as ordered gives a left-handed one; swapping two
    // axes rights it, so the pair feeding the frame is chosen by the
    // triple's own handedness.
    let (first, second) = if d[0].vector().cross(d[1].vector()).dot(d[2].vector()) >= 0.0 {
        (d[1], d[0])
    } else {
        (d[0], d[1])
    };
    let far = corner + (d[0].vector() + d[1].vector() + d[2].vector()) * radius;
    let frame = Frame::new(far, d[2].reversed(), first.reversed(), tol)?;
    debug_assert!(
        frame.y().vector().dot(second.reversed().vector()).abs() > 0.99,
        "the frame's derived axis is the remaining edge"
    );
    let block = ogeom_algo::make_box(model, frame, (radius, radius, radius), tol)?.shape;
    let ball = ogeom_algo::make_sphere(model, frame, radius, tol)?.shape;
    let tool = ogeom_bool::cut(model, &block, &ball, tol)?;
    let rounded = ogeom_bool::cut(model, solid, &tool.shape, tol)?;
    let mut built = Built {
        shape: rounded.shape,
        history: tool.history.then(&rounded.history),
    };
    built.history.modify(vertex, built.shape.clone());
    Ok(built)
}
