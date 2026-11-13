//! P3: each vertex closed by one method, chosen by geometry.
//!
//! Every vertex's hole is bounded by its faces' setback corners and the
//! strips' ends there, consistent by construction: each corner is shared
//! by exactly the two curves that meet at it. The method follows the ADR's
//! survey (docs/adr/0001-blend-architecture.md, §13): a chamfer corner is the
//! facet through its faces' setback corners where they lie on one plane,
//! and a patch through them where they do not; a round corner is a sphere
//! where one ball touches every face, and one n-sided patch elsewhere.

use ogeom_core::{OgeomResult, ogeom_bail};

use super::classify::VertexClass;
use super::request::Section;

/// How one vertex is closed. The extension point for vertex blends: the
/// facet here; the sphere and the n-sided patch with rounds in M2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum VertexMethod {
    /// The polygon through the faces' setback corners round the vertex, in
    /// turn: a facet where they span a plane, nothing where they coincide,
    /// part of a face where they lie in its plane. The assembly and its
    /// guards settle which, and refuse one off a plane by name.
    Facet,
}

/// The method for a vertex of `class` whose blends have `sections`.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for a
/// vertex whose blends this engine does not close yet.
pub(crate) fn method(_class: VertexClass, sections: &[&Section]) -> OgeomResult<VertexMethod> {
    if sections
        .iter()
        .all(|s| matches!(s, Section::Chamfer { .. }))
    {
        return Ok(VertexMethod::Facet);
    }
    ogeom_bail!(
        Construction,
        "this vertex's blends are not closed by the blend engine yet"
    )
}

/// How a vertex where rounds meet is closed, chosen by its geometry (ADR
/// 0001, §13, question 2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum RoundVertex {
    /// No blended edge: the vertex stays.
    Sharp,
    /// One round ending on the third face of a three-edge vertex, that face
    /// cutting the round's edge: the strip's end (its section's arc where
    /// the face stands square to the edge, an ellipse where it leans) lies
    /// in it, and the face takes the end for its corner, losing the corner
    /// of a convex round and gaining the fill of a concave one.
    Flush,
    /// Rounds of one radius at every edge of a vertex of three planes or
    /// more, all convex or all concave, where one ball touches every face
    /// at its setback corner: the spherical polygon between the strips' end
    /// arcs.
    Sphere,
    /// Rounds at every edge of a vertex of three planes, or of more of one
    /// sense, where no one ball touches the faces at their setback corners
    /// (radii differing, convex and concave meeting, or planes no ball
    /// touches together): the strips end through their faces' setback
    /// corners, skew where those stand at different stations, and one
    /// transfinite patch closes the hole tangent to all of them, on one
    /// face for three or four sides and on a face per side for more.
    Patch,
    /// A convex and a concave round beside an edge left sharp at a vertex of
    /// three planes: the strips end through their faces' setback corners,
    /// and one patch closes the corner between them, its fourth and third
    /// sides along the sharp edge's line through the vertex, tangent to the
    /// face whose corner stands on that line's extension below the vertex
    /// and to the other above it, its normal turning at the vertex (ADR
    /// 0001, §13, question 6).
    Signed,
    /// Rounds of one radius at every edge of a convex vertex of four planes
    /// where no one ball touches all four but one touches three at each end
    /// of a ridge (an elongated pyramid's apex): the ball rolls along the
    /// ridge between the two faces either side of it, a sphere at each end
    /// and the cylinder between, as the rolling ball rounds it.
    Ridge,
    /// Two rounds of one sense beside an edge left sharp. Of one radius with
    /// their axes meeting, the strips run on to the mitre plane, where the
    /// two cylinders cross on one ellipse, and meet along it. Otherwise the
    /// larger (its corner on the sharp edge the farther) runs on
    /// to the smaller's other face and ends flush in it, and the smaller ends
    /// on the larger along the curve the two cylinders cross in. Either way
    /// the sharp edge's crease carries on between them, and there is no
    /// vertex face.
    Mitre,
}
