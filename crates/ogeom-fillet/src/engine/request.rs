//! The request: every blended edge with its section, read as one set.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape};

use crate::chamfer::{Chamfer, Legs, legs};

/// Every blend of one operation. Its order carries no meaning: every phase
/// walks the solid's own elements in their key order.
#[derive(Debug, Clone)]
pub(crate) struct Request {
    pub blends: Vec<Blend>,
    /// Whether a vertex of three chamfers closes by their mitre, the point
    /// their planes share, rather than by the facet through their faces'
    /// setback corners ([`ChamferCorner`](crate::ChamferCorner)).
    pub mitre: bool,
}

/// One edge and the section swept along it.
#[derive(Debug, Clone)]
pub(crate) struct Blend {
    pub edge: Shape,
    pub section: Section,
}

/// A section, the only thing a blend says about itself; the engine's
/// phases past P1 never look inside one. A new kind of section (a round,
/// an asymmetric or conic profile, a law along the edge) is a new variant
/// here and its [`SectionLaw`](super::sections::SectionLaw).
#[derive(Debug, Clone)]
pub(crate) enum Section {
    /// A flat bevel: its distance in from the edge on each of the edge's
    /// two faces, each face named.
    Chamfer { on: [(Shape, f64); 2] },
    /// A round: the rolling ball of this radius.
    Round { radius: f64 },
}

impl Request {
    /// The chamfer specs as a request, each one's legs read on `solid` as
    /// given, or `None` where a spec's seat is not two planes.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// for an edge asked twice or a distance that bevels nothing.
    pub fn from_chamfers(
        model: &Model,
        solid: &Shape,
        specs: &[(Shape, Chamfer)],
        tol: Tolerances,
    ) -> OgeomResult<Option<Self>> {
        let mut blends: Vec<Blend> = Vec::with_capacity(specs.len());
        for (edge, spec) in specs {
            let Ok(Legs::Straight(seat, distances)) = legs(model, solid, edge, spec, tol) else {
                return Ok(None);
            };
            if blends.iter().any(|b| b.edge.is_same(edge)) {
                ogeom_bail!(Construction, "an edge is asked to be chamfered twice");
            }
            for d in distances {
                if !(d.is_finite() && d > tol.confusion()) {
                    ogeom_bail!(Construction, "a chamfer of {d} bevels nothing");
                }
            }
            blends.push(Blend {
                edge: edge.clone(),
                section: Section::Chamfer {
                    on: [
                        (seat.faces[0].clone(), distances[0]),
                        (seat.faces[1].clone(), distances[1]),
                    ],
                },
            });
        }
        Ok(Some(Self {
            blends,
            mitre: false,
        }))
    }

    /// Each edge rounded at its radius, as one request.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// for an edge asked twice or a radius that rounds nothing.
    pub fn from_rounds(rounds: &[(Shape, f64)], tol: Tolerances) -> OgeomResult<Self> {
        let mut blends: Vec<Blend> = Vec::with_capacity(rounds.len());
        for (edge, radius) in rounds {
            if blends.iter().any(|b| b.edge.is_same(edge)) {
                ogeom_bail!(Construction, "an edge is asked to be rounded twice");
            }
            if !(radius.is_finite() && *radius > tol.confusion()) {
                ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
            }
            blends.push(Blend {
                edge: edge.clone(),
                section: Section::Round { radius: *radius },
            });
        }
        Ok(Self {
            blends,
            mitre: false,
        })
    }
}
