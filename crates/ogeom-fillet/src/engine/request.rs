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
    /// Edges read as pieces: where a boolean split one edge of the form and
    /// its pieces are asked apart, each piece is an edge of its own, the
    /// edge read from the solid split at the vertices between them.
    pub splits: Vec<Split>,
    /// Vertices closed as setback corners (ADR 0001, §13 question 2), each
    /// with the distance along each of its edges at which that edge's strip
    /// stops: only where asked, never by default.
    pub corners: Vec<(Shape, crate::Setbacks)>,
}

/// One edge of the solid read as the pieces another reading of the part
/// holds: each piece's own edge and the vertices between them (shapes not
/// in the solid, named only to key the pieces and their history).
#[derive(Debug, Clone)]
pub(crate) struct Split {
    pub edge: Shape,
    pub pieces: Vec<Shape>,
    pub joints: Vec<Shape>,
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
    /// A round meeting each named face its distance from the edge, its
    /// section `profile` between the two contacts (an ellipse, or a conic
    /// of its own `rho`).
    Conic {
        on: [(Shape, f64); 2],
        profile: crate::asym::RoundProfile,
    },
    /// A round whose ball's radius follows `law` along the edge, its
    /// fraction nought at the vertex `start`.
    Varying {
        law: crate::law::RadiusLaw,
        start: Shape,
    },
}

impl Section {
    /// Whether the section is a round's (curved, built by the rounds'
    /// sheet) rather than a bevel's.
    pub(crate) const fn is_round(&self) -> bool {
        matches!(
            self,
            Self::Round { .. } | Self::Conic { .. } | Self::Varying { .. }
        )
    }
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
            splits: Vec::new(),
            corners: Vec::new(),
        }))
    }

    /// One edge rounded meeting `a.face` `a.on_face` from it and its other
    /// face `a.on_other`, its section `a.profile`, as a request of one;
    /// `None` where the edge is not straight between two planes or the face
    /// named is not one of them.
    ///
    /// # Errors
    ///
    /// As reading the edge's seat.
    pub fn from_round_distances(
        model: &Model,
        solid: &Shape,
        a: &crate::asym::Asymmetric,
        tol: Tolerances,
    ) -> OgeomResult<Option<Self>> {
        let Ok(seat) = crate::support::planar_seat(model, solid, &a.edge, tol) else {
            return Ok(None);
        };
        let other = if seat.faces[0].is_same(&a.face) {
            seat.faces[1].clone()
        } else if seat.faces[1].is_same(&a.face) {
            seat.faces[0].clone()
        } else {
            return Ok(None);
        };
        Ok(Some(Self {
            blends: vec![Blend {
                edge: a.edge.clone(),
                section: Section::Conic {
                    on: [(a.face.clone(), a.on_face), (other, a.on_other)],
                    profile: a.profile,
                },
            }],
            mitre: false,
            splits: Vec::new(),
            corners: Vec::new(),
        }))
    }

    /// Rounds of one radius and rounds of their own distances, as one
    /// request: each distance round's edge straight between two planes,
    /// `face` one of them; `None` where one is not.
    ///
    /// # Errors
    ///
    /// As [`Request::from_rounds`].
    pub fn from_mixed_rounds(
        model: &Model,
        solid: &Shape,
        rounds: &[(Shape, f64)],
        distances: &[crate::asym::Asymmetric],
        tol: Tolerances,
    ) -> OgeomResult<Option<Self>> {
        let mut request = Self::from_rounds(rounds, tol)?;
        for a in distances {
            if request.blends.iter().any(|b| b.edge.is_same(&a.edge)) {
                ogeom_bail!(Construction, "an edge is asked to be rounded twice");
            }
            let Some(one) = Self::from_round_distances(model, solid, a, tol)? else {
                return Ok(None);
            };
            request.blends.extend(one.blends);
        }
        Ok(Some(request))
    }

    /// One straight edge rounded at a radius following `law` from its first
    /// vertex, as a request of one; `None` for a law not linear along it.
    ///
    /// # Errors
    ///
    /// As reading the edge's ends.
    pub fn from_varying(
        model: &Model,
        edge: &Shape,
        law: &crate::law::RadiusLaw,
    ) -> OgeomResult<Option<Self>> {
        // A law curving along the edge bends the strip's contacts, which the
        // sheet runs straight: that is the current path's for now.
        if !law.is_linear() {
            return Ok(None);
        }
        let Some((start, _)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(None);
        };
        Ok(Some(Self {
            blends: vec![Blend {
                edge: edge.clone(),
                section: Section::Varying {
                    law: law.clone(),
                    start,
                },
            }],
            mitre: false,
            splits: Vec::new(),
            corners: Vec::new(),
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
            splits: Vec::new(),
            corners: Vec::new(),
        })
    }
}
