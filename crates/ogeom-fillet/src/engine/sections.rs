//! P1: each blend's section, solved against the solid as given.
//!
//! Nothing built here changes the solid another blend sees. A section law
//! answers two questions about its edge: how far in from the edge its
//! contact lies on each of the edge's faces (the setback the faces take in
//! P2), and what fills the strip between its contact lines once P2 has
//! found where the strip ends.

use ogeom_core::{OgeomResult, ogeom_bail};

use super::classify::{Planar, Seat};
use super::id;
use super::request::Section;

/// What a section law gives the engine for one blended edge between two
/// planes. The extension point for every kind of section: chamfers here;
/// rounds, conic and asymmetric profiles and laws along the edge in M2.
pub(crate) trait SectionLaw {
    /// The contact's distance in from the edge on each of the seat's two
    /// faces, in the seat's order.
    ///
    /// # Errors
    ///
    /// Where the section does not seat on the edge.
    fn setbacks(&self, seat: &Seat, planar: &Planar) -> OgeomResult<[f64; 2]>;

    /// What fills the strip between the contact lines, its corners given as
    /// points' indices: on the seat's first face at the edge's start and
    /// end, then on its second face at the start and end.
    ///
    /// # Errors
    ///
    /// Where the strip cannot be built.
    fn strip(&self, seat: &Seat, corners: StripCorners) -> OgeomResult<Strip>;
}

/// A strip's corners, as indices into the result's points.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StripCorners {
    pub first: [usize; 2],
    pub second: [usize; 2],
}

/// What fills a strip.
#[derive(Debug, Clone)]
pub(crate) enum Strip {
    /// The planar ring through the corners, wound out of the solid: a
    /// chamfer's bevel.
    Planar(Vec<usize>),
}

/// The law of `section`.
pub(crate) fn law(section: &Section) -> Box<dyn SectionLaw + '_> {
    match section {
        Section::Chamfer { on } => Box::new(ChamferLaw { on }),
        Section::Round { radius } => Box::new(RoundLaw { radius: *radius }),
        Section::Conic { on, .. } => Box::new(ConicLaw { on }),
    }
}

/// A flat bevel: the plane through its two contact lines.
struct ChamferLaw<'a> {
    on: &'a [(ogeom_topo::Shape, f64); 2],
}

impl SectionLaw for ChamferLaw<'_> {
    fn setbacks(&self, seat: &Seat, planar: &Planar) -> OgeomResult<[f64; 2]> {
        let mut out = [0.0; 2];
        for (k, face) in seat.faces.iter().enumerate() {
            let shape = &planar.faces[*face].shape;
            let Some((_, d)) = self.on.iter().find(|(f, _)| id(f) == id(shape)) else {
                ogeom_bail!(
                    Construction,
                    "a chamfer names a face its edge does not bound"
                );
            };
            out[k] = *d;
        }
        Ok(out)
    }

    fn strip(&self, _seat: &Seat, corners: StripCorners) -> OgeomResult<Strip> {
        // The first face walks the edge from its start to its end; out of
        // the solid the bevel goes round against it.
        Ok(Strip::Planar(vec![
            corners.first[1],
            corners.first[0],
            corners.second[0],
            corners.second[1],
        ]))
    }
}

/// A round: the rolling ball of `radius`, touching each face where its
/// contact line stands `r cot(theta / 2)` in from the edge, `theta` the
/// dihedral's opening on the ball's side (the material's at a convex edge,
/// the open side's at a concave one). Its strip is curved, built by
/// [`super::rounds`].
struct RoundLaw {
    radius: f64,
}

impl SectionLaw for RoundLaw {
    fn setbacks(&self, seat: &Seat, _planar: &Planar) -> OgeomResult<[f64; 2]> {
        let theta = if seat.convex {
            seat.opening
        } else {
            core::f64::consts::TAU - seat.opening
        };
        if !(theta > 0.0 && theta < core::f64::consts::PI) {
            ogeom_bail!(
                Construction,
                "a round's faces meet too flat or too sharply to seat it"
            );
        }
        let d = self.radius / (theta * 0.5).tan();
        Ok([d, d])
    }

    fn strip(&self, _seat: &Seat, _corners: StripCorners) -> OgeomResult<Strip> {
        ogeom_bail!(
            Invariant,
            "a round's strip is curved: the round sheet builds it, not the planar rings"
        )
    }
}

/// A conic round's: each face named with its distance, as a chamfer's; its
/// strip is the rounds' sheet's to build.
struct ConicLaw<'a> {
    on: &'a [(ogeom_topo::Shape, f64); 2],
}

impl SectionLaw for ConicLaw<'_> {
    fn setbacks(&self, seat: &Seat, planar: &Planar) -> OgeomResult<[f64; 2]> {
        ChamferLaw { on: self.on }.setbacks(seat, planar)
    }

    fn strip(&self, _seat: &Seat, _corners: StripCorners) -> OgeomResult<Strip> {
        ogeom_bail!(
            Construction,
            "a conic round's strip is built by the rounds' sheet"
        )
    }
}
