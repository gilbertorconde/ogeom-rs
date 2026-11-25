//! Exact section sweeps closed by their rails, crease and end sections.

use crate::RoundProfile;
use crate::law::Profile;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{SheetEdge, SheetFace, build_faces};
use super::fallback::{Sense, Tool};
use crate::round::section::PlanarSection;

/// Apply each tool on its seated faces, retaining only newly created
/// faces as the edge's generated geometry.
pub(super) fn applied(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    tools: &[Tool],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut built = Built::from_nothing(solid.clone());
    for tool in tools {
        let step = match tool.sense {
            Sense::Cut => ogeom_bool::cut_seated(model, &built.shape, &tool.solid, tol),
            Sense::Fuse => ogeom_bool::fuse_seated(model, &built.shape, &tool.solid, tol),
        }?;
        let inherited: Vec<_> = explore_unique(model, &built.shape, ShapeType::Face)?
            .iter()
            .flat_map(|f| step.history.trace(f).to_vec())
            .collect();
        let kept = explore_unique(model, &step.shape, ShapeType::Face)?;
        let mut history = built.history.then(&step.history);
        for (face, owner) in &tool.faces {
            for image in step.history.trace(face) {
                if kept.iter().any(|f| f.is_same(image))
                    && !inherited.iter().any(|f| f.is_same(image))
                {
                    history.generate(owner, image.clone());
                }
            }
        }
        built = Built::new(step.shape, history);
    }
    for edge in edges {
        built.history.delete(edge);
    }
    super::guards::solid_fits(model, &built.shape, tol)?;
    Ok(Some(built))
}

/// A straight planar round under a radius law, including a narrow strip
/// taken wholly at every radius. The law is read from the stored curve's
/// start, independently of the edge occurrence's orientation.
pub(crate) fn linear_round(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    law: &Profile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let consumed = |radius| -> OgeomResult<Option<crate::support::Seat>> {
        match crate::support::ball_fits(model, solid, edge, radius, false, tol) {
            Ok(_) => Ok(None),
            Err(ogeom_core::OgeomError::Construction(_)) => {
                crate::engine::consumed_seat::consumed_seat(model, solid, edge, radius, tol)
            }
            Err(e) => Err(e),
        }
    };
    let seat = match (consumed(law.least())?, consumed(law.largest())?) {
        (Some(_), Some(seat)) => seat,
        (None, None) => {
            let fractions: Vec<f64> = (1..32).map(|k| f64::from(k) / 32.0).collect();
            crate::support::ball_fits_by(
                model,
                solid,
                edge,
                &|s| law.at(s).0,
                &fractions,
                false,
                tol,
            )?;
            crate::support::ball_fits(model, solid, edge, law.largest(), false, tol)?;
            crate::support::planar_seat(model, solid, edge, tol)?
        }
        _ => ogeom_bail!(
            Construction,
            "a variable radius whose ball sits on a narrow face beside the edge at some radii of the law and takes it off whole at others is still owed: the constant round changes construction between them; see docs/CAPABILITIES.md, fillet.variable"
        ),
    };
    let unit = PlanarSection::new(&seat, &crate::round::section::RoundSpec::Circular(1.0), tol)?;
    let distance = law.scaled_by(unit.distances(0.0)[0], false);
    let section = PlanarSection::new(
        &seat,
        &crate::round::section::RoundSpec::Varying {
            laws: [distance.clone(), distance],
            profile: RoundProfile::Elliptical,
        },
        tol,
    )?;
    let tool = sweep(model, edge, &section, tol)?;
    applied(model, solid, std::slice::from_ref(edge), &[tool], tol)?.ok_or_else(|| {
        ogeom_core::ogeom_err!(
            Construction,
            "a varying round's tool does not close on its solid"
        )
    })
}

/// A straight round whose two distances vary, with its contacts checked
/// on the original faces; a concave one ending in the faces at its ends.
#[allow(clippy::too_many_arguments, reason = "one edge's two distance laws")]
pub(crate) fn varying_distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    laws: &[Profile; 2],
    profile: RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (seat, named) = crate::asym::straight_seat(model, solid, edge, face, tol)?;
    let laws = if named == 0 {
        laws.clone()
    } else {
        [laws[1].clone(), laws[0].clone()]
    };
    let spec = crate::round::section::RoundSpec::Varying { laws, profile };
    let section = PlanarSection::new(&seat, &spec, tol)?;
    // Short of the ends, where a face overhanging the corner may cut a
    // contact off by design.
    for k in 1..16 {
        let s = section.length() * f64::from(k) / 16.0;
        for (i, face) in seat.faces.iter().enumerate() {
            if !crate::support::on_face_or_boundary(model, solid, face, section.contact(i, s), tol)?
            {
                ogeom_bail!(
                    Construction,
                    "an asymmetric round sets back {} along a face of the edge, at {:.3} of it, past that face's side; see docs/CAPABILITIES.md, fillet.asymmetric-rounds",
                    section.distances(s)[i],
                    s / section.length()
                );
            }
        }
    }
    let tool =
        super::conic_tools::extended_sweep(model, solid, edge, &seat, &spec, [true; 2], tol)?;
    applied(model, solid, std::slice::from_ref(edge), &[tool], tol)?.ok_or_else(|| {
        ogeom_core::ogeom_err!(
            Construction,
            "a varying distance round's tool does not close on its solid"
        )
    })
}

/// Close a section's sweep over its seat into an engine tool. Every edge
/// is shared by its two faces; the planar legs and caps use the section's
/// own boundary curves.
pub(super) fn sweep(
    model: &mut Model,
    owner: &Shape,
    section: &PlanarSection,
    tol: Tolerances,
) -> OgeomResult<Tool> {
    let length = section.length();
    let [a, b] = section.order();
    let points = [
        section.crease(0.0),
        section.crease(length),
        section.contact(a, 0.0),
        section.contact(a, length),
        section.contact(b, 0.0),
        section.contact(b, length),
    ];
    let mut edges: Vec<SheetEdge> = [(0, 1), (0, 2), (1, 3), (0, 4), (1, 5)]
        .into_iter()
        .map(|ends| SheetEdge { ends, curve: None })
        .collect();
    for (ends, boundary) in [
        ((2, 3), section.contact_curve(a, (0.0, length), tol)?),
        ((4, 5), section.contact_curve(b, (0.0, length), tol)?),
        ((2, 4), section.arc_at(0.0, tol)?),
        ((3, 5), section.arc_at(length, tol)?),
    ] {
        edges.push(SheetEdge {
            ends,
            curve: Some(boundary),
        });
    }
    let rings = [
        vec![(0, true), (2, true), (5, false), (1, false)],
        vec![(0, false), (3, true), (6, true), (4, false)],
        vec![(1, true), (7, true), (3, false)],
        vec![(2, false), (4, true), (8, false)],
        vec![(5, true), (8, true), (6, false), (7, false)],
    ];
    let surface = section.surface((0.0, length), tol)?;
    // The strip's chart can run either way across the crease. Its outward
    // normal points from the crease through the middle of its arc.
    let (arc, range) = section.arc_at(length * 0.5, tol)?;
    let middle = arc.point_at(f64::midpoint(range.0, range.1), tol)?;
    let projected = ogeom_algo::project_on_surface(&surface, middle, 16, tol)?;
    let (du, dv) = surface.d1_at(projected.parameters.0, projected.parameters.1, tol)?;
    let reverse_strip = du.cross(dv).dot(middle - section.crease(length * 0.5)) < 0.0;
    let reverse = (points[2] - points[0])
        .cross(points[4] - points[0])
        .dot(section.seat().along)
        > 0.0;
    let sheet: Vec<SheetFace> = rings
        .into_iter()
        .enumerate()
        .map(|(i, mut ring)| {
            // The surface sense and the adjoining planes' senses can differ;
            // their final walks must still oppose along every shared edge.
            if i == 4 && reverse_strip != reverse {
                ring.reverse();
                for (_, forward) in &mut ring {
                    *forward = !*forward;
                }
            }
            SheetFace {
                surface: (i == 4).then(|| surface.clone()),
                ring,
                holes: Vec::new(),
            }
        })
        .collect();
    let mut faces = build_faces(model, &points, &edges, &sheet, tol)?;
    if reverse_strip {
        faces[4] = faces[4].reversed();
    }
    if reverse {
        for face in &mut faces[..4] {
            *face = face.reversed();
        }
    }
    let shell = model.add_shell(&faces)?;
    let solid = model.add_solid(&[shell])?;
    let faces = explore_unique(model, &solid, ShapeType::Face)?
        .into_iter()
        .map(|f| (f, owner.clone()))
        .collect();
    Ok(Tool {
        solid,
        sense: if section.seat().convex {
            Sense::Cut
        } else {
            Sense::Fuse
        },
        faces,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::round::section::RoundSpec;
    use ogeom_math::{Frame, Point};

    #[test]
    fn a_section_sweep_closes_with_consistent_outward_faces() -> OgeomResult<()> {
        let tol = Tolerances::default();
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (12.0, 12.0, 12.0), tol)?.shape;
        let mut selected = None;
        for edge in explore_unique(&model, &solid, ShapeType::Edge)? {
            let b = ogeom_algo::shape_bounds(&model, &edge, tol)?;
            if b.centre()
                .is_some_and(|p| p.distance(Point::new(0.0, 0.0, 6.0)) <= tol.confusion())
            {
                selected = Some(edge);
                break;
            }
        }
        let Some(edge) = selected else {
            ogeom_bail!(Construction, "the box has no edge along the Z axis");
        };
        let seat = crate::support::planar_seat(&model, &solid, &edge, tol)?;
        for spec in [
            RoundSpec::Circular(1.0),
            RoundSpec::Conic {
                d: [1.0, 1.0],
                profile: RoundProfile::Elliptical,
            },
        ] {
            let section = PlanarSection::new(&seat, &spec, tol)?;
            let tool = sweep(&mut model, &edge, &section, tol)?;
            let found = ogeom_algo::check(&model, &tool.solid, tol)?;
            assert!(found.is_valid(), "{spec:?}: {found}");
            let volume = ogeom_algo::volume_properties(
                &model,
                &tool.solid,
                ogeom_mesh::Deflection::default(),
                tol,
            )?
            .mass;
            let expected = 12.0 * (1.0 - std::f64::consts::FRAC_PI_4);
            // Both sections bound the unit circular spandrel, swept twelve mm.
            assert!(
                (volume - expected).abs() < 1e-6,
                "{spec:?}: {volume}, expected {expected}"
            );
        }
        Ok(())
    }
}
