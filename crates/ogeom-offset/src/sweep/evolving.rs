//! Pipe shells whose section scales and twists along the spine.
//!
//! The section is the profile placed in the law's frame at each point of
//! the spine, as for [`make_pipe_shell_with`](super::make_pipe_shell_with),
//! and before that scaled about the spine and turned about it, both in
//! proportion to the share of the spine's length run there.
//!
//! Scaled without a twist down a straight spine, every profile point runs
//! along a straight line, at a constant speed, to its scaled image at the
//! end: the solid is the ruled loft between the profile and that image,
//! exact (a pyramid's frustum over a polygon, a cone's over a circle, the
//! pyramid or cone itself where the scale reaches zero). Everywhere else
//! the walls are skinned through the placed sections and held to the
//! tolerance between them against the law itself.

use ogeom_algo::{Built, History, curve_length, find_plane};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve;
use ogeom_math::{Axis, Direction, Point, Transform};
use ogeom_topo::{Filter, Model, Shape, ShapeType, VertexData, explore};

use super::{
    PipeCorners, PipeLaw, WalkCurve, densified, evenly_by, has_sharp_corner, law_loft, make_loft,
    make_pipe_shell_with, moved_ring, profile_square_at_start, shell_stations, spine_curve_of,
    straight_legs,
};

/// How a pipe's section scales and twists along its spine, each in
/// proportion to the length run along the spine: none at its start, all
/// of it at its end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PipeEvolution {
    /// The section's scale at the spine's end, about the spine; it is one
    /// at the start. Zero closes the pipe on the spine's end point.
    pub scale: f64,
    /// The angle in radians the section has turned through about the
    /// spine by its end, counter-clockwise about the direction of travel.
    pub twist: f64,
}

impl PipeEvolution {
    /// The scale and turn a share `lambda` of the way along, in the start
    /// frame's coordinates: about its origin, the spine, and its `z`, the
    /// spine's tangent.
    pub(super) fn local(&self, lambda: f64, tol: Tolerances) -> OgeomResult<Transform> {
        let turn = Transform::rotation(
            Axis {
                location: Point::ORIGIN,
                direction: Direction::Z,
            },
            self.twist * lambda,
        );
        let scale = 1.0 + (self.scale - 1.0) * lambda;
        Ok(turn * Transform::scaling(Point::ORIGIN, scale, tol)?)
    }
}

/// The length of a spine's edges and of the run up to each.
pub(super) struct SpineLengths {
    before: Vec<f64>,
    total: f64,
}

impl SpineLengths {
    /// The lengths along `curves`, a spine's edges in order.
    pub(super) fn of(curves: &[WalkCurve], tol: Tolerances) -> OgeomResult<Self> {
        let mut before = Vec::with_capacity(curves.len());
        let mut total = 0.0;
        for (curve, range, _) in curves {
            before.push(total);
            total += curve_length(curve, *range, tol)?;
        }
        if total <= tol.confusion() {
            ogeom_bail!(Construction, "the spine has no length");
        }
        Ok(Self { before, total })
    }

    /// The share of the spine's length run to parameter `t` on edge `edge`.
    pub(super) fn fraction(
        &self,
        curves: &[WalkCurve],
        edge: usize,
        t: f64,
        tol: Tolerances,
    ) -> OgeomResult<f64> {
        let (curve, range, reversed) = &curves[edge];
        let run = if *reversed {
            curve_length(curve, (t, range.1), tol)?
        } else {
            curve_length(curve, (range.0, t), tol)?
        };
        Ok(((self.before[edge] + run) / self.total).clamp(0.0, 1.0))
    }
}

/// Sweep a planar profile along a spine under a frame law, the section
/// scaling and twisting as it goes.
///
/// The profile (a closed wire, or a face whose holes ride along) stands
/// square to the spine at its start, as for a pipe shell. At each point of
/// the spine the section is the profile scaled about the spine by
/// `1 + (scale - 1) λ` and turned about it by `twist λ`, `λ` the share of
/// the spine's length run there, then carried by `law`. With no scale and
/// no twist this is [`make_pipe_shell_with`] with mitred corners.
///
/// Scaled without a twist down one straight edge under the
/// rotation-minimizing, binormal or fixed law (none of which turns the
/// section along a line), the solid is the ruled loft from the profile to
/// its scaled image, exact; a scale of zero closes it on the spine's end,
/// the exact pyramid or cone. Otherwise each loop of the profile is
/// skinned through its placed sections and held to `tolerance` between
/// them against the law, along an open spine with no sharp corner.
///
/// Each hole sweeps its own solid under the same law, cut from the outer
/// one.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `scale` is negative or `twist` not finite; if the profile is not a
/// closed planar wire or face standing square to the spine at its start;
/// if a skinned sweep is asked along a closed spine or round a sharp
/// corner, of the fixed law, or to a scale of zero; and as
/// [`make_pipe_shell_with`] or the law refuses the spine.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the skin
/// cannot reach `tolerance`.
pub fn make_pipe_shell_evolving(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    evolution: PipeEvolution,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !evolution.scale.is_finite() || evolution.scale < 0.0 {
        ogeom_bail!(
            Construction,
            "a section scales to {} at the spine's end; a scale is naught or more",
            evolution.scale
        );
    }
    if !evolution.twist.is_finite() {
        ogeom_bail!(
            Construction,
            "a twist of {} is not an angle",
            evolution.twist
        );
    }
    if evolution.scale == 1.0 && evolution.twist == 0.0 {
        return make_pipe_shell_with(
            model,
            profile,
            spine,
            law,
            PipeCorners::Mitre,
            tolerance,
            tol,
        );
    }
    let rings: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => explore(model, profile, Filter::OfType(ShapeType::Wire))?,
        ShapeType::Wire => vec![profile.clone()],
        other => ogeom_bail!(
            Construction,
            "a pipe shell sweeps a planar wire or face, not a {other:?}"
        ),
    };
    let Some(plane) = find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let stations = shell_stations(model, spine, tol)?;
    profile_square_at_start(model, &plane, spine, stations[0].at, tol)?;
    if rings.len() > 1 && evolution.scale == 0.0 {
        ogeom_bail!(
            Construction,
            "a profile with holes cannot close on a point; its holes would meet"
        );
    }
    let mut solids = Vec::with_capacity(rings.len());
    for ring in &rings {
        solids.push(evolving_ring(
            model, ring, spine, law, evolution, tolerance, tol,
        )?);
    }
    let mut shape = solids[0].clone();
    for hole in &solids[1..] {
        shape = ogeom_bool::cut(model, &shape, hole, tol)?.shape;
    }
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// One closed ring of the profile swept under the evolving law, as a solid.
fn evolving_ring(
    model: &mut Model,
    ring: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    evolution: PipeEvolution,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let level = matches!(
        law,
        PipeLaw::RotationMinimizing | PipeLaw::Binormal(_) | PipeLaw::Fixed
    );
    if evolution.twist == 0.0
        && level
        && let Some(legs) = straight_legs(model, spine, tol)?
        && legs.len() == 1
    {
        let (start, run) = legs[0];
        let end = if evolution.scale == 0.0 {
            model.add_vertex(VertexData::new(start + run))
        } else {
            let motion =
                Transform::translation(run) * Transform::scaling(start, evolution.scale, tol)?;
            moved_ring(model, ring, &motion, tol)?
        };
        // A polygon or a circle lofts to its image as planes or a cone;
        // any other ring as the ruled walls between matched edges.
        let shape = if evolution.scale == 0.0 || analytic_ring(model, ring, start, tol)? {
            make_loft(model, ring, &end, tol)?.shape
        } else {
            crate::loft::make_loft_solid(model, &[ring.clone(), end], &[], true, tol)?.shape
        };
        return Ok(shape);
    }
    if matches!(law, PipeLaw::Fixed) {
        ogeom_bail!(
            Construction,
            "a fixed section scales or twists only down one straight edge, unturned"
        );
    }
    if evolution.scale == 0.0 {
        ogeom_bail!(
            Construction,
            "a section shrinks to a point only down one straight edge, unturned"
        );
    }
    let stations = shell_stations(model, spine, tol)?;
    let ring_spine =
        stations[0].at.distance(stations[stations.len() - 1].at) <= tol.confusion() * 10.0;
    if ring_spine || has_sharp_corner(model, spine, tol)? {
        ogeom_bail!(
            Construction,
            "a scaling or twisting section sweeps an open spine with no sharp corner"
        );
    }
    // As many stations as the law's turning or the twist asks, each edge
    // its share by length.
    let turned = densified(model, spine, stations.clone(), law, tolerance, tol)?.len();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let twisted = (evolution.twist.abs() / super::MOST_TURN).ceil() as usize + 1;
    let wanted = turned.max(twisted);
    let lengths = edge_lengths(model, spine, tol)?;
    let total: f64 = lengths.iter().sum();
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let share = |e: usize| -> u32 {
        let length = lengths.get(e).copied().unwrap_or(0.0);
        ((wanted as f64 * length / total).ceil() as u32).max(2)
    };
    let stations = evenly_by(model, spine, stations, &share, false, tol)?;
    Ok(law_loft(
        model,
        ring,
        spine,
        &stations,
        law,
        Some(evolution),
        tolerance,
        tol,
    )?
    .shape)
}

/// Whether a ring is all straight edges, or one circle about `axis_at`:
/// what [`make_loft`] lofts to its image scaled about `axis_at` exactly.
fn analytic_ring(
    model: &Model,
    ring: &Shape,
    axis_at: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let edges = explore(model, ring, Filter::OfType(ShapeType::Edge))?;
    let mut kinds = Vec::with_capacity(edges.len());
    for edge in &edges {
        kinds.push(spine_curve_of(model, edge)?.0);
    }
    Ok(kinds.iter().all(|c| matches!(c, Curve::Line(_)))
        || (kinds.len() == 1
            && matches!(&kinds[0], Curve::Circle(c)
                if c.circle().centre().distance(axis_at) <= tol.confusion())))
}

/// Each spine edge's length, in the spine's order.
fn edge_lengths(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<f64>> {
    let edges: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        _ => model.ordered_children_of(spine)?,
    };
    let mut out = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range) = spine_curve_of(model, edge)?;
        out.push(curve_length(&curve, range, tol)?);
    }
    Ok(out)
}
