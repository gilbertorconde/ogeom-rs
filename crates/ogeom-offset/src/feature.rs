//! Form features: the named operations a modeller thinks in.
//!
//! A pocket is a prism cut into a solid; a pad is the same prism fused onto
//! it; a rib is a thin pad; a slot is a pocket that runs out of both ends;
//! a revolved feature turns a profile instead of sweeping it. None of these
//! is a new geometric construction — each is a sweep and a boolean, and
//! what makes it a *feature* is that the operation says which it was and
//! carries the profile through the history.
//!
//! Building them here rather than leaving them to the caller is not
//! ceremony. It fixes the two things a caller gets wrong: which way the
//! sweep should run so the tool reaches the material it is meant to reach,
//! and what the history should say afterwards. Both are in one place, and
//! the vocabulary is the one drawings use.

use ogeom_algo::{Built, make_prism, make_revolution};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Axis, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

/// Which way a feature meets the material.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Feature {
    /// The swept tool is added: a pad, a boss, a rib.
    Added,
    /// The swept tool is removed: a pocket, a slot, a groove.
    Removed,
}

/// Sweep `profile` along `vector` and add or remove the result.
///
/// The profile is a face — a wire is not a tool, it is the boundary of one
/// — and the sweep is the ordinary prism, so the feature's walls are ruled
/// exactly as its profile's edges are. A pocket deeper than the material
/// simply cuts through; a pad shorter than nothing is refused.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// profile is not a face, the vector has no length, or the boolean refuses
/// the configuration.
pub fn feature_prism(
    model: &mut Model,
    solid: &Shape,
    profile: &Shape,
    vector: Vector,
    sense: Feature,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "a form feature sweeps a face; a wire is the boundary of one, \
             not a tool"
        );
    }
    let tool = make_prism(model, profile, vector, tol)?;
    applied(model, solid, &tool.shape, profile, sense, tol)
}

/// Turn `profile` about `axis` through `angle` and add or remove the result.
///
/// # Errors
///
/// As [`feature_prism`], plus whatever the revolution refuses — an angle
/// outside `(0, 2π]`, or a profile the axis passes through.
pub fn feature_revol(
    model: &mut Model,
    solid: &Shape,
    profile: &Shape,
    axis: Axis,
    angle: f64,
    sense: Feature,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "a form feature turns a face; a wire is the boundary of one, \
             not a tool"
        );
    }
    let tool = make_revolution(model, profile, axis, angle, tol)?;
    applied(model, solid, &tool.shape, profile, sense, tol)
}

/// A rib: a pad of stated thickness, swept from a profile face's own plane.
///
/// The rib is the profile thickened along `normal` by `thickness` and fused
/// on. It is `feature_prism` with the vector spelled for the case, and it
/// exists because a rib is a thing a drawing names and a caller should not
/// have to spell as a prism every time.
///
/// # Errors
///
/// As [`feature_prism`].
pub fn feature_rib(
    model: &mut Model,
    solid: &Shape,
    profile: &Shape,
    normal: Vector,
    thickness: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !thickness.is_finite() || thickness <= tol.confusion() {
        ogeom_bail!(Construction, "a rib of {thickness} thickness holds nothing");
    }
    let magnitude = normal.magnitude();
    if magnitude <= tol.confusion() {
        ogeom_bail!(Construction, "a rib needs a direction to stand in");
    }
    feature_prism(
        model,
        solid,
        profile,
        normal / magnitude * thickness,
        Feature::Added,
        tol,
    )
}

/// A slot: a pocket swept along a direction and cut clean through.
///
/// The prism runs `depth` each way from the profile, so the tool leaves the
/// material at both ends and the slot is open however the profile sits.
///
/// # Errors
///
/// As [`feature_prism`].
pub fn feature_slot(
    model: &mut Model,
    solid: &Shape,
    profile: &Shape,
    along: Vector,
    depth: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !depth.is_finite() || depth <= tol.confusion() {
        ogeom_bail!(Construction, "a slot of depth {depth} cuts nothing");
    }
    let magnitude = along.magnitude();
    if magnitude <= tol.confusion() {
        ogeom_bail!(Construction, "a slot needs a direction to run in");
    }
    let direction = along / magnitude;
    // Swept from behind the profile to past it: a slot is open at both
    // ends, and a tool that starts *on* the profile leaves a skin.
    let started = ogeom_algo::transformed(
        model,
        profile,
        ogeom_math::Transform::translation(-direction * depth),
    )?;
    feature_prism(
        model,
        solid,
        &started.shape,
        direction * (depth * 2.0),
        Feature::Removed,
        tol,
    )
}

/// The boolean half, with the profile carried into the history.
fn applied(
    model: &mut Model,
    solid: &Shape,
    tool: &Shape,
    profile: &Shape,
    sense: Feature,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut result = match sense {
        Feature::Added => ogeom_bool::fuse(model, solid, tool, tol)?,
        Feature::Removed => ogeom_bool::cut(model, solid, tool, tol)?,
    };
    // What the feature was made from is what a later edit will name, so the
    // profile generates the result rather than vanishing into the tool.
    result.history.generate(profile, result.shape.clone());
    Ok(result)
}
