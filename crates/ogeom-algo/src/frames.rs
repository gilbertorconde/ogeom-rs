//! Rotation-minimizing frames along a curve or a wire.
//!
//! A frame that rides a path should turn only as much as the path makes
//! it: the tangent turns with the path, and the normal turns with the
//! tangent and no more. That frame has no twist about the tangent, it is
//! defined where the curvature vanishes (on a straight run, through an
//! inflection), and unlike the Frenet frame it never flips.
//!
//! # How it is carried
//!
//! By double reflection: from one station to the next the frame is
//! reflected in the plane bisecting the chord between them, then in the
//! plane bisecting the reflected tangent and the next tangent. Two
//! reflections make a rotation, so the frame stays orthonormal and
//! right-handed; the result is exact on lines, circular arcs and every
//! planar curve, and fourth-order accurate elsewhere in the step.
//!
//! Stations stand at least eight to an edge and as many more as hold each
//! step to a sixty-fourth of a turn of the tangent, the turning measured
//! over seventeen samples of the edge. A helix of three turns carried
//! that way twists within about `1e-9` rad of its closed form; a cubic
//! spline knotted every few degrees of its turning, whose third
//! derivative jumps at each knot, within about `1e-5` rad. Between
//! stations the frame is one more double-reflection step from the station
//! behind, so a frame asked for at a station is the station's own.
//!
//! These stations and this carrying are the ones the pipe sweeps use, so a
//! pattern placed along a path with these frames agrees with a pipe swept
//! along the same path.
//!
//! # Corners, closed paths and the first frame
//!
//! - A sharp corner of a wire turns the frame about the corner's own axis
//!   (the cross product of the two tangents) by the corner's angle, the
//!   least rotation taking one tangent to the other; it is applied as a
//!   reflection in the plane halving the corner, which agrees on every
//!   vector square to the incoming tangent.
//! - A closed path brings the frame home turned about its start tangent by
//!   the path's holonomy. The frames are left as carried, and the turn is
//!   reported in [`CurveFrames::closing_twist`]; [`CurveFrames::closed_up`]
//!   pays it off along the path in proportion to length.
//! - The first normal is the caller's, made square to the start tangent,
//!   or by default the world `z` made square to it, or the world `x` where
//!   the tangent runs within thirty degrees of `z`.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _};
use ogeom_math::{Direction, Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType};

use crate::build::placed_edge_curve;
use crate::length::{curve_length, parameter_at_length};

/// One station a frame is carried through: where the path is, which way it
/// runs, and the edge and parameter it stands on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathStation {
    /// The point on the path.
    pub at: Point,
    /// The unit tangent, in the direction of travel.
    pub tangent: Vector,
    /// The edge this station stands on, by position along the path.
    pub edge: usize,
    /// The station's parameter on that edge's curve.
    pub t: f64,
}

/// The edges of a path (an edge or a wire), in the order it walks them.
fn path_edges(model: &Model, path: &Shape) -> OgeomResult<Vec<Shape>> {
    let edges = match model.kind_of(path)? {
        ShapeType::Edge => vec![path.clone()],
        ShapeType::Wire => model.ordered_children_of(path)?,
        other => ogeom_bail!(
            Construction,
            "a path runs along an edge or a wire, not a {other:?}"
        ),
    };
    if edges.is_empty() {
        ogeom_bail!(Construction, "the path has no edge to run along");
    }
    Ok(edges)
}

/// An edge's curve where it stands, its range, and whether the path walks
/// it backwards.
fn walked_curve(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64), bool)> {
    let Some((curve, range)) = placed_edge_curve(model, edge, tol)? else {
        ogeom_bail!(Construction, "a path edge has no curve");
    };
    Ok((curve, range, edge.orientation() == Orientation::Reversed))
}

/// The stations a frame is carried through along `path` (an edge or a
/// wire): each edge sampled evenly in its parameter, at least eight steps
/// and enough more that no step turns the tangent by more than a
/// sixty-fourth of a turn, the turning read from seventeen samples.
///
/// Where two edges meet the walk holds two stations on one point, the end
/// of one edge and the start of the next: heading different ways, that
/// pair is a sharp corner.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the path is neither an edge nor a wire, has an edge with no curve,
/// stands still somewhere (a zero tangent), or collapses to a point.
pub fn path_stations(
    model: &Model,
    path: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<PathStation>> {
    let edges = path_edges(model, path)?;
    let mut stations: Vec<PathStation> = Vec::new();
    for (ei, edge) in edges.iter().enumerate() {
        let (curve, range, reversed) = walked_curve(model, edge, tol)?;
        let turning = {
            let mut sum = 0.0_f64;
            let mut last: Option<Vector> = None;
            for i in 0..=16 {
                let t = range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
                let d = curve.d1_at(t, tol)?;
                let m = d.magnitude();
                if m <= tol.confusion() {
                    continue;
                }
                let u = d / m;
                if let Some(prev) = last {
                    sum += prev.dot(u).clamp(-1.0, 1.0).acos();
                }
                last = Some(u);
            }
            sum
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = (turning / (core::f64::consts::TAU / 64.0)).ceil().max(8.0) as usize;
        for i in 0..=count {
            #[allow(clippy::cast_precision_loss)]
            let f = (i as f64) / (count as f64);
            let t = if reversed {
                range.1 - (range.1 - range.0) * f
            } else {
                range.0 + (range.1 - range.0) * f
            };
            let p = curve.point_at(t, tol)?;
            let d = curve.d1_at(t, tol)?;
            let m = d.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "the path is degenerate at {t}");
            }
            let tangent = if reversed { -(d / m) } else { d / m };
            if let Some(prev) = stations.last()
                && prev.edge == ei
                && prev.at.distance(p) <= tol.confusion()
                && prev.tangent.cross(tangent).magnitude() <= tol.angular()
                && prev.tangent.dot(tangent) > 0.0
            {
                continue;
            }
            stations.push(PathStation {
                at: p,
                tangent,
                edge: ei,
                t,
            });
        }
    }
    if stations.len() < 2 {
        ogeom_bail!(Construction, "the path collapses to a point");
    }
    Ok(stations)
}

/// The first normal taken when none is given: the world `z` made square to
/// `tangent`, or the world `x` where `tangent` runs within thirty degrees
/// of `z`. `tangent` is a unit vector.
#[must_use]
pub fn default_start_normal(tangent: Vector) -> Vector {
    let seed = if tangent.cross(Vector::Z).magnitude() > 0.5 {
        Vector::Z
    } else {
        Vector::X
    };
    let v = seed - tangent * seed.dot(tangent);
    v / v.magnitude()
}

/// One rotation-minimizing step: the normal `n0` at `(p0, t0)` carried to
/// `(p1, t1)` by double reflection, reflected in the plane bisecting the
/// chord and then in the plane bisecting the reflected tangent and `t1`.
/// No travel means no change.
#[must_use]
pub fn carry_normal(p0: Point, t0: Vector, n0: Vector, p1: Point, t1: Vector) -> Vector {
    let v1 = p1 - p0;
    let c1 = v1.dot(v1);
    if c1 <= 1e-20 {
        return n0;
    }
    let nl = n0 - v1 * (2.0 / c1 * v1.dot(n0));
    let tl = t0 - v1 * (2.0 / c1 * v1.dot(t0));
    let v2 = t1 - tl;
    let c2 = v2.dot(v2);
    let next = if c2 > 1e-20 {
        nl - v2 * (2.0 / c2 * v2.dot(nl))
    } else {
        nl
    };
    next / next.magnitude()
}

/// The normal `n` at `from` carried on to `to`: one double-reflection step,
/// or across a corner (two stations on one point) the reflection in the
/// plane halving it.
fn carry_between(from: &PathStation, to: &PathStation, n: Vector) -> Vector {
    let v1 = to.at - from.at;
    if v1.dot(v1) > 1e-20 {
        return carry_normal(from.at, from.tangent, n, to.at, to.tangent);
    }
    // A corner's twin station: no travel to reflect through. For a vector
    // square to the incoming tangent, the reflection in the corner's mitre
    // plane is the turn about the corner's own axis, and its mirror
    // symmetry is what lands both legs' sheared sections on one ring.
    let bisector = from.tangent + to.tangent;
    let m = bisector.magnitude();
    if m <= 1e-12 {
        return n;
    }
    let b = bisector / m;
    n - b * (2.0 * n.dot(b))
}

/// Rotation-minimizing normals at each station, `start` at the first (a
/// unit vector square to its tangent), carried by [`carry_normal`] and
/// across each corner by the turn about the corner's axis.
#[must_use]
pub fn carry_normals(stations: &[PathStation], start: Vector) -> Vec<Vector> {
    let mut normals: Vec<Vector> = Vec::with_capacity(stations.len());
    if stations.is_empty() {
        return normals;
    }
    normals.push(start);
    for pair in stations.windows(2) {
        let n = normals[normals.len() - 1];
        normals.push(carry_between(&pair[0], &pair[1], n));
    }
    normals
}

/// A frame on a path: the point, the unit tangent in the direction of
/// travel and a unit normal square to it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveFrame {
    /// The point on the path.
    pub at: Point,
    /// The unit tangent, in the direction the path is walked.
    pub tangent: Direction,
    /// The unit normal, square to the tangent.
    pub normal: Direction,
    /// The share of the path's length the frame stands at, in `[0, 1]`.
    pub share: f64,
    /// The edge it stands on, by position along the path.
    pub edge: usize,
    /// Its parameter on that edge's curve.
    pub parameter: f64,
}

impl CurveFrame {
    /// The third axis, `tangent × normal`: the frame is right-handed.
    #[must_use]
    pub fn binormal(&self) -> Vector {
        self.tangent.vector().cross(self.normal.vector())
    }
}

/// Rotation-minimizing frames along a path, as
/// [`rotation_minimizing_frames`] returns them.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveFrames {
    /// One frame per station asked for, in the order asked.
    pub frames: Vec<CurveFrame>,
    /// On a closed path, the angle in radians, right-handed about the
    /// start tangent and in `(-π, π]`, from the start normal to the normal
    /// carried round the whole path and back across its closing join (the
    /// corner turn there included, if the path closes on a corner): the
    /// path's holonomy. `None` on an open path.
    pub closing_twist: Option<f64>,
}

impl CurveFrames {
    /// The frames with the closing twist paid off: each normal turned
    /// about its tangent by `-closing_twist` times its share of the path's
    /// length, so the frame at the end of a closed path is the frame at its
    /// start. The frames unchanged on an open path.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// only if a frame holds a non-finite normal.
    pub fn closed_up(&self, tol: Tolerances) -> OgeomResult<Vec<CurveFrame>> {
        let Some(twist) = self.closing_twist else {
            return Ok(self.frames.clone());
        };
        self.frames
            .iter()
            .map(|f| {
                let phi = -twist * f.share;
                let (t, n) = (f.tangent.vector(), f.normal.vector());
                let v = n * phi.cos() + t.cross(n) * phi.sin();
                let v = v - t * v.dot(t);
                Ok(CurveFrame {
                    normal: Direction::new(v, tol)?,
                    ..*f
                })
            })
            .collect()
    }
}

/// Rotation-minimizing frames along `path` (an edge or a wire) at each of
/// `stations`, given as shares of the path's length in `[0, 1]`.
///
/// The first normal is `start_normal` made square to the start tangent, or
/// by default the world `z` made square to it (the world `x` where the
/// tangent runs within thirty degrees of `z`). From there it is carried by
/// double reflection through the stations the pipe sweeps carry their
/// frames through (see [`path_stations`]), and to a share between two of
/// them by one more step from the one behind. So:
///
/// - On a straight run the frame is constant, and through an inflection
///   it carries on as anywhere else: it never needs the curvature.
/// - Along a planar path a normal starting square to the plane stays
///   square to it exactly, and one starting in the plane stays in it.
/// - At a sharp corner of a wire the frame turns about the corner's axis
///   by the corner's angle. A station exactly on a corner takes the frame
///   leaving it; the end of the path takes the frame arriving.
/// - A closed path (one whose end meets its start) brings the frame home
///   turned by its holonomy; that turn is
///   [`CurveFrames::closing_twist`] and is not spread into the frames.
///   [`CurveFrames::closed_up`] spreads it.
///
/// The frames agree, at the same stations, with the sections of a pipe
/// swept along the same path under the rotation-minimizing law, up to the
/// one fixed turn its profile's orientation sets; round a closed path the
/// pipe pays off the closing twist as [`CurveFrames::closed_up`] does.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a station is
/// not a finite share in `[0, 1]`;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `start_normal` is not finite or runs along the start tangent, or as
/// [`path_stations`]; as [`parameter_at_length`] if a length cannot be
/// measured.
pub fn rotation_minimizing_frames(
    model: &Model,
    path: &Shape,
    stations: &[f64],
    start_normal: Option<Vector>,
    tol: Tolerances,
) -> OgeomResult<CurveFrames> {
    if let Some(bad) = stations
        .iter()
        .find(|s| !s.is_finite() || **s < 0.0 || **s > 1.0)
    {
        ogeom_bail!(
            Domain,
            "a station is a share of the path's length in [0, 1], not {bad}"
        );
    }
    let walk = path_stations(model, path, tol)?;
    let t0 = walk[0].tangent;
    let n0 = match start_normal {
        None => default_start_normal(t0),
        Some(v) => {
            if !v.is_finite() {
                ogeom_bail!(Construction, "the start normal {v:?} is not finite");
            }
            let square = v - t0 * v.dot(t0);
            if v.magnitude() <= tol.confusion()
                || square.magnitude() <= tol.angular() * v.magnitude()
            {
                ogeom_bail!(
                    Construction,
                    "the start normal {v:?} runs along the path's start tangent; \
                     it sets no frame"
                );
            }
            square / square.magnitude()
        }
    };
    let normals = carry_normals(&walk, n0);
    let last = walk.len() - 1;
    let closing_twist = (walk[0].at.distance(walk[last].at) <= 10.0 * tol.confusion()).then(|| {
        let home = carry_between(&walk[last], &walk[0], normals[last]);
        n0.cross(home).dot(t0).atan2(n0.dot(home))
    });

    let edges = path_edges(model, path)?;
    let mut curves = Vec::with_capacity(edges.len());
    let mut lengths = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range, reversed) = walked_curve(model, edge, tol)?;
        lengths.push(curve_length(&curve, range, tol)?);
        curves.push((curve, range, reversed));
    }
    let total: f64 = lengths.iter().sum();
    if total <= tol.confusion() {
        ogeom_bail!(Construction, "the path collapses to a point");
    }
    let mut frames = Vec::with_capacity(stations.len());
    for &share in stations {
        let target = share * total;
        // The edge the share falls on: the first that runs past it, so a
        // share on a corner falls on the leg leaving it; the end of the
        // path on the last edge with any length.
        let mut start = 0.0;
        let mut found = None;
        for (e, &length) in lengths.iter().enumerate() {
            if length > 0.0 {
                found = Some((e, start));
                if start + length > target {
                    break;
                }
            }
            start += length;
        }
        let Some((e, from)) = found else {
            ogeom_bail!(Construction, "the path collapses to a point");
        };
        let (curve, range, reversed) = &curves[e];
        let local = (target - from).clamp(0.0, lengths[e]);
        let parameter = if *reversed {
            parameter_at_length(curve, *range, lengths[e] - local, tol)?
        } else {
            parameter_at_length(curve, *range, local, tol)?
        };
        let p = curve.point_at(parameter, tol)?;
        let d = curve.d1_at(parameter, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "the path is degenerate at {parameter}");
        }
        let tangent = if *reversed { -(d / m) } else { d / m };
        // The station behind on this edge: the last of its stations the
        // walk reaches at or before the parameter. Its first is the
        // edge's start, so one always is.
        let behind = |s: &PathStation| {
            s.edge == e
                && if *reversed {
                    s.t >= parameter
                } else {
                    s.t <= parameter
                }
        };
        let Some(k) = walk.iter().rposition(behind) else {
            ogeom_bail!(Construction, "the path has no station on edge {e}");
        };
        let n = carry_normal(walk[k].at, walk[k].tangent, normals[k], p, tangent);
        let n = n - tangent * n.dot(tangent);
        frames.push(CurveFrame {
            at: p,
            tangent: Direction::new(tangent, tol)?,
            normal: Direction::new(n, tol)?,
            share,
            edge: e,
            parameter,
        });
    }
    Ok(CurveFrames {
        frames,
        closing_twist,
    })
}
