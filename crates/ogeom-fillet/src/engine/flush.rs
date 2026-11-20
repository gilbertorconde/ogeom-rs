//! A round's free end at an inside corner, flush on the face it meets
//! there (ADR 0001, amendment "M5e", a free end; SMLib's mixed 3x1 corner:
//! "the end face will be extended out to the fillet surface"). The end
//! face's surface, a plane or a cylinder, cone or torus (an earlier
//! round's band among them), is carried on past the face, and the round's
//! tool loses what lies on its far side from the crease near that end.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// The surface a round's band ends flush on, and which side of it is the
/// crease's.
#[derive(Clone)]
pub(crate) enum EndSurface {
    /// A plane through `at`, its unit normal `keep` towards the crease.
    Plane { at: Point, keep: Vector },
    /// A cylinder, the crease's side outside it where `outside`.
    Cylinder {
        cylinder: ogeom_math::Cylinder,
        outside: bool,
    },
    /// A cone, the crease's side outside it where `outside`.
    Cone {
        cone: ogeom_math::Cone,
        outside: bool,
    },
    /// A torus, the crease's side outside its tube where `outside`, and
    /// then round the ring from the axis where `ringward` (farther from the
    /// axis than the tube's centre circle). Past the surface from outside
    /// lies the tube and, within the tube's height, the ring's hole on the
    /// far side of its centre circle: a way carried on through the tube
    /// does not come out again on the other side.
    Torus {
        torus: ogeom_math::Torus,
        outside: bool,
        ringward: bool,
    },
}

impl EndSurface {
    /// The end face's surface `surface` through `at`, the crease's side the
    /// side `crease` (a point just back along the crease from its end)
    /// stands on, `normal` the face's outward normal at `at`. `None` for a
    /// surface of another kind.
    pub(crate) fn of(
        surface: &SurfaceGeometry,
        at: Point,
        normal: Vector,
        crease: Point,
    ) -> Option<Self> {
        Some(match surface {
            SurfaceGeometry::Plane(_) => {
                let keep = if (crease - at).dot(normal) >= 0.0 {
                    normal
                } else {
                    -normal
                };
                Self::Plane { at, keep }
            }
            SurfaceGeometry::Cylinder(c) => {
                let cylinder = c.cylinder();
                let outside = off_axis(cylinder.axis(), crease) >= cylinder.radius();
                Self::Cylinder { cylinder, outside }
            }
            SurfaceGeometry::Cone(c) => {
                let cone = c.cone();
                let axis = cone.axis();
                let z = axis.direction.vector().dot(crease - axis.location);
                let outside = off_axis(axis, crease) >= cone.radius_at(z);
                Self::Cone { cone, outside }
            }
            SurfaceGeometry::Torus(t) => {
                let torus = t.torus();
                let outside = torus.signed_distance_to(crease) >= 0.0;
                let ringward = off_axis(torus.axis(), crease) >= torus.major_radius();
                Self::Torus {
                    torus,
                    outside,
                    ringward,
                }
            }
            _ => return None,
        })
    }
}

impl EndSurface {
    /// Whether `p` stands on the crease's side of the surface.
    pub(crate) fn keeps(&self, p: Point) -> bool {
        match self {
            Self::Plane { at, keep } => (p - *at).dot(*keep) >= 0.0,
            Self::Cylinder { cylinder, outside } => {
                (off_axis(cylinder.axis(), p) >= cylinder.radius()) == *outside
            }
            Self::Cone { cone, outside } => {
                let axis = cone.axis();
                let z = axis.direction.vector().dot(p - axis.location);
                (off_axis(axis, p) >= cone.radius_at(z)) == *outside
            }
            Self::Torus {
                torus,
                outside,
                ringward,
            } => {
                let in_tube = torus.signed_distance_to(p) < 0.0;
                if !*outside {
                    return in_tube;
                }
                let axis = torus.axis();
                let height = axis.direction.vector().dot(p - axis.location);
                let beside = height.abs() <= torus.minor_radius();
                let hole_side = (off_axis(axis, p) < torus.major_radius()) == *ringward;
                !(in_tube || (beside && hole_side))
            }
        }
    }
}

/// How far `p` stands from `axis`.
fn off_axis(axis: ogeom_math::Axis, p: Point) -> f64 {
    let d = axis.direction.vector();
    let off = p - axis.location;
    (off - d * d.dot(off)).magnitude()
}

/// `tool` with what lies on the far side of `end` from the crease, within
/// `span` of the vertex `at`, cut away, and `faces` of it traced through
/// the cut. A curved end face can leave the tool in pieces, its way
/// carried on through the surface and out again past it (an earlier
/// round's band, a tube): the piece holding `crease`, a point of the tool
/// on the crease's own stretch, is the band, and what lies past the
/// surface is past its end. `None` where a boolean refuses or no one piece
/// holds `crease`.
///
/// # Errors
///
/// As building the cutting solids and classifying `crease`.
#[allow(clippy::too_many_arguments, reason = "one cut, all its data")]
pub(crate) fn cut_far(
    model: &mut Model,
    tool: &Shape,
    faces: &[Shape],
    end: &EndSurface,
    at: Point,
    span: f64,
    crease: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let Some(far) = far_side(model, end, at, span, tol)? else {
        return Ok(None);
    };
    let Ok(cut) = ogeom_bool::cut(model, tool, &far, tol) else {
        return Ok(None);
    };
    let solids = explore_unique(model, &cut.shape, ShapeType::Solid)?;
    let kept = match (solids.as_slice(), crease) {
        ([one], _) => one,
        ([], _) | (_, None) => return Ok(None),
        (many, Some(p)) => {
            let mut holding = Vec::new();
            for s in many {
                if ogeom_algo::classify_in_solid_exact(model, s, p, tol)?
                    == ogeom_algo::Containment::In
                {
                    holding.push(s);
                }
            }
            let [one] = holding[..] else {
                return Ok(None);
            };
            one
        }
    };
    let own = explore_unique(model, kept, ShapeType::Face)?;
    let traced = faces
        .iter()
        .flat_map(|f| cut.history.trace(f).to_vec())
        .filter(|f| own.iter().any(|g| g.is_same(f)))
        .collect();
    Ok(Some((kept.clone(), traced)))
}

/// `tool` with the material of `solid` on the far side of `end`, within
/// `span` of `at`, cut away, and `faces` of it traced through the cut: a
/// tool standing off a wall past its surface takes nothing the solid holds
/// beyond that wall (a block fused against a drum, at an inside corner).
/// `None` where a boolean refuses or leaves other than one solid.
///
/// # Errors
///
/// As building the cutting solids.
#[allow(clippy::too_many_arguments, reason = "one cut, all its data")]
pub(crate) fn cut_material_far(
    model: &mut Model,
    tool: &Shape,
    faces: &[Shape],
    end: &EndSurface,
    at: Point,
    span: f64,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let Some(far) = far_side(model, end, at, span, tol)? else {
        return Ok(None);
    };
    let Ok(held) = ogeom_bool::common(model, &far, solid, tol) else {
        return Ok(None);
    };
    if explore_unique(model, &held.shape, ShapeType::Solid)?.is_empty() {
        return Ok(Some((tool.clone(), faces.to_vec())));
    }
    let Ok(cut) = ogeom_bool::cut(model, tool, &held.shape, tol) else {
        return Ok(None);
    };
    let solids = explore_unique(model, &cut.shape, ShapeType::Solid)?;
    let [kept] = solids.as_slice() else {
        return Ok(None);
    };
    let traced = faces
        .iter()
        .flat_map(|f| cut.history.trace(f).to_vec())
        .collect();
    Ok(Some((kept.clone(), traced)))
}

/// What lies on the far side of `end` from the crease within the cube of
/// side `2 span` about `at`.
fn far_side(
    model: &mut Model,
    end: &EndSurface,
    at: Point,
    span: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    if let EndSurface::Plane { at, keep } = end {
        return Ok(Some(behind(model, *at, *keep, span, tol)?));
    }
    let cube = {
        let corner = at - Vector::new(span, span, span);
        ogeom_algo::make_box(
            model,
            Frame::new(corner, Direction::Z, Direction::X, tol)?,
            (span * 2.0, span * 2.0, span * 2.0),
            tol,
        )?
        .shape
    };
    let reach = 2.0 * span;
    // The solid the surface bounds, through the cube.
    let (inside, outside) = match end {
        EndSurface::Plane { .. } => unreachable!("planes return above"),
        EndSurface::Cylinder { cylinder, outside } => {
            let axis = cylinder.axis();
            let d = axis.direction.vector();
            let foot = axis.location + d * d.dot(at - axis.location);
            let frame = Frame::new(foot - d * reach, axis.direction, across(d, tol)?, tol)?;
            let drum =
                ogeom_algo::make_cylinder(model, frame, cylinder.radius(), 2.0 * reach, tol)?.shape;
            (drum, *outside)
        }
        EndSurface::Cone { cone, outside } => {
            let axis = cone.axis();
            let d = axis.direction.vector();
            let z = d.dot(at - axis.location);
            let (lo, hi) = (z - reach, z + reach);
            let (r_lo, r_hi) = (cone.radius_at(lo), cone.radius_at(hi));
            // The cube must not reach the apex: a frustum is no solid there.
            if r_lo <= tol.confusion() || r_hi <= tol.confusion() {
                return Ok(None);
            }
            let frame = Frame::new(axis.location + d * lo, axis.direction, across(d, tol)?, tol)?;
            let frustum = ogeom_algo::make_cone(model, frame, r_lo, r_hi, hi - lo, tol)?.shape;
            (frustum, *outside)
        }
        EndSurface::Torus {
            torus,
            outside,
            ringward,
        } => {
            let frame = torus.frame();
            let (big, small) = (torus.major_radius(), torus.minor_radius());
            let ring = ogeom_algo::make_torus(model, frame, big, small, tol)?.shape;
            if !*outside {
                (ring, false)
            } else {
                // The tube with the slab of the tube's height on the far
                // side of its centre circle from the crease.
                let z = frame.z().vector();
                let foot = Frame::new(frame.origin() - z * small, frame.z(), frame.x(), tol)?;
                let puck = ogeom_algo::make_cylinder(model, foot, big, 2.0 * small, tol)?.shape;
                let slab = if *ringward {
                    puck
                } else {
                    let wide = off_axis(torus.axis(), at) + small + 2.0 * reach;
                    let disc =
                        ogeom_algo::make_cylinder(model, foot, wide, 2.0 * small, tol)?.shape;
                    let Ok(annulus) = ogeom_bool::cut(model, &disc, &puck, tol) else {
                        return Ok(None);
                    };
                    annulus.shape
                };
                let Ok(past) = ogeom_bool::fuse(model, &ring, &slab, tol) else {
                    return Ok(None);
                };
                (past.shape, true)
            }
        }
    };
    let far = if outside {
        ogeom_bool::common(model, &cube, &inside, tol)
    } else {
        ogeom_bool::cut(model, &cube, &inside, tol)
    };
    Ok(far.ok().map(|b| b.shape))
}

/// A unit direction square to the unit `d`.
fn across(d: Vector, tol: Tolerances) -> OgeomResult<Direction> {
    let x = if d.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    Direction::new(x - d * x.dot(d), tol)
}

/// The cube of side `2 span` on the far side of the plane through `at`
/// from `keep`, one face on that plane, centred on `at` across it.
fn behind(
    model: &mut Model,
    at: Point,
    keep: Vector,
    span: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let inward = Direction::new(-keep, tol)?;
    let across = across(keep, tol)?;
    let frame = Frame::new(at, inward, across, tol)?;
    let corner = at - (frame.x().vector() + frame.y().vector()) * span;
    Ok(ogeom_algo::make_box(
        model,
        Frame::new(corner, inward, across, tol)?,
        (span * 2.0, span * 2.0, span * 2.0),
        tol,
    )?
    .shape)
}

/// `built` with each of the faces `around` (new faces of it) unified with
/// the faces on its carrier it reaches, and the edge chains that leaves on
/// one curve joined: a band cut flush on a face's surface leaves a piece of
/// that surface beside the face, and a crease running on along both stays
/// one edge to its real end. Where merging leaves the solid unfit, `built`
/// as it is.
///
/// # Errors
///
/// As unifying and joining, and a cancelled guard.
pub(crate) fn rejoined(
    model: &mut Model,
    built: ogeom_algo::Built,
    around: &[Shape],
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    if around.is_empty() {
        return Ok(built);
    }
    let mut scratch = model.clone();
    let (unified, merged) =
        ogeom_heal::unify_same_domain_around(&mut scratch, &built.shape, around, tol)?;
    if merged == 0 {
        return Ok(built);
    }
    let (joined, _) = ogeom_heal::merge_edges(&mut scratch, &unified.shape, tol)?;
    let history = rejoined_edges(&scratch, &unified.shape, &joined)?;
    let shape = joined.shape;
    match super::guards::solid_fits(&scratch, &shape, tol) {
        Ok(()) => {
            *model = scratch;
            let history = built.history.then(&unified.history).then(&history);
            Ok(ogeom_algo::Built::new(shape, history))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(built),
    }
}

/// `joined`'s history with each edge of `before` it merged away traced to
/// the edge it joined: the neighbour across a vertex it shares, which
/// became a new edge reaching on over both.
fn rejoined_edges(
    model: &Model,
    before: &Shape,
    joined: &ogeom_algo::Built,
) -> OgeomResult<ogeom_algo::History> {
    let mut history = joined.history.clone();
    let edges = explore_unique(model, before, ShapeType::Edge)?;
    let ends = |e: &Shape| -> OgeomResult<Vec<Shape>> {
        Ok(ogeom_algo::edge_vertices(model, e)?
            .map(|(a, b)| vec![a, b])
            .unwrap_or_default())
    };
    let meets = |a: &[Shape], b: &[Shape]| a.iter().any(|v| b.iter().any(|w| w.is_same(v)));
    for gone in edges.iter().filter(|e| joined.history.is_deleted(e)) {
        let own = ends(gone)?;
        for next in &edges {
            if next.is_same(gone) || !meets(&ends(next)?, &own) {
                continue;
            }
            if let [image] = joined.history.trace(next)
                && !image.is_same(next)
                && meets(&ends(image)?, &own)
            {
                history.modify(gone, image.clone());
                break;
            }
        }
    }
    Ok(history)
}
