//! Rounds on seats only the march speaks (ADR 0001, amendment "marched
//! seats"): a crease between faces no closed-form section reads (a spline
//! edge, the seam of two cylinders crossing, a converted solid's every
//! edge). Each edge's strip is the march's band, closed into its wedge by
//! its legs on the two faces and its caps, built on the solid as given
//! with the request's other edges as its mates and handed to the engine
//! unapplied; the engine applies the wedges in one boolean each, as it
//! applies every tool of the one-boolean fallback.

use core::cell::RefCell;

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::fallback::{Sense, Tool};
use crate::support::Wedge;

/// What seats an edge's round, as the engine reads it for the march.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Seat {
    /// A crease no closed-form section reads: the march's.
    Marched,
    /// A circle between a planar cap and a coaxial cylinder or cone.
    Rim,
    /// A straight crease between planes, or between planes or drums
    /// parallel to it.
    Straight,
}

/// How `edge` of `solid` is seated: a straight edge between planes or a
/// ruled pair of hosts, a rim on a planar cap and its coaxial wall, and
/// any other crease the march's; a rim on a cone the march's too where
/// `fitted` asks its band fitted through the ball's arcs.
fn seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    fitted: bool,
    tol: Tolerances,
) -> OgeomResult<Seat> {
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    Ok(match curve {
        Curve::Line(_) if crate::support::hosts_planar(model, solid, edge, tol)? => Seat::Straight,
        Curve::Line(_) => match crate::ruled::ruled_seat(model, solid, edge, tol)? {
            Some(_) => Seat::Straight,
            None => Seat::Marched,
        },
        Curve::Circle(c) => match crate::support::revolved_seat(model, solid, edge, &c, tol) {
            Ok(seat) if fitted && seat.slope != 0.0 => Seat::Marched,
            Ok(_) => Seat::Rim,
            Err(OgeomError::Construction(_)) => Seat::Marched,
            Err(e) => return Err(e),
        },
        _ => Seat::Marched,
    })
}

/// What the march's bands came to.
enum Bands {
    /// The wedges, each set aside on the solid as given, or, for a lone
    /// edge whose leg lies on a sphere charted with its pole turned away,
    /// on the solid rebuilt on that chart: that solid, with its history
    /// from the solid as given.
    Wedges(Vec<Wedge>, Option<Box<Built>>),
    /// One edge's ball or band refused: a one-edge request may take the
    /// strip beside it off whole.
    Refused,
    /// Not the march's to build: the request goes on.
    Declined,
}

/// `once` rounded on its marched seats: every edge's wedge cut where the
/// edge is convex and fused where it is concave, each wedge's faces named
/// for the edge it rounds (or the edge its builder credits them to). A
/// one-edge request whose ball does not fit beside it, whose band does not
/// build or whose band the boolean refuses, takes the strip beside it off
/// whole where one comes off ([`super::strip_taken`]); so does a rim asked
/// alone. `None` where no edge is a marched seat and none is taken off
/// whole, a wedge does not build, a boolean refuses, or the result is not
/// valid, nor where a one-edge request's band is refused and no strip
/// comes off whole: the request goes on then. Where `fitted`, a rim on a
/// cone is a marched seat, its band fitted through the ball's arcs.
///
/// # Errors
///
/// As the march and the booleans, other than their refusals, which
/// decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    fitted: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    building_tools(|| seated(model, solid, once, fitted, tol))
}

/// Run an engine tool builder with the old-constructor back-edge guard.
pub(super) fn building_tools<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            IN_SEATS.with(|flag| flag.set(self.0));
        }
    }
    let _restore = Restore(IN_SEATS.with(|flag| flag.replace(true)));
    f()
}

thread_local! {
    /// Whether the engine's marched seats are at work on this thread.
    static IN_SEATS: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// Where an old marched constructor named `entry` (one that builds and
/// applies its own wedge, or reads a strip off whole for the one-edge
/// entries) is reached: a record for the old-path count, and in a debug
/// build a panic, where the engine's marched seats are at work. Their
/// tools are the engine's own; reaching one of these from there is a
/// regression.
pub(crate) fn old_marched_constructor(entry: &str) {
    if IN_SEATS.with(core::cell::Cell::get) {
        #[cfg(feature = "old-path-log")]
        let _old_path_entry = crate::fillet::to_current_path(entry);
        if cfg!(debug_assertions) {
            panic!("the engine's marched seats reached the old constructor {entry}");
        }
    }
}

/// `f` run as the engine's other routes run, outside the marched seats:
/// a crease restored under a strip taken off is rounded by the whole
/// engine, whose other routes are not the seats'.
pub(super) fn outside_seats<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            IN_SEATS.with(|flag| flag.set(self.0));
        }
    }
    let _restore = Restore(IN_SEATS.with(|flag| flag.replace(false)));
    f()
}

/// The wedges `once`'s edges on marched seats are rounded by, each built
/// on `solid` as given with the request's other edges as its mates and set
/// aside unapplied, for a route to apply with its own tools. `None` where
/// no edge is a marched seat, or the march declines or refuses an edge.
///
/// # Errors
///
/// As the march and the booleans, other than their refusals, which
/// decline.
pub(crate) fn wedges(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Wedge>>> {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            IN_SEATS.with(|flag| flag.set(self.0));
        }
    }
    let _restore = Restore(IN_SEATS.with(|flag| flag.replace(true)));
    let mut seats = Vec::with_capacity(once.len());
    for (edge, _) in once {
        seats.push(seat(model, solid, edge, false, tol)?);
    }
    if !seats.contains(&Seat::Marched) {
        return Ok(None);
    }
    Ok(match bands(model, solid, once, &seats, tol)? {
        // Wedges set aside against a rebuilt solid are not the solid as
        // given's to apply.
        Bands::Wedges(wedges, None) => Some(wedges),
        Bands::Wedges(_, Some(_)) | Bands::Refused | Bands::Declined => None,
    })
}

/// [`rounds`], its seats read.
fn seated(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    fitted: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut seats = Vec::with_capacity(once.len());
    for (edge, _) in once {
        seats.push(seat(model, solid, edge, fitted, tol)?);
    }
    let alone = match once {
        [(edge, radius)] => Some((edge, *radius)),
        _ => None,
    };
    if !seats.contains(&Seat::Marched) {
        // A rim asked alone whose own round refused: a ring beside it
        // taken off whole.
        return match (alone, seats.as_slice()) {
            (Some((edge, radius)), [Seat::Rim]) => {
                super::strip_taken::rounds(model, solid, edge, radius, tol)
            }
            _ => Ok(None),
        };
    }
    let mut scratch = model.clone();
    match bands(&mut scratch, solid, once, &seats, tol)? {
        Bands::Wedges(wedges, turned) => {
            let on = turned.as_ref().map_or(solid, |t| &t.shape);
            if let Some(built) = applied(&mut scratch, on, once, &wedges, tol)? {
                *model = scratch;
                return Ok(Some(match turned {
                    Some(turned) => Built::new(built.shape, turned.history.then(&built.history)),
                    None => built,
                }));
            }
            let Some((edge, radius)) = alone else {
                return Ok(None);
            };
            super::strip_taken::rounds(model, solid, edge, radius, tol)
        }
        Bands::Refused => {
            let Some((edge, radius)) = alone else {
                return Ok(None);
            };
            super::strip_taken::rounds(model, solid, edge, radius, tol)
        }
        Bands::Declined => Ok(None),
    }
}

/// The wedges the march rounds `once`'s edges by, each built on `solid` as
/// given with the request's other edges as its mates. Where two of the
/// edges' creases turn at a shared vertex, their bands meet there (see
/// [`super::marched_junctions`]).
fn bands(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    seats: &[Seat],
    tol: Tolerances,
) -> OgeomResult<Bands> {
    let all_marched = seats.iter().all(|s| *s == Seat::Marched);
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let radii: Vec<f64> = once.iter().map(|(_, r)| *r).collect();
    // A ball that does not fit beside its edge, or that rests on a far
    // edge, is not the march's band.
    match super::request_fits(model, solid, &edges, &radii, tol) {
        Ok(rolled) if !rolled.iter().any(|r| *r) => {}
        Err(crate::support::Blamed {
            error: e @ OgeomError::Cancelled,
            ..
        }) => return Err(e),
        Err(crate::support::Blamed {
            error: OgeomError::Construction(_),
            ..
        }) if all_marched && once.len() == 1 => return Ok(Bands::Refused),
        _ => return Ok(Bands::Declined),
    }
    let mut ends: Vec<Vec<Shape>> = Vec::with_capacity(edges.len());
    for edge in &edges {
        ends.push(
            ogeom_algo::edge_vertices(model, edge)?
                .map(|(a, b)| vec![a, b])
                .unwrap_or_default(),
        );
    }
    let at_vertex = |v: &Shape| ends.iter().flatten().filter(|w| w.is_same(v)).count();
    if !all_marched && ends.iter().flatten().any(|v| at_vertex(v) >= 3) {
        return Ok(Bands::Declined);
    }
    let met = match super::marched_junctions::met(
        model,
        solid,
        &edges,
        &radii,
        &ends,
        all_marched,
        tol,
    ) {
        Ok(Some(met)) => met,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Ok(None) | Err(_) => return Ok(Bands::Declined),
    };
    let mut set_aside: Vec<(usize, Wedge)> = Vec::new();
    let mut turned = None;
    for (index, ((edge, radius), seat)) in once.iter().zip(seats).enumerate() {
        let aside = RefCell::new(Vec::new());
        let built = match seat {
            Seat::Marched => crate::marched::marched_wedge(
                model,
                solid,
                edge,
                *radius,
                Some((index, &met.mates)),
                &aside,
                tol,
            ),
            Seat::Rim | Seat::Straight => super::marched_junctions::other_band(
                model,
                solid,
                edge,
                *radius,
                (index, &met.mates),
                &aside,
                tol,
            ),
        };
        let wedges = aside.into_inner();
        match built {
            Ok(built) if built.shape.is_same(solid) && !wedges.is_empty() => {}
            // A lone edge's leg on a sphere charted with its pole turned
            // away: the march hands back the solid rebuilt on that chart,
            // its faces on the sphere replaced and every edge kept, for
            // the wedge to melt against. The history it hands back ends in
            // the edge's own deletion, which the wedge's application
            // records again with the faces it generates.
            Ok(built) if once.len() == 1 && !wedges.is_empty() => {
                turned = Some(Box::new(Built::new(
                    built.shape,
                    without(&built.history, edge),
                )));
            }
            Err(e @ OgeomError::Cancelled) => return Err(e),
            Err(OgeomError::Construction(_) | OgeomError::NotDone(_)) if once.len() == 1 => {
                return Ok(Bands::Refused);
            }
            _ => return Ok(Bands::Declined),
        }
        for wedge in wedges {
            set_aside.push((index, wedge));
        }
    }
    match super::marched_junctions::joined(model, solid, once, &met, set_aside, tol) {
        Ok(joined) => Ok(Bands::Wedges(joined, turned)),
        Err(e @ OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(Bands::Declined),
    }
}

/// `history` with nothing recorded of `input`.
fn without(history: &ogeom_algo::History, input: &Shape) -> ogeom_algo::History {
    let mut kept = ogeom_algo::History::new();
    for shape in history.inputs() {
        if shape.is_same(input) {
            continue;
        }
        if history.is_deleted(&shape) {
            kept.delete(&shape);
        }
        if let Some(copy) = history.copy_of(&shape) {
            kept.copy(&shape, copy.clone());
        } else {
            for image in history.modified(&shape) {
                kept.modify(&shape, image.clone());
            }
        }
        for made in history.generated(&shape) {
            kept.generate(&shape, made.clone());
        }
    }
    kept
}

/// `wedges` applied to `solid` as the one-boolean fallback's tools, each
/// wedge's faces named for its edge or the edge it credits them to;
/// `None` where a boolean refuses, an edge of `once` is left nothing of,
/// or the result is not valid.
pub(super) fn applied(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    wedges: &[Wedge],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut tools = Vec::with_capacity(wedges.len());
    for wedge in wedges {
        let Some(own) = wedge.edges.first() else {
            return Ok(None);
        };
        let faces = explore_unique(model, &wedge.solid, ShapeType::Face)?
            .into_iter()
            .map(|face| {
                let from = wedge
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or_else(|| own.clone(), |(_, e)| e.clone());
                (face, from)
            })
            .collect();
        tools.push(Tool {
            solid: wedge.solid.clone(),
            sense: if wedge.additive {
                Sense::Fuse
            } else {
                Sense::Cut
            },
            faces,
        });
    }
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(model, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    Ok(Some(built))
}
