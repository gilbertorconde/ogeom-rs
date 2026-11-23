//! The largest size a blend builds at, read from its seat and confirmed
//! against the build itself (ADR 0001, amendment "largest-size readings
//! confirmed by the build").
//!
//! The seat reading (where a contact first leaves its face) is fast and
//! closed-form, but a request can be refused below it for a reason the
//! seat does not see (a corner that does not close, a bevel reaching a
//! concave floor, a result failing its check), and can build above it
//! (a bevel rolling over its face's far side as a notch). So the seat
//! reading only starts the search: the size is then confirmed by asking
//! the same call the request makes, on a copy of the model, at sizes
//! bracketing the boundary between the last size that built and the first
//! that was refused, bisected to [`LIMIT_TOLERANCE`].

use crate::law::Profile;
use crate::{Chamfer, ChamferCorner, FilletCorner, Overflow, RadiusLaw, RoundProfile};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape};

/// The relative tolerance a confirmed limit is bisected to: the smallest
/// size found refused is at most this share above the size returned.
pub const LIMIT_TOLERANCE: f64 = 1e-4;

/// The most builds one reading probes. The bracket and the bisection to
/// [`LIMIT_TOLERANCE`] take at most about 45 from any start in a range of
/// a billion; the seat reading usually brackets the limit in two, and
/// [`LIMIT_SAMPLES`] each side are spent looking for a size that builds
/// again.
pub const LIMIT_PROBES: u32 = 64;

/// How many sizes below the limit, and how many above its first refusal,
/// a reading probes for a range that is not one interval: a refusal below
/// a size that builds, or a size that builds above a refusal.
pub const LIMIT_SAMPLES: u32 = 3;

/// What a largest-size reading found: the size and how it was confirmed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlendLimit {
    /// The largest size confirmed to build, at most the bound asked: the
    /// size of the first boundary from below between sizes that build and
    /// sizes that are refused. None where no size probed above the
    /// confusion distance builds.
    pub size: Option<f64>,
    /// The smallest size probed above `size` that was refused, within
    /// [`LIMIT_TOLERANCE`] of it unless [`LIMIT_PROBES`] ran out first. None
    /// where the bound asked builds.
    pub refused: Option<f64>,
    /// A size above `refused` found to build again: the sizes that build
    /// are not one interval, and `size` is where they are first broken.
    /// None where none of the sizes probed for one builds, which does not
    /// prove there is none.
    pub builds_again: Option<f64>,
    /// The seat reading the search started from, where the seat could be
    /// read: the size where a contact first leaves its face.
    pub seat: Option<f64>,
    /// How many builds were probed.
    pub probes: u32,
}

/// A request whose largest size [`blend_limit`] reads, with the options
/// its blend call takes. The size scales what the request names.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum BlendSizing<'a> {
    /// [`fillet_edges_corner`](crate::fillet_edges_corner) on `edges`
    /// closed as `corner` says; the size is the radius. Every behaviour of
    /// [`Overflow`] is on; [`BlendSizing::RoundsCornerOverflow`] chooses
    /// among them.
    Rounds {
        /// The edges rounded together.
        edges: &'a [Shape],
        /// How a vertex where rounds meet is closed.
        corner: FilletCorner,
    },
    /// [`fillet_edge_distances`](crate::fillet_edge_distances) on `edge`,
    /// meeting `face` `size · on_face` from the edge and the other face
    /// `size · on_other` from it, its section `profile`; the size is the
    /// common scale. Every behaviour of [`Overflow`] is on;
    /// [`BlendSizing::RoundDistancesOverflow`] chooses among them.
    RoundDistances {
        /// The edge rounded.
        edge: &'a Shape,
        /// The face the first distance runs along.
        face: &'a Shape,
        /// The distance along `face` at a scale of one.
        on_face: f64,
        /// The distance along the other face at a scale of one.
        on_other: f64,
        /// The section between the contacts.
        profile: RoundProfile,
    },
    /// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered) on `edge`
    /// alone, closed as `corner` says; the size is the measure
    /// [`largest_chamfer_distance`](crate::largest_chamfer_distance) reads
    /// for `spec`'s spelling, the measure `spec` itself names not read.
    /// One edge has no vertex where three of its bevels meet, so either
    /// corner builds alike there; the corner is the call's own option.
    /// Every behaviour of [`Overflow`] is on; [`BlendSizing::Chamfers`]
    /// chooses among them.
    Chamfer {
        /// The edge bevelled.
        edge: &'a Shape,
        /// The bevel's spelling.
        spec: &'a Chamfer,
        /// How a vertex where three or more bevels meet is closed.
        corner: ChamferCorner,
    },
    /// [`fillet_edges_overflow`](crate::fillet_edges_overflow) on `edges`,
    /// a round reaching geometry other than its own two faces handled as
    /// `overflow` says; the size is the radius. Each vertex where rounds
    /// meet is closed by setbacks, as that call closes it. With
    /// [`Overflow::default`] this reads as [`BlendSizing::Rounds`] with
    /// [`FilletCorner::Setback`].
    RoundsOverflow {
        /// The edges rounded together.
        edges: &'a [Shape],
        /// What a round reaching other geometry does.
        overflow: Overflow,
    },
    /// [`chamfer_edges_overflow`](crate::chamfer_edges_overflow) on every
    /// bevel of `specs` at once, closed as `corner` says, a bevel reaching
    /// geometry other than its own two faces handled as `overflow` says;
    /// the size is the one measure every spec is scaled to, as
    /// [`BlendSizing::Chamfer`] scales one, the measures `specs` name not
    /// read. With [`Overflow::default`] the call is
    /// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered), and one
    /// spec reads as [`BlendSizing::Chamfer`]. The seat is
    /// [`largest_chamfer_distance_all`](crate::largest_chamfer_distance_all)'s
    /// seat reading: each bevel seated on its faces, and no bevel's contact
    /// within another's leg of its edge across a face they share.
    Chamfers {
        /// The bevels, each scaled by the size.
        specs: &'a [(Shape, Chamfer)],
        /// How a vertex where three or more bevels meet is closed.
        corner: ChamferCorner,
        /// What a bevel reaching other geometry does.
        overflow: Overflow,
    },
    /// [`fillet_edges_variable`](crate::fillet_edges_variable) on `chain`
    /// under `law`; the size is the law's value `size` frees, the rest of
    /// the law held as given. Every behaviour of [`Overflow`] is on;
    /// [`BlendSizing::RoundLawOverflow`] chooses among them.
    RoundLaw {
        /// The chain rounded, as the call takes it.
        chain: &'a [Shape],
        /// The law, its freed value as given ignored.
        law: &'a RadiusLaw,
        /// Which value of the law the size is.
        size: LawSize,
    },
    /// [`fillet_edge_variable`](crate::fillet_edge_variable) on `edge`
    /// from `start` to `end`, a linear law whose points are its start (0)
    /// and end (1); the size is the value `size` frees.
    RoundLinear {
        /// The edge rounded.
        edge: &'a Shape,
        /// The radius at the edge's start.
        start: f64,
        /// The radius at its end.
        end: f64,
        /// Which value of the law the size is.
        size: LawSize,
    },
    /// [`fillet_edge_distances_variable`](crate::fillet_edge_distances_variable)
    /// on `edge`, meeting `face` at the distances `on_face` gives and the
    /// other face at those `on_other` gives, its section `profile`; the size
    /// is the value `size` frees. That call takes no [`Overflow`], so this
    /// reading has no form under options.
    DistanceLaws {
        /// The edge rounded.
        edge: &'a Shape,
        /// The face the first law's distance runs along.
        face: &'a Shape,
        /// The distances along `face`.
        on_face: &'a RadiusLaw,
        /// The distances along the other face.
        on_other: &'a RadiusLaw,
        /// The section between the contacts.
        profile: RoundProfile,
        /// Which value of the laws the size is.
        size: DistanceLawSize,
    },
    /// [`fillet_edges_corner_overflow`](crate::fillet_edges_corner_overflow)
    /// on `edges` closed as `corner` says, a round reaching geometry other
    /// than its own two faces handled as `overflow` says; the size is the
    /// radius. The seat is [`BlendSizing::Rounds`]'. With
    /// [`Overflow::default`] this reads as [`BlendSizing::Rounds`], and with
    /// [`FilletCorner::Setback`] as [`BlendSizing::RoundsOverflow`].
    RoundsCornerOverflow {
        /// The edges rounded together.
        edges: &'a [Shape],
        /// How a vertex where rounds meet is closed.
        corner: FilletCorner,
        /// What a round reaching other geometry does.
        overflow: Overflow,
    },
    /// [`fillet_edge_distances_overflow`](crate::fillet_edge_distances_overflow)
    /// on `edge`, as [`BlendSizing::RoundDistances`] scales its distances,
    /// a round reaching geometry other than its own two faces handled as
    /// `overflow` says. The seat is [`BlendSizing::RoundDistances`]'. With
    /// [`Overflow::default`] this reads as [`BlendSizing::RoundDistances`].
    RoundDistancesOverflow {
        /// The edge rounded.
        edge: &'a Shape,
        /// The face the first distance runs along.
        face: &'a Shape,
        /// The distance along `face` at a scale of one.
        on_face: f64,
        /// The distance along the other face at a scale of one.
        on_other: f64,
        /// The section between the contacts.
        profile: RoundProfile,
        /// What a round reaching other geometry does.
        overflow: Overflow,
    },
    /// [`fillet_edges_variable_overflow`](crate::fillet_edges_variable_overflow)
    /// on `chain` under `law`, as [`BlendSizing::RoundLaw`] frees its value,
    /// a round reaching geometry other than its own two faces handled as
    /// `overflow` says. The seat is [`BlendSizing::RoundLaw`]'. With
    /// [`Overflow::default`] this reads as [`BlendSizing::RoundLaw`].
    RoundLawOverflow {
        /// The chain rounded, as the call takes it.
        chain: &'a [Shape],
        /// The law, its freed value as given ignored.
        law: &'a RadiusLaw,
        /// Which value of the law the size is.
        size: LawSize,
        /// What a round reaching other geometry does.
        overflow: Overflow,
    },
}

/// One value of a law: its first point, its last, or the point at an
/// index. A linear law's points are its start (0) and end (1). On a closed
/// chain a law with a point at the chain's start and one at its end gives
/// one point twice, and freeing either frees both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LawStation {
    /// The law's first point.
    Start,
    /// The law's last point.
    End,
    /// The law's point at this index, from 0.
    Point(usize),
}

/// Which size of a radius law a reading frees.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LawSize {
    /// One value, the rest held as given.
    At(LawStation),
    /// The scale of every radius of the law, as given at a scale of one.
    Scale,
}

/// Which size of a pair of distance laws a reading frees.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceLawSize {
    /// One value of the law along the named face, the rest held as given.
    OnFace(LawStation),
    /// One value of the law along the other face, the rest held as given.
    OnOther(LawStation),
    /// The scale of every distance of both laws, as given at a scale of one.
    Scale,
}

/// The largest size up to `up_to` at which `sizing`'s blend call builds on
/// `solid`, with how it was confirmed: the reading behind
/// [`largest_fillet_radius_all`](crate::largest_fillet_radius_all),
/// [`largest_fillet_distances`](crate::largest_fillet_distances),
/// [`largest_chamfer_distance`](crate::largest_chamfer_distance) and
/// [`largest_chamfer_distance_all`](crate::largest_chamfer_distance_all),
/// which return its [`BlendLimit::size`], with the options their calls
/// take.
///
/// The seat is read first, without building, as [`blend_seat`] reads it
/// alone; then the call is asked on a copy of the model, first at the
/// seat reading (at `up_to` where the seat cannot be read) and a
/// [`LIMIT_TOLERANCE`] past it, then doubling or halving the size until
/// one builds and one is refused, and bisecting between the largest size
/// that built and the smallest refused until they stand within
/// [`LIMIT_TOLERANCE`] of each other. [`LIMIT_SAMPLES`] sizes evenly
/// below the result, and as many spread geometrically from the refusal to
/// `up_to`, are then probed: a refusal below moves the answer down to the
/// boundary under it, a build above is reported as
/// [`BlendLimit::builds_again`]. At most [`LIMIT_PROBES`] builds are
/// probed.
///
/// A reading where the seat brackets the limit costs `2 + 2 ·`
/// [`LIMIT_SAMPLES`] builds of the request; one bisected from a factor of
/// two costs about 14 more.
///
/// Under a law ([`BlendSizing::RoundLaw`],
/// [`BlendSizing::RoundLawOverflow`], [`BlendSizing::RoundLinear`],
/// [`BlendSizing::DistanceLaws`]) the size is one value of the law, or the
/// scale of all of them, and each probe asks the call under the whole law
/// so set, over the whole chain. The seat is read at both ends of the edge
/// and every sixteenth between: the largest size at which the law's
/// measure at each station stays within that station's contact's reach.
/// On a chain of several edges it is the constant round's seat held at
/// every sixteenth of the chain.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to`, an edge list or a distance is not usable; where the seat
/// itself is refused and no size builds, that refusal.
pub fn blend_limit(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<BlendLimit> {
    read_seat(model, solid, sizing, up_to, tol, |seat, floor, builds| {
        settle(seat, up_to, floor, builds)
    })
}

/// The size where a contact of `sizing`'s blend first leaves its face, at
/// most `up_to`, read from the seats without building: the reading
/// [`blend_limit`] starts from, its [`BlendLimit::seat`]. None where no
/// size above the confusion distance stays seated.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to`, an edge list or a distance is not usable, as [`blend_limit`]
/// refuses them; otherwise the seat reading's own refusal, such as
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where the seat
/// is not read from the seat (a marched crease).
pub fn blend_seat(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    read_seat(model, solid, sizing, up_to, tol, |seat, _, _| seat)
}

/// The probe [`read_seat`] hands on: whether the request builds at a size.
type Builds<'b> = &'b mut dyn FnMut(f64) -> bool;

/// `sizing`'s seat reading, the smallest size worth probing and the probe,
/// handed to `then`: [`blend_limit`] confirms the seat with the probe,
/// [`blend_seat`] returns it. Nothing is built before `then` asks.
fn read_seat<R>(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    up_to: f64,
    tol: Tolerances,
    then: impl FnOnce(OgeomResult<Option<f64>>, f64, Builds<'_>) -> OgeomResult<R>,
) -> OgeomResult<R> {
    let conf = tol.confusion();
    if !up_to.is_finite() || up_to <= conf {
        ogeom_bail!(
            Construction,
            "a largest size up to {up_to} asks for nothing"
        );
    }
    let (seat, floor) = match *sizing {
        BlendSizing::RoundLaw { .. }
        | BlendSizing::RoundLawOverflow { .. }
        | BlendSizing::RoundLinear { .. }
        | BlendSizing::DistanceLaws { .. } => {
            return law_seat_reading(model, solid, sizing, up_to, tol, then);
        }
        BlendSizing::Rounds { edges, .. }
        | BlendSizing::RoundsOverflow { edges, .. }
        | BlendSizing::RoundsCornerOverflow { edges, .. } => {
            if edges.is_empty() {
                ogeom_bail!(Construction, "a request of no edges has no largest radius");
            }
            (
                crate::fits::seat_fillet_radius(model, solid, edges, up_to, tol),
                conf,
            )
        }
        BlendSizing::RoundDistances {
            edge,
            face,
            on_face,
            on_other,
            ..
        }
        | BlendSizing::RoundDistancesOverflow {
            edge,
            face,
            on_face,
            on_other,
            ..
        } => {
            for d in [on_face, on_other] {
                if !d.is_finite() || d <= conf {
                    ogeom_bail!(
                        Construction,
                        "a round set back {d} from its edge rounds nothing"
                    );
                }
            }
            (
                crate::asym::seat_fillet_distances(
                    model, solid, edge, face, on_face, on_other, up_to, tol,
                ),
                conf / on_face.min(on_other),
            )
        }
        BlendSizing::Chamfer { edge, spec, .. } => {
            let floor = chamfer_floor(spec, conf)?;
            (
                crate::chamfer::seat_chamfer_distance(model, solid, edge, spec, up_to, tol),
                floor,
            )
        }
        BlendSizing::Chamfers { specs, .. } => {
            if specs.is_empty() {
                ogeom_bail!(
                    Construction,
                    "a request of no edges has no largest distance"
                );
            }
            let mut floor = conf;
            for (_, spec) in specs {
                floor = floor.max(chamfer_floor(spec, conf)?);
            }
            (
                crate::chamfer::seat_chamfer_distance_all(model, solid, specs, up_to, tol),
                floor,
            )
        }
    };
    then(
        seat,
        floor,
        &mut counted(|size| probe(model, solid, sizing, size, tol)),
    )
}

/// The smallest measure worth asking of `spec`: the one at which its
/// shorter distance is a confusion distance.
fn chamfer_floor(spec: &Chamfer, conf: f64) -> OgeomResult<f64> {
    match spec {
        Chamfer::Distances {
            on_face, on_other, ..
        } => {
            for d in [*on_face, *on_other] {
                if !d.is_finite() || d <= conf {
                    ogeom_bail!(
                        Construction,
                        "a bevel set back {d} from its edge bevels nothing"
                    );
                }
            }
            Ok(conf / on_face.min(*on_other))
        }
        Chamfer::Symmetric(_) | Chamfer::Angle { .. } => Ok(conf),
    }
}

/// The confirmed reading from `seat`: searched from the bound where the
/// seat was refused, that refusal returned where no size builds.
fn settle(
    seat: OgeomResult<Option<f64>>,
    up_to: f64,
    floor: f64,
    builds: impl FnMut(f64) -> bool,
) -> OgeomResult<BlendLimit> {
    match seat {
        Ok(seat) => Ok(confirm(Some(seat), up_to, floor, builds)),
        Err(refusal) => {
            let found = confirm(None, up_to, floor, builds);
            if found.size.is_none() {
                return Err(refusal);
            }
            Ok(found)
        }
    }
}

/// The largest value up to `up_to` that the point `at` of `law` can take,
/// the rest of the law held, at which
/// [`fillet_edges_variable`](crate::fillet_edges_variable) on `chain` of
/// `solid` builds: [`blend_limit`] with [`BlendSizing::RoundLaw`]. Every
/// value probed below it builds and one [`LIMIT_TOLERANCE`] above it is
/// refused, the law checked over the whole chain by the build itself.
/// None where no value above the confusion distance builds.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `up_to` is not a usable length, `chain` is empty, `at` names no point
/// of the law, or the law at `up_to` states no usable radius.
pub fn largest_law_radius(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    at: LawStation,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let sizing = BlendSizing::RoundLaw {
        chain,
        law,
        size: LawSize::At(at),
    };
    Ok(blend_limit(model, solid, &sizing, up_to, tol)?.size)
}

/// The largest scale `k` up to `up_to` at which
/// [`fillet_edges_variable`](crate::fillet_edges_variable) on `chain` of
/// `solid` under `law` with every radius times `k` builds, as
/// [`largest_law_radius`] reads a value. A law holding one radius `r`
/// answers the constant round's largest radius over `r`.
///
/// # Errors
///
/// As [`largest_law_radius`].
pub fn largest_law_scale(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let sizing = BlendSizing::RoundLaw {
        chain,
        law,
        size: LawSize::Scale,
    };
    Ok(blend_limit(model, solid, &sizing, up_to, tol)?.size)
}

/// The largest value up to `up_to` of the distance `size` frees (one
/// point of either law, or the scale of both) at which
/// [`fillet_edge_distances_variable`](crate::fillet_edge_distances_variable)
/// builds with these laws, as [`largest_law_radius`] reads a radius.
///
/// # Errors
///
/// As [`largest_law_radius`], for a point of either law.
#[allow(
    clippy::too_many_arguments,
    reason = "fillet_edge_distances_variable's own, with the size freed and its bound"
)]
pub fn largest_law_distances(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: &RadiusLaw,
    on_other: &RadiusLaw,
    profile: RoundProfile,
    size: DistanceLawSize,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let sizing = BlendSizing::DistanceLaws {
        edge,
        face,
        on_face,
        on_other,
        profile,
        size,
    };
    Ok(blend_limit(model, solid, &sizing, up_to, tol)?.size)
}

/// Whether `sizing`'s blend call builds at `size`, asked on a copy of the
/// model as the request asks it.
fn probe(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    size: f64,
    tol: Tolerances,
) -> bool {
    let mut scratch = model.clone();
    match *sizing {
        BlendSizing::Rounds { edges, corner } => {
            crate::fillet_edges_corner(&mut scratch, solid, edges, size, corner, tol).is_ok()
        }
        BlendSizing::RoundsCornerOverflow {
            edges,
            corner,
            overflow,
        } => crate::fillet_edges_corner_overflow(
            &mut scratch,
            solid,
            edges,
            size,
            corner,
            overflow,
            tol,
        )
        .is_ok(),
        BlendSizing::RoundDistances {
            edge,
            face,
            on_face,
            on_other,
            profile,
        } => crate::fillet_edge_distances(
            &mut scratch,
            solid,
            edge,
            face,
            size * on_face,
            size * on_other,
            profile,
            tol,
        )
        .is_ok(),
        BlendSizing::RoundDistancesOverflow {
            edge,
            face,
            on_face,
            on_other,
            profile,
            overflow,
        } => crate::fillet_edge_distances_overflow(
            &mut scratch,
            solid,
            edge,
            face,
            size * on_face,
            size * on_other,
            profile,
            overflow,
            tol,
        )
        .is_ok(),
        BlendSizing::Chamfer { edge, spec, corner } => {
            let specs = [(edge.clone(), scaled(spec, size))];
            crate::chamfer_edges_cornered(&mut scratch, solid, &specs, corner, tol).is_ok()
        }
        BlendSizing::RoundsOverflow { edges, overflow } => {
            crate::fillet_edges_overflow(&mut scratch, solid, edges, size, overflow, tol).is_ok()
        }
        BlendSizing::Chamfers {
            specs,
            corner,
            overflow,
        } => {
            let specs: Vec<(Shape, Chamfer)> = specs
                .iter()
                .map(|(edge, spec)| (edge.clone(), scaled(spec, size)))
                .collect();
            crate::chamfer_edges_overflow(&mut scratch, solid, &specs, corner, overflow, tol)
                .is_ok()
        }
        // Read by `law_seat_reading`, which probes them itself.
        BlendSizing::RoundLaw { .. }
        | BlendSizing::RoundLawOverflow { .. }
        | BlendSizing::RoundLinear { .. }
        | BlendSizing::DistanceLaws { .. } => false,
    }
}

/// `spec` at the measure `size`, as
/// [`largest_chamfer_distance`](crate::largest_chamfer_distance) reads it:
/// the distance, the scale of both distances, or the distance at the
/// angle.
pub(crate) fn scaled(spec: &Chamfer, size: f64) -> Chamfer {
    match spec {
        Chamfer::Symmetric(_) => Chamfer::Symmetric(size),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => Chamfer::Distances {
            face: face.clone(),
            on_face: size * on_face,
            on_other: size * on_other,
        },
        Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
            face: face.clone(),
            distance: size,
            angle: *angle,
        },
    }
}

/// The fractions of a chain a law's seat is read at: both ends and every
/// sixteenth between.
const LAW_STATIONS: u32 = 16;

/// What a law reading frees in one law: a value by its index, or the
/// scale of all of them.
#[derive(Debug, Clone, Copy)]
enum Free {
    Index(usize),
    Scale,
}

/// The value `at` names in `law`, by its index among the law's points.
fn free_point(law: &RadiusLaw, at: LawStation) -> OgeomResult<Free> {
    let n = match law {
        RadiusLaw::Linear { .. } => 2,
        RadiusLaw::Points(points) => points.len(),
    };
    if n == 0 {
        ogeom_bail!(Construction, "a radius law of no points states no radius");
    }
    let index = match at {
        LawStation::Start => 0,
        LawStation::End => n - 1,
        LawStation::Point(i) => i,
    };
    if index >= n {
        ogeom_bail!(
            Construction,
            "a law of {n} points has no point {index} to free"
        );
    }
    Ok(Free::Index(index))
}

/// `law` with the value `free` names set to `value`, or every value times
/// `value`; on a `closed` chain a point at the start and one at the end
/// are one point, set together.
fn law_at(
    law: &RadiusLaw,
    free: Option<Free>,
    value: f64,
    closed: bool,
    tol: Tolerances,
) -> RadiusLaw {
    let set = |points: &mut [(f64, f64)]| match free {
        None => {}
        Some(Free::Scale) => points.iter_mut().for_each(|p| p.1 *= value),
        Some(Free::Index(i)) => {
            let n = points.len();
            points[i].1 = value;
            let joined = closed
                && n > 1
                && points[0].0 <= tol.parametric()
                && points[n - 1].0 >= 1.0 - tol.parametric();
            if joined && (i == 0 || i == n - 1) {
                points[0].1 = value;
                points[n - 1].1 = value;
            }
        }
    };
    match law {
        RadiusLaw::Linear { start, end } => {
            let mut points = [(0.0, *start), (1.0, *end)];
            set(&mut points);
            RadiusLaw::Linear {
                start: points[0].1,
                end: points[1].1,
            }
        }
        RadiusLaw::Points(points) => {
            let mut points = points.clone();
            set(&mut points);
            RadiusLaw::Points(points)
        }
    }
}

/// The values a law states.
fn law_values(law: &RadiusLaw) -> Vec<f64> {
    match law {
        RadiusLaw::Linear { start, end } => vec![*start, *end],
        RadiusLaw::Points(points) => points.iter().map(|p| p.1).collect(),
    }
}

/// Whether `chain` closes: every vertex of its edges is shared by two
/// of them, a closed edge's one vertex by itself.
fn chain_closes(model: &Model, chain: &[Shape]) -> OgeomResult<bool> {
    let mut open: Vec<Shape> = Vec::new();
    for edge in chain {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(false);
        };
        for v in [a, b] {
            if let Some(k) = open.iter().position(|o| o.is_same(&v)) {
                open.swap_remove(k);
            } else {
                open.push(v);
            }
        }
    }
    Ok(!chain.is_empty() && open.is_empty())
}

/// The largest measure up to `up_to`, above `floor`, at which `fits`
/// holds, bisected to [`LIMIT_TOLERANCE`]; None where it holds at none.
fn seat_by(fits: impl Fn(f64) -> bool, up_to: f64, floor: f64) -> Option<f64> {
    if fits(up_to) {
        return Some(up_to);
    }
    let (mut lo, mut hi) = (floor * 2.0, up_to);
    if lo >= hi || !fits(lo) {
        return None;
    }
    while !close(lo, hi) {
        let mid = f64::midpoint(lo, hi);
        if fits(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(lo)
}

/// A law's seat: the largest size at which, at every station read, the
/// law's measure there stays within its contact's reach. Each station is
/// its fraction of the chain, the law it reads (an index into the laws
/// `laws` gives for a size) and its reach.
fn law_seat(
    stations: &[(f64, usize, f64)],
    laws: impl Fn(f64) -> Option<Vec<Profile>>,
    up_to: f64,
    floor: f64,
) -> Option<f64> {
    let fits = |size: f64| {
        laws(size).is_some_and(|laws| {
            stations
                .iter()
                .all(|&(s, k, reach)| laws[k].at(s).0 <= reach)
        })
    };
    seat_by(fits, up_to, floor)
}

/// The fractions [`LAW_STATIONS`] names.
fn law_fractions() -> Vec<f64> {
    (0..=LAW_STATIONS)
        .map(|k| f64::from(k) / f64::from(LAW_STATIONS))
        .collect()
}

/// The largest measure any value of `laws` takes up to the bound: `up_to`
/// for a freed value, the largest value times `up_to` for a scale.
fn measure_bound(laws: &[(&RadiusLaw, Option<Free>)], up_to: f64) -> f64 {
    laws.iter()
        .flat_map(|(law, free)| {
            let k = if matches!(free, Some(Free::Scale)) {
                up_to
            } else {
                1.0
            };
            law_values(law).into_iter().map(move |v| v * k)
        })
        .fold(up_to, f64::max)
}

/// The smallest size worth asking: a confusion distance, or for a scale
/// the scale at which the least value is one.
fn law_floor(laws: &[(&RadiusLaw, Option<Free>)], conf: f64) -> f64 {
    let least = laws
        .iter()
        .filter(|(_, free)| matches!(free, Some(Free::Scale)))
        .flat_map(|(law, _)| law_values(law))
        .fold(f64::INFINITY, f64::min);
    if least.is_finite() && least > conf {
        conf / least
    } else {
        conf
    }
}

/// `builds`, each call counted under test, so a test can tell how many
/// builds a reading asked.
fn counted(mut builds: impl FnMut(f64) -> bool) -> impl FnMut(f64) -> bool {
    move |size| {
        #[cfg(test)]
        tests::BUILDS.with(|n| n.set(n.get() + 1));
        builds(size)
    }
}

/// [`read_seat`] for a request under a law: the seat read at
/// [`LAW_STATIONS`] along the chain under the law at each size, the probe
/// the call under the law at the size asked.
fn law_seat_reading<R>(
    model: &Model,
    solid: &Shape,
    sizing: &BlendSizing<'_>,
    up_to: f64,
    tol: Tolerances,
    then: impl FnOnce(OgeomResult<Option<f64>>, f64, Builds<'_>) -> OgeomResult<R>,
) -> OgeomResult<R> {
    let conf = tol.confusion();
    match *sizing {
        BlendSizing::RoundLaw { chain, law, size } => {
            let overflow = Overflow::default();
            round_law_reading(model, solid, chain, law, size, overflow, up_to, tol, then)
        }
        BlendSizing::RoundLawOverflow {
            chain,
            law,
            size,
            overflow,
        } => round_law_reading(model, solid, chain, law, size, overflow, up_to, tol, then),
        BlendSizing::RoundLinear {
            edge,
            start,
            end,
            size,
        } => {
            // One edge, read as an open law even where it closes.
            let law = RadiusLaw::Linear { start, end };
            let chain = std::slice::from_ref(edge);
            let free = Some(match size {
                LawSize::At(at) => free_point(&law, at)?,
                LawSize::Scale => Free::Scale,
            });
            let at = |v: f64| law_at(&law, free, v, false, tol);
            Profile::new(&at(up_to), false, tol)?;
            let laws = [(&law, free)];
            let seat = round_law_seat(model, solid, chain, false, &laws, &at, up_to, tol);
            let builds = |v: f64| {
                let RadiusLaw::Linear { start, end } = at(v) else {
                    return false;
                };
                let mut scratch = model.clone();
                crate::fillet_edge_variable(&mut scratch, solid, edge, start, end, tol).is_ok()
            };
            then(seat, law_floor(&laws, conf), &mut counted(builds))
        }
        BlendSizing::DistanceLaws {
            edge,
            face,
            on_face,
            on_other,
            profile,
            size,
        } => {
            let closed = chain_closes(model, std::slice::from_ref(edge))?;
            let (free_face, free_other) = match size {
                DistanceLawSize::OnFace(at) => (Some(free_point(on_face, at)?), None),
                DistanceLawSize::OnOther(at) => (None, Some(free_point(on_other, at)?)),
                DistanceLawSize::Scale => (Some(Free::Scale), Some(Free::Scale)),
            };
            let at = |v: f64| {
                [
                    law_at(on_face, free_face, v, closed, tol),
                    law_at(on_other, free_other, v, closed, tol),
                ]
            };
            for law in at(up_to) {
                Profile::new(&law, closed, tol)?;
            }
            let laws = [(on_face, free_face), (on_other, free_other)];
            let seat = (|| -> OgeomResult<Option<f64>> {
                let reaches = crate::fits::station_reaches(
                    model,
                    solid,
                    edge,
                    &law_fractions(),
                    measure_bound(&laws, up_to),
                    |_, _| 1.0,
                    tol,
                )?;
                let stations: Vec<(f64, usize, f64)> = reaches
                    .into_iter()
                    .map(|r| (r.fraction, usize::from(!r.face.is_same(face)), r.largest))
                    .collect();
                let profiles = |v: f64| {
                    at(v)
                        .iter()
                        .map(|law| Profile::new(law, closed, tol).ok())
                        .collect::<Option<Vec<_>>>()
                };
                Ok(law_seat(&stations, profiles, up_to, law_floor(&laws, conf)))
            })();
            let builds = |v: f64| {
                let [a, b] = at(v);
                let mut scratch = model.clone();
                crate::fillet_edge_distances_variable(
                    &mut scratch,
                    solid,
                    edge,
                    face,
                    &a,
                    &b,
                    profile,
                    tol,
                )
                .is_ok()
            };
            then(seat, law_floor(&laws, conf), &mut counted(builds))
        }
        _ => ogeom_bail!(Construction, "not a request under a law"),
    }
}

/// [`law_seat_reading`] for [`BlendSizing::RoundLaw`] and
/// [`BlendSizing::RoundLawOverflow`]: each probe builds under `overflow`.
#[allow(
    clippy::too_many_arguments,
    reason = "the request's parts and law_seat_reading's own"
)]
fn round_law_reading<R>(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    size: LawSize,
    overflow: Overflow,
    up_to: f64,
    tol: Tolerances,
    then: impl FnOnce(OgeomResult<Option<f64>>, f64, Builds<'_>) -> OgeomResult<R>,
) -> OgeomResult<R> {
    let conf = tol.confusion();
    if chain.is_empty() {
        ogeom_bail!(
            Construction,
            "a variable radius over no edges has no largest size"
        );
    }
    let closed = chain_closes(model, chain)?;
    let free = Some(match size {
        LawSize::At(at) => free_point(law, at)?,
        LawSize::Scale => Free::Scale,
    });
    let at = |v: f64| law_at(law, free, v, closed, tol);
    Profile::new(&at(up_to), closed, tol)?;
    let laws = [(law, free)];
    let seat = round_law_seat(model, solid, chain, closed, &laws, &at, up_to, tol);
    let builds = |v: f64| {
        let mut scratch = model.clone();
        crate::fillet_edges_variable_overflow(&mut scratch, solid, chain, &at(v), overflow, tol)
            .is_ok()
    };
    then(seat, law_floor(&laws, conf), &mut counted(builds))
}

/// A round's seat under a law: on one edge its contacts' reach at each
/// station; on a chain of several, where the stations' places along the
/// chain are not read, the constant round's seat held at every fraction.
#[allow(
    clippy::too_many_arguments,
    reason = "the request's parts, as law_seat_reading has them"
)]
fn round_law_seat(
    model: &Model,
    solid: &Shape,
    chain: &[Shape],
    closed: bool,
    laws: &[(&RadiusLaw, Option<Free>)],
    at: &impl Fn(f64) -> RadiusLaw,
    up_to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let bound = measure_bound(laws, up_to);
    let stations: Vec<(f64, usize, f64)> = if let [edge] = chain {
        crate::fits::station_reaches(
            model,
            solid,
            edge,
            &law_fractions(),
            bound,
            |_, factor| factor,
            tol,
        )?
        .into_iter()
        .map(|r| (r.fraction, 0, r.largest))
        .collect()
    } else {
        let Some(seat) = crate::fits::seat_fillet_radius(model, solid, chain, bound, tol)? else {
            return Ok(None);
        };
        law_fractions().into_iter().map(|s| (s, 0, seat)).collect()
    };
    let profiles = |v: f64| Profile::new(&at(v), closed, tol).ok().map(|p| vec![p]);
    Ok(law_seat(
        &stations,
        profiles,
        up_to,
        law_floor(laws, tol.confusion()),
    ))
}

/// The probes one reading has asked, each size once.
struct Probes<F> {
    builds: F,
    asked: Vec<(f64, bool)>,
}

impl<F: FnMut(f64) -> bool> Probes<F> {
    /// Whether the call builds at `size`; None once [`LIMIT_PROBES`] are
    /// spent.
    fn ask(&mut self, size: f64) -> Option<bool> {
        if let Some(&(_, built)) = self.asked.iter().find(|(s, _)| *s == size) {
            return Some(built);
        }
        if self.spent() >= LIMIT_PROBES {
            return None;
        }
        let built = (self.builds)(size);
        self.asked.push((size, built));
        Some(built)
    }

    fn spent(&self) -> u32 {
        u32::try_from(self.asked.len()).unwrap_or(u32::MAX)
    }
}

/// Whether `hi` stands within [`LIMIT_TOLERANCE`] of `lo`.
fn close(lo: f64, hi: f64) -> bool {
    hi <= lo * (1.0 + LIMIT_TOLERANCE)
}

/// Bisect between `lo`, which built, and `hi`, which was refused, until
/// they are [`close`] or the probes are spent.
fn bisect<F: FnMut(f64) -> bool>(probes: &mut Probes<F>, mut lo: f64, mut hi: f64) -> (f64, f64) {
    while !close(lo, hi) {
        let mid = f64::midpoint(lo, hi);
        match probes.ask(mid) {
            Some(true) => lo = mid,
            Some(false) => hi = mid,
            None => break,
        }
    }
    (lo, hi)
}

/// Search down from `hi`, refused, halving until a size above `floor`
/// builds, then bisect: the boundary below `hi`, or None with the smallest
/// size refused where none builds.
fn below<F: FnMut(f64) -> bool>(probes: &mut Probes<F>, hi: f64, floor: f64) -> (Option<f64>, f64) {
    let mut hi = hi;
    let mut t = hi * 0.5;
    loop {
        if t <= floor {
            return (None, hi);
        }
        match probes.ask(t) {
            Some(true) => break,
            Some(false) => {
                hi = t;
                t *= 0.5;
            }
            None => return (None, hi),
        }
    }
    let (lo, hi) = bisect(probes, t, hi);
    (Some(lo), hi)
}

/// The shared confirmation every largest-size reading ends in: the
/// boundary between sizes that build and sizes refused, by `builds`,
/// started from the `seat` reading where there is one, sizes at or below
/// `floor` never asked. See [`blend_limit`] for the search.
pub(crate) fn confirm(
    seat: Option<Option<f64>>,
    up_to: f64,
    floor: f64,
    builds: impl FnMut(f64) -> bool,
) -> BlendLimit {
    let mut probes = Probes {
        builds,
        asked: Vec::new(),
    };
    let seat_size = seat.flatten();
    // Where the seat holds no size, the smallest worth asking.
    let start = match seat {
        Some(None) => (floor * 2.0).min(up_to),
        Some(Some(s)) => s.min(up_to).max(floor * 2.0).min(up_to),
        None => up_to,
    };
    let (mut lo, mut hi): (Option<f64>, Option<f64>) = match probes.ask(start) {
        Some(true) => {
            // The seat can stand short of what builds (a bevel notching
            // past its face's far side): step just past it, then double.
            let mut lo = start;
            let mut hi = None;
            let mut t = (start * (1.0 + LIMIT_TOLERANCE)).min(up_to);
            while lo < up_to {
                match probes.ask(t) {
                    Some(true) => {
                        lo = t;
                        t = (t * 2.0).min(up_to);
                    }
                    Some(false) => {
                        hi = Some(t);
                        break;
                    }
                    None => break,
                }
            }
            match hi {
                Some(h) => {
                    let (l, h) = bisect(&mut probes, lo, h);
                    (Some(l), Some(h))
                }
                None => (Some(lo), None),
            }
        }
        Some(false) => {
            let (l, h) = below(&mut probes, start, floor);
            (l, Some(h))
        }
        None => (None, None),
    };
    let mut builds_again = None;
    if let Some(l) = lo {
        // A refusal below the answer: the sizes that build break there.
        let samples: Vec<f64> = (1..=LIMIT_SAMPLES)
            .map(|k| l * f64::from(k) / f64::from(LIMIT_SAMPLES + 1))
            .filter(|s| *s > floor)
            .collect();
        let mut first_refused = None;
        let mut built_below = None;
        for s in samples {
            match probes.ask(s) {
                Some(true) if first_refused.is_none() => built_below = Some(s),
                Some(false) if first_refused.is_none() => first_refused = Some(s),
                _ => {}
            }
        }
        if let Some(r) = first_refused {
            builds_again = Some(l);
            let (l2, h2) = match built_below {
                Some(b) => {
                    let (l2, h2) = bisect(&mut probes, b, r);
                    (Some(l2), h2)
                }
                None => below(&mut probes, r, floor),
            };
            lo = l2;
            hi = Some(h2);
        }
    }
    if let (Some(h), None) = (hi, builds_again) {
        // A size building again above the first refusal.
        if h < up_to {
            for k in 1..=LIMIT_SAMPLES {
                let s = h * (up_to / h).powf(f64::from(k) / f64::from(LIMIT_SAMPLES));
                if probes.ask(s) == Some(true) {
                    builds_again = Some(s);
                    break;
                }
            }
        }
    }
    BlendLimit {
        size: lo,
        refused: hi,
        builds_again,
        seat: seat_size,
        probes: probes.spent(),
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_geom::Curve3d as _;
    use ogeom_math::{Direction, Frame, Point};
    use ogeom_topo::{ShapeType, explore_unique};
    use std::cell::Cell;

    thread_local! {
        /// The builds asked on this thread, counted by [`counted`].
        pub(super) static BUILDS: Cell<u32> = const { Cell::new(0) };
    }

    const T: Tolerances = Tolerances::millimetres();

    /// A plate 20 wide (x), 4 deep (y), 20 tall (z), with its top front
    /// edge (y = 0, z = 20), its top face and its front face.
    fn plate() -> (Model, Shape, Shape, Shape, Shape) {
        let mut m = Model::new();
        let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
        let s = ogeom_algo::make_box(&mut m, frame, (20.0, 4.0, 20.0), T)
            .unwrap()
            .shape;
        let edge = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
                [a, b].iter().all(|&t| {
                    let p = c.point_at(t, T).unwrap();
                    p.y.abs() < 1e-9 && (p.z - 20.0).abs() < 1e-9
                })
            })
            .unwrap();
        let face_at = |on: &dyn Fn(Point) -> bool| {
            explore_unique(&m, &s, ShapeType::Face)
                .unwrap()
                .into_iter()
                .find(|f| {
                    let b = ogeom_algo::shape_bounds(&m, f, T).unwrap();
                    on(b.low().unwrap()) && on(b.high().unwrap())
                })
                .unwrap()
        };
        let top = face_at(&|p| (p.z - 20.0).abs() < 1e-6);
        let front = face_at(&|p| p.y.abs() < 1e-6);
        (m, s, edge, top, front)
    }

    /// `blend_seat` on `sizing`, checked to build nothing and to be
    /// `blend_limit`'s seat, which is checked to count its builds.
    fn seat_agrees(m: &Model, s: &Shape, sizing: &BlendSizing<'_>, up_to: f64) -> Option<f64> {
        BUILDS.with(|n| n.set(0));
        let seat = blend_seat(m, s, sizing, up_to, T).unwrap();
        assert_eq!(BUILDS.with(Cell::get), 0, "the seat reading built");
        let limit = blend_limit(m, s, sizing, up_to, T).unwrap();
        assert_eq!(BUILDS.with(Cell::get), limit.probes);
        assert!(limit.probes > 0);
        assert_eq!(seat, limit.seat);
        seat
    }

    #[test]
    fn a_round_s_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let edges = [edge];
        let sizing = BlendSizing::Rounds {
            edges: &edges,
            corner: FilletCorner::Rolling,
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn a_round_s_seat_under_overflow_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let edges = [edge];
        let sizing = BlendSizing::RoundsOverflow {
            edges: &edges,
            overflow: Overflow::default(),
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn seats_under_overflow_are_the_plain_seats_without_a_build() {
        let (m, s, edge, top, _) = plate();
        let edges = [edge.clone()];
        let law = RadiusLaw::Linear {
            start: 1.0,
            end: 2.0,
        };
        let off = Overflow::none();
        let pairs = [
            (
                BlendSizing::Rounds {
                    edges: &edges,
                    corner: FilletCorner::Rolling,
                },
                BlendSizing::RoundsCornerOverflow {
                    edges: &edges,
                    corner: FilletCorner::Rolling,
                    overflow: off,
                },
            ),
            (
                BlendSizing::RoundDistances {
                    edge: &edge,
                    face: &top,
                    on_face: 1.0,
                    on_other: 0.5,
                    profile: RoundProfile::Elliptical,
                },
                BlendSizing::RoundDistancesOverflow {
                    edge: &edge,
                    face: &top,
                    on_face: 1.0,
                    on_other: 0.5,
                    profile: RoundProfile::Elliptical,
                    overflow: off,
                },
            ),
            (
                BlendSizing::RoundLaw {
                    chain: &edges,
                    law: &law,
                    size: LawSize::Scale,
                },
                BlendSizing::RoundLawOverflow {
                    chain: &edges,
                    law: &law,
                    size: LawSize::Scale,
                    overflow: off,
                },
            ),
        ];
        for (plain, under) in &pairs {
            let seat = seat_agrees(&m, &s, under, 40.0);
            assert!(seat.is_some());
            assert_eq!(seat, blend_seat(&m, &s, plain, 40.0, T).unwrap());
        }
    }

    #[test]
    fn each_face_s_reach_is_read_without_a_build() {
        let (m, s, edge, top, front) = plate();
        for (face, reach) in [(&top, 4.0), (&front, 20.0)] {
            let sizing = BlendSizing::RoundDistances {
                edge: &edge,
                face,
                on_face: 1.0,
                on_other: 1e-5,
                profile: RoundProfile::Elliptical,
            };
            let seat = seat_agrees(&m, &s, &sizing, 40.0).unwrap();
            assert!((seat - reach).abs() < 1e-6, "{seat} against {reach}");
        }
    }

    #[test]
    fn a_bevel_s_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let sizing = BlendSizing::Chamfer {
            edge: &edge,
            spec: &Chamfer::Symmetric(1.0),
            corner: ChamferCorner::Point,
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn bevels_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let specs = [(edge, Chamfer::Symmetric(1.0))];
        let sizing = BlendSizing::Chamfers {
            specs: &specs,
            corner: ChamferCorner::Point,
            overflow: Overflow::default(),
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn a_law_s_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let chain = [edge];
        let law = RadiusLaw::Linear {
            start: 1.0,
            end: 2.0,
        };
        let sizing = BlendSizing::RoundLaw {
            chain: &chain,
            law: &law,
            size: LawSize::Scale,
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn a_linear_round_s_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, ..) = plate();
        let sizing = BlendSizing::RoundLinear {
            edge: &edge,
            start: 1.0,
            end: 2.0,
            size: LawSize::At(LawStation::End),
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn distance_laws_seat_is_the_limit_s_without_a_build() {
        let (m, s, edge, top, _) = plate();
        let on_face = RadiusLaw::Linear {
            start: 1.0,
            end: 2.0,
        };
        let on_other = RadiusLaw::Linear {
            start: 1.0,
            end: 1.0,
        };
        let sizing = BlendSizing::DistanceLaws {
            edge: &edge,
            face: &top,
            on_face: &on_face,
            on_other: &on_other,
            profile: RoundProfile::Elliptical,
            size: DistanceLawSize::Scale,
        };
        assert!(seat_agrees(&m, &s, &sizing, 40.0).is_some());
    }

    #[test]
    fn an_unusable_bound_is_refused_as_the_limit_refuses_it() {
        let (m, s, edge, ..) = plate();
        let edges = [edge];
        let sizing = BlendSizing::Rounds {
            edges: &edges,
            corner: FilletCorner::Setback,
        };
        let seat = blend_seat(&m, &s, &sizing, 0.0, T).unwrap_err();
        let limit = blend_limit(&m, &s, &sizing, 0.0, T).unwrap_err();
        assert_eq!(seat.to_string(), limit.to_string());
    }

    fn run(seat: Option<Option<f64>>, up_to: f64, builds: impl Fn(f64) -> bool) -> BlendLimit {
        confirm(seat, up_to, 1e-7, builds)
    }

    #[test]
    fn a_closed_chain_s_end_points_move_together_and_an_open_one_s_do_not() {
        let tol = Tolerances::millimetres();
        let law = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 2.0), (1.0, 1.0)]);
        let end = Some(Free::Index(2));
        assert_eq!(
            law_at(&law, end, 3.0, true, tol),
            RadiusLaw::Points(vec![(0.0, 3.0), (0.5, 2.0), (1.0, 3.0)])
        );
        assert_eq!(
            law_at(&law, end, 3.0, false, tol),
            RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 2.0), (1.0, 3.0)])
        );
        // A law short of the chain's end has no point there to join.
        let short = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 2.0)]);
        assert_eq!(
            law_at(&short, Some(Free::Index(1)), 3.0, true, tol),
            RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 3.0)])
        );
        assert_eq!(
            law_at(
                &RadiusLaw::Linear {
                    start: 1.0,
                    end: 2.0
                },
                Some(Free::Scale),
                3.0,
                false,
                tol
            ),
            RadiusLaw::Linear {
                start: 3.0,
                end: 6.0
            }
        );
    }

    #[test]
    fn an_exact_seat_is_confirmed_in_two_probes_and_the_samples() {
        let found = run(Some(Some(2.0)), 10.0, |s| s <= 2.0);
        assert_eq!(found.size, Some(2.0));
        assert!(close(2.0, found.refused.unwrap()));
        assert_eq!(found.builds_again, None);
        assert_eq!(found.probes, 2 + 2 * LIMIT_SAMPLES);
    }

    #[test]
    fn a_seat_short_of_what_builds_is_searched_past() {
        let found = run(Some(Some(2.0)), 100.0, |s| s <= 7.3);
        let size = found.size.unwrap();
        assert!(size <= 7.3 && close(size, 7.3), "{size}");
        assert!(found.refused.unwrap() > 7.3);
    }

    #[test]
    fn a_seat_past_what_builds_is_searched_below() {
        let found = run(Some(Some(2.0)), 100.0, |s| s <= 0.31);
        let size = found.size.unwrap();
        assert!(size <= 0.31 && close(size, 0.31), "{size}");
    }

    #[test]
    fn the_bound_building_answers_the_bound() {
        let found = run(Some(Some(0.5)), 0.5, |_| true);
        assert_eq!(found.size, Some(0.5));
        assert_eq!(found.refused, None);
    }

    #[test]
    fn nothing_building_answers_none_within_the_cap() {
        let found = run(None, 100.0, |_| false);
        assert_eq!(found.size, None);
        assert!(found.probes <= LIMIT_PROBES);
    }

    #[test]
    fn a_size_building_again_above_a_refusal_is_reported() {
        // Refused only between 2 and 3: the answer is the first boundary.
        let found = run(Some(Some(2.0)), 10.0, |s| !(s > 2.0 && s < 3.0));
        assert_eq!(found.size, Some(2.0));
        assert!(found.builds_again.unwrap() >= 3.0);
    }

    #[test]
    fn a_refusal_below_the_seat_s_boundary_moves_the_answer_down() {
        // Refused between 0.9 and 1.1, building again up to 2: the answer
        // is 0.9, and 2 is reported as building again.
        let found = run(Some(Some(2.0)), 10.0, |s| {
            s <= 0.9 || (1.1..=2.0).contains(&s)
        });
        let size = found.size.unwrap();
        assert!(size <= 0.9 && close(size, 0.9), "{size}");
        assert_eq!(found.builds_again, Some(2.0));
    }
}
