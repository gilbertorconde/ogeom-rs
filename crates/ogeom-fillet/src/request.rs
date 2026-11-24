//! One edge-blend request: every edge with its section and the options of
//! the operation, read as one set (ADR 0001, section 3, and amendment "the
//! public request"; docs/design/1265-fillet-request.md).
//!
//! The request picks the existing call that builds what it asks and hands
//! it the arguments that call takes, so each refusal text is that call's
//! own. A combination no call builds is refused by name.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_err};
use ogeom_topo::{Model, Shape};

use crate::{
    Blend, BlendLimit, BlendRefusal, BlendSizing, Chamfer, ChamferCorner, DistanceLawSize,
    FilletCorner, LawSize, LawStation, Overflow, RadiusLaw, RoundProfile, Setbacks,
};

/// Every edge blend of one operation, read as one set (ADR 0001, section 3).
/// The order edges are added in carries no meaning for the shape built; it
/// only numbers them for [`BlendRefusal::edges`], chains flattened into
/// their edges.
///
/// Each option is a method, and each has the default of the call without
/// it: setback round corners, facet chamfer corners, no setback corners at
/// vertices, every [`Overflow`] behaviour on. A corner option names no
/// blend of a request without that kind (the round corner of bevels alone)
/// and is not read there.
///
/// ```
/// # use ogeom_core::Tolerances;
/// # use ogeom_fillet::BlendRequest;
/// # use ogeom_math::Frame;
/// # use ogeom_topo::{Model, ShapeType, explore_unique};
/// let tol = Tolerances::millimetres();
/// let mut model = Model::new();
/// let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), tol)?.shape;
/// let edges = explore_unique(&model, &block, ShapeType::Edge)?;
/// let rounded = BlendRequest::new()
///     .round(edges[..1].iter().cloned(), 2.0)
///     .build(&mut model, &block, tol)?;
/// # let _ = rounded;
/// # Ok::<(), ogeom_core::OgeomError>(())
/// ```
#[derive(Debug, Clone, Default)]
#[must_use]
pub struct BlendRequest {
    entries: Vec<Entry>,
    round_corner: FilletCorner,
    chamfer_corner: ChamferCorner,
    setbacks: Vec<(Shape, Setbacks)>,
    overflow: Overflow,
}

/// What one blend sweeps along its edge or chain.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum BlendSection {
    /// A circular round: the rolling ball of this radius.
    Round(f64),
    /// A round whose radius follows the law along the chain.
    RoundLaw(RadiusLaw),
    /// A round meeting `face` `on_face` from the edge and the other face
    /// `on_other` from it, its section `profile` between the contacts.
    RoundDistances {
        /// The face the first distance runs along.
        face: Shape,
        /// The distance along `face`.
        on_face: f64,
        /// The distance along the edge's other face.
        on_other: f64,
        /// The section between the two contacts.
        profile: RoundProfile,
    },
    /// As [`BlendSection::RoundDistances`], each distance following its law
    /// along the chain.
    DistanceLaws {
        /// The face the first law's distance runs along.
        face: Shape,
        /// The distances along `face`.
        on_face: RadiusLaw,
        /// The distances along the other face.
        on_other: RadiusLaw,
        /// The section between the two contacts.
        profile: RoundProfile,
    },
    /// A flat bevel: symmetric, two distances, or a distance and an angle.
    Chamfer(Chamfer),
}

impl From<Chamfer> for BlendSection {
    fn from(chamfer: Chamfer) -> Self {
        Self::Chamfer(chamfer)
    }
}

impl From<Blend> for BlendSection {
    fn from(blend: Blend) -> Self {
        match blend {
            Blend::Round(radius) => Self::Round(radius),
            Blend::Chamfer(chamfer) => Self::Chamfer(chamfer),
            Blend::RoundDistances {
                face,
                on_face,
                on_other,
                profile,
            } => Self::RoundDistances {
                face,
                on_face,
                on_other,
                profile,
            },
        }
    }
}

/// Which size a largest-size reading of a [`BlendRequest`] varies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Vary {
    /// Every length of the request times one factor; angles, profiles and
    /// the shape of each law held. The reading is the factor: a round asked
    /// at radius `1.0` reads its largest radius.
    #[default]
    Scale,
    /// One station of the request's round law, the rest held. The reading
    /// is that station's radius.
    Law(LawStation),
    /// One station of the `on_face` law of a
    /// [`BlendSection::DistanceLaws`] section, the rest held.
    OnFace(LawStation),
    /// One station of its `on_other` law, the rest held.
    OnOther(LawStation),
}

/// Edges added by one builder call, under one section.
#[derive(Debug, Clone)]
struct Entry {
    edges: Vec<Shape>,
    section: BlendSection,
}

/// A request read as the call that builds it, with that call's arguments.
enum Lowered<'a> {
    /// Rounds of one radius.
    Rounds { edges: Vec<Shape>, radius: f64 },
    /// Rounds of several radii.
    RoundsWith(Vec<(Shape, f64)>),
    /// Bevels alone.
    Chamfers(Vec<(Shape, Chamfer)>),
    /// One round of two distances alone.
    Distances {
        edge: Shape,
        face: &'a Shape,
        on_face: f64,
        on_other: f64,
        profile: RoundProfile,
    },
    /// Constant sections of more than one kind.
    Mixed(Vec<(Shape, Blend)>),
    /// One chain under a round law.
    RoundLaw {
        chain: &'a [Shape],
        law: &'a RadiusLaw,
    },
    /// One chain under a pair of distance laws.
    DistanceLaws {
        chain: &'a [Shape],
        face: &'a Shape,
        on_face: &'a RadiusLaw,
        on_other: &'a RadiusLaw,
        profile: RoundProfile,
    },
}

impl Lowered<'_> {
    /// What the request asks, as a refusal names it.
    const fn what(&self) -> &'static str {
        match self {
            Self::Rounds { .. } => "rounds",
            Self::RoundsWith(_) => "rounds of several radii",
            Self::Chamfers(_) => "bevels",
            Self::Distances { .. } => "a round of two distances",
            Self::Mixed(_) => "blends of several kinds",
            Self::RoundLaw { .. } => "a round law",
            Self::DistanceLaws { .. } => "distance laws",
        }
    }
}

fn refused(text: String) -> OgeomError {
    ogeom_err!(Construction, "{text}")
}

/// `x` read back at the request's own measure: the call read sizes of the
/// measure `unit`.
fn per(x: Option<f64>, unit: f64) -> Option<f64> {
    x.map(|x| x / unit)
}

impl BlendRequest {
    /// A request of no blends, every option at its default.
    pub fn new() -> Self {
        Self::default()
    }

    /// Each edge rounded at `radius`.
    pub fn round(self, edges: impl IntoIterator<Item = Shape>, radius: f64) -> Self {
        self.chain(edges, BlendSection::Round(radius))
    }

    /// Each edge bevelled `distance` back along both faces.
    pub fn chamfer(self, edges: impl IntoIterator<Item = Shape>, distance: f64) -> Self {
        self.chain(edges, Chamfer::Symmetric(distance))
    }

    /// One edge with any section; a law on it is a chain of one edge.
    pub fn edge(self, edge: Shape, section: impl Into<BlendSection>) -> Self {
        self.chain([edge], section)
    }

    /// A chain under one section. A law runs along the whole chain, as
    /// [`fillet_edges_variable`](crate::fillet_edges_variable) takes it. A
    /// constant section on a chain is that section on each edge.
    pub fn chain(
        mut self,
        chain: impl IntoIterator<Item = Shape>,
        section: impl Into<BlendSection>,
    ) -> Self {
        self.entries.push(Entry {
            edges: chain.into_iter().collect(),
            section: section.into(),
        });
        self
    }

    /// How a vertex where rounds meet closes. Default
    /// [`FilletCorner::Setback`].
    pub const fn round_corner(mut self, corner: FilletCorner) -> Self {
        self.round_corner = corner;
        self
    }

    /// How a vertex where three or more bevels meet closes. Default
    /// [`ChamferCorner::Facet`].
    pub const fn chamfer_corner(mut self, corner: ChamferCorner) -> Self {
        self.chamfer_corner = corner;
        self
    }

    /// A vertex closed as a setback corner with these distances (ADR 0001,
    /// section 13, question 2), as
    /// [`fillet_edges_cornered`](crate::fillet_edges_cornered) closes it.
    pub fn setback(mut self, vertex: Shape, setbacks: Setbacks) -> Self {
        self.setbacks.push((vertex, setbacks));
        self
    }

    /// What a blend reaching other geometry does. Default: every behaviour
    /// on.
    pub const fn overflow(mut self, overflow: Overflow) -> Self {
        self.overflow = overflow;
        self
    }

    /// The request built on `solid`.
    ///
    /// | Request | Built as |
    /// |---|---|
    /// | rounds of one radius | [`fillet_edges_overflow`](crate::fillet_edges_overflow); [`fillet_edges_corner_overflow`](crate::fillet_edges_corner_overflow) when rolling; [`fillet_edges_cornered`](crate::fillet_edges_cornered) with setbacks |
    /// | rounds of several radii | [`fillet_edges_with`](crate::fillet_edges_with); [`blend_edges_overflow`](crate::blend_edges_overflow) under other options |
    /// | bevels alone | [`chamfer_edges_with`](crate::chamfer_edges_with) at the defaults; otherwise [`chamfer_edges_overflow`](crate::chamfer_edges_overflow) |
    /// | one round of two distances | [`fillet_edge_distances_overflow`](crate::fillet_edge_distances_overflow) |
    /// | constant sections of several kinds | [`blend_edges_overflow`](crate::blend_edges_overflow) |
    /// | one chain of a round law | [`fillet_edges_variable_overflow`](crate::fillet_edges_variable_overflow) |
    /// | one chain of distance laws | [`fillet_edges_distances_variable`](crate::fillet_edges_distances_variable) |
    ///
    /// A request with no edges is read by the kind of its first builder
    /// call (rounds where there is none) and takes that call's refusal.
    ///
    /// # Errors
    ///
    /// The call's own refusals; and
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction),
    /// by name, for a combination no call builds: a law with other blends
    /// or under rolling corners, setback corners with anything but rounds
    /// of one radius under the default overflow, rolling corners with
    /// anything but rounds of one radius, distance laws under overflow
    /// options.
    pub fn build(&self, model: &mut Model, solid: &Shape, tol: Tolerances) -> OgeomResult<Built> {
        let lowered = self.lower()?;
        let what = lowered.what();
        let rolling = self.round_corner == FilletCorner::Rolling;
        let defaults = self.overflow == Overflow::default();
        if !self.setbacks.is_empty() && !matches!(lowered, Lowered::Rounds { .. }) {
            return Err(refused(format!(
                "setback corners at vertices on {what} are not built"
            )));
        }
        if rolling && !matches!(lowered, Lowered::Rounds { .. } | Lowered::Chamfers(_)) {
            return Err(refused(format!("rolling corners on {what} are not built")));
        }
        match lowered {
            Lowered::Rounds { edges, radius } => {
                if !self.setbacks.is_empty() {
                    if rolling {
                        return Err(refused(
                            "setback corners at vertices with rolling corners are not built".into(),
                        ));
                    }
                    if !defaults {
                        return Err(refused(
                            "setback corners at vertices under overflow options are not built"
                                .into(),
                        ));
                    }
                    return crate::fillet_edges_cornered(
                        model,
                        solid,
                        &edges,
                        radius,
                        &self.setbacks,
                        tol,
                    );
                }
                if rolling {
                    crate::fillet::rounds_corner(
                        model,
                        solid,
                        &edges,
                        radius,
                        FilletCorner::Rolling,
                        self.overflow,
                        tol,
                    )
                } else {
                    crate::fillet::rounds(model, solid, &edges, radius, self.overflow, tol)
                }
            }
            Lowered::RoundsWith(rounds) => {
                if defaults {
                    crate::fillet_edges_with(model, solid, &rounds, tol)
                } else {
                    let pairs: Vec<(Shape, Blend)> = rounds
                        .into_iter()
                        .map(|(edge, radius)| (edge, Blend::Round(radius)))
                        .collect();
                    crate::blend_edges_overflow(
                        model,
                        solid,
                        &pairs,
                        self.chamfer_corner,
                        self.overflow,
                        tol,
                    )
                }
            }
            Lowered::Chamfers(specs) => {
                if defaults && self.chamfer_corner == ChamferCorner::Facet {
                    crate::chamfer::bevels(model, solid, &specs, tol)
                } else {
                    crate::chamfer_overflow::bevels_overflow(
                        model,
                        solid,
                        &specs,
                        self.chamfer_corner,
                        self.overflow,
                        tol,
                    )
                }
            }
            Lowered::Distances {
                edge,
                face,
                on_face,
                on_other,
                profile,
            } => crate::asym::distance_round(
                model,
                solid,
                &edge,
                face,
                on_face,
                on_other,
                profile,
                self.overflow,
                tol,
            ),
            Lowered::Mixed(pairs) => crate::blend_edges_overflow(
                model,
                solid,
                &pairs,
                self.chamfer_corner,
                self.overflow,
                tol,
            ),
            Lowered::RoundLaw { chain, law } => {
                crate::fillet::law_rounds(model, solid, chain, law, self.overflow, tol)
            }
            Lowered::DistanceLaws {
                chain,
                face,
                on_face,
                on_other,
                profile,
            } => {
                if !defaults {
                    return Err(refused(
                        "distance laws under overflow options are not built".into(),
                    ));
                }
                crate::asym_law::distance_law_rounds(
                    model, solid, chain, face, on_face, on_other, profile, tol,
                )
            }
        }
    }

    /// [`build`](Self::build), its refusal naming every refused edge, read
    /// in request order, so that the request without them builds (ADR 0001,
    /// amendment "checked blends under overflow"): as
    /// [`fillet_edges_corner_overflow_checked`](crate::fillet_edges_corner_overflow_checked)
    /// for rounds of one radius, and
    /// [`chamfer_edges_overflow_checked`](crate::chamfer_edges_overflow_checked)
    /// for bevels alone, under the request's corner and overflow.
    ///
    /// # Errors
    ///
    /// The checked call's refusal; and, naming no edge, an
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// by name for any other request, or one with setback corners at
    /// vertices.
    pub fn build_checked(
        &self,
        model: &mut Model,
        solid: &Shape,
        tol: Tolerances,
    ) -> Result<Built, BlendRefusal> {
        let whole = |error: OgeomError| BlendRefusal {
            edges: Vec::new(),
            error,
        };
        let lowered = self.lower().map_err(whole)?;
        if !self.setbacks.is_empty() {
            return Err(whole(refused(
                "a checked reading with setback corners at vertices is not built".into(),
            )));
        }
        match lowered {
            Lowered::Rounds { edges, radius } => crate::blame::rounds_checked(
                model,
                solid,
                &edges,
                radius,
                self.round_corner,
                self.overflow,
                tol,
            ),
            Lowered::Chamfers(specs) => crate::blame::bevels_checked(
                model,
                solid,
                &specs,
                self.chamfer_corner,
                self.overflow,
                tol,
            ),
            other => Err(whole(refused(format!(
                "a checked reading of {} is not built",
                other.what()
            )))),
        }
    }

    /// [`build`](Self::build) of the largest subset that builds, with each
    /// refused edge's refusal: as
    /// [`fillet_edges_partial`](crate::fillet_edges_partial) for rounds of
    /// one radius and
    /// [`chamfer_edges_partial`](crate::chamfer_edges_partial) for bevels
    /// alone, at the default options.
    ///
    /// # Errors
    ///
    /// The partial call's own; and
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// by name for any other request, or one with an option other than its
    /// default.
    pub fn build_partial(
        &self,
        model: &mut Model,
        solid: &Shape,
        tol: Tolerances,
    ) -> OgeomResult<(Built, Vec<BlendRefusal>)> {
        let lowered = self.lower()?;
        if !self.setbacks.is_empty() {
            return Err(refused(
                "a partial reading with setback corners at vertices is not built".into(),
            ));
        }
        if self.overflow != Overflow::default() {
            return Err(refused(
                "a partial reading under overflow options is not built".into(),
            ));
        }
        match lowered {
            Lowered::Rounds { edges, radius } => {
                if self.round_corner == FilletCorner::Rolling {
                    return Err(refused(
                        "a partial reading with rolling corners is not built".into(),
                    ));
                }
                crate::blame::rounds_partial(model, solid, &edges, radius, tol)
            }
            Lowered::Chamfers(specs) => {
                if self.chamfer_corner != ChamferCorner::Facet {
                    return Err(refused(
                        "a partial reading with point corners is not built".into(),
                    ));
                }
                crate::blame::bevels_partial(model, solid, &specs, tol)
            }
            other => Err(refused(format!(
                "a partial reading of {} is not built",
                other.what()
            ))),
        }
    }

    /// The largest size up to `up_to` at which the request builds, varying
    /// what `vary` names, with how it was confirmed: the reading of
    /// [`blend_limit`](crate::blend_limit), every size in the request's own
    /// measure (under [`Vary::Scale`], the factor).
    ///
    /// # Errors
    ///
    /// The reading's own refusals; and
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// by name for a request no reading takes: rounds of several radii,
    /// blends of several kinds, setback corners at vertices, bevels whose
    /// measures differ, a law station of a request with no such law.
    pub fn limit(
        &self,
        model: &Model,
        solid: &Shape,
        vary: Vary,
        up_to: f64,
        tol: Tolerances,
    ) -> OgeomResult<BlendLimit> {
        self.reading(vary, |sizing, unit| {
            let limit = crate::limit::limit_of(model, solid, sizing, up_to * unit, tol)?;
            Ok(BlendLimit {
                size: per(limit.size, unit),
                refused: per(limit.refused, unit),
                builds_again: per(limit.builds_again, unit),
                seat: per(limit.seat, unit),
                probes: limit.probes,
            })
        })
    }

    /// [`limit`](Self::limit)'s size alone.
    ///
    /// # Errors
    ///
    /// As [`limit`](Self::limit).
    pub fn largest(
        &self,
        model: &Model,
        solid: &Shape,
        vary: Vary,
        up_to: f64,
        tol: Tolerances,
    ) -> OgeomResult<Option<f64>> {
        Ok(self.limit(model, solid, vary, up_to, tol)?.size)
    }

    /// The seat reading [`limit`](Self::limit) starts from, without
    /// building: [`blend_seat`](crate::blend_seat).
    ///
    /// # Errors
    ///
    /// As [`limit`](Self::limit).
    pub fn seat(
        &self,
        model: &Model,
        solid: &Shape,
        vary: Vary,
        up_to: f64,
        tol: Tolerances,
    ) -> OgeomResult<Option<f64>> {
        self.reading(vary, |sizing, unit| {
            Ok(per(
                crate::limit::seat_of(model, solid, sizing, up_to * unit, tol)?,
                unit,
            ))
        })
    }

    /// The request read as the [`BlendSizing`] its reading takes, and the
    /// request's measure in that sizing's: the sizing reads `unit` where the
    /// request is at a factor of one.
    fn reading<R>(
        &self,
        vary: Vary,
        then: impl FnOnce(&BlendSizing<'_>, f64) -> OgeomResult<R>,
    ) -> OgeomResult<R> {
        let lowered = self.lower()?;
        let what = lowered.what();
        let rolling = self.round_corner == FilletCorner::Rolling;
        let defaults = self.overflow == Overflow::default();
        if !self.setbacks.is_empty() {
            return Err(refused(format!(
                "a reading of {what} with setback corners at vertices is not read"
            )));
        }
        let scale_only = || -> OgeomResult<()> {
            if vary == Vary::Scale {
                Ok(())
            } else {
                Err(refused(format!(
                    "a law reading of {what}, a request with no law, is not read"
                )))
            }
        };
        let no_rolling = || -> OgeomResult<()> {
            if rolling {
                Err(refused(format!(
                    "a reading of {what} with rolling corners is not read"
                )))
            } else {
                Ok(())
            }
        };
        match lowered {
            Lowered::Rounds { edges, radius } => {
                scale_only()?;
                let sizing = match (rolling, defaults) {
                    (_, true) => BlendSizing::Rounds {
                        edges: &edges,
                        corner: self.round_corner,
                    },
                    (false, false) => BlendSizing::RoundsOverflow {
                        edges: &edges,
                        overflow: self.overflow,
                    },
                    (true, false) => BlendSizing::RoundsCornerOverflow {
                        edges: &edges,
                        corner: self.round_corner,
                        overflow: self.overflow,
                    },
                };
                then(&sizing, usable(radius, "rounds of radius")?)
            }
            Lowered::Chamfers(specs) => {
                scale_only()?;
                let unit = bevel_measure(&specs)?;
                if let ([(edge, spec)], true) = (specs.as_slice(), defaults) {
                    then(
                        &BlendSizing::Chamfer {
                            edge,
                            spec,
                            corner: self.chamfer_corner,
                        },
                        unit,
                    )
                } else {
                    then(
                        &BlendSizing::Chamfers {
                            specs: &specs,
                            corner: self.chamfer_corner,
                            overflow: self.overflow,
                        },
                        unit,
                    )
                }
            }
            Lowered::Distances {
                edge,
                face,
                on_face,
                on_other,
                profile,
            } => {
                scale_only()?;
                no_rolling()?;
                let sizing = if defaults {
                    BlendSizing::RoundDistances {
                        edge: &edge,
                        face,
                        on_face,
                        on_other,
                        profile,
                    }
                } else {
                    BlendSizing::RoundDistancesOverflow {
                        edge: &edge,
                        face,
                        on_face,
                        on_other,
                        profile,
                        overflow: self.overflow,
                    }
                };
                then(&sizing, 1.0)
            }
            Lowered::RoundLaw { chain, law } => {
                no_rolling()?;
                let size = match vary {
                    Vary::Scale => LawSize::Scale,
                    Vary::Law(at) => LawSize::At(at),
                    Vary::OnFace(_) | Vary::OnOther(_) => {
                        return Err(refused(
                            "a distance-law reading of a round law is not read".into(),
                        ));
                    }
                };
                let sizing = if defaults {
                    BlendSizing::RoundLaw { chain, law, size }
                } else {
                    BlendSizing::RoundLawOverflow {
                        chain,
                        law,
                        size,
                        overflow: self.overflow,
                    }
                };
                then(&sizing, 1.0)
            }
            Lowered::DistanceLaws {
                chain,
                face,
                on_face,
                on_other,
                profile,
            } => {
                no_rolling()?;
                if !defaults {
                    return Err(refused(
                        "a reading of distance laws under overflow options is not read".into(),
                    ));
                }
                let [edge] = chain else {
                    return Err(refused(format!(
                        "a reading of distance laws over a chain of {} edges is not read",
                        chain.len()
                    )));
                };
                let size = match vary {
                    Vary::Scale => DistanceLawSize::Scale,
                    Vary::OnFace(at) => DistanceLawSize::OnFace(at),
                    Vary::OnOther(at) => DistanceLawSize::OnOther(at),
                    Vary::Law(_) => {
                        return Err(refused(
                            "a round-law reading of distance laws is not read".into(),
                        ));
                    }
                };
                then(
                    &BlendSizing::DistanceLaws {
                        edge,
                        face,
                        on_face,
                        on_other,
                        profile,
                        size,
                    },
                    1.0,
                )
            }
            Lowered::RoundsWith(_) | Lowered::Mixed(_) => {
                Err(refused(format!("a reading of {what} is not read")))
            }
        }
    }

    /// A request of `specs`' bevels in their order. With no bevels it is
    /// still a request of bevels, refused by their call's own text.
    pub(crate) fn bevels(specs: &[(Shape, Chamfer)]) -> Self {
        if specs.is_empty() {
            return Self::new().chamfer(std::iter::empty(), 1.0);
        }
        specs.iter().fold(Self::new(), |request, (edge, spec)| {
            request.edge(edge.clone(), spec.clone())
        })
    }

    /// The request as the call that builds it. A request with no edges
    /// takes the kind of its first builder call, rounds where there is
    /// none, so the call refuses it by its own text.
    fn lower(&self) -> OgeomResult<Lowered<'_>> {
        let asked: Vec<&Entry> = self
            .entries
            .iter()
            .filter(|e| !e.edges.is_empty())
            .collect();
        if asked.is_empty() {
            return Ok(match self.entries.first() {
                None => Lowered::Rounds {
                    edges: Vec::new(),
                    radius: 1.0,
                },
                Some(entry) => lowered_alone(entry, Vec::new())?,
            });
        }
        let laws = asked.iter().filter(|e| is_law(&e.section)).count();
        if laws > 0 {
            let [entry] = asked.as_slice() else {
                return Err(refused(
                    "a law with other blends in one request is not built".into(),
                ));
            };
            return lowered_alone(entry, Vec::new());
        }
        let pairs: Vec<(Shape, &BlendSection)> = asked
            .iter()
            .flat_map(|e| e.edges.iter().map(|edge| (edge.clone(), &e.section)))
            .collect();
        let radius = |s: &BlendSection| match s {
            BlendSection::Round(r) => Some(*r),
            _ => None,
        };
        if let Some(first) = radius(pairs[0].1)
            && pairs.iter().all(|(_, s)| radius(s).is_some())
        {
            // Bit for bit as well, so that rounds of one radius that is not
            // a number are still rounds of one radius.
            if pairs.iter().all(|(_, s)| {
                radius(s).is_some_and(|r| r == first || r.to_bits() == first.to_bits())
            }) {
                return Ok(Lowered::Rounds {
                    edges: pairs.into_iter().map(|(e, _)| e).collect(),
                    radius: first,
                });
            }
            return Ok(Lowered::RoundsWith(
                pairs
                    .into_iter()
                    .map(|(e, s)| (e, radius(s).unwrap_or(first)))
                    .collect(),
            ));
        }
        if pairs
            .iter()
            .all(|(_, s)| matches!(s, BlendSection::Chamfer(_)))
        {
            return Ok(Lowered::Chamfers(
                pairs
                    .into_iter()
                    .filter_map(|(e, s)| match s {
                        BlendSection::Chamfer(c) => Some((e, c.clone())),
                        _ => None,
                    })
                    .collect(),
            ));
        }
        if let [(edge, BlendSection::RoundDistances { .. })] = pairs.as_slice() {
            return lowered_alone(asked[0], vec![edge.clone()]);
        }
        Ok(Lowered::Mixed(
            pairs
                .into_iter()
                .filter_map(|(e, s)| constant(s).map(|b| (e, b)))
                .collect(),
        ))
    }
}

const fn is_law(section: &BlendSection) -> bool {
    matches!(
        section,
        BlendSection::RoundLaw(_) | BlendSection::DistanceLaws { .. }
    )
}

/// A constant section as [`blend_edges`](crate::blend_edges) takes it.
fn constant(section: &BlendSection) -> Option<Blend> {
    Some(match section {
        BlendSection::Round(r) => Blend::Round(*r),
        BlendSection::Chamfer(c) => Blend::Chamfer(c.clone()),
        BlendSection::RoundDistances {
            face,
            on_face,
            on_other,
            profile,
        } => Blend::RoundDistances {
            face: face.clone(),
            on_face: *on_face,
            on_other: *on_other,
            profile: *profile,
        },
        BlendSection::RoundLaw(_) | BlendSection::DistanceLaws { .. } => return None,
    })
}

/// `entry` alone in its request, its edges `edges` where its section is
/// constant (none for a request of no edges, its one edge for a round of
/// two distances), its chain as given under a law.
fn lowered_alone(entry: &Entry, edges: Vec<Shape>) -> OgeomResult<Lowered<'_>> {
    Ok(match &entry.section {
        BlendSection::Round(radius) => Lowered::Rounds {
            edges,
            radius: *radius,
        },
        BlendSection::Chamfer(c) => {
            Lowered::Chamfers(edges.into_iter().map(|e| (e, c.clone())).collect())
        }
        BlendSection::RoundDistances {
            face,
            on_face,
            on_other,
            profile,
        } => match edges.into_iter().next() {
            Some(edge) => Lowered::Distances {
                edge,
                face,
                on_face: *on_face,
                on_other: *on_other,
                profile: *profile,
            },
            None => {
                return Err(refused(
                    "a round of two distances on no edge rounds nothing".into(),
                ));
            }
        },
        BlendSection::RoundLaw(law) => Lowered::RoundLaw {
            chain: &entry.edges,
            law,
        },
        BlendSection::DistanceLaws {
            face,
            on_face,
            on_other,
            profile,
        } => Lowered::DistanceLaws {
            chain: &entry.edges,
            face,
            on_face,
            on_other,
            profile: *profile,
        },
    })
}

/// `measure` where it is a length a scale reading can multiply.
fn usable(measure: f64, what: &str) -> OgeomResult<f64> {
    if measure.is_finite() && measure > 0.0 {
        Ok(measure)
    } else {
        Err(refused(format!(
            "a scale reading of {what} {measure} reads nothing"
        )))
    }
}

/// The measure every bevel's reading scales from: the distance of a
/// symmetric or angle bevel, one for two distances (the reading scales
/// both). Bevels of different measures share no such unit.
fn bevel_measure(specs: &[(Shape, Chamfer)]) -> OgeomResult<f64> {
    let measure = |c: &Chamfer| match c {
        Chamfer::Symmetric(d) | Chamfer::Angle { distance: d, .. } => *d,
        Chamfer::Distances { .. } => 1.0,
    };
    let Some((_, first)) = specs.first() else {
        return Ok(1.0);
    };
    let unit = measure(first);
    if specs.iter().any(|(_, c)| measure(c) != unit) {
        return Err(refused(
            "a scale reading of bevels of different measures is not read".into(),
        ));
    }
    usable(unit, "bevels of measure")
}
