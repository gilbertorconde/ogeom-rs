//! Where a point sits relative to a shape.
//!
//! Three answers, never two: inside, outside, or *on* the boundary within
//! tolerance. The third is not a hedge. Geometry that meets is the normal case
//! in a kernel (a boolean's whole job is finding it), and a classifier that
//! forces every point to one side has to pick, silently, for exactly the points
//! where the choice matters most.
//!
//! # Accuracy
//!
//! [`classify_in_solid`] and [`classify_on_face`] work from the tessellation,
//! so a point nearer the boundary than the deflection cannot be told from one
//! on it. That is reported as [`Containment::On`] rather than guessed: the
//! band the answer is uncertain within is the deflection, and saying so is the
//! difference between an approximate answer and a wrong one.
//!
//! Tightening the deflection narrows the band. It never removes it: the
//! exact question needs ray/surface intersection, which is
//! [`classify_in_solid_exact`]: rays cast against the faces' true surfaces,
//! where the uncertain band shrinks from the deflection to the tolerance.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Aabb, Direction, Point, Point2, Vector};
use ogeom_mesh::{Deflection, face_boundary, inside_boundary, open_chart_ring, triangulate};
use ogeom_topo::{Model, NodeData, Shape, ShapeType};

use crate::measure::project_on_surface;

/// Where a point sits relative to a shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Containment {
    /// Strictly inside.
    In,
    /// On the boundary, within tolerance of it.
    On,
    /// Strictly outside.
    Out,
}

impl Containment {
    /// Whether the point is inside or on the boundary.
    #[must_use]
    pub const fn is_inside_or_on(self) -> bool {
        matches!(self, Self::In | Self::On)
    }

    /// The classification of the same point against the complement.
    ///
    /// `In` and `Out` swap; `On` is its own opposite, since a boundary is
    /// shared by both sides.
    #[must_use]
    pub const fn inverted(self) -> Self {
        match self {
            Self::In => Self::Out,
            Self::On => Self::On,
            Self::Out => Self::In,
        }
    }
}

/// Where a point sits relative to a face.
///
/// A point off the face's surface is [`Containment::Out`]: a face is a patch of
/// surface, so "inside" can only mean inside its trimming, and a point in space
/// that does not lie on the surface at all is not inside anything.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `face` is not a
/// face, or the deflection settings are unusable;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn classify_on_face(
    model: &Model,
    face: &Shape,
    point: Point,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Containment> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };

    // Into the surface's own frame first: the trimming lives in parameter
    // space, and the face may be placed anywhere.
    let placement = face.transform(model.datums())?;
    let local = placement.inverse()?.apply(point);

    // A grid dense enough to bracket a foot point on a surface that folds:
    // too coarse and Newton starts in the wrong basin and converges on a far
    // side of a cylinder.
    let projection = project_on_surface(surface, local, 32, tol)?;
    let reach = tol.confusion().max(data.tolerance.get());
    if projection.distance > reach {
        return Ok(Containment::Out);
    }

    // A plane's holes are drawn only near the point: a face with hundreds
    // of holes is asked about beside a few of them.
    let rings = Rings::of(model, face, surface, deflection, tol)?;
    let (u, v) = projection.parameters;
    let at = rings.place(surface, Point2::new(u, v), tol);
    // The uncertain band, converted from a distance in space into one in
    // parameter units through the surface's own scale. A fixed parameter
    // tolerance would be metres wide at a sphere's equator and nothing at its
    // pole. The rings are polylines drawn at the caller's deflection, so the
    // band carries that sag too: without it, a point between a tangent chord
    // and its arc (inside the true trim, outside the sampled one) would
    // read as Out when the honest answer at this resolution is On.
    let band = parametric_band(surface, (u, v), reach + deflection.chord, tol);
    if rings.within(model, face, deflection, at, band, tol)? {
        return Ok(Containment::On);
    }
    Ok(if rings.inside(model, face, deflection, at, tol)? {
        Containment::In
    } else {
        Containment::Out
    })
}

/// Where a point sits relative to a closed shell or solid.
///
/// Ray casting against the tessellation: a ray from the point crosses the
/// boundary an odd number of times if and only if it started inside.
///
/// # Errors
///
/// As [`triangulate()`], plus
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the boundary is
/// not closed (an open shell has no inside), and
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if every ray tried hit an
/// edge or a vertex, where the crossing count is ambiguous.
pub fn classify_in_solid(
    model: &Model,
    solid: &Shape,
    point: Point,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Containment> {
    SolidMesh::of(model, solid, deflection, tol)?.holds(point, tol)
}

/// A solid's boundary meshed once, to classify many points against.
///
/// [`classify_in_solid`] triangulates the whole solid for every point it is
/// asked about; a caller stepping a probe along a ray (a fillet feeling for
/// where its ball runs out of material) asks dozens of times about one
/// solid, and pays for dozens of meshes. This keeps the one.
#[derive(Debug, Clone)]
pub struct SolidMesh {
    triangles: Vec<[Point; 3]>,
    bound: ogeom_math::Aabb,
    reach: f64,
}

impl SolidMesh {
    /// Mesh `solid` at `deflection`.
    ///
    /// # Errors
    ///
    /// As [`classify_in_solid`].
    pub fn of(
        model: &Model,
        solid: &Shape,
        deflection: Deflection,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        deflection.validate()?;
        let mesh = triangulate(model, solid, deflection, tol)?;
        if mesh.is_empty() || !mesh.is_closed() {
            ogeom_bail!(
                Construction,
                "the boundary is not closed, so there is no inside to be in"
            );
        }
        let triangles: Vec<[Point; 3]> = mesh
            .triangles
            .iter()
            .map(|t| t.map(|i| mesh.positions[i as usize]))
            .collect();
        let bound = ogeom_math::Aabb::of_points(&mesh.positions);
        Ok(Self {
            triangles,
            bound,
            reach: tol.confusion() + deflection.chord,
        })
    }

    /// Where `point` stands against the meshed boundary.
    ///
    /// # Errors
    ///
    /// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if every ray
    /// tried hit an edge or a vertex, where the crossing count is ambiguous.
    pub fn holds(&self, point: Point, tol: Tolerances) -> OgeomResult<Containment> {
        // Outside the mesh's box by more than the boundary band, nothing is
        // near enough to be on it and no ray can cross it.
        if !self.bound.expanded(self.reach).contains(point) {
            return Ok(Containment::Out);
        }
        // On the boundary beats either side, and is decided in space rather
        // than along a ray: a point sitting on a face is on the boundary from
        // every direction, and no crossing count says so.
        for t in &self.triangles {
            if distance_to_triangle(point, *t) <= self.reach {
                return Ok(Containment::On);
            }
        }
        // A ray that grazes an edge or passes through a vertex is counted
        // once by one triangle and twice by its neighbour, or not at all.
        // Rather than patch the count, notice the near-miss and cast again
        // somewhere else.
        for direction in RAY_DIRECTIONS {
            let ray = Direction::new(Vector::new(direction[0], direction[1], direction[2]), tol)?;
            if let Some(crossings) = count_crossings(&self.triangles, point, ray, tol) {
                return Ok(if crossings % 2 == 1 {
                    Containment::In
                } else {
                    Containment::Out
                });
            }
        }
        ogeom_bail!(
            NotDone,
            "every ray tried met an edge or a vertex, where the crossing count is \
             ambiguous"
        )
    }
}

/// Where a point sits relative to a closed shell or solid, decided against
/// the true surfaces.
///
/// Rays are cast against each face's actual geometry through the
/// curve/surface intersector, so the band where the answer is *On* rather
/// than a side is the tolerance, not the deflection. A
/// point a micron off a sphere's wall classifies as the side it is on;
/// [`classify_in_solid`] at any practical deflection could only say *On*.
///
/// The crossing-parity argument is the same as the tessellated one, and so is
/// the discipline about degenerate hits. A ray that grazes a surface
/// tangentially, meets a face too near its boundary to be sure which side of
/// the trim it crossed, lies *in* a face's surface, or passes through a
/// pole or an apex, is not patched into a count; the ray is abandoned and
/// the next direction tried. Six directions, deterministic, none axis-aligned.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the shape is
/// not a solid or a shell, or its boundary is not closed;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if every ray met a
/// degeneracy, which a handful of deliberately skew directions makes an
/// engineered case rather than an encountered one.
pub fn classify_in_solid_exact(
    model: &Model,
    solid: &Shape,
    point: Point,
    tol: Tolerances,
) -> OgeomResult<Containment> {
    classify_in_solid_exact_banded(model, solid, point, tol.confusion() * 1e4, tol)
}

/// [`classify_in_solid_exact`] with the boundary band under the caller's
/// control.
///
/// The band is the ring polylines' chord tolerance, and with it the width of
/// the region that answers `On`. The default is generous (a boolean wants a
/// piece near a boundary called On and resolved against its partner), but a
/// caller that got On *without* a partner to resolve against needs to ask
/// again at a width where proximity stops impersonating coincidence.
///
/// # Errors
///
/// As [`classify_in_solid_exact`].
pub fn classify_in_solid_exact_banded(
    model: &Model,
    solid: &Shape,
    point: Point,
    ring_chord: f64,
    tol: Tolerances,
) -> OgeomResult<Containment> {
    SolidBoundary::of(model, solid, ring_chord, tol)?.holds(model, point, tol)
}

/// One face of a prepared boundary, with everything about it that does not
/// depend on the point being classified.
#[derive(Debug)]
struct PreparedFace {
    surface: ogeom_geom::SurfaceGeometry,
    /// The placement's inverse, for carrying a point into the surface's frame.
    inverse: ogeom_math::Transform,
    /// The face itself, for drawing its rings when first asked.
    face: Shape,
    /// The trimming rings, polylined at the boundary's stated chord, drawn
    /// the first time a point or a ray comes near the face: a boolean asks
    /// about a few faces of a large solid, and drawing every face's rings
    /// up front cost more than all its questions.
    rings: std::sync::OnceLock<OgeomResult<Rings>>,
    /// Where the face can be, padded past anything its bound could miss: a
    /// point outside is not on it, and a ray missing it does not cross it.
    bound: Aabb,
    /// How far off the surface a point still lies on the face: the face's
    /// own tolerance, never under the confusion distance.
    reach: f64,
}

/// A face's trimming rings at one chord.
#[derive(Debug)]
enum Rings {
    /// Every ring, drawn together.
    All(Vec<Vec<Point2>>),
    /// A plane's outer ring, and its holes each drawn the first time a
    /// point comes near it: a face with hundreds of holes is asked about
    /// near a few.
    Open {
        outer: Vec<Point2>,
        holes: Vec<Hole>,
    },
}

/// A hole of a plane face, its ring drawn when first needed.
#[derive(Debug)]
struct Hole {
    wire: Shape,
    /// A chart box holding every point its ring can have.
    low: Point2,
    high: Point2,
    ring: std::sync::OnceLock<OgeomResult<Option<Vec<Point2>>>>,
}

impl Hole {
    fn ring(
        &self,
        model: &Model,
        face: &Shape,
        deflection: Deflection,
        tol: Tolerances,
    ) -> OgeomResult<Option<&[Point2]>> {
        match self
            .ring
            .get_or_init(|| open_chart_ring(model, face, &self.wire, false, deflection, tol))
        {
            Ok(ring) => Ok(ring.as_deref()),
            Err(e) => Err(e.clone()),
        }
    }

    /// How far `p` stands from the box, nothing inside it.
    fn box_distance(&self, p: Point2) -> f64 {
        let dx = (self.low.x - p.x).max(p.x - self.high.x).max(0.0);
        let dy = (self.low.y - p.y).max(p.y - self.high.y).max(0.0);
        dx.hypot(dy)
    }
}

impl Rings {
    /// A face's rings as [`face_boundary`] draws them: a plane's with
    /// several wires one wire at a time, where each hole's box can be read
    /// from its pcurves and is smaller than the outer ring, so the outer
    /// ring is the largest, as drawing all of them would find.
    fn of(
        model: &Model,
        face: &Shape,
        surface: &ogeom_geom::SurfaceGeometry,
        deflection: Deflection,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let all = || face_boundary(model, face, deflection, tol).map(Rings::All);
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return all();
        };
        // Stored order, the outer wire first: the walk of a reversed face
        // lists its holes first.
        let wires = model.children_of(face)?;
        if !matches!(surface, ogeom_geom::SurfaceGeometry::Plane(_)) || wires.len() < 2 {
            return all();
        }
        let Some(outer) = open_chart_ring(model, face, &wires[0], true, deflection, tol)? else {
            return all();
        };
        let mut area = 0.0;
        for i in 0..outer.len() {
            let (p, q) = (outer[i], outer[(i + 1) % outer.len()]);
            area += p.x * q.y - q.x * p.y;
        }
        let area = area.abs() * 0.5;
        let mut holes = Vec::with_capacity(wires.len() - 1);
        for wire in &wires[1..] {
            let Some((low, high)) = chart_box(model, data.surface, wire, tol) else {
                return all();
            };
            if (high.x - low.x) * (high.y - low.y) >= area {
                return all();
            }
            holes.push(Hole {
                wire: wire.clone(),
                low,
                high,
                ring: std::sync::OnceLock::new(),
            });
        }
        Ok(Rings::Open { outer, holes })
    }

    /// [`place_on_rings`]: on a plane the point stands where it is.
    fn place(&self, surface: &ogeom_geom::SurfaceGeometry, at: Point2, tol: Tolerances) -> Point2 {
        match self {
            Rings::All(rings) => place_on_rings(surface, rings, at, tol),
            Rings::Open { .. } => at,
        }
    }

    /// Whether `at` lies within `band` of a ring.
    fn within(
        &self,
        model: &Model,
        face: &Shape,
        deflection: Deflection,
        at: Point2,
        band: ChartBand,
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        match self {
            Rings::All(rings) => Ok(band.meets_rings(rings, at)),
            Rings::Open { outer, holes } => {
                if band.meets_ring(outer, at) {
                    return Ok(true);
                }
                for hole in holes {
                    if hole.box_distance(at) <= band.radius()
                        && let Some(ring) = hole.ring(model, face, deflection, tol)?
                        && band.meets_ring(ring, at)
                    {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
        }
    }

    /// Whether `at` lies inside the rings, by the even-odd count.
    fn inside(
        &self,
        model: &Model,
        face: &Shape,
        deflection: Deflection,
        at: Point2,
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        match self {
            Rings::All(rings) => Ok(inside_boundary(rings, at)),
            Rings::Open { outer, holes } => {
                let mut inside = inside_boundary(core::slice::from_ref(outer), at);
                for hole in holes {
                    if hole.box_distance(at) <= 0.0
                        && let Some(ring) = hole.ring(model, face, deflection, tol)?
                        && inside_boundary(&[ring.to_vec()], at)
                    {
                        inside = !inside;
                    }
                }
                Ok(inside)
            }
        }
    }
}

/// A chart box holding every point a wire's pcurves on `surface` reach,
/// for lines and circles; `None` for a wire with another kind of pcurve.
fn chart_box(
    model: &Model,
    surface: ogeom_topo::SurfaceId,
    wire: &Shape,
    tol: Tolerances,
) -> Option<(Point2, Point2)> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::PlanarCurve;
    let mut low = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut high = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut take = |p: Point2, rx: f64, ry: f64| {
        low = Point2::new(low.x.min(p.x - rx), low.y.min(p.y - ry));
        high = Point2::new(high.x.max(p.x + rx), high.y.max(p.y + ry));
    };
    // How far a ring's ends may be moved onto the vertices they belong to:
    // a vertex lies within its tolerance of its edge's end, and the edge
    // within its own of its pcurve.
    let mut loose = tol.confusion();
    for edge in model.ordered_children_of(wire).ok()? {
        let data = model.node(&edge)?.data().as_edge()?;
        loose = loose.max(data.tolerance.get());
        for vertex in model.children_of(&edge).ok()? {
            if let Some(v) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                loose = loose.max(v.tolerance.get());
            }
        }
        let ogeom_topo::EdgeRepr::PCurve { curve, range, .. } =
            data.pcurve_for(surface, edge.location())?
        else {
            return None;
        };
        match model.geometry().pcurve(*curve)? {
            line @ PlanarCurve::Line(_) => {
                take(line.point_at(range.0, tol).ok()?, 0.0, 0.0);
                take(line.point_at(range.1, tol).ok()?, 0.0, 0.0);
            }
            PlanarCurve::Circle(arc) => {
                let circle = arc.circle();
                let (x, y) = (circle.frame().x().vector(), circle.frame().y().vector());
                let r = circle.radius();
                take(
                    circle.centre(),
                    r * (x.x.abs() + y.x.abs()),
                    r * (x.y.abs() + y.y.abs()),
                );
            }
            _ => return None,
        }
    }
    if !(low.x.is_finite() && low.y.is_finite() && high.x.is_finite() && high.y.is_finite()) {
        return None;
    }
    // Those moves twice over, and rounding in the points drawn and in a
    // line's points just past its ends.
    let slack = loose.mul_add(
        4.0,
        1e-9 * (1.0
            + low
                .x
                .abs()
                .max(low.y.abs())
                .max(high.x.abs())
                .max(high.y.abs())),
    );
    Some((
        Point2::new(low.x - slack, low.y - slack),
        Point2::new(high.x + slack, high.y + slack),
    ))
}

impl PreparedFace {
    /// The face's trimming rings at `deflection`, drawn once.
    fn rings(&self, model: &Model, deflection: Deflection, tol: Tolerances) -> OgeomResult<&Rings> {
        match self
            .rings
            .get_or_init(|| Rings::of(model, &self.face, &self.surface, deflection, tol))
        {
            Ok(rings) => Ok(rings),
            Err(e) => Err(e.clone()),
        }
    }

    /// Where a point sits against this face: [`classify_on_face`] on what
    /// was prepared, with no surface lookup, placement or ring walk per
    /// question.
    fn holds(
        &self,
        model: &Model,
        point: Point,
        deflection: Deflection,
        tol: Tolerances,
    ) -> OgeomResult<Containment> {
        let local = self.inverse.apply(point);
        let projection = project_on_surface(&self.surface, local, 32, tol)?;
        if projection.distance > self.reach {
            return Ok(Containment::Out);
        }
        // Against the rings drawn so far.
        let rings = self.rings(model, deflection, tol)?;
        let (u, v) = projection.parameters;
        let at = rings.place(&self.surface, Point2::new(u, v), tol);
        let band = parametric_band(&self.surface, (u, v), self.reach + deflection.chord, tol);
        if rings.within(model, &self.face, deflection, at, band, tol)? {
            return Ok(Containment::On);
        }
        Ok(if rings.inside(model, &self.face, deflection, at, tol)? {
            Containment::In
        } else {
            Containment::Out
        })
    }
}

/// A solid's boundary, prepared once and asked about many points.
///
/// Classifying a point casts rays and counts crossings, which is cheap. What
/// is not cheap is what the rays are cast *against*: every face's trimming
/// rings, polylined, plus its placement's inverse. That work depends on the
/// solid and the chord, never on the point, and a boolean asks once per face
/// piece. Preparing costs hundreds of times the ray casting it serves, so it
/// is done once.
#[derive(Debug)]
pub struct SolidBoundary {
    faces: Vec<PreparedFace>,
    bound: ogeom_math::Aabb,
    centre: Point,
    diagonal: f64,
    ring_chord: f64,
}

impl SolidBoundary {
    /// Prepare a solid's boundary for classification at a given ring chord.
    ///
    /// # Errors
    ///
    /// As [`classify_in_solid_exact`].
    pub fn of(model: &Model, solid: &Shape, ring_chord: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::prepare(model, solid, ring_chord, tol)
    }

    fn prepare(
        model: &Model,
        solid: &Shape,
        ring_chord: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let kind = model.kind_of(solid)?;
        // A compound of solids bounds their union: each lump's shells are
        // boundary, and a ray counts the crossings of all of them.
        let lumps = kind == ShapeType::Compound
            && model
                .children_of(solid)?
                .iter()
                .all(|part| model.kind_of(part).is_ok_and(|k| k == ShapeType::Solid));
        if !matches!(kind, ShapeType::Solid | ShapeType::Shell) && !lumps {
            ogeom_bail!(Construction, "expected a solid or a shell, got {kind:?}");
        }
        let shells = if kind == ShapeType::Shell {
            vec![solid.clone()]
        } else {
            ogeom_topo::explore_unique(model, solid, ShapeType::Shell)?
        };
        if shells.is_empty() {
            ogeom_bail!(Construction, "the shape has no shell, so no boundary");
        }
        for shell in &shells {
            if !crate::build::is_shell_closed(model, shell)? {
                ogeom_bail!(
                    Construction,
                    "the boundary is not closed, so there is no inside to be in"
                );
            }
        }

        let faces = ogeom_topo::explore_unique(model, solid, ShapeType::Face)?;
        // One prepared face per face, in face order, computed in parallel:
        // each preparation reads the model and writes nothing, and walking a
        // face's trimming rings is the whole cost of building a boundary.
        let prepared = ogeom_core::parallel::map_ordered(&faces, |_, face| {
            ogeom_core::progress::checkpoint()?;
            let Some(node) = model.node(face) else {
                ogeom_bail!(Dangling, "face is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "face node holds no face data");
            };
            let Some(surface) = model.geometry().surface(data.surface) else {
                ogeom_bail!(Dangling, "face refers to a surface not in this model");
            };
            let inverse = face.transform(model.datums())?.inverse()?;
            let own = crate::measure::shape_bounds(model, face, tol)?;
            let bound = own.expanded(
                ring_chord + data.tolerance.get() + tol.confusion() * 1e2 + own.diagonal() * 0.02,
            );
            Ok((
                PreparedFace {
                    surface: surface.clone(),
                    inverse,
                    face: face.clone(),
                    rings: std::sync::OnceLock::new(),
                    bound,
                    reach: tol.confusion().max(data.tolerance.get()),
                },
                own,
            ))
        })
        .into_iter()
        .collect::<OgeomResult<Vec<_>>>()?;
        // Anything outside the shape's bound is outside the shape, and the bound
        // also sets how long a ray must be to have left everything behind. A
        // solid holds shells and a shell faces, so the faces' bounds together
        // are the shape's.
        let bound = prepared
            .iter()
            .fold(Aabb::EMPTY, |bound, (_, own)| bound.union(own));
        let (Some(centre), diagonal) = (bound.centre(), bound.diagonal()) else {
            ogeom_bail!(Construction, "the boundary bounds nothing");
        };
        let prepared = prepared.into_iter().map(|(face, _)| face).collect();
        Ok(Self {
            faces: prepared,
            bound,
            centre,
            diagonal,
            ring_chord,
        })
    }

    /// Where a point stands against this boundary.
    ///
    /// # Errors
    ///
    /// As [`classify_in_solid_exact`].
    pub fn holds(&self, model: &Model, point: Point, tol: Tolerances) -> OgeomResult<Containment> {
        let ring_chord = self.ring_chord;
        let ring_deflection = Deflection {
            chord: ring_chord,
            angular: 0.05,
            ..Deflection::default()
        };
        let reach = tol.confusion();
        if !self.bound.expanded(reach).contains(point) {
            return Ok(Containment::Out);
        }
        let length = point.distance(self.centre) + self.diagonal + 1.0;

        // On the boundary beats either side, and each face answers exactly:
        // projection distance against the true surface, trimming in parameter
        // space.
        for prepared in &self.faces {
            if !prepared.bound.contains(point) {
                continue;
            }
            if prepared.holds(model, point, ring_deflection, tol)? != Containment::Out {
                return Ok(Containment::On);
            }
        }
        'directions: for direction in RAY_DIRECTIONS {
            let along = Vector::new(direction[0], direction[1], direction[2]);
            let far = point + along * length;
            let mut crossings = 0_usize;

            for prepared in &self.faces {
                let PreparedFace {
                    surface,
                    inverse,
                    bound,
                    ..
                } = prepared;
                if !segment_meets(bound, point, far) {
                    continue;
                }
                // Into the face's frame, as two points rather than a direction, so
                // a placement that scales still carries the ray faithfully.
                let from = inverse.apply(point);
                let to = inverse.apply(far);
                let ray: ogeom_geom::Curve = ogeom_geom::LineCurve::segment(from, to, tol)?.into();
                let found = ogeom_intersect::intersect_curve_surface(
                    &ray,
                    surface,
                    ogeom_intersect::CurveSurfaceOptions::default(),
                    tol,
                )?;
                if !found.lying.is_empty() {
                    // The ray runs in this face's surface: it crosses nothing and
                    // touches everything, which no parity expresses.
                    continue 'directions;
                }
                for hit in &found.crossings {
                    if hit.on_curve <= tol.confusion() {
                        // At the very start: the probe lies in this face's
                        // *surface*. Whether that matters depends on the trim:
                        // the boundary test above already said the point is off
                        // every face, so a start-crossing far from this face's
                        // rings is the unbounded surface talking, not the face,
                        // and it neither counts nor poisons the ray.
                        let (u, v) = hit.on_surface;
                        let rings = prepared.rings(model, ring_deflection, tol)?;
                        let at = rings.place(surface, Point2::new(u, v), tol);
                        let band = parametric_band(surface, (u, v), reach + ring_chord, tol);
                        let face = &prepared.face;
                        if rings.within(model, face, ring_deflection, at, band, tol)?
                            || rings.inside(model, face, ring_deflection, at, tol)?
                        {
                            continue 'directions;
                        }
                        continue;
                    }
                    let (u, v) = hit.on_surface;
                    use ogeom_geom::Surface as _;
                    let Ok((du, dv)) = surface.d1_at(u, v, tol) else {
                        continue 'directions;
                    };
                    let normal = du.cross(dv);
                    if normal.magnitude() <= tol.confusion() {
                        // A pole or an apex: no normal, no transversality.
                        continue 'directions;
                    }
                    let ray_direction = (to - from) / (to - from).magnitude();
                    if normal.dot(ray_direction).abs() <= GRAZING * normal.magnitude() {
                        // Tangential. A grazing contact counts once where parity
                        // needs zero or two; abandon the ray rather than guess.
                        continue 'directions;
                    }
                    let rings = prepared.rings(model, ring_deflection, tol)?;
                    let at = rings.place(surface, Point2::new(u, v), tol);
                    let band = parametric_band(surface, (u, v), reach + ring_chord, tol);
                    let face = &prepared.face;
                    if rings.within(model, face, ring_deflection, at, band, tol)? {
                        // Too near the face's boundary to know which side of the
                        // trim it crossed, and a shared edge would be counted by
                        // both faces or neither.
                        continue 'directions;
                    }
                    if rings.inside(model, face, ring_deflection, at, tol)? {
                        crossings += 1;
                    }
                }
            }
            return Ok(if crossings % 2 == 1 {
                Containment::In
            } else {
                Containment::Out
            });
        }
        ogeom_bail!(
            NotDone,
            "every ray tried met a tangency, a boundary, or a degenerate point, \
             where the crossing count is ambiguous"
        )
    }
}

/// The sine of the shallowest crossing angle a counted ray/surface crossing
/// may make.
///
/// Below this the hit is treated as a graze: a tangential contact is one
/// crossing where parity arithmetic needs zero or two, and the seed of the
/// threshold is the same as the marching intersector's `SHALLOWEST`: beneath
/// a microradian, rounding in the evaluated normal can no longer tell a
/// crossing from a touch.
const GRAZING: f64 = 1e-6;

/// Directions to cast rays along, tried in order.
///
/// Deterministic, not random: a classifier that gives different answers on
/// different runs is worse than one that fails, because the failure can be
/// handled and the inconsistency cannot. They are deliberately not axis-aligned
/// and share no common plane, so a mesh built on a regular grid (where an
/// axis-aligned ray runs along a whole row of edges) does not defeat all of
/// them at once.
const RAY_DIRECTIONS: [[f64; 3]; 6] = [
    [0.577_35, 0.577_35, 0.577_35],
    [-0.301_5, 0.904_5, 0.301_5],
    [0.727_6, -0.485_1, 0.485_1],
    [0.259_5, 0.259_5, -0.930_0],
    [-0.816_5, -0.408_2, 0.408_2],
    [0.132_5, -0.662_3, -0.737_5],
];

/// Count how many triangles a ray from `from` crosses, or `None` if any hit was
/// too close to an edge or vertex to count reliably.
fn count_crossings(
    triangles: &[[Point; 3]],
    from: Point,
    along: Direction,
    tol: Tolerances,
) -> Option<usize> {
    let mut crossings = 0;
    for t in triangles {
        match ray_hits_triangle(from, along, *t, tol) {
            Hit::Crosses => crossings += 1,
            Hit::Misses => {}
            Hit::Ambiguous => return None,
        }
    }
    Some(crossings)
}

/// What a ray did to a triangle.
enum Hit {
    /// Passed through its interior, ahead of the start.
    Crosses,
    /// Did not meet it.
    Misses,
    /// Met an edge, a vertex, or the plane edge-on, where counting it once is
    /// as defensible as counting it twice or not at all.
    Ambiguous,
}

/// Möller-Trumbore, with the degenerate cases separated out rather than
/// rounded away.
fn ray_hits_triangle(from: Point, along: Direction, t: [Point; 3], tol: Tolerances) -> Hit {
    let direction = along.vector();
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = direction.cross(e2);
    let determinant = e1.dot(h);

    // Scale the comparison by the triangle: a determinant is a volume, so a
    // fixed threshold rejects small triangles and accepts edge-on hits on
    // large ones.
    let scale = e1.magnitude() * e2.magnitude();
    let flat = tol.confusion() * scale;
    if determinant.abs() <= flat {
        // Edge-on. It cannot cross cleanly, but it may lie in the plane and
        // touch, which no crossing count expresses.
        let normal = e1.cross(e2);
        let reach = tol.confusion() * scale;
        return if normal.dot(from - t[0]).abs() <= reach {
            Hit::Ambiguous
        } else {
            Hit::Misses
        };
    }

    let inverse = 1.0 / determinant;
    let s = from - t[0];
    let u = inverse * s.dot(h);
    let q = s.cross(e1);
    let v = inverse * direction.dot(q);
    let w = 1.0 - u - v;

    // Barycentric coordinates near zero mean the ray passed along an edge, and
    // near one that it went through a vertex.
    let edge = tol.confusion();
    if [u, v, w].iter().any(|c| c.abs() <= edge) {
        // Only ambiguous if the ray would otherwise have hit: a grazing miss
        // well outside the triangle is a miss.
        return if u >= -edge && v >= -edge && w >= -edge {
            Hit::Ambiguous
        } else {
            Hit::Misses
        };
    }
    if u < 0.0 || v < 0.0 || w < 0.0 {
        return Hit::Misses;
    }

    let distance = inverse * e2.dot(q);
    if distance <= tol.confusion() {
        // Behind the start, or right at it, and right at it was already ruled
        // out by the on-boundary test before any ray was cast.
        return Hit::Misses;
    }
    Hit::Crosses
}

/// The distance from a point to a triangle.
fn distance_to_triangle(p: Point, t: [Point; 3]) -> f64 {
    // Clamp the projection onto the triangle's plane into the triangle, by
    // checking the three edge regions and the interior. Solving the 2×2 normal
    // equations directly and clamping is shorter than a region case analysis
    // and gives the same closest point.
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let d = t[0] - p;
    let (a, b, c) = (e1.dot(e1), e1.dot(e2), e2.dot(e2));
    let (dd, e) = (e1.dot(d), e2.dot(d));
    let determinant = b.mul_add(-b, a * c);

    if determinant.abs() <= f64::MIN_POSITIVE {
        // A degenerate triangle is a segment; its edges still answer.
        return edge_distance(p, t);
    }
    let mut s = b.mul_add(e, -(c * dd)) / determinant;
    let mut u = b.mul_add(dd, -(a * e)) / determinant;

    if s >= 0.0 && u >= 0.0 && s + u <= 1.0 {
        let closest = t[0] + e1 * s + e2 * u;
        return p.distance(closest);
    }
    // Outside: the closest point is on an edge.
    s = s.clamp(0.0, 1.0);
    u = u.clamp(0.0, 1.0);
    let _ = (s, u);
    edge_distance(p, t)
}

/// The distance from a point to the nearest of a triangle's three edges.
fn edge_distance(p: Point, t: [Point; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..3 {
        best = best.min(segment_distance(p, t[i], t[(i + 1) % 3]));
    }
    best
}

/// The distance from a point to a segment.
fn segment_distance(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let squared = d.dot(d);
    if squared <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(d) / squared).clamp(0.0, 1.0);
    p.distance(a + d * t)
}

/// The distance from a 2D point to a 2D segment.
fn segment_distance_2d(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let squared = d.dot(d);
    if squared <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(d) / squared).clamp(0.0, 1.0);
    p.distance(a + d * t)
}

/// Fold a chart point toward the rings' own window, one period at a time.
///
/// Projection and intersection answer parameters in a surface's principal
/// range, but a face's trim may live in any window of a periodic chart: a
/// band anchored where its rings happened to start. The trim tests compare
/// against the rings, so the point folds to them, not the other way round.
pub(crate) fn fold_toward_rings(
    surface: &ogeom_geom::SurfaceGeometry,
    rings: &[Vec<Point2>],
    mut at: Point2,
) -> Point2 {
    use ogeom_geom::Surface as _;
    // Whatever repeats: the canonical surfaces' angles, a whole turn of a
    // revolution, a periodic patch.
    let ((ua, ub), (va, vb)) = surface.domain();
    let u_period = (surface.is_periodic_u() && ub > ua).then_some(ub - ua);
    let v_period = (surface.is_periodic_v() && vb > va).then_some(vb - va);
    let fold = |x: f64, lo: f64, hi: f64, period: f64| -> f64 {
        let mut x = x;
        while x < lo && x + period <= hi + period {
            x += period;
            if x >= lo {
                break;
            }
        }
        while x > hi && x - period >= lo - period {
            x -= period;
            if x <= hi {
                break;
            }
        }
        x
    };
    if let Some(period) = u_period {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for ring in rings {
            for q in ring {
                lo = lo.min(q.x);
                hi = hi.max(q.x);
            }
        }
        if lo.is_finite() {
            at.x = fold(at.x, lo, hi, period);
        }
    }
    if let Some(period) = v_period {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for ring in rings {
            for q in ring {
                lo = lo.min(q.y);
                hi = hi.max(q.y);
            }
        }
        if lo.is_finite() {
            at.y = fold(at.y, lo, hi, period);
        }
    }
    // The window alone cannot say which side of a seam a point is on when
    // the seam runs diagonally round (a band opened along the widest gap
    // its mesh left): both a point and its copy a period over lie within
    // the rings' extent, only one inside them. The one inside is the face's.
    if !inside_boundary(rings, at) {
        let shifts = [
            u_period.map(|p| (p, 0.0)),
            u_period.map(|p| (-p, 0.0)),
            v_period.map(|p| (0.0, p)),
            v_period.map(|p| (0.0, -p)),
        ];
        if let Some(inside) = shifts
            .into_iter()
            .flatten()
            .map(|(du, dv)| Point2::new(at.x + du, at.y + dv))
            .find(|q| inside_boundary(rings, *q))
        {
            return inside;
        }
    }
    at
}

/// [`fold_toward_rings`], and where the point folds outside them on a
/// revolution, its other place in the chart.
///
/// A profile crossing the axis (a whole circle revolved) covers its surface
/// twice: the point at `(u, v)` is also at `u` a half turn over, on the
/// profile's far side. A projection or a crossing answers either, and a face
/// using one half of the profile holds only one of them.
fn place_on_rings(
    surface: &ogeom_geom::SurfaceGeometry,
    rings: &[Vec<Point2>],
    at: Point2,
    tol: Tolerances,
) -> Point2 {
    let folded = fold_toward_rings(surface, rings, at);
    if inside_boundary(rings, folded) {
        return folded;
    }
    let ogeom_geom::SurfaceGeometry::Revolution(revolution) = surface else {
        return folded;
    };
    let other = (|| -> Option<Point2> {
        use ogeom_geom::Surface as _;
        let point = surface.point_at(at.x, at.y, tol).ok()?;
        let ((u0, u1), _) = surface.domain();
        let half = core::f64::consts::PI;
        let u = [at.x + half, at.x - half]
            .into_iter()
            .find(|u| *u >= u0 - tol.angular() && *u <= u1 + tol.angular())?;
        let back = ogeom_math::Transform::rotation(revolution.axis(), -u).apply(point);
        let found = crate::measure::project_on_curve(revolution.curve(), back, 64, tol).ok()?;
        (found.distance <= tol.confusion()).then(|| Point2::new(u, found.parameter))
    })();
    other
        .map(|q| fold_toward_rings(surface, rings, q))
        .filter(|q| inside_boundary(rings, *q))
        .unwrap_or(folded)
}

/// A distance of `reach` in space about a chart point, read in the chart.
///
/// Each chart direction keeps its own scale, the length of the surface's
/// tangent along it: near a sphere's pole a step in `u` covers a sliver of
/// the space a step in `v` does, and one scale for both would either miss
/// rings along `v` or reach far past them across it. A chart step
/// `(du, dv)` is within the band when `(du * su, dv * sv)` is within
/// `reach`. Where both tangents vanish (a cone's apex) no chart distance
/// corresponds to a spatial one, and the band covers the whole
/// neighbourhood rather than closing to nothing.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ChartBand {
    reach: f64,
    su: f64,
    sv: f64,
}

impl ChartBand {
    /// The widest the band reaches along any chart direction, in chart
    /// units.
    pub(crate) fn radius(self) -> f64 {
        let scale = self.su.min(self.sv);
        if scale <= 0.0 {
            f64::INFINITY
        } else {
            self.reach / scale
        }
    }

    /// Whether a ring passes within the band about `p`.
    pub(crate) fn meets_ring(self, ring: &[Point2], p: Point2) -> bool {
        if self.su <= 0.0 && self.sv <= 0.0 {
            return true;
        }
        let scaled = |q: Point2| Point2::new((q.x - p.x) * self.su, (q.y - p.y) * self.sv);
        let origin = Point2::new(0.0, 0.0);
        (0..ring.len()).any(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            segment_distance_2d(origin, scaled(a), scaled(b)) <= self.reach
        })
    }

    /// Whether any of the rings passes within the band about `p`.
    pub(crate) fn meets_rings(self, rings: &[Vec<Point2>], p: Point2) -> bool {
        rings.iter().any(|ring| self.meets_ring(ring, p))
    }
}

/// The band a distance of `reach` in space makes about `(u, v)` in the
/// chart, scaled by the surface's tangents there.
pub(crate) fn parametric_band(
    surface: &ogeom_geom::SurfaceGeometry,
    at: (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> ChartBand {
    use ogeom_geom::Surface;
    let Ok((du, dv)) = surface.d1_at(at.0, at.1, tol) else {
        return ChartBand {
            reach,
            su: 1.0,
            sv: 1.0,
        };
    };
    // A tangent shorter than the confusion distance is a collapsed one.
    let scale = |t: ogeom_math::Vector| {
        let m = t.magnitude();
        if m <= tol.confusion() { 0.0 } else { m }
    };
    ChartBand {
        reach,
        su: scale(du),
        sv: scale(dv),
    }
}

/// Whether the segment from `a` to `b` passes through the box, by slabs.
fn segment_meets(bound: &Aabb, a: Point, b: Point) -> bool {
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        return false;
    };
    let (mut enter, mut leave) = (0.0_f64, 1.0_f64);
    for (from, to, lo, hi) in [
        (a.x, b.x, low.x, high.x),
        (a.y, b.y, low.y, high.y),
        (a.z, b.z, low.z, high.z),
    ] {
        let d = to - from;
        if d.abs() <= f64::EPSILON * (from.abs() + to.abs() + 1.0) {
            if from < lo || from > hi {
                return false;
            }
            continue;
        }
        let (t0, t1) = ((lo - from) / d, (hi - from) / d);
        let (t0, t1) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        enter = enter.max(t0);
        leave = leave.min(t1);
        if enter > leave {
            return false;
        }
    }
    true
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::make_box;
    use ogeom_math::Frame;
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_point_inside_a_box_is_inside_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(1.0, 1.0, 1.0),
            Point::new(0.1, 0.1, 0.1),
            Point::new(1.9, 1.9, 1.9),
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::In,
                "{p:?} should be inside"
            );
        }
    }

    /// A box with a spherical void bounded by a sphere's face turned inside
    /// out: a point in the void is out of the material, exactly and by mesh.
    #[test]
    fn a_point_in_a_spherical_void_is_outside_the_material() {
        let mut model = Model::new();
        let corner = Frame::new(
            Point::new(-20.0, -20.0, -20.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let block = make_box(&mut model, corner, (40.0, 40.0, 40.0), T).unwrap();
        let ball = crate::make_sphere(&mut model, Frame::WORLD, 10.0, T).unwrap();
        let outer = explore_unique(&model, &block.shape, ShapeType::Shell).unwrap()[0].clone();
        let skin = explore_unique(&model, &ball.shape, ShapeType::Face).unwrap()[0].clone();
        let void = model.add_shell(&[skin.reversed()]).unwrap();
        let holed = crate::make_solid(&mut model, &[outer, void]).unwrap().shape;
        for (p, want) in [
            (Point::new(0.0, 0.0, 0.0), Containment::Out),
            (Point::new(3.0, 3.0, 3.0), Containment::Out),
            (Point::new(-3.0, -3.0, -3.0), Containment::Out),
            (Point::new(15.0, 0.0, 0.0), Containment::In),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &holed, p, T).unwrap(),
                want,
                "{p:?}"
            );
            assert_eq!(
                classify_in_solid(&model, &holed, p, fine(), T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    /// A sphere's cap whose rim runs half a millimetre from the pole: near
    /// the pole a step round the chart covers a twentieth of the space a
    /// step toward it does, and the boundary band is read in space along
    /// each, so points a tenth of a millimetre either side of the rim are in
    /// and out, and only a point within the band is on it.
    #[test]
    fn points_beside_a_rim_near_a_sphere_s_pole_are_in_or_out() {
        use core::f64::consts::TAU;
        let mut model = Model::new();
        let radius = 10.0;
        let sphere = ogeom_math::Sphere::new(Frame::WORLD, radius, T).unwrap();
        let surface: ogeom_geom::SurfaceGeometry = ogeom_geom::SphereSurface::new(sphere).into();
        // The rim: the parallel at a polar distance of 0.5.
        let height = |off: f64| (radius * radius - off * off).sqrt();
        let rim_frame = Frame::new(
            Point::new(0.0, 0.0, height(0.5)),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let circle = ogeom_math::Circle::new(rim_frame, 0.5, T).unwrap();
        let start = crate::make_vertex(&mut model, Point::new(0.5, 0.0, height(0.5))).shape;
        let rim = crate::make_edge_between(
            &mut model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (0.0, TAU),
            &start,
            &start,
            T,
        )
        .unwrap()
        .shape;
        let top = crate::make_vertex(&mut model, Point::new(0.0, 0.0, radius)).shape;
        let mut data = ogeom_topo::EdgeData::new();
        data.degenerate = true;
        let pole = model.add_edge(data, &[top.clone(), top]).unwrap();
        let cap = crate::make_revolution_band(&mut model, &surface, &rim, &pole, T).unwrap();
        // The side probe's own resolution: a chord of a hundredth.
        let deflection = Deflection {
            chord: 1e-2,
            ..Deflection::default()
        };
        // Across the rim, on the side away from the seam.
        let at = |off: f64| Point::new(-off, 0.0, height(off));
        for (off, want) in [
            (0.4, Containment::In),
            (0.45, Containment::In),
            (0.5, Containment::On),
            (0.505, Containment::On),
            (0.55, Containment::Out),
            (0.6, Containment::Out),
        ] {
            assert_eq!(
                classify_on_face(&model, &cap, at(off), deflection, T).unwrap(),
                want,
                "{off} from the pole"
            );
        }
    }

    #[test]
    fn a_point_outside_a_box_is_outside_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(3.0, 1.0, 1.0),
            Point::new(-1.0, 1.0, 1.0),
            Point::new(1.0, 1.0, -0.5),
            Point::new(-5.0, -5.0, -5.0),
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::Out,
                "{p:?} should be outside"
            );
        }
    }

    #[test]
    fn a_point_on_a_boxs_face_is_on_it_rather_than_forced_to_a_side() {
        // The case a two-valued classifier has to guess at, and the case a
        // boolean spends all its time in.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for p in [
            Point::new(1.0, 1.0, 0.0), // face centre
            Point::new(0.0, 1.0, 1.0), // another face
            Point::new(2.0, 2.0, 1.0), // an edge
            Point::ORIGIN,             // a vertex
            Point::new(2.0, 2.0, 2.0), // the far vertex
        ] {
            assert_eq!(
                classify_in_solid(&model, &built.shape, p, fine(), T).unwrap(),
                Containment::On,
                "{p:?} should be on the boundary"
            );
        }
    }

    #[test]
    fn a_ray_along_a_grid_of_edges_does_not_defeat_the_classifier() {
        // A box tessellated into two triangles per face has a diagonal across
        // every face and an edge along every side. An axis-aligned ray from the
        // centre runs straight into a face centre; a diagonal one can run along
        // a triangle edge. The retry is what makes either survivable.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        // The centre of a cube: every axis-aligned ray hits a face centre, and
        // the main diagonal goes through a vertex.
        assert_eq!(
            classify_in_solid(&model, &built.shape, Point::new(1.0, 1.0, 1.0), fine(), T).unwrap(),
            Containment::In
        );
    }

    #[test]
    fn an_open_shell_has_no_inside() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();
        assert!(classify_in_solid(&model, &face, Point::ORIGIN, fine(), T).is_err());
    }

    #[test]
    fn a_point_on_a_face_is_inside_its_trimming_or_not() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        // The face at z = 0, spanning x in [0, 2] and y in [0, 3]. Found by
        // the role make_box gave it, which is what provenance is for.
        let bottom = explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                model
                    .provenance_of(f)
                    .and_then(ogeom_core::Provenance::role)
                    == Some(crate::primitive::roles::FACE_MIN_Z)
            })
            .expect("the box has a face at z = 0");

        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(1.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(5.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::Out,
            "on the surface's plane but outside the trimming"
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(1.0, 1.5, 1.0), fine(), T).unwrap(),
            Containment::Out,
            "off the surface entirely"
        );
        assert_eq!(
            classify_on_face(&model, &bottom, Point::new(0.0, 1.5, 0.0), fine(), T).unwrap(),
            Containment::On,
            "on the trimming boundary"
        );
    }

    #[test]
    fn the_answers_invert_the_way_a_complement_does() {
        assert_eq!(Containment::In.inverted(), Containment::Out);
        assert_eq!(Containment::Out.inverted(), Containment::In);
        // A boundary belongs to both sides, so complementing leaves it alone.
        assert_eq!(Containment::On.inverted(), Containment::On);

        assert!(Containment::In.is_inside_or_on());
        assert!(Containment::On.is_inside_or_on());
        assert!(!Containment::Out.is_inside_or_on());
    }

    #[test]
    fn a_translated_box_classifies_the_same_way_translated_points() {
        let offset = Vector::new(10.0, -20.0, 30.0);
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN + offset, Direction::Z, Direction::X, T).unwrap();
        let built = make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();

        assert_eq!(
            classify_in_solid(
                &model,
                &built.shape,
                Point::new(1.0, 1.0, 1.0) + offset,
                fine(),
                T
            )
            .unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid(&model, &built.shape, Point::new(1.0, 1.0, 1.0), fine(), T).unwrap(),
            Containment::Out,
            "the untranslated point is nowhere near the translated box"
        );
    }

    #[test]
    fn the_exact_classifier_agrees_with_the_tessellated_one_on_a_box() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

        for (p, want) in [
            (Point::new(1.0, 1.0, 1.0), Containment::In),
            (Point::new(0.1, 0.1, 0.1), Containment::In),
            (Point::new(3.0, 1.0, 1.0), Containment::Out),
            (Point::new(1.0, 1.0, -0.5), Containment::Out),
            (Point::new(-50.0, -50.0, -50.0), Containment::Out),
            (Point::new(1.0, 1.0, 0.0), Containment::On),
            (Point::new(2.0, 2.0, 1.0), Containment::On),
            (Point::ORIGIN, Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_resolves_what_the_deflection_band_cannot() {
        // The reason this function exists. A point a micron off a sphere's
        // wall is far inside any practical deflection band; the tessellated
        // classifier must say On, because against a mesh it genuinely cannot
        // tell. Against the true sphere the side is knowable, and known.
        let mut model = Model::new();
        let built = crate::make_sphere(&mut model, Frame::WORLD, 2.0, T).unwrap();

        let barely_in = Point::new(0.0, 0.0, 2.0 - 1e-5);
        let barely_out = Point::new(0.0, 0.0, 2.0 + 1e-5);

        assert_eq!(
            classify_in_solid(&model, &built.shape, barely_in, fine(), T).unwrap(),
            Containment::On,
            "the mesh cannot tell a micron from the wall"
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, barely_in, T).unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, barely_out, T).unwrap(),
            Containment::Out
        );
        // Exactly on the wall: On, decided by projection, not by a ray.
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(0.0, 0.0, 2.0), T).unwrap(),
            Containment::On
        );
    }

    #[test]
    fn the_exact_classifier_handles_a_cylinder_wall_and_caps() {
        let mut model = Model::new();
        let built = crate::make_cylinder(&mut model, Frame::WORLD, 1.5, 4.0, T).unwrap();

        for (p, want) in [
            (Point::new(0.0, 0.0, 2.0), Containment::In),
            (Point::new(1.5 - 1e-5, 0.0, 2.0), Containment::In),
            (Point::new(1.5 + 1e-5, 0.0, 2.0), Containment::Out),
            (Point::new(0.3, 0.4, 4.0 - 1e-5), Containment::In),
            (Point::new(0.3, 0.4, 4.0 + 1e-5), Containment::Out),
            (Point::new(1.5, 0.0, 2.0), Containment::On),
            (Point::new(0.3, 0.4, 0.0), Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_walks_the_general_path_through_a_torus() {
        // No analytic ray/torus case exists, so every crossing here came from
        // the seeded Newton path, and a torus also puts the hole in the
        // middle, where a ray to the outside crosses the tube wall twice.
        let mut model = Model::new();
        let built = crate::make_torus(&mut model, Frame::WORLD, 3.0, 1.0, T).unwrap();

        for (p, want) in [
            (Point::new(3.0, 0.0, 0.0), Containment::In),
            (Point::new(3.0, 0.0, 0.9), Containment::In),
            (Point::ORIGIN, Containment::Out),
            (Point::new(3.0, 0.0, 1.5), Containment::Out),
            (Point::new(5.0, 5.0, 0.0), Containment::Out),
            (Point::new(3.0, 0.0, 1.0), Containment::On),
        ] {
            assert_eq!(
                classify_in_solid_exact(&model, &built.shape, p, T).unwrap(),
                want,
                "{p:?}"
            );
        }
    }

    #[test]
    fn the_exact_classifier_respects_a_placed_solid() {
        let offset = Vector::new(10.0, -20.0, 30.0);
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN + offset, Direction::Z, Direction::X, T).unwrap();
        let built = make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();

        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(1.0, 1.0, 1.0) + offset, T)
                .unwrap(),
            Containment::In
        );
        assert_eq!(
            classify_in_solid_exact(&model, &built.shape, Point::new(1.0, 1.0, 1.0), T).unwrap(),
            Containment::Out
        );
    }

    #[test]
    fn the_exact_classifier_refuses_an_open_boundary() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();
        assert!(classify_in_solid_exact(&model, &face, Point::ORIGIN, T).is_err());
    }

    #[test]
    fn distance_to_a_triangle_is_measured_from_the_nearest_part_of_it() {
        let t = [
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        // Above the interior: the plane distance.
        assert!((distance_to_triangle(Point::new(0.25, 0.25, 2.0), t) - 2.0).abs() < 1e-12);
        // Beyond a vertex: the distance to that vertex.
        assert!((distance_to_triangle(Point::new(-3.0, 0.0, 0.0), t) - 3.0).abs() < 1e-12);
        // On it: nothing.
        assert!(distance_to_triangle(Point::new(0.25, 0.25, 0.0), t) < 1e-12);
    }

    #[test]
    fn an_unusable_deflection_is_refused() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let bad = Deflection {
            chord: f64::NAN,
            ..Deflection::default()
        };
        assert!(classify_in_solid(&model, &built.shape, Point::ORIGIN, bad, T).is_err());
    }
}
