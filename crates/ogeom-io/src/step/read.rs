//! From parsed exchange structure to a living model.
//!
//! The reader walks every `MANIFOLD_SOLID_BREP` and rebuilds it bottom-up:
//! points, placements, curves and surfaces into geometry; vertices, edges,
//! loops, faces and shells into topology, shared exactly as the file shares
//! them — a vertex referenced by eight edges is one vertex here too, which is
//! what lets a closed shell close. Edge ranges are re-derived on this
//! kernel's own parameterizations from the vertex geometry, because STEP's
//! parameterizations are its own business and carrying them over blind is
//! how off-by-a-period bugs are born.
//!
//! What the reader does not understand it *counts*: every instance never
//! visited lands in the report's skipped table by keyword, and every
//! compromise — a face without a pcurve, a shell that does not close — is a
//! warning with the instance number in it. An import that succeeded with
//! three warnings is a different thing from one that succeeded, and the
//! report is what keeps the difference visible.

use super::parse::{Arg, Exchange, Instance};
use ogeom_algo::{make_edge_between, make_face_on, make_shell, make_solid, make_vertex, make_wire};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{
    BSplineCurve, CircleCurve, ConeSurface, Curve, CylinderSurface, EllipseCurve, LineCurve,
    PlanarCurve, PlaneSurface, SphereSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{
    Axis, Circle, Cone, Cylinder, Direction, Ellipse, Frame, KnotVector, Plane, Point, Sphere,
    Torus, Transform, Transform2, Vector, Vector2,
};
use ogeom_topo::{Location, Model, Shape};
use std::collections::{BTreeMap, HashMap, HashSet};

/// How far a plane or a cylinder read from a file extends past what anything
/// in the file uses. A face's trim is its wires; the surface's domain is only
/// the parameter window, and this one is generous without being unbounded.
const SURFACE_EXTENT: f64 = 1e5;

/// What an import brought in, and what it left behind.
#[derive(Debug, Default)]
pub struct StepReport {
    /// Millimetres per file unit, as the file's own unit section states it.
    pub scale_mm: f64,
    /// Instance keywords the reader never visited, with counts. Presentation,
    /// product structure and annotation live here by design; geometry landing
    /// here is a gap worth reading about.
    pub skipped: BTreeMap<String, usize>,
    /// Everything that imported less than perfectly, one line each.
    pub warnings: Vec<String>,
    /// The warning flood, counted: one entry per *kind* of imperfection,
    /// with how often it happened, the worst measured value where the kind
    /// measures one, and one exemplar entity id. A 776-warning community
    /// file summarises to a handful of lines a consumer can actually show;
    /// `warnings` keeps the full prose. Sorted by count, largest first.
    pub summary: Vec<WarningSummary>,
    /// Faces that read without a complete trim: an edge's boundary sat too
    /// far from the surface for any honest pcurve (beyond the one-millimetre
    /// healing cap), so the face will refuse to triangulate. Deduplicated,
    /// in file order — the structured form of the warnings that name them,
    /// carrying the face itself so the instructed follow-up needs no search.
    /// `check` reports the same faces as broken from the model side. A
    /// refused id whose face never finished building has nothing to act on
    /// and stays in the warnings alone.
    pub untrimmed_faces: Vec<UntrimmedFace>,
}

/// One kind of imperfect import, counted rather than repeated.
#[derive(Debug, Clone)]
pub struct WarningSummary {
    /// The kind, stable across runs: `"vertex-miss"`, `"boundary-slop"`,
    /// `"fit-short"`, `"untrimmed"`.
    pub kind: &'static str,
    /// How many times it happened.
    pub count: usize,
    /// The worst measured value among them — a distance, for every kind
    /// that measures one; zero where none applies.
    pub worst: f64,
    /// One entity id to look at first.
    pub exemplar: u64,
}

/// A face the reader could not trim, with the shape to act on.
#[derive(Debug, Clone)]
pub struct UntrimmedFace {
    /// The file's id for the face — the same id the warnings name.
    pub entity: u64,
    /// The face as built; hand it to `ogeom_heal::fix_face_pcurves` with
    /// the cap the situation deserves.
    pub face: Shape,
}

/// A read exchange file: the model, the solids found, and the report.
#[derive(Debug)]
pub struct StepImport {
    /// The document everything was built into: the model, plus the file's
    /// product structure, names and colours.
    pub document: ogeom_doc::Document,
    /// One shape per `MANIFOLD_SOLID_BREP`, in file order.
    pub solids: Vec<Shape>,
    /// What happened along the way.
    pub report: StepReport,
}

/// Read a STEP exchange file's B-rep content.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the file does
/// not parse, contains no solid, or a solid's structure is broken in a way
/// topology cannot represent. Faces the reader cannot complete become
/// warnings, not errors — the report says exactly what was compromised.
pub fn read_step(text: &str, tol: Tolerances) -> OgeomResult<StepImport> {
    let exchange = super::parse::parse(text)?;
    let mut reader = Reader {
        exchange: &exchange,
        model: Model::new(),
        report: StepReport {
            scale_mm: 1.0,
            ..StepReport::default()
        },
        angle_scale: 1.0,
        visited: HashSet::new(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        callout_index: HashMap::new(),
        pcurves: HashMap::new(),
        untrimmed_ids: Vec::new(),
        tallies: HashMap::new(),
        cdsr_of_nauo: None,
        properties_of_definition: None,
        sdrs_of_property: None,
        tol,
    };
    reader.report.scale_mm = reader.unit_scale();
    reader.angle_scale = reader.angle_unit_scale();

    let mut solids = Vec::new();
    let mut by_msb: HashMap<u64, Shape> = HashMap::new();
    let mut ids: Vec<u64> = exchange
        .data
        .iter()
        .filter(|(_, inst)| inst.part("MANIFOLD_SOLID_BREP").is_some())
        .map(|(id, _)| *id)
        .collect();
    ids.sort_unstable();
    let total = ids.len() as u64;
    for (done, id) in ids.into_iter().enumerate() {
        ogeom_core::progress::checkpoint()?;
        ogeom_core::progress::stage_at("step: solid", done as u64 + 1, total);
        let solid = reader.solid(id)?;
        by_msb.insert(id, solid.clone());
        solids.push(solid);
    }
    // The tallies fold into the summary, largest first, ties by kind so the
    // order is the file's and not the map's.
    {
        let mut entries: Vec<WarningSummary> = reader
            .tallies
            .drain()
            .map(|(kind, (count, worst, exemplar))| WarningSummary {
                kind,
                count,
                worst,
                exemplar,
            })
            .collect();
        entries.sort_by(|a, b| b.count.cmp(&a.count).then(a.kind.cmp(b.kind)));
        reader.report.summary = entries;
    }
    // The noted refusals resolve to the faces themselves, now that they
    // exist: the same cache the shells were assembled from answers by the
    // very id the warnings name.
    for id in std::mem::take(&mut reader.untrimmed_ids) {
        if let Some(face) = reader.faces.get(&id) {
            reader.report.untrimmed_faces.push(UntrimmedFace {
                entity: id,
                face: face.clone(),
            });
        }
    }
    if solids.is_empty() {
        ogeom_bail!(
            Construction,
            "the exchange file contains no MANIFOLD_SOLID_BREP to read"
        );
    }
    let document = reader.document(&by_msb, &solids)?;

    // Everything never visited, counted by its leading keyword.
    for (id, instance) in &exchange.data {
        if !reader.visited.contains(id) {
            *reader
                .report
                .skipped
                .entry(instance.keyword().to_owned())
                .or_default() += 1;
        }
    }

    Ok(StepImport {
        document,
        solids,
        report: reader.report,
    })
}

/// An edge as built: the shape, its curve and range, and whether STEP's edge
/// direction runs against the curve.
type BuiltEdge = (Shape, Curve, (f64, f64), bool);

/// A pcurve derived ahead of the face that needs it.
///
/// Deriving one is a pure function of a curve, its range and a surface — no
/// model, no order — and on a real assembly it is 95% of the time spent
/// building solids. So it is done for a whole solid at once, off the walk
/// that attaches it.
enum PreparedPcurve {
    /// The projection had a closed form.
    Exact(PlanarCurve),
    /// It did not, and this is the fit, with what the fit cost.
    Fitted {
        curve: PlanarCurve,
        error: f64,
        met: bool,
        worst_off: f64,
        warning: Option<String>,
    },
    /// Neither worked; the face gets the warning the walk would have made.
    Refused(String),
}

struct Reader<'a> {
    exchange: &'a Exchange,
    model: Model,
    report: StepReport,
    /// Radians per file angle unit — degrees are common.
    angle_scale: f64,
    visited: HashSet<u64>,
    vertices: HashMap<u64, Shape>,
    edges: HashMap<u64, BuiltEdge>,
    faces: HashMap<u64, Shape>,
    /// STEP id → index into the callouts just built, for the views pass.
    callout_index: HashMap<u64, usize>,
    /// `(face, edge)` → the pcurve already derived for it, from the parallel
    /// pass at the head of each solid.
    pcurves: HashMap<(u64, u64), PreparedPcurve>,
    /// Faces noted untrimmed, by file id; resolved to shapes once the read
    /// is far enough along for the shapes to exist.
    untrimmed_ids: Vec<u64>,
    /// kind → (count, worst, exemplar), folded into the report's summary.
    tallies: HashMap<&'static str, (usize, f64, u64)>,
    /// Usage → its `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`, built once.
    /// The lookup used to rescan the whole exchange per assembly edge —
    /// O(usages × entities), 333 × 457 k on one reporting assembly, which
    /// was three quarters of the entire document build.
    cdsr_of_nauo: Option<HashMap<u64, u64>>,
    /// Definition → its `PROPERTY_DEFINITION`s, and property → its
    /// `SHAPE_DEFINITION_REPRESENTATION`s, built together once. The datum
    /// target lookup used to rescan the whole exchange per property *per
    /// target* — the same quadratic shape the assembly index retired, one
    /// storey deeper. Each list ascends by id, so whichever entry answers
    /// is the one the old scan would have reached first.
    properties_of_definition: Option<HashMap<u64, Vec<u64>>>,
    sdrs_of_property: Option<HashMap<u64, Vec<u64>>>,
    tol: Tolerances,
}

impl Reader<'_> {
    fn instance(&mut self, id: u64) -> OgeomResult<&'_ Instance> {
        self.visited.insert(id);
        self.exchange.data.get(&id).ok_or_else(|| {
            ogeom_core::ogeom_err!(
                Construction,
                "the file references #{id}, which does not exist"
            )
        })
    }

    fn args(&mut self, id: u64, keyword: &str) -> OgeomResult<Vec<Arg>> {
        let instance = self.instance(id)?;
        let Some(args) = instance.part(keyword) else {
            ogeom_bail!(
                Construction,
                "#{id} is {}, where {keyword} was needed",
                instance.keyword()
            );
        };
        Ok(args.to_vec())
    }

    // --- units ---------------------------------------------------------------

    /// The unit instances the representation context actually assigns.
    ///
    /// A file may carry both a radian and a degree — definition and
    /// conversion — and which applies is not a matter of existence but of
    /// assignment: `GLOBAL_UNIT_ASSIGNED_CONTEXT` lists the ones in force.
    fn assigned_units(&self) -> Vec<u64> {
        // A file may carry several unit contexts — inches for the model,
        // millimetres for its annotation sheet — and the shapes' own
        // representations name the one their coordinates mean. That context
        // goes first; the rest follow in entity order, so the answer never
        // depends on how a map happens to iterate.
        let mut cited: Vec<u64> = self
            .exchange
            .data
            .values()
            .filter(|inst| {
                inst.parts
                    .iter()
                    .any(|(k, _)| k.ends_with("SHAPE_REPRESENTATION"))
            })
            .filter_map(|inst| {
                inst.parts
                    .iter()
                    .find(|(k, _)| k.ends_with("SHAPE_REPRESENTATION"))
                    .and_then(|(_, args)| args.last())
                    .and_then(Arg::reference)
            })
            .collect();
        cited.sort_unstable();
        cited.dedup();

        let mut contexts: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.part("GLOBAL_UNIT_ASSIGNED_CONTEXT").is_some())
            .map(|(id, _)| *id)
            .collect();
        contexts.sort_unstable();
        contexts.sort_by_key(|id| !cited.contains(id));

        let mut out = Vec::new();
        for id in contexts {
            if let Some(instance) = self.exchange.data.get(&id)
                && let Some(args) = instance.part("GLOBAL_UNIT_ASSIGNED_CONTEXT")
            {
                for arg in args.iter().filter_map(Arg::list).flatten() {
                    if let Some(r) = arg.reference() {
                        out.push(r);
                    }
                }
            }
        }
        out
    }

    /// Millimetres per file length unit, from the units the context assigns.
    fn unit_scale(&mut self) -> f64 {
        let assigned = self.assigned_units();
        for id in assigned {
            let Some(instance) = self.exchange.data.get(&id) else {
                continue;
            };
            let instance = instance.clone();
            if instance.part("LENGTH_UNIT").is_none() {
                continue;
            }
            if let Some(args) = instance.part("SI_UNIT") {
                return match args.first() {
                    Some(Arg::Enum(prefix)) => match prefix.as_str() {
                        "MILLI" => 1.0,
                        "CENTI" => 10.0,
                        "DECI" => 100.0,
                        "KILO" => 1e6,
                        "MICRO" => 1e-3,
                        other => {
                            self.report
                                .warnings
                                .push(format!("#{id}: unknown SI prefix {other}; taking metres"));
                            1000.0
                        }
                    },
                    _ => 1000.0,
                };
            }
            if let Some(args) = instance.part("CONVERSION_BASED_UNIT")
                && let Some(measure) = args.get(1).and_then(Arg::reference)
                && let Some(inner) = self.exchange.data.get(&measure)
            {
                let factor = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(|a| match a {
                        Arg::Typed(k, v) if k == "LENGTH_MEASURE" => {
                            v.first().and_then(Arg::number)
                        }
                        _ => None,
                    });
                let base = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(Arg::reference)
                    .and_then(|b| self.exchange.data.get(&b))
                    .and_then(|u| u.part("SI_UNIT"))
                    .map_or(1000.0, |args| match args.first() {
                        Some(Arg::Enum(p)) if p == "MILLI" => 1.0,
                        Some(Arg::Enum(p)) if p == "CENTI" => 10.0,
                        _ => 1000.0,
                    });
                if let Some(f) = factor {
                    return f * base;
                }
            }
        }
        self.report
            .warnings
            .push("no length unit found; taking millimetres".to_owned());
        1.0
    }

    /// Radians per file angle unit — the same dance as length, for the files
    /// that measure their cones in degrees.
    fn angle_unit_scale(&mut self) -> f64 {
        let assigned = self.assigned_units();
        for id in assigned {
            let Some(instance) = self.exchange.data.get(&id) else {
                continue;
            };
            let instance = instance.clone();
            if instance.part("PLANE_ANGLE_UNIT").is_none() {
                continue;
            }
            if instance.part("SI_UNIT").is_some() {
                return 1.0;
            }
            if let Some(args) = instance.part("CONVERSION_BASED_UNIT")
                && let Some(measure) = args.get(1).and_then(Arg::reference)
                && let Some(inner) = self.exchange.data.get(&measure)
            {
                let factor = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(|a| match a {
                        Arg::Typed(k, v) if k == "PLANE_ANGLE_MEASURE" => {
                            v.first().and_then(Arg::number)
                        }
                        _ => None,
                    });
                if let Some(f) = factor {
                    return f;
                }
            }
        }
        1.0
    }

    // --- geometry ------------------------------------------------------------

    fn point(&mut self, id: u64) -> OgeomResult<Point> {
        let args = self.args(id, "CARTESIAN_POINT")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a point without coordinates");
        };
        let scale = self.report.scale_mm;
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0) * scale;
        Ok(Point::new(value(0), value(1), value(2)))
    }

    fn direction(&mut self, id: u64) -> OgeomResult<Direction> {
        let args = self.args(id, "DIRECTION")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            ogeom_bail!(Construction, "#{id}: a direction without components");
        };
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0);
        Direction::new(Vector::new(value(0), value(1), value(2)), self.tol)
    }

    /// An `AXIS2_PLACEMENT_3D` as a frame, with the standard's defaults for
    /// what the file leaves out.
    fn frame(&mut self, id: u64) -> OgeomResult<Frame> {
        let args = self.args(id, "AXIS2_PLACEMENT_3D")?;
        let origin = args
            .get(1)
            .and_then(Arg::reference)
            .map(|r| self.point(r))
            .transpose()?
            .unwrap_or(Point::ORIGIN);
        let z = args
            .get(2)
            .and_then(Arg::reference)
            .map(|r| self.direction(r))
            .transpose()?
            .unwrap_or(Direction::Z);
        let x = match args.get(3).and_then(Arg::reference) {
            Some(r) => self.direction(r)?,
            None => Direction::from_cross(z.vector(), Vector::new(0.31, 0.52, 0.8), self.tol)?,
        };
        Frame::new(origin, z, x, self.tol)
    }

    fn surface(&mut self, id: u64) -> OgeomResult<Option<SurfaceGeometry>> {
        let (keyword, args) = {
            let instance = self.instance(id)?;
            (
                instance.keyword().to_owned(),
                instance
                    .parts
                    .first()
                    .map(|(_, a)| a.clone())
                    .unwrap_or_default(),
            )
        };
        let scale = self.report.scale_mm;
        let radius_arg = |args: &[Arg], i: usize| args.get(i).and_then(Arg::number);
        // B-spline surfaces arrive two ways: a simple instance with every
        // attribute in one list, or a complex instance whose parts each
        // carry their own slice — the rational form always the latter.
        {
            let (base, knots_part, weights) = {
                let instance = self.instance(id)?;
                (
                    instance.part("B_SPLINE_SURFACE").map(<[Arg]>::to_vec),
                    instance
                        .part("B_SPLINE_SURFACE_WITH_KNOTS")
                        .map(<[Arg]>::to_vec),
                    instance
                        .part("RATIONAL_B_SPLINE_SURFACE")
                        .map(<[Arg]>::to_vec),
                )
            };
            if let Some(kp) = knots_part {
                let (degrees, grid_arg, mults_knots) = if let Some(base) = base {
                    (
                        (base.first().cloned(), base.get(1).cloned()),
                        base.get(2).cloned(),
                        (
                            kp.first().cloned(),
                            kp.get(1).cloned(),
                            kp.get(2).cloned(),
                            kp.get(3).cloned(),
                        ),
                    )
                } else {
                    (
                        (kp.get(1).cloned(), kp.get(2).cloned()),
                        kp.get(3).cloned(),
                        (
                            kp.get(8).cloned(),
                            kp.get(9).cloned(),
                            kp.get(10).cloned(),
                            kp.get(11).cloned(),
                        ),
                    )
                };
                return self
                    .bspline_surface(id, degrees, grid_arg, mults_knots, weights)
                    .map(Some);
            }
        }
        let out = match keyword.as_str() {
            "PLANE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                Some(
                    PlaneSurface::over(
                        Plane::new(frame),
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "CYLINDRICAL_SURFACE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let radius = radius_arg(&args, 2).unwrap_or(0.0) * scale;
                Some(
                    CylinderSurface::new(
                        Cylinder::new(frame, radius, self.tol)?,
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "CONICAL_SURFACE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let radius = radius_arg(&args, 2).unwrap_or(0.0) * scale;
                let angle = radius_arg(&args, 3).unwrap_or(0.0) * self.angle_scale;
                Some(
                    ConeSurface::new(
                        Cone::new(frame, radius, angle, self.tol)?,
                        (-SURFACE_EXTENT, SURFACE_EXTENT),
                    )?
                    .into(),
                )
            }
            "SPHERICAL_SURFACE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let radius = radius_arg(&args, 2).unwrap_or(0.0) * scale;
                Some(SphereSurface::new(Sphere::new(frame, radius, self.tol)?).into())
            }
            "TOROIDAL_SURFACE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let major = radius_arg(&args, 2).unwrap_or(0.0) * scale;
                let minor = radius_arg(&args, 3).unwrap_or(0.0) * scale;
                Some(TorusSurface::new(Torus::new(frame, major, minor, self.tol)?).into())
            }
            other => {
                self.report.warnings.push(format!(
                    "#{id}: surface kind {other} is not read yet; its face is skipped"
                ));
                None
            }
        };
        Ok(out)
    }

    /// Expand STEP's multiplicity-compressed knots.
    fn expand_knots(mults: Option<Arg>, knots: Option<Arg>) -> Vec<f64> {
        let mut out = Vec::new();
        let (Some(Arg::List(mults)), Some(Arg::List(knots))) = (mults, knots) else {
            return out;
        };
        for (m, k) in mults.iter().zip(&knots) {
            let count = m.number().unwrap_or(1.0);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let count = count as usize;
            for _ in 0..count {
                out.push(k.number().unwrap_or(0.0));
            }
        }
        out
    }

    #[allow(clippy::type_complexity)]
    fn bspline_surface(
        &mut self,
        id: u64,
        degrees: (Option<Arg>, Option<Arg>),
        grid_arg: Option<Arg>,
        mults_knots: (Option<Arg>, Option<Arg>, Option<Arg>, Option<Arg>),
        weights: Option<Vec<Arg>>,
    ) -> OgeomResult<SurfaceGeometry> {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let deg = |a: Option<Arg>| a.and_then(|x| x.number()).unwrap_or(1.0) as usize;
        let (u_degree, v_degree) = (deg(degrees.0), deg(degrees.1));
        let Some(Arg::List(rows)) = grid_arg else {
            ogeom_bail!(
                Construction,
                "#{id}: a b-spline surface without control points"
            );
        };
        let mut points = Vec::new();
        let (mut u_count, mut v_count) = (0, 0);
        for row in &rows {
            let Some(cells) = row.list() else { continue };
            u_count += 1;
            v_count = cells.len();
            for cell in cells {
                points.push(self.point(cell.reference().unwrap_or(0))?);
            }
        }
        let u_knots = KnotVector::new(Self::expand_knots(mults_knots.0, mults_knots.2), u_degree)?;
        let v_knots = KnotVector::new(Self::expand_knots(mults_knots.1, mults_knots.3), v_degree)?;
        let surface = if let Some(weights) = weights {
            let flat: Vec<f64> = weights
                .first()
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::list)
                .flatten()
                .filter_map(Arg::number)
                .collect();
            let weighted: Vec<ogeom_math::Weighted<Point>> = points
                .iter()
                .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                .map(|(p, w)| ogeom_math::Weighted::new(*p, *w, self.tol))
                .collect::<OgeomResult<_>>()?;
            ogeom_geom::BSplineSurface::rational(
                u_knots,
                v_knots,
                ogeom_math::ControlGrid::new(weighted, u_count, v_count)?,
            )?
        } else {
            ogeom_geom::BSplineSurface::new(
                u_knots,
                v_knots,
                &ogeom_math::ControlGrid::new(points, u_count, v_count)?,
                self.tol,
            )?
        };
        Ok(surface.into())
    }

    fn curve(&mut self, id: u64) -> OgeomResult<Option<Curve>> {
        let (keyword, args) = {
            let instance = self.instance(id)?;
            (
                instance.keyword().to_owned(),
                instance
                    .parts
                    .first()
                    .map(|(_, a)| a.clone())
                    .unwrap_or_default(),
            )
        };
        let scale = self.report.scale_mm;
        {
            let (base, kp, weights) = {
                let instance = self.instance(id)?;
                (
                    instance.part("B_SPLINE_CURVE").map(<[Arg]>::to_vec),
                    instance
                        .part("B_SPLINE_CURVE_WITH_KNOTS")
                        .map(<[Arg]>::to_vec),
                    instance
                        .part("RATIONAL_B_SPLINE_CURVE")
                        .map(<[Arg]>::to_vec),
                )
            };
            if let (Some(base), Some(kp), Some(weights)) = (base, kp, weights) {
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let degree = base.first().and_then(Arg::number).unwrap_or(1.0) as usize;
                let control: Vec<Point> = base
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                let knots = KnotVector::new(
                    Self::expand_knots(kp.first().cloned(), kp.get(1).cloned()),
                    degree,
                )?;
                let flat: Vec<f64> = weights
                    .first()
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::number)
                    .collect();
                let weighted: Vec<ogeom_math::Weighted<Point>> = control
                    .iter()
                    .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                    .map(|(p, w)| ogeom_math::Weighted::new(*p, *w, self.tol))
                    .collect::<OgeomResult<_>>()?;
                return Ok(Some(BSplineCurve::rational(knots, weighted)?.into()));
            }
        }
        let out: Option<Curve> = match keyword.as_str() {
            "LINE" => {
                let through = self.point(args[1].reference().unwrap_or(0))?;
                // The vector's magnitude scales STEP's parameter; ranges here
                // are re-derived from vertex geometry, so only the direction
                // matters.
                let vector = self.args(args[2].reference().unwrap_or(0), "VECTOR")?;
                let direction = self.direction(vector[1].reference().unwrap_or(0))?;
                Some(
                    LineCurve::new(Axis {
                        location: through,
                        direction,
                    })
                    .into(),
                )
            }
            "CIRCLE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let radius = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(CircleCurve::new(Circle::new(frame, radius, self.tol)?).into())
            }
            "ELLIPSE" => {
                let frame = self.frame(args[1].reference().unwrap_or(0))?;
                let a = args.get(2).and_then(Arg::number).unwrap_or(0.0) * scale;
                let b = args.get(3).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(EllipseCurve::new(Ellipse::new(frame, a, b, self.tol)?).into())
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = args.get(1).and_then(Arg::number).unwrap_or(0.0);
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let degree = degree as usize;
                let control: Vec<Point> = args
                    .get(2)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .map(|r| self.point(r))
                    .collect::<OgeomResult<_>>()?;
                let mults = args.get(6).and_then(Arg::list).unwrap_or(&[]).to_vec();
                let knots = args.get(7).and_then(Arg::list).unwrap_or(&[]).to_vec();
                let mut expanded = Vec::new();
                for (m, k) in mults.iter().zip(&knots) {
                    let count = m.number().unwrap_or(1.0);
                    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                    let count = count as usize;
                    for _ in 0..count {
                        expanded.push(k.number().unwrap_or(0.0));
                    }
                }
                Some(
                    BSplineCurve::new(KnotVector::new(expanded, degree)?, control, self.tol)?
                        .into(),
                )
            }
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
                // A curve dressed in its surface associations — what every
                // exporter derived from the reference kernel writes for
                // every edge. The 3D curve is the first argument; the
                // pcurve list is advisory and this reader re-derives its
                // own, so unwrapping is the whole job.
                self.curve(args.get(1).and_then(Arg::reference).unwrap_or(0))?
            }
            other => {
                self.report.warnings.push(format!(
                    "#{id}: curve kind {other} is not read yet; its edge is skipped"
                ));
                None
            }
        };
        Ok(out)
    }

    /// The parameter of a point on one of this kernel's curves.
    fn parameter_of(&self, curve: &Curve, p: Point) -> Option<f64> {
        match curve {
            Curve::Line(line) => {
                let axis = line.axis();
                Some((p - axis.location).dot(axis.direction.vector()))
            }
            Curve::Circle(c) => {
                let local = c.circle().frame().to_local(p);
                Some(local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU))
            }
            Curve::Ellipse(e) => {
                let local = e.ellipse().frame().to_local(p);
                Some(
                    (local.y / e.ellipse().minor_radius())
                        .atan2(local.x / e.ellipse().major_radius())
                        .rem_euclid(core::f64::consts::TAU),
                )
            }
            Curve::BSpline(_) => None,
            _ => None,
        }
    }

    // --- topology ------------------------------------------------------------

    fn vertex(&mut self, id: u64) -> OgeomResult<Shape> {
        if let Some(shape) = self.vertices.get(&id) {
            return Ok(shape.clone());
        }
        let args = self.args(id, "VERTEX_POINT")?;
        let point = self.point(args[1].reference().unwrap_or(0))?;
        let shape = make_vertex(&mut self.model, point).shape;
        self.vertices.insert(id, shape.clone());
        Ok(shape)
    }

    /// An `EDGE_CURVE`, built once, curve-forward.
    ///
    /// Returns the edge, its curve, its range, and whether STEP's edge
    /// direction runs *against* the curve — which each use folds into its
    /// own orientation.
    fn edge(&mut self, id: u64) -> OgeomResult<Option<BuiltEdge>> {
        if let Some(found) = self.edges.get(&id) {
            return Ok(Some(found.clone()));
        }
        let args = self.args(id, "EDGE_CURVE")?;
        let v1 = args[1].reference().unwrap_or(0);
        let v2 = args[2].reference().unwrap_or(0);
        let Some(curve) = self.curve(args[3].reference().unwrap_or(0))? else {
            return Ok(None);
        };
        let same_sense = !args.get(4).is_some_and(|a| a.is_enum("F"));

        let p1 = {
            let vargs = self.args(v1, "VERTEX_POINT")?;
            self.point(vargs[1].reference().unwrap_or(0))?
        };
        let p2 = {
            let vargs = self.args(v2, "VERTEX_POINT")?;
            self.point(vargs[1].reference().unwrap_or(0))?
        };

        // The edge is built along the curve's own parameter; a STEP edge
        // running the other way is flagged, and every use composes the flag
        // into its orientation.
        let (start, end, flipped) = if same_sense {
            (p1, p2, false)
        } else {
            (p2, p1, true)
        };
        let (range, closed) = match (
            self.parameter_of(&curve, start),
            self.parameter_of(&curve, end),
        ) {
            (Some(a), Some(b)) => {
                let period = if curve.is_periodic() {
                    let (lo, hi) = curve.domain();
                    hi - lo
                } else {
                    0.0
                };
                if v1 == v2 {
                    ((a, a + if period > 0.0 { period } else { 0.0 }), true)
                } else if period > 0.0 && b <= a + self.tol.parametric() {
                    ((a, b + period), false)
                } else {
                    ((a, b), false)
                }
            }
            _ => {
                // No closed-form inversion: take the curve's own domain and
                // hold the endpoints to it.
                let (lo, hi) = curve.domain();
                let head = curve.point_at(lo, self.tol)?;
                let tail = curve.point_at(hi, self.tol)?;
                if head.distance(start) > self.tol.confusion() * 10.0
                    || tail.distance(end) > self.tol.confusion() * 10.0
                {
                    self.report.warnings.push(format!(
                        "#{id}: edge endpoints sit {:.2e} and {:.2e} from its \
                         curve's ends; the curve's own domain was taken",
                        head.distance(start),
                        tail.distance(end)
                    ));
                }
                ((lo, hi), v1 == v2)
            }
        };
        let _ = closed;

        let (vlo, vhi) = if flipped {
            (self.vertex(v2)?, self.vertex(v1)?)
        } else {
            (self.vertex(v1)?, self.vertex(v2)?)
        };
        // Real files are imprecise, and NIST's own readme says so of these.
        // Where the curve's end misses its vertex by more than the default
        // tolerance, the vertex's tolerance *grows* to state the gap — the
        // data model's growing tolerances are exactly for this, and the
        // warning keeps the healing visible.
        for (vertex, t) in [(&vlo, range.0), (&vhi, range.1)] {
            let end = curve.point_at(t, self.tol)?;
            let stated = {
                let Some(node) = self.model.node(vertex) else {
                    continue;
                };
                let Some(data) = node.data().as_vertex() else {
                    continue;
                };
                (data.point, data.tolerance)
            };
            let gap = end.distance(stated.0);
            if gap > stated.1.get() {
                if let Some(node) = self.model.node_mut(vertex)
                    && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(gap + self.tol.confusion());
                }
                self.warn_vertex_miss(id, gap);
            }
        }
        let shape =
            make_edge_between(&mut self.model, curve.clone(), range, &vlo, &vhi, self.tol)?.shape;
        let entry = (shape, curve, range, flipped);
        self.edges.insert(id, entry.clone());
        Ok(Some(entry))
    }

    fn face(&mut self, id: u64) -> OgeomResult<Option<Shape>> {
        if let Some(shape) = self.faces.get(&id) {
            return Ok(Some(shape.clone()));
        }
        let args = self.args(id, "ADVANCED_FACE")?;
        let bounds: Vec<u64> = args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        let Some(surface) = self.surface(args[2].reference().unwrap_or(0))? else {
            return Ok(None);
        };
        let face_forward = !args.get(3).is_some_and(|a| a.is_enum("F"));
        let surface_id = self.model.geometry_mut().add_surface(surface.clone());

        // Outer bound first, so the face's first wire is its outer ring.
        let mut ordered = bounds.clone();
        ordered.sort_by_key(|b| {
            self.exchange
                .data
                .get(b)
                .map_or(1, |i| i32::from(i.part("FACE_OUTER_BOUND").is_none()))
        });

        // Which edges this face uses twice: those are seams, and get both
        // sides' pcurves.
        let mut edge_uses: HashMap<u64, usize> = HashMap::new();
        for &bound in &ordered {
            let bargs = self.bound_args(bound)?;
            if self.instance(bargs.0)?.part("VERTEX_LOOP").is_some() {
                continue;
            }
            let loop_args = self.args(bargs.0, "EDGE_LOOP")?;
            for oe in loop_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                if let Some(oe_id) = oe.reference() {
                    let oargs = self.args(oe_id, "ORIENTED_EDGE")?;
                    if let Some(e) = oargs.get(3).and_then(Arg::reference) {
                        *edge_uses.entry(e).or_default() += 1;
                    }
                }
            }
        }

        let mut wires = Vec::new();
        let mut annotated: HashSet<u64> = HashSet::new();
        for &bound in &ordered {
            let (loop_id, bound_forward) = self.bound_args(bound)?;
            if let Some(vertex_loop) = {
                let instance = self.instance(loop_id)?;
                instance.part("VERTEX_LOOP").map(<[Arg]>::to_vec)
            } {
                // A loop of one vertex: a pole or an apex. It has no edges,
                // but it still bounds the face in parameter space — as a
                // degenerate edge running across the chart at the row the
                // point collapses to, exactly as native cones and spheres
                // are built.
                let vertex_id = vertex_loop.get(1).and_then(Arg::reference).unwrap_or(0);
                let vargs = self.args(vertex_id, "VERTEX_POINT")?;
                let at = self.point(vargs[1].reference().unwrap_or(0))?;
                let vertex = self.vertex(vertex_id)?;
                let projection = ogeom_algo::project_on_surface(&surface, at, 32, self.tol)?;
                let ((ua, ub), _) = surface.domain();
                let row = projection.parameters.1;
                let mut data = ogeom_topo::EdgeData::new();
                data.degenerate = true;
                let edge = self
                    .model
                    .add_edge(data, &[vertex.clone(), vertex.clone()])?;
                let pcurve: PlanarCurve = ogeom_geom::Line2d::segment(
                    ogeom_math::Point2::new(ua, row),
                    ogeom_math::Point2::new(ub, row),
                    self.tol,
                )?
                .into();
                ogeom_algo::attach_pcurve(
                    &mut self.model,
                    &edge,
                    pcurve,
                    surface_id,
                    Location::identity(),
                    (0.0, ub - ua),
                )?;
                wires.push(make_wire(&mut self.model, &[edge], self.tol)?.shape);
                continue;
            }
            let loop_args = self.args(loop_id, "EDGE_LOOP")?;
            let mut uses: Vec<(Shape, u64)> = Vec::new();
            for oe in loop_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                let Some(oe_id) = oe.reference() else {
                    continue;
                };
                let oargs = self.args(oe_id, "ORIENTED_EDGE")?;
                let Some(edge_id) = oargs.get(3).and_then(Arg::reference) else {
                    continue;
                };
                let forward = !oargs.get(4).is_some_and(|a| a.is_enum("F"));
                let Some((shape, curve, range, flipped)) = self.edge(edge_id)? else {
                    self.report.warnings.push(format!(
                        "#{id}: a bound references unreadable edge #{edge_id}; \
                         the face is skipped"
                    ));
                    return Ok(None);
                };
                // The use's direction composes the loop's, the bound's and
                // the edge-against-curve flag.
                let mut use_forward = forward == bound_forward;
                if flipped {
                    use_forward = !use_forward;
                }
                let placed = if use_forward {
                    shape.clone()
                } else {
                    shape.reversed()
                };
                uses.push((placed, edge_id));

                if annotated.insert(edge_id) {
                    self.attach_pcurves(
                        id,
                        edge_id,
                        &shape,
                        &curve,
                        range,
                        &surface,
                        surface_id,
                        edge_uses.get(&edge_id).copied().unwrap_or(1) > 1,
                    )?;
                }
            }
            if !bound_forward {
                uses.reverse();
            }
            let edges: Vec<Shape> = uses.into_iter().map(|(s, _)| s).collect();
            if edges.is_empty() {
                continue;
            }
            wires.push(make_wire(&mut self.model, &edges, self.tol)?.shape);
        }
        if wires.is_empty() {
            self.report
                .warnings
                .push(format!("#{id}: a face with no readable bounds is skipped"));
            return Ok(None);
        }

        // A periodic face bound only by closed rings — a cylinder band
        // between two circles — arrives without a seam edge, which is a
        // legitimate STEP shape and an open rectangle in this kernel's
        // chart. The seam is synthesised the way native cylinders build it:
        // one edge at the period join, appearing in the wire twice.
        if wires.len() == 2
            && surface.is_periodic_u()
            && let [(e_lo, _v_lo), (e_hi, _v_hi)] =
                closed_ring_edges(&self.model, &wires)?.as_slice()
        {
            {
                // The band construction is ogeom-algo's make_revolution_band —
                // one authority shared with the healer. Anything it refuses
                // (a spline ring, ring vertices at different angles, a
                // surface with no closed-form iso-curve) becomes a warning
                // and the raw bounds, not an error.
                match ogeom_algo::make_revolution_band(
                    &mut self.model,
                    &surface,
                    e_lo,
                    e_hi,
                    self.tol,
                ) {
                    Ok(built) => {
                        let shape = if face_forward {
                            built
                        } else {
                            built.reversed()
                        };
                        self.faces.insert(id, shape.clone());
                        return Ok(Some(shape));
                    }
                    Err(e) => {
                        self.report.warnings.push(format!(
                            "#{id}: no seam could be synthesised ({e}); the \
                             face may not triangulate"
                        ));
                    }
                }
            }
        }

        // A cone face bounded by a *single* closed ring: the apex is the
        // other boundary, and some files simply never write it — no vertex
        // loop, nothing. The geometry leaves one choice of what the face
        // means, so the apex is synthesised and the band built as if the
        // file had said so.
        if wires.len() == 1
            && matches!(surface, SurfaceGeometry::Cone(_))
            && let [(ring, _)] = closed_ring_edges(&self.model, &wires)?.as_slice()
        {
            match ogeom_algo::make_apex_band(&mut self.model, &surface, ring, self.tol) {
                Ok(built) => {
                    let shape = if face_forward {
                        built
                    } else {
                        built.reversed()
                    };
                    self.faces.insert(id, shape.clone());
                    return Ok(Some(shape));
                }
                Err(e) => {
                    self.report.warnings.push(format!(
                        "#{id}: no apex could be synthesised ({e}); the face \
                         may not triangulate"
                    ));
                }
            }
        }

        let built = make_face_on(&mut self.model, surface_id, &wires, self.tol)?.shape;
        let shape = if face_forward {
            built
        } else {
            built.reversed()
        };
        self.faces.insert(id, shape.clone());
        Ok(Some(shape))
    }

    /// A bound's loop and orientation, whichever of the two bound kinds it is.
    fn bound_args(&mut self, id: u64) -> OgeomResult<(u64, bool)> {
        let instance = self.instance(id)?;
        let args = instance
            .part("FACE_OUTER_BOUND")
            .or_else(|| instance.part("FACE_BOUND"))
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "#{id} is not a face bound"))?
            .to_vec();
        let loop_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let forward = !args.get(2).is_some_and(|a| a.is_enum("F"));
        Ok((loop_id, forward))
    }

    /// Attach this face's pcurve — or both seam sides — to an edge.
    #[allow(clippy::too_many_arguments)]
    fn attach_pcurves(
        &mut self,
        face_id: u64,
        edge_id: u64,
        edge: &Shape,
        curve: &Curve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
        surface_id: ogeom_topo::SurfaceId,
        seam: bool,
    ) -> OgeomResult<()> {
        let widen = |p: PlanarCurve| -> PlanarCurve {
            // A line pcurve evaluates anywhere; its stated domain must still
            // cover the edge's range, which for a wrapped circle runs past
            // one period.
            if let PlanarCurve::Line(l) = &p {
                let (lo, hi) = (l.domain().0.min(range.0), l.domain().1.max(range.1));
                if let Ok(wider) = ogeom_geom::Line2d::over(l.axis(), lo, hi) {
                    return wider.into();
                }
            }
            p
        };
        // Claimed rather than derived, where the pass at the head of the solid
        // already did it. The fallbacks below stay exactly as they were, for
        // the faces no pass covered — a face reached outside a solid walk, or
        // one whose preparation refused.
        if let Some(prepared) = self.pcurves.remove(&(face_id, edge_id)) {
            match prepared {
                PreparedPcurve::Exact(exact) => {
                    return self.record_pcurve(edge, widen(exact), surface_id, seam, range);
                }
                PreparedPcurve::Fitted {
                    curve: fitted,
                    error,
                    met,
                    worst_off,
                    warning,
                } => {
                    if let Some(w) = warning {
                        self.warn_slop(w, worst_off, face_id);
                    }
                    if !met {
                        self.warn_fit_short(face_id, error);
                    }
                    if worst_off > self.tol.confusion()
                        && let Some(node) = self.model.node_mut(edge)
                        && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                    {
                        data.tolerance = data.tolerance.widen_to(worst_off + self.tol.confusion());
                    }
                    return self.record_pcurve(edge, fitted, surface_id, seam, range);
                }
                PreparedPcurve::Refused(why) => {
                    self.report.warnings.push(why);
                    self.note_untrimmed(face_id);
                    return Ok(());
                }
            }
        }
        let pcurve =
            match ogeom_intersect::exact_pcurve_over(curve, range, surface, self.tol).map(widen) {
                Some(exact) => exact,
                None => {
                    // No closed form — a spline surface, or a combination the
                    // projection table lacks. The pcurve is *fitted at the
                    // curve's own parameters*: sample the edge, project each
                    // sample into the chart, fit the trace with the parameters
                    // held fixed, so same-parameter is preserved by construction
                    // and the reported error is the true chart deviation.
                    match crate::pcurves::fit_projected_pcurve(curve, range, surface, self.tol) {
                        Ok((fitted, error, met, worst_off, slop_warning)) => {
                            if let Some(w) = slop_warning {
                                self.warn_slop(w, worst_off, face_id);
                            }
                            if !met {
                                self.warn_fit_short(face_id, error);
                            }
                            // The edge provably sits `worst_off` from the surface
                            // it bounds; its tolerance grows to cover that, the
                            // same honesty the vertex ends get.
                            if worst_off > self.tol.confusion()
                                && let Some(node) = self.model.node_mut(edge)
                                && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                            {
                                data.tolerance =
                                    data.tolerance.widen_to(worst_off + self.tol.confusion());
                            }
                            fitted
                        }
                        Err(e) => {
                            self.report.warnings.push(format!(
                                "face #{face_id}: no pcurve for an edge on this \
                             surface ({e}); the face may not triangulate"
                            ));
                            self.note_untrimmed(face_id);
                            return Ok(());
                        }
                    }
                }
            };
        self.record_pcurve(edge, pcurve, surface_id, seam, range)
    }

    /// Count one occurrence of a warning kind toward the summary.
    fn tally(&mut self, kind: &'static str, measured: f64, exemplar: u64) {
        let entry = self.tallies.entry(kind).or_insert((0, 0.0, exemplar));
        entry.0 += 1;
        if measured > entry.1 {
            entry.1 = measured;
            entry.2 = exemplar;
        }
    }

    /// A curve end missing its vertex: the prose, and the count.
    fn warn_vertex_miss(&mut self, id: u64, gap: f64) {
        self.report.warnings.push(format!(
            "#{id}: a curve end misses its vertex by {gap:.2e}; the vertex \
             tolerance grew to say so"
        ));
        self.tally("vertex-miss", gap, id);
    }

    /// A fitted trim that stopped short of its target: prose and count.
    fn warn_fit_short(&mut self, face_id: u64, error: f64) {
        self.report.warnings.push(format!(
            "face #{face_id}: a projected pcurve fit stopped at {error:.2e}; \
             the face's mesh may sit that far off along this edge"
        ));
        self.tally("fit-short", error, face_id);
    }

    /// The file's own boundary slop, carried and counted.
    fn warn_slop(&mut self, prose: String, worst: f64, exemplar: u64) {
        self.report.warnings.push(prose);
        self.tally("boundary-slop", worst, exemplar);
    }

    /// Note a face left without a complete trim, once.
    fn note_untrimmed(&mut self, face_id: u64) {
        if self.untrimmed_ids.last() != Some(&face_id) {
            self.untrimmed_ids.push(face_id);
            self.tally("untrimmed", 0.0, face_id);
        }
    }

    /// Attach a derived pcurve, seaming it where the edge bounds the chart
    /// twice. The tail of `attach_pcurves`, shared with the prepared path.
    fn record_pcurve(
        &mut self,
        edge: &Shape,
        pcurve: PlanarCurve,
        surface_id: ogeom_topo::SurfaceId,
        seam: bool,
        range: (f64, f64),
    ) -> OgeomResult<()> {
        if seam {
            let Some(surface) = self.model.geometry().surface(surface_id).cloned() else {
                ogeom_bail!(Dangling, "the surface is not in this model");
            };
            let ((ua, ub), _) = surface.domain();
            let span = ub - ua;
            // One side is where the projection landed; the other is one
            // period over, whichever way stays within the chart.
            let mid = pcurve.point_at(f64::midpoint(range.0, range.1), self.tol)?;
            let shift = if mid.x - ua < span * 0.5 { span } else { -span };
            let other =
                pcurve.transformed(&Transform2::translation(Vector2::new(shift, 0.0)), self.tol)?;
            ogeom_algo::attach_seam(
                &mut self.model,
                edge,
                pcurve,
                other,
                surface_id,
                Location::identity(),
                range,
            )?;
        } else {
            ogeom_algo::attach_pcurve(
                &mut self.model,
                edge,
                pcurve,
                surface_id,
                Location::identity(),
                range,
            )?;
        }
        Ok(())
    }

    /// Derive every pcurve this solid's faces will want, in parallel.
    ///
    /// Deriving one is a pure function of a curve, its range and a surface —
    /// it reads no model and depends on no order — and measured on a 330-solid
    /// assembly it is 95% of the time spent building solids. The walk that
    /// follows attaches them, in file order, exactly as it did when it derived
    /// them itself.
    ///
    /// Resolving *what* to derive still runs on the walk's own thread, since
    /// it builds edges and vertices into the model. That part is cheap; it is
    /// the projection and the fitting that are not.
    ///
    /// Best-effort by design: anything this cannot resolve is simply left out
    /// of the table, and `attach_pcurves` derives it the old way. So a face
    /// shape this does not anticipate costs time, never correctness.
    fn prepare_pcurves(&mut self, face_ids: &[u64]) {
        struct Job {
            face: u64,
            edge: u64,
            curve: Curve,
            range: (f64, f64),
            /// Index into `surfaces`. The surface is held once per face, not
            /// once per edge: a B-spline patch owns its whole control grid,
            /// and cloning that per edge costs more than the projection it
            /// was cloned for.
            surface: usize,
        }
        let mut surfaces: Vec<SurfaceGeometry> = Vec::new();
        let mut jobs: Vec<Job> = Vec::new();
        let mut seen: HashSet<(u64, u64)> = HashSet::new();
        for &fid in face_ids {
            let Ok(args) = self.args(fid, "ADVANCED_FACE") else {
                continue;
            };
            let Some(surface) = args
                .get(2)
                .and_then(Arg::reference)
                .and_then(|sid| self.surface(sid).ok().flatten())
            else {
                continue;
            };
            surfaces.push(surface);
            let at = surfaces.len() - 1;
            let bounds: Vec<u64> = args
                .get(1)
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::reference)
                .collect();
            for bound in bounds {
                let Ok((loop_id, _)) = self.bound_args(bound) else {
                    continue;
                };
                let Ok(loop_args) = self.args(loop_id, "EDGE_LOOP") else {
                    continue;
                };
                let uses: Vec<u64> = loop_args
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .collect();
                for oe_id in uses {
                    let Ok(oargs) = self.args(oe_id, "ORIENTED_EDGE") else {
                        continue;
                    };
                    let Some(edge_id) = oargs.get(3).and_then(Arg::reference) else {
                        continue;
                    };
                    if !seen.insert((fid, edge_id)) {
                        continue;
                    }
                    if let Ok(Some((_, curve, range, _))) = self.edge(edge_id) {
                        jobs.push(Job {
                            face: fid,
                            edge: edge_id,
                            curve,
                            range,
                            surface: at,
                        });
                    }
                }
            }
        }
        // Below a handful of edges the threads cost more than the work; the
        // sequential path through `attach_pcurves` is already correct, so the
        // table is simply left empty and the walk derives them itself.
        if jobs.len() < 16 {
            return;
        }
        let tol = self.tol;
        let derived = ogeom_core::parallel::map_ordered(&jobs, |_, job| {
            let surface = &surfaces[job.surface];
            match ogeom_intersect::exact_pcurve_over(&job.curve, job.range, surface, tol) {
                Some(exact) => PreparedPcurve::Exact(exact),
                None => {
                    match crate::pcurves::fit_projected_pcurve(&job.curve, job.range, surface, tol)
                    {
                        Ok((curve, error, met, worst_off, warning)) => PreparedPcurve::Fitted {
                            curve,
                            error,
                            met,
                            worst_off,
                            warning,
                        },
                        Err(e) => PreparedPcurve::Refused(format!(
                            "face #{}: no pcurve for an edge on this surface ({e}); \
                         the face may not triangulate",
                            job.face
                        )),
                    }
                }
            }
        });
        for (job, pcurve) in jobs.iter().zip(derived) {
            self.pcurves.insert((job.face, job.edge), pcurve);
        }
    }

    /// The solid a `MANIFOLD_SOLID_BREP` names: its shell's faces, sewn.
    ///
    /// (The comment that stood here described a fitted pcurve, which is not
    /// what this builds; it had been copied from elsewhere.)
    fn solid(&mut self, id: u64) -> OgeomResult<Shape> {
        let args = self.args(id, "MANIFOLD_SOLID_BREP")?;
        let shell_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let shell_instance = self.instance(shell_id)?;
        let shell_args = shell_instance
            .part("CLOSED_SHELL")
            .or_else(|| shell_instance.part("OPEN_SHELL"))
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "#{shell_id} is not a shell"))?
            .to_vec();
        let face_ids: Vec<u64> = shell_args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        self.prepare_pcurves(&face_ids);
        let mut faces = Vec::new();
        for fid in face_ids {
            if let Some(face) = self.face(fid)? {
                faces.push(face);
            }
        }
        if faces.is_empty() {
            ogeom_bail!(Construction, "#{id}: a solid with no readable faces");
        }
        let shell = make_shell(&mut self.model, &faces)?.shape;
        if !ogeom_algo::is_shell_closed(&self.model, &shell)? {
            self.report.warnings.push(format!(
                "#{id}: the shell does not close as read; measures needing an \
                 inside will refuse it"
            ));
        }
        Ok(make_solid(&mut self.model, std::slice::from_ref(&shell))?.shape)
    }

    // --- product structure, names, colours -----------------------------------

    /// Assemble the document: the model, plus everything the file says about
    /// products, assemblies, placements and appearance.
    ///
    /// Takes the model out of the reader — geometry reading is over by the
    /// time structure is read. Structure that resists becomes a warning and a
    /// flat document, never an error: the geometry is already good, and a
    /// mangled product tree should not take it down.
    fn document(
        &mut self,
        by_msb: &HashMap<u64, Shape>,
        solids: &[Shape],
    ) -> OgeomResult<ogeom_doc::Document> {
        // The graph is walked before the model moves, because frames scale
        // through the reader's own unit handling.
        let structure = self.product_structure(by_msb);
        let colours = self.colours(by_msb);
        let pmi = self.pmi_of();

        let mut document = ogeom_doc::Document::over(std::mem::take(&mut self.model));
        match structure {
            Some(products) => self.build_products(&mut document, products),
            None => {
                for (i, solid) in solids.iter().enumerate() {
                    document.add_part(format!("solid-{i}"), solid.clone());
                }
            }
        }
        for (shape, colour) in colours {
            document.set_colour(&shape, colour);
        }
        *document.pmi_mut() = pmi;
        for view in self.views() {
            document.add_view(view);
        }
        Ok(document)
    }

    /// The file's product graph, or `None` when it has none worth the name.
    fn product_structure(&mut self, by_msb: &HashMap<u64, Shape>) -> Option<Vec<PdEntry>> {
        // PRODUCT_DEFINITION -> name, via formation and product.
        let mut pds: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.part("PRODUCT_DEFINITION").is_some())
            .filter(|(_, inst)| inst.part("PRODUCT_DEFINITION_RELATIONSHIP").is_none())
            .map(|(id, _)| *id)
            .collect();
        pds.sort_unstable();
        if pds.is_empty() {
            return None;
        }

        // SHAPE_DEFINITION_REPRESENTATION: definition (a PRODUCT_DEFINITION_SHAPE
        // over a PD or a usage) -> shape representation.
        let mut sr_of_pd: HashMap<u64, u64> = HashMap::new();
        let sdrs: Vec<(u64, u64)> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.part("SHAPE_DEFINITION_REPRESENTATION").is_some())
            .filter_map(|(id, _)| {
                let args = self.args(*id, "SHAPE_DEFINITION_REPRESENTATION").ok()?;
                Some((args.first()?.reference()?, args.get(1)?.reference()?))
            })
            .collect();
        for (pds_id, sr) in sdrs {
            if let Some(definition) = self.definition_of_shape(pds_id) {
                sr_of_pd.insert(definition, sr);
            }
        }

        // A product's solids may live one representation over: AP203 files
        // routinely tie the product to a bare axis representation and hang
        // the B-rep off it through a plain SHAPE_REPRESENTATION_RELATIONSHIP.
        // Only the plain ones are followed — the transformation-carrying kind
        // is an assembly edge, and following it would leak one product's
        // geometry into another.
        let mut linked: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut srrs: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| {
                inst.part("SHAPE_REPRESENTATION_RELATIONSHIP").is_some()
                    && inst
                        .part("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
                        .is_none()
            })
            .map(|(id, _)| *id)
            .collect();
        srrs.sort_unstable();
        for srr in srrs {
            let args = {
                let Ok(instance) = self.instance(srr) else {
                    continue;
                };
                let Some(args) = instance
                    .part("SHAPE_REPRESENTATION_RELATIONSHIP")
                    .or_else(|| instance.part("REPRESENTATION_RELATIONSHIP"))
                else {
                    continue;
                };
                args.to_vec()
            };
            if let (Some(a), Some(b)) = (
                args.get(2).and_then(Arg::reference),
                args.get(3).and_then(Arg::reference),
            ) {
                linked.entry(a).or_default().push(b);
                linked.entry(b).or_default().push(a);
            }
        }

        let mut entries = Vec::new();
        for pd in pds {
            let name = self.product_name(pd).unwrap_or_else(|| format!("#{pd}"));
            let mut shapes: Vec<Shape> = Vec::new();
            if let Some(&sr) = sr_of_pd.get(&pd) {
                let mut reps = vec![sr];
                reps.extend(linked.get(&sr).into_iter().flatten().copied());
                for rep in reps {
                    if let Some(items) = self.representation_items(rep) {
                        shapes.extend(items.iter().filter_map(|item| by_msb.get(item).cloned()));
                    }
                }
            }
            entries.push(PdEntry {
                pd,
                name,
                shapes,
                children: Vec::new(),
            });
        }

        // NEXT_ASSEMBLY_USAGE_OCCURRENCE: parent -> child, with the placement
        // recovered from the CONTEXT_DEPENDENT_SHAPE_REPRESENTATION over it.
        let mut nauos: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.part("NEXT_ASSEMBLY_USAGE_OCCURRENCE").is_some())
            .map(|(id, _)| *id)
            .collect();
        nauos.sort_unstable();
        let entry_of_pd: HashMap<u64, usize> =
            entries.iter().enumerate().map(|(i, e)| (e.pd, i)).collect();
        for nauo in nauos {
            let Ok(args) = self.args(nauo, "NEXT_ASSEMBLY_USAGE_OCCURRENCE") else {
                continue;
            };
            let (Some(parent), Some(child)) = (
                args.get(3).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) else {
                continue;
            };
            let name = args
                .get(5)
                .and_then(|a| match a {
                    Arg::Str(s) if !s.is_empty() => Some(s.clone()),
                    _ => None,
                })
                .or_else(|| {
                    self.args(nauo, "NEXT_ASSEMBLY_USAGE_OCCURRENCE")
                        .ok()
                        .and_then(|args| match args.get(1) {
                            Some(Arg::Str(s)) if !s.is_empty() => Some(s.clone()),
                            _ => None,
                        })
                });
            let child_sr = sr_of_pd.get(&child).copied();
            let at = self
                .usage_transform(nauo, child_sr)
                .unwrap_or(Transform::IDENTITY);
            if let Some(&at_index) = entry_of_pd.get(&parent) {
                entries[at_index].children.push((child, at, name));
            }
        }
        Some(entries)
    }

    /// The `PRODUCT_DEFINITION` (or usage) a `PRODUCT_DEFINITION_SHAPE` is over.
    fn definition_of_shape(&mut self, pds_id: u64) -> Option<u64> {
        let args = self.args(pds_id, "PRODUCT_DEFINITION_SHAPE").ok()?;
        args.get(2).and_then(Arg::reference)
    }

    /// A product definition's product name.
    fn product_name(&mut self, pd: u64) -> Option<String> {
        let formation = self
            .args(pd, "PRODUCT_DEFINITION")
            .ok()?
            .get(2)
            .and_then(Arg::reference)?;
        let product = self
            .args(formation, "PRODUCT_DEFINITION_FORMATION")
            .ok()?
            .first()
            .and_then(Arg::reference)
            .or_else(|| {
                self.args(formation, "PRODUCT_DEFINITION_FORMATION")
                    .ok()?
                    .get(2)
                    .and_then(Arg::reference)
            })?;
        let args = self.args(product, "PRODUCT").ok()?;
        match args.get(1).or_else(|| args.first()) {
            Some(Arg::Str(name)) if !name.is_empty() => Some(name.clone()),
            _ => None,
        }
    }

    /// A representation's item references.
    fn representation_items(&mut self, sr: u64) -> Option<Vec<u64>> {
        let instance = self.instance(sr).ok()?;
        let args = instance
            .part("SHAPE_REPRESENTATION")
            .or_else(|| instance.part("ADVANCED_BREP_SHAPE_REPRESENTATION"))
            .or_else(|| instance.part("REPRESENTATION"))?
            .to_vec();
        Some(
            args.get(1)?
                .list()?
                .iter()
                .filter_map(Arg::reference)
                .collect(),
        )
    }

    /// The placement a usage's `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` states.
    ///
    /// The transformation aligns an axis placement in the child's space with
    /// one in the parent's; which item is which follows from which side of
    /// the representation relationship is the child's own shape
    /// representation, not from argument order, because real files disagree
    /// about the order.
    fn usage_transform(&mut self, nauo: u64, child_sr: Option<u64>) -> Option<Transform> {
        if self.cdsr_of_nauo.is_none() {
            // One pass over the CDSRs, each resolved to the usage it
            // describes; ascending id order so a usage described twice keeps
            // the same one the old lowest-id-first scan chose.
            let mut cdsrs: Vec<u64> = self
                .exchange
                .data
                .iter()
                .filter(|(_, inst)| {
                    inst.part("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
                        .is_some()
                })
                .map(|(id, _)| *id)
                .collect();
            cdsrs.sort_unstable();
            let mut index: HashMap<u64, u64> = HashMap::new();
            for id in cdsrs {
                let Some(owner) = self
                    .args(id, "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
                    .ok()
                    .and_then(|args| args.get(1).and_then(Arg::reference))
                    .and_then(|pds| self.definition_of_shape(pds))
                else {
                    continue;
                };
                index.entry(owner).or_insert(id);
            }
            self.cdsr_of_nauo = Some(index);
        }
        let cdsr = *self.cdsr_of_nauo.as_ref()?.get(&nauo)?;
        let rr = self
            .args(cdsr, "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
            .ok()?
            .first()
            .and_then(Arg::reference)?;
        let (rep_1, rep_2) = {
            let args = self.args(rr, "REPRESENTATION_RELATIONSHIP").ok()?;
            (
                args.get(2).and_then(Arg::reference)?,
                args.get(3).and_then(Arg::reference)?,
            )
        };
        let idt = self
            .args(rr, "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            .ok()?
            .first()
            .and_then(Arg::reference)?;
        let (item_1, item_2) = {
            let args = self.args(idt, "ITEM_DEFINED_TRANSFORMATION").ok()?;
            (
                args.get(2).and_then(Arg::reference)?,
                args.get(3).and_then(Arg::reference)?,
            )
        };
        let frame_1 = self.frame(item_1).ok()?;
        let frame_2 = self.frame(item_2).ok()?;
        // item_1 pairs with rep_1. When rep_1 is the child's representation,
        // the child's frame_1 lands on the parent's frame_2.
        let child_first = match child_sr {
            Some(sr) => rep_1 == sr || rep_2 != sr,
            None => true,
        };
        Some(if child_first {
            Transform::from_frame(&frame_2) * Transform::to_frame(&frame_1)
        } else {
            Transform::from_frame(&frame_1) * Transform::to_frame(&frame_2)
        })
    }

    /// Products into the document: assemblies for the parents, parts for the
    /// shaped, instances for the usage edges.
    fn build_products(&mut self, document: &mut ogeom_doc::Document, entries: Vec<PdEntry>) {
        let mut ids: HashMap<u64, ogeom_doc::ProductId> = HashMap::new();
        for entry in &entries {
            let shape = match entry.shapes.len() {
                0 => None,
                1 => Some(entry.shapes[0].clone()),
                _ => match ogeom_algo::build::make_compound(document.model_mut(), &entry.shapes) {
                    Ok(built) => Some(built.shape),
                    Err(_) => Some(entry.shapes[0].clone()),
                },
            };
            if entry.children.is_empty() {
                if let Some(shape) = shape {
                    ids.insert(entry.pd, document.add_part(&entry.name, shape));
                }
                // A product with neither shape nor children holds nothing a
                // document can say; it is left out.
            } else {
                let assembly = document.add_assembly(&entry.name);
                ids.insert(entry.pd, assembly);
                // An assembly with its own geometry keeps it as a body part
                // placed at identity — rare, but files do it.
                if let Some(shape) = shape {
                    let body = document.add_part(format!("{}-body", entry.name), shape);
                    let _ = document.add_instance(assembly, body, Transform::IDENTITY, None);
                }
            }
        }
        for entry in &entries {
            let Some(&parent) = ids.get(&entry.pd) else {
                continue;
            };
            for (child, at, name) in &entry.children {
                let Some(&child_id) = ids.get(child) else {
                    self.report.warnings.push(format!(
                        "#{child}: an assembly child holds nothing readable"
                    ));
                    continue;
                };
                if let Err(e) = document.add_instance(parent, child_id, *at, name.clone()) {
                    self.report
                        .warnings
                        .push(format!("assembly edge #{} -> #{child}: {e}", entry.pd));
                }
            }
        }
    }

    /// Colours from styled items, keyed to the shapes they style.
    fn colours(&mut self, by_msb: &HashMap<u64, Shape>) -> Vec<(Shape, ogeom_doc::Colour)> {
        let mut styled: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| {
                inst.part("STYLED_ITEM").is_some() || inst.part("OVER_RIDING_STYLED_ITEM").is_some()
            })
            .map(|(id, _)| *id)
            .collect();
        // Sorted so an item styled twice resolves the same way every run;
        // overriding styles carry higher instance numbers in practice, and a
        // later same-shape entry wins the map insertion downstream.
        styled.sort_unstable();
        let mut out = Vec::new();
        for id in styled {
            let args = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                let Some(args) = instance
                    .part("STYLED_ITEM")
                    .or_else(|| instance.part("OVER_RIDING_STYLED_ITEM"))
                else {
                    continue;
                };
                args.to_vec()
            };
            let Some(item) = args.get(2).and_then(Arg::reference) else {
                continue;
            };
            let Some(shape) = by_msb.get(&item).or_else(|| self.faces.get(&item)) else {
                continue;
            };
            let styles: Vec<u64> = args
                .get(1)
                .and_then(Arg::list)
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            let shape = shape.clone();
            if let Some(colour) = styles.iter().find_map(|&style| self.colour_in(style, 0)) {
                out.push((shape, colour));
            }
        }
        out
    }

    /// The first `COLOUR_RGB` reachable from a presentation style, depth-bounded.
    ///
    /// The styled-item chain has five links and real files rearrange them, so
    /// the walk follows references rather than the textbook path.
    fn colour_in(&mut self, id: u64, depth: usize) -> Option<ogeom_doc::Colour> {
        if depth > 6 {
            return None;
        }
        let (keyword, args) = {
            let instance = self.instance(id).ok()?;
            let all: Vec<Arg> = instance
                .parts()
                .flat_map(|(_, args)| args.iter().cloned())
                .collect();
            (instance.keyword().to_owned(), all)
        };
        if keyword == "COLOUR_RGB" {
            let channel = |i: usize| args.get(i).and_then(Arg::number);
            return Some(ogeom_doc::Colour::rgb(
                channel(1)?,
                channel(2)?,
                channel(3)?,
            ));
        }
        let mut refs: Vec<u64> = Vec::new();
        collect_refs(&args, &mut refs);
        refs.into_iter()
            .find_map(|next| self.colour_in(next, depth + 1))
    }

    // --- semantic PMI --------------------------------------------------------

    /// The file's semantic PMI: dimensions, geometric tolerances, datums.
    ///
    /// Annotations that resist stay out with a warning; PMI never takes the
    /// geometry down.
    fn pmi_of(&mut self) -> ogeom_doc::Pmi {
        let mut pmi = ogeom_doc::Pmi::new();
        // Which STEP instance each annotation came from, so the presentation
        // pass can name the annotation a callout draws exactly rather than by
        // matching a string two annotations may share.
        let mut annotation_ids: HashMap<u64, ogeom_doc::Annotated> = HashMap::new();

        // Which topology each shape aspect describes: directly through
        // GEOMETRIC_ITEM_SPECIFIC_USAGE, and one relationship step outward,
        // because composite aspects hold their pieces through relationships.
        let mut aspect_items: HashMap<u64, Vec<ogeom_topo::TShapeId>> = HashMap::new();
        let mut gisus = self.ids_with("GEOMETRIC_ITEM_SPECIFIC_USAGE");
        gisus.sort_unstable();
        for id in gisus {
            let Ok(args) = self.args(id, "GEOMETRIC_ITEM_SPECIFIC_USAGE") else {
                continue;
            };
            let (Some(aspect), Some(item)) = (
                args.get(2).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) else {
                continue;
            };
            let node = self
                .faces
                .get(&item)
                .map(Shape::node)
                .or_else(|| self.edges.get(&item).map(|(shape, ..)| shape.node()));
            if let Some(node) = node {
                aspect_items.entry(aspect).or_default().push(node);
            }
        }
        let mut adjacency: HashMap<u64, Vec<u64>> = HashMap::new();
        for id in self.ids_with("SHAPE_ASPECT_RELATIONSHIP") {
            let Ok(args) = self.args(id, "SHAPE_ASPECT_RELATIONSHIP") else {
                continue;
            };
            if let (Some(a), Some(b)) = (
                args.get(2).and_then(Arg::reference),
                args.get(3).and_then(Arg::reference),
            ) {
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
        let items_for = |aspect: u64| -> Vec<ogeom_topo::TShapeId> {
            // Three relationship steps: a composite aspect holds components,
            // a derived aspect sits behind a composite, and a datum one link
            // behind its features — the deepest chain the corpus exhibits.
            let mut reach = vec![aspect];
            for _ in 0..3 {
                let mut next = reach.clone();
                for a in &reach {
                    next.extend(adjacency.get(a).into_iter().flatten().copied());
                }
                next.sort_unstable();
                next.dedup();
                reach = next;
            }
            let mut out: Vec<ogeom_topo::TShapeId> = Vec::new();
            for a in reach {
                out.extend(aspect_items.get(&a).into_iter().flatten().copied());
            }
            out.sort_unstable();
            out.dedup();
            out
        };

        // Dimensions: characteristic -> representation, values from the
        // measure items, bounds from any plus/minus tolerance over the same
        // characteristic.
        let mut plus_minus: HashMap<u64, (Option<f64>, Option<f64>)> = HashMap::new();
        for id in self.ids_with("PLUS_MINUS_TOLERANCE") {
            let Ok(args) = self.args(id, "PLUS_MINUS_TOLERANCE") else {
                continue;
            };
            let (Some(tv), Some(dim)) = (
                args.first().and_then(Arg::reference),
                args.get(1).and_then(Arg::reference),
            ) else {
                continue;
            };
            let Ok(tv_args) = self.args(tv, "TOLERANCE_VALUE") else {
                continue;
            };
            let lower = tv_args
                .first()
                .and_then(Arg::reference)
                .and_then(|r| self.measure_value(r))
                .map(|(v, _)| v);
            let upper = tv_args
                .get(1)
                .and_then(Arg::reference)
                .and_then(|r| self.measure_value(r))
                .map(|(v, _)| v);
            plus_minus.insert(dim, (lower, upper));
        }
        let mut dcrs = self.ids_with("DIMENSIONAL_CHARACTERISTIC_REPRESENTATION");
        dcrs.sort_unstable();
        for id in dcrs {
            let Ok(args) = self.args(id, "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION") else {
                continue;
            };
            let (Some(dim), Some(sdr)) = (
                args.first().and_then(Arg::reference),
                args.get(1).and_then(Arg::reference),
            ) else {
                continue;
            };
            let mut values = Vec::new();
            let mut kind = ogeom_doc::MeasureKind::Length;
            if let Ok(sdr_args) = self.args(sdr, "SHAPE_DIMENSION_REPRESENTATION") {
                for item in sdr_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                    if let Some(r) = item.reference()
                        && let Some((v, k)) = self.measure_value(r)
                    {
                        values.push(v);
                        kind = k;
                    }
                }
            }
            let (name, location, aspects) = self.dimension_shape(dim);
            let features: Vec<Vec<ogeom_topo::TShapeId>> =
                aspects.into_iter().map(&items_for).collect();
            let (minus, plus) = plus_minus.get(&dim).copied().unwrap_or((None, None));
            annotation_ids.insert(dim, ogeom_doc::Annotated::Dimension(pmi.dimensions.len()));
            pmi.dimensions.push(ogeom_doc::Dimension {
                name,
                values,
                kind,
                plus,
                minus,
                features,
                location,
            });
        }

        // Geometric tolerances: any instance whose subtype names one. The
        // complex form keeps its attributes on the GEOMETRIC_TOLERANCE part;
        // the simple form flattens them into the subtype's own list.
        let is_subtype = |k: &str| {
            k.ends_with("_TOLERANCE")
                && !matches!(
                    k,
                    "GEOMETRIC_TOLERANCE"
                        | "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE"
                        | "GEOMETRIC_TOLERANCE_WITH_DEFINED_UNIT"
                        | "GEOMETRIC_TOLERANCE_WITH_MODIFIERS"
                        | "GEOMETRIC_TOLERANCE_WITH_MAXIMUM_TOLERANCE"
                        | "PLUS_MINUS_TOLERANCE"
                        | "TOLERANCE_VALUE"
                )
        };
        let mut gts: Vec<u64> = self
            .exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.parts().any(|(k, _)| is_subtype(k)))
            .map(|(id, _)| *id)
            .collect();
        gts.sort_unstable();
        for id in gts {
            let (name, magnitude_ref, aspect, kind, modifiers, datum_refs) = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                let subtype = instance
                    .parts()
                    .map(|(k, _)| k.to_owned())
                    .find(|k| is_subtype(k));
                let Some(subtype) = subtype else {
                    continue;
                };
                let base = instance
                    .part("GEOMETRIC_TOLERANCE")
                    .or_else(|| instance.part(&subtype))
                    .unwrap_or(&[])
                    .to_vec();
                let name = match base.first() {
                    Some(Arg::Str(s)) => s.clone(),
                    _ => String::new(),
                };
                let magnitude_ref = base.get(2).and_then(Arg::reference);
                let aspect = base.get(3).and_then(Arg::reference);
                let kind = Some(subtype.trim_end_matches("_TOLERANCE").to_lowercase());
                // Modifiers ride on their own part in the complex form, as
                // a list of enumeration words.
                let modifiers: Vec<String> = instance
                    .part("GEOMETRIC_TOLERANCE_WITH_MODIFIERS")
                    .and_then(|args| args.first())
                    .and_then(Arg::list)
                    .map(|list| {
                        list.iter()
                            .filter_map(|arg| match arg {
                                Arg::Enum(word) => Some(word.to_lowercase()),
                                _ => None,
                            })
                            .collect()
                    })
                    .unwrap_or_default();
                // The complex form keeps the datum list on its own part;
                // the simple form appends it as the subtype's fifth argument.
                let datum_refs: Vec<u64> = instance
                    .part("GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE")
                    .and_then(|args| args.first())
                    .and_then(Arg::list)
                    .or_else(|| base.get(4).and_then(Arg::list))
                    .map(|list| list.iter().filter_map(Arg::reference).collect())
                    .unwrap_or_default();
                (name, magnitude_ref, aspect, kind, modifiers, datum_refs)
            };
            let Some(kind) = kind else {
                continue;
            };
            let magnitude = magnitude_ref
                .and_then(|r| self.measure_value(r))
                .map_or(0.0, |(v, _)| v);
            let datums: Vec<String> = datum_refs
                .iter()
                .filter_map(|&r| self.datum_letter(r, 0))
                .collect();
            let items = aspect.map(items_for).unwrap_or_default();
            annotation_ids.insert(id, ogeom_doc::Annotated::Tolerance(pmi.tolerances.len()));
            pmi.tolerances.push(ogeom_doc::GeometricTolerance {
                kind,
                name,
                magnitude,
                modifiers,
                datums,
                items,
            });
        }

        // Datums: the letters, with the features they mark reached through
        // the aspect graph.
        let mut datums = self.ids_with("DATUM");
        datums.sort_unstable();
        for id in datums {
            let letter = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                if instance.part("DATUM_FEATURE").is_some()
                    || instance.part("DATUM_REFERENCE").is_some()
                    || instance.part("DATUM_REFERENCE_COMPARTMENT").is_some()
                    || instance.part("DATUM_SYSTEM").is_some()
                {
                    continue;
                }
                let Some(args) = instance.part("DATUM") else {
                    continue;
                };
                match args.get(4) {
                    Some(Arg::Str(s)) if !s.is_empty() => s.clone(),
                    _ => continue,
                }
            };
            annotation_ids.insert(id, ogeom_doc::Annotated::Datum(pmi.datums.len()));
            pmi.datums.push(ogeom_doc::Datum {
                label: letter,
                items: items_for(id),
            });
        }

        // Datum targets: the pads a datum is actually established at. The
        // target's identifier is the letter's number — `A1` is target 1 of
        // datum A — and its placement and size come through the shape
        // representation the feature is associated with.
        let mut targets = self.ids_with("PLACED_DATUM_TARGET_FEATURE");
        targets.sort_unstable();
        for id in targets {
            if let Some(target) = self.datum_target(id, &items_for) {
                pmi.targets.push(target);
            }
        }

        // Presentation: what a viewer draws. A callout holds tessellated
        // annotation occurrences, each an indexed set of polylines over a
        // coordinates list; an annotation plane says which plane they are
        // drawn in and which callouts it holds; and a model item association
        // says which semantic annotation a callout is the picture of.
        let (callouts, callout_index) = self.callouts(&annotation_ids);
        pmi.callouts = callouts;
        self.callout_index = callout_index;
        pmi
    }

    /// Saved views: every *named* draughting model is one — the unnamed one
    /// is the annotation-plane container this writer emits itself. The
    /// camera item gives the frame; the callout items give the subset.
    fn views(&mut self) -> Vec<ogeom_doc::View> {
        let mut out = Vec::new();
        for id in self.ids_with("DRAUGHTING_MODEL") {
            let Ok(args) = self.args(id, "DRAUGHTING_MODEL") else {
                continue;
            };
            let name = match args.first() {
                Some(Arg::Str(text)) => text.clone(),
                _ => String::new(),
            };
            if name.is_empty() {
                continue;
            }
            let mut frame = None;
            let mut callouts = Vec::new();
            for item in args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                let Some(item) = item.reference() else {
                    continue;
                };
                if let Ok(cam) = self.args(item, "CAMERA_MODEL_D3") {
                    if let Some(placement) = cam.get(1).and_then(Arg::reference) {
                        frame = self.frame(placement).ok();
                    }
                    continue;
                }
                if let Some(&index) = self.callout_index.get(&item) {
                    callouts.push(index);
                }
            }
            let Some(frame) = frame else { continue };
            out.push(ogeom_doc::View {
                name,
                frame,
                clipping: None,
                callouts,
            });
        }
        out
    }

    /// One placed datum target: which datum, which number, where and how big.
    fn datum_target(
        &mut self,
        id: u64,
        items_for: &dyn Fn(u64) -> Vec<ogeom_topo::TShapeId>,
    ) -> Option<ogeom_doc::DatumTarget> {
        let (target_id, description) = {
            let instance = self.instance(id).ok()?;
            let args = instance.part("PLACED_DATUM_TARGET_FEATURE")?;
            let target = match args.get(4) {
                Some(Arg::Str(s)) if !s.is_empty() => s.clone(),
                _ => return None,
            };
            let description = match args.get(1) {
                Some(Arg::Str(s)) => s.to_ascii_lowercase(),
                _ => String::new(),
            };
            (target, description)
        };
        // `A1`: the letter is the datum, the digits its number.
        let split = target_id
            .find(|c: char| c.is_ascii_digit())
            .unwrap_or(target_id.len());
        let (letter, number) = target_id.split_at(split);
        let index = number.parse::<u32>().unwrap_or(1);
        let datum = if letter.is_empty() {
            self.datum_letter(id, 0).unwrap_or_default()
        } else {
            letter.to_owned()
        };

        // The target's placement and size live in a shape representation the
        // target's own property definition names. The lengths come in the
        // file's own unit, as every length does. Both hops go through the
        // indexes: the old form rescanned every entity per property *per
        // target*, the assembly quadratic one storey deeper.
        self.ensure_property_indexes();
        let mut frame = None;
        let mut lengths: Vec<f64> = Vec::new();
        let properties: Vec<u64> = self
            .properties_of_definition
            .as_ref()
            .and_then(|m| m.get(&id).cloned())
            .unwrap_or_default();
        for property in properties {
            let sdrs: Vec<u64> = self
                .sdrs_of_property
                .as_ref()
                .and_then(|m| m.get(&property).cloned())
                .unwrap_or_default();
            for sdr in sdrs {
                let Ok(args) = self.args(sdr, "SHAPE_DEFINITION_REPRESENTATION") else {
                    continue;
                };
                let Some(rep) = args.get(1).and_then(Arg::reference) else {
                    continue;
                };
                for item in self.representation_items(rep).unwrap_or_default() {
                    let keyword = self
                        .instance(item)
                        .map(|i| i.keyword().to_owned())
                        .unwrap_or_default();
                    match keyword.as_str() {
                        "AXIS2_PLACEMENT_3D" => frame = self.frame(item).ok(),
                        "LENGTH_MEASURE_WITH_UNIT" | "MEASURE_REPRESENTATION_ITEM" => {
                            if let Some((value, _)) = self.measure_value(item) {
                                lengths.push(value);
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        // What the target *is* follows from the description the file gives
        // and how many sizes it carries: an area needs two, a circle one, a
        // line one, a point none.
        let kind = if description.contains("circle") || description.contains("circular") {
            ogeom_doc::DatumTargetKind::Circle {
                diameter: lengths.first().copied().unwrap_or(0.0),
            }
        } else if description.contains("rectangle") || lengths.len() >= 2 {
            ogeom_doc::DatumTargetKind::Rectangle {
                length: lengths.first().copied().unwrap_or(0.0),
                width: lengths.get(1).copied().unwrap_or(0.0),
            }
        } else if description.contains("line") || lengths.len() == 1 {
            ogeom_doc::DatumTargetKind::Line {
                length: lengths.first().copied().unwrap_or(0.0),
            }
        } else {
            ogeom_doc::DatumTargetKind::Point
        };
        Some(ogeom_doc::DatumTarget {
            datum,
            index,
            kind,
            at: frame.map_or(Point::ORIGIN, |f: Frame| f.origin()),
            frame,
            items: items_for(id),
        })
    }

    /// The drawn annotations: callouts, their polylines, their planes, and
    /// which semantic annotation each one draws.
    fn callouts(
        &mut self,
        annotation_ids: &HashMap<u64, ogeom_doc::Annotated>,
    ) -> (Vec<ogeom_doc::Callout>, HashMap<u64, usize>) {
        // Which callout each annotation plane holds, and the plane's frame.
        let mut plane_of: HashMap<u64, Frame> = HashMap::new();
        for id in self.ids_with("ANNOTATION_PLANE") {
            let Ok(args) = self.args(id, "ANNOTATION_PLANE") else {
                continue;
            };
            let frame = args
                .get(2)
                .and_then(Arg::reference)
                .and_then(|plane| {
                    let inner = self.args(plane, "PLANE").ok()?;
                    inner.get(1).and_then(Arg::reference)
                })
                .and_then(|placement| self.frame(placement).ok());
            let Some(frame) = frame else { continue };
            for element in args.get(3).and_then(Arg::list).unwrap_or(&[]) {
                if let Some(callout) = element.reference() {
                    plane_of.insert(callout, frame);
                }
            }
        }

        // Which semantic annotation each callout draws. The association names
        // the annotation's own STEP id, so the link is made by matching that
        // against the ids the semantic pass already resolved.
        // A callout may be associated more than once — once with the shape
        // aspect the annotation is *about*, once with the annotation itself —
        // so every association is kept and the first that resolves to a
        // semantic annotation is the answer.
        let mut draws: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut associations = self.ids_with("DRAUGHTING_MODEL_ITEM_ASSOCIATION");
        associations.sort_unstable();
        for id in associations {
            let Ok(args) = self.args(id, "DRAUGHTING_MODEL_ITEM_ASSOCIATION") else {
                continue;
            };
            if let (Some(definition), Some(item)) = (
                args.get(2).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) {
                draws.entry(item).or_default().push(definition);
            }
        }

        let mut out = Vec::new();
        let mut index_of: HashMap<u64, usize> = HashMap::new();
        let mut callouts = self.ids_with("DRAUGHTING_CALLOUT");
        callouts.sort_unstable();
        for id in callouts {
            let Ok(args) = self.args(id, "DRAUGHTING_CALLOUT") else {
                continue;
            };
            let name = match args.first() {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let mut polylines = Vec::new();
            for element in args.get(1).and_then(Arg::list).unwrap_or(&[]).to_vec() {
                let Some(occurrence) = element.reference() else {
                    continue;
                };
                polylines.extend(self.annotation_polylines(occurrence));
            }
            if polylines.is_empty() && !plane_of.contains_key(&id) {
                continue;
            }
            index_of.insert(id, out.len());
            out.push(ogeom_doc::Callout {
                name,
                plane: plane_of.get(&id).copied(),
                polylines,
                annotates: draws
                    .get(&id)
                    .cloned()
                    .unwrap_or_default()
                    .into_iter()
                    .find_map(|d| annotation_ids.get(&d).copied()),
            });
        }
        (out, index_of)
    }

    /// The polylines one annotation occurrence draws.
    ///
    /// A tessellated occurrence names a curve set, which names a coordinates
    /// list and gives each polyline as *one-based* indices into it. That is
    /// the whole of the drawn geometry, and it is read as it is written
    /// rather than resampled.
    fn annotation_polylines(&mut self, occurrence: u64) -> Vec<Vec<Point>> {
        let item = {
            let Ok(instance) = self.instance(occurrence) else {
                return Vec::new();
            };
            let args = instance
                .part("TESSELLATED_ANNOTATION_OCCURRENCE")
                .or_else(|| instance.part("ANNOTATION_OCCURRENCE"))
                .or_else(|| instance.part("STYLED_ITEM"));
            match args.and_then(|a| a.get(2).and_then(Arg::reference)) {
                Some(item) => item,
                None => return Vec::new(),
            }
        };
        self.tessellated_polylines(item, 0)
    }

    /// The polylines a tessellated item holds, however deep it nests them.
    ///
    /// The drawn geometry of one annotation is a *set*: a frame's box, its
    /// leader, its text strokes, each its own curve set over its own
    /// coordinates list, gathered under one item — which may itself be
    /// repositioned by a placement, and that placement is applied here rather
    /// than left for a consumer to discover.
    fn tessellated_polylines(&mut self, item: u64, depth: usize) -> Vec<Vec<Point>> {
        if depth > 4 {
            return Vec::new();
        }
        let (children, curve_set, placement) = {
            let Ok(instance) = self.instance(item) else {
                return Vec::new();
            };
            let children: Vec<u64> = instance
                .part("TESSELLATED_GEOMETRIC_SET")
                .and_then(|args| args.first().and_then(Arg::list))
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            let curve_set = instance.part("TESSELLATED_CURVE_SET").map(<[Arg]>::to_vec);
            let placement = instance
                .part("REPOSITIONED_TESSELLATED_ITEM")
                .and_then(|args| args.first().and_then(Arg::reference));
            (children, curve_set, placement)
        };

        let mut out = Vec::new();
        for child in children {
            out.extend(self.tessellated_polylines(child, depth + 1));
        }
        if let Some(args) = curve_set
            && let Some(list) = args.get(1).and_then(Arg::reference)
        {
            let points = self.coordinates_list(list);
            for line in args.get(2).and_then(Arg::list).unwrap_or(&[]) {
                let Some(indices) = line.list() else { continue };
                let mut polyline = Vec::with_capacity(indices.len());
                for index in indices {
                    // One-based, as the format states them.
                    let Some(k) = index.number() else { continue };
                    #[expect(
                        clippy::cast_possible_truncation,
                        clippy::cast_sign_loss,
                        reason = "an index into a list the file itself sized"
                    )]
                    let k = k as usize;
                    if k >= 1 && k <= points.len() {
                        polyline.push(points[k - 1]);
                    }
                }
                if polyline.len() >= 2 {
                    out.push(polyline);
                }
            }
        }
        if let Some(placement) = placement
            && let Ok(frame) = self.frame(placement)
        {
            for polyline in &mut out {
                for p in polyline.iter_mut() {
                    *p = frame.to_world(*p);
                }
            }
        }
        out
    }

    /// A `COORDINATES_LIST`'s points, in the document's own length unit.
    fn coordinates_list(&mut self, id: u64) -> Vec<Point> {
        let Ok(args) = self.args(id, "COORDINATES_LIST") else {
            return Vec::new();
        };
        let scale = self.report.scale_mm;
        args.get(2)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(|entry| {
                let coords = entry.list()?;
                let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(Point::new(value(0), value(1), value(2)))
            })
            .collect()
    }

    /// Every instance id carrying a part with this keyword.
    fn ids_with(&self, keyword: &str) -> Vec<u64> {
        self.exchange
            .data
            .iter()
            .filter(|(_, inst)| inst.part(keyword).is_some())
            .map(|(id, _)| *id)
            .collect()
    }

    /// One pass over the exchange builds both property-chain indexes:
    /// definition → its `PROPERTY_DEFINITION`s (by the definition argument),
    /// property → its `SHAPE_DEFINITION_REPRESENTATION`s (by the definition
    /// they represent). Ascending ids inside each list, so a lookup visits
    /// candidates in the same order the old full scan would have.
    fn ensure_property_indexes(&mut self) {
        if self.properties_of_definition.is_some() {
            return;
        }
        let mut properties: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut sdrs: HashMap<u64, Vec<u64>> = HashMap::new();
        let exchange = self.exchange;
        for (id, instance) in &exchange.data {
            if let Some(args) = instance.part("PROPERTY_DEFINITION") {
                // Read for the index is read: the skipped table should not
                // claim the reader never looked.
                self.visited.insert(*id);
                if let Some(definition) = args.get(2).and_then(Arg::reference) {
                    properties.entry(definition).or_default().push(*id);
                }
            }
            if let Some(args) = instance.part("SHAPE_DEFINITION_REPRESENTATION") {
                self.visited.insert(*id);
                if let Some(property) = args.first().and_then(Arg::reference) {
                    sdrs.entry(property).or_default().push(*id);
                }
            }
        }
        for list in properties.values_mut().chain(sdrs.values_mut()) {
            list.sort_unstable();
        }
        self.properties_of_definition = Some(properties);
        self.sdrs_of_property = Some(sdrs);
    }

    /// A measure item's value and kind, scaled into the document's units.
    fn measure_value(&mut self, id: u64) -> Option<(f64, ogeom_doc::MeasureKind)> {
        let args = {
            let instance = self.instance(id).ok()?;
            instance.part("MEASURE_WITH_UNIT")?.to_vec()
        };
        match args.first() {
            Some(Arg::Typed(kind, inner)) => {
                let value = inner.first().and_then(Arg::number)?;
                if kind.contains("ANGLE") {
                    Some((value * self.angle_scale, ogeom_doc::MeasureKind::Angle))
                } else {
                    Some((value * self.report.scale_mm, ogeom_doc::MeasureKind::Length))
                }
            }
            _ => None,
        }
    }

    /// A dimensional characteristic's name, kind and aspects.
    ///
    /// Sizes apply to one feature; locations — linear or angular — run
    /// between two. `ANGULAR_SIZE` and `ANGULAR_LOCATION` are the same
    /// shapes with an extra angle-selection argument at the end.
    fn dimension_shape(&mut self, dim: u64) -> (String, bool, Vec<u64>) {
        let size = {
            let Ok(instance) = self.instance(dim) else {
                return (String::new(), false, Vec::new());
            };
            instance
                .part("DIMENSIONAL_SIZE")
                .or_else(|| instance.part("ANGULAR_SIZE"))
                .map(<[Arg]>::to_vec)
        };
        if let Some(args) = size {
            let name = match args.get(1) {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let aspects = args.first().and_then(Arg::reference).into_iter().collect();
            return (name, false, aspects);
        }
        let location = {
            let Ok(instance) = self.instance(dim) else {
                return (String::new(), false, Vec::new());
            };
            instance
                .part("DIMENSIONAL_LOCATION")
                .or_else(|| instance.part("ANGULAR_LOCATION"))
                .map(<[Arg]>::to_vec)
        };
        if let Some(args) = location {
            let name = match args.first() {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let aspects = [args.get(2), args.get(3)]
                .into_iter()
                .flatten()
                .filter_map(Arg::reference)
                .collect();
            return (name, true, aspects);
        }
        (String::new(), false, Vec::new())
    }

    /// The datum letter reachable from a datum reference, depth-bounded: the
    /// reference chain runs through compartments and systems, and files
    /// arrange it differently.
    fn datum_letter(&mut self, id: u64, depth: usize) -> Option<String> {
        if depth > 4 {
            return None;
        }
        let (letter, members, refs) = {
            let instance = self.instance(id).ok()?;
            let letter = instance.part("DATUM").and_then(|args| match args.get(4) {
                Some(Arg::Str(s)) if !s.is_empty() => Some(s.clone()),
                _ => None,
            });
            // A compartment or common datum carrying a list of constituent
            // references is a composite: every constituent contributes a
            // letter, and the letters act as one datum.
            let members: Vec<u64> = ["DATUM_REFERENCE_COMPARTMENT", "COMMON_DATUM"]
                .iter()
                .find_map(|k| instance.part(k))
                .into_iter()
                .flatten()
                .find_map(Arg::list)
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            let mut refs = Vec::new();
            for (_, args) in instance.parts() {
                collect_refs(args, &mut refs);
            }
            (letter, members, refs)
        };
        if let Some(letter) = letter {
            return Some(letter);
        }
        if !members.is_empty() {
            let letters: Vec<String> = members
                .into_iter()
                .filter_map(|r| self.datum_letter(r, depth + 1))
                .collect();
            if !letters.is_empty() {
                return Some(letters.join("-"));
            }
        }
        refs.into_iter()
            .find_map(|r| self.datum_letter(r, depth + 1))
    }
}

/// A product definition gathered from the file: its shapes and its usage
/// edges, before anything is committed to the document.
struct PdEntry {
    pd: u64,
    name: String,
    shapes: Vec<Shape>,
    children: Vec<(u64, Transform, Option<String>)>,
}

/// Every reference in an argument tree, in order.
fn collect_refs(args: &[Arg], out: &mut Vec<u64>) {
    for arg in args {
        match arg {
            Arg::Ref(id) => out.push(*id),
            Arg::List(inner) | Arg::Typed(_, inner) => collect_refs(inner, out),
            _ => {}
        }
    }
}

/// The chart coordinates of a point on an analytic surface, by closed-form
/// inversion — `None` for surfaces that need iterative projection.
/// For a two-wire periodic face: each wire's single closed edge with its
/// vertex, empty when the shape is anything else.
fn closed_ring_edges(model: &Model, wires: &[Shape]) -> OgeomResult<Vec<(Shape, Shape)>> {
    let mut out = Vec::new();
    for wire in wires {
        let edges = ogeom_topo::explore(
            model,
            wire,
            ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Edge),
        )?;
        if edges.len() != 1 {
            return Ok(Vec::new());
        }
        let edge = edges[0].clone();
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &edge)? else {
            return Ok(Vec::new());
        };
        if !a.is_same(&b) {
            return Ok(Vec::new());
        }
        out.push((edge, a));
    }
    Ok(out)
}

/// Unused-import guard for kinds only touched via traits.
#[allow(dead_code)]
fn _keep(p: PlanarCurve) -> PlanarCurve {
    p
}
