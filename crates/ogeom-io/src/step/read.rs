//! From parsed exchange structure to a living model.
//!
//! The reader walks every `MANIFOLD_SOLID_BREP` (`BREP_WITH_VOIDS` is one,
//! by subtype) and every `SHELL_BASED_SURFACE_MODEL`, and rebuilds them
//! bottom-up:
//! points, placements, curves and surfaces into geometry; vertices, edges,
//! loops, faces and shells into topology, shared exactly as the file shares
//! them: a vertex referenced by eight edges is one vertex here too, which is
//! what lets a closed shell close. Edge ranges are re-derived on this
//! kernel's own parameterizations from the vertex geometry, because STEP's
//! parameterizations are its own business and carrying them over blind is
//! how off-by-a-period bugs are born. Geometric sets in geometrically
//! bounded representations are construction geometry: named items with no
//! topology, kept on the document under their product.
//!
//! What the reader does not understand it *counts*: every instance never
//! visited lands in the report's skipped table by keyword, and every
//! compromise (a face without a pcurve, a shell that does not close) is a
//! warning with the instance number in it. An import that succeeded with
//! three warnings is a different thing from one that succeeded, and the
//! report is what keeps the difference visible.

use super::parse::{Arg, Instance, Instances};
use ogeom_algo::{
    make_edge_between, make_face_on, make_shell, make_solid, make_vertex, make_wire,
    project_on_curve,
};
use ogeom_core::{FastMap, FastSet};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_doc::ConstructionGeometry;
use ogeom_geom::Curve2d as _;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{
    BSplineCurve, CircleCurve, ConeSurface, Curve, CylinderSurface, EllipseCurve, ExtrusionSurface,
    HyperbolaCurve, LineCurve, ParabolaCurve, PlanarCurve, PlaneSurface, SphereSurface,
    SurfaceGeometry, TorusSurface,
};
use ogeom_math::{
    Axis, Blend as _, Circle, Cone, Cylinder, Direction, Ellipse, Frame, Hyperbola, KnotVector,
    Matrix3, Parabola, Plane, Point, Sphere, Torus, Transform, Vector, Weighted,
};
use ogeom_topo::{Location, Model, Shape};
use std::collections::BTreeMap;

mod geometry;
mod pcurves;
mod pmi;
mod products;
mod splines;
mod topology;

use pcurves::{ChartMap, derive_pcurve, in_kernel_units, through_map, walked_back};
use splines::{Patch, expanded_knots, file_degree, implied_knots, join_patches, join_splines};

/// One use a bound makes of an edge: the edge as the loop walks it, the edge
/// as it was built, the file's id for it, and the curve and range it carries.
type BoundUse = (Shape, Shape, u64, Curve, (f64, f64));

/// How a point of space reads on a surface's unbounded directions, `None`
/// where the direction has none to read.
type ChartReading = dyn Fn(Point) -> (Option<f64>, Option<f64>);

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
    /// measures one, and one exemplar entity id. A file with hundreds of
    /// warnings summarises to a handful of lines a consumer can show;
    /// `warnings` keeps the full prose. Sorted by count, largest first.
    pub summary: Vec<WarningSummary>,
    /// Faces that read without a complete trim: an edge's boundary sat too
    /// far from the surface for any honest pcurve (beyond the one-millimetre
    /// healing cap), so the face will refuse to triangulate. Deduplicated,
    /// in file order: the structured form of the warnings that name them,
    /// carrying the face itself so the instructed follow-up needs no search.
    /// `check` reports the same faces as broken from the model side. A
    /// refused id whose face never finished building has nothing to act on
    /// and stays in the warnings alone.
    pub untrimmed_faces: Vec<UntrimmedFace>,
    /// Edges on faces whose pcurve the file states and the reader kept: it
    /// follows the edge's curve, once lifted, as closely as a derived one.
    /// Each face's use of an edge counts once.
    pub pcurves_kept: usize,
    /// Edges on faces whose pcurve the file states but the reader derived
    /// instead, the file's standing farther from the edge's curve.
    pub pcurves_derived: usize,
}

/// One kind of imperfect import, counted rather than repeated.
#[derive(Debug, Clone)]
pub struct WarningSummary {
    /// The kind, stable across runs: `"vertex-miss"`, `"boundary-slop"`,
    /// `"fit-short"`, `"untrimmed"`, `"pole-crossed"`.
    pub kind: &'static str,
    /// How many times it happened.
    pub count: usize,
    /// The worst measured value among them: a distance, for every kind
    /// that measures one; zero where none applies.
    pub worst: f64,
    /// One entity id to look at first.
    pub exemplar: u64,
}

/// A face the reader could not trim, with the shape to act on.
#[derive(Debug, Clone)]
pub struct UntrimmedFace {
    /// The file's id for the face: the same id the warnings name.
    pub entity: u64,
    /// The face as built; hand it to `ogeom_heal::fix_face_pcurves` with
    /// the cap the situation deserves.
    pub face: Shape,
}

/// A read exchange file: the model, the bodies found, and the report.
#[derive(Debug)]
pub struct StepImport {
    /// The document everything was built into: the model, plus the file's
    /// product structure, names, colours and construction geometry.
    pub document: ogeom_doc::Document,
    /// One shape per `MANIFOLD_SOLID_BREP` (`BREP_WITH_VOIDS` included,
    /// its cavities among its shells), in file order.
    pub solids: Vec<Shape>,
    /// One shape per `SHELL_BASED_SURFACE_MODEL`, in file order: its shell,
    /// or a compound of its shells when it names several.
    ///
    /// A surface body is what a modeller exports for a part built from faces
    /// rather than from a solid. It stays a shell here even when it happens
    /// to close: the file did not call it a solid, and the document carries
    /// it under its product exactly as the file placed it.
    pub shells: Vec<Shape>,
    /// Every face read, with its `ADVANCED_FACE` or `FACE_SURFACE` id, in id
    /// order: the face as the solids and shells above hold it, so that
    /// [`Shape::is_partner`] finds it in them. The id names the face across
    /// reads of the same file, as no handle does.
    pub faces: Vec<(u64, Shape)>,
    /// What happened along the way.
    pub report: StepReport,
}

/// Read a STEP exchange file's B-rep content.
///
/// The geometric sets and axes of its geometrically bounded wireframe and
/// surface representations read as construction geometry, named, under the
/// product each representation belongs to.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the file does
/// not parse, contains no solid, sheet or geometric set, or a solid's
/// structure is broken in a way topology cannot represent. Faces the reader cannot complete become
/// warnings, not errors; the report says exactly what was compromised.
pub fn read_step(text: &str, tol: Tolerances) -> OgeomResult<StepImport> {
    let (_, instances) = super::parse::parse_instances(text)?;
    let mut reader = Reader {
        instances: &instances,
        model: Model::new(),
        report: StepReport {
            scale_mm: 1.0,
            ..StepReport::default()
        },
        angle_scale: 1.0,
        visited: Visited::for_instances(&instances),
        vertices: FastMap::default(),
        edges: FastMap::default(),
        faces: FastMap::default(),
        callout_index: FastMap::default(),
        pcurves: FastMap::default(),
        pieces: FastMap::default(),
        band_faces: FastSet::default(),
        unnamed_colours: FastSet::default(),
        pole_vertices: Vec::new(),
        piece_key: u64::MAX,
        untrimmed_ids: Vec::new(),
        tallies: FastMap::default(),
        cdsr_of_nauo: None,
        properties_of_definition: None,
        sdrs_of_property: None,
        depth: 0,
        tol,
    };
    reader.report.scale_mm = reader.unit_scale();
    reader.angle_scale = reader.angle_unit_scale();

    let mut solids = Vec::new();
    let mut shells = Vec::new();
    let mut by_item: FastMap<u64, Shape> = FastMap::default();
    // Solids and surface models, one walk in file order; a body is a body
    // to the progress bar whichever kind it is.
    let mut ids: Vec<(u64, bool)> = instances
        .iter()
        .filter_map(|(id, inst)| {
            if inst.part("MANIFOLD_SOLID_BREP").is_some() || inst.part("BREP_WITH_VOIDS").is_some()
            {
                Some((id, true))
            } else if inst.part("SHELL_BASED_SURFACE_MODEL").is_some() {
                Some((id, false))
            } else {
                None
            }
        })
        .collect();
    ids.sort_unstable();
    let total = ids.len() as u64;
    for (done, (id, is_solid)) in ids.into_iter().enumerate() {
        ogeom_core::progress::checkpoint()?;
        ogeom_core::progress::stage_at("step: solid", done as u64 + 1, total);
        if is_solid {
            let solid = reader.solid(id)?;
            by_item.insert(id, solid.clone());
            solids.push(solid);
        } else {
            let shell = reader.surface_model(id)?;
            by_item.insert(id, shell.clone());
            shells.push(shell);
        }
    }
    // Every widening the reader made (an edge's tolerance recording how far
    // its pcurves sit from its curve) reaches the vertices that bound it,
    // and a face's reaches its edges: the containment rule the checker holds
    // a solid to, kept by the reader that built it.
    for body in solids.iter().chain(shells.iter()) {
        ogeom_algo::restore_containment(&mut reader.model, body)?;
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
    // The noted refusals resolve to the faces themselves, once they
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
    crate::bodies::rechart_beside_poles(
        &mut reader.model,
        solids.iter().chain(shells.iter()),
        &mut reader.report.warnings,
        tol,
    )?;
    let construction = reader.construction_sets()?;
    if solids.is_empty() && shells.is_empty() && construction.is_empty() {
        ogeom_bail!(
            Construction,
            "the exchange file contains no MANIFOLD_SOLID_BREP, \
             BREP_WITH_VOIDS, SHELL_BASED_SURFACE_MODEL or geometric set to read"
        );
    }
    // Each face read, by the file's id for it, as the bodies hold it: a face
    // the reader turned or rebuilt after first building it is the one its
    // map holds last, and the body's own occurrence is handed out.
    let mut held: FastMap<ogeom_topo::TShapeId, Shape> = FastMap::default();
    for body in solids.iter().chain(shells.iter()) {
        for face in ogeom_topo::explore_unique(&reader.model, body, ogeom_topo::ShapeType::Face)? {
            held.entry(face.node()).or_insert(face);
        }
    }
    let mut faces: Vec<(u64, Shape)> = reader
        .faces
        .iter()
        .filter_map(|(&id, face)| Some((id, held.get(&face.node())?.clone())))
        .collect();
    faces.sort_unstable_by_key(|&(id, _)| id);
    let document = reader.document(&by_item, &solids, &shells, construction)?;

    // Everything never visited, counted by its leading keyword.
    for (id, instance) in instances.iter() {
        if !reader.visited.seen(id) {
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
        shells,
        faces,
        report: reader.report,
    })
}

/// Named construction geometry by the representation holding it, in file
/// order.
type ConstructionSets = Vec<(u64, Vec<(String, ConstructionGeometry)>)>;

/// An edge as built: the shape, its curve and range, and whether STEP's edge
/// direction runs against the curve.
type BuiltEdge = (Shape, Curve, (f64, f64), bool);

/// A pcurve derived ahead of the face that needs it.
///
/// Deriving one is a pure function of a curve, its range and a surface (no
/// model, no order), and it is most of the time spent building solids. So
/// it is done for a whole solid at once, off the walk that attaches it.
/// An `off` within the confusion distance is only known to be within it:
/// the edge is widened only past it.
enum PreparedPcurve {
    /// The file stated it, and it stands `off` from the curve once lifted,
    /// within what the edge's own tolerance allows.
    Stated { curve: PlanarCurve, off: f64 },
    /// The projection had a closed form, standing `off` from the curve
    /// once lifted.
    Exact { curve: PlanarCurve, off: f64 },
    /// It did not, and this is the fit, with what the fit cost and how far
    /// it stands from the curve once lifted.
    Fitted {
        curve: PlanarCurve,
        error: f64,
        met: bool,
        worst_off: f64,
        off: f64,
        warning: Option<String>,
    },
    /// Neither worked, for this reason; the face gets the warning.
    Refused(String),
}

struct Reader<'a> {
    instances: &'a Instances,
    model: Model,
    report: StepReport,
    /// Radians per file angle unit; degrees are common.
    angle_scale: f64,
    /// Which instances the reader has looked at.
    visited: Visited,
    vertices: FastMap<u64, Shape>,
    edges: FastMap<u64, BuiltEdge>,
    faces: FastMap<u64, Shape>,
    /// STEP id → index into the callouts just built, for the views pass.
    callout_index: FastMap<u64, usize>,
    /// `(face, edge)` → the pcurve already derived for it, from the parallel
    /// pass at the head of each solid.
    /// The flag says whether the file stated a pcurve for it.
    pcurves: FastMap<(u64, u64), (PreparedPcurve, bool)>,
    /// Edge → the pieces it was cut into where its curve runs through a
    /// pole of a face it bounds, each under a key of its own that stands in
    /// for the edge's id wherever a face's bound names it.
    pieces: FastMap<u64, Vec<(u64, BuiltEdge)>>,
    /// Faces bounded by two whole circles of their surface whose rims share
    /// a vertex column only because one rim was split on the other's: each
    /// builds as a band of rings through its common column.
    band_faces: FastSet<u64>,
    /// The `DRAUGHTING_PRE_DEFINED_COLOUR`s whose name is none of the eight
    /// the standard defines, warned of once each.
    unnamed_colours: FastSet<u64>,
    /// The vertices made at poles, shared by every piece that meets one.
    pole_vertices: Vec<Shape>,
    /// The last key handed to a piece; keys count down from the top of the
    /// range, clear of every file id.
    piece_key: u64,
    /// Faces noted untrimmed, by file id; resolved to shapes once the read
    /// is far enough along for the shapes to exist.
    untrimmed_ids: Vec<u64>,
    /// kind → (count, worst, exemplar), folded into the report's summary.
    tallies: FastMap<&'static str, (usize, f64, u64)>,
    /// Usage → its `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`, built once,
    /// so an assembly edge is a lookup rather than a scan of the whole
    /// exchange, which is O(usages × entities).
    cdsr_of_nauo: Option<FastMap<u64, u64>>,
    /// Definition → its `PROPERTY_DEFINITION`s, and property → its
    /// `SHAPE_DEFINITION_REPRESENTATION`s, built together once, so a datum
    /// target's lookup is not a scan of the whole exchange per property
    /// *per target*. Each list ascends by id, so the first entry that
    /// answers is the lowest id, the one a scan in file order reaches first.
    properties_of_definition: Option<FastMap<u64, Vec<u64>>>,
    sdrs_of_property: Option<FastMap<u64, Vec<u64>>>,
    /// How many geometry or shell builders are open on the stack. Entities
    /// refer to entities of their own kind (an offset surface to its basis,
    /// an oriented shell to its shell), so a file that refers in a circle
    /// would otherwise recurse until the stack runs out.
    depth: u32,
    tol: Tolerances,
}

/// How deep entities of one kind may refer to one another. Real files nest
/// a handful of levels (a trimmed curve of an offset of a replica); a chain
/// this long is a cycle or a hostile file.
const MOST_NESTING: u32 = 64;

impl<'a> Reader<'a> {
    /// Run a builder one level deeper, refusing past [`MOST_NESTING`].
    fn nested<T>(&mut self, build: impl FnOnce(&mut Self) -> OgeomResult<T>) -> OgeomResult<T> {
        if self.depth >= MOST_NESTING {
            ogeom_bail!(
                Construction,
                "entities refer to one another more than {MOST_NESTING} deep; \
                 the file refers in a circle"
            );
        }
        self.depth += 1;
        let built = build(self);
        self.depth -= 1;
        built
    }

    /// The instance behind an id. Borrowed from the exchange, not from the
    /// reader, so its arguments can be read while the reader builds.
    fn instance(&mut self, id: u64) -> OgeomResult<&'a Instance> {
        self.visited.mark(id);
        let instances = self.instances;
        instances.get(id).ok_or_else(|| {
            ogeom_core::ogeom_err!(
                Construction,
                "the file references #{id}, which does not exist"
            )
        })
    }

    /// A face's arguments (name, bounds, surface, sense), whether the file
    /// wrote it as the `ADVANCED_FACE` every modern writer uses or as the
    /// plain `FACE_SURFACE` it specialises, which carries the same four.
    fn face_args(&mut self, id: u64) -> OgeomResult<&'a [Arg]> {
        let instance = self.instance(id)?;
        if let Some(args) = instance
            .part("ADVANCED_FACE")
            .or_else(|| instance.part("FACE_SURFACE"))
        {
            return Ok(args);
        }
        ogeom_bail!(
            Construction,
            "#{id} is {}, where a face (ADVANCED_FACE or FACE_SURFACE) was needed",
            instance.keyword()
        );
    }

    fn args(&mut self, id: u64, keyword: &str) -> OgeomResult<&'a [Arg]> {
        let instance = self.instance(id)?;
        let Some(args) = instance.part(keyword) else {
            ogeom_bail!(
                Construction,
                "#{id} is {}, where {keyword} was needed",
                instance.keyword()
            );
        };
        Ok(args)
    }

    // --- units ---------------------------------------------------------------

    /// The unit instances the representation context actually assigns.
    ///
    /// A file may carry both a radian and a degree (definition and
    /// conversion), and which applies is not a matter of existence but of
    /// assignment: `GLOBAL_UNIT_ASSIGNED_CONTEXT` lists the ones in force.
    fn assigned_units(&self) -> Vec<u64> {
        // A file may carry several unit contexts (inches for the model,
        // millimetres for its annotation sheet), and the shapes' own
        // representations name the one their coordinates mean. That context
        // goes first; the rest follow in entity order, so the answer never
        // depends on how a map happens to iterate.
        let mut cited: Vec<u64> = self
            .instances
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
            .instances
            .iter()
            .filter(|(_, inst)| inst.part("GLOBAL_UNIT_ASSIGNED_CONTEXT").is_some())
            .map(|(id, _)| id)
            .collect();
        contexts.sort_unstable();
        contexts.sort_by_cached_key(|id| !cited.contains(id));

        let mut out = Vec::new();
        for id in contexts {
            if let Some(instance) = self.instances.get(id)
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
            let Some(instance) = self.instances.get(id) else {
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
                && let Some(inner) = self.instances.get(measure)
            {
                let factor = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(|a| match a {
                        Arg::Typed(typed) if typed.0 == "LENGTH_MEASURE" => {
                            typed.1.first().and_then(Arg::number)
                        }
                        _ => None,
                    });
                let base = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(Arg::reference)
                    .and_then(|b| self.instances.get(b))
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

    /// Radians per file angle unit: the same dance as length, for the files
    /// that measure their cones in degrees.
    fn angle_unit_scale(&mut self) -> f64 {
        let assigned = self.assigned_units();
        for id in assigned {
            let Some(instance) = self.instances.get(id) else {
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
                && let Some(inner) = self.instances.get(measure)
            {
                let factor = inner
                    .parts
                    .iter()
                    .flat_map(|(_, a)| a.iter())
                    .find_map(|a| match a {
                        Arg::Typed(typed) if typed.0 == "PLANE_ANGLE_MEASURE" => {
                            typed.1.first().and_then(Arg::number)
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
}

/// Every reference in an argument tree, in order.
fn collect_refs(args: &[Arg], out: &mut Vec<u64>) {
    for arg in args {
        match arg {
            Arg::Ref(id) => out.push(*id),
            Arg::List(inner) => collect_refs(inner, out),
            Arg::Typed(typed) => collect_refs(&typed.1, out),
            _ => {}
        }
    }
}

/// Whether two points stand on one column about `circle`'s axis: the same
/// angle about it, to within the confusion distance along the circle
/// through `q`.
fn same_column(circle: &ogeom_math::Circle, p: Point, q: Point, tol: Tolerances) -> bool {
    let (lp, lq) = (circle.frame().to_local(p), circle.frame().to_local(q));
    let apart = (lp.y.atan2(lp.x) - lq.y.atan2(lq.x)).rem_euclid(core::f64::consts::TAU);
    let turn = apart.min(core::f64::consts::TAU - apart);
    turn * lq.x.hypot(lq.y) <= tol.confusion()
}

/// For a two-wire periodic face: each wire's single closed edge with its
/// vertex, empty when the shape is anything else.
/// The pole of `ball` a ring walked as `ring` turns toward: the side its
/// material lies on, left of the walk about the face's own normal, which is
/// the surface's (its chart's `du x dv`, inward on a ball framed
/// left-handed) where the face runs with its surface. `None` where the ring
/// has no curve to walk or stands square to neither pole.
fn pole_on_the_left(
    model: &Model,
    surface: &SurfaceGeometry,
    ball: ogeom_math::Sphere,
    ring: &Shape,
    at: &Shape,
    face_forward: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    let Some(data) = model.node(ring).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(None);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        return Ok(None);
    };
    let Some(start) = model
        .node(at)
        .and_then(|n| n.data().as_vertex())
        .map(|v| v.point)
    else {
        return Ok(None);
    };
    let reversed = ring.orientation() == ogeom_topo::Orientation::Reversed;
    let t = if reversed { range.1 } else { range.0 };
    let along = curve.d1_at(t, tol)? * if reversed { -1.0 } else { 1.0 };
    let local = ball.frame().to_local(start);
    let (u, v) = (
        local.y.atan2(local.x),
        (local.z / ball.radius()).clamp(-1.0, 1.0).asin(),
    );
    let (du, dv) = surface.d1_at(u, v, tol)?;
    let own = du.cross(dv);
    let normal = if face_forward { own } else { -own };
    let left = normal.cross(along);
    let axis = ball.frame().z().vector() * ball.radius();
    let lean = left.dot(axis);

    if lean.abs() <= tol.confusion() * left.magnitude() {
        return Ok(None);
    }
    Ok(Some(if lean > 0.0 {
        ball.centre() + axis
    } else {
        ball.centre() - axis
    }))
}

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

/// The window of `curve` between the feet of two vertices that stand on it
/// short of its ends: both within a hair of the curve, in the curve's own
/// order, and not the whole curve. `None` where the vertices are off the
/// curve, reversed against it, or already at its ends.
fn window_between_feet(
    curve: &Curve,
    start: Point,
    end: Point,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let (lo, hi) = curve.domain();
    let samples = match curve {
        Curve::BSpline(spline) => (spline.control_points().len() * 8).max(64),
        _ => 64,
    };
    let foot = |p: Point| -> Option<f64> {
        let found = project_on_curve(curve, p, samples, tol).ok()?;
        (found.distance <= tol.confusion() * 1e4).then_some(found.parameter)
    };
    let (a, b) = (foot(start)?, foot(end)?);
    let slack = tol.parametric().max((hi - lo) * 1e-9);
    if b <= a + slack {
        return None;
    }
    if (a - lo).abs() <= slack && (b - hi).abs() <= slack {
        return None;
    }
    Some((a, b))
}

/// Which instances the reader has looked at.
///
/// One slot per possible id where the instance table has one, as it does
/// for the dense ids exporters write: the reader touches instances millions
/// of times, and a direct index beats hashing every touch. A file whose ids
/// are sparse (one instance numbered in the billions) gets a set instead of
/// a table that size.
enum Visited {
    Dense(Vec<bool>),
    Sparse(ogeom_core::FastSet<u64>),
}

impl Visited {
    fn for_instances(instances: &Instances) -> Self {
        match instances {
            Instances::Dense(slots) => Self::Dense(vec![false; slots.len()]),
            Instances::Sparse(_) => Self::Sparse(ogeom_core::FastSet::default()),
        }
    }

    fn mark(&mut self, id: u64) {
        match self {
            Self::Dense(slots) => {
                if let Some(slot) = usize::try_from(id).ok().and_then(|i| slots.get_mut(i)) {
                    *slot = true;
                }
            }
            Self::Sparse(seen) => {
                seen.insert(id);
            }
        }
    }

    fn seen(&self, id: u64) -> bool {
        match self {
            Self::Dense(slots) => usize::try_from(id)
                .ok()
                .and_then(|i| slots.get(i).copied())
                .unwrap_or(false),
            Self::Sparse(seen) => seen.contains(&id),
        }
    }
}

/// The reference in argument `i`, or 0 (which names no entity, so the lookup
/// that follows refuses it) where the argument is missing or not a reference.
/// A representation item's name: its `REPRESENTATION_ITEM` part's in a
/// complex instance, else its own first argument's.
fn item_name(instance: &Instance) -> String {
    instance
        .part("REPRESENTATION_ITEM")
        .or_else(|| instance.parts.first().map(|(_, args)| args.as_slice()))
        .and_then(<[Arg]>::first)
        .and_then(|arg| match arg {
            Arg::Str(name) => Some(name.clone()),
            _ => None,
        })
        .unwrap_or_default()
}

fn ref_at(args: &[Arg], i: usize) -> u64 {
    args.get(i).and_then(Arg::reference).unwrap_or(0)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests;
