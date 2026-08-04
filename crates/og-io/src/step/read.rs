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
use og_algo::{make_edge_between, make_face_on, make_shell, make_solid, make_vertex, make_wire};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve2d as _;
use og_geom::Curve3d as _;
use og_geom::Surface as _;
use og_geom::{
    BSplineCurve, CircleCurve, ConeSurface, Curve, CylinderSurface, EllipseCurve, LineCurve,
    PlanarCurve, PlaneSurface, SphereSurface, SurfaceGeometry, TorusSurface,
};
use og_math::{
    Axis, Circle, Cone, Cylinder, Direction, Ellipse, Frame, KnotVector, Plane, Point, Sphere,
    Torus, Transform2, Vector, Vector2,
};
use og_topo::{Location, Model, Shape};
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
}

/// A read exchange file: the model, the solids found, and the report.
#[derive(Debug)]
pub struct StepImport {
    /// The model everything was built into.
    pub model: Model,
    /// One shape per `MANIFOLD_SOLID_BREP`, in file order.
    pub solids: Vec<Shape>,
    /// What happened along the way.
    pub report: StepReport,
}

/// Read a STEP exchange file's B-rep content.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the file does
/// not parse, contains no solid, or a solid's structure is broken in a way
/// topology cannot represent. Faces the reader cannot complete become
/// warnings, not errors — the report says exactly what was compromised.
pub fn read_step(text: &str, tol: Tolerances) -> OgResult<StepImport> {
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
        tol,
    };
    reader.report.scale_mm = reader.unit_scale();
    reader.angle_scale = reader.angle_unit_scale();

    let mut solids = Vec::new();
    let mut ids: Vec<u64> = exchange
        .data
        .iter()
        .filter(|(_, inst)| inst.part("MANIFOLD_SOLID_BREP").is_some())
        .map(|(id, _)| *id)
        .collect();
    ids.sort_unstable();
    for id in ids {
        solids.push(reader.solid(id)?);
    }
    if solids.is_empty() {
        og_bail!(
            Construction,
            "the exchange file contains no MANIFOLD_SOLID_BREP to read"
        );
    }

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
        model: reader.model,
        solids,
        report: reader.report,
    })
}

/// An edge as built: the shape, its curve and range, and whether STEP's edge
/// direction runs against the curve.
type BuiltEdge = (Shape, Curve, (f64, f64), bool);

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
    tol: Tolerances,
}

impl Reader<'_> {
    fn instance(&mut self, id: u64) -> OgResult<&'_ Instance> {
        self.visited.insert(id);
        self.exchange.data.get(&id).ok_or_else(|| {
            og_core::og_err!(
                Construction,
                "the file references #{id}, which does not exist"
            )
        })
    }

    fn args(&mut self, id: u64, keyword: &str) -> OgResult<Vec<Arg>> {
        let instance = self.instance(id)?;
        let Some(args) = instance.part(keyword) else {
            og_bail!(
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
        let mut out = Vec::new();
        for instance in self.exchange.data.values() {
            if let Some(args) = instance.part("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
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

    fn point(&mut self, id: u64) -> OgResult<Point> {
        let args = self.args(id, "CARTESIAN_POINT")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            og_bail!(Construction, "#{id}: a point without coordinates");
        };
        let scale = self.report.scale_mm;
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0) * scale;
        Ok(Point::new(value(0), value(1), value(2)))
    }

    fn direction(&mut self, id: u64) -> OgResult<Direction> {
        let args = self.args(id, "DIRECTION")?;
        let Some(coords) = args.get(1).and_then(Arg::list) else {
            og_bail!(Construction, "#{id}: a direction without components");
        };
        let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0);
        Direction::new(Vector::new(value(0), value(1), value(2)), self.tol)
    }

    /// An `AXIS2_PLACEMENT_3D` as a frame, with the standard's defaults for
    /// what the file leaves out.
    fn frame(&mut self, id: u64) -> OgResult<Frame> {
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

    fn surface(&mut self, id: u64) -> OgResult<Option<SurfaceGeometry>> {
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
    ) -> OgResult<SurfaceGeometry> {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let deg = |a: Option<Arg>| a.and_then(|x| x.number()).unwrap_or(1.0) as usize;
        let (u_degree, v_degree) = (deg(degrees.0), deg(degrees.1));
        let Some(Arg::List(rows)) = grid_arg else {
            og_bail!(
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
            let weighted: Vec<og_math::Weighted<Point>> = points
                .iter()
                .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                .map(|(p, w)| og_math::Weighted::new(*p, *w, self.tol))
                .collect::<OgResult<_>>()?;
            og_geom::BSplineSurface::rational(
                u_knots,
                v_knots,
                og_math::ControlGrid::new(weighted, u_count, v_count)?,
            )?
        } else {
            og_geom::BSplineSurface::new(
                u_knots,
                v_knots,
                &og_math::ControlGrid::new(points, u_count, v_count)?,
                self.tol,
            )?
        };
        Ok(surface.into())
    }

    fn curve(&mut self, id: u64) -> OgResult<Option<Curve>> {
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
                    .collect::<OgResult<_>>()?;
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
                let weighted: Vec<og_math::Weighted<Point>> = control
                    .iter()
                    .zip(flat.iter().chain(std::iter::repeat(&1.0)))
                    .map(|(p, w)| og_math::Weighted::new(*p, *w, self.tol))
                    .collect::<OgResult<_>>()?;
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
                    .collect::<OgResult<_>>()?;
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

    fn vertex(&mut self, id: u64) -> OgResult<Shape> {
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
    fn edge(&mut self, id: u64) -> OgResult<Option<BuiltEdge>> {
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
                    && let og_topo::NodeData::Vertex(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(gap + self.tol.confusion());
                }
                self.report.warnings.push(format!(
                    "#{id}: a curve end misses its vertex by {gap:.2e}; the                      vertex tolerance grew to say so"
                ));
            }
        }
        let shape =
            make_edge_between(&mut self.model, curve.clone(), range, &vlo, &vhi, self.tol)?.shape;
        let entry = (shape, curve, range, flipped);
        self.edges.insert(id, entry.clone());
        Ok(Some(entry))
    }

    fn face(&mut self, id: u64) -> OgResult<Option<Shape>> {
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
                let projection = og_algo::project_on_surface(&surface, at, 32, self.tol)?;
                let ((ua, ub), _) = surface.domain();
                let row = projection.parameters.1;
                let mut data = og_topo::EdgeData::new();
                data.degenerate = true;
                let edge = self
                    .model
                    .add_edge(data, &[vertex.clone(), vertex.clone()])?;
                let pcurve: PlanarCurve = og_geom::Line2d::segment(
                    og_math::Point2::new(ua, row),
                    og_math::Point2::new(ub, row),
                    self.tol,
                )?
                .into();
                og_algo::attach_pcurve(
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
                // The band construction is og-algo's make_revolution_band —
                // one authority shared with the healer. Anything it refuses
                // (a spline ring, ring vertices at different angles, a
                // surface with no closed-form iso-curve) becomes a warning
                // and the raw bounds, not an error.
                match og_algo::make_revolution_band(&mut self.model, &surface, e_lo, e_hi, self.tol)
                {
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
    fn bound_args(&mut self, id: u64) -> OgResult<(u64, bool)> {
        let instance = self.instance(id)?;
        let args = instance
            .part("FACE_OUTER_BOUND")
            .or_else(|| instance.part("FACE_BOUND"))
            .ok_or_else(|| og_core::og_err!(Construction, "#{id} is not a face bound"))?
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
        edge: &Shape,
        curve: &Curve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
        surface_id: og_topo::SurfaceId,
        seam: bool,
    ) -> OgResult<()> {
        let widen = |p: PlanarCurve| -> PlanarCurve {
            // A line pcurve evaluates anywhere; its stated domain must still
            // cover the edge's range, which for a wrapped circle runs past
            // one period.
            if let PlanarCurve::Line(l) = &p {
                let (lo, hi) = (l.domain().0.min(range.0), l.domain().1.max(range.1));
                if let Ok(wider) = og_geom::Line2d::over(l.axis(), lo, hi) {
                    return wider.into();
                }
            }
            p
        };
        let pcurve = match og_intersect::exact_pcurve_of(curve, surface, self.tol).map(widen) {
            Some(exact) => exact,
            None => {
                // No closed form — a spline surface, or a combination the
                // projection table lacks. The pcurve is *fitted at the
                // curve's own parameters*: sample the edge, project each
                // sample into the chart, fit the trace with the parameters
                // held fixed, so same-parameter is preserved by construction
                // and the reported error is the true chart deviation.
                match self.fit_projected_pcurve(curve, range, surface) {
                    Ok((fitted, error, met)) => {
                        if !met {
                            self.report.warnings.push(format!(
                                "face #{face_id}: a projected pcurve fit \
                                 stopped at {error:.2e}; the face's mesh may \
                                 sit that far off along this edge"
                            ));
                        }
                        fitted
                    }
                    Err(e) => {
                        self.report.warnings.push(format!(
                            "face #{face_id}: no pcurve for an edge on this \
                             surface ({e}); the face may not triangulate"
                        ));
                        return Ok(());
                    }
                }
            }
        };
        if seam {
            let ((ua, ub), _) = surface.domain();
            let span = ub - ua;
            // One side is where the projection landed; the other is one
            // period over, whichever way stays within the chart.
            let mid = pcurve.point_at(f64::midpoint(range.0, range.1), self.tol)?;
            let shift = if mid.x - ua < span * 0.5 { span } else { -span };
            let other =
                pcurve.transformed(&Transform2::translation(Vector2::new(shift, 0.0)), self.tol)?;
            og_algo::attach_seam(
                &mut self.model,
                edge,
                pcurve,
                other,
                surface_id,
                Location::identity(),
                range,
            )?;
        } else {
            og_algo::attach_pcurve(
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

    /// A pcurve fitted from projection at the curve's own parameters.
    fn fit_projected_pcurve(
        &mut self,
        curve: &Curve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
    ) -> OgResult<(PlanarCurve, f64, bool)> {
        const SAMPLES: usize = 96;
        let mut worst_off = 0.0_f64;
        let mut parameters = Vec::with_capacity(SAMPLES + 1);
        let mut trace = Vec::with_capacity(SAMPLES + 1);
        let mut space_run = 0.0;
        let mut parameter_run = 0.0;
        let mut previous: Option<(Point, og_math::Point2)> = None;
        for i in 0..=SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let t = range.0 + (range.1 - range.0) * i as f64 / SAMPLES as f64;
            let p = curve.point_at(t, self.tol)?;
            let (uv, off) = match chart_of(surface, p) {
                // Analytic surfaces invert in closed form — grid seeding
                // over a plane's or cylinder's enormous stated extents lands
                // microns off, and a fitted pcurve inherits every micron.
                Some(uv) => {
                    let lifted = surface.point_at(uv.x, uv.y, self.tol)?;
                    (uv, p.distance(lifted))
                }
                None => {
                    let projection = og_algo::project_on_surface(surface, p, 24, self.tol)?;
                    (
                        og_math::Point2::new(projection.parameters.0, projection.parameters.1),
                        projection.distance,
                    )
                }
            };
            if off > self.tol.confusion() * 1e4 {
                og_bail!(
                    Construction,
                    "the edge sits {off:.2e} from the surface it should bound"
                );
            }
            worst_off = worst_off.max(off);
            if let Some((lp, luv)) = previous {
                space_run += p.distance(lp);
                parameter_run += uv.distance(luv);
            }
            previous = Some((p, uv));
            parameters.push(t);
            trace.push(uv);
        }
        // A trace on a periodic chart may cross the seam mid-edge; unwrap it
        // pointwise so the fit sees a continuous curve.
        let ((ua, ub), (va, vb)) = surface.domain();
        let spans = (
            if surface.is_periodic_u() {
                ub - ua
            } else {
                0.0
            },
            if surface.is_periodic_v() {
                vb - va
            } else {
                0.0
            },
        );
        for i in 1..trace.len() {
            if spans.0 > 0.0 {
                while trace[i].x - trace[i - 1].x > spans.0 * 0.5 {
                    trace[i].x -= spans.0;
                }
                while trace[i].x - trace[i - 1].x < -spans.0 * 0.5 {
                    trace[i].x += spans.0;
                }
            }
            if spans.1 > 0.0 {
                while trace[i].y - trace[i - 1].y > spans.1 * 0.5 {
                    trace[i].y -= spans.1;
                }
                while trace[i].y - trace[i - 1].y < -spans.1 * 0.5 {
                    trace[i].y += spans.1;
                }
            }
        }
        // The tolerance carried into the chart through the trace's own
        // metric — the honest cheap version, refined by the fit's report.
        let scale = if space_run > self.tol.confusion() {
            parameter_run / space_run
        } else {
            1.0
        };
        let target = (self.tol.confusion() * 1e2 * scale).max(f64::MIN_POSITIVE);
        let fitted = og_geom::fit::fit_points_2d_at(&parameters, &trace, 3, target, self.tol)?;
        if worst_off > self.tol.confusion() * 1e3 {
            self.report.warnings.push(format!(
                "an edge sits up to {worst_off:.2e} from the surface it \
                 bounds; the file's own slop, carried into the chart"
            ));
        }
        Ok((fitted.curve.into(), fitted.error, fitted.met))
    }

    fn solid(&mut self, id: u64) -> OgResult<Shape> {
        let args = self.args(id, "MANIFOLD_SOLID_BREP")?;
        let shell_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let shell_instance = self.instance(shell_id)?;
        let shell_args = shell_instance
            .part("CLOSED_SHELL")
            .or_else(|| shell_instance.part("OPEN_SHELL"))
            .ok_or_else(|| og_core::og_err!(Construction, "#{shell_id} is not a shell"))?
            .to_vec();
        let face_ids: Vec<u64> = shell_args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        let mut faces = Vec::new();
        for fid in face_ids {
            if let Some(face) = self.face(fid)? {
                faces.push(face);
            }
        }
        if faces.is_empty() {
            og_bail!(Construction, "#{id}: a solid with no readable faces");
        }
        let shell = make_shell(&mut self.model, &faces)?.shape;
        if !og_algo::is_shell_closed(&self.model, &shell)? {
            self.report.warnings.push(format!(
                "#{id}: the shell does not close as read; measures needing an \
                 inside will refuse it"
            ));
        }
        Ok(make_solid(&mut self.model, std::slice::from_ref(&shell))?.shape)
    }
}

/// The chart coordinates of a point on an analytic surface, by closed-form
/// inversion — `None` for surfaces that need iterative projection.
fn chart_of(surface: &SurfaceGeometry, p: Point) -> Option<og_math::Point2> {
    let tau = core::f64::consts::TAU;
    match surface {
        SurfaceGeometry::Plane(s) => {
            let l = s.plane().frame().to_local(p);
            Some(og_math::Point2::new(l.x, l.y))
        }
        SurfaceGeometry::Cylinder(s) => {
            let l = s.cylinder().frame().to_local(p);
            Some(og_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), l.z))
        }
        SurfaceGeometry::Cone(s) => {
            let l = s.cone().frame().to_local(p);
            Some(og_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), l.z))
        }
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let l = sphere.frame().to_local(p);
            let lat = (l.z / sphere.radius()).clamp(-1.0, 1.0).asin();
            Some(og_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), lat))
        }
        SurfaceGeometry::Torus(s) => {
            let torus = s.torus();
            let l = torus.frame().to_local(p);
            let u = l.y.atan2(l.x).rem_euclid(tau);
            let radial = l.x.hypot(l.y) - torus.major_radius();
            let v = l.z.atan2(radial).rem_euclid(tau);
            Some(og_math::Point2::new(u, v))
        }
        _ => None,
    }
}

/// For a two-wire periodic face: each wire's single closed edge with its
/// vertex, empty when the shape is anything else.
fn closed_ring_edges(model: &Model, wires: &[Shape]) -> OgResult<Vec<(Shape, Shape)>> {
    let mut out = Vec::new();
    for wire in wires {
        let edges = og_topo::explore(
            model,
            wire,
            og_topo::Filter::OfType(og_topo::ShapeType::Edge),
        )?;
        if edges.len() != 1 {
            return Ok(Vec::new());
        }
        let edge = edges[0].clone();
        let Some((a, b)) = og_algo::edge_vertices(model, &edge)? else {
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
