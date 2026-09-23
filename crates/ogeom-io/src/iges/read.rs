//! From an IGES deck to a living model.
//!
//! Two kinds of file arrive under one extension. A *solid* file carries
//! manifold solid B-rep objects — entity 186 over shells, faces, loops, edge
//! lists and vertex lists — and reads bottom-up the way the STEP reader does,
//! sharing what the file shares. A *surface* file, the older and far more
//! common kind, is a loose collection of trimmed surfaces; those become
//! faces, the faces are sewn, and a shell that closes becomes a solid. Both
//! kinds re-derive edge ranges on this kernel's own parameterizations from
//! the endpoint geometry, because a 1980s file's parameterizations are its
//! own business.
//!
//! What the reader does not understand it *counts*: every entity never
//! visited lands in the report's skipped table under its type number, and
//! every compromise is a warning naming the directory entry. Refusals are by
//! name — a conic form this reader does not translate says which form and
//! where, and points at the parity ledger's `io.iges` row for the whole
//! picture.

use super::parse::{Entity, File};
use ogeom_algo::{make_edge_between, make_solid, make_vertex, sew};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{
    BSplineCurve, BSplineSurface, CircleCurve, ConeSurface, Curve, CylinderSurface, EllipseCurve,
    ExtrusionSurface, HyperbolaCurve, LineCurve, OffsetCurve, OffsetSurface, ParabolaCurve,
    PlaneSurface, RevolutionSurface, SphereSurface, SurfaceGeometry, TorusSurface, TrimmedCurve,
};
use ogeom_math::{
    Circle, Cone, ControlGrid, Cylinder, Direction, Ellipse, Frame, Hyperbola, KnotVector, Matrix3,
    Parabola, Plane, Point, Sphere, Torus, Transform, Vector, Weighted,
};
use ogeom_topo::{Model, Shape};
use std::collections::{BTreeMap, HashMap};

/// How far an unbounded plane or quadric extends past anything the file uses
/// — the same convention the STEP reader states: a face's trim is its wires,
/// and the surface's domain is only a parameter window.
const SURFACE_EXTENT: f64 = 1e5;

/// What an import brought in, and what it left behind.
#[derive(Debug, Default)]
pub struct IgesReport {
    /// Millimetres per file unit, as the global section states it.
    pub scale_mm: f64,
    /// Entity types the reader never visited, with counts, keyed as
    /// `"type NNN"` or `"type NNN form F"`. Annotation and drafting land
    /// here by design; geometry landing here is a gap worth reading about.
    pub skipped: BTreeMap<String, usize>,
    /// Everything that imported less than perfectly, one line each.
    pub warnings: Vec<String>,
}

/// A read IGES file: the document, the shapes found, and the report.
#[derive(Debug)]
pub struct IgesImport {
    /// The document everything was built into.
    pub document: ogeom_doc::Document,
    /// One shape per manifold solid, then one per surface group that sewed
    /// closed.
    pub solids: Vec<Shape>,
    /// Sewn shells and loose faces that do not enclose a volume.
    pub sheets: Vec<Shape>,
    /// What happened along the way.
    pub report: IgesReport,
}

/// The pieces an edge-list entry resolves to.
type BuiltEdge = (Shape, Curve, (f64, f64));

/// Read an IGES file's geometry.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// deck does not parse, its units are unreadable, or it contains nothing this
/// reader translates into shapes. Individual entities that fail to translate
/// become warnings and skipped counts rather than errors — the report says
/// exactly what was compromised.
pub fn read_iges(text: &str, tol: Tolerances) -> OgeomResult<IgesImport> {
    let file = super::parse::parse(text)?;
    let Some(scale_mm) = file.scale_mm() else {
        ogeom_bail!(
            Construction,
            "IGES global section names units this reader cannot convert to millimetres"
        );
    };
    let mut reader = Reader {
        file: &file,
        model: Model::new(),
        report: IgesReport {
            scale_mm,
            ..IgesReport::default()
        },
        visited: BTreeMap::new(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        vertex_misses: (0, 0.0),
        tol,
    };

    // Solids first: every 186, in directory order.
    let mut solids: Vec<(i64, Shape)> = Vec::new();
    let solid_des: Vec<i64> = file
        .entities
        .iter()
        .filter(|(_, e)| e.kind == 186)
        .map(|(de, _)| *de)
        .collect();
    let total = solid_des.len() as u64;
    for (done, de) in solid_des.into_iter().enumerate() {
        ogeom_core::progress::checkpoint()?;
        ogeom_core::progress::stage_at("iges: solid", done as u64 + 1, total);
        match reader.manifold_solid(de) {
            Ok(solid) => solids.push((de, solid)),
            Err(e) => reader
                .report
                .warnings
                .push(format!("D{de}: manifold solid failed to build: {e}")),
        }
    }

    // Then the surface file: every *independent* trimmed or bounded surface
    // becomes a face; the faces sew; closed shells become solids. Subordinate
    // entities belong to something else and are not top-level geometry — the
    // subordinate switch is the second two-digit field of the status word.
    let mut faces = Vec::new();
    let face_des: Vec<i64> = file
        .entities
        .iter()
        .filter(|(de, e)| {
            matches!(e.kind, 143 | 144)
                && (e.status / 10_000) % 100 == 0
                && !reader.visited.contains_key(de)
        })
        .map(|(de, _)| *de)
        .collect();
    let total = face_des.len() as u64;
    for (done, de) in face_des.into_iter().enumerate() {
        ogeom_core::progress::checkpoint()?;
        ogeom_core::progress::stage_at("iges: face", done as u64 + 1, total);
        match reader.face(de) {
            Ok(face) => faces.push(face),
            Err(e) => reader
                .report
                .warnings
                .push(format!("D{de}: trimmed surface failed to build: {e}")),
        }
    }
    let mut sheets = Vec::new();
    if !faces.is_empty() {
        let sewn = sew(&mut reader.model, &faces, tol)?;
        for shell in &sewn.shells {
            if ogeom_algo::is_shell_closed(&reader.model, shell)? {
                let solid = make_solid(&mut reader.model, std::slice::from_ref(shell))?.shape;
                solids.push((0, solid));
            } else {
                sheets.push(shell.clone());
            }
        }
    }

    if solids.is_empty() && sheets.is_empty() {
        // A deck that *had* candidates which all failed is a different
        // refusal from one with nothing to try, and the warnings say why.
        if reader.report.warnings.is_empty() {
            ogeom_bail!(
                Construction,
                "the IGES file contains no manifold solid and no independent \
                 trimmed surface this reader translates"
            );
        }
        ogeom_bail!(
            Construction,
            "every shape in the IGES file failed to build: {}",
            reader.report.warnings.join("; ")
        );
    }

    // Everything never visited, counted by type and form.
    for (de, entity) in &file.entities {
        if !reader.visited.contains_key(de) {
            let key = if entity.form == 0 {
                format!("type {}", entity.kind)
            } else {
                format!("type {} form {}", entity.kind, entity.form)
            };
            *reader.report.skipped.entry(key).or_default() += 1;
        }
    }

    if reader.vertex_misses.0 > 0 {
        let (count, worst) = reader.vertex_misses;
        reader.report.warnings.push(format!(
            "{count} vertices sat off the curve ends they bound, by up to {worst:.2e}; \
             their tolerances grew to say so"
        ));
    }
    let document = reader.document(&solids, &sheets);
    let solids = solids.into_iter().map(|(_, s)| s).collect();
    Ok(IgesImport {
        document,
        solids,
        sheets,
        report: reader.report,
    })
}

struct Reader<'a> {
    file: &'a File,
    model: Model,
    report: IgesReport,
    /// Every directory entry the reader consumed, for the skipped table.
    visited: BTreeMap<i64, ()>,
    /// Vertices by (vertex-list DE, 1-based index) — shared, which is what
    /// lets a closed shell close.
    vertices: HashMap<(i64, i64), Shape>,
    /// Edges by (edge-list DE, 1-based index), for the same reason.
    edges: HashMap<(i64, i64), BuiltEdge>,
    /// Vertices a curve end missed by more than the confusion tolerance
    /// and that widened to cover it: how many, and the widest miss.
    vertex_misses: (usize, f64),
    tol: Tolerances,
}

impl<'a> Reader<'a> {
    fn entity(&mut self, de: i64) -> OgeomResult<&'a Entity> {
        let Some(entity) = self.file.entity(de) else {
            ogeom_bail!(Construction, "IGES pointer D{de} names no entity");
        };
        self.visited.insert(de.abs(), ());
        Ok(entity)
    }

    /// The model-space transform an entity carries, identity when none. A
    /// transformation entity may itself be transformed; that composes.
    fn placement(&mut self, entity: &Entity) -> OgeomResult<Transform> {
        if entity.transform == 0 {
            return Ok(Transform::IDENTITY);
        }
        let de = entity.transform;
        let t = self.entity(de)?;
        if t.kind != 124 {
            ogeom_bail!(
                Construction,
                "D{de}: a transformation pointer names a type {} entity",
                t.kind
            );
        }
        let v = |i: usize| t.at(i).real();
        let s = self.report.scale_mm;
        let linear = Matrix3::new([[v(0), v(1), v(2)], [v(4), v(5), v(6)], [v(8), v(9), v(10)]]);
        let translation = Vector::new(v(3) * s, v(7) * s, v(11) * s);
        let m = Transform::from_parts(linear, 1.0, translation, self.tol.angular())?;
        if t.transform != 0 {
            let outer = self.placement(t)?;
            return Ok(outer * m);
        }
        Ok(m)
    }

    fn point3(&self, e: &Entity, i: usize) -> Point {
        let s = self.report.scale_mm;
        Point::new(
            e.at(i).real() * s,
            e.at(i + 1).real() * s,
            e.at(i + 2).real() * s,
        )
    }

    /// A model-space curve with the range its own definition covers.
    fn curve(&mut self, de: i64) -> OgeomResult<(Curve, (f64, f64))> {
        let entity = self.entity(de)?;
        let scale = self.report.scale_mm;
        let (curve, range) = match entity.kind {
            110 => {
                let a = self.point3(entity, 0);
                let b = self.point3(entity, 3);
                let line = LineCurve::segment(a, b, self.tol)?;
                (Curve::from(line), (0.0, a.distance(b)))
            }
            100 => {
                // Centre, start and end in the definition plane at z = zt;
                // the arc runs counter-clockwise from start to end.
                let zt = entity.at(0).real() * scale;
                let c = Point::new(entity.at(1).real() * scale, entity.at(2).real() * scale, zt);
                let s = Point::new(entity.at(3).real() * scale, entity.at(4).real() * scale, zt);
                let e = Point::new(entity.at(5).real() * scale, entity.at(6).real() * scale, zt);
                let radius = c.distance(s);
                let x = Direction::new(s - c, self.tol)?;
                let frame = Frame::new(c, Direction::Z, x, self.tol)?;
                let circle = Circle::new(frame, radius, self.tol)?;
                let to_end = e - c;
                let ang = to_end
                    .dot(frame.y().vector())
                    .atan2(to_end.dot(frame.x().vector()))
                    .rem_euclid(core::f64::consts::TAU);
                let sweep = if ang <= self.tol.parametric() {
                    core::f64::consts::TAU
                } else {
                    ang
                };
                (Curve::from(CircleCurve::new(circle)), (0.0, sweep))
            }
            104 => self.conic(de, entity)?,
            112 => self.spline_curve(de, entity)?,
            126 => self.nurbs_curve(de, entity)?,
            130 => self.offset_curve(de, entity)?,
            102 => ogeom_bail!(
                Construction,
                "D{de}: a composite curve is a sequence, not a curve; the \
                 caller walks its segments"
            ),
            kind => ogeom_bail!(
                Construction,
                "D{de}: curve entity type {kind}{} is not translated — \
                 docs/PARITY.md, io.iges",
                super::entity_name(kind).map_or_else(String::new, |n| format!(" ({n})"))
            ),
        };
        // The entity's placement moves the curve into model space.
        let placement = self.placement(entity)?;
        if placement == Transform::IDENTITY {
            Ok((curve, range))
        } else {
            Ok((curve.transformed(&placement, self.tol)?, range))
        }
    }

    /// Conic arc: `A x² + B xy + C y² + D x + E y + F = 0` in the definition
    /// plane, axis-aligned. The coefficients say which conic it is — both
    /// squares one sign an ellipse, opposite signs a hyperbola, one square
    /// missing a parabola — and each translates to its own curve, the arc's
    /// ends read off the start and terminate points. A rotated conic is
    /// refused by name until a file demands it.
    fn conic(&mut self, de: i64, entity: &Entity) -> OgeomResult<(Curve, (f64, f64))> {
        let scale = self.report.scale_mm;
        let (a, b, c, d, e, f) = (
            entity.at(0).real(),
            entity.at(1).real(),
            entity.at(2).real(),
            entity.at(3).real(),
            entity.at(4).real(),
            entity.at(5).real(),
        );
        if b.abs() > 1e-12 {
            ogeom_bail!(
                Construction,
                "D{de}: conic arc form {} turns its axes; only an axis-aligned \
                 conic is translated — docs/PARITY.md, io.iges",
                entity.form
            );
        }
        let zt = entity.at(6).real() * scale;
        let start = Point::new(entity.at(7).real() * scale, entity.at(8).real() * scale, zt);
        let end = Point::new(
            entity.at(9).real() * scale,
            entity.at(10).real() * scale,
            zt,
        );
        // Both squares present and of one sign: an ellipse, its coefficients
        // made positive.
        if a != 0.0 && c != 0.0 && (a > 0.0) == (c > 0.0) {
            let sign = if a > 0.0 { 1.0 } else { -1.0 };
            return self.ellipse_arc(de, [a, c, d, e, f].map(|k| k * sign), zt, start, end);
        }
        if a != 0.0 && c != 0.0 {
            return Self::hyperbola_arc(de, [a, c, d, e, f], zt, start, end, scale, self.tol);
        }
        if (a == 0.0) != (c == 0.0) {
            return Self::parabola_arc(de, [a, c, d, e, f], zt, start, end, scale, self.tol);
        }
        ogeom_bail!(Construction, "D{de}: conic arc coefficients close no conic")
    }

    /// The ellipse arm of [`Reader::conic`]: `A x² + C y² + D x + E y + F = 0`
    /// with `A` and `C` positive.
    fn ellipse_arc(
        &mut self,
        de: i64,
        [a, c, d, e, f]: [f64; 5],
        zt: f64,
        start: Point,
        end: Point,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        let scale = self.report.scale_mm;
        let cx = -d / (2.0 * a);
        let cy = -e / (2.0 * c);
        let rhs = a * cx * cx + c * cy * cy - f;
        let (ra2, rb2) = (rhs / a, rhs / c);
        if ra2 <= 0.0 || rb2 <= 0.0 {
            ogeom_bail!(
                Construction,
                "D{de}: conic arc coefficients close no ellipse"
            );
        }
        let centre = Point::new(cx * scale, cy * scale, zt);
        let (rx, ry) = (ra2.sqrt() * scale, rb2.sqrt() * scale);
        // The ellipse type wants major ≥ minor; when the x semi-axis is the
        // smaller one, a quarter-turn of the frame swaps the roles.
        let (frame, major, minor) = if rx >= ry {
            (
                Frame::new(centre, Direction::Z, Direction::X, self.tol)?,
                rx,
                ry,
            )
        } else {
            (
                Frame::new(centre, Direction::Z, Direction::Y, self.tol)?,
                ry,
                rx,
            )
        };
        let ellipse = Ellipse::new(frame, major, minor, self.tol)?;
        let angle_of = |p: Point| -> f64 {
            let local = frame.to_local(p);
            (local.y / minor)
                .atan2(local.x / major)
                .rem_euclid(core::f64::consts::TAU)
        };
        let t0 = angle_of(start);
        let mut t1 = angle_of(end);
        if t1 <= t0 + self.tol.parametric() {
            t1 += core::f64::consts::TAU;
        }
        Ok((Curve::from(EllipseCurve::new(ellipse)), (t0, t1)))
    }

    /// The hyperbola arm of [`Reader::conic`]: squares of opposite sign.
    /// The branch is the one the arc's ends stand on, its parameter the
    /// natural one — `(a cosh t, b sinh t)` — and the arc runs in increasing
    /// parameter, the frame turned over where the file's ends run the other
    /// way.
    fn hyperbola_arc(
        de: i64,
        [a, c, d, e, f]: [f64; 5],
        zt: f64,
        start: Point,
        end: Point,
        scale: f64,
        tol: Tolerances,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        let cx = -d / (2.0 * a);
        let cy = -e / (2.0 * c);
        let rhs = a * cx * cx + c * cy * cy - f;
        if rhs == 0.0 {
            ogeom_bail!(
                Construction,
                "D{de}: conic arc coefficients close no hyperbola; they cross"
            );
        }
        // Divided through by the right-hand side, the positive square names
        // the transverse axis.
        let (pa, pc) = (a / rhs, c / rhs);
        let (major, minor, along_x) = if pa > 0.0 {
            ((1.0 / pa).sqrt(), (-1.0 / pc).sqrt(), true)
        } else {
            ((1.0 / pc).sqrt(), (-1.0 / pa).sqrt(), false)
        };
        let centre = Point::new(cx * scale, cy * scale, zt);
        let axis = if along_x { Direction::X } else { Direction::Y };
        // The branch: where the ends stand along the transverse axis.
        let side = (start - centre).dot(axis.vector());
        let x = if side >= 0.0 {
            axis
        } else {
            Direction::new(-axis.vector(), tol)?
        };
        let build = |z: Direction| -> OgeomResult<(Curve, (f64, f64))> {
            let frame = Frame::new(centre, z, x, tol)?;
            let (major, minor) = (major * scale, minor * scale);
            let hyperbola = Hyperbola::new(frame, major, minor, tol)?;
            let t_of = |p: Point| (frame.to_local(p).y / minor).asinh();
            let (t0, t1) = (t_of(start), t_of(end));
            let extent = t0.abs().max(t1.abs()).max(1e-3) * 1.5;
            Ok((
                Curve::from(HyperbolaCurve::new(hyperbola, extent)?),
                (t0, t1),
            ))
        };
        let (curve, (t0, t1)) = build(Direction::Z)?;
        if t1 > t0 {
            return Ok((curve, (t0, t1)));
        }
        // Turned over: the same branch, the parameter running the other way.
        let (curve, (t0, t1)) = build(Direction::new(-Direction::Z.vector(), tol)?)?;
        Ok((curve, (t0, t1)))
    }

    /// The parabola arm of [`Reader::conic`]: one square missing. The axis
    /// is the missing square's direction; the parameter runs along the
    /// other, as this vocabulary's parabola does.
    fn parabola_arc(
        de: i64,
        [a, c, d, e, f]: [f64; 5],
        zt: f64,
        start: Point,
        end: Point,
        scale: f64,
        tol: Tolerances,
    ) -> OgeomResult<(Curve, (f64, f64))> {
        // `A x² + D x + E y + F = 0` opens along y; `C y² + D x + E y + F = 0`
        // along x. Either is `axis = apex + k (across - apex)²`.
        let (apex_across, apex_along, k, along) = if c == 0.0 {
            if e == 0.0 {
                ogeom_bail!(
                    Construction,
                    "D{de}: conic arc coefficients close no parabola"
                );
            }
            let x0 = -d / (2.0 * a);
            let y0 = -(a * x0 * x0 + d * x0 + f) / e;
            (x0, y0, -a / e, Direction::Y)
        } else {
            if d == 0.0 {
                ogeom_bail!(
                    Construction,
                    "D{de}: conic arc coefficients close no parabola"
                );
            }
            let y0 = -e / (2.0 * c);
            let x0 = -(c * y0 * y0 + e * y0 + f) / d;
            (y0, x0, -c / d, Direction::X)
        };
        let apex = if c == 0.0 {
            Point::new(apex_across * scale, apex_along * scale, zt)
        } else {
            Point::new(apex_along * scale, apex_across * scale, zt)
        };
        // The axis points the way the parabola opens.
        let x = if k >= 0.0 {
            along
        } else {
            Direction::new(-along.vector(), tol)?
        };
        // `axis = k across²` against `axis = across² / (4 f)`.
        let focal = scale / (4.0 * k.abs());
        let build = |z: Direction| -> OgeomResult<(Curve, (f64, f64))> {
            let frame = Frame::new(apex, z, x, tol)?;
            let parabola = Parabola::new(frame, focal, tol)?;
            let t_of = |p: Point| frame.to_local(p).y;
            let (t0, t1) = (t_of(start), t_of(end));
            let extent = t0.abs().max(t1.abs()).max(1e-3) * 1.5;
            Ok((Curve::from(ParabolaCurve::new(parabola, extent)?), (t0, t1)))
        };
        let (curve, (t0, t1)) = build(Direction::Z)?;
        if t1 > t0 {
            return Ok((curve, (t0, t1)));
        }
        let (curve, (t0, t1)) = build(Direction::new(-Direction::Z.vector(), tol)?)?;
        Ok((curve, (t0, t1)))
    }

    /// Offset curve (130): a base curve displaced a constant distance
    /// perpendicular to a reference direction. The file displaces along
    /// the reference crossed with the tangent; this vocabulary's offset runs
    /// along the tangent crossed with the reference, so the distance flips
    /// sign. A varying offset — a function of the parameter — is refused by
    /// name.
    fn offset_curve(&mut self, de: i64, entity: &Entity) -> OgeomResult<(Curve, (f64, f64))> {
        let scale = self.report.scale_mm;
        let kind = entity.at(1).int();
        let (d1, d2) = (entity.at(5).real(), entity.at(7).real());
        if kind != 1 || (d1 - d2).abs() > 1e-12 {
            ogeom_bail!(
                Construction,
                "D{de}: only a constant offset (type 1) is translated — \
                 docs/PARITY.md, io.iges"
            );
        }
        let (basis, range) = self.curve(entity.at(0).int())?;
        let reference = Direction::from_coords(
            entity.at(9).real(),
            entity.at(10).real(),
            entity.at(11).real(),
            self.tol,
        )?;
        let (tt1, tt2) = (entity.at(12).real(), entity.at(13).real());
        let range = if tt2 > tt1 { (tt1, tt2) } else { range };
        let offset = OffsetCurve::new(basis, -d1 * scale, reference)?;
        Ok((Curve::Offset(Box::new(offset)), range))
    }

    fn spline_curve(&mut self, de: i64, entity: &Entity) -> OgeomResult<(Curve, (f64, f64))> {
        let n = usize::try_from(entity.at(3).int()).unwrap_or(0);
        if n == 0 {
            ogeom_bail!(Construction, "D{de}: a spline curve with no segments");
        }
        let scale = self.report.scale_mm;
        // Break points T(1..=N+1), then 12 coefficients per segment; the
        // polynomial argument runs over [0, h] within each span.
        let t = |i: usize| entity.at(4 + i).real();
        let base = 4 + n + 1;
        let mut control: Vec<Point> = Vec::with_capacity(3 * n + 1);
        let mut knots: Vec<f64> = vec![t(0); 4];
        for seg in 0..n {
            let h = t(seg + 1) - t(seg);
            if h <= 0.0 {
                ogeom_bail!(Construction, "D{de}: spline segment {seg} has no span");
            }
            let co = |k: usize| entity.at(base + 12 * seg + k).real();
            let (ax, bx, cx, dx) = (co(0), co(1), co(2), co(3));
            let (ay, by, cy, dy) = (co(4), co(5), co(6), co(7));
            let (az, bz, cz, dz) = (co(8), co(9), co(10), co(11));
            let at = |s: f64| {
                Point::new(
                    (ax + s * (bx + s * (cx + s * dx))) * scale,
                    (ay + s * (by + s * (cy + s * dy))) * scale,
                    (az + s * (bz + s * (cz + s * dz))) * scale,
                )
            };
            // Bernstein form of a cubic on [0, h]: the endpoints, and one
            // third of the end derivatives standing off them.
            let p0 = at(0.0);
            let p3 = at(h);
            let d0 = Vector::new(bx, by, bz) * (h * scale / 3.0);
            let d1 = Vector::new(
                bx + 2.0 * cx * h + 3.0 * dx * h * h,
                by + 2.0 * cy * h + 3.0 * dy * h * h,
                bz + 2.0 * cz * h + 3.0 * dz * h * h,
            ) * (h * scale / 3.0);
            if seg == 0 {
                control.push(p0);
            }
            control.push(p0 + d0);
            control.push(p3 - d1);
            control.push(p3);
            if seg + 1 < n {
                knots.extend([t(seg + 1); 3]);
            }
        }
        knots.extend([t(n); 4]);
        let curve = BSplineCurve::new(KnotVector::new(knots, 3)?, control, self.tol)?;
        let range = ogeom_geom::Curve3d::domain(&curve);
        Ok((Curve::from(curve), range))
    }

    /// Rational B-spline curve, the direct translation.
    fn nurbs_curve(&mut self, de: i64, entity: &Entity) -> OgeomResult<(Curve, (f64, f64))> {
        let k = usize::try_from(entity.at(0).int()).unwrap_or(0);
        let degree = usize::try_from(entity.at(1).int()).unwrap_or(0);
        if degree == 0 {
            ogeom_bail!(Construction, "D{de}: a B-spline curve of degree zero");
        }
        let n_ctrl = k + 1;
        let n_knots = n_ctrl + degree + 1;
        let knots: Vec<f64> = (0..n_knots).map(|i| entity.at(6 + i).real()).collect();
        let w_base = 6 + n_knots;
        let p_base = w_base + n_ctrl;
        let scale = self.report.scale_mm;
        let mut control = Vec::with_capacity(n_ctrl);
        for i in 0..n_ctrl {
            let w = entity.at(w_base + i).real();
            let p = Point::new(
                entity.at(p_base + 3 * i).real() * scale,
                entity.at(p_base + 3 * i + 1).real() * scale,
                entity.at(p_base + 3 * i + 2).real() * scale,
            );
            control.push(Weighted::new(p, w, self.tol)?);
        }
        let curve = BSplineCurve::rational(KnotVector::new(knots, degree)?, control)?;
        let v0 = entity.at(p_base + 3 * n_ctrl).real();
        let v1 = entity.at(p_base + 3 * n_ctrl + 1).real();
        let domain = ogeom_geom::Curve3d::domain(&curve);
        let range = if v1 > v0 { (v0, v1) } else { domain };
        Ok((Curve::from(curve), range))
    }

    /// A model-space surface.
    fn surface(&mut self, de: i64) -> OgeomResult<SurfaceGeometry> {
        let entity = self.entity(de)?;
        let scale = self.report.scale_mm;
        let surface: SurfaceGeometry = match entity.kind {
            108 => {
                // A x + B y + C z = D, unbounded; any face's trim bounds it.
                let normal = Vector::new(
                    entity.at(0).real(),
                    entity.at(1).real(),
                    entity.at(2).real(),
                );
                let d = entity.at(3).real();
                let origin = Point::ORIGIN + normal * (d / normal.dot(normal)) * scale;
                let plane = Plane::through(origin, Direction::new(normal, self.tol)?);
                PlaneSurface::over(
                    plane,
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                )?
                .into()
            }
            190 => {
                let point = self.location_entity(entity.at(0).int())?;
                let normal = self.direction_entity(entity.at(1).int())?;
                let plane = Plane::through(point, normal);
                PlaneSurface::over(
                    plane,
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                )?
                .into()
            }
            120 => {
                // Axis line, generatrix, start and terminate angles.
                let axis = {
                    let (line, _) = self.curve(entity.at(0).int())?;
                    let Curve::Line(l) = line else {
                        ogeom_bail!(
                            Construction,
                            "D{de}: a surface of revolution's axis is not a line"
                        );
                    };
                    l.axis()
                };
                let (curve, range) = self.curve(entity.at(1).int())?;
                let mut curve = trimmed_to(curve, range, self.tol)?;
                let sa = entity.at(2).real();
                let ta = entity.at(3).real();
                let sweep = if ta > sa {
                    ta - sa
                } else {
                    core::f64::consts::TAU
                };
                if sa.abs() > self.tol.parametric() {
                    curve = curve.transformed(&Transform::rotation(axis, sa), self.tol)?;
                }
                RevolutionSurface::new(curve, axis, sweep)?.into()
            }
            122 => {
                // Directrix, plus the far end of a generator drawn through
                // the directrix's start point.
                let (curve, range) = self.curve(entity.at(0).int())?;
                let start = curve.point_at(range.0, self.tol)?;
                let far = self.point3(entity, 1);
                let vec = far - start;
                let curve = trimmed_to(curve, range, self.tol)?;
                ExtrusionSurface::new(curve, Direction::new(vec, self.tol)?, vec.magnitude())?
                    .into()
            }
            128 => self.nurbs_surface(de, entity)?,
            118 => self.ruled_surface(de, entity)?,
            140 => {
                // Normal, distance, base surface: displaced along the
                // file's direction, which is the base's own normal or its
                // opposite — the sign carries the difference.
                let given = Direction::from_coords(
                    entity.at(0).real(),
                    entity.at(1).real(),
                    entity.at(2).real(),
                    self.tol,
                )?;
                let distance = entity.at(3).real() * scale;
                let basis = self.surface(entity.at(4).int())?;
                let ((ua, ub), (va, vb)) = ogeom_geom::Surface::domain(&basis);
                let own = ogeom_geom::Surface::normal_at(
                    &basis,
                    f64::midpoint(ua, ub),
                    f64::midpoint(va, vb),
                    self.tol,
                )?;
                let signed = if own.vector().dot(given.vector()) < 0.0 {
                    -distance
                } else {
                    distance
                };
                SurfaceGeometry::Offset(Box::new(OffsetSurface::new(basis, signed)?))
            }
            192 => {
                let point = self.location_entity(entity.at(0).int())?;
                let dir = self.direction_entity(entity.at(1).int())?;
                let radius = entity.at(2).real() * scale;
                let frame = frame_about(point, dir, self.tol)?;
                CylinderSurface::new(
                    Cylinder::new(frame, radius, self.tol)?,
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                )?
                .into()
            }
            194 => {
                let point = self.location_entity(entity.at(0).int())?;
                let dir = self.direction_entity(entity.at(1).int())?;
                let radius = entity.at(2).real() * scale;
                let half_angle = entity.at(3).real().to_radians();
                let frame = frame_about(point, dir, self.tol)?;
                ConeSurface::new(
                    Cone::new(frame, radius, half_angle, self.tol)?,
                    (-SURFACE_EXTENT, SURFACE_EXTENT),
                )?
                .into()
            }
            196 => {
                let centre = self.location_entity(entity.at(0).int())?;
                let radius = entity.at(1).real() * scale;
                SphereSurface::new(Sphere::centred(centre, radius, self.tol)?).into()
            }
            198 => {
                let centre = self.location_entity(entity.at(0).int())?;
                let dir = self.direction_entity(entity.at(1).int())?;
                let major = entity.at(2).real() * scale;
                let minor = entity.at(3).real() * scale;
                let frame = frame_about(centre, dir, self.tol)?;
                TorusSurface::new(Torus::new(frame, major, minor, self.tol)?).into()
            }
            kind => ogeom_bail!(
                Construction,
                "D{de}: surface entity type {kind}{} is not translated — \
                 docs/PARITY.md, io.iges",
                super::entity_name(kind).map_or_else(String::new, |n| format!(" ({n})"))
            ),
        };
        let placement = self.placement(entity)?;
        if placement == Transform::IDENTITY {
            Ok(surface)
        } else {
            Ok(surface.transformed(&placement, self.tol)?)
        }
    }

    /// Ruled surface (118): the straight lines between points of equal
    /// parameter on two curves, as a patch of degree one across. Both
    /// curves in their exact spline form over `[0, 1]`, raised to one
    /// degree and refined to one knot vector, the second walked backward
    /// where the file's direction flag says so. The file's form 0 joins
    /// points of equal *arc length* fraction; that is read as equal
    /// parameter, exact for the uniform-speed curves and a hair off for the
    /// rest.
    fn ruled_surface(&mut self, de: i64, entity: &Entity) -> OgeomResult<SurfaceGeometry> {
        let (first, first_range) = self.curve(entity.at(0).int())?;
        let (second, second_range) = self.curve(entity.at(1).int())?;
        let mut a = first.to_bspline_over(first_range, self.tol)?;
        let mut b = second.to_bspline_over(second_range, self.tol)?;
        if entity.at(2).int() == 1 {
            let Curve::BSpline(turned) =
                ogeom_geom::Reversible::reversed(&Curve::BSpline(b.clone()))
            else {
                ogeom_bail!(Construction, "D{de}: a reversed spline is not a spline");
            };
            b = turned;
        }
        while a.degree() < b.degree() {
            a = a.elevated(self.tol)?;
        }
        while b.degree() < a.degree() {
            b = b.elevated(self.tol)?;
        }
        let tol = self.tol;
        let unify = |target: &mut BSplineCurve, source: &BSplineCurve| -> OgeomResult<()> {
            let (lo, hi) = source.knots().domain();
            for (value, multiplicity) in source.knots().distinct() {
                if value <= lo + 1e-12 || value >= hi - 1e-12 {
                    continue;
                }
                let held = target.knots().multiplicity_of(value);
                if multiplicity > held {
                    *target = target.with_knot_inserted(value, multiplicity - held, tol)?;
                }
            }
            Ok(())
        };
        unify(&mut a, &b)?;
        unify(&mut b, &a)?;
        let nu = a.control_points().len();
        if b.control_points().len() != nu {
            ogeom_bail!(
                Construction,
                "D{de}: the ruled surface's curves refine to different nets"
            );
        }
        let mut points: Vec<Weighted<Point>> = Vec::with_capacity(nu * 2);
        for i in 0..nu {
            points.push(a.control_points()[i]);
            points.push(b.control_points()[i]);
        }
        let grid = ControlGrid::new(points, nu, 2)?;
        let across = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
        Ok(BSplineSurface::rational(a.knots().clone(), across, grid)?.into())
    }

    fn nurbs_surface(&mut self, de: i64, entity: &Entity) -> OgeomResult<SurfaceGeometry> {
        let k1 = usize::try_from(entity.at(0).int()).unwrap_or(0);
        let k2 = usize::try_from(entity.at(1).int()).unwrap_or(0);
        let m1 = usize::try_from(entity.at(2).int()).unwrap_or(0);
        let m2 = usize::try_from(entity.at(3).int()).unwrap_or(0);
        if m1 == 0 || m2 == 0 {
            ogeom_bail!(Construction, "D{de}: a B-spline surface of degree zero");
        }
        let (nu, nv) = (k1 + 1, k2 + 1);
        let (nku, nkv) = (nu + m1 + 1, nv + m2 + 1);
        let base = 9;
        let u_knots: Vec<f64> = (0..nku).map(|i| entity.at(base + i).real()).collect();
        let v_knots: Vec<f64> = (0..nkv).map(|i| entity.at(base + nku + i).real()).collect();
        let w_base = base + nku + nkv;
        let p_base = w_base + nu * nv;
        let scale = self.report.scale_mm;
        // The file lists control points with the first (u) index varying
        // fastest; the grid stores row-major with u as the slow index, so
        // the read reorders.
        let raw = |u: usize, v: usize| -> OgeomResult<Weighted<Point>> {
            let i = v * nu + u;
            let w = entity.at(w_base + i).real();
            let p = Point::new(
                entity.at(p_base + 3 * i).real() * scale,
                entity.at(p_base + 3 * i + 1).real() * scale,
                entity.at(p_base + 3 * i + 2).real() * scale,
            );
            Weighted::new(p, w, self.tol)
        };
        let mut weighted = Vec::with_capacity(nu * nv);
        for u in 0..nu {
            for v in 0..nv {
                weighted.push(raw(u, v)?);
            }
        }
        let grid = ControlGrid::new(weighted, nu, nv)?;
        Ok(ogeom_geom::BSplineSurface::rational(
            KnotVector::new(u_knots, m1)?,
            KnotVector::new(v_knots, m2)?,
            grid,
        )?
        .into())
    }

    /// A point entity (116) or bare coordinate triple carrier.
    fn location_entity(&mut self, de: i64) -> OgeomResult<Point> {
        let e = self.entity(de)?;
        if e.kind != 116 {
            ogeom_bail!(
                Construction,
                "D{de}: expected a point entity, found type {}",
                e.kind
            );
        }
        Ok(self.point3(e, 0))
    }

    /// A direction entity (123) as a unit vector.
    fn direction_entity(&mut self, de: i64) -> OgeomResult<Direction> {
        let e = self.entity(de)?;
        if e.kind != 123 {
            ogeom_bail!(
                Construction,
                "D{de}: expected a direction entity, found type {}",
                e.kind
            );
        }
        Direction::from_coords(e.at(0).real(), e.at(1).real(), e.at(2).real(), self.tol)
    }

    /// The boundary of a trimmed face as model-space curve segments — from a
    /// curve-on-surface (142), a boundary entity (141), or a bare curve,
    /// walking composite curves flat.
    fn boundary_segments(&mut self, de: i64) -> OgeomResult<Vec<(Curve, (f64, f64))>> {
        let entity = self.entity(de)?;
        match entity.kind {
            // Curve on surface: the model-space curve is the trim's truth;
            // the file's pcurve is advisory, because faces recompute exact
            // pcurves and say when they cannot.
            142 => {
                let c = entity.at(2).int();
                if c == 0 {
                    ogeom_bail!(
                        Construction,
                        "D{de}: a curve-on-surface carries no model-space \
                         curve; pcurve-only trimming is not translated — \
                         docs/PARITY.md, io.iges"
                    );
                }
                self.curve_segments(c)
            }
            141 => {
                let n = usize::try_from(entity.at(3).int()).unwrap_or(0);
                let mut out = Vec::new();
                let mut i = 4;
                for _ in 0..n {
                    let cptr = entity.at(i).int();
                    let k = usize::try_from(entity.at(i + 2).int()).unwrap_or(0);
                    i += 3 + k;
                    // The sense flag is advisory here too: the wire builder
                    // chains segments by their geometry.
                    out.extend(self.curve_segments(cptr)?);
                }
                Ok(out)
            }
            _ => self.curve_segments(de),
        }
    }

    fn curve_segments(&mut self, de: i64) -> OgeomResult<Vec<(Curve, (f64, f64))>> {
        let entity = self.entity(de)?;
        if entity.kind == 102 {
            let n = usize::try_from(entity.at(0).int()).unwrap_or(0);
            let mut out = Vec::new();
            for i in 0..n {
                out.extend(self.curve_segments(entity.at(1 + i).int())?);
            }
            return Ok(out);
        }
        Ok(vec![self.curve(de)?])
    }

    /// A trimmed (144) or bounded (143) surface as a face.
    fn face(&mut self, de: i64) -> OgeomResult<Shape> {
        let entity = self.entity(de)?;
        let (surface_de, boundaries): (i64, Vec<i64>) = match entity.kind {
            144 => {
                let s = entity.at(0).int();
                let outer_given = entity.at(1).int() == 1;
                let n_inner = usize::try_from(entity.at(2).int()).unwrap_or(0);
                let mut bs = Vec::new();
                if outer_given && entity.at(3).int() != 0 {
                    bs.push(entity.at(3).int());
                }
                for i in 0..n_inner {
                    bs.push(entity.at(4 + i).int());
                }
                (s, bs)
            }
            143 => {
                let s = entity.at(1).int();
                let n = usize::try_from(entity.at(2).int()).unwrap_or(0);
                (s, (0..n).map(|i| entity.at(3 + i).int()).collect())
            }
            kind => ogeom_bail!(Construction, "D{de}: type {kind} is not a trimmed surface"),
        };
        let surface = self.surface(surface_de)?;
        if boundaries.is_empty() {
            // The surface's own natural boundary: a closed or bounded
            // surface can stand alone as a face.
            return Ok(ogeom_algo::make_natural_face(&mut self.model, surface)?.shape);
        }
        let mut wires = Vec::new();
        for boundary in boundaries {
            let segments = self.boundary_segments(boundary)?;
            wires.push(self.wire_edges(de, segments)?);
        }
        self.assemble_face(surface, wires)
    }

    /// Boundary segments into a closed chain of edges, head to tail.
    ///
    /// Surface files are loose about sense — a boundary's segments arrive in
    /// order but each may run either way — so the chain is stitched by
    /// geometry: each segment joins whichever of its ends sits at the chain's
    /// current head, and the last vertex is the first, which is what closes
    /// the wire.
    fn wire_edges(
        &mut self,
        face_de: i64,
        segments: Vec<(Curve, (f64, f64))>,
    ) -> OgeomResult<Vec<Shape>> {
        if segments.is_empty() {
            ogeom_bail!(Construction, "D{face_de}: a boundary with no curves");
        }
        let ends: Vec<(Point, Point)> = segments
            .iter()
            .map(|(c, r)| Ok((c.point_at(r.0, self.tol)?, c.point_at(r.1, self.tol)?)))
            .collect::<OgeomResult<_>>()?;
        let n = segments.len();
        let weld = self.tol.confusion() * 100.0;

        // One closed segment closes on a single vertex.
        if n == 1 {
            let (curve, range) = segments
                .into_iter()
                .next()
                .unwrap_or_else(|| unreachable!());
            let (s, e) = ends[0];
            if s.distance(e) > weld {
                ogeom_bail!(
                    Construction,
                    "D{face_de}: a one-curve boundary whose ends sit {:.2e} apart",
                    s.distance(e)
                );
            }
            let v = make_vertex(&mut self.model, s).shape;
            let edge = make_edge_between(&mut self.model, curve, range, &v, &v, self.tol)?.shape;
            return Ok(vec![edge]);
        }

        let head = make_vertex(&mut self.model, ends[0].0).shape;
        let mut at = ends[0].0;
        let mut at_vertex = head.clone();
        let mut edges = Vec::with_capacity(n);
        for (i, ((curve, range), (s, e))) in segments.into_iter().zip(ends).enumerate() {
            let forward = at.distance(s) <= at.distance(e);
            let (this_end, this_point) = if forward { (e, e) } else { (s, s) };
            let gap = at.distance(if forward { s } else { e });
            if gap > weld {
                ogeom_bail!(
                    Construction,
                    "D{face_de}: boundary segment {i} starts {gap:.2e} from \
                     where the previous one ended"
                );
            }
            let last = i + 1 == n;
            let next_vertex = if last {
                head.clone()
            } else {
                make_vertex(&mut self.model, this_end).shape
            };
            // The edge is built along its curve's own direction; a segment
            // running against the chain is used reversed, exactly as a
            // hand-built prism's rim edges are.
            let edge = if forward {
                make_edge_between(
                    &mut self.model,
                    curve,
                    range,
                    &at_vertex,
                    &next_vertex,
                    self.tol,
                )?
                .shape
            } else {
                make_edge_between(
                    &mut self.model,
                    curve,
                    range,
                    &next_vertex,
                    &at_vertex,
                    self.tol,
                )?
                .shape
                .reversed()
            };
            edges.push(edge);
            at = this_point;
            at_vertex = next_vertex;
        }
        Ok(edges)
    }

    /// A face from a surface and wires of already-built edges: the exact
    /// pcurve where the pair has a closed form, the fitted one where it does
    /// not, and both sides of the chart for an edge the wire uses twice —
    /// which is what a seam is.
    fn assemble_face(
        &mut self,
        surface: SurfaceGeometry,
        wires: Vec<Vec<Shape>>,
    ) -> OgeomResult<Shape> {
        // A single wire that is one edge used twice on a periodic surface is
        // a seam and nothing else, and a boundary that is nothing but the
        // seam encloses the whole chart: the face is the surface, and the
        // natural face carries its own degenerate boundary — a sphere's
        // poles — which the file had no edges for.
        if let [edges] = wires.as_slice()
            && let [a, b] = edges.as_slice()
            && a.node() == b.node()
        {
            use ogeom_geom::Surface as _;
            if surface.is_periodic_u() || surface.is_periodic_v() {
                return Ok(ogeom_algo::make_natural_face(&mut self.model, surface)?.shape);
            }
        }
        let surface_id = self.model.geometry_mut().add_surface(surface.clone());
        let mut wire_shapes = Vec::with_capacity(wires.len());
        for edges in &wires {
            wire_shapes.push(ogeom_algo::make_wire(&mut self.model, edges, self.tol)?.shape);
        }
        let face =
            ogeom_algo::make_face_on(&mut self.model, surface_id, &wire_shapes, self.tol)?.shape;

        for edges in &wires {
            self.chart_wire(edges, &surface, surface_id)?;
        }
        Ok(face)
    }

    /// One wire's pcurves, chained around the face's chart.
    ///
    /// Each edge's image is computed on its own — the exact projection where
    /// the pair has a closed form, the fitted one where it does not — and a
    /// periodic chart then has a branch to choose. Read one edge at a time
    /// the choice is arbitrary, and the wire comes apart: a bore's wall
    /// arrives with one rim written over `[−π, π]` and the other over
    /// `[−π/2, 3π/2]`, both the right circle and neither meeting the seam
    /// the wire closes on. So each image after the first is shifted by whole
    /// periods until its start meets where the last one ended, which is the
    /// rule `make_face_with_pcurves` follows for shapes this kernel builds
    /// itself.
    ///
    /// A seam falls out of the same walk. The wire uses it twice, up one
    /// column of the chart and down the other, and those columns *are* one
    /// image a period apart — so the chaining produces both, and the use
    /// that runs forward is the forward side.
    fn chart_wire(
        &mut self,
        edges: &[Shape],
        surface: &SurfaceGeometry,
        surface_id: ogeom_topo::SurfaceId,
    ) -> OgeomResult<()> {
        use ogeom_geom::Curve2d as _;
        let mut counts: HashMap<ogeom_topo::TShapeId, usize> = HashMap::new();
        for edge in edges {
            *counts.entry(edge.node()).or_default() += 1;
        }
        // One image per edge, however many times the wire walks it.
        let mut images: HashMap<ogeom_topo::TShapeId, (ogeom_geom::PlanarCurve, (f64, f64))> =
            HashMap::new();
        // The sides each walk of it left, by the direction that walk ran.
        let mut sides: HashMap<
            ogeom_topo::TShapeId,
            (
                Option<ogeom_geom::PlanarCurve>,
                Option<ogeom_geom::PlanarCurve>,
            ),
        > = HashMap::new();
        let mut order: Vec<ogeom_topo::TShapeId> = Vec::new();
        let mut previous: Option<ogeom_math::Point2> = None;
        for edge in edges {
            let (curve, range) = {
                let Some(data) = self.model.node(edge).and_then(|n| n.data().as_edge()) else {
                    continue;
                };
                let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d()
                else {
                    continue;
                };
                let Some(geometry) = self.model.geometry().curve(*curve) else {
                    continue;
                };
                (geometry.clone(), *range)
            };
            if let std::collections::hash_map::Entry::Vacant(slot) = images.entry(edge.node()) {
                let Some(image) = self.image_of(edge, &curve, range, surface)? else {
                    continue;
                };
                slot.insert((image, range));
                order.push(edge.node());
            }
            let Some((image, range)) = images.get(&edge.node()).cloned() else {
                continue;
            };
            let backwards = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let (start, end) = if backwards {
                (range.1, range.0)
            } else {
                (range.0, range.1)
            };
            let image =
                crate::pcurves::shifted_to_meet(&image, start, previous, surface, self.tol)?;
            previous = Some(image.point_at(end, self.tol)?);
            let walked = sides.entry(edge.node()).or_default();
            if backwards {
                walked.1 = Some(image);
            } else {
                walked.0 = Some(image);
            }
        }
        for node in order {
            let Some((image, range)) = images.get(&node).cloned() else {
                continue;
            };
            let Some(edge) = edges.iter().find(|e| e.node() == node) else {
                continue;
            };
            let walked = sides.remove(&node).unwrap_or((None, None));
            let seam = counts.get(&node).copied().unwrap_or(0) > 1;
            let columns = match &walked {
                (Some(forward), Some(reversed)) => {
                    let at = forward.point_at(range.0, self.tol)?;
                    let other = reversed.point_at(range.0, self.tol)?;
                    (at.distance(other) > self.tol.confusion()).then_some((forward, reversed))
                }
                _ => None,
            };
            match (seam, columns) {
                (true, Some((forward, reversed))) => {
                    ogeom_algo::attach_seam(
                        &mut self.model,
                        edge,
                        forward.clone(),
                        reversed.clone(),
                        surface_id,
                        ogeom_topo::Location::identity(),
                        range,
                    )?;
                }
                // The walk left the seam's two uses in one place: the chart
                // closes without being periodic — a skinned wall's is such a
                // chart, clamped and closed — and there is no period to
                // shift by. The other column goes a chart's width over,
                // toward the middle, which is where it went before there
                // was a walk to ask.
                (true, None) => {
                    let other = crate::pcurves::seam_other_side(&image, range, surface, self.tol)?;
                    ogeom_algo::attach_seam(
                        &mut self.model,
                        edge,
                        image,
                        other,
                        surface_id,
                        ogeom_topo::Location::identity(),
                        range,
                    )?;
                }
                (false, _) => {
                    let (forward, reversed) = walked;
                    ogeom_algo::attach_pcurve(
                        &mut self.model,
                        edge,
                        forward.or(reversed).unwrap_or(image),
                        surface_id,
                        ogeom_topo::Location::identity(),
                        range,
                    )?;
                }
            }
        }
        Ok(())
    }

    /// One edge's image on one surface, exact where the pair has a closed
    /// form and fitted where it does not — the same policy the STEP reader
    /// applies, through the shared machinery. Where the wire puts it on a
    /// periodic chart is [`Self::chart_wire`]'s business.
    fn image_of(
        &mut self,
        edge: &Shape,
        curve: &Curve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
    ) -> OgeomResult<Option<ogeom_geom::PlanarCurve>> {
        use ogeom_geom::PlanarCurve;
        let widen = |p: PlanarCurve| -> PlanarCurve {
            if let PlanarCurve::Line(l) = &p {
                use ogeom_geom::Curve2d as _;
                let (lo, hi) = (l.domain().0.min(range.0), l.domain().1.max(range.1));
                if let Ok(wider) = ogeom_geom::Line2d::over(l.axis(), lo, hi) {
                    return wider.into();
                }
            }
            p
        };
        let found = match ogeom_intersect::exact_pcurve_over(curve, range, surface, self.tol)
            .map(widen)
        {
            Some(exact) => exact,
            None => match crate::pcurves::fit_projected_pcurve(curve, range, surface, self.tol) {
                Ok((fitted, error, met, worst_off, slop)) => {
                    if let Some(w) = slop {
                        self.report.warnings.push(w);
                    }
                    if !met {
                        self.report.warnings.push(format!(
                            "a projected pcurve fit stopped at {error:.2e}; \
                             the face's mesh may sit that far off along this edge"
                        ));
                    }
                    if worst_off > self.tol.confusion()
                        && let Some(node) = self.model.node_mut(edge)
                        && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                    {
                        data.tolerance = data.tolerance.widen_to(worst_off + self.tol.confusion());
                    }
                    fitted
                }
                Err(e) => {
                    self.report.warnings.push(format!(
                        "no pcurve for an edge on this surface ({e}); the \
                         face may not triangulate"
                    ));
                    return Ok(None);
                }
            },
        };
        Ok(Some(found))
    }

    /// A manifold solid B-rep object: shell of faces of loops of edges.
    fn manifold_solid(&mut self, de: i64) -> OgeomResult<Shape> {
        let entity = self.entity(de)?;
        let shell = self.shell(entity.at(0).int())?;
        let n_voids = usize::try_from(entity.at(2).int()).unwrap_or(0);
        let mut shells = vec![shell];
        for i in 0..n_voids {
            shells.push(self.shell(entity.at(3 + 2 * i).int())?);
        }
        let solid = make_solid(&mut self.model, &shells)?.shape;
        // A fitted trim widens its edge to the offset it measured; the
        // edge's vertices come along, and the solid keeps the containment
        // rule the checker holds it to.
        ogeom_algo::restore_containment(&mut self.model, &solid)?;
        Ok(solid)
    }

    fn shell(&mut self, de: i64) -> OgeomResult<Shape> {
        let entity = self.entity(de)?;
        if entity.kind != 514 {
            ogeom_bail!(
                Construction,
                "D{de}: expected a shell, found type {}",
                entity.kind
            );
        }
        let n = usize::try_from(entity.at(0).int()).unwrap_or(0);
        let mut faces = Vec::with_capacity(n);
        for i in 0..n {
            let face_de = entity.at(1 + 2 * i).int();
            let same_sense = entity.at(2 + 2 * i).int() != 0;
            let face = self.brep_face(face_de)?;
            faces.push(if same_sense { face } else { face.reversed() });
        }
        Ok(ogeom_algo::make_shell(&mut self.model, &faces)?.shape)
    }

    fn brep_face(&mut self, de: i64) -> OgeomResult<Shape> {
        let entity = self.entity(de)?;
        if entity.kind != 510 {
            ogeom_bail!(
                Construction,
                "D{de}: expected a face, found type {}",
                entity.kind
            );
        }
        let surface = self.surface(entity.at(0).int())?;
        let n_loops = usize::try_from(entity.at(1).int()).unwrap_or(0);
        // Parameter 2 is the outer-loop flag; the loop pointers follow.
        let mut wires = Vec::with_capacity(n_loops);
        for i in 0..n_loops {
            wires.push(self.loop_edges(entity.at(3 + i).int())?);
        }
        self.assemble_face(surface, wires)
    }

    fn loop_edges(&mut self, de: i64) -> OgeomResult<Vec<Shape>> {
        let entity = self.entity(de)?;
        if entity.kind != 508 {
            ogeom_bail!(
                Construction,
                "D{de}: expected a loop, found type {}",
                entity.kind
            );
        }
        let n = usize::try_from(entity.at(0).int()).unwrap_or(0);
        let mut edges = Vec::with_capacity(n);
        let mut i = 1;
        for _ in 0..n {
            let is_vertex = entity.at(i).int() == 1;
            let list_de = entity.at(i + 1).int();
            let index = entity.at(i + 2).int();
            let orientation = entity.at(i + 3).int();
            let k = usize::try_from(entity.at(i + 4).int()).unwrap_or(0);
            i += 5 + 2 * k;
            if is_vertex {
                // A vertex entry marks a degenerate use; the face builder
                // rebuilds chart degeneracies from the surface itself.
                self.report.warnings.push(format!(
                    "D{de}: a loop lists a vertex entry, which this reader skips"
                ));
                continue;
            }
            let (edge, _, _) = self.list_edge(list_de, index)?;
            edges.push(if orientation != 0 {
                edge
            } else {
                edge.reversed()
            });
        }
        Ok(edges)
    }

    /// Edge `index` (1-based) of an edge list (504), built once and shared.
    fn list_edge(&mut self, list_de: i64, index: i64) -> OgeomResult<BuiltEdge> {
        let key = (list_de, index);
        if let Some(found) = self.edges.get(&key) {
            return Ok(found.clone());
        }
        let entity = self.entity(list_de)?;
        if entity.kind != 504 {
            ogeom_bail!(
                Construction,
                "D{list_de}: expected an edge list, found type {}",
                entity.kind
            );
        }
        let i = usize::try_from(index - 1).map_err(|_| {
            ogeom_core::ogeom_err!(Construction, "D{list_de}: edge index {index} out of range")
        })?;
        let base = 1 + 5 * i;
        let curve_de = entity.at(base).int();
        let (sv_list, sv_index) = (entity.at(base + 1).int(), entity.at(base + 2).int());
        let (tv_list, tv_index) = (entity.at(base + 3).int(), entity.at(base + 4).int());
        let (curve, mut range) = self.curve(curve_de)?;
        let vs = self.list_vertex(sv_list, sv_index)?;
        let ve = self.list_vertex(tv_list, tv_index)?;
        let (ps, pe) = (self.point_of(&vs), self.point_of(&ve));
        // Re-derive the range from the vertices on this kernel's own
        // parameterization, exactly as the STEP reader does and for the same
        // reason: the file's parameterization is its own business.
        if let (Some(a), Some(b)) = (parameter_on(&curve, ps), parameter_on(&curve, pe)) {
            let period = if curve.is_periodic() {
                let (lo, hi) = curve.domain();
                hi - lo
            } else {
                0.0
            };
            range = if ps.distance(pe) < self.tol.confusion() && period > 0.0 {
                (a, a + period)
            } else if period > 0.0 && b <= a + self.tol.parametric() {
                (a, b + period)
            } else {
                (a, b)
            };
            // A vertex a few nanometres past a bounded curve's end projects
            // past it: the range is held to the curve's own domain, and the
            // vertex widens below to cover the rest of the miss.
            if period == 0.0 {
                let (lo, hi) = curve.domain();
                range = (range.0.clamp(lo, hi), range.1.clamp(lo, hi));
            }
        }
        // IGES states no tolerances, so a vertex is built at the confusion
        // tolerance and a curve end that misses it by rounding — a writer's
        // last digit, a few tenths of a nanometre — would refuse the whole
        // solid. As the STEP reader does, the vertex's tolerance grows to
        // state the miss, up to the millimetre past which a boundary is not
        // this curve's at all; beyond that the edge still refuses by name.
        let cap = self.tol.confusion() * 1e7;
        for (vertex, t) in [(&vs, range.0), (&ve, range.1)] {
            let (Ok(end), Some(stated)) = (
                curve.point_at(t, self.tol),
                self.model
                    .node(vertex)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| (d.point, d.tolerance)),
            ) else {
                continue;
            };
            let gap = end.distance(stated.0);
            if gap > stated.1.get() && gap <= cap {
                self.model.widen(
                    vertex,
                    ogeom_core::Tolerance::new(gap + self.tol.confusion())?,
                )?;
                self.vertex_misses.0 += 1;
                self.vertex_misses.1 = self.vertex_misses.1.max(gap);
            }
        }
        let edge =
            make_edge_between(&mut self.model, curve.clone(), range, &vs, &ve, self.tol)?.shape;
        let built = (edge, curve, range);
        self.edges.insert(key, built.clone());
        Ok(built)
    }

    /// Vertex `index` (1-based) of a vertex list (502), built once and
    /// shared — sharing is what lets a closed shell close.
    fn list_vertex(&mut self, list_de: i64, index: i64) -> OgeomResult<Shape> {
        let key = (list_de, index);
        if let Some(found) = self.vertices.get(&key) {
            return Ok(found.clone());
        }
        let entity = self.entity(list_de)?;
        if entity.kind != 502 {
            ogeom_bail!(
                Construction,
                "D{list_de}: expected a vertex list, found type {}",
                entity.kind
            );
        }
        let i = usize::try_from(index - 1).map_err(|_| {
            ogeom_core::ogeom_err!(
                Construction,
                "D{list_de}: vertex index {index} out of range"
            )
        })?;
        let point = self.point3(entity, 1 + 3 * i);
        let vertex = make_vertex(&mut self.model, point).shape;
        self.vertices.insert(key, vertex.clone());
        Ok(vertex)
    }

    fn point_of(&self, vertex: &Shape) -> Point {
        self.model
            .node(vertex)
            .and_then(|n| n.data().as_vertex())
            .map_or(Point::ORIGIN, |d| d.point)
    }

    /// The document: parts named from entity labels, colours from the fixed
    /// palette and from 314 entities the directory colour fields point at.
    fn document(&mut self, solids: &[(i64, Shape)], sheets: &[Shape]) -> ogeom_doc::Document {
        let colours: Vec<(Shape, ogeom_doc::Colour)> = solids
            .iter()
            .filter_map(|(de, s)| self.colour_of(*de).map(|c| (s.clone(), c)))
            .collect();
        let mut document = ogeom_doc::Document::over(std::mem::take(&mut self.model));
        for (i, (de, solid)) in solids.iter().enumerate() {
            let label = self
                .file
                .entity(*de)
                .map(|e| e.label.clone())
                .filter(|l| !l.is_empty())
                .unwrap_or_else(|| format!("solid-{i}"));
            document.add_part(label, solid.clone());
        }
        for (i, sheet) in sheets.iter().enumerate() {
            document.add_part(format!("sheet-{i}"), sheet.clone());
        }
        for (shape, colour) in colours {
            document.set_colour(&shape, colour);
        }
        document
    }

    /// The colour a directory entry states: a negated field points at a 314
    /// entity's RGB percentages; a small positive names the fixed palette.
    fn colour_of(&mut self, de: i64) -> Option<ogeom_doc::Colour> {
        let entity = self.file.entity(de)?;
        let c = entity.colour;
        if c < 0 {
            let e = self.file.entity(-c)?;
            self.visited.insert(-c, ());
            if e.kind != 314 {
                return None;
            }
            return Some(ogeom_doc::Colour::rgb(
                (e.at(0).real() / 100.0).clamp(0.0, 1.0),
                (e.at(1).real() / 100.0).clamp(0.0, 1.0),
                (e.at(2).real() / 100.0).clamp(0.0, 1.0),
            ));
        }
        // 1 black, 2 red, 3 green, 4 blue, 5 yellow, 6 magenta, 7 cyan,
        // 8 white — the specification's own palette.
        let palette = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 1.0, 0.0),
            (1.0, 0.0, 1.0),
            (0.0, 1.0, 1.0),
            (1.0, 1.0, 1.0),
        ];
        let index = usize::try_from(c).ok()?.checked_sub(1)?;
        palette
            .get(index)
            .map(|&(r, g, b)| ogeom_doc::Colour::rgb(r, g, b))
    }
}

/// A trimmed carrier where the range is a strict part of the domain — a
/// generatrix used by a sweep is exactly its stated span.
fn trimmed_to(curve: Curve, range: (f64, f64), tol: Tolerances) -> OgeomResult<Curve> {
    let (lo, hi) = curve.domain();
    if (range.0 - lo).abs() < tol.parametric() && (range.1 - hi).abs() < tol.parametric() {
        return Ok(curve);
    }
    Ok(Curve::from(TrimmedCurve::new(
        curve, range.0, range.1, tol,
    )?))
}

use crate::inversion::parameter_on;

/// A frame with the given axis direction, reference direction chosen stably.
fn frame_about(origin: Point, axis: Direction, tol: Tolerances) -> OgeomResult<Frame> {
    let seed = if axis.vector().dot(Vector::X).abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = Direction::from_cross(axis.vector(), seed, tol)?;
    Frame::new(origin, axis, x, tol)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::parse::Value;
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    fn entity(kind: i64, form: i64, params: Vec<Value>) -> Entity {
        Entity {
            kind,
            form,
            transform: 0,
            colour: 0,
            level: 0,
            status: 0,
            label: String::new(),
            params,
        }
    }

    fn file(entities: Vec<(i64, Entity)>) -> File {
        let mut global = vec![Value::Default; 16];
        global[13] = Value::Int(2);
        File {
            global,
            entities: entities.into_iter().collect(),
        }
    }

    fn reader(file: &File) -> Reader<'_> {
        Reader {
            file,
            model: Model::new(),
            report: IgesReport {
                scale_mm: 1.0,
                ..IgesReport::default()
            },
            visited: BTreeMap::new(),
            vertices: HashMap::new(),
            edges: HashMap::new(),
            vertex_misses: (0, 0.0),
            tol: T,
        }
    }

    fn reals(values: &[f64]) -> Vec<Value> {
        values.iter().map(|v| Value::Real(*v)).collect()
    }

    fn line(from: [f64; 3], to: [f64; 3]) -> Entity {
        entity(
            110,
            0,
            reals(&[from[0], from[1], from[2], to[0], to[1], to[2]]),
        )
    }

    /// A conic arc whose squares differ in sign is a hyperbola, and one
    /// missing a square is a parabola; each reads as its own curve, the arc
    /// running from the file's start to its terminate point whichever way
    /// round the file lists them.
    #[test]
    fn hyperbola_and_parabola_conic_arcs_read_as_their_curves() {
        let (cosh1, sinh1) = (1.0_f64.cosh(), 1.0_f64.sinh());
        // x²/4 − y²/9 = 1, the arc from the vertex to parameter one.
        let hyperbola = |forward: bool| {
            let (s, e) = if forward {
                ([2.0, 0.0], [2.0 * cosh1, 3.0 * sinh1])
            } else {
                ([2.0 * cosh1, 3.0 * sinh1], [2.0, 0.0])
            };
            entity(
                104,
                2,
                reals(&[
                    0.25,
                    0.0,
                    -1.0 / 9.0,
                    0.0,
                    0.0,
                    -1.0,
                    0.0,
                    s[0],
                    s[1],
                    e[0],
                    e[1],
                ]),
            )
        };
        // x² − 4 y = 0 opening along y, and y² − 4 x = 0 opening along x.
        let parabola_y = entity(
            104,
            3,
            reals(&[1.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0]),
        );
        let parabola_x = entity(
            104,
            3,
            reals(&[0.0, 0.0, 1.0, -4.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0]),
        );
        let deck = file(vec![
            (1, hyperbola(true)),
            (3, hyperbola(false)),
            (5, parabola_y),
            (7, parabola_x),
        ]);
        let mut reader = reader(&deck);
        for de in [1, 3] {
            let (curve, (t0, t1)) = reader.curve(de).unwrap();
            assert!(matches!(curve, Curve::Hyperbola(_)), "D{de} is a hyperbola");
            assert!(t1 > t0, "the arc runs forward: {t0} .. {t1}");
            let entity = deck.entity(de).unwrap();
            let start = Point::new(entity.at(7).real(), entity.at(8).real(), 0.0);
            let end = Point::new(entity.at(9).real(), entity.at(10).real(), 0.0);
            assert!(curve.point_at(t0, T).unwrap().distance(start) < 1e-9);
            assert!(curve.point_at(t1, T).unwrap().distance(end) < 1e-9);
            let mid = curve.point_at(f64::midpoint(t0, t1), T).unwrap();
            assert!((mid.x * mid.x / 4.0 - mid.y * mid.y / 9.0 - 1.0).abs() < 1e-9);
        }
        for (de, along_y) in [(5, true), (7, false)] {
            let (curve, (t0, t1)) = reader.curve(de).unwrap();
            assert!(matches!(curve, Curve::Parabola(_)), "D{de} is a parabola");
            assert!(t1 > t0);
            let (start, end) = if along_y {
                (Point::ORIGIN, Point::new(4.0, 4.0, 0.0))
            } else {
                (Point::new(4.0, 4.0, 0.0), Point::ORIGIN)
            };
            assert!(curve.point_at(t0, T).unwrap().distance(start) < 1e-9);
            assert!(curve.point_at(t1, T).unwrap().distance(end) < 1e-9);
            let mid = curve.point_at(f64::midpoint(t0, t1), T).unwrap();
            let residual = if along_y {
                mid.x * mid.x - 4.0 * mid.y
            } else {
                mid.y * mid.y - 4.0 * mid.x
            };
            assert!(residual.abs() < 1e-9, "on the parabola: {mid:?}");
        }
    }

    /// A ruled surface between two lines is the bilinear patch between
    /// them, the second line walked backward where the direction flag says.
    #[test]
    fn a_ruled_surface_reads_as_the_patch_between_its_curves() {
        let deck = file(vec![
            (1, line([0.0, 0.0, 0.0], [10.0, 0.0, 0.0])),
            (3, line([0.0, 5.0, 2.0], [10.0, 5.0, 2.0])),
            (
                5,
                entity(
                    118,
                    1,
                    vec![Value::Int(1), Value::Int(3), Value::Int(0), Value::Int(0)],
                ),
            ),
            (
                7,
                entity(
                    118,
                    1,
                    vec![Value::Int(1), Value::Int(3), Value::Int(1), Value::Int(0)],
                ),
            ),
        ]);
        let mut reader = reader(&deck);
        let straight = reader.surface(5).unwrap();
        let turned = reader.surface(7).unwrap();
        use ogeom_geom::Surface as _;
        let middle = straight.point_at(0.5, 0.5, T).unwrap();
        assert!(
            middle.distance(Point::new(5.0, 2.5, 1.0)) < 1e-9,
            "{middle:?}"
        );
        let far = straight.point_at(0.0, 1.0, T).unwrap();
        assert!(far.distance(Point::new(0.0, 5.0, 2.0)) < 1e-9, "{far:?}");
        let far = turned.point_at(0.0, 1.0, T).unwrap();
        assert!(far.distance(Point::new(10.0, 5.0, 2.0)) < 1e-9, "{far:?}");
    }

    /// An offset surface displaces its base along the file's direction,
    /// and an offset curve along the file's reference crossed with the
    /// tangent — whichever way round this vocabulary spells either.
    #[test]
    fn offset_entities_displace_the_way_the_file_says() {
        let deck = file(vec![
            (1, entity(108, 0, reals(&[0.0, 0.0, 1.0, 0.0]))),
            (
                3,
                entity(
                    140,
                    0,
                    vec![
                        Value::Real(0.0),
                        Value::Real(0.0),
                        Value::Real(1.0),
                        Value::Real(3.0),
                        Value::Int(1),
                    ],
                ),
            ),
            (
                5,
                entity(
                    140,
                    0,
                    vec![
                        Value::Real(0.0),
                        Value::Real(0.0),
                        Value::Real(-1.0),
                        Value::Real(3.0),
                        Value::Int(1),
                    ],
                ),
            ),
            (7, line([0.0, 0.0, 0.0], [10.0, 0.0, 0.0])),
            (
                9,
                entity(
                    130,
                    0,
                    vec![
                        Value::Int(7),
                        Value::Int(1),
                        Value::Int(0),
                        Value::Int(3),
                        Value::Int(0),
                        Value::Real(2.0),
                        Value::Real(0.0),
                        Value::Real(2.0),
                        Value::Real(0.0),
                        Value::Real(0.0),
                        Value::Real(0.0),
                        Value::Real(1.0),
                        Value::Real(0.0),
                        Value::Real(10.0),
                    ],
                ),
            ),
        ]);
        let mut reader = reader(&deck);
        use ogeom_geom::Surface as _;
        let up = reader.surface(3).unwrap();
        assert!((up.point_at(1.0, 2.0, T).unwrap().z - 3.0).abs() < 1e-9);
        let down = reader.surface(5).unwrap();
        assert!((down.point_at(1.0, 2.0, T).unwrap().z + 3.0).abs() < 1e-9);
        let (offset, (t0, t1)) = reader.curve(9).unwrap();
        assert!(matches!(offset, Curve::Offset(_)));
        assert!((t0, t1) == (0.0, 10.0));
        let at = offset.point_at(5.0, T).unwrap();
        assert!(at.distance(Point::new(5.0, 2.0, 0.0)) < 1e-9, "{at:?}");
    }
}
