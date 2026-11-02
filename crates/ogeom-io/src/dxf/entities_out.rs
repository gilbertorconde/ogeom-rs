//! The R2000 writer: a drawing's entities as written, with its unit and
//! layers.
//!
//! An R2000 (AC1015) file is read by its handles: every table, table entry,
//! block, entity and object carries one, and each names its owner. The
//! file holds what a strict reader asks for and nothing besides:
//!
//! - HEADER: `$ACADVER`, `$DWGCODEPAGE`, `$INSUNITS`, `$MEASUREMENT` and
//!   `$HANDSEED`, the next free handle.
//! - CLASSES, empty: no entity or object here needs a class.
//! - TABLES: all nine, in the reference's order. VPORT, VIEW and UCS are
//!   empty; LTYPE holds `ByBlock`, `ByLayer`, `CONTINUOUS` and `DASHED`;
//!   LAYER holds layer `0` and the drawing's layers; STYLE holds
//!   `Standard`; APPID holds `ACAD`; DIMSTYLE holds `Standard`;
//!   BLOCK_RECORD holds `*Model_Space` and `*Paper_Space`.
//! - BLOCKS: the two spaces' BLOCK and ENDBLK.
//! - ENTITIES: each owned by `*Model_Space`.
//! - OBJECTS: the root dictionary, naming the empty `ACAD_GROUP`
//!   dictionary.

use super::{DxfCurve, DxfText, real};
use core::fmt::Display;
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point2;
use std::fmt::Write as _;

/// The unit a drawing's lengths are in, written as `$INSUNITS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DxfUnit {
    /// No unit declared (`$INSUNITS` 0).
    Unitless,
    /// Millimetres (`$INSUNITS` 4).
    Millimetre,
    /// Inches (`$INSUNITS` 1).
    Inch,
}

impl DxfUnit {
    /// The `$INSUNITS` code.
    #[must_use]
    pub fn insunits(self) -> i32 {
        match self {
            Self::Unitless => 0,
            Self::Millimetre => 4,
            Self::Inch => 1,
        }
    }
}

/// A layer: its name, and whether it draws dashed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DxfLayer {
    /// The name: printable ASCII, none of `< > / \ " : ; ? * | = ``, no
    /// space at either end, at most 255 characters. `0` names the layer
    /// every file has.
    pub name: String,
    /// Dashed (`DASHED` linetype) rather than continuous.
    pub dashed: bool,
}

/// What an item draws.
#[derive(Debug, Clone, PartialEq)]
pub enum DxfShape {
    /// A curve, as [`read_dxf_entities`](super::read_dxf_entities) reads it.
    Curve(DxfCurve),
    /// A line of text.
    Text(DxfText),
}

/// One entity to write, and the layer it goes on.
#[derive(Debug, Clone, PartialEq)]
pub struct DxfItem {
    /// Index into [`DxfOut::layers`].
    pub layer: usize,
    /// What it draws.
    pub shape: DxfShape,
}

/// A drawing to write with [`write_dxf_entities`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DxfOut<'a> {
    /// The unit of every length.
    pub unit: DxfUnit,
    /// The layers, each declared in the LAYER table whether drawn on or not.
    pub layers: &'a [DxfLayer],
    /// The entities, written in order.
    pub items: &'a [DxfItem],
}

/// Write a drawing as an R2000 (AC1015) ASCII DXF, each entity as itself.
///
/// Lines, arcs, circles, ellipses, `LWPOLYLINE`s with their bulges,
/// splines and `TEXT`, on named layers, with the unit as `$INSUNITS`.
/// Coordinates are written in their shortest exact form, so
/// [`read_dxf_entities`](super::read_dxf_entities) reads back the same
/// numbers; angles go out in degrees, as the format has them. A spline
/// with no knots is written by its fit points, as the reader returns one.
/// Text outside printable ASCII is written as `\U+XXXX` escapes, as are a
/// backslash, the first `%` of `%%`, and spaces at either end, so that the
/// string reads back as given in any code page.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// layer's name is not one DXF allows or two names match ignoring case, an
/// item names no layer, a number is not finite, a radius, a text height or
/// a spline weight is not positive, an ellipse's major axis is zero or its
/// ratio outside (0, 1], a polyline has fewer than two vertices or a bulge
/// per vertex, a spline's knots do not number its control points plus its
/// order or decrease, or a text holds a control character or a character
/// outside the Basic Multilingual Plane.
pub fn write_dxf_entities(drawing: &DxfOut<'_>) -> OgeomResult<String> {
    check_layers(drawing.layers)?;
    for item in drawing.items {
        if item.layer >= drawing.layers.len() {
            ogeom_bail!(
                Construction,
                "an item names layer {}, of {} layers",
                item.layer,
                drawing.layers.len()
            );
        }
        check_shape(&item.shape)?;
    }

    let mut w = Writer::default();
    let tables = w.tables(drawing.layers);
    let blocks = w.blocks(&tables);
    let entities = w.entities(drawing, tables.model_space);
    let objects = w.objects();

    let mut out = String::new();
    section(&mut out, "HEADER", &{
        let mut s = String::new();
        variable(&mut s, "$ACADVER", &[(1, "AC1015".into())]);
        variable(&mut s, "$DWGCODEPAGE", &[(3, "ANSI_1252".into())]);
        variable(
            &mut s,
            "$INSUNITS",
            &[(70, drawing.unit.insunits().to_string())],
        );
        let metric = u8::from(drawing.unit != DxfUnit::Inch);
        variable(&mut s, "$MEASUREMENT", &[(70, metric.to_string())]);
        variable(&mut s, "$HANDSEED", &[(5, format!("{:X}", w.next))]);
        s
    });
    section(&mut out, "CLASSES", "");
    section(&mut out, "TABLES", &tables.text);
    section(&mut out, "BLOCKS", &blocks);
    section(&mut out, "ENTITIES", &entities);
    section(&mut out, "OBJECTS", &objects);
    out.push_str("0\nEOF\n");
    Ok(out)
}

fn finite(values: &[f64], what: &str) -> OgeomResult<()> {
    if values.iter().any(|v| !v.is_finite()) {
        ogeom_bail!(Construction, "a {what} has a coordinate that is not finite");
    }
    Ok(())
}

fn positive(value: f64, what: &str) -> OgeomResult<()> {
    if !(value.is_finite() && value > 0.0) {
        ogeom_bail!(Construction, "a {what} must be positive; found {value}");
    }
    Ok(())
}

fn check_layers(layers: &[DxfLayer]) -> OgeomResult<()> {
    for (i, layer) in layers.iter().enumerate() {
        let name = &layer.name;
        let allowed = |c: char| c.is_ascii_graphic() || c == ' ';
        if name.is_empty()
            || name.len() > 255
            || name.starts_with(' ')
            || name.ends_with(' ')
            || !name.chars().all(allowed)
            || name.contains(['<', '>', '/', '\\', '"', ':', ';', '?', '*', '|', '=', '`'])
        {
            ogeom_bail!(Construction, "{name:?} is not a DXF layer name");
        }
        if layers[..i]
            .iter()
            .any(|other| other.name.eq_ignore_ascii_case(name))
        {
            ogeom_bail!(Construction, "layer {name:?} is declared twice");
        }
    }
    Ok(())
}

fn check_shape(shape: &DxfShape) -> OgeomResult<()> {
    match shape {
        DxfShape::Curve(DxfCurve::Line { start, end }) => {
            finite(&[start.x, start.y, end.x, end.y], "line")
        }
        DxfShape::Curve(DxfCurve::Circle { centre, radius }) => {
            finite(&[centre.x, centre.y], "circle")?;
            positive(*radius, "circle's radius")
        }
        DxfShape::Curve(DxfCurve::Arc {
            centre,
            radius,
            start_angle,
            end_angle,
        }) => {
            finite(&[centre.x, centre.y, *start_angle, *end_angle], "arc")?;
            positive(*radius, "arc's radius")
        }
        DxfShape::Curve(DxfCurve::Ellipse {
            centre,
            major,
            ratio,
            start_param,
            end_param,
        }) => {
            finite(
                &[
                    centre.x,
                    centre.y,
                    major.x,
                    major.y,
                    *start_param,
                    *end_param,
                ],
                "ellipse",
            )?;
            if major.x == 0.0 && major.y == 0.0 {
                ogeom_bail!(Construction, "an ellipse's major axis is zero");
            }
            if !(*ratio > 0.0 && *ratio <= 1.0) {
                ogeom_bail!(
                    Construction,
                    "an ellipse's minor-to-major ratio is in (0, 1]; found {ratio}"
                );
            }
            Ok(())
        }
        DxfShape::Curve(DxfCurve::Polyline { vertices, .. }) => {
            if vertices.len() < 2 {
                ogeom_bail!(
                    Construction,
                    "a polyline needs two vertices; found {}",
                    vertices.len()
                );
            }
            for (p, bulge) in vertices {
                finite(&[p.x, p.y, *bulge], "polyline")?;
            }
            Ok(())
        }
        DxfShape::Curve(DxfCurve::Spline {
            degree,
            knots,
            control_points,
            weights,
            ..
        }) => check_spline(*degree, knots, control_points, weights.as_deref()),
        DxfShape::Text(text) => {
            finite(&[text.at.x, text.at.y, text.angle], "text")?;
            positive(text.height, "text's height")?;
            if let Some(c) = text
                .text
                .chars()
                .find(|c| c.is_control() || u32::from(*c) > 0xFFFF)
            {
                ogeom_bail!(
                    Construction,
                    "a DXF text cannot hold {c:?}: a control character or one outside the \
                     Basic Multilingual Plane"
                );
            }
            Ok(())
        }
    }
}

fn check_spline(
    degree: usize,
    knots: &[f64],
    points: &[Point2],
    weights: Option<&[f64]>,
) -> OgeomResult<()> {
    if degree == 0 {
        ogeom_bail!(Construction, "a spline's degree is at least 1");
    }
    for p in points {
        finite(&[p.x, p.y], "spline")?;
    }
    if knots.is_empty() {
        // Fit points: a curve through them, of the degree given.
        if points.len() < 2 {
            ogeom_bail!(Construction, "a spline needs two fit points");
        }
        if weights.is_some() {
            ogeom_bail!(Construction, "a spline by fit points has no weights");
        }
        return Ok(());
    }
    if points.len() <= degree {
        ogeom_bail!(
            Construction,
            "a degree-{degree} spline needs {} control points; found {}",
            degree + 1,
            points.len()
        );
    }
    if knots.len() != points.len() + degree + 1 {
        ogeom_bail!(
            Construction,
            "a degree-{degree} spline with {} control points has {} knots; found {}",
            points.len(),
            points.len() + degree + 1,
            knots.len()
        );
    }
    finite(knots, "spline's knot vector")?;
    if knots.windows(2).any(|w| w[1] < w[0]) {
        ogeom_bail!(Construction, "a spline's knots decrease");
    }
    if let Some(weights) = weights {
        if weights.len() != points.len() {
            ogeom_bail!(
                Construction,
                "a spline with {} control points has {} weights",
                points.len(),
                weights.len()
            );
        }
        for w in weights {
            positive(*w, "spline's weight")?;
        }
    }
    Ok(())
}

/// A string as a DXF value that reads back as itself: anything outside
/// printable ASCII, a backslash, the first `%` of `%%` and a space at
/// either end become `\U+XXXX`.
fn escaped(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    for (i, &c) in chars.iter().enumerate() {
        let plain =
            c.is_ascii_graphic() && c != '\\' && !(c == '%' && chars.get(i + 1) == Some(&'%'))
                || c == ' ' && i != 0 && i + 1 != chars.len();
        if plain {
            out.push(c);
        } else {
            let _ = write!(out, "\\U+{:04X}", u32::from(c));
        }
    }
    out
}

/// The TABLES section, and the block records other records name as owner.
struct Tables {
    text: String,
    model_space: u64,
    paper_space: u64,
}

/// Group-coded text with a running handle count; handles start at 1.
struct Writer {
    next: u64,
}

impl Default for Writer {
    fn default() -> Self {
        Self { next: 1 }
    }
}

fn pair(out: &mut String, code: i32, value: impl Display) {
    let _ = writeln!(out, "{code}\n{value}");
}

fn point(out: &mut String, code: i32, p: Point2) {
    pair(out, code, real(p.x));
    pair(out, code + 10, real(p.y));
    pair(out, code + 20, "0.0");
}

fn section(out: &mut String, name: &str, body: &str) {
    pair(out, 0, "SECTION");
    pair(out, 2, name);
    out.push_str(body);
    pair(out, 0, "ENDSEC");
}

fn variable(out: &mut String, name: &str, values: &[(i32, String)]) {
    pair(out, 9, name);
    for (code, value) in values {
        pair(out, *code, value);
    }
}

impl Writer {
    fn handle(&mut self) -> u64 {
        let h = self.next;
        self.next += 1;
        h
    }

    /// Open a record with its handle and owner, returning the handle.
    fn record(&mut self, out: &mut String, kind: &str, owner: u64) -> u64 {
        self.record_coded(out, kind, 5, owner)
    }

    fn record_coded(&mut self, out: &mut String, kind: &str, code: i32, owner: u64) -> u64 {
        let h = self.handle();
        pair(out, 0, kind);
        pair(out, code, format_args!("{h:X}"));
        pair(out, 330, format_args!("{owner:X}"));
        h
    }

    /// A symbol table with its entries, each written by `entry` given the
    /// table's handle.
    fn table(
        &mut self,
        out: &mut String,
        name: &str,
        count: usize,
        entries: impl FnOnce(&mut Self, &mut String, u64),
    ) {
        let h = self.handle();
        pair(out, 0, "TABLE");
        pair(out, 2, name);
        pair(out, 5, format_args!("{h:X}"));
        pair(out, 330, 0);
        pair(out, 100, "AcDbSymbolTable");
        pair(out, 70, count);
        if name == "DIMSTYLE" {
            pair(out, 100, "AcDbDimStyleTable");
        }
        entries(self, out, h);
        pair(out, 0, "ENDTAB");
    }

    /// Open a table entry: handle, owner, the two subclass markers, name
    /// and flags.
    fn entry(
        &mut self,
        out: &mut String,
        kind: &str,
        owner: u64,
        subclass: &str,
        name: &str,
    ) -> u64 {
        let code = if kind == "DIMSTYLE" { 105 } else { 5 };
        let h = self.record_coded(out, kind, code, owner);
        pair(out, 100, "AcDbSymbolTableRecord");
        pair(out, 100, subclass);
        pair(out, 2, name);
        pair(out, 70, 0);
        h
    }

    fn tables(&mut self, layers: &[DxfLayer]) -> Tables {
        let mut out = String::new();
        self.table(&mut out, "VPORT", 0, |_, _, _| {});

        self.table(&mut out, "LTYPE", 4, |w, out, table| {
            let linetypes: [(&str, &str, &[f64]); 4] = [
                ("ByBlock", "", &[]),
                ("ByLayer", "", &[]),
                ("CONTINUOUS", "Solid line", &[]),
                ("DASHED", "Dashed line", &[0.5, -0.25]),
            ];
            for (name, description, pattern) in linetypes {
                w.entry(out, "LTYPE", table, "AcDbLinetypeTableRecord", name);
                pair(out, 3, description);
                pair(out, 72, 65);
                pair(out, 73, pattern.len());
                pair(
                    out,
                    40,
                    real(pattern.iter().fold(0.0, |sum, d| sum + d.abs())),
                );
                for dash in pattern {
                    pair(out, 49, real(*dash));
                    pair(out, 74, 0);
                }
            }
        });

        // Layer 0 is in every file; a drawing's own `0` replaces it.
        let mut all: Vec<(&str, bool)> = Vec::with_capacity(layers.len() + 1);
        if !layers.iter().any(|l| l.name == "0") {
            all.push(("0", false));
        }
        all.extend(layers.iter().map(|l| (l.name.as_str(), l.dashed)));
        self.table(&mut out, "LAYER", all.len(), |w, out, table| {
            for (name, dashed) in &all {
                w.entry(out, "LAYER", table, "AcDbLayerTableRecord", name);
                pair(out, 62, 7);
                pair(out, 6, if *dashed { "DASHED" } else { "CONTINUOUS" });
            }
        });

        self.table(&mut out, "STYLE", 1, |w, out, table| {
            w.entry(out, "STYLE", table, "AcDbTextStyleTableRecord", "Standard");
            pair(out, 40, "0.0");
            pair(out, 41, "1.0");
            pair(out, 50, "0.0");
            pair(out, 71, 0);
            pair(out, 42, "2.5");
            pair(out, 3, "txt");
            pair(out, 4, "");
        });
        self.table(&mut out, "VIEW", 0, |_, _, _| {});
        self.table(&mut out, "UCS", 0, |_, _, _| {});
        self.table(&mut out, "APPID", 1, |w, out, table| {
            w.entry(out, "APPID", table, "AcDbRegAppTableRecord", "ACAD");
        });
        self.table(&mut out, "DIMSTYLE", 1, |w, out, table| {
            w.entry(
                out,
                "DIMSTYLE",
                table,
                "AcDbDimStyleTableRecord",
                "Standard",
            );
        });

        let mut spaces = (0, 0);
        self.table(&mut out, "BLOCK_RECORD", 2, |w, out, table| {
            spaces.0 = w.entry(
                out,
                "BLOCK_RECORD",
                table,
                "AcDbBlockTableRecord",
                "*Model_Space",
            );
            spaces.1 = w.entry(
                out,
                "BLOCK_RECORD",
                table,
                "AcDbBlockTableRecord",
                "*Paper_Space",
            );
        });
        Tables {
            text: out,
            model_space: spaces.0,
            paper_space: spaces.1,
        }
    }

    fn blocks(&mut self, tables: &Tables) -> String {
        let mut out = String::new();
        for (name, owner, paper) in [
            ("*Model_Space", tables.model_space, false),
            ("*Paper_Space", tables.paper_space, true),
        ] {
            self.record(&mut out, "BLOCK", owner);
            pair(&mut out, 100, "AcDbEntity");
            if paper {
                pair(&mut out, 67, 1);
            }
            pair(&mut out, 8, "0");
            pair(&mut out, 100, "AcDbBlockBegin");
            pair(&mut out, 2, name);
            pair(&mut out, 70, 0);
            point(&mut out, 10, Point2::new(0.0, 0.0));
            pair(&mut out, 3, name);
            pair(&mut out, 1, "");
            self.record(&mut out, "ENDBLK", owner);
            pair(&mut out, 100, "AcDbEntity");
            if paper {
                pair(&mut out, 67, 1);
            }
            pair(&mut out, 8, "0");
            pair(&mut out, 100, "AcDbBlockEnd");
        }
        out
    }

    fn entities(&mut self, drawing: &DxfOut<'_>, owner: u64) -> String {
        let mut out = String::new();
        for item in drawing.items {
            let layer = drawing.layers[item.layer].name.as_str();
            let kind = match &item.shape {
                DxfShape::Curve(DxfCurve::Line { .. }) => "LINE",
                DxfShape::Curve(DxfCurve::Circle { .. }) => "CIRCLE",
                DxfShape::Curve(DxfCurve::Arc { .. }) => "ARC",
                DxfShape::Curve(DxfCurve::Ellipse { .. }) => "ELLIPSE",
                DxfShape::Curve(DxfCurve::Polyline { .. }) => "LWPOLYLINE",
                DxfShape::Curve(DxfCurve::Spline { .. }) => "SPLINE",
                DxfShape::Text(_) => "TEXT",
            };
            self.record(&mut out, kind, owner);
            pair(&mut out, 100, "AcDbEntity");
            pair(&mut out, 8, layer);
            shape(&mut out, &item.shape);
        }
        out
    }

    fn objects(&mut self) -> String {
        let mut out = String::new();
        let root = self.record(&mut out, "DICTIONARY", 0);
        pair(&mut out, 100, "AcDbDictionary");
        pair(&mut out, 281, 1);
        pair(&mut out, 3, "ACAD_GROUP");
        // The group dictionary is the next record, so takes the next handle.
        let group = self.next;
        pair(&mut out, 350, format_args!("{group:X}"));
        self.record(&mut out, "DICTIONARY", root);
        pair(&mut out, 100, "AcDbDictionary");
        pair(&mut out, 281, 1);
        out
    }
}

/// An entity's own groups, after its common ones.
fn shape(out: &mut String, shape: &DxfShape) {
    match shape {
        DxfShape::Curve(DxfCurve::Line { start, end }) => {
            pair(out, 100, "AcDbLine");
            point(out, 10, *start);
            point(out, 11, *end);
        }
        DxfShape::Curve(DxfCurve::Circle { centre, radius }) => {
            pair(out, 100, "AcDbCircle");
            point(out, 10, *centre);
            pair(out, 40, real(*radius));
        }
        DxfShape::Curve(DxfCurve::Arc {
            centre,
            radius,
            start_angle,
            end_angle,
        }) => {
            pair(out, 100, "AcDbCircle");
            point(out, 10, *centre);
            pair(out, 40, real(*radius));
            pair(out, 100, "AcDbArc");
            pair(out, 50, real(start_angle.to_degrees()));
            pair(out, 51, real(end_angle.to_degrees()));
        }
        DxfShape::Curve(DxfCurve::Ellipse {
            centre,
            major,
            ratio,
            start_param,
            end_param,
        }) => {
            pair(out, 100, "AcDbEllipse");
            point(out, 10, *centre);
            point(out, 11, Point2::new(major.x, major.y));
            normal(out);
            pair(out, 40, real(*ratio));
            pair(out, 41, real(*start_param));
            pair(out, 42, real(*end_param));
        }
        DxfShape::Curve(DxfCurve::Polyline { vertices, closed }) => {
            pair(out, 100, "AcDbPolyline");
            pair(out, 90, vertices.len());
            pair(out, 70, u8::from(*closed));
            pair(out, 43, "0.0");
            for (p, bulge) in vertices {
                pair(out, 10, real(p.x));
                pair(out, 20, real(p.y));
                if *bulge != 0.0 {
                    pair(out, 42, real(*bulge));
                }
            }
        }
        DxfShape::Curve(DxfCurve::Spline {
            degree,
            knots,
            control_points,
            weights,
            closed,
        }) => {
            pair(out, 100, "AcDbSpline");
            normal(out);
            // Planar, plus rational and closed as given.
            let flags = 8 | if weights.is_some() { 4 } else { 0 } | u8::from(*closed);
            pair(out, 70, flags);
            pair(out, 71, degree);
            if knots.is_empty() {
                pair(out, 72, 0);
                pair(out, 73, 0);
                pair(out, 74, control_points.len());
                pair(out, 44, "0.0");
                for p in control_points {
                    point(out, 11, *p);
                }
            } else {
                pair(out, 72, knots.len());
                pair(out, 73, control_points.len());
                pair(out, 74, 0);
                for k in knots {
                    pair(out, 40, real(*k));
                }
                for w in weights.iter().flatten() {
                    pair(out, 41, real(*w));
                }
                for p in control_points {
                    point(out, 10, *p);
                }
            }
        }
        DxfShape::Text(text) => {
            pair(out, 100, "AcDbText");
            point(out, 10, text.at);
            pair(out, 40, real(text.height));
            pair(out, 1, escaped(&text.text));
            pair(out, 50, real(text.angle.to_degrees()));
            pair(out, 7, "Standard");
            pair(out, 100, "AcDbText");
        }
    }
}

/// The extrusion of a drawing in the XY plane, seen from above.
fn normal(out: &mut String) {
    pair(out, 210, "0.0");
    pair(out, 220, "0.0");
    pair(out, 230, "1.0");
}
