//! PLY in all three of its encodings, for meshes and for point clouds.
//!
//! A PLY file is a text header declaring elements, each a count and a list
//! of typed properties, then the elements' data in that order: as text, or
//! as packed little- or big-endian binary. A scan keeps everything a
//! scanner measured at a vertex in it, so [`read_ply_bytes`] reads the
//! vertex element whole into a [`PointCloud`]: the position, a normal, a
//! colour, and every other scalar property by its name. Faces are
//! optional; a file without them is a cloud.
//!
//! The header's counts are promises the body must keep, and each is held
//! to the bytes left before anything is allocated for it, so a short or
//! hostile file is refused by name rather than read into memory it does
//! not back.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;
use std::fmt::Write as _;

use crate::cloud::{PointCloud, Scalar};

/// What [`read_ply_bytes`] found.
#[derive(Debug, Clone)]
pub struct PlyImport {
    /// The vertex element: positions and every other property read at
    /// each vertex, colours scaled into `[0, 1]`.
    pub vertices: PointCloud,
    /// The faces, each fanned from its first corner into triangles; empty
    /// for a point cloud.
    pub triangles: Vec<[u32; 3]>,
    /// An RGBA colour per triangle, where the face element carries one;
    /// a face fanned into several triangles gives each its colour.
    pub triangle_colours: Option<Vec<[f64; 4]>>,
    /// What was read with a caveat or passed over: elements other than
    /// vertices and faces, list properties of a vertex, faces of fewer
    /// than three corners.
    pub warnings: Vec<String>,
}

impl PlyImport {
    /// The triangles as a mesh over the vertices.
    ///
    /// A vertex keeps the file's normal where it has a unit one; any other
    /// takes the area-weighted normal of the triangles around it, and a
    /// vertex on none of them points along `+z`.
    #[must_use]
    pub fn mesh(&self) -> Triangulation {
        let positions = self.vertices.points.clone();
        let mut resolved = self
            .vertices
            .normals
            .clone()
            .unwrap_or_else(|| vec![Vector::ZERO; positions.len()]);
        for triangle in &self.triangles {
            let [a, b, c] = triangle.map(|i| positions[i as usize]);
            let face = (b - a).cross(c - a);
            for index in triangle {
                let slot = &mut resolved[*index as usize];
                if crate::mesh_formats::held_is_unset(slot) {
                    *slot += face;
                }
            }
        }
        let normals = resolved
            .into_iter()
            .map(|n| {
                let m = n.magnitude();
                if m > 0.0 { n / m } else { Vector::Z }
            })
            .collect();
        Triangulation {
            parameters: vec![(0.0, 0.0); positions.len()],
            positions,
            normals,
            triangles: self.triangles.clone(),
            deflection_met: true,
        }
    }
}

/// A scalar type the header names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
    /// A type the format does not define, readable only as text, where
    /// every value is a number whatever its type.
    Unknown,
}

impl Kind {
    fn parse(word: &str) -> Self {
        match word {
            "char" | "int8" => Self::I8,
            "uchar" | "uint8" => Self::U8,
            "short" | "int16" => Self::I16,
            "ushort" | "uint16" => Self::U16,
            "int" | "int32" => Self::I32,
            "uint" | "uint32" => Self::U32,
            "float" | "float32" => Self::F32,
            "double" | "float64" => Self::F64,
            _ => Self::Unknown,
        }
    }

    fn size(self) -> usize {
        match self {
            Self::I8 | Self::U8 => 1,
            Self::I16 | Self::U16 => 2,
            Self::I32 | Self::U32 | Self::F32 => 4,
            Self::F64 => 8,
            Self::Unknown => 0,
        }
    }

    fn is_integer(self) -> bool {
        !matches!(self, Self::F32 | Self::F64 | Self::Unknown)
    }

    /// What full intensity is for a colour channel of this type: the
    /// type's largest value for an integer, one for a float.
    fn full(self) -> f64 {
        match self {
            Self::I8 => f64::from(i8::MAX),
            Self::U8 => f64::from(u8::MAX),
            Self::I16 => f64::from(i16::MAX),
            Self::U16 => f64::from(u16::MAX),
            Self::I32 => f64::from(i32::MAX),
            Self::U32 => f64::from(u32::MAX),
            Self::F32 | Self::F64 | Self::Unknown => 1.0,
        }
    }
}

#[derive(Debug, Clone)]
enum Property {
    Scalar {
        name: String,
        kind: Kind,
    },
    List {
        name: String,
        count: Kind,
        item: Kind,
    },
}

impl Property {
    fn name(&self) -> &str {
        match self {
            Self::Scalar { name, .. } | Self::List { name, .. } => name,
        }
    }
}

#[derive(Debug)]
struct Element {
    name: String,
    count: usize,
    properties: Vec<Property>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Encoding {
    Ascii,
    Little,
    Big,
}

/// The header, and where the body starts.
fn header(bytes: &[u8]) -> OgeomResult<(Encoding, Vec<Element>, usize)> {
    let Some(end) = bytes.windows(10).position(|w| w == b"end_header") else {
        if bytes.starts_with(b"ply") {
            ogeom_bail!(Construction, "the PLY header never ends");
        }
        ogeom_bail!(Construction, "this is not a PLY file");
    };
    let mut body = end + 10;
    if bytes.get(body) == Some(&b'\r') {
        body += 1;
    }
    if bytes.get(body) == Some(&b'\n') {
        body += 1;
    }
    let text = String::from_utf8_lossy(&bytes[..end]);
    let mut lines = text.lines();
    if lines.next().map(str::trim) != Some("ply") {
        ogeom_bail!(Construction, "this is not a PLY file");
    }
    let mut encoding = None;
    let mut elements: Vec<Element> = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.as_slice() {
            ["format", kind, ..] => {
                encoding = Some(match *kind {
                    "ascii" => Encoding::Ascii,
                    "binary_little_endian" => Encoding::Little,
                    "binary_big_endian" => Encoding::Big,
                    _ => ogeom_bail!(
                        Construction,
                        "the PLY format {kind} is not one the format defines"
                    ),
                });
            }
            ["element", name, count] => {
                let Ok(count) = count.parse::<usize>() else {
                    ogeom_bail!(Construction, "the PLY element {name} has the count {count}");
                };
                elements.push(Element {
                    name: (*name).to_owned(),
                    count,
                    properties: Vec::new(),
                });
            }
            ["property", "list", count, item, name] => {
                // A property before any element belongs to nothing, and is
                // passed over.
                let Some(element) = elements.last_mut() else {
                    continue;
                };
                let count = Kind::parse(count);
                if !count.is_integer() && count != Kind::Unknown {
                    ogeom_bail!(Construction, "the PLY list {name} is counted by a float");
                }
                element.properties.push(Property::List {
                    name: (*name).to_owned(),
                    count,
                    item: Kind::parse(item),
                });
            }
            ["property", kind, name] => {
                let Some(element) = elements.last_mut() else {
                    continue;
                };
                element.properties.push(Property::Scalar {
                    name: (*name).to_owned(),
                    kind: Kind::parse(kind),
                });
            }
            ["element" | "property", ..] => {
                ogeom_bail!(Construction, "the PLY header line {line:?} is malformed");
            }
            _ => {}
        }
    }
    // A header without a format line is read as text, the one encoding
    // whose bytes say what they are.
    let encoding = encoding.unwrap_or(Encoding::Ascii);
    if encoding != Encoding::Ascii {
        for element in &elements {
            for property in &element.properties {
                let kinds = match property {
                    Property::Scalar { kind, .. } => [*kind, *kind],
                    Property::List { count, item, .. } => [*count, *item],
                };
                if kinds.contains(&Kind::Unknown) {
                    ogeom_bail!(
                        Construction,
                        "the PLY property {} has a type the format does not define, so its size is unknown",
                        property.name()
                    );
                }
            }
        }
    }
    Ok((encoding, elements, body))
}

/// The body, read one value at a time.
struct Body<'a> {
    bytes: &'a [u8],
    at: usize,
    encoding: Encoding,
}

impl Body<'_> {
    fn remaining(&self) -> usize {
        self.bytes.len() - self.at
    }

    fn value(&mut self, kind: Kind) -> OgeomResult<f64> {
        if self.encoding == Encoding::Ascii {
            while self.at < self.bytes.len() && self.bytes[self.at].is_ascii_whitespace() {
                self.at += 1;
            }
            let start = self.at;
            while self.at < self.bytes.len() && !self.bytes[self.at].is_ascii_whitespace() {
                self.at += 1;
            }
            if start == self.at {
                ogeom_bail!(
                    Construction,
                    "the PLY body ends before the header's counts are met"
                );
            }
            let word = String::from_utf8_lossy(&self.bytes[start..self.at]);
            let Ok(value) = word.parse::<f64>() else {
                ogeom_bail!(
                    Construction,
                    "the PLY body carries {word}, which is not a number"
                );
            };
            if kind.is_integer() && value.fract() != 0.0 {
                ogeom_bail!(
                    Construction,
                    "the PLY body carries {word} for an integer property"
                );
            }
            return Ok(value);
        }
        let size = kind.size();
        let Some(raw) = self.bytes.get(self.at..self.at + size) else {
            ogeom_bail!(
                Construction,
                "the PLY body ends before the header's counts are met"
            );
        };
        self.at += size;
        let mut b = [0_u8; 8];
        b[..size].copy_from_slice(raw);
        if self.encoding == Encoding::Big {
            b[..size].reverse();
        }
        Ok(match kind {
            Kind::I8 => f64::from(i8::from_le_bytes([b[0]])),
            Kind::U8 => f64::from(b[0]),
            Kind::I16 => f64::from(i16::from_le_bytes([b[0], b[1]])),
            Kind::U16 => f64::from(u16::from_le_bytes([b[0], b[1]])),
            Kind::I32 => f64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Kind::U32 => f64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Kind::F32 => f64::from(f32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            Kind::F64 | Kind::Unknown => f64::from_le_bytes(b),
        })
    }

    /// A list's length, refused where the bytes left could not hold it.
    fn length(&mut self, count: Kind, item: Kind) -> OgeomResult<usize> {
        let value = self.value(count)?;
        if value < 0.0 || value.fract() != 0.0 {
            ogeom_bail!(Construction, "a PLY list has the length {value}");
        }
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a non-negative integer of at most 32 bits"
        )]
        let n = value as usize;
        let least = match self.encoding {
            Encoding::Ascii => n.saturating_mul(2).saturating_sub(1),
            _ => n.saturating_mul(item.size()),
        };
        if least > self.remaining() {
            ogeom_bail!(
                Construction,
                "a PLY list of {n} items runs past the end of the file"
            );
        }
        Ok(n)
    }

    /// Refuse an element whose count the bytes left could not hold, at
    /// the least room each instance takes.
    fn holds(&self, element: &Element) -> OgeomResult<()> {
        let per: usize = element
            .properties
            .iter()
            .map(|p| match (self.encoding, p) {
                (Encoding::Ascii, _) => 2,
                (_, Property::Scalar { kind, .. }) => kind.size(),
                (_, Property::List { count, .. }) => count.size(),
            })
            .sum();
        let least = element.count.saturating_mul(per);
        let room = match self.encoding {
            Encoding::Ascii => self.remaining().saturating_add(1),
            _ => self.remaining(),
        };
        if least > room {
            ogeom_bail!(
                Construction,
                "the PLY header promises {} of {}, more than the {} bytes left could hold",
                element.count,
                element.name,
                self.remaining()
            );
        }
        Ok(())
    }
}

/// Where a colour channel sits among an element's properties, and its type.
fn channel(properties: &[Property], names: [&str; 3]) -> Option<(usize, Kind)> {
    properties.iter().enumerate().find_map(|(i, p)| match p {
        Property::Scalar { name, kind } if names.contains(&name.as_str()) => Some((i, *kind)),
        _ => None,
    })
}

/// The colour channels of an element: red, green, blue, and alpha where
/// it has one.
fn colour_channels(properties: &[Property]) -> Option<[Option<(usize, Kind)>; 4]> {
    let red = channel(properties, ["red", "r", "diffuse_red"]);
    let green = channel(properties, ["green", "g", "diffuse_green"]);
    let blue = channel(properties, ["blue", "b", "diffuse_blue"]);
    let alpha = channel(properties, ["alpha", "a", "diffuse_alpha"]);
    (red.is_some() && green.is_some() && blue.is_some()).then_some([red, green, blue, alpha])
}

fn colour_of(values: &[f64], channels: &[Option<(usize, Kind)>; 4]) -> [f64; 4] {
    channels.map(|c| c.map_or(1.0, |(i, kind)| values[i] / kind.full()))
}

/// Read a PLY, ASCII or binary of either byte order, mesh or point cloud.
///
/// Every element the header declares is read in its order. Of the vertex
/// element: `x`, `y` and `z` (required), `nx`, `ny` and `nz` as a normal,
/// `red`, `green`, `blue` and `alpha` as a colour (also spelled `r`, `g`,
/// `b`, `a` or `diffuse_red` and so on; an integer channel scaled by its
/// type's largest value, a float one taken as it is), and every other
/// scalar property as a [`Scalar`] of its name. Of the face element: the
/// list named `vertex_indices` (or `vertex_index`, or else its first
/// list), of any count and index type, fanned into triangles, and a colour
/// as for vertices. Other elements and properties are passed over.
///
/// PLY declares no unit; coordinates come back as written.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// bytes are not a PLY file, the header is malformed or names a type the
/// format does not define, a count promises more than the file holds, the
/// body ends early or carries a value that does not parse, the vertices
/// lack a coordinate or are absent, or a face names a vertex that is not
/// there.
pub fn read_ply_bytes(bytes: &[u8]) -> OgeomResult<PlyImport> {
    let (encoding, elements, start) = header(bytes)?;
    let mut body = Body {
        bytes,
        at: start,
        encoding,
    };
    let mut vertices: Option<PointCloud> = None;
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut triangle_colours: Vec<[f64; 4]> = Vec::new();
    let mut face_coloured = false;
    let mut warnings = Vec::new();
    let mut short_faces = 0_usize;
    let mut values: Vec<f64> = Vec::new();
    let mut corners: Vec<u32> = Vec::new();
    for element in &elements {
        if element.properties.is_empty() || element.count == 0 {
            continue;
        }
        body.holds(element)?;
        match element.name.as_str() {
            "vertex" => {
                let at = |name: &str| {
                    element
                        .properties
                        .iter()
                        .position(|p| matches!(p, Property::Scalar { name: n, .. } if n == name))
                };
                let (Some(x), Some(y), Some(z)) = (at("x"), at("y"), at("z")) else {
                    ogeom_bail!(Construction, "a vertex is missing a coordinate");
                };
                let normal = match (at("nx"), at("ny"), at("nz")) {
                    (Some(a), Some(b), Some(c)) => Some([a, b, c]),
                    _ => None,
                };
                let colour = colour_channels(&element.properties);
                let mut taken = vec![x, y, z];
                taken.extend(normal.iter().flatten());
                taken.extend(colour.iter().flatten().flatten().map(|(i, _)| *i));
                let scalar_at: Vec<usize> = (0..element.properties.len())
                    .filter(|i| {
                        !taken.contains(i)
                            && matches!(element.properties[*i], Property::Scalar { .. })
                    })
                    .collect();
                if element
                    .properties
                    .iter()
                    .any(|p| matches!(p, Property::List { .. }))
                {
                    warnings
                        .push("the vertices carry list properties, which are not read".to_owned());
                }
                let mut cloud = PointCloud {
                    points: Vec::with_capacity(element.count),
                    normals: normal.map(|_| Vec::with_capacity(element.count)),
                    colours: colour.map(|_| Vec::with_capacity(element.count)),
                    scalars: scalar_at
                        .iter()
                        .map(|&i| Scalar {
                            name: element.properties[i].name().to_owned(),
                            values: Vec::with_capacity(element.count),
                        })
                        .collect(),
                };
                for _ in 0..element.count {
                    read_instance(&mut body, &element.properties, &mut values, None)?;
                    cloud
                        .points
                        .push(Point::new(values[x], values[y], values[z]));
                    if let (Some(normals), Some([a, b, c])) = (cloud.normals.as_mut(), normal) {
                        normals.push(Vector::new(values[a], values[b], values[c]));
                    }
                    if let (Some(colours), Some(channels)) = (cloud.colours.as_mut(), &colour) {
                        colours.push(colour_of(&values, channels));
                    }
                    for (scalar, &i) in cloud.scalars.iter_mut().zip(&scalar_at) {
                        scalar.values.push(values[i]);
                    }
                }
                if vertices.is_some() {
                    ogeom_bail!(Construction, "the PLY declares its vertices twice");
                }
                vertices = Some(cloud);
            }
            "face" => {
                let lists = || {
                    element
                        .properties
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| matches!(p, Property::List { .. }))
                };
                let indices = lists()
                    .find(|(_, p)| matches!(p.name(), "vertex_indices" | "vertex_index"))
                    .or_else(|| lists().next())
                    .map(|(i, _)| i);
                let Some(indices) = indices else {
                    ogeom_bail!(Construction, "the PLY faces carry no list of vertices");
                };
                let colour = colour_channels(&element.properties);
                face_coloured = colour.is_some();
                for _ in 0..element.count {
                    read_instance(
                        &mut body,
                        &element.properties,
                        &mut values,
                        Some((indices, &mut corners)),
                    )?;
                    if corners.len() < 3 {
                        short_faces += 1;
                        continue;
                    }
                    for i in 1..corners.len() - 1 {
                        triangles.push([corners[0], corners[i], corners[i + 1]]);
                        if let Some(channels) = &colour {
                            triangle_colours.push(colour_of(&values, channels));
                        }
                    }
                }
            }
            other => {
                warnings.push(format!("the element {other} is not read"));
                for _ in 0..element.count {
                    read_instance(&mut body, &element.properties, &mut values, None)?;
                }
            }
        }
    }
    let Some(vertices) = vertices.filter(|v| !v.is_empty()) else {
        ogeom_bail!(Construction, "the file carries no vertices");
    };
    if let Some(past) = triangles
        .iter()
        .flatten()
        .find(|&&v| v as usize >= vertices.len())
    {
        ogeom_bail!(
            Construction,
            "a face names vertex {past}, past the {} the file has",
            vertices.len()
        );
    }
    if short_faces > 0 {
        warnings.push(format!(
            "{short_faces} faces of fewer than three corners were passed over"
        ));
    }
    Ok(PlyImport {
        vertices,
        triangles,
        triangle_colours: face_coloured.then_some(triangle_colours),
        warnings,
    })
}

/// Read one instance of an element: its scalars into `values` at their
/// property positions (a list's slot left at zero), and the list at
/// `corners.0`, if asked for, into `corners.1` as vertex indices.
fn read_instance(
    body: &mut Body<'_>,
    properties: &[Property],
    values: &mut Vec<f64>,
    mut corners: Option<(usize, &mut Vec<u32>)>,
) -> OgeomResult<()> {
    values.clear();
    for (k, property) in properties.iter().enumerate() {
        match property {
            Property::Scalar { kind, .. } => values.push(body.value(*kind)?),
            Property::List { count, item, .. } => {
                values.push(0.0);
                let n = body.length(*count, *item)?;
                match corners.as_mut() {
                    Some((at, out)) if *at == k => {
                        out.clear();
                        for _ in 0..n {
                            let v = body.value(*item)?;
                            if v < 0.0 || v > f64::from(u32::MAX) || v.fract() != 0.0 {
                                ogeom_bail!(Construction, "a face names vertex {v}");
                            }
                            #[allow(
                                clippy::cast_possible_truncation,
                                clippy::cast_sign_loss,
                                reason = "an integer within u32, checked above"
                            )]
                            out.push(v as u32);
                        }
                    }
                    _ => {
                        for _ in 0..n {
                            body.value(*item)?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// The property names the writer gives the position, normal and colour.
const RESERVED: [&str; 10] = [
    "x", "y", "z", "nx", "ny", "nz", "red", "green", "blue", "alpha",
];

/// Write a cloud, or a mesh over it, as binary little-endian PLY.
///
/// Positions as `double`, so they come back as the same `f64`; normals as
/// `float`; colours as `uchar` channels, alpha only where some point is
/// not opaque; each scalar as a `double` property of its name. With
/// triangles, a face element of `uchar`-counted `int` index lists follows;
/// without, the file is a point cloud. [`read_ply_bytes`] reads it back.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if an
/// attribute does not have one entry per point (see [`PointCloud::check`]);
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) for a triangle
/// naming a vertex that is not there, or past the `int` range;
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a scalar name
/// that is empty, holds white space, or is one the writer gives the
/// position, normal or colour.
pub fn write_ply_binary(vertices: &PointCloud, triangles: &[[u32; 3]]) -> OgeomResult<Vec<u8>> {
    vertices.check()?;
    let n = vertices.len();
    if let Some(past) = triangles
        .iter()
        .flatten()
        .find(|&&v| v as usize >= n || i32::try_from(v).is_err())
    {
        ogeom_bail!(Range, "a triangle names vertex {past} of a cloud of {n}");
    }
    for s in &vertices.scalars {
        if s.name.is_empty()
            || s.name.chars().any(char::is_whitespace)
            || RESERVED.contains(&s.name.as_str())
        {
            ogeom_bail!(
                Domain,
                "the scalar name {:?} cannot name a PLY property",
                s.name
            );
        }
    }
    let alpha = vertices
        .colours
        .as_ref()
        .is_some_and(|c| c.iter().any(|c| c[3] < 1.0));
    let mut head = String::from("ply\nformat binary_little_endian 1.0\ncomment ogeom\n");
    let _ = writeln!(head, "element vertex {n}");
    head.push_str("property double x\nproperty double y\nproperty double z\n");
    if vertices.normals.is_some() {
        head.push_str("property float nx\nproperty float ny\nproperty float nz\n");
    }
    if vertices.colours.is_some() {
        head.push_str("property uchar red\nproperty uchar green\nproperty uchar blue\n");
        if alpha {
            head.push_str("property uchar alpha\n");
        }
    }
    for s in &vertices.scalars {
        let _ = writeln!(head, "property double {}", s.name);
    }
    if !triangles.is_empty() {
        let _ = writeln!(head, "element face {}", triangles.len());
        head.push_str("property list uchar int vertex_indices\n");
    }
    head.push_str("end_header\n");
    let mut out = head.into_bytes();
    for (i, p) in vertices.points.iter().enumerate() {
        for c in [p.x, p.y, p.z] {
            out.extend_from_slice(&c.to_le_bytes());
        }
        if let Some(normals) = &vertices.normals {
            let v = normals[i];
            for c in [v.x, v.y, v.z] {
                #[allow(clippy::cast_possible_truncation, reason = "the file's float")]
                out.extend_from_slice(&(c as f32).to_le_bytes());
            }
        }
        if let Some(colours) = &vertices.colours {
            let channels = if alpha { 4 } else { 3 };
            for &c in &colours[i][..channels] {
                out.push(crate::cloud::channel_byte(c));
            }
        }
        for s in &vertices.scalars {
            out.extend_from_slice(&s.values[i].to_le_bytes());
        }
    }
    for t in triangles {
        out.push(3);
        for &v in t {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::cast_possible_truncation,
        reason = "test code"
    )]
    use super::*;

    fn tetrahedron_cloud() -> PointCloud {
        PointCloud {
            points: vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.25, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 0.0, 0.1),
            ],
            normals: Some(vec![
                Vector::new(-0.6, 0.0, -0.8),
                Vector::X,
                Vector::Y,
                Vector::Z,
            ]),
            colours: Some(vec![
                [1.0, 0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0, 1.0],
                [0.0, 0.0, 1.0, 51.0 / 255.0],
                [1.0, 1.0, 1.0, 1.0],
            ]),
            scalars: vec![Scalar {
                name: "confidence".into(),
                values: vec![0.1, 0.2, 0.3, 1e-300],
            }],
        }
    }

    const TETRA: [[u32; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];

    #[test]
    fn binary_ply_round_trips_a_mesh_and_a_cloud() {
        let cloud = tetrahedron_cloud();
        let bytes = write_ply_binary(&cloud, &TETRA).unwrap();
        let back = read_ply_bytes(&bytes).unwrap();
        assert!(back.warnings.is_empty(), "{:?}", back.warnings);
        assert_eq!(back.triangles, TETRA);
        assert!(back.triangle_colours.is_none());
        let mut expected = cloud.clone();
        // Normals travel as float.
        expected.normals = cloud.normals.as_ref().map(|n| {
            n.iter()
                .map(|v| {
                    Vector::new(
                        f64::from(v.x as f32),
                        f64::from(v.y as f32),
                        f64::from(v.z as f32),
                    )
                })
                .collect()
        });
        assert_eq!(back.vertices, expected);

        let points = write_ply_binary(&cloud, &[]).unwrap();
        let back = read_ply_bytes(&points).unwrap();
        assert!(back.triangles.is_empty());
        assert_eq!(back.vertices.points, cloud.points);
        assert!(!String::from_utf8_lossy(&points).contains("element face"));
    }

    /// The same mesh spelled three ways (text, little-endian with `ushort`
    /// counts and `uint` indices, big-endian with `char` counts and
    /// `short` indices) reads to the same vertices and triangles, a quad
    /// fanned into two of them, each with the face's colour.
    #[test]
    fn every_encoding_and_list_type_reads_alike() {
        let ascii = b"ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\n\
property float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty float quality\n\
element face 1\nproperty list uchar int vertex_indices\nproperty float red\nproperty float green\n\
property float blue\nend_header\n0 0 0 255 0 0 0.5\n1 0 0 0 255 0 1.5\n1 1 0 0 0 255 2.5\n0 1 0 0 0 0 3.5\n\
4 0 1 2 3 0.25 0.5 1\n";
        let reference = read_ply_bytes(ascii).unwrap();
        assert_eq!(reference.triangles, vec![[0, 1, 2], [0, 2, 3]]);
        assert_eq!(
            reference.triangle_colours,
            Some(vec![[0.25, 0.5, 1.0, 1.0]; 2])
        );
        assert_eq!(
            reference.vertices.scalar("quality"),
            Some(&[0.5, 1.5, 2.5, 3.5][..])
        );
        assert_eq!(
            reference.vertices.colours.as_ref().unwrap()[2],
            [0.0, 0.0, 1.0, 1.0]
        );

        for big in [false, true] {
            let (count, index) = if big {
                ("char", "short")
            } else {
                ("ushort", "uint")
            };
            let mut bytes = format!(
                "ply\nformat {} 1.0\nelement vertex 4\nproperty float x\nproperty float y\n\
property float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty float quality\n\
element face 1\nproperty list {count} {index} vertex_indices\nproperty float red\nproperty float green\n\
property float blue\nend_header\n",
                if big { "binary_big_endian" } else { "binary_little_endian" }
            )
            .into_bytes();
            let f = |out: &mut Vec<u8>, v: f32| {
                out.extend_from_slice(&if big {
                    v.to_be_bytes()
                } else {
                    v.to_le_bytes()
                });
            };
            let vertices = [
                ([0.0, 0.0, 0.0], [255, 0, 0], 0.5),
                ([1.0, 0.0, 0.0], [0, 255, 0], 1.5),
                ([1.0, 1.0, 0.0], [0, 0, 255], 2.5),
                ([0.0, 1.0, 0.0], [0, 0, 0], 3.5),
            ];
            for (p, c, q) in vertices {
                for v in p {
                    f(&mut bytes, v);
                }
                bytes.extend_from_slice(&c);
                f(&mut bytes, q);
            }
            if big {
                bytes.push(4);
                for i in 0..4_i16 {
                    bytes.extend_from_slice(&i.to_be_bytes());
                }
            } else {
                bytes.extend_from_slice(&4_u16.to_le_bytes());
                for i in 0..4_u32 {
                    bytes.extend_from_slice(&i.to_le_bytes());
                }
            }
            for v in [0.25, 0.5, 1.0] {
                f(&mut bytes, v);
            }
            let read = read_ply_bytes(&bytes).unwrap();
            assert_eq!(read.triangles, reference.triangles, "big endian: {big}");
            assert_eq!(read.triangle_colours, reference.triangle_colours);
            assert_eq!(read.vertices, reference.vertices);
        }
    }

    #[test]
    fn truncated_and_hostile_ply_is_refused_by_name() {
        let bytes = write_ply_binary(&tetrahedron_cloud(), &TETRA).unwrap();
        for cut in 0..bytes.len() {
            assert!(read_ply_bytes(&bytes[..cut]).is_err(), "cut at {cut}");
        }
        // A count no file of this size could back is refused before any
        // allocation is made for it.
        let huge = b"ply\nformat binary_little_endian 1.0\nelement vertex 4000000000000\n\
property double x\nproperty double y\nproperty double z\nend_header\n\0\0\0\0";
        let error = read_ply_bytes(huge).unwrap_err().to_string();
        assert!(error.contains("promises"), "{error}");
        let long_list = b"ply\nformat binary_little_endian 1.0\nelement vertex 1\nproperty uchar x\n\
property uchar y\nproperty uchar z\nelement face 1\nproperty list uint int vertex_indices\nend_header\n\
\0\0\0\xff\xff\xff\xff";
        assert!(read_ply_bytes(long_list).is_err());
        for bad in [
            &b"ply\nformat binary_little_endian 1.0\nelement vertex 1\nproperty quad x\nend_header\n0\n"[..],
            b"ply\nformat zip 1.0\nend_header\n",
            b"ply\nelement vertex 1\nproperty float x\nend_header\n0\n",
            b"ply\nformat ascii 1.0\nproperty float x\nend_header\n",
            b"ply\nformat ascii 1.0\nelement vertex x\nproperty float x\nend_header\n",
            b"ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nend_header\n1\n",
            b"ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nproperty float y\n\
property float z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n3 0 1 -2\n",
            b"\x00\x01garbage end_header",
            b"",
        ] {
            assert!(read_ply_bytes(bad).is_err(), "{}", String::from_utf8_lossy(bad));
        }
        // A type the format does not define is read as text, where a value
        // is a number whatever its type.
        let wide = read_ply_bytes(
            b"ply\nformat ascii 1.0\nproperty float stray\nelement vertex 1\nproperty int64 x\n\
property int64 y\nproperty int64 z\nend_header\n1 2 3\n",
        )
        .unwrap();
        assert_eq!(wide.vertices.points, vec![Point::new(1.0, 2.0, 3.0)]);
        let mut seed = 0x2545_f491_u32;
        for _ in 0..2000 {
            let mut noisy = bytes.clone();
            for _ in 0..4 {
                seed ^= seed << 13;
                seed ^= seed >> 17;
                seed ^= seed << 5;
                let at = seed as usize % noisy.len();
                noisy[at] = (seed >> 8) as u8;
            }
            let _ = read_ply_bytes(&noisy);
        }
    }

    #[test]
    fn the_writer_refuses_what_it_cannot_say() {
        let cloud = tetrahedron_cloud();
        assert!(write_ply_binary(&cloud, &[[0, 1, 4]]).is_err());
        let mut named = cloud.clone();
        named.scalars[0].name = "nx".into();
        assert!(write_ply_binary(&named, &[]).is_err());
        named.scalars[0].name = "two words".into();
        assert!(write_ply_binary(&named, &[]).is_err());
        let mut short = cloud;
        short.colours.as_mut().unwrap().pop();
        assert!(write_ply_binary(&short, &[]).is_err());
    }
}
