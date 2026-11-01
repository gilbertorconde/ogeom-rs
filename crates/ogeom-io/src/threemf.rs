//! 3MF: the mesh package, written with its own archive.
//!
//! A 3MF file is a ZIP holding three parts: the content-types declaration,
//! a relationship pointing at the model, and the model itself: an XML
//! document of meshes and the items that place them. The archive is written
//! and read here, with no library behind it.
//!
//! It writes *stored* entries (no compression), which the ZIP format allows
//! and every reader accepts: what writing a package needs is the container,
//! not the codec. A file written this way is larger than one deflated.
//!
//! Reading is the other way round. Every package a slicer or a modelling
//! tool writes is deflated, so [`read_3mf`] inflates, with the decoder half
//! of DEFLATE kept in this crate, and reads the archive through its central
//! directory, where the sizes are, since a streamed entry leaves them out
//! of its local header.

use ogeom_core::FastMap;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;
use std::fmt::Write as _;

use crate::image::{Image, Wrap};
use crate::mesh_formats::MeshColours;
use crate::xml;

/// One mesh in the package, with the name the item carries.
#[derive(Debug, Clone)]
pub struct Object<'a> {
    /// The tessellation.
    pub mesh: &'a Triangulation,
    /// The name the object is given, if any.
    pub name: Option<String>,
}

/// Write meshes as a 3MF package.
///
/// One object per mesh, one item per object, all in millimetres, the
/// format's own default unit and this kernel's.
#[must_use]
pub fn write_3mf(objects: &[Object<'_>]) -> Vec<u8> {
    let content_types = r#"<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>
</Types>
"#;
    let relationships = r#"<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>
</Relationships>
"#;

    let mut model = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<model unit=\"millimeter\" \
         xml:lang=\"en-US\" \
         xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n  \
         <resources>\n",
    );
    let mut items = String::new();
    for (i, object) in objects.iter().enumerate() {
        if object.mesh.triangles.is_empty() {
            continue;
        }
        let id = i + 1;
        match &object.name {
            Some(name) => {
                let _ = writeln!(
                    model,
                    "    <object id=\"{id}\" type=\"model\" name=\"{}\">",
                    escaped(name)
                );
            }
            None => {
                let _ = writeln!(model, "    <object id=\"{id}\" type=\"model\">");
            }
        }
        model.push_str("      <mesh>\n        <vertices>\n");
        for p in &object.mesh.positions {
            let _ = writeln!(
                model,
                "          <vertex x=\"{}\" y=\"{}\" z=\"{}\"/>",
                p.x, p.y, p.z
            );
        }
        model.push_str("        </vertices>\n        <triangles>\n");
        for [a, b, c] in &object.mesh.triangles {
            let _ = writeln!(
                model,
                "          <triangle v1=\"{a}\" v2=\"{b}\" v3=\"{c}\"/>"
            );
        }
        model.push_str("        </triangles>\n      </mesh>\n    </object>\n");
        let _ = writeln!(items, "    <item objectid=\"{id}\"/>");
    }
    model.push_str("  </resources>\n  <build>\n");
    model.push_str(&items);
    model.push_str("  </build>\n</model>\n");

    archive(&[
        ("[Content_Types].xml", content_types.as_bytes()),
        ("_rels/.rels", relationships.as_bytes()),
        ("3D/3dmodel.model", model.as_bytes()),
    ])
}

/// The XML text escapes, which are the five the specification names.
fn escaped(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// A ZIP archive of stored (uncompressed) entries.
///
/// Local header, data, then a central directory and its end record. Every
/// field the format requires and none it does not: no data descriptors (the
/// sizes are known before writing), no ZIP64 (a 3MF over four gigabytes is
/// a different problem), no extra fields.
fn archive(entries: &[(&str, &[u8])]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    let mut directory: Vec<u8> = Vec::new();
    let mut count = 0_u16;

    for (name, data) in entries {
        let offset = u32::try_from(out.len()).unwrap_or(u32::MAX);
        let crc = crc32(data);
        let size = u32::try_from(data.len()).unwrap_or(u32::MAX);
        let name_bytes = name.as_bytes();
        let name_len = u16::try_from(name_bytes.len()).unwrap_or(u16::MAX);

        // Local file header.
        out.extend_from_slice(&0x0403_4b50_u32.to_le_bytes());
        out.extend_from_slice(&20_u16.to_le_bytes()); // version needed
        out.extend_from_slice(&0_u16.to_le_bytes()); // flags
        out.extend_from_slice(&0_u16.to_le_bytes()); // stored
        out.extend_from_slice(&0_u16.to_le_bytes()); // time
        out.extend_from_slice(&0_u16.to_le_bytes()); // date
        out.extend_from_slice(&crc.to_le_bytes());
        out.extend_from_slice(&size.to_le_bytes());
        out.extend_from_slice(&size.to_le_bytes());
        out.extend_from_slice(&name_len.to_le_bytes());
        out.extend_from_slice(&0_u16.to_le_bytes()); // extra length
        out.extend_from_slice(name_bytes);
        out.extend_from_slice(data);

        // Central directory entry.
        directory.extend_from_slice(&0x0201_4b50_u32.to_le_bytes());
        directory.extend_from_slice(&20_u16.to_le_bytes()); // version made by
        directory.extend_from_slice(&20_u16.to_le_bytes()); // version needed
        directory.extend_from_slice(&0_u16.to_le_bytes());
        directory.extend_from_slice(&0_u16.to_le_bytes());
        directory.extend_from_slice(&0_u16.to_le_bytes());
        directory.extend_from_slice(&0_u16.to_le_bytes());
        directory.extend_from_slice(&crc.to_le_bytes());
        directory.extend_from_slice(&size.to_le_bytes());
        directory.extend_from_slice(&size.to_le_bytes());
        directory.extend_from_slice(&name_len.to_le_bytes());
        directory.extend_from_slice(&0_u16.to_le_bytes()); // extra
        directory.extend_from_slice(&0_u16.to_le_bytes()); // comment
        directory.extend_from_slice(&0_u16.to_le_bytes()); // disk
        directory.extend_from_slice(&0_u16.to_le_bytes()); // internal attrs
        directory.extend_from_slice(&0_u32.to_le_bytes()); // external attrs
        directory.extend_from_slice(&offset.to_le_bytes());
        directory.extend_from_slice(name_bytes);
        count += 1;
    }

    let directory_offset = u32::try_from(out.len()).unwrap_or(u32::MAX);
    let directory_size = u32::try_from(directory.len()).unwrap_or(u32::MAX);
    out.extend_from_slice(&directory);
    // End of central directory.
    out.extend_from_slice(&0x0605_4b50_u32.to_le_bytes());
    out.extend_from_slice(&0_u16.to_le_bytes()); // this disk
    out.extend_from_slice(&0_u16.to_le_bytes()); // directory's disk
    out.extend_from_slice(&count.to_le_bytes());
    out.extend_from_slice(&count.to_le_bytes());
    out.extend_from_slice(&directory_size.to_le_bytes());
    out.extend_from_slice(&directory_offset.to_le_bytes());
    out.extend_from_slice(&0_u16.to_le_bytes()); // comment length
    out
}

/// The ZIP checksum: CRC-32, reflected, polynomial `0xEDB8_8320`, a byte at
/// a time through its table.
pub(crate) fn crc32(data: &[u8]) -> u32 {
    const TABLE: [u32; 256] = {
        let mut table = [0_u32; 256];
        let mut i = 0_u32;
        while i < 256 {
            let mut crc = i;
            let mut k = 0;
            while k < 8 {
                crc = if crc & 1 != 0 {
                    (crc >> 1) ^ 0xEDB8_8320
                } else {
                    crc >> 1
                };
                k += 1;
            }
            table[i as usize] = crc;
            i += 1;
        }
        table
    };
    let mut crc = 0xFFFF_FFFF_u32;
    for byte in data {
        crc = (crc >> 8) ^ TABLE[((crc ^ u32::from(*byte)) & 0xFF) as usize];
    }
    !crc
}

// --- Reading -----------------------------------------------------------------

/// One entry of the archive's central directory.
struct Entry {
    name: String,
    method: u16,
    flags: u16,
    crc: u32,
    compressed: usize,
    size: usize,
    header: usize,
}

fn u16_at(bytes: &[u8], at: usize) -> OgeomResult<u16> {
    match at.checked_add(2).and_then(|end| bytes.get(at..end)) {
        Some(b) => Ok(u16::from_le_bytes([b[0], b[1]])),
        None => ogeom_bail!(Construction, "the archive is cut short"),
    }
}

fn u32_at(bytes: &[u8], at: usize) -> OgeomResult<u32> {
    match at.checked_add(4).and_then(|end| bytes.get(at..end)) {
        Some(b) => Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        None => ogeom_bail!(Construction, "the archive is cut short"),
    }
}

/// The central directory: every entry, with the sizes a streamed local
/// header leaves as zero.
fn directory(bytes: &[u8]) -> OgeomResult<Vec<Entry>> {
    const END: u32 = 0x0605_4b50;
    const LOCATOR: u32 = 0x0706_4b50;
    const END64: u32 = 0x0606_4b50;
    // The end record is the last thing in the file but for a comment of up
    // to 65 535 bytes.
    let floor = bytes.len().saturating_sub(22 + 0xFFFF);
    let mut end = None;
    let mut at = bytes.len().saturating_sub(22);
    while at >= floor && bytes.len() >= 22 {
        if u32_at(bytes, at)? == END {
            end = Some(at);
            break;
        }
        if at == 0 {
            break;
        }
        at -= 1;
    }
    let Some(end) = end else {
        ogeom_bail!(Construction, "these bytes are not a ZIP archive");
    };
    let zip64 = end >= 20 && u32_at(bytes, end - 20)? == LOCATOR;
    // In a ZIP64 archive the classic record's disk numbers may be
    // saturated like the rest of it; the ZIP64 record's are the ones read.
    if !zip64 && (u16_at(bytes, end + 4)? != 0 || u16_at(bytes, end + 6)? != 0) {
        ogeom_bail!(Construction, "the archive spans several disks");
    }
    let mut count = u64::from(u16_at(bytes, end + 10)?);
    let mut offset = u64::from(u32_at(bytes, end + 16)?);
    // A ZIP64 archive keeps the classic end record with its fields
    // saturated, and the true ones in a ZIP64 end record its locator, just
    // before the classic one, points to. Streaming writers emit it whatever
    // the archive's size.
    if zip64 {
        let at = offset_in(bytes, u64_at(bytes, end - 20 + 8)?)?;
        if u32_at(bytes, at)? != END64 {
            ogeom_bail!(
                Construction,
                "the archive's ZIP64 locator points at no ZIP64 end record"
            );
        }
        if u32_at(bytes, at + 16)? != 0 || u32_at(bytes, at + 20)? != 0 {
            ogeom_bail!(Construction, "the archive spans several disks");
        }
        count = u64_at(bytes, at + 32)?;
        offset = u64_at(bytes, at + 48)?;
    } else if count == 0xFFFF || offset == 0xFFFF_FFFF {
        ogeom_bail!(
            Construction,
            "the archive's end record is saturated and no ZIP64 record follows it"
        );
    }
    let mut entries = Vec::new();
    let mut at = offset_in(bytes, offset)?;
    for _ in 0..count {
        if u32_at(bytes, at)? != 0x0201_4b50 {
            ogeom_bail!(Construction, "the archive's central directory is damaged");
        }
        let name_len = usize::from(u16_at(bytes, at + 28)?);
        let extra_len = usize::from(u16_at(bytes, at + 30)?);
        let comment_len = usize::from(u16_at(bytes, at + 32)?);
        let Some(name) = bytes.get(at + 46..at + 46 + name_len) else {
            ogeom_bail!(Construction, "the archive is cut short");
        };
        let name = String::from_utf8_lossy(name).into_owned();
        let mut size = u64::from(u32_at(bytes, at + 24)?);
        let mut compressed = u64::from(u32_at(bytes, at + 20)?);
        let mut header = u64::from(u32_at(bytes, at + 42)?);
        // A saturated size or offset is carried in the ZIP64 extra field,
        // in the order uncompressed size, compressed size, header offset,
        // each present only where its fixed field is saturated.
        if [size, compressed, header].contains(&0xFFFF_FFFF) {
            let Some(extra) = extra_field(bytes, at + 46 + name_len, extra_len, 0x0001) else {
                ogeom_bail!(
                    Construction,
                    "the entry {name} has a saturated field and no ZIP64 extra"
                );
            };
            let mut read = 0;
            for field in [&mut size, &mut compressed, &mut header] {
                if *field == 0xFFFF_FFFF {
                    *field = u64_at(extra, read)?;
                    read += 8;
                }
            }
        }
        entries.push(Entry {
            name,
            flags: u16_at(bytes, at + 8)?,
            method: u16_at(bytes, at + 10)?,
            crc: u32_at(bytes, at + 16)?,
            compressed: offset_in(bytes, compressed)?,
            size: to_usize(size)?,
            header: offset_in(bytes, header)?,
        });
        at += 46 + name_len + extra_len + comment_len;
    }
    Ok(entries)
}

fn u64_at(bytes: &[u8], at: usize) -> OgeomResult<u64> {
    match at.checked_add(8).and_then(|end| bytes.get(at..end)) {
        Some(b) => Ok(u64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ])),
        None => ogeom_bail!(Construction, "the archive is cut short"),
    }
}

/// An offset or length within the archive, refused past its end: every
/// such field names bytes the file must hold.
fn offset_in(bytes: &[u8], value: u64) -> OgeomResult<usize> {
    match usize::try_from(value) {
        Ok(v) if v <= bytes.len() => Ok(v),
        _ => ogeom_bail!(Construction, "the archive names an offset past its own end"),
    }
}

fn to_usize(value: u64) -> OgeomResult<usize> {
    match usize::try_from(value) {
        Ok(v) => Ok(v),
        Err(_) => ogeom_bail!(
            Construction,
            "the archive names an offset past this machine's reach"
        ),
    }
}

/// The data of the extra field with header `id` among the `len` bytes of
/// extra fields starting at `at`.
fn extra_field(bytes: &[u8], at: usize, len: usize, id: u16) -> Option<&[u8]> {
    let fields = bytes.get(at..at + len)?;
    let mut i = 0;
    while i + 4 <= fields.len() {
        let tag = u16::from_le_bytes([fields[i], fields[i + 1]]);
        let size = usize::from(u16::from_le_bytes([fields[i + 2], fields[i + 3]]));
        let data = fields.get(i + 4..i + 4 + size)?;
        if tag == id {
            return Some(data);
        }
        i += 4 + size;
    }
    None
}

/// An entry's bytes: stored ones as they are, deflated ones inflated, and
/// either checked against the checksum the archive gives.
fn contents(bytes: &[u8], entry: &Entry) -> OgeomResult<Vec<u8>> {
    let name = &entry.name;
    if entry.flags & 1 != 0 {
        ogeom_bail!(Construction, "the entry {name} is encrypted");
    }
    if u32_at(bytes, entry.header)? != 0x0403_4b50 {
        ogeom_bail!(
            Construction,
            "the entry {name} has no local header where the directory puts it"
        );
    }
    let name_len = usize::from(u16_at(bytes, entry.header + 26)?);
    let extra_len = usize::from(u16_at(bytes, entry.header + 28)?);
    let data_at = entry.header + 30 + name_len + extra_len;
    let Some(data) = data_at
        .checked_add(entry.compressed)
        .and_then(|end| bytes.get(data_at..end))
    else {
        ogeom_bail!(
            Construction,
            "the entry {name} runs past the end of the file"
        );
    };
    let out = match entry.method {
        0 => data.to_vec(),
        8 => crate::inflate::inflate(data, entry.size)?,
        method => ogeom_bail!(
            Construction,
            "the entry {name} is compressed with method {method}; only stored and deflated entries are read"
        ),
    };
    if out.len() != entry.size || crc32(&out) != entry.crc {
        ogeom_bail!(Construction, "the entry {name} does not match its checksum");
    }
    Ok(out)
}

/// Read every part of a package: the entry names and their bytes.
///
/// Stored and deflated entries are read, each checked against its
/// checksum, from classic and ZIP64 archives alike; the central directory's
/// sizes are the ones trusted, since a streamed entry's local header leaves
/// them zero or saturated. Any other compression method, an encrypted
/// entry, or an archive spanning several disks is refused by name rather
/// than half-read.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// bytes are not a ZIP archive, an entry is damaged, or it is one of the
/// refusals above.
pub fn read_package(bytes: &[u8]) -> OgeomResult<Vec<(String, Vec<u8>)>> {
    directory(bytes)?
        .iter()
        .filter(|e| !e.name.ends_with('/'))
        .map(|e| Ok((e.name.clone(), contents(bytes, e)?)))
        .collect()
}

/// What a 3MF object is for, as its `type` attribute says.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectType {
    /// A part to be made: `model`, the default.
    Model,
    /// A support the designer modelled as solid: `solidsupport`.
    SolidSupport,
    /// Support structure, not part of the design: `support`.
    Support,
    /// An open surface, not a solid: `surface`.
    Surface,
    /// Anything else the file carries: `other`.
    Other,
}

/// One object of a package, placed as its build item places it.
#[derive(Debug, Clone)]
pub struct ThreeMfObject {
    /// The object's name, where it has one.
    pub name: Option<String>,
    /// What the object is for. Every type is read; the import's warnings
    /// say when one is not a part.
    pub object_type: ObjectType,
    /// The mesh in millimetres, components flattened into it and the build
    /// item's transform applied.
    pub mesh: Triangulation,
    /// The object's colour as the file writes it (sRGB, RGBA in `[0, 1]`),
    /// where every triangle carries the same one.
    pub colour: Option<[f64; 4]>,
    /// What the object is for, as the metadata of the package names its
    /// role; [`ObjectRole::Model`] where nothing does.
    pub role: ObjectRole,
    /// The colours of the mesh as the file writes them (sRGB, RGBA in
    /// `[0, 1]`): from colour groups, base materials or textures sampled
    /// at each corner's coordinates. Given where every triangle has one.
    pub colours: Option<MeshColours>,
    /// The volumes the mesh is made of, in order, each a run of its
    /// triangles: one per mesh the object's components flatten, or per
    /// triangle range the package's metadata divides a mesh into. Together
    /// they cover every triangle of [`ThreeMfObject::mesh`].
    pub volumes: Vec<ThreeMfVolume>,
}

/// What a volume is for, as the metadata a package's writer keeps names
/// its role.
///
/// The core specification has no word for this. Tools that prepare a model
/// for printing write it in their own vocabulary: as metadata of an object
/// in the model (`<metadata name="…">`), or in configuration parts under
/// `Metadata/` (`*.config`) that describe each object's volumes by an
/// attribute or a `key`/`value` pair. Those are read by name, whatever
/// namespace prefix the name carries: a `volume_type` or `subtype` of
/// `ModelPart`/`normal_part`, `ParameterModifier`/`modifier_part`,
/// `NegativeVolume`/`negative_part`, `SupportBlocker` or `SupportEnforcer`
/// (case, underscores and hyphens aside); and the flags `modifier`,
/// `infill_mesh` and `cutting_mesh` (a modifier), `anti_overhang_mesh` (a
/// support blocker) and `support_mesh` (a support enforcer) set true.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ObjectRole {
    /// Part of what is to be made.
    Model,
    /// A volume that changes how the parts it overlaps are made, and is
    /// not made itself.
    Modifier,
    /// A volume where no support is to be generated.
    SupportBlocker,
    /// A volume where support is to be generated.
    SupportEnforcer,
    /// A volume cut out of the parts it overlaps.
    NegativeVolume,
    /// A role the metadata names that is none of these, as written.
    Other(String),
}

impl ObjectRole {
    fn describe(&self) -> String {
        match self {
            Self::Model => "a part".to_owned(),
            Self::Modifier => "a modifier".to_owned(),
            Self::SupportBlocker => "a support blocker".to_owned(),
            Self::SupportEnforcer => "a support enforcer".to_owned(),
            Self::NegativeVolume => "a negative volume".to_owned(),
            Self::Other(role) => format!("of role {role}"),
        }
    }
}

/// One volume of an object: a run of the object's triangles, and what it
/// is for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThreeMfVolume {
    /// The volume's name, from the metadata or else its mesh object.
    pub name: Option<String>,
    /// What the volume is for.
    pub role: ObjectRole,
    /// Its triangles, as indices into the object's mesh.
    pub triangles: std::ops::Range<usize>,
}

/// What [`read_3mf`] found.
#[derive(Debug, Clone)]
pub struct ThreeMfImport {
    /// One object per build item, in build order.
    pub objects: Vec<ThreeMfObject>,
    /// What was read with a caveat, or skipped: required extensions this
    /// reader does not know, per-triangle colours dropped, objects that are
    /// not parts.
    pub warnings: Vec<String>,
}

const CORE: &str = "http://schemas.microsoft.com/3dmanufacturing/core/2015/02";
const MATERIAL: &str = "http://schemas.microsoft.com/3dmanufacturing/material/2015/02";
const PRODUCTION: &str = "http://schemas.microsoft.com/3dmanufacturing/production/2015/06";
const MODEL_RELATIONSHIP: &str = "http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel";
/// The extensions whose required use changes nothing this reader returns:
/// production (read; it is how multi-part packages reference their
/// objects), materials (read as far as colour), and the slicers' own.
const UNDERSTOOD: [&str; 2] = [PRODUCTION, MATERIAL];

/// Read a 3MF package into placed meshes.
///
/// The start part is the one `_rels/.rels` names as the 3D model. Its
/// objects are meshes or assemblies of components; a component, or a
/// build item, may name an object in another model part of the package
/// through the production extension's `path`, as slicers write every
/// object to its own part. Each build item comes back as one mesh, its
/// components flattened and every transform applied, scaled from the
/// model's unit to millimetres. A transform that mirrors has its triangles
/// rewound, so every mesh keeps its outward winding. Vertices the
/// flattening brings together within the confusion tolerance are welded.
///
/// A uniform object colour is kept, from a base material or a colour
/// group. Colours that vary across the triangles are read into
/// [`ThreeMfObject::colours`] where every triangle has one: from colour
/// groups and base materials, and from textures (PNG or JPEG parts of the
/// package) sampled by the nearest texel at each corner's coordinates,
/// tiled as the texture's `tilestyleu` and `tilestylev` say (`none` read as
/// `clamp`). Composite materials and multi-properties are passed over with
/// a warning, and so is a texture that cannot be decoded.
///
/// The roles of objects and volumes (modifiers, support blockers and
/// enforcers, negative volumes) are read from the metadata of the package
/// as [`ObjectRole`] lays out, and each object's [`ThreeMfObject::volumes`]
/// says which of its triangles each volume is. A configuration part that
/// cannot be read is a warning.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// archive cannot be read (see [`read_package`]), a model part is not
/// well-formed, a reference names an object that is not there or contains
/// itself, a triangle names a vertex that is not there, or the model
/// requires the secure-content extension.
pub fn read_3mf(bytes: &[u8], tol: Tolerances) -> OgeomResult<ThreeMfImport> {
    let entries = directory(bytes)?;
    let find = |path: &str| -> Option<&Entry> {
        let wanted = part_key(path);
        entries.iter().find(|e| part_key(&e.name) == wanted)
    };
    let mut warnings = Vec::new();
    let start = match find("_rels/.rels") {
        Some(rels) => {
            let text = String::from_utf8_lossy(&contents(bytes, rels)?).into_owned();
            start_part(&text)?
        }
        None => None,
    }
    .unwrap_or_else(|| "/3D/3dmodel.model".to_owned());

    let mut parts: FastMap<String, Part> = FastMap::default();
    let mut pending = vec![start.clone()];
    while let Some(path) = pending.pop() {
        let key = part_key(&path);
        if parts.contains_key(&key) {
            continue;
        }
        let Some(entry) = find(&path) else {
            ogeom_bail!(Construction, "the package has no model part {path}");
        };
        let raw = contents(bytes, entry)?;
        let Ok(text) = std::str::from_utf8(&raw) else {
            ogeom_bail!(Construction, "the model part {path} is not UTF-8");
        };
        let part = parse_part(text, &path, &mut warnings)?;
        pending.extend(part.referenced());
        parts.insert(key, part);
    }

    // Textures, each decoded once; one that cannot be is a warning, and
    // the triangles it colours go without.
    let mut budget = crate::image::Budget::new();
    for part in parts.values_mut() {
        for def in part.textures.values_mut() {
            let decoded = match find(&def.path) {
                None => Err(format!("the package has no texture part {}", def.path)),
                Some(entry) => contents(bytes, entry)
                    .and_then(|raw| budget.decode(&raw))
                    .map_err(|e| format!("the texture {} could not be decoded: {e}", def.path)),
            };
            match decoded {
                Ok(image) => def.image = Some(image),
                Err(warning) => warnings.push(warning),
            }
        }
    }
    // The configuration parts a writer keeps beside the model: read for the
    // roles of the volumes, and passed over with a warning if malformed.
    let mut config: FastMap<u32, ConfigObject> = FastMap::default();
    for entry in &entries {
        let key = part_key(&entry.name);
        if !(key.starts_with("metadata/") && key.ends_with(".config")) {
            continue;
        }
        let read = contents(bytes, entry).and_then(|raw| match std::str::from_utf8(&raw) {
            Ok(text) => parse_config(text),
            Err(_) => ogeom_bail!(Construction, "it is not UTF-8"),
        });
        match read {
            Ok(objects) => {
                for (id, object) in objects {
                    let held = config.entry(id).or_default();
                    held.pairs.extend(object.pairs);
                    held.volumes.extend(object.volumes);
                }
            }
            Err(e) => warnings.push(format!(
                "the configuration part {} could not be read: {e}",
                entry.name
            )),
        }
    }

    let context = Context {
        parts: &parts,
        start: part_key(&start),
        config: &config,
    };
    let root = &parts[&part_key(&start)];
    let scale = root.scale;
    let mut objects = Vec::with_capacity(root.build.len());
    for item in &root.build {
        let path = item.path.clone().unwrap_or_else(|| start.clone());
        let top = object(&parts, &path, item.object)?;
        let outer = Tag {
            name: None,
            role: ObjectRole::Model,
        };
        let role = context.tag(&path, item.object, top, &outer).role;
        let mut flat = Flat::default();
        flatten(
            &context,
            &path,
            item.object,
            item.transform,
            &mut flat,
            0,
            &outer,
        )?;
        let name = top.name.clone();
        let label = name
            .clone()
            .unwrap_or_else(|| format!("object {}", item.object));
        match top.object_type {
            ObjectType::Model | ObjectType::SolidSupport => {}
            ObjectType::Support => warnings.push(format!("{label} is a support structure")),
            ObjectType::Surface => warnings.push(format!("{label} is a surface, not a solid")),
            ObjectType::Other => warnings.push(format!("{label} is of type other, not a part")),
        }
        if role != ObjectRole::Model {
            warnings.push(format!("{label} is {}, not a part", role.describe()));
        }
        if flat.skipped_properties {
            warnings.push(format!(
                "{label} carries composite materials or multiple properties; they are not read"
            ));
        }
        if flat.undecoded {
            warnings.push(format!(
                "{label} is textured with an image that could not be decoded"
            ));
        }
        let several = flat.colours.len() > 1;
        let colour = match flat.colours.as_slice() {
            [] => None,
            [only] => *only,
            _ => resolve_colour(parts.get(&part_key(&path)), top.pid, top.pindex),
        };
        if flat.degenerate > 0 {
            warnings.push(format!(
                "{label}: {} triangles repeat a vertex and were dropped",
                flat.degenerate
            ));
        }
        let Welded {
            mesh,
            corners,
            volumes,
        } = flat.into_mesh(scale, tol.confusion());
        let coloured = corners.iter().filter(|c| c.is_some()).count();
        let colours = match corners.into_iter().collect::<Option<Vec<_>>>() {
            Some(corners) if !corners.is_empty() => Some(MeshColours::from_corners(
                &mesh.triangles,
                mesh.positions.len(),
                corners,
            )),
            _ => {
                if coloured > 0 || several {
                    warnings.push(format!(
                        "{label} carries several colours and not every triangle has one; the per-triangle colours are dropped"
                    ));
                }
                None
            }
        };
        for volume in &volumes {
            if volume.role != ObjectRole::Model && volume.role != role {
                let named = volume
                    .name
                    .as_ref()
                    .map_or_else(String::new, |n| format!(" {n}"));
                warnings.push(format!(
                    "{label} holds the volume{named}, {}; its triangles are in the mesh",
                    volume.role.describe()
                ));
            }
        }
        objects.push(ThreeMfObject {
            name,
            object_type: top.object_type,
            mesh,
            colour,
            role,
            colours,
            volumes,
        });
    }
    Ok(ThreeMfImport { objects, warnings })
}

/// A part name as the archive and every reference compare it: without the
/// leading slash, percent escapes decoded, and case folded, as the package
/// conventions make part names case-insensitive.
fn part_key(path: &str) -> String {
    let path = path.trim_start_matches('/');
    let bytes = path.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%'
            && let Some(hex) = path.get(i + 1..i + 3)
            && let Ok(byte) = u8::from_str_radix(hex, 16)
        {
            out.push(byte);
            i += 3;
            continue;
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).to_lowercase()
}

/// The target of the package's 3D-model relationship.
fn start_part(rels: &str) -> OgeomResult<Option<String>> {
    let mut reader = xml::Reader::new(rels);
    while let Some(event) = reader.next()? {
        if let xml::Event::Start(start) = event
            && start.local == "Relationship"
            && start.get("Type") == Some(MODEL_RELATIONSHIP)
            && let Some(target) = start.get("Target")
        {
            return Ok(Some(target.to_owned()));
        }
    }
    Ok(None)
}

/// A 3MF transform: the three rows the format writes, image of each axis,
/// then the translation; a point is carried as a row vector.
#[derive(Debug, Clone, Copy)]
struct Affine {
    axes: [Vector; 3],
    translation: Vector,
}

impl Affine {
    const IDENTITY: Self = Self {
        axes: [Vector::X, Vector::Y, Vector::Z],
        translation: Vector::ZERO,
    };

    fn parse(text: &str) -> OgeomResult<Self> {
        let values: Vec<f64> = text
            .split_ascii_whitespace()
            .map(str::parse)
            .collect::<Result<_, _>>()
            .or_else(|_| ogeom_bail!(Construction, "the transform {text} is not twelve numbers"))?;
        let [a, b, c, d, e, f, g, h, i, x, y, z] = values[..] else {
            ogeom_bail!(Construction, "the transform {text} is not twelve numbers");
        };
        Ok(Self {
            axes: [
                Vector::new(a, b, c),
                Vector::new(d, e, f),
                Vector::new(g, h, i),
            ],
            translation: Vector::new(x, y, z),
        })
    }

    fn vector(self, v: Vector) -> Vector {
        self.axes[0] * v.x + self.axes[1] * v.y + self.axes[2] * v.z
    }

    fn point(self, p: Point) -> Point {
        Point::ORIGIN + self.vector(p - Point::ORIGIN) + self.translation
    }

    /// `inner` first, then `self`.
    fn after(self, inner: Self) -> Self {
        Self {
            axes: inner.axes.map(|a| self.vector(a)),
            translation: self.vector(inner.translation) + self.translation,
        }
    }

    fn mirrors(self) -> bool {
        self.axes[0].dot(self.axes[1].cross(self.axes[2])) < 0.0
    }
}

/// A component or a build item: an object, where it lives, and where it
/// is placed.
struct Reference {
    path: Option<String>,
    object: u32,
    transform: Affine,
}

enum Content {
    Mesh {
        positions: Vec<Point>,
        triangles: Vec<[u32; 3]>,
        /// Each triangle's property, `(group, index at each corner)`,
        /// where it names one.
        properties: Vec<Option<(u32, [u32; 3])>>,
    },
    Components(Vec<Reference>),
}

struct ObjectDef {
    name: Option<String>,
    object_type: ObjectType,
    pid: Option<u32>,
    pindex: Option<u32>,
    content: Content,
    /// The `(name, value)` pairs of the object's metadata group.
    metadata: Vec<(String, String)>,
}

/// A property group of a model part.
enum Group {
    /// Colours by index: a colour group, or base materials' display
    /// colours.
    Colours(Vec<[f64; 4]>),
    /// Texture coordinates by index into the texture of that id.
    Texture {
        texture: u32,
        coords: Vec<(f64, f64)>,
    },
    /// A group this reader does not read: composites, multi-properties.
    Unread,
}

/// A `texture2d`: the image part it names, how it tiles, and the image
/// once decoded.
struct TextureDef {
    path: String,
    wrap: [Wrap; 2],
    image: Option<Image>,
}

/// One model part: its objects, property groups and build.
struct Part {
    scale: f64,
    objects: FastMap<u32, ObjectDef>,
    /// Property groups by id.
    groups: FastMap<u32, Group>,
    textures: FastMap<u32, TextureDef>,
    build: Vec<Reference>,
}

impl Part {
    /// The other parts this one's components and items name.
    fn referenced(&self) -> Vec<String> {
        let components = self.objects.values().flat_map(|o| match &o.content {
            Content::Components(c) => c.iter().collect::<Vec<_>>(),
            Content::Mesh { .. } => Vec::new(),
        });
        components
            .chain(&self.build)
            .filter_map(|r| r.path.clone())
            .collect()
    }
}

fn number<T: std::str::FromStr>(start: &xml::Start<'_>, name: &str) -> OgeomResult<Option<T>> {
    match start.get(name) {
        None => Ok(None),
        Some(text) => match text.trim().parse() {
            Ok(v) => Ok(Some(v)),
            Err(_) => ogeom_bail!(
                Construction,
                "the {} attribute {name}=\"{text}\" is not a number",
                start.local
            ),
        },
    }
}

fn required<T: std::str::FromStr>(start: &xml::Start<'_>, name: &str) -> OgeomResult<T> {
    match number(start, name)? {
        Some(v) => Ok(v),
        None => ogeom_bail!(Construction, "a {} has no {name}", start.local),
    }
}

fn parse_part(text: &str, path: &str, warnings: &mut Vec<String>) -> OgeomResult<Part> {
    let mut reader = xml::Reader::new(text);
    let mut part = Part {
        scale: 1.0,
        objects: FastMap::default(),
        groups: FastMap::default(),
        textures: FastMap::default(),
        build: Vec::new(),
    };
    // Where the reader is: the element names from the model down.
    let mut stack: Vec<(bool, String)> = Vec::new();
    let mut current: Option<(u32, ObjectDef)> = None;
    let mut group: Option<(u32, Group)> = None;
    while let Some(event) = reader.next()? {
        let start = match event {
            xml::Event::End => {
                if let Some((core, name)) = stack.pop() {
                    match (core, name.as_str()) {
                        (true, "object") => {
                            if let Some((id, def)) = current.take() {
                                part.objects.insert(id, def);
                            }
                        }
                        (_, "basematerials" | "colorgroup" | "texture2dgroup") => {
                            if let Some((id, g)) = group.take() {
                                part.groups.insert(id, g);
                            }
                        }
                        _ => {}
                    }
                }
                continue;
            }
            xml::Event::Start(start) => start,
        };
        let ns = start.namespace.as_ref();
        let core = ns == CORE;
        match (ns, start.local, stack.len()) {
            (CORE, "model", 0) => {
                part.scale = unit_scale(start.get("unit").unwrap_or("millimeter"))?;
                for prefix in start
                    .get("requiredextensions")
                    .unwrap_or("")
                    .split_whitespace()
                {
                    let uri = namespace_of(text, prefix);
                    if uri.contains("securecontent") {
                        ogeom_bail!(
                            Construction,
                            "the model {path} requires the secure-content extension; its meshes are encrypted"
                        );
                    }
                    if !UNDERSTOOD.contains(&uri.as_str()) {
                        warnings.push(format!(
                            "the model {path} requires the extension {uri}, which is not read"
                        ));
                    }
                }
            }
            (_, _, 0) => ogeom_bail!(Construction, "the part {path} is not a 3MF model"),
            (CORE, "object", _) => {
                let id = required(&start, "id")?;
                let object_type = match start.get("type").unwrap_or("model") {
                    "model" => ObjectType::Model,
                    "solidsupport" => ObjectType::SolidSupport,
                    "support" => ObjectType::Support,
                    "surface" => ObjectType::Surface,
                    _ => ObjectType::Other,
                };
                current = Some((
                    id,
                    ObjectDef {
                        name: start.get("name").map(str::to_owned),
                        object_type,
                        pid: number(&start, "pid")?,
                        pindex: number(&start, "pindex")?,
                        content: Content::Components(Vec::new()),
                        metadata: Vec::new(),
                    },
                ));
            }
            (CORE, "mesh", _) => {
                if let Some((_, def)) = current.as_mut() {
                    def.content = Content::Mesh {
                        positions: Vec::new(),
                        triangles: Vec::new(),
                        properties: Vec::new(),
                    };
                }
            }
            (CORE, "vertex", _) => {
                if let Some((_, def)) = current.as_mut()
                    && let Content::Mesh { positions, .. } = &mut def.content
                {
                    positions.push(Point::new(
                        required(&start, "x")?,
                        required(&start, "y")?,
                        required(&start, "z")?,
                    ));
                }
            }
            (CORE, "triangle", _) => {
                if let Some((_, def)) = current.as_mut() {
                    let pid = number(&start, "pid")?.or(def.pid);
                    let first = number(&start, "p1")?;
                    // The second and third corners take the first's property
                    // unless they name their own.
                    let corner = |name: &str, p1: u32| {
                        first
                            .and(start.get(name))
                            .and_then(|t| t.trim().parse().ok())
                            .unwrap_or(p1)
                    };
                    let index = first
                        .or(def.pindex)
                        .map(|p1| [p1, corner("p2", p1), corner("p3", p1)]);
                    if let Content::Mesh {
                        triangles,
                        properties,
                        ..
                    } = &mut def.content
                    {
                        triangles.push([
                            required(&start, "v1")?,
                            required(&start, "v2")?,
                            required(&start, "v3")?,
                        ]);
                        properties.push(pid.zip(index));
                    }
                }
            }
            (CORE, "component", _) => {
                let reference = reference(&start)?;
                if let Some((_, def)) = current.as_mut()
                    && let Content::Components(list) = &mut def.content
                {
                    list.push(reference);
                }
            }
            (CORE, "item", _) => part.build.push(reference(&start)?),
            (CORE, "basematerials", _) | (MATERIAL, "colorgroup", _) => {
                group = Some((required(&start, "id")?, Group::Colours(Vec::new())));
            }
            (CORE, "base", _) => {
                if let Some((_, Group::Colours(colours))) = group.as_mut() {
                    colours.push(parse_colour(
                        start.get("displaycolor").unwrap_or("#FFFFFF"),
                    )?);
                }
            }
            (MATERIAL, "color", _) => {
                if let Some((_, Group::Colours(colours))) = group.as_mut() {
                    colours.push(parse_colour(start.get("color").unwrap_or("#FFFFFF"))?);
                }
            }
            (MATERIAL, "texture2d", _) => {
                let wrap = |name: &str| match start.get(name).unwrap_or("wrap") {
                    "mirror" => Wrap::Mirror,
                    "clamp" | "none" => Wrap::Clamp,
                    _ => Wrap::Repeat,
                };
                if let Some(id) = number::<u32>(&start, "id")? {
                    part.textures.insert(
                        id,
                        TextureDef {
                            path: start.get("path").unwrap_or("").to_owned(),
                            wrap: [wrap("tilestyleu"), wrap("tilestylev")],
                            image: None,
                        },
                    );
                }
            }
            // A texture group that does not hold together is one this
            // reader does not read, not a reason to refuse the package.
            (MATERIAL, "texture2dgroup", _) => {
                if let Some(id) = number::<u32>(&start, "id")? {
                    let texture = start.get("texid").and_then(|t| t.trim().parse().ok());
                    group = Some((
                        id,
                        texture.map_or(Group::Unread, |texture| Group::Texture {
                            texture,
                            coords: Vec::new(),
                        }),
                    ));
                }
            }
            (MATERIAL, "tex2coord", _) => {
                if let Some((_, g @ Group::Texture { .. })) = group.as_mut() {
                    let at = |name: &str| {
                        start
                            .get(name)
                            .and_then(|t| t.trim().parse::<f64>().ok())
                            .filter(|x| x.is_finite())
                    };
                    match (at("u"), at("v"), &mut *g) {
                        (Some(u), Some(v), Group::Texture { coords, .. }) => coords.push((u, v)),
                        _ => *g = Group::Unread,
                    }
                }
            }
            (MATERIAL, _, _) => {
                // Textures, composites and multi-properties: their ids are
                // known, so a triangle naming one is reported, not failed.
                if let Some(id) = number::<u32>(&start, "id")? {
                    part.groups.insert(id, Group::Unread);
                }
                if !start.empty {
                    reader.skip()?;
                    continue;
                }
            }
            (CORE, "metadatagroup", _)
                if current.is_some() && stack.last().is_some_and(|(c, n)| *c && n == "object") => {}
            (CORE, "metadata", _)
                if current.is_some()
                    && stack
                        .last()
                        .is_some_and(|(c, n)| *c && n == "metadatagroup") =>
            {
                let name = start.get("name").unwrap_or("").to_owned();
                let value = reader.text().into_owned();
                if let Some((_, def)) = current.as_mut() {
                    def.metadata.push((name, value));
                }
                if !start.empty {
                    reader.skip()?;
                    continue;
                }
            }
            (CORE, "metadata" | "metadatagroup", _) if !start.empty => {
                reader.skip()?;
                continue;
            }
            _ => {}
        }
        // An empty element still ends, and its end pops it.
        stack.push((core, start.local.to_owned()));
    }
    Ok(part)
}

fn reference(start: &xml::Start<'_>) -> OgeomResult<Reference> {
    Ok(Reference {
        path: start.get_in(PRODUCTION, "path").map(str::to_owned),
        object: required(start, "objectid")?,
        transform: match start.get("transform") {
            Some(text) => Affine::parse(text)?,
            None => Affine::IDENTITY,
        },
    })
}

/// The URI a prefix is declared as on the model element, for naming a
/// required extension.
fn namespace_of(text: &str, prefix: &str) -> String {
    let declaration = format!("xmlns:{prefix}=");
    text.find(&declaration)
        .and_then(|i| {
            let rest = &text[i + declaration.len()..];
            let quote = rest.chars().next()?;
            let rest = &rest[1..];
            rest.find(quote).map(|end| rest[..end].to_owned())
        })
        .unwrap_or_else(|| prefix.to_owned())
}

fn unit_scale(unit: &str) -> OgeomResult<f64> {
    Ok(match unit {
        "micron" => 1e-3,
        "millimeter" => 1.0,
        "centimeter" => 10.0,
        "inch" => 25.4,
        "foot" => 304.8,
        "meter" => 1000.0,
        _ => ogeom_bail!(Construction, "the model unit {unit} is not one 3MF defines"),
    })
}

/// `#RRGGBB` or `#RRGGBBAA`, to RGBA in `[0, 1]`.
fn parse_colour(text: &str) -> OgeomResult<[f64; 4]> {
    let hex = text.trim().trim_start_matches('#');
    let channel = |i: usize| -> Option<f64> {
        let byte = u8::from_str_radix(hex.get(i..i + 2)?, 16).ok()?;
        Some(f64::from(byte) / 255.0)
    };
    match (hex.len(), channel(0), channel(2), channel(4)) {
        (6, Some(r), Some(g), Some(b)) => Ok([r, g, b, 1.0]),
        (8, Some(r), Some(g), Some(b)) => Ok([r, g, b, channel(6).unwrap_or(1.0)]),
        _ => ogeom_bail!(
            Construction,
            "the colour {text} is not #RRGGBB or #RRGGBBAA"
        ),
    }
}

fn object<'p>(parts: &'p FastMap<String, Part>, path: &str, id: u32) -> OgeomResult<&'p ObjectDef> {
    match parts.get(&part_key(path)).and_then(|p| p.objects.get(&id)) {
        Some(def) => Ok(def),
        None => ogeom_bail!(Construction, "the model {path} has no object {id}"),
    }
}

/// A colour a property names: `None` for none, or for a group this reader
/// does not read.
fn resolve_colour(part: Option<&Part>, pid: Option<u32>, index: Option<u32>) -> Option<[f64; 4]> {
    match part?.groups.get(&pid?)? {
        Group::Colours(colours) => colours.get(index? as usize).copied(),
        _ => None,
    }
}

/// The meshes of one build item, gathered as they flatten.
#[derive(Default)]
struct Flat {
    positions: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    /// Every distinct colour a triangle's first corner carries, `None`
    /// among them for a triangle that carries none.
    colours: Vec<Option<[f64; 4]>>,
    /// Each triangle's colour at its three corners, where it has one.
    corners: Vec<Option<[[f64; 4]; 3]>>,
    /// Each triangle's volume, an index into `volumes`.
    volume_of: Vec<usize>,
    volumes: Vec<(Option<String>, ObjectRole)>,
    skipped_properties: bool,
    /// Textures a triangle names that could not be decoded.
    undecoded: bool,
    degenerate: usize,
}

/// Components nest; past this depth a reference is taken to contain itself.
const DEPTH: usize = 64;

/// What the package's configuration parts say of one object of the start
/// part: its own `(key, value)` pairs, and its volumes'.
#[derive(Default)]
struct ConfigObject {
    pairs: Vec<(String, String)>,
    volumes: Vec<ConfigVolume>,
}

/// A volume a configuration part describes: a range of the object's mesh
/// triangles, first and last inclusive, or a component's object id.
#[derive(Default)]
struct ConfigVolume {
    triangles: Option<(usize, usize)>,
    object: Option<u32>,
    pairs: Vec<(String, String)>,
}

/// The role and name a volume carries down from what encloses it.
#[derive(Clone)]
struct Tag {
    name: Option<String>,
    role: ObjectRole,
}

impl Tag {
    /// This tag, overridden by what `pairs` say.
    fn with(&self, pairs: &[(String, String)]) -> Self {
        Self {
            name: value_of(pairs, "name")
                .map(str::to_owned)
                .or_else(|| self.name.clone()),
            role: role_of(pairs).unwrap_or_else(|| self.role.clone()),
        }
    }
}

/// A key without its namespace prefix, lower-cased.
fn bare_key(key: &str) -> String {
    key.rsplit(':').next().unwrap_or(key).to_ascii_lowercase()
}

fn value_of<'a>(pairs: &'a [(String, String)], key: &str) -> Option<&'a str> {
    pairs
        .iter()
        .find(|(k, _)| bare_key(k) == key)
        .map(|(_, v)| v.as_str())
}

/// The role `(key, value)` pairs name, by the vocabulary [`ObjectRole`]
/// lists: a type first, then a flag.
fn role_of(pairs: &[(String, String)]) -> Option<ObjectRole> {
    let squashed = |v: &str| {
        v.chars()
            .filter(|c| !matches!(c, '_' | '-' | ' '))
            .collect::<String>()
            .to_ascii_lowercase()
    };
    for key in ["volume_type", "subtype"] {
        if let Some(value) = value_of(pairs, key) {
            return Some(match squashed(value).as_str() {
                "modelpart" | "normalpart" | "model" | "part" => ObjectRole::Model,
                "parametermodifier" | "modifierpart" | "modifier" => ObjectRole::Modifier,
                "negativevolume" | "negativepart" => ObjectRole::NegativeVolume,
                "supportblocker" => ObjectRole::SupportBlocker,
                "supportenforcer" => ObjectRole::SupportEnforcer,
                _ => ObjectRole::Other(value.to_owned()),
            });
        }
    }
    let set = |key: &str| {
        value_of(pairs, key)
            .is_some_and(|v| matches!(v.trim().to_ascii_lowercase().as_str(), "1" | "true" | "yes"))
    };
    if set("modifier") || set("infill_mesh") || set("cutting_mesh") {
        return Some(ObjectRole::Modifier);
    }
    if set("anti_overhang_mesh") {
        return Some(ObjectRole::SupportBlocker);
    }
    if set("support_mesh") {
        return Some(ObjectRole::SupportEnforcer);
    }
    None
}

/// Read a configuration part: its `object` elements by id, each with the
/// `key`/`value` metadata it holds and the `volume` or `part` elements in
/// it, whose own attributes count as pairs too.
fn parse_config(text: &str) -> OgeomResult<FastMap<u32, ConfigObject>> {
    let mut reader = xml::Reader::new(text);
    let mut objects: FastMap<u32, ConfigObject> = FastMap::default();
    let mut depth = 0_usize;
    // The open object, the depth it opened at, and its open volume.
    let mut object: Option<(u32, usize, ConfigObject)> = None;
    let mut volume: Option<(usize, ConfigVolume)> = None;
    while let Some(event) = reader.next()? {
        let start = match event {
            xml::Event::End => {
                if volume.as_ref().is_some_and(|(d, _)| *d == depth)
                    && let (Some((_, v)), Some((_, _, o))) = (volume.take(), object.as_mut())
                {
                    o.volumes.push(v);
                }
                if object.as_ref().is_some_and(|(_, d, _)| *d == depth)
                    && let Some((id, _, o)) = object.take()
                {
                    let held = objects.entry(id).or_default();
                    held.pairs.extend(o.pairs);
                    held.volumes.extend(o.volumes);
                }
                depth = depth.saturating_sub(1);
                continue;
            }
            xml::Event::Start(start) => start,
        };
        depth += 1;
        let pairs = |start: &xml::Start<'_>| -> Vec<(String, String)> {
            start
                .attributes
                .iter()
                .map(|(_, k, v)| ((*k).to_owned(), v.clone().into_owned()))
                .collect()
        };
        match start.local {
            "object" if object.is_none() && depth == 2 => {
                if let Some(id) = number::<u32>(&start, "id")? {
                    object = Some((id, depth, ConfigObject::default()));
                }
            }
            "volume" | "part" if object.is_some() && volume.is_none() => {
                let first = number::<usize>(&start, "firstid")?;
                let last = number::<usize>(&start, "lastid")?;
                volume = Some((
                    depth,
                    ConfigVolume {
                        triangles: first.zip(last),
                        object: number(&start, "id")?,
                        pairs: pairs(&start),
                    },
                ));
            }
            "metadata" => {
                let key = start.get("key").or_else(|| start.get("name"));
                let value = match start.get("value") {
                    Some(v) => v.to_owned(),
                    None => reader.text().into_owned(),
                };
                if let Some(key) = key {
                    let pair = (key.to_owned(), value);
                    match (volume.as_mut(), object.as_mut()) {
                        (Some((_, v)), _) => v.pairs.push(pair),
                        (None, Some((_, _, o))) => o.pairs.push(pair),
                        (None, None) => {}
                    }
                }
            }
            _ => {}
        }
    }
    Ok(objects)
}

/// What flattening reads besides the parts: the package's configuration of
/// the start part's objects.
struct Context<'p> {
    parts: &'p FastMap<String, Part>,
    start: String,
    config: &'p FastMap<u32, ConfigObject>,
}

impl Context<'_> {
    /// The configuration of object `id` of the part at `path`: only the
    /// start part's objects are configured.
    fn config(&self, path: &str, id: u32) -> Option<&ConfigObject> {
        (part_key(path) == self.start)
            .then(|| self.config.get(&id))
            .flatten()
    }

    /// The tag of the object `id` at `path`, under `inherited`.
    fn tag(&self, path: &str, id: u32, def: &ObjectDef, inherited: &Tag) -> Tag {
        let mut tag = inherited.clone();
        if let Some(config) = self.config(path, id) {
            tag.role = role_of(&config.pairs).unwrap_or(tag.role);
        }
        if let Some(role) = role_of(&def.metadata) {
            tag.role = role;
        }
        tag
    }
}

/// The colour of a triangle at its corners, from the group its property
/// names: `Ok(None)` for no colour, `Err(())` for a group or texture that
/// is not read.
fn corner_colours(
    part: Option<&Part>,
    pid: u32,
    index: [u32; 3],
) -> Result<Option<[[f64; 4]; 3]>, bool> {
    let Some(part) = part else {
        return Ok(None);
    };
    match part.groups.get(&pid) {
        None => Ok(None),
        Some(Group::Unread) => Err(false),
        Some(Group::Colours(colours)) => Ok(index
            .iter()
            .map(|&i| colours.get(i as usize).copied())
            .collect::<Option<Vec<_>>>()
            .map(|c| [c[0], c[1], c[2]])),
        Some(Group::Texture { texture, coords }) => {
            let Some(def) = part.textures.get(texture) else {
                return Ok(None);
            };
            let Some(image) = def.image.as_ref() else {
                return Err(true);
            };
            Ok(index
                .iter()
                .map(|&i| {
                    coords
                        .get(i as usize)
                        .map(|&(u, v)| image.sample(u, v, def.wrap))
                })
                .collect::<Option<Vec<_>>>()
                .map(|c| [c[0], c[1], c[2]]))
        }
    }
}

fn flatten(
    context: &Context<'_>,
    path: &str,
    id: u32,
    placed: Affine,
    flat: &mut Flat,
    depth: usize,
    inherited: &Tag,
) -> OgeomResult<()> {
    if depth > DEPTH {
        ogeom_bail!(Construction, "the object {id} in {path} contains itself");
    }
    let parts = context.parts;
    let def = object(parts, path, id)?;
    let part = parts.get(&part_key(path));
    let tag = context.tag(path, id, def, inherited);
    let config = context.config(path, id);
    match &def.content {
        Content::Components(components) => {
            for c in components {
                let inner = c.path.as_deref().unwrap_or(path);
                let child = config
                    .and_then(|cfg| cfg.volumes.iter().find(|v| v.object == Some(c.object)))
                    .map_or_else(|| tag.clone(), |v| tag.with(&v.pairs));
                flatten(
                    context,
                    inner,
                    c.object,
                    placed.after(c.transform),
                    flat,
                    depth + 1,
                    &child,
                )?;
            }
        }
        Content::Mesh {
            positions,
            triangles,
            properties,
        } => {
            let base = u32::try_from(flat.positions.len()).unwrap_or(u32::MAX);
            let count = positions.len();
            flat.positions
                .extend(positions.iter().map(|p| placed.point(*p)));
            let mirrored = placed.mirrors();
            let ranges: Vec<(usize, usize, Tag)> = config
                .map(|cfg| {
                    cfg.volumes
                        .iter()
                        .filter_map(|v| v.triangles.map(|(a, b)| (a, b, tag.with(&v.pairs))))
                        .collect()
                })
                .unwrap_or_default();
            // The volume the last triangle went to, by range, in this mesh.
            let mut open: Option<Option<usize>> = None;
            for (k, (triangle, property)) in triangles.iter().zip(properties).enumerate() {
                if triangle.iter().any(|&v| v as usize >= count) {
                    ogeom_bail!(
                        Construction,
                        "a triangle of object {id} in {path} names a vertex past the {count} it has"
                    );
                }
                let [a, b, c] = *triangle;
                if a == b || b == c || c == a {
                    flat.degenerate += 1;
                    continue;
                }
                let triangle = if mirrored { [a, c, b] } else { [a, b, c] };
                flat.triangles.push(triangle.map(|v| base + v));
                let range = ranges.iter().position(|(a, b, _)| (*a..=*b).contains(&k));
                if open != Some(range) {
                    let volume = range.map_or(&tag, |r| &ranges[r].2);
                    flat.volumes.push((
                        volume.name.clone().or_else(|| def.name.clone()),
                        volume.role.clone(),
                    ));
                    open = Some(range);
                }
                flat.volume_of.push(flat.volumes.len() - 1);
                let colour = match property {
                    Some((pid, index)) => match part.and_then(|p| p.groups.get(pid)) {
                        Some(Group::Colours(group)) => group.get(index[0] as usize).copied(),
                        Some(Group::Unread) => {
                            flat.skipped_properties = true;
                            None
                        }
                        _ => None,
                    },
                    None => None,
                };
                if !flat.colours.contains(&colour) {
                    flat.colours.push(colour);
                }
                let corners = match property {
                    Some((pid, index)) => match corner_colours(part, *pid, *index) {
                        Ok(c) => c,
                        Err(undecoded) => {
                            flat.undecoded |= undecoded;
                            None
                        }
                    },
                    None => None,
                };
                // A mirrored triangle is rewound; its corners go with it.
                flat.corners
                    .push(corners.map(|[p, q, r]| if mirrored { [p, r, q] } else { [p, q, r] }));
            }
        }
    }
    Ok(())
}

/// A welded build item: the mesh, each kept triangle's corner colours, and
/// its volumes as runs of the kept triangles.
struct Welded {
    mesh: Triangulation,
    corners: Vec<Option<[[f64; 4]; 3]>>,
    volumes: Vec<ThreeMfVolume>,
}

impl Flat {
    /// The mesh, scaled to millimetres, coincident vertices welded, and a
    /// normal at each vertex from the triangles around it.
    fn into_mesh(self, scale: f64, weld: f64) -> Welded {
        let positions: Vec<Point> = self
            .positions
            .iter()
            .map(|p| Point::ORIGIN + (*p - Point::ORIGIN) * scale)
            .collect();
        // Weld on a grid of the weld distance, looking in the neighbouring
        // cells too, so two points either side of a cell wall still meet.
        // A cast saturates, so a coordinate past `i64`'s cells shares the
        // last one and is still compared by distance.
        #[allow(clippy::cast_possible_truncation, reason = "saturating")]
        let cell = |p: Point| {
            (
                (p.x / weld).floor() as i64,
                (p.y / weld).floor() as i64,
                (p.z / weld).floor() as i64,
            )
        };
        let mut grid: FastMap<(i64, i64, i64), Vec<u32>> = FastMap::default();
        let mut kept: Vec<Point> = Vec::with_capacity(positions.len());
        let mut remap = Vec::with_capacity(positions.len());
        for p in &positions {
            let (x, y, z) = cell(*p);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&(
                            x.saturating_add(dx),
                            y.saturating_add(dy),
                            z.saturating_add(dz),
                        )) && let Some(&k) = list
                            .iter()
                            .find(|&&k| kept[k as usize].distance(*p) <= weld)
                        {
                            found = Some(k);
                            break 'search;
                        }
                    }
                }
            }
            let index = found.unwrap_or_else(|| {
                let k = u32::try_from(kept.len()).unwrap_or(u32::MAX);
                kept.push(*p);
                grid.entry((x, y, z)).or_default().push(k);
                k
            });
            remap.push(index);
        }
        let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(self.triangles.len());
        let mut corners = Vec::with_capacity(self.triangles.len());
        let mut volume_of = Vec::with_capacity(self.triangles.len());
        for (k, t) in self.triangles.iter().enumerate() {
            let [a, b, c] = t.map(|v| remap[v as usize]);
            if a != b && b != c && c != a {
                triangles.push([a, b, c]);
                corners.push(self.corners[k]);
                volume_of.push(self.volume_of[k]);
            }
        }
        let mut cursor = 0;
        let volumes = self
            .volumes
            .into_iter()
            .enumerate()
            .map(|(v, (name, role))| {
                let start = cursor;
                while volume_of.get(cursor) == Some(&v) {
                    cursor += 1;
                }
                ThreeMfVolume {
                    name,
                    role,
                    triangles: start..cursor,
                }
            })
            .collect();
        let mut normals = vec![Vector::ZERO; kept.len()];
        for triangle in &triangles {
            let [a, b, c] = triangle.map(|i| kept[i as usize]);
            let face = (b - a).cross(c - a);
            for &i in triangle {
                normals[i as usize] += face;
            }
        }
        let normals = normals
            .into_iter()
            .map(|n| {
                let m = n.magnitude();
                if m > 0.0 { n / m } else { Vector::Z }
            })
            .collect();
        Welded {
            mesh: Triangulation {
                parameters: vec![(0.0, 0.0); kept.len()],
                positions: kept,
                normals,
                triangles,
                deflection_met: true,
            },
            corners,
            volumes,
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::cast_possible_truncation,
        reason = "test code"
    )]
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    const TYPES: &str = r#"<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"/>"#;
    const RELS: &str = r#"<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="r0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>"#;

    /// The vertices and twelve outward triangles of a unit cube at `x`,
    /// each triangle's extra attributes from `property(t)`.
    fn cube(x: f64, property: impl Fn(usize) -> String) -> String {
        let mut out = String::from("<mesh><vertices>");
        for k in 0..8 {
            let (dx, dy, dz) = (k & 1, (k >> 1) & 1, (k >> 2) & 1);
            out.push_str(&format!(
                r#"<vertex x="{}" y="{dy}" z="{dz}"/>"#,
                x + f64::from(dx)
            ));
        }
        out.push_str("</vertices><triangles>");
        let faces = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        for (f, [a, b, c, d]) in faces.iter().enumerate() {
            for (k, [p, q, r]) in [[a, b, c], [a, c, d]].iter().enumerate() {
                out.push_str(&format!(
                    r#"<triangle v1="{p}" v2="{q}" v3="{r}" {}/>"#,
                    property(2 * f + k)
                ));
            }
        }
        out.push_str("</triangles></mesh>");
        out
    }

    fn model(resources: &str, build: &str) -> String {
        format!(
            r#"<?xml version="1.0"?><model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:m="http://schemas.microsoft.com/3dmanufacturing/material/2015/02" xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06"><resources>{resources}</resources><build>{build}</build></model>"#
        )
    }

    fn package(parts: &[(&str, &[u8])]) -> Vec<u8> {
        let mut all: Vec<(&str, &[u8])> = vec![
            ("[Content_Types].xml", TYPES.as_bytes()),
            ("_rels/.rels", RELS.as_bytes()),
        ];
        all.extend_from_slice(parts);
        archive(&all)
    }

    /// One mesh of two cubes, the second a modifier by the triangle range
    /// the configuration part gives it.
    #[test]
    fn a_modifier_volume_reports_its_role_by_triangle_range() {
        let mut mesh = cube(0.0, |_| String::new());
        let second = cube(5.0, |_| String::new());
        // Merge the second cube's vertices and triangles into the first.
        let (v2, t2) = second
            .trim_start_matches("<mesh><vertices>")
            .trim_end_matches("</triangles></mesh>")
            .split_once("</vertices><triangles>")
            .unwrap();
        let t2 = (0..8).rev().fold(t2.to_owned(), |t, k| {
            t.replace(&format!("\"{k}\""), &format!("\"#{}\"", k + 8))
        });
        let t2 = t2.replace("\"#", "\"");
        mesh = mesh.replace("</vertices>", &format!("{v2}</vertices>"));
        mesh = mesh.replace("</triangles>", &format!("{t2}</triangles>"));
        let model = model(
            &format!(r#"<object id="1" type="model" name="bracket">{mesh}</object>"#),
            r#"<item objectid="1"/>"#,
        );
        let config = r#"<?xml version="1.0"?><config><object id="1"><metadata type="object" key="name" value="bracket"/>
<volume firstid="0" lastid="11"><metadata type="volume" key="name" value="body"/><metadata type="volume" key="volume_type" value="ModelPart"/></volume>
<volume firstid="12" lastid="23"><metadata type="volume" key="name" value="denser"/><metadata type="volume" key="volume_type" value="ParameterModifier"/></volume>
</object></config>"#;
        let bytes = package(&[
            ("3D/3dmodel.model", model.as_bytes()),
            ("Metadata/model.config", config.as_bytes()),
        ]);
        let import = read_3mf(&bytes, T).unwrap();
        let object = &import.objects[0];
        assert_eq!(object.role, ObjectRole::Model);
        assert_eq!(
            object.volumes,
            vec![
                ThreeMfVolume {
                    name: Some("body".into()),
                    role: ObjectRole::Model,
                    triangles: 0..12,
                },
                ThreeMfVolume {
                    name: Some("denser".into()),
                    role: ObjectRole::Modifier,
                    triangles: 12..24,
                },
            ]
        );
        assert!(
            object.mesh.triangles[12..]
                .iter()
                .flatten()
                .all(|&v| object.mesh.positions[v as usize].x >= 5.0)
        );
        assert!(
            import.warnings.iter().any(|w| w.contains("a modifier")),
            "{:?}",
            import.warnings
        );
        // Without the configuration, the same mesh is one volume, a part.
        let plain = package(&[("3D/3dmodel.model", model.as_bytes())]);
        let import = read_3mf(&plain, T).unwrap();
        assert_eq!(import.objects[0].volumes.len(), 1);
        assert_eq!(import.objects[0].volumes[0].triangles, 0..24);
        assert!(import.warnings.is_empty(), "{:?}", import.warnings);
    }

    /// An assembly whose components live in a part of their own, the
    /// configuration naming each component by its object id and the
    /// subtype attribute giving its role; and an object whose own metadata
    /// makes it a support enforcer.
    #[test]
    fn component_and_object_roles_read_by_name() {
        let parts = model(
            &format!(
                r#"<object id="2" type="model">{}</object><object id="3" type="model">{}</object>"#,
                cube(0.0, |_| String::new()),
                cube(0.25, |_| String::new()),
            ),
            "",
        );
        let root = model(
            &format!(
                r#"<object id="1" type="model" name="assembly"><components><component p:path="/3D/Objects/parts.model" objectid="2"/><component p:path="/3D/Objects/parts.model" objectid="3"/></components></object>
<object id="4" type="model" name="tower"><metadatagroup><metadata name="tool:support_mesh" preserve="1">True</metadata></metadatagroup>{}</object>"#,
                cube(9.0, |_| String::new())
            ),
            r#"<item objectid="1"/><item objectid="4"/>"#,
        );
        let config = r#"<?xml version="1.0"?><config><object id="1"><metadata key="name" value="assembly"/>
<part id="2" subtype="normal_part"><metadata key="name" value="housing"/></part>
<part id="3" subtype="negative_part"><metadata key="name" value="bore"/></part></object>
<plate><metadata key="plater_id" value="1"/></plate></config>"#;
        let bytes = package(&[
            ("3D/3dmodel.model", root.as_bytes()),
            ("3D/Objects/parts.model", parts.as_bytes()),
            ("Metadata/model_settings.config", config.as_bytes()),
        ]);
        let import = read_3mf(&bytes, T).unwrap();
        let (assembly, tower) = (&import.objects[0], &import.objects[1]);
        let roles: Vec<_> = assembly
            .volumes
            .iter()
            .map(|v| (v.name.as_deref(), v.role.clone(), v.triangles.clone()))
            .collect();
        assert_eq!(
            roles,
            vec![
                (Some("housing"), ObjectRole::Model, 0..12),
                (Some("bore"), ObjectRole::NegativeVolume, 12..24),
            ]
        );
        assert_eq!(tower.role, ObjectRole::SupportEnforcer);
        assert_eq!(tower.volumes[0].role, ObjectRole::SupportEnforcer);
        assert!(
            import
                .warnings
                .iter()
                .any(|w| w.contains("tower is a support enforcer"))
        );
        // An unknown role is kept as written.
        assert_eq!(
            role_of(&[("subtype".into(), "fuzzy_skin_part".into())]),
            Some(ObjectRole::Other("fuzzy_skin_part".into()))
        );
        assert_eq!(
            role_of(&[("modifier".into(), "1".into())]),
            Some(ObjectRole::Modifier)
        );
        assert_eq!(role_of(&[("name".into(), "x".into())]), None);
    }

    /// A colour group named at each corner, and a texture group sampled at
    /// each corner, come back as the colours the file gives.
    #[test]
    fn colour_groups_and_textures_read_as_mesh_colours() {
        let texels = [[255, 0, 0, 255], [0, 0, 255, 255]];
        let png = crate::image::encode_png(2, 1, &texels);
        let resources = format!(
            r##"<m:colorgroup id="5"><m:color color="#FF0000"/><m:color color="#00FF0080"/></m:colorgroup>
<m:texture2d id="6" path="/3D/Textures/two.png" contenttype="image/png" tilestyleu="clamp"/>
<m:texture2dgroup id="7" texid="6"><m:tex2coord u="0.25" v="0.5"/><m:tex2coord u="0.75" v="0.5"/></m:texture2dgroup>
<object id="1" type="model" name="painted">{}</object><object id="2" type="model" name="textured">{}</object>"##,
            cube(0.0, |t| format!(r#"pid="5" p1="{}" p2="1" p3="0""#, t % 2)),
            cube(3.0, |t| format!(r#"pid="7" p1="{}""#, t / 6)),
        );
        let model = model(&resources, r#"<item objectid="1"/><item objectid="2"/>"#);
        let bytes = package(&[
            ("3D/3dmodel.model", model.as_bytes()),
            ("3D/Textures/two.png", &png),
        ]);
        let import = read_3mf(&bytes, T).unwrap();
        assert!(import.warnings.is_empty(), "{:?}", import.warnings);
        let red = [1.0, 0.0, 0.0, 1.0];
        let green = [0.0, 1.0, 0.0, 128.0 / 255.0];
        let painted = &import.objects[0];
        let Some(MeshColours::PerCorner(corners)) = &painted.colours else {
            panic!("{:?}", painted.colours);
        };
        assert_eq!(corners[0], [red, green, red]);
        assert_eq!(corners[1], [green, green, red]);
        let textured = &import.objects[1];
        let Some(MeshColours::PerTriangle(faces)) = &textured.colours else {
            panic!("{:?}", textured.colours);
        };
        assert!(faces[..6].iter().all(|c| *c == red));
        assert!(faces[6..].iter().all(|c| *c == [0.0, 0.0, 1.0, 1.0]));
        assert_eq!(textured.colour, None);

        // The texture part missing, the mesh still reads, uncoloured, and
        // says why.
        let bytes = package(&[("3D/3dmodel.model", model.as_bytes())]);
        let import = read_3mf(&bytes, T).unwrap();
        assert!(import.objects[1].colours.is_none());
        assert!(
            import
                .warnings
                .iter()
                .any(|w| w.contains("no texture part"))
        );
    }

    /// Cut short or damaged anywhere, a package is refused or read, never
    /// a panic; a malformed configuration part is only a warning.
    #[test]
    fn hostile_packages_never_panic() {
        let model = model(
            &format!(
                r#"<object id="1" type="model">{}</object>"#,
                cube(0.0, |_| String::new())
            ),
            r#"<item objectid="1"/>"#,
        );
        let bytes = package(&[
            ("3D/3dmodel.model", model.as_bytes()),
            (
                "Metadata/a.config",
                b"<config><object id=\"x\"><volume firstid=\"q\"/>",
            ),
        ]);
        let import = read_3mf(&bytes, T).unwrap();
        assert!(
            import
                .warnings
                .iter()
                .any(|w| w.contains("configuration part"))
        );
        assert_eq!(import.objects[0].mesh.triangles.len(), 12);
        for cut in 0..bytes.len() {
            assert!(read_3mf(&bytes[..cut], T).is_err(), "cut at {cut}");
        }
        let mut seed = 0x9e37_79b9_u32;
        for _ in 0..2000 {
            let mut noisy = bytes.clone();
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            let at = seed as usize % noisy.len();
            noisy[at] ^= (seed >> 8) as u8 | 1;
            let _ = read_3mf(&noisy, T);
        }
        for garbage in [&b""[..], b"PK\x05\x06", b"\x89PNG\r\n\x1a\n"] {
            assert!(read_3mf(garbage, T).is_err());
        }
    }
}
