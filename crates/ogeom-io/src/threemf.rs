//! 3MF: the mesh package, written with its own archive.
//!
//! A 3MF file is a ZIP holding three parts: the content-types declaration,
//! a relationship pointing at the model, and the model itself — an XML
//! document of meshes and the items that place them. The mesh half is
//! ordinary; the archive is the part a kernel usually reaches for a library
//! to do, and this does not.
//!
//! It writes *stored* entries — no compression — which the ZIP format has
//! always allowed and every reader accepts. That is the whole reason the
//! archive can live here in two hundred lines instead of arriving as a
//! dependency with a compressor in it: what a CAD kernel needs from ZIP is
//! the container, not the codec. A file written this way is larger than one
//! deflated, and says so by being what it is.

use ogeom_topo::Triangulation;
use std::fmt::Write as _;

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
/// One object per mesh, one item per object, all in millimetres — the
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

/// A ZIP archive of stored — uncompressed — entries.
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

/// The ZIP checksum: CRC-32, reflected, polynomial `0xEDB8_8320`.
fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFF_u32;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            let carry = crc & 1;
            crc >>= 1;
            if carry != 0 {
                crc ^= 0xEDB8_8320;
            }
        }
    }
    !crc
}

/// Read the stored parts of a package: the entry names and their bytes.
///
/// Only *stored* entries come back, which is what this writes and what a
/// 3MF from a machine that cared about size will not be. An entry that was
/// deflated is reported by name in the error rather than silently missing —
/// a package half-read is worse than one refused.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// bytes are not a ZIP archive, or an entry is compressed.
pub fn read_package(bytes: &[u8]) -> ogeom_core::OgeomResult<Vec<(String, Vec<u8>)>> {
    let mut out = Vec::new();
    let mut at = 0_usize;
    while at + 30 <= bytes.len() {
        let signature =
            u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]]);
        if signature != 0x0403_4b50 {
            break;
        }
        let method = u16::from_le_bytes([bytes[at + 8], bytes[at + 9]]);
        let size = u32::from_le_bytes([
            bytes[at + 18],
            bytes[at + 19],
            bytes[at + 20],
            bytes[at + 21],
        ]) as usize;
        let name_len = u16::from_le_bytes([bytes[at + 26], bytes[at + 27]]) as usize;
        let extra_len = u16::from_le_bytes([bytes[at + 28], bytes[at + 29]]) as usize;
        let name_at = at + 30;
        let data_at = name_at + name_len + extra_len;
        if data_at + size > bytes.len() {
            ogeom_core::ogeom_bail!(
                Construction,
                "an archive entry runs past the end of the file"
            );
        }
        let name = String::from_utf8_lossy(&bytes[name_at..name_at + name_len]).into_owned();
        if method != 0 {
            ogeom_core::ogeom_bail!(
                Construction,
                "the entry {name} is compressed, and this archive reader stores only"
            );
        }
        out.push((name, bytes[data_at..data_at + size].to_vec()));
        at = data_at + size;
    }
    if out.is_empty() {
        ogeom_core::ogeom_bail!(Construction, "these bytes are not a stored ZIP archive");
    }
    Ok(out)
}
