//! Writing tessellated interchange: glTF 2.0, OBJ and PLY.
//!
//! Three formats, one philosophy carried over from the DXF writer: the
//! functions take bare tessellations rather than a model, so anything that
//! produces a [`Triangulation`] — a solid, one face, a healed import —
//! exports without this module knowing where it came from. glTF is written
//! as GLB, the single-file binary form, with positions, normals, indices
//! and an optional base colour per mesh; OBJ and PLY are the plain-text
//! dialects every downstream tool reads.

use ogeom_topo::Triangulation;
use std::fmt::Write as _;

/// One mesh to export, with what the format can say about it.
#[derive(Debug, Clone)]
pub struct ExportMesh<'a> {
    /// The tessellation.
    pub mesh: &'a Triangulation,
    /// An RGBA base colour in `[0, 1]`, where the format carries one.
    pub colour: Option<[f64; 4]>,
    /// A name, where the format carries one.
    pub name: Option<String>,
}

impl<'a> ExportMesh<'a> {
    /// A bare mesh, no colour, no name.
    #[must_use]
    pub fn plain(mesh: &'a Triangulation) -> Self {
        Self {
            mesh,
            colour: None,
            name: None,
        }
    }
}

// --- glTF 2.0 (GLB) ----------------------------------------------------------

/// Write meshes as a GLB — glTF 2.0's single-file binary form.
///
/// One buffer, one node per mesh under one scene; positions and normals as
/// `f32` vectors, indices as `u32`, and a metallic–roughness material with
/// the base colour where one was given. Empty meshes are skipped — a node
/// with nothing to draw is not something a viewer should be handed.
#[must_use]
pub fn write_glb(meshes: &[ExportMesh<'_>]) -> Vec<u8> {
    let mut binary: Vec<u8> = Vec::new();
    let mut accessors = String::new();
    let mut buffer_views = String::new();
    let mut mesh_json = String::new();
    let mut node_json = String::new();
    let mut material_json = String::new();
    let mut accessor_count = 0_usize;
    let mut view_count = 0_usize;
    let mut mesh_count = 0_usize;
    let mut material_count = 0_usize;

    for export in meshes {
        let mesh = export.mesh;
        if mesh.triangles.is_empty() || mesh.positions.is_empty() {
            continue;
        }
        let comma = |s: &mut String| {
            if !s.is_empty() {
                s.push(',');
            }
        };

        // Positions.
        let position_offset = binary.len();
        let (mut low, mut high) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &mesh.positions {
            for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                low[k] = low[k].min(v);
                high[k] = high[k].max(v);
                #[allow(clippy::cast_possible_truncation)]
                binary.extend_from_slice(&(v as f32).to_le_bytes());
            }
        }
        comma(&mut buffer_views);
        let _ = write!(
            buffer_views,
            r#"{{"buffer":0,"byteOffset":{position_offset},"byteLength":{}}}"#,
            binary.len() - position_offset
        );
        let position_view = view_count;
        view_count += 1;
        comma(&mut accessors);
        #[allow(clippy::cast_possible_truncation)]
        let (lo, hi) = (low.map(|v| v as f32), high.map(|v| v as f32));
        let _ = write!(
            accessors,
            r#"{{"bufferView":{position_view},"componentType":5126,"count":{},"type":"VEC3","min":[{},{},{}],"max":[{},{},{}]}}"#,
            mesh.positions.len(),
            lo[0],
            lo[1],
            lo[2],
            hi[0],
            hi[1],
            hi[2],
        );
        let position_accessor = accessor_count;
        accessor_count += 1;

        // Normals, normalized as glTF requires.
        let normal_offset = binary.len();
        for n in &mesh.normals {
            let magnitude = n.magnitude();
            let unit = if magnitude > 0.0 {
                [n.x / magnitude, n.y / magnitude, n.z / magnitude]
            } else {
                [0.0, 0.0, 1.0]
            };
            for v in unit {
                #[allow(clippy::cast_possible_truncation)]
                binary.extend_from_slice(&(v as f32).to_le_bytes());
            }
        }
        comma(&mut buffer_views);
        let _ = write!(
            buffer_views,
            r#"{{"buffer":0,"byteOffset":{normal_offset},"byteLength":{}}}"#,
            binary.len() - normal_offset
        );
        let normal_view = view_count;
        view_count += 1;
        comma(&mut accessors);
        let _ = write!(
            accessors,
            r#"{{"bufferView":{normal_view},"componentType":5126,"count":{},"type":"VEC3"}}"#,
            mesh.normals.len(),
        );
        let normal_accessor = accessor_count;
        accessor_count += 1;

        // Indices.
        let index_offset = binary.len();
        for t in &mesh.triangles {
            for &k in t {
                binary.extend_from_slice(&k.to_le_bytes());
            }
        }
        comma(&mut buffer_views);
        let _ = write!(
            buffer_views,
            r#"{{"buffer":0,"byteOffset":{index_offset},"byteLength":{}}}"#,
            binary.len() - index_offset
        );
        let index_view = view_count;
        view_count += 1;
        comma(&mut accessors);
        let _ = write!(
            accessors,
            r#"{{"bufferView":{index_view},"componentType":5125,"count":{},"type":"SCALAR"}}"#,
            mesh.triangles.len() * 3,
        );
        let index_accessor = accessor_count;
        accessor_count += 1;

        // Material, when a colour was given.
        let material = export.colour.map(|[r, g, b, a]| {
            comma(&mut material_json);
            let _ = write!(
                material_json,
                r#"{{"pbrMetallicRoughness":{{"baseColorFactor":[{r},{g},{b},{a}],"metallicFactor":0.1,"roughnessFactor":0.8}}}}"#,
            );
            material_count += 1;
            material_count - 1
        });

        comma(&mut mesh_json);
        let material_field = material.map_or(String::new(), |m| format!(r#","material":{m}"#));
        let _ = write!(
            mesh_json,
            r#"{{"primitives":[{{"attributes":{{"POSITION":{position_accessor},"NORMAL":{normal_accessor}}},"indices":{index_accessor}{material_field}}}]}}"#,
        );
        comma(&mut node_json);
        let name_field = export.name.as_ref().map_or(String::new(), |n| {
            format!(r#","name":"{}""#, escape_json(n))
        });
        let _ = write!(node_json, r#"{{"mesh":{mesh_count}{name_field}}}"#);
        mesh_count += 1;
    }

    let node_indices: Vec<String> = (0..mesh_count).map(|i| i.to_string()).collect();
    let json = format!(
        r#"{{"asset":{{"version":"2.0","generator":"ogeom"}},"scene":0,"scenes":[{{"nodes":[{}]}}],"nodes":[{node_json}],"meshes":[{mesh_json}],"materials":[{material_json}],"accessors":[{accessors}],"bufferViews":[{buffer_views}],"buffers":[{{"byteLength":{}}}]}}"#,
        node_indices.join(","),
        binary.len(),
    );
    // A colour-free file carries no materials array worth having.
    let json = json.replace(r#","materials":[],"#, ",");

    // GLB framing: 4-byte alignment, JSON padded with spaces, binary with
    // zeros.
    let mut json_bytes = json.into_bytes();
    while !json_bytes.len().is_multiple_of(4) {
        json_bytes.push(b' ');
    }
    while !binary.len().is_multiple_of(4) {
        binary.push(0);
    }
    let total = 12 + 8 + json_bytes.len() + 8 + binary.len();
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(&0x4654_6C67_u32.to_le_bytes()); // "glTF"
    out.extend_from_slice(&2_u32.to_le_bytes());
    out.extend_from_slice(&u32::try_from(total).unwrap_or(u32::MAX).to_le_bytes());
    out.extend_from_slice(
        &u32::try_from(json_bytes.len())
            .unwrap_or(u32::MAX)
            .to_le_bytes(),
    );
    out.extend_from_slice(&0x4E4F_534A_u32.to_le_bytes()); // "JSON"
    out.extend_from_slice(&json_bytes);
    out.extend_from_slice(
        &u32::try_from(binary.len())
            .unwrap_or(u32::MAX)
            .to_le_bytes(),
    );
    out.extend_from_slice(&0x004E_4942_u32.to_le_bytes()); // "BIN\0"
    out.extend_from_slice(&binary);
    out
}

fn escape_json(s: &str) -> String {
    s.chars()
        .flat_map(|c| match c {
            '"' => vec!['\\', '"'],
            '\\' => vec!['\\', '\\'],
            c if c.is_control() => format!("\\u{:04x}", c as u32).chars().collect(),
            c => vec![c],
        })
        .collect()
}

// --- OBJ ---------------------------------------------------------------------

/// Write meshes as Wavefront OBJ.
///
/// `v`/`vn`/`f` records, one `o` group per named mesh, indices one-based
/// and shared across the file as the format demands.
#[must_use]
pub fn write_obj(meshes: &[ExportMesh<'_>]) -> String {
    let mut out = String::from("# ogeom\n");
    let mut base = 1_usize;
    for (i, export) in meshes.iter().enumerate() {
        let mesh = export.mesh;
        if mesh.triangles.is_empty() {
            continue;
        }
        let name = export.name.clone().unwrap_or_else(|| format!("mesh-{i}"));
        let _ = writeln!(out, "o {name}");
        for p in &mesh.positions {
            let _ = writeln!(out, "v {} {} {}", p.x, p.y, p.z);
        }
        for n in &mesh.normals {
            let _ = writeln!(out, "vn {} {} {}", n.x, n.y, n.z);
        }
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|k| k as usize + base);
            let _ = writeln!(out, "f {a}//{a} {b}//{b} {c}//{c}");
        }
        base += mesh.positions.len();
    }
    out
}

// --- PLY ---------------------------------------------------------------------

/// Write one mesh as ASCII PLY.
///
/// Vertices with normals, faces as index lists; a colour, when given, is
/// carried per vertex as the `uchar` triple the format convention expects.
#[must_use]
pub fn write_ply(export: &ExportMesh<'_>) -> String {
    let mesh = export.mesh;
    let colour = export.colour.map(|[r, g, b, _]| {
        [
            (r.clamp(0.0, 1.0) * 255.0).round(),
            (g.clamp(0.0, 1.0) * 255.0).round(),
            (b.clamp(0.0, 1.0) * 255.0).round(),
        ]
    });
    let mut out = String::from("ply\nformat ascii 1.0\ncomment ogeom\n");
    let _ = writeln!(out, "element vertex {}", mesh.positions.len());
    out.push_str(
        "property float x\nproperty float y\nproperty float z\n\
         property float nx\nproperty float ny\nproperty float nz\n",
    );
    if colour.is_some() {
        out.push_str("property uchar red\nproperty uchar green\nproperty uchar blue\n");
    }
    let _ = writeln!(out, "element face {}", mesh.triangles.len());
    out.push_str("property list uchar uint vertex_indices\nend_header\n");
    for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
        let _ = write!(out, "{} {} {} {} {} {}", p.x, p.y, p.z, n.x, n.y, n.z);
        if let Some([r, g, b]) = colour {
            let _ = write!(out, " {r} {g} {b}");
        }
        out.push('\n');
    }
    for t in &mesh.triangles {
        let _ = writeln!(out, "3 {} {} {}", t[0], t[1], t[2]);
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_math::Frame;
    use ogeom_topo::Model;

    const T: Tolerances = Tolerances::millimetres();

    fn box_mesh() -> Triangulation {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
            .unwrap()
            .shape;
        ogeom_mesh::triangulate(&model, &solid, ogeom_mesh::Deflection::default(), T).unwrap()
    }

    #[test]
    fn a_glb_frames_its_chunks_and_counts_its_geometry() {
        let mesh = box_mesh();
        let glb = write_glb(&[ExportMesh {
            mesh: &mesh,
            colour: Some([0.8, 0.2, 0.1, 1.0]),
            name: Some("box".into()),
        }]);
        // Magic, version, and a total length that matches the file.
        assert_eq!(&glb[0..4], b"glTF");
        assert_eq!(u32::from_le_bytes(glb[4..8].try_into().unwrap()), 2);
        assert_eq!(
            u32::from_le_bytes(glb[8..12].try_into().unwrap()) as usize,
            glb.len()
        );
        // The JSON chunk parses far enough to carry the structure.
        let json_length = u32::from_le_bytes(glb[12..16].try_into().unwrap()) as usize;
        assert_eq!(&glb[16..20], b"JSON");
        let json = core::str::from_utf8(&glb[20..20 + json_length]).unwrap();
        assert!(json.contains(r#""version":"2.0""#));
        assert!(json.contains(r#""POSITION":0"#));
        assert!(json.contains(r#""indices":2"#));
        assert!(json.contains("baseColorFactor"));
        assert!(json.contains(r#""name":"box""#));
        assert!(json.contains(r#""min":["#), "POSITION carries bounds");
        // The binary chunk holds what the accessors promise: positions,
        // normals, indices.
        let expected =
            mesh.positions.len() * 12 + mesh.normals.len() * 12 + mesh.triangles.len() * 12;
        let padded = expected + (4 - expected % 4) % 4;
        let binary_length =
            u32::from_le_bytes(glb[20 + json_length..24 + json_length].try_into().unwrap());
        assert_eq!(binary_length as usize, padded);
    }

    #[test]
    fn an_obj_counts_its_records_and_shares_its_index_space() {
        let mesh = box_mesh();
        let one = ExportMesh::plain(&mesh);
        let text = write_obj(&[one.clone(), one]);
        assert_eq!(
            text.lines().filter(|l| l.starts_with("v ")).count(),
            mesh.positions.len() * 2
        );
        assert_eq!(
            text.lines().filter(|l| l.starts_with("f ")).count(),
            mesh.triangles.len() * 2
        );
        // The second mesh's faces reference the second mesh's vertices.
        let last_face = text.lines().rev().find(|l| l.starts_with("f ")).unwrap();
        let first_index: usize = last_face
            .split_whitespace()
            .nth(1)
            .unwrap()
            .split('/')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert!(first_index > mesh.positions.len());
    }

    #[test]
    fn a_ply_declares_what_it_carries() {
        let mesh = box_mesh();
        let text = write_ply(&ExportMesh {
            mesh: &mesh,
            colour: Some([0.0, 0.5, 1.0, 1.0]),
            name: None,
        });
        assert!(text.starts_with("ply\nformat ascii 1.0\n"));
        assert!(text.contains(&format!("element vertex {}", mesh.positions.len())));
        assert!(text.contains(&format!("element face {}", mesh.triangles.len())));
        assert!(text.contains("property uchar red"));
        let body_faces = text
            .lines()
            .skip_while(|l| *l != "end_header")
            .filter(|l| l.starts_with("3 "))
            .count();
        assert_eq!(body_faces, mesh.triangles.len());
    }
}
