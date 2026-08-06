//! §F3 of `docs/PLAN.md`: glTF read back.
//!
//! What this crate writes it reads, and what *another* writer might have
//! written it also reads: the indirection is the format, so the tests give the
//! reader documents built deliberately in the forms a writer chooses —
//! interleaved with a stride, indices as unsigned shorts, normals as
//! normalized signed bytes, a sparse accessor over a base, a node hierarchy
//! with a matrix and with translation–rotation–scale.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::io::mesh_formats::{ExportMesh, ImportedMesh, read_glb, read_gltf, write_glb};
use ogeom::math::{Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Triangulation};

const T: Tolerances = Tolerances::millimetres();

fn block_mesh() -> Triangulation {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 30.0), T)
        .unwrap()
        .shape;
    ogeom::mesh::triangulate(&model, &block, Deflection::default(), T).unwrap()
}

/// The signed volume a closed triangle soup encloses — the one claim a mesh
/// format really makes.
fn volume(m: &Triangulation) -> f64 {
    m.triangles
        .iter()
        .map(|[a, b, c]| {
            let (a, b, c) = (
                m.positions[*a as usize],
                m.positions[*b as usize],
                m.positions[*c as usize],
            );
            a.to_vector().dot(b.to_vector().cross(c.to_vector())) / 6.0
        })
        .sum::<f64>()
        .abs()
}

/// The round trip: what `write_glb` emits, `read_glb` reads.
#[test]
fn a_glb_comes_back_as_the_meshes_it_was() {
    let mesh = block_mesh();
    let written = write_glb(&[ExportMesh {
        mesh: &mesh,
        colour: Some([0.2, 0.4, 0.6, 1.0]),
        name: Some("block".to_owned()),
    }]);
    let read = read_glb(&written).unwrap();
    assert_eq!(read.len(), 1);
    let back = &read[0];
    assert_eq!(back.name.as_deref(), Some("block"));
    assert_eq!(back.colour, Some([0.2, 0.4, 0.6, 1.0]));
    assert_eq!(back.mesh.positions.len(), mesh.positions.len());
    assert_eq!(back.mesh.triangles.len(), mesh.triangles.len());
    for (a, b) in back.mesh.positions.iter().zip(&mesh.positions) {
        // The file carries `f32`, so the trip is exact to single precision
        // and says so rather than claiming more.
        assert!(a.distance(*b) < 1e-4, "{a:?} against {b:?}");
    }
    for (a, b) in back.mesh.normals.iter().zip(&mesh.normals) {
        assert!(a.cross(*b).magnitude() < 1e-5, "{a:?} against {b:?}");
    }
    assert!(
        (volume(&back.mesh) - 6000.0).abs() < 1e-2,
        "{}",
        volume(&back.mesh)
    );

    // Two meshes come back as two.
    let pair = write_glb(&[ExportMesh::plain(&mesh), ExportMesh::plain(&mesh)]);
    assert_eq!(read_glb(&pair).unwrap().len(), 2);
}

/// A base64 payload as a `.gltf` document's data uri.
fn data_uri(bytes: &[u8]) -> String {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::from("data:application/octet-stream;base64,");
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            chunk.get(1).copied().unwrap_or(0),
            chunk.get(2).copied().unwrap_or(0),
        ];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for k in 0..4 {
            if k <= chunk.len() {
                out.push(char::from(ALPHABET[((n >> (18 - 6 * k)) & 0x3F) as usize]));
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// A tetrahedron written the way a *different* writer would: positions and
/// normals interleaved behind one strided buffer view, indices as unsigned
/// shorts, normals as normalized signed bytes, and the whole thing in a
/// `.gltf` document with a data uri instead of a GLB.
#[test]
fn the_indirection_a_writer_chooses_is_read_as_it_finds_it() {
    // Four corners of a tetrahedron, and a normal per vertex pointing out
    // along its own position.
    let corners = [
        [0.0_f32, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ];
    // Interleaved: twelve bytes of position, then three signed bytes of
    // normal and one of padding — a sixteen-byte stride, which is the shape
    // an aligned writer picks.
    let normals: [[i8; 3]; 4] = [[0, 0, -127], [127, 0, 0], [0, 127, 0], [0, 0, 127]];
    let mut binary: Vec<u8> = Vec::new();
    for (corner, normal) in corners.iter().zip(&normals) {
        for v in corner {
            binary.extend_from_slice(&v.to_le_bytes());
        }
        for v in normal {
            binary.push(*v as u8);
        }
        binary.push(0);
    }
    let index_offset = binary.len();
    // Wound so the tetrahedron is closed and outward.
    let faces: [[u16; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
    for face in faces {
        for k in face {
            binary.extend_from_slice(&k.to_le_bytes());
        }
    }
    while !binary.len().is_multiple_of(4) {
        binary.push(0);
    }

    let document = format!(
        r#"{{
          "asset": {{"version": "2.0"}},
          "scene": 0,
          "scenes": [{{"nodes": [0]}}],
          "nodes": [{{"mesh": 0, "name": "tetra"}}],
          "meshes": [{{"primitives": [{{
              "attributes": {{"POSITION": 0, "NORMAL": 1}},
              "indices": 2, "mode": 4
          }}]}}],
          "accessors": [
            {{"bufferView": 0, "byteOffset": 0, "componentType": 5126, "count": 4, "type": "VEC3"}},
            {{"bufferView": 0, "byteOffset": 12, "componentType": 5120, "count": 4,
              "type": "VEC3", "normalized": true}},
            {{"bufferView": 1, "componentType": 5123, "count": 12, "type": "SCALAR"}}
          ],
          "bufferViews": [
            {{"buffer": 0, "byteOffset": 0, "byteLength": {index_offset}, "byteStride": 16}},
            {{"buffer": 0, "byteOffset": {index_offset}, "byteLength": 24}}
          ],
          "buffers": [{{"byteLength": {}, "uri": "{}"}}]
        }}"#,
        binary.len(),
        data_uri(&binary),
    );

    let read = read_gltf(&document).unwrap();
    assert_eq!(read.len(), 1);
    let mesh = &read[0].mesh;
    assert_eq!(read[0].name.as_deref(), Some("tetra"));
    assert_eq!(mesh.positions.len(), 4);
    assert_eq!(mesh.triangles.len(), 4);
    assert!(mesh.positions[1].distance(Point::new(4.0, 0.0, 0.0)) < 1e-6);
    // The normalized signed bytes read back as unit vectors, not as 127.
    for n in &mesh.normals {
        assert!((n.magnitude() - 1.0).abs() < 1e-6, "{n:?}");
    }
    // A tetrahedron with legs of four: a sixth of the box, which is 32/3.
    assert!(
        (volume(mesh) - 32.0 / 3.0).abs() < 1e-6,
        "the strided, short-indexed, byte-normalled tetrahedron: {}",
        volume(mesh)
    );
}

/// A sparse accessor: a base of zeros with three of its four elements
/// replaced. This is the form a writer uses for a mesh that is mostly one
/// thing, and the reader has to apply the override *after* the base.
#[test]
fn a_sparse_accessor_overrides_the_base_it_sits_on() {
    // The base: four positions, all the origin.
    let mut binary: Vec<u8> = vec![0; 4 * 3 * 4];
    // The sparse indices: elements 1, 2 and 3.
    let sparse_index_offset = binary.len();
    for k in [1_u16, 2, 3] {
        binary.extend_from_slice(&k.to_le_bytes());
    }
    binary.extend_from_slice(&0_u16.to_le_bytes()); // padding to four bytes
    // The sparse values.
    let sparse_value_offset = binary.len();
    for corner in [[4.0_f32, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]] {
        for v in corner {
            binary.extend_from_slice(&v.to_le_bytes());
        }
    }
    let index_offset = binary.len();
    for face in [[0_u8, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]] {
        binary.extend_from_slice(&face);
    }
    while !binary.len().is_multiple_of(4) {
        binary.push(0);
    }

    let document = format!(
        r#"{{
          "asset": {{"version": "2.0"}},
          "scene": 0,
          "scenes": [{{"nodes": [0]}}],
          "nodes": [{{"mesh": 0}}],
          "meshes": [{{"primitives": [{{"attributes": {{"POSITION": 0}}, "indices": 1}}]}}],
          "accessors": [
            {{"bufferView": 0, "componentType": 5126, "count": 4, "type": "VEC3",
              "sparse": {{"count": 3,
                "indices": {{"bufferView": 1, "componentType": 5123}},
                "values": {{"bufferView": 2}}}}}},
            {{"bufferView": 3, "componentType": 5121, "count": 12, "type": "SCALAR"}}
          ],
          "bufferViews": [
            {{"buffer": 0, "byteOffset": 0, "byteLength": 48}},
            {{"buffer": 0, "byteOffset": {sparse_index_offset}, "byteLength": 6}},
            {{"buffer": 0, "byteOffset": {sparse_value_offset}, "byteLength": 36}},
            {{"buffer": 0, "byteOffset": {index_offset}, "byteLength": 12}}
          ],
          "buffers": [{{"byteLength": {}, "uri": "{}"}}]
        }}"#,
        binary.len(),
        data_uri(&binary),
    );

    let read = read_gltf(&document).unwrap();
    let mesh = &read[0].mesh;
    assert_eq!(mesh.positions.len(), 4);
    assert!(mesh.positions[0].distance(Point::ORIGIN) < 1e-9, "the base");
    assert!(mesh.positions[2].distance(Point::new(0.0, 4.0, 0.0)) < 1e-9);
    assert!(
        (volume(mesh) - 32.0 / 3.0).abs() < 1e-6,
        "the same tetrahedron, three quarters of it sparse: {}",
        volume(mesh)
    );
    // No normals were given, so the triangles' own were worked out — and they
    // point outward, which is what the winding says.
    assert_eq!(mesh.normals.len(), 4);
    for n in &mesh.normals {
        assert!((n.magnitude() - 1.0).abs() < 1e-9);
    }
}

/// The node hierarchy places what it holds, stated either way, and an uneven
/// scale carries normals through the inverse transpose rather than through
/// the scale itself.
#[test]
fn nodes_place_their_meshes_and_normals_survive_an_uneven_scale() {
    // A ball, not a box: every normal of a box points along an axis, where a
    // diagonal scale and its inverse transpose agree after normalizing, and
    // the test would prove nothing. A sphere's normals are slanted.
    let mesh = {
        let mut model = Model::new();
        let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 5.0, T)
            .unwrap()
            .shape;
        ogeom::mesh::triangulate(&model, &ball, Deflection::with_chord(0.2).unwrap(), T).unwrap()
    };
    let base = write_glb(&[ExportMesh::plain(&mesh)]);
    let placed: Vec<ImportedMesh> = {
        // Re-frame the GLB's own JSON with a parent node that scales and
        // translates: the geometry is the same, the placement is not.
        let json_length = u32::from_le_bytes([base[12], base[13], base[14], base[15]]) as usize;
        let json = core::str::from_utf8(&base[20..20 + json_length]).unwrap();
        let rewritten = json.replace(r#""nodes":[0]"#, r#""nodes":[1]"#).replace(
            r#""nodes":[{"mesh":0}]"#,
            r#""nodes":[{"mesh":0},{"children":[0],"scale":[2,1,0.5],"translation":[100,0,0]}]"#,
        );
        assert!(rewritten.contains("children"), "the rewrite took: {json}");
        let mut out = base[..12].to_vec();
        let mut bytes = rewritten.into_bytes();
        while !bytes.len().is_multiple_of(4) {
            bytes.push(b' ');
        }
        out.extend_from_slice(&u32::try_from(bytes.len()).unwrap().to_le_bytes());
        out.extend_from_slice(&base[16..20]);
        out.extend_from_slice(&bytes);
        out.extend_from_slice(&base[20 + json_length..]);
        let total = u32::try_from(out.len()).unwrap();
        out[8..12].copy_from_slice(&total.to_le_bytes());
        read_glb(&out).unwrap()
    };
    assert_eq!(placed.len(), 1);
    let back = &placed[0].mesh;
    assert_eq!(back.positions.len(), mesh.positions.len());
    // The ball of radius five at the origin, scaled and moved: an ellipsoid
    // spanning x in [90, 110], y in [-5, 5], z in [-2.5, 2.5].
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &back.positions {
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[k] = lo[k].min(v);
            hi[k] = hi[k].max(v);
        }
    }
    assert!(
        (lo[0] - 90.0).abs() < 1e-2 && (hi[0] - 110.0).abs() < 1e-2,
        "{lo:?} {hi:?}"
    );
    assert!(
        (hi[1] - 5.0).abs() < 1e-2 && (hi[2] - 2.5).abs() < 1e-2,
        "{lo:?} {hi:?}"
    );
    // Volume scales by the determinant, which here is one — so the claim is
    // against the mesh that went in, not against the ball it approximates:
    // the reader changed the placement and nothing else.
    let want = volume(&mesh);
    assert!(
        (volume(back) - want).abs() < want * 1e-3,
        "{} against {want}",
        volume(back)
    );

    // And the normals came through the inverse transpose. For `diag(2, 1, ½)`
    // that is `diag(½, 1, 2)` up to a factor, which normalizing removes — and
    // it is a genuinely different direction from the scale's own, which is
    // what the second assertion holds it to.
    let mut differed = 0;
    for (before, after) in mesh.normals.iter().zip(&back.normals) {
        // A mesh may carry a vertex whose normal is undefined — a sphere's
        // pole is one — and a direction nothing stated is not one to check.
        if (before.magnitude() - 1.0).abs() > 1e-9 {
            continue;
        }
        let want = ogeom::math::Vector::new(before.x / 2.0, before.y, before.z * 2.0);
        let want = want / want.magnitude();
        assert!((after.magnitude() - 1.0).abs() < 1e-6, "unit: {after:?}");
        assert!(
            after.cross(want).magnitude() < 1e-5,
            "the inverse transpose: {after:?} against {want:?}"
        );
        let naive = ogeom::math::Vector::new(before.x * 2.0, before.y, before.z / 2.0);
        let naive = naive / naive.magnitude();
        if naive.cross(want).magnitude() > 1e-3 {
            differed += 1;
        }
    }
    assert!(
        differed > mesh.normals.len() / 4,
        "the scale itself would have been visibly wrong: {differed} of {}",
        mesh.normals.len()
    );
}

/// What the reader does not do, it says. Each refusal names the thing.
#[test]
fn what_is_not_read_is_refused_by_name() {
    let mesh = block_mesh();
    let base = write_glb(&[ExportMesh::plain(&mesh)]);

    // Not a GLB at all.
    let message = format!("{}", read_glb(b"not a glb at all").unwrap_err());
    assert!(message.contains("not a GLB"), "{message}");

    // A version this is not.
    let mut wrong = base.clone();
    wrong[4] = 3;
    let message = format!("{}", read_glb(&wrong).unwrap_err());
    assert!(message.contains("not glTF 2.0"), "{message}");

    // A buffer that points at a file.
    let document = r#"{"asset":{"version":"2.0"},"buffers":[{"byteLength":4,"uri":"scene.bin"}]}"#;
    let message = format!("{}", read_gltf(document).unwrap_err());
    assert!(message.contains("external file"), "{message}");

    // A primitive that is not triangles.
    let document = r#"{"asset":{"version":"2.0"},"scene":0,"scenes":[{"nodes":[0]}],
        "nodes":[{"mesh":0}],
        "meshes":[{"primitives":[{"attributes":{"POSITION":0},"mode":5}]}],
        "accessors":[{"componentType":5126,"count":3,"type":"VEC3"}]}"#;
    let message = format!("{}", read_gltf(document).unwrap_err());
    assert!(message.contains("not triangles"), "{message}");

    // Draco.
    let document = r#"{"asset":{"version":"2.0"},"scene":0,"scenes":[{"nodes":[0]}],
        "nodes":[{"mesh":0}],
        "meshes":[{"primitives":[{"attributes":{"POSITION":0},
          "extensions":{"KHR_draco_mesh_compression":{"bufferView":0}}}]}],
        "accessors":[{"componentType":5126,"count":3,"type":"VEC3"}]}"#;
    let message = format!("{}", read_gltf(document).unwrap_err());
    assert!(message.contains("Draco"), "{message}");

    // A node that is its own descendant.
    let document = r#"{"asset":{"version":"2.0"},"scene":0,"scenes":[{"nodes":[0]}],
        "nodes":[{"children":[1]},{"children":[0]}]}"#;
    let message = format!("{}", read_gltf(document).unwrap_err());
    assert!(message.contains("own descendant"), "{message}");

    // An index that names no vertex.
    let document = r#"{"asset":{"version":"2.0"},"scene":0,"scenes":[{"nodes":[0]}],
        "nodes":[{"mesh":0}],
        "meshes":[{"primitives":[{"attributes":{"POSITION":0},"indices":1}]}],
        "accessors":[{"componentType":5126,"count":3,"type":"VEC3"},
                     {"componentType":5125,"count":3,"type":"SCALAR",
                      "sparse":{"count":1,
                        "indices":{"bufferView":0,"componentType":5121},
                        "values":{"bufferView":0}}}],
        "bufferViews":[{"buffer":0,"byteOffset":0,"byteLength":8}],
        "buffers":[{"byteLength":8,"uri":"data:application/octet-stream;base64,AAAACgAAAAA="}]}"#;
    let message = format!("{}", read_gltf(document).unwrap_err());
    assert!(message.contains("names no vertex"), "{message}");
}
