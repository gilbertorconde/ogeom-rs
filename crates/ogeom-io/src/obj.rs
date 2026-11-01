//! Wavefront OBJ, read with its materials.
//!
//! The geometry is the file's `v`, `vn` and `f` records. The colour is in
//! other files: the material libraries `mtllib` names, whose `Kd` is a
//! material's diffuse colour and whose `map_Kd` names a texture image the
//! faces' `vt` coordinates sample. [`read_obj_with_materials`] is handed a
//! way to open those files by the names the OBJ gives them;
//! [`read_obj_file`] opens them beside the OBJ on disk.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;

use crate::image::{Image, Wrap};
use crate::mesh_formats::{MeshColours, held_is_unset};

/// What [`read_obj_with_materials`] found.
#[derive(Debug, Clone)]
pub struct ObjImport {
    /// The mesh, as [`read_obj`] reads it.
    pub mesh: Triangulation,
    /// The colours, where every triangle has one: from its material's
    /// texture at each corner's texture coordinate, or else from its
    /// material's `Kd`; or, where no material colours every triangle, the
    /// per-vertex colours some writers append to `v` records.
    pub colours: Option<MeshColours>,
    /// What was read with a caveat or passed over: a library or texture
    /// that could not be opened or decoded, a material the libraries do
    /// not define, triangles left without a colour.
    pub warnings: Vec<String>,
}

/// One face corner's indices: position, texture coordinate, normal.
type Corner = (usize, Option<usize>, Option<usize>);

/// The records of an OBJ, faces fanned into triangles.
struct Parsed {
    positions: Vec<Point>,
    /// The colour a `v x y z r g b` record appends.
    vertex_colours: Vec<Option<[f64; 4]>>,
    normals: Vec<Vector>,
    /// The `vt` records, `None` for one that does not parse.
    uvs: Vec<Option<(f64, f64)>>,
    triangles: Vec<[u32; 3]>,
    /// Each triangle's corners' texture coordinates.
    corner_uvs: Vec<[Option<usize>; 3]>,
    /// Each triangle's material, an index into `materials`.
    material_of: Vec<Option<usize>>,
    materials: Vec<String>,
    libraries: Vec<String>,
    normal_of: Vec<Option<usize>>,
}

fn parse(text: &str) -> OgeomResult<Parsed> {
    let mut p = Parsed {
        positions: Vec::new(),
        vertex_colours: Vec::new(),
        normals: Vec::new(),
        uvs: Vec::new(),
        triangles: Vec::new(),
        corner_uvs: Vec::new(),
        material_of: Vec::new(),
        materials: Vec::new(),
        libraries: Vec::new(),
        normal_of: Vec::new(),
    };
    let mut material: Option<usize> = None;
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        let mut fields = line.split_whitespace();
        match fields.next() {
            Some("v") => {
                let rest: Vec<&str> = fields.collect();
                let coords = triple(&mut rest.iter().copied(), "a vertex")?;
                p.positions
                    .push(Point::new(coords[0], coords[1], coords[2]));
                p.normal_of.push(None);
                // `v x y z r g b`: a colour in [0, 1] after the position.
                let colour = if rest.len() == 6 {
                    let rgb = rest[3..]
                        .iter()
                        .map(|f| f.parse::<f64>().ok().filter(|c| (0.0..=1.0).contains(c)))
                        .collect::<Option<Vec<f64>>>();
                    rgb.map(|c| [c[0], c[1], c[2], 1.0])
                } else {
                    None
                };
                p.vertex_colours.push(colour);
            }
            Some("vn") => {
                let coords = triple(&mut fields, "a normal")?;
                p.normals.push(Vector::new(coords[0], coords[1], coords[2]));
            }
            Some("vt") => {
                // `u`, then `v` (zero where only `u` is given). A record
                // that does not parse is kept as a gap, so the indices of
                // the rest still line up and the corners naming it read as
                // having no coordinate.
                let mut numbers = fields.map(|f| f.parse::<f64>().ok().filter(|x| x.is_finite()));
                let uv = match (numbers.next().flatten(), numbers.next()) {
                    (Some(u), None) => Some((u, 0.0)),
                    (Some(u), Some(Some(v))) => Some((u, v)),
                    _ => None,
                };
                p.uvs.push(uv);
            }
            Some("usemtl") => {
                let name = fields.collect::<Vec<_>>().join(" ");
                material = Some(match p.materials.iter().position(|m| *m == name) {
                    Some(i) => i,
                    None => {
                        p.materials.push(name);
                        p.materials.len() - 1
                    }
                });
            }
            Some("mtllib") => {
                // Names separated by spaces; a name ending in `.mtl` closes
                // one, so a name holding spaces still reads.
                let mut name = String::new();
                for field in fields {
                    if !name.is_empty() {
                        name.push(' ');
                    }
                    name.push_str(field);
                    if field.to_ascii_lowercase().ends_with(".mtl") {
                        p.libraries.push(std::mem::take(&mut name));
                    }
                }
                if !name.is_empty() {
                    p.libraries.push(name);
                }
            }
            Some("f") => {
                let mut corners: Vec<Corner> = Vec::new();
                for field in fields {
                    corners.push(corner(
                        field,
                        p.positions.len(),
                        p.uvs.len(),
                        p.normals.len(),
                    )?);
                }
                if corners.len() < 3 {
                    ogeom_bail!(
                        Construction,
                        "a face of {} corners is not a face",
                        corners.len()
                    );
                }
                for i in 1..corners.len() - 1 {
                    let fan = [corners[0], corners[i], corners[i + 1]];
                    let mut indices = [0_u32; 3];
                    let mut uvs = [None; 3];
                    for (k, (vertex, uv, normal)) in fan.iter().enumerate() {
                        if let Some(n) = normal {
                            p.normal_of[*vertex] = Some(*n);
                        }
                        indices[k] = u32::try_from(*vertex).unwrap_or(u32::MAX);
                        uvs[k] = *uv;
                    }
                    p.triangles.push(indices);
                    p.corner_uvs.push(uvs);
                    p.material_of.push(material);
                }
            }
            _ => {}
        }
    }
    if p.positions.is_empty() {
        ogeom_bail!(Construction, "the file carries no vertices");
    }
    Ok(p)
}

impl Parsed {
    fn mesh(&self) -> Triangulation {
        let positions = &self.positions;
        // A vertex whose normal the file did not give takes the average of
        // the triangles it belongs to, the same answer the writer would
        // have had.
        let mut resolved = vec![Vector::ZERO; positions.len()];
        for (i, held) in self.normal_of.iter().enumerate() {
            if let Some(n) = held.and_then(|n| self.normals.get(n)) {
                resolved[i] = *n;
            }
        }
        for triangle in &self.triangles {
            let [a, b, c] = triangle.map(|i| positions[i as usize]);
            let face = (b - a).cross(c - a);
            for index in triangle {
                let slot = &mut resolved[*index as usize];
                if held_is_unset(slot) {
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
            .collect::<Vec<_>>();
        Triangulation {
            positions: positions.clone(),
            normals,
            parameters: vec![(0.0, 0.0); positions.len()],
            triangles: self.triangles.clone(),
            deflection_met: true,
        }
    }
}

/// Read an OBJ into a triangulation.
///
/// Vertices, faces and vertex normals; everything else (materials, groups,
/// texture coordinates, smoothing) is a statement about *rendering* a mesh
/// rather than about the mesh, and is skipped rather than half-honoured;
/// [`read_obj_with_materials`] reads the colours. A face of more than
/// three vertices is fanned from its first, which is what a convex polygon
/// means and what OBJ writers emit.
///
/// Indices may be negative, which in OBJ counts back from the end, and may
/// carry the `v/vt/vn` triple, of which the first field is the position.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if an
/// index names a vertex that is not there, or a coordinate does not parse.
pub fn read_obj(text: &str) -> OgeomResult<Triangulation> {
    Ok(parse(text)?.mesh())
}

/// Three floats off an iterator.
fn triple<'a>(fields: &mut impl Iterator<Item = &'a str>, what: &str) -> OgeomResult<[f64; 3]> {
    let mut out = [0.0; 3];
    for slot in &mut out {
        let Some(field) = fields.next() else {
            ogeom_bail!(Construction, "{what} needs three coordinates");
        };
        let Ok(value) = field.parse::<f64>() else {
            ogeom_bail!(
                Construction,
                "{what} carries {field}, which is not a number"
            );
        };
        *slot = value;
    }
    Ok(out)
}

/// One OBJ face corner: `v`, `v/vt`, `v//vn` or `v/vt/vn`, one-based, and
/// negative counting back from the end.
fn corner(field: &str, vertices: usize, uvs: usize, normals: usize) -> OgeomResult<Corner> {
    let mut parts = field.split('/');
    let resolve = |text: &str, count: usize| -> Option<usize> {
        let index = text.parse::<isize>().ok()?;
        let zero = if index > 0 {
            usize::try_from(index).ok()?.checked_sub(1)?
        } else if index < 0 {
            count.checked_sub(usize::try_from(-index).ok()?)?
        } else {
            return None;
        };
        (zero < count).then_some(zero)
    };
    let Some(vertex) = parts.next().and_then(|t| resolve(t, vertices)) else {
        ogeom_bail!(
            Construction,
            "a face names vertex {field}, which is not there"
        );
    };
    let uv = parts.next().and_then(|t| resolve(t, uvs));
    let normal = parts.next().and_then(|t| resolve(t, normals));
    Ok((vertex, uv, normal))
}

/// A material as a library defines it.
#[derive(Default)]
struct Material {
    diffuse: Option<[f64; 3]>,
    opacity: f64,
    texture: Option<Texture>,
}

/// A `map_Kd`: the image, and the options that place it.
struct Texture {
    file: String,
    offset: (f64, f64),
    scale: (f64, f64),
    clamp: bool,
    image: Option<Image>,
}

/// The file a `map_Kd` names after its options, and the options it honours:
/// `-o` and `-s` (offset and scale of the coordinates) and `-clamp`.
fn texture_statement(fields: &[&str]) -> Texture {
    let mut texture = Texture {
        file: String::new(),
        offset: (0.0, 0.0),
        scale: (1.0, 1.0),
        clamp: false,
        image: None,
    };
    let mut i = 0;
    while i < fields.len() && fields[i].starts_with('-') {
        let option = fields[i];
        i += 1;
        // How many arguments the option takes: the three-number options
        // take up to three, as many as are numbers.
        let numbers = |i: &mut usize| -> Vec<f64> {
            let mut out = Vec::new();
            while out.len() < 3
                && let Some(v) = fields.get(*i).and_then(|f| f.parse::<f64>().ok())
            {
                out.push(v);
                *i += 1;
            }
            out
        };
        match option {
            "-o" => {
                let v = numbers(&mut i);
                texture.offset = (
                    v.first().copied().unwrap_or(0.0),
                    v.get(1).copied().unwrap_or(0.0),
                );
            }
            "-s" => {
                let v = numbers(&mut i);
                texture.scale = (
                    v.first().copied().unwrap_or(1.0),
                    v.get(1).copied().unwrap_or(1.0),
                );
            }
            "-t" => {
                numbers(&mut i);
            }
            "-clamp" => {
                texture.clamp = fields.get(i) == Some(&"on");
                i += 1;
            }
            "-mm" => i += 2,
            _ => i += 1,
        }
    }
    texture.file = fields.get(i..).unwrap_or_default().join(" ");
    texture
}

/// The materials a library defines, by name.
fn library(text: &str, into: &mut Vec<(String, Material)>, warnings: &mut Vec<String>) {
    let mut current: Option<usize> = None;
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some((&keyword, rest)) = fields.split_first() else {
            continue;
        };
        if keyword == "newmtl" {
            into.push((
                rest.join(" "),
                Material {
                    opacity: 1.0,
                    ..Material::default()
                },
            ));
            current = Some(into.len() - 1);
            continue;
        }
        let Some(material) = current.map(|i| &mut into[i]) else {
            continue;
        };
        let name = &material.0;
        match keyword {
            "Kd" => {
                let rgb: Option<Vec<f64>> = rest.iter().map(|f| f.parse::<f64>().ok()).collect();
                match rgb.as_deref() {
                    Some([r, g, b]) => material.1.diffuse = Some([*r, *g, *b]),
                    Some([v]) => material.1.diffuse = Some([*v; 3]),
                    _ => warnings.push(format!(
                        "the material {name} gives Kd as {}, which is not read",
                        rest.join(" ")
                    )),
                }
            }
            "d" => {
                if let Some(d) = rest.last().and_then(|f| f.parse::<f64>().ok()) {
                    material.1.opacity = d.clamp(0.0, 1.0);
                }
            }
            "Tr" => {
                if let Some(t) = rest.first().and_then(|f| f.parse::<f64>().ok()) {
                    material.1.opacity = (1.0 - t).clamp(0.0, 1.0);
                }
            }
            "map_Kd" => material.1.texture = Some(texture_statement(rest)),
            _ => {}
        }
    }
}

/// Read an OBJ with the colours its materials give it.
///
/// The mesh is the one [`read_obj`] reads. Each library `mtllib` names is
/// asked of `open` by the name the OBJ gives it, and each texture a
/// library's `map_Kd` names by the name the library gives it; `open`
/// answers with the file's bytes, or `None` where there is no such file.
///
/// A triangle takes the colour of the material `usemtl` last named before
/// it: with a texture and texture coordinates at its corners, the texel at
/// each corner's coordinate (after the `-o` offset and `-s` scale the
/// statement gives; repeating outside `[0, 1]`, or clamped under
/// `-clamp on`); otherwise the material's `Kd`. A texture stands for the
/// diffuse colour rather than scaling `Kd`, as scan textures are written
/// to be shown. Alpha is the material's `d` (or one less its `Tr`), times
/// the texel's. Textures are PNG or JPEG.
///
/// What cannot be opened or decoded (a library, a texture, a material the
/// libraries do not define) is a warning, and the triangles it would have
/// coloured fall back as above or go without. Colours are returned only
/// where every triangle has one; see [`ObjImport::colours`].
///
/// # Errors
///
/// As [`read_obj`].
pub fn read_obj_with_materials(
    text: &str,
    open: &mut dyn FnMut(&str) -> Option<Vec<u8>>,
) -> OgeomResult<ObjImport> {
    let parsed = parse(text)?;
    let mesh = parsed.mesh();
    let mut warnings = Vec::new();
    let mut defined: Vec<(String, Material)> = Vec::new();
    if !parsed.materials.is_empty() {
        for name in &parsed.libraries {
            match open(name) {
                Some(bytes) => library(
                    &String::from_utf8_lossy(&bytes),
                    &mut defined,
                    &mut warnings,
                ),
                None => warnings.push(format!("the material library {name} could not be opened")),
            }
        }
    }
    let mut budget = crate::image::Budget::new();
    for (name, material) in &mut defined {
        let Some(texture) = material.texture.as_mut() else {
            continue;
        };
        if !parsed
            .material_of
            .contains(&parsed.materials.iter().position(|m| m == name))
        {
            continue;
        }
        match open(&texture.file) {
            None => warnings.push(format!(
                "the texture {} of {name} could not be opened",
                texture.file
            )),
            Some(bytes) => match budget.decode(&bytes) {
                Ok(image) => texture.image = Some(image),
                Err(e) => warnings.push(format!(
                    "the texture {} of {name} could not be decoded: {e}",
                    texture.file
                )),
            },
        }
    }
    let used: Vec<Option<&Material>> = parsed
        .materials
        .iter()
        .map(|name| {
            let found = defined
                .iter()
                .rev()
                .find(|(n, _)| n == name)
                .map(|(_, m)| m);
            if found.is_none() {
                warnings.push(format!("the material {name} is not defined"));
            }
            found
        })
        .collect();

    let mut corners: Vec<[[f64; 4]; 3]> = Vec::with_capacity(parsed.triangles.len());
    let mut uncoloured = 0_usize;
    let mut untextured = 0_usize;
    for (t, material) in parsed.material_of.iter().enumerate() {
        let Some(material) = material.and_then(|m| used[m]) else {
            uncoloured += 1;
            continue;
        };
        let alpha = material.opacity;
        let sampled = material
            .texture
            .as_ref()
            .and_then(|tex| Some((tex, tex.image.as_ref()?)))
            .and_then(|(tex, image)| {
                let uvs = parsed.corner_uvs[t];
                let wrap = if tex.clamp { Wrap::Clamp } else { Wrap::Repeat };
                let mut out = [[0.0; 4]; 3];
                for (slot, uv) in out.iter_mut().zip(uvs) {
                    let (u, v) = parsed.uvs[uv?]?;
                    let (u, v) = (
                        u * tex.scale.0 + tex.offset.0,
                        v * tex.scale.1 + tex.offset.1,
                    );
                    let [r, g, b, a] = image.sample(u, v, [wrap; 2]);
                    *slot = [r, g, b, a * alpha];
                }
                Some(out)
            });
        if let Some(colours) = sampled {
            corners.push(colours);
            continue;
        }
        if material.texture.as_ref().is_some_and(|t| t.image.is_some()) {
            untextured += 1;
        }
        match material.diffuse {
            Some([r, g, b]) => corners.push([[r, g, b, alpha]; 3]),
            None => uncoloured += 1,
        }
    }
    if untextured > 0 {
        warnings.push(format!(
            "{untextured} textured triangles have no texture coordinates and take their material's Kd"
        ));
    }
    let colours = if uncoloured == 0 && !corners.is_empty() {
        Some(MeshColours::from_corners(
            &parsed.triangles,
            parsed.positions.len(),
            corners,
        ))
    } else if let Some(vertex) = parsed
        .vertex_colours
        .iter()
        .copied()
        .collect::<Option<Vec<[f64; 4]>>>()
    {
        Some(MeshColours::PerVertex(vertex))
    } else {
        if !corners.is_empty() {
            warnings.push(format!(
                "{uncoloured} triangles have no colour, so the colours of the rest are not returned"
            ));
        }
        None
    };
    Ok(ObjImport {
        mesh,
        colours,
        warnings,
    })
}

/// Read an OBJ file from disk with its materials, as
/// [`read_obj_with_materials`] reads them.
///
/// A library or texture is looked for at the relative path its statement
/// gives, under the OBJ's own directory; a path that is absolute, or that
/// climbs out with `..`, or that is not found there, is looked for by its
/// last component in that directory, and nowhere else. Backslashes
/// separate components as slashes do.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `path` cannot be read, and as [`read_obj_with_materials`].
pub fn read_obj_file(path: impl AsRef<std::path::Path>) -> OgeomResult<ObjImport> {
    let path = path.as_ref();
    let bytes = match std::fs::read(path) {
        Ok(b) => b,
        Err(e) => ogeom_bail!(Construction, "{} could not be read: {e}", path.display()),
    };
    let directory = path
        .parent()
        .map(std::path::Path::to_path_buf)
        .unwrap_or_default();
    let mut open = |name: &str| -> Option<Vec<u8>> {
        let name = name.replace('\\', "/");
        let parts: Vec<&str> = name
            .split('/')
            .filter(|p| !p.is_empty() && *p != ".")
            .collect();
        let contained =
            !name.starts_with('/') && !parts.iter().any(|p| *p == ".." || p.contains(':'));
        if contained {
            let mut nested = directory.clone();
            nested.extend(&parts);
            if let Ok(bytes) = std::fs::read(&nested) {
                return Some(bytes);
            }
        }
        std::fs::read(directory.join(parts.last()?)).ok()
    };
    read_obj_with_materials(&String::from_utf8_lossy(&bytes), &mut open)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_core::FastMap;

    /// A unit cube, each face its own quad with its own texture
    /// coordinates, laid out as six cells of a 3 × 2 texture.
    fn textured_cube(material: &str) -> String {
        let mut text = format!("mtllib cube.mtl\nusemtl {material}\n");
        for z in [0, 1] {
            for y in [0, 1] {
                for x in [0, 1] {
                    text.push_str(&format!("v {x} {y} {z}\n"));
                }
            }
        }
        // Cell (i, j) of the 3 × 2 grid, sampled at its centre.
        for i in 0..6 {
            let (u, v) = (
                (f64::from(i % 3) + 0.5) / 3.0,
                (f64::from(i / 3) + 0.5) / 2.0,
            );
            text.push_str(&format!("vt {u} {v}\n"));
        }
        // Six faces, outward, one-based; each corner the face's cell.
        let faces = [
            [1, 3, 4, 2],
            [5, 6, 8, 7],
            [1, 2, 6, 5],
            [3, 7, 8, 4],
            [1, 5, 7, 3],
            [2, 4, 8, 6],
        ];
        for (k, f) in faces.iter().enumerate() {
            let t = k + 1;
            text.push_str(&format!(
                "f {}/{t} {}/{t} {}/{t} {}/{t}\n",
                f[0], f[1], f[2], f[3]
            ));
        }
        text
    }

    const CELLS: [[u8; 4]; 6] = [
        [255, 0, 0, 255],
        [0, 255, 0, 255],
        [0, 0, 255, 255],
        [255, 255, 0, 255],
        [0, 255, 255, 255],
        [255, 0, 255, 255],
    ];

    /// The 3 × 2 texture: row 0 is the top, so cells 3..6 (v in the upper
    /// half) are the first row.
    fn texture() -> Vec<u8> {
        let mut texels = Vec::new();
        texels.extend_from_slice(&CELLS[3..6]);
        texels.extend_from_slice(&CELLS[0..3]);
        crate::image::encode_png(3, 2, &texels)
    }

    fn files() -> FastMap<String, Vec<u8>> {
        let mut files = FastMap::default();
        files.insert(
            "cube.mtl".to_owned(),
            b"newmtl painted\nKd 0.8 0.8 0.8\nmap_Kd -s 1 1 1 textures/cube.png\n\
newmtl plain\nKd 0.25 0.5 1\nd 0.5\n"
                .to_vec(),
        );
        files.insert("textures/cube.png".to_owned(), texture());
        files
    }

    #[test]
    fn a_textured_cube_reads_its_face_colours() {
        let files = files();
        let import =
            read_obj_with_materials(&textured_cube("painted"), &mut |n| files.get(n).cloned())
                .unwrap();
        assert!(import.warnings.is_empty(), "{:?}", import.warnings);
        assert_eq!(import.mesh.triangles.len(), 12);
        let Some(MeshColours::PerTriangle(colours)) = import.colours else {
            panic!("a colour per triangle: {:?}", import.colours);
        };
        for (t, colour) in colours.iter().enumerate() {
            assert_eq!(
                *colour,
                CELLS[t / 2].map(|c| f64::from(c) / 255.0),
                "triangle {t}"
            );
        }
    }

    #[test]
    fn kd_colours_faces_without_a_texture() {
        let files = files();
        let import =
            read_obj_with_materials(&textured_cube("plain"), &mut |n| files.get(n).cloned())
                .unwrap();
        assert_eq!(
            import.colours,
            Some(MeshColours::PerVertex(vec![[0.25, 0.5, 1.0, 0.5]; 8]))
        );
        // The texture missing, the painted material falls back to its Kd.
        let mut missing = files.clone();
        missing.remove("textures/cube.png");
        let import =
            read_obj_with_materials(&textured_cube("painted"), &mut |n| missing.get(n).cloned())
                .unwrap();
        assert_eq!(
            import.colours,
            Some(MeshColours::PerVertex(vec![[0.8, 0.8, 0.8, 1.0]; 8]))
        );
        assert!(
            import
                .warnings
                .iter()
                .any(|w| w.contains("could not be opened"))
        );
        // A library that is not there leaves the mesh uncoloured, and says so.
        let import = read_obj_with_materials(&textured_cube("painted"), &mut |_| None).unwrap();
        assert!(import.colours.is_none());
        assert_eq!(import.mesh, read_obj(&textured_cube("painted")).unwrap());
        assert!(!import.warnings.is_empty());
    }

    #[test]
    fn vertex_colours_read_where_no_material_colours() {
        let text = "v 0 0 0 1 0 0\nv 1 0 0 0 1 0\nv 0 1 0 0 0 1\nf 1 2 3\n";
        let import = read_obj_with_materials(text, &mut |_| None).unwrap();
        assert_eq!(
            import.colours,
            Some(MeshColours::PerVertex(vec![
                [1.0, 0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0, 1.0],
                [0.0, 0.0, 1.0, 1.0]
            ]))
        );
    }

    #[test]
    fn hostile_materials_and_textures_never_panic() {
        let text = textured_cube("painted");
        let png = texture();
        for cut in 0..png.len() {
            let mut files = files();
            files.insert("textures/cube.png".to_owned(), png[..cut].to_vec());
            let import = read_obj_with_materials(&text, &mut |n| files.get(n).cloned()).unwrap();
            assert!(!import.warnings.is_empty(), "cut at {cut}");
        }
        let mut files = files();
        files.insert(
            "cube.mtl".to_owned(),
            b"Kd 1 1 1\nnewmtl painted\nKd spectral x.rfl\nmap_Kd -o\nmap_Kd -clamp\nd\nTr x\n"
                .to_vec(),
        );
        let _ = read_obj_with_materials(&text, &mut |n| files.get(n).cloned()).unwrap();
        let gap = text.replace("vt 0.5 0.75", "vt nan 0.75");
        let fresh = super::tests::files();
        let import = read_obj_with_materials(&gap, &mut |n| fresh.get(n).cloned()).unwrap();
        assert!(
            import
                .warnings
                .iter()
                .any(|w| w.contains("no texture coordinates"))
        );
        assert!(read_obj_with_materials("v 0 0 0\nf 1 2 3\n", &mut |_| None).is_err());
    }
}
