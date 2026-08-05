//! STL: a triangle soup, in ASCII or binary.
//!
//! The format has no notion of a shared vertex — every triangle names its three
//! corners in full — and no notion of topology at all. Writing to it therefore
//! *loses* everything the kernel knows: which triangles belong to which face,
//! which edges were exact, what the surfaces actually were. That is the format's
//! nature and not a shortcoming of this writer.
//!
//! # Reading gives a mesh, not a shape
//!
//! [`read`] returns a [`Triangulation`], not a `Shape`, and that is deliberate.
//! Recovering a B-rep from a triangle soup means deciding which triangles are
//! coplanar enough to be one face, which chains of edges are one curve, and
//! what surface each face was cut from — that is surface reconstruction, a
//! research problem, not a file format concern. A function returning a `Shape`
//! here would have to guess, and the guess would be wrong in ways nothing
//! downstream could detect.
//!
//! # The normals are written and ignored
//!
//! Every STL triangle carries a facet normal. Most writers get it right, some
//! write zeros, and some write it inconsistent with the winding. [`read`]
//! therefore recomputes normals from the winding and discards what the file
//! said, which is what every robust reader does. [`write()`] emits the true
//! normal, because a reader that trusts it should not be punished for it.

use std::fmt::Write as _;

use og_core::{OgResult, Tolerances, og_bail};
use og_math::{Point, Vector};
use og_topo::Triangulation;

/// Which encoding to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// Human-readable. Roughly six times the size and exact only to the
    /// precision printed.
    Ascii,
    /// Compact and exact, but `f32` — see [`write()`].
    Binary,
}

/// The header a binary STL carries, and the name an ASCII one does.
const HEADER: &str = "ogeom";

/// Write a triangulation as STL.
///
/// # Precision
///
/// STL stores coordinates as `f32` in binary and as printed decimals in ASCII.
/// The kernel works in `f64`. A round trip therefore loses precision — about
/// seven significant digits in binary — and a model far from the origin loses
/// it where it matters most: a part at 1e6 units has a binary STL resolution of
/// about 0.06 units. This writes what the format can hold; it does not pretend
/// the result round-trips exactly, and [`read`] will not return what was
/// written.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the mesh has no
/// triangles, or holds a triangle index that is not a vertex.
pub fn write(mesh: &Triangulation, encoding: Encoding) -> OgResult<Vec<u8>> {
    if mesh.triangles.is_empty() {
        og_bail!(Construction, "an STL with no triangles describes nothing");
    }
    for triangle in &mesh.triangles {
        for index in triangle {
            if *index as usize >= mesh.positions.len() {
                og_bail!(
                    Construction,
                    "a triangle names vertex {index}, and the mesh has {}",
                    mesh.positions.len()
                );
            }
        }
    }
    Ok(match encoding {
        Encoding::Ascii => write_ascii(mesh).into_bytes(),
        Encoding::Binary => write_binary(mesh),
    })
}

/// Read STL, in either encoding.
///
/// The encoding is detected rather than asked for: a file that starts with
/// `solid` is *usually* ASCII, but binary files written by several well-known
/// programs start with it too, because their 80-byte header happens to. The
/// reliable test is whether the file's length matches what its triangle count
/// claims, so that is the test used.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the bytes are
/// not STL of either kind, or are truncated part-way through a triangle.
pub fn read(bytes: &[u8], tol: Tolerances) -> OgResult<Triangulation> {
    if looks_binary(bytes) {
        read_binary(bytes, tol)
    } else {
        read_ascii(bytes, tol)
    }
}

/// Whether the bytes are a binary STL.
///
/// Decided by arithmetic, not by the leading word. A binary file is 84 bytes of
/// header and count plus exactly 50 per triangle; nothing else lands on that
/// length by accident, and plenty of binary files begin with `solid`.
fn looks_binary(bytes: &[u8]) -> bool {
    const HEADER_BYTES: usize = 84;
    const PER_TRIANGLE: usize = 50;
    if bytes.len() < HEADER_BYTES {
        return false;
    }
    let Ok(count) = <[u8; 4]>::try_from(&bytes[80..84]) else {
        return false;
    };
    let claimed = u32::from_le_bytes(count) as usize;
    bytes.len() == HEADER_BYTES + claimed * PER_TRIANGLE
}

/// Render the ASCII form.
fn write_ascii(mesh: &Triangulation) -> String {
    let mut out = String::with_capacity(mesh.triangles.len() * 260);
    let _ = writeln!(out, "solid {HEADER}");
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
        let n = facet_normal(a, b, c);
        let _ = writeln!(out, "  facet normal {:e} {:e} {:e}", n.x, n.y, n.z);
        let _ = writeln!(out, "    outer loop");
        for p in [a, b, c] {
            let _ = writeln!(out, "      vertex {:e} {:e} {:e}", p.x, p.y, p.z);
        }
        let _ = writeln!(out, "    endloop");
        let _ = writeln!(out, "  endfacet");
    }
    let _ = writeln!(out, "endsolid {HEADER}");
    out
}

/// Render the binary form.
fn write_binary(mesh: &Triangulation) -> Vec<u8> {
    let mut out = Vec::with_capacity(84 + mesh.triangles.len() * 50);
    let mut header = [b' '; 80];
    let name = HEADER.as_bytes();
    header[..name.len()].copy_from_slice(name);
    out.extend_from_slice(&header);

    #[allow(clippy::cast_possible_truncation)]
    out.extend_from_slice(&(mesh.triangles.len() as u32).to_le_bytes());

    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
        let n = facet_normal(a, b, c);
        for v in [n, a.to_vector(), b.to_vector(), c.to_vector()] {
            #[allow(clippy::cast_possible_truncation)]
            for component in [v.x as f32, v.y as f32, v.z as f32] {
                out.extend_from_slice(&component.to_le_bytes());
            }
        }
        // The attribute-count field. Some tools smuggle colour through it;
        // writing anything but zero makes the file unreadable to the ones that
        // do not expect it.
        out.extend_from_slice(&0_u16.to_le_bytes());
    }
    out
}

/// Parse the binary form.
fn read_binary(bytes: &[u8], tol: Tolerances) -> OgResult<Triangulation> {
    let Ok(count) = <[u8; 4]>::try_from(&bytes[80..84]) else {
        og_bail!(Construction, "the binary STL header is truncated");
    };
    let count = u32::from_le_bytes(count) as usize;

    let mut mesh = Triangulation::new();
    for i in 0..count {
        let at = 84 + i * 50;
        // The facet normal is read past rather than used; see the module docs.
        let corners: Vec<Point> = (0..3)
            .map(|k| {
                let base = at + 12 + k * 12;
                Point::new(
                    f64::from(read_f32(bytes, base)),
                    f64::from(read_f32(bytes, base + 4)),
                    f64::from(read_f32(bytes, base + 8)),
                )
            })
            .collect();
        push(&mut mesh, corners[0], corners[1], corners[2]);
    }
    Ok(mesh.welded(tol))
}

/// One little-endian `f32`, or zero past the end.
///
/// A truncated file is caught by [`looks_binary`] before this runs, so the
/// fallback is unreachable in practice; returning zero rather than panicking
/// keeps a malformed file a bad *mesh* instead of a crash.
fn read_f32(bytes: &[u8], at: usize) -> f32 {
    <[u8; 4]>::try_from(bytes.get(at..at + 4).unwrap_or(&[0; 4]))
        .map(f32::from_le_bytes)
        .unwrap_or(0.0)
}

/// Parse the ASCII form.
///
/// Deliberately lenient about layout — indentation, blank lines, and the solid
/// name vary between writers, and the keywords are what carry the meaning.
/// Deliberately strict about a facet having exactly three vertices, because a
/// facet with four is a quad some writer emitted and silently dropping one
/// corner would put a hole in the mesh.
fn read_ascii(bytes: &[u8], tol: Tolerances) -> OgResult<Triangulation> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        og_bail!(
            Construction,
            "these bytes are neither binary STL nor valid UTF-8, so they are \
             not ASCII STL either"
        );
    };

    let mut mesh = Triangulation::new();
    let mut corners: Vec<Point> = Vec::with_capacity(3);
    let mut in_facet = false;
    let mut facets = 0_usize;

    for (line_number, line) in text.lines().enumerate() {
        let mut words = line.split_whitespace();
        let Some(keyword) = words.next() else {
            continue;
        };
        match keyword {
            "facet" => {
                in_facet = true;
                corners.clear();
            }
            "vertex" => {
                let values: Vec<f64> = words.filter_map(|w| w.parse().ok()).collect();
                if values.len() != 3 {
                    og_bail!(
                        Construction,
                        "line {}: a vertex needs three numbers, got {:?}",
                        line_number + 1,
                        line.trim()
                    );
                }
                corners.push(Point::new(values[0], values[1], values[2]));
            }
            "endfacet" => {
                if corners.len() != 3 {
                    og_bail!(
                        Construction,
                        "line {}: a facet has {} vertices; STL facets are \
                         triangles, and dropping a corner would put a hole in \
                         the mesh",
                        line_number + 1,
                        corners.len()
                    );
                }
                push(&mut mesh, corners[0], corners[1], corners[2]);
                in_facet = false;
                facets += 1;
            }
            _ => {}
        }
    }

    if in_facet {
        og_bail!(Construction, "the file ends part-way through a facet");
    }
    if facets == 0 {
        og_bail!(
            Construction,
            "no facets found; these bytes are not STL of either kind"
        );
    }
    Ok(mesh.welded(tol))
}

/// Append one triangle, with its vertices unshared.
///
/// Welding happens afterwards, over the whole mesh at once. Sharing as we go
/// would need a lookup per vertex against everything read so far, and get the
/// same answer more slowly.
fn push(mesh: &mut Triangulation, a: Point, b: Point, c: Point) {
    #[allow(clippy::cast_possible_truncation)]
    let base = mesh.positions.len() as u32;
    let normal = facet_normal(a, b, c);
    for p in [a, b, c] {
        mesh.positions.push(p);
        mesh.normals.push(normal);
        mesh.parameters.push((0.0, 0.0));
    }
    mesh.triangles.push([base, base + 1, base + 2]);
}

/// The unit normal a triangle's winding implies, or zero if it has no area.
fn facet_normal(a: Point, b: Point, c: Point) -> Vector {
    let n = (b - a).cross(c - a);
    let length = n.magnitude();
    if length <= f64::MIN_POSITIVE {
        Vector::ZERO
    } else {
        n / length
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    /// A unit tetrahedron, wound outward.
    fn tetrahedron() -> Triangulation {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        mesh.normals = vec![Vector::Z; 4];
        mesh.parameters = vec![(0.0, 0.0); 4];
        mesh.triangles = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        mesh
    }

    #[test]
    fn a_mesh_survives_a_round_trip_through_either_encoding() {
        let original = tetrahedron();
        for encoding in [Encoding::Ascii, Encoding::Binary] {
            let bytes = write(&original, encoding).unwrap();
            let back = read(&bytes, T).unwrap();

            assert_eq!(
                back.triangle_count(),
                original.triangle_count(),
                "{encoding:?}"
            );
            assert_eq!(
                back.vertex_count(),
                original.vertex_count(),
                "{encoding:?}: welding should recover the shared corners"
            );
            assert!(back.is_closed(), "{encoding:?}: the tetrahedron came apart");
            assert_relative_eq!(back.volume(), original.volume(), epsilon = 1e-6);
        }
    }

    #[test]
    fn the_encoding_is_detected_by_length_not_by_the_leading_word() {
        // Several well-known programs write binary files whose 80-byte header
        // begins with "solid". Sniffing the first word sends those down the
        // ASCII path and produces an empty mesh with no error.
        let mut bytes = write(&tetrahedron(), Encoding::Binary).unwrap();
        bytes[..5].copy_from_slice(b"solid");

        assert!(looks_binary(&bytes), "length is what decides");
        let back = read(&bytes, T).unwrap();
        assert_eq!(back.triangle_count(), 4);
    }

    #[test]
    fn the_written_normal_agrees_with_the_winding() {
        // A reader that trusts the facet normal should not be punished for it,
        // even though this one does not.
        let text = String::from_utf8(write(&tetrahedron(), Encoding::Ascii).unwrap()).unwrap();
        let normals: Vec<Vec<f64>> = text
            .lines()
            .filter(|l| l.trim_start().starts_with("facet normal"))
            .map(|l| {
                l.split_whitespace()
                    .filter_map(|w| w.parse().ok())
                    .collect()
            })
            .collect();
        assert_eq!(normals.len(), 4);

        let mesh = tetrahedron();
        for (written, triangle) in normals.iter().zip(&mesh.triangles) {
            let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
            let expected = facet_normal(a, b, c);
            assert_relative_eq!(written[0], expected.x, epsilon = 1e-9);
            assert_relative_eq!(written[1], expected.y, epsilon = 1e-9);
            assert_relative_eq!(written[2], expected.z, epsilon = 1e-9);
        }
    }

    #[test]
    fn a_read_mesh_takes_its_normals_from_the_winding_not_the_file() {
        // Files with zeroed or inverted facet normals are common. Trusting them
        // makes a mesh that renders inside out and whose volume comes out
        // negative, with nothing in the geometry to say why.
        let lying = "solid liar
  facet normal 0 0 -1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 0 1 0
    endloop
  endfacet
endsolid liar
";
        let mesh = read(lying.as_bytes(), T).unwrap();
        assert_eq!(mesh.triangle_count(), 1);
        for normal in &mesh.normals {
            assert_relative_eq!(normal.z, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn a_facet_that_is_not_a_triangle_is_refused_rather_than_trimmed() {
        // Some writers emit quads. Keeping the first three corners looks like it
        // works and leaves a hole where the fourth was.
        let quad = "solid q
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 1 1 0
      vertex 0 1 0
    endloop
  endfacet
endsolid q
";
        let refused = read(quad.as_bytes(), T);
        assert!(refused.is_err());
        assert!(format!("{}", refused.unwrap_err()).contains("triangles"));
    }

    #[test]
    fn a_truncated_or_empty_file_is_refused() {
        assert!(read(b"", T).is_err());
        assert!(read(b"solid x\nendsolid x\n", T).is_err(), "no facets");
        assert!(
            read(
                b"solid x\n  facet normal 0 0 1\n    outer loop\n      vertex 0 0 0\n",
                T
            )
            .is_err(),
            "ends inside a facet"
        );
    }

    #[test]
    fn a_mesh_with_nothing_in_it_is_not_written() {
        assert!(write(&Triangulation::new(), Encoding::Ascii).is_err());

        let mut broken = tetrahedron();
        broken.triangles.push([0, 1, 99]);
        assert!(write(&broken, Encoding::Binary).is_err());
    }

    #[test]
    fn binary_is_the_length_the_format_says_it_is() {
        let mesh = tetrahedron();
        let bytes = write(&mesh, Encoding::Binary).unwrap();
        assert_eq!(bytes.len(), 84 + 4 * 50);
        assert_eq!(&bytes[..5], b"ogeom");
        assert_eq!(u32::from_le_bytes(bytes[80..84].try_into().unwrap()), 4);
    }

    #[test]
    fn a_binary_round_trip_loses_precision_and_the_docs_say_so() {
        // Not a defect to fix — the format stores f32. The test exists so that
        // if someone later "fixes" a failing comparison by loosening a
        // tolerance, they meet this instead and learn where the loss comes from.
        let mut mesh = tetrahedron();
        mesh.positions[1] = Point::new(1.000_000_1, 0.0, 0.0);
        let back = read(&write(&mesh, Encoding::Binary).unwrap(), T).unwrap();

        let moved = back
            .positions
            .iter()
            .any(|p| (p.x - 1.000_000_1).abs() > 1e-9 && (p.x - 1.0).abs() < 1e-3);
        assert!(moved, "f32 should have rounded the seventh digit away");
    }
}
