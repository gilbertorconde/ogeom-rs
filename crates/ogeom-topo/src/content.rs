//! A face's content written out, for keeping what was computed from a face
//! across models in which it recurs.

use ogeom_geom::{Curve, Curve3d as _, PlanarCurve};
use ogeom_math::{KnotVector, Transform};

use crate::{EdgeRepr, Model, Shape, TShapeId};

/// Everything a face is made of, written out as bytes: its orientation,
/// placement, tolerance and surface, then each wire in order with each of
/// its edges' orientation, placement, tolerance, flags, every
/// representation (a curve, a pcurve, a seam's two pcurves, with their
/// placements and ranges, and whether a pcurve lies on the face's own
/// surface or another face's) and vertices. An edge met again in the face
/// is named by where it was first met, so a seam walked twice reads as one
/// edge. Geometry is written out, never named by its handle, so the same
/// face built again in another model reads the same.
///
/// Floats are written as their bits, or in their debug form, which reads
/// back to the same bits, and every piece of varying length carries its
/// length: two faces read the same only where every one of these is the
/// same. `None` where any part of the face does not resolve, or an edge
/// carries a mesh's indices, which name a triangulation of this model by
/// its handle.
#[must_use]
pub fn face_content(model: &Model, face: &Shape) -> Option<Vec<u8>> {
    let mut out = Vec::new();
    let data = model.node(face)?.data().as_face()?;
    let geometry = model.geometry();
    let datums = model.datums();
    text(&mut out, format_args!("F{:?}", face.orientation()));
    placement(&mut out, face.transform(datums).ok()?);
    text(
        &mut out,
        format_args!(
            "{:?}{}{:?}",
            data.tolerance,
            data.natural_restriction,
            geometry.surface(data.surface)?
        ),
    );
    let mut met: Vec<TShapeId> = Vec::new();
    let wires = model.ordered_children_of(face).ok()?;
    count(&mut out, wires.len());
    for wire in wires {
        text(&mut out, format_args!("W{:?}", wire.orientation()));
        placement(&mut out, wire.transform(datums).ok()?);
        let edges = model.ordered_children_of(&wire).ok()?;
        count(&mut out, edges.len());
        for edge in edges {
            let at = met
                .iter()
                .position(|&e| e == edge.node())
                .unwrap_or_else(|| {
                    met.push(edge.node());
                    met.len() - 1
                });
            let e = model.node(&edge)?.data().as_edge()?;
            text(&mut out, format_args!("E{at}{:?}", edge.orientation()));
            placement(&mut out, edge.transform(datums).ok()?);
            text(
                &mut out,
                format_args!(
                    "{:?}{}{}{}",
                    e.tolerance,
                    e.degenerate,
                    e.same_parameter(),
                    e.representations.len()
                ),
            );
            for repr in &e.representations {
                match repr {
                    EdgeRepr::Curve3d {
                        curve,
                        location,
                        range,
                    } => {
                        out.push(b'C');
                        space_curve(&mut out, geometry.curve(*curve)?);
                        placement(&mut out, location.composed(datums).ok()?);
                        bits(&mut out, &[range.0, range.1]);
                    }
                    EdgeRepr::PCurve {
                        curve,
                        surface,
                        location,
                        range,
                    } => {
                        out.push(b'P');
                        planar_curve(&mut out, geometry.pcurve(*curve)?);
                        out.push(u8::from(*surface == data.surface));
                        placement(&mut out, location.composed(datums).ok()?);
                        bits(&mut out, &[range.0, range.1]);
                    }
                    EdgeRepr::Seam {
                        forward,
                        reversed,
                        surface,
                        location,
                        range,
                    } => {
                        out.push(b'S');
                        planar_curve(&mut out, geometry.pcurve(*forward)?);
                        planar_curve(&mut out, geometry.pcurve(*reversed)?);
                        out.push(u8::from(*surface == data.surface));
                        placement(&mut out, location.composed(datums).ok()?);
                        bits(&mut out, &[range.0, range.1]);
                    }
                    EdgeRepr::PolygonOnTriangulation { .. } => return None,
                    other => {
                        out.push(b'R');
                        text(&mut out, format_args!("{other:?}"));
                    }
                }
            }
            let vertices = model.ordered_children_of(&edge).ok()?;
            count(&mut out, vertices.len());
            for vertex in vertices {
                let v = model.node(&vertex)?.data().as_vertex()?;
                text(&mut out, format_args!("V{:?}", vertex.orientation()));
                placement(&mut out, vertex.transform(datums).ok()?);
                bits(&mut out, &[v.point.x, v.point.y, v.point.z]);
                text(&mut out, format_args!("{:?}", v.tolerance));
            }
        }
    }
    Some(out)
}

/// A count or a length.
fn count(out: &mut Vec<u8>, n: usize) {
    out.extend_from_slice(&(n as u64).to_le_bytes());
}

/// Text, after its length.
fn text(out: &mut Vec<u8>, args: std::fmt::Arguments<'_>) {
    let written = args.to_string();
    count(out, written.len());
    out.extend_from_slice(written.as_bytes());
}

/// Floats as their bits.
fn bits(out: &mut Vec<u8>, values: &[f64]) {
    for x in values {
        out.extend_from_slice(&x.to_bits().to_le_bytes());
    }
}

/// A placement: the identity as a letter, any other written out.
fn placement(out: &mut Vec<u8>, transform: Transform) {
    if transform == Transform::IDENTITY {
        out.push(b'I');
    } else {
        out.push(b'T');
        text(out, format_args!("{transform:?}"));
    }
}

fn knots(out: &mut Vec<u8>, knots: &KnotVector) {
    count(out, knots.degree());
    count(out, knots.knots().len());
    bits(out, knots.knots());
}

/// A curve in space; a spline by its bits, its control points being the
/// bulk of what a face holds.
fn space_curve(out: &mut Vec<u8>, curve: &Curve) {
    match curve {
        Curve::BSpline(spline) => {
            out.push(b'B');
            out.push(u8::from(spline.is_rational()));
            out.push(u8::from(spline.is_periodic()));
            knots(out, spline.knots());
            count(out, spline.control_points().len());
            for c in spline.control_points() {
                bits(out, &[c.scaled.x, c.scaled.y, c.scaled.z, c.weight]);
            }
        }
        other => {
            out.push(b'D');
            text(out, format_args!("{other:?}"));
        }
    }
}

/// A curve in a chart, a spline by its bits.
fn planar_curve(out: &mut Vec<u8>, curve: &PlanarCurve) {
    match curve {
        PlanarCurve::BSpline(spline) => {
            out.push(b'B');
            out.push(u8::from(spline.is_rational()));
            knots(out, spline.knots());
            count(out, spline.control_points().len());
            for c in spline.control_points() {
                bits(out, &[c.scaled.x, c.scaled.y, c.weight]);
            }
        }
        other => {
            out.push(b'D');
            text(out, format_args!("{other:?}"));
        }
    }
}
