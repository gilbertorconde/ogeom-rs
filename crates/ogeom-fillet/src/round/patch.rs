//! P3's general vertex blend for rounds: one n-sided G1 patch through the
//! strips' end sections and the faces' setback curves (ADR 0001, §4.4).
//!
//! The hole's boundary and region are given to the patch: each side is
//! a strip's end section, which the patch meets tangent to the strip, or a
//! face's setback curve, which it meets tangent to the face. The patch is
//! fitted once, at the fixed net [`VERTEX_NET`], and built whatever it
//! achieved; its residual (each side's gap and angle, and how far it stands
//! past its region) is reported for P5 to judge. It is not refined until a
//! check passes.

use crate::fill_n::{FillBoundary, FillSide, Region, make_filling_n_fixed};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape};

/// The control count across a vertex patch's larger extent.
pub(crate) const VERTEX_NET: usize = 27;

/// A vertex patch and what it achieved.
#[derive(Debug, Clone)]
pub(crate) struct VertexPatch {
    /// The patch's face.
    pub face: Shape,
    /// Each side's report, in the order the sides were given.
    pub sides: Vec<FillSide>,
    /// How far the patch stands past its region at most, inside the hole.
    pub region_past: f64,
}

impl VertexPatch {
    /// The widest gap between the patch and a side's edge.
    pub(crate) fn worst_gap(&self) -> f64 {
        self.sides.iter().map(|s| s.gap).fold(0.0, f64::max)
    }

    /// The widest angle between the patch and a side's support.
    pub(crate) fn worst_angle(&self) -> f64 {
        self.sides
            .iter()
            .filter_map(|s| s.angle)
            .fold(0.0, f64::max)
    }
}

/// The vertex patch over the hole `boundary` bounds, held to `region`,
/// fitted once at [`VERTEX_NET`].
///
/// # Errors
///
/// As [`crate::make_filling_n`]'s reading of the boundary: a loop that does
/// not close, or that crosses itself seen along the plane it spans.
pub(crate) fn vertex_patch(
    model: &mut Model,
    boundary: &[FillBoundary],
    region: Region<'_>,
    tol: Tolerances,
) -> OgeomResult<VertexPatch> {
    let (filled, region_past) =
        make_filling_n_fixed(model, boundary, &[], VERTEX_NET, Some(region), tol)?;
    Ok(VertexPatch {
        face: filled.built.shape,
        sides: filled.sides,
        region_past,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::vertex_patch;
    use crate::fill_n::FillBoundary;
    use crate::round::section::{PlanarSection, RoundSpec};
    use crate::round::sphere::SphereCorner;
    use crate::support::Seat;
    use ogeom_core::Tolerances;
    use ogeom_geom::{Continuity, Surface as _};
    use ogeom_math::{Point, Vector};
    use ogeom_topo::{Model, Shape};

    const T: Tolerances = Tolerances::millimetres();

    /// One strip's face over `(s0, s1)` of its crease, its end arc at `s0`
    /// returned to bound the vertex's hole.
    fn strip(m: &mut Model, section: &PlanarSection, (s0, s1): (f64, f64)) -> (Shape, Shape) {
        let [first, second] = section.order();
        let at = |m: &mut Model, s: f64, i: usize| {
            ogeom_algo::make_vertex(m, section.contact(i, s)).shape
        };
        let (a0, b0, a1, b1) = (
            at(m, s0, first),
            at(m, s0, second),
            at(m, s1, first),
            at(m, s1, second),
        );
        let arc = |m: &mut Model, s: f64, from: &Shape, to: &Shape| {
            let (curve, range) = section.arc_at(s, T).unwrap();
            ogeom_algo::make_edge_between(m, curve, range, from, to, T)
                .unwrap()
                .shape
        };
        let near = arc(m, s0, &a0, &b0);
        let far = arc(m, s1, &a1, &b1);
        let rail = |m: &mut Model, i: usize, from: &Shape, to: &Shape| {
            crate::support::segment_between(
                m,
                (from, section.contact(i, s0)),
                (to, section.contact(i, s1)),
                T,
            )
            .unwrap()
        };
        let ra = rail(m, first, &a0, &a1);
        let rb = rail(m, second, &b0, &b1);
        let surface = section.surface((s0, s1), T).unwrap();
        let face = crate::support::face_from_edges(
            m,
            surface,
            &[near.clone(), rb, far.reversed(), ra.reversed()],
            T,
        )
        .unwrap();
        (face, near)
    }

    /// A box corner's three equal strips end on the corner ball's great
    /// circles. Through their end arcs, tangent to the strips, the fixed
    /// patch meets each strip within the tangency target and stands between
    /// the ball and the corner: outside the ball but for the fit near its
    /// sides, and inside the corner's planes. Being the patch of least
    /// bending, not the ball, it stands off the ball inside; the ball is
    /// the exact case the sphere vertex takes.
    #[test]
    fn the_fixed_patch_through_three_equal_strips_stands_between_ball_and_corner() {
        let r = 2.0;
        let mut m = Model::new();
        let v = ogeom_algo::make_vertex(&mut m, Point::ORIGIN).shape;
        let seat = |start: Point, end: Point, normals: [Vector; 2]| Seat {
            start,
            along: (end - start) / start.distance(end),
            end,
            normals,
            faces: [v.clone(), v.clone()],
            convex: true,
        };
        let (z, y, x) = (
            Vector::new(0.0, 0.0, 1.0),
            Vector::new(0.0, -1.0, 0.0),
            Vector::new(-1.0, 0.0, 0.0),
        );
        let corner = Point::new(0.0, 0.0, 10.0);
        let seats = [
            seat(corner, Point::new(20.0, 0.0, 10.0), [z, y]),
            seat(corner, Point::new(0.0, 20.0, 10.0), [z, x]),
            seat(corner, Point::new(0.0, 0.0, 0.0), [y, x]),
        ];
        let mut boundary = Vec::new();
        for s in &seats {
            let section = PlanarSection::new(s, &RoundSpec::Circular(r), T).unwrap();
            let (face, near) = strip(&mut m, &section, (r, 6.0));
            boundary.push(FillBoundary {
                edge: near,
                support: Some(face),
                continuity: Continuity::G1,
            });
        }
        let planes = [(z, corner), (y, corner), (x, corner)];
        let region = move |p: Point| -> (f64, Vector) {
            let mut worst = (f64::NEG_INFINITY, planes[0].0);
            for (n, on) in planes {
                let d = n.dot(p - on);
                if d > worst.0 {
                    worst = (d, n);
                }
            }
            worst
        };
        let patch = vertex_patch(&mut m, &boundary, &region, T).unwrap();
        let ball = SphereCorner::new(planes, r, true, T).unwrap().unwrap();
        let data = m.node(&patch.face).unwrap().data().as_face().unwrap();
        let surface = m.geometry().surface(data.surface).unwrap().clone();
        let ((u0, u1), (v0, v1)) = surface.domain();
        let fine = ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        };
        let mut off = 0.0_f64;
        let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
        let mut read = 0;
        for i in 0..=20 {
            for j in 0..=20 {
                let p = surface
                    .point_at(
                        u0 + (u1 - u0) * f64::from(i) / 20.0,
                        v0 + (v1 - v0) * f64::from(j) / 20.0,
                        T,
                    )
                    .unwrap();
                let inside = ogeom_algo::classify_on_face(&m, &patch.face, p, fine, T)
                    .is_ok_and(|c| c == ogeom_algo::Containment::In);
                if inside {
                    off = off.max((p.distance(ball.centre) - r).abs());
                    least = least.min(p.distance(ball.centre) - r);
                    most = most.max(region(p).0);
                    read += 1;
                }
            }
        }
        eprintln!(
            "PATCH gap {:.2e} angle {:.2e} past {:.2e} off the ball {off:.2e} least {least:.2e} most {most:.2e}",
            patch.worst_gap(),
            patch.worst_angle(),
            patch.region_past
        );
        assert!(
            patch.worst_gap() <= 2e-4,
            "a side stands {} off",
            patch.worst_gap()
        );
        assert!(
            patch.worst_angle() <= 2e-4,
            "a side meets its strip at {}",
            patch.worst_angle()
        );
        assert!(
            patch.region_past <= 2e-4,
            "the patch stands {} past the corner",
            patch.region_past
        );
        assert!(read > 50, "only {read} samples inside the patch");
        assert!(least >= -4e-4, "the patch stands {least} inside the ball");
        assert!(
            most <= 2e-4,
            "the patch stands {most} past the corner's planes"
        );
    }
}
