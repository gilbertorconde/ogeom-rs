//! Exact admission of interior nodes on a convex spherical corner. Its
//! great-circle trims are half-spaces through the sphere's centre, whereas
//! their polygonal UV chords can bulge outside the real face.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

#[derive(Debug)]
pub(crate) enum SphereTrim {
    /// Other faces keep the ordinary polygonal region test.
    Other,
    /// Each inward unit normal bounds the actual spherical polygon.
    Convex { centre: Point, inward: Vec<Vector> },
}

impl SphereTrim {
    pub(crate) fn of(
        model: &Model,
        face: &Shape,
        surface: &SurfaceGeometry,
        rings: &[Vec<Point2>],
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let SurfaceGeometry::Sphere(sphere) = surface else {
            return Ok(Self::Other);
        };
        if rings.len() != 1 {
            return Ok(Self::Other);
        }
        let centre = sphere.sphere().frame().origin();
        let to_local = face.transform(model.datums())?.inverse()?;
        let mut points = Vec::new();
        let mut normals = Vec::new();
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            let Some(NodeData::Edge(data)) = model.node(&edge).map(|n| n.data()) else {
                ogeom_bail!(Dangling, "a sphere boundary edge does not resolve");
            };
            if data.degenerate {
                continue;
            }
            let Some(EdgeRepr::Curve3d {
                curve,
                range,
                location,
            }) = data.curve3d()
            else {
                return Ok(Self::Other);
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "a sphere boundary curve does not resolve");
            };
            let curve = geometry
                .transformed(&location.composed(model.datums())?, tol)?
                .transformed(&edge.transform(model.datums())?, tol)?
                .transformed(&to_local, tol)?;
            let mut basis = &curve;
            while let Curve::Trimmed(trimmed) = basis {
                basis = trimmed.basis();
            }
            let Curve::Circle(circle) = basis else {
                return Ok(Self::Other);
            };
            let circle = circle.circle();
            // Minor great-circle arcs bound a convex polygon in one open
            // hemisphere. A complete sphere, a small circle, a hole or a
            // non-convex trim is not this construction.
            if circle.centre().distance(centre) > tol.confusion()
                || (circle.radius() - sphere.sphere().radius()).abs() > tol.confusion()
                || range.1 - range.0 > core::f64::consts::PI + tol.angular()
            {
                return Ok(Self::Other);
            }
            // The oriented wire selects the material side, not an arbitrary
            // hemisphere with the same boundary. Reversing a face reverses
            // both its edge walk and its surface normal.
            let middle = f64::midpoint(range.0, range.1);
            let mut tangent = curve.d1_at(middle, tol)?;
            if edge.orientation() == Orientation::Reversed {
                tangent = -tangent;
            }
            let mut radial = curve.point_at(middle, tol)? - centre;
            if face.orientation() == Orientation::Reversed {
                radial = -radial;
            }
            let n = circle.frame().z().vector();
            normals.push(if radial.cross(tangent).dot(n) >= 0.0 {
                n
            } else {
                -n
            });
            for t in [range.0, f64::midpoint(range.0, range.1), range.1] {
                points.push(curve.point_at(t, tol)?);
            }
        }
        if normals.len() < 3 {
            return Ok(Self::Other);
        }
        let mut inward = Vec::new();
        for n in normals {
            let (lo, hi) = points
                .iter()
                .map(|p| (*p - centre).dot(n))
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), h| {
                    (lo.min(h), hi.max(h))
                });
            if lo >= -tol.confusion() && hi > tol.confusion() {
                inward.push(n);
            } else {
                return Ok(Self::Other);
            }
        }
        Ok(Self::Convex { centre, inward })
    }

    pub(crate) fn holds(
        &self,
        surface: &SurfaceGeometry,
        p: Point2,
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        match self {
            Self::Other => Ok(true),
            Self::Convex { centre, inward } => {
                let p = surface.point_at(p.x, p.y, tol)?;
                Ok(inward
                    .iter()
                    .all(|n| (p - *centre).dot(*n) >= -tol.confusion()))
            }
        }
    }
}
