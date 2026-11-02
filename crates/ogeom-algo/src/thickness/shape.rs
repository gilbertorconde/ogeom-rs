//! Wall thickness against a shape's own surfaces.

use ogeom_core::{OgeomResult, parallel};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

use super::{Method, Thickness, inward, meets, spread, through};
use crate::classify::{SolidBoundary, exact_ring_chord};
use crate::distance::Target;

/// The wall thickness at each point of `shape`, a solid or a shell,
/// measured into the material against the point's outward normal on the
/// faces' own surfaces; in the points' order, `None` where it is thicker
/// than [`Thickness::max`] or the inward ray leaves the shape without
/// meeting a face.
///
/// A ray meets a face where it crosses the face's surface inside its trim,
/// or within the band of the trimming rings drawn at ten thousand times the
/// confusion distance, as [`classify_in_solid_exact`](crate::classify_in_solid_exact)
/// draws them. A sphere's nearest surface is
/// [`to_shape`](crate::distance::to_shape)'s nearest point, among the
/// faces, edges and vertices not within the facing limit of edge-on to the
/// point's normal. See the [module](super) for the rest.
///
/// A point that is not finite, or whose normal has no direction, answers
/// `None`.
///
/// # Errors
///
/// As [`Thickness::check`];
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if there is
/// not one normal per point;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the shape is not a solid, a shell or a compound of solids; and as the
/// intersections and projections report.
pub fn thickness_of_shape(
    model: &Model,
    shape: &Shape,
    points: &[Point],
    normals: &[Vector],
    options: &Thickness,
) -> OgeomResult<Vec<Option<f64>>> {
    let sin = options.checked(points.len(), normals.len())?;
    let tol = options.tol;
    let boundary = SolidBoundary::of_closed(model, shape, exact_ring_chord(tol), tol)?;
    let target = match options.method {
        Method::Sphere => Some(Target::of(model, shape, false, tol)?),
        Method::Ray { .. } => None,
    };
    let walls = Walls {
        model,
        boundary: &boundary,
        target: target.as_ref(),
        sin,
        options,
    };
    let answers = parallel::map_ordered(points, |i, &p| match inward(p, normals[i]) {
        Some(u) => walls.at(p, u),
        None => Ok(None),
    });
    answers.into_iter().collect()
}

/// A shape measured against, with what every point asks of it.
struct Walls<'a> {
    model: &'a Model,
    boundary: &'a SolidBoundary,
    /// The shape's elements, for a sphere's nearest points.
    target: Option<&'a Target>,
    sin: f64,
    options: &'a Thickness,
}

impl Walls<'_> {
    /// The wall at `p`, inward along the unit `u`.
    fn at(&self, p: Point, u: Vector) -> OgeomResult<Option<f64>> {
        let tol = self.options.tol;
        let on = self.boundary.faces_at(self.model, p, tol)?;
        let max = self.options.max;
        match (self.options.method, self.target) {
            (Method::Ray { cone, rays }, _) => {
                let mut least: Option<f64> = None;
                for k in 0..rays {
                    let v = spread(u, cone, rays, k);
                    if let Some(length) = self.ray(p, v, max, &on)?
                        && least.is_none_or(|held| length < held)
                    {
                        least = Some(length);
                    }
                }
                Ok(least)
            }
            (Method::Sphere, Some(target)) => {
                let ray = self.ray(p, u, max, &on)?;
                self.sphere(target, p, u, ray)
            }
            (Method::Sphere, None) => Ok(None),
        }
    }

    fn ray(&self, p: Point, v: Vector, max: f64, on: &[Shape]) -> OgeomResult<Option<f64>> {
        self.boundary
            .first_crossing(self.model, (p, v), max, on, self.sin, self.options.tol)
    }

    /// The diameter of the medial ball at `p` along the unit `u`, shrunk
    /// from half the inward ray's reading `ray`, or from half the thickest
    /// wall where the ray met nothing within it.
    ///
    /// Each step finds the surface nearest the centre within the ball and
    /// shrinks to the ball through it, until the nearest surface is within
    /// the confusion distance of the ball.
    fn sphere(
        &self,
        target: &Target,
        p: Point,
        u: Vector,
        ray: Option<f64>,
    ) -> OgeomResult<Option<f64>> {
        let tol = self.options.tol;
        let mut radius = match ray {
            Some(length) => length / 2.0,
            None if self.options.max.is_finite() => self.options.max / 2.0,
            None => return Ok(None),
        };
        let mut bounded = ray.is_some();
        let faces = |n: Vector| meets(n, u, self.sin);
        loop {
            let centre = p + u * radius;
            let Some(found) = target.nearest(self.model, centre, Some(&faces), radius, tol)? else {
                break;
            };
            bounded = true;
            if radius - found.distance <= tol.confusion() {
                break;
            }
            match through(p, u, found.point) {
                Some(r) if r < radius => radius = r,
                _ => break,
            }
        }
        Ok(bounded.then_some(2.0 * radius))
    }
}
