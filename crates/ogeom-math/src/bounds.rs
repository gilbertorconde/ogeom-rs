//! Axis-aligned bounding boxes.
//!
//! Every use of a bounding box in a kernel is a *rejection* test: cull this
//! pair, skip this subtree, prune this branch. So the one property that matters
//! is that a box genuinely contains what it claims to. A box that is too large
//! costs time; a box that is too small silently drops a real intersection, and
//! nothing downstream can tell.
//!
//! Everything that produces an [`Aabb`] here therefore errs outward, and says
//! how. Sampling a curve at a few parameters and taking the extremes would
//! *not* qualify: the curve bulges between the samples.

use core::fmt;

use ogeom_core::Tolerances;

use crate::{Point, Vector};

/// An axis-aligned box, or the empty box.
///
/// The empty box is a distinct state rather than a degenerate one with reversed
/// bounds: "no points at all" and "a box of zero size at the origin" answer
/// containment questions differently, and conflating them makes an empty shape
/// appear to be at the origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// `None` when empty.
    extent: Option<(Point, Point)>,
}

impl Default for Aabb {
    fn default() -> Self {
        Self::EMPTY
    }
}

impl Aabb {
    /// A box containing nothing.
    pub const EMPTY: Self = Self { extent: None };

    /// A box containing exactly one point.
    #[must_use]
    pub const fn of_point(p: Point) -> Self {
        Self {
            extent: Some((p, p)),
        }
    }

    /// A box spanning two corners, in either order.
    #[must_use]
    pub fn of_corners(a: Point, b: Point) -> Self {
        Self {
            extent: Some((a.min(b), a.max(b))),
        }
    }

    /// A box containing every point given.
    #[must_use]
    pub fn of_points(points: &[Point]) -> Self {
        points.iter().fold(Self::EMPTY, |acc, p| acc.with_point(*p))
    }

    /// Whether this box contains nothing.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.extent.is_none()
    }

    /// The lower corner, or `None` if empty.
    #[must_use]
    pub fn low(&self) -> Option<Point> {
        self.extent.map(|(low, _)| low)
    }

    /// The upper corner, or `None` if empty.
    #[must_use]
    pub fn high(&self) -> Option<Point> {
        self.extent.map(|(_, high)| high)
    }

    /// The centre, or `None` if empty.
    #[must_use]
    pub fn centre(&self) -> Option<Point> {
        self.extent.map(|(low, high)| low.midpoint(high))
    }

    /// The size along each axis, or zero if empty.
    #[must_use]
    pub fn size(&self) -> Vector {
        self.extent.map_or(Vector::ZERO, |(low, high)| high - low)
    }

    /// The length of the longest axis, or zero if empty.
    #[must_use]
    pub fn extent_max(&self) -> f64 {
        let s = self.size();
        s.x.max(s.y).max(s.z)
    }

    /// The length of the diagonal, or zero if empty.
    #[must_use]
    pub fn diagonal(&self) -> f64 {
        self.size().magnitude()
    }

    /// The enclosed volume, or zero if empty or flat.
    #[must_use]
    pub fn volume(&self) -> f64 {
        let s = self.size();
        s.x * s.y * s.z
    }

    /// This box grown to include `p`.
    #[must_use]
    pub fn with_point(&self, p: Point) -> Self {
        match self.extent {
            None => Self::of_point(p),
            Some((low, high)) => Self {
                extent: Some((low.min(p), high.max(p))),
            },
        }
    }

    /// This box grown to include `other`.
    #[must_use]
    pub fn union(&self, other: &Self) -> Self {
        match (self.extent, other.extent) {
            (None, _) => *other,
            (_, None) => *self,
            (Some((al, ah)), Some((bl, bh))) => Self {
                extent: Some((al.min(bl), ah.max(bh))),
            },
        }
    }

    /// This box grown by `margin` in every direction.
    ///
    /// The usual way to account for a tolerance before a rejection test: two
    /// shapes whose boxes miss by less than their tolerances may still touch.
    #[must_use]
    pub fn expanded(&self, margin: f64) -> Self {
        match self.extent {
            None => Self::EMPTY,
            Some((low, high)) => {
                let m = Vector::splat(margin.max(0.0));
                Self {
                    extent: Some((low - m, high + m)),
                }
            }
        }
    }

    /// This box grown by the confusion tolerance.
    #[must_use]
    pub fn with_tolerance(&self, tol: Tolerances) -> Self {
        self.expanded(tol.confusion())
    }

    /// Whether `p` lies inside, boundary included.
    #[must_use]
    pub fn contains(&self, p: Point) -> bool {
        self.extent.is_some_and(|(low, high)| {
            p.x >= low.x
                && p.x <= high.x
                && p.y >= low.y
                && p.y <= high.y
                && p.z >= low.z
                && p.z <= high.z
        })
    }

    /// Whether `other` lies entirely inside this box.
    #[must_use]
    pub fn contains_box(&self, other: &Self) -> bool {
        match other.extent {
            // Nothing is contained by anything, vacuously — including by the
            // empty box.
            None => true,
            Some((low, high)) => self.contains(low) && self.contains(high),
        }
    }

    /// Whether the two boxes share any point.
    ///
    /// The rejection test everything else is built on. An empty box intersects
    /// nothing.
    #[must_use]
    pub fn intersects(&self, other: &Self) -> bool {
        let (Some((al, ah)), Some((bl, bh))) = (self.extent, other.extent) else {
            return false;
        };
        al.x <= bh.x && ah.x >= bl.x && al.y <= bh.y && ah.y >= bl.y && al.z <= bh.z && ah.z >= bl.z
    }

    /// The overlap of two boxes, or empty if they do not meet.
    #[must_use]
    pub fn intersection(&self, other: &Self) -> Self {
        if !self.intersects(other) {
            return Self::EMPTY;
        }
        let (Some((al, ah)), Some((bl, bh))) = (self.extent, other.extent) else {
            return Self::EMPTY;
        };
        Self {
            extent: Some((al.max(bl), ah.min(bh))),
        }
    }

    /// The shortest distance from `p` to this box, zero if inside.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        let Some((low, high)) = self.extent else {
            return f64::INFINITY;
        };
        // Per axis, how far outside the slab the point lies. Zero inside.
        let outside = Vector::new(
            (low.x - p.x).max(p.x - high.x).max(0.0),
            (low.y - p.y).max(p.y - high.y).max(0.0),
            (low.z - p.z).max(p.z - high.z).max(0.0),
        );
        outside.magnitude()
    }

    /// A lower bound on the distance between two boxes, zero if they meet.
    ///
    /// A *bound*, not the distance between the shapes inside them, which is why
    /// it is only ever useful for rejection: if this exceeds a threshold the
    /// shapes certainly do too.
    #[must_use]
    pub fn distance_to_box(&self, other: &Self) -> f64 {
        let (Some((al, ah)), Some((bl, bh))) = (self.extent, other.extent) else {
            return f64::INFINITY;
        };
        let gap = Vector::new(
            (bl.x - ah.x).max(al.x - bh.x).max(0.0),
            (bl.y - ah.y).max(al.y - bh.y).max(0.0),
            (bl.z - ah.z).max(al.z - bh.z).max(0.0),
        );
        gap.magnitude()
    }

    /// The eight corners, or an empty list if the box is empty.
    #[must_use]
    pub fn corners(&self) -> Vec<Point> {
        let Some((low, high)) = self.extent else {
            return Vec::new();
        };
        let mut out = Vec::with_capacity(8);
        for i in 0..8 {
            out.push(Point::new(
                if i & 1 == 0 { low.x } else { high.x },
                if i & 2 == 0 { low.y } else { high.y },
                if i & 4 == 0 { low.z } else { high.z },
            ));
        }
        out
    }

    /// This box transformed, as the box of the transformed corners.
    ///
    /// The result contains the transformed box but is generally larger than the
    /// tightest one: a rotated box is not axis-aligned, and its bounding box
    /// must cover the rotation. Erring outward is the safe direction, and
    /// repeatedly transforming a box therefore inflates it — transform the
    /// geometry and re-bound instead of chaining this.
    #[must_use]
    pub fn transformed(&self, t: &crate::Transform) -> Self {
        Self::of_points(
            &self
                .corners()
                .iter()
                .map(|p| t.apply(*p))
                .collect::<Vec<_>>(),
        )
    }
}

impl fmt::Display for Aabb {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.extent {
            None => f.write_str("empty"),
            Some((low, high)) => write!(
                f,
                "[{:.6}, {:.6}, {:.6}] .. [{:.6}, {:.6}, {:.6}]",
                low.x, low.y, low.z, high.x, high.y, high.z
            ),
        }
    }
}

impl FromIterator<Point> for Aabb {
    fn from_iter<I: IntoIterator<Item = Point>>(iter: I) -> Self {
        iter.into_iter()
            .fold(Self::EMPTY, |acc, p| acc.with_point(p))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{Axis, Transform};
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    fn unit() -> Aabb {
        Aabb::of_corners(Point::ORIGIN, Point::new(1.0, 1.0, 1.0))
    }

    #[test]
    fn the_empty_box_is_distinct_from_a_zero_sized_one() {
        // Conflating them makes an empty shape appear to be at the origin, and
        // every containment question about it answers wrongly.
        let empty = Aabb::EMPTY;
        let degenerate = Aabb::of_point(Point::ORIGIN);

        assert!(empty.is_empty());
        assert!(!degenerate.is_empty());
        assert!(!empty.contains(Point::ORIGIN));
        assert!(degenerate.contains(Point::ORIGIN));
        assert_eq!(empty.centre(), None);
        assert_eq!(degenerate.centre(), Some(Point::ORIGIN));
        assert!(!empty.intersects(&degenerate));
        assert_eq!(empty.distance_to(Point::ORIGIN), f64::INFINITY);
    }

    #[test]
    fn corners_are_taken_in_either_order() {
        let a = Aabb::of_corners(Point::new(1.0, 2.0, 3.0), Point::ORIGIN);
        let b = Aabb::of_corners(Point::ORIGIN, Point::new(1.0, 2.0, 3.0));
        assert_eq!(a, b);
        assert_eq!(a.low(), Some(Point::ORIGIN));
        assert_eq!(a.high(), Some(Point::new(1.0, 2.0, 3.0)));
    }

    #[test]
    fn union_grows_and_the_empty_box_is_its_identity() {
        let a = Aabb::of_point(Point::ORIGIN);
        let b = Aabb::of_point(Point::new(1.0, 2.0, 3.0));
        let joined = a.union(&b);
        assert_eq!(joined.low(), Some(Point::ORIGIN));
        assert_eq!(joined.high(), Some(Point::new(1.0, 2.0, 3.0)));

        assert_eq!(a.union(&Aabb::EMPTY), a);
        assert_eq!(Aabb::EMPTY.union(&a), a);
        assert_eq!(Aabb::EMPTY.union(&Aabb::EMPTY), Aabb::EMPTY);
    }

    #[test]
    fn measurements_of_a_unit_box() {
        let b = unit();
        assert_eq!(b.size(), Vector::new(1.0, 1.0, 1.0));
        assert_relative_eq!(b.volume(), 1.0);
        assert_relative_eq!(b.extent_max(), 1.0);
        assert_relative_eq!(b.diagonal(), 3.0_f64.sqrt());
        assert_eq!(b.centre(), Some(Point::new(0.5, 0.5, 0.5)));
        assert_eq!(b.corners().len(), 8);
        assert!(Aabb::EMPTY.corners().is_empty());
        assert_relative_eq!(Aabb::EMPTY.volume(), 0.0);
    }

    #[test]
    fn containment_includes_the_boundary() {
        // A rejection test that excluded the boundary would drop shapes that
        // touch exactly, which is the case that matters most.
        let b = unit();
        assert!(b.contains(Point::ORIGIN));
        assert!(b.contains(Point::new(1.0, 1.0, 1.0)));
        assert!(b.contains(Point::new(0.5, 0.0, 1.0)));
        assert!(!b.contains(Point::new(1.000_001, 0.5, 0.5)));
    }

    #[test]
    fn boxes_that_touch_are_reported_as_intersecting() {
        // Two shapes whose boxes meet exactly may well touch, and rejecting
        // them would drop a real contact.
        let a = unit();
        let touching = Aabb::of_corners(Point::new(1.0, 0.0, 0.0), Point::new(2.0, 1.0, 1.0));
        let apart = Aabb::of_corners(Point::new(1.001, 0.0, 0.0), Point::new(2.0, 1.0, 1.0));

        assert!(a.intersects(&touching));
        assert!(!a.intersects(&apart));
        assert_relative_eq!(a.distance_to_box(&touching), 0.0);
        assert_relative_eq!(a.distance_to_box(&apart), 0.001, epsilon = 1e-12);
    }

    #[test]
    fn intersection_is_the_overlap_and_empty_when_there_is_none() {
        let a = unit();
        let b = Aabb::of_corners(Point::new(0.5, 0.5, 0.5), Point::new(2.0, 2.0, 2.0));
        let overlap = a.intersection(&b);
        assert_eq!(overlap.low(), Some(Point::new(0.5, 0.5, 0.5)));
        assert_eq!(overlap.high(), Some(Point::new(1.0, 1.0, 1.0)));

        let apart = Aabb::of_corners(Point::new(5.0, 5.0, 5.0), Point::new(6.0, 6.0, 6.0));
        assert!(a.intersection(&apart).is_empty());
        assert!(a.intersection(&Aabb::EMPTY).is_empty());
    }

    #[test]
    fn distance_is_zero_inside_and_measured_from_the_nearest_face() {
        let b = unit();
        assert_relative_eq!(b.distance_to(Point::new(0.5, 0.5, 0.5)), 0.0);
        assert_relative_eq!(
            b.distance_to(Point::new(0.5, 0.5, 1.0)),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(b.distance_to(Point::new(0.5, 0.5, 3.0)), 2.0);
        // Diagonally outside: the distance combines every axis it is outside on.
        assert_relative_eq!(
            b.distance_to(Point::new(-3.0, -4.0, 0.5)),
            5.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn expansion_errs_outward_and_never_shrinks() {
        let b = unit();
        let grown = b.expanded(0.5);
        assert_eq!(grown.low(), Some(Point::new(-0.5, -0.5, -0.5)));
        assert!(grown.contains_box(&b));

        // A negative margin would shrink the box, which for a rejection test
        // means dropping shapes that are really there.
        assert_eq!(b.expanded(-1.0), b);
        assert!(Aabb::EMPTY.expanded(1.0).is_empty());
        assert!(b.with_tolerance(T).contains_box(&b));
    }

    #[test]
    fn transforming_a_box_errs_outward() {
        // A rotated box is not axis-aligned, so its bound must cover the
        // rotation. That means the result is larger than the tightest box round
        // the rotated shape — which is the safe direction, and the reason to
        // re-bound the geometry rather than chain this.
        let b = unit();
        let rotated = b.transformed(&Transform::rotation(Axis::Z, core::f64::consts::FRAC_PI_4));

        for corner in b.corners() {
            let moved = Transform::rotation(Axis::Z, core::f64::consts::FRAC_PI_4).apply(corner);
            assert!(
                rotated.contains(moved),
                "corner {moved:?} escaped the bound"
            );
        }
        assert!(
            rotated.volume() > b.volume(),
            "a rotation cannot tighten a bound"
        );

        // A translation, by contrast, is exact.
        let moved = b.transformed(&Transform::translation(Vector::new(10.0, 0.0, 0.0)));
        assert_relative_eq!(moved.volume(), b.volume(), epsilon = 1e-12);
    }

    #[test]
    fn the_empty_box_is_contained_by_everything() {
        assert!(unit().contains_box(&Aabb::EMPTY));
        assert!(Aabb::EMPTY.contains_box(&Aabb::EMPTY));
        assert!(!Aabb::EMPTY.contains_box(&unit()));
    }

    #[test]
    fn collecting_points_builds_their_bound() {
        let points = [
            Point::new(1.0, 0.0, 0.0),
            Point::new(-2.0, 5.0, 1.0),
            Point::new(0.0, -1.0, 3.0),
        ];
        let from_iter: Aabb = points.iter().copied().collect();
        assert_eq!(from_iter, Aabb::of_points(&points));
        assert_eq!(from_iter.low(), Some(Point::new(-2.0, -1.0, 0.0)));
        assert_eq!(from_iter.high(), Some(Point::new(1.0, 5.0, 3.0)));
        for p in points {
            assert!(from_iter.contains(p));
        }
        assert!(Aabb::of_points(&[]).is_empty());
    }

    #[test]
    fn display_distinguishes_empty_from_degenerate() {
        assert_eq!(Aabb::EMPTY.to_string(), "empty");
        assert!(
            Aabb::of_point(Point::ORIGIN)
                .to_string()
                .contains("0.000000")
        );
    }
}
