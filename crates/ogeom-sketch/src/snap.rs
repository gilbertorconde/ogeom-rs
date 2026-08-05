//! Snapping: what the pointer means when it lands near something.
//!
//! A sketch is drawn with a pointer, and a pointer is never exactly on
//! anything. Snapping is the sketch answering "what did you mean" — and the
//! answer has to say *what* it snapped to, not just where, because the
//! caller's next move is usually to constrain against it. A point snapped
//! to a line's midpoint is a midpoint constraint waiting to happen; the
//! same coordinates snapped to nothing are just coordinates.
//!
//! The order is by kind before distance. A defined feature — a point, an
//! end, a centre, a midpoint — beats lying on a curve, and lying on a curve
//! beats the grid, because that is the order of how much the sketch knows
//! about each. Within one kind, the nearest wins.

use ogeom_math::Point2;

use crate::model::{ArcId, CircleId, LineId, PointId, Sketch};

/// What the pointer landed on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SnapKind {
    /// A sketch point itself.
    Point(PointId),
    /// The midpoint of a line.
    Midpoint(LineId),
    /// A circle's centre.
    CircleCentre(CircleId),
    /// An arc's centre.
    ArcCentre(ArcId),
    /// Somewhere along a line, between its ends.
    OnLine(LineId),
    /// Somewhere on a circle's rim.
    OnCircle(CircleId),
    /// A grid crossing.
    Grid,
}

impl SnapKind {
    /// How much the sketch knows about this kind: lower is more definite,
    /// and more definite wins.
    const fn rank(self) -> u8 {
        match self {
            Self::Point(_) => 0,
            Self::Midpoint(_) | Self::CircleCentre(_) | Self::ArcCentre(_) => 1,
            Self::OnLine(_) | Self::OnCircle(_) => 2,
            Self::Grid => 3,
        }
    }
}

/// Where the pointer was taken to mean, and what it meant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Snap {
    /// The snapped position.
    pub at: Point2,
    /// What it snapped to.
    pub kind: SnapKind,
    /// How far the pointer was from it.
    pub distance: f64,
}

/// How near counts as near, and what may be snapped to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SnapOptions {
    /// The reach, in sketch units. Nothing further away is considered.
    pub radius: f64,
    /// Grid spacing, if the grid is one of the things to snap to.
    pub grid: Option<f64>,
    /// Whether construction geometry counts. It usually should: that is
    /// what construction geometry is *for*.
    pub construction: bool,
}

impl Default for SnapOptions {
    fn default() -> Self {
        Self {
            radius: 5.0,
            grid: None,
            construction: true,
        }
    }
}

impl Sketch {
    /// What the pointer at `at` means, given how near it has to be.
    ///
    /// `None` when nothing — not even the grid, if one was offered — lies
    /// within the reach.
    #[must_use]
    pub fn snap(&self, at: Point2, options: SnapOptions) -> Option<Snap> {
        if !options.radius.is_finite() || options.radius <= 0.0 {
            return None;
        }
        let mut best: Option<Snap> = None;
        let mut offer = |candidate: Snap| {
            if candidate.distance > options.radius {
                return;
            }
            let better = best.as_ref().is_none_or(|held| {
                (candidate.kind.rank(), candidate.distance) < (held.kind.rank(), held.distance)
            });
            if better {
                best = Some(candidate);
            }
        };

        for (i, data) in self.points.iter().enumerate() {
            if data.construction && !options.construction {
                continue;
            }
            let p = Point2::new(self.params[data.at], self.params[data.at + 1]);
            offer(Snap {
                at: p,
                kind: SnapKind::Point(PointId(i)),
                distance: p.distance(at),
            });
        }

        for (i, data) in self.lines.iter().enumerate() {
            if data.construction && !options.construction {
                continue;
            }
            let (Ok(a), Ok(b)) = (self.point(data.a), self.point(data.b)) else {
                continue;
            };
            let middle = Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y));
            offer(Snap {
                at: middle,
                kind: SnapKind::Midpoint(LineId(i)),
                distance: middle.distance(at),
            });
            // The foot of the pointer on the segment, kept between the ends
            // — beyond them the line is not there to snap to.
            let along = b - a;
            let length = along.magnitude();
            if length > f64::MIN_POSITIVE {
                let t = ((at - a).dot(along) / (length * length)).clamp(0.0, 1.0);
                let foot = a + along * t;
                offer(Snap {
                    at: foot,
                    kind: SnapKind::OnLine(LineId(i)),
                    distance: foot.distance(at),
                });
            }
        }

        for (i, data) in self.circles.iter().enumerate() {
            if data.construction && !options.construction {
                continue;
            }
            let Ok(centre) = self.point(data.centre) else {
                continue;
            };
            offer(Snap {
                at: centre,
                kind: SnapKind::CircleCentre(CircleId(i)),
                distance: centre.distance(at),
            });
            let radius = self.params[data.radius_at];
            let out = at - centre;
            let reach = out.magnitude();
            if reach > f64::MIN_POSITIVE && radius > 0.0 {
                let rim = centre + out * (radius / reach);
                offer(Snap {
                    at: rim,
                    kind: SnapKind::OnCircle(CircleId(i)),
                    distance: rim.distance(at),
                });
            }
        }

        for (i, data) in self.arcs.iter().enumerate() {
            let Ok(centre) = self.point(data.centre) else {
                continue;
            };
            offer(Snap {
                at: centre,
                kind: SnapKind::ArcCentre(ArcId(i)),
                distance: centre.distance(at),
            });
        }

        if let Some(step) = options.grid
            && step.is_finite()
            && step > 0.0
        {
            let node = Point2::new((at.x / step).round() * step, (at.y / step).round() * step);
            offer(Snap {
                at: node,
                kind: SnapKind::Grid,
                distance: node.distance(at),
            });
        }

        best
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::model::Constraint;

    #[test]
    fn a_defined_feature_beats_a_curve_and_a_curve_beats_the_grid() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(Point2::new(0.0, 0.0));
        let b = sketch.add_point(Point2::new(10.0, 0.0));
        let line = sketch.add_line(a, b).unwrap();
        let options = SnapOptions {
            radius: 2.0,
            grid: Some(1.0),
            construction: true,
        };

        // Near the middle of the line: the midpoint is a feature, so it
        // wins over both lying on the line and the grid crossing, even
        // though all three are within reach.
        let snapped = sketch.snap(Point2::new(5.2, 0.3), options).unwrap();
        assert_eq!(snapped.kind, SnapKind::Midpoint(line));
        assert!((snapped.at.x - 5.0).abs() < 1e-12 && snapped.at.y.abs() < 1e-12);

        // Away from any feature but over the line: on the line.
        let snapped = sketch.snap(Point2::new(2.4, 0.3), options).unwrap();
        assert_eq!(snapped.kind, SnapKind::OnLine(line));
        assert!(snapped.at.y.abs() < 1e-12 && (snapped.at.x - 2.4).abs() < 1e-12);

        // Off the line entirely: the grid, since it was offered.
        let snapped = sketch.snap(Point2::new(2.4, 3.1), options).unwrap();
        assert_eq!(snapped.kind, SnapKind::Grid);
        assert!((snapped.at.x - 2.0).abs() < 1e-12 && (snapped.at.y - 3.0).abs() < 1e-12);

        // And nothing at all, out of reach, with no grid offered.
        assert!(
            sketch
                .snap(
                    Point2::new(40.0, 40.0),
                    SnapOptions {
                        grid: None,
                        ..options
                    }
                )
                .is_none()
        );
    }

    #[test]
    fn a_rim_snaps_to_the_rim_and_construction_can_be_left_out() {
        let mut sketch = Sketch::new();
        let centre = sketch.add_point(Point2::new(0.0, 0.0));
        let circle = sketch.add_circle(centre, 5.0).unwrap();
        let options = SnapOptions {
            radius: 1.0,
            grid: None,
            construction: true,
        };

        let snapped = sketch.snap(Point2::new(5.4, 0.0), options).unwrap();
        assert_eq!(snapped.kind, SnapKind::OnCircle(circle));
        assert!((snapped.at.x - 5.0).abs() < 1e-12);

        // Made construction and excluded, the same pointer finds nothing —
        // and its centre point goes with it, since the point is the
        // circle's own.
        sketch.set_circle_construction(circle, true).unwrap();
        sketch.set_point_construction(centre, true).unwrap();
        assert!(
            sketch
                .snap(
                    Point2::new(5.4, 0.0),
                    SnapOptions {
                        construction: false,
                        ..options
                    }
                )
                .is_none()
        );
        // Kept in, it is found again.
        assert_eq!(
            sketch.snap(Point2::new(5.4, 0.0), options).unwrap().kind,
            SnapKind::OnCircle(circle)
        );
        let _ = Constraint::Radius(circle, 5.0);
    }
}
