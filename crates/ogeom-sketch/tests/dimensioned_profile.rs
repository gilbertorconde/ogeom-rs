//! A dimensioned profile solves to its dimensions, and a contradictory one
//! is refused by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_math::Point2;
use ogeom_sketch::{Constraint, Sketch, SolveOptions};

/// A bracket profile solves to its dimensions and reports itself exactly
/// constrained — and when a contradictory dimension is added, the solver
/// refuses and *names* the two constraints that fight. Naming them is the
/// requirement: a solver that only says "did not converge" leaves the person
/// holding it to find the contradiction by bisection.
#[test]
fn a_sketch_solves_and_an_overconstrained_one_names_the_fight() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(Point2::new(0.1, -0.3));
    let p1 = sketch.add_point(Point2::new(52.0, 2.0));
    let p2 = sketch.add_point(Point2::new(49.0, 21.0));
    let p3 = sketch.add_point(Point2::new(-2.0, 19.0));
    let bottom = sketch.add_line(p0, p1).unwrap();
    let right = sketch.add_line(p1, p2).unwrap();
    let top = sketch.add_line(p3, p2).unwrap();
    let left = sketch.add_line(p0, p3).unwrap();
    let centre = sketch.add_point(Point2::new(20.0, 8.0));
    let bore = sketch.add_circle(centre, 4.0).unwrap();

    sketch
        .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
        .unwrap();
    sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
    sketch
        .constrain(Constraint::Distance(p0, p1, 50.0))
        .unwrap();
    sketch.constrain(Constraint::Vertical(right)).unwrap();
    sketch
        .constrain(Constraint::Distance(p1, p2, 20.0))
        .unwrap();
    sketch.constrain(Constraint::Horizontal(top)).unwrap();
    sketch.constrain(Constraint::Vertical(left)).unwrap();
    sketch.constrain(Constraint::Radius(bore, 5.0)).unwrap();
    sketch
        .constrain(Constraint::Fixed(centre, Point2::new(25.0, 10.0)))
        .unwrap();

    let solution = sketch.solve(SolveOptions::default()).unwrap();
    assert!(solution.converged, "residual {}", solution.residual);
    assert!(solution.diagnosis.is_well_constrained());
    assert_eq!(sketch.measure_distance(p0, p1).unwrap().round(), 50.0);
    assert_eq!(sketch.measure_radius(bore).unwrap().round(), 5.0);

    // Now the fight: a second, disagreeing width.
    let conflicting = sketch
        .constrain(Constraint::Distance(p0, p1, 55.0))
        .unwrap();
    let fought = sketch.solve(SolveOptions::default()).unwrap();
    assert!(!fought.converged, "a contradiction must not converge");
    assert_eq!(fought.diagnosis.conflicting.len(), 1);
    let group = &fought.diagnosis.conflicting[0];
    assert!(group.contains(&conflicting), "the new dimension is named");
    assert!(group.len() >= 2, "so is what it fights: {group:?}");
    for id in group {
        // Every named constraint prints as a sentence a person can act on.
        assert!(!sketch.describe(*id).unwrap().is_empty());
    }
}
