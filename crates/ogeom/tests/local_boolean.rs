//! A small cut into a large solid costs what the cut touches, not what the
//! solid holds.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A pin of radius 1.5 standing 20 tall from `at`.
fn pin(model: &mut Model, at: Point) -> Shape {
    let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
    make_cylinder(model, frame, 1.5, 20.0, T).unwrap().shape
}

/// A plate 10 thick drilled through with a `rows` by `rows` grid of holes
/// of radius 1.5, 10 apart and 10 from its sides, in one cut, with the pin
/// drilling one more hole at its corner.
fn plate_and_corner(rows: u32, faces: usize) -> (Model, Shape, Shape) {
    let mut model = Model::new();
    let size = 10.0 * f64::from(rows) + 10.0;
    let block = make_box(&mut model, Frame::WORLD, (size, size, 10.0), T)
        .unwrap()
        .shape;
    let mut pins = Vec::new();
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j), -5.0);
            pins.push(pin(&mut model, at));
        }
    }
    let pins = model.add_compound(&pins).unwrap();
    let plate = ogeom::boolean::cut(&mut model, &block, &pins, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .len(),
        faces
    );
    let corner = pin(&mut model, Point::new(5.0, 5.0, -5.0));
    (model, plate, corner)
}

/// The corner hole into a plate of 582 faces does the work the same hole
/// does in a plate of 22: the faces it does not reach are set aside, not
/// gathered, intersected, split or classified, so the faces gathered, the
/// face pairs asked for sections and the pieces classified are the same
/// counts for both plates.
#[test]
fn a_corner_hole_does_the_work_it_touches_not_what_the_plate_holds() {
    let work = |(model, plate, corner): &mut (Model, Shape, Shape)| {
        let before = explore_unique(model, plate, ShapeType::Face).unwrap().len();
        let (cut, work) =
            ogeom_bool::work::counted(|| ogeom::boolean::cut(model, plate, corner, T));
        let after = explore_unique(model, &cut.unwrap().shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(after, before + 1);
        work
    };
    let small = work(&mut plate_and_corner(4, 22));
    let large = work(&mut plate_and_corner(24, 582));
    assert!(small.faces_gathered > 0 && small.face_pairs > 0 && small.pieces > 0);
    assert_eq!(large, small, "582 faces against 22");
}

/// The processor time of the corner hole into the plate of 582 faces
/// against the plate of 22, reported: the least of many runs each, taken
/// in turns, the cut made again each time in the model holding its plate.
/// What the large plate still pays for every face it passes through
/// (`side`, the history copies, the shells put back together) is linear
/// and memory-bound, and moves with the machine.
#[test]
#[ignore = "heavy"]
fn a_corner_hole_s_time_against_the_plate_it_cuts_is_reported() {
    let mut small = plate_and_corner(4, 22);
    let mut large = plate_and_corner(24, 582);
    let least = |(model, plate, corner): &mut (Model, Shape, Shape)| {
        (0..6).fold(f64::INFINITY, |least, _| {
            let (cut, took) = cpu_time(|| ogeom::boolean::cut(model, plate, corner, T));
            cut.unwrap();
            least.min(took.as_secs_f64())
        })
    };
    let (mut least_small, mut least_large) = (f64::INFINITY, f64::INFINITY);
    for _ in 0..12 {
        least_small = least_small.min(least(&mut small));
        least_large = least_large.min(least(&mut large));
    }
    eprintln!(
        "corner hole: {least_large:.5} s into 582 faces, {least_small:.5} s into 22, ratio {:.2}",
        least_large / least_small
    );
}

/// How many faces of `result` are faces of `plate` passed through as they
/// stand.
fn passed_through(model: &Model, plate: &Shape, result: &Shape) -> usize {
    let own: ogeom::core::FastSet<_> = explore_unique(model, plate, ShapeType::Face)
        .unwrap()
        .iter()
        .map(Shape::node)
        .collect();
    explore_unique(model, result, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|face| own.contains(&face.node()))
        .count()
}

/// The vertices of the hole nearest the corner pin.
fn nearest_hole_vertices(model: &Model, plate: &Shape) -> Vec<Shape> {
    explore_unique(model, plate, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .filter(|v| {
            let at = model.node(v).unwrap().data().as_vertex().unwrap().point;
            at.x < 12.0 && at.y < 12.0 && at.x > 8.0 && at.y > 8.0
        })
        .collect()
}

/// Widen `vertex` through [`Model::widen`].
fn widen_vertex(model: &mut Model, vertex: &Shape) {
    model
        .widen(vertex, ogeom::core::Tolerance::new(1.1).unwrap())
        .unwrap();
}

/// Widen `vertex` through the node handed out for editing.
fn widen_in_node(model: &mut Model, vertex: &Shape) {
    if let Some(node) = model.node_mut(vertex)
        && let ogeom::topo::NodeData::Vertex(data) = node.data_mut()
    {
        data.tolerance = ogeom::core::Tolerance::new(1.1).unwrap();
    }
}

/// A boolean on a solid a boolean has read before reads it back as it was
/// read, and reads it as it stands after an edit in place: a vertex of the hole beside the corner pin,
/// widened until that hole's wall reaches the pin, takes the wall into the
/// cut, as the same cut into the plate edited before any boolean does.
#[test]
fn a_solid_edited_in_place_is_read_as_it_stands_by_the_next_boolean() {
    let editors: [fn(&mut Model, &Shape); 2] = [widen_vertex, widen_in_node];
    let (model, plate, corner) = plate_and_corner(8, 70);
    let mut unedited = model.clone();
    let cut = ogeom::boolean::cut(&mut unedited, &plate, &corner, T).unwrap();
    let untouched = passed_through(&unedited, &plate, &cut.shape);
    for edit in editors {
        // Read by a first cut and read back by a second, then edited.
        let mut read = model.clone();
        ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        let back = ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        assert_eq!(passed_through(&read, &plate, &back.shape), untouched);
        for vertex in nearest_hole_vertices(&read, &plate) {
            edit(&mut read, &vertex);
        }
        let again = ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        // Edited before any boolean read it.
        let mut fresh = model.clone();
        for vertex in nearest_hole_vertices(&fresh, &plate) {
            edit(&mut fresh, &vertex);
        }
        let first = ogeom::boolean::cut(&mut fresh, &plate, &corner, T).unwrap();
        let expected = passed_through(&fresh, &plate, &first.shape);
        assert!(expected < untouched, "{expected} of {untouched}");
        assert_eq!(passed_through(&read, &plate, &again.shape), expected);
    }
}

/// The cells of a 19 by 19 grid 10 apart on a plate 200 square, in an
/// order drawn from a fixed seed.
fn shuffled_cells() -> Vec<(u32, u32)> {
    let mut cells: Vec<(u32, u32)> = (0..19).flat_map(|i| (0..19).map(move |j| (i, j))).collect();
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    for k in (1..cells.len()).rev() {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let pick = usize::try_from(state >> 33).unwrap() % (k + 1);
        cells.swap(k, pick);
    }
    cells
}

/// A plate 200 square and 10 thick, then `steps` booleans on it in turn,
/// each on the result of the one before: a hole of radius 1.5 cut through
/// and a boss of radius 2.5 standing 5 above the top fused on, by turns,
/// each in a cell of its own. `each` is handed the model, the result and
/// the step's processor time after every step.
fn chain(steps: usize, mut each: impl FnMut(&Model, &Shape, usize, f64)) {
    let mut model = Model::new();
    let mut part = make_box(&mut model, Frame::WORLD, (200.0, 200.0, 10.0), T)
        .unwrap()
        .shape;
    for (step, &(i, j)) in shuffled_cells().iter().take(steps).enumerate() {
        let (x, y) = (10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j));
        let hole = step % 2 == 0;
        let tool = if hole {
            pin(&mut model, Point::new(x, y, -5.0))
        } else {
            let frame = Frame::new(Point::new(x, y, 5.0), Direction::Z, Direction::X, T).unwrap();
            make_cylinder(&mut model, frame, 2.5, 10.0, T)
                .unwrap()
                .shape
        };
        let (built, took) = cpu_time(|| {
            if hole {
                ogeom::boolean::cut(&mut model, &part, &tool, T)
            } else {
                ogeom::boolean::fuse(&mut model, &part, &tool, T)
            }
        });
        part = built.unwrap().shape;
        each(&model, &part, step, took.as_secs_f64());
    }
}

/// Each step's processor time along a chain of 200 holes and bosses,
/// reported every 20 steps with the faces the part holds then.
#[test]
#[ignore = "heavy"]
fn a_chain_of_booleans_s_time_per_step_is_reported() {
    let steps = std::env::var("OGEOM_CHAIN_STEPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    // Once untimed, so the timed run starts on a warm core.
    chain(steps, |_, _, _, _| {});
    let mut window = Vec::new();
    let mut total = 0.0;
    chain(steps, |model, part, step, took| {
        window.push(took);
        total += took;
        if window.len() == 20 || step + 1 == steps {
            window.sort_by(f64::total_cmp);
            let faces = explore_unique(model, part, ShapeType::Face).unwrap().len();
            eprintln!(
                "chain steps {:>3}..{:>3}: {faces:>4} faces, median {:.2} ms, least {:.2} ms",
                step + 1 - window.len(),
                step,
                window[window.len() / 2] * 1e3,
                window[0] * 1e3
            );
            window.clear();
        }
    });
    eprintln!("chain of {steps}: {:.3} s", total);
}

/// Along a chain of booleans, each result is read with the faces it shares
/// with its operand taken from what was kept of the operand's read, and
/// each read so taken is the one a fresh read gives: the reads are compared
/// on every step.
#[test]
fn a_chain_of_booleans_reads_each_result_as_a_fresh_read_would() {
    let (faces, taken) = ogeom::algo::derived::derived_reads_checked(|| {
        let mut faces = 0;
        chain(80, |model, part, _, _| {
            faces = explore_unique(model, part, ShapeType::Face).unwrap().len();
        });
        faces
    });
    // The plate's 6 faces, a wall per hole, a wall and a top per boss.
    assert_eq!(faces, 6 + 40 + 80);
    assert!(taken > 0);
}

/// A solid a boolean made, edited in place before the next boolean, is
/// read as it stands: a vertex of the hole beside the corner pin, a hole
/// the first cut passed through from the plate, widened until its wall
/// reaches the pin, takes the wall into the corner cut; nothing read of
/// the plate before the edit is taken for it.
#[test]
fn a_result_edited_in_place_is_read_as_it_stands_by_the_next_boolean() {
    let editors: [fn(&mut Model, &Shape); 2] = [widen_vertex, widen_in_node];
    let (mut model, plate, corner) = plate_and_corner(8, 70);
    let far = pin(&mut model, Point::new(75.0, 75.0, -5.0));
    let holed = ogeom::boolean::cut(&mut model, &plate, &far, T)
        .unwrap()
        .shape;
    let mut unedited = model.clone();
    let (cut, taken) = ogeom::algo::derived::derived_reads_checked(|| {
        ogeom::boolean::cut(&mut unedited, &holed, &corner, T).unwrap()
    });
    assert!(taken > 0);
    let untouched = passed_through(&unedited, &holed, &cut.shape);
    for edit in editors {
        let mut edited = model.clone();
        for vertex in nearest_hole_vertices(&edited, &holed) {
            edit(&mut edited, &vertex);
        }
        let (cut, _) = ogeom::algo::derived::derived_reads_checked(|| {
            ogeom::boolean::cut(&mut edited, &holed, &corner, T).unwrap()
        });
        let through = passed_through(&edited, &holed, &cut.shape);
        assert!(through < untouched, "{through} of {untouched}");
    }
}
