//! Splitting a face by section curves, in its own parameter space.
//!
//! This is the builder's 2D half, and it deliberately knows nothing about
//! geometry beyond polylines. The caller (the pave filler) has already done
//! the exact work: every boundary edge and every section curve is split at
//! every mutual crossing by the intersectors, so what arrives here is
//! *strands*: polyline scaffolding in parameter space, each carrying a tag
//! naming the exact sub-curve it stands for, meeting other strands only at
//! endpoints. This module decides the combinatorics (which strands bound
//! which region) and hands back pieces as sequences of directed tags, so the
//! rebuild works from the exact geometry and the polylines are never promoted
//! to an answer.

use ogeom_core::predicates::{Exact, Predicates, Sign};
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point2;

/// One pre-split piece of boundary or section, as scaffolding plus a name.
#[derive(Debug, Clone)]
pub(crate) struct Strand<T> {
    /// The curve's course through parameter space, finely enough sampled
    /// that its first few snaps from each end give the way it leaves.
    pub polyline: Vec<Point2>,
    /// Which exact sub-curve this stands for.
    pub tag: T,
    /// Whether this is a piece of the face's own boundary, as opposed to a
    /// section. Boundary strands define the material. A dangling boundary is
    /// an error where a dangling section is a pruning.
    pub boundary: bool,
}

/// One directed traversal of a strand inside a piece's ring.
#[derive(Debug, Clone)]
pub(crate) struct Traversal<T> {
    /// The strand's tag.
    pub tag: T,
    /// Whether the ring runs the strand backwards.
    pub reversed: bool,
}

/// One piece of a split face.
pub(crate) struct Piece<'s, T> {
    /// The boundary as directed strand traversals: ring `[0]` is the outer
    /// contour, counter-clockwise in parameter space. Further rings are
    /// holes, clockwise.
    pub rings: Vec<Vec<Traversal<T>>>,
    /// The same rings as parameter-space polylines, in the same order.
    ///
    /// The caller needs the region itself, not only its name: deciding
    /// whether two coincident faces' pieces stand for the *same* patch of one
    /// surface is a containment question, and containment is asked of an
    /// outline.
    pub outlines: Vec<Vec<Point2>>,
    /// The roomiest point strictly inside the piece: the first of
    /// [`Piece::interiors`].
    pub first: Point2,
    /// The scanlines the rest of the interior points are sought along, and
    /// the material they are held to, until they are asked for.
    probes: Probes<'s>,
}

/// What finds a piece's interior points, and the points once found.
struct Probes<'s> {
    scanlines: Scanlines,
    material: std::rc::Rc<Boxed<'s>>,
    snap: f64,
    found: Option<Vec<Point2>>,
}

impl<T> Piece<'_, T> {
    /// Points strictly inside the piece, best first.
    ///
    /// More than one, because a single probe can be unlucky: a piece that
    /// merely *touches* the other solid has a probe on that contact reading
    /// neither in nor out, and the way past it is to ask somewhere else in
    /// the same piece. Found when first asked for: most pieces are settled
    /// by [`Piece::first`] alone.
    pub(crate) fn interiors(&mut self) -> &[Point2] {
        let probes = &mut self.probes;
        probes.found.get_or_insert_with(|| {
            // Only probes inside the material this arrangement bounds are
            // of any use to a caller asking "where does this piece stand".
            probes
                .scanlines
                .points(probes.snap)
                .into_iter()
                .filter(|p| probes.material.inside(*p))
                .collect()
        })
    }
}

/// Assemble pre-split strands into the pieces they bound.
///
/// Dangling sections (chains that separate no material) are pruned;
/// regions outside the boundary strands' material (a hole's inside, say) are
/// dropped by an even-odd test against the boundary polylines.
///
/// `places` gives each strand's place among all its face's strands, and
/// `lone` the face's hole rings walked on their own (see [`Lone`]): each is
/// taken back into the piece it lies in, as if walked with the rest.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a boundary
/// strand dangles, or the graph yields no piece at all.
pub(crate) fn assemble<'s, T: Clone>(
    strands: &'s [Strand<T>],
    places: &[usize],
    snap: f64,
    lone: &'s [Lone<T>],
) -> OgeomResult<Vec<Piece<'s, T>>> {
    let walk = walk(strands, snap)?;
    let Walk {
        live,
        from,
        ends,
        cycles,
    } = &walk;
    let tail = |d: usize| -> usize {
        let (u, v) = ends[d / 2];
        if d.is_multiple_of(2) { u } else { v }
    };
    let outline = |cycle: &[usize]| outline_of(live, cycle);
    // Where a cycle stands among the face's: by its first dart's strand
    // and way, as the walk met them.
    let order = |cycle: &[usize]| (places[from[cycle[0] / 2]], cycle[0] % 2 == 1);

    let mut positives: Vec<(&Vec<usize>, Vec<Point2>)> = Vec::new();
    let mut negatives: Vec<HoleSide<'_, T>> = Vec::new();
    for cycle in cycles {
        let line = outline(cycle);
        let a = area(&line);
        if a > snap * snap {
            positives.push((cycle, line));
        } else if a < -(snap * snap) {
            negatives.push((order(cycle), Ok(cycle), line));
        }
    }
    if !lone.is_empty() {
        negatives.extend(
            lone.iter()
                .map(|ring| (ring.order, Err(ring), ring.outline.clone())),
        );
        negatives.sort_by_key(|(at, ..)| *at);
    }

    // The boundary strands' polylines, for the material test, each with
    // its box: a face with hundreds of holes asks it of every hole's disc.
    let material: Vec<&'s [Point2]> = live
        .iter()
        .filter(|s| s.boundary)
        .map(|s| s.polyline.as_slice())
        .chain(
            lone.iter()
                .flat_map(|ring| ring.lines.iter().map(Vec::as_slice)),
        )
        .collect();
    let material = std::rc::Rc::new(Boxed::new(material));

    // The nodes each cycle passes through: a hole that shares one with a
    // positive cycle is the same component, not a hole in it. Asked of the
    // nodes, not of the polylines' nearness: a hole can pass within the weld
    // of a loose boundary without meeting it, and read by distance it is
    // taken for part of the boundary and dropped.
    let nodes_of = |cycle: &[usize]| -> std::collections::BTreeSet<usize> {
        cycle.iter().map(|&d| tail(d)).collect()
    };
    // Each cycle's nodes, area and box, found once: the nesting below asks
    // them of every hole against every positive cycle, twice over.
    let boxed = |line: &[Point2]| -> (Point2, Point2) {
        line.iter().fold(
            (
                Point2::new(f64::INFINITY, f64::INFINITY),
                Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            ),
            |(lo, hi), p| {
                (
                    Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                    Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
                )
            },
        )
    };
    let within = |b: &(Point2, Point2), p: Point2| {
        p.x >= b.0.x && p.x <= b.1.x && p.y >= b.0.y && p.y <= b.1.y
    };
    let positive_nodes: Vec<_> = positives.iter().map(|(c, _)| nodes_of(c)).collect();
    let positive_area: Vec<f64> = positives.iter().map(|(_, l)| area(l).abs()).collect();
    let positive_box: Vec<_> = positives.iter().map(|(_, l)| boxed(l)).collect();
    // A lone ring meets no other cycle.
    let negative_nodes: Vec<_> = negatives
        .iter()
        .map(|(_, cycle, _)| {
            cycle.map_or_else(|_| std::collections::BTreeSet::new(), |c| nodes_of(c))
        })
        .collect();
    let meet = |x: &std::collections::BTreeSet<usize>, y: &std::collections::BTreeSet<usize>| {
        x.iter().any(|n| y.contains(n))
    };
    // A point outside a polygon's box is outside the polygon.
    let contains = |k: usize, p: Point2| within(&positive_box[k], p) && inside(&positives[k].1, p);
    // A hole belongs to the smallest positive cycle strictly containing
    // it. Sharing a node means same component, not a hole. Each hole's
    // containers are found once: a face with hundreds of holes has as many
    // positive cycles (each hole's own disc), and asking every pair about
    // every other cycle is the square of that again.
    //
    // Whether a hole lies in a cycle is asked at three of its chords'
    // midpoints, by majority. A hole may touch its cycle at isolated points
    // that are no node of either (a circle inscribed in a square touches it
    // at four), and a single point of the hole asked there lies on the
    // cycle and reads either way.
    let hole_inside = |k: usize, hole: &[Point2]| {
        let n = hole.len();
        let votes = [0, n / 3, 2 * n / 3]
            .into_iter()
            .filter(|&i| {
                let (a, b) = (hole[i % n], hole[(i + 1) % n]);
                contains(
                    k,
                    Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y)),
                )
            })
            .count();
        votes >= 2
    };
    // A cycle holding two of the three points holds the first or the
    // second, so only the cycles whose boxes hold one of those two are
    // asked, in their order.
    let grid = BoxGrid::new(&positive_box);
    let containers: Vec<Vec<usize>> = negatives
        .iter()
        .enumerate()
        .map(|(hi, (_, _, hole))| {
            let n = hole.len();
            let vote = |i: usize| {
                let (a, b) = (hole[i % n], hole[(i + 1) % n]);
                Point2::new(f64::midpoint(a.x, b.x), f64::midpoint(a.y, b.y))
            };
            grid.holding_either(vote(0), vote(n / 3))
                .into_iter()
                .filter(|&oi| {
                    hole_inside(oi, hole) && !meet(&negative_nodes[hi], &positive_nodes[oi])
                })
                .collect()
        })
        .collect();
    let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); positives.len()];
    for (hi, held) in containers.iter().enumerate() {
        for &pi in held {
            let direct = !held.iter().any(|&oi| {
                oi != pi
                    && contains(pi, positives[oi].1[0])
                    && positive_area[oi] < positive_area[pi]
            });
            if direct {
                holes_of[pi].push(hi);
            }
        }
    }
    let mut pieces = Vec::new();
    for (pi, (cycle, line)) in positives.iter().enumerate() {
        let mut rings = vec![traversals(cycle, live)];
        let mut rings_outline = vec![line.clone()];
        for &hi in &holes_of[pi] {
            let (_, hole_cycle, hole) = &negatives[hi];
            rings.push(match hole_cycle {
                Ok(cycle) => traversals(cycle, live),
                Err(ring) => ring.rings.clone(),
            });
            rings_outline.push(hole.clone());
        }
        // The roomiest probe decides whether the cycle bounds material: a
        // hole's own disc does not, and on a face with hundreds of holes
        // most cycles are such discs.
        let mut scanlines = Scanlines::new(&rings_outline, snap);
        let Some(interior) = scanlines.first_point(snap) else {
            continue;
        };
        if !material.inside(interior) {
            continue;
        }
        pieces.push(Piece {
            rings,
            outlines: rings_outline,
            first: interior,
            probes: Probes {
                scanlines,
                material: std::rc::Rc::clone(&material),
                snap,
                found: None,
            },
        });
    }
    if pieces.is_empty() {
        ogeom_bail!(Construction, "arrangement left no piece of the face");
    }
    Ok(pieces)
}

/// A hole side of a face: where it stands among the face's cycles, its
/// traversals (a cycle of the face's walk, or a lone ring's own), and its
/// outline.
type HoleSide<'a, T> = (
    (usize, bool),
    Result<&'a Vec<usize>, &'a Lone<T>>,
    Vec<Point2>,
);

/// A face's strands walked into cycles: the strands that take part, the
/// place of each among the strands given, the nodes each runs between, and
/// every cycle as darts (dart `2k` runs live strand `k` forward, `2k + 1`
/// backward), each starting at its least dart.
struct Walk<'s, T> {
    live: Vec<&'s Strand<T>>,
    from: Vec<usize>,
    ends: Vec<(usize, usize)>,
    cycles: Vec<Vec<usize>>,
}

/// Walk `strands` into cycles.
///
/// # Errors
///
/// As [`assemble`].
fn walk<T>(strands: &[Strand<T>], snap: f64) -> OgeomResult<Walk<'_, T>> {
    let (mut from, mut live): (Vec<usize>, Vec<&Strand<T>>) = strands
        .iter()
        .enumerate()
        .filter(|(_, s)| s.polyline.len() >= 2 && polyline_length(&s.polyline) > snap)
        .unzip();
    if live.is_empty() {
        ogeom_bail!(Construction, "a face with no boundary bounds nothing");
    }

    // Endpoints snapped to canonical nodes. Only endpoints: the pre-split
    // contract says strands meet nowhere else.
    let mut nodes: Vec<Point2> = Vec::new();
    // Nodes binned on a grid of the snap's own size, so an endpoint looks
    // only at the nine bins about it. Of the nodes in reach it takes the
    // first made, the one a scan in order would have met first.
    let cell = |p: Point2| -> (i64, i64) {
        let at = |x: f64| -> i64 {
            #[allow(clippy::cast_possible_truncation)]
            let k = (x / snap.max(f64::MIN_POSITIVE)).floor() as i64;
            k
        };
        (at(p.x), at(p.y))
    };
    let mut bins: ogeom_core::FastMap<(i64, i64), Vec<usize>> = ogeom_core::FastMap::default();
    let mut canon = |p: Point2, nodes: &mut Vec<Point2>| -> usize {
        let (cx, cy) = cell(p);
        let mut first: Option<usize> = None;
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                let key = (cx.saturating_add(dx), cy.saturating_add(dy));
                for &i in bins.get(&key).map_or(&[][..], Vec::as_slice) {
                    if nodes[i].distance(p) <= snap && first.is_none_or(|f| i < f) {
                        first = Some(i);
                    }
                }
            }
        }
        if let Some(i) = first {
            return i;
        }
        nodes.push(p);
        bins.entry((cx, cy)).or_default().push(nodes.len() - 1);
        nodes.len() - 1
    };
    let mut ends: Vec<(usize, usize)> = Vec::new();
    for strand in &live {
        let from = canon(strand.polyline[0], &mut nodes);
        let to = canon(
            *strand.polyline.last().unwrap_or(&strand.polyline[0]),
            &mut nodes,
        );
        ends.push((from, to));
    }

    // Dust welds through. A strand shorter than the snap is dropped from
    // the graph, but its two ends were still one junction of the boundary:
    // a run of consecutive dust pieces (the paving's crossing clusters on
    // a tolerant rail leave them) can jointly span more than the snap, and
    // dropping the pieces one by one would tear a gap no positional weld
    // reaches. Each dropped strand therefore aliases its endpoints, and a
    // chain of dust aliases end to end, so the surviving neighbours meet at
    // one node however long the run.
    {
        let mut parent: Vec<usize> = (0..nodes.len()).collect();
        let find = |parent: &mut Vec<usize>, mut i: usize| -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        };
        let mut any = false;
        for strand in strands {
            if strand.polyline.len() < 2 || polyline_length(&strand.polyline) > snap {
                continue;
            }
            let a = canon(strand.polyline[0], &mut nodes);
            let b = canon(
                *strand.polyline.last().unwrap_or(&strand.polyline[0]),
                &mut nodes,
            );
            while parent.len() < nodes.len() {
                parent.push(parent.len());
            }
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[rb] = ra;
                any = true;
            }
        }
        if any {
            for (from, to) in &mut ends {
                *from = find(&mut parent, *from);
                *to = find(&mut parent, *to);
            }
        }
    }

    // Prune dangling chains. A section that fails to separate material hangs
    // by an end. A *boundary* strand doing so means the face's own boundary
    // does not close, which no amount of pruning repairs.
    //
    // Peeled by a queue: a strand goes when either end is left on fewer
    // than two strands, and its going may leave its neighbours so, which
    // the queue takes up in turn. What goes is what peeling round after
    // round would take, and the survivors keep their order.
    {
        let mut degree = vec![0_usize; nodes.len()];
        let mut at_node: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
        for (k, (u, v)) in ends.iter().enumerate() {
            degree[*u] += 1;
            degree[*v] += 1;
            at_node[*u].push(k);
            at_node[*v].push(k);
        }
        let mut gone = vec![false; ends.len()];
        let mut queue: std::collections::VecDeque<usize> =
            (0..nodes.len()).filter(|n| degree[*n] < 2).collect();
        while let Some(n) = queue.pop_front() {
            for &k in &at_node[n] {
                if gone[k] {
                    continue;
                }
                let (u, v) = ends[k];
                if degree[u] >= 2 && degree[v] >= 2 {
                    continue;
                }
                let strand = live[k];
                if strand.boundary {
                    if std::env::var("OGEOM_ARRANGE_DEBUG").is_ok() {
                        eprintln!(
                            "DANGLE: boundary strand {:?} .. {:?} (deg {} / {})",
                            strand.polyline[0],
                            strand.polyline[strand.polyline.len() - 1],
                            degree[u],
                            degree[v]
                        );
                    }
                    ogeom_bail!(
                        Construction,
                        "a face boundary strand dangles; the boundary does not \
                         close in parameter space"
                    );
                }
                gone[k] = true;
                for end in [u, v] {
                    degree[end] -= 1;
                    if degree[end] < 2 {
                        queue.push_back(end);
                    }
                }
            }
        }
        let mut keep_live = Vec::with_capacity(live.len());
        let mut keep_from = Vec::with_capacity(live.len());
        let mut keep_ends = Vec::with_capacity(ends.len());
        for (k, ((strand, end), at)) in live.iter().zip(&ends).zip(&from).enumerate() {
            if !gone[k] {
                keep_live.push(*strand);
                keep_from.push(*at);
                keep_ends.push(*end);
            }
        }
        live = keep_live;
        from = keep_from;
        ends = keep_ends;
    }
    if live.is_empty() {
        ogeom_bail!(Construction, "the face's boundary vanished in arrangement");
    }

    // Darts: twins adjacent by construction: dart 2k runs strand k forward,
    // dart 2k + 1 backward.
    let dart_count = live.len() * 2;
    let head = |d: usize| -> usize {
        let (u, v) = ends[d / 2];
        if d.is_multiple_of(2) { v } else { u }
    };
    let tail = |d: usize| -> usize {
        let (u, v) = ends[d / 2];
        if d.is_multiple_of(2) { u } else { v }
    };
    // A dart's polyline point `i` steps from its tail, that way round.
    let point = |d: usize, i: usize| -> Point2 {
        let line = &live[d / 2].polyline;
        if d.is_multiple_of(2) {
            line[i]
        } else {
            line[line.len() - 1 - i]
        }
    };
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for d in 0..dart_count {
        around[tail(d)].push(d);
    }
    // The order of the darts round a node is read where each first leaves
    // a small circle about it, one circle for all of them. Strands that
    // meet only at their ends leave any such circle in the order they
    // stand round the node, where their first steps need not: a section
    // leaving a boundary strand all but tangentially (a side swept from a
    // face's outline crossing a round that meets that outline at a
    // tangent) can set out a hair to the wrong side of it before turning
    // in, and read by its first step it would cross the boundary. The
    // circle is a few snaps wide, and within half of every dart's reach
    // from the node, so each dart leaves it.
    let reach = |d: usize| -> f64 {
        let at = nodes[tail(d)];
        live[d / 2]
            .polyline
            .iter()
            .map(|p| p.distance(at))
            .fold(0.0_f64, f64::max)
    };
    let mut radius = vec![snap * 4.0; nodes.len()];
    for d in 0..dart_count {
        radius[tail(d)] = radius[tail(d)].min(reach(d) * 0.5);
    }
    let leaving = |d: usize| -> Point2 {
        let at = nodes[tail(d)];
        let r = radius[tail(d)];
        let count = live[d / 2].polyline.len();
        let mut previous = point(d, 0);
        for i in 1..count {
            let p = point(d, i);
            let (dp, dq) = (previous.distance(at), p.distance(at));
            if dq >= r {
                // Where the step crosses the circle, along the step.
                let t = if dq - dp > f64::MIN_POSITIVE {
                    ((r - dp) / (dq - dp)).clamp(0.0, 1.0)
                } else {
                    1.0
                };
                return Point2::new(
                    previous.x + (p.x - previous.x) * t,
                    previous.y + (p.y - previous.y) * t,
                );
            }
            previous = p;
        }
        point(d, count - 1)
    };
    // Each dart's angle once, not once per comparison.
    let heading: Vec<f64> = (0..dart_count)
        .map(|d| angle(nodes[tail(d)], leaving(d)))
        .collect();
    for ring in &mut around {
        ring.sort_by(|&x, &y| {
            heading[x]
                .partial_cmp(&heading[y])
                .unwrap_or(core::cmp::Ordering::Equal)
        });
    }
    let mut position = vec![0_usize; dart_count];
    for ring in &around {
        for (at, &d) in ring.iter().enumerate() {
            position[d] = at;
        }
    }
    let next = |d: usize| -> usize {
        let ring = &around[head(d)];
        let at = position[d ^ 1];
        // The dart before the twin in counter-clockwise order: the face on
        // the left continues there.
        ring[(at + ring.len() - 1) % ring.len()]
    };

    let mut seen = vec![false; dart_count];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for start in 0..dart_count {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut d = start;
        loop {
            seen[d] = true;
            cycle.push(d);
            d = next(d);
            if d == start {
                break;
            }
        }
        cycles.push(cycle);
    }

    Ok(Walk {
        live,
        from,
        ends,
        cycles,
    })
}

/// The polyline a cycle traces, for area, containment and interior points.
fn outline_of<T>(live: &[&Strand<T>], cycle: &[usize]) -> Vec<Point2> {
    let mut out: Vec<Point2> = Vec::new();
    for &d in cycle {
        let line = &live[d / 2].polyline;
        let mut points: Vec<Point2> = if d.is_multiple_of(2) {
            line.clone()
        } else {
            line.iter().rev().copied().collect()
        };
        points.pop();
        out.append(&mut points);
    }
    out
}

/// A hole's ring whose strands stand apart from every other strand of its
/// face: no other strand ends within a snap of them or runs inside them.
/// Walked alone it is the cycle the face's whole walk would find, so the
/// face's walk leaves it out and takes it back as it is.
pub(crate) struct Lone<T> {
    /// Where the ring stands among the face's cycles: the place of its
    /// first strand among the face's strands, and whether it runs that
    /// strand backward.
    order: (usize, bool),
    rings: Vec<Traversal<T>>,
    outline: Vec<Point2>,
    /// The ring's strands, which bound the face's material.
    lines: Vec<Vec<Point2>>,
}

/// The hole side of the one ring `strands` close, `places` their places
/// among their face's strands; `None` where they do not close exactly one
/// ring, every strand taking part, with a hole side.
pub(crate) fn lone_ring<T: Clone>(
    strands: &[Strand<T>],
    places: &[usize],
    snap: f64,
) -> Option<Lone<T>> {
    let walk = walk(strands, snap).ok()?;
    if walk.live.len() != strands.len() || walk.cycles.len() != 2 {
        return None;
    }
    let mut hole = None;
    let mut disc = false;
    for cycle in &walk.cycles {
        let line = outline_of(&walk.live, cycle);
        let a = area(&line);
        if a > snap * snap {
            disc = true;
        } else if a < -(snap * snap) {
            hole = Some((cycle, line));
        }
    }
    let (cycle, outline) = hole.filter(|_| disc)?;
    Some(Lone {
        order: (places[walk.from[cycle[0] / 2]], cycle[0] % 2 == 1),
        rings: traversals(cycle, &walk.live),
        outline,
        lines: walk.live.iter().map(|s| s.polyline.clone()).collect(),
    })
}

fn traversals<T: Clone>(cycle: &[usize], live: &[&Strand<T>]) -> Vec<Traversal<T>> {
    cycle
        .iter()
        .map(|&d| Traversal {
            tag: live[d / 2].tag.clone(),
            reversed: d % 2 == 1,
        })
        .collect()
}

fn polyline_length(line: &[Point2]) -> f64 {
    line.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The angle of the direction from `from` towards `to`.
fn angle(from: Point2, to: Point2) -> f64 {
    let v = to - from;
    v.y.atan2(v.x)
}

/// The signed area: positive for counter-clockwise.
fn area(ring: &[Point2]) -> f64 {
    let mut doubled = 0.0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        doubled += a.x.mul_add(b.y, -(b.x * a.y));
    }
    doubled * 0.5
}

/// Whether the ray from `p` toward `+x` crosses the segment `a`→`b`,
/// decided by the predicate rather than by a computed intersection.
///
/// The straddle test picks the segments the ray could cross. For those, the
/// crossing question is exactly "which side of the segment's line does `p`
/// lie on", which is `orient2d`'s question. A point exactly on the segment's
/// line reads as no crossing. The callers probe interior sample points,
/// never boundary ones.
fn ray_crosses_segment<P: Predicates>(a: Point2, b: Point2, p: Point2) -> bool {
    if (a.y > p.y) == (b.y > p.y) {
        return false;
    }
    let side = P::orient2d([a.x, a.y], [b.x, b.y], [p.x, p.y]);
    if b.y > a.y {
        side == Sign::Positive
    } else {
        side == Sign::Negative
    }
}

/// The leaning ray's slope. A tangent chain puts whole strands exactly
/// along an axis-aligned junction line, where a horizontal ray grazes
/// corner after corner inside rounding noise and counts them at random. No
/// real boundary runs along this slope.
const SLANT: f64 = 0.618_033_988_749_894_9;

/// As [`ray_crosses_segment`], with the ray leaning off the axes.
fn slanted_ray_crosses_segment<P: Predicates>(a: Point2, b: Point2, p: Point2) -> bool {
    let shifted = [p.x + 1.0, p.y + SLANT];
    let above = |q: Point2| P::orient2d([p.x, p.y], shifted, [q.x, q.y]) == Sign::Positive;
    if above(a) == above(b) {
        return false;
    }
    let side = P::orient2d([a.x, a.y], [b.x, b.y], [p.x, p.y]);
    if above(b) {
        side == Sign::Positive
    } else {
        side == Sign::Negative
    }
}

/// Whether the horizontal ray from `p` crosses the segment from `a` to `b`,
/// as [`inside_rings`] counts a crossing.
pub(crate) fn ray_crosses(a: Point2, b: Point2, p: Point2) -> bool {
    ray_crosses_segment::<Exact>(a, b, p)
}

/// Even-odd containment of a point in one closed polyline.
fn inside(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if ray_crosses_segment::<Exact>(a, b, p) {
            inside = !inside;
        }
    }
    inside
}

/// Even-odd containment of a point in a piece's rings.
///
/// A ring is a closed loop stored without repeating its first point, so each
/// one is closed as it is walked. Ring `[0]` is the outer contour and the
/// rest are holes, and the even-odd count over all of them answers for the
/// region they jointly bound.
///
/// This is the test for a *piece*. [`inside_many`] is the test for the
/// strands a piece is assembled from, which close only jointly and must not
/// be closed one by one.
pub(crate) fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    rings.iter().fold(false, |acc, ring| acc != inside(ring, p))
}

/// Open polylines that jointly close, held with their boxes, for asking
/// [`inside_many`] and its leaning twin of many points: each ray
/// passes over the polylines whose box it cannot meet.
pub(crate) struct Trim {
    lines: Vec<Vec<Point2>>,
    /// `(low x, high x, low y, high y)` of each line.
    boxes: Vec<(f64, f64, f64, f64)>,
}

impl Trim {
    pub(crate) fn new(lines: Vec<Vec<Point2>>) -> Self {
        let boxes = lines
            .iter()
            .map(|line| {
                line.iter().fold(
                    (
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                    ),
                    |(lx, hx, ly, hy), q| (lx.min(q.x), hx.max(q.x), ly.min(q.y), hy.max(q.y)),
                )
            })
            .collect();
        Self { lines, boxes }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    /// The middle of the box round every line.
    pub(crate) fn middle(&self) -> Point2 {
        let (lx, hx, ly, hy) = self.boxes.iter().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(lx, hx, ly, hy), b| (lx.min(b.0), hx.max(b.1), ly.min(b.2), hy.max(b.3)),
        );
        Point2::new(f64::midpoint(lx, hx), f64::midpoint(ly, hy))
    }

    /// [`inside_many`]: a line wholly above or below the point, or wholly
    /// left of it, is not crossed.
    pub(crate) fn inside(&self, p: Point2) -> bool {
        let mut inside = false;
        for (line, &(_, hx, ly, hy)) in self.lines.iter().zip(&self.boxes) {
            if ly > p.y || hy <= p.y || hx < p.x {
                continue;
            }
            for w in line.windows(2) {
                if ray_crosses_segment::<Exact>(w[0], w[1], p) {
                    inside = !inside;
                }
            }
        }
        inside
    }

    /// [`inside_many`] along a leaning ray: a line whose box stands wholly on one side
    /// of the leaning ray's line, by the same exact test the crossing
    /// asks of each end, is not crossed.
    pub(crate) fn inside_slanted(&self, p: Point2) -> bool {
        let shifted = [p.x + 1.0, p.y + SLANT];
        let above = |x: f64, y: f64| Exact::orient2d([p.x, p.y], shifted, [x, y]) == Sign::Positive;
        let mut inside = false;
        for (line, &(lx, hx, ly, hy)) in self.lines.iter().zip(&self.boxes) {
            if lx.partial_cmp(&hx).is_none_or(core::cmp::Ordering::is_gt) {
                continue;
            }
            let corners = [above(lx, ly), above(lx, hy), above(hx, ly), above(hx, hy)];
            if corners.iter().all(|&c| c == corners[0]) {
                continue;
            }
            for w in line.windows(2) {
                if slanted_ray_crosses_segment::<Exact>(w[0], w[1], p) {
                    inside = !inside;
                }
            }
        }
        inside
    }
}

/// Boxes filed on a uniform grid over their joint extent, for finding the
/// boxes that may hold a point without asking every one. A box spanning
/// many cells is kept apart and offered for every point.
struct BoxGrid {
    low: Point2,
    cell: (f64, f64),
    side: usize,
    cells: Vec<Vec<usize>>,
    wide: Vec<usize>,
}

impl BoxGrid {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "a count of cells"
    )]
    fn new(boxes: &[(Point2, Point2)]) -> Self {
        let finite = |b: &(Point2, Point2)| {
            b.0.x.is_finite() && b.0.y.is_finite() && b.1.x.is_finite() && b.1.y.is_finite()
        };
        let (mut low, mut high) = (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for b in boxes.iter().filter(|b| finite(b)) {
            low = Point2::new(low.x.min(b.0.x), low.y.min(b.0.y));
            high = Point2::new(high.x.max(b.1.x), high.y.max(b.1.y));
        }
        let side = ((boxes.len() as f64).sqrt().ceil() as usize).clamp(1, 64);
        let span = |lo: f64, hi: f64| {
            if hi > lo {
                (hi - lo) / side as f64
            } else {
                1.0
            }
        };
        let mut grid = Self {
            low,
            cell: (span(low.x, high.x), span(low.y, high.y)),
            side,
            cells: vec![Vec::new(); side * side],
            wide: Vec::new(),
        };
        let most = (side * side / 16).max(16);
        for (k, b) in boxes.iter().enumerate() {
            if !finite(b) {
                continue;
            }
            let (x0, y0) = grid.cell_of(b.0);
            let (x1, y1) = grid.cell_of(b.1);
            if (x1 - x0 + 1) * (y1 - y0 + 1) > most {
                grid.wide.push(k);
                continue;
            }
            for y in y0..=y1 {
                for x in x0..=x1 {
                    grid.cells[y * side + x].push(k);
                }
            }
        }
        grid
    }

    /// The cell a point falls in, clamped to the grid.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a cell index"
    )]
    fn cell_of(&self, p: Point2) -> (usize, usize) {
        let at = |x: f64, low: f64, size: f64| -> usize {
            (((x - low) / size).floor().max(0.0) as usize).min(self.side - 1)
        };
        (
            at(p.x, self.low.x, self.cell.0),
            at(p.y, self.low.y, self.cell.1),
        )
    }

    /// Every box that may hold `p` or `q`, ascending: all that do, and
    /// some that do not.
    fn holding_either(&self, p: Point2, q: Point2) -> Vec<usize> {
        let mut out = self.wide.clone();
        for point in [p, q] {
            let (x, y) = self.cell_of(point);
            out.extend_from_slice(&self.cells[y * self.side + x]);
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}

/// Polylines with their boxes, for asking [`inside_many`] of many points.
struct Boxed<'a> {
    lines: Vec<&'a [Point2]>,
    /// `(low x, high x, low y, high y)` of each line.
    boxes: Vec<(f64, f64, f64, f64)>,
}

impl<'a> Boxed<'a> {
    fn new(lines: Vec<&'a [Point2]>) -> Self {
        let boxes = lines
            .iter()
            .map(|line| {
                line.iter().fold(
                    (
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                    ),
                    |(lx, hx, ly, hy), q| (lx.min(q.x), hx.max(q.x), ly.min(q.y), hy.max(q.y)),
                )
            })
            .collect();
        Self { lines, boxes }
    }

    /// [`inside_many`] along the leaning ray, skipping the polylines whose
    /// box stands wholly on one side of its line. A probe stands midway
    /// between its own cycle's vertex heights, which a symmetric face puts
    /// level with a corner of the boundary; two strands ending at that
    /// corner a rounding apart in height straddle a level ray both or
    /// neither. The leaning ray passes no corner.
    fn inside(&self, p: Point2) -> bool {
        let shifted = [p.x + 1.0, p.y + SLANT];
        let above = |x: f64, y: f64| Exact::orient2d([p.x, p.y], shifted, [x, y]) == Sign::Positive;
        let mut inside = false;
        for (line, &(lx, hx, ly, hy)) in self.lines.iter().zip(&self.boxes) {
            if lx.partial_cmp(&hx).is_none_or(core::cmp::Ordering::is_gt) {
                continue;
            }
            let corners = [above(lx, ly), above(lx, hy), above(hx, ly), above(hx, hy)];
            if corners.iter().all(|&c| c == corners[0]) {
                continue;
            }
            for w in line.windows(2) {
                if slanted_ray_crosses_segment::<Exact>(w[0], w[1], p) {
                    inside = !inside;
                }
            }
        }
        inside
    }
}

/// Even-odd containment against open polylines that jointly close.
///
/// The strands are pieces of closed rings, so counting crossings segment by
/// segment over all of them gives the same even-odd answer the assembled
/// rings would.
pub(crate) fn inside_many(lines: &[&[Point2]], p: Point2) -> bool {
    let mut inside = false;
    for line in lines {
        for w in line.windows(2) {
            if ray_crosses_segment::<Exact>(w[0], w[1], p) {
                inside = !inside;
            }
        }
    }
    inside
}

/// Points strictly inside the region the rings bound, by scanline.
///
/// A horizontal line through the widest gap between distinct vertex heights
/// cannot pass through a vertex or run along a horizontal segment, so its
/// crossings are transversal and the midpoint of the first inside interval
/// is interior with room to spare.
///
/// The rest are for the caller who has to ask again. They vary in *both*
/// chart directions (other heights, and other positions along each) because
/// a piece that merely touches the other solid touches it somewhere, and a
/// second opinion taken from the same place is not one.
/// [`interior_points`] for rings the caller assembled itself.
pub(crate) fn interior_points_of(rings: &[Vec<Point2>], snap: f64) -> Vec<Point2> {
    interior_points(rings, snap)
}

/// Points strictly inside the region the rings bound, placed at golden-ratio
/// fractions of the widest scanline gaps and of each inside interval.
///
/// [`interior_points`] offers halves and quarters, and a solid touching a
/// piece only along lines or at points can stand on every one of those: a
/// cylinder inscribed in a box touches its wall at the quarter turns, which
/// are the quarter columns of the wall's chart. A golden fraction is no
/// simple ratio of the piece's extent, so a contact that met every regular
/// probe is not met here unless it covers a region.
pub(crate) fn off_contact_points(rings: &[Vec<Point2>], snap: f64) -> Vec<Point2> {
    const FRACTIONS: [f64; 3] = [
        0.381_966_011_250_105,
        0.618_033_988_749_895,
        0.145_898_033_750_315,
    ];
    let mut heights: Vec<f64> = rings.iter().flatten().map(|p| p.y).collect();
    heights.sort_unstable_by(f64::total_cmp);
    heights.dedup_by(|a, b| (*a - *b).abs() <= snap);
    let mut gaps: Vec<(f64, f64)> = heights
        .windows(2)
        .map(|pair| (pair[1] - pair[0], pair[0]))
        .filter(|(gap, _)| *gap > snap)
        .collect();
    gaps.sort_by(|a, b| b.0.total_cmp(&a.0));
    gaps.truncate(3);
    let segments = Bands::new(rings);
    let mut points = Vec::new();
    for (gap, low) in gaps {
        for along in FRACTIONS {
            let level = gap.mul_add(along, low);
            let mut crossings: Vec<f64> = Vec::new();
            for &(a, b) in segments.near(level) {
                if (a.y > level) != (b.y > level) {
                    crossings.push((b.x - a.x).mul_add((level - a.y) / (b.y - a.y), a.x));
                }
            }
            crossings.sort_by(f64::total_cmp);
            for pair in crossings.as_chunks::<2>().0 {
                let width = pair[1] - pair[0];
                if width > snap {
                    for across in FRACTIONS {
                        points.push(Point2::new(width.mul_add(across, pair[0]), level));
                    }
                }
            }
        }
    }
    points
}

/// The scanlines a region's interior points are sought along, cast as
/// asked and kept: the widest gap between distinct vertex heights first
/// with its quarter heights straight after it, then the gaps narrowing,
/// equal gaps lowest first.
struct Scanlines {
    gaps: Vec<(f64, f64)>,
    widest: std::collections::BinaryHeap<(u64, std::cmp::Reverse<usize>)>,
    levels: Vec<(f64, f64)>,
    segments: Bands,
}

impl Scanlines {
    fn new(rings: &[Vec<Point2>], snap: f64) -> Self {
        let mut heights: Vec<f64> = rings.iter().flatten().map(|p| p.y).collect();
        heights.sort_unstable_by(f64::total_cmp);
        heights.dedup_by(|a, b| (*a - *b).abs() <= snap);
        // Drawn from a heap as the search asks, since it mostly stops after
        // a few dozen of thousands.
        let gaps: Vec<(f64, f64)> = heights
            .windows(2)
            .map(|pair| (pair[1] - pair[0], f64::midpoint(pair[0], pair[1])))
            .collect();
        // A gap is never negative, and a non-negative float's bits order as
        // the float does.
        let widest = gaps
            .iter()
            .enumerate()
            .map(|(i, &(gap, _))| {
                (
                    if gap > 0.0 { gap.to_bits() } else { 0 },
                    std::cmp::Reverse(i),
                )
            })
            .collect();
        Self {
            gaps,
            widest,
            levels: Vec::new(),
            segments: Bands::new(rings),
        }
    }

    /// The `index`th scanline's gap and height.
    fn level(&mut self, index: usize) -> Option<(f64, f64)> {
        while self.levels.len() <= index {
            let (_, std::cmp::Reverse(i)) = self.widest.pop()?;
            self.levels.push(self.gaps[i]);
            // An unsplit chart (a whole sphere, pole to pole) has exactly
            // one gap and therefore one scanline, straight across its
            // middle. That is precisely where a solid seated on its equator
            // touches it, so the widest gap also offers its quarter heights:
            // same piece, different latitude.
            if self.levels.len() == 1 {
                let (gap, level) = self.gaps[i];
                self.levels.push((gap, gap.mul_add(-0.25, level)));
                self.levels.push((gap, gap.mul_add(0.25, level)));
            }
        }
        self.levels.get(index).copied()
    }

    /// Where the scanline at `level` crosses the rings, in order.
    fn crossings(&self, level: f64) -> Vec<f64> {
        let mut crossings: Vec<f64> = Vec::new();
        for &(a, b) in self.segments.near(level) {
            if (a.y > level) != (b.y > level) {
                crossings.push((b.x - a.x).mul_add((level - a.y) / (b.y - a.y), a.x));
            }
        }
        crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        crossings
    }

    /// The first of [`Scanlines::points`] alone: the roomiest candidate,
    /// the earliest among equals, found by casting scanlines only until no
    /// later one can offer more room.
    fn first_point(&mut self, snap: f64) -> Option<Point2> {
        let mut best: Option<(f64, Point2)> = None;
        for index in 0.. {
            let Some((gap, level)) = self.level(index) else {
                break;
            };
            // A candidate's room is never more than its scanline's gap, and
            // the gaps only narrow: from here on none is roomier.
            if best.is_some_and(|(room, _)| gap <= room) {
                break;
            }
            if gap <= snap {
                continue;
            }
            for pair in self.crossings(level).as_chunks::<2>().0 {
                let width = pair[1] - pair[0];
                if width > snap {
                    // A quarter position has half its midpoint's room, so
                    // only the midpoint can be the roomiest.
                    let room = width.min(gap);
                    if best.is_none_or(|(held, _)| room > held) {
                        best = Some((room, Point2::new(f64::midpoint(pair[0], pair[1]), level)));
                    }
                }
            }
        }
        best.map(|(_, p)| p)
    }

    /// Points strictly inside the region, best first.
    fn points(&mut self, snap: f64) -> Vec<Point2> {
        // Every inside interval of every scanline is a candidate, and they
        // are ranked by width. The first interval of the roomiest scanline
        // is not good enough: a piece with a cusp (the sliver beside a line
        // tangent to a circle, which is what a ball inscribed in a cylinder
        // leaves on any plane through both) puts that interval inside the
        // cusp, where the "interior" point is within rounding of the
        // boundary and reads as lying on it.
        let mut candidates: Vec<(f64, Point2)> = Vec::new();
        // A candidate's room is never more than its scanline's gap, and the
        // scanlines come widest first. Once the candidates roomier than the
        // next scanline's gap already settle every probe, no later scanline
        // can change the choice, and the rest (on a face with hundreds of
        // holes, thousands of scanlines each crossing every segment) are not
        // cast.
        let mut checked_at = 1;
        for index in 0.. {
            let Some((gap, level)) = self.level(index) else {
                break;
            };
            if index == checked_at {
                checked_at *= 2;
                let mut settled: Vec<(f64, Point2)> = candidates
                    .iter()
                    .copied()
                    .filter(|(room, _)| *room > gap)
                    .collect();
                settled.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(core::cmp::Ordering::Equal));
                if choose(&settled, snap).1 {
                    break;
                }
            }
            if gap <= snap {
                continue;
            }
            for pair in self.crossings(level).as_chunks::<2>().0 {
                let width = pair[1] - pair[0];
                if width > snap {
                    // Room is the lesser of the interval's width and its
                    // scanline's gap: a level in a thin gap may pass just
                    // beyond an arc's chords and inside the arc itself, where
                    // the chords report an interval as wide as the whole
                    // piece.
                    let room = width.min(gap);
                    candidates.push((room, Point2::new(f64::midpoint(pair[0], pair[1]), level)));
                    // And its quarter positions, ranked below the midpoint.
                    // Moving the scanline is not enough on its own: a piece
                    // symmetric about a chart-vertical line (a cylinder band,
                    // a revolved wall, a chart rectangle) has the same
                    // midpoint at every height, so a solid touching it along
                    // that line is met by every one of these "different"
                    // probes at once. The quarter heights above exist for
                    // the same reason in the other direction. This is that
                    // rule, applied to the width.
                    candidates.push((room * 0.5, Point2::new(width.mul_add(0.25, pair[0]), level)));
                    candidates.push((room * 0.5, Point2::new(width.mul_add(0.75, pair[0]), level)));
                }
            }
        }

        // Room first, but never room alone. Ranked purely by width, every
        // scanline of a piece symmetric about a chart-vertical line offers
        // its own midpoint before any of them offers a different position,
        // so a caller asking nine times asks the same column nine times.
        // Take the roomiest candidate at each distinct column first, then the
        // rest in width order, which leaves the other columns at the front
        // for a touch running down one, and the other heights right behind
        // them for a touch running across.
        //
        // Nine are returned, so nine columns end the search, and nine of the
        // rest are all that can follow them: a face with hundreds of holes
        // offers hundreds of thousands of candidates, and comparing each with
        // every column already chosen costs seconds for nine points.
        candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(core::cmp::Ordering::Equal));
        choose(&candidates, snap).0
    }
}

fn interior_points(rings: &[Vec<Point2>], snap: f64) -> Vec<Point2> {
    Scanlines::new(rings, snap).points(snap)
}

/// The rings' segments filed by height in bands, so a scanline is met
/// only by the segments of its own band: a face with hundreds of holes
/// casts thousands of scanlines.
struct Bands {
    low: f64,
    height: f64,
    bands: Vec<Vec<(Point2, Point2)>>,
}

impl Bands {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "a count of bands"
    )]
    fn new(rings: &[Vec<Point2>]) -> Self {
        let segments: Vec<(Point2, Point2)> = rings
            .iter()
            .flat_map(|ring| (0..ring.len()).map(move |i| (ring[i], ring[(i + 1) % ring.len()])))
            .collect();
        let (low, high) = segments
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), (a, b)| {
                (lo.min(a.y).min(b.y), hi.max(a.y).max(b.y))
            });
        let count = ((segments.len() as f64).sqrt() as usize).max(1);
        let height = if high > low {
            (high - low) / count as f64
        } else {
            1.0
        };
        let mut out = Self {
            low,
            height,
            bands: Vec::new(),
        };
        // Each band's share counted first, so every band is filled once
        // into room already its own.
        let spans: Vec<(usize, usize)> = segments
            .iter()
            .map(|(a, b)| {
                let (lo, hi) = if a.y <= b.y { (a.y, b.y) } else { (b.y, a.y) };
                (out.band_of(lo, count), out.band_of(hi, count))
            })
            .collect();
        let mut sizes = vec![0_usize; count];
        for &(from, to) in &spans {
            for size in &mut sizes[from..=to] {
                *size += 1;
            }
        }
        out.bands = sizes.into_iter().map(Vec::with_capacity).collect();
        for (segment, (from, to)) in segments.into_iter().zip(spans) {
            for band in &mut out.bands[from..=to] {
                band.push(segment);
            }
        }
        out
    }

    /// The band a height falls in, clamped to the bands there are.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a band index"
    )]
    fn band(&self, y: f64) -> usize {
        self.band_of(y, self.bands.len())
    }

    /// [`Bands::band`] among `count` bands.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a band index"
    )]
    fn band_of(&self, y: f64, count: usize) -> usize {
        // Truncation is the floor of a height at or above the lowest.
        let at = ((y - self.low) / self.height).max(0.0) as usize;
        at.min(count - 1)
    }

    /// The segments that may straddle `level`, in the rings' order.
    fn near(&self, level: f64) -> &[(Point2, Point2)] {
        &self.bands[self.band(level)]
    }
}

/// How many interior points a piece offers.
const PROBES: usize = 9;

/// The probes [`interior_points`] offers from its ranked candidates: the
/// roomiest at each distinct column first, then the rest in width order;
/// and whether the columns alone filled them, so that no candidate ranked
/// lower could change the choice.
fn choose(candidates: &[(f64, Point2)], snap: f64) -> (Vec<Point2>, bool) {
    let mut chosen: Vec<(f64, Point2)> = Vec::new();
    let mut rest: Vec<(f64, Point2)> = Vec::new();
    for &candidate in candidates {
        if chosen.len() == PROBES {
            break;
        }
        if chosen
            .iter()
            .any(|(_, p)| (p.x - candidate.1.x).abs() <= snap)
        {
            if rest.len() < PROBES {
                rest.push(candidate);
            }
        } else {
            chosen.push(candidate);
        }
    }
    let full = chosen.len() == PROBES;
    chosen.extend(rest);
    chosen.truncate(PROBES);
    (chosen.into_iter().map(|(_, p)| p).collect(), full)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Every strand walked together, in its own place.
    fn assemble<T: Clone>(strands: &[Strand<T>], snap: f64) -> OgeomResult<Vec<Piece<'_, T>>> {
        let places: Vec<usize> = (0..strands.len()).collect();
        super::assemble(strands, &places, snap, &[])
    }

    /// The even-odd containment is a sign decision, and it goes through the
    /// predicate: the crossing question near a slanted edge is answered by
    /// `orient2d` on both sides of the edge, at separations where a computed
    /// intersection abscissa would be pure rounding.
    #[test]
    fn containment_is_decided_by_the_predicate() {
        let ring = vec![
            Point2::new(0.0, -1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(-1.0, 0.0),
        ];
        assert!(inside(&ring, Point2::new(0.0, 0.0)));
        assert!(!inside(&ring, Point2::new(0.75, 0.5)));
        // Astride the slanted edge x + y = 1: the side flips with the sign.
        let eps = 1e-14;
        assert!(inside(&ring, Point2::new(0.5 - eps, 0.5)));
        assert!(!inside(&ring, Point2::new(0.5 + eps, 0.5)));
        // A point exactly on the edge's line, straddle satisfied, reads as
        // no crossing from either side: `Sign::Zero` is no crossing.
        assert!(!ray_crosses_segment::<Exact>(
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(0.5, 0.5),
        ));
    }

    /// A square's four sides as boundary strands, pre-split at the corners.
    fn square(side: f64) -> Vec<Strand<usize>> {
        let c = [
            Point2::new(0.0, 0.0),
            Point2::new(side, 0.0),
            Point2::new(side, side),
            Point2::new(0.0, side),
        ];
        (0..4)
            .map(|i| Strand {
                polyline: vec![c[i], c[(i + 1) % 4]],
                tag: i,
                boundary: true,
            })
            .collect()
    }

    fn strand(from: Point2, to: Point2, tag: usize) -> Strand<usize> {
        Strand {
            polyline: vec![from, to],
            tag,
            boundary: false,
        }
    }

    #[test]
    fn an_uncut_face_is_one_piece_with_its_tags_in_order() {
        let square = square(2.0);
        let pieces = assemble(&square, 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].rings.len(), 1);
        assert_eq!(pieces[0].rings[0].len(), 4);
    }

    #[test]
    fn a_rectangle_offers_probes_away_from_its_own_centre_line() {
        // A piece symmetric about a chart-vertical line has the same midpoint
        // at every height, so ranking probes by room alone offers one column
        // over and over. A solid touching this piece down that column would
        // meet every probe at once, and the classifier (which asks again
        // precisely because the first answer was "on the boundary") would
        // have nowhere else to ask.
        let square = square(2.0);
        let mut pieces = assemble(&square, 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);

        let interiors = pieces[0].interiors().to_vec();
        let columns = interiors.iter().map(|p| p.x);
        let spread = columns.clone().fold(f64::NEG_INFINITY, f64::max)
            - columns.fold(f64::INFINITY, f64::min);
        assert!(
            spread > 1e-6,
            "every probe stands in the same column: {interiors:?}"
        );
    }

    /// A section setting out all but along the boundary, its first step a
    /// hair to the outside before it turns in (a section starting at a
    /// tangent, sampled): read where the strands leave a circle about their
    /// node, it stands inside, and the face splits in two.
    #[test]
    fn a_section_leaving_along_the_boundary_still_splits_the_face() {
        let p = Point2::new;
        let boundary = |a: Point2, b: Point2, tag: usize| Strand {
            polyline: vec![a, b],
            tag,
            boundary: true,
        };
        let strands = vec![
            boundary(p(0.0, 0.0), p(0.5, 0.0), 0),
            boundary(p(0.5, 0.0), p(1.0, 0.0), 1),
            boundary(p(1.0, 0.0), p(1.0, 1.0), 2),
            boundary(p(1.0, 1.0), p(0.5, 1.0), 3),
            boundary(p(0.5, 1.0), p(0.0, 1.0), 4),
            boundary(p(0.0, 1.0), p(0.0, 0.0), 5),
            Strand {
                polyline: vec![
                    p(0.5, 0.0),
                    p(0.501, -1e-7),
                    p(0.51, 0.001),
                    p(0.55, 0.05),
                    p(0.6, 0.5),
                    p(0.5, 1.0),
                ],
                tag: 6,
                boundary: false,
            },
        ];
        let pieces = assemble(&strands, 1e-3).unwrap();
        assert_eq!(pieces.len(), 2);
    }

    /// A parallelogram whose left and right corners stand level at height
    /// 0, cut by two chords into three strips. The middle strip's probe
    /// stands midway between its corners' heights, level with the right
    /// corner, where the two boundary strands ending there disagree in
    /// height by a rounding: the strip is material all the same, and the
    /// three strips make the whole parallelogram.
    #[test]
    fn a_strip_whose_probe_is_level_with_a_corner_is_kept() {
        let p = Point2::new;
        let boundary = |a: Point2, b: Point2, tag: usize| Strand {
            polyline: vec![a, b],
            tag,
            boundary: true,
        };
        let strands = vec![
            boundary(p(0.0, 0.0), p(20.0, -5.0), 0),
            boundary(p(20.0, -5.0), p(20.4375, -3.25), 1),
            boundary(p(20.4375, -3.25), p(20.875, -1.5), 2),
            boundary(p(20.875, -1.5), p(21.25, 0.0), 3),
            boundary(p(21.25, 4.4e-16), p(1.25, 5.0), 4),
            boundary(p(1.25, 5.0), p(0.8125, 3.25), 5),
            boundary(p(0.8125, 3.25), p(0.375, 1.5), 6),
            boundary(p(0.375, 1.5), p(0.0, 0.0), 7),
            strand(p(0.8125, 3.25), p(20.875, -1.5), 8),
            strand(p(20.4375, -3.25), p(0.375, 1.5), 9),
        ];
        let pieces = assemble(&strands, 1e-6).unwrap();
        let areas: Vec<f64> = pieces.iter().map(|q| area(&q.outlines[0]).abs()).collect();
        assert_eq!(pieces.len(), 3, "{areas:?}");
        let whole: f64 = areas.iter().sum();
        assert!((whole - 106.25).abs() < 1e-9, "{areas:?}");
    }

    #[test]
    fn a_chord_split_at_the_boundary_makes_two_pieces() {
        // The caller's pre-split contract: the chord arrives as one strand
        // spanning wall to wall, and the walls arrive split at its feet.
        let mut strands = vec![
            Strand {
                polyline: vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)],
                tag: 0,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 0.0), Point2::new(2.0, 1.0)],
                tag: 1,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 1.0), Point2::new(2.0, 2.0)],
                tag: 2,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 2.0), Point2::new(0.0, 2.0)],
                tag: 3,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 2.0), Point2::new(0.0, 1.0)],
                tag: 4,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 1.0), Point2::new(0.0, 0.0)],
                tag: 5,
                boundary: true,
            },
        ];
        strands.push(strand(Point2::new(0.0, 1.0), Point2::new(2.0, 1.0), 6));
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        // The chord appears in both pieces, once each way round.
        let uses: Vec<bool> = pieces
            .iter()
            .flat_map(|p| p.rings[0].iter())
            .filter(|t| t.tag == 6)
            .map(|t| t.reversed)
            .collect();
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0], uses[1]);
    }

    #[test]
    fn a_dangling_section_is_pruned_not_walked() {
        let mut strands = square(2.0);
        strands.push(strand(Point2::new(1.0, 0.0), Point2::new(1.0, 1.0), 9));
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert!(pieces[0].rings[0].iter().all(|t| t.tag != 9));
    }

    #[test]
    fn a_dangling_boundary_is_an_error() {
        let mut strands = square(2.0);
        strands.pop();
        assert!(assemble(&strands, 1e-7).is_err());
    }

    #[test]
    fn a_closed_section_loop_gives_a_hole_with_its_tags() {
        // The loop arrives as two arcs (the closed-curve pre-split), wholly
        // inside the face: the inner region is a piece, and the rest is a
        // piece with a hole whose ring carries the arcs' tags.
        let mut strands = square(4.0);
        let top = Point2::new(2.0, 3.0);
        let bottom = Point2::new(2.0, 1.0);
        strands.push(Strand {
            polyline: vec![bottom, Point2::new(3.0, 2.0), top],
            tag: 10,
            boundary: false,
        });
        strands.push(Strand {
            polyline: vec![top, Point2::new(1.0, 2.0), bottom],
            tag: 11,
            boundary: false,
        });
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        let with_hole = pieces.iter().find(|p| p.rings.len() == 2).unwrap();
        let tags: Vec<usize> = with_hole.rings[1].iter().map(|t| t.tag).collect();
        assert!(tags.contains(&10) && tags.contains(&11));
        let inner = pieces
            .iter()
            .find(|p| p.rings.len() == 1 && p.rings[0].len() == 2);
        assert!(inner.is_some(), "the loop's inside is its own piece");
    }

    #[test]
    fn curved_strands_walk_like_straight_ones() {
        // A wavy section spanning the face: the angular sort reads where
        // each strand leaves a circle about its node, so curvature is
        // invisible to the walk.
        let strands = vec![
            Strand {
                polyline: vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
                tag: 0,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 0.0), Point2::new(4.0, 2.0)],
                tag: 1,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 2.0), Point2::new(4.0, 4.0)],
                tag: 2,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 4.0), Point2::new(0.0, 4.0)],
                tag: 3,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 4.0), Point2::new(0.0, 2.0)],
                tag: 4,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 2.0), Point2::new(0.0, 0.0)],
                tag: 5,
                boundary: true,
            },
            Strand {
                polyline: vec![
                    Point2::new(0.0, 2.0),
                    Point2::new(1.0, 2.4),
                    Point2::new(2.0, 2.0),
                    Point2::new(3.0, 1.6),
                    Point2::new(4.0, 2.0),
                ],
                tag: 6,
                boundary: false,
            },
        ];
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
    }
}
