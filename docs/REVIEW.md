# Kernel review: bugs, performance, memory

An adversarial read of the whole kernel (2026-09-27, at `cf744c4`), plus
measurements. Each item says where, what is wrong, the fix, and what guards
the fix. Nothing here changes behaviour that is right today: every fix is
either a pure speedup with identical output, or a change on inputs that are
wrong today.

How to use it: work down the order in [Plan](#plan). When an item lands,
change its status to `done <commit>`. Re-run `tools/ogeom-stress` with
`--check` after each performance item; a speedup that moves the baseline is
not a pure speedup.

Status keys: `open`, `done <commit>`, `dropped (reason)`.

Evidence keys:

- **reproduced**: a scratch program shows the failure.
- **read**: confirmed by reading the code path end to end.
- **plausible**: the mechanism is clear, the size is not measured.

Line numbers are as of `cf744c4` and will drift; the function names hold.

## Reproduced

These were run, not only read:

| Item | Input | Result |
|---|---|---|
| B1 | fuse an L-bracket with a disjoint box sitting in its notch | 1 solid, `check` valid; should be 2 solids |
| R1 | STEP entity with fewer arguments than expected (`VERTEX_POINT('')`) | panic, index out of bounds |
| R2 | `#5=ORIENTED_CLOSED_SHELL('',*,#5,.T.)` | stack overflow, process aborts |
| R2 | 200k nested parentheses in one entity | stack overflow, process aborts |
| R3 | one entity numbered `#999999999999` | aborts allocating the visited table |
| B2, B3 | extrema of two unit spheres 10 apart (want 8, 10, 10, 12) | 8, 9.05 (four times), 9.9999, 10; the farthest 12 is missing and 9.05 is not stationary |
| B4 | two planes through the origin, 1e-6 rad apart | `Meeting::Same` |

## Correctness bugs

### B1. Shells are nested by bounding box: disjoint solids become voids
`open` `reproduced`

- Where: `ogeom-bool/src/lib.rs` `assemble_result` (the "Nest" block) and
  `make_volume`.
- Problem: shell j is taken as a void of shell i when i's box contains j's
  box. A separate part sitting in a notch, or an island inside a cavity, is
  swallowed as a void of the wrong solid. `check` does not notice (see B8).
  The two functions also disagree: `make_volume` reverses the void,
  `assemble_result` does not.
- Fix: box test as a prefilter only; confirm by classifying a point of j
  against i; nest by depth parity (0 solid, 1 void, 2 solid again). One
  shared helper for both functions, one orientation rule.
- Guard: L plus notch box gives 2 solids of 750 and 8; hollow box plus
  island gives 2 solids.

### B2. Surface-surface extrema never polishes a farthest seed
`open` `reproduced`

- Where: `ogeom-intersect/src/extrema.rs`, the seed loops of
  surface/surface extrema and `thin`.
- Problem: near and far seeds alternate (even and odd indices); `thin` keeps
  every 8th, an even step, so every far seed is dropped and only a quarter
  of the near ones survive. Seeds are not local extrema either, so 256 4x4
  Newton solves always run.
- Fix: keep only local minima and maxima of each field over the grid (as
  curve/surface extrema already does), and thin near and far separately.
- Guard: two spheres report both nearest and farthest approach.

### B3. Extrema reports Newton results that never converged
`open` `reproduced`

- Where: `extrema.rs`, the three `newton_system(...).ok()?` sites.
- Problem: `newton_system` returns `Ok` when it gives up (exhausted, or no
  descent). The stall point is kept as a stationary approach. The walk and
  march re-check the gap; extrema does not.
- Fix: accept only a converged result, or a residual under a scale-aware
  bound (see N1).

### B4. Plane-plane parallel test is a million times too loose
`open` `reproduced`

- Where: `ogeom-intersect/src/surface.rs` `plane_plane`; the same formula in
  `plane_cylinder`.
- Problem: `(|cos| - 1).abs() <= angular` compares 1 - cos(t), about t^2/2,
  to an angle. Planes up to 1.4e-6 rad apart are "parallel", and sharing a
  point they are called `Same`: at 100 mm they are 1e-4 apart, 1000 times
  the confusion.
- Fix: `na.cross(nb).magnitude() <= angular`, as the rest of the crate does.
- Guard: two planes at 1e-6 rad through one point meet in a line; rerun the
  coplanar boolean tests.

### B5. Drafted faces about an oblique neutral can come out inside out
`open` `read`

- Where: `ogeom-offset/src/draft.rs` `general_draft`.
- Problem: the hinge chain runs in whatever direction the mesh segments
  chained; the new surface's normal follows it, while the face keeps the old
  face's orientation flag. Nothing compares the two. This is the inverted
  "drafted drum about an oblique neutral" already seen by the orientation
  probe.
- Fix: compare the new normal at the middle station with the old one
  (`ogeom_algo::normals_oppose`, as `heal/remodel.rs` does) and reverse the
  hinge if they oppose.
- Guard: oblique-neutral draft of a cylinder, both pull signs, valid and
  positive volume.

### B6. Canonical recognition can flip concave faces
`open` `read`

- Where: `ogeom-heal/src/canonical.rs`, where the recognised carrier
  replaces the spline.
- Problem: a cylinder, cone or sphere always has its natural normal away
  from the axis or centre. A bore stored as a spline with its normal toward
  the axis keeps its Forward flag and ends up inside out.
- Fix: `normals_oppose(old, carrier)` and reverse the rebuilt face when
  true (`remodel.rs` already does this).
- Guard: a spline bore through `simplify_to_canonical`: valid, volume
  unchanged.

### B7. Marched fillet decides its side at the wrong station
`open` `plausible`

- Where: `ogeom-fillet/src/marched.rs`, the inward/outward tests of the
  closed and open bands.
- Problem: the surface point is taken at the middle of the chord-length
  parameter, but compared with the ball centre of station `n/2`. Adaptive
  steps bunch stations where curvature is high, so the two can be far apart
  and the sign can flip. Likely source of the inverted marched fillet on a
  spline edge.
- Fix: take the point and normal at station `n/2` itself, or vote over a
  few stations.
- Guard: convex fillet on a closed elliptical rim: valid, volume matches
  the swept-ball estimate.

### B8. `check` misses inside-out faces and non-manifold shells
`open` `read` (known gap)

- Where: `ogeom-algo/src/check.rs`; `build.rs` `is_shell_closed`.
- Problem: nothing in `check` looks at orientation, and shells count as
  closed when every edge is used an even number of times, so an edge used 4
  times passes. B1's output passes `check`.
- Fix: run `mass::flags_agree` per closed shell (it answers true whenever it
  cannot tell, so no false alarms) and report `Broken` when false; report
  edges used more than twice. Plan from the stress notes: `flags_agree`
  false and the step-off probe confirms, then report.
- Guard: flip one box face; two cubes sharing one edge; B1's output.

### B9. Skinned walls orient against the centroid of the whole skin
`open` `read`

- Where: `ogeom-offset/src/sweep.rs` `skinned_wall`, `skinned_solid`,
  `apex_patch`, `cornered_loft`.
- Problem: the wall is kept or reversed by comparing its normal with the
  direction from the centroid of all rows. For a pipe on a semicircle in the
  XZ plane the centroid sits inside the bend and the test reverses a correct
  wall. XY spines hide it, which is why tests pass.
- Fix: test against the centre of the row nearest the middle, or orient the
  closed shell once by signed volume.
- Guard: `make_pipe_skinned` on semicircles in XY, XZ, YZ.

### B10. Unevenly spaced loft sections lose about 20% of the volume
`open` `plausible` (known bug, cause narrowed)

- Where: `cornered_loft` > `skinned_strip` > `fit_surface_grid`
  (`ogeom-geom/src/fit.rs`), centripetal parameters.
- Problem: a square twisting a quarter turn with steps 0.625 and 0.3125
  measures 64 instead of 80. The fit passes through every row but its error
  is measured only at the rows; centripetal parameters squash a 2:1 spacing
  change and the cubic sags between rows. Worked round in `4c2d481` by
  spacing stations evenly.
- Fix: first measure mid-span error to confirm; then chord-length
  parameters from section spacing, and a mid-span check in the fit.
- Guard: the reported case gives 80; all loft volume tests.

### B11. Smaller correctness bugs

| Id | Where | Problem | Fix |
|---|---|---|---|
| B11a | `ogeom-bool` `fill`, the `break` after `march_pair` | marching one curve drops the pair's other curves and can duplicate earlier ones | collect per pair, replace only non-tangent entries |
| B11b | `ogeom-bool` `contact_intervals` | trim outline lacks the second seam column and end welding, so tangency contacts on full drums read as outside | reuse `fill`'s outline (see P3) |
| B11c | `ogeom-bool` two contact paves | use the raw edge tolerance, not `honest()`; a loose imported edge merges every pave on it | `honest(contact.tolerance, tol)` |
| B11d | `ogeom-bool` `section()` | fresh vertices per sub-edge, so section edges never connect | weld by position |
| B11e | `ogeom-bool` `Rebuild::vertex`, `strands_of` | overlapping junctions resolve to the first by index, not the nearest | `min_by` distance |
| B11f | `ogeom-bool` `is_half_space` call in `fuse` | hardcodes `Tolerances::millimetres()` | pass `tol` |
| B11g | `ogeom-topo` `check_tolerances` | face-over-edge rule never checked (wires carry no tolerance, so the pair is skipped) | carry the nearest ancestor tolerance down |
| B11h | `ogeom-geom` trimmed surface and curve `transformed` | trim window not rescaled when the basis domain starts at 0, or when trims nest | map the trim affinely from old to new domain |
| B11i | `ogeom-geom` `traits.rs` periodic normalisation | NaN and infinity pass as `Ok(NaN)` on circles, cylinders, spheres | reject non-finite first |
| B11j | `ogeom-algo` `sew.rs` `for ... in merged.clone()` | iterates a `HashMap`: pcurve ids and repr order vary run to run | sort by `dropped.index()` |
| B11k | `ogeom-algo` mass | volume swallows exact-integration errors (mesh fallback), area propagates them; `.ok()` also swallows `Cancelled` and `Dangling` | one helper: geometric errors fall back, others propagate |
| B11l | `ogeom-hlr/src/exact.rs` `sampled` | closed silhouettes have span near 0, so they are drawn as octagons | use `ogeom_mesh::discretize` |
| B11m | `ogeom-mesh` `triangulate.rs`, pcurve-only edge fallback | passes `deflection` instead of the agreed `along` chord | pass `along` |
| B11n | `ogeom-io` STEP strings | writer emits raw UTF-8, reader reads Latin-1, `\X2\` escapes never decoded: "ø" round-trips as "Ã¸" | encode and decode Part 21 escapes |
| B11o | `ogeom-io` STEP writer `real()` | writes `NaN.0`, `inf.0` | refuse non-finite |
| B11p | `ogeom-fillet` `marched.rs` | overwrites the input face's shared surface before steps that can fail | add a new surface, commit after success |
| B11q | `ogeom-intersect` `walk.rs` loop closure | closes when within one step of the start; a hairpin or the next thread turn can close early | also require heading agreement |

## Numerics

| Id | Where | Problem | Fix |
|---|---|---|---|
| N1 | `extrema.rs` (`confusion^2`), `curves.rs` `polish_3d` (`confusion * 0.01`) | stopping thresholds are absolute on a residual scaled by the parametrisation speed: unreachable for a fast spline, loose for a slow one | normalise by the derivative norms |
| N2 | `curve_surface.rs` line/sphere, line/cylinder; `solve.rs` `quadratic_roots` | an exact tangency rounds the discriminant negative about half the time and is missed | treat a slightly negative discriminant as a double root; the existing polish and gap filter reject false ones |
| N3 | `walk.rs`, `march.rs`, `curve_surface.rs`, `curves.rs` Newton closures | a failed evaluation becomes a zero residual (`unwrap_or(ORIGIN)`), which Newton can accept | return infinite residuals so the line search backs off |

## Hostile input (a reader must never abort)

### R1. STEP reader indexes argument lists directly
`open` `reproduced`

- Where: `ogeom-io/src/step/read.rs`, about 25 `args[N]` sites in
  `surface`, `curve`, vertex, edge and face builders.
- Fix: `args.get(i)` like the rest of the reader, or a `ref_at` helper.
- Guard: every keyword with an empty argument list returns `Err`.

### R2. Unbounded recursion in STEP and IGES
`open` `reproduced`

- Where: STEP `parse.rs` `arguments`/`argument`; `read.rs` `surface`,
  `curve`, `shell` for offset, trimmed, replica and oriented entities. IGES
  `read.rs` `placement` (transform chain).
- Fix: an in-progress set or depth cap per builder; a nesting cap in the
  parser. `colour_in` and `datum_letter` already cap depth.

### R3. File-supplied numbers drive allocations
`open` `reproduced` (one case)

- STEP `visited` is sized by the largest entity id; negative ids wrap.
- STEP knot multiplicities and implied degree; IGES counts (`n_knots`,
  `3*n+1`, `nu*nv`); native and BREP counts; glTF accessor counts: all
  allocate or loop before checking against the data present.
- 3MF: ZIP64 sizes pre-allocate, inflate has no output cap (zip bomb),
  offset sums unchecked.
- Fix: validate each count against the bytes or parameters present,
  `checked_mul`/`checked_add`, cap `with_capacity`, bail as soon as inflate
  output passes the declared size. Dense `visited` only when ids are dense,
  else a set.

### R4. IGES slices text by byte
`open` `read`

- Where: `iges/parse.rs` fixed-column slicing; Hollerith length `j + 1 + n`.
- Problem: a non-ASCII character near column 64 or 72 panics (not a char
  boundary); a huge Hollerith count overflows.
- Fix: parse as bytes, decode Hollerith payloads lossily, checked adds.

### R5. Smaller input hazards

- `ogeom-topo` `Location::composed` multiplies once per unit of power, and
  `then` adds powers unchecked: a native file with `d0^2147483647` hangs or
  wraps. Use exponentiation by squaring and `checked_add`; reject extreme
  powers in the reader.
- `ogeom-topo/src/tessellation.rs` `welded` adds to grid keys unchecked;
  infinite STL coordinates overflow. Use `saturating_add` as 3MF does.
- STL binary `claimed * 50` overflows on 32-bit targets.
- `ogeom-doc` `contains_product` has no visited set: exponential on shared
  sub-assemblies, and recursive.

Add a `cargo fuzz` target per reader (`read_step`, `read_iges`, STL, 3MF)
once R1 to R4 are fixed, so they stay fixed.

## Performance

Grouped by where the time goes. "Identical" means the output should be
bit-for-bit the same; check with the stress baseline and the test suite.

### Walks and topology (every operation pays these)

| Id | Where | Problem | Fix | Output |
|---|---|---|---|---|
| P1 | `ogeom-topo` `explore` | descends below the wanted type: `OfType(Face)` walks every wire, edge and vertex too | descend only into compounds and kinds above the wanted one | identical |
| P2 | `ogeom-topo` `Shape::moved`, `composed` | three location clones per child per step | build the child shape field by field; identity fast path in `then` | identical |
| P3 | `ogeom-topo` `ancestors_of` | explores every candidate's subtree per call; offset `growing_edges` calls it per edge, O(E x F x subtree) | an ancestor map built in one walk | identical order if built in explore order |
| P4 | `ogeom-topo` `Model::widen` | allocates and walks the full subtree on every `add_edge` and `add_face`, almost always a no-op | stop where the tolerance already suffices (after B11g) | identical |
| P5 | `ogeom-core` `parallel::map_ordered` | queries `available_parallelism` and spawns threads on every call, even for 2 items | cache the count; minimum-work threshold | identical |

### Geometry evaluation (the inner loop of everything)

| Id | Where | Problem | Fix |
|---|---|---|---|
| P6 | `ogeom-geom` `transformed` on B-splines | no identity fast path; clones the whole grid or net then discards it | identity returns a clone; clone only knots; a `Cow` for read-only callers |
| P7 | `ogeom-math` `bspline::evaluate`, `derivatives` | a heap `Vec` per de Boor evaluation and per derivative call | `SmallVec<[P; 8]>` |
| P8 | `ogeom-geom` B-spline surface `d1_at`/`d2_at` | always takes the rational path; `curvature_at` evaluates twice | branch on `rational` as `jet_at` does |
| P9 | `ogeom-math` `surface_derivatives` | computes the full square table (k, l up to order) and repeats the inner v sums | only k + l <= order; factor the v sum |
| P10 | `ogeom-geom` `BSplineCurve` rebuild methods | `..self.clone()` clones the control net only to replace it; `segment` costs about 5 net clones | write the fields out |
| P11 | `ogeom-geom` periodic normalisation | every periodic evaluation goes through a `#[cold]` function | inline fast path when already in range |

### Intersection and solving

| Id | Where | Problem | Fix |
|---|---|---|---|
| P12 | `ogeom-math/src/solve.rs` `newton_system` | 15 to 25 heap allocations per evaluation (DMatrix from `Vec<Vec>`), used by walk, march, extrema, curves | const-generic fixed-size solver (`newton_system_2` already exists) |
| P13 | `ogeom-intersect` walk and march | 6 surface evaluations per Newton step where 2 do; extrema calls `point_at`, `d1_at`, `d2_at` instead of `jet_at` | one evaluation per surface per step; `jet_at` |
| P14 | `ogeom-intersect/src/walk.rs` | the tangent at each accepted point is computed twice | carry it to the next iteration (identical) |
| P15 | `ogeom-intersect` `walk.rs` null vector | recursive cofactor determinant with allocations; absolute degeneracy test | fixed-size, relative test |
| P16 | `ogeom-intersect` `march.rs` `seeds` | re-samples each surface once per spline border; 1.3M cell-pair box tests unindexed; `sample_by` evaluates shared corners 4 times | sample once, grid-evaluate once, sweep-and-prune |
| P17 | `ogeom-math` `solve.rs` damped Newton | near-singular Jacobians take huge steps, then up to 30 halvings, each a full residual and Jacobian | relative pivot test, step cap, residual-only line search |
| P18 | `ogeom-intersect` `curves.rs` seeding | one global reach for all 129 x 129 segment pairs | per-pair reach |

### Boolean

| Id | Where | Problem | Fix |
|---|---|---|---|
| P19 | `fill` face-pair loop | sequential, while the cheaper paving stage is parallel; `surface_surface` solved twice per closed-form pair | `map_ordered` over admitted pairs, concatenate in order |
| P20 | `strands_of` and the arrangement loop | sequential per face | `map_ordered` per face |
| P21 | `chart_point_of` | rebuilds the face outline on every call, inside nested piece loops | store the welded outline on `GFace` once (also fixes B11b) |
| P22 | contact edges x target edges | no bound check; up to 8,000 curve evaluations per pair | prefilter by edge bounds |
| P23 | hug test | `distance_to_edge_curve` on every edge without a box test; vote vectors allocated per candidate | box prefilter; index contacts by face |
| P24 | all-pairs scans: sections x sections, per-face scans of pieces, contacts and junctions, dust, hug dedupe | pre-bucket by face; `Bins` for junctions; flags instead of `contains` |
| P25 | `arrange.rs` | node snapping is O(n^2) and takes the first node, not the nearest; layer-by-layer pruning; hole nesting positives x negatives x positives with sets rebuilt per test; `atan2` inside the sort comparator | grid snap, queue peel, precomputed cycle data |
| P26 | `cells()` | runs the general fuse three times | fuse once, assemble three filters |
| P27 | `face_trim_lines` | recomputed per rebuilt sub-edge | once per face |
| P28 | face bounds of curved faces | padded by 0.75 of their own diagonal | intersect with the surface's own extent (control net hull, sphere box) |
| P29 | `section_face` | clones face data and surface; linear vertex search; `classify_on_face` per candidate | no clone when identity; bins; prepared face |

### Classification, mass, check, sewing

| Id | Where | Problem | Fix |
|---|---|---|---|
| P30 | `ogeom-algo` `SolidBoundary::holds` | calls `classify_on_face`, which rebuilds the rings `prepare` already stored (about 600 times the ray cast) for every point near a face | classify on the prepared face |
| P31 | `ogeom-algo` `classify_in_solid` | re-triangulates the whole solid per call; fillet, ruled, marched and draft call it in loops of up to 32 | a prepared `SolidMesh` with a box early-out |
| P32 | `ogeom-algo` `project_on_surface` | 33 x 33 grid even for planes, cylinders, cones, spheres, tori | closed-form seed, grid only as fallback |
| P33 | `ogeom-algo` `mass_chart` `integrate` | samples grow with fine^2 and are all stored (hundreds of MB at the last doubling); no cancel checkpoint | stream the comparison; knot lines once; checkpoint |
| P34 | `ogeom-algo` `mass.rs` `gauss2`, `Accumulator::add` | 22 quadrature passes to recover a rule; a heap `Vec` per simplex | `gauss_legendre_rule`; fixed array (identical) |
| P35 | `ogeom-algo` `flags_agree` | walks every face again and clones and transforms each edge curve per sample | cache per edge |
| P36 | `ogeom-algo` `check` | self-intersection rebuilds each face per pair, no box filter, adjacency by hash (a collision skips a pair); containment reports duplicates; six `explore_unique` passes | build once, sweep-and-prune, `TShapeId` sets, one walk |
| P37 | `ogeom-algo` `sew` | rebuilds all touched edges each round (orphans the previous), fingerprints clone curves, bins sized by the loosest edge | rebuild only changed, per-edge radius |
| P38 | `ogeom-algo` `shape_bounds` | recomputes shared sub-shapes; `SolidBoundary` computes face bounds twice | memo per call |
| P39 | `ogeom-algo` `mesh_solid` `split_across_slivers` | scans all triangles per sliver per round | edge-to-triangle map |

### Fillet, offset, sweep, heal

| Id | Where | Problem | Fix |
|---|---|---|---|
| P40 | `ogeom-geom/src/fit.rs` `fit_family` | rebuilds and inverts a dense normal matrix per row per round; 2000 stations is a 32 MB matrix inverted per control column | factor once per round, banded solve, no inverse |
| P41 | `sweep.rs` `pipe_shell_law` | asks `evenly` for the total station count per edge (E times too many); `densified` can double 12 times on a jumping law | split by length, cap, fail rather than grow |
| P42 | pipes, evolved, drill tools | pieces fused one at a time onto the growing result, O(n^2) | balanced fuse, or cut a compound once |
| P43 | `fillet/src/corner.rs` labellings | each attempt runs the full-solid cut | test the tool cheaply first, cut once |
| P44 | `marched.rs` chord `3e-6` | absolute; a large rim hits the 20,000 point cap and fails | scale with radius |
| P45 | `shape.rs` offset guard | retries meshing at absolute chord `1e-4` | scale to the part; keep the first error |
| P46 | `heal/fix_shape.rs` | meshes every edge to measure length | endpoint distance early out |
| P47 | `fillet.rs` `refind_edges` | full-solid projections, twice per target | box prefilter, one edge-face map |

### Exchange, mesh, HLR

| Id | Where | Problem | Fix |
|---|---|---|---|
| P48 | `ogeom-mesh` `add_interior_points` | tests each grid point against every ring edge while `RingBands` is built just before and unused | pass the bands (identical, one line) |
| P49 | `ogeom-mesh` refinement loop | re-measures unchanged triangles every round, 9 evaluations each | memo by edge, only faces near new vertices |
| P50 | `ogeom-mesh` shared edges | sampled once per adjacent face, again on redraw | a shared edge sample map |
| P51 | `ogeom-hlr/src/project.rs` `occluded` | every sample against every triangle, one thread | 2D grid of projected triangles, then parallel |
| P52 | `ogeom-hlr/src/exact.rs` | full curve-surface intersection against every face per sample; rings recomputed per candidate | projected box per face, rings once |
| P53 | STEP reader | `args()` deep-clones argument trees per access; face surfaces built twice; `BuiltEdge` cloned per cache hit; about 20 full scans by keyword | borrow args; cache surfaces by id; keyword index built once |
| P54 | STEP writer | clones and scans `written_nodes` per coloured product; whole file built as strings | map by node; `write_step_to(impl Write)` |
| P55 | `ogeom-mesh` `simplify` | O(T) scan per collapse | heap with lazy invalidation |

## Memory

| Id | Where | Problem | Fix |
|---|---|---|---|
| M1 | `ogeom-topo` model arenas | nothing is ever freed; failed attempts and retries (corner labellings, reaching retries, sew rounds) leave nodes and geometry behind | rollback on failure, or a reachability compaction into a fresh dense model |
| M2 | `mass_chart` | see P33, the largest transient allocation found | stream |
| M3 | STEP parse | an owned `String` per keyword, 56-byte `Arg`, `HashMap<u64, Instance>` with SipHash: several times the file size | borrow from the input, box `Typed`, dense index |
| M4 | walk states | `Vec<Vec<f64>>`, one allocation per point, plus a clone per push | flat stride storage |
| M5 | `ogeom-bool` `gather` | clones and transforms every surface and every edge curve per adjacent face, even at identity | skip identity; cache world curves by edge |
| M6 | `ogeom-topo` identity map | `set_derived` orphans the primitive provenance record; `shape_of` is a linear scan | `Vec` indexed by node, provenance up front |
| M7 | `ogeom-doc` undo | full deep copy per checkpoint, unbounded | cap history, share unchanged fields |

## Measured

Pending: the profiling run (stress `--times`, `perf` on the slowest cases,
peak RSS) is filled in here.

## Plan

In order. Each step is one or a few commits, each gated by `check.sh` and
the stress `--check`.

1. **Readers never abort**: R1, R2, R3, R4, R5, then fuzz targets. Small,
   isolated, reproduced.
2. **Wrong results**: B1, B4, B2 and B3 with N1, B11c, B11f, B11j
   (determinism), B11i, then N2, N3. Each with the guard test named.
3. **Orientation**: B8 first (so `check` can see the rest), then B5, B6,
   B7, B9. Re-run the orientation probe over the test models.
4. **Identical-output speedups, cheap**: P1, P2, P14, P30, P34, P48, P6,
   P7, P10, P21, P27. Expect the broadest gain from P1, P30 and P7.
5. **Parallel boolean**: P19, P20 (deterministic through `map_ordered`),
   after P5.
6. **Solver and evaluation**: P12, P13, P8, P9, P15, P16, P17.
7. **Fitting**: P40, then B10 with the chord-length change, then P41.
8. **Indexing**: P22 to P25, P36, P37, P39, P51, P52, P55.
9. **Memory**: M1 design first (rollback or compaction), then M3, M5.

Out of scope here and already tracked elsewhere: the tracer landing exactly
on a domain edge (the general form of the #74 fix), the boolean missing a
section along a patch boundary edge (#81), the drill scan classes.
