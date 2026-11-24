# M6: remove the old blend paths without losing regressions

## Status

Complete. On the final layer, with `RUST_TEST_THREADS=1` and both Cargo
invocations exiting 0, the counters recorded zero round and zero chamfer
requests over every `ogeom-fillet` and `ogeom` test, ignored ones included.
A workspace search found no caller of the old constructors. The old round and
chamfer construction, the `old-path-log` feature, `to_current_path` and the
`tools/old-path-*.sh` scripts are removed; the inventory below is kept as the
record of what was retired and is pinned to its own revision.

## Decision and baseline

**Do not delete the old round implementation yet.** Retire construction, not
public APIs or regression coverage. ADR [0001](adr/0001-blend-architecture.md)
§6 and its “what M6 refuses by name” amendment are the acceptance contract;
[#902](https://github.com/gilbertorconde/ogeom-rs/issues/902) is the tracker.
A formerly accepted solid becoming a named refusal is still a migration
blocker. Bit identity to an old fitted patch is not an acceptance oracle.

This inventory is pinned to `origin/main` at
`46d545d1c853d3934cf5d4f827ef1420427de429`, inspected on 2026-11-15. All
source locations below refer to that revision. Unless a full path is given,
`src/` and `tests/` mean `crates/ogeom-fillet/src/` and
`crates/ogeom-fillet/tests/`. A module location at line 1 names the whole file;
function and test locations name their definitions or dispatch sites.

### Measured reachability

Both repository scripts were run unchanged, including ignored tests, with
`RUST_TEST_THREADS=1`. Cargo's actual exit status and output were retained
separately: **all three underlying Cargo invocations exited 0**. The fillet
suite passed 532 tests on each invocation; the facade suite passed 1,244.

| Measurement | Requests | Distinct asking tests | Scope |
|---|---:|---:|---|
| `tools/old-path-rounds.sh -- --include-ignored` | **214** | **77** | All `ogeom-fillet` tests, including unit tests and ignored tests |
| `tools/old-path-chamfers.sh -- --include-ignored` | **0** | **0** | All `ogeom-fillet` and `ogeom` tests |
| Additional round records in the chamfer script's second Cargo run | **220** | **27** | `ogeom` tests, not included by the round script |

The additional 220 records are obtained by subtracting the first fillet run's
records from the chamfer script's cumulative log, not by counting chamfers as
rounds. They comprise 189 `fillet_edge` and 31 `fillet_edges` calls. They are
an additional M6 gate, not part of the 214-request baseline.

The supplied 198/74 expectation is not this main revision's count.
[Open PR #1092](https://github.com/gilbertorconde/ogeom-rs/pull/1092), head
`f9d61538af073364f76f08fd756fc6d766fd6f6d` when inspected, reports **214/77
on its baseline and 198/74 on its own branch**. Its three eliminated asking
tests are the narrow-face and side-cut-off cases in `tests/fillet/fillet.rs`, and
`rounds_at_inside_corners_off_square_or_leaning_both_ways_take_all_above_both`
in `tests/fillet/chain.rs`. This document neither imports nor duplicates its work.
Recount after it and other routing PRs merge; do not subtract branch claims
from a changing main.

Two limitations matter before trusting a zero:

1. At this baseline, both scripts discard Cargo output and use `|| true`
   (`tools/old-path-rounds.sh:19`, `tools/old-path-chamfers.sh:19–23`). A failed
   compile can report zero. The logger at `src/fillet.rs:102` can also
   interleave records between test threads. An initial two-thread run
   reported 214/78, including a concatenated, nonexistent test name. The
   serial 214/77 result is the usable baseline. PR #1092 contains counter
   reliability fixes; land those first or retain an independent exit-status
   record and use one test thread.
2. The counters log dispatch, not every old constructor. Engine-to-old
   constructor calls, direct vertex APIs, linear-variable APIs and private
   wedge tests bypass the logger. Zero dispatches is necessary, not sufficient.

## 1. Entry points and the exact branches to remove

The public facade is `crates/ogeom/src/lib.rs:17`, which re-exports the
**entire** fillet crate. Keep its API and the exports at `src/lib.rs:52–87`.
Production callers below include library-to-library calls; applications outside
this checkout cannot be enumerated and are not assumed absent. The test
appendices give every observed asking test, with its source location.

### Logged round entry points

| Entry point / fallback site | Current callers and tests | Removable implementation after migration |
|---|---|---|
| `fillet_edge`, `src/fillet.rs:80`; fallback `:91–92` | Public facade; `asym::fillet_edge_distances`' equal-round branch; facade vertex/feature tests. 18 logged requests in the fillet suite; 189 additional facade requests. | The post-engine call to `fillet_edge_wedge(:126)`, not the API. Then its private `cut_off_bands(:153)`, `fillet_edge_own(:208)`, `fillet_edge_meeting(:425)` and seat-routing construction once their other callers disappear. |
| `fillet_edges`, `src/fillet.rs:1062`; fallback `:1073` | `blend_symmetric`, checked/partial wrappers, setback composition, `facepair::curved_blend`; integration families in Appendix A. 101 logged fillet-suite requests, 31 additional facade requests. | Post-engine `on_own_nodes → cut_off_bands → fillet_edges_meeting(:1613)`. Preserve copy-on-write and history semantics in the engine adapter. |
| `fillet_edges_with`, `src/fillet.rs:1110`; fallback `:1135` | Mixed blend requests and checked/partial operations; `radius_per_edge.rs`, asymmetric corner fixtures. 2 logged requests. | The same meeting/corner/sequential band tail, not per-edge radius validation. |
| `fillet_edges_variable`, `src/fillet.rs:3423`; fallback `:3452` | Public API and `tests/fillet/variable.rs`. 35 logged requests. | `fillet_chain_law_own(:3677)`, `law_piece(:3781)`, `consumed_law(:3913)`, `chain_law(:4192)`, `run_on_past_turns(:4686)`, `straight_law_fillet(:4810)` and their private law-construction graph after variable routing covers the accepted domain. Keep law evaluation. |
| `fillet_edge_distances`, `src/asym.rs:364`; fallback `:437` | Public API; level laws in `asym_law.rs:74`; asymmetric tests. 47 logged requests. | `fillet_edge_distances_own(:452)`, rim-tool construction and the `marched_bevel_with`/`planar_fillet` tail after conic routes migrate. Keep `RoundProfile` and section arithmetic. |
| `blend_edges`, `src/blend.rs:51`; asymmetric fallback `:102` | Public mixed-round/bevel API; `tests/fillet/asymmetric.rs` and `blend_edges.rs`. 11 logged requests. | `asymmetric_corners`, `after_symmetric` and `close_asymmetric_corners` dispatch in `src/asym.rs:1272,679,1484`. Keep the API and the engine-owned bevel-then-round semantics of `blend_symmetric(:147)`. |

The 214 requests split exactly as **18 + 101 + 2 + 35 + 47 + 11**. These
are requests, not edges, successful builds, or one count per test. Refusal
requests also have to move to the engine with the same meaningful diagnosis.

### Logged chamfer entry points: zero observed, not zero static callers

| Entry point / fallback site | Current callers and tests | Removable implementation |
|---|---|---|
| `chamfer_edge`, `src/chamfer.rs:54`; `chamfer_edge_own(:70)`, fallback `:84–85` | Public facade; chamfer, curved-chamfer and feature tests. | `wedge_for(...).apply(...)` after engine decline. Keep seam no-op and ownership behavior. |
| `chamfer_edge_distances`, `src/chamfer.rs:101`; own `:121`, fallback `:142–143` | Public facade; asymmetric/distance chamfer tests. | Same wedge tail; keep named-face distance interpretation. |
| `chamfer_edge_angle`, `src/chamfer.rs:160`; own `:180`, fallback `:201–202` | Public facade; angle chamfer tests. | Same wedge tail; keep angle-to-legs conversion. |
| `chamfer_edges(:290)`, `chamfer_edges_with(:354)` → `chamfer_edges_blamed(:368)`; fallback `:400–406` | `blame::chamfer_edges_checked(:98)`, partial operations, mixed blends; all chamfer suites and facade blends. | `chamfer_edges_mitred_blamed(:1282)` → `chamfer_edges_with_own(:1313)` sequential application, preserving blame-index remapping. |
| `chamfer_edges_cornered`, `src/chamfer_corner.rs:65`; fallbacks `:86–87` and `:95–96` | `blend_symmetric`; corner and polyhedral chamfer tests. | Old point-mitre call `chamfer_edges_mitred(:1272)` and the old facets/corner-tool tail `facets(:128)`, `polygon_corner(:266)`, `tetrahedron(:494)`. Keep `ChamferCorner`, including Point and Facet behavior. |

The old chamfer graph below those tails includes `wedge_for(:1444)`,
`sphere_tool(:1510)`, `revolved_bevel(:1981)`, `revolved_bevel_over(:1996)`,
`sector_bevel(:2138)`, `bevel(:2224)` and `marched_step(:2332)` in
`src/chamfer.rs`. Its production callers are the tails above and its own
private builders. `support::tests::blend_wedges_walk_each_edge_once_each_way`
directly calls `chamfer_edges_mitred` at `src/support.rs:3127`, bypassing the
zero counter. Port that test's geometric assertion before removing this seam.

### Unlogged entry points and engine back-edges

| Entry point / location | What still reaches it | Required treatment |
|---|---|---|
| `fillet_edge_variable(:3063)` → `fillet_edge_variable_own(:3082)` in `src/fillet.rs` | `tests/fillet/fillet.rs:446`, several `tests/fillet/variable.rs` comparisons, `crates/ogeom/tests/facade/mechanical_part.rs:95`. Equal endpoints call old `fillet_edge_own` directly. | Route the linear law through the engine, including equal endpoints; delete the old linear solid builder only after those assertions pass. Not counted today. |
| `fillet_edge_distances_variable`, `src/asym_law.rs:56` → `varying(:92)` | `tests/fillet/asymmetric.rs:1375–1611`; constant laws call `fillet_edge_distances`, nonconstant laws build directly. | Route varying conic distances; keep `varying_net(:323)`, which the engine's section layer uses. |
| `round_vertex`, `src/corner.rs:68` → `round_vertex_own(:84)` → `round_vertex_with(:97)` | Public facade, #981 corner-history tests; old fillet, asymmetric, setback, rolled and sharp-corner builders; `crates/ogeom/tests/facade/blends.rs:2380–3385`. | Preserve the vertex-only API's result, including its planar caps, and generated-face history. It is not equivalent to rounding all incident edges. Replace its construction before deleting `corner.rs` / `corner_curved.rs`. |
| `round_vertex_setback`, `src/setback.rs:64`; `round_vertex_setback_sharp`, `src/sharp.rs:95` | Public facade; `tests/fillet/setback.rs` and `cell_guard.rs`; old corner orchestration. | Port public vertex-only semantics and continuity/refusal contracts. Do not remove them because edge-request counts reach zero. |
| `fillet_edges_cornered`, `src/setback.rs:256`; fallback in `cornered_own(:294)` at `:350` | Public facade; `tests/fillet/setback.rs`, `cell_guard.rs`. | Keep validation and engine-first `cornered_by_engine`; replace old `fillet_edges` followed by `close_corner`, which is not separately logged. |
| `request_fits`, `src/fillet.rs:1831` | `src/blame.rs:68,140`, `src/fits.rs:508`, old fillet construction, `marched_wedges(:638)`. | Retain/move the geometric fit and overlap guards and blame data until all production callers use an equivalent engine guard. These are not test-only. |
| `engine::marched_seats::rounds`, `src/engine/marched_seats.rs:29` | Last route from `rounds_by_engine_unplaced`, `src/fillet.rs:1407`; marched integration tests may return successfully **without a counter record**. | Calls old `marched_wedges(:36)`, `rim_ring_taken_whole(:45)`, `apply_set_aside(:54)`, `marched_strip_taken_whole(:83)`. Replace the old constructors with engine-owned tools before deleting their graph. Merely moving this whole case catalogue does not complete the redesign. |
| `marched_wedges`, `src/fillet.rs:618` | The engine back-edge above. It calls `route(:937)`, `request_fits`, `fillet_edge_meeting(:736)`, consumed-rim/marched constructors and `junction` joining. | Explicit M6 blocker even if every public fallback count becomes zero. Preserve accepted tools, fit targets, caps, ancestry and refusals while replacing it. |
| `engine::fallback`, `src/engine/fallback.rs:1` | `engine::run(:300,392,400,413,476)` and curved/marched engine routes. | **Keep.** This is the ADR's one-boolean engine assembly, not the old fallback dispatcher. `round_tool(:822)` and `group_tool(:1025)` call `super::run` at `:892,1158`, not `fillet_edge_wedge`. `straight_non_local(:1207)` calls shared `chamfer::legs` at `:1222`; move that reader, not this fallback. |
| Engine recursion / public routing adapter, `src/engine/mod.rs:340` | Engine defeaturing passes `fillet::rounds_by_engine` as callback. | Keep/move the engine-only adapter `src/fillet.rs:1152–1607`; it is not an old constructor. Remove its reverse dependency on an API implementation file when practical. |

Other APIs are not old-edge dispatchers: `blend_faces` (`src/facepair.rs:139`),
`fillet_sheet_edges` / `fillet_faces` (`src/sheet.rs:108,259`), marching,
bridging, filling, analysis and 2D corner operations remain supported. But
they retain shared pieces of old files: planar `blend_faces` calls
`fillet::seated_fillet` (`src/facepair.rs:30,297`), and curved face/sheet
blends use `pair_marched`, `sheet_curved` and `marched`. This prevents deleting
those files wholesale, even after the edge engine is independent.

## 2. Module deletion manifest

This is a file manifest, not a permission to delete a file whose retained
symbols still have callers. “Remove” means after the reachability, shared-helper
and test gates in this document. Each listed caller is production unless
explicitly marked test; the test column identifies the relevant owning suites.

| File / anchor | Callers / owning tests | Disposition |
|---|---|---|
| `src/fillet.rs:1` | Public round APIs; `asym`, `blame`, `fits`, `corner`, `facepair`, `chamfer`, `junction`, `support`, engine and section layer; round/variable/chain suites | Split API, engine routing and shared math first; remove old seat, sequential-band, run-out and junction orchestration. |
| `src/chamfer.rs:1` | Public APIs, `blame`, `blend`, `chamfer_corner`, `asym`, engine readers; all chamfer suites | Keep API/engine adapters and shared readers; remove wedge graph. |
| `src/chamfer_corner.rs:1` | Public API, mixed blends; `chamfer_corner.rs`, `polyhedral_chamfer.rs` tests | Keep enum/API; remove post-engine facets/mitres graph. |
| `src/asym.rs:1` | Public APIs, `blend`, `asym_law`, fillet, engine and `round/section`; `asymmetric.rs` | Keep profiles/section algebra and engine routing; remove old conic solids/corner orchestration. |
| `src/asym_law.rs:1` | Public varying-distance API; `round/section`; `asymmetric.rs` | Keep/move `varying_net`; remove direct solid builder only after routing. |
| `src/consumed.rs:51,404` | `fillet`, `rolled`, `stepped`, `carried`, `shelf`, `rolled_rim`, `support`, `setback`; consumed and narrow-face suites | Old consumed-seat/slanted tools removable with the old construction graph, after guard dependencies are discharged. Distinct from **retained** `src/engine/consumed.rs`. |
| `src/consumed_rim.rs:250` | `fillet`, including unlogged marched bridge; `consumed_rim.rs`, variable tests | Remove consumed-ring catalogue after replacing bridge and variable callers. |
| `src/consumed_marched.rs:48` | `fillet`, including unlogged bridge; `consumed_marched.rs` | Remove after engine tools replace consumed-marched construction. |
| `src/rolled.rs:465,1363,1385` | `fillet`, `support`, `setback`, `stepped`, `carried`, `shelf`, `rolled_rim`; concave/tangent/roll-over tests | Remove rolled construction after guard and unlogged bridge migration; extract reusable planar predicates if still called. |
| `src/carried.rs:641,1338` | `rolled`, `stepped`; consumed/roll-over/variable tests | Remove carried construction with its callers; preserve independent support-transfer tests. |
| `src/stepped.rs:123,746` | `rolled`; stepped/concave/narrow-face tests | Remove old stepped-run/band catalogue after rolled callers disappear. |
| `src/rolled_rim.rs:84,222` | `fillet::request_fits`, `rolled`; #981 carried-rim tests and roll-over suites | Retire builders only after fit/guard users and bridge are replaced; preserve bounded-contact regressions. |
| `src/shelf.rs:88,231` | `fillet::request_fits`, `rolled`; #981 shelf and rim-bounds regressions | Same treatment; zero public fallback logs do not prove this reader dead. |
| `src/ruled.rs:49,130,198` | `fillet::route`, `support`; contact/cap regressions | Remove old ruled solids once old routing and fit readers disappear. Engine's own ruled-seat construction remains. |
| `src/junction.rs:152,802,848,1006,1832,1869` | `fillet`, including `marched_wedges`; chain/tangent tests; unit test `:2311` | Remove junction/carry/join catalogue only after unlogged bridge and variable calls go. |
| `src/pinched.rs:75` | `marched`; tangent/pinched variable cases | Remove edge-specific pinched construction once old marched entry points disappear. |
| `src/turned_pole.rs:42` | `marched`; pole/mesh tests in `support` and facade | Remove edge-specific pole-leg construction only after migrated pole regressions pass. |
| `src/corner.rs:68,97` | Public vertex API and old corner callers; history and private labelling tests | Port vertex-only API; retire old tool/chart search, preserving #981 public ancestry. |
| `src/corner_curved.rs:38` | `corner`; curved-vertex facade tests | Remove when vertex-only curved behavior is ported, not at zero edge counters. |
| `src/setback.rs:64,256,1015,1726` | Public APIs, fillet, asymmetric, sharp/cusp, chamfer, engine; setback/mixed/cell-guard tests | Retain `Setbacks`, validation and engine route; move `point_of`/`trim_conic`; retire old fill/certification/corner closures after replacements. |
| `src/sharp.rs:95,250,438` | Public sharp-setback API, `setback`; sharp and cell-guard tests | Port API; retire old sharp corner sewing/fill. |
| `src/cusp.rs:90,234` | `setback`, `fillet`; mixed/tangent cusp tests | Remove after cusp behavior and unlogged construction are migrated. |
| `src/chamfer_marched.rs:151,164,485` | Old chamfer/asymmetric solids **and engine contact reader**; marched/crease chamfer tests | Keep/move `crease_contacts`; remove `marched_bevel[_with]` only after asymmetric callers move. |
| `src/chamfer_sphere.rs:175,309,346` | `chamfer` **and engine**; sphere-rim chamfer suites | Keep/move `SphereSeat`, seat and ring-surface reading; remove obsolete wedge application/sector construction only after symbol-level caller check. |
| `src/support.rs:1` | Nearly all retained geometry APIs, engine and old builders; private wedge tests | Never blanket-delete. Keep seats, geometry, distances and guards; remove collection/cut-off bookkeeping only after the unlogged bridge and tests stop needing it. |
| `src/marched.rs:54,69,3542` | Old edge builder; retained face/sheet, pair marching, fills, engine and round patch layers; private trim/contact tests | Split reusable fitting/chart/geometry from old edge construction. Keep `band_fit_target`, fitting and retained face/sheet callees. |
| `src/pair_marched.rs:77,242,1955,2097` | `facepair`, `sheet_curved`, `between`, old consumed-marched; contact-accuracy unit test | Keep face-pair and sheet functionality; only prune exclusively old-edge callers. |
| `src/between.rs:85` | `pair_marched`, `sheet_curved`, old shelf/rolled-rim; clearance unit test | Keep independent swept-ball clearance. |
| `src/fill_n.rs:195,269`; `src/fill_patch.rs:1` | Public `make_filling_n`, `round/patch`, old cusp/sharp/setback; filling and solver tests | Keep public filling and its solver; old corner call sites can go. |
| `src/facepair.rs:139`; `src/sheet.rs:108,259`; `src/sheet_curved.rs:1` | Public face/sheet APIs, `pair_marched`; facade `blends.rs`, `sheet_fillets.rs` | Keep. Extract shared old-file helpers before deleting any of their implementation. |
| `src/fits.rs:1`; `src/blame.rs:1`; `src/blend.rs:1` | Public largest-radius/distance, checked/partial and mixed APIs; largest/blame/blend suites | Keep API contracts; replace old guard/constructor dependencies, not the wrappers. |
| `src/law.rs:1`; `src/march.rs:1`; `src/bridge.rs:1`; `src/analyse.rs:1`; `src/corner2d.rs:1` | Public APIs; engine sections, contact construction and independent tests | Keep; not part of old-path deletion. |
| `src/round/{section,sphere,transfinite}.rs:1`; `src/round/mod.rs:1`; `src/engine/` | Engine section/vertex/assembly pipeline and tests | Keep. The unused fixed-fill prototype `src/round/patch.rs` is a separate pruning candidate, not the active transfinite vertex method; see §4. |
| `src/lib.rs:4–50`; `crates/ogeom-fillet/Cargo.toml:29`; `tools/old-path-{rounds,chamfers}.sh:1` | Module declarations, logger feature, migration checks | Remove declarations only for deleted files. Remove the feature, `to_current_path` and scripts last, after recording final zeros and the static no-old-constructor proof. |

### Shared-helper extraction checklist

These are concrete imports that make a file-name-only deletion fail to compile
or silently discard production behavior. Suggested destinations are ownership
boundaries, not a request to add compatibility shims.

| Helper in an old file | Retained caller | Survives / moves to |
|---|---|---|
| `fillet::{BallSeat(:5315), ball_seat(:5334), seat_centre(:5272), section_at(:5293)}` | `round/section.rs:23`; `asym.rs:95,114` | `round/section` or shared section math. Keep circular section arithmetic and its tolerances. |
| `setback::point_of(:90)`, `trim_conic(:852)` | `engine/classify.rs:13`, `engine/mod.rs:800`, `engine/setbacks.rs:868` | Shared topology/trim utilities or their engine owners. Keep world-placement handling. |
| `chamfer::{Chamfer, Legs, legs(:1690), wall_depth(:1426), rim_bevel_surface(:1939), same_spelling(:2419)}` | `engine/request.rs:6`, `fallback.rs:1222`, `trimmed.rs:47,8352`, `seams.rs:20` | Public request types remain; geometric readers move to engine seats/sections. `legs` also needs the relevant private seat/angle helpers. |
| `chamfer_sphere::{SphereSeat, sphere_seat(:175), sphere_bevel_rings(:346)}` | `chamfer::legs` and `engine/trimmed.rs:48` | Engine spherical-seat/section support; retain the transitive surface helpers needed by rings. |
| `chamfer_marched::crease_contacts(:485)` | `engine/trimmed.rs:2081` | Engine contact construction; keep its transitive section/station math, not old tool application. |
| `marched::band_fit_target(:3542)` and fitting helpers | `engine/{rims,rounds,round_strips,setbacks,trimmed}`, `round/transfinite.rs:1008`; retained face/sheet/fill paths | Shared fitting policy. Do not loosen fit targets during extraction. |
| `asym::{RoundProfile, legs(:133), conic_arc(:279), conic_pieces(:298)}`; `asym_law::varying_net(:323)` | `round/section.rs:22,177,360`, engine request/rim sections | Public profile type plus shared section geometry. |
| `support::{Seat, RevolvedSeat, planar_seat(:215), revolved_seat(:560)}` and geometry/guard utilities | `engine/{request,rounds,rims,trimmed}` and `round/section` | Shared support geometry, preserving checks as well as constructors. |
| `fillet::{whole_rims(:5669), rim_circle(:5652), same_circle(:5730)}` | `fillet::rounds_by_engine_routes(:1463)` before engine rim/curved assembly | Keep/move engine rim-chain normalization; it is not an obsolete band constructor. |
| `fillet::request_fits(:1831)`, `support::{ball_fits(:2016), ball_fits_by(:2047), bands_clear(:2529), Blamed}` | `blame.rs:68,140`, `fits.rs:508`, marched bridge | Engine/shared guard owner; preserve edge attribution, fit and overlap refusals. This can still call old shelf/rolled readers until ported. |
| `fillet::seated_fillet(:5393)` and its section/topology helpers | Public planar `blend_faces`, `facepair.rs:297`; vertex tools | Keep/extract until face-pair and vertex construction use replacements. Removing only edge fallback callers is insufficient. |
| `fillet::{marched_wedges, rim_ring_taken_whole, marched_strip_taken_whole, apply_set_aside}` plus `Mate(:342)`, junction/collection dependencies | `engine/marched_seats.rs:36–83`; `marched` takes `Mate` | Replace construction first; then remove the obsolete graph. This is not a small pure-helper move. |
| `fill_n::{FillBoundary, FillSide, Region}` and filling solver | Public filling API; the fixed-fill prototype `round/patch.rs:12` also calls `make_filling_n_fixed(:269)` | Keep public filling and its solver. `make_filling_n_fixed` has only the unused prototype as caller; it is not required by the active engine. |
| `bridge::compatible(:662)` | `engine/trimmed.rs:1928,1960,2135` | Keep; it is not old-path-specific. |

## 3. Preserve every regression introduced by #981

Discovery: `gh pr view 981 --json body,files,mergeCommit`, then
`git log origin/main --oneline --extended-regexp --grep='981|969|970|971|972|973|974|975|976'`
and `git show ddd31ffcf8b39c01e2a4e76f3705ff381d24ac5f`. The merged PR
contains **eleven test functions**: eight fillet integration tests, one fillet
unit test, one boolean integration test and one boolean unit test. The two
inline unit tests must not be missed by looking only in `contact_regressions.rs`.

All eleven pass at the pinned revision. A focused run of the eight fillet
integration tests with `old-path-log` created **no fallback log**. The same
public/API bypass caveat applies: the vertex-only half of `corner_history`
uses old construction even though it produces no log record.

| Issue / test at current location | Route today; old-specific assertions? | Port plan: preserve the regression, never delete it |
|---|---|---|
| #969: `tests/fillet/corner_history.rs:34`, `a_vertex_generates_its_corner_faces_through_composed_history` | Mixed: calls old public `round_vertex` at `:53` for all eight vertices; incident-edge and all-edge `fillet_edges` calls at `:68,72` use engine. `History::generated`, exact unit sphere/centre are **public geometry/ancestry**, not obsolete internals. | Keep both vertex-only and edge-request scenarios. Move the vertex-only API to its replacement when available; retain one sphere owned by the correct source vertex before/after composition, never replace it with “some generated face”. |
| #969: `src/corner.rs:1781`, `generated_corner_faces_belong_to_their_own_vertices` | Calls old private `round_vertex_with` and `round_vertices_together`, then engine `fillet_edges`. **Yes:** exactly four generated faces before band composition (sphere plus three planar caps), batching entry point and private chart/tool mechanics. | Move its geometric/ancestry coverage to the owning public vertex API test. Preserve eight distinct analytical centres/radii and final composed ownership. Replace the private four-face/tool-layout assertion with the public vertex-only cap/history contract; keep an explicit multiple-vertex composition scenario, not only the final all-edge result. No `pub` widening for tests. |
| #970: `crates/ogeom-bool/tests/bool/spline_extrusion.rs:10`, `a_coplanar_spline_cut_keeps_its_floor` | Boolean-only, no blend engine or old fillet call. No old-blend internals. | Keep unchanged: closed valid cut and independently converged volume **28840.851225099 mm³**. Its coplanar floor failure remains relevant after M6. |
| #971: `tests/fillet/extrusion_cap.rs:10`, `a_spline_cut_preserves_the_rounded_prisms_exact_cap` | Engine round of the quadrilateral prism, then boolean spline cut. No old-specific assertions. | Keep unchanged: cap vertices **and incident geometric curve endpoints** on **z=24.5** within model confusion. Checking only vertices would lose the 3.44303217e-6 mm endpoint regression. |
| #972: `tests/fillet/contact_regressions.rs:10`, `cylinder_contacts_leave_the_middle_wall_between_both_bands` | Engine ruled/meeting round route; no fallback records. Exact axes, positive middle-wall contact gap and volumes are independent geometry. | Keep both fixtures and closed-form volumes **1577.0636662254801 / 2610.337927095673 mm³**. Do not replace gap/tangency checks with validity or face counts. |
| Contact-radius support for #972/#974: `crates/ogeom-bool/src/lib.rs:17520`, `a_hug_measures_the_curve_rather_than_its_sampling_chords` | Boolean-only unit test of `hug_distance(:11227)`, not a fillet old path. Private helper, but **not** one removed by M6. | Keep unchanged: exact radial offsets on periodic circle ranges, and nearest endpoint on a bounded arc. If boolean distance ownership changes separately, port this oracle to that owner; M6 does not delete it. |
| #973: `tests/fillet/carried_rim.rs:55`, `a_stadium_foot_round_rolls_on_the_plate_edge` | Engine far-edge/support-transfer route; no fallback records. Exact envelope and quadrature, not `rolled_rim` internals. Marked heavy. | Keep unchanged, including `--include-ignored`: converged divergence-theorem volume at **1e-5 mm³**, and exact envelope probes through support switches on both stadium ends. |
| #974 positive + cavity negative: `tests/fillet/contact_regressions.rs:134`, `narrow_cylinder_shelf` | Engine far-edge route for the positive; engine void guard for the cavity negative. No old-path records. Error substring “a round reaches into a void” is a meaningful named refusal, not a private function name. | Keep all **187 exact envelope samples**, trimmed-face admission near transition endpoints and cavity-refusal check. Never compare against an old fitted shelf result or broaden to any error. |
| #974 bounded-rim negative: `tests/fillet/rim_bounds.rs:11`, `a_shelf_round_refuses_contact_beyond_the_circular_edges_end` | Engine refusal; no fallback record. “roll past a vertex” describes the unsupported bounded-support change. | Keep trimmed-rim fixture and named refusal. If that support change becomes supported, replace refusal only with an independently proved correct result, including the free plane-plane round near **y=2**, not an extrapolated circle. |
| #975: `tests/fillet/cap_history.rs:10`, `a_ruled_round_keeps_one_inherited_cap_at_each_end` | Engine ruled-seat route. Exactly two connected planar caps with inherited ancestry, not generated by the edge, is the regression's observable contract. | Keep unchanged. Preserve source-face tracing and exclusion from edge-generated history; merging/splitting implementation faces must not reclassify inherited caps as new blend faces. |
| #976: `tests/fillet/contact_regressions.rs:255`, `oblique_corner_mesh` | Engine planar round, then meshing; no old-path records. No private assertions. | Keep all **seven rotations**, every emitted mesh node within model confusion of the **trimmed BREP**, and closed mesh. Do not check only the supporting sphere or filter finished triangles to make it pass. |

No test is ported in this docs-only change. Most need no mechanical change;
the two #969 tests require care around vertex-only construction and the inline
private seam. The boolean and mesh fixes introduced by #981 are outside the
old blend deletion and must stay.

## 4. Tests coupled to old construction

An integration test using an old dispatcher is **not** an old-internals-only
test. Every Appendix A/B scenario stays and migrates. Filenames such as
`chain.rs`, `consumed_marched.rs`, `mixed_vertex.rs`, `setback.rs` and
`asymmetric.rs` do not authorize deletion. Update capability-ledger citations
if retained proof moves (`docs/capabilities.toml`, `docs/CAPABILITIES.md`).

### Rewrite the mechanism; keep the geometric or safety contract

| Test / location | Coupling and required rewrite |
|---|---|
| `support::tests::blend_wedges_walk_each_edge_once_each_way`, `src/support.rs:3104` | Forces `fillet_edge_wedge` and old mitred chamfer, then inspects collected wedges. Preserve independent spandrel/Pappus volumes and edge-use orientation using engine tools/result shells. The obsolete collection assertion can go only with its private seam. |
| `straight_seat_wedges_hold_their_pcurves_within_their_edges`, `src/support.rs:3162` | Replace forced old wedge with engine straight tools/results. Keep measured pcurve gap no greater than each edge's stated tolerance, for convex and concave fixtures. |
| `rim_arc_wedges_walk_each_edge_once_each_way`, `src/support.rs:3322`; `marched_wedges_walk_each_edge_once_each_way`, `:3601` | Use engine arc/marched tools or public result, preserving opposite edge uses, closedness and geometric volume. Do not replace with source-import assertions. |
| `a_pole_leg_meshes_to_its_area_as_the_bore_turns`, `src/support.rs:3421`; `a_pole_leg_blends_with_the_rim_s_middle_nearest_the_pole`, `:3439` | Shared fixture `pole_leg_blend(:3452)` collects old wedges. Port to engine pole-crossing output/tool; retain leg-area, mesh closure, volume and rotation/near-pole checks. |
| `a_band_round_crossing_bores_meshes_to_its_volume`, `src/support.rs:3535` | Replace old collected band; retain independently measured BREP/mesh volume and pole/seam coverage. |
| `a_band_carried_on_past_its_cap_checks_valid`, `src/junction.rs:2311` | Direct old wedge and `carried_into` call. Port the cap/run-out fixture to an engine meeting-tool or public chain request. Add an independent cap/volume property rather than keeping validity as its only proof. |
| `the_certified_check_agrees_with_the_fine_mesh`, `src/setback.rs:2943` | Old `SignedCorner`, `certify_within` and old corner fixtures; it logs two old requests. Preserve the safety property in the engine's patch-region guard with the fine-mesh cross-check **and negative strict margin**. Remove the old certifier test only in the same change that retains equivalent proof at the guard owner. |
| `generated_corner_faces_belong_to_their_own_vertices`, `src/corner.rs:1781` | Protected #981 regression: port as specified above; never delete its ancestry meaning. |
| `an_oblique_corner_rounds_the_same_under_every_labelling`, `src/corner.rs:1892`; `a_square_pyramid_apex_rounds_the_same_under_every_labelling`, `:1907` | Private forced-chart index. Preserve oblique/apex geometry and order/orientation invariance through public vertex API and analytical/independent volume oracle; drop forced internal chart enumeration. |
| `a_cone_s_rim_is_the_torus_the_current_path_cuts`, `src/engine/rims.rs:1020` | Despite its name, `:1029` calls engine-first public `fillet_edge`, producing no old log today: the “old” oracle is no longer independent. Rewrite expected removed volume from torus/section-centroid integration. Keep cone angle and rim geometry. |
| `a_seated_wedge_is_cut_by_trimming_only_where_nothing_reaches_into_it`, `src/support.rs:3225` | Tests **boolean** `cut_seated`/general-cut equivalence and work counters, but gets its tool through `fillet_edge_wedge`. Move ownership to `ogeom-bool` tests with an independent analytic tool, or use an engine tool without test-only exports. Keep drilled and undrilled cases and refusal of unsafe seated cutting. |

### Assertions that may disappear with an obsolete optimization

`corners_clear_of_each_other_round_in_one_cut_as_in_turn`
(`src/corner.rs:1843`) and `corners_whose_balls_come_close_round_in_turn`
(`:1880`) specifically test `round_vertices_together`, gathered-face work
counts, and its `None` choice near touching balls. These batching-selection
and work-count assertions can go when that optimization and its seam go.
**First retain** the public many-vertex/disjoint/near-touching geometry and
history cases, with an independent volume/region oracle; do not copy an
obsolete boolean-count target into an engine test. They are not an excuse to
remove #981's distinct per-vertex ownership checks.

The fixed-fill prototype test
`the_fixed_patch_through_three_equal_strips_stands_between_ball_and_corner`
(`src/round/patch.rs:133`) has no production caller of its `vertex_patch(:52)`
owner: the only call is in this test at `:177`. It tests the prototype's
least-bending departure from the exact ball, not the active engine's
transfinite method. It and `make_filling_n_fixed` may be pruned separately
with that unused prototype after checking documentation citations. Keep the
active engine's independent corner-region and tangency proofs; do not retain
this prototype solely to keep its test running. This does not authorize
removing public `make_filling_n` or its solver.

### Keep: shared algorithms, not deletion candidates

Keep `between::tests::a_run_clears_the_faces_beside_it_whichever_way_its_stations_run`
(`src/between.rs:577`),
`pair_marched::tests::a_line_of_contact_holds_its_stated_error_between_the_stations`
(`src/pair_marched.rs:2903`),
`support::tests::a_cone_s_farthest_exit_over_a_triangle_is_its_hull_exit`
(`src/support.rs:3013`), the two chart/continuation tests at
`src/marched.rs:5727,5771`, the two solver tests at
`src/fill_patch.rs:506,542`, and the law-window test at `src/law.rs:581`.
Their owners have retained face/sheet/fill/guard callers. Likewise keep the
engine and `round/*` tests. The rollback test
`chamfer::tests::a_request_the_engine_does_not_build_leaves_the_model_as_it_was`
(`src/chamfer.rs:2495`) already tests an engine attempt, not old wedges; move
it with the adapter and retain model-unchanged-on-decline behavior.

## 5. Landing sequence and completion gates

Estimates are engineering time excluding review queues and machine-dependent
cold builds. Stages 2–4 may need several family-sized PRs; do not merge a
partially green migration. Each PR keeps the workspace green and retains the
currently correct result/refusal contracts. Refresh this manifest against main
before editing; other routing work is active.

1. **Make the evidence trustworthy — one tooling PR, 0.5–1 day.** Land the
   independent counter reliability commit from #1092 (do not duplicate it).
   Preserve Cargo failure status/diagnostics and atomic records; include a
   failing-Cargo negative fixture. Extend round measurement to facade tests
   and identify direct/unlogged callers. Gate: full counts with successful
   Cargo exit statuses; no `unknown`/concatenated asking names.
2. **Finish the already-owned routing families — family PRs, 1–3 days each.**
   Integrate #1092 and other active migrations first. Remaining families are
   variable curved laws, asymmetric/conic curves and corners, carried and
   re-entrant chains/rims, signed/mixed/unequal corners, and named fit/overlap
   refusals (Appendix A). Treat the 27 facade asking tests in Appendix B as
   blockers too. Keep edge-order permutations, history, geometric oracles
   and the ADR's real-world Voron 2.4/DoomCube check at each routing change.
   Gate: each family's dispatch count reaches zero with no lost acceptance;
   final complete round and chamfer counts are zero, not just selected suites.
3. **Port unlogged APIs and constructors — small owner PRs, 1–3 days each.**
   Route linear and asymmetric-variable APIs. Replace the marched bridge's
   old constructors without moving the case catalogue intact. Port vertex-only,
   setback and sharp APIs, preserving caps, isolation and history semantics;
   retain face/sheet APIs. The marched bridge and vertex-only work are real
   implementation prerequisites, not a half-day deletion exercise. Port the
   inline #981 test and remaining private-seam tests as separate test changes
   coupled to their available engine owner. Gate: all eleven #981 tests pass;
   no engine route or retained API still reaches an old constructor.
4. **Extract shared survivors, then retire chamfer wedges — 2–3 small PRs,
   0.5–1 day each.** Move section/seat/fit/trim helpers by ownership without
   changing tolerance or geometry. Keep public enums, signatures, ownership,
   cancellation, blame and refusal behavior. Then remove the five chamfer
   entry families' fallback tails and the old facet/wedge graph, keeping shared
   spherical and marched contacts. Gate: all chamfer and mixed suites plus facade checks,
   pcurve/history proofs and strict Clippy pass; static old chamfer callers are
   absent, including formerly direct wedge tests.
5. **Execute the prepared round deletion — one mechanical PR, 1–2 days,
   after stages 1–4.** Remove the six logged round tails, old variable builders,
   dead corner/junction/consumed/rolled catalogue files and their `mod`
   declarations as permitted by the manifest. Remove obsolete tests only with
   the proof/ports above, not by mass deletion of integration files. Keep the
   engine fallback and all shared/public survivors. Run final zero counters
   before removing the logger feature/scripts in this same final retirement
   change. Gate: no residual old-constructor callers, no dead-code lint
   suppression added to keep an orphan graph compiling, full repeated workspace
   checks, capability-ledger validation and the real-world corpus pass.

Chamfer retirement can precede completion of every round family if its shared
survivors and mixed/asymmetric callers have been separated. It cannot delete
`chamfer_marched`, `chamfer_sphere` or `support` by filename. The final round
removal is the “one clean pass”; the preceding PRs make its prerequisites
explicit rather than bundling new geometry with deletion.

### Verification commands for the removal owner

Start with plain Cargo and `~/.local/share/mbx/bin` first in `PATH`; never set
`CARGO_TARGET_DIR`. Give each Cargo/build command a timeout of at least 900 s.
If mbx refuses a write with “Read-only file system”, use direct Cargo and the
worktree's own `target/`, then remove that worktree build directory at the end.
Do not remove a shared mbx cache through a `target` symlink. In this research
run the mbx wrapper timed out before producing diagnostics; direct
`/usr/bin/cargo` completed. No read-only refusal is claimed for that timeout.

```sh
export PATH="$HOME/.local/share/mbx/bin:$PATH"
export RUST_TEST_THREADS=1
# After the counter reliability change: a nonzero Cargo status must fail these.
tools/old-path-rounds.sh -- --include-ignored
tools/old-path-chamfers.sh -- --include-ignored
cargo test -p ogeom-fillet --features old-path-log -- --include-ignored
cargo test -p ogeom --features ogeom-fillet/old-path-log -- --include-ignored
cargo test -p ogeom-bool --test bool -- --include-ignored spline_extrusion::
cargo test -p ogeom-bool --lib a_hug_measures_the_curve_rather_than_its_sampling_chords
cargo fmt --all -- --check
cargo clippy -p ogeom-fillet --all-targets --all-features -- -D warnings
node tools/lint-comment-rot.mjs --all
./tools/check.sh
```

At the pinned baseline, retain independent Cargo statuses around the scripts;
their exit status alone cannot prove success. Before final deletion also run
`rg` over the whole workspace for the constructor names in §1, including
`fillet_edge_wedge`, `marched_wedges`, `apply_set_aside`, `round_vertex_with`,
`chamfer_edges_mitred` and the module manifest. Resolve each non-test caller;
resolve each test by retained contract, not by adding an allow-list. Inspect
multiline imports, re-exports and qualified calls, not just `use crate::...`.

Research validation is deliberately narrower than the final removal gate:
the two complete counter scripts passed their underlying fillet/facade Cargo
runs; focused #981 fillet integration, inline corner-history and boolean
regressions passed. No production code, test, tolerance, workflow or manifest
is changed by this plan. The whole repeated workspace gate and downstream
real-world corpus are requirements for the implementation PRs, not claims
made by this docs-only research.

## Appendix A: all 77 asking tests in the fillet count

Each row gives observed requests and logged entry names, not a prediction of
which helper branch is used. Every scenario stays; only the old mechanism
changes. These are also the test callers for the round entry-point table.

| Test location and name | Requests | Logged entries |
|---|---:|---|
| `src/setback.rs:2943` — `setback::tests::the_certified_check_agrees_with_the_fine_mesh` | 2 | fillet_edges × 2 |
| `tests/fillet/asymmetric.rs:264` — `a_slanted_edge_takes_the_ellipse_of_its_own_dihedral` | 2 | fillet_edge_distances × 2 |
| `tests/fillet/asymmetric.rs:456` — `what_is_not_built_is_refused_by_name` | 3 | fillet_edge_distances × 3 |
| `tests/fillet/asymmetric.rs:548` — `blend_edges_takes_an_asymmetric_round_beside_others_and_refuses_one_meeting_a_bevel_of_the_other_convexity` | 1 | blend_edges × 1 |
| `tests/fillet/asymmetric.rs:747` — `an_asymmetric_round_runs_out_along_its_wall_onto_a_ramp` | 8 | fillet_edge_distances × 8 |
| `tests/fillet/asymmetric.rs:795` — `an_asymmetric_round_running_up_under_a_roof_stops_on_it` | 8 | fillet_edge_distances × 8 |
| `tests/fillet/asymmetric.rs:991` — `a_bore_s_rim_takes_its_conic_corner_round_the_axis` | 4 | fillet_edge_distances × 4 |
| `tests/fillet/asymmetric.rs:1209` — `what_a_rim_does_not_build_is_refused_by_name` | 3 | fillet_edge_distances × 3 |
| `tests/fillet/asymmetric.rs:1268` — `blend_edges_takes_an_asymmetric_rim_beside_a_round_rim` | 1 | blend_edges × 1 |
| `tests/fillet/asymmetric.rs:1629` — `a_rim_arc_takes_its_conic_corner_over_its_own_turn` | 4 | fillet_edge_distances × 4 |
| `tests/fillet/asymmetric.rs:1689` — `a_tangent_loop_of_asymmetric_rounds_closes_round_a_pad` | 2 | blend_edges × 2 |
| `tests/fillet/asymmetric.rs:1758` — `a_tangent_loop_of_asymmetric_rounds_fills_round_a_pocket_s_floor` | 2 | blend_edges × 2 |
| `tests/fillet/asymmetric.rs:1838` — `asymmetric_rounds_meeting_other_than_as_a_tangent_chain_are_refused` | 2 | blend_edges × 2 |
| `tests/fillet/asymmetric.rs:2047` — `an_open_spline_crease_takes_its_conic_corner_along_its_length` | 4 | fillet_edge_distances × 4 |
| `tests/fillet/asymmetric.rs:2148` — `a_bore_s_rim_on_a_drum_takes_a_marched_conic_band` | 4 | fillet_edge_distances × 4 |
| `tests/fillet/asymmetric.rs:2253` — `asymmetric_rounds_close_a_box_corner_with_rounds_and_bevels` | 6 | blend_edges × 3, fillet_edge_distances × 3 |
| `tests/fillet/asymmetric.rs:2438` — `the_largest_pair_is_where_the_first_contact_leaves_its_face` | 4 | fillet_edge_distances × 4 |
| `tests/fillet/blame.rs:54` — `a_plate_s_overlap_names_an_edge_and_the_one_across_the_wall` | 1 | fillet_edges × 1 |
| `tests/fillet/blame.rs:194` — `a_notched_box_s_refusal_names_the_notch_edge` | 1 | fillet_edges × 1 |
| `tests/fillet/cell_guard.rs:294` — `a_convex_and_a_concave_round_with_a_void_in_their_cell_are_refused_by_name` | 3 | fillet_edges × 3 |
| `tests/fillet/chain.rs:545` — `a_short_arc_of_a_rim_that_stops_rounds_over_its_own_arc` | 2 | fillet_edges × 2 |
| `tests/fillet/chain.rs:785` — `blends_at_an_inside_corner_whose_vertex_stands_off_the_wall_meet_at_the_wall` | 2 | fillet_edges × 2 |
| `tests/fillet/chain.rs:1570` — `rounds_at_inside_corners_off_square_or_leaning_both_ways_take_all_above_both` | 4 | fillet_edges × 4 |
| `tests/fillet/chain.rs:1915` — `rounds_at_a_curved_crease_turning_nearly_back_on_itself` | 4 | fillet_edges × 4 |
| `tests/fillet/chain.rs:2650` — `a_cone_rim_that_meets_unresolved_is_rounded_fitted` | 4 | fillet_edge × 4 |
| `tests/fillet/chain.rs:2705` — `rounds_at_a_bore_meeting_a_drafted_chord_take_all_above_both` | 2 | fillet_edges × 2 |
| `tests/fillet/chain.rs:2739` — `a_pocket_of_a_disc_and_its_chord_meets_at_both_corners` | 1 | fillet_edges × 1 |
| `tests/fillet/chain.rs:2840` — `rounds_of_a_d_shaped_block_take_both_bands_as_alone` | 4 | fillet_edges × 4 |
| `tests/fillet/chain.rs:2911` — `rounds_of_a_d_shaped_block_with_an_overhanging_chord_close_in_both_orders` | 4 | fillet_edges × 4 |
| `tests/fillet/chain.rs:2939` — `rounds_of_a_d_shaped_block_on_a_drafted_drum_take_both_bands_as_alone` | 12 | fillet_edges × 12 |
| `tests/fillet/chain.rs:3089` — `rims_are_cut_by_trimming_to_the_general_fuse_s_solid` | 2 | fillet_edges × 2 |
| `tests/fillet/chain.rs:3430` — `rims_grouped_with_other_wedges_are_trimmed_before_the_rest_is_cut` | 1 | fillet_edges × 1 |
| `tests/fillet/chain.rs:3907` — `bands_running_into_a_rounded_corner_are_trimmed_to_the_general_fuse_s_solid` | 1 | fillet_edges × 1 |
| `tests/fillet/chain.rs:4325` — `rims_pieces_in_a_tangent_chain_are_trimmed_to_the_general_fuse_s_solid` | 1 | fillet_edges × 1 |
| `tests/fillet/concave.rs:810` — `a_round_rolling_on_a_far_edge_refuses_where_the_ball_leaves_it` | 1 | fillet_edge × 1 |
| `tests/fillet/concave.rs:1361` — `a_round_rolls_round_the_corner_of_a_narrow_triangle` | 1 | fillet_edge × 1 |
| `tests/fillet/consumed_marched.rs:127` — `a_bevel_wider_than_the_round_s_setback_is_still_refused` | 1 | fillet_edges × 1 |
| `tests/fillet/fillet.rs:418` — `a_rim_fillet_that_swallows_the_axis_is_refused` | 1 | fillet_edge × 1 |
| `tests/fillet/fillet.rs:623` — `a_fillet_the_faces_cannot_hold_is_refused` | 2 | fillet_edge × 2 |
| `tests/fillet/fillet.rs:638` — `blends_overlapping_across_a_shared_face_are_refused` | 1 | fillet_edges × 1 |
| `tests/fillet/fillet.rs:1119` — `seams_whose_ball_runs_off_a_narrow_face_round_beside_other_edges` | 10 | fillet_edges × 10 |
| `tests/fillet/fillet.rs:1359` — `a_seam_cut_off_at_a_side_inside_another_band_is_refused` | 1 | fillet_edges × 1 |
| `tests/fillet/fillet.rs:1426` — `every_edge_of_a_thin_plate_overlaps_alike_built_as_box_or_prism` | 4 | fillet_edges × 4 |
| `tests/fillet/fillet.rs:2100` — `the_top_loop_of_a_block_with_rounded_verticals_rounds_on_the_turns` | 1 | fillet_edges × 1 |
| `tests/fillet/fillet.rs:2238` — `a_rim_arc_ending_obliquely_on_a_flat_rounds_on_past_its_ends` | 3 | fillet_edge × 3 |
| `tests/fillet/largest.rs:63` — `a_plate_s_edge_rounds_at_most_its_thickness` | 1 | fillet_edges × 1 |
| `tests/fillet/largest.rs:131` — `a_plate_rounds_all_round_at_most_half_its_thickness` | 1 | fillet_edges × 1 |
| `tests/fillet/largest.rs:145` — `a_prism_built_plate_rounds_all_round_at_most_half_its_thickness` | 1 | fillet_edges × 1 |
| `tests/fillet/largest.rs:184` — `a_cone_rim_rounds_at_most_where_its_wall_ends` | 2 | fillet_edges × 2 |
| `tests/fillet/marched_fillet.rs:1156` — `rounds_of_a_rib_into_a_boss_take_the_union_of_each_alone` | 4 | fillet_edges × 4 |
| `tests/fillet/mixed_vertex.rs:87` — `a_convex_and_a_concave_round_beside_an_edge_left_sharp_close_one_corner_in_either_order` | 6 | fillet_edges × 6 |
| `tests/fillet/mixed_vertex.rs:194` — `three_rounds_of_both_kinds_at_a_vertex_close_one_corner_whatever_the_order` | 8 | fillet_edge × 2, fillet_edges × 6 |
| `tests/fillet/mixed_vertex.rs:311` — `a_signed_corner_whose_cell_reaches_through_the_wall_is_refused_by_name` | 2 | fillet_edges × 2 |
| `tests/fillet/radius_per_edge.rs:285` — `a_pyramid_apex_at_alternating_radii_closes_with_one_patch` | 2 | fillet_edges_with × 2 |
| `tests/fillet/tangent_chain.rs:140` — `a_foot_round_on_a_face_narrowing_to_a_cusp_is_cut_off_at_the_face_s_side` | 4 | fillet_edges × 4 |
| `tests/fillet/tangent_chain.rs:195` — `a_tangent_chain_through_the_cusp_rounds_either_side_of_it_in_either_order` | 2 | fillet_edges × 2 |
| `tests/fillet/tangent_chain.rs:654` — `the_arm_s_side_edge_and_its_arc_into_the_hub_turn_through_an_s_in_either_order` | 2 | fillet_edges × 2 |
| `tests/fillet/tangent_chain.rs:665` — `printcad_s_chain_asked_on_one_side_of_the_tangent_point_is_refused_by_name` | 2 | fillet_edges × 2 |
| `tests/fillet/variable.rs:230` — `a_law_of_points_rounds_a_drum_s_rim` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:294` — `a_law_of_points_rounds_a_marched_seat_tangent_to_both_faces` | 2 | fillet_edges_variable × 2 |
| `tests/fillet/variable.rs:419` — `a_law_of_points_rounds_an_open_arc_running_out_at_both_ends` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:481` — `a_law_of_points_rounds_an_open_chain_of_arcs_of_one_circle` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:559` — `a_law_of_points_rounds_a_flat_meeting_a_drum_along_its_ruling` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:648` — `a_law_of_points_rounds_a_seam_pinched_where_two_drums_turn_tangent` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:814` — `a_law_of_points_rounds_the_top_loop_of_a_rounded_rectangle_pad` | 2 | fillet_edges × 1, fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:1163` — `a_law_of_points_rounds_a_pad_s_top_loop_with_one_corner_left_sharp` | 2 | fillet_edges × 1, fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:1514` — `a_law_whose_ball_runs_off_a_narrow_face_is_cut_off_as_the_constant_round` | 6 | fillet_edge × 1, fillet_edges_variable × 5 |
| `tests/fillet/variable.rs:1577` — `a_law_whose_ball_runs_off_a_narrow_face_into_the_material_is_refused_by_name` | 4 | fillet_edge × 1, fillet_edges_variable × 3 |
| `tests/fillet/variable.rs:1630` — `a_law_takes_a_narrow_strip_beside_its_edge_off_whole` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:1718` — `a_law_whose_ball_sits_on_a_narrow_strip_rounds_the_edge_beside_it` | 2 | fillet_edges_variable × 2 |
| `tests/fillet/variable.rs:1780` — `a_law_of_points_rounds_a_drum_s_ruling_turning_into_its_rim` | 3 | fillet_edge × 2, fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:1845` — `a_variable_radius_refuses_by_name_what_it_does_not_round` | 6 | fillet_edges_variable × 6 |
| `tests/fillet/variable.rs:1905` — `a_law_steep_at_an_open_arc_s_end_caps_in_the_cut_it_ends_on` | 2 | fillet_edges_variable × 2 |
| `tests/fillet/variable.rs:2026` — `a_law_over_a_chain_beside_a_narrow_face_is_cut_off_across_a_tangent_junction` | 1 | fillet_edges_variable × 1 |
| `tests/fillet/variable.rs:2065` — `a_law_over_a_chain_beside_a_narrow_face_is_cut_off_across_a_convex_turn` | 3 | fillet_edges × 1, fillet_edges_variable × 2 |
| `tests/fillet/variable.rs:2185` — `a_law_round_a_marched_rim_beside_a_narrow_face_is_cut_off_as_the_constant_round` | 2 | fillet_edges_variable × 2 |
| `tests/fillet/variable.rs:2264` — `a_law_round_a_rim_crossing_a_bevel_s_width_is_refused_by_name` | 1 | fillet_edges_variable × 1 |

## Appendix B: the 27 additional facade asking tests

These 220 requests are visible in the chamfer script's log but excluded from
its chamfer count. They are not a separate API inventory: they enter through
the same public round functions. Direct vertex/face/sheet calls may make
additional unlogged old-constructor requests.

| Test location and name | Requests | Logged entries |
|---|---:|---|
| `crates/ogeom-fillet/tests/fillet/blend_reach.rs:73` — `a_fillet_wider_than_its_face_is_refused` | 1 | fillet_edges × 1 |
| `crates/ogeom-fillet/tests/fillet/blend_reach.rs:100` — `a_rim_fillet_wider_than_its_cap_is_refused` | 1 | fillet_edges × 1 |
| `crates/ogeom-fillet/tests/fillet/blend_reach.rs:128` — `a_concave_fillet_wider_than_its_floor_is_refused` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/facade/blends.rs:1916` — `an_edge_round_stops_where_an_upright_drum_s_foot_runs_off_the_block` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/facade/blends.rs:2099` — `an_upright_drum_s_face_blend_ends_where_the_crease_runs_off_the_block` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/facade/blends.rs:2292` — `b2_three_fillets_and_the_corner_tool_round_the_vertex` | 1 | fillet_edge × 1 |
| `crates/ogeom/tests/facade/blends.rs:2359` — `round_vertex_reproduces_the_b2_closed_form` | 1 | fillet_edge × 1 |
| `crates/ogeom/tests/facade/blends.rs:2410` — `round_vertex_rounds_the_corner_at_any_placement` | 10 | fillet_edge × 10 |
| `crates/ogeom/tests/facade/blends.rs:2485` — `round_vertex_rounds_an_oblique_corner` | 2 | fillet_edge × 2 |
| `crates/ogeom/tests/facade/blends.rs:2642` — `round_vertex_carries_a_band_stopped_by_a_run_out_back_to_the_corner` | 1 | fillet_edge × 1 |
| `crates/ogeom/tests/facade/blends.rs:2659` — `round_vertex_carries_two_bands_stopped_by_one_run_out` | 2 | fillet_edge × 2 |
| `crates/ogeom/tests/facade/blends.rs:2752` — `round_vertex_after_one_edge_at_a_time_is_the_corner_rounded_in_one_call` | 149 | fillet_edge × 149 |
| `crates/ogeom/tests/facade/blends.rs:3154` — `round_vertex_rounds_flat_and_oblique_apexes_no_ball_touches` | 10 | fillet_edge × 10 |
| `crates/ogeom/tests/facade/blends.rs:3361` — `a_curved_corner_in_one_call_has_one_ball_and_one_round_at_a_time_takes_the_union` | 3 | fillet_edge × 2, fillet_edges × 1 |
| `crates/ogeom/tests/facade/blends.rs:3408` — `a_curved_corner_closes_with_the_rim_s_band_first` | 5 | fillet_edge × 5 |
| `crates/ogeom/tests/facade/blends.rs:3507` — `a_curved_corner_s_straight_band_first_runs_out_across_the_drum` | 4 | fillet_edge × 4 |
| `crates/ogeom/tests/facade/blends.rs:3843` — `a_rim_blend_rolls_over_the_bracket_s_concave_blend` | 1 | fillet_edges × 1 |
| `crates/ogeom/tests/facade/blends.rs:3912` — `a_marched_blend_meets_a_straight_blend_at_its_corners` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/facade/blends.rs:4393` — `a_rim_whose_rail_passes_within_a_stride_of_the_pole_blends` | 3 | fillet_edge × 2, fillet_edges × 1 |
| `crates/ogeom/tests/facade/blends.rs:4474` — `blends_along_a_rounded_rim_turn_each_arc_by_its_quarter` | 1 | fillet_edges × 1 |
| `crates/ogeom/tests/facade/edge_chains.rs:331` — `edges_crossing_beside_a_drum_take_what_each_takes_alone` | 6 | fillet_edges × 6 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:1454` — `rounds_meeting_at_a_vertex_their_edge_overruns_are_capped_there` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:1501` — `rounds_meeting_at_a_vertex_off_one_edges_wall_are_capped_at_the_wall` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:1759` — `a_rib_s_chain_rounds_on_along_its_tangent_continuations` | 3 | fillet_edges × 3 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:1855` — `a_rib_s_lower_spline_crease_rounds_alone_and_with_its_chain` | 2 | fillet_edges × 2 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:2493` — `a_kidney_pocket_floor_loop_rounds_to_its_volume` | 1 | fillet_edges × 1 |
| `crates/ogeom/tests/corpus_suites/imported_corpus.rs:2517` — `a_kidney_pocket_floor_loop_rounds_in_time` | 1 | fillet_edges × 1 |
