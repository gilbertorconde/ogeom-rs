# Design note: fitting and mesh to solid out of `ogeom-algo` (#1264)

Tracked in #1272. Depends on #1260 (the explicit facade), which is merged.

## Verdict: partly valid

The problem is real and has grown. The proposed groups cannot be cut as
drawn: three of the "fitting" modules and `traced` are construction
internals, and the mesh-to-solid group has a user inside `ogeom-bool`.

Measured on `main` at `09979ec1`. Production lines only: inline
`#[cfg(test)] mod` blocks and `tests.rs` files are excluded. The issue's
numbers were measured at `456c3234`.

| | Issue | Now |
|---|---|---|
| `ogeom-algo`, production lines | ~78,000 | 88,707 |
| `.rs` files | 85 | 98 |
| mesh to solid (`mesh_solid`, `recognize*`, `register`, `traced`) | ~19,500 | 21,192 (`mesh_solid` alone 16,757) |
| fitting (`fit`, `held`, `network`, `pcurve_fit`, `cutoff`, `bins`, `length`) | ~12,300 | 15,580 |
| `ogeom_algo::` references in `ogeom-fillet/src` | 881 | 1,160 |

### Which crates use each group

Checked with `grep` over `crates/*/src` and `tools/`, test modules excluded.

- **Mesh to solid** (`solid_from_mesh`, `MeshSolidOptions`, ...):
  - `ogeom-bool`: `mesh_region.rs` (`cut_by_mesh_region` builds its slab with
    `solid_from_mesh`).
  - `ogeom-io`: `stl.rs`, only `single_precision_quantum`, a function of a
    `Triangulation` alone.
  - The facade, its tests, and `tools/ogeom-bench` and `ogeom-stress`.
- **Fitting**:
  - `ogeom-fillet`: `fit_network` and its types, `net_held_within`,
    `interpolate`, `Spacing`.
  - `ogeom-mesh`: tests only (`interpolate`, `Spacing`), through the existing
    path-only dev-dependency.
  - The facade.
- **Neither group** is used by `ogeom-heal`, `ogeom-hlr`, `ogeom-offset` or
  `ogeom-doc`. `traced` is used by `build`, `ogeom-bool`, `ogeom-heal` and
  `ogeom-io`.

### Are the groups separable?

The `crate::` references between top-level modules were read from source,
test modules excluded.

- **`pcurve_fit`, `length` and `bins` are construction internals.**
  `pcurve_fit` is used by `build` (7), `sew`, `fill_n`, `upgrade` and
  `pcurve_gap`. `length` is used by `check`, `extend` and `frames`. `bins` is
  used by `sew`. They stay in `ogeom-algo`.
- **`traced` is construction.** It is used by `build` and three other crates.
  It stays.
- **`recognize`, `orient` and `cutoff` belong with fitting.**
  - `fit/primitive.rs` uses `recognize::{centroid, covariance, fit_sphere,
    circle_through}`.
  - `fit/grow.rs` uses `orient::{check_indices, edge_uses, key}`.
  - Every fit uses `Cutoff`.
  - `recognize` itself only uses `measure`.
- **`outline`, `segments_2d` and `rest` use `fit::segments`.**
  `rest` also uses `orient::check_indices` and the private `hull3`, which
  nothing else uses. All four go with fitting.
- **Mesh to solid reaches into fitting** (`fit` 23 references, `recognize`
  22), so its crate sits above the fitting crate.
- **The core never reaches up into either group.**
  - The only reference from a module that stays to a module that moves is
    `classify.rs:3174`, in its test module. It builds a solid through
    `solid_from_mesh`.
  - From the moved modules into the core, seven `pub(crate)` items are used.
    They must become reachable (listed below). Everything else they use is
    already public: `build`, `check`, `mass`, `measure`, `pcurve_gap`,
    `pcurve_fit`, `primitive`, `tight`, `chart_cross`, `self_cross`.

### What the normative documents say

- `docs/DATA_MODEL.md` says nothing about crate boundaries. Line 210 lists the
  crates whose operations fill in history, and the new crates join that list.
- `docs/SCOPE.md` admits mesh to solid as "the one addition". It ends: "Never
  through a pull request that quietly adds a crate." So the first PR that adds
  a crate records it there.
- `docs/PLAN.md` ("Structure") already says: "Mesh conversion moved to its own
  crate would stop it rebuilding the whole graph."
- `CONTRIBUTING.md` says "users reach them only through [the facade]" and
  "An item that moves to another crate keeps its `ogeom::<module>` path".
  Under "Releasing" it counts fifteen crates published together, and records
  the `ogeom-mesh` → `ogeom-algo` path-only dev-dependency.

## Decision

### Target layout

Two new library crates above `ogeom-algo`, published with the others at the
workspace version. Module names do not change, so every facade module path
(`ogeom::algo::fit::zone`, `ogeom::algo::register::Icp`, ...) resolves
unchanged.

**`ogeom-fit`**: fitting and reading of measured points and meshes. About
20,800 lines.

| Module | Lines | Notes |
|---|---|---|
| `fit` | 8,043 | primitive fits, zones, region growing, `segments`, `tails` |
| `network` | 5,775 | surfaces fitted to curve networks |
| `recognize` | 1,555 | |
| `orient` | 1,311 | |
| `register` | 1,100 | |
| `outline` | 1,095 | |
| `segments_2d` | 864 | |
| `rest` | 442 | |
| `hull3` | 252 | private |
| `held` | 220 | |
| `cutoff` | 109 | private |

**`ogeom-recon`**: a B-rep solid from a triangle mesh. About 18,800 lines.

| Module | Lines | Notes |
|---|---|---|
| `mesh_solid` | 16,757 | with its 15 submodules |
| `recognize_patch` | 957 | private |
| `recognize_swept` | 560 | private |
| `mesh_region` | 491 | from `ogeom-bool` |

**`ogeom-algo` keeps** construction, sweeps and sewing, mass, check and
classify, and medial axis and distance: about 49,700 lines.

Two small pieces move down instead of up:

- **Curve fitting through points goes to `ogeom_geom::fit`.** That is
  `Spacing`, `interpolate`, `approximate` and `approximate_within`, plus the
  crate-internal `spaced` and `interpolate_at` as `#[doc(hidden)] pub`: about
  450 lines from `fit.rs`. They use only `ogeom-core`, `ogeom-math`,
  `ogeom-geom` and `nalgebra`. They wrap `ogeom_geom::fit::fit_points`, which
  their own docs already say "lives in `ogeom_geom::fit`". This keeps
  `ogeom-fillet` and the `ogeom-mesh` tests off the fitting crate for plain
  curve fitting, and adds no second dev-dependency cycle.
- **`single_precision_quantum` goes to `ogeom-mesh`.** It is a function of a
  `Triangulation`. `ogeom-io` then needs nothing from mesh to solid.

The crate names follow the issue. `ogeom-recon` is distinct from
`outside/crates/ogeom-recognize`, which is feature recognition and out of
scope.

### Public items that move

Each move keeps its facade path. "Module" means the `pub mod` path moves with
all its items.

**To `ogeom_fit`:**

- Module `fit`, with `AnyFit`, `Crease`, `FitKind`, `GrowOptions`, `Grown`,
  `Meeting`, `SeedReport`, `grow`, `Fit`, `FitFixed`, `FitOptions`, `Side`,
  `circle2`, `circle3`, `cone`, `cylinder`, `line2`, `line3`, `plane`,
  `sphere`, `torus`, and module `zone`. `ogeom_fit::fit` re-exports `Spacing`,
  `interpolate`, `approximate` and `approximate_within` from `ogeom_geom::fit`,
  so `ogeom::algo::fit::interpolate` still resolves.
- Module `network`, with `Across`, `ChainedPatch`, `CurveDef`,
  `IrregularCorner`, `LoopSide`, `Network`, `NetworkOptions`, `PatchDef`,
  `PatchFit`, `PatchSide`, `fit_loop_quads`, `fit_network`,
  `fit_network_chained`.
- Module `held`, with `NetOutside` and `net_held_within`.
- `Cutoff`.
- Module `recognize`, with `Canonical`, `Recognized` and `recognize_points`.
- Module `orient` (`Directions`, `FaceDirection`, `Revolution`, `Rim`,
  `face_directions`, `axis_of_revolution`, `principal_frame`, `open_rims`).
- Module `register` (`Icp`, `Starts` and the rest of its items).
- Module `outline`, with `Outline`, `Piece2`, `Wire2d`, `fit_outline`.
- Module `segments_2d`, with `Segment2`, `Segments`, `fit_segments_2d`.
- Module `rest`, with `Rest`, `resting_faces`, `resting_faces_of`.

**To `ogeom_recon`:**

- Module `mesh_solid`, with `FallbackReason`, `FitConstraints`, `MeshRegion`,
  `MeshRegions`, `MeshSolid`, `MeshSolidOptions`, `MeshSolidReport`,
  `RegionFallback`, `RegionId`, `RegionRefusal`, `SurfaceKind`,
  `refine_solid`, `solid_from_mesh`.
- `MeshRegions::find_as_listed`, behind a `test-hooks` feature of
  `ogeom-recon`, as it is behind `ogeom-algo`'s today.
- `ogeom_bool::cut_by_mesh_region` becomes `ogeom_recon::cut_by_mesh_region`.

**Down:**

- `ogeom_algo::{Spacing, interpolate, approximate, approximate_within}` →
  `ogeom_geom::fit::*`.
- `ogeom_algo::single_precision_quantum` → `ogeom_mesh::single_precision_quantum`.
  `ogeom_recon` re-exports it beside `MeshSolidOptions::quantum`.

**Made reachable across the boundary.** These are `pub(crate)` today and become
`#[doc(hidden)] pub`. Each is listed in `tools/facade.py` `INTERNAL`, the
pattern `FillRegion` and `merge_edges_at` already use.

- `check::check_wire` and `check::edge_spanned_past` (used by `mesh_solid`)
- `mass::volume_as_flagged` and `mass::VolumeKept` (used by `mesh_solid`)
- `self_cross::crossings` (used by `mesh_solid`)
- `distance::to_mesh_with` and `distance::Pseudonormals`, with
  `Pseudonormals::of` (used by `register`)

### Dependency graph after

Internal crates only. External dependencies (`nalgebra`, `smallvec`, `spade`)
go to whichever crate uses them.

```
core ← math ← geom ← topo
geom ← intersect
topo ← mesh
mesh, intersect ← algo
algo ← fit
algo, intersect, mesh ← bool
algo ← heal                  (unchanged; never bool)
fit, bool, mesh ← recon      (recon: core math geom topo mesh intersect algo fit bool)
fit, bool, heal ← fillet     (fillet gains fit)
fillet ← offset
algo, bool ← io, hlr         (io no longer needs mesh to solid)
everything ← ogeom           (facade gains fit and recon)
```

- There are no cycles. `ogeom-algo` and `ogeom-heal` still never depend on
  `ogeom-bool`.
- `ogeom-fit` does not depend on `ogeom-bool`.
- `ogeom-recon` depends on `ogeom-bool` only for `cut_by_mesh_region`. That is
  what lets `ogeom-bool`, and everything above it, stop rebuilding when mesh
  to solid changes.
- The only remaining cycle is the existing path-only dev-dependency
  `ogeom-mesh` → `ogeom-algo`. No new one is added.

**What each edit rebuilds:**

| Edited | Rebuilds |
|---|---|
| `ogeom-recon` | `ogeom-recon` and the facade |
| `ogeom-fit` | `ogeom-fit`, `ogeom-recon`, `ogeom-fillet`, `ogeom-offset` and the facade |
| `ogeom-algo` | as today, but `ogeom-algo` itself is 44% smaller, and `fit` and `recon` compile alongside `bool`, `heal`, `io` and `hlr` instead of before them |

### Deprecation path

- **Facade paths never change.** Every `ogeom::algo::*` and
  `ogeom::boolean::cut_by_mesh_region` path resolves before and after. The
  facade imports from the new crate, as `CONTRIBUTING.md` prescribes.
  `crates/ogeom/tests/downstream_paths.rs` enforces this: it lists
  `ogeom::algo::fit::*`, `ogeom::algo::register::*`, `solid_from_mesh`,
  `single_precision_quantum` and `ogeom::boolean::cut_by_mesh_region`.
- **Paths that move down keep their old crate path for one release** as a
  plain re-export with a doc line naming the new path and the release that
  removes it:
  - `ogeom_algo::{Spacing, interpolate, approximate, approximate_within,
    single_precision_quantum}`.
  - `rustc` rejects `#[deprecated]` on a `use` item (`useless_deprecated`), so
    the changelog fragment and the doc line carry the notice. This is how
    #1261 and #1266 kept theirs.
- **Paths that move up cannot be re-exported from below.** An old
  `ogeom_algo::solid_from_mesh` path would need `ogeom-algo` to depend on
  `ogeom-recon`, which depends on `ogeom-algo`: a cycle. The same holds for
  `ogeom_bool::cut_by_mesh_region`. These direct crate paths end in the release
  that moves them.
  - The facade's own documentation says the `ogeom-*` crate boundaries "will
    move".
  - The known downstream consumer uses only facade paths (#1272).
  - Each moved path is listed under "Downstream" in #1272 and in the
    changelog fragment, so a direct user can update in one pass.

### Migration for downstream

- **Facade users** (printCAD among them) change nothing. Every
  `ogeom::<module>::<item>` path resolves to the same item, with the same
  signature and behaviour.
- **Direct users of `ogeom-algo`** add `ogeom-fit` or `ogeom-recon` and rename
  `ogeom_algo::X` to `ogeom_fit::X` or `ogeom_recon::X` for the items listed
  above.
- **Direct users of `ogeom-bool`** rename `ogeom_bool::cut_by_mesh_region` to
  `ogeom_recon::cut_by_mesh_region`.
- **Users of the curve-fitting functions** may keep `ogeom_algo::interpolate`
  for one release, then use `ogeom_geom::fit::interpolate`. The same applies to
  `approximate`, `approximate_within`, `Spacing` and `single_precision_quantum`
  (→ `ogeom_mesh::single_precision_quantum`).
- **`outside/`** uses none of the moved items (checked). Its workspace still
  gets the version bump `CONTRIBUTING.md` asks for.

### Pull requests

Each PR moves code only and changes no behaviour. Each runs `./tools/check.sh`
and the full `cargo test -p ogeom -- --include-ignored`. Each touches
`changelog.d/1264-*.md` and the "Downstream" list in #1272.

1. **Prepare.**
   - `docs/SCOPE.md` records the two crates. Scope does not change: the
     "one addition" section names `ogeom-recon` as where it lives.
   - Widen the seven items above to `#[doc(hidden)] pub`, and add them to
     `INTERNAL`.
   - Move the `classify.rs` test that builds through `solid_from_mesh` to
     `crates/ogeom/tests/`. It tests classification of a converted solid, and
     the facade crate sees both.
   - Record the baseline build times (method under "Verification").
2. **Curve fitting down.**
   - Move the curve-fitting part of `fit.rs` to `ogeom_geom::fit`.
   - Add `ogeom_algo` root re-exports for one release.
   - Point the facade at `ogeom_geom`.
   - `ogeom-fillet` and the `ogeom-mesh` tests switch to `ogeom_geom::fit`
     paths.
3. **Quantum down.**
   - Move `single_precision_quantum` to `ogeom-mesh`.
   - `ogeom_algo` re-exports it for one release.
   - `ogeom-io` uses the new path.
4. **`ogeom-fit`.**
   - `git mv` the eleven modules. Add the crate manifest (metadata as the
     other crates) and the `[workspace.dependencies]` entry.
   - `ogeom-fillet` depends on it. The facade imports from it.
   - Move `crates/ogeom-algo/tests/{fit_grow, fit_segments_2d, orient,
     primitive_fits, register, rest, zone_fits}.rs`, or their modules in the
     single test program if #1267 has landed.
   - Update `docs/capabilities.toml` evidence paths and the ledger checksum.
   - Update the README crate table, `CONTRIBUTING.md` ("fifteen" becomes
     "seventeen" once both crates exist), and the `docs/DATA_MODEL.md` history
     list.
5. **`ogeom-recon`.**
   - `git mv` `mesh_solid` (with its directory), `recognize_patch` and
     `recognize_swept`, and `ogeom-bool/src/mesh_region.rs`.
   - Add a `test-hooks` feature. `crates/ogeom` enables it in its dev
     dependencies as it does for `ogeom-algo`.
   - Update `tools/converter-check.sh` (`-p ogeom-algo` becomes
     `-p ogeom-recon`), the capabilities paths, the README and
     `CONTRIBUTING.md`.
   - This PR meets the issue's acceptance: mesh to solid in its own crate, the
     facade path unchanged, build times before and after.

PRs 2 and 3 are independent of each other. PR 4 needs 1 and 2, and PR 5 needs
4. All five land in one release, so the old direct paths of PRs 4 and 5 end
together, and the re-exports of PRs 2 and 3 end one release later.

### Risks

- **Rebase conflicts with #1263 and #1270** (splitting large files). Those
  issues come first in #1272. Each move here is a `git mv` plus import edits,
  so it rebases on whatever split has landed. Splits inside `mesh_solid` or
  `fit` stay valid because module names are kept.
- **Lost `#[cfg(test)]` reach.** A test module inside a moved file may use
  core `pub(crate)` items beyond the seven above. The PR finds these at
  compile time. Each one is either widened as `#[doc(hidden)] pub` or the test
  moves to `crates/ogeom/tests`. No production `pub` is added for a test alone
  without saying so in the PR.
- **Rustdoc links.** Intra-doc links such as `crate::recognize_points` in
  core modules, and links into core from moved modules, break silently as
  warnings. `tools/check.sh` builds rustdoc with `-D warnings`.
- **Publishing order.** The workspace publish orders the crates from the
  graph. The dry run (`cargo publish --dry-run --workspace`) runs in PR 4 and
  PR 5.
- **Build time does not improve.** The measurement is part of the
  acceptance. If PR 5 does not cut the rebuild after a one-line edit in
  `mesh_solid` by at least half, the PR says so and #1272 records it. The code
  move still stands, because of the cohesion and the smaller `ogeom-algo`.
- **Hidden behaviour change.** None is expected from moves. The real-world
  STEP import check and `tools/converter-check.sh` run before PR 4 and PR 5
  merge, and must match the baseline face for face.

### Verification

Each PR is verified as follows:

- `./tools/check.sh`
- `python3 tools/facade.py check`
- `cargo test -p ogeom -- --include-ignored`
- `python3 tools/capabilities.py check` (PR 4 and PR 5; `./tools/check.sh`
  also runs it, after rustdoc with `-D warnings`)
- `cargo publish --dry-run --workspace` (PR 4 and PR 5)

Also:

- **Tests are moved, never edited.** Each moved test file's diff is a rename
  plus import lines. The count of tests run per crate before and after is
  listed in the PR, and the sum must match.
- **Build time** is measured with the same command and machine before and
  after:
  - a one-line comment edit in `mesh_solid/regions.rs` (PR 5) and in
    `fit/grow.rs` (PR 4), then `cargo build -p ogeom --tests`;
  - the wall time and the list of crates recompiled (`cargo build -v`) are
    recorded;
  - the same edit in `build.rs` shows the core's cost is unchanged.

## Questions decided

- **Where `recognize`, `orient` and `rest` go: `ogeom-fit`.** The fits depend
  on the first two, and `rest` needs `fit::segments` and `orient`. Keeping
  them in the core would mean moving those helpers down, which is a redesign,
  not a move.
- **`pcurve_fit`, `length`, `bins` and `traced` stay**: construction uses
  them.
- **`cut_by_mesh_region` goes up into `ogeom-recon`**, instead of keeping
  `ogeom-recon` below `ogeom-bool`. Keeping the reverse dependency would make
  every mesh-to-solid edit rebuild booleans, fillets, offsets, I/O and HLR,
  which is most of the gain. Its facade path stays.
- **Curve interpolation goes down to `ogeom-geom`**, not into `ogeom-fit`. It
  is curve construction with no dependency above geometry. It avoids a new
  dev-dependency cycle, and its old `ogeom_algo` path can be kept one release.
- **Direct crate paths of upward moves end in the release that moves them.**
  A re-export would be a dependency cycle. Facade paths, which are the API,
  never move.
- **The medial and distance group and the mass and check group stay** in
  `ogeom-algo`. Construction and query use them (`check` uses `proximity`,
  `register` uses `distance`), and nothing above the core gains from moving
  them.

Decided: implementation may start.
