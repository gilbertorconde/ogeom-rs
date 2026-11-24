# Design note: one blend request instead of option-named entry points (#1265)

Tracked in #1272. Depends on #1260 (the explicit facade), which is merged.
ADR 0001 governs the blend model. This note decides the public shape of a
blend request. An ADR 0001 amendment that points here lands with the first
code PR.

## Verdict: valid

Re-counted on `main` at `09979ec1`, from the `pub use` block in
`crates/ogeom-fillet/src/lib.rs` and the bodies behind it.

**Edge-blend entry points whose names spell options: 43 functions.** Each
new option still multiplies the set.

| Family | Functions |
|---|---|
| Uniform rounds (11) | `fillet_edge`, `fillet_edges`, `fillet_edges_with`, `fillet_edges_overflow`, `fillet_edges_corner`, `fillet_edges_corner_overflow`, `fillet_edges_cornered`, `fillet_edges_checked`, `fillet_edges_corner_checked`, `fillet_edges_corner_overflow_checked`, `fillet_edges_partial` |
| Laws (3) | `fillet_edge_variable`, `fillet_edges_variable`, `fillet_edges_variable_overflow` |
| Distances (4) | `fillet_edge_distances`, `fillet_edge_distances_overflow`, `fillet_edge_distances_variable`, `fillet_edges_distances_variable` |
| Chamfers (10) | `chamfer_edge`, `chamfer_edge_angle`, `chamfer_edge_distances`, `chamfer_edges`, `chamfer_edges_with`, `chamfer_edges_cornered`, `chamfer_edges_overflow`, `chamfer_edges_checked`, `chamfer_edges_overflow_checked`, `chamfer_edges_partial` |
| Mixed (2) | `blend_edges`, `blend_edges_overflow` |
| Largest-size readings (13) | `largest_fillet_radius{,_all,_overflow}`, `largest_fillet_distances`, `largest_chamfer_distance{,_all,_overflow,_all_overflow}`, `largest_law_radius`, `largest_law_scale`, `largest_law_distances`, and `blend_limit` / `blend_seat` over `BlendSizing` |

What changed since the issue was written:

- **`BlendSizing` repeats the problem for readings.** It is a public enum
  whose eleven variants are named by options: `Rounds`, `RoundsOverflow`,
  `RoundsCornerOverflow`, `RoundDistancesOverflow`, `RoundLawOverflow` and so
  on.
- **The checked calls read refusals two ways.**
  - `fillet_edges_checked`, `fillet_edges_corner_checked` and
    `chamfer_edges_checked` name the first edge whose prefix of the request
    fails, then the earlier edges it meets.
  - The `*_overflow_checked` calls name every refused edge, read in order.
    ADR 0001's amendment "checked blends under overflow" (#1240) argues this
    is what applications need.
- **Five one-edge functions are their own constructions, not wrappers.**
  `fillet_edge`, `fillet_edge_variable`, `chamfer_edge`, `chamfer_edge_angle`
  and `chamfer_edge_distances` have narrower domains (a straight edge between
  planes) and their own refusal texts. `fillet_edge` documents one real
  difference: it "stops flush" against a neighbouring blend instead of
  running on through it.

### Which functions are wrappers today

Read from the bodies.

**Pure delegation**, with defaults filled in:

- `fillet_edges` → `fillet_edges_overflow`
- `fillet_edges_corner` → `fillet_edges_corner_overflow`, which with `Setback`
  → `fillet_edges_overflow`
- `fillet_edges_variable` → `_variable_overflow`
- `fillet_edge_distances` → `_overflow`
- `fillet_edges_distances_variable` on one edge → `fillet_edge_distances_variable`
- `fillet_edges_corner_checked` with `Setback` → `fillet_edges_checked`
- `chamfer_edges` → `chamfer_edges_with`
- `chamfer_edges_cornered` → `chamfer_edges_overflow`
- `blend_edges` → `blend_edges_overflow`, which with bevels alone →
  `chamfer_edges_overflow`, and with rounds alone → `fillet_edges_with`
- every `largest_*` → `blend_limit` over one `BlendSizing` variant

**Distinct bodies: 23.**

- Rounds: `fillet_edges_overflow`, `fillet_edges_corner_overflow` (rolling),
  `fillet_edges_with`, `fillet_edges_cornered`,
  `fillet_edges_variable_overflow`.
- Distances: `fillet_edge_distances_overflow`, `fillet_edge_distances_variable`.
- Chamfers: `chamfer_edges_with` (through `chamfer_edges_blamed`),
  `chamfer_edges_overflow`.
- Mixed: `blend_edges_overflow`.
- Readings: `blend_limit` and `blend_seat`.
- Checked: the two prefix readings and the two in-order readings.
- Partial: the two partial readings.
- One edge: the five one-edge constructions.

So 43 public names sit over 23 bodies. The engine already reads one internal
request (`engine/request.rs`: blends of `(edge, Section)`, `mitre`, `corners`,
`rolling`), which is ADR 0001 §3's "set of `(edge, section)` pairs".

### What the ADR and the M6 plan settle

- **ADR 0001 §3 fixes the request's meaning**: a set of `(edge, section)`
  pairs, with sections circular, conic or distances, chamfer, or any of these
  under a law. The result is a function of the set, not of its order. Its
  amendments fix the options:
  - "overflow handling" defines `Overflow`, tried in order, every behaviour
    on by default;
  - "rolling corners" defines `FilletCorner`;
  - §13 question 1 makes the facet the default chamfer corner;
  - §13 question 2 adds setbacks only where asked.
- **§11 says every migration step keeps the public API.** Nothing in the ADR
  fixes function names.
- **`docs/m6-removal-plan.md`** keeps "API contracts" and "the wrappers" while
  old construction paths retire. It retired old constructions, never public
  names.

Neither document settles the API's shape, so this note does, inside ADR
0001's semantics.

### Downstream use

The issue names these, and `crates/ogeom/tests/downstream_paths.rs` lists
them: `fillet_edges`, `fillet_edges_variable`, `fillet_edges_corner`,
`fillet_edges_corner_checked`, `fillet_edge_distances`, `blend_edges`,
`chamfer_edges_checked`, `chamfer_edges_cornered`, `chamfer_edges_with`,
`largest_chamfer_distance`, `largest_fillet_distances`,
`largest_fillet_radius_all`, and the types `Blend`, `BlendRefusal`, `Chamfer`,
`FilletCorner`, `ChamferCorner`, `RadiusLaw` and `RoundProfile`.

`fillet_faces`, `fillet_sheet_edges`, `blend_faces`, `analyse_blend` and
`make_blend_surface` are other operations, and stay as they are.

## Decision

### The request type

A builder over a value. The fields are private, and every option is a method,
so a new option is one more method and breaks no caller. It lives in
`crates/ogeom-fillet/src/request.rs`, beside the engine's internal `Request`.

```rust
/// Every edge blend of one operation, read as one set (ADR 0001, section 3).
/// The order edges are added in carries no meaning for the shape built; it
/// only numbers them for `BlendRefusal::edges`.
#[derive(Debug, Clone, Default)]
#[must_use]
pub struct BlendRequest { /* private */ }

/// What one blend sweeps along its edge or chain.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum BlendSection {
    /// A circular round: the rolling ball of this radius.
    Round(f64),
    /// A round whose radius follows the law along the chain.
    RoundLaw(RadiusLaw),
    /// A round meeting `face` `on_face` from the edge and the other face
    /// `on_other` from it, its section `profile` between the contacts.
    RoundDistances { face: Shape, on_face: f64, on_other: f64, profile: RoundProfile },
    /// As `RoundDistances`, each distance following its law along the chain.
    DistanceLaws { face: Shape, on_face: RadiusLaw, on_other: RadiusLaw, profile: RoundProfile },
    /// A flat bevel: symmetric, two distances, or a distance and an angle.
    Chamfer(Chamfer),
}

impl From<Chamfer> for BlendSection { /* BlendSection::Chamfer */ }
impl From<Blend> for BlendSection { /* Round, Chamfer, RoundDistances map one to one */ }

/// Which size a largest-size reading varies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Vary {
    /// Every length of the request times one factor; angles, profiles and
    /// the shape of each law held. The reading is the factor.
    #[default]
    Scale,
    /// One station of the request's round law, the rest held.
    Law(LawStation),
    /// One station of the `on_face` law of a `DistanceLaws` section.
    OnFace(LawStation),
    /// One station of its `on_other` law.
    OnOther(LawStation),
}

impl BlendRequest {
    pub fn new() -> Self;

    // What is blended.
    /// Each edge rounded at `radius`.
    pub fn round(self, edges: impl IntoIterator<Item = Shape>, radius: f64) -> Self;
    /// Each edge bevelled `distance` back along both faces.
    pub fn chamfer(self, edges: impl IntoIterator<Item = Shape>, distance: f64) -> Self;
    /// One edge with any section.
    pub fn edge(self, edge: Shape, section: impl Into<BlendSection>) -> Self;
    /// A chain under one section. A law runs along the whole chain, as
    /// `fillet_edges_variable` takes it. A constant section on a chain is
    /// that section on each edge.
    pub fn chain(self, chain: impl IntoIterator<Item = Shape>, section: impl Into<BlendSection>) -> Self;

    // How.
    /// How a vertex where rounds meet closes. Default `Setback`.
    pub fn round_corner(self, corner: FilletCorner) -> Self;
    /// How a vertex where three or more bevels meet closes. Default `Facet`.
    pub fn chamfer_corner(self, corner: ChamferCorner) -> Self;
    /// A vertex closed as a setback corner with these distances (ADR 0001,
    /// section 13, question 2).
    pub fn setback(self, vertex: Shape, setbacks: Setbacks) -> Self;
    /// What a blend reaching other geometry does. Default: every behaviour on.
    pub fn overflow(self, overflow: Overflow) -> Self;

    // The operation and its readings.
    pub fn build(&self, model: &mut Model, solid: &Shape, tol: Tolerances)
        -> OgeomResult<Built>;
    pub fn build_checked(&self, model: &mut Model, solid: &Shape, tol: Tolerances)
        -> Result<Built, BlendRefusal>;
    pub fn build_partial(&self, model: &mut Model, solid: &Shape, tol: Tolerances)
        -> OgeomResult<(Built, Vec<BlendRefusal>)>;
    pub fn limit(&self, model: &Model, solid: &Shape, vary: Vary, up_to: f64, tol: Tolerances)
        -> OgeomResult<BlendLimit>;
    /// `limit(..)?.size`.
    pub fn largest(&self, model: &Model, solid: &Shape, vary: Vary, up_to: f64, tol: Tolerances)
        -> OgeomResult<Option<f64>>;
    /// The seat reading `limit` starts from, without building.
    pub fn seat(&self, model: &Model, solid: &Shape, vary: Vary, up_to: f64, tol: Tolerances)
        -> OgeomResult<Option<f64>>;
}
```

Also:

- **`ChamferCorner` gains `impl Default` (`Facet`)**, the ADR's default.
  Adding the impl breaks nothing.
- **These stay as they are:** `Overflow`, `FilletCorner`, `Chamfer`,
  `RadiusLaw`, `RoundProfile`, `Setbacks`, `BlendRefusal`, `BlendLimit`,
  `LawStation` and the `LIMIT_*` constants. `BlendRefusal::edges` indexes the
  request's edges in the order they were added, with chains flattened.

### How the old readings map onto the request

| Reading today | On the request |
|---|---|
| `_overflow` (`Overflow`) | `.overflow(o)`, for rounds and bevels alike, since both calls take the same struct |
| `_corner` (`FilletCorner`) | `.round_corner(c)` |
| `_cornered` (`ChamferCorner`) | `.chamfer_corner(c)` |
| `fillet_edges_cornered` (vertices with `Setbacks`) | `.setback(v, s)` per vertex |
| `_variable` (laws) | `BlendSection::RoundLaw` on `.chain(..)`. `fillet_edge_variable(e, a, b)` is `RoundLaw(RadiusLaw::Linear { start: a, end: b })` on that edge |
| `_distances` | `BlendSection::RoundDistances`, `DistanceLaws`, or `Chamfer::Distances` / `Chamfer::Angle` |
| `_with` (per-edge sizes) | one `.edge(e, section)` per edge |
| `_checked` | `.build_checked(..)`: every refused edge, read in order (#1240) |
| `_partial` | `.build_partial(..)` |
| `largest_*`, `blend_limit`, `blend_seat` | `.largest(..)`, `.limit(..)`, `.seat(..)` with a `Vary` |

**The readings are a scale.** `Vary::Scale` reads the largest factor `k` up to
`up_to` at which the request, every length times `k`, builds. This one rule
covers each old reading:

- the radius of a uniform round, asked at radius `1.0`;
- a symmetric or angle chamfer's distance, asked at distance `1.0`;
- the common scale of two distances, which `largest_fillet_distances` and
  `Chamfer::Distances` already read;
- `LawSize::Scale`.

The old functions pass unit sizes, so `k · 1.0` is `k` bit for bit. They probe
exactly the sizes they probe today.

### How a request reaches today's code

`build`, `limit` and the rest **choose an existing body and call it with the
arguments it takes today.** The dispatcher validates nothing itself. Every
refusal text stays the body's own: an empty request, a radius that cannot be
used, an edge asked twice. An empty `round(..)` or `chamfer(..)` call routes by
the kind its builder method recorded. Where the options are at their
defaults, the default call's body is used.

**`build`:**

| Request | Body |
|---|---|
| rounds, one radius, no setbacks, `Setback` | `fillet_edges_overflow(edges, r, overflow)` |
| rounds, one radius, `Rolling` | `fillet_edges_corner_overflow(.., Rolling, overflow)` |
| rounds, one radius, setbacks, `Setback`, default overflow | `fillet_edges_cornered` |
| rounds, several radii, `Setback`, default overflow | `fillet_edges_with` |
| rounds, several radii, `Setback`, other overflow | `blend_edges_overflow` |
| chamfers alone | `chamfer_edges_with` when `Facet` with default overflow (subject to audit (d) below), otherwise `chamfer_edges_overflow(specs, corner, overflow)` |
| chamfers with rounds or `RoundDistances`, `Setback` | `blend_edges_overflow(pairs, chamfer_corner, overflow)` |
| one `RoundDistances` edge | `fillet_edge_distances_overflow` |
| one chain of `RoundLaw` | `fillet_edges_variable_overflow` |
| one chain of `DistanceLaws`, default overflow | `fillet_edges_distances_variable`, which hands a chain of one to `fillet_edge_distances_variable` |

**`build_checked`:**

- uniform rounds, either corner, any overflow →
  `fillet_edges_corner_overflow_checked`;
- chamfers, either corner, any overflow → `chamfer_edges_overflow_checked`.

**`build_partial`:**

- uniform rounds with defaults → `fillet_edges_partial`;
- chamfers with defaults → `chamfer_edges_partial`.

**`limit` and `seat`** lower to today's `BlendSizing` variant, which becomes
crate-internal:

- uniform rounds → `Rounds`, `RoundsOverflow` or `RoundsCornerOverflow`;
- one `RoundDistances` edge → `RoundDistances` or `RoundDistancesOverflow`;
- chamfers → `Chamfer` for one edge with defaults, otherwise `Chamfers`;
- a `RoundLaw` chain → `RoundLaw` or `RoundLawOverflow`, with
  `Vary::Law(s)` → `LawSize::At(s)` and `Scale` → `LawSize::Scale`;
- one `DistanceLaws` edge → `DistanceLaws`, with `OnFace` / `OnOther` /
  `Scale` → `DistanceLawSize`.

**Any other combination is refused by name** with
`OgeomError::Construction`. Examples: a law under a rolling corner, setbacks
with overflow options, a checked mixed request, a law reading of a request
with no law. The text names what was asked and says no construction builds it
("a round law with rolling corners is not built"). One test per refused row
pins it, following "Refusals are by name, and pinned". Lifting a refusal later
is a PR with its own tests and no API change. That is the point of the
request.

### Tuples two old functions share

Several old functions lower to the same request with the same options.
Where their bodies are not one call chain, PR 1 audits them:

- **(a)** `fillet_edge(e, r)` vs `fillet_edges(&[e], r)`
- **(b)** `chamfer_edge`, `chamfer_edge_angle` and `chamfer_edge_distances` vs
  `chamfer_edges_with(&[(e, spec)])`
- **(c)** `fillet_edge_variable(e, a, b)` vs `fillet_edges_variable(&[e], Linear)`
- **(d)** `chamfer_edges_cornered(Facet)` vs `chamfer_edges_with`
- **(e)** the prefix-reading checked calls vs the in-order ones

A request always lowers to the **general** call's body: the multi-edge,
default-options one. For each pair, the audit runs the other function's own
test files through that body. If they pass, the function becomes a wrapper.
If not, it stays a verbatim deprecated function over its own body for one
release, and the difference is recorded in the ADR amendment and the
changelog.

Expected outcome, decided now:

- **The five one-edge constructions (a, b, c) stay verbatim.** They are the
  older straight-edge family that ADR 0001 §10 retires. Making them wrappers
  would change what they build or refuse.
- **The three prefix-reading checked calls (e) stay verbatim.**
  `build_checked` uses the #1240 reading.
- **(d) follows the audit.**

Every other function becomes a wrapper.

### Five entry points on the request

```rust
// fillet_edges(model, &solid, &edges, 2.0, tol)
BlendRequest::new()
    .round(edges.iter().cloned(), 2.0)
    .build(model, &solid, tol)

// fillet_edges_corner_overflow_checked(model, &solid, &edges, 2.0, FilletCorner::Rolling, overflow, tol)
BlendRequest::new()
    .round(edges.iter().cloned(), 2.0)
    .round_corner(FilletCorner::Rolling)
    .overflow(overflow)
    .build_checked(model, &solid, tol)

// fillet_edges_variable_overflow(model, &solid, &chain, &law, overflow, tol)
BlendRequest::new()
    .chain(chain.iter().cloned(), BlendSection::RoundLaw(law.clone()))
    .overflow(overflow)
    .build(model, &solid, tol)

// chamfer_edges_overflow(model, &solid, &specs, ChamferCorner::Point, overflow, tol)
specs
    .iter()
    .fold(BlendRequest::new(), |r, (e, c)| r.edge(e.clone(), c.clone()))
    .chamfer_corner(ChamferCorner::Point)
    .overflow(overflow)
    .build(model, &solid, tol)

// largest_law_distances(model, &solid, &edge, &face, &d1, &d2, profile,
//                       DistanceLawSize::OnFace(LawStation::Start), 8.0, tol)
BlendRequest::new()
    .edge(edge.clone(), BlendSection::DistanceLaws {
        face: face.clone(), on_face: d1.clone(), on_other: d2.clone(), profile,
    })
    .largest(model, &solid, Vary::OnFace(LawStation::Start), 8.0, tol)
```

And one each for setbacks and the radius reading:

```rust
// fillet_edges_cornered(model, &solid, &edges, 2.0, &[(v, setbacks)], tol)
BlendRequest::new().round(edges.iter().cloned(), 2.0).setback(v, setbacks).build(model, &solid, tol)

// largest_fillet_radius_all(model, &solid, &edges, 10.0, tol)
BlendRequest::new().round(edges.iter().cloned(), 1.0).largest(model, &solid, Vary::Scale, 10.0, tol)
```

### What stays out of the request

These are other operations, not option variants of an edge blend. They keep
their functions:

- vertex rounds: `round_vertex`, `round_vertex_setback`,
  `round_vertex_setback_sharp`, the local corner feature of the amendment
  "vertex-only rounding";
- face-pair blends: `blend_faces`, `fillet_faces`;
- sheet edges: `fillet_sheet_edges`, a shell rather than a solid;
- 2D corners: `fillet_corner_2d`, `chamfer_corner_2d`;
- bridges: `make_blend_curve`, `make_blend_surface`;
- the marching primitives: `march_blend*`;
- the analyses: `analyse_blend`, `face_curvature_*`;
- chain selection: `tangent_chain*`.

### Deprecation path

- **Release N** is the release the wrappers land in, expected 0.22.0.
  - The 43 functions and the types `Blend`, `BlendSizing`, `LawSize` and
    `DistanceLawSize` get
    `#[deprecated(since = "0.N.0", note = "use BlendRequest::...")]`, the
    note naming the exact call.
  - They keep compiling and behaving as before.
  - The facade keeps re-exporting them next to `BlendRequest`, `BlendSection`
    and `Vary`.
- **Release N + 1** removes them, with the verbatim one-edge and
  prefix-reading bodies.
  - Their tests are ported to the request, never deleted. That covers the
    regression tests of #969 to #976 named in `docs/m6-removal-plan.md`.
  - #1272's "Downstream" list records each removed path.

### Migration for downstream

The paths printCAD uses, and what each becomes:

| Today | On the request |
|---|---|
| `fillet_edges(m, s, &e, r, t)` | `BlendRequest::new().round(e, r).build(m, s, t)` |
| `fillet_edges_variable(m, s, &c, &law, t)` | `.chain(c, BlendSection::RoundLaw(law)).build(..)` |
| `fillet_edges_corner(.., corner, ..)` | `.round(e, r).round_corner(corner).build(..)` |
| `fillet_edges_corner_checked(..)` | `.round(e, r).round_corner(corner).build_checked(..)` (see note 1) |
| `fillet_edge_distances(m, s, e, f, d1, d2, p, t)` | `.edge(e, BlendSection::RoundDistances { face: f, on_face: d1, on_other: d2, profile: p }).build(..)` |
| `blend_edges(m, s, &pairs, corner, t)` | one `.edge(e, blend)` per pair (`Blend` converts), then `.chamfer_corner(corner).build(..)` |
| `chamfer_edges_with(m, s, &specs, t)` | one `.edge(e, spec)` per spec, then `.build(..)` |
| `chamfer_edges_cornered(.., corner, ..)` | the same, then `.chamfer_corner(corner).build(..)` |
| `chamfer_edges_checked(..)` | the same, then `.build_checked(..)` (see note 1) |
| `largest_fillet_radius_all(m, s, &e, up, t)` | `.round(e, 1.0).largest(m, s, Vary::Scale, up, t)` |
| `largest_fillet_distances(m, s, e, f, d1, d2, up, t)` | `.edge(e, RoundDistances { .. }).largest(.., Vary::Scale, ..)` |
| `largest_chamfer_distance(m, s, e, spec, up, t)` | `.edge(e, spec at measure 1.0).largest(.., Vary::Scale, ..)` |
| `Blend` | `BlendSection` (`From<Blend>`) |

Notes:

1. **The checked calls name edges differently.** The request names every
   edge it refuses, in request order, so that the request without them
   builds. The prefix-reading calls named the first culprit and its
   neighbours. An application that offers to drop the named edges gets one
   rebuild that succeeds, which is what #1240 argued for.
2. **These paths are unchanged:** `FilletCorner`, `ChamferCorner`,
   `RadiusLaw`, `RoundProfile`, `Chamfer`, `BlendRefusal`, `fillet_faces`,
   `fillet_sheet_edges`, `blend_faces`, `analyse_blend` and
   `make_blend_surface`.
3. **Inside this repository**, the users to migrate are
   `crates/ogeom-offset/src/shape.rs`, `tools/ogeom-stress`, `ogeom-bench`,
   `tools/ogeom-wasm/tests` and `outside/` (`fillet_edge`, `chamfer_edge`).

### Pull requests

1. **ADR 0001 amendment "the public request"** (docs only).
   - Points to this note.
   - Records the audit of tuples (a) to (e): which functions
     become wrappers and which stay verbatim, with the test files that
     decided it.
2. **`BlendRequest`, `BlendSection`, `Vary`, and `Default` for
   `ChamferCorner`.**
   - Add the dispatcher and the refusals by name. Old functions untouched.
   - Add to the facade and `tools/facade.py`, plus a changelog fragment.
   - New tests:
     - one per dispatch row: the request and the old call on the same
       corpus part give the same faces and bit-identical volume;
     - one per refused combination, pinning its text.
3. **Old functions become wrappers** (those the audit admits).
   - `BlendSizing` lowering moves behind the request.
   - `blend_limit` and `blend_seat` become wrappers.
   - No test file changes: the existing 75 fillet test files and the facade
     tests that call old names pass unchanged.
4. **Migrate in-repository callers** to the request (the list under
   "Migration", note 3). `outside/` takes the version bump.
5. **Deprecate.**
   - `#[deprecated]` on the 43 functions and the four types.
   - `#![allow(deprecated)]` at the top of each test file that still
     exercises them, including `downstream_paths.rs`. No test body changes.
   - Add the release N entry under #1272 "Downstream".
6. **Release N + 1: remove**, porting each removed function's tests to the
   request.

Later PRs may lift refusals: a checked or partial reading of mixed requests,
laws under rolling corners. Each comes with tests and needs no API change.

### Risks

- **The dispatcher picks a body the old call did not use.** The audit in PR 1
  and the per-row equivalence tests in PR 2 guard this, and PR 3 rests on
  unchanged tests passing. The real-world import check and the fillet stress
  scenario run before PR 3 merges and must match the baseline.
- **Refusal texts change.** Validation stays in the bodies, so no text
  changes. The equivalence tests compare errors as well as results.
- **`BlendRefusal::edges` indexes chains.** A chain is flattened into its
  edges. The checked readings that exist index per edge and do not take
  laws, so no existing index changes meaning.
- **Name clashes.** `fillet` and `chamfer` are public modules, so the request
  does not add free functions of those names. Its operations are methods.
- **Deprecation noise under `-D warnings`.** Internal callers stop using
  deprecated items in PR 4, before PR 5 adds the attribute. Test files get a
  crate-level allow, and wrappers that call each other carry
  `#[allow(deprecated)]`.
- **Performance.** The request owns its shapes, so building one clones the
  edge handles: cheap reference handles, negligible beside any blend.

### Verification

Every PR runs:

- `./tools/check.sh`
- `cargo test -p ogeom-fillet -- --include-ignored`
- `cargo test -p ogeom -- --include-ignored`
- `python3 tools/facade.py check`
- the MSRV and CI-stable clippy that `tools/check.sh` covers

Also:

- **PR 2:** an equivalence test per dispatch row, with result faces, volume
  to the bit and error text compared, and a pinned test per refusal.
- **PR 3:** the test tree's diff is empty. The count of tests run per crate
  before and after is equal.
- **PR 3 and PR 6:** the real-world import check runs.
- **PR 5:** the only test-tree diff is the added allow lines.

## Rationale

- **One request for fillets and chamfers**, not separate `fillet` and
  `chamfer` entry points. ADR 0001 is one blend model for both. A request
  mixing the two is already public (`blend_edges`), and a split by kind would
  rebuild that mixed entry point beside them. `round(..)` and `chamfer(..)`
  keep the common cases one call long.
- **A builder with private fields**, not a struct with public fields. This is
  the mainstream kernel API practice for blends: one blend operation object
  that takes edges or chains, each with its section (constant, two distances,
  distance and angle, or a law), plus global options (corner type, overflow
  behaviours, setbacks). It is then built, and asked which edges failed and
  how large a blend fits. In Rust, the guidelines' answer to many optional
  parameters is a builder, with `#[non_exhaustive]` so options and section
  kinds can be added without breaking callers. `Overflow` is already
  `#[non_exhaustive]`.
- **Checked and partial are methods, not a mode field.** Their return types
  and costs differ: a checked refusal builds once per edge. A mode field would
  make every caller unpack a sum type. Separate methods also keep the cost
  visible where the call is made.
- **Readings are a scale.** A scale is unit-consistent for every section
  kind, and the old readings are its unit-size cases, bit for bit.
- **Unsupported combinations are refused by name, never approximated.** They
  are not mapped to the nearest supported one. That is the project's rule
  ("No stubs"). Each such combination is now one refusal row instead of a
  missing function name.
- **The one-edge constructions and the prefix reading are not kept in the
  new API.** ADR 0001 retires the first, and its amendment #1240 replaced the
  second. Keeping them only as deprecated functions for one release lets
  callers move without a silent change.

Decided: implementation may start.
