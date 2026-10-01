# Kernel test methodology

The right tool per layer. Each contract has one owner at the strongest
boundary that is cheap and deterministic enough. There is no external
oracle: comparing against another kernel is not possible, because it would
mean vendoring one. Every oracle below can be computed here.

| Layer | Owner | Oracle and tool |
|---|---|---|
| Laws: transforms, predicates, knot vectors, basis functions, tolerance arithmetic | Property test over random inputs | `proptest`; the law stated in the test name (round trip, composition, antisymmetry, containment); see `crates/ogeom-core/tests/properties.rs` |
| Analytic geometry: primitives, extrusions, revolutions, offsets of analytic surfaces | Closed form | volume, area and centroid from the formula, at a stated tolerance; the same result built two ways |
| Curves and surfaces: evaluation, derivatives, conversion to splines | Agreement to a stated error | finite differences against derivatives; the converted spline sampled against the original; a jet against the separate accessors to rounding |
| Intersections and sections | Both surfaces | every point of the section on both surfaces within the stated tolerance; the loop closed; the count of sections |
| Booleans, blends, shells | Invariants plus a closed form where one exists | `check` passes; volume and area against the formula or against inclusion and exclusion (fuse plus common equals the sum); the history names every input |
| Refusals | The text | the error says which thing and why, asserted on the text, not on `is_err()` |
| Meshing | Agreement with the surface | every vertex within the deflection; shared edges agree; the volume of the mesh against the solid's |
| Exchange formats | Round trip and outside-authored input | write then read gives the same measurement; `tests/corpus/` files (every one with a manifest entry) read to the expected counts and volumes |
| Drawings | The closed form of the outline | an edge drawn once; the silhouette of an analytic solid against its formula |
| Robustness at the placements nobody wrote down | The stress harness | `tools/ogeom-stress`, seeded, judged per case by an oracle, gated against `baseline.json` nightly |
| Performance | Watched, not gated | `tools/ogeom-bench` ratios against `baseline.json` |
| Parity | The ledger's evidence | `docs/parity/parity.toml` cites the symbol and the test; `tools/parity.py check` holds it to them |
| The guide | The book's own code | `crates/ogeom/tests/book.rs` includes the examples by anchor; a book that builds is a book whose examples passed |

## Rules

- A test states its tolerance and why. A tolerance tightened to pass is a
  bug report; a tolerance loosened to pass needs its reason in the commit.
- A corpus file is copied in with a manifest entry, never referenced from
  outside the tree, and never taken from another kernel's test suite.
- Quality signal is whether a deliberate wrong answer goes red, not line
  coverage. When a test's value is in doubt, make one mutation of the owner
  (flip a sign, drop a face, skip the location) and confirm the test fails,
  then restore the source byte for byte.
- Property tests run twice in `tools/check.sh`. A property that fails on
  some seeds is a real case; add it by name.
- No file the user tested with is named in a test or a comment. The test
  names the geometric condition.
