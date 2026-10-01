# Plan

What is left to build, and the rules the work follows.

## Where the remaining work is tracked

- **GitHub issues** hold bugs and requests, one per issue.
- **The parity ledger** (`docs/PARITY.md`, generated from
  `docs/parity/parity.toml`) audits the reference kernel's modelling modules
  against named capabilities. The restriction on each `partial` row is the live
  list of what that capability does not do yet. `tools/check.sh` holds the
  ledger to its evidence.
- **This file** holds the open items that are neither: known limits the code
  refuses by name and points here for.

Finished work is recorded in `CHANGELOG.md` and in the tests that pin it, not
here.

## Open

**Edge and vertex contact in the boolean.** A piece lying on the other solid's
boundary with no coincident partner face is refused. What reaches that refusal
is genuine contact confined to a line or a point: there is no shared region to
find a partner in, and no pair of normals to compare. Resolving it needs the
contact classified from the neighbouring pieces instead of the piece itself.

**Coincident spline patches.** Two B-spline surfaces lying on one another (a
sheared copy sharing a plane with its original) are detected by sampling, and
each edge is carried into the other's chart over its stretch on the other's
window. An edge that leaves that window and comes back onto it is still
refused.

**Tangential contact as its own curve.** Two surfaces touching along a curve
without crossing it, where no closed form names the curve (a torus resting in
a free-form cradle), march into fragments on both surfaces that describe
nothing. `crates/ogeom-intersect/tests/hard.rs` pins the limit. The answer is a
tangential trace that follows the contact as a curve of its own kind.

**Blends between curved faces.** The face-to-face blend takes planar supports.
Curved ones need the marching seat: the spine where the two offset surfaces
meet, walked like the edge blend's.

**Three drills in `nist_ftc_06` are refused.** Replay with
`ogeom-stress --case drill/nist_ftc_06_asme1_rd/<n>`; a refusal prints the
drill's placement.

- Drill 4 ("the kept pieces did not close into a shell"): face 89 is a
  hemisphere bounded by one great circle through both poles. Its stored pcurve
  is a fit that strays far from the meridian (u = -0.30 where the circle is at
  u = 0, v = -0.65), so the drill's section ends 0.3 from where the boundary is
  paved and the strands do not join. Wanted: such a circle split at the poles
  into two exact half meridians, the neighbour rebuilt on the halves.
- Drills 0 and 3 ("arrangement left no piece of the face"): faces 4 and 86 are
  part cylinders whose chart window lies at negative `u` (-3.016 to -1.571).
  Their side edges are paved where the drill's section crosses them, but no
  section strand reaches the face, and the arrangement finds no cycle. Wanted:
  the section's pieces folded into the face's own window before they are kept
  or dropped.

**Facet walls drafted less than their chords sag.** A pad pushed down from
the top of a slab converted face for facet, whose rounded corner's rows hold
different numbers of facets and lean in by a few thousandths of a
millimetre (less than the corner's chords sag, so the rows' edges cross the
pad's walls in plan), refuses at some row counts: the kept pieces do not
close. The rows' planes and the pad's walls meet at small angles all round
the corner, the near-coplanar sliver band again.

**Speed, not correctness.**

- A per-face state cache would let the boolean's build phase read a piece's
  classification instead of probing for it. No case needs it for correctness.
- Coincidence of two patches is measured in the boolean because the projection
  it needs lives in `ogeom-algo`, above the intersector. Moving
  `project_on_surface` down into `ogeom-geom` would let `intersect_surfaces`
  answer `Same` itself, for every caller.
- Local operations, so an op costs what an edit touches, not what the solid
  holds:
  - Local boolean, in part. A face the tool leaves alone is neither split nor
    arranged, the other solid's trims and outlines are drawn only where
    asked, and an edge rebuilt whole is built once for its two faces. Still
    whole-solid: the pair filter (all pairs, by box), the rebuild of every
    kept face as new nodes, and the sew over all of them; a face the tool
    crosses is arranged whole, its untouched holes with it.
  - Unchanged faces in the history: done as `History::copy`, an exact copy
    on new nodes. Sharing the nodes themselves would need every operation
    that edits its result in place to copy on write first.
  - Local refine: done, `unify_same_domain_around`.
  - Fillets in one pass. Each blend piece and each corner is one whole-solid
    boolean today; wanted, all of a fillet's blends in one boolean, or
    replaced faces locally.
  - Face bounds in the model: not kept; `shape_bounds` and `tight_bounds`
    take about a millisecond on a four-hundred-face part, which a cache
    and its invalidation would not improve on.

**SAT, X\_T and JT are refused.** Their specifications are unpublished.
Implementing them would mean reverse engineering files instead of reading a
standard. If a specification is published, the refusal lifts.

## How this project works

These rules are enforced. Every one of them has caused a patch to be rejected.

- **Measured, not asserted.** A claim about geometry is backed by a test that
  measures it against a closed form, a published value, or an independent
  computation.
- **An approximation must say so.** Sampling, fitting and tolerances are fine
  when they are stated. A number that hides its own error is not.
- **Refusals are by name, and pinned.** When the kernel cannot do something,
  the error text says which thing and why, and a test holds it to that. A
  silent wrong answer is the only unacceptable outcome.
- **No stubs.** Either a thing is implemented, or it is refused by name and
  tracked as above.
- **Scope is parity with the reference kernel's modelling modules, plus what
  `docs/SCOPE.md` admits.** `docs/SCOPE.md` is normative and says how to decide
  a case. Parity is about capability, not structure.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against, bundles
  or imports another kernel, and the design is worked out here, not mirrored.
  The field's vocabulary is used throughout.

## Decisions, not gaps

These are settled. They are listed so nobody reopens them by accident.

- **A surface returns its point and derivatives in one evaluation, and states
  that the result agrees with the separate accessors only to rounding.**
  - A foot-point solve needs all six values at the same place. A
    tensor-product patch answering three separate accessors locates its spans,
    builds its basis functions and sums its control grid three times.
  - [`Surface::jet_at`] answers from one order-two table. This saves 3 to 25%
    of a spline-rich STEP read.
  - A patch sums its point by de Boor and its derivatives by basis functions.
    These reassociate differently, so a jet may differ from the separate
    accessors in the last ulp.
  - Consistency within a jet is what a Newton step needs, and what the type
    guarantees. A caller must not mix a jet with the accessors at the same
    parameters and expect identical bits.
- **A pcurve with no closed form is `None`, not a fit.** An exact curve with a
  fitted pcurve would have two descriptions that disagree by an amount nothing
  records. The consumer that needs one marches the pair instead.
- **A closed exact section partly outside a surface's extent is kept whole.**
  The restriction that matters is the face's trim, which is the boolean's own
  2D stage. The surface extent is only a parameterization window.
- **Scaled placements in the boolean are baked first.** A scale changes a
  surface's parameterization underneath its pcurves; `baked_shape` rebuilds
  the geometry, and the boolean calls it.
- **The crossing walker refuses tangential contact.** The tangential walker
  owns that case, and the section pipeline routes to it. The refusal stays
  pinned.
- **Bi-tangent construction is subsumed**: by the 2D repertoire in 2D, and by
  the blend family's own envelope in 3D.
- **Glue is subsumed** by the boolean's same-domain unification, which already
  skips nothing it needs and unifies what glue would.
