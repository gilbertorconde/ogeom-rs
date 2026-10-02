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
the corner, the near-coplanar sliver band again. A slab drafted 0.005 or
0.002 by the rows' depth to the power 1.5 passes at every row count tried;
at 0.001 some counts refuse, the open edges being the pad's corner walls
where they run along the corner's facets.

**Mesh conversion.** `solid_from_mesh` rebuilds planes, the four canonical
surfaces, extrusions and surfaces of revolution, and leaves every other
region faceted. Today it fits each region on its own, solves each seam point
by point along the mesh boundary between two regions, and takes its corners
at mesh vertices. The failures left come from those three steps:

- a fillet fitted apart from its supports meets them at a near tangency the
  seam solve cannot settle (C3 in `docs/REVIEW.md`);
- a seam that follows a jagged boundary crosses back over a small face and
  folds its trim, and the faces caught are found after the build and
  withdrawn whole;
- a whole turn hides a tilt of the axis from the profile fit.

Patching these case by case has stopped paying. The items below change the
steps themselves, in this order. Each is held to the converter's standard
(verified against every sample, edges placed on both surfaces, facets as the
fallback) and measured on `mesh_corpus` and the truth bench before and
after.

1. **Sweeps all the way round.** Done for a whole turn: the axis read
   from the normals (already the pitch-zero line complex) is brought close
   by the profile search and settled by Gauss-Newton on the samples, which
   lie on the surface where the normals only lean toward it. A sweep that
   goes round is laid out as a wrapped face, its seam seated between its
   rims. Still to do: the band layout for a sweep between two circles (its
   seam the profile itself), and a closed profile (an extrusion of a closed
   curve, a surface of revolution of a closed one). The full helical line
   complex would tell a thread from a surface of revolution; the mesh's
   own normals carry a small twist (a pitch of about 1e-3 on a turned
   wave meshed at 0.01), so its threshold has to sit above that. A thread
   goes to the patch (item 5).
2. **Fillets from their supports.** A constant-radius fillet between two
   recognized faces is built from them and a radius (Kos, Martin, Varady
   2000), not fitted freely. The rolling ball's centre runs where the two
   supports, offset by the radius, meet; the fillet is the cylinder or
   torus about that spine; its edges are the contact curves, in closed
   form. The tangency is exact by construction and no near-tangent seam is
   solved. Only the radius is fitted, and the fillet is verified against
   every sample of its region. Closed forms: two planes give a cylinder; a
   plane with a cylinder square to it, two coaxial cylinders, and a plane
   with a cone on its axis give a torus. A corner ball where three equal
   fillets meet is the sphere about the point their spines share. Other
   pairs keep today's path. Seating a corner sphere on its fillets alone was
   tried and opened two more slits: the seams also have to come from the
   corners (item 4), not from the mesh boundary.
3. **Loops checked in each face's chart.** Before the build, each face's
   loop, mapped into its surface's parameters, must be simple and turn the
   right way. A failure names the seam, so the culprit loop withdraws the
   faces of that seam instead of finding them afterwards by overlap and
   volume. This also measures how many folds item 4 has to remove.
4. **Topology first.** Before any seam is solved, the boundary graph on the
   mesh (corners where three or more regions meet, boundary paths between
   them) is cleaned:
   - corners closer than the tolerance are merged;
   - a facet region of one or two triangles at a junction goes to a
     neighbour that verifies it, or into the corner;
   - every face is bounded by simple cycles;
   - a junction stays only if its faces are pairwise adjacent (Benière et
     al. 2012).

   Then each corner is solved onto its surfaces (least squares with a
   residual check where four or more meet, refused by name otherwise), and
   each seam is traced as one branch of the intersection from corner to
   corner, guided by the mesh path, instead of point by point. A branch
   traced between fixed ends cannot fold. This is the general answer to the
   folded trims the crossed-seam, straight-seam and turned-over checks catch
   now, and to `nist_ctc_05`: two large planes and a row of small ones whose
   neighbours' seams cross back over them, which faceting their curved
   neighbours only shrinks. It replaces how the plan builds edges, so it
   starts with a design note and a corpus baseline.
5. **A fitted B-spline patch** for a smooth region nothing else fits:
   - a region that is not one disk with one loop stays faceted;
   - the chart comes from a canonical surface that nearly fits (within about
     ten times the tolerance) where there is one;
   - otherwise from a mean-value map onto a square (Floater 2003): one
     sparse linear solve, without folds when the boundary is convex, the
     corners where the neighbouring face changes. The sparse solver is
     chosen when the item starts;
   - `fit_surface_scattered` takes those parameters, with a fairing term and
     knots inserted at the worst span, then two or three rounds of
     parameter correction;
   - the patch is verified both ways, triangle interiors included (allowed
     the tolerance plus the triangle's own sag), and its Jacobian and
     normals checked.

   Tangency to a canonical neighbour comes after item 2.
6. **The steps as an API:** the regions found, merging two, splitting one
   along a vertex path, fitting a chosen surface type to one (optionally
   with a fixed axis or radius), then building the solid. The pieces exist
   inside `solid_from_mesh`; the work is a stable surface for them and tests
   that a corrected conversion builds what an automatic one would have. It
   follows item 4, whose corner and boundary graph is what an application
   edits.
7. **The free boundary of an open mesh.** A recognized region's edge with no
   neighbour has no second surface to be solved onto. Check whether it comes
   back as a polyline. If so, place it on the region's surface as a curve
   (a circle on a cylinder's rim, a fitted spline otherwise), split where it
   turns sharply.

**Faces whose trims fold in their chart tessellate open.** Some converted
solids (`nist_ftc_07` in `mesh_corpus`, `nist_ctc_02` meshed at a
thousandth of its diagonal) pass `check` and tessellate open: a face's pcurves stray from its 3D edges by more than the
tessellation's tolerance, or fold in its chart. Most should go with item 4
above; what is left after it is measured first. Then, per face: where a
pcurve strays, the boundary's chart positions are re-derived by projecting
the shared 3D edge points onto the surface, the edges staying as they are;
the loops are checked before meshing and the triangles' normals against the
surface after. Each face that needed it is named with its deviation, since
an open tessellation is the signal that a solid is suspect and must not be
hidden. A face that still fails could be meshed by an advancing front on
the surface itself, which needs no chart.

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
  - Fillets in few passes. A chain is taken in rounds, no two edges of a
    round sharing a vertex, and each round's blends are applied in one
    boolean each way, their wedges taken as a compound of disjoint solids.
    Each edge still builds its blend against what the earlier rounds left.
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
- **The converter uses no learned model.** Published learned reconstruction
  (2024 to 2026) gives a valid solid for 70 to 76% of parts under a hundred
  faces, with every face a B-spline, and places its own bottleneck in
  rebuilding the b-rep, which is the converter's verified build. A learned
  proposal is also a result nothing here measures. Global selection by an
  integer program is out for the same kernel reason: its time has no bound.
- **Glue is subsumed** by the boolean's same-domain unification, which already
  skips nothing it needs and unifies what glue would.
