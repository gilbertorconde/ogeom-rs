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

**Unpaired coincidence in the boolean.** Contact confined to lines or points
is classified off the contact. A piece that still reads on the other solid's
boundary at every point it is asked at, with no coincident partner face, is
refused: that is a contact over a region of faces the boolean did not pair as
coincident. No configuration known to reach it has been found; one that does
wants its two faces paired as same-domain.

**Blends between curved faces that share no edge.** `blend_faces` blends
curved faces that meet along edges of the solid through the edge blend (exact
or marched), and curved faces that share no edge where their surfaces share a
direction or an axis. Two gaps are refused by name:

- Faces sharing no edge whose surfaces share no direction or axis (a B-spline
  face, a cylinder at a slant to a plane). These want the march along the
  spine where the two offset surfaces meet, with no edge of the solid to
  guide it: the surfaces' own intersection as the guide.
- A round crossing a face that stands between the two (a ball smaller than a
  chamfer between the faces it rounds): the solid beside the middle of the
  round is open where the seat says material, or the reverse.

**Two drills in `nist_ftc_06` are refused.** Replay with
`ogeom-stress --case drill/nist_ftc_06_asme1_rd/<n>`; a refusal prints the
drill's placement. Both refuse the same way on the part with or without its
hemisphere's great circle split at the poles.

- Drill 0 ("a wire needs at least one edge"): the drill's wall hugs a face
  of the part, and the hug's junctions reach 72 mm. Merged, one junction
  reaches 148 mm, every edge of a piece of face 0 closes on its vertex as
  dust, and the piece is left with an empty ring.
- Drill 3 ("the kept pieces did not close into a shell"): not diagnosed.

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
2. **Fillets from their supports.** Done for the closed forms: rounds
   between two planes (already), a torus between a plane and a coaxial
   cylinder or cone, and corner balls, whose rounds now share one radius
   so all their axes meet at the ball's centre. Not done: a torus between
   two cylinders or a cylinder and a sphere, variable radius, and fillets
   between curved faces off a common axis, which keep their fitted
   surfaces. Below the data's noise the regions broke up before any fillet was
   recognized; such a distance is now raised to the noise (C3, closed).
3. **Loops checked in each face's chart.** Done, after the build rather
   than before it: a fold deeper than its two edges' tolerance, outside
   every corner's, and still there drawn a hundred times finer names its
   two seams, and the curved faces across them are withdrawn. Measured on
   the corpus before the check, final outputs held folds in Body28 (174
   faces, 7 of them past their edges' tolerance), 77777_1 (18), shelf_bracket
   (6) and handle-pickup (2); all but Body28's deep ones lie within their
   edges' or corners' tolerance, which a tolerant boundary allows. Body28
   gives up 25 more curved faces and still tessellates open, so its open
   seams are not these folds. A guard on the point-by-point section
   (refusing a curve that runs against its chain) was tried and dropped:
   it left Body28 invalid and moved 77777_1's volume away from its mesh.
4. **Topology first.** Measured 2026-10-02 at 093c0de (benches and
   counts kept with the scratch benches): of the corpus parts, Body28 and
   handle-pickup tessellate open, and of the NIST parts meshed at a
   thousandth of their diagonal, ctc_02, ftc_07 and ftc_10. Classified by
   the faces whose meshes leave an edge used other than twice:

   | Part | Bad mesh edges | At a planar facet of 3 or 4 edges | Elsewhere |
   |---|---|---|---|
   | Body28 | 436 | nearly all | a few on tori and spheres |
   | handle-pickup | 5 | 3 | 2 on a sphere |
   | ftc_07 | 212 | 69 | spheres 53, larger planes 62, cylinders 21 |
   | ctc_02 | 11 | 4 | tori 7 |

   Seams today: most are exact closed forms or point-by-point sections;
   the chord through the chain's vertices carries 868 (Body11), 1581
   (Body28) and 736 (77777_1). The work goes in stages, each measured on
   these counts, the corpus and the truth bench before the next:

   - **4a. Sliver facets at junctions.** A planar face of one or two
     triangles between curved faces, narrower than the tolerance of the
     seams that bound it, folds inside that tolerance and meshes over
     itself. Its triangles go to the curved neighbour across its widest
     seam when every vertex of it lies within that seam's tolerance of the
     neighbour's surface, the neighbour's tolerance raised to cover them;
     a facet within the tolerance of its corner collapses into the corner.
     Either way it changes nothing past what its seams already claim.
     Target: the first class in the table.

     Tried 2026-10-02 (patch kept with the scratch benches). Facets
     narrower than the reach are few: absorbing only those moves Body28
     from 436 bad mesh edges to 419 and ftc_07 from 212 to 209. The
     facets at the bad edges are ordinary triangles no surface claimed;
     their seams, often chords whose tolerance runs to a twentieth of
     their span, fold across the whole facet. Absorbing every facet of one
     or two triangles within the reach of a curved neighbour, whatever its
     width, does what the stage was for (Body28 to 68 bad mesh edges,
     ctc_02 closed, ftc_07 to 166, the truth bench's face ratio from 0.27
     to 0.17, Body11 from 5258 faces to 3621) and breaks what was exact: a
     facet on a cone's rim taken into the cone turns its exact rim circles
     into sections and the cone is withdrawn (boss_cone_chamfer, 9 faces to
     56), and shelf_bracket and 77777_1 lose their tori. It also tips
     Body28 invalid through a latent planar face of almost no area (three
     of its corners on a line) in a faceted fan far from any facet taken;
     `check`'s inside-out probe flags that face, and runs in about a second
     on the largest parts. Open: which facets to take without costing an
     exact seam. The zero-area planar faces are gone: a triangle no higher
     than the coplanar distance joins a neighbouring plane by its corners'
     distance (Body28 alone moves, from 436 bad mesh edges to 179).
     Absorbing every facet within the reach now gives a region's facets
     back before it is faceted, and makes the conversion again without
     them if a face ends up turned in: the cone and the tori stay exact,
     Body11 goes from 5258 faces to 3612 and linkage_bores_chamfer from 397
     to 223 with nothing faceted, Body28 to 122 bad mesh edges, and the
     truth bench's face ratio from 0.27 to 0.17 (noisy window-tube and
     rim-disc lose a little: 3.2e-6 to 6.1e-6, 7.2e-4 to 7.8e-4). Still
     open: handle-pickup's folding facet (its corners stand beyond the
     reach of both neighbours), Body28's remaining bad edges, and ctc_02
     and ftc_07, which go from 11 and 212 bad mesh edges to 23 and 223:
     some facets they take make worse seams than they removed.
   - **4b. Corners.** Measured: of the bad mesh edges, most sit at a
     vertex whose edges' curves end off it, each on its own side (gaps up
     to 6e-2, inside the vertex's tolerance). Solving the corners onto
     their surfaces in the converter was tried and dropped: four fitted
     surfaces seldom meet in a point, and moving the corners helped some
     parts and hurt others (ftc_07 from 221 bad edges to 307). The
     tessellator now draws each edge from its vertices, in space and in
     each face's chart, which serves imported solids too: 77777_1 from 64
     bad mesh edges to none, Body28 from 122 to 31, ctc_02 from 25 to 6,
     handle-pickup from 5 to 3, ftc_07 from 223 to 188.
   - **Closed by withdrawal.** What 4a and 4b left (long strips a row of
     facets wide whose sides cross inside the chord tolerance they claim,
     and a torus built with its seam through a hole) is caught where the
     solid meshes open: the curved faces whose own meshes do not meet
     their neighbours' fall back to facets. Every part measured now
     tessellates closed. The cost is those faces (ftc_07 gives up 20 of
     143) and time on parts that need it (ctc_02 from 27 s to 50 s).
   - **4c. Seams traced between solved corners**, as one branch of the
     intersection guided by the mesh path, replacing the point-by-point
     section where the two surfaces cross. Only if 4a and 4b leave folds,
     and nist_ctc_05's folded planes are the case to measure it on.
   - **4d. Junction rules** (Benière et al. 2012): a junction stays only
     if its faces are pairwise adjacent; dangling edges are dropped.
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
7. **A torus pierced across its outer equator.** The whole torus's seam
   runs round the outer equator, so a hole there cannot be an inner wire
   of it; the tube stays facets. Built whole, the hole would have to join
   the seam: the outer wire running along the equator to the hole, round
   it and on.
8. **The free boundary of an open mesh.** A recognized region's edge with no
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
