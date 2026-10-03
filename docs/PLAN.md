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
coincident. A facet leaning off a pad's wall by a ten-thousandth of a radian
or more, the two sharing an edge to their points' rounding, crosses the wall
on that edge rather than on the solved line a sliver off it. Leaning less,
the band where the two stand within the weld distance of each other is over
a tenth of the facet, and a slab drafted 0.00005 under a pad (below) still
reaches the refusal: the facet and the wall want pairing as same-domain.

**Blends between curved faces that share no edge.** `blend_faces` blends
curved faces that meet along edges of the solid through the edge blend (exact
or marched), curved faces that share no edge where their surfaces share a
direction or an axis, and otherwise marches the ball round where the
surfaces cross, capping the round in the ball's section where it leaves a
face, or through the end of the solid's edges along the crease where those
end first. `fillet_faces` marches between separate faces the same way, over
part of the seat or bridging a gap. A drum's foot leaning along a block's
edge closes at 1, 2 and 15 degrees with the axis on the edge or 0.5 past
it, as a face blend of the parted wall and as an edge round. With the axis
past the edge the cap stands in the ball's section through the crease's
end, which leans with the drum: one cap's wall corner stands past the
block's side and the other falls short, by 0.003 at 1 degree and 0.032 at
15. That is kept: every run-out caps in the ball's section, and a cap in
the side's plane would end the band off a section. Refused by name:

- Faces of a solid that the ball leaves at different places round its
  seat and never both at once: neither face gives the round both its ends.
- A marched round closing on itself between separate faces, and one whose
  line of contact crosses a face more than once.
- A round crossing a face that stands between the two (a ball smaller than a
  chamfer between the faces it rounds): the solid beside the middle of the
  round is open where the seat says material, or the reverse.

Not yet by name: an upright drum whose axis stands 0.5 past the block's
edge goes through the closed-form swept blend, which answers "does not
touch both faces" where the foot runs off the top.

**Section loops pinned at a sphere's pole.** Where a curved wall passes
through a pole of a sphere face, the marched section stalls and doubles
back at the pole, each loop starts and ends at the pole node of the chart,
and the arrangement can read one piece round both loops. Before the
boolean a whole ball is recharted about a clear axis, and so is a trimmed
sphere face of any solid (a dome, a dome on a drum, a cap): its trim edges
take fitted images in the new chart, its old seam and pole edges go, and
the new seam leaves the trim at a vertex. Drills through the pole of a
half ball, a cap and a half ball on a drum close at every wall offset
tried. Still pinned: a face whose trim the new chart would have to split
(two loops round the axis, as a band between two planes has, or a loop
the seam meridian crosses at no vertex), a face at a placement or on a
left-handed sphere, and an edge with no curve in space.

**A tilted half ball fused with a drum.** A half ball whose sphere is
charted about a tilted axis, fused with a drum, fails with "arrangement
left no piece of the face": an error, not a wrong answer, seen while
building the dome drills. Not yet diagnosed.

**A section whose fit misses widens its junctions as far.** A marched
section keeps its fit's error as its tolerance, and the paving's reach and
the junctions at the ends of a stretch it runs along an edge grow with that
tolerance, unbounded. The error is measured in space: the curve against
the trace, each pcurve lifted against the curve, and where the surfaces
cross at a shallow angle their gap over the angle's sine, capped by the
chart residual carried through the surface's stretch. A section across a
wide drum now states microns where it stated tenths, a branch missing by a
hundred times its tolerance is fitted in pieces, and a fitted section meets
an edge within its own reach plus the edge's radius. Two sources remain:

- A section through or beside a sphere's pole fits its samples but its
  sphere pcurve strays up to 0.09 between them. Split at the pole sample
  each piece measures under 1e-3, but the boolean then refuses a drill it
  passed whole (the kept pieces do not close), so the split waits on the
  pole arrangement. Recharting takes the pole away wherever it applies.
- A branch nothing fits (residuals of 0.1 to 10 near a tangency) still
  carries its miss, as does a near-tangent crossing whose capped angle term
  stands at 1e-3 or more.

**Facet walls drafted less than their chords sag.** A pad pushed down from
the top of a slab converted face for facet, whose rounded corner's rows lean
in by less than the corner's chords sag, fuses, cuts and is taken in common
at every row count tried when the slab is drafted 0.001 to 0.005 by the
rows' depth to the power 1.5 and each row stays a face of its own, and
also at the converter's own coplanar distance: the scatter estimate
compares each nearly flat edge's turn with the edges beside it and
continuing it, and a steady turn row after row is read as a curve drawn
finely, not scatter (84 of 84 row-count sets measured at draft 0.001,
against 15 before). Drafted 0.0001 to 0.0007 the top band leans off the
pad's walls by under a thousandth of a radian and crosses them on the top
edge both faces hold: of 17 row-count sets at those 4 drafts and 2 depths,
133 of 136 pass with each row a face of its own and 115 of 136 at the
converter's own distance (79 and 67 before). What remains:

- Drafted 0.00005 or less the band where a top facet stands within the
  weld distance of its wall is over a tenth of the facet, and the facet's
  edge is not taken for the crossing: the kept pieces are left open, or a
  piece is refused as an unpaired coincidence. The pair wants same-domain
  pairing.
- At the converter's own distance, rows that lean alike merge into one
  face on a fitted plane, its vertices up to 6e-4 off it within their
  tolerance. A pad whose floor crosses such a face came out off the mesh's
  volume by 2e-5 to 5e-4 (8 sets at a depth of 3). Not a boolean error:
  slicing the slab with a plain box at the floor moves cut plus common by
  the same amount to three digits, the pieces of a face tessellated on
  vertices off its plane differ from the face tessellated whole, and every
  case stays under a fiftieth of the crossed faces' doubt (vertex offset
  times area, 6e-3 to 4e-2); where the crossed faces hold their vertices
  the identities hold to 1e-12. The test's tolerance now adds that doubt.
  The converter merging rows that far off a plane is the cause to take up
  if exact volumes on such parts matter. Of the 15 committed row-count
  sets at 4 drafts and 2 depths, 107 of 120 pass at this distance; of the
  rest 10 are refused as kept pieces that do not close (9 at a depth of
  10), 2 as a fuse whose rays all meet a degenerate point, and one fuse
  is valid but its tessellation open (three to nine a row, drafted
  0.0003, depth 10).
- Each row a face of its own, three sets at a depth of 10 leave the kept
  pieces open or the fuse's tessellation open (eight to five or three to
  nine a row, drafted 0.0003 to 0.0005).

**Mesh conversion.** `solid_from_mesh` rebuilds planes, the four canonical
surfaces, extrusions, surfaces of revolution and fitted B-spline patches,
and leaves every other region faceted. Today it fits each region on its
own, solves each seam point by point along the mesh boundary between two
regions, and takes its corners at mesh vertices. The failures left come from those three steps:

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
5. **A fitted B-spline patch** for a smooth region nothing else fits.
   Done (`MeshSolidOptions::patches`, on by default): a region that is one
   disk with one loop, with vertices enough inside it, is charted by a
   canonical surface within ten times the distance or else by the
   mean-value map onto a square (one banded solve in reverse Cuthill-McKee
   order, no new dependency), fitted by `fit_surface_scattered_at` with a
   thin-plate term, knots split where feet stand off and parameters
   corrected, verified both ways, continued past its square and handed to
   the culprit loop. Refusals are counted (`patches_not_disk`,
   `patches_narrow`, `patches_unverified`). Measured: the corpus and the
   truth bench come back unchanged; no region there verifies (Body28 has
   13 unverified, the rest are not disks or are rows of facets).
   On a fine mesh, recognition cuts small spheres, cones, cylinders and
   tori out of a free-form surface, and the free triangles left round them
   are no disk. The pieces the smooth area encloses (curved regions and
   pockets of fewer than 32 free triangles that meet nothing outside it
   except across creases) are now tried with the smooth regions they join
   as one region, and kept apart where that patch does not verify. A round
   tangent to a face outside the smooth area is never taken. Measured on
   the bumped plate meshed at 80 and 120 cells across: before, 12140
   faces (4 curved, all faceted) and 24153 faces (70 curved, 46 faceted),
   one region not a disk each; now 6 faces, the top one patch, the exact
   surface within 6e-6 of it. The corpus, the NIST parts and the truth
   bench come back unchanged (faces, kinds, volumes, refusal counts).
   Still to do:
   - the fit holds vertices to half the distance, but between rows near a
     free edge it can stand off the true surface by about twice the
     distance on coarse meshes;
   - a patch meeting a neighbour tangentially falls to the chord seam or
     to facets.

   Tangency to a canonical neighbour comes after item 2.
6. **The steps as an API.** Done: `MeshRegions::find` gathers the
   regions (each with its triangles, surface, measured deviation and
   neighbours), `merge`, `split` along a vertex path and `fit` of a
   `SurfaceKind` with `FitConstraints` (a fixed axis, or a cylinder's,
   sphere's or torus tube's radius) change them, and `build` makes the
   solid. `solid_from_mesh` is `find` then `build`, so both go through one
   path; the corpus and the truth bench come back bitwise unchanged. Every
   surface a step puts on a region is verified at every vertex and by the
   triangles' sag, and a step that cannot be taken is a named
   `RegionRefusal`. Tests (`mesh_steps`): a cylinder drawn with 36 degree
   facets, faceted by the automatic conversion (volume 5877.85, the
   decagonal prism's), merged into one region builds the cylinder
   (6283.185307179579 against pi r^2 h = 6283.185307179587), the same
   solid the automatic conversion builds with the crease raised; a roof
   two hundredths high found as one plane builds 200.22 (a 3.9e-3 error),
   and split along its ridge 201.000000000001 of 201. Not done: an edited
   region keeps the surface it was given, and the tangent passes that
   recognition runs over its regions (rounds between planes, blends,
   corner balls) do not see it; a merged or split region is refitted as a
   plane, a canonical surface or a sweep, never as a patch, and `fit` has
   no sweep or patch kind. Item 4's corner and boundary graph, once there,
   is what an application would edit next.
7. **A torus pierced across its outer equator.** Done where a parallel
   is free: a whole torus with holes has its seam round the axis placed
   on the parallel the holes leave widest free, as its seam round the
   tube already was on the widest free meridian, and the holes are inner
   wires of it. Measured on a ring of radii 10 and 3 drilled with a bore
   of radius 1 blind through the outer equator, blind through the inner,
   through both, and blind through the top parallel: each comes back one
   toroidal face, valid, meshing closed, its exact volume within 4e-9 of
   the original's. The outer equator case first measured 5e-4 off: the
   exact integrator told the boundary from the holes by a single Gauss
   rule across each closed spline, which gave a hole the wrong winding;
   it now integrates knot span by knot span. Still open: holes that
   between them cross every parallel (a slot running all the way round
   the tube, or a ring of holes staggered round it) leave no seam
   position, and the tube stays facets. Built whole, a hole would have to
   join the seam: the outer wire running along the parallel to the hole,
   round it and on.
8. **The free boundary of an open mesh.** Done: a curved face's free
   boundary is cut where it turns by the crease angle and each run placed
   on the face's surface, as a parallel or ruling where one holds and
   otherwise as a curve fitted through its vertices' feet
   (`free_edges_fitted`). Measured on open tubes, sheets, a dome and a
   free-form sheet, which all came back faceted before. A mesh with no
   flat stretch no longer has its coplanar distance read off its gently
   curved triangle pairs: a bumped sheet alone measures 1e-4 (2e-3
   before; the closed plate 3e-5) and comes back one fitted patch. Where
   the turns change sign the estimate asks them to change linearly, so
   a coarse mesh of a surface whose curvature varies fast can still read
   its inflections as a little scatter.

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
- Coincidence of two patches: done. `project_on_surface` lives in
  `ogeom-geom`, and `intersect_surfaces` measures a pair without a closed
  form for coincidence and answers `Same` for every caller. The boolean
  still asks `coincide_as_stated` itself, since only the faces know their
  stated tolerance. No measurable time change on the stress run.
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
