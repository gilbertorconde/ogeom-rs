# Plan

What is left to build, and the rules the work follows.

## Where the remaining work is tracked

- **GitHub issues** hold bugs and requests, one per issue.
- **The parity ledger** (`docs/PARITY.md`, generated from
  `docs/parity/parity.toml`) audits the reference kernel's modelling modules
  against named capabilities. The restriction on each `partial` row is the
  list of what that capability does not do yet. `tools/check.sh` holds the
  ledger to its evidence.
- **This file** holds the open items that are neither: known limits of the
  code, refusals by name, design items and speed work.

Finished work is recorded in `CHANGELOG.md` and in the tests that pin it, not
here.

## Open

Within each group the items are ordered by risk: a possible wrong answer
first, then refusals, then cost.

### Checking and measurement

**Meshes of doubly curved faces stand off by up to twice the chord.** The
grid holds each row and each column to the chord, and a cell's middle
stands off by both sags together; the repair pass fires only at three
chords. A ball of radius 20 meshed at chord 0.1 stands off 0.17, at 0.5
0.68, and says so (`deflection_met` false). A cone's face, scaled by its
mean radius for Delaunay, joins points three columns apart where its
radius is far from the mean: the cone of #202 stands off 0.08 at chord
0.0125. What was tried (2026-10-12): holding each direction to half the
chord on balls and tori, with the repair asking each triangle's middle
and splitting its longest interior edge, puts balls within the chord at
every chord asked and passes the mesh conversion corpus, but the
inserted points break a ball's symmetry (`a_sphere_converges_on_the_volume_a_sphere_has`
reads its centre off the origin) and three `mesh_solid` conversions;
the repair at one chord on every edge, or at the middle alone, leaves
converted faces whose meshes do not meet (217 edges on `nist_ctc_02`) or
a button head of 57 faces. A grid fine enough by construction (rows and
columns from the surface's curvature both ways, the cone's chart scaled
per row) keeps the symmetry and the converter's inputs regular, and is
the next thing to try. Closes the rest of #202.

**Some fits through bare points are trusted between their points.** A fit
through points (`fit_points`, `fit_points_at`, the grid fits) states its
error at the points only. These callers still read it, or nothing, where
the shape between the points matters. Every edge the workspace's tests
build through them was measured densely against its stated tolerance, and
none stated less than its pcurves stand off:

- `ogeom-offset/src/sweep.rs` (`Skin::fit`, `skinned_wall`): a grid skin's
  error becomes the slack an adopted border is widened by, and
  `adopt_border` reads that border against the surface at 33 points.
  Between the sections the skin is the fit itself; along a section the
  section curve is the truth and is not read between the grid's points.
- `ogeom-fillet/src/pair_marched.rs` (`contact_cut`): the ball's line of
  contact is held to its target at the march's stations only; its pcurve
  is measured densely against it, but the true contact line between
  stations would take a march to read.
- `ogeom-hlr/src/exact.rs`: a silhouette is fitted through the walk's
  points and kept whether or not the fit met the chord. The fit's error
  bounds nothing there; the line is drawn.

Read and left: `ogeom-algo/src/convert.rs` (`fit_pcurve` stops refining
at 384 samples and checks middles only, but every re-derived edge is then
measured densely by `state_pcurve_gaps_of`), and
`ogeom-algo/src/recognize_swept.rs` (a mesh's vertices are all the truth
there is, and the profile is measured against every one). The edges an
offset rebuilds, a sheet offset or thickened makes, `split_face` makes and
a blend's wedge carries are measured densely the same way.

### Booleans

**Tools against a box rounded on every edge.** The stress `rounded`
scenario cuts a slanted box rounded on all twelve edges by an oblique
primitive and takes the common, the two adding up to Steiner's volume. At
seeds 1 to 4, 50 cases each, none is wrong, and two refuse once #199 is in:

- `--seed 3 --case rounded/primitives/118`: a cone through a box rounded
  at 0.62 refuses with "edge 0 ends where edge 1 does not begin", two
  vertices 1.07e-4 apart where the cone's sections cross the tangent seam
  between a box face and a band. The plane and the cylinder meet
  tangentially there, so the cone's two sections reach the seam's line
  that far apart. Exact conics (#200) leave it as it was.
- `--seed 2 --case rounded/primitives/134`: a drum of radius 3 through a
  box rounded at 1.39 refuses with "the kept pieces did not close into a
  shell", five open edges round a corner ball and its bands, as though a
  piece of one face were dropped.

**Nearly coplanar faces are not paired as coincident.** A piece that reads on
the other solid's boundary at every point asked, with no coincident partner
face, is refused by the boolean's piece classification ("a contact over a
region of faces not recognised as coincident is refused"). The case is a
faceted wall leaning off a pad's wall by under a ten-thousandth of a radian:
the band where the two stand within the weld distance covers most of the
facet. Two planar faces where one lies on one side of the other's plane, its
corners past it by rounding at most, no longer give their planes' line as a
section (`one_side_of` in `ogeom-bool/src/lib.rs`): a corner facet meeting
the wall only at a vertex, and a facet's side edge beside the next wall. Of
the 15 committed row-count sets at depths 3 and 10, all 30 pass with each
row a face of its own at 0.00001, 0.00002 and 0.00005; at the converter's
own distance 25, 30 and 27 of 30 pass. What remains:

- At 0.00001, [7,7,7,6,6,6,5] and [8,5,9,9,5,8,7] at the converter's
  distance, both depths: the kept pieces do not close. In the first, the
  slab's merged corner face meets the pad's flat wall in a section that
  leaves the corner edge 7.8 mm up and reaches the pad's floor two microns
  off it; the sliver between that section and the pad's corner edge is
  kept on one side and not the other, and the corner edge above the
  section's end is left open.
- At the converter's distance with depth 10 (the pad through the whole
  slab), the fuse comes out 4e-5 to 6e-5 cubic millimetres under the mesh's
  measure: 0.00001 [8;7] and three sets at 0.00005. Within the weld distance
  times the top band's area, past the test's 2e-5.
- The design for a band the two stand within the weld distance of each
  other over a region: pair it as same-domain and split it where the planes
  part by the weld distance.

**Section loops pinned at a sphere's pole.** Where a curved wall passes
through a pole of a sphere face, the marched section stalls and doubles back
at the pole, and the arrangement can read one piece round both loops. Before
the boolean, `ball_chart.rs` restates a whole ball or a trimmed sphere face on
a chart whose poles stand clear; a solid at a rigid placement is restated in
its own frame, and one scaled, mirrored or on a left-handed sphere is baked
first. It leaves as it is (and the pole stays): a face whose trim the new
chart would have to split (two loops round the axis, or a loop the seam
meridian crosses at no vertex), a face placed inside its solid, and a face
with an edge that has no curve in space. A ball cut to the slab |x| < 0.8 or
0.5 (radius 10) and drilled along z through its pole fails so: every clear
axis either lies within a tenth of the radius of the trim or has both loops
round it.

**A section whose fit misses widens its junctions as far.** A marched section
keeps its fit's measured error as its tolerance, and the paving's reach and
the junctions along an edge grow with it, unbounded. Two sources remain. A
section through or beside a sphere's pole fits its samples, but its sphere
pcurve strays up to 0.09 between them; split at the pole sample each piece
measures under 1e-3, but the boolean then refuses a drill it passes whole, so
the split waits on the pole item above. A branch nothing fits (residuals of
0.1 to 10 near a tangency), and a near-tangent crossing whose capped angle
term stands at 1e-3 or more, carry their miss.

**A drill lying in a face that runs into a round.** A drill whose lowest line
lies in a plate's underside meets the round at each end within a few microns
of tangency. The section is a figure eight crossing the round's tangent edge
twice at its double point, or two loops crossing it a hundredth of a
millimetre apart. With the drill's seam toward -y, the marcher returns the
section as six open pieces whose fits miss by up to 9e-3, and the kept pieces
do not close; the same placement moved by a few ulps leaves the drill less
the part unclosed. With the seam toward z it closes. The stress drill
scenario's seed 2 drill 5 on `nist_ctc_03` is the same: the drill's top
ruling lies in the plate's top face and runs into the round of radius 11 at
its edge, where the drill, the round and the plate's face are tangent at one
point; the marched drill-round section, fitted to 2.8e-4, crosses the round's
tangent edge twice 9.5e-3 apart, and the sliver between the crossings sews
to nothing (its twin edges stand 9.2e-4 apart).

**A drill's seam along the other part's edge.** A drill whose wall runs along
a round's tangent line with its own seam on that line: the section hugs an
edge on both sides and lays no split, and the boolean refuses. A drill tangent
to a plane along a face's edge, seam there, is the same. Carrying the seam
onto that face as a contact closes it, but the same contact breaks a drill
whose seam crosses a face's interior lying in its plane. Turned away, the
seam closes. The concave case: seed 3's drill 0 on `nist_ctc_03` has its seam
along the edge where a wall (y = 425.455) meets the lip's ceiling above it,
with the drill inside the material on both sides of the line, so the edge
turns nothing and is rightly carried as no contact (the touch along a
boundary edge is carried only where the solid turns across the touching
face). The drill's wall below the seam is then never split by the faces it
crosses behind the wall, and its band inside the material merges into the
outside piece, whose one ring runs along both sides of the seam.
A drill whose seam stands up to a few tenths of a micron above a block's top
face cuts through validly; a blind one in the same place is refused (the
piece where it ends lies on the block at every probe), and a drill with its
top line in a plate's top face where the face runs into a round is refused
for some seam directions and offsets of up to a tenth of a micron (radii 11
and 6, eight of 45 placements). Treating a plane and a drum that part by no
more than the confusion distance as tangent closes the blind block, but the
drum's section against the round still sees the drum above the plate, and
the pieces stop agreeing: twice as many of those placements refused.

**A section tangent to an edge at its end.** A drill whose circle on a rounded
box's base touches the round's tangent line at the corner vertex, where the
corner ball touches the base, is refused. At seeds 1 to 3 the stress drill
scenario refuses one drill each: seed 1's drill 2 on the rounded box is this
item, and seeds 2 and 3 each refuse one on `nist_ctc_03`, recorded under the
two items above. A drill along y whose top line touches `nist_ftc_11`'s bore
rim where the bore meets the torus above it (stress seed 11, drill 14) has
its section along the top stop a millimetre short of the touch; the face it
leaves does not close in its chart, and the boolean refuses it.

**A loop's pcurve stopped short of its seam column after a fuse.** A
section that crosses a torus's seam a hair from its own start is marched
as a long piece and one a few hundred-thousandths long across the seam;
the fuse drops the short piece, the loop closes in space within its
vertex's radius, and on the fused face its pcurve stops 1e-4 short of the
column it closes on. A later boolean on the fused solid meets that gap:
past its arrangement's node snap, the strand ending on the column dangles
and it refuses. The fuse should carry the short piece into the loop's
pcurve, or state the loop on the column it closes on. Of the chain stress
scenario's 900 cases at seeds 4, 9 and 10, 1, 2 and 1 are refused (3, 5
and 3 before the seam fixes of 2026-10-12): three marched fits beside a
sphere's pole ("a branch passing beside a chart's pole fits nothing
yet"), and seed 4's `chain/primitives/178`, two fused balls cut by a
cone, whose boundary strand dangles.

**A tilted half ball drilled on the diagonal chart.** A half ball charted
about the cube diagonal (1, 1, 1) with its seam toward z, on a drum and drilled
through, gives fuse + common and cut + common 1.2e-8 relative off the
operands, where other charts hold 5e-9. The rim and section are fitted on the
tilted chart (held to 7e-6 and 2e-5), and the measured volumes move by that
much.

### Fillets and blends

**Curved face blends refused by name.** `blend_faces` and `fillet_faces`
round curved faces that share no edge where their surfaces share a direction
or an axis, and otherwise march the ball round where the surfaces cross. They
refuse (`curved_face_blends_refuse_by_name` in `tests/blends.rs`):

- faces the ball leaves at different places round its seat and never both at
  once, so neither face gives the round both its ends;
- a marched round closing on itself between separate faces, and one whose
  line of contact crosses a face more than once;
- a round crossing a face that stands between the two (a ball crossing a
  chamfer between the faces it rounds), where the solid beside the round is
  open where the seat says material, or the reverse.

**A pole leg whose edge winds round the pole and whose rail does not.** A
marched band whose rail passes over a sphere's pole has a leg from the
winding loop to the pole with the other loop cut out as a hole. When the
winding loop is the rail, `pole_leg` in `marched.rs` splits it clear of the
hole and starts the seam there. When it is the edge, the seam leaves from the
edge's own vertex and can cross the hole. No measured part builds one.

### Construction and exchange

**The converter's result depends on how the mesh lists its triangles.**
Recognized regions grow from seeds taken in list order (a stride through
the triangle indices) and claim triangles first come, and the default
distance widens when any one fit presses against it. So the same mesh
listed otherwise finds other regions: on `nist_ftc_07` 378 to 666 faces
over 30 rotations (the distance widened in about half), and on every
corpus part some shuffled listing gives a different face count, a region
faceted (`grid_point_on_a_diagonal_boundary`, `nist_ftc_10`). All 30
rotations of `nist_ftc_07` meet, against 9, at 349 to 614 faces, one with a
region faceted; all 30 of `nist_ftc_10` meet, at 282 to 369 faces, three
with one to three regions faceted. Putting
the mesh in an order read off the geometry (points by position, triangles
by their corners; patch `canonical-order.patch` in the scratch crate)
makes every listing one result, but the one it picks facets two regions
of `nist_ftc_08` (a small cone at a cylinder's end) and six of
`nist_ftc_10` (one small cylinder's fit presses, the distance widens
twice, 1272 faces): it is a listing like any other. Before it lands, the
growth has to stop depending on who comes first: options are widening
only where many fits press (noise is everywhere, a run-on is local),
seeding from the most confident triangles first, or letting regions that
meet compete for their shared triangles by residual.

**STEP pcurves the file cannot state exactly.** The reader keeps a file's
pcurve where its lifted gap is within twice the projection's. A pcurve the
file cannot state in the curve's parameter (a reversed or left-handed conic,
a trig or offset trace, a curve written as a spline conversion) is derived on
reading. A foreign file's pcurves are read without its length or angle unit,
so one in inches or degrees fails the gap test and is derived.

**Edges through a pole on reading.** Both readers cut an edge through a pole
of a face it bounds (a sphere's pole, a cone's apex, a patch's collapsed side)
before any face is built. IGES trimmed surfaces (143, 144) cut every boundary
segment through a pole of any surface read with them, so neighbours sew piece
to piece. Not cut: a trim placed by one subfigure instance against a pole of
a surface placed by another, and pcurves other producers build through a
pole. The cut edge stays in the model unused, and a STEP style or shape
aspect naming its file id lands on every piece.

### Mesh conversion

`solid_from_mesh` rebuilds planes, the four canonical surfaces, extrusions,
surfaces of revolution and fitted B-spline patches, and leaves every other
region faceted; `MeshSolidReport::fallbacks` names each curved region it
faceted and why. Each change is held to the converter's standard (verified
against every sample, edges placed on both surfaces, facets as the fallback)
and measured on `mesh_corpus` and the truth bench before and after.

**Seams between fitted regions.** Each region is fitted on its own and each
seam solved point by point along the mesh boundary, or threaded as a chord
through the boundary's vertices. Folds and faces that mesh open are caught
after the build and their curved faces withdrawn to facets, so every
measured part tessellates closed at the cost of those faces. A chord where no
point carries onto both surfaces (a torus and a cylinder fitted to pieces of
one spline blend) can stand 0.2 off, stated as its tolerance. Design items,
in order: absorb single facets at junctions without costing an exact seam
(taking every facet within reach breaks exact rims); trace seams between
solved corners as one branch of the intersection guided by the mesh path;
junction rules (Benière et al. 2012: a junction stays only if its faces are
pairwise adjacent, dangling edges dropped).

**Fillets fitted apart from their supports.** Derived from their supports:
rounds between two planes, a torus between any two supports on its axis
(cylinders, cones and planes square to it, or one of them and a sphere
centred on it), corner balls, and a sphere put on the torus a round turns on
beside a plane. Not derived: a torus between two spheres, a piece of a torus
between curved supports met as a sphere, variable radius, and fillets between
curved faces off a common axis, which keep their fitted surfaces and meet
their supports at a near tangency. A torus corner a few facets round between
two fillet cylinders and a floor fits nothing (every sample reaches into the
cylinders) and stays facets.

**Sweeps.** A whole turn is recognized, its axis settled on the samples. Not
done: the band layout for a sweep between two circles (its seam the profile
itself), and a closed profile (an extrusion of a closed curve, a revolution
of one). The helical line complex would tell a thread from a surface of
revolution; its threshold has to sit above the twist a mesh's own normals
carry (a pitch of about 1e-3 on a turned wave meshed at 0.01). A thread goes
to the patch.

**Patches.** A smooth region nothing else fits becomes a fitted B-spline
patch where it is one disk with one loop and verifies
(`MeshSolidOptions::patches`). Limits:

- a tangent line inside a row of triangles (a scan, or a mesh without
  vertices along its edges): the curvature's jump falls inside one span, the
  fit runs out at 7e-5 to 4e-4 against 1.5e-5, and a near-flat foot gathers
  into narrow slanted planes;
- a run-out steeper than the fit can follow stays faceted (a hill on a round
  bar at 0.8 high, the fit at 1.24e-5 against 1.23e-5), and a mesh too coarse
  to leave a hill enough vertices keeps the sphere round the hill's facets;
- run-outs into cones and tori are threaded as into cylinders and spheres,
  but not measured;
- between vertices a fit is no closer to the true surface than its knots
  allow: on coarse meshes up to 4 times the distance in the row at a free
  edge, under 1 inside.

**A torus pierced with no circle free either way.** Holes in a whole torus
are inner wires, its seams placed on a free parallel and meridian or threaded
through the holes. Holes that leave no parallel and no meridian free, or that
overlap along the chain at every level tried, stay facets.

**Region edits.** `MeshRegions` `fit` has no sweep kind, and the tangent
passes do not revisit a neighbour an edit changes (a ball beside a merged
round keeps its fit). The corner and boundary graph of the seam work above is
what an application would edit next.

**Inflections on a coarse open mesh.** Where the turns change sign, the
coplanar distance estimate asks them to change linearly, so a coarse mesh of
a surface whose curvature varies fast can read its inflections as scatter.

**Converted parts of the real assemblies.** Fifty solids of the DoomCube 350
assembly and fifty-three of the Voron 2.4, each meshed at a thousandth of
its diagonal and converted back, come back valid within 3e-3 of their
volume but for four. Voron solid 990 comes back with a planar face of
0.09 square millimetres that `check` calls turned in; the cylinder beside
it meshes 0.023 off its surface at a chord of 1e-4 (the mesher says it
missed), so a winding number over that mesh is no witness either way. Not
diagnosed. DoomCube 630 and Voron 1170 and 1200 come back valid, 3.1e-3
and 3.2e-3 under their volume.

**A face whose rims touch within a micron meshes open.** One solid of the
DoomCube 350 assembly has a cylinder face between a circle and a spline edge
that dips to within half a micron of it and touches it: its ring crosses
itself in the chart at any chord, the crossing segments are refused as
constraints, and the face loses boundary points its neighbours keep (27
edges open). Clipping such a ring ear by ear closes the mesh but folds about
half the face's area (4 of 7 square millimetres) and hides from the mesh
conversion the faces it withdraws because they mesh open; it was not kept.
A face this degenerate wants healing (the spline edge projected onto the
cylinder, or the sliver removed) before it is meshed.

### Speed

**A per-face state cache.** The boolean's build phase could read a piece's
classification instead of probing for it. No case needs it for correctness.

**The rest of the local boolean.** Faces whose box misses the other solid's
are set aside (`ogeom-bool/src/aside.rs`): not gathered, split, classified or
rebuilt, and passed through as the same nodes; a plane's holes clear of the
other solid are left out of its arrangement and put back into the piece
holding them. What a boolean reads of a solid of 64 faces or more (the
setting aside's walk of every face, the classifier's prepared faces and
the rings drawn of them) is kept in the model (`Model::kept_read`) and
read back by the next boolean on the solid unchanged; any edit below it
forgets it. A corner hole into a plate of 582 faces costs about 1.5 times
the same hole into one of 22 (heavy test `local_boolean`, bound 2.0).
Still linear in the solid: `side`, the history copies and the shells of
the faces passed through, about 0.4 ms at one thread on the 582-face
plate, and the first boolean on a solid reads all of it. A solid made by a
boolean is read whole by the next one; deriving its read from the
operand's kept read and the faces made is the lever for a chain of
booleans. Where the assembly around the faces set aside does not close (the
sew welds a shared vertex, a seam join reaches a shared edge, a strand would
widen a shared vertex), the boolean runs again over every face; converted
slabs with tolerant vertices take that path.

**Operations that edit in place.** A boolean notes the nodes it passes
through as held by its operands too (`Model::note_held`). Heal, sew, the
fillets and chamfers, `remove_faces`, `unify_same_domain` and `tessellate`
copy what is held below the shape they are given before editing
(`Model::unshare`, `ogeom_algo::on_own_nodes`). Still edited where they
stand: a node held whole and handed in itself (a face set aside, given to
`fix_face_pcurves`); a profile's edges, which a sweep (`make_prism`,
`make_revolution`, the pipes and lofts) and `split_face` describe on the
surfaces they build; a boolean's own vertex widening and re-charting of a
ball. Each grows a tolerance, or adds a pcurve on a surface none of the
operand's faces lie on (clearing the edge's same-parameter flag). Copying
there instead parts the profile from the solid it came from: a decision
for those operations, not taken here. Sharing made by
anything but the boolean (`Reshape`, compounds, placed instances) is not
noted, and the notes are not saved with a document.

**Fillets in fewer passes.** A chain is taken in rounds, no two edges of a
round sharing a vertex, each round in one boolean each way (wedges whose
boxes meet in as few booleans as keep each one's boxes apart). Everything
else in a call reads only what is near each edge, so a call costs its
booleans plus a few milliseconds per edge: the four outer edges of a plate
with 256 holes take about 2.5 times the plain plate's time, every rim of
it about 25 ms a hole. What is left is the boolean's own cost over the
whole solid, and a polygon's outline still takes two rounds. Building the
later round against the solid as given, as though the earlier cuts were
made, and cutting every wedge at once (fused first, or as disjoint pieces)
gives a tool the boolean cannot close against even a plain box ("the kept
pieces did not close into a shell"): its legs lie on the box's faces all
round the outline.

### Quality and performance

From a read-only audit of every crate (2026-10-07). Nothing below is
measured yet: each item starts with a bench, and lands only where paired
runs beat the threshold it declares (the performance-optimizer skill).
Kernel threading already runs through `ogeom_core::parallel::map_ordered`
(scoped threads, results in item order, bit-identical at any thread count).

**Benches.** `tools/ogeom-bench` (see its README) reports min, median
and MAD per bench after a warm-up, loops short benches to 50 ms of
samples, takes `--threads N` and `--filter`, and `threads()` honours
`OGEOM_THREADS`. The baseline is recorded at one thread so its ratios
compare across core counts. Benches now cover mass properties, `check`
and `fix_shape` on the largest corpus part, sew of a 1536-face shell,
`solid_from_mesh` on a tessellated part, STEP write, IGES read (of the
kernel's own IGES of a corpus part: the corpus has no IGES file), weld and
STL read at 1M triangles, `project_exact` and the polygonal HLR, a
thickened spline face, a pipe through circular sections, a guided pipe, a
marched fillet and a box with all 12 edges filleted. Still missing: a STEP
file of 10 MB or more (the largest corpus file is 1.2 MB). On the
largest part `solid_from_mesh` does not finish a warm-up and three samples
in two minutes, so its bench runs on the small part. `project_exact` takes
about 30 s at one thread there (7 s at 20), nearly all of it Newton on its
spline faces, so its heavier bench runs on a 520-face part.

An instruction-count gate in CI would be deterministic on a loaded runner
and looks feasible: `iai-callgrind` needs valgrind on the runner (an apt
install on ubuntu-latest) and a bench target with a few of the benches
above at reduced size, since callgrind runs them 20 to 50 times slower;
`perf stat -e instructions` needs `perf_event_paranoid` at 1 or below,
which hosted runners do not promise. The work: an `iai` bench crate, a
recorded count per bench, and a gate at a few percent.

**Cheap and certain.**
- `map_ordered` starts fresh threads for every stage, two items included.
  Measured, and left alone: the items are heavy (faces, edges, face pairs)
  and few (a drill hands over 2, 2 and 5), so the default count already
  beats `--threads 1` on the small benches (boolean_drill 1.0 against
  1.45 ms, fillet_block 20 against 37 ms). A minimum batch per worker
  serialises those stages: 2 items per worker cost fillet_box_all 3%,
  4 cost it 10%, 16 nearly doubled fillet_block. A scoped spawn costs
  8 to 16 us per thread here, about 8% of boolean_drill and under 1% of
  the fillets. Letting the caller take items as one of the workers gained
  3% on boolean_drill and 7% on boolean_local, under the 10% bar; a
  persistent pool could win at most that 8% on the smallest booleans.
- Closed-form feet on curves: `project_on_curve` in
  `ogeom-algo/src/measure.rs` samples a line or a circle like any curve
  (`line_parameter` and `circle_parameter` exist in
  `ogeom-math/src/elementary.rs`), and refines with Brent's minimizer
  where Newton on `(C - P) . C'` after the bracket would do. Surfaces
  (cone, torus, extrusions and revolutions through their profile) already
  have closed-form feet. Reading the DoomCube and Voron 2.4 assemblies
  projects only onto splines this way (33,000 calls, 0.2 s of 27 s), so
  the gain is for other callers.
- The intersection walker: six or more unknowns still go to the heap
  solver, `Condition::system_at` allocates on every evaluation, and each
  step corrects from the last point rather than a tangent prediction.
- One jet call where separate point, first and second derivative calls
  are made: `CurveOnSurface::derivatives_at`, the mesh lift (point then
  normal), and the side probes of mass integration; allocation-free curve
  derivatives (a curve jet; `derivatives_at` returns a `Vec`).
- Box rejection: `general_3d` in the curve crossings samples both curves
  before any whole-curve box test, and `section_face.rs` does not filter
  edges by box; section-against-edge paving tests the other face's box,
  not the section's own; a spline's sub-range is bounded by its whole
  control polygon.
- The guided sweep law samples each guide at 257 points at every station
  and bisects each sign change 60 times.
- The mesh still draws an edge again for each face's rings and chords and
  in each settling round; its sealing passes build the edge-use map three
  times.
- Primitive refinement in the converter clones its parameters twice per
  point and builds a new matrix every iteration; piece nesting ray-casts
  every pair of pieces against every triangle with no box test.

**Algorithms worth the work.**
- Classify a boolean's pieces once per connected region, not per piece
  (Requicha and Voelcker 1985): pieces meeting across an unpaved sub-edge
  share a state. Measured: classification is about 5% of
  `boolean_many_faces` (36 probes over four drills, faces the tool misses
  already set aside) and 0.1% of the fillets, so it waits for a case where
  pieces are many; it needs the split and classify of `general_fuse_as`
  apart first. The boolean's time is the arrangement (about 2.9 ms of the
  four drills) and assembly (about 2.1 ms).
- Surface feet on a trimmed or offset spline still seed from a grid; the
  span search on a plain patch could clip its spans to a trim window, and
  bound an offset by its basis's boxes grown by the offset distance.
- Curve-surface intersection prepares its sample grid once per face and
  clips the ray to the face's box; point-in-solid rays pay it on every
  call now.
- March sections at the cubic's own interpolation error and check the fit
  between stations (Bajaj et al. 1988; Barnhill and Kersey 1990), instead
  of a straight-chord step about ten times finer than the fit needs. The
  fillet march takes stations eight times denser than its fit for the
  same reason. Measured on two crossed cylinders: the joint fit is half
  the fuse (5 to 7 ms a branch from one span through 1556 points); a fit
  started from knots placed by the cubic's error estimate is a third
  faster with a tighter true miss, but moves sections by sub-micron
  amounts and the shells of a few stress drills and two tests close or
  refuse the other way (patch kept as `hermite-stations-fit.patch` in the
  scratch crate). Shell closure that flips on sub-micron section changes
  is the item to settle first.
- March general pairs over the faces' parameter boxes, and cache a pair of
  surfaces marched once within one fill.
- Marched seeds search the cell pairs that come close without crossing
  for a crossing or for the stationary approach a small loop runs round
  (Sederberg and Meyers 1988). Where the surfaces lie close along a whole
  curve (a blend against its support) the squares around a search that
  found nothing are passed over, to keep the cost within the grid's, so a
  small loop beside such a band can still be missed. Normal cones over the
  spline patches' hulls would prove those squares empty instead.
- The winding number fallback of the exact classifier (where every ray
  meets a tangency) sums the faces' meshes triangle by triangle; a
  hierarchical far field (Barill et al. 2018) once a solid of many faces
  asks it often. Exact solid angles of trimmed planes would spare meshing
  them. A probe within a curved face's mesh band at the finest chord still
  refuses.
- Arrangement darts ordered by exact orientation, ties walked apart.
- Primitive fits with analytic gradients and Levenberg-Marquardt (Lukacs,
  Marshall and Martin 1998).
- Fitting a band or a pcurve: factor the normal matrix once per round for
  all rows, a cyclic system as band plus border; warm-start knots from
  the last round.
- Exact offsets of extrusions and revolutions through their profile's
  offset; spline offsets through the local-refinement fit.
- Mesh `simplify` by quadric error metrics with a heap and the link
  condition (Garland and Heckbert 1997); the triangulation keeps triangles
  by the exact flood rather than a centroid parity.
- Fillet corners cut together where their tools stand apart, not one at a
  time against the whole solid; the ruled fillet through the round's
  batching; the curved corner probes the exact boundary instead of meshing
  the whole solid.

**Structure.** `fill` (2,300 lines) and `general_fuse_as` (1,500) in the
boolean, `mesh_solid.rs` (12,900 lines), and the largest builders in the
fillet and offset crates split into named phases with typed outputs, which
the region classification and the parallel read-only phases need. Mesh
conversion moved to its own crate would stop it rebuilding the whole
graph.

**Consumers' builds.** Workspace profiles do not reach a dependent crate,
and without a target CPU every `mul_add` (515 of them) is a call into libm.
Measure `--release` against the default profile and `x86-64-v3` against
the baseline target, and document what a consumer should set.

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
- **No stubs.** A thing is implemented, or it is refused by name and tracked
  as above.
- **Scope is parity with the reference kernel's modelling modules, plus what
  `docs/SCOPE.md` admits.** `docs/SCOPE.md` is normative and says how to
  decide a case. Parity is about capability, not structure.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against,
  bundles or imports another kernel, and the design is worked out here, not
  mirrored. The field's vocabulary is used throughout.

## Decisions, not gaps

These are settled, and listed so nobody reopens them by accident.

- **A surface returns its point and derivatives in one evaluation, agreeing
  with the separate accessors only to rounding.** A foot-point solve needs
  all six values at one place, and a patch answering three accessors locates
  its spans, builds its basis and sums its grid three times.
  [`Surface::jet_at`] answers from one order-two table (3 to 25% of a
  spline-rich STEP read). A patch sums its point by de Boor and its
  derivatives by basis functions, which can differ in the last ulp.
  Consistency within a jet is what a Newton step needs and what the type
  guarantees; a caller must not mix a jet with the accessors at the same
  parameters and expect identical bits.
- **A pcurve with no closed form is `None`, not a fit.** An exact curve with a
  fitted pcurve would have two descriptions that disagree by an amount nothing
  records. The consumer that needs one marches the pair instead.
- **A closed exact section partly outside a surface's extent is kept whole.**
  The restriction that matters is the face's trim, the boolean's own 2D
  stage; the surface extent is only a parameterization window.
- **Scaled placements in the boolean are baked first.** A scale changes a
  surface's parameterization under its pcurves; `baked_shape` rebuilds the
  geometry, and the boolean calls it.
- **The crossing walker refuses tangential contact.** The tangential walker
  owns that case, and the section pipeline routes to it. The refusal stays
  pinned.
- **Bi-tangent construction is subsumed**: by the 2D repertoire in 2D, and by
  the blend family's own envelope in 3D.
- **Glue is subsumed** by the boolean's same-domain unification, which skips
  nothing it needs and unifies what glue would.
- **A reversed face walks its wires backward, and its outer wire is stored
  first.** `ordered_children_of` reverses a reversed parent's children, which
  a wire's walk needs and which lists a reversed face's holes before its
  outer wire. Code that wants the outer wire reads `Model::outer_wire`; code
  that wants every wire under the face's sense, outer first, reads
  `children_of`. Most readers pick the outer ring by area.
- **A fillet refuses an edge between tangent faces.** There is no corner to
  round ("the edge's faces are tangent"). The stress fillet scenario's
  refusals are all of this kind or an edge splitting two coplanar faces, and
  its crease list takes both.
- **A marched blend's run-out is capped in the ball's section.** Where a round
  runs off a face, the cap stands in the ball's section through the crease's
  end, which leans with the faces; a cap in a side's plane would end the band
  off a section.
- **Each face's box is kept in the model.** The model keeps a face's tight
  box in the face's own frame, found the first time it is asked for
  (`ogeom_algo::face_bounds`, `tight_bounds`, the boolean's pair filter) and
  read in constant time after (`Model::face_bounds`). It is forgotten when a
  node of the face is handed out by `node_mut` or a surface by
  `surface_mut`. `tight_bounds` of a 582-face plate a second time costs a
  twentieth of the first. Faces a boolean rebuilds are new nodes, so a
  result's boxes are found again (see "Sharing unchanged faces' nodes").
- **The converter uses no learned model.** Published learned reconstruction
  (2024 to 2026) gives a valid solid for 70 to 76% of parts under a hundred
  faces, every face a B-spline, and places its own bottleneck in rebuilding
  the b-rep, which is the converter's verified build. A learned proposal is
  also a result nothing here measures. Global selection by an integer program
  is out for the same reason: its time has no bound.
- **SAT, X\_T and JT are not read or written.** Their specifications are
  unpublished, and implementing them would mean reverse engineering files
  instead of reading a standard. If a specification is published, this lifts.
