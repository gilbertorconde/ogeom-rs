# Plan

What is left to build, and the rules the work follows.

## Where the remaining work is tracked

- **GitHub issues** hold bugs and requests, one per issue.
- **The capability ledger** (`docs/CAPABILITIES.md`, generated from
  `docs/capabilities.toml`) records the kernel's named capabilities. The restriction on each `partial` row is the
  list of what that capability does not do yet. `tools/check.sh` holds the
  ledger to its evidence.
- **This file** holds the open items that are neither: known limits of the
  code, refusals by name, design items and speed work.

Finished work is recorded in `CHANGELOG.md` and in the tests that pin it, not
here.

## Open

Within each group the items are ordered by risk: a possible wrong answer
first, then refusals, then cost.

### Booleans

**Sections through a sphere's pole are refused.** A marched section whose
image on a sphere face runs past the chart's pole row by more than its
tolerance, or passes within it of one of the face's poles between its ends,
is refused ("a section through a sphere chart's pole is not resolved"): the
image cannot follow the half turn of longitude at the pole, and the
arrangement read it on the far side of the chart. Before the boolean,
`ball_chart.rs` restates a whole ball or a trimmed sphere face, an operand
or a member of a compound one, on a chart whose poles stand clear of the
other operand's surfaces and of its planes other than one through a pole,
so the section meets no pole; a band whose two loops go round the new axis
takes a seam between them, its loop edges cut where the seam crosses them.
The face is left on its own chart, and the section refused, where the
rechart cannot run: a face with one new pole inside whose loop round it the
seam meridian crosses at no vertex, a face placed inside its solid, a face
with an edge that has no curve in space, and a face whose trim holds no
pole of its chart (a rounded box's corner whose chart's pole stands beside
the other operand's plane outside the corner). Recharting such a face on
the face path as it stands splits rounded corners differently and once
left a face crossing another.
A section passing a hair beside the pole of a face left on its chart is
refused the same way ("a branch passing beside a chart's pole fits nothing
yet"): its image fitted between the samples holds, but the arrangement does
not yet split and classify the piece it bounds round the pole rightly.
What would clear them: cut the loop at the seam as the band's is cut, or
march the section through the pole, breaking its image at the pole sample
so each piece has an image of its own.

**Nearly coplanar faces cross only where both are lines.** Two planar faces
leaning a hair apart cross at a small angle, and a strip of one between the
crossing and its edge, narrower than the confusion distance's reach divided
by the angle, stands on the other at every point. On the attempt that
settles such pieces from both sides (`settled` in `ogeom-bool/src/lib.rs`),
the edges bounding the strip are carried onto the other face from where a
plane-plane section line crosses them (`graze_contacts`), so the other face's
matching strip is split off and the two are kept as one contact: a cube
and a box whose wall leans 1e-5 or 1e-4 radians off the cube's, crossing a
micron inside its top or bottom edge, fuse to the closed form within the
weld distance over the band (`boolean_contact.rs`). The crossing is found
only between a section line and a line edge: a strip bounded by an arc,
or by a marched section, is not split this way, and the piece is refused
as before ("a contact over a region of faces not recognised as coincident
is refused") where nothing else splits it. Where the strip is wider
the ordinary crossing splits both faces exactly, and pairing the band as
same-domain would move the result by up to the weld distance over it, so
the band is not paired there.

**A section whose fit misses widens its junctions as far.** A marched
section keeps its fit's measured error as its tolerance, and the paving's
reach and the junctions along an edge grow with it, unbounded. Over the
workspace suite (heavy tier included) and the stress suite at seeds 1 to
4, a boolean's sections no longer state more than 2.1e-3, with a sphere
chart's pole met only where the chart is left unturned, which refuses.
Three still miss by more than a few chords. A section of two skinned
lens runs ends with its sample 5.8e-4 past the edge of one patch's chart,
and its pcurve, lifted there, stands 2.1e-3 off the curve
(`sections_off_the_spine_round_a_lens_glue_their_runs_on_the_mitre_planes`
in `crates/ogeom/tests/pipe_laws.rs`). A skin running tangent to a mitre
plane (`an_off_centre_section_down_an_l_mitres_between_skinned_runs` in
`crates/ogeom/tests/pipe_corners.rs`) states 1.6e-3, the surfaces' normals
there parallel to 1e-14. A drum meeting a
face blend's patch at a few ten-thousandths of a radian states up to
5.2e-4 (the stress face-blend scenarios), where the section solved
between two samples lands tens of chords off their chord and is not
filled in.

**A drill's seam beside the other part's face.** A drum touching a plane
of the other part to within the confusion distance is the touching drill
(see `drills_a_hair_off_a_plate_top_running_into_a_round_cut_as_the_touching_one`
in `crates/ogeom/tests/boolean_grazing_drill.rs`). What is left:
- The drill lying in `nist_ctc_03`'s plate underside between its rounds,
  its seam toward -z on the underside, is refused at every offset of up to
  a tenth of a micron, the touching one included ("the kept pieces did not
  close into a shell", raised by a tenth "marching resolved no branch").
  With its seam toward z, -y or y it cuts at all five offsets.
- Lowered under `nist_ctc_03`'s plate top past the confusion distance, by
  2e-7 to 1e-6, the stress drill seed 2, drill 5 placement is refused for
  every seam direction ("the kept pieces did not close into a shell"), as
  before: the drum then stands apart from the plate, and its section with
  the round crosses the drill's top line two to five microns past the
  round's tangent edge.
- A drum moved onto a plane carries its seams with it but not its rims,
  which its neighbours hold. Moved by the confusion distance itself (raised
  or lowered by 1e-7) with the seam on a diagonal, the results are valid
  but their volumes miss by about the move times the wall's length and
  radius: the common of the 2 by 30 drill through the block, lowered, by
  2.8e-6, and on the plate's drill of length 970 the drill less the part
  and the common add up to the drill short by 1.5e-3, against the
  touching drill's 2e-7.

**A loop's pcurve stops short of its seam column after a fuse.** A
section that crosses a torus's seam a hair from its own start is marched
as a long piece and one a few hundred-thousandths long across the seam.
The fuse drops the short piece, the loop closes in space within its
vertex's radius, and on the fused face its pcurve stops 7.6e-5 short of
the column it closes on (`a_fused_torus_face_whose_loop_stops_a_hair_short_of_its_seam_is_cut`
in `crates/ogeom/tests/boolean_chain.rs`). A later boolean reads the two
ends as meeting within their radii and cuts it, but the fused face's ring
is open in the chart. Ending the long piece's image on the column breaks
its agreement with its own curve by the dropped piece's length, and the
exact volume integral then gives the face up to the mesh. The fuse should
carry the short piece into the long one, across the closed section's
join, which a piece's range does not yet pass.

### Fillets and blends

**Curved face blends refused by name.** `blend_faces` and `fillet_faces`
round curved faces that share no edge where their surfaces share a direction
or an axis, and otherwise march the ball round where the surfaces cross;
between separate faces a march closed round its seat is a band fitted closed
round the loop. They refuse (`curved_face_blends_refuse_by_name` in
`tests/blends.rs`, `face_rounds_refuse_by_name` and `curved_rounds_refuse_by_name`
in `tests/sheet_fillets.rs`):

- faces the ball leaves at different places round its seat and never both at
  once. The round would be two runs, one over each overlap of the faces'
  stretches. Each face's end is found by cutting it along its line of
  contact over a run opened where the ball is off every face it leaves;
  with no such station, any opening runs one face's line out of it and back,
  and the cut gives two edges. Owed: each face's stretch read with the loop
  opened in that face's own gap, the stretches intersected round the loop,
  and one capped wedge per run applied in turn;
- a line of contact crossing a face more than once: several runs between the
  same two faces, and with `trim` a face cut along several lines whose piece
  away from every round may be more than one;
- corners where the surfaces cross in more than one place the ball rolls
  in (two drums fused, a drum's flat and its wall), rounded corner by
  corner by `blend_faces` and `fillet_faces`
  (`a_face_blend_fills_both_creases_over_older_smaller_rounds`,
  `two_drums_round_both_creases_their_walls_meet_in`), are told apart by
  the sum of their reaches, each the larger of its crease's distances to
  the contacts, or else by a line parting the convex polygons that hold
  each corner's section: its crease, its contacts and the meeting points
  of each arc's end tangents
  (`a_face_blend_takes_off_two_corners_whose_balls_overlap_but_whose_corners_stand_apart`).
  Owed: corners whose polygons overlap while the corners do not, told
  apart by their own boundaries; and with `trim` over several corners,
  each later cut keeps the earlier rails, which a cut rebuilding them is
  refused by name;
- a round crossing a face that stands between the two (a ball crossing a
  chamfer between the faces it rounds), where the solid beside the round is
  open where the seat says material, or the reverse. No ball rolls in that
  corner touching only the two faces: on the test's drum, chamfered 1, the
  ball of radius 1.5 touching the wall and the lid has its centre 1.41 from
  the chamfer's plane. The round would leave the material through the
  chamfer, giving round, chamfer, round meeting at creases. It needs a ball
  touching the chamfer as well (a three-face blend). On every path, two
  planes, swept and marched, `blend_faces` refuses by name where another
  face of the solid, on any surface, holds a point inside the ball at a
  station of the run (`between.rs`;
  `a_face_blend_refuses_a_ball_that_crosses_a_bevel_between_the_faces`,
  `a_face_blend_refuses_a_spline_bevel_or_a_bevel_over_part_of_the_edge`,
  `a_face_blend_refuses_a_wedge_its_ball_crosses_over_part_of_the_foot`,
  `a_face_blend_refuses_a_wedge_on_the_foot_s_edges_its_ball_crosses`):
  a chamfer the ball crosses, uneven or not, plane or B-spline, a round
  smaller than the one it would replace, a wedge fused into a drum's foot
  over a short stretch, the foot parted from the top or its edges shared
  either side of the wedge. The run is read whole, halved between stations
  until the centre's travel fits under the face's clearance. Between two
  planes and on the swept path a face inside the ball counts only inside
  its lune, the part between its centre and its arc
  (`a_face_blend_rounds_a_part_thinner_than_its_ball_where_the_lune_stays_inside`),
  and a face on a plane along a straight run is cleared of every lune at
  once from the run's two end sections
  (`a_face_blend_rounds_a_bar_whose_far_faces_pass_just_behind_the_ball_s_centre`).
  The marched path reads the ball's side inside its lune where no face
  stands in the ball
  (`a_marched_face_blend_fills_over_an_older_round_just_smaller`), but
  still asks the whole ball, so a part thinner than its ball is refused
  there. Owed: the marched seat's lunes, their planes turning along the
  run. An older marched round just smaller than the ball (within a few
  hundredths of it) is fitted to 2e-4 and reaches into the ball by a few
  micrometres beside a contact, where its border runs within its own
  tolerance of the face; the blend is refused there, as the fitted round
  does stand inside the ball. Faces the
  crease or a line of contact meets, and faces holding a capped run's
  ends, are where the round runs off (a crease running off a block's
  side) and are not asked. Owed: a face that ends the run without
  meeting either line or holding an end (a narrow face past the one a
  line of contact leaves by) is asked as a face between; and the face's
  distance is the kernel's point-to-face distance, which takes the foot
  on the face's surface only where it lands inside the trim, so a face
  whose nearest interior point is not its surface's is read from its
  edges.
- with `trim`, a closed line of contact across the seam of a closed B-spline
  face whose chart does not wrap round (a tube swept from a B-spline ring):
  `split_face` brings a trace back by a period only on a periodic chart, so
  the cut does not part the face into pieces covering it once. Owed: the
  line cut at the seam on such a chart, or the face restated on a periodic
  one.

**An edge round whose ball is wider than a face it rounds.** `fillet_edge`
refuses where the ball's contact on one of the edge's faces stands past
that face's far side at every station: the face is narrower than the
setback all along the edge (a sliver face where two oblique boxes are
fused, a triangle face whose apex stands under the setback). Where the
face is a planar strip whose far edge turns the same way as the edge
(convex both: a chamfer, a narrow top), the ball rolls between the faces
either side of the strip and, where the strip lies wholly in the corner
it takes off, the strip goes with the corner (`consumed.rs`). Owed there:
the concave mirror (a flat in an inside corner), whose fill's caps need
the solid past the strip's ends read as open, and ends not square to the
edge. Where the far edge turns the other way (a concave seam where two
solids are fused at a slant), the ball touches the other face and rests
on the far edge, its centre on the ellipse where the plane a radius off
that face meets the cylinder a radius about the edge, and where it rests
on that edge alone over the whole crease the round is the band through
its arcs (`rolled.rs`). Where the ball changes what it holds along the
crease, the round is a chain of bands meeting tangent where it holds
three things at once (`carried.rs`): a face on each side (the cylinder
between the planes), a face and the other face's edge (the ellipse), or a
face and the other face's corner (the circle over the face about the
corner). Of the 11 stress blend refusals of that kind at seeds 1 to 4,
five round alone: seed 1 cases 12 and 48 and seed 4 case 31 on the edge
all along, seed 3 case 33 from one strip's far edge down onto both faces
and onto the other's far edge, and seed 2 case 54 from the far edge round
the triangular strip's corner. Seed 8 case 46 rounds in its request, edge
to corner to edge. The rest are refused by name. In seed 3 cases 17, 32
and 54 the ball comes onto the face beyond the far edge, where it touches
neither of the edge's faces on that side and rounds the corner between
the face it touches and the face beyond: that is the round of the seam
between those two, which the request does not name, and needs that seam
in the request with a corner closing the vertex where the seams meet
(seed 3 case 32 also starts where its strip's end edge turns the ball's
sections back). In seed 1 case 36 and seed 3 case 4 the ball rests on
the far edges of both faces at once, its centre on the quartic where the
cylinders a radius about the two edges meet, which has no closed-form
piece here. In seed 4 case 48 the ball's contact on the face it touches
runs over a concave seam of that face onto the face beyond. Seed 1 case
12 rounds alone, but another edge of its request shares a face with it. Within a request, the ball also refuses where its contact
runs off a face anywhere along an edge that shares a face or a vertex
with another of the request (4 refusals): the band's contacts are checked
against its neighbours' and an off-face contact holds none to check. Owed:
the overlap check read on the stretch of the band the face holds, and the
contact beyond it carried onto the face the band is cut off by.

**A pole leg whose edge winds round the pole and whose rail does not.** A
marched band whose rail passes over a sphere's pole has a leg from the
winding loop to the pole with the other loop cut out as a hole. When the
winding loop is the rail, `pole_leg` in `marched.rs` splits it clear of the
hole and starts the seam there. When it is the edge, the seam leaves from the
edge's own vertex and can cross the hole. No measured part builds one.

### Construction and exchange

**The converter widens its distance for one fit.** The default coplanar
distance doubles, twice at most, when any one recognized fit presses
against it, as a fit that runs on tangentially into its neighbours does:
`nist_ftc_10` widens once and `nist_ctc_05` twice, each for one or two
fits. Regions do not depend on how the mesh lists its triangles, so
neither does the widening, but a wider distance lets facets of surfaces
nothing recognizes merge into planes: `nist_ftc_07` widened twice comes
back on 319 faces with its volume 6e-5 off, against 583 and 3e-7 at the
default. Exporter scatter is everywhere and a run-on is local, so the
distance should widen only where many fits press, or for the pressing
region alone.

**STEP pcurves the file cannot state exactly.** The reader keeps a file's
pcurve where its lifted gap is within twice the projection's, restated in
this kernel's units (the file's length and plane angle units, applied to
each surface parameter as ISO 10303-42 parameterizes the surface) and,
where the file states it at another pace or phase than its curve's, read
at the affine map of the curve's parameter that puts it on the edge's
ends: a conic either way round and turned, a line at any pace, a spline
over its own domain either way. Still derived: a trig or offset trace and
a pcurve of a curve written as a spline conversion, none of which STEP
states at a linear map of the curve's parameter; an ellipse turned other
than by half turns; a conic on a chart stretched unevenly; and pcurves on
a line swept into a plane or on a surface replica, whose charts here are
not the file's stretched along each direction. The writer still leaves out
a reversed conic pcurve: written forward, it would ask other readers to
take a pcurve at another pace than its curve.

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
faceted and why. It reads meshes drawn in strips along rulings, as
other kernels commonly draw drums and cones, as it reads ogeom's own. Each
change is held to the converter's standard (verified against every sample,
edges placed on both surfaces, facets as the fallback) and measured on
`mesh_corpus` and the truth bench before and after.

**Seams between fitted regions.** Each region is fitted on its own and each
seam solved point by point along the mesh boundary, or threaded as a chord
through the boundary's vertices. Folds and faces that mesh open are caught
after the build and their curved faces withdrawn to facets, so every
measured part tessellates closed at the cost of those faces. A chord where no
point carries onto both surfaces (a torus and a cylinder fitted to pieces of
one spline blend) can stand 0.2 off, stated as its tolerance. A single facet
whose seam with a curved face would be a chord wider than the reach is built
as a fan onto that face instead (a wedge split at its corner first), its
seam lifted onto the face: on `mesh_corpus` every edge between a curved face
and a single facet then holds within the reach except on `nist_ctc_02`,
whose widest (0.35) border two-triangle planes along a spline blend's tori.
Design items, in order: those two-facet planes, which no fan takes; seams
between two curved faces fitted to one spline (a torus beside a cylinder at
0.09, nearly tangent, no point on both); junction rules (Benière et al.
2012: a junction stays only if its faces are pairwise adjacent, dangling
edges dropped).

Tracing a seam as a branch of the two fitted surfaces' intersection has
nothing to trace on the measured parts. Of the 427 chords on `mesh_corpus`,
258 stand wider than the reach (twenty coplanar distances), and 254 of those
because one surface stands farther than the reach off the other at a point
of the mesh path. A signed distance moves no faster than the point it is
taken at, so no point on both surfaces lies within the reach of that point,
and Newton's step from the path lands 0.5 to 15 off it. The widest are all
of this kind: `nist_ctc_02`'s planes and tori fitted to pieces of B-spline
faces (0.63 over a span of 28, the torus's equator bowing off the plane by
its sag), its cylinders beside planes (0.46), and its tori beside cylinders
(0.11, meeting at a quarter of a degree to seven). The other four are cubics
through three vertices standing a little wider than the polyline would. Of
the 169 within the reach, 3 have a crossing within it all along. The fix is
upstream, in fits that share their neighbours' tangency or take the spline
face whole; the chord stays as the fallback. The planes beside those tori
and cylinders were single facets, now built as fans (above): their seams lie
on the curved faces, and of the main corpus with the step-conformance parts
the widest seam on a curved face is 0.35 (p99 0.046), against 0.63 (0.16).

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
  edge, under 1 inside;
- a ring (a spline blend round a pocket) is no disk, and cutting it into
  disks does not help on a coarse mesh. On `nist_ctc_02` its blends are
  rings four facets across (radius 10, chord 1.17). Halved, then halved
  again, the pieces fail: the mean-value map flattens triangles with all
  three corners on one side of the square, a half ring folds over its mean
  plane, and the quarter rings' patches stand 3 to 5 percent past a
  triangle's sag at its middle, or hold too few vertices inside. One piece
  of eight verified, and the facets left beside it widened the worst seam
  from 0.35 to 0.73. A ring needs a chart that goes round it (periodic in
  one direction) and enough rows across to verify against.

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
volume. Without recognition (`MeshSolidOptions::recognize` off) the build
asks nothing of the solid it returns, so a mesh passing through itself comes
back a solid crossing itself there; with recognition the conversion is
refused, naming where.

**A face whose ring runs back along itself is meshed, not healed.** Where a
face's boundary runs out along a path and back along it within the widest
tolerance of its edges (a spline edge rising to a circle and running on
along it to the vertex they share, a hair past it in the chart), the ring
crosses itself at any chord. On its last pass the mesher takes the tail out
of the ring and draws it as a strip at its edges' points
(`ogeom-mesh/src/tails.rs`), and the rest triangulates whole. The model
keeps the tail: the face still bounds a stretch of no width between two
edges its neighbours do not share with each other, and `fix_shape` strips
such a pair only where no other face walks either edge. A tail whose way
back leaves the way out across a segment wider than that tolerance is left
crossing, as is a ring that crosses itself anywhere else.

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
plate, and the first boolean on a solid reads all of it. A solid read for
the first time takes each face it shares with a solid whose read is kept
from that read (`Model::kept_reads`): a boolean's result takes the faces
it passed through from its operand's setting aside and prepared boundary,
rings drawn included, and reads only the faces made; the nodes below it
are walked down to the nodes a kept read holds. The setting aside asks
whether a face of the tool lies outside the solid of the boundary the
general fuse classifies against, prepared once. Along a chain of 200
holes and bosses into a plate (heavy test `local_boolean`) a step at 306
faces costs 1.31 ms against 1.48, the chain 8 percent fewer cycles and 12
percent at 361 steps, to the same result; at 66 faces nothing changes.
Taking the faces' reads is under 2 percent of that at 200 steps: the
edges' holders are rebuilt in the solid's own order from each face's walk
as a fresh read builds them. Still linear along the chain: those holders,
`side`, the history copies, the shells' closed check and each new face's
rim walk (`explore_unique` per face; a plain walk saves 2 percent).
`edge_pcurve_gap` answers a line, a circle or an ellipse against a line,
a circle or an ellipse on a plane, or an iso line of a cylinder, a cone,
a sphere or a torus, in closed form where that stands within the edge's
tolerance (`pcurve_gap/closed.rs`), within rounding of what its sampling
measures; the chain then costs 22 percent fewer cycles at 200 steps and
17 at 361, and 80 percent of the gaps the stress suite measures are
answered so. Still sampled: B-spline pcurves, the chart traces of
oblique sections, helices and any gap past the tolerance. A STEP read
measures the pcurves it derives with the same forms where it acts on the
gap only past the confusion distance
(`pcurve_fit::lifted_gap_past_confusion`): the Voron 2.4 assembly reads
in 23 percent fewer cycles at one thread. Still sampled there: the
derived pcurve's gap where the file states pcurves of its own, which
bounds the ones kept.
Measured and not landed: the classifier draws a
cylinder wall's rings at the full chord though its pcurves are straight
in the chart, about 15 percent of `boolean_drill` and of a step of the
chain and 40 percent of `offset_part`; drawn at the fewest points the
walk's period choices (rings placed by the mean of their points) put a
NIST part's wall a period apart in places, so the fewer points are not
the same rings. Boxing the section against each edge before crossing
them rejects almost no pair on the benches. Where the assembly around the faces set aside does not close (the
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
- **Scope is what `docs/SCOPE.md` admits.** It is normative and says how
  to decide a case. Scope is about capability, not structure.
- **Independence.** See `CONTRIBUTING.md`. Nothing here links against,
  bundles or imports another kernel, no one working here reads another
  kernel's implementation source, and the design is worked out here, not
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
