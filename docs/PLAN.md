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
the side's plane would end the band off a section. Upright, with the axis
on the edge or 0.5 past it, the closed-form blend is an exact torus ended
where the solid's edges along the crease end, capped in the meridian plane
there, and fills the ring's share of the foot's arc within 1e-10 relative,
as the edge round of the sharp foot does. Refused by name:

- Faces of a solid that the ball leaves at different places round its
  seat and never both at once: neither face gives the round both its ends.
- A marched round closing on itself between separate faces, and one whose
  line of contact crosses a face more than once.
- A round crossing a face that stands between the two (a ball smaller than a
  chamfer between the faces it rounds): the solid beside the middle of the
  round is open where the seat says material, or the reverse.

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

**A tilted half ball on a drum.** A half ball charted about a tilted axis
fuses, cuts and commons with a drum standing under its flat face, sharing
the face and its rim, and with a drum through it, also where the chart's
seam lies in the flat face (the axis level, the seam meridian level too)
and the common is empty. One limit remains: drilled through
on the chart about the cube diagonal (1, 1, 1) with its seam toward z, the
sums fuse + common and cut + common stand 1.2e-8 relative off the
operands, where the other charts tried hold 5e-9: the rim and the section
are fitted on the tilted chart, held to 7e-6 and 2e-5, and the measured
volumes move by that much.

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

**A drill lying in a face that runs into a round.** A drill whose lowest
line lies in a plate's underside meets the round at each end of the plate
within a few microns of tangency. Its section with the round is a figure
eight crossing the round's tangent edge twice at the double point, or two
loops crossing it a hundredth of a millimetre apart, which is one junction
where the edge cannot be told from the drill between the crossings. Both
close on the corpus drill along x with the drill's seam toward z. With its
seam toward -y the marcher returns the section with one round as six open
pieces whose fits miss by up to 9e-3, paving the round's edge at several
places with honesties up to 2e-2, and the kept pieces do not close. The
same placement moved by a few ulps leaves the drill less the part
unclosed.
The stress drills still refused: one on the rounded box at seed 1, six at
seed 2 (a converted plate, the frustum, `nist_ctc_01`, `nist_ftc_11`, the
rounded box, the shaved cube) and the rounded box at seed 3, some not yet
diagnosed. Five drills through `nist_ctc_03` once refused there too; they
cut since the STEP reader states how far its pcurves stand off their edges.

**Stress fillets still refused.** At seeds 1 to 3 the fillet scenario
refuses 23, 21 and 21 of 160, all by design: every edge of the rounded box
joins two tangent faces, and the stepped part's fuse leaves its base's and
its step's coplanar side faces split along an edge with no corner. The
harness's crease list takes both kinds.

**Edges walked the same way by both their faces.** The data model has each
edge between two faces walked once each way, senses composed: each face
keeps its material on the left of its rings about its surface's normal.
Prisms, tapered prisms, revolutions, baked shapes, revolution bands,
sewing, the STEP reader, the offset crate's lofts, pipes and sweeps,
the mesh converter, rebuilt faces whose chart turns (a reflecting bake,
a restated surface), seaming and ring re-anchoring, surface recognition,
face removal, the IGES writer, sheet rounds, faces rebuilt by draft,
offset or move, and every blend wedge do. `check` names an edge with two
faces in a shell, closed or open, that both walk the same way, as a
suspect finding.

A face whose ring is wound against its own surface's normal (a planar cap
built from a ring walked clockwise about the plane it was given) is
ambiguous: sewing trusts the ring and turns it with the faces beside it,
so a ruled sheet closed by such a cap can come back with that cap facing
in. Where the cap has no pcurves the shell cannot be meshed to settle it.
`make_face` given a plane and a ring could refuse or turn the plane to
agree with the ring.

- Sewing rebuilt a reversed face by reading its wires under its own sense
  and then reversing the rebuilt face again; it reads the stored wires
  now. That made a different rough box for the converter tests (the
  tessellator numbers a reversed face's vertices along its rings), and on
  it the converter bounds a fillet cylinder by an edge 0.018 loose with
  the straight iso line as its pcurve. STEP now carries each edge's
  pcurves where the file can state them exactly (`SURFACE_CURVE`,
  `SEAM_CURVE`, `PCURVE`), and the reader keeps a file's pcurve where its
  lifted gap is within twice the projection's (the iso line stands 0.0183
  off the curve, the projection 0.0096); the converted parts come back
  with every pcurve and their volume to the last digit. Left: an edge's
  pcurve the file cannot state in the curve's parameter (a reversed or
  left-handed conic, a trig or offset trace, a curve written as a spline
  conversion) is still derived on reading, and a foreign file's pcurves
  are read without its length or angle unit, so one in inches or degrees
  fails the gap test and is derived.
- The STEP reader stores a face turned against its surface with its
  loops walked back, so every third-party corpus file reads with each
  edge walked once each way (it was 166 edges on `nist_ctc_01`). The two
  ogeom-written screw files were written from faces read the old way and
  have their reversed faces' bounds flipped; any other file written so
  is repaired on reading: where a shell's edges pick out a consistent set
  of faces turned against their surfaces whose loops, walked back, leave
  every edge walked once each way, those faces are rebuilt and the report
  says so. No third-party corpus file sets it off.
- Reversed faces now tessellate in another order, which showed how much
  the mesh converter hangs on the order a mesh lists its triangles and
  their corners: with the old reader's meshes, shuffling the triangles or
  turning each triangle's corners gave `nist_ftc_07` up to 40 faceted
  regions and `nist_ftc_10` up to four. Three fixes hold: a first sample
  too small to be fitted as a torus grows to the next stage instead of
  ending the seed (`nist_ftc_07` comes back whatever corner each triangle
  lists first, and corpus parts come back on fewer faces: `nist_ctc_02`
  1807 to 1609, `nist_ctc_04` 581 to 483, `nist_ftc_10` 420 to 331); a
  curved fit whose distances round past the tolerance is refused (a
  sphere centred 1e17 off took 740 triangles of planes, cylinders and
  tori on `nist_ftc_10`); and facets beside a face that reaches past its
  triangles are fanned before the face is faceted (the sliver part, any
  triangle order). Two more held once the reader fix changed the
  meshes: a last-resort seam is threaded only through vertices lying on
  both faces (a cone moved 0.077 off a corner of `nist_ftc_08` was joined
  by a chord 0.12 loose instead of falling back), and a curved face that
  took facets in gives them back before a face beside a turned-in facet
  is blamed (on `nist_ftc_10` a coarse cylinder that absorbed a strip of
  facets covered its neighbours, and the corner ball beside them was
  faceted; the part now comes back on 314 faces with nothing faceted).
- Open: with the reader fix a torus corner of `nist_ftc_08` a few facets
  round (32 triangles between two fillet cylinders and the floor) is no
  longer recognized. Every sample of it reaches into the cylinders and
  fits nothing; a one-row sample of its bottom fits a cone, which takes
  eight of its triangles and three of the floor's, and the rest stay
  planar facets. Nothing falls back, but the corner is facets where it
  was a torus.
- The offset crate's 42 are fixed, every measured volume and centroid
  unchanged: a ruled loft's planar faces walked their corners in section
  order whatever side the material was on (the bottom cap always against
  the walls); the torus segment's half-tube patches ran clockwise in the
  chart, and its far cap walked the tube circle forward; the planar caps
  of skinned lofts and pipes, helical sweeps and pipe shells walked the
  border rings in the walls' order whatever way the outward normal
  pointed, and now walk each ring so the material is on its left (an
  outer ring turning positively about the normal, a hole negatively); and
  a pipe shell's flat strip read its normal off its first row's chord and
  first column's, which on a strip starting from a curved corner's join
  row lie along one line, so the strip's ring could run clockwise about
  its plane. It reads the turn of the whole sampled border now. The
  integration tests count same-way edges with `tests/support/walks.rs`.
  A shell check over the offset crate's tests and the sweep, loft, pipe
  and thread integration tests now flags only draft (`apply_draft`'s
  rebuilt faces) and the fillet wedges a thick solid's rounded corners go
  through.
- The converter stores every face's rings counter-clockwise in its
  chart, walked back from the triangles where the face is turned (it was
  134 edges on `nist_ctc_01`, 246 on `nist_ftc_10`, none left on the
  corpus or the truth bench). Faces and volumes are unchanged: an exact
  volume taken about a fixed point agrees to the last digits; about the
  first face's anchor, which the walk moves, a converted part that does
  not close exactly differs by up to 6e-6 (`nist_ftc_10`).
- A shell check run over the default tests (each closed shell an edge of
  which two faces walk the same way) flagged 80 tests and now flags the
  three that build a face turned inside out on purpose
  (`check_orientation`, `fix_shape`). Run on each blend wedge before it
  is applied as well, it flagged 105 more tests and now 26, all wedges
  whose results are sound: marched blends (converted edges and
  rims, spline edges, a branch cylinder's seam, the crossed bores' loop,
  pinched drums), face blends capped flush or running off a side, crease
  arcs split at a seam, and a tangent chain's wedges; each walked one to
  six edges the same way. None does now, over the whole workspace with
  the heavy tier: a rim arc's end caps walk back where their corners turn
  against their outward normals and its other faces run their rings
  anticlockwise in their charts; a marched band's run-out legs and caps
  walk back by the side the rail and the corners stand on, its annular
  legs run the outer ring anticlockwise and the hole clockwise, and a band
  between two rings walks its lower ring forward.
- `check` gained the rule. Over the workspace with the heavy tier and the
  stress run at seeds 1 to 3 it names only the faces the orientation and
  repair tests turn inside out on purpose, which expect it; no sheet is
  named, so open shells are asked too. The count costs under a
  millisecond on the corpus parts (0.3 ms on `nist_ctc_02`'s 663 faces,
  whose check takes 0.34 s). A marched band widens its rails and end arcs
  with their vertices, so its wedges and the blended solids pass `check`
  with nothing broken.
- Fixed with the same cause (rings read under a face's sense, then the
  face turned again): seaming and ring re-anchoring, surface
  recognition, face removal and the IGES writer (whose reader takes loops
  about the surface's normal and turns the face by the shell's flag).
  A rebuilt face whose new chart runs the other way (a reflecting bake of
  a plane, a recognised spline facing its axis) walks its rings back, its
  pole rows with them; a band of rings running clockwise stores its loop
  walked back; a sheet round's ring is made counter-clockwise before the
  round is turned to face its axis. Draft, face offset and move read the
  stored rings too. A blend wedge's planar face winds its corners about
  its outward normal, a fillet cap walks back where its corners turn
  against it, and a rim's annulus walks its wider ring with its normal.

- Sewing orients by the walk: in each sewn group the faces are turned
  breadth first over their shared edges until every edge two of them
  share is walked once each way, each connected run keeping the
  orientation most of its faces have (its first face's on a tie), open
  sheets too. A closed group whose faces were turned is then turned whole
  where its meshed volume comes out negative. A group with no consistent
  orientation (a half-twisted band) keeps its faces as given, and
  `Sewn::edges_walked_one_way` names the edges that show it. The probe
  and its two guards are gone. Over the workspace with the heavy tier,
  14 groups turn (17 faces), each face genuinely against its neighbours:
  boxes with a lid filled the wrong way round (two sew tests and two of
  the measured box's four builds; walls facing in turn the box whole), a
  sheet with one square wound the other way, a floor filled facing up
  into an extruded square, a blend surface whose ring runs against both
  strips it joins (six blends in the bridge tests), a floor and wall
  built to walk their crease the same way (the thick sheet test now turns
  the sewn wall back itself), and a ruled sheet's bottom cap whose ring
  is wound against its plane's normal (it now faces into the box; its
  caps carry no pcurves, so the volume cannot be meshed to settle it).
  No other result changed. The stress run at seeds 1 to 3 turns nothing
  and its outcomes are unchanged.

The fillet's side test no longer trusts the walk where the chart's chords
cannot settle it.

**Two convex fillets meeting where a curved face takes over.** On the
shaved cube a bottom edge and the upright edge at its end, rounded in that
order, take off the union of the two blends, the bottom's flush cap
leaving a sliver standing against the drum. Rounded the other way, the
bottom blend runs on through the upright's band and trims that sliver
too, 0.8% more at radius 0.18 and 4% at 1. Both are valid; one corner
should give one solid.

**Stress results on exact volumes.** The stress harness holds its volume
identities to 1e-6 on exact volumes, and at seeds 1 to 3 every drill is
measured exactly and holds (worst 1.8e-7, a converted corpus part;
generated parts within 7e-9). Nine drills once fell to a tessellation.
In eight (the torus drilled along z or y with its wall through a point of
the outer equator's seams, and drills through `nist_ftc_11` tangent to a
plane holding a circular edge) the boolean bounds the drill's wall, near
where the section hugs an existing circle, by that circle: the wall's
pcurve, fitted to the true section, stands up to 4e-5 (1.2e-4 on the
corpus part) off the edge's circle while the edge states 1e-7. The exact
integral now takes the strip between the lifted pcurve and the edge's
curve as a ruled surface with the face, up to a thousandth of a
millimetre. The ninth, on the rounded box, is a junction where two pieces'
ends stand 6e-5 apart inside a vertex of tolerance 4.2e-5: a junction may
now open to the vertex's diameter. The torus commons stand within 1.5e-9
of a ring-by-ring integral of the torus's thickness (2.6e-8 without the
strip), and the nine identities hold within 5.4e-9.
The sewing that takes the section's edge and the circle's for one edge
now raises the kept edge's tolerance to how far the carried pcurves stand
off its curve, sampled at 33 points (2.2e-5 to 4.0e-5 on the torus drills,
up to 1.2e-4 on `nist_ftc_11`, where the edges stated 1e-7), and the exact
integral still takes the strip on a curve with a closed form whose pcurve
stands off it by more than a hundred confusions, whatever its tolerance.
On such a curve the strip is now taken as wide as the edge states, not
only to a thousandth: mesh conversion bounded a planar face by a straight
chord stating its sag (0.49 on `nist_ctc_02`'s blends) while the torus
beside it kept the arc, and the sliver between was a hole in the
boundary. Five tori there came out 1.8% under their mesh area and the
converted volume 1.29e-3 under the part, against 1.7e-4 for its fine
tessellation; with the strips, 1.75e-4. (The converter no longer leaves
those chords: the converted part now measures 4.3e-5 under, 3.7e-5
without the strips.) A face the chart rectangle or disc would take goes
round its chart loops where a pcurve stands that wide of its edge, so it
gets the strip too. An area counts a strip as surface the face lacks, or
takes it away where it lies back over the face, decided once per lobe (a
lobe ends where the pcurve crosses the curve or a ruling's foot comes to
rest at the curve's end, and breaks the panels); it was counted with the
sign of the walk, so half the strips subtracted. Moving the reference
point by 1340 now moves the converted volumes by 1904 (`nist_ctc_02`,
4241 before), 15 (`nist_ctc_04`, 4888) and 134 (`nist_ftc_10`, 308): a
closed boundary would not move at all. Corpus areas move by up to the
strips' own area, 4.5e-5 of a face on the original `nist_ctc_02` and 0.4%
on small faces of `nist_ftc_07`.
Left: a fitted curve's pcurves within its stated tolerance are still
taken as along, so a converted face bounded by fitted edges stating
2e-2 to 5e-2 still leaves slits (two converted `nist_ctc_02` tori 0.2% and
0.7% under a fine mesh); whether they are what still moves with the
reference is not measured. The original `nist_ctc_02` moves by 923 with the reference, with
or without the strips.

`check` compares an edge's pcurves with its curve only where the edge
claims `same_parameter`, which almost no producer sets, so a pcurve
leaving its curve by more than the edge states reads valid. Comparing
every edge (the lifted pcurve against the nearest point of the curve's
stretch, 33 samples) flagged the readers, `fix_shape`, mesh conversion,
the geometry rebuild and the boolean. Each now raises its edges'
tolerance to where its pcurves stand, and with the stricter comparison
every workspace test passes, ignored ones included, as does the stress run
at seeds 1 to 3, so the check can follow them.
The boolean passes it. A fuzzy boolean's confusion is its fuzz, and the
pieces of an operand's edge, a contact and a section kept their source's
tolerance only above it, so edges of a faceted round stating 1e-6 to
4e-6 came out at 1e-7 with their pcurves where they were; the floor is
now the smallest tolerance. The sew measured carried pcurves beyond the
edge's tolerance or the confusion, whichever was wider, which under a
fuzz let twins up to 9e-6 apart through; it measures beyond the edge's
own. A contact's piece states how far its image on the target's surface
stands (a rim on a drum a micron off its own axis: 1e-6 at 1e-7). The
drill lying in a plate's underside cut a line of `nist_ctc_03` lying
4.6e-6 off its plane, which the reader now states and the pieces keep.
The STEP and IGES readers pass it: every pcurve they attach, exact or
fitted, is lifted against its edge's curve (`pcurve_fit::lifted_gap`: 257
samples, the nearest point within a step where the pace differs, and a
golden-section search for the summit about the widest peaks and both
ends, where a fit through a near-pole sweeps across the chart in a few
thousandths of the edge) and the edge widened to it. A line of
`nist_ctc_03` lies 4.6e-6 off its plane, a sphere's fitted pcurve strays
up to 2.6e-4 near a pole in `nist_ftc_07`, and converted solids and
sphere cuts read back 3e-6 to 7e-6 off, each at 1e-7 before. Past a
millimetre the reader warns instead of widening, since an edge that wide
swallows its neighbours. A curve through a pole of a face it bounds (a
sphere's pole, a cone's apex, a patch's collapsed side) has no single
chart image, so both readers cut such an edge at the pole before any face
is built, on the edge every face shares: a meridian circle of a sphere in
`nist_ftc_06` stood 3.9 off its fitted pcurve and now reads as two exact
meridians at 7e-9, the hemisphere's area exact, 2 fewer warnings, and the
drills through the part cut as before. Not cut: IGES trimmed surfaces
(144), whose faces sew whole edge to whole edge, and pcurves other
producers build through a pole. The edge that was cut stays in the model
with no face using it and nothing referring to it: a STEP style or shape
aspect naming its file id lands on every piece. IGES carries no attribute
on an edge of a 186 solid, whose edges are entries of an edge list.
With its B-spline edges' gaps (up to 3.9e-2) stated, the exact integral
takes `nist_ctc_02`, which moves its volume from the 0.1 mesh's 4.710940e7
to the exact 4.710056e7, against 4.710192e7 on a 0.005 mesh.
`fix_shape`'s tolerance reduction measured an edge at 9 samples, pcurves
only (no seams), and shrank fitted trims straying 0.019 to 0.029 between
them to microns; it now measures with the same helper.
Mesh conversion passes it: a plane face keeps a curve's closed-form
image only where it lands within the edge's tolerance (an arc fitted on a
neighbouring face may cross the plane, and its closed-form image is then a
circle of its radius, thousands off), and otherwise takes the curve's
projection; a last
pass raises each edge's tolerance to where its pcurves stand, sampled
denser than the check.
The geometry rebuild (`to_nurbs`, `to_nurbs_within`, `restate_geometry`)
passes it. Its pcurves on an offset surface came from a foot-point solve
that stopped at its seed, since an offset has no second derivative; the
solve now steps by Gauss-Newton there, and an offset drum's pcurves, 0.11
off at 1e-5, and an offset box's, 1.18e-7 at 1e-7 (the divide's case),
lie on their edges. An iso line on a fitted surface is kept only where it
follows the edge within the fit target, and every edge given an iso line
or a fit states how far its pcurves stand (256 samples, the nearest point
of the stretch).
A surface fit (`to_nurbs_within`, a patch's degree restriction) is made
at the surface's own parameters and measured against it at every span's
middle as well as at the samples; a miss splits the spans it lies in.
Against a 200x200 projection both ways, every fit of an offset drum,
plane, cone, torus and spline, and of a plain cylinder, sphere, torus,
cone and spline, now lies within 1e-4 and 1e-5 as asked (the drum at 1e-4
stood 1.29e-3 off between the old chord-length grid's samples, a sphere
0.14). A fit takes 4 to 95 ms where it took about 1 ms, and an offset
spline at 1e-5 takes 0.57 s. The offset drum's rims state 8.4e-5 at 1e-4.
Left: the stated error is the worst at the samples and span quarter
points, which can sit a little under the true miss; a fit whose surface
cannot be evaluated somewhere (an offset sphere's poles) refuses as
before.
Curve fits measured only at their samples follow the same rule through
`ogeom_algo::traced`: a fit starts at every span of its source (each cut
in degree + 1, at least the old count), is checked at each interval's
quarters, splits the intervals that miss, and reports the worst measured
(budget 8192 samples). Where a fit and its pcurve were made apart, the
edge states `lifted_gap` or `state_pcurve_gaps_of`. Measured before and
after, each pinned by a test that failed before:
- `reanchor_boundaries`: a 64-span edge weaving 0.1 moved onto a plane
  1.0e-1 off its shadow while reporting 3e-12; now 3.5e-5, reported 3.6e-5.
  `fix_face_pcurves` widens the vertices with the edge (they stayed at 1e-7
  under a 1e-2 edge, which `check` flagged) and states the lifted gap.
- `recognize_surface` verified on a 9x9 interior grid: a 0.059 bump in a
  plane patch's outer twentieth read as a plane with a certificate of 0. It
  now also verifies on a grid over the whole chart, two per span each way
  (16 to 256 a side), and the rebuilt face and its edges state the
  certificate and their pcurve gaps. `simplified_edge` verified at 17
  points: a 0.067 bump in an 80-span boundary became a line stating 1e-7;
  it samples every span now and a replaced edge states its deviation.
- IGES 142 with only a parameter-space trim: a 301-corner polyline trim
  lost corners by 5.9e-2 while warning 1.0e-5. A spline trim's chart is now
  exact (the maps are affine) and the lift is fitted piece by piece between
  its corners and joined C0: corners within 1e-6. IGES 130 types 2 and 3:
  on a 400-span base the offset strayed 4.3e-3 between samples while
  warning 3.5e-4; the law is evaluated at any parameter now, 7.1e-6.
- `make_face_with_pcurves`' projected fallback: a trim 2.5e-4 off a
  300-span edge stating 1e-7; the edge now states the lifted gap.
  `make_band_between`'s fitted connector on a sphere: 4.9e-5 at 1e-7, now
  fitted to 1e-6 between samples (1.5e-7) and stated. `project_edge_onto_plane`'s
  fallback (a curve on a surface, an offset): 2.0e-2 off a 400-span curve
  while stating under 1e-5; it is fitted at the edge's own parameter now
  (9.6e-6), so the result runs same-parameter with the edge.
- Boolean section pieces: images up to 1% past what the edge stated
  (2.38e-6 at 2.35e-6 in the pole drills); the piece states its gap as a
  contact piece does, and `pcurve_onto_ends` refits at every span.
Measured with no miss in the result, fixed anyway: a fillet leg's fitted
apex image stood up to 1.3e-3 off an apex edge stating 1e-7 (the boolean
that consumes the leg restates it); `attach_image` understated by 14%; a
half space's sides on a 60-span patch stood 3.5e-2 off at 1e-7 (they never
reach a result). No miss: the marched pair's contact line (gaps up to
9.98e-7 under its 1e-6), and the marched bands, whose points stand within
9.8e-5 of exact sections solved between the stations (target 2e-4); the
ball centred off the band's normal misses the hosts by up to 1.1e-3, a
normal-direction error no tolerance states. HLR's marched silhouettes fit
through on-curve points and sit 1.9e-6 off a torus's outline at a 1e-5
chord, closer than the walk; their comment no longer claims an error is
added. Left: `projected_into_shared_chart` ignores the owner fit's `met`
(the contact piece states the dense gap after), heal `split.rs`'s joint
fit and `divide.rs`'s iso fit (one case measured, 8e-15) are measured at
samples only.

The offset builders' fits are now measured the same way. Each builder
hands `fit_surface_sampled` (or `fit_curve_sampled`, or
`fit_trace_sampled` for a curve with its pcurve) the geometry it samples,
fitted at its own parameters and checked at every span's quarter points;
where the least-squares fit strays the spline through every sample
decides where to refine. Measured densely against the true geometry,
before and after: `make_pipe_skinned` along eight periods of a sine at
1e-3, 9.6e-2 to 3.1e-4; `make_filling` round a half disc's arc from four
samples at 1e-4, 2.1e-2 to 1.6e-5; `make_loft_skinned` through ellipses
10 by 0.2 at 1e-4, 1.4e-2 to 3.0e-5; `make_pipe_shell` of a half disc
down a bent spline at 1e-5, 6.5e-4 to 4.6e-6 on its arc;
`make_pipe_sections` between two such ellipses at 1e-3, 0.18 to 7.2e-4;
`normal_projection` of a circle onto a ball from 8 stations at 1e-3,
2.5e-2 (stating 2.2e-2, its edge 1e-7) to 9.3e-4 (stating 9.6e-4, the
edge widened to it). Covered: the skinned pipe, the skinned lofts and
their aligned, closed and apex forms and the loft a pipe law builds, the
strips of a faceted pipe shell, a cornered loft and a helical sweep, a
smooth profile's wall in a pipe shell, both closed pipe rings, a planar
strip's borders, the multisection pipe's rings, the filling, the
projection and the extruded draft's wall; an adopted strip border states
how far its image stands. A wavy wall drafted 0.1 at a 1e-4 target stood
2.5e-3 off its turned rulings and stands 3.1e-5 off; the wall's tangent
continuation past the profile's ends now reaches as far as the rulings
lean, for a neighbour standing across it at a slant. Its rebuild takes
0.27 s (about 1.5 s with the old fit).
The offset, sweep and stress runs cost about the same; `thread_groove`
takes about 8% longer.
A pipe shell's curved runs and closed rings are checked between their
stations too: the frame is carried from the station behind by one
double-reflection step and turned to meet the next station's normal, and
the skin's parameter across is even in the station index. A disc of
radius 0.3 down a spline turning a quarter within one span, at 1e-4,
stood 1.1e-2 off its tube and stands 1.5e-5 off (both laws); round an
ellipse 5 by 2, 1.3e-4 to 3.0e-5. Each wall takes about 20 to 100 ms
more. The law loft (auxiliary and binormal laws) is held to the first
section moved by the law's frame between its stations, the law asked for
its normal there: the same disc under a guide beside the spine, at 1e-4,
6.4e-3 to 3.6e-5; under a binormal at 1e-5, 5.4e-5 to 2.8e-6. A guide
starting on the first station's plane no longer has its crossing taken a
step along it (the start section turned by up to that step). The helical
strips are checked between their stations too (the screw is exact at any
turn): a disc of radius 1 a thousand out from the axis, once round, stood
4.5e-5 off at its 1e-5 target and stands 5.5e-6 off. A swept skin's
parameter across runs from naught to one over its stations. Cost: a
helical sweep builds about half as long again (a spring of 3 turns 0.49 s
to 0.75 s), a thread groove's pipe 0.13 s to 0.26 s; the thread groove
and helical sweep tests take 0 to 20% longer, the stress run the same.
`general_draft` fitted its hinge to 1e-4 and its rulings' tips (a unit
out) to 1e-4 over the reach, through 256 stations: the hinge's own error
grew along the rulings, and a drum drafted 0.1 about a plane tilted 0.7
stood 2.0e-3 off its exact rulings 20 out. The wall's two border rows are
now fitted on one knot vector (the surface fit through them) and checked
between the stations, where the crossing is found again from the chart
between the stations' feet, held on the chord between them so its pace
does not kink: 1.8e-5, and the draft about 30% faster there.
The patch closing a non-planar end of a skinned loft is the cone from the
end section to its centroid, now fitted to that cone between the samples
as well: rings waving 0.8 twelve times round, the end rising and falling
2, at 1e-4, stood up to 9.6e-2 off it (and failed `check`) and stands
6.0e-5 off; a plain saddle end at 1e-3, 1.8e-3 to 4.3e-5. The cap is
turned away from the next section's middle: the test from the section it
passes through is square to a nearly flat cap's normal, and could leave
it facing in.
A pipe shell's curved run ending on a curved corner's crossing (each
column along its own generator to the crossing, read by its own length)
is checked between its rows too, and so are a planar strip's end columns
wherever the skin is known along them. A disc of radius 2 down a quarter
arc of radius 20 turning square onto a straight leg, at 1e-4, stood
7.3e-3 off its tube near the corner and stands 7.7e-5 off; a square
section through the same corner was refused (an end column through its
rows only, 5.9e-3) and stands 6.6e-5 off. A column's length is read on
the stations' grid, so the arc's end, where an off-spine generator
changes pace onto the straight extension, is a step's end: read across
it, the column's parameter kinked and the fit refined to a 421 by 259
grid (11 s) at 1e-4 and failed at 1e-5. Now 0.1 s at 1e-4 and 0.3 s at
1e-5 (9.0e-6).

The exact area of a trimmed spline face grades an inner panel across
which `|n|` dips almost to nothing (the surface folding over) towards the
dip's bottom, found on the Legendre series of `n` at the panel's nodes;
only a face whose first run shows such a panel watches for it. A wavy
wall drafted 0.1 at a 1e-4 target (a 9 by 173 patch whose normal turns
over in a sliver near a crest) doubled its panels five times, about 5 s,
and stood 4.6e-10 off its area; it settles on the first doubling in
0.07 s, within 1e-10. Volumes and every other area the tests measure are
unchanged to the last bit.
The two fits held between their samples are one: `fit_traced` and
`fit_traced_2d` are the open case of `fit_curve_sampled` and its planar
twin `fit_curve_2d_sampled`, one refinement in either dimension. The
least-squares fit first and the spline through every sample where it
misses, every span measured at its eighths, an open fit refined to 8192
spans (`SAMPLED_SPANS`) and a closed one to 1024 (`SAMPLED_LOOP_SPANS`:
its join couples the ends, so its system is solved whole, and a closed
ring at 1e-7 that took 8192 spans cost 1 s a fit). Every caller was
measured against its own trace at 20000 points over the range and 16 in
each final span, over the default tier, the heavy helical sweep and
thread groove tests and the stress run: no fit stands further off than
before. The sweeps' fits are unchanged; the traced ones came closer and
mostly cheaper. A projected helix stood at 0.97 of its target and stands
at 0.20, on a third of the samples (2 s to 0.03 s); a re-anchored
boundary 0.42 to 1e-10; a varying IGES offset stood 0.1% past its target
between the quarter points and stands within it.
Left:
- `check` passed a skinned loft whose nearly flat end patch faced into the
  solid (the exact volume then came out 607 for 1278), and the mesh's
  volume did not see it either.
- A blended pipe section fitted to a tenth of the sweep's tolerance
  stands 0.8% past that tenth between its eighths on one ring (still a
  tenth of the tolerance); the eighths do not see it.
- `general_draft` on a closed face whose seam is closed only to position
  (a skinned loft's wall): the face's normal turns across the seam and the
  exact draft's rulings with it, so the drafted wall closes on the mean
  ruling and blends to each side's within a thirty-second of the hinge.
  There it stands off the exact draft by half the turn times the reach
  (6.5e-4 on a loft of circles of radius 10 skinned at 1e-3, reach 9);
  elsewhere within 1e-4. Such a draft rebuilds in about 1.3 s (0.23 s
  before), the time in the rebuild's exact volume over the wall.

Left: the rest of those commons' miss (up to 1.5e-9) is the fitted sections' own pcurves
standing off their edges within the 1e-5 taken as along, which every exact
volume carries: taking a strip on every edge that stands off at all
brings the torus commons within 1e-13 and their identities within 1e-14,
but makes the stress run 5.6 times slower and leaves seven drills on a
tessellation, so it is not done.

- **A drill all but touching inside a torus's outer equator.** A drill
  touching the torus inside its outer equator meets it in a figure eight;
  the marcher now finds the touch exactly and cuts the branches there, and
  these cut (at the seams' vertex and anywhere on the equator, every seam
  direction tried, each sum within 1e-8). Moved inside by 1e-7 or 1e-6
  with its own seam on the touching line, the section is two loops, above
  and below the equator, each turning sharply a thousandth or so from it
  across the drill's seam. Drills of radius 1.5, 3 and 4.35 at 26 places
  round the equator, both offsets, measure exactly: cut, common and drill
  less the torus within 4e-8 of the torus against a ray-by-ray integral of
  the common. The misses of 3e-4 to 1e-3 once read here were the
  measure's, not the boolean's: these results fell to a 2e-3 chord
  tessellation, and with chart loops left open at their junctions the
  exact integral was off by 1e-5 to 1e-4. The last two (radius 1.5, 0.8
  and 6.3 radians round, 1e-7 in) once refused (the kept pieces did not
  close): the drill's seam line was never found crossing a loop at its
  sharp turn, so the drill's wall kept no band between the loops. A loop
  fitted whole wobbled by 3.5e-4 there, more than the marcher's finest
  step, and is now fitted in pieces; another fit all but stopped along its
  parameter beside the crossing, holding the curve-curve Newton short of
  it, which now walks on by alternating feet. A classification guard
  refuses a cut, common or union
  whose sections leave every piece of both outside the other while the
  solids overlap, which is what the touching drill returned before (the
  torus whole).
- **A drill's seam along the other part's edge.** A drill whose wall runs
  along a round's tangent line with its own seam on that line: the section
  there hugs an edge on both sides and lays no split. Its seam turned away
  it closes. A drill tangent to a plane along a face's edge with its seam
  there is the same; carrying the seam onto that face as a contact closes
  it, but the same contact breaks a drill whose seam crosses a face's
  interior lying in its plane, so it is not carried.
- **A section tangent to an edge at its end.** A drill whose circle on a
  rounded box's base touches the round's tangent line at the corner vertex,
  where the corner ball touches the base.

The sections grazing an edge to the fourth order now close (a drill lying on
a frustum's base through its rim, or along a round's tangent line into the
corners). Drilled so at radii 0.3 to 2, cut and common add up to the part
within 1.2e-10 relative, and the drill less the part and the common add up
to the drill within 8e-9.

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
converter's own distance (79 and 67 before). The estimate also takes an
edge whose neighbours turn several times as much for a lull in a bend, and
end rows whose turns grow toward the edge for steady, comparing edges
continuing one past its ends only within five degrees of it: at its own
distance (1.65e-5 on every set, what single precision resolves) each row
stays a face of its own and the faces a pad's floor crosses hold their
vertices, doubt at most 1.6e-5 where rows merged 6e-4 off their planes
gave 4e-2. Of the 15 committed row-count sets at drafts 0.001 to 0.0001
and depths 3 and 10, 147 of 150 pass at that distance (137 before). A rung
between two rows leaning out of a wall by a couple of thousandths of a
radian is met by the sections either side a micron apart along it, and
the edge's doubt there runs past a micron where the crossing's gap is
within the weld distance, up to two microns, the most a joined crossing's
vertex then claims; and the mesher reads an inner ring's width
as twice its area over its perimeter, its mean width, so a ring of facets
poking a sixth of a micron through a wall, a micron and a third wide,
stays a hole. Of 17 row-count sets at drafts 0.001 to 0.0001 and depths 3
and 10, 340 of 340 pass at both distances (334 before). Drafted 0.00005
and 0.00002, a top facet's shared edge is taken for the crossing wherever
the facet's far side stands off the wall by more than the weld distance,
though the band within it runs a third of the way down the facet. Where a
section's end vertex is widened where it was welded, the sew compares each
end of two edges within that end's own tolerances, so the two sides of a
sliver a quarter of a micron wide stay two edges (the pad's floor at a
corner row): the one result found outside its doubt (2.1e-3 cubic
millimetres off, doubt 1.0e-3) now matches. Of 17 row-count sets at depths
3 and 10, each row a face of its own / at the converter's distance: draft
0.00005 34 / 28 of 34 (20 / 14 before), 0.00002 32 / 29 (18 / 15),
0.00001 24 / 18 (16 / 16); every other case refuses. What remains:

- Same-domain pairing for planar faces within the weld distance of each
  other over a region. Drafted 0.00001 the facet stands within the weld
  of its wall over most of its height; where a corner facet meets the wall
  only at a vertex, their planes' solve crosses it diagonally a tenth of a
  millimetre off and the kept pieces do not close (10 of 34 refuse with
  each row a face of its own, 16 at the converter's distance). Pairing the
  band as coincident, split where the planes part by the weld, is the
  design; dropping the solved line when each face lies on one side of the
  other's plane changed nothing measured. The other refusals at 0.00005
  and 0.00002 (depth 10 with merged rows, rising and falling counts at
  depth 3) are a lower corner row's facets crossing the pad's walls, or a
  ray meeting a tangency, not the top band.
- The mesher tells a slit from a hole by width alone: slits measured
  about six tenths of a micron wide, the ring above a micron and a third,
  either side of a one micron cutoff. A wall poked through more thinly
  would lose its hole again; telling them apart by what fills the ring
  would hold either way.

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
   A patch running out tangentially into a plane: done where the mesh
   has vertices along the line they meet on. No crease bounded such a
   region, so it took the plane with it and was no disk, or its fit failed
   and both went to facets. Smooth regions now stop at the planes the free
   triangles gather into, where a plane holds a smooth region's worth of
   triangles (32). Where the patch over the mean-value map does not verify,
   the region is charted by its projection onto the plane across its mean
   normal: the map pinches the rounded boundary of a run-out into the
   square's corners, the projection does not. The section solve finds
   nothing at a tangency and the seam is the chord threaded through the
   chain, which lies on the plane; its image on the patch is now sampled
   over each knot span alike (sampled evenly over the range, 1024 points
   over a chord of 190 spans left the image 6.2e-5 off it, per span
   8.8e-7). Measured on a lopsided twisted hill on a plate's top,
   running out with its slope (not its curvature) into the flat, meshed at
   40 cells across: before, 1153 faces (the top one smooth region with the
   plane, unverified, the hill faceted); now 7 faces, the seam a closed
   curve on the plane (within 2e-15 of it) at a tolerance of 7.7e-6
   against a distance of 2.9e-5, the volume 3.3e-5 off 2021.6, valid and
   tessellating closed. On the rough box one noisy corner, bounded now by
   planes and fillets alone, comes back a patch instead of facets. The
   truth bench and the stress baseline come back unchanged, and so do the
   corpus parts and the NIST parts in faces, kinds and volumes, but for
   thread_flank_narrower_than_a_chord (invalid before and after: 8213
   faces to 8242, 51 curved to 50). Their refusal counts move, as smooth
   regions stop at planes: fewer not disks (ctc_03 from 11 to 2), a few
   more narrow.

   A patch running out tangentially into a cylinder, a sphere or a round:
   done where the mesh has vertices along the line they meet on.
   Recognition runs first, so the canonical region already bounds the
   smooth region and the patch builds; what failed was the seam. The chord
   threaded through the chain's vertices stands off a curved surface
   between them and overshoots where the chain turns, and its image on the
   surface strays further. Measured before: a hill on a round bar's side
   (radius 5, 96 cells round, 40 along) came back 4 faces, the seam at a
   tolerance of 4.2e-4 against a distance of 2.5e-5 and 5.5e-5 off the
   cylinder; a hill on a ball (radius 10, 96 by 128) 2 faces, the seam at
   1.9e-3 against 1.4e-4 and as far off the sphere. Now, between a
   cylinder, cone, sphere or torus and a patch, the chain is threaded in
   the canonical surface's chart, cut where it turns by thirty degrees,
   and lifted onto the surface; the curve and both images come from the
   same points, and it is kept where it holds both faces within the reach
   (a loose seam keeps the chord). The bar's seam lies on the cylinder
   within 1e-9 at 4.2e-5 (the patch's own distance near the hill's
   corners), the volume 4.6e-5 off; the ball's within 7e-7 at 2.9e-4, the
   volume 3e-4 off; both valid and tessellating closed. A round along a
   block's edge running into a free-form top came back 7 faces before and
   after, its seam the round's ruling, exact. A curved region holding more
   triangles than the smooth regions it meets is no longer tried with them
   as one region: it bounds them (the ball was first tried whole, with the
   hill, as one patch). The truth bench, the stress baseline and the
   corpus and NIST parts come back unchanged in faces, kinds, volumes and
   refusal counts, but for thread_flank_narrower_than_a_chord (invalid
   before and after: 8242 faces to 8131, 50 tori to 52).

   A sphere cap with holes: done. A sphere one ring goes round and every
   other ring does not (a dome drilled off its pole, a hill on a
   hemisphere's side) is framed about an axis whose pole is on the face
   and five degrees clear of every ring, a ring's plane normal tried first
   so a flat rim is a latitude. Its seam runs from a rim vertex to the
   pole straight in the chart: a meridian where one keeps half the holes'
   mean step clear of them (a seam along a hole's side left the face
   meshing open), a slanted line otherwise; the pole is an edge of no
   length and the holes inner wires. A sphere two rings go round (a dome
   drilled through its pole off its axis) has its poles inside the two
   rims and is laid out wrapped, its seam between them. Measured on the
   hill on a hemisphere (radius 10, 48 rings by 128): 5304 faces before,
   the sphere faceted, the volume 5.9e-4 off; now 3 (plane, sphere,
   patch), valid, tessellating closed, 6.8e-8 off. At 40 by 96: 3364 faces
   to 3. A dome on a drum and a hemisphere on its base, drilled (radius
   1.5) off the pole and through it, meshed at 0.05 and 0.01: 299 to 3111
   faces before with the sphere faceted, now 4 and 3 faces, volume within
   3.2e-8. The truth bench and the stress baseline come back unchanged,
   and so do the corpus and NIST parts in faces, kinds and volumes, but
   for m5x16_bhcs, whose button head is such a cap: 1272 faces to 215, 2
   spheres to 3, volume 2.3e-3 off to 5.9e-5, valid and closed before and
   after.

   Still to do:
   - a tangent line inside a row of triangles (a scan, or a part meshed
     without vertices along its edges): the boundary row's vertices lie on
     the plane and the curvature's jump falls inside one span, which a
     cubic with knots no finer than two vertex steps cannot follow (the fit
     runs out at 7e-5 to 4e-4 against a target of 1.5e-5), and the hill's
     near flat foot gathers into narrow slanted planes beside the main one,
     each within the distance, the volume off by a few times the distance
     over the hill's area;
   - a hemisphere coarse enough to leave its hill too few vertices for a
     patch (24 rings by 64) keeps the sphere one face round the hill's
     facets. Its conversion took 7 to 10 s, nearly all in the exact volume
     of the body (each hole edge a cubic of some 160 spans, every span ten
     boundary points with a full inner integral each); now 1.4 s (2.1 s of
     CPU against 9), the volume within 1e-13 of before. A short boundary
     panel takes as few points as its error bound allows, the inner
     integrals on an analytic surface are not refined between runs, and a
     run is summed as it goes instead of held. What is left of the time is
     in planning and fitting, not in the checks;
   - a steeper run-out into a curved surface fails the fit as on the
     plane: the bar's hill at 0.8 high instead of 0.5 stays faceted, the
     fit running out at 1.24e-5 against a target of 1.23e-5;
   - cones and tori are threaded the same way, but not measured;
   - a chain of three vertices between a curved face and its neighbour:
     done, threaded with the parabola through them (a cubic needs four
     points, and the threading gave up). Carrying points onto both
     surfaces does not serve along a tangency: a point between two of the
     vertices moves along the plane far farther than the gap it closes.
     Measured (faces, curved faces faceted, volume error): ctc_05 378, 1,
     1.22e-4 to 318, 0, 1.25e-4 (a round of radius 6.35 running out into
     a plane); ctc_02 1997, 8, 4.3e-5 to 1861, 5, 2.9e-5;
     sliver_on_a_diagonal_of_the_grid 1191, 3, 4.1e-4, invalid to 1008, 0,
     5.5e-6, valid. ftc_07, the rest of the NIST and corpus parts, the
     truth bench and the stress baseline are unchanged. One edge of ctc_02 between a torus and a cylinder fitted
     to pieces of one spline blend now stands 0.21 off (45 coplanar
     distances), a line chord where no point carries onto both;
   - between vertices the fit is no closer to the true surface than its
     knots allow. Measured against the exact surface, by distance from a
     free edge in cells: the bumped sheet at 10 cells across 1.5 times the
     distance in the edge row, under 1 inside; at 20 cells 0.3 at most;
     the bowl at 10 cells 4.2 in the edge row, 2.1 in the next, under 1
     inside; the bowl at 16 cells, its distance read at 3.2e-5, up to 1.8
     inside; the tangent hill 1.1 near its top. Lowering the fairing a
     hundredfold made the edge row worse (6.3), holding the vertices to a
     quarter of the distance left the hill unverified, and knots one
     vertex step apart changed nothing: the vertices say no more, and no
     cheap change helps.
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
   and split along its ridge 201.000000000001 of 201.

   Edits now hold or say why not. `build` names each curved region it
   faceted and why (`MeshSolidReport::fallbacks`, a `RegionFallback` per
   region with a `FallbackReason`: boundary not placed, reaching past its
   triangles, turned in, body volume off, folded seam, overlap, meshing
   open). A region a merge, a split or a free fit leaves on a surface goes
   through the tangent passes for that region alone (a round between two
   planes onto the cylinder tangent to both, a fillet onto what its
   supports fix); a fit holding an axis or a radius keeps what it was
   given, and a patch takes none. A planar region of one or two triangles
   a step made is never absorbed into a curved neighbour. `fit` takes
   `SurfaceKind::Patch` (the automatic patch fit and verification, refused
   as `NotADisk`, `TooNarrowForAPatch` or `PatchDoesNotVerify`), and a
   merge or split that no plane, canonical surface, sweep or kept surface
   fits falls back to a patch where the union is a disk and it verifies
   (with `patches` on). Tests: a one-row band piece fitted as the sphere
   through both rims verifies and cannot be built beside the cylinder;
   both are named (`BoundaryNotPlaced`) and the solid is the 64-gon prism
   (6273.096981091773 of 6273.096981091878). A round drawn as two 45
   degree facets with one row lifted 1e-5 merges onto the tangent
   cylinder (radius 3 to 1e-12, where the free fit gives 2.999976, its
   axis 3.4e-5 off and crossing the faces at 8e-6 radians) and builds 480.68583470577 of 480.68583470577
   (9e-13 off). A facet split off a 64-facet cylinder stays a flat: the
   solid is 6283.027677084415 against the exact D-cut 6283.027677084466
   (the cylinder, 6283.185, without the protection). The top of a saddle
   block, a facet per triangle, merged into a patch builds 499.99970 of
   500 (3e-4 against an allowance of 1.5e-3, the distance over the base's
   area); with the patches off, `fit` of a patch builds the same.
   `solid_from_mesh` results are unchanged: the 34 corpus STEP parts give
   the same faces, curved and faceted counts and volume bits, the truth
   bench the same lines, and the stress baseline no regression. Not done: `fit`
   has no sweep kind, the tangent passes do not revisit a neighbour an
   edit changes (a ball beside a merged round keeps its fit), and item 4's
   corner and boundary graph, once there, is what an application would
   edit next.
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
   it now integrates knot span by knot span.

   Done where no parallel (or no meridian) is free: the seam that way is
   threaded through the holes. It runs straight in the chart from a point
   on the free circle the other way to the hole vertex furthest back, round
   the hole to the vertex furthest on, and on to the next hole; those
   vertices are kept as edge ends, the straight pieces are seam edges, one
   arc of each hole bounds the face along the chain and the other a whole
   turn across, and the seam the other way is the full circle through the
   chain's start. A torus whose rims go round one way only is a band
   whatever the gap they leave (a groove narrower than a quarter turn read
   as going round both ways before), and a band round the tube between
   rims that are no circles gets a seam of its own as one round the axis
   does. Measured, each one toroidal face, valid, meshing closed: three
   bores of radius 2 straight through the tube a third of a turn apart
   round the axis and a sixth round the tube (six holes over every
   parallel), 1563.1185470 of 1563.1185695 (1.4e-8); slots over the top
   on one side and under the bottom on the other (no meridian free),
   1418.14076 of 1418.14096 (1.4e-7); a groove round the axis cut by a
   second torus, 1e-15 off; a flat-sided groove round the tube, 5e-11 off.
   Still open: holes that leave no circle free either way (no free
   parallel and no free meridian) stay facets, as do holes that overlap
   along the chain at every level tried. The exact volume check of the
   six-hole ring took most of its conversion (1.5 s wall, 3 s of CPU);
   with short panels on fewer points it converts in 0.4 s (0.8 s of CPU).
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

**Converted faces whose edges fold in space.** Done for the NIST parts:
with the converter's closure check off, every one meshed at a thousandth
of its diagonal now tessellates closed (ftc_06 closes once welded, as
before). ftc_07's four slivers were huge spheres fitted to a few long
triangles of a thin fillet, two of their own edges crossing within those
edges' tolerance; `folded_seams` now refuses a curved face whose own
edges cross, and the fillet stays facets there. ctc_02's blend comes back
as tori side by side meeting at so slight an angle that where they cross
lies millimetres from the mesh's boundary, so the section solve is
refused and the edge falls to a chord through the boundary's vertices,
standing up to 1.1 mm off both tori across its long spans. Between two
curved faces such a chord is now threaded through points carried onto
both surfaces as well, and keeps within 0.07 mm. The closure check stays
as the backstop; no converter or corpus test needs it now. Every chord between two curved faces now carries points onto both, not
only one whose vertices alone stand past the reach; where the vertices
alone thread no curve within a twentieth of the longest span, the carried
points' curve is the one taken. Measured on top of the two changes below
(faces, curved faces faceted, volume error): ftc_06 345, 10, 8.5e-5 to
155, 0, 2.0e-6; ctc_02 2106, 10, 1.3e-3 to 1997, 8, 3.7e-5; ftc_07 695,
0, 1.8e-5 to 647, 0, 8.5e-6; ftc_10 513, 4 to 498, 3;
sliver_on_a_diagonal_of_the_grid (invalid before and after) 1242, 6 to
1191, 3. The rest of the NIST parts, the corpus, the truth bench and the
stress baseline are unchanged. Run alone, ftc_07 converts in about 4 s
and ctc_02 in 7 s, as before. The spheres the earlier trial faceted were
ftc_07's rounds, since fitted on their cylinders.

**Converted faces folded under a neighbour.** `check` probes every face
of a solid for orientation, stepping off it by shorter steps until both
sides answer, where it used to trust the faces' flags whenever they
agreed along shared edges. That finds planar faces in the converter's
results for ftc_10, ftc_07 and ctc_02 turned into the material within
about 0.03 of themselves: a fold, material on the face's outside and a
thin empty wedge on its inside, under a neighbour that runs past it. The
converter asks the same probe, so it catches these and facets the curved
faces beside them (ftc_07 1375 faces to 1460, ctc_02 4090 to 4138).
ftc_10 went from 704 to 977 because a threaded build with such a face
was thrown away whole for the unthreaded one, whose crossed seams then
facet 16 curved faces; the threading now stands and only the cylinder
beside the turned facets is faceted (726 faces, 144 curved, volume error
3.8e-4 from 8.1e-4).

What folds them, in every case a single-triangle facet meeting a curved
face all but tangentially, the facet a chord plane of a neighbour that
stayed faceted (a chamfer, a drafted wall, a blend):

- The exact seam wraps past the facet (ctc_02, two facets on a 57.5
  bore). The plane through the facet meets the cylinder in a long thin
  ellipse; the arc between the two corners on the cylinder bulges 8 in
  the plane for a 1.1 sag and passes 0.16 beyond the third corner, which
  sits 0.019 inside the cylinder. The face's loop then encloses the
  crescent between the arc and the facet's other two sides, wound the
  wrong way. No exact curve on both surfaces exists here: the facet lies
  wholly on the cylinder's material side.
- A chord seam leaves a slot (ftc_10 beside a radius 6 cylinder, ftc_07
  tori tangent to a drafted wall, ctc_02's third facet). The seam is a
  line threaded through the mesh corners with a tolerance of 0.09 to
  0.15, the curved face trimmed at its projection; between the line and
  that trim is an opening as wide as the tolerance. A probe point 0.003
  to 0.03 off the facet sits below the curved surface's level, and a
  ray through the opening reaches the curved face from inside: In. At
  0.1, past the tolerance, the probe answers right.

Done for single facets: a facet of one triangle with one side on a
curved face and two on planar ones is built as a fan when it is found
turned in, or when its loop as built winds the wrong way round its
plane (the seam bulging past its far corner). The fan is a B-spline
ruled from the far corner (its `v = 0` row, closed on an edge with no
length, as a cone's apex) to the seam drawn straight in the curved
face's chart and lifted onto it; every edge is shared exactly, and the
seam's tolerance is the lift's interpolation error (about 5e-7).
Measured against the build before (faces, curved faces faceted, volume
error):

- ftc_10: 726 to 687 faces, 11 to 10 faceted (the radius 6 cylinder
  stays curved; three fans), 3.8e-4 to 1.4e-4.
- ctc_02: 4138 to 4115 faces, 84 to 83 faceted, seven fans (six round
  the 57.5 bore, which stays curved, one on a cone). The volume now
  takes the closed form (before, one torus face fell to the mesh) and
  reads 1.3e-3 for 1.6e-4; the closed tessellation of the same result
  reads 1.7e-4. The difference is five R 112 r 10 tori whose chart
  integral comes out 1.8% under their fine mesh's area (2369 for 2412,
  the source triangles giving 2326), each 11500 short in volume: an
  error of the chart integration, open, not of the fans.
- ctc_04: 1213 to 1204 faces, its one Overlaps fallback gone.
- ftc_07 unchanged: its six tori turned in lie against faceted corners
  of drafted walls, runs of several facets, so no fan applies (built
  since as cones, below).

Every face meshed on its own still meets its neighbours, counting no
triangle with two corners at one point (an apex row's; a primitive
cone's apex meshes the same). Truth bench and stress baseline
unchanged; ftc_10 converts in about 7.8 s for 5.9 s (the extra
rebuilds), ctc_02 and ctc_04 as before.

ctc_02's third facet, between a cylinder and a torus (two curved
sides), is built as a wedge (below). The probe classifies two
points per face against the exact boundary, about 0.45 s on a 165-face
helical thread (0.05 s before), 0.16 s on ctc_02 (1.0 s before, when
the flag walk ran).

Done for facets whose seam does not place: a single facet beside a
curved face, a chord of a spline face left faceted, can cut across the
curved surface's curvature, its shared span standing off that surface by
just over a twentieth of its length, so no chord holds and the section
solve finds no crossing. Where the span is the facet's only side on a
curved face, the facet is built as a fan before the curved face is
faceted. Measured (faces, curved faces faceted, volume error):
grid_point_on_a_diagonal_boundary 823, 3, 6.3e-4 to 744, 1, 6.2e-4 (a
cylinder and a torus); ctc_02 1861, 5, 2.9e-5 to 1854, 4, 2.9e-5. The
other parts, the truth bench and the stress baseline are unchanged.
ctc_02's torus of radius 10 fitted to a spline blend has a facet whose two
sides lie on it; it is built as a wedge (below).
A facet whose plane meets a narrow curved face along a curve bulging past
that face's far side (a quarter round 0.17 long beside a single facet,
the seam 3e-4 over the far rim) folds the curved face; the facet is now
built as a fan before the curved face is faceted for the fold.
grid_point_on_a_diagonal_boundary 744, 1, 6.2e-4 to 741, 0, 6.2e-4; the
other parts, the truth bench and the stress baseline are unchanged.

A small region on another canonical kind than its larger neighbour joins
it where its vertices lie on the neighbour's surface within the distance,
as one of the same kind already did: a few facets of a fillet torus lie on
a sphere as exactly as on the torus, and the sphere could not meet its
neighbours. Measured (faces, curved faces, faceted, volume error): ftc_08
317, 150, 5, 2.8e-5 to 310, 151, 3, 2.8e-5; ctc_02 1854, 334, 4, 2.9e-5 to
1853, 333, 4, 2.7e-5; ctc_04 590, 273 to 588, 271 (volume unchanged);
m5x16_bhcs 215, 6, 2, 5.9e-5 to 212, 3, 2, 9.3e-5 (two spheres and a
second cylinder on its head's one cylinder join it, as in the source);
sliver_on_a_diagonal_of_the_grid 1008, 70 to 1006, 68, 5.5e-6 to 4.5e-6.
The truth bench's noisy rim-disc comes back 4 faces for 5, 4.1e-7 off for
8.4e-5 (mean volume error 3.5e-6 to 7.1e-7). A torus met only as such
pieces (ftc_08's other corner: two spheres, no torus region beside them)
has each sphere beside one plane and a cylinder or cone square to it put
on the torus tangent to both where that holds its vertices, and the
pieces then join: ftc_08 310, 151, 3, 2.8e-5 to 302, 153, 0, 2.2e-5, the
other parts, the truth bench and the stress baseline unchanged. ftc_10's
last fallback is a cylinder of seven triangles
fitted across a torus and a plane's edge (a vertex 0.07 off the torus).

A round along a plane's edge that turns round a rounded corner goes on as
a torus, and met a few facets round, that torus can come apart into
pieces whose vertices lie on spheres as exactly as on it, with no
cylinder square to the plane beside them to fix it (the corner's own
cylinder stayed facets). A sphere beside a plane and a round (a cylinder
whose axis runs parallel to the plane at its radius) is now put on the
torus the round turns on: its axis square to the plane, its tube the
round's, the circle of its tube's centres tangent to the round's axis
where the round ends, which the vertices the two regions share fix; the
major radius is solved from the sphere's vertices. Without the round's
end the points do not tell the torus from a sphere through them.
Measured (faces, curved faces faceted, volume error): ctc_04 588, 2,
1.06e-4 to 581, 0, 1.07e-4 (two spheres of radius 5.11 on a torus R 8
r 5). The other parts, the truth bench and the stress baseline are
unchanged.

A flat face much larger than a curved region's facets joins the region
whole only where each of its facets lies on the surface as the surface's
own would (sagging and leaning as a chord of it), not only its corners:
a disc capping a chamfer's cone has every corner on the cone's rim, and
the cone took it and could not close (BoundaryNotPlaced). m5x16_bhcs 212,
2, 9.3e-5 to 181, 1, 2.0e-4 (the chamfer's facets stood inside the cone
and offset the socket's faceted drill point, which stands outside its
cone); m5x16_bhcs_loops 470, 3, 2.5e-3 to 439, 2, 2.4e-3.

A drill point whose rim is arcs of one circle met by several faces (the
cone under a hex socket, its base circle touching the six walls) is laid
out as a cap as one with a single rim circle is: the arcs in turn along
the parallel's line in the chart, the seam a ruling up from the rim's
first vertex. m5x16_bhcs 181, 1, 2.0e-4 to 20, 0, 1.1e-4 (the source's
twenty faces); m5x16_bhcs_loops 439, 2 to 278, 1.

A smooth region counts as bending only on turns between its own
triangles. A flat face beside a recognized one across an edge smoother
than a crease (a screw head's flat top beside its sphere, 26 degrees)
was a smooth region of its own, was offered the sphere as a band, and
the two came back one sweep with a kink in its profile, faceted whole.
m5x16_bhcs_loops 278, 1, 2.5e-3 to 20, 0, 1.1e-4.

Done for facets with two curved sides: a single facet between two curved
faces and a planar one is built as a wedge, the B-spline surface ruled
between its two seams, each drawn straight in its curved face's chart and
lifted onto it at evenly spaced parameters so the two share their knots,
closing at their shared corner on an edge with no length; its third side
is the ruling at their far ends. ctc_02 1853, 4, 2.7e-5 to 1828, 2,
2.4e-5 (a torus and a cylinder faceted for Overlaps stay curved);
grid_point_on_a_diagonal_boundary 741 to 735 faces, 6.22e-4 to 6.25e-4;
the open thread flank 2035 to 2001 faces; the degenerate spline sliver
keeps its 13 regions faceted, one for a folded seam where it failed to
build.

A wedge whose two curved sides lie on one curved face shares one chain
with it, through their shared corner, and no curve through the chain's
three vertices holds (ctc_02's torus of radius 10 on a spline blend, the
parabola 0.24 off for a limit of 0.20). The shared corner is pinned, each
side is a seam of its own, and the facet a wedge. ctc_02 1828, 2 to 1799,
1, 2.5e-5; the open thread flank 2001 to 1916 faces.

A run of facets between two planes that are both single facets (chords
of a blend left faceted) is not put on the round tangent to them: on
ctc_02 a run of two such facets came back a cylinder of radius 5.07 that
reached past its triangles. ctc_02 1799, 1 to 1807, 0, 2.5e-5; the
degenerate spline sliver 243, 13 faceted to 228, 4 (two still for
`BuildFailed`). Over these seven changes every NIST and corpus solid
stays valid and tessellates closed, the truth bench's noisy rim-disc
moves from 4.1e-7 to 3.3e-7 and the rest of it is unchanged, and the
stress baseline is unchanged. Left: ftc_10's cylinder of seven triangles
fitted across a torus and a plane's edge, a misfit the fallback answers.

Done for ftc_07, in two steps:

- Its rounds of radius 0.43 along the drafted walls' foot are meshed
  a few rows across, with fans of long facets from single corners. A
  seed on such a strip fitted the round's cylinder exactly (1e-15) and a
  sphere of radius 2e3 to 1e5 within 1e-5 to 3e-4, and the sphere won
  the chord count (recognition prefers the fit laying most mesh edges on
  itself) by laying the strip's cross chords on itself as their plane
  does. A sphere or torus the samples tell from a closer fit (farther
  than twice it, or a thousandth of the distance) now takes no part in
  the count; a ruled fit still wins over a round one the samples cannot
  tell it from (a torus band whose two rows lie on a sphere exactly and
  on a cone at 2e-4, in sliver_on_a_diagonal_of_the_grid, stays a cone).
  ftc_07: 1460 to 1163 faces, 117 to 14 regions faceted, no sphere fits;
  ctc_05 gives one torus of radius 3e4 back to its cylinder.
- The drafted walls' corners are cones of 1 and 2 degrees meshed two to
  six facets round and one high: seven vertices, too few to fit a cone
  (nine) or for recognition to grow from. They stayed planar facets and
  the tori at their foot folded under them (Overlaps). `faceted_rounds`
  takes runs of planar facets joined across smooth edges between two
  planes, finds the rulings the run shares with each plane, and puts the
  run on the cone whose axis lies in each ruling's plane with its face's
  normal (a cylinder tangent to both where the rulings are parallel),
  verified at every vertex and against each facet's sag. A polygon whose
  flanking sides are chords too fails: the surface tangent to them along
  their edges misses its corners. ftc_07: 1163 to 927 faces, 14 to 8
  faceted (the countersinks' cones, BoundaryNotPlaced, untouched), all
  six tori built, volume error 2.2e-5 (1.8e-5 at the start), valid,
  tessellates closed, 11.0 s to 5.4 s. Elsewhere: ctc_01 223 to 220,
  ctc_02 4115 to 4104, ctc_03 201 to 169 (volume error 1.1e-4 to
  6.5e-5), ctc_04 1204 to 1184, ctc_05 658 to 648, ftc_08 346 to 328,
  ftc_09 145 to 132 (2.5e-5 to 2.2e-6). ctc_02 and ftc_08 each report
  one more faceted region: a run put on its cylinder whose seams do not
  place, built as the facets it was. Truth bench and stress baseline
  unchanged. A synthetic drafted boss (corners four facets round, a
  fillet at the foot) pins it: before, no cones and three tori folded.

Cones closing at their apex: done. ftc_07's eight cones of half angle 59
degrees are drill points, not countersinks: each region is bounded by one
circle (its bore's) and closes at its apex inside it, and only a sphere had
a layout for one rim. A cone whose one rim is a full circle starting on its
seam, with a vertex at its apex within the reach and every vertex on the
apex's nappe, is now laid out as a cap: the rim, a ruling up from the apex,
and the apex an edge of no length, the surface windowed down to the apex.
Measured (faces, curved faces faceted, volume error, conversion time with
six parts at once):

| Part | Before | After |
|---|---|---|
| ftc_07 | 927, 8, 2.2e-5, 7.1 s | 695, 0, 1.8e-5, 6.5 s |
| ctc_01 | 220, 2, 4.8e-6 | 150, 0, 4.0e-7 |
| ctc_02 | 4104, 84, 1.3e-3, 15.5 s | 2106, 10, 1.3e-3, 12.7 s |
| ctc_04 | 1184, 24, 1.5e-4, 6.0 s | 590, 2, 1.5e-4, 3.8 s |
| ctc_05 | 648, 11, 1.2e-4 | 378, 1, 1.3e-4 |
| ftc_06 | 403, 12, 8.6e-5 | 345, 10, 8.5e-5 |
| ftc_10 | 687, 10, 1.4e-4 | 513, 4, 1.4e-4 (one ReachesPast) |

All valid and tessellating closed; the other corpus parts, the truth
bench and the stress baseline come back unchanged. A blind hole with a
drill point, meshed and converted, pins it.

A rim circle need not start on the cone's seam: where it is shared with
the bore, its start comes from whichever surface placed it first. The cap
was refused there (ftc_10's two drill points, one faceted for its
boundary, the other for reaching past its triangles in the layout tried
instead). The cone is now turned about its axis so its seam runs through
the rim's start. ftc_10: 498 faces, 3 faceted, 2.2e-4 to 420, 1, 2.3e-4;
the other parts, the truth bench and the stress baseline are unchanged.
Its last fallback is a cylinder of seven triangles reaching past them.

Tori standing in for cylinders: done. ctc_03's rounds (radius 9.5 to
12.7) and two of ctc_05's walls (radius 279.4) came back as tori of major
radius 1e3 to 3e4 whose tube radius differs from it by the cylinder's: a
torus that close to its axis is the cylinder about it to within the
samples. The cylinder fit itself missed them (its axis read from a strip a
few rows wide, it stood 7 to 20 off), so no simpler fit was there to tie.
Now a torus that fits, with its tube radius (or its major radius) ten
times the samples' span or more, seeds the cylinder about its axis (or
along its tube), refined and kept where it ties the torus; requiring the
torus to fit keeps a strip across a real fillet, which neither fits, from
seeding a cylinder along it. Measured: ctc_03 7 tori to none (169 faces,
volume error 6.5e-5, unchanged), ctc_05 7 tori to 5 (378 faces, 1.3e-4 to
1.2e-4), ftc_11 4 tori to 2 tori and 2 cylinders (the source's kinds),
ctc_01 and m5x16_bhcs_loops each a torus to the source's cylinder,
ftc_08 328 faces to 317 with one more region faceted (4 to 5), and
grid_point_on_a_diagonal_boundary 826 to 823 with one region faceted for
a folded seam. The truth bench's rim-disc at 0.05 moves from 1.2e-6 to
6.3e-8 (its mean volume error 3.57e-6 to 3.50e-6); the stress baseline is
unchanged.

A face that cannot be built no longer fails the conversion. A face whose
build fails is withdrawn (for a fan, the curved face its seam lies on; for
another planar face, the curved faces beside it), a check after the build
that fails withdraws the faces that cannot be meshed or measured on their
own, or every recognized one where none is found, and a planning failure
withdraws the recognized regions; each is named with
`FallbackReason::BuildFailed`. A planar region whose plane does not hold
its vertices is gathered again from its triangles, and a planar face
withdrawn a second time is built a face per triangle; the error stands
only where nothing is left to withdraw. The corpus's degenerate spline
sliver failed outright (fans with no normal along their seam) and now
comes back, 243 faces, 13 faceted, 3 of them for `BuildFailed`. Every
other NIST and corpus part, the truth bench and the stress baseline are
unchanged, in the same time. Regions of five NIST parts and a screw put
one at a time on a surface moved off them (a plane offset, a curved
surface shifted, a radius grown or halved, through
`MeshRegions::put_surface_unverified`) all come back valid and closed,
the moved region faceted; before, an offset plane came back with faces
turned into the material. A conversion of ctc_02 that failed with "the
face's boundary enclosed no triangulable region" was not reproduced on the
current build.

The converter's last orientation probe now asks every face, as `check`
does (`inside_out_faces`), where it asked only curved faces and the faces
beside or across them: a planar face turned in among planes went
unprobed. A face it names with no curved face beside it or across its box
is gathered again from its triangles, then built a face per triangle
(`FallbackReason::TurnedIn`). No NIST or corpus part, the truth bench or
the stress baseline changes. The probe costs what it did: its time is
building the solid's boundary once, not the faces asked (ctc_02 0.29 s
to 0.33 s over its builds, ftc_10 0.64 s to 0.59 s, every other part
under 0.03 s). A box whose top is put on its plane turned over comes back
valid with the top gathered again; before, it came back turned in.

**Flags the exact volume trusts.** `volume_properties` takes the closed
form only where the faces' flags agree. The edge walks tie faces into sets
that agree among themselves; a face the walks cannot read (a cone apex, a
trim that is no loop) is a set of its own, and so is each closed shell.
Where there is more than one set, each set is probed as `check` probes a
face, against its own solid: one face of each first, and where that face
cannot tell (the broad face of a sheet thinner than the probe's shortest
step, as a cut 0.001 below a top face leaves), more of its faces at each
round, since any face of a set answers for it. A set facing in, or one no
face of which the probe settles, sends the volume to the mesh, which
mends a minority of turned faces. Placed faces (a prism's top) are
walked, an edge occurrence being its node at its placement, and the
shape's own placement is set aside. Counts: on the stress run 177 of 1129
exact volumes are probed, every set settled by its first face, none
falls back; over the workspace's tests (both tiers) 457 of 2148, one
falling back (a cone cap turned in on purpose) and one settled on a
second round (a 0.0005 sheet beside a drum). Stress time unchanged in
paired runs (14.5 and 15.1 s against 14.8 and 15.1 s). Open: a face
whose boundary middle lies outside it can still read its walk
backwards, and if it is also turned over the two errors cancel unseen.
A void turned inside out as a whole (every face of its shell reversed)
agrees with itself, and the mesh mends only a minority within each
connected piece: a 10 block with a 4 box void so turned weighed 1064 for
936, with a sphere of radius 2, 1031.4 for 966.5. A set facing in that is
a whole shell of a solid with several is now turned over in the closed
form (936 and 966.49, exact). Where such a solid goes to the mesh (a face
of the block turned as well), every face under a solid is probed as
`check` probes it, meshed on its own on the shape's agreed chords, and a
face facing in is counted turned: 936, and 966.52 at a 1e-3 chord. A face
the probe cannot settle there counts as its flag says, with no mending.
One turned face of the void or of the block, the other shell right, was
already mended. The mesh path pays one probe per face only for solids
with more than one shell.

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
    asked, an edge rebuilt whole is built once for its two faces, and each
    edge is read once for the faces on either side of it. The sew compares
    only the rebuilt pieces and the copies beside the other solid: a copy
    whose vertices only its own solid's faces hold passes through as it
    stands (`sew_around`). The classifier draws a plane's holes one at a
    time, as a point comes near one, and a hole of a crossed face that no
    other strand comes near is walked apart from the face's arrangement and
    put back in its piece. `boolean_local` in `tools/ogeom-bench` times a
    drill into the side of a plate with 225 holes and a slot across its top.
    Still whole-solid: gathering both solids, rebuilding every kept face as
    new nodes (the copies too), the history composed over all of them, the
    seam join and closure passes, and the classifier's face bounds; a
    crossed face's interior probes are still sought over all its holes. The
    pair filter by box costs nothing measurable. Profiled on a 900-hole
    plate (a side drill and a boss fuse, each touching one or two faces),
    the boolean's CPU goes to: rebuilding the kept faces 27 to 38% (copies
    included, through the general piece builder), gathering both solids
    14 to 17% (each edge sampled at 33 points for its box), the sew 8%, the
    two classifiers' face bounds 6%, the seam join and the membrane pass 2
    to 4% each, composing the history 3 to 8%. A drill into the top crosses two
    faces of 900 holes, and their arrangement is 45% of it: the scanline
    probe search over 200,000 outline points per face, and each hole
    walked on its own. No phase is above a fifth on its own, so no one
    change wins much: removing the quadratic edge lookups of a many-holed
    face, the allocations of `Model::widen` and the vertex bins, and the
    history clones when a step records nothing gave 8% on the side drill
    and the fuse and nothing measurable on the converted mesh or the top
    drill. A radix sort of the probe heights was slower than the
    comparison sort on outlines that arrive nearly ordered.
  - Unchanged faces in the history: done as `History::copy`, an exact copy
    on new nodes. Sharing the nodes themselves would remove most of the
    rebuild, a third at most, and needs: the rebuilt neighbours to name
    the shared edges and vertices (pcurves on the input's own surface ids,
    so the face keeps its placement), no weld or tolerance widening on a
    shared vertex (`Rebuild::vertex` and `Model::add_face` both widen in
    place), the sew to leave a shared face's derivation alone (it records
    every settled face as derived from itself), and every operation that
    edits its result in place (heal, fillet, offset, sew) to copy on write
    first, since the input would see the edit.
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
