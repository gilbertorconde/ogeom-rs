# Corpus manifest

One entry per committed file. A file without an entry does not get committed.

| File | Source (URL) | Licence / legal basis | Fetched |
|---|---|---|---|
| `nist_ctc_01_asme1_rd.stp` | [NIST-PMI-STEP-Files.zip](https://www.nist.gov/system/files/documents/noindex/2024/06/19/NIST-PMI-STEP-Files.zip), `AP203 geometry only/`, via the [NIST download page](https://www.nist.gov/ctl/smart-connected-systems-division/smart-connected-manufacturing-systems-group/mbe-pmi-0) | US government work, 17 U.S.C. §105; page states verbatim: "The test cases, CAD models, and STEP files can be used without any restrictions." NIST appreciates acknowledgement; NIST logo not to be used promotionally. | 2026-08-04 |
| `nist_ctc_02_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ctc_03_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ctc_04_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ctc_05_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_06_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_07_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_08_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ftc_09_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_10_asme1_rb.stp` | same | same | 2026-08-04 |
| `nist_ftc_11_asme1_rb.stp` | same | same | 2026-08-04 |
| `NIST-README.txt` | NIST's own readme from inside the zip, kept verbatim as provenance | same | 2026-08-04 |
| `nist_ctc_01_asme1_ap242-e1.stp` | Same NIST zip as the AP203 parts, `NIST-PMI-STEP-Files/` root: the CTC 1 part in AP242 edition 1 with full semantic PMI — dimensions, plus/minus bounds, geometric tolerances, datums. | Same NIST terms as above | 2026-08-05 |
| `ogeom_asm_bolted_plate.stp` | Authored for this project (generated, then committed): a plate with two bolts — AP214 product structure, three usage occurrences over two parts, per-product and per-face colours, reference designators. Ground truth for the assembly reader, every value chosen by hand. | This repository's own licence (MIT OR Apache-2.0) | 2026-08-05 |
| `threemf_cube_deflated.3mf` | Authored for this project (generated with Python's `zipfile`, then committed): one 10 mm cube with a base-material colour, placed at x = 10 by its build item, its start part named by the package relationships, every entry deflated and streamed with a data descriptor as slicers write them. | This repository's own licence (MIT OR Apache-2.0) | 2026-09-23 |
| `threemf_cube_zip64.3mf` | Authored for this project: the parts of `threemf_cube_deflated.3mf` repacked as streaming writers lay out a package whatever its size — ZIP64 end record and locator, sizes and offsets saturated in every header and carried in ZIP64 extra fields, entries at version 4.5. | This repository's own licence (MIT OR Apache-2.0) | 2026-09-23 |
| `threemf_components_inch.3mf` | Authored for this project, as above: model unit inch; one object of two components of a one-inch cube held in a model part of its own and reached through the production extension's `p:path`, the second component mirrored in x. | same | 2026-09-23 |
| `threemf_painted_sphere.3mf` | Authored for this project, as above: a 4 680-triangle sphere painted two colours from a colour group, large enough to span several dynamic-Huffman blocks, beside a cube of type `support`. | same | 2026-09-23 |

The eleven `.stp` files are the **AP203 geometry-only** exports of the CTC 1–5 and
FTC 6–11 test parts from the NIST MBE PMI Validation and Conformance Testing
project — part geometry with no PMI, which is the subset a geometry kernel reads
first. NIST's readme states plainly that these are **not** error-free reference
files: conformance checkers report syntax errors in them. That is part of their
value — a reader that only accepts clean files has not been tested. The
PMI-annotated AP242 variants exist in the same zip and can be added under the
same licence when PMI is in scope.

## m5x16_bhcs.step

An M5x16 button-head cap screw, five placed bodies, extracted standalone
from a community printer assembly (issue #37's fixture). Nothing but the
classic primitives — spherical button head, cylindrical shank, conical
chamfers — with the head's sphere zone slit along a meridian that sits at
the chart's own seam: the minimal reproducer for a doubly-used edge that
is a slit, not a period-wrapping seam.

## m5x16_bhcs_loops.step

The same M5x16 button-head cap screw as `m5x16_bhcs.step`, one body, but
extracted after the reader stopped slitting it: the head's sphere zone
keeps the file's own two closed rims — circles cut square to the screw on
a sphere whose chart runs along z, nested loops in the chart rather than
parallels of it. The reproducer for a two-ring periodic face that is a
face on its own bounds, not a band missing its seam (issue #37's
assembly-side half).

## nema17_coupler_faces.step

Six of the seventy-three single-face surface bodies a NEMA 17 motor
coupler was exported as, extracted standalone from the same community
printer assembly as the M5x16 screws: `SHELL_BASED_SURFACE_MODEL` over an
`OPEN_SHELL` of one `ADVANCED_FACE` each — two planes, two cylinders, two
tori — hung off one product's shape representation. The reproducer for a
part the file never calls a solid: a reader that walks only
`MANIFOLD_SOLID_BREP` leaves it invisible.

## nema17_coupler_hole.step

One more of the same coupler's surface bodies: a half-cylinder wall with
a two-edge hole loop straddling the drum's chart seam. The reproducer for
a hole whose two halves arrive on different branches of a periodic chart
— a loop that never closes in the chart and a hole the mesher draws but
cannot cut.

## eccentric_ellipse_edges.step

Eight faces round one corner of a community printer carriage, each face's
transitive closure extracted standalone and wrapped in a shell and a solid
of their own: planes and a cylinder bounded by short arcs of ellipses six
and a half metres by 1.8 millimetres, where a plane cuts a drum almost
along its axis. Their vertices sit a couple of microns off those curves.
The reproducer for a vertex inverted onto an eccentric ellipse by the
closed form alone, which read it nine millimetres along the curve and sent
the edge the long way round.

## spline_face_fit_runs_away.step

Sixty entities, one face's transitive closure extracted standalone from a
community printer assembly, with a shell and a solid wrapped round it so
the reader yields something: a degree 3×3 patch whose `u` direction is degenerate the whole
way across — a sliver four microns wide and a tenth of a millimetre long —
and the three spline trims that bound it. The reproducer for a fit with no
sound sample to lean on (issue #41).

## box_with_a_cavity.step

Authored for this project (generated, then committed): a 10 mm cube with a
4 mm cube hollowed out of its middle, written the way a real exporter
writes one — `BREP_WITH_VOIDS` over an outer `CLOSED_SHELL` and an
`ORIENTED_CLOSED_SHELL` naming the cavity the other way round. The
reproducer for a solid a reader matching on `MANIFOLD_SOLID_BREP` alone
never sees: three bodies of a community printer assembly are written this
way, and every one of them read as nothing at all. Its volume is exactly
1000 − 64 = 936 mm³, which is the assertion.

## sliver_face_falls_apart.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
quarter-arc forty-five millimetres long and eighteen microns wide, between
two nearly concentric circles. At a tenth of a millimetre the sagitta of
each bounding arc is twenty-nine microns — wider than the region — so the
inner polyline crosses the outer one, and the triangulator answers with
sixteen triangles in fifteen disconnected pieces. The reproducer for a
boundary drawn too coarsely to be a boundary.

## cone_apex_quarter_turn.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a cone
sector bounded by two rulings into the apex and one arc, the rulings
standing exactly a quarter turn apart. The run between them along the
apex row is boundary the ring needs — dropping its far end cuts the corner
through the face and loses the triangle at the apex. The reproducer for a
degenerate row that a quarter-period test read as no row at all.

## closed_tube_face_crosses_its_join.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a face
on a B-spline tube that closes on itself in `u` without being periodic,
whose trim crosses the join twice. The reproducer for a ring that must
fold across a *closed* chart's join — a fold that engaged on periodicity
alone never fired, the ring jumped the width of the chart twice, and the
face drew as six pieces.

## slit_loops_are_not_holes.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a plane
with twenty inner loops, five of which are slits — out along two arcs and
back along two splines fitted to the same arcs, three millimetres long and
a fifth of a micron wide, enclosing nothing. The reproducer for an inner
loop that is not a hole: read as one, the triangulator drew the face with
twelve holes it does not have.

## long_bore_between_cross_holes.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a bore
2.1 mm in radius and four hundred long, crossed by holes wider than
itself, so its one wire winds from rim to rim along the intersection
curves. The reproducer for a grid one row deep: a cylinder never sags
along its axis, so sag gave the bore a single interior row, and the
Delaunay triangulation bridged two hundred millimetres from each rim to
it with triangles a quarter turn wide — the bore drew as a square between
its holes.

## chart_far_from_its_origin.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
cylinder twenty millimetres tall whose axis point the file placed half a
metre away, so the face's chart spans `v` from −500 000 to −499 980. The
reproducer for a degenerate-triangle scale taken from where the ring sits
rather than how far it reaches: at that scale a quarter of a chart unit
was "degenerate", every cell of a finer grid was, and the face drew as
two triangles.

## annulus_narrower_than_its_rims_sag.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
plane annulus forty microns wide between rims of 2.845 and 2.805 mm. At
half a radian of angular deflection each rim is a sixteen-gon sagging
fifty-five microns, the two polygons cross, and the boundary encloses
nothing. The reproducer for an empty first pass refused outright instead
of drawn again finer, as a first pass that came back in fragments is.

## sliver_on_a_diagonal_of_the_grid.step

One solid extracted standalone from a community printer assembly: a
turned part of a hundred and fifty faces, among them a three-edged
B-spline patch whose bottom row collapses to a point. Meshed whole at
half a radian, three grid points on that patch sit on a diagonal, the
middle one a rounding off the line, and the sliver they make has its
centre at that middle point to the last bits; the sag repair inserted the
centre round after round, each a hair on the last, and the degenerate
filter dropped the hairs and left a hole. The reproducer for a repair
point that lands on a vertex already there. The whole solid, because the
grid that lines up is the one the solid's own edge refinement produces,
and the patch alone does not.

## grid_point_on_a_diagonal_boundary.step

One solid extracted standalone from a community printer assembly: a
turned part of two hundred and fifty faces, among them a B-spline patch
and the torus across one of its edges, an edge that runs diagonally
across the patch's chart. Meshed whole at half a radian, a grid point
falls exactly on that segment — the midpoint of two grid corners the ring
joins — and even-odd counting calls it inside; inserted, it split the
constraint on the patch alone, and the torus was drawn to the unsplit
edge. The reproducer for a T-junction born inside one face. The whole
solid, for the same reason as the one above.

## slots_cross_at_half_a_radian.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
cylinder wall with six slanted slots through it, each slot's two sides
different curves between the same two points. At half a radian of angular
deflection the two sides' polylines cross, and the face — with hundreds
of interior points — drew with two holes more than it has, while the
count of triangles against boundary points that was meant to catch a
crossing saw nothing: the crossing cost a handful of triangles and the
interior points buried the difference. The reproducer for a crossing that
only the boundary, triangulated on its own, can be asked about.

## a_spike_on_the_boundary.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
plane whose one wire runs out along an edge to a point and straight back
along a second edge over the same curve — a spike into the face that
bounds nothing. It triangulated to two hairs and one triangle too many,
and read as a crossing by the exact count. The reproducer for a spike
stripped from a ring before it is triangulated.

## closed_edge_vertex_off_the_seam.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
B-spline patch bounded by four edges, one of them a closed loop — a fitted
cubic of thirty-four controls that ends where it begins — whose one vertex
sits on the curve 2.18 mm along it from the seam. The reproducer for a
closed edge held to a seam its vertex misses: the reader widened the
vertex's tolerance to 2.18 mm to say so, and a solid meshed whole then
welded its borders at that reach.

## fillet_strip_with_a_narrow_chart.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: a
fillet strip a centimetre long and a couple of millimetres wide, turning
a quarter turn across its width, on a rational B-spline whose chart runs
`u` over a fiftieth of a unit and `v` over one. The reproducer for a
triangulation Delaunay in the chart but not on the surface: with the
chart sixty times narrower one way than the other, points were joined
along the strip across columns rather than to the row beside them, and
the triangles lifted folded — half the face shaded as creases.

## thread_flank_narrower_than_a_chord.step

One face extracted standalone from a community printer assembly, with a
shell and a solid wrapped round it so the reader yields something: the
flank of a machine thread, a helical strip a tenth of a millimetre wide
running twenty turns round a rod, on a B-spline with a long control net —
which is why the file is large. The reproducer for a face narrower than
the chord it is drawn at: its boundary sagged by three times the face's
width between points, and every triangle across the width stood off the
surface by that sag, so the thread shaded as a saw of fins.
