# Adversarial review, 2026-09-30

An adversarial review of the kernel at `ed3b155`, eight areas in parallel. Every item below was reproduced; the repro programs live outside the repository, under `/var/tmp/og-review/<area>/`. Items are worked in order and marked done with the commit that fixes them.

Status: `[ ]` open, `[x]` done (commit), `[-]` dropped (reason).

## Issues

- [x] **#101.** `fix_shape` turned a located face right way out on a new node with no location, so a void placed by a location (an assembly part read from STEP) moved out of its solid. The turned face now keeps its location. (12ef0e3)

## Tier 1: wrong results in everyday use

### M. Mirrors and negative scales

- [x] **M1.** A negative uniform scale is classified `Scale`, so `apply_direction` leaves directions unflipped (`ogeom-math/src/transform.rs`, also `Transform2`). Circles, planes and frames transformed by `scaling(c, -2)` keep the wrong orientation. (b36ce02)
- [x] **M2.** `RevolutionSurface::transformed` under a mirror evaluates at `M·R(-u)·c(v)`: the axis must be reversed when the transform does not preserve handedness (`ogeom-geom/src/surface.rs`). (b31289e)
- [x] **M3.** `OffsetCurve` and `OffsetSurface` land on the wrong side after a mirror: the offset distance must be negated when handedness flips (`curve.rs`, `surface.rs`). Planar offsets too; `Transform2::preserves_handedness` no longer counts a negative scale, a half turn in the plane, as a reflection. (b31289e)
- [x] **M4.** `OffsetSurface::analytic` rebuilds a plane on a right-handed frame, flipping a mirrored plane's normal; the cone case shifts the parameterization (`surface.rs`). (8d9ee48)
- [x] **M5.** Converting a mirrored cylinder or cone to a B-spline flips its normals: `section_at` rebuilds right-handed frames (`ogeom-geom/src/convert.rs`). (bd8e604)
- [x] **M6.** Primitives in a left-handed frame (`Frame::mirrored`) come out inside out; `make_sphere` fails (`ogeom-algo/src/primitive.rs`: box, cylinder, torus, cone, wedge, sphere). (3c9ee72)
- [x] **M7.** Booleans with mirrored curved solids: a mirrored cone gives wrong fuse, cut and common that pass `check`; mirrored cylinders, spheres and tori are refused (`ogeom-bool`). (46b2484)
- [x] **M8.** glTF reader ignores a mirroring node transform: winding is not swapped (`ogeom-io/src/mesh_formats.rs`). (ddfd682)

### S. Spheres through STEP and IGES

- [x] **S1.** A sphere face bounded only by its seam walked both ways (what our own writer emits) reads back with no triangulable region: whole spheres come back empty, spherical voids are lost, cut spheres keep the wrong side (`ogeom-io/src/step/read.rs` `face`; IGES the same). (eb6f270, 2fc8ba0)

### I. Intersections

- [x] **I1.** Closed sections of two crossing cylinders (an oblique pipe tee) come back 1 to 9 mm off both surfaces with a stated tolerance of 25 to 1200: the loop's fit fails (`ogeom-intersect/src/section.rs` `marched`, `approx.rs`). (c7eaa49)

### B. Splines

- [x] **B1.** `KnotVector::reversed` and `reparameterized` clamp knots outside the domain, changing unclamped and periodic splines by up to 9 mm (`ogeom-math/src/knots.rs`). (e2b4149)
- [x] **B2.** `to_bezier_segments` and `elevate_degree` assume clamped ends; wrong for unclamped splines (`ogeom-math/src/bspline.rs`). (bc456d3)
- [x] **B3.** `KnotVector::new` accepts end knots repeated `degree + 2` times; evaluation there gives NaN (`knots.rs`). (e30262b)

### O. Offsets and shells

- [x] **O1.** `offset_shape`, `offset_faces`, `move_faces` and `make_thick_solid` ignore a solid's Location: wrong geometry that passes `check`, or refusals (`ogeom-offset/src/shape.rs`). (572f441)
- [x] **O2.** `offset_shape` and `make_thick_solid` return inside-out solids when the offset passes through the part, instead of refusing (`shape.rs`). (2396893)

### F. Fillets

- [x] **F1.** `fillet_edge` never checks that the ball fits: a radius of 100 on a 2 mm box deletes the solid and reports success; concave edges add material far outside (`ogeom-fillet/src/fillet.rs`). (5074595, 0085b85)
- [x] **F2.** `fillet_edges` accepts blends whose setbacks overlap across a shared face (`fillet.rs`). (28f9106)

### N. Boolean fallbacks and small gaps

- [x] **N1.** The nesting fallback samples at 0.1 mm and accepts a solid that pokes 2 µm through as nested (`ogeom-bool/src/lib.rs` `or_nested`, `nested`). (da82038)
- [x] **N2.** Faces closer than about 1.5e-5 are welded into broken topology that passes `check` (weld floor `confusion · 1e2` in `assemble_result` and `outline_snap`). (7c09dfe)

### A. Other algorithms

- [x] **A1.** `make_revolution` of a face whose hole winds the same way as its outline adds the hole's volume (`ogeom-algo/src/sweep.rs` `revolution_over_face`). (c834382)
- [x] **A2.** `to_nurbs` and `general_transformed_shape` extend a pointed cone's window past the apex: 15 to 20% too much volume (`convert.rs` `bounded_to_face`). (bf5d6bd)

## Tier 2: narrower wrong results, and crashes

- [x] **C1.** `classify_in_solid_exact` says In for points in a spherical void: `shape_bounds` of a reversed whole sphere bounds only its seam (`ogeom-algo/src/measure.rs` `patch_bulge`). (17e0a28)
- [x] **C2.** Whole spheres and tori from f32 meshes fail `check`: `whole_face` widens the face's tolerance but not its seam, pole edges and vertices (`mesh_solid.rs`). (7bc5070)
- [ ] **C3.** A small valid `coplanar_distance` on f32 input gives a shape whose whole-shape tessellation is open while `check` passes (`mesh_solid.rs`). Diagnosed 2026-10-01: the slits run along recognized sphere faces (a rounded block's corners) where they meet single-facet planar faces; those sphere faces' triangles overrun their exact area. Withdrawing recognition on an open body tessellation in `body_culprits` was tried and rejected: inside the build loop sound bodies also tessellate open. Deferred to the end of Tier 2. Narrowed 2026-10-01: the overrunning sphere faces have self-intersecting trims. On the rounded block's corner sphere, the edge it shares with a neighbouring fillet cylinder is a fitted B-spline that dips 0.056 rad past its own corner vertex and crosses the next edge, which runs along the corner's planar facet, so the exact area (0.096) and the triangles (0.171) disagree. The pcurve follows the 3D curve to 1e-8, so the overshoot is in the fitted edge itself. The sphere and the cylinder meet tangentially along a great circle, which has a closed form, but from f32 input the recognized centre sits off the recognized axis and the radii differ, so the closed form does not apply and the edge is fitted along a tangential seam, where a fit is ill-conditioned. The remedy is to make neighbouring recognized surfaces consistent (centre onto axis, radii equal) before their shared edges are cut. Tried 2026-10-01: rejecting a section that strays from its chain by more than the chord's allowance (5% of its longest span), and threading the fallback curve through points midway between the two surfaces, closes the tessellation (0 slits, volume 3838.3 against about 3839.5). It also changes the rough rounded box: about 170 seams become curves loose by about 0.005 beside single-facet planes, the boolean then collapses the drill's short sections across those facets as dust, and every drill in `holes_beside_a_rough_corner_cut_and_fill_valid` fails to close. Absorbing lone facets into a neighbouring curved face was tried as well; it absorbs nothing on the f32 block and breaks `a_windowed_tube_s_wall_patches_each_become_a_face`. Parked until the boolean keeps short sections across small facets bounded by loose edges. Tried again 2026-10-01 after the boolean's facet fixes (summed crossing doubts, edges shared by planes at small angles): the rough box's drills still fail with the stray filter. Two other routes were tried and dropped: seating each corner sphere on its fillet cylinders (centre to their axes, radius theirs; the spheres move by microns and the tessellation gains two slits), and flooring a seam's reach at the data's single-precision noise when the asked coplanar distance is below it (no change). Measured 2026-10-02: the rounded block tessellates closed at a coplanar distance of 8e-7 and above (at 1e-5 it comes back 6 planes, 12 cylinders and 8 spheres) and still opens at 1e-7 and 5e-7, against a measured flat noise of 3.1e-7. The open tessellations were facets the solved seam folds: a single facet beside a fillet thinner than the curve its plane meets the fillet along bulges, or a sliver between two corner balls whose two seams cross in its plane, so the facet's triangles cover its fold twice or reach from one seam to the other. Those seams are now threaded straight through the facet's corners, carrying their stray as tolerance (06c7dd2, 1c30ea0); the rough rounded box, Body11, shelf_bracket and 77777_1 tessellate closed with their curved faces kept. Body28 still opens: threading its seams and faceting what follows leaves a face collapsed onto a line elsewhere, so its threaded build is set aside. Below 1.6 times the noise the remaining slits run between curved faces, where no facet lies to thread. Planned 2026-10-02: fillets and corner balls built from their supports, and seams traced from solved corners (`docs/PLAN.md`, Mesh conversion items 2 and 4), in place of further per-case fixes. Item 2 landed 2026-10-02 and leaves 1e-7 and 5e-7 as they were: below the data's noise the regions break into 625 and 189 faces before any fillet is recognized, so nothing is derived. A coplanar distance below the measured noise asks for more than the data holds; refusing it by name, or flooring it at the noise, is the open decision.
- [x] **C4.** `solid_from_mesh` accepts NaN, negative or zero `coplanar_distance`, and NaN or zero `crease`, silently (`mesh_solid.rs`). (2ca79a1)
- [x] **C5.** The converter's nesting ray counts a hit on a shared edge twice or not at all, turning a void into a separate solid (`mesh_solid.rs` `inside`). (bdb0a72)
- [x] **C6.** Flat caps of a fine f32 mesh, rotated and far from the origin, break into many faces: `coplanar_angle` ignores the quantum (`mesh_solid.rs` `coplanar_groups`). (9ca5089)
- [x] **R1.** PLY reader panics on an out-of-range face index (`mesh_formats.rs`). (676d426)
- [x] **R2.** Deeply nested glTF (JSON) and VRML overflow the stack and abort the process (`json.rs`, `vrml.rs`). (78a7599)
- [x] **X1.** Spindle tori: the folded branch is dropped by the plane, coaxial cylinder and coaxial torus closed forms (`ogeom-intersect/src/surface.rs`). (5cd70b7)
- [x] **X2.** Circle-circle overlaps across the periodic seam are lost or reported outside the window (`curves.rs` `clipped_to_windows`). (198e86e)
- [x] **X3.** Cones whose height range spans the apex: the far nappe's section is dropped by the coaxial cylinder and parallel cone closed forms (`surface.rs`). (71d3bc5)
- [x] **X4.** A fitted section near a cone apex breaks its stated tolerance about 200 times (`approx.rs`). (25dfc19)
- [x] **X5.** A tangential contact curve (torus in sphere) is 100 times off its stated tolerance and split in two (`section.rs`). (95196d7)
- [x] **X6.** The general curve/surface path reports each tangency up to 9 times (`curve_surface.rs` dedup radius). (c4c1db8)
- [x] **X7.** Planar curve/curve misses tangent and near-tangent crossings when a curve is trimmed (`curves.rs`). (3cf6676)
- [x] **K1.** `cone_parameters` is wrong past the apex and does not return the nearest point; `torus_parameters` does not round-trip on a spindle torus (`ogeom-math/src/elementary.rs`). (bc7c37b)
- [x] **K2.** A cone with a negative half angle is accepted but its `distance_to` is wrong (`quadric.rs`). (356a0a7)
- [x] **K3.** The circle-circle bisector returns the wrong branch when the first circle is smaller (`construct2d.rs`). (8add0f5)
- [x] **K4.** A near-circular ellipse (minor slightly over major) gives NaN (`conic.rs`). (7ad169b)
- [x] **K5.** Apollonius circles lose precision away from the origin (`construct2d.rs`). (fbc3008)
- [x] **G1.** A trimmed pcurve over a non-line basis is not rescaled by a scaling transform (`curve2d.rs`). (12d754e)
- [x] **G2.** `CurveOnSurface::transformed` ignores scaling of length-type surface parameters (`curve.rs`). (12d754e)
- [x] **G3.** Scaling a parabola does not rescale its parameter (`curve.rs`; `convert.rs` range adjustment misses it too). (12d754e)
- [x] **G4.** `HelixCurve::arc_length` ignores taper (`curve.rs`). (c11773a)
- [x] **G5.** A periodic B-spline reports `is_closed() == false` (`curve.rs`, `curve2d.rs`). (bb054ad)
- [x] **L1.** `curve_length` fails at a kink of a polyline spline, and `parameter_at_length` / `points_by_count` then return wrong points silently (`ogeom-algo/src/length.rs`). (922ab2f)
- [x] **L2.** `medial_axis` ignores the face's placement (`medial.rs`). (8283b48)
- [x] **L3.** `tight_bounds` can be up to 0.045 mm too small on spheres and tori in general frames (`tight.rs`). (6e5e189)
- [x] **L4.** Mass properties of a scaled shape are off by about 2%: the tessellation's chord scales with the transform (`mass.rs`, `ogeom-mesh`). (77e152f)
- [x] **T1.** Tessellation interiors sag up to twice the chord while `deflection_met` says true (`ogeom-mesh/src/triangulate.rs`). The flag is now honest; the mesh is unchanged. Holding each side of a doubly curved surface's grid to half the chord meets the chord in the interior, but measured it triples tessellation time (`tessellate_part` 12.6 to 38.3 ms, `tessellate_torus` 2.4 to 6.8 ms), so that trade is left to decide. (f013022)
- [x] **H1.** `divide_by_angle` and `divide_by_area` fail on cylinders shorter than 1 (`ogeom-heal/src/divide.rs`, probe at parameter 1). (77b5d92)
- [x] **H2.** `to_nurbs` / `to_bezier` of a torus leaves a face whose exact measures are off by 2 to 5% and that fails to triangulate 1000 mm out. (27db634)

## Tier 3: refusals and minor contracts

- [x] **P1.** Boolean refusals: near-coaxial cylinders (radius 5e-6 to 3e-3 apart), boxes offset 1e-6, a box rotated 1e-7 to 1e-4 rad on another, cut of two identical spheres, cylinder seam at some angles near a section line. Settled by the rule agreed with Gil on 2026-10-01: a sliver thinner than the weld distance (1e-5, or the faces' stated tolerance where wider) is welded, anything thicker kept as real geometry. Fixed: cut by an identical solid and the seam-angle refusal (f5f4ba7); drums 1e-3 apart (b29d3ae); faces within tolerance of the other solid's (51f131c); a cavity touching its wall (8b84617); drums 2e-5 to 1e-4 apart in radius or axis, read against exact boundaries and exact circle crossings (e09a500); slivers under the weld distance welded by a retry held to the weld distance (836f8c5); blocks tilted 1e-5 to 1e-4 rad, near-parallel planes crossing on their true line (d88c1f7). The near-degenerate sweep (`b3`, `p1`) now answers every case with no wrong result, the box slid exactly 1e-6 included, drums 1e-5 to 1.07e-5 apart in axis included (64f5b6f). The pad on the unrefined STL conversion (`facetpad`, 2 of 6 fuses) is fixed (2893d1a): a pad wall and the facet under its top edge crossed on a solved line a sliver below the edge they share, and facet-row edges lying all but in a wall's plane were split twice by the sections either side. (f5f4ba7, b29d3ae, 51f131c, 8b84617, e09a500, 836f8c5, d88c1f7, 64f5b6f, 2893d1a)
- [x] **P2.** Fillet refusals: a pipe tee junction loop, all edges of an L-bracket, a full round of a box. (ce954b0, 1cf4e3d, 1d9ba1b, b0b672d)
  - Full round of a 2 x 2 x 1 box at r = 0.5: fixed (b0b672d). The blends on each side already met along its middle; what refused was each short edge, which the balls of the corners at its two ends consume whole. Such an edge is now rounded by its corners, and the result has Steiner's volume, 1 + pi/2 + pi/6.
  - All 18 edges of an L-bracket at r = 0.1: fixed (ce954b0). One order failed because any block corner of three blends refused at r <= 0.1: the ball patch's rim, a trimmed circle, was held to the loose overlap width for fitted curves and read as running along the sharp edge it leaves tangentially. The other order failed because the corner tool took the step's re-entrant corner for a convex one; it now probes beside each corner edge.
  - The junction loop of a drum tee (r 1 and 0.6): fixed (1cf4e3d, 1d9ba1b). The blend already runs the whole loop from any one arc; its leg on the branch is a band whose two rings started on one chart column a period apart, so the seam wound a full turn as a helix. The arcs after the first then refused as consumed; an arc on a crease already rounded is now done.
- [x] **P3.** `make_prism_tapered` rejects a face straight from `make_face` (no pcurves attached first). (acb9aed)
- [x] **P4.** `make_polyhedron` takes a face's normal from its first three points and refuses a face starting with three collinear points; `medial_axis` refuses a polygon with a collinear vertex. (79b9658)
- [x] **P5.** `elevate_degree` leaves the knots it introduced (C0 where C2), contrary to its doc. (ba5ce5f)
- [x] **P6.** Quadratic and cubic solvers lose double roots; the doc's closed-form claim for quartics is untrue. (4354c03)
- [x] **P7.** `lines_tangent_to_two_circles` misses the third tangent of touching circles. (fd3ffab)
- [x] **P8.** Derivatives above the implemented order come back zero instead of refused (circle, ellipse, hyperbola, parabola, helix, 2D pad). (6fb58f6)
- [x] **P9.** `Triangulation::welded` panics on a mesh with no normals. (e3e6fea)
- [x] **P10.** HLR `project` draws every outline edge twice. (ee9c81a)
- [x] **P11.** `make_wedge` doc says a zero top extent is refused; it is accepted. (e8e0b44)
- [x] **P12.** `volume_properties` of a plain box takes seconds at 1000 mm and more at 3000 mm. (9ed3c96)
- [x] **P13.** Minor tolerance units: `Quaternion::between` near antiparallel, `SurfacePoint::normal` compares a squared length to a length. (7042e86)
- [x] **P14.** `make_periodic` records face history only at the solid level. (4cf4820)
- [x] **P15.** Intersector 2D edge cases: coincident 2D circles ignore phase; collinear segments touching end to end return nothing; a cone tangent to a plane along a ruling reads `Apart`. (5ff2437)
