# ADR 0001: One blend model for fillets and chamfers

- **Status:** proposed, for agreement before implementation (#902).
- **Scope:** `ogeom-fillet`: edge rounds, chamfers, asymmetric and variable
  sections, vertex blends. Face–face blends (`blend_faces`) and sheet rounds
  follow later and are out of this record's first steps.

## 1. Context

`ogeom-fillet` is 40 files and about 51,600 lines, with about 650 refusal
sites. It has had 187 commits since 2026-10-01, 77 of them `fix(fillet)`.

- **Wrong solids:** in the last two release cycles it has built wrong solids
  that `check` passed: #713/#839, #846, #868/#871, #843.
- **Order and radius dependence:** results depend on the order the edges are
  asked in, and on the radius (#838, #809, #878).
- **Speed:** it is slow (#870).
- **Refusals:** it refuses simple parts. The stepped bracket in #902 is one:
  chamfering every edge gives a concave bevel hanging past both side faces,
  unbevelled side profiles and inconsistent vertex facets, and filleting
  every edge doesn't build.

What the code does today. Each fact is checked against `main`; files are
under `crates/ogeom-fillet/src/`.

- **Edge by edge, against a changing solid.**
  - `fillet_edges_own` (`fillet.rs`) builds the bands of edges that share a
    vertex one after another, each against what the last one left. Its
    comment says so: the order the edges were asked in decides everything
    else.
  - Chamfers apply one boolean per wedge in request order
    (`chamfer_edges_with_own`, `chamfer.rs`).
  - The request is sorted into a canonical order only when a signed or cusp
    corner exists (`fillet_edges_meeting`).
- **Ends and corners as a catalogue.** Each configuration has its own
  construction, its own refusals and its own cell cut and re-sew:
  - corners: setback, mixed and signed (`setback.rs`), cusp (`cusp.rs`), one
    edge left sharp (`sharp.rs`), ball-and-block (`corner.rs`),
    chamfer-facet (`chamfer_corner.rs`);
  - junctions where creases turn (`junction.rs`);
  - a ball on a far edge (`rolled.rs`), a staircase (`stepped.rs`), a contact
    changing kind along the edge (`carried.rs`), a strip taken off whole
    (`consumed*.rs`), a section shrinking to a pole (`pinched.rs`);
  - cut-off-flush ends (`cut_off_bands`, `fillet.rs`).

  A configuration outside the catalogue is refused, or built wrong where a
  case is misapplied.
- **Many rebuilds per request.**
  - A request costs one corner boolean, at least two booleans per "round" of
    vertex-disjoint edges, about ten per junction corner, one cut and one
    full re-sew per special corner, and two booleans per cut-off.
  - Retry ladders multiply this by up to six: exact rims, then fitted rims;
    crossed junctions, then flush, then none.
- **Validation after the fact.** Corner patches are fitted by iterative net
  refinement, up to 135 controls, with hold passes and 202-control nets
  (`fill_n.rs`). Each is then measured against a region and refused if it
  strays. Soundness rests on those checks being complete, and they have not
  been.
- **Planar chamfers share the curved machinery.**
  - Planar chamfers use the shared wedge-and-boolean scaffolding.
  - Where a bevel meets a round, they go through a cell cut, an N-sided G1
    patch and a re-sew (`blend.rs`, `CornerEnd::Bevel`).
  - A chamfer of a polyhedral part is an exactly solvable plane problem, yet
    it can be slow, flaky or refused.

## 2. Decision

Replace the case catalogue with **one blend model**, built in fixed phases
from the **original solid and the whole request at once**:

```
request ─► P0 normalise & classify ─► P1 sections & contact curves (per edge, on the original solid)
        ─► P2 setbacks & ends (per face corner, by intersection) ─► P3 vertex blends (one method)
        ─► P4 one assembly (trimmed faces + strips + vertex patches, sewn by known adjacency)
        ─► P5 independent guards (check, region, overhang)
```

The model rests on one topological fact, the **face-setback structure**,
shared by every blend kind:

- Every original face shrinks to the region inside its **setback curves**:
  the contact curves of the blends on its edges, or the edge itself where
  the edge is not blended.
- Every blended edge becomes a **strip** between its contact curves on its
  two faces.
- Every vertex where blended edges meet becomes a **vertex blend** bounded
  by the strips' end sections and the faces' setback corners.

This is the polyhedral bevel and the setback vertex blend of Várady and
Rockwood [5, 6], made the backbone of the whole operation. The topology of
the result is known from the request alone, before any geometry is
computed. That pseudo-topology is the same idea as the SMLib fillet
executive's "fillet brep" built ahead of geometry [14].

**Exact polyhedral chamfers are the first phase delivered.** They prove the
backbone at no geometric risk (§5).

## 3. Semantics: what a request means

A request is a set of `(edge, section)` pairs. A section is one of:

- a circular round of radius `r`;
- a conic round with distances `d1`, `d2` and a shape parameter (#877);
- a chamfer with distances `d1`, `d2`, or a distance and an angle;
- any of these varying along the edge under a law.

The result is defined as a function of the set, not of a sequence, as
follows.

1. **Contact curves.**
   - A blend's contact curve on face `F` is where its section touches `F`:
     - for a round, the rolling ball's contact, so the spine is where the
       two faces' offsets meet [2, 4];
     - for a chamfer, the line or curve at the given distance from the edge,
       measured on `F`.
   - On a plane beside a straight edge it is a line parallel to the edge.
     On curved faces it is the marched contact the crate already computes.
2. **Setback corners.**
   - Where two edges of `F` meet at a vertex `v`, the face's setback corner
     `S(F, v)` is where their contact curves on `F` meet.
   - An unblended edge contributes the edge itself, at distance 0.
3. **Strips** run between `S(F1, v)` and `S(F1, w)` on one side and
   `S(F2, v)` and `S(F2, w)` on the other, for an edge `vw` between faces
   `F1` and `F2`. Their ends are where the section sweep meets the vertex
   blend.
4. **Vertex blends** fill the polygon of setback corners around `v`:
   - a planar facet for chamfers at a 3-valent vertex;
   - a sphere where an equal-radius convex or concave trihedral corner
     allows one;
   - otherwise one n-sided G1 patch (§4.4).
5. **Admissibility.** A request fits when:
   - every face's setback curves leave it a simple region, so no setback
     curves cross and no shrunk face inverts;
   - every strip's contact curves stay on their faces;
   - every vertex polygon is non-degenerate.

   These predicates are the domain. `largest_fillet_radius` and
   `largest_chamfer_distance` become their exact boundary, not a separate
   search. Outside the domain the request is refused by name, naming the
   face or vertex that collapses.

Rounds have a global reference semantics: rounding every convex edge is the
morphological opening by a ball of radius `r`, and filleting every concave
edge is the closing [3]. This is the **oracle** for tests on pure-convex and
pure-concave regions. It is not the construction. It has no meaning for
mixed vertices or partial requests, and building it is far costlier than the
local model.

## 4. Architecture

### 4.1 P0: normalise and classify

- **Chains.** Group tangent-continuous edges into chains, so a chain is one
  strip with one section law.
- **Vertex classes.** Classify every vertex as `N × M` (N edges, M blended),
  with the count of convex, concave and tangent edges. This is the
  classification commercial systems publish for vertex blends [7, 14]. It
  selects the vertex method in P3 and nothing else; no geometry is decided
  by class.
- **Stable keys.** Give every element a key derived from the geometry and
  topology of the input, not from request order. All iteration below runs in
  key order, which makes the result bit-identical under any permutation of
  the request.

### 4.2 P1: sections and contact curves, per chain, on the original solid

- Each chain's section is solved against the **original** faces. Nothing
  built in P1 changes the solid another chain sees, so order independence
  holds by construction, not by sorting.
- **Exact where a closed form exists.** Plane–plane is a cylinder (round) or
  a plane (chamfer); plane–cylinder and plane–cone coaxial are a torus; and
  the other closed forms in `support.rs` (`planar_seat`, `revolved_seat`),
  `fillet.rs` (`section_at`) and `asym.rs` (`seated_conic`, `conic_pieces`).
- **Marched otherwise.** `march.rs` (`march_blend`, `march_blend_law`,
  `seat_section`) and `pair_marched.rs` (`pair_seat`) already take two
  surfaces, a section law and a guide, and return the spine and contacts in
  each support's own parameters. They need no topology.
- **Strips are extended.** Each strip is computed past its edge's ends far
  enough to reach its vertex blends, as the SMLib sheets are made "way long"
  [14]. Its true ends come from P2, not from the section solver.

### 4.3 P2: setbacks and ends, by intersection

- On each face, intersect neighbouring contact curves to get the setback
  corners `S(F, v)`. On planes this is a 3×3 linear solve. On curved faces
  it is a curve–curve intersection in the face's parameter space.
- A blend ends where its contact curve meets:
  - another blend's contact curve on the same face (a setback corner);
  - the face's boundary (an unblended neighbour edge);
  - a face it was not seated on (roll-over: the contact leaves its face).

  This is one intersection query per end, not a case. The kind of end only
  labels the result.
- **Roll-over** is first refused by name, as today, then handled in a later
  phase. The strip continues on the new face with the contact re-seated
  there, the cliff-edge and roll-over geometry of [14, 15].
- **Admissibility predicates (§3.5)** are evaluated here, on curves and
  corners, before any surface is built. A request that fails them is refused
  with the colliding pair named. Nothing is built and then measured.

### 4.4 P3: vertex blends, one method with stated exact cases

Every vertex blend is a hole bounded by the strips' end sections and the
faces' setback curves near `v`. Its boundary is closed and its corner data
consistent by construction, because each corner is a setback corner shared
by exactly the two curves that meet there.

1. **Exact cases, chosen by geometry, not by catalogue:**
   - **Planar facet:** a chamfer at a 3-valent vertex. The setback corners
     are three points, so the facet is a plane.
   - **Sphere:** equal-radius rounds at a trihedral vertex where one ball
     touches all three faces. That is today's ball-and-block corner, kept as
     a special value of the general method.
   - **Torus or plane continuation:** a tangent `2 × 2` vertex, where one
     strip continues the other.
2. **General case:** one n-sided G1 patch interpolating the strips' end
   sections and the faces along the setback curves.
   - The method is a ribbon-based transfinite patch [8, 9], or Gregory's
     n-sided patch [10, 11] for n ≤ 6.
   - It is evaluated, or fitted once at a fixed resolution with the residual
     reported. It is not refined until a check passes.
   - `fill_n.rs` (`make_filling_n_region`) stays the fitter where a B-spline
     result is required. It is fed consistent, twist-compatible corner data
     and a region known in advance.

   Setbacks are the knob Várady and Rockwood give for patch quality: larger
   setbacks give rounder patches and fewer bad corner angles [5, 6].
3. **Chamfers at n > 3 vertices:**
   - The setback corners need not be coplanar.
   - The facet is planar where they are.
   - Otherwise it is triangulated from the polygon's centroid, or refused by
     name. That choice is an open decision (§13).

### 4.5 P4: one assembly

- Every original face is trimmed once, along its setback curves.
- The strips and vertex blends are built as faces whose boundary edges are
  the setback curves and section curves computed in P2.
- The shell is assembled by **known adjacency** from the P0 pseudo-topology:
  each new edge is created once and shared by its two faces. No tolerance
  sewing searches for matches, and no boolean runs on the local path.
- The one boolean that remains is for **non-local** blends [12]: a strip
  that, extended, crosses a face that is not its neighbour (a fillet
  reaching a boss). The strip set is closed into a blend solid and one
  boolean is applied (the "sheets then one boolean" design, §9 B). The
  first steps refuse it by name, as non-local; it is a later phase.

### 4.6 P5: independent guards

These are kept, and kept independent of the construction:

- `check`, and orientation consistency of every new face.
- **Region:**
  - every convex strip and vertex blend lies inside the original solid;
  - every concave one lies inside the solid's closing at the largest radius
    asked, within the fill tolerance;
  - every vertex blend lies inside the convex hull of its boundary
    ribbons' region.

  These bounds come from the semantics in §3, not from per-case cells.
- **No overhang:** no new face looks down past its strips' contact normals
  by more than the fill wobble.
- `analyse_blend` on every strip and vertex blend: gap, tangency within
  2e-4 rad, curvature where G2 is asked.

A guard that fails is an engine bug and refuses the request. It is never the
way the engine searches for a working construction.

## 5. Polyhedral fast path: exact chamfers (prototype result)

A throwaway prototype (Python, not committed) built the §3 structure for
planar solids. Per face and vertex it solves for `S(F, v)` (three planes),
then builds:

- shrunk faces;
- edge quads, which are planar because the two setback lines are parallel to
  the edge;
- vertex polygons.

It checked the result: closed (volume independent of the origin), no two
triangles crossing, and exact volumes. On a stepped bracket like #902's
(profile 60 deep, plate 6 thick, step to 16, block to 30, 80 long; 16
vertices, 24 edges, all 3-valent):

| request | result | time (Python) |
|---|---|---|
| every edge, d = 1, 2, 3 | valid, closed, 0 crossings; every edge bevelled, side profiles included; a triangular facet at every vertex, the mixed ones included | ≈ 4 ms |
| every edge, d = 4.9 | 13 crossings: past the domain, since the 6-thick plate admits d ≤ 3, and detected | ≈ 4 ms |
| the concave step edge alone, d = 2 | ends flush with both side faces; volume 69,760 = 69,600 + d²/2 × 80, exactly | ≈ 3 ms |
| asymmetric distances 1 / 2.5 | valid, closed | ≈ 4 ms |
| five random halves of the edges | all valid and closed | ≈ 3.5 ms each |

- **Order independence.** The construction reads a function of
  `(edge, face)`, never a sequence, so every permutation of the request
  gives the same solid.
- **No booleans and no fitting.** Every coordinate is a 3×3 solve.
- **The domain** is the face-polygon offset staying simple plus the facets
  staying non-degenerate. Both are exact predicates, and the largest distance
  is computable in closed form per face.

Planar-faced requests whose sections are all chamfers take this path
outright. That covers the #902 bracket's chamfer-all-edges case and every
chamfer of a prismatic part.

**Rounds on polyhedral parts** use the same backbone with exact geometry:

- **Strips:** cylinders, between tangent lines at `r · cot(α/2)` from the
  edge for an interior dihedral angle `α`.
- **Vertices:** spheres where one ball fits; otherwise the §4.4 patch.

Their oracle is the closed-form section area along each strip, plus the
opening and closing reference [3] on pure-convex or pure-concave parts.

## 6. Soundness by construction

| Failure seen | How the design rules it out |
|---|---|
| Order and radius dependence (#838, #809, #878) | P1 reads only the original solid. Every phase iterates in stable-key order. No phase depends on another chain's output except through P2's symmetric intersections. |
| Blend past a face's side over air (#868, #871) | A contact curve that leaves its face is an end found in P2. The strip is trimmed there, or roll-over is refused, before any surface exists. There is no "run off and cut later". |
| Corner patch outside the solid (#713, #839, #846) | The vertex hole's boundary and region are fixed before the patch. Exact cases need no fit. The general patch is evaluated or fitted once against known corner data, and P5 measures it independently. |
| Feature lost inside a corner cell (#843) | No cells are cut: faces are trimmed along curves, so nothing inside a box is discarded. |
| Concave bevel past both side faces, unbevelled profile edges, inconsistent facets (#902) | §5: the face-setback structure gives every edge its quad and every vertex its facet by one rule. |
| Fill nets refined to 202 controls, then checked (#870) | Exact cases need no fitting. The general patch has one fixed resolution with its residual reported. |

## 7. Complexity targets

- **P0–P2:** `O(E + V)` local work, plus one march per curved chain. The
  march's cost is the solver's, unchanged.
- **P3:** constant work per exact vertex. One patch per general vertex, at a
  fixed resolution.
- **P4:** one trim per face and one assembly: `O(F + E + V)`.

Targets, release build:

| Case | Target |
|---|---|
| Chamfer every edge of a 50-edge prismatic part | under 10 ms |
| Round every edge of a 50-edge prismatic part | under 1 s, as #902 asks; analytic strips and sphere or single-patch corners make well under 100 ms reachable |
| Curved chains | bounded by the march; no retry ladders and no whole-request reruns |

## 8. How #877 and variable radii fit

- **The section is the only per-edge law.**
  - #877's asymmetric distances and conic profiles (`asym.rs`,
    `asym_law.rs`, merged through #898) become section laws in P1. A conic
    section's contact curves sit at `d1` and `d2` on the two faces.
  - Asymmetric chamfers are a chamfer section with two distances. The
    prototype already builds them on the polyhedral path.
  - Variable radius is a section whose parameters follow `RadiusLaw`
    (`law.rs`) along the chain.
- **P2–P5 never look inside a section.** Setback corners are where contact
  curves meet, whatever made the curves. The vertex patch interpolates
  whatever ribbons the strips present.
- **One open semantic question:** what a chamfer meeting a round at one
  vertex should look like (§13).
- **Coordination:** the #729 agent's #877 work keeps landing on the current
  path. Its section math (`seated_conic`, `conic_pieces`, the law) moves
  into P1 unchanged. Its end and corner handling (`after_symmetric`, corner
  routing in `blend.rs`) is what retires.

## 9. Alternatives weighed

- **A. Keep extending the catalogue.**
  - Rejected. Each case adds constructions, refusals and checks. The wrong
    solids came from cases misapplied or checks incomplete, and order
    dependence is inherent to sequential application.
- **B. Blend sheets, then one boolean or merge** (the SMLib executive's
  global merge [14]).
  - Strong on generality: non-local blends fall out.
  - Booleans on blend sheets meet near-tangent faces at every contact curve,
    which is the hardest input a boolean has. The result's topology is
    discovered, not known.
  - Kept as the P4 fallback for non-local blends only.
- **C. Morphological offsets** (opening and closing [3]).
  - Exact rolling-ball semantics for whole regions.
  - Undefined for mixed vertices, partial requests, chamfers and asymmetric
    sections. Needs exact solid offsetting, which is costlier than the
    blends themselves.
  - Kept as the test oracle.
- **D. Implicit or potential blends** [16].
  - Smooth vertex blends for free, but they don't produce the trimmed
    parametric faces a b-rep needs. Rejected.
- **E. Direct topology surgery per edge** (Euler operators, edge by edge).
  - Same fragility as A: each step changes what the next sees.
  - The adopted design is surgery done once, from the precomputed
    pseudo-topology.

## 10. What is reused and what retires

Classes from a survey of `main`.

**Reused unchanged or nearly:**

| File or function | Role |
|---|---|
| `march.rs` (`march_blend`, `march_blend_law`, `seat_section`) | Marched sections; topology-free |
| `pair_marched.rs` (`pair_seat`) | Guide where the surfaces or their offsets cross |
| `support.rs` (`planar_seat`, `revolved_seat`), `fillet.rs` (`section_at`), `ruled.rs`, `asym.rs` (`seated_conic`, `conic_pieces`), `law.rs` | Closed-form sections and section laws |
| `marched.rs` (`conic_section`, `fit_open_band`, `fit_closed_band`) | Band fitting. Not the wedge or run-out builders. |
| `fill_n.rs` (`make_filling_n_region`) | Fitter for general vertex patches, fed consistent corner data |
| `fill_patch.rs`, `bridge.rs` | Hermite blend surfaces |
| `between.rs`, `support.rs` (`ball_fits_by`, `bands_clear`) | Admissibility predicates, read on the original solid |
| `fits.rs` | Domain boundary, reworked onto the §3.5 predicates |
| `analyse.rs` | P5 measurement |
| `setback.rs` (`certify_within`, `least_setbacks`) | Closed-form corner maths |
| `corner_curved.rs`, `turned_pole.rs`, `chamfer_marched.rs` | Section and chart utilities |

**Retired in steps:**

| File | Case |
|---|---|
| `setback.rs` (corner constructions), `cusp.rs`, `sharp.rs`, `corner.rs` (the corner booleans, not its frame helpers), `chamfer_corner.rs` | Corner catalogue |
| `junction.rs`, `carried.rs`, `rolled.rs`, `stepped.rs`, `consumed*.rs`, `pinched.rs` | End and junction catalogue |
| `fillet.rs` (`cut_off_bands`, the round/retry orchestration), `chamfer.rs` (sequential wedges), `blend.rs` | Sequential orchestration |

What these do correctly today becomes corpus cases for the new path.

## 11. Migration plan

Every step keeps the current public API. A step routes a class of request
to the new engine only when that class's corpus passes there. Everything
else stays on the current path.

- **M1: polyhedral chamfers.**
  - New engine for requests whose faces are all planar and whose sections
    are all chamfers, symmetric or asymmetric.
  - Oracles: exact volume (divergence over the built faces, against an
    independent plane-arrangement volume), closedness and no crossing, and
    every permutation of the request bit-identical.
  - Ships with the #902 bracket chamfer case.
  - The old chamfer path stays for curved faces.
- **M2: request model and P0–P2 for planar faces with rounds.** Exact
  cylinder strips. Sphere and planar exact vertices. The general vertex
  patch through `make_filling_n_region` with fixed resolution.
  - Corpus: box, L, T, U, stepped block, the #902 bracket fillet-all-edges,
    mixed rounds and chamfers.
- **M3: curved faces.** P1 on `march.rs` and `pair_marched.rs`. P2
  intersections in face parameter space. The plate with a boss, tangent and
  not (the #868 and #878 families).
- **M4: the general vertex patch** (ribbon-based transfinite or Gregory),
  replacing iterative fill refinement for vertex holes.
- **M5: roll-over and non-local blends.** Roll-over re-seats the contact on
  the new face. Non-local blends go through the one-boolean fallback.
- **M6: retire old paths,** file by file, once no routed class reaches them.
  The real-world check (Voron 2.4, DoomCube) runs at every routing change.

**"Keeps every currently correct result"** means that every case the current
tests accept passes on the new path within its oracle. Bit-identity to the
old path is not promised where the old construction was a fitted patch,
since the new one is exact or fitted differently. It is promised for
closed-form strips (cylinders, tori), which both paths build exactly. A case
the old path refused may now build. A case the old path built and the new
one refuses is a migration blocker.

## 12. Acceptance corpus

From #902:

- The stepped bracket, chamfer and fillet on every edge, matching the
  reference render.
- Every edge of a box, an L-bracket, a T and a U channel, a stepped block,
  and a plate with a boss (tangent and not), at several radii and distances.
- Mixed requests: rounds with bevels, asymmetric distances (#877), variable
  radii.

On every case:

- every permutation of the request gives the bit-identical solid;
- no refusal on planar polyhedral parts inside the §3.5 domain;
- exact oracles: polyhedral volume for planar chamfers, closed-form section
  areas for prismatic rounds, the morphological reference on pure-convex and
  pure-concave regions; plus the region and no-overhang guards;
- the §7 speed targets;
- the real-world check unchanged or better.

## 13. Open decisions

1. **Chamfer facets at vertices of more than three edges,** where the
   setback corners are not coplanar: one triangulated facet from the
   centroid, several facets, or a refusal by name.
2. **Default vertex blend for rounds:** the sphere wherever one ball fits,
   with the general patch elsewhere, or the setback patch everywhere for
   uniform shape. Also the default setback, if any: Várady–Rockwood
   setbacks trade corner size for patch quality.
3. **A chamfer and a round meeting at one vertex:** a planar facet, the
   general patch, or extending the round's strip into the chamfer plane.
4. **Non-local blends:** refuse by name until M5, or prioritise the
   one-boolean fallback.
5. **`*_partial` semantics** (greedy in request order today) under a
   set-valued request: drop them, or redefine them as the largest admissible
   subset by stable key.
6. **Special-case work on the current path until M2:** frozen (#900 parked,
   #877 continuing), or allowed for regressions only.

## References

1. J. Vida, R. R. Martin, T. Várady. "A survey of blending methods that use
   parametric surfaces." *Computer-Aided Design* 26(5):341–365, 1994.
   doi:10.1016/0010-4485(94)90023-X.
2. B. K. Choi, S. Y. Ju. "Constant-radius blending in surface modelling."
   *Computer-Aided Design* 21(4):213–220, 1989.
   doi:10.1016/0010-4485(89)90046-8.
3. J. R. Rossignac, A. A. G. Requicha. "Offsetting operations in solid
   modelling." *Computer Aided Geometric Design* 3:129–148, 1986. Rounding
   and filleting as compositions of offsets.
4. S. Y. Ju. *Construction of blend surfaces in surface modelling.*
   Doctoral thesis, KAIST, 1989. Constant- and variable-radius edge blends.
5. T. Várady, A. Rockwood. "Vertex blending based on the setback split." In
   *Mathematical Methods for Curves and Surfaces*, Vanderbilt University
   Press, 1995, pp. 1–16.
6. T. Várady, A. Rockwood. "Geometric construction for setback vertex
   blending." *Computer-Aided Design* 29(6):413–425, 1997.
7. T. Várady, C. M. Hoffmann. "Vertex blending: problems and solutions." In
   *Mathematical Methods for Curves and Surfaces II*, 1998, pp. 501–527.
8. P. Salvi, T. Várady, A. Rockwood. "Ribbon-based transfinite surfaces."
   *Computer Aided Geometric Design*, 2014.
9. P. Salvi, T. Várady. "A multi-sided generalization of the C⁰ Coons
   patch." arXiv:2002.11347, 2020. T. Várady, P. Salvi, Gy. Karikó. "A
   multi-sided Bézier patch with a simple control structure." *Computer
   Graphics Forum*, 2016, doi:10.1111/cgf.12833.
10. J. A. Gregory. "N-sided surface patches." In *The Mathematics of
    Surfaces*, Oxford University Press, 1986.
11. P. Charrot, J. A. Gregory. "A pentagonal surface patch for computer
    aided geometric design." *Computer Aided Geometric Design* 1(1), 1984.
12. I. C. Braid. "Non-local blending of boundary models." *Computer-Aided
    Design* 29(2):89–100, 1997.
13. I. Stroud. *Boundary Representation Modelling Techniques.* Springer,
    2006. Blending in b-rep modellers.
14. SMLib fillet design, read from the Apache-2.0 sources in `usd-brep`
    (`source/SMLib`: `SmFilletExecutive`, `SmFilletCorner`,
    `SmFilletSolver`, `SmFilletIntersector`) and their Rust port in
    `solidcore` (`crates/solidcore/src/blend/`). Concepts taken: a whole-model
    pseudo-topology built before geometry; phased corner and edge passes;
    the offset-surface intersection as the section solver; pluggable
    cross-section generators; N × M vertex classes with a bevel-or-blend
    choice; "way long" sheets trimmed by intersection; one merge at the end.
15. A. P. Rockwood, J. C. Owen. "Blending surfaces in solid modeling." In
    *Geometric Modeling: Algorithms and New Trends*, SIAM, 1987.
16. C. M. Hoffmann, J. Hopcroft. "The potential method for blending
    surfaces and corners." In *Geometric Modeling: Algorithms and New
    Trends*, SIAM, 1987.

No implementation source of any other kernel was read for this record.
