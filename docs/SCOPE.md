# Scope

This file says what belongs in the kernel, what does not, and how to decide a
case.

## The rule

ogeom targets parity with the reference kernel's modelling modules, and nothing
else. There is one deliberate addition: turning meshes into solids (see
[below](#the-one-addition-meshes-into-solids)).

Four modules are in scope:

| Module | What it covers |
|---|---|
| **FoundationClasses** | Arithmetic, primitives, solvers, tolerances, errors. |
| **ModelingData** | The geometry and topology vocabularies: curves, surfaces, the b-rep data model. |
| **ModelingAlgorithms** | Intersection, booleans, blending, offsets, sweeps, healing, tessellation, hidden-line removal. |
| **DataExchange** | STEP, IGES, STL, VRML, OBJ, glTF, PLY, and the document structure these formats carry. |

Three modules are permanently out of scope:

| Module | Why it is out |
|---|---|
| **Visualization** | Rendering, viewers, interactive selection. A kernel is not a renderer. |
| **ApplicationFramework** | The generic label-and-attribute document tree. The *exchange* document is in scope, because it is part of DataExchange. The framework beneath it belongs to the application. |
| **Draw** | A test harness with its own scripting language. |

Anything the reference kernel does not do at all is out of scope by default.
Examples: constraint solving, feature recognition and process planning. These
are real disciplines, but they are not this kernel.

## The one addition: meshes into solids

ogeom turns a mesh into a solid and recognizes its surfaces
(`solid_from_mesh`, `recognize_points`).

What the reference kernel does: its modelling algorithms build a shape on a
mesh (one planar face per triangle, with shared edges) and merge coplanar faces.

What ogeom adds: it finds which regions of triangles lie on a cylinder, a
cone, a sphere or a torus, and rebuilds those regions on those surfaces.

Why: the same reason the exchange module exists. Printers, slicers and model
sites exchange meshes. A kernel that reads STL, OBJ and 3MF but can only display
the result leaves the application to rebuild the geometry itself.

Recognition meets the kernel's standard, not a heuristic's:

- every surface is verified against every sample, at a stated tolerance;
- every edge between recognized faces is placed exactly on both surfaces;
- a region that cannot be built this way stays faceted, and is counted.

This is what makes recognition a construction the kernel can stand behind.

Still out of scope: fitting free-form surfaces to scans, and reading design
intent back out of topology.

## Code that is out of scope: `outside/`

Some out-of-scope disciplines were built here before this rule existed, and
they work. Instead of deleting working code, it lives in `outside/`.

- `outside/` is a separate workspace.
- The kernel's `Cargo.toml` excludes it by name.
- So no path dependency can pull it back in unless someone deliberately deletes
  that exclusion. This makes the rule structural, not just a statement of
  intent.

`outside/README.md` explains why each crate there is out of scope.

## How to decide a case

The answer comes from the reference tree's own files, not from opinion:

1. `adm/MODULES` maps each module to its toolkits.
2. `src/<Toolkit>/PACKAGES` maps each toolkit to its packages.
3. `src/<Package>/*.hxx` are the package's classes.

Across the four in-scope modules this gives **276 packages and 6,267 public
headers**. Then:

- If a capability's counterpart is in that set, it is in scope.
- If it is in Visualization, ApplicationFramework or Draw, it is out.
- If it has no counterpart at all, it is out.

`docs/parity/reference-index.tsv` is that set, committed, so you can answer the
question without a reference checkout. `docs/PARITY.md` records where ogeom
stands against it.

## What the scope rule does not mean

**It does not allow mirroring.** `CONTRIBUTING.md` forbids copying another
kernel's class hierarchy, decomposition or file layout, and that still holds.
Parity is about *capability*, not structure:

- the parity record is keyed on what a caller would ask for;
- each entry names the reference packages it accounts for;
- a capability we deliberately provide differently is recorded as `divergent`,
  with the reasoning. It is not a gap.

**It is not driven by usage data.** `docs/api_surface.json` profiles how one
large application uses the reference kernel.

- It is a **sequencing** input: it says what to get right first.
- It appears in the parity index as a column for that purpose only.
- It has never been a scope input. Its own generator says so: *"What it is
  emphatically not good for: deciding what to build."*
- A capability inside the four modules is in scope whether or not that
  application ever calls it.

**It is not a size target.** 6,267 headers are not 6,267 things to build. Most
are generic instantiations (`TColStd_Array1OfReal` and several hundred similar
headers), which Rust generics give for free. The triage rules in
`tools/apisurf/apisurf.py` reduce the count to the underlying capabilities.
Each rule is recorded with the headers it removed, so the reduction can be
audited.

## How the scope changes

Only by editing this file, with the reasoning written down. Never through a
pull request that quietly adds a crate.
