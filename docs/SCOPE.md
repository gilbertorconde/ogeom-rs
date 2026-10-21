# Scope

This file says what belongs in the kernel, what does not, and how to decide a
case.

## The rule

ogeom is a complete solid-modelling kernel with its data exchange, and nothing
else. One area goes beyond what a modelling kernel usually offers: turning
meshes into solids (see [below](#the-one-addition-meshes-into-solids)).

Four areas are in scope:

| Area | What it covers |
|---|---|
| **Foundations** | Arithmetic, primitives, solvers, tolerances, errors. |
| **Geometry and topology** | Curves, surfaces, and the b-rep data model. |
| **Modelling algorithms** | Intersection, booleans, blending, offsets, sweeps, healing, tessellation, hidden-line removal. |
| **Data exchange** | STEP, IGES, STL, VRML, OBJ, glTF, PLY, 3MF, DXF, and the document structure these formats carry. |

Three areas are permanently out of scope:

| Area | Why it is out |
|---|---|
| **Visualization** | Rendering, viewers, interactive selection. A kernel is not a renderer. |
| **Application frameworks** | A generic label-and-attribute document tree. The *exchange* document is in scope, because the formats carry it. The framework beneath it belongs to the application. |
| **Scripting harnesses** | Test shells with their own scripting languages. |

Disciplines that are not solid modelling are out of scope. Examples:
constraint solving, feature recognition and process planning. These are real
disciplines, but they are not this kernel.

## The one addition: meshes into solids

ogeom turns a mesh into a solid and recognizes its surfaces
(`solid_from_mesh`, `recognize_points`).

What a modelling kernel usually does: build a shape on a mesh (one planar
face per triangle, with shared edges) and merge coplanar faces.

What ogeom adds: it finds which regions of triangles lie on a cylinder, a
cone, a sphere or a torus, and rebuilds those regions on those surfaces. A
smooth region on none of them is rebuilt on a surface of revolution, an
extrusion, or a fitted B-spline patch, in that order of preference, before it
is left faceted.

Why: the same reason the exchange module exists. Printers, slicers and model
sites exchange meshes. A kernel that reads STL, OBJ and 3MF but can only display
the result leaves the application to rebuild the geometry itself.

Recognition meets the kernel's standard, not a heuristic's:

- every surface is verified against every sample, at a stated tolerance;
- every edge between recognized faces is placed exactly on both surfaces;
- a region that cannot be built this way stays faceted, and is counted;
- every step is a computation the kernel can repeat and measure. No learned
  model proposes a surface or a topology.

This is what makes recognition a construction the kernel can stand behind.

Neighbouring surfaces may be built from each other where a geometric
relation between them holds within the tolerance: a constant-radius fillet
is built from its two supports and a radius, so it meets them tangentially
by construction, and a corner is solved onto the surfaces that meet there.
The relation only places the surfaces. It is not kept, and the result is a
plain solid with no feature tree. A helical region (a thread) is told apart
from a surface of revolution so it can go to the fitted patch; no helical
surface type is added for it.

A fitted patch is held to the same standard. Facets are exact to the mesh
but nothing downstream can work with them: an offset, a shell or a fillet
across a faceted region gives facets back, or fails. A patch verified within
the tolerance (at the triangles' interiors as well as their corners) is the
better result, and the facets remain the fallback when no patch verifies.
The fitting itself is already in scope (`algo.fitting`); this only says
where the converter uses it.

The kernel also exposes the steps under the converter: the regions it found,
merging and splitting them, fitting a chosen surface type to one, and
building the solid from the result. An application can then correct an
automatic conversion. The selection and the interface for it belong to the
application.

Still out of scope: scans. Noisy or incomplete point data needs denoising,
hole filling and registration, which are a discipline of their own. Reading
design intent back out of topology is out of scope too.

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

- If a capability belongs to one of the four in-scope areas, it is in scope.
- If it belongs to visualization, an application framework or a scripting
  harness, it is out.
- If it is not solid modelling or data exchange at all, it is out.

`docs/capabilities.toml` lists the kernel's capabilities, one per thing a
caller would ask for, each with a verdict and the evidence for it.
`docs/CAPABILITIES.md` is its rendered form.

## What the scope rule does not mean

**It does not mean copying how another kernel does it.** `CONTRIBUTING.md`
forbids taking another kernel's code or reading its implementation source, and
forbids mirroring its class hierarchy, decomposition or file layout. Scope is
about *capability*, not structure:

- the ledger is keyed on what a caller would ask for;
- a capability we deliberately provide differently is recorded as `divergent`,
  with the reasoning. It is not a gap.

**It is not driven by any one application.** A capability inside the four
areas is in scope whether or not a given application calls it.

## How the scope changes

Only by editing this file, with the reasoning written down. Never through a
pull request that quietly adds a crate.
