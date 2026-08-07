# Scope

What belongs in this kernel, what does not, and how to tell without arguing.

## The rule

**ogeom targets parity with the reference kernel's modelling modules, and nothing
else.** Four modules are in scope:

| Module | What it is |
|---|---|
| **FoundationClasses** | Arithmetic, primitives, solvers, tolerances, errors. |
| **ModelingData** | The geometry and topology vocabularies — curves, surfaces, the b-rep data model. |
| **ModelingAlgorithms** | Intersection, booleans, blending, offsets, sweeps, healing, tessellation, hidden-line removal. |
| **DataExchange** | STEP, IGES, STL, VRML, OBJ, glTF, PLY, and the document structure they carry. |

Three are out, permanently:

| Module | Why |
|---|---|
| **Visualization** | Rendering, viewers, interactive selection. A kernel is not a renderer. |
| **ApplicationFramework** | The generic label-and-attribute document tree. The *exchange* document is in scope — it is part of DataExchange — but the framework beneath it is an application's concern. |
| **Draw** | A test harness with its own scripting language. |

Anything the reference does not do at all is out of scope by default. Reverse
engineering from meshes, constraint solving, feature recognition and process
planning are all real disciplines and none of them is this.

Each of those was built here before this rule was written, and each of them
works. Rather than delete working code to make a point, it lives in `outside/`,
which is a separate workspace the kernel's `Cargo.toml` excludes by name. The
exclusion is what makes the rule structural instead of aspirational: nothing
there can be pulled back in by a path dependency without someone deleting that
line on purpose. `outside/README.md` says why each crate is on the far side.

## How to decide a case

The question "is this in scope?" is answered mechanically, not by taste:

1. `adm/MODULES` in the reference tree maps each module to its toolkits.
2. `src/<Toolkit>/PACKAGES` maps each toolkit to its packages.
3. `src/<Package>/*.hxx` are that package's classes.

The union over the four in-scope modules is **276 packages / 6,267 public
headers**. If a capability's counterpart is in that set, it is in scope. If it is
in Visualization, ApplicationFramework or Draw, it is not. If it has no
counterpart at all, it is not.

`docs/parity/reference-index.tsv` is that set, committed, so the question can be
answered without a reference checkout. `docs/PARITY.md` records where we stand
against it.

## What the scope rule is *not*

**Not a licence to mirror.** `CONTRIBUTING.md` forbids reproducing another
kernel's class hierarchy, decomposition and file layout, and that still holds.
Parity is a claim about *capability*, not about structure: the parity record is
keyed on what a caller would ask for, and each entry names the reference
packages it accounts for. A capability we deliberately provide differently is
recorded as `divergent`, with the reasoning — not as a gap.

**Not driven by usage data.** `docs/api_surface.json` profiles how one large
application exercises the reference. It is a **sequencing** input — it says what
to get right first — and it appears in the parity index as a column for exactly
that purpose. It has never been a scope input and is not one now. Its own
generator says so: *"What it is emphatically not good for: deciding what to
build."* A capability inside the four modules is in scope whether or not that
application ever calls it.

**Not a size target.** 6,267 headers is not 6,267 things to build. Most of that
count is generic instantiation — `TColStd_Array1OfReal` and its several hundred
siblings — which Rust's generics give for free. The triage rules in
`tools/apisurf/apisurf.py` are what reduce the number to the capabilities
underneath it, and every rule is recorded with the headers it removed so the
reduction is auditable rather than asserted.

## Where the scope changes

Here, by editing this file, with the reasoning written down. Not in a pull
request that quietly adds a crate.
