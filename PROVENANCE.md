# Provenance

This file records where ogeom's code and design came from, so that users and
contributors can judge its licence for themselves. ogeom is LGPL-2.1-only,
and it is a derivative work of Open CASCADE Technology (OCCT) in the
licensing sense. It contains no code from OCCT's implementation, translated
or otherwise; the audit below is the evidence. What makes it derivative is
the reference-driven development recorded here.

`CONTRIBUTING.md` states the naming rule that goes with this: the licence
and the provenance name OCCT plainly, and the rest of the repository
writes in the field's vocabulary.

## The licence

ogeom is licensed under the GNU Lesser General Public License, version 2.1,
and no other version: the same licence OCCT itself carries. The reason is
the relationship this file records. The public API vocabulary, the
data-model conventions and the default tolerances follow OCCT; the
project's first goal was a kernel that could take OCCT's place in
OCCT-based software; and for one day its implementation source was read.
Reference-driven development makes the kernel a derivative work in the
licensing sense even where the audit below finds no translated code, and
no audit can bound what the AI assistants that wrote most of this code
learned from before this project began. LGPL-2.1 states that position
honestly, and it keeps the work usable: a library under it may be linked
from software under any licence, which is what an application kernel has
to allow.

Copies of releases published before the change were distributed under
MIT OR Apache-2.0 and remain under those terms; the licence above applies
to releases published after it.

## The rule, from the first commit

The first commit (31 July 2026) included `CONTRIBUTING.md` with this rule:
"Do not copy, translate, or transliterate another kernel's source into
`crates/`." The rule has only grown stricter since. Today it also forbids
reading another kernel's implementation source for work here, and giving it
to any tool, an AI assistant included, as context.

The implementation sources are:

- published algorithms and papers: Piegl & Tiller for NURBS, Shewchuk for
  robust predicates, Brent's method, Yuksel's polynomial roots, and the
  marching and surface/surface intersection literature;
- format standards: ISO 10303 for STEP, the IGES specification, and the
  published descriptions of the other formats;
- first principles, checked by tests against results computed independently:
  closed forms, round trips, invariants and properties over random inputs.

## Open CASCADE Technology (OCCT)

OCCT, which is LGPL-2.1, is the kernel ogeom was most often compared with.
The relationship is set out here in full.

**Why OCCT was a reference.** The project's first goal was a kernel that could
be plugged into FreeCAD, which is built on OCCT. To learn what that would
take:

- FreeCAD's use of OCCT's public API was profiled, from public headers only;
- OCCT's module and toolkit manifests and its public headers defined a
  feature-parity scope;
- capabilities were audited against that scope.

That goal has been dropped. ogeom is not meant to replace OCCT in FreeCAD or
in any other program, and the tooling and data built for it have been removed
from the repository. The scope is now stated in ogeom's own terms
(`docs/SCOPE.md`), and capabilities are recorded against it
(`docs/CAPABILITIES.md`).

**Public API names.** Several public names follow OCCT's vocabulary
(`make_prism`, `make_pipe_shell`, `fuse` / `cut` / `common`, `sew`,
`fix_shape`, `unify_same_domain`, `is_same` / `is_partner` / `is_equal`). They
were chosen so that code written against OCCT would map onto ogeom easily.
They are names only: the implementations behind them are ogeom's own. Names
that exist only to mirror OCCT may change in a future release.

**Data-model conventions.** `ogeom-topo` adopts B-rep conventions that OCCT
users know: a shape as a shared node plus a placement and an orientation, the
eight shape types, the four orientations and how they compose, three levels of
shape identity, placement as a chain of transforms, and several curve
representations on one edge. They were adopted so that models and `.brep`
files written by OCCT-based software map onto ogeom without loss. The code
that expresses them is independent: arenas with stable keys rather than
reference-counted handles, a placement chain built by appending rather than a
linked list, and an entity identity and provenance model
(`docs/DATA_MODEL.md` §8) that departs from OCCT's pointer identity on
purpose.

**Default tolerances.** The defaults in `ogeom-core` equal OCCT's `Precision`
values: confusion `1e-7`, angular `1e-12`, and intersection, parametric
confusion and approximation at the same ratios to confusion. The reason is
interoperability. STEP and `.brep` files written by OCCT-based software state
their vertex and edge tolerances on that basis, and a kernel that reads them
has to agree on when two points are the same, or valid models arrive invalid.
Nothing else in the tolerance code corresponds to OCCT's: ogeom scales linear
tolerances by the model's unit, adds a separate degeneracy threshold, and
keeps per-entity tolerances with containment rules.

**Vocabulary.** The boolean code uses the words "pave", "interference" and
"same domain". These words are also used in OCCT's public Boolean Operations
user guide. In ogeom they name ogeom's own structures: a `Pave` is a parameter
with a tolerance and spread model, there are no interference tables, and the
pipeline is arranged differently (face pairs found through a box tree and
intersected in parallel, then each face arranged in its own parameter space).

**AI-assisted development.** Most of ogeom was written with AI coding
assistants, under the same rule. OCCT's public headers and manifests were used
as input for the scope and API-surface work described above. Implementation
source was not supplied to these tools, with one exception. On 21 October
2026, agents researching two boolean robustness issues read parts of OCCT's
boolean implementation. One fix written afterwards (PR #388) was compared line
by line with the OCCT functions consulted; its methods differ from them. That
reading broke the spirit of this project's rule, and the rule now forbids it
explicitly. No one can say what public code the underlying AI models were
trained on. That is why the audit below examines the code itself.

## The audit

In October 2026, about two dozen of the most algorithm-specific functions and
structures were compared side by side with OCCT's implementation, across
mathematics, topology, intersection, booleans and tolerances. None is a
translation or a close paraphrase. No OCCT identifiers, naming conventions or
comments appear in the code, and no module follows OCCT's file layout or call
graph.

| Area | Finding |
|---|---|
| Mathematics | Published algorithms, written independently. Brent in inverse-quadratic form, Yuksel polynomial roots, de Boor evaluation and knot insertion from Piegl & Tiller, degree elevation by Bézier decomposition. |
| Topology | Shared design conventions (above), independent code. |
| Intersection | Independent. One walker shared by surface/surface intersection, silhouettes and blend marching, with a single spatial step and a sag bound; its own seeding and closed forms. |
| Booleans | Independent decomposition. No counterpart to OCCT's data structure, pave blocks, common blocks or interference tables. |
| Tolerances | Same default values, for interoperability (above). Nothing else shared. |

## Reporting a concern

If a function or file looks like a translation of another kernel's code,
open an issue naming it. It will be compared in public and, if it is
derived, either rewritten or recorded here with the rest of the
provenance: under the LGPL, a derivation in the open is compliance.
