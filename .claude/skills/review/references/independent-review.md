# Independent review panel

The lead scopes the review, picks the lenses, checks each claim and issues
the verdict. Lenses are evidence questions, not personas.

## Budget

- At most three lenses per round.
- At most two rounds for one task and one scope: one initial round and,
  only when corrections or unresolved evidence need it, one selective
  follow-up. Never a third round. After the budget, the lead checks directly
  and carries open evidence into the verdict.
- Pick only lenses whose question can change the verdict. Never launch a
  lens to fill a quota, and never defer an obviously required lens to round
  two.
- A follow-up reruns no lens only because it ran before. Pick the smallest
  set from the correction diff and the open findings, or none.

## Lenses

| Lens | Question |
|---|---|
| Facts | What changed, which operations and paths does it touch, and what evidence is missing? |
| Caution | How can this return a wrong shape that passes `check`, lose a void or a face, or fail at a degenerate placement? |
| Simplicity | Is this the smallest sufficient diff and the simplest correct algorithm, with superseded code gone? |

| Specialist | Trigger | Focus |
|---|---|---|
| Numerics | tolerances, predicates, solvers, marching, fitting | scale, conditioning, convergence, the stated error of an approximation |
| Topology | `ogeom-topo`, result assembly, `check`, healing | the invariants of `docs/DATA_MODEL.md`: orientation composition, location chains, identity levels, tolerance containment, history |
| Exchange | `ogeom-io`, `ogeom-doc` | round trips, both writings of a construct, units, voids, assemblies, colours and names |
| Tests and oracles | changed tests, or untested material behaviour | an oracle independent of the code, a property over examples, a refusal pinned by text (see test-audit) |
| Performance | hot paths, STEP reads, booleans, meshing | `ogeom-bench --check` ratios, amplification, allocation |

## Packet

Every reviewer gets the same frozen packet: the task and acceptance
criteria, base and head SHAs, changed and excluded scope, non-goals,
`CONTRIBUTING.md`, the risk tier, the allowed read-only commands, exactly
one lens, and the result schema below. Never include provisional or sibling
findings.

Reviewers may read, search and run non-mutating checks (`cargo test -p
<crate> <filter>`, a replay with `ogeom-stress --case`, a repro program
outside the repository). They may not edit tracked files, commit, push,
post, or launch nested reviewers. Retry a failed lens once, in the same
round, only when a concrete cause changed.

## Result schema

Each reviewer returns: coverage (files and paths read), candidate findings
(`file:line`, rule slug, failure scenario, change-causal evidence, smallest
fix), hypotheses it rejected, and open questions. "No findings" is valid.
The lead verifies every candidate against the code before it enters the
report.
