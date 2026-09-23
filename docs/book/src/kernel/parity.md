# The parity audit

Most projects describe how complete they are with adjectives. ogeom measures
it with an audit that the build enforces. This chapter explains how the audit
works, so that [the ledger](parity-ledger.md) can be read as a measurement.

## What is measured

The reference target is fixed by [Scope](scope.md):

- four modelling modules;
- 276 packages;
- 6,267 public headers.

This set is committed as `docs/parity/reference-index.tsv`, so no reference
checkout is needed to use it.

The audit is keyed on **capabilities, not headers**:

- `docs/parity/parity.toml` lists 98 capabilities.
- Each capability *claims* the reference headers it accounts for.
- The gate requires the claims to be **total and disjoint**: every kept header
  is claimed exactly once, or the build fails.

This turns "nothing was forgotten" into a checked property.

Some headers are not capabilities: generic instantiations that Rust generics
replace, containers, superseded internals. Written triage rules remove them.
Each rule has a stable id and records the headers it removed, so the reduction
can be audited.

## Verdicts and evidence

Every capability has one of six verdicts. The gate enforces the evidence each
verdict requires:

| Verdict | Means | Must cite |
|---|---|---|
| `covered` | built and tested | symbols that resolve in the built rustdoc, and tests that exist in the tree |
| `partial` | built, with a stated restriction | the restriction in words, plus symbols and tests |
| `divergent` | deliberately different | the reasoning |
| `absent` | not built | a plan entry |
| `n/a` | excluded by a triage rule | the rule id |
| `unreviewed` | not yet audited | counted against a ratchet that can only go down |

The citations are checked on every build. Rename a cited symbol or delete a
cited test, and the audit fails. `docs/PARITY.md` is generated from the index
and the ledger, committed, and checked for staleness, so the rendered page
cannot drift from the data.

## Current state

As of the audit's completion, **no capability is absent**. The remaining
work is exactly what the ledger lists: the `partial` restrictions, and the
`unreviewed` count, which the ratchet holds at its floor. The
[ledger chapter](parity-ledger.md) shows the current committed state and is
rebuilt with the book.
