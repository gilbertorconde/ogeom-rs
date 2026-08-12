# The parity audit

"How complete is this kernel?" is a question projects usually answer with
adjectives. ogeom answers it with a gated audit, and this chapter explains
the machinery so [the ledger](parity-ledger.md) reads as what it is: a
measurement.

## The shape of the claim

The reference target is fixed by [Scope](scope.md): four modelling
modules, 276 packages, 6,267 public headers, committed as
`docs/parity/reference-index.tsv` so the question is answerable without a
reference checkout. Against that index, **capabilities are the primary
key, not headers**: each of the 97 capabilities in
`docs/parity/parity.toml` *claims* the reference headers it accounts for,
and the gate requires the claims to be **total and disjoint** — every
kept header claimed exactly once, or the build fails. That is what makes
"nothing was forgotten" a checked property instead of a hope.

Headers that are not capabilities — generic instantiations Rust's
generics subsume, containers, superseded internals — are removed by
written triage rules, each with a stable id and its removed headers
recorded, so the reduction is auditable rather than asserted.

## Verdicts and evidence

Every capability carries one of six verdicts, and each verdict has an
evidence requirement the gate enforces:

| Verdict | Means | Must cite |
|---|---|---|
| `covered` | built and tested | symbols that resolve in the built rustdoc, tests that exist in the tree |
| `partial` | built, with a stated restriction | the restriction, in words, plus symbols and tests |
| `divergent` | deliberately different | the reasoning |
| `absent` | not built | a plan entry |
| `n/a` | excluded by a triage rule | the rule id |
| `unreviewed` | not yet audited | counted against a ratchet that only goes down |

The citations are live: rename a symbol and the audit fails the build;
delete a cited test, same. `docs/PARITY.md` is generated from the index
and the ledger, committed, and checked for staleness — the rendered form
cannot drift from the data.

## Where it stands

As of the audit's completion, **no capability is absent**. The residual
worklist is exactly what the ledger states: the `partial` restrictions,
and the `unreviewed` count the ratchet holds at its floor. The
[ledger chapter](parity-ledger.md) is the current committed state,
rebuilt with the book.
