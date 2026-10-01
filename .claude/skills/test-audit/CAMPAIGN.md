# Test-pruning campaign

Campaign mode prunes one crate's whole test surface in one batch, such as
`ogeom-io` or `ogeom-bool`. The value bar, retention bar, candidate
evidence and validation in [SKILL.md](SKILL.md) apply to every lane. This
file adds the order of work. Each step ends on its completion criterion; do
not start the next step early.

## 1. Baseline

Record the crate's test line count and every test file's pass or fail state
at a pinned `main` SHA. Keep baseline failures in their own list: treat each
one as a possible kernel bug, not a stale test.

Done when every in-scope test file has a recorded baseline result.

## 2. Lanes and inventory

Split the surface into lanes along production owner boundaries, not file
prefixes. For `ogeom-io` these could be the STEP reader, the STEP writer,
IGES, the mesh formats, and the corpus suites. Include the crate's cases in
`crates/ogeom/tests` and its citations in `docs/parity/parity.toml`.

Done when every test file belongs to exactly one lane.

## 3. Read-only ledger per lane

Give each lane to its own read-only agent. The agent reads every assigned
test in full, including tables and property strategies, the production
owners, their callers and history. Each test goes into a written ledger with
one mark:

- `R`: retain, naming the contract and the bug it catches;
- `F`: retain the contract but repair the assertion, such as a refusal test
  that accepts any error;
- `C`: consolidate, naming the owner that absorbs the assertion first: a
  sibling table case, a property, or a stronger suite;
- `D`: delete, naming the proof that remains, or why no contract exists.

Judge a test by its assertions, not its name.

Done when every test in the lane has a mark and an evidence line.

## 4. Layer plan per lane

A second read-only pass, starting from the ledger, looks for the redundant
layer: several suites that assert the same closed form through the same
primitive. Name the keeper for each contract. Correct any ledger errors
this pass finds.

Done when each lane plan names its retired files, its keeper per contract,
the assertions to carry into keepers, and the test-only items unlocked.

## 5. Cutover

Edit lane by lane. With each lane, remove the test-only items it unlocks.
Re-evidence any parity verdict that cited a retired test before retiring
it. Keep book anchors intact.

Done when every lane plan is applied, each lane's keepers pass, and
`tools/parity.py check` passes.

## 6. Preservation review

Before claiming completion, have independent reviewers compare deleted
coverage against the keepers, one reviewer per lane. They look for
contracts that lost their only proof, and for new assertions that cannot
fail.

For each restored contract, make one deliberate mutation of the production
owner and confirm the keeper goes red. Then restore the source byte for
byte.

Done when every reported gap is restored or rejected with source evidence,
and every restored contract has a caught mutation.

## 7. Kernel defects

A baseline failure that survives into a keeper is a bug. Fix it at its
owner as a separate commit, with a pinning test that fails on the pre-fix
code. Record unrelated discrepancies as follow-ups instead of fixing them in
the campaign.

Done when each repaired defect has a failing pre-fix run and a passing
post-fix run.

## 8. Hand off

Hand off with the [SKILL.md](SKILL.md) report, plus: baseline and final test
line counts, with production counted separately; lanes, retired layers and
keepers; preservation gaps found and their mutations; defects with their
pinning tests.
