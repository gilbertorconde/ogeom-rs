---
name: ci-iterate
description: "Loop on a PR or a pushed branch until the actionable checks pass and high and medium review feedback is fixed: classify each red check, reproduce locally, fix the cause, commit, repeat. Use for a red CI run, a failing or flaky check, a red nightly stress run, or review feedback on an open PR. Not for approval or merge gates, and not for a bug with no CI run behind it (use systematic-debugging)."
---

# ci-iterate

Fix actionable CI failures and high and medium review feedback on one PR or
branch. Stop and report when only human gates remain.

## The gates

| Workflow | Job | Local equivalent | Gating |
|---|---|---|---|
| `ci.yml` | `check` (fmt, comment rot, clippy, tests twice, docs, book, parity) | `tools/check.sh` | yes |
| `ci.yml` | `msrv` | `cargo +<rust-version> check --workspace --all-targets` (`check.sh` runs it when the toolchain is installed) | yes |
| `ci.yml` | `deny` | `cargo deny check` | yes |
| `ci.yml` | `outside` | `cargo test` in `outside/` | no: information, never a reason to block a kernel change |
| `stress.yml` | nightly and on demand | `cargo run --release -p ogeom-stress -- --seed 1 --check tools/ogeom-stress/baseline.json` | yes, on its own schedule |
| `docs.yml` | on push to `main` | `mdbook build docs/book` and `cargo doc` (both inside `check.sh`) | deploys only |

`tools/check.sh` is the whole gating set on one machine. Run it bare: piping
it into `tail` reports the pipe's status, not the script's.

## Rules

- Tie every diagnosis to the exact run, job and head SHA. Evidence from an
  older SHA is not evidence for this one.
- Classify before you fix: a kernel defect, a test defect, a property test
  seed that found a real case, a flake, infrastructure (runner, a registry
  outage, a toolchain download), or a baseline that is simply stale.
- Fix the root cause, with the systematic-debugging skill. No `--no-verify`,
  no skipped or commented-out test, no weakened assertion, no loosened
  tolerance without its reason in the commit, no `lint-comment-rot: ignore`
  to silence a comment that narrates a change.
- Never force push, rerun or cancel runs, edit workflow files, or post on
  the PR unless the user asks. Do not push at all unless the user has said
  to; commit, gate, and hand off.
- Same failure after two attempts: stop and ask.

## Loop

1. Identify the branch, PR and state:

   ```bash
   gh pr view --json number,url,headRefName,headRefOid,isDraft,reviewDecision
   gh pr checks --json name,state,bucket,link,workflow
   gh run list --branch <branch> --limit 10 --json databaseId,workflowName,conclusion,headSha
   ```

2. Feedback: read review threads and comments.

   ```bash
   gh api repos/{owner}/{repo}/pulls/{number}/comments
   gh api repos/{owner}/{repo}/pulls/{number}/reviews
   ```

   Fix high and medium items (wrong geometry, broken invariant, a rule from
   `CONTRIBUTING.md`) after verifying each one against the code. List low
   items (naming, taste) for the user to choose. A false positive gets a
   short reason, not a code change.

3. Checks:

   | State | Action |
   |---|---|
   | a check failed, none pending | fix failures |
   | actionable checks pending | wait and read feedback meanwhile |
   | only review or approval gates pending | report `BLOCKED_BY_REVIEW_GATE` |
   | no checks after a grace period | report `NO_CHECKS` |
   | all actionable checks passed | read feedback once more, then stop |

4. Fix each failure:

   ```bash
   gh run view <run-id> --log-failed
   gh run view <run-id> --json headSha,jobs
   ```

   Trace from the assertion, panic, lint or compiler error to its source.
   State the cause in one line before editing ("fails because X, reached by
   Y"). Search sibling call sites for the same defect and fix all of them.

   Failures with their own shape:
   - **A property test.** The log names the seed and the minimal case. Add
     the case to the regression file or as a named test, then fix the owner.
     The suite runs twice in `check.sh` for exactly this reason: one green
     run is not proof.
   - **The parity gate.** A verdict in `docs/parity/parity.toml` cites a
     symbol or test that no longer exists, or `docs/PARITY.md` is stale. Fix
     the citation or regenerate the ledger; never delete the verdict to pass.
   - **The comment-rot lint.** The comment narrates a change or cites an
     issue. Rewrite it to describe present behaviour.
   - **MSRV.** A construct the newer compiler accepts. Rewrite for the
     declared `rust-version`; raising it is a decision for the user.
   - **Stress.** A scenario has fewer `ok` or more `invalid`, `wrong` or
     `panic` than the baseline. Replay the case with `--case`, fix the
     kernel. A stress baseline is rewritten (`--write`) only when the
     improvement is intended and stated.
   - **The bench step** is informational. A drift is reported, not fixed
     here (use performance-optimizer).

5. Flake suspected: measure before believing it.

   ```bash
   gh run list --workflow=ci.yml --branch main --status completed \
     --limit 60 --json conclusion,headSha,createdAt,databaseId
   ```

   Confirm it is the same test each time from the failing job's log.
   Interleaved pass and fail on unchanged code is a flake; a long unbroken
   red streak is a regression, so find the last green and first red commit.
   In this kernel a flake is nearly always a property test finding a real
   case on some seeds, or a timing in a tool; reproduce it with the seed
   before calling it a flake.

6. Verify locally with `tools/check.sh`, commit in the repo convention (one
   scope per commit, `fix(<area>): <what the user gets>`, as in `fix(bool)`), and hand off
   for the push. Restart at step 3 once the run is in.

## Exit

| Result | When |
|---|---|
| `DONE` | actionable checks pass and no high or medium feedback is open |
| `ASK` | the same failure after two attempts, unclear feedback, or an infrastructure issue |
| `STOP` | no PR, only human gates remain, or the branch needs a conflict decision |

Report the branch, head SHA, each failure with its class, cause and fix,
the local checks run, and what remains.
