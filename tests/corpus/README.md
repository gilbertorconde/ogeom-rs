# The corpus

Outside-authored input for the verification strategy (`docs/SCOPE.md`,
*Verification*): files chosen by people who were not us, so the tests stop
measuring the authors' imagination. Consumed by the STEP importer when §17
lands; until then this directory holds the files and their provenance, ready.

## Rules

1. **Files are copied in, never referenced from outside the repository.**
   A test that reads a path outside the tree works on one machine.
2. **Every file has a manifest entry** in `MANIFEST.md`: filename, source URL,
   the license or legal basis for redistribution, and the date it was fetched.
   A file without an entry does not get committed.
3. **No file from another CAD kernel's own test suite.** The independence rule
   covers test assets that ship inside a kernel's repository; a neutral STEP
   file that happens to have passed through some kernel is fine — that is what
   STEP is for.

## Vetted sources

Licence notes below are a starting point for validation, not a substitute for
it — check the statement on the actual download page before committing files.

### NIST MBE PMI validation models (CTC / FTC test cases) — recommended

- What: STEP AP242 test models (also AP203/AP214 variants) designed
  specifically to exercise conformance — PMI, geometry, assemblies. Small,
  curated, and designed as *test cases*, which is exactly the role here.
- Where: NIST "MBE PMI Validation and Conformance Testing" project,
  nist.gov (search: NIST PMI CTC FTC STEP models).
- Legal basis: works of the U.S. federal government are not subject to
  copyright in the United States (17 U.S.C. § 105). NIST distributes them
  with a no-warranty disclaimer and a request to credit NIST as the source.
  **Verify on the download page** that the specific models carry the standard
  NIST statement and no contractor copyright notice, then record the exact
  wording in the manifest.

### STEPcode sample files — acceptable, verify per file

- What: assorted STEP files in the `stepcode` project's test data.
- Legal basis: the project is BSD-3-Clause; files inside the repo are
  presumptively under it, which is compatible with MIT OR Apache-2.0.
  Verify no per-file notices say otherwise.

## Sources considered and rejected

- **ABC dataset** (NYU, ~1M models scraped from Onshape): each model's
  copyright belongs to its Onshape author; the platform terms that made them
  public do not clearly grant redistribution in a permissively licensed
  repository. Widely used in academic ML work, which is a different legal
  posture from vendoring. Rejected as murky.
- **Fusion 360 Gallery dataset** (Autodesk): research licence with
  non-commercial terms. Incompatible with this repository's licence. Rejected.
- **GrabCAD and similar model-sharing sites**: per-model licences, generally
  no redistribution right. Rejected.
- **Any kernel's bundled test corpus**: independence rule. Rejected.
