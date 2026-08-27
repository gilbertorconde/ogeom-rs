# Corpus manifest

One entry per committed file. A file without an entry does not get committed.

| File | Source (URL) | Licence / legal basis | Fetched |
|---|---|---|---|
| `nist_ctc_01_asme1_rd.stp` | [NIST-PMI-STEP-Files.zip](https://www.nist.gov/system/files/documents/noindex/2024/06/19/NIST-PMI-STEP-Files.zip), `AP203 geometry only/`, via the [NIST download page](https://www.nist.gov/ctl/smart-connected-systems-division/smart-connected-manufacturing-systems-group/mbe-pmi-0) | US government work, 17 U.S.C. §105; page states verbatim: "The test cases, CAD models, and STEP files can be used without any restrictions." NIST appreciates acknowledgement; NIST logo not to be used promotionally. | 2026-08-04 |
| `nist_ctc_02_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ctc_03_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ctc_04_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ctc_05_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_06_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_07_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_08_asme1_rc.stp` | same | same | 2026-08-04 |
| `nist_ftc_09_asme1_rd.stp` | same | same | 2026-08-04 |
| `nist_ftc_10_asme1_rb.stp` | same | same | 2026-08-04 |
| `nist_ftc_11_asme1_rb.stp` | same | same | 2026-08-04 |
| `NIST-README.txt` | NIST's own readme from inside the zip, kept verbatim as provenance | same | 2026-08-04 |
| `nist_ctc_01_asme1_ap242-e1.stp` | Same NIST zip as the AP203 parts, `NIST-PMI-STEP-Files/` root: the CTC 1 part in AP242 edition 1 with full semantic PMI — dimensions, plus/minus bounds, geometric tolerances, datums. | Same NIST terms as above | 2026-08-05 |
| `ogeom_asm_bolted_plate.stp` | Authored for this project (generated, then committed): a plate with two bolts — AP214 product structure, three usage occurrences over two parts, per-product and per-face colours, reference designators. Ground truth for the assembly reader, every value chosen by hand. | This repository's own licence (MIT OR Apache-2.0) | 2026-08-05 |

The eleven `.stp` files are the **AP203 geometry-only** exports of the CTC 1–5 and
FTC 6–11 test parts from the NIST MBE PMI Validation and Conformance Testing
project — part geometry with no PMI, which is the subset a geometry kernel reads
first. NIST's readme states plainly that these are **not** error-free reference
files: conformance checkers report syntax errors in them. That is part of their
value — a reader that only accepts clean files has not been tested. The
PMI-annotated AP242 variants exist in the same zip and can be added under the
same licence when PMI is in scope.

## m5x16_bhcs.step

An M5x16 button-head cap screw, five placed bodies, extracted standalone
from a community printer assembly (issue #37's fixture). Nothing but the
classic primitives — spherical button head, cylindrical shank, conical
chamfers — with the head's sphere zone slit along a meridian that sits at
the chart's own seam: the minimal reproducer for a doubly-used edge that
is a slit, not a period-wrapping seam.
