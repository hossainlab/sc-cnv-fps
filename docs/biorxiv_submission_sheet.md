# bioRxiv submission sheet — null-control audit of malignant-cell calling

Prepared to accompany `preprint/scCNV-FP-Evaluation.pdf`. Everything below is either ready
to paste into the submission form or flagged as needing your input. Last reconciled against
the manuscript on 2026-10-27.

## Files to upload

| Form field | File | Note |
|---|---|---|
| Manuscript | `preprint/scCNV-FP-Evaluation.pdf` | Single PDF: text, 4 figures, 3 tables and Supplementary Tables S1–S6; 23 pages. Well within the 40 MB limit. |
| Supplemental material | `scCNV-FP-Evaluation_supplement.zip` | Result tables, analysis and figure scripts, figure files, LaTeX source. Built by `make supplement`. |

Line numbers are enabled in the current source to aid review (`\linenumbers`, flagged with a
comment in the preamble). Remove that line for final typesetting if you prefer a clean page.

The DOCX route is an alternative, not an addition: upload either the PDF, or a Word file plus
separate figure files for server-side conversion. Do not upload both. bioRxiv requires TeX to
be converted to PDF before submission; the `.tex` source travels in the supplement zip.

## Metadata fields

**Title**
Evaluating the false-positive rate of single-cell copy-number and malignant-cell classifiers

**Abstract for the metadata box** (about 320 words; the in-PDF abstract is longer)

Single-cell cancer studies split cells into malignant and non-malignant compartments using
copy-number inference or supervised classifiers, but these methods are rarely run on tissue
containing no tumour. We assembled 14 tumour-free scRNA-seq datasets (36,528 cells from 76
donor–dataset combinations; 10 healthy organs and 4 inflammatory-disease datasets) and ran
inferCNV (as infercnvpy), CopyKAT and ikarus under published decision rules; every positive
call is false by construction. infercnvpy and CopyKAT made false calls in every dataset and
ikarus in 12 of 14: 0.5–19.3% for the two-criterion rule of Tirosh and colleagues, 4.3–56.2%
for a forced two-group split with an immune reference and up to 78.0% without one, 3.0–63.0%
for CopyKAT and 0–100% for ikarus, with CopyKAT and ikarus each calling about 31% of cells
pooled. In-silico positive controls showed per-cell discrimination near chance for single arm
events (AUC 0.53), rising to 0.85 at 18 events. At a false-call rate matched on the same
diploid cells used to measure it, so that the factor is an upper bound, the continuous score
detected 2.2–5.3 times as many altered cells as the published rules in all eight estimable
rule comparisons, and up to 10.3 times as many as an external classifier. False calls share a
pseudo-karyotype that recurs across tissues (median Spearman
rho = 0.54) and tracks local gene density (rho = -0.53, P < 10^-4 against a within-chromosome
circular-shift null): gene-dense regions look deleted. A cell-type-conditional empirical null
with Benjamini–Hochberg control at alpha = 5% cut the per-cell false-call rate to a median of
0.03% in 13 of 14 datasets (28 false calls among 33,536 cells), although 7 of these 13 still
received at least one call, so extreme-tail control was not exact, and sensitivity was not
evaluated. In psoriatic skin from a separate study, a null from other tissues failed (26.8%);
in a separate donor-held-out evaluation the corresponding cross-tissue null called 27.1% of
cells while a null from held-out donors of the same study called 0.03%.

**Keywords**
single-cell RNA-seq; copy-number inference; malignant-cell classification; negative controls;
false-positive rate; calibration

**Suggested subject category**
Bioinformatics (primary). Cancer Biology or Genomics are defensible alternatives; the paper is
a methods audit, so the bioinformatics category matches the audience most likely to act on it.

**License**
CC BY 4.0 if you want the calibration layer reused and redistributed freely, which is the
point of releasing it. CC BY-NC-ND is the most common choice on the server but forbids
derivative works, which sits awkwardly with a paper whose contribution is a reusable layer.

## Declarations — drafted, confirm before posting

These are written into the manuscript. They are statements only you can make, so read them and
correct anything that is not true of your situation:

- **Competing interests**: "The authors declare no competing interests."
- **Funding**: "This work received no specific funding."
- **Author contributions**: now given as CRediT roles in the manuscript. Confirm the split
  between the two authors is accurate before posting.
- **Use of AI tools**: a dedicated section states that AI-assisted tools were used for data
  analysis, figure preparation, writing and critical review under author direction, and that
  the authors verified all outputs and take responsibility for the content.

## Still needed from you

1. **ORCID iD per author** (optional, recommended).
2. **Affiliation detail** for DeepBio (city, country) and confirmation of the CHIRAL
   Bangladesh line.
3. **Confirm the CRediT role assignment** now written into the manuscript.

## Before you click submit

- Posting is one-way: a preprint cannot be removed from bioRxiv once posted. Revisions are
  supported and keep the original DOI, but the first version stays visible.
- Every submission is screened in-house for completeness, scope and article type, and
  automatically for plagiarism and prior online appearance. Nothing here has been posted
  elsewhere.
- **The objection to expect.** The calibrated rule's specificity is demonstrated (about
  90-fold fewer false calls) but its sensitivity was never measured, so the proposed procedure
  is not yet shown to retain power. On tumour-free data a rule that calls almost nothing wins
  trivially. This is stated in the abstract, Results and Limitations; it is a scope decision,
  not a formatting one, and a journal would likely require it to be closed.
- The second objection to expect is that no allele-aware caller (Numbat, HoneyBADGER, CaSpER)
  was run, even though the paper predicts such methods should resist the compositional
  artefact it describes. Also disclosed.
