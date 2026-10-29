# Submission bundle — contents and how it was built

## Manuscript
- `preprint/scCNV-FP-Evaluation.pdf` — typeset preprint, 23 pages (bioRxiv upload file);
  built by `make preprint`
- `preprint/scCNV-FP-Evaluation.tex` — LaTeX source; figures resolve to `figures/Fig*.pdf`
  and generated tables to `tables/*.tex`
- `preprint/scCNV-FP-Evaluation.bib` — bibliography (31 references), typeset with BibTeX
  (`unsrtnat`, numbered in citation order)
- `docs/archive/` — dated pre-review copies of the source

## Figures
Fig. 1 false positives; Fig. 2 spike-in positive controls; Fig. 3 mechanism;
Fig. 4 calibration. PNG (330 dpi) and vector PDF versions of each are in `figures/`.
Rebuild with `make figures` (needs only `results/` and
`calls/TS_Lung_infercnv.parquet`).

## Tables
Main Tables 1–3 and Supplementary Tables S1–S6 are generated into `tables/`.
Table S6 lists the CELLxGENE dataset identifier for each panel dataset.

## Supplement
`scCNV-FP-Evaluation_supplement.zip` contains the analysis and figure scripts, the result
tables every manuscript number is read from (`results/T*.csv`, `results/*.json`), the figure
files and the LaTeX source. See `results/README.md` for which table backs which figure or
table, and for the two superseded calibration tables that no manuscript number uses.

## Reproducing
Callers: infercnvpy 0.6.1 (Python), ikarus 0.0.3 (Python), CopyKAT 1.2.5 (R 4.4.3).
Entry points: `run_infercnv.py`, `run_ikarus.py`, `export_for_r.py` + `run_copykat_all.R`.
Calibration layer: `calibrate.py` (cell-type-conditional empirical null with generalized-Pareto
tail and per-sample Benjamini–Hochberg control). Figures: `fig1.py` … `fig4.py`.

Not everything reruns end to end; see "Known gaps" in the top-level `README.md`. In
particular the table-aggregation driver, `chrom_arms.csv`, the recorded Ensembl release and
the per-cell spike-in scores were not retained.
