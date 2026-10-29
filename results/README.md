# Result tables: what backs what

Every number in `preprint/scCNV-FP-Evaluation.tex` is computed from the tables in this
directory. This file records which table backs which manuscript object, so that a reader
regenerating the figures does not have to guess.

| Table | Contents | Used for |
|---|---|---|
| `T0_datasets.csv` | panel composition (cells, donor--dataset counts, cell types, assay, disease) | Table 1, Table S6 |
| `T1_infercnv_fpr.csv` | infercnvpy false-call rate per dataset x reference mode x decision rule | Table 1, Table S1, Fig. 1a, 1b |
| `T2_ikarus_fpr.csv` | ikarus tumour-call rate per dataset | Table 1, Fig. 1a |
| `T2b_copykat_fpr.csv` | CopyKAT aneuploid and `not.defined` rate per dataset | Table 1, Fig. 1a |
| `T3_spike_auc.csv` | spike-in AUC by burden, tissue and reference mode | Table S2, Fig. 2a |
| `T3b_ikarus_burden.csv` | ikarus call rate by imposed burden | Table S3, Fig. 2b |
| `T3c_copykat_burden.csv` | CopyKAT call rate by imposed burden | Table S3, Fig. 2b |
| `T3d_copykat_spike.csv` | CopyKAT false-call and high-burden detection rate | Results text |
| `T4_matched_fpr.csv` | rule vs score at matched false-call rate | Table 2, Fig. 2c |
| `T5c_alpha_sweep.csv` | **leave-one-tissue-out (LOTO) calibration across alpha** | Table 3 (LOTO columns), Table S4, Fig. 4a, 4b, 4c |
| `T10_condition_matched.csv` | **donor-held-out calibration, cross-tissue vs same-dataset null** | Table 3 (donor-held-out columns), Fig. 4d |
| `T6_density_vs_bias.csv` | panel-averaged apparent shift and gene count per 10-Mb bin | Fig. 3b |
| `T6b_density_shift_test.csv` | circular-shift test of the gene-density association | Results, Methods, Fig. 3b |
| `T7_binned_profiles.csv` | per-dataset bin-level difference profiles | Fig. 3a, 3c; input to T8, T9b |
| `T8_cross_tissue_corr.csv` | pairwise Spearman correlation of difference profiles | Table S5, Fig. 3c |
| `T9_family_enrichment.csv` | clustered gene-family locus annotation (**annotation only, not a test**) | input to T9b |
| `T9b_family_enrichment_test.csv` | circular-shift test of family-locus enrichment | Results, Methods |
| `T9c_family_per_locus.csv` | per-family mean absolute difference and fold vs background | Results |

## Removed before posting: two superseded calibration tables

`T5_calibration_loto_labels.csv` and `T5b_calibration_loto_archetype.csv` were **earlier
calibration runs that no manuscript number was taken from**, and they were removed on
2026-10-27 so that nothing in this directory contradicts the published tables. Their file
names described the same quantity as the LOTO columns of Table 3, but their values differed
substantially from the reported ones (for example liver 4.31%, 129/2,995, against the
reported 0.00%, 0/2,995; and for `T5b` psoriatic skin 27.8% against the reported 57.0%).
`T5` carried no `alpha` column, so the level at which it was produced was never recorded,
and the driver that wrote both tables was not retained (see Limitations in the manuscript),
so neither could be regenerated or reconciled.

The calibration quantities the manuscript reports are those in `T5c` (LOTO) and `T10`
(donor-held-out), which are mutually consistent with Table 3, Table S4 and Fig. 4 as
published.

The known disagreement between `T5c` and the cross-tissue arm of `T10` affects three
datasets whose two arms score an identical set of cells (large intestine, liver, lung); it
is disclosed in the manuscript itself, and both are reported as computed. For 7 of the 13
datasets the two analyses do not score the same cells at all, because `T10` is restricted
to donors contributing at least 200 cells.
