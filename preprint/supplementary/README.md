# Supplementary Data

Supplementary Data 1-18 are the result tables underlying every figure, table and number in
the manuscript. Each file is a byte-identical copy of a table in `results/` (the "Source"
column), renamed for submission. Supplementary Tables S1-S6 are typeset in the manuscript PDF;
the files below are their machine-readable sources.

| File | Source | Contents | Used for |
|---|---|---|---|
| `Supplementary_Data_01_panel_datasets.csv` | `T0_datasets.csv` | Panel composition (cells, donor-dataset counts, cell types, assay, disease) | Table 1, Table S6 |
| `Supplementary_Data_02_infercnvpy_false_call_rates.csv` | `T1_infercnv_fpr.csv` | infercnvpy false-call rate per dataset x reference mode x decision rule | Table 1, Table S1, Fig. 1a, 1b |
| `Supplementary_Data_03_ikarus_false_call_rates.csv` | `T2_ikarus_fpr.csv` | ikarus tumor-call rate per dataset | Table 1, Fig. 1a |
| `Supplementary_Data_04_copykat_false_call_rates.csv` | `T2b_copykat_fpr.csv` | CopyKAT aneuploid and `not.defined` rate per dataset | Table 1, Fig. 1a |
| `Supplementary_Data_05_spikein_auc.csv` | `T3_spike_auc.csv` | Spike-in AUC by burden, tissue and reference mode | Table S2, Fig. 2a |
| `Supplementary_Data_06_ikarus_by_karyotype_burden.csv` | `T3b_ikarus_burden.csv` | ikarus call rate by imposed burden | Table S3, Fig. 2b |
| `Supplementary_Data_07_copykat_by_karyotype_burden.csv` | `T3c_copykat_burden.csv` | CopyKAT call rate by imposed burden | Table S3, Fig. 2b |
| `Supplementary_Data_08_copykat_spikein_detection.csv` | `T3d_copykat_spike.csv` | CopyKAT false-call and high-burden detection rate | Results text |
| `Supplementary_Data_09_matched_false_call_rate_comparison.csv` | `T4_matched_fpr.csv` | Published rule vs continuous score at matched false-call rate | Table 2, Fig. 2c |
| `Supplementary_Data_10_calibration_leave_one_tissue_out_alpha_sweep.csv` | `T5c_alpha_sweep.csv` | Leave-one-tissue-out calibration across alpha | Table 3 (LOTO columns), Table S4, Fig. 4a-c |
| `Supplementary_Data_11_calibration_donor_held_out.csv` | `T10_condition_matched.csv` | Donor-held-out calibration, cross-tissue vs same-dataset null | Table 3 (donor-held-out columns), Fig. 4d |
| `Supplementary_Data_12_gene_density_vs_apparent_shift.csv` | `T6_density_vs_bias.csv` | Panel-averaged apparent shift and gene count per 10-Mb bin | Fig. 3b |
| `Supplementary_Data_13_gene_density_circular_shift_test.csv` | `T6b_density_shift_test.csv` | Circular-shift test of the gene-density association | Results, Methods, Fig. 3b |
| `Supplementary_Data_14_binned_difference_profiles.csv` | `T7_binned_profiles.csv` | Per-dataset bin-level difference profiles | Fig. 3a, 3c |
| `Supplementary_Data_15_cross_tissue_profile_correlation.csv` | `T8_cross_tissue_corr.csv` | Pairwise Spearman correlation of difference profiles | Table S5, Fig. 3c |
| `Supplementary_Data_16_gene_family_locus_annotation.csv` | `T9_family_enrichment.csv` | Clustered gene-family locus annotation (annotation only, not a test) | Input to Supplementary Data 17 |
| `Supplementary_Data_17_gene_family_enrichment_test.csv` | `T9b_family_enrichment_test.csv` | Circular-shift test of family-locus enrichment | Results, Methods |
| `Supplementary_Data_18_gene_family_per_locus.csv` | `T9c_family_per_locus.csv` | Per-family mean absolute difference and fold vs background | Results |
