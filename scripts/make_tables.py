"""Manuscript tables (LaTeX) and derived numbers, built only from results/*.csv.

Writes tables/*.tex and results/derived_numbers.json. Every number quoted in
the preprint text that is not a direct table lookup is computed here, so the text,
tables and figures share one source.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path("results")
OUT = Path("tables")
OUT.mkdir(parents=True, exist_ok=True)

PRETTY = {"TS_Bone_Marrow": "Bone marrow", "TS_Kidney": "Kidney",
          "TS_Large_Intestine": "Large intestine", "TS_Liver": "Liver", "TS_Lung": "Lung",
          "TS_Pancreas": "Pancreas", "TS_Skin": "Skin", "TS_Small_Intestine": "Small intestine",
          "TS_Spleen": "Spleen", "TS_Stomach": "Stomach", "UC_colon": "Colon (UC cohort)",
          "PSO_skin": "Skin (psoriasis)", "CD_epithelial": "Colon (Crohn's, epithelial)",
          "CD_immune": "Colon (Crohn's, immune)"}
ORDER = list(PRETTY)
RULES = {"R1_twogroup": "Forced two-group split", "R1b_topcluster": "Highest-scoring cluster",
         "R3_refSD": "Reference SD", "R4_tirosh": "Two-criterion (Tirosh)"}
SPK = {"SPK_TS_Liver": "Liver", "SPK_TS_Lung": "Lung", "SPK_TS_Spleen": "Spleen"}


def wilson(k, n, a=0.05):
    if n == 0:
        return np.nan, np.nan
    z = stats.norm.ppf(1 - a / 2)
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def f1(x, nd="--"):
    return nd if pd.isna(x) else f"{x:.1f}"


def f2(x, nd="--"):
    return nd if pd.isna(x) else f"{x:.2f}"


def write(name, body):
    (OUT / f"{name}.tex").write_text(body)
    print("wrote", OUT / f"{name}.tex")


T0 = pd.read_csv(R / "T0_datasets.csv").set_index("dataset")
T1 = pd.read_csv(R / "T1_infercnv_fpr.csv")
T2 = pd.read_csv(R / "T2_ikarus_fpr.csv").set_index("dataset")
T2b = pd.read_csv(R / "T2b_copykat_fpr.csv").set_index("dataset")
T3 = pd.read_csv(R / "T3_spike_auc.csv")
T3b = pd.read_csv(R / "T3b_ikarus_burden.csv")
T3c = pd.read_csv(R / "T3c_copykat_burden.csv")
T4 = pd.read_csv(R / "T4_matched_fpr.csv")
T5c = pd.read_csv(R / "T5c_alpha_sweep.csv")
T6 = pd.read_csv(R / "T6_density_vs_bias.csv")
T8 = pd.read_csv(R / "T8_cross_tissue_corr.csv", index_col=0)
T10 = pd.read_csv(R / "T10_condition_matched.csv")

D = {}  # derived numbers

# ------------------------------------------------------------------ panel + callers
imm = T1[T1["mode"] == "refimm"].set_index("dataset")
ref = T1[T1["mode"] == "reffree"].set_index("dataset")
fig_mode = {d: ("refimm" if d in imm.index else "reffree") for d in ORDER}
tir = imm["R4_tirosh"].combine_first(ref["R4_tirosh"])
two = imm["R1_twogroup"].combine_first(ref["R1_twogroup"])

ck_n = T2b["n"]
ck_aneu = T2b["aneuploid_pct"] / 100 * ck_n
ck_nd = T2b["notdefined_pct"] / 100 * ck_n

rows = []
for d in ORDER:
    t = T0.loc[d]
    mark = "" if fig_mode[d] == "refimm" else "$^{\\dagger}$"
    rows.append(" & ".join([
        PRETTY[d], t["assay"].replace("'", "$'$"), f"{t.n_cells:,}", str(t.n_donors),
        str(t.n_cell_types), f1(tir[d]) + mark, f1(two[d]) + mark,
        f1(T2b.loc[d, "aneuploid_pct"]), f1(T2b.loc[d, "notdefined_pct"]),
        f"{int(ck_n[d]):,}", f1(T2.loc[d, "tumour_pct"])]) + r" \\")
    if d == "TS_Stomach":
        rows.append(r"\midrule")

n_all = T0.n_cells
pool = lambda s: (s * n_all[s.index]).sum() / n_all[s.index].sum()
D["pooled_tirosh_fig"] = round(pool(tir), 1)
D["pooled_twogroup_fig"] = round(pool(two), 1)
D["copykat_pooled"] = round(100 * ck_aneu.sum() / ck_n.sum(), 1)
D["copykat_pooled_defined"] = round(100 * ck_aneu.sum() / (ck_n - ck_nd).sum(), 1)
D["copykat_n_total"] = int(ck_n.sum())
D["ikarus_pooled"] = round((T2.tumour_pct * T2.n).sum() / T2.n.sum(), 1)
D["n_cells"] = int(n_all.sum())
D["n_donors"] = int(T0.n_donors.sum())
D["tirosh_range_fig"] = [round(tir.min(), 1), round(tir.max(), 1)]
D["twogroup_range_fig"] = [round(two.min(), 1), round(two.max(), 1)]
D["twogroup_range_reffree"] = [round(ref.R1_twogroup.min(), 1), round(ref.R1_twogroup.max(), 1)]
D["copykat_range"] = [round(T2b.aneuploid_pct.min(), 1), round(T2b.aneuploid_pct.max(), 1)]
D["copykat_n_over20"] = int((T2b.aneuploid_pct > 20).sum())
D["ikarus_n_over50"] = int((T2.tumour_pct > 50).sum())
for m, s in [("reffree", ref), ("refimm", imm)]:
    D[f"pooled_rules_{m}"] = {r: round(pool(s[r].dropna()), 1) for r in RULES if s[r].notna().any()}

pooled_row = " & ".join([r"\textbf{Pooled}", "", f"{D['n_cells']:,}", str(D["n_donors"]), "",
                         f1(D["pooled_tirosh_fig"]), f1(D["pooled_twogroup_fig"]),
                         f1(D["copykat_pooled"]), f1(100 * ck_nd.sum() / ck_n.sum()),
                         f"{D['copykat_n_total']:,}", f1(D["ikarus_pooled"])]) + r" \\"

write("table1_panel", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Tumour-free panel and false malignant-call rates.} Every call is a false positive by construction. infercnvpy columns give the percentage of cells called malignant by infercnvpy under the two-criterion (Tirosh) rule and the forced two-group split, using the immune reference where one was available ($^{\dagger}$reference-free mode: no immune cells in the Crohn's epithelial dataset). CopyKAT was run on a cell-type-proportional subsample ($n$ shown); the aneuploid percentage uses all submitted cells as the denominator, with \texttt{not.defined} cells counted as not aneuploid. Pooled rates are weighted by dataset size. The pooled donor count is a sum over datasets; Tabula Sapiens donors contribute to several organs and the two Crohn's datasets come from one cohort, so it overstates the number of unique donors. Rates for all rules and both reference modes are in Supplementary Table~S1.}
\label{tab:panel}
\scriptsize
\setlength{\tabcolsep}{3.2pt}
\begin{tabular}{@{}llrrrrrrrrr@{}}
\toprule
 & & & & & \multicolumn{2}{c}{infercnvpy (\%)} & \multicolumn{3}{c}{CopyKAT} & ikarus \\
\cmidrule(lr){6-7}\cmidrule(lr){8-10}\cmidrule(lr){11-11}
Dataset & Assay & Cells & Donors & Types & Tirosh & Two-group & Aneu.\ (\%) & n.d.\ (\%) & $n$ & Tumour (\%) \\
\midrule
""" + "\n".join(rows) + "\n\\midrule\n" + pooled_row + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ S1: all rules x modes
rows = []
for d in ORDER:
    for m, s in [("Reference-free", ref), ("Immune reference", imm)]:
        if d not in s.index:
            rows.append(" & ".join([PRETTY[d], m] + ["n.a."] * 4) + r" \\")
            continue
        rows.append(" & ".join([PRETTY[d], m] + [f1(s.loc[d, r], "n.d.") for r in RULES]) + r" \\")
for m, key in [("Reference-free", "reffree"), ("Immune reference", "refimm")]:
    p = D[f"pooled_rules_{key}"]
    rows.append(" & ".join([r"\textbf{Pooled}", m] + [f1(p.get(r, np.nan), "n.d.") for r in RULES]) + r" \\")
write("tableS1_rules", r"""\begin{table}[htbp]
\centering
\caption{\textbf{False malignant-call rate (\%) of infercnvpy under four published decision rules and two reference modes.} The per-cell quantile rule (90th percentile) is omitted because it calls 10\% of cells by construction. n.d.: the reference-SD rule is undefined without a reference. n.a.: no immune reference cells available. Pooled rates are weighted by dataset size over the datasets where the rule is defined.}
\label{tab:s1}
\scriptsize
\begin{tabular}{@{}llrrrr@{}}
\toprule
Dataset & Mode & Two-group & Top cluster & Ref.\ SD & Tirosh \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ spike-ins
burdens = sorted(T3.burden.unique())
rows = []
for m, lab in [("reffree", "Reference-free"), ("refimm", "Immune reference")]:
    for ds in SPK:
        s = T3[(T3["mode"] == m) & (T3.dataset == ds)].set_index("burden").auc_meansq
        rows.append(" & ".join([lab, SPK[ds]] + [f"{s[b]:.2f}" for b in burdens]) + r" \\")
    s = T3[T3["mode"] == m].groupby("burden").auc_meansq.mean()
    rows.append(" & ".join([lab, r"\textit{Mean}"] + [f"\\textit{{{s[b]:.2f}}}" for b in burdens]) + r" \\")
    if m == "reffree":
        rows.append(r"\midrule")
    D[f"auc_mean_{m}"] = {int(b): round(s[b], 3) for b in burdens}
n_spk = int(T3[T3["mode"] == "reffree"].n.sum())
n_dip = int(T3b[T3b.spike_burden == 0].n.sum())
D["n_spike_aneuploid"], D["n_spike_diploid"] = n_spk, n_dip
write("tableS2_auc", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Per-cell discrimination of in-silico aneuploid cells from diploid companions.} Area under the ROC curve for the continuous infercnvpy score (mean squared deviation of the window profile) by number of imposed arm-level events. Each entry compares 120 aneuploid cells with the 2,155--2,162 untouched diploid companions of the same tissue (""" + f"{n_spk:,}" + r""" aneuploid and """ + f"{n_dip:,}" + r""" diploid cells in total).}
\label{tab:s2}
\scriptsize
\begin{tabular}{@{}ll""" + "r" * len(burdens) + r"""@{}}
\toprule
 & & \multicolumn{""" + str(len(burdens)) + r"""}{c}{Arm-level events per cell} \\
\cmidrule(lr){3-""" + str(2 + len(burdens)) + r"""}
Mode & Tissue & """ + " & ".join(str(b) for b in burdens) + r""" \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

rows = []
for ds in SPK:
    ik = T3b[T3b.dataset == ds].set_index("spike_burden").tumour_pct
    ck = T3c[T3c.dataset == ds].set_index("spike_burden")
    rows.append(" & ".join([SPK[ds], "ikarus"] + [f1(ik[b]) for b in [0] + burdens]) + r" \\")
    rows.append(" & ".join(["", "CopyKAT"] + [f1(ck.aneu_pct[b]) for b in [0] + burdens]) + r" \\")
    rows.append(" & ".join(["", r"\quad CopyKAT $n$"] + [str(int(ck.n[b])) for b in [0] + burdens]) + r" \\")
write("tableS3_callers_burden", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Response of ikarus and CopyKAT to imposed karyotype burden.} Percentage of cells called tumour (ikarus) or aneuploid (CopyKAT); burden 0 denotes untouched diploid companion cells. ikarus was applied to all cells (""" + f"{n_dip:,}" + r""" diploid, 120 per burden per tissue); CopyKAT to a subsample (cell counts shown).}
\label{tab:s3}
\scriptsize
\begin{tabular}{@{}ll""" + "r" * (len(burdens) + 1) + r"""@{}}
\toprule
Tissue & Caller & """ + " & ".join(str(b) for b in [0] + burdens) + r""" \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ matched FPR (per tissue)
rows, n_better, n_cmp = [], 0, 0
for r in ["R4_tirosh", "R1b_topcluster", "R1_twogroup", "ikarus"]:
    for ds in SPK:
        a = T4[(T4.dataset == ds) & (T4.method == r)]
        b = T4[(T4.dataset == ds) & (T4.method == f"score@FPR({r})")]
        if a.empty:
            continue
        a = a.iloc[0]
        if b.empty:
            sc, fold = "n.e.", "--"
        else:
            b = b.iloc[0]
            sc = f1(b.sens_hi)
            fold = f"{b.sens_hi / a.sens_hi:.1f}$\\times$" if a.sens_hi > 0 else "--"
            n_cmp += 1
            n_better += int(b.sens_hi > a.sens_hi)
        name = RULES.get(r, "ikarus")
        rows.append(" & ".join([name, SPK[ds], f1(a.fpr), f1(a.sens_hi), sc, fold]) + r" \\")
    rows.append(r"\addlinespace")
D["matched_score_better"] = [n_better, n_cmp]
write("table2_matched", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Detection at matched false-call rate, per spike-in tissue.} For each rule, its false-call rate on diploid companion cells and its detection rate on cells carrying 12--18 arm events, compared with a threshold on the continuous infercnvpy score placed at exactly that false-call rate. The threshold is a quantile of the score among the same diploid companions used to measure the rule's false-call rate, so the ratios are an upper bound attainable only where such a reference exists. n.e.: not estimable; the rule made no false calls on spleen diploid cells, so no score threshold at 0\% false calls was evaluated.}
\label{tab:matched}
\small
\begin{tabular}{@{}llrrrr@{}}
\toprule
 & & False-call & \multicolumn{2}{c}{Detection, 12--18 events (\%)} & \\
\cmidrule(lr){4-5}
Rule & Tissue & rate (\%) & Rule & Score & Ratio \\
\midrule
""" + "\n".join(rows[:-1]) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ------------------------------------------------------------------ mechanism
rho, p = stats.spearmanr(T6.n_genes, T6.mean_diff)
D["density_rho"], D["density_p"], D["density_nbins"] = round(rho, 2), float(f"{p:.1e}"), len(T6)
M = T8.values
iu = np.triu_indices_from(M, 1)
D["cross_rho_median"] = round(float(np.median(M[iu])), 2)
D["cross_rho_npairs"] = int(len(iu[0]))
sub = T8.drop(index="CD_immune", columns="CD_immune").values
D["cross_rho_median_excl_cdimm"] = round(float(np.median(sub[np.triu_indices_from(sub, 1)])), 2)
per = T8.apply(lambda r: r.drop(r.name).median(), axis=1)
D["cross_rho_per_dataset"] = per.round(2).to_dict()
cd = T8.loc["CD_immune"].drop("CD_immune")
D["cdimm_rho_range"] = [round(cd.min(), 2), round(cd.max(), 2)]

# ------------------------------------------------------------------ calibration
def sweep(src):
    p = T5c[T5c.source == src].pivot(index="dataset", columns="alpha", values="realized")
    n = T5c[T5c.source == src].groupby("dataset").n.first()
    return p, n

lab, nlab = sweep("labels")
arc, _ = sweep("archetype")
alphas = list(lab.columns)
nfp = (lab.mul(nlab, axis=0)).round().astype(int)
D["cal_alpha05_median_13"] = round(100 * lab[0.05].drop("PSO_skin").median(), 2)
D["cal_alpha05_max_13"] = round(100 * lab[0.05].drop("PSO_skin").max(), 2)
D["cal_alpha05_pso"] = round(100 * lab.loc["PSO_skin", 0.05], 1)
D["cal_any_call_labels"] = {str(a): int((lab[a].drop("PSO_skin") > 0).sum()) for a in alphas}
D["cal_any_call_archetype"] = {str(a): int((arc[a].drop("PSO_skin") > 0).sum()) for a in alphas}
D["cal_nfp_alpha05"] = nfp[0.05].to_dict()


def _poibin_tail(probs, k):
    """P(X >= k) for a sum of independent Bernoulli(probs)."""
    dist = np.array([1.0])
    for q in probs:
        dist = np.convolve(dist, [1 - q, q])
    return float(dist[k:].sum())


# Excess of datasets with any call over the BH bound, under either BH grouping:
# per dataset, P(any call) <= alpha; per donor, P(any call) <= min(1, alpha * donors).
don = T0.n_donors.drop("PSO_skin")
excess = {}
for src, P in [("labels", lab), ("archetype", arc)]:
    for a in alphas:
        obs = int((P[a].drop("PSO_skin") > 0).sum())
        pd_ = np.full(len(don), a)
        pn_ = np.minimum(1.0, a * don.values)
        excess[f"{src}_{a}"] = {"observed": obs,
                                "bound_mean_per_dataset": round(float(pd_.sum()), 2),
                                "p_per_dataset": float(f"{_poibin_tail(pd_, obs):.1e}"),
                                "bound_mean_per_donor": round(float(pn_.sum()), 2),
                                "p_per_donor": float(f"{_poibin_tail(pn_, obs):.1e}")}
D["cal_any_call_excess"] = excess
D["cal_exceed_nominal_labels"] = [d for d in lab.index if (lab.loc[d] > np.array(alphas)).any()]
D["cal_exceed_nominal_archetype"] = [d for d in arc.index if (arc.loc[d] > np.array(alphas)).any()]

cm = T10.pivot(index="dataset", columns="null_source", values="realized")
cmn = T10.groupby("dataset").n.first()
D["donor_pso"] = [round(100 * cm.loc["PSO_skin", "healthy"], 1), round(100 * cm.loc["PSO_skin", "matched"], 2)]
D["donor_matched_max"] = round(100 * cm["matched"].max(), 2)
D["donor_n_datasets"] = int(len(cm))
disc = {}
for d in cm.index:
    disc[d] = [round(100 * lab.loc[d, 0.05], 2), round(100 * cm.loc[d, "healthy"], 2)]
D["t5c_vs_t10_crosstissue_alpha05"] = disc

rows = []
for d in ORDER:
    k = int(nfp.loc[d, 0.05]); n = int(nlab[d])
    lo, hi = wilson(k, n)
    t10h = f2(100 * cm.loc[d, "healthy"]) if d in cm.index else "--"
    t10m = f2(100 * cm.loc[d, "matched"]) if d in cm.index else "--"
    rows.append(" & ".join([PRETTY[d], f1(tir[d]), f"{k}/{n:,}", f2(100 * k / n),
                            f"{100 * lo:.2f}--{100 * hi:.2f}", f2(100 * arc.loc[d, 0.05]),
                            t10h, t10m]) + r" \\")
    if d == "TS_Stomach":
        rows.append(r"\midrule")
write("table3_calibration", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Calibrated decision rule at nominal level $\alpha = 5\%$.} Tirosh: false-call rate of the two-criterion rule (as in Table~\ref{tab:panel}). Leave-one-tissue-out (LOTO): null fitted on the healthy datasets other than the query, conditioned on annotated cell type (false calls / cells, per-cell false-call rate and Wilson 95\% interval) or on label-free archetype. Donor-held-out: a separate evaluation on cells of donors with at least 200 cells, comparing a null from other tissues with a null from the other donors of the same dataset. Because of that restriction its denominator is smaller than the LOTO denominator for 7 of 13 datasets, so the two sets of columns do not always score the same cells; where they do, its cross-tissue column still fails to reproduce the LOTO column for three datasets (large intestine, liver, lung; see Limitations). --: dataset not eligible (fewer than two qualifying donors).}
\label{tab:cal}
\scriptsize
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}lrrrrrrr@{}}
\toprule
 & & \multicolumn{4}{c}{LOTO null (\%)} & \multicolumn{2}{c}{Donor-held-out (\%)} \\
\cmidrule(lr){3-6}\cmidrule(lr){7-8}
Dataset & Tirosh (\%) & Calls & Rate & 95\% CI & Archetype & Other tissues & Same dataset \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

rows = []
for src, P in [("Cell type", lab), ("Archetype", arc)]:
    for d in ORDER:
        rows.append(" & ".join([src if d == ORDER[0] else "", PRETTY[d]] +
                               [f"{100 * P.loc[d, a]:.2f}" for a in alphas]) + r" \\")
    key = "cal_any_call_labels" if src == "Cell type" else "cal_any_call_archetype"
    rows.append(" & ".join(["", r"\textit{Datasets with $\geq$1 call (of 13)}"] +
                           [f"\\textit{{{D[key][str(a)]}}}" for a in alphas]) + r" \\")
    if src == "Cell type":
        rows.append(r"\midrule")
write("tableS4_alpha", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Realized per-cell false-call rate (\%) of the calibrated rule across nominal levels.} Leave-one-tissue-out null conditioned on cell type or on label-free archetype. The last row of each block counts datasets, excluding psoriatic skin, with at least one false call. Under Benjamini--Hochberg control at level $\alpha$ with valid $p$-values, the expected count is at most $13\alpha$ if each dataset is one BH family, or $\sum_i \min(1, \alpha d_i) = """ + ", ".join(f"{excess[f'labels_{a}']['bound_mean_per_donor']:.2f}" for a in alphas) + r"""$ at the four levels if each of the $d_i$ donors of dataset $i$ is a family.}
\label{tab:s4}
\scriptsize
\begin{tabular}{@{}llrrrr@{}}
\toprule
Conditioning & Dataset & """ + " & ".join(f"$\\alpha={100 * a:g}\\%$" for a in alphas) + r""" \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

rows = [" & ".join([PRETTY[d], f2(per[d]), str(int(T0.loc[d, "n_donors"]))]) + r" \\" for d in ORDER]
write("tableS5_crosstissue", r"""\begin{table}[htbp]
\centering
\caption{\textbf{Reproducibility of the false-call difference profile.} Median Spearman correlation of each dataset's 10-Mb-bin difference profile (falsely called minus uncalled cells, two-criterion rule, reference-free mode) with those of the other 13 datasets. The Crohn's immune profile is anti-correlated with every other dataset (Spearman $\rho$ from """ + f"{D['cdimm_rho_range'][0]:.2f} to {D['cdimm_rho_range'][1]:.2f}" + r""").}
\label{tab:s5}
\small
\begin{tabular}{@{}lrr@{}}
\toprule
Dataset & Median $\rho$ & Donors \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

def _py(o):
    if isinstance(o, dict):
        return {str(k): _py(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_py(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    return o

(R / "derived_numbers.json").write_text(json.dumps(_py(D), indent=1))
print(json.dumps(_py(D), indent=1))
