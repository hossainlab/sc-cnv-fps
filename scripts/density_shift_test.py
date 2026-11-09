"""Gene density vs apparent copy-number shift, tested against a circular-shift null.

Reads the panel-averaged binned difference profile (results/T6_density_vs_bias.csv)
and tests its Spearman correlation with gene count per 10-Mb bin. Adjacent bins share
smoothed windows, so the asymptotic P-value is not calibrated; instead the profile is
shifted circularly within each chromosome by an independent random offset, leaving
gene counts in place (same scheme as attribution.locus_enrichment). This keeps the
autocorrelation of the profile and the gene-count distribution while breaking their
positional correspondence.

Usage:  python scripts/density_shift_test.py [out.csv]   (default results/T6b_density_shift_test.csv)
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata


def shift_test(df, n_shifts=10_000, seed=0):
    df = df.sort_values(["chrom", "bin"]).reset_index(drop=True)
    # chromosomes are visited in string order (1, 10, 11, ..., 19, 2, 20, ...), as in the
    # original run; the order only affects which random offsets go to which chromosome
    chroms = sorted(pd.unique(df["chrom"]).tolist(), key=str)
    rho, p_asym = spearmanr(df["mean_diff"], df["n_genes"])
    groups = [np.flatnonzero(df["chrom"].values == c) for c in chroms]
    y = df["mean_diff"].values
    g_rank = rankdata(df["n_genes"].values)
    rng = np.random.default_rng(seed)
    null = np.empty(n_shifts)
    for i in range(n_shifts):
        ys = np.empty_like(y)
        for idx in groups:
            ys[idx] = np.roll(y[idx], rng.integers(0, len(idx)))
        null[i] = np.corrcoef(rankdata(ys), g_rank)[0, 1]
    # one-sided: association at least as negative as observed
    p_shift = (1 + int((null <= rho).sum())) / (1 + n_shifts)
    return {
        "test": "gene_density_vs_apparent_shift", "n_bins": len(df),
        "spearman_rho": round(float(rho), 4), "p_asymptotic": float(p_asym),
        "n_shifts": n_shifts, "p_circular_shift_onesided": p_shift,
        "shift_null_rho_mean": round(float(null.mean()), 4),
        "shift_null_rho_sd": round(float(null.std()), 4),
        "shift_null_rho_q2.5": round(float(np.quantile(null, 0.025)), 4),
        "shift_null_rho_q97.5": round(float(np.quantile(null, 0.975)), 4),
        "shift_null_rho_min": round(float(null.min()), 4),
        "z_vs_shift_null": round(float((null.mean() - rho) / null.std()), 2),
        "seed": seed,
    }


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "results/T6b_density_shift_test.csv"
    res = shift_test(pd.read_csv("results/T6_density_vs_bias.csv"))
    pd.DataFrame([res]).to_csv(out, index=False)
    print(res)
