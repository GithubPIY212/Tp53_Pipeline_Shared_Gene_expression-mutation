"""
Step 7 - Sensitivity analysis: does the overlap result depend on the cut-offs?

The primary analysis (settings in config.py) is fixed BEFORE looking at results. This step
re-runs the overlap test over a grid of stricter cut-offs and shows all of them side by side,
so the write-up can say honestly "the primary result is X; it changes/does not change if...".
Rows other than the primary one are post-hoc: report them as sensitivity analysis, never as
the headline.

Output: results/sensitivity_table.tsv
"""
import numpy as np
import pandas as pd
from scipy import stats
import config


def sensitivity_table(ad, ca, grid):
    common = ad.index.intersection(ca.index)
    M = len(common)
    primary = (config.AD_MIN_ABS_LOG2FC, config.CANCER_MIN_ABS_LOG2FC, config.FDR_THRESHOLD)
    rows = []
    for ad_cut, ca_cut, fdr in grid:
        A = set(ad.index[(ad.fdr < fdr) & (ad.log2FC.abs() >= ad_cut)]) & set(common)
        C = set(ca.index[(ca.fdr < fdr) & (ca.log2FC.abs() >= ca_cut)]) & set(common)
        shared = sorted(A & C)
        k, expected = len(shared), len(A) * len(C) / M
        p_overlap = stats.hypergeom.sf(k - 1, M, len(A), len(C)) if k else 1.0
        same = int((np.sign(ad.loc[shared, "log2FC"]) == np.sign(ca.loc[shared, "log2FC"])).sum()) if k else 0
        rows.append({
            "ad_min_abs_log2FC": ad_cut, "cancer_min_abs_log2FC": ca_cut, "fdr": fdr,
            "significant_AD": len(A), "significant_cancer": len(C), "shared": k,
            "expected_by_chance": round(expected, 1),
            "enrichment_fold": round(k / expected, 2) if expected else np.nan,
            "overlap_p": p_overlap, "same_direction": same,
            "fraction_same_direction": round(same / k, 3) if k else np.nan,
            "direction_p": stats.binomtest(same, k, 0.5).pvalue if k else np.nan,
            "analysis": "PRIMARY (pre-set)" if (ad_cut, ca_cut, fdr) == primary else "post-hoc sensitivity",
        })
    return pd.DataFrame(rows)


def run(ad=None, ca=None):
    if ad is None:
        ad = pd.read_csv(config.results_path("de_ad_annotated.tsv"), sep="\t", index_col=0)
    if ca is None:
        ca = pd.read_csv(config.results_path("de_cancer_annotated.tsv"), sep="\t", index_col=0)
    grid = list(config.SENSITIVITY_GRID)
    primary = (config.AD_MIN_ABS_LOG2FC, config.CANCER_MIN_ABS_LOG2FC, config.FDR_THRESHOLD)
    if primary not in grid:
        grid.insert(0, primary)
    table = sensitivity_table(ad, ca, grid)
    table.to_csv(config.results_path("sensitivity_table.tsv"), sep="\t", index=False)
    show = table[["ad_min_abs_log2FC", "cancer_min_abs_log2FC", "fdr", "shared", "expected_by_chance",
                  "enrichment_fold", "overlap_p", "fraction_same_direction", "analysis"]]
    print("[step7] Sensitivity of the overlap result to the cut-offs:")
    print(show.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    return table


if __name__ == "__main__":
    run()
