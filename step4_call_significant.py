"""
Step 4 - Multiple-testing correction and significance calls.

Testing ~20,000 genes at p<0.05 would give ~1,000 false hits by luck alone. Benjamini-Hochberg
FDR controls the expected share of false discoveries among the genes we call significant.
A gene is 'significant' only if FDR < threshold AND |log2FC| >= the arm's effect-size cut-off.

Output: results/de_ad_annotated.tsv, results/de_cancer_annotated.tsv,
        results/sig_ad.tsv, results/sig_cancer.tsv
"""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import config


def call_significant(de, min_abs_log2fc):
    de = de.copy()
    ok = de["p_value"].notna()
    de["fdr"] = np.nan
    de.loc[ok, "fdr"] = multipletests(de.loc[ok, "p_value"], method="fdr_bh")[1]
    de["significant"] = (de["fdr"] < config.FDR_THRESHOLD) & (de["log2FC"].abs() >= min_abs_log2fc)
    de["direction"] = np.where(de["log2FC"] > 0, "up", "down")
    return de


def run(de_ad=None, de_cancer=None):
    if de_ad is None:
        de_ad = pd.read_csv(config.results_path("de_ad.tsv"), sep="\t", index_col=0)
    if de_cancer is None:
        de_cancer = pd.read_csv(config.results_path("de_cancer.tsv"), sep="\t", index_col=0)

    ad = call_significant(de_ad, config.AD_MIN_ABS_LOG2FC)
    ca = call_significant(de_cancer, config.CANCER_MIN_ABS_LOG2FC)
    ad.to_csv(config.results_path("de_ad_annotated.tsv"), sep="\t")
    ca.to_csv(config.results_path("de_cancer_annotated.tsv"), sep="\t")
    ad[ad.significant].sort_values("fdr").to_csv(config.results_path("sig_ad.tsv"), sep="\t")
    ca[ca.significant].sort_values("fdr").to_csv(config.results_path("sig_cancer.tsv"), sep="\t")
    print(f"[step4] Significant genes: AD {int(ad.significant.sum())}, cancer {int(ca.significant.sum())}")
    return ad, ca


if __name__ == "__main__":
    run()
