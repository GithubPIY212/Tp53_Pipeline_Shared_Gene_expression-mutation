"""
Step 3 - Differential expression (case vs control) for every gene, in both arms.

Welch's t-test (does not assume equal variance) per gene.

FOLD CHANGE - the part that trips people up:
  * Cancer values are log2(RSEM+1)  -> log2FC = mean(case) - mean(control).
  * AD values are ALREADY log10(sample/reference). A difference of two means is therefore a
    log10 fold change; dividing by log10(2) converts it to log2 so both arms share one scale.
  Doing mean(case)-mean(control) on raw ratios, or treating the AD numbers as raw intensities,
  would give wrong effect sizes. P-values are unaffected by this rescaling.

Output: results/de_ad.tsv, results/de_cancer.tsv (all genes, no filtering yet)
"""
import numpy as np
import pandas as pd
from scipy import stats
import config
import subtype_adjustment as sa

LOG10_2 = np.log10(2)


def differential_expression(expr, groups, case, control, scale_to_log2=1.0):
    """expr: genes x samples. scale_to_log2: multiply the mean difference by this to get log2FC."""
    case_ids = groups.index[groups == case]
    ctrl_ids = groups.index[groups == control]
    if min(len(case_ids), len(ctrl_ids)) < config.MIN_SAMPLES_PER_GROUP:
        raise ValueError(f"Too few samples: {len(case_ids)} case / {len(ctrl_ids)} control")
    a = expr[case_ids].to_numpy(float)
    b = expr[ctrl_ids].to_numpy(float)
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    mean_a, mean_b = np.nanmean(a, axis=1), np.nanmean(b, axis=1)
    return pd.DataFrame({
        "gene": expr.index,
        "mean_case": mean_a,
        "mean_control": mean_b,
        "log2FC": (mean_a - mean_b) * scale_to_log2,
        "t_stat": t,
        "p_value": p,
    }).set_index("gene")


def _cancer_de(expr, groups):
    """Cancer arm DE under the model chosen in config.CANCER_MODEL (none / subtype / stratum)."""
    model = config.CANCER_MODEL
    case, ctrl = config.CANCER_CASE_LABEL, config.CANCER_CONTROL_LABEL
    if model == "none":
        return differential_expression(expr, groups, case, ctrl)
    subtypes = sa.load_subtypes(list(expr.columns))
    subtypes.to_frame().to_csv(config.results_path("cancer_subtype.tsv"), sep="\t")
    print(pd.crosstab(subtypes.fillna("unknown"), groups).to_string())   # sanity check: TP53 vs subtype
    if model == "subtype":
        return sa.adjusted_de(expr, groups, subtypes, case, ctrl)
    if model == "stratum":
        keep = list(expr.columns[(subtypes == config.CANCER_STRATUM).to_numpy()])
        print(f"[stratum] only {config.CANCER_STRATUM}: {len(keep)} samples")
        return differential_expression(expr[keep], groups[keep], case, ctrl)
    raise ValueError(f"Unknown CANCER_MODEL '{model}' (use none / subtype / stratum)")


def run(ad=None, cancer=None):
    if ad is None:
        ad = (pd.read_csv(config.results_path("ad_expression_clean.tsv"), sep="\t", index_col=0),
              pd.read_csv(config.results_path("ad_groups.tsv"), sep="\t", index_col=0)["group"])
    if cancer is None:
        cancer = (pd.read_csv(config.results_path("cancer_expression_clean.tsv"), sep="\t", index_col=0),
                  pd.read_csv(config.results_path("cancer_groups.tsv"), sep="\t", index_col=0)["group"])

    de_ad = differential_expression(ad[0], ad[1], config.AD_CASE_LABEL, config.AD_CONTROL_LABEL,
                                    scale_to_log2=1.0 / LOG10_2)
    de_cancer = _cancer_de(cancer[0], cancer[1])
    de_ad.to_csv(config.results_path("de_ad.tsv"), sep="\t")
    de_cancer.to_csv(config.results_path("de_cancer.tsv"), sep="\t")
    print(f"[step3] DE done: AD {len(de_ad)} genes, cancer {len(de_cancer)} genes")
    return de_ad, de_cancer


if __name__ == "__main__":
    run()
