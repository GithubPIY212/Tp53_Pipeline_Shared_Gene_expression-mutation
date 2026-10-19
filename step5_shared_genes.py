"""
Step 5 - Shared genes: the end point of Part 1.

1. Intersect significant genes from the two arms.
2. Is the overlap bigger than chance? Hypergeometric test, with the background = genes that
   were tested in BOTH datasets (using all human genes would exaggerate significance).
3. Do shared genes move the same way in both diseases? (binomial test vs 50/50)

Output: results/shared_genes.tsv, results/overlap_summary.tsv
"""
import pandas as pd
from scipy import stats
import config


def overlap_statistics(ad, ca):
    background = ad.index.intersection(ca.index)
    sig_ad = set(ad.index[ad.significant]) & set(background)
    sig_ca = set(ca.index[ca.significant]) & set(background)
    shared = sorted(sig_ad & sig_ca)

    M, n, N, k = len(background), len(sig_ad), len(sig_ca), len(shared)
    p_overlap = stats.hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
    expected = n * N / M if M else float("nan")

    table = pd.DataFrame({
        "gene": shared,
        "log2FC_AD": ad.loc[shared, "log2FC"].values,
        "fdr_AD": ad.loc[shared, "fdr"].values,
        "log2FC_cancer": ca.loc[shared, "log2FC"].values,
        "fdr_cancer": ca.loc[shared, "fdr"].values,
    })
    table["same_direction"] = (table.log2FC_AD * table.log2FC_cancer) > 0
    table["pattern"] = ("AD_" + ad.loc[shared, "direction"].values +
                        "__cancer_" + ca.loc[shared, "direction"].values)
    n_same = int(table.same_direction.sum())
    p_dir = stats.binomtest(n_same, k, 0.5).pvalue if k > 0 else float("nan")

    summary = pd.DataFrame([{
        "genes_tested_in_both": M, "significant_AD": n, "significant_cancer": N,
        "shared": k, "expected_by_chance": round(expected, 2),
        "enrichment_fold": round(k / expected, 2) if expected else float("nan"),
        "hypergeometric_p": p_overlap,
        "shared_same_direction": n_same, "direction_binomial_p": p_dir,
    }])
    return table.sort_values("fdr_cancer").reset_index(drop=True), summary


def run(ad=None, ca=None):
    if ad is None:
        ad = pd.read_csv(config.results_path("de_ad_annotated.tsv"), sep="\t", index_col=0)
    if ca is None:
        ca = pd.read_csv(config.results_path("de_cancer_annotated.tsv"), sep="\t", index_col=0)
    table, summary = overlap_statistics(ad, ca)
    table.to_csv(config.results_path("shared_genes.tsv"), sep="\t", index=False)
    summary.to_csv(config.results_path("overlap_summary.tsv"), sep="\t", index=False)
    s = summary.iloc[0]
    print(f"[step5] Shared genes: {int(s.shared)} (expected by chance {s.expected_by_chance}); "
          f"hypergeometric p = {s.hypergeometric_p:.3g}; same direction {int(s.shared_same_direction)}/{int(s.shared)}")
    return table, summary


if __name__ == "__main__":
    run()
