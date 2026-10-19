"""
Step 8 - Turn the shared-gene list into a short, ranked hand-off list for the next part
(CRISPRi/CRISPRa target choice), plus a status table for the genes the team is watching.

Tiers (all shared genes get one):
  A_strict_concordant : FDR < PRIORITY_FDR in BOTH arms, |log2FC| above the strict cut-offs in
                        BOTH arms, and the SAME direction in both diseases
  B_concordant        : same direction, passes only the primary cut-offs
  C_discordant        : opposite direction in the two diseases (kept: this is where the p53
                        'active in AD / lost in cancer' pattern would show up)
Rank inside tier A = geometric mean of the two |log2FC| values (large in both = high rank).

RUN THIS ON THE SUBTYPE-ADJUSTED RESULTS if you can: unadjusted cancer fold changes are dominated
by basal-like vs other-subtype biology.

Output: results/part2_priority_genes.tsv   (tier A, ranked; top_n_pick marks the first N)
        results/shared_genes_tiered.tsv    (every shared gene with its tier)
        results/part2_watchlist.tsv        (status of the watched genes in both arms)
"""
import numpy as np
import pandas as pd
import config


def assign_tiers(shared):
    s = shared.copy()
    strict = ((s.fdr_AD < config.PRIORITY_FDR) & (s.fdr_cancer < config.PRIORITY_FDR) &
              (s.log2FC_AD.abs() >= config.PRIORITY_AD_MIN_ABS_LOG2FC) &
              (s.log2FC_cancer.abs() >= config.PRIORITY_CANCER_MIN_ABS_LOG2FC))
    s["tier"] = np.where(~s.same_direction, "C_discordant",
                         np.where(strict, "A_strict_concordant", "B_concordant"))
    s["score"] = np.sqrt(s.log2FC_AD.abs() * s.log2FC_cancer.abs())
    return s.sort_values(["tier", "score"], ascending=[True, False]).reset_index(drop=True)


def watchlist_status(ad, ca, tiered):
    rows = []
    shared_tier = tiered.set_index("gene")["tier"].to_dict()
    for label, symbols in config.WATCHLIST.items():
        sym = next((s for s in symbols if s in ad.index or s in ca.index), None)
        row = {"watched_gene": label, "symbol_found": sym or "not in either dataset"}
        for arm, df in (("AD", ad), ("cancer", ca)):
            if sym is not None and sym in df.index:
                row[f"log2FC_{arm}"], row[f"fdr_{arm}"] = df.loc[sym, "log2FC"], df.loc[sym, "fdr"]
            else:
                row[f"log2FC_{arm}"], row[f"fdr_{arm}"] = np.nan, np.nan
        row["shared_tier"] = shared_tier.get(sym, "not shared")
        rows.append(row)
    return pd.DataFrame(rows)


def run(ad=None, ca=None, shared=None):
    if ad is None:
        ad = pd.read_csv(config.results_path("de_ad_annotated.tsv"), sep="\t", index_col=0)
    if ca is None:
        ca = pd.read_csv(config.results_path("de_cancer_annotated.tsv"), sep="\t", index_col=0)
    if shared is None:
        shared = pd.read_csv(config.results_path("shared_genes.tsv"), sep="\t")
    tiered = assign_tiers(shared)
    tiered.to_csv(config.results_path("shared_genes_tiered.tsv"), sep="\t", index=False)
    top = tiered[tiered.tier == "A_strict_concordant"].reset_index(drop=True)
    top.insert(0, "rank", range(1, len(top) + 1))
    top["top_n_pick"] = top["rank"] <= config.PRIORITY_TOP_N
    top.to_csv(config.results_path("part2_priority_genes.tsv"), sep="\t", index=False)
    watch = watchlist_status(ad, ca, tiered)
    watch.to_csv(config.results_path("part2_watchlist.tsv"), sep="\t", index=False)
    counts = tiered.tier.value_counts().to_dict()
    print(f"[step8] Tiers: {counts}; hand-off list = {len(top)} genes (top {config.PRIORITY_TOP_N} flagged)")
    return tiered, top, watch


if __name__ == "__main__":
    run()
