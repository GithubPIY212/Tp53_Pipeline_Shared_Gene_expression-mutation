"""
Step 9 - Pathway enrichment of the shared genes + a direct look at the p53 pathway.

A) Over-representation: are the shared genes enriched for GO / KEGG / Reactome / Hallmark terms?
   The background is the genes tested in BOTH datasets (not the whole genome), otherwise
   terms look enriched just because of which genes the platforms measure.
B) p53 axis: do p53 pathway genes move in OPPOSITE directions in AD and TP53-altered cancer
   (the pattern seen for CDKN1A and MDM2)? Judged against the genome-wide rate of opposite
   direction, and with a table of core p53 targets. HYPOTHESIS-GENERATING, not confirmatory.

Gene sets come from Enrichr library names (needs internet + `pip install gseapy`) or from local
.gmt files listed in config.ENRICHMENT_LIBRARIES.

Output: results/enrichment_<list>.tsv, results/p53_pathway_direction.tsv, results/p53_core_targets.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import config


# ---------- gene sets ----------
def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_gene_sets(source):
    p = Path(str(source))
    if not p.exists() and config.data_path(str(source)).exists():
        p = config.data_path(str(source))
    if str(source).lower().endswith(".gmt"):
        return read_gmt(p)
    import gseapy as gp                                   # imported here so offline use works
    return gp.get_library(name=str(source), organism="Human")


# ---------- A) over-representation ----------
def over_representation(genes, gene_sets, background):
    bg = set(background)
    genes = set(genes) & bg
    rows = []
    for term, members in gene_sets.items():
        in_bg = set(members) & bg
        k = len(genes & in_bg)
        if len(in_bg) < config.ENRICH_MIN_SET_SIZE or k < config.ENRICH_MIN_OVERLAP:
            continue
        p = stats.hypergeom.sf(k - 1, len(bg), len(in_bg), len(genes))
        rows.append({"term": term, "overlap": k, "term_size_in_background": len(in_bg),
                     "query_size": len(genes), "fold_enrichment": round(k / (len(in_bg) * len(genes) / len(bg)), 2),
                     "p_value": p, "genes": ";".join(sorted(genes & in_bg))})
    if not rows:
        return pd.DataFrame(columns=["term", "overlap", "term_size_in_background", "query_size",
                                     "fold_enrichment", "p_value", "genes", "fdr"])
    df = pd.DataFrame(rows).sort_values("p_value")
    df["fdr"] = multipletests(df.p_value, method="fdr_bh")[1]
    return df.reset_index(drop=True)


# ---------- B) p53 axis ----------
def p53_direction(ad, ca, p53_genes):
    common = ad.index.intersection(ca.index)
    genes = [g for g in p53_genes if g in common]
    if len(genes) < 5:
        raise ValueError(f"Only {len(genes)} p53 pathway genes found in both datasets")
    a, c = ad.loc[genes, "log2FC"], ca.loc[genes, "log2FC"]
    opposite = int((np.sign(a) != np.sign(c)).sum())
    baseline = float((np.sign(ad.loc[common, "log2FC"]) != np.sign(ca.loc[common, "log2FC"])).mean())
    return pd.DataFrame([{
        "p53_genes_tested_in_both": len(genes),
        "mean_log2FC_AD": a.mean(), "fraction_up_AD": float((a > 0).mean()),
        "wilcoxon_p_AD": stats.wilcoxon(a).pvalue,
        "mean_log2FC_cancer": c.mean(), "fraction_up_cancer": float((c > 0).mean()),
        "wilcoxon_p_cancer": stats.wilcoxon(c).pvalue,
        "opposite_direction": opposite, "fraction_opposite": opposite / len(genes),
        "genomewide_fraction_opposite": baseline,
        "opposite_vs_genomewide_p": stats.binomtest(opposite, len(genes), baseline).pvalue,
    }])


def core_targets(ad, ca):
    rows = []
    for g in config.P53_CORE_TARGETS:
        rows.append({"gene": g,
                     "log2FC_AD": ad.log2FC.get(g, np.nan), "fdr_AD": ad.fdr.get(g, np.nan),
                     "log2FC_cancer": ca.log2FC.get(g, np.nan), "fdr_cancer": ca.fdr.get(g, np.nan)})
    return pd.DataFrame(rows)


def run(ad=None, ca=None, tiered=None):
    if ad is None:
        ad = pd.read_csv(config.results_path("de_ad_annotated.tsv"), sep="\t", index_col=0)
    if ca is None:
        ca = pd.read_csv(config.results_path("de_cancer_annotated.tsv"), sep="\t", index_col=0)
    if tiered is None:
        tiered = pd.read_csv(config.results_path("shared_genes_tiered.tsv"), sep="\t")
    background = list(ad.index.intersection(ca.index))

    core = core_targets(ad, ca)
    core.to_csv(config.results_path("p53_core_targets.tsv"), sep="\t", index=False)
    print("[step9] Core p53 targets (log2FC AD vs cancer):")
    print(core.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    lists = {"shared_all": tiered.gene.tolist(),
             "priority_tierA": tiered.loc[tiered.tier == "A_strict_concordant", "gene"].tolist()}
    libs = {}
    for lib in config.ENRICHMENT_LIBRARIES:
        libs[str(lib)] = load_gene_sets(lib)
    for name, genes in lists.items():
        if len(genes) < 5:
            print(f"[step9] list '{name}' has {len(genes)} genes - skipped")
            continue
        parts = []
        for lib, sets in libs.items():
            res = over_representation(genes, sets, background)
            res.insert(0, "library", Path(lib).stem if lib.lower().endswith(".gmt") else lib)
            parts.append(res)
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(config.results_path(f"enrichment_{name}.tsv"), sep="\t", index=False)
        n_sig = int((out.fdr < 0.05).sum()) if len(out) else 0
        print(f"[step9] enrichment '{name}' ({len(genes)} genes): {n_sig} terms with FDR < 0.05")

    p53_lib = libs.get(config.P53_LIBRARY) or load_gene_sets(config.P53_LIBRARY)
    key = next((k for k in p53_lib if config.P53_SET_KEYWORD in k.lower()), None)
    if key is None:
        print(f"[step9] No set matching '{config.P53_SET_KEYWORD}' in {config.P53_LIBRARY}; p53 axis test skipped")
        return None
    direction = p53_direction(ad, ca, p53_lib[key])
    direction.insert(0, "gene_set", key)
    direction.to_csv(config.results_path("p53_pathway_direction.tsv"), sep="\t", index=False)
    d = direction.iloc[0]
    print(f"[step9] p53 pathway ('{key}', {int(d.p53_genes_tested_in_both)} genes): "
          f"opposite direction in {int(d.opposite_direction)} (genome-wide rate {d.genomewide_fraction_opposite:.2f}), "
          f"p = {d.opposite_vs_genomewide_p:.3g}")
    return direction


if __name__ == "__main__":
    run()
