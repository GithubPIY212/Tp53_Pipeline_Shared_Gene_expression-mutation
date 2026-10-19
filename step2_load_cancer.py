"""
Step 2 - Load TCGA-BRCA (cBioPortal), define TP53 status, log-transform.

TP53 status comes from the mutation file: any sample with a TP53 mutation = altered.
CAVEAT (state it in your write-up): copy-number and structural TP53 changes are not used here.

Output: results/cancer_expression_clean.tsv  (genes x samples, log2(RSEM+1))
        results/cancer_groups.tsv
"""
import numpy as np
import pandas as pd
import config


def _sample_key(s):
    # TCGA-XX-XXXX-01 (first 15 chars) is the shared sample ID across cBioPortal files.
    return str(s)[:15]


def run():
    raw = pd.read_csv(config.data_path(config.CANCER_EXPRESSION_FILE), sep="\t", comment="#")
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in raw.columns else raw.columns[0]
    raw = raw.drop(columns=[c for c in ("Entrez_Gene_Id",) if c in raw.columns])
    raw = raw.dropna(subset=[gene_col]).set_index(gene_col)
    raw = raw.groupby(level=0).mean()
    raw.columns = [_sample_key(c) for c in raw.columns]
    raw = raw.loc[:, ~raw.columns.duplicated()]

    # RSEM counts are heavily skewed; log2(x+1) makes them roughly bell-shaped for t-tests.
    expr = np.log2(raw.astype(float) + 1)
    expr = expr[expr.mean(axis=1) >= config.CANCER_MIN_MEAN_LOG2]

    mut = pd.read_csv(config.data_path(config.CANCER_MUTATION_FILE), sep="\t", comment="#",
                      usecols=["Hugo_Symbol", "Tumor_Sample_Barcode"])
    altered = {_sample_key(s) for s in mut.loc[mut.Hugo_Symbol == config.TP53_GENE, "Tumor_Sample_Barcode"]}
    groups = pd.Series(
        [config.CANCER_CASE_LABEL if s in altered else config.CANCER_CONTROL_LABEL for s in expr.columns],
        index=expr.columns, name="group")

    expr.to_csv(config.results_path("cancer_expression_clean.tsv"), sep="\t")
    groups.to_csv(config.results_path("cancer_groups.tsv"), sep="\t")
    print(f"[step2] Cancer: {expr.shape[0]} genes x {expr.shape[1]} samples | "
          f"{(groups == config.CANCER_CASE_LABEL).sum()} TP53-altered vs "
          f"{(groups == config.CANCER_CONTROL_LABEL).sum()} unaltered")
    return expr, groups


if __name__ == "__main__":
    run()
