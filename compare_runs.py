"""
Compare the primary run with the subtype-adjusted run: how much of the cancer signal (and of the
shared-gene list) was really breast-cancer subtype?

    python compare_runs.py results results_subtype
"""
import sys
from pathlib import Path
import pandas as pd


def compare(dir_a, dir_b):
    ca_a = pd.read_csv(Path(dir_a) / "de_cancer_annotated.tsv", sep="\t", index_col=0)
    ca_b = pd.read_csv(Path(dir_b) / "de_cancer_annotated.tsv", sep="\t", index_col=0)
    sh_a = set(pd.read_csv(Path(dir_a) / "shared_genes.tsv", sep="\t").gene)
    sh_b = set(pd.read_csv(Path(dir_b) / "shared_genes.tsv", sep="\t").gene)
    sig_a, sig_b = set(ca_a.index[ca_a.significant]), set(ca_b.index[ca_b.significant])
    both = sorted(sh_a & sh_b)
    return pd.DataFrame([{
        "cancer_significant_primary": len(sig_a), "cancer_significant_adjusted": len(sig_b),
        "cancer_significant_in_both": len(sig_a & sig_b),
        "fraction_of_primary_kept": round(len(sig_a & sig_b) / len(sig_a), 3) if sig_a else None,
        "shared_primary": len(sh_a), "shared_adjusted": len(sh_b), "shared_in_both": len(both),
    }]), both


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python compare_runs.py results results_subtype")
    table, genes = compare(sys.argv[1], sys.argv[2])
    print(table.T.to_string(header=False))
    pd.Series(genes, name="gene").to_csv(Path(sys.argv[2]) / "shared_genes_robust_to_subtype.txt", index=False)
    print(f"Saved {len(genes)} genes shared in BOTH runs to {sys.argv[2]}/shared_genes_robust_to_subtype.txt")
