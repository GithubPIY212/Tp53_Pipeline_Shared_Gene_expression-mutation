"""
Step 6 - Export for Cytoscape / STRING.

* cytoscape_node_table.tsv : import in Cytoscape (File > Import > Table from File) after building
  a network, or use it alone to style nodes. Columns carry fold changes and direction.
* string_gene_list.txt     : paste into STRING (or Cytoscape's stringApp) to fetch interactions.
"""
import pandas as pd
import config


def run(shared=None):
    if shared is None:
        shared = pd.read_csv(config.results_path("shared_genes.tsv"), sep="\t")
    nodes = shared.rename(columns={"gene": "name"})
    nodes["mean_log2FC"] = (nodes.log2FC_AD + nodes.log2FC_cancer) / 2
    nodes.to_csv(config.results_path("cytoscape_node_table.tsv"), sep="\t", index=False)
    with open(config.results_path("string_gene_list.txt"), "w") as fh:
        fh.write("\n".join(shared["gene"].astype(str)) + ("\n" if len(shared) else ""))
    print(f"[step6] Exported {len(shared)} shared genes for Cytoscape/STRING")
    return nodes


if __name__ == "__main__":
    run()
