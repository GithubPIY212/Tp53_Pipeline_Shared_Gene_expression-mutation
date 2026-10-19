# TP53 / LINE-1 pipeline - Part 1 (up to the shared gene list and its follow-up)

Question: do genes that change with TP53 alteration in breast cancer (TCGA-BRCA) overlap with
genes that change in Alzheimer's brain (GSE33000) more than chance predicts?

## Quick start
    pip install -r requirements.txt
    pytest tests                                   # 18 tests should pass
    python run_pipeline.py --demo                  # fake data with planted answers
    python run_pipeline.py                         # real data, PRIMARY analysis      -> results/
    python run_pipeline.py --cancer-model subtype  # adjusted for BRCA subtype        -> results_subtype/
    python compare_runs.py results results_subtype # how much of the signal was subtype?

## Real data in `data/`
| File | From |
|---|---|
| `GSE33000_expression.tsv.gz`, `GSE33000_metadata.tsv` | GEO (build metadata with `make_ad_metadata.py`) |
| `data_mrna_seq_v2_rsem.txt`, `data_mutations.txt` | cBioPortal BRCA PanCancer Atlas tarball |
| `data_clinical_patient.txt` (and/or `data_clinical_sample.txt`) | same tarball - needed ONLY for `--cancer-model subtype/stratum` |

## Steps
1. `step1_load_ad.py` - clean aliases, average duplicate probes, keep AD vs non-demented
2. `step2_load_cancer.py` - log2(RSEM+1), drop near-silent genes, TP53 mutated vs not
3. `step3_differential_expression.py` - Welch t-test (or subtype-adjusted model, `subtype_adjustment.py`); AD log10 -> log2FC
4. `step4_call_significant.py` - Benjamini-Hochberg FDR + effect-size cut-off
5. `step5_shared_genes.py` - intersection, hypergeometric overlap test, direction check
6. `step6_export_cytoscape.py` - node table + STRING gene list
7. `step7_sensitivity.py` - the overlap test repeated over a grid of cut-offs (`sensitivity_table.tsv`)
8. `step8_prioritize_for_part2.py` - tiered, ranked hand-off list + status of watched genes
9. `step9_pathway_enrichment.py` - GO/KEGG/Reactome/Hallmark enrichment + p53-axis direction test (needs internet for Enrichr names, or local .gmt files)

## Reading the results honestly
- The PRIMARY analysis is the one set in `config.py` BEFORE looking at results. Step 7 rows marked
  "post-hoc" are sensitivity analysis, not headline results.
- Hand the next part the list from `results_subtype/` (step 8), not from the unadjusted run.
- Top-ranked shared genes can reflect tissue composition (immune or neuronal genes), not p53 biology: check them.
- p53-axis test (step 9) is hypothesis-generating.
- Caveats to state: TP53 status = mutation calls only; microarray vs RNA-seq; brain vs breast;
  LINE-1 not measured directly; AD arm not adjusted for age/sex/batch.
