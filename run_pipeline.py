"""
Run Part 1 end to end: load -> clean -> DE -> FDR -> shared genes -> Cytoscape export
                       -> sensitivity table -> hand-off list -> pathway enrichment / p53 axis.

    python run_pipeline.py --demo                      # fake data with planted answers (do this first)
    python run_pipeline.py                             # real data, primary analysis  -> results/
    python run_pipeline.py --cancer-model subtype      # adjusted for BRCA subtype    -> results_subtype/
    python run_pipeline.py --cancer-model stratum      # one subtype only             -> results_stratum/
    python run_pipeline.py --skip-enrichment           # skip step 9 (it needs internet for Enrichr libraries)
"""
import argparse
import config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true", help="use synthetic data in data_demo/")
    ap.add_argument("--cancer-model", choices=["none", "subtype", "stratum"], default="none")
    ap.add_argument("--skip-enrichment", action="store_true")
    args = ap.parse_args()

    config.CANCER_MODEL = args.cancer_model
    suffix = "" if args.cancer_model == "none" else f"_{args.cancer_model}"
    config.RESULTS_DIR = config.ROOT / f"results{suffix}"

    if args.demo:
        import make_synthetic_data
        config.DATA_DIR = config.ROOT / "data_demo"
        config.RESULTS_DIR = config.ROOT / f"results_demo{suffix}"
        config.ENRICHMENT_LIBRARIES = ["demo_sets.gmt"]       # offline gene sets for the demo
        config.P53_LIBRARY = "demo_sets.gmt"
        config.P53_CORE_TARGETS = ["GENE0300", "GENE0301", "GENE0302"]
        config.WATCHLIST = {"demo gene": ["GENE0001"], "missing gene": ["NOT_A_GENE"]}
        make_synthetic_data.build(config.DATA_DIR)
        print("== DEMO MODE (synthetic data) ==")
    print(f"== Cancer model: {args.cancer_model} ==")

    import step1_load_ad, step2_load_cancer, step3_differential_expression
    import step4_call_significant, step5_shared_genes, step6_export_cytoscape
    import step7_sensitivity, step8_prioritize_for_part2, step9_pathway_enrichment

    ad = step1_load_ad.run()
    cancer = step2_load_cancer.run()
    de_ad, de_cancer = step3_differential_expression.run(ad, cancer)
    ad_ann, ca_ann = step4_call_significant.run(de_ad, de_cancer)
    shared, _ = step5_shared_genes.run(ad_ann, ca_ann)
    step6_export_cytoscape.run(shared)
    step7_sensitivity.run(ad_ann, ca_ann)
    tiered, _, _ = step8_prioritize_for_part2.run(ad_ann, ca_ann, shared)
    if not args.skip_enrichment:
        try:
            step9_pathway_enrichment.run(ad_ann, ca_ann, tiered)
        except Exception as err:            # most likely: no internet for Enrichr libraries
            print(f"[step9] SKIPPED - {type(err).__name__}: {err}\n"
                  "        (needs internet + `pip install gseapy`, or put .gmt files in "
                  "config.ENRICHMENT_LIBRARIES; rerun with --skip-enrichment to silence this)")
    print(f"\nDone. Outputs in {config.RESULTS_DIR}")


if __name__ == "__main__":
    main()
