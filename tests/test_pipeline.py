import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
import make_synthetic_data as synth
import step1_load_ad, step2_load_cancer, step3_differential_expression as s3
import step4_call_significant as s4, step5_shared_genes as s5, step6_export_cytoscape as s6
import step7_sensitivity as s7, step8_prioritize_for_part2 as s8, step9_pathway_enrichment as s9
import subtype_adjustment as sa


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    d = tmp_path_factory.mktemp("demo")
    config.DATA_DIR, config.RESULTS_DIR = d / "data", d / "results"
    truth = synth.build(config.DATA_DIR)
    ad, ca = step1_load_ad.run(), step2_load_cancer.run()
    de_ad, de_ca = s3.run(ad, ca)
    ad_a, ca_a = s4.run(de_ad, de_ca)
    shared, summary = s5.run(ad_a, ca_a)
    return dict(truth=truth, ad=ad, ca=ca, ad_a=ad_a, ca_a=ca_a, shared=shared, summary=summary)


def test_huntingtons_excluded(results):
    groups = results["ad"][1]
    assert set(groups) == {config.AD_CASE_LABEL, config.AD_CONTROL_LABEL}
    assert len(groups) == 110


def test_alias_symbols_cleaned(results):
    assert not any("///" in g for g in results["ad"][0].index)


def test_tp53_groups_counted(results):
    assert (results["ca"][1] == config.CANCER_CASE_LABEL).sum() == 90


def test_low_expression_genes_filtered(results):
    assert "GENE2950" not in results["ca"][0].index


def test_log10_to_log2_fold_change():
    expr = pd.DataFrame([[0.3010] * 5 + [0.0] * 5], index=["G"])   # 0.301 log10 = 2-fold = 1 log2
    groups = pd.Series(["case"] * 5 + ["ctrl"] * 5, index=expr.columns)
    out = s3.differential_expression(expr, groups, "case", "ctrl", scale_to_log2=1 / s3.LOG10_2)
    assert out.loc["G", "log2FC"] == pytest.approx(1.0, abs=1e-3)


def test_planted_shared_genes_recovered(results):
    found = set(results["shared"].gene)
    recall = len(found & results["truth"]) / len(results["truth"])
    assert recall >= 0.85
    assert len(found - results["truth"]) <= 3        # few false positives


def test_overlap_is_significant_and_directions_reported(results):
    s = results["summary"].iloc[0]
    assert s.hypergeometric_p < 1e-10
    assert 40 <= s.shared_same_direction <= s.shared


def test_fdr_never_below_pvalue(results):
    ad = results["ad_a"].dropna(subset=["fdr"])
    assert (ad.fdr >= ad.p_value - 1e-12).all()


def test_exports_exist(results):
    s6.run(results["shared"])
    assert config.results_path("cytoscape_node_table.tsv").exists()
    n = len(config.results_path("string_gene_list.txt").read_text().split())
    assert n == len(results["shared"])


# ---------- subtype adjustment ----------
@pytest.fixture(scope="module")
def adjusted(results, tmp_path_factory):
    old = (config.CANCER_MODEL, config.RESULTS_DIR)
    config.CANCER_MODEL = "subtype"
    config.RESULTS_DIR = tmp_path_factory.mktemp("adj")
    try:
        de_ad, de_ca = s3.run(results["ad"], results["ca"])
        _, ca_a = s4.run(de_ad, de_ca)
    finally:
        config.CANCER_MODEL, config.RESULTS_DIR = old
    return ca_a


def test_subtypes_load_from_cbioportal_style_file(results):
    sub = sa.load_subtypes(list(results["ca"][0].columns))
    assert sub.notna().all() and set(sub) == {"BRCA_Basal", "BRCA_LumA", "BRCA_LumB", "BRCA_Her2"}


def test_subtype_adjustment_removes_subtype_driven_genes(results, adjusted):
    only = synth.SUBTYPE_ONLY
    assert results["ca_a"].loc[only, "significant"].sum() >= 35      # confounded: look like TP53 hits
    assert adjusted.loc[only, "significant"].sum() <= 3              # adjustment removes them


def test_subtype_adjustment_keeps_real_tp53_genes(adjusted):
    assert adjusted.loc[list(synth.SHARED), "significant"].sum() >= 55


# ---------- step 7: sensitivity ----------
def test_sensitivity_has_one_primary_row_matching_step5(results):
    table = s7.sensitivity_table(results["ad_a"], results["ca_a"], config.SENSITIVITY_GRID +
                                 [(config.AD_MIN_ABS_LOG2FC, config.CANCER_MIN_ABS_LOG2FC, config.FDR_THRESHOLD)])
    prim = table[table.analysis.str.startswith("PRIMARY")]
    assert len(prim) >= 1
    assert prim.iloc[0].shared == results["summary"].iloc[0].shared


# ---------- step 8: hand-off ----------
def test_tiers_are_consistent(results):
    tiered = s8.assign_tiers(results["shared"])
    assert set(tiered.tier) <= {"A_strict_concordant", "B_concordant", "C_discordant"}
    assert (tiered.loc[tiered.tier == "C_discordant", "same_direction"] == False).all()
    a = tiered[tiered.tier == "A_strict_concordant"]
    assert (a.fdr_AD < config.PRIORITY_FDR).all() and (a.log2FC_AD.abs() >= config.PRIORITY_AD_MIN_ABS_LOG2FC).all()
    assert len(a) >= 40


def test_watchlist_handles_old_and_new_symbols_and_missing(results):
    old = config.WATCHLIST
    config.WATCHLIST = {"alias test": ["NOT_REAL", "GENE0001"], "missing": ["NOPE"]}
    try:
        tiered = s8.assign_tiers(results["shared"])
        w = s8.watchlist_status(results["ad_a"], results["ca_a"], tiered)
    finally:
        config.WATCHLIST = old
    assert w.loc[0, "symbol_found"] == "GENE0001"
    assert w.loc[1, "symbol_found"] == "not in either dataset" and w.loc[1, "shared_tier"] == "not shared"


# ---------- step 9: enrichment + p53 axis ----------
def test_over_representation_finds_enriched_set_and_rejects_random(results):
    bg = list(results["ad_a"].index.intersection(results["ca_a"].index))
    sets = {"enriched": list(synth.SHARED[:30]), "random": bg[1000:1040]}
    res = s9.over_representation(list(results["shared"].gene), sets, bg)
    assert res.iloc[0].term == "enriched" and res.iloc[0].fdr < 1e-6
    assert "random" not in set(res.term)             # overlap below the minimum -> not reported


def test_p53_axis_detects_planted_opposite_direction(results):
    out = s9.p53_direction(results["ad_a"], results["ca_a"], synth.P53_SET).iloc[0]
    assert out.fraction_opposite > 0.9 and out.opposite_vs_genomewide_p < 1e-6


def test_gmt_reader(tmp_path):
    f = tmp_path / "x.gmt"
    f.write_text("setA\tdesc\tG1\tG2\tG3\nsetB\tdesc\tG4\n")
    assert s9.read_gmt(f) == {"setA": ["G1", "G2", "G3"], "setB": ["G4"]}
