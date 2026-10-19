"""
config.py - ONE place for every path, label and threshold.

WHY: if a reviewer asks "what cut-off did you use?", the answer is in this file.
Steps read these values at run time, so run_pipeline.py can point them at demo data.
"""
from pathlib import Path

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"          # put your real input files here
RESULTS_DIR = ROOT / "results"    # every output is written here


def data_path(name):
    return Path(DATA_DIR) / name


def results_path(name):
    Path(RESULTS_DIR).mkdir(parents=True, exist_ok=True)
    return Path(RESULTS_DIR) / name


# ---------------- Alzheimer's arm (GEO GSE33000) ----------------
# Either file name is accepted (plain or gzipped).
AD_EXPRESSION_FILES = ["GSE33000_expression.tsv.gz", "GSE33000_expression.tsv"]
AD_METADATA_FILE = "GSE33000_metadata.tsv"
# Labels must match GEO's text EXACTLY (lowercase 'd' in disease!).
AD_CASE_LABEL = "Alzheimer's disease"
AD_CONTROL_LABEL = "non-demented"

# ---------------- Cancer arm (TCGA-BRCA from cBioPortal) ----------------
CANCER_EXPRESSION_FILE = "data_mrna_seq_v2_rsem.txt"   # genes x samples (RSEM)
CANCER_MUTATION_FILE = "data_mutations.txt"            # one row per mutation
TP53_GENE = "TP53"
# Samples with >=1 TP53 mutation = "TP53-altered"; the rest = "TP53-unaltered".
CANCER_CASE_LABEL = "TP53-altered"
CANCER_CONTROL_LABEL = "TP53-unaltered"
# Drop genes that are basically switched off (mean log2(RSEM+1) below this).
CANCER_MIN_MEAN_LOG2 = 1.0

# ---------------- Statistics ----------------
FDR_THRESHOLD = 0.05
# Effect-size cut-offs are per arm because microarray ratios are compressed
# compared with RNA-seq. Both are in log2 fold-change units.
AD_MIN_ABS_LOG2FC = 0.10
CANCER_MIN_ABS_LOG2FC = 0.50
MIN_SAMPLES_PER_GROUP = 3

# ---------------- Cancer model: how to deal with breast-cancer SUBTYPE ----------------
# TP53 mutations are far more common in basal-like tumours, and basal-like tumours have a very
# different expression profile. Without correction, "TP53-altered vs unaltered" partly measures
# "basal vs not basal" (mostly cell-division genes).
#   "none"    : plain TP53-altered vs unaltered (the primary, pre-set analysis)
#   "subtype" : same comparison, adjusted for subtype (linear model: gene ~ TP53 + subtype)
#   "stratum" : compare TP53-altered vs unaltered inside ONE subtype only (CANCER_STRATUM)
CANCER_MODEL = "none"
CANCER_CLINICAL_FILES = ["data_clinical_patient.txt", "data_clinical_sample.txt"]
CANCER_SUBTYPE_COLUMNS = ["SUBTYPE", "PAM50_SUBTYPE", "BRCA_SUBTYPE", "Subtype"]
CANCER_STRATUM = "BRCA_LumA"
MIN_SUBTYPE_SAMPLES = 20      # subtypes with fewer samples are left out of the adjusted model

# ---------------- Sensitivity analysis (step 7) ----------------
# (AD |log2FC| cut-off, cancer |log2FC| cut-off, FDR). The row matching the settings above is
# the pre-set PRIMARY analysis; every other row is post-hoc and must be reported as such.
SENSITIVITY_GRID = [
    (0.10, 0.50, 0.05), (0.20, 0.50, 0.05), (0.30, 0.50, 0.05), (0.30, 1.00, 0.05),
    (0.50, 1.00, 0.05), (0.20, 0.50, 0.001), (0.30, 0.50, 0.001), (0.40, 0.75, 0.001),
]

# ---------------- Hand-off list for the next part of the project (step 8) ----------------
PRIORITY_AD_MIN_ABS_LOG2FC = 0.30
PRIORITY_CANCER_MIN_ABS_LOG2FC = 1.00
PRIORITY_FDR = 0.001
PRIORITY_TOP_N = 50
# Genes the CRISPR part wants to watch. Several symbols = old/new names for the same gene.
WATCHLIST = {
    "SETDB1": ["SETDB1"], "EHMT2 (G9A)": ["EHMT2"], "cGAS": ["CGAS", "MB21D1"],
    "STING": ["STING1", "TMEM173"], "TREM2": ["TREM2"],
}

# ---------------- Pathway enrichment + p53 axis (step 9) ----------------
# Each entry is an Enrichr library name (needs internet) OR a path to a local .gmt file.
ENRICHMENT_LIBRARIES = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023",
                        "KEGG_2021_Human", "Reactome_2022"]
P53_LIBRARY = "MSigDB_Hallmark_2020"      # library that holds the p53 pathway gene set
P53_SET_KEYWORD = "p53 pathway"           # matched case-insensitively against set names
P53_CORE_TARGETS = ["TP53", "MDM2", "CDKN1A", "BAX", "GADD45A", "BBC3", "PMAIP1"]
ENRICH_MIN_OVERLAP = 3
ENRICH_MIN_SET_SIZE = 5
