"""
Subtype handling for the cancer arm (used by step 3 when config.CANCER_MODEL != "none").

Adjusted model, per gene:   expression ~ intercept + TP53_altered + subtype dummies
The TP53 coefficient is then the TP53 effect AFTER removing differences explained by subtype,
and it is reported as log2FC (expression is already log2). All genes are fitted at once with
matrix algebra, so it is fast.
"""
import numpy as np
import pandas as pd
from scipy import stats
import config


def _patient(sample_id):
    return str(sample_id)[:12]          # TCGA-XX-XXXX


def load_subtypes(sample_ids):
    """Return a Series (sample -> subtype). Looks for a SUBTYPE column in the cBioPortal clinical files."""
    tried = []
    for fname in config.CANCER_CLINICAL_FILES:
        path = config.data_path(fname)
        if not path.exists():
            tried.append(f"{fname} (missing)")
            continue
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
        col = next((c for c in df.columns if c in config.CANCER_SUBTYPE_COLUMNS), None)
        id_col = next((c for c in ("PATIENT_ID", "SAMPLE_ID") if c in df.columns), None)
        if col is None or id_col is None:
            tried.append(f"{fname} (no subtype column; columns: {list(df.columns)[:12]})")
            continue
        df["key"] = df[id_col].map(_patient)
        lookup = df.dropna(subset=[col]).drop_duplicates("key").set_index("key")[col]
        out = pd.Series([lookup.get(_patient(s)) for s in sample_ids], index=list(sample_ids), name="subtype")
        if out.notna().any():
            return out
        tried.append(f"{fname} (no sample IDs matched)")
    raise FileNotFoundError("Could not find tumour subtypes. Put the cBioPortal clinical files in data/ "
                            "or edit CANCER_CLINICAL_FILES / CANCER_SUBTYPE_COLUMNS in config.py. Tried: "
                            + "; ".join(tried))


def adjusted_de(expr, groups, subtypes, case, control):
    """Linear model for every gene: expr ~ TP53 + subtype. Returns the same columns as step 3's DE table."""
    usable = subtypes.notna()
    counts = subtypes[usable].value_counts()
    keep_types = counts[counts >= config.MIN_SUBTYPE_SAMPLES].index
    usable &= subtypes.isin(keep_types)
    samples = list(expr.columns[usable.to_numpy()])
    grp = groups[samples]
    in_two = grp.isin([case, control])
    samples = [s for s in samples if in_two[s]]

    y = expr[samples].to_numpy(float).T                       # samples x genes
    tp53 = (groups[samples] == case).astype(float).to_numpy()
    dummies = pd.get_dummies(subtypes[samples], drop_first=True).astype(float).to_numpy()
    X = np.column_stack([np.ones(len(samples)), tp53, dummies])
    n, p = X.shape
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ X.T @ y                                   # p x genes
    resid = y - X @ beta
    dof = n - p
    sigma2 = (resid ** 2).sum(axis=0) / dof
    se = np.sqrt(sigma2 * xtx_inv[1, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta[1] / se
    pval = 2 * stats.t.sf(np.abs(t), dof)

    is_case = tp53 == 1
    out = pd.DataFrame({
        "gene": expr.index,
        "mean_case": y[is_case].mean(axis=0),
        "mean_control": y[~is_case].mean(axis=0),
        "log2FC": beta[1],
        "t_stat": t,
        "p_value": pval,
    }).set_index("gene")
    print(f"[subtype model] {n} samples, subtypes used: {dict(subtypes[samples].value_counts())}")
    return out
