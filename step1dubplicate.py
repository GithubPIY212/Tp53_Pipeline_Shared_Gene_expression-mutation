"""
Step 1 - Load and clean the Alzheimer's (GSE33000) data.

Output: results/ad_expression_clean.tsv   (genes x samples, log10 ratio values)
        results/ad_groups.tsv             (sample -> AD / control)
"""
import re
import pandas as pd
import config


def clean_symbol(raw):
    """Microarray symbol fields can hold several aliases ('A1BG///ABG'). Keep the first."""
    if pd.isna(raw):
        return None
    first = re.split(r"\s*(?:///|//|;|,|\|)\s*", str(raw).strip())[0]
    return first.strip() or None


def _find_expression_file():
    for name in config.AD_EXPRESSION_FILES:
        if config.data_path(name).exists():
            return config.data_path(name)
    raise FileNotFoundError(
        f"No AD expression file found in {config.DATA_DIR}. Tried: {config.AD_EXPRESSION_FILES}")


def _detect_columns(meta, sample_ids):
    """Find (sample-id column, group column) without assuming their names."""
    id_col = max(meta.columns, key=lambda c: meta[c].astype(str).isin(sample_ids).sum())
    if meta[id_col].astype(str).isin(sample_ids).sum() == 0:
        raise ValueError(f"No metadata column matches expression sample names. Columns: {list(meta.columns)}")
    group_col = None
    for c in meta.columns:
        if (meta[c].astype(str) == config.AD_CASE_LABEL).any():
            group_col = c
            break
    if group_col is None:
        raise ValueError(
            f"Label '{config.AD_CASE_LABEL}' not found in any metadata column. "
            f"Check spelling/case against GEO. Columns: {list(meta.columns)}")
    return id_col, group_col


def run():
    expr = pd.read_csv(_find_expression_file(), sep="\t", index_col=0)
    expr.index = [clean_symbol(g) for g in expr.index]
    expr = expr[expr.index.notna()]
    # Several probes can map to one gene: average them so each gene appears once.
    expr = expr.groupby(level=0).mean()

    meta = pd.read_csv(config.data_path(config.AD_METADATA_FILE), sep="\t")
    id_col, group_col = _detect_columns(meta, set(expr.columns.astype(str)))
    meta[id_col] = meta[id_col].astype(str)
    meta = meta.set_index(id_col)

    # EXACT match keeps only AD and non-demented; Huntington's etc. fall away automatically.
    keep = meta[group_col].isin([config.AD_CASE_LABEL, config.AD_CONTROL_LABEL])
    meta = meta[keep]
    expr = expr[[s for s in meta.index if s in expr.columns]]
    groups = meta.loc[expr.columns, group_col].rename("group")

    expr.to_csv(config.results_path("ad_expression_clean.tsv"), sep="\t")
    groups.to_csv(config.results_path("ad_groups.tsv"), sep="\t")
    print(f"[step1] AD: {expr.shape[0]} genes x {expr.shape[1]} samples | "
          f"{(groups == config.AD_CASE_LABEL).sum()} AD vs {(groups == config.AD_CONTROL_LABEL).sum()} controls")
    return expr, groups


if __name__ == "__main__":
    run()
