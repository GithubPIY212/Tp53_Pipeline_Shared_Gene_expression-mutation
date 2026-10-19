
import re
import pandas as pd
import config


def clean_symbol(raw):
    """Extract the first gene symbol from an annotation field."""
    if pd.isna(raw):
        return None

    value = str(raw).strip()
    first = re.split(r"\s*(?:///|//|;|,)\s*", value)[0].strip()

    if not first or first.upper() in {"NA", "N/A", "---", "NONE"}:
        return None

    return first


def run():
    # Locate the configured expression file.
    path = config.data_path(config.AD_EXPRESSION_FILES[0])
    raw = pd.read_csv(path, sep="\t", low_memory=False)

    print("Input file:", path)
    print("Input columns:", raw.columns[:10].tolist())

    # Support raw and preprocessed expression matrices.
    gene_col = next(
        (c for c in ("gene_symbol", "Gene", "gene", "symbol")
         if c in raw.columns),
        None,
    )

    if gene_col is None:
        raise ValueError(
            f"No gene-symbol column found. "
            f"Columns: {raw.columns[:10].tolist()}"
        )

    annotation_cols = {
        "reporterID", "Gene", "Transcript", "gene_symbol", "gene", "symbol"
    }
    sample_cols = [c for c in raw.columns if c not in annotation_cols]

    if not sample_cols:
        raise ValueError("No expression sample columns were found.")

    # Build a genes x samples expression matrix.
    raw["gene_symbol"] = raw[gene_col].map(clean_symbol)
    raw = raw[raw["gene_symbol"].notna()].copy()

    raw[sample_cols] = raw[sample_cols].apply(
        pd.to_numeric, errors="coerce"
    )

    expr = raw.groupby("gene_symbol")[sample_cols].mean()

    # Metadata file: this must be a tabular file with sample/group columns.
    metadata_path = config.data_path(config.AD_METADATA_FILE)
    metadata = pd.read_csv(metadata_path, sep=None, engine="python")

    sample_id_col = config.AD_METADATA_SAMPLE_ID_COL
    group_col = config.AD_METADATA_GROUP_COL

    if sample_id_col not in metadata.columns or group_col not in metadata.columns:
        raise ValueError(
            f"Metadata needs columns '{sample_id_col}' and '{group_col}'. "
            f"Found: {metadata.columns.tolist()}"
        )

    metadata[sample_id_col] = metadata[sample_id_col].astype(str).str.strip()
    expr.columns = expr.columns.astype(str).str.strip()

    # Keep only samples that exist in both files and have the target labels.
    metadata = metadata[
        metadata[group_col].isin(
            [config.AD_DISEASE_LABEL, config.AD_CONTROL_LABEL]
        )
    ].copy()

    metadata = metadata.drop_duplicates(subset=sample_id_col).set_index(
        sample_id_col
    )

    matched = [s for s in expr.columns if s in metadata.index]

    if not matched:
        raise ValueError(
            "No metadata sample IDs match the expression columns. "
            f"Expression examples: {expr.columns[:8].tolist()}; "
            f"metadata examples: {metadata.index[:8].tolist()}. "
            "Check how PFC sample IDs map to the GEO metadata."
        )

    expr = expr[matched]
    groups = metadata.loc[matched, group_col]

    disease_samples = groups[
        groups == config.AD_DISEASE_LABEL
    ].index.tolist()
    control_samples = groups[
        groups == config.AD_CONTROL_LABEL
    ].index.tolist()

    if not disease_samples or not control_samples:
        raise ValueError(
            f"Matched {len(disease_samples)} disease samples and "
            f"{len(control_samples)} controls. Both groups must be present."
        )

    print("Expression matrix:", expr.shape)
    print("Matched samples:", len(matched))
    print("Alzheimer's disease samples:", len(disease_samples))
    print("Control samples:", len(control_samples))

    # Return a structure that includes the groups for downstream analysis.
    return {
        "expression": expr,
        "disease_samples": disease_samples,
        "control_samples": control_samples,
    }


if __name__ == "__main__":
    run()
