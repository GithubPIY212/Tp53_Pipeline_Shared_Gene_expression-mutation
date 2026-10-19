
import gzip
import pandas as pd
from pathlib import Path

data_dir = Path("data")
expression_file = data_dir / "GSE33000_raw_data.txt.gz"
metadata_file = Path("GSE33000_sample_metadata.tsv")

metadata = pd.read_csv(metadata_file, sep="\t")
metadata["title"] = metadata["title"].astype(str).str.strip()

# Keep only AD and non-demented controls
metadata["disease_status"] = (
    metadata["characteristics_ch2"]
    .astype(str)
    .str.replace("disease status: ", "", regex=False)
    .str.strip()
)

metadata = metadata[
    metadata["disease_status"].isin(
        ["Alzheimer's disease", "non-demented"]
    )
].copy()

metadata["group"] = metadata["disease_status"].map({
    "Alzheimer's disease": "AD",
    "non-demented": "Control"
})

keep_titles = set(metadata["title"])

# Read matrix and retain matching sample columns
with gzip.open(expression_file, "rt", encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_cols = [
        col for col in header[3:] if col in keep_titles
    ]
    selected_cols = header[:3] + sample_cols

    output_path = Path("GSE33000_AD_control_expression.tsv")

    with output_path.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(selected_cols) + "\n")

        for line in f:
            fields = line.rstrip("\n").split("\t")
            out.write("\t".join(
                fields[:3] + [
                    fields[header.index(col)] for col in sample_cols
                ]
            ) + "\n")

# Save aligned metadata in expression-column order
ordered_metadata = metadata.set_index("title").loc[sample_cols].reset_index()
ordered_metadata = ordered_metadata.rename(
    columns={"title": "expression_sample_id"}
)
ordered_metadata.to_csv(
    "GSE33000_AD_control_expression_metadata.tsv",
    sep="\t", index=False
)

print("Expression samples retained:", len(sample_cols))
print("AD:", (ordered_metadata["group"] == "AD").sum())
print("Controls:", (ordered_metadata["group"] == "Control").sum())
print("Expression file:", output_path)
print("Metadata: GSE33000_AD_control_expression_metadata.tsv")
