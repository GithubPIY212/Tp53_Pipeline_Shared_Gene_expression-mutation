
import pandas as pd
from pathlib import Path

file = Path(
    r"D:\Major Project\2nd Step\tp53_line1_pipeline\tp53_line1_pipeline\data\GSE33000_raw_data_Extract.txt"
)

output = Path(
    r"D:\Major Project\2nd Step\tp53_line1_pipeline\tp53_line1_pipeline\TNFSF12_TNFSF13_probe_rows.tsv"
)

targets = r"TNFSF12|TNFSF13"
matches = []

print("Searching for TNFSF12/TNFSF13...")

for chunk in pd.read_csv(
    file, sep="\t", chunksize=5000, dtype=str, low_memory=False
):
    # Search each row across all columns.
    mask = chunk.astype(str).apply(
        lambda col: col.str.contains(
            targets, case=False, na=False, regex=True
        )
    ).any(axis=1)

    if mask.any():
        matches.append(chunk.loc[mask])

if matches:
    result = pd.concat(matches, ignore_index=True)
    result.to_csv(output, sep="\t", index=False)

    print(f"Found {len(result)} matching rows.")
    print(f"Saved file to: {output}")
else:
    print("No matching gene symbols found.")
    print("The file may contain probe IDs rather than gene symbols.")

