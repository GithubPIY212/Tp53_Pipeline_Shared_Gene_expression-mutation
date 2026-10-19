"""
Build data/GSE33000_metadata.tsv from GSE33000_series_matrix.txt.gz (downloaded from GEO).

Usage (from the pipeline folder):
    python make_ad_metadata.py GSE33000_series_matrix.txt.gz
It finds the 'characteristics' row that contains "Alzheimer's disease" and saves
three columns: sample_id (GSM number), sample_title (e.g. PFC_1), disease_state.
"""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


def main(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    ids, titles, char_rows = None, None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                ids = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                char_rows.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    if ids is None:
        sys.exit("No '!Sample_geo_accession' line found - is this a GEO series matrix file?")

    states = None
    for row in char_rows:
        if any("alzheimer" in v.lower() for v in row):
            # values look like "disease state: Alzheimer's disease" -> keep text after the colon
            states = [v.split(":", 1)[1].strip() if ":" in v else v.strip() for v in row]
            break
    if states is None:
        sys.exit("No characteristics row mentions Alzheimer's. Open the file and check the wording.")

    out = Path("data")
    out.mkdir(exist_ok=True)
    with open(out / "GSE33000_metadata.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "sample_title", "disease_state"])
        w.writerows(zip(ids, titles or ids, states))
    print("Saved data/GSE33000_metadata.tsv")
    for label, n in Counter(states).items():
        print(f"  {n:4d}  {label}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python make_ad_metadata.py GSE33000_series_matrix.txt.gz")
    main(sys.argv[1])
