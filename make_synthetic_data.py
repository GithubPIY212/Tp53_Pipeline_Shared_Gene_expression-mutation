"""
Builds FAKE data in the same file formats as the real downloads, with known planted signals:
  * 60 genes truly changed in BOTH arms (50 same direction, 10 opposite)
  * 40 genes changed in AD only, 150 genes changed in cancer only
If the pipeline cannot find the 60 planted shared genes, something is broken.
Run on real data ONLY after it passes here.
"""
import numpy as np
import pandas as pd
from pathlib import Path

N_GENES = 3000
SHARED = [f"GENE{i:04d}" for i in range(0, 60)]
AD_ONLY = [f"GENE{i:04d}" for i in range(60, 100)]
CANCER_ONLY = [f"GENE{i:04d}" for i in range(100, 250)]
OPPOSITE = set(SHARED[50:])   # last 10 shared genes go opposite ways
SUBTYPE_ONLY = [f"GENE{i:04d}" for i in range(250, 290)]   # driven by basal subtype, NOT by TP53
P53_SET = [f"GENE{i:04d}" for i in range(300, 340)]        # up in AD, mildly down in TP53-altered cancer


def build(out_dir, seed=7):
    rng = np.random.default_rng(seed)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes = [f"GENE{i:04d}" for i in range(N_GENES)]

    # ---- AD arm: log10(sample/reference) values, 60 AD + 50 controls + 20 Huntington's
    ad_samples = [f"GSM{i}" for i in range(130)]
    labels = ["Alzheimer's disease"] * 60 + ["non-demented"] * 50 + ["Huntington's disease"] * 20
    X = rng.normal(0, 0.12, (N_GENES, 130))
    for g in SHARED + AD_ONLY + P53_SET:
        i = genes.index(g)
        X[i, :60] += 0.15 if g not in OPPOSITE else -0.15   # ~0.5 log2 units
    names = [f"{g}///ALIAS{n}" if n % 10 == 0 else g for n, g in enumerate(genes)]
    pd.DataFrame(X, index=names, columns=ad_samples).to_csv(out / "GSE33000_expression.tsv", sep="\t")
    pd.DataFrame({"sample_id": ad_samples, "disease_state": labels}).to_csv(
        out / "GSE33000_metadata.tsv", sep="\t", index=False)

    # ---- Cancer arm: RSEM-like counts, 90 TP53-altered + 130 unaltered
    samples = [f"TCGA-AA-{i:04d}-01" for i in range(220)]
    L = rng.normal(6, 1.0, (N_GENES, 1)) + rng.normal(0, 0.6, (N_GENES, 220))
    L[2900:] = rng.normal(-3, 0.2, (100, 220))           # near-silent genes (should be filtered)
    for g in SHARED + CANCER_ONLY:
        i = genes.index(g)
        L[i, :90] += 1.2
    for g in P53_SET:                                    # mild loss in TP53-altered tumours
        L[genes.index(g), :90] -= 0.35
    # Subtype confounding: TP53-altered tumours are mostly basal, and basal tumours have their own program.
    rng2 = np.random.default_rng(seed + 1)               # separate stream so other draws are unchanged
    basal = np.zeros(220, bool)
    basal[:72] = True                                    # 72 of 90 altered tumours are basal
    basal[90 + rng2.choice(130, 13, replace=False)] = True   # 13 of 130 unaltered are basal
    others = rng2.choice(["BRCA_LumA", "BRCA_LumB", "BRCA_Her2"], size=220, p=[0.5, 0.3, 0.2])
    subtype = np.where(basal, "BRCA_Basal", others)
    for g in SUBTYPE_ONLY:
        L[genes.index(g), basal] += 1.5
    rsem = np.clip(2 ** L - 1, 0, None)
    df = pd.DataFrame(rsem, index=genes, columns=samples)
    df.insert(0, "Entrez_Gene_Id", range(N_GENES))
    df.index.name = "Hugo_Symbol"
    df.reset_index().to_csv(out / "data_mrna_seq_v2_rsem.txt", sep="\t", index=False)

    mut = [("TP53", s + "") for s in samples[:90]] + [("PIK3CA", s) for s in samples[90:130]]
    pd.DataFrame(mut, columns=["Hugo_Symbol", "Tumor_Sample_Barcode"]).to_csv(
        out / "data_mutations.txt", sep="\t", index=False)
    clin = pd.DataFrame({"PATIENT_ID": [x[:12] for x in samples], "SUBTYPE": subtype})
    with open(out / "data_clinical_patient.txt", "w") as fh:     # cBioPortal-style: '#' comment lines first
        fh.write("#Patient Identifier\tSubtype\n#Patient identifier\tSubtype\n#STRING\tSTRING\n#1\t1\n")
        clin.to_csv(fh, sep="\t", index=False)
    sets = {"p53 Pathway": P53_SET, "Demo shared set": SHARED[:30],
            "Random set A": genes[1000:1040], "Random set B": genes[1500:1530]}
    with open(out / "demo_sets.gmt", "w") as fh:
        for name, members in sets.items():
            fh.write(name + "\tdemo\t" + "\t".join(members) + "\n")
    return set(SHARED)


if __name__ == "__main__":
    build("data_demo")
    print("Synthetic data written to data_demo/")
