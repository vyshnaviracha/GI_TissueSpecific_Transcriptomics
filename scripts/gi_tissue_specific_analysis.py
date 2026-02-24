"""
GI Tissue-Specific Tumor vs Normal Transcriptomic Analysis
Author: Your Name
Journal Target: npj Systems Biology and Applications

Description:
Generates tissue-specific tumor vs normal heatmap (GI cancers)
Extracts top 2% high/low expressed genes per cancer
Exports publication-ready SVG figure and CSV tables
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ==============================
# CONFIGURATION
# ==============================

EXPR_FILE = "GI_1396_Annotated_TPM.tsv"
PHENO_FILE = "TcgaTargetGTEX_phenotype.txt"

OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TOP_PERCENT = 0.02

# ==============================
# LOAD DATA
# ==============================

print("Loading expression data...")
expr = pd.read_csv(EXPR_FILE, sep="\t")
expr = expr.set_index("symbol")
expr = expr.select_dtypes(include=[np.number])

print("Loading phenotype data...")
pheno = pd.read_csv(PHENO_FILE, sep="\t", encoding="latin1")

# ==============================
# MATCH SAMPLES
# ==============================

common_samples = list(set(expr.columns) & set(pheno["sample"]))
expr = expr[common_samples]
pheno = pheno[pheno["sample"].isin(common_samples)]

print("Expression matrix:", expr.shape)
print("Phenotype matrix:", pheno.shape)

# ==============================
# DEFINE GI TISSUES
# ==============================

gi_sites = {
    "Colorectum": ["Colon", "Rectum"],
    "Stomach": ["Stomach"],
    "Liver": ["Liver"],
    "Pancreas": ["Pancreas"],
    "Esophagus": ["Esophagus"]
}

result = {}

for tissue, keywords in gi_sites.items():

    tissue_samples = pheno[
        pheno["_primary_site"].str.contains("|".join(keywords), case=False, na=False)
    ]

    tumor_ids = tissue_samples[
        tissue_samples["_sample_type"] == "Primary Tumor"
    ]["sample"].tolist()

    normal_ids = tissue_samples[
        tissue_samples["_study"] == "GTEX"
    ]["sample"].tolist()

    tumor_ids = list(set(tumor_ids) & set(expr.columns))
    normal_ids = list(set(normal_ids) & set(expr.columns))

    print(f"{tissue}: Tumor={len(tumor_ids)}, Normal={len(normal_ids)}")

    if len(tumor_ids) > 20 and len(normal_ids) > 20:
        result[f"{tissue}_Normal"] = expr[normal_ids].mean(axis=1)
        result[f"{tissue}_Tumor"] = expr[tumor_ids].mean(axis=1)

heatmap_df = pd.DataFrame(result)

# ==============================
# Z-SCORE NORMALIZATION
# ==============================

heatmap_scaled = heatmap_df.sub(
    heatmap_df.mean(axis=1), axis=0
).div(
    heatmap_df.std(axis=1), axis=0
)

heatmap_scaled = heatmap_scaled.replace([np.inf, -np.inf], np.nan)
heatmap_scaled = heatmap_scaled.dropna()

print("Final matrix shape:", heatmap_scaled.shape)

# ==============================
# GENERATE HEATMAP
# ==============================

sns.set(style="white")

g = sns.clustermap(
    heatmap_scaled,
    cmap="RdBu_r",
    center=0,
    vmin=-2,
    vmax=2,
    figsize=(12, 24),
    row_cluster=True,
    col_cluster=False,
    yticklabels=False
)

plt.title("GI Tissue-Specific Tumor vs Normal Expression", pad=120)

heatmap_path = os.path.join(OUTPUT_DIR, "Figure3_GI_TissueSpecific_TumorVsNormal.svg")
g.savefig(heatmap_path, format="svg", dpi=600)

print("Heatmap saved:", heatmap_path)

# ==============================
# DIFFERENTIAL ANALYSIS
# ==============================

delta_df = {}

for col in heatmap_scaled.columns:
    if "_Tumor" in col:
        tissue = col.replace("_Tumor","")
        normal_col = f"{tissue}_Normal"
        delta_df[tissue] = (
            heatmap_scaled[col] - heatmap_scaled[normal_col]
        )

delta_df = pd.DataFrame(delta_df)

n_genes = int(len(delta_df) * TOP_PERCENT)

print("Extracting top", n_genes, "genes per cancer")

for tissue in delta_df.columns:
    
    ranked = delta_df[tissue].sort_values(ascending=False)
    
    high = ranked.head(n_genes)
    low = ranked.tail(n_genes)
    
    high.to_csv(os.path.join(OUTPUT_DIR, f"{tissue}_Top2percent_High.csv"))
    low.to_csv(os.path.join(OUTPUT_DIR, f"{tissue}_Top2percent_Low.csv"))

print("Top 2% gene tables saved.")
print("Analysis complete.")
