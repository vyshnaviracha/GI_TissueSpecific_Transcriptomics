# GI Tissue-Specific Transcriptomic Remodeling in Gastrointestinal Cancers

This repository contains reproducible code for generating tissue-specific tumor vs normal mRNA expression heatmaps across gastrointestinal (GI) cancers using the TCGA-TARGET-GTEx cohort.

---

## Study Objective

To characterize transcriptomic remodeling in gastrointestinal malignancies by comparing tumor and matched normal tissues across:

- Colorectum
- Stomach
- Liver
- Pancreas
- Esophagus

Gene-level TPM values were used for a curated GI-relevant gene panel.

---

## Data Source

TCGA-TARGET-GTEx cohort  
Downloaded from UCSC Xena Browser.

Expression data type:
- RSEM normalized TPM

Sample annotation:
- `TcgaTargetGTEX_phenotype.txt`

---

## Repository Structure

```
analysis/
    heatmap_analysis.ipynb

scripts/
    gi_tissue_specific_analysis.py

results/
    Figure3_GI_TissueSpecific_TumorVsNormal.svg
```

---

## Methods Summary

1. Extract GI-relevant samples using phenotype annotation.
2. Match TCGA tumor and GTEx normal samples by tissue type.
3. Compute mean TPM per tissue-condition group.
4. Perform gene-wise Z-score normalization.
5. Apply hierarchical clustering (Euclidean distance, average linkage).
6. Visualize with consistent red–blue color scaling.

Color scale:
- Red = high expression (relative to gene mean)
- Blue = low expression

---

## Reproducibility

To reproduce:

```
pip install -r requirements.txt
python scripts/gi_tissue_specific_analysis.py
```

---

## Output

Primary figure:
- `Figure3_GI_TissueSpecific_TumorVsNormal.svg`
