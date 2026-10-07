# Processed data (safe to commit)

This folder holds **processed** intermediate results only.
Raw TCGA expression / CNV matrices and GEO series matrices are **not** stored here
(they are large and already public).

## Files in this folder

| File | Description |
|------|-------------|
| `ME4_module_genes.txt` | 454 genes in the ME4 WGCNA module |
| `ME4_CNV_driven_core_genes.txt` | 56 genes with \|ρ\| > 0.3 and p < 0.05 (CNV–mRNA) |
| `five_gene_signature.txt` | Final LASSO-Cox signature (MMP12, ASPM, OTX1, NCAPH, HELLS) |
| `clinical_characteristics.csv` | Summary table matching manuscript Table 1 |

## Raw data (download locally; do not commit)

| File | Source |
|------|--------|
| TCGA-ACC RNA-seq | [GDC](https://portal.gdc.cancer.gov/) |
| Gene-level CNV | GDC absolute copy-number |
| Clinical annotations | GDC / cBioPortal |
| GSE19750 series matrix | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19750) |

Place raw files in a local directory of your choice and set `DATA_DIR`
in `R/00_paths.R` (or pass the path as an environment variable).
