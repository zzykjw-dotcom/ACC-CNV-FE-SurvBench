# ACC-CNV-FE-SurvBench

Reproducible analysis pipeline for:

**CNV-constrained feature engineering and multi-algorithm benchmarking for high-dimensional survival data: a case study on adrenocortical carcinoma**

Target journal: *Journal of Bioinformatics and Computational Biology* (JBCB)

## Overview

1. **CNV-FE** (feature engineering)
   - Stage 1: signed WGCNA → ME4 module → CNV–mRNA filter (56 genes)
   - Stage 2: LASSO-Cox → five-gene signature (MMP12, ASPM, OTX1, NCAPH, HELLS)

2. **SurvBench** (algorithm comparison)
   - Classical models (CoxPH, LASSO-Cox, RSF) on 8-dimensional input
   - Cox-AAE on 59-dimensional input (asymmetric design)
   - Bootstrap optimism-corrected C-index, SurvSHAP(t), calibration, DCA
   - External validation on GSE19750

Biological characterisation (GO/KEGG, PPI, immune, hub expression) is in a separate module.

## Repository layout

```
ACC-CNV-FE-SurvBench/
├── README.md
├── data/  "Note: the ME4 module contains 454 genes in total; 410 of these have available gene-level CNV data and were used for CNV-expression correlation analysis."                        # processed lists only (safe to commit)
│   ├── ME4_module_genes.txt
│   ├── ME4_CNV_driven_core_genes.txt
│   ├── five_gene_signature.txt
│   ├── clinical_characteristics.csv
│   └── README_data.md
├── R/
│   ├── 00_paths.R                 # relative paths + seed
│   ├── 01_WGCNA.R                 # algorithm
│   ├── 02_CNV_selection.R         # algorithm
│   ├── 03_LASSO_multi_algorithm.R # algorithm
│   ├── 04_bootstrap_SHAP_calibration_DCA.R
│   ├── 05_external_validation_GSE19750.R
│   └── 06_biological_GO_KEGG_PPI_immune.R   # biological validation
├── Python/
│   └── Cox_AAE.py
└── results/
    ├── figures/
    └── tables/
```

**Raw TCGA / GEO matrices are not stored in this repository.**  
Download them locally and point `ACC_DATA_DIR` to that folder (see below).

## Requirements

**R** (>= 4.3): `WGCNA`, `glmnet`, `randomForestSRC`, `survival`, `survminer`, `timeROC`, `rms`, `pec`, `dcurves`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `GSVA`, `limma`, `data.table`, `dplyr`, `ggplot2`, `pheatmap`, `igraph`

**Python** (>= 3.8): `torch`, `numpy`, `pandas`, `scikit-learn`, `lifelines`

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("WGCNA", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "GSVA", "limma"))
install.packages(c("glmnet", "randomForestSRC", "survival", "survminer", "timeROC",
                   "rms", "pec", "dcurves", "data.table", "dplyr", "ggplot2", "pheatmap", "igraph"))
```

## How to run (no absolute paths)

1. Clone this repository.
2. Download raw data into a local folder (e.g. `~/ACC_raw/`).
3. From the **project root**:

```r
# Optional: point to raw data location
Sys.setenv(ACC_DATA_DIR = path.expand("~/ACC_raw"))

# Run in order (algorithm track)
source("R/00_paths.R")
source("R/01_WGCNA.R")
source("R/02_CNV_selection.R")
source("R/03_LASSO_multi_algorithm.R")
source("R/04_bootstrap_SHAP_calibration_DCA.R")
source("R/05_external_validation_GSE19750.R")

# Biological validation
source("R/06_biological_GO_KEGG_PPI_immune.R")
```

```bash
# Deep model (from project root)
export ACC_DATA_DIR=~/ACC_raw   # optional
python Python/Cox_AAE.py
```

All random seeds are fixed (`set.seed(42)` / `random_state=42`).

## Key results (manuscript)

| Metric | Value |
|--------|-------|
| Apparent C-index (5-gene Cox) | 0.906 |
| Optimism-corrected C-index (B=1000) | 0.882 (95% CI 0.867–0.968) |
| Cox-AAE (59-dim, 5-fold) | 0.831 ± 0.066 |
| GSE19750 log-rank p | 0.00068 |
| GSE19750 3-year AUC | 0.989 |

## Citation

[To be completed after acceptance]

## License

Scripts are provided for reproducibility of the associated manuscript.  
Source data remain subject to TCGA and GEO terms of use.
