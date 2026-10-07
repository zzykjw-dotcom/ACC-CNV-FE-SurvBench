############################################################
# 01_WGCNA.R
# ALGORITHM TRACK – Stage 1 of CNV-FE: signed WGCNA
# Identifies ME4 module (survival-linked mitotic module)
############################################################

source(file.path("R", "00_paths.R"))   # or paste 00_paths content if scripts run standalone

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(pheatmap)
})
cor <- WGCNA::cor
allowWGCNAThreads()
set.seed(42)

msg("========== MODULE: WGCNA ==========")

# Prefer existing intermediates for exact manuscript numbers
use_existing <- file.exists("module_genes.txt") &&
                file.exists("module_eigenes.csv")

if (use_existing) {
  msg("Loading existing WGCNA outputs...")
  module_genes <- readLines("module_genes.txt")
  module_genes <- module_genes[nzchar(module_genes)]
  MEs <- read.csv("module_eigenes.csv", row.names = 1, check.names = FALSE)
  softPower <- 7
  msg("ME4 genes:", length(module_genes), "(expect 454)")
  msg("CDC20 in module:", "CDC20" %in% module_genes)
} else {
  msg("Re-running WGCNA from top 3,000 variable genes...")
  dat <- load_expr_clin()
  expr_raw <- dat$expr
  gene_vars <- apply(expr_raw, 1, var, na.rm = TRUE)
  keep_3k <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
  datExpr <- t(expr_raw[keep_3k, ])

  powers <- 1:20
  sft <- pickSoftThreshold(datExpr, powerVector = powers,
                           networkType = "signed", verbose = 0)
  softPower <- sft$powerEstimate
  if (is.na(softPower)) softPower <- 7
  msg("Selected soft power:", softPower)

  net <- blockwiseModules(
    datExpr,
    power             = softPower,
    networkType       = "signed",
    TOMType           = "signed",
    minModuleSize     = 15,
    mergeCutHeight    = 0.25,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    verbose           = 2
  )
  MEs <- net$MEs
  moduleLabels <- net$colors

  # Align clinical for trait correlation
  clin <- dat$clin
  trait <- data.frame(
    OS.time   = as.numeric(clin$OS_time),
    OS.status = as.integer(clin$OS_status),
    row.names = colnames(expr_raw)
  )
  moduleTraitCor <- cor(MEs, trait, use = "pairwise.complete.obs")
  keyME <- names(which.min(moduleTraitCor[, "OS.time"]))
  keyLabel <- as.numeric(gsub("ME", "", keyME))
  module_genes <- colnames(datExpr)[moduleLabels == keyLabel]

  write.table(module_genes, file.path(output_dir, "module_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.csv(MEs, file.path(output_dir, "module_eigenes.csv"))
  write.csv(trait, file.path(output_dir, "trait_data.csv"))
  msg("Key module:", keyME, "genes =", length(module_genes))
}

cor <- stats::cor
msg("WGCNA module complete.")
