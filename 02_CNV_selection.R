############################################################
# 02_CNV_selection.R
# ALGORITHM TRACK – Stage 1 continued: CNV–mRNA filter
# Output: 56 CNV-driven core genes
############################################################

source(file.path("R", "00_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
})
set.seed(42)

msg("========== MODULE: CNV selection ==========")

# Expression
dat <- load_expr_clin()
expr_mat <- dat$expr
common_samples <- dat$common

# Gene-level CNV (Ensembl → Symbol)
cnv_raw <- fread("TCGA-ACC.gene-level_absolute.tsv",
                 data.table = FALSE, check.names = FALSE)
rownames(cnv_raw) <- cnv_raw[[1]]
cnv_raw[[1]] <- NULL
cnv_mat <- as.matrix(cnv_raw)
storage.mode(cnv_mat) <- "numeric"
colnames(cnv_mat) <- substr(colnames(cnv_mat), 1, 15)

ensembl_ids <- sub("\\..*$", "", rownames(cnv_mat))
symbol_map <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                     column = "SYMBOL", keytype = "ENSEMBL",
                     multiVals = "first")
keep <- !is.na(symbol_map)
cnv_mat <- cnv_mat[keep, ]
rownames(cnv_mat) <- symbol_map[keep]
cnv_mat <- cnv_mat[!duplicated(rownames(cnv_mat)), ]
msg("CNV matrix after Symbol mapping:", dim(cnv_mat))

# Align samples
common_samples <- intersect(colnames(cnv_mat), colnames(expr_mat))
cnv_mat  <- cnv_mat[, common_samples]
expr_mat <- expr_mat[, common_samples]
msg("Intersected samples:", length(common_samples))

# ME4 genes
me4_genes <- readLines("module_genes.txt")
me4_genes <- me4_genes[nzchar(me4_genes)]
me4_in_cnv <- intersect(me4_genes, rownames(cnv_mat))
me4_in_cnv <- intersect(me4_in_cnv, rownames(expr_mat))
msg("ME4 genes present in both CNV and expression:", length(me4_in_cnv))

# Spearman CNV–mRNA
cor_res <- data.frame(Gene = me4_in_cnv, Rho = NA_real_, Pval = NA_real_)
for (i in seq_along(me4_in_cnv)) {
  g <- me4_in_cnv[i]
  x <- as.numeric(cnv_mat[g, ])
  y <- as.numeric(expr_mat[g, ])
  if (sum(!is.na(x)) < 10 || sum(!is.na(y)) < 10) next
  sx <- sd(x, na.rm = TRUE); sy <- sd(y, na.rm = TRUE)
  if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) next
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  cor_res$Rho[i]  <- unname(ct$estimate)
  cor_res$Pval[i] <- ct$p.value
}
msg("Genes with valid correlation:", sum(!is.na(cor_res$Rho)))

core_genes <- cor_res$Gene[!is.na(cor_res$Rho) &
                             abs(cor_res$Rho) > 0.3 &
                             cor_res$Pval < 0.05]
msg("CNV-driven core genes:", length(core_genes), "(expect ~56)")

write.csv(cor_res, file.path(output_dir, "ME4_CNV_expression_correlation.csv"),
          row.names = FALSE)
writeLines(core_genes, file.path(output_dir, "ME4_CNV_driven_core_genes.txt"))
# also copy to working dir for downstream scripts
write.csv(cor_res, "ME4_CNV_expression_correlation.csv", row.names = FALSE)
writeLines(core_genes, "ME4_CNV_driven_core_genes.txt")

msg("CNV selection complete.")
