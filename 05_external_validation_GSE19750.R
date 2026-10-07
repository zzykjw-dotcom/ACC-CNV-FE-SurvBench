############################################################
# 05_external_validation_GSE19750.R
# ALGORITHM TRACK – External survival validation
# GSE19750: 5-gene risk score → KM + time-dependent ROC
############################################################

source(file.path("R", "00_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(timeROC)
  library(org.Hs.eg.db)
})
set.seed(42)

msg("========== MODULE: GSE19750 external validation ==========")

gse_file <- "GSE19750_series_matrix.txt"
if (!file.exists(gse_file)) {
  stop("GSE19750_series_matrix.txt not found in working directory.")
}

# Expression matrix (skip GEO metadata lines)
expr_gse <- read.delim(gse_file, comment.char = "!", header = TRUE,
                       row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE)
gse_mat <- as.matrix(expr_gse)
storage.mode(gse_mat) <- "numeric"

# Probe → Symbol (Affymetrix HG-U133 Plus 2.0)
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("hgu133plus2.db", update = FALSE, ask = FALSE)
}
library(hgu133plus2.db)
probe_ids <- rownames(gse_mat)
symbol_map <- mapIds(hgu133plus2.db, keys = probe_ids,
                     column = "SYMBOL", keytype = "PROBEID",
                     multiVals = "first")
keep <- !is.na(symbol_map)
gse_mat2 <- gse_mat[keep, ]
rownames(gse_mat2) <- symbol_map[keep]
gse_mat2 <- gse_mat2[order(-apply(gse_mat2, 1, var, na.rm = TRUE)), ]
gse_mat2 <- gse_mat2[!duplicated(rownames(gse_mat2)), ]
msg("GSE19750 matrix after Symbol mapping:", dim(gse_mat2))

core_genes <- c("MMP12", "ASPM", "OTX1", "NCAPH", "HELLS")
msg("Gene match:", paste(core_genes %in% rownames(gse_mat2), collapse = ", "))

# Survival annotations from series matrix
all_lines <- readLines(gse_file)
surv_time_line   <- grep("survival in years", all_lines, value = TRUE)[1]
surv_status_line <- grep("survival status", all_lines, value = TRUE)[1]

extract_field <- function(line, pattern) {
  parts <- strsplit(line, "\t")[[1]][-1]
  parts <- gsub('"', "", parts)
  gsub(pattern, "", parts)
}
surv_time   <- extract_field(surv_time_line, "survival in years: ")
surv_status <- extract_field(surv_status_line, "survival status: ")

surv_df <- data.frame(
  sample    = colnames(gse_mat2),
  OS_time   = as.numeric(surv_time),
  OS_status = ifelse(surv_status == "dead", 1L,
                     ifelse(surv_status == "alive", 0L, NA_integer_)),
  stringsAsFactors = FALSE
)

# Linear predictor using TCGA-trained coefficients if available;
# otherwise equal-weight sum of z-scored genes as fallback.
gse_sub <- gse_mat2[intersect(core_genes, rownames(gse_mat2)), , drop = FALSE]
if (file.exists("ACC_gene_level_model_comparison.RData")) {
  load("ACC_gene_level_model_comparison.RData")
  coefs <- coef(cox_fit_gene)
  genes_use <- intersect(names(coefs), rownames(gse_sub))
  risk_score <- as.numeric(t(gse_sub[genes_use, , drop = FALSE]) %*%
                             coefs[genes_use])
} else if (exists("cox_fit")) {
  coefs <- coef(cox_fit)
  genes_use <- intersect(names(coefs), rownames(gse_sub))
  risk_score <- as.numeric(t(gse_sub[genes_use, , drop = FALSE]) %*%
                             coefs[genes_use])
} else {
  msg("No TCGA Cox object found – using sum of scaled expression.")
  risk_score <- colSums(scale(t(gse_sub)))
}
names(risk_score) <- colnames(gse_sub)
surv_df$risk_score <- risk_score[match(surv_df$sample, names(risk_score))]
surv_df <- na.omit(surv_df)
msg("External validation samples:", nrow(surv_df),
    " Events:", sum(surv_df$OS_status))

# KM
surv_df$group <- ifelse(surv_df$risk_score > median(surv_df$risk_score),
                        "High", "Low")
fit <- survfit(Surv(OS_time, OS_status) ~ group, data = surv_df)
pdf(file.path(output_dir, "GSE19750_KM.pdf"), width = 6, height = 5)
print(ggsurvplot(fit, data = surv_df, pval = TRUE, risk.table = TRUE,
                 xlab = "Time (years)", ylab = "Survival probability",
                 title = "GSE19750 external validation"))
dev.off()
print(survdiff(Surv(OS_time, OS_status) ~ group, data = surv_df))

# Time-dependent ROC (times in years)
roc_res <- timeROC(T = surv_df$OS_time, delta = surv_df$OS_status,
                   marker = surv_df$risk_score, cause = 1,
                   times = c(1, 3, 5), iid = TRUE)
msg("GSE19750 time-dependent AUC:")
print(roc_res$AUC)

pdf(file.path(output_dir, "GSE19750_timeROC.pdf"), width = 6, height = 5)
plot(roc_res, time = 1, col = "blue", title = "GSE19750 Time-dependent ROC")
plot(roc_res, time = 3, col = "red", add = TRUE)
plot(roc_res, time = 5, col = "green", add = TRUE)
legend("bottomright",
       legend = sprintf("%d-year AUC = %.3f", c(1, 3, 5), roc_res$AUC),
       col = c("blue", "red", "green"), lwd = 2)
dev.off()

write.csv(surv_df, file.path(output_dir, "GSE19750_risk_score.csv"),
          row.names = FALSE)
msg("External validation complete.")
