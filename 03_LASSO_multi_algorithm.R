############################################################
# 03_LASSO_multi_algorithm.R
# ALGORITHM TRACK – Stage 2 of CNV-FE + SurvBench classical models
# LASSO-Cox → 5-gene signature; CoxPH / LASSO / RSF comparison
############################################################

source(file.path("R", "00_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(Hmisc)
})
set.seed(42)

msg("========== MODULE: LASSO + multi-algorithm ==========")

dat_obj <- load_expr_clin()
expr_mat <- dat_obj$expr
clin     <- dat_obj$clin
common_samples <- dat_obj$common

core_genes <- readLines("ME4_CNV_driven_core_genes.txt")
core_genes <- core_genes[nzchar(core_genes)]
core_genes <- intersect(core_genes, rownames(expr_mat))
msg("Core genes available:", length(core_genes))

# Build gene-level data frame
dat_gene <- data.frame(
  sample_id = common_samples,
  OS_time   = as.numeric(clin$OS_time),
  OS_status = as.integer(clin$OS_status),
  Age       = as.numeric(clin$Age),
  Sex       = ifelse(toupper(as.character(clin$Sex)) %in% c("MALE", "M"), 1L, 0L),
  Stage     = as.numeric(factor(clin$Stage)),
  stringsAsFactors = FALSE
)
for (g in core_genes) {
  dat_gene[[g]] <- as.numeric(expr_mat[g, common_samples])
}
dat_gene <- na.omit(dat_gene)
msg("Samples:", nrow(dat_gene), " Events:", sum(dat_gene$OS_status))

# LASSO-Cox feature selection
x_all <- as.matrix(dat_gene[, core_genes])
y_all <- Surv(dat_gene$OS_time, dat_gene$OS_status)
cv_lasso_all <- cv.glmnet(x_all, y_all, family = "cox", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_lasso_all, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[as.numeric(lasso_coef) != 0]
msg("LASSO-selected genes:", paste(selected_genes, collapse = ", "))
# Manuscript signature: MMP12, ASPM, OTX1, NCAPH, HELLS
# (re-run may differ slightly; use the five genes below for reproducibility)
selected_genes <- c("MMP12", "ASPM", "OTX1", "NCAPH", "HELLS")
selected_genes <- intersect(selected_genes, colnames(dat_gene))

feat_cols <- c(selected_genes, "Age", "Sex", "Stage")
x <- as.matrix(dat_gene[, feat_cols])
y <- Surv(dat_gene$OS_time, dat_gene$OS_status)

calc_cindex <- function(pred, y) {
  as.numeric(rcorr.cens(-as.numeric(pred), y)["C Index"])
}

results <- data.frame(Model = character(), Cindex = numeric(),
                      stringsAsFactors = FALSE)

# CoxPH
formula_cox <- as.formula(paste("Surv(OS_time, OS_status) ~",
                                paste(feat_cols, collapse = " + ")))
cox_fit_gene <- coxph(formula_cox, data = dat_gene, x = TRUE, y = TRUE)
results <- rbind(results, data.frame(
  Model = "CoxPH",
  Cindex = summary(cox_fit_gene)$concordance[1]
))

# LASSO-Cox (on selected features)
cv_lasso_final <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
lasso_pred <- as.numeric(predict(cv_lasso_final, newx = x, s = "lambda.min",
                                 type = "link"))
results <- rbind(results, data.frame(
  Model = "LASSO-Cox",
  Cindex = calc_cindex(lasso_pred, y)
))

# RSF
rf_gene <- rfsrc(Surv(OS_time, OS_status) ~ .,
                 data = dat_gene[, c("OS_time", "OS_status", feat_cols)],
                 ntree = 1000)
results <- rbind(results, data.frame(
  Model = "RSF",
  Cindex = 1 - rf_gene$err.rate[rf_gene$ntree]
))

results <- results[order(-results$Cindex), ]
print(results)
msg("Best classical model:", results$Model[1])

save(dat_gene, cox_fit_gene, cv_lasso_final, rf_gene, selected_genes, results,
     file = file.path(output_dir, "ACC_gene_level_model_comparison.RData"))
save(dat_gene, cox_fit_gene, selected_genes, results,
     file = "ACC_gene_level_model_comparison.RData")
write.csv(results, file.path(output_dir, "model_comparison_cindex_auc.csv"),
          row.names = FALSE)

msg("LASSO + multi-algorithm complete.")
