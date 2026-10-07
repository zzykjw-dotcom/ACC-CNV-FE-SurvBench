############################################################
# 04_bootstrap_SHAP_calibration_DCA.R
# ALGORITHM TRACK – SurvBench evaluation utilities
# Bootstrap optimism-corrected C-index, SurvSHAP(t),
# calibration (1/3/5 yr), DCA, NRI/IDI
############################################################

source(file.path("R", "00_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(rms)
  library(pec)
  library(dcurves)
  library(ggplot2)
})
set.seed(42)

msg("========== MODULE: Bootstrap / SHAP / Calibration / DCA ==========")

# ----- rebuild 5-gene + clinical data -----
dat_obj <- load_expr_clin()
expr <- dat_obj$expr
clin <- dat_obj$clin
common_samples <- dat_obj$common

genes5 <- c("MMP12", "ASPM", "OTX1", "NCAPH", "HELLS")
stopifnot(all(genes5 %in% rownames(expr)))

dat <- data.frame(
  sample_id = common_samples,
  OS_time   = as.numeric(clin$OS_time),
  OS_status = as.integer(clin$OS_status),
  Age       = as.numeric(clin$Age),
  Sex       = ifelse(toupper(as.character(clin$Sex)) %in% c("MALE", "M"), 1L, 0L),
  Stage     = as.numeric(factor(clin$Stage)),
  MMP12 = as.numeric(expr["MMP12", ]),
  ASPM  = as.numeric(expr["ASPM", ]),
  OTX1  = as.numeric(expr["OTX1", ]),
  NCAPH = as.numeric(expr["NCAPH", ]),
  HELLS = as.numeric(expr["HELLS", ]),
  stringsAsFactors = FALSE
)
dat <- na.omit(dat)
dat_gene <- dat
msg("Samples:", nrow(dat), " Events:", sum(dat$OS_status))

cox_formula <- Surv(OS_time, OS_status) ~
  MMP12 + ASPM + OTX1 + NCAPH + HELLS + Age + Sex + Stage
cox_fit <- coxph(cox_formula, data = dat, x = TRUE, y = TRUE)
app_c <- as.numeric(concordance(cox_fit)$concordance)
msg("Apparent C-index:", round(app_c, 4))

# ----- Bootstrap B = 1000 (LP sign corrected) -----
B <- 1000
c_boot <- numeric(B)
c_orig <- numeric(B)
n <- nrow(dat)
pb <- txtProgressBar(min = 0, max = B, style = 3)
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  d_b <- dat[idx, ]
  if (sum(d_b$OS_status) < 3L) {
    c_boot[b] <- NA_real_; c_orig[b] <- NA_real_; next
  }
  fit_b <- tryCatch(coxph(cox_formula, data = d_b, x = TRUE, y = TRUE),
                    error = function(e) NULL)
  if (is.null(fit_b)) {
    c_boot[b] <- NA_real_; c_orig[b] <- NA_real_; next
  }
  c_boot[b] <- as.numeric(concordance(fit_b)$concordance)
  lp_orig <- predict(fit_b, newdata = dat, type = "lp")
  # Higher LP = higher hazard → lower survival; use I(-lp) for concordance
  c_orig[b] <- as.numeric(
    concordance(Surv(dat$OS_time, dat$OS_status) ~ I(-lp_orig))$concordance
  )
  setTxtProgressBar(pb, b)
}
close(pb)

ok <- !is.na(c_boot) & !is.na(c_orig)
optimism    <- mean(c_boot[ok] - c_orig[ok])
c_corrected <- app_c - optimism
ci_boot     <- quantile(c_boot[ok], c(0.025, 0.975), na.rm = TRUE)

# rms confirmatory
dd <- datadist(dat)
options(datadist = "dd")
cph_fit <- cph(cox_formula, data = dat, x = TRUE, y = TRUE, surv = TRUE)
val <- validate(cph_fit, method = "boot", B = 500, dxy = TRUE)
c_rms_corr <- 0.5 + val["Dxy", "index.corrected"] / 2

cat("\n===== Bootstrap summary =====\n")
cat("Apparent C-index:           ", round(app_c, 4), "\n")
cat("Optimism:                   ", round(optimism, 4), "\n")
cat("Optimism-corrected C-index: ", round(c_corrected, 4), "\n")
cat("Bootstrap 95% CI:           ", round(ci_boot[1], 4), "–",
    round(ci_boot[2], 4), "\n")
cat("rms Corrected C:            ", round(c_rms_corr, 4), "\n")

sink(file.path(output_dir, "bootstrap_Cindex_summary.txt"))
cat("===== Bootstrap C-index Summary =====\n")
cat("Samples:", nrow(dat), "  Events:", sum(dat$OS_status), "\n")
cat("Apparent C-index:           ", round(app_c, 4), "\n")
cat("Optimism:                   ", round(optimism, 4), "\n")
cat("Optimism-corrected C-index: ", round(c_corrected, 4), "\n")
cat("Bootstrap 95% CI:           ", round(ci_boot[1], 4), " – ",
    round(ci_boot[2], 4), "\n")
cat("rms Corrected C:            ", round(c_rms_corr, 4), "\n")
sink()

# ----- Calibration 1/3/5 year -----
cal_times <- c(12, 36, 60)
pdf(file.path(output_dir, "Calibration_curves_1_3_5yr.pdf"),
    width = 12, height = 4)
par(mfrow = c(1, 3))
for (t in cal_times) {
  cal <- calibrate(cph_fit, cmethod = "KM", method = "boot",
                   B = 150, u = t, m = 15)
  plot(cal, xlab = paste0("Predicted ", t / 12, "-year survival"),
       ylab = paste0("Observed ", t / 12, "-year survival"),
       main = paste0(t / 12, "-year calibration"), subtitles = FALSE)
  abline(0, 1, col = "red", lty = 2)
}
dev.off()
msg("Calibration curves saved.")

# ----- DCA (5-gene vs clinical-only, 3-year) -----
cox_clin <- coxph(Surv(OS_time, OS_status) ~ Age + Sex + Stage,
                  data = dat, x = TRUE, y = TRUE)
times_dca <- 36
surv_5gene <- pec::predictSurvProb(cox_fit, newdata = dat, times = times_dca)
risk_5gene <- as.numeric(1 - surv_5gene)
surv_clin  <- pec::predictSurvProb(cox_clin, newdata = dat, times = times_dca)
risk_clin  <- as.numeric(1 - surv_clin)
dca_df <- data.frame(
  OS_time = dat$OS_time, OS_status = dat$OS_status,
  risk_5gene = risk_5gene, risk_clin = risk_clin
)
dca_res <- dca(Surv(OS_time, OS_status) ~ risk_5gene + risk_clin,
               data = dca_df, time = times_dca,
               thresholds = seq(0, 0.5, by = 0.01))
pdf(file.path(output_dir, "DCA_3year.pdf"), width = 7, height = 5)
print(plot(dca_res) +
        labs(title = "Decision Curve Analysis (3-year)",
             subtitle = "5-gene model vs Clinical-only model") +
        theme_bw(base_size = 12))
dev.off()
msg("DCA saved.")

# ----- NRI / IDI (optional package) -----
if (requireNamespace("survIDINRI", quietly = TRUE)) {
  library(survIDINRI)
  lp_old <- predict(cox_clin, type = "lp")
  lp_new <- predict(cox_fit, type = "lp")
  idi_nri <- IDI.INF(
    indata = cbind(dat$OS_time, dat$OS_status, lp_old, lp_new),
    t0 = 36, npert = 200
  )
  print(idi_nri)
  sink(file.path(output_dir, "NRI_IDI_results.txt"))
  print(idi_nri)
  sink()
} else {
  msg("survIDINRI not installed – NRI/IDI skipped.")
}

# ----- SurvSHAP(t) (optional: survex) -----
if (requireNamespace("survex", quietly = TRUE) &&
    requireNamespace("pec", quietly = TRUE)) {
  library(survex)
  feat_cols <- c(genes5, "Age", "Sex", "Stage")
  surv_pred <- function(model, newdata, times) {
    pec::predictSurvProb(model, newdata, times)
  }
  explainer <- explain_survival(
    model = cox_fit,
    data = dat[, feat_cols],
    y = Surv(dat$OS_time, dat$OS_status),
    predict_survival_function = surv_pred,
    label = "CoxPH (5-gene + clinical)"
  )
  set.seed(42)
  global_shap <- model_survshap(explainer = explainer, N = 50)
  pdf(file.path(output_dir, "SHAP_beeswarm_CoxPH.pdf"), width = 8, height = 6)
  print(plot(global_shap, geom = "beeswarm"))
  dev.off()
  msg("SurvSHAP beeswarm saved.")
} else {
  msg("survex not installed – SHAP skipped (use existing PDF if available).")
}

save(cox_fit, dat, dat_gene, app_c, optimism, c_corrected, ci_boot, c_rms_corr,
     file = file.path(output_dir, "bootstrap_Cindex_results.RData"))
msg("Bootstrap / calibration / DCA module complete.")
