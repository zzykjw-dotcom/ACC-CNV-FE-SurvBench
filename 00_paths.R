############################################################
# 00_paths.R
# Shared paths, seed and helpers (no machine-specific absolute paths)
############################################################

# Data directory: default = project root
# Override: Sys.setenv(ACC_DATA_DIR = "/path/to/raw/data")
.data_dir_env <- Sys.getenv("ACC_DATA_DIR", unset = "")
if (nzchar(.data_dir_env)) {
  data_dir <- .data_dir_env
} else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  data_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  data_dir <- getwd()
}

output_dir <- file.path(data_dir, "outputs")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (basename(getwd()) == "R") setwd("..")

set.seed(42)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
}

load_expr_clin <- function() {
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  library(data.table)

  expr_path <- file.path(data_dir, "ACC_expression_data.csv")
  clin_path <- file.path(data_dir, "ACC_clinical.csv")
  if (!file.exists(expr_path)) expr_path <- "ACC_expression_data.csv"
  if (!file.exists(clin_path)) clin_path <- "ACC_clinical.csv"

  expr <- fread(expr_path, data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  clin <- fread(clin_path)
  common <- intersect(colnames(expr), clin$sample_id)
  expr  <- expr[, common, drop = FALSE]
  clin  <- clin[match(common, clin$sample_id), ]

  list(expr = expr, clin = clin, common = common)
}

msg("Project / data directory:", data_dir)
msg("Output directory:", output_dir)
msg("Working directory:", getwd())
