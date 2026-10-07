############################################################
# 06_biological_GO_KEGG_PPI_immune.R
# BIOLOGICAL VALIDATION TRACK
# GO/KEGG of ME4, PPI hubs, GSE19750 hub expression, CIBERSORT
############################################################

source(file.path("R", "00_paths.R"))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(igraph)
})
set.seed(42)

msg("========== MODULE: Biological validation ==========")

module_genes <- readLines("module_genes.txt")
module_genes <- module_genes[nzchar(module_genes)]
msg("ME4 genes:", length(module_genes))

# ----- GO BP + KEGG -----
gene_df <- bitr(module_genes, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = gene_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                readable = TRUE)
ekegg <- enrichKEGG(gene = gene_df$ENTREZID, organism = "hsa",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ego),
          file.path(output_dir, "GO_BP_results.csv"), row.names = FALSE)
write.csv(as.data.frame(ekegg),
          file.path(output_dir, "KEGG_results.csv"), row.names = FALSE)

pdf(file.path(output_dir, "GO_bubble.pdf"), width = 10, height = 8)
print(dotplot(ego, showCategory = 15, title = "GO BP – ME4 module"))
dev.off()
pdf(file.path(output_dir, "KEGG_bar.pdf"), width = 9, height = 6)
print(barplot(ekegg, showCategory = 12, title = "KEGG – ME4 module"))
dev.off()
msg("Top GO terms:")
print(head(as.data.frame(ego)[, c("Description", "Count", "p.adjust")], 5))

# ----- PPI (local STRING TSV) -----
if (file.exists("string_interactions.tsv")) {
  ppi <- read.delim("string_interactions.tsv", stringsAsFactors = FALSE,
                    comment.char = "", check.names = FALSE)
  colnames(ppi)[1:2] <- c("from", "to")
  if (!"combined_score" %in% colnames(ppi)) {
    score_col <- grep("score", colnames(ppi), ignore.case = TRUE)[1]
    ppi$combined_score <- ppi[[score_col]]
  }
  ppi$combined_score <- as.numeric(ppi$combined_score)
  if (max(ppi$combined_score, na.rm = TRUE) <= 1) {
    ppi <- ppi[ppi$combined_score >= 0.4, ]
  } else {
    ppi <- ppi[ppi$combined_score >= 400, ]
  }
  g <- graph_from_data_frame(ppi[, c("from", "to")], directed = FALSE)
  E(g)$weight <- ppi$combined_score
  msg("PPI nodes:", vcount(g), "edges:", ecount(g))

  deg <- degree(g)
  hub_df <- data.frame(
    gene = names(deg),
    Degree = deg,
    Betweenness = betweenness(g, normalized = TRUE),
    Closeness = closeness(g, normalized = TRUE),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Degree))
  write.csv(hub_df, file.path(output_dir, "hub_genes_Degree.csv"),
            row.names = FALSE)
  msg("CDC20 degree rank:", which(hub_df$gene == "CDC20"))
  msg("Top 10:", paste(head(hub_df$gene, 10), collapse = ", "))

  top50 <- head(hub_df$gene, 50)
  top10 <- head(hub_df$gene, 10)
  sub_g <- induced_subgraph(g, intersect(top50, V(g)$name))
  V(sub_g)$color <- "grey70"
  V(sub_g)$size  <- 3
  V(sub_g)$color[V(sub_g)$name %in% top10] <- "orange"
  V(sub_g)$size[V(sub_g)$name %in% top10]  <- 8
  if ("CDC20" %in% V(sub_g)$name) {
    V(sub_g)$color[V(sub_g)$name == "CDC20"] <- "red"
    V(sub_g)$size[V(sub_g)$name == "CDC20"]  <- 12
  }
  pdf(file.path(output_dir, "PPI_top50.pdf"), width = 10, height = 10)
  set.seed(42)
  plot(sub_g, vertex.label = V(sub_g)$name, vertex.label.cex = 0.7,
       vertex.label.color = "black", vertex.label.dist = 0.5,
       layout = layout_with_fr(sub_g),
       main = "PPI network of top 50 hub genes (CDC20 highlighted)")
  dev.off()
} else {
  msg("string_interactions.tsv not found – PPI skipped.")
}

# ----- GSE19750 hub expression (ACC vs Normal) -----
gse_file <- "GSE19750_series_matrix.txt"
if (file.exists(gse_file)) {
  expr_gse <- read.delim(gse_file, comment.char = "!", header = TRUE,
                         row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)
  all_lines <- readLines(gse_file)
  title_line <- all_lines[grep("^!Sample_title", all_lines)][1]
  sample_titles <- gsub('^"|"$', "", strsplit(title_line, "\t")[[1]][-1])
  group <- ifelse(grepl("Normal", sample_titles, ignore.case = TRUE),
                  "Normal", "ACC")
  names(group) <- colnames(expr_gse)

  probe_map <- list(
    CDC20 = "202870_s_at", TOP2A = "201291_s_at", CCNA2 = "203418_at",
    CDK1  = "203213_at",   BUB1B = "203755_at",   AURKB = "209464_at",
    CDC45 = "204126_s_at", KIF23 = "204709_s_at", BUB1  = "209642_at",
    CDCA8 = "221591_s_at", FOXM1 = "202580_x_at"
  )
  genes <- names(probe_map)
  gene_expr <- sapply(probe_map, function(p) {
    if (p %in% rownames(expr_gse)) as.numeric(expr_gse[p, ])
    else rep(NA_real_, ncol(expr_gse))
  })
  rownames(gene_expr) <- colnames(expr_gse)
  gene_expr <- as.data.frame(gene_expr)
  gene_expr$Group <- group[rownames(gene_expr)]

  res_list <- lapply(genes, function(g) {
    acc  <- gene_expr[[g]][gene_expr$Group == "ACC"]
    norm <- gene_expr[[g]][gene_expr$Group == "Normal"]
    wt   <- wilcox.test(acc, norm)
    data.frame(
      Gene = g, n_ACC = sum(!is.na(acc)), n_Normal = sum(!is.na(norm)),
      median_ACC = median(acc, na.rm = TRUE),
      median_Normal = median(norm, na.rm = TRUE),
      log2FC = median(acc, na.rm = TRUE) - median(norm, na.rm = TRUE),
      wilcox_p = wt$p.value, stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, res_list)
  res$adj_p <- p.adjust(res$wilcox_p, method = "BH")
  res <- res[order(-res$log2FC), ]
  write.csv(res, file.path(output_dir, "GSE19750_hub_genes_expression.csv"),
            row.names = FALSE)
  msg("Hub genes upregulated (p<0.05 & log2FC>0):",
      sum(res$wilcox_p < 0.05 & res$log2FC > 0), "of", nrow(res))

  mat <- t(as.matrix(gene_expr[, genes]))
  mat_z <- t(scale(t(mat)))
  annotation_col <- data.frame(Group = group, row.names = colnames(mat_z))
  ann_colors <- list(Group = c(ACC = "#E64B35", Normal = "#4DBBD5"))
  pdf(file.path(output_dir, "GSE19750_hub_genes_heatmap.pdf"),
      width = 10, height = 5)
  pheatmap(mat_z, annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_cols = TRUE, cluster_rows = TRUE, show_colnames = FALSE,
           fontsize_row = 11,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = "Hub genes expression (GSE19750, row-scaled)")
  dev.off()
} else {
  msg("GSE19750_series_matrix.txt not found – hub expression skipped.")
}

# ----- Immune (CIBERSORT) -----
if (file.exists("CIBERSORT_fractions.csv") &&
    file.exists("module_eigenes.csv")) {
  immune <- read.csv("CIBERSORT_fractions.csv", check.names = FALSE)
  id_col <- colnames(immune)[1]
  cell_cols <- setdiff(colnames(immune),
                       c(id_col, "P-value", "Correlation", "RMSE",
                         "P.value", "Correlation.Coefficient"))
  rownames(immune) <- immune[[id_col]]
  MEs_i <- read.csv("module_eigenes.csv", row.names = 1, check.names = FALSE)
  dat_obj <- load_expr_clin()
  common <- intersect(rownames(immune), rownames(MEs_i))
  common <- intersect(common, colnames(dat_obj$expr))
  immune <- immune[common, cell_cols, drop = FALSE]
  immune[] <- lapply(immune, as.numeric)
  me4   <- MEs_i[common, "ME4"]
  cdc20 <- as.numeric(dat_obj$expr["CDC20", common])

  res_imm <- lapply(colnames(immune), function(cell) {
    x <- immune[[cell]]
    ct_me4 <- cor.test(me4, x, method = "spearman", exact = FALSE)
    ct_cdc <- cor.test(cdc20, x, method = "spearman", exact = FALSE)
    data.frame(
      Immune_Cell = cell,
      Rho_ME4 = unname(ct_me4$estimate), P_ME4 = ct_me4$p.value,
      Rho_CDC20 = unname(ct_cdc$estimate), P_CDC20 = ct_cdc$p.value,
      stringsAsFactors = FALSE
    )
  })
  res_imm <- do.call(rbind, res_imm)
  res_imm$P_ME4_adj   <- p.adjust(res_imm$P_ME4, method = "BH")
  res_imm$P_CDC20_adj <- p.adjust(res_imm$P_CDC20, method = "BH")
  res_imm <- res_imm[order(res_imm$Rho_ME4), ]
  write.csv(res_imm, file.path(output_dir, "Immune_Correlation_Results.csv"),
            row.names = FALSE)

  plot_df <- res_imm
  plot_df$label <- gsub("_", " ", gsub("_CIBERSORT", "", plot_df$Immune_Cell))
  plot_df$label <- factor(plot_df$label, levels = plot_df$label)
  long_df <- rbind(
    data.frame(Cell = plot_df$label, Rho = plot_df$Rho_ME4,
               Source = "ME4 module eigengene"),
    data.frame(Cell = plot_df$label, Rho = plot_df$Rho_CDC20,
               Source = "CDC20 expression")
  )
  p_lolli <- ggplot(long_df, aes(x = Rho, y = Cell, color = Source)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_segment(aes(x = 0, xend = Rho, y = Cell, yend = Cell),
                 linewidth = 0.6) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("ME4 module eigengene" = "#E64B35",
                                  "CDC20 expression" = "#4DBBD5")) +
    labs(x = "Spearman Rho", y = NULL,
         title = "Immune infiltration correlations") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(output_dir, "Figure10_immune_lollipop.pdf"),
         p_lolli, width = 8, height = 7)
  ggsave(file.path(output_dir, "Figure10_immune_lollipop.png"),
         p_lolli, width = 8, height = 7, dpi = 300)
  msg("Immune correlations written.")
} else {
  msg("CIBERSORT or module_eigenes missing – immune step skipped.")
}

msg("Biological validation module complete.")
