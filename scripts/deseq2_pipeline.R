# -----------------------------
# 0. Load packages
# -----------------------------
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(cowplot)
library(pheatmap)
library(stringr)

# -----------------------------
# 1. Load data
# -----------------------------
counts <- read.csv("readcount.csv")
rownames(counts) <- counts$geneID
counts$geneID <- NULL

coldata <- read.csv("metadata.csv")
rownames(coldata) <- coldata$Sample
coldata$Sample <- NULL

coldata <- coldata[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(coldata)))


coldata$Condition <- trimws(coldata$Condition)
coldata$Condition <- factor(coldata$Condition)
coldata$Condition <- relevel(coldata$Condition, ref = "Veh")

coldata$Batch <- factor(coldata$Batch)

# -----------------------------
# 2. Clean metadata
# -----------------------------
coldata$Condition <- trimws(coldata$Condition)
coldata$Condition <- factor(coldata$Condition)
coldata$Condition <- relevel(coldata$Condition, ref = "Veh")

coldata$Batch <- factor(coldata$Batch)

# -----------------------------
# 3. DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ Batch + Condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

# -----------------------------
# 4. Gene ID annotation
# Build master zebrafish annotation table once
# biomaRt first, org.Dr.eg.db fallback
# -----------------------------
gene_ids <- rownames(dds)

# biomaRt / Ensembl annotation
tryCatch({
  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "drerio_gene_ensembl",
    mirror = "useast"
  )
  
  bm_annot <- biomaRt::getBM(
    attributes = c("entrezgene_id", "external_gene_name", "description"),
    filters = "entrezgene_id",
    values = gene_ids,
    mart = mart
  )
}, error = function(e) {
  message("BioMart failed; using org.Dr.eg.db annotation only.")
  bm_annot <<- data.frame(
    entrezgene_id = character(),
    external_gene_name = character(),
    description = character()
  )
})

bm_annot <- bm_annot %>%
  dplyr::mutate(
    gene = as.character(entrezgene_id),
    bm_symbol = external_gene_name,
    bm_description = description
  ) %>%
  dplyr::select(
    gene,
    bm_symbol,
    bm_description
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# org.Dr.eg.db fallback annotation
org_annot <- AnnotationDbi::select(
  org.Dr.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL", "GENENAME", "ZFIN")
)

org_annot <- org_annot %>%
  dplyr::mutate(
    gene = as.character(ENTREZID),
    org_symbol = SYMBOL,
    org_description = GENENAME,
    org_zfin_id = ZFIN
  ) %>%
  dplyr::select(
    gene,
    org_symbol,
    org_description,
    org_zfin_id
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Combine annotation sources
gene_annot <- tibble(gene = as.character(gene_ids)) %>%
  dplyr::left_join(bm_annot, by = "gene") %>%
  dplyr::left_join(org_annot, by = "gene") %>%
  dplyr::mutate(
    symbol = dplyr::coalesce(bm_symbol, org_symbol, gene),
    gene_description = dplyr::coalesce(
      bm_description,
      org_description,
      "Not annotated"
    ),
    zfin_id = dplyr::coalesce(
      org_zfin_id,
      "Not annotated"
    ),
    
    placeholder_symbol = grepl("^LOC|^si:|^[0-9]+$", symbol),
    
    useful_description =
      !is.na(gene_description) &
      gene_description != "Not annotated" &
      !grepl(
        "hypothetical|uncharacterized|novel|predicted|LOC|si:",
        gene_description,
        ignore.case = TRUE
      ),
    
    display_label = dplyr::case_when(
      !placeholder_symbol ~ symbol,
      placeholder_symbol & useful_description ~ gene_description,
      TRUE ~ symbol
    ),
    
    display_label = stringr::str_remove(display_label, " \\[.*$"),
    display_label = stringr::str_trunc(display_label, width = 28)
  ) %>%
  dplyr::select(
    gene,
    zfin_id,
    symbol,
    display_label,
    gene_description)

gene_annot <- gene_annot %>%
  dplyr::mutate(
    symbol = dplyr::coalesce(symbol_override, symbol),
    display_label = dplyr::coalesce(symbol_override, display_label),
    gene_description = dplyr::coalesce(description_override, gene_description)
  ) %>%
  dplyr::select(
    gene,
    zfin_id,
    symbol,
    display_label,
    gene_description
  )

# -----------------------------
# 4b. PCA plot
# -----------------------------

vsd <- vst(dds)

pca_plot <- plotPCA(vsd, intgroup = "Condition") +
  theme_classic(base_size = 14)

ggsave(
  filename = "DEG_PCA_plot.pdf",
  plot = pca_plot,
  width = 6,
  height = 5
)

# -----------------------------
# 5. Comparisons
# -----------------------------
conds <- levels(dds$Condition)
comparisons <- conds[conds != "Veh"]

top_tables <- list()
full_tables <- list()

# -----------------------------
# 6. Loop through contrasts
# -----------------------------
for (cond in comparisons) {
  
  message("Processing: ", cond, " vs Veh")
  
  coef_name <- resultsNames(dds)[
    resultsNames(dds) == paste0("Condition_", cond, "_vs_Veh")
  ]
  
  if (length(coef_name) != 1) {
    stop(paste("Could not uniquely identify coefficient for", cond))
  }
  
  res <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene")
}

# -----------------------------
# 7. Add annotations
# -----------------------------
res_df <- res_df %>%
  dplyr::left_join(gene_annot, by = "gene")

# -----------------------------
# Classification
# -----------------------------
fc_cutoff <- 1
padj_cutoff <- 0.05
p_cutoff <- -log10(padj_cutoff)

res_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(ifelse(is.na(padj), pvalue, padj)),
    
    regulation = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange > fc_cutoff ~ "up",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -fc_cutoff ~ "down",
      is.na(padj) & pvalue < 0.001 & log2FoldChange > fc_cutoff ~ "up",
      is.na(padj) & pvalue < 0.001 & log2FoldChange < -fc_cutoff ~ "down",
      TRUE ~ "ns"
    )
  )

n_up <- sum(res_df$regulation == "up", na.rm = TRUE)
n_down <- sum(res_df$regulation == "down", na.rm = TRUE)
n_ns <- sum(res_df$regulation == "ns", na.rm = TRUE)
n_total <- n_up + n_down

# -----------------------------
# Complete DEG table for this contrast
# -----------------------------
full_table <- res_df %>%
  dplyr::mutate(
    contrast = paste0(cond, "_vs_Veh"),
    direction = regulation,
    entrez_id = gene
  ) %>%
  dplyr::select(
    contrast,
    direction,
    entrez_id,
    zfin_id,
    symbol,
    display_label,
    gene_description,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  dplyr::arrange(desc(abs(log2FoldChange)))

full_tables[[paste0(cond, "_vs_Veh")]] <- full_table

# -----------------------------
# 8. Top 5 up/down 
# -----------------------------
#top_up <- res_df %>%
#filter(regulation == "up") %>%
#arrange(desc(log2FoldChange)) %>%
#slice_head(n = 5)

#top_down <- res_df %>%
#filter(regulation == "down") %>%
#arrange(log2FoldChange) %>%
#slice_head(n = 5)

# top_genes <- bind_rows(top_up, top_down) %>%
#  mutate(
#   contrast = paste0(cond, "_vs_Veh"),
#  direction = regulation
#)
# -----------------------------
# 8. Label Top 25 DEGs
# -----------------------------
top_genes <- res_df %>%
  filter(regulation != "ns") %>%
  arrange(padj) %>%
  slice_head(n = 25)

top_fc_genes <- res_df %>%
  filter(!is.na(log2FoldChange), regulation != "ns") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 5)

top_table <- top_genes %>%
  dplyr::mutate(
    contrast = paste0(cond, "_vs_Veh"), 
    direction = regulation,
    entrez_id = gene
  ) %>%
  dplyr::select(
    contrast,
    direction,
    entrez_id,
    zfin_id,
    symbol,
    display_label,
    gene_description,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  )

top_tables[[paste0(cond, "_vs_Veh")]] <- top_table
# -----------------------------
# Heatmap of top 50 DEGs for this contrast
# -----------------------------
heatmap_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene)

heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(vsd)]

# Extract VST values
heatmap_mat <- assay(vsd)[heatmap_genes, ]

# Collapse replicates by treatment group
condition_groups <- colData(dds)$Condition

heatmap_mat_collapsed <- sapply(
  levels(condition_groups),
  function(condition) {
    rowMeans(heatmap_mat[, condition_groups == condition, drop = FALSE])
  }
)

# Use clean display labels
heatmap_symbols <- gene_annot$display_label[
  match(rownames(heatmap_mat_collapsed), gene_annot$gene)
]

heatmap_symbols[is.na(heatmap_symbols)] <- rownames(heatmap_mat_collapsed)

rownames(heatmap_mat_collapsed) <- make.unique(heatmap_symbols)

# Treatment-only annotation
annotation_col <- data.frame(
  Condition = colnames(heatmap_mat_collapsed)
)

rownames(annotation_col) <- colnames(heatmap_mat_collapsed)

pheatmap(
  heatmap_mat_collapsed,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 6,
  filename = paste0(cond, "_vs_Veh_top50_heatmap_collapsed.pdf"),
  width = 6,
  height = 10
)

# -----------------------------
# 9. Volcano plot
# -----------------------------
# Dynamic x-axis limit based on max absolute fold change
x_limit <- max(abs(res_df$log2FoldChange), na.rm = TRUE)
x_limit <- ceiling(x_limit + 0.5)

p <- ggplot(res_df, aes(log2FoldChange, neg_log10_padj)) +
  
  annotate("rect",
           xmin = fc_cutoff, xmax = Inf,
           ymin = p_cutoff, ymax = Inf,
           fill = "red", alpha = 0.08) +
  
  annotate("rect",
           xmin = -Inf, xmax = -fc_cutoff,
           ymin = p_cutoff, ymax = Inf,
           fill = "darkgreen", alpha = 0.08) +
  
  annotate("rect",
           xmin = -fc_cutoff, xmax = fc_cutoff,
           ymin = -Inf, ymax = p_cutoff,
           fill = "gray", alpha = 0.03) +
  
  geom_point(aes(color = regulation), size = 1, alpha = 0.9) +
  
  geom_hline(yintercept = p_cutoff, linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  
  geom_point(
    data = top_genes,
    shape = 21,
    size = 4,
    stroke = 1.2,
    color = "black",
    fill = NA
  ) +
  #  geom_point(
  #   data = top_fc_genes,
  #  shape = 21,
  # size = 5,
  #stroke = 1.2,
  #color = "purple",
  #fill = NA
  #) +
  geom_text_repel(
    data = top_genes,
    aes(label = display_label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.3
  ) +
  # geom_text_repel(
  #  data = top_fc_genes,
  # aes(label = display_label),
  #  size = 3.2,
  # color = "purple",
  #  max.overlaps = Inf,
  # box.padding = 0.5,
  #  point.padding = 0.3,
  # segment.size = 0.3
  #) +
  coord_cartesian(xlim = c(-x_limit, x_limit)) +
  
  scale_color_manual(values = c(
    "up" = "red",
    "down" = "darkgreen",
    "ns" = "gray"
  )) +
  
  labs(
    title = paste(cond, "vs Veh"),
    x = expression(log[2](fold~change)),
    y = expression(-log[10](padj))
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# -----------------------------
# 10. Separate annotation panel
# -----------------------------
p_annot <- ggplot() +
  annotate("text", x = 0, y = 0.90,
           label = paste0("DEGs (", n_total, ")"),
           hjust = 0, size = 5, fontface = "bold") +
  
  annotate("point", x = 0.05, y = 0.72,
           color = "red", size = 3.5) +
  annotate("text", x = 0.14, y = 0.72,
           label = paste0("up: ", n_up),
           hjust = 0, size = 4.2) +
  
  annotate("point", x = 0.05, y = 0.62,
           color = "darkgreen", size = 3.5) +
  annotate("text", x = 0.14, y = 0.62,
           label = paste0("down: ", n_down),
           hjust = 0, size = 4.2) +
  
  annotate("point", x = 0.05, y = 0.52,
           color = "gray", size = 3.5) +
  annotate("text", x = 0.14, y = 0.52,
           label = paste0("not significant: ", n_ns),
           hjust = 0, size = 4.2) +
  
  xlim(0, 1) +
  ylim(0.45, 0.95) +
  theme_void()

p_final <- plot_grid(
  p,
  p_annot,
  nrow = 1,
  rel_widths = c(3.6, 1.0)
)

ggsave(
  filename = paste0(cond, "_vs_Veh_volcano.pdf"),
  plot = p_final,
  width = 10,
  height = 6,
  dpi = 1200)

# -----------------------------
# 11. Export Excel workbooks
# -----------------------------
wb_top <- createWorkbook()

for (sheet in names(top_tables)) {
  addWorksheet(wb_top, sheetName = sheet)
  writeData(wb_top, sheet = sheet, top_tables[[sheet]])
}

combined_top <- bind_rows(top_tables)

addWorksheet(wb_top, sheetName = "All_Contrasts")
writeData(wb_top, sheet = "All_Contrasts", combined_top)

saveWorkbook(
  wb_top,
  file = "Top5DEGs_by_contrast.xlsx",
  overwrite = TRUE
)

# Complete DEG workbook
wb_full <- createWorkbook()

for (sheet in names(full_tables)) {
  addWorksheet(wb_full, sheetName = sheet)
  writeData(wb_full, sheet = sheet, full_tables[[sheet]])
}

combined_full <- bind_rows(full_tables)

addWorksheet(wb_full, sheetName = "All_Contrasts")
writeData(wb_full, sheet = "All_Contrasts", combined_full)

saveWorkbook(
  wb_full,
  file = "DEGs_by_contrast.xlsx",
  overwrite = TRUE
)
















