# 1. Dependencies and paths-----------------------
input_dir <- "data"
output_dir <- "res"
required_packages <- c("data.table", "dplyr", "ggplot2", "stringr", "clusterProfiler", "qvalue", "scRNAtoolVis", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) stop("Install these R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stage_labels <- paste0("S", 1:6)
load(file.path(input_dir, "kegg_and_class_annotations.RData"))
# 2. Data preparation-----------------------
transcriptome_files <- file.path(input_dir, paste0("transcriptome_matrix_part", 1:2, ".tsv.gz"))
transcriptome <- rbindlist(lapply(transcriptome_files, fread, showProgress = FALSE))
stopifnot(nrow(transcriptome) > 0L, !anyDuplicated(transcriptome$gene_id))
fpkm_columns <- grep("_fpkm$", names(transcriptome), value = TRUE)
Gemh_standardize_sample_id <- function(x) {
  x <- gsub("-", "_", x, fixed = TRUE)
  sub("_(count|fpkm)$", "", x)
}
sample_ids <- sort(Gemh_standardize_sample_id(fpkm_columns))
stopifnot(length(sample_ids) == 36L, !anyDuplicated(sample_ids))
sample_info <- data.table::data.table(
  Sample_ID = sample_ids,
  Cultivar = sub("^([A-Z]+)[1-6]_[1-3]$", "\\1", sample_ids),
  Stage = sub("^[A-Z]+([1-6])_[1-3]$", "S\\1", sample_ids),
  Replicate = as.integer(sub("^.*_([1-3])$", "\\1", sample_ids))
)
sample_info[, Cultivar := factor(Cultivar, levels = c("WS", "MS"))]
sample_info[, Stage := factor(Stage, levels = stage_labels)]
fpkm_matrix <- as.matrix(transcriptome[, ..fpkm_columns])
rownames(fpkm_matrix) <- transcriptome$gene_id
colnames(fpkm_matrix) <- Gemh_standardize_sample_id(fpkm_columns)
fpkm_matrix <- fpkm_matrix[, sample_ids, drop = FALSE]
transcriptome_contrasts <- c(
  paste0("MS", 2:6, "vsMS", 1:5),
  paste0("WS", 2:6, "vsWS", 1:5),
  paste0("MS", 1:6, "vsWS", 1:6)
)
stage_matched_contrasts <- paste0("MS", 1:6, "vsWS", 1:6)
expr_degs <- rbindlist(lapply(transcriptome_contrasts, function(contrast) {
  columns <- paste0(contrast, c(paste0("_", sub("vs.*$", "", contrast)), "_log2FoldChange", "_padj"))
  stopifnot(all(columns %in% names(transcriptome)))
  data.table(gene_id = transcriptome$gene_id, Case = transcriptome[[columns[1]]],
             log2FC = transcriptome[[columns[2]]], p.adj = transcriptome[[columns[3]]], group = contrast)
}))
# 3. Expression variation-----------------------
fpkm_sample_by_gene <- t(fpkm_matrix)
expressed_genes <- colMeans(fpkm_sample_by_gene > 1) > 0.1
gene_sd <- apply(fpkm_sample_by_gene[, expressed_genes, drop = FALSE], 2, stats::sd)
top_variable_genes <- names(sort(gene_sd, decreasing = TRUE))[seq_len(1000L)]
pca_model <- stats::prcomp(fpkm_sample_by_gene[, top_variable_genes, drop = FALSE], scale. = TRUE)
pca_data <- merge(sample_info, as.data.table(pca_model$x[, 1:2], keep.rownames = "Sample_ID"),
                  by = "Sample_ID", sort = FALSE)
pca_data[, Cultivar := factor(Cultivar, levels = c("MS", "WS"))]
pca_variance <- 100 * pca_model$sdev^2 / sum(pca_model$sdev^2)
stage_colors <- c("#332288", "#88CCEE", "#44AA99", "#DDCC77", "#882255", "#CC6677")
figure_3b <- ggplot(pca_data, aes(PC1, PC2, color = Stage, shape = Cultivar)) +
  geom_point(size = 3) + theme_bw() + scale_color_manual(values = stage_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_variance[1]), y = sprintf("PC2 (%.1f%%)", pca_variance[2]), shape = "Mutation Type")
ggsave(file.path(output_dir, "Figure_3B.pdf"), figure_3b, width = 5.5, height = 4.2)
# 4. Differential expression-----------------------
stage_degs <- expr_degs[group %in% stage_matched_contrasts & !is.na(Case)]
deg_counts <- stage_degs[p.adj < 0.05 & abs(log2FC) > 1,
  .(Up = sum(log2FC > 1), Down = sum(log2FC < -1)), by = group]
volcano <- stage_degs[, .(gene = gene_id, p_val = p.adj, avg_log2FC = log2FC,
                         cluster = sub("vs", " vs ", group, fixed = TRUE))]
volcano[, p_val_adj := stats::p.adjust(p_val, "BH")]
figure_3c <- scRNAtoolVis::jjVolcano(
  diffData = as.data.frame(volcano), pvalue.cutoff = 0.05, log2FC.cutoff = 1, col.type = "updown", topGeneN = 3,
  aesCol = c("#0072B2", "#A65B14"), tile.col = c("#ABB7DF", "#9FD3EA", "#A2D6CA", "#F2DE9B", "#D7A9CA", "#E5A98F"),
  celltypeSize = 3.6
) +
  labs(x = "Stage Comparison", y = "Average log2 fold change") + theme(legend.position = "top")
stage_label_data <- unique(as.data.frame(volcano[, .(cluster)]))
for (i in seq_along(figure_3c$layers)) {
  if (inherits(figure_3c$layers[[i]]$geom, c("GeomTile", "GeomText"))) {
    figure_3c$layers[[i]]$data <- stage_label_data
  }
}
set.seed(20261112)
ggsave(file.path(output_dir, "Figure_3C.pdf"), figure_3c, width = 12, height = 8)
# 5. Pathway enrichment-----------------------
s13_columns <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
transcriptome_enrichment <- rbindlist(lapply(transcriptome_contrasts, function(contrast) {
  differential <- expr_degs[group == contrast & !is.na(Case) & p.adj < 0.05 & abs(log2FC) > 1]
  ids <- unique(inner_join(differential[, .(gene.id = gene_id)], ko.id[, .(gene.id, KOID)], by = "gene.id")$KOID)
  stage_matched <- contrast %in% stage_matched_contrasts
  result <- clusterProfiler::enricher(
    gene = ids, universe = unique(ko.id$KOID),
    TERM2GENE = ko_pathway.all[, .(pathway, ko)], TERM2NAME = ko_pathway.all[, .(pathway, Description)],
    pvalueCutoff = 1, qvalueCutoff = 1
  )
  out <- as.data.table(as.data.frame(result))[, ..s13_columns]
  out[, qvalue := qvalue::qvalue(pvalue, lambda = 0.05, pi0.method = "bootstrap")$qvalues]
  out[, geneID := vapply(strsplit(geneID, "/", fixed = TRUE),
                        function(x) paste(intersect(ids, x), collapse = "/"), character(1))]
  if (!stage_matched) out <- out[qvalue <= 0.2]
  out[, group := contrast]
  out
}))
excluded <- "05171"
excluded_matched <- c(excluded, "04936", "05417", "05418")
transcriptome_enrichment <- transcriptome_enrichment[!sub("^(ko|map)", "", ID) %in% excluded]
matched_enrichment <- transcriptome_enrichment[
  group %in% stage_matched_contrasts & !sub("^(ko|map)", "", ID) %in% excluded_matched
]
selected_pathways <- matched_enrichment[p.adjust < 0.05] |>
  arrange(group, p.adjust) |> group_by(group) |> slice_head(n = 8)
pathway_data <- filter(matched_enrichment, ID %in% unique(selected_pathways$ID) & p.adjust < 0.2) %>%
  dplyr::select(ID,Description,p.adjust,Count,group) %>%
  mutate(num = str_extract(group,"\\d+"),
         Mutation.type = str_extract(group,"MS|WS"))
pathway_data <- arrange(pathway_data,Mutation.type,num)
pathway_data <- mutate(pathway_data,group = factor(group,unique(pathway_data$group)))
figure_3d <- ggplot(pathway_data, aes(group, Description)) +
  geom_point(aes(size = Count,color = p.adjust)) +
  scale_color_gradientn(colours =  c("#003c30",
                                     "#01665e", "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3",
                                     "#dfc27d", "#bf812d", "#8c510a", "#543005")) +
  theme_bw(base_rect_size = 1.5) +
  labs(x = "", y = "KEGG pathways")  +
  scale_y_discrete(labels=function(x) str_wrap(x, width=30),position = "left") +
  scale_x_discrete(position = "top",label = paste0("S",1:6)) +
  theme(axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 12))
ggsave(file.path(output_dir, "Figure_3D.pdf"), figure_3d, width = 6.5, height = 5)
table_s13 <- rbind(
  transcriptome_enrichment[!group %in% stage_matched_contrasts], matched_enrichment
)[startsWith(ID, "ko")]
setorder(table_s13, pvalue)
# 6. Results-----------------------
writexl::write_xlsx(list(S13 = table_s13), file.path(output_dir, "Table_S13.xlsx"), format_headers = FALSE)
save(pca_data, pca_variance, deg_counts, table_s13, pathway_data, file = file.path(output_dir, "Figure_3.RData"))
