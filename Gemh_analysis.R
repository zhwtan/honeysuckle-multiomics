# 1. Dependencies and paths-----------------------
required_packages <- c("data.table", "mixOmics", "clusterProfiler", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}
set.seed(2025)
transcriptome_files <- paste0("data/transcriptome_matrix_part", 1:2, ".tsv.gz")
metabolome_file <- "data/metabolome_matrix.tsv.gz"
results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- paste0(results_dir, "/figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
cultivar_colors <- c(WS = "#272E6A", MS = "#D51F26")
stage_colors <- c("#D51F26", "#F47D2B", "#FEE500", "#89C75F", "#2F8AC4", "#6A3E98")
font_family <- if (requireNamespace("systemfonts", quietly = TRUE) &&
                   any(systemfonts::system_fonts()$family == "Arial")) "Arial" else "sans"

Gemh_standardize_sample_id <- function(x) {
  x <- gsub("-", "_", x, fixed = TRUE)
  sub("_(count|fpkm)$", "", x)
}

Gemh_run_ora <- function(genes, universe, term_to_gene, term_to_name) {
  genes <- intersect(unique(genes), unique(term_to_gene$Gene_ID))
  universe <- intersect(unique(universe), unique(term_to_gene$Gene_ID))
  if (length(genes) == 0 || length(universe) == 0) return(data.table::data.table())
  result <- clusterProfiler::enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = term_to_gene[, c("Term_ID", "Gene_ID")],
    TERM2NAME = term_to_name[, c("Term_ID", "Term_Name")],
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (is.null(result)) return(data.table::data.table())
  data.table::as.data.table(as.data.frame(result))
}

Gemh_cross_cor <- function(x, y) {
  x_centered <- sweep(x, 2, colMeans(x), FUN = "-")
  y_centered <- sweep(y, 2, colMeans(y), FUN = "-")
  x_ss <- sqrt(colSums(x_centered^2))
  y_ss <- sqrt(colSums(y_centered^2))
  if (any(x_ss == 0) || any(y_ss == 0)) stop("Zero-variance features were detected.")
  crossprod(x_centered, y_centered) / outer(x_ss, y_ss)
}

Gemh_theme_publication <- function() {
  ggplot2::theme_bw(base_size = 14, base_family = font_family) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16, color = "black"),
      axis.text = ggplot2::element_text(size = 14, color = "black"),
      legend.title = ggplot2::element_text(size = 14, color = "black"),
      legend.text = ggplot2::element_text(size = 14, color = "black"),
      strip.text = ggplot2::element_text(size = 14, color = "black", face = "bold"),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
}

Gemh_save_pdf <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = paste0(figure_dir, "/", filename),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "transparent"
  )
}

Gemh_wrap_label <- function(x, width = 32) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

Gemh_p_stars <- function(p) {
  ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# 2. Data import and sample information-----------------------
stopifnot(all(file.exists(transcriptome_files)))
transcriptome <- data.table::rbindlist(lapply(
  transcriptome_files,
  data.table::fread,
  showProgress = FALSE
))
metabolome <- data.table::fread(metabolome_file, showProgress = FALSE)
stopifnot(nrow(transcriptome) == 37608L, nrow(metabolome) == 1326L)
data.table::fwrite(
  metabolome,
  file = paste0(results_dir, "/Table_S14_all_metabolites.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

count_columns <- grep("_count$", names(transcriptome), value = TRUE)
fpkm_columns <- grep("_fpkm$", names(transcriptome), value = TRUE)
metabolite_sample_columns <- grep("^(MS|WS)[1-6]-[1-3]$", names(metabolome), value = TRUE)
stopifnot(length(count_columns) == 36L, length(fpkm_columns) == 36L, length(metabolite_sample_columns) == 36L)

sample_ids <- sort(Gemh_standardize_sample_id(fpkm_columns))
stopifnot(setequal(sample_ids, Gemh_standardize_sample_id(count_columns)))
stopifnot(setequal(sample_ids, Gemh_standardize_sample_id(metabolite_sample_columns)))
sample_info <- data.table::data.table(
  Sample_ID = sample_ids,
  Cultivar = sub("^([A-Z]+)[1-6]_[1-3]$", "\\1", sample_ids),
  Stage = sub("^[A-Z]+([1-6])_[1-3]$", "S\\1", sample_ids),
  Replicate = as.integer(sub("^.*_([1-3])$", "\\1", sample_ids))
)
sample_info[, Cultivar := factor(Cultivar, levels = c("WS", "MS"))]
sample_info[, Stage := factor(Stage, levels = paste0("S", 1:6))]
data.table::fwrite(sample_info, file = paste0(results_dir, "/sample_information.tsv"), sep = "\t")

fpkm_matrix <- as.matrix(transcriptome[, ..fpkm_columns])
rownames(fpkm_matrix) <- transcriptome$gene_id
colnames(fpkm_matrix) <- Gemh_standardize_sample_id(fpkm_columns)
fpkm_matrix <- fpkm_matrix[, sample_ids, drop = FALSE]

metabolite_matrix <- t(as.matrix(metabolome[, ..metabolite_sample_columns]))
rownames(metabolite_matrix) <- Gemh_standardize_sample_id(metabolite_sample_columns)
colnames(metabolite_matrix) <- metabolome$Index
metabolite_matrix <- metabolite_matrix[sample_ids, , drop = FALSE]

# 3. Transcriptome analysis-----------------------
transcriptome_contrasts <- c(
  paste0("MS", 2:6, "vsMS", 1:5),
  paste0("WS", 2:6, "vsWS", 1:5),
  paste0("MS", 1:6, "vsWS", 1:6)
)

Gemh_extract_de_result <- function(contrast) {
  required_columns <- paste0(
    contrast,
    c("_log2FoldChange", "_pvalue", "_padj")
  )
  stopifnot(all(required_columns %in% names(transcriptome)))
  data.table::data.table(
    Gene_ID = transcriptome$gene_id,
    Contrast = contrast,
    log2FC = transcriptome[[required_columns[1]]],
    P_value = transcriptome[[required_columns[2]]],
    BH_FDR = transcriptome[[required_columns[3]]]
  )
}

transcriptome_de <- data.table::rbindlist(lapply(transcriptome_contrasts, Gemh_extract_de_result))
transcriptome_de[, Differential := !is.na(BH_FDR) & BH_FDR < 0.05 & abs(log2FC) > 1]
data.table::fwrite(
  transcriptome_de,
  file = paste0(results_dir, "/transcriptome_differential_expression.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

stage_matched_contrasts <- paste0("MS", 1:6, "vsWS", 1:6)
core_up_ids <- Reduce(
  intersect,
  lapply(stage_matched_contrasts, function(x) {
    transcriptome_de[Contrast == x & !is.na(BH_FDR) & BH_FDR < 0.05 & log2FC > 1, Gene_ID]
  })
)
core_down_ids <- Reduce(
  intersect,
  lapply(stage_matched_contrasts, function(x) {
    transcriptome_de[Contrast == x & !is.na(BH_FDR) & BH_FDR < 0.05 & log2FC < -1, Gene_ID]
  })
)
core_deg_ids <- c(core_up_ids, core_down_ids)
stopifnot(length(core_deg_ids) == 1248L)
core_degs <- unique(transcriptome[gene_id %in% core_deg_ids, .(Gene_ID = gene_id, gene_name, gene_description, KO_ID)])
core_degs[, Direction := ifelse(Gene_ID %in% core_up_ids, "Up", "Down")]
data.table::fwrite(core_degs, file = paste0(results_dir, "/core_differential_genes.tsv"), sep = "\t")

fpkm_sample_by_gene <- t(fpkm_matrix)
expressed_genes <- colMeans(fpkm_sample_by_gene > 1) > 0.1
gene_sd <- apply(fpkm_sample_by_gene[, expressed_genes, drop = FALSE], 2, stats::sd)
top_variable_genes <- names(sort(gene_sd, decreasing = TRUE))[seq_len(min(1000L, length(gene_sd)))]
pca_model <- stats::prcomp(fpkm_sample_by_gene[, top_variable_genes, drop = FALSE], scale. = TRUE)
pca_scores <- data.table::as.data.table(pca_model$x[, 1:5, drop = FALSE], keep.rownames = "Sample_ID")
pca_scores <- merge(sample_info, pca_scores, by = "Sample_ID", sort = FALSE)
data.table::fwrite(pca_scores, file = paste0(results_dir, "/transcriptome_pca_scores.tsv"), sep = "\t")

gene_pathway <- transcriptome[!is.na(KEGG_Pathways) & KEGG_Pathways != "", .(
  Pathway = unlist(strsplit(KEGG_Pathways, ";", fixed = TRUE))
), by = .(Gene_ID = gene_id)]
gene_pathway[, c("Term_ID", "Term_Name") := data.table::tstrsplit(Pathway, "::", fixed = TRUE)]
gene_term_to_gene <- unique(gene_pathway[, .(Term_ID, Gene_ID)])
gene_term_to_name <- unique(gene_pathway[, .(Term_ID, Term_Name)])
gene_universe <- unique(gene_term_to_gene$Gene_ID)
transcriptome_enrichment <- data.table::rbindlist(lapply(transcriptome_contrasts, function(contrast) {
  genes <- transcriptome_de[Contrast == contrast & Differential == TRUE, Gene_ID]
  result <- Gemh_run_ora(genes, gene_universe, gene_term_to_gene, gene_term_to_name)
  if (nrow(result) > 0) result[, Contrast := contrast]
  result
}), fill = TRUE)
data.table::fwrite(
  transcriptome_enrichment,
  file = paste0(results_dir, "/Table_S13_transcriptome_KEGG_enrichment.tsv"),
  sep = "\t"
)

# 4. Metabolome analysis-----------------------
Gemh_metabolite_de <- function(stage_number) {
  stage_label <- paste0("S", stage_number)
  stage_samples <- sample_info[Stage == stage_label, Sample_ID]
  x <- metabolite_matrix[stage_samples, , drop = FALSE]
  group <- sample_info[match(stage_samples, Sample_ID), Cultivar]
  group <- factor(group, levels = c("MS", "WS"))
  plsda_model <- suppressWarnings(mixOmics::plsda(x, group, ncomp = 3))
  vip_score <- suppressWarnings(mixOmics::vip(plsda_model)[, 1])
  test_result <- vapply(seq_len(ncol(x)), function(j) {
    ms_value <- x[group == "MS", j]
    ws_value <- x[group == "WS", j]
    fold_change <- mean(ms_value) / mean(ws_value)
    p_value <- if (length(unique(c(ms_value, ws_value))) > 1) {
      stats::t.test(ms_value, ws_value)$p.value
    } else {
      1
    }
    c(Fold_Change = fold_change, log2FC = log2(fold_change), P_value = p_value)
  }, numeric(3))
  result <- data.table::as.data.table(t(test_result))
  result[, Metabolite_ID := colnames(x)]
  result[, `:=`(
    Stage = stage_label,
    BH_FDR = stats::p.adjust(P_value, method = "BH"),
    VIP = vip_score[Metabolite_ID]
  )]
  score <- data.table::data.table(
    Sample_ID = stage_samples,
    Cultivar = group,
    Stage = stage_label,
    Component_1 = plsda_model$variates$X[, 1],
    Component_2 = plsda_model$variates$X[, 2]
  )
  list(
    differential = result[, .(Metabolite_ID, Stage, Fold_Change, log2FC, P_value, BH_FDR, VIP)],
    scores = score
  )
}

metabolite_stage_results <- lapply(1:6, Gemh_metabolite_de)
metabolite_de <- data.table::rbindlist(lapply(metabolite_stage_results, `[[`, "differential"))
plsda_scores <- data.table::rbindlist(lapply(metabolite_stage_results, `[[`, "scores"))
metabolite_de[, Differential := P_value < 0.05 & VIP > 1]
data.table::fwrite(
  metabolite_de,
  file = paste0(results_dir, "/Table_S15_metabolite_differential_abundance.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(plsda_scores, file = paste0(results_dir, "/metabolite_PLSDA_scores.tsv"), sep = "\t")

core_dam_ids <- Reduce(
  intersect,
  lapply(paste0("S", 1:6), function(x) metabolite_de[Stage == x & Differential == TRUE, Metabolite_ID])
)
stopifnot(length(core_dam_ids) == 259L)
metabolite_annotation_columns <- intersect(
  c("Index", "Compounds", "Class I", "Class II", "Formula", "Level", "CAS", "cpd_ID", "kegg_map"),
  names(metabolome)
)
core_dams <- unique(metabolome[Index %in% core_dam_ids, ..metabolite_annotation_columns])
data.table::fwrite(core_dams, file = paste0(results_dir, "/core_differential_metabolites.tsv"), sep = "\t")

metabolite_pathway <- metabolome[!is.na(kegg_map) & kegg_map != "-", .(
  Term_ID = trimws(unlist(strsplit(kegg_map, ",", fixed = TRUE)))
), by = .(Gene_ID = Index)]
metabolite_pathway <- unique(metabolite_pathway[grepl("^(ko|map)[0-9]{5}$", Term_ID)])
metabolite_pathway[, Term_ID := sub("^map", "ko", Term_ID)]
metabolite_term_to_name <- unique(gene_term_to_name[Term_ID %in% metabolite_pathway$Term_ID])
missing_term_ids <- setdiff(unique(metabolite_pathway$Term_ID), metabolite_term_to_name$Term_ID)
if (length(missing_term_ids) > 0) {
  metabolite_term_to_name <- data.table::rbindlist(list(
    metabolite_term_to_name,
    data.table::data.table(Term_ID = missing_term_ids, Term_Name = missing_term_ids)
  ))
}
metabolite_enrichment <- Gemh_run_ora(
  core_dam_ids,
  metabolome$Index,
  metabolite_pathway,
  metabolite_term_to_name
)
data.table::fwrite(
  metabolite_enrichment,
  file = paste0(results_dir, "/Table_S16_metabolite_KEGG_enrichment.tsv"),
  sep = "\t"
)

class_term_to_gene <- unique(metabolome[
  !is.na(get("Class II")) & get("Class II") != "-",
  .(Term_ID = get("Class II"), Gene_ID = Index)
])
class_term_to_name <- unique(class_term_to_gene[, .(Term_ID, Term_Name = Term_ID)])
class_enrichment <- Gemh_run_ora(
  core_dam_ids,
  metabolome$Index,
  class_term_to_gene,
  class_term_to_name
)
data.table::fwrite(
  class_enrichment,
  file = paste0(results_dir, "/core_metabolite_Class_II_enrichment.tsv"),
  sep = "\t"
)

peonidin_ids <- metabolome[grepl("peonidin", Compounds, ignore.case = TRUE) & as.character(Level) == "3", Index]
s6_samples <- sample_info[Stage == "S6", Sample_ID]
s6_group <- sample_info[match(s6_samples, Sample_ID), Cultivar]
peonidin_statistics <- data.table::rbindlist(lapply(peonidin_ids, function(id) {
  value <- metabolite_matrix[s6_samples, id]
  data.table::data.table(
    Metabolite_ID = id,
    Compound = metabolome[Index == id, Compounds][1],
    Fold_Change = mean(value[s6_group == "MS"]) / mean(value[s6_group == "WS"]),
    P_value = stats::t.test(value[s6_group == "MS"], value[s6_group == "WS"])$p.value
  )
}))
peonidin_sum <- rowSums(metabolite_matrix[s6_samples, peonidin_ids, drop = FALSE])
peonidin_statistics <- data.table::rbindlist(list(
  peonidin_statistics,
  data.table::data.table(
    Metabolite_ID = "Four_peonidin_metabolites",
    Compound = "Sum of four level 3 peonidin-related metabolites",
    Fold_Change = mean(peonidin_sum[s6_group == "MS"]) / mean(peonidin_sum[s6_group == "WS"]),
    P_value = stats::t.test(peonidin_sum[s6_group == "MS"], peonidin_sum[s6_group == "WS"])$p.value
  )
))
data.table::fwrite(peonidin_statistics, file = paste0(results_dir, "/peonidin_S6_statistics.tsv"), sep = "\t")

# 5. Transcriptome-metabolome integration-----------------------
expression_core <- t(fpkm_matrix[core_deg_ids, sample_ids, drop = FALSE])
metabolite_core <- metabolite_matrix[sample_ids, core_dam_ids, drop = FALSE]
stopifnot(nrow(expression_core) == 36L, ncol(expression_core) == 1248L)
stopifnot(nrow(metabolite_core) == 36L, ncol(metabolite_core) == 259L)

pearson_r <- Gemh_cross_cor(expression_core, metabolite_core)
pearson_df <- nrow(expression_core) - 2
pearson_t <- pearson_r * sqrt(pearson_df / pmax(1 - pearson_r^2, .Machine$double.eps))
pearson_p <- 2 * stats::pt(-abs(pearson_t), df = pearson_df)

design_matrix <- stats::model.matrix(~ Cultivar + Stage, data = sample_info)
design_qr <- qr(design_matrix)
expression_residual <- qr.resid(design_qr, expression_core)
metabolite_residual <- qr.resid(design_qr, metabolite_core)
partial_r <- Gemh_cross_cor(expression_residual, metabolite_residual)
partial_df <- nrow(expression_core) - design_qr$rank - 1
partial_t <- partial_r * sqrt(partial_df / pmax(1 - partial_r^2, .Machine$double.eps))
partial_p <- 2 * stats::pt(-abs(partial_t), df = partial_df)

correlation_results <- data.table::data.table(
  Gene_ID = rep(colnames(expression_core), times = ncol(metabolite_core)),
  Metabolite_ID = rep(colnames(metabolite_core), each = ncol(expression_core)),
  Pearson_r = as.vector(pearson_r),
  Pearson_P_value = as.vector(pearson_p),
  Pearson_BH_FDR = stats::p.adjust(as.vector(pearson_p), method = "BH"),
  Partial_Pearson_r = as.vector(partial_r),
  Partial_P_value = as.vector(partial_p),
  Partial_BH_FDR = stats::p.adjust(as.vector(partial_p), method = "BH")
)
stopifnot(nrow(correlation_results) == 323232L)
data.table::fwrite(
  correlation_results,
  file = paste0(results_dir, "/Table_S17_gene_metabolite_correlations.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)

candidate_map <- data.table::data.table(
  Gene_ID = c("Lj02G265400", "Lj02G307900"),
  Gene_Symbol = c("LjAN2", "LjGST")
)
candidate_correlations <- merge(correlation_results, candidate_map, by = "Gene_ID")
candidate_correlations <- merge(
  candidate_correlations,
  core_dams,
  by.x = "Metabolite_ID",
  by.y = "Index",
  all.x = TRUE,
  sort = FALSE
)
candidate_correlations[, Pearson_Selected := abs(Pearson_r) > 0.7 & Pearson_BH_FDR < 0.05]
data.table::fwrite(
  candidate_correlations,
  file = paste0(results_dir, "/LjAN2_LjGST_metabolite_correlations.tsv"),
  sep = "\t"
)

# 6. Figures-----------------------
pca_variance <- round(100 * pca_model$sdev^2 / sum(pca_model$sdev^2), 1)
figure_3b <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(PC1, PC2, color = Stage, shape = Cultivar)
) +
  ggplot2::geom_point(size = 3.2, alpha = 0.9) +
  ggplot2::scale_color_manual(values = stage_colors) +
  ggplot2::scale_shape_manual(values = c(WS = 16, MS = 17)) +
  ggplot2::labs(
    x = paste0("PC1 (", pca_variance[1], "%)"),
    y = paste0("PC2 (", pca_variance[2], "%)")
  ) +
  Gemh_theme_publication()
Gemh_save_pdf(figure_3b, "Figure_3B_transcriptome_PCA.pdf", 7.2, 5.4)

deg_counts <- transcriptome_de[
  Contrast %in% stage_matched_contrasts & Differential == TRUE,
  .N,
  by = .(Contrast, Direction = ifelse(log2FC > 0, "Up", "Down"))
]
deg_counts[, Stage := factor(paste0("S", sub("^MS([1-6]).*$", "\\1", Contrast)), levels = paste0("S", 1:6))]
data.table::fwrite(deg_counts, file = paste0(results_dir, "/Figure_3C_DEG_counts.tsv"), sep = "\t")
figure_3c <- ggplot2::ggplot(deg_counts, ggplot2::aes(Stage, N, fill = Direction)) +
  ggplot2::geom_col(position = "stack", width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = N),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 4,
    color = "white",
    family = font_family
  ) +
  ggplot2::scale_fill_manual(values = c(Up = "#D51F26", Down = "#2F8AC4")) +
  ggplot2::labs(x = "Developmental stage", y = "Number of DEGs") +
  Gemh_theme_publication()
Gemh_save_pdf(figure_3c, "Figure_3C_stage_matched_DEG_counts.pdf", 7.2, 5.4)

enrichment_plot_data <- transcriptome_enrichment[
  Contrast %in% stage_matched_contrasts & !is.na(p.adjust) & p.adjust < 0.05
]
data.table::setorder(enrichment_plot_data, Contrast, p.adjust)
enrichment_plot_data <- enrichment_plot_data[, head(.SD, 8), by = Contrast]
enrichment_plot_data[, Stage := factor(paste0("S", sub("^MS([1-6]).*$", "\\1", Contrast)), levels = paste0("S", 1:6))]
enrichment_plot_data[, Pathway := Gemh_wrap_label(Description)]
enrichment_plot_data[, Pathway := factor(Pathway, rev(unique(Pathway)))]
figure_3d <- ggplot2::ggplot(
  enrichment_plot_data,
  ggplot2::aes(Stage, Pathway, size = Count, color = p.adjust)
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#D51F26", high = "#2F8AC4", trans = "reverse") +
  ggplot2::scale_size_continuous(range = c(2, 7)) +
  ggplot2::labs(x = "Developmental stage", y = "KEGG pathway", color = "BH-FDR", size = "Gene count") +
  Gemh_theme_publication()
Gemh_save_pdf(figure_3d, "Figure_3D_stage_matched_KEGG_enrichment.pdf", 10.5, 8.5)

figure_s3 <- ggplot2::ggplot(
  enrichment_plot_data,
  ggplot2::aes(-log10(p.adjust), stats::reorder(Pathway, -p.adjust), size = Count, color = p.adjust)
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::facet_wrap(~Stage, scales = "free_y", ncol = 2) +
  ggplot2::scale_color_gradient(low = "#D51F26", high = "#2F8AC4", trans = "reverse") +
  ggplot2::scale_size_continuous(range = c(2, 7)) +
  ggplot2::labs(x = expression(-log[10]("BH-FDR")), y = "KEGG pathway", color = "BH-FDR", size = "Gene count") +
  Gemh_theme_publication()
Gemh_save_pdf(figure_s3, "Figure_S3_stage_specific_KEGG_enrichment.pdf", 14, 16)

volcano_data <- transcriptome_de[Contrast %in% stage_matched_contrasts & !is.na(P_value)]
volcano_data[, Status := ifelse(Differential & log2FC > 0, "Up", ifelse(Differential, "Down", "Not significant"))]
volcano_data[, `:=`(
  Stage = factor(paste0("S", sub("^MS([1-6]).*$", "\\1", Contrast)), levels = paste0("S", 1:6)),
  Minus_log10_P = -log10(pmax(P_value, .Machine$double.xmin))
)]
volcano_plot <- ggplot2::ggplot(volcano_data, ggplot2::aes(log2FC, Minus_log10_P, color = Status)) +
  ggplot2::geom_point(size = 0.45, alpha = 0.5) +
  ggplot2::facet_wrap(~Stage, ncol = 3) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
  ggplot2::scale_color_manual(values = c(Up = "#D51F26", Down = "#2F8AC4", `Not significant` = "grey75")) +
  ggplot2::labs(x = expression(log[2]("fold change")), y = expression(-log[10](P)), color = NULL) +
  Gemh_theme_publication()
Gemh_save_pdf(volcano_plot, "Supplementary_stage_matched_DEG_volcano.pdf", 12, 8)

plsda_plot <- ggplot2::ggplot(
  plsda_scores,
  ggplot2::aes(Component_1, Component_2, color = Cultivar, shape = Cultivar)
) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::facet_wrap(~Stage, scales = "free", ncol = 3) +
  ggplot2::scale_color_manual(values = cultivar_colors) +
  ggplot2::scale_shape_manual(values = c(WS = 16, MS = 17)) +
  ggplot2::labs(x = "PLS-DA component 1", y = "PLS-DA component 2") +
  Gemh_theme_publication()
Gemh_save_pdf(plsda_plot, "Supplementary_metabolite_PLSDA_scores.pdf", 12, 8)

dam_sets <- setNames(
  lapply(paste0("S", 1:6), function(x) metabolite_de[Stage == x & Differential == TRUE, Metabolite_ID]),
  paste0("S", 1:6)
)
dam_membership <- data.table::data.table(Metabolite_ID = unique(unlist(dam_sets)))
for (stage in names(dam_sets)) dam_membership[, (stage) := Metabolite_ID %in% dam_sets[[stage]]]
dam_membership[, Pattern := do.call(paste0, lapply(.SD, as.integer)), .SDcols = names(dam_sets)]
intersection_counts <- dam_membership[, .N, by = Pattern][order(-N)][seq_len(min(30L, .N))]
intersection_counts[, Intersection := factor(seq_len(.N), levels = seq_len(.N))]
intersection_matrix <- unique(dam_membership[Pattern %in% intersection_counts$Pattern, c("Pattern", names(dam_sets)), with = FALSE])
intersection_matrix <- merge(intersection_matrix, intersection_counts[, .(Pattern, Intersection)], by = "Pattern")
intersection_matrix <- data.table::melt(
  intersection_matrix,
  id.vars = c("Pattern", "Intersection"),
  variable.name = "Stage",
  value.name = "Present"
)
intersection_matrix[, Stage := factor(Stage, levels = rev(names(dam_sets)))]
intersection_matrix[, Stage_Index := as.numeric(Stage)]
intersection_segments <- intersection_matrix[Present == TRUE, .(
  Stage_min = min(Stage_Index),
  Stage_max = max(Stage_Index)
), by = Intersection]
set_sizes <- data.table::data.table(
  Stage = factor(rev(names(dam_sets)), levels = rev(names(dam_sets))),
  N = lengths(dam_sets)[rev(names(dam_sets))]
)
upset_top <- ggplot2::ggplot(intersection_counts, ggplot2::aes(Intersection, N)) +
  ggplot2::geom_col(fill = "#2F8AC4", width = 0.72) +
  ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 3.6, family = font_family) +
  ggplot2::labs(x = NULL, y = "Intersection size") +
  ggplot2::coord_cartesian(clip = "off") +
  Gemh_theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
upset_matrix <- ggplot2::ggplot(intersection_matrix, ggplot2::aes(Intersection, Stage_Index)) +
  ggplot2::geom_segment(
    data = intersection_segments,
    ggplot2::aes(Intersection, Stage_min, xend = Intersection, yend = Stage_max),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(ggplot2::aes(color = Present), size = 2.7) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") +
  ggplot2::scale_y_continuous(breaks = seq_along(levels(intersection_matrix$Stage)), labels = levels(intersection_matrix$Stage)) +
  ggplot2::labs(x = NULL, y = NULL) +
  Gemh_theme_publication() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
upset_sets <- ggplot2::ggplot(set_sizes, ggplot2::aes(-N, Stage, fill = Stage)) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = setNames(stage_colors, paste0("S", 1:6))) +
  ggplot2::scale_x_continuous(labels = abs) +
  ggplot2::labs(x = "Set size", y = NULL) +
  Gemh_theme_publication()
grDevices::cairo_pdf(
  filename = paste0(figure_dir, "/Figure_4A_core_metabolite_UpSet.pdf"),
  width = 9,
  height = 6,
  family = font_family,
  bg = "transparent"
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(
  nrow = 2,
  ncol = 2,
  widths = grid::unit(c(3, 9), "null"),
  heights = grid::unit(c(6, 4), "null")
)))
print(upset_top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2), newpage = FALSE)
print(upset_sets, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
print(upset_matrix, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2), newpage = FALSE)
grDevices::dev.off()

class_counts <- core_dams[!is.na(get("Class II")) & get("Class II") != "-", .N, by = .(Class_II = get("Class II"))]
data.table::setorder(class_counts, -N)
data.table::fwrite(class_counts, file = paste0(results_dir, "/Figure_4B_core_metabolite_classification.tsv"), sep = "\t")
class_plot_data <- head(class_counts, 15)
class_plot_data[, Class_II := factor(Class_II, rev(Class_II))]
figure_4b <- ggplot2::ggplot(class_plot_data, ggplot2::aes(N, Class_II)) +
  ggplot2::geom_col(fill = "#2F8AC4", width = 0.72) +
  ggplot2::labs(x = "Number of core differential metabolites", y = "Class II") +
  Gemh_theme_publication()
Gemh_save_pdf(figure_4b, "Figure_4B_core_metabolite_classification.pdf", 9, 7)

metabolite_enrichment_plot <- data.table::copy(metabolite_enrichment)
data.table::setorder(metabolite_enrichment_plot, pvalue)
metabolite_enrichment_plot <- head(metabolite_enrichment_plot, 15)
metabolite_enrichment_plot[, Pathway := factor(Gemh_wrap_label(Description), rev(Gemh_wrap_label(Description)))]
figure_4c <- ggplot2::ggplot(
  metabolite_enrichment_plot,
  ggplot2::aes(-log10(pvalue), Pathway, size = Count, color = pvalue)
) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_gradient(low = "#D51F26", high = "#2F8AC4", trans = "reverse") +
  ggplot2::scale_size_continuous(range = c(2.5, 8)) +
  ggplot2::labs(x = expression(-log[10](P)), y = "KEGG pathway", color = "P value", size = "Metabolite count") +
  Gemh_theme_publication()
Gemh_save_pdf(figure_4c, "Figure_4C_metabolite_KEGG_enrichment.pdf", 9, 7)

flavonoid_ids <- unique(metabolite_pathway[Term_ID == "ko00941", Gene_ID])
flavonoid_heatmap <- merge(
  metabolite_de[Metabolite_ID %in% intersect(core_dam_ids, flavonoid_ids)],
  metabolome[, .(Metabolite_ID = Index, Compound = Compounds)],
  by = "Metabolite_ID"
)
flavonoid_heatmap[, `:=`(
  Stage = factor(Stage, levels = paste0("S", 1:6)),
  Significance = Gemh_p_stars(P_value),
  Compound = factor(Compound, rev(unique(Compound)))
)]
data.table::fwrite(flavonoid_heatmap, file = paste0(results_dir, "/Figure_4D_flavonoid_heatmap_data.tsv"), sep = "\t")
figure_4d <- ggplot2::ggplot(flavonoid_heatmap, ggplot2::aes(Stage, Compound, fill = log2FC)) +
  ggplot2::geom_tile(color = "black", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = Significance), size = 4, family = font_family) +
  ggplot2::scale_fill_gradient2(low = "#2F8AC4", mid = "white", high = "#D51F26", midpoint = 0) +
  ggplot2::labs(x = "Developmental stage", y = NULL, fill = expression(log[2]("MS/WS"))) +
  Gemh_theme_publication()
Gemh_save_pdf(figure_4d, "Figure_4D_flavonoid_pathway_heatmap.pdf", 8.5, 5.5)

workflow_nodes <- data.table::data.table(
  x = c(2, 8, 2, 8, 8, 5, 5, 5),
  y = c(8, 8, 6.5, 6.5, 5, 3.5, 2, 0.5),
  label = c(
    "Transcriptome\n6 MS vs WS contrasts", "Metabolome\n6 MS vs WS contrasts",
    "1,248 core DEGs\nBH-FDR < 0.05; |log2FC| > 1", "259 core DAMs\nVIP > 1; P < 0.05",
    "Flavonoid biosynthesis\nKEGG ko00941", "Pearson correlation\n|r| > 0.7; BH-FDR < 0.05",
    "Candidate genes\nLjAN2, LjGST, RPS2, ENO, RPL44", "Partial Pearson correlation\nAdjusted for cultivar and stage"
  )
)
workflow_edges <- data.table::data.table(
  x = c(2, 8, 8, 2, 8, 5, 5),
  y = c(7.65, 7.65, 6.15, 6.15, 4.65, 3.15, 1.65),
  xend = c(2, 8, 8, 4.5, 5.5, 5, 5),
  yend = c(6.85, 6.85, 5.35, 3.85, 3.85, 2.35, 0.85)
)
figure_s4 <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = workflow_edges,
    ggplot2::aes(x, y, xend = xend, yend = yend),
    linewidth = 0.8,
    arrow = grid::arrow(length = grid::unit(0.18, "cm"))
  ) +
  ggplot2::geom_label(
    data = workflow_nodes,
    ggplot2::aes(x, y, label = label),
    size = 4.4,
    family = font_family,
    linewidth = 0.5,
    fill = "white"
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 10), ylim = c(0, 8.7), clip = "off") +
  ggplot2::theme_void(base_family = font_family) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent", color = NA))
Gemh_save_pdf(figure_s4, "Figure_S4_candidate_gene_screening_workflow.pdf", 10, 11)

all_candidate_map <- data.table::rbindlist(list(
  candidate_map,
  data.table::data.table(
    Gene_ID = c("Lj05G088500", "Lj07G096000", "Lj03G332700"),
    Gene_Symbol = c("RPS2", "ENO", "RPL44")
  )
))
candidate_expression <- data.table::as.data.table(
  t(fpkm_matrix[all_candidate_map$Gene_ID, sample_ids, drop = FALSE]),
  keep.rownames = "Sample_ID"
)
candidate_expression <- data.table::melt(
  candidate_expression,
  id.vars = "Sample_ID",
  variable.name = "Gene_ID",
  value.name = "FPKM"
)
candidate_expression <- merge(candidate_expression, sample_info, by = "Sample_ID")
candidate_expression <- merge(candidate_expression, all_candidate_map, by = "Gene_ID")
candidate_expression[, log2_FPKM := log2(FPKM + 1)]
expression_statistics <- candidate_expression[, .(
  P_value = stats::t.test(log2_FPKM[Cultivar == "MS"], log2_FPKM[Cultivar == "WS"])$p.value,
  Label_Y = max(log2_FPKM) * 1.08
), by = .(Gene_Symbol, Stage)]
expression_statistics[, Significance := Gemh_p_stars(P_value)]
data.table::fwrite(expression_statistics, file = paste0(results_dir, "/candidate_gene_expression_statistics.tsv"), sep = "\t")

Gemh_candidate_expression_plot <- function(symbols) {
  input_data <- candidate_expression[Gene_Symbol %in% symbols]
  input_stats <- expression_statistics[Gene_Symbol %in% symbols]
  ggplot2::ggplot(input_data, ggplot2::aes(Stage, log2_FPKM, fill = Cultivar)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(Stage, Cultivar)),
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Cultivar, group = Cultivar),
      position = ggplot2::position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
      size = 2
    ) +
    ggplot2::geom_text(
      data = input_stats,
      ggplot2::aes(Stage, Label_Y, label = Significance),
      inherit.aes = FALSE,
      size = 4.5,
      family = font_family
    ) +
    ggplot2::facet_wrap(~Gene_Symbol, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = cultivar_colors) +
    ggplot2::scale_color_manual(values = cultivar_colors) +
    ggplot2::labs(x = "Developmental stage", y = expression(log[2]("FPKM + 1"))) +
    Gemh_theme_publication()
}

figure_s5ab <- Gemh_candidate_expression_plot(c("LjAN2", "LjGST"))
Gemh_save_pdf(figure_s5ab, "Figure_S5AB_LjAN2_LjGST_expression.pdf", 10, 5.5)
figure_s6 <- Gemh_candidate_expression_plot(c("RPS2", "ENO", "RPL44"))
Gemh_save_pdf(figure_s6, "Figure_S6_RPS2_ENO_RPL44_expression.pdf", 14, 5.5)

network_edges <- candidate_correlations[
  Metabolite_ID %in% flavonoid_ids & Pearson_Selected == TRUE,
  .(Gene_Symbol, Metabolite_ID, Compound = Compounds, Pearson_r, Pearson_BH_FDR)
]
gene_nodes <- data.table::data.table(Gene_Symbol = c("LjAN2", "LjGST"), Gene_x = 0, Gene_y = c(0.7, 0.3))
metabolite_nodes <- unique(network_edges[, .(Metabolite_ID, Compound)])
metabolite_nodes[, `:=`(Metabolite_x = 1, Metabolite_y = seq(0.25, 0.75, length.out = .N))]
network_edges <- merge(network_edges, gene_nodes, by = "Gene_Symbol")
network_edges <- merge(network_edges, metabolite_nodes, by = c("Metabolite_ID", "Compound"))
figure_s5d <- ggplot2::ggplot() +
  ggplot2::geom_curve(
    data = network_edges,
    ggplot2::aes(Gene_x, Gene_y, xend = Metabolite_x, yend = Metabolite_y, linewidth = abs(Pearson_r)),
    color = "#D51F26",
    curvature = 0.15,
    alpha = 0.75
  ) +
  ggplot2::geom_point(data = gene_nodes, ggplot2::aes(Gene_x, Gene_y), size = 8, color = "#2F8AC4") +
  ggplot2::geom_point(
    data = metabolite_nodes,
    ggplot2::aes(Metabolite_x, Metabolite_y),
    size = 8,
    color = "#D51F26"
  ) +
  ggplot2::geom_text(
    data = gene_nodes,
    ggplot2::aes(Gene_x - 0.05, Gene_y, label = Gene_Symbol),
    hjust = 1,
    size = 4.5,
    family = font_family
  ) +
  ggplot2::geom_text(
    data = metabolite_nodes,
    ggplot2::aes(Metabolite_x + 0.05, Metabolite_y, label = Compound),
    hjust = 0,
    size = 4.5,
    family = font_family
  ) +
  ggplot2::scale_linewidth(range = c(0.8, 2), guide = "none") +
  ggplot2::coord_cartesian(xlim = c(-0.35, 1.7), ylim = c(0, 1), clip = "off") +
  ggplot2::theme_void(base_family = font_family) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent", color = NA))
Gemh_save_pdf(figure_s5d, "Figure_S5D_LjAN2_LjGST_flavonoid_network.pdf", 10, 5.5)

# 7. Reproducibility record-----------------------
run_summary <- data.table::data.table(
  Item = c("Samples", "Core_DEGs", "Core_DAMs", "Gene_metabolite_tests"),
  Value = c(nrow(sample_info), length(core_deg_ids), length(core_dam_ids), nrow(correlation_results))
)
data.table::fwrite(run_summary, file = paste0(results_dir, "/run_summary.tsv"), sep = "\t")
capture.output(utils::sessionInfo(), file = paste0(results_dir, "/sessionInfo.txt"))
cat("Analysis completed: 36 samples, 1,248 core DEGs, 259 core DAMs, and 323,232 correlations.\n")
