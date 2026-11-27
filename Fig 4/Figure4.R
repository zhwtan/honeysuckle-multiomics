# 1. Dependencies and paths-----------------------
input_dir <- "data"
output_dir <- "res"
required_packages <- c("data.table", "dplyr", "ggplot2", "mixOmics", "clusterProfiler", "UpSetR", "enrichplot", "igraph", "BioEnricher", "stringr", "pathview", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) stop("Install these R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(pathview)
})
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stage_labels <- paste0("S", 1:6)
set.seed(2025)
load(file.path(input_dir, "kegg_and_class_annotations.RData"))
# 2. Data preparation-----------------------
metabolome <- fread(file.path(input_dir, "metabolome_matrix.tsv.gz"), showProgress = FALSE)
stopifnot(nrow(metabolome) > 0L, !anyDuplicated(metabolome$Index))
metabolite_sample_columns <- grep("^(MS|WS)[1-6]-[1-3]$", names(metabolome), value = TRUE)
Gemh_standardize_sample_id <- function(x) {
  x <- gsub("-", "_", x, fixed = TRUE)
  sub("_(count|fpkm)$", "", x)
}
sample_ids <- sort(Gemh_standardize_sample_id(metabolite_sample_columns))
stopifnot(length(sample_ids) == 36L, !anyDuplicated(sample_ids))
sample_info <- data.table::data.table(
  Sample_ID = sample_ids,
  Cultivar = sub("^([A-Z]+)[1-6]_[1-3]$", "\\1", sample_ids),
  Stage = sub("^[A-Z]+([1-6])_[1-3]$", "S\\1", sample_ids),
  Replicate = as.integer(sub("^.*_([1-3])$", "\\1", sample_ids))
)
sample_info[, Cultivar := factor(Cultivar, levels = c("WS", "MS"))]
sample_info[, Stage := factor(Stage, levels = stage_labels)]
metabolite_matrix <- t(as.matrix(metabolome[, ..metabolite_sample_columns]))
rownames(metabolite_matrix) <- Gemh_standardize_sample_id(metabolite_sample_columns)
colnames(metabolite_matrix) <- metabolome$Index
metabolite_matrix <- metabolite_matrix[sample_ids, , drop = FALSE]
s14_columns <- c("Index", "Compounds", "Class I", "Class II", "Formula", "Level", "CAS",
                 sort(metabolite_sample_columns), "cpd_ID", "kegg_map", "mw_cpd_ID", "mw_map")
table_s14 <- metabolome[, ..s14_columns]
# 3. Differential abundance-----------------------
Gemh_stage_metabolite_tests <- function(i) {
  stage <- stage_labels[i]
  stage_samples <- sample_info[Stage == stage, Sample_ID]
  x <- metabolite_matrix[stage_samples, , drop = FALSE]
  group <- factor(sample_info[match(stage_samples, Sample_ID), Cultivar], levels = c("MS", "WS"))
  model <- suppressWarnings(mixOmics::plsda(x, group, ncomp = 3))
  vip <- suppressWarnings(mixOmics::vip(model)[, 1])
  variable <- colSums(sweep(x, 2, x[1, ], "-") != 0) > 0
  x <- x[, variable, drop = FALSE]
  tests <- vapply(seq_len(ncol(x)), function(j) {
    ms <- x[group == "MS", j]
    ws <- x[group == "WS", j]
    fc <- mean(ms) / mean(ws)
    c(FC = fc, log2FC = log2(fc), pvalue = stats::t.test(ms, ws)$p.value)
  }, numeric(3L))
  result <- as.data.table(t(tests))
  result[, .(ID = colnames(x), FC, log2FC, Group = "MS-WS", pvalue,
             FDR = stats::p.adjust(pvalue, "BH"), Stage = paste("Stage", i), VIP = vip[colnames(x)])]
}
s15 <- rbindlist(lapply(1:6, Gemh_stage_metabolite_tests))
setorder(s15, Stage, ID)
# 4. Shared metabolites-----------------------
stage_levels <- paste("Stage", 1:6)
dam_sets <- setNames(lapply(stage_levels, function(stage) {
  s15[Stage == stage & VIP > 1 & FDR < 0.05 & abs(log2FC) >= 0.58, ID]
}), stage_labels)
core_dam_ids <- Reduce(intersect, dam_sets)
kegg_term_to_gene <- as.data.frame(meatb_info_KEGG[, 1:2])
kegg_term_to_name <- as.data.frame(meatb_info_KEGG[, c(1, 3)])
kegg_universe <- intersect(s15$ID, kegg_term_to_gene$Index)
new_enrichment <- clusterProfiler::enricher(
  gene = intersect(core_dam_ids, kegg_term_to_gene$Index), universe = kegg_universe,
  TERM2GENE = kegg_term_to_gene, TERM2NAME = kegg_term_to_name, pvalueCutoff = 1, qvalueCutoff = 1
)
new_s16 <- as.data.frame(new_enrichment)[, c(
  "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
)]
names(new_s16)[5] <- "p value"
class_term_to_gene <- as.data.frame(meatb_info_ClassII[, 2:1])
class_enrichment <- clusterProfiler::enricher(
  gene = core_dam_ids, universe = unique(s15$ID), TERM2GENE = class_term_to_gene,
  pvalueCutoff = 1, qvalueCutoff = 1
)
class_results <- as.data.frame(class_enrichment)
significant_classes <- class_results |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(p.adjust, pvalue)
cluster_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500")
# 4.1. Set intersections-----------------------
upset_input <- UpSetR::fromList(dam_sets)
grDevices::pdf(paste0(output_dir, "/Figure_4A.pdf"), width = 5.5, height = 3.8, onefile = FALSE)
UpSetR::upset(
  upset_input, order.by = "freq", nsets = 6, nintersects = 30, matrix.color = "black",
  main.bar.color = "#377EB8", shade.color = "pink", point.size = 1.5, line.size = 0.6,
  sets.bar.color = cluster_colors, mainbar.y.label = "DEMs numbers", sets.x.label = "DEMs numbers"
)
grDevices::dev.off()
# 4.2. Chemical classes-----------------------
class_enrichment_plot <- clusterProfiler::enricher(
  gene = core_dam_ids, universe = unique(s15$ID), TERM2GENE = class_term_to_gene,
  pvalueCutoff = 0.05
)
figure_4b <- enrichplot::cnetplot(
  class_enrichment_plot, layout = igraph::layout_with_kk, color_edge = "category"
) +
  ggplot2::scale_size_continuous(name = "Size", range = c(3, 8)) +
  ggplot2::guides(
    size = ggplot2::guide_legend(override.aes = list(colour_ggnewscale_1 = "black"))
  ) +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 22),
    legend.title = ggplot2::element_text(size = 22)
  ) + ggplot2::labs(colour_ggnewscale_1 = "Category")
ggplot2::ggsave(
  filename = paste0(output_dir, "/Figure_4B.pdf"), plot = figure_4b,
  width = 13, height = 10, units = "in"
)
# 4.3. Pathway enrichment-----------------------
figure_4c <- BioEnricher::lzq_ORA.dotplot1(enrich.obj = new_enrichment, color.by = "pvalue") +
  ggplot2::labs(x = "Metabolite Ratio", y = "Pathways") +
  ggplot2::theme(axis.title.y = ggplot2::element_text()) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
  ggplot2::scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 30),
    position = "left"
  )
ggplot2::ggsave(
  filename = paste0(output_dir, "/Figure_4C.pdf"), plot = figure_4c,
  width = 6, height = 6, units = "in"
)
# 5. Pathway mapping-----------------------
compound_fc <- merge(s15[FC != 1 & ID %in% meatb_info_KEGG[pathway == "ko00941", Index]],
                     table_s14[, .(ID = Index, cpd_ID)], by = "ID")
compound_fc <- compound_fc[, .(Average_FC = mean(FC)), by = cpd_ID]
compound_values <- setNames(log2(compound_fc$Average_FC + 1), compound_fc$cpd_ID)
pathview_dir <- tempfile("pathview-", tmpdir = output_dir)
dir.create(pathview_dir)
old_dir <- getwd()
stopifnot(all(file.copy(file.path(input_dir, c("ko00941.png", "ko00941.xml")), pathview_dir)))
tryCatch({
  setwd(pathview_dir)
  invisible(pathview::pathview(cpd.data = compound_values, pathway.id = "00941", species = "ko",
                               cpd.idtype = "kegg", kegg.dir = ".", min.nnodes = 1,
                               low = list(gene = "#0a9396", cpd = "blue"),
                               mid = list(gene = "#e9d8a6", cpd = "gray"),
                               high = list(gene = "#bb3e03", cpd = "yellow")))
  setwd(old_dir)
  pathview_png <- file.path(pathview_dir, "ko00941.pathview.png")
  stopifnot(file.exists(pathview_png))
  stopifnot(file.copy(pathview_png, file.path(output_dir, "Figure_4D.png"), overwrite = TRUE))
}, finally = {
  setwd(old_dir)
  unlink(pathview_dir, recursive = TRUE)
})
# 6. Results-----------------------
writexl::write_xlsx(list(S14 = table_s14, S15 = s15, S16 = new_s16),
                    file.path(output_dir, "Tables_S14-S16.xlsx"), format_headers = FALSE)
save(s15, new_s16, core_dam_ids, significant_classes, compound_values, file = file.path(output_dir, "Figure_4.RData"))
