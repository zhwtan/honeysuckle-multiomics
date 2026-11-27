# UpSet plot 图A，同源基因簇（OrthoFinder）
# 输入文件：Orthogroups.GeneCount.tsv （OrthoFinder标准输出）
# 输出：upset_gene_family.pdf
library(UpSetR)
library(dplyr)

setwd("/Users/Users/Desktop/Project/金银花基因组/data/2_upset")

# 读取OrthoFinder基因家族计数表
og_count <- read.delim("Orthogroups.GeneCount.tsv", header = T, sep = "\t", row.names = 1)
# 去掉最后一列Total
og_count <- og_count[, -ncol(og_count)]
# 转为0/1矩阵：>0=存在(1)，0=不存在(0)
og_bin <- og_count %>% mutate(across(everything(), ~ifelse(.x > 0, 1, 0)))

# 转成UpSetR需要的输入格式
upset_input <- og_bin
# A.trichopoda, A.thaliana, V.vinifera, S.indicum, C.canephora, S.tuberosum, L.sativa, L.japonica, Z.mays, O.sativa
sp_order <- c("A.trichopoda","A.thaliana","V.vinifera","S.indicum","C.canephora","S.tuberosum","L.sativa","L.japonica","Z.mays","O.sativa")
upset_input <- upset_input[, sp_order]

pdf("upset_gene_family.pdf", width = 14, height = 7)
upset(upset_input,
      sets = sp_order,
      keep.order = TRUE,
      sets.bar.color = "black",
      order.by = "freq",
      decreasing = TRUE,
      show.numbers = TRUE,
      number.angle = 0,
      mainbar.y.label = "Gene family number",
      sets.x.label = "Gene family number",
      text.scale = c(1, 1, 1, 1, 1, 0.8))
dev.off()
