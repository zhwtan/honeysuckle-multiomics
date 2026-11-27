# GO气泡图（C） + KEGG条形图（D）
# 输入：go_enrich_result.tsv , kegg_enrich_result.tsv
# 输出：GO_bubble.pdf , KEGG_bar.pdf
library(ggplot2)
library(dplyr)
library(stringr)

setwd("/Users/Users/Desktop/Project/金银花基因组/data/3.GO.KEGG")

# ========== 图C GO气泡图 ==========
go <- read.delim("go_enrich_result.tsv", sep = "\t", header = T)
# 筛选BP/CC/MF，x=GeneRatio，y=Description，size=Count，color=-log10(p.adjust)
go <- go %>%
  mutate(
    GeneRatio = Count / BgRatio,
    negLogP = -log10(p.adjust)
  )
# 挑选top term（和示例图一样只画展示的这些）
go_plot <- go %>% arrange(desc(negLogP)) %>% slice(1:10)

pdf("GO_bubble.pdf", width = 8, height = 7)
ggplot(go_plot, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = negLogP)) +
  scale_color_gradient(low = "#87CEEB", high = "#EE6363") +
  scale_size_continuous(range = c(3, 10)) +
  labs(x = "Gene ratio", y = "GO terms", color = "-log(p value)", size = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

# ========== 图D KEGG横向条形图 ==========
kegg <- read.delim("kegg_enrich_result.tsv", sep = "\t", header = T)
kegg <- kegg %>%
  mutate(negLogP = -log10(p.adjust)) %>%
  arrange(negLogP)
# 取top10通路
kegg_plot <- kegg %>% slice_tail(n=10)

pdf("KEGG_bar.pdf", width = 10, height = 6)
ggplot(kegg_plot, aes(x = Count, y = reorder(Description, Count), fill = negLogP)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#800080", high = "#B22222") +
  labs(x = "Gene numbers", y = "Pathways", fill = "p.adj") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()
