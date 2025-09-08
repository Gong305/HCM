setwd("~/scRNA-heart-mitochodria")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure5/311/"
dir.create(output, recursive = TRUE)
########巨噬细胞5.5------
library(readxl)
metacell.obj <- qs::qread("~/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Macrophage")
sub_metacell.obj@meta.data
deg <- read_excel("~/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg1 <- subset(deg,subset = deg$`Cell Type` %in% "Macrophage")
deg1$change <- ""
deg1$change[deg1$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg1$change[deg1$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg1$change[deg1$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg1 <- deg1[,c(1,15,16,17,21)]
colnames(deg1) <- c("Gene","logFC","Pvalue","Adjusted_Pvalue","change")
deg1$log10Pvalue <- log10(deg1$Pvalue)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
down <- subset(deg1,subset = change %in% "down") 
down <- subset(down,subset = down$Gene %in% dbs$gene)
library(ggplot2)
library(dplyr)
library(ggrepel)
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
down_deg <-
  down[order(down$logFC, decreasing = F),] %>% pull("Gene")
down_bp <-
  enrichGO(
    down_deg,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
down_term <- down_bp@result
down_term$celltype <- "Macrophage"


write.table(
  down_term,
  file = "~/scRNA-heart-mitochodria/results/Macrophage/function/HCM-NF/down.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(down_bp, "~/scRNA-heart-mitochodria/results/Macrophage/function/HCM-NF/down.rds")
######GSEA网络图-----
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
install.packages("aPEAR")
library(aPEAR)
Idents(sub_metacell.obj) <- "group"
markers <- FindMarkers(sub_metacell.obj,
                       ident.1 ="HCM",
                       ident.2 = "NF",
                       min_pct = 0,
                       logfc.threshold = 0)
write.csv(markers,file = paste0(output,"gsea_deg.csv"))
markers$gene <- rownames(markers)
down_deg = markers[,2,drop = T]
names(down_deg) = as.character(markers[,6,drop = T])
down_deg = sort(down_deg, decreasing = TRUE)
enrich <- gseGO(down_deg, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",minGSSize = 1,maxGSSize = 50000)
term <- enrich@result
term <- subset(term,subset = term$Description %in% c("oxidative phosphorylation","proton motive force-driven ATP synthesis","aerobic respiration","carboxylic acid transport","aerobic electron transport chain","MHC class II protein complex assembly","MHC protein complex assembly","immune response-regulating signaling pathway","immune response-activating signaling pathway","immune response","regulation of immune system process"))
enrichmentNetwork(term,repelLabels = T, drawEllipses = TRUE)
p2
ggsave(paste0(output,"GO富集柱状图上调.pdf"), plot = p2, width = 4.45, height = 4.05)

####细胞比例的桑吉图-----
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
Ratio <- metacell.obj@meta.data %>%
  group_by(group, cell_type) %>% # 分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c('#efb306',
            '#7db954',
            '#852f88',
            '#4e54ac',
            '#0f8096',
            'pink',
            'green')

ggplot(Ratio, aes(x =group, y= relative_freq, fill = cell_type,
                  stratum=cell_type, alluvium=cell_type)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='group',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
####细胞比例的变化箱型图-----
####m6A酶的箱型图-----
library(Seurat)
library(tidyverse)
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
source("/home/xingwl/share/20220802_scrna-m6A/custom_plot_function.R")
library(tidyverse)

Ratio <- metacell.obj@meta.data %>%
  group_by(sample, cell_type) %>% # 分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
sub_Ratio <- subset(Ratio,Ratio$cell_type %in% "Macrophage")

sub_Ratio$group <- NA
sub_Ratio$group[grep("hcm",sub_Ratio$sample)] <- "HCM"
sub_Ratio$group[grep("nf",sub_Ratio$sample)] <- "NF"
sub_Ratio$group <- factor(sub_Ratio$group, levels = c("NF", "HCM"))
library(ggpubr)
f <- ggplot(sub_Ratio, aes(x = group, y = relative_freq)) +
  geom_boxplot(
    aes(fill = NA, color = group), # 边框颜色
    alpha = 0 # 填充透明
  ) +
  geom_jitter(
    aes(color = group), 
    width = 0.2, 
    size = 1.5
  ) +
  scale_fill_manual(values = c(NF = "#2e86c1", HCM = "#e32d32")) +
  scale_color_manual(values = c(NF = "#2e86c1", HCM = "#e32d32")) + # 边框颜色
  theme_classic() +
  NoLegend() +   
  ylim(0, 0.25) +
  theme(
    axis.text.x = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Macrophages") +
  # geom_signif(
  #   comparisons = list(c("NF", "HCM")), # 选择你想在哪2组上添加标签
  #   test = "wilcox.test", # "t 检验，比较两组（参数）" = "t.test","wilcox.test” 符号秩检验，比较两组（非参数）" = "wilcox.test"
  #   map_signif_level = F # 标签样式F为数字，T为*号
  # ) +
  xlab("Group") +
  ylab("Proportion of macrophages")

# 打印图形
print(f)
ggsave(
  file.path(output, "巨噬比例箱型图.pdf"),
  plot = f,
  height = 3.5,
  width = 3
)
####细胞比例与线粒体功能的相关性-----
meta <- sub_metacell.obj@meta.data
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
metascore.obj1<- AddModuleScore_UCell(sub_metacell.obj, features = dbs1,
                                      ncores = 20)
meta <- metascore.obj1@meta.data
colnames(meta)
meta1 <- meta[,c(1,99,120,62,63,144,94)]
meta1 <- aggregate(meta1, by=list(biosample=meta1$sample),mean)##求均值
meta1$sample <- meta1$biosample
adata1 <- merge(meta1,sub_Ratio,by = "sample")

cor_function_meta <- corr.test(adata1[,c(11)],adata1[,c(3:8)], method = "spearman", adjust = "fdr")
cp <- as.data.frame(cor_function_meta$p)
cr <- as.data.frame(cor_function_meta$r)
cr <- round(cr,2)

cp$names <- rownames(cp)
longp <- melt(cp,idvar = "names",v.names = "abd",direction = "long")
colnames(longp) <- c('na','va','p')

cr$names <- rownames(cr)
longr <- melt(cr,idvar = "names",v.names = "abd",direction = "long")
allnew <- as.data.frame(cbind(longr,longp$p))
colnames(allnew) <- c('feature_x','feature_y','r','p')


filtered_adata <- allnew %>%
  filter(p <= 0.05)

###相关性散点图----
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
# adata <- cbind(adata,t_results2)
# colnames(adata) <- gsub(" ","_",colnames(adata))
###去除0值需要两条通路单独来做-----
# data <- adata[,c("sample","Fatty_acid_oxidation_UCell","CELL_MIGRATION_INVOLVED_IN_HEART_DEVELOPMENT_UCell")]
# data1 <- data
# data1[data1 == 0] <- NA
# data1 <- na.omit(data1)
# # data1 <- aggregate(data1, by=list(biosample=data1$sample),mean)##求均值
# colnames(data1)
######相关性散点图带山峦的图-----
colnames(adata1) <- gsub(" ","_",colnames(adata1))
P <- adata1 %>%
  ggplot(aes(x = TCA_cycle_UCell, y = relative_freq)) +
  geom_point(size = 0.9,color="black",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = TCA_cycle_UCell, y = relative_freq),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Relative freq")+
  NoLegend()
library(ggExtra)
p12 <- ggMarginal(P, type = "density", 
                  xparams = list(fill = "pink"),
                  yparams = list(fill = "#7db954"))
p12
ggsave(paste0(output,"TCA和巨噬细胞比例相关性散点图.pdf"), plot = p12, width = 5.39, height = 3.35)
#####UCell打分热图展示------
########UCELL打分----柱状图可视化-----
library(UCell)
metascore.obj <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
metascore.obj <- subset(metascore.obj,subset = cell_type %in% "Macrophage")
meta<- metascore.obj@meta.data
meta <- meta[, c("group",
                 colnames(meta)[grep("_UCell", colnames(meta))])]   #提取出含GOBP的
meta <-
  aggregate(meta[,-1], list(meta$group), FUN = mean)
rownames(meta) <- meta[, 1]
meta <- meta[, -1]
meta <- t(meta)

#### 可视化 ----
library(ComplexHeatmap)
meta <- t(scale(t(meta), center = T, scale = T))

heat <- Heatmap(
  meta,
  cluster_rows = T,
  cluster_columns = T,
  col = colorRampPalette(c("#2f73bb", "white", "#e94644"))(10),
  heatmap_legend_param = list(grid_height = unit(10, 'mm')),
  show_row_names = T,
  show_column_dend = F,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 4),
  border = T,
  name = "Score",
  heatmap_width = unit(5, "cm") * ncol(meta),
  heatmap_height = unit(0.15, "cm") * nrow(meta),
  use_raster = F
)
heat

###############################################################
####功能通路差异分析 ----
library(limma)
meta<- metascore.obj@meta.data
es.matrix <- meta[,c(colnames(meta)[grep("_UCell", colnames(meta))])]

meta <-
  metascore.obj@meta.data[, c("group", "cell_type")]
meta$cell_id <- rownames(meta)
es.matrix <- t(es.matrix)
es.matrix <- es.matrix[,meta$cell_id]
clusters <- names(table(meta$cell_type))
library(purrr)
library(tidyverse)
case <- "HCM"
control <- "NF"
cluster <- "Macrophage"
t_results <- map_dfr(clusters, function(cluster) {
  cell_id.1 <- meta[meta$cell_type == cluster &
                      meta$group == case,]$cell_id
  cell_id.2 <- meta[meta$cell_type == cluster &
                      meta$group == control,]$cell_id
  es.matrix.1 <-
    as.data.frame(es.matrix[, cell_id.1],
                  row.names = row.names(es.matrix))
  es.matrix.2 <-
    as.data.frame(es.matrix[, cell_id.2],
                  row.names = row.names(es.matrix))
  
  es.matrix.f <- cbind(es.matrix.1, es.matrix.2)
  
  # 分组设计
  grouP <-
    c(rep("case", dim(es.matrix.1)[2]),
      rep("control", dim(es.matrix.2)[2])) %>% as.factor()
  
  desigN <- model.matrix( ~ grouP + 0)
  rownames(desigN) <-
    c(colnames(es.matrix.1), colnames(es.matrix.2))
  # desigN
  comparE <-
    makeContrasts(grouPcase - grouPcontrol, levels = desigN)
  fiT <- lmFit(es.matrix.f, desigN)
  fiT2 <- contrasts.fit(fiT, comparE)
  fiT3 <- eBayes(fiT2)
  diff <- topTable(fiT3, coef = 1, number = 200)
  temp_res <- diff %>% rownames_to_column("feature") %>%
    dplyr::select(feature, t) %>%   #报错
    mutate(cluster = cluster)
  return(temp_res)
})
t_results <- t_results[order(t_results[,2],decreasing=T),]
rownames(t_results) <- t_results[,1]
ggplot(t_results,aes(x=feature,y = t))+
  geom_bar(stat = "identity")+
  theme_classic()+
  ylim(-40,40)+
  coord_flip()


t_results$fill <- "no"
t_results$fill[t_results$t > 2.58 ] <- "up"
t_results$fill[t_results$t < -2.58 ] <- "down"

t_results1 <- subset(t_results,subset = feature %in% c("Creatine.metabolism_UCell","Complex.II_UCell","Complex.III_UCell","Complex.V_UCell","Complex.IV_UCell","Complex.I_UCell","Fatty.acid.oxidation_UCell","Mitophagy_UCell"))
library(ggplot2)
p3 <-
  ggplot(t_results1, aes(x = reorder(feature,t), y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#36648b",
      "up" = "#e94644",
      "no" = "#cccccc"
    ),
    guide = "none"
  ) +
  geom_hline(
    yintercept = c(-2.58, 2.58),
    colour = "white",
    linetype = "dashed",
    size = 0.5
  ) +
  coord_flip() +
  xlab("Macrophage  (HCM/NF)") +
  ylab("t value of UCell score") +
  geom_text(
    data = subset(t_results1, t < 0),
    aes(
      x = feature,
      y = 0.1,
      label = paste0(" ", feature),
      colour = "black"
    ),
    size = 2.5,
    hjust = "outward"
  ) +
  geom_text(
    data = subset(t_results1, t > 0),
    aes(
      x = feature,
      y = -0.1,
      label = paste0(" ", feature),
      colour = "black"
    ),
    size = 2.5,
    hjust = "inward"
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "black"),
                      guide = "none") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
p3
ggsave(paste0(output,"线粒体功能差异.pdf"), plot = p3, width = 7.25, height = 4.35)

library(UCell)
dbs2 <- read_excel("~/scRNA-heart-mitochodria/data/巨噬细胞下游功能通路.xlsx")
colnames(dbs2) <- c("pathway","gene")
dbs2 <- split(dbs2$gene, dbs2$pathway)
metascore.obj <- AddModuleScore_UCell(metascore.obj, features = dbs2,
                                    ncores = 50)

#####线粒体相关通路与炎症相关基因的相关性------

library(Seurat)
library(tidyverse)
library(psych)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
outdir <- output

# Ready to analyze cell type?
selected_cell_type <- "Macrophage"
###################   1 相关性部分可视化----
#### 气泡图展示关注基因与通路的相关性 ----
data <- metascore.obj@assays$RNA@data
data <- as.data.frame(data)
data1 <- data[c(dbs2[["GOBP_CELLULAR_HOMEOSTASIS"]],dbs2[["GOBP_INFLAMMATORY_RESPONSE"]]),]
data1 <- na.omit(data1)
data1 <- as.data.frame(t(data1))

adata <- metascore.obj@meta.data
adata1 <- cbind(adata,data1)
colnames(adata1)
feature1 <- c("A2M","AOAH","ASH1L","BLNK","CASP1","CASP4","CD200R1","CD28","CIITA","CSF1R","LYN") #inflammatory
feature2 <- c("VDR","UNC13B","TXNRD1","TXN","TTC7A","TRPM4","TCIRG1","SOD1","SLC39A8","CCR1","CCR2") #CELLULAR_HOMEOSTASIS
cor_function_meta <- corr.test(adata1[,c(62,63,94,99,120,144)],adata1[,c(feature1,feature2)], method = "spearman", adjust = "fdr")
cp <- as.data.frame(cor_function_meta$p)
cr <- as.data.frame(cor_function_meta$r)
cr <- round(cr,2)

cp$names <- rownames(cp)
longp <- melt(cp,idvar = "names",v.names = "abd",direction = "long")
colnames(longp) <- c('na','va','p')

cr$names <- rownames(cr)
longr <- melt(cr,idvar = "names",v.names = "abd",direction = "long")
allnew <- as.data.frame(cbind(longr,longp$p))
colnames(allnew) <- c('feature_x','feature_y','r','p')

filtered_adata <- allnew %>%
  filter(p <= 0.05)
filtered_adata$feature_y <- gsub("_UCell", "", filtered_adata$feature_y )
filtered_adata$feature_y <- gsub("GOBP_", "", filtered_adata$feature_y )
filtered_adata$feature_y  <- gsub("_", " ", filtered_adata$feature_y )
filtered_adata$feature_y  <- tolower(filtered_adata$feature_y )
filtered_adata$feature_x <- gsub("_UCell", "", filtered_adata$feature_x )

p <- filtered_adata %>%
  # filter(feature_x %in% selected_genes,
  #        feature_y %in% selected_features) %>%
  # dplyr:: mutate(
  #   feature_x = fct_relevel(feature_x, selected_genes),
  #   feature_y = fct_relevel(feature_y, selected_features)
  # ) %>%
  catdotplot(
    x = feature_x,
    y = feature_y,
    size = -log10(p),
    color = r,
    dot_scale = 7,
    title = selected_cell_type
  ) +
  coord_fixed() +
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text.x = element_text(face = "italic")
  )
p  
ggsave(
  file.path(output, "pathway_cor_dotplot.pdf"),
  plot = p,
  height = 8.43,
  width = 5.42,
)
#####通路与通路相关性散点图-----
expr <- metascore.obj@meta.data
expr <- as.data.frame(expr)
colnames(expr)
adata1 <- expr[,c(1,18:176)]
adata<- aggregate(adata1, by=list(sample=adata1$sample),mean)##求均值
rownames(adata) <- adata$sample
adata <- as.data.frame(adata)
filtered_adata <- adata|>
  as.data.frame()
filtered_adata <- filtered_adata[,-c(1,2)]
filtered_adata %>% mutate(across(where(is.character), as.numeric))  -> filtered_adata

###散点图----
#第1种画法
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
P <- filtered_adata %>%
  ggplot(aes(x = Mitophagy_UCell, y = GOBP_CELLULAR_HOMEOSTASIS_UCell)) +
  geom_point(size = 0.5,color="#7db954",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = Mitophagy_UCell, y = GOBP_CELLULAR_HOMEOSTASIS_UCell),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Cellular homeostasis")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"细胞稳态与线粒体自噬散点图.pdf"), plot = p4, width = 3.83, height = 3.35)

P <- filtered_adata %>%
  ggplot(aes(x = OXPHOS_UCell, y = GOBP_CELLULAR_HOMEOSTASIS_UCell)) +
  geom_point(size = 0.5,color="#7db954",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = OXPHOS_UCell, y = GOBP_CELLULAR_HOMEOSTASIS_UCell),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Cellular homeostasis")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"细胞稳态与OXPHOS散点图.pdf"), plot = p4, width = 3.83, height = 3.35)
#####相关性带显著程度的热图-通路与通路的相关性-----
library(vegan)
#install.packages('vegan')
library(psych)
library(reshape2)
library(ggplot2)
#data("varechem")

#### GSEA 
library(purrr)
deg_2$gene <- rownames(deg_2)
gsea_res <- deg_2 |>
  filter(group == "ComplexI") |>
  cat_gsea(category = "C5")
saveRDS(gsea_res, file.path(outdir, "ComplexI_C5_gsea.rds"))

TGF <- subset(TGFB,subset = TGFB$pathway %in% c("GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION","GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION"))
features <- TGF$gene
#####在HCM中组按照线粒体自噬分高低表达组-----
library(stringr)
meta <- metascore.obj@meta.data
selected_genes <- "Mitophagy_UCell"
median(meta[,selected_genes])
metascore.obj[[str_c(selected_genes, "_group")]] <-
  if_else(meta[,selected_genes] > median(meta[,selected_genes]),
          "high", "low")
table(metascore.obj$Mitophagy_UCell_group)
#######高低表达组的差异火山------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#差异分析--
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
### 高低表达组的差异火山图 ----
DefaultAssay(metascore.obj)
DefaultAssay(metascore.obj) <- "RNA"
Idents(metascore.obj) <- "Mitophagy_UCell_group"
markers <- FindMarkers(metascore.obj,
                       ident.1 ="high",
                       ident.2 = "low",
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
deg2 <- markers
deg2$group <- "Mitophagy_UCell"
deg2$change <- "not"
deg2$change[deg2$p_val <0.05 & deg2$avg_log2FC >0.25] <- "up"
deg2$change[deg2$p_val <0.05 & deg2$avg_log2FC < -0.25] <- "down"
deg2$log10Pvalue <- log10(deg2$p_val)
features <- c(feature1,feature2)
features <- c("LYN","CIITA","CASP4","TXNRD1")
a <- deg2[features,]
library(ggplot2)
library(dplyr)
library(ggrepel)
p1 <- ggplot(deg2, aes(x =avg_log2FC, y=-log10Pvalue, colour=change)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.9, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#7db954","grey","pink")) + xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围 #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = 2.853872, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="black",lwd=0.8)+
  labs(x="log2FC", y="-log10pvalue") +  #x、y轴标签
  ggtitle("Mitophagy-high vs Mitophagy-low") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()) 
P2 <- p1+ geom_label_repel(data = a,
                     aes(x =avg_log2FC, y=-log10Pvalue, label = features),
                     size = 2, fill="#CCFFFF",max.overlaps = 300,label.padding = 0.1)
ggsave(paste0(output,"Mitophagy火山图.pdf"), plot = P2, width = 5.29, height = 3.69)
write.csv(deg2,file = paste0(output,"Mitophagy_UCell差异基因.csv"))
###功能富集------
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
deg2$Gene <- rownames(deg2)
up_deg <- subset(deg2,subset = deg2$change %in% "up")
down_deg <- subset(deg2,subset = deg2$change %in% "down")
up_deg <-
  up_deg[order(up_deg$avg_log2FC, decreasing = T),] %>% pull("Gene")
down_deg <-
  down_deg[order(down_deg$avg_log2FC, decreasing = T),] %>% pull("Gene")
up_bp <-
  enrichGO(
    up_deg,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term1 <- up_bp@result
term1$celltype <- "Mitophagy_high"
down_bp <-
  enrichGO(
    down_deg,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term2 <- down_bp@result
term2$celltype <- "Mitophagy_low"
####功能富集结果----
#提取
subset_up_term <- subset(term1,subset = term1$Description %in% c("cell communication by electrical coupling","cellular carbohydrate metabolic process","negative regulation of interleukin-12 production","respiratory system development","organic hydroxy compound catabolic process","negative regulation of cytokine production","cell redox homeostasis"))
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
subset_up_term$text_x <- rep(0.03,7)
p6 <- ggplot(data = subset_up_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "pink", high = "pink") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p6
ggsave(paste0(output,"GO富集巨噬分组柱状图上调.pdf"), plot = p6, width = 4.45, height = 4.05) 
###下调---
subset_down_term <- subset(term2,subset = term2$Description %in% c("MHC class II protein complex assembly","MHC protein complex assembly","peptide antigen assembly with MHC class II protein complex","cell activation involved in immune response","calcium ion homeostasis","immune response-regulating cell surface receptor signaling pathway","positive regulation of interleukin-2 production"))
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
subset_down_term$text_x <- rep(0.03,7)
p6 <- ggplot(data = subset_down_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#7db954", high = "#7db954") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p6
ggsave(paste0(output,"GO富集巨噬分组柱状图下调.pdf"), plot = p6, width = 4.45, height = 4.05)
#### GSEA ----
Idents(metascore.obj) <- "Mitophagy_UCell_group"
markers1 <- FindMarkers(metascore.obj,
                       ident.1 ="high",
                       ident.2 = "low",
                       min.pct = 0,
                       logfc.threshold = 0)
deg1 <- markers
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(ggplot2)
gsea.input <-
  deg1[order(deg1$avg_log2FC, decreasing = T),]

category <- "C5"
genesets <- msigdbr(species = "Homo sapiens",#####Homo sapiens,Mus musculus
                    category = category)
genesets <- subset(genesets,
                   select = c("gs_name", "gene_symbol"))

genelist <-
  structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
genelist <- na.omit(genelist)
res <- GSEA(genelist, TERM2GENE = genesets, minGSSize = 1,maxGSSize = 5000,pvalueCutoff = 1)
library(enrichplot)
saveRDS(res, file.path(output, "Mitophagy_C5_gsea.rds"))
res <- readRDS(paste0(output, "Mitophagy_C5_gsea.rds"))
####可视化-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
a <- res@result
setid2 <- c("HP_MITOCHONDRIAL_INHERITANCE","GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN","GOBP_OXIDATIVE_PHOSPHORYLATION","GOCC_MHC_CLASS_II_PROTEIN_COMPLEX","GOBP_RESPONSE_TO_CYTOKINE")

p7 <- gseaNb(object = res,
       geneSetID = setid2,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.8,
       pvalY = 0.5,
       # pFill = "transparent",
       geneCol = "#009933",
       newCurveCol = c("#014f84","grey","#890102"))

p7
ggsave(paste0(output,"Mitophagy_GSEA.pdf"), plot = p7, width = 7.24, height = 5.13)
### 小提琴图-----
library(ggpubr)
my_comparisons <- c("low","high")
meta <- metascore.obj@meta.data
colnames(meta) <- gsub("_UCell|GOBP_", "", colnames(meta))
colnames(meta) <- tolower(colnames(meta))
p9 <- ggviolin(meta, x='mitophagy_group',y="cellular_homeostasis",
               xlab = 'Mitophagy group', ylab = "Cellular homeostasis",
               fill='mitophagy_group',palette = c("#7db954", "pink"),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$cellular_homeostasis),# 取得分的中值
             linetype = "dashed", size = 0.5, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("low","high")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"稳态小提琴.pdf"), plot = p9, width = 3.41, height = 5.92) 

library(msigdbr)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
dbs2 <- Dataset[grepl("CYTOKINE|MHC", Dataset$gs_name, ignore.case = TRUE), ]
dbs2 <- split(dbs2$gene_symbol, dbs2$gs_name)
metascore.obj <- AddModuleScore_UCell(metascore.obj, features = dbs2,
                                      ncores = 50)


meta <- metascore.obj@meta.data
colnames(meta) <- gsub("_UCell|GOBP_", "", colnames(meta))
colnames(meta) <- tolower(colnames(meta))
colnames(meta)
p9 <- ggviolin(meta, x='mitophagy_group',y="gocc_mhc_class_ii_protein_complex",
               xlab = 'Mitophagy group', ylab = "MHC class ii protein complex",
               fill='mitophagy_group',palette = c("#7db954", "pink"),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$gocc_mhc_class_ii_protein_complex),# 取得分的中值
             linetype = "dashed", size = 0.5, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("low","high")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"MHCII_小提琴.pdf"), plot = p9, width = 3.41, height = 5.92) 



#####在HCM中组按照OXPHOS分高低表达组-----
library(stringr)
meta <- metascore.obj1@meta.data
selected_genes <- "OXPHOS_UCell"
median(meta[,selected_genes])
metascore.obj1[[str_c(selected_genes, "_group")]] <-
  if_else(meta[,selected_genes] > median(meta[,selected_genes]),
          "high", "low")
table(metascore.obj1$OXPHOS_UCell_group)
#######高低表达组的差异火山------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#差异分析--
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
### 高低表达组的差异火山图 ----
DefaultAssay(metascore.obj1)
DefaultAssay(metascore.obj1) <- "RNA"
Idents(metascore.obj1) <- "OXPHOS_UCell_group"
markers <- FindMarkers(metascore.obj1,
                       ident.1 ="high",
                       ident.2 = "low",
                       min.pct = 0.1,
                       logfc.threshold = 0.1)
deg2 <- markers
deg2$group <- "OXPHOS_UCell"
deg2$change <- "not"
deg2$change[deg2$p_val <0.05 & deg2$avg_log2FC >0.25] <- "up"
deg2$change[deg2$p_val <0.05 & deg2$avg_log2FC < -0.25] <- "down"
deg2$log10Pvalue <- log10(deg2$p_val)
features <- c(feature1,feature2)
features <- c("LYN","CIITA","CASP4","TXNRD1")
a <- deg2[features,]
library(ggplot2)
library(dplyr)
library(ggrepel)
p1 <- ggplot(deg2, aes(x =avg_log2FC, y=-log10Pvalue, colour=change)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.9, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#7db954","grey","pink")) + xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围 #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = 2.853872, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="black",lwd=0.8)+
  labs(x="log2FC", y="-log10pvalue") +  #x、y轴标签
  ggtitle("OXPHOS-high vs OXPHOS-low") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()) 
P2 <- p1+ geom_label_repel(data = a,
                           aes(x =avg_log2FC, y=-log10Pvalue, label = features),
                           size = 2, fill="#CCFFFF",max.overlaps = 300,label.padding = 0.1)
ggsave(paste0(output,"OXPHOS火山图.pdf"), plot = P2, width = 5.29, height = 3.69)
write.csv(deg2,file = paste0(output,"OXPHOS_UCell差异基因.csv"))
###功能富集------
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
deg2$Gene <- rownames(deg2)
up_deg <- subset(deg2,subset = deg2$change %in% "up")
down_deg <- subset(deg2,subset = deg2$change %in% "down")
up_deg <-
  up_deg[order(up_deg$avg_log2FC, decreasing = T),] %>% pull("Gene")
down_deg <-
  down_deg[order(down_deg$avg_log2FC, decreasing = T),] %>% pull("Gene")
up_bp <-
  enrichGO(
    up_deg,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term1 <- up_bp@result
term1$celltype <- "OXPHOS_high"
down_bp <-
  enrichGO(
    down_deg,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term2 <- down_bp@result
term2$celltype <- "OXPHOS_low"
####功能富集结果----
#提取
subset_up_term <- subset(term1,subset = term1$Description %in% c("cell communication by electrical coupling","cellular carbohydrate metabolic process","negative regulation of interleukin-12 production","respiratory system development","organic hydroxy compound catabolic process","negative regulation of cytokine production","cell redox homeostasis"))
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
subset_up_term$text_x <- rep(0.03,7)
p6 <- ggplot(data = subset_up_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "pink", high = "pink") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p6
ggsave(paste0(output,"GO富集巨噬分组柱状图上调.pdf"), plot = p6, width = 4.45, height = 4.05) 
###下调---
subset_down_term <- subset(term2,subset = term2$Description %in% c("MHC class II protein complex assembly","MHC protein complex assembly","peptide antigen assembly with MHC class II protein complex","cell activation involved in immune response","calcium ion homeostasis","immune response-regulating cell surface receptor signaling pathway","positive regulation of interleukin-2 production"))
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
subset_down_term$text_x <- rep(0.03,7)
p6 <- ggplot(data = subset_down_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#7db954", high = "#7db954") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p6
ggsave(paste0(output,"GO富集巨噬分组柱状图下调.pdf"), plot = p6, width = 4.45, height = 4.05)
#### GSEA ----
markers1 <- FindMarkers(metascore.obj1,
                        ident.1 ="high",
                        ident.2 = "low",
                        min.pct = 0,
                        logfc.threshold = 0)
deg1 <- markers
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(ggplot2)
gsea.input <-
  deg1[order(deg1$avg_log2FC, decreasing = T),]

category <- "C5"
genesets <- msigdbr(species = "Homo sapiens",#####Homo sapiens,Mus musculus
                    category = category)
genesets <- subset(genesets,
                   select = c("gs_name", "gene_symbol"))

genelist <-
  structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
genelist <- na.omit(genelist)
res <- GSEA(genelist, TERM2GENE = genesets, minGSSize = 1,maxGSSize = 5000,pvalueCutoff = 1)
library(enrichplot)
saveRDS(res, file.path(output, "Mitophagy_C5_gsea.rds"))

####可视化-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
a <- res@result
setid2 <- c("GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE","GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE","GOBP_INFLAMMATORY_RESPONSE","GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE","GOBP_ENERGY_HOMEOSTASIS","GOBP_CELLULAR_HOMEOSTASIS")

p7 <- gseaNb(object = res,
             geneSetID = setid2,
             newGsea = T,
             addPval = T,
             rmHt = T,
             pvalX = 0.8,
             pvalY = 0.5,
             pFill = "transparent",
             geneCol = "#009933",
             newCurveCol = c("#7db954","grey","pink"))

p7
ggsave(paste0(output,"OXPHOS_GSEA.pdf"), plot = p7, width = 7.24, height = 5.13)
### 小提琴图-----
#通路打分---
dbs1 <- split(genesets$gene_symbol, genesets$gs_name)
metascore.obj1 <- AddModuleScore_UCell(metascore.obj1, features = dbs1,
                                    ncores = 50,maxRank = 1000000)
library(ggpubr)
my_comparisons <- c("low","high")
meta <- metascore.obj1@meta.data
colnames(meta) <- gsub("_UCell|GOBP_", "", colnames(meta))
colnames(meta) <- tolower(colnames(meta))
colnames(meta) <- gsub(" ","_",colnames(meta))
meta1 <- meta[,c("oxphos_group","cellular_homeostasis")]
p9 <- ggviolin(meta1, x='oxphos_group',y="cellular_homeostasis",
               xlab = 'OXPHOS group', ylab = "Cellular homeostasis",
               fill='oxphos_group',palette = c("#7db954", "pink"),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta1$cellular_homeostasis),# 取得分的中值
             linetype = "dashed", size = 0.5, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("low","high")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"稳态小提琴.pdf"), plot = p9, width = 3.41, height = 5.92) 
meta1 <- meta[,c("oxphos_group","inflammatory_response")]
p9 <- ggviolin(meta1, x='oxphos_group',y="inflammatory_response",
               xlab = 'OXPHOS group', ylab = "Inflammatory response",
               fill='oxphos_group',palette = c("#7db954", "pink"),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta1$inflammatory_response),# 取得分的中值
             linetype = "dashed", size = 0.5, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("low","high")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"OXPHOS炎症反应小提琴.pdf"), plot = p9, width = 3.41, height = 5.92) 


p9 <- ggviolin(meta, x='mitophagy_group',y="gocc_mhc_class_ii_protein_complex",
               xlab = 'Mitophagy group', ylab = "MHC class ii protein complex",
               fill='mitophagy_group',palette = c("#7db954", "pink"),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$gocc_mhc_class_ii_protein_complex),# 取得分的中值
             linetype = "dashed", size = 0.5, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("low","high")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"MHCII_小提琴.pdf"), plot = p9, width = 3.41, height = 5.92) 


#####转录因子可视化------
#####转录因子可视化-----
rm(list = ls())
#### 加载包 ----
library(AUCell)
library(SCENIC)
library(Seurat)
library(pheatmap)
library(patchwork)
metascore.obj <- NormalizeData(metascore.obj)
DefaultAssay(metascore.obj) <- "RNA"
regulonAUC <- importAUCfromText("~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

cellInfo <- metascore.obj@meta.data

# 将auc得分矩阵添加进seurat

metascore.obj[["scenic"]] <-
  CreateAssayObject(counts = t(getAUC(regulonAUC)))
Idents(metascore.obj) <- metascore.obj$group
DefaultAssay(metascore.obj) <- "scenic"

library(readr)
TFS_targets <- read_csv("~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/tfs_targer.csv")



library(tidyverse)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
# db <- subset(dbs,subset = dbs$pathway %in% "Mitophagy")
# tf_seleted <- subset(TFS_targets, target_gene%in% db$gene) %>% pull(tf)
# tf_seleted <- unique(tf_seleted)
# tf_seleted <- paste0(tf_seleted,"(+),")
tf_seleted <- c("NFKB1(+),","E2F1(+),","FLI1(+),","E2F3(+),","E2F2(+),","STAT1(+),","TAF1(+),")

# 导入auc数据
regulonAUC <- importAUCfromText("~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- metacell.obj@meta.data

# 将auc得分矩阵添加进seurat
# 
# cellInfo <- sub_metacell.obj@meta.data
# Idents(sub_metacell.obj) <- sub_metacell.obj$sub_cell_type
Idents(metascore.obj) <- metascore.obj$group
#### 1. 按照细胞类型进行热图可视化 ----
regulonActivity_byCellType <-
  sapply(split(rownames(cellInfo), Idents(metascore.obj)),
         function(cells)
           rowMeans(t(getAUC(regulonAUC))[, cells]))
# regulonActivity_byCellType <- regulonActivity_byCellType[,-c(6,7)]
a <- data.frame(names = row.names(regulonActivity_byCellType), regulonActivity_byCellType)


pheatmap(
  regulonActivity_byCellType[tf_seleted,],
  scale = "row",
  # annotation_col=meta,
  color =  colorRampPalette(c("navy", "white", "firebrick3"))(50),
  # breaks=seq(-1, 1, length.out = 100),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white",
  cellheight = 50,
  cellwidth = 50,
  fontsize = 10,
  angle_col = "45",
  cluster_cols = F,
  cluster_rows = F,
  filename = "~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/deg-heatmap_by_celltype.pdf"
)
#### 3. 转录因子表达水平可视化 ----
DefaultAssay(metascore.obj) <- "RNA"
tf <- "STAT1"
p1 <- FeaturePlot(
  metascore.obj,
  features = tf,
  reduction = "umap",
  pt.size = 2,
  order = T,
  min.cutoff = 0)
p1
p1 <- VlnPlot(sub_metacell.obj,
              features = tf,
              cols = c("lightgrey" , "#DE1F1F",
                       "navy"),
              pt.size = 0,
              group.by = "sub_cell_type") 
p1
#### 4. 转录因子活性在umap或tsne投影可视化 ----
DefaultAssay(metascore.obj) <- "scenic"
tf <- "STAT1(+),"
p2 <- FeaturePlot(
  metascore.obj,
  features = tf,
  reduction = 'umap',
  min.cutoff = 0,
  cols = c("lightgrey" , "#DE1F1F"),
  pt.size = 2,
  order = T
)
p2
library(ggpubr)
p2 <- VlnPlot(sub_metacell.obj,
              features = tf,
              cols = c("lightgrey" , "#DE1F1F",
                       "navy"),
              pt.size = 0,
              group.by = "group",sort = T) 
p2

p <- p1 + p2 + plot_layout(ncol = 2)
p


#######挑选转录因子去cytoscape做网络------
a <- read.csv("~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/tfs_targer.csv")
HCM_diff <- read_excel("~/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
HCM_diff <- HCM_diff[,c(1,3,15,16,17)]
HCM_diff <- subset(HCM_diff,HCM_diff$Gene %in% dbs$gene)
colnames(HCM_diff) <- c("Gene","cell_type","logFC","pvalue","adjusted_pvalue")
HCM_diff <- subset(HCM_diff,subset = HCM_diff$`adjusted_pvalue` < 0.01)
mac_diff <- subset(HCM_diff,subset = HCM_diff$cell_type %in% "Macrophage")
up <- subset(mac_diff,subset = mac_diff$logFC > 0) %>% pull(Gene)
down <- subset(mac_diff,subset = mac_diff$logFC < 0) %>% pull(Gene)
up_TF <- subset(a,subset = a$target_gene %in% up)
down_TF <- subset(a,subset = a$target_gene %in% down)
write.csv(up_TF,file = "~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/mito_updeg_TF.csv")
write.csv(down_TF,file = "~/scRNA-heart-mitochodria/results/Macrophage/pyscenic/HCM-NF/mito_downdeg_TF.csv")
######在HCM中根据线粒体自噬高低表达做cellchat-----
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
metascore.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
Idents(metascore.obj) <- "cell_type"
table(metascore.obj$cell_type)
metascore.obj0 <- subset(metascore.obj,idents = c("Adipocyte",
                                            "Cardiomyocyte","Endocardial_cell","Endothelial_cell","Fibroblast",
                                            "Lymphocyte","Mast_cell","Neuronal_cell","Pericyte","VSMC"))
metascore.obj1 <- subset(metascore.obj,idents = c("Macrophage"))
meta1 <- metascore.obj1@meta.data
selected_genes <- "Mitophagy_UCell"
median(meta1[,selected_genes])
metascore.obj1[[str_c(selected_genes, "_group")]] <-
  if_else(meta1[,selected_genes] > median(meta1[,selected_genes]),
          "high", "low")
table(metascore.obj1$Mitophagy_UCell_group)
group1 <- "low"
group2 <- "high"
# cellchat输出文件路径
metascore.obj2 <- metascore.obj1[, metascore.obj1[["Mitophagy_UCell_group"]] == group1] 
table(metascore.obj2$cell_type)
metascore.obj2$cell_type <- "Mac_Mitophagy_UCell_low"
metascore.obj3 <- metascore.obj1[, metascore.obj1[["Mitophagy_UCell_group"]] == group2] 
table(metascore.obj3$cell_type)
metascore.obj3$cell_type <- "Mac_Mitophagy_UCell_high"
merged <- merge(metascore.obj0, c(metascore.obj2, metascore.obj3))

metascore.obj <- merged

Idents(metascore.obj) <- "cell_type"
#### Preprocessing ----
group <- "HCM"
# group <- "NF"
# cellchat输出文件路径
output.dir <-
  paste0(output, group, "/") # must have "/"
dir.create(output.dir, recursive = T)
metascore.obj <- metascore.obj[, metascore.obj[["group"]] == group] #
names(table(Idents(metascore.obj)))
print(table(Idents(metascore.obj)))

expr <- metascore.obj@assays$RNA@data

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(metascore.obj))
colnames(meta) <- "labels"

unique(meta$labels) # check the cell labels
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <-
  as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-
  CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 20) # do parallel
cellchat <-
  identifyOverExpressedGenes(cellchat) # take a short time
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) # take a short time

#### Compute the communication probability and infer cellular communication network ----
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#### Extract the inferred cellular communication network as a data frame ----
#pass
#### Infer the cell-cell communication at a signaling pathway level ----
cellchat <- computeCommunProbPathway(cellchat)
#### Calculate the aggregated cell-cell communication network ----
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()

pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
  
)
dev.off()

pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"))
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()

pdf(file = paste0(output.dir, "netVisual_chord_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
saveRDS(cellchat, file = paste0(output.dir, "cellchat.rds"))

###细胞通讯具体配受体分析-----
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c(11,10), targets.use = c(1), remove.isolate = FALSE)
#> Comparing communications on a single object
#####multicellchat-----
library(CellChat)

HCM <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure5/311/HCM/cellchat.rds")
HCM@meta <- HCM@meta
NF <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure5/311/NF/cellchat.rds")
group.new = levels(HCM@idents)
NF <- liftCellChat(NF, group.new)



object.list <- list(NF = NF, 
                    HCM = HCM)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)) #因为有报错，原本是FALSE
cellchat

pos.dataset = "HCM"  #实验组，谁和谁比，即分子，case
features.name = pos.dataset
cellchat <-
  identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.1,
    thresh.p = 1
  )

#### 第一部分 气泡、热图、网络展示整体相互作用、通路 ----
#图一，呈现两组之间配受体数量的差异 
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2  

#图二，网络图绘制两组的，粗细-相互作用强度，红色-在上调的信号通路，蓝色-下调 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#图三，热图 
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#图四，网络图 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 5, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#图五，堆叠柱形图 
gg1 <- rankNet(cellchat,
               mode = "comparison", stacked = T,color.use = c("#266E94","#FF9797"), #原本没有comparison = c(1,2,3,4,5,6)，自己加上去的
               do.stat = FALSE) +coord_flip() #原本是TURE
gg1

#图七，气泡图
gg1 <- netVisual_bubble(cellchat,
                        sources.use = c(12,11),            #sources.use = 1,sources.use = c(1,2,3,4)配体的细胞类型选哪个
                        targets.use = c(1,2,3,4,5,6,7,8,9),       #targets.use = c(3:2),targets.use = c(1,2,3,4,5,6)受体的细胞类型选哪个到哪个，targets.use = c(1:4),c(4)是AT1
                        comparison = c(1,2),  
                        signaling = c("CD45","CD46","IGF","MHC-II","MHC-I","CD86","PDGF"),
                        #max.dataset = 2,
                        #title.name = "Increased signaling in LS", 
                        angle.x = 45, remove.isolate = T)
gg1
###图8 2D散点图-----
HCM <- netAnalysis_computeCentrality(
  object = HCM,
  slot.name = "netP",
  net = NULL,
  net.name = NULL,
  thresh = 0.05
)

NF <- netAnalysis_computeCentrality(
  object = NF,
  slot.name = "netP",
  net = NULL,
  net.name = NULL,
  thresh = 0.05
)
object.list <- list(NF = NF)
object.list <- list(HCM = HCM,
                    NF = NF)
gg1 <- netAnalysis_signalingRole_scatter(HCM) + scale_x_continuous(limits = c(0,60)) + scale_y_continuous(limits = c(10,50)) 
gg2 <- netAnalysis_signalingRole_scatter(NF) + scale_x_continuous(limits = c(0,60)) + scale_y_continuous(limits = c(10,50)) 
gg1+gg2 
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg) + scale_x_continuous(limits = c(0,5)) 

#####图9 功能相似性的散点图----
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
####WGCNA分析-----
################################# 加载包 ################################
# rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 
library(Seurat)
setwd(output)
dir.create('WGCNA')
setwd('./WGCNA')
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 
getwd()

################################# 0.输入数据准备 ################################
metascore.obj <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
metascore.obj1 <- subset(metascore.obj,subset = cell_type %in% "Macrophage")
# metascore.obj1 <- NormalizeData(metascore.obj1)
metascore.obj1@meta.data
#按照样本取平均，获得表达矩阵
expr <- AverageExpression(metascore.obj1, group.by = "sample", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
# 更改矩阵行名为样本＋疾病，为表型数据做准备
meta1 <- metascore.obj1@meta.data
selected_genes <- "Mitophagy_UCell"
median(meta1[,selected_genes])
metascore.obj1[[str_c(selected_genes, "_group")]] <-
  if_else(meta1[,selected_genes] > median(meta1[,selected_genes]),
          "high", "low")
table(metascore.obj1$Mitophagy_UCell_group)
expr2 <- as.data.frame(AverageExpression(metascore.obj1, group.by = c("Mitophagy_UCell_group","sample"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(expr2))
name2 <- str_split_fixed(name2$`colnames(expr2)`,"_",2)
name2 <- as.data.frame(name2)
table(name2$V1)
expr2 <- t(expr2)
expr2 <- as.data.frame(expr2)
rownames(expr2) <- paste0(rownames(expr2),"#",name2$V1)
expr2[1:4,1:4]
#                             MIR1302-2HG FAM138A OR4F5 AL627309.1
# high_LV_1422_1_hcm#high           0       0     0  0.0000000
# high_LV_1422_2_hcm#high           0       0     0  0.0000000
# 把命名赋值给expr
name1 <- as.data.frame(colnames(expr))
name1 <- str_split_fixed(name1$`colnames(expr)`,"_",4)
name1 <- as.data.frame(name1)
table(name1$V4)
expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name1$V4)
expr[1:4,1:4]
#                 MIR1302-2HG FAM138A 
# CK158#CHD           0       0   
# CK159#CHD           0       0

### 筛选数据前5000的基因
expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],]
dim(keep_data)
keep_data[1:4,1:4]
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
keep_data[1:4,1:4]
#            MALAT1    PRKG1    
# CK158#CHD  839.0129 64.49691 
# CK159#CHD  932.0754 55.67776

### 创建datTraits，包含分组、表型等信息
df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)

### 给分组加上编号
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)
#           group groupNo
# CK158#CHD   CHD       1
# CK159#CHD   CHD       1
table(datTraits$group)
table(datTraits$groupNo)

#使用输入文件
datExpr0 <- as.data.frame(keep_data)

############################## 1.判断数据质量 ################################
### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  pdf("step1_Sample dendrogram.pdf",width = 8,height = 6)
  p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
            cex.axis = 1, cex.main = 1,cex.lab=1)
  print(p)
  dev.off()
  # ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                  colors = rainbow(length(table(datTraits$group))),
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  p2 <- plotDendroAndColors(sampleTree, sample_colors,
                            groupLabels = "trait",
                            cex.dendroLabels = 0.8,
                            marAll = c(1, 4, 3, 1),
                            cex.rowText = 0.01,
                            main = "Sample dendrogram and trait"
  )
  print(p2)
  dev.off()
  # # Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  # dev.off()
}

##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

### 判断数据质量 : PCA进行分组查看
group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

##保存数据
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")


############################### 2.挑选最佳阈值power ###################################
rm(list = ls())  
load("step1_input.Rdata")
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
power = sft$powerEstimate
power 

# 若无向网络在power小于15或有向网络powerdx小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file = "step2_power_value.Rdata")


##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr), #默认5000
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 50,    ##越大模块越少
    mergeCutHeight = 0.5, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = T,
    verbose = 3
  )
  table(net$colors) 
  # 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
  # 401 1659  326  317  309  293  290  262  163  160  158  146  129  109  107   95   76
  
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}

## 模块可视化，层级聚类树展示各个模块
if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")



####################### 4.关联基因模块与表型 #####################################
rm(list = ls())
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")

## 模块与表型的相关性热图
if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) #get the group
  
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap2.pdf",
      width = 2*length(colnames(design)),
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 # colors = blueWhiteRed(50),
                 colors = colorRampPalette(c("#105eb7", "white","#d7131a"))(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}
#条形图
if(T){
  y = as.data.frame(design[ ,"hcm"]);
  GS1=as.numeric(cor(y,datExpr,use="p"))
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,moduleColors , mean , na.rm=T)
  pdf("step4-Module-single_trait-relationship_barplot_high.pdf",width = 8,height = 5)
  plotModuleSignificance(GeneSignificance,moduleColors)
  dev.off()
}

### 模块与表型的相关性boxplot图
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names")
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width = 8,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) #排布为每行2个
  dev.off()
}


### 基因与模块、表型的相关性散点图
#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数，
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

# 选择离散性状的表型
levels(datTraits$group)
choose_group <- "hcm"

if(T){
  modNames <- substring(names(MEs), 3)
  
  ### 计算模块与基因的相关性矩阵
  ## Module Membership: 模块内基因表达与模块特征值的相关性
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  ###  计算性状与基因的相关性矩阵
  ## Gene significance，GS：比较样本某个基因与对应表型的相关性
  ## 连续型性状
  # trait <- datTraits$groupNo
  ## 非连续型性状，需转为0-1矩阵, 已存于design中
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  ### 可视化基因与模块、表型的相关性.
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己想要的模块
  selectModule <- modNames  ## 全部模块批量作图
  pdf("step4_gene-Module-trait-significance_HCM.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}


#########################  5. WGCNA可视化：TOMplot  Eigengene-adjacency-heatmap ##################################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  ## draw all genes
  if(T){
    geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("step5_TOMplot_Network-heatmap.png",width = 800, height=600)
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(T){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


### 模块相关性展示 Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    # 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    # design
    HCM = as.data.frame(design[,1])
    names(HCM) = "HCM"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, HCM))
  }
  pdf("step5_module_cor_Eigengene-dendrogram_HCM.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="",
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

#################### 6. 选择感兴趣基因模块进行GO分析 ####################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

### 条件设置
OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
# black         blue        brown         cyan        green  greenyellow         grey    lightcyan      magenta midnightblue         pink       purple 
# 262          326          317          107          293          146          401           76          160           95          163          158 
# red       salmon          tan    turquoise       yellow 
# 290          109          129         1659          309 

choose_module <- c("black","blue","brown","green","grey","magenta","pink","purple","red",
                   "turquoise","yellow") 
# pink cyan green tan

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F)
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "BP",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",   #指定多重假设检验矫正的方法,一般选择 "BH" 或 "fdr"，BH较严格，fdr较温和（计算的q小些）
    pvalueCutoff = 0.05
    # qvalueCutoff = 0.25
    # minGSSize &maxGSSize：是富集的最小/大的基因集的大小（基因数目）
  )
  
  ###精简GO富集的结果,去冗余
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res,
    cutoff=0.5,
    by="p.adjust",
    select_fun=min
  )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  ### 绘制dotplot图
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE, #将有overlap的结果也展示出来
                  label_format=90)
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 12,
         height = 15)
}

############################### 7.感兴趣基因模块绘制热图 ######################################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step3_genes_modules.Rdata")
table(moduleColors)

module = "green"
### 感兴趣模块画热图
if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat)) #对基因做scale，并转置表达矩阵为行为基因、列为样本形式
  # n[n>2]=2
  # n[n< -2]= -2
  # n[1:4,1:4]
  
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames =T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col =ac,
                     width = 8,
                     height = 6,
                     angle_col=45,
                     main = paste0("module_",module,"-gene heatmap"),
                     filename = paste0("step7_module_",module,"_Gene-heatmap.pdf"))
  
}

################### 8.感兴趣模块基因导出 VisANT or cytoscape ######################
rm(list = ls())
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")
module = "green"  #turquoise
if(T){
  ### 提取感兴趣模块基因名
  gene <- colnames(datExpr)
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  # write.table(modgene,paste0("step8_",module,"_modgene.csv"))
  
  ### 模块对应的基因关系矩阵
  TOM <- TOMsimilarityFromExpr(datExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  ### 筛选连接度最大的top100基因（核心基因/hub基因）
  nTop = 100
  IMConn = softConnectivity(datExpr[, modgene]) #计算连接度
  top = (rank(-IMConn) <= nTop) #选取连接度最大的top100
  filter_modTOM <- modTOM[top, top]
  
  # for visANT
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  # for cytoscape
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste("step8_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,  #weighted权重筛选阈值，可调整
                                  nodeNames = modgene[top],
                                  nodeAttr = moduleColors[inModule][top])
}


#######WGCNA种提取模块gene进行打分-----
modulegene <- read.csv("~/scRNA-heart-mitochodria/figure/final/Figure5/311/WGCNA/step6_gene_moduleColors.csv")
library(ggplot2)
library(ggpubr)
####模块基因进行功能富集-----
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
######模块基因打分-----
library(Seurat)
library(readr)
library(UCell)
library(ggpubr)
metascore.obj <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
metascore.obj1 <- subset(metascore.obj,subset = cell_type %in% "Macrophage")
dbs1 <- split(modulegene$gene, modulegene$module)
metascore.obj1 <- AddModuleScore_UCell(metascore.obj1, features = dbs1,maxRank = 150000,
                                   ncores = 20)
library(msigdbr)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
dbs2 <- Dataset[grepl("CYTOKINE|MHC", Dataset$gs_name, ignore.case = TRUE), ]
dbs2 <- split(dbs2$gene_symbol, dbs2$gs_name)
metascore.obj1 <- AddModuleScore_UCell(metascore.obj1, features = dbs2,
                                      ncores = 50)
meta <- metascore.obj1@meta.data
colnames(meta)

####相关性热图绘制------
library(readxl)
library(reshape2)
library(UCell)
library(psych)
# adata <- aggregate(meta, by=list(saple=meta$sample),mean)
adata <- meta
colnames(adata)
colnames(adata) <- gsub("_UCell", "", colnames(adata) )
colnames(adata) <- gsub("GOBP_", "", colnames(adata) )
colnames(adata)  <- gsub("_", " ", colnames(adata) )
colnames(adata)  <- tolower(colnames(adata) )
cor_function_meta <- corr.test(adata[,c(165:171)],adata[,c(120,99,216,250,249)], method = "spearman", adjust = "fdr")#做相关性


####可视化-----
cp <- as.data.frame(cor_function_meta$p)
cr <- as.data.frame(cor_function_meta$r)
cr <- round(cr,2)

cp$names <- rownames(cp)
longp <- melt(cp,idvar = "names",v.names = "abd",direction = "long")
colnames(longp) <- c('na','va','p')

cr$names <- rownames(cr)
longr <- melt(cr,idvar = "names",v.names = "abd",direction = "long")

allnew <- as.data.frame(cbind(longr,longp$p))
colnames(allnew) <- c('names','variable','r','p')
#allnew$variable <- c(rep('Agriculture',9),rep('Forest',9),rep('Wetland',9),rep('Grass',9),rep('Desert',9))

allnew[which(allnew$p<0.001),'sig'] <- '***'
allnew[which(allnew$p<0.01 & allnew$p>0.001),'sig'] <- '**'
allnew[which(allnew$p<0.05 & allnew$p>0.01),'sig'] <- '*'
allnew[,1] <- factor(allnew[,1], levels = c("grey","brown","blue","turquoise", "red", "yellow","green"))
# 按照第一列的值进行排序
sorted_allnew <- allnew[order(allnew[,1]), ]

p13 <- ggplot(sorted_allnew, aes(variable,names)) +
  geom_tile(aes(fill=r),color="white",size=5) +
  geom_text(aes(label=r), color="black", size=4) + # 把相关性添加进去
  geom_text(aes(label=sig), color="black", size=4,vjust = 1.8)+ # 把星号添加进去
  scale_fill_gradient2(low='#2166AC', high='#B2182B',mid = 'white', limit=c(-1,1),  
                       name=paste0("*p < 0.05","\n\n","**p < 0.01","\n\n","***p < 0.001","\n\n","Correlation")) + # 把P值添加到图例，这一步非常巧妙
  labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
  theme(axis.text.x = element_text(size=10,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())
p13
ggsave(paste0(output,"WGCNA模块相关性热图.pdf"), plot = p13, width = 7.92, height = 4.72)


###散点图----
#第1种画法
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
P <- adata %>%
  ggplot(aes(x = mitophagy, y = green)) +
  geom_point(size = 0.5,color="#40a8c4",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#7db954", fill = "#9ACD32")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = mitophagy, y = green),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 15,  color = "black"),
        axis.title.y = element_text(size = 15,  color = "black")) +
  labs(y = "Green",x = "Mitophagy")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#40a8c4"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"WGCNA_green与线粒体自噬散点图.pdf"), plot = p4, width = 3.83, height = 3.35)

P <- adata %>%
  ggplot(aes(x = oxphos, y = green)) +
  geom_point(size = 0.5,color="#40a8c4",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#7db954", fill = "#9ACD32")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = oxphos, y = green),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 15,  color = "black"),
        axis.title.y = element_text(size = 15,  color = "black")) +
  labs(y = "Green",x = "OXPHOS")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#40a8c4"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"WGCNA_green与OXPHOS散点图.pdf"), plot = p4, width = 3.83, height = 3.35)
colnames(adata)
P <- adata %>%
  ggplot(aes(x = `gomf mhc class ii protein binding`, y = green)) +
  geom_point(size = 0.5,color="#40a8c4",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#7db954", fill = "#9ACD32")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = `gomf mhc class ii protein binding`, y = green),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 15,  color = "black"),
        axis.title.y = element_text(size = 15,  color = "black")) +
  labs(y = "Green",x = "MHC class ii protein binding")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#40a8c4"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"WGCNA_green与mhciiproteinbinding散点图.pdf"), plot = p4, width = 3.83, height = 3.35)

P <- adata %>%
  ggplot(aes(x = `positive regulation of macrophage cytokine production`, y = green)) +
  geom_point(size = 0.5,color="#40a8c4",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#7db954", fill = "#9ACD32")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = `positive regulation of macrophage cytokine production`, y = green),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 15,  color = "black"),
        axis.title.y = element_text(size = 15,  color = "black")) +
  labs(y = "Green",x = "Positive regulation of macrophage cytokine production")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#40a8c4"),
                 yparams = list(fill = "#7db954"))
p4
ggsave(paste0(output,"WGCNA_green与macrophagecytokine散点图.pdf"), plot = p4, width = 3.83, height = 3.35)





meta <- meta[,c(18,180,181,182)]
table(meta$Mitophagy_UCell_group)
###设置组间对比，排列组合
my_comparisons <- list( c("high", "low"), 
                        c("low", "zero"), 
                        c("zero", "high"))

#set color
green <- "green"
blue   <- "#5bc0eb"
grey   <- "#8693ab"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"

ggplot(data = meta,aes(x = Mitophagy_UCell_group, #分组列名
                       y = brown_UCell, #连续变量列名
                       fill = Mitophagy_UCell_group))+ #按分组填充颜色
  scale_fill_manual(values = c(green,darkred,blue, grey, lightred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  # geom_point(shape = 21, size=0, # 点的性状和大小
  #            position = position_jitterdodge(), # 让点散开
  #            color="black", alpha = 1) +
  theme_classic() + 
  ylab("brown_UCell") +
  xlab("Mitophagy_UCell_group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(method = "kruskal.test", label.y = min(meta$brown_UCell))

######基因表达的分裂小提琴图------
########小提琴图
dbs <- read_csv("/home/duansq/projects/20221014_oocytes-mitochondria/data/datasets/human_mito_v3.0.csv")
selected_pathway <- c("Complex I",
                      "Complex II",
                      "Complex III",
                      "Complex IV",
                      "Complex V",
                      "TCA cycle",
                      "OXPHOS",
                      "ROS and glutathione metabolism",
                      "Translation")
walk(pathways, function(x) {
  # selected_pathway <- "Complex I"
  selected_pathway <- x
  features <- dbs |>
    filter(pathway == selected_pathway) |>
    pull(gene)
  
  # selected_pathway <- "Ribosomal genes"
  # features <- grep("^RP[SL]", rownames(human_ribo_TPM), value = T)
  #
  matrix <- expr[features,c("Y","O") ]
  matrix <- na.omit(matrix)
  seurat_obj$group <- factor(seurat_obj$aging,
                             levels = c("Y", "O"))
  # Idents(seurat_obj) <- seurat_obj$aging
  # seurat_obj <- ScaleData(seurat_obj,split.by = "aging")
  p <- VlnPlot(
    seurat_obj,
    features = features,
    group.by = "group",
    split.by = "group",
    flip = T,
    stack = T,
    # slot = "scale.data",
  ) +
    theme_cat(aspect.ratio = 0.08) +
    theme(
      strip.text.y = element_text(hjust = 0, face = "italic"),
      axis.title.x = element_blank(),
      legend.position = "top",
      legend.margin = margin(b = -8)
    )
  p
  ggsave(
    file.path(outdir, str_c(
      str_replace_all(selected_pathway, " ", "_"), "_violin.pdf"
    )),
    plot = p,
    height = 10,
    width = 10
  )
})

library(ggpubr)
# metascore.obj <- readRDS("~/scRNA-heart-mitochodria/results/Macrophage/线粒体分解/高低0三组/seurat_obj.rds")
metascore.obj <- NormalizeData(metascore.obj)
expr <- as.matrix(metascore.obj@assays$RNA@data)
expr1 <- expr[c("CD14","DHX9","EPHB2","IL6R","NFKB1","TGFB1","PTPRC"),]
expr1 <- t(expr1)
meta <- metascore.obj@meta.data
meta1 <- meta[,c("disease","group")]
adata <- cbind(meta1,expr1)
library(reshape2)
library(ggunchained)
library(plyr)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
a_long<-melt(adata,
             id.vars = c("disease","group"),#需要保留不参与聚合的变量,
             variable.name='gene',
             value.name='expression')
Data_summary <- summarySE(a_long, measurevar="expression", groupvars=c("disease","gene"))
head(Data_summary)

p<-ggplot(a_long,aes(x=gene,y=expression,fill=disease))+geom_split_violin()+
  theme_bw()+
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.25),size= 1)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.05, 
                position= position_dodge(0.25), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c('#00BFC4','#eb990c'))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank())+
  stat_compare_means(aes(group = disease),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(a_long$expression),
                     hide.ns = T)+
  ylab('Pathway Score')
p
#######基因相关性计算-----
expr <-
  AverageExpression(metascore.obj,
                    group.by = "biosample_id", assays = "RNA")[["RNA"]]
expr1 <- expr[c("CD14","DHX9","EPHB2","IL6R","NFKB1","TGFB1","PTPRC"),]
expr1 <- t(expr1)
M <- cor(expr1)
corrplot(M,
         add = F,  #增加一个新图
         type = 'lower', #在左下角
         method = 'number', # 类型为数字
         tl.pos = 'n', #不添加文字标签
         cl.pos = 'n')
corrplot(M, type = 'lower', method = 'number', order = 'hclust',
         tl.col = 'black', cl.ratio = 0.2, tl.srt = 45, 
         col = rev(COL2('PuOr', 10)))

####打分相关的分析-----
library(reshape2)
library(UCell)
library(Seurat)
library(readr)
library(corrplot)
metascore.obj <- readRDS("~/scRNA-heart-mitochodria/data/HCM-NF.rds")
metascore.obj <- subset(metascore.obj,subset = cell_type %in% "Macrophage")
metascore.obj <- NormalizeData(metascore.obj)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
turqu <- as.data.frame(turquoise)
turqu$pathway <- "turquoise"
colnames(turqu) <- c("gene","pathway")
turqu <- turqu[,c(2,1)]
dbs <- rbind(dbs,turqu)
dbs1 <- split(dbs$gene, dbs$pathway)
metascore.obj <- AddModuleScore_UCell(metascore.obj, features = dbs1,maxRank = 150000,
                                   ncores = 20)
meta <- metascore.obj@meta.data
meta1 <- meta[,c("TCA.cycle_UCell","Calcium.cycle_UCell","turquoise_UCell")]
M <- cor(meta1)
corrplot(M,
         add = F,  #增加一个新图
         type = 'lower', #在左下角
         method = 'circle', # 类型为数字
         tl.pos = 'n', #不添加文字标签
         cl.pos = 'n')
corrplot(M,
         add = FALSE,
         type = 'lower',
         method = 'circle',
         col = colorRampPalette(c("green", "red"))(20))
corrplot(M,
         add = FALSE,
         type = 'lower',
         method = 'ellipse',
         col = colorRampPalette(c("blue", "white", "firebrick"))(20))
p4
print(p)

#######WGCNA种提取模块gene功能富集-----
modulegene <- read.csv("~/scRNA-heart-mitochodria/figure/final/Figure5/311/WGCNA/step6_gene_moduleColors.csv")
######clusterprofile转录因子靶基因功能富集网络图-------
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
turquoise <- subset(modulegene,subset = module %in% "turquoise") %>% pull(gene)
green <- subset(modulegene,subset = module %in% "green") %>% pull(gene)
blue <- subset(modulegene,subset = module %in% "blue") %>% pull(gene)
grey <- subset(modulegene,subset = module %in% "grey") %>% pull(gene)
red <- subset(modulegene,subset = module %in% "red") %>% pull(gene)
brown <- subset(modulegene,subset = module %in% "brown") %>% pull(gene)
yellow <- subset(modulegene,subset = module %in% "yellow") %>% pull(gene)


data <- list(turquoise,green,blue,grey,red,brown,yellow)
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
gene_cluster <- data
names(gene_cluster)=c("turquoise","green","blue","grey","red","brown","yellow")

xx <- compareCluster(gene_cluster, 
                     fun='enrichGO',
                     ont= 'BP',
                     OrgDb='org.Hs.eg.db' ,
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1
)


xx <- pairwise_termsim(xx)
a <- xx@compareClusterResult
a <- subset(a,subset = Description %in% c("phospholipid biosynthetic process","activation of innate immune response","immune response-activating signaling pathway","glycerophospholipid metabolic process","activation of immune response","regulation of cell-matrix adhesion","integrin-mediated signaling pathway","regulation of apoptotic signaling pathway","positive regulation of cytokine production","cellular response to calcium ion","fatty acid derivative metabolic process","fatty acid derivative biosynthetic process","fatty acid metabolic process","cellular response to ATP"))
xx@compareClusterResult <- a  
pdf("~/scRNA-heart-mitochodria/figure/final/Figure5/311/go_WGCNA2.pdf",height = 8,width = 8)
xx <- pairwise_termsim(xx)
p23 <- emapplot(
  xx,
  cex_category = 2,
  cex_line = 0.1,
  cex_label_category = 0.89,
  layout = "kk",
  # graphopt, fr
  # shadowtext = F,
  repel = T,
  nWords = 1,
  legend_n = 4,
  node_label = "category",
  showCategory = 5
) +
  scale_fill_manual(values = c('turquoise',
                               '#3cb371',
                               '#3498db',
                               "grey",
                               '#e63946','#ab7967',"#f1c40f"))
dev.off()
p23
ggsave(paste0(output,"WGCNA功能富集网络图.pdf"), plot = p23, width = 6.68, height = 5.20)
