setwd("~/scRNA-heart-mitochodria")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure3/311/"
dir.create(output, recursive = TRUE)
######差异基因的火山图----
library(readr)
library(Seurat)
seurat.obj <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
metacell.obj <- qs::qread("~/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
library(readxl)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg1 <- subset(deg,subset = deg$`Cell Type` %in% "Cardiomyocyte")
deg1$change <- ""
deg1$change[deg1$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg1$change[deg1$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg1$change[deg1$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg1 <- deg1[,c(1,15,16,17,21)]
colnames(deg1) <- c("Gene","logFC","Pvalue","Adjusted_Pvalue","change")
deg1$log10Pvalue <- log10(deg1$Pvalue)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
a <- subset(deg1,subset = deg1$Gene %in% dbs$gene)
features <- unique(a$Gene)
library(ggplot2)
library(dplyr)
library(ggrepel)
p1 <- ggplot(deg1, aes(x =logFC, y=-log10Pvalue, colour=change)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5","grey","#ff4757")) + xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围 #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = 2.853872, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="logFC", y="-log10pvalue") +  #x、y轴标签
  ggtitle("HCM/NF(Cardiomyocyte)") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()) 
p2 <- p1+ geom_label_repel(data = a,
                           aes(x =logFC, y=-log10Pvalue, label = features),
                           size = 1.5, fill="#CCFFFF",max.overlaps = 300,label.padding = 0.1)
p2
ggsave(paste0(output,"火山图.pdf"), plot = p2, width = 5.95, height = 4.22)
#####功能富集-----
up <- subset(deg1,deg1$change %in% "up")
down <- subset(deg1,deg1$change %in% "down")
down_deg <- subset(down,subset = down$Gene %in% dbs$gene)
up_deg <- subset(up,subset = up$Gene %in% dbs$gene)
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
down_deg <-
  down_deg[order(down_deg$logFC, decreasing = T),] %>% pull("Gene")

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
down_term$celltype <- "Cardiomyocyte"

write.table(
  down_term,
  file = "~/scRNA-heart-mitochodria/results/Cardiomyocyte/function/HCM-NF/down.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(down_bp, "~/scRNA-heart-mitochodria/results/Cardiomyocyte/function/HCM-NF/down.rds")
#up-----
up_deg <-
  up_deg[order(up_deg$logFC, decreasing = T),] %>% pull("Gene")

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
up_term <- up_bp@result
up_term$celltype <- "Cardiomyocyte"


write.table(
  up_term,
  file = "~/scRNA-heart-mitochodria/results/Cardiomyocyte/function/HCM-NF/up.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(up_bp, "~/scRNA-heart-mitochodria/results/Cardiomyocyte/function/HCM-NF/up.rds")

######功能富集柱状图-----
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
subset_down_term <- subset(down_term,subset = down_term$Description %in% c("cellular respiration","energy derivation by oxidation of organic compounds","respiratory electron transport chain","fatty acid beta-oxidation","fatty acid catabolic process","oxidative phosphorylation","ATP synthesis coupled electron transport","electron transport chain","aerobic respiration","tricarboxylic acid cycle"))
subset_down_term$text_x <- rep(0.03,10)
p3 <- ggplot(data = subset_down_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#EFEFEF", high = "#546de5") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p3
ggsave(paste0(output,"下调柱状图.pdf"), plot = p3, width = 4.63, height = 3.35)
subset_up_term <- subset(up_term,subset = up_term$Description %in% c("amino acid metabolic process","nucleotide biosynthetic process","nucleoside phosphate biosynthetic process","branched-chain amino acid metabolic process","purine nucleotide biosynthetic process","ribonucleotide metabolic process","purine-containing compound biosynthetic process","ribose phosphate metabolic process","cellular modified amino acid metabolic process"))
subset_up_term$text_x <- rep(0.03,9)
p4 <- ggplot(data = subset_up_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#F6edef", high = "#ff4757") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p4
ggsave(paste0(output,"上调柱状图.pdf"), plot = p4, width = 4.45, height = 3.76)
#######亚群降维图------
######细胞降维图
seurat.obj@meta.data
sub_seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")

library(RColorBrewer)
cell_type_cols <- brewer.pal(8, "Set2") 
cell_type_cols <- cell_type_cols[3:5]
sub_seurat.obj@meta.data$sub_cell_type <- gsub("Cardiomyocyte_I", "Cardiomyocyte I", sub_seurat.obj@meta.data$sub_cell_type)
sub_seurat.obj@meta.data$sub_cell_type <- gsub("Cardiomyocyte_II", "Cardiomyocyte II", sub_seurat.obj@meta.data$sub_cell_type)
sub_seurat.obj@meta.data$sub_cell_type <- gsub("Cardiomyocyte_III", "Cardiomyocyte III", sub_seurat.obj@meta.data$sub_cell_type)
Idents(sub_seurat.obj) <- "sub_cell_type"
p1  <- DimPlot(sub_seurat.obj, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Cardiomyocyte") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1
ggsave(paste0(output,"亚群umap.pdf"), plot = p1, width = 5.98, height = 4.65)
#######细胞比例的桑吉图------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
metacell.obj$cell <- rownames(metacell.obj@meta.data)
metacell.obj1 <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
seurat.obj1 <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")
Ratio <- metacell.obj1@meta.data %>%
  group_by(group, sub_cell_type) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c('#efb306',
            '#7db954',
            '#852f88',
            '#4e54ac',
            '#0f8096',
            'pink',
            'green')

p2 <- ggplot(Ratio, aes(x =group, y= relative_freq, fill = sub_cell_type,
                        stratum=sub_cell_type, alluvium=sub_cell_type)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Disease',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p2
ggsave(paste0(output,"亚群比例桑吉图.pdf"), plot = p2, width = 5.76, height = 3.62)
#####细胞比例变化饼图-----
library("SeuratObject")
table(seurat.obj1$sub_cell_type)

# 计算每个状态下每种细胞类型的数量
counts <- table(seurat.obj1$disease,seurat.obj1$sub_cell_type)
counts <- counts[c(2,3),c(3,4,5)]
# 计算每个状态下每种细胞类型的比例
prop <- prop.table(counts, margin = 2) * 100  # 转换为百分比

# 绘制饼图
par(mfrow=c(1,3))  # 将绘图区域划分为一行三列

# 遍历每个状态
for (i in 1:3) {
  # 提取当前状态下的细胞类型比例
  prop_sub <- prop[, i]
  # 绘制饼图
  pie(prop_sub, main = paste("Cardiomyocyte", i), col=mycolor)
}
######UCell心肌细胞亚型通路差异柱状图------
########UCELL打分----柱状图可视化-----
library(UCell)
library(limma)
data.metacell <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
data.metacell1 <- subset(data.metacell,subset = cell_type %in% "Cardiomyocyte")
# data.metacell <- NormalizeData(data.metacell)
meta<- data.metacell1@meta.data
es.matrix <- meta[,c(colnames(meta)[grep("_UCell", colnames(meta))])]
meta <- meta[, c("group",
                 colnames(meta)[grep("_UCell", colnames(meta))])]   #提取出含GOBP的
meta <-
  aggregate(meta[,-1], list(meta$group), FUN = mean)
rownames(meta) <- meta[, 1]
meta <- meta[, -1]
colnames(meta) <- gsub("GOBP_", "", colnames(meta))
colnames(meta) <- gsub("_", " ", colnames(meta))
meta <- t(meta)



#### 可视化 ----
library(ComplexHeatmap)
meta <- t(scale(t(meta), center = T, scale = T))
##将行名转换为第一列
meta <- data.frame(names = row.names(meta), meta)
meta <- meta[, -1]
heat <- Heatmap(
  meta,
  cluster_rows = T,
  cluster_columns = T,
  col = colorRampPalette(c("#436eee", "white", "#EE0000"))(10),
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
####线粒体通路与疾病通路的相关性-----
data.metacell1 <- NormalizeData(data.metacell1)
metacell.obj <- qs::qread("~/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
data.metacell1 <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
#### 气泡图展示关注基因与通路的相关性 ----
data <- data.metacell1@assays$RNA@data
data <- as.data.frame(data)
features <- c("IL12B","IL17B","IL2","IL20","IL23A","IL25","IL4","LTA","MDK","LIAS","KLRG1",  ###inflammatory
              "CLASP1","CLASP2","HAS3","NTN4","SMPD3","THSD4",    ####ECM
              "ACE2","ATP1A2","ATP2B1","ATP2B2","CAMK2D","CTNNA3","CXADR","ISL1","SCN1B",  ####CARDIAC_CONTRACTION
              "CRK","EID2","FNTA","GDF10","GIPC1","HDAC2","JUN","SMAD5"     #TGFb
)
data1 <- data[features,]
data1 <- na.omit(data1)
data1 <- as.data.frame(t(data1))
##打分----
library(UCell)
library(psych)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
data.metacell1<- AddModuleScore_UCell(data.metacell1, features = dbs1,
                                  ncores = 20)


adata <- data.metacell1@meta.data
adata1 <- cbind(adata,data1)
colnames(adata1)
###按照样本取均值-----
adata2 <- aggregate(adata1, by=list(sample=adata1$sample),mean)
colnames(adata2)
cor_function_meta <- corr.test(adata2[,c(87,121,145,64,48:52)],adata2[,c(166:199)], method = "spearman", adjust = "fdr")
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
filtered_adata$feature_x <- gsub("_UCell", "", filtered_adata$feature_x )
filtered_adata$feature_x <- gsub("GOBP_", "", filtered_adata$feature_x )
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
p <- filtered_adata %>%
  filter(
         feature_y %in% c("KLRG1","CLASP1","CLASP2","ACE2","GIPC1","JUN","LIAS")) %>%
  # dplyr:: mutate(
  #   feature_x = fct_relevel(feature_x, selected_genes),
  #   feature_y = fct_relevel(feature_y, selected_features)
  # ) %>%
  catdotplot(
    x = feature_y,
    y = feature_x,
    size = -log10(p),
    color = r,
    dot_scale = 7,
    title = "Cardiomyocyte"
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




####心肌细胞三/心肌细胞1和2的差异分析 ----
es.matrix <- t(es.matrix)
case <- "Cardiomyocyte_III"
control1 <- "Cardiomyocyte_I"
control2 <- "Cardiomyocyte_II"
meta <-
  data.metacell@meta.data[, c("sub_cell_type",
                            "group")]
meta$cell_id <- rownames(meta)
cell_id.1 <-
  meta[meta$sub_cell_type == case, ]$cell_id
cell_id.2 <-
  meta[meta$sub_cell_type == control1, ]$cell_id
cell_id.3 <-
  meta[meta$sub_cell_type == control2, ]$cell_id
es.matrix.1 <-
  as.data.frame(es.matrix[, cell_id.1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, c(cell_id.2,cell_id.3)],
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
es.matrix.f %>% mutate(across(where(is.character), as.numeric))  -> es.matrix.f
fiT <- lmFit(es.matrix.f, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
diff <- topTable(fiT3, coef = 1, number = 200000)
t_results <-
  as.data.frame(diff$t, row.names = rownames(diff))
head(t_results)
colnames(t_results) <- c("t")
# 保存t值结果
t_results <- na.omit(t_results)

# colnames(t_results) <- c("t", "Gobp")

# t_results$Gobp = with(t_results, reorder(Gobp, t))
t_results$fill <- ""
t_results[t_results$t >= 2.58,]$fill <-
  "up"
t_results[t_results$t <= -2.58,]$fill <-
  "down"
t_results[abs(t_results$t) < 2.58,]$fill <-
  "no"
t_results$color <- ""
t_results[abs(t_results$t) < 2.58,]$color <-
  "n"
t_results[abs(t_results$t) >= 2.58,]$color <-
  "y"
t_results$feature <- rownames(t_results)
library(ggplot2)

# t_results1 <- subset(t_results,subset = t_results$feature %in% c(rownames(t_results)[grep("INFLAMMA", rownames(t_results))],
#                                                                  rownames(t_results)[grep("CYTOKINES", rownames(t_results))],rownames(t_results)[grep("FIBRO", rownames(t_results))]))

t_results1 <- subset(t_results,subset = t_results$feature %in% c("Creatine.metabolism_UCell","Calcium.uniporter_UCell","Fatty.acid.oxidation_UCell","Lipid.metabolism_UCell","TCA.cycle_UCell","Calcium.cycle_UCell","ROS.and.glutathione.metabolism_UCell","OXPHOS_UCell"))
p7 <-
  ggplot(t_results1, aes(x = reorder(feature,t), y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      "down" = "#779fd3",
      "up" = "#a13037",
      "no" = "#cccccc"
    ),
    guide = "none"
  ) +
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dashed",
    size = 0.5
  ) +
  coord_flip() +
  xlab("Cardiomyocyte III/(I,II)") +
  ylab("t value of score") +
  geom_text(
    data = subset(t_results1, t < 0),
    aes(
      x = feature,
      y = 0.1,
      label = paste0(" ", feature),
      color = color
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
      color = color
    ),
    size = 2.5,
    hjust = "inward"
  ) +
  scale_colour_manual(values = c("y" = "black", "n" = "#cccccc"),
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
p7
ggsave(paste0(output,"通路差异图.pdf"), plot = p7, width = 5.97, height = 3.29)
#####monocle2轨迹-------
library(monocle)
#####monocle2-----
outdir <- "~/scRNA-heart-mitochodria/figure/final/轨迹/Cardiomyocyte/1e-38/"
# qvalues <- 1e-10
# cds <-
#   cat_monocle2(colonocyte,
#                group_by = "sub_cell_type",
#                qvalue = qvalue,
#                outdir = outdir)
cds1 <- qs::qread("~/scRNA-heart-mitochodria/figure/final/轨迹/Cardiomyocyte/1e-38/cds.qs")
# cds1 <- orderCells(cds1,root_state = 4)
p3 <-
  plot_cell_trajectory(cds1, color_by = "sub_cell_type", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p4 <- 
  plot_cell_trajectory(cds1, color_by = "Pseudotime", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p3
ggsave(paste0(output,"轨迹1.pdf"), plot = p3, width = 4.86, height = 3.85)
p4
ggsave(paste0(output,"轨迹2.pdf"), plot = p4, width = 4.86, height = 3.85)
###轨迹山峦图-----
#####轨迹山峦图------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=Biobase::pData(cds1)
plotdf <- plotdf[,c(17,18,19)]
a <- data.metacell1@meta.data
plotdf <- cbind(a,plotdf)
p5 <- ggplot(plotdf, aes(x=Pseudotime,y=group,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(3,31,47),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(c("#779fd3","white","firebrick"))(100))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
p5
ggsave(paste0(output,"轨迹山峦1.pdf"), plot = p5, width = 7.41, height = 5.54)

p6 <- ggplot(plotdf, aes(x=Pseudotime,y=State,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(1.7,31.3,46.5),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(c("#779fd3","white","firebrick"))(100))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
p6
ggsave(paste0(output,"轨迹山峦2.pdf"), plot = p6, width = 7.41, height = 5.54)
# ###线粒体通路打分阴阳表达组的山峦图-----
# data.metacell1@meta.data <- plotdf
# x <- "OXPHOS_UCell"
# median(plotdf[,"OXPHOS_UCell"])
# data.metacell1[[str_c(x, "_group")]] <-
#   if_else(a[,x] > median(plotdf[,"OXPHOS_UCell"]),
#           "high", "low")
# plotdf <- data.metacell1@meta.data
# 
# p6 <- ggplot(plotdf, aes(x=Pseudotime,y=OXPHOS_UCell_group,fill = stat(x))) +
#   geom_density_ridges_gradient(scale=1) +
#   geom_vline(xintercept = c(3.5,13.3),linetype=2)+
#   scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
#   scale_y_discrete("")+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank()
#   ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
# p6
# ggsave(paste0(output,"轨迹山峦2.pdf"), plot = p6, width = 7.41, height = 5.54)
######轨迹热图------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
library(dplyr)
library(tidytree)
library(viridis)
library(scales)
Time_diff <- read.csv("~/scRNA-heart-mitochodria/figure/final/轨迹/Cardiomyocyte/1e-38/diff_test_res.csv",row.names = 1)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
mito_gene <- unique(dbs$gene)
mito_gene1 <- subset(mito_gene,subset = mito_gene %in% rownames(cds1))
genes <- c("ACADSB","ATP5F1E","ATP5IF1","COXAL1","MT-ATP6","MT-CO1","MT-ND1","NDUFA1","UQCRB","CXCL12","IL33","CXCL1","IL7","C7","CFI","CFD","DCN","TGLN","ACE2","BMP10","TNNL1","COMP","INHBA","SULF1")
genes <- subset(genes,subset = genes %in% rownames(cds1))
p <- plot_pseudotime_heatmap(cds1[genes,],
                             num_cluster = 2,
                             show_rownames = T,
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
ggsave(paste0(output,"轨迹热图.pdf"), plot = p, width = 5.93, height = 5.17)
#####轨迹分支热图-----
#### 单细胞轨迹的“分支”分析 ----
BEAM_res <- BEAM(cds1[Time_diff$gene_short_name,], branch_point = 1, cores = 20)
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
pdf(file=file.path(output, "branched_diff_heatmap.pdf"), width=5, height=5)

cdsheat <- cds1[row.names(subset(BEAM_res,
                                qval < 0.05)),]

library(ClusterGVis)
p1 <- plot_genes_branched_heatmap2(cdsheat,
                                   branch_point = 1,
                                   num_clusters = 3,
                                   cores = 10,
                                   branch_colors = c("#979797", "#F05662", "#7990C8"),
                                   return_heatmap = T,
                                   hmcols = colorRampPalette(c("navy","white","firebrick3","firebrick3"))(100),
                                   use_gene_short_name = F,
                                   show_rownames = F)
p1

visCluster(object = p1,plot.type = "heatmap")

library(RColorBrewer)

visCluster(object = p1,
           plot.type = "heatmap",
           ht.col.list = list(col_range = seq(-2,2,length.out = 100),
                              col_color = colorRampPalette(brewer.pal(9,"PRGn"))(100)))

pdf(file=file.path(output, "multiple-branch.pdf"), height = 6,width = 8)
visCluster(object = p1,plot.type = "both")
dev.off()

# 修改热图
# pdf(file = "two-branch.pdf",height = 6,width = 7)
# visCluster(object = p1,plot.type = "both")
# dev.off()

# 添加上富集分析的内容
# enrich for clusters

library(org.Hs.eg.db)
enrich <- enrichCluster(object = p1,
                        OrgDb =org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 100,
                        seed = 5201314)

# check
head(enrich[1:3,])
enrich <- subset(enrich,subset = Description %in% c("extracellular matrix organization","cell adhesion mediated by integrin","calcium ion transport","cell-matrix adhesion","regulation of chemotaxis","regulation of mitochondrion organization","oxidative phosphorylation","phospholipid biosynthetic process","cellular respiration","ATP synthesis coupled electron transport","aerobic respiration","energy derivation by oxidation of organic compounds","histone modification","ATP metabolic process"))
enrich <- enrich[c(1:10,16:20),]
markGenes = c("TIMM17A","ACAA2","ATP5MG","ACOT13","ACAA1","ATP5F1B","ACADSB","ATP5F1E","ATP5IF1","COXAL1","MT-ATP6","MT-CO1","MT-ND1","NDUFA1","UQCRB","CXCL12","IL33","CXCL1","IL7","C7","CFI","CFD","DCN","TGLN","ACE2","BMP10","TNNL1","COMP","INHBA","SULF1")
# PLOT
pdf(file=file.path(output, "branch-enrich-2024043.pdf"),height = 9,width = 16,onefile = F)
visCluster(object = p1,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           go.col = rep(jjAnno::useMyCol("calm",n = 3),each = 5),
           add.bar = T,
           line.side = "left")


dev.off()
#' # 选择要在热图上添加的基因
#' selected_genes <- c("Gene1", "Gene2", "Gene3")
#' 
#' # 添加基因标签
#' heatmap <- add_heatmap_annotation(heatmap, 
#'                                   labels_col = selected_genes, 
#'                                   labels_col_side = "right",
#'                                   col = list(labels_col_side_colors = "black"))
#' 
#' #### 分组 ----
#' 
#' #'Normal'
#' case <- 'NF'
#' cds2 <- cds1[, cds1[["disease"]] == case]
#' # Time_diff <- differentialGeneTest(cds[ordergene,],
#' #                                   cores = 1,
#' #                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
#' # write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)
#' pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
#' Time_diff <- differentialGeneTest(cds2[ordergene,],
#'                                   cores = 10,
#'                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
#' p <- plot_pseudotime_heatmap(
#'   cds2[row.names(subset(Time_diff, qval < 0.01)),],
#'   cluster_rows = T,
#'   num_clusters = 3, cores = 10,
#'   use_gene_short_name = T,
#'   show_rownames = T,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
#' dev.off()
#' gene <- read_excel("~/scRNA-heart-mitochodria/results/Cardiomyocyte/monocle2/轨迹gene.xlsx") %>% pull(gene)
#' gene <- subset(gene,subset = gene %in% rownames(cds3))
#' #'HCM'
#' case <- 'NF'
#' cds2 <- cds1[, cds1[["disease"]] == case]
#' # Time_diff <- differentialGeneTest(cds[ordergene,],
#' #                                   cores = 1,
#' #                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
#' # write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)
#' pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
#' Time_diff <- differentialGeneTest(cds2[ordergene,],
#'                                   cores = 10,
#'                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
#' p <- plot_pseudotime_heatmap(
#'   cds3[gene,],
#'   cluster_rows = F,
#'   num_clusters = 3, cores = 10,
#'   use_gene_short_name = T,
#'   show_rownames = F,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
#' dev.off()
#' 
#' case <- 'HCM'
#' cds3 <- cds1[, cds1[["disease"]] == case]
#' # Time_diff <- differentialGeneTest(cds3[ordergene,],
#' #                                   cores = 10,
#' #                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
#' # write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)
#' pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
#' p <- plot_pseudotime_heatmap(
#'   cds3[row.names(subset(Time_diff, qval < 0.01)),],
#'   cluster_rows = T,
#'   num_clusters = 5, cores = 10,
#'   use_gene_short_name = T,
#'   show_rownames = T,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
#' p
#' dev.off()
#' 
#' 
#' p <- plot_pseudotime_heatmap(
#'   cds2[row.names(subset(Time_diff, qval < 0.01)),],
#'   num_clusters = 5,
#'   show_rownames = T,
#'   return_heatmap = T,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100)
#' )
#' p
#' #ggsave(file.path("~/Retina_cell/results/monocle2/figure1.png"),plot = p2,height = 2,width = 2)
#' 
#' # 提取热图的基因。
#' clusters <- cutree(p$tree_row, k = 3)
#' table(clusters)
#' genes <- names(clusters[clusters == 1])
#' dbs1 <- subset(dbs,subset = dbs$gene %in% genes)
#' write.table(genes,
#'             quote = F,
#'             row.names = F,
#'             col.names = F)

#####NMF非负矩阵分解-----
######按照线粒体基因进行非负矩阵分解分群------
####按照线粒体基因进行降维
library("RColorBrewer")
library(readr)
library(NMF)
library(Seurat)
data.metacell <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
gene_mito <- unique(dbs$gene)
DefaultAssay(data.metacell) <- "RNA" ##设置为integrated可理解为用整合后的data矩阵做下游分析，非原始值。

# # seurat.obj<-NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# scale.genes <-  rownames(data.metacell)
# data.metacell1 <- ScaleData(data.metacell, features = scale.genes)
# data.metacell1 <- FindVariableFeatures(data.metacell1)
# data.metacell1 <- RunPCA(object = seurat.obj1, pc.genes = gene_mito)
# # 肘部图确定最佳dims
# ElbowPlot(seurat.obj1, ndims = 50)
# seurat.obj1 <- FindNeighbors(seurat.obj1, 
#                              dims = 1:10)
# seurat.obj1 <- FindClusters(seurat.obj1, 
#                             resolution = 0.01)
sub_metacell.obj <- subset(data.metacell,subset = cell_type %in% "Cardiomyocyte")
# 
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
mito_gene <- unique(dbs$gene)
mito_gene <- subset(mito_gene,subset = mito_gene %in% rownames(sub_metacell.obj))
DefaultAssay(sub_metacell.obj) <- "RNA" ##设置为integrated可理解为用整合后的data矩阵做下游分析，非原始值。
sub_metacell.obj@assays$RNA@counts <- sub_metacell.obj@assays$RNA@counts[mito_gene, ]
sub_metacell.obj@assays$RNA@data <- sub_metacell.obj@assays$RNA@data[mito_gene, ]
sub_metacell.obj@assays$RNA@meta.features <- sub_metacell.obj@assays$RNA@meta.features[mito_gene,]
# sub_metacell.obj<-NormalizeData(sub_metacell.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# scale.genes <-  rownames(sub_metacell.obj)
# sub_metacell.obj <- ScaleData(sub_metacell.obj, features = scale.genes)
# sub_metacell.obj <- FindVariableFeatures(sub_metacell.obj)
sub_metacell.obj1 <-  CreateSeuratObject(counts = sub_metacell.obj@assays[["RNA"]]@counts, project = "metacell")
sub_metacell.obj1 <- NormalizeData(sub_metacell.obj1) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
sub_metacell.obj1@meta.data <- sub_metacell.obj@meta.data
vm <- sub_metacell.obj1@assays$RNA@scale.data
saveRDS(vm, file = "~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/vm.rds")
res <- nmf(vm, 3)  #很慢
save(res, file = "~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/nmf_res.rda")
load("~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/nmf_res.rda")
# 每个因子提取30个
fs <- extractFeatures(res, 50L)
fs <- lapply(fs, function(x) rownames(res)[x])
fs <- do.call("rbind", fs)
rownames(fs) <- paste0("cluster", 1:3)
write.csv(t(fs), "~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/pb,c_NMF_TopGenes.csv")
DT::datatable(t(fs))
### 选择用于后续分析的因子
s.f <- 1:3   # 因子 1 主要是线粒体和核糖体

## 降维
library(tidyverse)
cell1 <- colnames(sub_metacell.obj)
cell2 <- colnames(coef(res))
cells <- intersect(cell1, cell2)
sub_metacell.obj1 <- sub_metacell.obj1[,cells]
sub_metacell.obj1 <- ScaleData(sub_metacell.obj1)
sub_metacell.obj1 <- FindVariableFeatures(sub_metacell.obj1)
sub_metacell.obj1 <- RunPCA(object = sub_metacell.obj1, pc.genes = mito_gene)
sub_metacell.obj1@reductions$nmf <- sub_metacell.obj1@reductions$pca
sub_metacell.obj1@reductions$nmf@cell.embeddings <- t(coef(res)[,cells])
sub_metacell.obj1@reductions$nmf@feature.loadings <- basis(res)
sub_metacell.obj1 <- RunUMAP(sub_metacell.obj1, reduction='nmf', dims=s.f)
library(Seurat)
library(Matrix)
## 基于NMF降维矩阵的聚类
sub_metacell.obj1 <- FindNeighbors(sub_metacell.obj1, reduction='nmf', dims=s.f) %>% FindClusters()

## 基于因子最大载荷分类
sub_metacell.obj$cluster <- apply(NMF::coefficients(res)[s.f,], 2, which.max)
####聚类降维结果可视化----
p1 <- DimPlot(sub_metacell.obj, label = T) + ggtitle("Clustered by Louvain")
p2 <- DimPlot(sub_metacell.obj, group.by = "cluster", label = T) + ggtitle("Clustered by max loading")
pc <- p1|p2
ggsave(paste0(output,"sub_metacell.obj_NMF_Cluster.pdf"), pc, width = 10, height = 5)
p2
pc
saveRDS(sub_metacell.obj, file = "~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj.rds")
####NMF与心肌细胞亚群之间的相关性---
####基因降维热图-----
####基因聚类----
library(ggplot2)
seurat.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj.rds")
library(Seurat)
# data.metacell <- qs::qread("/home/gongfengcz/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
# sub_metacell.obj <- subset(data.metacell,subset = cell_type %in% "Cardiomyocyte")
# meta <- sub_metacell.obj@meta.data
# meta$cell <- rownames(meta) 
# meta1 <- seurat.obj@meta.data
# meta1$cell <- rownames(meta1)
# meta1 <- meta1[,c(169,170)]
# meta <- merge(meta,meta1,by = "cell")
# rownames(meta) <- meta$cell
# sub_metacell.obj@meta.data <- meta
library(tidyverse)
features <- mito_gene
# sub_metacell.obj1 <- readRDS("~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj1.rds")
sub_metacell.obj <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")
# sub_metacell.obj$cluster <- ""
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 1] <- "C1"
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 2] <- "C2"
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 3] <- "C3"
p2 <- DimPlot(sub_metacell.obj, group.by = "cluster", label = T) + ggtitle("Clustered by max loading")
p2
ggsave(paste0(output,"sub_metacell.obj_NMF_Cluster.pdf"), p2, width = 4.81, height = 5)
# sub_metacell.obj$cluster[sub_metacell.obj$cluster == 4] <- "C4"
all_expr <-
  AverageExpression(sub_metacell.obj,
                    assays = "RNA",
                    slot = "data",
                    group.by = "sub_cell_type")[["RNA"]]
all_expr <- all_expr[!is.infinite(rowSums(all_expr)),]
NMF_expr <-
  AverageExpression(
    sub_metacell.obj,
    assays = "RNA",
    slot = "data",
    group.by = "cluster"
  )[["RNA"]]
all_expr <- as.data.frame(all_expr)
NMF_expr <- as.data.frame(NMF_expr)
all_expr %>% mutate(across(where(is.character), as.numeric))  -> all_expr
NMF_expr %>% mutate(across(where(is.character), as.numeric))  -> NMF_expr
expr <- cbind(all_expr,NMF_expr)
expr %>% mutate(across(where(is.character), as.numeric))  -> expr
library(RColorBrewer)
bk <- c(seq(-1, 1, by = 0.01))
p9 <- cor(expr[,c(1:3)], expr[,c(4:6)]) %>%
  pheatmap::pheatmap(
    cellwidth = 8,
    cellheight = 8,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 8,
    cluster_rows = F,
    cluster_cols = F,
    scale = "row",
    color=colorRampPalette(c("#1E3163","navy", "#00C1D4", "#FFED99","#FF7600","firebrick3"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(-1, 1, 0.5),
    breaks = bk,
    main = "Correlation"
  )
p9
ggsave(paste0(output,"NMF相关性热图.pdf"), plot = p9, width = 4.49, height = 3.67)
#####NMF亚群比例----
#######细胞比例的桑吉图------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
Ratio <- sub_metacell.obj@meta.data %>%
  group_by(group, cluster) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c('#efb306',
            '#7db954',
            '#852f88',
            'firebrick',
            '#0f8096',
            'pink',
            'green')

p10 <- ggplot(Ratio, aes(x =group, y= relative_freq, fill = cluster,
                         stratum=cluster, alluvium=cluster)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Group',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p10
ggsave(paste0(output,"NMF比例.pdf"), plot = p10, width = 4.49, height = 3.67)
#####对c1进行功能富集-----Metascape----
##找到marker基因-----
Idents(sub_metacell.obj) <- "cluster"
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
markers <- FindAllMarkers(sub_metacell.obj,
                       min_pct = 0,
                       logfc.threshold = 0)
write.csv(markers,file = paste0(output,"marker_NMFcluster.csv"))
markers1 <- subset(markers,subset  = markers$cluster %in% "C1")
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
markers1 <-
  markers1[order(markers1$avg_log2FC, decreasing = T),] %>% pull("gene")

C1_bp <-
  enrichGO(
    markers1,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term1 <- C1_bp@result
term1$celltype <- "Cluster I"
#提取
subset_up_term <- subset(term1,subset = term1$Description %in% c("mitochondrial calcium ion homeostasis","mitochondrial calcium ion transmembrane transport","calcium import into the mitochondrion","mitochondrial fission","reactive oxygen species metabolic process","reactive oxygen species biosynthetic process","regulation of reactive oxygen species metabolic process"))
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
ggsave(paste0(output,"NMFC1_功能富集.pdf"), plot = p6, width = 4.45, height = 4.05) 






#####打分的小提琴图-----
# sub_metacell.obj <- readRDS("~/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj1.rds")
meta <- sub_metacell.obj@meta.data
colnames(meta)
meta1 <- meta[,c(167,120,24,25,26,63,99,144,135)]
library(reshape2)
data_long <- melt(meta1,
                  id.vars = c('cluster'),#需要保留不参与聚合的变量,
                  variable.name='Pathways',
                  value.name='Score')
                                                            
library(ggpubr)
my_comparisons <- c("C1","C2","C3")
###Calcium_uniporter通路
a1 <- subset(data_long,subset = Pathways %in% "ROS.and.glutathione.metabolism_UCell") 
comparisons <- list(c("C1","C2"),
                    c("C2","C3"),
                    c("C1","C3"))
a1$cluster <- factor(a1$cluster, levels = c("C1","C2","C3"))
p <- ggviolin(a1, 
              "cluster",
              "Score",
              color = "cluster",
              add = "boxplot", 
              palette = c("#efb306",  "#5797bc", "#e44349"),
              add.params = list(fill = "white"))+
  theme_classic2()+
  stat_compare_means(method = "t.test", 
                     label = "p.signif",##星号设置
                     comparisons = comparisons)+
  stat_compare_means(label.y = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cluster")+ylab("Calcium uniporter")+
  theme(plot.title = element_text(hjust = 0.3),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black") )
p
ggsave(file.path(output, "NMF小提琴calcium.pdf"),
       p,width = 4,height = 3.5)
####FAO----
a1 <- subset(data_long,subset = Pathways %in% "TCA.cycle_UCell") 
comparisons <- list(c("C1","C2"),
                    c("C2","C3"),
                    c("C1","C3"))
a1$cluster <- factor(a1$cluster, levels = c("C1","C2","C3"))
p <- ggviolin(a1, 
              "cluster",
              "Score",
              color = "cluster",
              add = "boxplot", 
              palette = c("#efb306",  "#5797bc", "#e44349"),
              add.params = list(fill = "white"))+
  theme_classic2()+
  stat_compare_means(method = "t.test", 
                     label = "p.signif",##星号设置
                     comparisons = comparisons)+
  stat_compare_means(label.y = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cluster")+ylab("TCA cycle")+
  theme(plot.title = element_text(hjust = 0.3),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black") ) 
p
ggsave(file.path(output, "NMF小提琴ROS.pdf"),
       p,width = 4,height = 3.5)
#######不同cluster功能的热图------
sub_metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/sub_seurat_obj.rds")
sub_NMF <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj.rds")
sub_metacell.obj <- NormalizeData(sub_metacell.obj)
pathway <- read_excel("~/scRNA-heart-mitochodria/data/心肌细胞功能.xlsx")
table(pathway$pathways)
features1 <- subset(pathway,subset = pathway$pathways %in% "GOBP_CARDIAC_CONDUCTION") %>% pull(gene)

expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
features <- c("IL12B","IL17B","IL2","IL20","IL23A","IL25","IL4","LTA","MDK","LIAS","KLRG1",  ###inflammatory
              "CLASP1","CLASP2","HAS3","NTN4","SMPD3","THSD4",    ####ECM
              "ACE2","ATP1A2","ATP2B1","ATP2B2","CAMK2D","CTNNA3","CXADR","ISL1","SCN1B",  ####CARDIAC_CONTRACTION
              "CRK","EID2","FNTA","GDF10","GIPC1","HDAC2","JUN","SMAD5"     #TGFb
)
expr1 <- expr[features,]
expr1 <- t(expr1)
expr1 <- as.data.frame(expr1)
expr1$cell <- rownames(expr1)
meta1 <- sub_NMF@meta.data
meta1$cell <- rownames(meta1)
adata <- merge(meta1,expr1,by = "cell")
colnames(adata)
adata <- adata[,c(168,169:202)]
adata1 <- aggregate(adata, by=list(cluster=adata$cluster),mean)
rownames(adata1) <- adata1[,1]
colnames(adata1)
adata1 <- adata1[,-c(1,2)]
rownames(adata1) <- c("C1","C2","C3")
bk <- c(seq(-1, 1, by = 0.01))
library(pheatmap)
p<-pheatmap(
  adata1,
  scale = "column",#'none', 'row' or 'column'
  # annotation_col=meta,
  color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white"))(50),
            colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),#红到蓝梯度
  # breaks=seq(-value.max, value.max, length.out = 100),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white",
  cellheight = 5,
  cellwidth = 10,
  fontsize = 5,
  angle_col = "45",
  main = 'Genes',
  cluster_cols = F,
  cluster_rows = F,
  legend_breaks = seq(-1, 1, 0.5)
)
####转录因子可视化----
#######转录因子分析 ------
library(readr)
library(SCENIC)
sub_metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/results/Cardiomyocyte/NMF/seurat_obj.rds")
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 1] <- "C1"
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 2] <- "C2"
sub_metacell.obj$cluster[sub_metacell.obj$cluster == 3] <- "C3"
tfs_targer <- read.csv("~/scRNA-heart-mitochodria/results/pyscenic/Cardiomyocytesub_metacell/tfs_targets.csv")
######转录因子可视化  #####
# 导入auc数据
# 导入auc数据
library(AUCell)
library(pheatmap)
regulonAUC <- importAUCfromText("~/scRNA-heart-mitochodria/results/pyscenic/Cardiomyocytesub_metacell/auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- sub_metacell.obj@meta.data
Idents(sub_metacell.obj) <- sub_metacell.obj$cluster
# 将auc得分矩阵添加进seurat
counts = t(getAUC(regulonAUC))
counts <- as.data.frame(counts)
counts$cell <- rownames(counts)
meta <- sub_metacell.obj@meta.data
meta$cell <- rownames(meta)
tf <- merge(meta,counts,by = "cell")
colnames(tf)
tf <- tf[,c(168,169:205)]
data1 <- aggregate(tf, by=list(cluster=tf$cluster),mean)##求均值
rownames(data1) <- data1[,1]
data1 <- data1[,-c(1,2)]
data1 <- t(data1)
library(readxl)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
tfs_targer1 <- subset(tfs_targer,subset = target %in% dbs$gene)
data1 <- as.data.frame(data1)
data1 <- subset(data1,subset = rownames(data1) %in% tfs_targer1$tfs)
p22 <-pheatmap(
  data1,
  scale = "row",#'none', 'row' or 'column'
  # annotation_col=meta,
  color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white"))(50),
            colorRampPalette(colors = c("white","#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
  # breaks=seq(-value.max, value.max, length.out = 100),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white",
  cellheight = 5,
  cellwidth = 10,
  fontsize = 5,
  angle_col = "45",
  main = 'Transcription Factor',
  cluster_cols = F,
  cluster_rows = T
)
print(p22)
ggsave(paste0(output,"转录因子热图.pdf"), plot = p22, width = 6.35, height = 4.29)
#####相关性分析-----
counts = t(getAUC(regulonAUC))
rownames(counts) <-gsub(",", "", rownames(counts))
counts <- t(counts)
counts <- counts[,c("EGR1(+)","CUX1(+)","SREBF2(+)","POUSF1(+)","MEF2D(+)","MEF2A(+)")]
target <- subset(tfs_targer1,subset = tfs %in% c("EGR1(+)","CUX1(+)","SREBF2(+)","POUSF1(+)","MEF2D(+)","MEF2A(+)"))
target <- unique(target$target)
data.metacell <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
data.metacell1 <- subset(data.metacell,subset = cell_type %in% "Cardiomyocyte")
expr <- as.data.frame(data.metacell1[["RNA"]]@data)
expr <- t(expr)
target <- c("MGST3","IDI1","MRPL33","PGS1","TIMM23","MRPL52","MUTYH","GRPEL2","TXN2","MCL1","ISCU")
expr <- expr[,target]
adata <- cbind(counts,expr)
meta <- data.metacell1@meta.data
meta <- meta[,c(1,2)]
adata <- cbind(adata,meta)
####相关性热图绘制
library(readxl)
library(UCell)
library(psych)
adata <- aggregate(adata, by=list(saple=adata$sample),mean)
colnames(adata)
cor_function_meta <- corr.test(adata[,c(9,10,22,23,32,28)],adata[,c(39:49)], method = "spearman", adjust = "fdr")#做相关性


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


p13 <- ggplot(allnew, aes(variable,names)) +
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
ggsave(paste0(output,"转录因子相关性热图.pdf"), plot = p13, width = 11.03, height = 3.81)
#####相关性散点图------
########线粒体相关功能与炎症功能
library(UCell)
library(readxl)
library(readr)
metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM-NF.rds")
metacell.obj <- NormalizeData(metacell.obj)
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
dbs2 <- read_excel("~/scRNA-heart-mitochodria/data/巨噬细胞下游功能通路.xlsx")
colnames(dbs2) <- c("pathway","gene")
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs <- rbind(dbs,dbs2)
dbs1 <- split(dbs$gene, dbs$pathway)
sub_metacell.obj <- AddModuleScore_UCell(sub_metacell.obj, features = dbs1,
                                       ncores = 50)
# data.HCM <- AddModuleScore_UCell(HCM, features = dbs1,
#                                  ncores = 20)
library(vegan)
#install.packages('vegan')
library(psych)
library(reshape2)
library(ggplot2)
#data("varechem")
library(Seurat)
adata <- sub_metacell.obj@meta.data
###提取出来之后把表达为0的值去掉-----
####相关性散点图-----
###相关性散点图----
###散点图----
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
colnames(adata)
data <- adata[,c(1,68,132,85,86,87)]
data1 <- data
data1[data1 == 0] <- NA
data1 <- na.omit(data1)
data1 <- aggregate(data1, by=list(sample=data1$biosample_id),mean)##求均值
colnames(data1)
p14 <- data1 %>%
  ggstatsplot::ggscatterstats( 
    x = Fatty.acid.oxidation_UCell,
    y = GOBP_INFLAMMATORY_RESPONSE_UCell,
    type =  "pearson",
    conf.level = 0.99, # confidence level 
    line.color = "#FFA300",
    bf.message = FALSE,  ##去除贝叶斯相关的统计值
    ggtheme = theme_bw(), # choosing a different theme
    ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
    xalpha = 0.6,
    yalpha = 0.6 ,
    centrality.para = "median",
    messages = F, # turn off messages and notes 
    title = "Relationship with FAO")
p14
ggsave(paste0(output,"FAO_炎症相关散点.pdf"), plot = p14, width = 6.41, height = 5.59)
p15 <- data1 %>%
  ggstatsplot::ggscatterstats( 
    x = OXPHOS_UCell,
    y = GOBP_INFLAMMATORY_RESPONSE_UCell,
    type =  "pearson",
    conf.level = 0.99, # confidence level 
    line.color = "#FFA300",
    bf.message = FALSE,  ##去除贝叶斯相关的统计值
    ggtheme = theme_bw(), # choosing a different theme
    ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
    xalpha = 0.6,
    yalpha = 0.6 ,
    centrality.para = "median",
    messages = F, # turn off messages and notes 
    title = "Relationship with OXPHOS")
p15
ggsave(paste0(output,"OXPHOS_炎症相关散点.pdf"), plot = p15, width = 6.41, height = 5.59)


######clusterprofile转录因子靶基因功能富集网络图-------
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
target <- tfs_targer
SREBF2 <- subset(target,subset = tfs %in% "SREBF2(+)") %>% pull(target)
POU2F1 <- subset(target,subset = tfs %in% "POU2F1(+)") %>% pull(target)
MEF2D <- subset(target,subset = tfs %in% "MEF2D(+)") %>% pull(target)
MEF2A <- subset(target,subset = tfs %in% "MEF2A(+)") %>% pull(target)
EGR1 <- subset(target,subset = tfs %in% "EGR1(+)") %>% pull(target)
CUX1 <- subset(target,subset = tfs %in% "CUX1(+)") %>% pull(target)

data <- list(SREBF2,POU2F1,MEF2D,MEF2A,EGR1,CUX1)
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
gene_cluster <- data
names(gene_cluster)=c("SREBF2","POU2F1","MEF2D","MEF2A","EGR1","CUX1")

xx <- compareCluster(gene_cluster, 
                     fun='enrichGO',
                     ont= 'BP',
                     OrgDb='org.Hs.eg.db' ,
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1
)


xx <- pairwise_termsim(xx)
a <- xx@compareClusterResult
a <- subset(a,subset = Description %in% c("fatty acid beta-oxidation","fatty acid catabolic process","fatty acid metabolic process","energy derivation by oxidation of organic compounds","energy homeostasis","ATP metabolic process","ATP biosynthetic process","ATP generation from ADP","cardiac muscle tissue development","ardiac muscle tissue morphogenesis","muscle tissue development","muscle organ development","striated muscle tissue development","skeletal muscle tissue development","skeletal muscle organ development","skeletal muscle cell differentiation","ventricular cardiac muscle tissue morphogenesis","oxidative phosphorylation","respiratory electron transport chain","ATP synthesis coupled electron transport","mitochondrial ATP synthesis coupled electron transport","electron transport chain","ATP synthesis coupled electron transport","negative regulation of transforming growth factor beta receptor signaling pathway"))
xx@compareClusterResult <- a  
pdf("~/scRNA-heart-mitochodria/results/Fibroblast/pyscenic/HCM_NF/go_net2.pdf",height = 8,width = 8)
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
  scale_fill_manual(values = c('#e74c3c',
                               '#3498db',
                               '#9b59b6',
                               "#f1c40f",
                               '#3cb371','#00688b'))
dev.off()
p23
ggsave(paste0(output,"转录因子功能富集网络图.pdf"), plot = p23, width = 6.68, height = 5.20)

###A功能富集和弦图----
###功能富集和弦图-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readxl)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg1 <- subset(deg,subset = deg$`Cell Type` %in% "Cardiomyocyte")
deg1$change <- ""
deg1$change[deg1$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg1$change[deg1$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg1$change[deg1$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg1 <- deg1[,c(1,15,16,17,21)]
colnames(deg1) <- c("Gene","logFC","Pvalue","Adjusted_Pvalue","change")
deg1$log10Pvalue <- log10(deg1$Pvalue)

down_deg <- subset(deg1,subset = Pvalue < 0.05 & logFC < -0.25)
down_deg$gene <- rownames(down_deg)
down <-
    down_deg[order(down_deg$logFC, decreasing = T),] %>% pull("Gene")

down_bp <-
  enrichGO(
    down,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
down_term <- down_bp@result
# up_term$celltype <- "Activated_fibroblast"

####图形绘制
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)

go <- data.frame(Category = "BP",
                 ID = down_term$ID,
                 Term = down_term$Description, 
                 Genes = gsub("/", ", ", down_term$geneID), 
                 adj_pval = down_term$p.adjust)
genelist <- data.frame(ID = down_deg$Gene, logFC = down_deg$logFC)
circ <- circle_dat(go, genelist)
pathways <- c("fatty acid metabolic process","fatty acid catabolic process","ATP metabolic process","oxidative phosphorylation")
head(circ)
n = 4 #圈图需要选定term，这里画前面5个

chord <- chord_dat(circ, genelist, pathways)
head(chord)
p3 <- GOChord(chord, 
              space = 0.02, # Gene block spacing
              gene.order = 'logFC', 
              # nlfc = ,
              lfc.min = -2,
              lfc.max = 0,
              lfc.col = c('#4CAF50', "white", '#006994'), # Adjusted colors for higher contrast
              gene.space = 0.26, # Distance of gene names from the circle
              gene.size = 5, # Font size of gene names 
              border.size = 0.1, # Thickness of the black border on the curves
              process.label = 8) # Font size of term labels

p3
ggsave(filename = file.path(output, "成纤维整体差异功能富集和弦图.pdf"),
       plot = p3,
       height = 8,
       width = 8)
