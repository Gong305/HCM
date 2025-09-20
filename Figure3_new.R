###########10月11日NC-----
setwd("~/scRNA-heart-mitochodria")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure2/311/"
dir.create(output, recursive = TRUE)
####Ucell打分----
library(tidyr)
library(dplyr)
library(Seurat)
######线粒体基因降维图----
#####本身降维信息-----
library(ggrepel)
library(AUCell)
library(clusterProfiler)
library(tidyverse)
library(msigdbr)
HCM <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
HCM <- NormalizeData(HCM)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
library(UCell)
data.seurat <- AddModuleScore_UCell(HCM, features = dbs1,maxRank = 15000,
                                    ncores = 50)
# saveRDS(data.seurat,file = "~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
data.seurat <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
meta<- data.seurat@meta.data
meta <- meta[, c("disease","cell_type",
                 colnames(meta)[grep("_UCell", colnames(meta))])]
t_results <- meta
# 保存t值结果
saveRDS(t_results,
        file =
          "~/scRNA-heart-mitochodria/results/整体/ALL_UCscore.rds")
t_results <- readRDS("~/scRNA-heart-mitochodria/results/整体/ALL_UCscore.rds")
#### 可视化 ----
colnames(t_results) <- gsub("_UCell", "", colnames(t_results))
colnames(t_results) <- gsub("_", " ", colnames(t_results))
t_results$cell_id <- rownames(t_results)
es.matrix <- t_results[,-c(1,2)]
es.matrix <- t(es.matrix)
#### 差异分析 ----
cell_types <- unique(t_results$`cell type`)
cell_types <- cell_types[-which(cell_types == "Epicardial cell")]
library(limma)
t_results1 <- vector("double", 5) 
for (i in cell_types) {
  case <- "HCM"
  control <- "NF"
  meta <-
    data.seurat@meta.data[, c("cell_type",
                              "disease")]
  meta$cell_id <- rownames(meta)
  cell_id.1 <-
    meta[c(meta$disease == case&meta$cell_type == i), ]$cell_id
  cell_id.2 <-
    meta[c(meta$disease == control&meta$cell_type == i), ]$cell_id
  
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
  
  
  
  # colnames(t_results) <- c("t", "Gobp")
  
  # t_results$Gobp = with(t_results, reorder(Gobp, t))
  t_results <- as.data.frame(t_results)
  t_results$fill <- ""
  t_results <- t_results[-150,]
  t_results$t <- as.numeric(t_results$t)
  if (any(t_results$t >= 2.58)) {
    t_results[t_results$t >= 2.58, ]$fill <- "up"
  }
  if (any(t_results$t <= -2.58)) {
    t_results[t_results$t <= -2.58,]$fill <- "down"
  }
  if (any(abs(t_results$t) < 2.58)) {
    t_results[abs(t_results$t) < 2.58,]$fill <-
      "no"
  }
  
  t_results$color <- ""
  if (any(abs(t_results$t) < 2.58)) {
    t_results[abs(t_results$t) < 2.58,]$color <-
      "n"
  }
  if (any(abs(t_results$t) >= 2.58)) {
    t_results[abs(t_results$t) >= 2.58,]$color <-
      "y"
  }
  t_results$feature <- rownames(t_results)
  t_results$cell_type <- i
  t_results1 <- rbind(t_results1,t_results)
}
t_results1 <- t_results1[-1,]
saveRDS(t_results1,file = "~/scRNA-heart-mitochodria/results/整体/all_celltype.rds")
t_results1 <- readRDS("~/scRNA-heart-mitochodria/results/整体/all_celltype.rds")
####大类对应----
library(readxl)
category <- read_excel("~/scRNA-heart-mitochodria/data/线粒体通路大类.xlsx")
category$MitoPathway <- gsub(" ", ".", category$MitoPathway)
category$MitoPathway <- gsub(",", ".", category$MitoPathway)
category$MitoPathway <- gsub("-", ".", category$MitoPathway)
results <- merge(t_results1,category,by.x = "feature",by.y = "MitoPathway",all.x = TRUE)
results1 <- subset(results,subset = color %in% "y")
#####饼状图绘制-----
df1 <- table(results1$`MitoPathways Hierarchy`)
df1 <- as.data.frame(df1)
df1$Percentage <- (df1$Freq / sum(df1$Freq)) * 100
library(ggforce)
library(dplyr)
library(forcats)
df2 <- df1 %>%
  mutate(perc = Freq/ sum(Freq)) %>%
  mutate(labels = scales::percent(perc)) %>%
  arrange(desc(perc)) %>%
  mutate(Var1 = fct_rev(fct_inorder(Var1))) %>%
  mutate(text_y = cumsum(Freq) - Freq/2)

df2
library(ggrepel)
library(ggsci)
library(paletteer)
colors <- c("#80c97f","#a68dc8","#7F7F7FFF","#8C564BFF",
            "#ffc000","#c00000",'#E5D2DD', '#F1BB72', '#F3B1A0', 
            '#D6E7A3', '#57C3F3','#476D87',"darkorange1", 
            "darksalmon","#DBDB8DFF")
p1 <- df2 %>%
  ggplot(aes(x = "", y = Freq, fill = Var1)) +
  # draw_key_point 使用一个点作为图例中每个类别的显示标记
  # color 增加白色边框
  geom_col(key_glyph = draw_key_point, color = "white") +
  coord_polar(theta = "y") +
  geom_label_repel(aes(label = labels, y = text_y),
                   nudge_x = 0.6,
                   nudge_y = 0.6,
                   size = 5,
                   show.legend = F,
                   segment.color = "grey50"
  ) +
  # override.aes 修改图例形状和大小
  guides(fill = guide_legend(title = "Category",
                             ncol = 1, # 图例3列
                             override.aes = list(shape = 21, size=5))) +
  scale_fill_manual(values = colors) +
  theme_void() +
  # 修改图例位置
  theme(legend.position = c(1.1,0.5),
        # 图例折叠为3列，增加边框
        legend.background = element_rect(fill = "white",
                                         size = 0.3,
                                         linetype="solid"))
ggsave(paste0(output,"通路紊乱饼图.pdf"), plot = p1, width = 10.23, height = 8.49)

#####细胞类型数量的柱状图-----
####数量统计-----
data <- results1
library(ggsci)
number <- data %>%
  group_by(feature) %>%
  summarize(NumYInCellTypes = sum(color == 'y'))
number$group <- ""
number$group[number$NumYInCellTypes == 0] <- "None"
number$group[number$NumYInCellTypes >=1 & number$NumYInCellTypes <=3] <- "1-3 celltypes"
number$group[number$NumYInCellTypes >=4 & number$NumYInCellTypes <=6] <- "4-6 celltypes"
number$group[number$NumYInCellTypes >=7 & number$NumYInCellTypes <=9] <- "7-9 celltypes"
number$group[number$NumYInCellTypes >=9 & number$NumYInCellTypes <=11] <- "9-11 celltypes"
number$group[number$NumYInCellTypes == 12] <- "All celltypes"
data1 <- table(number$group)
data1 <- as.data.frame(data1)

results2 <- subset(results,subset = color %in% "n")
data2 <- results2
number2 <- data2 %>%
  group_by(feature) %>%
  summarize(NumYInCellTypes = sum(color == 'n'))
data2 <- data.frame(
  Var1 = "None",
  Freq = 1
)
data <- rbind(data2,data1)
data$label <- rownames(data)
data$label <- as.numeric(data$label)
colnames(data) <- c("Group","Freq","label")
p1 <- ggplot(data,aes(reorder(Group,label), Freq,fill = Group,width = 0.8))  + theme(panel.grid.major=element_line(colour=NA),
                                                                                     panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                                                     panel.grid.minor = element_blank(),
                                                                                     axis.line = element_line(colour = "black"))+ labs(title = "Number of mitochondrial reated pathways", x = "Group" , y = "Number of statistically significant dysregulated pathways") + geom_bar(stat = 'identity') + geom_text(aes(label = Freq), size = 5, vjust = -0.3) + scale_y_continuous(expand = c(0,0), limits = c(0,75)) ## expand修改数据条的位置，limit修改y轴的标签限制

p2 <- p1 + scale_fill_jco()
ggsave(paste0(output,"Figure2B.pdf"), plot = p2, width = 6.95, height = 4.49)
library(reshape2)
######通路差异紊乱的热图------
#####脂类代谢-----
data <- results1
data1 <- subset(data,subset = data$`MitoPathways Hierarchy` %in% "Lipid metabolism")

data1 <- data1[,c(1,2,5)]
data_wide<-dcast(data1, cell_type~feature,
                 value.var = 't')
for (i in 1:ncol(data_wide)){
  data_wide[,i][is.na(data_wide[,i])] <- 0
}
rownames(data_wide) <- data_wide[,1]
data_wide <- data_wide[,-1]
data_wide %>% mutate(across(where(is.character), as.numeric))  -> data_wide
library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white"))(1000)

# breaks = seq(min(unlist(c(RPF3))), max(unlist(c(RPF3))), length.out=100)
# pheatmap(RPF3, color=colormap, breaks=breaks,border_color="black",
#          cutree_col = 2,cutree_row = 4)
breaks = seq(-10, 0, length.out=1000)
data_wide <- t(data_wide)
data_wide <- data_wide[c(1,5,6,8,9),c(1,2,3,4,5,7,8,10,11,12)]
rownames(data_wide) <- sub("\\.", " ", rownames(data_wide))
pheatmap(data_wide,
         cellheight=40, cellwidth=20,
         height = 10, width = 10,
         number_color="red", 
         number_format="%.2e",
         border="white",
         fontsize_number = 10, 
         fontsize = 10,
         show_rownames = T,
         scale = "none", # column
         main="Lipid metabolism",
         clustering_distance_rows = "minkowski",
         clustering_method="complete",
         cluster_cols = F, treeheight_col = 20,
         cluster_rows = F, treeheight_row = 20,
         color=colormap, breaks=breaks)

dev.off()
#####碳代谢-------
######通路差异紊乱的热图------
data <- results1
data1 <- subset(data,subset = data$`MitoPathways Hierarchy` %in% "Carbohydrate metabolism")
data1 <- data1[,c(1,2,5)]
data_wide<-dcast(data1, cell_type~feature,
                 value.var = 't')
for (i in 1:ncol(data_wide)){
  data_wide[,i][is.na(data_wide[,i])] <- 0
}
rownames(data_wide) <- data_wide[,1]
data_wide <- data_wide[,-1]
data_wide %>% mutate(across(where(is.character), as.numeric))  -> data_wide
library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(1000)
breaks = seq(-10, 10, length.out=1000)
data_wide <- t(data_wide)
data_wide <- data_wide[c(1,5,7,8,10),]
pheatmap(data_wide, 
         cellheight=20,cellwidth=10,
         height = 10,width = 10,
         number_color="red", 
         number_format="%.2e",
         border="white",
         fontsize_number = 10, 
         fontsize = 10,
         show_rownames = T,
         scale = "none",#column
         main="Carbohydrate metabolism",angle_col = 90,
         clustering_distance_rows = "minkowski",
         clustering_method="complete",
         cluster_cols = F,treeheight_col = 20,
         cluster_rows = F,treeheight_row = 20,
         color=colormap,breaks=breaks)
dev.off()
#######代谢----
#####碳代谢-------
######通路差异紊乱的热图------
data <- results1
data1 <- subset(data,subset = data$`MitoPathways Hierarchy` %in% "OXPHOS")
data1 <- data1[,c(1,2,5)]
data_wide<-dcast(data1, cell_type~feature,
                 value.var = 't')
for (i in 1:ncol(data_wide)){
  data_wide[,i][is.na(data_wide[,i])] <- 0
}
rownames(data_wide) <- data_wide[,1]
data_wide <- data_wide[,-1]
data_wide %>% mutate(across(where(is.character), as.numeric))  -> data_wide
library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(1000)
breaks = seq(-10, 10, length.out=1000)
data_wide <- t(data_wide)
data_wide <- data_wide[c(2,4,5,6,10,11,16,20),c(1,2,3,4,5,6,7,8,10,11,12)]
pheatmap(data_wide, 
         cellheight=20,cellwidth=10,
         height = 10,width = 10,
         number_color="red", 
         number_format="%.2e",
         border="white",
         fontsize_number = 10, 
         fontsize = 10,
         show_rownames = T,
         scale = "none",#column
         main="OXPHOS",angle_col = 90,
         clustering_distance_rows = "minkowski",
         clustering_method="complete",
         cluster_cols = F,treeheight_col = 20,
         cluster_rows = F,treeheight_row = 20,
         color=colormap,breaks=breaks)
dev.off()
# #########每个细胞类型差异基因的热图 ------
# ######差异基因计算------
# #### 寻找差异基因 ----
# # 加了一列细胞类型及分组信息
# Idents(HCM) <- HCM$cell_type
# HCM$celltype.group <-
#   paste(Idents(HCM), HCM$disease, sep = "_")
# # 加了一列细胞类型
# HCM$celltype <- Idents(HCM)
# # 将细胞类型及刺激状态作为分组
# Idents(HCM) <- "celltype.group"
# 
# # 查看当前分组
# table(Idents(HCM))
# 
# #使用FindMarkers函数寻找差异表达基因
# DefaultAssay(HCM) <- "RNA"
# celltypes <- unique(HCM$cell_type)
# 
# markers1 <- FindMarkers(HCM,
#                         ident.1 =paste0("Cardiomyocyte_HCM") ,
#                         ident.2 = paste0("Cardiomyocyte_NF"))
# markers1$celltype <- "Cardiomyocyte"
# markers1$gene <- rownames(markers1)
# 
# celltypes <- c("Adipocyte","Macrophage", "Endocardial","Fibroblast","Lymphatic_endothelial","VSMC",                 
#                "Endothelial","Pericyte","Neuronal","Lymphocyte","Mast_cell")
# for (i in celltypes) {
#   markers <- FindMarkers(HCM,
#                          ident.1 =paste0(i,"_HCM") ,
#                          ident.2 = paste0(i,"_NF"))
#   markers$celltype <- i
#   markers$gene <- rownames(markers)
#   markers1 <- rbind(markers,markers1)
#   
# }
# write_csv(markers1,file = "~/scRNA-heart-mitochodria/results/deg/Findmarker计算差异.csv")
#####差异基因提取-----
library(readxl)
library(reshape2)
library(tidyverse)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
HCM_diff <- read_excel("~/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
HCM_diff <- HCM_diff[,c(1,3,15,16,17)]
HCM_diff <- subset(HCM_diff,HCM_diff$Gene %in% dbs$gene)
colnames(HCM_diff) <- c("Gene","cell_type","logFC","pvalue","adjusted_pvalue")
HCM_diff <- subset(HCM_diff,subset = HCM_diff$`adjusted_pvalue` < 0.01)
# HCM_diff <- read_csv("~/scRNA-heart-mitochodria/results/deg/Findmarker计算差异.csv")
# colnames(HCM_diff) <- c("pvalue","logFC","pct.1","pct.2","adjusted_pvalue","cell_type","Gene")
# HCM_diff <- subset(HCM_diff,HCM_diff$Gene %in% dbs$gene)
######通路基因的雷达图-----
####线粒体差异基因玫瑰图展示-----
HCM_diff$group <- ""
HCM_diff$group[HCM_diff$logFC > 0] <- "Up"
HCM_diff$group[HCM_diff$logFC < 0] <- "Down"
a <- table(HCM_diff$cell_type,HCM_diff$group) %>% as.data.frame()
colnames(a) <- c("celltype","group","Freq")
library(Seurat)
library(tidyverse)
library(reshape2)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(readr)

#### /玫瑰图展示差异基因数量 ----
p2 <- ggplot(a, aes(x = celltype, y = Freq)) +
  geom_col(aes(fill = group),
           width = 0.9,
           size = 0,
           alpha = 0.8) +
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5), # 调整数值标签的位置
            size = 3, # 调整标签字体大小
            color = "black")  +
  labs(x = "The number of mitoDEG") +
  coord_polar() +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(
    breaks = c("Up", "Down"),
    values = c("#E64B35FF", "#4DBBD5FF")
  )
p2



#####差异基因的upset图------
#####upset-----
HCM_diff <- read_csv("~/scRNA-heart-mitochodria/results/deg/HCM_NF.csv")
HCM_diff <- read_excel("~/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
HCM_diff <- HCM_diff[,c(1,3,15,16,17)]
HCM_diff <- subset(HCM_diff,HCM_diff$Gene %in% dbs$gene)
colnames(HCM_diff) <- c("Gene","cell_type","logFC","pvalue","adjusted_pvalue")
HCM_diff <- subset(HCM_diff,subset = HCM_diff$`adjusted_pvalue` < 0.01)

dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")

library(gdata)
library(UpSetR)
up <- subset(HCM_diff,subset = logFC > 0 )
up_wide<-dcast(up, Gene~up$cell_type,
               value.var = 'Gene')
up_wide <- up_wide[,-1]
# up$group <- paste0(up$TF,up$celltype)
# up1 = distinct(up,group,.keep_all = T)
down <- subset(HCM_diff,subset = logFC < 0 )
down_wide<-dcast(down, Gene~down$cell_type,
                 value.var = 'Gene')
down_wide <- down_wide[,-1]
# down1 = distinct(down,group,.keep_all = T)
# write_csv(up1,file = "~/scRNA-heart-mitochodria/results/scenic1/deg/upall.csv")
# write_csv(down1,file = "~/scRNA-heart-mitochodria/results/scenic1/deg/downall.csv")
# up <- read_csv("~/scRNA-heart-mitochodria/results/scenic1/deg/uptfwide.csv")
# down <- read_csv("~/scRNA-heart-mitochodria/results/scenic1/deg/downtfwide.csv")
#默认
sets <- c("Pseudo-Bulk","Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast","Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal","Pericyte","VSMC")
upset(fromList(down_wide), 
      order.by = "freq",
      nsets = 100,
      nintersects = 30, #需要绘制的交集数目
      mb.ratio = c(0.7, 0.3),#柱状图与矩阵点图之间的比例大小
      # number.angles = 0,#柱子上方数字倾斜角度
      show.numbers = "yes",
      point.size = 1.8,#矩阵中圆圈的大小
      line.size = 0.8, #矩阵点图中点和线的大小
      sets.x.label = "Set Size", #柱状图的轴标签
      #main.bar.color = c("#73BAD6","#73BAD6","#73BAD6","#73BAD6","#73BAD6","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf"), #柱状图柱子颜色  #"#7dc3fe",
      main.bar.color = "#0780cf",
      sets.bar.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#6495ED",
                         "cyan1", "royalblue4", "darksalmon", 
                         "darkgoldenrod1"),
      # matrix.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F"),  #交点的颜色
      matrix.color = "#0780cf",
      mainbar.y.label = "Gene Intersections", 
      text.scale = 1.75,
      shade.color = "#0780cf",
      shade.alpha = 0.1
      
      # mainbar.y.max = 500
      # att.color = sets.bar.color
      
)
upset(fromList(up_wide), 
      order.by = "freq",
      nsets = 100,
      nintersects = 30, #需要绘制的交集数目
      mb.ratio = c(0.7, 0.3),#柱状图与矩阵点图之间的比例大小
      # number.angles = 0,#柱子上方数字倾斜角度
      show.numbers = "yes",
      point.size = 1.8,#矩阵中圆圈的大小
      line.size = 0.8, #矩阵点图中点和线的大小
      sets.x.label = "Set Size", #柱状图的轴标签
      #main.bar.color = c("#73BAD6","#73BAD6","#73BAD6","#73BAD6","#73BAD6","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf","#0780cf"), #柱状图柱子颜色  #"#7dc3fe",
      main.bar.color = "firebrick3",
      sets.bar.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#6495ED",
                         "cyan1"),
      #matrix.color = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#F0027F"),  #交点的颜色
      matrix.color = "firebrick3",
      mainbar.y.label = "Gene Intersections", 
      text.scale = 1.75,
      shade.color = "firebrick3",
      shade.alpha = 0.1
      
      # mainbar.y.max = 500
      # att.color = sets.bar.color
      
)

#####差异基因的气泡图----
library(readxl)
diff <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
HCM_diff <- diff[,c(1,2,3,15,16,17)]
colnames(HCM_diff) <- c("Gene","Ensembl_ID","cell_type","logFC","Pvalue","Adjusted_pvalue")
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
#####FAO---
mito_diff <- subset(HCM_diff,subset = HCM_diff$Gene %in% dbs$gene)
dbs1 <- subset(dbs,subset = dbs$pathway %in% "Fatty acid oxidation")
mito_diff <- subset(HCM_diff,subset = HCM_diff$Gene %in% dbs1$gene)
mito_diff$logFC[mito_diff$logFC >= 2] <- 2
mito_diff$logFC[mito_diff$logFC <= -2] <- -2
mito_diff <- subset(mito_diff,subset = !(Gene %in% c("ACADSB","PCCB","PCCA","MCEE","HSD17B10","HADH","ECHS1","CROT","ACSS1","ACSM5","ACOT11")))
p <- mito_diff %>%
  ggplot(aes(
    x = cell_type,
    y = Gene,  # 不再使用 expression(italic(Gene))
    size = -log10(Pvalue),
    color = logFC
  )) +
  geom_point(size = 3) +
  theme_cat() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "black",
      size = 10,
      # face = "italic"
    ),
    axis.text.y = element_text(
      colour = "black",
      face = "italic",  # 将 y 轴基因名称设置为斜体
      size = 10
    ),
    axis.title = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "pt"),
    legend.key = element_blank(),
    legend.margin = margin(r = 0, unit = "pt"),
    plot.title = element_text(
      size = 3,
      colour = "black",
      face = "plain",
      hjust = 0.5
    )
  ) +
  scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", limits = c(-2, 2)) +
  ggtitle("Fatty acid oxidation") +
  theme(plot.title = element_text(size = 14))  # 调整标题字体大小为14

p
####OXPHOS----
mito_diff <- subset(HCM_diff,subset = HCM_diff$Gene %in% dbs$gene)
dbs1 <- subset(dbs,subset = dbs$pathway %in% "Fatty acid oxidation")
mito_diff <- subset(HCM_diff,subset = HCM_diff$Gene %in% dbs1$gene)
mito_diff$logFC[mito_diff$logFC >= 2] <- 2
mito_diff$logFC[mito_diff$logFC <= -2] <- -2
mito_diff <- subset(mito_diff,subset = !(Gene %in% c("ACADSB","PCCB","PCCA","MCEE","HSD17B10","HADH","ECHS1","CROT","ACSS1","ACSM5","ACOT11")))
p <- mito_diff %>%
  ggplot(aes(
    x = cell_type,
    y = Gene,  # 不再使用 expression(italic(Gene))
    size = -log10(Pvalue),
    color = logFC
  )) +
  geom_point(size = 3) +
  theme_cat() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "black",
      size = 10,
      # face = "italic"
    ),
    axis.text.y = element_text(
      colour = "black",
      face = "italic",  # 将 y 轴基因名称设置为斜体
      size = 10
    ),
    axis.title = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "pt"),
    legend.key = element_blank(),
    legend.margin = margin(r = 0, unit = "pt"),
    plot.title = element_text(
      size = 3,
      colour = "black",
      face = "plain",
      hjust = 0.5
    )
  ) +
  scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", limits = c(-2, 2)) +
  ggtitle("OXPHOS") +
  theme(plot.title = element_text(size = 14))  # 调整标题字体大小为14

p


######差异的转录因子网络图--------
#####网络图准备-----
samples <- list.files("~/scRNA-heart-mitochodria/results/scenic/")
a <- vector("double", 0) 
for (x in samples) {
  # x <- 1
  data_dir <-
    str_c("~/scRNA-heart-mitochodria/results/scenic/",
          x,
          "/up/int/2.5_regulonTargetsInfo.Rds")
  counts <-
    readRDS(data_dir)
  # counts <- subset(counts,subset = counts$highConfAnnot %in% "TRUE")
  counts <- subset(counts,subset = counts$gene %in% dbs$gene)
  counts$celltype <- x
  a <- rbind(a,counts,fill=TRUE)
}
a1 <- a[,-1]
#####下调的-----
samples <- list.files("~/scRNA-heart-mitochodria/results/scenic/")
a <- vector("double", 0) 
for (x in samples) {
  # x <- 1
  data_dir <-
    str_c("~/scRNA-heart-mitochodria/results/scenic/",
          x,
          "/down/int/2.5_regulonTargetsInfo.Rds")
  counts <-
    readRDS(data_dir)
  # counts <- subset(counts,subset = counts$highConfAnnot %in% "TRUE")
  counts <- subset(counts,subset = counts$gene %in% dbs$gene)
  counts$celltype <- x
  a <- rbind(a,counts,fill=TRUE)
}
a2 <- a[,-1]

a1$Gene <- paste0(a1$gene,"_",a1$celltype)
a2$Gene <- paste0(a2$gene,"_",a2$celltype)
a1$group <- "up"
a2$group <- "down"
write_csv(a1,file = "~/scRNA-heart-mitochodria/results/scenic1/deg/up_mitoTF.csv")
write_csv(a2,file = "~/scRNA-heart-mitochodria/results/scenic1/deg/down_mitoTF.csv")
####转到cytoscape做网络图-----
####转录因子热图-----
#####转录因子交集热图-----
a1 <- read_csv("~/scRNA-heart-mitochodria/results/scenic1/deg/up_mitoTF.csv")
a2<- read_csv("~/scRNA-heart-mitochodria/results/scenic1/deg/down_mitoTF.csv")

heatmap <- read_excel("~/scRNA-heart-mitochodria/results/scenic1/deg/mitoTF_heatmap.xlsx")
rownames(heatmap) <- heatmap$'cell type'
heatmap <- as.data.frame(heatmap)
rownames(heatmap) <- heatmap$'cell type'
heatmap <- heatmap[,c(2:8)]
heatmap <- heatmap[-6,]
test2 <- as.matrix(heatmap)

#注释数据框'行名'必须和绘图矩阵'行名'一致
diagnosis <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast","Lymphocyte","Macrophage","Neuronal","Pericyte","VSMC"  )
annotation_row = data.frame(factor(diagnosis, levels = diagnosis))
rownames(annotation_row) = diagnosis
rownames(annotation_row) = factor(rownames(annotation_row) ,levels = diagnosis)
colnames(annotation_row) = 'diagnosis'
ann_colors = list(
  diagnosis = c('Adipocyte'='#0072B2','Cardiomyocyte'='#D55E00','Endocardial'='#F0E442',
                'Endothelial'='#009E73','Fibroblast'='#56B4E9',"Lymphocyte"= "#984EA3","Macrophage"= "#FF34B3","Neuronal"= "#FFFF33","Pericyte" = "#A65628","VSMC" = "#4DAF4A"))
library(pheatmap)
pdf(file.path("~/scRNA-heart-mitochodria/results/scenic1", 'inter_5TF_diagosis_heatmap.pdf'),width = 3,height = 3)
pheatmap(test2,cluster_rows = F,cluster_cols = F,legend = F,scale='none',
         main='UP/DOWN TFs',gaps_col =c(3),show_rownames = F,
         annotation_colors = ann_colors,annotation_row = annotation_row, 
         cellwidth = 12, cellheight = 12, border_color = "grey", fontsize = 8, 
         color = c("#ee6470","white", "#00a6e1"))
dev.off()
