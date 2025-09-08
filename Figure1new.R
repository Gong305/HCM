####设置工作路径-----
setwd("~/scRNA-heart-mitochodria/")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure1/311/"
dir.create(output, recursive = TRUE)
####加载包-----
library(tidyr)
library(dplyr)
library(Seurat)
######线粒体基因降维图----
#####本身降维信息-----
library(ggrepel)
# seurat.obj <- readRDS("/home/lixy/SCP1303/data/diseasetype/HCMvsNF.rds")
# table(seurat.obj$cell_type_leiden0.6)
# seurat.obj$cell_type <- ""
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Activated_fibroblast"] <- "Fibroblast"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Fibroblast_I"] <- "Fibroblast"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Fibroblast_II"] <- "Fibroblast"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Cardiomyocyte_I"] <- "Cardiomyocyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Cardiomyocyte_II"] <- "Cardiomyocyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Cardiomyocyte_III"] <- "Cardiomyocyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Endocardial"] <- "Endocardial cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Endothelial_I"] <- "Endothelial cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Endothelial_II"] <- "Endothelial cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Endothelial_III"] <- "Endothelial cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Epicardial"] <- "Epicardial cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Lymphatic_endothelial"] <- "LEC"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Mast_cell"] <- "Mast cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Neuronal"] <- "Neuronal cell"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Pericyte_I"] <- "Pericyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Pericyte_II"] <- "Pericyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Proliferating_macrophage"] <- "Macrophage"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Adipocyte"] <- "Adipocyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Lymphocyte"] <- "Lymphocyte"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "Macrophage"] <- "Macrophage"
# seurat.obj$cell_type[seurat.obj$cell_type_leiden0.6 == "VSMC"] <- "VSMC"
# table(seurat.obj$cell_type)
# seurat.obj <- NormalizeData(seurat.obj)
# saveRDS(seurat.obj,file = "~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
seurat.obj <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
##展示线粒体比例的小提琴图----
Idents(seurat.obj) <- "cell_type"
seurat.obj <- PercentageFeatureSet(seurat.obj,
                               pattern = "^MT-",
                               col.name = "percent.mt")
seurat.obj <- PercentageFeatureSet(seurat.obj, 
                               pattern = "^RP[SL]", 
                               col.name = "percent_ribo")
seurat.obj <- PercentageFeatureSet(seurat.obj, 
                               pattern = "^HB[^(P)]", 
                               col.name = "percent_hb")
VlnPlot(seurat.obj,cols = allcolour,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0)

table(seurat.obj$cell_type)

umap = seurat.obj@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = seurat.obj@meta.data$cell_type) # 注释后的label信息 ，改为cell_type
# allcolour=c("#FFA500","#1E90FF","#9370DB","#F08080","#FFFF00",
#             "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
#             "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
#             "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

allcolour <-c( "#9479ad",  # Astrocyte
               "#e77e2c",  # Bcell
               "#66c2a5",  # Endothelial 
               "#CD5C5C",#986156
               "#1271b4",  # Ependymal
               "#80B1D3",  # Granulocyte
               "#FF3030", # Macrophage#7fb687
               "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
               "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","#B09C85FF", 
               "#FFA500") 
p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =1 )  +
  scale_color_manual(values = allcolour)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm')) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) + #设置legend中 点的大小   
  geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                   xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                   xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "umap 1",
           color="black",size = 3, ) + 
  annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "umap 2",
           color="black",size = 3,angle=90)
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
p1 <- p +
  geom_label_repel(aes(label=cell_type), data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none") +
  ggtitle("All genes")+
  theme(plot.title = element_text(size = 25))

p1 <- p +
  theme(legend.position = "none") +
  ggtitle("All genes")+
  theme(plot.title = element_text(size = 25))
p1
ggsave(paste0(output,"no label umap.pdf"), plot = p1, width = 6.96, height = 5.85)

####按照线粒体基因进行降维
library("RColorBrewer")
library(readr)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
gene_mito <- unique(dbs$gene)
# seurat.obj <- readRDS("~/scRNA-heart-mitochodria/data/sample1.rds")
DefaultAssay(seurat.obj) <- "RNA" ##设置为integrated可理解为用整合后的data矩阵做下游分析，非原始值。

seurat.obj<-NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
scale.genes <-  rownames(seurat.obj)
seurat.obj1 <- ScaleData(seurat.obj, features = scale.genes)
seurat.obj1 <- FindVariableFeatures(seurat.obj1)
seurat.obj1 <- RunPCA(object = seurat.obj1, pc.genes = gene_mito)
# 肘部图确定最佳dims
ElbowPlot(seurat.obj1, ndims = 50)
seurat.obj1 <- FindNeighbors(seurat.obj1, 
                             dims = 1:10)
seurat.obj1 <- FindClusters(seurat.obj1, 
                            resolution = 0.6)
seurat.obj1 <- RunUMAP(seurat.obj1, 
                       dims = 1:10)
seurat.obj1 <- RunTSNE(seurat.obj1, 
                       #reduction = "harmony", 
                       dims = 1:10)
umap1 = seurat.obj1@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = seurat.obj1@meta.data$cell_type) # 注释后的label信息 ，改为cell_type
allcolour <-c( "#9479ad",  # Astrocyte
               "#e77e2c",  # Bcell
               "#66c2a5",  # Endothelial 
               "#CD5C5C",#986156
               "#1271b4",  # Ependymal
               "#80B1D3",  # Granulocyte
               "#FF3030", # Macrophage#7fb687
               "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
               "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","#B09C85FF", 
               "#FFA500") 
p <- ggplot(umap1,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =1 )  +
  scale_color_manual(values = allcolour)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm')) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) + #设置legend中 点的大小   
  geom_segment(aes(x = min(umap1$UMAP_1) , y = min(umap1$UMAP_2) ,
                   xend = min(umap1$UMAP_1) +3, yend = min(umap1$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap1$UMAP_1)  , y = min(umap1$UMAP_2)  ,
                   xend = min(umap1$UMAP_1) , yend = min(umap1$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap1$UMAP_1) +1.5, y = min(umap1$UMAP_2) -1, label = "umap_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap1$UMAP_1) -1, y = min(umap1$UMAP_2) + 1.5, label = "umap_2",
           color="black",size = 3, fontface="bold" ,angle=90)
cell_type_med <- umap1 %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p2 <- p +
  theme(legend.position = "none") +
  ggtitle("MitoCarta")+
  theme(plot.title = element_text(size = 25))
p2
ggsave(paste0(output,"nolabel_umap_mito.pdf"), plot = p2, width = 6.96, height = 5.85)
######基因聚类热图-----
####基因聚类全部基因与线粒体基因----
library(tidyverse)
features <- gene_mito
all_expr <-
  AverageExpression(seurat.obj,
                    assays = "RNA",
                    slot = "data",
                    group.by = "cell_type")[["RNA"]]
all_expr <- all_expr[!is.infinite(rowSums(all_expr)),]
mito_expr <-
  AverageExpression(
    seurat.obj,
    assays = "RNA",
    slot = "data",
    group.by = "cell_type",
    features = features
  )[["RNA"]]
all_expr <- as.data.frame(all_expr)
mito_expr <- as.data.frame(mito_expr)
all_expr %>% mutate(across(where(is.character), as.numeric))  -> all_expr
mito_expr %>% mutate(across(where(is.character), as.numeric))  -> mito_expr
bk <- c(seq(0.3, 0.69999, by = 0.01), seq(0.7, 1, by = 0.01))
p3 <- cor(all_expr) %>%
  pheatmap::pheatmap(
    cellwidth = 25,
    cellheight = 25,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 18,
    color=colorRampPalette(c("#1E3163","navy", "#00C1D4", "#FFED99","#FF7600","firebrick3"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(0, 1, 0.2),
    breaks = bk,
    main = "All genes")
ggsave(paste0(output,"Figure1B_1.pdf"), plot = p3, width = 8.93, height = 7.71)
library(RColorBrewer)
bk <- c(seq(0.3, 0.69999, by = 0.01), seq(0.7, 1, by = 0.01))
p4 <- cor(mito_expr) %>%
  pheatmap::pheatmap(
    cellwidth = 25,
    cellheight = 25,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 18,
    color=colorRampPalette(c("#1E3163","navy", "#00C1D4", "#FFED99","#FF7600","firebrick3"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(0, 1, 0.2),
    breaks = bk,
    main = "MitoCarta"
  )
ggsave(paste0(output,"Figure1B_2.pdf"), plot = p4, width = 8.93, height = 7.71)
####基因聚类线粒体基因的疾病组与正常组----
library(tidyverse)
features <- gene_mito
seurat.obj1 <- subset(seurat.obj,subset = group %in% "HCM")
seurat.obj2 <- subset(seurat.obj,subset = group %in% "NF")
mito_HCM <-
  AverageExpression(
    seurat.obj1,
    assays = "RNA",
    slot = "data",
    group.by = "cell_type",
    features = features
  )[["RNA"]]
mito_NF <-
  AverageExpression(
    seurat.obj2,
    assays = "RNA",
    slot = "data",
    group.by = "cell_type",
    features = features
  )[["RNA"]]
mito_HCM <- as.data.frame(mito_HCM)
mito_NF <- as.data.frame(mito_NF)
mito_NF %>% mutate(across(where(is.character), as.numeric))  -> mito_NF
mito_HCM %>% mutate(across(where(is.character), as.numeric))  -> mito_HCM
bk <- c(seq(0.1, 0.59999, by = 0.01), seq(0.6, 1, by = 0.01))
p3 <- cor(mito_HCM) %>%
  pheatmap::pheatmap(
    cellwidth = 25,
    cellheight = 25,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 18,
    color=colorRampPalette(c("#1E3163",
      "#3E5C84",
      "#6D8DB2",
      "white",
     "#BC7C89",
      "#D67474",
      "#8B0000"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(0, 1, 0.2),
    breaks = bk,
    main = "HCM")
p3
ggsave(paste0(output,"线粒体基因HCM聚类热图.pdf"), plot = p3, width = 8.93, height = 7.71)
bk <- c(seq(0.1, 0.59999, by = 0.01), seq(0.6, 1, by = 0.01))
p3 <- cor(mito_NF) %>%
  pheatmap::pheatmap(
    cellwidth = 25,
    cellheight = 25,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 18,
    color=colorRampPalette(c("#1E3163",
                             "#3E5C84",
                             "#6D8DB2",
                             "white",
                             "#BC7C89",
                             "#D67474",
                             "#8B0000"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(0, 1, 0.2),
    breaks = bk,
    main = "NF")
p3
ggsave(paste0(output,"线粒体基因NF聚类热图.pdf"), plot = p3, width = 8.93, height = 7.71)
####全部差异基因的瀑布图-----
library(readxl)
library(catplot)
library(dplyr)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg$change <- ""
deg$change[deg$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg$change[deg$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg$change[deg$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg <- deg[,c(1,3,15,16,17,21)]
colnames(deg) <- c("Gene","cell_type","logFC","Pvalue","Adjusted_Pvalue","change")
deg$log10Pvalue <- log10(deg$Pvalue)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
a <- subset(deg,subset = deg$Gene %in% dbs$gene)
up_deg <- subset(deg,subset = deg$change %in% "up")
down_deg <- subset(deg,subset = deg$change %in% "down")
heatData <- reshape2::dcast(deg[,c("cell_type", "Gene", "change")], 
                            Gene ~ cell_type)
heatData %<>%
  mutate(count = rowSums(is.na(heatData))) %>%
  arrange(count,'Pseudo-Bulk',Cardiomyocyte,Fibroblast,'Endothelial I',Pericyte,Macrophage,VSMC,'Endothelial II',Lymphocyte,Endocardial,Adipocyte,Neuronal,'Lymphatic Endothelial','Mast Cell','Endothelial III')
heatData <- heatData[,-17]
##瀑布图可视化
heatAnnot <- data.frame(number = 
                          case_when(
                            heatData$count < 15 ~ "all",
                            !is.na(heatData$Cardiomyocyte) ~ "Cardiomyocyte",
                            !is.na(heatData$'Pseudo-Bulk') ~ "Pseudo-Bulk",
                            !is.na(heatData$Fibroblast) ~ "Fibroblast",
                            !is.na(heatData$'Endothelial I') ~ "Endothelial I",
                            !is.na(heatData$'Endothelial II') ~ "Endothelial II",
                            !is.na(heatData$'Endothelial III') ~ "Endothelial III",
                            !is.na(heatData$Pericyte) ~ "Pericyte",
                            !is.na(heatData$Macrophage) ~ "Macrophage",
                            !is.na(heatData$VSMC) ~ "VSMC",
                            !is.na(heatData$Lymphocyte) ~ "Lymphocyte",
                            !is.na(heatData$Endocardial) ~ "Endocardial",
                            !is.na(heatData$Adipocyte) ~ "Adipocyte",
                            !is.na(heatData$Neuronal) ~ "Neuronal",
                            !is.na(heatData$'Lymphatic Endothelial') ~ "Lymphatic Endothelial",
                            !is.na(heatData$'Mast Cell') ~ "Mast Cell",
                          ))
####上调----
phData <- 0 + (heatData[rowSums(heatData[2:16] == "up", na.rm = T) > 0 & 
                          rowSums(heatData[2:16] == "down", na.rm = T) == 0, 
                        2:16] == "up")
phData[is.na(phData)] <- 0
dim(phData)
library(pheatmap)
sorted_index <- order(heatAnnot$number)
sorted_heatAnnot <- heatAnnot[sorted_index,,drop = F]
phData <-phData[,c(2,7,13,10,15,9,1,12,4,14,5,6,8,11,3)]
pheatmap(phData,
         cellwidth = 10, 
         cellheight = .05,
         border_color = NA,
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F, 
         color = c("grey97", "lightcoral"),
         annotation_row = heatAnnot, 
         annotation_colors = list(number = setNames(c("#9479ad",  # Astrocyte
                                                      "#e77e2c",  # Bcell
                                                      "#66c2a5",  # Endothelial 
                                                      "#CD5C5C",#986156
                                                      "#1271b4",  # Ependymal
                                                      "#80B1D3",  # Granulocyte
                                                      "#FF3030", # Macrophage#7fb687
                                                      "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
                                                      "#8491B4FF", "#91D1C2FF", "#DC0000FF","#FFA500","grey"),
                                                    c('Pseudo-Bulk',"Cardiomyocyte","Fibroblast",'Endothelial I',"Pericyte","Macrophage","VSMC",'Endothelial II',"Lymphocyte","Endocardial","Adipocyte","Neuronal",'Lymphatic Endothelial','Mast Cell','Endothelial III', "all"))),
         filename = "~/scRNA-heart-mitochodria/figure/final/Figure2/311/scRankHeatmp_up.pdf")
###下调-----
phData <- 0 + (heatData[rowSums(heatData[2:16] == "down", na.rm = T) > 0 & 
                          rowSums(heatData[2:16] == "up", na.rm = T) == 0, 
                        2:16] == "down")
phData[is.na(phData)] <- 0
dim(phData)
library(pheatmap)
phData <-phData[,c(2,7,13,10,15,9,3,1,12,4,14,5,6,8,11)]
pheatmap(phData,
         cellwidth = 10, 
         cellheight = .05, #为了让两个热图一样高，根据自己的基因数量调整
         border_color = NA,
         cluster_rows = F, cluster_cols = F, 
         show_rownames = F, 
         color = c("grey97", "cornflowerblue"),
         annotation_row = heatAnnot, 
         annotation_colors = list(number = setNames(c("#9479ad",  # Astrocyte
                                                      "#e77e2c",  # Bcell
                                                      "#66c2a5",  # Endothelial 
                                                      "#CD5C5C",#986156
                                                      "#1271b4",  # Ependymal
                                                      "#80B1D3",  # Granulocyte
                                                      "#FF3030", # Macrophage#7fb687
                                                      "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
                                                      "#8491B4FF", "#91D1C2FF", "#DC0000FF","#FFA500","grey"),
                                                    c('Pseudo-Bulk',"Cardiomyocyte","Fibroblast",'Endothelial I',"Pericyte","Macrophage","VSMC",'Endothelial II',"Lymphocyte","Endocardial","Adipocyte","Neuronal",'Lymphatic Endothelial','Mast Cell','Endothelial III', "all"))),
         filename = "~/scRNA-heart-mitochodria/figure/final/Figure2/311/scRankHeatmp_down.pdf")
####差异线粒体基因的瀑布图-----
library(readxl)
library(catplot)
library(dplyr)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg$change <- ""
deg$change[deg$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg$change[deg$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg$change[deg$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg <- deg[,c(1,3,15,16,17,21)]
colnames(deg) <- c("Gene","cell_type","logFC","Pvalue","Adjusted_Pvalue","change")
deg$log10Pvalue <- log10(deg$Pvalue)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
a <- subset(deg,subset = deg$Gene %in% dbs$gene)
up_deg <- subset(a,subset = a$change %in% "up")
down_deg <- subset(a,subset = a$change %in% "down")
heatData <- reshape2::dcast(a[,c("cell_type", "Gene", "change")], 
                            Gene ~ cell_type)
heatData %<>%
  mutate(count = rowSums(is.na(heatData))) %>%
  arrange(count,'Pseudo-Bulk',Cardiomyocyte,Fibroblast,'Endothelial I',Pericyte,Macrophage,VSMC,'Endothelial II',Lymphocyte,Endocardial,Adipocyte,Neuronal,'Lymphatic Endothelial','Mast Cell','Endothelial III')
# heatData <- heatData[,-17]
##瀑布图可视化
heatAnnot <- data.frame(number = 
                          case_when(
                            heatData$count < 14 ~ "all",
                            !is.na(heatData$Cardiomyocyte) ~ "Cardiomyocyte",
                            !is.na(heatData$'Pseudo-Bulk') ~ "Pseudo-Bulk",
                            !is.na(heatData$Fibroblast) ~ "Fibroblast",
                            !is.na(heatData$'Endothelial I') ~ "Endothelial I",
                            !is.na(heatData$'Endothelial II') ~ "Endothelial II",
                            !is.na(heatData$'Endothelial III') ~ "Endothelial III",
                            !is.na(heatData$Pericyte) ~ "Pericyte",
                            !is.na(heatData$Macrophage) ~ "Macrophage",
                            !is.na(heatData$VSMC) ~ "VSMC",
                            !is.na(heatData$Lymphocyte) ~ "Lymphocyte",
                            !is.na(heatData$Endocardial) ~ "Endocardial",
                            !is.na(heatData$Adipocyte) ~ "Adipocyte",
                            !is.na(heatData$Neuronal) ~ "Neuronal",
                            !is.na(heatData$'Lymphatic Endothelial') ~ "Lymphatic Endothelial",
                            !is.na(heatData$'Mast Cell') ~ "Mast Cell",
                          ))
####上调----
phData <- 0 + (heatData[rowSums(heatData[2:16] == "up", na.rm = T) > 0 & 
                          rowSums(heatData[2:16] == "down", na.rm = T) == 0, 
                        2:16] == "up")
phData[is.na(phData)] <- 0
dim(phData)
library(pheatmap)
sorted_index <- order(heatAnnot$number)
sorted_heatAnnot <- heatAnnot[sorted_index,,drop = F]
phData <-phData[,c(2,7,13,10,15,1,4,14,5,9,12,6,8,11,3)]
pheatmap(phData,
         cellwidth = 10, 
         cellheight = 1,
         border_color = NA,
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F, 
         color = c("grey97", "lightcoral"),
         annotation_row = heatAnnot, 
         annotation_colors = list(number = setNames(c("#9479ad",  # Astrocyte
                                                      "#e77e2c",  # Bcell
                                                      "#66c2a5",  # Endothelial 
                                                      "#CD5C5C",#986156
                                                      "#1271b4",  # Ependymal
                                                      "#80B1D3",  # Granulocyte
                                                      "#FF3030", # Macrophage#7fb687
                                                      "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
                                                      "#8491B4FF", "#91D1C2FF", "#DC0000FF","#FFA500","grey"),
                                                    c('Pseudo-Bulk',"Cardiomyocyte","Fibroblast",'Endothelial I',"Pericyte","Macrophage","VSMC",'Endothelial II',"Lymphocyte","Endocardial","Adipocyte","Neuronal",'Lymphatic Endothelial','Mast Cell','Endothelial III', "all"))),
         filename = "~/scRNA-heart-mitochodria/figure/final/Figure2/311/mito_scRankHeatmp_up.pdf")
###下调-----
phData <- 0 + (heatData[rowSums(heatData[2:16] == "down", na.rm = T) > 0 & 
                          rowSums(heatData[2:16] == "up", na.rm = T) == 0, 
                        2:16] == "down")
phData[is.na(phData)] <- 0
dim(phData)
library(pheatmap)
phData <-phData[,c(2,7,13,10,1,5,4,15,9,3,12,14,6,8,11)]
pheatmap(phData,
         cellwidth = 10, 
         cellheight = 0.6, #为了让两个热图一样高，根据自己的基因数量调整
         border_color = NA,
         cluster_rows = F, cluster_cols = F, 
         show_rownames = F, 
         color = c("grey97", "cornflowerblue"),
         annotation_row = heatAnnot, 
         annotation_colors = list(number = setNames(c("#9479ad",  # Astrocyte
                                                      "#e77e2c",  # Bcell
                                                      "#66c2a5",  # Endothelial 
                                                      "#CD5C5C",#986156
                                                      "#1271b4",  # Ependymal
                                                      "#80B1D3",  # Granulocyte
                                                      "#FF3030", # Macrophage#7fb687
                                                      "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
                                                      "#8491B4FF", "#91D1C2FF", "#DC0000FF","#FFA500","grey"),
                                                    c('Pseudo-Bulk',"Cardiomyocyte","Fibroblast",'Endothelial I',"Pericyte","Macrophage","VSMC",'Endothelial II',"Lymphocyte","Endocardial","Adipocyte","Neuronal",'Lymphatic Endothelial','Mast Cell','Endothelial III', "all"))),
         filename = "~/scRNA-heart-mitochodria/figure/final/Figure2/311/mito_scRankHeatmp_down.pdf")
#######功能富集气泡图------
library(readxl)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")
deg$change <- ""
deg$change[deg$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg$change[deg$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg$change[deg$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg <- deg[,c(1,3,15,16,17,21)]
colnames(deg) <- c("Gene","cell_type","logFC","Pvalue","Adjusted_Pvalue","change")
deg$log10Pvalue <- log10(deg$Pvalue)
write.csv(deg,file = "~/scRNA-heart-mitochodria/results/deg/degall.csv")
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
a <- subset(deg,subset = deg$Gene %in% dbs$gene)
up_deg <- subset(deg,subset = deg$change %in% "up")
down_deg <- subset(deg,subset = deg$change %in% "down")
#####for循环进行功能富集------
####上调---
cell_types <- unique(deg$cell_type)
for (i in cell_types) {
  #####功能富集-----
  up <- subset(up_deg,up_deg$cell_type %in% i)
  #### GO分析 ----
  library(clusterProfiler)
  library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
  library(dplyr)
  up <-
    up[order(up$logFC, decreasing = T),] %>% pull("Gene")
  
  up_bp <-
    enrichGO(
      up,
      OrgDb = org.Hs.eg.db,
      # 如果是鼠要对应进行更换
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  up_term <- up_bp@result
  up_term$celltype <- i

  write.table(
    up_term,
    file = paste0("~/scRNA-heart-mitochodria/results/function/",i,"_up.txt"),
    row.names = F,
    quote = F,
    sep = "\t"
  )
  saveRDS(up_bp, paste0("~/scRNA-heart-mitochodria/results/function/",i,"_up.rds"))
}
####下调-----
cell_types <- unique(deg$cell_type)
cell_types <- na.omit(cell_types)
for (i in cell_types) {
  #####功能富集-----
 down <- subset(down_deg,down_deg$cell_type %in% i)
  #### GO分析 ----
  library(clusterProfiler)
  library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
  library(dplyr)
 down <-
   down[order(down$logFC, decreasing = T),] %>% pull("Gene")
  
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
 down_term$celltype <- i
  
  write.table(
    down_term,
    file = paste0("~/scRNA-heart-mitochodria/results/function/down/",i,"_down.txt"),
    row.names = F,
    quote = F,
    sep = "\t"
  )
  saveRDS(down_bp, paste0("~/scRNA-heart-mitochodria/results/function/down/",i,"_down.rds"))
}
#####合并功能富集结果
txt_files <- list.files("~/scRNA-heart-mitochodria/results/function", pattern = "\\.txt$", full.names = TRUE)


# 读取并合并所有txt文件
data_list <- lapply(txt_files, function(file) {
  read.delim(file, header = TRUE)  # 假设文件中有标题行，根据实际情况修改
})
# 使用bind_cols()函数合并所有数据框
up_data <- do.call(rbind, data_list)
write_csv(up_data,file = "~/scRNA-heart-mitochodria/results/function/allup.csv")

####下调-------
txt_files <- list.files("~/scRNA-heart-mitochodria/results/function/down", pattern = "\\.txt$", full.names = TRUE)


# 读取并合并所有txt文件
data_list <- lapply(txt_files, function(file) {
  read.delim(file, header = TRUE)  # 假设文件中有标题行，根据实际情况修改
})
# 使用bind_cols()函数合并所有数据框
down_data <- do.call(rbind, data_list)
write_csv(down_data,file = "~/scRNA-heart-mitochodria/results/function/alldown.csv")
down_data <- read.csv("~/scRNA-heart-mitochodria/results/function/alldown.csv")
#####功能富集气泡图-----
#####上调------
library(forcats)
up_data$Description <- as.factor(up_data$Description)
up_data$Description <- fct_inorder(up_data$Description)
up_data1 <- subset(up_data,subset = pvalue < 0.05)
up_data1 <- subset(up_data1,subset = Description %in% c("axonogenesis","cell-matrix adhesion","regulation of cell-matrix adhesion","tissue migration","regulation of cell-substrate adhesion","regulation of actin cytoskeleton organization","response to transforming growth factor beta","regulation of chemotaxis","positive regulation of cytokine production","immune response-activating signaling pathway","positive regulation of Wnt signaling pathway","positive regulation of chemotaxis","interleukin-2 production","activation of immune response","calcium ion transport"))

p99 <- ggplot(up_data1, aes(celltype, Description)) +
  geom_point(aes(color=-log10(pvalue), size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#fe8b8b',high='#890102')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
p99
ggsave(paste0(output,"所有细胞类型功能富集上调气泡图.pdf"), plot = p99, width = 6.75, height = 6.69)
#####下调-----
library(forcats)
down_data$Description <- as.factor(down_data$Description)
down_data$Description <- fct_inorder(down_data$Description)
down_data1 <- subset(down_data,subset = pvalue < 0.05)
down_data1 <- subset(down_data1,subset = Description %in% c("response to oxygen levels","dicarboxylic acid metabolic process","glucose metabolic process","fatty acid metabolic process","response to fatty acid","cellular carbohydrate metabolic process","energy derivation by oxidation of organic compounds","mitochondrial transmembrane transport","fatty acid oxidation","lipid oxidation","pyruvate metabolic process"))

p98 <- ggplot(down_data1, aes(celltype, Description)) +
  geom_point(aes(color=-log10(pvalue), size=Count))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='#58bbfe',high='#014f84')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
p98
ggsave(paste0(output,"所有细胞类型功能富集下调气泡图.pdf"), plot = p98, width = 6.75, height = 6.69)
######得分的featureplot------
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
install.packages("scCustomize")
library(scCustomize)
library(patchwork) 
data.seurat <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
quantile(data.seurat$OXPHOS_UCell)
# 创建一个颜色渐变函数
#对Score进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}#定义矫正函数
data.seurat$OXPHOS_UCell=normalize(data.seurat$OXPHOS_UCell)
quantile(data.seurat$OXPHOS_UCell)

data.seurat1 <- subset(data.seurat,subset = disease %in% "HCM")
data.seurat2 <- subset(data.seurat,subset = disease %in% "NF")
####OXPHOS-Featureplot-----
p1=FeaturePlot(object = data.seurat2,features = c("OXPHOS_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$OXPHOS_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed() +
  theme(legend.position="none") 
p1
p2=FeaturePlot(object = data.seurat1,features = c("OXPHOS_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$OXPHOS_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed()
p2
library(cowplot)
p3 <- plot_grid(p1, p2, ncol = 2,align = "vh")
p3
ggsave(paste0(output,"OXPHOS_feature图.pdf"), plot = p3, width = 8.97, height = 5.13)
####FAO-Featureplot-----
data.seurat$Fatty.acid.oxidation_UCell=normalize(data.seurat$Fatty.acid.oxidation_UCell)
quantile(data.seurat$Fatty.acid.oxidation_UCell)
p1=FeaturePlot(object = data.seurat2,features = c("Fatty.acid.oxidation_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$Fatty.acid.oxidation_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed() +
  theme(legend.position="none") 
p1
p2=FeaturePlot(object = data.seurat1,features = c("Fatty.acid.oxidation_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$Fatty.acid.oxidation_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed()
p2
library(cowplot)
p3 <- plot_grid(p1, p2, ncol = 2,align = "vh")
p3
ggsave(paste0(output,"Fatty.acid.oxidation_feature图.pdf"), plot = p3, width = 8.97, height = 5.13)
#####TCA cycle-featureplot----
data.seurat$TCA.cycle_UCell=normalize(data.seurat$TCA.cycle_UCell)
quantile(data.seurat$TCA.cycle_UCell)
p1=FeaturePlot(object = data.seurat2,features = c("TCA.cycle_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$TCA.cycle_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed() +
  theme(legend.position="none") 
p1
p2=FeaturePlot(object = data.seurat1,features = c("TCA.cycle_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$TCA.cycle_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed()
p2
library(cowplot)
p3 <- plot_grid(p1, p2, ncol = 2,align = "vh")
p3
ggsave(paste0(output,"TCA_feature图.pdf"), plot = p3, width = 8.97, height = 5.13)
####ROS-----
data.seurat$ROS.and.glutathione.metabolism_UCell=normalize(data.seurat$ROS.and.glutathione.metabolism_UCell)
quantile(data.seurat$ROS.and.glutathione.metabolism_UCell)
p1=FeaturePlot(object = data.seurat2,features = c("ROS.and.glutathione.metabolism_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$ROS.and.glutathione.metabolism_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed() +
  theme(legend.position="none") 
p1
p2=FeaturePlot(object = data.seurat1,features = c("ROS.and.glutathione.metabolism_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$ROS.and.glutathione.metabolism_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed()
p2
library(cowplot)
p3 <- plot_grid(p1, p2, ncol = 2,align = "vh")
p3
ggsave(paste0(output,"ROS_feature图.pdf"), plot = p3, width = 8.97, height = 5.13)
#####Calcium uniporter----
data.seurat$Calcium.uniporter_UCell=normalize(data.seurat$Calcium.uniporter_UCell)
quantile(data.seurat$Calcium.uniporter_UCell)
p1=FeaturePlot(object = data.seurat2,features = c("Calcium.uniporter_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$Calcium.uniporter_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed() +
  theme(legend.position="none") 
p1
p2=FeaturePlot(object = data.seurat1,features = c("Calcium.uniporter_UCell"),raster=FALSE)+ scale_color_gradientn(values = quantile(data.seurat$Calcium.uniporter_UCell),colours = c('#779fd3','white','#a13037')) +
  tidydr::theme_dr()+theme(panel.grid = element_blank())+ggplot2::coord_fixed()
p2
library(cowplot)
p3 <- plot_grid(p1, p2, ncol = 2,align = "vh")
p3
ggsave(paste0(output,"Calcium_feature图.pdf"), plot = p3, width = 8.97, height = 5.13)
#####显示差异基因数量的火山图-----
# 计算上调和下调基因的线粒体数量----
library(ggprism)
for (i in cell_types) {
deg$log10AdjPvalue <- log10(deg$Adjusted_Pvalue)
deg1 <- subset(deg,subset = cell_type %in% i)
deg1 <- subset(deg1,subset = Gene %in% dbs$gene)
up_gene <- nrow(subset(deg1, change %in% "up"))
down_gene <- nrow(subset(deg1, change %in% "down"))

p4 <- ggplot(deg1, aes(x =logFC, y=-log10AdjPvalue, colour=change)) +
  geom_point(alpha=0.85, size=1.5) +
  scale_color_manual(values=c('steelblue','gray','brown')) +
  xlim(c(-8, 8)) +
  geom_vline(xintercept= 0,lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.01) , lty=4,col="black",lwd=0.8) + 
  labs(x="logFC", y="-log10adjPvalue") +
  ggtitle(paste0("Mitochondrial related DEGs in ",i)) + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank()) +
  theme_prism(border = T) +
  # # 在图形上添加文本
  geom_text(aes(x = 5, y=-log10(0.05), label = paste("Upregulated: ", up_gene)), color = "black") +
  geom_text(aes(x = -5, y=-log10(0.05), label = paste("Downregulated: ", down_gene)), color = "black")
p4
ggsave(paste0(output,i,"_差异线粒体基因火山图.pdf"), plot = p4, width = 8.59, height = 6.44)
}

# 计算上调和下调基因的全部数量----
library(ggprism)
cell_types <- unique(deg$cell_type)
for (i in cell_types) {
  deg$log10AdjPvalue <- log10(deg$Adjusted_Pvalue)
  deg1 <- subset(deg,subset = cell_type %in% i)
  up_gene <- nrow(subset(deg1, change %in% "up"))
  down_gene <- nrow(subset(deg1, change %in% "down"))
  
  p4 <- ggplot(deg1, aes(x =logFC, y=-log10AdjPvalue, colour=change)) +
    geom_point(alpha=0.85, size=1.5) +
    scale_color_manual(values=c('steelblue','gray','brown')) +
    xlim(c(-8, 8)) +
    geom_vline(xintercept= 0,lty=4,col="black",lwd=0.8)+
    geom_hline(yintercept = -log10(0.01) , lty=4,col="black",lwd=0.8) + 
    labs(x="logFC", y="-log10adjPvalue") +
    ggtitle(paste0("DEGs in ",i)) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank()) +
    theme_prism(border = T) +
    # # 在图形上添加文本
    geom_text(aes(x = 5, y=-log10(0.05), label = paste("Upregulated: ", up_gene)), color = "black") +
    geom_text(aes(x = -5, y=-log10(0.05), label = paste("Downregulated: ", down_gene)), color = "black")
  p4
  ggsave(paste0(output,i,"_差异基因火山图.pdf"), plot = p4, width = 8.59, height = 6.44)
}
######差异基因与线粒体基因的韦恩图------
library(VennDiagram)
library(RColorBrewer)
cell_types <- na.omit(cell_types)
for (i in cell_types) {
deg$log10AdjPvalue <- log10(deg$Adjusted_Pvalue)
deg1 <- subset(deg,subset = cell_type %in% i)
up_gene <- subset(deg1, change %in% "up") %>% pull(Gene)
down_gene <- subset(deg1, change %in% "down") %>% pull(Gene)
  
p5 <- venn.diagram(
  x = list(
    Upregulated = up_gene,
    Downregulated = down_gene,
    Mitocarta = dbs$gene
  ),
  filename = NULL,
  fill = c("#fb8072", "#80b1d3", "#fdb462"),
  cex = 1,
  cat.cex = 1.3, #标签字体大小
  cat.default.pos = "text",  # 标签位置, outer内;text 外
  cat.dist = c(0.1, 0.1, 0.1), # 标签距离圆圈的远近
  cat.pos = c(0,0,0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
  print.mode = c("percent", "raw")
)
grid.draw(p5)
ggsave(paste0(output,i,"_差异基因韦恩图.pdf"), plot = p5, width = 6.53, height = 6.13)}
####只有下调overlap的韦恩图----
cell <- c("Endocardial","Mast Cell","Neuronal")
for (i in cell) {
  deg$log10AdjPvalue <- log10(deg$Adjusted_Pvalue)
  deg1 <- subset(deg,subset = cell_type %in% i)
  down_gene <- subset(deg1, change %in% "down") %>% pull(Gene)
  
  p5 <- venn.diagram(
    x = list(
      Mitocarta = dbs$gene,
      Downregulated = down_gene
    ),
    filename = NULL,
    fill = c( "#fdb462","#80b1d3"),
    cex = 1,
    scaled = F,
    cat.cex = 1.3, #标签字体大小
    cat.default.pos = "text",  # 标签位置, outer内;text 外
    cat.dist = c(0.1, 0.1), # 标签距离圆圈的远近
    cat.pos = c(0,0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
    print.mode = c("percent", "raw")
  )
  grid.draw(p5)
  ggsave(paste0(output,i,"_差异基因韦恩图.pdf"), plot = p5, width = 6.53, height = 6.13)}
#####GSEA的山峦图------
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
seurat.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM_NF311.rds")
seurat.obj <- NormalizeData(seurat.obj)
DefaultAssay(seurat.obj) <- "RNA"
Idents(seurat.obj) <- "group"
markers <- FindMarkers(seurat.obj,
                       ident.1 ="HCM",
                       ident.2 = "NF",
                       min.pct = 0,
                       logfc.threshold = 0)
# 
deg <- markers
write.table(deg,file = "~/scRNA-heart-mitochodria/results/GSEA/HCM_NF_deg.csv")
category <- "C5"
genesets <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = category) %>%
  dplyr::select(gs_name, gene_symbol)

deg <- deg %>% column_to_rownames("gene")
deg <-
  deg[order(deg$log2FoldChange, decreasing = T), ]
genelist <-
  structure(deg$log2FoldChange, names = rownames(deg))
res <- clusterProfiler::GSEA(genelist,
                             TERM2GENE = genesets,
                             seed = 717)

######伪bulk线粒体通路热图------
######线粒体通路打分热图-----
library(readr)
library(readxl)
library(UCell)
# HCM <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
library(UCell)
library(Seurat)
category <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/线粒体通路大类.xlsx")
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
# dbs1 <- split(dbs$gene, dbs$pathway)
# data.HCM <- AddModuleScore_UCell(HCM, features = dbs1,
#                                  ncores = 20)
# dataHCM <- data.HCM@meta.data
# colnames(dataHCM)
seurat.obj <- NormalizeData(seurat.obj)
expr <-
  AverageExpression(seurat.obj,
                    group.by = "biosample_id", assays = "RNA")[["RNA"]]
# gene <- subset(dbs,subset = dbs$pathway %in% "ROS and glutathione metabolism") %>% pull(gene)
# expr1 <- expr[gene,]
expr1 <- expr[c("ACAA2","ACACB","ACADS","ACAT1","HADHB","ETFB","SLC25A20","HADHA","ACSL1",
                "ACO2","CS","DLD","DLST","OGDH",
                "AIFM1","ATP5MC1","COX10","MT-ND6","NDUFA10","NDUFA12",
                "ALDH7A1","MMAA","MMAB","MTHFD2L","MTHFS","PCCA","AMT"),]
library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# 
# data_wide[data_wide >= 0.2] <- 0.2
# breaks = seq(min(unlist(c(data_wide))), max(unlist(c(data_wide))), length.out=100)
breaks = seq(-1, 1, length.out=100)
data <- expr1[,c(1:14,35,36,41:44,47,48,51,52,55:60,15:34,37:40,45,46,49,50,53,54)]
# anno_row <- as.data.frame(rownames(data))
# anno_row$group <- c(rep("Cell type", 13), rep("Disease", 2))
# rownames(anno_row) <- anno_row[,1]

library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# expr1 <- log2(expr1 +1)
breaks = seq(-1, 1, length.out=100)
# expr1[expr1 >= 3] <- 3
data <- t(data)
data <- data[c(31:60,1:30),]
p5 <- pheatmap(data, 
               cellheight=10,cellwidth=40,
               height = 20,width = 20,
               number_color="red", 
               number_format="%.2e",
               border="white",
               fontsize_number = 10, 
               fontsize = 20,
               gaps_col = c(8,13,20),
               gaps_row = c(30),
               show_rownames = F,
               scale = "column",#column
               main="Fatty acid oxidation            OXPHOS                 TCA cycle                        ROS and     \n                                                                                                                      glutathione metabolism",angle_col = 90,
               clustering_distance_rows = "minkowski",
               clustering_method="complete",
               cluster_cols = F,treeheight_col = 20,
               cluster_rows = F,treeheight_row = 20,
               color=colormap,breaks=breaks)
ggsave(paste0(output,"Figure1C.pdf"), plot = p5, width = 25.6, height = 13.38)
#######全部细胞类型通路打分的小提琴图-----
#####Ucell打分----
seurat.obj <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
seurat.obj <- NormalizeData(seurat.obj)
library(UCell)
data.seurat <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
dataHCM <- data.seurat@meta.data
dataHCM$group <- dataHCM$disease
colnames(dataHCM)
dataHCM <- dataHCM[,c(17,28,49:53,59,65,88,101,122,137,146,167)]

library(reshape2)

#####通路带显著性小提琴图-----
df_long <-melt(dataHCM,
               id.vars = c('group','cell_type'),#需要保留不参与聚合的变量,
               value.name='value')
####可视化-----
library(ggpubr)
my_comparisons <- list(c("HCM","NF"))
df_long$group <- factor(df_long$group, levels = c("NF", "HCM"))
# colnames(df_long) <- c("group","pathways","value")

df_long1 <- subset(df_long,subset = variable %in% "Fatty.acid.oxidation_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  xlab("Cell type") + ylab("Fatty acid oxidation signature") +
  geom_hline(yintercept = mean(df_long1$value),  # 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
######分裂小提琴图-----
# library(ggplot2)
# library(ggunchained)
# colours <- c('#08519C',"#E7B800") #设置颜色
# df_long1 <- subset(df_long,subset = variable %in% "Fatty.acid.oxidation_UCell")
# df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
# p_vio <- ggplot(df_long1, aes(x = cell_type, y = value, fill = group)) +
#   geom_split_violin(trim = TRUE, draw_quantiles = NULL) +
#   geom_boxplot(width = .1, alpha = .6, fatten = NULL, colour = NA, outlier.shape = NA, show.legend = FALSE) + 
#   stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = FALSE, position = position_dodge(.175)) +
#   stat_summary(fun.data = "mean_se", geom = "crossbar", width = 0.2, show.legend = FALSE, position = position_dodge(.175)) +  # 将点变成细横线
#   xlab("Cell type") + ylab("Fatty acid oxidation signature") +
#   scale_fill_manual(values = colours, name = "group") +
#   theme_bw() +
#   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#         panel.border = element_rect(size=0.75, colour = "black"),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# print(p_vio)
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
df_long1$value <- as.numeric(df_long1$value)
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p6 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p6
ggsave(paste0(output,"Figure1D_FAO.pdf"), plot = p6, width = 10.91, height = 4.32)
#####OXPHOS小提琴图-----

df_long1 <- subset(df_long,subset = variable %in% "OXPHOS_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cell type")+ylab("OXPHOS signature")+
  geom_hline(yintercept = mean(df_long1$value),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial"))
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p7 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p7
ggsave(paste0(output,"Figure1D_OXPHOS.pdf"), plot = p7, width = 10.91, height = 4.32)
#####TCA cycle的小提琴图-----

df_long1 <- subset(df_long,subset = variable %in% "TCA.cycle_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cell type")+ylab("TCA cycle signature")+
  geom_hline(yintercept = mean(df_long1$value),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial"))
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p7 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p7
ggsave(paste0(output,"Figure1D_TCA.pdf"), plot = p7, width = 10.91, height = 4.32)
#####Calcium.uniporter的小提琴图-----
df_long1 <- subset(df_long,subset = variable %in% "Calcium.uniporter_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cell type")+ylab("Calcium uniporter signature")+
  geom_hline(yintercept = mean(df_long1$value),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial"))
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p7 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p7
ggsave(paste0(output,"Figure1D_Calcium.pdf"), plot = p7, width = 10.91, height = 4.32)
#####Mitophagy_UCell的小提琴图-----
df_long1 <- subset(df_long,subset = variable %in% "Mitophagy_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cell type")+ylab("Mitophagy signature")+
  geom_hline(yintercept = mean(df_long1$value),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial"))
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p7 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p7
ggsave(paste0(output,"Figure1D_Mitophagy.pdf"), plot = p7, width = 10.91, height = 4.32)
#####Lipid.metabolism_UCell_UCell的小提琴图-----
df_long1 <- subset(df_long,subset = variable %in% "Lipid.metabolism_UCell")
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial cell"))
p_vio <- ggviolin(df_long1, 
                  "cell_type",
                  "value",
                  color = "group",
                  add = "boxplot", 
                  palette = c('#08519C',"#E7B800","#FF6347"),
                  add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15),  # 调整字体大小为12
        axis.text.y = element_text(size = 12),  # 调整y轴标签字体大小
        axis.title = element_text(size = 16),   # 调整坐标轴标题字体大小
        legend.text = element_text(size = 12)) +  # 调整图例文字大小
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cell type")+ylab("Lipid metabolism signature")+
  geom_hline(yintercept = mean(df_long1$value),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")
p_vio
##使用 rstatix 包进行组内比较的统计分析
library(rstatix)
library(tidyverse)
library(ggpubr)
#如果期望使用参数检验，以 t 检验为例
df_long1 <- subset(df_long1,subset = !(df_long1$cell_type %in% "Epicardial"))
stat_t <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t1 <- t_test(group_by(df_long1,cell_type),value~group)  #获取组内两两比较的 p 值
stat_t <- stat_t[c(2,1,8,3,5,6,12,4,11,10,7,9),]
# stat_t <- as.data.frame(stat_t)
stat_t <- add_significance(stat_t, 'p') #根据 p 值添加显著性标记 * 符号
#根据 t 检验的结果，在柱形图中添加 p 值或者显著性标记 * 符号..
stat_t.test <-  add_xy_position(stat_t, x = "cell_type", dodge = 0.8)
stat_t.test1 <-  add_xy_position(stat_t1, x = "cell_type", dodge = 0.8)
stat_t.test$x <- stat_t.test1$x
stat_t.test$xmin <- stat_t.test1$xmin
stat_t.test$xmax <- stat_t.test1$xmax
p7 <- p_vio + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.05,
                                 color = 'black',label.size = 5,
                                 hide.ns = T,#是否隐藏NS值
                                 bracket.size = 0.2)  #显示为 * 符号
p7
ggsave(paste0(output,"Figure1D_Lipid.metabolism.pdf"), plot = p7, width = 10.91, height = 4.32)
######功能通路的雷达图-----
meta <- data.HCM@meta.data
colnames(meta)
data_wide <- aggregate(meta, by=list(cell_type=meta$cell_type),mean)##求均值
rownames(data_wide) <- data_wide[,1]
data_wide <- data_wide[,-1]
data_wide <- data_wide[,-1]
data_wide %>% mutate(across(where(is.character), as.numeric))  -> data_wide
library(tidyverse)
library(DESeq2)
library(combinat)
library(clusterProfiler)
library(ggplot2)
library(reshape2)
data_wide$cell_type <- rownames(data_wide)
data_long <- melt(data_wide,
                  id.vars = c('cell_type'),#需要保留不参与聚合的变量,
                  variable.name='Pathways',
                  value.name='Score')
library(tidyr)
ef1 <- subset(data_long,subset = Pathways %in% c("OXPHOS_UCell","Fatty.acid.oxidation_UCell","ROS.and.glutathione.metabolism_UCell","TCA.cycle_UCell"))
ef1<-spread(ef1,cell_type,Score)
rownames(ef1) <- ef1[,1]
ef1 <- ef1[,-1]

#构建雷达图数据
my.data <- matrix( c(rep(max(ef1),ncol(ef1)), 
                     rep(0, ncol(ef1)), 
                     rep(0, ncol(ef1))), nrow = 3, ncol = ncol(ef1), byrow=TRUE)


colnames(my.data) <- colnames(ef1)
rownames(my.data) <- c("max", "min", "score")

my.data <- rbind(my.data, ef1)
my.data <- as.data.frame(my.data)
#把p（显著性阈值）那一列移动到最后一列
my.data <- my.data[c(1,2,5,4,6,7),]
# ACM,             DCM,          HCM,        HF,        MI)
# "#E64B35CC","#4DBBD5CC","#00A087CC","#8491B4CC", "#F39B7FCC"
# my.data$OXPHOS_UCell <- as.numeric(my.data$OXPHOS_UCell)
# my.data <- t(my.data)
# my.data <- as.data.frame(my.data)
library(fmsb)
p8 <- radarchart(my.data, 
                 axistype = 1,
                 pcol = c("#F39B7FCC","#4DBBD5CC","#00A087CC","#E64B35CC","white"), #分组连线颜色
                 pfcol = c(scales::alpha(c("#F39B7FCC","#4DBBD5CC","#00A087CC","#E64B35CC", "white"), c(0.5, 0.5, 0.5, 0.5,0.5))),#每个组填充颜色以及透明度
                 plwd = c(3,3,3,3,3,3),#每个组的连线粗细
                 plty = 1,
                 cglcol = "grey60", #网线颜色
                 cglty = 1, 
                 cglwd = 1,
                 axislabcol = "black",
                 vlcex = 1.2, #外周label文字大小
                 vlabels = colnames(colnames(my.data)),#外周label
                 caxislabels = c(0, 0.06, 0.12, 0.18,0.24),#坐标轴标记
                 calcex=1.3#坐标轴数字大小
)
p8
ggsave(paste0(output,"Figure1E_1.pdf"), plot = p8, width = 9.19, height = 5.96)
# 创建一个空白的绘图画布
plot.new()
legend(x = "bottomright", legend = c("OXPHOS","ROS and glutathione metabolism","Fatty acid oxidation","TCA cycle"), horiz = F,
       bty = "n", pch = 15 , col = c("#E64B35CC","#4DBBD5CC","#00A087CC", "#F39B7FCC"),
       text.col = "black", cex = 1, pt.cex = 1.5)

##########环状热图-------
library(readr)
library(readxl)
library(UCell)
HCM <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311.rds")
library(UCell)
library(Seurat)
library(msigdbr)
dbs2 <- read_excel("~/scRNA-heart-mitochodria/data/seicence转化医学pathways.xlsx")
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
colnames(dbs2) <- c("gs_name", "gene_symbol")
colnames(dbs) <- c("gs_name", "gene_symbol")
deg <- read.table("~/scRNA-heart-mitochodria/results/Macrophage/GSEA/HCM_NF_deg.csv")
genesets <- rbind(dbs,dbs2)
HCM <- NormalizeData(HCM)
expr <-
  AverageExpression(HCM,
                    group.by = c("cell_type","group"), assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
colnames(expr)
expr1 <- expr[,c(1,2,3,4,7,8,10,11,16,17)]

library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
breaks = seq(-1, 1, length.out=100)
library(readxl)
row_an <- read_excel("~/scRNA-heart-mitochodria/data/ISR_legend.xlsx")
row_an <- as.data.frame(row_an)
rownames(row_an) <- row_an[,1]
row_an <- row_an[,-1,drop = FALSE]
col_an <- read_excel("~/scRNA-heart-mitochodria/data/col_an.xlsx")
col_an <- as.data.frame(col_an)
rownames(col_an) <- col_an[,1]
col_an <- col_an[,-1]
Complex <- rownames(row_an)
data <- expr1[Complex,]
data <- na.omit(data)
library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# 
# data_wide[data_wide >= 0.2] <- 0.2
# breaks = seq(min(unlist(c(data_wide))), max(unlist(c(data_wide))), length.out=100)
breaks = seq(-1, 1, length.out=100)

row_an <- subset(genesets,subset = genesets$gs_name %in% c("ATP/NAD+ Carriers","CoQ Synthesis","Cytosolic Protein Import","Type II fatty acid synthesis","mtFASII","TCA cycle"))

row_an <- subset(genesets,subset = genesets$gs_name %in% c("HIF regulators","Selected HIF target genes","HIF1α regulated mitochondrial genes","Hypoxia Inducible Factor family of proteins"))
# rownames(row_an) <- row_an$gene_symbol
col_an <- read_excel("~/scRNA-heart-mitochodria/data/col_an.xlsx")
col_an <- as.data.frame(col_an)
rownames(col_an) <- col_an[,1]
col_an <- col_an[,-1]
gene <- row_an$gene_symbol
data <- expr1[gene,]
data <- na.omit(data)
# anno_row <- as.data.frame(rownames(data))
# anno_row$group <- c(rep("Cell type", 13), rep("Disease", 2))
# rownames(anno_row) <- anno_row[,1]

library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# expr1 <- log2(expr1 +1)
breaks = seq(-1, 0.5, length.out=100)
# expr1[expr1 >= 3] <- 3
# data <- t(data)
ann_colors=list(
  Disease=c(HCM="#3B5998", NF="#F1C40F"),
  Cell_type=c(Cardiomyocyte="#e77e2c",Adipocyte="#9479ad",Endothelial="#CD5C5C",Fibroblast="#80B1D3",Macrophage="#FFC1C1"),
  Unit=c('Assembly Factors'="#F25C54",'Structural Subunits'= "#4A90E2",Others="#F2D5B3")
)
color_map <- c(Cardiomyocyte="#e77e2c",Adipocyte="#9479ad",Endothelial="#CD5C5C",Fibroblast="#80B1D3",Macrophage="#FFC1C1",HCM="#3B5998", NF="#F1C40F")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(dendsort)
library(gridBase)
data <- as.matrix(data)
cir1<-t(scale(t(data)))
head(cir1) #查看数据

### 定义热图颜色梯度：
#设置legend颜色，范围；可从https://www.58pic.com/peisebiao/网站进行配色
mycol1=colorRamp2(c(-2, 0, 2),c("#4575B4","white","#D73027"))
#mycol1 = colorRamp2(c(-2, 0, 2), c("green", "white", "yellow"))
mycol2=colorRamp2(c(-2, 0, 2),c("blue","white","red")) 

### 1 单组环形热图绘制
# 如果矩阵数据分组，可用split参数来指定分类变量
row_an <- as.data.frame(row_an)
row_an = distinct(row_an,gene_symbol,.keep_all = T)
row.names(row_an) = row_an[,2]
row_an <- as.matrix(row_an) # 在circlize函数中，需要为matrix
row_an <- row_an[,-2,drop = FALSE]
#分组绘图
rownames(cir1) <- sub("\\.1$", "", rownames(cir1))
cir1 <- na.omit(cir1)
cir1 <- as.matrix(cir1)
circos.clear()
circos.par(gap.after=c(2,2,2,2,30,2)) # circos.par()调整圆环首尾间的距离，数值越大，距离越宽;让分裂的一个口大一点，可以添加行信息
circos.heatmap(cir1,col=mycol1,
               # dend.side="inside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               rownames.side="inside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.5, # 轨道的高度，数值越大圆环越粗
               rownames.col="black",
               bg.border="black", # 背景边缘颜色
               split = row_an, # 用行注释分裂热图
               show.sector.labels = T,
               rownames.cex=0.5, # 字体大小
               rownames.font=1.2, # 字体粗细
               cluster=FALSE, # cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               # dend.track.height=0.18,#调整行聚类树的高度
               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
                 #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                 # color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               }
)
library(AnnotationDbi)
lg=Legend(title="Expression",col_fun=mycol1,direction = c("horizontal"))
p9 <- grid.draw(lg)
ggsave(paste0(output,"Figure1F_1.pdf"), plot = p9, width = 7.84, height = 5.65)
anno <- t(col_an)

####第二张热图---
row_an <- subset(genesets,subset = genesets$gs_name %in% c("Fatty acid oxidation","Folate and 1-C metabolism","Gluconeogenesis","Glutamate metabolism","Glycolysis","Nucleotide Synthesis"))
gene <- row_an$gene_symbol
data <- expr1[gene,]
data <- na.omit(data)
data <- as.matrix(data)
cir1<-t(scale(t(data)))
head(cir1) #查看数据

### 定义热图颜色梯度：
#设置legend颜色，范围；可从https://www.58pic.com/peisebiao/网站进行配色
mycol1=colorRamp2(c(-2, 0, 2),c("#4575B4","white","#D73027"))
#mycol1 = colorRamp2(c(-2, 0, 2), c("green", "white", "yellow"))
mycol2=colorRamp2(c(-2, 0, 2),c("blue","white","red")) 

### 1 单组环形热图绘制
# 如果矩阵数据分组，可用split参数来指定分类变量
row_an <- as.data.frame(row_an)
row_an = distinct(row_an,gene_symbol,.keep_all = T)
row_an$gs_name[row_an$gs_name == "Fatty acid oxidation"] <- "FAO"
row_an$gs_name[row_an$gs_name == "Folate and 1-C metabolism"] <- "F1C Met"
row_an$gs_name[row_an$gs_name == "Gluconeogenesis"] <- "GNG"
row.names(row_an) = row_an[,2]
row_an <- as.matrix(row_an) # 在circlize函数中，需要为matrix
row_an <- row_an[,-2,drop = FALSE]
row_an$
  #分组绘图
  rownames(cir1) <- sub("\\.1$", "", rownames(cir1))
cir1 <- na.omit(cir1)
cir1 <- as.matrix(cir1)
circos.clear()
circos.par(gap.after=c(2,2,2,2,30,2)) # circos.par()调整圆环首尾间的距离，数值越大，距离越宽;让分裂的一个口大一点，可以添加行信息
cir1 <- na.omit(cir1)
circos.heatmap(cir1,col=mycol1,
               # dend.side="inside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               rownames.side="inside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.4, # 轨道的高度，数值越大圆环越粗
               rownames.col="black",
               bg.border="black", # 背景边缘颜色
               split = row_an, # 用行注释分裂热图
               show.sector.labels = T,
               rownames.cex=0.5, # 字体大小
               rownames.font=10, # 字体粗细
               cluster=FALSE, # cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               # dend.track.height=0.18,#调整行聚类树的高度
               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
                 #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                 # color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               }
)
library(AnnotationDbi)
lg=Legend(title="Expression",col_fun=mycol1,direction = c("horizontal"))
p9 <- grid.draw(lg)
ggsave(paste0(output,"Figure1F_1.pdf"), plot = p9, width = 7.84, height = 5.65)

circos.clear()
rownames(anno) <- c("Group","Cell type")
circos.par(gap.after=c(350))
p10 <- circos.heatmap(anno,col=color_map,
                      dend.side="outside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
                      rownames.side="inside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
                      track.height = 0.5, # 轨道的高度，数值越大圆环越粗
                      rownames.col="black",
                      bg.border="black", # 背景边缘颜色
                      # split = row_an, # 用行注释分裂热图
                      show.sector.labels = T,
                      rownames.cex=0.4, # 字体大小
                      rownames.font=1, # 字体粗细
                      cluster=FALSE, # cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
                      # dend.track.height=0.18,#调整行聚类树的高度
                      dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
                        #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                        # color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
                      }
)
p10
#图例与列名设置
expr <-
  AverageExpression(HCM,
                    group.by = c("cell_type","disease"), assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
# gene <- subset(dbs,subset = dbs$pathway %in% "ROS and glutathione metabolism") %>% pull(gene)
# expr1 <- expr[gene,]
# expr1 <- expr[c("ACAA2","ACACB","ACADS","ACAT1","HADHB","ETFB","SLC25A20","HADHA","ACSL1",
#                 "ACO2","CS","DLD","DLST","OGDH",
#                 "AIFM1","ATP5MC1","COX10","MT-ND6","NDUFA10","NDUFA12","SDHAF1",
#                 "BCKDHB","DBT","ETFA","ETFDH","MCCC2",
#                 "ALDH7A1","MMAA","MMAB","MTHFD2L","MTHFS","PCCA","AMT"),]
colnames(expr)
expr1 <- expr[,c(1,2,3,4,7,8,10,11,16,17)]

library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# 
# data_wide[data_wide >= 0.2] <- 0.2
# breaks = seq(min(unlist(c(data_wide))), max(unlist(c(data_wide))), length.out=100)
breaks = seq(-1, 1, length.out=100)
# ComplexI <- c("NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFA5","NDUFAF3","NDUFAF4","NDUFAF5","NDUFAF6","NDUFAF7","MT-ND1","NDUFA3","NDUFA8","NDUFA9","NDUFA13","NDUFAF2","TIMMDC1","MT-ND2","MT-ND3","MT-ND6","MT-ND4L","NDUFA1","NDUFA10","NDUFC1","NDUFC2", 
# "NDUFS5","ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186","MT-ND4","NDUFB1","NDUFB4","NDUFB5","NDUFB6","NDUFB10","NDUFB11","DMAC2","ATP5SL","FOXRED1","TMEM70","MT-ND5","NDUFAB1","NDUFB2","NDUFB3","NDUFB7","NDUFB8","NDUFB9","DMAC1","TMEM261","NDUFA2","NDUFA6","NDUFA7","NDUFA11","NDUFA12","NDUFS1","NDUFS4","NDUFS6", 
# "NDUFV1","NDUFV2","NDUFV3","NUBPL","AIFM1")
library(readxl)
row_an <- read_excel("~/scRNA-heart-mitochodria/data/ISR_legend.xlsx")
row_an <- as.data.frame(row_an)
rownames(row_an) <- row_an[,1]
row_an <- row_an[,-1,drop = FALSE]
col_an <- read_excel("~/scRNA-heart-mitochodria/data/col_an.xlsx")
col_an <- as.data.frame(col_an)
rownames(col_an) <- col_an[,1]
col_an <- col_an[,-1]
Complex <- rownames(row_an)
data <- expr1[Complex,]
data <- na.omit(data)
# anno_row <- as.data.frame(rownames(data))
# anno_row$group <- c(rep("Cell type", 13), rep("Disease", 2))
# rownames(anno_row) <- anno_row[,1]

library(RColorBrewer)
library(pheatmap)
colormap <- colorRampPalette(c("#4575B4","white","#D73027"))(100)
# expr1 <- log2(expr1 +1)
breaks = seq(-1, 0.5, length.out=100)
# expr1[expr1 >= 3] <- 3
# data <- t(data)
ann_colors=list(
  Group=c(HCM="#3B5998", NF="#F1C40F"),
  'Cell type'=c(Cardiomyocyte="#e77e2c",Adipocyte="#9479ad",Endothelial="#CD5C5C",Fibroblast="#80B1D3",Macrophage="#FFC1C1"),
  Unit=c('Assembly Factors'="#F25C54",'Structural Subunits'="#4A90E2",Others="#F2D5B3"),
  ISRgroup=c('Cytosolic'="white",'Additional ATF4 target genes'="#4A90E2")
)
row_an$ISRgroup <- as.character(row_an$ISRgroup)
colnames(col_an) <- c("Group","Cell type")
pheatmap(data, 
         cellheight=10,cellwidth=10,
         height = 20,width = 20,
         number_color="red", 
         number_format="%.2e",
         border="grey",
         fontsize_number = 10, 
         fontsize = 10,
         # gaps_col = c(2,4,6,8),
         # gaps_row = c(3,7,8,23),
         # annotation_row=row_an,
         annotation_col=col_an,
         annotation_colors=ann_colors,
         show_rownames = T)

######细胞类型线粒体基因聚类网络图------
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
 gene_mito <- unique(dbs$gene)
 features <- gene_mito
 mito_expr <-
     AverageExpression(
         seurat.obj,
         assays = "RNA",
         slot = "data",
         group.by = "cell_type",
         features = features
       )[["RNA"]]
mito_expr <- as.data.frame(mito_expr)
corr <- cor(mito_expr)
###计算p值----
cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(mito_expr) 
head(p.corr[, 1:5])
#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
corpvlue <- corpvlue[corpvlue$pvalue < 0.0001,] #只保留pvalue < 0.0001的
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)
#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)
corpvlue$color <- ifelse(corpvlue$cor > 0, "#FB9A99", "#C6DBEF")
library(igraph)
net <- graph_from_data_frame(d=corpvlue, directed=T) 
V(net)$size <- (1 + corpvlue$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- E(net)$weight/70  #连接之间权重
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=corpvlue$color)


#####整体GSEA的棒棒糖图----
library(readr)
library(tidyverse)
deg <- read.csv("~/scRNA-heart-mitochodria/results/GSEA/HCM_NF_deg.csv",sep = " ")
category <- "C5"
genesets <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = category) %>%
  dplyr::select(gs_name, gene_symbol)

# deg <- deg %>% column_to_rownames("gene")
deg <-
  deg[order(deg$avg_log2FC, decreasing = T), ]
genelist <-
  structure(deg$avg_log2FC, names = rownames(deg))
res <- clusterProfiler::GSEA(genelist,
                             TERM2GENE = genesets,
                             seed = 717)
#### 保存结果 ---
saveRDS(res,
        file = file.path(output,"all_gsea.rds"
        ))
df <- res@result
df1 <- subset(df,subset = df$Description %in% c("GOBP_CELLULAR_RESPIRATION","GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS","GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX","GOBP_ATP_METABOLIC_PROCESS","GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN","GOBP_FATTY_ACID_METABOLIC_PROCESS","GOCC_MITOCHONDRIAL_MATRIX","GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY"))
df1$NES <- as.numeric(df1$NES)
library(ggridges)
cols=c('#779fd3','white','#a13037')
ggplot(df1, aes(x = NES, y = Description)) + 
  geom_density_ridges(
    alpha = 0.6,          # 透明度（0.9 表示 90% 不透明）
    color = 'white'      # 脊线的颜色为白色
    # rel_min_height = 0.01, # 相对最小高度，控制脊线的最小高度，避免过于密集
    # scale = 1,          # 缩放因子，调整脊线的宽度
    # quantile_lines = TRUE, # 是否显示分位数线
    # quantiles = 2,        # 显示的分位数数量（这里是显示中位数线）
    # size = 0.5,             # 修改线条宽度为1（包括中位数线）
    # show.legend = TRUE,   # 是否显示图例
    # bandwidth = 0.2
  ) +
  scale_fill_manual(values = cols) + # 手动设置填充颜色
  scale_color_manual(values = cols)+ # 手动设置线条颜色
  theme_bw() +
  ylab("Gene")+
  xlab("Expression")+
  scale_x_continuous(limits = c(0.5, 5), breaks = seq(0.5, 5, by = 1)) +
  # geom_vline(xintercept = c(0.5, 2.5), size = 0.5, color = 'grey', lty = 'dashed')+
  theme(axis.title.x = element_text(size = 5,color = "black"),  # 修改 x 轴标签字体大小
        axis.title.y = element_text(size = 5,color = "black"),  # 修改 y 轴标签字体大小
        axis.text.x = element_text(size = 5,color = "black"),   # 修改 x 轴刻度字体大小
        axis.text.y = element_text(size = 5,color = "black"))+
  theme(plot.title = element_text(hjust = 0.2, vjust = 2,size = 5)  # 调整标题的水平和垂直位置
  )
ggsave("all2.pdf",height = 3,width = 5)


####topmarker热图-----
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
Idents(seurat.obj) <- "cell_type"
markers <- FindAllMarkers(seurat.obj,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# get top 10 genes
top5markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
mean_gene_exp <- AverageExpression(seurat.obj,
                                   features = top5markers$gene,
                                   group.by = 'cell_type',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# add colnames
colnames(mean_gene_exp) <- gsub("RNA.","",colnames(mean_gene_exp))

# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))

# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
allcolour <-c( "#9479ad",  # Astrocyte
               "#e77e2c",  # Bcell
               "#66c2a5",  # Endothelial 
               "#CD5C5C",#986156
               "#1271b4",  # Ependymal
               "#80B1D3",  # Granulocyte
               "#FF3030", # Macrophage#7fb687
               "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
               "#8491B4FF", "#91D1C2FF") 
allcolour <- allcolour[c(2,1,9,3,6,4,11,7,13,12,10,8,5)]
# top annotation
htdf <- htdf[,c(2,1,9,3,6,4,11,7,13,12,10,8,5)]
names(allcolour) <- colnames(htdf)
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = allcolour))
Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_title = "Top 5 marker genes",
        column_title = "Cell type",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = column_ha,
        # column_split = paste('clsuter ',0:8,sep = ''),
        col = col_fun)
####细胞比例桑吉图-----
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
seurat.obj$cell <- rownames(seurat.obj@meta.data)
Ratio <- seurat.obj@meta.data %>%
  group_by(group, cell_type) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c("#9479ad",  # Astrocyte
            "#e77e2c",  # Bcell
            "#66c2a5",  # Endothelial 
            "#CD5C5C",#986156
            "#1271b4",  # Ependymal
            "#80B1D3",  # Granulocyte
            "#FF3030", # Macrophage#7fb687
            "#F39B7FFF","#FFC1C1", "#4DBBD5FF", "#3C5488FF", 
            "#8491B4FF", "#91D1C2FF")

p2 <- ggplot(Ratio, aes(x =group, y= relative_freq, fill = cell_type,
                        stratum=cell_type, alluvium=cell_type)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Group',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p2
ggsave(paste0(output,"细胞比例桑吉图.pdf"), plot = p2, width = 5.76, height = 3.62)
#####通路雷达图-----
metascore <- readRDS("~/scRNA-heart-mitochodria/data/HCM_NF311_score.rds")
metadata <- metascore@meta.data
data13 <- aggregate(cbind(metadata$OXPHOS_UCell), list(metadata$cell_type,metadata$disease), mean)

colnames(data13) <- c('cell_type',"group",
                      'OXPHOS')
data_wide<-dcast(data13, group~cell_type,
                  value.var = 'OXPHOS')
rownames(data_wide) <- data_wide[,1]
data_wide <- data_wide[,-1]
data_wide <- data_wide[,-5]
library(fmsb)
max(data_wide)
min(data_wide)
# 定义变量最大最小值
max_min <- data.frame(
  Adipocyte = c(0.16, 0), 
  Cardiomyocyte = c(0.16, 0), 
  'Endocardial cell' = c(0.16, 0),
  'Endothelial cell' = c(0.16, 0), 
  Fibroblast = c(0.16, 0), 
  LEC = c(0.16, 0),
  Lymphocyte = c(0.16, 0), 
  Macrophage = c(0.16, 0),
  'Mast cell' = c(0.16, 0),
  'Neuronal cell' = c(0.16, 0),
  Pericyte = c(0.16, 0),
  VSMC = c(0.16, 0)
)
rownames(max_min) <- c("Max", "Min")
colnames(max_min) <- colnames(data_wide)
# 合并数据
df <- rbind(max_min, data_wide)
df
radarchart(
  df, axistype = 1,
  # Customize the polygon
  pcol = c('#AFD796','#61BAAF','#EC74B1','#FFC179','#00C4E4'),
  pfcol = scales::alpha(c('#AFD796','#61BAAF','#EC74B1','#FFC179','#00C4E4'),0.3), plwd = 2, plty = 1,
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "grey", 
  # Variable labels
  vlcex = 0.7, vlabels = colnames(df),
  caxislabels = c(0, 0.04, 0.08, 0.12, 0.16))
# Add an horizontal legend
legend(
  x = "top", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c('#AFD796','#61BAAF','#EC74B1','#FFC179','#00C4E4'),
  text.col = "black", cex = 1, pt.cex = 1.5
)
