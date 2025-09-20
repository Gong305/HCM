setwd("/home/gongfengcz/scRNA-heart-mitochodria")
output <- "/home/gongfengcz/scRNA-heart-mitochodria/figure/final/Figure4/311/"
#######HCM-NF6.13--------
library(readxl)
deg <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/HCMvsNF.xlsx")

deg2 <- subset(deg,subset = deg$`Cell Type` %in% "Fibroblast")
deg2$change <- ""
deg2$change[deg2$`CellBender:\r\r\nlogFC` > 0 ] <- "up"
deg2$change[deg2$`CellBender:\r\r\nlogFC` < 0] <- "down"
deg2$change[deg2$`CellBender:\r\r\nAdjusted P-Value` >= 0.01] <- "not"
deg2 <- deg2[,c(1,15,16,17,21)]
colnames(deg2) <- c("Gene","logFC","Pvalue","Adjusted_Pvalue","change")
deg2$log10Pvalue <- log10(deg2$Pvalue)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
a <- subset(deg2,subset = deg2$Gene %in% dbs$gene)
features <- unique(a$Gene)
library(ggplot2)
library(dplyr)
library(ggrepel)
p1 <- ggplot(deg2, aes(x =logFC, y=-log10Pvalue, colour=change)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5","grey","#ff4757")) + xlim(c(-5, 5)) +  #调整点的颜色和x轴的取值范围 #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = 2.853872, lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="logFC", y="-log10pvalue") +  #x、y轴标签
  ggtitle("HCM/Control-Fibroblast") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank()) 
p1+ geom_label_repel(data = a,
                     aes(x =logFC, y=-log10Pvalue, label = features),
                     size = 2, fill="#CCFFFF",max.overlaps = 300,label.padding = 0.1)
#####功能富集-----
up <- subset(deg2,deg2$change %in% "up")
down <- subset(deg2,deg2$change %in% "down")
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
down_term$celltype <- "Fibroblast"

write.table(
  down_term,
  file = "~/scRNA-heart-mitochodria/results/Fibroblast/function/HCM-NF/down.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(down_bp, "~/scRNA-heart-mitochodria/results/Fibroblast/function/HCM-NF/down.rds")
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
up_term$celltype <- "Fibroblast"


write.table(
  up_term,
  file = "~/scRNA-heart-mitochodria/results/Fibroblast/function/HCM-NF/up.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(up_bp, "~/scRNA-heart-mitochodria/results/Fibroblast/function/HCM-NF/up.rds")

######功能富集柱状图-----
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
subset_down_term <- subset(down_term,subset = down_term$Description %in% c("mitochondrial gene expression","mitochondrial translation","organic acid catabolic process","carboxylic acid catabolic process","cellular respiration","ribose phosphate metabolic process","small molecule catabolic process","energy derivation by oxidation of organic compounds","electron transport chain","fatty acid oxidation"))
subset_down_term$text_x <- rep(0.03,10)
ggplot(data = subset_down_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#EFEFEF", high = "#546de5") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
subset_up_term <- subset(up_term,subset = up_term$Description %in% c("organic acid catabolic process","carboxylic acid catabolic process","amino acid metabolic process","small molecule catabolic process","mitochondrial transmembrane transport","calcium import into the mitochondrion","mitochondrial calcium ion transmembrane transport","mitochondrial calcium ion homeostasis","alpha-amino acid metabolic process"))
subset_up_term$text_x <- rep(0.03,9)
p1 <- ggplot(data = subset_up_term, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = Count)) +
  scale_fill_continuous(low = "#F6edef", high = "#ff4757") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p1
ggsave(paste0(output,"GO富集柱状图上调.pdf"), plot = p1, width = 4.45, height = 4.05)


########HCM-NF------
#######轨迹分析------
####cytotrace------
#### 导入包
library(Seurat)
library(dplyr)
library(ggplot2)
library(CytoTRACE)
set.seed(717)
#### 导入数据
#### 创建CDS对象并预处理数据
metacell.obj <- qs::qread("/home/gongfengcz/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
Fib <- subset(metacell.obj,subset = cell_type %in% "Fibroblast")
Fib <- NormalizeData(Fib)
Idents(Fib) <- "sub_cell_type"
DimPlot(Fib,
        label = T,        
        repel = T,
        reduction = "umap")

table(Idents(Fib))
table(is.na(Fib$group))

#### 提取表型文件 
table(Fib$sub_cell_type)
phe <- Fib$sub_cell_type
phe = as.character(phe)
names(phe) <- rownames(Fib@meta.data)
#### 提取表达矩阵
mat_3k <- as.matrix(Fib@assays$RNA@counts)
mat_3k[1:4,1:4]
#### 当数据集中的单元数超过3,000时，CytoTRACE将自动以快速模式运行，这是一种用于减少运行时和内存使用的子采样方法。
#此外用户还可以使用ncores(默认值为1)来多线程，或者使用subsamplingsize(默认值为1000 cells)来指示子采样大小。
#在快速模式下运行以下数据集，使用8个核心，子样本大小为1000。
results <- CytoTRACE(mat = mat_3k)       #可以问一下?CytoTRACE 看看help
plotCytoGenes(results, numOfGenes = 100)  #基因的数量可以改变
plotCytoTRACE(results, phenotype = phe, outputDir = "HCM_NF") #默认用tsne降维，可以输入自己的嵌入
#####monocle2-----
#######轨迹monocle2------
library(Seurat)
library(ggplot2)
library(patchwork)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
trace('project2MST', edit = T, where = asNamespace("monocle"))#版本问题加上这行代码，代开源代码，
# #去掉35行
# seurat_obj <-
#   qs::qread(
#     "./data/public/2021_human-space-time-gut-cell-atlas_RasaElmentaite/data/adult_colon.qs"
#   )
# stem_cells <- seurat_obj |> 
#   subset(sub_cell_type == "Stem cells")
# stem_cells$sub_cell_type <- as.character(stem_cells$sub_cell_type)

DimPlot(Fib)
Fib$sub_cell_type <- as.character(Idents(Fib))
table(Fib$sub_cell_type)
Fib$sub_cell_type
#colonocyte <- merge(Fib, stem_cells)
# sort(table(colonocyte$sub_cell_type))
colonocyte <- Fib

# dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
# dbs1 <- split(dbs$gene, dbs$pathway)
# Fib <- AddModuleScore_UCell(Fib, features = dbs1,
#                                        ncores = 20)
# Idents(Fib) <- "group"
# DimPlot(Fib)
# Fib$sub_cell_type <- as.character(Idents(Fib))
# table(Fib$sub_cell_type)
# Fib$sub_cell_type
# # colonocyte <- merge(Fib, stem_cells)
# # sort(table(colonocyte$sub_cell_type))
# colonocyte <- Fib
outdir <- "/home/gongfengcz/scRNA-heart-mitochodria/results/Fibroblast/monocle2/1e-10/"
qvalues <- 1e-10
cds <-
  cat_monocle2(colonocyte,
               group_by = "sub_cell_type",
               qvalue = qvalue,
               outdir = outdir)
cds <- qs::qread("~/scRNA-heart-mitochodria/results/Fibroblast/monocle2/1e-10/cds.qs")
cds <- orderCells(cds)
p1 <-
  plot_cell_trajectory(cds, color_by = "sub_cell_type", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p2 <-
  plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p <- p1 + p2
p
# 以细胞类型上色
plot_cell_trajectory(cds,
                     color_by = "disease",
                     size = 1,
                     show_backbone = TRUE)
# 以State上色
plot_cell_trajectory(cds,
                     color_by = "State",
                     size = 1,
                     show_backbone = TRUE)
#####对通路进行打分-----
plot_cell_trajectory(cds,
                     color_by = "OXPHOS_UCell",
                     size = 1,
                     show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 

#####轨迹山峦图------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=Biobase::pData(cds)

ggplot(plotdf, aes(x=Pseudotime,y=disease,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.6,9.7),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
ggsave("山峦图1.pdf",width = 13,height = 7,units = "cm")
ggplot(plotdf, aes(x=Pseudotime,y=sub_cell_type,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.7,9.7,15.5),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
ggsave("山峦图2.pdf",width = 13,height = 7,units = "cm")
###线粒体通路打分的山峦图-----
ggplot(plotdf, aes(x=Pseudotime,y='Fatty.acid.oxidation_UCell',fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.7,9.7,15.5),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
ggsave("山峦图4.pdf",width = 13,height = 7,units = "cm")

# library(ClusterGVis)
# Fib$state <- cds@phenoData@data[["State"]]
# Fib$State <- ""
# Fib$State[Fib$state == "1"] <- "state 1"
# Fib$State[Fib$state == "2"] <- "state 2"
# Fib$State[Fib$state == "3"] <- "state 3"
# Fib$State[Fib$state == "4"] <- "state 4"
# Fib$State[Fib$state == "5"] <- "state 5"
# Fib$Pseudotime <- cds@phenoData@data[["Pseudotime"]]
# saveRDS(Fib,"/home/gongfengcz/scRNA-heart-mitochodria/results/Fibroblast/monocle2/1e-04/sub_seurat_obj.rds")
# library(Seurat)
# expr <-
#   AverageExpression(Fib,
#                     group.by = "State", assays = "RNA")[["RNA"]]
# 
# df2 <- expr[Time_genes,]
# library(ClusterGVis)
# df2 <- df2[!is.infinite(rowSums(df2)),]
# df2 <- na.omit(df2)
# ck <- clusterData(exp = df2,
#                   cluster.method = "mfuzz",
#                   cluster.num = 3)
# markGenes = unique(dbs$gene)
# visCluster(object = ck,
#            plot.type = "heatmap",
#            column_names_rot = 45,
#            markGenes = markGenes)
# ####对每个cluster进行功能富集----
# clusters <- cutree(p$ph$tree_row, k = 1)
# table(clusters)
# a <- ck$wide.res
# library(org.Hs.eg.db)
# output <- vector("double", ncol(a)) 
# for (i in 1:5) {
#   a1 <- subset(a,subset = a$cluster %in% i)
#   genes <- a1$gene
#   bp <-
#     enrichGO(
#       genes,
#       OrgDb = org.Hs.eg.db,
#       # 如果是鼠要对应进行更换
#       keyType = 'SYMBOL',
#       ont = "BP",
#       pAdjustMethod = "BH",
#       pvalueCutoff = 0.05,
#       qvalueCutoff = 0.2
#     )
#   term <- bp@result
#   term$cluster <- i
#   output <- rbind(output,term)
# }
# output1 <- subset(output,subset = output$Description %in% c("positive regulation of cytokine production","MHC class II protein complex assembly","calcium ion homeostasis","response to ATP"###cluster1
#                                                             ,"mononuclear cell differentiation","cell activation involved in immune response"###c2
#                                                             ,"regulation of mononuclear cell proliferation","apoptotic mitochondrial changes","regulation of leukocyte proliferation"###C3
#                                                             ,"immune response-activating signal transduction","regulation of inflammatory response","positive regulation of cell activation"###C4
#                                                             ,"antigen processing and presentation of peptide antigen","myeloid leukocyte activation"))
# output1 <- output1[,c(2,5,10)]
# output1$cluster <- paste0("C",output1$cluster)
# termanno2 <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/巨噬细胞轨迹.xlsx")
# visCluster(object = ck,
#            plot.type = "both",
#            column_names_rot = 45,
#            markGenes = markGenes,
#            annoTerm.data = termanno2)
# visCluster(object = ck,
#            plot.type = "both",
#            column_names_rot = 45,
#            markGenes = markGenes,
#            markGenes.side = "left",
#            genes.gp = c('italic',fontsize = 12,col = "black"),
#            annoTerm.data = termanno2,
#            line.side = "left")


#####得分在轨迹上的变化----
library(UCell)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
seurat_obj<- AddModuleScore_UCell(seurat_obj, features = dbs1,
                                  ncores = 20)
a <- seurat_obj@meta.data
b <- plotdf[,c(19:23)]
a <- cbind(a,b)
ggplot(a, aes(x=OXPHOS_UCell,y=group,fill=group))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(0.465,0.467),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
#######HCM-NF 7.2-----
metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM-NF.rds")
Fib <- subset(metacell.obj,subset = cell_type %in% "Fibroblast")
Fib <- NormalizeData(Fib)
Idents(Fib) <- "sub_cell_type"
Fib2 <- subset(Fib,subset = disease %in% "HCM")
library(RColorBrewer)
cell_type_cols <- c('#7db954',
                    '#efb306',
                    '#852f88',"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999", 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00") 
p2 <- DimPlot(Fib2, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Fibroblast_NF") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2
ggsave(paste0(output,"umap成纤维HCM.pdf"), plot = p2, width = 4.99, height = 4.32)
Fib1 <- subset(Fib,subset = disease %in% "NF")
library(RColorBrewer)
cell_type_cols <- c('#7db954',
                    '#efb306',
                    '#852f88',"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999", 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00") 
p3 <- DimPlot(Fib1, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Fibroblast_NF") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p3
ggsave(paste0(output,"umap成纤维nf.pdf"), plot = p3, width = 7.01, height = 4.32)
#######细胞比例的桑吉图------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
Ratio <- Fib@meta.data %>%
  group_by(disease, sub_cell_type) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c(
  '#7db954',
  '#efb306',
  '#852f88',
  '#4e54ac',
  '#0f8096',
  'pink',
  'green')

p4 <- ggplot(Ratio, aes(x =disease, y= relative_freq, fill = sub_cell_type,
                        stratum=sub_cell_type, alluvium=sub_cell_type)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Disease',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p4
ggsave(paste0(output,"成纤维比例桑吉图.pdf"), plot = p4, width = 5.09, height = 2.84)
######UCell对线粒体通路在HCM以及NF组中进行打分-----
########UCELL打分----柱状图可视化-----
library(UCell)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
table(Fib$sub_cell_type)
data.seurat <- AddModuleScore_UCell(Fib, features = dbs1,
                                    ncores = 50)
meta<- data.seurat@meta.data
meta <- meta[, c("disease",
                 colnames(meta)[grep("_UCell", colnames(meta))])]   #提取出含GOBP的
meta <-
  aggregate(meta[,-1], list(meta$disease), FUN = mean)
rownames(meta) <- meta[, 1]
meta <- meta[, -1]
meta <- t(meta)

#### 差异分析 ----
library(limma)
meta<- data.seurat@meta.data
es.matrix <- meta[,c(colnames(meta)[grep("_UCell", colnames(meta))])]

meta <-
  data.seurat@meta.data[, c("disease", "cell_type")]
meta$cell_id <- rownames(meta)
es.matrix <- t(es.matrix)
es.matrix <- es.matrix[,meta$cell_id]
clusters <- names(table(meta$cell_type))

case <- "HCM"
control <- "NF"
cluster <- "Fibroblast"
t_results <- map_dfr(clusters, function(cluster) {
  cell_id.1 <- meta[meta$cell_type == cluster &
                      meta$disease == case,]$cell_id
  cell_id.2 <- meta[meta$cell_type == cluster &
                      meta$disease == control,]$cell_id
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

library(ggplot2)
p5 <-
  ggplot(t_results, aes(x = reorder(feature,t), y = t, fill = fill)) +
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
  xlab("Fibroblast  (HCM/NF)") +
  ylab("t value of UCell score") +
  geom_text(
    data = subset(t_results, t < 0),
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
    data = subset(t_results, t > 0),
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
      color = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black", size = 0.5),
  )
p5

t_results1 <- subset(t_results,subset = t_results$feature %in% c(rownames(t_results)[grep("Calcium", rownames(t_results))],
                                                                 rownames(t_results)[grep("OXPHOS", rownames(t_results))],rownames(t_results)[grep("Fatty", rownames(t_results))],rownames(t_results)[grep("Complex", rownames(t_results))]))
p5 <-
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
  xlab("Fibroblast  (HCM/NF)") +
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
      color = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black", size = 0.5),
  )
p5
ggsave(paste0(output,"线粒体功能打分差异柱状图.pdf"), plot = p5, width = 6.47, height = 3.49)
#######UCell在激活成纤维细胞以及未激活成纤维细胞之间的功能差异------
####Ucell打分----
library(AUCell)
library(clusterProfiler)
library(tidyverse)
library(msigdbr)
HCM <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM-NF.rds")
HCM <- NormalizeData(HCM)
Fib <- subset(HCM,subset = cell_type %in% "Fibroblast")
Idents(Fib) <- Fib$sub_cell_type
cells_rankings <- AUCell_buildRankings(Fib@assays$RNA@data)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)#
#Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name <- gsub('GOBP_', '', Dataset$gs_name)
length(unique(Dataset$gs_name)) #50
#View(Dataset)
#挑选通路
# feature_list <- vascular_GS$Geneset
# length(unique(feature_list)) #614
dataset <- Dataset #%>% filter(Dataset$gs_name %in% feature_list)
length(unique(dataset$gs_name)) #526
dbs1 <- split(dataset$gene_symbol, dataset$gs_name)
library(UCell)
data.seurat <- AddModuleScore_UCell(Fib, features = dbs1,maxRank = 15000,
                                    ncores = 100)
meta<- data.seurat@meta.data
meta <- meta[, c("sub_cell_type",
                 colnames(meta)[grep("_UCell", colnames(meta))])]
t_results <- meta
# 保存t值结果
saveRDS(t_results,
        file =
          "~/scRNA-heart-mitochodria/results/Fibroblast/UCell/C5_ALL_UCscore.rds")
t_results <- readRDS("~/scRNA-heart-mitochodria/results/Fibroblast/UCell/C5_ALL_UCscore.rds")
#### 可视化 ----
rownames(t_results) <- gsub("_UCell", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
t_results$cell_id <- rownames(t_results)
t_results1 <- t_results[,c(2:15475)]
# t_results1 <- t_results[,c(colnames(t_results)[grep("EXTRACELLULAR_MATRIX", colnames(t_results))])]
t_results1 <- t_results[,c(colnames(t_results)[grep("FIBROSIS", colnames(t_results))],
                           colnames(t_results)[grep("RANSFORMING_GROWTH_FACTOR", colnames(t_results))],colnames(t_results)[grep("CELL_PROLIFERATION", colnames(t_results))],colnames(t_results)[grep("CELL_MIGRATION", colnames(t_results))],colnames(t_results)[grep("EXTRACELLULAR_MATRIX", colnames(t_results))])]
# t_results2 <- t_results1[,c(colnames(t_results1)[grep("HP", colnames(t_results1))])]
es.matrix <- t(t_results1)
#### 差异分析 ----
case <- "Activated_fibroblast"
control1 <- "Fibroblast_I"
control2 <- "Fibroblast_II"
meta <-
  Fib@meta.data[, c("sub_cell_type",
                    "disease")]
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
# rownames(t_results)[grep("CYTOKINES", rownames(t_results))],rownames(t_results)[grep("FIBRO", rownames(t_results))]))

# t_results1 <- subset(t_results,subset = t_results$feature %in% c("REGULATION_OF_MONOCYTE_EXTRAVASATION_UCell","COMPLEMENT_ACTIVATION_ALTERNATIVE_PATHWAY_UCell","CHRONIC_INFLAMMATORY_RESPONSE_UCell","REGULATION_OF_FIBRINOLYSIS_UCell"))
p <-
  ggplot(t_results, aes(x = reorder(feature,t), y = t, fill = fill)) +
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
    color = "white",
    linetype = "dashed",
    size = 0.5
  ) +
  coord_flip() +
  xlab("Activated_fibroblast/Fibroblast_I") +
  ylab("t value of score") +
  geom_text(
    data = subset(t_results, t < 0),
    aes(
      x = feature,
      y = 0.1,
      label = paste0(" ", feature),
      color = color
    ),
    size = 2.5,
    hjust = "inward"
  ) +
  geom_text(
    data = subset(t_results, t > 0),
    aes(
      x = feature,
      y = -0.1,
      label = paste0(" ", feature),
      color = color
    ),
    size = 2.5,
    hjust = "outward"
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
p
ggsave(
  filename = "~/scRNA-heart-mitochodria/results/Fibroblast/UCell/activated/Fib_barplot.png",
  plot = p,
  height = 10,
  width = 15
)

####筛选通路绘制多颜色barplot-----
t_results1 <- subset(t_results,subset = t_results$feature %in% c("HP_MYOCARDIAL_FIBROSIS_UCell"))
t_results1$cluster <- "Fibrosis"
t_results2 <- subset(t_results,subset = t_results$feature %in% c("TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_UCell","RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_UCell","TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION_UCell","TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION_UCell"))
t_results2$cluster <- "TGF beta"
t_results3 <- subset(t_results,subset = t_results$feature %in% c("CELL_MIGRATION_UCell","CELL_MIGRATION_INVOLVED_IN_HEART_DEVELOPMENT_UCell"))
t_results3$cluster <- "Cell Migration"
t_results4 <- subset(t_results,subset = t_results$feature %in% c("EXTRACELLULAR_MATRIX_ASSEMBLY_UCell","REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY_UCell","POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY_UCell"))
t_results4$cluster <- "Extracellular Matrix"
t_results5 <- rbind(t_results1,t_results2,t_results3,t_results4)
library(ggpubr)
p6 <- ggbarplot(t_results5, x="feature", y="t", fill = "cluster", color = "white", 
                width = 0.5,
                orientation = "horiz",  #横向显示
                palette = "aaas",    #配色方案
                legend = "right",    #图例位置
                sort.val = "asc",    #上升排序，区别于desc
                sort.by.groups=TRUE,
                lab.size = 6,
)+    #按组排序
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))# 
p6
ggsave(paste0(output,"亚群功能富集条形图.pdf"), plot = p6, width = 13.15, height = 5.01)
######成纤维亚群中线粒体功能的横向正方体变化------
library(org.Hs.eg.db)
library(GOSemSim)
####UCell打分---
library(UCell)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
Fib <- AddModuleScore_UCell(Fib, features = dbs1,
                            ncores = 50)
meta <- Fib@meta.data
colnames(meta)
meta1 <- meta[,c(3,7,122,28,26)]
# meta2 <- aggregate(meta1, by=list(cell_type=meta1$sub_cell_type),mean)
f <- function(y) {
  r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  r[3] <- mean(y)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
p7 <- ggplot(meta1, aes(sub_cell_type, Calcium.cycle_UCell, fill = sub_cell_type)) + 
  scale_fill_brewer(palette="Set3") + #配色
  guides(fill=FALSE) + #不显示图例
  
  stat_summary(fun.data= f, geom='boxplot') + 
  geom_hline(aes(yintercept=0.2), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细

p7
ggsave(paste0(output,"附图成纤维亚群的钙离子功能柱状图.pdf"), plot = p7, width = 6.24, height = 4.01)
######细胞轨迹------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
cds <- qs::qread("~/scRNA-heart-mitochodria/results/Fibroblast/monocle2/循环轨迹/1e-10/cds.qs")
# cds <- orderCells(cds)
p8 <-
  plot_cell_trajectory(cds, color_by = "sub_cell_type", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p8
ggsave(paste0(output,"轨迹亚群图.pdf"), plot = p8, width = 6.70, height = 4.45)
p9 <-
  plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5) +
  theme(aspect.ratio = 1,
        legend.title = element_blank())
p9
ggsave(paste0(output,"轨迹拟时序图.pdf"), plot = p9, width = 6.70, height = 4.45)
# 以细胞类型上色
plot_cell_trajectory(cds,
                     color_by = "disease",
                     size = 1,
                     show_backbone = TRUE)# 以State上色
plot_cell_trajectory(cds,
                     color_by = "State",
                     size = 1,
                     show_backbone = TRUE)


###线粒体通路打分阴阳表达组的山峦图-----
library(UCell)
library(readxl)
library(readr)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
dbs1 <- split(dbs$gene, dbs$pathway)
sub_metacell.obj <- AddModuleScore_UCell(sub_metacell.obj, features = dbs1,
                                       ncores = 50)
plotdf=Biobase::pData(cds)
plotdf <- plotdf[,c(19,20,21)]
a <- sub_metacell.obj@meta.data
plotdf1 <- cbind(a,plotdf)
x <- "Calcium.cycle_UCell"
median(plotdf1[,"Calcium.cycle_UCell"])
sub_metacell.obj[[str_c(x, "_group")]] <-
  if_else(plotdf1[,x] > 0,
          "high", "low")
a <- sub_metacell.obj@meta.data
plotdf <- cbind(a,plotdf)
# sub_metacell.obj@meta.data <- plotdf

Biobase::pData(cds) <- plotdf
#####对通路进行打分-----
plot_cell_trajectory(cds,
                     color_by = "OXPHOS_UCell",
                     size = 1,
                     show_backbone = TRUE) + viridis::scale_color_viridis(option = "D") 
#####轨迹山峦图------
p10 <- ggplot(plotdf, aes(x=Pseudotime,y=disease,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.55,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
p10
ggsave(paste0(output,"疾病山峦图.pdf"), plot = p10, width = 5.05, height = 3.40)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
p11 <- ggplot(plotdf, aes(x=Pseudotime,y=Calcium.uniporter_UCell_group,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.7,13.2),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
p11
ggsave(paste0(outdir,"山峦图2.pdf"),width = 13,height = 7,units = "cm")
#####合并到一起的山峦图---
library(ggpubr)
library(ggsci)
library(tidyverse)

b <- ggplot(plotdf, aes(Pseudotime, colour = sub_cell_type, fill = sub_cell_type)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) +
  theme_bw() +
  scale_colour_npg() +
  scale_fill_npg() +
  theme(panel.background = element_blank())

b
#### 分组热图 ----
# 设置顺序 ----
ordergene <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
ordergene <- ordergene$gene
ordergene <- unique(ordergene)
ordergene <- subset(ordergene,subset = ordergene %in% rownames(cds)) 
#ordergene2 <- c('RFX5',ordergene)
Time_diff <- differentialGeneTest(cds[ordergene,],
                                  cores = 50,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
outdir <- "~/scRNA-heart-mitochodria/results/Fibroblast/monocle2/循环轨迹/1e-10/"
write.csv(Time_diff,file.path(outdir, "Time_diff_DEG2.csv"), row.names = T)
#Time_diff <- read.csv(file.path(outdir, "Time_diff_DEG2.csv"))

# 7.4绘图 ----
pdf(file=file.path(outdir, "non-branched_heatmap_all.pdf"), width= 15, height= 15)
plot_pseudotime_heatmap(
  cds[row.names(subset(Time_diff, qval < 0.01)),],
  cluster_rows = T,
  num_clusters = 5, cores = 20,
  use_gene_short_name = T,
  show_rownames = F
)
dev.off()

library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
library(dplyr)
library(tidytree)
library(viridis)
library(scales)

interested_genes <-  dbs %>% 
  filter(pathway %in% c("Calcium uniporter","Calcium cycle")) %>% pull(gene)
Time_diff <- Time_diff %>% filter(gene_short_name %in% interested_genes)
pdf(file=file.path(outdir, "non-branched_heatmap_all.pdf"), width= 15, height= 15)
p <- plot_pseudotime_heatmap(cds[row.names(subset(Time_diff, qval < 0.01)),], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
dev.off()
#### 分组 ----

#'Normal'
case <- 'HCM'
cds1 <- cds[, cds[["disease"]] == case]
# Time_diff <- differentialGeneTest(cds[ordergene,],
#                                   cores = 1,
#                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)

Time_diff <- differentialGeneTest(cds1[ordergene,],
                                  cores = 50,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
# interested_genes <-  dbs %>% 
#   filter(pathway %in% c("Calcium uniporter","Calcium cycle","Calcium homeostasis")) %>% pull(gene)
# Time_diff <- Time_diff %>% filter(gene_short_name %in% interested_genes)
pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
p <- plot_pseudotime_heatmap(cds1[row.names(subset(Time_diff, qval < 0.05)),], 
                             num_cluster = 3, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
dev.off()
gene <- read_excel("~/scRNA-heart-mitochodria/results/Fibdiomyocyte/monocle2/轨迹gene.xlsx") %>% pull(gene)
gene <- subset(gene,subset = gene %in% rownames(cds3))
#'HCM'
case <- 'NF'
cds2 <- cds[, cds[["disease"]] == case]
# Time_diff <- differentialGeneTest(cds[ordergene,],
#                                   cores = 1,
#                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)
pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
Time_diff <- differentialGeneTest(cds2[ordergene,],
                                  cores = 10,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
p <- plot_pseudotime_heatmap(cds2[row.names(subset(Time_diff, qval < 0.05)),], 
                             num_cluster = 3, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
dev.off()

case <- 'HCM'
cds3 <- cds[, cds[["disease"]] == case]
# Time_diff <- differentialGeneTest(cds3[ordergene,],
#                                   cores = 10,
#                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(Time_diff,file.path(outdir, paste0(case, "_Time_diff.csv")), row.names = T)
pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
Time_diff3 <- differentialGeneTest(cds3[ordergene,],
                                   cores = 50,
                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
interested_genes <-  dbs %>% 
  filter(pathway %in% c("Calcium uniporter","Calcium cycle")) %>% pull(gene)
Time_diff <- Time_diff %>% filter(gene_short_name %in% interested_genes)

p <- plot_pseudotime_heatmap(cds3[row.names(subset(Time_diff3, qval < 0.1)),], 
                             num_cluster = 3, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
dev.off()

genes <- rownames(cds2[row.names(subset(Time_diff, qval < 0.05)),])
genes1 <- subset(genes,subset = genes %in% row.names(subset(Time_diff3, qval < 0.5)))
p <- plot_pseudotime_heatmap(
  cds3[genes1,],
  cluster_rows = F,
  show_rownames = T,
  cores = 20,
  return_heatmap = T,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100)
)
p
#ggsave(file.path("~/Retina_cell/results/monocle2/figure1.png"),plot = p2,height = 2,width = 2)

# 提取热图的基因。
clusters <- cutree(p$tree_row, k = 5)
table(clusters)
cluster_genes <- names(clusters[clusters == 1])
dbs1 <- subset(dbs,subset = dbs$gene %in% genes)
write.table(genes,
            quote = F,
            row.names = F,
            col.names = F)
####轨迹整体热图的绘制----
pdf(file=file.path(outdir, paste0(case, "_non-branched_heatmap_lable.pdf")), width= 5, height= 15)
Time_diff <- differentialGeneTest(cds2[ordergene,],
                                  cores = 10,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
p <- plot_pseudotime_heatmap(cds2[row.names(subset(Time_diff, qval < 0.05)),], 
                             num_cluster = 3, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p

p <- plot_genes_in_pseudotime(cds[cluster_genes,], color_by = "group")+
  scale_color_manual(values = c("#ee6470","#00a6e1"))
df <- p$data
ggplot(df, aes(Pseudotime, adjusted_expression)) + theme_classic()+
  geom_smooth(aes(colour = disease),method = "loess", se = FALSE)+
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank())+
  scale_color_manual(values = c("#ee6470","#00a6e1"))

######相关性计算-----
#####相关性散点图-去掉0值-----
########线粒体相关功能与炎症功能
library(UCell)
library(readxl)
library(readr)
# metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM-NF.rds")
# sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Fibdiomyocyte")
# dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
# dbs1 <- split(dbs$gene, dbs$pathway)
# sub_metacell.obj <- AddModuleScore_UCell(sub_metacell.obj, features = dbs1,
#                                        ncores = 50)
# data.HCM <- AddModuleScore_UCell(HCM, features = dbs1,
#                                  ncores = 20)
library(vegan)
#install.packages('vegan')
library(psych)
library(reshape2)
library(ggplot2)
#data("varechem")
library(Seurat)
t_results <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/results/Fibroblast/UCell/C5_ALL_UCscore.rds")
t_results1 <- t_results[,c("sub_cell_type","HP_MYOCARDIAL_FIBROSIS_UCell","TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_UCell","RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_UCell","TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION_UCell","TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION_UCell","CELL_MIGRATION_UCell","CELL_MIGRATION_INVOLVED_IN_HEART_DEVELOPMENT_UCell","EXTRACELLULAR_MATRIX_ASSEMBLY_UCell","REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY_UCell","POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY_UCell")]
t_results2 <- subset(t_results1,subset = sub_cell_type %in% "Activated_fibroblast")
t_results2 <- t_results2[,-1]
sub_metascore.obj <- subset(metascore.obj,subset = cell_type %in% "Fibroblast")
Fib <- NormalizeData(sub_metascore.obj)
sub <- subset(Fib,subset = sub_cell_type %in% "Activated_fibroblast")
adata <- sub@meta.data
###提取出来之后把表达为0的值去掉-----
####相关性散点图-----
###相关性散点图----
###散点图----
#第1种画法
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
adata <- cbind(adata,t_results2)
colnames(adata)
# ###去除0值需要两条通路单独来做-----
# data <- adata[,c(1,26,168)]
# data1 <- data
# data1[data1 == 0] <- NA
# data1 <- na.omit(data1)
# data1 <- aggregate(data1, by=list(sample=data1$biosample_id),mean)##求均值
# colnames(data1)
######相关性散点图带山峦的图-----
P <- adata %>%
  ggplot(aes(x = Calcium.cycle_UCell, y = TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_UCell)) +
  geom_point(size = 1,color="#7db954",alpha = 0.5) + #"#efb306",  "#5797bc", "#e44349"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = Calcium.cycle_UCell, y = TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_UCell),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,  color = "black"),
        axis.title.y = element_text(size = 8,  color = "black")) +
  labs(y = "Transforming growth factor beta receptor signaling pathway")+
  NoLegend()
P
library(ggExtra)
p12 <- ggMarginal(P, type = "density", 
                  xparams = list(fill = "pink"),
                  yparams = list(fill = "#7db954"))
p12
ggsave(paste0(output,"Calcium.cycle和纤维化相关性散点图.pdf"), plot = p12, width = 5.39, height = 3.35)
####通路和基因的相关性-----
expr <- sub@assays$RNA@data
expr <- t(expr)
expr <- as.data.frame(expr)
adata1 <- cbind(adata,expr)
colnames(adata1)
# adata2 <- adata1[,c(26,24539)]
# adata2[adata2 == 0] <- NA
# adata2 <- na.omit(adata2)
P <- adata1 %>%
  ggplot(aes(x = Calcium.uniporter_UCell, y = SMAD2)) +
  geom_point(size = 1,color="#5797bc",alpha = 0.5) + #"#efb306",  "#5797bc", "#e44349"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = Calcium.uniporter_UCell, y = SMAD2),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,  color = "black"),
        axis.title.y = element_text(size = 8,  color = "black")) +
  labs(y = "NFKB1")+
  NoLegend()
P
library(ggExtra)
p12 <- ggMarginal(P, type = "density", 
                  xparams = list(fill = "pink"),
                  yparams = list(fill = "#5797bc"))
p12
ggsave(paste0(output,"Calcium.UNIPORTER和SMAD2相关性散点图.pdf"), plot = p12, width = 5.39, height = 3.35)
######根据钙离子转运分高低表达组------
library(UCell)
library(readxl)
library(readr)
#####按照基因集打分分高低组
library(stringr)
library(dplyr)
a <- sub_metascore.obj@meta.data
# a <- a[,-c(1,2,3,4,5,6,7,8)]
# a %>% mutate(across(where(is.character), as.numeric))  -> a
x <- "Calcium.uniporter_UCell"
median(a[,"Calcium.uniporter_UCell"])
sub_metascore.obj[[str_c(x, "_group")]] <-
  if_else(a[,x] > median(a[,x]),
          "high", "low")
sub_metascore.obj@meta.data
a1 <- sub@meta.data
x <- "Calcium.cycle_UCell"
median(a1[,"Calcium.cycle_UCell"])
sub[[str_c(x, "_group")]] <-
  if_else(a1[,x] > median(a1[,x]),
          "high", "low")
sub@meta.data
#####小提琴图-----
meta <- sub@meta.data
colnames(meta)
library(reshape2)
data_long <- melt(meta,
                  id.vars = c('Calcium.cycle_UCell_group'),#需要保留不参与聚合的变量,
                  variable.name='Pathways',
                  value.name='Score')
data_long$Score <- as.numeric(data_long$Score)
data_long <- na.omit(data_long)
a1 <- subset(data_long,subset = Pathways %in% "Complex.III_UCell") 
comparisons <- list(c("high","low"))
p <- ggviolin(a1, 
              "Calcium.cycle_UCell_group",
              "Score",
              color = "Calcium.cycle_UCell_group",
              add = "boxplot", 
              palette = c("#e44349","#5797bc" ),
              add.params = list(fill = "white"))+
  theme_classic2()+
  stat_compare_means(method = "t.test", 
                     label = "p.signif",##星号设置
                     comparisons = comparisons)+
  stat_compare_means(label.y = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab("Cluster")+ylab("FAO")+
  theme(plot.title = element_text(hjust = 0.3),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black") ) 
p
ggsave(paste0(output,"Calcium.UNIPORTER和SMAD2相关性散点图.pdf"), plot = p12, width = 5.39, height = 3.35)
#######细胞比例的桑吉图------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
Ratio <- sub_metascore.obj@meta.data %>%
  group_by(group, Calcium.uniporter_UCell_group) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


#堆叠柱状图
mycolor = c(
            '#e44349',
            '#5797bc',
            'pink',
            'green')

p15 <- ggplot(Ratio, aes(x =group, y= relative_freq, fill = Calcium.uniporter_UCell_group,
                         stratum=Calcium.uniporter_UCell_group, alluvium=Calcium.uniporter_UCell_group)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Disease',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p15
ggsave(paste0(output,"Calciumuniporter分组比例.pdf"), plot = p15, width = 5.17, height = 3.18)
#####高低表达组与细胞亚型聚类-----
####基因聚类----
library(tidyverse)
mito_gene <- unique(dbs$gene)
mito_gene <- subset(mito_gene,subset = mito_gene %in% rownames(sub_metascore.obj))
features <- mito_gene
sub_metascore.obj <- NormalizeData(sub_metascore.obj)
all_expr <-
  AverageExpression(sub_metascore.obj,
                    assays = "RNA",
                    slot = "data",
                    group.by = "sub_cell_type")[["RNA"]]
all_expr <- all_expr[!is.infinite(rowSums(all_expr)),]
group_expr <-
  AverageExpression(
    sub_metascore.obj,
    assays = "RNA",
    slot = "data",
    group.by = "Calcium.uniporter_UCell_group"
  )[["RNA"]]
all_expr <- as.data.frame(all_expr)
group_expr <- as.data.frame(group_expr)
all_expr %>% mutate(across(where(is.character), as.numeric))  -> all_expr
group_expr %>% mutate(across(where(is.character), as.numeric))  -> group_expr
expr <- cbind(all_expr,group_expr)
expr %>% mutate(across(where(is.character), as.numeric))  -> expr
library(RColorBrewer)
bk <- c(seq(-1, 1, by = 0.01))
p16 <- cor(expr[,c(1:3)], expr[,c(4:5)]) %>%
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
    color=colorRampPalette(c( 
                              '#5797bc',"#00C1D4", '#e44349',"firebrick3"))(length(bk)),
    # color = c(
    #   colorRampPalette(colors = c('#625D9E',"white"))(length(bk) / 2),
    #   colorRampPalette(colors = c("white",'#E95C59'))(length(bk) / 2)
    # ),
    legend_breaks = seq(-1, 1, 0.5),
    breaks = bk,
    main = "Correlation"
  )
p16
ggsave(paste0(output,"细胞类型相关性.pdf"), plot = p16, width = 3.53, height = 2.02)
#######不同cluster功能的热图------
library(pheatmap)
expr <- as.data.frame(sub_metascore.obj[["RNA"]]@data)
features <- c("BGN","FN1","TAGLN","ELN","LAMA1","ITGA1","ITGA2",     ####ECM
              "TGFB1","TGFB2","SMAD2","SMAD4","COL1A1","COL1A2","COL3A1"     #TGFb
)
expr1 <- expr[features,]
expr1 <- t(expr1)
meta <- sub_metascore.obj@meta.data
adata <- cbind(expr1,meta)
colnames(adata)
adata <- adata[,c(179,1:14)]
adata1 <- aggregate(adata, by=list(cluster=adata$Calcium.uniporter_UCell_group),mean)
rownames(adata1) <- adata1[,1]
adata1 <- adata1[,-c(1,2)]
bk <- c(seq(-1, 1, by = 0.01))
p17 <-pheatmap(
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
library(reshape2)
data_long <- melt(adata,
                  id.vars = c('Calcium.uniporter_UCell_group'),#需要保留不参与聚合的变量,
                  variable.name='Genes',
                  value.name='Score')
data_long$Score <- as.numeric(data_long$Score)
data_long <- na.omit(data_long)
for (i in features) {
  a1 <- subset(data_long,subset = Genes %in% i) 
  comparisons <- list(c("high","low"))
  p <- ggviolin(a1, 
                "Calcium.uniporter_UCell_group",
                "Score",
                color = "Calcium.uniporter_UCell_group",
                add = "boxplot", 
                palette = c("#e44349","#5797bc"),
                add.params = list(fill = "white"))+
    theme_classic2()+
    stat_compare_means(method = "t.test", 
                       label = "p.signif",##星号设置
                       comparisons = comparisons)+
    stat_compare_means(label.y = 2.2)+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab("Cluster")+ylab(i)+
    theme(plot.title = element_text(hjust = 0.3),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black") ) 
  p
  ggsave(paste0(output,i,"成纤维钙离子高低组小提琴.pdf"), plot = p, width = 6.47, height = 4.25)
}
####激活成纤维的小提琴图-----
expr <- as.data.frame(sub[["RNA"]]@data)
features <- c("BGN","FN1","TAGLN","ELN","LAMA1","ITGA1","ITGA2",     ####ECM
              "TGFB1","TGFB2","SMAD2","SMAD4","COL1A1","COL1A2","COL3A1"     #TGFb
)
expr1 <- expr[features,]
expr1 <- t(expr1)
meta <- sub@meta.data
adata <- cbind(expr1,meta)
colnames(adata)
adata <- adata[,c(179,1:14)]
library(reshape2)
data_long <- melt(adata,
                  id.vars = c('Calcium.uniporter_UCell_group'),#需要保留不参与聚合的变量,
                  variable.name='Genes',
                  value.name='Score')
data_long$Score <- as.numeric(data_long$Score)
data_long <- na.omit(data_long)
for (i in features) {
  a1 <- subset(data_long,subset = Genes %in% i) 
  comparisons <- list(c("high","low"))
  p <- ggviolin(a1, 
                "Calcium.uniporter_UCell_group",
                "Score",
                color = "Calcium.uniporter_UCell_group",
                add = "boxplot", 
                palette = c("#e44349","#5797bc"),
                add.params = list(fill = "white"))+
    theme_classic2()+
    stat_compare_means(method = "t.test", 
                       label = "p.signif",##星号设置
                       comparisons = comparisons)+
    stat_compare_means(label.y = 2.2)+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab("Cluster")+ylab(i)+
    theme(plot.title = element_text(hjust = 0.3),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black") ) 
  p
  ggsave(paste0(output,i,"激活成纤维钙离子高低组小提琴.pdf"), plot = p, width = 6.47, height = 4.25)
}
### 关注基因在轨迹上的表达-----
#####多组基因表达的曲线图-----
library(monocle)
cds <- qs::qread("/home/hutb/scmito/轨迹/Fibroblast/1e-31/cds.qs")
group_colors <- c('#efb306',
                  '#7db954',
                  '#852f88',
                  "firebrick",
                  '#4e54ac',
                  'pink'
)
s_gene <- 'ETV6'
p <- plot_genes_in_pseudotime(cds[s_gene,], color_by = "sub_cell_type")+
  scale_color_manual(values = group_colors)
df <- p$data
p18 <- ggplot(df, aes(Pseudotime, adjusted_expression)) + 
  geom_point(aes(colour = sub_cell_type), size = 0.5 )+
  geom_smooth(aes(colour = sub_cell_type),method = "loess", se = FALSE)+theme_bw()+#去除背景颜色
  theme(panel.grid = element_blank(),#去除背景网格
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 15))+
  scale_color_manual(values = group_colors) + labs(
    title = s_gene)
p18
ggsave(paste0(output,"ETV6轨迹表达图.pdf"), plot = p18, width = 5.46, height = 3.53)
s_gene <- 'EGR3'
p <- plot_genes_in_pseudotime(cds[s_gene,], color_by = "sub_cell_type")+
  scale_color_manual(values = group_colors)
df <- p$data
p19 <- ggplot(df, aes(Pseudotime, adjusted_expression)) + 
  geom_point(aes(colour = sub_cell_type), size = 0.5 )+
  geom_smooth(aes(colour = sub_cell_type),method = "loess", se = FALSE)+theme_bw()+#去除背景颜色
  theme(panel.grid = element_blank(),#去除背景网格
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 15))+
  scale_color_manual(values = group_colors) + labs(
    title = s_gene)
p19
ggsave(paste0(output,"EGR3轨迹表达图.pdf"), plot = p19, width = 5.46, height = 3.53)
s_gene <- 'SREBF2'
p <- plot_genes_in_pseudotime(cds[s_gene,], color_by = "sub_cell_type")+
  scale_color_manual(values = group_colors)
df <- p$data
p19 <- ggplot(df, aes(Pseudotime, adjusted_expression)) + 
  geom_point(aes(colour = sub_cell_type), size = 0.5 )+
  geom_smooth(aes(colour = sub_cell_type),method = "loess", se = FALSE)+theme_bw()+#去除背景颜色
  theme(panel.grid = element_blank(),#去除背景网格
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 15))+
  scale_color_manual(values = group_colors) + labs(
    title = s_gene)
p19
ggsave(paste0(output,"BPTF轨迹表达图.pdf"), plot = p19, width = 5.46, height = 3.53)
####在成纤维细胞中再按照Calcium高低表达分组-----
###分高低表达组----
library(readxl)
library(stringr)
library(dplyr)
a <- sub_metascore.obj@meta.data
x <- "Calcium.uniporter_UCell"
median(a[,"Calcium.uniporter_UCell"])
sub_metascore.obj[[str_c(x, "_group")]] <-
  if_else(a[,x] > median(a[,x]),
          "high", "low")
sub_metascore.obj@meta.data
####高低表达组中GSEA的分析----
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
### 高低表达组的差异火山图 ----
outdir <-
  file.path("/home/gongfengcz/scRNA-heart-mitochodria/results/Fibroblast/Calcium高低/")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}


DefaultAssay(sub_metascore.obj)
DefaultAssay(sub_metascore.obj) <- "RNA"
Idents(sub_metascore.obj) <- "Calcium.uniporter_UCell_group"
markers <- FindMarkers(sub_metascore.obj,
                       ident.1 ="high",
                       ident.2 = "low",
                       min_pct = 0,
                       logfc.threshold = 0)
deg_2 <- markers
deg_2$group <- "Calcium.uniporter_UCell"

#### GSEA 
library(purrr)
deg_2$gene <- rownames(deg_2)
gsea_res <- deg_2 |>
  cat_gsea(category = "C5")
saveRDS(gsea_res, file.path(output, "Calcium.uniporter_C5_gsea.rds"))
term <- gsea_res@result
####可视化-----
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
mygene <- c("POSTN","COL14A1","COL3A1","COL12A1","THSD4","THBS2")
p20 <- gseaNb(object = gsea_res,
              geneSetID = 'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT',
              newGsea = T,
              addGene = mygene,
              addPval = T)
p20
ggsave(paste0(output,"GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENTGSEA图.pdf"), plot = p20, width = 6.42, height = 4.57)
mygene <- c("POSTN","COL14A1","COL3A1","COL12A1","THSD4","THBS2")
p20 <- gseaNb(object = gsea_res,
              geneSetID = 'GOBP_COLLAGEN_FIBRIL_ORGANIZATION',
              newGsea = T,
              addGene = mygene,
              addPval = T)
p20
ggsave(paste0(output,"GOBP_COLLAGEN_FIBRIL_ORGANIZATIONGSEA图.pdf"), plot = p20, width = 6.42, height = 4.57)
mygene <- c("MCUB","MICU3","MCU","MICU2","MICU1","SMDT1","VDAC1")
p21 <- gseaNb(object = gsea_res,
              geneSetID = 'GOBP_CALCIUM_IMPORT_INTO_THE_MITOCHONDRION',
              newGsea = T,
              addGene = mygene,
              addPval = T)
p21
ggsave(paste0(output,"GOBP_CALCIUM_IMPORT_INTO_THE_MITOCHONDRIONGSEA图.pdf"), plot = p21, width = 6.42, height = 4.57)
p21 <- gseaNb(object = gsea_res,
              geneSetID = 'GOBP_MITOCHONDRIAL_CALCIUM_ION_HOMEOSTASIS',
              newGsea = T,
              addGene = mygene,
              addPval = T)
p21
ggsave(paste0(output,"GOBP_MITOCHONDRIAL_CALCIUM_ION_HOMEOSTASISGSEA图.pdf"), plot = p21, width = 6.42, height = 4.57)
####转录因子分析----
#######转录因子分析 ------
library(readr)
library(SCENIC)
library(Seurat)
library(tidyverse)
metacore.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
sub_metascore.obj <- subset(metacore.obj,subset = cell_type %in% "Fibroblast")
tfs_targer <- read.csv("/home/gongfengcz/scRNA-heart-mitochodria/results/pyscenic/Fibroblastsub_metacell/tfs_targets.csv")
######转录因子可视化  #####
# 导入auc数据
library(AUCell)
regulonAUC <- importAUCfromText("/home/gongfengcz/scRNA-heart-mitochodria/results/pyscenic/Fibroblastsub_metacell/auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- sub_metascore.obj@meta.data
Idents(sub_metascore.obj) <- sub_metascore.obj$sub_cell_type
# 将auc得分矩阵添加进seurat
sub_metascore.obj[["scenic"]] <-
  CreateAssayObject(counts = getAUC(regulonAUC))
#### 1. 按照细胞类型进行热图可视化 ----
regulonActivity_byCellType <-
  sapply(split(rownames(cellInfo), Idents(sub_metascore.obj)),
         function(cells)
           rowMeans(counts[,cells]))
a <- data.frame(names = row.names(regulonActivity_byCellType), regulonActivity_byCellType)
rownames(regulonActivity_byCellType) <-gsub(",", "", rownames(regulonActivity_byCellType))
library(readxl)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
tfs_targer1 <- subset(tfs_targer,subset = target %in% dbs$gene)

regulonActivity_byCellType <- subset(regulonActivity_byCellType,subset = rownames(regulonActivity_byCellType) %in% tfs_targer1$tfs)

library(pheatmap)
pheatmap(
  regulonActivity_byCellType,
  scale = "row",
  # annotation_col=meta,
  color =  colorRampPalette(c("#4575B4","white","#D73027"))(100),
  breaks=seq(-1, 1, length.out = 100),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = "white",
  cellheight = 10,
  cellwidth = 20,
  fontsize = 10,
  angle_col = "45",
  cluster_cols = F,
  cluster_rows = T,
  filename = "/home/gongfengcz/scRNA-heart-mitochodria/figure/final/Figure4/mitoheatmap_by_cluster.pdf"
)
pdf(file.path("/home/gongfengcz/scRNA-heart-mitochodria/figure/final/Figure4/mitoheatmap_by_subcelltype.pdf"), width = 3,height = 3)
p22 <-pheatmap(
  regulonActivity_byCellType,
  scale = "row",#'none', 'row' or 'column'
  # annotation_col=meta,
  color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white"))(50),
            colorRampPalette(colors = c("white","#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),#红到蓝梯度
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
dev.off()
######clusterprofile转录因子靶基因功能富集网络图-------
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
target <- tfs_targer
NR3C1 <- subset(target,subset = tfs %in% "NR3C1(+)") %>% pull(target)
ETV6 <- subset(target,subset = tfs %in% "ETV6(+)") %>% pull(target)
NFATC1 <- subset(target,subset = tfs %in% "NFATC1(+)") %>% pull(target)
EGR3 <- subset(target,subset = tfs %in% "EGR3(+)") %>% pull(target)
SREBF2 <- subset(target,subset = tfs %in% "SREBF2(+)") %>% pull(target)

data <- list(NR3C1,ETV6,NFATC1,EGR3,SREBF2)
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
gene_cluster <- data
names(gene_cluster)=c("NR3C1","ETV6","NFATC1","EGR3","SREBF2")

xx <- compareCluster(gene_cluster, 
                     fun='enrichGO',
                     ont= 'BP',
                     OrgDb='org.Hs.eg.db' ,
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1
)


xx <- pairwise_termsim(xx)
a <- xx@compareClusterResult
a <- subset(a,subset = Description %in% c("Wnt signaling pathway","cell-cell signaling by wnt","regulation of supramolecular fiber organization","transforming growth factor beta receptor signaling pathway","response to transforming growth factor beta","cellular response to transforming growth factor beta stimulus","calcium-mediated signaling","calcium ion transmembrane transport","mitochondrial calcium ion homeostasis","regulation of stress fiber assembly","positive regulation of cell-cell adhesion","response to fibroblast growth factor","response to hypoxia","positive regulation of Wnt signaling pathway","integrin-mediated signaling pathway"))
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
  scale_fill_manual(values = c('#efb306',
                               '#7db954',
                               '#852f88',
                               "firebrick",
                               '#4e54ac'))
dev.off()
p23
ggsave(paste0(output,"转录因子功能富集网络图.pdf"), plot = p23, width = 6.68, height = 5.20)

####表型相关的WGCNA------
#######表型的WGCNA分析------
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(Seurat)
library(randomcoloR)
if (!require('R.utils')) install.packages('R.utils')
library(UCell)
# saveRDS(sub_metacell.obj,file = "~/scRNA-heart-mitochodria/results/Fibroblast/UCell_seurat.rds")
data <- as.data.frame(sub_metascore.obj[["RNA"]]@data)
write.csv(data,file = "/home/gongfengcz/scRNA-heart-mitochodria/results/Fibroblast/WGCNA121/Mitodata.csv")
#参数设置
library(stringr)
GSE="Mitodata"    #表达矩阵文件名称，不用后缀
C="nf"              #正常控制组名称
P="hcm"              #疾病实验组的名称
rt <- as.data.frame(AverageExpression(sub_metacell.obj, group.by = c("sample"), assays = "RNA")[["RNA"]])
name2 <- as.data.frame(colnames(rt))
rownames(name2) <- name2[,1]
name2$group <- str_split_fixed(name2$`colnames(rt)`,"_",4)[,4]
name2 <- as.data.frame(name2)
# rt=read.csv(paste0(GSE,".csv"),comment = "/",header=T,check.names=F)
rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp <- rt
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

######################################################################################################################1.数据准备
#分组
sample=name2
rt=rt[,rownames(sample)]
afcon=sum(sample[,2]==C)
#判断原始数据是否去了log
max(rt)
if(max(rt)>50) rt=log2(rt+1)     #rt最大值大于50则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

#未标准化,mar调整画布范围下左上右，自行调整哈
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.raw.pdf",width=5,height = 4)
par(cex = 0.7,mar=c(8,8,8,8))
if(ncol(rt)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt,las=2,col =cols ) ###绘图
dev.off()

#标准化
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.nor.pdf",width=5,height = 4.5)
par(cex = 0.5,mar=c(8,8,8,8))
if(ncol(rt1)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt1,las=2,col =cols ) ###绘图
dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file=paste0("1.","norexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#保留原始结果
rt3=rbind(ID=colnames(rt),rt)
write.table(rt3,file=paste0("1.","rawexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#######################################################################################################################2.差异分析
#如果需要使用未经矫正的数据则将下方数据去除
data=rt1
#data=rt

conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#limma差异标准流程
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
#保存所有基因的差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("2.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC)>1 & adj.P.Val < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0("2.","DIFF_less.xls"),sep="\t",quote=F,col.names=F)

###3.GSEA分析
deg=Diff
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$symbol 
head(geneList)

#开始GSEA富集分析
library(readxl)
geneset1 <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
colnames(geneset1) <- c("term","gene")
geneset1 <- subset(geneset1,subset = term %in% c("Calcium uniporter","Calcium homeostasis","Calcium cycle"))
geneset2 <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/成纤维细胞功能通路.xlsx")
geneset2 <- subset(geneset2,subset = geneset2$pathways %in% c("GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA","GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY","GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION","GOBP_POSITIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION","GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION"))
colnames(geneset2) <- c("term","gene")
geneset <- rbind(geneset1,geneset2)
kk2 <- GSEA(geneList, TERM2GENE=geneset2, verbose=FALSE,pAdjustMethod = "fdr") 
#保存GSEA结果
GSEAOUT=as.data.frame(kk2@result)
write.table(GSEAOUT,file="4.GSEAOUT.xls",sep="\t",quote=F,col.names=T)
#看一看有没有富集到结果
head(GSEAOUT)

#排序后分别取GSEA结果的前5个和后5个
col=distinctColorPalette(100)        #随机颜色，不满意可以多跑几遍
#排序后取前5个和后5个一起展示
num=5
pdf(paste0("4.","all_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[1],color = "red",title = rownames(kk2@result)[1])
dev.off()

#生成目录存放批量GSEA单图，取前5个和后5个分别绘制
afname=rownames(kk2@result)
afdDir <- paste0(getwd(),"/4.GSEA")           #路径必须中文
dir.create(afdDir)
for (j in afname) {
  pdf(paste0(afdDir,"/",j,".pdf"),width = 10,height = 10)
  dd=gseaplot2(kk2, geneSetID = j,color = "red",title = j)
  print(dd)
  dev.off()
}

#############################################################################################################5.WGCNA
library(GSVA)
###准备
afdir <- paste0(getwd(),"/5.WGCNA")           #路径必须中文
dir.create(afdir)
afrt=read.table(paste0("/home/gongfengcz/scRNA-heart-mitochodria/1.norexp_Mitodata.txt"),sep="\t",header=T,check.names=F)
afrt=as.matrix(afrt)
# rownames(afrt)=afrt[,1]
# exp=afrt[,2:ncol(afrt)]
dimnames=list(rownames(exp),colnames(exp))
afrt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
afrt=avereps(afrt)
geneSet=split(geneset$gene,geneset$term)
ssgseaScore=gsva(afrt, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
afssgsea=normalize(ssgseaScore)
#保留原始结果
aaafssgsea=rbind(ID=colnames(afssgsea),afssgsea)
write.table(aaafssgsea,file=paste0("1.","ssgsea_",GSE,".xls"),sep="\t",quote=F,col.names = F)

traitData=sample
traitData[,1]=traitData[,2]
traitData[,1]=ifelse(traitData[,1]==C,1,0)
traitData[,2]=ifelse(traitData[,2]==P,1,0)
#修改性状名称
colnames(traitData)=c(C,P)
afssgsea=t(afssgsea)
afssname=intersect(rownames(afssgsea),rownames(traitData))
afssgsea=afssgsea[afssname,]
traitData=traitData[afssname,]
traitData=cbind(traitData,afssgsea)
###############正式分析------
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
fpkm = read.table("/home/gongfengcz/scRNA-heart-mitochodria/1.norexp_Mitodata.txt",header=T,sep = "\t",check.names=F)#########file name can be changed#####数据文件名，根据实际修改，如果工作路径不是实际数据路径，需要添加正确的数据路径
rownames(fpkm)=fpkm[,1]
dim(fpkm)
names(fpkm)
datExpr0 = as.data.frame(t(fpkm[,-1]))
names(datExpr0) = fpkm[,1];
rownames(datExpr0) = names(fpkm[,-1])

datExpr0

##check missing value
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

##filter
meanFPKM=0.1  ####the threshold can be changed---过滤标准
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]  # for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)


filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file=paste0(afdir,"/FPKM_filter.xls"),row.names=F, col.names=T,quote=FALSE,sep="\t")

sampleTree = hclust(dist(datExpr0), method = "average")


pdf(file =paste0(afdir,"/1_sampleClustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
##abline(h = 15, col = "red")##是否选择剪切
dev.off()


### Determine cluster under the line
##clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
##table(clust)


### clust 1 contains the samples we want to keep.
##keepSamples = (clust==1)
##datExpr0 = datExpr0[keepSamples, ]


####载入性状数据
#Loading clinical trait data
for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.

#sizeGrWindow(12,12)
pdf(file=paste0(afdir,"/2_Sample dendrogram and trait heatmap.pdf"),width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", 
                    marAll = c(1, 12, 3, 1))
dev.off()







#############################network constr########################################
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(1:30)

# Call the network topology analysis function
datExpr0 <- as.matrix(datExpr0)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file=paste0(afdir,"/3_Scale independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower
softPower =sft$powerEstimate
#显示软阈值
print(softPower)

adjacency = adjacency(datExpr0, power = 9)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste0(afdir,"/4_Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 0, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file=paste0(afdir,"/5_Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste0(afdir,"/6_Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######模块剪切高度
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file=paste0(afdir,"/7_merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")




##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste0(afdir,"/8_Module-trait relationships.pdf"),width=15,height=12)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(25, 10, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)

######

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      #sizeGrWindow(7, 7)
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####
names(datExpr0)
probes = names(datExpr0)


#################export GS and MM############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)



####功能富集-----
####打分相关的分析-----
library(reshape2)
library(UCell)
library(Seurat)
library(readr)
library(corrplot)
modulegene <- read.csv("~/scRNA-heart-mitochodria/results/Fibroblast/WGCNA121/WGCNA/step6_gene_moduleColors.csv")
turquoise <- subset(modulegene,subset = module %in% "turquoise") %>% pull(gene)
sub_metacell.obj <- readRDS("~/scRNA-heart-mitochodria/data/HCM-NF.rds")
sub_metacell.obj <- subset(sub_metacell.obj,subset = cell_type %in% "Fibroblast")
sub_metacell.obj <- NormalizeData(sub_metacell.obj)
dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
turqu <- as.data.frame(turquoise)
turqu$pathway <- "turquoise"
colnames(turqu) <- c("gene","pathway")
turqu <- turqu[,c(2,1)]
dbs <- rbind(dbs,turqu)
dbs1 <- split(dbs$gene, dbs$pathway)
sub_metacell.obj <- AddModuleScore_UCell(sub_metacell.obj, features = dbs1,maxRank = 150000,
                                   ncores = 20)
meta <- sub_metacell.obj@meta.data
meta1 <- meta[,c("Calcium.uniporter_UCell","Calcium.cycle_UCell","turquoise_UCell")]
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
p24 <- corrplot(M,
                add = FALSE,
                type = 'lower',
                method = 'ellipse',
                col = colorRampPalette(c("blue", "white", "firebrick"))(20))
p24
print(p24)
ggsave(paste0(output,"转录因子功能富集网络图.pdf"), plot = p24, width = 6.68, height = 5.20)
