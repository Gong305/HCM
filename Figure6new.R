setwd("~/scRNA-heart-mitochodria")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure6/311/"
dir.create(output, recursive = TRUE)
######细胞降维图
library(Seurat)
seurat.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/HCM_NF311.rds")
sub_seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Endothelial cell")
metacell.obj <- qs::qread("/home/gongfengcz/scRNA-heart-mitochodria/data/Nor_metacell_seuratobj.qs")
metacell.obj@meta.data
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Endothelial_cell")
Idents(sub_seurat.obj) <- "cell_type_leiden0.6"
library(RColorBrewer)
cell_type_cols <- c('#efb306','#7db954',
                    '#852f88',"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999", 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00") 
p1 <- DimPlot(sub_seurat.obj, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Endothelial") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1
ggsave(paste0(output,"内皮亚群umap.pdf"), plot = p1, width = 5.47, height = 4.11)

#堆叠柱状图
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据
Ratio <- sub_seurat.obj@meta.data %>%
  group_by(group, cell_type_leiden0.6) %>% # 分组
  dplyr::summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
mycolor = c('#efb306',
            '#7db954',
            '#852f88',
            '#4e54ac',
            '#0f8096',
            'pink',
            'green')

p2 <- ggplot(Ratio, aes(x =group, y= relative_freq, fill = cell_type_leiden0.6,
                        stratum=cell_type_leiden0.6, alluvium=cell_type_leiden0.6)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='group',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor)
p2
ggsave(paste0(output,"内皮亚群比例.pdf"), plot = p2, width = 5.22, height = 3.37)
####featureplot-----
library(ComplexHeatmap)
library(scCustomize)

library(circlize)
gene<-c("VEGFA","EFNB2","IL1R1","SMAD1","COMT","MT-CYB")
library(scCustomize)
library(Nebulosa)
DefaultAssay(sub_seurat.obj) <- 'RNA'
p1 <- Plot_Density_Custom(seurat_object = sub_seurat.obj, features = "VEGFA")
p1
ggsave(file.path(output,'VEGFA-featureplot.pdf'),
       p1,
       height=4,
       width=5)
p1 <- Plot_Density_Custom(seurat_object = sub_seurat.obj, features = "ENG")
p1
ggsave(file.path(output,'ENG-featureplot.pdf'),
       p1,
       height=4,
       width=5)
p1 <- Plot_Density_Custom(seurat_object = sub_seurat.obj, features = "IL1R1")
p1
ggsave(file.path(output,'IL1R1-featureplot.pdf'),
       p1,
       height=4,
       width=5)
p1 <- Plot_Density_Custom(seurat_object = sub_seurat.obj, features = "MT-CYB")
p1
ggsave(file.path(output,'MT-CYB-featureplot.pdf'),
       p1,
       height=4,
       width=5)

######top_marker基因的表达热图------
library(Seurat)
library(ggsci)
library(dplyr) # 用于数据操作
library(tidyverse)
library(scRNAtoolVis)
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
Idents(sub_metacell.obj) <- "sub_cell_type"
markers <- FindAllMarkers(
  sub_metacell.obj,
  only.pos = TRUE
)
# markers <- readRDS("~/IVL/data/all.markers.rds")
markers2 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)
library(msigdbr)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
dbs2 <- Dataset[grepl("INFLAMMATORY|APOPTOSIS|", Dataset$gs_name, ignore.case = TRUE), ]
table(dbs2$gs_name)
gene1 <- subset(dbs2,subset = gs_name %in% c("GOBP_INFLAMMATORY_RESPONSE")) %>% pull(gene_symbol)

dbs <- read_csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
gene1 <- unique(dbs$gene)
annoGene <- gene1
features <- c("PIK3R3","ENG","KDR","VEGFA","EFNB2", #生血相关
              "HDAC9","IL1R1","SMAD1","CIITA", #炎症相关
              "COMT","MACROD1","MT-CYB","SLC25A36"
)
Idents(sub_seurat.obj) <- "cell_type_leiden0.6"

p3 <- AverageHeatmap(object = sub_seurat.obj,
                     markerGene = markers$gene,column_names_rot = 45,
                     cluster.order = c(
                       'Endothelial_I','Endothelial_II','Endothelial_III'
                     ),
                     annoCol = TRUE,
                     myanCol = c('#efb306',
                                 '#7db954',
                                 '#852f88'
                     ),
                     colseed = 127,annoColType = 'dark',showRowNames = F,markGenes = features,fontsize = 5)
p3
ggsave(paste0(output,"内皮亚群marker基因.pdf"), plot = p3, width = 6.65, height = 5.70)

E1 <- subset(markers,subset = cluster %in% "Endothelial_I")
E2 <- subset(markers,subset = cluster %in% "Endothelial_II")
E3 <- subset(markers,subset = cluster %in% "Endothelial_III")
#### GO分析 ----
library(clusterProfiler)
library(org.Hs.eg.db) # 鼠：org.Mm.eg.db；人：org.Hs.eg.db
library(dplyr)
E1 <- E1[order(E1$avg_log2FC, decreasing = T),] %>% pull("gene")
E2 <- E2[order(E2$avg_log2FC, decreasing = T),] %>% pull("gene")
E3 <- E3[order(E3$avg_log2FC, decreasing = T),] %>% pull("gene")

E1_bp <-
  enrichGO(
    E1,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
E2_bp <-
  enrichGO(
    E2,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
E3_bp <-
  enrichGO(
    E3,
    OrgDb = org.Hs.eg.db,
    # 如果是鼠要对应进行更换
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
E1_term <- E1_bp@result
E2_term <- E2_bp@result
E3_term <- E3_bp@result
term1 <- subset(E1_term,subset = Description %in% c("endothelial cell apoptotic process","regulation of extrinsic apoptotic signaling pathway","regulation of endothelial cell apoptotic process"))
term1$cluster <- "Endothelial I"
term2 <- subset(E2_term,subset = Description %in% c("regulation of angiogenesis","sprouting angiogenesis","regulation of inflammatory response"))
term2$cluster <- "Endothelial II"
term3 <- subset(E3_term,subset = Description %in% c("phospholipid metabolic process","medium-chain fatty acid metabolic process","regulation of carbohydrate metabolic process"))
term3$cluster <- "Endothelial III"
library(ggplot2)
library(ggsci)
pathway <- rbind(term1,term2,term3)
pathway$cluster<-factor(pathway$cluster,
                        levels = c("Endothelial I","Endothelial II","Endothelial III"), ordered = TRUE)
pathway$Description<-factor(pathway$Description,
                            levels = rev(unique(pathway$Description)),
                            ordered = TRUE)
p1<-ggplot(pathway, aes(x=-log10(pvalue),
                        y=Description,
                        fill=cluster,na.rm = FALSE))+
  geom_bar(stat="identity",na.rm = FALSE)+ scale_fill_manual(values = c('#efb306',
                                                                        '#7db954',
                                                                        '#852f88'))+
  geom_text(aes(x=-log10(pvalue),
                y=Description,
                label=Count),
            size=2.5,
            hjust="right",
            nudge_x=0.1)+
  labs(x="-log10(pvalue)",
       y="",
       title="Enriched GO terms")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

p1
ggsave(paste0(output,"marker功能富集柱状图.pdf"), plot = p1, width = 6.38, height = 6.74)

######HCM_NF打分的小提琴图-----
########分裂小提琴图-----
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratDisk)
library(msigdbr)
library(AUCell)
library(UCell)
library(pheatmap)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
filtered_Dataset <- Dataset[grep("angiogenesis|migration|proliferation|inflammatory|extracellular|apopto", Dataset$gs_name), ]
metascore.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
sub_metascore.obj <- subset(metascore.obj,subset = cell_type %in% "Endothelial_cell")
dbs1 <- split(filtered_Dataset$gene_symbol, filtered_Dataset$gs_name)
sub_metascore.obj <- AddModuleScore_UCell(sub_metascore.obj, features = dbs1,maxRank = 150000,
                                    ncores = 50)
### 小提琴图
library(ggpubr)
my_comparisons <- c("NF","HCM")
meta <- sub_metascore.obj@meta.data
colnames(meta) <- gsub("_UCell|GOBP_", "", colnames(meta))
colnames(meta) <- tolower(colnames(meta))
meta$group <- factor(meta$group,levels = c("NF","HCM"))
p9 <- ggviolin(meta, x='group',y="oxphos",
               xlab = 'Group', ylab = "OXPHOS",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$oxphos),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"OXPHOS小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 
meta$fatty.acid.oxidation
p9 <- ggviolin(meta, x='group',y="fatty.acid.oxidation",
               xlab = 'Group', ylab = "FAO",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$fatty.acid.oxidation),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"FAO小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 
meta$`regulation of acute inflammatory response`
p9 <- ggviolin(meta, x='group',y="regulation of acute inflammatory response",
               xlab = 'Group', ylab = "regulation of acute inflammatory response",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$'regulation of acute inflammatory response'),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"inflammatory_response小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 

meta$`cytokine production involved in inflammatory response`
p9 <- ggviolin(meta, x='group',y="cytokine production involved in inflammatory response",
               xlab = 'Group', ylab = "cytokine production involved in inflammatory response",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$'cytokine production involved in inflammatory response'),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"cytokin_inflammatory小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 
meta$`negative regulation of endothelial cell apoptotic process`
p9 <- ggviolin(meta, x='group',y="negative regulation of endothelial cell apoptotic process",
               xlab = 'Group', ylab = "negative regulation of endothelial cell apoptotic process",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$'negative regulation of endothelial cell apoptotic process'),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"endothelial_cell_apoptotic小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 
meta$`negative regulation of extrinsic apoptotic signaling pathway`
p9 <- ggviolin(meta, x='group',y="negative regulation of extrinsic apoptotic signaling pathway",
               xlab = 'Group', ylab = "negative regulation of extrinsic apoptotic signaling pathway",
               fill='group',palette = c("#3288BD",  "#D53E4F" ),add = 'boxplot',
               add.params = list(fill='white', width=0.05))+
  geom_hline(yintercept = mean(meta$'negative regulation of extrinsic apoptotic signaling pathway'),# 取得分的中值
             linetype = "dashed", size = 0.2, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = my_comparisons) +
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("NF","HCM")), # 选择你想在哪2组上添加标签
    step_increase = 0.1,
    test="t.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=T)
p9
ggsave(paste0(output,"extrinsic_apoptotic小提琴.pdf"), plot = p9, width = 6.85, height = 4.92) 
####通路与通路相关性散点图------
colnames(meta)
adata1 <- meta[,c(1,16:626)]
adata<- aggregate(adata1, by=list(sample=adata1$sample),mean)##求均值
rownames(adata) <- adata$sample
adata <- as.data.frame(adata)
filtered_adata <- adata|>
  as.data.frame()
filtered_adata <- filtered_adata[,-c(1,2)]
filtered_adata %>% mutate(across(where(is.character), as.numeric))  -> filtered_adata

###相关性散点图----
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
meta$`chronic inflammatory response`
P <- filtered_adata %>%
  ggplot(aes(x = oxphos, y = `chronic inflammatory response`)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = oxphos, y = `chronic inflammatory response`),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Chronic inflammatory response",x = "OXPHOS")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#3288BD"),
                 yparams = list(fill = "#D53E4F"))
p4
ggsave(paste0(output,"OXPHOS与慢性炎症散点图.pdf"), plot = p4, width = 3.83, height = 3.35)
P <- filtered_adata %>%
  ggplot(aes(x = oxphos, y = `negative regulation of endothelial cell apoptotic process`)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = oxphos, y = `negative regulation of endothelial cell apoptotic process`),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Negative regulation of endothelial cell apoptotic process",x = "OXPHOS")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#3288BD"),
                 yparams = list(fill = "#D53E4F"))
p4
ggsave(paste0(output,"OXPHOS与endothelial_cell_apoptotic散点图.pdf"), plot = p4, width = 3.83, height = 3.35)

P <- filtered_adata %>%
  ggplot(aes(x = fatty.acid.oxidation, y = `chronic inflammatory response`)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = fatty.acid.oxidation, y = `chronic inflammatory response`),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Chronic inflammatory response",x = "FAO")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#3288BD"),
                 yparams = list(fill = "#D53E4F"))
p4
ggsave(paste0(output,"FAO与慢性炎症散点图.pdf"), plot = p4, width = 3.83, height = 3.35)

P <- filtered_adata %>%
  ggplot(aes(x = fatty.acid.oxidation, y = `negative regulation of extrinsic apoptotic signaling pathway`)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = fatty.acid.oxidation, y = `negative regulation of extrinsic apoptotic signaling pathway`),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Negative regulation of extrinsic apoptotic signaling pathway",x = "FAO")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#3288BD"),
                 yparams = list(fill = "#D53E4F"))
p4
ggsave(paste0(output,"FAO与凋亡散点图.pdf"), plot = p4, width = 3.83, height = 3.35)




#######轨迹分析------
#####monocle2轨迹-------
library(monocle)
#####monocle2-----
outdir <- output
cds1 <- qs::qread("~/scRNA-heart-mitochodria/figure/final/轨迹/Endothelial/1e-31/cds.qs")
cds1 <- orderCells(cds1,root_state = 3)
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
###轨迹饼状图------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=Biobase::pData(cds1)

# 计算每个状态下每种细胞类型的数量
counts <- table(plotdf$sub_cell_type,plotdf$State)
# 计算每个状态下每种细胞类型的比例
prop <- prop.table(counts, margin = 2) * 100  # 转换为百分比

# 绘制饼图
par(mfrow=c(1,3))  # 将绘图区域划分为一行三列

# 遍历每个状态
for (i in 1:3) {
  # 提取当前状态下的细胞类型比例
  prop_sub <- prop[, i]
  # 绘制饼图
  pie(prop_sub, main = paste("State", i), col=mycolor)
}

#####轨迹山峦图------
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=Biobase::pData(cds1)
plotdf <- plotdf[,c(17,18,19)]
a <- sub_metascore.obj@meta.data
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
  geom_vline(xintercept = c(3.2,20.6,21.4),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(c("#779fd3","white","firebrick"))(100))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
p6
ggsave(paste0(output,"轨迹山峦2.pdf"), plot = p6, width = 7.41, height = 5.54)

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
Time_diff <- read.csv("~/scRNA-heart-mitochodria/figure/final/轨迹/Endothelial/1e-31/diff_test_res.csv",row.names = 1)
dbs <- read.csv("/home/gongfengcz/oocytes-mitochondria/database/human_mito_v3.0.csv")
mito_gene <- unique(dbs$gene)
mito_gene1 <- subset(mito_gene,subset = mito_gene %in% rownames(cds1))
genes <- c("ACADSB","ATP5F1E","ATP5IF1","COXAL1","MT-ATP6","MT-CO1","MT-ND1","NDUFA1","UQCRB","CXCL12","IL33","CXCL1","IL7","C7","CFI","CFD","DCN","TGLN","ACE2","BMP10","TNNL1","COMP","INHBA","SULF1")
genes <- subset(genes,subset = genes %in% rownames(cds1))
p <- plot_pseudotime_heatmap(cds1[row.names(subset(Time_diff,
                                                   qval < 0.05)),],
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
                                   cores = 10,branch_colors = c("#979797", "firebrick3", "navy"),
                                   return_heatmap = T,
                                   use_gene_short_name = T,hmcols = colorRampPalette(c("navy","white","firebrick3","firebrick3"))(100),
                                   show_rownames = F)
print(p1)
dev.off()
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
                        topn = 100)

# check
head(enrich[1:3,])
######功能富集柱状图-----
mytheme <- theme(axis.title = element_text(size = 8),axis.text = element_text(size = 8),plot.title = element_text(size = 8,
                                                                                                                  hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 8),legend.text = element_text(size = 8))
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
term1 <- subset(enrich,subset = group %in% "C1")
term1 <- subset(term1,subset = Description %in% c("extracellular matrix constituent secretion","vasculogenesis","regulation of vasculogenesis","cytoskeleton-dependent intracellular transport"))
term1$text_x <- rep(0.03,4)
p1 <- ggplot(data = term1, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = 1)) +
  scale_fill_continuous(low = "#779fd3", high = "#779fd3") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
ggsave(paste0(output,"轨迹1功能富集柱状图.pdf"), plot = p1, width = 4.45, height = 4.05)

term2 <- subset(enrich,subset = group %in% "C2")
term2<- subset(term2,subset = Description %in% c("immune response-activating signaling pathway","immune response-regulating signaling pathway","activation of immune response","intrinsic apoptotic signaling pathway"))
term2$text_x <- rep(0.03,4)
p2 <- ggplot(data = term2, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = 1)) +
  scale_fill_continuous(low = "#7db954", high = "#7db954") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p2
ggsave(paste0(output,"轨迹2功能富集柱状图.pdf"), plot = p2, width = 4.45, height = 4.05)


term3 <- subset(enrich,subset = group %in% "C3")
term3<- subset(term3,subset = Description %in% c("ATP metabolic process","aerobic respiration","energy derivation by oxidation of organic compounds","oxidative phosphorylation"))
term3$text_x <- rep(0.03,4)
p3 <- ggplot(data = term3, aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill = 1)) +
  scale_fill_continuous(low = "pink", high = "pink") +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1,panel.grid = element_blank(),axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_blank()) +
  geom_text(aes(x = text_x, label = Description),hjust = 0,size = 3)+ mytheme2
p3
ggsave(paste0(output,"轨迹3功能富集柱状图.pdf"), plot = p3, width = 4.45, height = 4.05)
####轨迹上的得分------
plotdf=Biobase::pData(cds1)
plotdf <- plotdf[,c(18,19)]
meta <- sub_metascore.obj@meta.data
plotdf <- cbind(meta,plotdf)
#####对通路进行打分-----
Biobase::pData(cds1) <- plotdf
p27 <- plot_cell_trajectory(cds1,
                            color_by = "Fatty.acid.oxidation_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 
p27
ggsave(paste0(output,"轨迹Fatty.acid.oxidation图.pdf"), plot = p27, width = 4.63, height = 3.10)

p28 <- plot_cell_trajectory(cds1,
                            color_by = "OXPHOS_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 
p28
ggsave(paste0(output,"轨迹OXPHOS图.pdf"), plot = p28, width = 4.63, height = 3.10)

p28 <- plot_cell_trajectory(cds1,
                            color_by = "`Negative regulation of intrinsic apoptotic signaling pathway_UCell`",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 
p28
ggsave(paste0(output,"轨迹凋亡图.pdf"), plot = p28, width = 4.63, height = 3.10)

p28 <- plot_cell_trajectory(cds1,
                            color_by = "`Regulation of cell migration involved in sprouting angiogenesis_UCell`", #，，Negative regulation of inflammatory response to antigenic stimulus_UCell，Negative regulation of activated t cell proliferation_UCell
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 
p28
ggsave(paste0(output,"轨迹血管新生图.pdf"), plot = p28, width = 4.63, height = 3.10)

######内皮亚群中功能的打分小提琴图------
library(org.Hs.eg.db)
library(GOSemSim)
####UCell打分---
library(UCell)
library(ggpubr)
library(ggplot2)
meta <- sub_metacell.obj@meta.data
# meta <- subset(meta,subset = meta$group %in% "HCM")
meta1 <- meta[,c("sub_cell_type","group",colnames(meta)[grep("ATP_SYNTHESIS", colnames(meta))],
                 colnames(meta)[grep("ELECTRON_TRANSPORT", colnames(meta))],colnames(meta)[grep("FATTY", colnames(meta))],colnames(meta)[grep("APOPTOTIC", colnames(meta))],colnames(meta)[grep("CELL_DEATH", colnames(meta))])]
meta1 <- meta[,c("sub_cell_type","group",colnames(meta)[grep("FATTY", colnames(meta))])]
meta3 <- aggregate(meta1, by=list(type=meta1$sub_cell_type),mean)
meta3 <- t(meta3)
meta3 <- as.data.frame(meta3)
meta3 %>% mutate(across(where(is.character), as.numeric))  -> meta3
meta4 <- subset(meta3,subset = meta3$V1 > meta3$V3)
rownames(meta4)
colnames(meta1)

meta2 <- meta1[,c("sub_cell_type","GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_UCell")]
##ggviolin----
library(ggpubr)
type <- c("Endothelial_I","Endothelial_II","Endothelial_III")
compar <- list(c("Endothelial_I","Endothelial_II"),
               c("Endothelial_I","Endothelial_III"),
               c("Endothelial_III","Endothelial_II"))
p15<- ggviolin(meta, x='sub_cell_type',y="Negative regulation of intrinsic apoptotic signaling pathway_UCell",
               xlab = '', ylab = "Negative regulation of intrinsic apoptotic signaling pathway",
               fill='sub_cell_type',palette = c('#efb306',
                                                '#7db954',
                                                '#852f88'),add = 'boxplot',
               add.params = list(fill='white', width=0.05),order = type)+
  geom_hline(yintercept = mean(meta$`Negative regulation of intrinsic apoptotic signaling pathway_UCell`),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = compar)

p15
ggsave(paste0(output,"小提琴Negative regulation of intrinsic apoptotic signaling pathway图.pdf"), plot = p15, width = 7.02, height = 4.93) 
p16<- ggviolin(meta, x='sub_cell_type',y="OXPHOS_UCell",
               xlab = '', ylab = "OXPHOS",
               fill='sub_cell_type',palette = c('#efb306',
                                                '#7db954',
                                                '#852f88'),add = 'boxplot',
               add.params = list(fill='white', width=0.05),order = type)+
  geom_hline(yintercept = mean(meta$`OXPHOS_UCell`),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = compar)

p16
ggsave(paste0(output,"小提琴OXPHOS图.pdf"), plot = p16, width = 7.02, height = 4.93) 
p17<- ggviolin(meta, x='sub_cell_type',y="Regulation of cell migration involved in sprouting angiogenesis_UCell",
               xlab = '', ylab = "`Regulation of cell migration involved in sprouting angiogenesis_UCell`",
               fill='sub_cell_type',palette = c('#efb306',
                                                '#7db954',
                                                '#852f88'),add = 'boxplot',
               add.params = list(fill='white', width=0.05),order = type)+
  geom_hline(yintercept = mean(meta$`Regulation of cell migration involved in sprouting angiogenesis_UCell`),# 取得分的中值
             linetype = "dashed", size = 0.3, color = "black")+# 调整中间横线的粗细
  stat_compare_means(comparisons = compar)

p17
ggsave(paste0(output,"小提琴血管再生图.pdf"), plot = p17, width = 7.02, height = 4.93) 
#####Cytotrace------
#### 导入包
library(Seurat)
library(dplyr)
library(ggplot2)
library(CytoTRACE)
set.seed(717)
#### 导入数据
setwd("~/scRNA-heart-mitochodria/figure/final/Figure6/311/")
Idents(sub_metascore.obj) <- "sub_cell_type"
DimPlot(sub_metascore.obj,
        label = T,        
        repel = T,
        reduction = "umap")
#### 提取表型文件 
table(sub_metascore.obj$sub_cell_type)
phe <- sub_metascore.obj$sub_cell_type
phe = as.character(phe)
names(phe) <- rownames(sub_metascore.obj@meta.data)
#### 提取表达矩阵
mat_3k <- as.matrix(sub_metascore.obj@assays$RNA@counts)
mat_3k[1:4,1:4]
#### 当数据集中的单元数超过3,000时，CytoTRACE将自动以快速模式运行，这是一种用于减少运行时和内存使用的子采样方法。
#此外用户还可以使用ncores(默认值为1)来多线程，或者使用subsamplingsize(默认值为1000 cells)来指示子采样大小。
#在快速模式下运行以下数据集，使用8个核心，子样本大小为1000。
results <- CytoTRACE(mat = mat_3k)       #可以问一下?CytoTRACE 看看help
plotCytoGenes(results, numOfGenes = 10,outputDir = output)  #基因的数量可以改变
plotCytoTRACE(results, phenotype = phe, outputDir = paste0(output,"Macfour")) #默认用tsne降维，可以输入自己的嵌入
######转录因子网络------
library(readr)
library(SCENIC)
library(Seurat)
library(tidyverse)
tfs_targer <- read.csv("~/scRNA-heart-mitochodria/results/pyscenic/Endothelial_cellsub_metacell/tfs_targets.csv")
######转录因子可视化  #####
# 导入auc数据
library(AUCell)
regulonAUC <- importAUCfromText("~/scRNA-heart-mitochodria/results/pyscenic/Endothelial_cellsub_metacell/auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- sub_metascore.obj@meta.data
Idents(sub_metascore.obj) <- sub_metascore.obj$sub_cell_type
# 将auc得分矩阵添加进seurat
counts = t(getAUC(regulonAUC))
counts <- as.data.frame(counts)
counts$cell <- rownames(counts)
meta <- sub_metascore.obj@meta.data
meta$cell <- rownames(meta)
tf <- merge(meta,counts,by = "cell")
colnames(tf)
tf <- tf[,c(1,7,628:671)]
data1 <- aggregate(tf, by=list(sub_cell_type=tf$sub_cell_type),mean)##求均值
rownames(data1) <- data1[,1]
data1 <- data1[,-c(1,2,3)]
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
  color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
            colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),#红到蓝梯度
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
target <- read.csv("~/scRNA-heart-mitochodria/results/pyscenic/Endothelial_cellsub_metacell/tfs_targets.csv")
SREBF2 <- subset(target,subset = tfs %in% "SREBF2(+)") %>% pull(target)
SPIB <- subset(target,subset = tfs %in% "SPIB(+)") %>% pull(target)
IRF9 <- subset(target,subset = tfs %in% "IRF9(+)") %>% pull(target)
NR3C1 <- subset(target,subset = tfs %in% "NR3C1(+)") %>% pull(target)
ESRRA <- subset(target,subset = tfs %in% "ESRRA(+)") %>% pull(target)

data <- list(SREBF2,SPIB,IRF9,NR3C1,ESRRA)
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
gene_cluster <- data
names(gene_cluster)=c("SREBF2","SPIB","IRF9","NR3C1","ESRRA")

xx <- compareCluster(gene_cluster, 
                     fun='enrichGO',
                     ont= 'BP',
                     OrgDb='org.Hs.eg.db' ,
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1
)


xx <- pairwise_termsim(xx)
a <- xx@compareClusterResult
a <- subset(a,subset = Description %in% c("ATP synthesis coupled electron transport","mitochondrial ATP synthesis coupled electron transport","positive regulation of ATP-dependent activity","regulation of ATP-dependent activity","oxidative phosphorylation","negative regulation of innate immune response","negative regulation of cytokine-mediated signaling pathway","negative regulation of response to cytokine stimulus","negative regulation of extrinsic apoptotic signaling pathway"))
xx@compareClusterResult <- a  
pdf("~/scRNA-heart-mitochodria/results/Fibroblast/pyscenic/HCM_NF/go_net2.pdf",height = 8,width = 8)
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


### 关注基因在轨迹上的表达-----
#####多组基因表达的曲线图-----
cds <- qs::qread("/home/gongfengcz/scRNA-heart-mitochodria/results/Endothelial/monocle2/循环轨迹/1e-12/cds.qs")
group_colors <- c('#efb306',
                  '#7db954',
                  '#852f88',
                  "firebrick",
                  '#4e54ac',
                  'pink'
)
s_gene <- 'FLI1'
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
ggsave(paste0(output,"FLI1轨迹表达图.pdf"), plot = p18, width = 5.46, height = 3.53)


###线粒体通路打分的轨迹图-----
metascore.obj <- readRDS("~/scRNA-heart-mitochodria/data/metacell_score.rds")
metascore.obj <- subset(metascore.obj,subset = cell_type %in% "Endothelial_cell")
# outdir <- "~/scRNA-heart-mitochodria/results/Endothelial/monocle2/1e-12/"
cds <- qs::qread("~/scRNA-heart-mitochodria/figure/final/轨迹/Endothelial/1e-31/cds.qs")
# cds <- orderCells(cds1,root_state = 3)
plotdf=Biobase::pData(cds)
plotdf <- plotdf[,c(17,18,19)]
meta <- metascore.obj@meta.data
plotdf <- cbind(meta,plotdf)
sub_metacell.obj@meta.data <- plotdf
x <- "Fatty.acid.oxidation_UCell"
median(plotdf[,"Fatty.acid.oxidation_UCell"])
sub_metacell.obj[[str_c(x, "_group")]] <-
  if_else(plotdf[,x] > median(plotdf[,"Fatty.acid.oxidation_UCell"]),
          "high", "low")
plotdf <- sub_metacell.obj@meta.data
ggplot(plotdf, aes(x=Pseudotime,y=Fatty.acid.oxidation_UCell_group,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.5),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=20,face = "bold"))
#####对通路进行打分-----
Biobase::pData(cds) <- plotdf
p27 <- plot_cell_trajectory(cds,
                            color_by = "Fatty.acid.oxidation_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "F") 
p27
ggsave(paste0(output,"轨迹Fatty.acid.oxidation图.pdf"), plot = p27, width = 4.63, height = 3.10)

p28 <- plot_cell_trajectory(cds,
                            color_by = "OXPHOS_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "D") 
p28
ggsave(paste0(output,"轨迹OXPHOS图.pdf"), plot = p28, width = 4.63, height = 3.10)
p29 <- plot_cell_trajectory(cds,
                            color_by = "GOBP_NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "D") 
p29
ggsave(paste0(output,"轨迹NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE图.pdf"), plot = p29, width = 4.63, height = 3.10)
p30 <- plot_cell_trajectory(cds,
                            color_by = "GOBP_POSITIVE_REGULATION_OF_SPROUTING_ANGIOGENESIS_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "D") 
p30
ggsave(paste0(output,"轨迹POSITIVE_REGULATION_OF_SPROUTING_ANGIOGENESI图.pdf"), plot = p30, width = 4.63, height = 3.10)
p31 <- plot_cell_trajectory(cds,
                            color_by = "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_UCell",
                            size = 1,
                            show_backbone = TRUE) + viridis::scale_color_viridis(option = "D") 
p31
ggsave(paste0(output,"轨迹NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY图.pdf"), plot = p31, width = 4.63, height = 3.10)