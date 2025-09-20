setwd("~/scRNA-heart-mitochodria")
output <- "~/scRNA-heart-mitochodria/figure/final/Figure7/311/celltype/"
dir.create(output)
######不分高低表达的multicellchat------
library(CellChat)
HCM <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure7/311/subtype/HCM/cellchat.rds")
NF <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure7/311/subtype/NF/cellchat.rds")
NF <- updateCellChat(NF)
HCM <- updateCellChat(HCM)
object.list <- list(NF = NF, 
                    HCM = HCM)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)) #因为有报错，原本是FALSE
cellchat
# cellchat <- identifyCommunicationPatterns(HCM, pattern = "incoming")
# netAnalysis_river(HCM, pattern = "outgoing")
# HCM@net$count
# netAnalysis_dot(HCM, pattern = "outgoing")

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
p1 <- gg1 + gg2  
p1
ggsave(paste0(output,"通讯数量柱状图.pdf"), plot = p1, width = 5.52, height = 3.55)
#图二，网络图绘制两组的，粗细-相互作用强度，红色-在上调的信号通路，蓝色-下调 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
p2
ggsave(paste0(output,"通讯网络图.pdf"), plot = p2, width = 12.63, height = 8.77)
#图三，热图 
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
p3 <- gg1 + gg2
p3
ggsave(paste0(output,"通讯热图.pdf"), plot = p3, width = 9.85, height = 5.53)
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
p4 <- rankNet(cellchat,
              mode = "comparison", stacked = T,color.use = c("#266E94","#FF9797"), #原本没有comparison = c(1,2,3,4,5,6)，自己加上去的
              do.stat = FALSE,signaling = c("TIGIT","OCLN","AGT","NGF","TAC","PERIOSTIN","ANGPTL","CD86","MHC-I","CD45","MHC-II","IL16","COMPLEMENT","PTN","THBS","NOTCH","WNT","VEGF","TGFb","COLLAGEN")) +coord_flip() #原本是TURE
p4
ggsave(paste0(output,"配受体条形图.pdf"), plot = p4, width = 6.16, height = 6.97)
#图七，气泡图
p5 <- netVisual_bubble(cellchat,
                       sources.use = c(1,9,12),            #sources.use = 1,sources.use = c(1,2,3,4)配体的细胞类型选哪个
                       targets.use = c(11,2),       #targets.use = c(3:2),targets.use = c(1,2,3,4,5,6)受体的细胞类型选哪个到哪个，targets.use = c(1:4),c(4)是AT1
                       comparison = c(1,2), 
                       signaling = c("THBS","PERIOSTIN","MHC-II","IL16","CD45","COMPLEMENT","MHC-I"),
                       #max.dataset = 2,
                       #title.name = "Increased signaling in LS", 
                       angle.x = 45, remove.isolate = T)
p5
ggsave(paste0(output,"配受体气泡.pdf"), plot = p5, width = 4.50, height = 4.73)
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
# object.list <- list(NF = NF)
object.list <- list(HCM = HCM,
                    NF = NF)
gg1 <- netAnalysis_signalingRole_scatter(HCM) + scale_x_continuous(limits = c(0,100)) + scale_y_continuous(limits = c(0,75)) 
gg2 <- netAnalysis_signalingRole_scatter(NF) + scale_x_continuous(limits = c(0,100)) + scale_y_continuous(limits = c(0,75)) 
gg2+gg1 
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg) + scale_x_continuous(limits = c(0,5)) 
######带箱子的山峦图表达----
metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
library(readxl)
library(stringr)
library(tidyverse)
meta <- sub_metacell.obj@meta.data
x <- "OXPHOS_UCell"
median(meta[,x])
sub_metacell.obj[[str_c(x, "_group")]] <-
  if_else(meta[,x] > median(meta[,x]),
          "high", "low")
meta1 <- sub_metacell.obj@meta.data
meta3 <- meta1[,c("group","OXPHOS_UCell_group")]
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
expr <- expr[c("POSTN","THBS1","THBS4"),]
expr <- t(expr)
adata <- cbind(expr,meta3)
colnames(adata)
# [1] "COL4A1"                            "COL4A2"                           
# [3] "COL4A3"                            "COL4A4"                           
# [5] "disease"                           "Calcium.homeostasis_UCell_group"  
# [7] "Calcium.cycle_UCell"               "Calcium.homeostasis_UCell"        
# [9] "Calcium.uniporter_UCell"           "Calcium.homeostasis_UCell_group.1"
# 设置颜色

p.val <- t.test(POSTN ~ OXPHOS_UCell_group,
                data = adata)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab
######按照样本计算得到得分的山峦图------
######带箱子的山峦图表达----
#####按照基因集打分分高低组
metacell.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
library(readxl)
library(stringr)
library(tidyverse)
meta <- sub_metacell.obj@meta.data
x <- "Fatty.acid.oxidation_UCell"
median(meta[,x])
sub_metacell.obj[[str_c(x, "_group")]] <-
  if_else(meta[,x] > median(meta[,x]),
          "high", "low")
meta1 <- sub_metacell.obj@meta.data
meta3 <- meta1[,c("group","Fatty.acid.oxidation_UCell_group")]
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
expr <- expr[c("POSTN","THBS1","THBS4"),]
expr <- t(expr)
expr <- as.data.frame(expr)
expr$sample <- rownames(expr)
meta3$sample <- rownames(meta3)
adata <- merge(expr,meta3,by.x = "sample",by.y = "sample")
# adata <- cbind(expr1,meta3)
colnames(adata)
# [1] "COL4A1"                            "COL4A2"                           
# [3] "COL4A3"                            "COL4A4"                           
# [5] "disease"                           "Calcium.homeostasis_UCell_group"  
# [7] "Calcium.cycle_UCell"               "Calcium.homeostasis_UCell"        
# [9] "Calcium.uniporter_UCell"           "Calcium.homeostasis_UCell_group.1"
# 设置颜色
p.val <- t.test(THBS1 ~ Fatty.acid.oxidation_UCell_group,
                data = adata)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab
green <- "#ec3114"
cyan <- "#0f66b7"
blue <- "#1B90BE"

# 绘制上半部分密度图
p_top <- ggplot(adata, aes(x = THBS1, color = Fatty.acid.oxidation_UCell_group, fill = Fatty.acid.oxidation_UCell_group, alpha = 0.8)) +
  geom_density() +
  geom_vline(xintercept = c(0.3, 0.75), linetype = 2, colour = "white") +
  scale_color_manual(values = c("#ec3114", "#0f66b7")) +
  scale_fill_manual(values = c("#ec3114", "#0f66b7")) +
  theme_classic() +
  xlab("Expression of THBS1") + 
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()

p_top

p_bot <- ggplot(adata, aes(Fatty.acid.oxidation_UCell_group, THBS1, fill = Fatty.acid.oxidation_UCell_group)) + 
  geom_boxplot(aes(col = Fatty.acid.oxidation_UCell_group)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated AUC") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 0.5,
           y = max(adata$POSTN)-0.3,
           size = 3, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() # 翻转图像

# 用白色标记箱子的基本统计量
library(aplot)
dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot
p <- p_top %>% insert_bottom(p_bot, height = 0.4)
p
ggsave(paste0(output,"心肌——THBS1云雨图.pdf"), plot = p, width = 7.45, height = 4.83)
invisible(dev.off())
####THBS4
colnames(adata)
# [1] "COL4A1"                            "COL4A2"                           
# [3] "COL4A3"                            "COL4A4"                           
# [5] "disease"                           "Calcium.homeostasis_UCell_group"  
# [7] "Calcium.cycle_UCell"               "Calcium.homeostasis_UCell"        
# [9] "Calcium.uniporter_UCell"           "Calcium.homeostasis_UCell_group.1"
# 设置颜色
p.val <- t.test(THBS4 ~ Fatty.acid.oxidation_UCell_group,
                data = adata)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab
green <- "#ec3114"
cyan <- "#0f66b7"
blue <- "#1B90BE"

# 绘制上半部分密度图
p_top <- ggplot(adata, aes(x = THBS4, color = Fatty.acid.oxidation_UCell_group, fill = Fatty.acid.oxidation_UCell_group, alpha = 0.8)) +
  geom_density() +
  geom_vline(xintercept = c(0.3, 0.75), linetype = 2, colour = "white") +
  scale_color_manual(values = c("#ec3114", "#0f66b7")) +
  scale_fill_manual(values = c("#ec3114", "#0f66b7")) +
  theme_classic() +
  xlab("Expression of THBS1") + 
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()

p_top

p_bot <- ggplot(adata, aes(Fatty.acid.oxidation_UCell_group, THBS4, fill = Fatty.acid.oxidation_UCell_group)) + 
  geom_boxplot(aes(col = Fatty.acid.oxidation_UCell_group)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated AUC") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 0.5,
           y = max(adata$POSTN)-0.3,
           size = 3, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() # 翻转图像

# 用白色标记箱子的基本统计量
library(aplot)
dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot
p <- p_top %>% insert_bottom(p_bot, height = 0.4)
p
ggsave(paste0(output,"心肌——THBS4云雨图.pdf"), plot = p, width = 7.45, height = 4.83)
#######相关性的散点图--------
library(UCell)
library(readxl)
library(readr)
library(vegan)
#install.packages('vegan')
library(psych)
library(reshape2)
library(ggplot2)
#data("varechem")
library(Seurat)
###散点图----
library(UCell)
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
meta3 <- meta1[,c("group","sample","Fatty.acid.oxidation_UCell","OXPHOS_UCell")]
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
expr <- expr[c("POSTN","THBS1","THBS4"),]
expr <- t(expr)
expr <- as.data.frame(expr)
expr$cell <- rownames(expr)
meta3$cell <- rownames(meta3)
adata <- merge(expr,meta3,by.x = "cell",by.y = "cell")
colnames(adata)
data1 <- adata[,c(6,3,7)]
data1[data1 == 0] <- NA
data1 <- na.omit(data1)
# data1 <- aggregate(adata, by=list(sample=adata$sample),mean)##求均值
colnames(data1)
data1 <- data1[,-c(2,6,7)]
P <- data1 %>%
  ggplot(aes(x = THBS1, y = Fatty.acid.oxidation_UCell)) +
  geom_point(size = 1,color="#0f66b7",alpha = 0.2) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = THBS1, y = Fatty.acid.oxidation_UCell),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "FAO")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#0f66b7"),
                 yparams = list(fill = "#ec3114"))
p4
ggsave(paste0(output,"心肌FAO与THBS1相关性.pdf"), plot = p4, width = 3.83, height = 3.35)



####心肌-成纤维细胞云雨图-----
sub_metacell.obj <- subset(metacell.obj,subset = cell_type %in% "Fibroblast")
library(readxl)
library(stringr)
library(tidyverse)
meta <- sub_metacell.obj@meta.data
x <- "Calcium.uniporter_UCell"
median(meta[,x])
sub_metacell.obj[[str_c(x, "_group")]] <-
  if_else(meta[,x] > median(meta[,x]),
          "high", "low")
meta1 <- sub_metacell.obj@meta.data
meta3 <- meta1[,c("group","Calcium.uniporter_UCell_group")]
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
expr <- expr[c("ITGAV","ITGB1","ITGB3","ITGA3"),]
expr <- t(expr)
expr <- as.data.frame(expr)
expr$sample <- rownames(expr)
meta3$sample <- rownames(meta3)
adata <- merge(expr,meta3,by.x = "sample",by.y = "sample")
# adata <- cbind(expr1,meta3)
colnames(adata)
# [1] "COL4A1"                            "COL4A2"                           
# [3] "COL4A3"                            "COL4A4"                           
# [5] "disease"                           "Calcium.homeostasis_UCell_group"  
# [7] "Calcium.cycle_UCell"               "Calcium.homeostasis_UCell"        
# [9] "Calcium.uniporter_UCell"           "Calcium.homeostasis_UCell_group.1"
# 设置颜色
p.val <- t.test(ITGAV ~ Calcium.uniporter_UCell_group,
                data = adata)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab
green <- '#7db954'
cyan <- '#efb306'
blue <- "#1B90BE"

# 绘制上半部分密度图
p_top <- ggplot(adata, aes(x = ITGB1, color = Calcium.uniporter_UCell_group, fill = Calcium.uniporter_UCell_group, alpha = 0.8)) +
  geom_density() +
  # geom_vline(xintercept = c(0.3, 0.75), linetype = 2, colour = "white") +
  scale_color_manual(values = c(green, cyan)) +
  scale_fill_manual(values = c(green, cyan)) +
  theme_classic() +
  xlab("Expression of ITGB1") + 
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()

p_top

p_bot <- ggplot(adata, aes(Calcium.uniporter_UCell_group, ITGB1, fill = Calcium.uniporter_UCell_group)) + 
  geom_boxplot(aes(col = Calcium.uniporter_UCell_group)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated AUC") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 0.5,
           y = max(adata$POSTN)-0.3,
           size = 3, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() # 翻转图像

# 用白色标记箱子的基本统计量
library(aplot)
dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot
p <- p_top %>% insert_bottom(p_bot, height = 0.4)
p
ggsave(paste0(output,"成纤维——ITGAV云雨图.pdf"), plot = p, width = 7.45, height = 4.83)
###散点图----
library(UCell)
source("/home/gongfengcz/Retina/src/computing/custom_function.R")
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
sub_metacell.obj1 <- subset(metacell.obj,subset = cell_type %in% "Fibroblast")
sub_metacell.obj2 <- subset(metacell.obj,subset = cell_type %in% "Cardiomyocyte")
meta1 <- sub_metacell.obj1@meta.data
meta2 <- sub_metacell.obj2@meta.data
meta3 <- meta1[,c("sample","Calcium.uniporter_UCell","Calcium.cycle_UCell","Calcium.homeostasis_UCell","OXPHOS_UCell","Fatty.acid.oxidation_UCell")]
meta4 <- meta2[,c("sample","OXPHOS_UCell","Fatty.acid.oxidation_UCell")]
# sub_metacell.obj <- NormalizeData(sub_metacell.obj)
# expr <- as.data.frame(sub_metacell.obj[["RNA"]]@data)
# expr <- expr[c("ITGAV","ITGB1","ITGB3","ITGA3"),]
# expr <- t(expr)
# expr <- as.data.frame(expr)
# expr$cell <- rownames(expr)
meta3$cell <- rownames(meta3)
meta4$cell <- rownames(meta4)
meta3 <- aggregate(meta3, by=list(sample=meta3$sample),mean)##求均值
meta4 <- aggregate(meta4, by=list(sample=meta4$sample),mean)##求均值
meta3 <- meta3[,-c(2,8)]
meta4 <- meta4[,-c(2,5)]

adata <- merge(meta3,meta4,by.x = "sample",by.y = "sample")
colnames(adata)
# data1 <- adata[,c(5,7,8)]
# data1[data1 == 0] <- NA
# data1 <- na.omit(data1)

colnames(adata)
# data1 <- data1[,-c(3)]
P <- adata %>%
  ggplot(aes(x =OXPHOS_UCell.x , y = OXPHOS_UCell.y)) +
  geom_point(size = 1,color="#0f66b7",alpha = 0.2) + #color="#F1948A"
  # geom_smooth(method = "lm", formula = y~x, color = "#ee6470", fill = scales::alpha("#ee6470", 0.3))+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = OXPHOS_UCell.x, y = OXPHOS_UCell.y),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  labs(y = "Cardiomyocyte OXPHOS")+
  NoLegend()
library(ggExtra)
p4 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "#0f66b7"),
                 yparams = list(fill = "#ec3114"))
p4
ggsave(paste0(output,"心肌与成纤维OXPHOS.pdf"), plot = p4, width = 3.83, height = 3.35)
#########疾病数据库-----
disease <- read_excel("~/scRNA-heart-mitochodria/results/疾病数据库/心脏疾病基因.xlsx")
features <- c(disease$`Cardiac hypertrophy`,disease$`Coronary heart disease`,disease$`Heart failure`,disease$`Dilated cardiomyopathy`,disease$`Ventricular fibrillation`)
features <- na.omit(features)
deg <- read.csv("~/scRNA-heart-mitochodria/results/deg/HCM_NF.csv")
deg$change <- "down"
deg$change[deg$logFC >0] <-  "up"
deg1 <- subset(deg,subset = deg$Gene %in% features)
a <- deg1 %>% count(Gene,change,sort = TRUE)
deg.num <- a
######分比例饼状气泡图绘制------
library(ggplot2)
library(scatterpie)
library(RColorBrewer)
DATA <- read_excel("~/scRNA-heart-mitochodria/results/疾病数据库/疾病饼状图2.xlsx")
DATA$up <- as.numeric(DATA$up)
DATA$down <- as.numeric(DATA$down)
Color<-c("#4575B4","#D73027")
ggplot() +
  geom_scatterpie(data = DATA,aes(x=gene, y=disease,r=0.4),  
                  cols=colnames(DATA)[3:4]) +
  scale_fill_manual(name = "Change", values = Color)+
  theme(panel.background  = element_blank()
  )



# 首先，计算每个gene的出现次数
gene_counts <- table(DATA$gene)

# 将gene_counts转换为数据框，以便可以与原始数据合并
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("gene", "count")
DATA <- merge(DATA, gene_counts_df, by = "gene", all.x = TRUE)
# 使用ggplot2绘图
library(ggplot2)

p <- ggplot() +
  geom_scatterpie(data = DATA, aes(x = gene, y = disease, r = count/10), cols = colnames(DATA)[3:4]) + coord_fixed() +
  scale_fill_manual(name = "Change", values = Color) +
  theme(panel.background = element_blank()) +
  geom_vline(xintercept = seq_along(unique(DATA$gene)) - 0.5, linetype = "dashed", color = "gray", size = 0.5) +
  scale_x_discrete(labels = unique(DATA$gene))  # 设置x轴标签为唯一的gene值

p + geom_scatterpie_legend(DATA$count/10, x=20, y=1,n = 3)

########疾病基因与线粒体相关通路按照样本进行相关性分析------
library(vegan)
#install.packages('vegan')
library(psych)
library(reshape2)
library(ggplot2)
#data("varechem")

#####线粒体相关通路与纤维化相关通路的相关性------
library(readxl)
library(UCell)
metascore.obj <- readRDS("/home/gongfengcz/scRNA-heart-mitochodria/data/metacell_score.rds")
meta<- metascore.obj@meta.data
colnames(meta)
adata <- meta[,c(1,24:26,62,63,99,120,144)]
expr <- as.data.frame(metascore.obj[["RNA"]]@data)
expr <- log2(expr+1)
expr <- expr[c("ACE","ACTA1","ACTB","AGTR1","ALDH1A2","ATP2B1","CD36","GRK2","HIF1A","KCNE1","MAML3","MYH6","PPARG","SOD2","TCF7L2"),]
expr <- t(expr)
expr <- as.data.frame(expr)
adata$cell <- rownames(adata)
expr$cell <- rownames(expr)
adata1 <- merge(adata,expr,by = "cell")
adata1 <- aggregate(adata1, by=list(saple=adata1$sample),mean)
colnames(adata1)
cor_function_meta <- corr.test(adata1[,c(12:26)],adata1[,c(6,7,8,10,11)], method = "spearman", adjust = "fdr")#做相关性


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


p7 <- ggplot(allnew, aes(variable,names)) +
  geom_tile(aes(fill=r),color="black",size=0.7) +
  geom_text(aes(label=r), color="black", size=4) + # 把相关性添加进去
  geom_text(aes(label=sig), color="black", size=4,vjust = 1.8)+ # 把星号添加进去
  scale_fill_gradient2(low='#89b6d5', high='#e92b27',mid = '#faf9fa', limit=c(-1,1),  
                       name=paste0("*p < 0.05","\n\n","**p < 0.01","\n\n","***p < 0.001","\n\n","Correlation")) + # 把P值添加到图例，这一步非常巧妙
  labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
  theme(axis.text.x = element_text(size=10,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())
p7
ggsave(paste0(output,"致病基因与线粒体通路相关性.pdf"), plot = p7, width = 6.33, height = 5.10)





#######附图-----细胞之间的通讯配受体分析------
output <- "~/scRNA-heart-mitochodria/figure/final/Figure7/311/celltype/"
######multicellchat------
library(CellChat)
HCM <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure7/311/HCM/cellchat.rds")
NF <- readRDS("~/scRNA-heart-mitochodria/figure/final/Figure7/311/NF/cellchat.rds")
NF <- updateCellChat(NF)
HCM <- updateCellChat(HCM)
###保持通讯层面细胞类型数量一致-----
group.new = levels(HCM@idents)
NF <- liftCellChat(NF, group.new)

# NF <- liftCellChat(NF, group.new)
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
p1 <- gg1 + gg2  
ggsave(paste0(output,"通讯数量柱状图.pdf"), plot = p1, width = 4.66, height = 3.12)
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
                        sources.use = c(1),            #sources.use = 1,sources.use = c(1,2,3,4)配体的细胞类型选哪个
                        targets.use = c(3,4,5),       #targets.use = c(3:2),targets.use = c(1,2,3,4,5,6)受体的细胞类型选哪个到哪个，targets.use = c(1:4),c(4)是AT1
                        comparison = c(1,2),       
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
# object.list <- list(NF = NF)
object.list <- list(HCM = HCM,
                    NF = NF)
gg1 <- netAnalysis_signalingRole_scatter(HCM) + scale_x_continuous(limits = c(0,10)) + scale_y_continuous(limits = c(0,7)) 
gg2 <- netAnalysis_signalingRole_scatter(NF) + scale_x_continuous(limits = c(0,10)) + scale_y_continuous(limits = c(0,7)) 
gg2+gg1 
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
####配受体的散点图分上下调展示------
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3,5,7),  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat,sources.use = 1, targets.use = c(3,5,8),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat,sources.use = 1, targets.use = c(3,5,8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
# Chord diagram
pathways.show <- c("COLLAGEN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", signaling.name = paste(pathways.show, names(object.list)[i]))
}









