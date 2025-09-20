#### Load packages ----
library(Seurat)
library(tidyverse)
library(ggplot2)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")

up <- read.csv("~/scRNA-heart-mitochodria/results/function/allup.csv")
up$type <- "up"
up <- subset(up,up$pvalue < 0.05)
upterm <- c(# 上调
  "chronic inflammatory response",
  "inflammasome complex assembly",
  "cell-matrix adhesion",
  "regulation of actin cytoskeleton organization",
  "regulation of cell-substrate adhesion",
  "collagen fibril organization",
  "transforming growth factor beta activation",
  "transforming growth factor beta production")
updata <- subset(up,up$Description %in% upterm)
updata <- subset(updata,updata$celltype == "Pseudo-Bulk")
updata$value <- -log(updata$pvalue)

down <- read.csv("~/scRNA-heart-mitochodria/results/function/alldown.csv")
down$type <- "down"
down <- subset(down,down$pvalue < 0.05)
downterm = c(
  # 下调
  "mitochondrial transport",
  "mitochondrial transmembrane transport",
  "mitochondrial RNA processing",
  "mitochondrion morphogenesis",
  "dicarboxylic acid metabolic process",
  "glucose metabolic process",
  "fatty acid metabolic process"
)
downdata <- subset(down,down$Description %in% downterm)
downdata <- subset(downdata,downdata$celltype == "Pseudo-Bulk")
downdata$value <- log(downdata$pvalue)

# 合并
dat <- rbind(updata,downdata)
dat_plot <- data.frame(id = dat$Description,
                       value = dat$value,
                       GeneRatio = dat$Count)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
dat_plot$type <- dat$type

dat_plot <- dat_plot[order(dat_plot$value,decreasing = T),]
p2 <- ggplot(data = dat_plot,aes(x = reorder(id, order(value, decreasing=F)),y = value,fill = type)) +
  geom_col()+
  coord_flip() +
  geom_text(aes(label = GeneRatio),size =5) + 
  scale_fill_manual(values = c('down'= '#779fd3','up'='#a13037')) +
  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('GOBP Pathway Enrichment') + 
  ylab("-log(pvalue)") 
p2
# 依次从下到上添加标签

up_pathway <- length(which(dat_plot$value > 0))
down_pathway <- length(which(dat_plot$value < 0))
high <- nrow(dat_plot)
p3 <- p2 + geom_text(data = dat_plot[1:up_pathway,],aes(x = id,y = 0.1,label = id),
                     hjust=1.2,color = 'black',size = 4) +
  geom_text(data = dat_plot[(up_pathway +1):high,],aes(x = id,y = -0.1,label = id),
            hjust = -0.2,color = 'black',size = 4) +
  scale_x_discrete(labels = NULL) +
  theme_bw() +theme(panel.grid=element_blank())
p3



