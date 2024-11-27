###Fig. S2 (WT young 3wk CPZ)
library(Seurat)
library(ggplot2)

##Fig S2A
YoungCPZ <- readRDS("FCG_YoungCPZ_integrated_Annotation.rds")
DefaultAssay(YoungCPZ) <- "RNA"
Idents(YoungCPZ) <- "celltype"
DimPlot(YoungCPZ, reduction = "umap", label = TRUE)

##Fig S2B
YoungCPZ_OL <- readRDS("YoungCPZ_OL_7dims_res0.2.rds")

#renumber clusters to start at 1
n <- dim(table(YoungCPZ_OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
YoungCPZ_OL@active.ident <- plyr::mapvalues(x = YoungCPZ_OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
YoungCPZ_OL@active.ident <- factor(YoungCPZ_OL@active.ident, levels=1:n)

DimPlot(YoungCPZ_OL, split.by = "Condition", reduction = "umap", label = TRUE, ncol=2)

##Fig S2C
#calculate cell ratios
a<-as.data.frame(table(YoungCPZ_OL$seurat_clusters,YoungCPZ_OL$Sample_Name))
colnames(a)<-c("clusters","sample","cell.no")
agg<-aggregate(cell.no~sample,a,sum)
a$genotype.total <- agg$cell.no[match(a$sample,agg$sample)]
a$ratio<-a$cell.no/a$genotype.total
#make groups spreadsheet in Excel listing out genotype for each subcluster+sample (i.e. if 5 subclusters x 3 samples = 15 rows in that genotype group)
groups <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S2/Spreadsheets/OL_CellRatioGroups.csv", header = 1) #read in groups
a$groups <- groups$Group
write.csv(a, "OL_CellRatio.csv")
#calculate in Excel avg and sem, read back in
a <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S2/Spreadsheets/OL_CellRatio.csv", header=1)
#if need to force order bars
a$genotype <- factor(a$groups, levels=c('F_Ctrl', 'F_3wk', 'M_Ctrl', 'M_3wk'))
#plot
ggplot(a, aes(clusters, ratio, fill=genotype))+
  geom_errorbar( aes(ymin=avg-sem, ymax=avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#re-label X axis cluster number to start at 1 in Illustrator

##Fig S2D
YoungCPZ_MG <- readRDS("YoungCPZ_MG_3dims_res0.1.rds")

#renumber clusters to start at 1
n <- dim(table(YoungCPZ_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
YoungCPZ_MG@active.ident <- plyr::mapvalues(x = YoungCPZ_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
YoungCPZ_MG@active.ident <- factor(YoungCPZ_MG@active.ident, levels=1:n)

DimPlot(YoungCPZ_MG, reduction = "umap", label = TRUE)

##Fig S2E
#calculate cell ratios
a<-as.data.frame(table(YoungCPZ_MG$seurat_clusters,YoungCPZ_MG$Sample_Name))
colnames(a)<-c("clusters","sample","cell.no")
agg<-aggregate(cell.no~sample,a,sum)
a$genotype.total <- agg$cell.no[match(a$sample,agg$sample)]
a$ratio<-a$cell.no/a$genotype.total
groups <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S2/Spreadsheets/MG_CellRatioGroups.csv", header = 1) #read in groups
a$groups <- groups$Group
#make groups spreadsheet in Excel listing out genotype for each subcluster+sample (i.e. if 5 subclusters x 3 samples = 15 rows in that genotype group)
write.csv(a, "MG_CellRatio.csv")
#Excel calculate avg and sem, read back in
a <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S2/Spreadsheets/MG_CellRatio.csv", header=1)
#define order of bars in graph
a$genotype <- factor(a$groups, levels=c('F_Ctrl', 'F_3wk', 'M_Ctrl', 'M_3wk'))
#plot
ggplot(a, aes(clusters, ratio, fill=genotype))+
  geom_errorbar( aes(ymin=avg-sem, ymax=avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#re-label X axis cluster number to start at 1 in Illustrator

