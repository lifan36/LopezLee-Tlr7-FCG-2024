###Figure S12 (FCG P301S APOE4 snRNAseq)
library(ggplot2)
library(Seurat)

##Fig S12A, same RDS obj as Fig 4
LG93E4_OL <- readRDS('LG93E4_OL_clustered.rds')
DefaultAssay(LG93E4_OL) <- 'RNA'
DimPlot(LG93E4_OL, reduction = "umap", split.by = "Condition", label = TRUE, ncol = 4)

##Fig S12B
# calculate ratio of each cell type cluster in each sample
a<-as.data.frame(table(LG93E4_OL$seurat_clusters,LG93E4_OL$Sample_Name))
colnames(a)<-c("clusters","genotype","cell.no")
agg<-aggregate(cell.no~genotype,a,sum)
a$genotype.total <- agg$cell.no[match(a$genotype,agg$genotype)]
a$ratio<-a$cell.no/a$genotype.total
groups <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S12/Spreadsheets/OL_Groups.csv", header = 1) #read in groups
a$groups <- groups$Group
write.csv(a, "OL_CellRatio.csv")

#calculate avg and sem in Excel, read back in
a <- read.csv("/Users/lifan/Desktop/data_analysis/FCG/GitHub/Sequencing/R scripts/Figure S12/Spreadsheets/OL_CellRatio.csv", header=1)
#if need to force order bars
a$gclass <- factor(a$groups, levels=c('XXO_NTG_E4', 'XXO_P301S_E4', 'XYO_NTG_E4', 'XYO_P301S_E4', 'XXT_NTG_E4', 'XXT_P301S_E4', 'XYT_NTG_E4', 'XYT_P301S_E4'))
#plot
ggplot(a, aes(clusters, ratio, fill=gclass))+
  geom_errorbar( aes(ymin=avg-sem, ymax=avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))

##Fig S12D
#need OL3 markers
Idents(LG93E4_OL) <- "seurat_clusters"

#renumber
n <- dim(table(LG93E4_OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG93E4_OL@active.ident <- plyr::mapvalues(x = LG93E4_OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
LG93E4_OL@active.ident <- factor(LG93E4_OL@active.ident, levels=1:n)

clust_markers <- FindAllMarkers(LG93E4_OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(clust_markers, "OL_clustermarkers.csv")

#discard the last cluster, <2%
LG93E4_OL <- subset(LG93E4_OL, ident="6", invert=T)
#look through OL3 markers to identify upregulated markers
#dotplot by subcluster
DotPlot(LG93E4_OL, features = c("Fth1", "Scd2", "Actb", "Mal", "Plp1", "Cnp","Cst3", "Hsp90ab1", 
                                                                "Cryab", "Syt11", "Aplp1", 
                                                                "Lamp1", "Ckb", "Ptgds", "Cmss1", "Apod",
                                                                "Eef1a1", "Ptma", "Prnp", "Trf"), cols="RdBu", scale.by = "size") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))

