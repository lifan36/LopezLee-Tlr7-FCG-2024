
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/FCG/integration")
FCG_integrated <- readRDS("FCG_integrated_PCA_0.1.rds")
#Remove the last cluster: 15, unknown cell type, mixed, 0.17% of tatal cells. discard
FCG_integrated <- subset(FCG_integrated, idents="15", invert=T)
Idents(FCG_integrated) <- "seurat_clusters"
FCG_integrated <- RenameIdents(FCG_integrated,
                                 `0` = "excitatory neurons", `1`="oligodendrocytes", `2`="excitatory neurons", `3`="astrocytes",
                                 `4`="excitatory neurons", `5`="microglia", `6`="excitatory neurons", `7`="excitatory neurons",
                                 `8`="inhibitory neurons", `9`="inhibitory neurons", `10`="OPCs", `11`="excitatory neurons",
                                 `12`="vascular cells", `13`="ependymal cells", `14`="inhibitory neurons"
)

setwd("/athena/ganlab/scratch/lif4001/FCG/integration")

pdf("FCG_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(FCG_integrated, reduction = 'umap', label = F)
dev.off()

FCG_integrated$celltype.orig.ident <- paste(Idents(FCG_integrated), FCG_integrated$orig.ident, sep = "_")
FCG_integrated$celltype <- Idents(FCG_integrated)

saveRDS(FCG_integrated, file = "FCG_integrated_Annotation.rds")

data <- FCG_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Genotype.Food,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/FCG/integration",
       width=4,height=4,units="in")

#markers for annotation
pdf("annotation_1.pdf", width=9, height=3)
DotPlot(data, features = c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Gad1", "Gad2","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1","Clic6","Ttr")) + RotatedAxis()
dev.off()


data <- FCG_integrated
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$Sample_Name,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/FCG/integration",
       width=8,height=4,units="in")

Idents(FCG_integrated) <- "celltype"
DefaultAssay(FCG_integrated) <- 'RNA'
pdf("FCG_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(FCG_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(FCG_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(FCG_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(FCG_integrated, idents = "microglia")
Cluster_AST <- subset(FCG_integrated, idents = "astrocytes")
Cluster_OL <- subset(FCG_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(FCG_integrated, idents = "OPCs")
Cluster_VC <- subset(FCG_integrated, idents = "vascular cells")
Cluster_EP <- subset(FCG_integrated, idents = "ependymal cells")

saveRDS(Cluster_EN, file = "FCG_EN_subset.rds")
saveRDS(Cluster_IN, file = "FCG_IN_subset.rds")
saveRDS(Cluster_MG, file = "FCG_MG_subset.rds")
saveRDS(Cluster_AST, file = "FCG_AST_subset.rds")
saveRDS(Cluster_OL, file = "FCG_OL_subset.rds")
saveRDS(Cluster_OPC, file = "FCG_OPC_subset.rds")
saveRDS(Cluster_VC, file = "FCG_VC_subset.rds")
saveRDS(Cluster_EP, file = "FCG_EP_subset.rds")






