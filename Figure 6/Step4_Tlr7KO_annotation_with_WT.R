
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/Tlr7/integration_with_WT")
Tlr7_integrated <- readRDS("Tlr7_integrated_PCA_0.1.rds")
Idents(Tlr7_integrated) <- "seurat_clusters"
#Remove cluster 11 - unknown, cluster 17 - contaminated astrocytes, cluster 18 - doublets
Tlr7_integrated <- subset(Tlr7_integrated, idents=c("11","17","18"), invert=T)
Idents(Tlr7_integrated) <- "seurat_clusters"
Tlr7_integrated <- RenameIdents(Tlr7_integrated,
                                `0` = "EN_Slc17a7", `1`="OL_Mbp", `2`="EN_Slc17a7", `3`="AST_Atp1a2",
                                `4`="EN_Slc17a7", `5`="EN_Slc17a7", `6`="MG_Cx3cr1", `7`="IN_Gad1",
                                `8`="OPC_Vcan", `9`="IN_Gad1", `10`="EN_Slc17a7",
                                `12`="EN_Slc17a7",`13`="VC_Cped1", `14`="CHOR_Ttr", `15`="EN_Slc17a7", `16`="IN_Gad1"
)

pdf("Tlr7_integrated_umap_annotation_test.pdf", width=7, height=5)
DimPlot(Tlr7_integrated, reduction = 'umap', label = T)
dev.off()

#Tlr7_integrated$celltype.orig.ident <- paste(Idents(Tlr7_integrated), Tlr7_integrated$orig.ident, sep = "_")
Tlr7_integrated$celltype <- Idents(Tlr7_integrated)
Tlr7_integrated$celltype <- factor(Tlr7_integrated$celltype, levels = c("EN_Slc17a7",
                                                                        "IN_Gad1","OL_Mbp","AST_Atp1a2","MG_Cx3cr1","OPC_Vcan","CHOR_Ttr","VC_Cped1"))
Idents(Tlr7_integrated) <- "celltype"
#markers for annotation
pdf("annotation_1.pdf", width=13, height=4)
DotPlot(Tlr7_integrated, features = c("Snap25","Syt1","Slc17a7", "Nrgn","Cux2","Pdzrn3","Mlip","Foxp2","Hs3st4","Grik3",
                                      "Cpne4","Tshz2","Vwc2l","Ntng2","Nr4a2","Col11a1","Gad1","Gad2","Plp1", "Mbp", "Mobp", "Aqp4", "Clu",
                                      "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra","Ttr","Htr2c", 
                                      "Cped1","Flt1","Ebf1")) + RotatedAxis()
dev.off()


saveRDS(Tlr7_integrated, file = "Tlr7_integrated_Annotation.rds")

# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(Tlr7_integrated$Condition,Tlr7_integrated$celltype))
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

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(), width=4,height=4,units="in")


# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(Tlr7_integrated$Sample_Name,Tlr7_integrated$celltype))
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

ggsave("sample_celltype_distribution.pdf",plot=last_plot(), width=6,height=4,units="in")

Idents(Tlr7_integrated) <- "celltype"
DefaultAssay(Tlr7_integrated) <- 'RNA'
pdf("Tlr7_integrated_umap_annotation_noLabel.pdf", width=7, height=5)
DimPlot(Tlr7_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(Tlr7_integrated, idents = "EN_Slc17a7")
Cluster_IN <- subset(Tlr7_integrated, idents = "IN_Gad1")
Cluster_MG <- subset(Tlr7_integrated, idents = "MG_Cx3cr1")
Cluster_AST <- subset(Tlr7_integrated, idents = "AST_Atp1a2")
Cluster_OL <- subset(Tlr7_integrated, idents = "OL_Mbp")
Cluster_OPC <- subset(Tlr7_integrated, idents = "OPC_Vcan")
Cluster_VC <- subset(Tlr7_integrated, idents = c("VC_Cped1"))
Cluster_CHOR <- subset(Tlr7_integrated, idents = "CHOR_Ttr")

saveRDS(Cluster_EN, file = "Tlr7_EN_subset.rds")
saveRDS(Cluster_IN, file = "Tlr7_IN_subset.rds")
saveRDS(Cluster_MG, file = "Tlr7_MG_subset.rds")
saveRDS(Cluster_AST, file = "Tlr7_AST_subset.rds")
saveRDS(Cluster_OL, file = "Tlr7_OL_subset.rds")
saveRDS(Cluster_OPC, file = "Tlr7_OPC_subset.rds")
saveRDS(Cluster_VC, file = "Tlr7_VC_subset.rds")
saveRDS(Cluster_CHOR, file = "Tlr7_CHOR_subset.rds")


