
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/integration_with_female")
bothSex_integrated <- readRDS("bothSex_integrated_PCA_0.1.rds")
#Remove cluster 17 - Choroid plexus epithelial cells, cluster 9-mixed, cluster 18 - 7 cells only.
bothSex_integrated <- subset(bothSex_integrated, idents=c("9","17","18"), invert=T)
Idents(bothSex_integrated) <- "seurat_clusters"
bothSex_integrated <- RenameIdents(bothSex_integrated,
                                 `0` = "oligodendrocytes", `1`="excitatory neurons", `2`="astrocytes", `3`="microglia",
                                 `4`="excitatory neurons", `5`="excitatory neurons", `6`="excitatory neurons", `7`="excitatory neurons",
                                 `8`="inhibitory neurons", `10`="OPCs", `11`="inhibitory neurons", `12`="excitatory neurons",
                                `13`="vascular cells",`14`="excitatory neurons",`15`="excitatory neurons", `16`="excitatory neurons"
)

pdf("bothSex_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(bothSex_integrated, reduction = 'umap', label = F)
dev.off()

bothSex_integrated$celltype.orig.ident <- paste(Idents(bothSex_integrated), bothSex_integrated$orig.ident, sep = "_")
bothSex_integrated$celltype <- Idents(bothSex_integrated)

saveRDS(bothSex_integrated, file = "bothSex_integrated_Annotation.rds")

data <- bothSex_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Condition,data$celltype))
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

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),
       width=4,height=4,units="in")

#markers for annotation
pdf("annotation_1.pdf", width=10.5, height=2.7)
DotPlot(data, features = c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Gad1", "Gad2","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1","Clic6")) + RotatedAxis()
dev.off()


data <- bothSex_integrated
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

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),
       width=6,height=4,units="in")

Idents(bothSex_integrated) <- "celltype"
DefaultAssay(bothSex_integrated) <- 'RNA'
pdf("bothSex_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(bothSex_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(bothSex_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(bothSex_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(bothSex_integrated, idents = "microglia")
Cluster_AST <- subset(bothSex_integrated, idents = "astrocytes")
Cluster_OL <- subset(bothSex_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(bothSex_integrated, idents = "OPCs")
Cluster_VC <- subset(bothSex_integrated, idents = "vascular cells")


saveRDS(Cluster_EN, file = "bothSex_EN_subset.rds")
saveRDS(Cluster_IN, file = "bothSex_IN_subset.rds")
saveRDS(Cluster_MG, file = "bothSex_MG_subset.rds")
saveRDS(Cluster_AST, file = "bothSex_AST_subset.rds")
saveRDS(Cluster_OL, file = "bothSex_OL_subset.rds")
saveRDS(Cluster_OPC, file = "bothSex_OPC_subset.rds")
saveRDS(Cluster_VC, file = "bothSex_VC_subset.rds")


#######################################################################
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/integration_with_female/subclustering")
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.15)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'bothSex_MG_reclusted_res0.15.rds')
pdf("bothSex_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_MG_umap_Condition_res0.15.pdf", width=5, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_MG_umap_Sample_res0.15.pdf", width=15, height=5.5)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
bothSex_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_MG_markers, "bothSex_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.15.csv")
#######################################################################

MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.2)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'bothSex_MG_reclusted_res0.2.rds')
pdf("bothSex_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_MG_umap_Condition_res0.2.pdf", width=5, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_MG_umap_Sample_res0.2.pdf", width=15, height=5.5)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
bothSex_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_MG_markers, "bothSex_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.2.csv")
#######################################################################

MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.3)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'bothSex_MG_reclusted_res0.3.rds')
pdf("bothSex_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_MG_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_MG_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
bothSex_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_MG_markers, "bothSex_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.3.csv")
#######################################################################
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.15)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'bothSex_AST_reclusted_res0.15.rds')
pdf("bothSex_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_AST_umap_Condition.pdf", width=5, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_AST_umap_Sample.pdf", width=15, height=5.5)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
bothSex_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_AST_markers, "bothSex_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts.csv")
###############################################
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.2)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'bothSex_AST_reclusted_res0.2.rds')
pdf("bothSex_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_AST_umap_Condition_res0.2.pdf", width=5, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_AST_umap_Sample_res0.2.pdf", width=15, height=5.5)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
bothSex_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_AST_markers, "bothSex_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.2.csv")
#######################################################################

AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.3)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'bothSex_AST_reclusted_res0.3.rds')
pdf("bothSex_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_AST_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_AST_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
bothSex_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_AST_markers, "bothSex_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.3.csv")
#######################################################################
#######################################################################
OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.15)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'bothSex_OL_reclusted_res0.15.rds')
pdf("bothSex_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OL_umap_Condition.pdf", width=5, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OL_umap_Sample.pdf", width=15, height=5.5)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
bothSex_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OL_markers, "bothSex_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts.csv")
#######################################################################

OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.2)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'bothSex_OL_reclusted_res0.2.rds')
pdf("bothSex_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OL_umap_Condition_res0.2.pdf", width=5, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OL_umap_Sample_res0.2.pdf", width=15, height=5.5)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
bothSex_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OL_markers, "bothSex_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts_res0.2.csv")
#######################################################################

OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.3)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'bothSex_OL_reclusted_res0.3.rds')
pdf("bothSex_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OL_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OL_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
bothSex_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OL_markers, "bothSex_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts_res0.3.csv")
#######################################################################
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.15)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'bothSex_OPC_reclusted_res0.15.rds')
pdf("bothSex_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OPC_umap_Condition.pdf", width=5, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OPC_umap_Sample.pdf", width=15, height=5.5)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
bothSex_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OPC_markers, "bothSex_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts.csv")
#######################################################################

OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.2)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'bothSex_OPC_reclusted_res0.2.rds')
pdf("bothSex_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OPC_umap_Condition_res0.2.pdf", width=5, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OPC_umap_Sample_res0.2.pdf", width=15, height=5.5)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
bothSex_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OPC_markers, "bothSex_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.2.csv")
#######################################################################

OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.3)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'bothSex_OPC_reclusted_res0.3.rds')
pdf("bothSex_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_OPC_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_OPC_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
bothSex_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_OPC_markers, "bothSex_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.3.csv")
#######################################################################

IN <- Cluster_IN
DefaultAssay(IN) <- 'integrated'
IN <- ScaleData(IN, verbose = FALSE)
IN <- RunPCA(IN, features = VariableFeatures(object = IN), verbose = FALSE)
ElbowPlot(IN)
IN <- FindNeighbors(IN, dims = 1:15)
IN <- FindClusters(IN, resolution = 0.1)
IN <- RunUMAP(IN, dims = 1: 15)
# rename cluster
n <- dim(table(IN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
IN@active.ident <- plyr::mapvalues(x = IN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
IN@active.ident <- factor(IN@active.ident, levels=1:n)
saveRDS(IN, file = 'bothSex_IN_reclusted_res0.1.rds')
pdf("bothSex_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_IN_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_IN_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
bothSex_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_IN_markers, "bothSex_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subcluster_cell_counts_res0.3.csv")
#######################################################################

EN <- Cluster_EN
DefaultAssay(EN) <- 'integrated'
EN <- ScaleData(EN, verbose = FALSE)
EN <- RunPCA(EN, features = VariableFeatures(object = EN), verbose = FALSE)
ElbowPlot(EN)
EN <- FindNeighbors(EN, dims = 1:15)
EN <- FindClusters(EN, resolution = 0.1)
EN <- RunUMAP(EN, dims = 1: 15)
# rename cluster
n <- dim(table(EN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
EN@active.ident <- plyr::mapvalues(x = EN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
EN@active.ident <- factor(EN@active.ident, levels=1:n)
saveRDS(EN, file = 'bothSex_EN_reclusted_res0.1.rds')
pdf("bothSex_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_EN_umap_Condition_res0.3.pdf", width=5, height=5)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("bothSex_EN_umap_Sample_res0.3.pdf", width=15, height=5.5)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
bothSex_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(bothSex_EN_markers, "bothSex_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subcluster_cell_counts_res0.3.csv")
#######################################################################
