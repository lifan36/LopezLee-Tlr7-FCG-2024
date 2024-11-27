#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/Tlr7/integration_with_WT")
Cluster_EN <- readRDS("Tlr7_EN_subset.rds")
Cluster_IN <- readRDS("Tlr7_IN_subset.rds")
Cluster_MG <- readRDS("Tlr7_MG_subset.rds")
Cluster_AST <- readRDS("Tlr7_AST_subset.rds")
Cluster_OL <- readRDS("Tlr7_OL_subset.rds")
Cluster_OPC <- readRDS("Tlr7_OPC_subset.rds")
Cluster_VC <- readRDS("Tlr7_VC_subset.rds")
Cluster_CHOR <- readRDS("Tlr7_CHOR_subset.rds")

setwd("/athena/ganlab/scratch/lif4001/Tlr7/integration_with_WT/subclustering")

###################################################################################
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
saveRDS(MG, file = 'Tlr7_MG_reclusted_res0.15.rds')
pdf("Tlr7_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_MG_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_MG_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
Tlr7_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_MG_markers, "Tlr7_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "Tlr7_MG_subcluster_cell_counts_res0.15.csv")
##########
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
saveRDS(MG, file = 'Tlr7_MG_reclusted_res0.2.rds')
pdf("Tlr7_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_MG_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_MG_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
Tlr7_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_MG_markers, "Tlr7_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "Tlr7_MG_subcluster_cell_counts_res0.2.csv")
##########
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
saveRDS(MG, file = 'Tlr7_MG_reclusted_res0.3.rds')
pdf("Tlr7_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_MG_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_MG_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
Tlr7_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_MG_markers, "Tlr7_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "Tlr7_MG_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
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
saveRDS(AST, file = 'Tlr7_AST_reclusted_res0.15.rds')
pdf("Tlr7_AST_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_AST_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_AST_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
Tlr7_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_AST_markers, "Tlr7_AST_markers_res0.15.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "Tlr7_AST_subcluster_cell_counts_res0.15.csv")
##########
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
saveRDS(AST, file = 'Tlr7_AST_reclusted_res0.2.rds')
pdf("Tlr7_AST_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_AST_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_AST_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
Tlr7_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_AST_markers, "Tlr7_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "Tlr7_AST_subcluster_cell_counts_res0.2.csv")
##########
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
saveRDS(AST, file = 'Tlr7_AST_reclusted_res0.3.rds')
pdf("Tlr7_AST_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_AST_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_AST_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
Tlr7_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_AST_markers, "Tlr7_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "Tlr7_AST_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
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
saveRDS(OL, file = 'Tlr7_OL_reclusted_res0.15.rds')
pdf("Tlr7_OL_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OL_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OL_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
Tlr7_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OL_markers, "Tlr7_OL_markers_res0.15.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "Tlr7_OL_subcluster_cell_counts_res0.15.csv")
##########
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
saveRDS(OL, file = 'Tlr7_OL_reclusted_res0.2.rds')
pdf("Tlr7_OL_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OL_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OL_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
Tlr7_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OL_markers, "Tlr7_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "Tlr7_OL_subcluster_cell_counts_res0.2.csv")
##########
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
saveRDS(OL, file = 'Tlr7_OL_reclusted_res0.3.rds')
pdf("Tlr7_OL_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OL_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OL_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
Tlr7_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OL_markers, "Tlr7_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "Tlr7_OL_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
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
saveRDS(OPC, file = 'Tlr7_OPC_reclusted_res0.15.rds')
pdf("Tlr7_OPC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OPC_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OPC_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
Tlr7_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OPC_markers, "Tlr7_OPC_markers_res0.15.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "Tlr7_OPC_subcluster_cell_counts_res0.15.csv")
##########
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
saveRDS(OPC, file = 'Tlr7_OPC_reclusted_res0.2.rds')
pdf("Tlr7_OPC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OPC_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OPC_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
Tlr7_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OPC_markers, "Tlr7_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "Tlr7_OPC_subcluster_cell_counts_res0.2.csv")
##########
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
saveRDS(OPC, file = 'Tlr7_OPC_reclusted_res0.3.rds')
pdf("Tlr7_OPC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_OPC_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_OPC_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
Tlr7_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_OPC_markers, "Tlr7_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "Tlr7_OPC_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
CHOR <- Cluster_CHOR
DefaultAssay(CHOR) <- 'integrated'
CHOR <- ScaleData(CHOR, verbose = FALSE)
CHOR <- RunPCA(CHOR, features = VariableFeatures(object = CHOR), verbose = FALSE)
ElbowPlot(CHOR)
CHOR <- FindNeighbors(CHOR, dims = 1:15)
CHOR <- FindClusters(CHOR, resolution = 0.15)
CHOR <- RunUMAP(CHOR, dims = 1: 15)
# rename cluster
n <- dim(table(CHOR@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
CHOR@active.ident <- plyr::mapvalues(x = CHOR@active.ident, from = current.cluster.ids, to = new.cluster.ids)
CHOR@active.ident <- factor(CHOR@active.ident, levels=1:n)
saveRDS(CHOR, file = 'Tlr7_CHOR_reclusted_res0.15.rds')
pdf("Tlr7_CHOR_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_CHOR_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_CHOR_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
Tlr7_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_CHOR_markers, "Tlr7_CHOR_markers_res0.15.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "Tlr7_CHOR_subcluster_cell_counts_res0.15.csv")
##########
CHOR <- Cluster_CHOR
DefaultAssay(CHOR) <- 'integrated'
CHOR <- ScaleData(CHOR, verbose = FALSE)
CHOR <- RunPCA(CHOR, features = VariableFeatures(object = CHOR), verbose = FALSE)
ElbowPlot(CHOR)
CHOR <- FindNeighbors(CHOR, dims = 1:15)
CHOR <- FindClusters(CHOR, resolution = 0.2)
CHOR <- RunUMAP(CHOR, dims = 1: 15)
# rename cluster
n <- dim(table(CHOR@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
CHOR@active.ident <- plyr::mapvalues(x = CHOR@active.ident, from = current.cluster.ids, to = new.cluster.ids)
CHOR@active.ident <- factor(CHOR@active.ident, levels=1:n)
saveRDS(CHOR, file = 'Tlr7_CHOR_reclusted_res0.2.rds')
pdf("Tlr7_CHOR_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_CHOR_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_CHOR_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
Tlr7_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_CHOR_markers, "Tlr7_CHOR_markers_res0.2.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "Tlr7_CHOR_subcluster_cell_counts_res0.2.csv")
##########
CHOR <- Cluster_CHOR
DefaultAssay(CHOR) <- 'integrated'
CHOR <- ScaleData(CHOR, verbose = FALSE)
CHOR <- RunPCA(CHOR, features = VariableFeatures(object = CHOR), verbose = FALSE)
ElbowPlot(CHOR)
CHOR <- FindNeighbors(CHOR, dims = 1:15)
CHOR <- FindClusters(CHOR, resolution = 0.3)
CHOR <- RunUMAP(CHOR, dims = 1: 15)
# rename cluster
n <- dim(table(CHOR@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
CHOR@active.ident <- plyr::mapvalues(x = CHOR@active.ident, from = current.cluster.ids, to = new.cluster.ids)
CHOR@active.ident <- factor(CHOR@active.ident, levels=1:n)
saveRDS(CHOR, file = 'Tlr7_CHOR_reclusted_res0.3.rds')
pdf("Tlr7_CHOR_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_CHOR_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_CHOR_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
Tlr7_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_CHOR_markers, "Tlr7_CHOR_markers_res0.3.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "Tlr7_CHOR_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
VC <- Cluster_VC
DefaultAssay(VC) <- 'integrated'
VC <- ScaleData(VC, verbose = FALSE)
VC <- RunPCA(VC, features = VariableFeatures(object = VC), verbose = FALSE)
ElbowPlot(VC)
VC <- FindNeighbors(VC, dims = 1:15)
VC <- FindClusters(VC, resolution = 0.15)
VC <- RunUMAP(VC, dims = 1: 15)
# rename cluster
n <- dim(table(VC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
VC@active.ident <- plyr::mapvalues(x = VC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VC@active.ident <- factor(VC@active.ident, levels=1:n)
saveRDS(VC, file = 'Tlr7_VC_reclusted_res0.15.rds')
pdf("Tlr7_VC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(VC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_VC_umap_Condition_res0.15.pdf", width=10, height=5)
DimPlot(VC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_VC_umap_Sample_res0.15.pdf", width=12, height=8)
DimPlot(VC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(VC) <- 'RNA'
Tlr7_VC_markers <- FindAllMarkers(VC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_VC_markers, "Tlr7_VC_markers_res0.15.csv")
write.csv(table(VC$seurat_clusters, VC$Sample_Name), "Tlr7_VC_subcluster_cell_counts_res0.15.csv")
##########
VC <- Cluster_VC
DefaultAssay(VC) <- 'integrated'
VC <- ScaleData(VC, verbose = FALSE)
VC <- RunPCA(VC, features = VariableFeatures(object = VC), verbose = FALSE)
ElbowPlot(VC)
VC <- FindNeighbors(VC, dims = 1:15)
VC <- FindClusters(VC, resolution = 0.2)
VC <- RunUMAP(VC, dims = 1: 15)
# rename cluster
n <- dim(table(VC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
VC@active.ident <- plyr::mapvalues(x = VC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VC@active.ident <- factor(VC@active.ident, levels=1:n)
saveRDS(VC, file = 'Tlr7_VC_reclusted_res0.2.rds')
pdf("Tlr7_VC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(VC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_VC_umap_Condition_res0.2.pdf", width=10, height=5)
DimPlot(VC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_VC_umap_Sample_res0.2.pdf", width=12, height=8)
DimPlot(VC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(VC) <- 'RNA'
Tlr7_VC_markers <- FindAllMarkers(VC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_VC_markers, "Tlr7_VC_markers_res0.2.csv")
write.csv(table(VC$seurat_clusters, VC$Sample_Name), "Tlr7_VC_subcluster_cell_counts_res0.2.csv")
##########
VC <- Cluster_VC
DefaultAssay(VC) <- 'integrated'
VC <- ScaleData(VC, verbose = FALSE)
VC <- RunPCA(VC, features = VariableFeatures(object = VC), verbose = FALSE)
ElbowPlot(VC)
VC <- FindNeighbors(VC, dims = 1:15)
VC <- FindClusters(VC, resolution = 0.3)
VC <- RunUMAP(VC, dims = 1: 15)
# rename cluster
n <- dim(table(VC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
VC@active.ident <- plyr::mapvalues(x = VC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
VC@active.ident <- factor(VC@active.ident, levels=1:n)
saveRDS(VC, file = 'Tlr7_VC_reclusted_res0.3.rds')
pdf("Tlr7_VC_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(VC, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_VC_umap_Condition_res0.3.pdf", width=10, height=5)
DimPlot(VC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_VC_umap_Sample_res0.3.pdf", width=12, height=8)
DimPlot(VC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(VC) <- 'RNA'
Tlr7_VC_markers <- FindAllMarkers(VC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_VC_markers, "Tlr7_VC_markers_res0.3.csv")
write.csv(table(VC$seurat_clusters, VC$Sample_Name), "Tlr7_VC_subcluster_cell_counts_res0.3.csv")
###################################################################################
###################################################################################
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
saveRDS(EN, file = 'Tlr7_EN_reclusted_res0.1.rds')
pdf("Tlr7_EN_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_EN_umap_Condition_res0.1.pdf", width=10, height=5)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_EN_umap_Sample_res0.1.pdf", width=12, height=8)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
Tlr7_EN_markers <- FindAllMarkers(EN, logfc.threshold = 1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_EN_markers, "Tlr7_EN_markers_res0.1.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "Tlr7_EN_subcluster_cell_counts_res0.1.csv")
##########
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
saveRDS(IN, file = 'Tlr7_IN_reclusted_res0.1.rds')
pdf("Tlr7_IN_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_IN_umap_Condition_res0.1.pdf", width=10, height=5)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("Tlr7_IN_umap_Sample_res0.1.pdf", width=12, height=8)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
Tlr7_IN_markers <- FindAllMarkers(IN, logfc.threshold = 1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(Tlr7_IN_markers, "Tlr7_IN_markers_res0.1.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "Tlr7_IN_subcluster_cell_counts_res0.1.csv")
##########

Idents(Cluster_EN) <- "Condition"
Idents(Cluster_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"
Idents(Cluster_OPC) <- "Condition"

setwd("/athena/ganlab/scratch/lif4001/Tlr7/integration_with_WT/DEGs")

df <- FindMarkers(Cluster_EN, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "EN_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_IN, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "IN_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_MG, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OPC, ident.1 = "WT_XXO_CPZ", ident.2 = "WT_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OPC_WT_XXO_CPZ_vs_WT_XXO_Ctrl_DEGs.csv")

df <- FindMarkers(Cluster_EN, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "EN_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_IN, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "IN_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_MG, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OPC, ident.1 = "WT_XYT_CPZ", ident.2 = "WT_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OPC_WT_XYT_CPZ_vs_WT_XYT_Ctrl_DEGs.csv")

df <- FindMarkers(Cluster_EN, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "EN_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_IN, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "IN_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_MG, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OPC, ident.1 = "KO_XXO_CPZ", ident.2 = "KO_XXO_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OPC_KO_XXO_CPZ_vs_KO_XXO_Ctrl_DEGs.csv")

df <- FindMarkers(Cluster_EN, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "EN_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_IN, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "IN_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_MG, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "MG_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_AST, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "AST_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OL, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OL_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")
df <- FindMarkers(Cluster_OPC, ident.1 = "KO_XYT_CPZ", ident.2 = "KO_XYT_Ctrl", logfc.threshold = 0.1, min.pct = 0, only.pos = F,test.use = "MAST")
df <- df[df$p_val_adj < 0.05,]
write.csv(df, "OPC_KO_XYT_CPZ_vs_KO_XYT_Ctrl_DEGs.csv")

