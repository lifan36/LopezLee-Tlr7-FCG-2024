
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration")
FCG_YoungCPZ_integrated <- readRDS("FCG_YoungCPZ_integrated_PCA_0.1.rds")
#Remove cluster 13 - Choroid plexus epithelial cells, cluster 14-Reln+ CR cells, 12-mixed
FCG_YoungCPZ_integrated <- subset(FCG_YoungCPZ_integrated, idents=c("12","13","14"), invert=T)
Idents(FCG_YoungCPZ_integrated) <- "seurat_clusters"
FCG_YoungCPZ_integrated <- RenameIdents(FCG_YoungCPZ_integrated,
                                 `0` = "excitatory neurons", `1`="oligodendrocytes", `2`="excitatory neurons", `3`="astrocytes",
                                 `4`="excitatory neurons", `5`="inhibitory neurons", `6`="microglia", `7`="excitatory neurons",
                                 `8`="excitatory neurons", `9`="OPCs", `10`="vascular cells", `11`="excitatory neurons",
                                 `15`="vascular cells"
)

pdf("FCG_YoungCPZ_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(FCG_YoungCPZ_integrated, reduction = 'umap', label = T)
dev.off()

FCG_YoungCPZ_integrated$celltype.orig.ident <- paste(Idents(FCG_YoungCPZ_integrated), FCG_YoungCPZ_integrated$orig.ident, sep = "_")
FCG_YoungCPZ_integrated$celltype <- Idents(FCG_YoungCPZ_integrated)

saveRDS(FCG_YoungCPZ_integrated, file = "FCG_YoungCPZ_integrated_Annotation.rds")

data <- FCG_YoungCPZ_integrated
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
                           "Cped1","Clic6","Ttr")) + RotatedAxis()
dev.off()


data <- FCG_YoungCPZ_integrated
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

Idents(FCG_YoungCPZ_integrated) <- "celltype"
DefaultAssay(FCG_YoungCPZ_integrated) <- 'RNA'
pdf("FCG_YoungCPZ_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(FCG_YoungCPZ_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(FCG_YoungCPZ_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(FCG_YoungCPZ_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(FCG_YoungCPZ_integrated, idents = "microglia")
Cluster_AST <- subset(FCG_YoungCPZ_integrated, idents = "astrocytes")
Cluster_OL <- subset(FCG_YoungCPZ_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(FCG_YoungCPZ_integrated, idents = "OPCs")
Cluster_VC <- subset(FCG_YoungCPZ_integrated, idents = "vascular cells")


saveRDS(Cluster_EN, file = "FCG_YoungCPZ_EN_subset.rds")
saveRDS(Cluster_IN, file = "FCG_YoungCPZ_IN_subset.rds")
saveRDS(Cluster_MG, file = "FCG_YoungCPZ_MG_subset.rds")
saveRDS(Cluster_AST, file = "FCG_YoungCPZ_AST_subset.rds")
saveRDS(Cluster_OL, file = "FCG_YoungCPZ_OL_subset.rds")
saveRDS(Cluster_OPC, file = "FCG_YoungCPZ_OPC_subset.rds")
saveRDS(Cluster_VC, file = "FCG_YoungCPZ_VC_subset.rds")


#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/subclustering")
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.15.pdf", width=10, height=6)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.15.pdf", width=18, height=12)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.15.csv")
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.2.pdf", width=10, height=6)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.2.pdf", width=18, height=12)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.2.csv")
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.3.pdf", width=10, height=6)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.3.pdf", width=18, height=12)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.3.csv")
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition.pdf", width=10, height=6)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample.pdf", width=18, height=12)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts.csv")
##################################################
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.2.pdf", width=10, height=6)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.2.pdf", width=18, height=12)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.2.csv")
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.3.pdf", width=10, height=6)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.3.pdf", width=18, height=12)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.3.csv")
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition.pdf", width=10, height=6)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample.pdf", width=18, height=12)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers.csv")
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.2.pdf", width=10, height=6)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.2.pdf", width=18, height=12)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.2.csv")
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.3.pdf", width=10, height=6)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.3.pdf", width=18, height=12)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.3.csv")
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition.pdf", width=10, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample.pdf", width=18, height=12)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers.csv")
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.2.pdf", width=10, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.2.pdf", width=18, height=12)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.2.csv")
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.3.pdf", width=10, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.3.pdf", width=18, height=12)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.3.csv")
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
saveRDS(IN, file = 'FCG_YoungCPZ_IN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Condition_res0.3.pdf", width=10, height=6)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Sample_res0.3.pdf", width=18, height=12)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_YoungCPZ_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_IN_markers, "FCG_YoungCPZ_IN_markers_res0.3.csv")
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
saveRDS(EN, file = 'FCG_YoungCPZ_EN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Condition_res0.1.pdf", width=10, height=6)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Sample_res0.1.pdf", width=18, height=12)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_YoungCPZ_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_EN_markers, "FCG_YoungCPZ_EN_markers_res0.1.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subcluster_cell_counts_res0.3.csv")
#######################################################################


Idents(Cluster_EN) <- "Condition"
Idents(Cluster_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"
Idents(Cluster_OPC) <- "Condition"
Idents(Cluster_VC) <- "Condition"

setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/DEGs")
#FCG_YoungCPZ vs Ctrl
EN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_EN, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "EN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
IN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_IN, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "IN_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
MG_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_MG, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "MG_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
AST_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_AST, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "AST_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
OL_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_OL, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "OL_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
OPC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "OPC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")
VC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs <- FindMarkers(Cluster_VC, ident.1 = "XXO_P301S_E4", ident.2 = "XXO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(VC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs, "VC_XXO_P301S_E4_vs_XXO_NTG_E4_DEGs.csv")

EN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_EN, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "EN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
IN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_IN, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "IN_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
MG_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_MG, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "MG_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
AST_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_AST, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "AST_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
OL_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_OL, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "OL_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
OPC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "OPC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")
VC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs <- FindMarkers(Cluster_VC, ident.1 = "XYO_P301S_E4", ident.2 = "XYO_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(VC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs, "VC_XYO_P301S_E4_vs_XYO_NTG_E4_DEGs.csv")

EN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_EN, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "EN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
IN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_IN, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "IN_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
MG_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_MG, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "MG_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
AST_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_AST, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "AST_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
OL_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_OL, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "OL_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
OPC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "OPC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")
VC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs <- FindMarkers(Cluster_VC, ident.1 = "XXT_P301S_E4", ident.2 = "XXT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(VC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs, "VC_XXT_P301S_E4_vs_XXT_NTG_E4_DEGs.csv")

EN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_EN, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "EN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
IN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_IN, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "IN_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
MG_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_MG, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "MG_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
AST_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_AST, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "AST_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
OL_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_OL, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "OL_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
OPC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "OPC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")
VC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs <- FindMarkers(Cluster_VC, ident.1 = "XYT_P301S_E4", ident.2 = "XYT_NTG_E4", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(VC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs, "VC_XYT_P301S_E4_vs_XYT_NTG_E4_DEGs.csv")


Idents(Cluster_EN) <- "Condition"
Idents(Cluster_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"
Idents(Cluster_OPC) <- "Condition"
Idents(Cluster_VC) <- "Condition"

XXO_P301S_EN <- subset(Cluster_EN, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_IN <- subset(Cluster_IN, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_MG <- subset(Cluster_MG, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_AST <- subset(Cluster_AST, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_OL <- subset(Cluster_OL, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_OPC <- subset(Cluster_OPC, idents=c("XXO_NTG_E4","XXO_P301S_E4"))
XXO_P301S_VC <- subset(Cluster_VC, idents=c("XXO_NTG_E4","XXO_P301S_E4"))

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/subclustering_XXO_P301S")
MG <- XXO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.15.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.15.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXO_P301S_cell_counts_res0.15.csv")
#######################################################################
MG <- XXO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXO_P301S_cell_counts_res0.2.csv")
#######################################################################
MG <- XXO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################
AST <- XXO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXO_P301S_cell_counts.csv")
#####################################################
AST <- XXO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXO_P301S_cell_counts_res0.2.csv")
#######################################################################
AST <- XXO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################
OL <- XXO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXO_P301S_cell_counts.csv")
#######################################################################
OL <- XXO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXO_P301S_cell_counts_res0.2.csv")
#######################################################################
OL <- XXO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################
OPC <- XXO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXO_P301S_cell_counts.csv")
#######################################################################
OPC <- XXO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXO_P301S_cell_counts_res0.2.csv")
#######################################################################
OPC <- XXO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################
IN <- XXO_P301S_IN
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
saveRDS(IN, file = 'FCG_YoungCPZ_IN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_YoungCPZ_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_IN_markers, "FCG_YoungCPZ_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################
EN <- XXO_P301S_EN
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
saveRDS(EN, file = 'FCG_YoungCPZ_EN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_YoungCPZ_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_EN_markers, "FCG_YoungCPZ_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subXXO_P301S_cell_counts_res0.3.csv")
#######################################################################

XYO_P301S_EN <- subset(Cluster_EN, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_IN <- subset(Cluster_IN, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_MG <- subset(Cluster_MG, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_AST <- subset(Cluster_AST, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_OL <- subset(Cluster_OL, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_OPC <- subset(Cluster_OPC, idents=c("XYO_NTG_E4","XYO_P301S_E4"))
XYO_P301S_VC <- subset(Cluster_VC, idents=c("XYO_NTG_E4","XYO_P301S_E4"))

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/subclustering_XYO_P301S")
MG <- XYO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.15.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.15.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYO_P301S_cell_counts_res0.15.csv")
#######################################################################
MG <- XYO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYO_P301S_cell_counts_res0.2.csv")
#######################################################################
MG <- XYO_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################
AST <- XYO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYO_P301S_cell_counts.csv")
##################################################################
AST <- XYO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYO_P301S_cell_counts_res0.2.csv")
#######################################################################
AST <- XYO_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################
OL <- XYO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYO_P301S_cell_counts.csv")
#######################################################################
OL <- XYO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYO_P301S_cell_counts_res0.2.csv")
#######################################################################
OL <- XYO_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################
OPC <- XYO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYO_P301S_cell_counts.csv")
#######################################################################
OPC <- XYO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYO_P301S_cell_counts_res0.2.csv")
#######################################################################
OPC <- XYO_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################
IN <- XYO_P301S_IN
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
saveRDS(IN, file = 'FCG_YoungCPZ_IN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_YoungCPZ_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_IN_markers, "FCG_YoungCPZ_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################
EN <- XYO_P301S_EN
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
saveRDS(EN, file = 'FCG_YoungCPZ_EN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_YoungCPZ_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_EN_markers, "FCG_YoungCPZ_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subXYO_P301S_cell_counts_res0.3.csv")
#######################################################################

XXT_P301S_EN <- subset(Cluster_EN, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_IN <- subset(Cluster_IN, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_MG <- subset(Cluster_MG, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_AST <- subset(Cluster_AST, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_OL <- subset(Cluster_OL, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_OPC <- subset(Cluster_OPC, idents=c("XXT_NTG_E4","XXT_P301S_E4"))
XXT_P301S_VC <- subset(Cluster_VC, idents=c("XXT_NTG_E4","XXT_P301S_E4"))

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/subclustering_XXT_P301S")
MG <- XXT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.15.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.15.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXT_P301S_cell_counts_res0.15.csv")
#######################################################################
MG <- XXT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXT_P301S_cell_counts_res0.2.csv")
#######################################################################
MG <- XXT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################
AST <- XXT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXT_P301S_cell_counts.csv")
##################################################################
AST <- XXT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXT_P301S_cell_counts_res0.2.csv")
#######################################################################
AST <- XXT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################
OL <- XXT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXT_P301S_cell_counts.csv")
#######################################################################
OL <- XXT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXT_P301S_cell_counts_res0.2.csv")
#######################################################################
OL <- XXT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################
OPC <- XXT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXT_P301S_cell_counts.csv")
#######################################################################
OPC <- XXT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXT_P301S_cell_counts_res0.2.csv")
#######################################################################
OPC <- XXT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################
IN <- XXT_P301S_IN
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
saveRDS(IN, file = 'FCG_YoungCPZ_IN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_YoungCPZ_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_IN_markers, "FCG_YoungCPZ_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################
EN <- XXT_P301S_EN
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
saveRDS(EN, file = 'FCG_YoungCPZ_EN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_YoungCPZ_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_EN_markers, "FCG_YoungCPZ_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subXXT_P301S_cell_counts_res0.3.csv")
#######################################################################

XYT_P301S_EN <- subset(Cluster_EN, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_IN <- subset(Cluster_IN, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_MG <- subset(Cluster_MG, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_AST <- subset(Cluster_AST, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_OL <- subset(Cluster_OL, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_OPC <- subset(Cluster_OPC, idents=c("XYT_NTG_E4","XYT_P301S_E4"))
XYT_P301S_VC <- subset(Cluster_VC, idents=c("XYT_NTG_E4","XYT_P301S_E4"))

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration/subclustering_XYT_P301S")
MG <- XYT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.15.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.15.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYT_P301S_cell_counts_res0.15.csv")
#######################################################################
MG <- XYT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYT_P301S_cell_counts_res0.2.csv")
#######################################################################
MG <- XYT_P301S_MG
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
saveRDS(MG, file = 'FCG_YoungCPZ_MG_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_MG_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_YoungCPZ_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_MG_markers, "FCG_YoungCPZ_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################
AST <- XYT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYT_P301S_cell_counts.csv")
##################################################################
AST <- XYT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYT_P301S_cell_counts_res0.2.csv")
#######################################################################
AST <- XYT_P301S_AST
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
saveRDS(AST, file = 'FCG_YoungCPZ_AST_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_AST_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_YoungCPZ_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_AST_markers, "FCG_YoungCPZ_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################
OL <- XYT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYT_P301S_cell_counts.csv")
#######################################################################
OL <- XYT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYT_P301S_cell_counts_res0.2.csv")
#######################################################################
OL <- XYT_P301S_OL
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
saveRDS(OL, file = 'FCG_YoungCPZ_OL_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OL_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_YoungCPZ_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OL_markers, "FCG_YoungCPZ_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################
OPC <- XYT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.15.rds')
pdf("FCG_YoungCPZ_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYT_P301S_cell_counts.csv")
#######################################################################
OPC <- XYT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.2.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.2.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYT_P301S_cell_counts_res0.2.csv")
#######################################################################
OPC <- XYT_P301S_OPC
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
saveRDS(OPC, file = 'FCG_YoungCPZ_OPC_reclusted_res0.3.rds')
pdf("FCG_YoungCPZ_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_OPC_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_YoungCPZ_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_OPC_markers, "FCG_YoungCPZ_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################
IN <- XYT_P301S_IN
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
saveRDS(IN, file = 'FCG_YoungCPZ_IN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_IN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_YoungCPZ_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_IN_markers, "FCG_YoungCPZ_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################
EN <- XYT_P301S_EN
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
saveRDS(EN, file = 'FCG_YoungCPZ_EN_reclusted_res0.1.rds')
pdf("FCG_YoungCPZ_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("FCG_YoungCPZ_EN_umap_Sample_res0.3.pdf", width=11, height=2)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_YoungCPZ_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_YoungCPZ_EN_markers, "FCG_YoungCPZ_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subXYT_P301S_cell_counts_res0.3.csv")
#######################################################################


