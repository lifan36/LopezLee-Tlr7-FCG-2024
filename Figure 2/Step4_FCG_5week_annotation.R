
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration")
FCG_integrated <- readRDS("FCG_integrated_PCA_0.1.rds")
#Remove cluster 17-duplicates, 12-Reln+ CR cells, 5-unknown, mixed
FCG_integrated <- subset(FCG_integrated, idents=c("12"), invert=T)
Idents(FCG_integrated) <- "seurat_clusters"
FCG_integrated <- RenameIdents(FCG_integrated,
                                 `0` = "excitatory neurons", `1`="oligodendrocytes", `2`="astrocytes", `3`="microglia",
                                 `4`="excitatory neurons", `5`="excitatory neurons", `6`="inhibitory neurons", `7`="OPCs",
                                 `8`="excitatory neurons", `9`="vascular cells", `10`="excitatory neurons", `11`="CHOR",
                                `13`="vascular cells", `14`="vascular cells"
)

pdf("FCG_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(FCG_integrated, reduction = 'umap', label = F)
dev.off()

FCG_integrated$celltype.orig.ident <- paste(Idents(FCG_integrated), FCG_integrated$orig.ident, sep = "_")
FCG_integrated$celltype <- Idents(FCG_integrated)


FCG_integrated$Condition <- factor(x = FCG_integrated$Condition, levels = c("XXO","XXO_Cup", "XYO", "XYO_Cup","XXT","XXT_Cup", "XYT", "XYT_Cup"))
FCG_integrated$Sample_Name <- factor(x = FCG_integrated$Sample_Name, levels = c("XXO_1","XXO_2","XXO_3","XXO_Cup_1", "XXO_Cup_2", "XXO_Cup_3",
                                                                                      "XYO_1","XYO_2","XYO_3","XYO_Cup_1", "XYO_Cup_2", "XYO_Cup_3",
                                                                                      "XXT_1","XXT_2","XXT_3","XXT_Cup_1", "XXT_Cup_2", "XXT_Cup_3",
                                                                                      "XYT_1","XYT_2","XYT_3","XYT_Cup_1", "XYT_Cup_2", "XYT_Cup_3"))

saveRDS(FCG_integrated, file = "FCG_integrated_Annotation.rds")

pdf("FCG_QC.pdf", width=9, height=4)
Idents(FCG_integrated) <- "Condition"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


data <- FCG_integrated
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

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/FCG_5week/integration",
       width=4,height=4,units="in")

#markers for annotation
pdf("annotation_1.pdf", width=10.5, height=2.7)
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

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/FCG_5week/integration",
       width=6,height=4,units="in")

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
Cluster_CHOR <- subset(FCG_integrated, idents = "CHOR")

saveRDS(Cluster_EN, file = "FCG_EN_subset.rds")
saveRDS(Cluster_IN, file = "FCG_IN_subset.rds")
saveRDS(Cluster_MG, file = "FCG_MG_subset.rds")
saveRDS(Cluster_AST, file = "FCG_AST_subset.rds")
saveRDS(Cluster_OL, file = "FCG_OL_subset.rds")
saveRDS(Cluster_OPC, file = "FCG_OPC_subset.rds")
saveRDS(Cluster_VC, file = "FCG_VC_subset.rds")
saveRDS(Cluster_CHOR, file = "FCG_CHOR_subset.rds")

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(MG, file = 'FCG_MG_reclusted_res0.15.rds')
pdf("FCG_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_MG_umap_Condition_res0.15.pdf", width=11, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_MG_umap_Sample_res0.15.pdf", width=13, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_MG_markers, "FCG_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.15.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(MG, file = 'FCG_MG_reclusted_res0.2.rds')
pdf("FCG_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_MG_umap_Condition_res0.2.pdf", width=11, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_MG_umap_Sample_res0.2.pdf", width=13, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_MG_markers, "FCG_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(MG, file = 'FCG_MG_reclusted_res0.3.rds')
pdf("FCG_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("FCG_MG_umap_Condition_res0.3.pdf", width=11, height=5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_MG_umap_Sample_res0.3.pdf", width=13, height=8)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(MG) <- 'RNA'
FCG_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_MG_markers, "FCG_MG_markers_res0.3.csv")
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
saveRDS(AST, file = 'FCG_AST_reclusted_res0.15.rds')
pdf("FCG_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_AST_umap_Condition.pdf", width=11, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_AST_umap_Sample.pdf", width=13, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_AST_markers, "FCG_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts.csv")
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(AST, file = 'FCG_AST_reclusted_res0.2.rds')
pdf("FCG_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_AST_umap_Condition_res0.2.pdf", width=11, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_AST_umap_Sample_res0.2.pdf", width=13, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_AST_markers, "FCG_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(AST, file = 'FCG_AST_reclusted_res0.3.rds')
pdf("FCG_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("FCG_AST_umap_Condition_res0.3.pdf", width=11, height=5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_AST_umap_Sample_res0.3.pdf", width=13, height=8)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(AST) <- 'RNA'
FCG_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_AST_markers, "FCG_AST_markers_res0.3.csv")
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
saveRDS(OL, file = 'FCG_OL_reclusted_res0.15.rds')
pdf("FCG_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OL_umap_Condition.pdf", width=11, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OL_umap_Sample.pdf", width=13, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OL_markers, "FCG_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(OL, file = 'FCG_OL_reclusted_res0.2.rds')
pdf("FCG_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OL_umap_Condition_res0.2.pdf", width=11, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OL_umap_Sample_res0.2.pdf", width=13, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OL_markers, "FCG_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(OL, file = 'FCG_OL_reclusted_res0.3.rds')
pdf("FCG_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OL_umap_Condition_res0.3.pdf", width=11, height=5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OL_umap_Sample_res0.3.pdf", width=13, height=8)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OL) <- 'RNA'
FCG_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OL_markers, "FCG_OL_markers_res0.3.csv")
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
saveRDS(OPC, file = 'FCG_OPC_reclusted_res0.15.rds')
pdf("FCG_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OPC_umap_Condition.pdf", width=11, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OPC_umap_Sample.pdf", width=13, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OPC_markers, "FCG_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(OPC, file = 'FCG_OPC_reclusted_res0.2.rds')
pdf("FCG_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OPC_umap_Condition_res0.2.pdf", width=11, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OPC_umap_Sample_res0.2.pdf", width=13, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OPC_markers, "FCG_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(OPC, file = 'FCG_OPC_reclusted_res0.3.rds')
pdf("FCG_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("FCG_OPC_umap_Condition_res0.3.pdf", width=11, height=5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_OPC_umap_Sample_res0.3.pdf", width=13, height=8)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(OPC) <- 'RNA'
FCG_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_OPC_markers, "FCG_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.3.csv")
#######################################################################
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(EN, file = 'FCG_EN_reclusted_res0.1.rds')
pdf("FCG_EN_umap_res0.1.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_EN_umap_Condition_res0.1.pdf", width=11, height=5)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_EN_umap_Sample_res0.1.pdf", width=13, height=8)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(EN) <- 'RNA'
FCG_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_EN_markers, "FCG_EN_markers_res0.1.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subcluster_cell_counts_res0.1.csv")
#######################################################################
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(IN, file = 'FCG_IN_reclusted_res0.1.rds')
pdf("FCG_IN_umap_res0.1.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("FCG_IN_umap_Condition_res0.1.pdf", width=11, height=5)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_IN_umap_Sample_res0.1.pdf", width=13, height=8)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(IN) <- 'RNA'
FCG_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_IN_markers, "FCG_IN_markers_res0.1.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subcluster_cell_counts_res0.1.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(CHOR, file = 'FCG_CHOR_reclusted_res0.15.rds')
pdf("FCG_CHOR_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("FCG_CHOR_umap_Condition_res0.15.pdf", width=11, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_CHOR_umap_Sample_res0.15.pdf", width=13, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
FCG_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_CHOR_markers, "FCG_CHOR_markers_res0.15.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "CHOR_subcluster_cell_counts_res0.15.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(CHOR, file = 'FCG_CHOR_reclusted_res0.2.rds')
pdf("FCG_CHOR_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("FCG_CHOR_umap_Condition_res0.2.pdf", width=11, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_CHOR_umap_Sample_res0.2.pdf", width=13, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
FCG_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_CHOR_markers, "FCG_CHOR_markers_res0.2.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "CHOR_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration/subclustering")
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
saveRDS(CHOR, file = 'FCG_CHOR_reclusted_res0.3.rds')
pdf("FCG_CHOR_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(CHOR, reduction = 'umap', label = T)
dev.off()
pdf("FCG_CHOR_umap_Condition_res0.3.pdf", width=11, height=5)
DimPlot(CHOR, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()
pdf("FCG_CHOR_umap_Sample_res0.3.pdf", width=13, height=8)
DimPlot(CHOR, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
DefaultAssay(CHOR) <- 'RNA'
FCG_CHOR_markers <- FindAllMarkers(CHOR, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(FCG_CHOR_markers, "FCG_CHOR_markers_res0.3.csv")
write.csv(table(CHOR$seurat_clusters, CHOR$Sample_Name), "CHOR_subcluster_cell_counts_res0.3.csv")
#######################################################################
