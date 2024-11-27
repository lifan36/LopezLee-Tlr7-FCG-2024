
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/FCG_5week/DF_2ndRound")
XXO_1 <- readRDS(file = "XXO_1_singlets_PCA.rds")
XXO_2 <- readRDS(file = "XXO_2_singlets_PCA.rds")
XXO_3 <- readRDS(file = "XXO_3_singlets_PCA.rds")
XXO_Cup_1 <- readRDS(file = "XXO_Cup_1_singlets_PCA.rds")
XXO_Cup_2 <- readRDS(file = "XXO_Cup_2_singlets_PCA.rds")
XXO_Cup_3 <- readRDS(file = "XXO_Cup_3_singlets_PCA.rds")

XYO_1 <- readRDS(file = "XYO_1_singlets_PCA.rds")
XYO_2 <- readRDS(file = "XYO_2_singlets_PCA.rds")
XYO_3 <- readRDS(file = "XYO_3_singlets_PCA.rds")
XYO_Cup_1 <- readRDS(file = "XYO_Cup_1_singlets_PCA.rds")
XYO_Cup_2 <- readRDS(file = "XYO_Cup_2_singlets_PCA.rds")
XYO_Cup_3 <- readRDS(file = "XYO_Cup_3_singlets_PCA.rds")

XXT_1 <- readRDS(file = "XXT_1_singlets_PCA.rds")
XXT_2 <- readRDS(file = "XXT_2_singlets_PCA.rds")
XXT_3 <- readRDS(file = "XXT_3_singlets_PCA.rds")
XXT_Cup_1 <- readRDS(file = "XXT_Cup_1_singlets_PCA.rds")
XXT_Cup_2 <- readRDS(file = "XXT_Cup_2_singlets_PCA.rds")
XXT_Cup_3 <- readRDS(file = "XXT_Cup_3_singlets_PCA.rds")

XYT_1 <- readRDS(file = "XYT_1_singlets_PCA.rds")
XYT_2 <- readRDS(file = "XYT_2_singlets_PCA.rds")
XYT_3 <- readRDS(file = "XYT_3_singlets_PCA.rds")
XYT_Cup_1 <- readRDS(file = "XYT_Cup_1_singlets_PCA.rds")
XYT_Cup_2 <- readRDS(file = "XYT_Cup_2_singlets_PCA.rds")
XYT_Cup_3 <- readRDS(file = "XYT_Cup_3_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/FCG_5week/integration")
XXO <- c(XXO_1, XXO_2, XXO_3)
anchors_XXO <- FindIntegrationAnchors(object.list = XXO, dims = 1:30)
XXO_integrated <- IntegrateData(anchorset = anchors_XXO, dims = 1:30)
rm(XXO_1, XXO_2, XXO_3, XXO)
XXO_Cup <- c(XXO_Cup_1, XXO_Cup_2, XXO_Cup_3)
anchors_XXO_Cup <- FindIntegrationAnchors(object.list = XXO_Cup, dims = 1:30)
XXO_Cup_integrated <- IntegrateData(anchorset = anchors_XXO_Cup, dims = 1:30)
rm(XXO_Cup_1, XXO_Cup_2, XXO_Cup_3, XXO_Cup)

XYO <- c(XYO_1, XYO_2, XYO_3)
anchors_XYO <- FindIntegrationAnchors(object.list = XYO, dims = 1:30)
XYO_integrated <- IntegrateData(anchorset = anchors_XYO, dims = 1:30)
rm(XYO_1, XYO_2, XYO_3, XYO)
XYO_Cup <- c(XYO_Cup_1, XYO_Cup_2, XYO_Cup_3)
anchors_XYO_Cup <- FindIntegrationAnchors(object.list = XYO_Cup, dims = 1:30)
XYO_Cup_integrated <- IntegrateData(anchorset = anchors_XYO_Cup, dims = 1:30)
rm(XYO_Cup_1, XYO_Cup_2, XYO_Cup_3, XYO_Cup)

XXT <- c(XXT_1, XXT_2, XXT_3)
anchors_XXT <- FindIntegrationAnchors(object.list = XXT, dims = 1:30)
XXT_integrated <- IntegrateData(anchorset = anchors_XXT, dims = 1:30)
rm(XXT_1, XXT_2, XXT_3, XXT)
XXT_Cup <- c(XXT_Cup_1, XXT_Cup_2, XXT_Cup_3)
anchors_XXT_Cup <- FindIntegrationAnchors(object.list = XXT_Cup, dims = 1:30)
XXT_Cup_integrated <- IntegrateData(anchorset = anchors_XXT_Cup, dims = 1:30)
rm(XXT_Cup_1, XXT_Cup_2, XXT_Cup_3, XXT_Cup)

XYT <- c(XYT_1, XYT_2, XYT_3)
anchors_XYT <- FindIntegrationAnchors(object.list = XYT, dims = 1:30)
XYT_integrated <- IntegrateData(anchorset = anchors_XYT, dims = 1:30)
rm(XYT_1, XYT_2, XYT_3, XYT)
XYT_Cup <- c(XYT_Cup_1, XYT_Cup_2, XYT_Cup_3)
anchors_XYT_Cup <- FindIntegrationAnchors(object.list = XYT_Cup, dims = 1:30)
XYT_Cup_integrated <- IntegrateData(anchorset = anchors_XYT_Cup, dims = 1:30)
rm(XYT_Cup_1, XYT_Cup_2, XYT_Cup_3, XYT_Cup)


FCG <- c(XXO_integrated, XXO_Cup_integrated, XYO_integrated, XYO_Cup_integrated, XXT_integrated, XXT_Cup_integrated, XYT_integrated, XYT_Cup_integrated)
anchors_FCG <- FindIntegrationAnchors(object.list = FCG, dims = 1:30)
FCG_integrated <- IntegrateData(anchorset = anchors_FCG, dims = 1:30)
rm(XXO_integrated, XXO_Cup_integrated, XYO_integrated, XYO_Cup_integrated, XXT_integrated, XXT_Cup_integrated, XYT_integrated, XYT_Cup_integrated, FCG)

#saveRDS(FCG_integrated, file = "FCG_integrated.rds")

DefaultAssay(FCG_integrated) <- 'integrated'

# FCG_integrated <- NormalizeData(FCG_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# FCG_integrated <- FindVariableFeatures(FCG_integrated, selection.method = "vst", nfeatures = 3000)

FCG_integrated <- ScaleData(FCG_integrated, verbose = FALSE)
FCG_integrated <- RunPCA(FCG_integrated, features = VariableFeatures(object = FCG_integrated), verbose = FALSE)

FCG_integrated <- FindNeighbors(FCG_integrated, dims = 1:15)
FCG_integrated <- FindClusters(FCG_integrated, resolution = 0.1)
FCG_integrated <- RunUMAP(FCG_integrated, dims = 1: 15)

DefaultAssay(FCG_integrated) <- 'RNA'
FCG_integrated <- NormalizeData(FCG_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
FCG_integrated <- ScaleData(FCG_integrated, features = rownames(FCG_integrated))

#saveRDS(FCG_integrated, file = 'FCG_integrated_PCA_0.1.rds')
#FCG_integrated <- readRDS(file = "FCG_integrated_PCA_0.1.rds")

FCG_integrated$Condition <- factor(x = FCG_integrated$Condition, levels = c("XXO","XXO_Cup","XYO","XYO_Cup","XXT","XXT_Cup","XYT","XYT_Cup"))
FCG_integrated$Sample_Name <- factor(x = FCG_integrated$Sample_Name, levels = c("XXO_1","XXO_2","XXO_3","XXO_Cup_1","XXO_Cup_2","XXO_Cup_3",
                                                                                "XYO_1","XYO_2","XYO_3","XYO_Cup_1","XYO_Cup_2","XYO_Cup_3",
                                                                                "XXT_1","XXT_2","XXT_3","XXT_Cup_1","XXT_Cup_2","XXT_Cup_3",
                                                                                "XYT_1","XYT_2","XYT_3","XYT_Cup_1","XYT_Cup_2","XYT_Cup_3"
                                                                                ))

pdf("FCG_QC.pdf", width=9, height=4)
Idents(FCG_integrated) <- "Condition"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

pdf("FCG_QC_Sample.pdf", width=18, height=4)
Idents(FCG_integrated) <- "Sample_Name"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(FCG_integrated) <- "seurat_clusters"
pdf("FCG_integrated_umap.pdf", width=5, height=4)
DimPlot(FCG_integrated, reduction = 'umap', label = T)
dev.off()
pdf("FCG_integrated_umap_split_individual.pdf", width=12, height=7)
DimPlot(FCG_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
pdf("FCG_integrated_umap_split_Condition.pdf", width=10, height=4.5)
DimPlot(FCG_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()

write.csv(table(FCG_integrated$seurat_clusters, FCG_integrated$Sample_Name), "FCG_cell_counts_cluster_by_sample.csv")

saveRDS(FCG_integrated, file = 'FCG_integrated_PCA_0.1.rds')

DefaultAssay(FCG_integrated) <- 'RNA'

FCG_markers <- FindAllMarkers(FCG_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(FCG_markers, "FCG_markers.csv")

FCG_markers <- read.csv(file = "FCG_markers.csv", header=T,row.names =1)
top5 <- FCG_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("FCG_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(FCG_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("FCG_annotation_combine.pdf", width=10, height=5)
DotPlot(object = FCG_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

