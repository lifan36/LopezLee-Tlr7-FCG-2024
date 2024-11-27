
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
setwd("/athena/ganlab/scratch/lif4001/LG93E4/DF_2ndRound")

XXO_NTG_E4_1 <- readRDS(file = "XXO_NTG_E4_1_singlets_PCA.rds")
XXO_NTG_E4_2 <- readRDS(file = "XXO_NTG_E4_2_singlets_PCA.rds")
XXO_NTG_E4_3 <- readRDS(file = "XXO_NTG_E4_3_singlets_PCA.rds")

XXO_P301S_E4_1 <- readRDS(file = "XXO_P301S_E4_1_singlets_PCA.rds")
XXO_P301S_E4_2 <- readRDS(file = "XXO_P301S_E4_2_singlets_PCA.rds")
XXO_P301S_E4_3 <- readRDS(file = "XXO_P301S_E4_3_singlets_PCA.rds")
#XXO_P301S_E4_4 <- readRDS(file = "XXO_P301S_E4_4_singlets_PCA.rds")

XYO_NTG_E4_1 <- readRDS(file = "XYO_NTG_E4_1_singlets_PCA.rds")
XYO_NTG_E4_2 <- readRDS(file = "XYO_NTG_E4_2_singlets_PCA.rds")
XYO_NTG_E4_3 <- readRDS(file = "XYO_NTG_E4_3_singlets_PCA.rds")

#XYO_P301S_E4_1 <- readRDS(file = "XYO_P301S_E4_1_singlets_PCA.rds")
XYO_P301S_E4_2 <- readRDS(file = "XYO_P301S_E4_2_singlets_PCA.rds")
XYO_P301S_E4_3 <- readRDS(file = "XYO_P301S_E4_3_singlets_PCA.rds")
XYO_P301S_E4_4 <- readRDS(file = "XYO_P301S_E4_4_singlets_PCA.rds")

XXT_NTG_E4_1 <- readRDS(file = "XXT_NTG_E4_1_singlets_PCA.rds")
XXT_NTG_E4_2 <- readRDS(file = "XXT_NTG_E4_2_singlets_PCA.rds")
XXT_NTG_E4_3 <- readRDS(file = "XXT_NTG_E4_3_singlets_PCA.rds")

XXT_P301S_E4_1 <- readRDS(file = "XXT_P301S_E4_1_singlets_PCA.rds")
XXT_P301S_E4_2 <- readRDS(file = "XXT_P301S_E4_2_singlets_PCA.rds")
#XXT_P301S_E4_3 <- readRDS(file = "XXT_P301S_E4_3_singlets_PCA.rds")
XXT_P301S_E4_4 <- readRDS(file = "XXT_P301S_E4_4_singlets_PCA.rds")

XYT_NTG_E4_1 <- readRDS(file = "XYT_NTG_E4_1_singlets_PCA.rds")
XYT_NTG_E4_2 <- readRDS(file = "XYT_NTG_E4_2_singlets_PCA.rds")
XYT_NTG_E4_3 <- readRDS(file = "XYT_NTG_E4_3_singlets_PCA.rds")

XYT_P301S_E4_1 <- readRDS(file = "XYT_P301S_E4_1_singlets_PCA.rds")
XYT_P301S_E4_2 <- readRDS(file = "XYT_P301S_E4_2_singlets_PCA.rds")
#XYT_P301S_E4_3 <- readRDS(file = "XYT_P301S_E4_3_singlets_PCA.rds")
XYT_P301S_E4_4 <- readRDS(file = "XYT_P301S_E4_4_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/LG93E4/integration_without_outlier")

XXO_NTG_E4 <- c(XXO_NTG_E4_1, XXO_NTG_E4_2, XXO_NTG_E4_3)
anchors_XXO_NTG_E4 <- FindIntegrationAnchors(object.list = XXO_NTG_E4, dims = 1:30)
XXO_NTG_E4_integrated <- IntegrateData(anchorset = anchors_XXO_NTG_E4, dims = 1:30)
rm(XXO_NTG_E4_1, XXO_NTG_E4_2, XXO_NTG_E4_3, XXO_NTG_E4)

XXO_P301S_E4 <- c(XXO_P301S_E4_1, XXO_P301S_E4_2, XXO_P301S_E4_3)
anchors_XXO_P301S_E4 <- FindIntegrationAnchors(object.list = XXO_P301S_E4, dims = 1:30)
XXO_P301S_E4_integrated <- IntegrateData(anchorset = anchors_XXO_P301S_E4, dims = 1:30)
rm(XXO_P301S_E4_1, XXO_P301S_E4_2, XXO_P301S_E4_3, XXO_P301S_E4)

XYO_NTG_E4 <- c(XYO_NTG_E4_1, XYO_NTG_E4_2, XYO_NTG_E4_3)
anchors_XYO_NTG_E4 <- FindIntegrationAnchors(object.list = XYO_NTG_E4, dims = 1:30)
XYO_NTG_E4_integrated <- IntegrateData(anchorset = anchors_XYO_NTG_E4, dims = 1:30)
rm(XYO_NTG_E4_1, XYO_NTG_E4_2, XYO_NTG_E4_3, XYO_NTG_E4)

XYO_P301S_E4 <- c(XYO_P301S_E4_2, XYO_P301S_E4_3, XYO_P301S_E4_4)
anchors_XYO_P301S_E4 <- FindIntegrationAnchors(object.list = XYO_P301S_E4, dims = 1:30)
XYO_P301S_E4_integrated <- IntegrateData(anchorset = anchors_XYO_P301S_E4, dims = 1:30)
rm(XYO_P301S_E4_2, XYO_P301S_E4_3, XYO_P301S_E4, XYO_P301S_E4_4)

XXT_NTG_E4 <- c(XXT_NTG_E4_1, XXT_NTG_E4_2, XXT_NTG_E4_3)
anchors_XXT_NTG_E4 <- FindIntegrationAnchors(object.list = XXT_NTG_E4, dims = 1:30)
XXT_NTG_E4_integrated <- IntegrateData(anchorset = anchors_XXT_NTG_E4, dims = 1:30)
rm(XXT_NTG_E4_1, XXT_NTG_E4_2, XXT_NTG_E4_3, XXT_NTG_E4)

XXT_P301S_E4 <- c(XXT_P301S_E4_1, XXT_P301S_E4_2, XXT_P301S_E4_4)
anchors_XXT_P301S_E4 <- FindIntegrationAnchors(object.list = XXT_P301S_E4, dims = 1:30)
XXT_P301S_E4_integrated <- IntegrateData(anchorset = anchors_XXT_P301S_E4, dims = 1:30)
rm(XXT_P301S_E4_1, XXT_P301S_E4_2, XXT_P301S_E4, XXT_P301S_E4_4)

XYT_NTG_E4 <- c(XYT_NTG_E4_1, XYT_NTG_E4_2, XYT_NTG_E4_3)
anchors_XYT_NTG_E4 <- FindIntegrationAnchors(object.list = XYT_NTG_E4, dims = 1:30)
XYT_NTG_E4_integrated <- IntegrateData(anchorset = anchors_XYT_NTG_E4, dims = 1:30)
rm(XYT_NTG_E4_1, XYT_NTG_E4_2, XYT_NTG_E4_3, XYT_NTG_E4)

XYT_P301S_E4 <- c(XYT_P301S_E4_1, XYT_P301S_E4_2, XYT_P301S_E4_4)
anchors_XYT_P301S_E4 <- FindIntegrationAnchors(object.list = XYT_P301S_E4, dims = 1:30)
XYT_P301S_E4_integrated <- IntegrateData(anchorset = anchors_XYT_P301S_E4, dims = 1:30)
rm(XYT_P301S_E4_1, XYT_P301S_E4_2, XYT_P301S_E4, XYT_P301S_E4_4)

LG93E4 <- c(XXO_NTG_E4_integrated,XXO_P301S_E4_integrated,XYO_NTG_E4_integrated,XYO_P301S_E4_integrated,XXT_NTG_E4_integrated,XXT_P301S_E4_integrated,XYT_NTG_E4_integrated,XYT_P301S_E4_integrated)
anchors_LG93E4 <- FindIntegrationAnchors(object.list = LG93E4, dims = 1:30)
LG93E4_integrated <- IntegrateData(anchorset = anchors_LG93E4, dims = 1:30)
rm(XXO_NTG_E4_integrated,XXO_P301S_E4_integrated,XYO_NTG_E4_integrated,XYO_P301S_E4_integrated,XXT_NTG_E4_integrated,XXT_P301S_E4_integrated,XYT_NTG_E4_integrated,XYT_P301S_E4_integrated, LG93E4)

#saveRDS(LG93E4_integrated, file = "LG93E4_integrated.rds")

DefaultAssay(LG93E4_integrated) <- 'integrated'

# LG93E4_integrated <- NormalizeData(LG93E4_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG93E4_integrated <- FindVariableFeatures(LG93E4_integrated, selection.method = "vst", nfeatures = 3000)

LG93E4_integrated <- ScaleData(LG93E4_integrated, verbose = FALSE)
LG93E4_integrated <- RunPCA(LG93E4_integrated, features = VariableFeatures(object = LG93E4_integrated), verbose = FALSE)

LG93E4_integrated <- FindNeighbors(LG93E4_integrated, dims = 1:15)
LG93E4_integrated <- FindClusters(LG93E4_integrated, resolution = 0.1)
LG93E4_integrated <- RunUMAP(LG93E4_integrated, dims = 1: 15)

DefaultAssay(LG93E4_integrated) <- 'RNA'
LG93E4_integrated <- NormalizeData(LG93E4_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG93E4_integrated <- ScaleData(LG93E4_integrated, features = rownames(LG93E4_integrated))

#saveRDS(LG93E4_integrated, file = 'LG93E4_integrated_PCA_0.1.rds')
#LG93E4_integrated <- readRDS(file = "LG93E4_integrated_PCA_0.1.rds")

LG93E4_integrated$Condition <- factor(x = LG93E4_integrated$Condition, levels = c("XXO_NTG_E4","XXO_P301S_E4","XYO_NTG_E4","XYO_P301S_E4","XXT_NTG_E4","XXT_P301S_E4","XYT_NTG_E4","XYT_P301S_E4"))
LG93E4_integrated$Sample_Name <- factor(x = LG93E4_integrated$Sample_Name, levels = c("XXO_NTG_E4_1","XXO_NTG_E4_2","XXO_NTG_E4_3",
                                                                                      "XXO_P301S_E4_1","XXO_P301S_E4_2","XXO_P301S_E4_3",
                                                                                      "XYO_NTG_E4_1","XYO_NTG_E4_2","XYO_NTG_E4_3",
                                                                                      "XYO_P301S_E4_2","XYO_P301S_E4_3","XYO_P301S_E4_4",
                                                                                      "XXT_NTG_E4_1","XXT_NTG_E4_2","XXT_NTG_E4_3",
                                                                                      "XXT_P301S_E4_1","XXT_P301S_E4_2","XXT_P301S_E4_4",
                                                                                      "XYT_NTG_E4_1","XYT_NTG_E4_2","XYT_NTG_E4_3",
                                                                                      "XYT_P301S_E4_1","XYT_P301S_E4_2","XYT_P301S_E4_4"))

pdf("LG93E4_QC.pdf", width=9, height=4)
Idents(LG93E4_integrated) <- "Condition"
VlnPlot(object = LG93E4_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(LG93E4_integrated) <- "Sample_Name"
pdf("LG93E4_QC_Sample.pdf", width=12, height=4)

VlnPlot(object = LG93E4_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(LG93E4_integrated) <- "seurat_clusters"
pdf("LG93E4_integrated_umap.pdf", width=5, height=4)
DimPlot(LG93E4_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG93E4_integrated_umap_split_individual.pdf", width=18, height=12)
DimPlot(LG93E4_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
pdf("LG93E4_integrated_umap_split_Condition.pdf", width=10, height=6)
DimPlot(LG93E4_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()

write.csv(table(LG93E4_integrated$seurat_clusters, LG93E4_integrated$Sample_Name), "LG93E4_cell_counts_cluster_by_sample.csv")

DefaultAssay(LG93E4_integrated) <- 'RNA'

LG93E4_markers <- FindAllMarkers(LG93E4_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(LG93E4_markers, "LG93E4_markers.csv")


saveRDS(LG93E4_integrated, file = 'LG93E4_integrated_PCA_0.1.rds')

LG93E4_markers <- read.csv(file = "LG93E4_markers.csv", header=T,row.names =1)
top5 <- LG93E4_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("LG93E4_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(LG93E4_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG93E4_annotation_combine.pdf", width=10, height=5)
DotPlot(object = LG93E4_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

