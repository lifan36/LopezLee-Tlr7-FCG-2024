
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
setwd("/athena/ganlab/scratch/lif4001/FCG/DF_2ndRound")
XXO_1 <- readRDS(file = "XXO_1_singlets_PCA.rds")
XXO_2 <- readRDS(file = "XXO_2_singlets_PCA.rds")
XXO_3 <- readRDS(file = "XXO_3_singlets_PCA.rds")
XXO_4 <- readRDS(file = "XXO_4_singlets_PCA.rds")
XXO_5 <- readRDS(file = "XXO_5_singlets_PCA.rds")
XXO_6 <- readRDS(file = "XXO_6_singlets_PCA.rds")
XYO_1 <- readRDS(file = "XYO_1_singlets_PCA.rds")
XYO_2 <- readRDS(file = "XYO_2_singlets_PCA.rds")
XYO_3 <- readRDS(file = "XYO_3_singlets_PCA.rds")
XYO_4 <- readRDS(file = "XYO_4_singlets_PCA.rds")
XYO_5 <- readRDS(file = "XYO_5_singlets_PCA.rds")
XYO_6 <- readRDS(file = "XYO_6_singlets_PCA.rds")
XXT_1 <- readRDS(file = "XXT_1_singlets_PCA.rds")
XXT_2 <- readRDS(file = "XXT_2_singlets_PCA.rds")
XXT_3 <- readRDS(file = "XXT_3_singlets_PCA.rds")
XXT_4 <- readRDS(file = "XXT_4_singlets_PCA.rds")
XXT_5 <- readRDS(file = "XXT_5_singlets_PCA.rds")

XYT_1 <- readRDS(file = "XYT_1_singlets_PCA.rds")
XYT_2 <- readRDS(file = "XYT_2_singlets_PCA.rds")
XYT_3 <- readRDS(file = "XYT_3_singlets_PCA.rds")
XYT_4 <- readRDS(file = "XYT_4_singlets_PCA.rds")
XYT_5 <- readRDS(file = "XYT_5_singlets_PCA.rds")
XYT_6 <- readRDS(file = "XYT_6_singlets_PCA.rds")



setwd("/athena/ganlab/scratch/lif4001/FCG/integration")
XXOC <- c(XXO_1, XXO_5, XXO_6)
anchors_XXOC <- FindIntegrationAnchors(object.list = XXOC, dims = 1:30)
XXOC_integrated <- IntegrateData(anchorset = anchors_XXOC, dims = 1:30)
rm(XXO_1, XXO_5, XXO_6, XXOC)

XXOD <- c(XXO_2, XXO_3, XXO_4)
anchors_XXOD <- FindIntegrationAnchors(object.list = XXOD, dims = 1:30)
XXOD_integrated <- IntegrateData(anchorset = anchors_XXOD, dims = 1:30)
rm(XXO_2, XXO_3, XXO_4, XXOD)

XXO <- c(XXOC_integrated, XXOD_integrated)
anchors_XXO <- FindIntegrationAnchors(object.list = XXO, dims = 1:30)
XXO_integrated <- IntegrateData(anchorset = anchors_XXO, dims = 1:30)
rm(XXOC_integrated, XXOD_integrated, XXO)

XYOC <- c(XYO_1, XYO_2, XYO_5)
anchors_XYOC <- FindIntegrationAnchors(object.list = XYOC, dims = 1:30)
XYOC_integrated <- IntegrateData(anchorset = anchors_XYOC, dims = 1:30)
rm(XYO_1, XYO_2, XYO_5, XYOC)

XYOD <- c(XYO_3, XYO_4, XYO_6)
anchors_XYOD <- FindIntegrationAnchors(object.list = XYOD, dims = 1:30)
XYOD_integrated <- IntegrateData(anchorset = anchors_XYOD, dims = 1:30)
rm(XYO_3, XYO_4, XYO_6, XYOD)

XYO <- c(XYOC_integrated, XYOD_integrated)
anchors_XYO <- FindIntegrationAnchors(object.list = XYO, dims = 1:30)
XYO_integrated <- IntegrateData(anchorset = anchors_XYO, dims = 1:30)
rm(XYOC_integrated, XYOD_integrated, XYO)

XXTC <- c(XXT_1, XXT_5)
anchors_XXTC <- FindIntegrationAnchors(object.list = XXTC, dims = 1:30)
XXTC_integrated <- IntegrateData(anchorset = anchors_XXTC, dims = 1:30)
rm(XXT_1, XXT_5, XXTC)

XXTD <- c(XXT_2, XXT_3, XXT_4)
anchors_XXTD <- FindIntegrationAnchors(object.list = XXTD, dims = 1:30)
XXTD_integrated <- IntegrateData(anchorset = anchors_XXTD, dims = 1:30)
rm(XXT_2, XXT_3, XXT_4, XXTD)

XXT <- c(XXTC_integrated, XXTD_integrated)
anchors_XXT <- FindIntegrationAnchors(object.list = XXT, dims = 1:30)
XXT_integrated <- IntegrateData(anchorset = anchors_XXT, dims = 1:30)
rm(XXTC_integrated, XXTD_integrated, XXT)

XYTC <- c(XYT_1, XYT_2, XYT_4)
anchors_XYTC <- FindIntegrationAnchors(object.list = XYTC, dims = 1:30)
XYTC_integrated <- IntegrateData(anchorset = anchors_XYTC, dims = 1:30)
rm(XYT_1, XYT_2, XYT_4, XYTC)

XYTD <- c(XYT_3, XYT_5, XYT_6)
anchors_XYTD <- FindIntegrationAnchors(object.list = XYTD, dims = 1:30)
XYTD_integrated <- IntegrateData(anchorset = anchors_XYTD, dims = 1:30)
rm(XYT_3, XYT_5, XYT_6, XYTD)

XYT <- c(XYTC_integrated, XYTD_integrated)
anchors_XYT <- FindIntegrationAnchors(object.list = XYT, dims = 1:30)
XYT_integrated <- IntegrateData(anchorset = anchors_XYT, dims = 1:30)
rm(XYTC_integrated, XYTD_integrated, XYT)

FCG <- c(XXO_integrated, XYO_integrated, XXT_integrated, XYT_integrated)
anchors_FCG <- FindIntegrationAnchors(object.list = FCG, dims = 1:30)
FCG_integrated <- IntegrateData(anchorset = anchors_FCG, dims = 1:30)
rm(XXO_integrated, XYO_integrated, XXT_integrated, XYT_integrated, FCG)

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

saveRDS(FCG_integrated, file = 'FCG_integrated_PCA_0.1.rds')
FCG_integrated <- readRDS(file = "FCG_integrated_PCA_0.1.rds")

Idents(FCG_integrated) <- "orig.ident"
FCG_integrated <- RenameIdents(FCG_integrated,
                              `FCG_787` = "Ctrl", `FCG_803`="Demye", `FCG_815`="Demye", `FCG_837`="Demye",`FCG_878`="Ctrl",`FCG_881`="Ctrl",
                              `FCG_839` = "Ctrl", `FCG_859`="Ctrl", `FCG_866`="Demye", `FCG_882`="Demye",`FCG_895`="Ctrl",`FCG_896`="Demye",
                              `FCG_832` = "Ctrl", `FCG_845`="Demye", `FCG_846`="Demye", `FCG_899`="Demye",`FCG_902`="Ctrl",
                              `FCG_854` = "Ctrl", `FCG_863`="Ctrl", `FCG_869`="Demye", `FCG_891`="Ctrl",`FCG_900`="Demye",`FCG_901`="Demye"
)

FCG_integrated$Food <- Idents(FCG_integrated)
FCG_integrated$Genotype.Food <- paste(FCG_integrated$Condition, Idents(FCG_integrated), sep = "_")

FCG_integrated$Food <- factor(x = FCG_integrated$Food, levels = c("Ctrl","Demye"))
FCG_integrated$Condition <- factor(x = FCG_integrated$Condition, levels = c("XXO","XYO","XXT","XYT"))
FCG_integrated$Genotype.Food <- factor(x = FCG_integrated$Genotype.Food, levels = c("XXO_Ctrl","XXO_Demye","XYO_Ctrl","XYO_Demye",
                                                                                    "XXT_Ctrl","XXT_Demye","XYT_Ctrl","XYT_Demye"))

Idents(FCG_integrated) <- "orig.ident"
FCG_integrated <- RenameIdents(FCG_integrated,
                               `FCG_787` = "XXO_Ctrl_1", `FCG_803`="XXO_Demye_1", `FCG_815`="XXO_Demye_2", `FCG_837`="XXO_Demye_3",`FCG_878`="XXO_Ctrl_2",`FCG_881`="XXO_Ctrl_3",
                               `FCG_839` = "XYO_Ctrl_1", `FCG_859`="XYO_Ctrl_2", `FCG_866`="XYO_Demye_1", `FCG_882`="XYO_Demye_2",`FCG_895`="XYO_Ctrl_3",`FCG_896`="XYO_Demye_3",
                               `FCG_832` = "XXT_Ctrl_1", `FCG_845`="XXT_Demye_1", `FCG_846`="XXT_Demye_2", `FCG_899`="XXT_Demye_3",`FCG_902`="XXT_Ctrl_2",
                               `FCG_854` = "XYT_Ctrl_1", `FCG_863`="XYT_Ctrl_2", `FCG_869`="XYT_Demye_1", `FCG_891`="XYT_Ctrl_3",`FCG_900`="XYT_Demye_2",`FCG_901`="XYT_Demye_3"
)
FCG_integrated$Sample <- Idents(FCG_integrated)
FCG_integrated$Sample <- factor(x = FCG_integrated$Sample, levels = c("XXO_Ctrl_1","XXO_Ctrl_2","XXO_Ctrl_3","XXO_Demye_1","XXO_Demye_2","XXO_Demye_3",
                                                                             "XYO_Ctrl_1","XYO_Ctrl_2","XYO_Ctrl_3","XYO_Demye_1","XYO_Demye_2","XYO_Demye_3",
                                                                             "XXT_Ctrl_1","XXT_Ctrl_2","XXT_Demye_1","XXT_Demye_2","XXT_Demye_3",
                                                                             "XYT_Ctrl_1","XYT_Ctrl_2","XYT_Ctrl_3","XYT_Demye_1","XYT_Demye_2","XYT_Demye_3"))

pdf("FCG_QC.pdf", width=12, height=4)
Idents(FCG_integrated) <- "Condition"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(FCG_integrated) <- "Food"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(FCG_integrated) <- "Genotype.Food"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

pdf("FCG_QC_Sample.pdf", width=18, height=4)
Idents(FCG_integrated) <- "Sample"
VlnPlot(object = FCG_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(FCG_integrated) <- "seurat_clusters"
pdf("FCG_integrated_umap.pdf", width=5, height=4)
DimPlot(FCG_integrated, reduction = 'umap', label = T)
dev.off()
pdf("FCG_integrated_umap_split_individual.pdf", width=12, height=7)
DimPlot(FCG_integrated, reduction = "umap", split.by = "Sample", label = T, ncol = 6)
dev.off()
pdf("FCG_integrated_umap_split_Condition.pdf", width=10, height=4.5)
DimPlot(FCG_integrated, reduction = "umap", split.by = "Genotype.Food", label = T, ncol = 4)
dev.off()

write.csv(table(FCG_integrated$seurat_clusters, FCG_integrated$Sample), "FCG_cell_counts_cluster_by_sample.csv")

saveRDS(FCG_integrated, file = 'FCG_integrated_PCA_0.1.rds')

DefaultAssay(FCG_integrated) <- 'RNA'

FCG_markers <- FindAllMarkers(FCG_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(FCG_markers, "FCG_markers.csv")

FCG_markers <- read.csv(file = "FCG_markers.csv", header=T,row.names =1)
top5 <- FCG_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
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

