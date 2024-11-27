
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
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/DF_2ndRound")

M_Ctrl_1 <- readRDS(file = "M_Ctrl_1_singlets_PCA.rds")
M_Ctrl_2 <- readRDS(file = "M_Ctrl_2_singlets_PCA.rds")
M_Ctrl_3 <- readRDS(file = "M_Ctrl_3_singlets_PCA.rds")

M_3wk_1 <- readRDS(file = "M_3wk_1_singlets_PCA.rds")
M_3wk_2 <- readRDS(file = "M_3wk_2_singlets_PCA.rds")
M_3wk_3 <- readRDS(file = "M_3wk_3_singlets_PCA.rds")

F_Ctrl_1 <- readRDS(file = "F_Ctrl_1_singlets_PCA.rds")
F_Ctrl_2 <- readRDS(file = "F_Ctrl_2_singlets_PCA.rds")
F_Ctrl_3 <- readRDS(file = "F_Ctrl_3_singlets_PCA.rds")

F_3wk_1 <- readRDS(file = "F_3wk_1_singlets_PCA.rds")
F_3wk_2 <- readRDS(file = "F_3wk_2_singlets_PCA.rds")
F_3wk_3 <- readRDS(file = "F_3wk_3_singlets_PCA.rds")


setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/integration")

M_Ctrl <- c(M_Ctrl_1, M_Ctrl_2, M_Ctrl_3)
anchors_M_Ctrl <- FindIntegrationAnchors(object.list = M_Ctrl, dims = 1:30)
M_Ctrl_integrated <- IntegrateData(anchorset = anchors_M_Ctrl, dims = 1:30)
rm(M_Ctrl_1, M_Ctrl_2, M_Ctrl_3, M_Ctrl)

M_3wk <- c(M_3wk_1, M_3wk_2, M_3wk_3)
anchors_M_3wk <- FindIntegrationAnchors(object.list = M_3wk, dims = 1:30)
M_3wk_integrated <- IntegrateData(anchorset = anchors_M_3wk, dims = 1:30)
rm(M_3wk_1, M_3wk_2, M_3wk_3, M_3wk)

F_Ctrl <- c(F_Ctrl_1, F_Ctrl_2, F_Ctrl_3)
anchors_F_Ctrl <- FindIntegrationAnchors(object.list = F_Ctrl, dims = 1:30)
F_Ctrl_integrated <- IntegrateData(anchorset = anchors_F_Ctrl, dims = 1:30)
rm(F_Ctrl_1, F_Ctrl_2, F_Ctrl_3, F_Ctrl)

F_3wk <- c(F_3wk_1, F_3wk_2, F_3wk_3)
anchors_F_3wk <- FindIntegrationAnchors(object.list = F_3wk, dims = 1:30)
F_3wk_integrated <- IntegrateData(anchorset = anchors_F_3wk, dims = 1:30)
rm(F_3wk_1, F_3wk_2, F_3wk_3, F_3wk)

FCG_YoungCPZ <- c(M_Ctrl_integrated,M_3wk_integrated,F_Ctrl_integrated,F_3wk_integrated)
anchors_FCG_YoungCPZ <- FindIntegrationAnchors(object.list = FCG_YoungCPZ, dims = 1:30)
FCG_YoungCPZ_integrated <- IntegrateData(anchorset = anchors_FCG_YoungCPZ, dims = 1:30)
rm(M_Ctrl_integrated,M_3wk_integrated,F_Ctrl_integrated,F_3wk_integrated, FCG_YoungCPZ)

#saveRDS(FCG_YoungCPZ_integrated, file = "FCG_YoungCPZ_integrated.rds")

DefaultAssay(FCG_YoungCPZ_integrated) <- 'integrated'

# FCG_YoungCPZ_integrated <- NormalizeData(FCG_YoungCPZ_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# FCG_YoungCPZ_integrated <- FindVariableFeatures(FCG_YoungCPZ_integrated, selection.method = "vst", nfeatures = 3000)

FCG_YoungCPZ_integrated <- ScaleData(FCG_YoungCPZ_integrated, verbose = FALSE)
FCG_YoungCPZ_integrated <- RunPCA(FCG_YoungCPZ_integrated, features = VariableFeatures(object = FCG_YoungCPZ_integrated), verbose = FALSE)

FCG_YoungCPZ_integrated <- FindNeighbors(FCG_YoungCPZ_integrated, dims = 1:15)
FCG_YoungCPZ_integrated <- FindClusters(FCG_YoungCPZ_integrated, resolution = 0.1)
FCG_YoungCPZ_integrated <- RunUMAP(FCG_YoungCPZ_integrated, dims = 1: 15)

DefaultAssay(FCG_YoungCPZ_integrated) <- 'RNA'
FCG_YoungCPZ_integrated <- NormalizeData(FCG_YoungCPZ_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
FCG_YoungCPZ_integrated <- ScaleData(FCG_YoungCPZ_integrated, features = rownames(FCG_YoungCPZ_integrated))

#saveRDS(FCG_YoungCPZ_integrated, file = 'FCG_YoungCPZ_integrated_PCA_0.1.rds')
#FCG_YoungCPZ_integrated <- readRDS(file = "FCG_YoungCPZ_integrated_PCA_0.1.rds")

FCG_YoungCPZ_integrated$Condition <- factor(x = FCG_YoungCPZ_integrated$Condition, levels = c("M_Ctrl","M_3wk","F_Ctrl","F_3wk"))
FCG_YoungCPZ_integrated$Sample_Name <- factor(x = FCG_YoungCPZ_integrated$Sample_Name, levels = c("M_Ctrl_1","M_Ctrl_2","M_Ctrl_3",
                                                                                      "M_3wk_1","M_3wk_2","M_3wk_3",
                                                                                      "F_Ctrl_1","F_Ctrl_2","F_Ctrl_3",
                                                                                      "F_3wk_1","F_3wk_2","F_3wk_3"))

pdf("FCG_YoungCPZ_QC.pdf", width=9, height=4)
Idents(FCG_YoungCPZ_integrated) <- "Condition"
VlnPlot(object = FCG_YoungCPZ_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(FCG_YoungCPZ_integrated) <- "Sample_Name"
pdf("FCG_YoungCPZ_QC_Sample.pdf", width=12, height=4)

VlnPlot(object = FCG_YoungCPZ_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(FCG_YoungCPZ_integrated) <- "seurat_clusters"
pdf("FCG_YoungCPZ_integrated_umap.pdf", width=5, height=4)
DimPlot(FCG_YoungCPZ_integrated, reduction = 'umap', label = T)
dev.off()
pdf("FCG_YoungCPZ_integrated_umap_split_individual.pdf", width=15, height=5.5)
DimPlot(FCG_YoungCPZ_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
pdf("FCG_YoungCPZ_integrated_umap_split_Condition.pdf", width=6, height=6)
DimPlot(FCG_YoungCPZ_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

write.csv(table(FCG_YoungCPZ_integrated$seurat_clusters, FCG_YoungCPZ_integrated$Sample_Name), "FCG_YoungCPZ_cell_counts_cluster_by_sample.csv")

DefaultAssay(FCG_YoungCPZ_integrated) <- 'RNA'

FCG_YoungCPZ_markers <- FindAllMarkers(FCG_YoungCPZ_integrated, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1.5, test.use = "MAST")
write.csv(FCG_YoungCPZ_markers, "FCG_YoungCPZ_markers.csv")


saveRDS(FCG_YoungCPZ_integrated, file = 'FCG_YoungCPZ_integrated_PCA_0.1.rds')

FCG_YoungCPZ_markers <- read.csv(file = "FCG_YoungCPZ_markers.csv", header=T,row.names =1)
top5 <- FCG_YoungCPZ_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("FCG_YoungCPZ_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(FCG_YoungCPZ_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("FCG_YoungCPZ_annotation_combine.pdf", width=10, height=5)
DotPlot(object = FCG_YoungCPZ_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

