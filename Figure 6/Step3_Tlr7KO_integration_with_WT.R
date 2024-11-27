
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====
setwd("/athena/ganlab/scratch/lif4001/FCG/DF_2ndRound")

WT_XXO_Ctrl_1 <- readRDS(file = "XXO_1_singlets_PCA.rds")
WT_XXO_Ctrl_2 <- readRDS(file = "XXO_5_singlets_PCA.rds")
WT_XXO_Ctrl_3 <- readRDS(file = "XXO_6_singlets_PCA.rds")
WT_XXO_CPZ_1 <- readRDS(file = "XXO_2_singlets_PCA.rds")
WT_XXO_CPZ_2 <- readRDS(file = "XXO_3_singlets_PCA.rds")
WT_XXO_CPZ_3 <- readRDS(file = "XXO_4_singlets_PCA.rds")

WT_XYT_Ctrl_1 <- readRDS(file = "XYT_1_singlets_PCA.rds")
WT_XYT_Ctrl_2 <- readRDS(file = "XYT_2_singlets_PCA.rds")
WT_XYT_Ctrl_3 <- readRDS(file = "XYT_4_singlets_PCA.rds")
WT_XYT_CPZ_1 <- readRDS(file = "XYT_3_singlets_PCA.rds")
WT_XYT_CPZ_2 <- readRDS(file = "XYT_5_singlets_PCA.rds")
WT_XYT_CPZ_3 <- readRDS(file = "XYT_6_singlets_PCA.rds")

WT_XXO_Ctrl_1[["Condition"]] = c('WT_XXO_Ctrl')
WT_XXO_Ctrl_1[["Sample_Name"]] = c('WT_XXO_Ctrl_1')
WT_XXO_Ctrl_2[["Condition"]] = c('WT_XXO_Ctrl')
WT_XXO_Ctrl_2[["Sample_Name"]] = c('WT_XXO_Ctrl_2')
WT_XXO_Ctrl_3[["Condition"]] = c('WT_XXO_Ctrl')
WT_XXO_Ctrl_3[["Sample_Name"]] = c('WT_XXO_Ctrl_3')
WT_XXO_CPZ_1[["Condition"]] = c('WT_XXO_CPZ')
WT_XXO_CPZ_1[["Sample_Name"]] = c('WT_XXO_CPZ_1')
WT_XXO_CPZ_2[["Condition"]] = c('WT_XXO_CPZ')
WT_XXO_CPZ_2[["Sample_Name"]] = c('WT_XXO_CPZ_2')
WT_XXO_CPZ_3[["Condition"]] = c('WT_XXO_CPZ')
WT_XXO_CPZ_3[["Sample_Name"]] = c('WT_XXO_CPZ_3')

WT_XYT_Ctrl_1[["Condition"]] = c('WT_XYT_Ctrl')
WT_XYT_Ctrl_1[["Sample_Name"]] = c('WT_XYT_Ctrl_1')
WT_XYT_Ctrl_2[["Condition"]] = c('WT_XYT_Ctrl')
WT_XYT_Ctrl_2[["Sample_Name"]] = c('WT_XYT_Ctrl_2')
WT_XYT_Ctrl_3[["Condition"]] = c('WT_XYT_Ctrl')
WT_XYT_Ctrl_3[["Sample_Name"]] = c('WT_XYT_Ctrl_3')
WT_XYT_CPZ_1[["Condition"]] = c('WT_XYT_CPZ')
WT_XYT_CPZ_1[["Sample_Name"]] = c('WT_XYT_CPZ_1')
WT_XYT_CPZ_2[["Condition"]] = c('WT_XYT_CPZ')
WT_XYT_CPZ_2[["Sample_Name"]] = c('WT_XYT_CPZ_2')
WT_XYT_CPZ_3[["Condition"]] = c('WT_XYT_CPZ')
WT_XYT_CPZ_3[["Sample_Name"]] = c('WT_XYT_CPZ_3')


setwd("/athena/ganlab/scratch/lif4001/Tlr7/DF_2ndRound")

LG70_NY_23 <- readRDS(file = "LG70_NY_23_singlets_PCA.rds")
LG70_NY_24 <- readRDS(file = "LG70_NY_24_singlets_PCA.rds")
LG70_NY_25 <- readRDS(file = "LG70_NY_25_singlets_PCA.rds")
LG70_NY_26 <- readRDS(file = "LG70_NY_26_singlets_PCA.rds")
LG70_NY_33 <- readRDS(file = "LG70_NY_33_singlets_PCA.rds")
LG70_NY_34 <- readRDS(file = "LG70_NY_34_singlets_PCA.rds")
LG70_NY_36 <- readRDS(file = "LG70_NY_36_singlets_PCA.rds")
LG70_NY_37 <- readRDS(file = "LG70_NY_37_singlets_PCA.rds")
LG70_NY_38 <- readRDS(file = "LG70_NY_38_singlets_PCA.rds")
LG70_NY_39 <- readRDS(file = "LG70_NY_39_singlets_PCA.rds")
LG70_NY_40 <- readRDS(file = "LG70_NY_40_singlets_PCA.rds")
LG70_NY_41 <- readRDS(file = "LG70_NY_41_singlets_PCA.rds")


LG70_NY_39[["Condition"]] = c('KO_XXO_Ctrl')
LG70_NY_39[["Sample_Name"]] = c('KO_XXO_Ctrl_1')
LG70_NY_40[["Condition"]] = c('KO_XXO_Ctrl')
LG70_NY_40[["Sample_Name"]] = c('KO_XXO_Ctrl_2')
LG70_NY_41[["Condition"]] = c('KO_XXO_Ctrl')
LG70_NY_41[["Sample_Name"]] = c('KO_XXO_Ctrl_3')

LG70_NY_23[["Condition"]] = c('KO_XXO_CPZ')
LG70_NY_23[["Sample_Name"]] = c('KO_XXO_CPZ_1')
LG70_NY_33[["Condition"]] = c('KO_XXO_CPZ')
LG70_NY_33[["Sample_Name"]] = c('KO_XXO_CPZ_2')
LG70_NY_34[["Condition"]] = c('KO_XXO_CPZ')
LG70_NY_34[["Sample_Name"]] = c('KO_XXO_CPZ_3')

LG70_NY_36[["Condition"]] = c('KO_XYT_Ctrl')
LG70_NY_36[["Sample_Name"]] = c('KO_XYT_Ctrl_1')
LG70_NY_37[["Condition"]] = c('KO_XYT_Ctrl')
LG70_NY_37[["Sample_Name"]] = c('KO_XYT_Ctrl_2')
LG70_NY_38[["Condition"]] = c('KO_XYT_Ctrl')
LG70_NY_38[["Sample_Name"]] = c('KO_XYT_Ctrl_3')

LG70_NY_24[["Condition"]] = c('KO_XYT_CPZ')
LG70_NY_24[["Sample_Name"]] = c('KO_XYT_CPZ_1')
LG70_NY_25[["Condition"]] = c('KO_XYT_CPZ')
LG70_NY_25[["Sample_Name"]] = c('KO_XYT_CPZ_2')
LG70_NY_26[["Condition"]] = c('KO_XYT_CPZ')
LG70_NY_26[["Sample_Name"]] = c('KO_XYT_CPZ_3')


setwd("/athena/ganlab/scratch/lif4001/Tlr7/integration_with_WT")

WT_XXO_Ctrl <- c(WT_XXO_Ctrl_1, WT_XXO_Ctrl_2, WT_XXO_Ctrl_3)
anchors_WT_XXO_Ctrl <- FindIntegrationAnchors(object.list = WT_XXO_Ctrl, dims = 1:30)
WT_XXO_Ctrl_integrated <- IntegrateData(anchorset = anchors_WT_XXO_Ctrl, dims = 1:30)
rm(WT_XXO_Ctrl_1, WT_XXO_Ctrl_2, WT_XXO_Ctrl_3, WT_XXO_Ctrl)

WT_XXO_CPZ <- c(WT_XXO_CPZ_1, WT_XXO_CPZ_2, WT_XXO_CPZ_3)
anchors_WT_XXO_CPZ <- FindIntegrationAnchors(object.list = WT_XXO_CPZ, dims = 1:30)
WT_XXO_CPZ_integrated <- IntegrateData(anchorset = anchors_WT_XXO_CPZ, dims = 1:30)
rm(WT_XXO_CPZ_1, WT_XXO_CPZ_2, WT_XXO_CPZ_3, WT_XXO_CPZ)

WT_XYT_Ctrl <- c(WT_XYT_Ctrl_1, WT_XYT_Ctrl_2, WT_XYT_Ctrl_3)
anchors_WT_XYT_Ctrl <- FindIntegrationAnchors(object.list = WT_XYT_Ctrl, dims = 1:30)
WT_XYT_Ctrl_integrated <- IntegrateData(anchorset = anchors_WT_XYT_Ctrl, dims = 1:30)
rm(WT_XYT_Ctrl_1, WT_XYT_Ctrl_2, WT_XYT_Ctrl_3, WT_XYT_Ctrl)

WT_XYT_CPZ <- c(WT_XYT_CPZ_1, WT_XYT_CPZ_2, WT_XYT_CPZ_3)
anchors_WT_XYT_CPZ <- FindIntegrationAnchors(object.list = WT_XYT_CPZ, dims = 1:30)
WT_XYT_CPZ_integrated <- IntegrateData(anchorset = anchors_WT_XYT_CPZ, dims = 1:30)
rm(WT_XYT_CPZ_1, WT_XYT_CPZ_2, WT_XYT_CPZ_3, WT_XYT_CPZ)

F_Ctrl <- c(LG70_NY_39, LG70_NY_40, LG70_NY_41)
anchors_F_Ctrl <- FindIntegrationAnchors(object.list = F_Ctrl, dims = 1:30)
F_Ctrl_integrated <- IntegrateData(anchorset = anchors_F_Ctrl, dims = 1:30)
rm(LG70_NY_39, LG70_NY_40, LG70_NY_41, F_Ctrl)

F_Cup_3week <- c(LG70_NY_23, LG70_NY_33, LG70_NY_34)
anchors_F_Cup_3week <- FindIntegrationAnchors(object.list = F_Cup_3week, dims = 1:30)
F_Cup_3week_integrated <- IntegrateData(anchorset = anchors_F_Cup_3week, dims = 1:30)
rm(LG70_NY_23, LG70_NY_33, LG70_NY_34, F_Cup_3week)

M_Ctrl <- c(LG70_NY_36, LG70_NY_37, LG70_NY_38)
anchors_M_Ctrl <- FindIntegrationAnchors(object.list = M_Ctrl, dims = 1:30)
M_Ctrl_integrated <- IntegrateData(anchorset = anchors_M_Ctrl, dims = 1:30)
rm(LG70_NY_36, LG70_NY_37, LG70_NY_38, M_Ctrl)

M_Cup_3week <- c(LG70_NY_24, LG70_NY_25, LG70_NY_26)
anchors_M_Cup_3week <- FindIntegrationAnchors(object.list = M_Cup_3week, dims = 1:30)
M_Cup_3week_integrated <- IntegrateData(anchorset = anchors_M_Cup_3week, dims = 1:30)
rm(LG70_NY_24, LG70_NY_25, LG70_NY_26, M_Cup_3week)

Tlr7 <- c(WT_XXO_Ctrl_integrated, WT_XXO_CPZ_integrated, WT_XYT_Ctrl_integrated, WT_XYT_CPZ_integrated, F_Ctrl_integrated, F_Cup_3week_integrated,
          M_Ctrl_integrated, M_Cup_3week_integrated)
anchors_Tlr7 <- FindIntegrationAnchors(object.list = Tlr7, dims = 1:30)
Tlr7_integrated <- IntegrateData(anchorset = anchors_Tlr7, dims = 1:30)
rm(WT_XXO_Ctrl_integrated, WT_XXO_CPZ_integrated, WT_XYT_Ctrl_integrated, WT_XYT_CPZ_integrated, F_Ctrl_integrated, F_Cup_3week_integrated,
   M_Ctrl_integrated, M_Cup_3week_integrated, Tlr7)

saveRDS(Tlr7_integrated, file = "Tlr7_integrated.rds")

DefaultAssay(Tlr7_integrated) <- 'integrated'

# Tlr7_integrated <- NormalizeData(Tlr7_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# Tlr7_integrated <- FindVariableFeatures(Tlr7_integrated, selection.method = "vst", nfeatures = 3000)

Tlr7_integrated <- ScaleData(Tlr7_integrated, verbose = FALSE)
Tlr7_integrated <- RunPCA(Tlr7_integrated, features = VariableFeatures(object = Tlr7_integrated), verbose = FALSE)

Tlr7_integrated <- FindNeighbors(Tlr7_integrated, dims = 1:15)
Tlr7_integrated <- FindClusters(Tlr7_integrated, resolution = 0.1)
Tlr7_integrated <- RunUMAP(Tlr7_integrated, dims = 1: 15)

DefaultAssay(Tlr7_integrated) <- 'RNA'
Tlr7_integrated <- NormalizeData(Tlr7_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Tlr7_integrated <- ScaleData(Tlr7_integrated, features = rownames(Tlr7_integrated))

#saveRDS(Tlr7_integrated, file = 'Tlr7_integrated_PCA_0.1.rds')
#Tlr7_integrated <- readRDS(file = "Tlr7_integrated_PCA_0.1.rds")

Tlr7_integrated$Condition <- factor(x = Tlr7_integrated$Condition, levels = c("WT_XXO_Ctrl","WT_XXO_CPZ","WT_XYT_Ctrl","WT_XYT_CPZ",
                                                                              "KO_XXO_Ctrl","KO_XXO_CPZ","KO_XYT_Ctrl","KO_XYT_CPZ"))
Tlr7_integrated$Sample_Name <- factor(x = Tlr7_integrated$Sample_Name, levels = c("WT_XXO_Ctrl_1","WT_XXO_Ctrl_2","WT_XXO_Ctrl_3",
                                                                                  "WT_XXO_CPZ_1","WT_XXO_CPZ_2","WT_XXO_CPZ_3",
                                                                                  "WT_XYT_Ctrl_1","WT_XYT_Ctrl_2","WT_XYT_Ctrl_3",
                                                                                  "WT_XYT_CPZ_1","WT_XYT_CPZ_2","WT_XYT_CPZ_3",
                                                                                  "KO_XXO_Ctrl_1","KO_XXO_Ctrl_2","KO_XXO_Ctrl_3",
                                                                                  "KO_XXO_CPZ_1","KO_XXO_CPZ_2","KO_XXO_CPZ_3",
                                                                                  "KO_XYT_Ctrl_1","KO_XYT_Ctrl_2","KO_XYT_Ctrl_3",
                                                                                  "KO_XYT_CPZ_1","KO_XYT_CPZ_2","KO_XYT_CPZ_3"))

pdf("Tlr7_QC.pdf", width=9, height=4)
Idents(Tlr7_integrated) <- "Condition"
VlnPlot(object = Tlr7_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(Tlr7_integrated) <- "Sample_Name"
pdf("Tlr7_QC_Sample.pdf", width=18, height=4)

VlnPlot(object = Tlr7_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(Tlr7_integrated) <- "seurat_clusters"
pdf("Tlr7_integrated_umap.pdf", width=5, height=4)
DimPlot(Tlr7_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Tlr7_integrated_umap_split_individual.pdf", width=12, height=8)
DimPlot(Tlr7_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 6)
dev.off()
pdf("Tlr7_integrated_umap_split_Condition.pdf", width=10, height=5)
DimPlot(Tlr7_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()

write.csv(table(Tlr7_integrated$seurat_clusters, Tlr7_integrated$Sample_Name), "Tlr7_cell_counts_cluster_by_sample.csv")

DefaultAssay(Tlr7_integrated) <- 'RNA'
#Add marker genes
sig_EN<-c("Snap25","Vsnl1","Dnm1","Pde1a","Cplx1","Nap1l5","Erbb4","Slc6a1","Hpca","Grin2b","Rasgrp1","Ppp3r1","Camk4","Brinp1",
          "Nrgn","C1ql3","Cplx2","Rbfox1","Prox1","Slc17a7", "Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", 
          "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r","Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("Tlr7_annotation_combine.pdf", width=20, height=5)
DotPlot(object = Tlr7_integrated, features = markers.to.plot) + RotatedAxis()
dev.off()


Tlr7_markers <- FindAllMarkers(Tlr7_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(Tlr7_markers, "Tlr7_markers.csv")


saveRDS(Tlr7_integrated, file = 'Tlr7_integrated_PCA_0.1.rds')

#Tlr7_integrated <- readRDS(file = "Tlr7_integrated_PCA_0.1.rds")
#Tlr7_markers <- read.csv(file = "Tlr7_markers.csv", header=T,row.names =1)
#top5 <- Tlr7_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top5$gene <- as.character(top5$gene)
#pdf("Tlr7_HeatMapTop5_0.1_new.pdf", width=24, height=16)
#DoHeatmap(Tlr7_integrated, features = top5$gene) + NoLegend()
#dev.off()




