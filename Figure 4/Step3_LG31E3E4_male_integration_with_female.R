
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
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/DF_2ndRound")

LG31E4_R47H_78 <- readRDS(file = "LG31E4_R47H_78_singlets_PCA.rds")
LG31E4_R47H_80 <- readRDS(file = "LG31E4_R47H_80_singlets_PCA.rds")
LG31E4_R47H_94 <- readRDS(file = "LG31E4_R47H_94_singlets_PCA.rds")
LG31E4_R47H_60 <- readRDS(file = "LG31E4_R47H_60_singlets_PCA.rds")
LG31E4_R47H_61 <- readRDS(file = "LG31E4_R47H_61_singlets_PCA.rds")
LG31E4_R47H_82 <- readRDS(file = "LG31E4_R47H_82_singlets_PCA.rds")

LG31E4_R47H_78[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_78[["Sample_Name"]] = c('M_E4_NTG_WT_1')
LG31E4_R47H_80[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_80[["Sample_Name"]] = c('M_E4_NTG_WT_2')
LG31E4_R47H_94[["Condition"]] = c('M_E4_NTG_WT')
LG31E4_R47H_94[["Sample_Name"]] = c('M_E4_NTG_WT_3')

LG31E4_R47H_60[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_60[["Sample_Name"]] = c('M_E4_P301S_WT_1')
LG31E4_R47H_61[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_61[["Sample_Name"]] = c('M_E4_P301S_WT_2')
LG31E4_R47H_82[["Condition"]] = c('M_E4_P301S_WT')
LG31E4_R47H_82[["Sample_Name"]] = c('M_E4_P301S_WT_3')


setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/DF_2ndRound")

E4_NTG_WT_1 <- readRDS(file = "E4_NTG_WT_1_singlets_PCA.rds")
E4_NTG_WT_2 <- readRDS(file = "E4_NTG_WT_2_singlets_PCA.rds")
E4_NTG_WT_3 <- readRDS(file = "E4_NTG_WT_3_singlets_PCA.rds")

#E4_P301S_WT_1 <- readRDS(file = "E4_P301S_WT_1_singlets_PCA.rds")
E4_P301S_WT_1 <- readRDS(file = "E4_P301S_WT_2_singlets_PCA.rds")
E4_P301S_WT_2 <- readRDS(file = "E4_P301S_WT_3_singlets_PCA.rds")
E4_P301S_WT_3 <- readRDS(file = "E4_P301S_WT_4_singlets_PCA.rds")

E4_NTG_WT_1[["Condition"]] = c('F_E4_NTG_WT')
E4_NTG_WT_1[["Sample_Name"]] = c('F_E4_NTG_WT_1')
E4_NTG_WT_2[["Condition"]] = c('F_E4_NTG_WT')
E4_NTG_WT_2[["Sample_Name"]] = c('F_E4_NTG_WT_2')
E4_NTG_WT_3[["Condition"]] = c('F_E4_NTG_WT')
E4_NTG_WT_3[["Sample_Name"]] = c('F_E4_NTG_WT_3')

E4_P301S_WT_1[["Condition"]] = c('F_E4_P301S_WT')
E4_P301S_WT_1[["Sample_Name"]] = c('F_E4_P301S_WT_1')
E4_P301S_WT_2[["Condition"]] = c('F_E4_P301S_WT')
E4_P301S_WT_2[["Sample_Name"]] = c('F_E4_P301S_WT_2')
E4_P301S_WT_3[["Condition"]] = c('F_E4_P301S_WT')
E4_P301S_WT_3[["Sample_Name"]] = c('F_E4_P301S_WT_3')

setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/integration_with_female")

M_E4_NTG_WT <- c(LG31E4_R47H_78, LG31E4_R47H_80, LG31E4_R47H_94)
anchors_M_E4_NTG_WT <- FindIntegrationAnchors(object.list = M_E4_NTG_WT, dims = 1:30)
M_E4_NTG_WT_integrated <- IntegrateData(anchorset = anchors_M_E4_NTG_WT, dims = 1:30)
rm(LG31E4_R47H_78, LG31E4_R47H_80, LG31E4_R47H_94, M_E4_NTG_WT)

M_E4_P301S_WT <- c(LG31E4_R47H_60, LG31E4_R47H_61, LG31E4_R47H_82)
anchors_M_E4_P301S_WT <- FindIntegrationAnchors(object.list = M_E4_P301S_WT, dims = 1:30)
M_E4_P301S_WT_integrated <- IntegrateData(anchorset = anchors_M_E4_P301S_WT, dims = 1:30)
rm(LG31E4_R47H_60, LG31E4_R47H_61, LG31E4_R47H_82, M_E4_P301S_WT)

E4_NTG_WT <- c(E4_NTG_WT_1, E4_NTG_WT_2, E4_NTG_WT_3)
anchors_E4_NTG_WT <- FindIntegrationAnchors(object.list = E4_NTG_WT, dims = 1:30)
E4_NTG_WT_integrated <- IntegrateData(anchorset = anchors_E4_NTG_WT, dims = 1:30)
rm(E4_NTG_WT_1, E4_NTG_WT_2, E4_NTG_WT_3, E4_NTG_WT)

E4_P301S_WT <- c(E4_P301S_WT_1, E4_P301S_WT_2, E4_P301S_WT_3)
anchors_E4_P301S_WT <- FindIntegrationAnchors(object.list = E4_P301S_WT, dims = 1:30)
E4_P301S_WT_integrated <- IntegrateData(anchorset = anchors_E4_P301S_WT, dims = 1:30)
rm(E4_P301S_WT_1, E4_P301S_WT_2, E4_P301S_WT_3, E4_P301S_WT)

bothSex <- c(M_E4_NTG_WT_integrated, M_E4_P301S_WT_integrated,E4_NTG_WT_integrated, E4_P301S_WT_integrated)
anchors_bothSex <- FindIntegrationAnchors(object.list = bothSex, dims = 1:30)
bothSex_integrated <- IntegrateData(anchorset = anchors_bothSex, dims = 1:30)
rm(M_E4_NTG_WT_integrated, M_E4_P301S_WT_integrated,E4_NTG_WT_integrated, E4_P301S_WT_integrated, bothSex)

#saveRDS(bothSex_integrated, file = "bothSex_integrated.rds")

DefaultAssay(bothSex_integrated) <- 'integrated'

# bothSex_integrated <- NormalizeData(bothSex_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# bothSex_integrated <- FindVariableFeatures(bothSex_integrated, selection.method = "vst", nfeatures = 3000)

bothSex_integrated <- ScaleData(bothSex_integrated, verbose = FALSE)
bothSex_integrated <- RunPCA(bothSex_integrated, features = VariableFeatures(object = bothSex_integrated), verbose = FALSE)

bothSex_integrated <- FindNeighbors(bothSex_integrated, dims = 1:15)
bothSex_integrated <- FindClusters(bothSex_integrated, resolution = 0.1)
bothSex_integrated <- RunUMAP(bothSex_integrated, dims = 1: 15)

DefaultAssay(bothSex_integrated) <- 'RNA'
bothSex_integrated <- NormalizeData(bothSex_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
bothSex_integrated <- ScaleData(bothSex_integrated, features = rownames(bothSex_integrated))

#saveRDS(bothSex_integrated, file = 'bothSex_integrated_PCA_0.1.rds')
#bothSex_integrated <- readRDS(file = "bothSex_integrated_PCA_0.1.rds")

bothSex_integrated$Condition <- factor(x = bothSex_integrated$Condition, levels = c("M_E4_NTG_WT","M_E4_P301S_WT","F_E4_NTG_WT","F_E4_P301S_WT"))
bothSex_integrated$Sample_Name <- factor(x = bothSex_integrated$Sample_Name, levels = c("M_E4_NTG_WT_1","M_E4_NTG_WT_2","M_E4_NTG_WT_3",
                                                                                  "M_E4_P301S_WT_1","M_E4_P301S_WT_2","M_E4_P301S_WT_3",
                                                                                  "F_E4_NTG_WT_1","F_E4_NTG_WT_2","F_E4_NTG_WT_3",
                                                                                  "F_E4_P301S_WT_1","F_E4_P301S_WT_2","F_E4_P301S_WT_3"))

pdf("bothSex_QC.pdf", width=9, height=4)
Idents(bothSex_integrated) <- "Condition"
VlnPlot(object = bothSex_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(bothSex_integrated) <- "Sample_Name"
pdf("bothSex_QC_Sample.pdf", width=22, height=4)

VlnPlot(object = bothSex_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(bothSex_integrated) <- "seurat_clusters"
pdf("bothSex_integrated_umap.pdf", width=5, height=4)
DimPlot(bothSex_integrated, reduction = 'umap', label = T)
dev.off()
pdf("bothSex_integrated_umap_split_individual.pdf", width=15, height=12)
DimPlot(bothSex_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 5)
dev.off()
pdf("bothSex_integrated_umap_split_Condition.pdf", width=9, height=6)
DimPlot(bothSex_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

write.csv(table(bothSex_integrated$seurat_clusters, bothSex_integrated$Sample_Name), "bothSex_cell_counts_cluster_by_sample.csv")

DefaultAssay(bothSex_integrated) <- 'RNA'

bothSex_markers <- FindAllMarkers(bothSex_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(bothSex_markers, "bothSex_markers.csv")


saveRDS(bothSex_integrated, file = 'bothSex_integrated_PCA_0.1.rds')

bothSex_markers <- read.csv(file = "bothSex_markers.csv", header=T,row.names =1)
top5 <- bothSex_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("bothSex_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(bothSex_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("bothSex_annotation_combine.pdf", width=10, height=5)
DotPlot(object = bothSex_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()
