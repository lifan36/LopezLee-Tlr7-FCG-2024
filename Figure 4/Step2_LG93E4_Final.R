#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG93E4/DF_2ndRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8238) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_503", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_NTG_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_NTG_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_503" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_NTG_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_NTG_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_NTG_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8656) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_528", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_NTG_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_NTG_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_528" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_NTG_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_NTG_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_NTG_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8886) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_542", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_NTG_E4_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_NTG_E4_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_542" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_NTG_E4_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_NTG_E4_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_NTG_E4_3_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9001) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_621", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_P301S_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_P301S_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_621" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_P301S_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_P301S_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_P301S_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7581) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_409", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_P301S_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_P301S_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_409" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_P301S_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_P301S_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_P301S_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8838) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_539", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_NTG_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_NTG_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_539" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_NTG_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_NTG_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_NTG_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8729) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_532", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_NTG_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_NTG_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_532" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_NTG_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_NTG_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_NTG_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8394) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_512", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_NTG_E4_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_NTG_E4_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_512" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_NTG_E4_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_NTG_E4_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_NTG_E4_3_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8548) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_521", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_P301S_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_P301S_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_521" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_P301S_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_P301S_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_P301S_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_P301S_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9152) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_631", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_P301S_E4_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_P301S_E4_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_631" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_P301S_E4_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_P301S_E4_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_P301S_E4_3_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8357) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_510", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_NTG_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_NTG_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_510" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_NTG_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_NTG_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_NTG_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9396) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_648", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_NTG_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_NTG_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_648" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_NTG_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_NTG_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_NTG_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9372) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_647", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_NTG_E4_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_NTG_E4_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_647" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_NTG_E4_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_NTG_E4_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_NTG_E4_3_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8909) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_543", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_P301S_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_P301S_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_543" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_P301S_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_P301S_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_P301S_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9859) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_680", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_P301S_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_P301S_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_680" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_P301S_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_P301S_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_P301S_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7822) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_422", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_NTG_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_NTG_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_422" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_NTG_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_NTG_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_NTG_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9349) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_645", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_NTG_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_NTG_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_645" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_NTG_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_NTG_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_NTG_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9437) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_651", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_NTG_E4_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_NTG_E4_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_651" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_NTG_E4_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_NTG_E4_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_NTG_E4_3_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9024) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_623", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_P301S_E4_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_P301S_E4_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_623" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_P301S_E4_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_P301S_E4_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_P301S_E4_1_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8711) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_531", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_P301S_E4_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_P301S_E4_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_531" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_P301S_E4_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_P301S_E4_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_P301S_E4_2_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXO_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*4551) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_141", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXO_P301S_E4_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXO_P301S_E4_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.03_141" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXO_P301S_E4_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXO_P301S_E4_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXO_P301S_E4_4_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XXT_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5577) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_218", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XXT_P301S_E4_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XXT_P301S_E4_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_218" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XXT_P301S_E4_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XXT_P301S_E4_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XXT_P301S_E4_4_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYO_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5859) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_229", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYO_P301S_E4_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYO_P301S_E4_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_229" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYO_P301S_E4_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYO_P301S_E4_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYO_P301S_E4_4_singlets_PCA.rds")
###########################################################################
#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound/XYT_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5256) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_205", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("XYT_P301S_E4_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("XYT_P301S_E4_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_205" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"XYT_P301S_E4_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("XYT_P301S_E4_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"XYT_P301S_E4_4_singlets_PCA.rds")
###########################################################################










