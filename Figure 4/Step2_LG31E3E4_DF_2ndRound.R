#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/DF_2ndRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_NTG_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_NTG_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6790) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_312", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_NTG_WT_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_NTG_WT_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_312" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_NTG_WT_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_NTG_WT_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_NTG_WT_1_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_P301S_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7865) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_425", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_P301S_WT_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_P301S_WT_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_425" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_P301S_WT_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_P301S_WT_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_P301S_WT_1_singlets_PCA.rds")
###########################################################################
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_P301S_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7146) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_386", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_P301S_WT_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_P301S_WT_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_386" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_P301S_WT_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_P301S_WT_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_P301S_WT_2_singlets_PCA.rds")
###########################################################################
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_P301S_R47H_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10932) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_831", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_P301S_R47H_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_P301S_R47H_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_831" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_P301S_R47H_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_P301S_R47H_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_P301S_R47H_1_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_P301S_R47H_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*4894) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_152", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_P301S_R47H_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_P301S_R47H_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_152" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_P301S_R47H_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_P301S_R47H_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_P301S_R47H_2_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E3_P301S_R47H_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8171) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_498", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E3_P301S_R47H_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E3_P301S_R47H_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_498" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E3_P301S_R47H_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E3_P301S_R47H_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E3_P301S_R47H_3_singlets_PCA.rds")

###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_NTG_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_NTG_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6975) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_321", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_NTG_WT_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_NTG_WT_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_321" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_NTG_WT_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_NTG_WT_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_NTG_WT_1_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_NTG_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_NTG_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7715) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_417", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_NTG_WT_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_NTG_WT_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_417" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_NTG_WT_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_NTG_WT_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_NTG_WT_2_singlets_PCA.rds")
###########################################################################
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.121*16668) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.3, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_2017", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_WT_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_WT_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.3_2017" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_WT_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_WT_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_WT_1_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8227) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_502", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_WT_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_WT_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_502" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_WT_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_WT_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_WT_2_singlets_PCA.rds")
###########################################################################
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_WT_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6983) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_321", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_WT_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_WT_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_321" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_WT_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_WT_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_WT_3_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_WT_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7737) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_418", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_WT_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_WT_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_418" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_WT_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_WT_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_WT_4_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_R47H_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7927) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_428", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_R47H_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_R47H_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_428" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_R47H_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_R47H_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_R47H_1_singlets_PCA.rds")
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_R47H_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6889) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_317", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_R47H_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_R47H_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_317" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_R47H_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_R47H_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_R47H_2_singlets_PCA.rds")
###########################################################################
###########################################################################
all <- readRDS(file = "/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound/E4_P301S_R47H_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9215) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_636", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("E4_P301S_R47H_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("E4_P301S_R47H_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_636" #visualizing the singlet vs doublet cells
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"E4_P301S_R47H_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
singlets <- ScaleData(object = singlets)
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
pdf("E4_P301S_R47H_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
saveRDS(singlets,"E4_P301S_R47H_3_singlets_PCA.rds")
###########################################################################













