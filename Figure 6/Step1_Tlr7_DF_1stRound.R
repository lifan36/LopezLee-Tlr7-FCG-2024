#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Tlr7/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_11/outs')
sc = autoEstCont(sc)
LG70_NY_11.counts = adjustCounts(sc)
LG70_NY_11 <- CreateSeuratObject(counts = LG70_NY_11.counts, project = "LG70_NY_11", min.cells = 3, min.features = 200)
rm(LG70_NY_11.counts)
#vizualize QC metrics and filtering====
LG70_NY_11[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_11, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_11
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_11_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_11_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_11_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_11_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_12/outs')
sc = autoEstCont(sc)
LG70_NY_12.counts = adjustCounts(sc)
LG70_NY_12 <- CreateSeuratObject(counts = LG70_NY_12.counts, project = "LG70_NY_12", min.cells = 3, min.features = 200)
rm(LG70_NY_12.counts)
#vizualize QC metrics and filtering====
LG70_NY_12[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_12, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_12
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_12_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_12_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_12_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_12_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_13/outs')
sc = autoEstCont(sc)
LG70_NY_13.counts = adjustCounts(sc)
LG70_NY_13 <- CreateSeuratObject(counts = LG70_NY_13.counts, project = "LG70_NY_13", min.cells = 3, min.features = 200)
rm(LG70_NY_13.counts)
#vizualize QC metrics and filtering====
LG70_NY_13[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_13, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_13
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_13_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_13_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_13_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_13_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_16/outs')
sc = autoEstCont(sc)
LG70_NY_16.counts = adjustCounts(sc)
LG70_NY_16 <- CreateSeuratObject(counts = LG70_NY_16.counts, project = "LG70_NY_16", min.cells = 3, min.features = 200)
rm(LG70_NY_16.counts)
#vizualize QC metrics and filtering====
LG70_NY_16[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_16, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_16
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_16_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_16_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_16_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_16_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_17/outs')
sc = autoEstCont(sc)
LG70_NY_17.counts = adjustCounts(sc)
LG70_NY_17 <- CreateSeuratObject(counts = LG70_NY_17.counts, project = "LG70_NY_17", min.cells = 3, min.features = 200)
rm(LG70_NY_17.counts)
#vizualize QC metrics and filtering====
LG70_NY_17[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_17, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_17
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_17_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_17_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_17_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_17_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_18/outs')
sc = autoEstCont(sc)
LG70_NY_18.counts = adjustCounts(sc)
LG70_NY_18 <- CreateSeuratObject(counts = LG70_NY_18.counts, project = "LG70_NY_18", min.cells = 3, min.features = 200)
rm(LG70_NY_18.counts)
#vizualize QC metrics and filtering====
LG70_NY_18[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_18, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_18
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_18_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_18_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_18_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_18_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_19/outs')
sc = autoEstCont(sc)
LG70_NY_19.counts = adjustCounts(sc)
LG70_NY_19 <- CreateSeuratObject(counts = LG70_NY_19.counts, project = "LG70_NY_19", min.cells = 3, min.features = 200)
rm(LG70_NY_19.counts)
#vizualize QC metrics and filtering====
LG70_NY_19[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_19, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_19
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_19_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_19_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_19_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_19_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_20/outs')
sc = autoEstCont(sc)
LG70_NY_20.counts = adjustCounts(sc)
LG70_NY_20 <- CreateSeuratObject(counts = LG70_NY_20.counts, project = "LG70_NY_20", min.cells = 3, min.features = 200)
rm(LG70_NY_20.counts)
#vizualize QC metrics and filtering====
LG70_NY_20[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_20, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_20
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_20_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_20_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_20_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_20_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_21/outs')
sc = autoEstCont(sc)
LG70_NY_21.counts = adjustCounts(sc)
LG70_NY_21 <- CreateSeuratObject(counts = LG70_NY_21.counts, project = "LG70_NY_21", min.cells = 3, min.features = 200)
rm(LG70_NY_21.counts)
#vizualize QC metrics and filtering====
LG70_NY_21[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_21, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_21
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_21_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_21_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_21_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_21_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_23/outs')
sc = autoEstCont(sc)
LG70_NY_23.counts = adjustCounts(sc)
LG70_NY_23 <- CreateSeuratObject(counts = LG70_NY_23.counts, project = "LG70_NY_23", min.cells = 3, min.features = 200)
rm(LG70_NY_23.counts)
#vizualize QC metrics and filtering====
LG70_NY_23[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_23, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_23
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_23_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_23_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_23_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_23_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_24/outs')
sc = autoEstCont(sc)
LG70_NY_24.counts = adjustCounts(sc)
LG70_NY_24 <- CreateSeuratObject(counts = LG70_NY_24.counts, project = "LG70_NY_24", min.cells = 3, min.features = 200)
rm(LG70_NY_24.counts)
#vizualize QC metrics and filtering====
LG70_NY_24[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_24, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_24
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_24_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_24_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_24_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_24_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_25/outs')
sc = autoEstCont(sc)
LG70_NY_25.counts = adjustCounts(sc)
LG70_NY_25 <- CreateSeuratObject(counts = LG70_NY_25.counts, project = "LG70_NY_25", min.cells = 3, min.features = 200)
rm(LG70_NY_25.counts)
#vizualize QC metrics and filtering====
LG70_NY_25[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_25, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_25
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_25_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_25_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_25_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_25_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_26/outs')
sc = autoEstCont(sc)
LG70_NY_26.counts = adjustCounts(sc)
LG70_NY_26 <- CreateSeuratObject(counts = LG70_NY_26.counts, project = "LG70_NY_26", min.cells = 3, min.features = 200)
rm(LG70_NY_26.counts)
#vizualize QC metrics and filtering====
LG70_NY_26[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_26, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_26
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_26_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_26_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_26_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_26_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_30/outs')
sc = autoEstCont(sc)
LG70_NY_30.counts = adjustCounts(sc)
LG70_NY_30 <- CreateSeuratObject(counts = LG70_NY_30.counts, project = "LG70_NY_30", min.cells = 3, min.features = 200)
rm(LG70_NY_30.counts)
#vizualize QC metrics and filtering====
LG70_NY_30[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_30, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_30
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_30_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_30_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_30_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_30_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_31/outs')
sc = autoEstCont(sc)
LG70_NY_31.counts = adjustCounts(sc)
LG70_NY_31 <- CreateSeuratObject(counts = LG70_NY_31.counts, project = "LG70_NY_31", min.cells = 3, min.features = 200)
rm(LG70_NY_31.counts)
#vizualize QC metrics and filtering====
LG70_NY_31[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_31, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_31
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_31_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_31_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_31_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_31_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_32/outs')
sc = autoEstCont(sc)
LG70_NY_32.counts = adjustCounts(sc)
LG70_NY_32 <- CreateSeuratObject(counts = LG70_NY_32.counts, project = "LG70_NY_32", min.cells = 3, min.features = 200)
rm(LG70_NY_32.counts)
#vizualize QC metrics and filtering====
LG70_NY_32[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_32, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_32
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_32_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_32_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_32_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_32_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_33/outs')
sc = autoEstCont(sc)
LG70_NY_33.counts = adjustCounts(sc)
LG70_NY_33 <- CreateSeuratObject(counts = LG70_NY_33.counts, project = "LG70_NY_33", min.cells = 3, min.features = 200)
rm(LG70_NY_33.counts)
#vizualize QC metrics and filtering====
LG70_NY_33[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_33, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_33
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_33_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_33_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_33_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_33_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_34/outs')
sc = autoEstCont(sc)
LG70_NY_34.counts = adjustCounts(sc)
LG70_NY_34 <- CreateSeuratObject(counts = LG70_NY_34.counts, project = "LG70_NY_34", min.cells = 3, min.features = 200)
rm(LG70_NY_34.counts)
#vizualize QC metrics and filtering====
LG70_NY_34[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_34, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_34
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_34_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_34_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_34_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_34_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_36/outs')
sc = autoEstCont(sc)
LG70_NY_36.counts = adjustCounts(sc)
LG70_NY_36 <- CreateSeuratObject(counts = LG70_NY_36.counts, project = "LG70_NY_36", min.cells = 3, min.features = 200)
rm(LG70_NY_36.counts)
#vizualize QC metrics and filtering====
LG70_NY_36[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_36, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_36
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_36_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_36_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_36_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_36_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_37/outs')
sc = autoEstCont(sc)
LG70_NY_37.counts = adjustCounts(sc)
LG70_NY_37 <- CreateSeuratObject(counts = LG70_NY_37.counts, project = "LG70_NY_37", min.cells = 3, min.features = 200)
rm(LG70_NY_37.counts)
#vizualize QC metrics and filtering====
LG70_NY_37[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_37, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_37
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_37_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_37_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_37_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_37_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_38/outs')
sc = autoEstCont(sc)
LG70_NY_38.counts = adjustCounts(sc)
LG70_NY_38 <- CreateSeuratObject(counts = LG70_NY_38.counts, project = "LG70_NY_38", min.cells = 3, min.features = 200)
rm(LG70_NY_38.counts)
#vizualize QC metrics and filtering====
LG70_NY_38[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_38, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_38
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_38_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_38_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_38_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_38_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_39/outs')
sc = autoEstCont(sc)
LG70_NY_39.counts = adjustCounts(sc)
LG70_NY_39 <- CreateSeuratObject(counts = LG70_NY_39.counts, project = "LG70_NY_39", min.cells = 3, min.features = 200)
rm(LG70_NY_39.counts)
#vizualize QC metrics and filtering====
LG70_NY_39[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_39, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_39
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_39_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_39_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_39_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_39_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_40/outs')
sc = autoEstCont(sc)
LG70_NY_40.counts = adjustCounts(sc)
LG70_NY_40 <- CreateSeuratObject(counts = LG70_NY_40.counts, project = "LG70_NY_40", min.cells = 3, min.features = 200)
rm(LG70_NY_40.counts)
#vizualize QC metrics and filtering====
LG70_NY_40[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_40, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_40
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_40_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_40_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_40_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_40_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Tlr7/cellranger/LG70_NY_41/outs')
sc = autoEstCont(sc)
LG70_NY_41.counts = adjustCounts(sc)
LG70_NY_41 <- CreateSeuratObject(counts = LG70_NY_41.counts, project = "LG70_NY_41", min.cells = 3, min.features = 200)
rm(LG70_NY_41.counts)
#vizualize QC metrics and filtering====
LG70_NY_41[["percent.mt"]] <- PercentageFeatureSet(object = LG70_NY_41, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG70_NY_41
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG70_NY_41_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG70_NY_41_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG70_NY_41_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG70_NY_41_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
