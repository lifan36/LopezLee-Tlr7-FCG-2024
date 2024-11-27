#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4_male/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_121/outs')
sc = autoEstCont(sc)
LG31E3_R47H_121.counts = adjustCounts(sc)
LG31E3_R47H_121 <- CreateSeuratObject(counts = LG31E3_R47H_121.counts, project = "LG31E3_R47H_121", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_121.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_121[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_121, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_121
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_121_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_121_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_121_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_121_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_122/outs')
sc = autoEstCont(sc)
LG31E3_R47H_122.counts = adjustCounts(sc)
LG31E3_R47H_122 <- CreateSeuratObject(counts = LG31E3_R47H_122.counts, project = "LG31E3_R47H_122", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_122.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_122[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_122, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_122
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_122_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_122_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_122_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_122_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_22/outs')
sc = autoEstCont(sc)
LG31E3_R47H_22.counts = adjustCounts(sc)
LG31E3_R47H_22 <- CreateSeuratObject(counts = LG31E3_R47H_22.counts, project = "LG31E3_R47H_22", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_22.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_22[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_22, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_22
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_22_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_22_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_22_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_22_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_23/outs')
sc = autoEstCont(sc)
LG31E3_R47H_23.counts = adjustCounts(sc)
LG31E3_R47H_23 <- CreateSeuratObject(counts = LG31E3_R47H_23.counts, project = "LG31E3_R47H_23", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_23.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_23[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_23, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_23
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_23_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_23_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_23_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_23_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_28/outs')
sc = autoEstCont(sc)
LG31E3_R47H_28.counts = adjustCounts(sc)
LG31E3_R47H_28 <- CreateSeuratObject(counts = LG31E3_R47H_28.counts, project = "LG31E3_R47H_28", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_28.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_28[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_28, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_28
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_28_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_28_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_28_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_28_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_30/outs')
sc = autoEstCont(sc)
LG31E3_R47H_30.counts = adjustCounts(sc)
LG31E3_R47H_30 <- CreateSeuratObject(counts = LG31E3_R47H_30.counts, project = "LG31E3_R47H_30", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_30.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_30[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_30, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_30
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_30_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_30_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_30_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_30_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_45/outs')
sc = autoEstCont(sc)
LG31E3_R47H_45.counts = adjustCounts(sc)
LG31E3_R47H_45 <- CreateSeuratObject(counts = LG31E3_R47H_45.counts, project = "LG31E3_R47H_45", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_45.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_45[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_45, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_45
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_45_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_45_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_45_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_45_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_93/outs')
sc = autoEstCont(sc)
LG31E3_R47H_93.counts = adjustCounts(sc)
LG31E3_R47H_93 <- CreateSeuratObject(counts = LG31E3_R47H_93.counts, project = "LG31E3_R47H_93", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_93.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_93[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_93, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_93
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_93_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_93_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_93_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_93_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
##############################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_R47H_95/outs')
sc = autoEstCont(sc)
LG31E3_R47H_95.counts = adjustCounts(sc)
LG31E3_R47H_95 <- CreateSeuratObject(counts = LG31E3_R47H_95.counts, project = "LG31E3_R47H_95", min.cells = 3, min.features = 200)
rm(LG31E3_R47H_95.counts)
#vizualize QC metrics and filtering====
LG31E3_R47H_95[["percent.mt"]] <- PercentageFeatureSet(object = LG31E3_R47H_95, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E3_R47H_95
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E3_R47H_95_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E3_R47H_95_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E3_R47H_95_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E3_R47H_95_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_101/outs')
sc = autoEstCont(sc)
LG31E4_R47H_101.counts = adjustCounts(sc)
LG31E4_R47H_101 <- CreateSeuratObject(counts = LG31E4_R47H_101.counts, project = "LG31E4_R47H_101", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_101.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_101[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_101, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_101
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_101_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_101_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_101_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_101_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_32/outs')
sc = autoEstCont(sc)
LG31E4_R47H_32.counts = adjustCounts(sc)
LG31E4_R47H_32 <- CreateSeuratObject(counts = LG31E4_R47H_32.counts, project = "LG31E4_R47H_32", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_32.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_32[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_32, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_32
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_32_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_32_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_32_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_32_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_60/outs')
sc = autoEstCont(sc)
LG31E4_R47H_60.counts = adjustCounts(sc)
LG31E4_R47H_60 <- CreateSeuratObject(counts = LG31E4_R47H_60.counts, project = "LG31E4_R47H_60", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_60.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_60[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_60, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_60
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_60_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_60_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_60_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_60_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_61/outs')
sc = autoEstCont(sc)
LG31E4_R47H_61.counts = adjustCounts(sc)
LG31E4_R47H_61 <- CreateSeuratObject(counts = LG31E4_R47H_61.counts, project = "LG31E4_R47H_61", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_61.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_61[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_61, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_61
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_61_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_61_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_61_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_61_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_68/outs')
sc = autoEstCont(sc)
LG31E4_R47H_68.counts = adjustCounts(sc)
LG31E4_R47H_68 <- CreateSeuratObject(counts = LG31E4_R47H_68.counts, project = "LG31E4_R47H_68", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_68.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_68[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_68, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_68
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_68_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_68_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_68_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_68_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_78/outs')
sc = autoEstCont(sc)
LG31E4_R47H_78.counts = adjustCounts(sc)
LG31E4_R47H_78 <- CreateSeuratObject(counts = LG31E4_R47H_78.counts, project = "LG31E4_R47H_78", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_78.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_78[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_78, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_78
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_78_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_78_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_78_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_78_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_80/outs')
sc = autoEstCont(sc)
LG31E4_R47H_80.counts = adjustCounts(sc)
LG31E4_R47H_80 <- CreateSeuratObject(counts = LG31E4_R47H_80.counts, project = "LG31E4_R47H_80", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_80.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_80[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_80, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_80
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_80_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_80_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_80_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_80_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_82/outs')
sc = autoEstCont(sc)
LG31E4_R47H_82.counts = adjustCounts(sc)
LG31E4_R47H_82 <- CreateSeuratObject(counts = LG31E4_R47H_82.counts, project = "LG31E4_R47H_82", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_82.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_82[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_82, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_82
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_82_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_82_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_82_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_82_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_R47H_94/outs')
sc = autoEstCont(sc)
LG31E4_R47H_94.counts = adjustCounts(sc)
LG31E4_R47H_94 <- CreateSeuratObject(counts = LG31E4_R47H_94.counts, project = "LG31E4_R47H_94", min.cells = 3, min.features = 200)
rm(LG31E4_R47H_94.counts)
#vizualize QC metrics and filtering====
LG31E4_R47H_94[["percent.mt"]] <- PercentageFeatureSet(object = LG31E4_R47H_94, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG31E4_R47H_94
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG31E4_R47H_94_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG31E4_R47H_94_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG31E4_R47H_94_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG31E4_R47H_94_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
