#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG93E4/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_32/outs')
sc = autoEstCont(sc)
XXO_NTG_E4_1.counts = adjustCounts(sc)
XXO_NTG_E4_1 <- CreateSeuratObject(counts = XXO_NTG_E4_1.counts, project = "LG93E4_32", min.cells = 3, min.features = 200)
XXO_NTG_E4_1[["Condition"]] = c('XXO_NTG_E4')
XXO_NTG_E4_1[["Sample_Name"]] = c('XXO_NTG_E4_1')
rm(XXO_NTG_E4_1.counts)
#vizualize QC metrics and filtering====
XXO_NTG_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XXO_NTG_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_NTG_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_NTG_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_NTG_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_37/outs')
sc = autoEstCont(sc)
XXO_NTG_E4_2.counts = adjustCounts(sc)
XXO_NTG_E4_2 <- CreateSeuratObject(counts = XXO_NTG_E4_2.counts, project = "LG93E4_37", min.cells = 3, min.features = 200)
XXO_NTG_E4_2[["Condition"]] = c('XXO_NTG_E4')
XXO_NTG_E4_2[["Sample_Name"]] = c('XXO_NTG_E4_2')
rm(XXO_NTG_E4_2.counts)
#vizualize QC metrics and filtering====
XXO_NTG_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XXO_NTG_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_NTG_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_NTG_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_NTG_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_38/outs')
sc = autoEstCont(sc)
XXO_NTG_E4_3.counts = adjustCounts(sc)
XXO_NTG_E4_3 <- CreateSeuratObject(counts = XXO_NTG_E4_3.counts, project = "LG93E4_38", min.cells = 3, min.features = 200)
XXO_NTG_E4_3[["Condition"]] = c('XXO_NTG_E4')
XXO_NTG_E4_3[["Sample_Name"]] = c('XXO_NTG_E4_3')
rm(XXO_NTG_E4_3.counts)
#vizualize QC metrics and filtering====
XXO_NTG_E4_3[["percent.mt"]] <- PercentageFeatureSet(object = XXO_NTG_E4_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_NTG_E4_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_NTG_E4_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_NTG_E4_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_22/outs')
sc = autoEstCont(sc)
XXO_P301S_E4_1.counts = adjustCounts(sc)
XXO_P301S_E4_1 <- CreateSeuratObject(counts = XXO_P301S_E4_1.counts, project = "LG93E4_22", min.cells = 3, min.features = 200)
XXO_P301S_E4_1[["Condition"]] = c('XXO_P301S_E4')
XXO_P301S_E4_1[["Sample_Name"]] = c('XXO_P301S_E4_1')
rm(XXO_P301S_E4_1.counts)
#vizualize QC metrics and filtering====
XXO_P301S_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XXO_P301S_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_P301S_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_P301S_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_P301S_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_24/outs')
sc = autoEstCont(sc)
XXO_P301S_E4_2.counts = adjustCounts(sc)
XXO_P301S_E4_2 <- CreateSeuratObject(counts = XXO_P301S_E4_2.counts, project = "LG93E4_24", min.cells = 3, min.features = 200)
XXO_P301S_E4_2[["Condition"]] = c('XXO_P301S_E4')
XXO_P301S_E4_2[["Sample_Name"]] = c('XXO_P301S_E4_2')
rm(XXO_P301S_E4_2.counts)
#vizualize QC metrics and filtering====
XXO_P301S_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XXO_P301S_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_P301S_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_P301S_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_P301S_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_23/outs')
sc = autoEstCont(sc)
XYO_NTG_E4_1.counts = adjustCounts(sc)
XYO_NTG_E4_1 <- CreateSeuratObject(counts = XYO_NTG_E4_1.counts, project = "LG93E4_23", min.cells = 3, min.features = 200)
XYO_NTG_E4_1[["Condition"]] = c('XYO_NTG_E4')
XYO_NTG_E4_1[["Sample_Name"]] = c('XYO_NTG_E4_1')
rm(XYO_NTG_E4_1.counts)
#vizualize QC metrics and filtering====
XYO_NTG_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XYO_NTG_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_NTG_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_NTG_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_NTG_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_3/outs')
sc = autoEstCont(sc)
XYO_NTG_E4_2.counts = adjustCounts(sc)
XYO_NTG_E4_2 <- CreateSeuratObject(counts = XYO_NTG_E4_2.counts, project = "LG93E4_3", min.cells = 3, min.features = 200)
XYO_NTG_E4_2[["Condition"]] = c('XYO_NTG_E4')
XYO_NTG_E4_2[["Sample_Name"]] = c('XYO_NTG_E4_2')
rm(XYO_NTG_E4_2.counts)
#vizualize QC metrics and filtering====
XYO_NTG_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XYO_NTG_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_NTG_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_NTG_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_NTG_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_41/outs')
sc = autoEstCont(sc)
XYO_NTG_E4_3.counts = adjustCounts(sc)
XYO_NTG_E4_3 <- CreateSeuratObject(counts = XYO_NTG_E4_3.counts, project = "LG93E4_41", min.cells = 3, min.features = 200)
XYO_NTG_E4_3[["Condition"]] = c('XYO_NTG_E4')
XYO_NTG_E4_3[["Sample_Name"]] = c('XYO_NTG_E4_3')
rm(XYO_NTG_E4_3.counts)
#vizualize QC metrics and filtering====
XYO_NTG_E4_3[["percent.mt"]] <- PercentageFeatureSet(object = XYO_NTG_E4_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_NTG_E4_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_NTG_E4_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_NTG_E4_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_39/outs')
sc = autoEstCont(sc)
XYO_P301S_E4_2.counts = adjustCounts(sc)
XYO_P301S_E4_2 <- CreateSeuratObject(counts = XYO_P301S_E4_2.counts, project = "LG93E4_39", min.cells = 3, min.features = 200)
XYO_P301S_E4_2[["Condition"]] = c('XYO_P301S_E4')
XYO_P301S_E4_2[["Sample_Name"]] = c('XYO_P301S_E4_2')
rm(XYO_P301S_E4_2.counts)
#vizualize QC metrics and filtering====
XYO_P301S_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XYO_P301S_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_P301S_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_P301S_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_P301S_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_42/outs')
sc = autoEstCont(sc)
XYO_P301S_E4_3.counts = adjustCounts(sc)
XYO_P301S_E4_3 <- CreateSeuratObject(counts = XYO_P301S_E4_3.counts, project = "LG93E4_42", min.cells = 3, min.features = 200)
XYO_P301S_E4_3[["Condition"]] = c('XYO_P301S_E4')
XYO_P301S_E4_3[["Sample_Name"]] = c('XYO_P301S_E4_3')
rm(XYO_P301S_E4_3.counts)
#vizualize QC metrics and filtering====
XYO_P301S_E4_3[["percent.mt"]] <- PercentageFeatureSet(object = XYO_P301S_E4_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_P301S_E4_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_P301S_E4_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_P301S_E4_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_P301S_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_30/outs')
sc = autoEstCont(sc)
XXT_NTG_E4_1.counts = adjustCounts(sc)
XXT_NTG_E4_1 <- CreateSeuratObject(counts = XXT_NTG_E4_1.counts, project = "LG93E4_30", min.cells = 3, min.features = 200)
XXT_NTG_E4_1[["Condition"]] = c('XXT_NTG_E4')
XXT_NTG_E4_1[["Sample_Name"]] = c('XXT_NTG_E4_1')
rm(XXT_NTG_E4_1.counts)
#vizualize QC metrics and filtering====
XXT_NTG_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XXT_NTG_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_NTG_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_NTG_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_NTG_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_31/outs')
sc = autoEstCont(sc)
XXT_NTG_E4_2.counts = adjustCounts(sc)
XXT_NTG_E4_2 <- CreateSeuratObject(counts = XXT_NTG_E4_2.counts, project = "LG93E4_31", min.cells = 3, min.features = 200)
XXT_NTG_E4_2[["Condition"]] = c('XXT_NTG_E4')
XXT_NTG_E4_2[["Sample_Name"]] = c('XXT_NTG_E4_2')
rm(XXT_NTG_E4_2.counts)
#vizualize QC metrics and filtering====
XXT_NTG_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XXT_NTG_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_NTG_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_NTG_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_NTG_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_47/outs')
sc = autoEstCont(sc)
XXT_NTG_E4_3.counts = adjustCounts(sc)
XXT_NTG_E4_3 <- CreateSeuratObject(counts = XXT_NTG_E4_3.counts, project = "LG93E4_47", min.cells = 3, min.features = 200)
XXT_NTG_E4_3[["Condition"]] = c('XXT_NTG_E4')
XXT_NTG_E4_3[["Sample_Name"]] = c('XXT_NTG_E4_3')
rm(XXT_NTG_E4_3.counts)
#vizualize QC metrics and filtering====
XXT_NTG_E4_3[["percent.mt"]] <- PercentageFeatureSet(object = XXT_NTG_E4_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_NTG_E4_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_NTG_E4_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_NTG_E4_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_45/outs')
sc = autoEstCont(sc)
XXT_P301S_E4_1.counts = adjustCounts(sc)
XXT_P301S_E4_1 <- CreateSeuratObject(counts = XXT_P301S_E4_1.counts, project = "LG93E4_45", min.cells = 3, min.features = 200)
XXT_P301S_E4_1[["Condition"]] = c('XXT_P301S_E4')
XXT_P301S_E4_1[["Sample_Name"]] = c('XXT_P301S_E4_1')
rm(XXT_P301S_E4_1.counts)
#vizualize QC metrics and filtering====
XXT_P301S_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XXT_P301S_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_P301S_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_P301S_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_P301S_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_48/outs')
sc = autoEstCont(sc)
XXT_P301S_E4_2.counts = adjustCounts(sc)
XXT_P301S_E4_2 <- CreateSeuratObject(counts = XXT_P301S_E4_2.counts, project = "LG93E4_48", min.cells = 3, min.features = 200)
XXT_P301S_E4_2[["Condition"]] = c('XXT_P301S_E4')
XXT_P301S_E4_2[["Sample_Name"]] = c('XXT_P301S_E4_2')
rm(XXT_P301S_E4_2.counts)
#vizualize QC metrics and filtering====
XXT_P301S_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XXT_P301S_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_P301S_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_P301S_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_P301S_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_29/outs')
sc = autoEstCont(sc)
XYT_NTG_E4_1.counts = adjustCounts(sc)
XYT_NTG_E4_1 <- CreateSeuratObject(counts = XYT_NTG_E4_1.counts, project = "LG93E4_29", min.cells = 3, min.features = 200)
XYT_NTG_E4_1[["Condition"]] = c('XYT_NTG_E4')
XYT_NTG_E4_1[["Sample_Name"]] = c('XYT_NTG_E4_1')
rm(XYT_NTG_E4_1.counts)
#vizualize QC metrics and filtering====
XYT_NTG_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XYT_NTG_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_NTG_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_NTG_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_NTG_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_NTG_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_79/outs')
sc = autoEstCont(sc)
XYT_NTG_E4_2.counts = adjustCounts(sc)
XYT_NTG_E4_2 <- CreateSeuratObject(counts = XYT_NTG_E4_2.counts, project = "LG93E4_79", min.cells = 3, min.features = 200)
XYT_NTG_E4_2[["Condition"]] = c('XYT_NTG_E4')
XYT_NTG_E4_2[["Sample_Name"]] = c('XYT_NTG_E4_2')
rm(XYT_NTG_E4_2.counts)
#vizualize QC metrics and filtering====
XYT_NTG_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XYT_NTG_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_NTG_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_NTG_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_NTG_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_NTG_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_92/outs')
sc = autoEstCont(sc)
XYT_NTG_E4_3.counts = adjustCounts(sc)
XYT_NTG_E4_3 <- CreateSeuratObject(counts = XYT_NTG_E4_3.counts, project = "LG93E4_92", min.cells = 3, min.features = 200)
XYT_NTG_E4_3[["Condition"]] = c('XYT_NTG_E4')
XYT_NTG_E4_3[["Sample_Name"]] = c('XYT_NTG_E4_3')
rm(XYT_NTG_E4_3.counts)
#vizualize QC metrics and filtering====
XYT_NTG_E4_3[["percent.mt"]] <- PercentageFeatureSet(object = XYT_NTG_E4_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_NTG_E4_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_NTG_E4_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_NTG_E4_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_NTG_E4_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_NTG_E4_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_18/outs')
sc = autoEstCont(sc)
XYT_P301S_E4_1.counts = adjustCounts(sc)
XYT_P301S_E4_1 <- CreateSeuratObject(counts = XYT_P301S_E4_1.counts, project = "LG93E4_18", min.cells = 3, min.features = 200)
XYT_P301S_E4_1[["Condition"]] = c('XYT_P301S_E4')
XYT_P301S_E4_1[["Sample_Name"]] = c('XYT_P301S_E4_1')
rm(XYT_P301S_E4_1.counts)
#vizualize QC metrics and filtering====
XYT_P301S_E4_1[["percent.mt"]] <- PercentageFeatureSet(object = XYT_P301S_E4_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_P301S_E4_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_P301S_E4_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_P301S_E4_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_P301S_E4_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_28/outs')
sc = autoEstCont(sc)
XYT_P301S_E4_2.counts = adjustCounts(sc)
XYT_P301S_E4_2 <- CreateSeuratObject(counts = XYT_P301S_E4_2.counts, project = "LG93E4_28", min.cells = 3, min.features = 200)
XYT_P301S_E4_2[["Condition"]] = c('XYT_P301S_E4')
XYT_P301S_E4_2[["Sample_Name"]] = c('XYT_P301S_E4_2')
rm(XYT_P301S_E4_2.counts)
#vizualize QC metrics and filtering====
XYT_P301S_E4_2[["percent.mt"]] <- PercentageFeatureSet(object = XYT_P301S_E4_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_P301S_E4_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_P301S_E4_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_P301S_E4_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_P301S_E4_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_84/outs')
sc = autoEstCont(sc)
XYO_P301S_E4_4.counts = adjustCounts(sc)
XYO_P301S_E4_4 <- CreateSeuratObject(counts = XYO_P301S_E4_4.counts, project = "LG93E4_84", min.cells = 3, min.features = 200)
XYO_P301S_E4_4[["Condition"]] = c('XYO_P301S_E4')
XYO_P301S_E4_4[["Sample_Name"]] = c('XYO_P301S_E4_4')
rm(XYO_P301S_E4_4.counts)
#vizualize QC metrics and filtering====
XYO_P301S_E4_4[["percent.mt"]] <- PercentageFeatureSet(object = XYO_P301S_E4_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYO_P301S_E4_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYO_P301S_E4_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYO_P301S_E4_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYO_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYO_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_36/outs')
sc = autoEstCont(sc)
XXO_P301S_E4_4.counts = adjustCounts(sc)
XXO_P301S_E4_4 <- CreateSeuratObject(counts = XXO_P301S_E4_4.counts, project = "LG93E4_36", min.cells = 3, min.features = 200)
XXO_P301S_E4_4[["Condition"]] = c('XXO_P301S_E4')
XXO_P301S_E4_4[["Sample_Name"]] = c('XXO_P301S_E4_4')
rm(XXO_P301S_E4_4.counts)
#vizualize QC metrics and filtering====
XXO_P301S_E4_4[["percent.mt"]] <- PercentageFeatureSet(object = XXO_P301S_E4_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXO_P301S_E4_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXO_P301S_E4_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXO_P301S_E4_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXO_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXO_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_81/outs')
sc = autoEstCont(sc)
XYT_P301S_E4_4.counts = adjustCounts(sc)
XYT_P301S_E4_4 <- CreateSeuratObject(counts = XYT_P301S_E4_4.counts, project = "LG93E4_81", min.cells = 3, min.features = 200)
XYT_P301S_E4_4[["Condition"]] = c('XYT_P301S_E4')
XYT_P301S_E4_4[["Sample_Name"]] = c('XYT_P301S_E4_4')
rm(XYT_P301S_E4_4.counts)
#vizualize QC metrics and filtering====
XYT_P301S_E4_4[["percent.mt"]] <- PercentageFeatureSet(object = XYT_P301S_E4_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XYT_P301S_E4_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XYT_P301S_E4_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XYT_P301S_E4_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XYT_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XYT_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG93E4/cellranger/LG93E4_80/outs')
sc = autoEstCont(sc)
XXT_P301S_E4_4.counts = adjustCounts(sc)
XXT_P301S_E4_4 <- CreateSeuratObject(counts = XXT_P301S_E4_4.counts, project = "LG93E4_80", min.cells = 3, min.features = 200)
XXT_P301S_E4_4[["Condition"]] = c('XXT_P301S_E4')
XXT_P301S_E4_4[["Sample_Name"]] = c('XXT_P301S_E4_4')
rm(XXT_P301S_E4_4.counts)
#vizualize QC metrics and filtering====
XXT_P301S_E4_4[["percent.mt"]] <- PercentageFeatureSet(object = XXT_P301S_E4_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- XXT_P301S_E4_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("XXT_P301S_E4_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("XXT_P301S_E4_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"XXT_P301S_E4_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("XXT_P301S_E4_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################










