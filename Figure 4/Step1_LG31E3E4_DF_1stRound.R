#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG31E3E4/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_17/outs')
sc = autoEstCont(sc)
E3_NTG_WT_1.counts = adjustCounts(sc)
E3_NTG_WT_1 <- CreateSeuratObject(counts = E3_NTG_WT_1.counts, project = "LG31E3_17", min.cells = 3, min.features = 200)
E3_NTG_WT_1[["Condition"]] = c('E3_NTG_WT')
E3_NTG_WT_1[["Sample_Name"]] = c('E3_NTG_WT_1')
rm(E3_NTG_WT_1.counts)
#vizualize QC metrics and filtering====
E3_NTG_WT_1[["percent.mt"]] <- PercentageFeatureSet(object = E3_NTG_WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_NTG_WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_NTG_WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_NTG_WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_NTG_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_NTG_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_1/outs')
sc = autoEstCont(sc)
E3_P301S_WT_1.counts = adjustCounts(sc)
E3_P301S_WT_1 <- CreateSeuratObject(counts = E3_P301S_WT_1.counts, project = "LG31E3_1", min.cells = 3, min.features = 200)
E3_P301S_WT_1[["Condition"]] = c('E3_P301S_WT')
E3_P301S_WT_1[["Sample_Name"]] = c('E3_P301S_WT_1')
rm(E3_P301S_WT_1.counts)
#vizualize QC metrics and filtering====
E3_P301S_WT_1[["percent.mt"]] <- PercentageFeatureSet(object = E3_P301S_WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_P301S_WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_P301S_WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_P301S_WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_P301S_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_18/outs')
sc = autoEstCont(sc)
E3_P301S_WT_2.counts = adjustCounts(sc)
E3_P301S_WT_2 <- CreateSeuratObject(counts = E3_P301S_WT_2.counts, project = "LG31E3_18", min.cells = 3, min.features = 200)
E3_P301S_WT_2[["Condition"]] = c('E3_P301S_WT')
E3_P301S_WT_2[["Sample_Name"]] = c('E3_P301S_WT_2')
rm(E3_P301S_WT_2.counts)
#vizualize QC metrics and filtering====
E3_P301S_WT_2[["percent.mt"]] <- PercentageFeatureSet(object = E3_P301S_WT_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_P301S_WT_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_P301S_WT_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_P301S_WT_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_P301S_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_4/outs')
sc = autoEstCont(sc)
E3_P301S_R47H_1.counts = adjustCounts(sc)
E3_P301S_R47H_1 <- CreateSeuratObject(counts = E3_P301S_R47H_1.counts, project = "LG31E3_4", min.cells = 3, min.features = 200)
E3_P301S_R47H_1[["Condition"]] = c('E3_P301S_R47H')
E3_P301S_R47H_1[["Sample_Name"]] = c('E3_P301S_R47H_1')
rm(E3_P301S_R47H_1.counts)
#vizualize QC metrics and filtering====
E3_P301S_R47H_1[["percent.mt"]] <- PercentageFeatureSet(object = E3_P301S_R47H_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_P301S_R47H_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_P301S_R47H_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_P301S_R47H_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_P301S_R47H_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_11/outs')
sc = autoEstCont(sc)
E3_P301S_R47H_2.counts = adjustCounts(sc)
E3_P301S_R47H_2 <- CreateSeuratObject(counts = E3_P301S_R47H_2.counts, project = "LG31E3_11", min.cells = 3, min.features = 200)
E3_P301S_R47H_2[["Condition"]] = c('E3_P301S_R47H')
E3_P301S_R47H_2[["Sample_Name"]] = c('E3_P301S_R47H_2')
rm(E3_P301S_R47H_2.counts)
#vizualize QC metrics and filtering====
E3_P301S_R47H_2[["percent.mt"]] <- PercentageFeatureSet(object = E3_P301S_R47H_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_P301S_R47H_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_P301S_R47H_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_P301S_R47H_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_P301S_R47H_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E3_12/outs')
sc = autoEstCont(sc)
E3_P301S_R47H_3.counts = adjustCounts(sc)
E3_P301S_R47H_3 <- CreateSeuratObject(counts = E3_P301S_R47H_3.counts, project = "LG31E3_12", min.cells = 3, min.features = 200)
E3_P301S_R47H_3[["Condition"]] = c('E3_P301S_R47H')
E3_P301S_R47H_3[["Sample_Name"]] = c('E3_P301S_R47H_3')
rm(E3_P301S_R47H_3.counts)
#vizualize QC metrics and filtering====
E3_P301S_R47H_3[["percent.mt"]] <- PercentageFeatureSet(object = E3_P301S_R47H_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E3_P301S_R47H_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E3_P301S_R47H_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E3_P301S_R47H_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E3_P301S_R47H_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E3_P301S_R47H_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_25/outs')
sc = autoEstCont(sc)
E4_NTG_WT_1.counts = adjustCounts(sc)
E4_NTG_WT_1 <- CreateSeuratObject(counts = E4_NTG_WT_1.counts, project = "LG31E4_25", min.cells = 3, min.features = 200)
E4_NTG_WT_1[["Condition"]] = c('E4_NTG_WT')
E4_NTG_WT_1[["Sample_Name"]] = c('E4_NTG_WT_1')
rm(E4_NTG_WT_1.counts)
#vizualize QC metrics and filtering====
E4_NTG_WT_1[["percent.mt"]] <- PercentageFeatureSet(object = E4_NTG_WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_NTG_WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_NTG_WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_NTG_WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_NTG_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_NTG_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_27/outs')
sc = autoEstCont(sc)
E4_NTG_WT_2.counts = adjustCounts(sc)
E4_NTG_WT_2 <- CreateSeuratObject(counts = E4_NTG_WT_2.counts, project = "LG31E4_27", min.cells = 3, min.features = 200)
E4_NTG_WT_2[["Condition"]] = c('E4_NTG_WT')
E4_NTG_WT_2[["Sample_Name"]] = c('E4_NTG_WT_2')
rm(E4_NTG_WT_2.counts)
#vizualize QC metrics and filtering====
E4_NTG_WT_2[["percent.mt"]] <- PercentageFeatureSet(object = E4_NTG_WT_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_NTG_WT_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_NTG_WT_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_NTG_WT_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_NTG_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_NTG_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_21/outs')
sc = autoEstCont(sc)
E4_P301S_WT_1.counts = adjustCounts(sc)
E4_P301S_WT_1 <- CreateSeuratObject(counts = E4_P301S_WT_1.counts, project = "LG31E4_21", min.cells = 3, min.features = 200)
E4_P301S_WT_1[["Condition"]] = c('E4_P301S_WT')
E4_P301S_WT_1[["Sample_Name"]] = c('E4_P301S_WT_1')
rm(E4_P301S_WT_1.counts)
#vizualize QC metrics and filtering====
E4_P301S_WT_1[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_28/outs')
sc = autoEstCont(sc)
E4_P301S_WT_2.counts = adjustCounts(sc)
E4_P301S_WT_2 <- CreateSeuratObject(counts = E4_P301S_WT_2.counts, project = "LG31E4_28", min.cells = 3, min.features = 200)
E4_P301S_WT_2[["Condition"]] = c('E4_P301S_WT')
E4_P301S_WT_2[["Sample_Name"]] = c('E4_P301S_WT_2')
rm(E4_P301S_WT_2.counts)
#vizualize QC metrics and filtering====
E4_P301S_WT_2[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_WT_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_WT_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_WT_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_WT_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_54/outs')
sc = autoEstCont(sc)
E4_P301S_WT_3.counts = adjustCounts(sc)
E4_P301S_WT_3 <- CreateSeuratObject(counts = E4_P301S_WT_3.counts, project = "LG31E4_54", min.cells = 3, min.features = 200)
E4_P301S_WT_3[["Condition"]] = c('E4_P301S_WT')
E4_P301S_WT_3[["Sample_Name"]] = c('E4_P301S_WT_3')
rm(E4_P301S_WT_3.counts)
#vizualize QC metrics and filtering====
E4_P301S_WT_3[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_WT_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_WT_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_WT_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_WT_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_WT_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_73/outs')
sc = autoEstCont(sc)
E4_P301S_WT_4.counts = adjustCounts(sc)
E4_P301S_WT_4 <- CreateSeuratObject(counts = E4_P301S_WT_4.counts, project = "LG31E4_73", min.cells = 3, min.features = 200)
E4_P301S_WT_4[["Condition"]] = c('E4_P301S_WT')
E4_P301S_WT_4[["Sample_Name"]] = c('E4_P301S_WT_4')
rm(E4_P301S_WT_4.counts)
#vizualize QC metrics and filtering====
E4_P301S_WT_4[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_WT_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_WT_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_WT_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_WT_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_WT_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_WT_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_22/outs')
sc = autoEstCont(sc)
E4_P301S_R47H_1.counts = adjustCounts(sc)
E4_P301S_R47H_1 <- CreateSeuratObject(counts = E4_P301S_R47H_1.counts, project = "LG31E4_22", min.cells = 3, min.features = 200)
E4_P301S_R47H_1[["Condition"]] = c('E4_P301S_R47H')
E4_P301S_R47H_1[["Sample_Name"]] = c('E4_P301S_R47H_1')
rm(E4_P301S_R47H_1.counts)
#vizualize QC metrics and filtering====
E4_P301S_R47H_1[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_R47H_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_R47H_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_R47H_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_R47H_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_R47H_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_33/outs')
sc = autoEstCont(sc)
E4_P301S_R47H_2.counts = adjustCounts(sc)
E4_P301S_R47H_2 <- CreateSeuratObject(counts = E4_P301S_R47H_2.counts, project = "LG31E4_33", min.cells = 3, min.features = 200)
E4_P301S_R47H_2[["Condition"]] = c('E4_P301S_R47H')
E4_P301S_R47H_2[["Sample_Name"]] = c('E4_P301S_R47H_2')
rm(E4_P301S_R47H_2.counts)
#vizualize QC metrics and filtering====
E4_P301S_R47H_2[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_R47H_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_R47H_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_R47H_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_R47H_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_R47H_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG31E3E4/cellranger/LG31E4_45/outs')
sc = autoEstCont(sc)
E4_P301S_R47H_3.counts = adjustCounts(sc)
E4_P301S_R47H_3 <- CreateSeuratObject(counts = E4_P301S_R47H_3.counts, project = "LG31E4_45", min.cells = 3, min.features = 200)
E4_P301S_R47H_3[["Condition"]] = c('E4_P301S_R47H')
E4_P301S_R47H_3[["Sample_Name"]] = c('E4_P301S_R47H_3')
rm(E4_P301S_R47H_3.counts)
#vizualize QC metrics and filtering====
E4_P301S_R47H_3[["percent.mt"]] <- PercentageFeatureSet(object = E4_P301S_R47H_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- E4_P301S_R47H_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("E4_P301S_R47H_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("E4_P301S_R47H_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"E4_P301S_R47H_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("E4_P301S_R47H_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()









