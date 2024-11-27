#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_1/outs')
sc = autoEstCont(sc)
M_Ctrl_1.counts = adjustCounts(sc)
M_Ctrl_1 <- CreateSeuratObject(counts = M_Ctrl_1.counts, project = "Jax_1", min.cells = 3, min.features = 200)
M_Ctrl_1[["Condition"]] = c('M_Ctrl')
M_Ctrl_1[["Sample_Name"]] = c('M_Ctrl_1')
rm(M_Ctrl_1.counts)
#vizualize QC metrics and filtering====
M_Ctrl_1[["percent.mt"]] <- PercentageFeatureSet(object = M_Ctrl_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_Ctrl_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_Ctrl_1_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_Ctrl_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_Ctrl_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_Ctrl_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_2/outs')
sc = autoEstCont(sc)
M_Ctrl_2.counts = adjustCounts(sc)
M_Ctrl_2 <- CreateSeuratObject(counts = M_Ctrl_2.counts, project = "Jax_2", min.cells = 3, min.features = 200)
M_Ctrl_2[["Condition"]] = c('M_Ctrl')
M_Ctrl_2[["Sample_Name"]] = c('M_Ctrl_2')
rm(M_Ctrl_2.counts)
#vizualize QC metrics and filtering====
M_Ctrl_2[["percent.mt"]] <- PercentageFeatureSet(object = M_Ctrl_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_Ctrl_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_Ctrl_2_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_Ctrl_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_Ctrl_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_Ctrl_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_3/outs')
sc = autoEstCont(sc)
M_Ctrl_3.counts = adjustCounts(sc)
M_Ctrl_3 <- CreateSeuratObject(counts = M_Ctrl_3.counts, project = "Jax_3", min.cells = 3, min.features = 200)
M_Ctrl_3[["Condition"]] = c('M_Ctrl')
M_Ctrl_3[["Sample_Name"]] = c('M_Ctrl_3')
rm(M_Ctrl_3.counts)
#vizualize QC metrics and filtering====
M_Ctrl_3[["percent.mt"]] <- PercentageFeatureSet(object = M_Ctrl_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_Ctrl_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_Ctrl_3_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_Ctrl_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_Ctrl_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_Ctrl_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_14/outs')
sc = autoEstCont(sc)
M_3wk_1.counts = adjustCounts(sc)
M_3wk_1 <- CreateSeuratObject(counts = M_3wk_1.counts, project = "Jax_14", min.cells = 3, min.features = 200)
M_3wk_1[["Condition"]] = c('M_3wk')
M_3wk_1[["Sample_Name"]] = c('M_3wk_1')
rm(M_3wk_1.counts)
#vizualize QC metrics and filtering====
M_3wk_1[["percent.mt"]] <- PercentageFeatureSet(object = M_3wk_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_3wk_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_3wk_1_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_3wk_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_3wk_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_3wk_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_15/outs')
sc = autoEstCont(sc)
M_3wk_2.counts = adjustCounts(sc)
M_3wk_2 <- CreateSeuratObject(counts = M_3wk_2.counts, project = "Jax_15", min.cells = 3, min.features = 200)
M_3wk_2[["Condition"]] = c('M_3wk')
M_3wk_2[["Sample_Name"]] = c('M_3wk_2')
rm(M_3wk_2.counts)
#vizualize QC metrics and filtering====
M_3wk_2[["percent.mt"]] <- PercentageFeatureSet(object = M_3wk_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_3wk_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_3wk_2_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_3wk_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_3wk_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_3wk_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_16/outs')
sc = autoEstCont(sc)
M_3wk_3.counts = adjustCounts(sc)
M_3wk_3 <- CreateSeuratObject(counts = M_3wk_3.counts, project = "Jax_16", min.cells = 3, min.features = 200)
M_3wk_3[["Condition"]] = c('M_3wk')
M_3wk_3[["Sample_Name"]] = c('M_3wk_3')
rm(M_3wk_3.counts)
#vizualize QC metrics and filtering====
M_3wk_3[["percent.mt"]] <- PercentageFeatureSet(object = M_3wk_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- M_3wk_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("M_3wk_3_FeatureScatter.pdf", width=12, height=4)
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
pdf("M_3wk_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"M_3wk_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("M_3wk_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_30/outs')
sc = autoEstCont(sc)
F_Ctrl_1.counts = adjustCounts(sc)
F_Ctrl_1 <- CreateSeuratObject(counts = F_Ctrl_1.counts, project = "Jax_30", min.cells = 3, min.features = 200)
F_Ctrl_1[["Condition"]] = c('F_Ctrl')
F_Ctrl_1[["Sample_Name"]] = c('F_Ctrl_1')
rm(F_Ctrl_1.counts)
#vizualize QC metrics and filtering====
F_Ctrl_1[["percent.mt"]] <- PercentageFeatureSet(object = F_Ctrl_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_Ctrl_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_Ctrl_1_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_Ctrl_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_Ctrl_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_Ctrl_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_31/outs')
sc = autoEstCont(sc)
F_Ctrl_2.counts = adjustCounts(sc)
F_Ctrl_2 <- CreateSeuratObject(counts = F_Ctrl_2.counts, project = "Jax_31", min.cells = 3, min.features = 200)
F_Ctrl_2[["Condition"]] = c('F_Ctrl')
F_Ctrl_2[["Sample_Name"]] = c('F_Ctrl_2')
rm(F_Ctrl_2.counts)
#vizualize QC metrics and filtering====
F_Ctrl_2[["percent.mt"]] <- PercentageFeatureSet(object = F_Ctrl_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_Ctrl_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_Ctrl_2_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_Ctrl_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_Ctrl_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_Ctrl_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_32/outs')
sc = autoEstCont(sc)
F_Ctrl_3.counts = adjustCounts(sc)
F_Ctrl_3 <- CreateSeuratObject(counts = F_Ctrl_3.counts, project = "Jax_32", min.cells = 3, min.features = 200)
F_Ctrl_3[["Condition"]] = c('F_Ctrl')
F_Ctrl_3[["Sample_Name"]] = c('F_Ctrl_3')
rm(F_Ctrl_3.counts)
#vizualize QC metrics and filtering====
F_Ctrl_3[["percent.mt"]] <- PercentageFeatureSet(object = F_Ctrl_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_Ctrl_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_Ctrl_3_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_Ctrl_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_Ctrl_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_Ctrl_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_38/outs')
sc = autoEstCont(sc)
F_3wk_1.counts = adjustCounts(sc)
F_3wk_1 <- CreateSeuratObject(counts = F_3wk_1.counts, project = "Jax_38", min.cells = 3, min.features = 200)
F_3wk_1[["Condition"]] = c('F_3wk')
F_3wk_1[["Sample_Name"]] = c('F_3wk_1')
rm(F_3wk_1.counts)
#vizualize QC metrics and filtering====
F_3wk_1[["percent.mt"]] <- PercentageFeatureSet(object = F_3wk_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_3wk_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_3wk_1_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_3wk_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_3wk_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_3wk_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_36/outs')
sc = autoEstCont(sc)
F_3wk_2.counts = adjustCounts(sc)
F_3wk_2 <- CreateSeuratObject(counts = F_3wk_2.counts, project = "Jax_36", min.cells = 3, min.features = 200)
F_3wk_2[["Condition"]] = c('F_3wk')
F_3wk_2[["Sample_Name"]] = c('F_3wk_2')
rm(F_3wk_2.counts)
#vizualize QC metrics and filtering====
F_3wk_2[["percent.mt"]] <- PercentageFeatureSet(object = F_3wk_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_3wk_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_3wk_2_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_3wk_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_3wk_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_3wk_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/FCG_YoungCPZ/cellranger/Jax_37/outs')
sc = autoEstCont(sc)
F_3wk_3.counts = adjustCounts(sc)
F_3wk_3 <- CreateSeuratObject(counts = F_3wk_3.counts, project = "Jax_37", min.cells = 3, min.features = 200)
F_3wk_3[["Condition"]] = c('F_3wk')
F_3wk_3[["Sample_Name"]] = c('F_3wk_3')
rm(F_3wk_3.counts)
#vizualize QC metrics and filtering====
F_3wk_3[["percent.mt"]] <- PercentageFeatureSet(object = F_3wk_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- F_3wk_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("F_3wk_3_FeatureScatter.pdf", width=12, height=4)
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
pdf("F_3wk_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"F_3wk_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("F_3wk_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################










