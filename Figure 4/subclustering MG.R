##MG subclustering

LG93E4 <- readRDS("LG93E4_integrated_Annotation.rds")
Idents(LG93E4) <- "celltype"
LG93E4_MG <- subset(LG93E4, idents = "microglia")

DefaultAssay(LG93E4_MG) <- 'integrated'
LG93E4_MG <- ScaleData(LG93E4_MG, verbose = FALSE)
LG93E4_MG <- RunPCA(LG93E4_MG, features = VariableFeatures(object = LG93E4_MG), verbose = FALSE)
ElbowPlot(LG93E4_MG)
LG93E4_MG <- FindNeighbors(LG93E4_MG, dims = 1:8)
LG93E4_MG <- FindClusters(LG93E4_MG, resolution = 0.2)
LG93E4_MG <- RunUMAP(LG93E4_MG, dims = 1:8)
DimPlot(LG93E4_MG, reduction = "umap", label = TRUE)

#renumber
n <- dim(table(LG93E4_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG93E4_MG@active.ident <- plyr::mapvalues(x = LG93E4_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
LG93E4_MG@active.ident <- factor(LG93E4_MG@active.ident, levels=1:n)
saveRDS(LG93E4_MG, file = 'LG93E4_MG_clustered.rds')
