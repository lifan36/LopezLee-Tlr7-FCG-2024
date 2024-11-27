##OL subclustering

LG93E4 <- readRDS("LG93E4_integrated_Annotation.rds")
Idents(LG93E4) <- "celltype"
LG93E4_OL <- subset(LG93E4, idents = "oligodendrocytes")

DefaultAssay(LG93E4_OL) <- 'integrated'
LG93E4_OL <- ScaleData(LG93E4_OL, verbose = FALSE)
LG93E4_OL <- RunPCA(LG93E4_OL, features = VariableFeatures(object = LG93E4_OL), verbose = FALSE)
ElbowPlot(LG93E4_OL)
LG93E4_OL <- FindNeighbors(LG93E4_OL, dims = 1:8)
LG93E4_OL <- FindClusters(LG93E4_OL, resolution = 0.2)
LG93E4_OL <- RunUMAP(LG93E4_OL, dims = 1:8)
DimPlot(LG93E4_OL, reduction = "umap", label = TRUE)

#renumber
n <- dim(table(LG93E4_OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG93E4_OL@active.ident <- plyr::mapvalues(x = LG93E4_OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
LG93E4_OL@active.ident <- factor(LG93E4_OL@active.ident, levels=1:n)

#save
saveRDS(LG93E4_OL, file = 'LG93E4_OL_clustered.rds')