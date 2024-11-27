##OL subclustering
YoungCPZ_OL <- subset(YoungCPZ, idents = "oligodendrocytes")

DefaultAssay(YoungCPZ_OL) <- 'integrated'
YoungCPZ_OL <- ScaleData(YoungCPZ_OL, verbose = FALSE)
YoungCPZ_OL <- RunPCA(YoungCPZ_OL, features = VariableFeatures(object = YoungCPZ_OL), verbose = FALSE)
ElbowPlot(YoungCPZ_OL)
YoungCPZ_OL <- FindNeighbors(YoungCPZ_OL, dims = 1:7)
YoungCPZ_OL <- FindClusters(YoungCPZ_OL, resolution = 0.2)
YoungCPZ_OL <- RunUMAP(YoungCPZ_OL, dims = 1:7)
DimPlot(YoungCPZ_OL, split.by = "Condition", reduction = "umap", label = TRUE, ncol=2)

#renumber to start at 1 instead of 0
n <- dim(table(YoungCPZ_OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
YoungCPZ_OL@active.ident <- plyr::mapvalues(x = YoungCPZ_OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
YoungCPZ_OL@active.ident <- factor(YoungCPZ_OL@active.ident, levels=1:n)

saveRDS(YoungCPZ_OL, "YoungCPZ_OL_7dims_res0.2.rds")