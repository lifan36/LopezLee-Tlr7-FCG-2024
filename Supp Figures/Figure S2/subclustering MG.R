##subcluster MG
#isolate microglia
YoungCPZ_MG <- subset(YoungCPZ, idents = "microglia")
saveRDS(YoungCPZ_MG, "MG.rds")

DimPlot(YoungCPZ_MG, reduction = "umap", label = TRUE)

#need to subcluster
DefaultAssay(YoungCPZ_MG) <- 'integrated'
YoungCPZ_MG <- ScaleData(YoungCPZ_MG, verbose = FALSE)
YoungCPZ_MG <- RunPCA(YoungCPZ_MG, features = VariableFeatures(object = YoungCPZ_MG), verbose = FALSE)
ElbowPlot(YoungCPZ_MG)
YoungCPZ_MG <- FindNeighbors(YoungCPZ_MG, dims = 1:15)
YoungCPZ_MG <- FindClusters(YoungCPZ_MG, resolution = 0.1)
YoungCPZ_MG <- RunUMAP(YoungCPZ_MG, dims = 1:15)
DimPlot(YoungCPZ_MG, split.by = "Condition", reduction = "umap", label = TRUE, ncol=2)

saveRDS(YoungCPZ_MG, "YoungCPZ_MG_15dims_res0.1.rds")