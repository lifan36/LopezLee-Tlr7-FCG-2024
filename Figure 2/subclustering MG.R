##subcluster FCG 5wk CPZ MG dataset
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:5)
MG <- FindClusters(MG, resolution = 0.2)
MG <- RunUMAP(MG, dims = 1: 5)

DimPlot(MG, reduction = "umap", split.by = "Condition", ncol = 4)

n <- dim(table(Demye_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
Demye_MG@active.ident <- plyr::mapvalues(x = Demye_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
Demye_MG@active.ident <- factor(Demye_MG@active.ident, levels=1:n)
saveRDS(MG, file = 'Chloe_FCG_MG_reclusted_res0.2.rds')