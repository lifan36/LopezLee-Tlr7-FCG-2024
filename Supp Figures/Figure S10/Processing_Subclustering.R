##processing
Spat <- readRDS("FCG_Cup_5week_Visium_June2023.rds")
Idents(Spat) <- "Condition"
DefaultAssay(Spat) <- 'Spatial'

Spat_XXO_XYT_clean <- subset(Spat, idents = c("XXO_Ctrl", "XXO_Cup", "XYT_Ctrl", "XYT_Cup"))
Spat_XXO_XYT_clean <- Spat_XXO_XYT_clean %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform()
Spat_XXO_XYT_clean <- RunPCA(Spat_XXO_XYT_clean, assay = "SCT", npcs = 50)
NormalizeData(Spat_XXO_XYT_clean, verbose = FALSE)
FindVariableFeatures(Spat_XXO_XYT_clean, selection.method = "vst", nfeatures = 2000)
ScaleData(Spat_XXO_XYT_clean, verbose = FALSE)
Spat_XXO_XYT_clean <- SCTransform(Spat_XXO_XYT_clean, assay = "SCT", verbose = FALSE)

RunPCA(Spat_XXO_XYT_clean, pc.genes = Spat_XXO_XYT_clean@var.genes, npcs = 20, verbose = FALSE)
#harmonize/subcluster
pc_mat <- Spat_XXO_XYT_clean@reductions$pca@cell.embeddings
meta_data <- harmonized_seurat@meta.data
harmony_embeddings <- HarmonyMatrix(pc_mat, meta_data = meta_data, 
                                    vars_use = "orig.ident", plot_convergence = TRUE)
Spat_XXO_XYT_clean <- CreateDimReducObject(
  embeddings = Spat_XXO_XYT_clean$rotation,
  loadings = pcs$x,
  stdev = pcs$sdev,
  key = "PC",
  assay = "harmony"
)
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:6)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
DimPlot(object = harmonized_seurat, reduction = "harmony", pt.size = .1, group.by = "orig.ident")

#rename cluster so starts with 1
n <- dim(table(harmonized_seurat@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
harmonized_seurat@active.ident <- plyr::mapvalues(x = harmonized_seurat@active.ident, from = current.cluster.ids, to = new.cluster.ids)
harmonized_seurat@active.ident <- factor(harmonized_seurat@active.ident, levels=1:n)
saveRDS(harmonized_seurat, file = 'XXO_XYT_clean_harm.rds')

