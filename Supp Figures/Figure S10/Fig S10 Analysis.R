##Figure S10 (spatial transcriptomics)
library(BayesSpace)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(loomR)
library(patchwork)
library(SeuratData)
library(UpSetR)
library(RColorBrewer)
library(gplots)
library(pals)
library(ggrepel)
library(ggpubr)
library(dplyr)

##Fig S10A
harmonized_seurat <- readRDS("XXO_XYT_clean_harm.rds")
DimPlot(harmonized_seurat, reduction = "umap", split.by = "Condition", label = TRUE, pt.size = .1, cols = "alphabet2")

##Fig S10B
#find cluster markers
markers <- FindAllMarkers(harmonized_seurat, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(markers, "markers.csv")

#Bayespace to increase image resolution
#pick representative image
Idents(harmonized_seurat) <- "orig.ident"
FCG_193 <- subset(harmonized_seurat, idents = c("FCG_193"))

#Convert to SCE
FCG_193_Bayes = as.SingleCellExperiment(FCG_193) #convert seurat to SCE
colData(FCG_193_Bayes) = cbind(colData(FCG_193_Bayes), FCG_193@images$FCG_193@coordinates) #add spatial info to SCE

#BayesSpace Workflow
FCG_193_Bayes = spatialPreprocess(FCG_193_Bayes, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
FCG_193_Bayes = spatialCluster(FCG_193_Bayes, platform = "Visium", nrep = 10000, burn.in = 100, q = 24) # cluster via BayesSpace
clusterPlot(FCG_193_Bayes) #plot via BayesSpace
FCG_193_Bayes.enhanced <- spatialEnhance(FCG_193_Bayes, q=24, platform="Visium", d=7,
                                         model="t", gamma=2,
                                         jitter_prior=0.3, jitter_scale=3.5,
                                         nrep=10000, burn.in=100,
                                         save.chain=TRUE)

clusterPlot(FCG_193_Bayes.enhanced)

#now need to enhance gene expression, in this case just doing all genes
# extract gene names
all_genes = rownames(FCG_193_Bayes.enhanced)
FCG_193_Bayes.enhanced.allfeatures <- enhanceFeatures(FCG_193_Bayes.enhanced, FCG_193_Bayes, nrounds=0, feature_names = all_genes) 

saveRDS(FCG_193_Bayes.enhanced.allfeatures, "FCG_193_Bayes.enhanced.allfeatures.rds")

#show rep section colored by expression of marker gene for subclusters of interest
#clust 7
featurePlot(FCG_193_Bayes.enhanced.allfeatures, "Plp1", color=NA)+scale_fill_viridis_c(option = "inferno")+coord_flip() + scale_x_reverse()

#clust 13
featurePlot(FCG_193_Bayes.enhanced.allfeatures, "Wipf3", color=NA)+scale_fill_viridis_c(option = "inferno")+coord_flip() + scale_x_reverse()

#clust 1
featurePlot(FCG_193_Bayes.enhanced.allfeatures, "Egr1", color=NA)+scale_fill_viridis_c(option = "inferno")+coord_flip() + scale_x_reverse()

#clust4
featurePlot(FCG_193_Bayes.enhanced.allfeatures, "Ccn3", color=NA)+scale_fill_viridis_c(option = "inferno")+coord_flip() + scale_x_reverse()

##Fig S10C
#compare XXO CPZ and XYT CPZ within cluster 7
DefaultAssay(harmonized_seurat) <- 'SCT'
Idents(harmonized_seurat) <- "seurat_clusters"

#identify cluster 7 markers
harmonized_seurat_7 <- FindMarkers(harmonized_seurat, ident.1 = "7", min.pct = 0.1, logfc.threshold = 0.1)

#subset cluster 7 to compare btwn sexes
harmonized_seurat_clust7 <- subset(harmonized_seurat, idents = "7")
#Identify DEGs
Idents(harmonized_seurat_clust7) <- "Condition"
clust7_XXOCPZvsXYTCPZ <- FindMarkers(harmonized_seurat_clust7, ident.1 = "XXO_Cup", ident.2 = "XYT_Cup", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(clust7_XXOCPZvsXYTCPZ, "clust7_XXOCPZvsXYTCPZ.csv")

#Volcano plot
clust7_XXOCPZvsXYTCPZ$color[clust7_XXOCPZvsXYTCPZ$avg_log2FC > 0.1 & clust7_XXOCPZvsXYTCPZ$p_val_adj < 0.05] <- "red"
clust7_XXOCPZvsXYTCPZ$color[clust7_XXOCPZvsXYTCPZ$avg_log2FC < -0.1 & clust7_XXOCPZvsXYTCPZ$p_val_adj < 0.05] <- "blue"
clust7_XXOCPZvsXYTCPZ <- subset(clust7_XXOCPZvsXYTCPZ, clust7_XXOCPZvsXYTCPZ$p_val_adj<0.05)
ggplot(clust7_XXOCPZvsXYTCPZ, aes(x = avg_log2FC,y = -log(p_val_adj), color=color, label = Gene))+ 
  geom_point()+
  theme_classic(base_size = 20)+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey","black"="black"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+  
  geom_text_repel(show.legend = "none")+
  theme(legend.position="none")+
  geom_vline(xintercept=c(0.1, -0.1),linetype=3)

##Fig S10D
#read in gene lists
DAM <- read.csv("DAM_KerenShaul_0.1.csv")
DAM$Gene <- DAM$gene
DAO <- read.csv("DAO_0.1.csv")
DAO$Gene <- DAO$symbol
DAA <- read.csv("DAA_Habib_0.1.csv")
RemyOL <- read.csv("RemyeOL_Voskuhl_0.1.csv")

#identify overlapping genes of each list with cluster 7
Clust7_DAM_overlap <- merge(clust7_XXOCPZvsXYTCPZ, DAM, by = "Gene")
Clust7_DAO_overlap <- merge(clust7_XXOCPZvsXYTCPZ, DAO, by = "Gene")
Clust7_DAA_overlap <- merge(clust7_XXOCPZvsXYTCPZ, DAA, by = "Gene")
Clust7_RemyeOL_overlap <- merge(clust7_XXOCPZvsXYTCPZ, RemyOL, by = "Gene")

#upset plot
listInput <- list(clust7_XXOCPZvsXYTCPZ = clust7_XXOCPZvsXYTCPZ$Gene, DAM = Clust7_DAM_overlap$Gene,
                  DAO = Clust7_DAO_overlap$Gene,RemyeOL = Clust7_RemyeOL_overlap$Gene, DAA = Clust7_DAA_overlap$Gene)
fromList(listInput)
upset(fromList(listInput), order.by = "freq")

#add in rep section with gene module score
# create enhanced seurat object for module score calculation
FCG_193_Bayes.enhanced.seurat <- as.Seurat(FCG_193_Bayes.enhanced.allfeatures, data = "logcounts", counts = "logcounts") 
saveRDS(FCG_193_Bayes.enhanced.seurat, "FCG_193_Bayes.enhanced.seurat.rds")

# calculate gene module scores
DefaultAssay(FCG_193_Bayes.enhanced.seurat) <- "originalexp"
NormalizeData(FCG_193_Bayes.enhanced.seurat)
#need to find overlap btwn seurat object and each list
#first need list of genes in seurat object
FCG_193_Bayes.enhanced.seurat_genes <- rownames(FCG_193_Bayes.enhanced.seurat@assays$originalexp@counts)
FCG_193_Bayes.enhanced.seurat_genes <- as.data.frame(FCG_193_Bayes.enhanced.seurat_genes)
FCG_193_Bayes.enhanced.seurat_genes$Gene <- FCG_193_Bayes.enhanced.seurat_genes$FCG_193_Bayes.enhanced.seurat_genes
FCG_193_Bayes.enhanced.seurat_genes <- FCG_193_Bayes.enhanced.seurat_genes[ -c(1) ]

FCG_193_DAM_overlap <- merge(FCG_193_Bayes.enhanced.seurat_genes, DAM, by = "Gene")

FCG_193_Bayes.enhanced.seurat <- AddModuleScore(FCG_193_Bayes.enhanced.seurat, features = list(FCG_193_DAM_overlap), name = "DAM_module", assay = "originalexp", nbin = 10)
#will show any genes in not sig detected in seurat object
#can call metadata to check module score is in there
meta <- FCG_193_Bayes.enhanced.seurat@meta.data

#spatial enrichment plot
featurePlot(FCG_193_Bayes.enhanced.allfeatures, feature = FCG_193_Bayes.enhanced.seurat$DAM_module1) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 1") + NoLegend()+coord_flip() + scale_x_reverse()

#plot legend, only need 1 since using same color scheme and goes min-max 
leg <- as_ggplot(get_legend(featurePlot(FCG_193_Bayes.enhanced.allfeatures, feature = FCG_193_Bayes.enhanced.seurat$DAM_module1) + 
                              scale_fill_viridis_c(option = "inferno", guide = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = FALSE, label = FALSE)) + 
                              theme(aspect.ratio = 1, legend.text=element_text(size=8), legend.position = "top") + labs(fill = "Relative Gene\nModule Enrichment")))
postscript("relative_gene_module_enrichment_legend_clust7.ps", width = 6, height = 4)

#repeat for other DAO
FCG_193_DAO_overlap <- merge(FCG_193_Bayes.enhanced.seurat_genes, DAO, by = "Gene")
FCG_193_Bayes.enhanced.seurat <- AddModuleScore(FCG_193_Bayes.enhanced.seurat, features = list(FCG_193_DAO_overlap), name = "DAO_module", assay = "originalexp", nbin = 10)
DAOplot <- featurePlot(FCG_193_Bayes.enhanced.allfeatures, feature = FCG_193_Bayes.enhanced.seurat$DAO_module1) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 1") + NoLegend()+coord_flip() + scale_x_reverse()
DAOplot

#repeat for RemyeOL
FCG_193_RemyeOL_overlap <- merge(FCG_193_Bayes.enhanced.seurat_genes, RemyOL, by = "Gene")
FCG_193_Bayes.enhanced.seurat <- AddModuleScore(FCG_193_Bayes.enhanced.seurat, features = list(FCG_193_RemyeOL_overlap), name = "RemyeOL_module", assay = "originalexp", nbin = 10)
RemyeOLplot <- featurePlot(FCG_193_Bayes.enhanced.allfeatures, feature = FCG_193_Bayes.enhanced.seurat$RemyeOL_module1) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 1") + NoLegend()+coord_flip() + scale_x_reverse()
RemyeOLplot

##Fig S10E
#need DEG lists between XXO CPZ and XYT CPZ within cluster 1
#subset cluster 1
Idents(harmonized_seurat) <- "seurat_clusters"
harmonized_seurat_clust1 <- subset(harmonized_seurat, idents = "1")
#identify DEGs (data in Table S6)
Idents(harmonized_seurat_clust1) <- "Condition"
clust1_XXOCPZvsXYTCPZ <- FindMarkers(harmonized_seurat_clust1, ident.1 = "XXO_Cup", ident.2 = "XYT_Cup", min.pct = 0.1, logfc.threshold = 0.1)
clust1_XXOCPZvsXYTCPZ_sig <- subset(clust1_XXOCPZvsXYTCPZ, clust1_XXOCPZvsXYTCPZ$p_val_adj<0.05)
write.csv(clust1_XXOCPZvsXYTCPZ, "clust1_XXOCPZvsXYTCPZ.csv")

#repeat for cluster 13
#subset
Idents(harmonized_seurat) <- "seurat_clusters"
harmonized_seurat_clust13 <- subset(harmonized_seurat, idents = "13")

#DEG list
Idents(harmonized_seurat_clust13) <- "Condition"
clust13_XXOCPZvsXYTCPZ <- FindMarkers(harmonized_seurat_clust13, ident.1 = "XXO_Cup", ident.2 = "XYT_Cup", min.pct = 0.1, logfc.threshold = 0.1)
clust13_XXOCPZvsXYTCPZ <- subset(clust13_XXOCPZvsXYTCPZ, clust13_XXOCPZvsXYTCPZ$p_val_adj<0.05)

#clust 1 and 7
Clust1_Clust7_SexCPZ_DEGs <- merge(clust1_XXOCPZvsXYTCPZ_sig, clust7_XXOCPZvsXYTCPZ, by = "Gene")
write.csv(Clust1_Clust7_SexCPZ_DEGs, "Overlap_Clust1and7_XXOCPZvsXYTCPZ.csv")

#now clust13 and clust 7 overlap
Clust13_Clust7_SexCPZ_DEGs <- merge(clust7_XXOCPZvsXYTCPZ, clust13_XXOCPZvsXYTCPZ, by = "Gene")
write.csv(Clust13_Clust7_SexCPZ_DEGs, "Overlap_Clust7and13_XXOCPZvsXYTCPZ.csv")

#clust 1 and 13
Clust1_Clust13_SexCPZ_DEGs <- merge(clust1_XXOCPZvsXYTCPZ_sig, clust13_XXOCPZvsXYTCPZ, by = "Gene")
write.csv(Clust1_Clust13_SexCPZ_DEGs, "Overlap_Clust1and13_XXOCPZvsXYTCPZ.csv")

#now overlap all 3
Clust13_Clust1_Clust7_SexCPZ_DEGs <- merge(Clust1_Clust7_SexCPZ_DEGs, clust13_XXOCPZvsXYTCPZ, by = "Gene")
write.csv(Clust13_Clust1_Clust7_SexCPZ_DEGs, "Overlap_AllClust_XXOCPZvsXYTCPZ.csv")

#clust7+13 vs all
All_Clust7vs13_overlap <- merge(Clust13_Clust1_Clust7_SexCPZ_DEGs, Clust13_Clust7_SexCPZ_DEGs, by = "Gene")
write.csv(All_Clust7vs13_overlap, "Overlap_All_Clust7and13_XXOCPZvsXYTCPZ.csv")

#clust1+13 vs all
All_Clust1vs13_overlap <- merge(Clust13_Clust1_Clust7_SexCPZ_DEGs, Clust1_Clust13_SexCPZ_DEGs, by = "Gene")
write.csv(All_Clust1vs13_overlap, "Overlap_All_Clust1and13_XXOCPZvsXYTCPZ.csv")

#clust7+1 vs all
All_Clust1vs7_overlap <- merge(Clust13_Clust1_Clust7_SexCPZ_DEGs, Clust1_Clust7_SexCPZ_DEGs, by = "Gene")
write.csv(All_Clust1vs7_overlap, "Overlap_All_Clust1and7_XXOCPZvsXYTCPZ.csv")

##Fig S10F
#heatmaps
#read in GO list of Myelination
myelination <- read.csv("GO_myelination.csv")
myelination$Gene <- myelination$Symbol

#overlap for each cluster
overlap.mye.clust7 <- clust7_XXOCPZvsXYTCPZ %>%
  inner_join(myelination, by = "Gene")
write.csv(overlap.mye.clust7, "overlap.mye.clust7.csv")

overlap.mye.clust1 <- clust1_XXOCPZvsXYTCPZ_sig %>%
  inner_join(myelination, by = "Gene")
write.csv(overlap.mye.clust1, "overlap.mye.clust1.csv")

overlap.mye.clust13 <- clust13_XXOCPZvsXYTCPZ %>%
  inner_join(myelination, by = "Gene")
write.csv(overlap.mye.clust13, "overlap.mye.clust13.csv")

#format and read in
mye_heatmap <- read.csv("Myelination_Cluster_Heatmap.csv", row.names = 1)
mye_heatmap <- as.matrix(mye_heatmap)

col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(20))
heatmap.2(mye_heatmap, trace = "none", scale = "none", dendrogram = "none", Rowv = T, density.info = "none", col =  col, breaks=seq(-1,1,0.1), Colv = F)

#repeat for synaptic transmission
synaptictrans <- read.csv("GO_SynTrans.csv")
synaptictrans$Gene <- synaptictrans$Symbol

overlap.syntrans.clust7 <- clust7_XXOCPZvsXYTCPZ %>%
  inner_join(synaptictrans, by = "Gene")
write.csv(overlap.syntrans.clust7, "overlap.syntrans.clust7.csv")

overlap.syntrans.clust13 <- clust13_XXOCPZvsXYTCPZ %>%
  inner_join(synaptictrans, by = "Gene")
overlap.syntrans.clust13 <- subset(overlap.syntrans.clust13, overlap.syntrans.clust13$p_val_adj<0.05)
write.csv(overlap.syntrans.clust13, "overlap.syntrans.clust13.csv")

overlap.syntrans.clust1 <- clust1_XXOCPZvsXYTCPZ %>%
  inner_join(synaptictrans, by = "Gene")
overlap.syntrans.clust1 <- subset(overlap.syntrans.clust1, overlap.syntrans.clust11$p_val_adj<0.05)
write.csv(overlap.syntrans.clust1, "overlap.syntrans.clust1.csv")

#read formatted back in
syntrans_heatmap <- read.csv("SynTrans_Cluster_Heatmap.csv", row.names = 1)
syntrans_heatmap <- as.matrix(syntrans_heatmap)
#heatmap
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(20))
heatmap.2(syntrans_heatmap, trace = "none", scale = "none", dendrogram = "none", Rowv = T, density.info = "none", col =  col, breaks=seq(-1,1,0.1), Colv = F, labCol = T)

#again for cytosolic ribo genes
cytoribo <- read.csv("GO_CytosolicRibo.csv")
cytoribo$Gene <- cytoribo$Symbol
overlap.ribo.clust7 <- clust7_XXOCPZvsXYTCPZ %>%
  inner_join(cytoribo, by = "Gene")
write.csv(overlap.ribo.clust7, "overlap.cytoribo.clust7.csv")

overlap.ribo.clust1 <- clust1_XXOCPZvsXYTCPZ_sig %>%
  inner_join(cytoribo, by = "Gene")
write.csv(overlap.ribo.clust1, "overlap.cytoribo.clust1.csv")

overlap.ribo.clust13 <- clust13_XXOCPZvsXYTCPZ %>%
  inner_join(cytoribo, by = "Gene")
write.csv(overlap.ribo.clust13, "overlap.cytoribo.clust13.csv")
#read in formatted
ribo_heatmap <- read.csv("CytoRibo_Cluster_Heatmap.csv", row.names = 1)
ribo_heatmap <- as.matrix(ribo_heatmap)
#heatmap
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(20))
heatmap.2(ribo_heatmap, trace = "none", scale = "none", dendrogram = "none", Rowv = T, density.info = "none", col =  col, breaks=seq(-1,1,0.1), Colv = F)
