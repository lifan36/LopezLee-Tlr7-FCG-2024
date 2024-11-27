###Fig. S6 (FCG 3wk CPZ dataset)
library(Seurat)
library(enrichplot)
library(DOSE)
library(KEGGREST)
library(clusterProfiler)
library(org.Mm.eg.db)

#Fig. S6A
#RDS objects are the same as used for Fig 1
OL <- readRDS("FCG_OL_final_res0.2.rds") 
DefaultAssay(OL) <- 'RNA'
DimPlot(OL, reduction = "umap", split.by = "Genotype.Food", ncol = 4, raster = F, label = T)

#Fig S6B
DotPlot(OL, split.by = "seurat_clusters", features = c("Mag", "Mog", "Plp1", "Mbp", "Pcdh9", "Chn2",
                                                        "Tcf7l2", "Slc1a1", "Bcan", "Ssh3", "Gng4", "Prickle1", 
                                                        "Grik2", "Grik1", "Itpr2", "Sox6", "Vcan", "Cspg4", "Pdgfra", 
                                                       "Dscam", "Neu4"), 
        cols = c("coral1","darkgoldenrod", "chartreuse4", "cyan3","cadetblue","darkorchid"))

##Fig. S6C
#identify OL subcluster 4 markers
FCG_OL4_DEGs <- FindMarkers(OL, ident.1= "4",logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
FCG_OL4_DEGs <- subset(FCG_OL4_DEGs, FCG_OL4_DEGs$p_val_adj<0.05)
write.csv(FCG_OL4_DEGs, "FCG_OL4_DEGs.csv")
FCG_OL4_DEGs <- read.csv("FCG_OL4_DEGs.csv")
#run through KEGG
kegg_organism = "mmu"
FCG_OL4_DEGs <- bitr(FCG_OL4_DEGs$Gene, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)

kegg_enrich <-  enrichKEGG(gene  = FCG_OL4_DEGs$ENTREZID,
                           organism     = 'mmu',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'fdr')
kegg <- setReadable(kegg_enrich, 'org.Mm.eg.db', 'ENTREZID')
write.csv(kegg, "OL4_KEGG.csv")
kegg <- read.csv("OL4_KEGG.csv")

#plot
categorys <- c("Axon guidance - Mus musculus (house mouse)", "Cell adhesion molecules - Mus musculus (house mouse)",
               "Endocytosis - Mus musculus (house mouse)", "Regulation of actin cytoskeleton - Mus musculus (house mouse)", "Phospholipase D signaling pathway - Mus musculus (house mouse)")
cnetplot(kegg, circular=T, colorEdge=T, showCategory=categorys)

##Fig. 6D
MG <- readRDS("FCG_MG_reclusted_res0.2.rds") 
DefaultAssay(MG) <- 'RNA'
Idents(MG) <- "seurat_clusters"

#renumber
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
#plot
DimPlot(MG, reduction = "umap", split.by = "Genotype.Food", ncol = 4, raster = F, label=T)


