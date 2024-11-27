###Fig S15
library(ggpubr)
library(ggrepel)

##Fig S15A, same dataset as Fig 6
LG70_integr_MG <- readRDS("Tlr7_MG_reclusted_res0.3.rds")
DefaultAssay(LG70_integr_MG) <- "RNA"

VlnPlot(
  object = LG70_integr_MG,
  idents = c("1","2","3","4","5","6"),
  log = FALSE,
  features = c("Lpl", "Myo1e", "Axl", "Csf1"),
  pt.size = 0, ncol=2)

##Fig S15B
OL <- readRDS("Tlr7_OL_reclusted_res0.15.rds")
DefaultAssay(OL) <- "RNA"

DimPlot(OL, reduction = "umap", split.by = "Condition", label = TRUE, ncol=4)

##Fig S15C (in Table S9)
OL_clust4_markers <- FindMarkers(OL, ident.1 = "4", min.pct = 0.25, logfc.threshold = 0.1, test.use = "MAST")
OL_clust4_markers <- subset(OL_clust4_markers, OL_clust4_markers$p_val_adj<0.05)
OL_clust4_markers$gene <- row.names(OL_clust4_markers)
write.csv(OL_clust4_markers, "OL_clust4_markers.csv")

OL_clust4_markers$color[OL_clust4_markers$avg_log2FC > 0.1 & OL_clust4_markers$p_val_adj < 0.05] <- "red"
OL_clust4_markers$color[OL_clust4_markers$avg_log2FC < -0.1 & OL_clust4_markers$p_val_adj < 0.05] <- "blue"

OL_clust4_markers$Gene <- row.names(OL_clust4_markers)

ggplot(OL_clust4_markers, aes(x = avg_log2FC,y = -log(p_val_adj), color=color, label = Gene))+ 
  geom_point()+
  theme_classic(base_size = 20)+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey","black"="black"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+  
  geom_text_repel(data = OL_clust4_markers, 
                  aes((label = OL_clust4_markers$Gene), color = "black", size=5) + geom_point(size=5), max.overlaps = 10)+
  theme(legend.position="none")+
  geom_vline(xintercept=c(0.1, -0.1),linetype=3)

##Fig S15D
DAO <- read.csv("DAO_0.1.csv")
DAO$gene <- DAO$symbol
#overlap
OL4_DAO_overlap <- merge(OL_clust4_markers, DAO, by = "gene")

#write out and format
write.csv(OL4_DAO_overlap, "OL4_DAO_overlap.csv")
OL4_DAO_Corr <- read.csv("OL4_DAO_Corr.csv")
#plot
ggscatter(OL4_DAO_Corr, x = "OL4", y = "DAO", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "OL4", ylab = "DAO Pandey et al.",
          repel=TRUE)+
  geom_text_repel(label = OL4_DAO_Corr$gene, max.overlaps = 10)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")

