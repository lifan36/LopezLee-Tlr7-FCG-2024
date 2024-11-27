###Figure 1 snRNAseq Analysis
remotes::install_version("Seurat", version = "4.0.0")
library(Seurat)
library(dplyr)
library(ggpubr)
library(ggrepel)

##Fig 1H, cell ratios plotted in Prism
##Fig 1I
OL <- readRDS("FCG_OL_final_res0.2.rds")
DefaultAssay(OL) <- "RNA"
Idents(OL) <- "seurat_clusters"

#rename clusters to start with 1
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)

#DEG list
OL4 <- FindMarkers(OL, ident.1= "4",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
OL4 <- subset(OL4, OL4$p_val_adj<0.05)
OL4$gene <- row.names(OL4)

#identify overlapping genes with DAO
DAO <- read.csv("DAO_0.1.csv")
DAO$gene <- DAO$symbol
overlap.OL4_DAO <- OL4 %>%
  inner_join(DAO, by = "gene")

#export and format, one col for gene name, one col for OL4 log2FC, one col for DAM log2FC
write.csv(overlap.OL4_DAO, "FCG_OL4_DEGs_overlap2.csv")

#read formatted back in
overlap.OL4_DAO <- read.csv("FCG_OL4_DEGs_overlap.csv", row.names = 1)
overlap.OL4_DAO.df<- as.data.frame(overlap.OL4_DAO)
overlap.OL4_DAO.df$gene <- row.names(overlap.OL4_DAO.df)

#plot
ggscatter(overlap.OL4_DAO.df, x = "OL4", y = "DAO", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "OL4", ylab = "DAO",
          repel=TRUE)+
  geom_text_repel(label = overlap.OL4_DAO.df$gene, max.overlaps = 10)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")

##Fig 1J, cell ratios plotted in Prism
##Figure 1K
MG <- readRDS("FCG_MG_reclusted_res0.2.rds")
DefaultAssay(MG) <- "RNA"
Idents(MG) <- "seurat_clusters"

#rename clusters to start with 1
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)

#DEG list ################FAN please check###################
MG3 <- FindMarkers(MG, ident.1 = "3", min.pct = 0.1, logfc.threshold = 0.1, test.use = "MAST")
MG3 <- subset(MG3, MG3$p_val_adj<0.05)
MG3$gene <- row.names(MG3)
#####
#DEGs were from Fan's list, thresholded avg_logFC by 0.1 in Excel and read in
MG3 <- read.csv("FCG_MG3_vs_others_DEGs_0.1.csv")

#get genes shared between the 2 datasets
DAM <- read.csv("DAM_KerenShaul_0.1.csv")
overlap.MG3_DAM <- MG3 %>%
  inner_join(DAM, by = "gene")

#export and write out to csv to format, one col for gene name, one col for MG3 log2FC, one col for DAM log2FC
write.csv(overlap.MG3_DAM, "overlap.MG3_DAM>0.1.csv")

#read formatted back in
overlap.MG3_DAM <- read.csv("DAM_MG3_overlapdf2.csv")
overlap.MG3_DAM.df<- as.data.frame(overlap.MG3_DAM)

#graph as correlation
ggscatter(overlap.MG3_DAM.df, x = "MG3", y = "DAM", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MG3", ylab = "DAM_KerenShaul",
          repel=TRUE)+
  geom_text_repel(label = overlap.MG3_DAM.df$gene, max.overlaps = 10)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")

##Fig. 1L Dotplot
DotPlot(MG, idents = "3", split.by = "Genotype.Food", features = "Lpl", 
        cols = "RdBu", scale.by = "size")

##Figure 1O, MG4 DEGs + volcano
MG4 <- FindMarkers(MG, ident.1 = "4", min.pct = 0.1, logfc.threshold = 0.1, test.use = "MAST")
MG4 <- subset(MG4, MG4$p_val_adj<0.05)
############################
#MG4 DEGs made by Fan, read in
MG4 <- read.csv("FCG_MG4_vs_others_DEGs.csv")
MG4$color[MG4$avg_logFC > 0.1 & MG4$p_val_adj < 0.05] <- "red"
MG4$color[MG4$avg_logFC < -0.1 & MG4$p_val_adj < 0.05] <- "blue"

ggplot(MG4, aes(x = avg_logFC,y = -log(p_val_adj), color=color, label = gene))+ 
  geom_point()+
  theme_classic(base_size = 20)+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey","black"="black"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+  
  geom_text_repel(data = MG4, 
                  aes((label = MG4$gene), color = "black", size=5) + geom_point(size=5), max.overlaps = 50)+
  theme(legend.position="none")+
  geom_vline(xintercept=c(0.1, -0.1),linetype=3)

