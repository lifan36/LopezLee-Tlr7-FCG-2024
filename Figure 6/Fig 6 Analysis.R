###Fig 6 (FCG integrated with TLR7KO snRNAseq & bulk RNA seq of PSUK mice with TLR7 inhibitor treatment)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(data.table)
install.packages("biomaRt")
library(biomaRt)
library(reshape2)
library(circlize)
library(stringr)
library(Seurat)
library(GeneOverlap)
library(ggpubr)
library(ggrepel)

##Fig 6B
LG70_integr_MG <- readRDS("Tlr7_MG_reclusted_res0.3.rds")
DefaultAssay(LG70_integr_MG) <- "RNA"
Idents(LG70_integr_MG) <- "Condition"
VlnPlot(LG70_integr_MG, features = "Tlr7", pt.size = 0)
#changed colors in Illustrator

##Fig 6C
#calculate cell ratios for all groups
a<-as.data.frame(table(LG70_integr_MG$seurat_clusters,LG70_integr_MG$Sample_Name))
colnames(a)<-c("clusters","genotype","cell.no")
agg<-aggregate(cell.no~genotype,a,sum)
a$genotype.total <- agg$cell.no[match(a$genotype,agg$genotype)]
a$ratio<-a$cell.no/a$genotype.total
#manually make .csv listing Condition for each subcluster+sample (e.g. 3 samples x 5 subclusters = 15 rows of that Condition) 
groups <- read.csv("MG_CellRatioGroups.csv", header = 1) #read in groups
a$groups <- groups$Group
write.csv(a, "MG_CellRatio.csv")
#calculate avg and sem calculated in Excel
#split into 2 spreadsheets, one for FCG and one for TLR7KO groups
#read in FCG to graph first
a <- read.csv("MG_FCG_CellRatio.csv", header=1)
#if need to force order bars
a$gclass <- factor(a$groups, levels=c('XXO_Ctrl', 'XXO_CPZ', 'XYT_Ctrl', 'XYT_CPZ'))
#plot
ggplot(a, aes(clusters, ratio, fill=gclass))+
  geom_errorbar( aes(ymin=Avg-sem, ymax=Avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#cluster numbers on x-axis changed to start at "1" in Illustrator
#last 2 clusters were not graphed due to low ratio (<2% in all samples)

##Fig 6D
#read in TLR7KO Excel with cell ratios
a <- read.csv("MG_TLR7KO_CellRatio.csv", header=1)
#if need to force order bars
a$gclass <- factor(a$groups, levels=c('KO_F_Ctrl', 'KO_F_CPZ', 'KO_M_Ctrl', 'KO_M_CPZ'))
#plot
ggplot(a, aes(clusters, ratio, fill=gclass))+
  geom_errorbar( aes(ymin=Avg-sem, ymax=Avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="PuOr")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#last 2 clusters were not graphed due to low ratio (<2% in all samples)

##Fig 6E
OL <- readRDS("Tlr7_OL_reclusted_res0.15.rds")
DefaultAssay(OL) <- "RNA"

#cell ratio
a<-as.data.frame(table(OL$seurat_clusters,OL$Sample_Name))
colnames(a)<-c("clusters","sample","cell.no")
agg<-aggregate(cell.no~sample,a,sum)
a$genotype.total <- agg$cell.no[match(a$sample,agg$sample)]
a$ratio<-a$cell.no/a$genotype.total
groups <- read.csv("OL_CellRatioGroups.csv", header = 1) #read in groups
a$groups <- groups$Group
write.csv(a, "OL_CellRatio.csv")
#read back in with avg and sem, FCG only
a <- read.csv("OL_CellRatio_FCG.csv", header=1)
#if need to force order bars
a$gclass <- factor(a$groups, levels=c('XXO_Ctrl', 'XXO_CPZ', 'XYT_Ctrl', 'XYT_CPZ'))
#plot
ggplot(a, aes(clusters, ratio, fill=gclass))+
  geom_errorbar( aes(ymin=Avg-sem, ymax=Avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#change cluster #s to start at 0 in Illustrator

##Fig 6F
#read in Excel with avg and sem, TLR7KO only
a <- read.csv("OL_CellRatio_LG70.csv", header=1)
#if need to force order bars
a$gclass <- factor(a$groups, levels=c('KO_F_Ctrl', 'KO_F_CPZ', "KO_M_Ctrl", "KO_M_CPZ"))
##
ggplot(a, aes(clusters, ratio, fill=gclass))+
  geom_errorbar( aes(ymin=Avg-sem, ymax=Avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="PuOr")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))
#change cluster #s to start at 0 in Illustrator

##Fig 6M
data <- read.csv("gene_count.csv", header=T, row.names=1)
#remove duplicates
data <- subset(data, !duplicated(data$gene_name))
row.names(data) <- data$gene_name
#remove metadata cols at end
data <- data[1:(length(data)-9)]
filtered = data[rowSums(data)>15,] #removes any genes with <15 counts
write.csv(data, "counts_genes_sorted_PSUK.csv")

#metadata
sample <- colnames(data)
sex <- c(rep("F_",3),rep("M_",3), rep("F_",5), rep("M_",6))
treatment <- c(rep("Ctrl",6),rep("Inhib",11))
sex_treatment <- paste(sex,treatment,sep="_")
meta <- data.frame(sample=sample, sex=sex, treatment=treatment, sex_treatment=sex_treatment)
all(colnames(filtered) %in% meta$sample)

#make DEG lists
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment)
dds <- DESeq(dds)

# 1. M Inhib vs. M Ctrl
contrast_MInhibvsMCtrl <- c("sex_treatment","M__Inhib","M__Ctrl")
res_MInhibvsMCtrl_unshrunken <- results(dds,contrast=contrast_MInhibvsMCtrl,alpha=0.05)
res_MInhibvsMCtrl <- lfcShrink(dds,contrast=contrast_MInhibvsMCtrl,res=res_MInhibvsMCtrl_unshrunken, type="normal")
write.csv(res_MInhibvsMCtrl, "DE_MInhibvsMCtrl.csv")
MInhibvsMCtrl <- read.csv("DE_MInhibvsMCtrl.csv")
MInhibvsMCtrl <- subset(MInhibvsMCtrl, MInhibvsMCtrl$pvalue<0.05)

# 2. F Inhib vs. F Ctrl
contrast_FInhibvsFCtrl <- c("sex_treatment","F__Inhib","F__Ctrl")
res_FInhibvsFCtrl_unshrunken <- results(dds,contrast=contrast_FInhibvsFCtrl,alpha=0.05)
res_FInhibvsFCtrl <- lfcShrink(dds,contrast=contrast_FInhibvsFCtrl,res=res_FInhibvsFCtrl_unshrunken, type="normal")
write.csv(res_FInhibvsFCtrl, "DE_FInhibvsFCtrl.csv")
FInhibvsFCtrl <- read.csv("DE_FInhibvsFCtrl.csv")
FInhibvsFCtrl <- subset(FInhibvsFCtrl, FInhibvsFCtrl$pvalue<0.05)

#now do overlap 
go.obj.InhibvsCtrl <- newGeneOverlap(FInhibvsFCtrl$X, MInhibvsMCtrl$X, genome.size = 21988)
go.obj.InhibvsCtrl <- testGeneOverlap(go.obj.InhibvsCtrl)
print(go.obj.InhibvsCtrl)
#gives you # of genes in overlap (A+B), unique to F (A), and unique to M (B)

##Fig 6N
#get unique DEGs
#setdiff(x, y) tells you elements in x not in y
Unique_MInhibvsMCtrl <- setdiff(MInhibvsMCtrl$X, FInhibvsFCtrl$X)
Unique_MInhibvsMCtrl <- as.data.frame(Unique_MInhibvsMCtrl)
#need to map back with logfc values, so merge back with OG dataset
Unique_MInhibvsMCtrl$X <- Unique_MInhibvsMCtrl$Unique_MInhibvsMCtrl
Unique_MInhibvsMCtrl <- merge(MInhibvsMCtrl, Unique_MInhibvsMCtrl, by = "X")
write.csv(Unique_MInhibvsMCtrl, "Unique_MInhibvsMCtrl.csv")

GO_Unique_MInhibvsMCtrl_dn <- read.csv("GO_Unique_InhibvsCtrl_M.csv")
colnames(GO_Unique_MInhibvsMCtrl_dn) <- c("name", "genes_in_gsea", "genes_in_data", "k/K", "p-value", "FDR")
GO_Unique_MInhibvsMCtrl_dn$FDR <- as.numeric(GO_Unique_MInhibvsMCtrl_dn$FDR)
GO_Unique_MInhibvsMCtrl_dn$logFDR <- -log10(GO_Unique_MInhibvsMCtrl_dn$FDR)
GO_Unique_MInhibvsMCtrl_dn$enrich <- paste(GO_Unique_MInhibvsMCtrl_dn$genes_in_data, "/", GO_Unique_MInhibvsMCtrl_dn$genes_in_gsea, sep=" ")
GO_Unique_MInhibvsMCtrl_dn <- GO_Unique_MInhibvsMCtrl_dn[order(-GO_Unique_MInhibvsMCtrl_dn$logFDR),]
GO_Unique_MInhibvsMCtrl_dn$Name <-  gsub("GO_", "", GO_Unique_MInhibvsMCtrl_dn$name)
GO_Unique_MInhibvsMCtrl_dn$Name <- gsub("*_", " ", GO_Unique_MInhibvsMCtrl_dn$Name)
GO_Unique_MInhibvsMCtrl_dn$Name <-  factor(GO_Unique_MInhibvsMCtrl_dn$Name, levels=rev(GO_Unique_MInhibvsMCtrl_dn$Name))
area.color <- c("blue", "blue", "blue", "blue","blue", "blue")
#plot
ggplot(data=GO_Unique_MInhibvsMCtrl_dn, aes(x=reorder(name, logFDR), y=logFDR)) +
  theme_classic() +
  ylab("-log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity", fill=area.color) +
  geom_text(aes(label = enrich), vjust = 0.5, hjust = 1.25, colour = "white")+
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  theme(plot.title = element_text(hjust = -0.5))+
  ggtitle("GO_Unique_KOMyevsWTMye_M_dn")

##Fig 6O
#want overlapping genes
FInhibvsFCtrl_MInhibvsMCtrl_overlap <- merge(FInhibvsFCtrl, MInhibvsMCtrl, by = "X")
write.csv(FInhibvsFCtrl_MInhibvsMCtrl_overlap, "FInhibvsFCtrl_MInhibvsMCtrl_overlap.csv")

#correlate
FInhibvsFCtrl_MInhibvsMCtrl_Corr <- read.csv("FInhibvsFCtrl_MInhibvsMCtrl_corr.csv")
ggscatter(FInhibvsFCtrl_MInhibvsMCtrl_Corr, x = "FEMALElog2FoldChange.x", y = "MALElog2FoldChange.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "FEMALElog2FoldChange.x", ylab = "MALElog2FoldChange.x", 
          repel=TRUE)+
  geom_text_repel(label = FInhibvsFCtrl_MInhibvsMCtrl_Corr$Gene, max.overlaps = 30)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")
