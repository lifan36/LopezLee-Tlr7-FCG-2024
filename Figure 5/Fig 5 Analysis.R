##Fig 5 (bulk RNA sequencing of TLR7 Agonist in vivo & TLR7KO Primary MG datasets)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(data.table)
library(biomaRt)
library(reshape2)
library(circlize)
library(stringr)
library(GeneOverlap)

#Fig 5B
data <- read.csv("gene_count.csv", header=T)
data <- subset(data, !duplicated(data$gene_name))
row.names(data) <- data$gene_name
data <- data[,-1]
#remove 9 cols at the end that are non-numeric
data <- data[1:(length(data)-9)]
#save as new spreadsheet
write.csv(data, "counts_genes_sorted.csv")

filtered = data[rowSums(data)>15,] #removes any genes with <15 counts
#add in metadata
sample <- colnames(data)
sex <- c(rep("M_",14),rep("F_",14))
treatment <- c(rep("Saline",7),rep("GDQ",7), rep("Saline", 7), rep("GDQ",7))
sex_treatment <- paste(sex,treatment,sep="_")
meta <- data.frame(sample=sample, sex=sex, treatment=treatment, sex_treatment=sex_treatment)
all(colnames(filtered) %in% meta$sample)
#make DEG lists
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment)
dds <- DESeq(dds)

# 1. F GDQ vs. F Saline
contrast_FGDQvsFSaline <- c("sex_treatment","F__GDQ","F__Saline")
res_FGDQvsFSaline_unshrunken <- results(dds,contrast=contrast_FGDQvsFSaline,alpha=0.05)
res_FGDQvsFSaline <- lfcShrink(dds,contrast=contrast_FGDQvsFSaline,res=res_FGDQvsFSaline_unshrunken, type="normal")
write.csv(res_FGDQvsFSaline, "DE_FGDQvsFSaline.csv")
FGDQvsFSaline <- read.csv("DE_FGDQvsFSaline.csv")
FGDQvsFSaline <- subset(FGDQvsFSaline,FGDQvsFSaline$pvalue<0.05)

#2. M GDQ vs. M Saline
contrast_MGDQvsMSaline <- c("sex_treatment","M__GDQ","M__Saline")
res_MGDQvsMSaline_unshrunken <- results(dds,contrast=contrast_MGDQvsMSaline,alpha=0.05)
res_MGDQvsMSaline <- lfcShrink(dds,contrast=contrast_MGDQvsMSaline,res=res_MGDQvsMSaline_unshrunken, type="normal")
write.csv(res_MGDQvsMSaline, "DE_MGDQvsMSaline.csv")
MGDQvsMSaline <- read.csv("DE_MGDQvsMSaline.csv")
MGDQvsMSaline <- subset(MGDQvsMSaline,MGDQvsMSaline$pvalue<0.05)

#want to see # of DEGs that are overlapping/unique
go.obj.GDQvsSaline <- newGeneOverlap(FGDQvsFSaline$X, MGDQvsMSaline$X, genome.size = 21988)
go.obj.GDQvsSaline <- testGeneOverlap(go.obj.GDQvsSaline)
print(go.obj.GDQvsSaline)
#tells you number of genes in both (A+B), unique to FGDQvsFSaline (A), and unique to MGDQvsMSaline (B)

#get the unique DEGs for each list
#setdiff(x, y) tells you elements in x not in y
Unique_F_GDQvsSaline <- setdiff(FGDQvsFSaline$X, MGDQvsMSaline$X)
Unique_F_GDQvsSaline <- as.data.frame(Unique_F_GDQvsSaline)

#need to map back with logfc values, so merge back with OG dataset
Unique_F_GDQvsSaline$X <- Unique_F_GDQvsSaline$Unique_F_GDQvsSaline
Unique_F_GDQvsSaline <- merge(FGDQvsFSaline, Unique_F_GDQvsSaline, by = "X")
write.csv(Unique_F_GDQvsSaline, "Unique_F_GDQvsSaline.csv")

#now for males
#setdiff(x, y) tells you elements in x not in y
Unique_M_GDQvsSaline <- setdiff(MGDQvsMSaline$X, FGDQvsFSaline$X)
Unique_M_GDQvsSaline <- as.data.frame(Unique_M_GDQvsSaline)
Unique_M_GDQvsSaline$X <- Unique_M_GDQvsSaline$Unique_M_GDQvsSaline

Unique_M_GDQvsSaline <- merge(MGDQvsMSaline, Unique_M_GDQvsSaline, by = "X")
write.csv(Unique_M_GDQvsSaline, "Unique_M_GDQvsSaline.csv")

##Fig 5C
#put unique genes for M that were upregulated into GSEA
#graph GSEA results
GO_M_GDQvsSaline_Unique_up <- read.csv("GO_M_GDQvsSaline_Unique_up_Annotated.csv",header=T)

colnames(GO_M_GDQvsSaline_Unique_up) <- c("name", "genes_in_gsea", "genes_in_data", "k/K", "p-value", "FDR")
GO_M_GDQvsSaline_Unique_up$FDR <- as.numeric(GO_M_GDQvsSaline_Unique_up$FDR)
GO_M_GDQvsSaline_Unique_up$logFDR <- -log10(GO_M_GDQvsSaline_Unique_up$FDR)
GO_M_GDQvsSaline_Unique_up$enrich <- paste(GO_M_GDQvsSaline_Unique_up$genes_in_data, "/", GO_M_GDQvsSaline_Unique_up$genes_in_gsea, sep=" ")
GO_M_GDQvsSaline_Unique_up <- GO_M_GDQvsSaline_Unique_up[order(-GO_M_GDQvsSaline_Unique_up$logFDR),]
GO_M_GDQvsSaline_Unique_up$Name <-  gsub("GO_", "", GO_M_GDQvsSaline_Unique_up$name)
GO_M_GDQvsSaline_Unique_up$Name <- gsub("*_", " ", GO_M_GDQvsSaline_Unique_up$Name)
GO_M_GDQvsSaline_Unique_up$Name <-  factor(GO_M_GDQvsSaline_Unique_up$Name, levels=rev(GO_M_GDQvsSaline_Unique_up$Name))
area.color <- c("red1", "red1", "red1", "red1","red1", "red1", "red1", "red1")
#plot
ggplot(data=GO_M_GDQvsSaline_Unique_up, aes(x=reorder(name, logFDR), y=logFDR)) +
  theme_classic() +
  ylab("log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity", fill=area.color) +
  geom_text(aes(label = enrich), vjust = 0.5, hjust = 1.25, colour = "white")+
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  theme(plot.title = element_text(hjust = -0.5))+
  ggtitle("GO_M_GDQvsSaline_Unique_up")+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "lightgrey")

##Fig 5E, F (TLR7KO dataset)
data <- read.csv("gene_count.csv", header=T)
#remove duplicates
data <- subset(data, !duplicated(data$gene_name))
row.names(data) <- data$gene_name
data <- data[,-1]
#remove 9 cols at the end that are non-numeric
data <- data[1:(length(data)-9)]
#save as new spreadsheet
write.csv(data, "counts_genes_sorted.csv")

filtered = data[rowSums(data)>15,] #removes any genes with <15 counts
#metadata
sample <- colnames(data)
sex <- c(rep("M_",12),rep("F_",12))
treatment <- c(rep("DMSO",3),rep("Mye",3), rep("DMSO", 3), rep("Mye",3), rep("DMSO", 3), rep("Mye", 3), rep("DMSO", 3), rep("Mye", 3))
genotype <- c(rep("WT",6),rep("KO",6), rep("WT", 6), rep("KO",6))
sex_treatment_genotype <- paste(sex,treatment,genotype,sep="_")
meta <- data.frame(sample=sample, sex=sex, genotype = genotype, treatment=treatment, sex_treatment_genotype=sex_treatment_genotype)

all(colnames(filtered) %in% meta$sample)

#DEG lists
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment_genotype)
dds <- DESeq(dds)

# 1. KO M Mye vs. WT M Mye
contrast_KOMMyevsWTMMye <- c("sex_treatment_genotype","M__Mye_KO","M__Mye_WT")
res_KOMMyevsWTMMye_unshrunken <- results(dds,contrast=contrast_KOMMyevsWTMMye,alpha=0.05)
res_KOMMyevsWTMMye <- lfcShrink(dds,contrast=contrast_KOMMyevsWTMMye,res=res_KOMMyevsWTMMye_unshrunken, type="normal")
write.csv(res_KOMMyevsWTMMye, "DE_KOMMyevsWTMMye.csv")
KOMMyevsWTMMye <- read.csv("DE_KOMMyevsWTMMye.csv")
KOMMyevsWTMMye <- subset(KOMMyevsWTMMye, KOMMyevsWTMMye$pvalue<0.05)

# 2. KO F Mye vs. WT F Mye
contrast_KOFMyevsWTFMye <- c("sex_treatment_genotype","F__Mye_KO","F__Mye_WT")
res_KOFMyevsWTFMye_unshrunken <- results(dds,contrast=contrast_KOFMyevsWTFMye,alpha=0.05)
res_KOFMyevsWTFMye <- lfcShrink(dds,contrast=contrast_KOFMyevsWTFMye,res=res_KOFMyevsWTFMye_unshrunken, type="normal")
write.csv(res_KOFMyevsWTFMye, "DE_KOFMyevsWTFMye.csv")
KOFMyevsWTFMye <- read.csv("DE_KOFMyevsWTFMye.csv")
KOFMyevsWTFMye <- subset(KOFMyevsWTFMye, KOFMyevsWTFMye$pvalue<0.05)

#calculate DEG overlap
go.obj.KOvsWT_Mye <- newGeneOverlap(KOFMyevsWTFMye$X, KOMMyevsWTMMye$X, genome.size = 21988)
go.obj.KOvsWT_Mye <- testGeneOverlap(go.obj.KOvsWT_Mye)
print(go.obj.KOvsWT_Mye)

#get unique genes for each comparison
#setdiff(x, y) tells you elements in x not in y
Unique_KOMyevsWTMye_M <- setdiff(KOMMyevsWTMMye$X, KOFMyevsWTFMye$X)
Unique_KOMyevsWTMye_M <- as.data.frame(Unique_KOMyevsWTMye_M)
#need to map back with logfc values, so merge back with OG dataset
Unique_KOMyevsWTMye_M$X <- Unique_KOMyevsWTMye_M$Unique_KOMyevsWTMye_M
Unique_KOMyevsWTMye_M <- merge(KOMMyevsWTMMye, Unique_KOMyevsWTMye_M, by = "X")
write.csv(Unique_KOMyevsWTMye_M, "Unique_KOMyevsWTMye_M.csv")

#now for F
#setdiff(x, y) tells you elements in x not in y
Unique_KOMyevsWTMye_F <- setdiff(KOFMyevsWTFMye$X, KOMMyevsWTMMye$X)
Unique_KOMyevsWTMye_F <- as.data.frame(Unique_KOMyevsWTMye_F)
#need to map back with logfc values, so merge back with OG dataset
Unique_KOMyevsWTMye_F$X <- Unique_KOMyevsWTMye_F$Unique_KOMyevsWTMye_F
Unique_KOMyevsWTMye_F <- merge(KOFMyevsWTFMye, Unique_KOMyevsWTMye_F, by = "X")
write.csv(Unique_KOMyevsWTMye_F, "Unique_KOMyevsWTMye_F.csv")

#volcano the unique for males
Unique_KOMyevsWTMye_M$color[Unique_KOMyevsWTMye_M$log2FoldChange > 0.1 & Unique_KOMyevsWTMye_M$pvalue < 0.05] <- "red"
Unique_KOMyevsWTMye_M$color[Unique_KOMyevsWTMye_M$log2FoldChange < -0.1 & Unique_KOMyevsWTMye_M$pvalue < 0.05] <- "blue"

#label specific genes
Unique_KOMyevsWTMye_M$genelabels <- ""
Unique_KOMyevsWTMye_M$genelabels <- ifelse(Unique_KOMyevsWTMye_M$X == "Ifit1" 
                                           |Unique_KOMyevsWTMye_M$X == "Ifit2"
                                           |Unique_KOMyevsWTMye_M$X == "Ifit3"
                                           |Unique_KOMyevsWTMye_M$X == "Ifit3b"
                                           |Unique_KOMyevsWTMye_M$X == "Irf7"
                                           |Unique_KOMyevsWTMye_M$X == "Irf9"
                                           |Unique_KOMyevsWTMye_M$X == "Irf2"
                                           |Unique_KOMyevsWTMye_M$X == "Irf5"
                                           |Unique_KOMyevsWTMye_M$X == "Irf3"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi205"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi3b"
                                           |Unique_KOMyevsWTMye_M$X == "Ifit3"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi27l2a"
                                           |Unique_KOMyevsWTMye_M$X == "Ifit1"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi206"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi211"
                                           |Unique_KOMyevsWTMye_M$X == "Ifit2"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi44"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi204"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi203-ps"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi35"
                                           |Unique_KOMyevsWTMye_M$X == "Ifi207"
                                           |Unique_KOMyevsWTMye_M$X == "Ifnar1", TRUE,FALSE)

ggplot(Unique_KOMyevsWTMye_M) +
  theme_classic(base_size = 20)+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  geom_point(aes(log2FoldChange, -log10(pvalue), col = color)) +
  geom_text_repel(
    aes(log2FoldChange, -log10(pvalue)),
    label = ifelse(Unique_KOMyevsWTMye_M$genelabels, Unique_KOMyevsWTMye_M$X, ""),
    box.padding = unit(0.45, "lines"),
    hjust = 1,
    max.overlaps = 1000
  ) +
  theme(legend.title = element_blank(), text = element_text(size = 20))+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey","black"="black"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  theme(legend.position="none")+
  geom_vline(xintercept=c(0.1, -0.1),linetype=3)

##Fig 5G 
#graph pathways
GO_Unique_KOMyevsWTMye_M_dn <- read.csv("GO_Unique_KOMyevsWTMye_M.csv",header=T)

colnames(GO_Unique_KOMyevsWTMye_M_dn) <- c("name", "genes_in_gsea", "genes_in_data", "k/K", "p-value", "FDR")
GO_Unique_KOMyevsWTMye_M_dn$FDR <- as.numeric(GO_Unique_KOMyevsWTMye_M_dn$FDR)
GO_Unique_KOMyevsWTMye_M_dn$logFDR <- -log10(GO_Unique_KOMyevsWTMye_M_dn$FDR)
GO_Unique_KOMyevsWTMye_M_dn$enrich <- paste(GO_Unique_KOMyevsWTMye_M_dn$genes_in_data, "/", GO_Unique_KOMyevsWTMye_M_dn$genes_in_gsea, sep=" ")
GO_Unique_KOMyevsWTMye_M_dn <- GO_Unique_KOMyevsWTMye_M_dn[order(-GO_Unique_KOMyevsWTMye_M_dn$logFDR),]
GO_Unique_KOMyevsWTMye_M_dn$Name <-  gsub("GO_", "", GO_Unique_KOMyevsWTMye_M_dn$name)
GO_Unique_KOMyevsWTMye_M_dn$Name <- gsub("*_", " ", GO_Unique_KOMyevsWTMye_M_dn$Name)
GO_Unique_KOMyevsWTMye_M_dn$Name <-  factor(GO_Unique_KOMyevsWTMye_M_dn$Name, levels=rev(GO_Unique_KOMyevsWTMye_M_dn$Name))
area.color <- c("blue", "blue", "blue", "blue","blue", "blue", "blue", "blue", "blue", "blue")
#plot
ggplot(data=GO_Unique_KOMyevsWTMye_M_dn, aes(x=reorder(name, logFDR), y=logFDR)) +
  theme_classic() +
  ylab("log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity", fill=area.color) +
  geom_text(aes(label = enrich), vjust = 0.5, hjust = 1.25, colour = "white")+
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  theme(plot.title = element_text(hjust = -0.5))+
  ggtitle("GO_Unique_KOMyevsWTMye_M_dn")

##Fig 5H
Unique_KOMyevsWTMye_F$color[Unique_KOMyevsWTMye_F$log2FoldChange > 0.1 & Unique_KOMyevsWTMye_F$pvalue < 0.05] <- "red"
Unique_KOMyevsWTMye_F$color[Unique_KOMyevsWTMye_F$log2FoldChange < -0.1 & Unique_KOMyevsWTMye_F$pvalue < 0.05] <- "blue"

Unique_KOMyevsWTMye_F$genelabels <- ""
Unique_KOMyevsWTMye_F$genelabels <- ifelse(Unique_KOMyevsWTMye_F$X == "Psm13" 
                                           |Unique_KOMyevsWTMye_F$X == "Mcm6"
                                           |Unique_KOMyevsWTMye_F$X == "Top2a"
                                           |Unique_KOMyevsWTMye_F$X == "Fanca"
                                           |Unique_KOMyevsWTMye_F$X == "Washc1"
                                           |Unique_KOMyevsWTMye_F$X == "Rmi1"
                                           |Unique_KOMyevsWTMye_F$X == "Ddb1"
                                           |Unique_KOMyevsWTMye_F$X == "Topbp1"
                                           |Unique_KOMyevsWTMye_F$X == "Rad50"
                                           |Unique_KOMyevsWTMye_F$X == "Mcm5"
                                           |Unique_KOMyevsWTMye_F$X == "Tubg1"
                                           |Unique_KOMyevsWTMye_F$X == "Ncapd2"
                                           |Unique_KOMyevsWTMye_F$X == "Mlh1"
                                           |Unique_KOMyevsWTMye_F$X == "Mcm7"
                                           |Unique_KOMyevsWTMye_F$X == "Ppp4c"
                                           |Unique_KOMyevsWTMye_F$X == "Yy1"
                                           |Unique_KOMyevsWTMye_F$X == "Kat5"
                                           |Unique_KOMyevsWTMye_F$X == "Rmi2"
                                           |Unique_KOMyevsWTMye_F$X == "Rad51ap1"
                                           |Unique_KOMyevsWTMye_F$X == "Exd2"
                                           |Unique_KOMyevsWTMye_F$X == "Mcm4"
                                           |Unique_KOMyevsWTMye_F$X == "Fzr1"
                                           |Unique_KOMyevsWTMye_F$X == "Xrn2"
                                           , TRUE,FALSE)
ggplot(Unique_KOMyevsWTMye_F) +
  theme_classic(base_size = 20)+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  geom_point(aes(log2FoldChange, -log10(pvalue), col = color)) +
  geom_text_repel(
    aes(log2FoldChange, -log10(pvalue)),
    label = ifelse(Unique_KOMyevsWTMye_F$genelabels, Unique_KOMyevsWTMye_F$X, ""),
    box.padding = unit(0.45, "lines"),
    hjust = 1,
    max.overlaps = 18000
  ) +
  theme(legend.title = element_blank(), text = element_text(size = 20))+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey","black"="black"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  theme(legend.position="none")+
  geom_vline(xintercept=c(0.1, -0.1),linetype=3)

##Fig 5I
GO_Unique_KOMyevsWTMye_F_dn <- read.csv("GO_KOMyevsWTMye_F_Unique.csv",header=T)

colnames(GO_Unique_KOMyevsWTMye_F_dn) <- c("name", "genes_in_gsea", "genes_in_data", "k/K", "p-value", "FDR")
GO_Unique_KOMyevsWTMye_F_dn$FDR <- as.numeric(GO_Unique_KOMyevsWTMye_F_dn$FDR)
GO_Unique_KOMyevsWTMye_F_dn$logFDR <- -log10(GO_Unique_KOMyevsWTMye_F_dn$FDR)
GO_Unique_KOMyevsWTMye_F_dn$enrich <- paste(GO_Unique_KOMyevsWTMye_F_dn$genes_in_data, "/", GO_Unique_KOMyevsWTMye_F_dn$genes_in_gsea, sep=" ")
GO_Unique_KOMyevsWTMye_F_dn <- GO_Unique_KOMyevsWTMye_F_dn[order(-GO_Unique_KOMyevsWTMye_F_dn$logFDR),]
GO_Unique_KOMyevsWTMye_F_dn$Name <-  gsub("GO_", "", GO_Unique_KOMyevsWTMye_F_dn$name)
GO_Unique_KOMyevsWTMye_F_dn$Name <- gsub("*_", " ", GO_Unique_KOMyevsWTMye_F_dn$Name)
GO_Unique_KOMyevsWTMye_F_dn$Name <-  factor(GO_Unique_KOMyevsWTMye_F_dn$Name, levels=rev(GO_Unique_KOMyevsWTMye_F_dn$Name))
area.color <- c("blue", "blue", "blue", "blue","blue", "blue", "blue", "blue", "blue", "blue")
#plot
ggplot(data=GO_Unique_KOMyevsWTMye_F_dn, aes(x=reorder(name, logFDR), y=logFDR)) +
  theme_classic() +
  ylab("log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity", fill=area.color) +
  geom_text(aes(label = GO_Unique_KOMyevsWTMye_F_dn$enrich), vjust = 0.5, hjust = 1.25, colour = "white")+
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  theme(plot.title = element_text(hjust = -0.5))+
  ggtitle("GO_Unique_KOMyevsWTMye_F_dn")

