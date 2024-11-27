##Fig S13
library(DESeq2)
library(pheatmap)

#Fig S13A, same dataset as Fig 5A-C (TLR7 Agonist)
data <- read.csv("counts_genes_sorted.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,] #removes any genes with <15 counts
sample <- colnames(data)
sex <- c(rep("M_",14),rep("F_",14))
treatment <- c(rep("Saline",7),rep("GDQ",7), rep("Saline", 7), rep("GDQ",7))
sex_treatment <- paste(sex,treatment,sep="_")
meta <- data.frame(sample=sample, sex=sex, treatment=treatment, sex_treatment=sex_treatment)
all(colnames(filtered) %in% meta$sample)
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment)
rld <- rlog(dds, blind = T)
#plot PCA
plotPCA(rld, intgroup = "sex_treatment", ntop = 50)+
  theme_classic()

##Fig S13B
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sex_treatment)
colnames(sampleDistMatrix) <- paste(vsd$sample, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(20)
#plot
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

##Fig S13C, same dataset as Fig 5E-I, (TLR7KO)
data <- read.csv("counts_genes_sorted.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,] #removes any genes with <15 counts

sample <- colnames(data)
sex <- c(rep("M_",12),rep("F_",12))
treatment <- c(rep("DMSO",3),rep("Mye",3), rep("DMSO", 3), rep("Mye",3), rep("DMSO", 3), rep("Mye", 3), rep("DMSO", 3), rep("Mye", 3))
genotype <- c(rep("WT",6),rep("KO",6), rep("WT", 6), rep("KO",6))
sex_treatment_genotype <- paste(sex,treatment,genotype,sep="_")
meta <- data.frame(sample=sample, sex=sex, genotype = genotype, treatment=treatment, sex_treatment_genotype=sex_treatment_genotype)

all(colnames(filtered) %in% meta$sample)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment_genotype)
rld <- rlog(dds, blind = T)
#plot
plotPCA(rld, intgroup = "sex_treatment_genotype", ntop = 500)+
  theme_classic()

##Fig S13D
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sex_treatment_genotype)
colnames(sampleDistMatrix) <- paste(vsd$sample, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(20)
#plot
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


