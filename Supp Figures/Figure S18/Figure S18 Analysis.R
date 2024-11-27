###Fig S18

##Fig S18A, same dataset as Fig 6 (TLR7 inhibitor)
data <- read.csv("counts_genes_sorted_PSUK.csv", header=T, row.names=1)

#metadata
filtered = data[rowSums(data)>15,] #removes any genes with <15 counts
sample <- colnames(data)
sex <- c(rep("F_",3),rep("M_",3), rep("F_",5), rep("M_",6))
treatment <- c(rep("Ctrl",6),rep("Inhib",11))
sex_treatment <- paste(sex,treatment,sep="_")
meta <- data.frame(sample=sample, sex=sex, treatment=treatment, sex_treatment=sex_treatment)
all(colnames(filtered) %in% meta$sample)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~sex_treatment)
rld <- rlog(dds, blind = T)

#plot
plotPCA(rld, intgroup = "sex_treatment", ntop = 50)+
  theme_classic()

##Fig S18B
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
