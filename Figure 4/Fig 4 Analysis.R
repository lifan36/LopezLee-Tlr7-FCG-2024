###Figure 4 (FCG;P301S;APOE4 and P301S;APOE4)
library(ggpubr)
library(Seurat)

##Fig. 4E
LG93E4_MG <- readRDS("LG93E4_MG_clustered.rds")
DefaultAssay(LG93E4_MG) <- "RNA"

# calculate ratio of each cell type cluster in each sample
a<-as.data.frame(table(LG93E4_MG$seurat_clusters,LG93E4_MG$Sample_Name))
colnames(a)<-c("clusters","genotype","cell.no")
agg<-aggregate(cell.no~genotype,a,sum)
a$genotype.total <- agg$cell.no[match(a$genotype,agg$genotype)]
a$ratio<-a$cell.no/a$genotype.total
#make groups spreadsheet in Excel listing out genotype for each subcluster+sample (i.e. if 5 subclusters x 3 samples = 15 rows in that genotype group)
groups <- read.csv("Groups.csv", header = 1) #read in groups
a$groups <- groups$Group
write.csv(a, "MG_CellRatio.csv")
#calculate in Excel avg and sem, read back in
a <- read.csv("MG_CellRatio.csv", header=1)
#order bars on graph
a$gclass <- factor(a$groups, levels=c('XXO_NTG_E4', 'XXO_P301S_E4', 'XYO_NTG_E4', 'XYO_P301S_E4', 'XXT_NTG_E4', 'XXT_P301S_E4', 'XYT_NTG_E4', 'XYT_P301S_E4'))
#plot
ggplot(a, aes(clusters, ratio, fill=groups))+
  geom_errorbar( aes(ymin=Avg-sem, ymax=Avg+sem), colour="black", position = "dodge")+
  geom_bar(position = "dodge", stat="summary", fun = "mean")+
  theme_classic()+
  scale_fill_brewer(palette="RdBu")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  xlab("Cluster")+
  ylab("Cell no./Condition Total")+
  geom_point(position=position_dodge(width=0.9))

##Fig 4F
Idents(LG93E4_MG) <- "seurat_clusters"

#renumber to start with 1
n <- dim(table(LG93E4_MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG93E4_MG@active.ident <- plyr::mapvalues(x = LG93E4_MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
LG93E4_MG@active.ident <- factor(LG93E4_MG@active.ident, levels=1:n)

#markers
clust_markers <- FindAllMarkers(LG93E4_MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(clust_markers, "markers_cluster.csv")

#copy clust 3 markers from above, paste in new spreadsheet in Excel, filter for only p_val_adj<0.05, read in
MG3 <- read.csv("clust3_markers.csv")
#correlate with DAM
DAM <- read.csv("DAM_KerenShaul_0.1.csv")
#overlap
MG3_DAM_overlap <- merge(MG3, DAM, by = "gene")
#write out and format so that 1 col = gene, 1 col = MG2 avg_logFC, 1 col=DAM
write.csv(MG3_DAM_overlap, "MG3_DAM_overlap.csv")
MG3_DAM_Corr <- read.csv("MG3_DAM_Corr.csv")
ggscatter(MG3_DAM_Corr, x = "MG3", y = "DAM", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MG3", ylab = "DAM Keren Shaul et al.",
          repel=TRUE)+
  geom_text_repel(label = MG3_DAM_Corr$gene, max.overlaps = 10)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")

##Fig 4I
#Identify pseudobulk DEGs (DEG list in Table S8)
Idents(LG93E4_MG) <- "Condition"
XXOTauvsXYTTau <- FindMarkers(LG93E4_MG, ident.1 = "XXO_P301S_E4", ident.2 = "XYT_P301S_E4", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(XXOTauvsXYTTau, "XXOTauvsXYTTau_0.1.csv")

#copy and paste DEGs into IPA software (set logFC cutoff=0.1, adj_p_val<0.05)
#graph IPA upstream reg results
#read in output
ipa<-read.csv("IPA_UpstreamRegs_MG_XYTPS19vsXXOPS19.csv",header=T)
#order by z score
ipa <- ipa[order(ipa$Activation.z.score), ]
#plot
ggplot(data=ipa, aes(x=reorder(Upstream.Regulator,-log(p.value.of.overlap)), y= -log(p.value.of.overlap), label=round(Activation.z.score, digits=1))) + 
  geom_point(stat='identity', aes(col=Predicted.Activation.State, size=abs(Activation.z.score)))  +
  scale_size_continuous(range = c(5, 8))+
  scale_color_manual(name="Predicted Activation State", 
                     labels = c("Activated", "Inhibited"), 
                     values = c("Activated"="tomato1", "Inhibited"="steelblue3")) + 
  geom_text(color="white", size=2.5) +
  facet_wrap(~Predicted.Activation.State, scales='free_y',ncol=2 )+
  labs(title="XYT PS19+ vs. XXO PS19+ Upstream Regulators") +
  theme_bw()+
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()

##Fig 4N
#read in dataset of E4 P301S (no FCG)
E4_noFCG_MG <- readRDS("MG.rds")
DefaultAssay(E4_noFCG_MG) <- "RNA"
Idents(E4_noFCG_MG) <- "Condition"

#dotplot
DotPlot(E4_noFCG_MG, split.by = "Condition", features = c("Trim30a", "Herc6", "Ddx60", "Apobec3", "Gm4951", "Rnf213","Parp14", "Sp100", 
                                                          "Trim30d", "Stat1", "Slfn8", 
                                                          "Nlrc5", "Eif2ak2", "Cd300lf","Stat2", "Irf9", "Isg15",
                                                          "Ifit1", "Oas1a"), cols="RdBu", scale.by = "size") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))


