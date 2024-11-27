###Figure 2 Analysis
#5wk CPZ snRNAseq

##Fig 2A-B
library(Seurat)
library(MAST)
library(RColorBrewer)
MG <- readRDS("Chloe_FCG_MG_reclusted_res0.2.rds")

DefaultAssay(MG) <- 'RNA'
Idents(MG) <- "Condition"
XXODemyevsXYTDemye <- FindMarkers(MG, ident.1 = "XXO_Cup", ident.2 = "XYT_Cup", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
XXODemyevsXYTDemye <- subset(XXODemyevsXYTDemye, XXODemyevsXYTDemye$p_val_adj<0.05)
write.csv(XXODemyevsXYTDemye, "MG_XXODemyevsXYTDemye_0.1.csv")
XXODemyevsXYTDemye$Gene <- row.names(XXODemyevsXYTDemye)

#overlap with 3 wk MG to see if sex differences are conserved
EarlyMG_XXODemyevsXYTDemye <- read.csv("3wk_MG_XXODemyevsXYTDemye_0.1.csv")
MG_XXODemyevsXYTDemye_Overlap <- XXODemyevsXYTDemye %>%
  inner_join(EarlyMG_XXODemyevsXYTDemye, by = "Gene")
write.csv(MG_XXODemyevsXYTDemye_Overlap, "MG_XXODemyevsXYTDemye_3wk5wk_Overlap.csv")

#reformat for heatmap in Excel so that genes=rows and cols = log2FC in 3 wk or 5wk
#read in and make heatmap
OverlapHeatmap <- read.csv("MG_3wk5wkOverlap_Hformat.csv", row.names = 1)
#had to remove Tsix, skewing heatmap 
OverlapHeatmap = OverlapHeatmap[-1,]
#also Xist, skewing heatmap
OverlapHeatmap = OverlapHeatmap[-1,]

OverlapHeatmap.matrix <- as.matrix(OverlapHeatmap)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
heatmap(OverlapHeatmap.matrix, scale = "none", col =  col, Colv = NA)

##Fig. 2E 
Idents(MG) <- "seurat_clusters"

#rename clusters to start with 1
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)

MG2 <- FindMarkers(MG, ident.1= "2", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
MG2 <- subset(MG2, MG2$p_val_adj<0.05)
write.csv(MG2, "MG_Clust2DEGs_0.1.csv")
MG2$gene <- row.names(MG2)

DAM <- read.csv("DAM_KerenShaul_0.1.csv")
overlap.MG2_DAM <- MG2 %>%
  inner_join(DAM, by = "gene")

#export and format
write.csv(overlap.MG2_DAM, "overlap.MG2_DAM.csv")

#format so that 1 col=gene, 1 col=MG2 logFC, 1 col=DAM
overlap.MG2_DAM.df <- read.csv("overlap.MG2_DAM_df2.csv", row.names = 1)
overlap.MG2_DAM.df<- as.data.frame(overlap.MG2_DAM.df)
overlap.MG2_DAM.df$gene <- rownames(overlap.MG2_DAM.df)

#correlate
ggscatter(overlap.MG2_DAM.df, x = "MG2", y = "DAM", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MG2", ylab = "DAM",
          repel=TRUE)+
  geom_text_repel(label = overlap.MG2_DAM.df$gene, max.overlaps = 10)+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey")

##Fig 2F-G
OL <- readRDS("Chloe_FCG_OL_reclusted_res0.2.rds")
DefaultAssay(OL) <- 'RNA'
Idents(OL) <- "Condition"

OL_XXODemyevsXYTDemye <- FindMarkers(OL, ident.1 = "XXO_Cup", ident.2 = "XYT_Cup", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
OL_XXODemyevsXYTDemye <- subset(OL_XXODemyevsXYTDemye, OL_XXODemyevsXYTDemye$p_val_adj<0.05)
OL_XXODemyevsXYTDemye$Gene <- rownames(OL_XXODemyevsXYTDemye)
write.csv(OL_XXODemyevsXYTDemye, "OL_XXODemyevsXYTDemye_0.1.csv")

#overlap with same comparison in 3 wk OL to see overlapping changes
#read in 3wk DEG list and overlap
EarlyOL_XXODemyevsXYTDemye <- read.csv("3wk_OL_XXODemyevsXYTDemye.csv")
EarlyOL_XXODemyevsXYTDemye$Gene <- EarlyOL_XXODemyevsXYTDemye$X
XXODemyevsXYTDemye_Overlap <- OL_XXODemyevsXYTDemye %>%
  inner_join(EarlyOL_XXODemyevsXYTDemye, by = "Gene")
write.csv(XXODemyevsXYTDemye_Overlap, "XXODemyevsXYTDemye_3wk5wk_Overlap.csv")

#reformat for heatmap in Excel so that genes=rows and cols = log2FC in 3 wk or 5wk
#read in and make heatmap
OverlapHeatmap <- read.csv("OL_3wk5wkOverlap_Hformat.csv", row.names = 1)
#had to remove Xist again, was skewing heatmap
OverlapHeatmap = OverlapHeatmap[-49,]
OverlapHeatmap.matrix <- as.matrix(OverlapHeatmap)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
heatmap(OverlapHeatmap.matrix, scale = "none", col =  col, Colv = NA)

##Fig. 2H
#identify XX vs. XY DEGs in CPZ
#subset out so have object with CPZ groups only
Idents(OL) <- "Condition"
Demye_OL_CPZ <- subset(OL, idents = c("XXO_Cup", "XXT_Cup", "XYO_Cup", "XYT_Cup"))

#now find XX vs XY DEGs
Idents(Demye_OL_CPZ) <- "Sex.Chromosome"
Demye_OL_CPZ_XXvsXY <- FindMarkers(Demye_OL_CPZ, ident.1 = "XX", ident.2 = "XY", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
Demye_OL_CPZ_XXvsXY <- subset(Demye_OL_CPZ_XXvsXY, Demye_OL_CPZ_XXvsXY$p_val_adj<0.05)
Demye_OL_CPZ_XXvsXY$Gene <- rownames(Demye_OL_CPZ_XXvsXY)

#repeat for Ovaries vs. Testes
Idents(Demye_OL_CPZ) <- "Gonad"
Demye_OL_CPZ_OvsT <- FindMarkers(Demye_OL_CPZ, ident.1 = "Ovaries", ident.2 = "Testes", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
Demye_OL_CPZ_OvsT <- subset(Demye_OL_CPZ_OvsT, Demye_OL_CPZ_OvsT$p_val_adj<0.05)
Demye_OL_CPZ_OvsT$Gene <- rownames(Demye_OL_CPZ_OvsT)

#overlap XXOCPZvsXYTCPZ with sex chromo DEGs and gonad DEGs
OL_XXOCPZvsXYTCPZ_SexChromoOverlap <- OL_XXODemyevsXYTDemye %>%
  inner_join(Demye_OL_CPZ_XXvsXY, by = "Gene")
write.csv(OL_XXOCPZvsXYTCPZ_SexChromoOverlap, "OL_XXOCPZvsXYTCPZ_SexChromoOverlap.csv")

OL_XXOCPZvsXYTCPZ_GonadOverlap <- OL_XXODemyevsXYTDemye %>%
  inner_join(Demye_OL_CPZ_OvsT, by = "Gene")
write.csv(OL_XXOCPZvsXYTCPZ_GonadOverlap, "OL_XXOCPZvsXYTCPZ_GonadOverlap.csv")

#now overlap chromo overlap and gonad overlap to find DEGs overlapping between all 3 (sex chromo, gonad, and XXOCPZvsXYTCPZ)
OL_XXOCPZvsXYTCPZ_SexChromoGonadOverlap <- OL_XXOCPZvsXYTCPZ_SexChromoOverlap %>%
  inner_join(OL_XXOCPZvsXYTCPZ_GonadOverlap, by = "Gene")
write.csv(OL_XXOCPZvsXYTCPZ_SexChromoGonadOverlap, "OL_XXOCPZvsXYTCPZ_SexChromoGonadOverlap.csv")

##Fig. 2I
#Cellchat
install.packages("remotes")
remotes::install_github("sqjin/CellChat", force = TRUE)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
#need to use RDS obj with all cell types
All <- readRDS("FCG_integrated_Annotation.rds")
#only want CPZ groups
Idents(All) <- "Condition"
Demye_All <- subset(All, idents = c("XXO_Cup", "XYO_Cup", "XXT_Cup", "XYT_Cup"))

#proceed with Cellchat
Demye_All[[]]
Idents(Demye_All) <- "celltype"

data.input <- GetAssayData(Demye_All, assay = "RNA", layer = "data") # normalized data matrix
labels <- Idents(Demye_All)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
##below is optional, only if you run into memory error from computeCommunProb
options(future.globals.maxSize = 8.074e+8)
##continue 
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#create df with all inferred cell-cell communications
df.net <- subsetCommunication(cellchat)

#infer cell-cell comm at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#calculate aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

#show all the significant interactions (L-R pairs) coming from OL only (defined by 'sources.use') to other cell types (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:8), remove.isolate = FALSE)

#manually colored the ligand-receptor pairs that named at least 1 DEG from (OL) XXO vs XYT DEG list or Ovaries vs Testes DEG list
#created graph with only sex-biased DEGs in Illustrator
#noticed that Cdh2-Cdh2 signaled to multiple cell types and was gonad DEG
#visualize CDH signaling in XXO CPZ and XYT CPZ, subset out into their own RDS obj
Idents(Demye_All) <- "Condition"
Demye_All_XXOCPZ <- subset(Demye_All, idents = "XXO_Cup")
Demye_All_XXOCPZ[[]]

Idents(Demye_All_XXOCPZ) <- "celltype"
data.input <- GetAssayData(Demye_All_XXOCPZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Demye_All_XXOCPZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#visualize CDH signaling chord pathway
pathways.show <- "CDH" 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#repeat for XYT CPZ
Demye_All_XYTCPZ <- subset(Demye_All, idents = "XYT_Cup")
Demye_All_XYTCPZ[[]]
Idents(Demye_All_XYTCPZ) <- "celltype"
data.input <- GetAssayData(Demye_All_XYTCPZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Demye_All_XYTCPZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#visualize CDH signaling chord pathway
pathways.show <- c("CDH") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#repeat for gonad comparison, Ovaries first
Idents(Demye_All) <- "Condition"
Demye_All_OvariesCPZ <- subset(Demye_All, idents = c("XXO_Cup", "XYO_Cup"))

Demye_All_OvariesCPZ[[]]
Idents(Demye_All_OvariesCPZ) <- "celltype"
data.input <- GetAssayData(Demye_All_OvariesCPZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Demye_All_OvariesCPZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#visualize CDH signaling chord pathway
pathways.show <- c("CDH") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#now testes
Demye_All_TestesCPZ <- subset(Demye_All, idents = c("XXT_Cup", "XYT_Cup"))

Demye_All_TestesCPZ[[]]
Idents(Demye_All_TestesCPZ) <- "celltype"
data.input <- GetAssayData(Demye_All_TestesCPZ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Demye_All_TestesCPZ)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#visualize CDH signaling chord pathway
pathways.show <- c("CDH") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
