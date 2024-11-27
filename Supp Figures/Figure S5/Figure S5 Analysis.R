###Figure S5 Analysis
library(Seurat)
library(dplyr)
library(gplots)
library(RColorBrewer)

##Fig S5A
data <- readRDS("FCG_final.rds")
#UMAP by cell type
DimPlot(data, reduction = "umap",raster = F)

##Fig S5B 
#subset each cell type
Idents(data) = "celltype"
AST <- subset(data, idents = "astrocytes")
Epen <- subset(data, idents = "ependymal cells")
ExcitNeur <- subset(data, idents = "excitatory neurons")
InhibNeur <- subset(data, idents = "inhibitory neurons")
MG <- subset(data, idents = "microglia")
OL <- subset(data, idents = "oligodendrocytes")
OPC <- subset(data, idents = "OPCs")
Vasc <- subset(data, idents = "vascular cells")

#calculate # DEGs btwn CPZ vs. Ctrl for each cell type
##################read in RDS obj
AST <- readRDS("FCG_AST_final_res0.2.rds")
Epen <- readRDS("FCG_Epen_final_res0.2.rds") 
ExcitNeur <- readRDS("FCG_ExcitNeur_final_res0.2.rds") 
InhibNeur <- readRDS("FCG_InhibNeur_final_res0.2.rds") 
MG <- readRDS("FCG_MG_final_res0.2.rds") 
OL <- readRDS("FCG_OL_final_res0.2.rds") 
OPC <- readRDS("FCG_OPC_final_res0.2.rds") 
Vasc <- readRDS("FCG_Vasc_final_res0.2.rds") 
###########################
#set assay and idents for each cell type
DefaultAssay(AST) <- 'RNA'
DefaultAssay(Epen) <- 'RNA'
DefaultAssay(ExcitNeur) <- 'RNA'
DefaultAssay(InhibNeur) <- 'RNA'
DefaultAssay(MG) <- 'RNA'
DefaultAssay(OL) <- 'RNA'
DefaultAssay(OPC) <- 'RNA'
DefaultAssay(Vasc) <- 'RNA'

Idents(AST) <- "Food"
Idents(Epen) <- "Food"
Idents(ExcitNeur) <- "Food"
Idents(InhibNeur) <- "Food"
Idents(MG) <- "Food"
Idents(OL) <- "Food"
Idents(OPC) <- "Food"
Idents(Vasc) <- "Food"

#identify DEGs, export to .csv, subset only sig for each cell type
AST_DemyevsCtrl <- FindMarkers(AST, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(AST_DemyevsCtrl, "AST.DemyevsCtrl.csv")
AST_DemyevsCtrl <- subset(AST_DemyevsCtrl, AST_DemyevsCtrl$p_val_adj <0.05)

Epen_DemyevsCtrl <- FindMarkers(Epen, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(Epen_DemyevsCtrl, "Epen.DemyevsCtrl.csv")
Epen_DemyevsCtrl <- subset(Epen_DemyevsCtrl, Epen_DemyevsCtrl$p_val_adj <0.05)

ExcitNeur_DemyevsCtrl <- FindMarkers(ExcitNeur, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(ExcitNeur_DemyevsCtrl, "ExcitNeur.DemyevsCtrl.csv")
ExcitNeur_DemyevsCtrl <- subset(ExcitNeur_DemyevsCtrl, ExcitNeur_DemyevsCtrl$p_val_adj <0.05)

InhibNeur_DemyevsCtrl <- FindMarkers(InhibNeur, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(InhibNeur_DemyevsCtrl, "InhibNeur.DemyevsCtrl.csv")
InhibNeur_DemyevsCtrl <- subset(InhibNeur_DemyevsCtrl, InhibNeur_DemyevsCtrl$p_val_adj <0.05)

MG_DemyevsCtrl <- FindMarkers(MG, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(MG_DemyevsCtrl, "MG.DemyevsCtrl.csv")
MG_DemyevsCtrl <- subset(MG_DemyevsCtrl, MG_DemyevsCtrl$p_val_adj <0.05)

OL_DemyevsCtrl <- FindMarkers(OL, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(OL_DemyevsCtrl, "OL.DemyevsCtrl.csv")
OL_DemyevsCtrl <- subset(OL_DemyevsCtrl, OL_DemyevsCtrl$p_val_adj <0.05)

OPC_DemyevsCtrl <- FindMarkers(OPC, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(OPC_DemyevsCtrl, "OPC.DemyevsCtrl.csv")
OPC_DemyevsCtrl <- subset(OPC_DemyevsCtrl, OPC_DemyevsCtrl$p_val_adj <0.05)

Vasc_DemyevsCtrl <- FindMarkers(Vasc, ident.1 = "Demye", ident.2 = "Ctrl", logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = F)
write.csv(Vasc_DemyevsCtrl, "Vasc.DemyevsCtrl.csv")
Vasc_DemyevsCtrl <- subset(Vasc_DemyevsCtrl, Vasc_DemyevsCtrl$p_val_adj <0.05)
#can see DEG #s in R viewer, graphed in Prism

##Fig S5D
#identify DEGs in female CPZ vs. male CPZ for MG, AST, and OL
Idents(MG) <- "Genotype.Food"
Idents(OL) <- "Genotype.Food"
Idents(AST) <- "Genotype.Food"

MG_XXODemyevsXYTDemye <- FindMarkers(MG, ident.1 = "XXO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
MG_XXODemyevsXYTDemye <- subset(MG_XXODemyevsXYTDemye, MG_XXODemyevsXYTDemye$p_val_adj<0.05)
write.csv(MG_XXODemyevsXYTDemye, "MG_XXODemyevsXYTDemye_0.1.csv")

OL_XXODemyevsXYTDemye <- FindMarkers(OL, ident.1 = "XXO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
OL_XXODemyevsXYTDemye <- subset(OL_XXODemyevsXYTDemye, OL_XXODemyevsXYTDemye$p_val_adj<0.05)
write.csv(OL_XXODemyevsXYTDemye, "OL_XXODemyevsXYTDemye_0.1.csv")

AST_XXODemyevsXYTDemye <- FindMarkers(AST, ident.1 = "XXO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = F)
AST_XXODemyevsXYTDemye <- subset(AST_XXODemyevsXYTDemye, AST_XXODemyevsXYTDemye$p_val_adj<0.05)
write.csv(AST_XXODemyevsXYTDemye, "AST_XXODemyevsXYTDemye_0.1.csv")

#open each .csv in Excel, sort by log2FC, copy and paste names of top 25 up- and top 25 down-regulated DEGs into separate .csv
#did not include Xist (skewed)
MG_top25 <- read.csv("MG_top25DEGs.csv")
AST_top25 <- read.csv("AST_top25DEGs.csv")
OL_top25 <- read.csv("OL_top25DEGs.csv")

#now identify the DEGs in each gonad and sex chromosome comparison (within each cell type)
#then overlap with corresponding top25 for that cell type
#in this case we do not logfc.threshold bc we need to know the fc contribution of each comparison for each DEG (see calculation in S5C)
#MG gonad comparisons
MG_XYODemyevsXYTDemye <- FindMarkers(MG, ident.1 = "XYO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(MG_XYODemyevsXYTDemye, "MG_XYODemyevsXYTDemye.csv")
MG_XYODemyevsXYTDemye$Gene <- row.names(MG_XYODemyevsXYTDemye)
MG_XYODemyevsXYTDemye_Overlap <- MG_XYODemyevsXYTDemye %>%
  inner_join(MG_top25, by = "Gene")
write.csv(MG_XYODemyevsXYTDemye_Overlap, "MG_XYODemyevsXYTDemye_top25.csv")

MG_XXTDemyevsXXODemye <- FindMarkers(MG, ident.1 = "XXT_Demye", ident.2 = "XXO_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(MG_XXTDemyevsXXODemye, "MG_XXTDemyevsXXODemye.csv")
MG_XXTDemyevsXXODemye$Gene <- row.names(MG_XXTDemyevsXXODemye)
MG_XXTDemyevsXXODemye_Overlap <- MG_XXTDemyevsXXODemye %>%
  inner_join(MG_top25, by = "Gene")
write.csv(MG_XXTDemyevsXXODemye_Overlap, "MG_XXTDemyevsXXODemye_top25.csv")

#MG sex chromosome comparisons
MG_XYODemyevsXXODemye <- FindMarkers(MG, ident.1 = "XYO_Demye", ident.2 = "XXO_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(MG_XYODemyevsXXODemye, "MG_XYODemyevsXXODemye.csv")
MG_XYODemyevsXXODemye$Gene <- row.names(MG_XYODemyevsXXODemye)
MG_XYODemyevsXXODemye_Overlap <- MG_XYODemyevsXXODemye %>%
  inner_join(MG_top25, by = "Gene")
write.csv(MG_XYODemyevsXXODemye_Overlap, "MG_XXODemyevsXYODemye_top25.csv")

MG_XXTDemyevsXYTDemye <- FindMarkers(MG, ident.1 = "XXT_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(MG_XXTDemyevsXYTDemye, "MG_XXTDemyevsXYTDemye.csv")
MG_XXTDemyevsXYTDemye$Gene <- row.names(MG_XXTDemyevsXYTDemye)
MG_XXTDemyevsXYTDemye_Overlap <- MG_XXTDemyevsXYTDemye %>%
  inner_join(MG_top25, by = "Gene")
write.csv(MG_XXTDemyevsXYTDemye_Overlap, "MG_XXTDemyevsXYTDemye_top25.csv")

#repeat for OL
OL_XYODemyevsXYTDemye <- FindMarkers(OL, ident.1 = "XYO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(OL_XYODemyevsXYTDemye, "OL_XYODemyevsXYTDemye.csv")
OL_XYODemyevsXYTDemye$Gene <- row.names(OL_XYODemyevsXYTDemye)
OL_XYODemyevsXYTDemye_Overlap <- OL_XYODemyevsXYTDemye %>%
  inner_join(OL_top25, by = "Gene")
write.csv(OL_XYODemyevsXYTDemye_Overlap, "OL_XYODemyevsXYTDemye_top25.csv")

OL_XXODemyevsXXTDemye <- FindMarkers(OL, ident.1 = "XXO_Demye", ident.2 = "XXT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(OL_XXODemyevsXXTDemye, "OL_XXODemyevsXXTDemye.csv")
OL_XXODemyevsXXTDemye$Gene <- row.names(OL_XXODemyevsXXTDemye)
OL_XXODemyevsXXTDemye_Overlap <- OL_XXODemyevsXXTDemye %>%
  inner_join(OL_top25, by = "Gene")
write.csv(OL_XXODemyevsXXTDemye_Overlap, "OL_XXODemyevsXXTDemye_top25.csv")

OL_XXODemyevsXYODemye <- FindMarkers(OL, ident.1 = "XXO_Demye", ident.2 = "XYO_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(OL_XXODemyevsXYODemye, "OL_XXODemyevsXYODemye.csv")
OL_XXODemyevsXYODemye$Gene <- row.names(OL_XXODemyevsXYODemye)
OL_XXODemyevsXYODemye_Overlap <- OL_XXODemyevsXYODemye %>%
  inner_join(OL_top25, by = "Gene")
write.csv(OL_XXODemyevsXYODemye_Overlap, "OL_XXODemyevsXYODemye_top25.csv")

OL_XXTDemyevsXYTDemye <- FindMarkers(OL, ident.1 = "XXT_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(OL_XXTDemyevsXYTDemye, "OL_XXTDemyevsXYTDemye.csv")
OL_XXTDemyevsXYTDemye$Gene <- row.names(OL_XXTDemyevsXYTDemye)
OL_XXTDemyevsXYTDemye_Overlap <- OL_XXTDemyevsXYTDemye %>%
  inner_join(OL_top25, by = "Gene")
write.csv(OL_XXTDemyevsXYTDemye_Overlap, "OL_XXTDemyevsXYTDemye_top25.csv")

#repeat for AST
AST_XYODemyevsXYTDemye <- FindMarkers(AST, ident.1 = "XYO_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(AST_XYODemyevsXYTDemye, "AST_XYODemyevsXYTDemye.csv")
AST_XYODemyevsXYTDemye$Gene <- row.names(AST_XYODemyevsXYTDemye)
AST_XYODemyevsXYTDemye_Overlap <- AST_XYODemyevsXYTDemye %>%
  inner_join(AST_top25, by = "Gene")
write.csv(AST_XYODemyevsXYTDemye_Overlap, "AST_XYODemyevsXYTDemye_top25.csv")

AST_XXODemyevsXXTDemye <- FindMarkers(AST, ident.1 = "XXO_Demye", ident.2 = "XXT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(AST_XXODemyevsXXTDemye, "AST_XXODemyevsXXTDemye.csv")
AST_XXODemyevsXXTDemye$Gene <- row.names(AST_XXODemyevsXXTDemye)
AST_XXODemyevsXXTDemye_Overlap <- AST_XXODemyevsXXTDemye %>%
  inner_join(AST_top25, by = "Gene")
write.csv(AST_XXODemyevsXXTDemye_Overlap, "AST_XXODemyevsXXTDemye_top25.csv")

AST_XXODemyevsXYODemye <- FindMarkers(AST, ident.1 = "XXO_Demye", ident.2 = "XYO_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(AST_XXODemyevsXYODemye, "AST_XXODemyevsXYODemye.csv")
AST_XXODemyevsXYODemye$Gene <- row.names(AST_XXODemyevsXYODemye)
AST_XXODemyevsXYODemye_Overlap <- AST_XXODemyevsXYODemye %>%
  inner_join(AST_top25, by = "Gene")
write.csv(AST_XXODemyevsXYODemye_Overlap, "AST_XXODemyevsXYODemye_top25.csv")

AST_XXTDemyevsXYTDemye <- FindMarkers(AST, ident.1 = "XXT_Demye", ident.2 = "XYT_Demye", logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
write.csv(AST_XXTDemyevsXYTDemye, "AST_XXTDemyevsXYTDemye.csv")
AST_XXTDemyevsXYTDemye$Gene <- row.names(AST_XXTDemyevsXYTDemye)
AST_XXTDemyevsXYTDemye_Overlap <- AST_XXTDemyevsXYTDemye %>%
  inner_join(AST_top25, by = "Gene")
write.csv(AST_XXTDemyevsXYTDemye_Overlap, "AST_XXTDemyevsXYTDemye_top25.csv")

#made a spreadsheet in Excel where copied+pasted the log2fc of each comparison for overlapping DEGs into a tab for each cell type (Table S1)
#used Excel to calculate % sex chromo and % gonad expression composing XXO CPZ vs XYT CPZ
#then Excel to calculate ratio of sex chromosome:gonad influence, did stringent cutoff of ratio<0.5 = gonad dominant, ratio>1.1 = chromo dominant, 0.5<ratio<1.1 = chromo-gonad interaction 
#Then in Excel, ordered DEGs by cell type, 1 col= gonad dominant, 1 col=sex chromo dominant, 1 col=sex chromo-gonad interaction. Used binary matrix where in each row/DEG, 0=not dominant, 1=dominant
#read in to make heatmap 
SexChromoGonadHeatmap <- read.csv("SexChromoGonad_Heatmap.csv", row.names = 1)
SexChromoGonadHeatmap <- as.matrix(SexChromoGonadHeatmap)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
heatmap.2(SexChromoGonadHeatmap, cellnote = SexChromoGonadHeatmap, notecol = "black", trace = "none", scale = "none", dendrogram = "none", Rowv = F, density.info = "none", col =  col, Colv = T, labCol = T)
#changed colors in Illustrator, green for not score dominant, orange = score dominant

##Fig S5E
#Identify marker genes for each cell type and keep only X-linked genes 
All <- readRDS("FCG_final.rds")
DefaultAssay(All) <- 'RNA'
Idents(All) <- "celltype"

Xgenes <- read.csv("XlinkedGenes_Mouse.csv")
Xgenes$gene <- Xgenes$Gene_ID

mg <- FindMarkers(All, ident.1 = "microglia",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
mg$gene <- row.names(mg)
overlap.MG_Xgenes <- mg %>%
  inner_join(Xgenes, by = "gene")

ol <- FindMarkers(All, ident.1 = "oligodendrocytes",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
ol$gene <- row.names(ol)
overlap.OL_Xgenes <- ol %>%
  inner_join(Xgenes, by = "gene")

opc <- FindMarkers(All, ident.1 = "OPCs",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
opc$gene <- row.names(opc)
overlap.OPC_Xgenes <- opc %>%
  inner_join(Xgenes, by = "gene")

ast <- FindMarkers(All, ident.1 = "astrocytes",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
ast$gene <- row.names(ast)
overlap.AST_Xgenes <- ast %>%
  inner_join(Xgenes, by = "gene")

epen <- FindMarkers(All, ident.1 = "ependymal cells",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
epen$gene <- row.names(epen)
overlap.EPEN_Xgenes <- epen %>%
  inner_join(Xgenes, by = "gene")

IN <- FindMarkers(All, ident.1 = "inhibitory neurons",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
IN$gene <- row.names(IN)
overlap.IN_Xgenes <- IN %>%
  inner_join(Xgenes, by = "gene")

VC <- FindMarkers(All, ident.1 = "vascular cells",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
VC$gene <- row.names(VC)
overlap.VC_Xgenes <- VC %>%
  inner_join(Xgenes, by = "gene")

EN <- FindMarkers(All, ident.1 = "excitatory neurons",logfc.threshold = 0, test.use = "MAST",min.pct = 0.25, only.pos = F)
EN$gene <- row.names(EN)
overlap.EN_Xgenes <- EN %>%
  inner_join(Xgenes, by = "gene")

#change col title for log2FC to cell type so you can compile all cell types into 1 spreadsheet 
overlap.AST_Xgenes$AST <- overlap.AST_Xgenes$avg_log2FC
overlap.OL_Xgenes$OL <- overlap.OL_Xgenes$avg_log2FC
overlap.MG_Xgenes$MG <- overlap.MG_Xgenes$avg_log2FC
overlap.IN_Xgenes$IN <- overlap.IN_Xgenes$avg_log2FC
overlap.EN_Xgenes$EN <- overlap.EN_Xgenes$avg_log2FC
overlap.EPEN_Xgenes$EPEN <- overlap.EPEN_Xgenes$avg_log2FC

#read out to Excel to format into 1 spreadsheet with rows = Gene names, cols = cell type log2FC
write.csv(overlap.OL_Xgenes, "overlap.AST_Xgenes_CellType.csv")
write.csv(overlap.AST_Xgenes, "overlap.OL_Xgenes_CellType.csv")
write.csv(overlap.MG_Xgenes, "overlap.MG_Xgenes_CellType.csv")
write.csv(overlap.IN_Xgenes, "overlap.IN_Xgenes_CellType.csv")
write.csv(overlap.EN_Xgenes, "overlap.EN_Xgenes_CellType.csv")
write.csv(overlap.EPEN_Xgenes, "overlap.EPEN_Xgenes_CellType.csv")

#read formatted spreadsheet back in, if gene was not a cell type marker for a given cell type, log2fc=0
XgenesCellType <- read.csv("Xgenes_CellType_CellMarkers.csv", row.names = 1)
XgenesCellType.matrix <- as.matrix(XgenesCellType)
#heatmap
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(60))
heatmap.2(XgenesCellType.matrix, trace = "none", scale = "none", dendrogram = "none", Rowv = F, density.info = "none", col =  col, breaks=seq(-3,3,0.1), Colv = F)

#now manually search in Excel if each DEG is also an XXO CPZ vs. XYT CPZ DEG (Table S1)
#if yes, add black dot for corresponding cell type (Illustrator)
