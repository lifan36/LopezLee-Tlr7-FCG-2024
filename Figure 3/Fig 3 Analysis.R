###Fig. 3 (protein-level, no sequencing)
library(edgeR)
library(limma)

##Fig. 3G-J pipeline repeated for each panel, bar graphs were made in Prism
#calculate avg negative concentration of ctrl samples within each cytokine in Excel
#read in
data <- read.csv("MagPix_RList.csv")
data$PureFoldChange <- (data$Concentration)/data$Avg.Cupri.Neg.All
data$log2transform <- log2(data$PureFoldChange)
#Compare btwn CPZ groups only
data <- subset(data, data$Treatment=="Cuprizone")
#write out and format in Excel so that each row=cytokine, col=each sample
write.csv(data, "MagPix_RList_FC.csv")
#read in formatted
data <- read.csv("MagPix_RList_Heatmap.csv", row.names = 1)
#read in sample names + genotype/any treatment info (metadata)
names <- read.csv("Names.csv",header=T) 
#differential expression
mm <- model.matrix(~0 + names$Gonad, data = data) #build model matrix 
colnames(mm) <- gsub("^.*\\$","",colnames(mm))
fit <- lmFit(data, mm) #linearfit
head(coef(fit))
contr <- makeContrasts(Condition="GonadOvaries - GonadTestes", levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
demg<-topTable(tmp, sort.by = "P", n = "inf")
write.csv(demg,"DE_OvsTCPZ_log2FC.csv") #output differential expression analysis
#start from beginning and repeat for Fig 3H, 3I, 3J

##Fig. 3K
#format the data for heatmap, rows = cytokines, cols = each sample
#read in only XYO, XYT samples
data2 <- read.csv("DemyeOnly_XYOvsXYT.csv", row.names=1)
data2.matrix <- as.matrix(data2)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
#heatmap
heatmap(data2.matrix, scale = "row", col =  col,  Colv = NA)

##Fig. 3L
#format the data for heatmap, rows = cytokines, cols = each sample
#only XXO, XXT samples
data2 <- read.csv("DemyeOnly_XXOvsXXT.csv", row.names = 1)
data2.matrix <- as.matrix(data2)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
#heatmap
heatmap(data2.matrix, scale = "row", col =  col, Colv = NA)

