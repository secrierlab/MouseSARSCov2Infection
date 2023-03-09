library(DESeq2)
library(gdata)
library(ggplot2)
library(pheatmap)

gset <- read.delim("RNAseq Data C52/raw_counts.tsv")
gset <- gset[which(!duplicated(gset$gene)),]
df <- gset[,-c(1:2)]
rownames(df) <- gset$gene

ids <- read.xls("RNAseq Data C52/S5962_QC_results_27.6.22_bmer.xlsx")
ids$Day <- sapply(ids$Externe.ID, function(x) strsplit(x,"_")[[1]][2])
ids$Condition <- sapply(ids$Externe.ID, function(x) ifelse(grepl("E",x),"Control","Infected"))
ids$Condition <- factor(ids$Condition)
ids$Group <- apply(ids[,c("Day","Condition")],1,
                   function(x) ifelse(x[2] == "Control","Control",
                                      paste0(x[2],"_",x[1])))
ids$Group <- factor(ids$Group)
colnames(ids)[1] <- "Sample"
rownames(ids) <- ids$Sample

ids$Group <- relevel(ids$Group, "Control")


ids <- ids[which(ids$Sample %in% colnames(df)),]
df <- df[,ids$Sample] 
#Next, construct a DESeqDataSet data object:
  dds <- DESeqDataSetFromMatrix(countData = df, colData = ids, design = ~ Condition)



#  1. Visual exploration of the data

#Check counts on dds:
head(counts(dds))

#Remove genes with less than 10 reads in all conditions. This will help DeSeq2 to perform faster. Furthermore, DeSeq2 DGE on lowly and highly expressed genes doesn't perform that well so it's better to eliminate them.

dds <- dds[ rowSums(counts(dds)) > 10, ]

# Check counts on dds after filtering
head(counts(dds))

## Normalisation

### Sequencing depth normalisation

#Investigate different sequencing library sizes. DGE software will use this sequencing depth information to normalise between samples.
colSums(counts(dds)) 


#Calculate the size factor and add it to the data set.
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## Data transformation

library(vsn)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Colors palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)


# Remove outlier and do PCA :

df.keep <- df[,-which(colnames(df) %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34"))]
ids.keep <- ids[which(!(ids$Sample %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34"))),]
## keep only day 2:
ids.keep <- ids.keep[which(ids.keep$Group %in% c("Infected_d2","Control")),]
ids.keep$Group <- factor(ids.keep$Group)
df.keep <- df.keep[,ids.keep$Sample]
dds.keep <- DESeqDataSetFromMatrix(countData = df.keep, colData = ids.keep, design = ~ Group)
dds.keep <- dds.keep[ rowSums(counts(dds.keep)) > 10, ]
dds.keep <- estimateSizeFactors(dds.keep)
vsd.keep <- vst(dds.keep, blind = FALSE)
#Sample distances with VST:
sampleDists_vsd <- dist(t(assay(vsd.keep)))
#Transform sample distances to matrix:
sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
#Modify row names joining treatment and cell type:
rownames(sampleDistMatrix_vsd) <- paste0(vsd.keep$Group,":",vsd.keep$Sample)
#Remove column names:
colnames(sampleDistMatrix_vsd) <- NULL
# Draw heatmap
heatmap <- pheatmap(sampleDistMatrix_vsd,
                    clustering_distance_rows = sampleDists_vsd,
                    clustering_distance_cols = sampleDists_vsd,
                    col = colors)
pdf("plots.c52.feb2023/heatmap.QC.outliersRemoved2.Day2vsControl.pdf")
heatmap
dev.off()
pdf("plots.c52.feb2023/PCA.sampleID.outlierRemoved2.Day2vsControl.pdf",w=10,h=10)
plotPCA(vsd.keep, intgroup=c("Sample"))+geom_label(aes(label = Sample), position = position_nudge(y = 1))
dev.off()
pdf("plots.c52.feb2023/PCA.bygroup.outlierRemoved2.Day2vsControl.pdf")
plotPCA(vsd.keep, intgroup="Group")
dev.off()

#  2. Differential gene expression


# Run the DESeq function on the data:
dds_DGE <- DESeq(dds.keep)

## Building the results table


dds_DGE_results <- results(dds_DGE)
head(dds_DGE_results)


#The resulting object can be filtered like a data frame. 
table(dds_DGE_results$padj < 0.05)

summary(dds_DGE_results)

save(dds_DGE_results, file="DE.c52.4outliersremoved.RData")


#Select significant genes with a p adjusted lower than 0.05:
resSig <- subset(dds_DGE_results, padj < 0.05)

# Markers of interest:
allgenes <- c("Tnf","Fasl","Tnfsf8","Tnfsf4","Tnfsf9","Tnfsf10","Tnfsf14","Tnfsf13b",
             "Tnfsf15","Tnfsf18","Lta","Tnfsf11","Ltb","Tnfsf12","Tnfsf13","Eda",
             "Cd40lg","Tnfsf13os","Cd70")

resSig.selected <- dds_DGE_results[which(rownames(dds_DGE_results) %in% allgenes),]
resSig.selected <- resSig.selected[order(resSig.selected$padj),]
write.csv(resSig.selected, file="plots.c52.feb2023/DESeq2_Sig_results_InfectedDay2VsControl.MarkersOfInterest.csv")


## Volcano plot

library("dplyr")

# Add labels for significance:
results_order <- as.data.frame(dplyr::mutate(as.data.frame(resSig.selected), 
                                             sig=ifelse(resSig.selected$padj<0.05, 
                                                        ifelse(resSig.selected$log2FoldChange>1,
                                                               "FDR<0.05 & logFC>1", 
                                                               ifelse(resSig.selected$log2FoldChange<(-1),
                                                                      "FDR<0.05 & logFC<-1","Not Sig")),"Not Sig")), 
                               row.names=rownames(resSig.selected))
head(results_order)


library(ggrepel)

# Select genes that have log2 fold change higher than 2 and p adjusted lower than 0.05:
DEgenes_DESeq <- results_order[which(abs(results_order$log2FoldChange) > 1 & results_order$padj < 0.05),]
DEgenes_DESeq <- DEgenes_DESeq[order(DEgenes_DESeq$padj),]

# Add this annotation to the plot:
pdf("plots.c52.feb2023/volcanoPlot.InfectedDay2VsControl.markersOfInterest.pdf")
volcanoP + 
  ggrepel::geom_text_repel(data=data.frame(resSig.selected), 
                           aes(label=rownames(results_order)),
                           max.overlaps=100)+
  xlim(c(max(abs(min(DEgenes_DESeq$log2FoldChange)),max(DEgenes_DESeq$log2FoldChange))*(-1)-1,
         max(abs(min(DEgenes_DESeq$log2FoldChange)),max(DEgenes_DESeq$log2FoldChange))+1))

dev.off() 


## Heat maps

# define the colours
annoCol<-list(Day=c(d2="#E7D2E1", d3="#7E668C", d4="#8A0089", d5="#330055"),
              Condition=c(Control="#C5C5C8",Infected="#6A004A"))

annoCol2<-list(Condition=c(Control="#C5C5C8",Infected="#6A004A"))


#Select the markers of interest:
mat  <- assay(vsd.keep)[ allgenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd.keep)[, c("Condition","Day")])

#sorted genes:
resSig.selected <- resSig.selected[order(resSig.selected$log2FoldChange),]

pdf("plots.c52.feb2023/heatmap.InfectedDay2VsControl.markersOfInterest.vsdNormalised.pdf",w=5)
pheatmap(mat[rev(rownames(resSig.selected)),], cluster_rows = FALSE,
         annotation_col = anno, 
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

resSig.selected <- resSig.selected[order(resSig.selected$padj),]

pdf("plots.c52.feb2023/heatmap.InfectedDay2VsControl.markersOfInterest.vsdNormalised.sortedByPvalue.pdf",w=5)
pheatmap(mat[rownames(resSig.selected),], cluster_rows = FALSE,
         annotation_col = anno, 
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()


## Also get the full heat map with all conditions:

df.keep2 <- df[,-which(colnames(df) %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34"))]
ids.keep2 <- ids[which(!(ids$Sample %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34"))),]
dds.full <- DESeqDataSetFromMatrix(countData = df.keep2, colData = ids.keep2, design = ~ Group)
dds.full <- dds.full[ rowSums(counts(dds.full)) > 10, ]
dds.full <- estimateSizeFactors(dds.full)
vsd.full <- vst(dds.full, blind = FALSE)

#Select the markers of interest:
mat  <- assay(vsd.full)[ allgenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd.full)[, c("Condition","Day")])

anno.sort <- anno[order(anno$Day),]
anno.sort <- anno[order(anno.sort$Condition),]

pdf("plots.c52.feb2023/heatmap.InfectedAllDaysVsControl.markersOfInterest.vsdNormalised.sortedByDay.pdf",w=5,h=5)
pheatmap(mat[rev(rownames(resSig.selected)),rownames(anno.sort)], cluster_rows = FALSE,
         cluster_col = FALSE,
         annotation_col = anno, 
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

pdf("plots.c52.feb2023/heatmap.InfectedAllDaysVsControl.markersOfInterest.vsdNormalised.sortedByDayAndPvalue.pdf",w=5,h=5)
pheatmap(mat[rownames(resSig.selected),rownames(anno.sort)], cluster_rows = FALSE,
         cluster_col = FALSE,
         annotation_col = anno, 
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()


# also show normalised expr values:
## Read in mouse data:
gset.all <- read.delim("RNAseq Data C52/normalized_counts.tsv")

# Remove outliers:
gset.all <- gset.all[,which(!(colnames(gset.all) %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34")))]
gset <- gset.all[which(gset.all$gene %in% allgenes),]
rownames(gset) <- gset$gene
gset <- gset[,-c(1:2)]
gset <- log2(gset+1)


anno1 <- as.data.frame(colData(vsd.keep)[, c("Condition","Day")])
gset1  <- gset - rowMeans(gset)

pdf("plots.c52.feb2023/heatmap.InfectedDay2VsControl.markersOfInterest.normalisedValues.pdf",w=5)
pheatmap(gset[rev(rownames(resSig.selected)),rownames(anno)], cluster_rows = FALSE,
         annotation_col = anno1, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

pdf("plots.c52.feb2023/heatmap.InfectedDay2VsControl.markersOfInterest.normalisedValues.zscore.pdf",w=5)
pheatmap(gset1[rev(rownames(resSig.selected)),rownames(anno)], cluster_rows = FALSE,
         annotation_col = anno1, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

anno <- as.data.frame(colData(vsd.full)[, c("Condition","Day")])

anno.sort <- anno[order(anno$Day),]
anno.sort <- anno[order(anno.sort$Condition),]


pdf("plots.c52.feb2023/heatmap.InfectedAllDaysVsControl.markersOfInterest.normalisedValues.pdf",w=5,h=5)
pheatmap(gset[rev(rownames(resSig.selected)),rownames(anno.sort)], cluster_rows = FALSE,
         cluster_col = FALSE,
         annotation_col = anno, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

pdf("plots.c52.feb2023/heatmap.InfectedAllDaysVsControl.markersOfInterest.normalisedValues.zscore.pdf",w=5,h=5)
pheatmap(gset1[rev(rownames(resSig.selected)),rownames(anno.sort)], cluster_rows = FALSE,
         cluster_col = FALSE,
         annotation_col = anno, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()


#### Other lists of interest:
celldeath <- read.xls("gene list Maria-Anti, pro, IFNs.xlsx", sheet=1)
celldeath$Cell.death <- toupper(celldeath$Cell.death)
rownames(celldeath) <- celldeath$Cell.death
celldeath$Type <- 1
ifn <- read.xls("ifn_marie.xlsx", header=FALSE)$V1
ifnr <- read.xls("infreceptors_marie.xlsx", header=FALSE)$V1

gset.all <- read.delim("RNAseq Data C52/normalized_counts.tsv")
gset.all <- gset.all[,which(!(colnames(gset.all) %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34")))]
gset.celldeath <- gset.all[which(toupper(gset.all$gene) %in% celldeath$Cell.death),]
gset.ifn <- gset.all[which(toupper(gset.all$gene) %in% toupper(ifn)),]
gset.ifnr <- gset.all[which(toupper(gset.all$gene) %in% toupper(ifnr)),]

rownames(gset.celldeath) <- toupper(gset.celldeath$gene)
gset.celldeath <- gset.celldeath[,-c(1:2)]
gset.celldeath <- log2(gset.celldeath+1)
rownames(gset.ifn) <- gset.ifn$gene
gset.ifn <- gset.ifn[,-c(1:2)]
gset.ifn <- log2(gset.ifn+1)
rownames(gset.ifnr) <- gset.ifnr$gene
gset.ifnr <- gset.ifnr[,-c(1:2)]
gset.ifnr <- log2(gset.ifnr+1)

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFN.pdf",w=6,h=5)
pheatmap(gset.ifn[,rownames(anno.sort)], cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_col = anno, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFNreceptors.pdf",w=6,h=2)
pheatmap(gset.ifnr[,rownames(anno.sort)], cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_col = anno, 
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.CellDeath.pdf",w=6,h=12)
pheatmap(gset.celldeath[,rownames(anno.sort)], cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_col = anno, 
         annotation_row = celldeath[,c(2:3)],
         annotation_colors = annoCol,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100))
dev.off()


######################
### Now average by day and normalise to control:

gset.all <- read.delim("RNAseq Data C52/normalized_counts.tsv")
gset.all <- gset.all[,which(!(colnames(gset.all) %in% c("S5962Nr3","S5962Nr6","S5962Nr72","S5962Nr34")))]
gset.celldeath <- gset.all[which(toupper(gset.all$gene) %in% celldeath$Cell.death),]
gset.ifn <- gset.all[which(toupper(gset.all$gene) %in% toupper(ifn)),]
gset.ifnr <- gset.all[which(toupper(gset.all$gene) %in% toupper(ifnr)),]
gset.markersOfInterest <- gset.all[which(gset.all$gene %in% allgenes),]

rownames(gset.celldeath) <- toupper(gset.celldeath$gene)
gset.celldeath <- gset.celldeath[,-c(1:2)]
rownames(gset.ifn) <- gset.ifn$gene
gset.ifn <- gset.ifn[,-c(1:2)]
rownames(gset.ifnr) <- gset.ifnr$gene
gset.ifnr <- gset.ifnr

rownames(gset.markersOfInterest) <- gset.markersOfInterest$gene
gset.markersOfInterest <- gset.markersOfInterest[,-c(1:2)]

anno <- data.frame(anno)
ctrl <- rownames(anno[which(anno$Condition=="Control"),])
d2 <- rownames(anno[which((anno$Condition=="Infected")&(anno$Day=="d2")),])
d3 <- rownames(anno[which((anno$Condition=="Infected")&(anno$Day=="d3")),])
d4 <- rownames(anno[which((anno$Condition=="Infected")&(anno$Day=="d4")),])
d5 <- rownames(anno[which((anno$Condition=="Infected")&(anno$Day=="d5")),])
ct <- rowMeans(gset.ifn[,ctrl])
ct[ct == 0] <- 0.0001
gset.ifn.avg <- cbind(ct/ct,
                      rowMeans(gset.ifn[,d2])/ct,
                      rowMeans(gset.ifn[,d3])/ct,
                      rowMeans(gset.ifn[,d4])/ct,
                      rowMeans(gset.ifn[,d5])/ct)
colnames(gset.ifn.avg) <- c("Control","D2","D3","D4","D5") 
gset.ifnr.avg <- cbind(rowMeans(gset.ifnr[,ctrl])/rowMeans(gset.ifnr[,ctrl]),
                      rowMeans(gset.ifnr[,d2])/rowMeans(gset.ifnr[,ctrl]),
                      rowMeans(gset.ifnr[,d3])/rowMeans(gset.ifnr[,ctrl]),
                      rowMeans(gset.ifnr[,d4])/rowMeans(gset.ifnr[,ctrl]),
                      rowMeans(gset.ifnr[,d5])/rowMeans(gset.ifnr[,ctrl]))
colnames(gset.ifnr.avg) <- c("Control","D2","D3","D4","D5") 
gset.celldeath.avg <- cbind(rowMeans(gset.celldeath[,ctrl])/rowMeans(gset.celldeath[,ctrl]),
                      rowMeans(gset.celldeath[,d2])/rowMeans(gset.celldeath[,ctrl]),
                      rowMeans(gset.celldeath[,d3])/rowMeans(gset.celldeath[,ctrl]),
                      rowMeans(gset.celldeath[,d4])/rowMeans(gset.celldeath[,ctrl]),
                      rowMeans(gset.celldeath[,d5])/rowMeans(gset.celldeath[,ctrl]))
colnames(gset.celldeath.avg) <- c("Control","D2","D3","D4","D5")


gset.markersOfInterest.avg <- cbind(rowMeans(gset.markersOfInterest[,ctrl])/rowMeans(gset.markersOfInterest[,ctrl]),
                            rowMeans(gset.markersOfInterest[,d2])/rowMeans(gset.markersOfInterest[,ctrl]),
                            rowMeans(gset.markersOfInterest[,d3])/rowMeans(gset.markersOfInterest[,ctrl]),
                            rowMeans(gset.markersOfInterest[,d4])/rowMeans(gset.markersOfInterest[,ctrl]),
                            rowMeans(gset.markersOfInterest[,d5])/rowMeans(gset.markersOfInterest[,ctrl]))
colnames(gset.markersOfInterest.avg) <- c("Control","D2","D3","D4","D5") 


pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFN.foldChange.pdf",w=3,h=4)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.ifn.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.ifn.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.ifn.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFN.Log2foldChange.pdf",w=3,h=4)
gset.ifn.avg <- log2(gset.ifn.avg+1)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.ifn.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.ifn.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.ifn.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFN.Log2foldChange.without6_7_11.pdf",w=3,h=4)
gset.ifn.avg <- gset.ifn.avg[which(!(rownames(gset.ifn.avg) %in% c("Ifna7","Ifna11","Ifna6"))),]
gset.ifn.avg <- log2(gset.ifn.avg+1)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.ifn.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.ifn.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.ifn.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.IFNreceptors.foldChange.pdf",w=4,h=3)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.ifnr.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.ifnr.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.ifnr.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.CellDeath.foldChangeUncapped.pdf",w=4,h=12)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.celldeath.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.celldeath.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.celldeath.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_row = celldeath[,c(2:3)],
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.CellDeath.foldChangeCappedAt4.pdf",w=4,h=12)
gset.celldeath.avg[gset.celldeath.avg>4] <- 4 
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.celldeath.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.celldeath.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.celldeath.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_row = celldeath[,c(2:3)],
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.CellDeath.Log2foldChange.pdf",w=4,h=12)
gset.celldeath.avg <- log2(gset.celldeath.avg+1) 
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.celldeath.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.celldeath.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.celldeath.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         annotation_row = celldeath[,c(2:3)],
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()


pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.TNFMarkersOfInterest.foldChange.pdf",w=3,h=4)
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.markersOfInterest.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.markersOfInterest.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.markersOfInterest.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

pdf("plots.c52.4outliersRemoved.feb2023/heatmapsByDayAndCondition/heatmap.InfectedAllDaysVsControl.normalisedValues.TNFMarkersOfInterest.Log2foldChange.pdf",w=3,h=4)
gset.markersOfInterest.avg <- log2(gset.markersOfInterest.avg+1) 
paletteLength <- 100
myColor <- colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength)
myBreaks <- c(seq(min(gset.markersOfInterest.avg), 1, length.out=ceiling(paletteLength/2) ), 
              seq(1.01, max(gset.markersOfInterest.avg), length.out=floor(paletteLength/2)))
pheatmap(gset.markersOfInterest.avg, cluster_rows = TRUE,
         cluster_col = FALSE,
         color = colorRampPalette(c("#17255A","white","#BD1E1E"))(paletteLength),
         breaks = myBreaks)
dev.off()

## 4. Annotate results

library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)


# Add functional path column:
resSig$PATH = mapIds(org.Mm.eg.db,
                     keys=rownames(resSig), 
                     column="PATH",
                     keytype="SYMBOL",
                     multiVals="first")

# Add Entrez ID column:
resSig$ENTREZID = mapIds(org.Mm.eg.db,
                         keys=rownames(resSig), 
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")

# Add ENSEMBL column:
resSig$ENSEMBL <- rownames(resSig)

resSig_DF <- as.data.frame(resSig)


## 4. Functional Gene Ontology (GO) analysis

library("clusterProfiler")

OrgDb <- org.Mm.eg.db # this is the mouse database, but you can also use other organisms

# Get ENTREZID as vector:
genes <- as.character(resSig$ENTREZID)

genes.upreg <- unique(as.character(resSig[which(resSig$log2FoldChange>0),]$ENTREZID))
genes.downreg <- unique(as.character(resSig[which(resSig$log2FoldChange<0),]$ENTREZID))

save(resSig, file="resSig.c52.RData")

mydf <- data.frame(Entrez=resSig$ENTREZID, FC=resSig)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

resSig$Group <- sapply(resSig$log2FoldChange, function(x) ifelse(x>0,"upregulated","downregulated"))

resSig.filter <- resSig[which(!is.na(resSig$ENTREZID)),]
formula_res <- compareCluster(ENTREZID~Group, 
                              data=data.frame(resSig.filter), 
                              fun="enrichKEGG",
                              OrgDb='org.Mm.eg.db')

pdf("plots.c52.feb2023/clusterProfiler.GOenrichment.pdf")
dotplot(formula_res)
dev.off()

pdf("plots.c52.feb2023/clusterProfiler.GOnetwork.pdf")
cnetplot(formula_res)
dev.off()

# Next, load mouse to human mapping:
load("../covid/mappingMusHu.GSE33266.RData")

resSig.filter.hu <- merge(data.frame(resSig.filter),
                          mappingMusHu[,c("MGI.symbol","HGNC.symbol")],
                          by.x="ENSEMBL", by.y="MGI.symbol",
                          all.x=FALSE, all.y=FALSE)
resSig.filter.hu$ENTREZID = mapIds(org.Hs.eg.db,
                         keys=resSig.filter.hu$HGNC.symbol, 
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")

formula_res_hu <- compareCluster(ENTREZID~Group, 
                              data=data.frame(resSig.filter.hu), 
                              fun="enrichKEGG")

pdf("plots.c52.feb2023/clusterProfiler.KEGGenrichment.pdf")
dotplot(formula_res_hu)
dev.off()

pdf("plots.c52.feb2023/clusterProfiler.KEGGnetwork.pdf")
cnetplot(formula_res_hu)
dev.off()



library(pathfindR)
p <- "KEGG"
  pathname <- paste0("pathfindR_Results_c52_InfectedDay2VsControl")
  #dir.create(file.path(pathname), showWarnings = FALSE)
  df <- data.frame(resSig_DF[which(abs(resSig_DF$log2FoldChange)>2),c("ENSEMBL","log2FoldChange","padj")])
  output_df <- run_pathfindR(df,output_dir=pathname,
                             gene_sets=p) 
  od <- data.frame(output_df)
  write.csv(od, file=paste0(pathname,"/enrichment.csv"))
  pdf(paste0(pathname,"/topEnrichedPathways.InfectedVsControl.pdf"),w=15)
  print(enrichment_chart(result_df = output_df, top_terms = 20) )
  dev.off()
  
  pdf(paste0(pathname,"/geneGraph.InfectedVsControl.pdf"),w=15,h=15)
  print(term_gene_graph(RA_output))
  dev.off()

  od.keep <- od[which(od$Term_Description %in%
                        c("Necroptosis","Cytokine-cytokine receptor interaction",
                          "Viral protein interaction with cytokine and cytokine receptor",
                          "TNF signaling pathway","Antigen processing and presentation",
                          "Osteoclast differentiation","Proteasome",
                          "Coronavirus disease - COVID-19",
                          "JAK-STAT signaling pathway",
                          "NOD-like receptor signaling pathway",
                          "PD-L1 expression and PD-1 checkpoint pathway in cancer",
                          "MAPK signaling pathway",
                          "NF-kappa B signaling pathway",
                          "Prolactin signaling pathway",
                          "RIG-I-like receptor signaling pathway",
                          "Complement and coagulation cascades",
                          "C-type lectin receptor signaling pathway",
                          "Phagosome",
                          "Apoptosis",
                          "Growth hormone synthesis, secretion and action",
                          "IL-17 signaling pathway",
                          "Cytosolic DNA-sensing pathway",
                          "Estrogen signaling pathway",
                          "ECM-receptor interaction",
                          "Ubiquitin mediated proteolysis",
                          "Neutrophil extracellular trap formation",
                          "Th17 cell differentiation",
                          "Th1 and Th2 cell differentiation",
                          "Cell cycle",
                          "Cellular senescence",
                          "TGF-beta signaling pathway",
                          "p53 signaling pathway",
                          "Adipocytokine signaling pathway",
                          "Natural killer cell mediated cytotoxicity",
                          "ErbB signaling pathway",
                          "B cell receptor signaling pathway",
                          "Endocytosis",
                          "Hippo signaling pathway",
                          "HIF-1 signaling pathway")),]

 df.p <-  od.keep[rev(order(od.keep$Fold_Enrichment)),]
df.p$Term_Description <- factor(df.p$Term_Description, levels=rev(df.p$Term_Description))

pdf(paste0(pathname,"/pathwayEnrichment.pdf"))
ggplot(df.p, aes(x = Term_Description, y = Fold_Enrichment)) +
  geom_segment(aes(x = Term_Description, xend = Term_Description, y = 0, yend = Fold_Enrichment),
               color = "gray", lwd = 1.5) +
  geom_point(size = df.p$support*50, pch = 21, bg = 4, col = 1) +
  geom_text(aes(label = paste0(round(df.p$support*100,1),"%")), 
            color = "#13315C", size = 2.5,
            hjust = -0.6) +
  coord_flip()+
  geom_hline(yintercept=1, linetype="dashed", 
             color = "black", size=0.2)+
  ylim(c(0,max(df.p$Fold_Enrichment)+1))+
  ylab("Fold enrichment")+
  xlab("")
dev.off() 

### finally, make heat maps for selected pathways:

pathways <- c("TNF signaling pathway","Necroptosis","Apoptosis")

for (pt in pathways) {
  
  p1 <- df.p[which(df.p$Term_Description == pt),]
  p1.up <- strsplit(p1$Up_regulated,", ")[[1]]
  p1.down <- strsplit(p1$Down_regulated,", ")[[1]]
  p1.genes <- c(p1.up, p1.down)
  
  wd <- length(p1.genes)/3+1
   
  gset.p1 <- gset.all[which(toupper(gset.all$gene) %in% p1.genes),]
  rownames(gset.p1) <- gset.p1$gene
  gset.p1 <- gset.p1[,-c(1:2)]
  gset.p1 <- log2(gset.p1+1)
  
  pdf(paste0(pathname,"/heatmap.InfectedDay2vsControl.",pt,".pdf"),w=5,h=wd)
  print(pheatmap(gset.p1[,rownames(anno1)], cluster_rows = FALSE,
           annotation_col = anno1, 
           annotation_colors = annoCol,
           color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100)))
  dev.off() 
  
  pdf(paste0(pathname,"/heatmap.AllDaysvsControl.",pt,".pdf"),w=5,h=wd)
  print(pheatmap(gset.p1[,rownames(anno.sort)], cluster_rows = FALSE,
                 cluster_cols = FALSE,
           annotation_col = anno, 
           annotation_colors = annoCol,
           color = colorRampPalette(c("#17255A","white","#BD1E1E"))(100)))
  dev.off() 

}

