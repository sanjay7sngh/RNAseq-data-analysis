##create plots heatmap using fpkm value
##Set colour and heatmap scaling breaks
## https://www.biostars.org/p/284721/
## Use this code (below).

## Your FPKM values will be stored in MyFPKMValues
## DiffExpressedGenes will comprise a single vector of genes that are differentially expressed
## zFPKM package will be used to convert your FPKM values to the Z-scale prior to clustering.
require("RColorBrewer")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)
#############import data######################
library(readxl)
ballgown_cotton_fpkm <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/ballgown_cotton_fpkm.xlsx", sheet = "inp")
                                                                         

cotton_goi <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/goi2.xlsx", sheet = "inp")

###
MyFPKMValues <- as.data.frame(ballgown_cotton_fpkm[, 2:36])
row.names(MyFPKMValues) <- ballgown_cotton_fpkm$gene_id
###
DiffExpressedGenes <- cotton_goi$symbol

##Scale the FPKM values to the Z scale

library(zFPKM)
heat <- zFPKM(MyFPKMValues)
activeGenes <- which(rowMeans(heat) > -3)
## Filter your dataset to include your differentially expressed genes:
  
heat <- heat[which(rownames(heat) %in% DiffExpressedGenes), ]
write.csv(heat, file = "fpkm_goi_data.csv")
###############################################################################
## general pheatmap script
library(readxl)
cotton_heat <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/final_dataset/metadata", sheet = "final_data")
library(pheatmap)
##cotton_heat2 <- dplyr::left_join(cotton_heat, cotton_goi, by = "gene_id")
##write.csv(cotton_heat2, file = "cotton_heat_data_joined.csv")
#############################################

cotton_meta <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/goi2.xlsx", sheet = "metadata")
##cotton_heat <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/goi2.xlsx", sheet = "inp2")
##############################################
cotton_heat2 <- as.data.frame(cotton_data[, 2:36])
rownames(cotton_heat2) <- cotton_data$symbol
class(cotton_heat2) ##make sure it is a matrix or dataframe
## column annotation
##ann_col <- data.frame(ann_col = c("G12", "G24", "G48", "G72", "S12", "S24", "S48", "S72"))
##rownames(ann_col) <- names(cotton_heat2)
##or create in excel and import
## this is modification to make heatmap elegant
to_ann <- as.data.frame(cotton_meta[,2:3]) 
rownames(to_ann) <- cotton_meta$symbol
##to_ann2 <- dplyr::left_join(gbs_heat_data, to_ann, by = "gene_id")


##ann_row <- as.data.frame(to_ann2[,c(4,5)])
##row.names(ann_row) <- to_ann2$gene_id
### ann_col_m$HPI <- as.factor(ann_col_m$HPI)
head(ann_row)
##row annotation
cotton_meta <- read_excel("C:/Users/IASRI/Dropbox/IASRI/cotton_117/cotton_gbs/rnaseq/goi2.xlsx", sheet = "label")
ann_col <- as.data.frame(cot_label[,3])
row.names(ann_col) <- cot_label$sample_name
## Define color and break--- https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
paletteLength <- 50
myColor <- colorRampPalette(c("magenta", "white", "blue"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heat40), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heat40)/paletteLength, max(heat40), length.out=floor(paletteLength/2)))

# Plot the heatmap
##pheatmap(test, color=myColor, breaks=myBreaks)
##colors = c(seq(-8.0,-1.0,length=70),seq(-0.99,0.99,length=20),seq(1.0,5.0,length=40))
##colors <- c(-8, -1, 1, 5)
##my_palette <- colorRampPalette(c("blue","green","black","magenta"))(n = 100)

### to remove rownames use --- show_rownames = F ----
library(RColorBrewer)
library(pheatmap)
pheatmap::pheatmap(cotton_heat2, cellwidth = 12, annotation_col = ann_col, color=myColor, breaks=myBreaks)
pheatmap::pheatmap(heat40, cellwidth = 12, annotation_row = to_ann, annotation_col = ann_col, color=myColor, breaks=myBreaks, show_rownames = F)

## to export cluster data folloew this link --- https://www.biostars.org/p/287512/

#####################################################################
## Generate heatmaps with Euclidean distance (first) and '1 - Pearson correlation' distance (second) (both use Ward's linkage)
###################################################################
# require("gplots")

#Euclidean distance
# heatmap.2(heat,
  col=myCol,
  breaks=myBreaks,
  main="Title",
  key=T, keysize=1.0,
  scale="none",
  density.info="none",
  reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
  trace="none",
  cexRow=0.2,
  cexCol=0.8,
  distfun=function(x) dist(x, method="euclidean"),
  hclustfun=function(x) hclust(x, method="ward.D2"))

#1-cor distance
# heatmap.2(heat,
  col=myCol,
  breaks=myBreaks,
  main="Title",
  key=T, keysize=1.0,
  scale="none",
  density.info="none",
  reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
  trace="none",
  cexRow=0.2,
  cexCol=0.8,
  distfun=function(x) as.dist(1-cor(t(x))),
  hclustfun=function(x) hclust(x, method="ward.D2"))