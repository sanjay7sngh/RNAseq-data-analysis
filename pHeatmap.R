## general pheatmap script
library(readxl)
osa_photosyn <- read_excel("osa_DEGs_photosynthesis.xlsx", sheet = "input")
library(pheatmap)
photosyn_heat <- as.data.frame(osa_photosyn[, 2:9])
rownames(photosyn_heat) <- osa_photosyn$Gene_symbol
class(photosyn_heat) ##make sure it is a matrix or dataframe
## column annotation
##ann_col <- data.frame(ann_col = c("G12", "G24", "G48", "G72", "S12", "S24", "S48", "S72"))
##rownames(ann_col) <- names(photosyn_heat)
##or create in excel and import
library(readxl)
ann_col <- read_excel("osa_DEGs_photosynthesis.xlsx", sheet = "ann_col")
head(ann_col)
ann_col_m <- as.data.frame(ann_col[,2:3])
row.names(ann_col_m) <- ann_col$sample
ann_col_m$HPI <- as.factor(ann_col_m$HPI)
head(ann_col_m)
##row annotation
ann_row <- as.data.frame(osa_photosyn[, 14])
row.names(ann_row) <- osa_photosyn$Gene_symbol
## Define color and break--- https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
paletteLength <- 50
myColor <- colorRampPalette(c("magenta", "white", "blue"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(photosyn_heat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(photosyn_heat)/paletteLength, max(photosyn_heat), length.out=floor(paletteLength/2)))

# Plot the heatmap
##pheatmap(test, color=myColor, breaks=myBreaks)
##colors = c(seq(-8.0,-1.0,length=70),seq(-0.99,0.99,length=20),seq(1.0,5.0,length=40))
##colors <- c(-8, -1, 1, 5)
##my_palette <- colorRampPalette(c("blue","green","black","magenta"))(n = 100)
library(RColorBrewer)
pheatmap(photosyn_heat, cellwidth = 12, annotation_col = ann_col_m, color=myColor, breaks=myBreaks)
pheatmap(photosyn_heat, cellwidth = 12, annotation_col = ann_col_m, annotation_row = ann_row, color=myColor, breaks=myBreaks)



##colors = c(seq(0,0.0001,length=100),seq(0.0099,0.05,length=100),seq(0.051,0.15,length=100), seq(0.159,1,length=100))
##pheatmap(enrFunHeat, cluster_cols = FALSE, annotation_row = enrFunAnn, color = my_palette, breaks = colors)
##my_palette <- colorRampPalette(c("blue","magenta", "black", "green", "red"))(n = 299)
##pheatmap(enrFunHeat, cluster_cols = FALSE, annotation_row = enrFunAnn, color = my_palette, breaks = colors)
##pheatmap(enrFunHeat, cluster_cols = FALSE, annotation_row = enrFunAnn, color = my_palette, breaks = colors,  symm=F,symkey=F,symbreaks=T, scale="none")
##my_palette <- colorRampPalette(c("blue","magenta", "black", "green"))(n = 299)
##pheatmap(enrFunHeat, cluster_cols = FALSE, annotation_row = enrFunAnn, color = my_palette, breaks = colors,  symm=F,symkey=F,symbreaks=T, scale="none")
##pheatmap(enrFunHeat, cluster_cols = FALSE, annotation_row = enrFunAnn, color = my_palette, breaks = colors,  symm=F,symkey=F,symbreaks=F, scale="none")