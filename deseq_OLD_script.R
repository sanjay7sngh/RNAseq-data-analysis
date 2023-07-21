## This is the code of old deseq packege
## here samples without replicated can be analyzed
##analysis of Al and bacteria sample of rice
##no biological replicates, using DESeq package

# load required library
library(readxl)
# library(DESeq)
library(edgeR)

# Deseq is not working for recent version of R, so using edgeR
# import data
count_data <- read_excel("C:/Users/IASRI/Dropbox/IASRI/osa/count.xlsx", sheet = "inp")

metadata <- read_excel("C:/Users/IASRI/Dropbox/IASRI/osa/count.xlsx", sheet = "metadata")

# clean count data in the format of deseq 
inp_micro <-  as.data.frame(count_data[, 7:11])
rownames(inp_micro) <- count_data$Geneid

# clean meta data in the format of deseq 
meta_micro <- as.data.frame(metadata[, 2:3])
rownames(meta_micro) <- metadata$sample
# check everything is formatted
head(meta_micro)
head(inp_micro)
## for aluminium
inp_micro_al <- inp_micro[, -c(4:5)]
##define condition
condition_al <- meta_micro[-c(4:5), ]

# al check
head(inp_micro_al)
head(condition_al)
##############################################################################
# edgeR specific code for DEG
# yedge <- DGEList(counts = inp_micro_al)
# head(yedge$samples)
# group <- condition_al$condition
# yedge <- DGEList(counts = inp_micro_al, group = group)
# yedge$samples
# yedge <- calcNormFactors(yedge)
# plotMDS(yedge)
# yedge <- estimateDisp(yedge)
# bcv <- 0.2
# 
# alnb_res <- exactTest(yedge, pair = c("control", "Al_NB"), dispersion=bcv^2)
# head(alnb_res)
##############################################################################
# create cds
# cds <- newCountDataSet( countData = data, conditions = as.factor(c("treated","untreated")))
cds_al = newCountDataSet(inp_micro_al, conditions = as.factor(c("control","Al_NB","Al_WB")))
# estimate size factor
cds_al = estimateSizeFactors(cds_al)
#dispersion estimate
cds_al = estimateDispersions( cds_al,  method='blind', sharingMode="fit-only")
#differential analysis
res_alwb = nbinomTest(cds_al, "control", "Al_WB")
res_alnb = nbinomTest(cds_al, "control", "Al_NB")
#export result
write.csv(res_alnb, "alnb.csv")
write.csv(res_alwb, "alwb.csv")
vsdFull = varianceStabilizingTransformation( cds_al)
select = order(rowMeans(counts(cds_al)), decreasing=TRUE)[1:30]
pheatmap(exprs(vsdFull)[select,])
dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cds_al), paste(condition))
pheatmap(mat)
print(plotPCA(vsdFull))
# 