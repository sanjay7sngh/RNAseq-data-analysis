
## Differential expression analysis using ballgown
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

## Read phenotype sample data

pheno_data <- read.csv("geuvadis_phenodata.csv")

## Read in expression data
bg_cotton <- ballgown(dataDir = "ballgown/count_cotton/", samplePattern="SRR")

## Filter low abundance genes
bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)

## DE by transcript
results_transcripts <-  stattest(bg_chrX_filt, feature='transcript', covariate='sex', 
                                 adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## DE by gene
results_genes <-  stattest(bg_chrX_filt, feature='gene', covariate='sex', 
                           adjustvars=c('population'), getFC=TRUE, meas='FPKM')

## Add gene name
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),
                                  geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, "chrX_transcripts_results.csv", row.names=FALSE)
write.csv(results_genes, "chrX_genes_results.csv", row.names=FALSE)

## Filter for genes with q-val <0.05
results_transcripts_sub <- subset(results_transcripts, results_transcripts$qval <=0.05)
results_genes <- subset(results_genes, results_genes$qval <=0.05)

## Plotting setup
tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

## Plotting gene abundance distribution
fpkm <- texpr(bg_chrX, meas='FPKM')
fpkm <- log2(fpkm +1)
boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')


## Plot gene of transcript 1729
plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))


