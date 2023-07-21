##analysis of Al and bacteria sample of rice
##no biological replicates, using DESeq package
inp_micro_al <- inp_micro[, -c(2:3)]
##define condition
condition_al <- meta_micro[-c(2:3), ]
## create cds
cds_al = newCountDataSet(inp_micro_al, condition_al)
## estimate size factor
cds_al = estimateSizeFactors( cds_al )
##dispersion estimate
cds_al = estimateDispersions( cds_al, method="blind", sharingMode="fit-only" )
##differential analysis
res_alnb = nbinomTest( cds_al, "cont", "alnb" )
res_alnb = nbinomTest( cds_al, "cont", "alwb" )
##export result
write.csv(res_alnb, "alnb.csv")
write.csv(res_alwb, "alwb.csv")
vsdFull = varianceStabilizingTransformation( cds_al)
select = order(rowMeans(counts(cds_al)), decreasing=TRUE)[1:30]
pheatmap(exprs(vsdFull)[select,])
dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cds), paste(condition_al))
pheatmap(mat)
print(plotPCA(vsdFull))
