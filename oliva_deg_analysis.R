##analysis of Al and bacteria sample of rice
##no biological replicates, using DESeq package
inp_micro_fe <- inp_micro[, -c(4:5)]
##define condition
condition_fe <- meta_micro[1:3, ]
## create cds
cds_fe = newCountDataSet(inp_micro_fe, condition_fe)
## estimate size factor
cds_fe = estimateSizeFactors( cds_fe )
##dispersion estimate
cds_fe = estimateDispersions( cds_fe, method="blind", sharingMode="fit-only" )
##differential analysis
res_fenb = nbinomTest( cds_fe, "cont", "fenb" )
res_fenb = nbinomTest( cds_fe, "cont", "fewb" )
##export result
write.csv(res_fenb, "fenb.csv")
write.csv(res_fewb, "fewb.csv")
vsdFull = varianceStabilizingTransformation( cds_fe)
select = order(rowMeans(counts(cds_fe)), decreasing=TRUE)[1:30]
pheatmap(exprs(vsdFull)[select,])
dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cds_fe), paste(condition_fe))
pheatmap(mat)
print(plotPCA(vsdFull))
