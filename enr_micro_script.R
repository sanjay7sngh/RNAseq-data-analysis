## enrichment with clusterProfiler of microbiology data
# creating vector of significant gene 
alnb_v <- alnb$GeneID
alwb_v <- alwb$GeneID
fenb_v <- fenb$GeneID
fewb_v <- fewb$GeneID
##enrichment with clusterProfiler
alnb_enr <- enricher(alnb_v, TERM2GENE = rice_go, pvalueCutoff = 1)
alwb_enr <- enricher(alwb_v, TERM2GENE = rice_go, pvalueCutoff = 1)
fenb_enr <- enricher(fenb_v, TERM2GENE = rice_go, pvalueCutoff = 1)
fewb_enr <- enricher(fewb_v, TERM2GENE = rice_go, pvalueCutoff = 1)

#join enr_data to GO annotation
alnb_enr <- as.data.frame(alnb_enr) %>% left_join(., rice_goAnno, by = "ID")
alwb_enr <- as.data.frame(alwb_enr) %>% left_join(., rice_goAnno, by = "ID")
fenb_enr <- as.data.frame(fenb_enr) %>% left_join(., rice_goAnno, by = "ID")
fewb_enr <- as.data.frame(fewb_enr) %>% left_join(., rice_goAnno, by = "ID")
# export result as csv file
write.csv(alnb_enr, file = "alnb_enr.csv")
write.csv(alwb_enr, file = "alwb_enr.csv")
write.csv(fenb_enr, file = "fenb_enr.csv")
write.csv(fewb_enr, file = "fewb_enr.csv")