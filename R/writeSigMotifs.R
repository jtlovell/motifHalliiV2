writeSigMotifs<-function(calcMotifEnrichmentOutput,csv.file="motifs.csv", threshold=0.05){
  stats<-calcMotifEnrichmentOutput$enrichmentStats
  write.csv(stats[stats$fdr.pvalue<=threshold,], file=csv.file, row.names=F)
}
