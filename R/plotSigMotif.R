plotSigMotif<-function(calcMotifEnrichmentOutput,pdf.file="motifs.pdf", threshold=0.05){
  stats<-calcMotifEnrichmentOutput$enrichmentStats
  results<-calcMotifEnrichmentOutput$enrichmentResults
  n.sig<-sum(stats$fdr.pvalue<=threshold)
  pdf(pdf.file, width=20, height=n.sig)
  plot(results[1:n.sig])
  dev.off()
}
