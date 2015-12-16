calcMotifEnrichment<-function(v2geneID, dropDuplicates=TRUE, verbose=TRUE, FDR.method="qvalue",...){
  require(Biostrings)
  require(PWMEnrich)
  require(qvalue)

  if(verbose) cat("loading promoter and motif information","\n")
  data(all.info)
  data(pwms)
  data(promoters)
  data(bg)

  if(verbose) cat("grabbing genes from database of", length(promoters), "promoter sequences\n", sep=" ")
  seq<-promoters[v2geneID]
  if(dropDuplicates){
    tokeep<-which(!duplicated(names(seq)))
    seq<-seq[tokeep]
  }

  if(verbose) cat("running motif enrichment analysis using ", length(seq), "promoter sequences\n", sep=" ")
  res = motifEnrichment(seq, score="affinity", pwms=bg, group.only=T, bg="logn")
  r = groupReport(res)
  report<-as.data.frame(r@d)
  report<-report[report$p.value>=0 & report$p.value<=1,]
  if(FDR.method=="qvalue"){
    report$fdr.pvalue<-qvalue(report$p.value)$qvalues
  }else{
    report$fdr.pvalue<-p.adjust(report$p.value, ...)
  }
  colnames(report)[2]<-"ID"

  enrich.output<-merge(report, all.info, by="ID", all.x=T)
  enrich.output<-enrich.output[order(enrich.output$fdr.pvalue),]
  return(list(enrichmentStats=enrich.output, enrichmentResults=r))
}
