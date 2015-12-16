\name{plotSigMotif}
\alias{plotSigMotif}
\title{Plotting function for motifHalliiV2}
\description{Writes a plot of significant motifs data (including logo) to file}
\usage{
plotSigMotif(calcMotifEnrichmentOutput, pdf.file = "motifs.pdf", threshold = 0.05)
}

\arguments{
  \item{calcMotifEnrichmentOutput}{output from calcMotifEnrichment
}
  \item{pdf.file}{the name (possibly including the directory) of the output file)
}
  \item{threshold}{the threshold to plot - uses mutliple testing corrected pvalue from calcMotifEnrichmentOutput
}
}

\value{A .pdf file is writen. Nothing is returned
}

\author{
J.T. Lovell
}


\examples{
library(motifHalliiV2)
library(devtools)

v2geneID<-c("Pahal.J02145", "Pahal.J01516","Pahal.H01249","Pahal.D01315","Pahal.J00094","Pahal.D00330","Pahal.D02670",
"Pahal.D01649", "Pahal.D02828", "Pahal.J00146", "Pahal.C04014", "Pahal.H01253", "Pahal.D00532","Pahal.J00293",
"Pahal.F02021", "Pahal.F02121", "Pahal.F01847", "Pahal.F01958", "Pahal.F00534", "Pahal.B02910")

#bonferroni multiple correction
test<-calcMotifEnrichment(v2geneID[1:10], FDR.method="p.adjust", method="bonferroni")

plotSigMotif(test, pdf.file="~/Desktop/test.pdf", threshold=0.5)


\keyword{ ~Motif }
\keyword{ ~Enrichment }
