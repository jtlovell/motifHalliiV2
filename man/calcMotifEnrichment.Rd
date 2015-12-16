\name{calcMotifEnrichment}
\alias{calcMotifEnrichment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simple inference of motif enrichment for P. hallii}
\description{Calculates Motif Enrichment for a vector of P. hallii V2 gene IDs}
\usage{calcMotifEnrichment(v2geneID, dropDuplicates = TRUE, verbose = TRUE, FDR.method = "qvalue", ...)}

\arguments{
  \item{v2geneID}{A Vector of V2 gene ID e.g. Pahal.J02145}
  \item{dropDuplicates}{Logical, If there are duplicate geneIDs in your list, should they be dropped?}
  \item{verbose}{Logical, should stepwise progress be reported?}
  \item{FDR.method}{If "qvalue", use the qvalue package, any other specification uses p.adjust}
  \item{\dots}{If FDR.method != "qvalue", specify the parameters of p.adjust}
}
\details{
%%  ~~ Nothing yet ~~
}
\value{
A two element list with the motif enrichment output
\item{enrichmentStats}{A dataframe of statistics and annotations for the motifs}
\item{enrichmentResults}{A PWMEnrich object of enrichment results, can be used for plotting}
}
\references{
%% ~ Nothing yet ~
}
\author{
%%  ~~J.T. Lovell~~
}
\note{
%%  ~~Nothing Yet~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(motifHalliiV2)
library(devtools)

v2geneID<-c("Pahal.J02145", "Pahal.J01516","Pahal.H01249","Pahal.D01315","Pahal.J00094","Pahal.D00330","Pahal.D02670",
"Pahal.D01649", "Pahal.D02828", "Pahal.J00146", "Pahal.C04014", "Pahal.H01253", "Pahal.D00532","Pahal.J00293",
"Pahal.F02021", "Pahal.F02121", "Pahal.F01847", "Pahal.F01958", "Pahal.F00534", "Pahal.B02910")

#no multiple correction
test<-calcMotifEnrichment(v2geneID[1:10], FDR.method="p.adjust", method="none")

#qvalue multiple correction
test<-calcMotifEnrichment(v2geneID[1:10], FDR.method="qvalue")

#bonferroni multiple correction
test<-calcMotifEnrichment(v2geneID[1:10], FDR.method="p.adjust", method="bonferroni")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~motif }
\keyword{ ~enrichment }
