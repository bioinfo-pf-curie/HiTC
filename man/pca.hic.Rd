\name{pca.hic}
\Rdversion{1.1}
\alias{pca.hic}
\title{Perform Principle Component Analysis on Hi-C contact map}
\description{
  Perform Principle Component Analysis on Hi-C contact map
}
\usage{pca.hic(x, normPerExpected=TRUE, npc=2, asGRangesList=TRUE,
gene.gr=NULL, ...)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{normPerExpected}{normalized by expected interaction using the loess
    calculation of distance dependency. see normPerExpected}
  \item{npc}{numeric; number of first principal component to return}
  \item{asGRangesList}{if TRUE a GRangesList object is returned where the scores represent the eigenvector}
  \item{gene.gr}{object of class \code{GenomicRanges} describing the
    genes position. If used, the A/B compartments classes are defined
    based on gene density}
  \item{...}{additional parameters passed to \code{normPerExpected} function}
}

\details{
 This method was apply by Lieberman-Aiden et al. 2009 to correlate the
 annotation profiles (genes, ChIP-seq, etc.) with the topological
 domains observed in Hi-C (see Fig3G of Lieberman-Aiden et al. 2009)
}

\value{
  A list with the eigen vector(s) of the \code{npc} first principal
  component(s), and their importance
}
\seealso{\code{\link{normPerExpected}}}
\author{N. Servant, B. Lajoie, R. McCord}
\examples{
## Get Lieberman-Aiden Hi-C data
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
import.my5C)
hiC <- HTClist(l)

## Performed PCA
pr<-pca.hic(hiC$chr14chr14, npc=1, asGRangesList=TRUE)

}
\keyword{manip}
