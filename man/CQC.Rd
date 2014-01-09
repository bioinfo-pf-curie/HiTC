\name{CQC}
\Rdversion{1.1}
\alias{CQC}
\title{Quality Control for high-throughput 'C' experiment}
\description{
  Quality Control for high-throughput 'C' experiment
}
\usage{CQC(x,  cis.trans.ratio = TRUE, hist.interac=TRUE, scat.interac.dist=TRUE,
hist.dist=TRUE, trim.range=0.98, winsize=NA)}
\arguments{
  \item{x}{object that inherits from class
    \code{HTCexp} or \code{HTClist}}
  \item{cis.trans.ratio}{logical; barplot of percentage of
    inter-intrachromsomal interactions}
  \item{hist.interac}{logical; histogram of the interaction frequency}
  \item{scat.interac.dist}{logical; scatter plot of interaction count
    versus the genomic distance between two elements}
  \item{hist.dist}{logical; histogram of the distance between the 'x'
    and 'y' intervals}
  \item{trim.range}{remove the extreme values by trimming the counts. Only use for plotting functions. [0,1]}
  \item{winsize}{used for the scat.interac.dist. If specify, the data
    are windowed before plotting}
}
\details{
If \code{x} is a \code{HTClist} object, all \code{HTCexp} objects are merged.
The zero values are not used to compute the descriptive statistics and to display the data.
If \code{trim.range} are lower than 1. The highest values (quantile probability is equal to \code{trim.range}) are discarded.
}
\value{
  Create quality plots and return a \code{matrix} with some simple statistics on all, cis and trans data.
}
\author{N. Servant, B. Lajoie}
\seealso{\code{\link{HTCexp-class}}}
\examples{
data(Nora_5C)

## Quality Control
CQC(E14)
}
\keyword{manip}
