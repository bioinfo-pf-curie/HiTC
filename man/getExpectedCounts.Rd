\name{getExpectedCounts}
\Rdversion{1.1}
\alias{getExpectedCounts}
\title{Estimate expected interaction counts of a High-Throughput C
  intrachromsomal map based on the genomic distance between two loci
}
\description{
The expected interaction is defined as the linear relationship between
the interaction counts and the distance between two loci. See details for additional informations.
}
\usage{getExpectedCounts(x, span=0.01, bin=0.005, stdev=FALSE, plot=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{span}{fraction of the data used for smoothing at each x point.}
  \item{bin}{interpolation parameter}
  \item{stdev}{logical, calculate the variance}
  \item{plot}{logical, display lowess smoothing and variance estimation points}
}

\details{
The expected value is the interaction frequency between two loci that one would expect
based on a sole dependency on the genomic proximity of these fragments in the linear
genome. This can be estimated using a Lowess regression model.
The lowess smoothing has two parameters : span and bin.
The span corresponds to the fraction of the data used for smoothing.
Instead of computing the local polynomial fitting at each data point, a
window of size delta (bin parameter) is applied on the data and a linear interpolation
is used to fill in the fitted values within the window. 
The default is 1\% of the range of x. If delta=0 all but identical x values are estimated independently.
The variance is then estimated using the same span and bin parameter, at each interpolation points.
The points inside a window are weighted so that nearby points get the most weight (tricube weight function). 
}

\value{
  A list with the expected interaction map and the estimated variance
}
\seealso{\code{\link{HTCexp-class}},\code{\link{normPerExpected}},
  \code{\link{normPerExpected}}, \code{\link{lowess}}}
\author{N. Servant, B. Lajoie}
\examples{
data(Nora_5C)

## Estimate expected interaction from distance between intervals
E14.exp<-getExpectedCounts(E14$chrXchrX, stdev=TRUE, plot=FALSE)
mapC(E14.exp$exp.interaction)
}
\keyword{manip}
