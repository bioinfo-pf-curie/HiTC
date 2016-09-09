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
\usage{getExpectedCounts(x, method=c("mean","loess"), asList=FALSE, ...)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{method}{method used to estimate the expected counts based on the
    genomic distance. See details}
  \item{asList}{return the results as a list. Otherwise, return an
    object of class \code{HTCexp}}
  \item{...}{arguments for mean or loess method, see below}
}

\section{For 'mean' method}{
  \describe{
    \item{logbin}{if true logarithm based bins are used. In practice, it
  means that the bin size will change as we move away from the diagonal}
    \item{step}{multiplicative factor between each bin. Use if logbin is
  true}
    \item{filter.low}{fraction of low counts to filter before
  normalizing. default: 0.05.}
  }
}
\section{For 'loess' method}{
  \describe{
    \item{span}{fraction of the data used for smoothing at each x point.}
    \item{bin}{interpolation parameter}
    \item{stdev}{logical, calculate the variance}
    \item{plot}{logical, display lowess smoothing and variance estimation points}
  }
}

\details{
The expected value is the interaction frequency between two loci that one would expect
based on a sole dependency on the genomic proximity of these fragments in the linear
genome. This can be estimated using two different methods, mean or loess.

The first method (default) is simply based on the mean counts of each
diagonal. If logbin is false, the expected counts will be estimated for
each bin of the contact maps. If logbin is true, the binsize will change
according to the distance to the diagonal. Short (resp. long) distance will be
estimated with smaller (larger) bins. This method works for all resolutions.

The second method is based on a Lowess regression model and works
usually fine with low resolution data (250Kb to 1Mb). At higher
resolution, the lowess regression might be difficult to fit.
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
E14.exp<-getExpectedCounts(E14$chrXchrX, method="loess", stdev=TRUE, plot=FALSE)
mapC(E14.exp)
}
\keyword{manip}
