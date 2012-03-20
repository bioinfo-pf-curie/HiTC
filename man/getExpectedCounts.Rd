\name{getExpectedCounts}
\Rdversion{1.1}
\alias{getExpectedCounts}
\title{Estimate expected interaction counts of a High-Throughput C
  experiment based on the genomic distance between two loci}
\description{
The expected interaction is defined as the linear relationship between
the interaction counts and the distance between two primers. See details for addtional informations.
}
\usage{getExpectedCounts(x, span=0.01, bin=0.005, stdev=FALSE, plot=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{span}{fraction of the data used for smoothing at each x point.}
  \item{bin}{interpolation parameter}
  \item{stdev}{logical, calculate the variance}
  \item{plot}{logical, display loess smoothing and variance estimation points}
}

\details{
  The estimation of the background is based on the linear interpolation
  of the counts with the primers distances.  A lowess smoothing is
  used to estimate this linear relationship. Lowess uses robust locally
  linear fits. 
  A window is placed about each x value; points that are inside the
  window are weighted so that nearby points get the most weight (tricube
  weight function). 
  The lowess smoothing has two parameters : span (alpha) and bin (beta).
  The span corresponds to the fraction of the data used to for smoothing
  at each x point, i.e. to define the neighboring used for the local
  smoothing.
  The bin is the interpolotion parameter, and define the interval size
  in units corresponding to x. If lowess estimates at two x values
  within delta of one another, it fits any points between them by linear
  interpolation. The default is 1\% of the range of x. If delta=0 all but
  identical x values are estimated independently.
  The bin is used to speed up computation: instead of computing the local
  polynomial fit at each data point it is not computed for points within
  delta of the last computed point, and linear interpolation is used to fill in the fitted values for the skipped points. 
  This function may be slow for large numbers of points. Increasing bin
  should speed things up, as will decreasing span.

  The variance is then estimated using the same span and bin parameter,
  at each interpolation points.
}

\value{
  A list with the expected interaction map and the estimated variance
}
\seealso{\code{\link{HTCexp-class}},\code{\link{normPerZscore}},
  \code{\link{normPerExpected}}, \code{\link{lowess}}}
\author{N. Servant, B. Lajoie}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Estimate expected interaction from distance between intervals
GM12878.exp<-getExpectedCounts(GM12878$chr16chr16, stdev=TRUE, plot=FALSE)
mapC(GM12878.exp$exp.interaction)
}
\keyword{manip}
