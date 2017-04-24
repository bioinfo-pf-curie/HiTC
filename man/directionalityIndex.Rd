\name{directionalityIndex}
\Rdversion{1.1}
\alias{directionalityIndex}
\title{Directionality index calculation}
\description{
  Calculate the directionality index as proposed by Dixon et al. 2012
}

\usage{directionalityIndex(x, winup = 2e+06, windown = 2e+06)}
\arguments{
  \item{x}{\code{HTClist} object}
  \item{winup}{size of upstream window}
  \item{windown}{size of downstrean window}
}
\details{
  Calculate the directionality index as proposed by Dixon et al. This
  index is then used to call topological domains in Hi-C/5C data.
}

\value{
  A \code{numeric} vector
}
\author{N. Servant}
\seealso{\code{\link{HTClist-class}}}
\examples{
require(HiCDataHumanIMR90)
data(Dixon2012_IMR90)
hox <- extractRegion(hic_imr90_40$chr6chr6, chr="chr6", from=50e6, to=58e6)
di<-directionalityIndex(hox)
}
\keyword{IO}
