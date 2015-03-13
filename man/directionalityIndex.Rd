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
\examples{}
\keyword{IO}
