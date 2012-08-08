\name{removeIntervals}
\Rdversion{1.1}
\alias{removeIntervals}
\title{Remove intervals from HTC object}
\description{
  Remove primers intervals from HTC object
}
\usage{removeIntervals(x, ids)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{ids}{character; vector of primers Ids to remove from the object}
}
\value{
  A \code{HTCexp} object without the discarded intervals
}
\author{N. Servant}
\seealso{\code{\link[genomeIntervals]{Genome_intervals-class}}}
\examples{
data(Nora_5C)

## Remove intervals from a HTCexp object 
removeIntervals(E14$chrXchrX, ids=c("5C_938_XIC-3_REV_2", "5C_938_XIC-3_REV_4"))
}
\keyword{manip}
