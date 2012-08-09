\name{extractRegion}
\Rdversion{1.1}
\alias{extractRegion}
\title{Extract a subset of the HTCexp object}
\description{
  Extract a subset of the \code{HTCexp} object based on genomic ranges
}
\usage{extractRegion(x, chr, from, to, exact=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{chr}{character; the chromosome of the genomic region}
  \item{from}{numeric;  start of the genomic region}
  \item{to}{numeric; end of the genomic region }
  \item{exact}{logical; exact genomic region}

}
\details{
  By default, only the intervals fully included in the genomic ranges are returned.
  If exact is true, the overlapping intervals are also used, and
  forced to start/end at the specified position. If no intervals are
  overlapping, an interval with NA values is added.
}
\value{
  A \code{HTCexp} object
}
\author{N. Servant}
\seealso{\code{\link[genomeIntervals]{Genome_intervals-class}}, \code{\link[girafe]{fracOverlap}}}
\examples{
data(Nora_5C)

## Focus on the genomic region chrX:98000000-100000000
E14sub<-extractRegion(E14$chrXchrX, chr="chrX", from=98000000, to=100000000)
E14sub
}
\keyword{manip}
