\name{extractRegion}
\Rdversion{1.1}
\alias{extractRegion}
\title{Extract a subset of the HTCexp object}
\description{
  Extract a subset of the \code{HTCexp} object based on genomic ranges
}
\usage{extractRegion(x, from, to, exact=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
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
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Focus on the genomic region chrX:98000000-100000000
GM12878sub<-extractRegion(GM12878$chr16chr16, from=100000, to=300000)
GM12878sub
}
\keyword{manip}
