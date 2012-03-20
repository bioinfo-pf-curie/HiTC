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
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 


## Remove intervals from a HTCexp object 
removeIntervals(GM12878$chr16chr16, ids=c("5C_305_ENm008_FOR_7", "5C_305_ENm008_REV_60"))
}
\keyword{manip}
