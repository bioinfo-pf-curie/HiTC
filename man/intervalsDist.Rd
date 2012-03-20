\name{intervalsDist}
\Rdversion{1.1}
\alias{intervalsDist}
\title{intervalsDist}
\description{
 Compute the distance between the intervals of a 'C' experiment
}
\usage{intervalsDist(x)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
}

\details{
  If \eqn{A} and \eqn{B} are the two sets of primers and
  \eqn{s} and \eqn{e}, the start and end of a primer, the distance is calculated as :\cr
  \deqn{\min(|A_e - B_s|, |A_s - B_e|)}
}

\value{
  A matrix of distances between primers
}
\seealso{\code{\link{HTCexp-class}}}
\author{N. Servant}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Calculate distances between primers/intervals
intervalsDist(GM12878$chr16chr16)
}
\keyword{manip}
