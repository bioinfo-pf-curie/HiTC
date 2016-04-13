\name{getPearsonMap}
\Rdversion{1.1}
\alias{getPearsonMap}
\title{Pearson correlation map}
\description{
  Generate pearson correlation map, usually used to call
    chromosomal compartments
}
\usage{getPearsonMap(x, normPerExpected=TRUE, center=TRUE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{normPerExpected}{normalized by expected interaction using the loess
    calculation of distance dependency. see normPerExpected}
  \item{center}{default=true. center the observed/expected map before calculating the
    Pearson correlation}
}

\details{
  The function returns an \code{HTCexp} object with Pearson correlation
  map. This is usually the first step of the Principal Component
  Analysis (see \code{pca.hic}).
  Centering the rows of the observed/expected matrix allows to avoid
  bias to due to ranges of interaction counts.
  If true, the correlation of small values should be as valuable as
  correlation of large values
}

\value{
  A \code{HTCexp} object
}
\seealso{\code{\link{normPerExpected}}}
\author{N. Servant, B. Lajoie, R. McCord}
\examples{
## Get Lieberman-Aiden Hi-C data
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
import.my5C)
hiC <- HTClist(l)

## get Pearson correlation map
pm <- getPearsonMap(hiC$chr14chr14)
mapC(HTClist(pm), maxrange=1, col.pos=c("black","red"), col.neg=c("black","blue"))
}
\keyword{manip}
