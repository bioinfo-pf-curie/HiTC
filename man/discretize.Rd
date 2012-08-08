\name{discretize}
\Rdversion{1.1}
\alias{discretize}
\title{Transform matrix of counts data into discrete matrix}
\description{
  Transform matrix of counts data into discrete matrix
}
\usage{discretize(x, nb.lev=4, quant=TRUE)}
\arguments{
  \item{x}{data matrix}
  \item{nb.lev}{number of discretization level}
  \item{quant}{logical; use quantile distribution or split data into equals 'nb.lev' levels}
}
\value{
  A discrete matrix
}
\author{N. Servant}
\seealso{\code{quantile}}
\examples{
data(Nora_5C)

## Data binning
E14bin<-binningC(E14$chrXchrX)

## Discretize matrix
dismat<-discretize(intdata(E14bin))
mapC(dismat)
}
\keyword{manip}
