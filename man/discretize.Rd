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
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Data binning
GM12878bin<-binningC(GM12878$chr16chr16)

## Discretize matrix
dismat<-discretize(intdata(GM12878bin))
mapC(dismat)
}
\keyword{manip}
