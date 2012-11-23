\name{exportC}
\Rdversion{1.1}
\alias{exportC}
\title{Export \code{HTCexp} object}
\description{
  Export \code{HTCexp} object to csv format
}
\usage{exportC(x, file)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{file}{character; the name of the output file}
}
\value{
  A csv file, with :\cr
  chrA,startA,endA,nameA,strandA,chrB,startB,endB,nameB,strandB,countAB
}
\author{N. Servant}
\seealso{\code{\link{export.my5C}}, \code{\link{importC}}}
\examples{
data(Nora_5C)

## Data binning
## E14.bin<-binningC(E14$chrXchrX)

## Export the new intervals definition
## exportC(E14.bin, file="E14.csv")
}
\keyword{IO}
