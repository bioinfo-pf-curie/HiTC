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
  \item{file}{character; the basename of the output file}
}
\value{
 Three output files will be created ; 2 BED files for each genomic
 intervals, and one matrix file
}
\author{N. Servant}
\seealso{\code{\link{export.my5C}}, \code{\link{importC}}}
\examples{
\dontrun{
data(Nora_5C)

## Data binning
E14.bin<-binningC(E14$chrXchrX)

## Export the new intervals definition
exportC(E14.bin, file="E14")
}
}
\keyword{IO}
