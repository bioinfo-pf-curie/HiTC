\name{importC}
\Rdversion{1.1}
\alias{importC}
\title{Import high-htroughput 'C' data}
\description{
  Import 5C or Hi-C data from csv file
}

\usage{importC(con, xgi.bed, ygi.bed = NULL, all.pairwise=TRUE)}

\arguments{
  \item{con}{input csv file. See details}
  \item{xgi.bed}{BED file describing the 'x' Intervals (i.e. column
    names) of the contact map. Required for the my5C list format}
  \item{ygi.bed}{BED file describing the 'y' intervals (i.e. row names) of the
    contact map. Required for the my5C list format}
  \item{all.pairwise}{logical; generate all pairwise chromosomal contact maps, i.e chr1-chr2, chr2-chr1}
}
\details{
  This function import high-throughput data from a matrix file.
}

\value{
  A list of \code{HTCexp} object(s)
}
\author{N. Servant}
\seealso{\code{\link{exportC}},\code{\link{import.my5C}}, \code{\link{HTCexp-class}}}
\examples{
\dontrun{
data(Nora_5C)

## Data binning
E14.bin<-binningC(E14$chrXchrX)

## Export the new intervals definition
exportC(E14.bin, file="E14")

##Import
importC("E14.csv")
}
}
\keyword{IO}
