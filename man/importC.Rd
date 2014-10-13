\name{importC}
\Rdversion{1.1}
\alias{importC}
\title{Import high-htroughput 'C' data}
\description{
  Import 5C or Hi-C data from list file
}

\usage{importC(con, xgi.bed, ygi.bed=NULL, allPairwise=FALSE, forceSymmetric=FALSE)}
\arguments{
  \item{con}{input csv file. See details}
  \item{xgi.bed}{BED file describing the 'x' Intervals (i.e. column
    names) of the contact map}
  \item{ygi.bed}{BED file describing the 'y' intervals (i.e. row names) of the
    contact map}
  \item{allPairwise}{logical; generate all pairwise chromosomal contact
    maps, i.e chr1-chr2, chr2-chr1}
  \item{forceSymmetric}{logical; force the intra-chromosomal contact maps to be symmetrical}
}
\details{
  This function import high-throughput data from a tab list file.\cr
  The standard format for 5C/Hi-C data is the following :\cr
  ** One list file (tab delimited)
  bin1 bin2 x12
  bin1 bin3 x13
  ...
  ** The BED file(s) describing the intervals ('xgi.bed' and 'ygi.bed'
  are usually the same for Hi-C but can be different for 5C data)
  chr1 1 1000000 bin1
  chr1 1000001 2000000 bin2
  ...

  Note that this format is particularly interesting for sparse
  data as only non null values are stored.
}

\value{
  A \code{HTClist} object
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
importC(con="E14.count", xgi="E14_xgi.bed", ygi="E14_ygi.bed")
}
}
\keyword{IO}
