\name{exportC}
\Rdversion{1.1}
\alias{exportC}
\title{Export \code{HTCexp} object}
\description{
  Export \code{HTCexp} object to tab format
}
\usage{exportC(x, file, per.chromosome=FALSE, use.names=FALSE, header=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{file}{character; the basename of the output file}
  \item{per.chromosome}{logical; export each contact maps in a different files (i.e one per chromosome pair)}
  \item{use.names}{if TRUE, keep the original row/colnames of the contact matrix}
  \item{header}{if TRUE, add an header with the package version and the date}
}
\value{
 Three output files will be created ; 2 BED files for each genomic
 intervals, and one tab file.
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
 If per.chromosome=FALSE, the data will be exported in one genome scaled file.
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
