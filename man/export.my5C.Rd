\name{export.my5C}
\Rdversion{1.1}
\alias{export.my5C}
\title{Export \code{HTCexp} object to my5C website format}
\description{
  Export \code{HTCexp} object to my5C website format
}
\usage{export.my5C(x, outputfile)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{outputfile}{character; the name of the output file}
}
\value{
  A my5C tabbed delimited file (BED format), with :\cr
  Y_INTERVAL_NAME/X_INTERVAL_NAME/INTERACTION_COUNT
}
\author{N. Servant}
\seealso{\code{\link{export}}}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Data binning
GM12878.bin<-binningC(GM12878$chr16chr16, binsize=50000, step=5)

## Export the new intervals definition
export.my5C(GM12878.bin, outputfile="GM12878my5C.csv")
}
\keyword{plot}
