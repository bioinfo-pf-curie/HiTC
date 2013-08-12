\name{export.my5C}
\Rdversion{1.1}
\alias{export.my5C}
\title{Export \code{HTCexp} object to my5C website format}
\description{
  Export \code{HTCexp} object to my5C website format
}
\usage{export.my5C(x, file, format=c("mat","list"), genome="mm9", header=TRUE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{file}{character; the name of the output file}
  \item{format}{Either 'list' or 'mat'. See details}
  \item{genome}{The genome version. This information is only used for
    the 'mat' export format. See details}
  \item{header}{if true; a header is added in the output file with the
    package version and the date}
  
}
\details{
  If 'format=list',a my5C tabbed delimited file is created, with :\cr
  Y_INTERVAL_NAME/X_INTERVAL_NAME/INTERACTION_COUNT\cr
  Otherwise, a tab-delimited matrix file is generated with the row and colnames
  defined as follow as in the my5C web tool :\cr
  REV_2|mm9|chrX:98831149-98834145
}
\author{N. Servant}
\seealso{\code{\link{exportC}}}
\examples{
\dontrun{
data(Nora_5C)

## Data binning
E14.bin<-binningC(E14$chrXchrX)

## Export the new intervals definition
export.my5C(E14.bin, file="E14my5C.csv")
}
}
\keyword{plot}
