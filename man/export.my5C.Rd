\name{export.my5C}
\Rdversion{1.1}
\alias{export.my5C}
\title{Export \code{HTCexp} object to my5C website format}
\description{
  Export \code{HTCexp} object to my5C website format
}
\usage{export.my5C(x, file)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{file}{character; the name of the output file}
}
\value{
  A my5C tabbed delimited file, with :\cr
  Y_INTERVAL_NAME/X_INTERVAL_NAME/INTERACTION_COUNT
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
