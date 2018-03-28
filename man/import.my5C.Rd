\name{import.my5C}
\Rdversion{1.1}
\alias{import.my5C}
\title{Import data from my5C webtool}
\description{
  Import data from my5C webtool
}

\usage{import.my5C(file, allPairwise=FALSE, rm.trans=FALSE, lazyload=FALSE)}

\arguments{
  \item{file}{input file from the my5C webtool}
  \item{allPairwise}{logical; generate all pairwise chromosomal
    contact maps, i.e chr1-chr2, chr2-chr1}
  \item{rm.trans}{if true, only intra-chromosomal maps are loaded}
  \item{lazyload}{logical; force the intra-chromosomal contact maps to
  be stored as triangular matrix}
}
\details{
  This function allows data import from the \href{http://my5c.umassmed.edu}{the my5C webtool}.\cr
  The matrix format is a tab-delimited format, corresponding to the
  contact map. The rownames and columnames are splitted to created
  the genome intervals (example : REV_2|mm9|chrX:98831149-98834145). 
  
  The \code{allPairwise} option is not necessary in case of symetric design. Otherwise, it will return all the pairwise contact maps.

  The matrix will be stored as a matrix inheriting from \code{Matrix}
  class. If forcesymmetrical=TRUE, the intrachromosomal matrix as coerced to
  \code{symmetricMatrix} class allowing a much more efficient memory usage.
}

\value{
  A \code{HTClist} object
}
\author{N. Servant}
\seealso{\code{\link{HTClist-class}},\code{\link{HTCexp-class}},\code{\link{importC}},
  \code{\link{Matrix-class}}, \code{\link[Matrix]{symmetricMatrix-class}}
}
\examples{
exDir <- system.file("extdata", package="HiTC")
## Load my5C matrix format
hiC<-import.my5C(file.path(exDir,"HIC_gm06690_chr14_chr14_1000000_obs.txt"))
detail(hiC)
}
\keyword{IO}
