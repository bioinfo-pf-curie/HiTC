\name{Nora_5C}
\Rdversion{1.1}
\alias{Nora_5C}
\alias{E14}
\alias{MEF}
\docType{data}
\title{
HiTC - 5C data
}
\description{
  5C data described by Nora et al. (2012)
}
\usage{data(Nora_5C)}
\format{
  Contains two list of \code{HTCexp} objects (E14 and MEF). Data from the
  chromosome X are available.
}
\details{
This 5C dataset published by Nora et al
(\href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35721}{GSE35721}),
contains two different samples, a male undifferentiated ES cells (E14,
GSM873935) and a mouse embryonic fibroblasts (MEF, GSM873924). This
dataset is mainly used to describe the available functionalities of the \code{HiTC} package
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35721
}
\references{
Nora E, Lajoie B,
}
\examples{
data(Nora_5C)
show(E14)
show(MEF)
}
\keyword{datasets}