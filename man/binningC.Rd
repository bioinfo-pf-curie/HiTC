\name{binningC}
\Rdversion{1.1}
\alias{binningC}
\title{Windowing of high-throughput 'C' interaction matrix}
\description{
  Windowing of 'C' interaction map
}
\usage{binningC(x, binsize=100000, bin.adjust=TRUE, upa=TRUE,
method="median", use.zero=TRUE, step=1, bnorm=FALSE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{binsize}{size of the bin to consider for windowing}
  \item{bin.adjust}{logical; adjust the size of the bin to the
   size of the genomic region}
  \item{upa}{logical; unique primer assignment. Allow one primer to
    belong to one or several bins}
  \item{method}{the method used to combine the counts. Must be
    \sQuote{mean}, \sQuote{median} or \sQuote{sum}}
  \item{use.zero}{logical; use the zero values in the \code{method} calculation}
  \item{step}{numeric; binning step size in \code{n} coverage \emph{i.e.} window step}
  \item{bnorm}{logical; normalise each combined counts by the number
    of primers in the bin}
}

\details{
 \code{bin.adjust} allows to work with bin of the same size. Otherwise,
 the last bin will has a size different from \code{binsize}.
 A primer is assigned to a bin, if there is at least one base overlap
 between the bin and the primer region.\cr
 The \code{method} used to combine the counts in a bin, must be \sQuote{mea}, \sQuote{median} or \sQuote{sum}.
 The \code{step} parameter allows to choose the overlap between the
 bins. A \code{step} of 2 means a 50\% overlap between two bins, a step
 of 3 means a 60\% overlap between two bins, \emph{etc.}
}

\value{
  An \code{HTCexp-class} object with binned intraction data. In this
  case, the primers are converted into bins, and the reverse or forward
  intervals are similar.
  The interaction matrix is symetric.
}
\seealso{\code{\link{HTCexp-class}}}
\author{N. Servant, B. Lajoie}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Data binning 50kb, with a 4/5 overlap
GM12878.bin <- binningC(GM12878$chr16chr16, binsize=50000, step=5)
show(GM12878.bin)
}
\keyword{manip}
