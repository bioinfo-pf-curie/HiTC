\name{setIntervalScale}
\Rdversion{1.1}
\alias{setIntervalScale}
\title{Set x and y interval of the HTCexp object}
\description{Set x and y interval of the HTCexp object and update the
  interaction map accordingly}
\usage{setIntervalScale(x, xgi, ygi, upa=TRUE, method=c("median","mean","sum"),
use.zero=TRUE, optimize.by = c("speed", "memory"))}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{ygi}{y intervals;
    see class \code{GRanges} for details}
  \item{xgi}{x intervals;
    see class \code{GRanges} for details}
   \item{upa}{logical; unique primer assignment. Allow one primer to
    belong to one or several bins}
  \item{method}{the method used to combine the counts. Must be
    \sQuote{mean}, \sQuote{median} or \sQuote{sum}}
  \item{use.zero}{logical; use the zero values in the \code{method}
    calculation}
  \item{optimize.by}{"speed" will use faster methods but more RAM, and
      "memory" will be slower, but require less RAM}
}

\details{
  Define new interaction map based on the specified xgi and ygi
  intervals.

  This function has to be used carefully and can has important impact on
  the interaction map. 
  It is important to note that the \code{setIntervalScale} function is different from the \code{binningC} function in the way that the output
  is not symetrical.
}

\value{
  A \code{HTCexp} object
}
\seealso{\code{\link{HTCexp-class}}}
\author{N. Servant}
\examples{
data(Nora_5C)

E14.bin<-binningC(E14$chrXchrX)

## I have two HTCexp samples defined with different intervals.
show(E14.bin)
show(MEF$chrXchrX)

## How to compare them ?
## One idea is to force the intervals definition of one object using the
## intervals of the other.

setIntervalScale(MEF$chrXchrX, xgi=x_intervals(E14.bin), ygi=y_intervals(E14.bin))
}
\keyword{manip}
