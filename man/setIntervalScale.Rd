\name{setIntervalScale}
\Rdversion{1.1}
\alias{setIntervalScale}
\title{Set x and y interval of the HTCexp object}
\description{Set x and y interval of the HTCexp object and update the
  interaction map accordingly}
\usage{setIntervalScale(x, xgi, ygi, upa=TRUE, method="mean", use.zero=TRUE)}
\arguments{
  \item{x}{object that inherits from class \code{HTCexp}}
  \item{ygi}{y intervals;
    see class \code{genome_intervals} for details}
  \item{xgi}{x intervals;
    see class \code{genome_intervals} for details}
   \item{upa}{logical; unique primer assignment. Allow one primer to
    belong to one or several bins}
  \item{method}{the method used to combine the counts. Must be
    \sQuote{mean}, \sQuote{median} or \sQuote{sum}}
  \item{use.zero}{logical; use the zero values in the \code{method} calculation}
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
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 
GM12878.bin<-binningC(GM12878$chr16chr16)

K562<-import.my5C(file.path(exDir,"nsmb.1936-S6.txt"), 
                 xgi.bed=file.path(exDir,"Bau_K562_REV.bed"), ygi.bed=file.path(exDir,"Bau_K562_FOR.bed")) 

## I have two HTCexp samples defined with different intervals.
show(GM12878.bin)
show(K562$chr16chr16)

## How to compare them ?
## One idea is to force the intervals definition of one object using the
## intervals of the other.

setIntervalScale(K562$chr16chr16, xgi=x_intervals(GM12878.bin), ygi=y_intervals(GM12878.bin))
}
\keyword{manip}
