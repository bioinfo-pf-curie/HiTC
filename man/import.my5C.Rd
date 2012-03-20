\name{import.my5C}
\Rdversion{1.1}
\alias{import.my5C}
\title{Import data from my5C webtool}
\description{
  Import data from my5C webtool
}

\usage{import.my5C(my5C.datafile, xgi.bed, ygi.bed, all.pairwise=TRUE)}

\arguments{
  \item{my5C.datafile}{input file from the my5C webtool}
  \item{xgi.bed}{BED file describing the 'x' Intervals (i.e. column
  names) of the interaction map}
  \item{ygi.bed}{BED file describing the 'y' intervals (i.e. row names) of the
  interaction map}
  \item{all.pairwise}{logical; generate all pairwise chromosomal interaction maps}
}
\details{
  The list format from \href{http://my5c.umassmed.edu}{the my5C webtool}
  is a tabbed delimited format (BED format), with :\cr
  FORWARD_PRIMER_NAME/REVERSE_PRIMER_NAME/INTERACTION_COUNT

  The matrix format is tabbed delimited format, corresponding to the
  interaction map. The rownames and columnames are splitted using the "|" separator
  (example : REV_2|mm9|chrX:98831149-98834145). The rownames and
  colnames are then intersected with the ids of the intervals defined in
  the BED files.
  
  The BED format is a standard format provided by the \href{http://genome.ucsc.edu/goldenPath/help/customTrack.html}{the UCSC Genome Browser}.

  The \code{all.pairwise} option is not necessary in case of symetric design. Otherwise, it will return all the pairwise interaction maps.
}

\value{
  A list of \code{HTCexp} object(s)
}
\author{N. Servant}
\seealso{\code{\link[genomeIntervals]{Genome_intervals-class}}, \code{\link{HTCexp-class}}}
\examples{
exDir <- system.file("extdata", package="HiTC")
hiC<-import.my5C(file.path(exDir,"HIC_gm06690_chr14_chr14_1000000_obs.txt"), xgi.bed=file.path(exDir,"GSE18199_gm06690_chr14_chr14_1Mb.bed"))
hiC

}
\keyword{IO}
