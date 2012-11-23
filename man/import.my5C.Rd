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
  names) of the interaction map. Required for the my5C list format}
  \item{ygi.bed}{BED file describing the 'y' intervals (i.e. row names) of the
  interaction map. Required for the my5C list format}
  \item{all.pairwise}{logical; generate all pairwise chromosomal
  interaction maps, i.e chr1-chr2, chr2-chr1}
}
\details{
  This function allows data import from the \href{http://my5c.umassmed.edu}{the my5C webtool}.
  Two input formats can be used :
  -  The list format is composed of three files; two BED files
  describing the genomic intervals (i.e. primers); and a tabbed
  delimited format to specify the interaction between each genomic regions, with :\cr
  FORWARD_PRIMER_NAME/REVERSE_PRIMER_NAME/INTERACTION_COUNT
    
  - The matrix format is a tab-delimited format, corresponding to the
  interaction map. The rownames and columnames are splitted to created
  the genome intervals (example : REV_2|mm9|chrX:98831149-98834145). 
  
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
## Load my5C matrix format
hiC<-import.my5C(file.path(exDir,"HIC_gm06690_chr14_chr14_1000000_obs.txt"))
detail(hiC)
}
\keyword{IO}
