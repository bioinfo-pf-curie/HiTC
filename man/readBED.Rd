\name{readBED}
\Rdversion{1.1}
\alias{readBED}
\title{readBED}
\description{
  read BED files and convert tracks in genomeIntervals objects
}
\usage{readBED(con)}
\arguments{
  \item{con}{BED file to read}
}

\details{
  If a score column is specified in the BED file, it will be saved as a
  'score' annotation slot.
}

\value{
  A list of GenomeIntervals object(s).
  Each element corresponds to a track.
}
\seealso{\code{\link{Genome_intervals-class}}}
\author{N. Servant}
\examples{
exDir <- system.file("extdata", package="HiTC")
gene <- readBED(file.path(exDir,"refseq_mm9_chrX_98831149_103425150.bed"))
}
\keyword{manip}
