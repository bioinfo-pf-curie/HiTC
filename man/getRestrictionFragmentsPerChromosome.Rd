\name{getRestrictionFragmentsPerChromosome}
\Rdversion{1.1}
\alias{getRestrictionFragmentsPerChromosome}

\title{Get a list of DNA restriction fragments}
\description{
  Performs the detection of restriction sites on a given genome and
  convert this information as a list of restriction fragments.}

\usage{getRestrictionFragmentsPerChromosome(resSite="AAGCTT", overhangs5=1,
chromosomes=NULL, genomePack="BSgenome.Mmusculus.UCSC.mm9")}

\arguments{
  \item{resSite}{the sequence of the restriction site to look
  for. Default is HindIII restriction site}
  \item{overhangs5}{5' overhangs on the DNA resulted from the cutting}
  \item{chromosomes}{vector of chromosome number to focus on. Default
  all the chromosomes for the specified genome.}
  \item{genomePack}{name of the genome package to access the DNA sequence}
}

\value{
  Returns a \code{GRanges} object with all restriction fragments for a
  given genome/chromosome.
}
\author{N. Servant}
\seealso{\code{\link{normLGF}}, \code{\link{setGenomicFeatures}}, \code{\link{getAnnotatedRestrictionSites}}}
\examples{

\dontrun{
## Mappability data From http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/
map_hg18<- import("wgEncodeCrgMapabilityAlign100mer.bw",format="BigWig", asRangedData=FALSE)

## 1- Get the list of restriction fragments for Human hg18 after HindIII digestion
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", overhangs5=1, chromosomes="chrX", genomePack="BSgenome.Hsapiens.UCSC.hg18")
}

}
\keyword{hplot}
