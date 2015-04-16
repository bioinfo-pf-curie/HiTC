\name{getAnnotatedRestrictionSites}
\Rdversion{1.1}
\alias{getAnnotatedRestrictionSites}

\title{Annotation of restriction sites}
\description{
  Performs the annotation of all restriction sites of a given genome
  (i.e. GC content, mappability, effective fragment length)
}

\usage{getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1,
chromosomes=NULL, genomePack="BSgenome.Mmusculus.UCSC.mm9", mappability=NULL, wingc=200, winmap=500)}

\arguments{
  \item{resSite}{the sequence of the restriction site to look
  for. Default is HindIII restriction site}
  \item{overhangs5}{5' overhangs on the DNA resulted from the cutting}
  \item{chromosomes}{vector of chromosome number to focus on. Default
  all the chromosomes for the specified genome.}
  \item{genomePack}{name of the genome package to access the DNA sequence}
  \item{wingc}{size of the window upstream and downstream the restriction site
  used to calculate the GC content}
  \item{mappability}{a \code{GRanges} with a 'score' describing the
  mappability of the genome}
  \item{winmap}{size of the window upstream and downstream the restriction site
  used to calculate the mappability}
}


\details{
  This function automatically annotate all the restriction sites of a
  given chromosome. The mappability is optional but strongly advice for
  Hi-C contact map normalization. This information can be easily
  download from public ressources like ftp://hgdownload.cse.ucsc.edu/gbdb/mm9/bbi/.
}
\value{
  Returns a \code{GRanges} object annotation data upstream (U) and
  downstream (D) the restriction sites.
}
\author{N. Servant}
\seealso{\code{\link{normLGF}}, \code{\link{setGenomicFeatures}}}
\examples{

\dontrun{
## Mappability data From http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/
map_hg18<- import("wgEncodeCrgMapabilityAlign100mer.bw",format="BigWig", asRangedData=FALSE)

## 1- Example of restriction sites annnotation
cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1, chromosomes="chr1", genomePack="BSgenome.Hsapiens.UCSC.hg18", wingc=200, mappability=map_hg18, winmap=500)
}

}
\keyword{hplot}
