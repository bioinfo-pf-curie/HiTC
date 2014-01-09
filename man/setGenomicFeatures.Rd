\name{setGenomicFeatures}
\Rdversion{1.1}
\alias{setGenomicFeatures}

\title{Annotation of Hi-C contact map}
\description{
  Annotate a Hi-C contact map with the genomic local features (i.e. GC content, mappability, effective fragment length)
}

\usage{setGenomicFeatures(x, ...)}

\arguments{
  \item{x}{\code{HTCexp} object to annotate}
  \item{...}{any option from the \code{getAnnotatedRestrictionSites} function
  (expect the chromosome)}
}


\details{
  The function will first call the \code{getAnnotatedRestrictionSites}
  to download the local genomic annotation, and will then used them to
  annotate the \code{HTCexp} object. The upstream and downstream
  annotation are combined (average or sum) according to their property.
}
\value{
  Returns a \code{HTCexp} object with local genomic features annotations.}
\author{N. Servant}
\seealso{\code{\link{normLGF}}, \code{\link{setGenomicFeatures}}}
\examples{

\dontrun{
require(BSgenome.Hsapiens.UCSC.hg18)

##Lieberman data
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
            import.my5C)
hiC <- HTClist(l)
hiC <- hiC[isIntraChrom(hiC)]
names(hiC)

## Mappability data From http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/
map_hg18<- import("wgEncodeCrgMapabilityAlign100mer.bw", format="BigWig", asRangedData=FALSE)

## Get the genomic feature of the HiC chr12 data
chr12_annot <- setGenomicFeatures(hiC$chr12chr12, resSite="AAGCTT", overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg18", wingc=200, mappability=map_hg18, winmap=500)
}
}

