\name{setGenomicFeatures}
\Rdversion{1.1}
\alias{setGenomicFeatures}

\title{Annotation of Hi-C contact map}
\description{
  Annotate a Hi-C contact map with the genomic local features (i.e. GC content, mappability, effective fragment length)
}

\usage{setGenomicFeatures(x, cutSites, minFragMap=.5, effFragLen=1000)}

\arguments{
  \item{x}{\code{HTCexp} object to annotate}
  \item{cutSites}{\code{GRangesList} or \code{GRanges} object with restriction sites
    annotation obtained using the \code{getAnnotatedRestrictionSites}
    function}
  \item{minFragMap}{Minimum Fragment Mappbility. All fragments with a lower
    mappability are not used for the annotation.}
  \item{effFragLen}{Efective Fragment Length. Size of specific fragment ligation}
}


\details{
  The function require the restriction sites annotation as provided by
  the \code{getAnnotatedRestrictionSites} function.
  The restriction sites are first filtered according to their
  mappability. This threshold has to be defined according to the data
  pre-processing. All remaining restriction sites are then intersected
  with the genomic bins of the contact map. All restriction sites
  included within a bin are averaged.
  The effective fragment length is defined as the size of specific
  ligation product. (See Yaffe and Tanay, 2011). In this paper, the
  authors define specific ligation as sum of distance to cutter sites
  (d1+d2) <= 500 bp. Such criterion implies that d1<=500 bp and d2 <=
  500 bp. So for each fragment end, only reads mapped within 500 bp to
  cutter sites are used for downstream analysis.
  All defults paramters correspond to the ones used in the HiCNorm method.
  
}
\value{
  Returns a \code{HTCexp} object with local genomic features annotations.}
\author{N. Servant}
\seealso{\code{\link{normLGF}}, \code{\link{setGenomicFeatures}}}
\examples{

\dontrun{
require(BSgenome.Hsapiens.UCSC.hg18)
require(rtracklayer)

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
cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1, chromosomes=seqlevels(hiC), genomePack="BSgenome.Hsapiens.UCSC.hg18", wingc=200, mappability=map_hg18, winmap=500)

chr12_annot <- setGenomicFeatures(hiC$chr12chr12, cutSites)
}
}

