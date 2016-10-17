\name{normLGF}
\Rdversion{1.1}
\alias{normLGF}


\title{Local Genomic Feature (LGF) normalization}
\description{
  Parametric model to remove systematic biases in the raw contact maps
}

\usage{normLGF(x, family=c("poisson", "nb"))}

\arguments{
    \item{x}{object that inherits from class \code{HTCexp}}
    \item{family}{parametric model to fit (poisson or nb)}
}

\details{
  This function implements the HiCNorm method proposed by Hu et al.
  Briefly, the method uses a generalized linear model to correct the
  systematic biases (effective fragment length, GC content, mappability)
  in a Hi-C contact map.
}
\value{
  Returns a \code{HTCexp} object with a normalized contact map.
}
\references{
Hu M, Deng K, Selvaraj S, Qin Z, Ren B, Liu JS. HiCNorm: removing biases in Hi-C data via Poisson regression. Bioinformatics. 2012;28(23):3131-3.
}
\author{N. Servant, M. Hu, S. Selvaraj}
\seealso{\code{\link{getAnnotatedRestrictionSites}}, \code{\link{setGenomicFeatures}}}
\examples{
\dontrun{
require(HiTC)
require(BSgenome.Hsapiens.UCSC.hg18)

##Lieberman data
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
            import.my5C)
hiC <- HTClist(l)
hiC <- hiC[isIntraChrom(hiC)]
names(hiC)

## Mappability data From http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/
map_hg18<- import("wgEncodeCrgMapabilityAlign100mer.bw", format="BigWig")

## Get the genomic feature of the chromosome 12
hiC_annot <- HTClist(lapply(hiC, setGenomicFeatures, resSite="AAGCTT", overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg18", wingc=200, mappability=map_hg18, winmap=500))
hiC_annot$chr12chr12

## Normalize the data
hiCnorm <- HTClist(lapply(hiC_annot, normLGF))

}
}
