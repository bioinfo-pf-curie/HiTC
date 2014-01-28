\name{normICE}
\Rdversion{1.1}
\alias{normICE}

\title{Iterative Correction of Hi-C data (ICE)}
\description{
 Iterative correction leverages the unique pairwise genome-wide structure of Hi-C data to decompose the data into a set of biases and a map of relative contact probabilities between any two genomic loci, achieving equal visibility across all genomic regions.
}

\usage{normICE(x,max_iter=200, eps=1e-4, sparse.filter=0.02)}

\arguments{
    \item{x}{object that inherits from class \code{HTCexp}}
    \item{max_iter}{maximum number of iteration}
    \item{eps}{the relative increment in the results before declaring
      convergence}
    \item{sparse.filter}{Define which pourcentage of bins with high
    sparsity should be force to zero}
}

\details{
  The imakaev normalization of Hi-C data consists of iteratively
  estimating the bias using the l1 norm.
  Note that the original method is applied on the genome-wide Hi-C
  map. For computation reason, we will applied the normalization to
  each chromosome separatly.
}
\value{
  Returns a \code{HTCexp} object with a corrected contact map.
}
\references{Imakaev M, Fudenberg G, McCord RP, Naumova N, Goloborodko A, Lajoie BR, Dekker J, Mirny LA. Iterative correction of Hi-C data reveals hallmarks of chromosome organization.Nat Methods. 2012 Oct;9(10):999-1003.}
\author{N. Servant, N. Varoqaux}
\seealso{\code{\link{normLGF}}}
\examples{
\dontrun{
##Lieberman data
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
            import.my5C)
hiC <- HTClist(l)
hiC <- hiC[isIntraChrom(hiC)]

## Run ICE
hiC_iced <- HTClist(lapply(hiC, normICE))
}
}

