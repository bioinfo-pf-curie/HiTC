\name{HTClist-class}
\Rdversion{1.1}
\docType{class}
\alias{HTClist-class}
\alias{HTClist}
\alias{c,HTClist-method}
\alias{detail,HTClist-method}
\alias{as.list,HTClist-method}
\alias{forcePairwise}
\alias{forcePairwise,HTClist-method}
\alias{forceSymmetric,HTClist,character-method}
\alias{forceSymmetric,HTClist,missing-method}
\alias{getCombinedContacts}
\alias{getCombinedContacts-method}
\alias{getCombinedContacts,HTClist-method}
\alias{getCombinedIntervals}
\alias{getCombinedIntervals-method}
\alias{getCombinedIntervals,HTClist-method}
\alias{isComplete}
\alias{isComplete-methods}
\alias{isComplete,HTClist-method}
\alias{isPairwise}
\alias{isPairwise-method}
\alias{isPairwise,HTClist-method}
\alias{isBinned,HTClist-method}
\alias{isIntraChrom,HTClist-method}
\alias{plot,HTClist,ANY-method}
\alias{range,HTClist-method}
\alias{ranges}
\alias{ranges-methods}
\alias{ranges,HTClist-method}
\alias{reduce}
\alias{reduce-methods}
\alias{reduce,HTClist-method}
\alias{show,HTClist-method}
\alias{seqlevels,HTClist-method}
\alias{summary,HTClist-method}
\alias{[,HTClist,ANY,ANY-method}
\alias{[,HTClist,ANY,ANY,ANY-method}


\title{Class 'HTClist'}
\description{
  A class for representing a list of high throughput Chromosome Conformation Capture data from next-generation sequencing experiments.
}

\section{Constructor}{
  \describe{
    The \code{HTClist} represents a list of \code{HTCexp} objects and can
    be created as follow :
    \item{}{
      \code{HTClist(...)} : Creates a HTClist object using HTCexp objects supplied in '...'
    }
  }
}

\section{Methods}{
  \describe{
    \item{c(x, ...)}{Combines a \code{signature("HTClist")} object 'x'
      with \code{signature("HTClist")} or \code{signature("HTCexp")} 
      objects in '...'. The results is an object of class \code{signature("HTCList")}}
    \item{detail(x)}{\code{signature("HTClist")}: a more
      detailed output of the experiment than provided by \code{show}.}
    \item{forcePairwise(x)}{return a \code{signature("HTClist")} with
      all the pairwise contact maps}
    \item{forceSymmetric(x)}{return a \code{signature("HTClist")} with
       half of the pairwise contact maps}
     \item{isComplete(x)}{Logical; true if 'x' contains all intra and
      interchromosomal maps}
    \item{isPairwise(x)}{Logical; true if 'x' contains all interchromosomal
      pairs, i.e. chr1chr2 and chr2chr1}
    \item{isBinned(x)}{applies 'isBinned' to each element in 'x'}
    \item{isIntraChrom(x)}{applies 'isIntraChrom' to each element in
      'x'}
    \item{getCombinedContacts(x)}{merge all contact maps in a single big matrix}
    \item{getCombinedIntervals(x, merge=FALSE)}{merge all x and y intervals in single
  GRangesList object, or in a single GRanges object if merge=TRUE}
    \item{ranges(x)}{applies 'range' to each element in 'x'}
    \item{range(x)}{return the reduce range of all elements in 'x'}
    \item{reduce(x, chr, cis=TRUE, trans=TRUE)}{reduce a HTClist object
  to the list of provided chromosomes. Intra/Interchromosomal maps are
  returned accoring to the cis and trans args}
    \item{seqlevels(x)}{return the sequence levels of all elements in 'x'}
    \item{as.list(x)}{coercion to simple \code{list} object}
    \item{names(x)}{get the names of the elements}
    \item{show(x)}{summarized output of the experiment, with informations
      about the data dimension}
    \item{summary(x)}{return descriptive summary statistics for each interaction map}
    \item{x[i]}{Get elements \code{i} from \code{x}. Can be the
    positional index or its name.} 
  }
}
\author{Nicolas Servant}

\details{
  A \code{signature("HTClist")} is composed of a list of contact maps,
  representing the chromosomal interactions between pair of
  chromosomes.\cr
  The expected number of maps for a complete
  \code{signature("HTClist")} object should be equal to
  'lchrs+(lchrs*(lchrs-1)/2)'. In this case, the chr1-chr2 map is stored
  once, but the dataset is complete.\cr
  If a \code{signature("HTClist")} object is composed of all pairwise
  interaction maps, it means that the chr1-chr2 AND the chr2-chr1 maps
  will be stored. This way of storing the data is less memory efficient
  but can ease the use of some genome-wide algorithm.\cr
  Note that the \code{getCombinedContact} method should be used
  carefully. This method merges all single contact maps in one
  genome-wide map, therefore creating a very big matrix requiring memory
  space. However, this method is useful for many genome-wide analysis.
}

\seealso{
  \code{\link[GenomicRanges]{GRangesList-class}}, \code{\link{HTCexp-class}}
}

\examples{
exDir <- system.file("extdata", package="HiTC")
l <- sapply(list.files(exDir, pattern=paste("HIC_gm06690_"), full.names=TRUE),
            import.my5C)
hiC <- HTClist(l)
names(hiC)

## Methods
ranges(hiC)
range(hiC)
isComplete(hiC)
isPairwise(hiC)
isBinned(hiC)
isIntraChrom(hiC)
seqlevels(hiC)
}
\keyword{classes}
