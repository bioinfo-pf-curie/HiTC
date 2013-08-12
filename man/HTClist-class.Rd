\name{HTClist-class}
\Rdversion{1.1}
\docType{class}
\alias{HTClist-class}
\alias{HTClist}
\alias{c,HTClist-method}
\alias{detail,HTClist-method}
\alias{as.list,HTClist-method}
\alias{isBinned,HTClist-method}
\alias{isIntraChrom,HTClist-method}
\alias{plot,HTClist,ANY-method}
\alias{range,HTClist-method}
\alias{ranges}
\alias{ranges-methods}
\alias{ranges,HTClist-method}
\alias{show,HTClist-method}
\alias{seqlevels,HTClist-method}
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
    \item{isBinned(x)}{applies 'isBinned' to each element in 'x'}
    \item{isIntraChrom(x)}{applies 'isIntraChrom' to each element in 'x'}
    \item{ranges(x)}{applies 'range' to each element in 'x'}
    \item{range(x)}{return the reduce range of all elements in 'x'}
    \item{seqlevels(x)}{return the sequence levels of all elements in 'x'}
    \item{as.list(x)}{coercion to simple \code{list} object}
    \item{names(x)}{get the names of the elements}
    \item{show(x)}{summarized output of the experiment, with informations
      about the data dimension}
    \item{x[i]}{Get elements \code{i} from \code{x}. Can be the
    positional index or its name.} 
  }
}
\author{Nicolas Servant}

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
isBinned(hiC)
isIntraChrom(hiC)
seqlevels(hiC)
}
\keyword{classes}
