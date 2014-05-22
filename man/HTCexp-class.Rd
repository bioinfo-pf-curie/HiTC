\name{HTCexp-class}
\Rdversion{1.1}
\docType{class}
\alias{HTCexp-class}
\alias{HTCexp}
\alias{c,HTCexp-method}
\alias{detail}
\alias{detail-methods}
\alias{detail,HTCexp-method}
\alias{divide}
\alias{divide-methods}
\alias{divide,HTCexp,HTCexp-method}
\alias{intdata}
\alias{intdata-methods}
\alias{intdata,HTCexp-method}
\alias{intdata<--methods}
\alias{intdata<-}
\alias{intdata<-,HTCexp,Matrix-method}
\alias{isBinned}
\alias{isBinned-methods}
\alias{isBinned,HTCexp-method}
\alias{isIntraChrom}
\alias{isIntraChrom-methods}
\alias{isIntraChrom,HTCexp-method}
\alias{isSymmetric}
\alias{isSymmetric-methods}
\alias{isSymmetric,HTCexp-method}
\alias{id}
\alias{id-methods}
\alias{id,GRanges-method}
\alias{normPerReads}
\alias{normPerReads-methods}
\alias{normPerReads,HTCexp-method}
\alias{normPerExpected}
\alias{normPerExpected-methods}
\alias{normPerExpected,HTCexp-method}
\alias{normPerTrans}
\alias{normPerTrans-methods}
\alias{normPerTrans,HTCexp,HTCexp,HTCexp-method}
\alias{plot,HTCexp,ANY-method}
\alias{plot,HTCexp,HTCexp-method}
\alias{range}
\alias{range-method}
\alias{range,HTCexp-method}
\alias{seqlevels}
\alias{seqlevels-method}
\alias{seqlevels,HTCexp-method}
\alias{show,HTCexp-method}
\alias{substract}
\alias{substract-methods}
\alias{substract,HTCexp,HTCexp-method}
\alias{summary,HTCexp-method}
\alias{x_intervals}
\alias{x_intervals-methods}
\alias{x_intervals,HTCexp-method}
\alias{x_intervals<-}
\alias{x_intervals<--methods}
\alias{x_intervals<-,HTCexp,GRanges-method}
\alias{y_intervals,HTCexp-method}
\alias{y_intervals}
\alias{y_intervals-methods}
\alias{y_intervals,HTCexp-method}
\alias{y_intervals<-}
\alias{y_intervals<--methods}
\alias{y_intervals<-,HTCexp,GRanges-method}
\alias{xy_intervals,HTCexp-method}
\alias{xy_intervals}
\alias{xy_intervals-methods}

% Old stuff (deprecated or defunct).
\alias{export-defunct}
\alias{normPerZscore-defunct}
\alias{seq_name-deprecated}

\title{Class 'HTCexp'}
\description{
  A class for representing high throughput Chromosome Conformation Capture data from next-generation sequencing experiments.
}
\section{Objects from the Class}{
  Objects can be created either by:
  \enumerate{
    \item calls of the form
    \code{new("HTCexp", intdata, GRanges, GRanges)}.
    \item using the auxiliary function \code{HTCexp} and supplying
    contact Matrix with x and y intervals definition.
    The frceSymmetric option can used to force intra-chromosomal
    contact maps to be stored as symmetrical matrix.
  }
}
\section{Slots}{
	 \describe{
    \item{\code{intdata}:}{Dense or Sparse Matrix, holding the
      interaction level between each pairs of 'x-y' intervals. The 'y' intervals must be in rows, and the 'x' in columns.}
    \item{\code{ygi}:}{Genomic ranges of y intervals;
      see class \code{granges} for details}
    \item{\code{xgi}:}{Genomic ranges of x intervals;
      see class \code{granges} for details}
  }
}

\section{Methods}{
  \describe{
    \item{c(x, ...)}{Combines 'x' and the \code{signature("HTCexp")}
  objects in '...' together. The results is an object of class  \code{signature("HTCList")}}
    \item{detail(x)}{\code{signature("HTCexp")}: a more
      detailed output of the experiment than provided by \code{show}.}
    \item{divide(x)}{comparison of two \code{signature("HTCexp")}
      objects. Perform the division of the two contact matrices on the
      common 'x' and 'y' intervals. The operation is done only on the
      common intervals of both objects. If one of the two objects has a count
      to zero, the divided value will be NA}
    \item{intdata(x)}{return the \code{intdata} Matrix counts}
    \item{export(x)}{Defunct. See exportC method}
    \item{isBinned}{return TRUE if the data are binned. The method tests
      if the 'x' and 'y' genome intervals are the same, if 90\% of the bins have the
      same size and if the full genomic range is covered}
    \item{isIntraChrom(x)}{return TRUE if the current
      \code{signature("HTCexp")} object contains intrachromosomal contact data}
    \item{isSymmetric(x)}{return TRUE if the contact map is
      symmetrical, i.e inherits the \code{symmetricMatrix} class}
    \item{normPerReads(x)}{normalize the contact matrix by the total
      number of reads of the matrix.}
    \item{normPerExpected(x, stdev=TRUE)}{normalize the contact matrix by the expected
      number of reads based on the distance between two loci. See details.}
    \item{normPerZscore(x)}{Defunct. See normPerExpected method}
    \item{normPerTrans(x, xtrans, ytrans, method="sum")}{Normalize cis contact map based on the trans
      interactions. See details}
    \item{plot(x)}{visualization method; Display an heatmap of the
      contact data. Refer to the documentation of
      \code{\link{mapC}} for more details of the plotting function}
    \item{range(x)}{return the genomic range of the
      \code{signature("HTCexp")} object} 
    \item{seq_name(x)}{Defunct. See seqlevels method} 
    \item{seqlevels(x)}{return the sequence levels of the \code{signature("HTCexp")} object} 
    \item{show(x)}{summarized output of the experiment, with informations
      about the data dimension and the genomic region studied}
    \item{substract(x)}{comparison of two \code{signature("HTCexp")}
      objects. Perform the substraction of the two contact matrices on the
      common 'x' and 'y' intervals. The operation is done only on the
      common intervals of both objects. If one of the two objects has a count
      to zero, the divided value will be NA}
    \item{summary(x)}{return descriptive summary statistics about the contact map}
    \item{x_intervals(x)}{return the \code{xgi} GRanges object defining the
      x intervals}
    \item{y_intervals(x)}{return the \code{ygi} GRanges object defining the
      y intervals}
    \item{xy_intervals(x)}{return both \code{xgi} and \code{ygi} objects as a
      \code{GRangesList} object}
  }
}
\author{Nicolas Servant}
\details{
  The \code{normPerExpected} method estimates the expected interactions based on a the dependency on the genomic
  proximity between two loci. If \code{stdev} is false, the ratio observed/expected is returned,
  otherwise, the zscore ((observed-expected)/stdev) is returned
  Look at the   \code{\link{getExpectedCounts}} function for details.\cr
  
  The \code{normPerTrans} method is based on the assumption that all trans contacts should be the same. Thus, the cis contacts can be normalized by the interaction level of trans data.
  The \code{xtrans} trans map has to share its 'xgi' ranges with the
  cis map, and the \code{ytrans} has to share its 'ygi' ranges with
  the cismap.
  The \code{method} is used to combine the normalization factor from x
  and y ranges. Must be \sQuote{sum}, \sQuote{mult} or \sQuote{mean}.

}

\seealso{
  \code{\link[GenomicRanges]{GRanges-class}},\code{\link[GenomicRanges]{GRangesList-class}},\code{\link[Matrix]{Matrix-class}}
}

\examples{
data(Nora_5C)

## HTCexp descriptio
show(E14)
detail(E14)

## Is binned data ?
isBinned(E14$chrXchrX)

## Is a inter or intrachromsomal experiment ?
isIntraChrom(E14$chrXchrX)

## Divide by expected interaction counts
E14norm<-normPerExpected(E14$chrXchrX)

## Operation on HTCexp object
E14_d_MEF<-divide(normPerReads(E14$chrXchrX), normPerReads(MEF$chrXchrX))
E14_s_MEF<-substract(normPerReads(E14$chrXchrX), normPerReads(MEF$chrXchrX))

## Overlap with genomic annotation
require(rtracklayer)
gene <- import(file.path(system.file("extdata", package="HiTC"),"refseq_mm9_chrX_98831149_103425150.bed"), format="bed", asRangedData=FALSE)
plot(E14$chrXchrX, tracks=list(RefSeqGene=gene))

\dontrun{
## normPerTrans data normalization applied on \href{http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUmassDekker5C}{ENCODE data}.
ENCODE=import.my5C("./ENM-GM12878-R1.matrix")

## Look at raw contact map
mapC(ENCODE$chr7chr7)

## look at normalize by trans contact map
mapC(normPerTrans(ENCODE$chr7chr7, xtrans=ENCODE$chr7chr5, ytrans=ENCODE$chr5chr7))
}

\dontrun{
## Export
exportC(E14$chrXchrX, con="E14.csv")
}
}
\keyword{classes}
