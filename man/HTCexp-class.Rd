\name{HTCexp-class}
\Rdversion{1.1}
\docType{class}
\alias{HTCexp-class}
\alias{HTCexp}
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
\alias{intdata<-,HTCexp,matrix-method}
\alias{isBinned}
\alias{isBinned-methods}
\alias{isBinned,HTCexp-method}
\alias{isIntraChrom}
\alias{isIntraChrom-methods}
\alias{isIntraChrom,HTCexp-method}
\alias{id}
\alias{id-methods}
\alias{id,GRanges-method}
\alias{id<-}
\alias{id<--methods}
\alias{id<-,GRanges,factor-method}
\alias{normPerReads}
\alias{normPerReads-methods}
\alias{normPerReads,HTCexp-method}
\alias{normPerExpected}
\alias{normPerExpected-methods}
\alias{normPerExpected,HTCexp-method}
\alias{normPerTrans}
\alias{normPerTrans-methods}
\alias{normPerTrans,HTCexp,HTCexp,HTCexp-method}
\alias{plot}
\alias{plot-methods}
\alias{plot,HTCexp,ANY-method}
\alias{plot,HTCexp,HTCexp-method}
\alias{range}
\alias{range-methods}
\alias{range,HTCexp-method}
\alias{seq_name}
\alias{seq_name-methods}
\alias{seq_name,HTCexp-method}
\alias{show}
\alias{show-methods}
\alias{show,HTCexp-method}
\alias{substract}
\alias{substract-methods}
\alias{substract,HTCexp,HTCexp-method}
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

% Old stuff (deprecated or defunct).
\alias{export}
\alias{export-methods}
\alias{export,HTCexp,character,ANY-method}
\alias{normPerZscore}
\alias{normPerZscore-methods}
\alias{normPerZscore,HTCexp-method}

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
    interaction matrix with x and y intervals definition.
  }
}
\section{Slots}{
	 \describe{
    \item{\code{intdata}:}{Integer matrix, holding the
      interaction level between each pairs of 'x-y' intervals. The 'y' intervals must be in rows, and the 'x' in columns.}
    \item{\code{ygi}:}{Genomic ranges of y intervals;
      see class \code{granges} for details}
    \item{\code{xgi}:}{Genomic ranges of x intervals;
      see class \code{granges} for details}
  }
}

\section{Methods}{
  \describe{
    \item{detail}{\code{signature("HTCexp")}: a more
      detailed output of the experiment than provided by \code{show}.}
    \item{divide}{comparison of two \code{signature("HTCexp")}
      objects. Perform the division of the two interaction matrices on the
      common 'x' and 'y' intervals. The operation is done only on the
      common intervals of both objects. If one of the two objects has a count
      to zero, the divided value will be NA.}
    \item{export}{Deprecated. See exportC function}
    \item{isBinned}{return TRUE if the data are binned. The method tests
      if the 'x' and 'y' genome intervals are the same, if each bin has the
      same size and if the full genomic range is covered}
    \item{isIntraChrom}{return TRUE if the current
      \code{signature("HTCexp")} object contains intrachromosomal interaction data}
    \item{normPerReads}{normalize the interaction matrix by the total
      number of reads of the matrix.}
    \item{normPerExpected}{normalize the interaction matrix by the expected
      number of reads based on the distance between two loci.}
    \item{normPerZscore}{Depracted. See normPerExpected}
    \item{normPerTrans}{Normalize cis interaction map based on the trans
      interactions. see details.}
    \item{plot}{visualization method; Display an heatmap of the
      interaction data. Refer to the documentation of
      \code{\link{mapC}} for more details of the plotting function.}
    \item{range}{return the genomic range of the \code{signature("HTCexp")} object} 
    \item{show}{summarized output of the experiment, with informations
      about the data dimension and the genomic region studied.}
    \item{substract}{comparison of two \code{signature("HTCexp")}
      objects. Perform the substraction of the two interaction matrices on the
      common 'x' and 'y' intervals. The operation is done only on the
      common intervals of both objects. If one of the two objects has a count
      to zero, the divided value will be NA.}
  }
}
\author{Nicolas Servant}
\details{
  The \code{normPerExpected} method estimates the expected interactions based on a the dependency on the genomic
  proximity between two loci. Look at the
  \code{\link{getExpectedCounts}} function for details.\cr
  
  The \code{normPerTrans} method is based on the assumption that all trans interactions should be the same. Thus, the cis interactions can be normalized by the interaction level of trans data.
  The \code{xtrans} trans map has to share its 'xgi' ranges with the
  cis map, and the \code{ytrans} has to share its 'ygi' ranges with
  the cismap.
  The \code{method} is used to combine the normalization factor from x
  and y ranges. Must be \sQuote{sum}, \sQuote{mult} or \sQuote{mean}.

}

\seealso{
  \code{\link[GenomicRanges]{GRanges-class}},
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

## ## Plotting
## plot(E14$chrXchrX)
## plot(E14$chrXchrX, view=2)
## plot(binningC(E14$chrXchrX), binningC(MEF$chrXchrX), maxrange=20)

## Divide by expected interaction counts
E14norm<-normPerExpected(E14$chrXchrX)

## Operation on HTCexp object
E14_d_MEF<-divide(normPerReads(E14$chrXchrX), normPerReads(MEF$chrXchrX))
E14_s_MEF<-substract(normPerReads(E14$chrXchrX), normPerReads(MEF$chrXchrX))

## Overlap with genomic annotation
require(rtracklayer)
gene <- import(file.path(system.file("extdata", package="HiTC"),"refseq_mm9_chrX_98831149_103425150.bed"), format="bed", asRangedData=FALSE)
plot(E14$chrXchrX, giblocs=GRangesList(RefSeqGene=gene))

\dontrun{
## normPerTrans data normalization applied on \href{http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUmassDekker5C}{ENCODE data}.
ENCODE=import.my5C("./ENM-GM12878-R1.matrix")

## Look at raw interaction map
mapC(ENCODE$chr7chr7)

## look at normalize by trans interaction map
mapC(normPerTrans(ENCODE$chr7chr7, xtrans=ENCODE$chr7chr5, ytrans=ENCODE$chr5chr7))
}

\dontrun{
## Export
exportC(E14$chrXchrX, con="E14.csv")
}
}
\keyword{classes}
