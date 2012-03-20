\name{HTCexp-class}
\Rdversion{1.1}
\docType{class}
\alias{HTCexp-class}
\alias{HTCexp}
\alias{export}
\alias{export-methods}
\alias{export,HTCexp,character,ANY-method}
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
\alias{id,Genome_intervals-method}
\alias{id<-}
\alias{id<--methods}
\alias{id<-,Genome_intervals,factor-method}
\alias{normPerReads}
\alias{normPerReads-methods}
\alias{normPerReads,HTCexp-method}
\alias{normPerExpected}
\alias{normPerExpected-methods}
\alias{normPerExpected,HTCexp-method}
\alias{normPerZscore}
\alias{normPerZscore-methods}
\alias{normPerZscore,HTCexp-method}
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
\alias{x_intervals<-,HTCexp,Genome_intervals-method}
\alias{y_intervals,HTCexp-method}
\alias{y_intervals}
\alias{y_intervals-methods}
\alias{y_intervals,HTCexp-method}
\alias{y_intervals<-}
\alias{y_intervals<--methods}
\alias{y_intervals<-,HTCexp,Genome_intervals-method}

\title{Class 'HTCexp'}
\description{
  A class for representing high throughput Chromosome Conformation Capture data from next-generation sequencing experiments.
}
\section{Objects from the Class}{
  Objects can be created either by:
  \enumerate{
    \item calls of the form
    \code{new("HTCexp", intdata, Genome_intervals, Genome_intervals)}.
    \item using the auxiliary function \code{HTCexp} and supplying
    interaction matrix with x and y intervals definition.
  }
}
\section{Slots}{
	 \describe{
    \item{\code{intdata}:}{Integer matrix, holding the
      interaction level between each pairs of 'x-y' intervals. The 'y' intervals must be in rows, and the 'x' in columns.}
    \item{\code{ygi}:}{Genomic interval of y intervals;
      see class \code{genome_intervals} for details}
    \item{\code{xgi}:}{Genomic interval of x intervals;
      see class \code{genome_intervals} for details}
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
\item{export}{create a BED file with the 'x' and 'y' intervals
  information.}
\item{isBinned}{return TRUE if the data are binned. The method tests
  if the 'x' and 'y' genome intervals are the same, if each bin has the
  same size and if the full genomic range is covered.}
\item{isIntraChrom}{return TRUE if the current
  \code{signature("HTCexp")} object contains intrachromosomal interaction data}
  \item{normPerReads}{normalize the interaction matrix by the total
    number of reads of the matrix.}
  \item{normPerExpected}{normalize the interaction matrix by the expected
    number of reads based on the distance between two loci.}
  \item{normPerZscore}{normalize the interaction matrix by the zscore
    calculation, which take into account the expected number of counts
  and the variance.}
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
\seealso{
  \code{\link[genomeIntervals]{Genome_intervals-class}},
  \code{\link[girafe]{AlignedGenomeIntervals-class}},
}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 
K562<-import.my5C(file.path(exDir,"nsmb.1936-S6.txt"), 
                 xgi.bed=file.path(exDir,"Bau_K562_REV.bed"), ygi.bed=file.path(exDir,"Bau_K562_FOR.bed")) 

## HTCexp description
show(GM12878)
detail(GM12878)

## Is binned data ?
isBinned(GM12878$chr16chr16)

## Is a inter or intrachromsomal experiment ?
isIntraChrom(GM12878$chr16chr16)

## Plotting
plot(GM12878$chr16chr16)
plot(GM12878$chr16chr16, view=2)
plot(binningC(GM12878$chr16chr16), binningC(K562$chr16chr16), maxrange=20)

## Zscore Normalization
GM12878norm<-normPerZscore(GM12878$chr16chr16)

## Operation on HTCexp object
GM12878_d_K562<-divide(normPerReads(GM12878$chr16chr16), normPerReads(K562$chr16chr16))
GM12878_s_K562<-substract(normPerReads(GM12878$chr16chr16), normPerReads(K562$chr16chr16))

## Overlap with genomic annotation
Refgene <- readBED(file.path(system.file("extdata", package="HiTC"),"refseq_hg19_chr16_1_500000.bed"))
plot(GM12878$chr16chr16, giblocs=list(RefSeqGene=Refgene$Refseq_Gene))
}
\keyword{classes}
