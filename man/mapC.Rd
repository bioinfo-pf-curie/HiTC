\name{mapC}
\Rdversion{1.1}
\alias{mapC}
\title{Visualize 'C' interaction map}
\description{
  Visualize 'C' interaction counts matrix
}

\usage{mapC(x, y=NULL, view=1, giblocs=NULL, minrange=NA, maxrange=NA, trim.range=0.98, names=FALSE, value=FALSE, show.na=FALSE, log.data=FALSE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", mask.data=NULL, grid=FALSE, title=NULL)}

\arguments{
  \item{x}{object that inherits from class \code{HTCexp} or from class
    \code{matrix}}
  \item{y}{optional. object that inherits from class \code{HTCexp} or
    from class \code{matrix}. If specified, view is set to 2}
  \item{view}{interaction map representation. See details}
  \item{giblocs}{genomeIntervals object of blocks to display as
  annotation track(s)}
  \item{minrange}{the minimum range of values used to define the color palette}
  \item{maxrange}{the maximum range of values used to define the color palette}
  \item{trim.range}{define the maxrange and minrange values using the
    percentile of the interaction matrix.}
  \item{names}{logical; display the names of the intervals. Useful for
    small matrices}
  \item{value}{logical; display the interaction values on the
    matrix. Useful for small matrices}
  \item{show.na}{logical; show the NA values in gray}
  \item{log.data}{logical; do you want to log the data before plotting
    the heatmap}
  \item{col.pos}{color for (low,mid,high) positive interaction counts. Must
  be a vectore of size 3. mid can be NA}
  \item{col.neg}{color for (low,mid,high) negative interaction counts. Must
  be a vectore of size 3. mid can be NA}
  \item{col.na}{color for NA values}
  \item{mask.data}{matrix to add to the heatmap as a mask. Must have the
  same dimension as the interaction matrix}
  \item{grid}{logical; add a grid on the heatmap}
  \item{title}{character; add a title to the heatmap}

}
\details{
  This function implements the \code{plot} method for objects
  of class \code{HTCexp}.
  
  By default, the maxrange and minrange values are fixed as the 98th
  percentile (resp. 2th percentile) of the interaction matrix. These
  values are useful to play with the contrast and remove the extreme
  values from the matrix.

  Two different views are available. The heatmap view (view=1) display
  the data in two dimension. The triangle view (view=2) only represent
  the top-right part the interaction matrix. If two \code{HTCexp}
  objects are specified the view is force to 2, in order to compare both
  interaction maps. The two maps have to be binned to ensure comparison
  between genomic ranges.

  Annotation tracks can be added to both views. In case of binned data,
  the exact genomic positions of each features are takken into
  account. Otherwise, the 'C' intervals which overlap with the
  annotation features are colored.
  
}
\value{
  Returns \code{NULL}; this function is called for the side-effect of
  creating the plot.
}
\author{N. Servant, B. Lajoie}
\seealso{\code{\link[genomeIntervals]{interval_overlap}}}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Interaction map
mapC(GM12878$chr16chr16)

## Play with contrast and color
mapC(GM12878$chr16chr16, maxrange=100, col.pos=c("black","red","yellow"))

## Add annotation and change view
exDir <- system.file("extdata", package="HiTC")
gene <- readBED(file.path(exDir,"refseq_hg19_chr16_1_500000.bed"))
mapC(GM12878$chr16chr16, giblocs=list(Refseq=gene$Refseq_Gene), view=2)

## Compare two samples
K562<-import.my5C(file.path(exDir,"nsmb.1936-S6.txt"), 
                 xgi.bed=file.path(exDir,"Bau_K562_REV.bed"), ygi.bed=file.path(exDir,"Bau_K562_FOR.bed")) 
mapC(binningC(GM12878$chr16chr16, binsize=50000, step=5),
binningC(K562$chr16chr16, binsize=50000, step=5), giblocs=list(Refseq=gene$Refseq_Gene))

}
\keyword{hplot}
