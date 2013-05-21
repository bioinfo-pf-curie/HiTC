\name{mapC}
\Rdversion{1.1}
\alias{mapC}
\alias{mapC-methods}
\alias{mapC,HTClist,ANY-method}
\alias{mapC,HTCexp,ANY-method}
\alias{mapC,HTCexp,HTCexp-method}

\title{Visualize 'C' interaction map}
\description{
  Visualize 'C' interaction map
}

\section{For \code{HTCexp} and \code{HTClist} objects}{
  \describe{
    \item{x}{object that inherits from class \code{HTCexp} or \code{HTClist}}
    \item{tracks}{List of GRanges objects of data to display as annotation track(s)}
    \item{minrange}{the minimum range of values used to define the color palette}
    \item{maxrange}{the maximum range of values used to define the color palette}
    \item{trim.range}{define the maxrange and minrange values using the
      percentile of the interaction matrix.}
    \item{show.na}{logical; show the NA values in gray}
    \item{log.data}{logical; do you want to log the data before plotting
      the heatmap}
    \item{col.pos}{color for (low,mid,high) positive interaction counts. Must
      be a vectore of size 3. mid can be NA}
    \item{col.neg}{color for (low,mid,high) negative interaction counts. Must
      be a vectore of size 3. mid can be NA}
    \item{col.na}{color for NA values}
      \item{grid}{logical; add a grid on the heatmap}
      \item{title}{character; add a title to the \code{HTCexp} plot(s)}
      \item{value}{logical; display the interaction values on the
	matrix. Useful for small matrices}
      \item{mask.data}{matrix to add to the heatmap as a mask. Must have the
      same dimension as the interaction matrix}
  }
}
\section{For \code{HTCexp} objects only}{
\describe{
  \item{y}{optional. object that inherits from class \code{HTCexp}.}
 }
}
\section{For \code{HTClist} objects only}{
\describe{
  \item{names}{logical; display the names of the intervals. Useful for
      small matrices}
  }
}


\details{
  This function implements the \code{plot} method for objects
  of class \code{HTCexp} and \code{HTClist}.
  
  By default, the maxrange and minrange values are fixed as the 98th
  percentile (resp. 2th percentile) of the interaction matrix. These
  values are useful to play with the contrast and remove the extreme
  values from the matrix.

  The \code{HTCexp} and \code{HTClist} are not represented in the same
  way. The heatmap view is used to display the \code{HTClist} objects in
  two dimension. This view is mainly useful to have an overview of the
  data, as Hi-C data.
  The triangle view is used for \code{HTCexp} only and represent
  the top-right part the interaction matrix. If two \code{HTCexp}
  objects are specified, they will be displayed in order to compare both
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
\seealso{\code{\link{HTCexp-class}}, \code{\link{HTClist-class}}}
\examples{
data(Nora_5C)

## Interaction map
## HTClist view
mapC(E14)

## HTCexp view
mapC(E14$chrXchrX)

## Play with contrast and color
mapC(E14$chrXchrX, maxrange=100, col.pos=c("black","red","yellow"))

## Add annotation and change view
require(rtracklayer)
exDir <- system.file("extdata", package="HiTC")
gene <- import(file.path(exDir,"refseq_mm9_chrX_98831149_103425150.bed"), format="bed", asRangedData=FALSE)
mapC(E14$chrXchrX, tracks=list(Refseq=gene))

## Compare two samples
mapC(binningC(E14$chrXchrX), binningC(MEF$chrXchrX), tracks=list(Refseq=gene))

}
\keyword{hplot}
