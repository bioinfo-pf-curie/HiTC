\name{CQC}
\Rdversion{1.1}
\alias{CQC}
\title{Quality Control for high-throughput 'C' experiment}
\description{
  Quality Control for high-throughput 'C' experiment
}
\usage{CQC(x,  cis.trans.ratio = TRUE, hist.interac=TRUE, scat.interac.dist=TRUE,
hist.dist=TRUE, dev.new=FALSE)}
\arguments{
  \item{x}{object or list of objects that inherits from class
    \code{HTCexp}}
  \item{cis.trans.ratio}{logical; barplot of percentage of
    inter-intrachromsomal interactions}
  \item{hist.interac}{logical; histogram of the interaction frequency}
  \item{scat.interac.dist}{logical; scatter plot of interaction count
    versus the genomic distance between two elements}
  \item{hist.dist}{logical; histogram of the distance between the 'x'
    and 'y' intervals}
  \item{dev.new}{logical; specifying if each plots must be in a separate
    graphical device}
}
\value{
    Return \code{NULL}; Create quality plots and print some additional informations
}
\author{N. Servant, B. Lajoie}
\seealso{\code{\link{HTCexp-class}}}
\examples{
exDir <- system.file("extdata", package="HiTC")
GM12878<-import.my5C(file.path(exDir,"nsmb.1936-S5.txt"), 
                 xgi.bed=file.path(exDir,"Bau_GM12878_REV.bed"), ygi.bed=file.path(exDir,"Bau_GM12878_FOR.bed")) 

## Quality Control
CQC(GM12878)
}
\keyword{manip}
