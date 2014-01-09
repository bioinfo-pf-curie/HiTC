################
## New S4 class to represent list of HTCexp objects
################

## class definition
setClass("HTClist",
         prototype = prototype(elementType = "HTCexp"),
         contains = "list")


setValidity("HTClist",
            function(object){
                fails <- character(0)
                
                ## Check Chromosome Pairs Data
                 chrom <- sapply(object, function(x){
                     paste(seqlevels(x_intervals(x)), seqlevels(y_intervals(x)), sep="")
                 })
                 if (any(duplicated(chrom))){
                     fails <- c(fails, "Multiple 'HTCexp' found from the same chromosome pairs")
                 }

                ## Check Chromosome bounderies
                ##if (length(unique(unlist(ranges(object)))) > length(seqlevels(object))){
                ##    fails <- c(fails, "Same chromosome found with different ranges")
                ##}
                
                if (length(fails) > 0) return(fails)
                return(TRUE)
            }
)

## constructor
HTClist <- function(...)
{
  listData <- list(...)
  stopifnot(length(listData) > 0L)
  
  if (length(listData) == 1L && is.list(listData[[1L]]))
      listData <- listData[[1L]]
  if (!all(sapply(listData, is, "HTCexp")))
      stop("all elements in '...' must be HTCexp objects")
  listData <- unname(listData)

  names(listData) <- sapply(listData, function(x){
    paste(seqlevels(x_intervals(x)), seqlevels(y_intervals(x)), sep="")
  })
  
  new("HTClist", listData)
}


################
##
## Methods
##
################

setMethod("c", "HTClist", function(x, ...){
    if (missing(x))
        args <- unname(list(...))
    else{
        l <- list(...)
        l <- sapply(l, function(x){
            if (inherits(x,"HTClist"))
                x <- as.list(x)
            x
        })
        args <- unname(c(as.list(x), l))
    }
    HTClist(args)
})

setMethod("detail",signature(x="HTClist"),
          function(x){
              stopifnot(validObject(x))
              unlist(x)
          }
)

setMethod(f="isBinned", signature(x="HTClist"),
          function(x){
              sapply(x, isBinned)
          }
)


setMethod(f="isIntraChrom", signature(x="HTClist"),
          function(x){
              sapply(x, isIntraChrom) 
          }
)

setMethod(f="ranges", signature(x="HTClist"),
          function(x){
              GRangesList(sapply(x, range)) 
          }
)

setMethod(f="range", signature(x="HTClist"),
          function(x){
            reduce(unlist(ranges(x)))
          }
)

setMethod(f="seqlevels", signature(x="HTClist"),
          function(x){
              unique(unlist(unname(sapply(x, seqlevels))))
          }
)

setMethod("show",signature="HTClist",
          function(object){
              stopifnot(validObject(object))
              cat("HTClist object of length",length(object),"\n")
              im <- isIntraChrom(object)
              cat(length(which(im)),"intra /", length(which(!im)), "inter-chromosomal maps\n")
              invisible(NULL)
          }
)

## setMethod("summary", signature=c(object="HTClist"),
##           function(object){
              
##               xdata.intra <- xdata.inter <- NULL
##               nnzero.intra <- nnzero.inter <- NA_integer_
##               l.intra <- l.inter <- NA_integer_
              
##               ## Connot apply of summary for HTCexp - averge of averages different from global average
##               if (length(which(isIntraChrom(object)))>0){
##                   x.intra <- object[isIntraChrom(object)]
##                   xdata.intra <- unlist(lapply(x.intra, function(x) {
##                       xdata <- as(as(intdata(x), "sparseMatrix"), "dgTMatrix")
##                       return(xdata@x)
##                   }))
##                   nnzero.intra <- sapply(x.intra, function(x){nnzero(intdata(x), na.counted=TRUE)})
##                   l.intra <- sapply(x.intra, function(x){length(intdata(x))})
##               }
              
##               if (length(which(!isIntraChrom(object)))>0){
##                   x.inter <- object[!isIntraChrom(object)]
##                   xdata.inter <- unlist(lapply(x.inter, function(x) {
##                       xdata <- as(as(intdata(x), "sparseMatrix"),"dgTMatrix")
##                       return(xdata@x)
##                   }))
##                   nnzero.inter <- sapply(x.inter, function(x){nnzero(intdata(x), na.counted=TRUE)})
##                   l.inter <- sapply(x.inter, function(x){length(intdata(x))})
##               }
              
##               dstat <- matrix(0,ncol=5,nrow=3)
##               rownames(dstat) <- c("all","intra","inter")
##               colnames(dstat) <- c("nbreads","nbinteraction", "averagefreq","medfreq","meansparsity")
              
##               if (!is.null(xdata.intra)){
##                   mx.intra <- round(mean(xdata.intra, na.rm=TRUE),3)
##                   mdx.intra <- round(median(xdata.intra, na.rm=TRUE),3)
##                   spars.intra <- round(mean(nnzero.intra/l.intra, na.rm=TRUE),3)
##                   dstat[2,] <- c(sum(xdata.intra, na.rm=TRUE), length(xdata.intra), mx.intra, mdx.intra, spars.intra)
##               }
##               if (!is.null(xdata.inter)){
##                   mx.inter <- round(mean(xdata.inter, na.rm=TRUE),3)
##                   mdx.inter <- round(median(xdata.inter, na.rm=TRUE),3)
##                   spars.inter <- round(mean(nnzero.inter/l.inter, na.rm=TRUE),3)
##                   dstat[3,] <- c(sum(xdata.inter, na.rm=TRUE), length(xdata.inter), mx.inter, mdx.inter, spars.inter)
##               }
##               mx <- round(mean(c(xdata.intra, xdata.inter), na.rm=TRUE),3)
##               mdx <- round(median(c(xdata.intra, xdata.inter), na.rm=TRUE),3)
##               spars <- round(sum(c(nnzero.intra, nnzero.inter), na.rm=TRUE)/sum(c(l.intra, l.inter), na.rm=TRUE),3)
##               dstat[1,] <- c(sum(c(xdata.intra, xdata.inter), na.rm=TRUE), length(c(xdata.intra, xdata.inter)), mx, mdx, spars)
##               dstat
##           })

setMethod("summary", signature=c(object="HTClist"),
           function(object){
               t(sapply(object, summary))
           })


setMethod("[", "HTClist",
    function(x, i, ...)
    {
        if (is.character(i)){
            i <- match( i, names(x))
          }        
        HTClist(unlist(x)[i])
    }
)          

setMethod("as.list", "HTClist",
    function(x)
    {
        as.list(unlist(x))
    }
)          
