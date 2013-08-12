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
