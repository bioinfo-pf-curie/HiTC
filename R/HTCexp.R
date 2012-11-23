################
## New S4 class for representing 5C or HiC data
## This class is mainly based on the GRanges class
## The interaction matrix is added to the object 
################

## class definition
setClass("HTCexp",
         representation = representation(
         "intdata"="matrix",
         "xgi"="GRanges",
         "ygi"="GRanges"),
         validity = function(object) {
             fails <- character(0)

             ## Format
             if(class(object@intdata)!="matrix") {
                 fails <- c(fails, "Expected matrix for the data slot!\n")
             }
             if(class(object@ygi)!="GRanges" || class(object@xgi)!="GRanges") {
                 fails <- c(fails, "Intervals 'xgi' and 'ygi' objects must be two GRanges objects!\n")
             }
             if (dim(object@intdata)[1] != length(object@ygi)){
                 fails <- c(fails, "The 'ygi' GRanges object must have exactly the same length as the number of rows of the data matrix")
             }
             if (dim(object@intdata)[2] != length(object@xgi)){
                 fails <- c(fails, "The 'xgi' GRanges object must have exactly the same length as the number of column of the data matrix")
             }
             if (is.null(colnames(object@intdata)) || is.null(id(object@xgi))){
                 fails <- c(fails,"No ids for 'xgi' intervals")
             }
             if (is.null(rownames(object@intdata)) || is.null(id(object@ygi))){
                 fails <- c(fails,"No ids for 'ygi' intervals")
             }
             
             ## Content
             if (length(object@xgi) == 0 || length(object@ygi) == 0){
                 fails <- c(fails, "Inervals xgi or ygi are of size 0")
             }
             if (length(chromosome(object@xgi)) > 1 || length(chromosome(object@ygi)) > 1){
                 fails <- c(fails, "Multiple chromosome found in xgi or ygi object")
             }
             if (length(setdiff(rownames(object@intdata), id(object@ygi)))>0){
                 fails <- c(fails, "The 'ygi' intervals from the interaction matrix are different from those defined in the GRanges object")
             }
             if (length(setdiff(colnames(object@intdata), id(object@xgi)))>0){
                 fails <- c(fails, "The 'xgi' intervals from the interaction matrix are different from those defined in the GRanges object")
             }
             
             if (length(fails) > 0) return(fails)
             return(TRUE)
         }
         )


## constructor
HTCexp <- function(intdata, xgi, ygi){
  
    ## check xgi/ygi length
    if (length(xgi)==0 || length(ygi)==0){
        stop("Cannot create HTCexp object of size 0")
    }
    
    ## check format
    stopifnot(is.matrix(intdata))
    stopifnot(inherits(ygi,"GRanges"))
    stopifnot(inherits(xgi,"GRanges"))

    ## sort xgi and ygi data
    xgi <- sort(xgi)
    ygi <- sort(ygi)

    ## deal with seqlevels in GRanges objects to avoid warnings
    seqlevels(xgi) <- unique(c(seqlevels(xgi), seqlevels(ygi)))
    seqlevels(ygi) <- unique(c(seqlevels(xgi), seqlevels(ygi)))

    intdata <- intdata[id(ygi), id(xgi)]
    new("HTCexp", intdata, xgi, ygi)
}

################
##
## Functions
##
################
removeIntervals <- function(x, ids){
    stopifnot(inherits(x, "HTCexp"))
    xgi <- x_intervals(x)
    ygi <- y_intervals(x)
    xgi.id <- id(xgi)
    ygi.id <- id(ygi)
    stopifnot(!is.null(xgi.id) | !is.null(ygi.id))

    xgi.id.n <- setdiff(xgi.id, ids)
    if (length(xgi.id.n) != length(xgi.id)){
        message("Discard ",(length(xgi.id)-length(xgi.id.n)), " 'x' intervals")
        x_intervals(x) <- xgi[which(id(xgi)%in%xgi.id.n),]
    }
    
    ygi.id.n <- setdiff(ygi.id, ids)
    if (length(ygi.id.n) != length(ygi.id)){
        message("Discard ",(length(ygi.id)-length(ygi.id.n)), " 'y' intervals")
        y_intervals(x) <- ygi[which(id(ygi)%in%ygi.id.n),]
    }
    x
}

extractRegion <- function(x, MARGIN=c(1,2), chr, from, to, exact=FALSE){
    stopifnot(inherits(x, "HTCexp"))
    stopifnot(is.element(MARGIN, c(1,2)))

    fromto <- GRanges(seqnames = chr, ranges = IRanges(start=from, end=to))
    xgi <- x_intervals(x)
    ygi <- y_intervals(x)
    data <-  x@intdata

    ## Force exact overlap by adding flanking regions
    if (exact){
        ## if no interval with from, add NA values
        ## else force overlapping coordinates
        ## x_intervals
        if (is.element(c(1), MARGIN)){
            xgi <- subsetByOverlaps(xgi, fromto, ignore.strand=TRUE)
            
            if (min(start(xgi))>from){
                ##flank(xgi[which.min(start(xgi))], width=min(start(xgi))-from, start=TRUE, ignore.strand=TRUE)
                ir <- IRanges(start=from, end=min(start(xgi))-1)
                edata <- elementMetadata(xgi)[1,]
                if (is.element("score",colnames(edata)))
                    edata$score <- 0
                if (is.element("name",colnames(edata)))
                    edata$name <- paste(chr,":",from,"-",min(start(xgi))-1, sep="")
                if (is.element("thick",colnames(edata)))
                    edata$thick <- ir
                
                ngr <- GRanges(seqnames = chr, ranges = ir, strand = "*", edata)
                xgi <- c(ngr,xgi)
                data <- cbind(rep(NA, nrow(data)), data)
            }else if (min(start(xgi))<=from){
                maxind <- max(which(start(xgi)<=from & end(xgi)>=from))
                start(xgi[maxind,]) <- from
                data <- data[,maxind:length(xgi)]
                xgi <- xgi[maxind:length(xgi),]
            }
            if (max(end(xgi))<to){
                ir <- IRanges(start=max(end(xgi))+1, end=to)
                edata <- elementMetadata(xgi)[length(xgi),]
                if (is.element("score",colnames(edata)))
                    edata$score <- 0
                if (is.element("name",colnames(edata)))
                    edata$name <- paste(chr,":",max(end(xgi))+1,"-",to, sep="")
                if (is.element("thick",colnames(edata)))
                    edata$thick <- ir
                
                ngr <- GRanges(seqnames = chr, ranges = ir, strand = "*", edata)
                xgi <- c(xgi, ngr)
                data <- cbind(data, rep(NA, nrow(data)))
            }else if (max(end(xgi))>=to){
                maxind <- min(which(start(xgi)<=to & end(xgi)>=to))
                end(xgi[maxind,]) <- to
                data <- data[,1:maxind]
                xgi <- xgi[1:maxind,]
            }
        }
        
        ## y_intervals
        if (is.element(c(2), MARGIN)){
            ygi <- subsetByOverlaps(ygi, fromto, ignore.strand=TRUE)
            
            if (min(start(ygi))>from){
                ir <- IRanges(start=from, end=min(start(ygi))-1)
                edata <- elementMetadata(ygi)[1,]
                if (is.element("score",colnames(edata)))
                    edata$score <- 0
                if (is.element("name",colnames(edata)))
                    edata$name <- paste(chr,":",from,"-",min(start(ygi))-1, sep="")
                if (is.element("thick",colnames(edata)))
                    edata$thick <- ir
             
                ngr <- GRanges(seqnames = chr, ranges = ir, strand = "*", edata)
                ygi <- c(ngr,ygi)
                data <- rbind(rep(NA, ncol(data)), data)
            }else if (min(start(ygi))<=from){
                maxind <- max(which(start(ygi)<=from & end(ygi)>=from))
                start(ygi[maxind,]) <- from
             data <- data[maxind:length(ygi),]
                ygi <- ygi[maxind:length(ygi),]
            }        
            if (max(end(ygi))<to){
                ir <- IRanges(start=max(end(ygi))+1, end=to)
                edata <- elementMetadata(ygi)[length(ygi),]
                if (is.element("score",colnames(edata)))
                    edata$score <- 0
                if (is.element("name",colnames(edata)))
                    edata$name <- paste(chr,":",max(end(ygi))+1,"-",to, sep="")
                if (is.element("thick",colnames(edata)))
                    edata$thick <- ir
                
                ngr <- GRanges(seqnames = chr, ranges = ir, strand = "*", edata)
                ygi <- c(ygi, ngr)
                data <- rbind(data, rep(NA, ncol(data)))
            }else if (max(end(ygi))>=to){
                maxind <- min(which(start(ygi)<=to & end(ygi)>=to))
                end(ygi[maxind,]) <- to
                data <- data[1:maxind,]
                ygi <- ygi[1:maxind,]
            }
        }
        
        data <-  x@intdata[id(ygi), id(xgi)]
        colnames(data) <- id(xgi)
        rownames(data) <- id(ygi)
        HTCexp(data, xgi, ygi)
     }else{
         ## Full overlap
         if (is.element(c(1), MARGIN))
             xgi <- subsetByOverlaps(xgi, fromto, type="within", ignore.strand=TRUE)
         
         if (is.element(c(2), MARGIN))
             ygi <- subsetByOverlaps(ygi, fromto, type="within", ignore.strand=TRUE)
         
         data <-  x@intdata[id(ygi), id(xgi)]
         HTCexp(data, xgi, ygi)
     }
}##extractRegion


################
##
## Methods
##
################
setMethod("chromosome",signature(x="GRanges"), function(x){
    as.vector(unique(seqnames(x)))
})


setMethod("detail",signature(x="HTCexp"), function(x){
    stopifnot(validObject(x))
    cat("HTC object\n")

    r <- range(x)
    if (length(r)==1){
        cat("Focus on genomic region [", chromosome(r), ":",
                start(r),"-",end(r),"]\n", sep="")
    }else{
        cat("Focus on genomic regions [",
                chromosome(r)[1], ":", start(r)[1], "-", end(r)[1],"] vs [",
                chromosome(r)[2], ":", start(r)[2], "-", end(r)[2],"]\n", sep="")
    }
   
    if(isIntraChrom(x)){
        cat("CIS Interaction Map\n")
    }else{
        cat("TRANS Interaction Map\n")
    }
    cat("Matrix of Interaction data: [", dim(intdata(x))[1],"-",dim(intdata(x))[2], "]\n", sep="")
    if (isBinned(x)){
        cat("Binned data - window size =", width(x_intervals(x))[1],"\n")
        cat(length(x_intervals(x)),"genome intervals\n")
    }
    else{
        cat(length(x_intervals(x))," genomic ranges from 'xgi' object\n")
        cat(length(y_intervals(x))," genomic ranges from 'ygi' object\n")
    }
    data <- as.vector(intdata(x))
    cat("Total Reads = ",sum(data, na.rm=TRUE),"\n")
    cat("Number of Interactions = ",length(data[which(data>0)]),"\n")
    cat("Median Frequency = ",median(data[which(data>0)]),"\n")
    invisible(NULL)
}) 

## Divide method
setMethod("divide", signature=c("HTCexp","HTCexp"), definition=function(x,y){
    xgi <- subsetByOverlaps(x@xgi, y@xgi, type="equal")
    ygi <- subsetByOverlaps(x@ygi, y@ygi, type="equal")

    if (length(xgi) == 0 || length(ygi) == 0){
        stop("No equal intervals found in x and y objects")
    }
    
    a <- x@intdata[id(ygi),id(xgi)]
    b <- y@intdata[id(ygi), id(xgi)]
    a[which(a==0 | b==0)]<-NA
    b[which(a==0 | b==0)]<-NA
    data <- a/b
    HTCexp(data, xgi , ygi)
})


## Interaction Data
setMethod(f="intdata", signature(x="HTCexp"), definition=function(x){
    x@intdata
})
  
setReplaceMethod(f="intdata", signature(x="HTCexp",value="matrix"),function(x, value) {
    x@intdata <- value
    return(x)
})

## Initialize method
setMethod("initialize",signature=c("HTCexp"), function(.Object, interaction.data, xgi, ygi){

    if (!missing(interaction.data) && !missing(ygi) && !missing(xgi)){
        .Object@intdata<-interaction.data
        .Object@ygi <-ygi
        .Object@xgi <-xgi
    }
    return(.Object)
})

## Binned data
setMethod("isBinned", signature(x="HTCexp"), function(x){
    stopifnot(validObject(x))
    ret <- TRUE

    ## Same bins for intra chromosomal
    if(isIntraChrom(x)){
        if (length(setdiff(ranges(x@xgi),ranges(x@ygi))) != 0)
            ret <- FALSE
    }
    ## Same width (exept for last bins - checked in 90% of bins)
    if (length(unique(width(x@xgi)[1:round(length(x@xgi)*.9)])) != 1 || length(unique(width(x@ygi)[1:round(length(x@ygi)*.9)])) != 1)
        ret <- FALSE
    ## no gaps
    if (length(reduce(x@xgi)) != 1 || length(reduce(x@ygi)) != 1)
        ret <- FALSE
    
    ret
})

## Inter/Intra chromosomal interaction
setMethod("isIntraChrom", signature(x="HTCexp"), function(x){
    stopifnot(validObject(x))
    ret <- FALSE
    if (chromosome(x@xgi) == chromosome(x@ygi))
        ret <- TRUE
    ret
})

## Intervals names
setMethod(f="id", signature(x="GRanges"), definition=function(x){
    ret <- NULL
    if (!is.null(names(x)))
        ret <- names(x)
    else if (!is.null(elementMetadata(x)$name))
        ret <- elementMetadata(x)$name
    ret
})

#setReplaceMethod(f="id", signature(x="GRanges", value="character"),function(x, value) {
#    elementMetadata(x)$name <- value
#    return (x)
#})

## Plot method
setMethod("plot", signature="HTCexp", definition=function(x, ...) {
    mapC(x, ...)
})

setMethod("plot", signature=c("HTCexp","HTCexp"), definition=function(x, y, ...) {
    mapC(x, y, ...)
})


## setMethod(f="range", signature(x="HTCexp"), definition=function(x){
##     if (chromosome(x@xgi) == chromosome(x@ygi)){
##         r <- range(c(ranges(x@ygi), ranges(x@xgi)))
##         ret <- c(chromosome(x@xgi), start(r), end(r))
##     }else{
##         ry <- range(x@ygi)
##         rx <- range(x@xgi)
##         ret <- as.data.frame(matrix(c(chromosome(x@ygi), start(ry), end(ry),  chromosome(x@xgi), start(rx), end(rx)), nrow=2, byrow=TRUE))
##         colnames(ret) <- c("chrom","start","end")
##     }
##     ret
## })

setMethod(f="range", signature(x="HTCexp"), definition=function(x){
    range(c(x@ygi, x@xgi), ignore.strand=TRUE)
})


## Seqlevels
setMethod(f="seqlevels", signature(x="HTCexp"), definition=function(x){
    unique(c(seqlevels(x@xgi), seqlevels(x@ygi)))
})

## Show method
setMethod("show",signature="HTCexp", function(object){
    stopifnot(validObject(object))
    cat("HTC object\n")

    r <- range(object)
    if (length(r)==1){
        cat("Focus on genomic region [", chromosome(r), ":",
                start(r),"-",end(r),"]\n", sep="")
    }else{
        cat("Focus on genomic regions [",
                chromosome(r)[1], ":", start(r)[1], "-", end(r)[1],"] vs [",
                chromosome(r)[2], ":", start(r)[2], "-", end(r)[2],"]\n", sep="")
    }
    #if(isIntraChrom(object)){
    #    cat("CIS Interaction Map")
    #}else{
    #    cat("TRANS Interaction Map")
    #}
    cat("Matrix of Interaction data: [", dim(intdata(object))[1],"-",dim(intdata(object))[2], "]\n", sep="")
    #if (isBinned(object)){
    #    cat("Binned data - window size =", y_intervals(object)[1,2]-y_intervals(object)[1,1]+1)
    #    cat(nrow(y_intervals(object)),"genome intervals")
    #}
    #else{
    #    cat(length(x_intervals(object))," genome intervals from 'xgi' object")
    #    cat(length(y_intervals(object))," genome intervals from 'ygi' object")
    #}
    invisible(NULL)
})


## Substraction
setMethod("substract", signature(x="HTCexp",y="HTCexp"), definition=function(x,y){
    xgi <- subsetByOverlaps(x@xgi, y@xgi, type="equal")
    ygi <- subsetByOverlaps(x@ygi, y@ygi, type="equal")
 
    a <- x@intdata[id(ygi),id(xgi)]
    b <- y@intdata[id(ygi), id(xgi)]
    a[which(a==0 | b==0)]<-NA
    b[which(a==0 | b==0)]<-NA
    data <- a-b

    HTCexp(data, xgi , ygi)
})



## X Intervals
setMethod(f="x_intervals", signature(x="HTCexp"), definition=function(x){
                return(x@xgi)
        })

setReplaceMethod(f="x_intervals", signature(x="HTCexp",value="GRanges"),function(x, value) {
    x@xgi <- value
    x@intdata<-x@intdata[,id(x@xgi)]
    return (x)
})


## Y Intervals
setMethod(f="y_intervals", signature(x="HTCexp"), definition=function(x){
                return(x@ygi)
        })

setReplaceMethod(f="y_intervals", signature(x="HTCexp",value="GRanges"),function(x, value) {
    x@ygi <- value
    x@intdata<-x@intdata[id(x@ygi),]
    return (x)
})
