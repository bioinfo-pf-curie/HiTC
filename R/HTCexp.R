################
## New S4 class for representing 5C or HiC data
## This class is mainly based on the genomeIntervals class
## The interaction matrix is added to the object 
################

## class definition
setClass("HTCexp",
         representation = representation(
         "intdata"="matrix",
         "xgi"="Genome_intervals",
         "ygi"="Genome_intervals"),
         validity = function(object) {
             fails <- character(0)
              if(class(object@intdata)!="matrix") {
                  fails <- c(fails, "Expected matrix for the data slot!\n")
              }
              if(class(object@ygi)!="Genome_intervals" || class(object@xgi)!="Genome_intervals") {
                  fails <- c(fails, "Intervals 'xgi' and 'ygi' objects must be two GenomeIntervals object!\n")
              }
              if (dim(object@intdata)[1] != dim(object@ygi)[1]){
                  fails <- c(fails, "The 'ygi' Genome_Intervals object must have exactly the same length as the number of rows of the data matrix")
              }
             if (dim(object@intdata)[2] != dim(object@xgi)[1]){
                   fails <- c(fails, "The 'xgi' Genome_Intervals object must have exactly the same length as the number of column of the data matrix")
              }
             if (length(setdiff(rownames(object@intdata), as.character(object@ygi$id)))>0){
                 fails <- c(fails, "The 'ygi' intervals from the interaction matrix are different from those defined in the Genome_Intervals object")
             }
             if (length(setdiff(colnames(object@intdata), as.character(object@xgi$id)))>0){
                 fails <- c(fails, "The 'xgi' intervals from the interaction matrix are different from those defined in the Genome_Intervals object")
             }
             if (length(levels(as.factor(chromosome(object@xgi)))) > 1 || length(levels(as.factor(chromosome(object@ygi)))) > 1){
                 fails <- c(fails, "Multiple chromosome found in xgi or ygi object")
             }
             if (is.null(colnames(object@intdata)) || is.null(object@xgi$id)){
                 fails <- c(fails,"No ids for 'xgi' intervals")
             }
             if (is.null(rownames(object@intdata)) || is.null(object@ygi$id)){
                 fails <- c(fails,"No ids for 'ygi' intervals")
             }
             if (length(fails) > 0) return(fails)
             return(TRUE)
         }
         )


## constructor
HTCexp <- function(intdata, xgi, ygi){
    stopifnot(is.matrix(intdata))
    stopifnot(inherits(ygi,"Genome_intervals"))
    stopifnot(inherits(xgi,"Genome_intervals"))

    ## sort xgi and ygi data
    xgis <- sort(xgi)
    ygis <- sort(ygi)

    if (is.null(colnames(intdata)) || is.null(xgi$id)){
        stop("No ids for 'xgi' intervals")
    }
    if (is.null(rownames(intdata)) || is.null(ygi$id)){
        stop("No ids for 'ygi' intervals")
    }
    intdata <- intdata[as.character(ygi$id), as.character(xgi$id)]

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

extractRegion <- function(x, chr, from, to, exact=FALSE){
    stopifnot(inherits(x, "HTCexp"))
    fromto <- new("Genome_intervals", matrix(c(from,to), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=data.frame(seq_name=chr, inter_base=FALSE))

    ## Deal with partial overlap
    xr <- fracOverlap(x@xgi, fromto, min.frac=0, both=FALSE)$Index1
    yr <- fracOverlap(x@ygi, fromto, min.frac=0, both=FALSE)$Index1
    
    if (length(xr)>0)
        xgi <- x@xgi[xr,]
    else
        xgi <- x@xgi
    if (length(yr)>0)
        ygi <- x@ygi[yr,]
    else
        ygi <- x@ygi
    
    data <-  x@intdata[as.character(ygi$id), as.character(xgi$id)]
    
    if (exact){
        ## if no interval with from, add NA values
        ## else force overlapping coordinates
        if (length(xr)>0){
            if (min(xgi)>from){
                annot <- data.frame(seq_name=chr, inter_base=FALSE, id=paste(chr,":",from,"-",x@xgi[min(xr),1]-1, sep=""))
                                    
                if (is.element("score",colnames(annotation(xgi))))
                    annot$score <- 0
                if (is.element("strand",colnames(annotation(xgi))))
                    annot$strand <- "+"
                    
                xgi <- c(new("Genome_intervals", matrix(c(from,x@xgi[min(xr),1]-1), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=annot), xgi)
                data <- cbind(rep(NA, nrow(data)), data)
            }else if (min(xgi)<from){
                maxind <- max(which(xgi[,1]<from & xgi[,2]>from))
                xgi[maxind,1] <- from
                data <- data[,maxind:nrow(xgi)]
                xgi <- xgi[maxind:nrow(xgi),]
            }
            
            if (max(xgi)<to){
                annot <- data.frame(seq_name=chr, inter_base=FALSE, id=paste(chr,":",x@xgi[max(xr),2]+1,"-",to, sep=""))

                if (is.element("score",colnames(annotation(xgi))))
                    annot$score <- 0
                if (is.element("strand",colnames(annotation(xgi))))
                    annot$strand <- 0
              
                xgi <- c(xgi, new("Genome_intervals", matrix(c(x@xgi[max(xr),2]+1, to), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=annot))
                data <- cbind(data, rep(NA, nrow(data)))
            }else if (max(xgi)>to){
                maxind <- min(which(xgi[,1]<to & xgi[,2]>to))
                xgi[maxind,2] <- to
                data <- data[,1:maxind]
                xgi <- xgi[1:maxind,]
            }
        }

        if (length(yr)>0){
            if (min(ygi)>from){
                annot <- data.frame(seq_name=chr, inter_base=FALSE, id=paste(chr,":",from,"-",x@ygi[min(yr),1]-1, sep=""))
                
                if (is.element("score",colnames(annotation(ygi))))
                    annot$score <- 0
                if (is.element("strand",colnames(annotation(ygi))))
                    annot$strand <- "+"

                ygi <- c(new("Genome_intervals", matrix(c(from,x@ygi[min(yr),1]-1), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=annot),ygi)
                data <- rbind(rep(NA, ncol(data)), data)
            }else if (min(ygi)<from){
                maxind <- max(which(ygi[,1]<from & ygi[,2]>from))
                ygi[maxind,1] <- from
                data <- data[maxind:nrow(ygi),]
                ygi <- ygi[maxind:nrow(ygi),]
            }
            
            if (max(ygi)<to){
                annot <- data.frame(seq_name=chr, inter_base=FALSE, id=paste(chr,":",x@ygi[max(yr),2]+1,"-",to, sep=""))
                
                if (is.element("score",colnames(annotation(ygi))))
                    annot$score <- 0
                if (is.element("strand",colnames(annotation(ygi))))
                    annot$strand <- "+"
                    
                ygi <- c(ygi,new("Genome_intervals", matrix(c(x@ygi[max(yr),2]+1, to), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=annot))
                data <- rbind(data, rep(NA, ncol(data)))
            }else if (max(ygi)>to){
                maxind <- min(which(ygi[,1]<to & ygi[,2]>to))
                ygi[maxind,2] <- to
                data <- data[1:maxind,]
                ygi <- ygi[1:maxind,]
            }
        }
        colnames(data) <- as.character(id(xgi))
        rownames(data) <- as.character(id(ygi))
    }
    new("HTCexp", data, xgi, ygi)
}


################
##
## Methods
##
################


## Detail method
setMethod("detail",signature(x="HTCexp"), function(x){
    cat("HTC object\n")
    if (length(seq_name(x))==1){
        cat("Focus on genomic region [", as.character(seq_name(x)), ":",
            min(c(x_intervals(x)[,1:2],y_intervals(x)[,1:2])),"-",max(c(x_intervals(x)[,1:2],y_intervals(x)[,1:2])),"]\n", sep="")
    }else{
        cat("Focus on genomic regions [",
            unique(as.character(seq_name(y_intervals(x)))), ": ",min(y_intervals(x)[,1:2]),"-",max(y_intervals(x)[,1:2]),"] vs [",
            unique(as.character(seq_name(x_intervals(x)))), ":",min(x_intervals(x)[,1:2]),"-",max(x_intervals(x)[,1:2]),"]\n", sep="")
    }
    cat("Matrix of Interaction data: [", dim(intdata(x))[1],"-",dim(intdata(x))[2], "]\n", sep="")
    if (isBinned(x)){
        cat("Binned data - window size =", y_intervals(x)[1,2]-y_intervals(x)[1,1]+1,"\n")
        cat(nrow(y_intervals(x)),"genome intervals\n")
    }
    else{
        cat(nrow(x_intervals(x))," genome intervals from 'xgi' object\n")
        cat(nrow(y_intervals(x))," genome intervals from 'ygi' object\n")
    }
    data <- as.vector(intdata(x))
    cat("Total Reads =",sum(data, na.rm=TRUE),"\n")
    cat("Number of Interactions = ",length(data[which(data>0)]),"\n")
    cat("Median Frequency = ",median(data[which(data>0)]),"\n")
    invisible(NULL)
}) 


## Divide method
setMethod("divide", signature=c("HTCexp","HTCexp"), definition=function(x,y){
    ygi <- x@ygi[which(chromosome(x@ygi)==chromosome(y@ygi) & x@ygi[,1]==y@ygi[,1] & x@ygi[,2]==y@ygi[,2]),]
    xgi <- x@xgi[which(chromosome(x@xgi)==chromosome(y@xgi) & x@xgi[,1]==y@xgi[,1] & x@xgi[,2]==y@xgi[,2]),]
 
    a <- x@intdata[as.character(annotation(ygi)$id), as.character(annotation(xgi)$id)]
    b <- y@intdata[as.character(annotation(ygi)$id), as.character(annotation(xgi)$id)]

    a[which(a==0 | b==0)]<-NA
    b[which(a==0 | b==0)]<-NA
    #a[which(a==0)]<-b[which(a==0)]
    #b[which(b==0)]<-a[which(b==0)]
 
    data <- a/b
    new("HTCexp", data, ygi , xgi)
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
    ## Symetrical and full overlap
    if (length(setdiff(x@xgi,x@ygi))==0 && dim(interval_union(x@xgi))[1] == 1 && length(unique(x@xgi[-nrow(x@xgi),1]-x@xgi[-1,1]))==1){
        return(TRUE)
    }else{
        return(FALSE)
    }
})

## Inter/Intra chromosomal interaction
setMethod("isIntraChrom", signature(x="HTCexp"), function(x){
    stopifnot(validObject(x))
    if (levels(as.factor(chromosome(x_intervals(x)))) == levels(as.factor(chromosome(y_intervals(x))))){
        return(TRUE)
    }else{
        return(FALSE)
    }
})

## Intervals names
setMethod(f="id", signature(object="Genome_intervals"), definition=function(object){
    object$id
})

setReplaceMethod(f="id", signature(object="Genome_intervals", value="factor"),function(object, value) {
    annotation(object)$id <- value
    return (object)
})


## Plot method
setMethod("plot", signature="HTCexp", definition=function(x, ...) {
    mapC(x, ...)
})

setMethod("plot", signature=c("HTCexp","HTCexp"), definition=function(x, y, ...) {
    mapC(x, y, ...)
})

setMethod(f="range", signature(x="HTCexp"), definition=function(x){
    return(c(min(x_intervals(x), y_intervals(x)), max(x_intervals(x), y_intervals(x))))
})


## Seq_name
setMethod(f="seq_name", signature(x="HTCexp"), definition=function(x){
    seqname <- unique(c(as.character(unique(seq_name(x@xgi))), as.character(unique(seq_name(x@ygi)))))
    ##stopifnot(length(seqname)==1)
    factor(seqname)
    })

## Show method
setMethod("show",signature="HTCexp", function(object){
    cat("HTC object\n")
    if (length(seq_name(object))==1){
        cat("Focus on genomic region [", as.character(seq_name(object)), ":",
            min(c(x_intervals(object)[,1:2],y_intervals(object)[,1:2])),"-",max(c(x_intervals(object)[,1:2],y_intervals(object)[,1:2])),"]\n", sep="")
    }else{
        cat("Focus on genomic regions [",
            unique(as.character(seq_name(y_intervals(object)))), ":",min(y_intervals(object)[,1:2]),"-",max(y_intervals(object)[,1:2]),"] vs [",
            unique(as.character(seq_name(x_intervals(object)))), ":",min(x_intervals(object)[,1:2]),"-",max(x_intervals(object)[,1:2]),"]\n", sep="")
    }
    if(isIntraChrom(object)){
        cat("CIS Interaction Map\n")
    }else{
        cat("TRANS Interaction Map\n")
    }
    cat("Matrix of Interaction data: [", dim(intdata(object))[1],"-",dim(intdata(object))[2], "]\n", sep="")
    if (isBinned(object)){
        cat("Binned data - window size =", y_intervals(object)[1,2]-y_intervals(object)[1,1]+1,"\n")
        cat(nrow(y_intervals(object)),"genome intervals\n")
    }
    else{
        cat(nrow(y_intervals(object))," genome intervals from ",unique(as.character(seq_name(y_intervals(object)))), " ('ygi' object)\n")
        cat(nrow(x_intervals(object))," genome intervals from ",unique(as.character(seq_name(x_intervals(object)))), " ('xgi' object)\n")
     }
    invisible(NULL)
})


## Substraction
setMethod("substract", signature(x="HTCexp",y="HTCexp"), definition=function(x,y){
    ygi <- x@ygi[which(chromosome(x@ygi)==chromosome(y@ygi) & x@ygi[,1]==y@ygi[,1] & x@ygi[,2]==y@ygi[,2]),]
    xgi <- x@xgi[which(chromosome(x@xgi)==chromosome(y@xgi) & x@xgi[,1]==y@xgi[,1] & x@xgi[,2]==y@xgi[,2]),]
    
    a <- x@intdata[as.character(annotation(ygi)$id), as.character(annotation(xgi)$id)]
    b <- y@intdata[as.character(annotation(ygi)$id), as.character(annotation(xgi)$id)]
  
    a[which(a==0 | b==0)]<-NA
    b[which(a==0 | b==0)]<-NA

    data <- a-b
    new("HTCexp", data, ygi , xgi)
})



## X Intervals
setMethod(f="x_intervals", signature(x="HTCexp"), definition=function(x){
                return(x@xgi)
        })

setReplaceMethod(f="x_intervals", signature(x="HTCexp",value="Genome_intervals"),function(x, value) {
    x@xgi <- value
    x@intdata<-x@intdata[,as.character(x@xgi$id)]
    return (x)
})


## Y Intervals
setMethod(f="y_intervals", signature(x="HTCexp"), definition=function(x){
                return(x@ygi)
        })

setReplaceMethod(f="y_intervals", signature(x="HTCexp",value="Genome_intervals"),function(x, value) {
    x@ygi <- value
    x@intdata<-x@intdata[as.character(x@ygi$id),]
    return (x)
})
