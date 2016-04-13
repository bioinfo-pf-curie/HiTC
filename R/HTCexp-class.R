################
## New S4 class to represent 5C or HiC data
## This class is mainly based on the GRanges class
## The interaction matrix is added to the object 
################

## class definition
setClass("HTCexp",
         representation = representation(
         "intdata"="Matrix",
         "xgi"="GRanges",
         "ygi"="GRanges")
         )


setValidity("HTCexp",
            function(object){
                fails <- character(0L)
                
                ## Format
                if(!inherits(object@intdata,"Matrix")) {
                    fails <- c(fails, "Expected Matrix for the data slot!\n")
                }
                if(!inherits(object@ygi,"GRanges") || !inherits(object@xgi,"GRanges")) {
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
                if (length(object@xgi) == 0L || length(object@ygi) == 0L){
                    fails <- c(fails, "Inervals xgi or ygi are of size 0")
                }
                if (length(seqlevels(object@xgi)) > 1L || length(seqlevels(object@ygi)) > 1L){
                    fails <- c(fails, "Multiple chromosome found in xgi or ygi object")
                }
                if (length(setdiff(rownames(object@intdata), id(object@ygi)))>0L){
                  fails <- c(fails, "The 'ygi' intervals from the interaction matrix are different from those defined in the GRanges object")
                }
                if (length(setdiff(colnames(object@intdata), id(object@xgi)))>0L){
                    fails <- c(fails, "The 'xgi' intervals from the interaction matrix are different from those defined in the GRanges object")
                }
                
                if (length(fails) > 0L) return(fails)
                return(TRUE)
            }
)

## constructor
HTCexp <- function(intdata, xgi, ygi, lazyload=FALSE)
{
    ## check xgi/ygi length
    stopifnot(inherits(ygi,"GRanges"))
    stopifnot(inherits(xgi,"GRanges"))
    
    if (length(xgi)==0L || length(ygi)==0L)
        stop("Cannot create HTCexp object of size 0")
    
    if(length(ygi) != nrow(intdata))
        stop("Length 'ygi' is not equal to number of rows in 'intdata'")
    
    if(length(xgi) != ncol(intdata))
        stop("Length 'xgi' is not equal to number of columns in 'intdata'")

    
    ## sort xgi and ygi data
    xgi <- sort(xgi)
    ygi <- sort(ygi)
  
    ## Check matrix format
    if (!inherits(intdata, "Matrix")){
        intdata <- as(as.matrix(intdata), "Matrix")
    }

    ## Check matrix names
    colnames(intdata) <- id(xgi)
    rownames(intdata) <- id(ygi)

    ## Create HTCexp
    obj <- new("HTCexp", intdata, xgi, ygi)

    
    if (lazyload & seqlevels(xgi)==seqlevels(ygi)){
        obj <- forceSymmetric(obj)
        obj <- forceTriangular(obj)
    }
    obj
}


################
##
## Functions
##
################
removeIntervals <- function(x, ids)
{
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

extractRegion <- function(x, MARGIN=c(1,2), chr, from, to, exact=FALSE)
{
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
                    edata$score <- 0L
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
                    edata$score <- 0L
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
                    edata$score <- 0L
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
                    edata$score <- 0L
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

setMethod("c", "HTCexp", function(x, ...){
    if (missing(x))
        args <- unname(list(...))
    else
        args <- unname(list(x, ...))
    HTClist(args)
})


setMethod("detail",signature(x="HTCexp"),
          function(x){
              stopifnot(validObject(x))
              cat("HTC object\n")
              
              r <- range(x)
              if (length(r)==1){
                  cat("Focus on genomic region [", seqlevels(r), ":",
                      start(r),"-",end(r),"]\n", sep="")
              }else{
                  cat("Focus on genomic regions [",
                      seqlevels(r)[1], ":", start(r)[1], "-", end(r)[1],"] vs [",
                      seqlevels(r)[2], ":", start(r)[2], "-", end(r)[2],"]\n", sep="")
              }
              
              if(isIntraChrom(x)){
                  cat("CIS Interaction Map\n")
              }else{
                  cat("TRANS Interaction Map\n")
              }
              cat("Matrix of Interaction data: [", dim(x@intdata)[1],"-",dim(x@intdata)[2], "]\n", sep="")
              if (isBinned(x)){
                  cat("Binned data - window size =", width(x_intervals(x))[1],"\n")
                  cat(length(x_intervals(x)),"genome intervals\n")
              }
              else{
                  cat(length(x_intervals(x))," genomic ranges from 'xgi' object\n")
                  cat(length(y_intervals(x))," genomic ranges from 'ygi' object\n")
              }
	      
           
              xdata <- x@intdata
              
              cat("Total Reads = ",sum(xdata, na.rm=TRUE),"\n")
              cat("Number of Interactions = ",nnzero(xdata),"\n")
              cat("Median Frequency = ",median(xdata@x[xdata@x>0L]),"\n")
              #cat("Sparsity = ",round((length(xdata)-nnzero(xdata))/length(xdata),3),"\n")
              cat("Sparsity = ",round(length(xdata@x)/length(xdata),3),"\n")

	      invisible(NULL)
          }
) 

## Divide method
setMethod("divide", signature=c("HTCexp","HTCexp"),
          function(x,y){

              xgi.ov <- findOverlaps(x@xgi, y@xgi, type="equal")
              ygi.ov <- findOverlaps(x@ygi, y@ygi, type="equal")

              ##xgi <- subsetByOverlaps(x@xgi, y@xgi, type="equal")
              ##ygi <- subsetByOverlaps(x@ygi, y@ygi, type="equal")
              
              if (length(queryHits(xgi.ov)) == 0L || length(queryHits(ygi.ov)) == 0L){
                  stop("No equal intervals found in x and y objects")
              }

              x.xgi <- x@xgi[queryHits(xgi.ov)]
              x.ygi <- x@ygi[queryHits(ygi.ov)]
              y.xgi <- y@xgi[subjectHits(xgi.ov)]
              y.ygi <- y@ygi[subjectHits(ygi.ov)]

              a <- x@intdata[id(x.ygi), id(x.xgi)]
              b <- y@intdata[id(y.ygi), id(y.xgi)]
              ##a <- x@intdata[id(ygi), id(xgi)]
              ##b <- y@intdata[id(ygi), id(xgi)]
              a[which(a==0L | b==0L)]<-NA
              b[which(a==0L | b==0L)]<-NA

              ## Use x as reference for intervals description
              colnames(b) <- colnames(a)
              rownames(b) <- rownames(a)
              data <- a/b
              
              HTCexp(data, x.xgi , x.ygi)
          }
)

## forceSymmetric
setMethod(f="forceSymmetric", signature(x="HTCexp", uplo="missing"),
          function(x, uplo){forceSymmetric(x, uplo="U")})

setMethod(f="forceSymmetric", signature(x="HTCexp", uplo="character"),
          function(x, uplo){
              stopifnot(uplo=="U" || uplo=="L")# || uplo=="S")
              idata <- x@intdata
              #if (uplo == "U" || uplo == "L"){
              #  idata <- forceSymmetric(idata, uplo)
              #}else if (uplo== "S"){
              #   idata <- idata + t(idata)
              #   diag(idata) <- diag(idata)/2
              #   idata <- as(idata, "symmetricMatrix")
              #}
              #intdata(x) <- idata
              intdata(x) <- forceSymmetric(idata, uplo)
              x
          }
)

setMethod(f="forceTriangular", signature(x="HTCexp"),
          function(x){
            stopifnot(isSymmetric(x))
            idata <- intdata(x)
            idata[lower.tri(idata)] <- 0
            intdata(x)<-as(idata, "sparseMatrix")
            x
})

#setMethod(f="forcePairwise", signature(x="HTCexp"),
#          function(x){
#            idata <- intdata(x)
#            idata <- idata + t(idata)
#            diag(idata) <- diag(idata)/2
#            idata <- as(idata, "sparseMatrix")
#            intdata(x) <- idata
#            x
#          })

## Interaction Data
setMethod(f="intdata", signature(x="HTCexp"),
          function(x){
            if (isTriangular(x) && isIntraChrom(x)){
              x <- forceSymmetric(x)
            }
             x@intdata
          }
)

setReplaceMethod(f="intdata", signature(x="HTCexp",value="Matrix"),
                 function(x, value){
                     x@intdata <- value
                     stopifnot(validObject(x))
                     return(x)
                 }
)

## Symetric interaction maps
setMethod("isSymmetric", signature(object="HTCexp"),
	  function(object){
              ret <- FALSE
              idata <- object@intdata
              if (Matrix::isSymmetric(idata))
                ret <- TRUE
              else if (isTriangular(idata))
                ret <- TRUE
              else if (seqlevels(object@xgi)==seqlevels(object@ygi)){
                if (ncol(idata)==nrow(idata))
                  ret <- suppressWarnings(all(idata-t(idata)==0, na.rm=TRUE))
              }
              ret
          }
)

## is Triangular Matrix
setMethod("isTriangular", signature(object="HTCexp"),
	  function(object){
            object@intdata[which(is.na(object@intdata))] <- 0
            ret <- Matrix::isTriangular(object@intdata)
          }
)



## Initialize method
setMethod("initialize",signature=c("HTCexp"),
          function(.Object, interaction.data, xgi, ygi){
              if (!missing(interaction.data) && !missing(ygi) && !missing(xgi)){
                  .Object@intdata<-interaction.data
                  .Object@ygi <-ygi
                  .Object@xgi <-xgi
              }
              return(.Object)
          }
)

## Binned data
setMethod("isBinned", signature(x="HTCexp"),
          function(x){
              stopifnot(validObject(x))
              ret <- TRUE
              
              ## Same bins for intra chromosomal
              ##if(isIntraChrom(x) && !all.equal(x@xgi, x@ygi)){
              if (isIntraChrom(x) && sum(countMatches(x@xgi, x@ygi)) != length(x@xgi)){
                  ret <- FALSE
              }
              ## Same width (exept for last bins - checked in 95% of bins)
              x.excl <- round(length(x@xgi)*.95)
              y.excl <- round(length(x@ygi)*.95)
              if (length(duplicated(width(x@xgi)))< x.excl || length(duplicated(width(x@ygi)))< y.excl)
                  ret <- FALSE
              
              ## no gaps
              ##if (isDisjoint(x@xgi) || isDisjoint(x@ygi))
              if (length(reduce(x@xgi)) != 1 || length(reduce(x@ygi)) != 1)
                ret <- FALSE

              ret
          }
)

## Inter/Intra chromosomal interaction
setMethod("isIntraChrom", signature(x="HTCexp"),
          function(x){
              stopifnot(validObject(x))
              ret <- FALSE
              if (seqlevels(x@xgi) == seqlevels(x@ygi))
                  ret <- TRUE
              ret
          }
)

## Intervals names
setMethod(f="id", signature(x="GRanges"),
          function(x){
              ret <- NULL
              if (!is.null(names(x)))
                  ret <- names(x)
              else if (!is.null(elementMetadata(x)$name))
                  ret <- elementMetadata(x)$name
              as.character(ret)
          }
)

## Range
setMethod(f="range", signature(x="HTCexp"),
          function(x){
              suppressWarnings(range(c(x@ygi, x@xgi), ignore.strand=TRUE))
          }
)

## Seqlevels
setMethod(f="seqlevels", signature(x="HTCexp"),
          function(x){
              unique(c(seqlevels(x@ygi), seqlevels(x@xgi)))
          }
)

## Show method
setMethod("show",signature="HTCexp",
          function(object){
              stopifnot(validObject(object))
              cat("HTCexp object\n")
          
              r <- range(object)
              if (length(r)==1){
                  cat("Focus on genomic region [", seqlevels(r), ":",
                      start(r),"-",end(r),"]\n", sep="")
              }else{
                  cat("Focus on genomic regions [",
                      seqlevels(r)[1], ":", start(r)[1], "-", end(r)[1],"] vs [",
                      seqlevels(r)[2], ":", start(r)[2], "-", end(r)[2],"]\n", sep="")
              }
              cat("Matrix of Interaction data: [", dim(object@intdata)[1],"-",dim(object@intdata)[2], "]\n", sep="")
              gf <- union(names(elementMetadata(object@xgi)), names(elementMetadata(object@ygi)))
              if (length(gf)>0)
                  cat("Genomic Features : ", paste(gf, collapse="-"),"\n")
              invisible(NULL)
          }
)


## Substraction
setMethod("substract", signature(x="HTCexp",y="HTCexp"),
          function(x,y){
              xov <- findOverlaps(x@xgi, y@xgi, type="equal")
              yov <- findOverlaps(x@ygi, y@ygi, type="equal")

              a <- x@intdata[id(x@ygi)[queryHits(yov)],
                        id(x@xgi)[queryHits(xov)]]
              b <- y@intdata[id(y@ygi)[subjectHits(yov)],
                        id(y@xgi)[subjectHits(xov)]]
              
              #a <- x@intdata[id(x@ygi),id(x@xgi)]
              #b <- y@intdata[id(y@ygi),id(y@xgi)]
              #a[which(a==0L | b==0L)]<-NA
              #b[which(a==0L | b==0L)]<-NA
              data <- a-b
              xgi <- x@xgi[queryHits(xov)]
              ygi <- x@ygi[queryHits(yov)]
              
              HTCexp(data, xgi , ygi)
          }
)

## Summary
setMethod("summary", signature=c(object="HTCexp"),
          function(object){
              interactors <- seqlevels(object)
              if (isIntraChrom(object))
                  interactors <- rep(seqlevels(object),2)
              else
                  interactors <- seqlevels(object)

              if(max(object@intdata, na.rm=TRUE)==0){
                  dstat <- c(interactors, 0, 0, 0, 0, 1)                 
              }else{
                  ## Force cohersion to dgTMatrix to have @x equal of non zero values
                  xdata <- as(as(object@intdata, "sparseMatrix"),"dgTMatrix")
                  mx <- round(mean(xdata@x, na.rm=TRUE),4)
                  mdx <- round(median(xdata@x, na.rm=TRUE),4)
                  dstat <- c(interactors, sum(xdata@x, na.rm=TRUE), nnzero(xdata, na.counted=TRUE), mx, mdx, round((length(xdata)-length(xdata@x))/length(xdata),4))
              }
              names(dstat) <- c("seq1", "seq2", "nbreads","nbinteraction","averagefreq","medfreq","sparsity")  
              data.frame(summary=dstat, stringsAsFactors=FALSE)
          })



## X Intervals
setMethod(f="x_intervals", signature(x="HTCexp"),
          function(x){x@xgi})

setReplaceMethod(f="x_intervals", signature(x="HTCexp",value="GRanges"),
                 function(x, value){
                     x@xgi <- value
                     intdata <- as(as.matrix(x@intdata[,as.vector(id(x@xgi))]), "Matrix")
                     colnames(intdata) <-as.vector(id(x@xgi))
                     rownames(intdata) <- rownames(x@intdata)
                     x@intdata <- intdata
                     x
                 }
)


## Y Intervals
setMethod(f="y_intervals", signature(x="HTCexp"),
          function(x){x@ygi})

setReplaceMethod(f="y_intervals", signature(x="HTCexp",value="GRanges"),
                 function(x, value){
                     x@ygi <- value
                     intdata <- as(as.matrix(x@intdata[as.vector(id(x@ygi)),]), "Matrix")
                     rownames(intdata) <-as.vector(id(x@ygi))
                     colnames(intdata) <- colnames(x@intdata)
                      x@intdata <- intdata
                     x
                 }
)


## X and Y Intervals
setMethod(f="xy_intervals", signature(x="HTCexp"),
          function(x){
              GRangesList(xgi=x@xgi, ygi=x@ygi)
          }
)
