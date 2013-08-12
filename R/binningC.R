##################################
## binningC
##
## Windowing of 'C' data
##
## x = an object of class HTCexp
## binsize = the window size
## bin.adjust = logical, the window size is adjusted taking into account the genome size
## upa = unique primer assignement. If true, a primer is assign to only one window
## method = the method used to summarize the window information (mean, median, sum)
## use.zero = use zero value for the summarization
## step = overlap of window. step = 1 for no overlap window, step = 2 50% overlap, etc
## bnorm = normalize by the log2 of number of primers in the window
###################################


binningC <- function(x, binsize=100000, bin.adjust=TRUE, upa=TRUE, method=c("median", "mean", "sum"), use.zero=TRUE, step=1, optimize.by=c("speed","memory")){

    stopifnot(inherits(x,"HTCexp"))

    #met.agglo <- c("mean", "median", "sum")
    #method <- met.agglo[pmatch(method, met.agglo)]
    met <- match.arg(method)
    optim <- match.arg(optimize.by)
        
    #if (is.na(method)) 
    #    stop("Error :Unknown method.")
  
    ygi <- y_intervals(x)
    xgi <- x_intervals(x)

    ## Optimization
    mat.data <- intdata(x)
    if (optim=="speed")
        mat.data <- as.matrix(mat.data)

    xmin <- start(range(xgi))
    xmax <- end(range(xgi))
    ymin <- start(range(ygi))
    ymax <- end(range(ygi))

    ## For cis data - set the same ranges
    if (isIntraChrom(x)){
        xmin <- ymin <- start(range(c(xgi,ygi), ignore.strand=TRUE))
        xmax <- ymax <- end(range(c(xgi,ygi), ignore.strand=TRUE))
    }
       
    x.nb.bin <- floor((xmax - xmin)/binsize)
    y.nb.bin <- floor((ymax - ymin)/binsize)
     
    if (bin.adjust){
        x.size.bin <- ceiling(binsize+((xmax-xmin+1)-(x.nb.bin*binsize))/x.nb.bin)
        y.size.bin <- ceiling(binsize+((ymax-ymin+1)-(y.nb.bin*binsize))/y.nb.bin)
    } else{ 
        x.size.bin=binsize
        y.size.bin=binsize
    }

    x.pas <- seq(from=xmin, to=xmax, by=floor(x.size.bin/step))
    x.pas <- x.pas[1:(length(x.pas)-(step-1))]
    y.pas <- seq(from=ymin, to=ymax, by=floor(y.size.bin/step))
    y.pas <- y.pas[1:(length(y.pas)-(step-1))]

    message("Bin size 'xgi' =",floor(x.size.bin/step)*step," [",step,"x",floor(x.size.bin/step),"]", sep="")
    message("Bin size 'ygi' =",floor(y.size.bin/step)*step," [",step,"x",floor(y.size.bin/step),"]", sep="")

    x.nb.bin <- length(x.pas)
    y.nb.bin <- length(y.pas)
    
    ## Unique Primer Assignment
    if (upa){
        start(ranges(xgi)) <- end(ranges(xgi)) <- mean(ranges(xgi))
        start(ranges(ygi)) <- end(ranges(ygi)) <- mean(ranges(ygi))
    }

    x.se.bin <- matrix(NA, ncol=2, nrow=x.nb.bin, byrow=TRUE)
    x.se.bin[,1] <- x.pas
    x.se.bin[,2] <- ifelse(x.pas+x.size.bin>xmax,xmax,x.pas+x.size.bin)
        
    y.se.bin <- matrix(NA, ncol=2, nrow=y.nb.bin, byrow=TRUE)
    y.se.bin[,1] <- y.pas
    y.se.bin[,2] <- ifelse(y.pas+y.size.bin>ymax,ymax,y.pas+y.size.bin)
 
    x.bin.set <- GRanges(seqnames=seqlevels(xgi), ranges = IRanges(start=x.se.bin[,1], end=x.se.bin[,2], names=paste(seqlevels(xgi),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")))
    y.bin.set <- GRanges(seqnames=seqlevels(ygi), ranges = IRanges(start=y.se.bin[,1], end=y.se.bin[,2], names=paste(seqlevels(ygi),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")))

    ## Overlap with both xgi and ygi (for 5C)
    xx.bin.over <- as.list(findOverlaps(x.bin.set, xgi))
    xy.bin.over <- as.list(findOverlaps(x.bin.set, ygi))
    yy.bin.over <- as.list(findOverlaps(y.bin.set, ygi))
    yx.bin.over <- as.list(findOverlaps(y.bin.set, xgi))

    if (optim=="speed")
        mat.bin <- matrix(0, ncol=x.nb.bin, nrow=y.nb.bin)
    else
        mat.bin <- Matrix(0, ncol=x.nb.bin, nrow=y.nb.bin)

    colnames(mat.bin) <- paste(seqlevels(xgi),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")
    rownames(mat.bin) <- paste(seqlevels(ygi),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")

    for (i in 1:(y.nb.bin)){
        fA <-yy.bin.over[[i]]
        rA <-yx.bin.over[[i]]
        
        for (j in 1:(x.nb.bin)){
            fB <-xy.bin.over[[j]]
            rB <-xx.bin.over[[j]]

            if ((length(fA>0) || length(rA)>0) && (length(fB>0) || length(rB)>0)){
                if (met=="sum"){
                    mat.bin[i,j] <- sum(mat.data[fA,rB], na.rm=TRUE)+sum(mat.data[fB,rA], na.rm=TRUE)
                    if (rownames(mat.bin)[i]==colnames(mat.bin)[j]){
                        mat.bin[i,j] <- mat.bin[i,j]/2
                    }
                }
                else if(met=="mean"){
                    sdata <- c(as.vector(mat.data[fA,rB]), as.vector(mat.data[fB,rA]))
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <- mean(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[i,j] <- mean(sdata, na.rm=TRUE)
                    }
                }
                else if (met=="median"){
                    sdata <- c(as.vector(mat.data[fA,rB]), as.vector(mat.data[fB,rA]))
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <-  median(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else{
                        mat.bin[i,j] <-  median(sdata, na.rm=TRUE)
                    }
                    a <- median(sdata, na.rm=TRUE)
                }
            }
        }
    }

    if (optim=="speed")
        mat.bin <- as(mat.bin, "Matrix")
 
    return(HTCexp(mat.bin, x.bin.set, y.bin.set))  
}

###################################
## setIntervalScale
##
## Force xgi and ygi intervals of a HTCexp object
##
## x = HTCexp object
## xgi = xgi Genome_interval object to use to define the HTCexp object
## ygi = ygi Genome_interval object to use to define the HTCexp object
## upa = unique primer assignement. If true, a primer is assign to only one window
## method = the method used to summarize the new interval information (mean, median, sum)
## use.zero = use zero value for the summarization
##
##################################

setIntervalScale <- function(x, xgi, ygi, upa=TRUE, method=c("median","mean","sum"), use.zero=TRUE, optimize.by=c("speed","memory")){
    
    stopifnot(inherits(x,"HTCexp"))

    #met.agglo <- c("mean", "median", "sum")
    #method <- met.agglo[pmatch(method, met.agglo)]
    #if (is.na(method)) 
    #    stop("Error :Unknown method.")
    met <- match.arg(method)
    optim <- match.arg(optimize.by)

    x.ygi <- y_intervals(x)
    x.xgi <- x_intervals(x)

    mat.data <- intdata(x)
    if (optim=="speed")
       mat.data <- as.matrix(mat.data)
    
    ## Unique Primer Assignment
    if (upa){
        start(ranges(x.xgi)) <- end(ranges(x.xgi)) <- mean(ranges(x.xgi))
        start(ranges(x.ygi)) <- end(ranges(x.ygi)) <- mean(ranges(x.ygi))
    }

    x.nb.bin <- length(xgi)
    y.nb.bin <- length(ygi)
    
    bin.over.y <- as.list(findOverlaps(ygi, x.ygi))
    bin.over.x <- as.list(findOverlaps(xgi, x.xgi))
  
    if (optim=="speed")
        mat.bin <- matrix(0, ncol=x.nb.bin, nrow=y.nb.bin)
    else
        mat.bin <- Matrix(0, ncol=x.nb.bin, nrow=y.nb.bin)

    colnames(mat.bin) <- id(xgi)
    rownames(mat.bin) <- id(ygi)
 
    for (i in 1:(x.nb.bin)){
        rA <-bin.over.x[[i]]
        for (j in 1:(y.nb.bin)){
            fB <-bin.over.y[[j]]
            if (length(rA)>0 && length(fB)>0){
                if (met=="sum"){
                    mat.bin[j,i] <- sum(mat.data[fB,rA], na.rm=TRUE)
                } else if(met=="mean"){
                    sdata <- c(mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[j,i] <- mean(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[j,i] <- mean(sdata, na.rm=TRUE)
                    }
                }
                else if (met=="median"){
                    sdata <- c(mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[j,i] <- median(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[j,i] <- median(sdata, na.rm=TRUE)
                    }
                }
            }
        }
    }

    if (optim=="speed")
        mat.bin <- as(mat.bin,"Matrix")
    
    return(HTCexp(mat.bin, xgi, ygi))
}

