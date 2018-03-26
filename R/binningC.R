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
###################################

binningC <- function(x, binsize=100000, bin.adjust=TRUE, upa=TRUE, method=c("sum", "median", "mean"), use.zero=TRUE, step=1, optimize.by=c("speed","memory")){

    stopifnot(inherits(x,"HTCexp"))
    met <- match.arg(method)
    optim <- match.arg(optimize.by)        
  
    ygi <- y_intervals(x)
    xgi <- x_intervals(x)

    ## Optimization
    mat.data <- intdata(x)
    if (optim=="speed")
        mat.data <- as.matrix(mat.data)

    rxi <- range(xgi, ignore.strand=TRUE)
    ryi <- range(ygi, ignore.strand=TRUE)

    ## For cis data - set the same ranges
    if (isIntraChrom(x)){
        xmin <- ymin <- start(range(c(xgi,ygi), ignore.strand=TRUE))
        xmax <- ymax <- end(range(c(xgi,ygi), ignore.strand=TRUE))
    }else{
        xmin <- start(rxi)
        xmax <- end(rxi)
        ymin <- start(ryi)
        ymax <- end(ryi)
    }

    if (max(xmax, ymax)<binsize){
        binsize=xmax
        warning("Data smaller than binsize. Cannot be binned at such resolution.")
    }
           
    if (bin.adjust){
        x.nb.bin <- floor((xmax - xmin)/binsize)
        y.nb.bin <- floor((ymax - ymin)/binsize)
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
        start(xgi) <- end(xgi) <- start(xgi)+(width(xgi)/2)
        start(ygi) <- end(ygi) <- start(ygi)+(width(ygi)/2)
    }

    x.se.bin <- matrix(NA_integer_, ncol=2, nrow=x.nb.bin, byrow=TRUE)
    x.se.bin[,1] <- x.pas
    x.se.bin[,2] <- ifelse(x.pas+x.size.bin>xmax,xmax,x.pas+x.size.bin)
    x.se.bin[,2] <- x.se.bin[,2] - 1
    
    y.se.bin <- matrix(NA_integer_, ncol=2, nrow=y.nb.bin, byrow=TRUE)
    y.se.bin[,1] <- y.pas
    y.se.bin[,2] <- ifelse(y.pas+y.size.bin>ymax,ymax,y.pas+y.size.bin)
    y.se.bin[,2] <- y.se.bin[,2] - 1

    x.bin.set <- GRanges(seqnames=seqlevels(xgi), ranges = IRanges(start=x.se.bin[,1], end=x.se.bin[,2], names=paste(seqlevels(xgi),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")))
    y.bin.set <- GRanges(seqnames=seqlevels(ygi), ranges = IRanges(start=y.se.bin[,1], end=y.se.bin[,2], names=paste(seqlevels(ygi),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")))

    ## Binning for 5C data - Overlap with both xgi and ygi
    if (!isBinned(x)){
        xx.bin.over <- as.list(findOverlaps(x.bin.set, xgi))
        xy.bin.over <- as.list(findOverlaps(x.bin.set, ygi))
        yy.bin.over <- as.list(findOverlaps(y.bin.set, ygi))
        yx.bin.over <- as.list(findOverlaps(y.bin.set, xgi))
        
        ## vector of pairs to combine
        if (isIntraChrom(x)){
            p <- matrix(rep(1:y.nb.bin,2), ncol=2)
            p <- rbind(p, t(combn(1:y.nb.bin, 2)))
        }else{
            p <- matrix(c(rep(1:y.nb.bin,x.nb.bin), rep(1:x.nb.bin, each=y.nb.bin)), ncol=2)
        }
                
        out<- apply(p, 1, function(idx){
            fA <-yy.bin.over[[idx[1]]]
            rA <-yx.bin.over[[idx[1]]]
            fB <-xy.bin.over[[idx[2]]]
            rB <-xx.bin.over[[idx[2]]]
               
             if ((length(fA>0) || length(rA)>0) && (length(fB>0) || length(rB)>0)){
                 if (met=="sum"){
                     if (idx[1]==idx[2])
                         (sum(mat.data[fA,rB], na.rm=TRUE)+sum(mat.data[fB,rA], na.rm=TRUE))/2
                     else
                        sum(mat.data[fA,rB], na.rm=TRUE)+sum(mat.data[fB,rA], na.rm=TRUE)
                 }
                else if(met=="mean"){
                    sdata <- c(as.vector(mat.data[fA,rB]), as.vector(mat.data[fB,rA]))
                    if (!use.zero && length(sdata[sdata!=0]>0))
                        mean(sdata[sdata!=0], na.rm=TRUE)
                    else
                        mean(sdata, na.rm=TRUE)
                }
                else if (met=="median"){
                    sdata <- c(as.vector(mat.data[fA,rB]), as.vector(mat.data[fB,rA]))
                    if (!use.zero && length(sdata[sdata!=0]>0))
                        median(sdata[sdata!=0], na.rm=TRUE)
                    else
                        median(sdata, na.rm=TRUE)
                }
            }else{
                return(0)
            }
        })
        out <- unlist(out)
        ## Create new HTCexp object
        mat.bin <- Matrix::sparseMatrix(i=p[which(out!=0),1], j=p[which(out!=0),2], x=out[which(out!=0)], dims=c(y.nb.bin, x.nb.bin))
        
        if (isIntraChrom(x))
            mat.bin <- forceSymmetric(mat.bin)
        
    }else{
        ## Binning for Hi-C data
        if (isSymmetric(x)){
            xx.bin.over <- as.list(findOverlaps(x.bin.set, xgi))
            xy.bin.over <- xx.bin.over
            
            ## vector of pairs to combine
            if (y.nb.bin>1){
                p <- matrix(rep(1:y.nb.bin,2), ncol=2)
                p <- rbind(p, t(combn(1:y.nb.bin, 2)))
            }else{
                p <- matrix(c(1,1), ncol=2, byrow=2)
            }
        }else{
            xx.bin.over <- as.list(findOverlaps(x.bin.set, xgi))
            xy.bin.over <- as.list(findOverlaps(y.bin.set, ygi))
            
            ## vector of pairs to combine
            p <- matrix(c(rep(1:y.nb.bin,x.nb.bin), rep(1:x.nb.bin, each=y.nb.bin)), ncol=2)
        }
        ## sum of pairs of overlap
        out<- apply(p, 1, function(idx){
            i <-xy.bin.over[[idx[1]]]
            j <-xx.bin.over[[idx[2]]]

            if (met=="sum"){
                if (idx[1]==idx[2])
                    sum(mat.data[i,j], na.rm=TRUE) ##/2
                else
                    sum(mat.data[i,j], na.rm=TRUE)
            }
            else if(met=="mean"){
                sdata <- as.vector(mat.data[i,j])
                if (!use.zero && length(sdata[sdata!=0]>0))
                    mean(sdata[sdata!=0], na.rm=TRUE)
                else
                    mean(sdata, na.rm=TRUE)
            }
            else if (met=="median"){
                sdata <- as.vector(mat.data[i,j])
                if (!use.zero && length(sdata[sdata!=0]>0))
                    median(sdata[sdata!=0], na.rm=TRUE)
                else
                    median(sdata, na.rm=TRUE)
            }
        })
        ## Create new HTCexp object
        mat.bin <- Matrix::sparseMatrix(i=p[which(out!=0),1], j=p[which(out!=0),2], x=out[which(out!=0)], dims=c(y.nb.bin, x.nb.bin))

        if (isSymmetric(x))
            mat.bin <- forceSymmetric(mat.bin)
    }
    colnames(mat.bin) <- paste(seqlevels(xgi),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")
    rownames(mat.bin) <- paste(seqlevels(ygi),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")
  
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

setIntervalScale <- function(x, xgi, ygi, upa=TRUE, method=c("sum", "median", "mean"), use.zero=TRUE, optimize.by=c("speed","memory")){
    
    stopifnot(inherits(x,"HTCexp"))
    met <- match.arg(method)
    optim <- match.arg(optimize.by)

    x.ygi <- y_intervals(x)
    x.xgi <- x_intervals(x)

    mat.data <- intdata(x)
    if (optim=="speed")
       mat.data <- as.matrix(mat.data)
    
    ## Unique Primer Assignment
    if (upa){
        start(ranges(x.xgi)) <- end(ranges(x.xgi)) <- start(ranges(x.xgi)) + floor(width(ranges(x.xgi)))/2
        start(ranges(x.ygi)) <- end(ranges(x.ygi)) <- start(ranges(x.ygi)) + floor(width(ranges(x.ygi)))/2
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

