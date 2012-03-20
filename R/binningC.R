###################################
## binning.5C
##
## Windowing of 5C data
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


binningC <- function(x, binsize=100000, bin.adjust=TRUE, upa=TRUE, method="median", use.zero=TRUE, step=1, bnorm=FALSE){
    
    stopifnot(inherits(x,"HTCexp"))

    met.agglo <- c("mean", "median", "sum")
    method <- met.agglo[pmatch(method, met.agglo)]
    if (is.na(method)) 
        stop("Error :Unknown method.")
  
    ygi <- y_intervals(x)
    xgi <- x_intervals(x)
    mat.data <- intdata(x)
    
    mn <- min(c(ygi[,1],xgi[,1]))
    mx <- max(c(ygi[,2],xgi[,2]))
    nb.bin <- floor((mx - mn)/binsize)

    if (bin.adjust){
        size.bin <- floor(binsize+((mx-mn)-(nb.bin*binsize))/nb.bin)
    }else {
        size.bin=binsize
    }
    pas <- seq(from=mn, to=mx, by=floor(size.bin/step))
    message("Bin size=",floor(size.bin/step)*step," [",step,"x",floor(size.bin/step),"]", sep="")
    nb.bin <- length(pas)-1

    if (pas[length(pas)]<mx){
        pas<-c(pas, mx)
        nb.bin <- nb.bin+1
    }

    ##Genome Intervals classes
    if (upa){
        ygi[,1] <-ygi[,2] <- round(apply(ygi[,c(1,2)],1, mean))
        xgi[,1] <-xgi[,2] <- round(apply(xgi[,c(1,2)],1, mean))
    }
    se.bin <- matrix(NA, ncol=2, nrow=nb.bin, byrow=TRUE)
    se.bin[,1] <- pas[-length(pas)]
    se.bin[,2] <- c(pas[-c(1:step)],rep(pas[length(pas)], step-1))
  
    bin.set <- new("Genome_intervals",se.bin, closed=c(TRUE,TRUE), annotation=data.frame(seq_name=rep(unique(seq_name(ygi)),nb.bin), inter_base=FALSE, strand="+", id=paste(unique(seq_name(ygi)),":",se.bin[,1],"-",se.bin[,2], sep="")))
    bin.over.FOR <- interval_overlap(bin.set, ygi)
    bin.over.REV <- interval_overlap(bin.set, xgi)
    mat.bin <- matrix(NA, ncol=nb.bin, nrow=nb.bin)
    colnames(mat.bin) <- rownames(mat.bin) <- paste(unique(seq_name(ygi)),":",se.bin[,1],"-",se.bin[,2], sep="")
    
    for (i in 1:(nb.bin)){
        fA <-bin.over.FOR[[i]]
        rA <-bin.over.REV[[i]]
        for (j in 1:(nb.bin)){
            fB <-bin.over.FOR[[j]]
            rB <-bin.over.REV[[j]]
            
            if ((length(fA>0) || length(rA)>0) && (length(fB>0) || length(rB)>0)){
                if (method=="sum"){
                    mat.bin[i,j] <- mat.bin[j,i] <- sum(mat.data[fA,rB], na.rm=TRUE)+sum(mat.data[fB,rA], na.rm=TRUE)
                    if (i==j){
                        mat.bin[i,j] <- mat.bin[j,i] <- mat.bin[i,j]/2
                    }
                }
                else if(method=="mean"){
                    sdata <- c(mat.data[fA,rB], mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <- mat.bin[j,i] <- mean(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[i,j] <- mat.bin[j,i] <- mean(sdata, na.rm=TRUE)
                    }
                }
                else if (method=="median"){
                    sdata <- c(mat.data[fA,rB], mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <- mat.bin[j,i] <-  median(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else{
                        mat.bin[i,j] <- mat.bin[j,i] <-  median(sdata, na.rm=TRUE)
                    }
                }
                
                ##Normalisation
                if (bnorm){
                    mat.bin[i,j] <- mat.bin[j,i] <- mat.bin[i,j]*log2(length(c(mat.data[fA,rB], mat.data[fB,rA])))
                }
                
            }
        }
    }
    return(HTCexp(mat.bin, bin.set, bin.set))
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

setIntervalScale <- function(x, xgi, ygi, upa=TRUE, method="mean", use.zero=TRUE){
    
    stopifnot(inherits(x,"HTCexp"))

    met.agglo <- c("mean", "median", "sum")
    method <- met.agglo[pmatch(method, met.agglo)]
    if (is.na(method)) 
        stop("Error :Unknown method.")
  
    x.ygi <- y_intervals(x)
    x.xgi <- x_intervals(x)
    mat.data <- intdata(x)
    
    ##Genome Intervals classes
    if (upa){
        x.ygi[,1] <- x.ygi[,2] <- round(apply(x.ygi[,c(1,2)],1, mean))
        x.xgi[,1] <- x.xgi[,2] <- round(apply(x.xgi[,c(1,2)],1, mean))
    }

    x.nb.bin <- dim(xgi)[1]
    y.nb.bin <- dim(ygi)[1]
    
    bin.over.y <- interval_overlap(ygi, x.ygi)
    bin.over.x <- interval_overlap(xgi, x.xgi)
  
    mat.bin <- matrix(NA, ncol=x.nb.bin, nrow=y.nb.bin)
    colnames(mat.bin) <- id(xgi)
    rownames(mat.bin) <- id(ygi)
     
    for (i in 1:(x.nb.bin)){
        rA <-bin.over.x[[i]]
        for (j in 1:(y.nb.bin)){
            fB <-bin.over.y[[j]]
            if (length(rA)>0 && length(fB)>0){
                if (method=="sum"){
                    mat.bin[j,i] <- sum(mat.data[fB,rA], na.rm=TRUE)
                } else if(method=="mean"){
                    sdata <- c(mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[j,i] <- mean(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[j,i] <- mean(sdata, na.rm=TRUE)
                    }
                }
                else if (method=="median"){
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
    return(HTCexp(mat.bin, xgi, ygi))
}

