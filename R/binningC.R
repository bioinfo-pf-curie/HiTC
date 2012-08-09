###################################
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


binningC <- function(x, binsize=100000, bin.adjust=TRUE, upa=TRUE, method="median", use.zero=TRUE, step=1){
    
    stopifnot(inherits(x,"HTCexp"))

    met.agglo <- c("mean", "median", "sum")
    method <- met.agglo[pmatch(method, met.agglo)]
    if (is.na(method)) 
        stop("Error :Unknown method.")
  
    ygi <- y_intervals(x)
    xgi <- x_intervals(x)
    mat.data <- intdata(x)

    xmin=min(xgi[,1])
    xmax=max(xgi[,2])
    ymin=min(ygi[,1])
    ymax=max(ygi[,2])

    ## For cis data - set the same ranges
    if (isIntraChrom(x)){
        xmin=min(c(ygi[,1],xgi[,1]))
        ymin=min(c(ygi[,1],xgi[,1]))
        xmax=max(c(xgi[,2], ygi[,2]))
        ymax=max(c(xgi[,2], ygi[,2]))
    }
       
    x.nb.bin <- floor((xmax - xmin)/binsize)
    y.nb.bin <- floor((ymax - ymin)/binsize)
     
    if (bin.adjust){
        x.size.bin <- floor(binsize+((xmax-xmin)-(x.nb.bin*binsize))/x.nb.bin)
        y.size.bin <- floor(binsize+((ymax-ymin)-(y.nb.bin*binsize))/y.nb.bin)
    } else{ 
        x.size.bin=binsize
        y.size.bin=binsize
    }

    x.pas <- seq(from=xmin, to=xmax, by=floor(x.size.bin/step))
    y.pas <- seq(from=ymin, to=ymax, by=floor(y.size.bin/step))
    message("Bin size 'xgi' =",floor(x.size.bin/step)*step," [",step,"x",floor(x.size.bin/step),"]", sep="")
    message("Bin size 'ygi' =",floor(y.size.bin/step)*step," [",step,"x",floor(y.size.bin/step),"]", sep="")

    
    x.nb.bin <- length(x.pas)-1
    y.nb.bin <- length(y.pas)-1
    
    if (x.pas[length(x.pas)]<xmax){
        x.pas<-c(x.pas, xmax)
        x.nb.bin <- x.nb.bin+1
    }
    if (y.pas[length(y.pas)]<ymax){
        y.pas<-c(y.pas, ymax)
        y.nb.bin <- y.nb.bin+1
    }

    ## Genome Intervals classes
    if (upa){
        ygi[,1] <-ygi[,2] <- round(apply(ygi[,c(1,2)],1, mean))
        xgi[,1] <-xgi[,2] <- round(apply(xgi[,c(1,2)],1, mean))
    }

    x.se.bin <- matrix(NA, ncol=2, nrow=x.nb.bin, byrow=TRUE)
    x.se.bin[,1] <- x.pas[-length(x.pas)]
    x.se.bin[,2] <- c(x.pas[-c(1:step)],rep(x.pas[length(x.pas)], step-1))
        
    y.se.bin <- matrix(NA, ncol=2, nrow=y.nb.bin, byrow=TRUE)
    y.se.bin[,1] <- y.pas[-length(y.pas)]
    y.se.bin[,2] <- c(y.pas[-c(1:step)],rep(y.pas[length(y.pas)], step-1))

    x.bin.set <- new("Genome_intervals",x.se.bin, closed=c(TRUE,TRUE), annotation=data.frame(seq_name=rep(unique(seq_name(xgi)),x.nb.bin), inter_base=FALSE, strand="+", id=paste(unique(seq_name(xgi)),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")))
    y.bin.set <- new("Genome_intervals",y.se.bin, closed=c(TRUE,TRUE), annotation=data.frame(seq_name=rep(unique(seq_name(ygi)),y.nb.bin), inter_base=FALSE, strand="+", id=paste(unique(seq_name(ygi)),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")))

    ## Overlap with both xgi and ygi (for 5C)
    xx.bin.over <- interval_overlap(x.bin.set, xgi)
    xy.bin.over <- interval_overlap(x.bin.set, ygi)
    yy.bin.over <- interval_overlap(y.bin.set, ygi)
    yx.bin.over <- interval_overlap(y.bin.set, xgi)

    
    mat.bin <- matrix(NA, ncol=x.nb.bin, nrow=y.nb.bin)
    colnames(mat.bin) <- paste(unique(seq_name(xgi)),":",x.se.bin[,1],"-",x.se.bin[,2], sep="")
    rownames(mat.bin) <- paste(unique(seq_name(ygi)),":",y.se.bin[,1],"-",y.se.bin[,2], sep="")

    for (i in 1:(y.nb.bin)){
        fA <-yy.bin.over[[i]]
        rA <-yx.bin.over[[i]]
        
        for (j in 1:(x.nb.bin)){
            fB <-xy.bin.over[[j]]
            rB <-xx.bin.over[[j]]

            #print ("--------------------")
            #print (bin.set[i,])
            #print ("---------fa--------")
            #print (id(ygi)[fA])
            #print ("---------ra--------")
            #print (id(xgi)[rA])
            #print ("--------------------")
            #print (bin.set[j,])
            #print ("---------fb--------")
            #print (id(xgi)[fB])
            #print ("---------rb--------")
            #print (id(ygi)[rB])
            
            if ((length(fA>0) || length(rA)>0) && (length(fB>0) || length(rB)>0)){
                if (method=="sum"){
                    mat.bin[i,j] <- sum(mat.data[fA,rB], na.rm=TRUE)+sum(mat.data[fB,rA], na.rm=TRUE)
                    if (rownames(mat.bin)[i]==colnames(mat.bin)[j]){
                        mat.bin[i,j] <- mat.bin[i,j]/2
                    }
                }
                else if(method=="mean"){
                    sdata <- c(mat.data[fA,rB], mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <- mean(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else {
                        mat.bin[i,j] <- mean(sdata, na.rm=TRUE)
                    }
                }
                else if (method=="median"){
                    sdata <- c(mat.data[fA,rB], mat.data[fB,rA])
                    if (!use.zero && length(sdata[sdata!=0]>0)){
                        mat.bin[i,j] <-  median(sdata[sdata!=0], na.rm=TRUE)
                    }
                    else{
                        mat.bin[i,j] <-  median(sdata, na.rm=TRUE)
                    }
                }
            }
        }
    }   
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

