###################################
## primersDist
##
## Calculate the distances between primers/intervals for cis data only
##
## x = an object of class HTCexp.
##
###################################

intervalsDist<-function(x){

    stopifnot(inherits(x,"HTCexp"))
    if (! isIntraChrom(x)){
       stop("The distance can only be calculated for intrachromosomal distances")
    }

    dist.mat <- apply(x_intervals(x), 1, function(xgi){
        apply(y_intervals(x), 1, function(ygi,xgi){
            ## Check if overlap
            if ((xgi[2]>=ygi[1] && xgi[2]<=ygi[2]) || (ygi[2]>=xgi[1] && ygi[2]<=xgi[2])){
                z <- abs(mean(ygi)-mean(xgi))
            }else{
                ##abs(ygi[2]-xgi[2])
                z <- min(abs(ygi[2]-xgi[1]),abs(ygi[1]- xgi[2]))
            }
            return(z)
        }, xgi=xgi)
    })
        
    rownames(dist.mat) <- id(y_intervals(x))
    colnames(dist.mat) <- id(x_intervals(x))

    dist.mat
}

###################################
## CQC
##
## 'C' Quality Control
##
## x = an object (or list) of class HTCexp. In case of list, the x are merged as one.
## trans.ratio = if true, plot the histogram of inter/intrachromosomal interactions 
## hist.interac = if true, plot the interaction frequency. How many interaction have n reads 
## scat.interac.dist = if true, plot the scatter plot of counts vs genomic distances. Interaction Distance vs interaction frequency
## hist.dist = if true, plot an histogram of distances between y and x intervals. How many interaction have n distance
## trim.range = remove the extreme values by trimming the counts. Only use for plotting functions
## dev.new = if true, draw each plot in a new view
##
###################################


CQC <- function(x, cis.trans.ratio=TRUE, hist.interac=TRUE, scat.interac.dist=TRUE, hist.dist=TRUE, trim.range=0.98, dev.new=FALSE){

    if (is.list(x)){
        ## Separate Intra/Inter chromosomal interaction
        x.intra <- x[unlist(lapply(x,isIntraChrom))]
        x.inter <- x[!unlist(lapply(x,isIntraChrom))]
        
        xdata.intra <- unlist(lapply(x.intra,function(x){
            stopifnot(inherits(x,"HTCexp"))
            return(as.vector(intdata(x)))
        }))
        xdata.intra.dist <- unlist(lapply(x.intra,function(x){
            stopifnot(inherits(x,"HTCexp"))
            return(as.vector(intervalsDist(x)))}))
        
        xdata.inter <- unlist(lapply(x.inter,function(x){
            stopifnot(inherits(x,"HTCexp"))
            return(as.vector(intdata(x)))
        }))
    } else if(inherits(x,"HTCexp")){
        x.intra <- x.inter <- NULL
        xdata.intra <- xdata.inter <- NULL
        if (isIntraChrom(x)){
            x.intra <- x
            xdata.intra <- as.vector(intdata(x.intra))
            xdata.intra.dist <- as.vector(intervalsDist(x.intra))
        }else {
            x.inter <- x
            xdata.inter <- as.vector(intdata(x.inter))
        }
    }else{
        stop("Wrong input type. 'HTCexp' or list of 'HTCexp' objects expected.")
    }

    ## Remove zero-values
    xdata.intra.dist <- xdata.intra.dist[xdata.intra>0]
    xdata.intra <- xdata.intra[xdata.intra>0]
    xdata.inter <- xdata.inter[xdata.inter>0]
    xdata <- c(xdata.inter, xdata.intra)

    ## Basic Descriptive Statistics
    dstat <- matrix(0,ncol=4,nrow=3)
    colnames(dstat) <- c("nbreads","nbinteraction","averagefreq","medfreq")
    rownames(dstat) <- c("all","cis","trans")

    dstat["all",]<-c(sum(xdata), sum(length(xdata)), round(mean(xdata),3), median(xdata))
    dstat["cis",]<-c(sum(xdata.intra), sum(length(xdata.intra)),  round(mean(xdata.intra),3), median(xdata.intra))
    if (!is.null(xdata.inter))
        dstat["trans",]<-c(sum(xdata.inter), sum(length(xdata.inter)),  round(mean(xdata.inter),3), median(xdata.inter))
    
    nbplot <- length(which(c(cis.trans.ratio, scat.interac.dist, hist.dist, hist.interac)))
    if (length(x.inter)>0){
       nbplot <- nbplot+1
    }
    if (!dev.new){
             par(mfrow=c(ceiling(nbplot/2),2), mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
     }

    ## Cis/Trans ratio
    if (cis.trans.ratio && length(xdata)>0){
        if (dev.new){
            dev.new()
            par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
        }
        barplot(c(sum(xdata.intra)/sum(xdata), sum(xdata.inter)/sum(xdata)), ylab="Fraction of Reads", names.arg=c("CIS","TRANS"),main="Intra/Inter Chromosomal Interaction", col=c("#FBB4AE","#B3CDE3"), cex.lab=0.7,cex.axis=0.7, cex.main=0.9)
    }

    ## Scatter plot of interaction distance for intrachromosomal interactions
    if(scat.interac.dist && length(x.intra)>0){
        if (dev.new){
            dev.new()
            par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
        }
        y1 <- sort(xdata.intra)
        y1 <-  quantile(y1[which(y1>1)], probs=trim.range)
        plot(x=xdata.intra.dist, y=xdata.intra,  xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="Interaction Counts", main="Scatter Plot (Frequency(Y) vs Distance(X))\nCIS Interaction Counts", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, cex.main=0.9, frame=FALSE)
    }

    ## Histogram of interaction counts for intra and/or interchromosomal interaction
    if (hist.interac){
        if (length(x.intra)>0){
            if (dev.new){
                dev.new()
                par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
            }
            x <- xdata.intra
            x <- x[which(x<quantile(x, probs=trim.range))]
            histdens<-sort(hist(x,breaks=500, plot=FALSE)$density)
            ymax <- histdens[floor(length(histdens)*.995)]
            hist(x,freq=FALSE,breaks=500,col="grey", ylim=c(0,ymax),main="Interaction Frequency Histogram\nCIS Interaction Counts",ylab="Probability Density",xlab="Interaction Frequency",cex.lab=0.7, cex.axis=0.7, pch=20, cex.main=0.9)
            lines(density(x), col="red")
        }
       if (length(x.inter)>0){
           if (dev.new){
               dev.new()
               par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
           }
           x <- xdata.inter
           x <- x[which(x<quantile(x, probs=trim.range))]
            histdens<-sort(hist(x,breaks=500, plot=FALSE)$density)
            ymax <- histdens[floor(length(histdens)*.995)]
            hist(x,freq=FALSE,breaks=500,col="grey", ylim=c(0,ymax),main="Interaction Frequency Histogram\nTRANS Interaction Counts",ylab="Probability Density",xlab="Interaction Frequency",cex.lab=0.7, cex.axis=0.7, pch=20, cex.main=0.9)
            lines(density(x), col="red")
        }
    }

    ## Histogram of Distance for intrachromosomal interactions
    if (hist.dist && length(x.intra)>0){
        if (dev.new){
             dev.new()
             par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
         }
        x <- xdata.intra.dist
        x <- xdata.intra.dist[which(x<quantile(x, probs=trim.range))]
        histdens<-sort(hist(x,breaks=500, plot=FALSE)$density)
        ymax <- histdens[floor(length(histdens)*.995)]
        hist(x,freq=FALSE,breaks=500,col="grey", ylim=c(0,ymax),main="Histogram of Distances Between \nX & Y intervals",ylab="Probability Density",xlab="CIS Interaction Distance",cex.lab=0.7, cex.axis=0.7, pch=20, cex.main=0.9)
        lines(density(x), col="red")
    }


    
    
  return(dstat)
}
