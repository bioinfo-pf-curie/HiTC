###################################
## primersDist
##
## Calculate the distances between primers/intervals for cis data only
##
## x = an object of class HTCexp.
##
##
###################################

intervalsDist<-function(x, use.zero=TRUE){

    stopifnot(inherits(x,"HTCexp"))
    if (! isIntraChrom(x)){
       stop("The distances can only be calculated for intrachromosomal interactions")
    }

    xgi <- x_intervals(x)
    ygi <- y_intervals(x)
    xdata <- intdata(x)
    
    dist.mat <- apply(cbind(as.matrix(ranges(xgi)), id(xgi)), 1, function(xi, yi){
        sxi <- as.numeric(xi[1])
        exi <- sxi+as.numeric(xi[2])-1
        nxi <- xi[3]
        
        mxi <- (exi-sxi)/2L
        sygi <- start(yi)
        eygi <- end(yi)

        ##over <- yi %over% IRanges(start=xi[1], end=xi[2]) ## !!
        over <- (exi>=sygi & exi<=eygi) | (sxi>=sygi & sxi<eygi)
        id.over <- which(over)
        id.nover <- which(!over)
                
        d <- rep(NA_integer_, length(ygi))
        if (length(id.nover)>0){
            aft <- abs(eygi[id.nover]-sxi)
            bef <- abs(sygi[id.nover]- exi)
            di <- bef
            di[which(bef>aft)] <- aft[which(bef>aft)]
            d[id.nover] <- di
            ##d[id.nover] <- apply(cbind(abs(eygi[id.nover]-sxi), abs(sygi[id.nover]- exi)), 1, min) ## !!
        }
        if (length(id.over)>0){
            myi <- (eygi[id.over] - sygi[id.over])/2
            d[id.over] <- abs(myi - mxi)
        }
        d
    }, yi=ranges(ygi))

    rownames(dist.mat) <- id(ygi)
    colnames(dist.mat) <- id(xgi)

    ## report only non-zero values
    if (!use.zero){
        dist.mat <- dist.mat + 1 ## zero distance are different from zero values
        sl.xdata <- xdata!=0
        dist.mat[which(!sl.xdata)] <- 0
        dist.mat <- as(dist.mat, "dgTMatrix")
    }
    dist.mat
}##intervalsDist

## Cis/Trans ratio
plotInterIntraRatio <- function(xdata.intra, xdata.inter, ...){
    sintra <- sum(xdata.intra, na.rm=TRUE)
    sinter <- sum(xdata.inter, na.rm=TRUE)
    sall <- sintra+sinter
    
    barplot(c(sintra/sall, sinter/sall), ylab="Fraction of Reads", names.arg=c("Intra","Inter"),main="Intra/Inter Chromosomal Interaction", ...)
}

## Scatter plot of interaction distance for intrachromosomal interactions
slidingWindow <- function(mydata=consV, win=c(1, 100), start=1, end=length(consV)) {
        mydata <- mydata[start:end]
        windex <- t(sapply(0:length(mydata), function(x) win+x))
        mywind <- sapply(seq(along=mydata), function(x) mean(mydata[windex[x,1]:windex[x,2]],
        na.rm=TRUE))
        mywind <- mywind[1:(length(mywind)-win[2])]
        return(mywind)
}


plotIntraDist <- function(xdata.intra, xdata.intra.dist, trim.range=.98, winsize=NA, add=FALSE, log=FALSE, fit=FALSE, fit.lim=NA, fit.out=1, ...){
    stopifnot(length(xdata.intra)==length(xdata.intra.dist))
    r <- range(xdata.intra.dist)

    ## Remove zeros
    idx <- which(xdata.intra>0)
    xdata.intra <- xdata.intra[idx]
    xdata.intra.dist <- xdata.intra.dist[idx]
    
    ## Windowing
    ## To optimize - see tapply function
    if (!is.na(winsize)){
        win <- seq.int(from=r[1], to=r[2], by=winsize)
        xp <- yp <- rep(NA, length(win)-1)
      
        tmp <- sapply(1L:(length(win)-1L), function(i){
            idx <- which(xdata.intra.dist>=win[i] & xdata.intra.dist<win[i+1])

            sub <- xdata.intra.dist[idx]
            x <- mean(xdata.intra.dist[idx], na.rm=TRUE)
            y <- mean(xdata.intra[idx], na.rm=TRUE)
            return(c(x, y))
        })
       
        xp <- tmp[1,]
        yp <- tmp[2,]
        ptype <- "b"
    }else{
        xp <- xdata.intra.dist
        yp <- xdata.intra
        ptype <- "p"
    }
    names(xp) <- paste("n",1:length(xp), sep="")
    names(yp) <- paste("n",1:length(yp), sep="")
    
    ## Trim the interaction counts
    if (trim.range<1){
        qt <- quantile(yp, probs=c((1-trim.range),trim.range), na.rm=TRUE)
        idx <- which(yp>=qt[1] & yp<=qt[2])
        yp <- yp[idx]
        xp <- xp[idx]
    }
    
    ## Log
    if (log){
        xp <- log10(xp)
        yp <- log10(yp)
    }

    if (fit){
        ## remove extreme distances
        #qt <- quantile(xp, probs=c(0.1,0.8), na.rm=TRUE)
        #idx <- which(xp>=qt[1] & xp<=qt[2])
        #xp <- xp[idx]
        #yp <- yp[idx]

        ## Define the scaling region
        if (length(fit.lim)==2){
            yp.fit <- yp[which(xp>fit.lim[1] & xp<fit.lim[2])]
            xp.fit <- xp[which(xp>fit.lim[1] & xp<fit.lim[2])]
        }else{
            xp.fit <- xp
            yp.fit <- yp
        }
        res.fit <- lm(yp.fit~xp.fit)

        ## remove outliers on the fitted region
        ## outliers are defined as the points with the higher residual values (distance to the regression curve)
        if (fit.out<1){
            r <- residuals(res.fit)
            th<-quantile(r, probs=fit.out)
            out <- names(which(r>th))
            
            xp.fit <- xp.fit[setdiff(names(xp.fit), out)]
            yp.fit <- yp.fit[setdiff(names(xp.fit), out)]
            res.fit <- lm(yp.fit~xp.fit)

            xp <- xp[setdiff(names(xp), out)]
            yp <- yp[setdiff(names(xp), out)]
        }
    }

    ## Plotting function
    if (!add){
        plot(x=c(min(xp), max(xp)), y=c(min(yp), max(yp)),  xlab="Genomic Distance (log10)", ylab="Interaction Counts (log10)",
              frame=FALSE, type="n", ...)
    }
    
    if (fit){
        pcol <- RColorBrewer::brewer.pal(8, "Pastel2")
        if (length(fit.lim)==2)
            rect(fit.lim[1], min(yp)-1, fit.lim[2], max(yp), col=pcol[5], border=pcol[5])
        abline(res.fit, ...)
        #text(x=max(xp), y=max(yp), labels=paste("a=", round(res.fit$coefficients[2],6)), font=2, cex=.7)
    }
    points(x=xp, y=yp, type=ptype, ...)
    if (fit)
        return(res.fit)
    else
        invisible(NULL)
}

plotHistCounts <- function(x, trim.range=.98, title, ...){
    if (trim.range<1)
        x <- x[which(x<quantile(x, probs=trim.range, na.rm=TRUE))]
    histdens<-sort(hist(x,breaks=500, plot=FALSE)$density)
    ymax <- histdens[floor(length(histdens)*.995)]
    hist(x,freq=FALSE,breaks=500,col="grey", ylim=c(0,ymax),main=paste("Interaction Frequency Histogram\n",title," Interaction Counts",sep=""),ylab="Probability Density",xlab="Interaction Frequency")
    lines(density(x), col="red")
}

plotHistDist <- function(xdata.intra.dist, trim.range=.98, ...){
    x <- xdata.intra.dist
    if (trim.range<1)
        xdata.intra.dist <- xdata.intra.dist[which(x<quantile(x, probs=trim.range, na.rm=TRUE))]
    
    histdens<-sort(hist(xdata.intra.dist,breaks=500, plot=FALSE)$density)
    ymax <- histdens[floor(length(histdens)*.995)]
    hist(xdata.intra.dist,freq=FALSE,breaks=500,col="grey", ylim=c(0,ymax),main="Histogram of Distances Between \nX & Y intervals",ylab="Probability Density",xlab="CIS Interaction Distance",...)
    lines(density(x), col="red")
}

extractCounts <- function(x){
  x.intra <- x.inter <- NULL
  xdata.intra <- xdata.inter <- NULL
  
  if (inherits(x,"HTClist")){
      ## Separate Intra/Inter chromosomal interaction
      if (length(which(isIntraChrom(x)))>0){
          x.intra <- x[isIntraChrom(x)]
          xdata.intra <- unlist(lapply(x.intra,function(x){
              return(as(as(intdata(x), "sparseMatrix"),"dgTMatrix")@x)
          }))
          xdata.intra.dist <- unlist(lapply(x.intra,function(x){
              return(intervalsDist(x, use.zero=FALSE)@x)}))
      }
      if (length(which(!isIntraChrom(x)))>0){
          x.inter <- x[!isIntraChrom(x)]
          xdata.inter <- unlist(lapply(x.inter,function(x){
              return(intdata(x)@x)
          }))
      }
  } else if(inherits(x,"HTCexp")){
      if (isIntraChrom(x)){
          x.intra <- x
          xdata.intra <- as(as(intdata(x.intra), "sparseMatrix"),"dgTMatrix")@x
          xdata.intra.dist <- intervalsDist(x.intra, use.zero=FALSE)@x
      }else {
          x.inter <- x
          xdata.inter <- intdata(x.inter)@x
          xdata.intra.dist <- NULL
      }
  }else{
      stop("Wrong input type. 'HTCexp' or 'HTClist' objects expected.")
  }
  return(list(intra=xdata.intra, intra.dist=xdata.intra.dist, inter=xdata.inter))
}


###################################
## CQC
##
## 'C' Quality Control
##
## x = an object of class HTCexp/HTClist. In case of list, the x are merged as one.
## trans.ratio = if true, plot the histogram of inter/intrachromosomal interactions 
## hist.interac = if true, plot the interaction frequency. How many interaction have n reads 
## scat.interac.dist = if true, plot the scatter plot of counts vs genomic distances. Interaction Distance vs interaction frequency
## hist.dist = if true, plot an histogram of distances between y and x intervals. How many interaction have n distance
## trim.range = remove the extreme values by trimming the counts. Only used for plotting functions
## dev.new = if true, draw each plot in a new view
##
###################################

CQC<- function(x, cis.trans.ratio=TRUE, hist.interac=TRUE, scat.interac.dist=TRUE, hist.dist=TRUE, trim.range=0.98, winsize=NA, dev.new=TRUE){
  
    message("Get data ...")
    data <- extractCounts(x)
    xdata.inter <- data$inter
    xdata.intra <- data$intra
    xdata <- c(xdata.inter, xdata.intra)
    xdata.intra.dist <- data$intra.dist
    
    message("Generate quality control plots ...")
    nbplot <- length(which(c(cis.trans.ratio, scat.interac.dist, hist.dist, hist.interac)))
    if (!is.null(xdata.inter)){
        nbplot <- nbplot+1
    }
    if (!dev.new)
        par(mfrow=c(ceiling(nbplot/2),2), mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)

  ## Cis/Trans ratio
  if (cis.trans.ratio && length(xdata)>0){
      if (dev.new){
          dev.new()
          par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
      }
      plotInterIntraRatio(xdata.intra, xdata.inter,
                          cex.lab=0.7, cex.axis=0.7, cex.main=0.9,
                          col=c("#FBB4AE","#B3CDE3"))
  }
   
  ## Scatter plot of interaction distance for intrachromosomal interactions
  if(scat.interac.dist && !is.null(xdata.intra)>0){
      if (dev.new){
          dev.new()
          par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
      }
      plotIntraDist(xdata.intra, xdata.intra.dist, trim.range, winsize=winsize,
                    cex=0.5, cex.lab=0.7, pch=20,
                    cex.axis=0.7, cex.main=0.9, main="Scatter Plot (Frequency(Y) vs Distance(X))\nCIS Interaction Counts")
  }
 
  ## Histogram of interaction counts for intra and/or interchromosomal interaction
  if (hist.interac){
      if (dev.new){
          dev.new()
          par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
      }
      plotHistCounts(xdata.intra, trim.range,
                     cex.lab=0.7, cex.axis=0.7,
                     pch=20, cex.main=0.9, title="CIS")
      if (!is.null(xdata.inter)){
          if (dev.new){
              dev.new()
              par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
          }
          plotHistCounts(xdata.inter, trim.range,
                         ,cex.lab=0.7, cex.axis=0.7,
                         pch=20, cex.main=0.9, title="TRANS")
      }
  }
 
  ## Histogram of Distance for intrachromosomal interactions
  if (hist.dist && !is.null(xdata.intra)){
      if (dev.new){
          dev.new()
          par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
      }
      plotHistDist(xdata.intra.dist, trim.range,
                   cex.lab=0.7, cex.axis=0.7,
                   pch=20, cex.main=0.9)
  }
}
