

## Normalized per reads number
setMethod("normPerReads", signature=c("HTCexp"), definition=function(x){
    x@intdata <- x@intdata/sum(x@intdata,na.rm=TRUE)
    x
})

## Normalized per expected number of count
setMethod("normPerExpected", signature=c("HTCexp"), definition=function(x, stdev=FALSE, ...){
    expCounts <- getExpectedCounts(x, stdev=stdev, ...)
    if (stdev){
        x@intdata <- (x@intdata-expCounts$exp.interaction)/expCounts$stdev.estimate
    }else{
        x@intdata <- x@intdata/expCounts$exp.interaction
    }
    x
})


## Normalized using TRANS data
setMethod("normPerTrans", signature=c("HTCexp","HTCexp","HTCexp"), definition=function(x, xtrans, ytrans, method="sum"){

    ## Match the good objects
    ## TODO check whether the two objects are really the same
    ## x and trans share the same x_intervals
    if(length(intersect(id(x_intervals(x)),id(x_intervals(xtrans)))) != length(id(x_intervals(x)))){
   	stop("No match between xgi cis and trans objects")
    }   
    
    ## x and trans share the same y_intervals
    if(length(intersect(id(y_intervals(x)),id(y_intervals(ytrans)))) != length(id(y_intervals(x)))){
 	stop("No match between ygi cis and trans objects")
    }   

    ix <- intdata(x)
    iy <- intdata(ytrans)
    iz <- intdata(xtrans)

    ## Weigth matrices from trans data
    normfacRow<-apply(iy, 1, "mean")
    wY <- matrix(rep(normfacRow,ncol(ix)), ncol=ncol(ix), byrow=FALSE)
    
    normfacCol<-apply(iz, 2, "mean")
    wZ <- matrix(rep(normfacCol,nrow(ix)), nrow=nrow(ix), byrow=TRUE)

    ## Normalize cis data
    ## Method
    met.agglo <- c("mult", "sum", "mean")
    method <- met.agglo[pmatch(method, met.agglo)]
    if (is.na(method)) 
        stop("Error :Unknown method ['mult','sum','mean' are expected].")

    if (method == "sum"){
        xnorm <- ix/(wY+wZ)
    }else if (method == "mult"){
        xnorm <- ix/(wY*wZ)
    }else if (method == "mean"){
        wM <- matrix(NA, ncol=ncol(wY), nrow=nrow(wY))
        for (i in 1:nrow(wY)){wM[i,]<-pmax(wY[i,], wZ[i,])}
        xnorm <- ix/(wM)
    }
    
    intdata(x) <- xnorm
    return(x)
})

###################################
## getExpectedCounts
##
## Estimate the expected interaction counts from a HTCexp object based on the interaction distances
##
## x = a HTCexp object
## span=fraction of the data used for smoothing at each x point. The larger the f value, the smoother the fit. 
## bin=interval size (in units corresponding to x). If lowess estimates at two x values within delta of one another, it fits any points between them by linear interpolation. The default is 1% of the range of x. If delta=0 all but identical x values are estimated independently.
## stdev = logical,the standard deviation is estimated for each interpolation point
## plot = logical, display lowess and variance smoothing
##
## bin is used to speed up computation: instead of computing the local polynomial fit at each data point it is not computed for points within delta of the last computed point, and linear interpolation is used to fill in the fitted values for the skipped points. 
## This function may be slow for large numbers of points. Increasing bin should speed things up, as will decreasing span. 
## Lowess uses robust locally linear fits. A window, dependent on f, is placed about each x value; points that are inside the window are weighted so that nearby points get the most weight. 
##
## NOTES
## All variances are calculated (even identical values) because of the parallel implementation.
## Cannot use rle object because the neighboring have to be calculated on the wall dataset. Or have to find a way to convert rle index to real index ...
## Easy to do with a 'for' loop but the parallel advantages are much bigger


###################################

getExpectedCounts<- function(x, span=0.01, bin=0.005, stdev=FALSE, plot=FALSE){
    stopifnot(inherits(x,"HTCexp"))
    if ("package:parallel" %in% search()){
        lFun <- mclapply
    } else {
        lFun <- lapply
    }
    
    ydata <- as.vector(intdata(x))
    ydata[which(is.na(ydata))] <- 0
    xdata.dist <- as.vector(intervalsDist(x))
    o<- order(xdata.dist)
    xdata.dist <- xdata.dist[o]
    ydata <- ydata[o]
    
    delta <- bin*diff(range(xdata.dist))
    ######################
    ## Lowess Fit
    ######################
    message("Lowess fit ...")
    #lowess.fit <- .C("lowess", x = as.double(xdata.dist), as.double(ydata), 
    #                length(ydata), as.double(span), as.integer(3), as.double(delta), 
    #                y = double(length(ydata)), double(length(ydata)), double(length(ydata)), PACKAGE = "stats")$y

    lowess.fit <-lowess(x=xdata.dist, y=ydata, f=span, delta=delta)$y
      
    y1 <- sort(ydata)
    y1 <-  quantile(y1[which(y1>1)], probs=0.99)

    if (plot){
        par(font.lab=2, mar=c(4,4,1,1))
        plot(x=xdata.dist, y=ydata,  xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="5C counts", main="", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, col="gray", frame=FALSE)
        points(x=xdata.dist[order(lowess.fit)], y=sort(lowess.fit), type="l", col="red")
    }
    lowess.mat <- matrix(lowess.fit[order(o)], nrow=length(y_intervals(x)), byrow=FALSE)
    rownames(lowess.mat) <- id(y_intervals(x))
    colnames(lowess.mat) <- id(x_intervals(x))

    ######################
    ## Variance estimation
    ######################
    stdev.mat <- NULL
    if (stdev){
        message("Standard deviation calculation ...")
        ##interpolation
        ind <- getDeltaRange(delta, xdata.dist)
        lx <- length(xdata.dist)
        Q <- floor(lx*span)
        stdev.delta <- unlist(lFun(1:length(ind), function(k){
            i <- ind[k]
            x1 <- xdata.dist[i]
            
            ## Neighbors selection 2*Q
            ll <- i-Q-1
            lr <- i+Q-1
            if (ll<0) ll=0
            if (lr>lx) lr=lx
            xdata.dist.sub <- xdata.dist[ll:lr]
            ydata.sub <- ydata[ll:lr]
            ## Select the Q closest distances
            d <- abs(x1-xdata.dist.sub)
            o2 <- order(d)[1:Q]
            x2 <- xdata.dist.sub[o2]
            y2 <- ydata.sub[o2]
            ## Distance between x and other points
            dref <- d[o2] 
            drefs <- dref/max(abs(dref-x1)) ##max(dref) - NS
            ## Tricube weigths and stdev calculation
            w <- tricube(drefs)
            sqrt <- w*(y2-lowess.fit[i])^2
            
            stdev <- sqrt(sum(sqrt)/
                          (((length(sqrt)-1) * sum(w))/length(sqrt)))
        }))

        if (plot){
            points(x=xdata.dist[ind], y=lowess.fit[ind], col="black", cex=.8, pch="+")
            legend(x="topright", lty=c(1,NA), pch=c(NA,"+"), col=c("red","black"),legend=c("Lowess fit","Interpolation points"), cex=.8, bty="n")
        }
        
        ## Approximation according to delta
        stdev.estimate <- approx(x=xdata.dist[ind], y=stdev.delta, method="linear", xout=xdata.dist)$y
        stdev.mat <- matrix(stdev.estimate[order(o)], nrow=length(y_intervals(x)), byrow=FALSE)
        rownames(stdev.mat) <- id(y_intervals(x))
        colnames(stdev.mat) <- id(x_intervals(x))
    }    
    return(list(exp.interaction=lowess.mat,stdev.estimate=stdev.mat))
}
    
###################################
## getDeltaRange
## INTERNAL FUNCTION
## Calculate the interpolation points from the delta value
##
## delta = lowess delta parameter
## xdata = Intervals distances matrix
###################################

getDeltaRange <- function(delta, xdata){
    message("Delta=",delta)
    
    if (delta>0){
        ind <- 1
        for (i in 1:length(xdata)-1){
            if (xdata[i+1]>=delta){
                ind <- c(ind,i)
                delta=delta+delta
            }
        }
        if (max(ind)<length(xdata)){
            ind <- c(ind, length(xdata))
        }
        message("Calculating stdev for",xdata[ind],"bps")
    }else{
        ind <- 1:length(xdata)
    }
    ind
}

###################################
## tricube
## INTERNAL FUNCTION
## tricube distance weigth
##
## x = numeric. A distance
###################################

tricube <- function(x) {
    ifelse (abs(x) < 1, (1 - (abs(x))^3)^3, 0)
}
