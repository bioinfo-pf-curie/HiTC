###################################
## colorC
## INTERNAL FUNCTION
## Define a palette for the heatmap
##
## low = color for minimal value
## high = color for extreme values (up and down)
## mid = color for medium values
## k = number of color generated
##
##################################

colorC <- function (low="white", high=c("red"), mid=NA, k = 50){
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    if (is.na(mid)) {
        r <- seq(low[1], high[1], len = k)
        g <- seq(low[2], high[2], len = k)
        b <- seq(low[3], high[3], len = k)
    }
    if (!is.na(mid)) {
        k2 <- round(k/2)
        mid <- col2rgb(mid)/255
        r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1], len = k2))
        g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2], len = k2))
        b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3], len = k2))
    }
    rgb(r, g, b)
    ##colorRampPalette(colors=c(low, mid, high))(k)
}


###################################
## plottingfunc
## INTERNAL FUNCTION
## alias to the function used to plot the data
##
## x = color for minimal value
## high = color for extreme values (up and down)
## mid = color for medium values
## k = number of color generated
##
##################################

plottingfunc <- function(x, y, z, show.zero=FALSE, col, ...){
    if (!show.zero){
        zi <- which(!is.na(z) & z == 0) ##zi <- which(Matrix:::is0(z))
        z <- as.matrix(z)
        z[zi] <- NA
        col <- col[-1]
    }else{
        z <- as.matrix(z)
    }
    image(x, y, z, useRaster = TRUE, col=col, ...)
}


###################################
## heatmapC
## INTERNAL FUNCTION
## Draw heatmap of the C data
##
## xdata = 'C' interaction map as a matrix
## names = logical, add axis with the bin/interval's names
## value = logical, add values on the matrix
## show.na = logical, na are shown in gray
## show.zero = logical, zero are plotted 
## col.neg = color used in mypalette for neg count (low, mid, high)
## col.pos = color used in mypalette for pos count (low, mid, high)
## col.na = color used for NA values
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################

heatmapC <- function(xdata,  names=FALSE, value=FALSE, show.zero=FALSE,  show.na=TRUE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", grid=FALSE, title=NULL){

    xdata <- t(xdata)
    xdata <- xdata[,ncol(xdata):1]

    xdata.pos <- xdata.neg <- NULL
   
    if (max(xdata, na.rm=TRUE)>0){
        xdata.pos <- xdata
        xdata.pos[xdata<0] <- 0
    }
    if (min(xdata, na.rm=TRUE)<0){
        xdata.neg <- xdata
        xdata.neg[xdata>=0] <- 0
    }

    ## k is the same for positive or negative values, to ensure that the same range of colors are used
    k <- length(unique(as.vector(abs(xdata))))
    
    if (!is.null(xdata.neg)){   
        col.neg <- colorC(col.neg[3],col.neg[1],mid=col.neg[2], k=k)
        plottingfunc(x=1:nrow(xdata),y=1:ncol(xdata),z=xdata.neg, show.zero=show.zero, axes=FALSE,ylab="",xlab="",col=col.neg)
        par(new = TRUE)
    }
    
    if (!is.null(xdata.pos)){   
        col.pos <- colorC(col.pos[1],col.pos[3],mid=col.pos[2],k=k)
        xdata <- as(xdata,"sparseMatrix")
        plottingfunc(x=1:nrow(xdata),y=1:ncol(xdata),z=xdata.pos,show.zero=show.zero,axes=FALSE,ylab="",xlab="",col=col.pos)
    }

    if (names){
        axis(side=2, at=1:ncol(xdata), labels=colnames(xdata), lwd=0.5, las=1, cex.axis=0.7)  
        axis(side=1, at=1:nrow(xdata), labels=rownames(xdata), lwd=0.5, las=2, cex.axis=0.7)
    }
    if (value){
        for (i in nrow(xdata):1){
            text(x=i,y=1:ncol(xdata),labels=round(xdata[i,],2), cex=0.7)
        }
    }
    
    ## ###################
    ## Heatmap options
    ## ###################
    if (show.na){
        if (length(which(is.na(xdata))>0)){
            par(new = TRUE)
            na.xdata <- matrix(NA, ncol=ncol(xdata), nrow=nrow(xdata))
            na.xdata[which(is.na(xdata))] <- 1
            plottingfunc(x=1:nrow(na.xdata),y=1:ncol(na.xdata),z=na.xdata,show.zero=show.zero,axes=FALSE,ylab="",xlab="",col=col.na,add=TRUE)
        }
    }
    
    if (grid){
        abline(h=seq(0.5,ncol(xdata)+0.5, by=1), mar=c(0,0))
        abline(v=seq(0.5,nrow(xdata)+0.5, by=1), mar=c(0,0))
    }
    
    if (!is.null(title)){
        text(x=ncol(xdata)-nchar(title), y=nrow(xdata)-5, title, col="darkgray", font=2)
    }
}

###################################
## triViewC
## INTERNAL FUNCTION
## Draw triangle view of the C data
##
## xdata = 'C' interaction map as a matrix
## value = logical, add values on the matrix
## show.na = logical, na are shown in gray
## col.low = color used in mypalette for low count (low, mid, high)
## col.high = color used in mypalette for high count (low, mid, high)
## col.na = color used for NA values
## plot.zero = plot the zero value or not
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################

triViewC <- function(xdata, flip=FALSE, value=FALSE, show.zero=FALSE, show.na=TRUE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", title=NULL){

    d <- min(dim(xdata))

    ## Start with a matrix for speed implementation - optimize.by option ???
    trimat <- matrix(NA, ncol=d*2, nrow=d)
    xdata <- as.matrix(xdata)
    for (w in 0:(d-1)){
        s <- gdiag(xdata, w=w)
        ls <- length(s)
        ss <- c(s,s)[as.vector(sapply(1:ls,function(x){return(c(x,x+ls))}))]
        trimat[w+1, (w+1):(d*2-w)] <- ss
    }
    trimat <- as(trimat, "Matrix")
    
    if (flip){
        trimat <- trimat[nrow(trimat):1,]
    }

    trimat.pos <- trimat.neg <- NULL
    if (max(trimat, na.rm=TRUE)>0){
        trimat.pos <- trimat
        trimat.pos[trimat<0] <- 0
    }
    if (min(trimat, na.rm=TRUE)<0){
        trimat.neg <- trimat
        trimat.neg[trimat>0] <- 0
    }
    
    k <- length(unique(as.vector(abs(trimat))))
    
    if (!is.null(trimat.neg)){     
        col.neg <- colorC(col.neg[3],col.neg[1],mid=col.neg[2],k=k)
        plottingfunc(y=1:nrow(trimat),x=1:ncol(trimat),z=t(trimat.neg),show.zero=show.zero,axes=FALSE,ylab="",xlab="",col=col.neg)
        par(new = TRUE)
    }

    if (!is.null(trimat.pos)){     
        col.pos <- colorC(col.pos[1],col.pos[3],mid=col.pos[2],k=k)
        plottingfunc(y=1:nrow(trimat),x=1:ncol(trimat),z=t(trimat.pos),show.zero=show.zero,axes=FALSE,ylab="",xlab="",col=col.pos)
    }
    if (show.na){
        if (length(which(is.na(trimat))>0)){
            par(new = TRUE)
            na.trimat <- matrix(NA, ncol=ncol(trimat), nrow=nrow(trimat))
            na.trimat[which(is.na(trimat))] <- 1
            plottingfunc(y=1:nrow(na.trimat),x=1:ncol(na.trimat),z=t(na.trimat),show.zero=show.zero,axes=FALSE,ylab="",xlab="",col=col.na)
        }
    }

    if (value){
        coordodd <- seq(.5,ncol(trimat), by=2)
        coord <- seq(1.5,ncol(trimat), by=2)
        for (i in 1:nrow(trimat)){
            if (i%%2)
                text(coord,y=i,labels=round(trimat[i,seq(1,ncol(trimat),2)],2), cex=.7)
            else
                text(coordodd,y=i,labels=round(trimat[i,seq(1,ncol(trimat),2)],2), cex=.7)
        }
    }
    if (!is.null(title)){
        text(x=1, y=nrow(trimat), adj=c(0,1), title, col="darkgray", font=2)
    }
}

###################################
## gdiag
## INTERNAL FUNCTION
## Extract the diagonal of a matrix - Generalisation of diag function
##
## x = a matrix
## w = step from the real diagonal. w=0 return the diagonal of the matrix.
##
##################################

gdiag <- function(x, w=0){
    if ((m <- min(dim(x))) == 0L) 
        return(vector(typeof(x), 0L))
    if (w<m){
        if (w>0)
            y <- x[1L + (0L+w):(m - 1L)* (dim(x)[1L] + 1L) - w]
        else
            y <- x[1L + 0L:(m - 1L - abs(w))* (dim(x)[1L] + 1L) + abs(w)]
        return(y)
    }else{
        return(vector(typeof(x), 0L)) 
    }
}


###################################
## setEnvDisplay
## INTERNAL FUNCTION
## Define the environment display for visualization
## Based on the graphics package
##
## x = HTCexp/HTClist object
## y = HTCexp object
## tracks = List of GRanges objects. Each object represent a genome track information
##
##################################
setEnvDisplay <- function(x, y=NULL, view, tracks=NULL){

  chrom <- seqlevels(x)
  lc <- length(chrom)

  if(!is.null(tracks)){
      stopifnot(unlist(lapply(tracks, inherits,"GRanges")))
      ntrack <- length(tracks)
      sizeblocs <- .05
  }

  ## HTClist object
  if (view==1){
      rx <- range(x)
      w <- as.numeric(width(rx))
      names(w) <- seqlevels(rx)
      w <- w[chrom]
    
    if(length(tracks) > 0){
      design <- matrix(NA, lc+1, lc+1)
      design[1,] <- c(1,seq(lc^2+2,(lc+1)^2-1,2))
      design[,1] <- c(1,seq(lc^2+3,(lc+1)^2,2))
      design[2:(lc+1), 2:(lc+1)] <- matrix(2:(lc^2+1), lc, lc, byrow=FALSE)

      heatspace <- 1-sizeblocs*ntrack
      layout(design, widths=c(sizeblocs*ntrack,round(w/sum(w)*heatspace,3)), heights=c(sizeblocs*ntrack,round(w/sum(w)*heatspace,3)))
      
      ##blank plot at position 1
      par(mar=c(0,0,0,0))
      plot(1, type="n", axes=FALSE, xlab="", ylab="")     
    }else{
      design <- matrix(1:lc^2, lc, lc, byrow=FALSE)
      layout(design, widths=round(w/sum(w),3), heights=round(w/sum(w),3))
    }
    ## HTCexp object
  }else if (view==2L){
    rx <- range(x)
    ## Annotation
    if(length(tracks) > 0L){
      if (!is.null(y)){
        ry <- range(y)
        if (width(rx)>=width(ry)){
          design <- rbind(rep(2L,3),rep(1L,3),c(4L,3L,5L),rep(6L,3))
          lmar <- (start(ry)-start(rx))/width(rx)
          rmar <- (end(rx)-end(ry))/width(rx)
          h <- (1-sizeblocs*ntrack)/2
          hy <- width(ry)/width(rx)
          layout(design, heights=c(h, sizeblocs*ntrack,h*hy,h*(1-hy)), widths=c(lmar,(1-rmar-lmar), rmar))
        }else{
          design <- rbind(rep(6L,3), c(4L,2L,5L), rep(1L,3),rep(3L,3)) 
          lmar <- (start(rx)-start(ry))/width(ry)
          rmar <- (end(ry)-end(rx))/width(ry)
          h <- (1-sizeblocs*ntrack)/2
          hx <- width(rx)/width(ry)
          layout(design, heights=c(h*(1-hx), h*hx, sizeblocs*ntrack,h), widths=c(lmar,(1-rmar-lmar), rmar))
        }
      }else{
        design <- matrix(1:2, 2, 1, byrow=TRUE)
        layout(design, heights=c(1-sizeblocs*ntrack, sizeblocs*ntrack))
      }
    }else{
      if (!is.null(y)){
        ry <- range(y)
        ## Adjust position of x and y if with are not the same
        if (width(rx)>=width(ry)){
          design <- matrix(c(1L,1L,1L,3L,2L,4L,5L,5L,5L), ncol=3, byrow=TRUE)
          lmar <- (start(ry)-start(rx))/width(rx)
          rmar <- (end(rx)-end(ry))/width(rx)
          hy <- width(ry)/width(rx)
          layout(design, widths=c(lmar,(1-rmar-lmar), rmar), heights=c(1L, hy, 1-hy))
          
        }else{
          design <- matrix(c(5L,5L,5L,3L,1L,4L,2L,2L,2L), ncol=3, byrow=TRUE)
          lmar <- (start(rx)-start(ry))/width(ry)
          rmar <- (end(ry)-end(rx))/width(ry)
          hx <- width(rx)/width(ry)
          layout(design, widths=c(lmar,(1-rmar-lmar), rmar), heights=c(1-hx, hx ,1L))
        }
      }else{
        design <- matrix(1:lc^2, lc, lc, byrow=TRUE)
        layout(design, widths=rep(1/lc, lc), heights=rep(1/lc, lc))
      }
    }
  }
}

###################################
## getMapData
## INTERNAL FUNCTION
## Operations to apply to the data before the visualization
##
## x = HTCexp object
## y = HTCexp object
## tracks = List of GRanges objects. Each object represent a genome track information
##
##################################

getData2Map <- function(x, minrange, maxrange, trim.range, log.data){
    stopifnot(inherits(x,"HTCexp"))
    xdata <- intdata(x)
    
    ## #####################
    ## Data Transformation
    ## #####################
    if (log.data){
        xdata <- log2(xdata)
    }

    ## #################
    ## Play with contrast
    ## #################
    if (trim.range <1 && is.na(maxrange) && is.na(minrange)){
        if (inherits(xdata, "sparseMatrix")){
            xmaxrange <- quantile(abs(xdata@x), probs=trim.range, na.rm=TRUE)
            xminrange <- quantile(abs(xdata@x), probs=1-trim.range, na.rm=TRUE)
        }else{
            xmaxrange <- quantile(abs(xdata@x[which(xdata@x!=0)]), probs=trim.range, na.rm=TRUE)
            xminrange <- quantile(abs(xdata@x[which(xdata@x!=0)]), probs=1-trim.range, na.rm=TRUE)
        }
    }
    else{
        if (is.na(maxrange))
            xmaxrange <- max(abs(xdata@x), na.rm=TRUE)
        else
            xmaxrange=maxrange
        if (is.na(minrange))
            xminrange <- min(abs(xdata@x), na.rm=TRUE)         
        else
            xminrange=minrange
    }

    ## The minrange/maxrange are simetrical around zero. The same threshold are used for positive and negative values
    xdata@x[which(xdata@x<=xminrange & xdata@x>0)] <- xminrange
    xdata@x[which(xdata@x>=-xminrange & xdata@x<0)] <- -xminrange
    xdata@x[which(xdata@x>=xmaxrange & xdata@x>0)] <- xmaxrange
    xdata@x[which(xdata@x<=-xmaxrange & xdata@x<0)] <- -xmaxrange

    print(paste("Minrange=",min(xdata, na.rm=TRUE)," - Maxrange=", max(xdata, na.rm=TRUE)))

    print(paste("minrange=",round(xminrange,6)," - maxrange=", round(xmaxrange,6)))

    xdata
}
       
###################################
## mapC methods
## 
## Visualization of HTCexp or HTClist objects
##
## x = HTCexp/HTClist object
## y = optional. HTCexp/HTClist object or matrix data
## tracks = List of GRanges objects. Each object represent a genome track information
## minrange = minimum value to draw
## maxrange = maximum value to draw
## trim.range = remove the outliers values by trimming the maxrange (quantile)
## names = logical, add axis with the bin/interval's names
## value = logical, add values on the matrix
## show.na = logical, na are shown in gray
## col.low = color used in mypalette for low count (low, mid, high)
## col.high = color used in mypalette for high count (low, mid, high)
## col.na = color used for NA values
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################

setMethod("mapC", signature="HTClist",
          function(x, tracks=NULL,
                   minrange=NA, maxrange=NA, trim.range=0.98, show.zero=FALSE, show.na=FALSE, log.data=FALSE, names=FALSE, value=FALSE,
                   col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", grid=FALSE){
              
              ## Set Graphical Environment
              setEnvDisplay(x, tracks=tracks, view=1)
              
              ## Get data to map and plots
              tmp <- sapply(names(pair.chrom(seqlevels(x))), function(i){
                  if (is.element(i,names(x))){
                      obj <- x[[i]]
                      xdata <- getData2Map(obj, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)
                      message("Plotting ",i,"...")

                      if (!names)
                          par(mar=c(0,0,0,0))
                      else
                          par(mar=c(mean(sapply(rownames(xdata),nchar))/2,mean(sapply(colnames(xdata),nchar))/2,0,0))

                      heatmapC(xdata, names=names, value=value, show.zero=show.zero,show.na=show.na, col.pos=col.pos,
                               col.neg=col.neg, col.na=col.na, grid=grid)
                  }else{
                      plot(1, type="n", axes=FALSE, xlab="", ylab="")
                  }
              })
              
              ## Add annotation based on intrachromosomal maps
              if (!is.null(tracks)){
                  tmp <- sapply(paste(seqlevels(x), seqlevels(x), sep=""), function(i){
                      if (is.element(i,names(x))){
                          obj <- x[[i]]
                          addImageTracks(obj, tracks, orientation="h", names=FALSE)
                          addImageTracks(obj, tracks, orientation="v", names=FALSE)
                }else{
                    warning("Intrachromosomal map for ",i," not found. Annotation skipped.")
                    plot(1, type="n", axes=FALSE, xlab="", ylab="")
                    plot(1, type="n", axes=FALSE, xlab="", ylab="")
                }
                  })
              }
          }
)


setMethod("mapC", signature="HTCexp",
          function(x, tracks=NULL,
                   minrange=NA, maxrange=NA, trim.range=0.98, value=FALSE, show.zero=FALSE, show.na=FALSE, log.data=FALSE,
                   col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", grid=FALSE, title=NULL){

              if (!isIntraChrom(x))
                  stop("The triangle view is available for intrachromosomal data only")

              ## Set Graphical Environment
              setEnvDisplay(x, tracks=tracks, view=2)
              
              ## Get data to map
              xdata <- getData2Map(x, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)

              ## Plots tracks and C map
              par(mar=c(0,0,0,0))
              triViewC(xdata, show.zero=show.zero, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, value=value, title=title)
              if (!is.null(tracks))
                addImageTracks(x, tracks, orientation="h")
          }
)


setMethod("mapC", signature=c("HTCexp","HTCexp"),
          function(x, y, tracks=NULL,
                   minrange=NA, maxrange=NA, trim.range=0.98, value=FALSE, show.zero=FALSE, show.na=FALSE, log.data=FALSE,
                   col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", grid=FALSE, title=NULL){

              if (!isBinned(x) || !isBinned(y))
                  stop("x and y have to be binned to plot them on the same scale")

              if (seqlevels(x) != seqlevels(y))
                  stop("x and y have to come from the same chromosome")
              
              ## Set Graphical Environment
              setEnvDisplay(x, y, tracks=tracks, view=2)
           
              ## Get data to map and plots
              xdata <- getData2Map(x, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)
              ydata <- getData2Map(y, minrange=minrange, maxrange=maxrange, trim.range=trim.range, log.data=log.data)

              ## Plots tracks and C map
              if (!is.null(tracks)){
                if (width(range(x))>=width(range(y))){
                  addImageTracks(x, tracks, orientation="h")
                }else{
                  addImageTracks(y, tracks, orientation="h")
                }
              }
              par(mar=c(.5,0,0,0))
              triViewC(xdata, value=value, show.zero=show.zero, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, title=title[1])
              par(mar=c(0,0,.5,0))
              triViewC(ydata, flip=TRUE, value=value, show.zero=show.zero, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, title=title[2])
            }
)

###################################
## plot Alias
##################################

setMethod("plot", signature="HTClist",
          function(x, ...){
              mapC(x, ...)
          }
)

setMethod("plot", signature="HTCexp",
          function(x, ...){
              mapC(x, ...)
          }
)

setMethod("plot", signature=c("HTCexp","HTCexp"),
          function(x, y, ...){
              mapC(x, y, ...)
          }
)
