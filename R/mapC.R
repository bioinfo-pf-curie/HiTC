###################################
## myPalette
## INTERNAL FUNCTION
## Define a palette for the heatmap
##
## low = color for minimal value
## high = color for extreme values (up and down)
## mid = color for medium values
## k = number of color generated
##
##################################

myPalette <- function (low="white", high=c("green", "red"), mid=NA, k = 50){
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
}

###################################
## gdiag
## INTERNAL FUNCTION
## Extract the diagonal of a matrix - Generatlisation of diag function
##
## x = a matrix
## w = step from the real diagonal. w=0 return the diagonal of the matrix.
##
##################################

gdiag <- function(x, w=0){
    if (is.matrix(x)) {
        if ((m <- min(dim(x))) == 0L) 
            return(vector(typeof(x), 0L))
         if (w<m){
            if (w>0)
                y <- c(x)[1L + (0L+w):(m - 1L)* (dim(x)[1L] + 1L) - w]
            else
                y <- c(x)[1L + 0L:(m - 1L - abs(w))* (dim(x)[1L] + 1L) + abs(w)]
            return(y)
        }else{
            return(vector(typeof(x), 0L)) 
        }
    }
}
###################################
## addImageTracks
## INTERNAL FUNCTION
## Add a genome track information on interaction map
##
## cset = 'C' xgi for vertical view or ygi for horizontal view. Used for to define the tracks positions
## giblocs = List of Genome_intervals objects, i.e. tracks information 
## orientation = plot vertical or horizontal tracks
##
##################################
     
addImageTracks <- function(x, giblocs, orientation=c("h","v")){
  stopifnot(inherits(x,"HTCexp"))
  
  ntrack <- length(giblocs)
  stopifnot(unlist(lapply(giblocs, inherits,"Genome_intervals")))
  
  colblocs <- brewer.pal(ntrack*2,"Paired")
  colblocs.minus <- colblocs[seq(1, ntrack*2, by=2)]
  colblocs.plus <- colblocs[seq(2, ntrack*2, by=2)]
  blocnames <- names(giblocs)
  
  if (orientation=="h"){     
    ypos <- 1
    if (isBinned(x)){
      par(mar=c(0,0,0,0))
      plot(c(range(x)[1],range(x)[2]), c(-1, ntrack*3+1),type="n",axes=FALSE, frame=FALSE, xlab="", ylab="", xlim=range(x), xaxs="i", yaxs="i")
      
      for (t in 1:ntrack){
        blocs <- giblocs[[t]]
        blocs.plus <- blocs.minus <- NULL
        if (is.element("strand", names(annotation(blocs)))){
          blocs.plus <- blocs[which(annotation(blocs)$strand=="+")]
          blocs.minus <- blocs[which(annotation(blocs)$strand=="-")]
        }else{
          blocs.plus <- blocs
        }
        
        if (length(blocs.plus)>0)
          rect(blocs.plus[,1], ypos+.1, blocs.plus[,2], ypos+.6, col=colblocs.plus[t], border=colblocs.plus[t])
        
        if (length(blocs.minus)>0)
          rect(blocs.minus[,1], ypos-.1, blocs.minus[,2], ypos-0.6, col=colblocs.minus[t], border=colblocs.minus[t])
        
        text(x=range(x)[1]+diff(range(x))/2, y=ypos+1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t])
        ypos <- ypos+3
        
      }
    }else{
      cset <- x_intervals(x)
      par(mar=c(0,0,0,0))
      plot(c(0,dim(cset)[1]), c(-1, ntrack*3+1),type="n", axes=FALSE, xlab="", ylab="", frame=FALSE, xlim=c(0,dim(cset)[1]), xaxs="i", yaxs="i")
      for (t in 1:ntrack){
        blocs <- giblocs[[t]]
        blocs.plus <- blocs.minus <- NULL
        ##Which intervals in 'to' overlap with 'from'
        if (is.element("strand", names(annotation(blocs)))){
          ovplus <- unlist(lapply(interval_overlap(from=blocs[which(annotation(blocs)$strand=="+")], to=cset), function(x){if(length(x)>0){return(range(x))}}))
          if (!is.null(ovplus))
            blocs.plus <- matrix(ovplus, ncol=2, byrow=TRUE)
          ovminus <- unlist(lapply(interval_overlap(from=blocs[which(annotation(blocs)$strand=="-")], to=cset), function(x){if(length(x)>0){return(range(x))}}))
          if (!is.null(ovminus))
            blocs.minus <- matrix(ovminus, ncol=2, byrow=TRUE)
        }else{
          ov <- unlist(lapply(interval_overlap(from=blocs, to=cset), function(x){if(length(x)>0){return(range(x))}}))
          blocs.plus <- matrix(ov, ncol=2, byrow=TRUE)
        }
        
        if (length(blocs.plus)>0)
          rect(blocs.plus[,1]-1, ypos+.1, blocs.plus[,2], ypos+.6, col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
          rect(blocs.minus[,1]-1, ypos-.1, blocs.minus[,2], ypos-0.6, col=colblocs.minus[t], border=colblocs.minus[t])
        text(x=dim(cset)[1]/2, y=ypos+1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t])
        ypos <- ypos+3
      }
    }
  }else{
    ypos <- -1
    if (isBinned(x)){
      par(mar=c(0,0,0,0))
      plot(c(1, -ntrack*3-1),c(range(x)[1],range(x)[2]), type="n",axes=FALSE, xlab="", ylab="", frame=FALSE, ylim=range(x), xaxs="i", yaxs="i")
      
      for (t in 1:ntrack){
        blocs <- giblocs[[t]]
        blocs.plus <- blocs.minus <- NULL
        if (is.element("strand", names(annotation(blocs)))){
          blocs.plus <- blocs[which(annotation(blocs)$strand=="+")]
          blocs.minus <- blocs[which(annotation(blocs)$strand=="-")]
        }else{
          blocs.plus <- blocs
        }
           
        if (length(blocs.plus)>0)
          rect(ypos-.1, sum(range(x))-blocs.plus[,1], ypos-.6, sum(range(x))-blocs.plus[,2], col=colblocs.plus[t], border=colblocs.plus[t])
        
        if (length(blocs.minus)>0)
          rect(ypos+.1, sum(range(x))-blocs.minus[,1], ypos+.6, sum(range(x))-blocs.minus[,2], col=colblocs.minus[t], border=colblocs.minus[t])
        
        text(y=range(x)[1]+diff(range(x))/2, x=ypos-1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t], srt=90)
        ypos <- ypos-3
      }
    }else{
      cset <- y_intervals(x)
      par(mar=c(0,0,0,0))
      plot(c(0, -ntrack*3-1),c(0,dim(cset)[1]),type="n", axes=FALSE, xlab="", ylab="", frame=FALSE, ylim=c(0,dim(cset)[1]), xaxs="i", yaxs="i")
      
      for (t in 1:ntrack){
        blocs <- giblocs[[t]]
        blocs.plus <- blocs.minus <- NULL

        ##Which intervals in 'to' overlap with 'from'
        if (is.element("strand", names(annotation(blocs)))){
          ovplus <- unlist(lapply(interval_overlap(from=blocs[which(annotation(blocs)$strand=="+")], to=cset), function(x){if(length(x)>0){return(range(x))}}))
          if (!is.null(ovplus))
            blocs.plus <- matrix(ovplus, ncol=2, byrow=TRUE)
          ovminus <- unlist(lapply(interval_overlap(from=blocs[which(annotation(blocs)$strand=="-")], to=cset), function(x){if(length(x)>0){return(range(x))}}))
          if (!is.null(ovminus))
            blocs.minus <- matrix(ovminus, ncol=2, byrow=TRUE)
        }else{
          ov <- unlist(lapply(interval_overlap(from=blocs, to=cset), function(x){if(length(x)>0){return(range(x))}}))
          blocs.plus <- matrix(ov, ncol=2, byrow=TRUE)
        }
        
        if (length(blocs.plus)>0)
          rect(ypos-.1, dim(cset)[1]-blocs.plus[,1]+1, ypos-.6, dim(cset)[1]-blocs.plus[,2], col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
          rect( ypos+.1, dim(cset)[1]-blocs.minus[,1]+1,  ypos+.6,dim(cset)[1]- blocs.minus[,2], col=colblocs.minus[t], border=colblocs.minus[t])
        text(y=dim(cset)[1]/2, x=ypos-1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t], srt=90)
        ypos <- ypos-3
      }
    }
  }
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
## col.low = color used in mypalette for low count (low, mid, high)
## col.high = color used in mypalette for high count (low, mid, high)
## col.na = color used for NA values
## mask.data = add a mask on the matrix. Has to have the same dimension as the input matrix
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################

heatmapC <- function(xdata,  names=FALSE,  value=FALSE,  show.na=TRUE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", mask.data=NULL, grid=FALSE, title=NULL){
    xdata <- t(xdata)
    xdata <- xdata[,ncol(xdata):1]

    xdata.pos <- xdata.neg <- xdata
    xdata.pos[xdata.neg<0] <- NA
    xdata.neg[xdata.pos>=0] <- NA

    k <- length(unique(as.vector(abs(xdata))))
    if (length(unique(xdata.neg[xdata.neg<0 & !is.na(xdata.neg)]))>0){
        col.neg <- myPalette(col.neg[3],col.neg[1],mid=col.neg[2], k=k)
        image(x=1:nrow(xdata),y=1:ncol(xdata),z=xdata.neg,axes=FALSE,ylab="",xlab="",col=col.neg)
        par(new = TRUE)
    }
 
    if (length(unique(xdata.pos[xdata.pos>0 & !is.na(xdata.pos)]))>0){
        col.pos <- myPalette(col.pos[1],col.pos[3],mid=col.pos[2],k=k)
        image(x=1:nrow(xdata),y=1:ncol(xdata),z=xdata.pos,axes=FALSE,ylab="",xlab="",col=col.pos)
    }
 
    if (names){
        axis(side=2, at=1:ncol(xdata), labels=colnames(xdata), lwd=0.5, las=1, cex.axis=0.7)  
        axis(side=1, at=1:nrow(xdata), labels=rownames(xdata), lwd=0.5, las=2, cex.axis=0.7)
    }
    if (value){
        if (is.null(mask.data)){
            for (i in nrow(xdata):1){
                text(x=i,y=1:nrow(xdata),labels=round(xdata[i,],2), cex=0.7)
            }
        }
        else{
            if (ncol(xdata) != ncol(mask.data) || nrow(xdata) != nrow(mask.data)){
                stop("Mask data have wrong dimension")
            }
            mask.data <- t(mask.data)
            mask.data <- mask.data[,ncol(mask.data):1]
            for (i in ncol(mask.data):1){
                text(x=i,y=1:nrow(mask.data),labels=round(mask.data[i,],2), cex=0.7)
            }
        }
    }
 
    ###################
    ## Heatmap options
    ###################
    if (show.na){
        if (length(which(is.na(xdata))>0)){
            par(new = TRUE)
            na.xdata <- matrix(NA, ncol=ncol(xdata), nrow=nrow(xdata))
            na.xdata[which(is.na(xdata))] <- 1
            image(x=1:nrow(na.xdata),y=1:ncol(na.xdata),z=na.xdata,axes=FALSE,ylab="",xlab="",col=col.na,add=TRUE)
        }
    }
    
    if (grid){
        abline(h=seq(0.5,nrow(xdata)+0.5, by=1), mar=c(0,0))
        abline(v=seq(0.5,ncol(xdata)+0.5, by=1), mar=c(0,0))
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
## names = logical, add axis with the bin/interval's names
## value = logical, add values on the matrix
## show.na = logical, na are shown in gray
## col.low = color used in mypalette for low count (low, mid, high)
## col.high = color used in mypalette for high count (low, mid, high)
## col.na = color used for NA values
## mask.data = add a mask on the matrix. Has to have the same dimension as the input matrix
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################

triViewC <- function(xdata, flip=FALSE, show.na=TRUE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80"){
    
    d <- min(dim(xdata))
    trimat <- matrix(NA, ncol=d*2, nrow=d)
    for (w in 0:(d-1)){
        s <- gdiag(xdata, w=w)
        ss <- c(s,s)[as.vector(sapply(1:length(s),function(x){return(c(x,x+length(s)))}))]
        trimat[w+1, (w+1):(d*2-w)] <- ss
    }

    if (flip){
        trimat <- trimat[nrow(trimat):1,]
    }
    
    trimat.pos <- trimat.neg <- trimat
    trimat.pos[trimat.neg<0] <- NA
    trimat.neg[trimat.pos>=0] <- NA
    
    k <- length(unique(as.vector(abs(trimat))))
    
    if (length(unique(trimat.neg[trimat.neg<0 & !is.na(trimat.neg)]))>0){
        col.neg <- myPalette(col.neg[3],col.neg[1],mid=col.neg[2],k=k)
        image(y=1:nrow(trimat),x=1:ncol(trimat),z=t(trimat.neg),axes=FALSE,ylab="",xlab="",col=col.neg)
        par(new = TRUE)
    }

    if (length(unique(trimat.pos[trimat.pos>0 & !is.na(trimat.pos)]))>0){
        col.pos <- myPalette(col.pos[1],col.pos[3],mid=col.pos[2],k=k)
        image(y=1:nrow(trimat),x=1:ncol(trimat),z=t(trimat.pos),axes=FALSE,ylab="",xlab="",col=col.pos)
    }
 
    if (show.na){
        if (length(which(is.na(trimat))>0)){
            par(new = TRUE)
            na.trimat <- matrix(NA, ncol=ncol(trimat), nrow=nrow(trimat))
            na.trimat[which(is.na(trimat))] <- 1
            image(y=1:nrow(na.trimat),x=1:ncol(na.trimat),z=t(na.trimat),axes=FALSE,ylab="",xlab="",col=col.na)
        }
    }
    
}


###################################
## mapC
##
## Draw heatmap of the C data
##
## x = HTCexp object or matrix data
## y = optional. HTCexp object or matrix data
## view = 1 - heatmapC, 2 - triViewC
## giblocs = list of Genome_intervals objects. Each object represent a genome track information
## minrange = minimum value to draw
## maxrange = maximum value to draw
## trim.range = remove the outliers values by trimming the maxrange (quantile)
## names = logical, add axis with the bin/interval's names
## value = logical, add values on the matrix
## show.na = logical, na are shown in gray
## col.low = color used in mypalette for low count (low, mid, high)
## col.high = color used in mypalette for high count (low, mid, high)
## col.na = color used for NA values
## mask.data = add a mask on the matrix. Has to have the same dimension as the input matrix
## grid = logical, the grid is shown
## title = add a title to the plot
##
##################################


mapC <- function(x, y=NULL, view=1, giblocs=NULL, minrange=NA, maxrange=NA, trim.range=0.98, names=FALSE, value=FALSE, show.na=FALSE, log.data=FALSE, col.pos=c("white",NA,"red"), col.neg=c("white",NA,"blue"), col.na="gray80", mask.data=NULL, grid=FALSE, title=NULL){

    if(missing(x)){
        stop("Missing x input")
    }

    if (inherits(x,"HTCexp"))
        xdata <- intdata(x)
    else if(is.matrix(x))
        xdata <- x
    else
        stop("Input has to belong to 'HTCexp' or 'matrix' classes")
    
    if (!is.null(y)){
        view=2
        if (inherits(y,"HTCexp")){
            ydata <- intdata(y)
        }
        else if(is.matrix(y))
            ydata <- y
        else
            stop("Input2 has to belong to 'HTCexp' or 'matrix' classes")
    }
    
    ##Logged data
    if (log.data){
        xdata[which(xdata<0)] <- NA
        xdata[which(xdata>0)] <- log2(xdata[which(xdata>0)])
        if (!is.null(y)){
            ydata[which(ydata<0)] <- NA
            ydata[which(ydata>0)] <- log2(ydata[which(ydata>0)])
        }
    }
    
    ## Check binned data for sample comparison
    if (!is.null(y)){
      if (!isBinned(x) || !isBinned(y))
        stop("x and y have to be binned to plot them on the same scale")
    }
    
    ###################
    ## Play with contrast
    ###################
    if (trim.range <1 && is.na(maxrange) && is.na(minrange)){
        xmaxrange <- quantile(abs(xdata[xdata!=0]), probs=trim.range, na.rm=TRUE)
        xminrange <- quantile(abs(xdata[xdata!=0]), probs=1-trim.range, na.rm=TRUE)
        if (!is.null(y)){
            ymaxrange <- quantile(abs(ydata[ydata!=0]), probs=trim.range, na.rm=TRUE)
            yminrange <- quantile(abs(ydata[ydata!=0]), probs=1-trim.range, na.rm=TRUE)
        }
        
    }
    else{
        if (is.na(maxrange)){
            xmaxrange <- max(abs(xdata[xdata!=0]), na.rm=TRUE)
            if (!is.null(y))
                ymaxrange <- max(abs(ydata[ydata!=0]), na.rm=TRUE)
        }
        else{
            xmaxrange=maxrange
            if (!is.null(y))
                ymaxrange=maxrange
        }
        if (is.na(minrange)){
            xminrange <- min(abs(xdata[xdata!=0]), na.rm=TRUE)
            if (!is.null(y))
                yminrange <- min(abs(ydata[ydata!=0]), na.rm=TRUE)
        }
        else{
            xminrange=minrange
            if (!is.null(y))
                yminrange=minrange
        }
    }
    print(paste("minrange=",round(xminrange,6)," - maxrange=", round(xmaxrange,6)))
    xdata[which(xdata<=xminrange & xdata>0)] <- xminrange
    xdata[which(xdata>=-xminrange & xdata<0)] <- -xminrange
    xdata[which(xdata>=xmaxrange & xdata>0)] <- xmaxrange
    xdata[which(xdata<=-xmaxrange & xdata<0)] <- -xmaxrange

   if (!is.null(y)){
       print(paste("minrange=",round(yminrange,6)," - maxrange=", round(ymaxrange,6)))
       ydata[which(ydata<=yminrange & ydata>0)] <- yminrange
       ydata[which(ydata>=-yminrange & ydata<0)] <- -yminrange
       ydata[which(ydata>=ymaxrange & ydata>0)] <- ymaxrange
       ydata[which(ydata<=-ymaxrange & ydata<0)] <- -ymaxrange
   }

    ###################
    ##Graphical design
    ###################
    if(!is.null(giblocs)){
        if (!inherits(x,"HTCexp")){
            warning("Cannot diplay genomeIntervals blocs. 'x' has to be a HTCexp object.")
        }else{
            if (names){
                names <- FALSE
                warning("Cannot diplay names and genomeIntervals blocs.")
            }
            if (!is.list(giblocs))
                giblocs <- list(giblocs)
            
            stopifnot(unlist(lapply(giblocs, inherits,"Genome_intervals")))
            ntrack <- length(giblocs)
            sizeblocs <- .1#.05
            ygi <- y_intervals(x)
            xgi <- x_intervals(x)
        }
    }
    
    if (view==1){
        if(!is.null(giblocs)){
            design <- matrix(1:4, 2, 2, byrow=TRUE)
                        layout(design, widths=c(sizeblocs*ntrack,1-sizeblocs*ntrack), heights=c(sizeblocs*ntrack,1-sizeblocs*ntrack))
            
            ##blank plot at position 1
            par(mar=c(0,0,0,0))
            plot(1, type="n", axes=FALSE, xlab="", ylab="")
            
            addImageTracks(x, giblocs, orientation="h")
            addImageTracks(x, giblocs, orientation="v")
          }else{
            layout(matrix(1, 1, 1, byrow=TRUE), heights=c(1))
          }
      }else if (view == 2){
        rx <- range(x)
        
        if (!is.null(y)){
          ##Check if non overlap between x and y
          if (range(x)[1]<range(y)[1] && range(x)[2]<range(y)[2])
            y <- extractRegion(y, from=range(y)[1], to=range(x)[2], exact=TRUE)
          else if (range(x)[1]>range(y)[1] && range(x)[2]>range(y)[2])
            y <- extractRegion(y, from=range(x)[1], to=range(y)[2], exact=TRUE)
          else if (range(y)[1]<range(x)[1] && range(y)[2]<range(x)[2])
            y <- extractRegion(y, from=range(x)[1], to=range(y)[2], exact=TRUE)
          else if(range(y)[1]>range(x)[1] && range(y)[2]>range(x)[2])
            y <- extractRegion(y, from=range(y)[1], to=range(x)[2], exact=TRUE)
          
          ry <- range(y)
        }
        
        if(!is.null(giblocs)){
            if (!is.null(y)){
                if (diff(rx)>=diff(ry)){
                    design <- rbind(rep(2,3),rep(1,3), c(4,3,5)) 
                    lmar <- (ry[1]-rx[1])/diff(rx)
                    rmar <- (rx[2]-ry[2])/diff(rx)
                    layout(design, heights=c((1-sizeblocs*ntrack)/2, sizeblocs*ntrack,(1-sizeblocs*ntrack)/2), widths=c(lmar,(1-rmar-lmar), rmar))
                    
                    addImageTracks(x, giblocs, orientation="h")
                }else{
                    design <- rbind( c(4,2,5), rep(1,3),rep(3,3)) 
                    
                    lmar <- (rx[1]-ry[1])/diff(ry)
                    rmar <- (ry[2]-rx[2])/diff(ry)
                    layout(design, heights=c((1-sizeblocs*ntrack)/2, sizeblocs*ntrack,(1-sizeblocs*ntrack)/2), widths=c(lmar,(1-rmar-lmar), rmar))
                    addImageTracks(x, giblocs, orientation="h")
                }
            }else{
                design <- matrix(2:1, 2, 1, byrow=TRUE)
                layout(design, heights=c(1-sizeblocs*ntrack, sizeblocs*ntrack))
                addImageTracks(x, giblocs, orientation="h")
            }
        }else{
            if (!is.null(y)){
                if (diff(rx)>=diff(ry)){
                    design <- matrix(c(1,1,1,3,2,4), ncol=3, byrow=TRUE)
                    lmar <- (ry[1]-rx[1])/diff(rx)
                    rmar <- (rx[2]-ry[2])/diff(rx)
                    layout(design, widths=c(lmar,(1-rmar-lmar), rmar), heights=c(1,1))
                }else{
                    design <- matrix(c(3,1,4,2,2,2), ncol=3, byrow=TRUE)
                    lmar <- (rx[1]-ry[1])/diff(ry)
                    rmar <- (ry[2]-rx[2])/diff(ry)
                    layout(design, widths=c(lmar,(1-rmar-lmar), rmar), heights=c(1,1))
                }
            }else{
                layout(matrix(1, 1, 1, byrow=TRUE), heights=c(1))
            }
        }
    }
    
    ###################
    ##Graphical view
    ###################
    
    if (!names)
        par(mar=c(0,0,0,0))
    else
        par(mar=c(mean(sapply(rownames(xdata),nchar))/2,mean(sapply(colnames(xdata),nchar))/2,0,0))
    
    if (view==1)
        heatmapC(xdata, names=names, value=value, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na, mask.data=mask.data, grid=grid, title=title)
    else{
        triViewC(xdata, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na)
        if (!is.null(y))
            triViewC(ydata, flip=TRUE, show.na=show.na, col.pos=col.pos, col.neg=col.neg, col.na=col.na) 
    }
}

###################################
## discretize
##
## Transform matrix of counts data into discrete matrix
##
## x = data matrix with 5C interactions (counts)
## quant = if true, use quantile, else just split into equals 'nb.lev' levels
## nb.lev = number of level
##
###################################

discretize <- function (x, nb.lev = 4, quant=TRUE) {

    stopifnot(is.matrix(x))
    out <- x
    if (quant) {
        lev <- quantile(x, seq(0, 1, by = 1/nb.lev), na.rm=TRUE)
    }else{
        mind <- min(x, na.rm=TRUE)
        maxd <- max(x, na.rm=TRUE)
        diff <- (maxd - mind)/nb.lev
        lev=seq(0, maxd, by=diff)
    }
    
    if (length(which(duplicated(lev)))>0){
        lev <- unique(lev)
        warning(paste("Find only ",length(lev)-1," levels (nb.lev=",nb.lev,")", sep=""))
    }

    for (k in 1:length(lev)){
        out[which(x>=lev[k] & x<=lev[k+1])] <- k
    }
    out
}

