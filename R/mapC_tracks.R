###################################
## getBlocsIndex
## INTERNAL FUNCTION
## Overlap between a GRanges objects (HTCexp vs tracks)
## Return the index (start/end) of overlap region based on the HTCexp ranges
##
## gr = GRanges object, i.e. 'C' xgi for vertical view or ygi for horizontal view. 
## track = list of GRanges objects, i.e. tracks information 
##
##################################

getBlocsIndex <- function(gr, track){
    suppressWarnings(r <- countOverlaps(gr, track, ignore.strand=TRUE))
    ov <- cumsum(width(Rle(ifelse(r==0,0,1))))

    if (!is.null(ov)){
        if (r[length(r)]==0)
            ov <- ov[-length(ov)]
        if (r[1]!=0)
            ov <- c(0,ov)
        blocs <- matrix(ov, ncol=2, byrow=TRUE)
        blocs[,1] <- blocs[,1]+1
    }
    blocs
}

###################################
## addImageTracks
## INTERNAL FUNCTION
## Add a genome track information on interaction map
##
## x = 'HTCexp' object
## tracks = list of GRanges objects, i.e. tracks information 
## orientation = plot vertical or horizontal tracks
##
##################################

addImageTracks <- function(x, tracks, orientation=c("h","v"), names=TRUE){

  stopifnot(inherits(x,"HTCexp"))
  if (inherits(x, "GRanges")){
    x <- list(x)
  }
  stopifnot(unlist(lapply(tracks, inherits,"GRanges")))

  ntrack <- length(tracks)
  suppressWarnings(colblocs <- brewer.pal(ntrack*2,"Paired"))
  colblocs.minus <- colblocs[seq(1, ntrack*2, by=2)]
  colblocs.plus <- colblocs[seq(2, ntrack*2, by=2)]
  blocnames <- names(tracks)
  
  if (orientation=="h"){
    ypos <- 1
    if (isBinned(x)){
      par(mar=c(0,0,0,0))
      plot(c(start(range(x_intervals(x))),end(range(x_intervals(x)))), c(-1, ntrack*3+1),type="n",axes=FALSE, frame=FALSE, xlab="", ylab="", xlim=c(start(range(x_intervals(x))), end(range(x_intervals(x)))), xaxs="i", yaxs="i")
      for (t in 1:ntrack){
        blocs <- tracks[[t]]
        ## keep only overlaping features between
        blocs <- subsetByOverlaps(blocs,range(x_intervals(x)), ignore.strand=TRUE)
        ## get strand information
        blocs.plus <- blocs[which(strand(blocs)=="+" | strand(blocs)=="*")]
        blocs.minus <- blocs[which(strand(blocs)=="-")]
        ## draw features
        if (length(blocs.plus)>0)
          rect(start(blocs.plus), ypos+.1, end(blocs.plus), ypos+.6, col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
          rect(start(blocs.minus), ypos-.1, end(blocs.minus), ypos-0.6, col=colblocs.minus[t], border=colblocs.minus[t])
        if (names)
          text(x=start(range(x_intervals(x)))+width(range(x_intervals(x)))/2, y=ypos+1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t])
        ypos <- ypos+3
      }
    }else{
      warning("The data are not binned, thus the scale is not linear.\nOnly the HTCexp intervals overlapping with the track's feature are displayed.")
      cset <- x_intervals(x)
      par(mar=c(0,0,0,0))
      plot(c(0,length(cset)), c(-1, ntrack*3+1),type="n", axes=FALSE, xlab="", ylab="", frame=FALSE, xlim=c(0,length(cset)), xaxs="i", yaxs="i")
      for (t in 1:ntrack){
        blocs <- tracks[[t]]
        ## keep only overlaping features between
        suppressWarnings(blocs <- subsetByOverlaps(blocs,range(x_intervals(x)), ignore.strand=TRUE))
        ## get strand information and index of cset overlapping with the annotation features
        blocs.plus <- getBlocsIndex(cset, blocs[which(strand(blocs)=="+" | strand(blocs)=="*")])
        blocs.minus <- getBlocsIndex(cset, blocs[which(strand(blocs)=="-")])
        ## draw features
        if (length(blocs.plus)>0)
          rect(blocs.plus[,1]-1, ypos+.1, blocs.plus[,2], ypos+.6, col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
                    rect(blocs.minus[,1]-1, ypos-.1, blocs.minus[,2], ypos-0.6, col=colblocs.minus[t], border=colblocs.minus[t])
        if (names)
          text(x=length(cset)/2, y=ypos+1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t])
        ypos <- ypos+3
      }
        }
  }else{
    ypos <- -1
    if (isBinned(x)){
      par(mar=c(0,0,0,0))
      plot(c(1, -ntrack*3-1),c(start(range(y_intervals(x))),end(range(y_intervals(x)))), type="n",axes=FALSE, xlab="", ylab="", frame=FALSE, ylim=c(start(range(y_intervals(x))),end(range(y_intervals(x)))), xaxs="i", yaxs="i")
      
      for (t in 1:ntrack){
        blocs <- tracks[[t]]
        ## keep only overlaping features between
        blocs <- subsetByOverlaps(blocs,range(y_intervals(x)), ignore.strand=TRUE)
        ## get strand information
        blocs.plus <- blocs[which(strand(blocs)=="+" | strand(blocs)=="*")]
        blocs.minus <- blocs[which(strand(blocs)=="-")]
        ## draw features    
        if (length(blocs.plus)>0)
          rect(ypos-.1, (start(range(y_intervals(x)))+end(range(y_intervals(x))))-start(blocs.plus),
               ypos-.6, (start(range(y_intervals(x)))+end(range(y_intervals(x))))-end(blocs.plus),
               col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
          rect(ypos+.1, (start(range(y_intervals(x)))+end(range(y_intervals(x))))-start(blocs.minus),
               ypos+.6, (start(range(y_intervals(x)))+end(range(y_intervals(x))))-end(blocs.minus),
               col=colblocs.minus[t], border=colblocs.minus[t])
        if (names)
          text(y=start(range(y_intervals(x)))+width(range(y_intervals(x)))/2, x=ypos-1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t], srt=90)
        ypos <- ypos-3
      }
    }else{
      warning("The data are not binned, and the scale is not linear.\nOnly the HTCexp intervals overlapping with the track's feature are displayed.")
      cset <- y_intervals(x)
      par(mar=c(0,0,0,0))
            plot(c(0, -ntrack*3-1),c(0,length(cset)),type="n", axes=FALSE, xlab="", ylab="", frame=FALSE, ylim=c(0,length(cset)), xaxs="i", yaxs="i")
      
      for (t in 1:ntrack){
        blocs <- tracks[[t]]
        ## keep only overlaping features between
        suppressWarnings(blocs <- subsetByOverlaps(blocs,range(y_intervals(x)), ignore.strand=TRUE))
        ## get strand information and index of cset overlapping with the annotation features
        blocs.plus <- getBlocsIndex(cset, blocs[which(strand(blocs)=="+" | strand(blocs)=="*")])
        blocs.minus <- getBlocsIndex(cset, blocs[which(strand(blocs)=="-")])
        ##draw features
        if (length(blocs.plus)>0)
          rect(ypos-.1, length(cset)-blocs.plus[,1]+1, ypos-.6, length(cset)-blocs.plus[,2], col=colblocs.plus[t], border=colblocs.plus[t])
        if (length(blocs.minus)>0)
          rect( ypos+.1, length(cset)-blocs.minus[,1]+1,  ypos+.6,length(cset)- blocs.minus[,2], col=colblocs.minus[t], border=colblocs.minus[t])
        if (names)
          text(y=length(cset)/2, x=ypos-1, labels=blocnames[t], cex=.7, font=2, col=colblocs.plus[t], srt=90)
                ypos <- ypos-3
      }
    }
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

