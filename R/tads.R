

###########################################
##
## Directionality Index
##
############################################

directionalityIndex <- function(x, winup=2e6, windown=2e6){
    stopifnot(isBinned(x))
    xdata <- intdata(x)
    rx <- ranges(x_intervals(x))
    upranges <- IRanges(start(rx)-winup+1, start(rx))
    A <- rep(NA, length(upranges))
    
    dwranges <- IRanges(end(rx), end(rx)+windown-1)
    B <- rep(NA, length(upranges))

    overup <- findOverlaps(rx, upranges,  type="within")
    overdown <- findOverlaps(rx, dwranges,  type="within")

    ## For all bins
    for (i in 1:length(rx)){
        ## What are the interactions in the up/downstream windows
        idx <- which(subjectHits(overup)==i)
        if (length(idx)>0){
            idx.rx <- queryHits(overup)[idx]
            A[i] <- sum(xdata[idx.rx, i], na.rm=TRUE)
        }
        idx <- which(subjectHits(overdown)==i)
        if (length(idx)>0){
            idx.rx <- queryHits(overdown)[idx]
            B[i] <- sum(xdata[idx.rx, i], na.rm=TRUE)
        }
        #print(paste0(A[i], "-", B[i], "-", A[i]+B[i]/2))
        #if(i>=100){stop()}
    }

    E <- (A+B)/2
    DI <- ((B-A)/abs(B-A))*(((A-E)^2)/E + ((B-E)^2)/E)
    DI[which(is.na(DI))] <- 0
    DI
}

