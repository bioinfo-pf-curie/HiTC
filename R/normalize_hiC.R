## Nicolas Servant
## HiTC BioConductor package
##**********************************************************************************************************##
##
## HIC Normalization procedure from Lieberman-Aiden et al. 2009
##
##**********************************************************************************************************##

## Normalized per expected number of count
setMethod("normPerExpected", signature=c("HTCexp"), definition=function(x, stdev=FALSE, ...){
    expCounts <- getExpectedCounts(x, stdev=stdev, ...)
    if (stdev){
        x@intdata <- (x@intdata-expCounts$exp.interaction)/expCounts$stdev.estimate
    }else{
        x@intdata <- x@intdata/expCounts$exp.interaction
    }
    ## Remove NaN or Inf values for further analyses
    intdata(x)[which(is.na(x@intdata) | is.infinite(x@intdata))]<-NA
    x
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
    lowess.mat <- Matrix(lowess.fit[order(o)], nrow=length(y_intervals(x)), byrow=FALSE)
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
        stdev.delta <- unlist(mclapply(1:length(ind), function(k){
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
        message("Calculating stdev ... ")
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



##**********************************************************************************************************##
##
## HIC Normalization procedure from HiCNorm package, Hu et al. 2012
##
##**********************************************************************************************************##

###################################
## normLGF
## Local Genomic Features normalization
## 
##
## x = HTCexp/HTClist object
## family = regression model Poisson or Neg Binon
##
## TODO
## - trans data
##
##################################

normLGF <- function(x,  family=c("poisson", "nb")){
    family <- match.arg(family)
    message("start ", seqlevels(x))
    counts <- intdata(x)
    
    ##remove rowCounts=0
    rc <- which(rowSums(counts)==0)
    if (length(rc)>0){
        counts.rc <- counts[-rc,-rc]
        elt <- elementMetadata(x_intervals(x)[-rc])
    }
    else{
        counts.rc <- counts
        elt <- elementMetadata(x_intervals(x))
    }

    if(dim(counts.rc)[1]==0){
        warning("Unable to normalize ",seqlevels(x_intervals(x))," interaction map")
        invisible(NULL)
    }else{
        len <- elt$len
        gcc <- elt$GC
        map <- elt$map

        if(is.null(len) || is.null(gcc) || is.null(map))
            stop("Genomic features are missing. Effective fragments length, GC content and mappability are required.")
        
        ##get cov matrix
        len_m<-as.matrix(log(1+len%o%len))
        gcc_m<-as.matrix(log(1+gcc%o%gcc))
        
        ##error for regions with 0 mappability
        map[which(map==0)] <- 10e-4
        map_m<-as.matrix(log(map%o%map))      

        ##centralize cov matrix of enz, gcc
        len_m<-(len_m-mean(len_m, na.rm=TRUE))/apply(len_m, 2, sd, na.rm=TRUE)
        gcc_m<-(gcc_m-mean(gcc_m, na.rm=TRUE))/apply(gcc_m, 2, sd, na.rm=TRUE)
        
        ##change matrix into vector
        counts_vec<-counts.rc[which(upper.tri(counts.rc,diag=FALSE))]
        len_vec<-len_m[upper.tri(len_m,diag=FALSE)]
        gcc_vec<-gcc_m[upper.tri(gcc_m,diag=FALSE)]
        map_vec<-map_m[upper.tri(map_m,diag=FALSE)]
        
        print("fit ...")
        if (family=="poisson"){
            ##fit Poisson regression: u~len+gcc+offset(map)
            fit<-glm(counts_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")
            ##fit<-bigglm(counts_vec~len_vec+gcc_vec+offset(map_vec),family="poisson", data=cbind(counts_vec, len_vec, gcc_vec, map_vec))
        }else{
            fit<-glm.nb(counts_vec~len_vec+gcc_vec+offset(map_vec))
        }
        
        coeff<-round(fit$coeff,4)
        counts.cor<-round(counts.rc/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
        counts[rownames(counts.rc), colnames(counts.rc)]<-counts.cor
        
        intdata(x) <- counts
        return(x)
    }
}

###################################
## setGenomicFeatures
## Annotate a HTCexp or HTClist object with the GC content and the mappability features
## 
##
## x = HTCexp/HTClist object
## cutSites = GRanges object ir GRangesList from getAnnotatedRestrictionSites function
## minFragMap = Discard restriction with mappability lower the this threshold (and NA)
## effFragLen = Effective fragment length
##################################

setGenomicFeatures <- function(x, cutSites, minFragMap=.5, effFragLen=1000){
    stopifnot(inherits(x,"HTCexp"))
    if (inherits(cutSites, "GRangesList")){
        cutSites <- cutSites[[seqlevels(x)]]
        seqlevels(cutSites) <- seqlevels(x)
    }
    stopifnot(seqlevels(x)==seqlevels(cutSites))
    obj <- x
    xgi <- x_intervals(x)
    xgi <- annotateIntervals(xgi, cutSites, minfragmap=minFragMap, efffraglen=effFragLen)
    x_intervals(obj) <- xgi
    if (isIntraChrom(x) && isBinned(x)){
        y_intervals(obj) <- xgi
    }else{
        ygi <- y_intervals(x)
        ygi <- annotateIntervals(ygi, cutSites, minfragmap=minFragMap, efffraglen=effFragLen)
        y_intervals(obj) <- ygi
    }
   obj
}

###################################
## annotateIntervals
## INTERNAL FUNCTION
## 
##
## gi = GRanges object from x_intervals or y_intervals methods
## annot = GRanges object from getAnnotatedRestrictionSites function
##
##################################

annotateIntervals <- function(gi, annot, minfragmap=.5, efffraglen=1000){

    ## Preprocess, keep fragments ends with mappability score larger than .5.
    ## Depends on the data processing (see Yaffe and Tanay, 2011). These fragments will be exclude from the analysis.
    if (!is.na(minfragmap)){
        idxmap <- which(annot$map_U<minfragmap | is.na(annot$map_U))
        elementMetadata(annot)[idxmap,c("len_U", "GC_U", "map_U")]<-NA_real_
        idxmap <- which(annot$map_D<minfragmap | is.na(annot$map_D))
        elementMetadata(annot)[idxmap,c("len_D", "GC_D", "map_D")]<-NA_real_
    }

    ## Get all restriction sites which overlap with the bins
    ## Slip upstream and downstream bins to deal with restriction sites which overlap the start or end of a fragment
    annot_up <- annot
    end(annot_up)<-start(annot)
    elementMetadata(annot_up) <- NULL
    annot_up$len=as.numeric(annot$len_U)
    annot_up$GC=as.numeric(annot$GC_U)
    annot_up$map=as.numeric(annot$map_U)
    annot_down <- annot
    start(annot_down) <- end(annot)
    elementMetadata(annot_down) <- NULL
    annot_down$len=as.numeric(annot$len_D)
    annot_down$GC=as.numeric(annot$GC_D)
    annot_down$map=as.numeric(annot$map_D)

    outl_up<- as.list(findOverlaps(gi, annot_up))
    outl_dw<- as.list(findOverlaps(gi, annot_down))
    
    annotscores <- t(sapply(1:length(outl_up), function(i){
        id_up <- outl_up[[i]]
        id_dw <- outl_dw[[i]]
        temp <-  c(annot_up[id_up], annot_down[id_dw])

        ## len - effective length" is the fragment length truncated by 1000 bp, which is the number of bases with specific ligation.
        ## In Yaffe & Tanay's paper Figure 1b, they define specific ligation as sum of distance to cutter sites (d1+d2) <= 500 bp. Such criterion implies that d1<=500 bp and d2 <= 500 bp. So for each fragment end, only reads mapped within 500 bp to cutter sites are used for downstream analysis. 
        lenv <- unique(temp$len)
        if (!is.na(efffraglen))
            lenscore <- sum(lenv>efffraglen, na.rm=TRUE)*efffraglen + sum(lenv[lenv<efffraglen], na.rm=TRUE)
        else
            lenscore <- sum(lenv, na.rm=TRUE)
        
        ##GC
        gcscore <- mean(temp$GC, na.rm=TRUE)

        ##map
        mapscore <- mean(temp$map, na.rm=TRUE)

        c(lenscore, gcscore, mapscore)
    }))
    
    colnames(annotscores) <- c("len", "GC", "map")
    elementMetadata(gi)<-round(annotscores,3)
    gi
}


###################################
## getAnnotatedRestrictionFragments
## Return the restriction fragments for a given enzyme, annotated with the GC content and the mappability
## 
##
## resSite = Cutting site of the restriction enzyme used (default HindIII)
## overhangs5 =  Cleavage 5 overhang
## chromosome = chromosomes list to focus on. If NULL, all genome chromosome are investigated
## genomePack = name of the BSgenome package to load
## w = size of the downstream/upstream window to use around the restriction site to calculate the GC content. Default is 200. See Yaffe and Tanay for more details
## mappability = GRanges object of the mappability (see the ENCODE mappability tracks)
##
## D = downstream / U = upstream the restriction site
##################################

getAnnotatedRestrictionSites <- function(resSite="AAGCTT", overhangs5=1, chromosomes=NA, genomePack="BSgenome.Mmusculus.UCSC.mm9",  mappability=NULL, wingc=200, winmap=500){

    if(genomePack %in% loadedNamespaces()==FALSE){
        stopifnot(require(genomePack, character.only=TRUE))
    }
    genome <- eval(as.name(genomePack))

    if (length(chromosomes)==0){
        chromosomes <- seqlevels(genome)
    }
    
    genomeCutSites <- mclapply(chromosomes, function(chr){
        message("Get restriction sites for ", chr, " ...")
        cutSites <- getRestrictionSitesPerChromosome(resSite, overhangs5, genome, chr)
        
        message("Calculate fragment length ...")
        ## Add chromosome start/end
        len_D <- c(end(cutSites)[-1], length(genome[[chr]])) - start(cutSites)
        len_U <- end(cutSites) - c(0, start(cutSites)[-length(cutSites)])
        cutSites$len_U <- len_U
        cutSites$len_D <- len_D

        message("Calculate GC content ...")
        ## Upstream GC content
        win <- start(cutSites)-wingc
        win[win<0] <- 1
        seq <- Biostrings::getSeq(genome, chr, start=win, end=start(cutSites)-1)
        ##cutSites$seq_U <- seq
        cutSites$GC_U<- round(Biostrings::letterFrequency(seq, as.prob=FALSE, letters="CG")/Biostrings::letterFrequency(seq, as.prob=FALSE, letters="ACGT"),3)
        
        ## Downstream GC content
        win <- start(cutSites)+wingc-1
        win[win>length(genome[[chr]])] <- length(genome[[chr]])
        seq <- Biostrings::getSeq(genome, chr, start(cutSites), win)
        cutSites$GC_D<- round(Biostrings::letterFrequency(seq, as.prob=FALSE, letters="CG")/Biostrings::letterFrequency(seq, as.prob=FALSE, letters="ACGT"),3)
        ##cutSites$seq_D <- seq

        if (!is.null(mappability)){
            message("Calculate mappability ...")
            stopifnot(inherits(mappability,"GRanges"))
            
            mappability <- mappability[seqnames(mappability)==chr]
            win <- start(cutSites)-winmap+1
            win[win<0] <- 1
            gr <- GRanges(seqnames = chr, ranges = IRanges(start=win, end=start(cutSites)))
            overl <- as.list(findOverlaps(gr, mappability))
            mscore <- mappability$score
            cutSites$map_U<- unlist(lapply(overl, function(idx){
                round(mean(mscore[idx], na.rm=TRUE),3)
            }))
        
            win <- start(cutSites)+winmap
            win[win>length(genome[[chr]])] <- length(genome[[chr]])
            gr <- GRanges(seqnames = chr, ranges = IRanges(start=start(cutSites)+1, end=win))
            overl <- as.list(findOverlaps(gr, mappability))
            mscore <- mappability$score
            cutSites$map_D<- unlist(lapply(overl, function(idx){
                round(mean(mscore[idx], na.rm=TRUE),3)
            }))
        }else{
            cutSites$map_U<-NA_real_
            cutSites$map_D<-NA_real_
        }
        cutSites
    })
    grl <- GRangesList(genomeCutSites)
    names(grl) <- chromosomes
    grl
}


###################################
## getRestrictionFragmentsPerChromosome
## INTERNAL FUNCTION
## 
##
## resSite = Cutting site of the restriction enzyme used
## overhangs5 =  Cleavage 5 overhang
## genome = BSgenome object of the reference genome
## chromosome = chromosome to focus on
##
##################################

getRestrictionSitesPerChromosome <- function(resSite, overhangs5, genome, chromosome){

    stopifnot(inherits(genome,"BSgenome"))

    restrictionSites<-Biostrings::matchPattern(resSite, genome[[chromosome]], fixed=FALSE)
 
    ## Deal with restriction enzyme 5' overhangs
    s <- start(restrictionSites) + overhangs5
    e <- end(restrictionSites) - overhangs5

    ir <- IRanges(start=s, end=e)    
    restrictionFrag <- GRanges(seqnames = chromosome, ranges = ir, strand = "*")
    return(restrictionFrag)
}

##**********************************************************************************************************##
##
## ICE Normalization procedure from Imakaev et al .2012
##
##**********************************************************************************************************##

###################################
## IterativeCorNormalization
## INTERNAL FUNCTION
## 
##
## x = HTCexp object
## max_iter = maximum number of iteration to converge
## eps = threshold to converge
##
##################################

IterativeCorNormalization <- function(x, max_iter=200, eps=1e-4){
    m <- dim(x)[1]
    ## for single end reads
    sum_ss <- matrix(rep(0, m), ncol=1)

    #bias <- matrix(rep(1, m), ncol=1)
    old_dbias <- NULL
    
    for (it in 1:max_iter){
        sum_ds <- rowSums(x, na.rm=TRUE)
        ##sum_ds <- sqrt(rowSums(x^2))
        dbias <- as.matrix(sum_ds, ncol=1) + sum_ss
        dbias <- dbias/mean(dbias[dbias!=0])
        dbias[dbias==0] <- 1
        ##bias <- bias * dbias
        
        ## normalize by the dbias matrix product
        x <- x/dbias %*% t(dbias)
        
        if (!is.null(old_dbias) && sum(abs(old_dbias - dbias))<eps){
            message("Break at iteration ", it)
            break
        }
        old_dbias <- dbias 
    }
    if (it == max_iter){
        warning("Did not converged. Stop at iteration ",max_iter)
    }
    return(x)
}

###################################
## IterativeCorNormalization
## ICE normlization
## 
##
## x = HTCexp object
## max_iter = maximum number of iteration to converge
## eps = threshold to converge
## spars.filter = filter out the (default 2%) bins the more sparse
##
##################################

normICE <- function(x, max_iter=200, eps=1e-4, sparse.filter=0.02){

    stopifnot(isSymmetric(x))
    idata <- intdata(x)
    srow <- rowSums(idata, na.rm=TRUE)

    spars <- apply(idata, 1, function(x){length(which(x==0))})
    f <- quantile(spars[spars!=dim(idata)[1]], probs=sparse.filter)
    idx <- spars[which(spars<as.numeric(f))]
    idata[idx,] <- 0
    
    message("start ", seqlevels(x))
    xmat <- IterativeCorNormalization(idata, max_iter=max_iter, eps=eps)
    intdata(x) <- xmat
    message("end ", seqlevels(x))
    x
}
