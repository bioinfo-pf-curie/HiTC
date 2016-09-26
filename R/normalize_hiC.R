## Nicolas Servant
## HiTC BioConductor package
##**********************************************************************************************************##
##
## HIC Normalization procedure from Lieberman-Aiden et al. 2009
##
##**********************************************************************************************************##

## Normalized per expected number of count
setMethod("normPerExpected", signature=c("HTCexp"), definition=function(x, ...){
    
    expCounts <- getExpectedCounts(forceSymmetric(x), asList=TRUE, ...)
    if (! is.null(expCounts$stdev.estimate)){
        x@intdata <- (x@intdata-expCounts$exp.interaction)/expCounts$stdev.estimate
    }else{
        x@intdata <- x@intdata/(expCounts$exp.interaction) 
    }
    ## Remove NaN or Inf values for further analyses
    #x@intdata[which(is.na(x@intdata) | is.infinite(x@intdata))]<-NA
    x@intdata[Matrix::which(is.infinite(x@intdata))]<-0
    x
})

## Normalized per expected number of counts across all cis maps
setMethod("normPerExpected", signature=c("HTClist"), definition=function(x, ...){

  xintra <- x[isIntraChrom(x)]

  ## estimated expected counts for all cis maps
  exp <- lapply(xintra, function(xx){
      r <- getExpectedCounts(forceSymmetric(xx), method="mean", asList=TRUE, ...)
      r$exp.interaction
  })

  ## combined all cis expected counts
  N <- max(sapply(exp, dim))
  counts <- matrix(0, ncol=N, nrow=N)
  ss <- matrix(0, ncol=N, nrow=N)
  for (i in 1:length(exp)){
    n <- dim(exp[[i]])[1]
    counts[1:n, 1:n] <- counts[1:n, 1:n]+1
    ss[1:n, 1:n] <- ss[1:n, 1:n] + as.matrix(exp[[i]])
  }
 
  ## Mean over all expected matrices
  ss <- ss / counts
  
  xintranorm <- lapply(xintra, function(xx){
    n <- dim(xx@intdata)[1]
    xx@intdata <- xx@intdata/ss[1:n,1:n]
    xx@intdata[which(is.na(xx@intdata) | is.infinite(xx@intdata))]<-0
    xx
  })
  x[isIntraChrom(x)] <- xintranorm
  x
})

###################################
## getExpectedCountsMean
##
## This way of calculate expected counts was used in Naumova et al.
## The idea is just to look at all diagonals and to calculate their mean
##
## x = a HTCexp object
##
##
## NOTES
## Migth be interesting to add an isotonic regression on the mean to force the expected value to decrease with the distance
###################################


getExpectedCounts <- function(x, method=c("mean","loess"), asList=FALSE, ...){
  met <- match.arg(method)

  if (dim(intdata(x))[1]>500 & met=="loess"){
    warning("Contact map looks big. Use mean method instead or be sure that the loess fit gives good results.")
  }
  
  if (met=="mean"){
    ret <- getExpectedCountsMean(x, ...)
  }else if (met=="loess"){
    ret <- getExpectedCountsLoess(x, ...)
  }else{
    stop("Unknown method")
  }

  if (asList){
    return(ret)
  }else{
    intdata(x) <- ret$exp.interaction
    return(x)
  }
}


logbins<- function(from, to, step=1.05, N=NULL) {
  if (is.null(N)){
    unique(round(c(from, exp(seq(log(from), log(to), by=log(step))), to)))
  }else{
    unique(round(c(from, exp(seq(log(from), log(to), length.out=N)), to)))
  }
}

getExpectedCountsMean <- function(x, logbin=TRUE, step=1.05, filter.low=0.05){

  xdata <- intdata(x)
  N <- dim(xdata)[1]
  if (logbin){
    bins <- logbins(from=1,to=N, step=step)
    bins <- as.vector(Rle(values=bins, lengths=c(diff(bins),1)))
    stopifnot(length(bins)==N)
  }else{
    bins <- 1:N
  }

  message("Estimate expected using mean contact frequency per genomic distance ...")
  
  xdata <- as.matrix(xdata)
  rc <- colSums(xdata, na.rm=TRUE)
  ##rc <- which(rc==0)
  rc <- which(rc < ceiling(quantile(rc[which(rc>0)], probs=filter.low)))
  rr <- rowSums(xdata, na.rm=TRUE)
  ##rr <- which(rr==0)
  rr <- which(rr <  ceiling(quantile(rr[which(rr>0)], probs=filter.low)))

  ## rm line with only zeros
  xdata[rr,] <- NA
  xdata[,rc] <- NA
           
  ## create an indicator for all diagonals in the matrix
  rows <- matrix(rep.int(bins, N), nrow=N)
  ##d <- rows - t(rows)

  d <- matrix(bins[1+abs(col(rows) - row(rows))],nrow=N) - 1
  d[lower.tri(d)] <- -d[upper.tri(d)]
 
  if (isSymmetric(xdata)){
    ## remove half of the matrix
    d[lower.tri(d)] <- NA
  }
  ## use split to group on these values
  mi <- split(xdata, d)
  milen <- lapply(mi, length)
  mimean <- lapply(mi, mean, na.rm=TRUE)
  miexp <- lapply(1:length(milen), function(i){rep(mimean[[i]], milen[[i]])})
  names(miexp) <- names(mi)
  expmat <- as(matrix(unsplit(miexp, d), nrow=nrow(xdata), ncol=ncol(xdata)), "Matrix")
  if (isSymmetric(xdata)){
    expmat <- forceSymmetric(expmat, uplo="U")
  }
  
  colnames(expmat) <- colnames(xdata)
  rownames(expmat) <- rownames(xdata)

  ## Put NA at rc and cc
  expmat[rr,] <- NA
  expmat[,rc] <- NA
  return(list(exp.interaction=expmat, stdev.estimate=NULL))
}


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

getExpectedCountsLoess<- function(x, span=0.01, bin=0.005, stdev=FALSE, plot=FALSE){
    stopifnot(inherits(x,"HTCexp"))

    xdata <- as.matrix(intdata(x))
    rc <- which(colSums(xdata, na.rm=TRUE)==0)
    rr <- which(rowSums(xdata, na.rm=TRUE)==0)

    ## rm line with only zeros
    xdata[rr,] <- NA
    xdata[,rc] <- NA
      
    ydata <- as.vector(xdata)
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

        ##plotIntraDist(ydata, xdata.dist, xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="Counts", main="", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, col="gray", frame=FALSE)
        
        plot(x=xdata.dist, y=ydata,  xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="Counts", main="", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, col="gray", frame=FALSE)
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
    
    ## Put NA at rc and cc
    lowess.mat[rr,] <- NA
    lowess.mat[,rc] <- NA
 
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
##
##################################

normLGF <- function(x,  family=c("poisson", "nb")){
    family <- match.arg(family)
    message("Starting LGF normalization on ", seqlevels(x), " ...")
    counts <- intdata(x)
    
    ## Remove rowCounts=0 & colCounts=0
    rc <- which(rowSums(counts)>0)

    ## Intrachromosomal maps
    if (isIntraChrom(x)){
        cc <- rc
        stopifnot(length(rc)>0)
        counts.rc <- counts[rc,rc]

        elt <- elementMetadata(y_intervals(x)[rc])
        len <- elt$len
        gcc <- elt$GC
        map <- elt$map
        
        if(all(is.na(len)) || all(is.na(gcc)) || all(is.na(map)))
            stop("Genomic features are missing. Effective fragments length, GC content and mappability are required.")
        
        ##get cov matrix
        len_m<-as.matrix(log(1+len%o%len))
        gcc_m<-as.matrix(log(1+gcc%o%gcc))
        
        ##error for regions with 0 mappability
        map[which(map==0)] <- 10e-4
        map_m<-as.matrix(log(map%o%map))      

    }else{

    ## Interchromosomal maps
        cc <- which(colSums(counts)>0)
        stopifnot(length(rc)>0 & length(cc)>0)
        
        counts.rc <- counts[rc,cc]
        yelt <- elementMetadata(y_intervals(x)[rc])
        xelt <- elementMetadata(x_intervals(x)[cc])
        
        ylen <- yelt$len
        xlen <- xelt$len
        ygcc <- yelt$GC
        xgcc <- xelt$GC
        ymap <- yelt$map
        xmap <- xelt$map

        if(all(is.na(ylen)) || all(is.na(ygcc)) || all(is.na(ymap)) || all(is.na(xlen)) || all(is.na(xgcc)) || all(is.na(xmap)))
            stop("Genomic features are missing. Effective fragments length, GC content and mappability are required.")
      
        ##get cov matrix
        len_m<-as.matrix(log(1+ylen%o%xlen))
        gcc_m<-as.matrix(log(1+ygcc%o%xgcc))
        
        ##error for regions with 0 mappability
        ymap[which(ymap==0)] <- 10e-4
        xmap[which(xmap==0)] <- 10e-4
        map_m<-as.matrix(log(ymap%o%xmap))      
    }
    
    ##centralize cov matrix of enz, gcc
    len_m<-(len_m-mean(len_m, na.rm=TRUE))/apply(len_m, 2, sd, na.rm=TRUE)
    gcc_m<-(gcc_m-mean(gcc_m, na.rm=TRUE))/apply(gcc_m, 2, sd, na.rm=TRUE)
    
    ##change matrix into vector
    if (isIntraChrom(x)){
        counts_vec<-counts.rc[which(upper.tri(counts.rc,diag=FALSE))]
        len_vec<-len_m[upper.tri(len_m,diag=FALSE)]
        gcc_vec<-gcc_m[upper.tri(gcc_m,diag=FALSE)]
        map_vec<-map_m[upper.tri(map_m,diag=FALSE)]
    }else{
        counts_vec<-as.vector(counts.rc)
        len_vec<-as.vector(len_m)
        gcc_vec<-as.vector(gcc_m)
        map_vec<-as.vector(map_m)
    }
    
    print("fit ...")
    if (family=="poisson"){
        ##fit Poisson regression: u~len+gcc+offset(map)
        fit<-glm(counts_vec~ len_vec+gcc_vec+offset(map_vec),family="poisson")
        ##fit<-bigglm(counts_vec~len_vec+gcc_vec+offset(map_vec),family="poisson", data=cbind(counts_vec, len_vec, gcc_vec, map_vec))
    }else{
        fit<-glm.nb(counts_vec~len_vec+gcc_vec+offset(map_vec))
    }

    coeff<-fit$coeff

    ## The corrected values (residuals) can be seen as a observed/expected correction.
    ## So I will compare the normalized counts with one: the observed count is higher or lower than the expected count. We may not want to compare the range of the normalized count with the range of the raw count. They have different interpretations.
    counts.cor<-round(counts.rc/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
    counts[rownames(counts.rc), colnames(counts.rc)]<-counts.cor
  
    intdata(x) <- counts
    return(x)
}##normLGF

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
    stopifnot(seqlevels(x) %in% seqlevels(cutSites))
    obj <- x
    xgi <- x_intervals(x)
    message("Annotation of ", seqlevels(x), " ...")
    xgi <- annotateIntervals(xgi, cutSites[[seqlevels(xgi)]], minfragmap=minFragMap, efffraglen=effFragLen)
    x_intervals(obj) <- xgi
    if (isIntraChrom(x) & isBinned(x)){
        y_intervals(obj) <- xgi
    }else{
        ygi <- y_intervals(x)
        ygi <- annotateIntervals(ygi, cutSites[[seqlevels(ygi)]], minfragmap=minFragMap, efffraglen=effFragLen)
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
    if (!is.na(minfragmap) & !all(is.na(annot$map_U)) & !all(is.na(annot$map_D))){
        idxmap <- which(annot$map_U<minfragmap | is.na(annot$map_U))
        elementMetadata(annot)[idxmap,c("len_U", "GC_U", "map_U")]<-NA_real_
        idxmap <- which(annot$map_D<minfragmap | is.na(annot$map_D))
        elementMetadata(annot)[idxmap,c("len_D", "GC_D", "map_D")]<-NA_real_
    }

    ## Get all restriction sites which overlap with the bins
    ## Split upstream and downstream bins to deal with restriction sites which overlap the start or end of a fragment
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
  
    annotscores <- lapply(1:length(outl_up), function(i){
        id_up <- outl_up[[i]]
        id_dw <- outl_dw[[i]]
        ##temp <-  c(annot_up[id_up], annot_down[id_dw])
        temp_up <- annot_up[id_up]
        temp_dw <- annot_down[id_dw]

        ## len - effective length" is the fragment length truncated by 1000 bp, which is the number of bases with specific ligation.
        ## In Yaffe & Tanay's paper Figure 1b, they define specific ligation as sum of distance to cutter sites (d1+d2) <= 500 bp. Such criterion implies that d1<=500 bp and d2 <= 500 bp. So for each fragment end, only reads mapped within 500 bp to cutter sites are used for downstream analysis. 
        lenv <- unique(c(temp_up$len, temp_dw$len))
        if (!is.na(efffraglen))
            lenscore <- sum(lenv>efffraglen, na.rm=TRUE)*efffraglen + sum(lenv[lenv<efffraglen], na.rm=TRUE)
        else
            lenscore <- sum(lenv, na.rm=TRUE)
        
        ##GC
        gcscore <- mean(c(temp_up$GC, temp_dw$GC), na.rm=TRUE)

        ##map
        mapscore <- mean(c(temp_up$map, temp_dw$map), na.rm=TRUE)

        c(lenscore, gcscore, mapscore)
    })
    annotscores <- matrix(unlist(annotscores), ncol=3, byrow=TRUE)
    colnames(annotscores) <- c("len", "GC", "map")
    elementMetadata(gi)$len <- round(annotscores[,"len"],3)
    elementMetadata(gi)$GC <- round(annotscores[,"GC"],3)
    elementMetadata(gi)$map <- round(annotscores[,"map"],3)

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

getAnnotatedRestrictionSites <- function(resSite="AAGCTT", overhangs5=1, chromosomes=NULL, genomePack="BSgenome.Mmusculus.UCSC.mm9",  mappability=NULL, wingc=200, winmap=500){

    if(genomePack %in% loadedNamespaces()==FALSE){
        stopifnot(require(genomePack, character.only=TRUE))
    }
    genome <- eval(as.name(genomePack))

    if (is.null(chromosomes)){
        chromosomes <- seqlevels(genome)
    }
    
    genomeCutSites <- mclapply(chromosomes, function(chr){
        message("Get restriction sites for ", chr, " ...")
        cutSites <- getRestrictionSitesPerChromosome(resSite, overhangs5, genome, chr)
        message(length(cutSites), " sites")
        
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
        message("done ...")
        cutSites
    })
    grl <- GRangesList(genomeCutSites)
    names(grl) <- chromosomes
    grl
}


###################################
## getRestrictionSitesPerChromosome
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
    stopifnot(length(chromosome)==1)
    
    restrictionSites<-Biostrings::matchPattern(resSite, genome[[chromosome]])
 
    ## Deal with restriction enzyme 5' overhangs
    s <- start(restrictionSites) + overhangs5
    e <- end(restrictionSites) - overhangs5

    ir <- IRanges(start=s, end=e)    
    restrictionSites <- GRanges(seqnames = chromosome, ranges = ir, strand = "*")
    return(restrictionSites)
}

###################################
## getRestrictionFragmentsPerChromosome
## 
##
## resSite = Cutting site of the restriction enzyme used
## overhangs5 =  Cleavage 5 overhang
## genome = BSgenome object of the reference genome
## chromosome = chromosome to focus on
##
##################################

getRestrictionFragmentsPerChromosome <- function(resSite="AAGCTT", overhangs5=1,
chromosomes=NULL, genomePack="BSgenome.Mmusculus.UCSC.mm9"){
    
    if(genomePack %in% loadedNamespaces()==FALSE){
        stopifnot(require(genomePack, character.only=TRUE))
    }
    genome <- eval(as.name(genomePack))
    stopifnot(inherits(genome,"BSgenome"))

    if (is.null(chromosomes)){
        chromosomes <- seqlevels(genome)
    }
  
    genomeResFrag <- mclapply(chromosomes, function(chromosome){
        message("Get restriction fragments for ", chromosome, " ...")
        restrictionSites<-getRestrictionSitesPerChromosome(resSite, overhangs5, genome, chromosome)        
        restrictionFrag <- GRanges(seqnames=chromosome, ranges=IRanges(
                                                         start=c(1,start(restrictionSites)),
                                                         end=c(start(restrictionSites)-1, seqlengths(genome)[chromosome])), strand="+")
    })
    return(genomeResFrag)
}


##**********************************************************************************************************##
##
## ICE Normalization procedure from Imakaev et al .2012
##
##**********************************************************************************************************##

###################################
## balancingSK
## INTERNAL FUNCTION
##
## Matrix balancing used in ICE normalization
## Based on the Sinkhorn-Knopp algorithm
##
## x = HTCexp object
## max_iter = maximum number of iteration to converge
## eps = threshold to converge
##
##################################

balancingSK<- function(x, max_iter=50, eps=1e-4){
    m <- dim(x)[1]

    ## Initialization    
    sum_ss <- matrix(rep(0, m), ncol=1)
    bias <- matrix(rep(1, m), ncol=1)
    old_dbias <- NULL
    ## Remove Diagonal ?
    
    for (it in 1:max_iter){
        message("it=",it," ", Sys.time())

        ## 1- calculate sum of W over all rows ++
        sum_ds <- rowSums(x, na.rm=TRUE)
        ##sum_ds <- sqrt(rowSums(x^2))

        ## 2- Calculate a vector of corrected ss reads
        ## NOT DONE
        
        ## 3- Calculate vector of bias
        dbias <- as.matrix(sum_ds, ncol=1) + sum_ss

        ## 4 - Renormalize bias by its mean valude over non-zero bins to avoid numerical instabilities
        dbias <- dbias/mean(dbias[dbias!=0])

        ## 5- Set zero values of bias to 1 to avoir 0/0 error
        dbias[dbias==0] <- 1
               
        ## 6- Divide W by bias BiBj for all (i,j) ++++
        x <- x/(dbias %*% t(dbias))

        ## 7- Multiple total vector of bias by additional biases
        ##bias <- bias * dbias

        if (!is.null(old_dbias) && sum(abs(old_dbias - dbias))<eps){
            message("Break at iteration ", it)
            break
        }
        old_dbias <- dbias 
    }
    if (it == max_iter){
        message("Did not converged. Stop at iteration ",max_iter)
    }else{
        message("Converge in ",it," iteration")
    }
    return(x)
}


###################################
## IterativeCorNormalization
## ICE normlization
## 
##
## x = HTCexp object or HTClist object
## max_iter = maximum number of iteration to converge
## eps = threshold to converge
## spars.filter = Percentage of row and column to discard based on sparsity (default=0.02)
##
##################################

normICE <- function(x, max_iter=50, eps=1e-4, sparse.filter=0.02){

    if (inherits(x, "HTCexp")){
        stopifnot(isSymmetric(x))
        idata <- intdata(x)
        gr <- y_intervals(x)
    }else if (inherits(x, "HTClist")){
        idata <- getCombinedContacts(x)
        gr <- getCombinedIntervals(x)
    }

    if (!is.na(sparse.filter)){
        message("Start filtering  ...", Sys.time())
        spars <- apply(idata, 1, function(x){length(which(x==0))}) 
        spars.t <- quantile(spars[spars!=dim(idata)[1]], probs=(1-sparse.filter))
        idx <- which(spars>as.numeric(spars.t))
        idata[idx,] <- 0
        idata[,idx] <- 0
        message("Filter out ",length(idx)," rows and columns ...")
    }
    
    message("Start Iterative Correction ...")
    xmat <- balancingSK(idata, max_iter=max_iter, eps=eps)
    
    if (inherits(x, "HTCexp")){
        intdata(x) <- xmat
    }else if (inherits(x, "HTClist")){
        ##     gr <- dimnames2gr(xmat, pattern="\\||\\:|\\-", feat.names=c("name","chr","start", "end"))
        ##     xgi <- gr[[1]]
        ##     ygi <- gr[[2]]
        ##     rownames(xmat) <- id(ygi)
        ##     colnames(xmat) <- id(xgi)
        if (is.null(gr$xgi))
            x <- splitCombinedContacts(xmat, xgi=gr$ygi, ygi=gr$ygi)
        else
            x <- splitCombinedContacts(xmat, xgi=gr$xgi, ygi=gr$ygi)
    }
    x
}
