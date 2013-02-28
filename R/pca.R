###################################
## pca.hic
##
## PCA analysis on interactio map as proposed in Lieberman et al.
##
## x = HTCexp object
## norm = if true, apply the normPerExpected normalization
## npc = number of component to report
## asGRangesList = if TRUE a GRangesList object is returned where the score represents the eigenvector
##################################

pca.hic <- function(x, normPerExpected=TRUE, npc=2, asGRangesList=TRUE){
    stopifnot(inherits(x,"HTCexp"))
    stopifnot(isBinned(x))
    
    ## observed/expected map using the loess calculation of distance dependency
    if (normPerExpected){
        x <- normPerExpected(x, stdev=FALSE)
    }
    
    ## Remove bins with sd=0 to avoid correlation failure
    idx <- union(names(which(apply(intdata(x), c(1), sd, na.rm=TRUE)==0)),
                 names(which(apply(intdata(x), c(2), sd, na.rm=TRUE)==0)))
    x <- removeIntervals(x,idx)
    
    ## Remove bins with NAs
    idx <- union(names(which(apply(intdata(x), 1, function(x) all(is.na(x))))),
                 names(which(apply(intdata(x), 2, function(x) all(is.na(x))))))
    x <- removeIntervals(x,idx)


    ## Calculate correlation
    ## use="pairwise.complete.obs" means that the correlation between two vectors is calculated using only the paired non-NAs values.
    ## It means that two correlations can be calculated using a different set of points ...
    ## In theory this never appends because bins with NAs values should be removed earlier
    xdata.cor <- cor(intdata(x), use="pairwise.complete.obs")

    ## Perform PCA
    pca <- prcomp(xdata.cor, scale=TRUE)
    
    ## Get results
    if (!asGRangesList){
        pca.res <- list()
        pca.res$varianceProp <- summary(pca)$importance[2,1:npc]
        for (i in 1:npc){
            pca.res[[eval(paste("PC",i, sep=""))]] <- pca$rotation[,i]
        }
    }else{
        pca.res <- GRangesList()
        xgi <- x_intervals(x)
        for (i in 1:npc){
            pca.res[[eval(paste("PC",i, sep=""))]] <- GRanges(seqnames(xgi), ranges=ranges(xgi), strand=strand(xgi), score=round(pca$rotation[,i],3))
        }
    }
    
    return(pca.res)
} 
