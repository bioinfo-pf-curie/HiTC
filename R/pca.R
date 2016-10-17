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

pca.hic <- function(x, normPerExpected=TRUE, npc=2, asGRangesList=TRUE, gene.gr=NULL, ...){
    stopifnot(inherits(x,"HTCexp"))
    stopifnot(isBinned(x))
    
    xcor <- getPearsonMap(x, normPerExpected, ...)
    xdata.cor <- intdata(xcor)
    xgi <- x_intervals(x)
    
    ## remove NA if still there
    idx <- which(apply(xdata.cor, 1, function(x){length(which(is.na(x)))}) != dim(xdata.cor)[1])
    xdata.cor <- xdata.cor[idx,idx]

    ## Perform PCA
    if (length(xdata.cor) == 0){
      warning("Empty correlation matrix for ", seqlevels(x))
      pca.res <- NULL
    }
    else{
      pca <- prcomp(xdata.cor, scale=TRUE)

      ## Gene density
      if (!is.null(gene.gr)){
          stopifnot(inherits(gene.gr,"GenomicRanges"))
          gene.density <- countOverlaps(xgi, gene.gr)
      }
     
      ## Get results
      if (!asGRangesList){
        pca.res <- list()
        pscore <- rep(NA, length(xgi))
        pca.res$varianceProp <- summary(pca)$importance[2,1:npc]
        for (i in 1:npc){
            pscore[idx] <- round(pca$rotation[,i],3)
                       
            if (!is.null(gene.gr)){
              gd.pos <- sum(gene.density[which(pscore>=0)])
              gd.neg <- sum(gene.density[which(pscore<0)])
              message("Gene density per bin - >0 = ", gd.pos)
              message("Gene density per bin - <0 = ", gd.neg)
              
              ## Add Chromosome compartment information
              if (gd.pos > gd.neg)
                cc <- ifelse(pscore>0, "A", "B")
              else{
                cc <- ifelse(pscore<0, "A", "B")
                pscore <- -pscore
              }
              pca.res[[eval(paste("PC",i, sep=""))]] <- rbind(pscore, cc, gene.density)
            }else{
              pca.res[[eval(paste("PC",i, sep=""))]] <- pscore
            }
          }
      }else{
        pca.res <- GRangesList()      
        for (i in 1:npc){
          pscore <- rep(NA, length(xgi))
          pscore[idx] <- round(pca$rotation[,i],3)
            
          if (!is.null(gene.gr)){
            gd.pos <- sum(gene.density[which(pscore>=0)])
            gd.neg <- sum(gene.density[which(pscore<0)])
            message("Gene density per bin - >0 = ", gd.pos)
            message("Gene density per bin - <0 = ", gd.neg)
            
            ## Add Chromosome compartment information
            if (gd.pos > gd.neg)
              cc <- ifelse(pscore>0, "A", "B")
            else{
              pscore <- -pscore
              cc <- ifelse(pscore>0, "A", "B")
            }
            pca.res[[eval(paste("PC",i, sep=""))]] <- GRanges(seqnames(xgi), ranges=ranges(xgi), strand=strand(xgi), score=pscore, genedens=gene.density, ccompartments=cc)
          }else{
            pca.res[[eval(paste("PC",i, sep=""))]] <- GRanges(seqnames(xgi), ranges=ranges(xgi), strand=strand(xgi), score=pscore)
          }         
        }
      }
    }
    return(pca.res)
  }

###################################
## sparseCor
##
## Correlation methods for sparse Matrix
##
## x = Matrix object
##################################

sparseCor <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x, na.rm=TRUE)
    covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
    sdvec <- sqrt(diag(covmat)) 
    cormat <- covmat/tcrossprod(sdvec)
       
    colnames(cormat) <- colnames(x)
    rownames(cormat) <- rownames(x)
    as(cormat, "Matrix")
}


getPearsonMap <- function(x, normPerExpected=TRUE, center=TRUE, ...){
  stopifnot(inherits(x,"HTCexp"))
  stopifnot(isBinned(x))

  ## force symmetric matrix
  x <- forceSymmetric(x)
  x.save <- x
 
  ## observed/expected map
  if (normPerExpected){
    x <- normPerExpected(x, ...)
  }
  
  ## Remove bins with sd=0 to avoid correlation failure
  idx.sd <- union(names(which(apply(intdata(x), c(1), sd, na.rm=TRUE)==0)),
               names(which(apply(intdata(x), c(2), sd, na.rm=TRUE)==0)))
  x <- removeIntervals(x,idx.sd)
  
  ## Remove bins with NAs
  idx.na <- union(names(which(apply(intdata(x), 1, function(x) all(is.na(x))))),
               names(which(apply(intdata(x), 2, function(x) all(is.na(x))))))
  x <- removeIntervals(x,idx.na)
  
  ## Calculate correlation
  ## use="pairwise.complete.obs" means that the correlation between two vectors is calculated using only the paired non-NAs values.
  ## It means that two correlations can be calculated using a different set of points ...
  ## In theory this never appends because bins with NAs values should be removed earlier

  ## scale on rows
  ## The Pearson's is calculated after subtracting the mean of the row from the O/E matrix.  This is because correlation of small values should be as valuable as correlation of big values.  We subtract the row mean from every entry in the row, in effect recentering the distribution.  O/E ranges from something like 1/5 to 5, and values below 1 that are correlated/anti-correlated with values above 1 need to count that way.
  if (center){
    xdata <- t(scale(t(intdata(x)), center=TRUE, scale=FALSE))
  }else{
    xdata <- as.matrix(intdata(x))
  }

  ## calculate correlation of columns
  xdata.cor <- sparseCor(xdata)

  ## reput NA
  ##idx <- setdiff(colnames(intdata()), c(idx.sd, idx.na))
  N <- dim(intdata(x.save))[1]
  cids <- colnames(intdata(x.save))
  rids <- rownames(intdata(x.save))
 
  xdata.final <- matrix(NA, ncol=N, nrow=N)
  colnames(xdata.final) <- cids
  rownames(xdata.final) <- rids
  xdata.final[rownames(xdata.cor),colnames(xdata.cor)] <- as.vector(xdata.cor)
  
  x.save@intdata <- as(xdata.final, "Matrix")
  x.save
}
