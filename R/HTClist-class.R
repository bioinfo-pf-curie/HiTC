################
## New S4 class to represent list of HTCexp objects
################

## class definition
setClass("HTClist",
         prototype = prototype(elementType = "HTCexp"),
         contains = "list")


setValidity("HTClist",
            function(object){
                fails <- character(0)                
                ## Check Chromosome Pairs Data
                 chrom <- sapply(object, function(x){
                     paste(seqlevels(x_intervals(x)), seqlevels(y_intervals(x)), sep="")
                 })
                 if (any(duplicated(chrom))){
                     fails <- c(fails, "Multiple 'HTCexp' found from the same chromosome pairs")
                 }

                ## Check Chromosome bounderies
                ##if (length(unique(unlist(ranges(object)))) > length(seqlevels(object))){
                ##    fails <- c(fails, "Same chromosome found with different ranges")
                ##}
                
                if (length(fails) > 0) return(fails)
                return(TRUE)
            })

## constructor
HTClist <- function(...)
{
  listData <- list(...)
  stopifnot(length(listData) > 0L)
  
  if (length(listData) == 1L && is.list(listData[[1L]]))
      listData <- listData[[1L]]
  if (!all(sapply(listData, is, "HTCexp")))
      stop("all elements in '...' must be HTCexp objects")
  listData <- unname(listData)

  names(listData) <- sapply(listData, function(x){
    paste(seqlevels(y_intervals(x)), seqlevels(x_intervals(x)), sep="")
  })
  
  new("HTClist", listData)
}


################
##
## Methods
##
################

setMethod("[", "HTClist",
    function(x, i, ...)
    {
        if (is.character(i)){
            i <- match( i, names(x))
          }        
        HTClist(unlist(x)[i])
    })          


setMethod("as.list", "HTClist",
    function(x)
    {
        as.list(unlist(x))
    })          


setMethod("c", "HTClist", function(x, ...){
    if (missing(x))
        args <- unname(list(...))
    else{
        l <- list(...)
        l <- sapply(l, function(x){
            if (inherits(x,"HTClist"))
                x <- as.list(x)
            x
        })
        args <- unname(c(as.list(x), l))
    }
    HTClist(args)
})


setMethod("detail",signature(x="HTClist"),
          function(x){
              stopifnot(validObject(x))
              lapply(x, show)
              invisible(NULL)
          })


setMethod(f="forcePairwise", signature(x="HTClist"),
          function(x){
            ##stopifnot(isComplete(x))
            ## intra sym
            #if (length(x[isIntraChrom(x)])>0){
            #  x[isIntraChrom(x)] <- HTClist(lapply(x[isIntraChrom(x)], forcePairwise))
            #}
            ## inter maps
            if (!isPairwise(x)){
              chrs <- seqlevels(x)
              pchr <- pair.chrom(chrs, rm.cis=TRUE)
              isin <- rep(0, length(pchr))
              names(isin) <- names(pchr)
              isin[names(x)] <- 1
              
              ptoadd <- pchr[names(which(isin==0))]
              nmaps <- mclapply(ptoadd, function(obj){
                symobj <- x[[paste0(obj[2], obj[1])]]
                HTCexp(intdata=t(intdata(symobj)), xgi=y_intervals(symobj), ygi=x_intervals(symobj))
              })
              
              c(x, nmaps)
            }else{
              x
            }
          })

setMethod(f="forceSymmetric", signature(x="HTClist", uplo="missing"),
          function(x, uplo){forceSymmetric(x, uplo="U")})

setMethod(f="forceSymmetric", signature(x="HTClist", uplo="character"),
          function(x, uplo){
              stopifnot(uplo=="U" || uplo=="L")
              stopifnot(isComplete(x))
              ## cis
              x[isIntraChrom(x)] <- lapply(x[isIntraChrom(x)], forceSymmetric)
              ## trans
              x <- sortSeqlevels(x)
              chrs <- seqlevels(x)
              pch <- pair.chrom(chrs, use.order=FALSE)
              if (uplo=="L"){
                  pch <- lapply(pch, "[", 2:1)
                  names(pch) <- sapply(pch, paste0, collapse="")
              }
              x[names(pch)]
          })


setMethod(f="isPairwise", signature(x="HTClist"),
          function(x){
              chrs <- seqlevels(x)
              all.pairs <- names(pair.chrom(chrs, rm.cis=TRUE))
              ret <- FALSE
              if (length(setdiff(all.pairs, names(x)))==0)
                  ret <- TRUE
              ret
          })


setMethod(f="getCombinedIntervals", signature(x="HTClist"),
          function(x, merge=FALSE){

            misintra <- length(setdiff(paste0(seqlevels(x), seqlevels(x)), names(x)))
            ## If all intrachromosomal maps are available
            if (misintra == 0){
              x <- x[isIntraChrom(x)]
              ygr <- lapply(x, y_intervals)
              ygi <- suppressWarnings(do.call(c,unname(ygr)))
              
              if(length(which(!sapply(x, isSymmetric)))>0){
                xgr <- lapply(x, x_intervals)
                xgi <- suppressWarnings(do.call(c,unname(xgr)))
              }else{
                xgi <- ygi
              }
            }else{
              ygr <- lapply(x, y_intervals)
              ygi <- suppressWarnings(do.call(c,unname(ygr)))
              ygi <- unique(ygi)

              xgr <- lapply(x, x_intervals)
              xgi <- suppressWarnings(do.call(c,unname(xgr)))
              xgi <- unique(xgi)
            }
            if (merge){
              gr <- suppressWarnings(unique(c(xgi, ygi)))
              gr <- sortSeqlevels(gr)
              gr <- sort(gr)
            }else{
              gr <- list(ygi=ygi, xgi=xgi)
            }
            gr
          })


setMethod(f="getCombinedContacts", signature(x="HTClist"),
          function(x){
              
            message("Start combining HTCexp objects ...")
            pchr <- t(as.data.frame(pair.chrom(seqlevels(x))))
            gr <- getCombinedIntervals(x, merge=TRUE)
            ## Need to use both x and y intervals for symmetric data
            dimchr <- sapply(split(gr, seqnames(gr)), length)

            col.merged <- lapply(seqlevels(x), function(chr){
              allpairs <- pchr[which(pchr[,1]==chr),]
              ##do.call(cBind,lapply(x[allpairs], function(xx){as(intdata(xx), "sparseMatrix")}})         
              do.call(cBind, apply(allpairs, 1, function(chrs){
                mapname <- paste0(chrs, collapse="")
                if (is.element(mapname, names(x))){
                    as(intdata(x[[mapname]]), "sparseMatrix")
                }else{
                  as(matrix(0, ncol=dimchr[chrs[2]], nrow=dimchr[chrs[1]]), "sparseMatrix")
                }
              }))
            })
                        
            bigMat <- do.call(rBind, col.merged)
            colnames(bigMat) <- rownames(bigMat) <- names(gr)
            message("Object size: ",object.size(bigMat))
            bigMat
          })


setMethod(f="isComplete", signature(x="HTClist"),
          function(x){
              chrs <- seqlevels(x)
              lchrs <- length(chrs)
              ret <- TRUE

              if (lchrs != length(which(isIntraChrom(x)))){
                ret <- FALSE
              }else{
                mat.pairs <- matrix(0, ncol=lchrs, nrow=lchrs, dimnames=list(chrs, chrs))
                for (m in x){
                  mat.pairs[seqlevels(m@ygi), seqlevels(m@xgi)] <- 1
                }         
                maps_per_chr <- apply(mat.pairs, 1, sum) + apply(mat.pairs, 2, sum)
                if (any(maps_per_chr<(length(seqlevels(x))+1)))
                  ret <- FALSE
              }
              ret
          })


setMethod(f="isBinned", signature(x="HTClist"),
          function(x){
              sapply(x, isBinned)
          })


setMethod(f="isIntraChrom", signature(x="HTClist"),
          function(x){
              sapply(x, isIntraChrom) 
          })


setMethod(f="ranges", signature(x="HTClist"),
          function(x){
              GRangesList(sapply(x, range)) 
          })


setMethod(f="range", signature(x="HTClist"),
          function(x){
              reduce(unlist(ranges(x)))
        })


setMethod(f="reduce", signature(x="HTClist"),
          function(x, chr, cis=TRUE, trans=TRUE, extend=FALSE){
            pc <- as.matrix(data.frame(pair.chrom(seqlevels(x), rm.cis=!cis)))
            if (extend){
                sel <- colnames(pc)[union(which(pc[1,] %in% chr), which(pc[2,] %in% chr))]
                ##sel <- colnames(pc)[which(pc[1,] %in% chr | pc[2,] %in% chr)]
            }else{
                sel <- colnames(pc)[intersect(which(pc[1,] %in% chr), which(pc[2,] %in% chr))]
                ##sel <- colnames(pc)[which(pc[1,] %in% chr & pc[2,] %in% chr)]
            }
            xr <- x[intersect(names(x), sel)]
            if (!trans){
              xr <- xr[isIntraChrom(xr)]
            }
            xr
          })


setMethod(f="seqlevels", signature(x="HTClist"),
          function(x){
              unique(unlist(as.vector(sapply(x, seqlevels))))
          })


setMethod("show",signature="HTClist",
          function(object){
              stopifnot(validObject(object))
              cat("HTClist object of length",length(object),"\n")
              if (length(object)>0){
                im <- isIntraChrom(object)
                cat(length(which(im)),"intra /", length(which(!im)), "inter-chromosomal maps\n")
              }
              invisible(NULL)
          })


setMethod("sortSeqlevels", signature="HTClist",
          function(x){
              stopifnot(validObject(x))
              chrs <- seqlevels(x)
              chrs.sort <- sortSeqlevels(chrs)
              out <- x[intersect(names(pair.chrom(chrs.sort)), names(x))]
              out
          })


setMethod("summary", signature=c(object="HTClist"),
           function(object){
               sy <- as.data.frame(t(data.frame(lapply(object, summary),stringsAsFactors=FALSE)), stringsAsFactors=FALSE)
               rownames(sy) <- paste0(sy$seq1, sy$seq2)
               sy
           })

