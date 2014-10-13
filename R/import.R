###################################
## importC
##
## Import HTCexp object to standard format
## rows/col/counts
##
## file = name of intput file to read
###################################

importC <- function(con, xgi.bed, ygi.bed=NULL, allPairwise=FALSE, forceSymmetric=FALSE){
    
    stopifnot(!missing(con))

    if(is.null(xgi.bed) && is.null(ygi.bed))
        stop("BED files of x/y intervals are required")
    
    message("Loading Genomic intervals ...")
    xgi <- rtracklayer::import(xgi.bed, format="bed", asRangedData=FALSE)
    names(xgi) <- id(xgi)
    if (!is.null(ygi.bed)){
        ygi <- rtracklayer::import(ygi.bed, format="bed", asRangedData=FALSE)
        names(ygi) <- id(ygi)
    }else{
        ygi <- xgi
    }
    
    message("Reading file ...")
    cdata <- read.table(con,comment.char = "#", colClasses=c("character","character","numeric"), check.names=FALSE)
    stopifnot(ncol(cdata)==3)

    id1 <- cdata[,1]
    id2 <- cdata[,2]
    pos1 <- match(id1, id(ygi))
    pos2 <- match(id2, id(xgi))  
    ipos1 <- which(!is.na(pos1))
    
    
    ## -1 is performed in the sparseMatrix function
    bigMat <- Matrix::sparseMatrix(i=pos1, j=pos2, x=cdata[,3], dims=c(length(ygi), length(xgi)), dimnames=list(id(ygi), id(xgi)))
    rm(cdata)
    
    message("Convert 'C' file in HTCexp object(s)")
    x <- splitCombinedContacts(bigMat, xgi, ygi, allPairwise, forceSymmetric)
}##importC


###################################
## import.my5C
##
## Create a HTCexp object from my5C's data files (matrix data file).
## Intervals not define in the BED files are not taken into account
## If multiple chromosomes are available in the intervals files, several HTCexp object are created, one per chromosome pair.
##
## my5C.datafile: data file at the matrix format
##
##################################

import.my5C <- function(file, allPairwise=FALSE, forceSymmetric=FALSE){
    
    ## Read data
    stopifnot(!missing(file))
    message("Reading file ...")
    my5Cdata <- read.table(file,comment.char = "#", check.names=FALSE)

    message("Convert my5C matrix file in HTCexp object(s)")
    my5Cdata <- as(as.matrix(my5Cdata),"Matrix")

    ## Create xgi and ygi object
    gr <- dimnames2gr(my5Cdata, pattern="\\||\\:|\\-", feat.names=c("name","org","chr","start", "end"))
    ygi <- gr[[1]]
    xgi <- gr[[2]]
    
    ## Create HTClist object from my5Cdata
    rownames(my5Cdata) <- id(ygi)
    colnames(my5Cdata) <- id(xgi)
    obj <- splitCombinedContacts(my5Cdata, xgi, ygi, allPairwise, forceSymmetric)
    
    return(HTClist(unlist(obj[which(!unlist(lapply(obj, is.null)))])))
}##import.my5C



###################################
##
## INTERNAL FUNCTION
## Split my5C row/colnames matrix to create intervals objects
## Generalized function of my5C2gr - now deprecated
## 
## x: Matrix data
## pattern: regular expression to split the colnames/rownames of x
## feat.names: features names associated with the regular expression
##
##################################

dimnames2gr <- function(x, pattern="\\||\\:|\\-", feat.names=c("name","chr","start", "end")){
    rdata <- strsplit(rownames(x), split=pattern)
    dr <- do.call(rbind.data.frame, rdata)
    stopifnot(dim(dr)[2]==length(feat.names))

    colnames(dr)<-feat.names
    rgr <- GRanges(seqnames=dr$chr, ranges = IRanges(start=as.numeric(as.character(dr$start)), end=as.numeric(as.character(dr$end)), names=as.character(dr$name)))

    if (length(setdiff(colnames(x), rownames(x))) > 0){
        cdata <- strsplit(colnames(x), pattern)
        cr <- do.call(rbind.data.frame, cdata)
        colnames(cr)<-feat.names
        cgr <- GRanges(seqnames=cr$chr, ranges = IRanges(start=as.numeric(as.character(cr$start)), end=as.numeric(as.character(cr$end)), names=as.character(cr$name)))
    }else{
        cgr <- rgr
    }
    list(rgr, cgr)
}##dimnames2gr

###################################
## splitCombinedContacts
## INTERNAL FUNCTION
## Split a genome-wide Matrix into HTClist
## Selection is done by ids from xgi and ygi objects
## 
## x: Matrix data
## xgi: GenomicRanges of x_intervals
## ygi: GenomicRanges of y_intervals
## allPairwise: see pair.chrom
## forceSymmetric: see HTCexp
##
##################################

splitCombinedContacts <- function(x, xgi, ygi, allPairwise=TRUE, forceSymmetric=FALSE){
    
    chromPair <- pair.chrom(sortSeqlevels(c(seqlevels(xgi), seqlevels(ygi))), use.order = allPairwise)
    obj <- mclapply(chromPair, function(chr) {
        ygi.subset <- ygi[which(seqnames(ygi) == chr[1]),] 
        seqlevels(ygi.subset) <- as.character(unique(seqnames(ygi.subset)))
        xgi.subset <- xgi[which(seqnames(xgi) == chr[2]),]
        seqlevels(xgi.subset) <- as.character(unique(seqnames(xgi.subset)))
     
        if (length(xgi.subset) > 0 && length(ygi.subset) > 0) {
            message("Creating ", chr[2], "-", chr[1], " Contact Map ...")
            if (length(ygi.subset)==1 || length(xgi.subset)==1){
                intdata <- Matrix(x[id(ygi.subset), id(xgi.subset)], nrow=length(ygi.subset), ncol=length(xgi.subset))
            }else{
                intdata <- x[id(ygi.subset), id(xgi.subset)]
            }
            colnames(intdata) <- id(xgi.subset)
            rownames(intdata) <- id(ygi.subset)
            HTCexp(intdata, xgi.subset, ygi.subset, forceSymmetric = forceSymmetric)
        }
    })
    ##obj
    HTClist(unlist(obj))
}##splitCombinedContacts
    
###################################
## pair.chrom
## INTERNAL FUNCTION
## Compute chromsome pairwise combinaison
##
## chrom: character i.e vector with chromosome names
## use.order: if TRUE, all the pairwise combinaison are returned (i.e. chr1chr2 AND chr2chr1)
##
##################################

pair.chrom <- function(chrom, use.order=TRUE){
    v <- unique(chrom)
    if (use.order)
      z <- unlist(sapply(1:length(v),function(i){paste(v[i], v)}))
    else
      z <- unlist(sapply(1:length(v),function(i){paste(v[i], v[i:length(v)])}))
    lz <- strsplit(z," ")
    names(lz) <- gsub(" ","",z)
    lz
}##pair.chrom
