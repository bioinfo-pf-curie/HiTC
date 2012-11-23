###################################
## importC
##
## Import HTCexp object to standard format
##
## con = name of outputfile to create
## chrA,startA,endA,nameA,strandA,chrB,startB,endB,nameB,strandB,countAB
###################################
## Import Method
importC <- function(con, all.pairwise=TRUE){

    alldata<-read.table(con, comment.char = "#", check.names=FALSE, sep=",")

    rgr <- alldata[,1:5]
    rgr <- rgr[which(!duplicated(rgr)),]
    rownames(rgr) <- 1:nrow(rgr)
    cgr <- alldata[,6:10]
    cgr <- cgr[which(!duplicated(cgr)),]
    rownames(cgr) <- 1:nrow(cgr)
    data <- alldata[,c(4,9,11)]

    ## Creating GRanges objects
    ygi <- GRanges(seqnames=rgr[,1], ranges = IRanges(start=rgr[,2], end=rgr[,3], names=rgr[,4]), strand =rgr[,5])
    xgi <- GRanges(seqnames=cgr[,1], ranges = IRanges(start=cgr[,2], end=cgr[,3], names=cgr[,4]), strand =cgr[,5])

    chromPair <- pair.chrom(chromosome(c(xgi,ygi)), use.order=all.pairwise)
    
    obj <- lapply(chromPair, function(chr){
        xgi.subset <- subset(xgi, as.vector(seqnames(xgi)==chr[1]))
        ygi.subset <- subset(ygi, as.vector(seqnames(ygi)==chr[2]))
  
      outs <- sapply(id(ygi.subset), function(fp){
                 out <- rep(0,length(xgi.subset))
                 subdata <- data[which(data[,1]==fp),]
        
                 if(dim(subdata)[1] == 0)
                     warning(fp," not found in the interaction data",call.=FALSE, immediate.=TRUE)
                 else {
                     ind.subdata <- which(subdata[,2]%in%id(xgi.subset))
                     ind.intdata <- match(subdata[,2],id(xgi.subset))
                     out[ind.intdata[which(!is.na(ind.intdata))]] <- subdata[ind.subdata,3]
                 }
                 out
      }, USE.NAMES=TRUE)
      intdata <- as.matrix(t(outs))
      colnames(intdata) <- id(xgi.subset)
      HTCexp(intdata, xgi.subset, ygi.subset)
    })
   return(obj[which(!unlist(lapply(obj, is.null)))])
}##ImportC TO VALIDATE



###################################
## import.my5C
##
## Create a HTCexp object from my5C's data files (line or matrix data file).
## The object is created from the intervals (BED) definition.
## Intervals not define in the BED files are not taken into account
## if multiple chromosomes are available in the intervals files, several HTCexp object are created, one per chromosome pair.
##
## my5C.datafile: data file. Tree columns (priFOR, priREV, counts)
## xgi.bed: BED file of x intervals to consider
## ygi.bed: BED file of y intervals to consider
##
##################################

import.my5C <- function(my5C.datafile, xgi.bed=NULL, ygi.bed=NULL, all.pairwise=TRUE){
    
    ## Read data
    stopifnot(!missing(my5C.datafile))
    my5Cdata <- read.table(my5C.datafile,comment.char = "#", check.names=FALSE)
    
    if (ncol(my5Cdata)==3){
        if(is.null(xgi.bed) || is.null(ygi.bed))
            stop("BED files of x/y intervals are required for my5C list format")
        
        xgi <- import(xgi.bed, format="bed", asRangedData=FALSE)
        if (!is.null(ygi.bed)){
            ygi <- import(ygi.bed, format="bed", asRangedData=FALSE)
        }else{
            ygi <- xgi
        }
        seqlevels(xgi) <- unique(c(seqlevels(xgi), seqlevels(ygi)))
        seqlevels(ygi) <- unique(c(seqlevels(xgi), seqlevels(ygi)))
        message("Genomic intervals loaded")
        
        my5Cdata <- read.table(my5C.datafile,comment.char = "#", check.names=FALSE)
        my5Cdata[,1] <- unlist(lapply(strsplit(as.character(my5Cdata[,1]),"|", fixed=TRUE),function(x){x[1]}))
        my5Cdata[,2] <- unlist(lapply(strsplit(as.character(my5Cdata[,2]),"|", fixed=TRUE),function(x){x[1]}))
        
        message("Convert my5C list file in HTCexp object(s)")
        chromPair <- pair.chrom(chromosome(c(xgi,ygi)), use.order=all.pairwise)
        
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- subset(xgi, as.vector(seqnames(xgi)==chr[1]))
            ygi.subset <- subset(ygi, as.vector(seqnames(ygi)==chr[2]))
            
            if (length(xgi.subset)>0 && length(ygi.subset)>0){
                outs <- sapply(id(ygi.subset), function(fp){
                    out <- rep(0,length(xgi.subset))
                    subdata <- my5Cdata[which(my5Cdata[,1]==fp),]
                    
                    if(dim(subdata)[1] == 0)
                        warning(fp," not found in the interaction data",call.=FALSE, immediate.=TRUE)
                    else {
                        ind.subdata <- which(subdata[,2]%in%id(xgi.subset))
                        ind.intdata <- match(subdata[,2],id(xgi.subset))
                        out[ind.intdata[which(!is.na(ind.intdata))]] <- subdata[ind.subdata,3]
                    }
                    out
                }, USE.NAMES=TRUE)
                intdata <- as.matrix(t(outs))
                colnames(intdata) <- as.character(id(xgi.subset))
                HTCexp(intdata, xgi.subset, ygi.subset)
            }
        })
    } else if (is.data.frame(my5Cdata)){
        message("Convert my5C matrix file in HTCexp object(s)")
        ## Create xgi and ygi object
        gr <- my5C2gr(my5Cdata)
        xgi <- gr[[1]]
        ygi <- gr[[2]]
        
        ## Create objects
        rownames(my5Cdata) <- unlist(lapply(strsplit(rownames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
        colnames(my5Cdata) <- unlist(lapply(strsplit(colnames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
        
        chromPair <- pair.chrom(chromosome(c(xgi, ygi)), use.order=all.pairwise)
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- xgi[which(seqnames(xgi)==chr[1]),]
            ygi.subset <- ygi[which(seqnames(ygi)==chr[2]),]
            
            if (length(xgi.subset)>0 && length(ygi.subset)>0){
                message("Loading ",chr[1],"-",chr[2],"...")
                intdata <- as.matrix(my5Cdata[as.vector(id(ygi.subset)), as.vector(id(xgi.subset))])
                HTCexp(intdata, xgi.subset, ygi.subset)
            }
        })
    }else{
        stop("Unknown my5C input")
    }
    return(obj[which(!unlist(lapply(obj, is.null)))])
} ##import.my5C

###################################
## pair.chrom
## INTERNAL FUNCTION
## Split my5C row/colnames matrix to create intervals objects
##
## my5Cdata: my5C matrix data
##
##################################

my5C2gr <- function(my5Cdata){
    if (is.null(rownames(my5Cdata)) || is.null(colnames(my5Cdata)))
        stop("rownames and colnames in the my5C format are required")
    
    ## Split row/colnames
    rdata <- matrix(unlist(lapply(strsplit(rownames(my5Cdata),"|", fixed=TRUE),function(x){
        tmp<-unlist(strsplit(x[3],":", fixed=TRUE));
        tmp2<-unlist(strsplit(tmp[2],"-", fixed=TRUE));
        return(c(tmp[1], tmp2[1], tmp2[2]))
    })), ncol=3, byrow=TRUE)
    rname <- unlist(lapply(strsplit(rownames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
    
    cdata <- matrix(unlist(lapply(strsplit(colnames(my5Cdata),"|", fixed=TRUE),function(x){
        tmp<-unlist(strsplit(x[3],":", fixed=TRUE));
        tmp2<-unlist(strsplit(tmp[2],"-", fixed=TRUE));
        return(c(tmp[1], tmp2[1], tmp2[2]))
    })), ncol=3, byrow=TRUE)
    cname <- unlist(lapply(strsplit(colnames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
    
    ## Create GRanges objects
    rgr <- GRanges(seqnames=rdata[,1], ranges = IRanges(start=as.numeric(rdata[,2]), end=as.numeric(rdata[,3]), names=rname))
    cgr <- GRanges(seqnames=cdata[,1], ranges = IRanges(start=as.numeric(cdata[,2]), end=as.numeric(cdata[,3]), names=cname))
    seqlevels(rgr) <- unique(c(seqlevels(rgr), seqlevels(cgr)))
    seqlevels(cgr) <- unique(c(seqlevels(rgr), seqlevels(cgr)))
 
    return(list(cgr, rgr))
}


    
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
