###################################
## importC
##
## Import HTCexp object to standard format
##
## con = name of intput file to read
## chrA,startA,endA,nameA,strandA,chrB,startB,endB,nameB,strandB,countAB
###################################
## Import Method

importC <- function(con, xgi.bed=NULL, ygi.bed=NULL, all.pairwise=TRUE){
    ## Read data
    stopifnot(!missing(con))
    data <- read.table(con,comment.char = "#", check.names=FALSE)
    
    xgi <- import(xgi.bed, format="bed", asRangedData=FALSE)
    if (!is.null(ygi.bed)){
        ygi <- import(ygi.bed, format="bed", asRangedData=FALSE)
    }else{
        ygi <- xgi
    }
    message("Genomic intervals loaded")
    
    chromPair <- pair.chrom(seqlevels(GRangesList(xgi, ygi)), use.order=all.pairwise)
    obj <- lapply(chromPair, function(chr){
        xgi.subset <- xgi[which(seqnames(xgi)==chr[1]),]
        seqlevels(xgi.subset)<-as.character(unique(seqnames(xgi.subset)))
        ygi.subset <- ygi[which(seqnames(ygi)==chr[2]),]
        seqlevels(ygi.subset)<-as.character(unique(seqnames(ygi.subset)))
        
        if (length(xgi.subset)>0 && length(ygi.subset)>0){
            message("Loading ",chr[1],"-",chr[2],"...")
            intdata <- as.matrix(data[as.vector(id(ygi.subset)), as.vector(id(xgi.subset))])
            HTCexp(intdata, xgi.subset, ygi.subset)
        }
    })    
    HTClist(unlist(obj))
}##importC



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

        message("Genomic intervals loaded")
        
        my5Cdata <- read.table(my5C.datafile,comment.char = "#", check.names=FALSE)
        my5Cdata[,1] <- unlist(lapply(strsplit(as.character(my5Cdata[,1]),"|", fixed=TRUE),function(x){x[1]}))
        my5Cdata[,2] <- unlist(lapply(strsplit(as.character(my5Cdata[,2]),"|", fixed=TRUE),function(x){x[1]}))
        
        message("Convert my5C list file in HTCexp object(s)")
        chromPair <- pair.chrom(seqlevels(GRangesList(xgi,ygi)), use.order=all.pairwise)
        
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- subset(xgi, as.vector(seqnames(xgi)==chr[1]))
            seqlevels(xgi.subset)<-as.character(unique(seqnames(xgi.subset)))
            ygi.subset <- subset(ygi, as.vector(seqnames(ygi)==chr[2]))
            seqlevels(ygi.subset)<-as.character(unique(seqnames(ygi.subset)))

            
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
        
        chromPair <- pair.chrom(seqlevels(GRangesList(xgi, ygi)), use.order=all.pairwise)
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- xgi[which(seqnames(xgi)==chr[1]),]
            seqlevels(xgi.subset)<-as.character(unique(seqnames(xgi.subset)))
            ygi.subset <- ygi[which(seqnames(ygi)==chr[2]),]
            seqlevels(ygi.subset)<-as.character(unique(seqnames(ygi.subset)))

            if (length(xgi.subset)>0 && length(ygi.subset)>0){
                message("Loading ",chr[1],"-",chr[2],"...")
                intdata <- as.matrix(my5Cdata[id(ygi.subset), id(xgi.subset)])
                colnames(intdata) <- id(xgi.subset)
                rownames(intdata) <- id(ygi.subset)
                HTCexp(intdata, xgi.subset, ygi.subset)
            }
        })
    }else{
        stop("Unknown my5C input")
    }
    return(HTClist(unlist(obj[which(!unlist(lapply(obj, is.null)))])))
}##import.my5C

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
    return(list(cgr, rgr))
}##my5C2gr


    
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
