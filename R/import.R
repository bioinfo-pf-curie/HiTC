###################################
## import.my5C
## Create a HiTC object from my5C's data files (line or matrix data file).
## The object is created from the intervals (BED) definition.
## Intervals not define in the BED files are not taken into account
## if multiple chromosomes are available in the intervals files, several HTCexp object are created, one per chromosome pair.
##
## my5C.datafile : data file. Tree columns (priFOR, priREV, counts)
## xgi.bed : BED file of x intervals to consider
## ygi.bed : BED file of y intervals to consider
##
##################################

import.my5C <- function(my5C.datafile, xgi.bed, ygi.bed=NULL, all.pairwise=TRUE){

    ##Load primers
    xgi <- readBED(xgi.bed)[[1]]
    if (!is.null(ygi.bed)){
        ygi <- readBED(ygi.bed)[[1]]
    }else{
        ygi <- xgi
    }

    message("Intervals files loaded")
    my5Cdata <- read.table(my5C.datafile,comment.char = "#", check.names=FALSE)

    if (ncol(my5Cdata)==3){
        my5Cdata[,1] <- unlist(lapply(strsplit(as.character(my5Cdata[,1]),"|", fixed=TRUE),function(x){x[1]}))
        my5Cdata[,2] <- unlist(lapply(strsplit(as.character(my5Cdata[,2]),"|", fixed=TRUE),function(x){x[1]}))

        message("Convert my5C list file in HiTC object(s)")
        chromPair <- pair.chrom(chromosome(xgi), use.order=all.pairwise)
   
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- xgi[which(chromosome(xgi)==chr[1]),]
            ygi.subset <- ygi[which(chromosome(ygi)==chr[2]),]
            
            outs <- sapply(as.character(id(ygi.subset)), function(fp){
                out <- rep(0,nrow(xgi.subset))
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
            new("HTCexp", intdata, xgi.subset, ygi.subset)
        })
    } else if (ncol(my5Cdata)==dim(xgi)[1] && nrow(my5Cdata)==dim(ygi)[1]){
        message("Convert my5C matrix file in HiTC object(s)")
        rownames(my5Cdata) <- unlist(lapply(strsplit(rownames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
        colnames(my5Cdata) <- unlist(lapply(strsplit(colnames(my5Cdata),"|", fixed=TRUE),function(x){x[1]}))
  
        chromPair <- pair.chrom(chromosome(xgi), use.order=all.pairwise)
        obj <- lapply(chromPair, function(chr){
            xgi.subset <- xgi[which(chromosome(xgi)==chr[1]),]
            ygi.subset <- ygi[which(chromosome(ygi)==chr[2]),]
            
            intdata <- as.matrix(my5Cdata[as.vector(id(ygi.subset)), as.vector(id(xgi.subset))])
            new("HTCexp", intdata, xgi.subset, ygi.subset)
        })
        
    }else{
        stop("Unknown my5C input")
    }
    return(obj)
}

###################################
## pair.chrom
## INTERNAL FUNCTION
## Compute chromsome pairwise combinaison
##
## chrom : character i.e vector with chromosome names
## use.order : if TRUE, all the pairwise combinaison are returned (i.e. chr1chr2 AND chr2chr1)
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
}


###################################
## readBED
## 
## Read a BED file
## If several tracks are found, return a list of genomeIntervals objects
##
## con : The BED file
##
##################################

readBED <- function(con){
  
    BEDheader <- c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
        
    r <- readLines(con, warn=FALSE)   
    r <- r[grep("^#",r,invert=TRUE)]

    tracks <- grep("^track",r)
    if (length(tracks)>0){
        desctrack <- sapply(r[tracks], function(x){strsplit(x,split =" ")})
        tracknames <- lapply(1:length(desctrack), function(x){
            indname <- grep("^name",desctrack[[x]])
            if (length(indname)>0){
                name <- strsplit(desctrack[[x]][indname],split="=")[[1]][2]
                return(gsub("\"","",name))
            }else{
                return(paste("track",x,sep=""))
            }
        })
    }else{
        tracks <- 1
        tracknames <- "track1"
    }
    
    tracks <- c(tracks, length(r)+1)
    outl <- lapply (1:(length(tracks)-1), function(x){
        sp <- as.character("\t")
        track.start <- tracks[x]
        track.end <- tracks[x+1]-1
        beddata <- r[track.start:track.end]
        beddata <- beddata[grep("^ch",beddata)]
     
        if (length(strsplit(beddata[1], split=sp)[[1]])==1)
            sp <- as.character(" ")

        out <- t(sapply(beddata, function(x){strsplit(x,split=sp)[[1]]}))
     
        if (ncol(out)>length(BEDheader)){
            stop("Too many columns in your BED file. See http://genome.ucsc.edu/FAQ/FAQformat#format1")
        }
        if (ncol(out)<3){
            stop("Columns Chr/Start/End are needed in your BED file. See http://genome.ucsc.edu/FAQ/FAQformat#format1")
        }
        
        colnames(out)<-BEDheader[1:ncol(out)]
        rownames(out)<-NULL
        out <- as.data.frame(out)
        ##Sort data by chromosome location
        out <- out[order(out[,"chrom"],as.numeric(as.character(out[,"chromStart"])),as.numeric(as.character(out[,"chromEnd"])), decreasing=FALSE),]

        if (!is.element("strand",colnames(out))){
            out <- cbind(out,strand=rep("+",nrow(out)))
        }
        if (!is.element("name",colnames(out))){
            out <- cbind(out,name=rep("-",nrow(out)))
        }

        ##Annotation slot for gi object
        if (is.element("score",colnames(out))){
            annotation <- data.frame(seq_name=out[,"chrom"], inter_base=FALSE, strand=out[,"strand"], id=out[,"name"], score=out[,"score"])
        }else{
            annotation <- data.frame(seq_name=out[,"chrom"], inter_base=FALSE, strand=out[,"strand"], id=out[,"name"])
        }
        new("Genome_intervals", matrix(as.numeric(as.matrix(out[,c("chromStart","chromEnd")])), ncol=2, byrow=FALSE), closed=c(TRUE,TRUE), annotation=annotation)     
    })
 
    message("Reading",length(tracknames),"tracks from",basename(con))
    names(outl) <- tracknames
    
    outl
}
