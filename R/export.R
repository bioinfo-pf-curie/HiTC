###################################
## exportC
##
## Export HiTC object to standard format
##
## x = an object of class HTCexp
## con = basename of files to create
##
###################################
exportC <- function(x, file="HiTCexport"){
    stopifnot(inherits(x,"HTCexp"))

    xgi <- x_intervals(x)
    ygi <- y_intervals(x)

    message("Export genomic ranges as BED files ...")
    export(xgi, con=paste(file,"_xgi.bed", sep=""), format="bed")
    export(ygi, con=paste(file,"_ygi.bed", sep=""), format="bed")

    message("Export interaction maps as matrix file ...")
    data2export <- t(intdata(x))
    write.table(as.data.frame(as.matrix(data2export)), file=paste(file,".mat", sep=""), quote=FALSE, sep="\t")
}##exportC



###################################
## export.my5C
##
## Export HiTC object to my5C
##
## x = an object of class HTCexp
## file = name of file to create
##
###################################

export.my5C <- function(x, file, format=c("mat", "list"), genome="mm9", header=TRUE){

    stopifnot(inherits(x,"HTCexp"))
    format <- match.arg(format)
    
    if (header){
        write(paste("##HiTC - v", packageVersion("HiTC"), sep=""),file=file)
        write(paste("##",date(), sep=""),file=file, append=TRUE)
    }
    
    data2export <- t(intdata(x))
    xgi <- x_intervals(x)
    xnames <- paste(id(xgi),"|",genome,"|",seqnames(xgi),":",start(xgi),"-",end(xgi), sep="")
    ygi <- y_intervals(x)
    ynames <- paste(id(ygi),"|",genome,"|",seqnames(ygi),":",start(ygi),"-",end(ygi), sep="")

    ## Inversion because of t(intdata(x))
    colnames(data2export) <- ynames
    rownames(data2export) <- xnames
    
    if (format=="list"){
        primers <- sapply(colnames(data2export), function(x){paste(x,rownames(data2export), sep="\t")})
        primers <- primers[which(data2export>0 & !is.na(data2export))]
        write(paste(primers,as.vector(data2export[which(data2export>0)]),sep="\t"), file=file, append=TRUE)
    }else{
        write.table(as.data.frame(as.matrix(data2export)), file=file, quote=FALSE, sep="\t")
    }
}##export.my5C
