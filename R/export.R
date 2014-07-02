###################################
## exportC
##
## Export HiTC object to standard format
##
## x = an object of class HTCexp
## con = basename of files to create
##
###################################

writeC <- function(data2export, xgi, ygi, file, genome="mm9", header=TRUE){
    message("Export genomic ranges as BED files ...")

    ## Write header
    bed.xgi <- paste(file,"_xgi.bed", sep="")  
    write(paste("##HiTC - v", packageVersion("HiTC"), sep=""),file=bed.xgi)
    write(paste("##",date(), sep=""),file=bed.xgi, append=TRUE)
    rtracklayer::export(xgi, con=paste(file,"_xgi.bed", sep=""), format="bed", append=TRUE)

    if(!is.null(ygi)){
        bed.ygi <- paste(file,"_ygi.bed", sep="")  
        write(paste("##HiTC - v", packageVersion("HiTC"), sep=""),file=bed.ygi)
        write(paste("##",date(), sep=""),file=bed.ygi, append=TRUE)
        rtracklayer::export(ygi, con=paste(file,"_ygi.bed", sep=""), format="bed", append=TRUE)
    }
    count.out <- paste(file,".count", sep="")
    message("Export interaction map in '",count.out,"' ...")
    write(paste("##HiTC - v", packageVersion("HiTC"), sep=""),file=count.out)
    write(paste("##",date(), sep=""),file=count.out, append=TRUE)
    
    data2export <- as(data2export, "TsparseMatrix")
    write(paste(rownames(data2export)[data2export@i+1], colnames(data2export)[data2export@j+1], data2export@x, sep="\t"), file=count.out, append=TRUE)
    invisible(NULL)
}
    


exportC <- function(x, file, per.chromosome=FALSE){
    if (inherits(x, "HTCexp")){
        stopifnot(isSymmetric(x))
        data2export <- intdata(x)
        xgi <- x_intervals(x)
        if(!isSymmetric(x))
            ygi <- y_intervals(x)
        else
            ygi <- NULL
        writeC(data2export, xgi, ygi, file)
    }else if (inherits(x, "HTClist")){
        if (per.chromosome){
            lapply(x, function(xx){
                data2export <- intdata(xx)
                xgi <- x_intervals(xx)
                if(!isSymmetric(xx))
                    ygi <- y_intervals(xx)
                else
                    ygi <- NULL
                chrname <- paste0(seqlevels(xx), collapse="")
                writeC(data2export, xgi, ygi, file=paste0(file, "_", chrname))
            })
        }else{
            if (isComplete(x)){
                data2export <- getCombinedContacts(x)
                combi <- getCombinedIntervals(x)
                writeC(data2export, combi$xgi, combi$ygi, file=paste0(file, "_fullgen"))
            }else{
                stop("isComplete(x) is not TRUE. Cannot export incomplete HTClist object as one file. Please use per.chromosome=TRUE.")
            }
        }
    }
    invisible(NULL)
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

write.my5C <- function(data2export, xgi, ygi, file, genome="mm9", header=TRUE){

    file <- paste0(file, ".mat")
    if (header){
        write(paste0("##HiTC - v", packageVersion("HiTC")),file=file)
        write(paste0("##",date()),file=file, append=TRUE)
        write("##my5C export matrix format",file=file, append=TRUE)
    }

    xnames <- paste(id(xgi),"|",genome,"|",seqnames(xgi),":",start(xgi),"-",end(xgi), sep="")
    ynames <- paste(id(ygi),"|",genome,"|",seqnames(ygi),":",start(ygi),"-",end(ygi), sep="")
    colnames(data2export) <- xnames
    rownames(data2export) <- ynames

    message("Exporting data in '",file,"' ...")
    suppressWarnings(write.table(as.data.frame(as.matrix(data2export)), file=file, quote=FALSE, sep="\t", append=TRUE))
}
    
export.my5C <- function(x, file, genome="mm9", per.chromosome=FALSE){
    if (inherits(x, "HTCexp")){
        stopifnot(isSymmetric(x))
        data2export <- intdata(x)
        xgi <- x_intervals(x)
        ygi <- y_intervals(x)
        write.my5C(data2export, xgi, ygi, file=file)
      
    }else if (inherits(x, "HTClist")){
        if (per.chromosome){
            lapply(x, function(xx){
                data2export <- intdata(xx)
                xgi <- x_intervals(xx)
                ygi <- y_intervals(xx)
                chrname <- paste0(seqlevels(xx), collapse="")
                write.my5C(data2export, xgi, ygi, file=paste0(file, "_", chrname))
            })
        }else{
            if (isComplete(x)){
                data2export <- getCombinedContacts(x)
                combi <- getCombinedIntervals(x)
                ygi <- combi$ygi
                if (!is.null(combi$xgi))
                    xgi <- combi$xgi
                else
                    xgi <- ygi
                write.my5C(data2export, xgi, ygi, file=paste0(file, "_", genome))
            }else{
                stop("isComplete(x) is not TRUE. Cannot export incomplete HTClist object as one file. Please use per.chromosome=TRUE.")
            }
        }
    }
    invisible(NULL)
}##export.my5C
