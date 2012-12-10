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
    write.table(data2export, file=paste(file,".mat", sep=""), quote=FALSE, sep="\t")
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

export.my5C <- function(x, file){

    stopifnot(inherits(x,"HTCexp"))

    write("##HiTC - export.my5C",file=file)
    write(paste("##",date(), sep=""),file=file, append=TRUE)

    data2export <- t(intdata(x))
    primers <- sapply(colnames(data2export), function(x){paste(x,rownames(data2export), sep="\t")})
    primers <- primers[which(data2export>0)]
    write(paste(primers,as.vector(data2export[which(data2export>0)]),sep="\t"), file=file, append=TRUE)
}
