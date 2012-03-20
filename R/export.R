###################################
## export.my5C
##
## Export HiTC object to my5C
##
## x = an object of class HTCexp
## outputfile = name of outputfile to create
##
###################################

export.my5C <- function(x, outputfile){

    stopifnot(inherits(x,"HTCexp"))

    write("##HiTC - export.my5C",file=outputfile)
    write(paste("##",date(), sep=""),file=outputfile, append=TRUE)

    data2export <- t(intdata(x))
    primers <- sapply(colnames(data2export), function(x){paste(x,rownames(data2export), sep="\t")})
    primers <- primers[which(data2export>0)]
    write(paste(primers,as.vector(data2export[which(data2export>0)]),sep="\t"), file=outputfile, append=TRUE)
}


## Export Method
setMethod("export", signature("HTCexp", "character"), function(object, con){
    stopifnot(inherits(object,"HTCexp"))

    if (length(grep("bed",con))==0){
        stop("File extension has to be .bed")
    }
    
    export(as(y_intervals(object),"Genome_intervals_stranded"), con=con, name="y Intervals", description="HiTC - y Intervals", color="205,102,102", nameColumn="id")
    export(as(x_intervals(object),"Genome_intervals_stranded"), con=con, name="x Intervals", description="HiTC - x Intervals", color="102,255,153", nameColumn="id", append=TRUE)

})
