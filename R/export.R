###################################
## exportC
##
## Export HiTC object to standard format
##
## x = an object of class HTCexp
## con = name of file to create
##
## chrA,startA,endA,nameA,strandA,chrB,startB,endB,nameB,strandB,countAB
###################################
## Export Method
exportC <- function(x, file){
    stopifnot(inherits(x,"HTCexp"))

    xgi <- x_intervals(x)
    ygi <- x_intervals(x)

    data <- as.data.frame(matrix(NA, ncol=11, nrow=nrow(xgi)*nrow(ygi)))
    colnames(data) <- c("chrA","startA","endA","nameA","strandA","chrB","startB","endB","nameB","strandB","countAB")

    
    data[,1:10] <- do.call("rbind", lapply(1:nrow(xgi), function(i){
        sA <- "."
        if (is.element("strand", names(annotation(xgi)))){
            sA <- annotation(xgi[i,])$strand
        }
        sB <- rep(".", nrow(ygi))
        if (is.element("strand", names(annotation(ygi)))){
            sB <- annotation(ygi)$strand
        }
        data.frame(seq_name(xgi)[i], xgi[i,], id(xgi)[i], sA, seq_name(ygi), ygi, id(ygi), sB)
    }))
    data[,11] <- as.vector(intdata(x))

    write("##HiTC data export", file=file)
    write(paste("##",date(), sep=""),file=file, append=TRUE)
    write("##chrA,startA,endA,nameA,strandA,chrB,startB,endB,nameB,strandB,countAB", file=file, append=TRUE)
    suppressWarnings(write.table(data, file=file, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE, append=TRUE))
}




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
