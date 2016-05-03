## Normalized per reads number
setMethod("normPerReads", signature=c("HTCexp"), definition=function(x){
    x@intdata <- x@intdata/sum(x@intdata,na.rm=TRUE)
    x
})

## Normalized using TRANS data
setMethod("normPerTrans", signature=c("HTCexp","HTCexp","HTCexp"), definition=function(x, xtrans, ytrans, method="max"){

    ## Match the good objects
    ## TODO check whether the two objects are really the same
    ## x and trans share the same x_intervals
    if(length(intersect(id(x_intervals(x)),id(x_intervals(xtrans)))) != length(id(x_intervals(x)))){
   	stop("No match between xgi cis and trans objects")
    }   
    
    ## x and trans share the same y_intervals
    if(length(intersect(id(y_intervals(x)),id(y_intervals(ytrans)))) != length(id(y_intervals(x)))){
 	stop("No match between ygi cis and trans objects")
    }   

    ix <- intdata(x)
    iy <- intdata(ytrans)
    iz <- intdata(xtrans)

    ## Weigth matrices from trans data
    normfacRow<-apply(iy, 1, "mean")
    wY <- matrix(rep(normfacRow,ncol(ix)), ncol=ncol(ix), byrow=FALSE)
    
    normfacCol<-apply(iz, 2, "mean")
    wZ <- matrix(rep(normfacCol,nrow(ix)), nrow=nrow(ix), byrow=TRUE)

    ## Normalize cis data
    ## Method
    met.agglo <- c("mult", "sum", "max")
    method <- met.agglo[pmatch(method, met.agglo)]
    if (is.na(method)) 
        stop("Error :Unknown method ['mult','sum','max' are expected].")

    if (method == "sum"){
        w <- (wY+wZ)
    }else if (method == "mult"){
        w <- (wY*wZ)
    }else if (method == "max"){
        w <- matrix(NA, ncol=ncol(wY), nrow=nrow(wY))
        for (i in 1:nrow(wY)){w[i,]<-pmax(wY[i,], wZ[i,])}
    }
    w[which(w==0)] <- 1
    xnorm <- ix/w

    intdata(x) <- xnorm
    return(x)
})
