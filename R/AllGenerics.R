##Generics
setGeneric(name="intdata", def=function(x) standardGeneric("intdata"))
setGeneric(name="intdata<-", def=function(x,value) standardGeneric("intdata<-"))

setGeneric(name="y_intervals", def=function(x) standardGeneric("y_intervals"))
setGeneric(name="y_intervals<-", def=function(x, value) standardGeneric("y_intervals<-"))

setGeneric(name="x_intervals", def=function(x) standardGeneric("x_intervals"))
setGeneric(name="x_intervals<-", def=function(x, value) standardGeneric("x_intervals<-"))

setGeneric(name="isBinned", def=function(x) standardGeneric("isBinned"))
setGeneric(name="isIntraChrom", def=function(x) standardGeneric("isIntraChrom"))

setGeneric(name="divide", def=function(x,y) standardGeneric("divide"))
setGeneric(name="substract", def=function(x,y) standardGeneric("substract"))

setGeneric(name="normPerReads", def=function(x) standardGeneric("normPerReads"))
setGeneric(name="normPerExpected", def=function(x, ...) standardGeneric("normPerExpected"))
setGeneric(name="normPerTrans", def=function(x, xtrans, ytrans, ...) standardGeneric("normPerTrans"))

