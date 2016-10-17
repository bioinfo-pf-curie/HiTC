##Generics

## Data access

setGeneric(name="intdata", def=function(x) standardGeneric("intdata"))
setGeneric(name="intdata<-", def=function(x,value) standardGeneric("intdata<-"))

setGeneric(name="id", def=function(x) standardGeneric("id"))
##setGeneric(name="id<-", def=function(x,value) standardGeneric("id<-"))

setGeneric(name="y_intervals", def=function(x) standardGeneric("y_intervals"))
setGeneric(name="y_intervals<-", def=function(x, value) standardGeneric("y_intervals<-"))
setGeneric(name="x_intervals", def=function(x) standardGeneric("x_intervals"))
setGeneric(name="x_intervals<-", def=function(x, value) standardGeneric("x_intervals<-"))
setGeneric(name="xy_intervals", def=function(x) standardGeneric("xy_intervals"))

## Update objects
setGeneric(name="getCombinedIntervals", def=function(x, ...) standardGeneric("getCombinedIntervals"))
setGeneric(name="getCombinedContacts", def=function(x) standardGeneric("getCombinedContacts"))
setGeneric(name="forcePairwise", def=function(x) standardGeneric("forcePairwise"))
setGeneric(name="forceTriangular", def=function(x) standardGeneric("forceTriangular"))

## is 
setGeneric(name="isComplete", def=function(x) standardGeneric("isComplete"))
setGeneric(name="isPairwise", def=function(x) standardGeneric("isPairwise"))

setGeneric(name="isBinned", def=function(x) standardGeneric("isBinned"))
setGeneric(name="isIntraChrom", def=function(x) standardGeneric("isIntraChrom"))

## Operation
#setGeneric(name="directionalityIndex", def=function(x) standardGeneric("directionalityIndex"))

setGeneric(name="divide", def=function(x,y) standardGeneric("divide"))
setGeneric(name="substract", def=function(x,y) standardGeneric("substract"))

setGeneric(name="normPerReads", def=function(x) standardGeneric("normPerReads"))
setGeneric(name="normPerExpected", def=function(x, ...) standardGeneric("normPerExpected"))
setGeneric(name="normPerTrans", def=function(x, xtrans, ytrans, ...) standardGeneric("normPerTrans"))

setGeneric(name="mapC", def=function(x, ...) standardGeneric("mapC"))
setGeneric(name="mapC", def=function(x, y, ...) standardGeneric("mapC"))

## Deprecated/Defunct
setGeneric(name="normPerZscore", def=function(x, ...) standardGeneric("normPerZscore"))
setGeneric(name="seq_name", def=function(x,...) standardGeneric("seq_name"))

## standard generics
setGeneric("summary")
