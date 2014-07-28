#####################
## Deprecated
#####################
setMethod("seq_name", signature=c("HTCexp"), definition=function(x, ...){
    .Deprecated("seqlevels")
})


#####################
## Defunct
#####################

## Normalized per zscore
setMethod("normPerZscore", signature=c("HTCexp"), definition=function(x, ...){
    .Defunct("normPerExpected")
})

## Export Method
setMethod("export", signature("HTCexp", "character"), function(object, con){
    .Defunct("exportC")
})


