#####################
## Deprecated
#####################
#setMethod("forcePairwise", signature=c("HTCexp"), definition=function(x){
#    .Deprecated("forcePairwise")
#})


#####################
## Defunct
#####################

##seq_name
setMethod("seq_name", signature=c("HTCexp"), definition=function(x, ...){
    .Defunct("seqlevels")
})

## Normalized per zscore
setMethod("normPerZscore", signature=c("HTCexp"), definition=function(x, ...){
    .Defunct("normPerExpected")
})

## Export Method
setMethod("export", signature("HTCexp", "character"), function(object, con){
    .Defunct("exportC")
})


