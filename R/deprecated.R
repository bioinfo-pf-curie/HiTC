## Normalized per zscore
setMethod("normPerZscore", signature=c("HTCexp"), definition=function(x, ...){
    .Deprecated("normPerExpected")
})

## Export Method
#setMethod("export", signature("HTCexp", "character"), function(object, con){
#    .Deprecated("exportC")
#})

##seq_name
