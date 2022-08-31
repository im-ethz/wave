#
# MeanRank.q
#
#----------------------------------------------------------------------
#
# PURPOSE:  create data with unique failure times and mean rank
#
# AUTHOR:  P. Heagerty
#
# DATE:  28 Feb 2011
#
#----------------------------------------------------------------------
#
MeanRank <- function( survival.time, survival.status, marker, start=NULL ){
#
# INPUT:  survival.time = follow-up time
#         survival.status = indicator of observed event
#         marker = quantitative marker
#         start = start time for interval (optional)
#
# RETURN:  time = unique event times
#          mean.rank = mean rank of case within risk set
#
#####
if( is.null(start) ) start <- rep( 0, length(survival.time) )
#####
#
##### drop any missing marker values...
#
keep <- !is.na(survival.time) & 
        !is.na(survival.status) & 
        !is.na(marker) &
        !is.na(start)
survival.time <- survival.time[ keep ]
survival.status <- survival.status[ keep ]
marker <- marker[ keep ]
start <- start[ keep ]
#
utimes <- unique( survival.time[ survival.status==1 ] )
utimes <- utimes[ order(utimes) ]
#
nonparmAUC <- rep( NA, length(utimes) )
nControls <- rep( NA, length(utimes) )
TheMarker <- marker
#
for( j in 1:length(utimes) ){
  dead.guy <- TheMarker[ (survival.time==utimes[j]) & (survival.status==1) ]
  is.started <- start < utimes[j]
  control.set <- TheMarker[ is.started & (survival.time > utimes[j]) ]
  set.size <- length(control.set)
  nControls[j] <- set.size
  ndead <- length(dead.guy)
  if( ndead==1 ){
    nonparmAUC[j] <- sum( rep(dead.guy,set.size) > control.set )/set.size
  }else{
    mean.rank <- 0
    for( k in 1:ndead ){
      this.rank <- sum( rep(dead.guy[k],set.size) > control.set )/set.size
      mean.rank <- mean.rank + this.rank/ndead
    }
    nonparmAUC[j] <- mean.rank
  }
}
out <- list( time = utimes, mean.rank = nonparmAUC, nControls=nControls )
out
}
#----------------------------------------------------------------------