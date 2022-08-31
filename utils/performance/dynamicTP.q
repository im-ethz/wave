#
# dynamicTP.q
#
#----------------------------------------------------------------------
#
# PURPOSE:  create data with unique failure times and mean rank
#
# AUTHOR:  A.Bansal, P. Heagerty
#
# DATE:  7 Oct 2015
#
#----------------------------------------------------------------------
#
##Helper functions for Nonparametric estimation

# Calculate percentile values using empirical method
pcval <- function(n_d, n_db, y_j, Y_db, t) {
#	n_db <- length(d[d==0])
#	Y_db <- sort(y[1:n_db])
#	y_j <- sort(y[(n_db+1):length(y)])

	# Calculate P(Y_db <= y_j)
	F_le <- ecdf(Y_db)

	# Calculate P(Y_db < y_j)
	y_j_neg <- -(y_j)
	Y_db_neg <- -(Y_db)
	F_l_neg <- ecdf(Y_db_neg)
	F_l <- rep(1, length(y_j)) - F_l_neg(y_j_neg)
	
	pcv_j <- F_l		# P(Y_db < y_j)

	return(pcv_j)
}

roct <- function(n_d, n_db, y_j, Y_db, t) {
	pcvals <- pcval(n_d, n_db, y_j, Y_db)
	F_pl <- ecdf(1-pcvals)
	return(F_pl(t))
}


dynamicTP <- function( p, survival.time, survival.status, marker, start=NULL ){
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
nCases <- rep( NA, length(utimes) )
TheMarker <- marker
#
for( j in 1:length(utimes) ){
  dead.guy <- TheMarker[ (survival.time==utimes[j]) & (survival.status==1) ]
  is.started <- start < utimes[j]
  control.set <- TheMarker[ is.started & (survival.time > utimes[j]) ]
  set.size <- length(control.set)
  nControls[j] <- set.size
  nCases[j] <- ndead <- length(dead.guy)
  
#  if(nControls[j] > 0)
#     nonparmAUC[j] <- roct(n_d = ndead, n_db = set.size, y_j = dead.guy, Y_db = control.set, t=p) 

  if( ndead==1 ){
    nonparmAUC[j] <- (sum( rep(dead.guy,set.size) > control.set )/set.size) > (1-p)
  }else{
    mean.rank <- 0
    for( k in 1:ndead ){
      this.rank <- (sum( rep(dead.guy[k],set.size) > control.set )/set.size) > (1-p)
      mean.rank <- mean.rank + this.rank/ndead
    } 
    nonparmAUC[j] <- mean.rank 
  }
}
out <- list( time = utimes, mean.rank = nonparmAUC, nControls=nControls, nCases=nCases )
out
}
#----------------------------------------------------------------------