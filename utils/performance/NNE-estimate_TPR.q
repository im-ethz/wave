#
# NNE-estimate.q
#
#----------------------------------------------------------------------
#
# PURPOSE:  compute NNE estimator
#
# AUTHOR:  P. Heagerty
#
# DATE:  06 April 2011
#
#----------------------------------------------------------------------
make.pair.min <- function( y ){
  ooo <- order(y)
  y <- y[ooo]
  n <- length(y)
  weight <- (n - c(1:n) + 1) + (n - c(1:n))
  weight <- weight/sum(weight)
  mu <- sum( y * weight )
  mu
}
make.slow <- function( y ){
  n <- length( y )
  y1 <- rep( y, n )
  y2 <- rep( y, rep(n,n) )
  mu <- mean( apply( cbind( y1, y2 ), 1, min ) )
  mu
}

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

#----------------------------------------------------------------------
#
nne_TPR <- function( x, y, lambda, nControls=NULL, nCases=NULL, 
   p, survival.time, survival.status, marker, start=NULL ){
#
if( is.null(start) ) start <- rep( 0, length(survival.time) )

# INPUT:  x      = vector
#         lambda = half-bandwidth
#
### drop any missing
#
bad.egg <- is.na(x) | is.na(y)
#
x <- x[ !bad.egg ]
y <- y[ !bad.egg ]
#
keep <- !is.na(survival.time) & 
        !is.na(survival.status) & 
        !is.na(marker) &
        !is.na(start)
survival.time <- survival.time[ keep ]
survival.status <- survival.status[ keep ]
marker <- marker[ keep ]
start <- start[ keep ]

##### order observations
#
ooo <- order(x)
#
x <- x[ ooo ]
y <- y[ ooo ]
#
n <- length(x)
#
half.width <- n * lambda / 2
#
nne.estimate <- rep( NA, n )
var.estimate <- rep( NA, n )
#
for( j in 1:n ){
  iLower <- j - half.width
  iUpper <- j + half.width
  keep <- ( c(1:n) >= iLower ) & (  c(1:n) <= iUpper )
  keep <- which(keep==TRUE)
  nne.estimate[j] <- mean( y[keep] )
  
  dead.guy <- NULL
  for(k in 1:length(keep)) {
     currTime <- x[(keep[k])]    
     dead.guy <- c(dead.guy, marker[ (survival.time==currTime) & (survival.status==1) ])
  }
  windowMaxTime <- max(x[keep])
  is.started <- start < windowMaxTime
  control.set <- marker[ is.started & (survival.time > windowMaxTime) ]
  set.size <- length(control.set)
#  nControls[j] <- set.size
#  nCases[j] <- ndead <- length(dead.guy)
  ndead <- length(dead.guy)

  if( ndead==1 ){
    roc <- (sum( rep(dead.guy,set.size) > control.set )/set.size) > (1-p)
  }else{
    mean.rank <- 0
    for( k in 1:ndead ){
      this.rank <- (sum( rep(dead.guy[k],set.size) > control.set )/set.size) > (1-p)
      mean.rank <- mean.rank + this.rank/ndead
    } 
    roc <- mean.rank 
  }

  # if(nControls[j] > 0) {
     # roc <- roct(n_d = ndead, n_db = set.size, y_j = dead.guy, Y_db = control.set, t=p) 
     var.estimate[j] <- (1 - roc)*(roc)/ndead
  #}
 }
out <- list( nne = nne.estimate, var=var.estimate, x = x, lambda=lambda )
out
}
#
#
# end-of-file...

