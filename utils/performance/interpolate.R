#
# interpolate.q
#
#----------------------------------------------------------------------
#
# PURPOSE:  interpolation
#
# DATE:  11 May 2010
#
# AUTHOR:  P. Heagerty
#
#----------------------------------------------------------------------
interpolate <- function( x, y, target ){
  ooo <- order(x)
  x <- x[ooo]
  y <- y[ooo]
  n <- length(target)
  m <- length(x)
  iii <- c( 1:length(x) )
  output <- rep( NA, n )
  for( j in 1:n ){
    index.low <- max( iii[ x <= target[j] ] )
    index.high <- min( iii[  target[j] <= x ] )
    if( target[j] > max(x) ){
	index.low <- m
        index.high <- m
    }
    if( target[j] < min(x) ){
	  index.low <- 1
        index.high <- 1
    }
#
    if( index.low == index.high ){
	output[j] <- y[ index.low ]
     }else{
        dx <- ( target[j] - x[index.low] )/( x[index.high]-x[index.low] )
        dy <- ( y[ index.high ] - y [ index.low ] )
	output[j] <- y[ index.low ] + dx * dy
     }
}
output
}
#----------------------------------------------------------------------
