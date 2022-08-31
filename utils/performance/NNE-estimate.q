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
#----------------------------------------------------------------------
#
nne <- function( x, y, lambda, nControls=NULL ){
#
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
  nne.estimate[j] <- mean( y[keep] )
  var.estimate[j] <- var( y[keep] )/sum(keep)
  if( !is.null(nControls) ){
#
# 08 Jan 2013
# make n0 bigger...
#
#  n0 <- mean( nControls[keep] )
   n0 <- max( nControls[keep] )
   m <- sum(keep)
#
# 08 Jan 2013
# drop "diagnonal" here...
#
   pMin <- make.pair.min( y[keep] )
   p <- nne.estimate[j]
# old...
#   var.estimate[j] <- var.estimate[j] + ( pMin - p^2)/n0
#
# new...
   add.term <- ( pMin - p^2)
   add.term <- add.term - (1/m)*( mean( y[keep] ) - mean( y[keep]^2 ) )
   var.estimate[j] <- var.estimate[j] + add.term/n0
   }
 }
out <- list( nne = nne.estimate, var=var.estimate, x = x, lambda=lambda )
out
}
#
#
# end-of-file...

