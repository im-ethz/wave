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
#
#nne <- function( x, y, lambda ){
nne_simple <- function( x, y, lambda ){
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
  #var.estimate[j] <- var( y[keep] )/sum(keep)
  #Change made by Aasthaa (Oct 29, 2012):
  var.estimate[j] <- ifelse(sum(keep)==1, 0, var( y[keep] )/sum(keep))
 }
out <- list( nne = nne.estimate, var=var.estimate, x = x, lambda=lambda )
out
}
#
#
# end-of-file...

