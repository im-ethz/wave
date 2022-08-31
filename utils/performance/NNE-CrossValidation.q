#
# NNE-CrossValidation.q
#
#----------------------------------------------------------------------
#
# PURPOSE: Estimate IMSE
#
# AUTHOR:  P. Heagerty
#
# DATE:  06 April 2011
#
#----------------------------------------------------------------------
#
nne.CrossValidate <- function( x, y, lambda, trim=0.05 ){
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
half.width <- ( n+1 ) * lambda / 2
#
cv.estimate <- rep( NA, n )
#
for( j in 1:n ){
  iLower <- j - half.width
  iUpper <- j + half.width
  keep <- ( c(1:n) >= iLower ) & (  c(1:n) <= iUpper )
  keep <- keep & ( c(1:n) != j )
  cv.estimate[j] <- mean( y[keep] )
 }
if( !is.na(trim) ){
  keep <- ( c(1:n) >= trim*n ) & ( c(1:n) <= (1-trim)*n )
  IMSE <- mean( (cv.estimate[keep]-y[keep])^2 )
}else{
  IMSE <- mean( (cv.estimate-y)^2 )
}
out <- list( IMSE = IMSE, lambda=lambda )
out
}
#
#
# end-of-file...

