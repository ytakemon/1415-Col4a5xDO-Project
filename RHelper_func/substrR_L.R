## simplified version of subst(), giving one numerical input along with string object to subset. 
## This pre-specifies starting point of either left or right.

#	Subset string from the right, x = character string, n = number of characters from the right
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#	Subset string from the left, x = character string, n = number of characters from the left
substrLeft <- function(x, n){
  substr(x, 1, n)
}