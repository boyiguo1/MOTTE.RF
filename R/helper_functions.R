# TODO: 1) figuring out what does this calculation do.
#       2) Change the description. Description is outragously wrong now.


# Calculation of population variance
#' Population variance calculation
#'
#' @param x a numeric vaector
#'
#' @return a numeric value
#'
#' @export
#' @examples
#' pop.var(c(1,2,3))
pop.var <- function(x){
  res <- outer(x,x,FUN="-")
  res <- sum(res^2)/2/(length(x)^2)
  return(res)
}

#' Checking if elements of x belongs to a exclusive interval from min to max
#'
#' @param x The elemnt want to check. Must be an integer
#' @param min the lower bound
#' @param max the upper bound
#'
#' @return True or False
#'
#' @export
#' @examples
#' a <- 5
#' is.between(a, 1, 6)
is.between <- function(x,min,max) {
  if(min > max) {
    #return(0)
    stop("is.between::Invalid parameters: Minimum is larger than Maximum")
  }
  return((x-min)*(max-x)>0)
}


#' Title
#'
#' @param x is the matrix
#' @param epsilon some threshold
#'
#' @return
#' @export
#'
#' @examples
findMaxIndex <- function(x , epsilon){
  x <- as.matrix(x)
  thres <- abs(x[1,1])*epsilon

  index <- max(which((abs(diag(x))>thres)==TRUE))

  return(index)

}

