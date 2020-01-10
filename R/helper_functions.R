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

#' Checking if elements of a vector belongs to a exclusive interval from min to max
#'
#' @param x A numeric vector/scalar
#' @param min the lower bound of the exclusive interval
#' @param max the upper bound of the exclusive interval
#'
#' @return True or False for valid comparison result. If \code{min} > \code{max}, the function return 0.
#'
#' #@export
#' #@examples
#' #is.between(3, 1, 6)  #TRUE
#' #is.between(7, 1, 6)  #FALSE
#' #is.between(5, 4, 1)  #0
#'
#' #is.between(c(1,3,2,7), 1, 6)
is.between <- function(x,min,max) {
  if(min > max) {
    return(0)
    #stop("is.between::Invalid parameters: Minimum is larger than Maximum")
  }
  return((x-min)*(max-x)>0)
}


#' Didn't find this function in tree ir forest yet. DOn't know if this is necessary
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

