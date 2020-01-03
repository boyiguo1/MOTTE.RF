library(CCA)              # R Library supporting CCA
library(data.tree)        # R Library supporting the tree results

# Calculation of population variance
#' Population mean calculation
#' 
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
pop.var <- function(x){
  res <- outer(x,x,FUN="-")
  res <- sum(res^2)/2/(length(x)^2)
  return(res)
}

#' Checking if elements of x belongs to a exclusive interval from min to max
#'
#' @param x The elemnt want to check. Must be an integer
#' @param min
#' @param max 
#'
#' @return True or False
#' @export
#'
#' @examples
is.between <- function(x,min,max) {
  if(min > max)
    stop("Invalid parameters: Minimum is larger than Maximum")
  
  return((x-min)*(max-x)>0)
}
