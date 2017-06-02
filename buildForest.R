library(parallel)
library(doParallel)
library(foreach)            # R library which supports parallel computing
source("buildTree.R")    # Tree 

#' buildForest
#' 
#' The function to fitting tree/trees
#'
#' @param x.b Before treatment covariates, a n by p matrix
#' @param x.e After treatment covariates, a n by p matrix
#' @param treat Treatment received, a n by 1 vector
#' @param y.b Before treatment outcomes, a n by q matrix
#' @param y.e After treatment outcomes, a n by q matrix
#' #@param method Method to use when choosing split value. Two options: "Exhaust" and "Random"
#' @param ntree Number of trees want to construct. By default it is 1; however, when Random method used, recommand setting it as 200
#' @param nodesize Parameter to control the node size. When the number of observations in Node smaller than nodesize, stop splitting
#' @param nsplits The number of split condidate want to examine when constructing split rule
#' @param nCore The number of cores use for forest contruction when doing parallel computation
#'
#' @return a list of data.tree. Even when only one tree construct, it is a list containing the single tree
#' @export
#'
#' @examples
buildForest <- function(
  x.b, x.e,
  treat,
  y.b, y.e,
  #method = "Exhaust",
  nsplits = NULL,
  nodesize = 2*(ncol(x.b)+1),
  left.out = 0.1,
  ntree = ifelse(is.null(nsplits),1,200), 
  nCore=ifelse(is.null(nsplits),1,detectCores()-1)) {
  
  ### Error Prevention Code
  
  if(nCore==1)
  {
    forest <- lapply(1:ntree,FUN = function(x){
      return(buildTree(
        x.b=x.b, x.e=x.e, treat=treat, y.b=y.b, y.e=y.e,
        nodesize=nodesize, nsplits=nsplits, left.out = left.out
      )
      )
    })
  }
  else{
    ### Set up cluster enviroment
    cl <- makeCluster(min(nCore,detectCores()-1) , outfile="log.out")
    registerDoParallel(cl)
    
    ### Construct forest with 
    forest <- foreach(i = 1:ntree,
                      .combine = c,
                      .export=c("buildTree","is.between"),
                      .multicombine = TRUE,
                      .verbose=TRUE,
                      .packages = c("CCA","data.tree"))  %dopar%
                      {
                        tree <- buildTree(
                          x.b=x.b, x.e=x.e, treat=treat, y.b=y.b, y.e=y.e,
                          nodesize=nodesize, nsplits=nsplits, left.out = left.out
                        )
                        print("finished building tree")
                        return(tree)
                      }
    stopCluster(cl)
    if(ntree==1) return(list(forest))
  }
  return(forest)
}