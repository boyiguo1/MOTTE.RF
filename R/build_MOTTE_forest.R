# library(parallel)
# library(doParallel)
# library(foreach)            # R library which supports parallel computing


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
#' @param left.out left.out is ensure at least left.out*2 sample for either treated or untreated sample in the group
# left.out is used for how many treated or untreated are left out when selecting split value
# e.g. if left.out= 1 choosing max(min(treated x),min( untreated x))
#'
#'
#' @return a list of data.tree. Even when only one tree construct, it is a list containing the single tree
#' @export
#'
#' @import parallel
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @examples
#' #' set.seed(1)
#' B <- create.B(10)
#' Z <- create.Z(10, 3)
#'
#' tmp.dat <- sim_MOTTE_data( n.train = 500, n.test = 200,
#' p = 10, q = 3, ratio = 0.5,
#' B = B, Z = Z)
#'
#' train.dat <- tmp.dat$train
#'
#' x.b <- scale(train.dat$x.b, center = FALSE, scale = TRUE)
#' x.e <- scale(train.dat$x.e, center = FALSE, scale = TRUE)
#' y.b <- scale(train.dat$y.b, center = FALSE, scale = TRUE)
#' y.e <- scale(train.dat$y.e, center = FALSE, scale = TRUE)
#' treat <- train.dat$trt
#'
#' # with(train.dat,
#'     build_MOTTE_forest(x.b, x.e, treat, y.b, y.e)
#' #)

build_MOTTE_forest <- function(
  x.b, x.e,
  treat,
  y.b, y.e,
  #method = "Exhaust",
  nsplits = NULL,
  nodesize = 2*(ncol(x.b)+1),
  left.out = 0.1,
  ntree = ifelse(is.null(nsplits),1,200),
  nCore=ifelse(is.null(nsplits),1,detectCores()-1)) {

  if(!is.null(nsplits))
    nsplits <- as.integer(nsplits)

  ### Error Prevention Code
  if(!is.matrix(x.b) || !is.matrix(x.e) || !is.matrix(y.b) || !is.matrix(y.e))
    stop("Error Message: x.b, x.e, y.b, y.e must be matrices")

# TODO: sort out the data structures in the future
# if(!is.vector(treat))
#    stop("Error Message: treat must be a vector")

  if(length(unique(
    nrow(x.b),nrow(x.e),length(treat),nrow(y.b), nrow(y.b)
  ))!=1)
    stop("Error Message: incosistent observation numbers in x.b, x.e, treat, y.b, y.e")

  if(!is.numeric(left.out) || !is.between(left.out,0,0.5))
    stop("Error Message: left.out must be a numeric value between 0 and 0.5")

  if(!is.null(nsplits) & !is.integer(nsplits))
    stop("Error Message: nsplits must be NULL or numeric")

  nCore <- round(nCore)
  if(!is.numeric(nCore) || nCore <= 0)
    stop("Error Message: nCore must be a positive integer")

  ntree <- round(ntree)
  if(!is.numeric(ntree) || ntree <= 0)
    stop("Error Message: ntree must be a positive integer")

  if(is.null(nsplits))
    cat("Exhaustively searching for split value\n")

  treat <- factor(treat)
  if(length(levels(treat)) != 2)
    stop("Error Message: Incorrect number of treatment groups. Must be 2 groups")

  if(nCore==1)
  {
    forest <- lapply(1:ntree,FUN = function(x){
      trt_tmp <- data.frame(treat) %>% mutate(id=1:n())

      idx <- purrr::map(levels(treat), function(lvl){
        trt_tmp %>% filter(treat==lvl) %>% pull(id) %>% sample(replace=T)
      }) %>% unlist


      return(build_MOTTE_tree(
        x.b=x.b[idx,], x.e=x.e[idx,], treat=treat[idx,], y.b=y.b[idx,], y.e=y.e[idx,],
        nodesize=nodesize, nsplits=nsplits, left.out = left.out
      )
      )
    })
  }
  else{
    ### Set up cluster enviroment
    cl <- makeCluster(min(nCore,detectCores()-1))
    registerDoParallel(cl)

    ### Construct forest with
    forest <- foreach(i = 1:ntree,
                      .combine = c,
                      .export=c("build_MOTTE_tree","is.between"),
                      .multicombine = TRUE,
                      .verbose=TRUE,
                      .packages = c("CCA","data.tree", "purrr", "dplyr"))  %dopar%
                      {

                        trt_tmp <- data.frame(treat) %>% mutate(id=1:n())

                        idx <- purrr::map(levels(treat), function(lvl){
                          trt_tmp %>% filter(treat==lvl) %>% pull(id) %>% sample(replace=T)
                        }) %>% unlist

                        tree <- build_MOTTE_tree(
                          x.b=x.b[idx,], x.e=x.e[idx,], treat=treat[idx,], y.b=y.b[idx,], y.e=y.e[idx,],
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
