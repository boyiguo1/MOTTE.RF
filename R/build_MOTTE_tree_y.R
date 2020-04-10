# source("global.R")

# library(CCA)              # R Library supporting CCA
# library(data.tree)        # R Library supporting the tree results

# Try everything method
# CCA diff.x and diff.y. One variable as split rule
# Currently using correlation to detect the association between x and projection of x
# left.out is ensure at least left.out*2 sample for either treated or untreated sample in the group
# left.out is used for how many treated or untreated are left out when selecting split value
# e.g. if left.out= 1 choosing max(min(treated x),min( untreated x))
# impurity has two option: "MeanTreatmentEffect", "varReduce"


#' Fitting MOTTE Tree in the MOTTE.RF
#'
#' @param x.b Pre-treatment covariates, i.e. microbiomes
#' @param x.e Post-treatment covariates, i.e. microbiomes
#' @param treat A vector of binary value, the arm of treatment
#' @param y.b Pre-treatment response, i.e. biomarkers
#' @param y.e Post-treatment response, i.e. biomarkers
# TODO: change the language, it is awkward.
#' @param nodesize An integer value. The threshold that control the maximum size of a node
#' @param nsplits A numeric value, the number of maximum splits
#' @param left.out left.out is ensure at least left.out*2 sample for either treated or untreated sample in the group
# left.out is used for how many treated or untreated are left out when selecting split value
# e.g. if left.out= 1 choosing max(min(treated x),min( untreated x))
#' @param seed a seed number to generate the random subsets of split candidates. Doesn't work when applying nsplit==NULL
#'
#' @return A data.tree object, node
#' @export
#'
#' @importFrom stats quantile var
#' @import data.tree
#'
# TODO: add import function here
# TODO: add description to setting reference trt.lvl. as R convention, the first level from levels function are use as the refence group
#       i.e. treatment 0/ treatment control
#' @examples
#' set.seed(1)
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
#' #with(train.dat,
#'     build_MOTTE_tree(x.b, x.e, factor(treat), y.b, y.e,
#'                      nodesize=30, nsplits=NULL, left.out = 0.1)
#'  #)
#'

build_MOTTE_tree_CO <- function(x.b, x.e.1, x.e.2, y.b, y.e.1, y.e.2,
                        nodesize, nsplits, left.out#) {
                        #, seed = 1
  ) {

  #set.seed(seed)

  # trt.lvl <- levels(treat)
  # if(length(trt.lvl) != 2)
  #   stop("Error Message: trt.lvl !=2 in build_MOTTE_tree")

  # Dimension
  n <- nrow(x.b)
  # n.treat.1 <- sum(treat==trt.lvl[1])
  p <- ncol(x.b)
  q <- ncol(y.b)

  ### Base cases:
  # a)
  # Only one treatment group in terminal node
  # Comment: extremely unlikely to happen
  # if(length(unique(treat))==1){
  #   return(
  #     data.tree::Node$new(
  #       paste("Terminal Node: ", n ," members", "Unique"),
  #       #xcenter = NULL,
  #       split.comb=NULL, split.value=NULL,
  #       Outcome=y.e, Treatment=treat)
  #   )}

  ### Base cases:
  # b)
  # The number of observations is smaller than threshold
  if(n <= nodesize){
    return(
      data.tree::Node$new(paste("Terminal Node: ", n," members"),
               #xcenter = NULL,
               split.comb=NULL, split.value=NULL,
               Outcome.1=y.e.1,# Treatment=treat)
               Outcome.2 = y.e.2
               )
    )}

  ### Base cases:
  # c)
  # The number of observations in each group is not larger than dimensions
  # Enforced to prevent error in CCA

  # if(min(n-n.treat.1,n.treat.1) <= max(p,q)){
    if(n <= max(p,q)){
    return(
      data.tree::Node$new(paste("Terminal Node: ", n ," members"),
               #xcenter = NULL,
               split.comb=NULL, split.value=NULL,
               Outcome.1=y.e.1,# Treatment=treat)
               Outcome.2 = y.e.2
      )
    )}

  ### Recursive cases:

  # Local level CCA using the subset of variables
  diff.x.1 <- x.e.1 - x.b
  diff.x.2 <- x.e.2 - x.b
  diff.y.1 <- y.e.1 - y.b
  diff.y.2 <- y.e.2 - y.b

  # Create the augmented matrices for both treatment arm
  # Each of the matrix have p(X^b) + p(\Delta X) + q (\Delta Y) columns
  # Trtment lvl 1

  Left.matrix <-
    rbind(
      cbind(x.b,matrix(0,nrow=n,ncol=2*p+2*q)),
      cbind(x.b,matrix(0,nrow=n,ncol=2*p+2*q)),
      cbind(matrix(0,nrow=n,ncol=3*p), diff.y.1, matrix(0, nrow=0, ncol=q)),
      cbind(matrix(0,nrow=n,ncol=3*p+q), diff.y.2),
      cbind(matrix(0,nrow=n,ncol=3*p), diff.y.1, -1*diff.y.2),
      cbind(matrix(0,nrow=n,ncol=p),diff.x.1, -1*diff.x.2, matrix(0,nrow=n,ncol=2*q)),
      cbind(matrix(0,nrow=n,ncol=p),diff.x.1, matrix(0,nrow=n,ncol=p+2*q)),
      cbind(matrix(0,nrow=n,ncol=2*p),diff.x.2, matrix(0,nrow=n,ncol=2*q))
      )

  Right.matrix <-
    rbind(
      cbind(matrix(0,nrow=n,ncol=3*p), diff.y.1, -1*diff.y.2),
      cbind(matrix(0,nrow=n,ncol=p),diff.x.1, -1*diff.x.2, matrix(0,nrow=n,ncol=2*q)),
      cbind(matrix(0,nrow=n,ncol=p),diff.x.1, matrix(0,nrow=n,ncol=p+2*q)),
      cbind(matrix(0,nrow=n,ncol=2*p),diff.x.2, matrix(0,nrow=n,ncol=2*q)),
      cbind(x.b,matrix(0,nrow=n,ncol=2*p+2*q)),
      cbind(x.b,matrix(0,nrow=n,ncol=2*p+2*q)),
      cbind(matrix(0,nrow=n,ncol=3*p), diff.y.1, matrix(0, nrow=0, ncol=q)),
      cbind(matrix(0,nrow=n,ncol=3*p+q), diff.y.2)
      )

  # Conduct CCA
  cancor.res <- CCA::cc(Left.matrix, Right.matrix)

  x.loading <- cancor.res$xcoef[1:p, 1]
  diff.y.1.loading <- cancor.res$xcoef[(3*p+1):(3*p+q),1]
  diff.y.2.loading <- cancor.res$xcoef[(3*p+q+1):(ncol(Left.matrix)),1]

  # Calculate the canonical variates
  x.proj <- x.b%*%x.loading
  diff.y.proj.1 <- diff.y.1 %*% diff.y.1.loading
  diff.y.proj.2 <- diff.y.2 %*% diff.y.2.loading

  # split.value.cand.treat1 <- unique(x.proj[treat==trt.lvl[1]])
  # split.value.cand.treat2 <- unique(x.proj[treat==trt.lvl[2]])

  # treat1.boundry <- stats::quantile(split.value.cand.treat1, c(left.out, 1-left.out))
  # treat2.boundry <- stats::quantile(split.value.cand.treat2, c(left.out, 1-left.out))

  split.value.cand <- unique(x.proj)
  # Here it used the internal funciton is.between

  if(length(split.value.cand)==0) {
    return(
      data.tree::Node$new(paste("Terminal Node: ", n," members. No split"),
                          #xcenter = NULL,
                          split.comb=NULL, split.value=NULL,
                          Outcome.1=y.e.1,# Treatment=treat)
                          Outcome.2 = y.e.2
      )
    )
  }

  # Only checking on subset of candidates
  # Reduce calculation load
  if(!is.null(nsplits)) {
    split.value.cand <- sample(split.value.cand,
                               min(nsplits, length(split.value.cand)),
                               replace=FALSE)
  }

  # Calculate impurity scores regarding variance reduction
  # impurity.score is a 2*n' matrix, the first row is the split value candidates, and second row is the
  # corresponding variance reduction
  impurity.score <- sapply(split.value.cand,FUN=function(x) {

    L.node.indices <- x.proj >= x
    R.node.indices <- x.proj < x

    # t1.indices <- treat == trt.lvl[1]
    # t2.indices <- treat == trt.lvl[2]

    L.length <- sum(L.node.indices)
    R.length <- sum(R.node.indices)
    # revision to split based on treatment difference reflected on X^b
    total.var <- (n-1)/n*var(diff.y.proj.1-diff.y.proj.2)
    left.var <- (L.length-1)/n*(var(diff.y.proj.1[L.node.indices]-diff.y.proj.2[L.node.indices]))
    right.var <- (L.length-1)/n*(var(diff.y.proj.1[R.node.indices]-diff.y.proj.2[R.node.indices]))
    return(
      matrix(
        c(x,total.var-left.var-right.var),
        nrow=2,ncol=1
      )
    )
  })

  split.value <- impurity.score[1,which.max(impurity.score[2,])]

  # Create a new node
  node <-data.tree::Node$new(
    paste("split.value = ", round(split.value, digits=3)),                        # Node name: must be unique to siblings
    #xcenter = x.center,
    split.comb=x.loading, split.value=split.value,
    Outcome.1 = NULL,# Treatment=treat)
    Outcome.2 = NULL
  )

  greater.indices <- which(x.proj>=split.value)
  less.indices <- which(x.proj<split.value)

  if(length(greater.indices)<=0) stop("greater indices <=0")
  if(length(less.indices)<=0) stop("less indices < 0")

  child.ge <- build_MOTTE_tree_CO(
    x.b = x.b[greater.indices,,drop=FALSE],
    x.e.1 = x.e.1[greater.indices,,drop=FALSE],
    x.e.2 = x.e.2[greater.indices,,drop=FALSE],
    y.b = y.b[greater.indices,,drop=FALSE],
    y.e.1 = y.e.1[greater.indices,,drop=FALSE],
    y.e.2 = y.e.2[greater.indices,,drop=FALSE],
    nodesize = nodesize, nsplits = nsplits, left.out = left.out
  )

  #
  child.l <- build_MOTTE_tree_CO(
    x.b = x.b[less.indices,,drop=FALSE],
    x.e.1 = x.e.1[less.indices,,drop=FALSE],
    x.e.2 = x.e.2[less.indices,,drop=FALSE],
    y.b = y.b[less.indices,,drop=FALSE],
    y.e.1 = y.e.1[less.indices,,drop=FALSE],
    y.e.2 = y.e.2[less.indices,,drop=FALSE],
    nodesize = nodesize, nsplits = nsplits, left.out = left.out
  )
  # side == TRUE means greater or equal than
  # side ==FALSE means less than
  # For error prevention when traverse the tree
  child.ge$side <- TRUE
  child.l$side <- FALSE

  # Rename child nodes in case of duplicate names between silblings
  child.ge$name <- paste("Greater Child: ",child.ge$name)
  child.l$name <- paste("Smaller Child: ",child.l$name)

  # Insert Children in the tree
  # TODO::add data.tree:: to somewhere
  node$AddChildNode(child.ge)
  node$AddChildNode(child.l)
  return(node)
}


#TODO: add a recursive function, which contains standardization of the data.
# Record the centers
# x.center <- attr(scale(x.b, center = T, scale = FALSE),"scaled:center")
# diff.y.0.center <- attr(scale(diff.y[treat==trt.lvl[1],,drop=FALSE],center = T, scale = FALSE),"scaled:center")
# diff.y.1.center <- attr(scale(diff.y[treat==trt.lvl[2],,drop=FALSE],center = T, scale = FALSE),"scaled:center")
