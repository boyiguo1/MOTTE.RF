#' Fitting MOTTE Tree in the MOTTE.RF
#'
#' @param x.b Pre-treatment covariates, i.e. microbiomes
#' @param x.e Post-treatment covariates, i.e. microbiomes
#' @param treat A vector of binary value, the arm of treatment
#' @param y.b Pre-treatment response, i.e. biomarkers
#' @param y.e Post-treatment response, i.e. biomarkers
#' @param nodesize An integer value. The threshold that control the maximum size of a node
#' @param nsplits A numeric value, the number of maximum splits
#' @param left.out left.out is ensure at least left.out*2 sample for either treated or untreated sample in the group
# left.out is used for how many treated or untreated are left out when selecting split value
# e.g. if left.out= 1 choosing max(min(treated x),min( untreated x))
#'
#' @return A data.tree object, node
#' @export
#'
#' @importFrom stats quantile var
#' @import data.tree
#'
#' @examples
#' set.seed(1)
#' B <- create.B(10)
#' Z <- create.Z(10, 3)
#' tmp.dat <- sim_MOTTE_data( n.train = 500, n.test = 200,
#' p = 10, q = 3, ratio = 0.5,
#' B = B, Z = Z)
#'
#' train.dat <- tmp.dat$train
#'
#' with(train.dat,
#'     build_MOTTE_tree(x.b, x.e, factor(treat), y.b, y.e,
#'                      nodesize=30, nsplits=NULL, left.out = 0.1)
#'  )
#'
#'  train.dat <- tmp.dat$train
#'
#'    x.b <- train.dat$x.b
#'    x.e <- train.dat$x.e
#'    treat <- train.dat$trt
#'    y.b <- train.dat$y.b
#'    y.e <- train.dat$y.e
#'

build_MOTTE_tree <- function(x.b, x.e, treat, y.b, y.e,
                             nodesize, nsplits, left.out) {

  trt.lvl <- levels(treat)
  if(length(trt.lvl) != 2)
    stop("Error Message: trt.lvl !=2 in build_MOTTE_tree")

  # Dimension
  n <- nrow(x.b)
  n.treat.1 <- sum(treat==trt.lvl[1])
  p <- ncol(x.b)
  q <- ncol(y.b)

  ### Base cases:
  # a)
  # Only one treatment group in terminal node
  # Comment: extremely unlikely to happen
  if(length(unique(treat))==1){
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n ," members", "Unique"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
    )}

  ### Base cases:
  # b)
  # The number of observations is smaller than threshold
  if(n <= nodesize){
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n," members"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
    )}

  ### Base cases:
  # c)
  # The number of observations in each group is not larger than dimensions
  # Enforced to prevent error in CCA

  if(min(n-n.treat.1,n.treat.1) <= max(p,q)){
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n ," members"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
    )}

  ### Recursive cases:

  # Local level CCA using the subset of variables
  treat.code <- case_when(
    treat==trt.lvl[1] ~ 1,
    treat==trt.lvl[2] ~ -1
  )

  .x.b <- scale(x.b, scale=F)
  delta.x <- x.e-x.b
  T.delta.x <- scale(treat.code*delta.x, scale=F)
  delta.y <- y.e-y.b
  T.delta.y <- scale(treat.code*delta.y, scale=F)

  # x.b.1 <- .x.b[treat==trt.lvl[1],,drop=F]
  # x.b.2 <- .x.b[treat!=trt.lvl[1],,drop=F]
  # x.e.1 <- delta.x[treat==trt.lvl[1],,drop=F] %>% scale(scale=F)
  # x.e.2 <- delta.x[treat!=trt.lvl[1],,drop=F] %>% scale(scale=F)
  # y.e.1 <- delta.y[treat==trt.lvl[1],,drop=F] %>% scale(scale=F)
  # y.e.2 <- delta.y[treat!=trt.lvl[1],,drop=F] %>% scale(scale=F)




  Left.matrix <-
    rbind(
      # cbind(x.b.1,matrix(0,nrow=n.treat.1,ncol=p+q)),
      # cbind(x.b.2,matrix(0,nrow=n-n.treat.1,ncol=p+q)),
      cbind(.x.b,matrix(0,nrow=n,ncol=p+q)),
      # cbind(x.b.1,matrix(0,nrow=n.treat.1,ncol=p+q)),
      # cbind(x.b.2,matrix(0,nrow=n-n.treat.1,ncol=p+q)),
      cbind(.x.b,matrix(0,nrow=n,ncol=p+q)),
      #cbind(matrix(0,nrow=n,ncol=p),treat.code*delta.x, matrix(0,nrow=n,ncol=q))
      cbind(matrix(0,nrow=n,ncol=p),T.delta.x, matrix(0,nrow=n,ncol=q))
    )

  Right.matrix <-
    rbind(
      # cbind(matrix(0,nrow=n.treat.1,ncol=2*p),y.e.1),
      # cbind(matrix(0,nrow=n-n.treat.1,ncol=2*p),-1*y.e.2),
      cbind(matrix(0,nrow=n,ncol=2*p),T.delta.y),
      # cbind(matrix(0,nrow=n.treat.1,ncol=p),x.e.1,matrix(0,nrow=n.treat.1,ncol=q)),
      # cbind(matrix(0,nrow=n-n.treat.1,ncol=p),-1*x.e.2,matrix(0,nrow=n-n.treat.1,ncol=q)),
      cbind(matrix(0,nrow=n,ncol=p),T.delta.x,matrix(0,nrow=n,ncol=q)),
      # cbind(matrix(0,nrow=n,ncol=2*p),treat.code*delta.y)
      cbind(matrix(0,nrow=n,ncol=2*p),T.delta.y)
    )

  # cca.res <- CCA::cc(rbind(Left.matrix, Right.matrix),
  #                    rbind(Right.matrix, Left.matrix))


  cca.res <- cancor(rbind(Left.matrix, Right.matrix),
                    rbind(Right.matrix, Left.matrix),
                    xcenter=F, ycenter=F)

  # cca.res <- stable.CCA(rbind(Left.matrix, Right.matrix),
  #                    rbind(Right.matrix, Left.matrix),
  #                   xcenter=F, ycenter=F)

  if(any(dim(cca.res$xcoef) < rep(2*p+q,2))) {
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n," members. No split"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
    )
  }

  # Calculate the canonical variates
  x.loading <- cca.res$xcoef[1:p,1]
  y.loading <- cca.res$xcoef[((2*p+1):(2*p+q)), 1]
  x.proj <- .x.b %*% x.loading
  # y.proj <- delta.y %*% y.loading
  y.proj <- T.delta.y %*% y.loading


  # Generate a vector consists of split value candidates
  #split.value.cand <- unique(x.proj)

  split.value.cand.treat1 <- unique(x.proj[treat==trt.lvl[1]])
  split.value.cand.treat2 <- unique(x.proj[treat==trt.lvl[2]])

  treat1.boundry <- stats::quantile(split.value.cand.treat1, c(left.out, 1-left.out))
  treat2.boundry <- stats::quantile(split.value.cand.treat2, c(left.out, 1-left.out))

  split.value.cand <- unique(x.proj)
  split.value.cand <- split.value.cand[is.between(unique(x.proj),
                                                  min = max(treat1.boundry[1],treat2.boundry[1]),
                                                  max = min(treat1.boundry[2],treat2.boundry[2]))]

  if(length(split.value.cand)==0) {
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n," members. No split"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
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

    L.length <- sum(L.node.indices)
    R.length <- sum(R.node.indices)

    treat.bool <- treat==trt.lvl[1]

    # # revision to split based on treatment difference reflected on X^b
    # total.var <- var(y.proj[treat.bool]) + var(y.proj[!treat.bool])
    # left.var <- var(y.proj[L.node.indices & treat.bool]) + var(y.proj[L.node.indices & (!treat.bool)])
    # right.var <- var(y.proj[R.node.indices & treat.bool]) + var(y.proj[R.node.indices & (!treat.bool)])

    total.var <- (n-1)/n*var(y.proj)
    left.var <- (L.length-1)/n*var(y.proj[L.node.indices])
    right.var <- (R.length-1)/n*var(y.proj[R.node.indices])

    return(
      matrix(
        c(x,total.var-left.var-right.var),
        nrow=2,ncol=1
      )
    )
  })

  split.value <- impurity.score[1,which.max(impurity.score[2,])]

  if(all(is.na(impurity.score[2,])))
    return(
      data.tree::Node$new(
        paste("Terminal Node: ", n," members. No split"),
        sample_size = n,
        xcenter = NULL, split.comb=NULL, split.value=NULL,
        Outcome.1=(y.e-y.b)[treat==trt.lvl[1],,drop=F],
        Outcome.2=(y.e-y.b)[treat==trt.lvl[2],,drop=F])
    )


  # Create a new node
  node <-data.tree::Node$new(
    paste("split.value = ", round(split.value, digits=3)),
    sample_size = n,# Node name: must be unique to siblings
    xcenter = attr(.x.b, "scaled:center"), #x.center,
    split.comb=x.loading, split.value=split.value,
    # Outcome=NULL, Treatment=NULL
    Outcome.1=NULL,
    Outcome.2=NULL
  )

  greater.indices <- which(x.proj>=split.value)
  less.indices <- which(x.proj<split.value)

  if(length(greater.indices)<=0) stop("greater indices <=0")
  if(length(less.indices)<=0) stop("less indices < 0")

  child.ge <- build_MOTTE_tree(
    x.b = x.b[greater.indices,,drop=FALSE], x.e = x.e[greater.indices,,drop=FALSE],
    treat = treat[greater.indices],
    y.b = y.b[greater.indices,,drop=FALSE], y.e = y.e[greater.indices,,drop=FALSE],
    nodesize = nodesize, nsplits = nsplits, left.out = left.out
  )

  #
  child.l <- build_MOTTE_tree(
    x.b = x.b[less.indices,,drop=FALSE], x.e = x.e[less.indices,,drop=FALSE],
    treat = treat[less.indices],
    y.b = y.b[less.indices,,drop=FALSE], y.e = y.e[less.indices,,drop=FALSE],
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
  node$AddChildNode(child.ge)
  node$AddChildNode(child.l)
  return(node)
}
