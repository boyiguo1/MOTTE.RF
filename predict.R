recommendResult <- function(tree.list,x.b,w) {
  # Check root is a list of tree
  # x.b can be a matrix
  # w is a weight function matches with dimension of y.e

  if(!is.matrix(x.b))
    stop("Error Message: x.b must be a matrix")

  return(
    apply(x.b,1,FUN=function(x,tree.list) {
      comp <- recommendResult.single(tree.list,x)
      return(
        ifelse(w%*%comp$Treat.1 > w%*%comp$Treat.2,
               comp$Levels[1], comp$Levels[2])
      )
    },
    tree.list=tree.list)
  )
}

predictResult <- function(tree.list, x.b){
  
  if(!is.matrix(x.b))
    stop("Error Message: x.b must be a matrix")
  
  return(
    apply(x.b,1,FUN=function(x,tree.list){
      recommendResult.single(tree.list,x)
    },
    tree.list=tree.list)
  )
}

recommendResult.single <- function(tree.list, x.b){
  # Check root is a list of trees
  # Check x.b is one observation
  result.list <- traverseForest(tree.list,x.b)
  out <- NULL
  treat <- NULL
  for(i in 1:length(result.list)){
    out <- rbind(out,result.list[[i]]$OUTCOME)
    treat <- c(treat,result.list[[i]]$TREAT)
  }
  
  trt.lvl <- levels(treat)
  
  treat1.means <- colMeans(out[treat==trt.lvl[1],,drop=F])
  treat2.means <- colMeans(out[treat==trt.lvl[2],,drop=F])
  #if(sum(is.na(untreat.means))>0) untreat.means <- rep(0,ncol(out))
  #if(sum(is.na(treat.means))>0) treat.means <- rep(0,ncol(out)) 
  return(list(Treat.1=treat1.means,Treat.2=treat2.means, Levels= levels(treat)))
}


#' Title
#'  A wrapper function using traverseTree to traverse the forest
#'
#' @param forest A list of trees
#' @param x.b a single row vector contains the baseline variates of one observation
#'
#' @return A list of list containing OUTCOME and TREATMENT for each tree 
#' @export
#'
#' @examples
traverseForest <- function(forest, x.b) {
  # Check root is a list of trees
  # Check x.b is one observation
  return(lapply(forest,traverseTree,x.b=x.b))
}


#' Title
#'
#' @param root data.tree
#' @param x.b a single row vector contains the baseline variates of one observation
#'
#' @return A list contain two matrices: OUTCOME and TREATMENT
#' @export
#'
#' @examples
traverseTree <- function(root, x.b){
  # Check the validity of root
  
  # Check x.b is a vector instead of a matrix
  
  x.b <- matrix(x.b, nrow=1)
  
  split.comb <- root$split.comb
  split.value <- root$split.value
  
  if(root$isLeaf){
    tmp <- 1
    return(list(OUTCOME=root$Outcome, TREATMENT=root$Treatment))}
  else{
    # Test cases
    if(length(root$children) > 2)
      stop("Invalid node: More than 2 children")
    if(length(root$children) == 0)
      stop("Invalid node: Inner node has no child")
    
    ge.node.index <- ifelse(root$children[[1]]$side,1,2)
    l.node.index <- ifelse(root$children[[1]]$side,2,1)
    if((x.b-root$xcenter)%*%split.comb>=split.value)
      return(traverseTree(root$children[[ge.node.index]],x.b))
    else
      return(traverseTree(root$children[[l.node.index]],x.b))
  }
}