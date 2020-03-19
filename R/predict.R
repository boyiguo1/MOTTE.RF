#' Recommend treatment based on pre-treatment covariates and a weight
#'
#' @param tree.list MOTTE forest
#' @param x.b a maitrx, a matrix that contains pre-treatment covarites
#' @param w a numeric vector, indicating the direction of desired outcome
#'
#' @return a vector that contains the recommanded treatment level for each individual
#' @export
#'
#' @examples
recommendResult <- function(tree.list,x.b,w) {
  # Check root is a list of tree
  # x.b can be a matrix
  # w is a weight function matches with dimension of y.e

  if(!is.matrix(x.b))
    stop("Error Message: x.b must be a matrix")
  if(!is.matrix(w))
    w.matrix <- matrix(w, ncol=length(w), nrow=nrow(x.b), byrow=TRUE)
  else {
    w.matrix <- w
    if(nrow(w.matrix)!=nrow(x.b))
      stop("Error Message: Inconsistent dimension between x.b and w")
  }

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

#' Predicting post-treatment responses
#'
#' @param tree.list MOTTE forest
#' @param x.b a matrix, contain all observation that requires a prediction
#'
#' @return Predicted post-treatment responses for both treatment arm
#' @export
#'
#' @examples
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

#' Create predicted post-treatment responses for a single obervation
#'
#' @param tree.list MOTTE.RF object
#' @param x.b a vector of pre-treatment covariates for one observation
#'
#' @return a list contain the post-treament response mean for both treatment arms
#' @export
#'
#' @examples
recommendResult.single <- function(tree.list, x.b){
  # Check root is a list of trees
  # Check x.b is one observation
  result.list <- traverseForest(tree.list,x.b)
  out <- NULL
  treat <- NULL
  for(i in 1:length(result.list)){
    out <- rbind(out,result.list[[i]]$OUTCOME)
    treat <- c(treat,as.character(result.list[[i]]$TREAT))
  }

  treat <- as.factor(treat)
  trt.lvl <- levels(treat)

  treat1.means <- colMeans(out[treat==trt.lvl[1],,drop=F])
  treat2.means <- colMeans(out[treat==trt.lvl[2],,drop=F])
  #if(sum(is.na(untreat.means))>0) untreat.means <- rep(0,ncol(out))
  #if(sum(is.na(treat.means))>0) treat.means <- rep(0,ncol(out))
  return(list(Treat.1=treat1.means,Treat.2=treat2.means, Levels= levels(treat)))
}


#' Title
#'
#' @param forest A MOTTE.RF object
#' @param x.b a data matrix for the testing data
#'
#' @return A data matrix that contains the treatment difference
#' @export
#'
#' @examples
calcTrtDiff <- function(forest, x.b){
  apply(x.b, 1, FUN = function(x, forest)
    calcTrtDiff.single(forest = forest, x.b=x),
    forest = forest) %>%
    dplyr::bind_rows()
}

#' @importFrom magrittr "%>%"
#' @import dplyr
calcTrtDiff.single <- function(forest, x.b){
  res <- traverseForest(forest, x.b) %>%
    group_by(TREATMENT) %>%
    summarize_all(.funs = mean, na.rm=TRUE) %>%
    #TODO: Improve the treatment naming part for general function use
    ungroup %>% select(-TREATMENT)
  #TODO: this is Bad
  res[1,]-res[2,]
}

#' Traverse Forest
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
  return(purrr::map_dfr(forest,traverseTree,x.b=x.b))
}


#' Traverse Tree
#'
#' Traverse MOTTE.Tree
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

  if(root$isLeaf)
    return(data.frame(OUTCOME=root$Outcome, TREATMENT=root$Treatment))
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
