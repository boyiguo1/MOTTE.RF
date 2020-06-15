#' Calcualte feature importance
#'
#' @param forest
#'
#' @return
#' @export
#'
#' @examples
MOTTE_feature_importance <- function(forest){

  VI_tree <- function(root){
    if(root$isLeaf){
      return(0)
    }
    else{
      # Test cases
      if(length(root$children) > 2)
        stop("Invalid node: More than 2 children")
      if(length(root$children) == 0)
        stop("Invalid node: Inner node has no child")

      ge.node.index <- ifelse(root$children[[1]]$side,1,2)
      l.node.index <- ifelse(root$children[[1]]$side,2,1)


      return(
        tcrossprod(root$split.comb)*(root$sample_size) +
        VI_tree(root$children[[l.node.index]]) +
        VI_tree(root$children[[l.node.index]])
      )
    }
  }

  imp_mat_trees <- purrr::map(forest, VI_tree)

  imp_mat_forest <- Reduce("+", imp_mat_trees)

  eigen(imp_mat_forest)

}
