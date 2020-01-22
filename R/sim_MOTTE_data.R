#' Simulate MOTTE formated data
#'
#' @param n.train An integer value, sample size for training data
#' @param n.test An integer value, sample size for testing data
#' @param p An integer value, the dimensionality of covaraites (X)
#' @param q An integer value, the dimensionality of responses (Y)
#' @param pi An fraction value, the ratio between treatment groups
#' #@param seed An integer value, set the seed number
#'
#' @return A nested list that contain training data and testing data. TODO: add more introdcution to whats in the list
#' @export
#' @importFrom MASS mvrnorm
#'
#' @examples
#' set.seed(1)
#' sim_MOTTE_data( n.train = 500, n.test = 200,
#' p = 10, q = 3, pi = 0.5)


sim_MOTTE_data <- function(
  n.train = 500, # <-  500
  n.test = 200, # <- 200
  p = 10, #  <- 10
  q = 3, #  <- 3
  # pi is the ratio between treatment groups
  ratio = 0.5#, # <- 0.5,
  #seed = 1
){
  #set.seed(seed)
  # Simulate the binary treatment assignment for training data
  Treat.train <- rbinom(n.train, 1, ratio)

  # Set up for the design matrix
  # TODO: extract this as an argument
  # figure out how this has been in the paper
  Z <- matrix(
    cbind(
      c(rep(1,5),rep(0,5)),
      c(rep(0,2),rep(1,5),rep(0,3)),
      c(rep(0,5),rep(1,5))
    ),nrow=p,ncol=q)

  # Simulate the autoregressive covariance matrix of X
  x.sig <- 0.8^(abs(outer(1:p,1:p,"-")))

  # Simulate x.b
  X.train.base <- MASS::mvrnorm(n.train,rep(0,p),diag(p))

  # TODO: why the treatment effect have the same i, and j
  # TODO: Extract treat.effect as a parameter
  X.train.end <- X.train.base
  X.train.end[,5] <- X.train.end[,5] + treat.effect(X.train.base,Treat.train,1,3)
  X.train.end[,6] <- X.train.end[,6] + treat.effect(X.train.base,Treat.train,1,3)
  X.train.end[,7] <- X.train.end[,7] + treat.effect(X.train.base,Treat.train,1,3)
  X.train.end <- X.train.end + mvrnorm(n.train, rep(0,p), 0.01*diag(p))

  #TODO: fix the sd here. Too small
  Y.train.base <- (X.train.base)%*%Z + mvrnorm(n.train, rep(0,q), 0.01*diag(q))
  Y.train.end <- (X.train.end)%*%Z + mvrnorm(n.train, rep(0,q), 0.01*diag(q))


  ####################################
  # Simulate Testing data
  #
  ####################################

  X.test.base  <- MASS::mvrnorm(n.test,rep(0,p),diag(p))
  # With/Without treatment X.end
  X.test.case.end <-  X.test.base
  X.test.control.end <- X.test.base
  X.test.case.end[,5] <- X.test.case.end[,5] + treat.effect(X.test.base,rep(1,n.test),1,3)
  X.test.case.end[,6] <- X.test.case.end[,6] + treat.effect(X.test.base,rep(1,n.test),1,3)
  X.test.case.end[,7] <- X.test.case.end[,7] + treat.effect(X.test.base,rep(1,n.test),1,3)
  X.test.control.end[,5] <- X.test.control.end[,5] + treat.effect(X.test.base,rep(0,n.test),1,3)
  X.test.control.end[,6] <- X.test.control.end[,6] + treat.effect(X.test.base,rep(0,n.test),1,3)
  X.test.control.end[,7] <- X.test.control.end[,7] + treat.effect(X.test.base,rep(0,n.test),1,3)
  X.test.case.end <- X.test.case.end
  X.test.control.end <- X.test.control.end
  # With/Without treatment Y.base
  Y.test.case.end <- (X.test.case.end)%*%Z
  Y.test.control.end <- (X.test.control.end)%*%Z

  return(
    list(
      train = list(x.b = X.train.base,
                  x.e = X.train.end,
                  treat = Treat.train,
                  y.b = Y.train.base,
                  y.e = Y.train.end),
      test = list(
        x.b = X.test.base,
        y.e.case = Y.test.case.end,
        y.e.control = Y.test.control.end
      )
    )
  )
}

treat.effect <- function(x.base,treat,i,j)
{
  # Create treatment effective population
  treat1.pop <- x.base[,i]+x.base[,j]>0
  treat0.pop <- rep(TRUE,nrow(x.base))

  # Create treatment effective magnitude
  no.treat.effect <- 0
  treat1.effect <- ifelse(treat1.pop,6,no.treat.effect)
  treat0.effect <- ifelse(treat0.pop,3,no.treat.effect)

  return(ifelse(treat,treat1.effect,treat0.effect))
}
