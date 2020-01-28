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
#' #TODO: extract making B to a function
#' B <- matrix(
#'     c(rep(1,3), rep(0, 7),
#'     rep(0,3), rep(1,3), rep(0, 4),
#'     rep(0,6), rep(1,3), 0),
#'     nrow = 10, ncol = 3) %>%
#'     cbind(matrix(0,nrow = 10, ncol = 7))
#'
#' #TODO: extract making Z to a function
#' Z <- matrix(
#'    c(rep(1,3), rep(0, 7),
#'     rep(0,3), rep(1,3), rep(0, 4),
#'     rep(0,6), rep(1,3), 0),
#'     nrow = 10, ncol = 3
#' )
#' sim_MOTTE_data( n.train = 500, n.test = 200,
#' p = 10, q = 3, ratio = 0.5,
#' B = B, Z = Z)

# TODO: setup the magnitude for the errors

sim_MOTTE_data <- function(
  n.train = 500, # <-  500
  n.test = 200, # <- 200
  p = 10, #  <- 10
  q = 3, #  <- 3
  # pi is the ratio between treatment groups
  ratio = 0.5, # <- 0.5,
  cov.mat = diag(p),
  # TODO: incorporate the linear and polynormial for both treat.f and link.f in the code
  treat.f = c("Linear", "Polynomial"),
  link.f = c("Linear", "Polynomial"),
  B ,
  Z
  #seed = 1
){
  #set.seed(seed)
  # Simulate the binary treatment assignment for training data
  Trt.lvls <- c("Trt 1", "Trt 2")
  Trt.train <- factor(Trt.lvls[rbinom(n.train, 1, ratio)+1])

  # Simulate the autoregressive covariance matrix of X
  # TODO: replace all x.sig with cov.mat
  x.sig <- 0.8^(abs(outer(1:p,1:p,"-")))

  # Simulate x.b
  X.train.base <- MASS::mvrnorm(n.train,rep(0,p),diag(p))

  # TODO: add treatment effect in it by making treat.train in it.
  X.train.end <- X.train.base + X.train.base %*% B + mvrnorm(n.train, rep(0,p), 0.01*diag(p))

  Y.train.base <- (X.train.base)%*%Z + mvrnorm(n.train, rep(0,q), 0.01*diag(q))
  Y.train.end <- (X.train.end)%*%Z + mvrnorm(n.train, rep(0,q), 0.01*diag(q))


  ####################################
  # Simulate Testing data
  #
  ####################################

  X.test.base  <- MASS::mvrnorm(n.test,rep(0,p),diag(p))
  # With/Without treatment X.end
  # TODO: change the matrix calc to add treatment effect on top
  X.test.case.end <-  X.test.base + X.test.base %*% B
  X.test.control.end <-  X.test.base + X.test.base %*% B

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
