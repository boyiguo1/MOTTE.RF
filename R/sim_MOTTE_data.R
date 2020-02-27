#' Simulate MOTTE formated data
#'
#' @param n.train An integer value, sample size for training data
#' @param n.test An integer value, sample size for testing data
#' @param p An integer value, the dimensionality of covaraites (X)
#' @param q An integer value, the dimensionality of responses (Y)
#' @param ratio An fraction value, the ratio between treatment groups
#' @param cov.mat Covariance matrix for X.b
#' @param trt.f character, the name of treatment effect funciton
#' @param link.f character, the name of treatment effect funciton
#' @param B matrix, the coefficient matrix used in the trt.f
#' @param Z matrix, the coefficient matrix used in the link.f
#'
#' @return A nested list that contain training data and testing data. TODO: add more introdcution to whats in the list
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom
#'
#' @examples
#' set.seed(1)
#' B <- create.B(10)
#' Z <- create.Z(10, 3)
#'
#' sim_MOTTE_data( n.train = 500, n.test = 200,
#' p = 10, q = 3, ratio = 0.5,
#' B = B, Z = Z)
#'


# n.train = 200
# n.test = 100
# p = 10
# q = 3
# ratio = 0.5
# cov.mat = diag(p)
# trt.f = "Polynomial"
# link.f = "Polynomial"
# B = create.B(p)
# Z = create.Z(p,q)

# TODO: setup the magnitude for the errors
sim_MOTTE_data <- function(
  n.train = 500, # <-  500
  n.test = 200, # <- 200
  p = 10, #  <- 10
  q = 3, #  <- 3
  ratio = 0.5, # <- 0.5,
  # TODO: set up error prevention for cov.mat, i.e. measuring the dimension
  cov.mat = diag(p),
  # TODO: incorporate the linear and polynormial for both treat.f and link.f in the code
  trt.f = c("Linear", "Polynomial"),
  link.f = c("Linear", "Polynomial"),
  B ,
  Z
  #seed = 1
){

  # TODO: extract this as argument for the parameter
  c.x <- sum((1:3)^2)/2
  c.y <- sum((1:3)^2)/2

  trt.f <- trt.f[[1]]
  link.f <- link.f[[1]]

  # TODO: Improve the warning language
  .link.f <- switch(link.f,
                    "Linear" = function(x){x%*%Z},
                    "Polynomial" = function(x){(x^2)%*%Z},
                    stop("Link function doesn't exist, choose from 'Linear' or 'Polynomial'"))
  .trt.f <- switch(trt.f,
                   "Linear" = function(x, trt){sweep(x, 1, trt, "*")%*%B},
                   "Polynomial" = function(x, trt){(sweep(x^2, 1, trt, "*"))%*%B},
                   stop("Trt.f doesn't exist, choose from 'Linear' or 'Polynomial'"))


  # Simulate the binary treatment assignment for training data
  # Clarification, when trt1, X.e = X.b + X.b%*%B and trt2, X.2 =X.b - X.b%*%B
  Trt.lvls <- c("Trt 1", "Trt 2")
  Trt.train <- factor(Trt.lvls[rbinom(n.train, 1, ratio)+1])

  # Simulate the autoregressive covariance matrix of X
  # cov.mat <- 0.8^(abs(outer(1:p,1:p,"-")))

  # Simulate x.b
  X.train.base <- MASS::mvrnorm(n.train,rep(0,p), cov.mat)


  X.train.end <- X.train.base + .trt.f(X.train.base, ifelse(Trt.train=="Trt 1", 1, -1)) + mvrnorm(n.train, rep(0,p), c.x*diag(p))

  Y.train.base <- .link.f(X.train.base) + mvrnorm(n.train, rep(0,q), c.y*diag(q))
  Y.train.end <- .link.f(X.train.end) + mvrnorm(n.train, rep(0,q), c.y*diag(q))


  ####################################
  # Simulate Testing data
  #
  ####################################

  X.test.base  <- MASS::mvrnorm(n.test,rep(0,p),diag(p))
  # With/Without treatment X.end
  X.test.trt1.end <-  X.test.base + .trt.f(X.test.base, 1)
  X.test.trt2.end <-  X.test.base + .trt.f(X.test.base, -1)


  # With/Without treatment Y.base
  Y.test.trt1.end <- .link.f(X.test.trt1.end)
  Y.test.trt2.end <- .link.f(X.test.trt2.end)

  return(
    list(
      train = list(x.b = X.train.base,
                  x.e = X.train.end,
                  trt = Trt.train,
                  y.b = Y.train.base,
                  y.e = Y.train.end),
      test = list(
        x.b = X.test.base,
        y.e.case = Y.test.trt1.end,
        y.e.control = Y.test.trt2.end
      )
    )
  )
}

#' Title
#' #TODO: improve the documentation
#' @param p The dimension of the design matrix
#'
#' @return A matrix
#' @export
#'
#' @examples
#' create.B(10)
create.B <- function(p){
  if(p < 9)
    stop("Minimum value for p is 9")
  cbind(
    matrix(
      c(1:3, rep(0, p-3),
        rep(0,3), 1:3, rep(0, p-6),
        rep(0,6), 1:3, rep(0, p-9)),
      nrow = p, ncol = 3),
    matrix(0,nrow = p, ncol = p-3)
  )
}

#' Title
#' #TODO: improve the documentation
#' @param p The number of rows
#' @param q The number of columns
#'
#' @return A design matrix
#' @export
#'
#' @examples
#' create.Z(10,3)
create.Z <- function(p, q){
  matrix(
    c(1:3, rep(0, p-3),
      rep(0,3), 1:3, rep(0, p-6),
      rep(0,6), 1:3, rep(0, p-9),
      rep(0, p*(q-3))
      ),
    nrow = p, ncol = q
  )
}
