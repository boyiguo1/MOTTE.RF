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
#' sim_MOTTE_data_cross_over( n.train = 500, n.test = 200,
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
sim_MOTTE_data_cross_over <- function(
  n.train = 500, # <-  500
  n.test = 200, # <- 200
  p = 10, #  <- 10
  q = 3, #  <- 3
  ratio = 0.5, # <- 0.5,
  # TODO: set up error prevention for cov.mat, i.e. measuring the dimension
  cov.mat = diag(p),
  # TODO: incorporate the linear and polynormial for both treat.f and link.f in the code
  trt.f = c("Linear", "Polynomial", "Box"),
  link.f = c("Linear", "Polynomial"),
  B ,
  Z
){

  # TODO: extract this as argument for the parameter
  # c.x <- sum((1:3)^2)/2
  # c.y <- sum((1:3)^2)/2

  c.x <- 0.5
  c.y <- 0.5
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
                   "Box" = function(x, trt){
                     .x <- x
                     for (i in 1: nrow(.x)) {
                       if(abs(x[i,1])<1 & abs(x[i,2])<1) .x[i, 1:3] <- 0
                       if(abs(x[i,4])<1 & abs(x[i,5])<1) .x[i, 4:6] <- 0
                       if(abs(x[i,7])<1 & abs(x[i,8])<1) .x[i, 7:9] <- 0
                     }
                     sweep(.x, 1, trt, "*") %*% B
                   },
                   stop("Trt.f doesn't exist, choose from 'Linear' or 'Polynomial' or 'Box'")
  )


  # Simulate the binary treatment assignment for training data
  # Clarification, when trt1, X.e = X.b + X.b%*%B and trt2, X.2 =X.b - X.b%*%B
  Trt.lvls <- c("Trt 1", "Trt 2")
  Trt.train <- factor(Trt.lvls[rbinom(n.train, 1, ratio)+1])

  # Simulate x.b
  X.train.base.1 <- MASS::mvrnorm(n.train,rep(0,p), cov.mat)
  X.train.base.2 <- X.train.base.1 + MASS::mvrnorm(n.train,rep(0,p), 0.1*diag(p))
  X.train.end.1 <- X.train.base.1 + .trt.f(X.train.base.1, 1)+ mvrnorm(n.train, rep(0,p), c.x*diag(p))
  X.train.end.2 <- X.train.base.2 + .trt.f(X.train.base.2, -1)+ mvrnorm(n.train, rep(0,p), c.x*diag(p))
  # X.train.end <- X.train.base + .trt.f(X.train.base, ifelse(Trt.train=="Trt 1", 1, -1)) + mvrnorm(n.train, rep(0,p), c.x*diag(p))

  Y.train.base.1 <- .link.f(X.train.base.1) + mvrnorm(n.train, rep(0,q), c.y*diag(q))
  Y.train.base.2 <- .link.f(X.train.base.2) + mvrnorm(n.train, rep(0,q), c.y*diag(q))
  Y.train.end.1 <- .link.f(X.train.end.1) + mvrnorm(n.train, rep(0,q), c.y*diag(q))
  Y.train.end.2 <- .link.f(X.train.end.2) + mvrnorm(n.train, rep(0,q), c.y*diag(q))

  ####################################
  # Simulate Testing data
  #
  ####################################

  #X.test.base  <- mvrnorm(n.test,rep(0,p),x.sig)
  X.test.base  <- MASS::mvrnorm(n.test,rep(0,p), cov.mat)
  # X.test.base.2  <- X.test.base.1 + MASS::mvrnorm(n.test,rep(0,p), 0.1*diag(p))
  # With/Without treatment X.end
  X.test.trt1.end <-  X.test.base + .trt.f(X.test.base, 1)
  X.test.trt2.end <-  X.test.base + .trt.f(X.test.base, -1)


  # With/Without treatment Y.base
  Y.test.trt1.end <- .link.f(X.test.trt1.end)
  Y.test.trt2.end <- .link.f(X.test.trt2.end)

  # TODO: add return list
  return(
    list(
      train = list(x.b.1 = X.train.base.1,
                   x.b.2 = X.train.base.2,
                   x.e.1 = X.train.end.1,
                   x.e.2 = X.train.end.2,
                   #trt = Trt.train,
                   y.b.1 = Y.train.base.1,
                   y.b.2 = Y.train.base.2,
                   y.e.1 = Y.train.end.1,
                   y.e.2 = Y.train.end.2),
      test = list(
        x.b = X.test.base,
        # x.b.2 = X.test.base.2,
        y.e.1 = Y.test.trt1.end,
        y.e.2 = Y.test.trt2.end
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
    matrix(
      c(1:3, rep(0, p-3),
        rep(0,3), 1:3, rep(0, p-6),
        rep(0,6), 1:3, rep(0, p-9)),
      nrow = p, ncol = 3),
    matrix(
      c(1:3, rep(0, p-3),
        rep(0,3), 1:3, rep(0, p-6),
        rep(0,6), 1:3, rep(0, p-9)),
      nrow = p, ncol = 3),
    matrix(0,nrow = p, ncol = p-9)
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
    c((1:3)/3, rep(0, p-3),
      rep(0,3), (1:3)/3, rep(0, p-6),
      rep(0,6), (1:3)/3, rep(0, p-9),
      rep(0, p*(q-3))
    ),
    nrow = p, ncol = q
  )
}
