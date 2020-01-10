source("global.R")

stable.CCA <- function (x, y, xcenter = TRUE, ycenter = TRUE, epsilon = 1e-04) 
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  if ((nr <- nrow(x)) != nrow(y)) 
    stop("unequal number of rows in 'cancor'")
  ncx <- ncol(x)
  ncy <- ncol(y)
  if (!nr || !ncx || !ncy) 
    stop("dimension 0 in 'x' or 'y'")
  if (is.logical(xcenter)) {
    if (xcenter) {
      xcenter <- colMeans(x, )
      x <- x - rep(xcenter, rep.int(nr, ncx))
    }
    else xcenter <- rep.int(0, ncx)
  }
  else {
    xcenter <- rep_len(xcenter, ncx)
    x <- x - rep(xcenter, rep.int(nr, ncx))
  }
  if (is.logical(ycenter)) {
    if (ycenter) {
      ycenter <- colMeans(y)
      y <- y - rep(ycenter, rep.int(nr, ncy))
    }
    else ycenter <- rep.int(0, ncy)
  }
  else {
    ycenter <- rep_len(ycenter, ncy)
    y <- y - rep(ycenter, rep.int(nr, ncy))
  }
  
  ## Centering finished
  
  qr.x <- qr(x)
  qr.y <- qr(y)
  
  p.x <- qr.x$pivot
  p.y <- qr.y$pivot
  
  ## QR decomposition finished
  R.x <- qr.R(qr.x, complete=TRUE)
  R.y <- qr.R(qr.y, complete=TRUE)
  
  Xi.x <- findMaxIndex(R.x, epsilon)
  Xi.y <- findMaxIndex(R.y, epsilon)

  Q.x <- qr.Q(qr.x, complete=TRUE)[,1:Xi.x,drop=FALSE]
  Q.y <- qr.Q(qr.y, complete=TRUE)[,1:Xi.y,drop=FALSE]

  R.x <- R.x[1:Xi.x, 1:Xi.x,drop=FALSE]
  R.y <- R.y[1:Xi.y, 1:Xi.y,drop=FALSE]
  
  v <- min(Xi.x, Xi.y)
  
  if(Xi.x >= Xi.y){
    svd.res <- svd(crossprod(Q.x,Q.y))
    u <- svd.res$u
    d <- svd.res$d
    z <- svd.res$v
  }
  else {
    svd.res <- svd(crossprod(Q.y,Q.x))
    u <- svd.res$z
    d <- svd.res$d
    z <- svd.res$u
  }
  
  u <- u[,1:v,drop=FALSE]
  z <- z[,1:v,drop=FALSE]
  
  A <- backsolve(R.x, u)
  B <- backsolve(R.y, z)
  
  xcoef <- matrix(0,ncol=ncol(A),nrow=ncol(x))
  ycoef <- matrix(0,ncol=ncol(B),nrow=ncol(y))
  
  xcoef[p.x[1:nrow(A)],] <- A
  ycoef[p.y[1:nrow(B)],] <- B
  
  # dx <- qx$rank
  # if (!dx) 
  #   stop("'x' has rank 0")
  # dy <- qy$rank
  # if (!dy) 
  #   stop("'y' has rank 0")
  # z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , 
  #                                                 drop = FALSE], dx, dy)
  # xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
  # rownames(xcoef) <- colnames(x)[qx$pivot][1L:dx]
  # ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)
  # rownames(ycoef) <- colnames(y)[qy$pivot][1L:dy]
  # list(cor = z$d, xcoef = xcoef, ycoef = ycoef, xcenter = xcenter, 
  #      ycenter = ycenter)
  
  return(
    list(cor = d, xcoef = xcoef, ycoef=ycoef,
         xcenter=xcenter, ycenter=ycenter)
  )
}