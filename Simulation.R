library(MASS)     # R library supporting multivariate normal simulation

####################################
# Simulate Setting
#
####################################

n.train <-  500
n.test <- 200
p <- 10
q <- 3
# pi is the ratio between treatment groups
pi <- 0.5

# Treatment effect function
treat.effect <- function(x.base)
{
  ifelse(x.base[,2]<0.25 & x.base[,5]<0.25,sqrt(x.base[,2]^2+x.base[,5]^2),0)
}

# Beta matrix
Z <- matrix(
  cbind(
    c(rep(1,5),rep(0,5)),
    c(rep(0,2),rep(1,5),rep(0,3)),
    c(rep(0,5),rep(1,5))
  ),nrow=p,ncol=q)


####################################
# Simulate Training data
#
####################################
# Check if the portion size is integer

# Simulate the treatment
Treat.train <- sample(
  c(
    rep(1, floor(n.train*pi)),           # Treatment Group 1
    rep(0, ceiling(n.train*(1-pi)))      # Treatment Group 2
  )
)

# Simulate the autoregressive covariance matrix of X
x.sig <- 0.8^(abs(outer(1:p,1:p,"-")))

# Simulate x.b
X.train.base <- mvrnorm(n.train,rep(0,p),x.sig)


X.train.end <- X.train.base
X.train.end[,5] <- X.train.end[,5] + Treat.train*treat.effect(X.train.base)
X.train.end[,6] <- X.train.end[,6] + Treat.train*treat.effect(X.train.base)
X.train.end[,7]<- X.train.end[,7] + Treat.train*treat.effect(X.train.base)
X.train.end <- X.train.end + mvrnorm(n.train, rep(0,p), 0.1*diag(p))   

Y.train.base <- (X.train.base)%*%Z + mvrnorm(n.train, rep(0,q), 0.1*diag(q))
Y.train.end <- (X.train.end)%*%Z + mvrnorm(n.train, rep(0,q), 0.1*diag(q))


####################################
# Simulate Testing data
#
####################################

X.test.base  <- mvrnorm(n.test,rep(0,p),x.sig)
# With/Without treatment X.end
X.test.case.end <-  X.test.base
X.test.control.end <- X.test.base
X.test.case.end[,5] <- X.test.case.end[,5] + rep(1,n.test)*treat.effect(X.test.base)
X.test.case.end[,6] <- X.test.case.end[,6] + rep(1,n.test)*treat.effect(X.test.base)
X.test.case.end[,7] <- X.test.case.end[,7] + rep(1,n.test)*treat.effect(X.test.base)
X.test.control.end[,5] <- X.test.control.end[,5] + rep(0,n.test)*treat.effect(X.test.base)
X.test.control.end[,6] <- X.test.control.end[,6] + rep(0,n.test)*treat.effect(X.test.base)
X.test.control.end[,7] <- X.test.control.end[,7] + rep(0,n.test)*treat.effect(X.test.base)
X.test.case.end <- X.test.case.end + mvrnorm(n.test,rep(0,p),0.1*diag(p)) 
X.test.control.end <- X.test.control.end + mvrnorm(n.test,rep(0,p),0.1*diag(p))
# With/Without treatment Y.end
Y.test.case.end <- (X.test.case.end)%*%Z+mvrnorm(n.test,rep(0,q),0.1*diag(q))
Y.test.control.end <- (X.test.control.end)%*%Z+mvrnorm(n.test,rep(0,q),0.1*diag(q))

# Generating Rdata file for GUIDE
# train.lab variable in this file is for the prediction purpose for Loh's method
train.lab <- c(rep(1,n.train),rep(0, 2*n.test))
Xs <- rbind(X.train.base, X.test.base, X.test.base)
Ys <- rbind(Y.train.end, matrix(NA,nrow=2*n.test,ncol=ncol(Y.train.end)))
Treat.all <- c(Treat.train,rep(1,n.test),rep(0,n.test))
all.data <- cbind(train.lab, Treat.all, Xs, Ys)
colnames(all.data) <- c("Train","Treated",paste0("X",1:p),paste0("Y",1:q))
# Save data for Loh's method
write.csv(all.data,paste0(args[1],"input_data.rdata"),row.names=F)



