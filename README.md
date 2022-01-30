# bspcov
An R package for Bayesian Sparse Estimation of a Covariance Matrix

# Reference
Lee, Jo, and Lee (2021+). The beta-mixture shrinkage prior for sparse covariances with posteriornear-minimax rate.

# An example in the manuscript
### Dimension
p = 12
n = 250

### True covariance structure
TrSig <- matrix(0, nrow = p, ncol = p) 

diag(TrSig) <- c(0.239,1.554,0.362,0.199,0.349,0.295,0.715,0.164,0.518,0.379,0.159,0.207)  
TrSig[c(2,8),1] <- TrSig[1,c(2,8)] <- c(0.117,0.031)  
TrSig[4,3] <- TrSig[3,4] <- 0.002  
TrSig[5,4] <- TrSig[4,5] <- 0.094  
TrSig[p,5] <- TrSig[5,p] <- -0.036  
TrSig[c(7,8),6] <- TrSig[6,c(7,8)] <- c(-0.229,0.002)  
TrSig[9:11,8] <- TrSig[8,9:11] <- c(0.112,-0.028,-0.008)  
TrSig[10:11,9] <- TrSig[9,10:11] <- c(-0.193,-0.090)  
TrSig[11,10] <- TrSig[10,11] <- 0.167  
iTrSig <- solve(TrSig)  

### Generate datasets
set.seed(7)  
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = TrSig)

### mcmc parameters
burnin <- 5000  
nmc <- 5000

### hyper-parameters for the beta-mixture shrinkage prior
lambda <- 1  
a <- 1/2  
b <- 1/2  

### run the algorithm
foo <- bspcov::bmspcov(S=crossprod(X), n=n, Sig=cov(X), a=a, b=b, lambda=lambda, burnin=burnin, nmc=nmc)

est.sig <- apply(foo$Sig, c(1,2), mean)  
Matrix::norm(TrSig - est.sig, type = 'F')/p

### variable selection performance
bspcov::postmcmc_vs(foo)
