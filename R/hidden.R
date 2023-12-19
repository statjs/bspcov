#' @importFrom MASS mvrnorm
#' @importFrom stats cor quantile
#' @importFrom BayesFactor correlationBF
#' @importFrom purrr map
#' @importFrom matrixStats logSumExp
#' @importFrom mvnfast dmvn
#' @importFrom RSpectra eigs
#' @importFrom Matrix Diagonal
#' @importFrom plyr alply
#' @import dplyr
#'

# Calculate the log-predictive likelihoods given posterior samples (postsample) and test data (Xnew).
pred_loglik <- function(postsample, Xnew){
  p <- dim(postsample[[1]])[1]
  logpdf <- postsample %>%
    purrr::map(~mvnfast::dmvn(Xnew, mu=rep(0,p), .x, log = TRUE)) %>%
    unlist
  matrixStats::logSumExp(logpdf) - log(length(postsample))
}



# Wrapper function for the LOOCV to calculate the log-predictive likelihoods (pred_loglik).
# This function splits data X into the train data and test observation and calculate the log-predictive likelihood.
loo_pred_loglik <- function(X, k, eps=NULL, nsample=100, prior=list()){
  n <- dim(X)[1]

  1:n %>% purrr::map(~pred_loglik(
    bandPPP(X[-.x,],k=k,eps=eps,nsample=nsample,prior=prior)$Sigma %>%
      plyr::alply(1,ks::invvech),X[.x,])) %>% unlist

}



# This function adjusts Sigma to be positive-definite by adding a diagonal matrix. See the case "lambda_{min}(B_k(Sigma)) <eps"
# in equation (1) of Lee, Lee and Lee (2023+).
pd_adjustment_Matrix <- function(Sigma, eps=1e-04){
  egmin <- RSpectra::eigs(Sigma, 1, sigma = -RSpectra::eigs(Sigma, 1, which = "LM")$values)$values
  if(egmin < eps){
    return(Sigma + Matrix::Diagonal(dim(Sigma)[1],-egmin+eps))
  }
  return (Sigma)
}





# Select a cutoff
select_cutoff <- function(n, rho, FNR, nsimdata=1000, seed) {
  # N: Null
  # A: Alternative

  p <- 2
  tSig.A <- matrix(c(1,rho,rho,1),p,p)

  BF.A <- numeric(nsimdata)
  if (!missing(seed)) {
    set.seed(seed=seed)
  } else {
    set.seed(sample(.Machine$integer.max,1))
  }
  for (isim in 1:nsimdata) {
    X.A = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = tSig.A)
    BF.A[isim] <- exp(BayesFactor::correlationBF(X.A[,1], X.A[,2], rscale = "ultrawide")@bayesFactor$bf)
  }

  # selected cutoff satisfying FNR criterion
  cutoff = quantile(BF.A, probs = FNR)

  return(cutoff)
}



# Sure Screening
pairwise.Jeffreys <- function(X){
  n = nrow(X)
  p = ncol(X)
  BF10.mat = matrix(NA, p, p) # whose (i,j)th entry is BF_{10}(i,j) for i>j (lower triangular matrix)

  for(i in 2:p){
    for(j in 1:(i-1)){
      BF10.mat[i,j] = exp(BayesFactor::correlationBF(X[,i], X[,j], rscale = "ultrawide")@bayesFactor$bf)
    }
  }

  return(BF10.mat)
}


# Find non-zero indices via sure screening
BayesCGM.SS <- function(X, p, C.th) {
  foo.BF <- pairwise.Jeffreys(X)
  support.mat <- (foo.BF > C.th)  # should use this object
  support.mat[upper.tri(support.mat)] <- t(support.mat)[upper.tri(support.mat)]

  INDzero <- list()
  for (i in 1:p) {
    INDzero[[i]] <- which(!support.mat[,i])
  }
  INDzero
}


# Find non-zero indices via sample correlations
BayesCGM.SS.CORR <- function(X, p, thr) {
  foo.cor <- stats::cor(X)
  foo.cor[upper.tri(foo.cor,diag=TRUE)] <- NA
  support.mat <- (abs(foo.cor) > stats::quantile(abs(foo.cor), prob = 1-thr, na.rm=TRUE))  # should use this object
  support.mat[upper.tri(support.mat)] <- t(support.mat)[upper.tri(support.mat)]

  INDzero <- list()
  for (i in 1:p) {
    INDzero[[i]] <- which(!support.mat[,i])
  }
  INDzero
}

