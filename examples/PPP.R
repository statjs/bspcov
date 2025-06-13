# fix random seed
set.seed(123)

n <- 25
p <- 50
Sigma0 <- diag(1, p)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)

# post-processed posterior (PPP) for Bayesian inference of a banded covariance matrix
res_bandPPP <- bspcov::bandPPP(X, 2, 0.01, nsample = 100)

# post-processed posterior (PPP) for Bayesian inference of a sparse covariance matrix
thresvec <- c(0.01, 0.1)
epsvec <- c(0.01, 0.1)
res_thresPPP <- bspcov::thresPPP(X, eps = 0.01, thres = list(value = 0.5, fun = "hard"), nsample = 100)

# results
est.sig_band <- bspcov::estimate(res_bandPPP)
est.sig_thres <- bspcov::estimate(res_thresPPP)
Matrix::norm(Sigma0 - est.sig_band, type = "F") / p
Matrix::norm(Sigma0 - est.sig_thres, type = "F") / p
