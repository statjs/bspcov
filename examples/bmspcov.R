# Dimension
p <- 12
n <- 250

# True covariance
TrSig <- matrix(0, nrow = p, ncol = p)
diag(TrSig) <- c(0.239, 1.554, 0.362, 0.199, 0.349, 0.295, 0.715, 0.164, 0.518, 0.379, 0.159, 0.207)
TrSig[c(2, 8), 1] <- TrSig[1, c(2, 8)] <- c(0.117, 0.031)
TrSig[4, 3] <- TrSig[3, 4] <- 0.002
TrSig[5, 4] <- TrSig[4, 5] <- 0.094
TrSig[p, 5] <- TrSig[5, p] <- -0.036
TrSig[c(7, 8), 6] <- TrSig[6, c(7, 8)] <- c(-0.229, 0.002)
TrSig[9:11, 8] <- TrSig[8, 9:11] <- c(0.112, -0.028, -0.008)
TrSig[10:11, 9] <- TrSig[9, 10:11] <- c(-0.193, -0.090)
TrSig[11, 10] <- TrSig[10, 11] <- 0.167
iTrSig <- solve(TrSig)

# Generate datasets
set.seed(7)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = TrSig)

# Bayesian sparse and positive definite estimate of a covariance matrix via the beta-mixture shrinkage prior
foo_bm <- bspcov::bmspcov(X, Sigma = cov(X))

# Bayesian sparse and positive definite estimate of a covariance matrix via screened beta-mixture prior
cutoff <- list(method = "FNR", rho = 0.25, FNR = 0.05)
foo_sbm <- bspcov::sbmspcov(X = X, Sigma = cov(X), cutoff = cutoff)

# results
est.sig_bm <- bspcov::estimate(foo_bm)
est.sig_sbm <- bspcov::estimate(foo_sbm)
Matrix::norm(TrSig - est.sig_bm, type = "F") / p
Matrix::norm(TrSig - est.sig_sbm, type = "F") / p
