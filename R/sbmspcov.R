#' Bayesian Sparse Covariance Estimation using Sure Screening
#'
#' Provides a Bayesian sparse and positive definite estimate of a covariance matrix via screened beta-mixture prior.
#'
#' Lee, Jo, Lee, and Lee (2023+) proposed the screened beta-mixture shrinkage prior for estimating a sparse and positive definite covariance matrix.
#' The screened beta-mixture shrinkage prior for \eqn{\Sigma = (\sigma_{jk})} is defined as
#' \deqn{
#'  \pi(\Sigma) = \frac{\pi^{u}(\Sigma)I(\Sigma \in C_p)}{\pi^{u}(\Sigma \in C_p)}, ~ C_p = \{\mbox{ all } p \times p \mbox{ positive definite matrices }\},
#' }
#' where \eqn{\pi^{u}(\cdot)} is the unconstrained prior given by
#' \deqn{
#' \pi^{u}(\sigma_{jk} \mid \psi_{jk}) = N\left(\sigma_{jk} \mid 0, \frac{\psi_{jk}}{1 - \psi_{jk}}\tau_1^2\right), ~ \psi_{jk} = 1 - 1/(1 + \phi_{jk})}
#' \deqn{
#' \pi^{u}(\psi_{jk}) = Beta(\psi_{jk} \mid a, b) ~ \mbox{for } (j, k) \in S_r(\hat{R})}
#' \deqn{
#' \pi^{u}(\sigma_{jj}) = Exp(\sigma_{jj} \mid \lambda),
#' }
#' where \eqn{S_r(\hat{R}) = \{(j,k) : 1 < j < k \le p, |\hat{\rho}_{jk}| > r\}}, \eqn{\hat{\rho}_{jk}} are sample correlations, and \eqn{r} is a threshold given by user.
#'
#' For more details, see Lee, Jo, Lee and Lee (2022+).
#'
#' @param X a n \eqn{\times} p data matrix with column mean zero.
#' @param Sigma an initial guess for Sigma.
#' @param cutoff a list giving the information for the threshold.
#' The list includes the following parameters (with default values in parentheses):
#' \code{method ('FNR')} giving the method for determining the threshold value
#' (method='FNR' uses the false negative rate (FNR)-based approach, method='corr' chooses the threshold value by sample correlations),
#' \code{rho} a lower bound of meaningfully large correlations whose the defaults values are 0.25 and 0.2 for method = 'FNR' and method = 'corr', respectively.
#' Note. If method='corr', \code{rho} is used as the threshold value.
#' \code{FNR (0.05)} giving the prespecified FNR level when method = 'FNR'.
#' \code{nsimdata (1000)} giving the number of simulated datasets for calculating Jeffreys’ default Bayes factors when method = 'FNR'.
#' @param prior a list giving the prior information.
#' The list includes the following parameters (with default values in parentheses):
#' \code{a (1/2)} and \code{b (1/2)} giving the shape parameters for beta distribution,
#' \code{lambda (1)} giving the hyperparameter for the diagonal elements,
#' \code{tau1sq (log(p)/(p^2*n))} giving the hyperparameter for the shrinkage prior of covariance.
#' @param nsample a list giving the MCMC parameters.
#' The list includes the following integers (with default values in parentheses):
#' \code{burnin (1000)} giving the number of MCMC samples in transition period,
#' \code{nmc (1000)} giving the number of MCMC samples for analysis.
#'
#' @return \item{Sigma}{a nmc \eqn{\times} p(p+1)/2 matrix including lower triangular elements of covariance matrix.}
#' \item{p}{dimension of covariance matrix.}
#' \item{Phi}{a nmc \eqn{\times} p(p+1)/2 matrix including shrinkage parameters corresponding to lower triangular elements of covariance matrix.}
#' \item{INDzero}{a list including indices of off-diagonal elements screened by sure screening.}
#' \item{cutoff}{the cutoff value specified by FNR-approach.}
#' @author Kyoungjae Lee and Seongil Jo
#' @seealso bmspcov
#' @keywords sparse covariance
#'
#' @references Lee, K., Jo, S., Lee, K., and Lee, J. (2023+), "Scalable and optimal Bayesian inference for sparse covariance matrices via screened beta-mixture prior",
#' 	arXiv:2206.12773.
#'
#' @importFrom stats rnorm rgamma
#' @importFrom mvnfast rmvn
#' @importFrom GIGrvg rgig
#' @importFrom progress progress_bar
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' n <- 20
#' p <- 5
#'
#' # generate a sparse covariance matrix:
#' True.Sigma <- matrix(0, nrow = p, ncol = p)
#' diag(True.Sigma) <- 1
#' Values <- -runif(n = p*(p-1)/2, min = 0.2, max = 0.8)
#' nonzeroIND <- which(rbinom(n=p*(p-1)/2,1,prob=1/p)==1)
#' zeroIND = (1:(p*(p-1)/2))[-nonzeroIND]
#' Values[zeroIND] <- 0
#' True.Sigma[lower.tri(True.Sigma)] <- Values
#' True.Sigma[upper.tri(True.Sigma)] <- t(True.Sigma)[upper.tri(True.Sigma)]
#' if(min(eigen(True.Sigma)$values) <= 0){
#'   delta <- -min(eigen(True.Sigma)$values) + 1.0e-5
#'   True.Sigma <- True.Sigma + delta*diag(p)
#' }
#'
#' # generate a data
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = True.Sigma)
#'
#' # compute sparse, positive covariance estimator:
#' fout <- bspcov::sbmspcov(X = X, Sigma = diag(diag(cov(X))))
#' post.est.m <- bspcov::estimate(fout)
#' sqrt(mean((post.est.m - True.Sigma)^2))
#' sqrt(mean((cov(X) - True.Sigma)^2))
#'
sbmspcov <- function(X, Sigma, cutoff = list(), prior = list(), nsample = list()) {
  # Estimate a sparse covariance matrix using the screened beta-mixture shrinkage prior

  n <- dim(X)[1]
  p <- dim(X)[2]
  stopifnot(p > 1)

  # Select a cutoff
  if (is.null(cutoff$method)) {
    cutoff$method <- 'FNR'
  }
  if (cutoff$method == 'FNR') {
    cutoff_vals <- list(rho=0.25, FNR=0.05, nsimdata=1000)
    cutoff_vals[names(cutoff)] <- cutoff
    C.th <- select_cutoff(n, rho=cutoff_vals$rho, FNR=cutoff_vals$FNR, nsimdata=cutoff_vals$nsimdata)

    # Sure Screening
    INDzero <- BayesCGM.SS(X, p, C.th)
    for (i in 1:p) {
      if (length(INDzero[[i]]) > 0) {
        Sigma[INDzero[[i]],i] <- 0
      }
    }
  } else if (cutoff$method == 'corr') {
    cutoff_vals <- list(thr = 0.2)
    cutoff_vals[names(cutoff)] <- cutoff

    # Sure Screening
    INDzero <- BayesCGM.SS.CORR(X, p, cutoff_vals$thr)
    for (i in 1:p) {
      if (length(INDzero[[i]]) > 0) {
        Sigma[INDzero[[i]],i] <- 0
      }
    }
  }

  # prior information
  privals <- list(a=1/2, b=1/2, lambda=1, tau1sq=log(p)/(p^2*n))
  privals[names(prior)] <- prior

  # mcmc parameters
  mcvals <- list(burnin=1000, nmc=1000)
  mcvals[names(nsample)] <- nsample

  # Estimating covariance using SBM prior
  sbmout <- sbm.covest(X=X, S=crossprod(X), n=n, p=p, Sigma=Sigma, a=privals$a, b=privals$b,
                       lambda=privals$lambda, tau1sq=privals$tau1sq,
                       INDzero, burnin=mcvals$burnin, nmc=mcvals$nmc)

  out <- sbmout
  out$INDzero <- INDzero
  if (cutoff$method == 'FNR') {
    out$cutoff <- C.th
  }
  class(out) <- 'bspcov'
  out
}




# Estimating covariance using SBM prior
sbm.covest <- function(X, S, n, p, Sigma, a, b, lambda, tau1sq, INDzero, burnin, nmc) {

  ind_noi_all <- list()
  ind_nozeroi_all <- list()
  ind_noi_add_all <- list()
  for (i in 1:p) {
    ind_noi_all[[i]] <- (1:p)[-i]
    if (length(INDzero[[i]])>0) {
      ind_nozeroi_all[[i]] <- (1:p)[-c(i,INDzero[[i]])]
      ind_tmp <- c()
      for (j in 1:length(INDzero[[i]])) {
        ind_tmp <- c(ind_tmp, which(INDzero[[i]][j] == ind_noi_all[[i]]))
      }
      ind_noi_add_all[[i]] <- (1:(p-1))[-ind_tmp]
    }else {
      ind_nozeroi_all[[i]] <- ind_noi_all[[i]]
      ind_noi_add_all[[i]] <- 1:(p-1)
    }
  }

  # an array for posterior samples
  Sigma_save <- Phi_save <- matrix(0, nrow = nmc, ncol = p*(p+1)/2)

  # initial values
  min_eig <- min(eigen(Sigma, only.values = TRUE)$values)
  if(min_eig <= 1e-15) {
    Sigma <- Sigma + (-min_eig + 0.001)*diag(p)
  }
  C <- solve(Sigma)
  tau <- matrix(tau1sq, p, p)
  Phi <- Psi <- matrix(1,p,p)

  # fixed parameters
  lam <- 1-n/2

  # blocked gibbs sampler
  pb <- progress::progress_bar$new(total = burnin+nmc)
  elapsed.time <- system.time({
  for (iter in 1:(burnin+nmc)) {
    pb$tick()

    for (i in 1:p) {
      ind_noi <- ind_noi_all[[i]]
      ind_nozeroi <- ind_nozeroi_all[[i]]
      ind_noi_add <- ind_noi_add_all[[i]]

      if(length(ind_nozeroi_all[[i]]) != 0) {
        tau_temp <- tau[ind_nozeroi, i]

        C11 <- C[ind_noi, ind_noi]; C12 <- C[ind_noi, i]
        S11 <- S[ind_noi, ind_noi]; S12 <- S[ind_noi, i]

        invSig11 <- C11 - tcrossprod(C12)/C[i,i]
        invSig11.reduced = invSig11[ind_noi_add, , drop=FALSE]
        invSig11S12.reduced <- invSig11.reduced %*% S12

        if(n >= p){
          W1 = tcrossprod(invSig11.reduced %*% S11, invSig11.reduced)
        }else{
          X.invSig11.reduced = tcrossprod(X[,ind_noi, drop=FALSE], invSig11.reduced)
          W1 = crossprod(X.invSig11.reduced)
        }

        # Sample gamma (= v in the paper)
        beta <- Sigma[ind_nozeroi, i, drop = FALSE]
        chi <- drop(crossprod(beta, W1 %*% beta) - 2*crossprod(beta, invSig11S12.reduced) + S[i,i])
        psi <- lambda
        gam <- GIGrvg::rgig(n=1, lambda = lam, chi = chi, psi = psi)
        if (gam <= 1e-06) {
          gam <- 1e-06
        }

        # Sample beta (= u in the paper)
        W <- W1/gam + diag(1/tau_temp,length(ind_noi_add)) + lambda*invSig11.reduced[ ,ind_noi_add, drop=FALSE]
        W <- (W + t(W))/2
        tryCatch(
          {
            W_chol <- chol(W)
          },
          error = function(cond) {
            W <- W + (-min(eigen(W)$values) + 0.001)*diag(ncol(W))
            W_chol <<- chol(W)
          }
        )
        mu_i <- backsolve(W_chol, forwardsolve(t(W_chol), invSig11S12.reduced)) / gam
        beta <- mu_i + backsolve(W_chol, rnorm(length(mu_i)))

        Sigma[ind_nozeroi, i] <- beta
        Sigma[i, ind_nozeroi] <- beta
        Sigma[i,i] <- gam + crossprod(beta, invSig11[ind_noi_add,ind_noi_add] %*% beta)

        for(j in ind_nozeroi) {
          chi.phi <- Sigma[j,i]^2/tau1sq
          if (chi.phi <= 1e-06) {
            chi.phi <- 1e-06
          }
          Phi[j,i] <- GIGrvg::rgig(n=1, lambda = a-1/2, chi = chi.phi, psi = 2*Psi[j,i])
          Phi[i,j] <- Phi[j,i]

          Psi[j,i] <- stats::rgamma(n=1, a+b, Phi[j,i]+1)
          Psi[i,j] <- Psi[j,i]
        }
        tau[ind_nozeroi,i] <- Phi[ind_nozeroi,i] * tau1sq
        tau[i,ind_nozeroi] <- tau[ind_nozeroi,i]

        # Below updating Precision matrix according to one-column change of precision matrix
        invSig11beta <- crossprod(invSig11.reduced, beta)
        C[ind_noi, ind_noi] <- invSig11 + tcrossprod(invSig11beta, invSig11beta)/gam
        C12 <- -invSig11beta/gam
        C[ind_noi, i] <- C12
        C[i, ind_noi] <- t(C12)
        C[i, i] <- 1/gam
      } else {
        gam <- GIGrvg::rgig(n=1, lambda = lam, chi = S[i,i], psi = lambda)
        if (gam <= 1e-06) {
          gam <- 1e-06
        }

        Sigma[i,i] <- gam
        C[i, i] <- 1/gam
      }
    }

    if (iter > burnin) {
      Sigma_save[iter-burnin,] <- Sigma[!upper.tri(Sigma)]
      Phi_save[iter-burnin,] <- Phi[!upper.tri(Sigma)]
    }
  }})[3]

  out <- list()
  out$prior <- 'sbmsp'
  out$p <- p
  out$Sigma <- Sigma_save
  out$Phi <- Phi_save
  out$mcmctime <- elapsed.time
  out
}
