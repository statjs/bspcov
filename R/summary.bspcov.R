#' Summary of Posterior Distribution
#'
#' Provides the summary statistics for posterior samples of covariance matrix.
#'
#' @param object an object from \strong{bandPPP}, \strong{bmspcov}, \strong{sbmspcov}, and \strong{thresPPP}.
#' @param cols a scalar or a vector including specific column indices.
#' @param rows a scalar or a vector including specific row indices greater than or equal to columns indices.
#' @param quantiles a numeric vector of quantiles to compute. Default is c(0.025, 0.25, 0.5, 0.75, 0.975).
#' @param ... additional arguments for the summary function.
#' @return \item{summary}{a table of summary statistics including empirical mean, standard deviation, and quantiles for posterior samples. For multiple chains, also includes effective sample size (n_eff) and R-hat diagnostics.}
#' @note If both \code{cols} and \code{rows} are vectors, they must have the same length.
#'
#' @author Seongil Jo and Kyeongwon Lee
#'
#' @importFrom coda as.mcmc mcmc.list effectiveSize gelman.diag
#' @importFrom ks invvech
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
#' summary(fout, cols = c(1, 3, 4), rows = c(1, 3, 4))
#' summary(fout, cols = 1, rows = 1:p)
#'
#' # Custom quantiles:
#' summary(fout, cols = 1, rows = 1:3, quantiles = c(0.05, 0.5, 0.95))
#'
#' # Multiple chains with diagnostics:
#' # fout_multi <- bspcov::sbmspcov(X = X, Sigma = diag(diag(cov(X))), nchain = 3)
#' # summary(fout_multi, cols = 1, rows = 1:3)  # Shows n_eff, R-hat, and nchain
#'
summary.bspcov <- function(object, cols, rows, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  stopifnot(!is.null(object$Sigma))
  p <- object$p
  if (missing(cols)) {
    cols <- 1
  }
  if (missing(rows)) {
    rows <- 1:p
  }
  stopifnot(all(cols <= p))
  stopifnot(all(rows <= p))
  if ((length(cols)>1) && (length(rows)>1)) {
    stopifnot(length(cols) == length(rows))
  }

  # Handle multiple chains vs single chain
  if (is.list(object$Sigma)) {
    # Multiple chains case - create mcmc.list for proper summary and diagnostics
    nchain <- length(object$Sigma)
    mcmc_list <- list()
    
    cat("Number of chains:", nchain, "\n\n")
    
    for (i in 1:nchain) {
      chain_sample <- t(apply(object$Sigma[[i]], 1, ks::invvech)[(cols-1)*p+rows,,drop=FALSE])
      colnames(chain_sample) <- paste('sigma[',rows,',',cols,']',sep='')
      mcmc_list[[i]] <- coda::as.mcmc(chain_sample)
    }
    
    # Convert to mcmc.list
    mcmc_list <- coda::mcmc.list(mcmc_list)
    
    # Get basic summary
    basic_summary <- summary(mcmc_list, quantiles = quantiles, ...)
    
    # Add MCMC diagnostics
    cat("MCMC Summary Statistics:\n")
    print(basic_summary)
    
    cat("\nMCMC Diagnostics:\n")
    
    # Effective sample size
    n_eff <- coda::effectiveSize(mcmc_list)
    cat("\nEffective Sample Size (n_eff):\n")
    print(round(n_eff, 1))
    
    # R-hat (Gelman-Rubin diagnostic)
    if (nchain > 1) {
      tryCatch({
        rhat <- coda::gelman.diag(mcmc_list, multivariate = FALSE)
        cat("\nPotential Scale Reduction Factor (R-hat):\n")
        print(round(rhat$psrf, 3))
        
        if (any(rhat$psrf[,1] > 1.1)) {
          cat("\nWarning: Some R-hat values > 1.1, indicating potential convergence issues.\n")
        }
      }, error = function(e) {
        cat("\nNote: R-hat calculation failed. This may occur with very short chains or identical starting values.\n")
      })
    }
    
    # Return invisibly for programmatic access
    invisible(list(
      summary = basic_summary,
      n_eff = n_eff,
      rhat = if(nchain > 1) tryCatch(coda::gelman.diag(mcmc_list, multivariate = FALSE), error = function(e) NULL) else NULL,
      nchain = nchain
    ))
    
  } else {
    # Single chain case
    post.sample <- coda::as.mcmc(t(apply(object$Sigma, 1, ks::invvech)[(cols-1)*p+rows,,drop=FALSE]))
    colnames(post.sample) <- paste('sigma[',rows,',',cols,']',sep='')
    
    # Get summary with custom quantiles
    basic_summary <- summary(post.sample, quantiles = quantiles, ...)
    print(basic_summary)
    
    # Return invisibly for programmatic access
    invisible(list(
      summary = basic_summary,
      nchain = 1
    ))
  }
}
