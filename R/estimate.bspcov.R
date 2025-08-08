#' Point-estimate of posterior distribution
#'
#' Compute the point estimate (mean) to describe posterior distribution.
#' For multiple chains, combines all chains to compute a more robust estimate.
#'
#' @param object an object from \strong{bandPPP}, \strong{bmspcov}, \strong{sbmspcov}, and \strong{thresPPP}.
#' @param ... additional arguments for estimate.
#'
#' @return \item{Sigma}{the point estimate (mean) of covariance matrix. For multiple chains, uses combined samples from all chains.}
#' @author Seongil Jo and Kyeongwon Lee
#' @seealso plot.postmean.bspcov
#'
#' @importFrom ks invvech
#' @export
#'
#' @examples
#'
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::bandPPP(X,2,0.01,nsample=100)
#' est <- bspcov::estimate(res)
#'
estimate <- function(object, ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
estimate.bspcov <- function(object, ...) {
  stopifnot(!is.null(object$Sigma))

  # Check if we have multiple chains (nchain > 1)
  if (is.list(object$Sigma)) {
    # Multiple chains case
    # Combine all chains into a single matrix
    combined_samples <- do.call(rbind, object$Sigma)
    
    # posterior mean
    post.est.m <- ks::invvech(colMeans(combined_samples))
    
  } else {
    # Single chain case (original code)
    post.sample <- object$Sigma
    
    # posterior mean
    post.est.m <- ks::invvech(colMeans(post.sample))
  }

  class(post.est.m) <- 'postmean.bspcov'
  post.est.m
}
