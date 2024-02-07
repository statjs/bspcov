#' Point-estimate of posterior distribution
#'
#' Compute the point estimate (mean) to describe posterior distribution.
#'
#' @param object an object from \strong{bandPPP}, \strong{bmspcov}, \strong{sbmspcov}, and \strong{thresPPP}.
#' @param ... additional arguments for estimate.
#'
#' @return \item{Sigma}{the point estimate (mean) of covariance matrix.}
#' @author Seongil Jo
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

  # posterior samples
  post.sample <- object$Sigma
  p <- object$p
  nsample <- nrow(post.sample)

  # posterior mean
  post.est.m <- ks::invvech(colMeans(post.sample))

  class(post.est.m) <- 'postmean.bspcov'
  post.est.m
}
