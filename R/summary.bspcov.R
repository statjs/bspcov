#' Summary of Posterior Distribution
#'
#' Provides the summary statistics for posterior samples of covariance matrix.
#'
#' @param object an object from \strong{bandPPP}, \strong{bmspcov}, \strong{sbmspcov}, and \strong{thresPPP}.
#' @param cols a scalar or a vector including specific column indices.
#' @param rows a scalar or a vector including specific row indices greater than or equal to columns indices.
#' @param ... additional arguments for the summary function.
#'
#' @note If both \code{cols} and \code{rows} are vectors, they must have the same length.
#'
#' @author Seongil Jo
#'
#' @importFrom coda as.mcmc
#' @importFrom ks invvech
#' @export
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' n <- 100
#' p <- 20
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
#' #summary(fout, cols = 1, rows = 1:p)}
#'
summary.bspcov <- function(object, cols, rows, ...) {
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
  if ((length(cols)>1) & (length(rows)>1)) {
    stopifnot(length(cols) == length(rows))
  }

  # posterior samples
  post.sample <- coda::as.mcmc(t(apply(object$Sigma, 1, ks::invvech)[(cols-1)*p+rows,,drop=FALSE]))
  colnames(post.sample) <- paste('sigma[',rows,',',cols,']',sep='')
  summary(post.sample, ...)
}
