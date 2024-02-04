#' Plot Diagnostics of Posterior Samples and Cross-Validation
#'
#' Provides a trace plot of posterior samples and a plot of a learning curve for cross-validation
#'
#' @param x an object from \strong{bmspcov}, \strong{sbmspcov}, \strong{cv.bandPPP}, and \strong{cv.thresPPP}.
#' @param ... additional arguments for ggplot2.
#' @param cols a scalar or a vector including specific column indices for the trace plot.
#' @param rows a scalar or a vector including specific row indices greater than or equal to columns indices for the trace plot.
#' @return \item{plot}{a plot for diagnostics of posterior samples \strong{x}.}
#' @author Seongil Jo
#'
#' @importFrom coda mcmc
#' @importFrom ggmcmc ggs ggs_traceplot
#' @importFrom magrittr `%>%`
#' @importFrom mvtnorm rmvnorm
#' @export
#'
#' @examples
#'
#' \donttest{
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
#' plot(fout, cols = c(1, 3, 4), rows = c(1, 3, 4))
#' plot(fout, cols = 1, rows = 1:3)
#'
#' # Cross-Validation for Banded Structure
#' Sigma0 <- diag(1,50)
#' X <- mvtnorm::rmvnorm(25,sigma = Sigma0)
#' kvec <- 1:2
#' epsvec <- c(0.01,0.05)
#' res <- bspcov::cv.bandPPP(X,kvec,epsvec,nsample=10,ncores=4)
#' plot(res)}
#'
plot.bspcov <- function(x, ..., cols, rows) {
  if (is.null(x$ppp)) {
    p <- x$p
    if (missing(cols)) {
      cols <- 1
    }
    if (missing(rows)) {
      rows <- 1:3
    }
    stopifnot(all(cols <= p))
    stopifnot(all(rows <= p))
    if ((length(cols)>1) & (length(rows)>1)) {
      stopifnot(length(cols) == length(rows))
    }

    # posterior samples
    post.sample <- t(apply(x$Sigma, 1, ks::invvech)[(cols-1)*p+rows,,drop=FALSE])
    colnames(post.sample) <- paste('sigma[',rows,',',cols,']',sep='')
    post.sample %>% mcmc %>% ggs %>% ggs_traceplot
  } else {
    if (x$ppp == 'cv.band') {
      ggplot2::qplot(y = x$elpd$logpdf, ylab = 'Expected Log Predictive Density (ELPD)', ...) +
        ggplot2::geom_line(size = 1) + ggplot2::theme_bw()
    } else {
      ggplot2::qplot(y = x$error$err, ylab = 'Spectral Norm Error', ...) +
        ggplot2::geom_line(size = 1) + ggplot2::theme_bw()
    }
  }
}
