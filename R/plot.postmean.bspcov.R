#' Draw a Heat Map for Point Estimate of Covariance Matrix
#'
#' Provides a heat map for posterior mean estimate of sparse covariance matrix
#'
#' @param x an object from \strong{estimate}.
#' @param ... additional arguments for ggplot2.
#' @return \item{plot}{a heatmap for point estimate of covariance matrix \strong{x}.}
#' @author Seongil Jo
#' @seealso estimate
#'
#' @importFrom ggplot2 ggplot aes geom_tile
#' @export
#'
#' @examples
#'
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::thresPPP(X, eps=0.01, thres=list(value=0.5,fun='hard'), nsample=100)
#' est <- bspcov::estimate(res)
#' plot(est)
#'
plot.postmean.bspcov <- function(x, ...) {
  p <- ncol(x)
  x_grid <- 1:p
  y_grid <- 1:p
  dat <- expand.grid(x_grid = x_grid, y_grid = y_grid)
  z <- as.vector(x)
  dat$z <- z
  ggplot2::ggplot(dat, ggplot2::aes(x = x_grid, y = y_grid, fill = z), ...) +
    ggplot2::geom_tile() + ggplot2::theme_classic() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(title="cov"))
}
