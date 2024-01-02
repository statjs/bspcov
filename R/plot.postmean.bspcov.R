#' Draw a Heat Map for Point Estimate of Covariance Matrix
#'
#' Provides a heat map for posterior mean estimate of sparse covariance matrix
#'
#' @param object an object from \strong{estimate}.
#'
#' @author Seongil Jo
#' @seealso estimate
#'
#' @importFrom ggplot2 ggplot aes geom_tile
#' @export
#'
#' @examples
#'
#' \dontrun{
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::thresPPP(X, eps=0.01, thres=list(value=0.5,fun='hard'), nsample=100)
#' est <- bspcov::estimate(res)
#' plot(est)}
#'
plot.postmean.bspcov <- function(object, ...) {
  p <- ncol(object)
  x <- 1:p
  y <- 1:p
  dat <- expand.grid(x=x, y=y)
  dat$z <- as.vector(object)
  ggplot2::ggplot(dat, ggplot2::aes(x, y, fill= z), ...) +
    ggplot2::geom_tile() + ggplot2::theme_classic() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(title="cov"))
}
