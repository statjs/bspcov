#' Bayesian Estimation of a Sparse Covariance Matrix
#'
#' Provides a post-processed posterior (PPP) for Bayesian inference of a sparse covariance matrix.
#'
#' Lee and Lee (2023) proposed a two-step procedure generating samples from the post-processed posterior for Bayesian inference of a sparse covariance matrix:
#' \itemize{
#'  \item Initial posterior computing step: Generate random samples from the following initial posterior obtained by using the inverse-Wishart prior \eqn{IW_p(B_0, \nu_0)}
#'  \deqn{
#'  \Sigma \mid X_1, \ldots, X_n \sim IW_p(B_0 + nS_n, \nu_0 + n),
#'  }
#'  where \eqn{S_n = n^{-1}\sum_{i=1}^{n}X_iX_i^\top}.
#'  \item Post-processing step: Post-process the samples generated from the initial samples
#'  \deqn{
#'  \Sigma_{(i)} := \left\{\begin{array}{ll}H_{\gamma_n}(\Sigma^{(i)}) + \left[\epsilon_n - \lambda_{\min}\{H_{\gamma_n}(\Sigma^{(i)})\}\right]I_p, &
#'  \mbox{ if } \lambda_{\min}\{H_{\gamma_n}(\Sigma^{(i)})\} < \epsilon_n, \\
#'  H_{\gamma_n}(\Sigma^{(i)}), & \mbox{ otherwise },
#'  \end{array}\right.
#'  }
#' }
#' where \eqn{\Sigma^{(1)}, \ldots, \Sigma^{(N)}} are the initial posterior samples,
#' \eqn{\epsilon_n} is a positive real number, and \eqn{H_{\gamma_n}(\Sigma)} denotes the generalized threshodling operator given as
#' \deqn{
#' (H_{\gamma_n}(\Sigma))_{ij} = \left\{\begin{array}{ll}\sigma_{ij}, & \mbox{ if } i = j, \\
#' h_{\gamma_n}(\sigma_{ij}), & \mbox{ if } i \neq j, \end{array}\right.
#' }
#' where \eqn{\sigma_{ij}} is the \eqn{(i,j)} element of \eqn{\Sigma} and \eqn{h_{\gamma_n}(\cdot)} is a generalized thresholding function.
#'
#' For more details, see Lee and Lee (2023).
#'
#' @param X a n \eqn{\times} p data matrix with column mean zero.
#' @param eps a small positive number decreasing to \eqn{0}.
#' @param thres a list giving the information for thresholding PPP procedure.
#' The list includes the following parameters (with default values in parentheses):
#' \code{value (0.1)} giving the positive real number for the thresholding PPP procedure,
#' \code{fun ('hard')} giving the thresholding function ('hard' or 'soft') for the thresholding PPP procedure.
#' @param prior a list giving the prior information.
#' The list includes the following parameters (with default values in parentheses):
#' \code{A (I)} giving the positive definite scale matrix for the inverse-Wishart prior,
#' \code{nu (p + 1)} giving the degree of freedom of the inverse-Wishar prior.
#' @param nsample a scalar value giving the number of the post-processed posterior samples.
#' @return \item{Sigma}{a nsample \eqn{\times} p(p+1)/2 matrix including lower triangular elements of covariance matrix.}
#' \item{p}{dimension of covariance matrix.}
#' @author Kwangmin Lee
#' @seealso cv.thresPPP
#' @keywords sparse covariance
#'
#' @references Lee, K. and Lee, J. (2023), "Post-processes posteriors for sparse covariances", \emph{Journal of Econometrics}.
#'
#' @importFrom purrr map
#' @importFrom CholWishart rInvWishart
#' @importFrom matrixcalc vech
#' @importFrom FinCovRegularization hard.thresholding soft.thresholding
#' @import magrittr
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
#'
thresPPP <- function(X, eps, thres = list(), prior = list(), nsample = 2000){
  p <- ncol(X)
  n <- nrow(X)

  # thresholding parameters
  thresvals <- list(value = 0.1, fun = 'hard')
  thresvals[names(thres)] <- thres
  thrval <- thresvals$value
  thrfun <- thresvals$fun
  if (thrfun == 'hard') {
    thresfun <- FinCovRegularization::hard.thresholding
  } else {
    thresfun <- FinCovRegularization::soft.thresholding
  }

  # prior
  privals <- list(nu = p + 1, A = diag(1,p))
  privals[names(prior)] <- prior
  A = privals$A
  nu = privals$nu

  # initial posterior step
  parameter_initposterior <- list(A = crossprod(X) + A, nu = nu + n)

  # post-processing step
  Sigma_save <- purrr::map(1:nsample,
                           ~(CholWishart::rInvWishart(1, parameter_initposterior$nu, parameter_initposterior$A)[,,1]) %>%
                             thresfun(thres = thrval) %>%
                             pd_adjustment_Matrix(eps = eps)) %>%
    purrr::map(~matrixcalc::vech(as.matrix(.x)))

  out <- list()
  out$ppp <- 'sparse'
  out$Sigma <- matrix(unlist(Sigma_save), nrow = nsample, ncol = p*(p+1)/2, byrow = TRUE)
  out$p <- p
  class(out) <- 'bspcov'
  out
}

