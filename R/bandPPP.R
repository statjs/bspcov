#' Bayesian Estimation of a Banded Covariance Matrix
#'
#' Provides a post-processed posterior for Bayesian inference of a banded covariance matrix.
#'
#' Lee, Lee, and Lee (2023+) proposed a two-step procedure generating samples from the post-processed posterior for Bayesian inference of a banded covariance matrix:
#' \itemize{
#'  \item Initial posterior computing step: Generate random samples from the following initial posterior obtained by using the inverse-Wishart prior \eqn{IW_p(B_0, \nu_0)}
#'  \deqn{
#'  \Sigma \mid X_1, \ldots, X_n \sim IW_p(B_0 + nS_n, \nu_0 + n),
#'  }
#'  where \eqn{S_n = n^{-1}\sum_{i=1}^{n}X_iX_i^\top}.
#'  \item Post-processing step: Post-process the samples generated from the initial samples
#'  \deqn{
#'  \Sigma_{(i)} := \left\{\begin{array}{ll}B_{k}(\Sigma^{(i)}) + \left[\epsilon_n - \lambda_{\min}\{B_{k}(\Sigma^{(i)})\}\right]I_p, &
#'  \mbox{ if } \lambda_{\min}\{B_{k}(\Sigma^{(i)})\} < \epsilon_n, \\
#'  B_{k}(\Sigma^{(i)}), & \mbox{ otherwise },
#'  \end{array}\right.
#'  }
#' }
#' where \eqn{\Sigma^{(1)}, \ldots, \Sigma^{(N)}} are the initial posterior samples,
#' \eqn{\epsilon_n} is a small positive number decreasing to \eqn{0} as \eqn{n \rightarrow \infty},
#' and \eqn{B_k(B)} denotes the \eqn{k}-band operation given as
#' \deqn{
#' B_{k}(B) = \{b_{ij}I(|i - j| \le k)\} \mbox{ for any } B = (b_{ij}) \in R^{p \times p}.
#' }
#' For more details, see Lee, Lee and Lee (2023+).
#'
#' @param X a n \eqn{\times} p data matrix with column mean zero.
#' @param k a scalar value (natural number) specifying the bandwidth of covariance matrix.
#' @param eps a small positive number decreasing to \eqn{0} with default value \eqn{(log(k))^2 * (k + log(p))/n}.
#' @param prior a list giving the prior information.
#' The list includes the following parameters (with default values in parentheses):
#' \code{A (I)} giving the positive definite scale matrix for the inverse-Wishart prior,
#' \code{nu (p + k)} giving the degree of freedom of the inverse-Wishar prior.
#' @param nsample a scalar value giving the number of the post-processed posterior samples.
#' @return \item{Sigma}{a nsample \eqn{\times} p(p+1)/2 matrix including lower triangular elements of covariance matrix.}
#' \item{p}{dimension of covariance matrix.}
#' @author Kwangmin Lee
#' @seealso cv.bandPPP estimate
#' @keywords banded covariance
#'
#' @references Lee, K., Lee, K., and Lee, J. (2023+), "Post-processes posteriors for banded covariances",
#' 	\emph{Bayesian Analysis}, DOI: 10.1214/22-BA1333.
#'
#' @importFrom purrr rerun
#' @importFrom CholWishart rInvWishart
#' @importFrom Matrix Matrix band
#' @importFrom matrixcalc vech
#' @import magrittr
#' @export
#'
#' @examples
#'
#' n <- 25
#' p <- 50
#' Sigma0 <- diag(1, p)
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
#' res <- bspcov::bandPPP(X,2,0.01,nsample=100)
#'
bandPPP <- function(X, k, eps, prior = list(), nsample = 2000){
  p <- dim(X)[2]
  n <- dim(X)[1]
  stopifnot(!missing(k))

  # prior
  privals <- list(nu = p + k, A = diag(1,p))
  privals[names(prior)] <- prior
  A = privals$A
  nu = privals$nu

  if(missing(eps)){
    eps <- (log(k))^2 * (k + log(p))/n
  }

  # initial posterior step
  param_initpost <- list(A = crossprod(X) + A, nu = nu + n)

  # post-processing step
  Sigma_save <- purrr::rerun(nsample, CholWishart::rInvWishart(1, param_initpost$nu, param_initpost$A)[,,1] %>%
                               Matrix::Matrix(sparse = T) %>%
                               Matrix::band(k1 = -k, k2 = k) %>% pd_adjustment_Matrix(eps = eps)) %>%
    purrr::map(~matrixcalc::vech(as.matrix(.x)))

  out <- list()
  out$ppp <- 'band'
  out$Sigma <- matrix(unlist(Sigma_save), nrow = nsample, ncol = p*(p+1)/2, byrow = TRUE)
  out$p <- p
  class(out) <- 'bspcov'
  out
}




