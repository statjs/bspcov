#' CV for Bayesian Estimation of a Banded Covariance Matrix
#'
#' Performs leave-one-out cross-validation (LOOCV) to calculate the predictive log-likelihood for a post-processed posterior of a banded covariance matrix and selects the optimal parameters.
#'
#' The predictive log-likelihood for each \eqn{k} and \eqn{\epsilon_n} is estimated as follows:
#'  \deqn{
#'  \sum_{i=1}^n \log S^{-1} \sum_{s=1}^S p(X_i \mid B_k^{(\epsilon_n)}(\Sigma_{i,s})),
#'  }
#' where \eqn{X_i} is the ith observation, \eqn{\Sigma_{i,s}} is the sth posterior sample based on \eqn{(X_1,\ldots,X_{i-1},X_{i+1},\ldots,X_n)}, and \eqn{B_k^{(\epsilon_n)}} represents the banding post-processing function.
#' For more details, see (3) in Lee, Lee and Lee (2023+).
#'
#' @param X a n \eqn{\times} p data matrix with column mean zero.
#' @param kvec a vector of natural numbers specifying the bandwidth of covariance matrix.
#' @param epsvec a vector of small positive numbers decreasing to \eqn{0}.
#' @param prior a list giving the prior information.
#' The list includes the following parameters (with default values in parentheses):
#' \code{A (I)} giving the positive definite scale matrix for the inverse-Wishart prior,
#' \code{nu (p + k)} giving the degree of freedom of the inverse-Wishar prior.
#' @param nsample a scalar value giving the number of the post-processed posterior samples.
#' @param ncores a scalar value giving the number of CPU cores.
#' @return \item{elpd}{a M \eqn{\times} 3 dataframe having the expected log predictive density (ELPD) for each k and eps, where M = length(kvec) * length(epsvec).}
#' @author Kwangmin Lee
#' @seealso bandPPP
#' @keywords banded covariance
#'
#' @references Lee, K., Lee, K., and Lee, J. (2023+), "Post-processes posteriors for banded covariances",
#' 	\emph{Bayesian Analysis}, DOI: 10.1214/22-BA1333.
#'
#' 	Gelman, A., Hwang, J., and Vehtari, A. (2014). "Understanding predictive information criteria for Bayesian models." \emph{Statistics and computing},
#' 	24(6), 997-1016.
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession
#' @importFrom dplyr arrange desc
#' @import magrittr
#' @export
#'
#' @examples
#'
#' \donttest{
#' Sigma0 <- diag(1,50)
#' X <- mvtnorm::rmvnorm(25,sigma = Sigma0)
#' kvec <- 1:2
#' epsvec <- c(0.01,0.05)
#' res <- bspcov::cv.bandPPP(X,kvec,epsvec,nsample=10,ncores=4)
#' plot(res)}
#'
cv.bandPPP <- function(X, kvec, epsvec, prior = list(), nsample = 2000, ncores = 2){
  stopifnot(!missing(kvec))
  stopifnot(!missing(epsvec))
  future::plan(future::multisession, workers = ncores)

  tuneparams <- expand.grid(k = kvec, eps = epsvec)
  elpd <- furrr::future_map(1:dim(tuneparams)[1],
                            ~mean(loo_pred_loglik(X = X, k = tuneparams$k[.x],
                                                  eps = tuneparams$eps[.x],
                                                  nsample = nsample, prior = prior)),
                            .options = furrr::furrr_options(seed = TRUE)) %>% unlist %>%
    cbind(tuneparams, logpdf=.) %>% dplyr::arrange(dplyr::desc(logpdf))

  out <- list()
  out$ppp <- 'cv.band'
  out$elpd <- elpd
  class(out) <- 'bspcov'
  out
}
