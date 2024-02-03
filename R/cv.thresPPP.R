#' CV for Bayesian Estimation of a Sparse Covariance Matrix
#'
#' Performs cross-validation to estimate spectral norm error for a post-processed posterior of a sparse covariance matrix.
#'
#' Given a set of train data and validation data, the spectral norm error for each \eqn{\gamma} and \eqn{\epsilon_n} is estimated as follows:
#'  \deqn{
#'  ||\hat{\Sigma}(\gamma,\epsilon_n)^{(train)} - S^{(val)}||_2
#'  }
#' where \eqn{\hat{\Sigma}(\gamma,\epsilon_n)^{(train)}} is the estimate for the covariance based on the train data and \eqn{S^{(val)}} is the sample covariance matrix based on the validation data.
#' The spectral norm error is estimated by the \eqn{10}-fold cross-validation.
#' For more details, see the first paragraph on page 9 in Lee and Lee (2023).
#'
#' @param X a n \eqn{\times} p data matrix with column mean zero.
#' @param thresvec a vector of real numbers specifying the parameter of the threshold function.
#' @param epsvec a vector of small positive numbers decreasing to \eqn{0}.
#' @param prior a list giving the prior information.
#' The list includes the following parameters (with default values in parentheses):
#' \code{A (I)} giving the positive definite scale matrix for the inverse-Wishart prior,
#' \code{nu (p + k)} giving the degree of freedom of the inverse-Wishar prior.
#' @param thresfun a string to specify the type of threshold function. \code{fun ('hard')} giving the thresholding function ('hard' or 'soft') for the thresholding PPP procedure.
#' @param nsample a scalar value giving the number of the post-processed posterior samples.
#' @param ncores a scalar value giving the number of CPU cores.
#' @return \item{CVdf}{a M \eqn{\times} 3 dataframe having the estimated spectral norm error for each thres and eps, where M = length(thresvec) * length(epsvec)}
#' @author Kwangmin Lee
#' @seealso thresPPP
#' @keywords sparse covariance
#'
#' @references Lee, K. and Lee, J. (2023), "Post-processes posteriors for sparse covariances", \emph{Journal of Econometrics}, 236(3), 105475.
#'
#' @importFrom purrr map reduce
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom FinCovRegularization hard.thresholding soft.thresholding
#' @importFrom ks invvech
#' @importFrom caret createMultiFolds
#' @importFrom dplyr group_by summarise
#' @importFrom plyr alply
#' @import magrittr
#' @export
#'
#' @examples
#'
#' \dontrun{
#' Sigma0 <- diag(1,50)
#' X <- mvtnorm::rmvnorm(25,sigma = Sigma0)
#' thresvec <- c(0.01,0.1)
#' epsvec <- c(0.01,0.1)
#' res <- bspcov::cv.thresPPP(X,thresvec,epsvec,nsample=100)
#' plot(res)}
#'
cv.thresPPP <- function(X, thresvec, epsvec, prior = NULL,
                        thresfun = 'hard',
                        nsample = 2000, ncores = 2){

  thresPPPCV <- function(Xtrain, Xnew, tuneparams, prior=NULL,
                         thresfun = 'hard', num.mcmc=100){
    n <- nrow(Xtrain)
    p <- ncol(Xtrain)
    if(is.null(prior)){
      prior <- list(A=diag(mean(diag(cov(Xtrain))),p),nu=p+1)
    }
    parameter_initposterior <- list(A=crossprod(Xtrain)+ prior$A,
                                    nu=prior$nu + n)

    Snew <- crossprod(Xnew)/nrow(Xnew)


    PPPmean_list <- 1:nrow(tuneparams) %>%
      purrr::map(~(thresPPP(Xtrain,thres= list(value = tuneparams$thres[.x], fun = thresfun),
                            eps=tuneparams$eps[.x],prior=prior,nsample)$Sigma %>%
                     plyr::alply(1,ks::invvech) %>%
                     purrr::reduce(`+`))/nsample)

    PPPmean_list %>% purrr::map(~norm(.x-Snew,type="2")) %>%
      unlist %>% data.frame(tuneparams,err=.)

  }

  future::plan(future::multisession, workers = ncores)

  tuneparams <- expand.grid(thres=thresvec,eps=epsvec)
  list.trainind <- caret::createMultiFolds(1:nrow(X),times = 1)

  CVdf <- list.trainind %>%
    furrr::future_map(~thresPPPCV(X[.x,],X[-.x,],tuneparams,prior,
                                  thresfun=thresfun),
                      .options = furrr::furrr_options(seed=TRUE)) %>%
    do.call("rbind",.) %>% dplyr::group_by(thres,eps) %>%
    dplyr::summarise(err=mean(err)) %>% dplyr::arrange(err)

  out <- list()
  out$ppp <- 'cv.sparse'
  out$error <- CVdf
  class(out) <- 'bspcov'
  out

}
