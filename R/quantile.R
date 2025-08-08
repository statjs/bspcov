#' Quantiles of posterior distribution
#'
#' Compute quantiles to describe posterior distribution.
#' For multiple chains, combines all chains to compute more robust quantiles.
#'
#' @param object an object from \strong{bandPPP}, \strong{bmspcov}, \strong{sbmspcov}, and \strong{thresPPP}.
#' @param probs numeric vector of probabilities with values in [0,1]. Default is c(0.025, 0.5, 0.975).
#' @param ... additional arguments for quantile.
#'
#' @return \item{quantiles}{a list containing quantile matrices for each probability level. For multiple chains, uses combined samples from all chains.}
#' @author Kyeongwon Lee
#' @seealso estimate, plot.postmean.bspcov
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
#' quant <- bspcov::quantile(res)
#' # Get 95% credible intervals
#' quant <- bspcov::quantile(res, probs = c(0.025, 0.975))
#'
quantile <- function(x, ...) {
  UseMethod("quantile")
}

#' @rdname quantile
#' @export
quantile.bspcov <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  stopifnot(!is.null(object$Sigma))
  stopifnot(all(probs >= 0 & probs <= 1))

  # Check if we have multiple chains (nchain > 1)
  if (is.list(object$Sigma)) {
    # Multiple chains case
    # Combine all chains into a single matrix
    combined_samples <- do.call(rbind, object$Sigma)
    
    # Compute quantiles for each column (vech element)
    quantile_matrix <- apply(combined_samples, 2, quantile, probs = probs)
    
  } else {
    # Single chain case (original code)
    post.sample <- object$Sigma
    
    # Compute quantiles for each column (vech element)
    quantile_matrix <- apply(post.sample, 2, quantile, probs = probs)
  }

  # Convert each quantile level to a matrix
  quantiles <- list()
  for (i in seq_along(probs)) {
    if (length(probs) == 1) {
      quantiles[[paste0("q", probs[i])]] <- ks::invvech(quantile_matrix)
    } else {
      quantiles[[paste0("q", probs[i])]] <- ks::invvech(quantile_matrix[i, ])
    }
  }
  
  # Add probability levels as names for easier access
  names(quantiles) <- paste0("q", probs)
  
  class(quantiles) <- 'quantile.bspcov'
  quantiles
}
