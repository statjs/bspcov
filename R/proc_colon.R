#' Preprocess Colon Gene Expression Data
#'
#' @description
#' The \code{proc_colon} function preprocesses colon gene expression data by:
#' \enumerate{
#'   \item Log transforming the raw counts.
#'   \item Performing two-sample t-tests for each gene between normal and tumor samples.
#'   \item Selecting the top 50 genes by absolute t-statistic.
#'   \item Returning the filtered expression matrix and sample indices/groups.
#' }
#'
#' @param colon
#'   A numeric matrix of raw colon gene expression values (genes × samples).
#'   Rows are genes; columns are samples.
#' @param tissues
#'   A numeric vector indicating tissue type per sample:
#'   positive for normal, negative for tumor.
#'
#' @return A list with components:
#'   \describe{
#'     \item{X}{A numeric matrix (samples x 50 genes) of selected, log‐transformed expression values.}
#'     \item{normal_idx}{Integer indices of normal‐tissue columns in the original data.}
#'     \item{tumor_idx}{Integer indices of tumor‐tissue columns in the original data.}
#'     \item{group}{Integer vector of length \code{ncol(colon)}, with 1 = normal, 2 = tumor.}
#'   }
#'
#' @importFrom stats t.test
#' @export
#'
#' @examples
#' data("colon")
#' data("tissues")
#' set.seed(1234)
#' colon_data <- proc_colon(colon, tissues)
#' X <- colon_data$X
#' \donttest{
#' foo <- bmspcov(X, Sigma = cov(X))
#' sigmah <- estimate(foo)
#' }

proc_colon <- function(colon, tissues) {
  colon <- log10(colon)
  N <- ncol(colon)
  p <- nrow(colon)

  # preprocessing the data
  normal_idx <- which(tissues > 0)
  N.n <- length(normal_idx)
  tumor_idx <- which(tissues < 0)
  N.t <- length(tumor_idx)
  group <- ifelse(tissues > 0, 1, 2)

  normal <- colon[, normal_idx]
  tumor <- colon[, tumor_idx]

  # select genes by two-sample t-test
  stats <- pvalues <- numeric(p)
  for (i in 1:p) {
    out <- t.test(normal[i, ], tumor[i, ])
    stats[i] <- abs(out$statistic)
    pvalues[i] <- out$p.value
  }
  o <- order(stats, decreasing = TRUE)
  colon.o <- colon[o, ]

  q <- 50
  X <- t(colon.o[1:q, ])

  return(list(
    X = X,
    normal_idx = normal_idx,
    tumor_idx = tumor_idx,
    group = group
  ))
}
