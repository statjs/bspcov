#'
#'
#' @examples
#'
#' data("colon")
#' data("tissues")
#' set.seed(1234)
#' colon_data <- proc_colon(colon, tissues)
#' X <- colon_data$X
#' \donttest{
#' foo <- bmspcov(X, Sigma = cov(X))
#' sigmah <- estimate(foo)}
#'
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
