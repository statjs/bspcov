# attach library
library(bspcov)
library(gridExtra)
# load data
data("colon")
colon <- log10(colon)
N <- ncol(colon)
p <- nrow(colon)
data("tissues")

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

# beta-mixture prior
foo <- bmspcov(X, Sigma = cov(X)) # X: n x p matrix

# an option for choosing a sure-screening threshold using false negative rate (FNR)
cutoff <- list(method = "FNR", rho = 0.25, FNR = 0.1)
# screened beta-mixture prior
foo.two <- sbmspcov(X = X, Sigma = cov(X), cutoff)

# posterior mean estimate
sigmah <- estimate(foo)
sigmah.two <- estimate(foo.two)

# heatmap in Figure 2
g <- grid.arrange(plot(sigmah), plot(sigmah.two), ncol = 2)
# ggsave("figs/colon.png", plot = g, width = 12, height = 6)
