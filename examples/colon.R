# attach library
library(bspcov)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)

# fix random seed
set.seed(123)

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

# visulization

## function to visualize posterior mean covariance matrix
vis_postmean <- function(postmean, title) {
  p <- ncol(postmean)
  sigmah_mat <- matrix(as.numeric(postmean), nrow = p, ncol = p)
  df <- melt(sigmah_mat, varnames = c("x", "y"), value.name = "cov")
  ggplot(df, aes(x, y, fill = cov)) +
    geom_raster(interpolate = FALSE) +
    labs(title = title, x = "Gene", y = "Gene") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
}

## compute range for color scale
min_sigmah <- min(c(as.numeric(sigmah), as.numeric(sigmah.two)))
max_sigmah <- max(c(as.numeric(sigmah), as.numeric(sigmah.two)))
cont_scale <- scale_fill_gradient(
  name   = "",
  low    = "black", 
  high   = "white",
  limits = c(min_sigmah, max_sigmah),
  guide  = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(0.5, "cm"))
)

## plots
p1 <- vis_postmean(sigmah, "bspcov")
p2 <- vis_postmean(sigmah.two, "sbmspcov")
g <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  cont_scale &
  theme(legend.position = "bottom")
ggsave("figs/colon.png", g, width = 12, height = 7)

# LDA with the posterior mean estimate (sigmah)
mu.1 <- as.matrix(colMeans(X[normal_idx, ]))
mu.2 <- as.matrix(colMeans(X[tumor_idx, ]))
omega.1 <- N.n / (N.n + N.t)
omega.2 <- N.t / (N.n + N.t)
lda <- function(X, mu.1, mu.2, sigmah, omega.1, omega.2) {
  delat.1 <- as.numeric(X %*% solve(sigmah) %*% mu.1) - as.numeric(0.5 * t(mu.1) %*% solve(sigmah) %*% mu.1 + log(omega.1))
  delat.2 <- as.numeric(X %*% solve(sigmah) %*% mu.2) - as.numeric(0.5 * t(mu.2) %*% solve(sigmah) %*% mu.2 + log(omega.2))
  delta <- ifelse(delat.1 > delat.2, 1, 2)
  return(delta)
}
# predict using the posterior mean estimate
pred <- lda(X, mu.1, mu.2, sigmah, omega.1, omega.2)
pred.two <- lda(X, mu.1, mu.2, sigmah.two, omega.1, omega.2)
# misclassification rate
misclass <- mean(pred != group)
misclass.two <- mean(pred.two != group)
cat("Misclassification rate for bspcov:", misclass, "\n")
cat("Misclassification rate for sbmspcov:", misclass.two, "\n")
