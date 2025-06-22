# attach library
library(bspcov)

# fix random seed
set.seed(1234)

# load and preprocess data
data("colon")
data("tissues")
colon_data <- proc_colon(colon, tissues = tissues)
X <- colon_data$X

# beta-mixture prior
foo <- bmspcov(X, Sigma = cov(X))

# screened beta-mixture prior
cutoff <- list(method = "corr", rho = 0.25)
foo.two <- sbmspcov(X = X, Sigma = cov(X), cutoff)

# posterior mean estimate
sigmah <- estimate(foo)
sigmah.two <- estimate(foo.two)

# ------- Analysis -------
library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)
library(dplyr)
library(caret)
library(mvnfast)

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
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
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

## Cross-validation for LDA
k_folds <- 5

# Create folds
X_normal <- X[colon_data$normal_idx, ]
X_tumor <- X[colon_data$tumor_idx, ]
folds1 <- caret::createFolds(seq_len(nrow(X_normal)), k = k_folds, returnTrain = TRUE)
folds2 <- caret::createFolds(seq_len(nrow(X_tumor)), k = k_folds, returnTrain = TRUE)
folds <- list(folds1, folds2)

## testacc function, LDA implementation and evaluation
testacc <- function(ind1, ind2) {
  train1 <- X_normal[ind1, ]
  train2 <- X_tumor[ind2, ]
  test1 <- X_normal[-ind1, ]
  test2 <- X_tumor[-ind2, ]

  train <- rbind(train1, train2)
  foo <- bmspcov(train, Sigma = cov(train))
  foo.two <- sbmspcov(train, Sigma = cov(train), cutoff = cutoff)
  sigmah <- estimate(foo)
  sigmah.two <- estimate(foo.two)

  count1 <- sum(mvnfast::dmvn(test1, colMeans(train1), sigmah, log = TRUE) >
                mvnfast::dmvn(test1, colMeans(train2), sigmah, log = TRUE))
  count2 <- sum(mvnfast::dmvn(test2, colMeans(train1), sigmah, log = TRUE) <
                mvnfast::dmvn(test2, colMeans(train2), sigmah, log = TRUE))                      
  count <- count1 + count2
  prob <- count / (nrow(test1) + nrow(test2))

  count1.two <- sum(mvnfast::dmvn(test1, colMeans(train1), sigmah.two, log = TRUE) >
                   mvnfast::dmvn(test1, colMeans(train2), sigmah.two, log = TRUE))
  count2.two <- sum(mvnfast::dmvn(test2, colMeans(train1), sigmah.two, log = TRUE) <
                   mvnfast::dmvn(test2, colMeans(train2), sigmah.two, log = TRUE))
  count.two <- count1.two + count2.two
  prob.two <- count.two / (nrow(test1) + nrow(test2))

  c(prob, prob.two)
}
res <- purrr::map2(folds[[1]], folds[[2]], ~testacc(.x, .y))
acc <- do.call("rbind", res) %>% colMeans() %>% round(4)
print(1 - acc)
# [1] 0.1146 0.1802