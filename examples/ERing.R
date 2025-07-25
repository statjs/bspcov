# attach library
library(bspcov)
library(dplyr)

# fix random seed
set.seed(1234)

# load and preprocess data
library(mlmts)
data("ERing")
Ering_class1 <- ERing$data[which(ERing$classes == 1)] %>%
  purrr::map(~ as.vector(t(.x[, 1]))) %>%
  do.call("rbind", .)
Ering_class2 <- ERing$data[which(ERing$classes == 2)] %>%
  purrr::map(~ as.vector(t(.x[, 1]))) %>%
  do.call("rbind", .)


# estimation by band PPP
## gesture 1
cov1 <- bandPPP(scale(Ering_class1, scale = FALSE, center = TRUE), k = 3, eps = 0.05)
postmean1 <- estimate(cov1)
## gesture 2
cov2 <- bandPPP(scale(Ering_class2, scale = FALSE, center = TRUE), k = 3, eps = 0.05)
postmean2 <- estimate(cov2)

# ----- Analysis -----
library(purrr)
library(readr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(caret)
library(mvnfast)
library(patchwork)

# visualization

## function to visualize posterior mean covariance matrix
vis_postmean <- function(postmean, title) {
  p <- ncol(postmean)
  sigmah <- matrix(as.numeric(postmean), nrow = p, ncol = p)
  df <- melt(sigmah, varnames = c("x", "y"), value.name = "cov")
  ggplot(df, aes(x, y, fill = cov)) +
    geom_raster(interpolate = FALSE) +
    labs(title = title, x = "v1", y = "v2") +
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
min_postmean <- min(c(as.numeric(postmean1), as.numeric(postmean2)))
max_postmean <- max(c(as.numeric(postmean1), as.numeric(postmean2)))
cont_scale <- scale_fill_gradient(
  name   = "",
  low    = "black",
  high   = "white",
  limits = c(min_postmean, max_postmean),
  guide  = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(0.5, "cm"))
)

## plots
p1 <- vis_postmean(postmean1, "Gesture 1")
p2 <- vis_postmean(postmean2, "Gesture 2")
g <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  cont_scale &
  theme(legend.position = "bottom")
ggsave("figs/cov_gestures.png", g, width = 12, height = 7)

## Cross-validation for QDA
k_folds <- 5

# Create folds
folds1 <- caret::createFolds(seq_len(nrow(Ering_class1)), k = k_folds, returnTrain = TRUE)
folds2 <- caret::createFolds(seq_len(nrow(Ering_class2)), k = k_folds, returnTrain = TRUE)
folds <- list(folds1, folds2)

## testacc function, QDA implementation and evaluation
testacc <- function(ind1, ind2) {
  train1 <- Ering_class1[ind1, ]
  train2 <- Ering_class2[ind2, ]
  test1 <- Ering_class1[-ind1, ]
  test2 <- Ering_class2[-ind2, ]

  cov1 <- bandPPP(scale(train1, scale = FALSE, center = TRUE), k = 3, eps = 0.05)
  postmean1 <- estimate(cov1)
  cov2 <- bandPPP(scale(train2, scale = FALSE, center = TRUE), k = 3, eps = 0.05)
  postmean2 <- estimate(cov2)
  prob1 <- mean(mvnfast::dmvn(test1, colMeans(train1), postmean1, log = T) > mvnfast::dmvn(test1, colMeans(train2), postmean2, log = T))
  prob2 <- mean(mvnfast::dmvn(test2, colMeans(train1), postmean1, log = T) < mvnfast::dmvn(test2, colMeans(train2), postmean2, log = T))

  c(prob1, prob2)
}

res <- purrr::map2(folds[[1]], folds[[2]], testacc)
do.call("rbind", res) %>% colMeans()
# [1] 1.0000000 0.7420202