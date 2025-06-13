# attach library
library(bspcov)
library(mlmts)
library(purrr)
library(readr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(caret)
library(dplyr)

set.seed(123)

# load and preprocess data
data("ERing")
Ering_class3 <- ERing$data[which(ERing$classes == 3)] %>%
  purrr::map(~ as.vector(t(.x[, 1]))) %>%
  do.call("rbind", .) %>%
  scale()
Ering_class2 <- ERing$data[which(ERing$classes == 2)] %>%
  purrr::map(~ as.vector(t(.x[, 1]))) %>%
  do.call("rbind", .) %>%
  scale()

# estimation by band PPP
p <- ncol(Ering_class2)
cov1 <- bandPPP(Ering_class3, k = 10, eps = 0.1)
postmean1 <- estimate(cov1)
cov2 <- bandPPP(Ering_class2, k = 10, eps = 0.1)
postmean2 <- estimate(cov2)

# visulization

## function to visualize posterior mean covariance matrix
vis_postmean <- function(postmean, title) {
  p <- ncol(postmean)
  sigmah <- matrix(as.numeric(postmean), nrow = p, ncol = p)
  df <- melt(sigmah, varnames = c("x", "y"), value.name = "cov")
  ggplot(df, aes(x, y, fill = cov)) +
    geom_raster(interpolate = FALSE) +
    labs(title = title, x = "v1", y = "v2") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())
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


# QDA implementation and evaluation

## Calculate means for each class
mu1 <- colMeans(Ering_class3)
mu2 <- colMeans(Ering_class2)

## Function to calculate log-likelihood for multivariate normal
log_likelihood <- function(x, mu, Sigma) {
  p <- length(mu)
  dev <- x - mu
  -0.5 * (p * log(2 * pi) + log(det(Sigma)) + t(dev) %*% solve(Sigma) %*% dev)
}

## QDA classifier function
qda_classify <- function(x, mu1, Sigma1, mu2, Sigma2, pi1 = 0.5, pi2 = 0.5) {
  delta1 <- log_likelihood(x, mu1, Sigma1) + log(pi1)
  delta2 <- log_likelihood(x, mu2, Sigma2) + log(pi2)
  ifelse(delta1 > delta2, 1, 2)
}

## Cross-validation for QDA
k_folds <- 10

## Combine data and create labels
class1_data <- Ering_class3
class2_data <- Ering_class2
X <- rbind(class1_data, class2_data)
y <- c(rep(1, nrow(class1_data)), rep(2, nrow(class2_data)))
data_with_labels <- data.frame(X, class = y)

## Create folds
folds <- createFolds(y = data_with_labels$class, k = k_folds, list = TRUE, returnTrain = FALSE)

## Initialize confusion matrix
cm <- matrix(0, nrow = 2, ncol = 2)

## Perform k-fold CV
for (i in 1:k_folds) {
  # Split data into training and testing
  test_indices <- folds[[i]]
  train_data <- data_with_labels[-test_indices, ]
  test_data <- data_with_labels[test_indices, ]
  
  # Separate by class
  train_class1 <- train_data[train_data$class == 1, -ncol(train_data)] %>% as.matrix()
  train_class2 <- train_data[train_data$class == 2, -ncol(train_data)] %>% as.matrix()
  
  # Estimate means and covariances
  mu1_train <- colMeans(train_class1)
  mu2_train <- colMeans(train_class2)
  
  cov1_train <- bandPPP(train_class1, k = 10, eps = 0.1)
  Sigma1_train <- matrix(as.numeric(estimate(cov1_train)), nrow = p, ncol = p)
  
  cov2_train <- bandPPP(train_class2, k = 10, eps = 0.1)
  Sigma2_train <- matrix(as.numeric(estimate(cov2_train)), nrow = p, ncol = p)
  
  # Classify test data
  test_X <- as.matrix(test_data[, -ncol(test_data)])
  true_class <- test_data$class
  
  # Apply QDA to each test sample
  pred_class <- apply(test_X, 1, function(x) {
    qda_classify(x, mu1_train, Sigma1_train, mu2_train, Sigma2_train)
  })
  
  # Update confusion matrix
  for (j in 1:length(pred_class)) {
    cm[true_class[j], pred_class[j]] <- cm[true_class[j], pred_class[j]] + 1
  }
}

## Normalize confusion matrix to get proportions
cm_norm <- cm / rowSums(cm)

## Create a formatted contingency table
contingency_table <- data.frame(
  "True_Class" = c("Gesture 1", "Gesture 2"),
  "Gesture_1" = round(cm_norm[,1], 4),
  "Gesture_2" = round(cm_norm[,2], 4),
  "Total" = c(1, 1)
)
print(contingency_table)