# attach library
library(bspcov)
library(mlmts)
library(purrr)
library(readr)
library(reshape2)
library(ggplot2)
library(gridExtra)


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


# visualization
vis_Sigma <- function(Sigma0, ind = 1) {
  ggSigma <- melt(Sigma0) %>% ggplot(aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    labs(title = paste0("Gesture ", ind), x = "i", y = "j") +
    scale_fill_gradient2(mid = "white")
  return(ggSigma)
}

g <- grid.arrange(
  plot(postmean1) + 
    scale_fill_gradient2(mid = "white") +
    labs(title = paste0("Gesture 1")),
  plot(postmean2) + scale_fill_gradient2(mid = "white") +
    labs(title = paste0("Gesture 2")),
  ncol = 2
)
ggsave(filename = "figs/cov_gestures.png", g, width = 12, height = 6)
