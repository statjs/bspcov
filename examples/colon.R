# attach library
library(bspcov)

# fix random seed
set.seed(1234)

# load and preprocess data
data("colon")
data("tissues")
colon_data <- proc_colon(colon, tissues = tissues)
X <- colon_data$X

# beta-mixture prior (parallel processing automatically set up when do.parallel = TRUE)
foo.bmspcov <- bmspcov(X,
  Sigma = cov(X), nchain = 4,
  nsample = list(burnin = 1000, nmc = 2000),
  seed = 1234, do.parallel = TRUE
)

# screened beta-mixture prior (parallel processing automatically set up when do.parallel = TRUE)
cutoff <- list(method = "corr", rho = 0.25)
foo.sbmspcov <- sbmspcov(
  X = X, Sigma = cov(X), cutoff, nchain = 4,
  nsample = list(burnin = 1000, nmc = 2000),
  seed = 1234, do.parallel = TRUE
)

# thresPPP
foo.thres.cv <- cv.thresPPP(X,
  thresvec = c(0.001, 0.1, 0.2, 0.3, 0.5, 1.0),
  epsvec = c(0.01, 0.05), nsample = 100, ncores = 4
)
foo.thres <- thresPPP(X = X, thres = foo.thres.cv$error$thres[1], eps = foo.thres.cv$error$eps[1], nsample = 100)

# posterior mean estimate
sigmah.bmspcov <- estimate(foo.bmspcov)
sigmah.sbmspcov <- estimate(foo.sbmspcov)
sigmah.thres <- estimate(foo.thres)


# ------- Analysis -------
saveRDS(list(foo.bmspcov, foo.sbmspcov, foo.thres), "examples/results/colon_bspcov.rds")
colon_bspcov <- readRDS("examples/results/colon_bspcov.rds")
foo.bmspcov <- colon_bspcov[[1]]
foo.sbmspcov <- colon_bspcov[[2]]
foo.thres <- colon_bspcov[[3]]
sigmah.bmspcov <- estimate(foo.bmspcov)
sigmah.sbmspcov <- estimate(foo.sbmspcov)
sigmah.thres <- estimate(foo.thres)
quantile.bmspcov <- quantile(foo.bmspcov, probs = c(0.025, 0.5, 0.975))
quantile.sbmspcov <- quantile(foo.sbmspcov, probs = c(0.025, 0.5, 0.975))
quantile.thres <- quantile(foo.thres, probs = c(0.025, 0.5, 0.975))

library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)
library(dplyr)
library(caret)
library(mvnfast)
library(rstan)
library(rjags)
library(R2jags)
library(tictoc)

# summary for diagnostics - save summary to plain text file
summary(foo.bmspcov, cols = 24:26, rows = 24:26) %>%
  capture.output() %>%
  writeLines("examples/results/colon_bmspcov_summary.txt")
summary(foo.sbmspcov, cols = 24:26, rows = 24:26) %>%
  capture.output() %>%
  writeLines("examples/results/colon_sbmspcov_summary.txt")

# traceplots
plot(foo.bmspcov, cols = 24:26, rows = 24:26)
ggsave("figs/colon_bmspcov_trace.png", width = 12, height = 7)

plot(foo.sbmspcov, cols = 24:26, rows = 24:26)
ggsave("figs/colon_sbmspcov_trace.png", width = 12, height = 7)

# Set options for rstan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# visualization

## compute range for color scale
min_sigmah <- min(c(as.numeric(sigmah.bmspcov), as.numeric(sigmah.sbmspcov), as.numeric(sigmah.thres),
  sapply(quantile.bmspcov, min),
  sapply(quantile.sbmspcov, min),
  sapply(quantile.thres, min)
))
max_sigmah <- max(c(as.numeric(sigmah.bmspcov), as.numeric(sigmah.sbmspcov), as.numeric(sigmah.thres),
  sapply(quantile.bmspcov, max),
  sapply(quantile.sbmspcov, max),
  sapply(quantile.thres, max)
))
cont_scale <- scale_fill_gradient(
  name   = "",
  low    = "black",
  high   = "white",
  limits = c(min_sigmah, max_sigmah),
  guide  = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(0.5, "cm"))
)

## plots
p1 <- plot(sigmah.bmspcov, title = "bmspcov", x_label = "Gene", y_label = "Gene")
p2 <- plot(sigmah.sbmspcov, title = "sbmspcov", x_label = "Gene", y_label = "Gene")
p3 <- plot(sigmah.thres, title = "thresPPP", x_label = "Gene", y_label = "Gene")
g <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect") &
  cont_scale &
  theme(legend.position = "bottom")
ggsave("figs/colon.png", g, width = 15, height = 6)

## plots without thres
g_no_thres <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  cont_scale &
  theme(legend.position = "bottom")
ggsave(filename = "figs/colon_no_thres.png", width = 12, height = 7)

## visualize the uncertainty using new quantile function
g_bmspcov_u <- plot(quantile.bmspcov,
  titles = c("bmspcov (2.5%)", "bmspcov (50%)", "bmspcov (97.5%)"),
  x_label = "Gene", y_label = "Gene",
  ncol = 3
) & cont_scale
ggsave("figs/colon_bmspcov_u.png", g_bmspcov_u, width = 14, height = 6)

g_sbmspcov_u <- plot(quantile.sbmspcov,
  titles = c("sbmspcov (2.5%)", "sbmspcov (50%)", "sbmspcov (97.5%)"),
  x_label = "Gene", y_label = "Gene",
  ncol = 3
) & cont_scale
ggsave("figs/colon_sbmspcov_u.png", g_sbmspcov_u, width = 14, height = 6)

g_thres_u <- plot(quantile.thres,
  titles = c("thresPPP (2.5%)", "thresPPP (50%)", "thresPPP (97.5%)"),
  x_label = "Gene", y_label = "Gene",
  ncol = 3
) & cont_scale
ggsave("figs/colon_thresPPP_u.png", g_thres_u, width = 14, height = 6)

# ------- Stan and JAGS Implementation -------

# Prepare data dimensions
n <- nrow(X)
p <- ncol(X)
# Stan implementation
stan_code <- "
data {
  int<lower=0> n; // number of observations
  int<lower=0> p; // number of variables
  matrix[n, p] X; // data matrix
  vector[p] mu; // mean vector
}
parameters {
  cov_matrix[p] Sigma; // covariance matrix
}
model {
  Sigma ~ inv_wishart(p, diag_matrix(rep_vector(1.0, p)));
  for (i in 1:n) {
    X[i, ] ~ multi_normal(mu, Sigma);
  }
}
"
stan_data <- list(n = n, p = p, X = X, mu = colMeans(X))
stan_model <- rstan::stan_model(model_code = stan_code)

# JAGS implementation
jags_code <- "
model {
  # Priors
  Omega ~ dwish(R, p)  # Wishart prior for precision matrix with p+1 degrees of freedom
  Sigma <- inverse(Omega)    # Covariance matrix

  # Likelihood
  for (i in 1:n) {
    X[i, 1:p] ~ dmnorm(mu, Omega)
  }
}
"
# Write JAGS model to file
jags_model_file <- tempfile(fileext = ".txt")
writeLines(jags_code, jags_model_file)
# Prepare data for JAGS
jags_data <- list(
  n = n,
  p = p,
  X = X,
  mu = colMeans(X),
  R = diag(p)
)
jags_params <- c("Sigma")

## Cross-validation for LDA with runtime tracking
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
  tic("bmspcov start")
  foo.bmspcov <- bmspcov(train, Sigma = cov(train), nchain = 4, seed = 1234, do.parallel = TRUE)
  bmspcov.elapsed <- toc(quiet = TRUE)
  bmspcov.elapsed <- bmspcov.elapsed$toc - bmspcov.elapsed$tic
  tic("sbmspcov start")
  foo.sbmspcov <- sbmspcov(train, Sigma = cov(train), cutoff = cutoff, nchain = 4, seed = 1234, do.parallel = TRUE)
  sbmspcov.elapsed <- toc(quiet = TRUE)
  sbmspcov.elapsed <- sbmspcov.elapsed$toc - sbmspcov.elapsed$tic
  tic("thresPPP start")
  foo.thres <- thresPPP(train, thres = foo.thres.cv$error$thres[1], eps = foo.thres.cv$error$eps[1], nsample = 100)
  thresPPP.elapsed <- toc(quiet = TRUE)
  thresPPP.elapsed <- thresPPP.elapsed$toc - thresPPP.elapsed$tic
  tic("Stan start")
  foo.stan <- rstan::sampling(
    stan_model,
    data = list(n = nrow(train), p = ncol(train), X = train, mu = colMeans(train)),
    iter = 3000,
    warmup = 1000,
    chains = 4,
    seed = 1234,
    init = function() {
      list(Sigma = cov(train) + 1e-4 * diag(ncol(train))) # Initial value for Sigma
    }
  )
  stan.elapsed <- toc(quiet = TRUE)
  stan.elapsed <- stan.elapsed$toc - stan.elapsed$tic
  tic("JAGS start")
  foo.jags <- R2jags::jags(
    data = list(n = nrow(train), p = ncol(train), X = train, mu = colMeans(train), R = diag(ncol(train))),
    parameters.to.save = jags_params,
    model.file = jags_model_file,
    n.chains = 4,
    n.iter = 3000,
    n.burnin = 1000,
    n.thin = 1
  )
  jags.elapsed <- toc(quiet = TRUE)
  jags.elapsed <- jags.elapsed$toc - jags.elapsed$tic

  sigmah.bmspcov <- estimate(foo.bmspcov)
  sigmah.sbmspcov <- estimate(foo.sbmspcov)
  sigmah.thres <- estimate(foo.thres)
  est.sig_stan <- apply(rstan::extract(foo.stan, "Sigma")$Sigma, c(2, 3), mean)
  est.sig_jags <- apply(foo.jags$BUGSoutput$sims.list$Sigma, c(2, 3), mean)

  count1 <- sum(mvnfast::dmvn(test1, colMeans(train1), sigmah.bmspcov, log = TRUE) >
    mvnfast::dmvn(test1, colMeans(train2), sigmah.bmspcov, log = TRUE))
  count2 <- sum(mvnfast::dmvn(test2, colMeans(train1), sigmah.bmspcov, log = TRUE) <
    mvnfast::dmvn(test2, colMeans(train2), sigmah.bmspcov, log = TRUE))
  count <- count1 + count2
  prob <- count / (nrow(test1) + nrow(test2))

  count1.two <- sum(mvnfast::dmvn(test1, colMeans(train1), sigmah.sbmspcov, log = TRUE) >
    mvnfast::dmvn(test1, colMeans(train2), sigmah.sbmspcov, log = TRUE))
  count2.two <- sum(mvnfast::dmvn(test2, colMeans(train1), sigmah.sbmspcov, log = TRUE) <
    mvnfast::dmvn(test2, colMeans(train2), sigmah.sbmspcov, log = TRUE))
  count.two <- count1.two + count2.two
  prob.two <- count.two / (nrow(test1) + nrow(test2))

  count1.thres <- sum(mvnfast::dmvn(test1, colMeans(train1), sigmah.thres, log = TRUE) >
    mvnfast::dmvn(test1, colMeans(train2), sigmah.thres, log = TRUE))
  count2.thres <- sum(mvnfast::dmvn(test2, colMeans(train1), sigmah.thres, log = TRUE) <
    mvnfast::dmvn(test2, colMeans(train2), sigmah.thres, log = TRUE))
  count.thres <- count1.thres + count2.thres
  prob.thres <- count.thres / (nrow(test1) + nrow(test2))

  count1.stan <- sum(mvnfast::dmvn(test1, colMeans(train1), est.sig_stan, log = TRUE) >
    mvnfast::dmvn(test1, colMeans(train2), est.sig_stan, log = TRUE))
  count2.stan <- sum(mvnfast::dmvn(test2, colMeans(train1), est.sig_stan, log = TRUE) <
    mvnfast::dmvn(test2, colMeans(train2), est.sig_stan, log = TRUE))
  count.stan <- count1.stan + count2.stan
  prob.stan <- count.stan / (nrow(test1) + nrow(test2))

  count1.jags <- sum(mvnfast::dmvn(test1, colMeans(train1), est.sig_jags, log = TRUE) >
    mvnfast::dmvn(test1, colMeans(train2), est.sig_jags, log = TRUE))
  count2.jags <- sum(mvnfast::dmvn(test2, colMeans(train1), est.sig_jags, log = TRUE) <
    mvnfast::dmvn(test2, colMeans(train2), est.sig_jags, log = TRUE))
  count.jags <- count1.jags + count2.jags
  prob.jags <- count.jags / (nrow(test1) + nrow(test2))

  list(
    method = c("bmspcov", "sbmspcov", "thresPPP", "Stan", "JAGS"),
    accuracy = c(prob, prob.two, prob.thres, prob.stan, prob.jags),
    elapsed = c(
      bmspcov.elapsed, sbmspcov.elapsed, thresPPP.elapsed,
      stan.elapsed, jags.elapsed
    )
  )
}
res <- purrr::map2(folds[[1]], folds[[2]], ~ testacc(.x, .y))
saveRDS(res, "examples/results/colon_cv.rds")

res <- readRDS("examples/results/colon_cv.rds")
summary_table <- lapply(res, data.frame) %>%
  do.call("rbind", .) %>%
  group_by(method) %>%
  mutate(misclass = 1 - round(accuracy, 4)) %>%
  summarise(
    misclass_mean = mean(misclass), misclass_sd = sd(misclass),
    elapsed_mean = mean(elapsed), elapsed_sd = sd(elapsed)
  ) %>%
  arrange(elapsed_mean)
summary_table

# print summary table in latex format with bold formatting for best values
summary_table_formatted <- summary_table %>%
  mutate(
    misclass_combined = ifelse(misclass_mean == min(misclass_mean), 
                              paste0("\\textbf{", sprintf("%.3f", misclass_mean), " (", sprintf("%.3f", misclass_sd), ")}"), 
                              paste0(sprintf("%.3f", misclass_mean), " (", sprintf("%.3f", misclass_sd), ")")),
    elapsed_combined = ifelse(elapsed_mean == min(elapsed_mean), 
                             paste0("\\textbf{", sprintf("%.3f", elapsed_mean), " (", sprintf("%.3f", elapsed_sd), ")}"), 
                             paste0(sprintf("%.3f", elapsed_mean), " (", sprintf("%.3f", elapsed_sd), ")"))
  ) %>%
  select(method, misclass_combined, elapsed_combined)
# Rename columns for better presentation
colnames(summary_table_formatted) <- c("Method", "Misclassification Rate", "Elapsed Time (sec)")

summary_kable <- knitr::kable(summary_table_formatted, format = "latex", booktabs = TRUE, escape = FALSE)
# save to file examples/results/colon_summary.tex
writeLines(summary_kable, "examples/results/colon_summary.tex")

# Clean up temporary file
unlink(jags_model_file)
