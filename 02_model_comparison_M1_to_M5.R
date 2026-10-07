#!/usr/bin/env Rscript

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

log_msg("Starting 02_model_comparison_M1_to_M5.R")

suppressPackageStartupMessages({
  library(rstan)
  library(bridgesampling)
  library(dplyr)
  library(tibble)
})

rstan_options(auto_write = TRUE)

log_msg("Loading fit_results_M1_to_M5.rds")

fit_results <- readRDS("fit_results_M1_to_M5.rds")

fdr_name <- fit_results$fdr_name
fdr_label <- fit_results$fdr_label
d_j <- fit_results$d_j
d_k <- fit_results$d_k
delta_kj <- fit_results$delta_kj

fit_m1 <- fit_results$fit_m1
fit_m2 <- fit_results$fit_m2
fit_m3 <- fit_results$fit_m3
fit_m4 <- fit_results$fit_m4
fit_m5 <- fit_results$fit_m5

stan_degree_data <- list(
  J = length(d_j),
  K = length(d_k),
  d_j = as.integer(d_j),
  d_k = as.integer(d_k)
)

storage.mode(delta_kj) <- "integer"

stan_delta_data <- list(
  K = nrow(delta_kj),
  J = ncol(delta_kj),
  delta = delta_kj
)

# Compile the same Stan models used to create fit_m1-fit_m5
log_msg("Compiling Stan models M1-M5")

stan_m1 <- stan_model(
  "stan_models/beta_binomial.stan"
)

stan_m2 <- stan_model(
  "stan_models/truncated_geometric.stan"
)

stan_m3 <- stan_model(
  "stan_models/edge_bernoulli_psi_j.stan"
)

stan_m4 <- stan_model(
  "stan_models/edge_bernoulli_psi_k.stan"
)

stan_m5 <- stan_model(
  "stan_models/edge_bernoulli_psi_k_times_psi_j.stan"
)

# Create zero-chain Stan objects for evaluating the log posterior
log_msg("Creating zero-chain Stan model objects")

stanfit_m1 <- suppressWarnings(
  sampling(
    stan_m1,
    data = stan_degree_data,
    chains = 0
  )
)

stanfit_m2 <- suppressWarnings(
  sampling(
    stan_m2,
    data = stan_degree_data,
    chains = 0
  )
)

stanfit_m3 <- suppressWarnings(
  sampling(
    stan_m3,
    data = stan_delta_data,
    chains = 0
  )
)

stanfit_m4 <- suppressWarnings(
  sampling(
    stan_m4,
    data = stan_delta_data,
    chains = 0
  )
)

stanfit_m5 <- suppressWarnings(
  sampling(
    stan_m5,
    data = stan_delta_data,
    chains = 0
  )
)

# Marginal likelihoods p(Y | M)
log_msg("Bridge sampling M1")

bridge_m1 <- bridge_sampler(
  fit_m1,
  stanfit_model = stanfit_m1,
  silent = TRUE
)

log_msg("Bridge sampling M2")

bridge_m2 <- bridge_sampler(
  fit_m2,
  stanfit_model = stanfit_m2,
  silent = TRUE
)

log_msg("Bridge sampling M3")

bridge_m3 <- bridge_sampler(
  fit_m3,
  stanfit_model = stanfit_m3,
  silent = TRUE
)

log_msg("Bridge sampling M4")

bridge_m4 <- bridge_sampler(
  fit_m4,
  stanfit_model = stanfit_m4,
  silent = TRUE
)

# Try normal bridge sampling for M5
log_msg("Trying normal bridge sampling M5")

bridge_m5 <- try(
  bridge_sampler(
    fit_m5,
    stanfit_model = stanfit_m5,
    method = "normal",
    silent = TRUE
  ),
  silent = TRUE
)

m5_method <- "normal"

# If normal fails, use Warp-III with 10 repetitions
if (
  inherits(bridge_m5, "try-error") ||
  !is.finite(bridge_m5$logml)
) {
  
  log_msg(
    "Normal bridge sampling failed for M5; ",
    "running Warp-III with 10 repetitions"
  )
  
  bridge_m5 <- bridge_sampler(
    fit_m5,
    stanfit_model = stanfit_m5,
    method = "warp3",
    repetitions = 10,
    silent = TRUE
  )
  
  m5_method <- "warp3, 10 repetitions"
}

# Use the normal estimate or median of the 10 Warp-III estimates
if (m5_method == "normal") {
  m5_logml <- bridge_m5$logml
} else {
  m5_logml <- logml(bridge_m5)
}

# log p(Y | M)
logml_values <- c(
  M1_BetaBinomial =
    bridge_m1$logml,
  
  M2_TruncatedGeometric =
    bridge_m2$logml,
  
  M3_EdgeBernoulli_psi_j =
    bridge_m3$logml,
  
  M4_EdgeBernoulli_psi_k =
    bridge_m4$logml,
  
  M5_EdgeBernoulli_psi_k_times_psi_j =
    m5_logml
)

# Bridge sampling percentage error
bridge_error_m1 <-
  error_measures(bridge_m1)$percentage

bridge_error_m2 <-
  error_measures(bridge_m2)$percentage

bridge_error_m3 <-
  error_measures(bridge_m3)$percentage

bridge_error_m4 <-
  error_measures(bridge_m4)$percentage

if (m5_method == "normal") {
  
  bridge_error_m5 <-
    error_measures(bridge_m5)$percentage
  
} else {
  
  m5_error <- error_measures(bridge_m5)
  
  bridge_error_m5 <- paste0(
    "min = ", m5_error$min,
    "; max = ", m5_error$max,
    "; IQR = ", m5_error$IQR
  )
}

bridge_error <- c(
  M1_BetaBinomial =
    as.character(bridge_error_m1),
  
  M2_TruncatedGeometric =
    as.character(bridge_error_m2),
  
  M3_EdgeBernoulli_psi_j =
    as.character(bridge_error_m3),
  
  M4_EdgeBernoulli_psi_k =
    as.character(bridge_error_m4),
  
  M5_EdgeBernoulli_psi_k_times_psi_j =
    as.character(bridge_error_m5)
)

# log p(M) with equal prior probability
log_prior <- log(
  rep(
    1 / length(logml_values),
    length(logml_values)
  )
)

names(log_prior) <- names(logml_values)

# log(p(Y | M)p(M))
log_bayes_numerator <- logml_values + log_prior

# log p(Y)
log_bayes_denominator <-
  max(log_bayes_numerator) +
  log(
    sum(
      exp(
        log_bayes_numerator -
          max(log_bayes_numerator)
      )
    )
  )

# p(M | Y)
posterior_model_probability <- exp(
  log_bayes_numerator -
    log_bayes_denominator
)

# One comparison table
model_comparison <- tibble(
  model = names(logml_values),
  
  data_used = c(
    "(d_j, d_k)",
    "(d_j, d_k)",
    "D",
    "D",
    "D"
  ),
  
  model_family = c(
    "M1: Beta-Binomial",
    "M2: Truncated geometric",
    "M3: Bernoulli psi_j",
    "M4: Bernoulli psi_k",
    "M5: Bernoulli psi_k * psi_j"
  ),
  
  bridge_method = c(
    "normal",
    "normal",
    "normal",
    "normal",
    m5_method
  ),
  
  log_marginal_likelihood =
    as.numeric(logml_values),
  
  bridge_percentage_error =
    as.character(bridge_error),
  
  prior_model_probability =
    exp(log_prior),
  
  posterior_model_probability =
    as.numeric(posterior_model_probability)
) %>%
  arrange(
    desc(log_marginal_likelihood)
  )

print(
  model_comparison,
  n = Inf,
  width = Inf
)

write.csv(
  model_comparison,
  "model_comparison_M1_to_M5.csv",
  row.names = FALSE
)

log_msg("Saved model_comparison_M1_to_M5.csv")
log_msg("02_model_comparison_M1_to_M5.R finished")