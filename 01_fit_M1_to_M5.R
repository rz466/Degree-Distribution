#!/usr/bin/env Rscript

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

log_msg("Starting 01_fit_M1_to_M5.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(rstan)
})

rstan_options(auto_write = TRUE)

all_de <- readRDS("norman_limma_all_perturbations.rds")

genes <- sort(unique(as.character(all_de$gene)))
perts <- setdiff(
  sort(unique(as.character(all_de$perturbation))),
  "control"
)

fdr_name <- "FDR_0.01"
fdr_label <- "0.01"
fdr_cutoff <- 0.01

log_msg("Constructing edges and degrees for FDR_0.01")

edge_df <- all_de %>%
  filter(
    q_value <= fdr_cutoff,
    gene != perturbation
  ) %>%
  distinct(gene, perturbation)

log_msg(
  "FDR 0.01: number of significant non-self edges = ",
  nrow(edge_df)
)

edge_mat <- sparseMatrix(
  i = match(edge_df$gene, genes),
  j = match(edge_df$perturbation, perts),
  x = rep(TRUE, nrow(edge_df)),
  dims = c(length(genes), length(perts)),
  dimnames = list(genes, perts)
)

d_j <- as.integer(Matrix::rowSums(edge_mat != 0))
d_k <- as.integer(Matrix::colSums(edge_mat != 0))

delta_kj <- t(as.matrix(edge_mat))
storage.mode(delta_kj) <- "integer"

stan_degree_data <- list(
  J = length(d_j),
  K = length(d_k),
  d_j = d_j,
  d_k = d_k
)

stan_delta_data <- list(
  K = nrow(delta_kj),
  J = ncol(delta_kj),
  delta = delta_kj
)

log_msg("Length of d_j: ", length(d_j))
log_msg("Length of d_k: ", length(d_k))
log_msg(
  "Dimension of delta_kj: ",
  paste(dim(delta_kj), collapse = " x ")
)

log_msg("Writing Stan model files for M1-M5")

dir.create("stan_models", showWarnings = FALSE)

writeLines(
  c(
'data {
  int<lower=1> J;
  int<lower=1> K;
  
  array[J] int<lower=0, upper=K> d_j;
  array[K] int<lower=0, upper=J> d_k;
}

transformed data {
  int n_unique_j;
  int n_unique_k;
  
  array[J] int unique_dj;
  array[J] int count_dj;
  
  array[K] int unique_dk;
  array[K] int count_dk;
  
  n_unique_j = 0;
  
  for (j in 1:J) {
    int found;
    found = 0;
    
    if (n_unique_j > 0) {
      for (u in 1:n_unique_j) {
        if (d_j[j] == unique_dj[u]) {
          count_dj[u] += 1;
          found = 1;
        }
      }
    }
    
    if (found == 0) {
      n_unique_j += 1;
      unique_dj[n_unique_j] = d_j[j];
      count_dj[n_unique_j] = 1;
    }
  }
  
  n_unique_k = 0;
  
  for (k in 1:K) {
    int found;
    found = 0;
    
    if (n_unique_k > 0) {
      for (u in 1:n_unique_k) {
        if (d_k[k] == unique_dk[u]) {
          count_dk[u] += 1;
          found = 1;
        }
      }
    }
    
    if (found == 0) {
      n_unique_k += 1;
      unique_dk[n_unique_k] = d_k[k];
      count_dk[n_unique_k] = 1;
    }
  }
}

parameters {
  real<lower=1e-6, upper=10> alpha_j;
  real<lower=1e-6, upper=10> beta_j;
  
  real<lower=1e-6, upper=10> alpha_k;
  real<lower=1e-6, upper=10> beta_k;
}

model {
  target += uniform_lpdf(alpha_j | 1e-6, 10);
  target += uniform_lpdf(beta_j  | 1e-6, 10);
  target += uniform_lpdf(alpha_k | 1e-6, 10);
  target += uniform_lpdf(beta_k  | 1e-6, 10);
  
  for (u in 1:n_unique_j) {
    target += count_dj[u] *
      beta_binomial_lpmf(unique_dj[u] | K, alpha_j, beta_j);
  }
  
  for (u in 1:n_unique_k) {
    target += count_dk[u] *
      beta_binomial_lpmf(unique_dk[u] | J, alpha_k, beta_k);
  }
}',
    ""
  ),
  "stan_models/beta_binomial.stan"
)

writeLines(
  c(
'data {
  int<lower=1> J;
  int<lower=1> K;
  
  array[J] int<lower=0, upper=K> d_j;
  array[K] int<lower=0, upper=J> d_k;
}

transformed data {
  int sum_dj;
  int sum_dk;
  
  sum_dj = 0;
  sum_dk = 0;
  
  for (j in 1:J) {
    sum_dj += d_j[j];
  }
  
  for (k in 1:K) {
    sum_dk += d_k[k];
  }
}

parameters {
  real<lower=1e-9, upper=0.999999999> r_j;
  real<lower=1e-9, upper=0.999999999> r_k;
}

model {
  real log_r_j;
  real log_r_k;
  real log_norm_dj;
  real log_norm_dk;
  
  target += uniform_lpdf(r_j | 1e-9, 0.999999999);
  target += uniform_lpdf(r_k | 1e-9, 0.999999999);
  
  log_r_j = log(r_j);
  log_r_k = log(r_k);
  
  log_norm_dj =
    log1m_exp((K + 1) * log_r_j) -
    log1m_exp(log_r_j);
  
  log_norm_dk =
    log1m_exp((J + 1) * log_r_k) -
    log1m_exp(log_r_k);
  
  target += sum_dj * log_r_j - J * log_norm_dj;
  target += sum_dk * log_r_k - K * log_norm_dk;
}',
    ""
  ),
  "stan_models/truncated_geometric.stan"
)

writeLines(
  c(
'data {
  int<lower=1> K;
  int<lower=1> J;
  array[K, J] int<lower=0, upper=1> delta;
}

transformed data {
  array[J] int target_degree;
  array[K + 1] int target_degree_count;
  
  for (k in 1:(K + 1)) {
    target_degree_count[k] = 0;
  }
  
  for (j in 1:J) {
    target_degree[j] = 0;
    
    for (k in 1:K) {
      target_degree[j] += delta[k, j];
    }
    
    target_degree_count[target_degree[j] + 1] += 1;
  }
}

parameters {
  real<lower=1e-6, upper=10> alpha;
  real<lower=1e-6, upper=10> beta;
}

model {
  target += uniform_lpdf(alpha | 1e-6, 10);
  target += uniform_lpdf(beta | 1e-6, 10);
  
  for (k in 0:K) {
    target += target_degree_count[k + 1] *
      (
        lbeta(k + alpha, K - k + beta) -
          lbeta(alpha, beta)
      );
  }
}',
    ""
  ),
  "stan_models/edge_bernoulli_psi_j.stan"
)

writeLines(
  c(
'data {
  int<lower=1> K;
  int<lower=1> J;
  array[K, J] int<lower=0, upper=1> delta;
}

transformed data {
  array[K] int perturbation_degree;
  int n_unique;
  array[K] int unique_degree;
  array[K] int degree_count;
  
  n_unique = 0;
  
  for (k in 1:K) {
    int found;
    
    perturbation_degree[k] = 0;
    
    for (j in 1:J) {
      perturbation_degree[k] += delta[k, j];
    }
    
    found = 0;
    
    if (n_unique > 0) {
      for (u in 1:n_unique) {
        if (
          perturbation_degree[k] ==
          unique_degree[u]
        ) {
          degree_count[u] += 1;
          found = 1;
        }
      }
    }
    
    if (found == 0) {
      n_unique += 1;
      unique_degree[n_unique] =
        perturbation_degree[k];
      degree_count[n_unique] = 1;
    }
  }
}

parameters {
  real<lower=1e-6, upper=10> alpha;
  real<lower=1e-6, upper=10> beta;
}

model {
  target += uniform_lpdf(alpha | 1e-6, 10);
  target += uniform_lpdf(beta | 1e-6, 10);
  
  for (u in 1:n_unique) {
    target += degree_count[u] *
      (
        lbeta(
          unique_degree[u] + alpha,
          J - unique_degree[u] + beta
        ) -
          lbeta(alpha, beta)
      );
  }
}',
    ""
  ),
  "stan_models/edge_bernoulli_psi_k.stan"
)

writeLines('
data {
  int<lower=1> K;
  int<lower=1> J;
  array[K, J] int<lower=0, upper=1> delta;
}

parameters {
  real<lower=1e-6, upper=10> eta;
  real<lower=1e-6, upper=10> zeta;

  real<lower=1e-6, upper=10> chi;
  real<lower=1e-6, upper=10> omega;

  vector<lower=0, upper=1>[J] psi_j;
  vector<lower=0, upper=1>[K] psi_k;
}

model {
  target += uniform_lpdf(eta  | 1e-6, 10);
  target += uniform_lpdf(zeta | 1e-6, 10);

  target += uniform_lpdf(chi   | 1e-6, 10);
  target += uniform_lpdf(omega | 1e-6, 10);

  target += beta_lpdf(psi_j | eta, zeta);
  target += beta_lpdf(psi_k | chi, omega);

  for (k in 1:K) {
    target += bernoulli_lpmf(
      delta[k] | psi_k[k] * psi_j
    );
  }
}
', "stan_models/edge_bernoulli_psi_k_times_psi_j.stan")

log_msg("Compiling Stan models M1-M5")

stan_m1 <- stan_model("stan_models/beta_binomial.stan")

stan_m2 <- stan_model("stan_models/truncated_geometric.stan")

stan_m3 <- stan_model("stan_models/edge_bernoulli_psi_j.stan")

stan_m4 <- stan_model("stan_models/edge_bernoulli_psi_k.stan")

stan_m5 <- stan_model("stan_models/edge_bernoulli_psi_k_times_psi_j.stan")

log_msg("Sampling M1")
fit_m1 <- sampling(stan_m1, data = stan_degree_data, iter = 50000, warmup = 1000, chains = 3, seed = 123, cores = 3)

log_msg("Sampling M2")
fit_m2 <- sampling(stan_m2, data = stan_degree_data, iter = 50000, warmup = 1000, chains = 3, seed = 123, cores = 3)

log_msg("Sampling M3")
fit_m3 <- sampling(stan_m3, data = stan_delta_data, iter = 50000, warmup = 1000, chains = 3, seed = 123, cores = 3)

log_msg("Sampling M4")
fit_m4 <- sampling(stan_m4, data = stan_delta_data, iter = 50000, warmup = 1000, chains = 3, seed = 123, cores = 3)

log_msg("Sampling M5")
fit_m5 <- sampling(stan_m5,data = stan_delta_data, iter = 50000, warmup = 1000, chains = 3, seed = 123, cores = 3)

saveRDS(
  list(
    fdr_name = fdr_name,
    fdr_label = fdr_label,
    d_j = d_j,
    d_k = d_k,
    delta_kj = delta_kj,
    fit_m1 = fit_m1,
    fit_m2 = fit_m2,
    fit_m3 = fit_m3,
    fit_m4 = fit_m4,
    fit_m5 = fit_m5
  ),
  "fit_results_M1_to_M5.rds"
)

log_msg("Saved fit_results_M1_to_M5.rds")
log_msg("01_fit_M1_to_M5.R finished")


