#!/usr/bin/env Rscript

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

log_msg("Starting 03_plot_ppc_M1_to_M5.R")

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

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

storage.mode(delta_kj) <- "integer"

output_dir <- paste0("ppc_M1_to_M5_", fdr_name)
dir.create(output_dir, showWarnings = FALSE)

pdf_file <- file.path(
  output_dir,
  paste0("PPC_M1_to_M5_", fdr_name, ".pdf")
)

# Posterior predictive checks
# For M1-M2: simulate posterior predictive degree vectors for d_j and d_k
set.seed(123)

# Posterior samples used for degree
ppc_degree_sample <- 500

# Plotting function for both d_j and d_k ppc
plot_ppc_degree <- function(obs_degree, sim_degree, model_name, degree_type,
                            max_degree, bins = 30) {
  
  # Stack posterior predictive degrees across posterior samples
  stacked_sim_degree <- as.vector(sim_degree)
  
  if (degree_type == "Perturbation degree d_k") {
    
    # Create common bins
    breaks <- seq(-0.5, max_degree + 0.5, length.out = bins + 1)
    
    # Calculate observed histogram counts
    obs_hist <- hist(obs_degree, breaks = breaks, plot = FALSE, include.lowest = TRUE)
    
    # Calculate posterior predictive histogram counts
    sim_hist <- hist(stacked_sim_degree, breaks = breaks, plot = FALSE, include.lowest = TRUE)
    
    # Convert histogram counts into probabilities
    plot_df <- bind_rows(
      tibble(
        bin_mid = obs_hist$mids,
        prob = obs_hist$counts / sum(obs_hist$counts),
        source = "Observed"
      ),
      tibble(
        bin_mid = sim_hist$mids,
        prob = sim_hist$counts / sum(sim_hist$counts),
        source = "Posterior predictive"
      )
    )
    
    # Plot probability histogram for d_k
    ggplot(plot_df, aes(x = bin_mid, y = prob, fill = source)) +
      geom_col(position = "identity", alpha = 0.45, width = diff(breaks)[1]) +
      labs(
        title = paste0(model_name, ", FDR = ", fdr_label),
        x = "Degree",
        y = "Probability",
        fill = NULL
      ) +
      coord_cartesian(
        xlim = c(-0.5, max_degree + 0.5),
        ylim = c(0, max(plot_df$prob, na.rm = TRUE))
      ) +
      theme_minimal()
    
  } else {
    
    # Convert degrees into probabilities
    plot_df <- bind_rows(
      tibble(degree = obs_degree, source = "Observed"),
      tibble(degree = stacked_sim_degree, source = "Posterior predictive")
    ) %>%
      count(source, degree, name = "freq") %>%
      group_by(source) %>%
      mutate(prob = freq / sum(freq)) %>%
      ungroup()
    
    # Plot exact degree probability bars for d_j
    ggplot(plot_df, aes(x = degree, y = prob, fill = source)) +
      geom_col(position = "identity", alpha = 0.45, width = 1) +
      labs(
        title = paste0(model_name, ", FDR = ", fdr_label),
        x = "Degree",
        y = "Probability",
        fill = NULL
      ) +
      coord_cartesian(
        xlim = c(-0.5, max_degree + 0.5),
        ylim = c(0, max(plot_df$prob, na.rm = TRUE))
      ) +
      theme_minimal()
  }
}

degree_model <- list(
  list(model_name = "M1(d_j): Beta-Binomial", model_type = "M1", fit = fit_m1,
       obs_degree = d_j, n_obs = length(d_j), max_degree = length(d_k),
       degree_type = "Target degree d_j"),
  
  list(model_name = "M2(d_j): Truncated geometric", model_type = "M2", fit = fit_m2,
       obs_degree = d_j, n_obs = length(d_j), max_degree = length(d_k),
       degree_type = "Target degree d_j"),
  
  list(model_name = "M1(d_k): Beta-Binomial", model_type = "M1", fit = fit_m1,
       obs_degree = d_k, n_obs = length(d_k), max_degree = length(d_j),
       degree_type = "Perturbation degree d_k"),
  
  list(model_name = "M2(d_k): Truncated geometric", model_type = "M2", fit = fit_m2,
       obs_degree = d_k, n_obs = length(d_k), max_degree = length(d_j),
       degree_type = "Perturbation degree d_k")
)

log_msg("Creating PPC PDF: ", pdf_file)
pdf(pdf_file, width = 10, height = 7)

for (model in degree_model) {
  
  log_msg("Posterior predictive simulation for ", model$model_name)
  
  # Define possible degree values for this model
  degree_grid <- 0:model$max_degree
  
  if (model$model_type == "M1") {
    
    # Extract posterior draws of alpha and beta
    if (model$degree_type == "Target degree d_j") {
      draws <- rstan::extract(model$fit, pars = c("alpha_j", "beta_j"))
      alpha_draws <- draws$alpha_j
      beta_draws <- draws$beta_j
    } else {
      draws <- rstan::extract(model$fit, pars = c("alpha_k", "beta_k"))
      alpha_draws <- draws$alpha_k
      beta_draws <- draws$beta_k
    }
    
    # Randomly select posterior samples
    idx <- sample(seq_along(alpha_draws), ppc_degree_sample)
    
    # Create matrix to store simulated degree vectors
    sim_degree <- matrix(NA, nrow = length(idx), ncol = model$n_obs)
    
    # Loop over selected posterior samples
    for (s in seq_along(idx)) {
      
      # Get alpha and beta from posterior sample
      alpha_s <- alpha_draws[idx[s]]
      beta_s <- beta_draws[idx[s]]
      
      # Simulate edge probabilities and degrees
      pi_s <- rbeta(model$n_obs, shape1 = alpha_s, shape2 = beta_s)
      sim_degree[s, ] <- rbinom(model$n_obs, size = model$max_degree, prob = pi_s)
    }
  }
  
  if (model$model_type == "M2") {
    
    # Extract posterior draws of r
    if (model$degree_type == "Target degree d_j") {
      draws <- rstan::extract(model$fit, pars = "r_j")
      r_draws <- draws$r_j
    } else {
      draws <- rstan::extract(model$fit, pars = "r_k")
      r_draws <- draws$r_k
    }
    
    # Randomly select posterior samples
    idx <- sample(seq_along(r_draws), ppc_degree_sample)
    
    # Create matrix to store simulated degree vectors
    sim_degree <- matrix(NA, nrow = length(idx), ncol = model$n_obs)
    
    # Loop over selected posterior samples
    for (s in seq_along(idx)) {
      
      # Simulate degrees
      r_s <- r_draws[idx[s]]
      log_prob_s <- degree_grid * log(r_s)
      prob_s <- exp(log_prob_s - max(log_prob_s))
      prob_s <- prob_s / sum(prob_s)
      sim_degree[s, ] <- sample(degree_grid, size = model$n_obs, replace = TRUE, prob = prob_s)
    }
  }
  
  # Plot observed vs stacked posterior predictive degree distribution
  print(plot_ppc_degree(
    obs_degree = model$obs_degree,
    sim_degree = sim_degree,
    model_name = model$model_name,
    degree_type = model$degree_type,
    max_degree = model$max_degree
  ))
}

# For M3-M5: simulate posterior predictive edge matrices D, then compute d_j and d_k
# Posterior sample used
ppc_edge_sample <- 2000

# Edge matrix dimensions
K_edge <- nrow(delta_kj)
J_edge <- ncol(delta_kj)

# Compute observed d_j and d_k
observed_dj <- as.integer(colSums(delta_kj))
observed_dk <- as.integer(rowSums(delta_kj))

edge_model <- list(
  list(model_name = "M3", model_type = "M3", model_label = "Bernoulli psi_j", fit = fit_m3),
  list(model_name = "M4", model_type = "M4", model_label = "Bernoulli psi_k", fit = fit_m4),
  list(model_name = "M5", model_type = "M5", model_label = "Multiplicative Bernoulli psi_k × psi_j", fit = fit_m5)
)

for (model in edge_model) {
  
  log_msg("Posterior predictive simulation for ", model$model_name)
  
  # Extract posterior draws
  if (model$model_type == "M5") {
    
    draws <- rstan::extract(model$fit, pars = c("psi_j", "psi_k"))
    idx <- sample(seq_len(nrow(draws$psi_j)), size = ppc_edge_sample, replace = FALSE)
    
  } else {
    
    draws <- rstan::extract(model$fit, pars = c("alpha", "beta"))
    idx <- sample(seq_along(draws$alpha), size = ppc_edge_sample, replace = FALSE)
    
  }
  
  # Create matrices to store simulated d_j and d_k
  sim_dj <- matrix(NA, nrow = length(idx), ncol = J_edge)
  sim_dk <- matrix(NA, nrow = length(idx), ncol = K_edge)
  
  # Loop over selected posterior samples
  for (s in seq_along(idx)) {
    
    if (model$model_type != "M5") {
      
      # Get alpha and beta from posterior sample
      alpha_s <- draws$alpha[idx[s]]
      beta_s <- draws$beta[idx[s]]
    }
    
    if (model$model_type == "M3") {
      
      # Simulate target-specific edge probabilities
      psi_j_s <- rbeta(J_edge, shape1 = alpha_s, shape2 = beta_s)
      D_s <- matrix(rbinom(K_edge * J_edge, size = 1, prob = rep(psi_j_s, each = K_edge)), nrow = K_edge, ncol = J_edge)
    }
    
    if (model$model_type == "M4") {
      
      # Simulate perturbation-specific edge probabilities
      psi_k_s <- rbeta(K_edge, shape1 = alpha_s, shape2 = beta_s)
      D_s <- matrix(rbinom(K_edge * J_edge, size = 1, prob = rep(psi_k_s, times = J_edge)), nrow = K_edge, ncol = J_edge)
    }
    
    if (model$model_type == "M5") {
      
      # Get one posterior draw of psi_j and psi_k
      psi_j_s <- draws$psi_j[idx[s], ]
      psi_k_s <- draws$psi_k[idx[s], ]
      
      # P(delta_kj = 1) = psi_k * psi_j
      prob_s <- outer(psi_k_s, psi_j_s, "*")
      
      # Simulate edge matrix
      D_s <- matrix(rbinom(K_edge * J_edge, size = 1, prob = as.vector(prob_s)), nrow = K_edge, ncol = J_edge)
    }
    
    # Compute simulated d_j
    sim_dj[s, ] <- colSums(D_s)
    
    # Compute simulated d_k
    sim_dk[s, ] <- rowSums(D_s)
  }
  
  # Plot observed vs posterior predictive d_j
  print(
    plot_ppc_degree(
      obs_degree = observed_dj,
      sim_degree = sim_dj,
      model_name = paste0(model$model_name, "(d_j): ", model$model_label),
      degree_type = "Target degree d_j",
      max_degree = K_edge
    )
  )
  
  # Plot observed vs posterior predictive d_k
  print(
    plot_ppc_degree(
      obs_degree = observed_dk,
      sim_degree = sim_dk,
      model_name = paste0(model$model_name, "(d_k): ", model$model_label),
      degree_type = "Perturbation degree d_k",
      max_degree = J_edge
    )
  )
}

dev.off()

log_msg("Saved ", pdf_file)
log_msg("03_plot_ppc_M1_to_M5.R finished")