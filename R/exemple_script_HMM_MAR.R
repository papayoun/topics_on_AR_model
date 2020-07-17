rm(list = ls())
library(tidyverse)
source("R/functions_ornstein_uhlenbeck.R")
source("R/functions_simulation.R")
source("R/functions_estimation.R")
set.seed(123)
true_ou_parameters <- list(list(mean = c(0, 0),
                           recall = matrix(c(0.31, -0.3, 0.3, 0.31),
                                           ncol = 2, byrow = T),
                           diffusion = diag(c(0.1, 0.8))),
                           list(mean = c(5, 5),
                                recall = matrix(c(0.5, 0, 0, 0.5),
                                                ncol = 2, byrow = T),
                                diffusion = diag(c(0.1, 0.8))))
true_ar_parameters <- purrr::map(true_ou_parameters,
                          set_ar_from_ou,
                          delta_t = 0.5)

transition_matrix <- matrix(c(0.90, 0.1,
                              0.1, 0.9), byrow = TRUE,
                            nrow = 2)
simulated_tibble <- simulate_HMM_MAR(500, Y0 = c(-10, -10, 0, 0),
                                     transition_matrix = transition_matrix,
                                     parameters = true_ar_parameters,
                                     delta_t = 0.5)
ggplot(simulated_tibble) +
  aes(Y3, Y4) +
  geom_path() +
  geom_point(aes(color = state))
simulated_tibble %>%
  tidyr::gather(-t, -state, key = "Metric", value = "Value") %>%
  ggplot(aes(x = t, y = Value)) +
  geom_path() +
  geom_point(aes(color = state)) +
  facet_wrap(~Metric, scales = "free_y", nrow = 4)

simulated_tibble %>%
  dplyr::select(-Y3,-Y4) %>% 
  mutate(weight = as.numeric(state == 1)) %>% 
  tidyr::gather(-t, -weight, -state, key = "Metric", value = "Value") %>%
  ggplot(aes(x = t, y = Value)) +
  geom_path() +
  geom_point(aes(color = state, alpha = weight)) +
  facet_wrap(~Metric, scales = "free_y", nrow = 4)

# Estimation --------------------------------------------------------------

obs_matrix <- simulated_tibble %>%
  dplyr::select(-t, -state) %>%
  as.matrix()
# weight <- as.numeric(simulated_tibble$state == 1)
# get_mle_estimates(obs_matrix, weight[-1]) %>% set_ou_from_ar(delta_t = 0.5)
# log_densities <- get_log_densites_HMM(obs_matrix, true_ar_parameters)
# alphas_log_lik <- compute_alphas(obs_matrix, true_ar_parameters, transition_matrix,
#                          log_densites = log_densities, return_likelihood = TRUE)
# alphas <- alphas_log_lik$alpha
# loglik <- alphas_log_lik$loglik
# betas <- compute_betas(obs_matrix, true_ar_parameters, transition_matrix,
#                        log_densites = log_densities)
# posterior_probabilities <- compute_gammas(obs_matrix,
#                                           true_ar_parameters, 
#                                           transition_matrix, 
#                                           log_densites = log_densities)
# posterior_marginal_prediction <- apply(posterior_probabilities, 1, which.max)
# 
# xis <- compute_xis(observations = obs_matrix, dynamics_param = true_ar_parameters,
#                    transition_mat = transition_matrix, alphas = alphas, betas = betas, 
#                    log_lik = loglik, 
#                    log_densites = log_densities, get_sum = T)
# 
# rm(list = ls())
start_ou_parameters <- list(list(mean = c(0, 0),
                                recall = matrix(c(1, 0, 0, 1),
                                                ncol = 2, byrow = T),
                                diffusion = diag(c(1, 1))),
                           list(mean = c(10, 10),
                                recall = matrix(c(0.5, 0, 0, 0.5),
                                                ncol = 2, byrow = T),
                                diffusion = diag(c(1, 0.8))))
start_ar_parameters <- purrr::map(start_ou_parameters,
                                 set_ar_from_ou,
                                 delta_t = 0.5)
start_transition_matrix <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)

estimation <- EM_algorithm(observations = obs_matrix, 
                           dynamics_param_0 = start_ar_parameters, 
                           transition_mat_0 = start_transition_matrix, n_iter = 10)
estimation$lls %>% 
  enframe(name = "Iteration") %>% 
  ggplot(aes(x = Iteration, y = value)) +
  geom_line() +
  geom_point()
