rm(list = ls())
source("R/functions_ornstein_uhlenbeck.R")
source("R/functions_simulation.R")
source("R/functions_estimation.R")

true_ou_parameters <- list(mean = c(0, 0), 
                        recall = matrix(c(2, -0.7, 0.7, 1),
                                        ncol = 2, byrow = T),
                        diffusion = diag(c(0.1, 0.8)))
true_ar_parameters <- set_ar_from_ou(ou_parameter_list =  true_ou_parameters,
                                     delta_t = 0.1)

simulated_matrix <- simulate_MAR(100, Y0 = c(-10, -10, 0, 0), 
                                parameters = true_ar_parameters)
get_mle_estimates(simulated_matrix)

simulated_matrix %>% 
  as_tibble() %>% 
  rename(X1 = V3, X2 = V4) %>% 
  ggplot(aes(X1, X2)) +
  geom_path()


