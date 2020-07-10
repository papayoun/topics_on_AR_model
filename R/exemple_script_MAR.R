rm(list = ls())
library(tidyverse)
source("R/functions_ornstein_uhlenbeck.R")
source("R/functions_simulation.R")
source("R/functions_estimation.R")

true_ou_parameters <- list(mean = c(0, 0), 
                        recall = matrix(c(0.31, -0.3, 0.3, 0.31),
                                        ncol = 2, byrow = T),
                        diffusion = diag(c(0.1, 0.8)))
true_ar_parameters <- set_ar_from_ou(ou_parameter_list =  true_ou_parameters,
                                     delta_t = 0.5)

simulated_matrix <- simulate_MAR(1000, Y0 = c(-10, -10, 0, 0), 
                                parameters = true_ar_parameters)
get_mle_estimates(simulated_matrix)

simulated_tibble <- simulated_matrix %>% 
  as_tibble() %>% 
  rename(X1 = V3, X2 = V4) %>% 
  mutate(time = seq(0, length.out = nrow(simulated_matrix),
                    by = true_ar_parameters$delta_t))
ggplot(simulated_tibble) + 
  aes(X1, X2) +
  geom_path()
simulated_tibble %>% 
  tidyr::gather(-time, key = "Metric", value = "Value") %>% 
  ggplot(aes(x = time, y = Value)) +
  geom_path() +
  facet_wrap(~Metric, scales = "free_y")
