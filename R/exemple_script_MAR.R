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

foo <- function(...){
  simulated_matrix <- simulate_MAR(300, Y0 = c(-10, -10, 0, 0), 
                                   parameters = true_ar_parameters)
  mle <- get_mle_estimates(simulated_matrix)
  iou <- get_iou_estimates(simulated_matrix)
  mle$A <- mle$A[,c(1,2)]
  truth <- true_ar_parameters
  truth$A <- true_ar_parameters$A[,c(1,2)]
  res <- mapply(function(x, y, z, name_) {
    tibble(parameter = name_,
           mle = sum((x - z)^2),
           iou = sum((y - z)^2))
  },
  mle, iou, truth[-4],
  names(true_ar_parameters[-4]), SIMPLIFY = FALSE) %>% 
    bind_rows()
  res 
}
essais <- mclapply(1:(1e4), foo, mc.cores = 10) %>% 
  bind_rows()
essais %>% 
  gather(-parameter, key = "method", value = "erreur") %>%
  ggplot(aes(y = erreur, x = method)) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free_y")
essais %>% 
  mutate(diff = mle - iou) %>% 
  ggplot(aes(x = parameter, y = diff)) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free")
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
