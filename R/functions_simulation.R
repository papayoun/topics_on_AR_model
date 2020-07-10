# Simple multivariate autoregressive --------------------------------------

simulate_MAR <- function(n, Y0, parameters){
  dimension <- length(Y0)
  # Initialization
  output <- matrix(NA, ncol =  dimension, nrow = n + 1)
  output[1, ] <- Y0
  m <- parameters$m
  A <- parameters$A
  Sigma <- parameters$Sigma
  m_innov <- rep(0, dimension)
  for(t in 1:n){
    output[t + 1, ] <- m + A %*% output[t, ] + 
      mixtools::rmvnorm(1, m_innov, Sigma)[1, ]
  }
  return(output)
}


# HMM ---------------------------------------------------------------------

simulate_markov_chain <- function(n, transition_matrix, first_state = NA){
  K <- ncol(transition_matrix)
  if(is.na(first_state)){# On fixe le premier état à 1 si non spécifié
    first_state <- 1
  }
  state_sequence <- rep(NA, n + 1)# Initialisation du vecteur de sorties
  state_sequence[1] <- first_state
  for(k in 1:n){
    state_sequence[k + 1] <- sample(K, 
                                prob = transition_matrix[state_sequence[k], ],
                                size = 1)
    # On tire dans une loi dont les poids sont données par P
  }
  return(state_sequence)
}

simulate_HMM_MAR <- function(n, parameters, transition_matrix, Y0, 
                             first_state = NA, delta_t = 1){
  dimension <- length(Y0)
  m_innov <- rep(0, dimension)
  m_intercepts <- purrr::map(parameters, "m") # K-elements list
  A_matrices <- purrr::map(parameters, "A") # K-elements list
  Sigma_matrices <- purrr::map(parameters, "Sigma") # K-elements list
  state_sequence <- simulate_markov_chain(n, transition_matrix, first_state)
  output <- matrix(NA, ncol =  dimension, nrow = n + 1, dimnames = 
                     list(row = NULL, col = paste0("Y", 1:dimension)))
  output[1, ] <- Y0
  for(t in 1:n){
    k <- state_sequence[t + 1]
    output[t + 1, ] <- m_intercepts[[k]] + 
      A_matrices[[k]] %*% output[t, ] + 
      mixtools::rmvnorm(1, m_innov, Sigma_matrices[[k]])[1, ]
  }
  output <- mutate(as_tibble(output),
                   state = factor(state_sequence),
                   t = seq(from = 0, length.out = n + 1, 
                           by = delta_t))
  return(output)
}
