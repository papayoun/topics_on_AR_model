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
