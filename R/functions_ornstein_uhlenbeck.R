set_ar_from_ou <- function(ou_parameter_list, delta_t = 1, integrated = TRUE){
  mu_ <- ou_parameter_list$mean
  gamma_ <- ou_parameter_list$recall
  gamma_inv <- solve(gamma_)
  sigma_ <- ou_parameter_list$diffusion
  dimension <- length(mu_)
  Id <- diag(1, dimension)
  # First, m and A
  # Dealing with velocity component
  A_vel <- expm::expm(-gamma_ * delta_t)
  m_vel <- as.numeric((Id - A_vel) %*% mu_)
  # Dealing with position components
  A_pos <- NULL
  m_pos <- NULL
  A_vel_pos <- NULL # How former position influence next vel
  A_pos_vel <- NULL # How former vel influence next pos
  if(integrated){
    A_pos <- Id
    recall_mat <- gamma_inv %*% (Id - A_vel)
    m_pos <- as.numeric((Id * delta_t - recall_mat) %*% mu_)
    A_vel_pos <- matrix(0, ncol = dimension, nrow = dimension)
    A_pos_vel <- recall_mat
  }
  m <- c(m_vel, m_pos)
  A <- rbind(cbind(A_vel, A_vel_pos),
             cbind(A_pos_vel, A_pos))
  # Then covariance matrix
  vec_S <- as.numeric(solve(kronecker(gamma_, Id) + kronecker(Id, gamma_)) %*% 
    as.numeric(sigma_ %*% t(sigma_)))
  S <- matrix(vec_S, nrow = dimension, ncol = dimension)
  covariance_vel <- S - A_vel %*% S %*% t(A_vel)
  covariance_pos <- NULL
  covariance_pos_vel <- NULL
  covariance_vel_pos <- NULL
  if(integrated){
    M <- S %*% t(gamma_inv) + gamma_inv %*% S
    covariance_pos_vel <- S %*% t(gamma_inv) %*% (Id - t(A_vel)) +
      gamma_inv %*% (A_vel - Id) %*% S %*% t(A_vel)
    covariance_vel_pos <- t(covariance_pos_vel)
    covariance_pos <- M * delta_t - 
      (Id - A_vel) %*% gamma_inv %*% M -
      M %*% t(gamma_inv) %*% (Id - t(A_vel)) +
      gamma_inv %*% S %*% t(gamma_inv) -
      A_vel %*% gamma_inv %*% S %*% t(gamma_inv) %*% t(A_vel)
  }
  Sigma <- rbind(cbind(covariance_vel, covariance_vel_pos), 
                 cbind(covariance_pos_vel, covariance_pos))
  list(m = m, A = A, Sigma = Sigma, delta_t = delta_t)
}
