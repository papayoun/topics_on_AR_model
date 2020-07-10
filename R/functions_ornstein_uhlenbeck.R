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

set_ou_from_ar <- function(ar_parameter_list, delta_t = NULL, integrated = TRUE){
  m <- ar_parameter_list$m
  A <- ar_parameter_list$A
  Sigma <- ar_parameter_list$Sigma
  delta_t <- ifelse(is.null(delta_t), 1,
                    delta_t)
  if(integrated){
    if((ncol(A) %% 2) != 0){
      stop("In the integrated case, the number of columns of A must be even")
    }
    dimension <- ncol(A) * 0.5
  }
  Id <- diag(1, dimension)
  A_vel <- A[1:dimension, 1:dimension]
  gamma_ <- -expm::logm(A_vel) / delta_t 
  m_vel <- m[1:dimension]
  mu_ <- as.numeric(solve(Id - A_vel) %*% m_vel)
  cov_vel_vec <- Sigma[1:dimension, 1:dimension] %>% as.numeric()
  S_vec <- solve(diag(1, 2 * dimension) - kronecker(A_vel, A_vel)) %*% cov_vel_vec
  diffusion_ <- ((kronecker(gamma_, Id) + kronecker(Id, gamma_)) %*% S_vec) %>% 
    matrix(nrow = dimension, ncol = dimension) %>% chol() %>% round(12)
  list(mean = round(mu_, 15), recall = round(gamma_, 15), 
       diffusion = round(diffusion_, 15))
}
