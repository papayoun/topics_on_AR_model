get_mle_estimates <- function(observations, weights = NULL){
  n <- nrow(observations) - 1
  if(is.null(weights)){
    weights <- rep(1, n)
  }
  W <- diag(weights)
  if(nrow(W) != n){
    stop("weights must be a vector of size nrow(observations) - 1")
  }
  Y <- observations[-1, ]
  X <- cbind(rep(1, n), observations[-(n + 1), ])
  B_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  m_hat <- B_hat[1, ]
  A_hat <- t(B_hat[-1, ])
  pred_mat <- Y - X %*% B_hat 
  Sigma_hat <- t(pred_mat) %*% diag(weights / sum(weights)) %*% pred_mat
  list(m = m_hat, A = A_hat, Sigma = Sigma_hat)
}

get_log_densites_HMM <- function(observations, parameters){
  n <- nrow(observations) - 1
  sapply(parameters, 
         function(param){ 
           c(0, # Initial distribution, arbitrary in our case
             parallel::mclapply(1:n, # From library parallel
                                function(t){
                                  mean = as.numeric(param$m + 
                                                      param$A %*% observations[t, ])
                                  var = param$Sigma
                                  mixtools::logdmvnorm(observations[t + 1, ], 
                                                       mean, var)
                                },
                                mc.cores = parallel::detectCores() - 1) %>% 
               unlist())
         }
  )
}

compute_alphas <- function(observations, 
                           dynamics_param,
                           transition_mat,
                           initial_distribution = NULL,
                           return_likelihood = F,
                           log_densites = NULL){
  observations <- as.matrix(observations)
  n <- nrow(observations) - 1
  K <- length(dynamics_param)
  if(is.null(initial_distribution)){
    initial_distribution <- rep(1 / K, K)
  }
  if(is.null(log_densites)){
    log_densites <- get_log_densites_HMM(observations, dynamics_param)
  }
  mat_log_alpha <- matrix(nrow = n + 1, ncol = K)
  mat_log_alpha[1, ] <- log(initial_distribution) + log_densites[1, ]
  for(t in 1:n){
    max_log_alpha_t <- max(mat_log_alpha[t, ])
    old_alpha <- exp(mat_log_alpha[t, ] - max_log_alpha_t)
    mat_log_alpha[t + 1, ] <- sapply(1:K, function(k){
      somme <- max_log_alpha_t + log(sum(transition_mat[, k] * old_alpha))
      return(log_densites[t + 1, k] + somme)
    }) 
  }
  max_log_alpha_n <- max(mat_log_alpha[n, ])
  likelihood <- max_log_alpha_n + log(sum(exp(mat_log_alpha[n, ] - max_log_alpha_n)))
  if(return_likelihood){
    return(list(alpha = mat_log_alpha, loglik = likelihood))
  }
  return(mat_log_alpha)
}

compute_betas <- function(observations, dynamics_param, 
                          transition_mat,
                          initial_distribution = NULL,
                          log_densites = NULL){
  n <- nrow(observations) - 1
  K <- ncol(transition_mat)
  if(is.null(log_densites)){
    log_densites <- get_log_densites_HMM(observations, dynamics_param)
  }
  mat_beta <- matrix(nrow = n + 1, ncol = ncol(transition_mat))
  mat_beta[n + 1, ] <- 0 # log scale
  for(t in n:1){
    max_beta <- max(mat_beta[t + 1, ])
    p_obs <- log_densites[t + 1, ]
    mat_beta[t, ] <- sapply(1:K, function(k){
      max_beta + log(sum(transition_mat[k, ] * exp(p_obs) * 
                           exp(mat_beta[t + 1, ] - max_beta)))
    })
  }
  return(mat_beta)
}


compute_gammas <- function(observations, 
                           dynamics_param, 
                           transition_mat,
                           alphas = NULL, betas = NULL, log_lik = NULL,
                           log_densites = NULL){
  if(is.null(alphas) | is.null(log_lik)){
    new_alphas <- compute_alphas(observations, dynamics_param, 
                                 transition_mat, log_densites = log_densites, 
                                 return_likelihood = T)
    alphas <- new_alphas$alpha
    log_lik <- new_alphas$loglik
  }
  if(is.null(betas)){
    betas <- compute_betas(observations, dynamics_param, 
                           transition_mat, log_densites = log_densites)
  }
  res_nn <- exp(betas + alphas - log_lik)
  res <- t(apply(res_nn, 1, function(x) x / sum(x)))
  res
}

compute_xis <- function(observations, dynamics_param, 
                        transition_mat,
                        alphas, betas, log_lik, log_densites = NULL,
                        get_sum = F){
  observations <- as.matrix(observations)
  n <- nrow(observations) - 1
  K <- ncol(transition_mat)
  if(is.null(log_densites)){
    log_densites <- get_log_densites_HMM(observations, dynamics_param)
  }
  res <- array(dim = c(K, K, n + 1))
  for(t in 1:n){
    foo <- function(j, k){ # z_{t - 1} = j, z_t = k
      exp(alphas[t, j] + betas[t + 1, k] + log(transition_mat[j, k]) + 
            log_densites[t + 1, k] - log_lik)
    }
    res[ , ,t] <- sapply(1:K, function(j) sapply(1:K, 
                                                 function(k) foo(j, k)))
    res[,,t] <- res[,,t] / sum(res[,,t])
  }
  if(get_sum){
    return(Reduce("+", lapply(1:n, function(t) res[,,t])))
  }
  res
} 

viterbi_algorithm <- function(observations, dynamics_param, transition_mat,
                              initial_distribution = NULL,
                              log_densites = NULL){
  observations <- as.matrix(observations)
  n <- nrow(observations) - 1
  K <- ncol(transition_mat)
  if(is.null(initial_distribution)){
    initial_distribution <- rep(1 / K, K) 
  }
  if(is.null(log_densites)){
    log_densites <- get_log_densites_HMM(observations, dynamics_param)
  }
  mat_delta <- matrix(nrow = n + 1, ncol = K)
  mat_delta[1, ] <- log(initial_distribution) + log_densites[1, ]
  mat_argmax <- matrix(NA, nrow = n, ncol = K)
  for(t in 1:n){
    max_argmax <- sapply(1:K, function(k){
      p_obs <- log_densites[t + 1, k]
      back_kern <- log(transition_mat[, k]) + mat_delta[t, ]
      c(max = max(back_kern) + p_obs, argmax =  which.max(back_kern))
    })
    mat_delta[t + 1, ] <- max_argmax["max", ]
    mat_argmax[t, ] <- max_argmax["argmax",] 
  }
  viter_sequence <- rep(NA, n + 1)
  viter_sequence[n + 1] <- which.max(mat_delta[n + 1, ])
  for(t in n:1){
    viter_sequence[t] <- mat_argmax[t, viter_sequence[t + 1]]
  }
  return(viter_sequence)
}

Q_theta <- function(observations, dynamics_param, transition_mat, log_densities, gammas, xis,
                    initial_distribution = NULL){
  K <- ncol(transition_mat)
  n <- nrow(observations) - 1
  if(is.null(initial_distribution)){
    initial_distribution <- rep(1 / K, K) 
  }
  term1 <- sum(gammas[1, ] * (log(initial_distribution) + log_densities[1, ]))
  term2 <- sapply(1:n, function(t){
    sum(t(xis[,,t]) * log(transition_mat)) 
  }) %>% sum()
  term3 <- sum(gammas[-1, ] * log_densities[-1, ])
  c(init = term1, trans = term2, obs = term3, tot = term1 + term2 + term3)
}


EM_algorithm <- function(observations, dynamics_param_0, transition_mat_0,
                         n_iter = 5){
  likelihoods <- rep(NA, n_iter)
  n <- nrow(observations) - 1
  K <- ncol(transition_mat_0)
  for(iter in 1:n_iter){
    print(iter)
    log_densites <- get_log_densites_HMM(observations, dynamics_param_0)
    alphas <- compute_alphas(observations = observations, dynamics_param = dynamics_param_0, 
                             transition_mat = transition_mat_0,
                             return_likelihood = T,
                             log_densites = log_densites)
    likelihoods[iter] <- alphas$loglik
    alphas_cour <- alphas$alpha
    betas_cour <- compute_betas(observations = observations, 
                                dynamics_param = dynamics_param_0, 
                                transition_mat = transition_mat_0, 
                                log_densites = log_densites)
    gammas_cour <- compute_gammas(observations = observations, 
                                  dynamics_param = dynamics_param_0, 
                                  transition_mat = transition_mat_0,  
                                  alphas = alphas_cour, 
                                  betas = betas_cour, 
                                  log_lik = likelihoods[iter],
                                  log_densites = log_densites)
    xis_cour <- compute_xis(observations = observations, 
                            dynamics_param = dynamics_param_0, 
                            transition_mat = transition_mat_0,
                            alphas = alphas_cour, 
                            betas = betas_cour, 
                            log_lik = likelihoods[iter], 
                            log_densites = log_densites,
                            get_sum = F)
    dynamics_param_0 <- lapply(1:K,
                               function(k)
                                 get_mle_estimates(observations, 
                                                   weights = gammas_cour[-1, k]))
    xis_cour_sum <- Reduce("+", lapply(1:n, function(t) xis_cour[,,t]))
    transition_mat_0 <- (t(xis_cour_sum) / colSums(gammas_cour[-(n + 1), ])) %>%
    {. / rowSums(.)}
  }
  likelihoods <- c(likelihoods, compute_alphas(observations = observations, 
                                               dynamics_param = dynamics_param_0,
                                               transition_mat = transition_mat_0,
                                               return_likelihood = T)$loglik)
  est_sequence <- viterbi_algorithm(observations = observations, 
                                    dynamics_param = dynamics_param_0, 
                                    transition_mat = transition_mat_0)
  return(list(estimate = dynamics_param_0, lls  = likelihoods, states = est_sequence))
}


