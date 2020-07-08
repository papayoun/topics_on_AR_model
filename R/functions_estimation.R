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
  Sigma_hat <- Reduce("+", lapply(1:n, function(i){
    residual <- observations[i + 1, ] - m_hat - A_hat %*% observations[i, ]
    weights[i] * residual %*% t(residual)
  })) / sum(weights)
  list(m_hat = m_hat, A_hat = A_hat, Sigma_hat = Sigma_hat)
}



