AMUSE_diagonalizer <- function(data, lags = 1) {
  lag_tau_sigma <- autocov_matrices(data, lags = lags)[[1]]
  
  # Symmetrize the lagged matrices
  lag_tau_sigma_sym <- (lag_tau_sigma + t(lag_tau_sigma)) / 2
  
  U_tau <- eigen(lag_tau_sigma_sym, symmetric = T)$vectors
  
  # Calculate the AMUSE unmixing matrix
  Lambda <- diag(t(U_tau) %*% lag_tau_sigma_sym %*% U_tau)
  Lambda <- Lambda^2

  list(
    u_tau = U_tau,
    lambda = Lambda
  )
}

SOBI_diagonalizer <- function(data, lags){
  Sigma_lags <- autocov_matrices(data, lags)
  
  # Store matrices in a p x p x K array
  p <- ncol(Z)
  K <- length(Sigma_lags)
  Sigma_array <- array(0, dim = c(p, p, K))
  
  for (k in seq_len(K)) {
    Sigma_array[, , k] <- Sigma_lags[[k]]
  }
  
  # Joint diagonalization
  U_TAU <- JADE::rjd(Sigma_array)$V
  
  # Calculate the unmixing matrix
  Lambda <- numeric(p)
  
  for (k in seq_along(Sigma_lags)) {
    Sigma_tau <- Sigma_lags[[k]]
    Lambda <- Lambda + 
      diag(t(U_TAU) %*% Sigma_tau %*% U_TAU)^2
  }

  list(
    u_tau = U_TAU,
    lambda = Lambda
  )
}


diagonalizer <- function(data,
                         lags = 1,
                         n_boot = 200,
                         bootstrap_method = block_bootstrap,
                         diagonalizer = AMUSE_diagonalizer) {
  res <- diagonalizer(data, lags)
  
  # then we calculate B_star 
  # B_star is an array
  B_star_array <- replicate(n_boot, {
    X_j <- bootstrap_method(data = data, B = 10)
    B_hat <- diagonalizer(data = data, lags)$u_tau
  })
  
  list(
    B_hat = res$u_tau, 
    B_star = B_star_array,
    lambda = Lambda
  )
}








