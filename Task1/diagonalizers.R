AMUSE_diagonalizer <- function(data, tau = 1) {
  lag_tau_sigma <- autocov_matrices(data, lags = tau)[[1]]
  
  # Symmetrize the lagged matrices
  lag_tau_sigma_sym <- lapply(Sigma_lags, function(Sigma) {
    (Sigma + t(Sigma)) / 2
  })
  
  U_tau <- eigen(lag_tau_sigma_sym, symmetric = T)$vectors
  B_hat <- U_tau  # doesn't do anything, but I think it might be good for clarity for when we code the SOBI estimator B_hat
  return(B_hat)
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
  JADE::rjd(Sigma_array)$V
}