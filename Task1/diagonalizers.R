AMUSE_diagonalizer <- function(data, tau = 1) {
  lag_tau_sigma <- autocov_matrices(data, lags = tau)[[1]]
  U_tau <- eigen(lag_tau_sigma, symmetric = T)$vectors
  B_hat <- U_tau  # doesn't do anything, but I think it might be good for clarity for when we code the SOBI estimator B_hat
  return(B_hat)
}

SOBI_diagonalizer <- function(data, tau = 1){
  
}