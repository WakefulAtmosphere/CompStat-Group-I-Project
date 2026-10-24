AMUSE_diagonalizer <- function(data, lags = 1) {
  lag_tau_sigma <- autocov_matrices(data, lags = lags)[[1]]
  
  # Symmetrize the lagged matrices
  lag_tau_sigma_sym <- (lag_tau_sigma + t(lag_tau_sigma)) / 2
  
  eig <- eigen(lag_tau_sigma_sym, symmetric = T)
  U_tau <- eig$vectors
  Lambda <- eig$values^2
  
  ord <- order(Lambda, decreasing = TRUE)
  Lambda <- Lambda[ord]
  U_tau <- U_tau[, ord, drop = FALSE]

  list(
    u_tau = U_tau,
    lambda = Lambda
  )
}

SOBI_diagonalizer <- function(data, lags){
  Sigma_lags <- autocov_matrices(data, lags)
  
  Sigma_sym <- lapply(Sigma_lags, function(Sigma) {
    (Sigma + t(Sigma)) / 2
  })
  
  # Store matrices in a p x p x K array
  p <- ncol(data)
  K <- length(Sigma_sym)
  Sigma_array <- array(0, dim = c(p, p, K))
  
  for (k in seq_len(K)) {
    Sigma_array[, , k] <- Sigma_sym[[k]]
  }
  
  # Joint diagonalization
  U_TAU <- JADE::rjd(
    Sigma_array,
    maxiter = 1000,
    eps = 1e-6
    )$V
  
  # Calculate the unmixing matrix
  Lambda <- numeric(p)
  
  for (k in seq_along(Sigma_sym)) {
    D_tau <- crossprod(
      U_TAU,
      Sigma_sym[[k]] %*% U_TAU
    )
    Lambda <- Lambda + diag(D_tau)^2
  }
  
  # Sort lambda and reorder the corresponding directions
  ord <- order(Lambda, decreasing = TRUE)
  Lambda <- Lambda[ord]
  U_TAU <- U_TAU[, ord, drop = FALSE]
  
  list(
    u_tau = U_TAU,
    lambda = Lambda
  )
}


diagonalizer_full <- function(data,
                         lags = 1,
                         n_boot = 200,
                         bootstrap_method = block_bootstrap,
                         block_size = 10,
                         diagonalizer = AMUSE_diagonalizer) {
  res <- diagonalizer(data, lags)
  
  # then we calculate B_star
  # B_star is an array
  B_star_array <- replicate(n_boot, {
    X_j <- bootstrap_method(data = data, block_size)
    B_hat <- diagonalizer(data = X_j, lags)$u_tau
  })
  
  list(
    B_hat = res$u_tau, 
    B_star = B_star_array,
    lambda = res$lambda
  )
}
