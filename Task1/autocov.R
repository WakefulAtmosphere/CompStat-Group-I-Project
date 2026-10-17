autocov_matrices <- function(Z, lags) {
  Z <- as.matrix(Z)
  n <- nrow(Z)

  if (any(lags < 1 | lags >= n)) {
    stop("Each lag must be between 1 and n - 1.")
  }
  
  Sigma_lags <- vector("list", length(lags))
  names(Sigma_lags) <- paste0("lag_", lags)
  
  for (i in seq_along(lags)) {
    tau <- lags[i]
    
    # Pair observations at t and t + tau
    Z_t <- Z[1:(n - tau), , drop = FALSE]
    Z_future <- Z[(tau + 1):n, , drop = FALSE]
    
    # Estimate E[z_t z_(t+tau)^T]
    # crossprod calculates Z_t^T Z_t+tau, which is a p x p matrix
    # Corresponds to the estimate presented in the paper for matrices with
    # observations as rows
    Sigma_lags[[i]] <- crossprod(Z_t, Z_future) /
      (n - tau)
  }
  
  Sigma_lags
}