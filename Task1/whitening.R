whiten_data <- function (X) {
  X <- as.matrix(X)
  
  # Center each variable
  X_c <- scale(X, center = TRUE, scale = FALSE)
  
  # Estimate sample covariance
  S <- cov(X_c)
  
  # Eigendecomposition of covariance matrix.
  # Used to calculate the inverse square root.
  eig <- eigen(S, symmetric = TRUE)
  
  # Check that whitening is numerically possible
  if (any(eig$values <= 0)) {
    stop("Covariance matrix must be positive definite.")
  }
  
  # Calculate the whitening matrix which corresponds to Sigma^(-1/2) in the paper
  W <- eig$vectors %*%
    diag(1 / sqrt(eig$values)) %*%
    t(eig$vectors)
  
  # Whiten the data. X_c has observations in rows.
  Z <- X_c %*% W
  Z
}