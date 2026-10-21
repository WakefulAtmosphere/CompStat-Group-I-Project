stationary_bootstrap <- function(X, p) {
  n <- nrow(X)
  
  if (n < 1) {
    stop("X must contain at least one observation.")
  }
  
  if (p < 0 || p > 1) {
    stop("p must be between 0 and 1.")
  }
  
  # Store the original time indices selected
  index <- integer(n)
  
  # First observation: uniform random starting index
  index[1] <- sample.int(n, size = 1)
  
  # Generate the remaining bootstrap indices
  if (n > 1) {
    for (t in 2:n) {
      if (runif(1) < p) {
        # Start a new block
        index[t] <- sample.int(n, size = 1)
      } else {
        # Continue the current block, wrapping around
        index[t] <- (index[t - 1] %% n) + 1
      }
    }
  }
  # Resample entire rows
  X_boot <- X[index, , drop = FALSE]
  X_boot
}