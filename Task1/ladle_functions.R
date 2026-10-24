phi <- function(lambda, k) {
  lambdas <- sort(as.numeric(lambda), decreasing = TRUE)
  p <- length(lambdas)

  if (length(k) != 1 || is.na(k) ||
      k != as.integer(k) || k < 0 || k >= p) {
    stop("k must be an integer between 0 and p - 1.")
  }

  if (any(!is.finite(lambdas)) || any(lambdas < 0)) {
    stop("Squared eigenvalues must be finite and non-negative.")
  }
  # Mathematical k = 0 corresponds to the first eigenvalue
  lambdas[k + 1] / (1 + sum(lambdas))
}

f_0_k <- function(B_hat, B_star_array, k) {
 p <- nrow(B_hat)
 
 if (length(k) != 1 || is.na(k) ||
     k != as.integer(k) || k < 0 || k >= p) {
   stop("k must be an integer between 0 and p - 1.")
 }
 
 if (k == 0) return(0) # f_0(0) = 0
 
 n_boot <- dim(B_star_array)[3]
 B_hat_k <- B_hat[, 1:k, drop = FALSE]
 
 f0_values <- numeric(n_boot)
 for (b in seq_len(n_boot)) {
   # Keep the expected dimensions
   B_star_j_k <- matrix(
     B_star_array[, 1:k, b],
     nrow = p,
     ncol = k
   )
   f0_values[b] <- 1 - abs(
     det(t(B_hat_k) %*% B_star_j_k)
   )
 }
 mean(f0_values)
}

f_k <- function(B_hat, B_star_array, k) {
  # f_k calls f_0_k inside it
  # we calculate f_0's for all k = 1, ..., p-1
  # and combine them to a vector f_0
  # where p is the number of variables
  p = nrow(B_hat)
  
  if (length(k) != 1 || is.na(k) ||
      k != as.integer(k) || k < 0 || k >= p) {
    stop("k must be an integer between 0 and p - 1.")
  }
  
  if (k == 0) return(0) # f_0 (0) = 0 so f_k is also 0
  
  f_0 <- numeric(p)   # not p-1 because the sum starts form zero
  f_0[1] <- 0 # R index 1 corresponds to mathematical index 0
  for (j in seq_len(p - 1)) {
    f_0[j + 1] <- f_0_k(B_hat, B_star_array, j)
  }
  # R index k+1 corresponds to mathematical index k
  f_0[k + 1] / (1 + sum(f_0))
}

g_k <- function(lambda, B_hat, B_star_array, k) {
  p = nrow(B_hat)

  if (length(k) != 1 || is.na(k) ||
      k != as.integer(k) || k < 0 || k >= p) {
    stop("k must be an integer between 0 and p - 1.")
  }

  f_k(B_hat, B_star_array, k) + phi(lambda, k)
}
