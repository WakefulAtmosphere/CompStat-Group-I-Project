phi <- function (lambda, k) {
  lambdas <- diag(lambda)
  lambdas <- sort(lambdas, decreasing = TRUE)
  
  lambdas[k]^2 / (1 + sum(lambdas^2))
}




f_0_k <- function(B_hat, B_star, k) {
  # function for one f_0(k)
  B_hat_k <- B_hat[, 1:k, drop = FALSE]
  B_star_k <- B_star[,1:k, ,drop = FALSE]
  # apply goes through bootstrap samples and calculates mean across them 
  f_0_k <- mean(apply(B_star_k, 3, function(B_star_j_k) {
    1- abs(det(t(B_hat_k)%*% B_star_j_k))
  }))
  return(f_0_k)
}



f_k <- function(B_hat, B_star, k) {
  # f_k calls f_0_k inside it
  # we calculate f_0's for all k = 1, ..., p-1
  # and combine them to a vector f_0
  # where p is the number of variables
  p = nrow(B_hat)
  f_0 <- numeric(p)   # not p-1 because the sum starts form zero
  for (k in seq(1, p-1)){
    f_0[k] <- f_0_k(B_hat, B_star_array, k)
  }
  # calculate the sum of f_0-vector 
  sum_f_0 <- sum(f_0)
  
  # final step
  f_0_k(B_hat, B_star_array, k) /(1 + sum_f_0)
}

g_k <- function(lambda, B_hat, B_star, k) {
  f_k(B_hat, B_star, k) + phi(lambda, k)
}