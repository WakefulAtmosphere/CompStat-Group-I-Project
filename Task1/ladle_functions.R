phi <- function (lambda, k) {
  lambdas <- diag(lambda)
  lambdas <- sort(lambdas, decreasing = TRUE)
  
  lambdas[k]^2 / (1 + sum(lambdas^2))
}