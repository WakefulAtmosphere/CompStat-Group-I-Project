generate_ladle_data <- function (n) {
  p <- 10
  q <- 5
  
  # Generate the five signal components
  Z_signal <- cbind(
    arima.sim(n = n, model = list(ar = 0.6, mr = c(0, -0.4))),
    arima.sim(n = n, model = list(ar = c(0.4, 0.1, 0.3), mr = c(0.2, 0.4))),
    arima.sim(n = n, model = list(ar = 0.4)),
    arima.sim(n = n, model = list(ma = c(0.5, 0.3, 0.1, 0.2))),
    arima.sim(n = n, model = list(ar = c(0.1, 0.1, 0.3)))
  )
  
  Z_noise <- matrix(rnorm(n * (p - q)), nrow = n, ncol = p - q)
  
  # Combine signal and noise and standardize each component
  Z <- scale(cbind(Z_signal, Z_noise))
  
  Omega <- diag(p)
  X <- Z %*% t(Omega)
  X
}