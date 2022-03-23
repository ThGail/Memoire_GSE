print("Hello GSE")
print("Second update")


A_fct <- function(u, param, dj) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * dj * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_fct <- function(u, param, dj) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * dj * sigma**2)
  return((- 2 * dj * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

lambdat_fct2 <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda_t <- param[4]
  for (i in 1:t) {
    lambda_t <- lambda_t + (1 / N) * k * (mu - lambda_t) + sigma * sqrt(lambda_t * (lambda_t > 0))* rnorm(N, 0, 1)
  }
  return(lambda_t)
}

proba_defaut_i <- function(t, TT, param, M, D, i){
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda0 <- param[4]
  
  invM <- solve(M)
  K <- length(D[1,])
  
  sum = 0
  for (j in (1:(K-1))){
    esp <- A_fct(TT-t, param, D[j,j])*exp(-B_fct(TT-t, param, D[j,j])*pi(t))
    sum = sum + M[i,j]*invM[i,j]*(esp-1)
  }
  return(sum)
}