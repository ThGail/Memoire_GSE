puissanceMatrice <- function(A, n){
  mat = diag(length(A[1,]))
  for (i in 1:n){mat = mat%*%A}
  return(mat)
}

logMatrice <- function(Q, n){
  QI = Q-diag(length(Q[1,])) # Q-Identite
  mat = matrix(rep(0,length(Q[1,])^2),nrow=length(Q[1,]))
  for (i in 1:n){
    mat = mat + (-1)**(i+1)*puissanceMatrice(QI,i)/i
  }
  return(mat)
}

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

# On suppose que le mouvement brownien est standard
pit_fct <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  pi0 <- param[4]
  pi_t <- pi0
  if (t==0) {return(pi_t)}
  for (i in 1:t) {
    pi_t <- abs(pi_t + k * (mu - pi_t) + sigma * sqrt(pi_t) * rnorm(N, 0, 1))
  }
  return(pi_t)
}

proba_defaut_i <- function(N, t, TT, param, M, D, i){
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda0 <- param[4]
  
  invM <- solve(M) # le mettre en param ?
  K <- length(D[1,])
  
  sum = 0
  for (j in (1:(K-1))){
    esp <- A_fct(TT-t, param, D[j,j])*exp(-B_fct(TT-t, param, D[j,j])*pit_fct(N,t,param))
    sum = sum + M[i,j]*invM[j,K]*(esp-1)
  }
  
  return(sum)
}