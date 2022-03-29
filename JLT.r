rm(list = ls())

# exemple de matrice de transition de rating
Q = matrix(c(91.06,8.25,0.6,0.07,0.02,0,0,0,0,
             0.86,89.58,8.96,0.45,0.07,0.04,0.02,0,0.02,
             0.05,2.58,91.02,5.63,0.51,0.11,0.04,0.01,0.05,
             0.03,0.15,4.34,90.41,3.99,0.72,0.16,0.02,0.18,
             0.01,0.04,0.46,6.66,83.17,7.81,0.77,0.12,0.96,
             0.01,0.03,0.16,0.5,5.35,82.3,7.41,0.58,3.66,
             0,0.01,0.03,0.1,0.4,7.6,79.24,3.33,9.29,
             0,0,0.06,0,0.72,2.94,11.48,50.57,34.23,
             0,0,0,0,0,0,0,0,100)/100,ncol=9,byrow=T)

# transformation en matrice log
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

# transformation pour la positivité hors diag
Ltemp=logMatrice(Q,5000)
Lpos = pmax(Ltemp,0)
Lneg = pmin(Ltemp,0)
L = Lpos + diag(rowSums(Lneg))

# diagonalisation de la matrice L
M = eigen(L)$vectors
D = diag(eigen(L)$values)
# round(M%*%D%*%solve(M),15) permet de retrouver la matrice de base (modulo erreur numérique)

# formule de la proba de défaut
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
pit_fct <- function(N, t, param)
{
  #k <- param[,1]
  #mu <- param[,2]
  #sigma <- param[,3]
  #pi_t <- param[,4]
  
  k <- matrix(rep(param_test[,1], N), nrow = 8)
  mu <- matrix(rep(param_test[,2], N), nrow = 8)
  sigma <- matrix(rep(param_test[,3], N), nrow = 8)
  pi_t <- matrix(rep(param_test[,4], N), nrow = 8)
  
  if (t==0) {return(pi_t)}
  for (i in 1:t)
  {
    #pi_t <- abs(pi_t + k * (mu - pi_t) + sigma * sqrt(pi_t) * rnorm(N, 0, 1))
    #pi_t <- abs(pi_t + k * (mu - pi_t) + (sigma * sqrt(pi_t)) %o% rnorm(N, 0, 1))
    pi_t <- abs(pi_t + k * (mu - pi_t) + t(t(sigma * sqrt(pi_t)) * rnorm(N, 0, 1)))
    print(i)
    print(pi_t)
  }
  return(pi_t)
}

calcul_esp <- function(pi_t)
{
  
}

proba_defaut_i <- function(N, t, TT, param, M, D, i){
  
  invM <- solve(M) # le mettre en param ?
  K <- length(D[1,])
  
  sum = 0
  for (j in (1:(K-1))){
    esp <- (1/N)*rowSums(A_fct(TT-t, param, D[j,j])*exp(-B_fct(TT-t, param, D[j,j])*pit_fct(N,t,param)))
    print('esp')
    print(esp)
    sum = sum + M[i,j]*invM[j,K]*(esp-1)
  }
  
  return(sum)
}
<<<<<<< HEAD

######### CALIBRAGE#########
#On considere param comme une matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 8 classes de rating et 4 parametres pour chacune d entre elle,
#La dimension de la matrice param est 8 x 4
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, kCaa, kCa-C, muAAA, muAA, muA, muBaa, muBa, muB, muCaa, muCaC, sigma, pi0)

#### TEST ( a garder temporairement) ####
param_test <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.3, 0.2, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.5, 0.5, 0.5, 0.7, 0.5, 0.5, 0.5, 0.8, 0.5, 0.5, 0.5), nrow =8)
print(param_test)

#temp
n <- 2
k_test <- matrix(rep(param_test[,1], n), nrow = 8)
mu_test <- matrix(rep(param_test[,2], n), nrow = 8)
sigma_test <- matrix(rep(param_test[,3], n), nrow = 8)
pi_t_test <- matrix(rep(param_test[,4], n), nrow = 8)

print(pi_t_test)

#pi_t <- abs(pi_t + k * (mu - pi_t) + sigma * sqrt(pi_t) * rnorm(N, 0, 1))
random <- rnorm(n,0,1)
print(random)
pi_t_test <- abs(pi_t_test + k_test * (mu_test - pi_t_test) + t(t(sigma_test * sqrt(pi_t_test)) * random))
random <- rnorm(n,0,1)
pi_t_test <- abs(pi_t_test + k_test * (mu_test - pi_t_test) + t(t(sigma_test * sqrt(pi_t_test)) * random))
random <- rnorm(n,0,1)
pi_t_test <- abs(pi_t_test + k_test * (mu_test - pi_t_test) + t(t(sigma_test * sqrt(pi_t_test)) * random))

print(temp)

k_test2 <- param_test[,1]
mu_test2 <- param_test[,2]
sigma_test2 <- param_test[,3]
pi_t_test2 <- param_test[,4]
print(abs(pi_t_test2) + k_test2 * (mu_test2 - pi_t_test2))
print(sigma_test2)
temp2 <- abs(pi_t_test2) + k_test2 * (mu_test2 - pi_t_test2) + (sigma_test2 * sqrt(pi_t_test2)) %o% random
print(temp2)








##### Calibrage des paramètres ####
hist(pit_fct(1000, 3, param_test), freq=FALSE, main='Histogramme des simulations (param random)', xlab='Simulation des lambda t = 3')

Ecart_JLT <- function(param){
  e <- 0
  for (t in 0:TT){
    spread <- Prix_Spread_fct(0,t,param,LGD)
    print(t)
    #print(spread)
    print()
    print(spread-spread_mkt)
    e <- e + (spread-spread_mkt)^2
    #print(e)
  }
  return(e)
}

test <- c(1,2,3,4)
print(test[1:3])
=======
>>>>>>> 608f574455448e132c160d839d5f0a371be34100
