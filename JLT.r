rm(list = ls())
set.seed(1)

library(readxl) 
library(dfoptim)

setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC    <- as.numeric(sheetDGlo$...2 [7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

LGD = 0.3

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
Ltemp=logMatrice(Q,100)
Lpos = pmax(Ltemp,0)
Lneg = pmin(Ltemp,0)
L = Lpos + diag(rowSums(Lneg))

# diagonalisation de la matrice L
M = eigen(L)$vectors
D = diag(eigen(L)$values)
# round(M%*%D%*%solve(M),15) permet de retrouver la matrice de base (modulo erreur numérique)

# formule de la proba de défaut
A_fct <- function(u, param, dj) {
  k <- param[,1]
  mu <- param[,2]
  sigma <- param[,3]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_fct <- function(u, param, dj) {
  k <- param[,1]
  mu <- param[,2]
  sigma <- param[,3]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  
  return((- 2 * dj * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# On suppose que le mouvement brownien est standard
pit_fct <- function(N, t, param) {
  k <- matrix(rep(param_test[,1], N), nrow = 8)
  mu <- matrix(rep(param_test[,2], N), nrow = 8)
  sigma <- matrix(rep(param_test[,3], N), nrow = 8)
  pi_t <- matrix(rep(param_test[,4], N), nrow = 8)
  
  if (t==0) {return(pi_t)}
  for (i in 1:t){
    pi_t <- abs(pi_t + k * (mu - pi_t) + t(t(sigma * sqrt(pi_t)) * rnorm(N, 0, 1)))
  }
  return(pi_t)
}

proba_defaut_i <- function(N, t, TT, param, M, D, i){
  # pour un rating donnée, on fait N simulations
  
  invM <- solve(M) # le mettre en param ?
  K <- length(D[1,])
  
  sum = 0
  for (j in (1:(K-1))){
    esp <- rowMeans(A_fct(TT-t, param, D[j,j])*exp(-B_fct(TT-t, param, D[j,j])*pit_fct(N,t,param)))[i] # enlever rowMeans pour avoir les trajectoires ?
    sum <- sum + M[i,j]*invM[j,K]*(esp-1)
  }
  return(sum)
}

proba_defaut_i_calibrage <- function(TT, param, M, D, i){
  # pour un rating donnée, et une maturité TT donnée, on a juste besion de t=0
  
  invM <- solve(M)
  K <- length(D[1,])
  pi0 <- param_test[,4]
  
  sum <- 0
  for (j in (1:(K-1))){
    esp <- A_fct(TT, param, D[j,j])[i]*exp(-B_fct(TT, param, D[j,j])[i]*pi0[i])
    sum <- sum + M[i,j]*invM[j,K]*(esp-1)
  }
  return(sum)
}

spread_i_fct <- function(N, t, TT, param, M, D, i){
  # if (t==TT){return(0)} # a verifier au temps t==TT
  # if (1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)<=0){return(0)} # à vérifier la condition
  return(-(1/(TT-t))*log(1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)))
}

spread_i_calibrage <- function(TT, param, M, D, i, LGD){
  return(-log(1 - LGD * proba_defaut_i_calibrage(TT, param, M, D, i))/TT)
}

# fonction objective à minimiser
# à minimiseer sur la maturité (ici 1 an ?) et sur les ratings
# un seul appel à cette fonction necessaire pour le calibrage 
# à calibrer K fois les paramètres k mu et pi0
# et 1 fois le paramètre sigma (car la volatilité est supposé constante)
Ecart_JLT <- function(param){
  e <- 0
  for (t in 1:TT){ for (i in 1:6){
    # spread <- spread_i_fct(1000, t, TT, param_test, M, D, i,LGD)
    spread <- spread_i_calibrage(t, param, M, D, i, LGD)
    e <- e + (spread-SpreadMarket[i])^2
  }}
  return(e)
} ## parametre sur le rating


######### CALIBRAGE#########
#On considere param comme une matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 8 classes de rating et 4 parametres pour chacune d entre elle,
#La dimension de la matrice param est 8 x 4
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, kCaa, kCa-C, muAAA, muAA, muA, muBaa, muBa, muB, muCaa, muCaC, sigma, pi0)

# param a garder temporairement
# les params du CIR
param_test <- matrix(c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003, 0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000, 8.727, 9.614, 9.367, 1.411, 0.001, 0.001, 0.001, 0.001,0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055), nrow =8)
# les params hasards
param_test <- matrix(c(2, 1.5, 2, 1.5, 2, 1.5, 2, 1.5, 0.3, 0.9, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.5, 0.6, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow =8)

N <- 1000
TT <- 5 # quelle maturité prendre ?
LB <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 8)
UB <- matrix(c(9, 9, 9, 9, 9, 9, 9, 9, 20, 20, 20, 20, 20, 20, 20, 20, 15, 15, 15, 15, 15, 15, 15, 15, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 8)
param_tt <- hjkb(param_test,Ecart_JLT,lower=LB,upper=UB)$par

# nlminb(start = param_test,Ecart_JLT,lower=LB,upper=UB)
# optim(par=param_test,fn=Ecart_JLT,lower=LB,upper=UB,method="L-BFGS-B")

spreadAAA <- c()
spreadAA <- c()
spreadA <- c()
spreadBBB <- c()
spreadBB <- c()
spreadB <- c()
for (t in 0:20) {
  spreadAAA <- c(spreadAAA, spread_i_calibrage(t, param_tt, M, D, 1,LGD))
  spreadAA <- c(spreadAA, spread_i_calibrage(t, param_tt, M, D, 2,LGD))
  spreadA <- c(spreadA, spread_i_calibrage(t, param_tt, M, D, 3,LGD))
  spreadBBB <- c(spreadBBB, spread_i_calibrage(t, param_tt, M, D, 4,LGD))
  spreadBB <- c(spreadBB, spread_i_calibrage(t, param_tt, M, D, 5,LGD))
  spreadB <- c(spreadB, spread_i_calibrage(t, param_tt, M, D, 6,LGD))
}
plot(spreadB, main='Spread', type='l', ylab="Spread", ylim=c(0,0.03), xlab="Maturité", col='purple')
lines(spreadAAA, col='red')
lines(spreadAA, col='orange')
lines(spreadA, col='brown')
lines(spreadBBB, col='lightblue')
lines(spreadBB, col='blue')


######### TEST #########
n <- 2
hist(pit_fct(100, 3, param_test), freq=FALSE, main='Histogramme des simulations (param random)', xlab='Simulation des lambda t = 3')
random <- rnorm(n,0,1)
print(random)

k_test <- matrix(rep(param_test[,1], n), nrow = 8)
mu_test <- matrix(rep(param_test[,2], n), nrow = 8)
sigma_test <- matrix(rep(param_test[,3], n), nrow = 8)
pi_t_test <- matrix(rep(param_test[,4], n), nrow = 8)

Ecart_JLT <- function(param){
  e <- 0
  for (t in 0:TT){
    spread <- spread_fct(1000, t, TT, param_test, M, D, i)
    e <- e + (spread-SpreadMarket)^2
  }
  return(e)
}

# Calibrage des paramètres
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



