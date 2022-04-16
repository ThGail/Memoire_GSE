rm(list = ls())
set.seed(1)

library(readxl) 
library(dfoptim)

setwd("~/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC    <- as.numeric(sheetDGlo$...2 [7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

LGD = 0.3

# exemple de matrice de transition de rating 2017
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

expMatrice <- function(Q, n){
  mat = matrix(rep(0,length(Q[1,])^2),nrow=length(Q[1,]))
  for (i in 0:n){
    mat = mat + puissanceMatrice(Q,i)/factorial(i)
  }
  return(mat)
}

# transformation pour la positivité hors diag
Ltemp=logMatrice(Q,100)
Lpos = pmax(Ltemp,0)
Lneg = pmin(Ltemp,0)
L = Lpos + diag(rowSums(Lneg))

# expMatrice(L,100)
# 
# LL <- matrix(rep(0,9*9),nrow=9)
# for (i in 1:9){LL[i,i] <- log(Q[i,i])}
# for (i in 1:9){for (j in 1:9){if(j!=i){LL[i,j] <- Q[i,j]*log(Q[i,i])/(Q[i,i]-1)}}}

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
    pi_t <- abs(pi_t + k * (mu - pi_t) + t(t(sigma * sqrt(pi_t)) * rnorm(N, 0, 1)))}
  return(pi_t)
}

proba_defaut_i <- function(N, t, TT, param, M, D, i){
# pour un rating donnée, on fait N simulations
  
  invM <- solve(M)
  K <- length(D[1,])

  sum <- 0
  for (j in (1:(K-1))){
    esp <- (A_fct(TT-t, param, D[j,j])*exp(-B_fct(TT-t, param, D[j,j])*pit_fct(N,t,param)))[i,]
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

spread_i_fct <- function(N, t, TT, param, M, D, i, LGD){
  # if (t==TT){return(0)}
  # if (1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)<=0){return(0)}
  return(-(1/(TT-t))*log(1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)))
}

spread_i_calibrage <- function(TT, param, M, D, i, LGD){
  if (TT==0){return(0)}
  return(-log(1 - LGD * proba_defaut_i_calibrage(TT, param, M, D, i))/TT)
}


######### CALIBRAGE#########
#On considere param comme une matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 8 classes de rating et 4 parametres pour chacune d entre elle,
#La dimension de la matrice param est 8 x 4
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, kCaa, kCa-C, muAAA, muAA, muA, muBaa, muBa, muB, muCaa, muCaC, sigma, pi0)

#### TEST ( a garder temporairement) ####
# poser sigma constante !!!!
param_test <- matrix(c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003, 
                       0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000, 
                       8.727, 9.614, 9.367, 1.411, 0.001, 0.001, 0.001, 0.001,
                       0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055), nrow =8)
#param_test <- matrix(c(2, 1.5, 2, 1.5, 2, 1.5, 2, 1.5, 0.3, 0.9, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.5, 0.6, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow =8)
#param_test <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), nrow =8)


Ecart_JLT <- function(param){
  e <- 0
  for (t in 1:TT){ for (i in 1:6){
    # spread <- spread_i_fct(1000, t, TT, param_test, M, D, i,LGD)
    spread <- spread_i_calibrage(t, param, M, D, i, LGD)
    e <- e + (spread-SpreadMarket[i])^2
  }}
  return(e)
} ## parametre sur le rating

TT <- 2
LB <- matrix(rep(0,8*4), nrow = 8)
UB <- matrix(c(9, 9, 9, 9, 9, 9, 9, 9, 20, 20, 20, 20, 20, 20, 20, 20, 15, 15, 15, 15, 15, 15, 15, 15, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 8)
# UB <- matrix(rep(50,8*4), nrow = 8)
# UB <- matrix(c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), nrow = 8) 
param_tt <- hjkb(param_test,Ecart_JLT,lower=LB,upper=UB)$par # 2.189641e-05

# nlminb(start = param_test,Ecart_JLT,lower=LB,upper=UB)
# optim(par=param_test,fn=Ecart_JLT,lower=LB,upper=UB,method="L-BFGS-B")

spreadAAA <- c()
spreadAA <- c()
spreadA <- c()
spreadBBB <- c()
spreadBB <- c()
spreadB <- c()
for (t in 1:20) {
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

# simulation de spread pour une maturité de 5 ans 
tt <- seq(0,150)
Spread_i_fct <- Vectorize(spread_i_fct,"t")
SP = Spread_i_fct(500, tt, 5, param_tt, M, D, 2, LGD)
matplot(tt,t(SP),type="l",main="Simulation spread sur une maturité de 5ans pour AA",
        ylab="spread", xlab="temps t")
# on a parfois des valeurs négatives, 
# mais en faisant plot(colMeans(Spread_i_fct(500, tt, 5, param_tt, M, D, 2, LGD))) c'est ok

plot(tt,colMeans(Spread_i_fct(1000, tt, 5, param_tt, M, D, 3, LGD)),'l',
     main="Test de spread pour les actions AA de maturité 5 ans",
     xlab="Temps",
     ylab="Taux spread")



#### Test de martingalité ####  pb here
# prixZC_Vas <- function(t, TT, param, rt) {
#   a <- param[1]
#   b <- param[2]
#   sigma <- param[3]
#   return(exp(-b * (TT - t)) * exp(-(rt - b) * (1 - exp(-a * (TT - t))) / a + 0.5 * (sigma^2 * (TT - t) / a^2 - sigma^2 / a^3 * (1 - exp(-a * (TT - t))) - sigma^2 / (2 * a^3) * (1 - exp(-a * (TT - t)))^2)))
# }
# 
# PZC_JLT_act <- function(N,t,TT,rt,param_taux,param_JLT,M,D,i,LGD){
#   return(exp(-rt*t)*prixZC_Vas (t, TT, param_taux, rt)/(1-spread_i_fct(N, t, TT, param_JLT, M, D, i, LGD)^(TT-t)))
# }
# 
# tt <- seq(0,30)
# test <- c()
# for (t in tt){
#   test <- c(test,mean(PZC_JLT_act(N,0,t,TauxZC[t],param_free,param_tt,M,D,4,LGD)))
# }
# plot(test)
# param_free est parametre de Vasicek

# N=1000;i=1
# spread.act <- c()
# for (t in 0:5){spread.act <- c(spread.act, mean(spread_i_fct(N, t, TT, param_tt, M, D, i,LGD)))}
# plot(spread.act,ylim=c(0.9,1.1),pch=20,col="darkgrey",
#      main="Test de martingalité pour les actions",
#      xlab="Maturité",
#      ylab="Moyenne de l'indice actualisé ")
# abline(h=S0,col="red",lty=3,lwd=2)
# 
# abline(h=S0+0.025,col="blue",lty=4)
# abline(h=S0-0.025,col="blue",lty=4)
# legend("topleft",legend=c("Prix initial","Borne à 5%"),
#        col=c("red","blue"),pch=20,
#        cex=0.8)
