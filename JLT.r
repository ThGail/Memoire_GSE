rm(list = ls())

library(readxl) 
library(dfoptim)

setwd("~/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
Maturite30 <- Maturite[1:30]
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
  k = param[1:8]
  mu = param[9:16]
  sigma = param[17]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_fct <- function(u, param, dj) {
  k = param[1:8]
  mu = param[9:16]
  sigma = param[17]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  return((- 2 * dj * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# On suppose que le mouvement brownien est standard
pit_fct <- function(N, t, param) {
  k <- matrix(rep(param[1:8], N), nrow = 8)
  mu <- matrix(rep(param[(1:8)+8], N), nrow = 8)
  sigma <- matrix(rep(param[17], N*8), nrow = 8)
  pi_t <- matrix(rep(param[18:25], N), nrow = 8)
  
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
proba_defaut_i.T <- Vectorize(proba_defaut_i,"TT")

proba_defaut_i_calibrage <- function(TT, param, M, D, i){
# pour un rating donnée, et une maturité TT donnée, on a juste besion de t=0
  
  invM <- solve(M)
  K <- length(D[1,])
  pi0 = param[18:25]
  
  sum <- 0
  for (j in (1:(K-1))){
    esp <- A_fct(TT, param, D[j,j])[i]*exp(-B_fct(TT, param, D[j,j])[i]*pi0[i])
    sum <- sum + M[i,j]*invM[j,K]*(esp-1)
  }
  return(sum)
}
proba_defaut_i_calibrage.T <- Vectorize(proba_defaut_i_calibrage,"TT")

spread_i_fct <- function(N, t, TT, param, M, D, i, LGD){
  if (t==TT){return(rep(0,N))}
  # if (1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)<=0){return(0)}
  return(-(1/(TT-t))*log(1- LGD * proba_defaut_i(N, t, TT, param, M, D, i)))
}
spread_i_fct.t <- Vectorize(spread_i_fct,"t")
spread_i_fct.T <- Vectorize(spread_i_fct,"TT")

spread_i_calibrage <- function(TT, param, M, D, i, LGD){
  if (TT==0){return(0)}
  return(-log(1 - LGD * proba_defaut_i_calibrage(TT, param, M, D, i))/TT)
}
spread_i_calibrage.T <- Vectorize(spread_i_calibrage,"TT")





########## CALIBRAGE ##########

#On considere param comme une matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 8 classes de rating et 4 parametres pour chacune d entre elle,
#La dimension de la matrice param est 8 x 4
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, kCaa, kCa-C, muAAA, muAA, muA, muBaa, muBa, muB, muCaa, muCaC, sigma, pi0)

#### TEST ( a garder temporairement)
# poser sigma constante !!!!
# param_init <- matrix(c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003,
#                        0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000,
#                        8.727, 9.614, 9.367, 1.411, 0.001, 0.001, 0.001, 0.001,
#                        0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055), nrow =8)
param_init <- c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003,
                0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000,
                0.01,
                0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055)

Ecart_JLT <- function(param){
  e <- 0
  for (t in 1:TT){ for (i in 1:6){
    # spread <- spread_i_fct(1000, t, TT, param_test, M, D, i,LGD)
    spread <- spread_i_calibrage(t, param, M, D, i, LGD)
    e <- e + (spread-SpreadMarket[i])^2
  }}
  return(e)
} ## parametre sur le rating

TT <- 1 # les données fornis sont les spreads de maturité 1 an
LB <- c(rep(0,8),rep(0,8),0,rep(0,8))
UB <- c(rep(30,8),rep(20,8),10,rep(30,8))
#paramJLT = hjkb(param_init,Ecart_JLT,lower=LB,upper=UB)$par # 2.731286e-05
paramJLT = nlminb(start = param_init,Ecart_JLT,lower=LB,upper=UB)$par
#paramJLT = optim(par=param_init,fn=Ecart_JLT,lower=LB,upper=UB,method="L-BFGS-B")$par
(listJLT <- list(
  k = paramJLT[1:8],
  mu = paramJLT[9:16],
  sigma = paramJLT[17],
  pi0 = paramJLT[18:25]))
(paramJLT[17] = min(sqrt(2*listJLT$k*listJLT$mu)[sqrt(2*listJLT$k*listJLT$mu) > 0]))
# il faut que ce parametre soit non nul



########## Some plot ##########

# Plot spread
spreadAAA <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 1, LGD)
spreadAA <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 2, LGD)
spreadA <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 3, LGD)
spreadBBB <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 4, LGD)
spreadBB <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 5, LGD)
spreadB <- spread_i_calibrage.T(Maturite30, paramJLT, M, D, 6, LGD)

plot(spreadAAA, main='Spread', type='l', ylab="Spread", ylim=c(0,0.04), xlab="Maturité", col='red')
lines(spreadAA, col='orange')
lines(spreadA, col='brown')
lines(spreadBBB, col='lightblue')
lines(spreadBB, col='blue')
lines(spreadB, col='purple')

#Essai avec proba de défaut calibrage :
probaDefAAA <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 1)
probaDefAA <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 2)
probaDefA <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 3)
probaDefBBB <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 4)
probaDefBB <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 5)
probaDefB <- proba_defaut_i_calibrage.T(Maturite30, paramJLT, M, D, 6)

plot(probaDefAAA, main='Probabilité de défaut calibrage', ylim=c(0,1), type='l', ylab="Spread", xlab="Maturité", col='red')
lines(probaDefAA, col='orange')
lines(probaDefA, col='brown')
lines(probaDefBBB, col='lightblue')
lines(probaDefBB, col='blue')
lines(probaDefB, col='purple')

#Essai avec proba de défaut normal :
probaDefAAAnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 1))
probaDefAAnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 2))
probaDefAnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 3))
probaDefBBBnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 4))
probaDefBBnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 5))
probaDefBnorm <- colMeans(proba_defaut_i.T(1000, 0, Maturite30, paramJLT, M, D, 6))

plot(probaDefAAAnorm, main='Probabilité de défaut calibrage', ylim=c(0,1), type='l', ylab="Spread", xlab="Maturité", col='red')
lines(probaDefAAnorm, col='orange')
lines(probaDefAnorm, col='brown')
lines(probaDefBBBnorm, col='lightblue')
lines(probaDefBBnorm, col='blue')
lines(probaDefBnorm, col='purple')

#simulation de spread AA pour une maturité de 10 ans
tt <- seq(0,10,0.1)[-length(seq(0,10,0.1))]
SP <- spread_i_fct.t(10, tt, 10, paramJLT, M, D, 2, LGD)
matplot(tt,t(SP),type="l",main="Simulation spread sur une maturité de 10 ans pour AA",
        ylab="spread", xlab="temps t")
plot(tt,colMeans(spread_i_fct.t(1000, tt, 10, paramJLT, M, D, 2, LGD)),'l',
     main="Test de spread pour les actions AA de maturité 10 ans",
     xlab="Temps",
     ylab="Taux spread")





########## MARTINGALITE ##########

#### Les fonctions utiles du Vasicek
# paramVas <- c(0.0095787241,0.0382841814,0.0002632241)
paramVas <- c(0.004084473,0.292354584,0.003962402) # calibrer à partir de TZC_Vas_FF_cali

taux_Vas_sim <- function(N, t, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(rnorm(N, r0 * exp(-a * t) + b * (1 - exp(-a * t)), sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a))))
}
taux_Vas_sim.t <- Vectorize(taux_Vas_sim,"t")

TZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  return(Ri - ((Ri-r0)*(1-exp(-a*TT))-sigma^2/(4*a^2)*(1-exp(-a*TT))^2)/(a*TT))
}

PZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b*TT) * exp(-(r0-b)*(1-exp(-a*TT))/a + 0.5 * (sigma^2 * TT / a^2 - sigma^2 / a^3 * (1 - exp(-a * TT)) - sigma^2 / (2 * a^3) * (1 - exp(-a * TT))^2)))
}
PZC_Vas_FF_calibrage.T <- Vectorize(PZC_Vas_FF_calibrage,"TT")

# resultat calibration sur le TZC
plot(TauxZC)
lines(TZC_Vas_FF_calibrage(Maturite,paramVas))

# resultat sur PZC
plot(exp(-Maturite*TauxZC))
lines(exp(-Maturite*TZC_Vas_FF_calibrage(Maturite,paramVas)))

# test pour voir erreur de calibrage du Vas
plot(PZC_Vas_FF_calibrage(Maturite,paramVas)/exp(-Maturite*TauxZC))


########## Test de martingalité : CASH FLOW ##########

# prix risqué à t=0, calculé à partir du spread FF et PZC FF
PZCr_i_CF_JLT_FF <- function(TT, param_taux, param_JLT, M, D, i, LGD){
  return(PZC_Vas_FF_calibrage.T(TT, param_taux)/(1+spread_i_calibrage.T(TT, param_JLT, M, D, i, LGD))^TT)
}
plot(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 1, LGD))

# moyenne prix ZC risqué actualisé au taux sans risque
PZCr_i_CF_JLT_sim <- function(N, TT, param_taux, param_JLT, M, D, i, LGD){
  pd_i <- proba_defaut_i.T(N, 0, TT, param_JLT, M, D, i)
  #return(exp(-TT*TauxZC[TT])*t(1*(1-pd_i)+(1-LGD)*pd_i))
  return(exp(-TT*TZC_Vas_FF_calibrage(TT, param_taux))*t(1*(1-pd_i)+(1-LGD)*pd_i))
}
plot(Maturite, rowMeans(PZCr_i_CF_JLT_sim(N=100, Maturite, paramVas, paramJLT, M, D, 1, LGD)))

# les plots
plot(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 1, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 1, LGD)),
     ylim=c(0.985,1.015),"l",col="red",ylab="CashFlow Actualisé",main="test de martingalité sur CASHFLOW")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 2, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 3, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 4, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 5, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 6, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, Maturite, paramVas, paramJLT, M, D, 6, LGD)),col="purple")

########## Test de martingalité : PZCr avec rating ##########

# PZCr_i_JLT blabla reste à coder

