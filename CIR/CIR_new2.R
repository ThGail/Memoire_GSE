rm(list = ls())

library(readxl) 
library(dfoptim)

setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC    <- as.numeric(sheetDGlo$...2 [7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

N <- 1000

A_fct <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_fct <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return((2 * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# On suppose que le mouvement brownien est standard
lambdat_fct <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda_t <- param[4]
  if (t==0) {return(rep(lambda_t, N))}
  for (i in 1:t) {
    lambda_t <- abs(lambda_t + k * (mu - lambda_t) + sigma * sqrt(lambda_t) * rnorm(N, 0, 1))
  }
  return(lambda_t)
}

#Bien une Chi-2 ? 
param_test <- c(0.5, 0.5, 0.5, 0.5)
hist(lambdat_fct(1000, 3, param_test), freq=FALSE, main='Histogramme des simulations (param random)', xlab='Simulation des lambda t = 3')

survie_FF_fct <- function(N, t, TT, param) {
  # On a calculé la moyenne en sortie de fonction et non directement dans la fonction
  lambda_t <- mean(lambdat_fct(N, t, param))
  B <- B_fct(TT - t, param)
  A <- A_fct(TT - t, param)
  return(A * exp(-B * lambda_t))
}

survie_FF_fct_sim <- function(N, t, TT, param) {
  # On a calculé la moyenne en sortie de fonction et non directement dans la fonction
  lambda_t <- lambdat_fct(N, t, param)
  B <- B_fct(TT - t, param)
  A <- A_fct(TT - t, param)
  return(A * exp(-B * lambda_t))
}

Prix_Spread_fct <- function(t,TT,param,LGD){
  return((1 - LGD + LGD * survie_FF_fct(N,t,TT,param) )**(-1/(TT-t))-1)
}

#Attention il n'y a pas de N en paramètre ...
Prix_Spread_fct_sim <- function(t,TT,param,LGD){
  return((1 - LGD + LGD * survie_FF_fct_sim(N,t,TT,param) )**(-1/(TT-t))-1)
}


#### Fonction ?cart ? la moyenne pour le calibrage ####
Ecart_CIR <- function(param){
  e <- 0
  for (t in 0:TT){
    spread <- Prix_Spread_fct(0,t,param,LGD)
    e <- e + (spread-spread_mkt)^2
  }
  return(e)
}

#====== ESSAI =====
LGD <- 0.3 # car 1-R=LGD
spread_mkt <- SpreadMarket[6]
TT <- 15
Ecart_CIR(param_test)

LB = c(-2,0,0,0)
UB = c(2,5,1,1)
param <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
print(param)

#param T = 150 : k = 2.257538e-03, mu = 3.841823e+00, sigma = 3.814697e-06, lambda0 = 5.140305e-02
param_150 <- c(2.257538e-03, 3.841823e+00, 3.814697e-06, 5.140305e-02)
plot(1:15,Prix_Spread_fct(0, 1:15, param_150, LGD), main='Spread T = 150', ylab = 'Spread', ylim=c(0.005,0.03))
abline(h=spread_mkt, col='red')

#param T = 15 : k = 5.923462e-03, mu = 9.856308e-01, sigma = 3.814697e-06, lambda0 = 5.830078e-02
param_15 <- c(5.923462e-03, 9.856308e-01, 3.814697e-06, 5.830078e-02)
plot(1:15,Prix_Spread_fct(0, 1:15, param_15, LGD), main='Spread T = 15', ylab = 'Spread', ylim=c(0.005,0.03))
abline(h=spread_mkt, col='red')

#============= RESULTATS ==============
#On va entrainer notre modèle uniquement sur T = 15
TT <- 15
#AAA :
spread_mkt <- SpreadMarket[1]
param_AAA <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#AA :
spread_mkt <- SpreadMarket[2]
param_AA <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#A
spread_mkt <- SpreadMarket[3]
param_A <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#BBB
spread_mkt <- SpreadMarket[4]
param_BBB <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#BB
spread_mkt <- SpreadMarket[5]
param_BB <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#B
spread_mkt <- SpreadMarket[6]
param_B <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par
#NR 
spread_mkt <- SpreadMarket[7]
param_NR <- hjkb(c(0.2,0.2,0.2,0.2),Ecart_CIR,lower=LB,upper=UB)$par

plot(1:15,Prix_Spread_fct(0, 1:15, param_AAA, LGD), main='Spread reproduits par le modèle CIR (T = 15)', type='l', ylab="Spread", xlab="Maturité", ylim=c(0.0005, 0.025), col='red')
lines(Prix_Spread_fct(0, 1:15, param_AA, LGD), col='orange')
lines(Prix_Spread_fct(0, 1:15, param_A, LGD), col='brown')
lines(Prix_Spread_fct(0, 1:15, param_BBB, LGD), col='lightblue')
lines(Prix_Spread_fct(0, 1:15, param_BB, LGD), col='blue')
lines(Prix_Spread_fct(0, 1:15, param_B, LGD), col='purple')
#lines(Prix_Spread_fct(0, 1:15, param_NR, LGD), col='blue')

legend("topright", c("AAA", "AA", "A", "BBB", "BB", "B"),
       col = c("red", "orange", "brown", "lightblue", "blue", "purple"), lty = c(1, 1, 1), cex = 0.6)
#Retour sur le AAA :
plot(1:15,Prix_Spread_fct(0, 1:15, param_AA, LGD), main='Spread AA (par CIR) comparé à la données (T = 15)', ylab="Spread", xlab="Maturité", ylim=c(0.001, 0.006), col='red')
abline(h = SpreadMarket[2])

print("    k            mu         sigma       lambda0")
print(param)

#### Tests de martingalité ####
spread.act <- c()
N <- 1000
print(Prix_Spread_fct_sim(0,2,param_AAA,LGD))
print(survie_FF_fct_sim(N, 0, 2, param_AAA))
for (t in 0:150){spread.act <- c(spread.act, mean(Prix_Spread_fct_sim(t,TT,param,LGD)))}
plot(spread.act,pch=20,col="darkgrey",
     main="Test de martingalité pour les actions",
     xlab="Maturité",
     ylab="Moyenne de l'indice actualisé ")
abline(h=S0,col="red",lty=3,lwd=2)

abline(h=S0+0.025,col="blue",lty=4)
abline(h=S0-0.025,col="blue",lty=4)
legend("topleft",legend=c("Prix initial","Borne à 5%"),
       col=c("red","blue"),pch=20,
       cex=0.8)

plot(TauxZC)
