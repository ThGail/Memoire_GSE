rm(list = ls())
library(readxl)

#### Importation des donnees ###################################################
setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])

# parametres taux Vasicek
a.vas <- 0.01298906
b.vas <- 0.03872287
sigma.vas <- 0.002397156
param <- c(a.vas,b.vas,sigma.vas)
r0 <- TauxZC[1]
# parametres actifs Black-Scholes
S0 <- 1
K <- 1
sigma.action <- 0.152396
sigma.immo <- 0.10
# nombre de simulation de Monte Carlo
N <- 50
TT <- length(TauxZC)
# matrice de correlaion entre taux_court, action_glo et immo
Corr <- matrix(c(1,0.1,0.05,0.1,1,-0.3,0.05,-0.3,1),nrow=3)

# fonctions interm?diaires, qui ? chaque Brownien donn?, ? chaque maturit? donn?e,  
# nous donne des sc?narios de taux Vasicek
Vas_W <- function(t, r0, param, W) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  X <- r0 * exp(-a * t) + b * (1 - exp(-a * t)) + sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a)) * W
  return(X)
}

# nous donne des sc?narios de taux d'actif par Black-Scholes
BS_W <- function(t, S0, r0, param, sigma.act, W.vas, W.act) {
  X <- S0 * exp((Vas_W(t, r0, param, W.vas) - 0.5 * sigma.act^2) * t + sigma.act * sqrt(t) * W.act)
  return(X)
}

# la fonction principale qui donne le GSE t?moin avec les corr?lations,
# pour une maturit? donn?e
GSE_mat_temoin <- function(N, t, r0, param, sigma.action, sigma.immo, C){
  sigma.vas <- param[3]
  
  # simulation de mouvements browniens corr?l?s
  Z <- matrix(rnorm(3 * N), nrow = 3)
  L <- t(chol(C, pivot = FALSE))
  W <- L %*% Z
  W.vas <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  vas <- Vas_W(t, r0, param, W.vas)
  action <- BS_W(t, S0, r0, param, sigma.action, W.vas, W.action)
  immo <- BS_W(t, S0, r0, param, sigma.immo, W.vas, W.immo)
  
  return(cbind(vas,action,immo))
}

vas.test <- c()
action.test <- c()
immo.test <- c()
for (t in 1:TT){
  test <- GSE_mat_temoin(N, t, r0, param, sigma.action, sigma.immo, Corr)
  vas.test <- cbind(vas.test, test[,1])
  action.test <- cbind(action.test, test[,2])
  immo.test <- cbind(immo.test, test[,3])
}
# plot du taux sans risque Vasicek avec matrice de corr?lation
matplot(t(vas.test), type='l', main='Simulation Taux ZC corrélé vs Taux ZC réel', ylab='Taux ZC', xlab='Maturité')
lines(TauxZC, col='black', lw=4)
lines(colMeans(vas.test),  type = "l", col="red", lw=3)
legend("topleft", legend=c("Moyenne des simulations", 'Taux Réel'),
       col=c("red" ,"black"),pch=20,
       cex=0.8)

#Plot du prix action :
matplot(t(action.test), type='l', xlim=c(0,100), ylim=c(0,1000),  main='Simulation Prix Action', ylab='Prix Action', xlab='Maturité')
lines(colMeans(action.test), type = 'l', col="black", xlim=c(0,100), ylim=c(0,1000), lw = 4)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("black"),pch=20,
       cex=1)

#Plot du prix Immobilier :
matplot(t(immo.test), type='l', xlim=c(0,100), ylim=c(0,1000), main='Simulation Prix Immobilier', ylab='Prix Immobilier', xlab='Maturité')
lines(colMeans(immo.test), type = 'l', col="black", xlim=c(0,100), ylim=c(0,1000), lw = 4)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("black"),pch=20,
       cex=1)

#Transformation du taux zéro coupon en prix zéro coupon :

#Le prix Zéro-Coupon réel : 
#On transforme les données issues du GSE de la même manière :
matplot(exp(- Maturite * t(vas.test)), type='l', main='Prix ZC simulé', ylab='Prix ZC', xlab='Maturité', ylim=c(0,15))
lines(exp(-Maturite * TauxZC), type = "l", xlab = "Maturité", ylab = "Prix ZC",  col = 'black', lw=4)
legend("topleft", legend=c("Prix réel via Formule de Transition"),
       col=c("black"),pch=20,
       cex=0.8)

