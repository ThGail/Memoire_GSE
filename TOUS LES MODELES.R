rm(list = ls())
library(readxl)
library(dfoptim)

### Importation des donnees ###################################################
setwd("~/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)
sheetEIOPA <- read_excel("RFR_for_publication_RFR_curves.xlsx", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])
Maturite10 <- Maturite[1:10]
Maturite30 <- Maturite[1:30]





### PARTIE VASICEK #############################################################
#### Définition des fonctions Vasicek ###################################################
# simulation taux instantané Vasicek
taux_Vas_sim <- function(N, t, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(rnorm(N, r0 * exp(-a * t) + b * (1 - exp(-a * t)), sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a))))
}
taux_Vas_sim.t <- Vectorize(taux_Vas_sim,"t")

# simulation TZC
TZC_Vas_FF_sim <- function(N, t, TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  nb_mat <- length(TT)
  theta <- t(matrix(rep(TT,N),nrow=nb_mat))-t
  taux_vas <- matrix(rep(taux_Vas_sim(N,t,param,r0),nb_mat),ncol=nb_mat)
  return(Ri - ((Ri-taux_vas)*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}

# simulation TZC des trajectoires pour une maturité donnée
TZC_Vas_FF_sim.t <- Vectorize(TZC_Vas_FF_sim,"t")

# formule fermée du TZC pour le calibrage
TZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  return(Ri - ((Ri-r0)*(1-exp(-a*TT))-sigma^2/(4*a^2)*(1-exp(-a*TT))^2)/(a*TT))
}

# simulation PZC sur les trajectoires
PZC_Vas_FF_sim <- function(N, t, TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  nb_mat <- length(TT)
  theta <- t(matrix(rep(TT,N),nrow=nb_mat))-t
  taux_vas <- matrix(rep(taux_Vas_sim(N,t,param,r0),nb_mat),ncol=nb_mat)
  return(exp(-b * theta) * exp(-(taux_vas - b) * (1 - exp(-a * theta)) / a + 0.5 * (sigma^2 * theta / a^2 - sigma^2 / a^3 * (1 - exp(-a * theta)) - sigma^2 / (2 * a^3) * (1 - exp(-a * theta))^2)))
}

# simulation PZC sur les trajectoires pour une maturité donnée
PZC_Vas_FF_sim.t <- Vectorize(PZC_Vas_FF_sim,"t")

# formule fermée du PZC (pour le calibrage)
PZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b*TT) * exp(-(r0-b)*(1-exp(-a*TT))/a + 0.5 * (sigma^2 * TT / a^2 - sigma^2 / a^3 * (1 - exp(-a * TT)) - sigma^2 / (2 * a^3) * (1 - exp(-a * TT))^2)))
}

# fonction objective à minimiser
ecart_Vas <- function(param){
  e <- 0
  for (t in 1:150){
    tzp<-TZC_Vas_FF_calibrage(t, param, TauxZC[1])
    # tzp <- mean(taux_Vas_sim.t(1000,t,paramVas))
    e <- e + (tzp-TauxZC[t])^2
  }
  return(e)
}





#### Calibrage du modèle Vasicek ###################################################
param_init = c(0.005,0.005,0.005)
LB = c(0,0,1e-16)
UB = c(1,1,2)
(paramVas = hjkb(param_init,ecart_Vas,lower=LB,upper=UB)$par)
# 0.004084473 0.292354584 0.003962402
# erreur = 0.0005183257





#### Some plot Vasicek ###################################################
# plot calibrage TZC
plot(TauxZC, type="l", main="Calibrage du taux zéro-coupon",col="red",lty=3,lwd=2,
     xlab="Maturité",ylab="Taux zéro-coupon")
lines(TZC_Vas_FF_calibrage(Maturite,paramVas),lwd=2)
legend("bottomright",legend=c("Courbe EIOPA (input)","Courbe TZC du modèle"),
       col=c("red","black"),pch=20,
       cex=0.6)

# plot PZC
plot(exp(-Maturite*TauxZC), type="l", main="Vérification avec le prix zéro-coupon",col="red",lty=3,lwd=2,
     xlab="Maturité",ylab="Prix zéro-coupon")
lines(PZC_Vas_FF_calibrage(Maturite,paramVas),lwd=2)
legend("topright",legend=c("Courbe issue de la EIOPA","Courbe PZC du modèle"),
       col=c("red","black"),pch=20,
       cex=0.6)

# simulation TZC de maturité 10 ans (1000 simulations)
tt <- seq(0,10,0.1)[-length(seq(0,10,0.1))]
plotTZC <- TZC_Vas_FF_sim.t(1000, tt, 10, paramVas)
matplot(tt,t(plotTZC[1:10,]),type="l",
        main="Scénarios de taux zéro-coupon de maturité 10 ans",
        xlab="Temps",ylab="TZC")
plot(colMeans(plotTZC),type="l",lwd=2,col="red",
     main="Moyenne de taux zéro-coupon de maturité 10 ans",
     xlab="Temps",ylab="TZC")

# simulation TZC de maturité 10 ans (1000 simulations)
plotPZC <- PZC_Vas_FF_sim.t(1000, tt, 10, paramVas)
matplot(tt,t(plotPZC[1:10,]),type="l",
        main="Scénarios du prix zéro-coupon de maturité 10 ans",
        xlab="Temps",ylab="PZC")
plot(colMeans(plotPZC),type="l",lwd=2,col="red",
     main="Moyenne du prix zéro-coupon de maturité 10 ans",
     xlab="Temps",ylab="TZC")





### PARTIE BLACK AND SCHOLES #############################################################
#### Définition des fonctions BS ###################################################
S0 = 1; K = 1
paramBS <- list(action=0.152396,immo=0.10)

# simulation prix actif
PA_BS_FF_sim <- function(N, t, S0, param_taux, param_actif, r0=TauxZC[1]){
  nb_mat <- length(t)
  mat <- t(matrix(rep(t,N),nrow=nb_mat))
  r <- taux_Vas_sim(N, t, param_taux, r0)
  return(S0 * exp((r - 0.5 * param_actif^2) * t + param_actif * rnorm(N, 0, sqrt(t))))
}

# simulation prix actif sur plusieurs temps
PA_BS_FF_sim.t <- Vectorize(PA_BS_FF_sim,"t")

# simulation prix de l'option Call actualisé
# call_BS_sim <- function(N,t,TT,S0,param_taux,param_actif,r0=TauxZC[1]){
#   tau <- TT - t
#   nb_mat <- length(tau)
#   theta <- t(matrix(rep(tau,N),nrow=nb_mat))
#   r <- taux_Vas_sim.t(N, tau, param_taux, r0)
#   return(exp(-r*theta)*pmax(PA_BS_FF_sim(N, TT, S0, param_taux, param_actif, r0)-K,0))
# }
# (z = call_BS_sim(5,0,Maturite10,S0,paramVas,paramBS$action))





#### Some plot BS ###################################################
# moyenne simulation indice action pour 1000 simulations
(plotPAaction = PA_BS_FF_sim.t(1000,Maturite30,S0,paramVas,paramBS$action))
matplot(t(plotPAaction[1:10,]),type="l",
        main="Scénarios de l'indice action sur 30 ans",
        xlab="Temps",ylab="indice action")
lines(colMeans(plotPAaction),type="l",lwd=2,col="red")

# moyenne simulation indice immo pour 1000 simulations
(plotPAimmo = PA_BS_FF_sim.t(1000,Maturite30,S0,paramVas,paramBS$immo))
matplot(t(plotPAimmo[1:10,]),type="l",
        main="Scénarios de l'indice immobilier sur 30 ans",
        xlab="Temps",ylab="indice immobilier")
lines(colMeans(plotPAimmo),type="l",lwd=2,col="red")



