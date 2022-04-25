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
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
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
  for (t in 1:length(TauxZC)){
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

# simulation prix actif pour le test de martingalité
PAact_BS_FF_sim <- function(N,t,S0,param_taux, param_actif, r0=TauxZC[1]){
  r <- TZC_Vas_FF_calibrage(t, param_taux, r0)
  return(exp(-r*t)*S0*exp((r-0.5*param_actif^2)*t+param_actif*rnorm(N,0,sqrt(t))))
}
PAact_BS_FF_sim.t <- Vectorize(PAact_BS_FF_sim,"t")

# simulation prix de l'option Call actualisé pour une maturité donnée
call_BS_sim <- function(N,t,TT,S0,param_taux,param_actif,r0=TauxZC[1]){
  tau <- TT - t
  nb_mat <- length(tau)
  theta <- t(matrix(rep(tau,N),nrow=nb_mat))
  r <- taux_Vas_sim.t(N, t, param_taux, r0)
  return(exp(-r*theta)*pmax(PA_BS_FF_sim(N, TT, S0, param_taux, param_actif, r0)-K,0))
}
call_BS_sim.T <- Vectorize(call_BS_sim,"TT")

call_BS_FF <- function(TT,S0,K,param_taux,param_actif,r0=TauxZC[1]) {
  r <- TZC_Vas_FF_calibrage(TT,param_taux,r0)
  d1 <- (log(S0/K) + (r + 0.5*param_actif^2)*TT)/(param_actif*sqrt(TT))
  d2 <- d1 - param_actif*sqrt(TT)
  return(S0*pnorm(d1) - K*exp(-r*TT)*pnorm(d2))
}





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

# convergence du call issu du prix action pour maturité 1 an
plotcallaction1an <- c()
for (n in seq(1, 5000, 100)){
  plotcallaction1an <- c(plotcallaction1an, mean(call_BS_sim(n,0,1,S0,paramVas,paramBS$action)))
}
plot(seq(1,5000,100),plotcallaction1an, 
     main=c("Convergence du call issu du prix action pour maturité 1 an"),pch=20,
     xlab="Nombre de simulations",ylab="Prix call")
abline(h=call_BS_FF(1,S0,K,paramVas,paramBS$action),lty=2,lwd=3,col="red")

# convergence du call issu du prix immo pour maturité 1 an
plotcallimmo1an <- c()
for (n in seq(1, 5000, 100)){
  plotcallimmo1an <- c(plotcallimmo1an, mean(call_BS_sim(n,0,1,S0,paramVas,paramBS$immo)))
}
plot(seq(1, 5000, 100),plotcallimmo1an, 
     main=c("Convergence du call issu du prix immobilier pour maturité 1 an"),pch=20,
     xlab="",ylab="Prix call")
abline(h=call_BS_FF(1,S0,K,paramVas,paramBS$immo),lty=2,lwd=3,col="red")

# test de martingalité sur indice action
plotaction = PAact_BS_FF_sim.t(3e5,Maturite,S0,paramVas,paramBS$action)
plot(colMeans(plotaction),ylim=c(0.9,1.1),pch=20,col="darkgrey",
     main="Test de martingalité pour les actions",
     xlab="Maturité",
     ylab="Moyenne de l'indice actualisé ")
abline(h=S0,col="red",lty=3,lwd=2)
abline(h=S0+0.025,col="blue",lty=4)
abline(h=S0-0.025,col="blue",lty=4)
legend("topleft",legend=c("Prix initial","Borne à 5%"),
       col=c("red","blue"),pch=20,
       cex=0.8)

# test de martingalité sur indice immo
plotimmo = PAact_BS_FF_sim.t(1e5,Maturite,S0,paramVas,paramBS$immo)
plot(colMeans(plotimmo),ylim=c(0.9,1.1),pch=20,col="darkgrey",
     main="Test de martingalité pour les immobiliers",
     xlab="Maturité",
     ylab="Moyenne de l'indice actualisé ")
abline(h=S0,col="red",lty=3,lwd=2)
abline(h=S0+0.025,col="blue",lty=4)
abline(h=S0-0.025,col="blue",lty=4)
legend("topleft",legend=c("Prix initial","Borne à 5%"),
       col=c("red","blue"),pch=20,
       cex=0.8)





### PARTIE CIR #############################################################
#### Définition des fonctions CIR ###################################################
A_CIR <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_CIR <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return((2 * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# calcule de l'intensite de défaut
lambdat_CIR_sim <- function(N, t, param) {
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

# simulation de la proba de survie
survie_CIR_sim <- function(N, t, TT, param) {
  lambda_t <- lambdat_CIR_sim(N, t, param)
  B <- B_CIR(TT - t, param)
  A <- A_CIR(TT - t, param)
  return(matrix(A*exp(- lambda_t %*% t(B)), nrow = N))
}
survie_CIR_sim.t <- Vectorize(survie_CIR_sim,"t")

# survie en version FF pour la calibration
survie_CIR_calibrage <- function(TT, param) {
  # le calibrage se fait en t=0
  lambda0 <- param[4]
  B <- B_CIR(TT, param)
  A <- A_CIR(TT, param)
  return(A*exp(-B*lambda0))
}

# le spread pour chaque rating
spread_CIR_sim <- function(N,t,TT,param,LGD){
  return((1 - LGD + LGD * survie_CIR_sim(N,t,TT,param) )**(-1/(TT-t))-1)
}
spread_CIR_sim.t <- Vectorize(spread_CIR_sim,"t")


# spread en version FF pour la calibration
spread_CIR_FF_calibrage <- function(TT,param,LGD){
  return((1-LGD+LGD*survie_CIR_calibrage(TT,param))**(-1/TT)-1)
}


# fonction objective à minimiser pour le modèle CIR
ecart_CIR <- function(param){
  e <- 0
  for (t in Maturite[1:TT]){
    spread <- spread_CIR_FF_calibrage(t,param,LGD)
    e <- e + (spread-spread_mkt)^2
  }
  return(e)
}




#### Calibrage du modèle CIR ###################################################
# on a uniquement les spreads de maturité 1 an
TT = 1; LGD = 0.3

paramCIR <- list()
param_init = c(0.2,0.2,0.2,0.2)
LB = c(0,0,0,0)
UB = c(10,10,10,10)

#AAA : 0.00781250 0.21624298 0.05812742 0.06718750
spread_mkt <- SpreadMarket[1]
paramCIR$AAA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AAA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AAA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AAA[3] <- min(paramCIR$AAA[3],sqrt(2*paramCIR$AAA[1]*paramCIR$AAA[2]))

#AA : 0.125000 0.137500 0.185405 0.075000
spread_mkt <- SpreadMarket[2]
paramCIR$AA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AA[3] <- min(paramCIR$AA[3],sqrt(2*paramCIR$AA[1]*paramCIR$AA[2]))

#A : 0.0312500 1.2156250 0.2756385 0.0781250
spread_mkt <- SpreadMarket[3]
paramCIR$A <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$A <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$A <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$A[3] <- min(paramCIR$A[3],sqrt(2*paramCIR$A[1]*paramCIR$A[2]))

#BBB : 0.00390625 0.16676254 0.03609477 0.19609375
spread_mkt <- SpreadMarket[4]
paramCIR$BBB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BBB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BBB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BBB[3] <- min(paramCIR$BBB[3],sqrt(2*paramCIR$BBB[1]*paramCIR$BBB[2]))

#BB : 1.2000000 0.2000000 0.6928203 0.2000000
spread_mkt <- SpreadMarket[5]
paramCIR$BB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BB[3] <- min(paramCIR$BB[3],sqrt(2*paramCIR$BB[1]*paramCIR$BB[2]))

#B : 1.6687500 0.2000000 0.8170067 0.2000000
spread_mkt <- SpreadMarket[6]
paramCIR$B <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$B <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$B <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$B[3] <- min(paramCIR$B[3],sqrt(2*paramCIR$B[1]*paramCIR$B[2]))





#### Some plot CIR ###################################################
# loi de lambda est suit une Chi-2 ? 
hist(lambdat_CIR_sim(1000, 3, param_init), freq=FALSE, main='Histogramme des simulations (param random)', xlab='Simulation des lambda t = 3')

# plot du spread pour tous les rating sur 30 ans
plot(Maturite30,spread_CIR_FF_calibrage(Maturite30, paramCIR$AAA, LGD), main='Les spread de crédit reproduits par le modèle CIR', type='l', ylab="Spread", xlab="Maturité", ylim=c(0.001, 0.03), col='red')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$AA, LGD), col='orange')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$A, LGD), col='brown')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$BBB, LGD), col='lightblue')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$BB, LGD), col='blue')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$B, LGD), col='purple')
legend("topright", c("AAA", "AA", "A", "BBB", "BB", "B"),
       col = c("red", "orange", "brown", "lightblue", "blue", "purple"), lty = c(1, 1, 1), cex = 0.6)

# simulation spread AAA pour une maturité de 30 ans
tt <- seq(1,30,0.5)[-length(seq(1,30,0.5))]
plotCIRsurvieAAA <- spread_CIR_sim.t(1000,tt,30,paramCIR$AAA,LGD)
matplot(tt,t(plotCIRsurvieAAA[1:50,]),type="l",
        main="Scénarios de spread AAA de maturité 30 ans",
        xlab="Temps",ylab="Spread")
lines(colMeans(plotCIRsurvieAAA),type="l",lwd=2,col="red")



