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
PrixZC <- exp(-Maturite*TauxZC)
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
TZC_Vas_FF_calibrage.T <- Vectorize(TZC_Vas_FF_calibrage,"TT")

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
PZC_Vas_FF_calibrage.T <- Vectorize(PZC_Vas_FF_calibrage,"TT")

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
N = 1000

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
plotTZC <- TZC_Vas_FF_sim.t(N, tt, 10, paramVas)
matplot(tt,t(plotTZC[1:10,]),type="l",
        main="Scénarios de taux zéro-coupon de maturité 10 ans",
        xlab="Temps",ylab="TZC")
plot(tt,colMeans(plotTZC),type="l",lwd=2,col="red",
     main="Moyenne de taux zéro-coupon de maturité 10 ans",
     xlab="Temps",ylab="TZC")

# simulation TZC de maturité 10 ans (1000 simulations)
plotPZC <- PZC_Vas_FF_sim.t(N, tt, 10, paramVas)
matplot(tt,t(plotPZC[1:10,]),type="l",
        main="Scénarios du prix zéro-coupon de maturité 10 ans",
        xlab="Temps",ylab="PZC")
plot(tt,colMeans(plotPZC),type="l",lwd=2,col="red",
     main="Moyenne du prix zéro-coupon de maturité 10 ans",
     xlab="Temps",ylab="TZC")





### PARTIE HULL AND WHITE #############################################################





### PARTIE BLACK AND SCHOLES #############################################################
#### Définition des fonctions BS ###################################################
S0 = 1; K = 1
paramBS <- list(action=0.152396,immo=0.10)

# simulation prix actif
PA_BS_FF_sim <- function(N, t, S0, param_Vas, param_actif, r0=TauxZC[1]){
  nb_mat <- length(t)
  mat <- t(matrix(rep(t,N),nrow=nb_mat))
  r <- taux_Vas_sim(N, t, param_Vas, r0)
  return(S0 * exp((r - 0.5 * param_actif^2) * t + param_actif * rnorm(N, 0, sqrt(t))))
}
# simulation prix actif sur plusieurs temps
PA_BS_FF_sim.t <- Vectorize(PA_BS_FF_sim,"t")

# simulation prix actif pour le test de martingalité
PAact_BS_FF_sim <- function(N,t,S0,param_Vas, param_actif, r0=TauxZC[1]){
  r <- TZC_Vas_FF_calibrage(t, param_Vas, r0)
  return(exp(-r*t)*S0*exp((r-0.5*param_actif^2)*t+param_actif*rnorm(N,0,sqrt(t))))
}
PAact_BS_FF_sim.t <- Vectorize(PAact_BS_FF_sim,"t")

# simulation prix de l'option Call actualisé pour une maturité donnée
call_BS_sim <- function(N,t,TT,S0,param_Vas,param_actif,r0=TauxZC[1]){
  tau <- TT - t
  nb_mat <- length(tau)
  theta <- t(matrix(rep(tau,N),nrow=nb_mat))
  r <- taux_Vas_sim.t(N, t, param_Vas, r0)
  return(exp(-r*theta)*pmax(PA_BS_FF_sim(N, TT, S0, param_Vas, param_actif, r0)-K,0))
}
call_BS_sim.T <- Vectorize(call_BS_sim,"TT")

call_BS_FF <- function(TT,S0,K,param_Vas,param_actif,r0=TauxZC[1]) {
  r <- TZC_Vas_FF_calibrage(TT,param_Vas,r0)
  d1 <- (log(S0/K) + (r + 0.5*param_actif^2)*TT)/(param_actif*sqrt(TT))
  d2 <- d1 - param_actif*sqrt(TT)
  return(S0*pnorm(d1) - K*exp(-r*TT)*pnorm(d2))
}





#### Some plot BS ###################################################
N = 1000

# moyenne simulation indice action pour 1000 simulations
(plotPAaction = PA_BS_FF_sim.t(N,Maturite30,S0,paramVas,paramBS$action))
matplot(t(plotPAaction[1:10,]),type="l",
        main="Scénarios de l'indice action sur 30 ans",
        xlab="Temps",ylab="indice action")
lines(colMeans(plotPAaction),type="l",lwd=2,col="red")

# moyenne simulation indice immo pour 1000 simulations
(plotPAimmo = PA_BS_FF_sim.t(N,Maturite30,S0,paramVas,paramBS$immo))
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
survie_CIR_sim.T <- Vectorize(survie_CIR_sim,"TT")

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
spread_CIR_sim.T <- Vectorize(spread_CIR_sim,"TT")

# spread en version FF pour la calibration
spread_CIR_FF_calibrage <- function(TT,param,LGD){
  return((1-LGD+LGD*survie_CIR_calibrage(TT,param))**(-1/TT)-1)
}

# prix zéro-coupon risqué en version FF pour chaque rating à la maturité
PZCr_CIR_FF <- function(TT, param_Vas, param_CIR, LGD){
  return(PZC_Vas_FF_calibrage(TT, param_Vas)/(1+spread_CIR_FF_calibrage(TT,param_CIR,LGD))^TT)
}

# simulation zéro-coupon risqué
PZCr_CIR_sim <- function(N,t,TT,param_Vas,param_CIR,LGD){
  surv <- survie_CIR_sim(N, t, TT, param_CIR)
  TZC <- TZC_Vas_FF_sim(N, t, TT, param_Vas)
  return(exp(-TT*t(TZC))*t(surv+(1-LGD)*(1-surv)))
}
# pour une maturité T fixé
PZCr_CIR_sim.t <- Vectorize(PZCr_CIR_sim,"t")
# pour t=0, à maturité variable
PZCr_CIR_sim.T <- Vectorize(PZCr_CIR_sim,"TT")

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

#AAA : 1.599722e-01 1.288132e-01 2.030100e-01 2.126438e-13
spread_mkt <- SpreadMarket[1]
#paramCIR$AAA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AAA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$AAA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AAA[3] <- min(paramCIR$AAA[3],sqrt(2*paramCIR$AAA[1]*paramCIR$AAA[2]))

#AA : 0.1793781 0.1526860 0.2028968 0.0000000
spread_mkt <- SpreadMarket[2]
#paramCIR$AA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$AA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AA[3] <- min(paramCIR$AA[3],sqrt(2*paramCIR$AA[1]*paramCIR$AA[2]))

#A : 0.1960534 0.1731844 0.2027408 0.0000000
spread_mkt <- SpreadMarket[3]
#paramCIR$A <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$A <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$A <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$A[3] <- min(paramCIR$A[3],sqrt(2*paramCIR$A[1]*paramCIR$A[2]))

#BBB : 0.19780523 0.17654727 0.20257493 0.01274617
spread_mkt <- SpreadMarket[4]
#paramCIR$BBB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BBB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$BBB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BBB[3] <- min(paramCIR$BBB[3],sqrt(2*paramCIR$BBB[1]*paramCIR$BBB[2]))

#BB : 0.19887557 0.18437281 0.20176218 0.04921644
spread_mkt <- SpreadMarket[5]
#paramCIR$BB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$BB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BB[3] <- min(paramCIR$BB[3],sqrt(2*paramCIR$BB[1]*paramCIR$BB[2]))

#B : 0.19890474 0.18459973 0.20173836 0.05140548
spread_mkt <- SpreadMarket[6]
#paramCIR$B <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$B <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$B <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$B[3] <- min(paramCIR$B[3],sqrt(2*paramCIR$B[1]*paramCIR$B[2]))





#### Some plot CIR ###################################################
N = 1000

# loi de lambda est suit une Chi-2 ? 
hist(lambdat_CIR_sim(N, 3, param_init), freq=FALSE, main='Histogramme des simulations (param random)', xlab='Simulation des lambda t = 3')

# plot du spread pour tous les rating sur 30 ans
plot(Maturite30,spread_CIR_FF_calibrage(Maturite30, paramCIR$AAA, LGD), main='Les spread de crédit reproduits par le modèle CIR', type='l', ylab="Spread", xlab="Maturité", ylim=c(0.001, 0.03), col='red')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$AA, LGD), col='orange')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$A, LGD), col='brown')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$BBB, LGD), col='lightblue')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$BB, LGD), col='blue')
lines(spread_CIR_FF_calibrage(Maturite30, paramCIR$B, LGD), col='purple')
legend("topright", c("AAA", "AA", "A", "BBB", "BB", "B"),
       col = c("red", "orange", "brown", "lightblue", "blue", "purple"), lty = c(1, 1, 1), cex = 0.6)
# remarque : on trouve pareil avec 
# plot(colMeans(spread_CIR_sim.T(N,0,Maturite30,paramCIR$AAA,LGD)))

# simulation spread AAA pour une maturité de 30 ans
tt <- seq(1,30,0.5)[-length(seq(1,30,0.5))]
plotCIRsurvieAAA <- spread_CIR_sim.t(N,tt,30,paramCIR$AAA,LGD)
matplot(tt,t(plotCIRsurvieAAA[1:10,]),type="l",
        main="Scénarios de spread AAA de maturité 30 ans",
        xlab="Temps",ylab="Spread")
# plot moyenne spread de maturité 30 ans avec Monte Carlo (sur 1000 simulations)
lines(tt,colMeans(plotCIRsurvieAAA),type="l",lwd=2,col="red")





#### (Test de martingalité CIR) ###################################################
plot(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$AAA, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$AAA, LGD)),"l",
     ylab="CashFlow Actualisé",ylim=c(0.99,1.01),
     main="Test de martingalité CIR", 
     col = 'red')
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$AA, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$AA, LGD)),"l",col="orange")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$A, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$A, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$BBB, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$BBB, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$BB, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$BB, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$B, LGD)/colMeans(PZCr_CIR_sim.T(N, 0, Maturite, paramVas, paramCIR$B, LGD)),"l",col="brown")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)





### PARTIE JLT #############################################################
#### Prépration donnée historique JLT #########################################
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

# diagonalisation de la matrice L
M = eigen(L)$vectors
D = diag(eigen(L)$values)
# round(M%*%D%*%solve(M),15) permet de retrouver la matrice de base (modulo erreur numérique)





#### Définition des fonctions JLT ###################################################
A_JLT <- function(u, param, dj) {
  k = param[1:8]
  mu = param[9:16]
  sigma = param[17]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_JLT <- function(u, param, dj) {
  k = param[1:8]
  mu = param[9:16]
  sigma = param[17]
  ga <- sqrt(k**2 - 2 * dj * sigma**2)
  return((- 2 * dj * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# la dynamique pour la proba de défaut
pit_JLT <- function(N, t, param) {
  k <- matrix(rep(param[1:8], N), nrow = 8)
  mu <- matrix(rep(param[(1:8)+8], N), nrow = 8)
  sigma <- matrix(rep(param[17], N*8), nrow = 8)
  pi_t <- matrix(rep(param[18:25], N), nrow = 8)
  
  if (t==0) {return(pi_t)}
  for (i in 1:t){
    pi_t <- abs(pi_t + k * (mu - pi_t) + t(t(sigma * sqrt(pi_t)) * rnorm(N, 0, 1)))}
  return(pi_t)
}

# la proba de faire défaut à partir du rating i
proba_defaut_i_JLT <- function(N, t, TT, param, M, D, i){
  # pour un rating donnée, on fait N simulations
  
  invM <- solve(M)
  K <- length(D[1,])
  
  sum <- 0
  for (j in (1:(K-1))){
    esp <- (A_JLT(TT-t, param, D[j,j])*exp(-B_JLT(TT-t, param, D[j,j])*pit_JLT(N,t,param)))[i,]
    sum <- sum + M[i,j]*invM[j,K]*(esp-1)
  }
  return(sum)
}
proba_defaut_i_JLT.T <- Vectorize(proba_defaut_i_JLT,"TT")

# pareil en version calibrage
proba_defaut_i_JLT_calibrage <- function(TT, param, M, D, i){
  # pour un rating donnée, et une maturité TT donnée, on a juste besion de t=0
  
  invM <- solve(M)
  K <- length(D[1,])
  pi0 = param[18:25]
  
  sum <- 0
  for (j in (1:(K-1))){
    esp <- A_JLT(TT, param, D[j,j])[i]*exp(-B_JLT(TT, param, D[j,j])[i]*pi0[i])
    sum <- sum + M[i,j]*invM[j,K]*(esp-1)
  }
  return(sum)
}
proba_defaut_i_JLT_calibrage.T <- Vectorize(proba_defaut_i_JLT_calibrage,"TT")

# formule générale de la proba de passage
proba_passage_ij_JLT <- function(N, t, TT, param, M, D, i, j){
  # pour un rating donnée, on fait N simulations
  
  invM <- solve(M)
  K <- length(D[1,])
  
  sum <- 0
  for (k in (1:(K-1))){
    esp <- (A_JLT(TT-t, param, D[k,k])*exp(-B_JLT(TT-t, param, D[k,k])*pit_JLT(N,t,param)))[i,]
    sum <- sum + M[i,k]*invM[k,j]*esp
  }
  
  sum <- sum + M[i,K]*invM[K,j]
  return(sum)
}

# calcule du spread de rating i
spread_i_JLT <- function(N, t, TT, param, M, D, i, LGD){
  if (t==TT){return(rep(0,N))}
  # if (1- LGD * proba_defaut_i_JLT(N, t, TT, param, M, D, i)<=0){return(0)}
  return(-(1/(TT-t))*log(1- LGD * proba_defaut_i_JLT(N, t, TT, param, M, D, i)))
}
spread_i_JLT.t <- Vectorize(spread_i_JLT,"t")
spread_i_JLT.T <- Vectorize(spread_i_JLT,"TT")

# spread pour la calibration
spread_i_JLT_calibrage <- function(TT, param, M, D, i, LGD){
  if (TT==0){return(0)}
  return(-log(1 - LGD * proba_defaut_i_JLT_calibrage(TT, param, M, D, i))/TT)
}
spread_i_JLT_calibrage.T <- Vectorize(spread_i_JLT_calibrage,"TT")





#### Calibrage du modèle CIR ###################################################
#On considere param comme une liste/matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 6 classes de rating et 4 parametres pour chacune d entre elle,
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, 
#                      muAAA, muAA, muA, muBaa, muBa, muB, 
#                      sigma, 
#                      pi0AAA, pi0AA, pi0A, pi0BBB, pi0BB, pi0B)

param_init <- c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003,
                0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000,
                0.01,
                0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055)

ecart_JLT <- function(param){
  e <- 0
  for (t in 1:TT){ for (i in 1:6){
    # spread <- spread_i_fct(1000, t, TT, param_test, M, D, i,LGD)
    spread <- spread_i_JLT_calibrage(t, param, M, D, i, LGD)
    e <- e + (spread-SpreadMarket[i])^2
  }}
  return(e)
} ## parametre sur le rating

TT <- 1 # les données fornis sont les spreads de maturité 1 an
LB <- c(rep(0,8),rep(0,8),0,rep(0,8))
UB <- c(rep(30,8),rep(20,8),10,rep(30,8))
#paramJLT = hjkb(param_init,Ecart_JLT,lower=LB,upper=UB)$par # 2.731286e-05
paramJLT = nlminb(start = param_init,ecart_JLT,lower=LB,upper=UB)$par
#paramJLT = optim(par=param_init,fn=Ecart_JLT,lower=LB,upper=UB,method="L-BFGS-B")$par
(listJLT <- list(
  k = paramJLT[1:8],
  mu = paramJLT[9:16],
  sigma = paramJLT[17],
  pi0 = paramJLT[18:25]))
(paramJLT[17] = min(sqrt(2*listJLT$k*listJLT$mu)[sqrt(2*listJLT$k*listJLT$mu) > 0]))
# il faut que sigma soit non nul
listJLT





#### Some plot JLT ###################################################
N = 1000

# Plot spread
spreadAAA <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 1, LGD)
spreadAA <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 2, LGD)
spreadA <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 3, LGD)
spreadBBB <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 4, LGD)
spreadBB <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 5, LGD)
spreadB <- spread_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 6, LGD)

plot(spreadAAA, 
     main='Spread de crédit reproduits par le modèle JLT', type='l', 
     ylab="Spread", 
     ylim=c(0,0.04), xlab="Maturité", col='red')
lines(spreadAA, col='orange')
lines(spreadA, col='brown')
lines(spreadBBB, col='lightblue')
lines(spreadBB, col='blue')
lines(spreadB, col='purple')
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)

#Essai avec proba de défaut calibrage :
probaDefAAA <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 1)
probaDefAA <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 2)
probaDefA <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 3)
probaDefBBB <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 4)
probaDefBB <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 5)
probaDefB <- proba_defaut_i_JLT_calibrage.T(Maturite30, paramJLT, M, D, 6)

plot(1-probaDefAAA, 
     main='Probabilité de Survie selon la notation', 
     ylim=c(0,1), type='l', 
     ylab="Probabilité de Survie", xlab="Maturité", col='red')
lines(1-probaDefAA, col='orange')
lines(1-probaDefA, col='brown')
lines(1-probaDefBBB, col='lightblue')
lines(1-probaDefBB, col='blue')
lines(1-probaDefB, col='purple')
legend("bottomleft",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)

#Essai avec proba de défaut normal :
probaDefAAAnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 1))
probaDefAAnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 2))
probaDefAnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 3))
probaDefBBBnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 4))
probaDefBBnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 5))
probaDefBnorm <- colMeans(proba_defaut_i_JLT.T(1000, 0, Maturite30, paramJLT, M, D, 6))

plot(probaDefAAAnorm, 
     main='Probabilité de défaut selon la notation', 
     ylim=c(0,1), type='l', 
     ylab="Probabilité de défaut", xlab="Maturité", col='red')
lines(probaDefAAnorm, col='orange')
lines(probaDefAnorm, col='brown')
lines(probaDefBBBnorm, col='lightblue')
lines(probaDefBBnorm, col='blue')
lines(probaDefBnorm, col='purple')
legend("topleft",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.7)

#simulation de spread AA pour une maturité de 10 ans
tt <- seq(0,10,0.1)[-length(seq(0,10,0.1))]
SP <- spread_i_JLT.t(N, tt, 10, paramJLT, M, D, 2, LGD)
matplot(tt,t(SP[1:10,]),type="l",
        main="Simulation spread sur maturité 10 ans pour AA",
        ylab="Spread", 
        xlab="Maturité")
plot(tt,colMeans(SP),'l',col="red",lwd=1.5,
     main=c("Moyenne de spread actions AA maturité 10 ans", "(sur 1000 simulations)"),
     xlab="Maturité",
     ylab="Spread")





#### Test de martingalité sur Cash Flow JLT ###################################################
# prix risqué à t=0, calculé à partir du spread FF et PZC FF
PZCr_i_CF_JLT_FF <- function(TT, param_Vas, param_JLT, M, D, i, LGD){
  return(PZC_Vas_FF_calibrage.T(TT, param_Vas)/(1+spread_i_JLT_calibrage.T(TT, param_JLT, M, D, i, LGD))^TT)
}

# moyenne prix ZC risqué CF actualisé au taux sans risque
PZCr_i_CF_JLT_sim <- function(N, t, TT, param_Vas, param_JLT, M, D, i, LGD){
  pd_i <- proba_defaut_i_JLT.T(N, t, TT, param_JLT, M, D, i)
  return(exp(-TT*TZC_Vas_FF_calibrage(TT, param_Vas))*t(1*(1-pd_i)+(1-LGD)*pd_i))
}

plot(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 1, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 1, LGD)),
     ylim=c(0.99,1.01),"l",col="red",
     ylab="CashFlow Actualisé",
     main="Test Martingalité par CashFlow du modèle JLT")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 2, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 3, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 4, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 5, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite,PZCr_i_CF_JLT_FF(Maturite, paramVas, paramJLT, M, D, 6, LGD)/rowMeans(PZCr_i_CF_JLT_sim(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)





#### Test de martingalité sur PZCr avec rating JLT ###################################################
# adaptation pour JLT
PZC_Vas_sim <- function(N, t, TT, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b * (TT - t)) * exp(-(t(taux_Vas_sim.t(N, TT, param, r0)) - b) * (1 - exp(-a * (TT - t))) / a + 0.5 * (sigma^2 * (TT - t) / a^2 - sigma^2 / a^3 * (1 - exp(-a * (TT - t))) - sigma^2 / (2 * a^3) * (1 - exp(-a * (TT - t)))^2)))
}

# moyenne prix ZC risqué actualisé au taux sans risque
PZCr_i_JLT_sim <- function(N, t, TT, param_Vas, param_JLT, M, D, i, LGD){
  sum = proba_defaut_i_JLT_calibrage(t, param_JLT, M, D, i)*(1-LGD)*PZC_Vas_sim(N, t, TT, param_Vas)
  for (l in 1:(8-1)){
    p_il = proba_passage_ij_JLT(N, 0, t, param_JLT, M, D, i, l)
    PZCr_l = PZCr_i_CF_JLT_sim(N, t, TT, param_Vas, param_JLT, M, D, l, LGD)
    sum = sum + p_il*PZCr_l
  }
  return(exp(-t*TZC_Vas_FF_calibrage(t, param_Vas))*sum)
}

t = 0.1
plot(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 1, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 1, LGD)),
     ylim=c(0.9,1.1),"l",col="red",
     xlab ="Maturité",
     ylab="PZC Actualisé",
     main="Test Martingalité sur PZC risqué du modèle JLT")
lines(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 2, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 3, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 4, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 5, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite30,PZCr_i_CF_JLT_FF(Maturite30, paramVas, paramJLT, M, D, 6, LGD)/rowMeans(PZCr_i_JLT_sim(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)





### END CALIBRATION ET VALIDATION #############################################################
