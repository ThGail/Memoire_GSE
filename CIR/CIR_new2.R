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
  
  return(matrix(A*exp(- lambda_t %*% t(B)), nrow = N))
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

plot(1:15,Prix_Spread_fct(0, 1:15, param_AAA, LGD), 
     main='Spread reproduits par le modèle CIR (T = 15)', 
     type='l', ylab="Spread", xlab="Maturité", 
     ylim=c(0.0005, 0.05), col='red')
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

#survie sim :
survie_FF_sim_AAA <- c()
survie_FF_sim_AA <- c()
survie_FF_sim_A <- c()
survie_FF_sim_BBB <- c()
survie_FF_sim_BB <- c()
survie_FF_sim_B <- c()

#Survie  normale
survie_FF_AAA <- c()
survie_FF_AA <- c()
survie_FF_A <- c()
survie_FF_BBB <- c()
survie_FF_BB <- c()
survie_FF_B <- c()

N <- 1000

for (t in 0:150){
  #Survie sim : 
  survie_FF_sim_AAA <- c(survie_FF_sim_AAA, mean(survie_FF_fct_sim(N, 0, t, param_AAA)))
  survie_FF_sim_AA <- c(survie_FF_sim_AA, mean(survie_FF_fct_sim(N, 0, t, param_AA)))
  survie_FF_sim_A <- c(survie_FF_sim_A, mean(survie_FF_fct_sim(N, 0, t, param_A)))
  survie_FF_sim_BBB <- c(survie_FF_sim_BBB, mean(survie_FF_fct_sim(N, 0, t, param_BBB)))
  survie_FF_sim_BB <- c(survie_FF_sim_BB, mean(survie_FF_fct_sim(N, 0, t, param_BB)))
  survie_FF_sim_B <- c(survie_FF_sim_B, mean(survie_FF_fct_sim(N, 0, t, param_B)))
  
  
  #Survie normale : 
  survie_FF_AAA <- c(survie_FF_AAA, survie_FF_fct(N, 0, t, param_AAA))
  survie_FF_AA <- c(survie_FF_AA, survie_FF_fct(N, 0, t, param_AA))
  survie_FF_A <- c(survie_FF_A, survie_FF_fct(N, 0, t, param_A))
  survie_FF_BBB <- c(survie_FF_BBB, survie_FF_fct(N, 0, t, param_BBB))
  survie_FF_BB <- c(survie_FF_BB, survie_FF_fct(N, 0, t, param_BB))
  survie_FF_B <- c(survie_FF_B, survie_FF_fct(N, 0, t, param_B))
}

#Plot avec Sim
plot(survie_FF_sim_AAA,col="red",type="l",
     main = 'Probabilité de survie selon la notation',
     ylim = c(0, 1),
     xlab="Maturité",
     ylab="Probilité de survie")
lines(survie_FF_sim_AA, col = 'orange')
lines(survie_FF_sim_A, col = 'brown')
lines(survie_FF_sim_BBB, col = 'lightblue')
lines(survie_FF_sim_BB , col = 'blue')
lines(survie_FF_sim_B, col = 'purple')
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)

#Plot avec normal :
plot(survie_FF_AAA,col="red",type="l",
     main = 'Probabilité de survie selon la notation',
     ylim = c(0, 1),
     xlab="Maturité",
     ylab="Probabilité de survie ")
lines(survie_FF_AA, col = 'orange')
lines(survie_FF_A, col = 'brown')
lines(survie_FF_BBB, col = 'lightblue')
lines(survie_FF_BB , col = 'blue')
lines(survie_FF_B, col = 'purple')
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)






########## MARTINGALITE ##########

#### Les fonctions utiles du Vasicek
# paramVas <- c(0.0095787241,0.0382841814,0.0002632241)
paramVas <- c(0.004084473,0.292354584,0.003962402) # calibrer à partir de TZC_Vas_FF_cali
paramVas2 <- c(0.008738708, 0.041931152, 0.001186371)
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
plot(TauxZC, main='courbe zéro coupon EIOPA')
lines(TZC_Vas_FF_calibrage(Maturite,paramVas), col ='red')
lines(TZC_Vas_FF_calibrage(Maturite,paramVas2), col= 'blue')

# resultat sur PZC
plot(exp(-Maturite*TauxZC))
lines(exp(-Maturite*TZC_Vas_FF_calibrage(Maturite,paramVas)), col = 'red')

# test pour voir erreur de calibrage du Vas
plot(PZC_Vas_FF_calibrage(Maturite,paramVas)/exp(-Maturite*TauxZC))
abline(h=1)
abline(h=1.05)
abline(h=0.95)


########## Test de martingalité : CASH FLOW ##########
paramCIR <- list(AAA = param_AAA,
                 AA = param_AA, 
                 A = param_A,
                 BBB = param_BBB, 
                 BB = param_BB, 
                 B = param_B)
paramCIR$AAA

# prix risqué à t=0, calculé à partir du spread FF et PZC FF
PZCr_i_CF_CIR_FF <- function(TT, param_taux, param_CIR, LGD){
  return(PZC_Vas_FF_calibrage.T(TT, param_taux)/(1+Prix_Spread_fct(0,TT,param_CIR,LGD))^TT)
}
plot(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$AAA, LGD), main='Prix ZC risqué par Formule Fermée', ylab = 'Prix')
lines(exp(-Maturite*TauxZC))

#survie_FF_fct_sim.T <- Vectorize(survie_FF_fct_sim,"TT")
# moyenne prix ZC risqué actualisé au taux sans risque




PZCr_i_CF_CIR_sim <- function(N, TT, param_taux, param_CIR, LGD){
  surv_i <- survie_FF_fct_sim(N, 0.1, TT, param_CIR)
  #return(exp(-TT*TauxZC[TT])*t(1*(surv_i)+(1-LGD)*(1-surv_i))
  return(exp(-TT*TZC_Vas_FF_calibrage(TT, param_taux))*t(surv_i+(1-LGD)*(1-surv_i)))
}

#
plot(Maturite, rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$AAA, LGD)), main='Prix empirique ZC risqué', ylab = 'Prix')

#Comparaison des deux courbes :
plot(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$AAA, LGD), main='Prix ZC risqué par Formule Fermée', ylab = 'Prix', col = 'black', type = 'l')
lines(Maturite, rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$AAA, LGD)), main='Prix empirique ZC risqué', ylab = 'Prix', col = 'red')
lines(exp(-Maturite*TauxZC), col = 'blue')
legend("topright", legend=c("Prix ZCR FF","Prix ZC Empirique", "Prix ZC EIOPA"),
       col=c("black","red", "blue"),pch=20,
       cex=0.8)

# Les plots :
plot(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$AAA, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$AAA, LGD)),"l",
     ylab="CashFlow Actualisé",
     main="Test de martingalité CIR", 
     col = 'red')
lines(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$AA, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$AA, LGD)),"l",col="orange")
lines(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$A, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$A, LGD)),"l",col="brown")
lines(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$BBB, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$BBB, LGD)),"l",col="lightblue")
lines(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$BB, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$BB, LGD)),"l",col="blue")
lines(Maturite,PZCr_i_CF_CIR_FF(Maturite, paramVas, paramCIR$B, LGD)/rowMeans(PZCr_i_CF_CIR_sim(N, Maturite, paramVas, paramCIR$B, LGD)),"l",col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)
