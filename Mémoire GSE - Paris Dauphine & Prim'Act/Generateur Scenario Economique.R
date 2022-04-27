rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Fonction"
chemin_cal = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Calibrage"
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\BS_calibrage.R",sep=""))
source(paste(chemin_cal,"\\CIR_calibrage.R",sep=""))
source(paste(chemin_cal,"\\JLT_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
PrixZC <- exp(-Maturite*TauxZC)
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]

### Parametres des modèles #####################################
# Vasicek
paramVas
# BS
paramBS
# CIR
paramCIR
# JLT
paramJLT

# matrice de correlaion entre taux_court, action_glo et immo
Corr <- matrix(c(1,0.1,0.05,0.1,1,-0.3,0.05,-0.3,1),nrow=3)
# décompostion de Cholesky :
L <- t(chol(Corr, pivot = FALSE))

### Définition fonctions utiles ###################################
# Generateur de mouvement brownien N * n_type
# n_type est le nombre de donnée intéret
# ex : TZC pour modèle de taux, spread pour modèle de crédit, etc
# MvtBrow <- function(N,n_type){
#   return(matrix(rnorm(N * n_type), nrow = N))
# }

# la fonction qui permet de simuler les taux instantanées de Vasicek
Vas_W <- function(t, param, W, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(r0 * exp(-a * t) + b * (1 - exp(-a * t)) + sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a)) * W)
}

# transformation du taux instantanée en TZC de Vasicek
tauxToTZC_Vas <- function(t,TT,taux_t,param){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  theta <- TT-t
  return(Ri - ((Ri-taux_t)*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}

# la fonction qui simule le prix des actifs du BS
BS_W <- function(t, param_Vas, param_actif, W_vas, W_act, S0=1, r0=TauxZC[1]) {
  return(S0 * exp((Vas_W(t, param_Vas, W_vas, r0) - 0.5 * param_actif^2) * t + param_actif * sqrt(t) * W_act))
}

### GSE temoin en utilisant Vasicek pour le generateur de taux sans risque ##############
GSE_Vas_temoin <- function(N, t, TT, param_Vas, param_action, param_immo, L, S0=1,r0=TauxZC[1]){

  # simulation de gaussiennes corrélées
  # on simule 3 car : 1 taux sans risque + 2 actifs + 6 crédits
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.vas <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  # calcul des grandeurs de base
  vas <- Vas_W(t, param_Vas, W.vas, r0)
  action <- BS_W(t, param_Vas, param_action, W.vas, W.action,S0,r0)
  immo <- BS_W(t, param_Vas, param_immo, W.vas, W.immo,S0,r0)
  
  TZCvas <- tauxToTZC_Vas(t,TT,vas,param_Vas)
  
  return(list(TZC_Vas=TZCvas,ind_action=action,ind_immo=immo))
}

# Some plot GSE temoin avec Vasicek dans le cas maturite T=50
N = 1000; TT=50
Tvas.test <- c()
action.test <- c()
immo.test <- c()

for (t in (1:TT)){
  GSE.test <- GSE_Vas_temoin(N, t, TT, paramVas, paramBS$action, paramBS$immo, L)
  Tvas.test <- cbind(Tvas.test, GSE.test$TZC_Vas)
  action.test <- cbind(action.test, GSE.test$ind_action)
  immo.test <- cbind(immo.test, GSE.test$ind_immo)
}
# plot du taux sans risque Vasicek avec matrice de correlation
matplot(t(Tvas.test[1:10,]), type='l', main="Simulation du taux zéro-coupon Vasicek", ylab='Taux ZC', xlab='Maturité')
lines(colMeans(Tvas.test),  type = "l", col="red",lwd=2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)
# plot de l'indice action :
matplot(t(action.test[1:10,]), type='l',  main="Simulation de l'indice action", ylab='Prix Action', xlab='Maturité')
lines(colMeans(action.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)
# plot de l'indice immobilier :
matplot(t(immo.test[1:10,]), type='l', main="Simulation de l'indice immobilier", ylab='Prix Immobilier', xlab='Maturité')
lines(colMeans(immo.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

### GSE complet avec risque de crédit en utilisant Vasicek pour le generateur de taux sans risque ####
GSE_Vas_JLT_complet <- function(N, t, TT, param_Vas, param_action, param_immo, param_JLT, L, LGD, S0=1,r0=TauxZC[1]){
  
  # simulation de gaussiennes corrélées
  # on simule 3 car : 1 taux sans risque + 2 actifs + 6 crédits
  Z <- matrix(rnorm(N*3), nrow = N)
  W <- L %*% t(Z)
  W.vas <- W[1,]
  W.action <- W[2,]
  W.immo <- W[3,]
  
  # calcul des grandeurs de base
  vas <- Vas_W(t, param_Vas, W.vas, r0)
  action <- BS_W(t, param_Vas, param_action, W.vas, W.action,S0,r0)
  immo <- BS_W(t, param_Vas, param_immo, W.vas, W.immo,S0,r0)
  
  TZCvas <- tauxToTZC_Vas(t,TT,vas,param_Vas)
  PZCvas <- exp(-(TT-t)*TZCvas)
  
  # calcul des spreads de crédit / proba de survie / proba de défaut / PZCr ?
  # on simule 6 autres pour les chacunes de rating (transformé en proba de survie)
  # les crédits sont indépendants les uns des autres (hypothses)
  ps.AAA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 1)
  ps.AA <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 2)
  ps.A <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 3)
  ps.BBB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 4)
  ps.BB <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 5)
  ps.B <- proba_survie_i_JLT(N, t, TT, param_JLT, M, D, 6)
  
  # calcul du PZR risqué
  pzcr.AAA <- PZCvas * (1-LGD+LGD*ps.AAA)
  pzcr.AA <- PZCvas * (1-LGD+LGD*ps.AA)
  pzcr.A <- PZCvas * (1-LGD+LGD*ps.A)
  pzcr.BBB <- PZCvas * (1-LGD+LGD*ps.BBB)
  pzcr.BB <- PZCvas * (1-LGD+LGD*ps.BB)
  pzcr.B <- PZCvas * (1-LGD+LGD*ps.B)
  
  return(list(TZC_Vas=TZCvas,ind_action=action,ind_immo=immo,
              ps_AAA=ps.AAA,ps_AA=ps.AA,ps_A=ps.A,ps_BBB=ps.BBB,ps_BB=ps.BB,ps_B=ps.B,
              PZCr_AAA=pzcr.AAA,PZCr_AA=pzcr.AA,PZCr_A=pzcr.A,PZCr_BBB=pzcr.BBB,PZCr_BB=pzcr.BB,PZCr_B=pzcr.B))
}

# Some plot GSE complet avec Vasicek
N = 1000; TT=50
PZCrAAAt1.test <- c()
PZCrAAAt1sc.test <- c()
PZCrAAAT50.test <- c()
for (t in (1:TT)){
  GSET50.test <- GSE_Vas_JLT_complet(N, t, TT, paramVas, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  GSEt1.test <- GSE_Vas_JLT_complet(N, 1, t, paramVas, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  PZCrAAAT50.test <- cbind(PZCrAAAT50.test, GSET50.test$PZCr_AAA)
  PZCrAAAt1.test <- cbind(PZCrAAAt1.test, GSEt1.test$PZCr_AAA)
  PZCrAAAt1sc.test <- c(PZCrAAAt1sc.test, mean(PZCr_i_JLT_sim(N, 1, t, paramVas, paramJLT, M, D, 1, LGD)))
}
matplot(t(PZCrAAAt1.test[1:10,]), type='l',
        main="Projection PZCr AAA au temps t=1",xlab="Maturite",ylab="PZCr")
lines(colMeans(PZCrAAAt1.test),  type = "l", col="red",lwd=2)
lines(PZCrAAAt1sc.test)
legend("bottomleft", legend=c("PZC risqué AAA avec corrélation des taux (issu GSE)","PZC risqué AAA sans corrélation des taux"),
       col=c("red","black"),pch=20,
       cex=0.7)

matplot(t(PZCrAAAT50.test[1:10,]), type='l',
        main="Evolution PZCr AAA de maturité 50 ans",xlab="Temps",ylab="PZCr")
lines(colMeans(PZCrAAAT50.test),  type = "l", col="red",lwd=2)


