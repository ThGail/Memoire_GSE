rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\HW_calibrage.R",sep=""))
source(paste(chemin_cal,"\\BS_calibrage.R",sep=""))
source(paste(chemin_cal,"\\CIR_calibrage.R",sep=""))
source(paste(chemin_cal,"\\JLT_calibrage.R",sep=""))
source(paste(chemin_fct,"\\GSE_fonction.R",sep=""))

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
Maturite50 <- Maturite[1:50]

### Parametres des modèles #####################################
# Vasicek
paramVas
# HW
paramHW
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

N = 1000; TT = 50; LGD = 0.3





#### Some plot GSE Vas CIR ########################################
Tvas.test <- c()
action.test <- c()
immo.test <- c()
PZCrAAAt1.test <- c()
PZCrAAAt1sc.test <- c()
PZCrAAAT50.test <- c()

for (t in (1:TT)){
  GSET50.test <- GSE_Vas_CIR(N, t, t+10, paramVas, paramBS$action, paramBS$immo, paramCIR, L, LGD)
  GSEt1.test <- GSE_Vas_CIR(N, 1, t, paramVas, paramBS$action, paramBS$immo, paramCIR, L, LGD)
  Tvas.test <- cbind(Tvas.test, GSET50.test$TZC_Vas)
  action.test <- cbind(action.test, GSET50.test$ind_action)
  immo.test <- cbind(immo.test, GSET50.test$ind_immo)
  PZCrAAAT50.test <- cbind(PZCrAAAT50.test, GSET50.test$PZCr_AAA)
  PZCrAAAt1.test <- cbind(PZCrAAAt1.test, GSEt1.test$PZCr_AAA)
  PZCrAAAt1sc.test <- c(PZCrAAAt1sc.test, mean(PZCr_CIR_Vas_sim(N, 1, t, paramVas, paramCIR$AAA, LGD)))
}

# plot du taux sans risque Vasicek avec matrice de correlation
matplot(t(Tvas.test[1:10,]), type='l', main="Simulation du taux zéro-coupon Vasicek", ylab='Taux ZC', xlab='Temps')
lines(colMeans(Tvas.test),  type = "l", col="red",lwd=2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot de l'indice action :
matplot(t(action.test[1:20,]), type='l',  main="Evolution de la performance action sur 50 ans", ylab='Indice action', xlab='Temps')
lines(colMeans(action.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot de l'indice immobilier :
matplot(t(immo.test[1:20,]), type='l', main="Evolution de la performance immobilière sur 50 ans", ylab='Indice immobilier', xlab='Temps')
lines(colMeans(immo.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot PZCr AAA de maturite 10 ans
matplot(t(PZCrAAAT50.test[1:10,]), type='l',
        main=c("Evolution PZCr AAA de maturité 10 ans","Vasicek CIR"),xlab="Temps",ylab="PZCr")
lines(colMeans(PZCrAAAT50.test),  type = "l", col="red",lwd=2)

# plot PZCr AAA au temps t=1
matplot(t(PZCrAAAt1.test[1:10,]), type='l',
        main=c("Projection PZCr AAA au temps t=1","Vasicek CIR"),xlab="Maturite",ylab="PZCr")
lines(colMeans(PZCrAAAt1.test),  type = "l", col="red",lwd=2)
lines(PZCrAAAt1sc.test)
legend("bottomleft", legend=c("PZC risqué AAA avec corrélation des taux (issu GSE)","PZC risqué AAA sans corrélation des taux"),
       col=c("red","black"),pch=20,
       cex=0.7)





#### Some plot GSE Vas JLT ########################################
PZCrAAAt1.test <- c()
PZCrAAAt1sc.test <- c()
PZCrAAAT50.test <- c()

for (t in (1:TT)){
  GSET50.test <- GSE_Vas_JLT(N, t, t+10, paramVas, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  GSEt1.test <- GSE_Vas_JLT(N, 1, t, paramVas, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  PZCrAAAT50.test <- cbind(PZCrAAAT50.test, GSET50.test$PZCr_AAA)
  PZCrAAAt1.test <- cbind(PZCrAAAt1.test, GSEt1.test$PZCr_AAA)
  PZCrAAAt1sc.test <- c(PZCrAAAt1sc.test, mean(PZCr_i_JLT_Vas_sim(N, 1, t, paramVas, paramJLT, M, D, 1, LGD)))
}

# plot PZCr AAA de maturite 10 ans
matplot(t(PZCrAAAT50.test[1:10,]), type='l',
        main=c("Evolution PZCr AAA de maturité 10 ans","Vasicek JLT"),xlab="Temps",ylab="PZCr")
lines(colMeans(PZCrAAAT50.test),  type = "l", col="red",lwd=2)

# plot PZCr AAA au temps t=1
matplot(t(PZCrAAAt1.test[1:10,]), type='l',
        main=c("Projection PZCr AAA au temps t=1","Vasicek JLT"),xlab="Maturite",ylab="PZCr")
lines(colMeans(PZCrAAAt1.test),  type = "l", col="red",lwd=2)
lines(PZCrAAAt1sc.test)
legend("bottomleft", legend=c("PZC risqué AAA avec corrélation des taux (issu GSE)","PZC risqué AAA sans corrélation des taux"),
       col=c("red","black"),pch=20,
       cex=0.7)






#### Some plot GSE HW CIR  ########################################
Thw.test <- c()
action.test <- c()
immo.test <- c()
PZCrAAAt1.test <- c()
PZCrAAAt1sc.test <- c()
PZCrAAAT50.test <- c()

for (t in (1:TT)){
  GSET50.test <- GSE_HW_CIR(N, t, t+10, paramHW, paramBS$action, paramBS$immo, paramCIR, L, LGD)
  GSEt1.test <- GSE_HW_CIR(N, 1, t, paramHW, paramBS$action, paramBS$immo, paramCIR, L, LGD)
  Thw.test <- cbind(Thw.test, GSET50.test$TZC_HW)
  action.test <- cbind(action.test, GSET50.test$ind_action)
  immo.test <- cbind(immo.test, GSET50.test$ind_immo)
  PZCrAAAT50.test <- cbind(PZCrAAAT50.test, GSET50.test$PZCr_AAA)
  PZCrAAAt1.test <- cbind(PZCrAAAt1.test, GSEt1.test$PZCr_AAA)
  PZCrAAAt1sc.test <- c(PZCrAAAt1sc.test, mean(PZCr_CIR_HW_sim(N, 1, t, paramHW, paramCIR$AAA, LGD)))
}

# plot du taux sans risque Hull&White avec matrice de correlation
matplot(t(Thw.test[1:10,]), type='l', main="Simulation du taux zéro-coupon Hull&White", ylab='Taux ZC', xlab='Maturité')
lines(colMeans(Thw.test),  type = "l", col="red",lwd=2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot de l'indice action :
matplot(t(action.test[1:20,]), ylim=c(0,150), type='l',  main="Evolution de la performance action sur 50 ans", ylab='Indice action', xlab='Temps')
lines(colMeans(action.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot de l'indice immobilier :
matplot(t(immo.test[1:20,]), ylim=c(0,150), type='l', main="Evolution de la performance immobilière sur 50 ans", ylab='Indice immobilier', xlab='Temps')
lines(colMeans(immo.test), type = 'l', col="red", lwd = 2)
legend("topleft", legend=c("Moyenne des simulations"),
       col=c("red"),pch=20,
       cex=0.7)

# plot PZCr AAA de maturite 10 ans
matplot(t(PZCrAAAT50.test[1:10,]), type='l',
        main=c("Evolution PZCr AAA de maturité 10 ans","HW CIR"),xlab="Temps",ylab="PZCr")
lines(colMeans(PZCrAAAT50.test),  type = "l", col="red",lwd=2)

# plot PZCr AAA au temps t=1
matplot(t(PZCrAAAt1.test[1:10,]), type='l',
        main=c("Projection PZCr AAA au temps t=1","HW CIR"),xlab="Maturite",ylab="PZCr")
lines(colMeans(PZCrAAAt1.test),  type = "l", col="red",lwd=2)
lines(PZCrAAAt1sc.test)
legend("bottomleft", legend=c("PZC risqué AAA avec corrélation des taux (issu GSE)","PZC risqué AAA sans corrélation des taux"),
       col=c("red","black"),pch=20,
       cex=0.7)





#### Some plot GSE HW JLT ########################################
PZCrAAAt1.test <- c()
PZCrAAAt1sc.test <- c()
PZCrAAAT50.test <- c()

for (t in (1:TT)){
  GSET50.test <- GSE_HW_JLT(N, t, t+10, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  GSEt1.test <- GSE_HW_JLT(N, 1, t, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  PZCrAAAT50.test <- cbind(PZCrAAAT50.test, GSET50.test$PZCr_AAA)
  PZCrAAAt1.test <- cbind(PZCrAAAt1.test, GSEt1.test$PZCr_AAA)
  PZCrAAAt1sc.test <- c(PZCrAAAt1sc.test, mean(PZCr_i_JLT_HW_sim(N, 1, t, paramHW, paramJLT, M, D, 1, LGD)))
}

# plot PZCr AAA de maturite 10 ans
matplot(t(PZCrAAAT50.test[1:10,]), type='l',
        main=c("Evolution PZCr AAA de maturité 10 ans","HW JLT"),xlab="Temps",ylab="PZCr")
lines(colMeans(PZCrAAAT50.test),  type = "l", col="red",lwd=2)

# plot PZCr AAA au temps t=1
matplot(t(PZCrAAAt1.test[1:10,]), type='l',
        main=c("Projection PZCr AAA au temps t=1","HW JLT"),xlab="Maturite",ylab="PZCr")
lines(colMeans(PZCrAAAt1.test),  type = "l", col="red",lwd=2)
lines(PZCrAAAt1sc.test)
legend("bottomleft", legend=c("PZC risqué AAA avec corrélation des taux (issu GSE)","PZC risqué AAA sans corrélation des taux"),
       col=c("red","black"),pch=20,
       cex=0.7)



