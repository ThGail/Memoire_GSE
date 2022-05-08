rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))
source(paste(chemin_fct,"\\BS_fonction.R",sep=""))
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\BS_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]
Maturite50 <- Maturite[1:50]

N = 1000; S0 = 1; K = 1

### Plot du Black and Scholes ##############################################
# moyenne simulation indice action pour 1000 simulations
(plotPAaction = PA_BS_Vas_FF_sim.t(N,Maturite30,S0,paramVas,paramBS$action))
matplot(t(plotPAaction[1:10,]),type="l",
        main="Scénarios de l'indice action sur 30 ans",
        xlab="Temps",ylab="indice action")
lines(colMeans(plotPAaction),type="l",lwd=2,col="red")

# moyenne simulation indice immo pour 1000 simulations
(plotPAimmo = PA_BS_Vas_FF_sim.t(N,Maturite30,S0,paramVas,paramBS$immo))
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
