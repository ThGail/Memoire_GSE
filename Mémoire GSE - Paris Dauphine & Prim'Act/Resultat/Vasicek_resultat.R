rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]

N = 1000

### Plot du Vasicek ###################################################
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

# simulation PZC de maturité 10 ans (1000 simulations)
plotPZC <- PZC_Vas_FF_sim.t(N, tt, 10, paramVas)
matplot(tt,t(plotPZC[1:10,]),type="l",
        main="Scénarios du prix zéro-coupon de maturité 10 ans",
        xlab="Temps",ylab="PZC")
plot(tt,colMeans(plotPZC),type="l",lwd=2,col="red",
     main="Moyenne du prix zéro-coupon de maturité 10 ans",
     xlab="Temps",ylab="TZC")

# projection à 5 ans du PZC pour différente maturité
PZC5 = PZC_Vas_FF_sim(1000, 5, Maturite, paramVas)
matplot(t(PZC5[1:10,]),type='l',
        main = "Projection du PZC dans 5 ans",xlab="Maturite",ylab="PZC")
lines(PZC_Vas_FF_calibrage(Maturite,paramVas),lwd=2,lty=20)
lines(colMeans(PZC5),lwd=2,col="red")
legend("topright",legend=c("Courbe PZC du modèle","Moyenne PZC avec t=5"),
       col=c("black","red"),pch=20,
       cex=0.6)

