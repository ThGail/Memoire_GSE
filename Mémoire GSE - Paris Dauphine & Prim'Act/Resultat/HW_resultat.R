rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\HW_fonction.R",sep=""))
source(paste(chemin_cal,"\\HW_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]
Maturite50 <- Maturite[1:50]

N = 1000

### Plot du HW ###################################################
# test pour le PZC
plot(PZC_HW_sim(N=1,0,Maturite,paramHW),type="l",lwd=2,
     main = 'Vérification avec le prix zéro-coupon', 
     xlab = 'Maturité', ylab = 'Prix zéro-coupon')
lines(PrixZC,col="red",lty=3,lwd=2)
# pour t=0, PZC est fermée, il suffit de faire "1" simulation
legend("topright",legend=c("PZC donnée (input)","PZC du modèle HW"),
       col=c("red","black"),pch=20,
       cex=0.6)

# test pour le TZC
plot(TZC_HW_sim(N=1,0,Maturite,paramHW),type="l",lwd=2,
     main = 'Vérification avec le taux zéro-coupon', 
     xlab = 'Maturité', ylab = 'Taux ZC')
lines(TauxZC,col="red",lty=3,lwd=2)
legend("bottomright",legend=c("TZC donnée (input)","TZC du modèle HW"),
       col=c("red","black"),pch=20,
       cex=0.6)

# simulation PZC dans 1 an pour différentes maturités
PZCt1 <- PZC_HW_sim.T(N,1,Maturite,paramHW)
matplot(Maturite,t(PZCt1[1:10,]),type='l',
        main = 'Prix zéro-coupon dans 1 an simulé par HW',
        xlab = 'Maturité', ylab = 'Prix ZC')
lines(colMeans(PZCt1),type="l",lwd=2,col="red")
lines(PrixZC)
legend("topright",legend=c("PZC en t=1","PZC en t=0"),
       col=c("red","black"),pch=20,
       cex=0.6)

# trajectoire du PZC de maturité 10 ans
PZCT10 <- c()
for (t in c(0,Maturite)){PZCT10 <- cbind(PZCT10,PZC_HW_sim(N,t,t+10,paramHW))}
matplot(t(PZCT10[1:10,]),type='l',
        main = 'Trajectoire prix zéro-coupon de maturité 10 ans',
        xlab = 'Temps', ylab = 'Prix ZC')
lines(colMeans(PZCT10),type="l",lwd=2,col="red")



