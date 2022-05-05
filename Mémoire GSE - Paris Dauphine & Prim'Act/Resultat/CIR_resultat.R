rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))
source(paste(chemin_fct,"\\CIR_fonction.R",sep=""))
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\CIR_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]

N = 1000

### Plot CIR ####################################################
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

### (Test de martingalité CIR) ###################################################
plot(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$AAA, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$AAA, LGD)),"l",
     ylab="CashFlow Actualisé",ylim=c(0.99,1.01),
     main="Test de martingalité CIR", 
     col = 'red')
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$AA, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$AA, LGD)),"l",col="orange")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$A, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$A, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$BBB, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$BBB, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$BB, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$BB, LGD)),"l",col="brown")
lines(Maturite,PZCr_CIR_FF(Maturite, paramVas, paramCIR$B, LGD)/colMeans(PZCr_CIR_Vas_sim.T(N, 0, Maturite, paramVas, paramCIR$B, LGD)),"l",col="brown")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.8)
