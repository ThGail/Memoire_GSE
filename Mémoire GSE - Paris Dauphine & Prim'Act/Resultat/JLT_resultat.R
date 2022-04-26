rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Fonction"
chemin_cal = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Calibrage"
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))
source(paste(chemin_fct,"\\JLT_fonction.R",sep=""))
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\JLT_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]

N = 1000; LGD = 0.3

### Plot JLT ########################################################
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

### Test de martingalité sur Cash Flow JLT ###################################################
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

### Test de martingalité sur PZCr avec rating JLT ###################################################
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

