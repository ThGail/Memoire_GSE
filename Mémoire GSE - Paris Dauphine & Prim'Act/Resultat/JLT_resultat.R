rm(list = ls())
library(readxl)
library(dfoptim)
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
chemin_cal = paste(chemin,"/Calibrage", sep="")
setwd(chemin)

### Importation des donnees ###################################################
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))
source(paste(chemin_fct,"\\HW_fonction.R",sep=""))
source(paste(chemin_fct,"\\JLT_fonction.R",sep=""))
source(paste(chemin_cal,"\\Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"\\HW_calibrage.R",sep=""))
source(paste(chemin_cal,"\\JLT_calibrage.R",sep=""))

sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])
Maturite10 <- Maturite[1:10]
Maturite15 <- Maturite[1:15]
Maturite30 <- Maturite[1:30]
Maturite50 <- Maturite[1:50]

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

# simulation de spread B pour une maturité de 3 ans
SP <- c()
for (t in c(0,Maturite50)){SP <- cbind(SP,spread_i_JLT(N, t, t+3, paramJLT, M, D, 6, LGD))}
matplot(t(SP[1:10,]),type="l",
        main=c("Evolution de spread B de maturité 3 ans sur 50 ans","JLT"),
        ylab="Spread",xlab="Temps")
lines(colMeans(SP),col="red",lwd=1.5)

# Simulation spread AA de maturité 10 ans
SP <- c()
for (t in c(0,Maturite50)){SP <- cbind(SP,spread_i_JLT(N, t, t+10, paramJLT, M, D, 2, LGD))}
matplot(t(SP[1:10,]),type="l",
        main="Simulation spread sur maturité 10 ans pour AA",
        ylab="Spread", 
        xlab="Temps")
lines(colMeans(SP),col="red",lwd=1.5)

# comparaison selon notation pour maturité 5 ans
SPAAA <- c()#spread_i_JLT(N, 0, 0+5, paramJLT, M, D, 1, LGD)
for (t in Maturite50){SPAAA <- cbind(SPAAA,spread_i_JLT(N, t, t+5, paramJLT, M, D, 1, LGD))}
SPAA <- c()#spread_i_JLT(N, 0, 0+5, paramJLT, M, D, 2, LGD)
for (t in Maturite50){SPAA <- cbind(SPAA,spread_i_JLT(N, t, t+5, paramJLT, M, D, 2, LGD))}
SPA <- c()#spread_i_JLT(N, 0, 0+5, paramJLT, M, D, 3, LGD)
for (t in Maturite50){SPA <- cbind(SPA,spread_i_JLT(N, t, t+5, paramJLT, M, D, 3, LGD))}
SPBBB <- c()#spread_i_JLT(N, 0, 0+5, paramJLT, M, D, 4, LGD)
for (t in Maturite50){SPBBB <- cbind(SPBBB,spread_i_JLT(N, t, t+5, paramJLT, M, D, 4, LGD))}

plot(colMeans(SPAAA),type="l",ylim=c(-0.0001,0.0055),col="red",xlab="Temps",ylab="Spread",
     main=c("Evolution des spreads par rating de", "maturité 5 ans sur 50 ans (JLT)"))
lines(colMeans(SPAA),col="orange")
lines(colMeans(SPA),col="brown")
lines(colMeans(SPBBB),col="lightblue")
legend("topright",legend=c("AAA","AA", "A","BBB"),
       col=c("red","orange", "brown", "lightblue"),pch=20,
       cex=0.7)

# comparaison de spread A pour différente maturité 
SP1 <- c()
for (t in Maturite50){SP1 <- cbind(SP1,spread_i_JLT(N, t, t+1, paramJLT, M, D, 3, LGD))}
SP5 <- c()
for (t in Maturite50){SP5 <- cbind(SP5,spread_i_JLT(N, t, t+5, paramJLT, M, D, 3, LGD))}
SP10 <- c()
for (t in Maturite50){SP10 <- cbind(SP10,spread_i_JLT(N, t, t+10, paramJLT, M, D, 3, LGD))}
SP30 <- c()
for (t in Maturite50){SP30 <- cbind(SP30,spread_i_JLT(N, t, t+30, paramJLT, M, D, 3, LGD))}

plot(colMeans(SP1),type="l",ylim=c(0,0.005),xlab="Temps",ylab="Spread",
     main="Evolution du spread A pour différentes maturités (JLT)")
lines(colMeans(SP5),col="purple")
lines(colMeans(SP10),col="blue")
lines(colMeans(SP30),col="lightblue")
legend("topright",legend=c("maturité 1 an","maturité 5 ans", "maturité 10 ans","maturité 30 ans"),
       col=c("black","purple", "blue", "lightblue"),pch=20,
       cex=0.7)

# comparaison de spread B pour différente maturité 
SP1 <- c()
for (t in Maturite50){SP1 <- cbind(SP1,spread_i_JLT(N, t, t+1, paramJLT, M, D, 6, LGD))}
SP5 <- c()
for (t in Maturite50){SP5 <- cbind(SP5,spread_i_JLT(N, t, t+5, paramJLT, M, D, 6, LGD))}
SP10 <- c()
for (t in Maturite50){SP10 <- cbind(SP10,spread_i_JLT(N, t, t+10, paramJLT, M, D, 6, LGD))}
SP30 <- c()
for (t in Maturite50){SP30 <- cbind(SP30,spread_i_JLT(N, t, t+30, paramJLT, M, D, 6, LGD))}

plot(colMeans(SP1),type="l",ylim=c(0,0.07),xlab="Temps",ylab="Spread",
     main="Evolution du spread B pour différentes maturités (JLT)")
lines(colMeans(SP5),col="purple")
lines(colMeans(SP10),col="blue")
lines(colMeans(SP30),col="lightblue")
legend("topright",legend=c("maturité 1 an","maturité 5 ans", "maturité 10 ans","maturité 30 ans"),
       col=c("black","purple", "blue", "lightblue"),pch=20,
       cex=0.7)

### Test de martingalité sur Cash Flow JLT Vas ###################################################
plot(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 1, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 1, LGD)),
     ylim=c(0.99,1.01),"l",col="red",
     ylab="CashFlow actualisé",
     main=c("Test de martingalité sur CashFlow", "du modèle JLT (avec Vasicek)"))
lines(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 2, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 3, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 4, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 5, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite,PZCr_i_CF_JLT_Vas_FF(Maturite, paramVas, paramJLT, M, D, 6, LGD)/colMeans(PZCr_i_CF_JLT_Vas_sim.T(N=1000, 0, Maturite, paramVas, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)

### Test de martingalité sur PZCr avec rating JLT Vas ###################################################
t = 0.2
plot(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 1, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 1, LGD)),
     ylim=c(0.9,1.1),"l",col="red",
     xlab ="Maturité",
     ylab="PZC risqué actualisé",
     main=c("Test de martingalité sur PZC risqué", "du modèle JLT (avec Vasicek)"))
lines(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 2, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 3, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 4, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 5, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite30,PZCr_i_CF_JLT_Vas_FF(Maturite30, paramVas, paramJLT, M, D, 6, LGD)/colMeans(PZCr_i_JLT_Vas_sim.T(N=1000, t, Maturite30, paramVas, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)

### Test de martingalité sur Cash Flow JLT HW ###################################################
plot(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 1, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 1, LGD)),
     ylim=c(0.99,1.01),"l",col="red",
     ylab="CashFlow actualisé",
     main=c("Test de martingalité sur CashFlow", "du modèle JLT (avec HW)"))
lines(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 2, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 3, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 4, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 5, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite,PZCr_i_CF_JLT_HW_FF(Maturite, paramHW, paramJLT, M, D, 6, LGD)/colMeans(PZCr_i_CF_JLT_HW_sim.T(N=1000, 0, Maturite, paramHW, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)

### Test de martingalité sur PZCr avec rating JLT HW ###################################################
t = 0.1
plot(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 1, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 1, LGD)),
     ylim=c(0.9,1.1),"l",col="red",
     xlab ="Maturité",
     ylab="PZC risqué actualisé",
     main=c("Test de martingalité sur PZC risqué", "du modèle JLT (avec HW)"))
lines(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 2, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 2, LGD)),col="orange")
lines(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 3, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 3, LGD)),col="brown")
lines(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 4, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 4, LGD)),col="lightblue")
lines(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 5, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 5, LGD)),col="blue")
lines(Maturite30,PZCr_i_CF_JLT_HW_FF(Maturite30, paramHW, paramJLT, M, D, 6, LGD)/colMeans(PZCr_i_JLT_HW_sim.T(N=1000, t, Maturite30, paramHW, paramJLT, M, D, 6, LGD)),col="purple")
legend("topright",legend=c("AAA","AA", "A","BBB","BB","B"),
       col=c("red","orange", "brown", "lightblue", "blue", "purple"),pch=20,
       cex=0.6)

