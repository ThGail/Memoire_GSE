rm(list = ls())
library(readxl)
library(dfoptim)
library(writexl)

chemin = "/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire/memoire_GSE/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = "/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire/memoire_GSE/Mémoire GSE - Paris Dauphine & Prim'Act/Fonction"
chemin_cal = "/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire/memoire_GSE/Mémoire GSE - Paris Dauphine & Prim'Act/Calibrage"
setwd(chemin)

chemin_export_courbe_tx = "/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire/memoire_GSE/Best Estimate - Simbel/input/parametres/esg/ESG/Simulation_CourbeDesTaux"
chemin_export_credit = "/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire/memoire_GSE/Best Estimate - Simbel/input/parametres/esg/ESG/Simulation_Credits"


### Importation des donnees ###################################################
source(paste(chemin_cal,"/Vasicek_calibrage.R",sep=""))
source(paste(chemin_cal,"/HW_calibrage.R",sep=""))
source(paste(chemin_cal,"/BS_calibrage.R",sep=""))
source(paste(chemin_cal,"/CIR_calibrage.R",sep=""))
source(paste(chemin_cal,"/JLT_calibrage.R",sep=""))
source(paste(chemin_fct,"/GSE_fonction.R",sep=""))

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

### Parametrisation #######################################
# matrice de correlaion entre taux_court, action_glo et immo
Corr <- matrix(c(1,0.1,0.05,0.1,1,-0.3,0.05,-0.3,1),nrow=3)
# d?compostion de Cholesky :
L <- t(chol(Corr, pivot = FALSE))

N = 20; TT = 50; LGD = 0.3
projection <- 50
### Export GSE HW JLT - Taux Zéro-Coupon #####################################


for(t in 1:projection)
{
  print(paste("Année de projection : ", t))
  CourbeZC <- c()
  
  #Maturité 0.083333
  GSET50.test <- GSE_HW_JLT(N, t, 0.0833333, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.25
  GSET50.test <- GSE_HW_JLT(N, t, 0.25, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.5
  GSET50.test <- GSE_HW_JLT(N, t, 0.5, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 0.75
  GSET50.test <- GSE_HW_JLT(N, t, 0.75, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  
  
  for(Tt in 1:TT)
  {
    GSET50.test <- GSE_HW_JLT(N, t, Tt, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
    CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  }
  
  #Maturité 40
  GSET50.test <- GSE_HW_JLT(N, t, 40, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  #Maturité 50
  GSET50.test <- GSE_HW_JLT(N, t, 50, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  CourbeZC <- cbind(CourbeZC, GSET50.test$TZC_HW)
  
  CourbeZC.dataframe <- as.data.frame(CourbeZC)
  
  #On renomme les colonnes
  names(CourbeZC.dataframe)[1] <- 0.0833333
  names(CourbeZC.dataframe)[2] <- 0.25
  names(CourbeZC.dataframe)[3] <- 0.5
  names(CourbeZC.dataframe)[4] <- 0.75
  for(i in 5:TT)
  {
    names(CourbeZC.dataframe)[i] <- i
  }
  write_xlsx(CourbeZC.dataframe,paste(chemin_export_courbe_tx,"/courbeZC", t, ".xlsx", sep=""))
}

### Export GSE HW JLT - Risque de crédit #####################################

Proba_survie_AAA <- c()
Proba_survie_AA <- c()
Proba_survie_A <- c()
Proba_survie_BBB <- c()
Proba_survie_BB <- c()
Proba_survie_B <- c()

GSET50.test <- GSE_HW_JLT(N, 0, 50, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
Proba_survie_AAA <- cbind(Proba_survie_AAA, GSET50.test$ps_AAA)
Proba_survie_AA <- cbind(Proba_survie_AA, GSET50.test$ps_AA)
Proba_survie_A <- cbind(Proba_survie_A, GSET50.test$ps_A)
Proba_survie_BBB <- cbind(Proba_survie_BBB, GSET50.test$ps_BBB)
Proba_survie_BB <- cbind(Proba_survie_BB, GSET50.test$ps_BB)
Proba_survie_B <- cbind(Proba_survie_B, GSET50.test$ps_B)



for(i in 1:TT)
{
  print(paste("Année de projection : ", i))
  GSET50.test <- GSE_HW_JLT(N, i, 50, paramHW, paramBS$action, paramBS$immo, paramJLT, L, LGD)
  Proba_survie_AAA <- cbind(Proba_survie_AAA, GSET50.test$ps_AAA)
  Proba_survie_AA <- cbind(Proba_survie_AA, GSET50.test$ps_AA)
  Proba_survie_A <- cbind(Proba_survie_A, GSET50.test$ps_A)
  Proba_survie_BBB <- cbind(Proba_survie_BBB, GSET50.test$ps_BBB)
  Proba_survie_BB <- cbind(Proba_survie_BB, GSET50.test$ps_BB)
  Proba_survie_B <- cbind(Proba_survie_B, GSET50.test$ps_B)
}

Proba_survie_AAA.dataframe <- as.data.frame(Proba_survie_AAA)
Proba_survie_AA.dataframe <- as.data.frame(Proba_survie_AA)
Proba_survie_A.dataframe <- as.data.frame(Proba_survie_A)
Proba_survie_BBB.dataframe <- as.data.frame(Proba_survie_BBB)
Proba_survie_BB.dataframe <- as.data.frame(Proba_survie_BB)
Proba_survie_B.dataframe <- as.data.frame(Proba_survie_B)

names(Proba_survie_AAA.dataframe)[1] <- 0
names(Proba_survie_AA.dataframe)[1] <- 0
names(Proba_survie_A.dataframe)[1] <- 0
names(Proba_survie_BBB.dataframe)[1] <- 0
names(Proba_survie_BB.dataframe)[1] <- 0
names(Proba_survie_B.dataframe)[1] <- 0

for(i in 2:TT)
{
  names(Proba_survie_AAA.dataframe)[i] <- i-1
}
write_xlsx(Proba_survie_AAA.dataframe,paste(chemin_export_credit,"/Proba_survie_AAA", t, ".xlsx", sep=""))
write_xlsx(Proba_survie_AA.dataframe,paste(chemin_export_credit,"/Proba_survie_AA", t, ".xlsx", sep=""))
write_xlsx(Proba_survie_A.dataframe,paste(chemin_export_credit,"/Proba_survie_A", t, ".xlsx", sep=""))
write_xlsx(Proba_survie_BBB.dataframe,paste(chemin_export_credit,"/Proba_survie_BBB", t, ".xlsx", sep=""))
write_xlsx(Proba_survie_BB.dataframe,paste(chemin_export_credit,"/Proba_survie_BB", t, ".xlsx", sep=""))
write_xlsx(Proba_survie_B.dataframe,paste(chemin_export_credit,"/Proba_survie_B", t, ".xlsx", sep=""))




