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

### Generateur de mouvement brownien N * n_type ######################
# n_type est le nombre de donnée intéret 
# ex : TZC pour modèle de taux, spread pour modèle de crédit, etc
MvtBrow <- function(N,n_type){}
