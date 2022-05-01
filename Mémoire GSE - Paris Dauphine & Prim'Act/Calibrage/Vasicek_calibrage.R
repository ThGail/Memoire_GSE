# procédure de calibrage du modèle Vasicek

# définition de l'environnement
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Fonction"
setwd(chemin)

# librairie utile
library(readxl)
library(dfoptim)

# importation fonctions utiles du Vasicek
source(paste(chemin_fct,"\\Vasicek_fonction.R",sep=""))

# importation données utiles pour le calibrage
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])

### CALIBRAGE ####################################################
param_init = c(0.005,0.005,0.005)
LB = c(0,0,1e-3)
UB = c(1,1,1)
(paramVas = hjkb(param_init,ecart_Vas,lower=LB,upper=UB)$par)
ecart_Vas(paramVas)
