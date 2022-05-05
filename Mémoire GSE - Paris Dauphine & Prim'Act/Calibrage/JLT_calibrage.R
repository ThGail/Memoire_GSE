# procédure de calibrage du modèle JLT

# définition de l'environnement
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = paste(chemin,"/Fonction", sep="")
setwd(chemin)

# librairie utile
library(readxl)
library(dfoptim)

# importation fonctions utiles du JLT
source(paste(chemin_fct,"\\JLT_fonction.R",sep=""))

# importation données utiles pour le calibrage
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

### PREPARATION MATRICE HISTORIQUE ##################################
# exemple de matrice de transition de rating 2017
Q = matrix(c(91.06,8.25,0.6,0.07,0.02,0,0,0,0,
             0.86,89.58,8.96,0.45,0.07,0.04,0.02,0,0.02,
             0.05,2.58,91.02,5.63,0.51,0.11,0.04,0.01,0.05,
             0.03,0.15,4.34,90.41,3.99,0.72,0.16,0.02,0.18,
             0.01,0.04,0.46,6.66,83.17,7.81,0.77,0.12,0.96,
             0.01,0.03,0.16,0.5,5.35,82.3,7.41,0.58,3.66,
             0,0.01,0.03,0.1,0.4,7.6,79.24,3.33,9.29,
             0,0,0.06,0,0.72,2.94,11.48,50.57,34.23,
             0,0,0,0,0,0,0,0,100)/100,ncol=9,byrow=T)

# transformation pour la positivité hors diag
Ltemp=logMatrice(Q,100)
Lpos = pmax(Ltemp,0)
Lneg = pmin(Ltemp,0)
L = Lpos + diag(rowSums(Lneg))

# diagonalisation de la matrice L
M = eigen(L)$vectors
D = diag(eigen(L)$values)
# round(M%*%D%*%solve(M),15) permet de retrouver la matrice de base (modulo erreur numérique)

### CALIBRAGE #############################################
TT <- 1 # les données fornis sont les spreads de maturité 1 an
LGD <- 0.3

#On considere param comme une liste/matrice, chaque ligne represente les parametres d une classe de rating
#Sachant qu'il y a 6 classes de rating et 4 parametres pour chacune d entre elle,
#Les parametres sont c(kAAA, kAA, kA, kBaa, kBa, kB, 
#                      muAAA, muAA, muA, muBaa, muBa, muB, 
#                      sigma, 
#                      pi0AAA, pi0AA, pi0A, pi0BBB, pi0BB, pi0B)

param_init <- c(0.327, 0.291, 0.200, 0.060, 0.004, 0.003, 0.003, 0.003,
                0.197, 0.325, 0.574, 0.543, 1.628, 2.018, 2.000, 2.000,
                0.01,
                0.003, 0.001, 0.001, 0.015, 0.054, 0.055, 0.055, 0.055)

LB <- c(rep(0,8),rep(0,8),0,rep(0,8))
UB <- c(rep(30,8),rep(20,8),10,rep(30,8))
#paramJLT = hjkb(param_init,Ecart_JLT,lower=LB,upper=UB)$par # 2.731286e-05
paramJLT = nlminb(start = param_init,ecart_JLT,lower=LB,upper=UB)$par
#paramJLT = optim(par=param_init,fn=Ecart_JLT,lower=LB,upper=UB,method="L-BFGS-B")$par
(listJLT <- list(
  k = paramJLT[1:8],
  mu = paramJLT[9:16],
  sigma = paramJLT[17],
  pi0 = paramJLT[18:25]))
paramJLT[17] = min(sqrt(2*listJLT$k*listJLT$mu)[sqrt(2*listJLT$k*listJLT$mu) > 0])
# il faut que sigma soit non nul
(listJLT$sigma <- paramJLT[17])

