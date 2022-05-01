# procédure de calibrage du modèle CIR

# définition de l'environnement
chemin = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act"
chemin_fct = "C:/Users/mengr/OneDrive/Documents/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1/Mémoire GSE - Paris Dauphine & Prim'Act/Fonction"
setwd(chemin)

# librairie utile
library(readxl)
library(dfoptim)

# importation fonctions utiles du CIR
source(paste(chemin_fct,"\\CIR_fonction.R",sep=""))

# importation données utiles pour le calibrage
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

### CALIBRAGE #####################################################
TT = 1; LGD = 0.3

paramCIR <- list()
param_init = c(0.2,0.2,0.2,0.2)
LB = c(0,0,0,0)
UB = c(10,10,10,10)

#AAA : 1.599722e-01 1.288132e-01 2.030100e-01 2.126438e-13
spread_mkt <- SpreadMarket[1]
#paramCIR$AAA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AAA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$AAA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AAA[3] <- min(paramCIR$AAA[3],sqrt(2*paramCIR$AAA[1]*paramCIR$AAA[2]))

#AA : 0.1793781 0.1526860 0.2028968 0.0000000
spread_mkt <- SpreadMarket[2]
#paramCIR$AA <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$AA <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$AA <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$AA[3] <- min(paramCIR$AA[3],sqrt(2*paramCIR$AA[1]*paramCIR$AA[2]))

#A : 0.1960534 0.1731844 0.2027408 0.0000000
spread_mkt <- SpreadMarket[3]
#paramCIR$A <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$A <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$A <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$A[3] <- min(paramCIR$A[3],sqrt(2*paramCIR$A[1]*paramCIR$A[2]))

#BBB : 0.19780523 0.17654727 0.20257493 0.01274617
spread_mkt <- SpreadMarket[4]
#paramCIR$BBB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BBB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$BBB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BBB[3] <- min(paramCIR$BBB[3],sqrt(2*paramCIR$BBB[1]*paramCIR$BBB[2]))

#BB : 0.19887557 0.18437281 0.20176218 0.04921644
spread_mkt <- SpreadMarket[5]
#paramCIR$BB <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$BB <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$BB <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$BB[3] <- min(paramCIR$BB[3],sqrt(2*paramCIR$BB[1]*paramCIR$BB[2]))

#B : 0.19890474 0.18459973 0.20173836 0.05140548
spread_mkt <- SpreadMarket[6]
#paramCIR$B <- hjkb(param_init,ecart_CIR,lower=LB,upper=UB)$par
paramCIR$B <- nlminb(start = param_init,ecart_CIR,lower=LB,upper=UB)$par
#paramCIR$B <- optim(par=param_init,fn=ecart_CIR,lower=LB,upper=UB,method="L-BFGS-B")$par
paramCIR$B[3] <- min(paramCIR$B[3],sqrt(2*paramCIR$B[1]*paramCIR$B[2]))


