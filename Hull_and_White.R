rm(list = ls())

library(readxl)
library(dfoptim)

setwd("C:/Users/dgahe/Documents/thibault_cours")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
Maturite30 <- Maturite[1:30]
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
Caps <- as.numeric(sheetCali$...5[3:22])
SpreadMarket <- as.numeric(sheetDGlo$...15[3:9])

r0 <- TauxZC[1]
TT <- length(TauxZC)
n <- 1000

# Affichage TauxZC et PrixZC :
plot(TauxZC)
PrixZC <- exp(-Maturite * TauxZC)

# Fonction passage au Prix Zéro coupon :
prixZC_Vas <- function(t, TT, param, rt) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  theta <- TT - t
  return(exp(-b * theta) * exp(-(rt - b) * (1 - exp(-a * theta)) / a + 0.5 * (sigma^2 * theta / a^2 - sigma^2 / a^3 * (1 - exp(-a * theta)) - sigma^2 / (2 * a^3) * (1 - exp(-a * theta))^2)))
}

#### Les paramètres ####
# Dans ce modèle, il y a deux paramètres : a et sigma
# param <- c(a, sigma)
paramHW <- c(0.5, 0.5)

# Calcul de la dérivée par interpolation :
# On note TxFI le tx_foward_instantane

logPrixZC_fct <- splinefun(x = 1:150, y = log(PrixZC), method = "natural")

TxFI <- -logPrixZC_fct(1:150, 1)

plot(-log(PrixZC), main = "-log(Prix_ZC)")
plot(TxFI, main = "Taux Forward Instantané")

TxFI_fct <- splinefun(x = 1:150, y = TxFI, method = "natural")

TxFI_deriv <- TxFI_fct(1:150, 1)
plot(TxFI_deriv, main = "Dérivée du taux forward instantané")


#### Fonction utiles ####

beta_fct <- function(t, param) {
  a <- param[1]
  sigma <- param[2]
  if(t == 0)
  {
    TxFI <- 0
  }
  else
  {
    TxFI <- TxFI[t]
  }
  return(TxFI + ((sigma^2) / (2 * a)) * (1 - exp(-a * (t)))^2)
}

tauxHW <- function(n, t, param, r0 = TauxZC[1]) {
  a <- param[1]
  sigma <- param[2]
  return(r0 * exp(-a * (t)) + beta_fct(t, param) - beta_fct(0, param) * exp(-a * (t)) + sqrt(((sigma^2) / (2 * a)) * (1 - exp(-2 * a * (t)))) * rnorm(n, 0, 1))
}
# simulation PZC des trajectoires pour une maturité donnée
tauxHW.T <- Vectorize(tauxHW, "TT")


B_fct <- function(t, TT, param) {
  a <- param[1]
  sigma <- param[2]
  return((1 / a) * (1 - exp(-a * (TT - t))))
}

A_fct <- function(t, TT, param) {
  a <- param[1]
  sigma <- param[2]
  B <- B_fct(t, TT, param)
  
  if(t == 0)
  {
    PrixZC_t = 1
    TxFI <- 0
  }
  else
  {
    PrixZC_t <- PrixZC[t]
    TxFI <- TxFI[t]
  }
  
  return((PrixZC[TT] / PrixZC_t) * exp(B * TxFI - ((sigma^2) / (4 * a)) * (1 - exp(-2 * a * t)) * B^2))
}

PrixZC_HW_sim <- function(n, t, TT, param) {
  if(t == 0 & TT == 0)
  {
    return(1)
  }
  return(A_fct(t, TT, param) * exp(-B_fct(t, TT, param) * tauxHW(n, t, param)))
}
# simulation PZC des trajectoires pour une maturité donnée
PrixZC_HW_sim.T <- Vectorize(PrixZC_HW_sim, "TT")

matplot(t(PrixZC_HW_sim.T(10, 1, 1:5, paramHW)))

PrixZC_HW <- function(n, t, TT, param) {
  return(mean(PrixZC_HW_sim(n, t, TT, param)))
}
PrixZC_HW.T <- Vectorize(PrixZC_HW, "TT")

# Transformation en Caps :

Caps_HW <- function(n, K, t, TT, param) {
  a <- param[1]
  sigma <- param[2]

  cap <- 0
  for (ti in 1:TT)
  {
    print(ti)
    PZCti <- PrixZC_HW_sim(n, t, ti, param)
    print(PZCti)
    PZCti_1 <- PrixZC_HW_sim(n, t, ti - 1, param)
    print(PZCti_1)
    v <- sigma * sqrt((1 - exp(-2 * a * (ti - t))) / (2 * a)) * B_fct(ti - 1, ti, param)
    print(v)
    d1 <- (1 / v) * log((PZCti * (1 + K)) / PZCti_1) + v / 2
    print(d1)
    print(pnorm(d1))
    cap <- cap + (1 + K) * PZCti * pnorm(d1) - PZCti_1 * pnorm(d1 - v)
    print(cap)
  }

  return(cap)
}


# Calibrage des paramètres :




Ecart_HW <- function(param) {
  e <- 0
  for (t in 1:TT) {
    e <- e + (PrixZC_HW(n, 0, TT, param) - PrixZC[t])^2
  }
  return(e)
}

param_init <- c(0.5, 0.5)
LB <- c(-1, 0)
UB <- c(1, 1)
paramHW <- hjkb(param_init, Ecart_HW, lower = LB, upper = UB)$par
# 0.004084473 0.292354584 0.003962402
# erreur = 0.0005183257


