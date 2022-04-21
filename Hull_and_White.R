rm(list = ls())

library(readxl) 
library(dfoptim)

setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
Maturite30 <- Maturite[1:30]
TauxZC    <- as.numeric(sheetDGlo$...2 [7:156])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

#### Les paramètres ####
#Dans ce modèle, il y a deux paramètres : a et sigma
#param <- c(a, sigma)

#### Fonction utiles ####
#On note TxFI le tx_foward_instantane
TxFI_fct <- function(PZC, t, TT)
{
  #code pour ce taux
}

TxFI_deriv_fct <- function(TxFI)
{
  #On dérive le TxFI
}

b_t_fct <- function(TxFI, derivTxFI, t, TT, param)
{
  #code b_t
}

Tx_HW_fct <- function(n, t, r0, param)
{
  #code Taux instantané Hull & White
  #Simulation à l'aide d'une loi normale
}

A_tT_fct <- function(TxFI,derivTxFI, t, TT, param)
{
  #calcul de A
}

B_tT_fct <- function(t, TT, param)
{
  #calcul de B
}

Tx_HW_vers_PZC <- function(TxFI, derivTxFI, t, TT, Tx_HW, param)
{
  #calcul du prix ZC
}


