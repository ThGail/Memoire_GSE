rm(list = ls())
library(readxl) 
library(dfoptim)

#### Importation des données ###################################################
setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet=1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet=3)

Maturite  <- as.numeric(sheetDGlo$`Courbe des taux ZC` [7:156])
TauxZC    <- as.numeric(sheetDGlo$...2 [7:156])
TauxIflZC <- as.numeric(sheetDGlo$...5 [7:106])
PrixCaps  <- as.numeric(sheetCali$...5 [3:22])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

N <- 1000
r0 <- TauxZC[1]



A_fct2 <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_fct2 <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return((2 * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# On suppose que le mouvement brownien est standard
lambdat_fct2 <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda_t <- param[4]
  for (i in 1:t) {
    lambda_t <- lambda_t + (1 / N) * k * (mu - lambda_t) + sigma * sqrt(lambda_t * (lambda_t > 0))* rnorm(N, 0, 1)
  }
  return(lambda_t)
}

survie_FF_fct2 <- function(N, t, TT, param) {
  # On a calculé la moyenne en sortie de fonction et non directement dans la fonction
  lambda_t <- mean(lambdat_fct(N, t, param))
  B <- B_fct(TT - t, param)
  A <- A_fct(TT - t, param)
  return(A * exp(-B * lambda_t))
}


Prix_Spread_fct_bis <- function(N, t, LGD, TT,param){
  return((-1/(TT-t))*(1-LGD + LDG * survie_FF_fct(N, t, TT, param)))
}

k <- 0.5
mu <- 0.5
sigma <- 0.5
TT <- 50
LGD <- 0.6
N <- 100
lambda0 <- 0.5



Ecart <- function(N, LGD, TT, param){
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <-sqrt(k**2 + 2 * sigma**2)
  e <- 0
  for (i in 1:length(SpreadMarket)){
    tzp <- Prix_Spread_fct(N, i, LGD, TT, param)
    e <- e + (tzp-SpreadMarket[i])^2
  }
  return(e)
}
Prix_Spread_fct <- function(N, t, LGD, TT,param){
  return((1-LGD + LGD * survie_FF_fct(N, t, TT, param))**(-1/(TT-t)))
}








rm(list = ls())

library(readxl) 
library(dfoptim)

setwd("~/Université/Master 1_Semestre 1_2021-2022/Mémoire de M1")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])
SpreadMarket  <- as.numeric(sheetDGlo$...15[3:9])

N <- 1000

A_fct <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k*2 + 2 * sigma*2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))*((2 * k * mu) / (sigma*2)))
}

B_fct <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k*2 + 2 * sigma*2)
  return((2 * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# On suppose que le mouvement brownien est standard
lambdat_fct <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda_t <- param[4]
  if (t==0) {return(lambda_t)}
  for (i in 1:t) {
    lambda_t <- lambda_t + (1 / N) * k * (mu - lambda_t) + sigma * sqrt(lambda_t* (lambda_t>0) ) * rnorm(N, 0, 1)
  }
  return(lambda_t)
}

survie_FF_fct <- function(N, t, TT, param) {
  # On a calculÃ©️ la moyenne en sortie de fonction et non directement dans la fonction
  lambda_t <- mean(lambdat_fct(N, t, param))
  B <- B_fct(TT - t, param)
  A <- A_fct(TT - t, param)
  return(A * exp(-B * lambda_t))
}

Prix_Spread_fct <- function(t,TT,param,LGD){
  return((1 - LGD + LGD * survie_FF_fct(N,t,TT,param) )^(-1/(TT-t))-1)
}

#### Fonction écart à la moyenne pour le calibrage ####
Ecart_CIR <- function(param){
  e <- 0
  for (t in 0:TT){
    spread <- Prix_Spread_fct(0,t,param,LGD)
    print(t)
    #print(spread)
    print()
    print(spread-spread_mkt)
    e <- e + (spread-spread_mkt)^2
    #print(e)
  }
  return(e)
}

LGD <- 0.4
param <- c(0.5, 0.5, 0.5, 0.5)
spread_mkt <- SpreadMarket[6]
TT <- 150

Ecart_CIR(param)

LB = c(0,0,0,0)
UB = c(1,1,1,1)
param <- hjkb(c(0.5, 0.5, 0.5, 0.5),Ecart_CIR,lower=LB,upper=UB)$par
param

tqt <- numeric(TT)
for(i in 1:TT){
  tqt[i] <- Prix_Spread_fct(0, i, param, LGD) 
}
plot(tqt, ylim= c(0,0.05))
abline(h= spread_mkt, col= 'red')

