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

N <- 10000
r0 <- TauxZC[1]

Calculrt_Vas <- function(N, t, param, r0){
  # N nb de simulation
  # t année de projection
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(mean(rnorm(N,r0*exp(-a*t)+b*(1-exp(-a*t)), sigma*sqrt((1-exp(-2*a*t))/(2*a)))))
}

prixZC_Vas <- function(t, TT, param, rt){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  theta <- TT-t
  return(exp(-b*theta)*exp(-(rt-b)*(1-exp(-a*theta))/a + 0.5*(sigma^2*theta/a^2-sigma^2/a^3*(1-exp(-a*theta))-sigma^2/(2*a^3)*(1-exp(-a*theta))^2)))
}

TauxZCFF_fct <- function(N, t, TT, param){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  theta <- TT-t
  return(Ri - ((Ri-Calculrt_Vas(N,TT,param,r0))*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}

Ecart <- function(param){
  e <- 0
  for (i in 1:150){
    tzp<-TauxZCFF_fct(N, 0, i, param)
    e <- e + (tzp-TauxZC[i])^2
  }
  return(e)
}

#======= Régression linéaire ======== 
# Calibrage des parametres a b avec la courbe taux ZC
# Rt+1 = Rt*e(-a) + b*(1-e(-a)) + epsilon

TT <- length(TauxZC) # la date de maturite
Tx_t1 <- TauxZC[1:(TT-1)]
Tx_t2 <- TauxZC[2:(TT)]
reg <- lm(Tx_t2~Tx_t1)
summary(reg)
a_reg <- -log(coef(reg)[2])
b_reg <- coef(reg)[1]/(1-coef(reg)[2])
#On choisi les valeurs suivantes ici d'une certaine simulation : 
#a_reg vaut 0.01298906
#b_reg vaut 0.03872287

#On détermine sigma en fiasant un MCO en fixant a et b :
LB = c(a_reg, b_reg, 0.00000001)
UB = c(a_reg, b_reg, 0.5)

#On fait un MCO avec les valeurs 
param <- hjkb(c(a_reg, b_reg, 0.1), Ecart, lower=LB,upper=UB)$par

test <- c()
for (i in 1:150){
  test<-c(test,TauxZCFF_fct(N, 0, i, param))
}

#En répétant cette sans fixer a et b on obtient des valeurs proches de : 
# param a = 0.008738708 b = 0.041931152 sigma = 0.001186371
LB_free = c(0, 0, 0.00000001)
UB_free = c(1, 1, 0.5)
param_free <- hjkb(c(a_reg, b_reg, 0.0015), Ecart, lower=LB_free,upper=UB_free)$par
#param_free <- c(0.008738708, 0.041931152, 0.001186371)
test_free <- c()
for (i in 1:150){
  test_free<-c(test_free,TauxZCFF_fct(N, 0, i, param_free))
}



plot(test, main="Comparaison simulation TZC et Tx Réel", type='l', col='red', xlab='Maturité', ylab='Taux Zéro-Coupon')
lines(test_free, col="blue")
lines(TauxZC)
legend("bottomright", c("RFR réel", "RFR sim a-b fixés", "RFR sim a-b libres"),
       col = c("black", "red", "blue"), lty = c(1, 1), cex = 0.8)
