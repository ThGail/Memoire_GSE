rm(list = ls())
library(readxl)

#### Inportation des donn?es ###################################################
setwd("/Users/thibaultgaillard/Documents/M1 Actuariat/Mémoire")
sheetDGlo <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 1)
sheetCali <- read_excel("Input_20210118_18h41m33s.xlsm", sheet = 3)
sheetEIOPA <- read_excel("RFR_for_publication_RFR_curves.xlsx", sheet = 3)

Maturite <- as.numeric(sheetDGlo$`Courbe des taux ZC`[7:156])
TauxZC <- as.numeric(sheetDGlo$...2[7:156])
TauxIflZC <- as.numeric(sheetDGlo$...5[7:106])
PrixCaps <- as.numeric(sheetCali$...5[3:22])
CourbeRFC <- as.numeric(sheetEIOPA$Euro[9:158])


#### Vasicek ###################################################################
######## Paramètres #########
# Ceux issu par moindre carrés ordinaires :
a <- 0.008738708
b <- 0.041931152
sigma <- 0.001186371

# Méthode par discrétisation, on fixe a t b puis on fait un MCO sur le sigma
a2 <- 0.01298906
b2 <- 0.03872287
#Obtenu après avoir fait un MCO en fixant a et b
sigma2 <- 0.002397156

R0 <- TauxZC[1]

N <- 10000
TT <- length(TauxZC)
###########################

plot(CourbeRFC, main="Comparaison Taux sans risque EIOPA & TZC", xlab="Maturité", ylab="Taux ZC/Sans risque")
lines(TauxZC, col = "red")
legend("bottomright", c("Taux sans risque données", "Taux zéro coupon données"),
  col = c("black", "red"), lty = c(1, 1), cex = 0.8
)

# Fonction taux Vasicek avec formule ferm?e
# Permet de simuler le taux d'intérêt sans risque (ou taux ZC)
TauxVasicek_fct <- function(n, t, r0, param) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(rnorm(n, r0 * exp(-a * t) + b * (1 - exp(-a * t)), sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a))))
}

TauxZCFF_fct <- function(N, t, TT, param, r0){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  theta <- TT-t
  return(Ri - ((Ri-Calculrt_Vas(N,TT,param,r0))*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}

# Même fonction, on calcule directement la moyenne dedans, permet de simplifier le code par moments :
Calculrt_Vas <- function(N, t, param, r0) {
  # N nb de simulation
  # t année de projection
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(mean(rnorm(N, r0 * exp(-a * t) + b * (1 - exp(-a * t)), sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a)))))
}

# ============ TEST viabilité taux ZC ==============
# On compare la valeur empirique du taux sans risque avec la courbe EIOPA
rep <- c()
rep2 <- c()
for (i in 0:TT) {
  rep <- c(rep, mean(TauxVasicek_fct(N, i, R0, c(a2, b2, sigma2))))
  rep2 <-c(rep2, TauxZCFF_fct(N, 0, i, c(a, b, sigma), R0))
}
plot(TauxZC, main="Comparaison simulation Taux ZC/Taux ZC Réel", ylab = "Taux ZC", xlab = "Maturité", type = "l")
#lines(rep, col = "red")
lines(rep2, col = "red")
legend("bottomright", c("Taux ZC réel", "Taux ZC sim"),
  col = c("black", "red"), lty = c(1, 1), cex = 0.8)
# =================================
# Fonction prix ZC avec formule fermée
PrixZC_fct <- function(n, t, TT, r0, param) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b * (TT - t)) * exp(-(TauxVasicek_fct(N, TT, R0, param) - b) * (1 - exp(-a * (TT - t))) / a + 0.5 * (sigma^2 * (TT - t) / a^2 - sigma^2 / a^3 * (1 - exp(-a * (TT - t))) - sigma^2 / (2 * a^3) * (1 - exp(-a * (TT - t)))^2)))
}
prixZC_Vas <- function(t, TT, param, rt) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b * (TT - t)) * exp(-(rt - b) * (1 - exp(-a * (TT - t))) / a + 0.5 * (sigma^2 * (TT - t) / a^2 - sigma^2 / a^3 * (1 - exp(-a * (TT - t))) - sigma^2 / (2 * a^3) * (1 - exp(-a * (TT - t)))^2)))
}

#========== TEST Prix zéro COUPON ==========
Ptest <- c()
Ptest2 <- c()
for (i in 0:TT) {
  Ptest <- c(Ptest, mean(PrixZC_fct(N, 0, i, R0, c(a, b, sigma))))
  Ptest2 <- c(Ptest2, prixZC_Vas(0, i, c(a, b, sigma), mean(TauxVasicek_fct(N, i, R0, c(a, b, sigma)))))
}
plot(exp(-Maturite * TauxZC), main="Prix ZC réel vs Prix ZC simulé", type = "l", xlab = "Maturité", ylab = "Prix ZC")
#lines(Ptest, col = "blue", pch = 20)
lines(Ptest2, col = "red", pch = 20)
legend("topright", c("PZC données", "Simulation PZC"),
  col = c("black", "red"), lty = c(1, 1), cex = 0.8)
#============================================
# Comparaison des taux
TauxZC_FPassage_fct <- function(n, t, TT, r0, param) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(-log(PrixZC_fct(n, t, TT, r0, param)) / (TT - t))
}
TauxZC_FF_fct <- function(n, t, TT, r0, param) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b - sigma^2 / (2 * a^2)
  theta <- TT - t
  return(Ri - ((Ri - mean(TauxVasicek_fct(N, TT, R0, param))) * (1 - exp(-a * theta)) - sigma^2 / (4 * a^2) * (1 - exp(-a * theta))^2) / (a * theta))
}

# ===== TEST Taux zéro coupon ==========
TZCtest <- c()
TZCtest2 <- c()
for (i in 1:150) {
  TZCtest <- c(TZCtest, mean(TauxZC_FPassage_fct(N, 0, i, R0, c(a, b, sigma))))
  TZCtest2 <- c(TZCtest2, mean(TauxZC_FF_fct(N, 0, i, R0, c(a, b, sigma))))
}
plot(TauxZC, type = 'l', xlab = "Maturité", main="Comparaison Simulation TZC")
lines(TZCtest, col = "blue", lw=4)
lines(TZCtest2, col = "red")
legend("bottomright", c("TZC données", "TZC via Formule de Passage", "TZC via Formule Fermée"),
  col = c("black", "blue", "red"), lty = c(1, 1, 1), cex = 0.6)
# ========================================

# ======================================
# ===== MODELE BLACK & SCHOLES ========
# ======================================

# tau = T-t est la duree restante jusqu'à la maturité
S0 <- 1
K <- 1
TT <- 1
t <- 0
tau <- 1 # Par défaut
N <- 100000 # Nombre de simulation Monte Carlo
# Paramètres issus des données :
sigma.action <- 0.152396
sigma.immo <- 0.10

# Prix de l'actif par simulation
PrixactifSt_fct <- function(n, t, S0, r, sigma) {
  return(S0 * exp((r - 0.5 * sigma^2) * t + sigma * rnorm(n, 0, sqrt(t))))
}

# ====================================
# Exemple : Prix action pour une maturité de 1 pour 100 simulations
plot(PrixactifSt_fct(n = 100, t = 1, S0, mean(TauxVasicek_fct(n = 100, t = 1, R0, c(a, b, sigma.action))), sigma.action), main="Simulation de prix de l'actif", xlab="Temps", ylab='Prix')

# Evolution du prix de l'action en fonction de la maturité sur 50 ans
Prix.action <- c()
for (t in 0:50) {
  Prix.action <- c(Prix.action, mean(PrixactifSt_fct(N, t, S0, Calculrt_Vas(N, t, c(a, b, sigma.action), R0), sigma.action)))
}
plot(Prix.action, main = "Simulation du prix du sous-jacent action", xlab = "Maturité", ylab = "Prix sous-jacent Action")
# =====================================
# Evolution du prix de l'immo en fonction de la maturité sur 50 ans
Prix.immo <- c()
for (t in 0:50) {
  Prix.immo <- c(Prix.immo, mean(PrixactifSt_fct(N, t, S0, Calculrt_Vas(N, t, c(a, b, sigma.immo), R0), sigma.immo)))
}
plot(Prix.immo,
  main = "Simulation du prix du sous-jacent immobilier", xlab = "Maturité", ylab = "Prix sous-jacent Immobilier"
)

#### Vérification ##############################################################
# Prix du call par simulation
bscallsimul_fct <- function(n,tau,S0,K,r,sigma){
  return(exp(-r*tau)*pmax(PrixactifSt_fct(n,tau,S0,r,sigma)-K,0))
}

# Prix call théorique
bscall_fct <- function(tau,S0,K,r,sigma) {
  d1 <- (log(S0/K) + (r + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
  d2 <- d1 - sigma*sqrt(tau)
  if (tau < 1) {return(max(S0-K,0))} 
  else {return(S0*pnorm(d1) - K*exp(-r*tau)*pnorm(d2))}
}

# Convergence du call issu du prix action pour maturité 1 an
call.act.simul <- c()
for (i in seq(0, N, 100)){
  call.act.simul <- c(call.act.simul, mean(bscallsimul_fct(i,tau,S0,K,TauxZC[tau],sigma.action)))
}
plot(seq(0,N,100),call.act.simul, 
     main=c("Convergence du call issu du prix action pour maturité 1 an"),pch=20,
     xlab="Nombre de simulations",ylab="Prix call")
abline(h=bscall_fct(tau,S0,K,TauxZC[tau],sigma.action),lty=2,lwd=3,col="red")

# Convergence du call issu du prix immo pour maturité 1 an
call.immo.simul <- c()
for (i in seq(0,N,100)){
  call.immo.simul <- c(call.immo.simul, mean(bscallsimul_fct(i,tau,S0,K,TauxZC[tau],sigma.immo)))
}
plot(seq(0,n,100),call.immo.simul, 
     main=c("Convergence du call issu du prix immobilier pour maturité 1 an"),pch=20,
     xlab="",ylab="Prix call")
abline(h=bscall_fct(tau,S0,K,TauxZC[tau],sigma.immo),lty=2,lwd=3,col="red")

#### Validation du modèle : Test de martingalité ###############################
Prixaction.act_fct <- function(n,t,S0,r,sigma){
  return(exp(-r*t)*S0*exp((r-0.5*sigma^2)*t+sigma*rnorm(n,0,sqrt(t))))
}

# Test pour action
Prix.action.act <- c()
for (t in 0:150){Prix.action.act <- c(Prix.action.act, mean(Prixaction.act_fct(N,t,S0,TauxZC[t],sigma.action)))}
plot(Prix.action.act,ylim=c(0.9,1.1),pch=20,col="darkgrey",
     main="Test de martingalité pour les actions",
     xlab="Maturité",
     ylab="Moyenne de l'indice actualisé ")
abline(h=S0,col="red",lty=3,lwd=2)

abline(h=S0+0.025,col="blue",lty=4)
abline(h=S0-0.025,col="blue",lty=4)
legend("topleft",legend=c("Prix initial","Borne à 5%"),
       col=c("red","blue"),pch=20,
       cex=0.8)

# Test pour immo
Prix.immo.act <- c()
for (t in 0:150){Prix.immo.act <- c(Prix.immo.act, mean(Prixaction.act_fct(N,t,S0,TauxZC[t],sigma.immo)))}
plot(Prix.immo.act,ylim=c(0.9,1.1),pch=20,col="darkgrey",
     main="Test de martingalité pour l'immobilier",
     xlab="Maturité",
     ylab="Moyenne de l'indice actualisé ")
abline(h=S0,col="red",lty=3,lwd=2)

abline(h=S0+0.025,col="blue",lty=4)
abline(h=S0-0.025,col="blue",lty=4)
legend("topleft",legend=c("Prix initial","Borne à 5%"),
       col=c("red","blue"),pch=20,
       cex=0.8)
#========= FIN MODELE SANS RISQUE DE CREDIT ===============

