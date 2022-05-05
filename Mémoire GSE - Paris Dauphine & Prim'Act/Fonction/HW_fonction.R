# définition des fonctions utiles du modèle Hull&White

# la fonction theta(t) de la dynamique issu par interpolation
theta_HW <- function(t,param){
  a = param[1]
  sigma = param[2]
  return(dTfI_fct(t)+a*TfI_fct(t)+sigma^2/(2*a)*(1-exp(-2*a*t)))
}

beta_HW <- function(t,param){
  a = param[1]
  sigma = param[2]
  return(TfI_fct(t)+((sigma^2)/(2*a)^2)*(1-exp(-a*t))^2)
}

taux_HW_fs_sim <- function(N,t,s,param,r0=TauxZC[1]){
  a <- param[1]
  sigma <- param[2]
  return(r0*exp(-a*(t-s))+beta_HW(t,param)-beta_HW(s,param)*exp(-a*(t-s)) + sqrt((sigma^2)/(2*a)*(1-exp(-2*a*(t-s))))*rnorm(N))
}

taux_HW_sim <- function(N,t,param,r0=TauxZC[1]){
  return(taux_HW_fs_sim(N,t,0,param,r0=TauxZC[1]))
}

B_HW <- function(t, TT, param) {
  a <- param[1]
  sigma <- param[2]
  return((1 / a) * (1 - exp(-a * (TT - t))))
}

A_HW <- function(t, TT, param) {
  a <- param[1]
  sigma <- param[2]
  B <- B_HW(t,TT,param)
  return((PZC_fct(TT)/PZC_fct(t))*exp(B*TfI_fct(t)-((sigma^2)/(4*a))*(1-exp(-2*a*t))*B^2))
}

# simulation PZC avec la formule fermée
PZC_HW_sim <- function(N,t,TT,param,r0=TauxZC[1]){
  return(A_HW(t,TT,param)*exp(-B_HW(t,TT,param)*taux_HW_sim(N, t, param,r0)))
}
PZC_HW_sim.T <- Vectorize(PZC_HW_sim,"TT")
PZC_HW_sim.t <- Vectorize(PZC_HW_sim,"t")

# PZC pour la partie calibrage (t=0)
PZC_HW_FF <- function(TT,param,r0=TauxZC[1]){
  return(PZC_HW_sim(1,0,TT,param,r0))
}

# simulation TZC à partir de la formule de passage
TZC_HW_sim <- function(N,t,TT,param,r0=TauxZC[1]){
  return(-log(PZC_HW_sim(N,t,TT,param,r0))/(TT-t))
}

# TZC à partie de la formule de passage (t=0)
TZC_HW_FF <- function(TT,param,r0=TauxZC[1]){
  return(TZC_HW_sim(1,0,TT,param,r0))
}

# deflateur 
deflateur_HW <- function(N,t,TT,param,r0=TauxZC[1]){
  return(exp(-(TT-t)*TZC_HW_sim(N,t,TT,param,r0)))
}
deflateur_HW.T <- Vectorize(deflateur_HW,"TT")
deflateur_HW.t <- Vectorize(deflateur_HW,"t")

# prix d'un cap de maturite T au temps t
caps_HW_sim <- function(N,t,TT=20,param,K=K_ATM,r0=TauxZC[1]){
  # nominal = 1
  # pas de temps annuel
  # remarque t < TT
  a <- param[1]
  sigma <- param[2]
  
  cap <- 0
  t0 <- ceiling(t) # la partir entière supérieure
  tn <- TT
  for (ti in ((t0+1):tn)){
    PZCti_1 <- PZC_HW_sim(N, t, ti-1, param, r0)
    PZCti <- PZC_HW_sim(N, t, ti, param, r0)
    v <- sigma*sqrt((1-exp(-2*a * (ti-1-t)))/(2*a))*B_HW(ti-1,ti,param)
    b <- (1/v)*log(PZCti*(1 + K)/PZCti_1)+v/2
    cap <- cap + PZCti_1*pnorm(-b+v) - (1+K)*PZCti*pnorm(-b) 
  }
  return(cap)
}

caps_HW_calibrage <- function(TT,param,K=K_ATM,r0=TauxZC[1]){
  # calibrage en t=0
  # en t=0, c'est une formule FF (pas d'aléa)
  return(caps_HW_sim(N=1,t=0,TT,param,K,r0))
}

# la fonction objective pour le calibrage
ecart_HW_cap <- function(param) {
  e <- 0
  for (t in 1:length(PrixCaps)) {
    e <- e + (caps_HW_calibrage(t,param) - PrixCaps[t])^2
  }
  return(e)
}

# pour comparer l'erreur au TZC du Vasicek
ecart_HW_TZC <- function(param) {
  e <- 0
  for (t in 1:length(Maturite)) {
    e <- e + (TZC_HW_sim(1,0,Maturite[t],param) - TauxZC[t])^2
  }
  return(e)
}
