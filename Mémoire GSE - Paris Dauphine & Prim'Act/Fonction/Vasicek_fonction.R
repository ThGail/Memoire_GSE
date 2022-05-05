# définition des fonctions utiles du modèle Vasicek

# simulation taux instantané Vasicek
taux_Vas_sim <- function(N, t, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(rnorm(N, r0 * exp(-a * t) + b * (1 - exp(-a * t)), sigma * sqrt((1 - exp(-2 * a * t)) / (2 * a))))
}
taux_Vas_sim.t <- Vectorize(taux_Vas_sim,"t")

# simulation TZC
TZC_Vas_FF_sim <- function(N, t, TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  nb_mat <- length(TT)
  theta <- t(matrix(rep(TT,N),nrow=nb_mat))-t
  taux_vas <- matrix(rep(taux_Vas_sim(N,t,param,r0),nb_mat),ncol=nb_mat)
  return(Ri - ((Ri-taux_vas)*(1-exp(-a*theta))-sigma^2/(4*a^2)*(1-exp(-a*theta))^2)/(a*theta))
}
# simulation TZC des trajectoires pour une maturité donnée
TZC_Vas_FF_sim.t <- Vectorize(TZC_Vas_FF_sim,"t")

# formule fermée du TZC pour le calibrage
TZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  Ri <- b-sigma^2/(2*a^2)
  return(Ri - ((Ri-r0)*(1-exp(-a*TT))-sigma^2/(4*a^2)*(1-exp(-a*TT))^2)/(a*TT))
}
TZC_Vas_FF_calibrage.T <- Vectorize(TZC_Vas_FF_calibrage,"TT")

# deflateur 
deflateur_Vas <- function(N,t,TT,param,r0=TauxZC[1]){
  return(exp(-TT*TZC_Vas_FF_sim(N,t,TT,param,r0)))
}
deflateur_Vas.T <- Vectorize(deflateur_Vas,"TT")

# simulation PZC sur les trajectoires
PZC_Vas_FF_sim <- function(N, t, TT, param, r0=TauxZC[1]){
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  nb_mat <- length(TT)
  theta <- t(matrix(rep(TT,N),nrow=nb_mat))-t
  taux_vas <- matrix(rep(taux_Vas_sim(N,t,param,r0),nb_mat),ncol=nb_mat)
  return(exp(-b * theta) * exp(-(taux_vas - b) * (1 - exp(-a * theta)) / a + 0.5 * (sigma^2 * theta / a^2 - sigma^2 / a^3 * (1 - exp(-a * theta)) - sigma^2 / (2 * a^3) * (1 - exp(-a * theta))^2)))
}
# simulation PZC sur les trajectoires pour une maturité donnée
PZC_Vas_FF_sim.t <- Vectorize(PZC_Vas_FF_sim,"t")

# formule fermée du PZC (pour le calibrage)
PZC_Vas_FF_calibrage <- function(TT, param, r0=TauxZC[1]) {
  a <- param[1]
  b <- param[2]
  sigma <- param[3]
  return(exp(-b*TT) * exp(-(r0-b)*(1-exp(-a*TT))/a + 0.5 * (sigma^2 * TT / a^2 - sigma^2 / a^3 * (1 - exp(-a * TT)) - sigma^2 / (2 * a^3) * (1 - exp(-a * TT))^2)))
}
PZC_Vas_FF_calibrage.T <- Vectorize(PZC_Vas_FF_calibrage,"TT")

# fonction objective à minimiser
ecart_Vas <- function(param){
  e <- 0
  for (t in 1:length(TauxZC)){
    tzp<-TZC_Vas_FF_calibrage(t, param, TauxZC[1])
    # tzp <- mean(taux_Vas_sim.t(1000,t,paramVas))
    e <- e + (tzp-TauxZC[t])^2
  }
  return(e)
}