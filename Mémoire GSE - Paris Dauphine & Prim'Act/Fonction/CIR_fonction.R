# définition des fonctions utiles du modèle CIR

A_CIR <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return(((2 * ga * exp((k + ga) * (u / 2))) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))**((2 * k * mu) / (sigma**2)))
}

B_CIR <- function(u, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  ga <- sqrt(k**2 + 2 * sigma**2)
  return((2 * (exp(ga * u) - 1)) / ((k + ga) * (exp(ga * u) - 1) + 2 * ga))
}

# calcule de l'intensite de défaut
lambdat_CIR_sim <- function(N, t, param) {
  k <- param[1]
  mu <- param[2]
  sigma <- param[3]
  lambda_t <- param[4]
  if (t==0) {return(rep(lambda_t, N))}
  for (i in 1:t) {
    lambda_t <- abs(lambda_t + k * (mu - lambda_t) + sigma * sqrt(lambda_t) * rnorm(N, 0, 1))
  }
  return(lambda_t)
}

# simulation de la proba de survie
survie_CIR_sim <- function(N, t, TT, param) {
  lambda_t <- lambdat_CIR_sim(N, t, param)
  B <- B_CIR(TT - t, param)
  A <- A_CIR(TT - t, param)
  return(matrix(A*exp(- lambda_t %*% t(B)), nrow = N))
}
survie_CIR_sim.t <- Vectorize(survie_CIR_sim,"t")
survie_CIR_sim.T <- Vectorize(survie_CIR_sim,"TT")

# survie en version FF pour la calibration
survie_CIR_calibrage <- function(TT, param) {
  # le calibrage se fait en t=0
  lambda0 <- param[4]
  B <- B_CIR(TT, param)
  A <- A_CIR(TT, param)
  return(A*exp(-B*lambda0))
}

# le spread pour chaque rating
spread_CIR_sim <- function(N,t,TT,param,LGD){
  return((1 - LGD + LGD * survie_CIR_sim(N,t,TT,param) )**(-1/(TT-t))-1)
}
spread_CIR_sim.t <- Vectorize(spread_CIR_sim,"t")
spread_CIR_sim.T <- Vectorize(spread_CIR_sim,"TT")

# spread en version FF pour la calibration
spread_CIR_FF_calibrage <- function(TT,param,LGD){
  return((1-LGD+LGD*survie_CIR_calibrage(TT,param))**(-1/TT)-1)
}

# prix zéro-coupon risqué en version FF pour chaque rating à la maturité
PZCr_CIR_FF <- function(TT, param_Vas, param_CIR, LGD){
  return(PZC_Vas_FF_calibrage(TT, param_Vas)/(1+spread_CIR_FF_calibrage(TT,param_CIR,LGD))^TT)
}

# simulation zéro-coupon risqué
PZCr_CIR_Vas_sim <- function(N,t,TT,param_Vas,param_CIR,LGD){
  surv <- survie_CIR_sim(N, t, TT, param_CIR)
  return(PZC_Vas_FF_sim(N,t,TT,param_Vas)*(surv+(1-LGD)*(1-surv)))
}
# pour une maturité T fixé
PZCr_CIR_Vas_sim.t <- Vectorize(PZCr_CIR_Vas_sim,"t")
# pour t=0, à maturité variable
PZCr_CIR_Vas_sim.T <- Vectorize(PZCr_CIR_Vas_sim,"TT")

PZCr_CIR_HW_sim <- function(N,t,TT,param_HW,param_CIR,LGD){
  surv <- survie_CIR_sim(N, t, TT, param_CIR)
  return(PZC_HW_sim(N,t,TT,param_HW)*(surv+(1-LGD)*(1-surv)))
}
PZCr_CIR_HW_sim.T <- Vectorize(PZCr_CIR_HW_sim,"TT")

# fonction objective à minimiser pour le modèle CIR
ecart_CIR <- function(param){
  e <- 0
  for (t in Maturite[1:TT]){
    spread <- spread_CIR_FF_calibrage(t,param,LGD)
    e <- e + (spread-spread_mkt)^2
  }
  return(e)
}
